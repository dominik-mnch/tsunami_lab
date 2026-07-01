GPU Parallelization with CUDA
=============================

**Authors:** Magdalena Schwarzkopf, Dominik Münch

For the individual phase of our project, we want to parallelize our existing tsunami simulation code using CUDA. 
We've previously implemented a (CPU) parallelized version of our solver using OpenMP but GPU parallelization is much more powerful.
To test just how much more powerful we want to implement a CUDA version of our solver and compare the performance of the two implementations to find out just how fast we can get.

Implementation Details
----------------------

GPU Wrapper and Memory Management
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The CUDA implementation is encapsulated in the ``WavePropagation2dCuda`` GPU wrapper class (located in ``src/cuda/WavePropagation2d_cuda.h`` and ``src/cuda/WavePropagation2d_cuda.cpp``), which manages all GPU-CPU interactions. The wrapper is responsible for:

- **GPU memory allocation**: Allocating device memory for water heights (``h``), x-momentum (``hu``), and y-momentum (``hv``) with double-buffering (current and next step buffers). The bathymetry (``b``) is allocated once and remains constant throughout the simulation.
- **Device initialization**: The static method ``initialize()`` queries GPU capabilities and sets up the CUDA device, returning success/failure status for error handling.
- **Host-device transfers**: Methods ``copyToGpu()`` and ``copyFromGpu()`` handle bidirectional data movement between host and device memory using ``cudaMemcpy()``.
- **Buffer swapping**: The ``swapBuffers()`` method toggles between current and next step buffers using an index (0 or 1), allowing efficient time-stepping without repeated allocations.

Memory layout uses a row-major format with padding for ghost cells on all boundaries. The stride (row size including ghost cells) is ``nCellsX + 2``, and the total number of stored values per variable is ``(nCellsX + 2) × (nCellsY + 2)``. This padding accommodates ghost cells for boundary conditions and enables adjacent cells in the same row to be accessed with predictable offsets, improving memory coalescing on CUDA.

**Double-buffering strategy**: Each variable (h, hu, hv) maintains two device arrays indexed as ``m_d_h[0]`` and ``m_d_h[1]``. After each time step, instead of copying new values back to the old buffer, only the index ``m_step`` is toggled between 0 and 1. This avoids expensive ``cudaMemcpy()`` operations and leverages the GPU's ability to handle multiple buffers efficiently. The operator-split kernels read from buffers indexed by ``1 - m_step`` (the old state) and write to buffers indexed by ``m_step`` (the new state), then ``swapBuffers()`` flips the index.

Kernel Design and Thread Organization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The CUDA implementation uses a **lock-free, operator-split approach** without atomic operations:

- **One thread per cell**: Each thread is assigned a single cell (indexed by thread and block IDs) and computes all updates for that cell. The thread ID is computed as ``threadIdx.x + blockIdx.x * blockDim.x`` in the x-direction and ``threadIdx.y + blockIdx.y * blockDim.y`` in the y-direction, with bounds checking to handle grids that don't divide evenly into blocks. This mapping eliminates write conflicts because each cell is updated by exactly one thread.

- **Operator-split x-sweep and y-sweep kernels**: 
  - The **x-sweep kernel** (``xSweepKernel``) iterates over vertical edges (those connecting cells to the left and right). For each edge at position (x, y), the kernel:
    
    1. Retrieves the left cell's state (h, hu) and the right cell's state
    2. Calls the Riemann solver (either Roe or F-Wave, selected at runtime) to compute net updates
    3. Accumulates the updates into the new-step buffer: the left cell receives the left net-update, and the right cell receives the right net-update
    4. Preserves the hv values from the old buffer (unchanged in the x-sweep)
    5. Applies wet/dry checking: if either adjacent cell is dry (h ≤ 0), special handling mirrors the water surface across the interface

  - The **y-sweep kernel** (``ySweepKernel``) is identical in structure but operates on horizontal edges. It reads h and hv values, computes y-momentum updates, and adds them to the h values that were already modified by the x-sweep.

- **Block dimensions**: Kernels use 16×16 blocks (256 threads per block). This configuration balances L2 cache efficiency (each 16×16 region of cells has predictable memory access patterns) with register pressure and occupancy. The block width is configurable at wrapper instantiation.

- **Grid and block indexing**: The grid is laid out as 2D blocks of 16×16 threads each, covering the computational domain. For a domain of size ``(nCellsX + 1) × nCellsY`` (x-sweep edges) or ``nCellsX × (nCellsY + 1)`` (y-sweep edges), the grid dimensions are calculated as ``((n + blockWidth - 1) / blockWidth, (m + blockWidth - 1) / blockWidth)``.

- **Synchronization**: Kernel launches are placed sequentially in the default CUDA stream, guaranteeing that all threads in the x-sweep complete before the y-sweep begins. This stream-serialized approach provides implicit synchronization without explicit device-side barriers. The cudaDeviceSynchronize() call (and error checking via cudaGetLastError()) is performed after each kernel launch to catch runtime errors and ensure data consistency before the next kernel.

- **Wet/dry handling**: When a cell is dry (h ≤ 0 from previous steps), special logic is applied:
  
  - The solver is still called, but instead of standard flux computation, the water surface is mirrored across the interface (so that both sides have the same water surface elevation).
  - Updates are zeroed for dry cells, preventing spurious momentum from accumulating.
  - This preserves the simulation's physical validity by preventing water from materializing in dry regions.

Alternative: Atomic-Based Approach
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An earlier version of the implementation used **atomic operations** (specifically ``atomicAdd``) to accumulate Riemann solver updates. While conceptually simpler in some ways, this approach proved problematic and was replaced by the lock-free design currently in use. Understanding the atomic approach and its limitations provides valuable insight into why the lock-free method was chosen.

**Atomic implementation structure**:

1. **Initialization kernel** (``initNewCells``): A separate kernel copies the old state (h, hu, hv) to the new buffers, effectively resetting the new buffer to the previous timestep's values.

2. **X-sweep kernel**: For each x-direction edge, the kernel computes net updates using the Riemann solver and uses ``atomicAdd()`` to accumulate the updates into the new buffer:
    
    .. code-block:: cuda
    
        atomicAdd(&d_h_new[idx], netUpdateH);
        atomicAdd(&d_hu_new[idx], netUpdateHu);
    
   This approach allows multiple threads to safely modify the same cell by serializing writes at the hardware level.

3. **Y-sweep kernel**: Similarly uses ``atomicAdd()`` to accumulate y-direction updates onto the already-modified new buffer.

**Why atomic operations are problematic for this application**:

- **Inherent serialization**: Each ``atomicAdd()`` operation serializes updates to a single memory location. When multiple threads in the same warp attempt to update the same cell (which is common at grid points where four cells meet), only one thread executes the atomic operation while others stall. This serialization defeats much of the parallelism gained by using a GPU.

- **Memory bandwidth bottleneck**: Atomic operations are extremely expensive on modern GPUs. While a regular global memory write costs ~400-500 cycles (hidden by other warp activity), an ``atomicAdd()`` can cost 1000+ cycles. In a finite-volume solver where each cell receives multiple updates (left, right, bottom, top), the overhead multiplies.

- **Race condition hazards**: The initialization-then-atomicAdd pattern creates potential inter-kernel data races that are **intermittent and nondeterministic**:
    
    - **Write-after-write (WAW) hazard**: If ``initNewCells`` completes but its writes haven't fully propagated to all cache levels when the x-sweep begins, threads may see partially-updated data. However, this doesn't happen deterministically every run.
    - **Read-after-write (RAW) hazard**: The x-sweep reads from the old buffers and writes to the new buffer with atomics. If the GPU's memory system doesn't guarantee visibility between kernels, the y-sweep might read stale values from the x-sweep's updates.
    - **Intermittent manifestation**: The race conditions manifest only probabilistically across multiple timesteps. A particular simulation run may complete successfully (same binary produces correct output), but other runs of the exact same code may produce different, incorrect results. This unpredictability makes the bug extremely difficult to track down—individual test runs often pass, but running the same test hundreds of times reveals failures. Testing under compute-sanitizer (which serializes kernel execution by inserting explicit synchronization) always passes, masking the underlying race condition. Only after running many timesteps does the race become apparent and reproducible.

- **Floating-point atomics are not universally supported**: ``atomicAdd()`` on 32-bit floats (``t_real``) is supported on modern architectures, but performance varies. 64-bit double atomics are not natively supported on many GPUs, requiring software fallbacks or workarounds.

**Comparison with lock-free approach**:

The lock-free, operator-split design avoids these pitfalls entirely:

- **No serialization**: Each thread is the sole writer of its assigned cell, so there are no conflicts, no atomic operations, and no thread stalls waiting for other threads.

- **Deterministic execution**: The sequential kernel launches (x-sweep, then y-sweep in a CUDA stream) provide implicit synchronization and eliminate inter-kernel race conditions. All x-sweep threads complete before any y-sweep thread begins.

- **Higher throughput**: Without atomic overhead, memory bandwidth is used efficiently. Reads and writes to the same cell from different kernels are ordered by stream semantics, avoiding the memory consistency issues that plagued the atomic version.

- **Simpler debugging**: Deterministic output makes correctness easier to verify. If GPU and CPU results differ, the difference is reproducible and localized.

**Trade-off: Extra passes over the domain**:

The lock-free approach requires that the entire x-sweep completes before the y-sweep begins, whereas the atomic version could theoretically overlap sweeps (in practice, it couldn't due to atomics' latency). The lock-free approach uses slightly more device memory (because all intermediate cell updates must be held in the new buffers) and makes two full passes over the grid instead of one. However, the elimination of atomic overhead and the improvement in cache locality far outweigh these costs.

**Summary**: The atomic-based approach was abandoned because it introduced nondeterminism, suffered from massive atomic serialization overhead, and created memory ordering bugs. The lock-free, operator-split approach is both faster and correct by design.

Solver Integration with Header-Only Design
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To support both CPU and GPU execution, the Riemann solvers (Roe and F-Wave) are implemented as **header-only classes** with ``__host__ __device__`` annotations. This design avoids a common pitfall in CUDA heterogeneous computing:

- **Dual-compilation and code generation**: The header files (``src/solvers/Roe.h`` and ``src/solvers/F_wave.h``) use a preprocessor macro ``TSUNAMI_CUDA_HOST_DEVICE`` that expands to ``__host__ __device__`` when compiling with NVCC (CUDA's compiler). When compiling host-only code with a standard C++ compiler (e.g., g++), the macro expands to an empty string. This allows the same source code to be compiled twice:
    
    1. Once by g++ for CPU execution (producing host-only binaries)
    2. Once by NVCC for GPU execution (producing device binaries that are linked into the CUDA runtime)

- **Inline solver methods**: All public and private static methods in the solver classes are marked ``inline`` and decorated with ``TSUNAMI_CUDA_HOST_DEVICE``. For example, the Roe solver's ``netUpdates()`` method computes the net updates by evaluating characteristic speeds and eigenstructure—all operations are inlined, so the compiler can optimize them together with the calling kernel code.

- **No function pointer overhead**: By inlining solvers directly into kernel code at compile time, there is no need for function pointers or indirect function calls on the GPU. If a separate compiled solver binary were used, the GPU would incur latency and branch divergence from indirect calls. Instead, the compiler's ptxas assembler generates device code with the full solver logic inline and optimized within each kernel.

- **Compilation pipeline**: When building with NVCC, a single source file (e.g., ``WavePropagation2d_kernels.cu``) includes both the solver headers and kernel definitions. NVCC compiles the entire translation unit, and the inline solver code is automatically instantiated in every kernel that calls it. The resulting PTX or SASS code contains the full solver logic, optimized for the specific architecture.

- **Solvers supported**: Both the Roe solver and F-Wave solver are available and can be selected at instantiation time via the ``useFWaveSolver`` flag in the wrapper constructor. A runtime branch inside the kernel checks this flag and dispatches to the appropriate solver. Because both solvers have identical interfaces, this dispatch is simple and can be optimized by the compiler to avoid redundant computation.

- **Header-only benefits**:
  
  - **Single source of truth**: The solver logic is defined once, in a header, and used by both CPU and GPU code.
  - **No linking issues**: There are no unresolved function symbols because the entire solver is compiled into every translation unit that includes it.
  - **Compiler optimizations**: The inlining and aggressive optimization (e.g., constant folding, loop unrolling) are possible because the entire solver is visible at compile time.

Regression Testing and Validation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The CUDA implementation includes comprehensive regression tests (in ``src/cuda/CudaRegressionTest.h``, ``.cu``, and ``.cpp``) to verify correctness by comparing GPU and CPU results:

- **Test infrastructure**: The ``CudaRegressionTest`` class provides static helper methods for allocating GPU memory, copying data to/from the device, comparing grids, and releasing resources. Key utilities include:
    
    - ``allocateGPUMemory()``: Allocates device arrays for h, hu, hv, and b using ``cudaMalloc()``.
    - ``copyToGPU()`` and ``copyFromGPU()``: Perform bidirectional host-device transfers with ``cudaMemcpy()``.
    - ``compareGrids()``: Compares two 2D arrays (CPU and GPU results) element-wise.
    - ``compareSingleCells()``: Validates a single cell, useful for debugging localized discrepancies.
    - ``freeGPUMemory()``: Releases all allocated device memory with ``cudaFree()``.

- **Floating-point comparison methodology**: Direct floating-point equality is unreliable due to rounding differences between CPU and GPU execution. Instead, comparisons use a **magnitude-scaled relative tolerance**:
    
    $$\text{allowed error} = \text{tolerance} \times \max(1.0, |\text{CPU value}|)$$
    
    The default tolerance is 1e-5. This accounts for the fact that large values can accumulate rounding errors, while small values near zero should be nearly exact. For multi-step tests, the tolerance is increased to 0.02 to accommodate error accumulation over time.

- **Test workflow**: Each regression test follows a standard pattern:
    
    1. Initialize a test problem (e.g., a dam-break scenario with specific bathymetry and initial conditions).
    2. Copy the initial state to the GPU device.
    3. Execute a single x-sweep and y-sweep on both GPU and CPU independently.
    4. Copy GPU results back to the host.
    5. Call ``compareGrids()`` to check h, hu, and hv arrays.
    6. If all comparisons pass within tolerance, the test passes; otherwise, it fails and reports the maximum error found.

- **Deterministic correctness across time steps**: The current lock-free implementation produces **deterministic results**—the same simulation run repeatedly yields identical output to machine precision. This is verified by running multiple time steps (e.g., 10-100 steps) and confirming that GPU and CPU results diverge minimally (~1e-5 relative error). The lock-free, operator-split design eliminates race conditions and atomics, ensuring that the computation follows a well-defined order.

- **Test coverage**: All **53 regression test cases** pass, covering:
    
    - **Small-scale problems** (1×1 to 10×10 grids): Verify correctness for minimal domains.
    - **Medium-scale problems** (32×32 to 64×64): Test typical problem sizes.
    - **Different bathymetry configurations**: Flat bottom, sloped, island, shelf.
    - **Different initial conditions**: Dry domain, wet domain, dam-break, Gaussian bump.
    - **Both solver types**: Tests for Roe and F-Wave solvers.
    - **Multiple time steps**: Validates error accumulation over simulations.

- **Test file organization**: The test code is split to handle CUDA compilation complexity:
    
    - ``CudaRegressionTest.cu``: Compiled by NVCC, contains GPU-specific functions like ``allocateGPUMemory()``, ``copyToGPU()``, kernel calls, and GPU synchronization.
    - ``CudaRegressionTest_host.cpp``: Compiled by g++, contains test orchestration, CPU reference computations, and comparison logic.
    - ``CudaRegressionTest.h``: Shared header declaring the ``CudaRegressionTest`` class interface.
    
    This split prevents symbol conflicts (e.g., duplicate definitions of ``freeGPUMemory()``) that would otherwise occur if both CPU and GPU code were in the same translation unit.


Lock-free vs Atomic Benchmark Comparison
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Both implementations were benchmarked at three spatial resolutions on the same GPU.
The tables below summarise the most relevant figures. Cell counts are derived
from the reported total time, timestep count, and time-per-cell-per-iteration.

Simulation-Level Performance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Wall-clock performance, lock-free (LF) vs atomic (A)
   :header-rows: 1
   :widths: 12 10 12 12 12 10 12

   * - Resolution
     - Cells
     - Timesteps
     - LF time (s)
     - Atomic time (s)
     - Speedup
     - Time saved
   * - 2000 m
     - ~1.01 M
     - 6,614
     - 0.354
     - 0.419
     - 1.18×
     - 15.5 %
   * - 1000 m
     - ~4.05 M
     - 13,308
     - 3.126
     - 3.596
     - 1.15×
     - 13.1 %
   * - 500 m
     - ~16.2 M
     - 26,608
     - 43.92
     - 60.19
     - 1.37×
     - 27.0 %

*Speedup* is atomic time divided by lock-free time; *time saved* is the relative
reduction in wall-clock time. The lock-free advantage grows with problem size:
each halving of the resolution quadruples the cell count and increases atomic
contention, widening the gap from ~15 % to 27 %.

GPU Kernel Time Breakdown
^^^^^^^^^^^^^^^^^^^^^^^^^^

Total kernel time summed over the whole run, in milliseconds. The atomic approach
adds a dedicated ``initNewCells`` kernel that copies the old state into the new
buffers before every step — pure overhead that the lock-free design eliminates.

.. list-table:: Total GPU kernel time (ms), summed over all timesteps
   :header-rows: 1
   :widths: 22 11 11 11 11 12 12

   * - Kernel
     - 2000 LF
     - 2000 A
     - 1000 LF
     - 1000 A
     - 500 LF
     - 500 A
   * - X-Sweep
     - 87.2
     - 95.6
     - 1089.2
     - 1149.5
     - 16931.3
     - 19055.0
   * - Y-Sweep
     - 80.4
     - 81.1
     - 1356.0
     - 1202.4
     - 17772.0
     - 18302.2
   * - ``initNewCells``
     - —
     - 49.9
     - —
     - 597.8
     - —
     - 13621.0
   * - Wet/Dry threshold
     - 32.3
     - 32.5
     - 352.3
     - 311.5
     - 8503.1
     - 8474.9
   * - Ghost cells
     - 29.8
     - 29.8
     - 64.6
     - 58.3
     - 175.0
     - 175.0
   * - **Total**
     - **229.8**
     - **288.8**
     - **2862.1**
     - **3319.5**
     - **43381.3**
     - **59628.0**

The ``initNewCells`` kernel alone accounts for 17–23 % of the atomic version's
GPU time and does no useful computation. At 500 m it costs 13.6 s — larger than
the entire wall-clock difference between the two approaches.

Synchronization and Launch Overhead
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Share of total time spent in the CUDA API
   :header-rows: 1
   :widths: 30 12 12 12 12 12 12

   * - CUDA API call
     - 2000 LF
     - 2000 A
     - 1000 LF
     - 1000 A
     - 500 LF
     - 500 A
   * - ``cudaDeviceSynchronize``
     - 55.2 %
     - 56.3 %
     - 90.4 %
     - 90.2 %
     - 98.6 %
     - 98.9 %
   * - ``cudaLaunchKernel``
     - 28.3 %
     - 29.4 %
     - 7.4 %
     - 7.8 %
     - 1.2 %
     - 1.0 %
   * - Kernels per timestep
     - 5
     - 6
     - 5
     - 6
     - 5
     - 6

At low resolution, launch and synchronization overhead dominate (the GPU is
underutilized). As the domain grows, ``cudaDeviceSynchronize`` approaches 99 %,
confirming both implementations are memory-bandwidth bound rather than
compute bound — as expected for a finite-volume solver.

Memory Transfers
^^^^^^^^^^^^^^^^^

Host-to-device transfers are essentially identical between the two approaches
(same data layout), so the speedup comes entirely from kernel efficiency, not I/O.

.. list-table:: Host-to-device memory transfers
   :header-rows: 1
   :widths: 30 16 16 16

   * - Metric
     - 2000 m
     - 1000 m
     - 500 m
   * - Data transferred (MB)
     - 16.27
     - 64.93
     - 259.47
   * - Transfer time, LF (ms)
     - 0.69
     - 3.14
     - 14.32
   * - Transfer time, Atomic (ms)
     - 0.68
     - 3.14
     - 14.25

Conclusion
^^^^^^^^^^

The lock-free approach is faster at every resolution and correct by design
(deterministic, no atomics). Its advantage scales with problem size — reaching
27 % at 500 m — driven mainly by eliminating the ``initNewCells`` copy kernel and
avoiding atomic serialization. Since finite-volume tsunami simulation is
memory-bound, removing redundant memory traffic is the most effective
optimization available.