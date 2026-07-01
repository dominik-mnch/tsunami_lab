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

The CUDA implementation is encapsulated in the ``WavePropagation2dCuda`` wrapper
class (``src/cuda/WavePropagation2d_cuda.h`` and ``.cpp``), which owns all GPU
state and mediates every host–device interaction. It allocates device memory for
the water height ``h``, the two momentum components ``hu`` and ``hv``, and the
bathymetry ``b``; the three evolving fields are double-buffered, while the static
bathymetry is allocated once. A static ``initialize()`` method queries the device
and sets it up, ``copyToGpu()`` and ``copyFromGpu()`` move data across the PCIe
boundary with ``cudaMemcpy()``, and ``swapBuffers()`` flips the active buffer
index after each step.

The fields are stored in row-major order with a one-cell halo of ghost cells on
every side. The row stride is therefore ``nCellsX + 2`` and each array holds
``(nCellsX + 2) × (nCellsY + 2)`` values. Keeping neighbouring cells contiguous
within a row lets adjacent threads access memory with a predictable stride, which
is what allows the hardware to coalesce their loads.

Rather than copying results back over the old state at the end of every step,
each field keeps two device arrays (``m_d_h[0]`` and ``m_d_h[1]``) and the wrapper
simply toggles the index ``m_step``. The kernels read the old state from the
buffers indexed by ``1 - m_step`` and write the new state into the buffers indexed
by ``m_step``; ``swapBuffers()`` then flips the index for the next step. This
turns an expensive per-step copy into a single pointer swap.

Kernel Design and Thread Organization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The solver uses a lock-free, operator-split scheme with no atomic operations.
Work is distributed one thread per cell: a thread finds its cell from
``threadIdx`` and ``blockIdx`` in both dimensions, with a bounds check so that the
partial blocks on the trailing edges return early. Because every cell is written
by exactly one thread, there are no write conflicts to guard against.

Each step is split into an x-sweep and a y-sweep kernel. The x-sweep
(``computeXSweepKernel``) solves the Riemann problem across the two vertical edges
of its cell — using either the Roe or F-Wave solver, chosen at launch — subtracts
the resulting net-updates from ``h`` and ``hu``, and carries ``hv`` through
unchanged. The y-sweep (``computeYSweepKernel``) does the same across the two
horizontal edges: it reads its neighbours from the old buffers, then adds its
net-updates onto the height the x-sweep already wrote for that same cell and
writes the final ``hv``. Both sweeps read their neighbours only from the
read-only old buffers and each writes only its own cell, so the two passes never
race.

Kernels launch with 16×16 blocks (256 threads), and the grid is sized as
``ceil(n / blockWidth) × ceil(m / blockWidth)`` to cover the domain; the block
width is configurable at construction. Ordering between the passes comes for free
from the CUDA stream: because the kernels are enqueued sequentially on the default
stream, the x-sweep is guaranteed to finish — with its writes visible — before the
y-sweep begins, so no explicit device-side barrier is needed. Error checking with
``cudaGetLastError`` and ``cudaDeviceSynchronize`` after each launch surfaces
failures immediately.

Dry cells receive special treatment. At a wet/dry interface the water surface is
mirrored across the edge so that both sides share the same elevation, and any
update flowing into a dry cell is zeroed. This keeps water from spuriously
appearing on land while preserving the physical state of the wet region.

Alternative: Atomic-Based Approach
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An earlier version of the implementation used **atomic operations** (specifically ``atomicAdd``) to accumulate Riemann solver updates. While conceptually simpler in some ways, this approach proved problematic and was replaced by the lock-free design currently in use. Understanding the atomic approach and its limitations provides valuable insight into why the lock-free method was chosen.

**Atomic implementation structure**:

The atomic approach is *edge-parallel*: it launches one thread per cell edge, closely mirroring the structure of the CPU reference. Each time step runs three kernels, all in the same CUDA stream:

1. **Initialization kernel** (``initNewCells``): copies the old state (h, hu, hv) into the new buffers.

2. **X-sweep kernel**: one thread per vertical edge. It reads the old state, calls the Riemann solver, and accumulates the scaled net-updates into the new buffer for *both* adjacent cells using ``atomicAdd``:

   .. code-block:: cuda

       atomicAdd( &io_h[l_ce],   -i_scaling * l_nu[0][0] );  // left cell
       atomicAdd( &io_hu[l_ce],  -i_scaling * l_nu[0][1] );
       atomicAdd( &io_h[l_ceR],  -i_scaling * l_nu[1][0] );  // right cell
       atomicAdd( &io_hu[l_ceR], -i_scaling * l_nu[1][1] );

3. **Y-sweep kernel**: the same pattern for horizontal edges, accumulating onto the buffer already updated by the x-sweep.

**Why the initialization kernel is required**:

The edge-parallel design is exactly what forces ``initNewCells`` to exist. Every cell lies between two edges in each sweep, so several different threads contribute to the same cell. A plain write (``h_new[i] = h_old[i] + delta``) would let those threads overwrite one another (a lost update), so the contributions must instead be *accumulated* with ``atomicAdd`` — a read-modify-write. But accumulation only adds the deltas; it needs the base value ``h_old[i]`` to already be present in the destination. Since the new buffer is one half of a double buffer, before initialization it still holds stale data from two steps ago, not ``h_old``. ``initNewCells`` seeds the destination with the old state so the subsequent atomic accumulations land on the correct base:

.. math::

   Q_{\text{new}}[i] = \underbrace{Q_{\text{old}}[i]}_{\text{initNewCells}} \; - \; \frac{\Delta t}{\Delta x} \sum_{\text{edges around } i} \text{netUpdate}

The lock-free design removes this kernel because it is *cell-parallel*: one thread owns one cell, gathers the contributions from all surrounding edges itself, and writes the complete result — base term included — exactly once. With no accumulation into memory, there is nothing to pre-seed.

**Why the atomic approach is slower**

The atomic approach is slower for three reasons. First, the ``initNewCells`` pass
reads and writes the entire grid every step while doing no physics; profiling
shows it consumes 17–23 % of the atomic version's GPU kernel time (13.6 s at 500 m
resolution), by itself larger than the total wall-clock gap between the two
approaches. Second, because adjacent edges share a cell, many threads issue
``atomicAdd`` to the same address; the hardware serializes these conflicting
updates, and the extra read-modify-write traffic to global memory lowers effective
bandwidth, so the atomic sweep kernels are measurably slower than their lock-free
counterparts even though they perform the same arithmetic. Third, the results are
not reproducible: ``atomicAdd`` accumulates in whatever order the scheduler runs
the threads, and floating-point addition is not associative, so the same binary
can produce slightly different results from one run to the next. Those ULP-scale
differences are harmless in a single step but accumulate over thousands of steps,
particularly near shocks. The lock-free kernel instead sums each cell in a fixed
order inside a single thread, making its results bit-for-bit deterministic.

**A note on correctness (what is *not* a problem)**:

The atomic kernels are **not** affected by classic inter-kernel race conditions. All kernels are launched into the same (default) CUDA stream, and CUDA guarantees that each kernel completes — with its memory writes visible — before the next one begins. So ``initNewCells`` always finishes before the x-sweep reads the buffers, and the x-sweep always finishes before the y-sweep. Within a single sweep, concurrent ``atomicAdd`` calls to a shared cell are also safe: atomics guarantee no lost updates. The only genuine correctness concern is the floating-point non-determinism described above — not memory-visibility hazards.

**Comparison with the lock-free approach**

Compared with the atomic version, the lock-free, operator-split design avoids
these pitfalls entirely. Each thread is the sole writer of its cell, so there are
no conflicts, no atomics, and no threads stalling on one another. Each cell's
net-updates are summed in a fixed order within one thread, making the result
bit-for-bit reproducible — a genuine advantage over the atomic version, whose
accumulation order, and therefore floating-point result, varies between runs. And
without the redundant ``initNewCells`` pass or atomic read-modify-write traffic,
global-memory bandwidth (the dominant cost in this memory-bound solver) is used
more efficiently. Deterministic output also makes debugging easier, since any
divergence from the CPU reference is reproducible and localized.

**Trade-off: redundant edge computation**:

The lock-free (cell-parallel) design is not free of costs. Because each thread owns one cell and reads all of its surrounding edges, every interior edge's Riemann problem is solved **twice** — once by the cell on each side — whereas the edge-parallel atomic kernel solves it only once. The lock-free version therefore does roughly twice the flux arithmetic. In this solver that trade is clearly worth it: the redundant arithmetic is cheap compared with the memory traffic it saves (no ``initNewCells`` copy, no atomic read-modify-write), and the problem is memory-bound, so eliminating memory operations matters far more than reducing compute. Both approaches use the same double-buffered memory footprint.

**Summary**: The atomic-based approach was abandoned because it required a redundant per-step initialization kernel, incurred atomic-contention overhead, and produced non-reproducible (run-to-run varying) results due to non-associative floating-point accumulation. It was *not*, however, affected by inter-kernel memory races — same-stream launches are safely ordered. The lock-free, operator-split approach is faster and deterministic by design.

Solver Integration with Header-Only Design
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To support both CPU and GPU execution, the Riemann solvers (Roe and F-Wave) are
implemented as header-only classes annotated for both host and device. The
headers (``src/solvers/Roe.h`` and ``src/solvers/F_wave.h``) define a macro
``TSUNAMI_CUDA_HOST_DEVICE`` that expands to ``__host__ __device__`` under NVCC and
to nothing under a plain host compiler, so the same source compiles twice — once
by g++ for the CPU build and once by NVCC for the device build.

Every solver method is ``inline`` and carries that annotation, which lets the
compiler inline the full solver body directly into each kernel. There are no
function pointers or indirect device calls to pay for; instead ptxas emits the
solver logic inline and optimizes it together with the surrounding kernel code.
Because a single ``.cu`` translation unit (such as
``WavePropagation2d_kernels.cu``) includes both the solver headers and the
kernels, the solver is instantiated wherever it is used and the resulting device
code carries the complete, optimized logic.

Both solvers are available and selected at construction through the
``useFWaveSolver`` flag; a runtime branch in the kernel dispatches to the chosen
one, and since the two share an identical interface the branch is trivial. Keeping
the solver header-only gives a single source of truth shared by the CPU and GPU
builds, avoids unresolved-symbol link errors, and exposes the whole solver to the
compiler for aggressive inlining and constant folding.

Regression Testing and Validation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The CUDA implementation is validated by regression tests (in
``src/cuda/CudaRegressionTest.h``, ``.cu`` and ``.cpp``) that compare the GPU
results against the CPU reference. The ``CudaRegressionTest`` class collects the
supporting utilities: ``allocateGPUMemory()`` and ``freeGPUMemory()`` manage
device buffers for h, hu, hv and b, ``copyToGPU()`` / ``copyFromGPU()`` move data
across the host\u2013device boundary, and ``compareGrids()`` / ``compareSingleCells()``
check the results element-wise (the latter is handy for pinning down a localized
discrepancy).

Because floating-point rounding differs between CPU and GPU, the comparison uses a
magnitude-scaled relative tolerance rather than exact equality:

.. math::

   \text{allowed error} = \text{tolerance} \times \max(1.0,\ |\text{CPU value}|)

The default tolerance is 1e-5, so large values may absorb some rounding while
values near zero must match closely; for multi-step tests it is relaxed to 0.02 to
allow for error accumulation over time. Each test seeds a problem (for example a
dam-break with a given bathymetry), copies it to the device, advances one x- and
y-sweep on both the GPU and the CPU independently, copies the GPU result back, and
compares the h, hu and hv fields — reporting the maximum error on failure.

Because the lock-free design sums each cell in a fixed order and uses no atomics,
its output is deterministic: repeated runs are identical to machine precision, and
multi-step runs stay within ~1e-5 of the CPU. All 53 test cases pass, spanning
small grids (1\u00d71 to 10\u00d710) and larger ones (32\u00d732 to 64\u00d764), a range of
bathymetries (flat, sloped, island, shelf) and initial conditions (dry, wet,
dam-break, Gaussian bump), both the Roe and F-Wave solvers, and multi-step
integration.

The test sources are split to keep the CUDA and host compilers apart:
``CudaRegressionTest.cu`` (NVCC) holds the GPU primitives and kernel launches,
``CudaRegressionTest_host.cpp`` (g++) holds the test orchestration and CPU
reference, and ``CudaRegressionTest.h`` declares the shared interface. The split
avoids symbol clashes (such as a duplicate ``freeGPUMemory()``) that would arise if
CPU and GPU code shared one translation unit.


Thread-Count (Block-Size) Benchmarking
--------------------------------------

To justify the default launch configuration, the lock-free solver was benchmarked
across every supported CUDA block width on a single fixed grid. All runs use the
Tohoku setup at 1000 m resolution (2700×1500 ≈ 4.05 M cells), the lock-free
strategy, and identical inputs — only the block size changes.

A block width ``N`` launches ``N × N`` threads per block. Because a CUDA warp is
32 threads wide, the block width directly controls how many of those 32 lanes are
actually used.

Results
~~~~~~~

.. list-table:: Block-size sweep (1000 m grid, lock-free solver)
   :header-rows: 1
   :widths: 10 12 12 16 16 12

   * - Block width
     - Threads/block
     - Warps/block
     - ns / cell / iter
     - Kernel ms / run
     - Speedup vs 1×1
   * - 1×1
     - 1
     - 1/32 of a warp
     - 1.5982
     - 85911.9
     - 1.0×
   * - 2×2
     - 4
     - 1/8 of a warp
     - 0.4051
     - 21611.5
     - 3.9×
   * - 4×4
     - 16
     - ½ warp
     - 0.1155
     - 6086.7
     - 13.8×
   * - 8×8
     - 64
     - 2 warps
     - 0.0669
     - 3458.1
     - 23.9×
   * - **16×16**
     - **256**
     - **8 warps**
     - **0.0565**
     - **2941.6**
     - **28.3×**
   * - 32×32
     - 1024
     - 32 warps
     - 0.0568
     - 2982.1
     - 28.1×

Per-kernel time and the CUDA-API overhead split from the same runs:

.. list-table:: Per-kernel time (ms/run) and API overhead
   :header-rows: 1
   :widths: 10 12 12 12 12 12

   * - Block width
     - X-sweep ms
     - Y-sweep ms
     - Wet/Dry ms
     - Sync %
     - Launch %
   * - 1×1
     - 28613.8
     - 28623.6
     - 28602.6
     - 99.6
     - 0.3
   * - 2×2
     - 7188.6
     - 7201.4
     - 7169.4
     - 98.7
     - 1.2
   * - 4×4
     - 1995.8
     - 2227.7
     - 1814.5
     - 95.6
     - 3.8
   * - 8×8
     - 1268.5
     - 1627.7
     - 506.3
     - 92.4
     - 6.5
   * - 16×16
     - 1123.6
     - 1400.4
     - 350.8
     - 91.4
     - 7.4
   * - 32×32
     - 1178.9
     - 1225.9
     - 439.1
     - 91.5
     - 7.3

Analysis
~~~~~~~~

**Warp utilization dominates the low end.** The jump from 1×1 to 16×16 is a
**28×** speedup — almost exactly the 32-thread warp width. This is not a
coincidence: a 1×1 block launches a single thread, so 31 of the warp's 32 lanes
sit idle and the GPU runs at ~1/32 efficiency. Widening the block fills the warp:
2×2 uses 4 lanes (≈4× faster), 4×4 fills half a warp (≈14×), and 8×8 is the first
size that fully occupies whole warps (≈24×). Below one full warp, performance
tracks lane occupancy almost linearly.

**The sweet spot is 16×16.** Once blocks are several warps wide the curve
flattens, and 16×16 (256 threads, 8 warps) gives the best time at
0.0565 ns/cell/iter. Going further to 32×32 (1024 threads, the hardware maximum
per block) gives no improvement and a slight regression (0.0568): the larger
blocks reduce the number of blocks that can be co-resident on each SM, so
occupancy and latency-hiding no longer improve. Interestingly the **y-sweep**
alone keeps getting faster at 32×32 (1225.9 ms vs 1400.4 ms at 16×16), but the
**x-sweep** and **wet/dry** kernels regress, so the overall step is marginally
slower.

**Memory-bound, as expected.** ``cudaDeviceSynchronize`` dominates the API time
at every block size (99.6 % → 91.4 %), confirming the solver is bandwidth-bound
rather than compute-bound. As the kernels get faster with wider blocks, launch
overhead grows in relative terms (0.3 % → 7.4 %) simply because the useful work
shrinks while the fixed per-launch cost stays constant. The host-to-device
transfer volume (~130 MB over the two profiled runs) and bandwidth (~21.6 GB/s)
are constant across all block sizes — a sanity check that only the launch
configuration changed.

**Recommendation.** 16×16 (256 threads/block) is the best configuration and is
used as the default throughout the solver. 8×8 is within ~18 % and is a
reasonable fallback on GPUs with tighter register limits; anything below 4×4
should be avoided because it cannot fill a warp.

Lock-free vs Atomic Benchmark Comparison
----------------------------------------

Both implementations were benchmarked at three spatial resolutions on the same GPU.
The tables below summarise the most relevant figures. Cell counts are derived
from the reported total time, timestep count, and time-per-cell-per-iteration.

.. note::

   These runs use the **default block size** of 16×16, i.e. **256 threads per
   block**. The number of blocks scales with the domain: the grid is launched with
   ``ceil((nCellsX + 2) / 16) × ceil((nCellsY + 2) / 16)`` blocks, so the total
   thread count is roughly one thread per cell (plus a partial block of padding
   threads on the trailing edges, which return early). The block size was not
   tuned for these measurements — a different block width could shift the absolute
   timings, but the lock-free vs atomic *comparison* holds since both variants use
   the identical launch configuration.

Simulation-Level Performance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Wall-clock performance, lock-free (LF) vs atomic (A)
   :header-rows: 1
   :widths: 10 8 9 10 10 12 12 8 10

   * - Resolution
     - Cells
     - Timesteps
     - LF time (s)
     - Atomic time (s)
     - LF ns/cell/iter
     - A ns/cell/iter
     - Speedup
     - Time saved
   * - 2000 m
     - ~1.01 M
     - 6,614
     - 0.354
     - 0.419
     - 0.0528
     - 0.0625
     - 1.18×
     - 15.5 %
   * - 1000 m
     - ~4.05 M
     - 13,308
     - 3.126
     - 3.596
     - 0.0580
     - 0.0667
     - 1.15×
     - 13.1 %
   * - 500 m
     - ~16.2 M
     - 26,608
     - 43.92
     - 60.19
     - 0.1019
     - 0.1396
     - 1.37×
     - 27.0 %

*Speedup* is atomic time divided by lock-free time; *time saved* is the relative
reduction in wall-clock time; *ns/cell/iter* is the time per cell per iteration.
The lock-free advantage grows with problem size:
each halving of the resolution quadruples the cell count and increases atomic
contention, widening the gap from ~15 % to 27 %.

GPU Kernel Time Breakdown
~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
   * - Other (one-time setup/teardown)
     - 16.5 %
     - 14.3 %
     - 2.2 %
     - 2.0 %
     - 0.2 %
     - 0.1 %
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
compute bound — as expected for a finite-volume solver. The *other* row is
almost entirely the one-time ``cudaDeviceReset`` at shutdown (plus ``cudaMalloc``
/ ``cudaFree`` / ``cudaMemcpy``); it is a fixed cost that does not scale with the
simulation, which is why it looks large at 2000 m (a 0.35 s run) but vanishes to
~0.1 % at 500 m. The first three rows sum to 100 %; the last row is a count, not
a percentage.

Memory Transfers
~~~~~~~~~~~~~~~~

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
~~~~~~~~~~

The lock-free approach is faster at every resolution and correct by design
(deterministic, no atomics). Its advantage scales with problem size — reaching
27 % at 500 m — driven mainly by eliminating the ``initNewCells`` copy kernel and
avoiding atomic serialization. Since finite-volume tsunami simulation is
memory-bound, removing redundant memory traffic is the most effective
optimization available.