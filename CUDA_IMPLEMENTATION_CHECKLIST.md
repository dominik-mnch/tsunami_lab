# CUDA Implementation Checklist - 3-Week Plan

## Week 1: CUDA Infrastructure & Baseline

### Task 1: CUDA Infrastructure Setup
- [x] Install CUDA Toolkit (match GPU driver version)
  - [x] Download CUDA toolkit from NVIDIA
  - [x] Install nvcc compiler
  - [x] Set CUDA_PATH environment variable
  - [x] Verify installation: `nvcc --version`
- [x] Install profiling tools
  - [x] Install NVIDIA Nsys
  - [x] Install nvprof
  - [x] Install Nsight Compute
  - [x] Verify tools are in PATH
- [x] Set up SCons CUDA tool support
  - [x] Modify `SConstruct` to detect CUDA installation
  - [x] Add CUDA compiler variables (NVCC, CUDA_ARCH)
  - [x] Add NVCC flags (C++17 standard, optimization, compute capability)
  - [x] Test SCons can find and use nvcc
- [x] Integrate with build system
  - [x] Modify `src/SConscript` to compile .cu files
  - [x] Create CUDA object file linking
  - [x] Build succeeds: `scons mode=release`

**Validation Checklist:**
- [x] `nvcc --version` returns valid version
- [x] `scons mode=release` builds without CUDA errors

---

### Task 2: Regression Testing Framework
- [x] Create GPU vs CPU comparison harness
  - [x] Create `src/cuda/CudaRegressionTest.h/cpp`
  - [x] Implement GPU memory allocation for solver arrays (h, hu, hv, b)
  - [x] Copy data from CPU to GPU
  - [x] Copy results back from GPU to CPU
  - [x] Create comparison function with tolerance
- [x] Implement tolerance checking
  - [x] Define floating-point tolerance constant (e.g., 1e-5 relative error)
  - [x] Create tolerance checking function for single cells
  - [x] Create comparison function for entire grids
  - [x] Log mismatches with cell coordinates and values
- [x] Test infrastructure setup
  - [x] Create test cases for 500×500 grid
  - [x] Create test cases for 1000×1000 grid
  - [x] Create test cases for 4000×4000 grid
  - [x] Ensure tests run on all grid sizes
- [x] Create benchmark logging
  - [x] Create CSV output format: `timestamp, grid_size, kernel_name, time_ms, occupancy, memory_bw`
  - [x] Create `scripts/log_benchmark.py` or logging class
  - [x] Implement append-to-CSV functionality
  - [x] Create `data/benchmark_results/` directory for outputs
- [ ] Unit tests
  - [x] Create `src/cuda/CudaRegressionTest.test.cpp`
  - [x] Test memory allocation/deallocation
  - [x] Test data copy CPU→GPU→CPU
  - [x] Test tolerance checking with known differences
  - [ ] Tests compile and pass

**Validation Checklist:**
- [x] Regression test framework compiles
- [x] All 3 grid sizes can be tested
- [x] CSV logging works and files created
- [ ] Tolerance checking catches significant differences
- [ ] Unit tests pass for memory and copying

---

### Task 3: Baseline Implementation
- [x] Implement fused single-step kernel
  - [x] Create `src/cuda/WavePropagation2d_kernels.cu`
  - [x] Implement `computeStepKernel`: one thread per cell, reads only old buffers, writes new cell once
  - [x] Accumulates net-updates from all four surrounding edges (left/right x, bottom/top y) in one kernel launch
  - [x] Handle wet/dry boundary conditions via `edgeNetUpdate` device function
  - [x] Race-free by design — no atomics needed
- [x] Implement ghost cell boundary kernels
  - [x] `setGhostOutflowLeftRightKernel`: one thread per row, copies first/last interior column to ghost columns
  - [x] `setGhostOutflowBottomTopKernel`: one thread per column, copies first/last interior row to ghost rows (runs after left/right to set corners correctly)
- [x] Implement GPU wrapper
  - [x] Create `src/cuda/WavePropagation2d_cuda.h` and `src/cuda/WavePropagation2d_cuda.cu`
  - [x] Allocate double-buffered GPU memory for h, hu, hv; single buffer for b
  - [x] `copyToGpu` / `copyToHost` for initialization and result retrieval
  - [x] `setGhostOutflow()`: launches left/right then bottom/top kernels with one `cudaDeviceSynchronize`
  - [x] `computeStep(scaling)`: launches `computeStepKernel` with one `cudaDeviceSynchronize`
  - [x] `swapBuffers()`: CPU-only flip of the active buffer index (zero-cost)
  - [x] `initialize(deviceId)` / `finalize()` static methods for CUDA device lifecycle
- [x] Integrate with build system
  - [x] `src/SConscript` compiles `.cu` files via `CUDAObject` builder under `CUDA_ENABLED` guard
  - [x] `build/benchmark_cuda` executable added to `SConstruct`, built only when CUDA is enabled
  - [x] `src/benchmark_cuda.cpp` implements 3-run benchmark using `computeStep` path, same grid (2160×1200) and input data as `benchmark`
- [x] Verify correctness
  - [x] Run kernel on small test cases via `CudaRegressionTest`
  - [x] Compare GPU output with CPU `WavePropagation2d` reference within tolerance
  - [x] Check numerical accuracy within tolerance
  - [x] Verify no segfaults or CUDA memory errors
- [x] Unit tests
  - [x] `src/cuda/CudaRegressionTest.test.cpp` tests memory allocation, copy round-trip, and multi-step regression
  - [x] Regression test passes on multiple grid sizes and timesteps

**Validation Checklist:**
- [x] All kernels compile without warnings
- [x] Kernels execute without CUDA runtime errors
- [x] GPU results match CPU results within tolerance
- [x] Numerical correctness verified across multiple timesteps
- [x] `benchmark_cuda` executable builds and runs end-to-end
- [x] Per-time-step workflow: `setGhostOutflow` → `computeStep` → `swapBuffers` (3 kernel launches, 2 syncs, 0 host↔device transfers)

---

## Week 2: Performance Benchmarking & Strategies

### Task 4: Thread Count Benchmarking
- [ ] Set up benchmark infrastructure
  - [ ] Create benchmark driver program
  - [x] Implement timer with microsecond precision
  - [x] Ensure warm-up runs before measurement
  - [x] Run multiple iterations (e.g., 10) and average
- [x] Test thread configurations
  - [x] Test block size 1
  - [x] Test block size 2
  - [x] Test block size 4
  - [x] Test block size 8
  - [x] Test block size 16
  - [x] Test block size 32
- [x] Measure on 1000m grid
  - [x] Run each configuration 10 times
  - [x] Record kernel execution time
  - [x] Record SM occupancy with nsys
  - [x] Record registers per thread
  - [x] Record memory bandwidth (GB/s)
- [x] Analyze with profiler
  - [x] Use nsys to capture timeline for each block size
  - [x] Identify register pressure threshold
  - [x] Identify occupancy limitations
  - [x] Identify memory bandwidth saturation
- [ ] Create comparison table
  - [ ] Columns: Block Size | Time (ms) | Occupancy (%) | Registers/Thread | Bandwidth (GB/s)
  - [ ] Sort by performance
  - [ ] Highlight optimal configuration
  - [ ] Document tradeoffs
- [ ] Document findings
  - [ ] Write analysis of results
  - [ ] Explain when occupancy becomes limiting
  - [ ] Explain when registers become limiting
  - [ ] Document warp efficiency for each configuration

**Validation Checklist:**
- [ ] All 5 block sizes tested successfully
- [ ] Benchmark table complete with all metrics
- [ ] Best-performing block size identified
- [ ] Profiler output collected for analysis
- [ ] Clear recommendation provided
- [ ] Results reproducible

---

### Task 5: Memory Layout Experimentation
- [ ] Implement row-major layout
  - [ ] Keep existing natural layout
  - [ ] Document layout: `m_h[row * nCellsX + col]`
  - [ ] Verify correctness with regression tests
- [ ] Implement column-major layout
  - [ ] Create alternative memory allocation
  - [ ] Adjust indexing for column-major: `m_h[col * nCellsY + row]`
  - [ ] Create kernel that works with column-major
  - [ ] Verify correctness with regression tests
- [ ] Implement blocked/tiled layout
  - [ ] Decide on tile size (32×32 or 64×64)
  - [ ] Create memory allocation for tiled format
  - [ ] Adjust indexing for tile calculations
  - [ ] Create kernel that works with tiled layout
  - [ ] Verify correctness with regression tests
- [ ] Benchmark each layout
  - [ ] Test on 500×500 grid
  - [ ] Test on 1000×1000 grid
  - [ ] Test on 4000×4000 grid
  - [ ] For each: measure time, memory bandwidth, L1/L2 cache hit rates
  - [ ] Run 10 iterations per configuration, average results
- [ ] Profile memory behavior
  - [ ] Use nsys to analyze cache hit rates
  - [ ] Analyze memory coalescing efficiency
  - [ ] Identify cache misses
  - [ ] Calculate achieved memory bandwidth
- [ ] Account for overhead
  - [ ] Measure transpose/layout conversion time
  - [ ] Measure memory padding overhead
  - [ ] Include in total time calculation
- [ ] Create visualization
  - [ ] Performance vs grid size for each layout
  - [ ] Cache hit rate comparison
  - [ ] Memory bandwidth comparison
  - [ ] Recommendation for each grid size range
- [ ] Document results
  - [ ] Write analysis of layout performance
  - [ ] Explain why certain layouts prefer certain grid sizes
  - [ ] Provide tier recommendations

**Validation Checklist:**
- [ ] All 3 layouts tested on all 3 grid sizes
- [ ] Performance table complete with all metrics
- [ ] Cache hit rates measured and documented
- [ ] Overhead accounted for in timing
- [ ] Visualization created and clear
- [ ] Recommendation provided for each grid size

---

### Task 6: Atomic vs Non-Atomic Synchronization
- [x] Implement atomic-based version
  - [x] Create kernel using `atomicAdd()` for edge updates
  - [x] Handle each neighbor cell update with atomic
  - [x] Test correctness with regression tests
  - [x] Verify all values update correctly
- [x] Implement alternative (lock-free operator-split)
  - [x] Create x-sweep kernel: one thread per cell, reads OLD buffers, writes to NEW buffer
  - [x] Create y-sweep kernel: one thread per cell, reads OLD h/hv, accumulates y-updates onto x-swept h
  - [x] Use double-buffered GPU memory for h, hu, hv (two buffers each)
  - [x] Stream-serialized kernel launches guarantee sequential consistency (no sync needed between sweeps)
  - [x] Test correctness with regression tests (all 53 test cases pass deterministically)
- [x] Benchmark both approaches
  - [x] Test on block sizes from Task 4
  - [x] Test on 2000m resolution grid
  - [x] Test on 1000m resolution grid
  - [x] Test on 500m resolution grid
  - [x] Record execution time for each
- [x] Create comparison table
  - [x] Strategy | Block Size | Grid Size | Time (ms) 
  - [x] Identify winner for each configuration
  - [x] Document performance ratio
- [x] Document tradeoffs
  - [x] Explain atomic overhead
  - [x] Explain merge overhead
  - [x] Identify grid size thresholds
  - [x] Provide recommendations

**Validation Checklist:**
- [x] Both synchronization strategies tested
- [x] Correctness verified for both
- [x] All grid sizes benchmarked
- [x] Performance winner identified for each scenario
- [x] Clear documentation of tradeoffs

---

## Week 3: Advanced Strategies & Final Report

### Task 7: Red-Black Ordering Strategy
- [ ] Implement red-black partitioning
  - [ ] Create function to classify cells as red or black (checkerboard)
  - [ ] Red: `(i + j) % 2 == 0`
  - [ ] Black: `(i + j) % 2 == 1`
  - [ ] Verify partitioning is correct
- [ ] Implement Phase 1 kernel (red cells)
  - [ ] Create kernel that updates only red cells
  - [ ] Each thread handles one edge
  - [ ] No conflicts within red phase
  - [ ] Verify correctness with regression tests
- [ ] Implement Phase 2 kernel (black cells)
  - [ ] Create kernel that updates only black cells
  - [ ] Each thread handles one edge
  - [ ] No conflicts within black phase
  - [ ] Verify correctness with regression tests
- [ ] Benchmark red-black approach
  - [ ] Test on 1000×1000 grid
  - [ ] Test on 2000×2000 grid
  - [ ] Test on 4000×4000 grid
  - [ ] Measure time for red phase
  - [ ] Measure time for black phase
  - [ ] Measure synchronization overhead between phases
  - [ ] Record total time (red + sync + black)
- [ ] Compare with atomic version
  - [ ] Run atomic version on same grids
  - [ ] Compare total kernel time
  - [ ] Compare SM efficiency
  - [ ] Measure occupancy for both
- [ ] Analyze synchronization overhead
  - [ ] Measure host-GPU synchronization cost
  - [ ] Measure kernel launch overhead (two vs one)
  - [ ] Identify if overhead outweighs benefits
- [ ] Determine winning strategy
  - [ ] Quantify performance difference
  - [ ] Identify grid size thresholds
  - [ ] Provide clear recommendation
- [ ] Document findings
  - [ ] Explain when red-black wins
  - [ ] Explain when atomic wins
  - [ ] Write analysis of synchronization bottlenecks

**Validation Checklist:**
- [ ] Red-black partitioning implemented correctly
- [ ] Both phases produce correct results
- [ ] Numerical correctness maintained
- [ ] All 3 grid sizes benchmarked
- [ ] Comparison with Task 6 results complete
- [ ] Clear winner identified for each scenario
- [ ] Overhead quantified and documented

---

### Task 8: Final Benchmarking & Optimization Report
- [ ] Run complete simulations
  - [ ] Implement y-sweep kernel (or use baseline)
  - [ ] Create benchmark that runs x-sweep + y-sweep
  - [ ] Run on 500×500 grid for 10 timesteps
  - [ ] Run on 1000×1000 grid for 10 timesteps
  - [ ] Run on 2000×2000 grid for 5 timesteps
  - [ ] Run on 4000×4000 grid for 3 timesteps
- [ ] Profile full pipeline
  - [ ] Use nsys to capture entire simulation
  - [ ] Identify remaining bottlenecks
  - [ ] Analyze kernel launch overhead
  - [ ] Analyze memory transfer overhead
  - [ ] Check for GPU stalls
- [ ] Create comprehensive table
  - [ ] Rows: all tested configurations
  - [ ] Columns: Grid Size | Block Size | Layout | Sync Strategy | Time (ms) | Speedup vs CPU
  - [ ] Calculate speedup ratio
  - [ ] Highlight best configuration
- [ ] Summarize performance impact
  - [ ] Which parameters had biggest impact?
  - [ ] Rank by performance contribution: thread count, memory layout, synchronization
  - [ ] Quantify performance gain from each decision
- [ ] Document best configuration
  - [ ] Best thread count: ___
  - [ ] Best memory layout: ___
  - [ ] Best synchronization: ___
  - [ ] Achieved speedup: ___x
  - [ ] Register usage: ___ per thread
  - [ ] SM occupancy: ___%
- [ ] Identify future optimization opportunities
  - [ ] Y-sweep kernel fusion
  - [ ] Boundary condition fusion with main computation
  - [ ] Shared memory optimization
  - [ ] Double buffering for overlapping computation/transfer
  - [ ] Kernel fusion opportunities
  - [ ] Other optimization ideas
- [ ] Write final report
  - [ ] Executive summary: key findings and recommendations
  - [ ] Methodology: what was tested and why
  - [ ] Results: comprehensive performance table
  - [ ] Analysis: which parameters matter most
  - [ ] Recommendations: best configuration for cluster
  - [ ] Future work: opportunities for further optimization
- [ ] Create deployment guide
  - [ ] Best parameters for different grid sizes
  - [ ] Performance expectations
  - [ ] Build instructions for cluster
  - [ ] Profiling commands for cluster verification
  - [ ] Troubleshooting guide

**Validation Checklist:**
- [ ] Complete simulations run successfully
- [ ] All configurations benchmarked
- [ ] Performance table accurate and complete
- [ ] Speedup calculations verified
- [ ] Best configuration clearly identified
- [ ] Future optimization opportunities documented
- [ ] Final report is comprehensive and clear
- [ ] Deployment guide is actionable

---

## Cross-Cutting Concerns

### Version Control & Documentation
- [ ] Initial commit with task 1 files
- [ ] Daily commits with progress updates
- [ ] Code review checklist for each task
- [ ] README updated with CUDA status
- [ ] GitHub/GitLab issues tracking for blockers

### Performance Metrics Tracking
- [ ] Create `data/benchmarks/` directory structure
- [ ] Week 1 baseline results saved
- [ ] Week 2 thread benchmarks saved
- [ ] Week 2 memory layout results saved
- [ ] Week 2 atomic vs lock-free results saved
- [ ] Week 3 red-black results saved
- [ ] Week 3 final comprehensive results saved
- [ ] All profiler outputs archived

### Risk Mitigation
- [ ] Fallback to CPU version always available
- [ ] Original CPU code not modified
- [ ] CUDA code is optional (build flag)
- [ ] Regression tests prevent silent failures
- [ ] Profiler output captured for debugging

---

## Summary

**Total Checklist Items:** ~150+ actionable items

**Estimated Time per Week:**
- Week 1: 15-20 hours (infrastructure heavy)
- Week 2: 20-25 hours (experimentation, iteration)
- Week 3: 15-20 hours (analysis, reporting)

**Key Success Criteria:**
- [ ] All tasks completed on schedule
- [ ] GPU speedup achieved (target: 5-10x)
- [ ] Numerical correctness maintained
- [ ] Code is production-ready
- [ ] Documentation is comprehensive
