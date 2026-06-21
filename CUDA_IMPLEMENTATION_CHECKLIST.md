# CUDA Implementation Checklist - 3-Week Plan

## Week 1: CUDA Infrastructure & Baseline

### Task 1: CUDA Infrastructure Setup
- [ ] Install CUDA Toolkit (match GPU driver version)
  - [ ] Download CUDA toolkit from NVIDIA
  - [ ] Install nvcc compiler
  - [ ] Set CUDA_PATH environment variable
  - [ ] Verify installation: `nvcc --version`
- [ ] Install profiling tools
  - [ ] Install NVIDIA Nsys
  - [ ] Install nvprof
  - [ ] Install Nsight Compute
  - [ ] Verify tools are in PATH
- [ ] Set up SCons CUDA tool support
  - [ ] Modify `SConstruct` to detect CUDA installation
  - [ ] Add CUDA compiler variables (NVCC, CUDA_ARCH)
  - [ ] Add NVCC flags (C++17 standard, optimization, compute capability)
  - [ ] Test SCons can find and use nvcc
- [ ] Create baseline CUDA kernel
  - [ ] Create `src/cuda/CudaUtils.h` with error checking macros
  - [ ] Create `src/cuda/CudaTestKernel.cu` (simple vector addition)
  - [ ] Create `src/cuda/CudaTestKernel.h` with declarations
  - [ ] Verify kernel compiles without errors
- [ ] Integrate with build system
  - [ ] Modify `src/SConscript` to compile .cu files
  - [ ] Create CUDA object file linking
  - [ ] Build succeeds: `scons mode=release`
- [ ] Documentation
  - [ ] Write `docs/CUDA_SETUP.md` (installation steps)
  - [ ] Write `docs/CUDA_BUILD_INSTRUCTIONS.md` (build with CUDA)
  - [ ] Document environment variables needed
  - [ ] Document target compute capabilities

**Validation Checklist:**
- [ ] `nvcc --version` returns valid version
- [ ] `scons mode=release` builds without CUDA errors
- [ ] Test kernel compiles and runs
- [ ] Simple baseline kernel produces correct output
- [ ] Documentation is clear and reproducible

---

### Task 2: Regression Testing Framework
- [ ] Create GPU vs CPU comparison harness
  - [ ] Create `src/cuda/CudaRegressionTest.h/cpp`
  - [ ] Implement GPU memory allocation for solver arrays (h, hu, hv, b)
  - [ ] Copy data from CPU to GPU
  - [ ] Copy results back from GPU to CPU
  - [ ] Create comparison function with tolerance
- [ ] Implement tolerance checking
  - [ ] Define floating-point tolerance constant (e.g., 1e-5 relative error)
  - [ ] Create tolerance checking function for single cells
  - [ ] Create comparison function for entire grids
  - [ ] Log mismatches with cell coordinates and values
- [ ] Test infrastructure setup
  - [ ] Create test cases for 500×500 grid
  - [ ] Create test cases for 1000×1000 grid
  - [ ] Create test cases for 4000×4000 grid
  - [ ] Ensure tests run on all grid sizes
- [ ] Create benchmark logging
  - [ ] Create CSV output format: `timestamp, grid_size, kernel_name, time_ms, occupancy, memory_bw`
  - [ ] Create `scripts/log_benchmark.py` or logging class
  - [ ] Implement append-to-CSV functionality
  - [ ] Create `data/benchmark_results/` directory for outputs
- [ ] Unit tests
  - [ ] Create `src/cuda/CudaRegressionTest.test.cpp`
  - [ ] Test memory allocation/deallocation
  - [ ] Test data copy CPU→GPU→CPU
  - [ ] Test tolerance checking with known differences
  - [ ] Tests compile and pass

**Validation Checklist:**
- [ ] Regression test framework compiles
- [ ] All 3 grid sizes can be tested
- [ ] CSV logging works and files created
- [ ] Tolerance checking catches significant differences
- [ ] Unit tests pass for memory and copying

---

### Task 3: Baseline Implementation
- [ ] Implement baseline x-sweep kernel
  - [ ] Create `src/patches/WavePropagation2d_kernels.cu`
  - [ ] Implement "one thread per edge" x-sweep kernel
  - [ ] Use simple, clear algorithm (no optimizations yet)
  - [ ] Handle boundary conditions
- [ ] Implement GPU wrapper
  - [ ] Create `src/patches/WavePropagation2d_cuda.h/cpp`
  - [ ] Allocate GPU memory for state arrays
  - [ ] Implement host functions to launch x-sweep kernel
  - [ ] Copy results back to host
- [ ] Create minimal host code
  - [ ] Create `src/cuda/cuda_driver.h/cpp`
  - [ ] Implement kernel launch wrapper
  - [ ] Implement memory management interface
  - [ ] Document kernel launch parameters
- [ ] Verify correctness
  - [ ] Run baseline kernel on small test case
  - [ ] Compare output with CPU version
  - [ ] Check numerical accuracy within tolerance
  - [ ] Verify no segfaults or memory issues
- [ ] Unit tests
  - [ ] Create baseline kernel tests
  - [ ] Test on 500×500 grid
  - [ ] Test on multiple timesteps
  - [ ] Regression test passes

**Validation Checklist:**
- [ ] Baseline kernel compiles without warnings
- [ ] Kernel executes without CUDA runtime errors
- [ ] GPU results match CPU results within tolerance
- [ ] Numerical correctness verified
- [ ] Can run for multiple timesteps
- [ ] Documentation describes baseline algorithm

---

## Week 2: Performance Benchmarking & Strategies

### Task 4: Thread Count Benchmarking
- [ ] Set up benchmark infrastructure
  - [ ] Create benchmark driver program
  - [ ] Implement timer with microsecond precision
  - [ ] Ensure warm-up runs before measurement
  - [ ] Run multiple iterations (e.g., 10) and average
- [ ] Test thread configurations
  - [ ] Test block size 64 threads
  - [ ] Test block size 128 threads
  - [ ] Test block size 256 threads
  - [ ] Test block size 512 threads
  - [ ] Test block size 1024 threads
- [ ] Measure on 1000×1000 grid
  - [ ] Run each configuration 10 times
  - [ ] Record kernel execution time
  - [ ] Record SM occupancy with nsys
  - [ ] Record registers per thread
  - [ ] Record memory bandwidth (GB/s)
- [ ] Analyze with profiler
  - [ ] Use nsys to capture timeline for each block size
  - [ ] Use nvprof to get detailed metrics
  - [ ] Identify register pressure threshold
  - [ ] Identify occupancy limitations
  - [ ] Identify memory bandwidth saturation
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
- [ ] Implement atomic-based version
  - [ ] Create kernel using `atomicAdd()` for edge updates
  - [ ] Handle each neighbor cell update with atomic
  - [ ] Test correctness with regression tests
  - [ ] Verify all values update correctly
- [ ] Implement alternative (lock-free buffer)
  - [ ] Create kernel that uses thread-local buffers
  - [ ] Each thread computes updates to local array
  - [ ] After kernel, merge local buffers to global memory
  - [ ] Test correctness with regression tests
- [ ] Benchmark both approaches
  - [ ] Test on block sizes from Task 4
  - [ ] Test on 500×500 grid
  - [ ] Test on 1000×1000 grid
  - [ ] Test on 4000×4000 grid
  - [ ] Record execution time for each
- [ ] Measure atomic contention
  - [ ] Use nsys to profile atomic operations
  - [ ] Count atomic operations per kernel launch
  - [ ] Measure serialization/stalls caused by atomics
  - [ ] Identify hotspots in the code
- [ ] Analyze memory stalls
  - [ ] Profile memory stall cycles
  - [ ] Compare stall rates between approaches
  - [ ] Identify threshold where atomics become bottleneck
- [ ] Create comparison table
  - [ ] Strategy | Block Size | Grid Size | Time (ms) | Atomic Ops | Contention
  - [ ] Identify winner for each configuration
  - [ ] Document performance ratio
- [ ] Document tradeoffs
  - [ ] Explain atomic overhead
  - [ ] Explain merge overhead
  - [ ] Identify grid size thresholds
  - [ ] Provide recommendations

**Validation Checklist:**
- [ ] Both synchronization strategies tested
- [ ] Correctness verified for both
- [ ] All grid sizes benchmarked
- [ ] Atomic contention measured
- [ ] Memory stalls profiled
- [ ] Performance winner identified for each scenario
- [ ] Clear documentation of tradeoffs

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
