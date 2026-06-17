GPU Parallelization: 3-Week Timeline
=====================================================

Week 1: CUDA - Thread Configuration Experimentation
====================================================

**Task 1:** Set up CUDA project with infrastructure
  - Install CUDA toolkit and profiling tools (nsys, nvprof)
  - Create regression testing framework (GPU vs CPU comparison)
  - Implement baseline kernel for x-sweep flux computation

**Task 2:** Benchmark different thread counts
  - Test block sizes: 64, 128, 256, 512, 1024 threads
  - Run on 1000×1000 grid for each configuration
  - Measure: execution time, SM occupancy, register usage, memory bandwidth
  - Create comparison table and identify best configuration

**Task 3:** Compare atomic vs lock-free strategies
  - Implement atomic-based flux accumulation
  - Test on same block sizes as Task 2
  - Measure atomic contention and compare performance

Week 2: CUDA - Memory and Algorithmic Strategies
=================================================

**Task 4:** Experiment with memory layouts
  - Test row-major, column-major, and tiled data arrangements
  - Benchmark each on 500×500, 1000×1000, 4000×4000 grids
  - Measure memory bandwidth and cache hit rates
  - Document tradeoffs

**Task 5:** Optimize with shared memory
  - Implement shared memory caching for bathymetry/state
  - Test different shared memory allocations (32KB, 64KB, 96KB)
  - Compare against baseline global-memory-only kernel
  - Measure occupancy and bandwidth improvements

**Task 6:** Evaluate red-black ordering vs atomics
  - Implement red-black checkerboard synchronization strategy
  - Compare execution time and SM stalls vs atomic version
  - Test on multiple grid sizes
  - Determine which strategy is faster and why

Week 3: Apple GPU - Metal Framework Evaluation
===============================================

**Task 7:** Prototype Metal Performance Shaders (MPS)
  - Set up Xcode Metal project
  - Implement flux computation using MPS or custom Metal shader
  - Test on 1000×1000 grid
  - Measure GPU utilization, bandwidth, kernel time

**Task 8:** Test Metal compute kernels with threadgroup experimentation
  - Write raw Metal shading language kernel
  - Test threadgroup sizes: 16, 32, 64, 128, 256
  - Measure execution time and occupancy for each size
  - Compare with MPS results

**Task 9:** Comparative analysis and documentation
  - Create final performance comparison table (CUDA vs Apple across grid sizes)
  - Document which configurations/strategies performed best
  - Summarize key findings: which parameters matter most, unexpected results
  - Finalize code with CPU fallback and profiler data
