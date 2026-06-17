GPU Parallelization: Bullet Point Plan
=======================================

Week 1: CUDA - Setup and Baseline Implementation
--------------------------------------------------

**Task 1:** Set up CUDA project with infrastructure
  - Install CUDA toolkit and profiling tools (nsys, nvprof)
  - Create regression testing framework (GPU vs CPU comparison)
  - Implement baseline kernel for x-sweep flux computation

**Task 2:** Compare atomic vs lock-free strategies
  - Implement atomic-based flux accumulation
  - Test on same block sizes as Task 2
  - Measure atomic contention and compare performance

Week 2: CUDA - Benchmark and memory layouts
-------------------------------------------

**Task 3:** Benchmark different thread counts
  - Test block sizes: 64, 128, 256, 512, 1024 threads
  - Run on 1000×1000 grid for each configuration
  - Measure: execution time, SM occupancy, register usage, memory bandwidth
  - Create comparison table and identify best configuration

**Task 4:** Experiment with memory layouts
  - Test row-major, column-major, and tiled data arrangements
  - Benchmark each on 500×500, 1000×1000, 4000×4000 grids
  - Measure memory bandwidth and cache hit rates
  - Document tradeoffs

Week 3: Apple GPU - Metal Framework Evaluation
----------------------------------------------

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

**Task 9:** Compare atomic vs lock-free strategies
  - Implement atomic-based flux accumulation
  - Test on same block sizes as Task 2
  - Measure atomic contention and compare performance

**Task 10:** Comparative analysis and documentation
  - Create final performance comparison table (CUDA vs Apple across grid sizes)
  - Document which configurations/strategies performed best
  - Summarize key findings: which parameters matter most, unexpected results
  - Finalize code with CPU fallback and profiler data
