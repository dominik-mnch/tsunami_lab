GPU Parallelization: CUDA Implementation Plan
=====================================================

Overview
--------

Week 1 establishes infrastructure and baseline performance. Week 2 explores memory layouts and optimization strategies. Week 3 implements advanced techniques and completes the y-sweep kernel. The emphasis throughout is on systematic experimentation with different configurations to understand performance tradeoffs.

Week 1 — CUDA: Setup and Infrastructure
---------------------------------------

Goal
~~~~
Establish a working CUDA development environment with infrastructure, build system, and regression testing framework. This week focuses on setup and preparation for subsequent weeks of optimization work.

Planned activities
~~~~~~~~~~~~~~~~~~

**Task 1: CUDA Infrastructure Setup**

   - Install CUDA toolkit matching target GPU drivers
   - Install profiling tools (NVIDIA Nsys, nvprof/Nsight Compute)
   - Set up SCons CUDA tool support in SConstruct/SConscript files
   - Create baseline CUDA kernel that compiles and runs with existing build system
   - Document installation process and environment variables
   - Verify CUDA installation with simple test kernel

**Task 2: Regression Testing Framework**

   - Create automated test harness that compares GPU vs CPU outputs
   - Implement tolerance checking for floating-point differences
   - Set up test infrastructure for multiple grid sizes (500×500, 1000×1000, 4000×4000)
   - Document test procedures and acceptance criteria

**Task 3: Baseline Implementation**

   - Implement baseline x-sweep kernel using "one thread per edge" approach
   - Ensure numerical correctness with simple algorithm first
   - Create minimal host code for kernel launching
   - Verify kernel compiles and runs without errors
   - Document baseline kernel design and assumptions

Validation criteria
~~~~~~~~~~~~~~~~~~~
- CUDA toolkit and profilers successfully installed and verified
- SCons build system compiles CUDA code without errors
- Regression test framework runs and compares GPU vs CPU outputs within tolerance
- Baseline kernel executes and produces correct results
- All infrastructure documented and reproducible

Deliverables
~~~~~~~~~~~~
- Installation guide with version information and environment setup
- SConstruct/SConscript modifications with CUDA tool integration
- Regression test harness with test cases for multiple grid sizes
- Baseline x-sweep kernel implementation
- Documentation of test procedures and tolerance thresholds

Week 2 — CUDA: Performance Benchmarking and Strategy Comparison
---------------------------------------------------------------

Goal
~~~~
Systematically benchmark different thread configurations and memory layouts. Evaluate atomic versus alternative synchronization strategies. Identify which parameters have the largest impact on performance.

Planned activities
~~~~~~~~~~~~~~~~~~

**Task 4: Thread Count Benchmarking**

   - Test block sizes: 64, 128, 256, 512, 1024 threads/block
   - Run each configuration on a 1000×1000 grid for consistency
   - For each configuration measure: kernel time (ms), SM occupancy (%), registers/thread, memory bandwidth (GB/s)
   - Use nsys profiler to identify register pressure and memory bottlenecks at each block size
   - Create comparison table identifying optimal thread count
   - Document when register pressure or occupancy becomes limiting factor
   - Analyze warp efficiency and identify stall reasons

**Task 5: Memory Layout Experimentation**

   - Implement three data layouts: row-major (natural), column-major, blocked/tiled
   - Benchmark each layout on 500×500, 1000×1000, 4000×4000 grids
   - For each configuration measure: kernel time, memory bandwidth, L1/L2 cache hit rates
   - Analyze which layout has best memory coalescing for different grid sizes
   - Account for any data layout overhead (transpose, padding, etc.)
   - Create visualization showing performance vs grid size for each layout
   - Identify preferred layout for small, medium, and large grids

**Task 6: Atomic vs Non-Atomic Synchronization Comparison**

   - Implement atomic-based wave propagation accumulation (atomicAdd to neighboring cells)
   - Implement alternative: each thread computes to local buffer, then merged
   - Compare execution time on all block sizes identified in Task 4
   - Measure atomic operation frequency and contention rate
   - Profile memory stalls caused by atomic serialization
   - Identify grid size threshold where atomics become bottleneck
   - Document which approach is faster for different scenarios

Validation criteria
~~~~~~~~~~~~~~~~~~~
- All 5 thread count configurations successfully benchmarked
- Clear best-performing block size identified based on occupancy and bandwidth
- All 3 memory layouts tested and documented with performance data
- Atomic contention is quantified and compared to alternatives
- Clear recommendations for which configuration to use

Deliverables
~~~~~~~~~~~~
- Benchmark table: Block Size | Time | Occupancy (%) | Registers | Bandwidth (GB/s)
- Memory layout comparison: Layout | Grid Size | Time | Bandwidth | Cache Hit Rate (%)
- Synchronization strategy comparison: Strategy | Contention | Performance Impact
- Profiler output and analysis documents
- Recommendation document for best Week 2 configuration

Week 3 — CUDA: Advanced Strategies and Complete Implementation
--------------------------------------------------------------

Goal
~~~~
Implement advanced parallelization strategies, complete the y-sweep kernel, and finalize optimizations. Deliver a fully optimized GPU solver ready for cluster deployment.

Planned activities
~~~~~~~~~~~~~~~~~~

**Task 7: Red-Black Ordering Strategy**

   - Implement red-black checkerboard partitioning (cells split into two non-conflicting sets)
   - Phase 1 kernel: update all red cells in parallel (no conflicts)
   - Phase 2 kernel: update all black cells in parallel (no conflicts)
   - Compare execution time vs atomic version from Task 6
   - Measure synchronization overhead of two kernel launches
   - Test on multiple grid sizes: 1000×1000, 2000×2000, 4000×4000
   - Determine if synchronization overhead or atomic contention is the limiting factor

**Task 8: Final Benchmarking and Optimization Report**

   - Run complete simulation (x-sweep + y-sweep) on multiple grid sizes
   - Profile full kernel pipeline to identify remaining bottlenecks
   - Create comprehensive comparison table of all tested configurations
   - Summarize which parameters had biggest performance impact
   - Document best overall configuration (thread count, memory layout, synchronization)

Validation criteria
~~~~~~~~~~~~~~~~~~~
- Red-black and atomic strategies are quantitatively compared
- Full simulation timesteps execute efficiently on GPU
- Complete configuration tested documented with reasoning

Deliverables
~~~~~~~~~~~~
- Red-black ordering kernel implementation with comparison analysis
- Complete benchmark data across all tested configurations
- Final performance report: Strategy | Grid Size | Time | Speedup vs CPU

Timeline and Success Criteria
-----------------------------

+--------+----------------------------+---------------------------------------------+
| Week   | Milestone                  | Acceptance Criteria                         |
+========+============================+=============================================+
| 1      | Setup and infrastructure   | - CUDA toolkit installed and verified       |
|        | complete                   | - Build system configured                   |
|        |                            | - Regression tests implemented              |
+--------+----------------------------+---------------------------------------------+
| 2      | Memory and synchronization | - 5 thread sizes benchmarked                |
|        | strategies documented      | - 3 memory layouts tested                   |
|        |                            | - Atomic vs alternatives quantified         |
+--------+----------------------------+---------------------------------------------+
| 3      | Complete GPU solver        | - Red-black strategy evaluated              |
|        | ready for production       | - 15+ configurations documented             |
|        |                            | - Final report delivered                    |
+--------+----------------------------+---------------------------------------------+

Expected Outcomes
~~~~~~~~~~~~~~~~~
- **Performance gain:** 5-50x speedup depending on grid size (small grids: 5-10x, large grids: 20-50x)
- **Identified bottleneck:** Memory bandwidth for most grids, register pressure for very small grids
- **Best configuration:** Likely 256-512 threads/block with row-major or tiled layout and red-black or atomic ordering
- **Configurations tested:** 15+ different thread count, memory layout, and synchronization combinations


