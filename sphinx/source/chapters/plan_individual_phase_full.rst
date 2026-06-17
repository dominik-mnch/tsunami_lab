GPU Parallelization: Formalized Plan
=====================================

Each week focuses on a distinct goal: establishing a working baseline on CUDA hardware, refining performance through benchmarking and memory-layout experiments, and finally evaluating an equivalent implementation on Apple GPUs using Metal. 

Week 1 — CUDA: Project Setup and Baseline Implementation
--------------------------------------------------------

Goal
~~~~
Establish a reliable CUDA-based development environment and provide a working baseline kernel for the wave propagation computations. The baseline will enable fair comparisons with later optimizations and with CPU results.

Planned activities
~~~~~~~~~~~~~~~~~~
1. Tooling and environment: Install a suitable CUDA toolkit and profiling utilities. At minimum, install and verify a CUDA SDK that matches the target GPU drivers. Add profiling tools such as NVIDIA Nsight Systems and Nsight Compute (or legacy `nvprof` where appropriate) so that kernel execution and memory behavior can be measured precisely.

2. Regression testing: Create an automated test harness that runs selected solver work units on both CPU and GPU, compares numerical outputs, and reports deviations. The harness should be runnable from the project root and produce deterministic results for the configured input cases.

3. Baseline kernel: Implement a straightforward GPU kernel that computes the x-sweep flux updates. The goal is correctness first; performance tuning follows in Week 2. Document the kernel launch configuration and any simplifying assumptions.

Validation criteria
~~~~~~~~~~~~~~~~~~~
- The regression harness runs GPU vs CPU and reports identical or acceptably close floating-point results for representative grids.
- The baseline kernel compiles and executes on the target GPU and is integrated into the project build and test flow.

Deliverables
~~~~~~~~~~~~
- Installation notes for CUDA and profilers.
- A repeatable test script that validates GPU results against CPU outputs.
- Source file(s) containing the baseline wave propagation kernel and the minimal host-side driver code for running it.

Week 2 — CUDA: Benchmarking and Memory-Layout Experiments
---------------------------------------------------------

Goal
~~~~
Explore configuration space and memory layouts to identify kernel launch parameters and data arrangements that maximize throughput and resource utilization on CUDA hardware.

Planned activities
~~~~~~~~~~~~~~~~~~
1. Thread-block sizing: Benchmark a range of block sizes (for example 64, 128, 256, 512, 1024 threads per block, respecting hardware limits) on representative grid sizes, initially using a 1000×1000 test grid. Record raw execution time and secondary metrics such as SM occupancy, register usage per thread, and achieved memory bandwidth.

2. Memory-layout experiments: Evaluate row-major, column-major, and tiled (blocked) memory layouts to determine which arrangement provides the best coalesced memory access patterns and cache behavior for the sweep kernels. Run these layouts on multiple grid sizes such as 500×500, 1000×1000, and 4000×4000 to reveal scaling trends.

3. Measurement and reporting: For each configuration record time-to-solution, achieved occupancy, registers per thread, and measured bandwidth. Aggregate the results into comparison tables and include a short analysis explaining why the top configurations performed best.

Validation criteria
~~~~~~~~~~~~~~~~~~~
- A clear best-performing block size or small set of block sizes is identified for each kernel variant.
- Memory-layout tradeoffs are documented, showing which layout is preferred for small, medium, and large grids.
- Benchmark scripts and data tables are stored with the project for reproducibility.

Deliverables
~~~~~~~~~~~~
- Benchmark scripts for running block-size sweeps and memory-layout comparisons.
- A comparison table summarizing performance across configurations.
- Notes explaining the observed behavior and recommended defaults.

Week 3 — Apple GPU: Metal Framework Evaluation
----------------------------------------------

Goal
~~~~
Assess the feasibility and performance characteristics of porting the wave propagation computation to Apple GPUs using Metal. This week focuses on producing a prototype, exploring threadgroup sizing, and comparing raw Metal kernels with higher-level Metal Performance Shaders (MPS) where applicable.

Planned activities
~~~~~~~~~~~~~~~~~~
1. Metal prototype: Create an Xcode Metal project that implements the wave propagation computation either via MPS primitives or a custom Metal shading language kernel. Ensure the same numerical algorithm is used so results remain comparable to CUDA and CPU versions.

2. Threadgroup experiments: Test a set of threadgroup (workgroup) sizes such as 16, 32, 64, 128, and 256 to determine which sizing provides the best execution time and occupancy on Apple hardware. Track kernel execution time and any available occupancy or pipeline metrics.

3. Atomic vs lock-free strategies: As with the CUDA experiments, implement both atomic-based accumulation and a lock-free alternative (if feasible) and measure contention and throughput for the same block/threadgroup sizes.

4. Comparative analysis: Produce a final table that juxtaposes CUDA and Metal results across the same grid sizes and configurations. Summarize which strategies and parameters appear to be hardware-dependent and which are algorithmic.

Validation criteria
~~~~~~~~~~~~~~~~~~~
- A working Metal prototype that produces results consistent with the CPU baseline.
- Threadgroup size recommendations and measured kernel times for the tested Apple GPU.
- A comparative report that highlights key differences between CUDA and Metal implementations.

Deliverables
~~~~~~~~~~~~
- Xcode project or standalone Metal kernel source implementing the wave propagation computations.
- Benchmark data and a short report comparing CUDA vs Metal performance.
- Recommendations for a production-quality implementation or a decision to favor one platform for further investment.

Cross-cutting concerns and measurement guidance
----------------------------------------------

Reproducibility and test inputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Use the same deterministic input datasets across all runs so that timing results are comparable. Document the exact input files, grid sizes, and any pre-processing steps required to generate the test cases.

Metrics to capture
~~~~~~~~~~~~~~~~~~
- Wall-clock kernel execution time (mean and standard deviation over multiple runs).
- Achieved memory bandwidth (where measurable) and SM or GPU occupancy.
- Register usage per thread and any notable spill to local memory.
- Any profiler-reported hotspots or memory-bound markers.

Performance reporting
~~~~~~~~~~~~~~~~~~~~~
Generate concise tables that highlight the best configurations. Include short commentary that explains tradeoffs and anomalies (for example, reductions in performance due to register pressure or poor memory coalescing).

Risks and mitigations
~~~~~~~~~~~~~~~~~~~~~
- Numerical drift: Floating-point results on GPU may differ slightly from CPU due to different reduction orders. Establish acceptable tolerance thresholds and verify correctness via unit/regression tests.
- Platform differences: Hardware limits (maximum threads per block, register file size) differ between GPU vendors. Keep configuration code modular so platform-specific launch parameters can be tuned easily.
- Contention costs: Atomic operations may create contention; where contention is high, assess lock-free reductions or hierarchical accumulation techniques.

Timeline and acceptance
-----------------------
- End of Week 1: Baseline CUDA kernel and working regression tests.
- End of Week 2: Benchmark results and a recommended CUDA configuration.
- End of Week 3: Metal prototype and a comparative report.

