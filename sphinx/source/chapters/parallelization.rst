Parallelization
===============

**Authors:** Magdalena Schwarzkopf, Dominik Münch

Runtime Benchmarking
--------------------

This week we implemented a dedicated benchmarking executable to measure the runtime performance and parallelization efficiency of the 2D tsunami simulator. This benchmark is particularly useful for evaluating how well the OpenMP parallelization scales with different thread counts and for comparing performance across different hardware configurations.

Purpose
~~~~~~~

The benchmark executable (``build/benchmark``) is designed to isolate the computational core of the tsunami simulation—specifically the time-stepping loop where most of the computation occurs. Unlike the main simulator which includes setup overhead, netCDF I/O, and station output, the benchmark measures only the essential numerical computation: the ghost cell updates (``setGhostOutflow()``) and the wave propagation solver (``timeStep()``).

This isolation makes the benchmark ideal for:

- Measuring solver performance without I/O interference
- Evaluating OpenMP scaling efficiency with different thread counts
- Comparing relative performance improvements across code optimizations
- Profiling computational kernels in the wave propagation

Configuration
~~~~~~~~~~~~~~

The benchmark is configured with fixed parameters to provide consistent, reproducible measurements:

- **Grid resolution**: Fixed at 2160 × 1200 cells
- **Physical domain**: The fixed Tohoku domain extent ``[-200000, 2500000] × [-750000, 750000]`` meters
- **Solver**: F-Wave solver for 2D wave propagation
- **Setup**: Tohoku 2011 tsunami event (loaded from netCDF bathymetry and displacement files)

Building and Running
~~~~~~~~~~~~~~~~~~~~

Build the benchmark target::

    scons

Run the benchmark with OpenMP thread control::

    ./build/benchmark RUNS OMP_NUM_THREADS [END_TIME] [BATHY_NC] [DISPL_NC]

Where:

- ``RUNS``: Number of repeated benchmark runs (timing values are averaged)
- ``OMP_NUM_THREADS``: Number of OpenMP threads for the parallel solver
- ``END_TIME``: Optional simulation end time in seconds (default: 10800)
- ``BATHY_NC``: Optional bathymetry netCDF file path
- ``DISPL_NC``: Optional displacement netCDF file path

Example: Run 3 repetitions with 4 OpenMP threads::

    ./build/benchmark 3 4

Output
~~~~~~

For each run, the benchmark prints:

- Setup information (domain extents, cell sizes, time-step size, estimated number of steps)
- Progress updates during the time loop showing current simulation time and step count
- Per-run statistics:

  - ``time stepping seconds``: Total wall-clock time spent in time-stepping (excluding setup)
  - ``time steps``: Number of completed time steps
  - ``time per cell and iteration``: Average time per cell per iteration (in seconds)
  - ``time per cell and iteration in ns``: Same metric in nanoseconds for easier interpretation

- Averaged statistics across all runs

Results Summary
~~~~~~~~~~~~~~~

The following table summarizes the performance metrics across different OpenMP thread counts:

.. table:: Runtime Performance with Varying Thread Counts
   :widths: auto

   +---------+----------------+------------------------+-------------------------+---------------------+
   | Threads | Time (seconds) | Time per Cell/Iter (s) | Time per Cell/Iter (ns) | Speedup vs 1 thread |
   +=========+================+========================+=========================+=====================+
   | 1       | 433.16         | 1.57107e-08            | 15.7107                 | 1.00×               |
   +---------+----------------+------------------------+-------------------------+---------------------+
   | 8       | 76.9283        | 2.79018e-09            | 2.79018                 | 5.63×               |
   +---------+----------------+------------------------+-------------------------+---------------------+
   | 16      | 38.7993        | 1.40725e-09            | 1.40725                 | 11.17×              |
   +---------+----------------+------------------------+-------------------------+---------------------+
   | 32      | 21.0277        | 7.62671e-10            | 0.762671                | 20.61×              |
   +---------+----------------+------------------------+-------------------------+---------------------+
   | 64      | 11.3072        | 4.10112e-10            | 0.410112                | 38.33×              |
   +---------+----------------+------------------------+-------------------------+---------------------+
   | 72      | 21.2947        | 7.72356e-10            | 0.772356                | 20.34×              |
   +---------+----------------+------------------------+-------------------------+---------------------+
   | 144     | 110.373        | 4.00322e-09            | 4.00322                 | 3.92×               |
   +---------+----------------+------------------------+-------------------------+---------------------+

We can observe that the parallelization increases the speed significantly up to 64 threads after which the simulation seems to slow down again.
The NVIDIA Grace cluster we ran this on is a dual socket system with 72 cores per socket, so the performance drop at 144 threads is likely 
because staying on one socket is more efficient than using both sockets.

The drop at 72 threads could be due to a cache limit or a NUMA effect.

Grace Benchmark
~~~~~~~~~~~~~~~

The Grace benchmark (``build/benchmark_grace``) is designed to meet all the requirements for the optional task 5. It uses the 2011 Tohoku input data with 250m resolution.
It uses 10800 x 6000 cells which means a 250m resolution for the cells as well. It runs 10000 time steps with netCDF output every 100 steps.

Grace Benchmark Results
^^^^^^^^^^^^^^^^^^^^^^^

**Team:** TBD

.. rubric:: Run Configurations

:Grid Resolution: :math:`10800 \times 6000 \text{ cells}`
:Total Cells: :math:`6.48 \times 10^7`
:Time Steps: :math:`10000`

.. rubric:: Time Allocation

.. math::

   T_{\text{total}} = 464.941\,\text{s}

.. math::

   T_{\text{comp}} = 405.035\,\text{s} \quad (87.12\%)

.. math::

   T_{\text{io}} = 59.8984\,\text{s} \quad (12.88\%)

.. rubric:: Core Performance Metrics

.. math::

   \text{Total Cell Updates} = 10800 \times 6000 \times 10000 = 6.48 \times 10^{11}

.. math::

   \text{Throughput} = 1,393,726,589\,\text{Updates/s} \approx 1393.73\,\text{MCUps}

.. math::

   \text{Latency per Cell/Iteration} = 6.25054 \times 10^{-10}\,\text{s} = 0.625054\,\text{ns}


Outer vs. Inner Loop
--------------------

The 2D solver uses a single ``#pragma omp parallel`` region per time step. 
The question is whether to place the ``#pragma omp for`` directive on the inner loop (over X) or the outer loop (over Y). 

this is what the code looks like with the inner loop parallelized:
.. code-block:: cpp

   // current implementation — directive on inner loop
   #pragma omp parallel
   {
     for( t_idx l_cy = 1; l_cy <= m_nCellsY; l_cy++ ) {
       #pragma omp for schedule(static)
       for( t_idx l_cx = 0; l_cx <= m_nCellsX; l_cx++ ) { ... }
     }
   }

and this is what the code looks like with the outer loop parallelized:

.. code-block:: cpp

   // outer loop variant
   #pragma omp parallel
   {
     #pragma omp for schedule(static)
     for( t_idx l_cy = 1; l_cy <= m_nCellsY; l_cy++ ) {
       for( t_idx l_cx = 0; l_cx <= m_nCellsX; l_cx++ ) { ... }
     }
   }

Every #pragma omp for has an implicit barrier at the end - all threads must finish
their chunk of the inner loop before any thread can proceed to the next iteration.
So with the outer parallelization, there is one barrier per nested loop - threads wait
once at the end, when all threads have finished all their assigned rows. 
With the inner loop parallelization, there is 1 barrier per outer iteration per nested loop -
threads wait after every single row, because the #pragma omp for and it's implicit barrier
are inside the outer loop. For the X-sweep that is 1200 waits, for the Y-sweep 2160, so it
is way more overhead then the single barrier in the outer loop variant.

Furthermore, parallelizing the inner loop leads to race conditions between adjacent columns
assigned to differen threads. l_ceR of thread 0 may be the l_ce of thread 1, and both threads
could silmutaneously write to the same cell. 

The inner loop results lower average (2.97 ns vs 4.79 ns for outer ``static``)
is misleading: the outer-loop results are skewed by the NUMA warmup penalty
in run 0 (7.15 ns), a consequence of the serial initialization in
``benchmark.cpp``. Once pages are warm, outer ``static`` reaches 2.59 ns by
run 2 — already faster than the inner-loop average.
The init loop is serial when comparing inner and outer loop, so one thread touches
all pages, placing everything on one NUMA node. During the time loop all other threads 
access remote memory, which is slow. This leads to initialisation overhead which slows
down the first run of the outer loop. The inner loop's barrier granularity limits 
how much remote memory any thread accesses at once, accidentally softening the NUMA penalty in run 0.

Scheduling Strategies
---------------------

The ``schedule()`` clause on ``#pragma omp for`` was varied on the outer-loop
variant. Schedules tested on the outer loop:

``static``: iterations are divided into equal-sized chunks and assigned to threads
before the loop starts, therefore no runtime coordinaton is needed. 

``dynamic``: threads grab one iteration at a time from a shared work queue at runtime.
When a thread finishes its iteration it immediatly grabs the next available one. This 
leads to good load balancing (which we do not need in this part of code), but shared 
queue requires atomic operations on every grab.

``guided``: Early on threads grab large chunks, towards the end chunks shrink so that 
faster threads can steal remaining work from slower ones.

``adding the argument 4``: divides iterations into chunks of 4, thread 0 gets rows 0-3,
thread 1 gets rows 4-7, etc. This interleaves work across threads more finely.


.. code-block:: cpp

   #pragma omp for schedule(static)      // equal chunks, assigned upfront 
   #pragma omp for schedule(static, 4)   // chunks of 4 rows, round-robin
   #pragma omp for schedule(dynamic)     // threads claim 1 chunk on-demand
   #pragma omp for schedule(guided)      // decreasing chunk sizes on-demand
   #pragma omp for schedule(guided, 4)   // guided, minimum chunk size 4


Pinning Strategies
------------------

Thread pinning was tested with ``schedule(static)`` on the outer loop so that
differences reflect placement alone. The Grace cluster has 2 CPU sockets. We changed the 
OMP_PLACES={0:16} to OMP_PLACES={0:32}, allowing Grace to use both sockets, otherwise comparing
pinning strategies would not make sense, because the threads could only "spread" across the
same socket. 

.. code-block:: bash

   OMP_PLACES={0:32} OMP_PROC_BIND=close  ./build/benchmark 3 16 
   OMP_PLACES={0:32} OMP_PROC_BIND=spread ./build/benchmark 3 16
   OMP_PLACES={0:32} OMP_PROC_BIND=master ./build/benchmark 3 16

``close`` packs threads on one socket (shared cache, one memory controller);
``spread`` distributes across both sockets (doubled bandwidth, cross-socket
traffic, threads on different sockets accessing the same data pay a NUMA penalty).  
``master`` all threads are placed as close as possible to the master thread, thread 0

NUMA and First-Touch Initialization
------------------------------------

Because the benchmark's init loop is serial, all pages are allocated on one
socket. Remote threads pay a NUMA penalty on every time step — visible as
run 0 being ~2.8× slower than run 2 (7.15 vs 2.59 ns) in the ``static``
results. Enabling parallalelization of the init loop fixes this:

.. code-block:: cpp

   // benchmark.cpp — enable for NUMA-aware init
   #pragma omp parallel for schedule(static)   // must match compute schedule
   for( t_idx l_cy = 0; l_cy < i_ny; l_cy++ ) {
     for( t_idx l_cx = 0; l_cx < i_nx; l_cx++ ) {
       // first write → OS places page on this thread's local NUMA node
     }
   }

Results Summary
---------------

All runs: 16 threads, ``OMP_PLACES={0:16}``, 2160×1200 cells, 3 repetitions.

.. list-table:: Benchmark Results (ns per cell per iteration)
   :header-rows: 1
   :widths: 36 11 11 11 11

   * - Configuration
     - Run 0
     - Run 1
     - Run 2
     - Avg
   * - Inner loop, ``static`` (as implemented)
     - 3.12
     - 2.89
     - 2.89
     - 2.97
   * - Outer loop, ``static``
     - 7.15
     - 4.62
     - 2.59
     - 4.79
   * - Outer loop, ``dynamic``
     - 81.80
     - 90.89
     - 90.92
     - 87.87
   * - Outer loop, ``guided``
     - 1.53
     - 1.30
     - 1.28
     - **1.37**
   * - Outer loop, ``static,4``
     - 2.08
     - 1.66
     - 1.67
     - 1.80
   * - Outer loop, ``guided,4``
     - 1.54
     - 1.29
     - 1.29
     - 1.37
   * - ``static`` + ``close``
     - 1.78
     - 1.31
     - 1.33
     - 1.47
   * - ``static`` + ``spread``
     - 2.88
     - 2.36
     - 2.36
     - 2.53
   * - ``static`` + ``master``
     - 1.62
     - 1.32
     - 1.33
     - 1.41
   * - First-touch init, ``static``
     - 1.60
     - 1.30
     - 1.31  
     - 1.40

Conclusions
-----------

The outer loop parallelization is actually faster than the inner loop once the NUMA warmup penalty is mitigated,
which can be seen when comparing the inner loop's 2.97 ns average to the outer loop's 2.59 ns in run 2. 
The guided scheduling performs best because it balances the load while keeping scheduling overhead low. Balancing the
load does not seem important at first when looking at our code, but because of the 
.. code-block:: cpp

        if( l_leftDry && l_rightDry ) {
          continue;
        }
the work loaded is actually not perfectly balanced, since some iterations can skip part of the code. 
The static schedule with small chunks (``static,4``) also performs well, but not as good as guided.     
The dynamic schedule performs very poorly due to high scheduling overhead and lack of locality. 
The scheduling ``dynamic`` is catastrophic here (~88 ns avg): threads must atomically claim
work from a shared queue on every iteration of a tight numerical loop —
~18× slower than ``static``. The advantage of a dynamic schedule can not be used here, since
the work load is equally distributed.  ``guided`` wins: large initial chunks minimise 
scheduling overhead while shrinking tail chunks improve load balance.        

The close and master splitter strategies yield similar performance, yet the spread strategie is significantly slower.
The close and master strategies both pack all 16 threads onto the same socket, so they are euqally fast.
With the spread strategie, the threads are being distributed across both sockets, but the init loop runs on socket 0, the 8 threads 
on socket 1 have to fetch all their data across the socket interconnect on every single time step.

With NUMA's first touch policy thread 0 first-touches rows ``0…N/p``, thread 1 rows ``N/p…2N/p``, etc.,
matching the static decomposition used at compute time. With this change, 
run 0 drops to ~1.55 ns — nearly equal to run 2. (``dynamic`` must not be
used here as random chunk assignment destroys the locality guarantee.) Initialization with and without first touch policy
is being compared with the static schedule. 

The best configuration is **outer loop with** ``guided`` **scheduling** (1.37 ns avg).  
Enabling NUMA-aware initialization closes the run-0 warmup gap from 7.15 ns in the same conditions(static scheduling, outer loop) except the 
NUMA-aware initialization to 1.55 ns, demonstrating that the serial init is the main source of variance
between runs.
