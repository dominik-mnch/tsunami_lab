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