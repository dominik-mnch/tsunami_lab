###########
Tsunami Lab
###########

This is the code of the Tsunami Lab Project taught at Friedrich Schiller University Jena.

This repository is the work of Dominik Münch and Magdalena Schwarzkopf based on the initial code by Alexander Breuer which can be found `here <https://github.com/breuera/tsunami_lab>`_.

Further information is available from: https://scalable.uni-jena.de/opt/tsunami/

To view the user documentation and project reports, click `here <https://dominik-mnch.github.io/tsunami_lab/sphinx/>`_.

To view the code documentation, click `here <https://dominik-mnch.github.io/tsunami_lab/doxygen/>`_.

Usage
=====

Command Line Interface
----------------------

The simulator is invoked with the following syntax::

   ./build/tsunami_lab NX NY X_LOWER X_UPPER Y_LOWER Y_UPPER K SOLVER_MODE END_TIME [GPU_FLAG] PROPAGATION [PROP_PARAMS] SETUP [SETUP_PARAMS]
   ./build/tsunami_lab <RES>m  X_LOWER X_UPPER Y_LOWER Y_UPPER K SOLVER_MODE END_TIME [GPU_FLAG] PROPAGATION [PROP_PARAMS] SETUP [SETUP_PARAMS]

**General Arguments:**

- ``NX``: Number of cells in x-direction (positive integer).
- ``NY``: Number of cells in y-direction (positive integer); forced to 1 for 1D propagation.
- ``<RES>m``: Resolution in meters (e.g. ``250m``). Replaces ``NX``/``NY``; cell counts are derived as ``NX = round((X_UPPER - X_LOWER) / RES)`` and ``NY = round((Y_UPPER - Y_LOWER) / RES)``.
- ``X_LOWER``: Lower bound of the domain in x-direction (real).
- ``X_UPPER``: Upper bound of the domain in x-direction (real, must be greater than ``X_LOWER``).
- ``Y_LOWER``: Lower bound of the domain in y-direction (real).
- ``Y_UPPER``: Upper bound of the domain in y-direction (real, must be greater than ``Y_LOWER``).
- ``K``: Coarsening factor; ``k*k`` cells are averaged into one output cell (positive integer, default: ``1`` for no coarsening).
- ``SOLVER_MODE``: ``1`` for F-Wave solver, ``0`` for Roe solver.
- ``END_TIME``: Simulation end time in seconds (positive real).
- ``GPU_FLAG`` (optional, 2D only): Enable GPU acceleration for wave propagation computation. Options: ``0``, ``1``, ``true``, or ``false`` (default: ``false``). GPU acceleration is only available for 2D simulations.

**Propagation Modes:**

- ``1d <boundary_mode>``: One-dimensional wave propagation. Boundary modes:
  
  - ``outflow``: Ghost outflow boundary conditions (default).
  - ``right``: Reflecting boundary on the right, outflow on the left.
  - ``left``: Reflecting boundary on the left, outflow on the right.
  - ``both``: Reflecting boundaries on both sides.

- ``2d [boundary_mode]``: Two-dimensional wave propagation with rectangular grids (NX × NY).

   - ``outflow``: Outflow boundaries on all four sides (default).
   - ``bothx``: Reflecting boundaries on left and right sides.
   - ``bothy``: Reflecting boundaries on bottom and top sides.
   - ``all``: Reflecting boundaries on all four sides.
   - ``left,right,bottom,top``: Any comma-separated subset of sides can be reflective.

**Setups:**

1. ``dam_break_1d hL huL hR huR xDam``
   
   Discontinuous dam break with water heights and momenta on left and right sides.
   
   - ``hL``, ``hR``: Water heights (left, right).
   - ``huL``, ``huR``: Momenta in x-direction (left, right).
   - ``xDam``: x-coordinate of the dam location.

2. ``shock_shock_1d h hu xDisc [useParabola [offset center scale]]``
   
   Shock-shock Riemann problem with optional parabolic bathymetry.
   
   - ``h``, ``hu``: Initial water height and momentum.
   - ``xDisc``: x-coordinate of discontinuity.
   - ``useParabola``: ``1`` or ``true`` to enable bathymetry; ``0`` or ``false`` (optional, default: false).
   - ``offset``, ``center``, ``scale``: Bathymetry parabola parameters (optional, only with useParabola=true).

3. ``rare_rare_1d h hu xDisc [useParabola [offset center scale]]``
   
   Rare-rare (fan-fan) Riemann problem with optional parabolic bathymetry.
   
   - Same parameters as shock_shock_1d.

4. ``subcritical_1d``
   
   Subcritical shallow water flow (predefined setup, no parameters).

5. ``supercritical_1d``
   
   Supercritical shallow water flow (predefined setup, no parameters).

6. ``tsunami_event_1d``
   
   Tsunami event from bathymetry file (predefined setup, no parameters).

7. ``circular_dam_break_2d hIn hOut huIn hvIn huOut hvOut xMid yMid radius``
   
   Two-dimensional circular dam break in a rectangular domain.
   
   - ``hIn``, ``hOut``: Water heights inside and outside the dam.
   - ``huIn``, ``hvIn``: Momenta inside the dam (x and y directions).
   - ``huOut``, ``hvOut``: Momenta outside the dam (x and y directions).
   - ``xMid``, ``yMid``: Center coordinates of the circular dam.
   - ``radius``: Radius of the dam.

8. ``artificial_tsunami_2d``

   Two-dimensional artificial tsunami setup (predefined setup, no parameters).

9. ``tsunami2d i_bathymetryFile [i_displacementFile]``

   Two-dimensional tsunami event setup with netCDF input files.
   
   - ``i_bathymetryFile``: Path to the bathymetry netCDF file (required, must have ``.nc`` extension).
   - ``i_displacementFile``: Path to the displacement netCDF file (optional, must have ``.nc`` extension if provided).

**Examples:**

1. 1D dam break with F-Wave solver::

   ./build/tsunami_lab 100 1 0 10000 0 0 1 1 2.5 1d outflow dam_break_1d 10 0 5 0 5000

2. 1D shock-shock with reflecting boundary on the right::

   ./build/tsunami_lab 50 1 0 5000 0 0 1 1 1.0 1d right shock_shock_1d 5 2.5 2500

3. 1D rare-rare with parabolic bathymetry::

   ./build/tsunami_lab 100 1 0 10000 0 0 1 1 5.0 1d outflow rare_rare_1d 4 -1 5000 1 0.5 5000 -0.0001

4. 1D subcritical flow (Roe solver)::

   ./build/tsunami_lab 200 1 0 50000 0 0 1 0 10.0 1d outflow subcritical_1d

5. 2D circular dam break::

   ./build/tsunami_lab 50 40 0 5000 0 4000 1 1 2.0 2d circular_dam_break_2d 10 5 0 0 0 0 2500 2000 1000

6. 2D circular dam break with reflecting left/right boundaries::

   ./build/tsunami_lab 50 40 0 5000 0 4000 1 1 2.0 2d bothx circular_dam_break_2d 10 5 0 0 0 0 2500 2000 1000

7. 2D artificial tsunami::

   ./build/tsunami_lab 100 100 -1000 1000 -1000 1000 1 1 5.0 2d artificial_tsunami_2d

8. 2D tsunami event with bathymetry file only::

   ./build/tsunami_lab 500 400 0 100000 0 80000 1 1 10.0 2d tsunami2d res/artificialtsunami_bathymetry_1000.nc

9. 2D tsunami event with bathymetry and displacement files::

   ./build/tsunami_lab 500 400 0 100000 0 80000 1 1 10.0 2d tsunami2d res/artificialtsunami_bathymetry_1000.nc res/artificialtsunami_displ_1000.nc

10. Tohoku 2011 at 1 km resolution using resolution mode::

   ./build/tsunami_lab 1000m -200000 2500000 -750000 750000 1 1 7200 2d outflow tsunami2d data/output/tohoku_gebco20_ucsb3_250m_bath.nc data/output/tohoku_gebco20_ucsb3_250m_displ.nc

11. 2D circular dam break with GPU acceleration::

   ./build/tsunami_lab 50 40 0 5000 0 4000 1 1 2.0 true 2d circular_dam_break_2d 10 5 0 0 0 0 2500 2000 1000

12. 2D Tohoku tsunami with GPU acceleration::

   ./build/tsunami_lab 500 400 0 100000 0 80000 1 1 10.0 1 2d tsunami2d res/artificialtsunami_bathymetry_1000.nc res/artificialtsunami_displ_1000.nc

Runtime Benchmark
-----------------

The benchmark executable is intended for runtime measurements of the 2D Tohoku tsunami setup without writing simulation output. It uses a fixed grid of ``2160 x 1200`` cells, the fixed Tohoku domain ``[-200000, 2500000] x [-750000, 750000]``, the F-Wave solver, and applies the same time-step timing as the main executable: only ``setGhostOutflow()`` and ``timeStep()`` are measured inside the simulation loop.

Build the benchmark target with::

   scons build/benchmark

Run it with::

   ./build/benchmark RUNS OMP_NUM_THREADS [END_TIME] [BATHY_NC] [DISPL_NC]

**Arguments:**

- ``RUNS``: Number of repeated benchmark runs.
- ``OMP_NUM_THREADS``: Number of OpenMP threads used by the solver.
- ``END_TIME``: Optional simulation end time in seconds. Defaults to ``10800``.
- ``BATHY_NC``: Optional bathymetry netCDF file. Defaults to ``data/output/tohoku_gebco20_usgs_250m_bath.nc``.
- ``DISPL_NC``: Optional displacement netCDF file. Defaults to ``data/output/tohoku_gebco20_usgs_250m_displ.nc``.

Example using the default dataset and end time for three repeated runs with four OpenMP threads::

   ./build/benchmark 3 4

Example with a shorter simulation end time and four OpenMP threads::

   ./build/benchmark 3 4 3600

Before entering the time loop, the benchmark prints the benchmark domain, cell sizes, time-step size, and estimated number of time steps. During the time loop it prints progress updates. For each run, the benchmark prints ``time stepping seconds``, ``time steps``, ``time per cell and iteration``, and ``time per cell and iteration in ns``. It also prints the average of these values over all runs.

Grace System Benchmark
----------------------

The Grace benchmark is designed for performance evaluation on high-performance computing systems. It benchmarks the Tohoku 2011 M 9.1 tsunami event at 250 meter resolution with periodic file I/O.

**Benchmark Configuration:**

- **Event:** Tohoku 2011 M 9.1
- **Resolution:** 250 meters
- **Grid:** 10,800 × 6,000 cells
- **Domain:** [-200000, 2500000] × [-750000, 750000] meters
- **Time Steps:** Fixed at 10,000 steps
- **File I/O:** NetCDF output written every 100 time steps
- **Solver:** F-Wave with OpenMP parallelization

Build the Grace benchmark with::

   scons build/benchmark_grace

Run it with::

   ./build/benchmark_grace [BATHY_NC] [DISPL_NC] [OUTPUT_NC]

**Optional Arguments:**

- ``BATHY_NC``: Bathymetry netCDF file. Defaults to ``data/output/tohoku_gebco20_usgs_250m_bath.nc``.
- ``DISPL_NC``: Displacement netCDF file. Defaults to ``data/output/tohoku_gebco20_usgs_250m_displ.nc``.
- ``OUTPUT_NC``: Output netCDF file. Defaults to ``solutions/grace_benchmark.nc``.

Example using default settings::

   ./build/benchmark_grace

The benchmark outputs:

- **Total Elapsed Time:** Total wall-clock time for the entire benchmark
- **Computation Time:** Time spent in solver computations (setGhostOutflow and timeStep)
- **I/O Time:** Time spent writing NetCDF output
- **Total Cell Updates:** Grid size × time steps = 648,000,000 updates
- **Cell Updates per Second:** Performance metric in cell updates per second (CUps)
- **Mega Cell Updates per Second:** Performance metric in MCUps (millions of cell updates per second)