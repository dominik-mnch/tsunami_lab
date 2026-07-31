Optimization
============

**Authors:** Magdalena Schwarzkopf, Dominik Münch

Executing code on the Draco Cluster
-----------------------------------

To execute our code on the Draco cluster, we first needed to get the code there
in the first place. This was done via git by cloning the repository,
initializing the submodules, and loading a Python module to install SCons.
A netCDF module is already available on the cluster by default.
Detailed setup, build, Slurm, and troubleshooting instructions are collected in
the :doc:`Draco cluster guide <draco_cluster_guide>`.

Results of the cluster simulations
----------------------------------

To verify that the code runs correctly on the cluster, we executed the same simulations as in
chapter 6 (Tsunami Simulations) on the cluster. We chose a 1000m and 500m resolution for both cases.
We also tried to run a 250m resolution simulation, but it did not finish within the allocated time of 24h.
The results of the 1000m and 500m resolution simulations are shown below.

Unfortunately, the 500m resolution simulation in the Chile case did not finish within the allocated time of 12h.
We therefore only show it until about 6000s of simulation time when the job was cancelled.


**Tohoku with 1000m resolution and outflow boundary conditions:**

.. video:: ../../../res/optimization/tohoku_1000m_outflow.mp4
   :align: center
   :width: 100%

**Tohoku with 500m resolution and outflow boundary conditions:**

.. video:: ../../../res/optimization/tohoku_500m_outflow.mp4
   :align: center
   :width: 100%

**Chile with 1000m resolution and outflow boundary conditions:**

.. video:: ../../../res/optimization/chile_1000m_outflow.mp4
   :align: center
   :width: 100%

**Chile with 500m resolution and outflow boundary conditions:**

.. video:: ../../../res/optimization/chile_500m_outflow.mp4
   :align: center
   :width: 100%

Comparing and evaluating the results of the local runs vs the cluster runs
--------------------------------------------------------------------------

To compare the local machine and the Draco cluster independently of the chosen
grid resolution and simulation length, we added a performance metric to the
solver. The timer only measures the numerical time stepping part of the
simulation, i.e., the repeated update of the wave propagation patch. Setup work,
input loading, and output writing are not included in this measurement.

At the end of each run, the code prints the accumulated time spent in the time
stepping loop, the number of executed time steps, and the normalized metric
``time per cell and iteration``. This value is computed as

.. math::

   \frac{t_\text{time stepping}}{N_x \cdot N_y \cdot N_\text{steps}}

where :math:`t_\text{time stepping}` is measured in seconds, :math:`N_x` and
:math:`N_y` are the numbers of grid cells in x- and y-direction, and
:math:`N_\text{steps}` is the number of time steps. The resulting unit is
seconds per cell update. A lower value means that, on average, the machine needs
less time to update one cell during one time step.

We are also counting how much real time the simulation takes and outputting that 
at the end of the run as well. We call this ``time stepping seconds``. 
This has nothing to do with the simulation time.

We can now compare which machine is faster for the same simulation.

Local runs
~~~~~~~~~~
Since we didn't want to run all the simulations again locally as well,
we only ran the Tohoku 1000m and 500m resolution simulations locally.

**Tohoku 1000m:**

- time stepping seconds: :math:`1909.94s = 31.81min`
- time steps: :math:`8872`
- time per cell and iteration: :math:`5.31549 \times 10^{-8} = 53.1549ns`

**Tohoku 500m:**

- time stepping seconds: :math:`14694.1s = 244.9min = 4.08h`
- time steps: :math:`17744`
- time per cell and iteration: :math:`5.11182 \times 10^{-8} = 51.1182ns`

Draco cluster runs
~~~~~~~~~~~~~~~~~~
Since the 500m resolution simulation in the Chile case did not finish within the allocated time of 12h,
we don't have the final performance metrics for that case.

**Tohoku 1000m:**

- time stepping seconds: :math:`1609.1s = 26.82min = 0.45h`
- time steps: :math:`8872`
- time per cell and iteration: :math:`4.47822 \times 10^{-8} = 44.7822ns`

**Tohoku 500m:**

- time stepping seconds: :math:`12989s = 216.48min = 3.61h`
- time steps: :math:`17744`
- time per cell and iteration: :math:`4.51864 \times 10^{-8} = 45.1864ns`

**Chile 1000m:**

- time stepping seconds: :math:`3885.37s = 64.76min = 1.08h`
- time steps: :math:`8206`
- time per cell and iteration: :math:`4.50933 \times 10^{-8} = 45.0933ns`

Comparison and evaluation
~~~~~~~~~~~~~~~~~~~~~~~~~

We can observe that the time per cell and iteration is lower on the Draco Cluster than on the local machine for all cases.
This is, of course, expected. However, the difference is on the order of a few nanoseconds, which is not a huge difference.

We want to make a note here about the time per cell and iteration of the Chile 1000m simulation on the cluster. Even though
it has about the same amount of time steps as the Tohoku 1000m simulation but a much larger total runtime, it still has about the same
time per cell and iteration. This is because the Chile case has a larger grid size thant the Tohoku case, which means that more cells
are updated during each time step. 

The Draco Cluster didn't give us a huge improvement but this is also expected since we are effectively only using 1 core on one of the nodes
because our code hasn't been parallelized yet. Apart from that, it should be mentioned that the local machine that ran the simulation
is already a pretty powerful machine with a good CPU (AMD Ryzen 5 9600X) so we shouldn't expect a huge improvement from the cluster.

The next chapter will be about parallelization and we should expect to see a much bigger improvement from the cluster then because we 
will actually be using more of the cores.

Compiler and Optimization Experiments
--------------------------------------

The experiments were performed on the GPU partition because the short
partition was fully allocated at the time. All measurements were performed
on compute nodes and not on login nodes.

The benchmark case used for all compiler experiments was a one-dimensional
shock-shock problem:

* resolution: ``50m``
* domain: ``0 <= x <= 5000``
* end time: ``1.0``
* solver: ``FWave``

The measured quantity is the time per cell and iteration.

Support for generic compilers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The build system was extended to support different C++ compilers using the
``CXX`` environment variable. This allows selecting the compiler during the
SCons build process, for example:

.. code-block:: bash

   CXX=g++ scons

or

.. code-block:: bash

   CXX=clang++ scons

The compiler is read from the environment in the SCons build script and
forwarded to the local SCons environment.

Comparison of GCC and Clang
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The solver was compiled with both GCC and Clang using the ``-O2``
optimization level. Each configuration was executed three times.

GCC ``-O2`` results:

+------+-----------------------------+
| Run  | Time per cell and iteration |
+======+=============================+
| 1    | 2.2988e-07                 |
+------+-----------------------------+
| 2    | 5.5373e-07                 |
+------+-----------------------------+
| 3    | 2.76615e-07                |
+------+-----------------------------+

The average runtime for GCC was:

::

   (2.2988e-07 + 5.5373e-07 + 2.76615e-07) / 3
   = 3.53e-07 seconds per cell and iteration


Clang ``-O2`` results:

+------+-----------------------------+
| Run  | Time per cell and iteration |
+======+=============================+
| 1    | 6.37281e-06                |
+------+-----------------------------+
| 2    | 5.87092e-06                |
+------+-----------------------------+
| 3    | 6.77842e-06                |
+------+-----------------------------+

The average runtime for Clang was:

::

   (6.37281e-06 + 5.87092e-06 + 6.77842e-06) / 3
   = 6.34e-06 seconds per cell and iteration


Comparison:

+----------+-----------------------------+
| Compiler | Average time                |
+==========+=============================+
| GCC      | 3.53e-07 s/cell/iteration   |
+----------+-----------------------------+
| Clang    | 6.34e-06 s/cell/iteration   |
+----------+-----------------------------+

For this benchmark, GCC generated faster code. The GCC executable was
approximately 18 times faster than the Clang executable.

The benchmark problem is relatively small, so some variation between runs
can occur due to system effects such as scheduling and cache effects.
Nevertheless, the difference between GCC and Clang was consistent across
the three Clang measurements.

Optimization flags and runtime comparison
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The effect of compiler optimization flags was investigated using GCC.
The optimization levels ``-O2`` and ``-Ofast`` were compared.

GCC ``-O2`` results:

+------+-----------------------------+
| Run  | Time per cell and iteration |
+======+=============================+
| 1    | 3.49127e-07                |
+------+-----------------------------+
| 2    | 2.61535e-07                |
+------+-----------------------------+
| 3    | 1.26583e-07                |
+------+-----------------------------+

Average:

::

   (3.49127e-07 + 2.61535e-07 + 1.26583e-07) / 3
   = 2.46e-07 seconds per cell and iteration


GCC ``-Ofast`` results:

+------+-----------------------------+
| Run  | Time per cell and iteration |
+======+=============================+
| 1    | 2.43328e-07                |
+------+-----------------------------+
| 2    | 1.31695e-07                |
+------+-----------------------------+
| 3    | 4.83725e-07                |
+------+-----------------------------+

Average:

::

   (2.43328e-07 + 1.31695e-07 + 4.83725e-07) / 3
   = 2.86e-07 seconds per cell and iteration


Comparison:

+----------------+---------------------------+
| Configuration  | Average time              |
+================+===========================+
| GCC ``-O2``    | 2.46e-07 s/cell/iteration |
+----------------+---------------------------+
| GCC ``-Ofast`` | 2.86e-07 s/cell/iteration |
+----------------+---------------------------+

For this benchmark, ``-Ofast`` did not improve the measured runtime compared
to ``-O2``. The difference between the runs is small and the benchmark size
causes noticeable measurement variation.

Additional Clang optimization tests could not be performed on this cluster
because the ``clang++`` compiler was not available in the environment. The
Clang ``-O2`` results above were obtained from an environment where Clang was
available.

Including report to find vectorization opportunities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To investigate possible vectorization opportunities, GCC optimization reports
were enabled using ``-fopt-info-vec-all``. The report was generated during
compilation with:

.. code-block:: bash

   scons report=1

The generated report analyzed the wave propagation implementation in
``WavePropagation1d``. No successful vectorizations were reported for the
inspected sections. The main reasons given by GCC for failed vectorization were:

* non-consecutive memory accesses,
* insufficient independent data references,
* missing grouped stores.

These limitations are caused by the current memory access patterns and data
dependencies in the implementation.

Including report to find vectorization opportunities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
I enabled GCC optimization reports using -fopt-info-vec-all. The report
shows that GCC analyzes the wave propagation implementation in WavePropagation1d. 
No loops were successfully vectorized in the inspected sections. The main reasons
for failed vectorization were non-consecutive memory accesses, insufficient
independent data references, and missing grouped stores. These limitations are
related to the current memory layout and data dependencies of the algorithm.

Conclusion
~~~~~~~~~~

The experiments show that compiler choice and optimization settings can have
a significant influence on runtime.

For this benchmark:

* GCC produced faster code than Clang at ``-O2``.
* GCC ``-Ofast`` did not provide a measurable speedup over ``-O2``.
* Aggressive optimization flags should be evaluated carefully in numerical
  simulations because they may affect floating-point accuracy.

.. _vtune-profiling:

VTune Profiling of the Tsunami Simulator
----------------------------------------

Build Configuration for Profiling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To obtain meaningful, line-level profiling data, the build system
(``SConstruct``) was extended with a new boolean build variable, ``profile``,
which appends the required compiler flags:

.. code-block:: python

   BoolVariable('profile',
                'add -g and disable inlining (-fno-inline) for profiling with VTune',
                False)

   # later in the file:
   if env['profile']:
     env.Append( CXXFLAGS = [ '-g', '-fno-inline' ] )

* ``-g`` embeds debug symbols, allowing VTune to map hotspots back to
  source-level functions and lines.
* ``-fno-inline`` disables function inlining, so that individually small
  but hot functions (e.g. flux/Riemann-solver kernels) show up as distinct
  entries in the profiler instead of being merged into their callers.

The benchmark target was then rebuilt with these flags enabled, alongside
the project's normal release optimizations:

.. code-block:: bash

   scons build/benchmark mode=release opt=2 profile=1

The resulting binary was verified to contain debug information:

.. code-block:: text

   $ file build/benchmark
   build/benchmark: ELF 64-bit LSB executable, x86-64, ...,
   with debug_info, not stripped

Running the Hotspots Analysis via the Command Line (Batch Job)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rather than running the analysis interactively on the login node, the
``Hotspots`` collection was launched from a Slurm batch job so that the
actual profiling work executes on an allocated compute node.

**Job script** (``profile.sh``):

.. code-block:: bash

   #!/bin/bash
   #SBATCH --job-name=vtune-hotspots
   #SBATCH -N 1
   #SBATCH -n 8
   #SBATCH -p short
   #SBATCH -t 00:30:00
   #SBATCH -o vtune_job.out
   #SBATCH -e vtune_job.err

   export PATH=/cluster/intel/oneapi/2025.0.0/vtune/2025.0/bin64:$PATH
   export OMP_NUM_THREADS=8

   cd $SLURM_SUBMIT_DIR

   vtune -collect hotspots -result-dir r001hs -- ./build/benchmark 1 8 600

The job was submitted with:

.. code-block:: bash

   $ sbatch profile.sh
   Submitted batch job 8522676

and executed successfully on ``node010`` (partition ``short``), producing a
result directory ``r001hs`` (~6.8 MB) containing the collected profiling
data. Aside from expected warnings about missing debug information for
system libraries (``libgomp``, ``libpthread``, ``libc``, ``libnetcdf``,
``libhdf5`` — none of which were compiled with ``-g``), the collection and
finalization completed without errors:

.. code-block:: text

   vtune: Collection started.
   vtune: Collection stopped.
   vtune: Executing actions 100 % done

Hotspot Analysis Results
~~~~~~~~~~~~~~~~~~~~~~~~

The benchmark executable was run on the fixed 2160 × 1200-cell Tohoku 2011
tsunami setup for a simulated end time of 600 s, using 8 OpenMP threads.

.. list-table:: Top Hotspots (Hotspots analysis, ``benchmark`` executable)
   :header-rows: 1
   :widths: 45 15 20 20

   * - Function
     - Module
     - CPU Time
     - % of CPU Time
   * - ``tsunami_lab::solvers::F_wave::netUpdates``
     - benchmark
     - 53.112 s
     - 44.7 %
   * - ``tsunami_lab::patches::WavePropagation2d::timeStep._omp_fn.0``
     - benchmark
     - 46.328 s
     - 39.0 %
   * - ``std::sqrt``
     - benchmark
     - 10.780 s
     - 9.1 %
   * - OpenMP runtime (``func@0x1dfb0``)
     - libgomp.so.1
     - 5.349 s
     - 4.5 %
   * - ``NC_get_vara`` (netCDF I/O)
     - libnetcdf.so.15
     - 1.610 s
     - 1.4 %

Additionally, VTune reported low core utilization for this run:

* Effective Physical Core Utilization: **6.5 %** (3.1 of 48 physical cores)
* Effective Logical Core Utilization: **6.4 %** (6.2 of 96 logical cores)

Which parts are compute-intensive? Was this expected?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The dominant cost, combining ``netUpdates`` and the enclosing ``timeStep``
loop, accounts for roughly **84 %** of total CPU time. This matches
expectations: these functions implement the F-Wave Riemann solver's flux
computation and the per-cell time-stepping loop, which is the numerical
core of a finite-volume shallow-water solver. Seeing the vast majority of
time spent here — rather than in I/O or runtime overhead — indicates the
profiled workload is genuinely compute-bound in the expected place.

One less expected, but notable, finding is that ``std::sqrt`` alone
accounts for **9.1 %** of total CPU time. This points to the wave-speed
computations (e.g. ``sqrt(g * h)``) inside ``netUpdates`` as a specific,
addressable cost within the dominant hotspot.

Separately, the low physical/logical core utilization (6.5 % / 6.4 %)
suggests the OpenMP parallelization is not effectively using the 8
requested threads across physical cores — a parallel-efficiency issue
distinct from the per-cell algorithmic cost discussed above.

Optimization Ideas
~~~~~~~~~~~~~~~~~~

Based on the above hotspot data, the following concrete optimization
directions were identified:

1. **Reuse of computed values.** Check whether wave-speed terms such as
   ``sqrt(g * h)`` are recomputed multiple times for the same cell across
   neighboring interface evaluations inside ``netUpdates``; caching such a
   value once per cell per time step would directly reduce the 9.1 %
   ``std::sqrt`` cost.

2. **Avoiding unnecessary square roots.** Where only a *comparison* of wave
   speeds is required (e.g. for CFL/time-step control), squared quantities
   may be usable instead of an explicit ``sqrt`` call.

3. **Re-enabling inlining.** ``-fno-inline`` was used only to get
   fine-grained profiling data. ``netUpdates`` is called once per
   cell-interface per time step, so it is exactly the kind of small, hot
   function that benefits from inlining (or being placed in a header) in
   the actual performance build; ``-fno-inline`` should be removed again
   once profiling is complete.

4. **Templates instead of runtime branching.** If ``netUpdates`` currently
   branches at runtime on solver mode (Roe vs. F-Wave) or dimensionality
   (1D vs. 2D), templating on these as compile-time parameters would let
   the compiler generate a fully specialized, branch-free version per case
   — directly inside the hottest function of the program.

5. **Investigating parallel efficiency.** Given the low core utilization
   observed, it is worth separately checking whether the OpenMP loop in
   ``timeStep`` parallelizes over a sufficiently fine-grained dimension,
   whether the grid decomposition (2160 × 1200 cells) splits evenly across
   threads, and whether thread affinity/oversubscription (8 requested
   threads on a node with 2 hardware threads per physical core) affects the
   measured utilization.

Individual Contributions
~~~~~~~~~~~~~~~~~~~~~~~~
Dominik Münch did task 8.1 this week. Magdalena Schwarzkopf did task 8.2, and 8.3 but did them at a later time, due to time constraints in this week. 