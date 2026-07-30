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
=====================================

The experiments were performed on the GPU partition because the short
partition was fully allocated at the time. All measurements were performed
on compute nodes and not on login nodes.

The benchmark case used was a one-dimensional shock-shock problem:

* resolution: ``50m``
* domain: ``0 <= x <= 5000``
* end time: ``1.0``
* solver: ``FWave``

The measured quantity is the time per cell and iteration.

Task 1: Support for generic compilers
-------------------------------------

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

Task 2: Comparison of GCC and Clang
-----------------------------------

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

Task 3: Optimization flags and runtime comparison
-------------------------------------------------

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

Effect of optimization flags on numerical accuracy
--------------------------------------------------

The optimization level can influence both runtime and numerical behaviour.

``-O2`` enables many performance optimizations while preserving relatively
strict floating-point behaviour. It generally keeps the numerical results
more predictable.

``-Ofast`` enables more aggressive optimizations and includes fast floating
point optimizations. These optimizations may:

* reorder floating-point operations,
* use approximations that are faster but less precise,
* assume that NaNs and infinities do not occur,
* violate strict IEEE floating-point rules.

For scientific simulations such as tsunami propagation, these changes can
affect reproducibility and may introduce small differences in the numerical
solution. Therefore, ``-Ofast`` should only be used after verifying that the
resulting accuracy is acceptable.

Conclusion
----------

The experiments show that compiler choice and optimization settings can have
a significant influence on runtime.

For this benchmark:

* GCC produced faster code than Clang at ``-O2``.
* GCC ``-Ofast`` did not provide a measurable speedup over ``-O2``.
* Aggressive optimization flags should be evaluated carefully in numerical
  simulations because they may affect floating-point accuracy.

Individual Contributions
------------------------
Dominik Münch did task 8.1 this week. Magdalena Schwarzkopf did task 8.2, and 8.3 but did them at a later time, due to time constraints in this week. 