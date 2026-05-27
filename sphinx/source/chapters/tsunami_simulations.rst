Tsunami Simulations
===================

**Authors:** Magdalena Schwarzkopf, Dominik Münch

This week we did some simulations of real events. To do this we used the ``TsunamiEvent2d`` setup, which initializes the simulation using bathymetry displacement data loaded from the 
provided netCDF files. We simulated the 2011 Tohoku earthquake which occured off the coast of Japan as well as the 2010 Chile earthquake.
The results are summarized in the following project report.

Unfortunately, it was not possible for us to do a 250 m resolution simulation. The 500 m resolution prouduced output files on the order of 300 GB.
It was therefore not an option to run this simulation locally on our machines. When we tried to run the simulation on the machines in the lab using ssh, 
the simulation was unfortunately killed after a few hours of runtime because the client disconnected. Therefore we only ran 1000 m and 500 m resolution simulations.

Chile Event
-----------
To understand the Chile event, we want to first understand the input data we are working with.
We can see the bathymetry visualised from the file ``chile_gebco20_usgs_250m_bath.nc`` in the following figure:

.. image:: ../../res/tsunami_simulations/chile_bathymetry.png
   :align: center

To see where the earthquake occurred, we can overlay the bathymetry file and displacement file (``chile_gebco20_usgs_250m_displ.nc``), warped by the bathymetry and colored by displacement:

.. image:: ../../res/tsunami_simulations/chile_displacement.png
   :align: center

Domain and simulation setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The 2010 Chile earthquake (:math:`M_w` 8.8) occurred on 27 February 2010 off the Maule coast at
approximately 35.8°S, 72.7°W.  The simulation domain uses a local Cartesian coordinate system
centred on the earthquake epicenter:

.. list-table::
   :header-rows: 1
   :widths: 20 25 25 20

   * - Axis
     - Lower bound
     - Upper bound
     - Extent
   * - :math:`x` (west–east)
     - −3 000 km
     - +500 km
     - 3 500 km
   * - :math:`y` (south–north)
     - −1 500 km
     - +1 500 km
     - 3 000 km

Bathymetry and displacement data are provided by GEBCO 2020 (USGS source model) at 250 m
native resolution:

* ``data/output/chile_gebco20_usgs_250m_bath.nc``
* ``data/output/chile_gebco20_usgs_250m_displ.nc``

Computational requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~

The cell-update count per time step follows from the 2D sweep (four Riemann problems per
interior cell plus one per ghost-cell edge on each boundary row/column):

.. math::

   N_\text{updates} = 4\,n_x\,n_y + 2\,n_x + 2\,n_y

The time step is CFL-limited by the shallowest wet cell.  At 500 m resolution the
observed step size is :math:`\Delta t \approx 10.97\,\text{s}`. At 1 000 m it is approximately
:math:`\Delta t \approx 21.94\,\text{s}` (CFL scales linearly with :math:`\Delta x`).

.. list-table::
   :header-rows: 1
   :widths: 15 22 20 25 18

   * - Resolution
     - Grid (:math:`n_x \times n_y`)
     - Cells
     - Updates / step
     - Steps (7 200 s)
   * - 1 000 m
     - 3 500 × 3 000
     - 10 500 000
     - :math:`4 \times 3500 \times 3000 + 2 \times 3500 + 2 \times 3000 = 42\,013\,000`
     - ≈ 328
   * - 500 m
     - 7 000 × 6 000
     - 42 000 000
     - :math:`4 \times 7000 \times 6000 + 2 \times 7000 + 2 \times 6000 = 168\,026\,000`
     - ≈ 657

Total cell updates (steps × updates/step):

.. list-table::
   :header-rows: 1
   :widths: 20 35 45

   * - Resolution
     - Expression
     - Total cell updates
   * - 1 000 m
     - :math:`328 \times 42\,013\,000`
     - :math:`\approx 1.38 \times 10^{10}`
   * - 500 m
     - :math:`657 \times 168\,026\,000`
     - :math:`\approx 1.10 \times 10^{11}`

Simulation commands
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # 1000 m resolution, outflow boundaries, 7200 s
   ./build/tsunami_lab 1000m -3000000 500000 -1500000 1500000 1 7200 2d outflow tsunami2d \
       data/output/chile_gebco20_usgs_250m_bath.nc \
       data/output/chile_gebco20_usgs_250m_displ.nc

   # 1000 m resolution, reflecting boundaries, 7200 s
   ./build/tsunami_lab 1000m -3000000 500000 -1500000 1500000 1 7200 2d all tsunami2d \
       data/output/chile_gebco20_usgs_250m_bath.nc \
       data/output/chile_gebco20_usgs_250m_displ.nc

   # 500 m resolution, outflow boundaries, 7200 s
   ./build/tsunami_lab 500m -3000000 500000 -1500000 1500000 1 7200 2d outflow tsunami2d \
       data/output/chile_gebco20_usgs_250m_bath.nc \
       data/output/chile_gebco20_usgs_250m_displ.nc

Visualisations
~~~~~~~~~~~~~~

**1 000 m – outflow boundaries**

.. video:: ../../res/tsunami_simulations/chile_1000m_outflow.mp4
   :align: center
   :width: 100%

**1 000 m – reflecting boundaries**

.. video:: ../../res/tsunami_simulations/chile_1000m_reflecting.mp4
   :align: center
   :width: 100%

**500 m – outflow boundaries**

.. video:: ../../res/tsunami_simulations/chile_500m_outflow.mp4
   :align: center
   :width: 100%


Tohoku Event
------------

Once again, we want to inspect the input data first. The bathymetry (file:
``tohoku_gebco20_ucsb3_250m_bath.nc``) is visualised in the following figure:

.. image:: ../../res/tsunami_simulations/japan_bathymetry.png
   :align: center

To see where the earthquake occurred, we can overlay the bathymetry file and displacement file (``tohoku_gebco20_ucsb3_250m_displ.nc``), warped by the bathymetry and colored by displacement:

.. image:: ../../res/tsunami_simulations/tohoku_displacement.png
   :align: center

Domain and simulation setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The 2011 Tohoku earthquake (:math:`M_w` 9.1) occurred on 11 March 2011 off the
Pacific coast of Japan at approximately 38.3°N, 142.4°E.  The simulation domain
uses a local Cartesian coordinate system centred on the earthquake epicenter:

.. list-table::
   :header-rows: 1
   :widths: 20 25 25 20

   * - Axis
     - Lower bound
     - Upper bound
     - Extent
   * - :math:`x` (west–east)
     - −200 km
     - +2 500 km
     - 2 700 km
   * - :math:`y` (south–north)
     - −750 km
     - +750 km
     - 1 500 km

Bathymetry and displacement data are provided by GEBCO 2020 (UCSB3 source model)
at 250 m native resolution:

* ``data/output/tohoku_gebco20_ucsb3_250m_bath.nc``
* ``data/output/tohoku_gebco20_ucsb3_250m_displ.nc``

The first wave leaves the computational domain after 210 frames (30 fps).  Since
each time step is approximately 10.1457 s of simulated time, this corresponds to
about 2 131 s (≈ 35.5 min).

Computational requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the same cell-update formula as for Chile:

.. math::

   N_\text{updates} = 4\,n_x\,n_y + 2\,n_x + 2\,n_y

The observed time step at 1 000 m resolution is :math:`\Delta t \approx 10.15\,\text{s}`;
at 500 m it is approximately :math:`\Delta t \approx 5.07\,\text{s}`.

.. list-table::
   :header-rows: 1
   :widths: 15 22 20 25 18

   * - Resolution
     - Grid (:math:`n_x \times n_y`)
     - Cells
     - Updates / step
     - Steps
   * - 1 000 m
     - 2 700 × 1 500
     - 4 050 000
     - :math:`4 \times 2700 \times 1500 + 2 \times 2700 + 2 \times 1500 = 16\,208\,400`
     - ≈ 710 (7 200 s)
   * - 500 m
     - 5 400 × 3 000
     - 16 200 000
     - :math:`4 \times 5400 \times 3000 + 2 \times 5400 + 2 \times 3000 = 64\,816\,800`
     - ≈ 2 130 (10 800 s)

Total cell updates (steps × updates/step):

.. list-table::
   :header-rows: 1
   :widths: 20 35 45

   * - Resolution
     - Expression
     - Total cell updates
   * - 1 000 m
     - :math:`710 \times 16\,208\,400`
     - :math:`\approx 1.15 \times 10^{10}`
   * - 500 m
     - :math:`2130 \times 64\,816\,800`
     - :math:`\approx 1.38 \times 10^{11}`

Simulation commands
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # 1000 m resolution, outflow boundaries, 7200 s
   ./build/tsunami_lab 1000m -200000 2500000 -750000 750000 1 7200 2d outflow tsunami2d \
       data/output/tohoku_gebco20_ucsb3_250m_bath.nc \
       data/output/tohoku_gebco20_ucsb3_250m_displ.nc

   # 500 m resolution, outflow boundaries, 10800 s
   ./build/tsunami_lab 500m -200000 2500000 -750000 750000 1 10800 2d outflow tsunami2d \
       data/output/tohoku_gebco20_ucsb3_250m_bath.nc \
       data/output/tohoku_gebco20_ucsb3_250m_displ.nc

   # 1000 m resolution, reflecting boundaries, 7200 s
   ./build/tsunami_lab 1000m -200000 2500000 -750000 750000 1 7200 2d all tsunami2d \
       data/output/tohoku_gebco20_ucsb3_250m_bath.nc \
       data/output/tohoku_gebco20_ucsb3_250m_displ.nc

   # 500 m resolution, reflecting boundaries, 10800 s
   ./build/tsunami_lab 500m -200000 2500000 -750000 750000 1 10800 2d all tsunami2d \
       data/output/tohoku_gebco20_ucsb3_250m_bath.nc \
       data/output/tohoku_gebco20_ucsb3_250m_displ.nc

Visualisations
~~~~~~~~~~~~~~

**1 000 m – outflow boundaries**

.. video:: ../../res/tsunami_simulations/tohoku_1000m_outflow.mp4
   :align: center
   :width: 100%

**1 000 m – reflecting boundaries**

.. video:: ../../res/tsunami_simulations/tohoku_1000m_reflecting.mp4
   :align: center
   :width: 100%

**500 m – outflow boundaries**

.. video:: ../../res/tsunami_simulations/tohoku_500m_outflow.mp4
   :align: center
   :width: 100%

**500 m – reflecting boundaries**

.. video:: ../../res/tsunami_simulations/tohoku_500m_reflecting.mp4
   :align: center
   :width: 100%

Sōma Travel Time Estimate
~~~~~~~~~~~~~~~~~~~~~~~~~~

Sōma is a coastal city in Fukushima Prefecture, approximately 55 km south and
128 km west of the March 11, 2011, :math:`M_w` 9.1 Tohoku earthquake epicenter
(38.297°N, 142.373°E, 05:46:24 UTC).  The straight-line distance from the
epicenter to Sōma is therefore

.. math::

   d = \sqrt{128^2 + 55^2} \approx 139 \text{ km}.

Measured data
^^^^^^^^^^^^^

Tide-gauge records for the 2011 Tohoku event are available through the
`NOAA NCEI Global Historical Tsunami Database
<https://www.ngdc.noaa.gov/hazel/view/hazards/tsunami/runup-search>`_.
Searching for location *Soma*, Japan, year 2011 yields a measured first-wave
travel time of **9 minutes** and a maximum run-up water height of ~9.3 m.
The maximum wave arrives after **1 hour and 4 minutes**.

Bathymetry cut analysis
^^^^^^^^^^^^^^^^^^^^^^^

The file ``soma_epicenter_tohoku_2011_bathymetry.csv`` contains a
ParaView cross-section of the GEBCO08 bathymetry along the straight line
connecting the epicenter and Sōma. The columns are:

* ``Points:0`` — :math:`x`-coordinate (m, positive east)
* ``Points:1`` — :math:`y`-coordinate (m, positive north)
* ``z`` — bathymetry value (m; positive = above sea level, negative = ocean depth)

The coordinate origin (0, 0) coincides with the epicenter.  Two key points
are identified from the data:

.. list-table::
   :header-rows: 1
   :widths: 25 15 20 20 20

   * - Point
     - Index
     - :math:`x` (m)
     - :math:`y` (m)
     - :math:`z` (m)
   * - Sōma coastline
     - 267
     - −125 000
     - −53 487
     - ≈ 0
   * - Epicenter
     - 356
     - −1 000
     - −428
     - −969

The coastline is identified as the first sign change in :math:`z` (land →
ocean).  The epicenter is the cut point closest to the coordinate origin; it
lies only 1.1 km from (0, 0), confirming the coordinate system is centred on
the epicenter.  The ocean path length along the cut is **134.9 km**, consistent
with the geometric distance above.

Wave-speed rule of thumb
^^^^^^^^^^^^^^^^^^^^^^^^

For a shallow-water tsunami the wave speed is approximated by

.. math::

   \lambda \approx \sqrt{g \cdot h},

where :math:`h` is the water depth and :math:`g = 9.81\,\text{m/s}^2`.
Using the depth at the epicenter (:math:`h = 969\,\text{m}`, read from the
CSV at the point closest to the origin) as a single representative speed:

.. math::

   \lambda \approx \sqrt{9.81 \cdot 969} \approx 97.5\,\text{m/s}

.. math::

   T \approx \frac{d}{\lambda} = \frac{134\,900\,\text{m}}{97.5\,\text{m/s}}
           \approx 1384\,\text{s} \approx \mathbf{23\,\text{min}}

This estimate noticeably exceeds the NOAA-measured first-wave arrival of **9 minutes**.
The discrepancy arises because the calculation uses the epicenter depth (969 m) as a
single representative speed, whereas the real tsunami radiates outward and its
fastest-arriving energy travels through the **Japan Trench**, which reaches depths of
7 000–9 000 m:

.. math::

   c_\text{trench} = \sqrt{9.81 \times 7000} \approx 262\,\text{m/s},\qquad
   T_\text{trench} = \frac{139\,000}{262} \approx 531\,\text{s} \approx \mathbf{9\,\text{min}}

The straight-line cross-section only samples the shallower continental slope near the
epicenter, so it systematically underestimates the wave speed and overestimates the
travel time.  In two dimensions, energy takes the fastest available path — here through
the deep trench — rather than the straight-line path.

Simulated station output
^^^^^^^^^^^^^^^^^^^^^^^^

The figure below shows the water height and momenta recorded at the Sōma
station from the 1 000 m resolution simulation.  The dashed vertical line
marks the simulated first arrival at :math:`t = 464.3\,\text{s}` (≈ 7.7 min),
after which the wave drops for a while (wave trough) before rising to a maximum
height of around 9 m, reached at around 70 minutes.

These values coincide well with the measured data (first arrival at 9 minutes,
maximum height of 9.3 m after 1 hour and 4 minutes).  Note that the plot shows
the total water column height :math:`h`; the surface elevation anomaly
:math:`\eta = h + b` at the station cell (bathymetry :math:`b \approx -20\,\text{m}`)
yields a peak of :math:`\approx 9\,\text{m}`, matching the measured run-up closely.

.. image:: ../../res/tsunami_simulations/Soma_plot_1000m.png
   :align: center
   :width: 100%

Individual contributions
------------------------
Unfortunately, Magdalena Schwarzkopf was sick this week which is why the work was done by Dominik Münnch this week.