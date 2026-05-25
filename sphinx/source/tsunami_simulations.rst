Tsunami Simulations
===================


Tohoku Event
------------

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
travel time of **28 minutes** and a maximum run-up water height of ~9.3 m.

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

This single-speed estimate is close to the measured 28 minutes.

Using the epicenter depth gives good agreement because it captures the
open-ocean wave speed that governs most of the travel time in reality.  In
the actual event the tsunami propagates two-dimensionally; energy arriving via
deeper offshore routes reaches the coast much faster than a shallow-shelf path
would suggest.

Domain and simulation setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~

First wave leaves the computational domain after 210 frames (30 fps). Since
each time step is about 10.1457 seconds of simulated time, the first wave
leaves the computational domain after about 2130.6 seconds (~35.5 minutes).

The domain spans −200 km to 2500 km in :math:`x` and −750 km to 750 km in
:math:`y`.  Grid sizes and cell-update counts:

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Resolution
     - Grid size
     - Cell updates per time step
   * - 1000 m
     - 2700 × 1500
     - :math:`4 \times 2700 \times 1500 + 2 \times 2700 + 2 \times 1500 = 16\,208\,400`
   * - 500 m
     - 5400 × 3000
     - :math:`4 \times 5400 \times 3000 + 2 \times 5400 + 2 \times 3000 = 64\,816\,800`

Simulation commands:

.. code-block:: bash

   # 1000 m resolution, outflow boundaries, 7200 s
   ./build/tsunami_lab 1000m -200000 2500000 -750000 750000 1 7200 2d outflow tsunami2d \
       data/output/tohoku_gebco20_ucsb3_250m_bath.nc \
       data/output/tohoku_gebco20_ucsb3_250m_displ.nc

   # 500 m resolution, outflow boundaries, 10800 s
   ./build/tsunami_lab 500m -200000 2500000 -750000 750000 1 10800 2d outflow tsunami2d \
       data/output/tohoku_gebco20_ucsb3_250m_bath.nc \
       data/output/tohoku_gebco20_ucsb3_250m_displ.nc