Large Data Input and Output
===========================

**Authors:** Magdalena Schwarzkopf, Dominik Münch

Overview
--------
In the previous week, we extended our code to support two-dimensional simulations.
This week, we added support for reading and writing large-scale data using netCDF files.


COARDS Convention
-----------------

The netCDF files follow the COARDS conventions and contain:

* coordinate variables ``x`` and ``y``
* a data variable ``z``

The bathymetry file represents:

.. math::

   b_{in}(x,y)

The displacement file represents:

.. math::

   d(x,y)

Both are stored on a regular grid.


TsunamiEvent2d
--------------

The ``TsunamiEvent2d`` setup initializes the simulation using bathymetry and
displacement data loaded from netCDF files.

Initial conditions are:

Water height:

.. math::

   h =
   \begin{cases}
   \max(-b_{in}, \delta), & b_{in} < 0 \\
   0, & \text{otherwise}
   \end{cases}

Momentum:

.. math::

   hu = 0, \quad hv = 0

Bathymetry:

.. math::

   b =
   \begin{cases}
   \min(b_{in}, -\delta) + d, & b_{in} < 0 \\
   \max(b_{in}, \delta) + d, & \text{otherwise}
   \end{cases}


Nearest-Neighbor Sampling
-------------------------

Simulation coordinates usually do not match the grid points in the netCDF files.
Therefore, values are obtained using nearest-neighbor lookup:

* closest x-grid point is selected
* closest y-grid point is selected
* corresponding value is returned

This may introduce small aliasing effects on coarse grids.


NetCDF Output
-------------

The ``NetCdf`` class writes simulation results to netCDF files.

It stores:

* bathymetry
* water height
* momentum in x and y direction
* simulation time

The output files can be visualized with tools such as ParaView.


Stored Variables
----------------

Coordinate variables:

* ``x`` — x-coordinate of cell centers
* ``y`` — y-coordinate of cell centers
* ``time`` — simulation time

Field variables:

* ``bathymetry``
* ``height``
* ``momentum_x``
* ``momentum_y``


Coordinate Generation
---------------------

Cell-center coordinates are computed as:

.. math::

   x_i = x_0 + \left(i + \frac{1}{2}\right)\Delta x

and similarly in y-direction.


Reading netCDF Data
-------------------

The ``read`` method:

#. opens the bathymetry file
#. reads grid dimensions and coordinates
#. loads bathymetry values
#. optionally loads displacement data
#. stores everything in a ``NetCdf::Data`` object

Data is stored in flattened arrays for efficiency.


Writing Time Steps
------------------

Each call to ``writeTimeStep`` appends one simulation step:

#. write simulation time
#. write water height
#. write momentum fields

Bathymetry is written only once since it is static.


Artificial Tsunami Setup
------------------------

The ``ArtificialTsunami2d`` setup is used for testing and validation of the
netCDF-based implementation.

It models a simple swimming pool with:

* constant bathymetry of ``-100 m``
* constant water height of ``100 m``
* zero initial momentum
* a localized displacement in the center

The displacement is defined as:

.. math::

   d(x,y) = 5 \cdot f(x) \cdot g(y)

with:

.. math::

   f(x) = \sin\left(\left(\frac{x}{500}+1\right)\pi\right), \quad
   g(y) = -\left(\frac{y}{500}\right)^2 + 1

for:

.. math::

   x,y \in [-500,500]

Outside this region, the displacement is zero.

The bathymetry is:

.. math::

   b(x,y) = -100 + d(x,y)

and the initial conditions are:

.. math::

   h(x,y) = 100, \quad hu = 0, \quad hv = 0


Input-Based vs Artificial Setup
-------------------------------

To validate the input-based ``TsunamiEvent2d`` setup, its results are compared
against the artificial setup. Since both describe the same physical scenario
(using different representations), their outputs approximately match.

Individual contributions
------------------------
Dominik Münch installed netCDF and implemented it's reading and writing functionality, as well as the
``TsunamiEvent2d`` setup. Magdalena Schwarzkopf implemented the ``ArtificialTsunami2d`` setup and 
performed the validation against the input-based setup. Both contributed to testing and documentation.