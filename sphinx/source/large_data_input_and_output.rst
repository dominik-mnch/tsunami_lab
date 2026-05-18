Two-Dimensional Solver
======================

**Authors:** Magdalena Schwarzkopf, Dominik Münch

Overview
--------
Last week, we extended our code to support two-dimensional simulations. This week, we have added support for large data input and output to our code.

artificial two-dimensional tsunami setup
----------------------------------------

The ``ArtificialTsunami2d`` setup implements an artificial tsunami scenario
used for testing two-dimensional tsunami simulations and validating the
netCDF-based ``TsunamiEvent2d`` setup.

The setup models a swimming pool with:

* constant bathymetry of ``-100 m``,
* constant water height of ``100 m``,
* zero initial momentum,
* a localized vertical displacement in the center of the domain.

The vertical displacement is defined as:

.. math::

   d(x,y) = 5 \cdot f(x) \cdot g(y)

with

.. math::

   f(x) = \sin\left(\left(\frac{x}{500}+1\right)\pi\right)

and

.. math::

   g(y) = -\left(\frac{y}{500}\right)^2 + 1

for:

.. math::

   x,y \in [-500,500]

Outside this region the displacement is zero.

The bathymetry is defined as:

.. math::

   b(x,y) = -100 + d(x,y)

The initial water height is constant:

.. math::

   h(x,y) = 100

The initial momenta are zero:

.. math::

   hu(x,y) = 0

.. math::

   hv(x,y) = 0



input-based Tsunami Setup versus the artificial Setup
-----------------------------------------------------
To check the 