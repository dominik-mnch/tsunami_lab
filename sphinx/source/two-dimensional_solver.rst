Two-Dimensional Solver
======================

**Authors:** Magdalena Schwarzkopf, Dominik Münch

Overview
--------
So far, we have only been able to simulate one-dimensional settings. This week, we have extended our code to also support two-dimensional simulations.

The ``WavePropagation.h`` header file was always designed to support both one-dimensional and two-dimensional simulations.
The implementation ``WavePropagation1d.cpp`` simply chose to not support two-dimensional simulations by ignoring the y-direction entirely and making functions like ``getMomentumY`` return dummy values.

This week we have also improved our command line interface (in ``main.cpp``) to accept parameters to control the wave propagation mode (1d or 2d), the boundary conditions for 1d simulations and the setup along with their respective parameters.
A summary of the command line interface can be found in the ``README.rst`` file.

Two-Dimensional Wave Propagation
--------------------------------

The new implementation ``WavePropagation2d.cpp`` now fully supports two-dimensional simulations. To do this, we are using the unsplit method:

.. math::

   \begin{aligned}
   Q_{i,j}^{n+1} = Q_{i,j}^n 
   &- \frac{\Delta t}{\Delta x} \left( A^+ \Delta Q_{i-1/2,j} + A^- \Delta Q_{i+1/2,j} \right) \\
   &- \frac{\Delta t}{\Delta y} \left( B^+ \Delta Q_{i,j-1/2} + B^- \Delta Q_{i,j+1/2} \right) \\
   &\quad \forall i \in \{ 1, \dots, n \}, \; j \in \{ 0, \dots, n \}.
   \end{aligned}

This means that a cell now has 4 different updates which each of its neighbours, two in the x-direction and two in the y-direction.
We can calculate the net updates by applying our F-Wave solver (or Roe solver) to the edges between the cell and its neighbours and then applying the net updates to the cell itself.

We can see this in action in the following simulation of a dam break in a 2d domain:

.. video:: ../../res/Dam_Break_2d.mp4
   :align: center
   :width: 100%
   :autoplay:
   :loop:

This looks very similar to the 1d Dam Break simulation but it doesn't really add any new information since it's basically just a copy of the 1d simulation for each y-value.

To illustrate the 2d simulation when the y-direction actually has an effect, we implemented the ``CircularDamBreak2d.cpp`` setup which simulates a circular dam break in the middle of a :math:`100 \times 100` domain.

.. video:: ../../res/Circular_Dam_Break_2d_no_bath.mp4
   :align: center
   :width: 100%
   :autoplay:
   :loop:

We can now see the wave propagating in a circular manner in all directions. 

The above simulation also didn't have any bathymetry. To show the effects that bathymetry can have in two dimensions, we added bathymetry as follows:

.. math::

    b(x, y) = \begin{cases}
        \frac{x^2}{400} + \frac{y^2}{400} & \text{if } \sqrt{(x + 20)^2 + y^2} < 7 \\
        0 & \text{otherwise}
    \end{cases}

That means that there is now a small hill in the water. To show the effects this has on the simulation, we can take a look at the following simulation:

.. video:: ../../res/Circular_Dam_Break_2d_bath.mp4
   :align: center
   :width: 100%
   :autoplay:
   :loop:

As we can see there is now a small water valley right above the bathymetry. We can compare this to our simulation of the :ref:`subcritical setup <subcritical_setup>` where we can also see this exact phenomenon.
A small water valley forming right above the hill in the bathymetry.