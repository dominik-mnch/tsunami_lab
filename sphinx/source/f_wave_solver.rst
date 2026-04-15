F-Wave Solver
=============

Overview
--------

The F-Wave Solver is a numerical solver for the one-dimensional shallow water equations. It implements the F-Wave method, which is based on computing wave speeds (eigenvalues) and wave strengths to compute net updates for fluid height and momentum.

**Authors:** Magdalena Schwarzkopf, Dominik Münch

Mathematical Background
-----------------------

The F-Wave solver addresses the shallow water equations system:

.. math::

    \frac{\partial h}{\partial t} + \frac{\partial hu}{\partial x} = 0
    
    \frac{\partial hu}{\partial t} + \frac{\partial}{\partial x}(hu^2 + \frac{1}{2}gh^2) = 0

where:

- :math:`h` is the height of the water surface
- :math:`u` is the horizontal velocity
- :math:`g` is the gravitational acceleration (approximated as :math:`\sqrt{g} \approx 3.1316`)

The solver uses a Riemann solver approach to compute net updates at cell interfaces based on left and right states.

Key Components
--------------

Wave Speeds (Eigenvalues)
~~~~~~~~~~~~~~~~~~~~~~~~~
This corresponds to the `eigenvalues` method in the `F_wave.cpp` source file.

**Algorithm:**

The method computes wave speeds using Roe averages:

1. Calculate square roots of left and right heights
2. Compute Roe-averaged height as :math:`\bar{h} = \frac{1}{2}\cdot(h_L + h_R)`
3. Compute Roe-averaged velocity: :math:`\bar{u} = \frac{\sqrt{h_L}u_L + \sqrt{h_R}u_R}{\sqrt{h_L} + \sqrt{h_R}}`
4. Compute wave speeds:
   
   - :math:`s_L = \bar{u} - \sqrt{g\bar{h}}`
   - :math:`s_R = \bar{u} + \sqrt{g\bar{h}}`

Wave Strengths
~~~~~~~~~~~~~~~
This corresponds to the `waveStrengths` method in the `F_wave.cpp` source file.

**Algorithm:**

The method computes wave strengths by solving the linearized system:

1. Compute the inverse of the right eigenvector matrix :math:`R`
2. Compute the flux jump:

.. math::

   \Delta f = f(q_R) - f(q_L), \quad
   q = \begin{bmatrix} h \\ hu \end{bmatrix}, \quad
   f(q) = \begin{bmatrix} hu \\ hu^2 + \frac{1}{2} g h^2 \end{bmatrix}

This is one of the differences between the F-wave and the Roe solver. The Roe uses the jump in quantities directly (:math:`q_R - q_L`) whereas the F-wave solver uses the jump in flux.

3. Solve :math:`\alpha = R^{-1} \Delta f` to get wave strengths

Net Updates
~~~~~~~~~~~
This corresponds to the `netUpdates` method in the `F_wave.cpp` source file.

**Algorithm:**

The method follows these steps:

1. Compute particle velocities from height and momentum: :math:`u = \frac{hu}{h}`
2. Call ``eigenvalues`` to compute wave speeds
3. Call ``waveStrengths`` to compute wave amplitudes
4. Compute net updates based on wave contributions:
   
   - Positive waves contribute to the right cell
   - Negative waves contribute to the left cell
   
This is another difference to the Roe solver. Here we add up the components of multiple different waves (if needed, see supersonic problem). This does not happen in the Roe solver.

5. Return updates scaled by wave speeds and strengths


How to build and test the code
------------------------------

The project is built using SCons:

.. code-block:: console
    
    scons

After building, you can run the tests using:

.. code-block:: console

   ./build/tests

Testing
-------

The F-Wave solver includes comprehensive unit tests in ``F_wave.test.cpp`` that verify:

- Correct computation of eigenvalues, (test case: "Eigenvalues")
- Accurate wave strength calculations, (test case: "RoeStrengths")
- net update computation: standard, supersonic problem, steady state, (test case: "FWaveUpdates")

Individual Contributions
-------------------------

This week we did every part of the project as a team, since it was the first assignment. 
