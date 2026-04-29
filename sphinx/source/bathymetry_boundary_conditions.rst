Bathymetry and Boundary Conditions
==================================

**Authors:** Magdalena Schwarzkopf, Dominik Münch

Adding bathymetry
-----------------

Previously, we have assumed that the bathymetry (the underwater topography) is flat and has no effect on the wave propagation.
This is obviously not realistic, so we need a way to include the bathymetry in our model.

To do this we modified the F-Wave Solver to include the bathymetry in the flux calculations.
More specifically, when calculating the jump in flux :math:`\Delta F`, we have to subtract the bathymetry term from the flux jump,
which is given by the following formula:

:math::

    \begin{split}\Delta x \Psi_{i-1/2} := \begin{bmatrix}
                                0 \\
                                -g (b_r - b_l) \frac{h_l+h_r}{2}
                              \end{bmatrix}.\end{split}

This was implemented in the F-Wave solver like this:

..code-block:: cpp

    // compute jump in flux
    t_real l_deltaF[2] = {0};
    l_deltaF[0] = i_huR - i_huL;
    l_deltaF[1] = (std::pow(i_huR,2) / i_hR  + 0.5 * std::pow(m_gSqrt,2) * std::pow(i_hR,2))
                    - (std::pow(i_huL,2) / i_hL  + 0.5 * std::pow(m_gSqrt,2) * std::pow(i_hL,2));

    //respect bathymetry
    t_real bathymetryTerm = -std::pow(m_gSqrt, 2)*(i_bR - i_bL) * ((i_hL + i_hR) / 2);
    l_deltaF[1] = l_deltaF[1] - bathymetryTerm;

Since this was previously unaccounted for, a lot of other classes had to be modified to include the bathymetry as well:

``WavePropagation.cpp`` now also has an array for the bathymetry which it passes on to the F-Wave solver when calculating the net updates.
Additionally, when setting the ghostOutflow it also sets this for the bathymetry.

``DamBreak1d.cpp`` now also has an array for the bathymetry and sets it to 0 for the whole domain.

``RareRare1d.cpp`` now also has an array for the bathymetry and sets it to 0 for the whole domain.
There is an additional option to set the bathymetry to a parabola which is controlled through additional parameters passed into the constructor:

..code-block:: cpp

    tsunami_lab::setups::RareRare1d::RareRare1d( t_real i_height,
                                                t_real i_momentum,
                                                t_real i_locationDiscontinuity,
                                                bool i_useBathymetryParabola,
                                                t_real i_bathymetryOffset,
                                                t_real i_bathymetryCenter,
                                                t_real i_bathymetryScale ) {
        m_height = i_height;
        m_momentum = i_momentum;
        m_locationDiscontinuity = i_locationDiscontinuity;
        m_useBathymetryParabola = i_useBathymetryParabola;
        m_bathymetryOffset = i_bathymetryOffset;
        m_bathymetryCenter = i_bathymetryCenter;
        m_bathymetryScale = i_bathymetryScale;
    }

The parabolic bathymetry is defined mathematically as:

.. math::

    b(x) = \begin{cases}
    b_{\text{offset}} + b_{\text{scale}} \left( x - x_c \right)^2 & \text{if } |x - x_c| \leq w \\
    0 & \text{otherwise}
    \end{cases}

where the half-width :math:`w` is determined by the constraint that the parabola touches zero at its extent:

.. math::

    w = \sqrt{-\frac{b_{\text{offset}}}{b_{\text{scale}}}}

The three controlling parameters are:

- :math:`b_{\text{offset}}`: vertical offset of the parabola (must be non-positive for the parabola to have real roots)
- :math:`b_{\text{scale}}`: quadratic coefficient controlling the curvature (must be positive for a valley-shaped parabola)
- :math:`x_c`: x-coordinate of the parabola's center (axis of symmetry)

This implementation ensures that the bathymetry smoothly transitions from the parabolic profile to zero bathymetry outside the parabola's extent.

..code-block:: cpp

    tsunami_lab::t_real tsunami_lab::setups::RareRare1d::getBathymetry( t_real i_x,
                                                                    t_real ) const {
        if( m_useBathymetryParabola ) {
            if( m_bathymetryScale != 0 ) {
            t_real l_rootSq = -m_bathymetryOffset / m_bathymetryScale;

            if( l_rootSq > 0 ) {
                t_real l_halfWidth = std::sqrt( l_rootSq );

                if( std::abs( i_x - m_bathymetryCenter ) <= l_halfWidth ) {
                return m_bathymetryOffset + m_bathymetryScale * std::pow( i_x - m_bathymetryCenter, 2 );
                }
            }
            }

            return 0;
        }

        return 0;
    }

``ShockShock1d.cpp`` now also has an array for the bathymetry and sets it to 0 for the whole domain.
It also provides the additional option to set the bathymetry to a parabola which works identically to the one in the Rare Rare setup.

To illustrate the effect of the bathymetry, we can look at the following simulation of the Shock Shock setup with a parabolic bathymetry:

- Shock Shock with parameters :math:`(h = 0.5, m = 0.3, location_discontinuity = 50, b_{\text{offset}} = 0.25, b_{\text{center}} = 50.0, b_{\text{scale}} = -0.08)`:

.. image:: ../../scripts/bathymetry_effect.gif
   :align: center

We can clearly see that instead of the wave in the middle having equal height everywhere, there are two peaks on either side.
This would not have happened without the bathymetry.

Adding boundary conditions
--------------------------

Until now, at the edges of the domain we have been using ghost cells that simply copy the values of the adjacent inner cells, which corresponds to an outflow boundary condition.
We now also want to be able to implement a reflective boundary condition, which corresponds to a wall at the edge of the domain where the fluid cannot flow through.
This means that the momentum at the ghost cell should be the negative of the adjacent inner cell, while the height should be the same as the adjacent inner cell.

To do this we implemented an enum in the ``WavePropagation1d.h`` header file to specify the type of boundary condition we want to use:

..code-block:: cpp

    //! ghost outflow boundary mode
        enum class BoundaryCondition {
        GhostOutflow,
        BoundaryRight,
        BoundaryLeft,
        BoundaryBoth
        };

This defines the 4 possible ways to set the boundary conditions. According to the specified mode, the ghost cells on the left and right edge
are set accordingly in the ``setGhostCells`` method of the ``WavePropagation1d.cpp`` file:

.. code-block: cpp

    void tsunami_lab::patches::WavePropagation1d::setGhostOutflow() {
        t_real * l_h = m_h[m_step];
        t_real * l_hu = m_hu[m_step];

        if (m_boundaryCondition == BoundaryCondition::GhostOutflow) {
            // set left boundary
            l_h[0] = l_h[1];
            l_hu[0] = l_hu[1];
            m_b[0] = m_b[1];

            // set right boundary
            l_h[m_nCells+1] = l_h[m_nCells];
            l_hu[m_nCells+1] = l_hu[m_nCells];
            m_b[m_nCells+1] = m_b[m_nCells];
        } else if (m_boundaryCondition == BoundaryCondition::BoundaryRight) {
            //set left boundary (ghost outflow)
            l_h[0] = l_h[1];
            l_hu[0] = l_hu[1];
            m_b[0] = m_b[1];

            // set right boundary (reflecting boundary)
            l_h[m_nCells+1] = l_h[m_nCells];
            l_hu[m_nCells+1] = - l_hu[m_nCells];
            m_b[m_nCells+1] = m_b[m_nCells];
        } else if (m_boundaryCondition == BoundaryCondition::BoundaryLeft) {
            // set left boundary (reflecting boundary)
            l_h[0] = l_h[1];
            l_hu[0] = - l_hu[1];
            m_b[0] = m_b[1];

            // set right boundary (ghost outflow)
            l_h[m_nCells+1] = l_h[m_nCells];
            l_hu[m_nCells+1] = l_hu[m_nCells];
            m_b[m_nCells+1] = m_b[m_nCells];
        } else if (m_boundaryCondition == BoundaryCondition::BoundaryBoth) {
            // set left boundary (reflecting boundary)
            l_h[0] = l_h[1];
            l_hu[0] = - l_hu[1];
            m_b[0] = m_b[1];

            // set right boundary (reflecting boundary)
            l_h[m_nCells+1] = l_h[m_nCells];
            l_hu[m_nCells+1] = - l_hu[m_nCells];
            m_b[m_nCells+1] = m_b[m_nCells];
        } 
    }

We can see the effect of the reflecting boundary conditions in the following Dam Break simulation with a reflecting boundary on the right edge.

- Dam Break with parameters :math:`(h_L = 5, m_L = 0, h_R = 2, location_Dam = 50)`:

.. image:: ../../scripts/reflecting_boundary_conditions.gif
   :align: center

When the shock wave comes in contact with the right wall, it is reflected back and starts traveling to the left.
This is equivalent to the left side of a shock shock setup.


3.3.1
- changes: main(after dambreak()), setup header, setup cpp, supercritical and subcritical get header and cpp setups, bathymetry now works with t_real
- location 3.3.1: 10.01
- Froude number 3.3.1: 0.584356
- location 3.3.2: 10.01
- Froude number 3.3.2: 1.22602
./build/tsunami_lab 1000 25 1 20

3.3 
[nix-shell:~/Project/tsunami_lab]$ ./build/tsunami_lab 1000 25 1 200
Hydraulic jump is where flow switches from Fr > 1 to Fr < 1, so where Froude number crosses 1
Question: should I be looking for the strongest jump? Or at each jump at every timestep?
Jump detected at x = 11.475 at time t = 3.82103
19.9735,11.4625,0.0125,0.0600494,-0.235125,0.124749
19.9735,11.4875,0.0125,0.146527,-0.238781,0.153369
19.9735,11.5125,0.0125,0.208031,-0.2425,0.124748
Die anderen Werte unterscheiden sich nur in der sechsten Nachkommastelle, diese findet sich in solutions 115