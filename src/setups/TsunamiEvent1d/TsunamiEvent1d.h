/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * One-dimensional Tsunami problem.
 **/
#ifndef TSUNAMI_LAB_SETUPS_TSUNAMI_EVENT_1D_H
#define TSUNAMI_LAB_SETUPS_TSUNAMI_EVENT_1D_H

#include "../Setup.h"
#include <vector>

namespace tsunami_lab {
  namespace setups {
    class TsunamiEvent1d;
  }
}

/**
 * 1d Tsunami Event setup.
 **/
class tsunami_lab::setups::TsunamiEvent1d: public Setup {

  private:

    //! bathymetry data
    std::vector<t_real> m_x;
    std::vector<t_real> m_b;

  public:

        /**
        * Constructor for TsunamiEvent1d setup.
        **/
    TsunamiEvent1d( );
    /**
     * Get's constant d(x).
     *
     * @param i_x x-coordinate.
     * @param i_y y-coordinate.
     * @return d(x).
     **/
    t_real getDisplacement( t_real i_x,
                      t_real i_y ) const;

    /**
     * Gets the water height at a given point.
     *
     * @param i_x x-coordinate.
     * @param i_y y-coordinate.
     * @return water height.
     **/
    t_real getHeight( t_real i_x,
                      t_real i_y ) const;

    /**
     * Gets the bathymetry at a given point.
     *
     * @param i_x x-coordinate.
     * @param i_y y-coordinate.
     * @return bathymetry.
     **/
    t_real getBathymetry( t_real i_x,
                          t_real i_y ) const;

    /**
     * Gets the momentum in x-direction at a given point.
     *
     * @param i_x x-coordinate.
     * @param i_y y-coordinate.
     * @return momentum in x-direction.
     **/
    t_real getMomentumX( t_real i_x,
                         t_real i_y ) const;

    /**
     * Gets the momentum in y-direction at a given point.
     *
     * @param i_x x-coordinate.
     * @param i_y y-coordinate.
     * @return momentum in y-direction.
     **/
    t_real getMomentumY( t_real i_x,
                         t_real i_y ) const;
};

#endif
