/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * One-dimensional Shock Shock problem.
 **/
#ifndef TSUNAMI_LAB_SETUPS_SHOCK_SHOCK_1D_H
#define TSUNAMI_LAB_SETUPS_SHOCK_SHOCK_1D_H

#include "../Setup.h"

namespace tsunami_lab {
  namespace setups {
    class ShockShock1d;
  }
}

/**
 * 1d Rare Rare setup.
 **/
class tsunami_lab::setups::ShockShock1d: public Setup {
  private:
    //! inital height of the water streams
    t_real m_height = 0;

    //! momentum of the water streams.
    t_real m_momentum = 0;

    //! location of the discontinuity
    t_real m_locationDiscontinuity = 0;

  public:
    /**
     * Constructor.
     *
     * @param i_height initial water height of the water streams.
     * @param i_momentum momentum of the water streams.
     * @param i_locationDiscontinuity location (x-coordinate) of the discontinuity.
     **/
    ShockShock1d( t_real i_height,
                t_real i_momentum,
                t_real i_locationDiscontinuity);

    /**
     * Gets the water height at a given point.
     *
     * @param i_x x-coordinate of the queried point.
     * @return height at the given point.
     **/
    t_real getHeight( t_real i_x,
                      t_real      ) const;

    /**
     * Gets the bathymetry at a given point.
     *
     * @return bathymetry at the given point.
     **/
    t_real getBathymetry( t_idx,
                          t_idx   ) const;

    /**
     * Gets the momentum in x-direction.
     *
     * @param i_x x-coordinate of the queried point.
     * @return momentum in x-direction.
     **/
    t_real getMomentumX( t_real i_x,
                         t_real ) const;

    /**
     * Gets the momentum in y-direction.
     *
     * @param i_x x-coordinate of the queried point.
     * @return momentum in y-direction.
     **/
    t_real getMomentumY( t_real i_x,
                         t_real ) const;

};

#endif