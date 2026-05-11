/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Two-dimensional Tsunami problem.
 **/
#ifndef TSUNAMI_LAB_SETUPS_TSUNAMI_EVENT_2D_H
#define TSUNAMI_LAB_SETUPS_TSUNAMI_EVENT_2D_H

#include "../Setup.h"
#include <vector>

namespace tsunami_lab {
  namespace setups {
    class ArtificialTsunami2d;
  }
}

/**
 * 2d Tsunami Event setup.
 **/
class tsunami_lab::setups::ArtificialTsunami2d: public Setup {

  private:



  public:

    /**
     * Constructor for TsunamiEvent2d setup.
    **/
    ArtificialTsunami2d();
    
    /**
     * Gets water displacement at a given point.
     *
     * @param i_x x-coordinate.
     * @param i_y y-coordinate.
     * @return water displacement at the given point.
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
