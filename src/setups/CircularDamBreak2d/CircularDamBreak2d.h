/**
 * @author Madalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Two-dimensional dam break problem.
 **/
#ifndef TSUNAMI_LAB_SETUPS_CIRCULAR_DAM_BREAK_2D_H
#define TSUNAMI_LAB_SETUPS_CIRCULAR_DAM_BREAK_2D_H

#include "../Setup.h"

namespace tsunami_lab {
  namespace setups {
    class CircularDamBreak2d;
  }
}

/**
 * 2d circular dam break setup.
 **/
class tsunami_lab::setups::CircularDamBreak2d: public Setup {
  private:
    //! height inside the dam
    t_real m_height_inside = 0;

    //! height outside the dam
    t_real m_height_outside = 0;

    //! momentum inside the dam in x-direction
    t_real m_momentumInsideX = 0;

    //! momentum inside the dam in y-direction
    t_real m_momentumInsideY = 0;

    //! momentum outside the dam in x-direction
    t_real m_momentumOutsideX = 0;

    //! momentum outside the dam in y-direction
    t_real m_momentumOutsideY = 0;

    //! location of the dam (radius)
    t_real m_middlePointDam_x = 0;
    t_real m_middlePointDam_y = 0;
    t_real m_radiusDam = 0;

  public:
    /**
     * Constructor.
     *
     * @param i_height_inside water height inside the dam.
     * @param i_height_outside water height outside the dam.
     * @param i_momentumInsideX momentum on the inside of the dam in x-direction.
     * @param i_momentumInsideY momentum on the inside of the dam in y-direction.
     * @param i_momentumOutsideX momentum on the outside of the dam in x-direction.
     * @param i_momentumOutsideY momentum on the outside of the dam in y-direction.
     * @param i_momentumRight momentum on the right side of the dam.
     * @param i_momentumUpwards momentum on the upwards side of the dam.
     * @param i_momentumDownwards momentum on the downwards side of the dam.
     * @param i_middlePointDam_x middle point (x-coordinate) of the dam.
     * @param i_middlePointDam_y middle point (y-coordinate) of the dam.
     * @param i_radiusDam radius of the dam.
     **/
    CircularDamBreak2d( t_real i_height_inside,
                        t_real i_height_outside,
                        t_real i_momentumInsideX,
                        t_real i_momentumInsideY,
                        t_real i_momentumOutsideX,
                        t_real i_momentumOutsideY,
                        t_real i_middlePointDam_x,
                        t_real i_middlePointDam_y,
                        t_real i_radiusDam);

    /**
     * Gets the water height at a given point.
     *
     * @param i_x x-coordinate of the queried point.
     * @return height at the given point.
     **/
    t_real getHeight( t_real i_x,
                      t_real i_y     ) const;

    /**
     * Gets the bathymetry at a given point.
     *
     * @return bathymetry at the given point.
     **/
    t_real getBathymetry( t_real i_x,
                          t_real i_y   ) const;

    /**
     * Gets the momentum in x-direction.
     *
     * @return momentum in x-direction.
     **/
    t_real getMomentumX( t_real i_x,
                         t_real i_y) const;

    /**
     * Gets the momentum in y-direction.
     *
     * @return momentum in y-direction.
     **/
    t_real getMomentumY( t_real i_x,
                         t_real i_y ) const;

};

#endif