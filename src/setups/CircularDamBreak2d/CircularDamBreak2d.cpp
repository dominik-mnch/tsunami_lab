/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Two-dimensional circular dam break problem.
 **/
#include "CircularDamBreak2d.h"

tsunami_lab::setups::CircularDamBreak2d::CircularDamBreak2d( t_real i_height_inside,
                                                             t_real i_height_outside,       
                                                             t_real i_momentumInsideX,
                                                             t_real i_momentumInsideY,
                                                             t_real i_momentumOutsideX,
                                                             t_real i_momentumOutsideY,
                                                             t_real i_middlePointDam_x,
                                                             t_real i_middlePointDam_y,
                                                             t_real i_radiusDam) {
  m_height_inside = i_height_inside;
  m_height_outside = i_height_outside;
  m_momentumInsideX = i_momentumInsideX;
  m_momentumInsideY = i_momentumInsideY;
  m_momentumOutsideX = i_momentumOutsideX;
  m_momentumOutsideY = i_momentumOutsideY;
  m_locationDam = i_middlePointDam_x;
  m_locationDam = i_middlePointDam_y;
  m_locationDam = i_radiusDam;
      
}

tsunami_lab::t_real tsunami_lab::setups::CircularDamBreak2d::getHeight( t_real i_x,
                                                                        t_real i_y ) const {
  if( sqrt( pow(i_x, 2) + pow(i_y, 2) ) < m_locationDam ) {
    return 10;
  }
  else {
    return 5;
  }
}

tsunami_lab::t_real tsunami_lab::setups::CircularDamBreak2d::getBathymetry( t_real,
                                                                            t_real ) const {
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::CircularDamBreak2d::getMomentumX( t_real,
                                                                           t_real ) const {
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::CircularDamBreak2d::getMomentumY( t_real,
                                                                           t_real ) const {
  return 0;
}