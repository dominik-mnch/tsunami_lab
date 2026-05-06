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
  m_middlePointDam_x = i_middlePointDam_x;
  m_middlePointDam_y = i_middlePointDam_y;
  m_radiusDam = i_radiusDam;
      
}

tsunami_lab::t_real tsunami_lab::setups::CircularDamBreak2d::getHeight( t_real i_x,
                                                                        t_real i_y ) const {
  t_real l_xDistance = i_x - m_middlePointDam_x;
  t_real l_yDistance = i_y - m_middlePointDam_y;
  t_real l_distance = std::sqrt( l_xDistance * l_xDistance + l_yDistance * l_yDistance );

  if( l_distance < m_radiusDam ) {
    return m_height_inside;
  }

  return m_height_outside;
}

tsunami_lab::t_real tsunami_lab::setups::CircularDamBreak2d::getBathymetry( t_real i_x,
                                                                            t_real i_y ) const {
  if( sqrt( pow(i_x + 20, 2) + pow(i_y, 2) ) < 7 ) {
    return (pow((i_x/20), 2) + pow((i_y/20), 2));
  }
  else {
    return 0;
  }
}

tsunami_lab::t_real tsunami_lab::setups::CircularDamBreak2d::getMomentumX( t_real i_x,
                                                                           t_real i_y ) const {
  t_real l_xDistance = i_x - m_middlePointDam_x;
  t_real l_yDistance = i_y - m_middlePointDam_y;
  t_real l_distance = std::sqrt( l_xDistance * l_xDistance + l_yDistance * l_yDistance );

  if( l_distance < m_radiusDam ) {
    return m_momentumInsideX;
  }

  return m_momentumOutsideX;
}

tsunami_lab::t_real tsunami_lab::setups::CircularDamBreak2d::getMomentumY( t_real i_x,
                                                                           t_real i_y ) const {
  t_real l_xDistance = i_x - m_middlePointDam_x;
  t_real l_yDistance = i_y - m_middlePointDam_y;
  t_real l_distance = std::sqrt( l_xDistance * l_xDistance + l_yDistance * l_yDistance );

  if( l_distance < m_radiusDam ) {
    return m_momentumInsideY;
  }

  return m_momentumOutsideY;
}