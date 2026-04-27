/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * One-dimensional dam break problem.
 **/
#include "DamBreak1d.h"

tsunami_lab::setups::DamBreak1d::DamBreak1d( t_real i_heightLeft,
                                             t_real i_momentumLeft,
                                             t_real i_heightRight,
                                             t_real i_momentumRight,
                                             t_real i_locationDam) {
  m_heightLeft = i_heightLeft;
  m_momentumLeft = i_momentumLeft;
  m_heightRight = i_heightRight;
  m_momentumRight = i_momentumRight;
  m_locationDam = i_locationDam;
}

tsunami_lab::t_real tsunami_lab::setups::DamBreak1d::getHeight( t_real i_x,
                                                                t_real      ) const {
  if( i_x < m_locationDam ) {
    return m_heightLeft;
  }
  else {
    return m_heightRight;
  }
}

tsunami_lab::t_real tsunami_lab::setups::DamBreak1d::getBathymetry( t_idx,
                                                                    t_idx ) const {
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::DamBreak1d::getMomentumX( t_real i_x,
                                                                   t_real ) const {
  if (i_x < m_locationDam ) {
    return m_momentumLeft;
  }
  else {
    return m_momentumRight;
  }
}

tsunami_lab::t_real tsunami_lab::setups::DamBreak1d::getMomentumY( t_real,
                                                                   t_real ) const {
  return 0;
}