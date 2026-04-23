/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * One-dimensional Rare Rare problem.
 **/
#include "ShockShock1d.h"

tsunami_lab::setups::ShockShock1d::ShockShock1d( t_real i_height,
                                             t_real i_momentum,
                                             t_real i_locationDiscontinuity,
                                             t_real * i_b ) {
  m_height = i_height;
  m_momentum = i_momentum;
  m_locationDiscontinuity = i_locationDiscontinuity;
  m_b = i_b;
}

tsunami_lab::t_real tsunami_lab::setups::ShockShock1d::getHeight( t_real,
                                                                t_real      ) const {
  //height is constant in this setup
  return m_height;
}

tsunami_lab::t_real tsunami_lab::setups::ShockShock1d::getBathymetry( t_idx i_ix,
                                                                      t_idx ) const {
  return m_b[i_ix];
}

tsunami_lab::t_real tsunami_lab::setups::ShockShock1d::getMomentumX( t_real i_x,
                                                                   t_real ) const {
  if(i_x < m_locationDiscontinuity){
    return m_momentum; // Left side moves to the right
  }
  else{
    return -m_momentum; // Right side moves to the left
  }
}

tsunami_lab::t_real tsunami_lab::setups::ShockShock1d::getMomentumY( t_real,
                                                                   t_real ) const {
  //y-direction is not relevant in the 1d case
  return 0;
}