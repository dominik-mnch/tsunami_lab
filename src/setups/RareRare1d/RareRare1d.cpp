/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * One-dimensional Rare Rare problem.
 **/
#include "RareRare1d.h"

tsunami_lab::setups::RareRare1d::RareRare1d( t_real i_height,
                                             t_real i_momentum,
                                             t_real i_locationDiscontinuity ) {
  m_height = i_height;
  m_momentum = i_momentum;
  m_locationDiscontinuity = i_locationDiscontinuity;
}

tsunami_lab::t_real tsunami_lab::setups::RareRare1d::getHeight( t_real i_x,
                                                                t_real      ) const {
  //height is constant in this setup
  return m_height;
}



tsunami_lab::t_real tsunami_lab::setups::RareRare1d::getMomentumX( t_real i_x,
                                                                   t_real ) const {
  if(i_x < m_locationDiscontinuity){
    return -m_momentum; // Left side moves to the left
  }
  else{
    return m_momentum; // Right side moves to the right
  }
}

tsunami_lab::t_real tsunami_lab::setups::RareRare1d::getMomentumY( t_real,
                                                                   t_real ) const {
  //y-direction is not relevant in the 1d case
  return 0;
}