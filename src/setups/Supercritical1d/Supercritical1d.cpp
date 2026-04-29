/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * One-dimensional supercritical problem.
 **/
#include "Supercritical1d.h"
#include <cmath>

tsunami_lab::setups::Supercritical1d::Supercritical1d() {
}

tsunami_lab::t_real tsunami_lab::setups::Supercritical1d::getHeight( t_real i_x,
                                                                t_real      ) const {
    return -1 * getBathymetry( i_x, 0 );
}

tsunami_lab::t_real tsunami_lab::setups::Supercritical1d::getBathymetry( t_real i_x,
                                                                    t_real ) const {
  if( i_x > 8 && i_x < 12 ) {
    return -0.13 -0.05 * std::pow( i_x - 10, 2 );
  }
  else {
    return -0.33;
  }
}

tsunami_lab::t_real tsunami_lab::setups::Supercritical1d::getMomentumX( t_real ,
                                                                   t_real ) const {
    return 0.18;                                                                
}

tsunami_lab::t_real tsunami_lab::setups::Supercritical1d::getMomentumY( t_real,
                                                                   t_real ) const {
  return 0;
}