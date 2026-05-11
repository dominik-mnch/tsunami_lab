/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Two-dimensional Tsunami problem.
 **/
#include "ArtificialTsunami2d.h"
#include <cmath>
#include <algorithm>
#include "../../io/Csv.h"



tsunami_lab::setups::ArtificialTsunami2d::ArtificialTsunami2d() {
}


tsunami_lab::t_real tsunami_lab::setups::ArtificialTsunami2d::getDisplacement( t_real i_x,
                                                                t_real i_y ) const {
    if(i_x < -500 || i_x > 500 || i_y < -500 || i_y > 500) {
        return 0;
    }
    else {
        t_real f_x = std::sin ( ((i_x/500) + 1) * M_PI);
        t_real f_y = -1 * std::pow((i_y/500),2) + 1;
        return 5 * f_x * f_y;
    }
}



tsunami_lab::t_real tsunami_lab::setups::ArtificialTsunami2d::getHeight( t_real ,
                                                                t_real  ) const {  
 return 100;

}

tsunami_lab::t_real tsunami_lab::setups::ArtificialTsunami2d::getBathymetry( t_real ,
                                                                      t_real ) const {
 return -100;
}


tsunami_lab::t_real tsunami_lab::setups::ArtificialTsunami2d::getMomentumX( t_real ,
                                                                 t_real ) const {
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::ArtificialTsunami2d::getMomentumY( t_real ,
                                                                   t_real  ) const {
  //y-direction is not relevant in the 1d case
  return 0;
}