#include "Setup.h"
#include "../constants.h"
#include <cmath>

float tsunami_lab::setups::Setup::getFroudeNumber( float i_x,
                              float i_y ) const {

    float h = getHeight(i_x,i_y);
    float hu = getMomentumX(i_x,i_y);
    
    if (h <= 0) return 0.0;

    float u = hu / h;

    return std::abs(u) / std::sqrt(9.81 * h);
}
