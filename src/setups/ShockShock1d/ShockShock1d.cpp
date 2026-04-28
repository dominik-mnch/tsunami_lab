/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * One-dimensional Rare Rare problem.
 **/
#include "ShockShock1d.h"
#include <cmath>

tsunami_lab::setups::ShockShock1d::ShockShock1d( t_real i_height,
                                                 t_real i_momentum,
                                                 t_real i_locationDiscontinuity,
                                                 bool i_useBathymetryParabola,
                                                 t_real i_bathymetryOffset,
                                                 t_real i_bathymetryCenter,
                                                 t_real i_bathymetryScale ) {
  m_height = i_height;
  m_momentum = i_momentum;
  m_locationDiscontinuity = i_locationDiscontinuity;
  m_useBathymetryParabola = i_useBathymetryParabola;
  m_bathymetryOffset = i_bathymetryOffset;
  m_bathymetryCenter = i_bathymetryCenter;
  m_bathymetryScale = i_bathymetryScale;
}

tsunami_lab::t_real tsunami_lab::setups::ShockShock1d::getHeight( t_real,
                                                                t_real      ) const {
  //height is constant in this setup
  return m_height;
}

tsunami_lab::t_real tsunami_lab::setups::ShockShock1d::getBathymetry( t_real i_x,
                                                                      t_real ) const {
  if( m_useBathymetryParabola ) {
    if( m_bathymetryScale != 0 ) {
      t_real l_rootSq = -m_bathymetryOffset / m_bathymetryScale;

      if( l_rootSq > 0 ) {
        t_real l_halfWidth = std::sqrt( l_rootSq );

        if( std::abs( i_x - m_bathymetryCenter ) <= l_halfWidth ) {
          return m_bathymetryOffset + m_bathymetryScale * std::pow( i_x - m_bathymetryCenter, 2 );
        }
      }
    }

    return 0;
  }

  return 0;
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