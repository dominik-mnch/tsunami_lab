/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * One-dimensional Tsunami problem.
 **/
#include "TsunamiEvent1d.h"
#include <cmath>
#include <algorithm>
#include "../../io/Csv.h"

#ifndef BATHYMETRY_CSV
#define BATHYMETRY_CSV "bathymetry.csv"
#endif

tsunami_lab::setups::TsunamiEvent1d::TsunamiEvent1d() {

    bool result = tsunami_lab::io::Csv::readBathymetry( BATHYMETRY_CSV, m_x, m_b );  

    if (!result ) {
        std::cerr << "Error: Failed to read bathymetry data for TsunamiEvent1d setup." << std::endl;
    }
}



tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent1d::getDisplacement( t_real i_x,
                                                      t_real ) const {
    if(175000 < i_x && i_x < 250000) {
      return 10 * std::sin((i_x - 175000) / 37500.0 * M_PI + M_PI);
    }
    else {
      return 0;
    }
}


tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent1d::getHeight( t_real i_x,
                                                                t_real      ) const {

  t_idx cell_index = (i_x / 440000) * m_b.size();                                                   
  if( m_b[cell_index] < 0 ) {
    return std::max( -1 * m_b[cell_index], static_cast<t_real>(20.0) );
  }
  else {
    return 0;
  }
}

tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent1d::getBathymetry( t_real i_x,
                                                                      t_real ) const {

  t_idx cell_index = (i_x / 440000) * m_b.size();

  if(m_b[cell_index] < 0) {
    return  std::min( m_b[cell_index], static_cast<t_real>(-20.0)) + getDisplacement(i_x, 0);
  }
  else {
    return  std::max( m_b[cell_index], static_cast<t_real>(20.0)) + getDisplacement(i_x, 0);
  }
}


tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent1d::getMomentumX( t_real,
                                                                 t_real ) const {
    return 0;
}

tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent1d::getMomentumY( t_real,
                                                                   t_real ) const {
  //y-direction is not relevant in the 1d case
  return 0;
}