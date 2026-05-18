/**
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Two-dimensional Tsunami Event setup with netCDF input.
 **/
#include "TsunamiEvent2d.h"
#include <cmath>

tsunami_lab::setups::TsunamiEvent2d::TsunamiEvent2d( std::string const & i_bathymetryFile,
                                                      std::string const & i_displacementFile )
	: m_data( tsunami_lab::io::NetCdf::read( i_bathymetryFile, i_displacementFile ) ) {
}

tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent2d::getHeight( t_real i_x,
                                                                    t_real i_y ) const {
	// Read bathymetry from netCDF file
	t_real l_bIn = m_data.getBathymetry( i_x, i_y );

	if( l_bIn < 0 ) {
		return std::max( -l_bIn, m_delta );
	}
	else {
		return 0;
	}
}

tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent2d::getBathymetry( t_real i_x,
                                                                        t_real i_y ) const {
	// Read bathymetry and displacement from netCDF files
	t_real l_bIn = m_data.getBathymetry( i_x, i_y );
	t_real l_d = m_data.getDisplacement( i_x, i_y );

	if( l_bIn < 0 ) {
		return std::min( l_bIn, -m_delta ) + l_d;
	}
	else {
		return std::max( l_bIn, m_delta ) + l_d;
	}
}

tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent2d::getMomentumX( t_real,
                                                                       t_real ) const {
	return 0;
}

tsunami_lab::t_real tsunami_lab::setups::TsunamiEvent2d::getMomentumY( t_real,
                                                                       t_real ) const {
	return 0;
}
