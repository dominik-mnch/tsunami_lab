/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Interface for NetCDF.
 **/

#include "NetCdf.h"
#include <filesystem>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
	constexpr tsunami_lab::t_idx c_invalidId = std::numeric_limits<tsunami_lab::t_idx>::max();

	// convert from t_idx to regular int for id
	int toNcId( tsunami_lab::t_idx i_id ) {
		if( i_id == c_invalidId ) return -1;
		return static_cast<int>( i_id );
	}

	// convert from regular int to t_idx for id
	tsunami_lab::t_idx fromNcId( int i_id ) {
		if( i_id < 0 ) return c_invalidId;
		return static_cast<tsunami_lab::t_idx>( i_id );
	}

	// check if the current operation was succesful, if not throw an error with context information
	void checkNc( int i_status, std::string const & i_context ) {
		if( i_status != NC_NOERR ) {
			std::ostringstream l_stream;
			l_stream << i_context << ": " << nc_strerror( i_status );
			throw std::runtime_error( l_stream.str() );
		}
	}

	// helper function to add text to a variable
	void putAttText(int i_ncId,
					int i_varId,
					char const * i_attName,
					std::string const & i_attValue ) {
		checkNc( nc_put_att_text(   i_ncId,
									i_varId,
									i_attName,
									i_attValue.size(),
									i_attValue.c_str() ),
						 			"nc_put_att_text" );
	}

	// helper function to compute and write the coordinates of the cell centers
	void putCoordinates(int i_ncId,
						int i_varId,
						tsunami_lab::t_idx i_n,
						tsunami_lab::t_real i_origin,
						tsunami_lab::t_real i_delta) {

		std::vector<tsunami_lab::t_real> l_coords( i_n );
		for( tsunami_lab::t_idx l_id = 0; l_id < i_n; l_id++ ) {
			l_coords[l_id] = i_origin +
											 ( static_cast<tsunami_lab::t_real>( l_id ) +
												 static_cast<tsunami_lab::t_real>( 0.5 ) ) * i_delta;
		}
		checkNc( nc_put_var_float( i_ncId, i_varId, l_coords.data() ), "nc_put_var_float" );
	}
}

tsunami_lab::io::NetCdf::NetCdf(t_real i_dx,
								t_real i_dy,
								t_real i_originX,
								t_real i_originY,
								t_idx i_nx,
								t_idx i_ny,
								t_idx i_stride,
								t_real i_time,
								t_real const * i_h,
								t_real const * i_b,
								t_real const * i_hu,
								t_real const * i_hv,
								std::string const & i_filePath )
	: m_dx( i_dx ),
		m_dy( i_dy ),
		m_originX( i_originX ),
		m_originY( i_originY ),
		m_nx( i_nx ),
		m_ny( i_ny ),
		m_stride( i_stride ),
		m_time( i_time ),
		m_h( i_h ),
		m_b( i_b ),
		m_hu( i_hu ),
		m_hv( i_hv ),
		m_timeStep( 0 ),
		m_ncId( c_invalidId ),
		m_varTimeId( c_invalidId ),
		m_varBathyId( c_invalidId ),
		m_varHeightId( c_invalidId ),
		m_varMomentumXId( c_invalidId ),
		m_varMomentumYId( c_invalidId ) {
	std::filesystem::path l_path( i_filePath );
	if( l_path.has_parent_path() ) {
		std::filesystem::create_directories( l_path.parent_path() );
	}

	bool l_isNewFile = !std::filesystem::exists( l_path );

	if( l_isNewFile ) {
		int l_ncId = -1;
		checkNc( nc_create( i_filePath.c_str(), NC_NETCDF4 | NC_CLOBBER, &l_ncId ), "nc_create" );
		m_ncId = fromNcId( l_ncId );

		// index of time dimension
		int l_dimTimeId = -1;
		// index of y dimension
		int l_dimYId = -1;
		// index of x dimension
		int l_dimXId = -1;
		// list of indices for y and x dimensions
		int l_dimsYX[2] = { -1, -1 };
		// list of indices for time, y and x dimensions
		int l_dimsTYX[3] = { -1, -1, -1 };

		// define dimensions (time (unlimited), y (size of y-direction), x (size of x-direction))
		checkNc( nc_def_dim( m_ncId, "time", NC_UNLIMITED, &l_dimTimeId ), "nc_def_dim(time)" );
		checkNc( nc_def_dim( m_ncId, "y", m_ny, &l_dimYId ), "nc_def_dim(y)" );
		checkNc( nc_def_dim( m_ncId, "x", m_nx, &l_dimXId ), "nc_def_dim(x)" );

		l_dimsYX[0] = l_dimYId;
		l_dimsYX[1] = l_dimXId;

		l_dimsTYX[0] = l_dimTimeId;
		l_dimsTYX[1] = l_dimYId;
		l_dimsTYX[2] = l_dimXId;

		int l_varXId = -1;
		int l_varYId = -1;
		int l_varTimeId = -1;
		// define variables for time and coordinates
		checkNc( nc_def_var( toNcId( m_ncId ), "time", NC_FLOAT, 1, &l_dimTimeId, &l_varTimeId ), "nc_def_var(time)" );
		m_varTimeId = fromNcId( l_varTimeId );
		checkNc( nc_def_var( toNcId( m_ncId ), "y", NC_FLOAT, 1, &l_dimYId, &l_varYId ), "nc_def_var(y)" );
		checkNc( nc_def_var( toNcId( m_ncId ), "x", NC_FLOAT, 1, &l_dimXId, &l_varXId ), "nc_def_var(x)" );

		// add text attributes to time and coordinate variables for paraview
		putAttText( toNcId( m_ncId ), l_varTimeId, "long_name", "time" );
		// paraview needs seconds since for time dimension
		putAttText( toNcId( m_ncId ), l_varTimeId, "units", "seconds since the earthquake event" );
		putAttText( toNcId( m_ncId ), l_varYId, "long_name", "y-coordinate of cell center" );
		putAttText( toNcId( m_ncId ), l_varYId, "units", "m" );
		putAttText( toNcId( m_ncId ), l_varXId, "long_name", "x-coordinate of cell center" );
		putAttText( toNcId( m_ncId ), l_varXId, "units", "m" );

		// define bathymetry variable (needs dimensions y and x)
		if( m_b != nullptr ) {
			int l_varBathyId = -1;
			checkNc( nc_def_var( toNcId( m_ncId ), "bathymetry", NC_FLOAT, 2, l_dimsYX, &l_varBathyId ), "nc_def_var(bathymetry)" );
			m_varBathyId = fromNcId( l_varBathyId );
			putAttText( toNcId( m_ncId ), l_varBathyId, "long_name", "sea floor depth" );
			putAttText( toNcId( m_ncId ), l_varBathyId, "units", "m" );
		}

		// define height variable (needs dimensions time, y and x)
		if( m_h != nullptr ) {
			int l_varHeightId = -1;
			checkNc( nc_def_var( toNcId( m_ncId ), "height", NC_FLOAT, 3, l_dimsTYX, &l_varHeightId ), "nc_def_var(height)" );
			m_varHeightId = fromNcId( l_varHeightId );
			putAttText( toNcId( m_ncId ), l_varHeightId, "long_name", "water height" );
			putAttText( toNcId( m_ncId ), l_varHeightId, "units", "m" );
		}

		// define momentum x variable (needs dimensions time, y and x)
		if( m_hu != nullptr ) {
			int l_varMomentumXId = -1;
			checkNc( nc_def_var( toNcId( m_ncId ), "momentum_x", NC_FLOAT, 3, l_dimsTYX, &l_varMomentumXId ), "nc_def_var(momentum_x)" );
			m_varMomentumXId = fromNcId( l_varMomentumXId );
			putAttText( toNcId( m_ncId ), l_varMomentumXId, "long_name", "x momentum" );
			putAttText( toNcId( m_ncId ), l_varMomentumXId, "units", "m^2 s^-1" );
		}

		// define momentum y variable (needs dimensions time, y and x)
		if( m_hv != nullptr ) {
			int l_varMomentumYId = -1;
			checkNc( nc_def_var( toNcId( m_ncId ), "momentum_y", NC_FLOAT, 3, l_dimsTYX, &l_varMomentumYId ), "nc_def_var(momentum_y)" );
			m_varMomentumYId = fromNcId( l_varMomentumYId );
			putAttText( toNcId( m_ncId ), l_varMomentumYId, "long_name", "y momentum" );
			putAttText( toNcId( m_ncId ), l_varMomentumYId, "units", "m^2 s^-1" );
		}

		putAttText( toNcId( m_ncId ), NC_GLOBAL, "Conventions", "COARDS" );
		checkNc( nc_enddef( toNcId( m_ncId ) ), "nc_enddef" );

		// write cell positions for x and y dimensions
		putCoordinates( toNcId( m_ncId ), l_varXId, m_nx, m_originX, m_dx );
		putCoordinates( toNcId( m_ncId ), l_varYId, m_ny, m_originY, m_dy );

		// write bathymetry data once as it does not change over time
		if( m_varBathyId != c_invalidId ) {
			std::vector<t_real> l_buffer( m_nx * m_ny );
			for( t_idx l_iy = 0; l_iy < m_ny; l_iy++ ) {
				for( t_idx l_ix = 0; l_ix < m_nx; l_ix++ ) {
					l_buffer[l_iy * m_nx + l_ix] = m_b[l_iy * m_stride + l_ix];
				}
			}
			checkNc( nc_put_var_float( toNcId( m_ncId ), toNcId( m_varBathyId ), l_buffer.data() ), "nc_put_var_float(bathymetry)" );
		}
	}
}

void tsunami_lab::io::NetCdf::writeTimeStep(t_real simTime) {
	// netCDF uses pointers for the start of the data and amount of data to write
	// write time dimension (start at current time step and write 1 value)
	std::size_t l_startT[1] = { m_timeStep };
	std::size_t l_countT[1] = { 1 };
	checkNc( nc_put_vara_float( toNcId( m_ncId ), toNcId( m_varTimeId ), l_startT, l_countT, &simTime ), "nc_put_vara_float(time)" );

	// buffer for data to be written to the netCDF file
	std::vector<t_real> l_buffer( m_nx * m_ny );

	// begin at current time step and write a block of size nx*ny
	std::size_t l_startData[3] = { m_timeStep, 0, 0 };
	std::size_t l_countData[3] = { 1, m_ny, m_nx };

	// write height data
	if( m_h != nullptr && m_varHeightId != c_invalidId ) {
		for( t_idx l_iy = 0; l_iy < m_ny; l_iy++ ) {
			for( t_idx l_ix = 0; l_ix < m_nx; l_ix++ ) {
				l_buffer[l_iy * m_nx + l_ix] = m_h[l_iy * m_stride + l_ix];
			}
		}
		checkNc( nc_put_vara_float( toNcId( m_ncId ), toNcId( m_varHeightId ), l_startData, l_countData, l_buffer.data() ), "nc_put_vara_float(height)" );
	}

	// write x momentum data
	if( m_hu != nullptr && m_varMomentumXId != c_invalidId ) {
		for( t_idx l_iy = 0; l_iy < m_ny; l_iy++ ) {
			for( t_idx l_ix = 0; l_ix < m_nx; l_ix++ ) {
				l_buffer[l_iy * m_nx + l_ix] = m_hu[l_iy * m_stride + l_ix];
			}
		}
		checkNc( nc_put_vara_float( toNcId( m_ncId ), toNcId( m_varMomentumXId ), l_startData, l_countData, l_buffer.data() ), "nc_put_vara_float(momentum_x)" );
	}

	// write y momentum data
	if( m_hv != nullptr && m_varMomentumYId != c_invalidId ) {
		for( t_idx l_iy = 0; l_iy < m_ny; l_iy++ ) {
			for( t_idx l_ix = 0; l_ix < m_nx; l_ix++ ) {
				l_buffer[l_iy * m_nx + l_ix] = m_hv[l_iy * m_stride + l_ix];
			}
		}
		checkNc( nc_put_vara_float( toNcId( m_ncId ), toNcId( m_varMomentumYId ), l_startData, l_countData, l_buffer.data() ), "nc_put_vara_float(momentum_y)" );
	}

	m_timeStep++;
}

tsunami_lab::io::NetCdf::~NetCdf() {
	if( m_ncId != c_invalidId ) {
		nc_close( toNcId( m_ncId ) );
		m_ncId = c_invalidId;
	}
}

