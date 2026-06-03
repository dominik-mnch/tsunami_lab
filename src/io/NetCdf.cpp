/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Interface for NetCDF.
 **/

#include "NetCdf.h"
#include <algorithm>
#include <cmath>
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

void tsunami_lab::io::NetCdf::helperWritingData(t_real const * m_some, t_idx m_someId) const {
	// begin at current time step and write a block of size nx*ny
	std::size_t l_startData[3] = { m_timeStep, 0, 0 };
	std::size_t l_countData[3] = { 1, m_ny/m_k, m_nx/m_k };

	// buffer for data to be written to the netCDF file
	std::vector<t_real> l_buffer( (m_nx + m_k - 1) / m_k * (m_ny + m_k - 1) / m_k );
	
	if( m_some != nullptr && m_someId != c_invalidId ) {
		for( t_idx l_iy = 0; l_iy < m_ny; l_iy+= m_k ) {
			for( t_idx l_ix = 0; l_ix < m_nx; l_ix+= m_k ) {
				t_real momentumYAvg = 0;

				// Handle blocks that might be smaler at the edges
				t_idx blockHeight = std::min(m_k, m_ny - l_iy);
				t_idx blockWidth = std::min(m_k, m_nx - l_ix);

				for( t_idx l_j = 0; l_j < blockHeight; l_j++ ) {
					for( t_idx l_i = 0; l_i < blockWidth; l_i++ ) {
						momentumYAvg += m_some[(l_iy+l_j) * m_stride + (l_ix+l_i)];
					}
				}
				t_idx l_oy = l_iy / m_k;
				t_idx l_ox = l_ix / m_k;
				l_buffer[l_oy * (m_nx/m_k) + l_ox] = momentumYAvg / (blockHeight * blockWidth);
			}
		}
		checkNc( nc_put_vara_float( toNcId( m_ncId ), toNcId( m_varMomentumYId ), l_startData, l_countData, l_buffer.data() ), "nc_put_vara_float(momentum_y)" );
	}
}

tsunami_lab::io::NetCdf::NetCdf(t_real i_dx,
								t_real i_dy,
								t_real i_originX,
								t_real i_originY,
								t_idx i_nx,
								t_idx i_ny,
								t_idx i_k,
								t_idx i_stride,
								t_real i_time,
								t_real i_endTime,
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
		m_k( i_k ),
		m_stride( i_stride ),
		m_time( i_time ),
		m_endTime( i_endTime ),
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
		m_varMomentumYId( c_invalidId ),
		m_checkpointNcId( c_invalidId ),
		m_lastSimTime( 0 ) {
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
		checkNc( nc_def_dim( m_ncId, "y", m_ny/m_k, &l_dimYId ), "nc_def_dim(y)" );
		checkNc( nc_def_dim( m_ncId, "x", m_nx/m_k, &l_dimXId ), "nc_def_dim(x)" );

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

		// enable compression on time variable
		checkNc( nc_def_var_deflate( toNcId( m_ncId ), l_varTimeId, 1, 1, 4 ), "nc_def_var_deflate(time)" );

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
			// enable compression on height variable
			checkNc( nc_def_var_deflate( toNcId( m_ncId ), l_varHeightId, 1, 1, 4 ), "nc_def_var_deflate(height)" );
			putAttText( toNcId( m_ncId ), l_varHeightId, "long_name", "water height" );
			putAttText( toNcId( m_ncId ), l_varHeightId, "units", "m" );
		}

		// define momentum x variable (needs dimensions time, y and x)
		if( m_hu != nullptr ) {
			int l_varMomentumXId = -1;
			checkNc( nc_def_var( toNcId( m_ncId ), "momentum_x", NC_FLOAT, 3, l_dimsTYX, &l_varMomentumXId ), "nc_def_var(momentum_x)" );
			m_varMomentumXId = fromNcId( l_varMomentumXId );
			// enable compression on momentum_x variable
			checkNc( nc_def_var_deflate( toNcId( m_ncId ), l_varMomentumXId, 1, 1, 4 ), "nc_def_var_deflate(momentum_x)" );
			putAttText( toNcId( m_ncId ), l_varMomentumXId, "long_name", "x momentum" );
			putAttText( toNcId( m_ncId ), l_varMomentumXId, "units", "m^2 s^-1" );
		}

		// define momentum y variable (needs dimensions time, y and x)
		if( m_hv != nullptr ) {
			int l_varMomentumYId = -1;
			checkNc( nc_def_var( toNcId( m_ncId ), "momentum_y", NC_FLOAT, 3, l_dimsTYX, &l_varMomentumYId ), "nc_def_var(momentum_y)" );
			m_varMomentumYId = fromNcId( l_varMomentumYId );
			// enable compression on momentum_y variable
			checkNc( nc_def_var_deflate( toNcId( m_ncId ), l_varMomentumYId, 1, 1, 4 ), "nc_def_var_deflate(momentum_y)" );
			putAttText( toNcId( m_ncId ), l_varMomentumYId, "long_name", "y momentum" );
			putAttText( toNcId( m_ncId ), l_varMomentumYId, "units", "m^2 s^-1" );
		}

		putAttText( toNcId( m_ncId ), NC_GLOBAL, "Conventions", "COARDS" );
		checkNc( nc_enddef( toNcId( m_ncId ) ), "nc_enddef" );

		// write cell positions for x and y dimensions
		putCoordinates( toNcId( m_ncId ), l_varXId, m_nx/m_k, m_originX, m_dx );
		putCoordinates( toNcId( m_ncId ), l_varYId, m_ny/m_k, m_originY, m_dy );

		// write bathymetry data once as it does not change over time 
		helperWritingData(m_b, m_varBathyId);
	}
	else {
		// File already exists (checkpoint resume): open for appending and look up variable IDs.
		int l_ncId = -1;
		checkNc( nc_open( i_filePath.c_str(), NC_WRITE, &l_ncId ), "nc_open(existing)" );
		m_ncId = fromNcId( l_ncId );

		int l_varTimeId = -1;
		checkNc( nc_inq_varid( l_ncId, "time", &l_varTimeId ), "nc_inq_varid(time)" );
		m_varTimeId = fromNcId( l_varTimeId );

		int l_varBathyId = -1;
		if( m_b != nullptr && nc_inq_varid( l_ncId, "bathymetry", &l_varBathyId ) == NC_NOERR )
			m_varBathyId = fromNcId( l_varBathyId );

		int l_varHeightId = -1;
		if( m_h != nullptr && nc_inq_varid( l_ncId, "height", &l_varHeightId ) == NC_NOERR )
			m_varHeightId = fromNcId( l_varHeightId );

		int l_varMomentumXId = -1;
		if( m_hu != nullptr && nc_inq_varid( l_ncId, "momentum_x", &l_varMomentumXId ) == NC_NOERR )
			m_varMomentumXId = fromNcId( l_varMomentumXId );

		int l_varMomentumYId = -1;
		if( m_hv != nullptr && nc_inq_varid( l_ncId, "momentum_y", &l_varMomentumYId ) == NC_NOERR )
			m_varMomentumYId = fromNcId( l_varMomentumYId );

		// Set m_timeStep to the existing number of steps so new writes append correctly.
		int         l_dimTimeId = -1;
		std::size_t l_nSteps    = 0;
		if( nc_inq_dimid( l_ncId, "time", &l_dimTimeId ) == NC_NOERR )
			nc_inq_dimlen( l_ncId, l_dimTimeId, &l_nSteps );
		m_timeStep = static_cast<t_idx>( l_nSteps );
	}
	defineCheckpoint( (l_path.parent_path() / "checkpoint.nc").string() );
}

void tsunami_lab::io::NetCdf::writeTimeStep(t_real simTime) {
	// netCDF uses pointers for the start of the data and amount of data to write
	// write time dimension (start at current time step and write 1 value)
	std::size_t l_startT[1] = { m_timeStep };
	std::size_t l_countT[1] = { 1 };
	checkNc( nc_put_vara_float( toNcId( m_ncId ), toNcId( m_varTimeId ), l_startT, l_countT, &simTime ), "nc_put_vara_float(time)" );


	// write height data
	helperWritingData(m_h, m_varHeightId);


	// write x momentum data
	helperWritingData(m_hu, m_varMomentumXId);

	// write y momentum data
	helperWritingData(m_hv, m_varMomentumYId);

	m_lastSimTime = simTime;
	overwriteCheckpointSimTime();
	checkNc( nc_sync( toNcId( m_ncId ) ), "nc_sync(solution)" );
	m_timeStep++;
}

void tsunami_lab::io::NetCdf::defineCheckpoint( std::string const & i_filePath ) {
	std::filesystem::path l_path( i_filePath );
	if( l_path.has_parent_path() ) {
		std::filesystem::create_directories( l_path.parent_path() );
	}

	int l_ncId = -1;
	checkNc( nc_create( i_filePath.c_str(), NC_NETCDF4 | NC_CLOBBER, &l_ncId ), "nc_create(checkpoint)" );
	m_checkpointNcId = fromNcId( l_ncId );

	int   l_nx      = static_cast<int>( m_nx );
	int   l_ny      = static_cast<int>( m_ny );
	float l_dx      = static_cast<float>( m_dx );
	float l_dy      = static_cast<float>( m_dy );
	float l_originX = static_cast<float>( m_originX );
	float l_originY = static_cast<float>( m_originY );
	float l_endTime = static_cast<float>( m_endTime );
	float l_simTime = 0.0f;


	checkNc( nc_put_att_int(   l_ncId, NC_GLOBAL, "nx",       NC_INT,   1, &l_nx      ), "nc_put_att_int(nx)"       );
	checkNc( nc_put_att_int(   l_ncId, NC_GLOBAL, "ny",       NC_INT,   1, &l_ny      ), "nc_put_att_int(ny)"       );
	checkNc( nc_put_att_float( l_ncId, NC_GLOBAL, "dx",       NC_FLOAT, 1, &l_dx      ), "nc_put_att_float(dx)"     );
	checkNc( nc_put_att_float( l_ncId, NC_GLOBAL, "dy",       NC_FLOAT, 1, &l_dy      ), "nc_put_att_float(dy)"     );
	checkNc( nc_put_att_float( l_ncId, NC_GLOBAL, "origin_x", NC_FLOAT, 1, &l_originX ), "nc_put_att_float(origin_x)" );
	checkNc( nc_put_att_float( l_ncId, NC_GLOBAL, "origin_y", NC_FLOAT, 1, &l_originY ), "nc_put_att_float(origin_y)" );
	checkNc( nc_put_att_float( l_ncId, NC_GLOBAL, "end_time", NC_FLOAT, 1, &l_endTime ), "nc_put_att_float(end_time)" );
	checkNc( nc_put_att_float( l_ncId, NC_GLOBAL, "sim_time", NC_FLOAT, 1, &l_simTime ), "nc_put_att_float(sim_time)" );

	checkNc( nc_sync( l_ncId ), "nc_sync(checkpoint)" );
}

void tsunami_lab::io::NetCdf::overwriteCheckpointSimTime() {
	if( m_checkpointNcId == c_invalidId ) return;

	float l_simTimeF = static_cast<float>( m_lastSimTime );
	checkNc( nc_put_att_float( toNcId( m_checkpointNcId ), NC_GLOBAL, "sim_time", NC_FLOAT, 1, &l_simTimeF ),
	         "nc_put_att_float(sim_time)" );
	checkNc( nc_sync( toNcId( m_checkpointNcId ) ), "nc_sync(checkpoint)" );
}

tsunami_lab::io::NetCdf::CheckpointData tsunami_lab::io::NetCdf::readCheckpoint( std::string const & i_checkpointFile ) {
	if( !std::filesystem::exists( i_checkpointFile ) ) {
		throw std::invalid_argument( "checkpoint file does not exist: " + i_checkpointFile );
	}

	CheckpointData l_cp;

	// ── 1. Read simulation parameters from checkpoint.nc ──────────────────────
	int l_cpId = -1;
	checkNc( nc_open( i_checkpointFile.c_str(), NC_NOWRITE, &l_cpId ), "nc_open(checkpoint)" );

	int   l_nx = 0, l_ny = 0;
	float l_dx = 0, l_dy = 0, l_originX = 0, l_originY = 0, l_endTime = 0, l_simTime = 0;

	checkNc( nc_get_att_int(   l_cpId, NC_GLOBAL, "nx",       &l_nx      ), "nc_get_att(nx)"       );
	checkNc( nc_get_att_int(   l_cpId, NC_GLOBAL, "ny",       &l_ny      ), "nc_get_att(ny)"       );
	checkNc( nc_get_att_float( l_cpId, NC_GLOBAL, "dx",       &l_dx      ), "nc_get_att(dx)"       );
	checkNc( nc_get_att_float( l_cpId, NC_GLOBAL, "dy",       &l_dy      ), "nc_get_att(dy)"       );
	checkNc( nc_get_att_float( l_cpId, NC_GLOBAL, "origin_x", &l_originX ), "nc_get_att(origin_x)" );
	checkNc( nc_get_att_float( l_cpId, NC_GLOBAL, "origin_y", &l_originY ), "nc_get_att(origin_y)" );
	checkNc( nc_get_att_float( l_cpId, NC_GLOBAL, "end_time", &l_endTime ), "nc_get_att(end_time)" );
	checkNc( nc_get_att_float( l_cpId, NC_GLOBAL, "sim_time", &l_simTime ), "nc_get_att(sim_time)" );

	nc_close( l_cpId );

	l_cp.nx      = static_cast<t_idx>( l_nx );
	l_cp.ny      = static_cast<t_idx>( l_ny );
	l_cp.dx      = static_cast<t_real>( l_dx );
	l_cp.dy      = static_cast<t_real>( l_dy );
	l_cp.originX = static_cast<t_real>( l_originX );
	l_cp.originY = static_cast<t_real>( l_originY );
	l_cp.endTime = static_cast<t_real>( l_endTime );
	l_cp.simTime = static_cast<t_real>( l_simTime );

	// ── 2. Derive solution.nc path (sibling file in the same directory) ────────
	std::filesystem::path l_solutionFile =
		std::filesystem::path( i_checkpointFile ).parent_path() / "solution.nc";

	if( !std::filesystem::exists( l_solutionFile ) ) {
		throw std::runtime_error( "solution file not found: " + l_solutionFile.string() );
	}

	// ── 3. Open solution.nc and find the last valid time step ──────────────────
	int l_ncId = -1;
	checkNc( nc_open( l_solutionFile.string().c_str(), NC_NOWRITE, &l_ncId ), "nc_open(solution)" );

	int         l_dimTimeId = -1;
	std::size_t l_nSteps    = 0;
	checkNc( nc_inq_dimid(  l_ncId, "time", &l_dimTimeId ), "nc_inq_dimid(time)" );
	checkNc( nc_inq_dimlen( l_ncId, l_dimTimeId, &l_nSteps ), "nc_inq_dimlen(time)" );

	if( l_nSteps == 0 ) {
		nc_close( l_ncId );
		throw std::runtime_error( "solution.nc contains no time steps" );
	}

	int l_varTimeId = -1, l_varHId = -1, l_varHuId = -1, l_varHvId = -1, l_varBId = -1;
	checkNc( nc_inq_varid( l_ncId, "time",       &l_varTimeId ), "nc_inq_varid(time)"       );
	checkNc( nc_inq_varid( l_ncId, "height",     &l_varHId    ), "nc_inq_varid(height)"     );
	checkNc( nc_inq_varid( l_ncId, "momentum_x", &l_varHuId   ), "nc_inq_varid(momentum_x)" );
	checkNc( nc_inq_varid( l_ncId, "bathymetry", &l_varBId    ), "nc_inq_varid(bathymetry)" );

	int  l_varHvId2     = -1;
	bool l_hasHv        = ( nc_inq_varid( l_ncId, "momentum_y", &l_varHvId2 ) == NC_NOERR );
	l_varHvId           = l_varHvId2;

	// scan backwards for last step whose time == end_time (commit marker)
	std::size_t l_lastStep = std::numeric_limits<std::size_t>::max();
	for( std::size_t l_i = l_nSteps; l_i-- > 0; ) {
		float       l_t   = 0;
		std::size_t l_s[] = { l_i };
		std::size_t l_c[] = { 1   };
		nc_get_vara_float( l_ncId, l_varTimeId, l_s, l_c, &l_t );
		if( l_t == l_endTime ) { l_lastStep = l_i; break; }
	}

	// fallback: last step with a valid finite time value
	if( l_lastStep == std::numeric_limits<std::size_t>::max() ) {
		for( std::size_t l_i = l_nSteps; l_i-- > 0; ) {
			float       l_t   = 0;
			std::size_t l_s[] = { l_i };
			std::size_t l_c[] = { 1   };
			nc_get_vara_float( l_ncId, l_varTimeId, l_s, l_c, &l_t );
			if( std::isfinite( l_t ) && l_t < 9e+36f ) { l_lastStep = l_i; break; }
		}
	}

	if( l_lastStep == std::numeric_limits<std::size_t>::max() ) {
		nc_close( l_ncId );
		throw std::runtime_error( "solution.nc contains no valid time step" );
	}

	// ── 4. Read h, hu, hv, b at the last valid time step ──────────────────────
	t_idx l_size = l_cp.nx * l_cp.ny;
	l_cp.h.resize(  l_size );
	l_cp.hu.resize( l_size );
	l_cp.hv.resize( l_size, 0 );
	l_cp.b.resize(  l_size );

	std::size_t l_startTYX[3] = { l_lastStep, 0, 0 };
	std::size_t l_countTYX[3] = { 1, l_cp.ny, l_cp.nx };

	checkNc( nc_get_vara_float( l_ncId, l_varHId,  l_startTYX, l_countTYX, l_cp.h.data()  ), "nc_get_vara(height)"     );
	checkNc( nc_get_vara_float( l_ncId, l_varHuId, l_startTYX, l_countTYX, l_cp.hu.data() ), "nc_get_vara(momentum_x)" );

	if( l_hasHv ) {
		checkNc( nc_get_vara_float( l_ncId, l_varHvId, l_startTYX, l_countTYX, l_cp.hv.data() ), "nc_get_vara(momentum_y)" );
	}

	std::size_t l_startYX[2] = { 0, 0 };
	std::size_t l_countYX[2] = { l_cp.ny, l_cp.nx };
	checkNc( nc_get_vara_float( l_ncId, l_varBId, l_startYX, l_countYX, l_cp.b.data() ), "nc_get_vara(bathymetry)" );

	nc_close( l_ncId );
	return l_cp;
}

tsunami_lab::io::NetCdf::Data tsunami_lab::io::NetCdf::read( std::string const & i_bathymetryFile,
														   std::string const & i_displacementFile ) {
	Data l_data;

	// Open bathymetry file
	int l_ncId = -1;
	checkNc( nc_open( i_bathymetryFile.c_str(), NC_NOWRITE, &l_ncId ), "nc_open (bathymetry)" );

	// Query dimensions
	int l_dimXId = -1;
	int l_dimYId = -1;
	t_idx l_nx = 0;
	t_idx l_ny = 0;

	// Get dimension ids for x and y dimensions and their lengths
	checkNc( nc_inq_dimid( l_ncId, "x", &l_dimXId ), "nc_inq_dimid (x)" );
	checkNc( nc_inq_dimid( l_ncId, "y", &l_dimYId ), "nc_inq_dimid (y)" );
	checkNc( nc_inq_dimlen( l_ncId, l_dimXId, &l_nx ), "nc_inq_dimlen (x)" );
	checkNc( nc_inq_dimlen( l_ncId, l_dimYId, &l_ny ), "nc_inq_dimlen (y)" );

	l_data.gridNx = l_nx;
	l_data.gridNy = l_ny;

	// Get variable ids for x and y coordinates
	int l_varXId = -1;
	int l_varYId = -1;
	checkNc( nc_inq_varid( l_ncId, "x", &l_varXId ), "nc_inq_varid (x)" );
	checkNc( nc_inq_varid( l_ncId, "y", &l_varYId ), "nc_inq_varid (y)" );

	// Set vectors to have the correct size
	l_data.xCoords.resize( l_data.gridNx );
	l_data.yCoords.resize( l_data.gridNy );

	// Read coordinate data for x and y dimensions
	checkNc( nc_get_var_float( l_ncId, l_varXId, l_data.xCoords.data() ), "nc_get_var_float (x)" );
	checkNc( nc_get_var_float( l_ncId, l_varYId, l_data.yCoords.data() ), "nc_get_var_float (y)" );

	// Get variable id for bathymetry data (named z according to COARDS convention)
	int l_varZId = -1;
	checkNc( nc_inq_varid( l_ncId, "z", &l_varZId ), "nc_inq_varid (z)" );

	// Determine dimension ordering of z variable
	int l_zNdims = 0;
	int l_zDimIds[NC_MAX_VAR_DIMS];
	checkNc( nc_inq_var( l_ncId, l_varZId, nullptr, nullptr, &l_zNdims, l_zDimIds, nullptr ), "nc_inq_var (z)" );

	// Set correct size for bathymetry data vector
	l_data.bathymetryData.resize( l_data.gridNx * l_data.gridNy );

	if( l_zNdims == 2 && l_zDimIds[0] == l_dimXId ) {
		// (x, y) ordering - read raw and transpose to (y, x)
		std::vector<t_real> l_raw( l_data.gridNx * l_data.gridNy );
		checkNc( nc_get_var_float( l_ncId, l_varZId, l_raw.data() ), "nc_get_var_float (z)" );
		for( t_idx l_iy = 0; l_iy < l_data.gridNy; l_iy++ ) {
			for( t_idx l_ix = 0; l_ix < l_data.gridNx; l_ix++ ) {
				l_data.bathymetryData[l_iy * l_data.gridNx + l_ix] = l_raw[l_ix * l_data.gridNy + l_iy];
			}
		}
	} else {
		// (y, x) ordering - read directly
		checkNc( nc_get_var_float( l_ncId, l_varZId, l_data.bathymetryData.data() ), "nc_get_var_float (z)" );
	}

	// Close bathymetry file
	nc_close( l_ncId );

	// Read displacement data if file is provided
	if( !i_displacementFile.empty() ) {
		// Open displacement file
		int l_ncDispId = -1;
		checkNc( nc_open( i_displacementFile.c_str(), NC_NOWRITE, &l_ncDispId ), "nc_open (displacement)" );

		// Query dimension IDs and sizes from displacement file
		int l_dispDimXId = -1;
		int l_dispDimYId = -1;
		t_idx l_dispNx = 0;
		t_idx l_dispNy = 0;
		checkNc( nc_inq_dimid( l_ncDispId, "x", &l_dispDimXId ), "nc_inq_dimid (displacement x)" );
		checkNc( nc_inq_dimid( l_ncDispId, "y", &l_dispDimYId ), "nc_inq_dimid (displacement y)" );
		checkNc( nc_inq_dimlen( l_ncDispId, l_dispDimXId, &l_dispNx ), "nc_inq_dimlen (displacement x)" );
		checkNc( nc_inq_dimlen( l_ncDispId, l_dispDimYId, &l_dispNy ), "nc_inq_dimlen (displacement y)" );

		// Read displacement file's own x and y coordinates
		int l_varDispXId = -1;
		int l_varDispYId = -1;
		checkNc( nc_inq_varid( l_ncDispId, "x", &l_varDispXId ), "nc_inq_varid (displacement x)" );
		checkNc( nc_inq_varid( l_ncDispId, "y", &l_varDispYId ), "nc_inq_varid (displacement y)" );
		std::vector<t_real> l_dispXCoords( l_dispNx );
		std::vector<t_real> l_dispYCoords( l_dispNy );
		checkNc( nc_get_var_float( l_ncDispId, l_varDispXId, l_dispXCoords.data() ), "nc_get_var_float (displacement x)" );
		checkNc( nc_get_var_float( l_ncDispId, l_varDispYId, l_dispYCoords.data() ), "nc_get_var_float (displacement y)" );

		// Get variable id for displacement data (named z according to COARDS convention)
		int l_varDispId = -1;
		checkNc( nc_inq_varid( l_ncDispId, "z", &l_varDispId ), "nc_inq_varid (displacement z)" );

		// Determine dimension ordering of displacement z variable
		int l_dispNdims = 0;
		int l_dispDimIds[NC_MAX_VAR_DIMS];
		checkNc( nc_inq_var( l_ncDispId, l_varDispId, nullptr, nullptr, &l_dispNdims, l_dispDimIds, nullptr ), "nc_inq_var (displacement z)" );
		bool l_dispIsXY = ( l_dispNdims == 2 && l_dispDimIds[0] == l_dispDimXId );

		// Read raw displacement data with size
		std::vector<t_real> l_dispRaw( l_dispNx * l_dispNy );
		checkNc( nc_get_var_float( l_ncDispId, l_varDispId, l_dispRaw.data() ), "nc_get_var_float (displacement z)" );

		// Store displacement with its own grid dimensions
		l_data.dispNx = l_dispNx;
		l_data.dispNy = l_dispNy;
		l_data.dispXCoords = std::move( l_dispXCoords );
		l_data.dispYCoords = std::move( l_dispYCoords );

		// Store raw displacement data in its native layout, normalised to (y, x) row-major
		l_data.displacementData.resize( l_dispNx * l_dispNy );
		if( l_dispIsXY ) {
			// Transpose from (x, y) to (y, x)
			for( t_idx l_iy = 0; l_iy < l_dispNy; l_iy++ ) {
				for( t_idx l_ix = 0; l_ix < l_dispNx; l_ix++ ) {
					l_data.displacementData[l_iy * l_dispNx + l_ix] = l_dispRaw[l_ix * l_dispNy + l_iy];
				}
			}
		} else {
			// Already (y, x) — copy directly
			l_data.displacementData = std::move( l_dispRaw );
		}

		// Close displacement file
		nc_close( l_ncDispId );
	}

	return l_data;
}

tsunami_lab::t_real tsunami_lab::io::NetCdf::Data::getBathymetry( t_real i_x, t_real i_y ) const {
	if( bathymetryData.empty() || xCoords.empty() || yCoords.empty() ) {
		return 0;
	}

	// Find closest grid indices using nearest-neighbor
	t_idx l_ix = 0;
	t_real l_minDistX = std::abs( i_x - xCoords[0] );
	for( t_idx l_i = 1; l_i < gridNx; l_i++ ) {
		t_real l_dist = std::abs( i_x - xCoords[l_i] );
		if( l_dist < l_minDistX ) {
			l_minDistX = l_dist;
			l_ix = l_i;
		}
	}

	t_idx l_iy = 0;
	t_real l_minDistY = std::abs( i_y - yCoords[0] );
	for( t_idx l_i = 1; l_i < gridNy; l_i++ ) {
		t_real l_dist = std::abs( i_y - yCoords[l_i] );
		if( l_dist < l_minDistY ) {
			l_minDistY = l_dist;
			l_iy = l_i;
		}
	}

	// Return value from grid
	return bathymetryData[l_iy * gridNx + l_ix];
}

tsunami_lab::t_real tsunami_lab::io::NetCdf::Data::getDisplacement( t_real i_x, t_real i_y ) const {
	if( displacementData.empty() || dispXCoords.empty() || dispYCoords.empty() ) {
		return 0;
	}

	// Return 0 for any query point outside the displacement file's coordinate range
	t_real l_xMin = std::min( dispXCoords.front(), dispXCoords.back() );
	t_real l_xMax = std::max( dispXCoords.front(), dispXCoords.back() );
	t_real l_yMin = std::min( dispYCoords.front(), dispYCoords.back() );
	t_real l_yMax = std::max( dispYCoords.front(), dispYCoords.back() );
	if( i_x < l_xMin || i_x > l_xMax || i_y < l_yMin || i_y > l_yMax ) {
		return 0;
	}

	// Find closest grid indices using nearest-neighbor in the displacement grid
	t_idx l_ix = 0;
	t_real l_minDistX = std::abs( i_x - dispXCoords[0] );
	for( t_idx l_i = 1; l_i < dispNx; l_i++ ) {
		t_real l_dist = std::abs( i_x - dispXCoords[l_i] );
		if( l_dist < l_minDistX ) {
			l_minDistX = l_dist;
			l_ix = l_i;
		}
	}

	t_idx l_iy = 0;
	t_real l_minDistY = std::abs( i_y - dispYCoords[0] );
	for( t_idx l_i = 1; l_i < dispNy; l_i++ ) {
		t_real l_dist = std::abs( i_y - dispYCoords[l_i] );
		if( l_dist < l_minDistY ) {
			l_minDistY = l_dist;
			l_iy = l_i;
		}
	}

	// Return value from displacement grid (row-major y, x)
	return displacementData[l_iy * dispNx + l_ix];
}

tsunami_lab::io::NetCdf::~NetCdf() {
	if( m_ncId != c_invalidId ) {
		nc_close( toNcId( m_ncId ) );
		m_ncId = c_invalidId;
	}
	if( m_checkpointNcId != c_invalidId ) {
		nc_close( toNcId( m_checkpointNcId ) );
		m_checkpointNcId = c_invalidId;
	}
}

