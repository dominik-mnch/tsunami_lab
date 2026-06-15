#include <catch2/catch.hpp>
#include "NetCdf.h"
#include <cstdio>
#include <netcdf.h>
#include <string>
#include <vector>

TEST_CASE( "NetCDF writer appends timesteps and ignores ghost cells.", "[NetCdf]" ) {
	std::string l_filePath = "test_solution.nc";
	std::remove( l_filePath.c_str() );

	tsunami_lab::t_real l_h[16] = {
			-100, -101, -102, -103,
			-104,    1,    2, -105,
			-106,    3,    4, -107,
			-108, -109, -110, -111
	};
	tsunami_lab::t_real l_b[16] = {
			-200, -201, -202, -203,
			-204,   10,   20, -205,
			-206,   30,   40, -207,
			-208, -209, -210, -211
	};
	tsunami_lab::t_real l_hu[16] = {
			-300, -301, -302, -303,
			-304,   11,   12, -305,
			-306,   13,   14, -307,
			-308, -309, -310, -311
	};
	tsunami_lab::t_real l_hv[16] = {
			-400, -401, -402, -403,
			-404,   21,   22, -405,
			-406,   23,   24, -407,
			-408, -409, -410, -411
	};

	// Create object once before writing multiple timesteps
	tsunami_lab::io::NetCdf l_writer( 1.0,
											1.0,
											0.0,
											0.0,
											2,
											2,
											1,
											4,
											0.0,
											1.0,
											1,
											"2d",
											"artificial_tsunami_2d",
											"nx=2;ny=2;k=1;solver=1;propagation=2d;setup=artificial_tsunami_2d",
											l_h + 5,
											l_b + 5,
											l_hu + 5,
											l_hv + 5,
											l_filePath,
											true );

	// Write first timestep
	l_writer.writeTimeStep( 0.0 );

	// Modify data and write second timestep
	l_h[5] = 101;
	l_h[6] = 102;
	l_h[9] = 103;
	l_h[10] = 104;
	l_writer.writeTimeStep( 1.0 );

	// Verify the output file
	int l_ncId = -1;
	REQUIRE( nc_open( l_filePath.c_str(), NC_NOWRITE, &l_ncId ) == NC_NOERR );

	int l_dimXId = -1;
	int l_dimYId = -1;
	int l_dimTimeId = -1;
	std::size_t l_dimX = 0;
	std::size_t l_dimY = 0;
	std::size_t l_dimTime = 0;

	REQUIRE( nc_inq_dimid( l_ncId, "x", &l_dimXId ) == NC_NOERR );
	REQUIRE( nc_inq_dimid( l_ncId, "y", &l_dimYId ) == NC_NOERR );
	REQUIRE( nc_inq_dimid( l_ncId, "time", &l_dimTimeId ) == NC_NOERR );
	REQUIRE( nc_inq_dimlen( l_ncId, l_dimXId, &l_dimX ) == NC_NOERR );
	REQUIRE( nc_inq_dimlen( l_ncId, l_dimYId, &l_dimY ) == NC_NOERR );
	REQUIRE( nc_inq_dimlen( l_ncId, l_dimTimeId, &l_dimTime ) == NC_NOERR );

	REQUIRE( l_dimX == 2 );
	REQUIRE( l_dimY == 2 );
	REQUIRE( l_dimTime == 2 );

	int l_varTimeId = -1;
	int l_varHeightId = -1;
	int l_varBathyId = -1;
	REQUIRE( nc_inq_varid( l_ncId, "time", &l_varTimeId ) == NC_NOERR );
	REQUIRE( nc_inq_varid( l_ncId, "height", &l_varHeightId ) == NC_NOERR );
	REQUIRE( nc_inq_varid( l_ncId, "bathymetry", &l_varBathyId ) == NC_NOERR );

	std::vector<tsunami_lab::t_real> l_times( 2 );
	REQUIRE( nc_get_var_float( l_ncId, l_varTimeId, l_times.data() ) == NC_NOERR );
	REQUIRE( l_times[0] == Approx( 0.0 ) );
	REQUIRE( l_times[1] == Approx( 1.0 ) );

	std::vector<tsunami_lab::t_real> l_height( 2 * 2 * 2 );
	REQUIRE( nc_get_var_float( l_ncId, l_varHeightId, l_height.data() ) == NC_NOERR );

	REQUIRE( l_height[0] == Approx( 1.0 ) );
	REQUIRE( l_height[1] == Approx( 2.0 ) );
	REQUIRE( l_height[2] == Approx( 3.0 ) );
	REQUIRE( l_height[3] == Approx( 4.0 ) );

	REQUIRE( l_height[4] == Approx( 101.0 ) );
	REQUIRE( l_height[5] == Approx( 102.0 ) );
	REQUIRE( l_height[6] == Approx( 103.0 ) );
	REQUIRE( l_height[7] == Approx( 104.0 ) );

	std::vector<tsunami_lab::t_real> l_bathy( 2 * 2 );
	REQUIRE( nc_get_var_float( l_ncId, l_varBathyId, l_bathy.data() ) == NC_NOERR );
	REQUIRE( l_bathy[0] == Approx( 10.0 ) );
	REQUIRE( l_bathy[1] == Approx( 20.0 ) );
	REQUIRE( l_bathy[2] == Approx( 30.0 ) );
	REQUIRE( l_bathy[3] == Approx( 40.0 ) );

	REQUIRE( nc_close( l_ncId ) == NC_NOERR );

	int l_checkpointNcId = -1;
	REQUIRE( nc_open( "checkpoint.nc", NC_NOWRITE, &l_checkpointNcId ) == NC_NOERR );

	int l_k = 0;
	int l_solverMode = -1;
	REQUIRE( nc_get_att_int( l_checkpointNcId, NC_GLOBAL, "k", &l_k ) == NC_NOERR );
	REQUIRE( nc_get_att_int( l_checkpointNcId, NC_GLOBAL, "solver_mode", &l_solverMode ) == NC_NOERR );
	REQUIRE( l_k == 1 );
	REQUIRE( l_solverMode == 1 );

	size_t l_propagationLen = 0;
	size_t l_setupLen = 0;
	size_t l_inputSignatureLen = 0;
	REQUIRE( nc_inq_attlen( l_checkpointNcId, NC_GLOBAL, "propagation", &l_propagationLen ) == NC_NOERR );
	REQUIRE( nc_inq_attlen( l_checkpointNcId, NC_GLOBAL, "setup", &l_setupLen ) == NC_NOERR );
	REQUIRE( nc_inq_attlen( l_checkpointNcId, NC_GLOBAL, "input_signature", &l_inputSignatureLen ) == NC_NOERR );
	std::string l_propagation( l_propagationLen, '\0' );
	std::string l_setup( l_setupLen, '\0' );
	std::string l_inputSignature( l_inputSignatureLen, '\0' );
	REQUIRE( nc_get_att_text( l_checkpointNcId, NC_GLOBAL, "propagation", l_propagation.data() ) == NC_NOERR );
	REQUIRE( nc_get_att_text( l_checkpointNcId, NC_GLOBAL, "setup", l_setup.data() ) == NC_NOERR );
	REQUIRE( nc_get_att_text( l_checkpointNcId, NC_GLOBAL, "input_signature", l_inputSignature.data() ) == NC_NOERR );
	REQUIRE( l_propagation == "2d" );
	REQUIRE( l_setup == "artificial_tsunami_2d" );
	REQUIRE( l_inputSignature == "nx=2;ny=2;k=1;solver=1;propagation=2d;setup=artificial_tsunami_2d" );
	REQUIRE( nc_close( l_checkpointNcId ) == NC_NOERR );
	std::remove( l_filePath.c_str() );
	std::remove( "checkpoint.nc" );
}

TEST_CASE( "NetCDF writer keeps partial downsample blocks at domain edges.", "[NetCdf]" ) {
	std::string l_filePath = "test_solution_coarse_edges.nc";
	std::remove( l_filePath.c_str() );
	std::remove( "checkpoint.nc" );

	tsunami_lab::t_real l_h[6]  = { 1, 2, 3,
	                               4, 5, 6 };
	tsunami_lab::t_real l_b[6]  = { 10, 20, 30,
	                               40, 50, 60 };
	tsunami_lab::t_real l_hu[6] = { 0, 0, 0,
	                               0, 0, 0 };
	tsunami_lab::t_real l_hv[6] = { 0, 0, 0,
	                               0, 0, 0 };

	{
		tsunami_lab::io::NetCdf l_writer( 1.0,
		                                      1.0,
		                                      0.0,
		                                      0.0,
		                                      3,
		                                      2,
		                                      2,
		                                      3,
		                                      0.0,
		                                      1.0,
		                                      1,
		                                      "2d",
		                                      "artificial_tsunami_2d",
		                                      "nx=3;ny=2;k=2",
		                                      l_h,
		                                      l_b,
		                                      l_hu,
		                                      l_hv,
			                                      l_filePath,
			                                      false );
		l_writer.writeTimeStep( 0.0 );
	}

	int l_ncId = -1;
	REQUIRE( nc_open( l_filePath.c_str(), NC_NOWRITE, &l_ncId ) == NC_NOERR );

	int l_dimXId = -1;
	int l_dimYId = -1;
	std::size_t l_dimX = 0;
	std::size_t l_dimY = 0;
	REQUIRE( nc_inq_dimid( l_ncId, "x", &l_dimXId ) == NC_NOERR );
	REQUIRE( nc_inq_dimid( l_ncId, "y", &l_dimYId ) == NC_NOERR );
	REQUIRE( nc_inq_dimlen( l_ncId, l_dimXId, &l_dimX ) == NC_NOERR );
	REQUIRE( nc_inq_dimlen( l_ncId, l_dimYId, &l_dimY ) == NC_NOERR );
	REQUIRE( l_dimX == 2 );
	REQUIRE( l_dimY == 1 );

	int l_varXId = -1;
	int l_varYId = -1;
	int l_varHeightId = -1;
	int l_varBathyId = -1;
	REQUIRE( nc_inq_varid( l_ncId, "x", &l_varXId ) == NC_NOERR );
	REQUIRE( nc_inq_varid( l_ncId, "y", &l_varYId ) == NC_NOERR );
	REQUIRE( nc_inq_varid( l_ncId, "height", &l_varHeightId ) == NC_NOERR );
	REQUIRE( nc_inq_varid( l_ncId, "bathymetry", &l_varBathyId ) == NC_NOERR );

	std::vector<tsunami_lab::t_real> l_x( 2 );
	std::vector<tsunami_lab::t_real> l_y( 1 );
	REQUIRE( nc_get_var_float( l_ncId, l_varXId, l_x.data() ) == NC_NOERR );
	REQUIRE( nc_get_var_float( l_ncId, l_varYId, l_y.data() ) == NC_NOERR );
	REQUIRE( l_x[0] == Approx( 1.0 ) );
	REQUIRE( l_x[1] == Approx( 2.5 ) );
	REQUIRE( l_y[0] == Approx( 1.0 ) );

	std::vector<tsunami_lab::t_real> l_height( 2 );
	std::vector<tsunami_lab::t_real> l_bathy( 2 );
	REQUIRE( nc_get_var_float( l_ncId, l_varHeightId, l_height.data() ) == NC_NOERR );
	REQUIRE( nc_get_var_float( l_ncId, l_varBathyId, l_bathy.data() ) == NC_NOERR );
	REQUIRE( l_height[0] == Approx( 3.0 ) );
	REQUIRE( l_height[1] == Approx( 4.5 ) );
	REQUIRE( l_bathy[0] == Approx( 30.0 ) );
	REQUIRE( l_bathy[1] == Approx( 45.0 ) );

	REQUIRE( nc_close( l_ncId ) == NC_NOERR );
	std::remove( l_filePath.c_str() );
	std::remove( "checkpoint.nc" );
}

TEST_CASE( "NetCDF writer can disable checkpoint metadata.", "[NetCdf]" ) {
	std::string l_filePath = "test_solution_no_checkpoint.nc";
	std::remove( l_filePath.c_str() );
	std::remove( "checkpoint.nc" );

	tsunami_lab::t_real l_h[1] = { 1 };
	tsunami_lab::t_real l_b[1] = { -10 };
	tsunami_lab::t_real l_hu[1] = { 0 };
	tsunami_lab::t_real l_hv[1] = { 0 };

	{
		tsunami_lab::io::NetCdf l_writer( 1.0,
		                                      1.0,
		                                      0.0,
		                                      0.0,
		                                      1,
		                                      1,
		                                      1,
		                                      1,
		                                      0.0,
		                                      1.0,
		                                      1,
		                                      "2d",
		                                      "artificial_tsunami_2d",
		                                      "nx=1;ny=1;k=1",
		                                      l_h,
		                                      l_b,
		                                      l_hu,
		                                      l_hv,
		                                      l_filePath,
		                                      false );
		l_writer.writeTimeStep( 0.0 );
	}

	int l_ncId = -1;
	REQUIRE( nc_open( l_filePath.c_str(), NC_NOWRITE, &l_ncId ) == NC_NOERR );
	REQUIRE( nc_close( l_ncId ) == NC_NOERR );
	REQUIRE( nc_open( "checkpoint.nc", NC_NOWRITE, &l_ncId ) == NC_ENOTFOUND );

	std::remove( l_filePath.c_str() );
	std::remove( "checkpoint.nc" );
}

TEST_CASE( "NetCDF reader reads bathymetry and displacement data.", "[NetCdf]" ) {
	std::string l_bathyFile = "test_reader_bathymetry.nc";
	std::string l_dispFile = "test_reader_displacement.nc";

	// Clean up any existing test files
	std::remove( l_bathyFile.c_str() );
	std::remove( l_dispFile.c_str() );

	// Create test bathymetry file
	int l_bathyNcId = -1;
	REQUIRE( nc_create( l_bathyFile.c_str(), NC_NETCDF4 | NC_CLOBBER, &l_bathyNcId ) == NC_NOERR );

	int l_dimXId = -1;
	int l_dimYId = -1;
	REQUIRE( nc_def_dim( l_bathyNcId, "x", 3, &l_dimXId ) == NC_NOERR );
	REQUIRE( nc_def_dim( l_bathyNcId, "y", 2, &l_dimYId ) == NC_NOERR );

	int l_varXId = -1;
	int l_varYId = -1;
	int l_varZId = -1;
	REQUIRE( nc_def_var( l_bathyNcId, "x", NC_FLOAT, 1, &l_dimXId, &l_varXId ) == NC_NOERR );
	REQUIRE( nc_def_var( l_bathyNcId, "y", NC_FLOAT, 1, &l_dimYId, &l_varYId ) == NC_NOERR );

	int l_dimsYX[2] = { l_dimYId, l_dimXId };
	REQUIRE( nc_def_var( l_bathyNcId, "z", NC_FLOAT, 2, l_dimsYX, &l_varZId ) == NC_NOERR );

	REQUIRE( nc_enddef( l_bathyNcId ) == NC_NOERR );

	// Write coordinate data
	std::vector<float> l_xCoords = { 0.0f, 100.0f, 200.0f };
	std::vector<float> l_yCoords = { 0.0f, 100.0f };
	REQUIRE( nc_put_var_float( l_bathyNcId, l_varXId, l_xCoords.data() ) == NC_NOERR );
	REQUIRE( nc_put_var_float( l_bathyNcId, l_varYId, l_yCoords.data() ) == NC_NOERR );

	// Write bathymetry data
	std::vector<float> l_bathyData = { -100.0f, -100.0f, -100.0f, -100.0f, -100.0f, -100.0f };
	REQUIRE( nc_put_var_float( l_bathyNcId, l_varZId, l_bathyData.data() ) == NC_NOERR );

	REQUIRE( nc_close( l_bathyNcId ) == NC_NOERR );

	// Create test displacement file
	int l_dispNcId = -1;
	REQUIRE( nc_create( l_dispFile.c_str(), NC_NETCDF4 | NC_CLOBBER, &l_dispNcId ) == NC_NOERR );

	REQUIRE( nc_def_dim( l_dispNcId, "x", 3, &l_dimXId ) == NC_NOERR );
	REQUIRE( nc_def_dim( l_dispNcId, "y", 2, &l_dimYId ) == NC_NOERR );

	REQUIRE( nc_def_var( l_dispNcId, "x", NC_FLOAT, 1, &l_dimXId, &l_varXId ) == NC_NOERR );
	REQUIRE( nc_def_var( l_dispNcId, "y", NC_FLOAT, 1, &l_dimYId, &l_varYId ) == NC_NOERR );

	REQUIRE( nc_def_var( l_dispNcId, "z", NC_FLOAT, 2, l_dimsYX, &l_varZId ) == NC_NOERR );

	REQUIRE( nc_enddef( l_dispNcId ) == NC_NOERR );

	REQUIRE( nc_put_var_float( l_dispNcId, l_varXId, l_xCoords.data() ) == NC_NOERR );
	REQUIRE( nc_put_var_float( l_dispNcId, l_varYId, l_yCoords.data() ) == NC_NOERR );

	// Write displacement data
	std::vector<float> l_dispData = { 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f };
	REQUIRE( nc_put_var_float( l_dispNcId, l_varZId, l_dispData.data() ) == NC_NOERR );

	REQUIRE( nc_close( l_dispNcId ) == NC_NOERR );

	SECTION( "Read bathymetry and displacement from netCDF files" ) {
		auto l_data = tsunami_lab::io::NetCdf::read( l_bathyFile, l_dispFile );

		REQUIRE( l_data.gridNx == 3 );
		REQUIRE( l_data.gridNy == 2 );

		REQUIRE( l_data.xCoords[0] == Approx( 0.0f ) );
		REQUIRE( l_data.xCoords[1] == Approx( 100.0f ) );
		REQUIRE( l_data.xCoords[2] == Approx( 200.0f ) );

		REQUIRE( l_data.yCoords[0] == Approx( 0.0f ) );
		REQUIRE( l_data.yCoords[1] == Approx( 100.0f ) );

		// Test bathymetry values at grid points
		REQUIRE( l_data.getBathymetry( 0.0f, 0.0f ) == Approx( -100.0f ) );
		REQUIRE( l_data.getBathymetry( 100.0f, 100.0f ) == Approx( -100.0f ) );
		REQUIRE( l_data.getBathymetry( 200.0f, 0.0f ) == Approx( -100.0f ) );
	}

	SECTION( "Test nearest-neighbor interpolation for bathymetry" ) {
		auto l_data = tsunami_lab::io::NetCdf::read( l_bathyFile, l_dispFile );

		// Test at off-grid points - should use nearest neighbor
		// Point (50, 50) is closest to (0, 0) or (100, 100)
		// Distance to (0, 0) = sqrt(50^2 + 50^2) = 70.7
		// Distance to (100, 100) = sqrt(50^2 + 50^2) = 70.7
		// Both equidistant, should get one of them
		float l_bathy = l_data.getBathymetry( 50.0f, 50.0f );
		REQUIRE( l_bathy == Approx( -100.0f ) );

		// Point (25, 25) is closest to (0, 0)
		l_bathy = l_data.getBathymetry( 25.0f, 25.0f );
		REQUIRE( l_bathy == Approx( -100.0f ) );
	}

	SECTION( "Test displacement values at grid points" ) {
		auto l_data = tsunami_lab::io::NetCdf::read( l_bathyFile, l_dispFile );

		// Displacement data layout: [0, 1, 2] for y=0, [3, 4, 5] for y=100
		REQUIRE( l_data.getDisplacement( 0.0f, 0.0f ) == Approx( 0.0f ) );
		REQUIRE( l_data.getDisplacement( 100.0f, 0.0f ) == Approx( 1.0f ) );
		REQUIRE( l_data.getDisplacement( 200.0f, 0.0f ) == Approx( 2.0f ) );
		REQUIRE( l_data.getDisplacement( 0.0f, 100.0f ) == Approx( 3.0f ) );
		REQUIRE( l_data.getDisplacement( 100.0f, 100.0f ) == Approx( 4.0f ) );
		REQUIRE( l_data.getDisplacement( 200.0f, 100.0f ) == Approx( 5.0f ) );
	}

	// Clean up test files
	std::remove( l_bathyFile.c_str() );
	std::remove( l_dispFile.c_str() );
}
