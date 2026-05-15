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
											4,
											0.0,
											l_h + 5,
											l_b + 5,
											l_hu + 5,
											l_hv + 5,
											l_filePath );

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
	std::remove( l_filePath.c_str() );
}
