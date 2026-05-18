/**
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Tests the 2D tsunami event setup with netCDF input.
 **/
#include <catch2/catch.hpp>
#include "TsunamiEvent2d.h"
#include "../../io/NetCdf.h"
#include <cstdio>
#include <netcdf.h>
#include <vector>
#include <cmath>

namespace {
	using tsunami_lab::t_idx;
	using tsunami_lab::t_real;

	// Helper function to create a simple test netCDF file with bathymetry data
	void createTestBathymetryFile( std::string const & i_filePath,
	                                 t_idx i_nx,
	                                 t_idx i_ny,
	                                 t_real i_xMin,
	                                 t_real i_xMax,
	                                 t_real i_yMin,
	                                 t_real i_yMax ) {
		int l_ncId = -1;
		int l_retVal = nc_create( i_filePath.c_str(), NC_NETCDF4 | NC_CLOBBER, &l_ncId );
		if( l_retVal != NC_NOERR ) {
			throw std::runtime_error( "Failed to create netCDF file" );
		}

		// Define dimensions
		int l_dimXId = -1;
		int l_dimYId = -1;
		nc_def_dim( l_ncId, "x", i_nx, &l_dimXId );
		nc_def_dim( l_ncId, "y", i_ny, &l_dimYId );

		// Define coordinate variables
		int l_varXId = -1;
		int l_varYId = -1;
		nc_def_var( l_ncId, "x", NC_FLOAT, 1, &l_dimXId, &l_varXId );
		nc_def_var( l_ncId, "y", NC_FLOAT, 1, &l_dimYId, &l_varYId );

		// Define bathymetry variable
		int l_dimsYX[2] = { l_dimYId, l_dimXId };
		int l_varZId = -1;
		nc_def_var( l_ncId, "z", NC_FLOAT, 2, l_dimsYX, &l_varZId );

		// End define mode
		nc_enddef( l_ncId );

		// Write x coordinates
		std::vector<float> l_xCoords( i_nx );
		for( t_idx l_i = 0; l_i < i_nx; l_i++ ) {
			l_xCoords[l_i] = i_xMin + ( i_xMax - i_xMin ) * static_cast<float>( l_i ) / static_cast<float>( i_nx - 1 );
		}
		nc_put_var_float( l_ncId, l_varXId, l_xCoords.data() );

		// Write y coordinates
		std::vector<float> l_yCoords( i_ny );
		for( t_idx l_i = 0; l_i < i_ny; l_i++ ) {
			l_yCoords[l_i] = i_yMin + ( i_yMax - i_yMin ) * static_cast<float>( l_i ) / static_cast<float>( i_ny - 1 );
		}
		nc_put_var_float( l_ncId, l_varYId, l_yCoords.data() );

		// Write bathymetry data (constant depth)
		std::vector<float> l_bathyData( i_nx * i_ny, -100.0f );
		nc_put_var_float( l_ncId, l_varZId, l_bathyData.data() );

		nc_close( l_ncId );
	}

	// Helper function to create a test displacement file
	void createTestDisplacementFile( std::string const & i_filePath,
	                                   t_idx i_nx,
	                                   t_idx i_ny,
	                                   t_real i_xMin,
	                                   t_real i_xMax,
	                                   t_real i_yMin,
	                                   t_real i_yMax ) {
		int l_ncId = -1;
		int l_retVal = nc_create( i_filePath.c_str(), NC_NETCDF4 | NC_CLOBBER, &l_ncId );
		if( l_retVal != NC_NOERR ) {
			throw std::runtime_error( "Failed to create netCDF file" );
		}

		// Define dimensions
		int l_dimXId = -1;
		int l_dimYId = -1;
		nc_def_dim( l_ncId, "x", i_nx, &l_dimXId );
		nc_def_dim( l_ncId, "y", i_ny, &l_dimYId );

		// Define coordinate variables
		int l_varXId = -1;
		int l_varYId = -1;
		nc_def_var( l_ncId, "x", NC_FLOAT, 1, &l_dimXId, &l_varXId );
		nc_def_var( l_ncId, "y", NC_FLOAT, 1, &l_dimYId, &l_varYId );

		// Define displacement variable
		int l_dimsYX[2] = { l_dimYId, l_dimXId };
		int l_varZId = -1;
		nc_def_var( l_ncId, "z", NC_FLOAT, 2, l_dimsYX, &l_varZId );

		// End define mode
		nc_enddef( l_ncId );

		// Write x coordinates
		std::vector<float> l_xCoords( i_nx );
		for( t_idx l_i = 0; l_i < i_nx; l_i++ ) {
			l_xCoords[l_i] = i_xMin + ( i_xMax - i_xMin ) * static_cast<float>( l_i ) / static_cast<float>( i_nx - 1 );
		}
		nc_put_var_float( l_ncId, l_varXId, l_xCoords.data() );

		// Write y coordinates
		std::vector<float> l_yCoords( i_ny );
		for( t_idx l_i = 0; l_i < i_ny; l_i++ ) {
			l_yCoords[l_i] = i_yMin + ( i_yMax - i_yMin ) * static_cast<float>( l_i ) / static_cast<float>( i_ny - 1 );
		}
		nc_put_var_float( l_ncId, l_varYId, l_yCoords.data() );

		// Write displacement data (sinusoidal pattern in the center)
		std::vector<float> l_dispData( i_nx * i_ny, 0.0f );
		
		// Define center and size of displacement
		t_real l_centerX = ( i_xMin + i_xMax ) / 2.0;
		t_real l_centerY = ( i_yMin + i_yMax ) / 2.0;
		t_real l_sizeX = 500.0;  // +/- 500m
		t_real l_sizeY = 500.0;

		for( t_idx l_iy = 0; l_iy < i_ny; l_iy++ ) {
			for( t_idx l_ix = 0; l_ix < i_nx; l_ix++ ) {
				t_real l_x = l_xCoords[l_ix];
				t_real l_y = l_yCoords[l_iy];
				t_real l_dx = l_x - l_centerX;
				t_real l_dy = l_y - l_centerY;

				// Only apply displacement within the specified region
				if( std::abs( l_dx ) <= l_sizeX && std::abs( l_dy ) <= l_sizeY ) {
					t_real l_fx = std::sin( ( l_dx / l_sizeX + 1.0 ) * M_PI );
					t_real l_gy = -( l_dy / l_sizeY ) * ( l_dy / l_sizeY ) + 1.0;
					l_dispData[l_iy * i_nx + l_ix] = 5.0f * l_fx * l_gy;
				}
			}
		}
		
		nc_put_var_float( l_ncId, l_varZId, l_dispData.data() );

		nc_close( l_ncId );
	}
}

TEST_CASE( "Test TsunamiEvent2d setup reads netCDF files correctly.", "[TsunamiEvent2d]" ) {
	std::string l_bathyFile = "test_bathymetry_2d.nc";
	std::string l_dispFile = "test_displacement_2d.nc";
	
	// Clean up any existing test files
	std::remove( l_bathyFile.c_str() );
	std::remove( l_dispFile.c_str() );

	// Create test files
	createTestBathymetryFile( l_bathyFile, 10, 10, 0.0, 10000.0, 0.0, 10000.0 );
	createTestDisplacementFile( l_dispFile, 10, 10, 0.0, 10000.0, 0.0, 10000.0 );

	SECTION( "Read bathymetry and displacement from netCDF files" ) {
		tsunami_lab::setups::TsunamiEvent2d l_setup( l_bathyFile, l_dispFile );

		// Test at a location where bathymetry is negative (underwater)
		// Height should be max(-b_in, delta) = max(100, 20) = 100
		REQUIRE( l_setup.getHeight( 5000.0, 5000.0 ) == Approx( 100.0 ) );

		// Momentum should be zero
		REQUIRE( l_setup.getMomentumX( 5000.0, 5000.0 ) == Approx( 0.0 ) );
		REQUIRE( l_setup.getMomentumY( 5000.0, 5000.0 ) == Approx( 0.0 ) );

		// Bathymetry should be min(b_in, -delta) + d = min(-100, -20) + d = -100 + d
		// At center (5000, 5000), displacement d >= 0, so bathymetry should be >= -100
		REQUIRE( l_setup.getBathymetry( 5000.0, 5000.0 ) >= -100.0f );
	}

	SECTION( "Test nearest-neighbor interpolation" ) {
		tsunami_lab::setups::TsunamiEvent2d l_setup( l_bathyFile, l_dispFile );

		// Test at grid points
		REQUIRE( l_setup.getHeight( 0.0, 0.0 ) == Approx( 100.0 ) );
		REQUIRE( l_setup.getHeight( 10000.0, 10000.0 ) == Approx( 100.0 ) );

		// Test between grid points (should use nearest neighbor)
		t_real l_h1 = l_setup.getHeight( 1234.5, 5678.9 );
		REQUIRE( l_h1 == Approx( 100.0 ) );  // Should still be 100 (constant bathymetry)
	}

	// Clean up test files
	std::remove( l_bathyFile.c_str() );
	std::remove( l_dispFile.c_str() );
}
