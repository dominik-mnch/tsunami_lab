/**
 * @author Madalena Schwarzkopf
 * @author Dominik Münch
 * 
 * @section DESCRIPTION
 * Unit tests for CUDA refression test infrastructure across grid sizes.  
 **/

#include <catch2/catch.hpp>
#include "../constants.h"
#include "CudaRegressionTest.h"
#include "../setups/CircularDamBreak2d/CircularDamBreak2d.h"

/**
 * Helper: allocate, fill, and run regression check for a given grid size.
 **/
static void runGridTest( tsunami_lab::t_idx i_nCellsX,
                         tsunami_lab::t_idx i_nCellsY ) {
    // --- CPU side --- 
    tsunami_lab::patches::WavePropagation2d l_waveProp( i_nCellsX,
                                                        i_nCellsY,
                                                        true);
    tsunami_lab::setups::CircularDamBreak2d l_setup( 10.0, 5.0,
                                                     0, 0, 0, 0,
                                                     0, 0, 10 );

    tsunami_lab::t_real l_domainSize = 50.0;
    tsunami_lab::t_real l_dx = l_domainSize / i_nCellsX;
    tsunami_lab::t_real l_dy = l_domainSize / i_nCellsY;

    for( tsunami_lab::t_idx l_iy = 0; l_iy < i_nCellsY; l_iy++ ) {
        for( tsunami_lab::t_idx l_ix = 0; l_ix < i_nCellsX; l_ix++ ) {
            tsunami_lab::t_real l_x = -l_domainSize / 2.0 + (l_ix + 0.5) * l_dx;
            tsunami_lab::t_real l_y = -l_domainSize / 2.0 + (l_iy + 0.5) * l_dy;

            l_waveProp.setHeight(    l_ix, l_iy, l_setup.getHeight(    l_x, l_y ) );
            l_waveProp.setMomentumX( l_ix, l_iy, l_setup.getMomentumX( l_x, l_y ) );
            l_waveProp.setMomentumY( l_ix, l_iy, l_setup.getMomentumY( l_x, l_y ) );
            l_waveProp.setBathymetry( l_ix, l_iy, l_setup.getBathymetry( l_x, l_y ) );
        }
    }

        // --- GPU side: allocate ---
    tsunami_lab::t_real* l_d_h  = nullptr;
    tsunami_lab::t_real* l_d_hu = nullptr;
    tsunami_lab::t_real* l_d_hv = nullptr;
    tsunami_lab::t_real* l_d_b  = nullptr;

    tsunami_lab::cuda::CudaRegressionTest::allocateGPUMemory( i_nCellsX,
                                                              i_nCellsY,
                                                              &l_d_h,
                                                              &l_d_hu,
                                                              &l_d_hv,
                                                              &l_d_b );
    // allocation succeeded if pointers are non-null
    REQUIRE( l_d_h  != nullptr );
    REQUIRE( l_d_hu != nullptr );
    REQUIRE( l_d_hv != nullptr );
    REQUIRE( l_d_b  != nullptr );  
    
    // --- copy CPU → GPU ---
    tsunami_lab::cuda::CudaRegressionTest::copyToGPU( i_nCellsX,
                                                      i_nCellsY,
                                                      l_waveProp.getHeight(),
                                                      l_waveProp.getMomentumX(),
                                                      l_waveProp.getMomentumY(),
                                                      l_waveProp.getBathymetry(),
                                                      l_d_h,
                                                      l_d_hu,
                                                      l_d_hv,
                                                      l_d_b );
    // --- copy GPU → CPU into a fresh buffer and compare ---
    tsunami_lab::t_idx l_nCells = i_nCellsX * i_nCellsY;
    tsunami_lab::t_real* l_h_back  = new tsunami_lab::t_real[l_nCells];
    tsunami_lab::t_real* l_hu_back = new tsunami_lab::t_real[l_nCells];
    tsunami_lab::t_real* l_hv_back = new tsunami_lab::t_real[l_nCells];

    tsunami_lab::cuda::CudaRegressionTest::copyFromGPU( i_nCellsX,
                                                         i_nCellsY,
                                                         l_d_h,
                                                         l_d_hu,
                                                         l_d_hv,
                                                         l_h_back,
                                                         l_hu_back,
                                                         l_hv_back );

    // round-trip check: CPU → GPU → CPU should be identical
    REQUIRE( tsunami_lab::cuda::CudaRegressionTest::compareGrids( i_nCellsX,
                                                                   i_nCellsY,
                                                                   l_waveProp.getHeight(),
                                                                   l_h_back ) );

    REQUIRE( tsunami_lab::cuda::CudaRegressionTest::compareGrids( i_nCellsX,
                                                                   i_nCellsY,
                                                                   l_waveProp.getMomentumX(),
                                                                   l_hu_back ) );

    REQUIRE( tsunami_lab::cuda::CudaRegressionTest::compareGrids( i_nCellsX,
                                                                   i_nCellsY,
                                                                   l_waveProp.getMomentumY(),
                                                                   l_hv_back ) );

    // --- free ---
    delete[] l_h_back;
    delete[] l_hu_back;
    delete[] l_hv_back;

    tsunami_lab::cuda::CudaRegressionTest::freeGPUMemory( l_d_h,
                                                           l_d_hu,
                                                           l_d_hv,
                                                           l_d_b );
}

TEST_CASE( "CUDA regression infrastructure: 500x500 grid", "[CudaRegression500]" ) {
    runGridTest( 500, 500 );
}

TEST_CASE( "CUDA regression infrastructure: 1000x1000 grid", "[CudaRegression1000]" ) {
    runGridTest( 1000, 1000 );
}

TEST_CASE( "CUDA regression infrastructure: 4000x4000 grid", "[CudaRegression4000]" ) {
    runGridTest( 4000, 4000 );
}
