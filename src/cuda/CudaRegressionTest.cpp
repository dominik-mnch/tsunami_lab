#include "CudaRegressionTest.h"
#include <cuda_runtime.h>
#include <cmath>
#include <iostream>
#include "../patches/WavePropagation2d.h"
#include "../patches/WavePropagation2d_cuda.h"
#include "../setups/CircularDamBreak2d/CircularDamBreak2d.h"

namespace tsunami_lab {
namespace cuda {

/**
 * @brief Calculate maximum absolute error across grid cells.
 */
t_real CudaRegressionTest::getMaxError( t_idx i_nCellsX,
                                        t_idx i_nCellsY,
                                        t_real const* i_cpu,
                                        t_real const* i_gpu ) {
    t_real l_maxError = 0.0;
    for( t_idx l_y = 0; l_y < i_nCellsY; l_y++ ) {
        for( t_idx l_x = 0; l_x < i_nCellsX; l_x++ ) {
            t_idx l_idx = l_y * i_nCellsX + l_x;
            t_real l_diff = std::fabs( i_cpu[l_idx] - i_gpu[l_idx] );
            if( l_diff > l_maxError ) {
                l_maxError = l_diff;
            }
        }
    }
    return l_maxError;
}

/**
 * @brief Compare GPU and CPU solver results after one timestep.
 */
bool CudaRegressionTest::compareKernelResults( t_idx i_nCellsX,
                                               t_idx i_nCellsY,
                                               bool i_useFWave,
                                               t_real* o_maxError ) {
    // --- Setup ---
    patches::WavePropagation2d l_wavePropCpu( i_nCellsX, i_nCellsY, i_useFWave );
    patches::WavePropagation2dCuda l_wavePropGpu( i_nCellsX, i_nCellsY, i_useFWave );

    setups::CircularDamBreak2d l_setup( 10.0, 5.0, 0, 0, 0, 0, 0, 0, 10 );

    t_real l_domainSize = 50.0;
    t_real l_dx = l_domainSize / i_nCellsX;
    t_real l_dy = l_domainSize / i_nCellsY;

    // Initialize both solvers with identical data
    for( t_idx l_iy = 0; l_iy < i_nCellsY; l_iy++ ) {
        for( t_idx l_ix = 0; l_ix < i_nCellsX; l_ix++ ) {
            t_real l_x = -l_domainSize / 2.0 + (l_ix + 0.5) * l_dx;
            t_real l_y = -l_domainSize / 2.0 + (l_iy + 0.5) * l_dy;

            t_real l_h  = l_setup.getHeight( l_x, l_y );
            t_real l_hu = l_setup.getMomentumX( l_x, l_y );
            t_real l_hv = l_setup.getMomentumY( l_x, l_y );
            t_real l_b  = l_setup.getBathymetry( l_x, l_y );

            l_wavePropCpu.setHeight( l_ix, l_iy, l_h );
            l_wavePropCpu.setMomentumX( l_ix, l_iy, l_hu );
            l_wavePropCpu.setMomentumY( l_ix, l_iy, l_hv );
            l_wavePropCpu.setBathymetry( l_ix, l_iy, l_b );
        }
    }

    // Copy to GPU
    l_wavePropGpu.copyToGpu( l_wavePropCpu.getHeight(),
                             l_wavePropCpu.getMomentumX(),
                             l_wavePropCpu.getMomentumY(),
                             l_wavePropCpu.getBathymetry() );

    // --- Run one timestep on both ---
    t_real l_scaling = 0.4;  // CFL scaling factor

    // GPU step
    l_wavePropGpu.initNewCells();
    l_wavePropGpu.xSweep( l_scaling );
    l_wavePropGpu.ySweep( l_scaling );
    l_wavePropGpu.swapBuffers();

    // CPU step (equivalent operations)
    l_wavePropCpu.setGhostOutflow();
    l_wavePropCpu.computeNetUpdates( l_scaling );

    // Copy GPU results back to host
    t_real* l_h_gpu_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hu_gpu_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hv_gpu_result = new t_real[i_nCellsX * i_nCellsY];

    l_wavePropGpu.copyToHost( l_h_gpu_result,
                              l_hu_gpu_result,
                              l_hv_gpu_result );

    // --- Compare ---
    bool l_hMatch = compareGrids( i_nCellsX, i_nCellsY,
                                   l_wavePropCpu.getHeight(), l_h_gpu_result );
    bool l_huMatch = compareGrids( i_nCellsX, i_nCellsY,
                                    l_wavePropCpu.getMomentumX(), l_hu_gpu_result );
    bool l_hvMatch = compareGrids( i_nCellsX, i_nCellsY,
                                    l_wavePropCpu.getMomentumY(), l_hv_gpu_result );

    if( o_maxError != nullptr ) {
        t_real l_maxH = getMaxError( i_nCellsX, i_nCellsY,
                                      l_wavePropCpu.getHeight(), l_h_gpu_result );
        t_real l_maxHu = getMaxError( i_nCellsX, i_nCellsY,
                                       l_wavePropCpu.getMomentumX(), l_hu_gpu_result );
        t_real l_maxHv = getMaxError( i_nCellsX, i_nCellsY,
                                       l_wavePropCpu.getMomentumY(), l_hv_gpu_result );
        *o_maxError = std::max( { l_maxH, l_maxHu, l_maxHv } );
    }

    // Cleanup
    delete[] l_h_gpu_result;
    delete[] l_hu_gpu_result;
    delete[] l_hv_gpu_result;

    return l_hMatch && l_huMatch && l_hvMatch;
}

/**
 * @brief Compare GPU and CPU solver results over multiple timesteps.
 */
bool CudaRegressionTest::compareMultipleTimesteps( t_idx i_nCellsX,
                                                   t_idx i_nCellsY,
                                                   t_idx i_nTimesteps,
                                                   bool i_useFWave,
                                                   t_idx i_checkInterval,
                                                   t_real* o_maxError ) {
    // --- Setup ---
    patches::WavePropagation2d l_wavePropCpu( i_nCellsX, i_nCellsY, i_useFWave );
    patches::WavePropagation2dCuda l_wavePropGpu( i_nCellsX, i_nCellsY, i_useFWave );

    setups::CircularDamBreak2d l_setup( 10.0, 5.0, 0, 0, 0, 0, 0, 0, 10 );

    t_real l_domainSize = 50.0;
    t_real l_dx = l_domainSize / i_nCellsX;
    t_real l_dy = l_domainSize / i_nCellsY;

    // Initialize both solvers with identical data
    for( t_idx l_iy = 0; l_iy < i_nCellsY; l_iy++ ) {
        for( t_idx l_ix = 0; l_ix < i_nCellsX; l_ix++ ) {
            t_real l_x = -l_domainSize / 2.0 + (l_ix + 0.5) * l_dx;
            t_real l_y = -l_domainSize / 2.0 + (l_iy + 0.5) * l_dy;

            t_real l_h  = l_setup.getHeight( l_x, l_y );
            t_real l_hu = l_setup.getMomentumX( l_x, l_y );
            t_real l_hv = l_setup.getMomentumY( l_x, l_y );
            t_real l_b  = l_setup.getBathymetry( l_x, l_y );

            l_wavePropCpu.setHeight( l_ix, l_iy, l_h );
            l_wavePropCpu.setMomentumX( l_ix, l_iy, l_hu );
            l_wavePropCpu.setMomentumY( l_ix, l_iy, l_hv );
            l_wavePropCpu.setBathymetry( l_ix, l_iy, l_b );
        }
    }

    // Copy to GPU
    l_wavePropGpu.copyToGpu( l_wavePropCpu.getHeight(),
                             l_wavePropCpu.getMomentumX(),
                             l_wavePropCpu.getMomentumY(),
                             l_wavePropCpu.getBathymetry() );

    // --- Time stepping ---
    t_real l_scaling = 0.4;
    bool l_allMatch = true;
    t_real l_maxErrorAccum = 0.0;
    t_real* l_h_gpu_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hu_gpu_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hv_gpu_result = new t_real[i_nCellsX * i_nCellsY];

    for( t_idx l_step = 0; l_step < i_nTimesteps; l_step++ ) {
        // GPU step
        l_wavePropGpu.initNewCells();
        l_wavePropGpu.xSweep( l_scaling );
        l_wavePropGpu.ySweep( l_scaling );
        l_wavePropGpu.swapBuffers();

        // CPU step
        l_wavePropCpu.setGhostOutflow();
        l_wavePropCpu.computeNetUpdates( l_scaling );

        // Check at intervals
        if( (l_step + 1) % i_checkInterval == 0 ) {
            l_wavePropGpu.copyToHost( l_h_gpu_result,
                                      l_hu_gpu_result,
                                      l_hv_gpu_result );

            bool l_hMatch = compareGrids( i_nCellsX, i_nCellsY,
                                           l_wavePropCpu.getHeight(), l_h_gpu_result,
                                           TOLERANCE * 10.0 );  // Relaxed tolerance for multi-step
            bool l_huMatch = compareGrids( i_nCellsX, i_nCellsY,
                                            l_wavePropCpu.getMomentumX(), l_hu_gpu_result,
                                            TOLERANCE * 10.0 );
            bool l_hvMatch = compareGrids( i_nCellsX, i_nCellsY,
                                            l_wavePropCpu.getMomentumY(), l_hv_gpu_result,
                                            TOLERANCE * 10.0 );

            if( !(l_hMatch && l_huMatch && l_hvMatch) ) {
                std::cout << "Mismatch at timestep " << (l_step + 1) << std::endl;
                l_allMatch = false;
            }

            t_real l_maxH = getMaxError( i_nCellsX, i_nCellsY,
                                         l_wavePropCpu.getHeight(), l_h_gpu_result );
            t_real l_maxHu = getMaxError( i_nCellsX, i_nCellsY,
                                          l_wavePropCpu.getMomentumX(), l_hu_gpu_result );
            t_real l_maxHv = getMaxError( i_nCellsX, i_nCellsY,
                                          l_wavePropCpu.getMomentumY(), l_hv_gpu_result );
            t_real l_stepMaxError = std::max( { l_maxH, l_maxHu, l_maxHv } );
            if( l_stepMaxError > l_maxErrorAccum ) {
                l_maxErrorAccum = l_stepMaxError;
            }
        }
    }

    if( o_maxError != nullptr ) {
        *o_maxError = l_maxErrorAccum;
    }

    // Cleanup
    delete[] l_h_gpu_result;
    delete[] l_hu_gpu_result;
    delete[] l_hv_gpu_result;

    return l_allMatch;
}

/**
 * @brief Release device memory allocated for regression test arrays.
 *
 * Frees the device buffers allocated for height, x-momentum, y-momentum, and bathymetry.
 */
void CudaRegressionTest::freeGPUMemory( t_real* i_h,
                                        t_real* i_hu,
                                        t_real* i_hv,
                                        t_real* i_b ) {
    cudaFree( i_h  );
    cudaFree( i_hu );
    cudaFree( i_hv );
    cudaFree( i_b  );
}

} // namespace cuda
} // namespace tsunami_lab