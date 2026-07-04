#include "CudaRegressionTest.h"
#include <cuda_runtime.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "../patches/WavePropagation2d.h"
#include "WavePropagation2d_cuda.h"
#include "../setups/CircularDamBreak2d/CircularDamBreak2d.h"

namespace tsunami_lab {
namespace cuda {

namespace {
/**
 * @brief Copy the interior cells of a ghost-celled, strided grid into a flat
 *        contiguous array of size nCellsX * nCellsY.
 *
 * @param i_nCellsX number of interior cells in x-direction.
 * @param i_nCellsY number of interior cells in y-direction.
 * @param i_stride  row stride of the source grid (nCellsX + 2).
 * @param i_interior pointer to interior cell (0,0) of the source grid.
 * @param o_flat    destination flat array (row-major, stride nCellsX).
 */
void flattenInterior( t_idx i_nCellsX,
                      t_idx i_nCellsY,
                      t_idx i_stride,
                      t_real const* i_interior,
                      t_real* o_flat ) {
    for( t_idx l_iy = 0; l_iy < i_nCellsY; l_iy++ ) {
        for( t_idx l_ix = 0; l_ix < i_nCellsX; l_ix++ ) {
            o_flat[l_iy * i_nCellsX + l_ix] = i_interior[l_iy * i_stride + l_ix];
        }
    }
}

void flattenInteriorColumnMajor( t_idx i_nCellsX,
                                 t_idx i_nCellsY,
                                 t_idx i_strideY,
                                 t_real const* i_interior,
                                 t_real* o_flat ) {
    for( t_idx l_iy = 0; l_iy < i_nCellsY; l_iy++ ) {
        for( t_idx l_ix = 0; l_ix < i_nCellsX; l_ix++ ) {
            // column-major source index: col * strideY + row
            o_flat[l_iy * i_nCellsX + l_ix] = i_interior[l_ix * i_strideY + l_iy];
        }
    }
}

} // anonymous namespace

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
    patches::cuda::WavePropagation2dCuda l_wavePropGpu( i_nCellsX, i_nCellsY, i_useFWave );

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

    // Apply outflow boundary on the CPU BEFORE copying so the GPU starts from
    // an identical full (ghost-celled) grid. A single GPU step then only relies
    // on ghost values that were fixed here, so every interior cell must match.
    l_wavePropCpu.setGhostOutflow();

    t_idx l_stride = l_wavePropCpu.getStride();        // nCellsX + 2
    t_idx l_nFull  = (i_nCellsY + 2) * l_stride;

    // Base pointers to the full ghost-celled CPU buffers.
    const t_real* l_hBase  = l_wavePropCpu.getHeight()     - (l_stride + 1);
    const t_real* l_huBase = l_wavePropCpu.getMomentumX()  - (l_stride + 1);
    const t_real* l_hvBase = l_wavePropCpu.getMomentumY()  - (l_stride + 1);
    const t_real* l_bBase  = l_wavePropCpu.getBathymetry() - (l_stride + 1);

    // Copy the full grid to GPU.
    l_wavePropGpu.copyToGpu( l_hBase, l_huBase, l_hvBase, l_bBase );

    // --- Run one timestep on both ---
    t_real l_scaling = 0.4;  // CFL scaling factor

    // GPU step (single fused, race-free kernel)
    l_wavePropGpu.computeStep( l_scaling );
    l_wavePropGpu.swapBuffers();

    // CPU step (ghost cells already applied above)
    l_wavePropCpu.timeStep( l_scaling );

    // Copy the full GPU grid back to host.
    t_real* l_h_gpu_full  = new t_real[l_nFull];
    t_real* l_hu_gpu_full = new t_real[l_nFull];
    t_real* l_hv_gpu_full = new t_real[l_nFull];

    l_wavePropGpu.copyToHost( l_h_gpu_full,
                              l_hu_gpu_full,
                              l_hv_gpu_full );

    // Flatten interiors of both CPU and GPU grids for comparison.
    t_real* l_h_gpu_result  = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hu_gpu_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hv_gpu_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_h_cpu_result  = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hu_cpu_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hv_cpu_result = new t_real[i_nCellsX * i_nCellsY];

    flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_h_gpu_full  + l_stride + 1, l_h_gpu_result );
    flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_hu_gpu_full + l_stride + 1, l_hu_gpu_result );
    flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_hv_gpu_full + l_stride + 1, l_hv_gpu_result );
    flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_wavePropCpu.getHeight(),    l_h_cpu_result );
    flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_wavePropCpu.getMomentumX(), l_hu_cpu_result );
    flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_wavePropCpu.getMomentumY(), l_hv_cpu_result );

    // --- Compare ---
    bool l_hMatch = compareGrids( i_nCellsX, i_nCellsY,
                                   l_h_cpu_result, l_h_gpu_result );
    bool l_huMatch = compareGrids( i_nCellsX, i_nCellsY,
                                    l_hu_cpu_result, l_hu_gpu_result );
    bool l_hvMatch = compareGrids( i_nCellsX, i_nCellsY,
                                    l_hv_cpu_result, l_hv_gpu_result );

    if( o_maxError != nullptr ) {
        t_real l_maxH = getMaxError( i_nCellsX, i_nCellsY,
                                      l_h_cpu_result, l_h_gpu_result );
        t_real l_maxHu = getMaxError( i_nCellsX, i_nCellsY,
                                       l_hu_cpu_result, l_hu_gpu_result );
        t_real l_maxHv = getMaxError( i_nCellsX, i_nCellsY,
                                       l_hv_cpu_result, l_hv_gpu_result );
        *o_maxError = std::max( { l_maxH, l_maxHu, l_maxHv } );
    }

    // Cleanup
    delete[] l_h_gpu_full;
    delete[] l_hu_gpu_full;
    delete[] l_hv_gpu_full;
    delete[] l_h_gpu_result;
    delete[] l_hu_gpu_result;
    delete[] l_hv_gpu_result;
    delete[] l_h_cpu_result;
    delete[] l_hu_cpu_result;
    delete[] l_hv_cpu_result;

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
    patches::cuda::WavePropagation2dCuda l_wavePropGpu( i_nCellsX, i_nCellsY, i_useFWave );

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

    // Apply outflow boundary on the CPU before copying so the GPU starts from
    // an identical full (ghost-celled) grid.
    l_wavePropCpu.setGhostOutflow();

    t_idx l_stride = l_wavePropCpu.getStride();   // nCellsX + 2
    t_idx l_nFull  = (i_nCellsY + 2) * l_stride;

    // Copy the full ghost-celled grid to the GPU.
    l_wavePropGpu.copyToGpu( l_wavePropCpu.getHeight()     - (l_stride + 1),
                             l_wavePropCpu.getMomentumX()  - (l_stride + 1),
                             l_wavePropCpu.getMomentumY()  - (l_stride + 1),
                             l_wavePropCpu.getBathymetry() - (l_stride + 1) );

    // --- Time stepping ---
    // Use a CFL-stable scaling (dt/dx). The maximum wave speed for this setup is
    // ~sqrt(g*h) with h up to 10, i.e. ~10, so dt/dx must stay below ~0.1 to keep
    // the explicit scheme stable. A larger value makes BOTH the CPU and GPU
    // simulations diverge to non-physical magnitudes, which then cannot match.
    t_real l_scaling = 0.04;
    bool l_allMatch = true;
    t_real l_maxErrorAccum = 0.0;

    // Full-grid GPU result buffers (sized for the ghost-celled grid).
    t_real* l_h_gpu_full  = new t_real[l_nFull];
    t_real* l_hu_gpu_full = new t_real[l_nFull];
    t_real* l_hv_gpu_full = new t_real[l_nFull];

    // Flattened interior buffers for comparison.
    t_real* l_h_gpu_result  = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hu_gpu_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hv_gpu_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_h_cpu_result  = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hu_cpu_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hv_cpu_result = new t_real[i_nCellsX * i_nCellsY];

    for( t_idx l_step = 0; l_step < i_nTimesteps; l_step++ ) {
        // GPU step
        l_wavePropGpu.setGhostOutflow();
        l_wavePropGpu.computeStep( l_scaling );
        l_wavePropGpu.swapBuffers();

        // CPU step
        l_wavePropCpu.setGhostOutflow();
        l_wavePropCpu.timeStep( l_scaling );

        // Check at intervals
        if( (l_step + 1) % i_checkInterval == 0 ) {
            l_wavePropGpu.copyToHost( l_h_gpu_full,
                                      l_hu_gpu_full,
                                      l_hv_gpu_full );

            // Flatten interiors of both CPU and GPU grids for comparison.
            flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_h_gpu_full  + l_stride + 1, l_h_gpu_result );
            flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_hu_gpu_full + l_stride + 1, l_hu_gpu_result );
            flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_hv_gpu_full + l_stride + 1, l_hv_gpu_result );
            flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_wavePropCpu.getHeight(),    l_h_cpu_result );
            flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_wavePropCpu.getMomentumX(), l_hu_cpu_result );
            flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_wavePropCpu.getMomentumY(), l_hv_cpu_result );

            // Relative tolerance for multi-step comparison. Both CPU and GPU use
            // the identical per-edge solver and the same accumulation order, but
            // floating-point round-off still differs slightly and is amplified
            // near shocks/steep gradients over many explicit time steps, so CPU
            // and GPU differ by a small relative amount (well below 2%).
            // compareGrids scales this tolerance by magnitude.
            const t_real l_multiStepTol = static_cast<t_real>(0.02);
            bool l_hMatch = compareGrids( i_nCellsX, i_nCellsY,
                                           l_h_cpu_result, l_h_gpu_result,
                                           l_multiStepTol );
            bool l_huMatch = compareGrids( i_nCellsX, i_nCellsY,
                                            l_hu_cpu_result, l_hu_gpu_result,
                                            l_multiStepTol );
            bool l_hvMatch = compareGrids( i_nCellsX, i_nCellsY,
                                            l_hv_cpu_result, l_hv_gpu_result,
                                            l_multiStepTol );

            if( !(l_hMatch && l_huMatch && l_hvMatch) ) {
                std::cout << "Mismatch at timestep " << (l_step + 1) << std::endl;
                l_allMatch = false;
            }

            t_real l_maxH = getMaxError( i_nCellsX, i_nCellsY,
                                         l_h_cpu_result, l_h_gpu_result );
            t_real l_maxHu = getMaxError( i_nCellsX, i_nCellsY,
                                          l_hu_cpu_result, l_hu_gpu_result );
            t_real l_maxHv = getMaxError( i_nCellsX, i_nCellsY,
                                          l_hv_cpu_result, l_hv_gpu_result );
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
    delete[] l_h_gpu_full;
    delete[] l_hu_gpu_full;
    delete[] l_hv_gpu_full;
    delete[] l_h_gpu_result;
    delete[] l_hu_gpu_result;
    delete[] l_hv_gpu_result;
    delete[] l_h_cpu_result;
    delete[] l_hu_cpu_result;
    delete[] l_hv_cpu_result;

    return l_allMatch;
}

/**
 * @brief Compare GPU and CPU solver results after one timestep (atomic kernels).
 */
bool CudaRegressionTest::compareKernelResultsAtomic( t_idx i_nCellsX,
                                                     t_idx i_nCellsY,
                                                     bool i_useFWave,
                                                     t_real* o_maxError ) {
    // --- Setup ---
    patches::WavePropagation2d l_wavePropCpu( i_nCellsX, i_nCellsY, i_useFWave );
    patches::cuda::WavePropagation2dCuda l_wavePropGpu( i_nCellsX, i_nCellsY, i_useFWave );

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

    // Apply outflow boundary on the CPU BEFORE copying so the GPU starts from
    // an identical full (ghost-celled) grid. A single GPU step then only relies
    // on ghost values that were fixed here, so every interior cell must match.
    l_wavePropCpu.setGhostOutflow();

    t_idx l_stride = l_wavePropCpu.getStride();        // nCellsX + 2
    t_idx l_nFull  = (i_nCellsY + 2) * l_stride;

    // Base pointers to the full ghost-celled CPU buffers.
    const t_real* l_hBase  = l_wavePropCpu.getHeight()     - (l_stride + 1);
    const t_real* l_huBase = l_wavePropCpu.getMomentumX()  - (l_stride + 1);
    const t_real* l_hvBase = l_wavePropCpu.getMomentumY()  - (l_stride + 1);
    const t_real* l_bBase  = l_wavePropCpu.getBathymetry() - (l_stride + 1);

    // Copy the full grid to GPU.
    l_wavePropGpu.copyToGpu( l_hBase, l_huBase, l_hvBase, l_bBase );

    // --- Run one timestep on both ---
    t_real l_scaling = 0.4;  // CFL scaling factor

    // GPU step (atomic kernels)
    l_wavePropGpu.computeStepAtomic( l_scaling );
    l_wavePropGpu.swapBuffers();

    // CPU step (ghost cells already applied above)
    l_wavePropCpu.timeStep( l_scaling );

    // Copy the full GPU grid back to host.
    t_real* l_h_gpu_full  = new t_real[l_nFull];
    t_real* l_hu_gpu_full = new t_real[l_nFull];
    t_real* l_hv_gpu_full = new t_real[l_nFull];

    l_wavePropGpu.copyToHost( l_h_gpu_full,
                              l_hu_gpu_full,
                              l_hv_gpu_full );

    // Flatten interiors of both CPU and GPU grids for comparison.
    t_real* l_h_gpu_result  = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hu_gpu_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hv_gpu_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_h_cpu_result  = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hu_cpu_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hv_cpu_result = new t_real[i_nCellsX * i_nCellsY];

    flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_h_gpu_full  + l_stride + 1, l_h_gpu_result );
    flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_hu_gpu_full + l_stride + 1, l_hu_gpu_result );
    flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_hv_gpu_full + l_stride + 1, l_hv_gpu_result );
    flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_wavePropCpu.getHeight(),    l_h_cpu_result );
    flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_wavePropCpu.getMomentumX(), l_hu_cpu_result );
    flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_wavePropCpu.getMomentumY(), l_hv_cpu_result );

    // --- Compare ---
    bool l_hMatch = compareGrids( i_nCellsX, i_nCellsY,
                                   l_h_cpu_result, l_h_gpu_result );
    bool l_huMatch = compareGrids( i_nCellsX, i_nCellsY,
                                    l_hu_cpu_result, l_hu_gpu_result );
    bool l_hvMatch = compareGrids( i_nCellsX, i_nCellsY,
                                    l_hv_cpu_result, l_hv_gpu_result );

    if( o_maxError != nullptr ) {
        t_real l_maxH = getMaxError( i_nCellsX, i_nCellsY,
                                      l_h_cpu_result, l_h_gpu_result );
        t_real l_maxHu = getMaxError( i_nCellsX, i_nCellsY,
                                       l_hu_cpu_result, l_hu_gpu_result );
        t_real l_maxHv = getMaxError( i_nCellsX, i_nCellsY,
                                       l_hv_cpu_result, l_hv_gpu_result );
        *o_maxError = std::max( { l_maxH, l_maxHu, l_maxHv } );
    }

    // Cleanup
    delete[] l_h_gpu_full;
    delete[] l_hu_gpu_full;
    delete[] l_hv_gpu_full;
    delete[] l_h_gpu_result;
    delete[] l_hu_gpu_result;
    delete[] l_hv_gpu_result;
    delete[] l_h_cpu_result;
    delete[] l_hu_cpu_result;
    delete[] l_hv_cpu_result;

    return l_hMatch && l_huMatch && l_hvMatch;
}

/**
 * @brief Compare GPU and CPU solver results over multiple timesteps (atomic kernels).
 */
bool CudaRegressionTest::compareMultipleTimestepsAtomic( t_idx i_nCellsX,
                                                         t_idx i_nCellsY,
                                                         t_idx i_nTimesteps,
                                                         bool i_useFWave,
                                                         t_idx i_checkInterval,
                                                         t_real* o_maxError ) {
    // --- Setup ---
    patches::WavePropagation2d l_wavePropCpu( i_nCellsX, i_nCellsY, i_useFWave );
    patches::cuda::WavePropagation2dCuda l_wavePropGpu( i_nCellsX, i_nCellsY, i_useFWave );

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

    // Apply outflow boundary on the CPU before copying so the GPU starts from
    // an identical full (ghost-celled) grid.
    l_wavePropCpu.setGhostOutflow();

    t_idx l_stride = l_wavePropCpu.getStride();   // nCellsX + 2
    t_idx l_nFull  = (i_nCellsY + 2) * l_stride;

    // Copy the full ghost-celled grid to the GPU.
    l_wavePropGpu.copyToGpu( l_wavePropCpu.getHeight()     - (l_stride + 1),
                             l_wavePropCpu.getMomentumX()  - (l_stride + 1),
                             l_wavePropCpu.getMomentumY()  - (l_stride + 1),
                             l_wavePropCpu.getBathymetry() - (l_stride + 1) );

    // --- Time stepping ---
    // Use a CFL-stable scaling (dt/dx). The maximum wave speed for this setup is
    // ~sqrt(g*h) with h up to 10, i.e. ~10, so dt/dx must stay below ~0.1 to keep
    // the explicit scheme stable. A larger value makes BOTH the CPU and GPU
    // simulations diverge to non-physical magnitudes, which then cannot match.
    t_real l_scaling = 0.04;
    bool l_allMatch = true;
    t_real l_maxErrorAccum = 0.0;

    // Full-grid GPU result buffers (sized for the ghost-celled grid).
    t_real* l_h_gpu_full  = new t_real[l_nFull];
    t_real* l_hu_gpu_full = new t_real[l_nFull];
    t_real* l_hv_gpu_full = new t_real[l_nFull];

    // Flattened interior buffers for comparison.
    t_real* l_h_gpu_result  = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hu_gpu_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hv_gpu_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_h_cpu_result  = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hu_cpu_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hv_cpu_result = new t_real[i_nCellsX * i_nCellsY];

    for( t_idx l_step = 0; l_step < i_nTimesteps; l_step++ ) {
        // GPU step (atomic kernels)
        l_wavePropGpu.setGhostOutflow();
        l_wavePropGpu.computeStepAtomic( l_scaling );
        l_wavePropGpu.swapBuffers();

        // CPU step
        l_wavePropCpu.setGhostOutflow();
        l_wavePropCpu.timeStep( l_scaling );

        // Check at intervals
        if( (l_step + 1) % i_checkInterval == 0 ) {
            l_wavePropGpu.copyToHost( l_h_gpu_full,
                                      l_hu_gpu_full,
                                      l_hv_gpu_full );

            // Flatten interiors of both CPU and GPU grids for comparison.
            flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_h_gpu_full  + l_stride + 1, l_h_gpu_result );
            flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_hu_gpu_full + l_stride + 1, l_hu_gpu_result );
            flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_hv_gpu_full + l_stride + 1, l_hv_gpu_result );
            flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_wavePropCpu.getHeight(),    l_h_cpu_result );
            flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_wavePropCpu.getMomentumX(), l_hu_cpu_result );
            flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_wavePropCpu.getMomentumY(), l_hv_cpu_result );

            // Relative tolerance for multi-step comparison. Both CPU and GPU use
            // the identical per-edge solver and the same accumulation order, but
            // floating-point round-off still differs slightly and is amplified
            // near shocks/steep gradients over many explicit time steps, so CPU
            // and GPU differ by a small relative amount (well below 2%).
            // compareGrids scales this tolerance by magnitude.
            const t_real l_multiStepTol = static_cast<t_real>(0.02);
            bool l_hMatch = compareGrids( i_nCellsX, i_nCellsY,
                                           l_h_cpu_result, l_h_gpu_result,
                                           l_multiStepTol );
            bool l_huMatch = compareGrids( i_nCellsX, i_nCellsY,
                                            l_hu_cpu_result, l_hu_gpu_result,
                                            l_multiStepTol );
            bool l_hvMatch = compareGrids( i_nCellsX, i_nCellsY,
                                            l_hv_cpu_result, l_hv_gpu_result,
                                            l_multiStepTol );

            if( !(l_hMatch && l_huMatch && l_hvMatch) ) {
                std::cout << "Mismatch at timestep " << (l_step + 1) << std::endl;
                l_allMatch = false;
            }

            t_real l_maxH = getMaxError( i_nCellsX, i_nCellsY,
                                         l_h_cpu_result, l_h_gpu_result );
            t_real l_maxHu = getMaxError( i_nCellsX, i_nCellsY,
                                          l_hu_cpu_result, l_hu_gpu_result );
            t_real l_maxHv = getMaxError( i_nCellsX, i_nCellsY,
                                          l_hv_cpu_result, l_hv_gpu_result );
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
    delete[] l_h_gpu_full;
    delete[] l_hu_gpu_full;
    delete[] l_hv_gpu_full;
    delete[] l_h_gpu_result;
    delete[] l_hu_gpu_result;
    delete[] l_hv_gpu_result;
    delete[] l_h_cpu_result;
    delete[] l_hu_cpu_result;
    delete[] l_hv_cpu_result;

    return l_allMatch;
}

/**
 * @brief Compare GPU row-major kernel results against GPU column-major kernel
 *        results after one timestep. Both instances use the same public API
 *        (copyToGpu/computeStep/copyToHost) — layout dispatch happens
 *        internally via m_memoryLayout. Only the host-side reinterpretation
 *        of the returned buffer differs, since copyToHost does not
 *        un-transpose column-major data (unlike Tile32, which does un-tile).
 */
bool CudaRegressionTest::compareKernelResultsColumnMajor( t_idx i_nCellsX,
                                                          t_idx i_nCellsY,
                                                          bool i_useFWave,
                                                          t_real* o_maxError ) {
    patches::WavePropagation2d l_wavePropCpu( i_nCellsX, i_nCellsY, i_useFWave );
    patches::cuda::WavePropagation2dCuda l_wavePropGpuRow( i_nCellsX, i_nCellsY, i_useFWave,
                                                           patches::cuda::MemoryLayout::RowMajor );
    patches::cuda::WavePropagation2dCuda l_wavePropGpuCol( i_nCellsX, i_nCellsY, i_useFWave,
                                                           patches::cuda::MemoryLayout::ColumnMajor );

    setups::CircularDamBreak2d l_setup( 10.0, 5.0, 0, 0, 0, 0, 0, 0, 10 );

    t_real l_domainSize = 50.0;
    t_real l_dx = l_domainSize / i_nCellsX;
    t_real l_dy = l_domainSize / i_nCellsY;

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

    // Identical, already-correct ghost cells for both GPU layouts.
    l_wavePropCpu.setGhostOutflow();

    t_idx l_stride = l_wavePropCpu.getStride();   // nCellsX + 2 (row-major)
    t_idx l_height = i_nCellsY + 2;               // nCellsY + 2 (column-major)
    t_idx l_nFull  = (i_nCellsY + 2) * l_stride;

    const t_real* l_hBase  = l_wavePropCpu.getHeight()     - (l_stride + 1);
    const t_real* l_huBase = l_wavePropCpu.getMomentumX()  - (l_stride + 1);
    const t_real* l_hvBase = l_wavePropCpu.getMomentumY()  - (l_stride + 1);
    const t_real* l_bBase  = l_wavePropCpu.getBathymetry() - (l_stride + 1);

    // Same copyToGpu for both — row-major host buffer in both cases;
    // copyToGpu internally transposes for ColumnMajor (mirroring the Tile32 branch).
    l_wavePropGpuRow.copyToGpu( l_hBase, l_huBase, l_hvBase, l_bBase );
    l_wavePropGpuCol.copyToGpu( l_hBase, l_huBase, l_hvBase, l_bBase );

    t_real l_scaling = 0.4;

    l_wavePropGpuRow.computeStep( l_scaling );
    l_wavePropGpuRow.swapBuffers();

    l_wavePropGpuCol.computeStep( l_scaling );
    l_wavePropGpuCol.swapBuffers();

    t_real* l_h_row_full  = new t_real[l_nFull];
    t_real* l_hu_row_full = new t_real[l_nFull];
    t_real* l_hv_row_full = new t_real[l_nFull];
    t_real* l_h_col_full  = new t_real[l_nFull];
    t_real* l_hu_col_full = new t_real[l_nFull];
    t_real* l_hv_col_full = new t_real[l_nFull];

    l_wavePropGpuRow.copyToHost( l_h_row_full, l_hu_row_full, l_hv_row_full );
    // NOTE: copyToHost does NOT un-transpose column-major layout — the
    // returned buffer is still laid out as col * l_height + row.
    l_wavePropGpuCol.copyToHost( l_h_col_full, l_hu_col_full, l_hv_col_full );

    t_real* l_h_row_result  = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hu_row_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hv_row_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_h_col_result  = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hu_col_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hv_col_result = new t_real[i_nCellsX * i_nCellsY];

    flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_h_row_full  + l_stride + 1, l_h_row_result );
    flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_hu_row_full + l_stride + 1, l_hu_row_result );
    flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_hv_row_full + l_stride + 1, l_hv_row_result );

    flattenInteriorColumnMajor( i_nCellsX, i_nCellsY, l_height, l_h_col_full  + l_height + 1, l_h_col_result );
    flattenInteriorColumnMajor( i_nCellsX, i_nCellsY, l_height, l_hu_col_full + l_height + 1, l_hu_col_result );
    flattenInteriorColumnMajor( i_nCellsX, i_nCellsY, l_height, l_hv_col_full + l_height + 1, l_hv_col_result );

    bool l_hMatch  = compareGrids( i_nCellsX, i_nCellsY, l_h_row_result,  l_h_col_result );
    bool l_huMatch = compareGrids( i_nCellsX, i_nCellsY, l_hu_row_result, l_hu_col_result );
    bool l_hvMatch = compareGrids( i_nCellsX, i_nCellsY, l_hv_row_result, l_hv_col_result );

    if( o_maxError != nullptr ) {
        t_real l_maxH  = getMaxError( i_nCellsX, i_nCellsY, l_h_row_result,  l_h_col_result );
        t_real l_maxHu = getMaxError( i_nCellsX, i_nCellsY, l_hu_row_result, l_hu_col_result );
        t_real l_maxHv = getMaxError( i_nCellsX, i_nCellsY, l_hv_row_result, l_hv_col_result );
        *o_maxError = std::max( { l_maxH, l_maxHu, l_maxHv } );
    }

    delete[] l_h_row_full;
    delete[] l_hu_row_full;
    delete[] l_hv_row_full;
    delete[] l_h_col_full;
    delete[] l_hu_col_full;
    delete[] l_hv_col_full;
    delete[] l_h_row_result;
    delete[] l_hu_row_result;
    delete[] l_hv_row_result;
    delete[] l_h_col_result;
    delete[] l_hu_col_result;
    delete[] l_hv_col_result;

    return l_hMatch && l_huMatch && l_hvMatch;
}

/**
 * @brief Compare GPU row-major vs GPU column-major kernel results over
 *        multiple timesteps, checking at a fixed interval.
 */
bool CudaRegressionTest::compareMultipleTimestepsColumnMajor( t_idx i_nCellsX,
                                                              t_idx i_nCellsY,
                                                              t_idx i_nTimesteps,
                                                              bool i_useFWave,
                                                              t_idx i_checkInterval,
                                                              t_real* o_maxError ) {
    patches::WavePropagation2d l_wavePropCpu( i_nCellsX, i_nCellsY, i_useFWave );
    patches::cuda::WavePropagation2dCuda l_wavePropGpuRow( i_nCellsX, i_nCellsY, i_useFWave,
                                                           patches::cuda::MemoryLayout::RowMajor );
    patches::cuda::WavePropagation2dCuda l_wavePropGpuCol( i_nCellsX, i_nCellsY, i_useFWave,
                                                           patches::cuda::MemoryLayout::ColumnMajor );

    setups::CircularDamBreak2d l_setup( 10.0, 5.0, 0, 0, 0, 0, 0, 0, 10 );

    t_real l_domainSize = 50.0;
    t_real l_dx = l_domainSize / i_nCellsX;
    t_real l_dy = l_domainSize / i_nCellsY;

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

    l_wavePropCpu.setGhostOutflow();

    t_idx l_stride = l_wavePropCpu.getStride();
    t_idx l_height = i_nCellsY + 2;
    t_idx l_nFull  = (i_nCellsY + 2) * l_stride;

    const t_real* l_hBase  = l_wavePropCpu.getHeight()     - (l_stride + 1);
    const t_real* l_huBase = l_wavePropCpu.getMomentumX()  - (l_stride + 1);
    const t_real* l_hvBase = l_wavePropCpu.getMomentumY()  - (l_stride + 1);
    const t_real* l_bBase  = l_wavePropCpu.getBathymetry() - (l_stride + 1);

    l_wavePropGpuRow.copyToGpu( l_hBase, l_huBase, l_hvBase, l_bBase );
    l_wavePropGpuCol.copyToGpu( l_hBase, l_huBase, l_hvBase, l_bBase );

    t_real l_scaling = 0.04;
    bool l_allMatch = true;
    t_real l_maxErrorAccum = 0.0;

    t_real* l_h_row_full  = new t_real[l_nFull];
    t_real* l_hu_row_full = new t_real[l_nFull];
    t_real* l_hv_row_full = new t_real[l_nFull];
    t_real* l_h_col_full  = new t_real[l_nFull];
    t_real* l_hu_col_full = new t_real[l_nFull];
    t_real* l_hv_col_full = new t_real[l_nFull];

    t_real* l_h_row_result  = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hu_row_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hv_row_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_h_col_result  = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hu_col_result = new t_real[i_nCellsX * i_nCellsY];
    t_real* l_hv_col_result = new t_real[i_nCellsX * i_nCellsY];

    for( t_idx l_step = 0; l_step < i_nTimesteps; l_step++ ) {
        l_wavePropGpuRow.setGhostOutflow();
        l_wavePropGpuRow.computeStep( l_scaling );
        l_wavePropGpuRow.swapBuffers();

        l_wavePropGpuCol.setGhostOutflow();
        l_wavePropGpuCol.computeStep( l_scaling );
        l_wavePropGpuCol.swapBuffers();

        if( (l_step + 1) % i_checkInterval == 0 ) {
            l_wavePropGpuRow.copyToHost( l_h_row_full, l_hu_row_full, l_hv_row_full );
            l_wavePropGpuCol.copyToHost( l_h_col_full, l_hu_col_full, l_hv_col_full );

            flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_h_row_full  + l_stride + 1, l_h_row_result );
            flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_hu_row_full + l_stride + 1, l_hu_row_result );
            flattenInterior( i_nCellsX, i_nCellsY, l_stride, l_hv_row_full + l_stride + 1, l_hv_row_result );

            flattenInteriorColumnMajor( i_nCellsX, i_nCellsY, l_height, l_h_col_full  + l_height + 1, l_h_col_result );
            flattenInteriorColumnMajor( i_nCellsX, i_nCellsY, l_height, l_hu_col_full + l_height + 1, l_hu_col_result );
            flattenInteriorColumnMajor( i_nCellsX, i_nCellsY, l_height, l_hv_col_full + l_height + 1, l_hv_col_result );

            const t_real l_multiStepTol = static_cast<t_real>(0.02);
            bool l_hMatch  = compareGrids( i_nCellsX, i_nCellsY, l_h_row_result,  l_h_col_result,  l_multiStepTol );
            bool l_huMatch = compareGrids( i_nCellsX, i_nCellsY, l_hu_row_result, l_hu_col_result, l_multiStepTol );
            bool l_hvMatch = compareGrids( i_nCellsX, i_nCellsY, l_hv_row_result, l_hv_col_result, l_multiStepTol );

            if( !(l_hMatch && l_huMatch && l_hvMatch) ) {
                std::cout << "Row-major/column-major mismatch at timestep " << (l_step + 1) << std::endl;
                l_allMatch = false;
            }

            t_real l_maxH  = getMaxError( i_nCellsX, i_nCellsY, l_h_row_result,  l_h_col_result );
            t_real l_maxHu = getMaxError( i_nCellsX, i_nCellsY, l_hu_row_result, l_hu_col_result );
            t_real l_maxHv = getMaxError( i_nCellsX, i_nCellsY, l_hv_row_result, l_hv_col_result );
            t_real l_stepMaxError = std::max( { l_maxH, l_maxHu, l_maxHv } );
            if( l_stepMaxError > l_maxErrorAccum ) {
                l_maxErrorAccum = l_stepMaxError;
            }
        }
    }

    if( o_maxError != nullptr ) {
        *o_maxError = l_maxErrorAccum;
    }

    delete[] l_h_row_full;
    delete[] l_hu_row_full;
    delete[] l_hv_row_full;
    delete[] l_h_col_full;
    delete[] l_hu_col_full;
    delete[] l_hv_col_full;
    delete[] l_h_row_result;
    delete[] l_hu_row_result;
    delete[] l_hv_row_result;
    delete[] l_h_col_result;
    delete[] l_hu_col_result;
    delete[] l_hv_col_result;

    return l_allMatch;
}

} // namespace cuda
} // namespace tsunami_lab