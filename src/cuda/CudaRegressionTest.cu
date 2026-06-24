#include "CudaRegressionTest.h"
#include <cuda_runtime.h>
#include <cmath>
#include <iostream>

namespace tsunami_lab {
namespace cuda {

void CudaRegressionTest::allocateGPUMemory( t_idx    i_nCellsX,
                                            t_idx    i_nCellsY,
                                            t_real** o_h,
                                            t_real** o_hu,
                                            t_real** o_hv,
                                            t_real** o_b ) {
    std::size_t l_size = i_nCellsX * i_nCellsY * sizeof(t_real);
    cudaMalloc( (void**)o_h,  l_size );
    cudaMalloc( (void**)o_hu, l_size );
    cudaMalloc( (void**)o_hv, l_size );
    cudaMalloc( (void**)o_b,  l_size );
}

void CudaRegressionTest::copyToGPU( t_idx          i_nCellsX,
                                    t_idx          i_nCellsY,
                                    t_real const*  i_h,
                                    t_real const*  i_hu,
                                    t_real const*  i_hv,
                                    t_real const*  i_b,
                                    t_real*        o_h_gpu,
                                    t_real*        o_hu_gpu,
                                    t_real*        o_hv_gpu,
                                    t_real*        o_b_gpu ) {
    std::size_t l_size = i_nCellsX * i_nCellsY * sizeof(t_real);
    cudaMemcpy( o_h_gpu,  i_h,  l_size, cudaMemcpyHostToDevice );
    cudaMemcpy( o_hu_gpu, i_hu, l_size, cudaMemcpyHostToDevice );
    cudaMemcpy( o_hv_gpu, i_hv, l_size, cudaMemcpyHostToDevice );
    cudaMemcpy( o_b_gpu,  i_b,  l_size, cudaMemcpyHostToDevice );
}

void CudaRegressionTest::copyFromGPU( t_idx          i_nCellsX,
                                      t_idx          i_nCellsY,
                                      t_real const*  i_h_gpu,
                                      t_real const*  i_hu_gpu,
                                      t_real const*  i_hv_gpu,
                                      t_real*        o_h,
                                      t_real*        o_hu,
                                      t_real*        o_hv ) {
    std::size_t l_size = i_nCellsX * i_nCellsY * sizeof(t_real);
    cudaMemcpy( o_h,  i_h_gpu,  l_size, cudaMemcpyDeviceToHost );
    cudaMemcpy( o_hu, i_hu_gpu, l_size, cudaMemcpyDeviceToHost );
    cudaMemcpy( o_hv, i_hv_gpu, l_size, cudaMemcpyDeviceToHost );
}

bool CudaRegressionTest::compareGrids( t_idx          i_nCellsX,
                                       t_idx          i_nCellsY,
                                       t_real const*  i_cpu,
                                       t_real const*  i_gpu,
                                       t_real         i_tolerance ) {
    bool l_match = true;
    for( t_idx l_y = 0; l_y < i_nCellsY; l_y++ ) {
        for( t_idx l_x = 0; l_x < i_nCellsX; l_x++ ) {
            t_idx  l_idx  = l_y * i_nCellsX + l_x;
            t_real l_diff = std::fabs( i_cpu[l_idx] - i_gpu[l_idx] );
            t_real l_allowed = i_tolerance * std::max( static_cast<t_real>(1),
                                                       std::fabs( i_cpu[l_idx] ) );
            if( l_diff > l_allowed ) {
                std::cerr << "MISMATCH at ("
                          << l_x << ", " << l_y << "): "
                          << "CPU=" << i_cpu[l_idx]
                          << " GPU=" << i_gpu[l_idx]
                          << " diff=" << l_diff
                          << std::endl;
                l_match = false;
            }
        }
    }
    return l_match;
}

bool CudaRegressionTest::compareSingleCells( t_idx         i_x,
                                             t_idx         i_y,
                                             t_idx         i_nCellsX,
                                             t_real const* i_cpu,
                                             t_real const* i_gpu,
                                             t_real        i_tolerance ) {
    t_idx  l_idx    = i_y * i_nCellsX + i_x;
    t_real l_gpuVal = 0;
    cudaMemcpy( &l_gpuVal, i_gpu + l_idx, sizeof(t_real), cudaMemcpyDeviceToHost );

    t_real l_diff = std::fabs( i_cpu[l_idx] - l_gpuVal );
    t_real l_allowed = i_tolerance * std::max( static_cast<t_real>(1),
                                               std::fabs( i_cpu[l_idx] ) );
    if( l_diff > l_allowed ) {
        std::cerr << "MISMATCH at ("
                  << i_x << ", " << i_y << "): "
                  << "CPU=" << i_cpu[l_idx]
                  << " GPU=" << l_gpuVal
                  << " diff=" << l_diff
                  << std::endl;
        return false;
    }
    return true;
}

void CudaRegressionTest::freeGPUMemory( t_real* i_h,
                                        t_real* i_hu,
                                        t_real* i_hv,
                                        t_real* i_b ) {
    cudaFree(i_h);
    cudaFree(i_hu);
    cudaFree(i_hv);
    cudaFree(i_b);
}

} // namespace cuda
} // namespace tsunami_lab