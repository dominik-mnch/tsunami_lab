#include "CudaRegressionTest.h"
#include <cuda_runtime.h>
#include <cmath>
#include <iostream>

namespace tsunami_lab {
namespace cuda {

void CudaRegressionTest::allocateGPUMemory( int    i_nCellsX,
                                            int    i_nCellsY,
                                            float** o_h,
                                            float** o_hu,
                                            float** o_hv,
                                            float** o_b ) {
    int l_size = i_nCellsX * i_nCellsY * sizeof(float);
    cudaMalloc( (void**)o_h,  l_size );
    cudaMalloc( (void**)o_hu, l_size );
    cudaMalloc( (void**)o_hv, l_size );
    cudaMalloc( (void**)o_b,  l_size );
}

void CudaRegressionTest::copyToGPU( int          i_nCellsX,
                                    int          i_nCellsY,
                                    float const* i_h,
                                    float const* i_hu,
                                    float const* i_hv,
                                    float const* i_b,
                                    float*       o_h_gpu,
                                    float*       o_hu_gpu,
                                    float*       o_hv_gpu,
                                    float*       o_b_gpu ) {
    int l_size = i_nCellsX * i_nCellsY * sizeof(float);
    cudaMemcpy( o_h_gpu,  i_h,  l_size, cudaMemcpyHostToDevice );
    cudaMemcpy( o_hu_gpu, i_hu, l_size, cudaMemcpyHostToDevice );
    cudaMemcpy( o_hv_gpu, i_hv, l_size, cudaMemcpyHostToDevice );
    cudaMemcpy( o_b_gpu,  i_b,  l_size, cudaMemcpyHostToDevice );
}

void CudaRegressionTest::copyFromGPU( int          i_nCellsX,
                                      int          i_nCellsY,
                                      float const* i_h_gpu,
                                      float const* i_hu_gpu,
                                      float const* i_hv_gpu,
                                      float const* i_b_gpu,
                                      float*       o_h,
                                      float*       o_hu,
                                      float*       o_hv,
                                      float*       o_b ) {
    int l_size = i_nCellsX * i_nCellsY * sizeof(float);
    cudaMemcpy( o_h,  i_h_gpu,  l_size, cudaMemcpyDeviceToHost );
    cudaMemcpy( o_hu, i_hu_gpu, l_size, cudaMemcpyDeviceToHost );
    cudaMemcpy( o_hv, i_hv_gpu, l_size, cudaMemcpyDeviceToHost );
    cudaMemcpy( o_b,  i_b_gpu,  l_size, cudaMemcpyDeviceToHost );
}

//generic: h, hu, hv, b
bool CudaRegressionTest::compareGrids( int          i_nCellsX,
                                       int          i_nCellsY,
                                       float const* i_cpu,
                                       float const* i_gpu,
                                       float        i_tolerance ) {
    bool l_match = true;
    for( int l_y = 0; l_y < i_nCellsY; l_y++ ) {
        for( int l_x = 0; l_x < i_nCellsX; l_x++ ) {
            int   l_idx  = l_y * i_nCellsX + l_x;
            float l_diff = std::fabs( i_cpu[l_idx] - i_gpu[l_idx] );
            if( l_diff > i_tolerance ) {
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

void CudaRegressionTest::freeGPUMemory( float* i_h,
                                        float* i_hu,
                                        float* i_hv,
                                        float* i_b ) {
    cudaFree( i_h  );
    cudaFree( i_hu );
    cudaFree( i_hv );
    cudaFree( i_b  );
}

} // namespace cuda
} // namespace tsunami_lab