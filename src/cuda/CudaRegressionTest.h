#ifndef CUDA_REGRESSION_TEST_H
#define CUDA_REGRESSION_TEST_H

namespace tsunami_lab {
namespace cuda {

class CudaRegressionTest {
public:
    // tolerance for floating point comparison
    static constexpr float TOLERANCE = 1e-5f;

    // allocate GPU memory for grid of given size
    static void allocateGPUMemory( int i_nCellsX,
                                   int i_nCellsY,
                                   float** o_h,
                                   float** o_hu,
                                   float** o_hv,
                                   float** o_b );

    // copy data from CPU arrays to GPU arrays
    static void copyToGPU( int          i_nCellsX,
                           int          i_nCellsY,
                           float const* i_h,
                           float const* i_hu,
                           float const* i_hv,
                           float const* i_b,
                           float*       o_h_gpu,
                           float*       o_hu_gpu,
                           float*       o_hv_gpu,
                           float*       o_b_gpu );

    // copy results back from GPU to CPU
    static void copyFromGPU( int          i_nCellsX,
                             int          i_nCellsY,
                             float const* i_h_gpu,
                             float const* i_hu_gpu,
                             float const* i_hv_gpu,
                             float*       o_h,
                             float*       o_hu,
                             float*       o_hv );

    // compare two grids within tolerance, returns true if they match
    static bool compareGrids( int          i_nCellsX,
                              int          i_nCellsY,
                              float const* i_cpu,
                              float const* i_gpu,
                              float        i_tolerance = TOLERANCE );

    // free GPU memory
    static void freeGPUMemory( float* i_h,
                               float* i_hu,
                               float* i_hv,
                               float* i_b );
};

} // namespace cuda
} // namespace tsunami_lab

#endif