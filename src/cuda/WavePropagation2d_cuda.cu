/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * GPU wrapper implementation for two-dimensional wave propagation solver.
 **/

#include "WavePropagation2d_cuda.h"
#include <cuda_runtime.h>
#include <cstring>

// Forward declarations of kernel functions
namespace tsunami_lab {
  namespace patches {
    namespace cuda {
      __global__
      void xSweep( const t_real *i_h_old,
                   const t_real *i_hu_old,
                   const t_real *i_hv_old,
                   const t_real *i_b,
                   t_real *o_h_new,
                   t_real *o_hu_new,
                   t_idx i_nCellsX,
                   t_idx i_nCellsY,
                   t_idx i_stride,
                   t_real i_scaling,
                   bool i_useFWave );

      __global__
      void ySweep( const t_real *i_h_old,
                   const t_real *i_hu_old,
                   const t_real *i_hv_old,
                   const t_real *i_b,
                   t_real *o_h_new,
                   t_real *o_hv_new,
                   t_idx i_nCellsX,
                   t_idx i_nCellsY,
                   t_idx i_stride,
                   t_real i_scaling,
                   bool i_useFWave );

      __global__
      void initNewCells( const t_real *i_h_old,
                        const t_real *i_hu_old,
                        const t_real *i_hv_old,
                        t_real *o_h_new,
                        t_real *o_hu_new,
                        t_real *o_hv_new,
                        t_idx i_nValues );
    }
  }
}

tsunami_lab::patches::cuda::WavePropagation2dCuda::WavePropagation2dCuda(
    t_idx i_nCellsX,
    t_idx i_nCellsY,
    bool i_useFWaveSolver ) : m_nCellsX( i_nCellsX ),
                              m_nCellsY( i_nCellsY ),
                              m_useFWaveSolver( i_useFWaveSolver ) {

  // Stride includes ghost cells on both sides: nCellsX + 2
  m_stride = i_nCellsX + 2;

  // Total number of values: (nCellsY + 2) * stride
  m_nValues = (i_nCellsY + 2) * m_stride;

  // Allocate GPU memory for all arrays
  // Each array needs space for two time steps (current and next)
  size_t l_size = m_nValues * sizeof(t_real);

  for( unsigned short l_st = 0; l_st < 2; l_st++ ) {
    cudaMalloc( (void**)&m_d_h[l_st], l_size );
    cudaMalloc( (void**)&m_d_hu[l_st], l_size );
    cudaMalloc( (void**)&m_d_hv[l_st], l_size );

    // Initialize to zero
    cudaMemset( m_d_h[l_st], 0, l_size );
    cudaMemset( m_d_hu[l_st], 0, l_size );
    cudaMemset( m_d_hv[l_st], 0, l_size );
  }

  // Allocate bathymetry (single copy, constant throughout simulation)
  cudaMalloc( (void**)&m_d_b, l_size );
  cudaMemset( m_d_b, 0, l_size );
}

tsunami_lab::patches::cuda::WavePropagation2dCuda::~WavePropagation2dCuda() {
  for( unsigned short l_st = 0; l_st < 2; l_st++ ) {
    if( m_d_h[l_st] != nullptr ) {
      cudaFree( m_d_h[l_st] );
      m_d_h[l_st] = nullptr;
    }
    if( m_d_hu[l_st] != nullptr ) {
      cudaFree( m_d_hu[l_st] );
      m_d_hu[l_st] = nullptr;
    }
    if( m_d_hv[l_st] != nullptr ) {
      cudaFree( m_d_hv[l_st] );
      m_d_hv[l_st] = nullptr;
    }
  }
  if( m_d_b != nullptr ) {
    cudaFree( m_d_b );
    m_d_b = nullptr;
  }
}

void tsunami_lab::patches::cuda::WavePropagation2dCuda::copyToGpu(
    const t_real *i_h,
    const t_real *i_hu,
    const t_real *i_hv,
    const t_real *i_b ) {

  size_t l_size = m_nValues * sizeof(t_real);

  // Copy to active buffers
  cudaMemcpy( m_d_h[m_step], i_h, l_size, cudaMemcpyHostToDevice );
  cudaMemcpy( m_d_hu[m_step], i_hu, l_size, cudaMemcpyHostToDevice );
  cudaMemcpy( m_d_hv[m_step], i_hv, l_size, cudaMemcpyHostToDevice );

  // Copy bathymetry (constant)
  cudaMemcpy( m_d_b, i_b, l_size, cudaMemcpyHostToDevice );
}

void tsunami_lab::patches::cuda::WavePropagation2dCuda::copyToHost(
    t_real *o_h,
    t_real *o_hu,
    t_real *o_hv ) {

  size_t l_size = m_nValues * sizeof(t_real);

  // Copy from active buffers
  cudaMemcpy( o_h, m_d_h[m_step], l_size, cudaMemcpyDeviceToHost );
  cudaMemcpy( o_hu, m_d_hu[m_step], l_size, cudaMemcpyDeviceToHost );
  cudaMemcpy( o_hv, m_d_hv[m_step], l_size, cudaMemcpyDeviceToHost );
}

void tsunami_lab::patches::cuda::WavePropagation2dCuda::initNewCells() {
  // Get pointers to active and next buffers
  t_real *l_h_old = m_d_h[m_step];
  t_real *l_hu_old = m_d_hu[m_step];
  t_real *l_hv_old = m_d_hv[m_step];

  t_real *l_h_new = m_d_h[1 - m_step];
  t_real *l_hu_new = m_d_hu[1 - m_step];
  t_real *l_hv_new = m_d_hv[1 - m_step];

  // Launch initialization kernel
  // Use 256 threads per block (or 512 for better occupancy)
  int l_threadsPerBlock = 256;
  int l_blocksPerGrid = (m_nValues + l_threadsPerBlock - 1) / l_threadsPerBlock;

  initNewCells<<<l_blocksPerGrid, l_threadsPerBlock>>>(
      l_h_old, l_hu_old, l_hv_old,
      l_h_new, l_hu_new, l_hv_new,
      m_nValues );

  cudaDeviceSynchronize();
}

void tsunami_lab::patches::cuda::WavePropagation2dCuda::xSweep( t_real i_scaling ) {
  // Get pointers to active and next buffers
  t_real *l_h_old = m_d_h[m_step];
  t_real *l_hu_old = m_d_hu[m_step];
  t_real *l_hv_old = m_d_hv[m_step];

  t_real *l_h_new = m_d_h[1 - m_step];
  t_real *l_hu_new = m_d_hu[1 - m_step];

  // Configure grid and block dimensions
  // Each edge in x-sweep has coordinates (cx, cy) where:
  //   cx ranges from 0 to nCellsX (nCellsX+1 edges per row)
  //   cy ranges from 1 to nCellsY (interior rows only)

  dim3 l_blockDim( 16, 16 );  // 256 threads per block
  dim3 l_gridDim(
      (m_nCellsX + 1 + l_blockDim.x - 1) / l_blockDim.x,
      (m_nCellsY + l_blockDim.y - 1) / l_blockDim.y
  );

  xSweep<<<l_gridDim, l_blockDim>>>(
      l_h_old, l_hu_old, l_hv_old, m_d_b,
      l_h_new, l_hu_new,
      m_nCellsX, m_nCellsY, m_stride,
      i_scaling, m_useFWaveSolver );

  cudaDeviceSynchronize();
}

void tsunami_lab::patches::cuda::WavePropagation2dCuda::ySweep( t_real i_scaling ) {
  // Get pointers to active and next buffers
  t_real *l_h_old = m_d_h[m_step];
  t_real *l_hu_old = m_d_hu[m_step];
  t_real *l_hv_old = m_d_hv[m_step];

  t_real *l_h_new = m_d_h[1 - m_step];
  t_real *l_hv_new = m_d_hv[1 - m_step];

  // Configure grid and block dimensions
  // Each edge in y-sweep has coordinates (cx, cy) where:
  //   cx ranges from 1 to nCellsX (interior columns only)
  //   cy ranges from 0 to nCellsY (nCellsY+1 edges per column)

  dim3 l_blockDim( 16, 16 );
  dim3 l_gridDim(
      (m_nCellsX + l_blockDim.x - 1) / l_blockDim.x,
      (m_nCellsY + 1 + l_blockDim.y - 1) / l_blockDim.y
  );

  ySweep<<<l_gridDim, l_blockDim>>>(
      l_h_old, l_hu_old, l_hv_old, m_d_b,
      l_h_new, l_hv_new,
      m_nCellsX, m_nCellsY, m_stride,
      i_scaling, m_useFWaveSolver );

  cudaDeviceSynchronize();
}

void tsunami_lab::patches::cuda::WavePropagation2dCuda::swapBuffers() {
  m_step = 1 - m_step;
}
