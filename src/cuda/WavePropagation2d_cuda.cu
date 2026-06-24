/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * GPU wrapper implementation for two-dimensional wave propagation solver.
 **/

#include "WavePropagation2d_cuda.h"
#include <cuda_runtime.h>
#include <cstring>
#include <cstdio>

// Static variable to track initialization
static bool g_cuda_initialized = false;

// Forward declarations of kernel functions
namespace tsunami_lab {
  namespace patches {
    namespace cuda {
      __global__
      void xSweepKernel( const t_real *i_h_old,
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
      void ySweepKernel( const t_real *i_h_old,
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
      void initNewCellsKernel( const t_real *i_h_old,
                        const t_real *i_hu_old,
                        const t_real *i_hv_old,
                        t_real *o_h_new,
                        t_real *o_hu_new,
                        t_real *o_hv_new,
                        t_idx i_nValues );

      __global__
      void setGhostOutflowLeftRightKernel( t_real *io_h,
                                           t_real *io_hu,
                                           t_real *io_hv,
                                           t_idx i_nCellsX,
                                           t_idx i_nCellsY,
                                           t_idx i_stride );

      __global__
      void setGhostOutflowBottomTopKernel( t_real *io_h,
                                           t_real *io_hu,
                                           t_real *io_hv,
                                           t_idx i_nCellsX,
                                           t_idx i_nCellsY,
                                           t_idx i_stride );

      __global__
      void computeStepKernel( const t_real *i_h,
                              const t_real *i_hu,
                              const t_real *i_hv,
                              const t_real *i_b,
                              t_real *o_h,
                              t_real *o_hu,
                              t_real *o_hv,
                              t_idx i_nCellsX,
                              t_idx i_nCellsY,
                              t_idx i_stride,
                              t_real i_scaling,
                              bool i_useFWave );
    }
  }
}

bool tsunami_lab::patches::cuda::WavePropagation2dCuda::initialize( int i_deviceId ) {
  // Set device
  cudaError_t l_error = cudaSetDevice( i_deviceId );
  if( l_error != cudaSuccess ) {
    fprintf( stderr, "Error setting CUDA device: %s\n", cudaGetErrorString( l_error ) );
    return false;
  }

  g_cuda_initialized = true;
  return true;
}

void tsunami_lab::patches::cuda::WavePropagation2dCuda::finalize() {
  if( g_cuda_initialized ) {
    cudaDeviceReset();
    g_cuda_initialized = false;
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

  initNewCellsKernel<<<l_blocksPerGrid, l_threadsPerBlock>>>(
      l_h_old, l_hu_old, l_hv_old,
      l_h_new, l_hu_new, l_hv_new,
      m_nValues );

  cudaError_t l_launchErr = cudaGetLastError();
  if( l_launchErr != cudaSuccess ) {
    fprintf( stderr, "initNewCellsKernel launch error: %s\n", cudaGetErrorString( l_launchErr ) );
  }
  cudaError_t l_syncErr = cudaDeviceSynchronize();
  if( l_syncErr != cudaSuccess ) {
    fprintf( stderr, "initNewCellsKernel sync error: %s\n", cudaGetErrorString( l_syncErr ) );
  }
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
      (m_nCellsY + 1 + l_blockDim.y - 1) / l_blockDim.y
  );

  xSweepKernel<<<l_gridDim, l_blockDim>>>(
      l_h_old, l_hu_old, l_hv_old, m_d_b,
      l_h_new, l_hu_new,
      m_nCellsX, m_nCellsY, m_stride,
      i_scaling, m_useFWaveSolver );

  cudaError_t l_launchErr = cudaGetLastError();
  if( l_launchErr != cudaSuccess ) {
    fprintf( stderr, "xSweepKernel launch error: %s\n", cudaGetErrorString( l_launchErr ) );
  }
  cudaError_t l_syncErr = cudaDeviceSynchronize();
  if( l_syncErr != cudaSuccess ) {
    fprintf( stderr, "xSweepKernel sync error: %s\n", cudaGetErrorString( l_syncErr ) );
  }
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
      (m_nCellsX + 1 + l_blockDim.x - 1) / l_blockDim.x,
      (m_nCellsY + 1 + l_blockDim.y - 1) / l_blockDim.y
  );

  ySweepKernel<<<l_gridDim, l_blockDim>>>(
      l_h_old, l_hu_old, l_hv_old, m_d_b,
      l_h_new, l_hv_new,
      m_nCellsX, m_nCellsY, m_stride,
      i_scaling, m_useFWaveSolver );

  cudaError_t l_launchErr = cudaGetLastError();
  if( l_launchErr != cudaSuccess ) {
    fprintf( stderr, "ySweepKernel launch error: %s\n", cudaGetErrorString( l_launchErr ) );
  }
  cudaError_t l_syncErr = cudaDeviceSynchronize();
  if( l_syncErr != cudaSuccess ) {
    fprintf( stderr, "ySweepKernel sync error: %s\n", cudaGetErrorString( l_syncErr ) );
  }
}

void tsunami_lab::patches::cuda::WavePropagation2dCuda::swapBuffers() {
  m_step = 1 - m_step;
}

void tsunami_lab::patches::cuda::WavePropagation2dCuda::setGhostOutflow() {
  // Operate on the active buffer (the one the step reads this iteration).
  t_real *l_h  = m_d_h[m_step];
  t_real *l_hu = m_d_hu[m_step];
  t_real *l_hv = m_d_hv[m_step];

  int l_threadsPerBlock = 256;

  // Left/right ghost columns first (one thread per row), then bottom/top ghost
  // rows (one thread per column). The two launches are serialized by the stream,
  // so the corner ghost cells take their value from the bottom/top edge - this
  // matches the sequential two-loop CPU reference and avoids the data race that
  // a single fused kernel would have (overlapping ghost read/writes).
  int l_blocksLeftRight = ( (m_nCellsY + 2) + l_threadsPerBlock - 1 ) / l_threadsPerBlock;
  setGhostOutflowLeftRightKernel<<<l_blocksLeftRight, l_threadsPerBlock>>>(
      l_h, l_hu, l_hv,
      m_nCellsX, m_nCellsY, m_stride );

  cudaError_t l_lrErr = cudaGetLastError();
  if( l_lrErr != cudaSuccess ) {
    fprintf( stderr, "setGhostOutflowLeftRightKernel launch error: %s\n", cudaGetErrorString( l_lrErr ) );
  }

  int l_blocksBottomTop = ( (m_nCellsX + 2) + l_threadsPerBlock - 1 ) / l_threadsPerBlock;
  setGhostOutflowBottomTopKernel<<<l_blocksBottomTop, l_threadsPerBlock>>>(
      l_h, l_hu, l_hv,
      m_nCellsX, m_nCellsY, m_stride );

  cudaError_t l_launchErr = cudaGetLastError();
  if( l_launchErr != cudaSuccess ) {
    fprintf( stderr, "setGhostOutflowBottomTopKernel launch error: %s\n", cudaGetErrorString( l_launchErr ) );
  }
  cudaError_t l_syncErr = cudaDeviceSynchronize();
  if( l_syncErr != cudaSuccess ) {
    fprintf( stderr, "setGhostOutflow sync error: %s\n", cudaGetErrorString( l_syncErr ) );
  }
}

void tsunami_lab::patches::cuda::WavePropagation2dCuda::computeStep( t_real i_scaling ) {
  // One fused kernel reads only the old buffers and writes each new cell once,
  // so there are no atomics and no inter-kernel write hazards (race-free).
  t_real *l_h_old  = m_d_h[m_step];
  t_real *l_hu_old = m_d_hu[m_step];
  t_real *l_hv_old = m_d_hv[m_step];

  t_real *l_h_new  = m_d_h[1 - m_step];
  t_real *l_hu_new = m_d_hu[1 - m_step];
  t_real *l_hv_new = m_d_hv[1 - m_step];

  // One thread per cell over the full ghost-celled grid.
  dim3 l_blockDim( 16, 16 );
  dim3 l_gridDim(
      (m_nCellsX + 2 + l_blockDim.x - 1) / l_blockDim.x,
      (m_nCellsY + 2 + l_blockDim.y - 1) / l_blockDim.y
  );

  computeStepKernel<<<l_gridDim, l_blockDim>>>(
      l_h_old, l_hu_old, l_hv_old, m_d_b,
      l_h_new, l_hu_new, l_hv_new,
      m_nCellsX, m_nCellsY, m_stride,
      i_scaling, m_useFWaveSolver );

  cudaError_t l_launchErr = cudaGetLastError();
  if( l_launchErr != cudaSuccess ) {
    fprintf( stderr, "computeStepKernel launch error: %s\n", cudaGetErrorString( l_launchErr ) );
  }
  cudaError_t l_syncErr = cudaDeviceSynchronize();
  if( l_syncErr != cudaSuccess ) {
    fprintf( stderr, "computeStepKernel sync error: %s\n", cudaGetErrorString( l_syncErr ) );
  }
}
