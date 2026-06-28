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
      void computeXSweepKernel( const t_real *i_h,
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

      __global__
      void computeYSweepKernel( const t_real *i_h,
                                 const t_real *i_hv,
                                 const t_real *i_b,
                                 t_real *io_h,
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
    bool i_useFWaveSolver,
    int  i_blockWidth ) : m_nCellsX( i_nCellsX ),
                          m_nCellsY( i_nCellsY ),
                          m_useFWaveSolver( i_useFWaveSolver ),
                          m_blockWidth( i_blockWidth ) {

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

  // Copy bathymetry (constant throughout the simulation).
  // The caller is responsible for having set ghost cells on the source array
  // before this call (e.g. via WavePropagation2d::setGhostOutflow()), so that
  // the bathymetry ghost cells copied here are already correct.
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

void tsunami_lab::patches::cuda::WavePropagation2dCuda::swapBuffers() {
  m_step = 1 - m_step;
}

void tsunami_lab::patches::cuda::WavePropagation2dCuda::setGhostOutflow() {
  // Operate on the active buffer (the one the step reads this iteration).
  t_real *l_h  = m_d_h[m_step];
  t_real *l_hu = m_d_hu[m_step];
  t_real *l_hv = m_d_hv[m_step];

  int l_threadsPerBlock = m_blockWidth * m_blockWidth;

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
  // Operator-split implementation: x-sweep then y-sweep, matching the CPU
  // reference exactly. Each kernel is one thread per cell and requires no
  // atomics. The two kernels are launched into the same (default) stream so
  // the x-sweep is guaranteed to finish before the y-sweep starts.
  t_real *l_h_old  = m_d_h[m_step];
  t_real *l_hu_old = m_d_hu[m_step];
  t_real *l_hv_old = m_d_hv[m_step];

  t_real *l_h_new  = m_d_h[1 - m_step];
  t_real *l_hu_new = m_d_hu[1 - m_step];
  t_real *l_hv_new = m_d_hv[1 - m_step];

  dim3 l_blockDim( m_blockWidth, m_blockWidth );
  dim3 l_gridDim(
      (m_nCellsX + 2 + l_blockDim.x - 1) / l_blockDim.x,
      (m_nCellsY + 2 + l_blockDim.y - 1) / l_blockDim.y
  );

  // X-sweep: apply x-direction net-updates to h and hu; copy hv through.
  computeXSweepKernel<<<l_gridDim, l_blockDim>>>(
      l_h_old, l_hu_old, l_hv_old, m_d_b,
      l_h_new, l_hu_new, l_hv_new,
      m_nCellsX, m_nCellsY, m_stride,
      i_scaling, m_useFWaveSolver );

  cudaError_t l_xErr = cudaGetLastError();
  if( l_xErr != cudaSuccess ) {
    fprintf( stderr, "computeXSweepKernel launch error: %s\n", cudaGetErrorString( l_xErr ) );
  }

  // Y-sweep: reads OLD h/hv for solver inputs; accumulates y-updates onto
  // the x-swept h_new and writes the final hv_new.
  // The default stream serialises launches, so x-sweep output is visible here.
  computeYSweepKernel<<<l_gridDim, l_blockDim>>>(
      l_h_old, l_hv_old, m_d_b,
      l_h_new, l_hv_new,
      m_nCellsX, m_nCellsY, m_stride,
      i_scaling, m_useFWaveSolver );

  cudaError_t l_yErr = cudaGetLastError();
  if( l_yErr != cudaSuccess ) {
    fprintf( stderr, "computeYSweepKernel launch error: %s\n", cudaGetErrorString( l_yErr ) );
  }
  cudaError_t l_syncErr = cudaDeviceSynchronize();
  if( l_syncErr != cudaSuccess ) {
    fprintf( stderr, "computeStep sync error: %s\n", cudaGetErrorString( l_syncErr ) );
  }
}
