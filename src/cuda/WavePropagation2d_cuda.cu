/**
 *
 * @section DESCRIPTION
 * GPU wrapper implementation for two-dimensional wave propagation solver.
 **/

#include "WavePropagation2d_cuda.h"
#include "../constants.h"
#include <cuda_runtime.h>
#include <cstring>
#include <cstdio>
#include <vector>

// Static variable to track initialization
static bool g_cuda_initialized = false;

// Forward declarations of kernel functions
namespace tsunami_lab {
  namespace patches {
    namespace cuda {
      __global__
      void setGhostOutflowLeftRightKernelRowMajor( t_real *io_h,
                                                   t_real *io_hu,
                                                   t_real *io_hv,
                                                   t_real *io_b,
                                                   t_idx i_nCellsX,
                                                   t_idx i_nCellsY,
                                                   t_idx i_stride );

      __global__
      void setGhostOutflowBottomTopKernelRowMajor( t_real *io_h,
                                                    t_real *io_hu,
                                                    t_real *io_hv,
                                                    t_real *io_b,
                                                    t_idx i_nCellsX,
                                                    t_idx i_nCellsY,
                                                    t_idx i_stride );

      __global__
      void computeXSweepKernelRowMajor( const t_real *i_h,
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
      void computeYSweepKernelRowMajor( const t_real *i_h,
                                         const t_real *i_hv,
                                         const t_real *i_b,
                                         t_real *io_h,
                                         t_real *o_hv,
                                         t_idx i_nCellsX,
                                         t_idx i_nCellsY,
                                         t_idx i_stride,
                                         t_real i_scaling,
                                         bool i_useFWave );

      __global__
      void initNewCellsKernelRowMajor( const t_real *i_h,
                                        const t_real *i_hu,
                                        const t_real *i_hv,
                                        t_real *o_h,
                                        t_real *o_hu,
                                        t_real *o_hv,
                                        t_idx i_nCellsX,
                                        t_idx i_nCellsY,
                                        t_idx i_stride );

      __global__
      void computeXSweepAtomicKernelRowMajor( const t_real *i_h,
                                               const t_real *i_hu,
                                               const t_real *i_hv,
                                               const t_real *i_b,
                                               t_real *io_h,
                                               t_real *io_hu,
                                               t_idx i_nCellsX,
                                               t_idx i_nCellsY,
                                               t_idx i_stride,
                                               t_real i_scaling,
                                               bool i_useFWave );

      __global__
      void computeYSweepAtomicKernelRowMajor( const t_real *i_h,
                                               const t_real *i_hv,
                                               const t_real *i_b,
                                               t_real *io_h,
                                               t_real *io_hv,
                                               t_idx i_nCellsX,
                                               t_idx i_nCellsY,
                                               t_idx i_stride,
                                               t_real i_scaling,
                                               bool i_useFWave );

      __global__
      void applyWetDryThresholdKernelRowMajor( t_real *io_h,
                                                t_real *io_hu,
                                                t_real *io_hv,
                                                t_idx i_nCellsX,
                                                t_idx i_nCellsY,
                                                t_idx i_stride );

      __global__
      void setGhostOutflowLeftRightKernelColumnMajor( t_real *io_h,
                                                      t_real *io_hu,
                                                      t_real *io_hv,
                                                      t_idx i_nCellsX,
                                                      t_idx i_nCellsY,
                                                      t_idx i_height );

      __global__
      void setGhostOutflowBottomTopKernelColumnMajor( t_real *io_h,
                                                       t_real *io_hu,
                                                       t_real *io_hv,
                                                       t_real *io_b,
                                                       t_idx i_nCellsX,
                                                       t_idx i_nCellsY,
                                                       t_idx i_height );

      __global__
      void computeXSweepKernelColumnMajor( const t_real *i_h,
                                            const t_real *i_hu,
                                            const t_real *i_hv,
                                            const t_real *i_b,
                                            t_real *o_h,
                                            t_real *o_hu,
                                            t_real *o_hv,
                                            t_idx i_nCellsX,
                                            t_idx i_nCellsY,
                                            t_idx i_height,
                                            t_real i_scaling,
                                            bool i_useFWave );

      __global__
      void computeYSweepKernelColumnMajor( const t_real *i_h,
                                            const t_real *i_hv,
                                            const t_real *i_b,
                                            t_real *io_h,
                                            t_real *o_hv,
                                            t_idx i_nCellsX,
                                            t_idx i_nCellsY,
                                            t_idx i_height,
                                            t_real i_scaling,
                                            bool i_useFWave );

      __global__
      void initNewCellsKernelColumnMajor( const t_real *i_h,
                                           const t_real *i_hu,
                                           const t_real *i_hv,
                                           t_real *o_h,
                                           t_real *o_hu,
                                           t_real *o_hv,
                                           t_idx i_nCellsX,
                                           t_idx i_nCellsY,
                                           t_idx i_height );

      __global__
      void computeXSweepAtomicKernelColumnMajor( const t_real *i_h,
                                                 const t_real *i_hu,
                                                 const t_real *i_hv,
                                                 const t_real *i_b,
                                                 t_real *io_h,
                                                 t_real *io_hu,
                                                 t_idx i_nCellsX,
                                                 t_idx i_nCellsY,
                                                 t_idx i_height,
                                                 t_real i_scaling,
                                                 bool i_useFWave );

      __global__
      void computeYSweepAtomicKernelColumnMajor( const t_real *i_h,
                                                 const t_real *i_hv,
                                                 const t_real *i_b,
                                                 t_real *io_h,
                                                 t_real *io_hv,
                                                 t_idx i_nCellsX,
                                                 t_idx i_nCellsY,
                                                 t_idx i_height,
                                                 t_real i_scaling,
                                                 bool i_useFWave );

      __global__
      void applyWetDryThresholdKernelColumnMajor( t_real *io_h,
                                                   t_real *io_hu,
                                                   t_real *io_hv,
                                                   t_idx i_nCellsX,
                                                   t_idx i_nCellsY,
                                                   t_idx i_height );

__global__
void setGhostOutflowLeftRightKernelTile(
    t_real *io_h,
    t_real *io_hu,
    t_real *io_hv,
    t_real *io_b,
    t_idx i_nCellsX,
    t_idx i_nCellsY,
    t_idx i_stride,
    t_idx i_height,
    t_idx i_nTilesX,
    t_idx i_nTilesY );

__global__
void setGhostOutflowBottomTopKernelTile(
    t_real *io_h,
    t_real *io_hu,
    t_real *io_hv,
    t_real *io_b,
    t_idx i_nCellsX,
    t_idx i_nCellsY,
    t_idx i_tileSizeX, //used to be stride
    t_idx i_tileSizeY, //used to be height
    t_idx i_nTilesX,
    t_idx i_nTilesY );

__global__
void computeXSweepKernelTile(
    const t_real *i_h,
    const t_real *i_hu,
    const t_real *i_hv,
    const t_real *i_b,
    t_real *o_h,
    t_real *o_hu,
    t_real *o_hv,
    t_idx i_tileSizeX,
    t_idx i_tileSizeY,
    t_idx i_stride,
    t_idx i_height,
    t_idx i_nCellsX,
    t_idx i_nCellsY,
    t_real i_scaling,
    bool i_useFWave );

__global__
void computeYSweepKernelTile(
    const t_real *i_h,
    const t_real *i_hv,
    const t_real *i_b,
    t_real *io_h,
    t_real *o_hv,
    t_idx i_nCellsX,
    t_idx i_nCellsY,
    t_idx i_stride,
    t_idx i_height,
    t_idx i_nTilesX,
    t_idx i_nTilesY,
    t_real i_scaling,
    bool i_useFWave );

__global__
void initNewCellsKernelTile(
    const t_real *i_h,
    const t_real *i_hu,
    const t_real *i_hv,
    t_real *o_h,
    t_real *o_hu,
    t_real *o_hv,
    t_idx i_nCellsX,
    t_idx i_nCellsY,
    t_idx i_stride,
    t_idx i_height,
    t_idx i_nTilesX,
    t_idx i_nTilesY );

__global__
void computeXSweepAtomicKernelTile(
    const t_real *i_h,
    const t_real *i_hu,
    const t_real *i_hv,
    const t_real *i_b,
    t_real *io_h,
    t_real *io_hu,
    t_idx i_nCellsX,
    t_idx i_nCellsY,
    t_idx i_stride,
    t_idx i_height,
    t_idx i_nTilesX,
    t_idx i_nTilesY,
    t_real i_scaling,
    bool i_useFWave );

__global__
void computeYSweepAtomicKernelTile(
    const t_real *i_h,
    const t_real *i_hv,
    const t_real *i_b,
    t_real *io_h,
    t_real *io_hv,
    t_idx i_nCellsX,
    t_idx i_nCellsY,
    t_idx i_stride,
    t_idx i_height,
    t_idx i_nTilesX,
    t_idx i_nTilesY,
    t_real i_scaling,
    bool i_useFWave );

__global__
void applyWetDryThresholdKernelTile(
    t_real *io_h,
    t_real *io_hu,
    t_real *io_hv,
    t_idx i_nCellsX,
    t_idx i_nCellsY,
    t_idx i_stride,
    t_idx i_height,
    t_idx i_nTilesX,
    t_idx i_nTilesY
 );

namespace {
  inline tsunami_lab::t_idx tileIndex( tsunami_lab::t_idx i_row,
                                       tsunami_lab::t_idx i_col,
                                       tsunami_lab::t_idx i_nTilesX,
                                       tsunami_lab::t_idx i_tileSize = 32 ) {
    const tsunami_lab::t_idx l_tileX = i_col / i_tileSize;
    const tsunami_lab::t_idx l_tileY = i_row / i_tileSize;
    const tsunami_lab::t_idx l_localX = i_col % i_tileSize;
    const tsunami_lab::t_idx l_localY = i_row % i_tileSize;
    return ((l_tileY * i_nTilesX + l_tileX) * i_tileSize * i_tileSize) +
           (l_localY * i_tileSize + l_localX);
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
    int  i_blockWidth,
    MemoryLayout i_memoryLayout ) : m_nCellsX( i_nCellsX ),
                                    m_nCellsY( i_nCellsY ),
                                    m_useFWaveSolver( i_useFWaveSolver ),
                                    m_blockWidth( i_blockWidth ),
                                    m_memoryLayout( i_memoryLayout ) {

  // Stride includes ghost cells on both sides: nCellsX + 2
  m_stride = i_nCellsX + 2;
  m_height = i_nCellsY + 2;

  // Total number of values depends on the selected memory layout.
  if( m_memoryLayout == MemoryLayout::Tile32 ) {
    m_nTilesX = (m_stride + 31) / 32;
    m_nTilesY = (m_height + 31) / 32;
    m_nValues = m_nTilesX * m_nTilesY * 32 * 32;
  }
  else {
    m_nValues = m_height * m_stride;
  }

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

  if( m_memoryLayout == MemoryLayout::Tile32 ) {
    std::vector<t_real> l_hTiled( m_nValues, 0.0 );
    std::vector<t_real> l_huTiled( m_nValues, 0.0 );
    std::vector<t_real> l_hvTiled( m_nValues, 0.0 );
    std::vector<t_real> l_bTiled( m_nValues, 0.0 );

    for( t_idx l_row = 0; l_row < m_height; ++l_row ) {
      for( t_idx l_col = 0; l_col < m_stride; ++l_col ) {
        const t_idx l_src = l_row * m_stride + l_col;
        const t_idx l_dst = tileIndex( l_row, l_col, m_nTilesX );
        l_hTiled[l_dst] = i_h[l_src];
        l_huTiled[l_dst] = i_hu[l_src];
        l_hvTiled[l_dst] = i_hv[l_src];
        l_bTiled[l_dst] = i_b[l_src];
      }
    }

    cudaMemcpy( m_d_h[m_step], l_hTiled.data(), l_size, cudaMemcpyHostToDevice );
    cudaMemcpy( m_d_hu[m_step], l_huTiled.data(), l_size, cudaMemcpyHostToDevice );
    cudaMemcpy( m_d_hv[m_step], l_hvTiled.data(), l_size, cudaMemcpyHostToDevice );
    cudaMemcpy( m_d_b, l_bTiled.data(), l_size, cudaMemcpyHostToDevice );
  }
  else {
    cudaMemcpy( m_d_h[m_step], i_h, l_size, cudaMemcpyHostToDevice );
    cudaMemcpy( m_d_hu[m_step], i_hu, l_size, cudaMemcpyHostToDevice );
    cudaMemcpy( m_d_hv[m_step], i_hv, l_size, cudaMemcpyHostToDevice );
    cudaMemcpy( m_d_b, i_b, l_size, cudaMemcpyHostToDevice );
  }
}

void tsunami_lab::patches::cuda::WavePropagation2dCuda::copyToHost(
    t_real *o_h,
    t_real *o_hu,
    t_real *o_hv ) {

  size_t l_size = m_nValues * sizeof(t_real);

  if( m_memoryLayout == MemoryLayout::Tile32 ) {
    std::vector<t_real> l_hTiled( m_nValues, 0.0 );
    std::vector<t_real> l_huTiled( m_nValues, 0.0 );
    std::vector<t_real> l_hvTiled( m_nValues, 0.0 );

    cudaMemcpy( l_hTiled.data(), m_d_h[m_step], l_size, cudaMemcpyDeviceToHost );
    cudaMemcpy( l_huTiled.data(), m_d_hu[m_step], l_size, cudaMemcpyDeviceToHost );
    cudaMemcpy( l_hvTiled.data(), m_d_hv[m_step], l_size, cudaMemcpyDeviceToHost );

    for( t_idx l_row = 0; l_row < m_height; ++l_row ) {
      for( t_idx l_col = 0; l_col < m_stride; ++l_col ) {
        const t_idx l_src = tileIndex( l_row, l_col, m_nTilesX );
        const t_idx l_dst = l_row * m_stride + l_col;
        o_h[l_dst] = l_hTiled[l_src];
        o_hu[l_dst] = l_huTiled[l_src];
        o_hv[l_dst] = l_hvTiled[l_src];
      }
    }
  }
  else {
    cudaMemcpy( o_h, m_d_h[m_step], l_size, cudaMemcpyDeviceToHost );
    cudaMemcpy( o_hu, m_d_hu[m_step], l_size, cudaMemcpyDeviceToHost );
    cudaMemcpy( o_hv, m_d_hv[m_step], l_size, cudaMemcpyDeviceToHost );
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
  t_real *l_b  = m_d_b;

  int l_threadsPerBlock = m_blockWidth * m_blockWidth;

  // Left/right ghost columns first (one thread per row), then bottom/top ghost
  // rows (one thread per column). The two launches are serialized by the stream,
  // so the corner ghost cells take their value from the bottom/top edge - this
  // matches the sequential two-loop CPU reference and avoids the data race that
  // a single fused kernel would have (overlapping ghost read/writes).
  int l_blocksLeftRight = ( (m_nCellsY + 2) + l_threadsPerBlock - 1 ) / l_threadsPerBlock;
  if( m_memoryLayout == MemoryLayout::ColumnMajor ) {
    setGhostOutflowLeftRightKernelColumnMajor<<<l_blocksLeftRight, l_threadsPerBlock>>>(
        l_h, l_hu, l_hv, 
        m_nCellsX, m_nCellsY, m_height );
  }
  else if( m_memoryLayout == MemoryLayout::Tile32 ) {
    setGhostOutflowLeftRightKernelTile<<<l_blocksLeftRight, l_threadsPerBlock>>>(
        l_h, l_hu, l_hv, l_b,
        m_nCellsX, m_nCellsY, m_stride, m_height, m_nTilesX, m_nTilesY );
  }
  else {
    setGhostOutflowLeftRightKernelRowMajor<<<l_blocksLeftRight, l_threadsPerBlock>>>(
        l_h, l_hu, l_hv, l_b,
        m_nCellsX, m_nCellsY, m_stride );
  }

  cudaError_t l_lrErr = cudaGetLastError();
  if( l_lrErr != cudaSuccess ) {
    fprintf( stderr, "setGhostOutflowLeftRightKernel launch error: %s\n", cudaGetErrorString( l_lrErr ) );
  }

  int l_blocksBottomTop = ( (m_nCellsX + 2) + l_threadsPerBlock - 1 ) / l_threadsPerBlock;
  if( m_memoryLayout == MemoryLayout::ColumnMajor ) {
setGhostOutflowBottomTopKernelColumnMajor<<<l_blocksBottomTop, l_threadsPerBlock>>>(
    l_h, l_hu, l_hv, l_b,
    m_nCellsX, m_nCellsY, m_height );
  }
  else if( m_memoryLayout == MemoryLayout::Tile32 ) {
    setGhostOutflowBottomTopKernelTile<<<l_blocksBottomTop, l_threadsPerBlock>>>(
        l_h, l_hu, l_hv, l_b,
        m_nCellsX, m_nCellsY, m_stride, m_height, m_nTilesX, m_nTilesY );
  }
  else {
    setGhostOutflowBottomTopKernelRowMajor<<<l_blocksBottomTop, l_threadsPerBlock>>>(
        l_h, l_hu, l_hv, l_b,
        m_nCellsX, m_nCellsY, m_stride );
  }

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
  setGhostOutflow();

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
  if( m_memoryLayout == MemoryLayout::ColumnMajor ) {
    computeXSweepKernelColumnMajor<<<l_gridDim, l_blockDim>>>(
        l_h_old, l_hu_old, l_hv_old, m_d_b,
        l_h_new, l_hu_new, l_hv_new,
        m_nCellsX, m_nCellsY, m_height,
        i_scaling, m_useFWaveSolver );
  }
  else if( m_memoryLayout == MemoryLayout::Tile32 ) {
    computeXSweepKernelTile<<<l_gridDim, l_blockDim>>>(
        l_h_old, l_hu_old, l_hv_old, m_d_b,
        l_h_new, l_hu_new, l_hv_new,
        m_nCellsX, m_nCellsY, m_stride, m_height, m_nTilesX, m_nTilesY,
        i_scaling, m_useFWaveSolver );
  }
  else {
    computeXSweepKernelRowMajor<<<l_gridDim, l_blockDim>>>(
        l_h_old, l_hu_old, l_hv_old, m_d_b,
        l_h_new, l_hu_new, l_hv_new,
        m_nCellsX, m_nCellsY, m_stride,
        i_scaling, m_useFWaveSolver );
  }

  cudaError_t l_xErr = cudaGetLastError();
  if( l_xErr != cudaSuccess ) {
    fprintf( stderr, "computeXSweepKernel launch error: %s\n", cudaGetErrorString( l_xErr ) );
  }

  // Y-sweep: reads OLD h/hv for solver inputs; accumulates y-updates onto
  // the x-swept h_new and writes the final hv_new.
  // The default stream serialises launches, so x-sweep output is visible here.
  if( m_memoryLayout == MemoryLayout::ColumnMajor ) {
    computeYSweepKernelColumnMajor<<<l_gridDim, l_blockDim>>>(
        l_h_old, l_hv_old, m_d_b,
        l_h_new, l_hv_new,
        m_nCellsX, m_nCellsY, m_height,
        i_scaling, m_useFWaveSolver );
  }
  else if( m_memoryLayout == MemoryLayout::Tile32 ) {
    computeYSweepKernelTile<<<l_gridDim, l_blockDim>>>(
        l_h_old, l_hv_old, m_d_b,
        l_h_new, l_hv_new,
        m_nCellsX, m_nCellsY, m_stride, m_height, m_nTilesX, m_nTilesY,
        i_scaling, m_useFWaveSolver );
  }
  else {
    computeYSweepKernelRowMajor<<<l_gridDim, l_blockDim>>>(
        l_h_old, l_hv_old, m_d_b,
        l_h_new, l_hv_new,
        m_nCellsX, m_nCellsY, m_stride,
        i_scaling, m_useFWaveSolver );
  }

  cudaError_t l_yErr = cudaGetLastError();
  if( l_yErr != cudaSuccess ) {
    fprintf( stderr, "computeYSweepKernel launch error: %s\n", cudaGetErrorString( l_yErr ) );
  }

  // Apply wet/dry threshold to eliminate NaNs in dry cells
  if( m_memoryLayout == MemoryLayout::ColumnMajor ) {
    applyWetDryThresholdKernelColumnMajor<<<l_gridDim, l_blockDim>>>(
        l_h_new, l_hu_new, l_hv_new,
        m_nCellsX, m_nCellsY, m_height );
  }
  else if( m_memoryLayout == MemoryLayout::Tile32 ) {
    applyWetDryThresholdKernelTile<<<l_gridDim, l_blockDim>>>(
        l_h_new, l_hu_new, l_hv_new,
        m_nCellsX, m_nCellsY, m_stride, m_height, m_nTilesX, m_nTilesY );
  }
  else {
    applyWetDryThresholdKernelRowMajor<<<l_gridDim, l_blockDim>>>(
        l_h_new, l_hu_new, l_hv_new,
        m_nCellsX, m_nCellsY, m_stride );
  }

  cudaError_t l_threshErr = cudaGetLastError();
  if( l_threshErr != cudaSuccess ) {
    fprintf( stderr, "applyWetDryThresholdKernel launch error: %s\n", cudaGetErrorString( l_threshErr ) );
  }

  cudaError_t l_syncErr = cudaDeviceSynchronize();
  if( l_syncErr != cudaSuccess ) {
    fprintf( stderr, "computeStep sync error: %s\n", cudaGetErrorString( l_syncErr ) );
  }
}

void tsunami_lab::patches::cuda::WavePropagation2dCuda::computeStepAtomic( t_real i_scaling ) {
  setGhostOutflow();

  // Atomic-based implementation: one thread per edge, uses atomicAdd to accumulate
  // updates into the new buffers. Follows the CPU reference structure more closely.
  // Requires manual initialization of new buffers before each sweep.
  
  t_real *l_h_old  = m_d_h[m_step];
  t_real *l_hu_old = m_d_hu[m_step];
  t_real *l_hv_old = m_d_hv[m_step];

  t_real *l_h_new  = m_d_h[1 - m_step];
  t_real *l_hu_new = m_d_hu[1 - m_step];
  t_real *l_hv_new = m_d_hv[1 - m_step];

  dim3 l_blockDim( m_blockWidth, m_blockWidth );
  
  // Grid dimensions for 2D kernels (cell initialization)
  dim3 l_gridDimCells(
      (m_nCellsX + 2 + l_blockDim.x - 1) / l_blockDim.x,
      (m_nCellsY + 2 + l_blockDim.y - 1) / l_blockDim.y
  );

  // Initialize new buffers with values from old buffers
  if( m_memoryLayout == MemoryLayout::ColumnMajor ) {
    initNewCellsKernelColumnMajor<<<l_gridDimCells, l_blockDim>>>(
        l_h_old, l_hu_old, l_hv_old,
        l_h_new, l_hu_new, l_hv_new,
        m_nCellsX, m_nCellsY, m_height );
  }
  else if( m_memoryLayout == MemoryLayout::Tile32 ) {
    initNewCellsKernelTile<<<l_gridDimCells, l_blockDim>>>(
        l_h_old, l_hu_old, l_hv_old,
        l_h_new, l_hu_new, l_hv_new,
        m_nCellsX, m_nCellsY, m_stride, m_height, m_nTilesX, m_nTilesY );
  }
  else {
    initNewCellsKernelRowMajor<<<l_gridDimCells, l_blockDim>>>(
        l_h_old, l_hu_old, l_hv_old,
        l_h_new, l_hu_new, l_hv_new,
        m_nCellsX, m_nCellsY, m_stride );
  }

  cudaError_t l_initErr = cudaGetLastError();
  if( l_initErr != cudaSuccess ) {
    fprintf( stderr, "initNewCellsKernel launch error: %s\n", cudaGetErrorString( l_initErr ) );
  }

  // X-sweep: one thread per edge (nCellsX+1) * nCellsY edges
  // Grid: x dimension = 0..nCellsX (nCellsX+1 threads), y dimension = 1..nCellsY
  dim3 l_gridDimXSweep(
      (m_nCellsX + 1 + l_blockDim.x - 1) / l_blockDim.x,
      (m_nCellsY + l_blockDim.y - 1) / l_blockDim.y
  );

  if( m_memoryLayout == MemoryLayout::ColumnMajor ) {
    computeXSweepAtomicKernelColumnMajor<<<l_gridDimXSweep, l_blockDim>>>(
        l_h_old, l_hu_old, l_hv_old, m_d_b,
        l_h_new, l_hu_new,
        m_nCellsX, m_nCellsY, m_height,
        i_scaling, m_useFWaveSolver );
  }
  else if( m_memoryLayout == MemoryLayout::Tile32 ) {
    computeXSweepAtomicKernelTile<<<l_gridDimXSweep, l_blockDim>>>(
        l_h_old, l_hu_old, l_hv_old, m_d_b,
        l_h_new, l_hu_new,
        m_nCellsX, m_nCellsY, m_stride, m_height, m_nTilesX, m_nTilesY,
        i_scaling, m_useFWaveSolver );
  }
  else {
    computeXSweepAtomicKernelRowMajor<<<l_gridDimXSweep, l_blockDim>>>(
        l_h_old, l_hu_old, l_hv_old, m_d_b,
        l_h_new, l_hu_new,
        m_nCellsX, m_nCellsY, m_stride,
        i_scaling, m_useFWaveSolver );
  }

  cudaError_t l_xErr = cudaGetLastError();
  if( l_xErr != cudaSuccess ) {
    fprintf( stderr, "computeXSweepAtomicKernel launch error: %s\n", cudaGetErrorString( l_xErr ) );
  }

  // Y-sweep: one thread per edge nCellsX * (nCellsY+1) edges
  // Grid: x dimension = 1..nCellsX, y dimension = 0..nCellsY (nCellsY+1 threads)
  dim3 l_gridDimYSweep(
      (m_nCellsX + l_blockDim.x - 1) / l_blockDim.x,
      (m_nCellsY + 1 + l_blockDim.y - 1) / l_blockDim.y
  );

  if( m_memoryLayout == MemoryLayout::ColumnMajor ) {
    computeYSweepAtomicKernelColumnMajor<<<l_gridDimYSweep, l_blockDim>>>(
        l_h_old, l_hv_old, m_d_b,
        l_h_new, l_hv_new,
        m_nCellsX, m_nCellsY, m_height,
        i_scaling, m_useFWaveSolver );
  }
  else if( m_memoryLayout == MemoryLayout::Tile32 ) {
    computeYSweepAtomicKernelTile<<<l_gridDimYSweep, l_blockDim>>>(
        l_h_old, l_hv_old, m_d_b,
        l_h_new, l_hv_new,
        m_nCellsX, m_nCellsY, m_stride, m_height, m_nTilesX, m_nTilesY,
        i_scaling, m_useFWaveSolver );
  }
  else {
    computeYSweepAtomicKernelRowMajor<<<l_gridDimYSweep, l_blockDim>>>(
        l_h_old, l_hv_old, m_d_b,
        l_h_new, l_hv_new,
        m_nCellsX, m_nCellsY, m_stride,
        i_scaling, m_useFWaveSolver );
  }

  cudaError_t l_yErr = cudaGetLastError();
  if( l_yErr != cudaSuccess ) {
    fprintf( stderr, "computeYSweepAtomicKernel launch error: %s\n", cudaGetErrorString( l_yErr ) );
  }

  // Apply wet/dry threshold to eliminate NaNs in dry cells
  dim3 l_gridDimThresh(
      (m_nCellsX + 2 + l_blockDim.x - 1) / l_blockDim.x,
      (m_nCellsY + 2 + l_blockDim.y - 1) / l_blockDim.y
  );
  if( m_memoryLayout == MemoryLayout::ColumnMajor ) {
    applyWetDryThresholdKernelColumnMajor<<<l_gridDimThresh, l_blockDim>>>(
        l_h_new, l_hu_new, l_hv_new,
        m_nCellsX, m_nCellsY, m_height );
  }
  else if( m_memoryLayout == MemoryLayout::Tile32 ) {
    applyWetDryThresholdKernelTile<<<l_gridDimThresh, l_blockDim>>>(
        l_h_new, l_hu_new, l_hv_new,
        m_nCellsX, m_nCellsY, m_stride, m_height, m_nTilesX, m_nTilesY );
  }
  else {
    applyWetDryThresholdKernelRowMajor<<<l_gridDimThresh, l_blockDim>>>(
        l_h_new, l_hu_new, l_hv_new,
        m_nCellsX, m_nCellsY, m_stride );
  }

  cudaError_t l_threshErr = cudaGetLastError();
  if( l_threshErr != cudaSuccess ) {
    fprintf( stderr, "applyWetDryThresholdKernel launch error: %s\n", cudaGetErrorString( l_threshErr ) );
  }

  cudaError_t l_syncErr = cudaDeviceSynchronize();
  if( l_syncErr != cudaSuccess ) {
    fprintf( stderr, "computeStepAtomic sync error: %s\n", cudaGetErrorString( l_syncErr ) );
  }
}

    } // namespace cuda
  } // namespace patches
} // namespace tsunami_lab
