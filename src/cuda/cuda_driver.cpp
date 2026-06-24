/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * CUDA driver implementation for tsunami simulation.
 **/

#include "cuda_driver.h"
#include <cuda_runtime.h>
#include <cstdio>
#include <cstring>

// Static variable to track initialization
static bool g_initialized = false;

bool tsunami_lab::cuda::CudaDriver::initialize( int i_deviceId ) {
  // Set device
  cudaError_t l_error = cudaSetDevice( i_deviceId );
  if( l_error != cudaSuccess ) {
    fprintf( stderr, "Error setting CUDA device: %s\n", cudaGetErrorString( l_error ) );
    return false;
  }

  g_initialized = true;
  return true;
}

bool tsunami_lab::cuda::CudaDriver::getDeviceInfo( char *o_deviceName,
                                                   size_t i_nameSize,
                                                   int o_computeCapability[2],
                                                   int &o_maxThreadsPerBlock ) {
  if( !g_initialized ) {
    fprintf( stderr, "CUDA not initialized\n" );
    return false;
  }

  cudaDeviceProp l_props;
  cudaError_t l_error = cudaGetDeviceProperties( &l_props, 0 );
  if( l_error != cudaSuccess ) {
    fprintf( stderr, "Error getting device properties: %s\n", cudaGetErrorString( l_error ) );
    return false;
  }

  strncpy( o_deviceName, l_props.name, i_nameSize - 1 );
  o_deviceName[i_nameSize - 1] = '\0';

  o_computeCapability[0] = l_props.major;
  o_computeCapability[1] = l_props.minor;

  o_maxThreadsPerBlock = l_props.maxThreadsPerBlock;

  return true;
}

void tsunami_lab::cuda::CudaDriver::finalize() {
  if( g_initialized ) {
    cudaDeviceReset();
    g_initialized = false;
  }
}

void tsunami_lab::cuda::CudaDriver::launchXSweepKernel(
    const t_real *i_h,
    const t_real *i_hu,
    const t_real *i_hv,
    const t_real *i_b,
    t_real *o_h,
    t_real *o_hu,
    t_idx i_nCellsX,
    t_idx i_nCellsY,
    t_idx i_stride,
    t_real i_scaling ) {
  // This is a wrapper function
  // In a full implementation, this would launch the kernel
  // For now, this is documented but not fully implemented
  // The actual kernel launch is in WavePropagation2d_cuda.cpp
}

void tsunami_lab::cuda::CudaDriver::launchYSweepKernel(
    const t_real *i_h,
    const t_real *i_hu,
    const t_real *i_hv,
    const t_real *i_b,
    t_real *o_h,
    t_real *o_hv,
    t_idx i_nCellsX,
    t_idx i_nCellsY,
    t_idx i_stride,
    t_real i_scaling ) {
  // This is a wrapper function
  // In a full implementation, this would launch the kernel
  // For now, this is documented but not fully implemented
  // The actual kernel launch is in WavePropagation2d_cuda.cpp
}

bool tsunami_lab::cuda::CudaDriver::allocateMemory(
    t_real **o_h,
    t_real **o_hu,
    t_real **o_hv,
    t_real **o_b,
    t_idx i_nValues ) {

  size_t l_size = i_nValues * sizeof(t_real);

  // Allocate two buffers for h, hu, hv (current and next step)
  for( int i = 0; i < 2; i++ ) {
    cudaError_t l_error = cudaMalloc( (void**)&o_h[i], l_size );
    if( l_error != cudaSuccess ) {
      fprintf( stderr, "Error allocating h buffer %d: %s\n", i, cudaGetErrorString( l_error ) );
      return false;
    }
    cudaMemset( o_h[i], 0, l_size );

    l_error = cudaMalloc( (void**)&o_hu[i], l_size );
    if( l_error != cudaSuccess ) {
      fprintf( stderr, "Error allocating hu buffer %d: %s\n", i, cudaGetErrorString( l_error ) );
      return false;
    }
    cudaMemset( o_hu[i], 0, l_size );

    l_error = cudaMalloc( (void**)&o_hv[i], l_size );
    if( l_error != cudaSuccess ) {
      fprintf( stderr, "Error allocating hv buffer %d: %s\n", i, cudaGetErrorString( l_error ) );
      return false;
    }
    cudaMemset( o_hv[i], 0, l_size );
  }

  // Allocate bathymetry buffer (single copy)
  cudaError_t l_error = cudaMalloc( (void**)o_b, l_size );
  if( l_error != cudaSuccess ) {
    fprintf( stderr, "Error allocating bathymetry buffer: %s\n", cudaGetErrorString( l_error ) );
    return false;
  }
  cudaMemset( *o_b, 0, l_size );

  return true;
}

void tsunami_lab::cuda::CudaDriver::freeMemory(
    t_real **i_h,
    t_real **i_hu,
    t_real **i_hv,
    t_real *i_b ) {

  // Free h, hu, hv buffers
  for( int i = 0; i < 2; i++ ) {
    if( i_h[i] != nullptr ) {
      cudaFree( i_h[i] );
      i_h[i] = nullptr;
    }
    if( i_hu[i] != nullptr ) {
      cudaFree( i_hu[i] );
      i_hu[i] = nullptr;
    }
    if( i_hv[i] != nullptr ) {
      cudaFree( i_hv[i] );
      i_hv[i] = nullptr;
    }
  }

  // Free bathymetry buffer
  if( i_b != nullptr ) {
    cudaFree( i_b );
  }
}

void tsunami_lab::cuda::CudaDriver::copyHostToDevice(
    const t_real *i_host,
    t_real *o_device,
    t_idx i_nValues ) {

  size_t l_size = i_nValues * sizeof(t_real);
  cudaError_t l_error = cudaMemcpy( o_device, i_host, l_size, cudaMemcpyHostToDevice );

  if( l_error != cudaSuccess ) {
    fprintf( stderr, "Error copying data to device: %s\n", cudaGetErrorString( l_error ) );
  }
}

void tsunami_lab::cuda::CudaDriver::copyDeviceToHost(
    const t_real *i_device,
    t_real *o_host,
    t_idx i_nValues ) {

  size_t l_size = i_nValues * sizeof(t_real);
  cudaError_t l_error = cudaMemcpy( o_host, i_device, l_size, cudaMemcpyDeviceToHost );

  if( l_error != cudaSuccess ) {
    fprintf( stderr, "Error copying data to host: %s\n", cudaGetErrorString( l_error ) );
  }
}

void tsunami_lab::cuda::CudaDriver::synchronize() {
  cudaDeviceSynchronize();
}

bool tsunami_lab::cuda::CudaDriver::checkError( const char *i_message ) {
  cudaError_t l_error = cudaGetLastError();
  if( l_error != cudaSuccess ) {
    fprintf( stderr, "%s: %s\n", i_message, cudaGetErrorString( l_error ) );
    return true;
  }
  return false;
}
