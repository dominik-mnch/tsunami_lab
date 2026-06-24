/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * CUDA driver for tsunami simulation.
 * Provides high-level interface for GPU acceleration.
 **/

#ifndef TSUNAMI_LAB_CUDA_DRIVER_H
#define TSUNAMI_LAB_CUDA_DRIVER_H

#include "../constants.h"

namespace tsunami_lab {
  namespace cuda {
    class CudaDriver;
  }
}

/**
 * CUDA driver for GPU-accelerated wave propagation.
 * 
 * Manages the interface between CPU solver and GPU implementation.
 * Handles memory transfers, kernel launches, and synchronization.
 **/
class tsunami_lab::cuda::CudaDriver {

  public:
    /**
     * Initialize CUDA device and query capabilities.
     * 
     * @param i_deviceId GPU device ID to use (default: 0)
     * @return true if initialization successful, false otherwise
     **/
    static bool initialize( int i_deviceId = 0 );

    /**
     * Get information about selected CUDA device.
     * 
     * @param o_deviceName buffer to store device name
     * @param i_nameSize maximum size of device name buffer
     * @param o_computeCapability major and minor compute capability
     * @param o_maxThreadsPerBlock maximum threads per block on device
     * @return true if device info retrieved successfully
     **/
    static bool getDeviceInfo( char *o_deviceName,
                               size_t i_nameSize,
                               int o_computeCapability[2],
                               int &o_maxThreadsPerBlock );

    /**
     * Finalize CUDA device (cleanup).
     **/
    static void finalize();

    /**
     * Launch x-sweep kernel for wave propagation.
     * 
     * Parameters and behavior:
     * - Processes all horizontal edges in the grid
     * - One thread per edge for "one thread per edge" mapping
     * - Uses atomic operations to handle write conflicts
     * - Block size: 16x16 threads (256 threads per block)
     * - Grid size: computed to cover all edges
     * 
     * @param i_h water heights (read-only)
     * @param i_hu x-momentum (read-only)
     * @param i_hv y-momentum (read-only)
     * @param i_b bathymetry (read-only, constant)
     * @param o_h updated water heights
     * @param o_hu updated x-momentum
     * @param i_nCellsX number of cells in x-direction
     * @param i_nCellsY number of cells in y-direction
     * @param i_stride row stride in elements
     * @param i_scaling time step scaling (dt/dx)
     **/
    static void launchXSweepKernel( const t_real *i_h,
                                    const t_real *i_hu,
                                    const t_real *i_hv,
                                    const t_real *i_b,
                                    t_real *o_h,
                                    t_real *o_hu,
                                    t_idx i_nCellsX,
                                    t_idx i_nCellsY,
                                    t_idx i_stride,
                                    t_real i_scaling );

    /**
     * Launch y-sweep kernel for wave propagation.
     * 
     * Parameters and behavior:
     * - Processes all vertical edges in the grid
     * - One thread per edge for "one thread per edge" mapping
     * - Uses atomic operations to handle write conflicts
     * - Block size: 16x16 threads (256 threads per block)
     * - Grid size: computed to cover all edges
     * 
     * @param i_h water heights (read-only)
     * @param i_hu x-momentum (read-only)
     * @param i_hv y-momentum (read-only)
     * @param i_b bathymetry (read-only, constant)
     * @param o_h updated water heights
     * @param o_hv updated y-momentum
     * @param i_nCellsX number of cells in x-direction
     * @param i_nCellsY number of cells in y-direction
     * @param i_stride row stride in elements
     * @param i_scaling time step scaling (dt/dx)
     **/
    static void launchYSweepKernel( const t_real *i_h,
                                    const t_real *i_hu,
                                    const t_real *i_hv,
                                    const t_real *i_b,
                                    t_real *o_h,
                                    t_real *o_hv,
                                    t_idx i_nCellsX,
                                    t_idx i_nCellsY,
                                    t_idx i_stride,
                                    t_real i_scaling );

    /**
     * Allocate GPU memory for state arrays.
     * 
     * Allocates memory for h, hu, hv arrays with two time step buffers.
     * Also allocates memory for bathymetry array.
     * 
     * @param o_h GPU pointers to water height buffers (size 2)
     * @param o_hu GPU pointers to x-momentum buffers (size 2)
     * @param o_hv GPU pointers to y-momentum buffers (size 2)
     * @param o_b GPU pointer to bathymetry buffer
     * @param i_nValues total number of values per array
     * @return true if allocation successful
     **/
    static bool allocateMemory( t_real **o_h,
                                t_real **o_hu,
                                t_real **o_hv,
                                t_real **o_b,
                                t_idx i_nValues );

    /**
     * Free GPU memory previously allocated.
     * 
     * @param i_h GPU pointers to water height buffers (size 2)
     * @param i_hu GPU pointers to x-momentum buffers (size 2)
     * @param i_hv GPU pointers to y-momentum buffers (size 2)
     * @param i_b GPU pointer to bathymetry buffer
     **/
    static void freeMemory( t_real **i_h,
                            t_real **i_hu,
                            t_real **i_hv,
                            t_real *i_b );

    /**
     * Copy data from host to GPU.
     * 
     * @param i_h_host host array with water heights
     * @param o_h_device GPU array for water heights
     * @param i_nValues number of values to copy
     **/
    static void copyHostToDevice( const t_real *i_host,
                                  t_real *o_device,
                                  t_idx i_nValues );

    /**
     * Copy data from GPU to host.
     * 
     * @param i_device GPU array to read from
     * @param o_host host array to receive data
     * @param i_nValues number of values to copy
     **/
    static void copyDeviceToHost( const t_real *i_device,
                                  t_real *o_host,
                                  t_idx i_nValues );

    /**
     * Synchronize GPU to ensure all kernels complete.
     **/
    static void synchronize();

    /**
     * Check for CUDA errors and print message if error occurred.
     * 
     * @param i_message message to print if error detected
     * @return true if error occurred, false otherwise
     **/
    static bool checkError( const char *i_message );
};

#endif
