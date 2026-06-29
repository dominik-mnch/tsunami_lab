#ifndef CUDA_REGRESSION_TEST_H
#define CUDA_REGRESSION_TEST_H

#include "../constants.h"

/**
 * @file CudaRegressionTest.h
 * @brief Helper utilities for CUDA regression tests.
 */

namespace tsunami_lab {
namespace cuda {

/**
 * @brief Utilities for allocating, copying, comparing, and freeing CUDA test buffers.
 *
 * This class provides static helper functions used by CUDA regression tests to move
 * data between host and device memory, validate output, and release allocated GPU memory.
 */
class CudaRegressionTest {
public:
    /**
     * @brief Default tolerance for floating point comparisons.
     */
    static constexpr t_real TOLERANCE = static_cast<t_real>(1e-5f);

    /**
     * @brief Allocate device memory for the shallow water variables and bathymetry.
     *
     * @param i_nCellsX Number of cells in the x direction.
     * @param i_nCellsY Number of cells in the y direction.
     * @param o_h Output pointer for device height array.
     * @param o_hu Output pointer for device x-momentum array.
     * @param o_hv Output pointer for device y-momentum array.
     * @param o_b Output pointer for device bathymetry array.
     */
    static void allocateGPUMemory( t_idx i_nCellsX,
                                   t_idx i_nCellsY,
                                   t_real** o_h,
                                   t_real** o_hu,
                                   t_real** o_hv,
                                   t_real** o_b );

    /**
     * @brief Copy host data arrays to device memory.
     *
     * @param i_nCellsX Number of cells in the x direction.
     * @param i_nCellsY Number of cells in the y direction.
     * @param i_h Host height array.
     * @param i_hu Host x-momentum array.
     * @param i_hv Host y-momentum array.
     * @param i_b Host bathymetry array.
     * @param o_h_gpu Destination device height array.
     * @param o_hu_gpu Destination device x-momentum array.
     * @param o_hv_gpu Destination device y-momentum array.
     * @param o_b_gpu Destination device bathymetry array.
     */
    static void copyToGPU( t_idx          i_nCellsX,
                           t_idx          i_nCellsY,
                           t_real const* i_h,
                           t_real const* i_hu,
                           t_real const* i_hv,
                           t_real const* i_b,
                           t_real*       o_h_gpu,
                           t_real*       o_hu_gpu,
                           t_real*       o_hv_gpu,
                           t_real*       o_b_gpu );

    /**
     * @brief Copy device data arrays back to host memory.
     *
     * @param i_nCellsX Number of cells in the x direction.
     * @param i_nCellsY Number of cells in the y direction.
     * @param i_h_gpu Device height array.
     * @param i_hu_gpu Device x-momentum array.
     * @param i_hv_gpu Device y-momentum array.
     * @param i_b_gpu Device bathymetry array.
     * @param o_h Destination host height array.
     * @param o_hu Destination host x-momentum array.
     * @param o_hv Destination host y-momentum array.
     * @param o_b Destination host bathymetry array.
     */
static void copyFromGPU( t_idx          i_nCellsX,
                         t_idx          i_nCellsY,
                         t_real const*  i_h_gpu,
                         t_real const*  i_hu_gpu,
                         t_real const*  i_hv_gpu,
                         t_real*        o_h,
                         t_real*        o_hu,
                         t_real*        o_hv );

    /**
     * @brief Compare two 2D real grids with a tolerance.
     *
     * @param i_nCellsX Number of cells in the x direction.
     * @param i_nCellsY Number of cells in the y direction.
     * @param i_cpu Reference CPU array.
     * @param i_gpu Observed GPU array.
     * @param i_tolerance Maximum allowed absolute difference.
     * @return true if all grid values match within tolerance, false otherwise.
     */
    static bool compareGrids( t_idx          i_nCellsX,
                              t_idx          i_nCellsY,
                              t_real const* i_cpu,
                              t_real const* i_gpu,
                              t_real        i_tolerance = TOLERANCE );

    /**
     * @brief Compare one cell between CPU and GPU arrays.
     *
     * @param i_x X index of the cell.
     * @param i_y Y index of the cell.
     * @param i_nCellsX Number of cells in the x direction.
     * @param i_cpu Reference CPU array.
     * @param i_gpu Observed GPU array.
     * @param i_tolerance Maximum allowed absolute difference.
     * @return true if the cell values match within tolerance, false otherwise.
     */
    static bool compareSingleCells( t_idx         i_x,
                                    t_idx         i_y,
                                    t_idx         i_nCellsX,
                                    t_real const* i_cpu,
                                    t_real const* i_gpu,
                                    t_real        i_tolerance = TOLERANCE );

    /**
     * @brief Free device memory allocated for regression test arrays.
     *
     * @param i_h Device height array to free.
     * @param i_hu Device x-momentum array to free.
     * @param i_hv Device y-momentum array to free.
     * @param i_b Device bathymetry array to free.
     */
    static void freeGPUMemory( t_real* i_h,
                               t_real* i_hu,
                               t_real* i_hv,
                               t_real* i_b );

    /**
     * @brief Compare GPU and CPU solver results after one timestep (lock-free kernels).
     *
     * Runs a single x-sweep and y-sweep on both GPU and CPU, then compares results.
     *
     * @param i_nCellsX Number of cells in x direction.
     * @param i_nCellsY Number of cells in y direction.
     * @param i_useFWave Use F-Wave solver (true) or Roe solver (false).
     * @param o_maxError Optional output: maximum absolute error across all cells.
     * @return true if GPU and CPU results match within tolerance, false otherwise.
     */
    static bool compareKernelResults( t_idx i_nCellsX,
                                      t_idx i_nCellsY,
                                      bool i_useFWave,
                                      t_real* o_maxError = nullptr );

    /**
     * @brief Compare GPU and CPU solver results after one timestep (atomic kernels).
     *
     * Runs a single x-sweep and y-sweep using atomic operations on GPU and CPU, then compares results.
     *
     * @param i_nCellsX Number of cells in x direction.
     * @param i_nCellsY Number of cells in y direction.
     * @param i_useFWave Use F-Wave solver (true) or Roe solver (false).
     * @param o_maxError Optional output: maximum absolute error across all cells.
     * @return true if GPU and CPU results match within tolerance, false otherwise.
     */
    static bool compareKernelResultsAtomic( t_idx i_nCellsX,
                                            t_idx i_nCellsY,
                                            bool i_useFWave,
                                            t_real* o_maxError = nullptr );

    /**
     * @brief Compare GPU and CPU solver results over multiple timesteps (lock-free kernels).
     *
     * Runs both solvers for the specified number of timesteps and periodically
     * compares their results. Allows accumulated error from floating-point operations.
     *
     * @param i_nCellsX Number of cells in x direction.
     * @param i_nCellsY Number of cells in y direction.
     * @param i_nTimesteps Number of timesteps to run.
     * @param i_useFWave Use F-Wave solver (true) or Roe solver (false).
     * @param i_checkInterval Compare results every N timesteps (default: 1).
     * @param o_maxError Optional output: maximum error across all checks.
     * @return true if all comparisons are within tolerance, false otherwise.
     */
    static bool compareMultipleTimesteps( t_idx i_nCellsX,
                                          t_idx i_nCellsY,
                                          t_idx i_nTimesteps,
                                          bool i_useFWave,
                                          t_idx i_checkInterval = 1,
                                          t_real* o_maxError = nullptr );

    /**
     * @brief Compare GPU and CPU solver results over multiple timesteps (atomic kernels).
     *
     * Runs both solvers for the specified number of timesteps using atomic operations
     * and periodically compares their results. Allows accumulated error from floating-point operations.
     *
     * @param i_nCellsX Number of cells in x direction.
     * @param i_nCellsY Number of cells in y direction.
     * @param i_nTimesteps Number of timesteps to run.
     * @param i_useFWave Use F-Wave solver (true) or Roe solver (false).
     * @param i_checkInterval Compare results every N timesteps (default: 1).
     * @param o_maxError Optional output: maximum error across all checks.
     * @return true if all comparisons are within tolerance, false otherwise.
     */
    static bool compareMultipleTimestepsAtomic( t_idx i_nCellsX,
                                                t_idx i_nCellsY,
                                                t_idx i_nTimesteps,
                                                bool i_useFWave,
                                                t_idx i_checkInterval = 1,
                                                t_real* o_maxError = nullptr );

    /**
     * @brief Calculate maximum absolute error across all cells in three arrays.
     *
     * Useful for diagnostic output and error tracking.
     *
     * @param i_nCellsX Number of cells in x direction.
     * @param i_nCellsY Number of cells in y direction.
     * @param i_cpu Reference CPU array.
     * @param i_gpu Comparison GPU array.
     * @return Maximum absolute error found.
     */
    static t_real getMaxError( t_idx i_nCellsX,
                               t_idx i_nCellsY,
                               t_real const* i_cpu,
                               t_real const* i_gpu );
};

} // namespace cuda
} // namespace tsunami_lab

#endif