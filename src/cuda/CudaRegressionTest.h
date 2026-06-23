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
                             t_real const* i_h_gpu,
                             t_real const* i_hu_gpu,
                             t_real const* i_hv_gpu,
                             t_real const* i_b_gpu,
                             t_real*       o_h,
                             t_real*       o_hu,
                             t_real*       o_hv,
                             t_real*       o_b );

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
     * @brief Compare one value of a single cell between CPU and GPU results.
     *
     * @param i_x X index of the cell.
     * @param i_y Y index of the cell.
     * @param i_stride Number of cells in the x direction.
     * @param i_cpu Reference CPU array.
     * @param i_gpu Observed GPU array.
     * @param i_tolerance Maximum allowed absolute difference.
     * @return true if the cell values match within tolerance, false otherwise.
     */
    static bool compareSingleCellsHelper( t_idx         i_x,
                                    t_idx         i_y,
                                    t_idx         i_stride,
                                    t_real const* i_cpu,
                                    t_real const* i_gpu,
                                    t_real        i_tolerance = TOLERANCE );
    /**
     * @brief Compare all value of a single cell between CPU and GPU results.
     *
bool CudaRegressionTest::compareSingleCells( t_idx         i_x,
                                             t_idx         i_y,
                                             t_idx         i_stride,
                                             t_real const* i_h_cpu,
                                             t_real const* i_hu_cpu,
                                             t_real const* i_hv_cpu,
                                             t_real const* i_b_cpu,
                                             t_real const* i_h_gpu,
                                             t_real const* i_hu_gpu,
                                             t_real const* i_hv_gpu,
                                             t_real const* i_b_gpu,
                                             t_real        i_tolerance );

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
};

} // namespace cuda
} // namespace tsunami_lab

#endif