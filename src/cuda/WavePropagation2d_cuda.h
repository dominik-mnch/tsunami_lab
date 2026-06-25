/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * GPU wrapper for two-dimensional wave propagation solver.
 **/

#ifndef TSUNAMI_LAB_PATCHES_WAVE_PROPAGATION_2D_CUDA_H
#define TSUNAMI_LAB_PATCHES_WAVE_PROPAGATION_2D_CUDA_H

#include "../constants.h"

namespace tsunami_lab {
  namespace patches {
    namespace cuda {
      class WavePropagation2dCuda;
    }
  }
}

/**
 * GPU-accelerated 2D wave propagation solver.
 * Manages GPU memory and provides interface for time stepping kernels.
 **/
class tsunami_lab::patches::cuda::WavePropagation2dCuda {

  private:
    //! GPU memory for water heights (two buffers for current/next step)
    t_real *m_d_h[2] = {nullptr, nullptr};

    //! GPU memory for x-momentum (two buffers for current/next step)
    t_real *m_d_hu[2] = {nullptr, nullptr};

    //! GPU memory for y-momentum (two buffers for current/next step)
    t_real *m_d_hv[2] = {nullptr, nullptr};

    //! GPU memory for bathymetry (constant throughout simulation)
    t_real *m_d_b = nullptr;

    //! Current step index (0 or 1) indicating active buffers
    unsigned short m_step = 0;

    //! Number of cells in x-direction
    t_idx m_nCellsX = 0;

    //! Number of cells in y-direction
    t_idx m_nCellsY = 0;

    //! Row stride (includes ghost cells)
    t_idx m_stride = 0;

    //! Total number of values per array
    t_idx m_nValues = 0;

    //! Whether to use F-Wave solver (true) or Roe solver (false)
    bool m_useFWaveSolver = false;

  public:
    /**
     * Initialize CUDA device and query capabilities.
     * 
     * @param i_deviceId GPU device ID to use (default: 0)
     * @return true if initialization successful, false otherwise
     **/
    static bool initialize( int i_deviceId = 0 );

    /**
     * Finalize CUDA device (cleanup).
     **/
    static void finalize();

    /**
     * Constructor: allocates GPU memory for rectangular domain.
     *
     * @param i_nCellsX number of cells in x-direction
     * @param i_nCellsY number of cells in y-direction
     * @param i_useFWaveSolver if true, uses F-Wave solver; otherwise Roe solver
     **/
    WavePropagation2dCuda( t_idx i_nCellsX,
                           t_idx i_nCellsY,
                           bool i_useFWaveSolver = false );

    /**
     * Destructor: frees GPU memory.
     **/
    ~WavePropagation2dCuda();

    /**
     * Copy host data to GPU memory.
     *
     * @param i_h host array with water heights
     * @param i_hu host array with x-momentum
     * @param i_hv host array with y-momentum
     * @param i_b host array with bathymetry
     **/
    void copyToGpu( const t_real *i_h,
                    const t_real *i_hu,
                    const t_real *i_hv,
                    const t_real *i_b );

    /**
     * Copy GPU results back to host memory.
     *
     * @param o_h host array to receive water heights
     * @param o_hu host array to receive x-momentum
     * @param o_hv host array to receive y-momentum
     **/
    void copyToHost( t_real *o_h,
                     t_real *o_hu,
                     t_real *o_hv );

    /**
     * Apply outflow boundary conditions on the active buffers by copying the
     * outermost interior cells into the surrounding ghost cells. Must be called
     * before computeStep() each time step to mirror the CPU solver.
     **/
    void setGhostOutflow();

    /**
     * Swap active buffers for next time step.
     **/
    void swapBuffers();

    /**
     * Perform one full time step in a single fused kernel launch: one thread per
     * cell reads only the old buffers and writes its new cell once, computing the
     * net-updates from all four surrounding edges. Race-free and deterministic.
     * Call swapBuffers() afterwards to make the result the new active state.
     *
     * @param i_scaling time-step scaling dt/dx.
     **/
    void computeStep( t_real i_scaling );

    /**
     * Get pointer to current active water height array on GPU.
     * @return GPU pointer to active h array
     **/
    t_real* getActiveH() { return m_d_h[m_step]; }

    /**
     * Get pointer to next water height array on GPU.
     * @return GPU pointer to next h array
     **/
    t_real* getNextH() { return m_d_h[1 - m_step]; }

    /**
     * Get stride (row width in elements).
     * @return stride value
     **/
    t_idx getStride() const { return m_stride; }

    /**
     * Get number of cells in x-direction.
     * @return nCellsX
     **/
    t_idx getNcellsX() const { return m_nCellsX; }

    /**
     * Get number of cells in y-direction.
     * @return nCellsY
     **/
    t_idx getNcellsY() const { return m_nCellsY; }

    /**
     * Get total number of values.
     * @return nValues
     **/
    t_idx getNValues() const { return m_nValues; }
};

#endif
