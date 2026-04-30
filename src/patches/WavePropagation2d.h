/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Two-dimensional wave propagation patch.
 **/
#ifndef TSUNAMI_LAB_PATCHES_WAVE_PROPAGATION_2D
#define TSUNAMI_LAB_PATCHES_WAVE_PROPAGATION_2D

#include "WavePropagation.h"

namespace tsunami_lab {
  namespace patches {
    class WavePropagation2d;
  }
}

class tsunami_lab::patches::WavePropagation2d: public WavePropagation {

  private:
    //! boolean that determines which solver to use
    bool m_useFWaveSolver;

    //! current step which indicates the active values in the arrays below
    unsigned short m_step = 0;

    //! number of cells discretizing the computational domain
    t_idx m_nCells = 0;

    //! water heights for the current and next time step for all cells
    t_real * m_h[2] = {nullptr, nullptr};

    //! momenta for the current and next time step for all cells in x-direction
    t_real * m_hu[2] = { nullptr, nullptr };

    //! momenta for the current and next time step for all cells in y-direction
    t_real * m_hv[2] = { nullptr, nullptr };

    //! bathymetry for all cells
    t_real * m_b = nullptr;

  public:
    /**
     * Constructs the 2d wave propagation solver.
     *
     * @param i_nCells number of cells.
     * @param i_useFWaveSolver if true, the class uses the f wave solver, otherwise it uses the Roe solver
     * @param i_boundaryCondition boundary condition for the ghost cells.
    **/
    WavePropagation2d( t_idx i_nCells, bool i_useFWaveSolver);

    /**
     * Destructor which frees all allocated memory.
     **/
    ~WavePropagation2d();

    /**
     * Performs a time step.
     *
     * @param i_scaling scaling of the time step (dt / dx) and (dt /).
     **/
    void timeStep( t_real i_scaling );

    /**
     * Sets the values of the ghost cells according to outflow boundary conditions.
     **/
    void setGhostOutflow();

    /**
     * Gets the stride in y-direction. x-direction is stride-1.
     *
     * @return stride in y-direction.
     **/
    t_idx getStride(){
      return m_nCells+2;
    }

    /**
     * Gets cells' water heights.
     *
     * @return water heights.
     */
    t_real const * getHeight(){
      return m_h[m_step]+1;
    }

    /**
     * Gets cells' bathymetry.
     *
     * @return bathymetry.
     */
    t_real const * getBathymetry(){
      return m_b;
    }

    /**
     * Gets the cells' momenta in x-direction.
     *
     * @return momenta in x-direction.
     **/
    t_real const * getMomentumX(){
      return m_hu[m_step]+1;
    }

    /**
     * Gets the cells' momenta in y-direction.
     *
     * @return momenta in y-direction.
     **/
    t_real const * getMomentumY(){
      return m_hv[m_step]+1;
    }

    /**
     * Sets the height of the cell to the given value.
     *
     * @param i_ix id of the cell in x-direction.
     * @param i_iy id of the cell in y-direction.
     * @param i_h water height.
     **/
    void setHeight( t_idx  i_ix,
                    t_idx  i_iy,
                    t_real i_h ) {
      m_h[m_step][i_iy * getStride() + i_ix + 1] = i_h;
    }

    /**
     * Sets the bathymetry of the cell to the given value.
     *
     * @param i_ix id of the cell in x-direction.
     * @param i_iy id of the cell in y-direction.
     * @param i_b bathymetry.
     **/
    void setBathymetry( t_idx  i_ix,
                        t_idx i_iy,
                        t_real i_b ) {
      m_b[i_iy * getStride() + i_ix + 1] = i_b;
    }

    /**
     * Sets the momentum in x-direction to the given value.
     *
     * @param i_ix id of the cell in x-direction.
     * @param i_iy id of the cell in y-direction.
     * @param i_hu momentum in x-direction.
     **/
    void setMomentumX( t_idx  i_ix,
                       t_idx  i_iy,
                       t_real i_hu ) {
      m_hu[m_step][i_iy * getStride() + i_ix + 1] = i_hu;
    }

    /**
     * Sets the momentum in y-direction to the given value.
     * 
     * @param i_ix id of the cell in x-direction.
     * @param i_iy id of the cell in y-direction.
     * @param i_hv momentum in y-direction.
     **/
    void setMomentumY( t_idx  i_ix,
                       t_idx  i_iy,
                       t_real i_hv ) {
      m_hv[m_step][i_iy * getStride() + i_ix + 1] = i_hv;
    }
};

#endif