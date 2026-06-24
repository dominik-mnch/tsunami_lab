/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * CUDA kernels for two-dimensional wave propagation.
 **/

#include "../constants.h"
#include "../solvers/Roe.h"
#include "../solvers/F_wave.h"

// CUDA annotations for kernel functions
#ifdef __CUDACC__
  #define TSUNAMI_CUDA_GLOBAL __global__
#else
  #define TSUNAMI_CUDA_GLOBAL
#endif

namespace tsunami_lab {
  namespace patches {
    namespace cuda {
      /**
       * X-sweep kernel: computes net-updates over horizontal edges.
       * 
       * One thread per edge in x-direction.
       * Thread (bx, by) processes edge at column bx+1, row by+1.
       *
       * @param i_h_old water heights from previous time step
       * @param i_hu_old x-momentum from previous time step
       * @param i_hv_old y-momentum from previous time step
       * @param i_b bathymetry
       * @param o_h_new updated water heights
       * @param o_hu_new updated x-momentum
       * @param i_nCellsX number of interior cells in x-direction
       * @param i_nCellsY number of interior cells in y-direction
       * @param i_stride row stride in elements
       * @param i_scaling time step scaling factor (dt/dx)
       * @param i_useFWave if true, use F-Wave solver; otherwise Roe solver
       **/
      TSUNAMI_CUDA_GLOBAL
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
                   bool i_useFWave ) {

        // Get thread indices
        t_idx l_cx = blockIdx.x * blockDim.x + threadIdx.x;  // column (edge index)
        t_idx l_cy = blockIdx.y * blockDim.y + threadIdx.y;  // row

        // Check bounds: process edges for columns 0 to nCellsX (inclusive)
        // and rows 1 to nCellsY (interior only)
        if( l_cx > i_nCellsX || l_cy < 1 || l_cy > i_nCellsY ) {
          return;
        }

        // Compute linear indices
        t_idx l_ce  = l_cy * i_stride + l_cx;      // left cell
        t_idx l_ceR = l_ce + 1;                     // right cell

        // Check if cells are dry (bathymetry above water surface)
        bool l_leftDry = (i_b[l_ce] > 0);
        bool l_rightDry = (i_b[l_ceR] > 0);

        // Skip if both cells are dry
        if( l_leftDry && l_rightDry ) {
          return;
        }

        bool l_reflectAtShore = (l_leftDry != l_rightDry);

        // Load cell data
        t_real l_hL  = i_h_old[l_ce];
        t_real l_hR  = i_h_old[l_ceR];
        t_real l_huL = i_hu_old[l_ce];
        t_real l_huR = i_hu_old[l_ceR];
        t_real l_bL  = i_b[l_ce];
        t_real l_bR  = i_b[l_ceR];

        // Handle wet/dry interfaces: reflect at shore
        if( l_reflectAtShore ) {
          if( l_leftDry ) {
            l_hL  = l_hR;
            l_huL = -l_huR;
            l_bL  = l_bR;
          }
          else {
            l_hR  = l_hL;
            l_huR = -l_huL;
            l_bR  = l_bL;
          }
        }

        // Compute net updates using actual solver
        t_real l_netUpdates[2][2];

        if( i_useFWave ) {
          solvers::F_wave::netUpdates( l_hL,
                                       l_hR,
                                       l_huL,
                                       l_huR,
                                       l_bL,
                                       l_bR,
                                       l_netUpdates[0],
                                       l_netUpdates[1] );
        }
        else {
          solvers::Roe::netUpdates( l_hL,
                                    l_hR,
                                    l_huL,
                                    l_huR,
                                    l_netUpdates[0],
                                    l_netUpdates[1] );
        }

        t_real l_netUpdateL_h  = l_netUpdates[0][0];
        t_real l_netUpdateL_hu = l_netUpdates[0][1];
        t_real l_netUpdateR_h  = l_netUpdates[1][0];
        t_real l_netUpdateR_hu = l_netUpdates[1][1];

        // Handle reflecting boundary at shore
        if( l_reflectAtShore ) {
          if( l_leftDry ) {
            l_netUpdateL_h  = 0;
            l_netUpdateL_hu = 0;
          }
          else {
            l_netUpdateR_h  = 0;
            l_netUpdateR_hu = 0;
          }
        }

        // Atomically update cell values
        // Note: In a production implementation, we would use shared memory to reduce atomics
        atomicAdd( &o_h_new[l_ce],   -i_scaling * l_netUpdateL_h );
        atomicAdd( &o_hu_new[l_ce],  -i_scaling * l_netUpdateL_hu );
        atomicAdd( &o_h_new[l_ceR],  -i_scaling * l_netUpdateR_h );
        atomicAdd( &o_hu_new[l_ceR], -i_scaling * l_netUpdateR_hu );
      }

      /**
       * Y-sweep kernel: computes net-updates over vertical edges.
       * 
       * One thread per edge in y-direction.
       * Thread (bx, by) processes edge at column bx+1, row by.
       *
       * @param i_h_old water heights from previous time step
       * @param i_hu_old x-momentum from previous time step
       * @param i_hv_old y-momentum from previous time step
       * @param i_b bathymetry
       * @param o_h_new updated water heights
       * @param o_hv_new updated y-momentum
       * @param i_nCellsX number of interior cells in x-direction
       * @param i_nCellsY number of interior cells in y-direction
       * @param i_stride row stride in elements
       * @param i_scaling time step scaling factor (dt/dx)
       * @param i_useFWave if true, use F-Wave solver; otherwise Roe solver
       **/
      TSUNAMI_CUDA_GLOBAL
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
                   bool i_useFWave ) {

        // Get thread indices
        t_idx l_cx = blockIdx.x * blockDim.x + threadIdx.x;  // column
        t_idx l_cy = blockIdx.y * blockDim.y + threadIdx.y;  // row (edge index)

        // Check bounds: process edges for columns 1 to nCellsX (interior only)
        // and rows 0 to nCellsY (inclusive)
        if( l_cx < 1 || l_cx > i_nCellsX || l_cy > i_nCellsY ) {
          return;
        }

        // Compute linear indices
        t_idx l_ceB = l_cy * i_stride + l_cx;           // bottom cell
        t_idx l_ceT = (l_cy + 1) * i_stride + l_cx;     // top cell

        // Check if cells are dry
        bool l_bottomDry = (i_b[l_ceB] > 0);
        bool l_topDry = (i_b[l_ceT] > 0);

        // Skip if both cells are dry
        if( l_bottomDry && l_topDry ) {
          return;
        }

        bool l_reflectAtShore = (l_bottomDry != l_topDry);

        // Load cell data
        // Note: In y-sweep, hv acts like hu in the solver
        t_real l_hB  = i_h_old[l_ceB];
        t_real l_hT  = i_h_old[l_ceT];
        t_real l_hvB = i_hv_old[l_ceB];
        t_real l_hvT = i_hv_old[l_ceT];
        t_real l_bB  = i_b[l_ceB];
        t_real l_bT  = i_b[l_ceT];

        // Handle wet/dry interfaces
        if( l_reflectAtShore ) {
          if( l_bottomDry ) {
            l_hB  = l_hT;
            l_hvB = -l_hvT;
            l_bB  = l_bT;
          }
          else {
            l_hT  = l_hB;
            l_hvT = -l_hvB;
            l_bT  = l_bB;
          }
        }

        // Compute net updates using actual solver
        t_real l_netUpdates[2][2];

        if( i_useFWave ) {
          solvers::F_wave::netUpdates( l_hB,
                                       l_hT,
                                       l_hvB,
                                       l_hvT,
                                       l_bB,
                                       l_bT,
                                       l_netUpdates[0],
                                       l_netUpdates[1] );
        }
        else {
          solvers::Roe::netUpdates( l_hB,
                                    l_hT,
                                    l_hvB,
                                    l_hvT,
                                    l_netUpdates[0],
                                    l_netUpdates[1] );
        }

        t_real l_netUpdateB_h  = l_netUpdates[0][0];
        t_real l_netUpdateB_hv = l_netUpdates[0][1];
        t_real l_netUpdateT_h  = l_netUpdates[1][0];
        t_real l_netUpdateT_hv = l_netUpdates[1][1];

        // Handle reflecting boundary at shore
        if( l_reflectAtShore ) {
          if( l_bottomDry ) {
            l_netUpdateB_h  = 0;
            l_netUpdateB_hv = 0;
          }
          else {
            l_netUpdateT_h  = 0;
            l_netUpdateT_hv = 0;
          }
        }

        // Atomically update cell values
        atomicAdd( &o_h_new[l_ceB],  -i_scaling * l_netUpdateB_h );
        atomicAdd( &o_hv_new[l_ceB], -i_scaling * l_netUpdateB_hv );
        atomicAdd( &o_h_new[l_ceT],  -i_scaling * l_netUpdateT_h );
        atomicAdd( &o_hv_new[l_ceT], -i_scaling * l_netUpdateT_hv );
      }

      /**
       * Initialize new cell quantities from old values.
       *
       * @param i_h_old water heights from previous time step
       * @param i_hu_old x-momentum from previous time step
       * @param i_hv_old y-momentum from previous time step
       * @param o_h_new initialized water heights
       * @param o_hu_new initialized x-momentum
       * @param o_hv_new initialized y-momentum
       * @param i_nValues total number of values to initialize
       **/
      TSUNAMI_CUDA_GLOBAL
      void initNewCells( const t_real *i_h_old,
                        const t_real *i_hu_old,
                        const t_real *i_hv_old,
                        t_real *o_h_new,
                        t_real *o_hu_new,
                        t_real *o_hv_new,
                        t_idx i_nValues ) {

        t_idx l_idx = blockIdx.x * blockDim.x + threadIdx.x;

        if( l_idx < i_nValues ) {
          o_h_new[l_idx]  = i_h_old[l_idx];
          o_hu_new[l_idx] = i_hu_old[l_idx];
          o_hv_new[l_idx] = i_hv_old[l_idx];
        }
      }
    }
  }
}
