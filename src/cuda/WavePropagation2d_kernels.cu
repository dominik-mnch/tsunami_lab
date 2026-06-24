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
  #define TSUNAMI_CUDA_DEVICE __device__
#else
  #define TSUNAMI_CUDA_GLOBAL
  #define TSUNAMI_CUDA_DEVICE
#endif

namespace tsunami_lab {
  namespace patches {
    namespace cuda {
      /**
       * Compute the net-updates for a single edge between cells A and B,
       * replicating the wet/dry reflection logic of the CPU reference.
       *
       * The returned updates are already zeroed for a dry cell at a wet/dry
       * interface, so the caller can simply subtract them. If both cells are
       * dry, all four outputs are zero.
       *
       * @param i_hA water height of cell A.
       * @param i_hB water height of cell B.
       * @param i_quA momentum of cell A in the sweep direction.
       * @param i_quB momentum of cell B in the sweep direction.
       * @param i_bA bathymetry of cell A.
       * @param i_bB bathymetry of cell B.
       * @param i_useFWave if true, use the F-Wave solver; otherwise Roe.
       * @param o_updA_h net-update on the height of cell A.
       * @param o_updA_q net-update on the momentum of cell A.
       * @param o_updB_h net-update on the height of cell B.
       * @param o_updB_q net-update on the momentum of cell B.
       **/
      TSUNAMI_CUDA_DEVICE
      inline void edgeNetUpdate( t_real i_hA,
                                 t_real i_hB,
                                 t_real i_quA,
                                 t_real i_quB,
                                 t_real i_bA,
                                 t_real i_bB,
                                 bool i_useFWave,
                                 t_real &o_updA_h,
                                 t_real &o_updA_q,
                                 t_real &o_updB_h,
                                 t_real &o_updB_q ) {
        o_updA_h = 0; o_updA_q = 0;
        o_updB_h = 0; o_updB_q = 0;

        bool l_aDry = (i_bA > 0);
        bool l_bDry = (i_bB > 0);

        // Both dry: no flux across this edge.
        if( l_aDry && l_bDry ) {
          return;
        }

        bool l_reflect = (l_aDry != l_bDry);

        t_real l_hA = i_hA, l_hB = i_hB;
        t_real l_quA = i_quA, l_quB = i_quB;
        t_real l_bA = i_bA, l_bB = i_bB;

        // Mirror the wet state into the dry side to emulate a reflecting wall.
        if( l_reflect ) {
          if( l_aDry ) {
            l_hA = l_hB; l_quA = -l_quB; l_bA = l_bB;
          }
          else {
            l_hB = l_hA; l_quB = -l_quA; l_bB = l_bA;
          }
        }

        t_real l_nu[2][2];
        if( i_useFWave ) {
          solvers::F_wave::netUpdates( l_hA, l_hB, l_quA, l_quB, l_bA, l_bB,
                                       l_nu[0], l_nu[1] );
        }
        else {
          solvers::Roe::netUpdates( l_hA, l_hB, l_quA, l_quB,
                                    l_nu[0], l_nu[1] );
        }

        o_updA_h = l_nu[0][0]; o_updA_q = l_nu[0][1];
        o_updB_h = l_nu[1][0]; o_updB_q = l_nu[1][1];

        // Only update the wet cell at a wet/dry interface.
        if( l_reflect ) {
          if( l_aDry ) { o_updA_h = 0; o_updA_q = 0; }
          else         { o_updB_h = 0; o_updB_q = 0; }
        }
      }

      /**
       * Fused single-step kernel: one thread per cell. Each thread reads only
       * the OLD buffers and writes its NEW cell value exactly once, so there are
       * no atomics and no inter-kernel write hazards. Interior cells receive the
       * net-updates from their four surrounding edges (left/right in x, and
       * bottom/top in y, all evaluated from the old state, matching the CPU's
       * dimensional splitting). Ghost cells are copied through unchanged.
       *
       * The accumulation order (old - leftX - rightX - bottomY - topY) mirrors
       * the CPU reference so the floating-point results agree to round-off.
       *
       * @param i_h old water heights.
       * @param i_hu old x-momentum.
       * @param i_hv old y-momentum.
       * @param i_b bathymetry.
       * @param o_h new water heights.
       * @param o_hu new x-momentum.
       * @param o_hv new y-momentum.
       * @param i_nCellsX number of interior cells in x-direction.
       * @param i_nCellsY number of interior cells in y-direction.
       * @param i_stride row stride in elements (nCellsX + 2).
       * @param i_scaling time-step scaling dt/dx.
       * @param i_useFWave if true, use the F-Wave solver; otherwise Roe.
       **/
      TSUNAMI_CUDA_GLOBAL
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
                              bool i_useFWave ) {
        t_idx l_cx = blockIdx.x * blockDim.x + threadIdx.x;
        t_idx l_cy = blockIdx.y * blockDim.y + threadIdx.y;

        if( l_cx >= i_nCellsX + 2 || l_cy >= i_nCellsY + 2 ) {
          return;
        }

        t_idx l_ce = l_cy * i_stride + l_cx;

        // Ghost cells: copy through unchanged (they are refreshed by the
        // outflow boundary update at the start of the next step anyway).
        if( l_cx < 1 || l_cx > i_nCellsX || l_cy < 1 || l_cy > i_nCellsY ) {
          o_h[l_ce]  = i_h[l_ce];
          o_hu[l_ce] = i_hu[l_ce];
          o_hv[l_ce] = i_hv[l_ce];
          return;
        }

        t_real l_hNew  = i_h[l_ce];
        t_real l_huNew = i_hu[l_ce];
        t_real l_hvNew = i_hv[l_ce];

        t_real l_uA_h, l_uA_q, l_uB_h, l_uB_q;

        // Left x-edge (ce-1, ce): this cell is the right (B) cell.
        edgeNetUpdate( i_h[l_ce - 1], i_h[l_ce],
                       i_hu[l_ce - 1], i_hu[l_ce],
                       i_b[l_ce - 1], i_b[l_ce],
                       i_useFWave,
                       l_uA_h, l_uA_q, l_uB_h, l_uB_q );
        l_hNew  -= i_scaling * l_uB_h;
        l_huNew -= i_scaling * l_uB_q;

        // Right x-edge (ce, ce+1): this cell is the left (A) cell.
        edgeNetUpdate( i_h[l_ce], i_h[l_ce + 1],
                       i_hu[l_ce], i_hu[l_ce + 1],
                       i_b[l_ce], i_b[l_ce + 1],
                       i_useFWave,
                       l_uA_h, l_uA_q, l_uB_h, l_uB_q );
        l_hNew  -= i_scaling * l_uA_h;
        l_huNew -= i_scaling * l_uA_q;

        // Bottom y-edge (ce-stride, ce): this cell is the top (B) cell.
        edgeNetUpdate( i_h[l_ce - i_stride], i_h[l_ce],
                       i_hv[l_ce - i_stride], i_hv[l_ce],
                       i_b[l_ce - i_stride], i_b[l_ce],
                       i_useFWave,
                       l_uA_h, l_uA_q, l_uB_h, l_uB_q );
        l_hNew  -= i_scaling * l_uB_h;
        l_hvNew -= i_scaling * l_uB_q;

        // Top y-edge (ce, ce+stride): this cell is the bottom (A) cell.
        edgeNetUpdate( i_h[l_ce], i_h[l_ce + i_stride],
                       i_hv[l_ce], i_hv[l_ce + i_stride],
                       i_b[l_ce], i_b[l_ce + i_stride],
                       i_useFWave,
                       l_uA_h, l_uA_q, l_uB_h, l_uB_q );
        l_hNew  -= i_scaling * l_uA_h;
        l_hvNew -= i_scaling * l_uA_q;

        o_h[l_ce]  = l_hNew;
        o_hu[l_ce] = l_huNew;
        o_hv[l_ce] = l_hvNew;
      }
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
      void initNewCellsKernel( const t_real *i_h_old,
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

      /**
       * Apply outflow boundary conditions on the left and right edges: copy the
       * first/last interior column into the left/right ghost columns. One thread
       * per row (0 .. nCellsY+1). Writes only ghost columns, reads only interior
       * columns, so there is no data race.
       *
       * @param io_h water heights (modified in place).
       * @param io_hu x-momentum (modified in place).
       * @param io_hv y-momentum (modified in place).
       * @param i_nCellsX number of interior cells in x-direction.
       * @param i_nCellsY number of interior cells in y-direction.
       * @param i_stride row stride in elements (nCellsX + 2).
       **/
      TSUNAMI_CUDA_GLOBAL
      void setGhostOutflowLeftRightKernel( t_real *io_h,
                                           t_real *io_hu,
                                           t_real *io_hv,
                                           t_idx i_nCellsX,
                                           t_idx i_nCellsY,
                                           t_idx i_stride ) {

        t_idx l_idx = blockIdx.x * blockDim.x + threadIdx.x;
        if( l_idx >= i_nCellsY + 2 ) {
          return;
        }

        t_idx l_row = l_idx * i_stride;

        // left ghost <- first interior column
        io_h[l_row]  = io_h[l_row + 1];
        io_hu[l_row] = io_hu[l_row + 1];
        io_hv[l_row] = io_hv[l_row + 1];

        // right ghost <- last interior column
        io_h[l_row + i_nCellsX + 1]  = io_h[l_row + i_nCellsX];
        io_hu[l_row + i_nCellsX + 1] = io_hu[l_row + i_nCellsX];
        io_hv[l_row + i_nCellsX + 1] = io_hv[l_row + i_nCellsX];
      }

      /**
       * Apply outflow boundary conditions on the bottom and top edges: copy the
       * first/last interior row into the bottom/top ghost rows. One thread per
       * column (0 .. nCellsX+1). Must run AFTER the left/right kernel so that the
       * corner ghost cells take their value from the bottom/top edge, matching
       * the CPU reference implementation.
       *
       * @param io_h water heights (modified in place).
       * @param io_hu x-momentum (modified in place).
       * @param io_hv y-momentum (modified in place).
       * @param i_nCellsX number of interior cells in x-direction.
       * @param i_nCellsY number of interior cells in y-direction.
       * @param i_stride row stride in elements (nCellsX + 2).
       **/
      TSUNAMI_CUDA_GLOBAL
      void setGhostOutflowBottomTopKernel( t_real *io_h,
                                           t_real *io_hu,
                                           t_real *io_hv,
                                           t_idx i_nCellsX,
                                           t_idx i_nCellsY,
                                           t_idx i_stride ) {

        t_idx l_idx = blockIdx.x * blockDim.x + threadIdx.x;
        if( l_idx >= i_nCellsX + 2 ) {
          return;
        }

        // bottom ghost <- first interior row
        io_h[l_idx]  = io_h[i_stride + l_idx];
        io_hu[l_idx] = io_hu[i_stride + l_idx];
        io_hv[l_idx] = io_hv[i_stride + l_idx];

        // top ghost <- last interior row
        t_idx l_top      = (i_nCellsY + 1) * i_stride + l_idx;
        t_idx l_topInner = i_nCellsY * i_stride + l_idx;
        io_h[l_top]  = io_h[l_topInner];
        io_hu[l_top] = io_hu[l_topInner];
        io_hv[l_top] = io_hv[l_topInner];
      }
    }
  }
}
