/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * CUDA kernels for two-dimensional wave propagation.
 **/

#include "../constants.h"
#include "../solvers/Roe.h"
#include "../solvers/F_wave.h"
#include "CudaLayout.h"

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

        // Check bathymetry-based dryness AND catch NaN values
        bool l_aDry = (i_bA > 0) || (i_hA != i_hA);
        bool l_bDry = (i_bB > 0) || (i_hB != i_hB);

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
       * @param i_height column height in elements (nCellsX + 2).
       **/
      TSUNAMI_CUDA_GLOBAL
      void setGhostOutflowLeftRightKernelColumnMajor( t_real *io_h,
                                                      t_real *io_hu,
                                                      t_real *io_hv,
                                                      t_idx i_nCellsX,
                                                      t_idx i_nCellsY,
                                                      t_idx i_height ) {

        t_idx l_row = blockIdx.x * blockDim.x + threadIdx.x;
        if( l_row >= i_nCellsY + 2 ) {
          return;
        }

        // left ghost <- first interior column
        io_h[l_row]  = io_h[i_height + l_row];
        io_hu[l_row] = io_hu[i_height + l_row];
        io_hv[l_row] = io_hv[i_height + l_row];

        // right ghost <- last interior column
        const t_idx l_rightGhost = (i_nCellsX + 1) * i_height + l_row;
        const t_idx l_rightInner = i_nCellsX * i_height + l_row;
        io_h[l_rightGhost]  = io_h[l_rightInner];
        io_hu[l_rightGhost] = io_hu[l_rightInner];
        io_hv[l_rightGhost] = io_hv[l_rightInner];
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
      void setGhostOutflowBottomTopKernelColumnMajor( t_real *io_h,
                                                       t_real *io_hu,
                                                       t_real *io_hv,
                                                       t_real *io_b,
                                                       t_idx i_nCellsX,
                                                       t_idx i_nCellsY,
                                                       t_idx i_height ) {

        t_idx l_col = blockIdx.x * blockDim.x + threadIdx.x;
        if( l_col >= i_nCellsX + 2 ) {
          return;
        }

        const t_idx l_base = l_col * i_height;

        // bottom ghost <- first interior row
        io_h[l_base]  = io_h[l_base + 1];
        io_hu[l_base] = io_hu[l_base + 1];
        io_hv[l_base] = io_hv[l_base + 1];

        // top ghost <- last interior row
        const t_idx l_top = l_base + i_nCellsY + 1;
        const t_idx l_topInner = l_base + i_nCellsY;
        io_h[l_top]  = io_h[l_topInner];
        io_hu[l_top] = io_hu[l_topInner];
        io_hv[l_top] = io_hv[l_topInner];
      }

      /**
       * X-sweep kernel: one thread per cell.
       *
       * Reads the OLD buffers. For ghost cells, copies all three fields
       * through unchanged. For interior cells, applies x-direction net-updates
       * to h and hu, and copies hv unchanged to the NEW buffer.
       * Mirrors the x-sweep pass of the CPU operator-split reference.
       **/
      TSUNAMI_CUDA_GLOBAL
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
                                            bool i_useFWave ) {
        t_idx l_cx = blockIdx.x * blockDim.x + threadIdx.x;
        t_idx l_cy = blockIdx.y * blockDim.y + threadIdx.y;

        if( l_cx >= i_nCellsX + 2 || l_cy >= i_nCellsY + 2 ) return;

        t_idx l_ce = COLUMN_MAJOR(l_cy, l_cx, i_height);

        // Ghost cells: copy through unchanged.
        if( l_cx < 1 || l_cx > i_nCellsX || l_cy < 1 || l_cy > i_nCellsY ) {
          o_h[l_ce]  = i_h[l_ce];
          o_hu[l_ce] = i_hu[l_ce];
          o_hv[l_ce] = i_hv[l_ce];
          return;
        }

        t_real l_hNew  = i_h[l_ce];
        t_real l_huNew = i_hu[l_ce];

        t_real l_uA_h, l_uA_q, l_uB_h, l_uB_q;

        // Left x-edge: neighbor is in the previous column.
        const t_idx l_left = l_ce - i_height;
        edgeNetUpdate( i_h[l_left], i_h[l_ce],
                       i_hu[l_left], i_hu[l_ce],
                       i_b[l_left], i_b[l_ce],
                       i_useFWave,
                       l_uA_h, l_uA_q, l_uB_h, l_uB_q );
        l_hNew  -= i_scaling * l_uB_h;
        l_huNew -= i_scaling * l_uB_q;

        // Right x-edge: neighbor is in the next column.
        const t_idx l_right = l_ce + i_height;
        edgeNetUpdate( i_h[l_ce], i_h[l_right],
                       i_hu[l_ce], i_hu[l_right],
                       i_b[l_ce], i_b[l_right],
                       i_useFWave,
                       l_uA_h, l_uA_q, l_uB_h, l_uB_q );
        l_hNew  -= i_scaling * l_uA_h;
        l_huNew -= i_scaling * l_uA_q;

        o_h[l_ce]  = l_hNew;
        o_hu[l_ce] = l_huNew;
        o_hv[l_ce] = i_hv[l_ce];
      }

      /**
       * Y-sweep kernel: one thread per cell.
       *
       * Uses the OLD h and hv buffers for Riemann-problem inputs (matching the
       * CPU reference which reads h_old in the y-sweep). Accumulates the
       * y-direction net-update on top of the x-swept h already in io_h, and
       * writes the final hv value to o_hv. hu is left untouched (already
       * written by the x-sweep kernel).
       *
       * No atomics are required: each thread is the sole writer of its cell.
       **/
      TSUNAMI_CUDA_GLOBAL
      void computeYSweepKernelColumnMajor( const t_real *i_h,
                                           const t_real *i_hv,
                                           const t_real *i_b,
                                           t_real *io_h,
                                           t_real *o_hv,
                                           t_idx i_nCellsX,
                                           t_idx i_nCellsY,
                                           t_idx i_height,
                                           t_real i_scaling,
                                           bool i_useFWave ) {
        t_idx l_cx = blockIdx.x * blockDim.x + threadIdx.x;
        t_idx l_cy = blockIdx.y * blockDim.y + threadIdx.y;

        if( l_cx >= i_nCellsX + 2 || l_cy >= i_nCellsY + 2 ) return;

        t_idx l_ce = COLUMN_MAJOR(l_cy, l_cx, i_height);

        // Ghost cells: copy hv through (h was already handled by x-sweep).
        if( l_cx < 1 || l_cx > i_nCellsX || l_cy < 1 || l_cy > i_nCellsY ) {
          o_hv[l_ce] = i_hv[l_ce];
          return;
        }

        t_real l_uA_h, l_uA_q, l_uB_h, l_uB_q;

        // Bottom y-edge: neighbor is in the previous row.
        const t_idx l_bottom = l_ce - 1;
        edgeNetUpdate( i_h[l_bottom], i_h[l_ce],
                       i_hv[l_bottom], i_hv[l_ce],
                       i_b[l_bottom], i_b[l_ce],
                       i_useFWave,
                       l_uA_h, l_uA_q, l_uB_h, l_uB_q );
        io_h[l_ce] -= i_scaling * l_uB_h;
        t_real l_hvNew = i_hv[l_ce] - i_scaling * l_uB_q;

        // Top y-edge: neighbor is in the next row.
        const t_idx l_top = l_ce + 1;
        edgeNetUpdate( i_h[l_ce], i_h[l_top],
                       i_hv[l_ce], i_hv[l_top],
                       i_b[l_ce], i_b[l_top],
                       i_useFWave,
                       l_uA_h, l_uA_q, l_uB_h, l_uB_q );
        io_h[l_ce] -= i_scaling * l_uA_h;
        l_hvNew    -= i_scaling * l_uA_q;

        o_hv[l_ce] = l_hvNew;
      }

      /**
       * Initialize new cell quantities by copying from old buffers.
       * One thread per cell.
       *
       * @param i_h old water heights
       * @param i_hu old x-momentum
       * @param i_hv old y-momentum
       * @param o_h new water heights
       * @param o_hu new x-momentum
       * @param o_hv new y-momentum
       * @param i_nCellsX number of interior cells in x-direction
       * @param i_nCellsY number of interior cells in y-direction
       * @param i_height column height in elements
       **/
      TSUNAMI_CUDA_GLOBAL
      void initNewCellsKernelColumnMajor( const t_real *i_h,
                                          const t_real *i_hu,
                                          const t_real *i_hv,
                                          t_real *o_h,
                                          t_real *o_hu,
                                          t_real *o_hv,
                                          t_idx i_nCellsX,
                                          t_idx i_nCellsY,
                                          t_idx i_height ) {
        t_idx l_cx = blockIdx.x * blockDim.x + threadIdx.x;
        t_idx l_cy = blockIdx.y * blockDim.y + threadIdx.y;

        if( l_cx >= i_nCellsX + 2 || l_cy >= i_nCellsY + 2 ) return;

        t_idx l_ce = COLUMN_MAJOR(l_cy, l_cx, i_height);
        o_h[l_ce]  = i_h[l_ce];
        o_hu[l_ce] = i_hu[l_ce];
        o_hv[l_ce] = i_hv[l_ce];
      }

      /**
       * X-sweep kernel with atomic updates: one thread per edge.
       *
       * Computes net-updates from the Riemann solver and uses atomicAdd to
       * accumulate them into the new buffers. Similar structure to CPU version.
       *
       * @param i_h old water heights
       * @param i_hu old x-momentum
       * @param i_hv old y-momentum
       * @param i_b bathymetry
       * @param io_h new water heights (modified in place with atomicAdd)
       * @param io_hu new x-momentum (modified in place with atomicAdd)
       * @param i_nCellsX number of interior cells in x-direction
       * @param i_nCellsY number of interior cells in y-direction
       * @param i_height column height in elements
       * @param i_scaling time-step scaling dt/dx
       * @param i_useFWave whether to use F-Wave solver
       **/
      TSUNAMI_CUDA_GLOBAL
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
                                                 bool i_useFWave ) {
        // One thread per x-edge: cy in [1, nCellsY], cx in [0, nCellsX]
        t_idx l_cx = blockIdx.x * blockDim.x + threadIdx.x;  // edge index 0..nCellsX
        t_idx l_cy = blockIdx.y * blockDim.y + threadIdx.y;  // row index 1..nCellsY

        if( l_cx > i_nCellsX || l_cy < 1 || l_cy > i_nCellsY ) return;

        t_idx l_ce = COLUMN_MAJOR(l_cy, l_cx, i_height);   // left cell
        t_idx l_ceR = l_ce + i_height;               // right cell

        bool l_leftDry = (i_b[l_ce] > 0);
        bool l_rightDry = (i_b[l_ceR] > 0);

        if( l_leftDry && l_rightDry ) {
          return;
        }

        bool l_reflectAtShore = (l_leftDry != l_rightDry);

        t_real l_hL = i_h[l_ce];
        t_real l_hR = i_h[l_ceR];
        t_real l_huL = i_hu[l_ce];
        t_real l_huR = i_hu[l_ceR];
        t_real l_bL = i_b[l_ce];
        t_real l_bR = i_b[l_ceR];

        // Mirror at wet/dry interface
        if( l_reflectAtShore ) {
          if( l_leftDry ) {
            l_hL = l_hR;
            l_huL = -l_huR;
            l_bL = l_bR;
          }
          else {
            l_hR = l_hL;
            l_huR = -l_huL;
            l_bR = l_bL;
          }
        }

        t_real l_nu[2][2];
        if( i_useFWave ) {
          solvers::F_wave::netUpdates( l_hL, l_hR, l_huL, l_huR, l_bL, l_bR,
                                       l_nu[0], l_nu[1] );
        }
        else {
          solvers::Roe::netUpdates( l_hL, l_hR, l_huL, l_huR,
                                    l_nu[0], l_nu[1] );
        }

        if( l_reflectAtShore ) {
          if( l_leftDry ) {
            l_nu[0][0] = 0;
            l_nu[0][1] = 0;
          }
          else {
            l_nu[1][0] = 0;
            l_nu[1][1] = 0;
          }
        }

        // Use atomicAdd to update both cells
        atomicAdd( &io_h[l_ce], -i_scaling * l_nu[0][0] );
        atomicAdd( &io_hu[l_ce], -i_scaling * l_nu[0][1] );
        atomicAdd( &io_h[l_ceR], -i_scaling * l_nu[1][0] );
        atomicAdd( &io_hu[l_ceR], -i_scaling * l_nu[1][1] );
      }

      /**
       * Y-sweep kernel with atomic updates: one thread per edge.
       *
       * Computes net-updates from the Riemann solver and uses atomicAdd to
       * accumulate them into the new buffers. Similar structure to CPU version.
       *
       * @param i_h old water heights
       * @param i_hv old y-momentum
       * @param i_b bathymetry
       * @param io_h new water heights (modified in place with atomicAdd)
       * @param io_hv new y-momentum (modified in place with atomicAdd)
       * @param i_nCellsX number of interior cells in x-direction
       * @param i_nCellsY number of interior cells in y-direction
       * @param i_height column height in elements
       * @param i_scaling time-step scaling dt/dx
       * @param i_useFWave whether to use F-Wave solver
       **/
      TSUNAMI_CUDA_GLOBAL
      void computeYSweepAtomicKernelColumnMajor( const t_real *i_h,
                                                 const t_real *i_hv,
                                                 const t_real *i_b,
                                                 t_real *io_h,
                                                 t_real *io_hv,
                                                 t_idx i_nCellsX,
                                                 t_idx i_nCellsY,
                                                 t_idx i_height,
                                                 t_real i_scaling,
                                                 bool i_useFWave ) {
        // One thread per y-edge: cx in [1, nCellsX], cy in [0, nCellsY]
        t_idx l_cx = blockIdx.x * blockDim.x + threadIdx.x;  // column index 1..nCellsX
        t_idx l_cy = blockIdx.y * blockDim.y + threadIdx.y;  // edge index 0..nCellsY

        if( l_cx < 1 || l_cx > i_nCellsX || l_cy > i_nCellsY ) return;

        t_idx l_ceB = COLUMN_MAJOR(l_cy, l_cx, i_height);      // bottom cell
        t_idx l_ceT = COLUMN_MAJOR(l_cy + 1, l_cx, i_height); // top cell

        bool l_bottomDry = (i_b[l_ceB] > 0);
        bool l_topDry = (i_b[l_ceT] > 0);

        if( l_bottomDry && l_topDry ) {
          return;
        }

        bool l_reflectAtShore = (l_bottomDry != l_topDry);

        t_real l_hB = i_h[l_ceB];
        t_real l_hT = i_h[l_ceT];
        t_real l_hvB = i_hv[l_ceB];
        t_real l_hvT = i_hv[l_ceT];
        t_real l_bB = i_b[l_ceB];
        t_real l_bT = i_b[l_ceT];

        // Mirror at wet/dry interface
        if( l_reflectAtShore ) {
          if( l_bottomDry ) {
            l_hB = l_hT;
            l_hvB = -l_hvT;
            l_bB = l_bT;
          }
          else {
            l_hT = l_hB;
            l_hvT = -l_hvB;
            l_bT = l_bB;
          }
        }

        t_real l_nu[2][2];
        if( i_useFWave ) {
          solvers::F_wave::netUpdates( l_hB, l_hT, l_hvB, l_hvT, l_bB, l_bT,
                                       l_nu[0], l_nu[1] );
        }
        else {
          solvers::Roe::netUpdates( l_hB, l_hT, l_hvB, l_hvT,
                                    l_nu[0], l_nu[1] );
        }

        if( l_reflectAtShore ) {
          if( l_bottomDry ) {
            l_nu[0][0] = 0;
            l_nu[0][1] = 0;
          }
          else {
            l_nu[1][0] = 0;
            l_nu[1][1] = 0;
          }
        }

        // Use atomicAdd to update both cells
        atomicAdd( &io_h[l_ceB], -i_scaling * l_nu[0][0] );
        atomicAdd( &io_hv[l_ceB], -i_scaling * l_nu[0][1] );
        atomicAdd( &io_h[l_ceT], -i_scaling * l_nu[1][0] );
        atomicAdd( &io_hv[l_ceT], -i_scaling * l_nu[1][1] );
      }

      /**
       * Post-processing kernel: Clear NaN values in state variables.
       * Cells with NaN in any state variable (h, hu, hv) are zeroed to prevent
       * NaN propagation through subsequent timesteps.
       *
       * @param io_h device height array (full grid with ghost cells)
       * @param io_hu device x-momentum array
       * @param io_hv device y-momentum array
       * @param i_nCellsX number of interior cells in x-direction
       * @param i_nCellsY number of interior cells in y-direction
       * @param i_height column height in elements
       **/
      TSUNAMI_CUDA_GLOBAL
      void applyWetDryThresholdKernelColumnMajor( t_real *io_h,
                                                   t_real *io_hu,
                                                   t_real *io_hv,
                                                   t_idx i_nCellsX,
                                                   t_idx i_nCellsY,
                                                   t_idx i_height ) {
        // One thread per cell (including ghost cells)
        t_idx l_cx = blockIdx.x * blockDim.x + threadIdx.x;
        t_idx l_cy = blockIdx.y * blockDim.y + threadIdx.y;

        if( l_cx >= i_nCellsX + 2 || l_cy >= i_nCellsY + 2 ) return;

        t_idx l_idx = COLUMN_MAJOR( l_cy, l_cx, i_height );

        // Clear cells with NaN values to prevent propagation
        if( io_h[l_idx] != io_h[l_idx] ||  // NaN check for h
            io_hu[l_idx] != io_hu[l_idx] ||  // NaN check for hu
            io_hv[l_idx] != io_hv[l_idx] ) {  // NaN check for hv
          io_h[l_idx] = 0.0;
          io_hu[l_idx] = 0.0;
          io_hv[l_idx] = 0.0;
        }
      }
    }
  }
}
