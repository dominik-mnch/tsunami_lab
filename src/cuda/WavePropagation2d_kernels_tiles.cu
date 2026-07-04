/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * CUDA kernels for tiled 32x32 wave propagation.
 **/

#include "../constants.h"
#include "../solvers/Roe.h"
#include "../solvers/F_wave.h"
#include "CudaLayout.h"

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
      using tsunami_lab::t_idx;
      using tsunami_lab::t_real;

      namespace {
        TSUNAMI_CUDA_DEVICE
        inline t_idx tileLinearIndex( t_idx i_row,
                                      t_idx i_col,
                                      t_idx i_stride,
                                      t_idx i_height,
                                      t_idx i_nTilesX,
                                      t_idx i_tileSize = 32 ) {
          const t_idx l_tileX = i_col / i_tileSize;
          const t_idx l_tileY = i_row / i_tileSize;
          const t_idx l_localX = i_col % i_tileSize;
          const t_idx l_localY = i_row % i_tileSize;
          const t_idx l_tileOffset = ((l_tileY * i_nTilesX + l_tileX) * i_tileSize * i_tileSize);
          return l_tileOffset + (l_localY * i_tileSize + l_localX);
        }

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

          bool l_aDry = (i_bA > 0) || (i_hA != i_hA);
          bool l_bDry = (i_bB > 0) || (i_hB != i_hB);

          if( l_aDry && l_bDry ) {
            return;
          }

          bool l_reflect = (l_aDry != l_bDry);
          t_real l_hA = i_hA, l_hB = i_hB;
          t_real l_quA = i_quA, l_quB = i_quB;
          t_real l_bA = i_bA, l_bB = i_bB;

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

          if( l_reflect ) {
            if( l_aDry ) { o_updA_h = 0; o_updA_q = 0; }
            else         { o_updB_h = 0; o_updB_q = 0; }
          }
        }
      }

      TSUNAMI_CUDA_GLOBAL
      void setGhostOutflowLeftRightKernelTile( t_real *io_h,
                                               t_real *io_hu,
                                               t_real *io_hv,
                                               t_real *io_b,
                                               t_idx i_nCellsX,
                                               t_idx i_nCellsY,
                                               t_idx i_stride,
                                               t_idx i_height,
                                               t_idx i_nTilesX,
                                               t_idx i_nTilesY ) {
        const t_idx l_row = blockIdx.y * blockDim.y + threadIdx.y;
        if( l_row >= i_height ) {
          return;
        }
        const t_idx l_col = blockIdx.x * blockDim.x + threadIdx.x;
        if( l_col >= i_stride ) {
          return;
        }
        const t_idx l_idx = tileLinearIndex( l_row, l_col, i_stride, i_height, i_nTilesX );
        if( l_col == 0 ) {
          const t_idx l_src = tileLinearIndex( l_row, 1, i_stride, i_height, i_nTilesX );
          io_h[l_idx] = io_h[l_src];
          io_hu[l_idx] = io_hu[l_src];
          io_hv[l_idx] = io_hv[l_src];
          io_b[l_idx] = io_b[l_src];
        }
        else if( l_col == i_nCellsX + 1 ) {
          const t_idx l_src = tileLinearIndex( l_row, i_nCellsX, i_stride, i_height, i_nTilesX );
          io_h[l_idx] = io_h[l_src];
          io_hu[l_idx] = io_hu[l_src];
          io_hv[l_idx] = io_hv[l_src];
          io_b[l_idx] = io_b[l_src];
        }
      }

      TSUNAMI_CUDA_GLOBAL
      void setGhostOutflowBottomTopKernelTile( t_real *io_h,
                                                t_real *io_hu,
                                                t_real *io_hv,
                                                t_real *io_b,
                                                t_idx i_nCellsX,
                                                t_idx i_nCellsY,
                                                t_idx i_stride,
                                                t_idx i_height,
                                                t_idx i_nTilesX,
                                                t_idx i_nTilesY ) {
        const t_idx l_row = blockIdx.y * blockDim.y + threadIdx.y;
        if( l_row >= i_height ) {
          return;
        }
        const t_idx l_col = blockIdx.x * blockDim.x + threadIdx.x;
        if( l_col >= i_stride ) {
          return;
        }
        const t_idx l_idx = tileLinearIndex( l_row, l_col, i_stride, i_height, i_nTilesX );
        if( l_row == 0 ) {
          const t_idx l_src = tileLinearIndex( 1, l_col, i_stride, i_height, i_nTilesX );
          io_h[l_idx] = io_h[l_src];
          io_hu[l_idx] = io_hu[l_src];
          io_hv[l_idx] = io_hv[l_src];
          io_b[l_idx] = io_b[l_src];
        }
        else if( l_row == i_nCellsY + 1 ) {
          const t_idx l_src = tileLinearIndex( i_nCellsY, l_col, i_stride, i_height, i_nTilesX );
          io_h[l_idx] = io_h[l_src];
          io_hu[l_idx] = io_hu[l_src];
          io_hv[l_idx] = io_hv[l_src];
          io_b[l_idx] = io_b[l_src];
        }
      }

      TSUNAMI_CUDA_GLOBAL
      void computeXSweepKernelTile( const t_real *i_h,
                                    const t_real *i_hu,
                                    const t_real *i_hv,
                                    const t_real *i_b,
                                    t_real *o_h,
                                    t_real *o_hu,
                                    t_real *o_hv,
                                    t_idx i_nCellsX,
                                    t_idx i_nCellsY,
                                    t_idx i_stride,
                                    t_idx i_height,
                                    t_idx i_nTilesX,
                                    t_idx i_nTilesY,
                                    t_real i_scaling,
                                    bool i_useFWave ) {
        const t_idx l_row = blockIdx.y * blockDim.y + threadIdx.y;
        const t_idx l_col = blockIdx.x * blockDim.x + threadIdx.x;
        if( l_row >= i_height || l_col >= i_stride ) {
          return;
        }
        const t_idx l_ce = tileLinearIndex( l_row, l_col, i_stride, i_height, i_nTilesX );
        if( l_col < 1 || l_col > i_nCellsX || l_row < 1 || l_row > i_nCellsY ) {
          o_h[l_ce] = i_h[l_ce];
          o_hu[l_ce] = i_hu[l_ce];
          o_hv[l_ce] = i_hv[l_ce];
          return;
        }

        t_real l_hNew = i_h[l_ce];
        t_real l_huNew = i_hu[l_ce];
        t_real l_uA_h, l_uA_q, l_uB_h, l_uB_q;

        const t_idx l_left = tileLinearIndex( l_row, l_col - 1, i_stride, i_height, i_nTilesX );
        edgeNetUpdate( i_h[l_left], i_h[l_ce],
                       i_hu[l_left], i_hu[l_ce],
                       i_b[l_left], i_b[l_ce],
                       i_useFWave,
                       l_uA_h, l_uA_q, l_uB_h, l_uB_q );
        l_hNew -= i_scaling * l_uB_h;
        l_huNew -= i_scaling * l_uB_q;

        const t_idx l_right = tileLinearIndex( l_row, l_col + 1, i_stride, i_height, i_nTilesX );
        edgeNetUpdate( i_h[l_ce], i_h[l_right],
                       i_hu[l_ce], i_hu[l_right],
                       i_b[l_ce], i_b[l_right],
                       i_useFWave,
                       l_uA_h, l_uA_q, l_uB_h, l_uB_q );
        l_hNew -= i_scaling * l_uA_h;
        l_huNew -= i_scaling * l_uA_q;

        o_h[l_ce] = l_hNew;
        o_hu[l_ce] = l_huNew;
        o_hv[l_ce] = i_hv[l_ce];
      }

      TSUNAMI_CUDA_GLOBAL
      void computeYSweepKernelTile( const t_real *i_h,
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
                                    bool i_useFWave ) {
        const t_idx l_row = blockIdx.y * blockDim.y + threadIdx.y;
        const t_idx l_col = blockIdx.x * blockDim.x + threadIdx.x;
        if( l_row >= i_height || l_col >= i_stride ) {
          return;
        }
        const t_idx l_ce = tileLinearIndex( l_row, l_col, i_stride, i_height, i_nTilesX );
        if( l_col < 1 || l_col > i_nCellsX || l_row < 1 || l_row > i_nCellsY ) {
          o_hv[l_ce] = i_hv[l_ce];
          return;
        }

        t_real l_uA_h, l_uA_q, l_uB_h, l_uB_q;
        const t_idx l_bottom = tileLinearIndex( l_row - 1, l_col, i_stride, i_height, i_nTilesX );
        edgeNetUpdate( i_h[l_bottom], i_h[l_ce],
                       i_hv[l_bottom], i_hv[l_ce],
                       i_b[l_bottom], i_b[l_ce],
                       i_useFWave,
                       l_uA_h, l_uA_q, l_uB_h, l_uB_q );
        io_h[l_ce] -= i_scaling * l_uB_h;
        t_real l_hvNew = i_hv[l_ce] - i_scaling * l_uB_q;

        const t_idx l_top = tileLinearIndex( l_row + 1, l_col, i_stride, i_height, i_nTilesX );
        edgeNetUpdate( i_h[l_ce], i_h[l_top],
                       i_hv[l_ce], i_hv[l_top],
                       i_b[l_ce], i_b[l_top],
                       i_useFWave,
                       l_uA_h, l_uA_q, l_uB_h, l_uB_q );
        io_h[l_ce] -= i_scaling * l_uA_h;
        l_hvNew -= i_scaling * l_uA_q;
        o_hv[l_ce] = l_hvNew;
      }

      TSUNAMI_CUDA_GLOBAL
      void computeXSweepAtomicKernelTile( const t_real *i_h,
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
                                          bool i_useFWave ) {
        const t_idx l_cx = blockIdx.x * blockDim.x + threadIdx.x;
        const t_idx l_row = blockIdx.y * blockDim.y + threadIdx.y;

        if( l_cx > i_nCellsX || l_row < 1 || l_row > i_nCellsY ) {
          return;
        }

        const t_idx l_left = tileLinearIndex( l_row, l_cx, i_stride, i_height, i_nTilesX );
        const t_idx l_right = tileLinearIndex( l_row, l_cx + 1, i_stride, i_height, i_nTilesX );

        bool l_leftDry = (i_b[l_left] > 0);
        bool l_rightDry = (i_b[l_right] > 0);

        if( l_leftDry && l_rightDry ) {
          return;
        }

        bool l_reflectAtShore = (l_leftDry != l_rightDry);

        t_real l_hL = i_h[l_left];
        t_real l_hR = i_h[l_right];
        t_real l_huL = i_hu[l_left];
        t_real l_huR = i_hu[l_right];
        t_real l_bL = i_b[l_left];
        t_real l_bR = i_b[l_right];

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

        atomicAdd( &io_h[l_left], -i_scaling * l_nu[0][0] );
        atomicAdd( &io_hu[l_left], -i_scaling * l_nu[0][1] );
        atomicAdd( &io_h[l_right], -i_scaling * l_nu[1][0] );
        atomicAdd( &io_hu[l_right], -i_scaling * l_nu[1][1] );
      }

      TSUNAMI_CUDA_GLOBAL
      void computeYSweepAtomicKernelTile( const t_real *i_h,
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
                                          bool i_useFWave ) {
        const t_idx l_col = blockIdx.x * blockDim.x + threadIdx.x;
        const t_idx l_row = blockIdx.y * blockDim.y + threadIdx.y;

        if( l_col < 1 || l_col > i_nCellsX || l_row > i_nCellsY ) {
          return;
        }

        const t_idx l_bottom = tileLinearIndex( l_row, l_col, i_stride, i_height, i_nTilesX );
        const t_idx l_top = tileLinearIndex( l_row + 1, l_col, i_stride, i_height, i_nTilesX );

        bool l_bottomDry = (i_b[l_bottom] > 0);
        bool l_topDry = (i_b[l_top] > 0);

        if( l_bottomDry && l_topDry ) {
          return;
        }

        bool l_reflectAtShore = (l_bottomDry != l_topDry);

        t_real l_hB = i_h[l_bottom];
        t_real l_hT = i_h[l_top];
        t_real l_hvB = i_hv[l_bottom];
        t_real l_hvT = i_hv[l_top];
        t_real l_bB = i_b[l_bottom];
        t_real l_bT = i_b[l_top];

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

        atomicAdd( &io_h[l_bottom], -i_scaling * l_nu[0][0] );
        atomicAdd( &io_hv[l_bottom], -i_scaling * l_nu[0][1] );
        atomicAdd( &io_h[l_top], -i_scaling * l_nu[1][0] );
        atomicAdd( &io_hv[l_top], -i_scaling * l_nu[1][1] );
      }

      TSUNAMI_CUDA_GLOBAL
      void initNewCellsKernelTile( const t_real *i_h,
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
                                   t_idx i_nTilesY ) {
        const t_idx l_row = blockIdx.y * blockDim.y + threadIdx.y;
        const t_idx l_col = blockIdx.x * blockDim.x + threadIdx.x;
        if( l_row >= i_height || l_col >= i_stride ) {
          return;
        }
        const t_idx l_idx = tileLinearIndex( l_row, l_col, i_stride, i_height, i_nTilesX );
        o_h[l_idx] = i_h[l_idx];
        o_hu[l_idx] = i_hu[l_idx];
        o_hv[l_idx] = i_hv[l_idx];
      }

      TSUNAMI_CUDA_GLOBAL
      void applyWetDryThresholdKernelTile( t_real *io_h,
                                            t_real *io_hu,
                                            t_real *io_hv,
                                            t_idx i_nCellsX,
                                            t_idx i_nCellsY,
                                            t_idx i_stride,
                                            t_idx i_height,
                                            t_idx i_nTilesX,
                                            t_idx i_nTilesY ) {
        const t_idx l_row = blockIdx.y * blockDim.y + threadIdx.y;
        const t_idx l_col = blockIdx.x * blockDim.x + threadIdx.x;
        if( l_row >= i_height || l_col >= i_stride ) {
          return;
        }
        const t_idx l_idx = tileLinearIndex( l_row, l_col, i_stride, i_height, i_nTilesX );
        if( io_h[l_idx] != io_h[l_idx] || io_hu[l_idx] != io_hu[l_idx] || io_hv[l_idx] != io_hv[l_idx] ) {
          io_h[l_idx] = 0.0;
          io_hu[l_idx] = 0.0;
          io_hv[l_idx] = 0.0;
        }
      }
    }
  }
}
