/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Two-dimensional wave propagation patch.
 **/
#include "WavePropagation2d.h"
#include "../solvers/Roe.h"
#include "../solvers/F_wave.h"

tsunami_lab::patches::WavePropagation2d::WavePropagation2d( t_idx i_nCells,
                                                            bool i_useFWaveSolver,
                                                            BoundaryCondition i_boundaryCondition ):
  WavePropagation2d( i_nCells,
                     i_nCells,
                     i_useFWaveSolver,
                     i_boundaryCondition ) {
}

tsunami_lab::patches::WavePropagation2d::WavePropagation2d( t_idx i_nCellsX,
                                                            t_idx i_nCellsY,
                                                            bool  i_useFWaveSolver,
                                                            BoundaryCondition i_boundaryCondition ) {
  m_nCellsX = i_nCellsX;
  m_nCellsY = i_nCellsY;
  m_useFWaveSolver = i_useFWaveSolver;
  m_boundaryCondition = i_boundaryCondition;

  t_idx l_nValues = (m_nCellsY + 2) * getStride();

  // allocate memory including a single ghost cell on each side
  for( unsigned short l_st = 0; l_st < 2; l_st++ ) {
    m_h[l_st] = new t_real[ l_nValues ]{};
    m_hu[l_st] = new t_real[ l_nValues ]{};
    m_hv[l_st] = new t_real[ l_nValues ]{};
  }
  m_b = new t_real[ l_nValues ]{};
}

tsunami_lab::patches::WavePropagation2d::~WavePropagation2d() {
  for( unsigned short l_st = 0; l_st < 2; l_st++ ) {
    delete[] m_h[l_st];
    delete[] m_hu[l_st];
    delete[] m_hv[l_st];
  }
  delete[] m_b;
}

void tsunami_lab::patches::WavePropagation2d::timeStep( t_real i_scaling ) {
  t_idx l_stride = getStride();
  t_idx l_nValues = (m_nCellsY + 2) * l_stride;

  // pointers to old and new data
  t_real * l_hOld = m_h[m_step];
  t_real * l_huOld = m_hu[m_step];
  t_real * l_hvOld = m_hv[m_step];

  m_step = 1 - m_step;
  t_real * l_hNew =  m_h[m_step];
  t_real * l_huNew = m_hu[m_step];
  t_real * l_hvNew = m_hv[m_step];

  // use one parallel region to avoid creating threads multiple times
#pragma omp parallel
  {
    // init new cell quantities
#pragma omp for schedule(static)
    for( t_idx l_ce = 0; l_ce < l_nValues; l_ce++ ) {
      l_hNew[l_ce] = l_hOld[l_ce];
      l_huNew[l_ce] = l_huOld[l_ce];
      l_hvNew[l_ce] = l_hvOld[l_ce];
    }

// ensure all threads have finished initializing new cell quantities before starting to update them
#pragma omp barrier

    // X-sweep: net-updates over horizontal edges (between left and right cells).
    #pragma omp for schedule(static)
    for( t_idx l_cy = 1; l_cy <= m_nCellsY; l_cy++ ) {
      t_real l_netUpdates[2][2];
      for( t_idx l_cx = 0; l_cx <= m_nCellsX; l_cx++ ) {
        t_idx l_ce  = l_cy * l_stride + l_cx;
        t_idx l_ceR = l_ce + 1;

        bool l_leftDry = (m_b[l_ce] > 0);
        bool l_rightDry = (m_b[l_ceR] > 0);

        if( l_leftDry && l_rightDry ) {
          continue;
        }

        bool l_reflectAtShore = (l_leftDry != l_rightDry);

        t_real l_hL  = l_hOld[l_ce];
        t_real l_hR  = l_hOld[l_ceR];
        t_real l_huL = l_huOld[l_ce];
        t_real l_huR = l_huOld[l_ceR];
        t_real l_bL  = m_b[l_ce];
        t_real l_bR  = m_b[l_ceR];

        if( l_reflectAtShore ) {
          // Mirror wet state at wet/dry interfaces to emulate a reflecting wall.
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

        if( m_useFWaveSolver ) {
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

        if( l_reflectAtShore ) {
          // Only update the wet cell; keep the dry cell unchanged.
          if( l_leftDry ) {
            l_netUpdates[0][0] = 0;
            l_netUpdates[0][1] = 0;
          }
          else {
            l_netUpdates[1][0] = 0;
            l_netUpdates[1][1] = 0;
          }
        }

        l_hNew[l_ce]   -= i_scaling * l_netUpdates[0][0];
        l_huNew[l_ce]  -= i_scaling * l_netUpdates[0][1];
        l_hNew[l_ceR]  -= i_scaling * l_netUpdates[1][0];
        l_huNew[l_ceR] -= i_scaling * l_netUpdates[1][1];
      }
    }

// ensure X-sweep is complete before starting Y-sweep
#pragma omp barrier

    // Y-sweep: net-updates over vertical edges (between bottom and top cells).
    #pragma omp for schedule(static)
    for( t_idx l_cx = 1; l_cx <= m_nCellsX; l_cx++ ) {
      t_real l_netUpdates[2][2];
      for( t_idx l_cy = 0; l_cy <= m_nCellsY; l_cy++ ) {
        t_idx l_ceB = l_cy * l_stride + l_cx;
        t_idx l_ceT = (l_cy + 1) * l_stride + l_cx;

        bool l_bottomDry = (m_b[l_ceB] > 0);
        bool l_topDry = (m_b[l_ceT] > 0);

        if( l_bottomDry && l_topDry ) {
          continue;
        }

        bool l_reflectAtShore = (l_bottomDry != l_topDry);

        t_real l_hL  = l_hOld[l_ceB];
        t_real l_hR  = l_hOld[l_ceT];
        t_real l_huL = l_hvOld[l_ceB];
        t_real l_huR = l_hvOld[l_ceT];
        t_real l_bL  = m_b[l_ceB];
        t_real l_bR  = m_b[l_ceT];

        if( l_reflectAtShore ) {
          // Mirror wet state at wet/dry interfaces to emulate a reflecting wall.
          if( l_bottomDry ) {
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

        if( m_useFWaveSolver ) {
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

        if( l_reflectAtShore ) {
          // Only update the wet cell; keep the dry cell unchanged.
          if( l_bottomDry ) {
            l_netUpdates[0][0] = 0;
            l_netUpdates[0][1] = 0;
          }
          else {
            l_netUpdates[1][0] = 0;
            l_netUpdates[1][1] = 0;
          }
        }

        l_hNew[l_ceB]  -= i_scaling * l_netUpdates[0][0];
        l_hvNew[l_ceB] -= i_scaling * l_netUpdates[0][1];
        l_hNew[l_ceT]  -= i_scaling * l_netUpdates[1][0];
        l_hvNew[l_ceT] -= i_scaling * l_netUpdates[1][1];
      }
    }
  }
}

void tsunami_lab::patches::WavePropagation2d::setGhostOutflow() {
  t_idx l_stride = getStride();

  t_real * l_h = m_h[m_step];
  t_real * l_hu = m_hu[m_step];
  t_real * l_hv = m_hv[m_step];

  // read boundary conditions from bitmask
  unsigned short l_boundaryMask = static_cast<unsigned short>( m_boundaryCondition );
  bool l_reflectLeft = (l_boundaryMask & static_cast<unsigned short>( BoundaryCondition::BoundaryLeft )) != 0;
  bool l_reflectRight = (l_boundaryMask & static_cast<unsigned short>( BoundaryCondition::BoundaryRight )) != 0;
  bool l_reflectBottom = (l_boundaryMask & static_cast<unsigned short>( BoundaryCondition::BoundaryBottom )) != 0;
  bool l_reflectTop = (l_boundaryMask & static_cast<unsigned short>( BoundaryCondition::BoundaryTop )) != 0;

  // left and right boundary for all rows
  for( t_idx l_cy = 0; l_cy < m_nCellsY + 2; l_cy++ ) {
    t_idx l_row = l_cy * l_stride;

    l_h[l_row] = l_h[l_row + 1];
    if( l_reflectLeft ) {
      l_hu[l_row] = -l_hu[l_row + 1];
    }
    else {
      l_hu[l_row] = l_hu[l_row + 1];
    }
    l_hv[l_row] = l_hv[l_row + 1];
    m_b[l_row] = m_b[l_row + 1];

    l_h[l_row + m_nCellsX + 1] = l_h[l_row + m_nCellsX];
    if( l_reflectRight ) {
      l_hu[l_row + m_nCellsX + 1] = -l_hu[l_row + m_nCellsX];
    }
    else {
      l_hu[l_row + m_nCellsX + 1] = l_hu[l_row + m_nCellsX];
    }
    l_hv[l_row + m_nCellsX + 1] = l_hv[l_row + m_nCellsX];
    m_b[l_row + m_nCellsX + 1] = m_b[l_row + m_nCellsX];
  }

  // bottom and top boundary for all columns
  for( t_idx l_cx = 0; l_cx < m_nCellsX + 2; l_cx++ ) {
    l_h[l_cx] = l_h[l_stride + l_cx];
    l_hu[l_cx] = l_hu[l_stride + l_cx];
    if( l_reflectBottom ) {
      l_hv[l_cx] = -l_hv[l_stride + l_cx];
    }
    else {
      l_hv[l_cx] = l_hv[l_stride + l_cx];
    }
    m_b[l_cx] = m_b[l_stride + l_cx];

    t_idx l_top = (m_nCellsY + 1) * l_stride + l_cx;
    t_idx l_topInner = m_nCellsY * l_stride + l_cx;
    l_h[l_top] = l_h[l_topInner];
    l_hu[l_top] = l_hu[l_topInner];
    if( l_reflectTop ) {
      l_hv[l_top] = -l_hv[l_topInner];
    }
    else {
      l_hv[l_top] = l_hv[l_topInner];
    }
    m_b[l_top] = m_b[l_topInner];
  }
}