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

tsunami_lab::patches::WavePropagation2d::WavePropagation2d( t_idx i_nCells , bool i_useFWaveSolver) {
  m_nCells = i_nCells;
  m_useFWaveSolver = i_useFWaveSolver;

  // allocate memory including a single ghost cell on each side
  for( unsigned short l_st = 0; l_st < 2; l_st++ ) {
    m_h[l_st] = new t_real[  getStride() * getStride() ];
    m_hu[l_st] = new t_real[ getStride() ];
    m_hv[l_st] = new t_real[ getStride() ];
  }
  m_b = new t_real[ (getStride()) * (getStride()) ];

  // init to zero
  for( unsigned short l_st = 0; l_st < 2; l_st++ ) {
    for( t_idx l_cx = 0; l_cx < getStride(); l_cx++ ) {
      m_hu[l_st][l_cx] = 0;
      m_hv[l_st][l_cx] = 0;
      for ( t_idx l_cy = 0; l_cy < getStride(); l_cy++ ) {
        m_h[l_st][l_cy * getStride() + l_cx + 1] = 0;
      }
    }
  }
  for( t_idx l_ce = 0; l_ce < getStride(); l_ce++ ) {
    m_b[l_ce] = 0;
  }
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
  // pointers to old and new data
  t_real * l_hOld = m_h[m_step];
  t_real * l_huOld = m_hu[m_step];
  t_real * l_hvOld = m_hv[m_step];

  m_step = (m_step+1) % 2;
  t_real * l_hNew =  m_h[m_step];
  t_real * l_huNew = m_hu[m_step];
  t_real * l_hvNew = m_hv[m_step];

  // init new cell quantities
  for( t_idx l_cx = 1; l_cx < m_nCells+1; l_cx++ ) {
    l_huNew[l_cx] = l_huOld[l_cx];
    l_hvNew[l_cx] = l_hvOld[l_cx];
    for ( t_idx l_cy = 1; l_cy < m_nCells+1; l_cy++ ) {
        l_hNew[l_cy * getStride() + l_cx + 1] = l_hOld[l_cy * getStride() + l_cx + 1];
    }
  }

  // iterate over edges and update with Riemann solutions
  for( t_idx l_cx = 0; l_cx < m_nCells + 1; l_cx++ ) {
    for ( t_idx l_cy = 0; l_cy < m_nCells + 1; l_cy++ ) {

    // !left and right update
    // determine left, right and top cell-id
    t_idx l_ce = l_cy * getStride() + l_cx;
    t_idx l_ceR = l_cy * getStride() + l_cx + 1;
    t_idx l_ceT = (l_cy + 1) * getStride() + l_cx;

    t_real l_hL = l_hOld[l_ce];
    t_real l_hR = l_hOld[l_ceR];
    t_real l_huL = l_huOld[l_ce];
    t_real l_huR = l_huOld[l_ceR];
    t_real l_bL = m_b[l_ce];
    t_real l_bR = m_b[l_ceR];


    // compute net-updates
    t_real l_netUpdates[2][2];

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

    // update the cells' quantities
    l_hNew[l_ce]  -= i_scaling * l_netUpdates[0][0];
    l_huNew[l_ce] -= i_scaling * l_netUpdates[0][1];

    l_hNew[l_ceR]  -= i_scaling * l_netUpdates[1][0];
    l_huNew[l_ceR] -= i_scaling * l_netUpdates[1][1];


    // !top and bottom update
    t_real l_hB = l_hOld[l_ce];
    t_real l_hT = l_hOld[l_ceT];
    t_real l_hvB = l_hvOld[l_ce];
    t_real l_hvT = l_hvOld[l_ceT];
    t_real l_bB = m_b[l_ce];
    t_real l_bT = m_b[l_ceT];

    // compute net updates
    if( m_useFWaveSolver ) {
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
                                    l_hvT,
                                    l_hvB,
                                    l_netUpdates[0],
                                    l_netUpdates[1] );
    }

    // update the cells' quantities
    l_hNew[l_ce]  -= i_scaling * l_netUpdates[0][0];
    l_hvNew[l_ce] -= i_scaling * l_netUpdates[0][1];

    l_hNew[l_ceT]  -= i_scaling * l_netUpdates[1][0];
    l_hvNew[l_ceT] -= i_scaling * l_netUpdates[1][1];
  }

    }
}

void tsunami_lab::patches::WavePropagation2d::setGhostOutflow() {
  t_real * l_h = m_h[m_step];
  t_real * l_hu = m_hu[m_step];
  t_real * l_hv = m_hv[m_step];

  // set left boundary
  l_hu[0] = l_hu[1];
  l_hv[0] = l_hv[1];

  // set right boundary
  l_hu[m_nCells+1] = l_hu[m_nCells];
  l_hv[m_nCells+1] = l_hu[m_nCells];

  // set all 4 boundaries for bathymetry and height

  // bottom boundary
  for ( t_idx l_ce = 0; l_ce < m_nCells + 1; l_ce++ ) {
    l_h[l_ce] = l_h[l_ce + getStride()];
    m_b[l_ce] = m_b[l_ce + getStride()];
  }

  // left boundary
  for ( t_idx l_ce = 0; l_ce < getStride() * getStride(); l_ce += getStride() ) {
    l_h[l_ce] = l_h[l_ce + 1];
    m_b[l_ce] = m_b[l_ce + 1];
  }

  // right boundary
  for ( t_idx l_ce = m_nCells + 1; l_ce < getStride() * getStride(); l_ce += getStride() ) {
    l_h[l_ce] = l_h[l_ce - 1];
    m_b[l_ce] = m_b[l_ce - 1];
  }

  // top boundary
  for ( t_idx l_ce = (getStride() - 1) * getStride(); l_ce < getStride() * getStride(); l_ce++) {
    l_h[l_ce] = l_h[l_ce - getStride()];
    m_b[l_ce] = m_b[l_ce - getStride()];
  }
}