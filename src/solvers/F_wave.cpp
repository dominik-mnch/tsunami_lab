/**
 * @author Magdalena Schwarzkopf 
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * f-wave-solver for the shallow water equations.
 **/
#include "F_wave.h"
#include <cmath>

void tsunami_lab::solvers::F_wave::netUpdates( t_real i_hL,
                                t_real i_hR,
                                t_real i_huL,
                                t_real i_huR,
                                t_real i_bL,
                                t_real i_bR,
                                t_real o_netUpdateL[2],
                                t_real o_netUpdateR[2] ) {
  t_real l_g = m_gSqrt * m_gSqrt;
  t_real l_gHalf = 0.5f * l_g;

 // compute particle velocities
  t_real l_uL = i_huL / i_hL;
  t_real l_uR = i_huR / i_hR;

  // compute wave speeds
  t_real l_hSqrtL = std::sqrt( i_hL );
  t_real l_hSqrtR = std::sqrt( i_hR );
  t_real l_uRoe = (l_hSqrtL * l_uL + l_hSqrtR * l_uR) / (l_hSqrtL + l_hSqrtR);
  t_real l_ghSqrtRoe = std::sqrt( l_gHalf * (i_hL + i_hR) );
  t_real l_sL = l_uRoe - l_ghSqrtRoe;
  t_real l_sR = l_uRoe + l_ghSqrtRoe;

  // compute jump in flux and bathymetry source term
  t_real l_deltaF0 = i_huR - i_huL;
  t_real l_deltaF1 = (i_huR * l_uR + l_gHalf * i_hR * i_hR)
                   - (i_huL * l_uL + l_gHalf * i_hL * i_hL);
  t_real l_bathymetryTerm = -l_gHalf * (i_bR - i_bL) * (i_hL + i_hR);
  l_deltaF1 = l_deltaF1 - l_bathymetryTerm;

  // compute wave strengths
  t_real l_detInv = 1 / (l_sR - l_sL);
  t_real l_aL = (l_sR * l_deltaF0 - l_deltaF1) * l_detInv;
  t_real l_aR = (l_deltaF1 - l_sL * l_deltaF0) * l_detInv;

  t_real l_zL0 = l_aL;
  t_real l_zL1 = l_aL * l_sL;
  t_real l_zR0 = l_aR;
  t_real l_zR1 = l_aR * l_sR;

  o_netUpdateL[0] = 0;
  o_netUpdateL[1] = 0;
  o_netUpdateR[0] = 0;
  o_netUpdateR[1] = 0;

  // set net-updates depending on wave speeds
  if( l_sL < 0 ) {
    o_netUpdateL[0] += l_zL0;
    o_netUpdateL[1] += l_zL1;
  }
  else if( l_sL > 0 ) {
    o_netUpdateR[0] += l_zL0;
    o_netUpdateR[1] += l_zL1;
  }

  if( l_sR > 0 ) {
    o_netUpdateR[0] += l_zR0;
    o_netUpdateR[1] += l_zR1;
  }
  else if( l_sR < 0 ) {
    o_netUpdateL[0] += l_zR0;
    o_netUpdateL[1] += l_zR1;
  }

}