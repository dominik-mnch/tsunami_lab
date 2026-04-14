/**
 * @author Magdalena Schwarzkopf 
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * f-wave-solver for the shallow water equations.
 **/
#include "F_wave.h"
#include <cmath>

void tsunami_lab::solvers::F_wave::eigenvalues( t_real   i_hL,
                                            t_real   i_hR,
                                            t_real   i_uL,
                                            t_real   i_uR,
                                            t_real & o_waveSpeedL,
                                            t_real & o_waveSpeedR ) {
  // pre-compute square-root ops
  t_real l_hSqrtL = std::sqrt( i_hL );
  t_real l_hSqrtR = std::sqrt( i_hR );

  // compute Roe averages
  t_real l_hRoe = 0.5f * ( i_hL + i_hR );
  t_real l_uRoe = l_hSqrtL * i_uL + l_hSqrtR * i_uR;
  l_uRoe /= l_hSqrtL + l_hSqrtR;

  // compute wave speeds
  t_real l_ghSqrtRoe = m_gSqrt * std::sqrt( l_hRoe );
  o_waveSpeedL = l_uRoe - l_ghSqrtRoe;
  o_waveSpeedR = l_uRoe + l_ghSqrtRoe;
}

void tsunami_lab::solvers::F_wave::waveStrengths( t_real   i_hL,
                                               t_real   i_hR,
                                               t_real   i_huL,
                                               t_real   i_huR,
                                               t_real   i_waveSpeedL,
                                               t_real   i_waveSpeedR,
                                               t_real & o_strengthL,
                                               t_real & o_strengthR ) {
  // compute inverse of right eigenvector-matrix
  t_real l_detInv = 1 / (i_waveSpeedR - i_waveSpeedL);

  t_real l_rInv[2][2] = {0};
  l_rInv[0][0] =  l_detInv * i_waveSpeedR;
  l_rInv[0][1] = -l_detInv;
  l_rInv[1][0] = -l_detInv * i_waveSpeedL;
  l_rInv[1][1] =  l_detInv;

  // compute jump in flux
  t_real l_deltaF[2] = {0};
  l_deltaF[0] = i_huR - i_huL;
  l_deltaF[1] = (std::pow(i_huR,2) / i_hR  + 0.5 * std::pow(m_gSqrt,2) * std::pow(i_hR,2))
                - (std::pow(i_huL,2) / i_hL  + 0.5 * std::pow(m_gSqrt,2) * std::pow(i_hL,2));

  // compute wave strengths
  o_strengthL  = l_rInv[0][0] * l_deltaF[0];
  o_strengthL += l_rInv[0][1] * l_deltaF[1];

  o_strengthR  = l_rInv[1][0] * l_deltaF[0];
  o_strengthR += l_rInv[1][1] * l_deltaF[1];
}

void tsunami_lab::solvers::F_wave::netUpdates( t_real i_hL,
                                t_real i_hR,
                                t_real i_huL,
                                t_real i_huR,
                                t_real o_netUpdateL[2],
                                t_real o_netUpdateR[2] ) {

 // compute particle velocities
  t_real l_uL = i_huL / i_hL;
  t_real l_uR = i_huR / i_hR;

  // compute wave speeds
  t_real l_sL = 0;
  t_real l_sR = 0;

  eigenvalues( i_hL,
              i_hR,
              l_uL,
              l_uR,
              l_sL,
              l_sR );

  // compute wave strengths
  t_real l_aL = 0;
  t_real l_aR = 0;

  waveStrengths( i_hL,
                 i_hR,
                 i_huL,
                 i_huR,
                 l_sL,
                 l_sR,
                 l_aL,
                 l_aR );

   //initialize waves
   t_real l_Z1[2] = {0}; 
   t_real l_Z2[2] = {0};

   //calculate first wave
   l_Z1[0] = l_aL;
   l_Z1[1] = l_aL * l_sL;

   //calculate second wave
   l_Z2[0] = l_aR;
   l_Z2[1] = l_aR * l_sR;

   // set net-updates depending on wave speeds
  for( unsigned short l_qt = 0; l_qt < 2; l_qt++ ) {
    // init
    o_netUpdateL[l_qt] = 0;
    o_netUpdateR[l_qt] = 0;

    // 1st wave
    if( l_sL < 0 ) {
      o_netUpdateL[l_qt] += l_Z1[l_qt];
    }
    else if (l_sL > 0) {
      o_netUpdateR[l_qt] += l_Z1[l_qt];
    }

    // 2nd wave
    if( l_sR > 0 ) {
      o_netUpdateR[l_qt] += l_Z2[l_qt];
    }
    else if (l_sR < 0) {
      o_netUpdateL[l_qt] += l_Z2[l_qt];
    }
  }

}