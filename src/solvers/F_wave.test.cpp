/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Unit tests of the F-wave solver.
 **/
#include <catch2/catch.hpp>
#define private public
#include "F_wave.h"
#undef public

TEST_CASE( "Test the derivation of the eigenvalues.", "[Eigenvalues]" ) {
   /*
    * Test case:
    *  h: 4 | 10.3
    *  u: 8 | -2
    *
    * roe height: 7.15
    * roe velocity: (sqrt(4) * 8 + (-2) * sqrt(10.3)) / ( sqrt(4) + sqrt(10.3) )
    *               = 1.8392422450039443
    * roe speeds: s1 = 1.8392422450039443 - sqrt(9.80665 * 7.15) = -6.5343798805
    *             s2 = 1.8392422450039443 + sqrt(9.80665 * 7.15) = 10.2128643705
    */
  float l_waveSpeedL = 0;
  float l_waveSpeedR = 0;
  tsunami_lab::solvers::F_wave::eigenvalues( 4,
                                         10.3,
                                         8,
                                         -2,
                                         l_waveSpeedL,
                                         l_waveSpeedR );

  REQUIRE( l_waveSpeedL == Approx(  -6.5343798805 ) );
  REQUIRE( l_waveSpeedR == Approx(  10.2128643705 ) );
}

TEST_CASE( "Test the derivation of the Roe wave strengths.", "[RoeStrengths]" ) {
  /*
   * Test case:
   *  h:   4 | 10.3
   *  u:   8 | -2
   *  hu: 32 |-20.6
   *
   * The derivation of the Roe speeds (s1, s2) is given above.
   *
   *  Matrix of right eigenvectors:
   *
   *      | 1   1 |
   *  R = |       |
   *      | s1 s2 |
   *
   * Inversion yields:
   *
   * wolframalpha.com query: invert {{1, 1}, {-6.5343798805, 10.2128643705}}
   *
   *         |  0.60982357559  -0.059711316382 |
   * R_inv = |                                 |
   *         |  0.39017642441   0.059711316382 |
   *
   * Multiplicaton with the jump in flux gives the wave strengths:
   *
   * wolframalpha.com query: {{0.60982357559, -0.059711316382},
   *                          {0.39017642441, 0.059711316382}} 
   *                        *  {-20.6 - 32,
                                (-20.6)^2 / 10.3  + 0.5 * 3.131557121^2 * 10.3^2
                                - (32^2 /  4    + 0.5 * 3.131557121^2 * 4^2 ) }
   *
   *        |                                               -20.6 - 32                                       |   | -45.6276 |
   * Rinv * |                                                                                                | = |          |
   *        | (-20.6)^2 / 10.3 + 0.5 * 3.131557121^2 * 10.3^2 - (32^2) /  4    + 0.5 * 3.131557121^2 * 4^2 ) |   | -6.97236 |
   */
  float l_strengthL = 0;
  float l_strengthR = 0;

  tsunami_lab::solvers::F_wave::waveStrengths( 4,
                                            10.3,
                                            32,
                                            -20.6,
                                            0,
                                            0,
                                            -6.5343798805,
                                            10.2128643705,
                                            l_strengthL,
                                            l_strengthR );

  REQUIRE( l_strengthL == Approx( -45.6276 ) );
  REQUIRE( l_strengthR == Approx( -6.97236 ) );
}

TEST_CASE( "Test the derivation of the F wave net-updates.", "[FWaveUpdates]" ) {
  /*
   * Test case:
   *  h:   4 | 10.3
   *  u:   8 | -2
   *  hu: 32 |-20.6
   *
   * The derivation of the Roe speeds (s1, s2) and wave strengths (a1, a1) is given above.
   *
   * The net-updates are given through the scaled eigenvectors.
   * 
   *
   *                      |  1 |   |  -45.6276  |
   * wave #1:   z1 = a1 * |    | = |            |
   *                      | s1 |   |  298.14807 |
   *
   *                      |  1 |   | -6.97236  |
   * wave #2:   z2 = a2 * |    | = |           |
   *                      | s2 |   | -71.20777 |
   * 
   * update #1: z1 (because s1 is negative)
   * 
   * update #2: z2 (because s2 is positive)
   * 
   */
  float l_netUpdatesL[2] = {0};
  float l_netUpdatesR[2] = {0};

  tsunami_lab::solvers::F_wave::netUpdates( 4,
                                         10.3,
                                         32,
                                         -20.6,
                                         0,
                                         0,
                                         l_netUpdatesL,
                                         l_netUpdatesR );

  REQUIRE( l_netUpdatesL[0] == Approx( -45.6276 ) );
  REQUIRE( l_netUpdatesL[1] == Approx( 298.14807 ) );

  REQUIRE( l_netUpdatesR[0] == Approx( -6.97236 ) );
  REQUIRE( l_netUpdatesR[1] == Approx( -71.20777 ) );

  /*
   * Test case (supersonic problem):
   *
   *      left | right
   *   h:  10  | 10
   *   hu: 80  | 120
   *
   * Roe speeds are given as:
   *
   *   s1 =  0.097
   *   s2 = 19.903
   *
   * Inversion of the matrix of right Eigenvectors:
   * 
   *   wolframalpha.com query: invert {{1, 1}, {0.097, 19.903}}
   *
   *          |  1.0049     -0.0504898  |
   *   Rinv = |                         |
   *          | -0.00489751   0.0504898 |
   * 
   * 
   *
   * Multiplicaton with the jump in flux gives the wave strengths:
   *
   *        |                                     120 - 80                                  |   | -0.19584  |   | a1 |
   * Rinv * |                                                                               | = |           | = |    |
   *        |  (120)^2 / 10 + 0.5 * 9.80665 * 10^2 - ((80^2) / 10 + 0.5 * 9.80665 * 10^2 )) |   |  40.1959  |   | a2 |
   *
   * The net-updates are given through the scaled eigenvectors.
   *
   *                      |  1 |   | -0.1959  |
   * wave #1: z1 = a1 *   |    | = |          |
   *                      | s1 |   |  -0.0190 |
   *
   *                      |  1 |   | 40.1959  |
   * wave #2: z2 = a2 *   |    | = |          |
   *                      | s2 |   | 800.0190 |
   * 
   * update #1 = 0 (because both s-values are positive)
   * 
   *                        | 40  |
   * update #2 = z1 + z2 =  |     |
   *                        | 800 |
   */

  tsunami_lab::solvers::F_wave::netUpdates( 10,
                                         10,
                                         80,
                                         120,
                                         0,
                                         0,
                                         l_netUpdatesL,
                                         l_netUpdatesR ); 

  REQUIRE( l_netUpdatesL[0] ==  Approx(0) );
  REQUIRE( l_netUpdatesL[1] ==  Approx(0) );

  REQUIRE( l_netUpdatesR[0] ==  Approx(40) );
  REQUIRE( l_netUpdatesR[1] ==  Approx(800) );

  /*
   * Test case (trivial steady state):
   *
   *     left | right
   *   h:  10 | 10
   *  hu:   0 |  0
   */

  tsunami_lab::solvers::F_wave::netUpdates( 10,
                                         10,
                                         0,
                                         0,
                                         0,
                                         0,
                                         l_netUpdatesL,
                                         l_netUpdatesR );

  REQUIRE( l_netUpdatesL[0] == Approx(0) );
  REQUIRE( l_netUpdatesL[1] == Approx(0) );

  REQUIRE( l_netUpdatesR[0] == Approx(0) );
  REQUIRE( l_netUpdatesR[1] == Approx(0) );
}