/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Unit tests of the F-wave solver.
 **/
#include <catch2/catch.hpp>
#include "F_wave.h"

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

TEST_CASE( "Test the derivation of the F-wave net-updates with bathymetry.", "[FWaveUpdatesBathymetry]" ) {
  /*
   * Test case (bathymetry effect on net-updates):
   *
   *     left | right
   *   h:  10 | 10
   *  hu:   0 |  0
   *   b:  -2 | -4
   *
  * The Roe speeds and wave strengths are:
   *   s1 ≈ -9.902853,  s2 ≈ 9.902853
   *   alpha1 ≈ 9.902853, alpha2 ≈ -9.902853
   *
   * Net-updates from scaled eigenvectors:
   *
   *                        |  1 |   |  9.902853                    |   |  9.902853 |
   *   wave #1:  z1 = a1 *  |    | = |                              | = |           |
   *                        | s1 |   |  9.902853 * (-9.902853)      |   | -98.0665  |
   *
   *                        |  1 |   | -9.902853                    |   | -9.902853 |
   *   wave #2:  z2 = a2 *  |    | = |                              | = |           |
   *                        | s2 |   | -9.902853 * 9.902853         |   | -98.0665  |
   *
   *   s1 < 0 → wave #1 goes left:  netUpdateL = z1
   *   s2 > 0 → wave #2 goes right: netUpdateR = z2
   *
   *   wolframalpha.com query: 9.902853 * 9.902853
   *   result: ≈ 98.0665 (= 9.80665 * 10)
   */
  float l_netUpdatesL[2] = {0};
  float l_netUpdatesR[2] = {0};

  tsunami_lab::solvers::F_wave::netUpdates( 10,
                                            10,
                                            0,
                                            0,
                                            -2,
                                            -4,
                                            l_netUpdatesL,
                                            l_netUpdatesR );

  REQUIRE( l_netUpdatesL[0] == Approx(  9.902853 ) );
  REQUIRE( l_netUpdatesL[1] == Approx( -98.0665  ) );

  REQUIRE( l_netUpdatesR[0] == Approx( -9.902853 ) );
  REQUIRE( l_netUpdatesR[1] == Approx( -98.0665  ) );
}