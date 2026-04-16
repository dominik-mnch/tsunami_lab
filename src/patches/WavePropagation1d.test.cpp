/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Unit tests for the one-dimensional wave propagation patch.
 **/
#include <catch2/catch.hpp>
#include "WavePropagation1d.h"

TEST_CASE( "Test the 1d wave propagation with Roe solver.", "[WavePropRoe1d]" ) {
  /*
   * Test case:
   *
   *   Single dam break problem between cell 49 and 50.
   *     left | right
   *       10 | 8
   *        0 | 0
   *
   *   Elsewhere steady state.
   *
   * The net-updates at the respective edge are given as
   * (see derivation in Roe solver):
   *    left          | right
   *      9.394671362 | -9.394671362
   *    -88.25985     | -88.25985
   */

  // construct solver and setup a dambreak problem
  tsunami_lab::patches::WavePropagation1d m_waveProp( 100 , false );

  for( std::size_t l_ce = 0; l_ce < 50; l_ce++ ) {
    m_waveProp.setHeight( l_ce,
                          0,
                          10 );
    m_waveProp.setMomentumX( l_ce,
                             0,
                             0 );
  }
  for( std::size_t l_ce = 50; l_ce < 100; l_ce++ ) {
    m_waveProp.setHeight( l_ce,
                          0,
                          8 );
    m_waveProp.setMomentumX( l_ce,
                             0,
                             0 );
  }

  // set outflow boundary condition
  m_waveProp.setGhostOutflow();

  // perform a time step
  m_waveProp.timeStep( 0.1 );

  // steady state
  for( std::size_t l_ce = 0; l_ce < 49; l_ce++ ) {
    REQUIRE( m_waveProp.getHeight()[l_ce]   == Approx(10) );
    REQUIRE( m_waveProp.getMomentumX()[l_ce] == Approx(0) );
  }

  // dam-break
  REQUIRE( m_waveProp.getHeight()[49]   == Approx(10 - 0.1 * 9.394671362) );
  REQUIRE( m_waveProp.getMomentumX()[49] == Approx( 0 + 0.1 * 88.25985) );

  REQUIRE( m_waveProp.getHeight()[50]   == Approx(8 + 0.1 * 9.394671362) );
  REQUIRE( m_waveProp.getMomentumX()[50] == Approx(0 + 0.1 * 88.25985) );

  // steady state
  for( std::size_t l_ce = 51; l_ce < 100; l_ce++ ) {
    REQUIRE( m_waveProp.getHeight()[l_ce]   == Approx(8) );
    REQUIRE( m_waveProp.getMomentumX()[l_ce] == Approx(0) );
  }
}

TEST_CASE( "Test the 1d wave propagation with F-Wave solver.", "[WaveProp1d]" ) {
  /*
   *
   *
   * 
   * Test case (1st test case from csv file):
   *  h:   8899.326826472694  | 8899.326826472694
   *  hu:   122.0337839252433 | -122.0337839252433
   *
   *   Elsewhere steady state.
   *
   * The net-updates at the respective edge are given as
   * (see derivation in Roe solver):
   *    left          | right
   *    -45.6276      | -6.97236
   *    298.14807     | -71.20777
   */

  // construct solver and setup a dambreak problem
  tsunami_lab::patches::WavePropagation1d m_waveProp( 100 , true );

  for( std::size_t l_ce = 0; l_ce < 50; l_ce++ ) {
    m_waveProp.setHeight( l_ce,
                          0,
                          8899.326826472694 );
    m_waveProp.setMomentumX( l_ce,
                             0,
                             122.0337839252433 );
  }
  for( std::size_t l_ce = 50; l_ce < 100; l_ce++ ) {
    m_waveProp.setHeight( l_ce,
                          0,
                          8899.326826472694 );
    m_waveProp.setMomentumX( l_ce,
                             0,
                             -122.0337839252433 );
  }

  // derive maximum wave speed in setup; the momentum is ignored
  tsunami_lab::t_real m_speedMax = std::sqrt( 9.81 * 8899.326826472694 );

  // derive constant time step; changes at simulation time are ignored
  tsunami_lab::t_real m_scaling = 0.5 / m_speedMax; 


  for (int i = 0; i < 100; i++){
      // set outflow boundary condition
    m_waveProp.setGhostOutflow();

    // perform a time step
    m_waveProp.timeStep( m_scaling );
  }


  REQUIRE( m_waveProp.getHeight()[49]   == Approx(8899.739847378269) );
}