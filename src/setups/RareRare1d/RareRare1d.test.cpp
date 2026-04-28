/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Tests the rare rare setup.
 **/
#include <catch2/catch.hpp>
#include "RareRare1d.h"

TEST_CASE( "Test the one-dimensional rare rare setup.", "[RareRare1d]" ) {

  tsunami_lab::setups::RareRare1d l_rareRare( 10,
                                              20,
                                               3);

  // left side
  REQUIRE( l_rareRare.getHeight( 2, 0 ) == 10 );

  REQUIRE( l_rareRare.getMomentumX( 2, 0 ) == -20 );

  REQUIRE( l_rareRare.getMomentumY( 2, 0 ) == 0 );

  REQUIRE( l_rareRare.getBathymetry( 2, 0 ) == 0 );

  REQUIRE( l_rareRare.getHeight( 2, 5 ) == 10 );

  REQUIRE( l_rareRare.getMomentumX( 2, 5 ) == -20 );

  REQUIRE( l_rareRare.getMomentumY( 2, 5 ) == 0 );

  REQUIRE( l_rareRare.getBathymetry( 2, 5 ) == 0 );

  // right side
  REQUIRE( l_rareRare.getHeight( 4, 0 ) == 10 );

  REQUIRE( l_rareRare.getMomentumX( 4, 0 ) == 20 );

  REQUIRE( l_rareRare.getMomentumY( 4, 0 ) == 0 );

  REQUIRE( l_rareRare.getBathymetry( 4, 0 ) == 0 );

  REQUIRE( l_rareRare.getHeight( 4, 5 ) == 10 );

  REQUIRE( l_rareRare.getMomentumX( 4, 5 ) == 20 );

  REQUIRE( l_rareRare.getMomentumY( 4, 5 ) == 0 );
  
  REQUIRE( l_rareRare.getBathymetry( 4, 5 ) == 0 );
}

TEST_CASE( "Test parabolic bathymetry for one-dimensional rare rare setup.", "[RareRare1d]" ) {
  // b(x) = 1 - 0.25 * (x - 4)^2 for |x - 4| <= sqrt(4) = 2, else b(x) = 0
  tsunami_lab::setups::RareRare1d l_rareRare( 10,
                                              20,
                                              3,
                                              true,
                                              1,
                                              4,
                                              -0.25 );

  REQUIRE( l_rareRare.getBathymetry( 4, 0 ) == Approx( 1 ) );
  REQUIRE( l_rareRare.getBathymetry( 3.5, 0 ) == Approx( 0.9375 ) );
  REQUIRE( l_rareRare.getBathymetry( 2, 0 ) == Approx( 0 ) );
  REQUIRE( l_rareRare.getBathymetry( 6, 0 ) == Approx( 0 ) );
  REQUIRE( l_rareRare.getBathymetry( 7, 0 ) == Approx( 0 ) );
}