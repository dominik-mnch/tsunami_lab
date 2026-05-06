/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Tests the two-dimensional circular dam break setup.
 **/
#include <catch2/catch.hpp>
#include "CircularDamBreak2d.h"

TEST_CASE( "Test the two-dimensional circular dam break setup.", "[CircularDamBreak2dNoVelocity]" ) {
  tsunami_lab::setups::CircularDamBreak2d l_damBreak( 10,
                                                      5,
                                                      3,
                                                      4,
                                                      1,
                                                      2,
                                                      50,
                                                      0,
                                                      10 );
    // inside the dam
  REQUIRE( l_damBreak.getHeight( 50, 0 ) == 10 );

  REQUIRE( l_damBreak.getMomentumX( 50, 0 ) == 3 );

  REQUIRE( l_damBreak.getMomentumY( 50, 0 ) == 4 );

  REQUIRE( l_damBreak.getBathymetry ( 50, 0 ) == 0 );

  REQUIRE( l_damBreak.getHeight( 2, 3 ) == 5 );


  // outside the dam

  REQUIRE( l_damBreak.getHeight( 61, 0 ) == 5 );

  REQUIRE( l_damBreak.getMomentumX( 61, 0 ) == 1 );

  REQUIRE( l_damBreak.getMomentumY( 61, 0 ) == 2 );

  REQUIRE( l_damBreak.getBathymetry ( 61, 0 ) == 0 );

  // outside the dam, but in the part of the hill
  REQUIRE( l_damBreak.getHeight( -22, 1 ) == 5 );

  REQUIRE( l_damBreak.getMomentumX( -22, 1 ) == 1 );

  REQUIRE( l_damBreak.getMomentumY( -22, 1 ) == 2 );

  REQUIRE( l_damBreak.getBathymetry ( -22, 1 ) == Approx(1.2125) );
  }

