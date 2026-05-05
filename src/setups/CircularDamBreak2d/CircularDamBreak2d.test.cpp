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
                                                      0,
                                                      0,
                                                      0,
                                                      0,
                                                      0,
                                                      0,
                                                     10);
    // inside the dam
  REQUIRE( l_damBreak.getHeight( 2, 3 ) == 10 );

  REQUIRE( l_damBreak.getMomentumX( 2, 3 ) == 0 );

  REQUIRE( l_damBreak.getMomentumY( 2, 3 ) == 0 );

  REQUIRE( l_damBreak.getBathymetry ( 2, 3 ) == 0 );


  // outside the dam

  REQUIRE( l_damBreak.getHeight( 10, 15 ) == 5 );

  REQUIRE( l_damBreak.getMomentumX( 10, 15 ) == 0 );

  REQUIRE( l_damBreak.getMomentumY( 10, 15 ) == 0 );

  REQUIRE( l_damBreak.getBathymetry ( 10, 15 ) == 0 );

  // outside the dam, but in the part of the hill
  REQUIRE( l_damBreak.getHeight( -22, 1 ) == 5 );

  REQUIRE( l_damBreak.getMomentumX( -22, 1 ) == 0 );

  REQUIRE( l_damBreak.getMomentumY( -22, 1 ) == 0 );

  REQUIRE( l_damBreak.getBathymetry ( -22, 1 ) == Approx(1.2125) );
  }

