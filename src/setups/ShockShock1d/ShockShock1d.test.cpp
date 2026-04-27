/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Tests the shock shock setup.
 **/
#include <catch2/catch.hpp>
#include "ShockShock1d.h"

TEST_CASE( "Test the one-dimensional shock shock setup.", "[ShockShock1d]" ) {

  tsunami_lab::setups::ShockShock1d l_shockShock( 10,
                                              20,
                                              3);

  // left side
  REQUIRE( l_shockShock.getHeight( 2, 0 ) == 10 );

  REQUIRE( l_shockShock.getMomentumX( 2, 0 ) == 20 );

  REQUIRE( l_shockShock.getMomentumY( 2, 0 ) == 0 );

  REQUIRE( l_shockShock.getBathymetry( 2, 0 ) == 0 );

  REQUIRE( l_shockShock.getHeight( 2, 5 ) == 10 );

  REQUIRE( l_shockShock.getMomentumX( 2, 5 ) == 20 );

  REQUIRE( l_shockShock.getMomentumY( 2, 5 ) == 0 );

  REQUIRE( l_shockShock.getBathymetry( 2, 5 ) == 0 );

  // right side
  REQUIRE( l_shockShock.getHeight( 4, 0 ) == 10 );

  REQUIRE( l_shockShock.getMomentumX( 4, 0 ) == -20 );

  REQUIRE( l_shockShock.getMomentumY( 4, 0 ) == 0 );

  REQUIRE( l_shockShock.getBathymetry( 4, 0 ) == 0 );

  REQUIRE( l_shockShock.getHeight( 4, 5 ) == 10 );

  REQUIRE( l_shockShock.getMomentumX( 4, 5 ) == -20 );

  REQUIRE( l_shockShock.getMomentumY( 4, 5 ) == 0 );  

  REQUIRE( l_shockShock.getBathymetry( 4, 5 ) == 0 );
}