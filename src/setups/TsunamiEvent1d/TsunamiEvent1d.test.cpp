/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Tests the tsunami event setup.
 **/
#include <catch2/catch.hpp>
#include "TsunamiEvent1d.h"

TEST_CASE( "Test the tsunami event setup.", "[TsunamiEvent1d]" ) {
  tsunami_lab::setups::TsunamiEvent1d l_tsunami_event;

  // case b_in < 0, delta > -b_in
  REQUIRE( l_tsunami_event.getHeight( 231.656660144, 0 ) == Approx(20) );

  REQUIRE( l_tsunami_event.getMomentumX( 231.656660144, 0 ) == Approx(0) );

  REQUIRE( l_tsunami_event.getMomentumY( 231.656660144, 0 ) == Approx(0) );

  REQUIRE( l_tsunami_event.getBathymetry( 231.656660144, 0 ) == Approx(-20) );

  // case b_in < 0, delta  < -b_in
  REQUIRE( l_tsunami_event.getHeight( 4169.8197377, 0 ) == Approx(22.2054342573) );

  REQUIRE( l_tsunami_event.getMomentumX( 4169.8197377, 0 ) == Approx(0) );

  REQUIRE( l_tsunami_event.getMomentumY( 4169.8197377, 0 ) == Approx(0) );

  REQUIRE( l_tsunami_event.getBathymetry( 4169.8197377, 0 ) == Approx(-22.2054342573) );

    // case b_in >= 0
  REQUIRE( l_tsunami_event.getHeight( 0, 0 ) == Approx(0) );

  REQUIRE( l_tsunami_event.getMomentumX( 0, 0 ) == Approx(0) );

  REQUIRE( l_tsunami_event.getMomentumY( 0, 0 ) == Approx(0) );

  REQUIRE( l_tsunami_event.getBathymetry( 0, 0 ) == Approx(20) );


}
