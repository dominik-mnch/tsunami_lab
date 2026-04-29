/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Tests the subcritical setup.
 **/
#include <catch2/catch.hpp>
#include "Supercritical1d.h"

TEST_CASE( "Test the subcritical setup.", "[Subcritical1d]" ) {
  tsunami_lab::setups::Subcritical1d l_subcritical;

// testing within the hill
  REQUIRE( l_subcritical.getHeight( 10, 0 ) == Approx(-0.13) );

  REQUIRE( l_subcritical.getMomentumX( 10, 0 ) == Approx(0.18) );

  REQUIRE( l_subcritical.getMomentumY( 10, 0 ) == Approx(0) );

  REQUIRE( l_subcritical.getBathymetry( 10, 0 ) == Approx(-0.13) );

  REQUIRE( l_subcritical.getHeight( 10, 5 ) == Approx(0.13) );

  REQUIRE( l_subcritical.getMomentumX( 10, 5 ) == Approx(0.18) );

  REQUIRE( l_subcritical.getMomentumY( 10, 5 ) == Approx(0) );

  REQUIRE( l_subcritical.getBathymetry( 10, 5 ) == Approx(-0.13) );

//testing at the part of flat bathymetry
  REQUIRE( l_subcritical.getHeight( 2, 0 ) == Approx(0.33) );

  REQUIRE( l_subcritical.getMomentumX( 2, 0 ) == Approx(0.18) );

  REQUIRE( l_subcritical.getMomentumY( 2, 0 ) == Approx(0) );

  REQUIRE( l_subcritical.getBathymetry( 2, 0 ) == Approx(-0.33) );

  REQUIRE( l_subcritical.getHeight( 2, 5 ) == Approx(0.33) );

  REQUIRE( l_subcritical.getMomentumX( 2, 5 ) == Approx(-0.18) );

  REQUIRE( l_subcritical.getMomentumY( 2, 5 ) == Approx(0) );

  REQUIRE( l_subcritical.getBathymetry( 2, 5 ) == Approx(-0.33) );
}