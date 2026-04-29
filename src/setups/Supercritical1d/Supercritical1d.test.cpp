/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Tests the supercritical setup.
 **/
#include <catch2/catch.hpp>
#include "Supercritical1d.h"

TEST_CASE( "Test the supercritical setup.", "[Supercritical1d]" ) {
  tsunami_lab::setups::Supercritical1d l_supercritical;

// testing within the hill
  REQUIRE( l_supercritical.getHeight( 10, 0 ) == Approx(0.13) );

  REQUIRE( l_supercritical.getMomentumX( 10, 0 ) == Approx(0.18) );

  REQUIRE( l_supercritical.getMomentumY( 10, 0 ) == Approx(0) );

  REQUIRE( l_supercritical.getBathymetry( 10, 0 ) == Approx(-0.13) );

  REQUIRE( l_supercritical.getHeight( 10, 5 ) == Approx(0.13) );

  REQUIRE( l_supercritical.getMomentumX( 10, 5 ) == Approx(0.18) );

  REQUIRE( l_supercritical.getMomentumY( 10, 5 ) == Approx(0) );

  REQUIRE( l_supercritical.getBathymetry( 10, 5 ) == Approx(-0.13) );

//testing at the part of flat bathymetry
  REQUIRE( l_supercritical.getHeight( 2, 0 ) == Approx(0.33) );

  REQUIRE( l_supercritical.getMomentumX( 2, 0 ) == Approx(0.18) );

  REQUIRE( l_supercritical.getMomentumY( 2, 0 ) == Approx(0) );

  REQUIRE( l_supercritical.getBathymetry( 2, 0 ) == Approx(-0.33) );

  REQUIRE( l_supercritical.getHeight( 2, 5 ) == Approx(0.33) );

  REQUIRE( l_supercritical.getMomentumX( 2, 5 ) == Approx(0.18) );

  REQUIRE( l_supercritical.getMomentumY( 2, 5 ) == Approx(0) );

  REQUIRE( l_supercritical.getBathymetry( 2, 5 ) == Approx(-0.33) );
}