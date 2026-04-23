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

  tsunami_lab::t_real i_b[6] = {-5, -5, -5, -5, -5, -5};

  tsunami_lab::setups::RareRare1d l_rareRare( 10,
                                              20,
                                               3,
                                              i_b);

  // left side
  REQUIRE( l_rareRare.getHeight( 2, 0 ) == 10 );

  REQUIRE( l_rareRare.getMomentumX( 2, 0 ) == -20 );

  REQUIRE( l_rareRare.getMomentumY( 2, 0 ) == 0 );

  REQUIRE( l_rareRare.getBathymetry( 2, 0 ) == -5 );

  REQUIRE( l_rareRare.getHeight( 2, 5 ) == 10 );

  REQUIRE( l_rareRare.getMomentumX( 2, 5 ) == -20 );

  REQUIRE( l_rareRare.getMomentumY( 2, 5 ) == 0 );

  REQUIRE( l_rareRare.getBathymetry( 2, 5 ) == -5 );

  // right side
  REQUIRE( l_rareRare.getHeight( 4, 0 ) == 10 );

  REQUIRE( l_rareRare.getMomentumX( 4, 0 ) == 20 );

  REQUIRE( l_rareRare.getMomentumY( 4, 0 ) == 0 );

  REQUIRE( l_rareRare.getBathymetry( 4, 0 ) == -5 );

  REQUIRE( l_rareRare.getHeight( 4, 5 ) == 10 );

  REQUIRE( l_rareRare.getMomentumX( 4, 5 ) == 20 );

  REQUIRE( l_rareRare.getMomentumY( 4, 5 ) == 0 );
  
  REQUIRE( l_rareRare.getBathymetry( 4, 5 ) == -5 );
}