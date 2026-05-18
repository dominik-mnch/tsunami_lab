/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Tests the 2D tsunami event setup with netCDF input.
 **/
#include <catch2/catch.hpp>
#include "TsunamiEvent2d.h"

namespace {
	using tsunami_lab::t_real;
}

TEST_CASE( "Test TsunamiEvent2d setup reads netCDF files correctly.", "[TsunamiEvent2d]" ) {
	tsunami_lab::setups::TsunamiEvent2d l_setup( ARTIFICIAL_TSUNAMI_BATHY_NC,
	                                             ARTIFICIAL_TSUNAMI_DISP_NC );

	SECTION( "Water height is 100 everywhere (bathymetry uniformly -100)" ) {
		// b_in = -100 everywhere; getHeight = max(-b_in, delta) = max(100, 20) = 100
		REQUIRE( l_setup.getHeight(    0.0,    0.0 ) == Approx( 100.0 ) );
		REQUIRE( l_setup.getHeight( 1000.0, 1000.0 ) == Approx( 100.0 ) );
		REQUIRE( l_setup.getHeight(-4000.0, 3000.0 ) == Approx( 100.0 ) );
	}

	SECTION( "Momenta are zero everywhere" ) {
		REQUIRE( l_setup.getMomentumX( 0.0, 0.0 ) == Approx( 0.0 ) );
		REQUIRE( l_setup.getMomentumY( 0.0, 0.0 ) == Approx( 0.0 ) );
	}

	SECTION( "Bathymetry equals base depth plus displacement" ) {
		// b_in = -100 everywhere, so getBathymetry = min(-100, -20) + d = -100 + d
		// Displacement is bounded by amplitude 5, so bathymetry lies in [-105, -95]
		REQUIRE( l_setup.getBathymetry(    0.0,    0.0 ) >= -105.0 );
		REQUIRE( l_setup.getBathymetry(    0.0,    0.0 ) <= -95.0  );
		REQUIRE( l_setup.getBathymetry( -495.0, -495.0 ) >= -105.0 );
		REQUIRE( l_setup.getBathymetry( -495.0, -495.0 ) <= -95.0  );
		REQUIRE( l_setup.getBathymetry( 2000.0, 3000.0 ) >= -105.0 );
		REQUIRE( l_setup.getBathymetry( 2000.0, 3000.0 ) <= -95.0  );
	}
}

