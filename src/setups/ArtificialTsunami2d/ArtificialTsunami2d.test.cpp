/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Tests the artificial tsunami 2d setup.
 **/
#include <catch2/catch.hpp>
#include "ArtificialTsunami2d.h"

TEST_CASE("ArtificialTsunami2d: constants are correct", "[setup]") {
    tsunami_lab::setups::ArtificialTsunami2d setup;

    REQUIRE(setup.getHeight(0.0, 0.0) == Approx(100.0));
    REQUIRE(setup.getHeight(123.0, -45.0) == Approx(100.0));

    REQUIRE(setup.getBathymetry(0.0, 0.0) == Approx(-100.0));
    REQUIRE(setup.getBathymetry(999.0, 999.0) == Approx(-100.0));

    REQUIRE(setup.getMomentumX(0.0, 0.0) == Approx(0.0));
    REQUIRE(setup.getMomentumY(0.0, 0.0) == Approx(0.0));
}

TEST_CASE("ArtificialTsunami2d: displacement is zero outside region", "[setup]") {
    tsunami_lab::setups::ArtificialTsunami2d setup;

    REQUIRE(setup.getDisplacement(600.0, 0.0) == Approx(0.0));
    REQUIRE(setup.getDisplacement(-600.0, 0.0) == Approx(0.0));
    REQUIRE(setup.getDisplacement(0.0, 600.0) == Approx(0.0));
    REQUIRE(setup.getDisplacement(0.0, -600.0) == Approx(0.0));
}

TEST_CASE("ArtificialTsunami2d: displacement known value check", "[setup]") {
    tsunami_lab::setups::ArtificialTsunami2d setup;

    REQUIRE(setup.getDisplacement(250.0, 250.0) == Approx(-3.75).margin(1e-12));
}
