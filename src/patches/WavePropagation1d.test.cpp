/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Unit tests for the one-dimensional wave propagation patch.
 **/
#include <catch2/catch.hpp>
#include "WavePropagation1d.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#ifndef MIDDLE_STATES_CSV
#define MIDDLE_STATES_CSV "middle_states.csv"
#endif

TEST_CASE( "Test the 1d wave propagation with Roe solver.", "[WavePropRoe1d]" ) {
  /*
   * Test case:
   *
   *   Single dam break problem between cell 49 and 50.
   *     left | right
   *       10 | 8
   *        0 | 0
   *
   *   Elsewhere steady state.
   *
   * The net-updates at the respective edge are given as
   * (see derivation in Roe solver):
   *    left          | right
   *      9.394671362 | -9.394671362
   *    -88.25985     | -88.25985
   */

  // construct solver and setup a dambreak problem
  tsunami_lab::patches::WavePropagation1d m_waveProp( 100 , false );

  for( std::size_t l_ce = 0; l_ce < 50; l_ce++ ) {
    m_waveProp.setHeight( l_ce,
                          0,
                          10 );
    m_waveProp.setMomentumX( l_ce,
                             0,
                             0 );
  }
  for( std::size_t l_ce = 50; l_ce < 100; l_ce++ ) {
    m_waveProp.setHeight( l_ce,
                          0,
                          8 );
    m_waveProp.setMomentumX( l_ce,
                             0,
                             0 );
  }

  // set outflow boundary condition
  m_waveProp.setGhostOutflow();

  // perform a time step
  m_waveProp.timeStep( 0.1 );

  // steady state
  for( std::size_t l_ce = 0; l_ce < 49; l_ce++ ) {
    REQUIRE( m_waveProp.getHeight()[l_ce]   == Approx(10) );
    REQUIRE( m_waveProp.getMomentumX()[l_ce] == Approx(0) );
  }

  // dam-break
  REQUIRE( m_waveProp.getHeight()[49]   == Approx(10 - 0.1 * 9.394671362) );
  REQUIRE( m_waveProp.getMomentumX()[49] == Approx( 0 + 0.1 * 88.25985) );

  REQUIRE( m_waveProp.getHeight()[50]   == Approx(8 + 0.1 * 9.394671362) );
  REQUIRE( m_waveProp.getMomentumX()[50] == Approx(0 + 0.1 * 88.25985) );

  // steady state
  for( std::size_t l_ce = 51; l_ce < 100; l_ce++ ) {
    REQUIRE( m_waveProp.getHeight()[l_ce]   == Approx(8) );
    REQUIRE( m_waveProp.getMomentumX()[l_ce] == Approx(0) );
  }
}

TEST_CASE( "Test the 1d wave propagation with F-wave solver and bathymetry step.", "[WavePropFWaveBathymetry1d]" ) {
  /*
   * Test case:
   *
   *   Single dam break problem with bathymetry between cell 49 and 50.
   *     left | right
   *   h:  10 | 10
   *  hu:   0 |  0
   *   b:  -2 | -4
   *
   *   Elsewhere steady state.
   *
   * At the interface, the F-wave solver yields:
   *   netUpdateL = (  9.902853, -98.0665 )
   *   netUpdateR = ( -9.902853, -98.0665 )
   *
   * With dt/dx = 0.1:
   *   h49  = 10 - 0.1 *  9.902853 =  9.0097147
   *   hu49 =  0 - 0.1 * -98.0665  =  9.80665
   *   h50  = 10 - 0.1 * -9.902853 = 10.9902853
   *   hu50 =  0 - 0.1 * -98.0665  =  9.80665
   */

  tsunami_lab::patches::WavePropagation1d l_waveProp( 100, true );

  for( std::size_t l_ce = 0; l_ce < 50; l_ce++ ) {
    l_waveProp.setHeight( l_ce,
                          0,
                          10 );
    l_waveProp.setMomentumX( l_ce,
                             0,
                             0 );
    l_waveProp.setBathymetry( l_ce,
                              0,
                              -2 );
  }

  for( std::size_t l_ce = 50; l_ce < 100; l_ce++ ) {
    l_waveProp.setHeight( l_ce,
                          0,
                          10 );
    l_waveProp.setMomentumX( l_ce,
                             0,
                             0 );
    l_waveProp.setBathymetry( l_ce,
                              0,
                              -4 );
  }

  l_waveProp.setGhostOutflow();
  l_waveProp.timeStep( 0.1 );

  for( std::size_t l_ce = 0; l_ce < 49; l_ce++ ) {
    REQUIRE( l_waveProp.getHeight()[l_ce] == Approx(10) );
    REQUIRE( l_waveProp.getMomentumX()[l_ce] == Approx(0) );
  }

  REQUIRE( l_waveProp.getHeight()[49] == Approx( 9.0097147 ) );
  REQUIRE( l_waveProp.getMomentumX()[49] == Approx( 9.80665 ) );

  REQUIRE( l_waveProp.getHeight()[50] == Approx( 10.9902853 ) );
  REQUIRE( l_waveProp.getMomentumX()[50] == Approx( 9.80665 ) );

  for( std::size_t l_ce = 51; l_ce < 100; l_ce++ ) {
    REQUIRE( l_waveProp.getHeight()[l_ce] == Approx(10) );
    REQUIRE( l_waveProp.getMomentumX()[l_ce] == Approx(0) );
  }
}

//Helper for CSV reading
struct TestCaseData {
  double hLeft, hRight;
  double huLeft, huRight;
  double hStar;
};

std::vector<TestCaseData> readCsv(const std::string& path) {
  std::vector<TestCaseData> data;
  std::ifstream file(path);

  if (!file.is_open()) {
    throw std::runtime_error("Could not open CSV file: " + path);
  }

  std::string line;

  // skip header
  std::getline(file, line);

  while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string value;
    TestCaseData row;

    std::getline(ss, value, ','); row.hLeft = std::stod(value);
    std::getline(ss, value, ','); row.hRight = std::stod(value);
    std::getline(ss, value, ','); row.huLeft = std::stod(value);
    std::getline(ss, value, ','); row.huRight = std::stod(value);
    std::getline(ss, value, ','); row.hStar = std::stod(value);

    data.push_back(row);
  }

  return data;
}

TEST_CASE("WavePropagation1d CSV validation", "[WavePropFWave1d]") {
  /*
   * Test case:
   * Reads the csv file and then runs 100 time steps of the simulation.
   * Checks if the height of the middle cell is with 10e-3 window of the value provided in middle_states.csv 
   */

  auto testData = readCsv(MIDDLE_STATES_CSV);

  REQUIRE_FALSE(testData.empty());

  for (const auto& tc : testData) {

    tsunami_lab::patches::WavePropagation1d waveProp(100, true);

    // left half
    for (std::size_t i = 0; i < 50; i++) {
      waveProp.setHeight(i, 0, tc.hLeft);
      waveProp.setMomentumX(i, 0, tc.huLeft);
    }

    // right half
    for (std::size_t i = 50; i < 100; i++) {
      waveProp.setHeight(i, 0, tc.hRight);
      waveProp.setMomentumX(i, 0, tc.huRight);
    }

    // compute timestep scaling
    tsunami_lab::t_real speedMax = std::sqrt(9.81 * std::max(tc.hLeft, tc.hRight));
    tsunami_lab::t_real scaling = 0.5 / speedMax;

    // run simulation
    for (int step = 0; step < 100; step++) {
      waveProp.setGhostOutflow();
      waveProp.timeStep(scaling);
    }

    // check result at interface cell
    REQUIRE(waveProp.getHeight()[49] == Approx(tc.hStar).epsilon(1e-3));
  }
}