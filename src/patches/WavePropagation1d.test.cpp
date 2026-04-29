/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Unit tests for the one-dimensional wave propagation patch.
 **/
#include <catch2/catch.hpp>
#include "WavePropagation1d.h"
#include "../solvers/Roe.h"
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
  tsunami_lab::patches::WavePropagation1d m_waveProp( 100 , false, tsunami_lab::patches::WavePropagation1d::BoundaryCondition::GhostOutflow );

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

  tsunami_lab::patches::WavePropagation1d l_waveProp( 100, true, tsunami_lab::patches::WavePropagation1d::BoundaryCondition::GhostOutflow );

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

TEST_CASE( "Test non-ghost-outflow boundary condition changes boundary behavior.", "[WavePropBoundaryBehaviourChange1d]" ) {
  // Same initial state, but different boundary conditions on the right side.
  tsunami_lab::patches::WavePropagation1d l_outflow( 20,
                                                     false,
                                                     tsunami_lab::patches::WavePropagation1d::BoundaryCondition::GhostOutflow );
  tsunami_lab::patches::WavePropagation1d l_reflectRight( 20,
                                                          false,
                                                          tsunami_lab::patches::WavePropagation1d::BoundaryCondition::BoundaryRight );

  for( std::size_t l_ce = 0; l_ce < 20; l_ce++ ) {
    l_outflow.setHeight( l_ce,
                         0,
                         5 );
    l_outflow.setMomentumX( l_ce,
                            0,
                            2 );

    l_reflectRight.setHeight( l_ce,
                              0,
                              5 );
    l_reflectRight.setMomentumX( l_ce,
                                 0,
                                 2 );
  }

  l_outflow.setGhostOutflow();
  l_reflectRight.setGhostOutflow();

  l_outflow.timeStep( 0.1 );
  l_reflectRight.timeStep( 0.1 );

  // Left boundary uses outflow in both configurations: first cell should match.
  REQUIRE( l_reflectRight.getHeight()[0] == Approx( l_outflow.getHeight()[0] ) );
  REQUIRE( l_reflectRight.getMomentumX()[0] == Approx( l_outflow.getMomentumX()[0] ) );

  // Right boundary differs (outflow vs reflecting): last cell should differ.
  REQUIRE( l_reflectRight.getHeight()[19] != Approx( l_outflow.getHeight()[19] ) );
  REQUIRE( l_reflectRight.getMomentumX()[19] != Approx( l_outflow.getMomentumX()[19] ) );
}

TEST_CASE( "Test BoundaryRight computes correct right-edge update with Roe solver.", "[WavePropBoundaryCalculation1d]" ) {
  tsunami_lab::patches::WavePropagation1d l_waveProp( 3,
                                                      false,
                                                      tsunami_lab::patches::WavePropagation1d::BoundaryCondition::BoundaryRight );

  for( std::size_t l_ce = 0; l_ce < 3; l_ce++ ) {
    l_waveProp.setHeight( l_ce,
                          0,
                          5 );
    l_waveProp.setMomentumX( l_ce,
                             0,
                             2 );
  }

  l_waveProp.setGhostOutflow();

  tsunami_lab::t_real const l_hLast = 5;
  tsunami_lab::t_real const l_huLast = 2;
  tsunami_lab::t_real const l_scaling = 0.1;

  tsunami_lab::t_real l_netUpdates[2][2];
  tsunami_lab::solvers::Roe::netUpdates( l_hLast,
                                         l_hLast,
                                         l_huLast,
                                         -l_huLast,
                                         l_netUpdates[0],
                                         l_netUpdates[1] );

  l_waveProp.timeStep( l_scaling );

  // Interior edge 1|2 has identical states, so only the right boundary edge contributes.
  REQUIRE( l_waveProp.getHeight()[2] == Approx( l_hLast - l_scaling * l_netUpdates[0][0] ) );
  REQUIRE( l_waveProp.getMomentumX()[2] == Approx( l_huLast - l_scaling * l_netUpdates[0][1] ) );
}

TEST_CASE( "Test wet-dry shoreline reflection updates only wet cell.", "[WavePropWetDryReflect1d]" ) {
  tsunami_lab::patches::WavePropagation1d l_waveProp( 3,
                                                      false,
                                                      tsunami_lab::patches::WavePropagation1d::BoundaryCondition::GhostOutflow );

  // Cell 1 is dry, cells 0 and 2 are wet to avoid dry-dry interfaces.
  l_waveProp.setHeight( 0, 0, 5 );
  l_waveProp.setMomentumX( 0, 0, 2 );
  l_waveProp.setBathymetry( 0, 0, -1 );

  l_waveProp.setHeight( 1, 0, 0 );
  l_waveProp.setMomentumX( 1, 0, 0 );
  l_waveProp.setBathymetry( 1, 0, 1 );

  l_waveProp.setHeight( 2, 0, 7 );
  l_waveProp.setMomentumX( 2, 0, -1 );
  l_waveProp.setBathymetry( 2, 0, -2 );

  l_waveProp.setGhostOutflow();

  tsunami_lab::t_real const l_scaling = 0.1;
  tsunami_lab::t_real l_expectedNetUpdatesL[2][2];
  tsunami_lab::t_real l_expectedNetUpdatesR[2][2];

  // Left wet/dry interface: mirrored state of cell 0 against dry cell 1.
  tsunami_lab::solvers::Roe::netUpdates( 5,
                                         5,
                                         2,
                                         -2,
                                         l_expectedNetUpdatesL[0],
                                         l_expectedNetUpdatesL[1] );

  // Right wet/dry interface: mirrored state of cell 2 against dry cell 1.
  tsunami_lab::solvers::Roe::netUpdates( 7,
                                         7,
                                         1,
                                         -1,
                                         l_expectedNetUpdatesR[0],
                                         l_expectedNetUpdatesR[1] );

  l_waveProp.timeStep( l_scaling );

  REQUIRE( l_waveProp.getHeight()[0] == Approx( 5 - l_scaling * l_expectedNetUpdatesL[0][0] ) );
  REQUIRE( l_waveProp.getMomentumX()[0] == Approx( 2 - l_scaling * l_expectedNetUpdatesL[0][1] ) );

  // Dry cell stays unchanged because dry-side net updates are suppressed.
  REQUIRE( l_waveProp.getHeight()[1] == Approx( 0 ) );
  REQUIRE( l_waveProp.getMomentumX()[1] == Approx( 0 ) );

  REQUIRE( l_waveProp.getHeight()[2] == Approx( 7 - l_scaling * l_expectedNetUpdatesR[1][0] ) );
  REQUIRE( l_waveProp.getMomentumX()[2] == Approx( -1 - l_scaling * l_expectedNetUpdatesR[1][1] ) );
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

    tsunami_lab::patches::WavePropagation1d waveProp(100, true, tsunami_lab::patches::WavePropagation1d::BoundaryCondition::GhostOutflow);

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