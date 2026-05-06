#include <catch2/catch.hpp>
#include "../constants.h"
#include "../patches/WavePropagation1d.h"
#include <sstream>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <filesystem>
#include "Stations.h"

#ifndef STATIONS_CSV
#define STATIONS_CSV "./stations/Stations.csv"
#endif

#ifndef STATIONS_OUTPUT_DIR
#define STATIONS_OUTPUT_DIR "./stations/output"
#endif

TEST_CASE("Stations reads mock configuration and writes CSV output with frequency control.", "[Stations]") {
  // Save original stations configuration file
  std::string stations_config = STATIONS_CSV;
  std::string stations_backup = stations_config + ".backup";
  
  // Back up the original config if it exists
  if (std::filesystem::exists(stations_config)) {
    std::filesystem::copy_file(stations_config, stations_backup, std::filesystem::copy_options::overwrite_existing);
  }

  // Create mock configuration file with test stations
  {
    std::ofstream config(stations_config);
    config << "TestStation1,30,0\n";
    config << "TestStation2,70,0\n";
    config.close();
  }

  // Create a 1D wave propagation solver with test data
  tsunami_lab::patches::WavePropagation1d wave_prop(100, false, tsunami_lab::patches::WavePropagation1d::BoundaryCondition::GhostOutflow);

  // Set up initial heights and momenta for testing
  for (tsunami_lab::t_idx i = 0; i < 100; i++) {
    wave_prop.setHeight(i, 0, 10.0f + i * 0.1f);
    wave_prop.setMomentumX(i, 0, 2.0f * i);
  }

  // Create Stations object with output frequency of 1.0 second
  // This will now read the mock configuration
  tsunami_lab::io::Stations stations(1.0);

  // Test 1: Write at time 1.0 (first write should happen since time >= frequency)
  stations.writeToCSV(1.0, 1.0, 1.0, &wave_prop);

  // Verify that output directory exists
  std::string output_dir = STATIONS_OUTPUT_DIR;
  REQUIRE(std::filesystem::exists(output_dir));

  // Verify that output files were created from mock configuration
  std::string station1_file = std::string(output_dir) + "/TestStation1.csv";
  std::string station2_file = std::string(output_dir) + "/TestStation2.csv";
  
  REQUIRE(std::filesystem::exists(station1_file));
  REQUIRE(std::filesystem::exists(station2_file));

  // Test 2: Read TestStation1.csv and verify format and content
  {
    std::ifstream file(station1_file);
    REQUIRE(file.is_open());

    std::string line;
    std::getline(file, line);
    // First line should be the header
    REQUIRE(line.find("time,water_height,momentum_x,momentum_y") != std::string::npos);

    std::getline(file, line);
    // Second line should contain data from time 1.0
    REQUIRE(line.find("1") != std::string::npos);

    file.close();
  }

  // Test 3: Verify frequency control - write at time 1.5 should NOT happen (< 2.0)
  stations.writeToCSV(1.5, 1.0, 1.0, &wave_prop);

  // Count lines in TestStation1.csv (should still be 2: header + one data line)
  {
    std::ifstream file(station1_file);
    int line_count = 0;
    std::string line;
    while (std::getline(file, line)) {
      line_count++;
    }
    file.close();
    REQUIRE(line_count == 2); // header + data from time 1.0
  }

  // Test 4: Write at time 2.0 should happen (>= 1.0 since last write at 1.0)
  stations.writeToCSV(2.0, 1.0, 1.0, &wave_prop);

  // Count lines in TestStation1.csv (should now be 3: header + two data lines)
  {
    std::ifstream file(station1_file);
    int line_count = 0;
    std::string line;
    while (std::getline(file, line)) {
      line_count++;
    }
    file.close();
    REQUIRE(line_count == 3); // header + data from time 1.0 + data from time 2.0
  }

  // Test 5: Verify the data is written with correct time values
  {
    std::ifstream file(station1_file);
    std::string line;
    
    // Skip header
    std::getline(file, line);
    
    // Read first data line
    std::getline(file, line);
    // Extract time value (should be 1)
    int comma_pos = line.find(',');
    std::string time_str = line.substr(0, comma_pos);
    REQUIRE(time_str == "1");

    // Read second data line
    std::getline(file, line);
    // Extract time value (should be 2)
    comma_pos = line.find(',');
    time_str = line.substr(0, comma_pos);
    REQUIRE(time_str == "2");

    file.close();
  }

  // Test 6: Verify TestStation2.csv also has the same data structure
  {
    std::ifstream file(station2_file);
    REQUIRE(file.is_open());

    std::string line;
    std::getline(file, line);
    REQUIRE(line.find("time,water_height,momentum_x,momentum_y") != std::string::npos);

    int line_count = 1; // header
    while (std::getline(file, line)) {
      line_count++;
    }
    file.close();
    REQUIRE(line_count == 3); // header + 2 data lines
  }

  // Cleanup: Restore original stations configuration file
  if (std::filesystem::exists(stations_backup)) {
    std::filesystem::copy_file(stations_backup, stations_config, std::filesystem::copy_options::overwrite_existing);
    std::filesystem::remove(stations_backup);
  }
}