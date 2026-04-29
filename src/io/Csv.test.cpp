/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Unit tests for the CSV-interface.
 **/
#include <catch2/catch.hpp>
#include "../constants.h"
#include <sstream>
#include <fstream>
#include <cstdio>
#define private public
#include "Csv.h"
#undef public

TEST_CASE( "Test the CSV-writer for 1D settings.", "[CsvWrite1d]" ) {
  // define a simple example
  tsunami_lab::t_real l_h[7]  = { 0, 1, 2, 3, 4, 5, 6 };
  tsunami_lab::t_real l_b[7]  = { 0, -1, -2, -3, -4, -5, -6 };
  tsunami_lab::t_real l_hu[7] = { 6, 5, 4, 3, 2, 1, 0 };

  std::stringstream l_stream0;
  tsunami_lab::io::Csv::write( 0.5,
                               5,
                               1,
                               7,
                               1.0,
                               l_h+1,
                               l_b+1,
                               l_hu+1,
                               nullptr,
                               l_stream0 );

  std::string l_ref0 = R"V0G0N(time,x,y,height,bathymetry,momentum_x
1,0.25,0.25,1,-1,5
1,0.75,0.25,2,-2,4
1,1.25,0.25,3,-3,3
1,1.75,0.25,4,-4,2
1,2.25,0.25,5,-5,1
)V0G0N";

  REQUIRE( l_stream0.str().size() == l_ref0.size() );
  REQUIRE( l_stream0.str() == l_ref0 );
}

TEST_CASE( "Test the CSV-writer for 2D settings.", "[CsvWrite2d]" ) {
  // define a simple example
  tsunami_lab::t_real l_h[16]  = {  0,  1,  2,  3,
                                    4,  5,  6,  7,
                                    8,  9, 10, 11,
                                   12, 13, 14, 15 };
  tsunami_lab::t_real l_hu[16] = { 15, 14, 13, 12,
                                   11, 10,  9,  8,
                                    7,  6,  5,  4,
                                    3,  2,  1,  0 };
  tsunami_lab::t_real l_b[16]  = {  0,  1,  2,  3,
                                    4,  5,  6,  7,
                                    8,  9, 10, 11,
                                   12, 13, 14, 15 };
  tsunami_lab::t_real l_hv[16] = {  0,  4,  8, 12,
                                    1,  5,  9, 13,
                                    2,  6, 10, 14,
                                    3,  7, 11, 15 };

  std::stringstream l_stream1;
  tsunami_lab::io::Csv::write( 10,
                               2,
                               2,
                               4,
                               2.0,
                               l_h+4+1,
                               l_b+4+1,
                               l_hu+4+1,
                               l_hv+4+1,
                               l_stream1 );

  std::string l_ref1 = R"V0G0N(time,x,y,height,bathymetry,momentum_x,momentum_y
2,5,5,5,5,10,5
2,15,5,6,6,9,9
2,5,15,9,9,6,6
2,15,15,10,10,5,10
)V0G0N";

  REQUIRE( l_stream1.str().size() == l_ref1.size() );
  REQUIRE( l_stream1.str() == l_ref1 );
}

TEST_CASE( "Test the CSV-reader for bathymetry data.", "[CsvReadBathymetry]" ) {
  // create a temporary CSV file with bathymetry data
  std::string l_testFile = "test_bathymetry.csv";
  std::ofstream l_outFile( l_testFile );
  
  l_outFile << "lat,long,x,height\n";
  l_outFile << "141.024949,37.316569,0,14.7254650696\n";
  l_outFile << "141.02756337,37.3166237288,0.231656660144,-6.93627624916\n";
  l_outFile << "141.030177745,37.3166784,0.463313319336,-7.48690520552\n";
  
  l_outFile.close();
  
  std::vector<tsunami_lab::t_real> l_x;
  std::vector<tsunami_lab::t_real> l_b;
  
  bool l_success = tsunami_lab::io::Csv::readBathymetry( l_testFile, l_x, l_b );
  
  REQUIRE( l_success == true );
  REQUIRE( l_x.size() == 3 );
  REQUIRE( l_b.size() == 3 );
  
  // check the values are read correctly
  REQUIRE( l_x[0] == Approx( 0.0 ) );
  REQUIRE( l_x[1] == Approx( 0.231656660144 ) );
  REQUIRE( l_x[2] == Approx( 0.463313319336 ) );
  
  REQUIRE( l_b[0] == Approx( 14.7254650696 ) );
  REQUIRE( l_b[1] == Approx( -6.93627624916 ) );
  REQUIRE( l_b[2] == Approx( -7.48690520552 ) );
  
  // clean up
  std::remove( l_testFile.c_str() );
}

TEST_CASE( "Test the CSV-reader for non-existent file.", "[CsvReadBathymetryNonExistent]" ) {
  std::vector<tsunami_lab::t_real> l_x;
  std::vector<tsunami_lab::t_real> l_b;
  
  bool l_success = tsunami_lab::io::Csv::readBathymetry( "nonexistent_file.csv", l_x, l_b );
  
  REQUIRE( l_success == false );
  REQUIRE( l_x.size() == 0 );
  REQUIRE( l_b.size() == 0 );
}

TEST_CASE( "Test the CSV-reader for invalid format.", "[CsvReadBathymetryInvalidFormat]" ) {
  // create a temporary CSV file with invalid data
  std::string l_testFile = "test_bathymetry_invalid.csv";
  std::ofstream l_outFile( l_testFile );
  
  l_outFile << "lat,long,x,height\n";
  l_outFile << "141.024949,37.316569,0,14.7254650696\n";
  l_outFile << "invalid,37.3166237288,0.231656660144,-6.93627624916\n";
  
  l_outFile.close();
  
  std::vector<tsunami_lab::t_real> l_x;
  std::vector<tsunami_lab::t_real> l_b;
  
  bool l_success = tsunami_lab::io::Csv::readBathymetry( l_testFile, l_x, l_b );
  
  REQUIRE( l_success == false );
  
  // clean up
  std::remove( l_testFile.c_str() );
}