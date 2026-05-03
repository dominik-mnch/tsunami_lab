/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Entry-point for simulations.
 **/
#include "patches/WavePropagation1d.h"
#include "patches/WavePropagation2d.h"
#include "setups/CircularDamBreak2d/CircularDamBreak2d.h"
#include "setups/DamBreak1d/DamBreak1d.h"
#include "setups/RareRare1d/RareRare1d.h"
#include "setups/ShockShock1d/ShockShock1d.h"
#include "setups/Subcritical1d/Subcritical1d.h"
#include "setups/Supercritical1d/Supercritical1d.h"
#include "setups/TsunamiEvent1d/TsunamiEvent1d.h"
#include "io/Csv.h"
#include <cctype>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cmath>
#include <fstream>
#include <limits>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

int main( int i_argc, char *i_argv[] ) {
  // number of cells in x- and y-direction
  tsunami_lab::t_idx l_nx = 0;
  tsunami_lab::t_idx l_ny = 1;

  // set cell size, domain, and end time
  tsunami_lab::t_real l_dxy = 1;
  tsunami_lab::t_real l_domainStart = 0.0;
  tsunami_lab::t_real l_domainEnd = 50000.0;
  tsunami_lab::t_real l_domainSize = l_domainEnd - l_domainStart;
  tsunami_lab::t_real l_endTime = 1.25;
  bool l_useFWaveSolver = true;
  bool l_useWavePropagation1d = false;
  std::string l_setupName;

  std::cout << "####################################" << std::endl;
  std::cout << "### Tsunami Lab                  ###" << std::endl;
  std::cout << "###                              ###" << std::endl;
  std::cout << "### https://scalable.uni-jena.de ###" << std::endl;
  std::cout << "####################################" << std::endl;

  auto l_printUsage = []() {
    std::cerr << "usage:" << std::endl;
    std::cerr << "  ./build/tsunami_lab NX NY DOMAIN_LOWER DOMAIN_UPPER SOLVER_MODE END_TIME PROPAGATION [PROP_PARAMS] SETUP [SETUP_PARAMS]" << std::endl;
    std::cerr << "" << std::endl;
    std::cerr << "general arguments:" << std::endl;
    std::cerr << "  NX, NY       number of cells in x- and y-direction." << std::endl;
    std::cerr << "  DOMAIN_LOWER lower bound of the domain in x-direction." << std::endl;
    std::cerr << "  DOMAIN_UPPER upper bound of the domain in x-direction." << std::endl;
    std::cerr << "  SOLVER_MODE  1 for FWaveSolver, 0 for RoeSolver." << std::endl;
    std::cerr << "  END_TIME     simulation end time in seconds." << std::endl;
    std::cerr << "" << std::endl;
    std::cerr << "propagation modes:" << std::endl;
    std::cerr << "  1d <boundary_mode>" << std::endl;
    std::cerr << "    boundary_mode: outflow | right | left | both" << std::endl;
    std::cerr << "  2d" << std::endl;
    std::cerr << "" << std::endl;
    std::cerr << "setups:" << std::endl;
    std::cerr << "  dam_break_1d hL huL hR huR xDam" << std::endl;
    std::cerr << "  shock_shock_1d h hu xDisc [useParabola [offset center scale]]" << std::endl;
    std::cerr << "  rare_rare_1d  h hu xDisc [useParabola [offset center scale]]" << std::endl;
    std::cerr << "  subcritical_1d" << std::endl;
    std::cerr << "  supercritical_1d" << std::endl;
    std::cerr << "  tsunami_event_1d" << std::endl;
    std::cerr << "  circular_dam_break_2d hIn hOut huIn hvIn huOut hvOut xMid yMid radius" << std::endl;
  };

  auto l_parseReal = []( std::string const & i_value, std::string const & i_name ) {
    char * l_endPtr = nullptr;
    double l_val = std::strtod( i_value.c_str(), &l_endPtr );
    if( l_endPtr == i_value.c_str() || *l_endPtr != '\0' ) {
      throw std::invalid_argument( "invalid real value for '" + i_name + "': " + i_value );
    }
    return static_cast<tsunami_lab::t_real>( l_val );
  };

  auto l_parseIdx = []( std::string const & i_value, std::string const & i_name ) {
    char * l_endPtr = nullptr;
    long long l_val = std::strtoll( i_value.c_str(), &l_endPtr, 10 );
    if( l_endPtr == i_value.c_str() || *l_endPtr != '\0' || l_val < 1 ) {
      throw std::invalid_argument( "invalid positive integer for '" + i_name + "': " + i_value );
    }
    return static_cast<tsunami_lab::t_idx>( l_val );
  };

  auto l_parseBool = []( std::string i_value, std::string const & i_name ) {
    for( char & l_ch: i_value ) l_ch = static_cast<char>( std::tolower( static_cast<unsigned char>( l_ch ) ) );
    if( i_value == "1" || i_value == "true" ) return true;
    if( i_value == "0" || i_value == "false" ) return false;
    throw std::invalid_argument( "invalid boolean for '" + i_name + "': " + i_value );
  };

  if( i_argc < 9 ) {
    std::cerr << "invalid number of arguments" << std::endl;
    l_printUsage();
    return EXIT_FAILURE;
  }

  tsunami_lab::setups::Setup *l_setup = nullptr;
  tsunami_lab::patches::WavePropagation *l_waveProp = nullptr;

  try {
    l_nx = l_parseIdx( i_argv[1], "NX" );
    l_ny = l_parseIdx( i_argv[2], "NY" );
    l_domainStart = l_parseReal( i_argv[3], "DOMAIN_LOWER" );
    l_domainEnd = l_parseReal( i_argv[4], "DOMAIN_UPPER" );
    l_domainSize = l_domainEnd - l_domainStart;

    int l_solverMode = static_cast<int>( l_parseIdx( i_argv[5], "SOLVER_MODE" ) );
    if( l_solverMode != 0 && l_solverMode != 1 ) {
      throw std::invalid_argument( "SOLVER_MODE must be 0 or 1" );
    }
    l_useFWaveSolver = (l_solverMode == 1);

    l_endTime = l_parseReal( i_argv[6], "END_TIME" );

    if( l_domainSize <= 0 || l_endTime <= 0 ) {
      throw std::invalid_argument( "DOMAIN_UPPER must be > DOMAIN_LOWER and END_TIME must be > 0" );
    }

    std::string l_propagationMode = i_argv[7];
    for( char & l_ch: l_propagationMode ) l_ch = static_cast<char>( std::tolower( static_cast<unsigned char>( l_ch ) ) );

    int l_setupArgStart = 0;
    tsunami_lab::patches::WavePropagation1d::BoundaryCondition l_boundaryCondition =
      tsunami_lab::patches::WavePropagation1d::BoundaryCondition::GhostOutflow;

    if( l_propagationMode == "1d" ) {
      l_useWavePropagation1d = true;
      if( i_argc < 10 ) {
        throw std::invalid_argument( "propagation mode '1d' requires a boundary mode and setup" );
      }

      std::string l_boundaryMode = i_argv[8];
      for( char & l_ch: l_boundaryMode ) l_ch = static_cast<char>( std::tolower( static_cast<unsigned char>( l_ch ) ) );

      if( l_boundaryMode == "outflow" ) {
        l_boundaryCondition = tsunami_lab::patches::WavePropagation1d::BoundaryCondition::GhostOutflow;
      }
      else if( l_boundaryMode == "right" ) {
        l_boundaryCondition = tsunami_lab::patches::WavePropagation1d::BoundaryCondition::BoundaryRight;
      }
      else if( l_boundaryMode == "left" ) {
        l_boundaryCondition = tsunami_lab::patches::WavePropagation1d::BoundaryCondition::BoundaryLeft;
      }
      else if( l_boundaryMode == "both" ) {
        l_boundaryCondition = tsunami_lab::patches::WavePropagation1d::BoundaryCondition::BoundaryBoth;
      }
      else {
        throw std::invalid_argument( "unknown 1d boundary mode: " + l_boundaryMode );
      }

      l_setupArgStart = 9;
      l_waveProp = new tsunami_lab::patches::WavePropagation1d( l_nx,
                                                                 l_useFWaveSolver,
                                                                 l_boundaryCondition );

      if( l_ny != 1 ) {
        std::cout << "warning: forcing NY=1 for 1d wave propagation" << std::endl;
        l_ny = 1;
      }
    }
    else if( l_propagationMode == "2d" ) {
      l_useWavePropagation1d = false;
      l_setupArgStart = 8;
      l_waveProp = new tsunami_lab::patches::WavePropagation2d( l_nx,
                                                                 l_ny,
                                                                 l_useFWaveSolver );
    }
    else {
      throw std::invalid_argument( "unknown propagation mode: " + l_propagationMode );
    }

    if( l_setupArgStart >= i_argc ) {
      throw std::invalid_argument( "missing setup name" );
    }

    l_setupName = i_argv[l_setupArgStart];
    for( char & l_ch: l_setupName ) l_ch = static_cast<char>( std::tolower( static_cast<unsigned char>( l_ch ) ) );

    std::vector<std::string> l_setupArgs;
    for( int l_ar = l_setupArgStart + 1; l_ar < i_argc; l_ar++ ) {
      l_setupArgs.emplace_back( i_argv[l_ar] );
    }

    bool l_setupIs2d = (l_setupName == "circular_dam_break_2d");
    if( l_useWavePropagation1d && l_setupIs2d ) {
      throw std::invalid_argument( "setup circular_dam_break_2d requires 2d propagation" );
    }

    if( l_setupName == "dam_break_1d" ) {
      if( l_setupArgs.size() != 5 ) {
        throw std::invalid_argument( "dam_break_1d expects 5 parameters" );
      }

      l_setup = new tsunami_lab::setups::DamBreak1d( l_parseReal( l_setupArgs[0], "hL" ),
                                                      l_parseReal( l_setupArgs[1], "huL" ),
                                                      l_parseReal( l_setupArgs[2], "hR" ),
                                                      l_parseReal( l_setupArgs[3], "huR" ),
                                                      l_parseReal( l_setupArgs[4], "xDam" ) );
    }
    else if( l_setupName == "shock_shock_1d" ) {
      if( l_setupArgs.size() == 3 ) {
        l_setup = new tsunami_lab::setups::ShockShock1d( l_parseReal( l_setupArgs[0], "h" ),
                                                          l_parseReal( l_setupArgs[1], "hu" ),
                                                          l_parseReal( l_setupArgs[2], "xDisc" ) );
      }
      else if( l_setupArgs.size() == 4 ) {
        l_setup = new tsunami_lab::setups::ShockShock1d( l_parseReal( l_setupArgs[0], "h" ),
                                                          l_parseReal( l_setupArgs[1], "hu" ),
                                                          l_parseReal( l_setupArgs[2], "xDisc" ),
                                                          l_parseBool( l_setupArgs[3], "useParabola" ) );
      }
      else if( l_setupArgs.size() == 7 ) {
        l_setup = new tsunami_lab::setups::ShockShock1d( l_parseReal( l_setupArgs[0], "h" ),
                                                          l_parseReal( l_setupArgs[1], "hu" ),
                                                          l_parseReal( l_setupArgs[2], "xDisc" ),
                                                          l_parseBool( l_setupArgs[3], "useParabola" ),
                                                          l_parseReal( l_setupArgs[4], "offset" ),
                                                          l_parseReal( l_setupArgs[5], "center" ),
                                                          l_parseReal( l_setupArgs[6], "scale" ) );
      }
      else {
        throw std::invalid_argument( "shock_shock_1d expects 3, 4 or 7 parameters" );
      }
    }
    else if( l_setupName == "rare_rare_1d" ) {
      if( l_setupArgs.size() == 3 ) {
        l_setup = new tsunami_lab::setups::RareRare1d( l_parseReal( l_setupArgs[0], "h" ),
                                                        l_parseReal( l_setupArgs[1], "hu" ),
                                                        l_parseReal( l_setupArgs[2], "xDisc" ) );
      }
      else if( l_setupArgs.size() == 4 ) {
        l_setup = new tsunami_lab::setups::RareRare1d( l_parseReal( l_setupArgs[0], "h" ),
                                                        l_parseReal( l_setupArgs[1], "hu" ),
                                                        l_parseReal( l_setupArgs[2], "xDisc" ),
                                                        l_parseBool( l_setupArgs[3], "useParabola" ) );
      }
      else if( l_setupArgs.size() == 7 ) {
        l_setup = new tsunami_lab::setups::RareRare1d( l_parseReal( l_setupArgs[0], "h" ),
                                                        l_parseReal( l_setupArgs[1], "hu" ),
                                                        l_parseReal( l_setupArgs[2], "xDisc" ),
                                                        l_parseBool( l_setupArgs[3], "useParabola" ),
                                                        l_parseReal( l_setupArgs[4], "offset" ),
                                                        l_parseReal( l_setupArgs[5], "center" ),
                                                        l_parseReal( l_setupArgs[6], "scale" ) );
      }
      else {
        throw std::invalid_argument( "rare_rare_1d expects 3, 4 or 7 parameters" );
      }
    }
    else if( l_setupName == "subcritical_1d" ) {
      if( !l_setupArgs.empty() ) {
        throw std::invalid_argument( "subcritical_1d expects no parameters" );
      }
      l_setup = new tsunami_lab::setups::Subcritical1d();
    }
    else if( l_setupName == "supercritical_1d" ) {
      if( !l_setupArgs.empty() ) {
        throw std::invalid_argument( "supercritical_1d expects no parameters" );
      }
      l_setup = new tsunami_lab::setups::Supercritical1d();
    }
    else if( l_setupName == "tsunami_event_1d" ) {
      if( !l_setupArgs.empty() ) {
        throw std::invalid_argument( "tsunami_event_1d expects no parameters" );
      }
      l_setup = new tsunami_lab::setups::TsunamiEvent1d();
    }
    else if( l_setupName == "circular_dam_break_2d" ) {
      if( l_setupArgs.size() != 9 ) {
        throw std::invalid_argument( "circular_dam_break_2d expects 9 parameters" );
      }
      l_setup = new tsunami_lab::setups::CircularDamBreak2d( l_parseReal( l_setupArgs[0], "hIn" ),
                                                              l_parseReal( l_setupArgs[1], "hOut" ),
                                                              l_parseReal( l_setupArgs[2], "huIn" ),
                                                              l_parseReal( l_setupArgs[3], "hvIn" ),
                                                              l_parseReal( l_setupArgs[4], "huOut" ),
                                                              l_parseReal( l_setupArgs[5], "hvOut" ),
                                                              l_parseReal( l_setupArgs[6], "xMid" ),
                                                              l_parseReal( l_setupArgs[7], "yMid" ),
                                                              l_parseReal( l_setupArgs[8], "radius" ) );
    }
    else {
      throw std::invalid_argument( "unknown setup: " + l_setupName );
    }

    // Calculate cell size using the domain argument
    l_dxy = l_domainSize / l_nx;
  }
  catch( std::exception const & i_ex ) {
    std::cerr << "invalid configuration: " << i_ex.what() << std::endl;
    l_printUsage();
    delete l_setup;
    delete l_waveProp;
    return EXIT_FAILURE;
  }

  if( l_setup == nullptr || l_waveProp == nullptr ) {
    std::cerr << "internal error: setup or propagation was not created" << std::endl;
    delete l_setup;
    delete l_waveProp;
    return EXIT_FAILURE;
  }

  if( l_nx < 1 || l_ny < 1 || l_domainSize <= 0 || l_endTime <= 0 ) {
      std::cerr << "invalid configuration: check cell counts, domain size, and end time" << std::endl;
      delete l_setup;
      delete l_waveProp;
      return EXIT_FAILURE;
  }

  std::cout << "runtime configuration" << std::endl;
  std::cout << "  number of cells in x-direction: " << l_nx << std::endl;
  std::cout << "  number of cells in y-direction: " << l_ny << std::endl;
  std::cout << "  domain lower bound x:          " << l_domainStart << std::endl;
  std::cout << "  domain upper bound x:          " << l_domainEnd << std::endl;
  std::cout << "  domain size:                    " << l_domainSize << std::endl;
  std::cout << "  cell size:                      " << l_dxy << std::endl;
  std::cout << "  propagation:                    " << (l_useWavePropagation1d ? "1d" : "2d") << std::endl;
  std::cout << "  setup:                          " << l_setupName << std::endl;
  std::cout << "  solver:                         " << (l_useFWaveSolver ? "FWave" : "Roe") << std::endl;
  std::cout << "  end time:                       " << l_endTime << std::endl;

  // clean up solutions directory
  std::cout << "cleaning solutions directory" << std::endl;
  std::string l_solutionsDir = "./solutions";
  if( std::filesystem::exists( l_solutionsDir ) ) {
    for( const auto& l_entry : std::filesystem::directory_iterator( l_solutionsDir ) ) {
      if( l_entry.is_regular_file() && l_entry.path().extension() == ".csv" ) {
        std::filesystem::remove( l_entry.path() );
      }
    }
  }

  // maximum observed height in the setup
  tsunami_lab::t_real l_hMax = std::numeric_limits< tsunami_lab::t_real >::lowest();

  // set up solver
  for( tsunami_lab::t_idx l_cy = 0; l_cy < l_ny; l_cy++ ) {
    tsunami_lab::t_real l_y = l_domainStart + l_cy * l_dxy;

    for( tsunami_lab::t_idx l_cx = 0; l_cx < l_nx; l_cx++ ) {
      tsunami_lab::t_real l_x = l_domainStart + l_cx * l_dxy;

      // get initial values of the setup
      tsunami_lab::t_real l_h = l_setup->getHeight( l_x,
                                                    l_y );
      l_hMax = std::max( l_h, l_hMax );

      tsunami_lab::t_real l_hu = l_setup->getMomentumX( l_x,
                                                        l_y );
      tsunami_lab::t_real l_hv = l_setup->getMomentumY( l_x,
                                                        l_y );
      tsunami_lab::t_real l_b = l_setup->getBathymetry( l_x,
                                                        l_y );

      // set initial values in wave propagation solver
      l_waveProp->setHeight( l_cx,
                             l_cy,
                             l_h );

      l_waveProp->setMomentumX( l_cx,
                                l_cy,
                                l_hu );

      l_waveProp->setMomentumY( l_cx,
                                l_cy,
                                l_hv );

      l_waveProp->setBathymetry( l_cx,
                                 l_cy,
                                 l_b );

    }
  }

  // derive maximum wave speed in setup; the momentum is ignored
  tsunami_lab::t_real l_speedMax = std::sqrt( 9.81 * std::max( l_hMax, static_cast<tsunami_lab::t_real>( 0 ) ) );
  if( l_speedMax <= 0 ) {
    std::cerr << "invalid setup state: maximum wave speed is zero" << std::endl;
    delete l_setup;
    delete l_waveProp;
    return EXIT_FAILURE;
  }

  // derive constant time step; changes at simulation time are ignored
  tsunami_lab::t_real l_dt = 0.5 * l_dxy / l_speedMax;

  // derive scaling for a time step
  tsunami_lab::t_real l_scaling = l_dt / l_dxy;

  // set up time and print control (endTime is set at the top)
  tsunami_lab::t_idx  l_timeStep = 0;
  tsunami_lab::t_idx  l_nOut = 0;
  tsunami_lab::t_real l_simTime = 0;


  std::cout << "entering time loop" << std::endl;

  // iterate over time
  while( l_simTime < l_endTime ){
    if( l_timeStep % 25 == 0 ) {
      std::cout << "  simulation time / #time steps: "
                << l_simTime << " / " << l_timeStep << std::endl;

      std::string l_path = "./solutions/solution_" + std::to_string(l_nOut) + ".csv";
      std::cout << "  writing wave field to " << l_path << std::endl;

      std::ofstream l_file;
      l_file.open( l_path  );  

      tsunami_lab::t_real const * l_bath = l_waveProp->getBathymetry();
      if( l_useWavePropagation1d ) l_bath = l_bath + 1;

      tsunami_lab::io::Csv::write( l_dxy,
                                   l_domainStart,
                                   l_domainStart,
                                   l_nx,
                                   l_ny,
                                   l_useWavePropagation1d ? 1 : l_waveProp->getStride(),
                                   l_simTime,
                                   l_waveProp->getHeight(),
                                   l_bath,
                                   l_waveProp->getMomentumX(),
                                   l_useWavePropagation1d ? nullptr : l_waveProp->getMomentumY(),
                                   l_file );
      l_file.close();
      l_nOut++;
    }
    l_waveProp->setGhostOutflow();
    l_waveProp->timeStep( l_scaling );

    l_timeStep++;
    l_simTime += l_dt;
  }

  std::cout << "finished time loop" << std::endl;

  // free memory
  std::cout << "freeing memory" << std::endl;
  delete l_setup;
  delete l_waveProp;

  std::cout << "finished, exiting" << std::endl;
  return EXIT_SUCCESS;
}
