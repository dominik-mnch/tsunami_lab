/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Entry-point for simulations.
 **/
#include "patches/WavePropagation1d.h"
#include "patches/WavePropagation2d.h"
#include "setups/ArtificialTsunami2d/ArtificialTsunami2d.h"
#include "setups/CircularDamBreak2d/CircularDamBreak2d.h"
#include "setups/DamBreak1d/DamBreak1d.h"
#include "setups/RareRare1d/RareRare1d.h"
#include "setups/ShockShock1d/ShockShock1d.h"
#include "setups/Subcritical1d/Subcritical1d.h"
#include "setups/Supercritical1d/Supercritical1d.h"
#include "setups/TsunamiEvent1d/TsunamiEvent1d.h"
#include "setups/TsunamiEvent2d/TsunamiEvent2d.h"
#include "setups/Checkpoint/Checkpoint.h"
#include "io/Csv.h"
#include "io/NetCdf.h"
#include "io/Stations.h"
#include <chrono>
#include <cctype>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <limits>
#include <filesystem>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

int main( int i_argc, char *i_argv[] ) {
  // number of cells in x- and y-direction
  tsunami_lab::t_idx l_nx = 0;
  tsunami_lab::t_idx l_ny = 1;

  // set cell size, domain, and end time
  tsunami_lab::t_real l_dx = 1;
  tsunami_lab::t_real l_dy = 1;
  tsunami_lab::t_real l_domainStartX = 0.0;
  tsunami_lab::t_real l_domainEndX   = 50000.0;
  tsunami_lab::t_real l_domainStartY = 0.0;
  tsunami_lab::t_real l_domainEndY   = 50000.0;
  tsunami_lab::t_real l_domainSizeX  = l_domainEndX - l_domainStartX;
  tsunami_lab::t_real l_domainSizeY  = l_domainEndY - l_domainStartY;
  tsunami_lab::t_idx l_k = 1;
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
    std::cerr << "  ./build/tsunami_lab NX NY X_LOWER X_UPPER Y_LOWER Y_UPPER K SOLVER_MODE END_TIME PROPAGATION [PROP_PARAMS] SETUP [SETUP_PARAMS]" << std::endl;
    std::cerr << "  ./build/tsunami_lab <RES>m  X_LOWER X_UPPER Y_LOWER Y_UPPER K SOLVER_MODE END_TIME PROPAGATION [PROP_PARAMS] SETUP [SETUP_PARAMS]" << std::endl;
    std::cerr << "" << std::endl;
    std::cerr << "general arguments:" << std::endl;
    std::cerr << "  NX, NY       number of cells in x- and y-direction." << std::endl;
    std::cerr << "  <RES>m       resolution in meters (e.g. 250m); derives NX/NY from domain extents." << std::endl;
    std::cerr << "  X_LOWER      lower bound of the domain in x-direction." << std::endl;
    std::cerr << "  X_UPPER      upper bound of the domain in x-direction." << std::endl;
    std::cerr << "  Y_LOWER      lower bound of the domain in y-direction." << std::endl;
    std::cerr << "  Y_UPPER      upper bound of the domain in y-direction." << std::endl;
    std::cerr << "  K            for a coarse output, k*k cells get averaged." << std::endl;
    std::cerr << "  SOLVER_MODE  1 for FWaveSolver, 0 for RoeSolver." << std::endl;
    std::cerr << "  END_TIME     simulation end time in seconds." << std::endl;
    std::cerr << "" << std::endl;
    std::cerr << "propagation modes:" << std::endl;
    std::cerr << "  1d <boundary_mode>" << std::endl;
    std::cerr << "    boundary_mode: outflow | right | left | both" << std::endl;
    std::cerr << "  2d [boundary_mode]" << std::endl;
    std::cerr << "    boundary_mode (optional): outflow | bothx | bothy | all | left,right,bottom,top (any comma-separated subset)" << std::endl;
    std::cerr << "" << std::endl;
    std::cerr << "setups:" << std::endl;
    std::cerr << "  dam_break_1d hL huL hR huR xDam" << std::endl;
    std::cerr << "  shock_shock_1d h hu xDisc [useParabola [offset center scale]]" << std::endl;
    std::cerr << "  rare_rare_1d  h hu xDisc [useParabola [offset center scale]]" << std::endl;
    std::cerr << "  subcritical_1d" << std::endl;
    std::cerr << "  supercritical_1d" << std::endl;
    std::cerr << "  tsunami_event_1d" << std::endl;
    std::cerr << "  circular_dam_break_2d hIn hOut huIn hvIn huOut hvOut xMid yMid radius" << std::endl;
    std::cerr << "  artificial_tsunami_2d" << std::endl;
    std::cerr << "  tsunami_event_2d i_bathymetryFile i_displacementFile" << std::endl;
    std::cerr << "  tsunami2d i_bathymetryFile [i_displacementFile]" << std::endl;
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

  auto l_parseSolverMode = []( std::string const & i_value ) {
    char * l_endPtr = nullptr;
    long l_val = std::strtol( i_value.c_str(), &l_endPtr, 10 );
    if( l_endPtr == i_value.c_str() || *l_endPtr != '\0' || (l_val != 0 && l_val != 1) ) {
      throw std::invalid_argument( "SOLVER_MODE must be 0 or 1" );
    }
    return static_cast<int>( l_val );
  };

  auto l_parseBool = []( std::string i_value, std::string const & i_name ) {
    for( char & l_ch: i_value ) l_ch = static_cast<char>( std::tolower( static_cast<unsigned char>( l_ch ) ) );
    if( i_value == "1" || i_value == "true" ) return true;
    if( i_value == "0" || i_value == "false" ) return false;
    throw std::invalid_argument( "invalid boolean for '" + i_name + "': " + i_value );
  };

  auto l_parseResolution = []( std::string const & i_value, std::string const & i_name ) {
    if( i_value.empty() || i_value.back() != 'm' ) {
      throw std::invalid_argument( "resolution must end with 'm' (e.g. 250m): " + i_value );
    }
    std::string l_numStr = i_value.substr( 0, i_value.size() - 1 );
    char * l_endPtr = nullptr;
    double l_val = std::strtod( l_numStr.c_str(), &l_endPtr );
    if( l_endPtr == l_numStr.c_str() || *l_endPtr != '\0' || l_val <= 0 ) {
      throw std::invalid_argument( "invalid resolution for '" + i_name + "': " + i_value );
    }
    return static_cast<tsunami_lab::t_real>( l_val );
  };

  auto l_parseBoundaryCondition2d = []( std::string i_value ) {
    using BoundaryCondition2d = tsunami_lab::patches::WavePropagation2d::BoundaryCondition;

    for( char & l_ch: i_value ) {
      l_ch = static_cast<char>( std::tolower( static_cast<unsigned char>( l_ch ) ) );
      if( l_ch == '+' || l_ch == '|' ) l_ch = ',';
    }

    if( i_value == "outflow" || i_value == "none" ) {
      return BoundaryCondition2d::GhostOutflow;
    }
    if( i_value == "all" ) {
      return BoundaryCondition2d::BoundaryAll;
    }
    if( i_value == "bothx" ) {
      return BoundaryCondition2d::BoundaryBothX;
    }
    if( i_value == "bothy" || i_value == "vertical" ) {
      return BoundaryCondition2d::BoundaryBothY;
    }

    unsigned short l_mask = 0;
    std::size_t l_pos = 0;
    while( l_pos <= i_value.size() ) {
      std::size_t l_next = i_value.find( ',', l_pos );
      std::string l_token = i_value.substr( l_pos, l_next == std::string::npos ? std::string::npos : l_next - l_pos );

      if( !l_token.empty() ) {
        if( l_token == "left" ) {
          l_mask |= static_cast<unsigned short>( BoundaryCondition2d::BoundaryLeft );
        }
        else if( l_token == "right" ) {
          l_mask |= static_cast<unsigned short>( BoundaryCondition2d::BoundaryRight );
        }
        else if( l_token == "bottom" ) {
          l_mask |= static_cast<unsigned short>( BoundaryCondition2d::BoundaryBottom );
        }
        else if( l_token == "top" ) {
          l_mask |= static_cast<unsigned short>( BoundaryCondition2d::BoundaryTop );
        }
        else {
          throw std::invalid_argument( "unknown 2d boundary token: " + l_token );
        }
      }

      if( l_next == std::string::npos ) break;
      l_pos = l_next + 1;
    }

    return static_cast<BoundaryCondition2d>( l_mask );
  };

  auto l_isKnownSetupName = []( std::string i_name ) {
    for( char & l_ch: i_name ) l_ch = static_cast<char>( std::tolower( static_cast<unsigned char>( l_ch ) ) );
    return i_name == "dam_break_1d" ||
           i_name == "shock_shock_1d" ||
           i_name == "rare_rare_1d" ||
           i_name == "subcritical_1d" ||
           i_name == "supercritical_1d" ||
           i_name == "tsunami_event_1d" ||
          i_name == "circular_dam_break_2d" ||
          i_name == "artificial_tsunami_2d" ||
          i_name == "tsunami_event_2d" ||
          i_name == "tsunami2d";
  };

  // Detect resolution mode: first positional arg ends with 'm' (e.g. "250m")
  bool l_resolutionMode = (i_argc >= 2 && std::string(i_argv[1]).back() == 'm');
  int  l_minArgc        = l_resolutionMode ? 10 : 11;

  if( i_argc < l_minArgc ) {
    std::cerr << "invalid number of arguments" << std::endl;
    l_printUsage();
    return EXIT_FAILURE;
  }

  tsunami_lab::setups::Setup *l_setup = nullptr;
  tsunami_lab::patches::WavePropagation *l_waveProp = nullptr;

  // Check whether a checkpoint exists. The actual resume decision is made after
  // parsing the current command so we can compare the full input signature.
  static constexpr char const * l_checkpointPath = "solutions/checkpoint.nc";
  std::string l_solutionPath = "solutions/solution.nc";
  std::string l_inputSignature;
  bool l_defaultCheckpointExists = std::filesystem::exists( l_checkpointPath );
  bool l_checkpointReadable = false;
  bool l_useCheckpoint = false;
  tsunami_lab::io::NetCdf::CheckpointData l_checkpointData;
  if( l_defaultCheckpointExists ) {
    try {
      l_checkpointData = tsunami_lab::io::NetCdf::readCheckpoint( l_checkpointPath );
      l_checkpointReadable = true;
    }
    catch( std::exception const & i_ex ) {
      std::cout << "Failed to read checkpoint data: " << i_ex.what() << std::endl;
    }
  }

  auto l_makeSeparateSolutionPath = []() {
    auto const l_ticks = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::system_clock::now().time_since_epoch()
    ).count();
    std::filesystem::path l_dir = std::filesystem::path( "solutions" ) /
      ( "run_" + std::to_string( l_ticks ) );
    tsunami_lab::t_idx l_suffix = 0;
    while( !std::filesystem::create_directories( l_dir ) ) {
      l_dir = std::filesystem::path( "solutions" ) /
        ( "run_" + std::to_string( l_ticks ) + "_" + std::to_string( l_suffix++ ) );
    }
    return ( l_dir / "solution.nc" ).string();
  };

  try {
    // l_argBase is the index of the SOLVER_MODE argument
    int l_argBase;
    if( l_resolutionMode ) {
      tsunami_lab::t_real l_res = l_parseResolution( i_argv[1], "RESOLUTION" );
      l_domainStartX = l_parseReal( i_argv[2], "X_LOWER" );
      l_domainEndX   = l_parseReal( i_argv[3], "X_UPPER" );
      l_domainStartY = l_parseReal( i_argv[4], "Y_LOWER" );
      l_domainEndY   = l_parseReal( i_argv[5], "Y_UPPER" );
      l_k            = l_parseIdx( i_argv[6], "K" );
      l_domainSizeX  = l_domainEndX - l_domainStartX;
      l_domainSizeY  = l_domainEndY - l_domainStartY;
      if( l_domainSizeX <= 0 || l_domainSizeY <= 0 ) {
        throw std::invalid_argument( "X_UPPER must be > X_LOWER and Y_UPPER must be > Y_LOWER" );
      }
      l_nx = std::max( (tsunami_lab::t_idx)1,
                       (tsunami_lab::t_idx)std::round( l_domainSizeX / l_res ) );
      l_ny = std::max( (tsunami_lab::t_idx)1,
                       (tsunami_lab::t_idx)std::round( l_domainSizeY / l_res ) );
      if( l_ny < l_k || l_nx < l_k ) {
        throw std::invalid_argument( "K must be a positive integer." );
      }
      l_argBase = 7;
    }
    else {
      l_nx           = l_parseIdx(  i_argv[1], "NX" );
      l_ny           = l_parseIdx(  i_argv[2], "NY" );
      l_domainStartX = l_parseReal( i_argv[3], "X_LOWER" );
      l_domainEndX   = l_parseReal( i_argv[4], "X_UPPER" );
      l_domainStartY = l_parseReal( i_argv[5], "Y_LOWER" );
      l_domainEndY   = l_parseReal( i_argv[6], "Y_UPPER" );
      l_k            = l_parseIdx( i_argv[7], "K" );
      l_domainSizeX  = l_domainEndX - l_domainStartX;
      l_domainSizeY  = l_domainEndY - l_domainStartY;
      if( l_domainSizeX <= 0 || l_domainSizeY <= 0 ) {
        throw std::invalid_argument( "X_UPPER must be > X_LOWER and Y_UPPER must be > Y_LOWER" );
      }
      if( l_ny < l_k || l_nx < l_k ) {
        throw std::invalid_argument( "K must be a positive integer divisor of NX and NY" );
      }
      l_argBase = 8;
    }

    int l_solverMode = l_parseSolverMode( i_argv[l_argBase] );
    l_useFWaveSolver = (l_solverMode == 1);

    l_endTime = l_parseReal( i_argv[l_argBase + 1], "END_TIME" );
    if( l_endTime <= 0 ) {
      throw std::invalid_argument( "END_TIME must be > 0" );
    }

    int l_propArgIdx = l_argBase + 2;
    std::string l_propagationMode = i_argv[l_propArgIdx];
    for( char & l_ch: l_propagationMode ) l_ch = static_cast<char>( std::tolower( static_cast<unsigned char>( l_ch ) ) );

    int l_setupArgStart = 0;
    tsunami_lab::patches::WavePropagation1d::BoundaryCondition l_boundaryCondition =
      tsunami_lab::patches::WavePropagation1d::BoundaryCondition::GhostOutflow;
    tsunami_lab::patches::WavePropagation2d::BoundaryCondition l_boundaryCondition2d =
      tsunami_lab::patches::WavePropagation2d::BoundaryCondition::GhostOutflow;
    std::string l_boundarySignature = "0";

    if( l_propagationMode == "1d" ) {
      l_useWavePropagation1d = true;
      if( i_argc < l_propArgIdx + 3 ) {
        throw std::invalid_argument( "propagation mode '1d' requires a boundary mode and setup" );
      }

      std::string l_boundaryMode = i_argv[l_propArgIdx + 1];
      for( char & l_ch: l_boundaryMode ) l_ch = static_cast<char>( std::tolower( static_cast<unsigned char>( l_ch ) ) );
      l_boundarySignature = l_boundaryMode;

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

      l_setupArgStart = l_propArgIdx + 2;
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

      std::string l_nextArg = i_argv[l_propArgIdx + 1];
      for( char & l_ch: l_nextArg ) l_ch = static_cast<char>( std::tolower( static_cast<unsigned char>( l_ch ) ) );

      if( l_isKnownSetupName( l_nextArg ) ) {
        l_setupArgStart = l_propArgIdx + 1;
        l_boundarySignature = "0";
      }
      else {
        if( i_argc < l_propArgIdx + 3 ) {
          throw std::invalid_argument( "propagation mode '2d' boundary mode requires a setup" );
        }
        l_boundaryCondition2d = l_parseBoundaryCondition2d( l_nextArg );
        l_boundarySignature = std::to_string( static_cast<unsigned short>( l_boundaryCondition2d ) );
        l_setupArgStart = l_propArgIdx + 2;
      }

      l_waveProp = new tsunami_lab::patches::WavePropagation2d( l_nx,
                                                                 l_ny,
                                                                 l_useFWaveSolver,
                                                                 l_boundaryCondition2d );
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

    bool l_setupIs2d = (l_setupName == "circular_dam_break_2d" ||
                        l_setupName == "artificial_tsunami_2d" ||
                        l_setupName == "tsunami_event_2d" ||
                        l_setupName == "tsunami2d");
    if( l_useWavePropagation1d && l_setupIs2d ) {
      throw std::invalid_argument( "2d setup requires 2d propagation" );
    }

    std::ostringstream l_signatureStream;
    l_signatureStream << std::setprecision( std::numeric_limits<tsunami_lab::t_real>::max_digits10 );
    l_signatureStream << "argc=" << ( i_argc - 1 );
    for( int l_ar = 1; l_ar < i_argc; l_ar++ ) {
      std::string const l_arg = i_argv[l_ar];
      l_signatureStream << ";arg_len=" << l_arg.size() << ";arg=" << l_arg;
    }
    l_signatureStream << ";nx=" << l_nx
                      << ";ny=" << l_ny
                      << ";x_lower=" << l_domainStartX
                      << ";x_upper=" << l_domainEndX
                      << ";y_lower=" << l_domainStartY
                      << ";y_upper=" << l_domainEndY
                      << ";k=" << l_k
                      << ";solver_mode=" << ( l_useFWaveSolver ? 1 : 0 )
                      << ";end_time=" << l_endTime
                      << ";propagation=" << l_propagationMode
                      << ";boundary=" << l_boundarySignature
                      << ";setup=" << l_setupName;
    for( std::string const & l_arg: l_setupArgs ) {
      l_signatureStream << ";setup_arg=" << l_arg;
    }
    l_inputSignature = l_signatureStream.str();

    if( l_defaultCheckpointExists ) {
      if( l_checkpointReadable && l_checkpointData.simTime > 0 &&
          !l_checkpointData.inputSignature.empty() &&
          l_checkpointData.inputSignature == l_inputSignature ) {
        l_useCheckpoint = true;
        std::cout << "checkpoint detected — resuming from t = " << l_checkpointData.simTime << std::endl;
      }
      else {
        l_solutionPath = l_makeSeparateSolutionPath();
        std::cout << "checkpoint exists but does not match current input; writing separate output to "
                  << l_solutionPath << std::endl;
      }
    }

    if( l_useCheckpoint ) {
      l_setup = new tsunami_lab::setups::Checkpoint( l_checkpointPath );
      auto * l_cpSetup = static_cast<tsunami_lab::setups::Checkpoint *>( l_setup );
      l_nx         = l_cpSetup->getNx();
      l_ny         = l_cpSetup->getNy();
      l_dx         = l_cpSetup->getDx();
      l_dy         = l_cpSetup->getDy();
      l_k          = l_cpSetup->getK();
      l_domainStartX = l_cpSetup->getOriginX();
      l_domainStartY = l_cpSetup->getOriginY();
      l_domainSizeX  = l_dx * l_nx;
      l_domainSizeY  = l_dy * l_ny;
      l_domainEndX   = l_domainStartX + l_domainSizeX;
      l_domainEndY   = l_domainStartY + l_domainSizeY;
      l_endTime      = l_cpSetup->getEndTime();
      l_useFWaveSolver = (l_cpSetup->getSolverMode() == 1);
      l_useWavePropagation1d = (l_cpSetup->getPropagation() == "1d");
      l_setupName  = l_cpSetup->getSetup();

      delete l_waveProp;
      if( l_useWavePropagation1d ) {
        tsunami_lab::patches::WavePropagation1d::BoundaryCondition l_bc =
          tsunami_lab::patches::WavePropagation1d::BoundaryCondition::GhostOutflow;
        l_waveProp = new tsunami_lab::patches::WavePropagation1d( l_nx,
                                                                   l_useFWaveSolver,
                                                                   l_bc );
      }
      else {
        tsunami_lab::patches::WavePropagation2d::BoundaryCondition l_bc2d =
          tsunami_lab::patches::WavePropagation2d::BoundaryCondition::GhostOutflow;
        l_waveProp = new tsunami_lab::patches::WavePropagation2d( l_nx,
                                                                   l_ny,
                                                                   l_useFWaveSolver,
                                                                   l_bc2d );
      }
    }

    if( !l_useCheckpoint ) {
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
      else if( l_setupName == "artificial_tsunami_2d" ) {
        if( !l_setupArgs.empty() ) {
          throw std::invalid_argument( "artificial_tsunami_2d expects no parameters" );
        }
        l_setup = new tsunami_lab::setups::ArtificialTsunami2d();
      }
      else if( l_setupName == "tsunami_event_2d" ) {
        if( l_setupArgs.size() != 2 ) {
          throw std::invalid_argument( "tsunami_event_2d expects 2 parameters" );
        }
        l_setup = new tsunami_lab::setups::TsunamiEvent2d( l_setupArgs[0], l_setupArgs[1] );
      }
      else if( l_setupName == "tsunami2d" ) {
        if( l_setupArgs.size() == 1 ) {
          l_setup = new tsunami_lab::setups::TsunamiEvent2d( l_setupArgs[0], "" );
        }
        else if( l_setupArgs.size() == 2 ) {
          l_setup = new tsunami_lab::setups::TsunamiEvent2d( l_setupArgs[0], l_setupArgs[1] );
        }
        else {
          throw std::invalid_argument( "tsunami2d expects 1 or 2 parameters (bathymetry file and optional displacement file)" );
        }
      }
      else {
        throw std::invalid_argument( "unknown setup: " + l_setupName );
      }
    }

    // Calculate cell sizes from independent x and y domain extents
    l_dx = l_domainSizeX / l_nx;
    l_dy = l_domainSizeY / l_ny;
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

  if( l_nx < 1 || l_ny < 1 || l_domainSizeX <= 0 || l_domainSizeY <= 0 || l_endTime <= 0 ) {
      std::cerr << "invalid configuration: check cell counts, domain size, and end time" << std::endl;
      delete l_setup;
      delete l_waveProp;
      return EXIT_FAILURE;
  }

  std::cout << "runtime configuration" << std::endl;
  std::cout << "  number of cells in x-direction: " << l_nx << std::endl;
  std::cout << "  number of cells in y-direction: " << l_ny << std::endl;
  std::cout << "  domain x:                       [" << l_domainStartX << ", " << l_domainEndX << "]" << std::endl;
  std::cout << "  domain y:                       [" << l_domainStartY << ", " << l_domainEndY << "]" << std::endl;
  std::cout << "  domain size x:                  " << l_domainSizeX << std::endl;
  std::cout << "  domain size y:                  " << l_domainSizeY << std::endl;
  std::cout << "  cell size x:                    " << l_dx << std::endl;
  std::cout << "  cell size y:                    " << l_dy << std::endl;
  std::cout << "  coarse output factor k:         " << l_k << std::endl;
  std::cout << "  propagation:                    " << (l_useWavePropagation1d ? "1d" : "2d") << std::endl;
  std::cout << "  setup:                          " << l_setupName << std::endl;
  std::cout << "  solver:                         " << (l_useFWaveSolver ? "FWave" : "Roe") << std::endl;
  std::cout << "  end time:                       " << l_endTime << std::endl;
  
  tsunami_lab::io::Stations stations(10.0); // output frequency

  // maximum observed height in the setup
  tsunami_lab::t_real l_hMax = std::numeric_limits< tsunami_lab::t_real >::lowest();

  // set up solver
  for( tsunami_lab::t_idx l_cy = 0; l_cy < l_ny; l_cy++ ) {
    tsunami_lab::t_real l_y = l_domainStartY + l_cy * l_dy;

    for( tsunami_lab::t_idx l_cx = 0; l_cx < l_nx; l_cx++ ) {
      tsunami_lab::t_real l_x = l_domainStartX + l_cx * l_dx;

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
  // use minimum cell size for CFL stability condition
  tsunami_lab::t_real l_dMin = std::min( l_dx, l_dy );
  tsunami_lab::t_real l_dt = 0.25 * l_dMin / l_speedMax;

  // derive scaling for a time step
  tsunami_lab::t_real l_scaling = l_dt / l_dMin;

  // set up time and print control (endTime is set at the top)
  tsunami_lab::t_idx  l_timeStep = 0;
  tsunami_lab::t_idx  l_nOut = 0;
  tsunami_lab::t_real l_simTime = 0;
  double l_timeSteppingSeconds = 0.0;


  std::cout << "entering time loop" << std::endl;

  std::filesystem::create_directories( "solutions" );

  // When resuming from a checkpoint, keep the existing solution.nc so the
  // writer can append to it; otherwise start fresh.
  if( !l_useCheckpoint ) {
    std::filesystem::remove( l_solutionPath );
  }
  else {
    // restore simTime to the last committed checkpoint time
    l_simTime = static_cast<tsunami_lab::setups::Checkpoint *>( l_setup )->getSimTime();
    std::cout << "  resuming simulation at t = " << l_simTime << std::endl;
  }

  // create netCDF writer
  tsunami_lab::io::NetCdf l_netCdf( l_dx,
                                    l_dy,
                                    l_domainStartX,
                                    l_domainStartY,
                                    l_nx,
                                    l_ny,
                                    l_k,
                                    l_useWavePropagation1d ? 1 : l_waveProp->getStride(),
                                    l_simTime,
                                    l_endTime,
                                    l_useFWaveSolver ? 1 : 0,
                                    l_useWavePropagation1d ? "1d" : "2d",
                                    l_setupName,
                                    l_inputSignature,
                                    l_waveProp->getHeight(),
                                    l_waveProp->getBathymetry(),
                                    l_waveProp->getMomentumX(),
                                    l_useWavePropagation1d ? nullptr : l_waveProp->getMomentumY(),
                                    l_solutionPath );

  // iterate over time
  while( l_simTime < l_endTime ){
    if( l_timeStep % 25 == 0 ) {
      std::cout << "  simulation time / #time steps: "
                << l_simTime << " / " << l_timeStep << std::endl;

      std::cout << "  appending wave field to " << l_solutionPath << std::endl;

      // write time step to netCDF file
      l_netCdf.writeTimeStep( l_simTime );

      l_nOut++;
    }
    // write station data to CSV files
    stations.writeToCSV(  l_simTime,
                          l_dx,
                          l_dy,
                          l_domainStartX,
                          l_domainStartY,
                          l_waveProp
                        );

    auto const l_stepStart = std::chrono::steady_clock::now();
    l_waveProp->setGhostOutflow();
    l_waveProp->timeStep( l_scaling );
    auto const l_stepEnd = std::chrono::steady_clock::now();

    l_timeSteppingSeconds += std::chrono::duration<double>(
      l_stepEnd - l_stepStart
    ).count();

    l_timeStep++;
    l_simTime += l_dt;
  }

  double const l_timePerCellIteration = l_timeSteppingSeconds /
    ( static_cast<double>( l_nx ) * static_cast<double>( l_ny ) *
      static_cast<double>( l_timeStep ) );

  std::cout << "time stepping seconds: " << l_timeSteppingSeconds << std::endl;
  std::cout << "time steps: " << l_timeStep << std::endl;
  std::cout << "time per cell and iteration: "
            << l_timePerCellIteration << std::endl;

  std::cout << "finished time loop" << std::endl;

  // free memory
  std::cout << "freeing memory" << std::endl;
  delete l_setup;
  delete l_waveProp;

  std::cout << "finished, exiting" << std::endl;
  return EXIT_SUCCESS;
}
