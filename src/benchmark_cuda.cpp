/**
 * @section DESCRIPTION
 * CUDA runtime benchmark for the two-dimensional tsunami setup.
 * Benchmarks the Tohoku 2011 event using the GPU-accelerated solver.
 **/
#include "cuda/WavePropagation2d_cuda.h"
#include "patches/WavePropagation2d.h"
#include "setups/TsunamiEvent2d/TsunamiEvent2d.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
  struct RunResult {
    double timeSteppingSeconds = 0.0;
    tsunami_lab::t_idx timeSteps = 0;
  };

  struct VerificationResult {
    double maxErrorH  = 0.0;
    double maxErrorHu = 0.0;
    double maxErrorHv = 0.0;
    tsunami_lab::t_idx timeSteps = 0;
  };

  unsigned int parsePositiveInt( char const * i_value,
                                 std::string const & i_name ) {
    char * l_endPtr = nullptr;
    unsigned long l_val = std::strtoul( i_value, &l_endPtr, 10 );
    if( l_endPtr == i_value || *l_endPtr != '\0' || l_val < 1 ) {
      throw std::invalid_argument( i_name + " must be a positive integer" );
    }
    return static_cast<unsigned int>( l_val );
  }

  tsunami_lab::t_real parsePositiveReal( char const * i_value,
                                         std::string const & i_name ) {
    char * l_endPtr = nullptr;
    double l_val = std::strtod( i_value, &l_endPtr );
    if( l_endPtr == i_value || *l_endPtr != '\0' || l_val <= 0 ) {
      throw std::invalid_argument( i_name + " must be a positive real value" );
    }
    return static_cast<tsunami_lab::t_real>( l_val );
  }

  std::string resolvePath( char const * i_argv0,
                           std::string const & i_path ) {
    namespace fs = std::filesystem;

    if( fs::exists( i_path ) ) return i_path;

    try {
      fs::path l_root = fs::canonical( i_argv0 ).parent_path().parent_path();
      fs::path l_path = l_root / i_path;
      if( fs::exists( l_path ) ) return l_path.string();
    }
    catch( fs::filesystem_error const & ) {}

    return i_path;
  }

  RunResult runBenchmark( tsunami_lab::t_idx i_nx,
                          tsunami_lab::t_idx i_ny,
                          tsunami_lab::t_real i_endTime,
                          std::string const & i_bathymetryPath,
                          std::string const & i_displacementPath ) {
    std::unique_ptr<tsunami_lab::setups::TsunamiEvent2d> l_setup(
      new tsunami_lab::setups::TsunamiEvent2d( i_bathymetryPath,
                                               i_displacementPath )
    );

    tsunami_lab::t_real l_domainStartX = -200000.0;
    tsunami_lab::t_real l_domainEndX   = 2500000.0;
    tsunami_lab::t_real l_domainStartY = -750000.0;
    tsunami_lab::t_real l_domainEndY   =  750000.0;
    tsunami_lab::t_real l_domainSizeX  = l_domainEndX - l_domainStartX;
    tsunami_lab::t_real l_domainSizeY  = l_domainEndY - l_domainStartY;
    if( l_domainSizeX <= 0 || l_domainSizeY <= 0 ) {
      throw std::runtime_error( "invalid hard-coded benchmark domain" );
    }

    tsunami_lab::t_real l_dx = l_domainSizeX / i_nx;
    tsunami_lab::t_real l_dy = l_domainSizeY / i_ny;

    tsunami_lab::t_real l_hMax = std::numeric_limits<tsunami_lab::t_real>::lowest();

    std::cout << "Starting setup..." << std::endl;

    // Use a temporary CPU solver to stage the initial state for copyToGpu
    std::unique_ptr<tsunami_lab::patches::WavePropagation2d> l_cpuInit(
      new tsunami_lab::patches::WavePropagation2d( i_nx, i_ny, true )
    );

    for( tsunami_lab::t_idx l_cy = 0; l_cy < i_ny; l_cy++ ) {
      tsunami_lab::t_real l_y = l_domainStartY + l_cy * l_dy;

      for( tsunami_lab::t_idx l_cx = 0; l_cx < i_nx; l_cx++ ) {
        tsunami_lab::t_real l_x = l_domainStartX + l_cx * l_dx;
        tsunami_lab::t_real l_h  = l_setup->getHeight( l_x, l_y );
        tsunami_lab::t_real l_hu = l_setup->getMomentumX( l_x, l_y );
        tsunami_lab::t_real l_hv = l_setup->getMomentumY( l_x, l_y );
        tsunami_lab::t_real l_b  = l_setup->getBathymetry( l_x, l_y );

        l_hMax = std::max( l_h, l_hMax );

        l_cpuInit->setHeight( l_cx, l_cy, l_h );
        l_cpuInit->setMomentumX( l_cx, l_cy, l_hu );
        l_cpuInit->setMomentumY( l_cx, l_cy, l_hv );
        l_cpuInit->setBathymetry( l_cx, l_cy, l_b );
      }
    }

    // Create GPU solver and upload initial state
    std::unique_ptr<tsunami_lab::patches::cuda::WavePropagation2dCuda> l_wavePropGpu(
      new tsunami_lab::patches::cuda::WavePropagation2dCuda( i_nx, i_ny, true )
    );
    l_wavePropGpu->copyToGpu( l_cpuInit->getHeight(),
                              l_cpuInit->getMomentumX(),
                              l_cpuInit->getMomentumY(),
                              l_cpuInit->getBathymetry() );
    l_cpuInit.reset(); // free CPU staging memory

    tsunami_lab::t_real l_speedMax = std::sqrt(
      static_cast<tsunami_lab::t_real>( 9.81 ) *
      std::max( l_hMax, static_cast<tsunami_lab::t_real>( 0 ) )
    );
    if( l_speedMax <= 0 ) {
      throw std::runtime_error( "maximum wave speed is zero" );
    }

    tsunami_lab::t_real l_dMin = std::min( l_dx, l_dy );
    tsunami_lab::t_real l_dt = static_cast<tsunami_lab::t_real>( 0.25 ) *
      l_dMin / l_speedMax;
    tsunami_lab::t_real l_scaling = l_dt / l_dMin;
    tsunami_lab::t_idx l_estimatedTimeSteps = static_cast<tsunami_lab::t_idx>(
      std::ceil( static_cast<double>( i_endTime ) / static_cast<double>( l_dt ) )
    );
    tsunami_lab::t_idx l_progressInterval = std::max(
      static_cast<tsunami_lab::t_idx>( 1 ),
      l_estimatedTimeSteps / static_cast<tsunami_lab::t_idx>( 20 )
    );

    RunResult l_result;
    tsunami_lab::t_real l_simTime = 0;

    std::cout << "Setup complete. Entering time loop..." << std::endl;
    std::cout << "  domain x:           [" << l_domainStartX << ", " << l_domainEndX << "]" << std::endl;
    std::cout << "  domain y:           [" << l_domainStartY << ", " << l_domainEndY << "]" << std::endl;
    std::cout << "  cell size x:        " << l_dx << std::endl;
    std::cout << "  cell size y:        " << l_dy << std::endl;
    std::cout << "  time step size:     " << l_dt << std::endl;
    std::cout << "  estimated steps:    " << l_estimatedTimeSteps << std::endl;

    while( l_simTime < i_endTime ) {
      auto const l_stepStart = std::chrono::steady_clock::now();

      l_wavePropGpu->setGhostOutflow();
      l_wavePropGpu->computeStep( l_scaling );
      l_wavePropGpu->swapBuffers();

      auto const l_stepEnd = std::chrono::steady_clock::now();

      l_result.timeSteppingSeconds += std::chrono::duration<double>(
        l_stepEnd - l_stepStart
      ).count();

      l_result.timeSteps++;
      l_simTime += l_dt;

      if( l_result.timeSteps % l_progressInterval == 0 || l_simTime >= i_endTime ) {
        std::cout << "  time loop progress: "
                  << l_simTime << " / " << i_endTime
                  << " (#time steps: " << l_result.timeSteps << ")" << std::endl;
      }
    }

    return l_result;
  }

  VerificationResult runVerification( tsunami_lab::t_idx i_nx,
                                       tsunami_lab::t_idx i_ny,
                                       tsunami_lab::t_real i_endTime,
                                       std::string const & i_bathymetryPath,
                                       std::string const & i_displacementPath ) {
    std::unique_ptr<tsunami_lab::setups::TsunamiEvent2d> l_setup(
      new tsunami_lab::setups::TsunamiEvent2d( i_bathymetryPath,
                                               i_displacementPath )
    );

    tsunami_lab::t_real l_domainStartX = -200000.0;
    tsunami_lab::t_real l_domainEndX   = 2500000.0;
    tsunami_lab::t_real l_domainStartY = -750000.0;
    tsunami_lab::t_real l_domainEndY   =  750000.0;
    tsunami_lab::t_real l_domainSizeX  = l_domainEndX - l_domainStartX;
    tsunami_lab::t_real l_domainSizeY  = l_domainEndY - l_domainStartY;

    tsunami_lab::t_real l_dx = l_domainSizeX / i_nx;
    tsunami_lab::t_real l_dy = l_domainSizeY / i_ny;

    tsunami_lab::t_real l_hMax = std::numeric_limits<tsunami_lab::t_real>::lowest();

    std::cout << "Setting up CPU and GPU solvers for verification..." << std::endl;

    std::unique_ptr<tsunami_lab::patches::WavePropagation2d> l_cpuSolver(
      new tsunami_lab::patches::WavePropagation2d( i_nx, i_ny, true )
    );

    for( tsunami_lab::t_idx l_cy = 0; l_cy < i_ny; l_cy++ ) {
      tsunami_lab::t_real l_y = l_domainStartY + l_cy * l_dy;
      for( tsunami_lab::t_idx l_cx = 0; l_cx < i_nx; l_cx++ ) {
        tsunami_lab::t_real l_x = l_domainStartX + l_cx * l_dx;
        tsunami_lab::t_real l_h  = l_setup->getHeight( l_x, l_y );
        tsunami_lab::t_real l_hu = l_setup->getMomentumX( l_x, l_y );
        tsunami_lab::t_real l_hv = l_setup->getMomentumY( l_x, l_y );
        tsunami_lab::t_real l_b  = l_setup->getBathymetry( l_x, l_y );
        l_hMax = std::max( l_h, l_hMax );
        l_cpuSolver->setHeight( l_cx, l_cy, l_h );
        l_cpuSolver->setMomentumX( l_cx, l_cy, l_hu );
        l_cpuSolver->setMomentumY( l_cx, l_cy, l_hv );
        l_cpuSolver->setBathymetry( l_cx, l_cy, l_b );
      }
    }

    std::unique_ptr<tsunami_lab::patches::cuda::WavePropagation2dCuda> l_gpuSolver(
      new tsunami_lab::patches::cuda::WavePropagation2dCuda( i_nx, i_ny, true )
    );
    // copyToGpu accepts interior pointers (stride = nCellsX + 2)
    l_gpuSolver->copyToGpu( l_cpuSolver->getHeight(),
                            l_cpuSolver->getMomentumX(),
                            l_cpuSolver->getMomentumY(),
                            l_cpuSolver->getBathymetry() );

    tsunami_lab::t_real l_speedMax = std::sqrt(
      static_cast<tsunami_lab::t_real>( 9.81 ) *
      std::max( l_hMax, static_cast<tsunami_lab::t_real>( 0 ) )
    );
    if( l_speedMax <= 0 ) {
      throw std::runtime_error( "maximum wave speed is zero" );
    }
    tsunami_lab::t_real l_dMin    = std::min( l_dx, l_dy );
    tsunami_lab::t_real l_dt      = static_cast<tsunami_lab::t_real>( 0.25 ) * l_dMin / l_speedMax;
    tsunami_lab::t_real l_scaling = l_dt / l_dMin;

    std::cout << "Running CPU + GPU in lock-step..." << std::endl;

    VerificationResult l_result;
    tsunami_lab::t_real l_simTime = 0;
    while( l_simTime < i_endTime ) {
      l_gpuSolver->setGhostOutflow();
      l_gpuSolver->computeStep( l_scaling );
      l_gpuSolver->swapBuffers();

      l_cpuSolver->setGhostOutflow();
      l_cpuSolver->timeStep( l_scaling );

      l_result.timeSteps++;
      l_simTime += l_dt;
    }

    // Copy final GPU state back to host (full ghost-celled array)
    tsunami_lab::t_idx const l_stride = i_nx + 2;
    tsunami_lab::t_idx const l_nFull  = ( i_ny + 2 ) * l_stride;
    std::vector<tsunami_lab::t_real> l_hGpu( l_nFull, 0 );
    std::vector<tsunami_lab::t_real> l_huGpu( l_nFull, 0 );
    std::vector<tsunami_lab::t_real> l_hvGpu( l_nFull, 0 );
    l_gpuSolver->copyToHost( l_hGpu.data(), l_huGpu.data(), l_hvGpu.data() );

    // CPU interior pointer has pitch l_stride; interior cell (cx,cy) is at cy*l_stride+cx
    tsunami_lab::t_real const* l_hCpu  = l_cpuSolver->getHeight();
    tsunami_lab::t_real const* l_huCpu = l_cpuSolver->getMomentumX();
    tsunami_lab::t_real const* l_hvCpu = l_cpuSolver->getMomentumY();

    for( tsunami_lab::t_idx l_cy = 0; l_cy < i_ny; l_cy++ ) {
      for( tsunami_lab::t_idx l_cx = 0; l_cx < i_nx; l_cx++ ) {
        tsunami_lab::t_idx l_gpuIdx = ( l_cy + 1 ) * l_stride + ( l_cx + 1 );
        tsunami_lab::t_idx l_cpuIdx =   l_cy        * l_stride +   l_cx;
        l_result.maxErrorH  = std::max( l_result.maxErrorH,
            static_cast<double>( std::fabs( l_hGpu[l_gpuIdx]  - l_hCpu[l_cpuIdx]  ) ) );
        l_result.maxErrorHu = std::max( l_result.maxErrorHu,
            static_cast<double>( std::fabs( l_huGpu[l_gpuIdx] - l_huCpu[l_cpuIdx] ) ) );
        l_result.maxErrorHv = std::max( l_result.maxErrorHv,
            static_cast<double>( std::fabs( l_hvGpu[l_gpuIdx] - l_hvCpu[l_cpuIdx] ) ) );
      }
    }

    return l_result;
  }

  void printUsage() {
    std::cerr << "usage:" << std::endl;
    std::cerr << "  ./build/benchmark_cuda RUNS [END_TIME] [BATHY_NC] [DISPL_NC] [--verify]" << std::endl;
    std::cerr << "" << std::endl;
    std::cerr << "RUNS:     number of repeated benchmark runs" << std::endl;
    std::cerr << "END_TIME: optional simulation end time in seconds, default 10800" << std::endl;
    std::cerr << "BATHY_NC: optional path to bathymetry NetCDF file" << std::endl;
    std::cerr << "DISPL_NC: optional path to displacement NetCDF file" << std::endl;
    std::cerr << "--verify: after GPU runs, also run once on CPU and compare final state" << std::endl;
  }
}

int main( int i_argc,
          char * i_argv[] ) {
  std::cout << "####################################" << std::endl;
  std::cout << "### Tsunami Lab                  ###" << std::endl;
  std::cout << "###                              ###" << std::endl;
  std::cout << "### https://scalable.uni-jena.de ###" << std::endl;
  std::cout << "####################################" << std::endl;

  // Scan for --verify flag and remove it from the positional argument list
  bool l_verify = false;
  std::vector<char*> l_args;
  l_args.reserve( i_argc );
  for( int l_i = 0; l_i < i_argc; l_i++ ) {
    if( std::string( i_argv[l_i] ) == "--verify" ) {
      l_verify = true;
    } else {
      l_args.push_back( i_argv[l_i] );
    }
  }
  int    l_argc = static_cast<int>( l_args.size() );
  char** l_argv = l_args.data();

  if( l_argc < 2 || l_argc > 5 ) {
    printUsage();
    return EXIT_FAILURE;
  }

  try {
    unsigned int l_runs = parsePositiveInt( l_argv[1], "RUNS" );

    tsunami_lab::t_idx l_nx = 2160;
    tsunami_lab::t_idx l_ny = 1200;

    tsunami_lab::t_real l_endTime = 10800.0;
    if( l_argc >= 3 ) {
      l_endTime = parsePositiveReal( l_argv[2], "END_TIME" );
    }

    std::string l_bathymetryPath = resolvePath(
      l_argv[0],
      l_argc >= 4 ? l_argv[3] : "data/output/tohoku_gebco20_usgs_250m_bath.nc"
    );
    std::string l_displacementPath = resolvePath(
      l_argv[0],
      l_argc >= 5 ? l_argv[4] : "data/output/tohoku_gebco20_usgs_250m_displ.nc"
    );

    std::cout << "benchmark_cuda configuration" << std::endl;
    std::cout << "  number of cells in x-direction: " << l_nx << std::endl;
    std::cout << "  number of cells in y-direction: " << l_ny << std::endl;
    std::cout << "  runs:                           " << l_runs << std::endl;
    std::cout << "  end time:                       " << l_endTime << std::endl;
    std::cout << "  bathymetry:                     " << l_bathymetryPath << std::endl;
    std::cout << "  displacement:                   " << l_displacementPath << std::endl;
    std::cout << "  verify (CPU vs GPU):            " << (l_verify ? "yes" : "no") << std::endl;

    if( !tsunami_lab::patches::cuda::WavePropagation2dCuda::initialize( 0 ) ) {
      throw std::runtime_error( "CUDA initialization failed" );
    }

    std::vector<RunResult> l_results( l_runs );
    double l_totalTimeSteppingSeconds = 0.0;
    tsunami_lab::t_idx l_totalTimeSteps = 0;

    for( unsigned int l_run = 0; l_run < l_runs; l_run++ ) {
      std::cout << "starting run " << l_run << std::endl;
      l_results[l_run] = runBenchmark( l_nx,
                                       l_ny,
                                       l_endTime,
                                       l_bathymetryPath,
                                       l_displacementPath );

      double const l_timePerCellIteration = l_results[l_run].timeSteppingSeconds /
        ( static_cast<double>( l_nx ) * static_cast<double>( l_ny ) *
          static_cast<double>( l_results[l_run].timeSteps ) );
      double const l_timePerCellIterationNs = l_timePerCellIteration * 1.0e9;

      std::cout << "time stepping seconds: " << l_results[l_run].timeSteppingSeconds << std::endl;
      std::cout << "time steps: " << l_results[l_run].timeSteps << std::endl;
      std::cout << "time per cell and iteration: "
                << l_timePerCellIteration << std::endl;
      std::cout << "time per cell and iteration in ns: "
                << l_timePerCellIterationNs << std::endl;
      std::cout << "done with run " << l_run << std::endl;

      l_totalTimeSteppingSeconds += l_results[l_run].timeSteppingSeconds;
      l_totalTimeSteps += l_results[l_run].timeSteps;
    }

    double const l_averageTimeSteppingSeconds = l_totalTimeSteppingSeconds /
      static_cast<double>( l_runs );
    double const l_averageTimeSteps = static_cast<double>( l_totalTimeSteps ) /
      static_cast<double>( l_runs );
    double const l_averageTimePerCellIteration = l_averageTimeSteppingSeconds /
      ( static_cast<double>( l_nx ) * static_cast<double>( l_ny ) *
        l_averageTimeSteps );
    double const l_averageTimePerCellIterationNs = l_averageTimePerCellIteration * 1.0e9;

    std::cout << "benchmark average" << std::endl;
    std::cout << "time stepping seconds: " << l_averageTimeSteppingSeconds << std::endl;
    std::cout << "time steps: " << l_averageTimeSteps << std::endl;
    std::cout << "time per cell and iteration: "
              << l_averageTimePerCellIteration << std::endl;
    std::cout << "time per cell and iteration in ns: "
              << l_averageTimePerCellIterationNs << std::endl;

    if( l_verify ) {
      std::cout << "starting CPU vs GPU verification run" << std::endl;
      VerificationResult l_vr = runVerification( l_nx,
                                                  l_ny,
                                                  l_endTime,
                                                  l_bathymetryPath,
                                                  l_displacementPath );
      std::cout << "verification complete" << std::endl;
      std::cout << "  time steps:        " << l_vr.timeSteps << std::endl;
      std::cout << "  max error h:       " << l_vr.maxErrorH  << std::endl;
      std::cout << "  max error hu:      " << l_vr.maxErrorHu << std::endl;
      std::cout << "  max error hv:      " << l_vr.maxErrorHv << std::endl;
    }

    tsunami_lab::patches::cuda::WavePropagation2dCuda::finalize();

    std::cout << "finished, exiting" << std::endl;
  }
  catch( std::exception const & i_ex ) {
    std::cerr << "benchmark_cuda failed: " << i_ex.what() << std::endl;
    printUsage();
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
