/**
 * @section DESCRIPTION
 * Grace system benchmark for the two-dimensional tsunami setup.
 * Benchmarks the Tohoku 2011 M 9.1 event at 250m resolution with file I/O every 100 time steps.
 **/
#include "patches/WavePropagation2d.h"
#include "setups/TsunamiEvent2d/TsunamiEvent2d.h"
#include "io/NetCdf.h"
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
#include <omp.h>

namespace {
  struct BenchmarkResult {
    double totalSeconds = 0.0;
    double computationSeconds = 0.0;
    double ioSeconds = 0.0;
    tsunami_lab::t_idx timeSteps = 0;
    tsunami_lab::t_idx cellsX = 0;
    tsunami_lab::t_idx cellsY = 0;
  };

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

  BenchmarkResult runGraceBenchmark(
    std::string const & i_bathymetryPath,
    std::string const & i_displacementPath,
    std::string const & i_outputPath ) {
    
    // Tohoku 2011 M 9.1 at 250m resolution
    static constexpr tsunami_lab::t_real l_domainStartX = -200000.0;
    static constexpr tsunami_lab::t_real l_domainEndX = 2500000.0;
    static constexpr tsunami_lab::t_real l_domainStartY = -750000.0;
    static constexpr tsunami_lab::t_real l_domainEndY = 750000.0;
    static constexpr tsunami_lab::t_real l_dx = 250.0;
    static constexpr tsunami_lab::t_real l_dy = 250.0;
    static constexpr tsunami_lab::t_idx l_nx = 10800;  // (2500000 - (-200000)) / 250 = 10800
    static constexpr tsunami_lab::t_idx l_ny = 6000;   // (750000 - (-750000)) / 250 = 6000
    static constexpr tsunami_lab::t_idx l_numTimeSteps = 10000;
    static constexpr tsunami_lab::t_idx l_ioInterval = 100;

    std::unique_ptr<tsunami_lab::setups::TsunamiEvent2d> l_setup(
      new tsunami_lab::setups::TsunamiEvent2d( i_bathymetryPath,
                                               i_displacementPath )
    );

    std::unique_ptr<tsunami_lab::patches::WavePropagation> l_waveProp(
      new tsunami_lab::patches::WavePropagation2d( l_nx,
                                                   l_ny,
                                                   true )
    );

    tsunami_lab::t_real l_hMax = std::numeric_limits<tsunami_lab::t_real>::lowest();

    std::cout << "Setting up Tohoku 2011 M 9.1 event at 250m resolution..." << std::endl;
    std::cout << "  Grid: " << l_nx << " x " << l_ny << " cells" << std::endl;
    std::cout << "  Domain: [" << l_domainStartX << ", " << l_domainEndX << "] x "
              << "[" << l_domainStartY << ", " << l_domainEndY << "]" << std::endl;
    std::cout << "  Cell size: " << l_dx << " x " << l_dy << " meters" << std::endl;
    std::cout << "  Total cells: " << (l_nx * l_ny) << std::endl;

    auto const l_setupStart = std::chrono::steady_clock::now();

    for( tsunami_lab::t_idx l_cy = 0; l_cy < l_ny; l_cy++ ) {
      tsunami_lab::t_real l_y = l_domainStartY + l_cy * l_dy;

      for( tsunami_lab::t_idx l_cx = 0; l_cx < l_nx; l_cx++ ) {
        tsunami_lab::t_real l_x = l_domainStartX + l_cx * l_dx;
        tsunami_lab::t_real l_h = l_setup->getHeight( l_x, l_y );
        tsunami_lab::t_real l_hu = l_setup->getMomentumX( l_x, l_y );
        tsunami_lab::t_real l_hv = l_setup->getMomentumY( l_x, l_y );
        tsunami_lab::t_real l_b = l_setup->getBathymetry( l_x, l_y );

        l_hMax = std::max( l_h, l_hMax );

        l_waveProp->setHeight( l_cx, l_cy, l_h );
        l_waveProp->setMomentumX( l_cx, l_cy, l_hu );
        l_waveProp->setMomentumY( l_cx, l_cy, l_hv );
        l_waveProp->setBathymetry( l_cx, l_cy, l_b );
      }
    }

    auto const l_setupEnd = std::chrono::steady_clock::now();
    double l_setupTime = std::chrono::duration<double>(
      l_setupEnd - l_setupStart
    ).count();

    std::cout << "Setup complete in " << l_setupTime << " seconds." << std::endl;

    tsunami_lab::t_real l_speedMax = std::sqrt(
      static_cast<tsunami_lab::t_real>( 9.81 ) *
      std::max( l_hMax, static_cast<tsunami_lab::t_real>( 0 ) )
    );
    if( l_speedMax <= 0 ) {
      throw std::runtime_error( "maximum wave speed is zero" );
    }

    tsunami_lab::t_real l_dMin = std::min( l_dx, l_dy );
    tsunami_lab::t_real l_dt = static_cast<tsunami_lab::t_real>( 0.5 ) *
      l_dMin / l_speedMax;
    tsunami_lab::t_real l_scaling = l_dt / l_dMin;

    std::cout << "Time step size: " << l_dt << " seconds" << std::endl;
    std::cout << "Maximum wave speed: " << l_speedMax << " m/s" << std::endl;

    // Initialize NetCDF output
    std::unique_ptr<tsunami_lab::io::NetCdf> l_io(
      new tsunami_lab::io::NetCdf(
        l_dx,
        l_dy,
        l_domainStartX,
        l_domainStartY,
        l_nx,
        l_ny,
        1,  // k=1, no averaging
        l_nx,  // stride
        0.0,  // initial simulation time
        l_numTimeSteps * l_dt,  // end time
        1,  // solver mode: 1 for F-Wave
        "2d",  // propagation mode
        "tsunami2d",  // setup name
        "Grace HPC Benchmark",  // input signature
        l_waveProp->getHeight(),
        l_waveProp->getBathymetry(),
        l_waveProp->getMomentumX(),
        l_waveProp->getMomentumY(),
        i_outputPath
      )
    );

    BenchmarkResult l_result;
    l_result.cellsX = l_nx;
    l_result.cellsY = l_ny;
    l_result.timeSteps = l_numTimeSteps;

    auto const l_totalStart = std::chrono::steady_clock::now();

    std::cout << "\nEntering time loop: " << l_numTimeSteps << " time steps..." << std::endl;

    for( tsunami_lab::t_idx l_step = 0; l_step < l_numTimeSteps; l_step++ ) {
      auto const l_stepStart = std::chrono::steady_clock::now();

      l_waveProp->setGhostOutflow();
      l_waveProp->timeStep( l_scaling );

      auto const l_stepEnd = std::chrono::steady_clock::now();
      l_result.computationSeconds += std::chrono::duration<double>(
        l_stepEnd - l_stepStart
      ).count();

      // Output every 100 time steps
      if( (l_step + 1) % l_ioInterval == 0 || l_step == l_numTimeSteps - 1 ) {
        auto const l_ioStart = std::chrono::steady_clock::now();

        tsunami_lab::t_real l_simTime = (l_step + 1) * l_dt;
        l_io->writeTimeStep( l_simTime );

        auto const l_ioEnd = std::chrono::steady_clock::now();
        l_result.ioSeconds += std::chrono::duration<double>(
          l_ioEnd - l_ioStart
        ).count();

        std::cout << "  Step " << (l_step + 1) << " / " << l_numTimeSteps
                  << " (t = " << l_simTime << " s)" << std::endl;
      }
    }

    auto const l_totalEnd = std::chrono::steady_clock::now();
    l_result.totalSeconds = std::chrono::duration<double>(
      l_totalEnd - l_totalStart
    ).count();

    return l_result;
  }

  void printUsage() {
    std::cerr << "usage:" << std::endl;
    std::cerr << "  ./build/benchmark_grace [BATHY_NC] [DISPL_NC] [OUTPUT_NC]" << std::endl;
    std::cerr << "" << std::endl;
    std::cerr << "BATHY_NC:   optional bathymetry netCDF file path" << std::endl;
    std::cerr << "DISPL_NC:   optional displacement netCDF file path" << std::endl;
    std::cerr << "OUTPUT_NC:  optional output netCDF file path" << std::endl;
  }
}

int main( int i_argc,
          char * i_argv[] ) {
  std::cout << "####################################" << std::endl;
  std::cout << "### Tsunami Lab - Grace Benchmark ###" << std::endl;
  std::cout << "###                              ###" << std::endl;
  std::cout << "### https://scalable.uni-jena.de ###" << std::endl;
  std::cout << "####################################" << std::endl;
  std::cout << std::endl;

  if( i_argc > 4 ) {
    printUsage();
    return EXIT_FAILURE;
  }

  try {
    std::string l_bathymetryPath = resolvePath(
      i_argv[0],
      i_argc >= 2 ? i_argv[1] : "data/output/tohoku_gebco20_usgs_250m_bath.nc"
    );
    std::string l_displacementPath = resolvePath(
      i_argv[0],
      i_argc >= 3 ? i_argv[2] : "data/output/tohoku_gebco20_usgs_250m_displ.nc"
    );
    std::string l_outputPath = resolvePath(
      i_argv[0],
      i_argc >= 4 ? i_argv[3] : "solutions/grace_benchmark.nc"
    );

    std::cout << "Benchmark Configuration:" << std::endl;
    std::cout << "  Event: Tohoku 2011 M 9.1" << std::endl;
    std::cout << "  Resolution: 250 meters" << std::endl;
    std::cout << "  Time steps: 10000" << std::endl;
    std::cout << "  File I/O interval: Every 100 time steps" << std::endl;
    std::cout << "  Bathymetry: " << l_bathymetryPath << std::endl;
    std::cout << "  Displacement: " << l_displacementPath << std::endl;
    std::cout << "  Output: " << l_outputPath << std::endl;
    std::cout << std::endl;

    BenchmarkResult l_result = runGraceBenchmark(
      l_bathymetryPath,
      l_displacementPath,
      l_outputPath
    );

    // Calculate metrics
    double l_totalCellUpdates = static_cast<double>( l_result.cellsX ) *
                                static_cast<double>( l_result.cellsY ) *
                                static_cast<double>( l_result.timeSteps );
    double l_cellUpdatesPerSecond = l_totalCellUpdates / l_result.totalSeconds;
    double l_megaCellUpdatesPerSecond = l_cellUpdatesPerSecond / 1e6;

    std::cout << std::endl;
    std::cout << "=====================================" << std::endl;
    std::cout << "GRACE BENCHMARK RESULTS" << std::endl;
    std::cout << "=====================================" << std::endl;
    std::cout << std::endl;
    std::cout << "Team: Waves of Change Collective" << std::endl;
    std::cout << std::endl;
    std::cout << "Performance Metrics:" << std::endl;
    std::cout << "  Total Elapsed Time:         " << l_result.totalSeconds << " seconds" << std::endl;
    std::cout << "  Computation Time:           " << l_result.computationSeconds << " seconds ("
              << (100.0 * l_result.computationSeconds / l_result.totalSeconds) << "%)" << std::endl;
    std::cout << "  I/O Time:                   " << l_result.ioSeconds << " seconds ("
              << (100.0 * l_result.ioSeconds / l_result.totalSeconds) << "%)" << std::endl;
    std::cout << std::endl;
    std::cout << "  Grid Resolution:            " << l_result.cellsX << " x " << l_result.cellsY << " cells" << std::endl;
    std::cout << "  Total Cells:                " << l_totalCellUpdates / l_result.timeSteps << std::endl;
    std::cout << "  Time Steps:                 " << l_result.timeSteps << std::endl;
    std::cout << "  Total Cell Updates:         " << static_cast<long long>( l_totalCellUpdates ) << std::endl;
    std::cout << std::endl;
    std::cout << "  Cell Updates per Second:    " << static_cast<long long>( l_cellUpdatesPerSecond ) << std::endl;
    std::cout << "  Mega Cell Updates/Second:   " << l_megaCellUpdatesPerSecond << " MCUps" << std::endl;
    std::cout << std::endl;
    std::cout << "=====================================" << std::endl;
    std::cout << "Benchmark completed successfully!" << std::endl;
    std::cout << "=====================================" << std::endl;
  }
  catch( std::exception const & i_ex ) {
    std::cerr << "benchmark failed: " << i_ex.what() << std::endl;
    printUsage();
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
