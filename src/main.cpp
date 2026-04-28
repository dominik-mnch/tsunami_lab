/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Entry-point for simulations.
 **/
#include "patches/WavePropagation1d.h"
#include "setups/DamBreak1d/DamBreak1d.h"
#include "setups/RareRare1d/RareRare1d.h"
#include "setups/ShockShock1d/ShockShock1d.h"
#include "setups/Subcritical1d/Subcritical1d.h"
#include "setups/Supercritical1d/Supercritical1d.h"
#include "io/Csv.h"
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <limits>
#include <filesystem>

int main( int i_argc, char *i_argv[] ) {
  // number of cells in x- and y-direction
  tsunami_lab::t_idx l_nx = 0;
  tsunami_lab::t_idx l_ny = 1;

  // set cell size, domain, and end time
  tsunami_lab::t_real l_dxy = 1;
  tsunami_lab::t_real l_domainSize = 50000.0;
  tsunami_lab::t_real l_endTime = 1.25;
  bool l_useFWaveSolver = true;

  std::cout << "####################################" << std::endl;
  std::cout << "### Tsunami Lab                  ###" << std::endl;
  std::cout << "###                              ###" << std::endl;
  std::cout << "### https://scalable.uni-jena.de ###" << std::endl;
  std::cout << "####################################" << std::endl;

  // Expecting 5 arguments
  if( i_argc != 5 ) {
    std::cerr << "invalid number of arguments, usage:" << std::endl;
    std::cerr << "  ./build/tsunami_lab N_CELLS_X DOMAIN_SIZE SOLVER_MODE END_TIME" << std::endl;
    std::cerr << "where:" << std::endl;
    std::cerr << "  N_CELLS_X   is the number of cells in x-direction." << std::endl;
    std::cerr << "  DOMAIN_SIZE is the total length of the domain." << std::endl;
    std::cerr << "  SOLVER_MODE is 1 for FWaveSolver or 0 for RoeSolver." << std::endl;
    std::cerr << "  END_TIME    is the simulation time in seconds." << std::endl;
    return EXIT_FAILURE;
  }
  else {
    l_nx = atoi( i_argv[1] );
    l_domainSize = atof( i_argv[2] );
    l_useFWaveSolver = (atoi( i_argv[3] ) == 1);
    l_endTime = atof( i_argv[4] );

    if( l_nx < 1 || l_domainSize <= 0 || l_endTime <= 0 ) {
      std::cerr << "invalid configuration: check cell count, domain size, or end time" << std::endl;
      return EXIT_FAILURE;
    }
    
    // Calculate cell size using the domain argument
    l_dxy = l_domainSize / l_nx;
  }

  std::cout << "runtime configuration" << std::endl;
  std::cout << "  number of cells in x-direction: " << l_nx << std::endl;
  std::cout << "  number of cells in y-direction: " << l_ny << std::endl;
  std::cout << "  domain size:                    " << l_domainSize << std::endl;
  std::cout << "  cell size:                      " << l_dxy << std::endl;
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

  // construct setup
  tsunami_lab::setups::Setup *l_setup;
  l_setup = new tsunami_lab::setups::DamBreak1d( 5,
                                                  0,
                                                  2,
                                                  0,
                                                  50);

  /*l_setup = new tsunami_lab::setups::ShockShock1d( 5,
                                                  0,
                                                  50 );*/

  //  l_setup = new tsunami_lab::setups::Subcritical1d();
  //l_setup = new tsunami_lab::setups::Supercritical1d();

  // variable for maximum Froude number
  tsunami_lab::t_real maxF = -1.0;
  // variable for position x where maximum Froude number occurs
  tsunami_lab::t_real maxX = 0.0;

  //number of sample points for Froude number evaluation
  const int N = 1000; 
  // define interval for x values
  const tsunami_lab::t_real x0 = 0.0;
  const tsunami_lab::t_real x1 = 25.0;

  // loop over sample points
  for (int i = 0; i < N; i++) {
      // calculate x value for current sample point
      tsunami_lab::t_real x = x0 + i * (x1 - x0) / (N - 1);
      
      // evaluate Froude number at current x value

      tsunami_lab::t_real F = l_setup->getFroudeNumber(x, 0.0);

 
      
      // update maximum Froude number and corresponding x value if necessary
      if (F > maxF) {
          maxF = F;
          maxX = x;
      }
  }

std::cout << "\n=== Froude analysis at t = 0 ===" << std::endl;
std::cout << "Max Froude number: " << maxF << std::endl;
std::cout << "Position of max Froude number: x = " << maxX << std::endl;





  tsunami_lab::patches::WavePropagation *l_waveProp;
  l_waveProp = new tsunami_lab::patches::WavePropagation1d(l_nx , l_useFWaveSolver, tsunami_lab::patches::WavePropagation1d::BoundaryCondition::BoundaryRight);

  // maximum observed height in the setup
  tsunami_lab::t_real l_hMax = std::numeric_limits< tsunami_lab::t_real >::lowest();

  // set up solver
  for( tsunami_lab::t_idx l_cy = 0; l_cy < l_ny; l_cy++ ) {
    tsunami_lab::t_real l_y = l_cy * l_dxy; 

    for( tsunami_lab::t_idx l_cx = 0; l_cx < l_nx; l_cx++ ) {
      tsunami_lab::t_real l_x = l_cx * l_dxy; 

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
  tsunami_lab::t_real l_speedMax = std::sqrt( 9.81 * l_hMax );

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

      tsunami_lab::t_real F_left = 0.0; // Initialize Froude number for the left cell

      for (tsunami_lab::t_idx i = 0; i < l_nx; i++) {
   

        tsunami_lab::t_real x = i * l_dxy;
        tsunami_lab::t_real h = l_waveProp->getHeight()[i];
        tsunami_lab::t_real hu = l_waveProp->getMomentumX()[i];
        tsunami_lab::t_real u = hu / h;
        tsunami_lab::t_real F_right = std::abs(u) / std::sqrt(9.81 * h);

        if(i > 0) {
          if (F_left > 1.0 && F_right < 1.0) {
            std::cout << "Jump detected at x = " << x << " at time t = " << l_simTime << std::endl;
          }
        }

        F_left = F_right; // Update Froude number for the next iteration
      }

      tsunami_lab::io::Csv::write( l_dxy,
                                   l_nx,
                                   1,
                                   1,
                                   l_simTime,
                                   l_waveProp->getHeight(),
                                   l_waveProp->getBathymetry()+1,
                                   l_waveProp->getMomentumX(),
                                   nullptr,
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
