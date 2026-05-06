/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * IO-class that summarizes a collection of stations to write water height and momenta values to a designated CSV file
 * at specified time intervals. The station locations and names are read from a CSV file "Stations.csv" (runtime configuration).
 */

#ifndef TSUNAMI_LAB_IO_STATIONS_H
#define TSUNAMI_LAB_IO_STATIONS_H

#include "../constants.h"
#include "../patches/WavePropagation.h"
#include <string>
#include <fstream>
#include <vector>

namespace tsunami_lab {
  namespace io {
    class Stations;
  }
}

/**
 * @brief Represents an output station for simulation data.
 *
 * @var name Name of the station.
 *
 * @var x X-coordinate in simulation space.
 *
 * @var y Y-coordinate in simulation space.
 *
 * @var output File stream used for writing station measurements.
 */
typedef struct {
  std::string name;
  tsunami_lab::t_real x;
  tsunami_lab::t_real y;
  std::ofstream output;
} Station;


class tsunami_lab::io::Stations {
  private:

    //! list of stations that write to CSV files
    std::vector<Station> stations;

    //! frequency at which the stations write to their CSV files
    t_real output_frequency;

    //! last time at which the stations wrote to their CSV files
    t_real lastTime = 0.0;

  public:

    /**
     * Constructor.
     * 
     * @param output_frequency frequency at which the stations write to their CSV files.
     */
    Stations(t_real output_frequency);

    /**
     * Writes the current water height and momenta values at the stations' locations to the stations' CSV files
     * if the time since the last write exceeds the output frequency.
     * 
     * @param time current simulation time.
     * @param dx cell width in x-direction.
     * @param dy cell width in y-direction.
     * @param waveProp pointer to the wave propagation solver to get the current water height and momenta values at the stations' locations.
     */
    void writeToCSV(t_real time,
                    t_real dx,
                    t_real dy,
                    tsunami_lab::patches::WavePropagation *waveProp);
};
#endif