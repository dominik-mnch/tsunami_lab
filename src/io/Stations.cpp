/**
 * @author Magalena Schwarzkopf
 * @author Dominik Münch
 * 
 * @section DESCRIPTION
 * Class for handling stations that write water height and momenta values at their locations to CSV files.
 **/

#include "Stations.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <filesystem>


tsunami_lab::io::Stations::Stations(t_real output_frequency) {
    
    std::ifstream file("./stations/Stations.csv");
    if (!file.is_open()) {
        std::cerr << "Could not open file\n";
        return;
    }

  // clean up output directory
  std::cout << "cleaning output directory" << std::endl;
  std::string l_outputDir = "./stations/output";
  if( std::filesystem::exists( l_outputDir ) ) {
    for( const auto& l_entry : std::filesystem::directory_iterator( l_outputDir ) ) {
      if( l_entry.is_regular_file() && l_entry.path().extension() == ".csv" ) {
        std::filesystem::remove( l_entry.path() );
      }
    }
  }

    std::string line; 
    std::istringstream iss(line);
    stations = std::vector<Station>();

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string xStr, yStr, name;

        std::getline(ss, name, ',');
        std::getline(ss, xStr, ',');
        std::getline(ss, yStr, ',');


        t_real x = 0, y = 0;

        try {
            x = std::stof(xStr);
            y = std::stof(yStr);
        } catch ( ... ) {
            std::cerr << "Error: Could not parse value '" << x << "' or '" << y << "' as a number" << std::endl;
            return;
        }

        stations.push_back(Station{
            name,
            x,
            y,
            std::ofstream("./stations/output/" + name + ".csv")
        });
    }
    file.close();
    this->output_frequency = output_frequency;
}

void writeCSVLine(tsunami_lab::t_real time,
                  tsunami_lab::t_real water_height,
                  tsunami_lab::t_real momentum_x,
                  tsunami_lab::t_real momentum_y,
                  std::ostream & io_stream) {

    io_stream << time << ","
              << water_height << ","
              << momentum_x << ","
              << momentum_y << "\n";
}



void tsunami_lab::io::Stations::writeToCSV( t_real time,
                                            t_real l_dx,
                                            t_real l_dy,
                                            tsunami_lab::patches::WavePropagation* l_waveProp) {


    if(time - lastTime < output_frequency) {
        return;
    }

    for(const Station& s : stations) {
        std::string filename = s.name + ".csv";

        // check if file exists / is empty
        std::ifstream check("./stations/output/" + filename);
        bool is_empty = !check.good() ||
                        check.peek() == std::ifstream::traits_type::eof();
        check.close();

        // open file
        std::ofstream csv_file("./stations/output/" + filename, std::ios::app);

        // write header if needed
        if (is_empty) {
            csv_file << "time,water_height,momentum_x,momentum_y\n";
        }

        // calculate the index of the location
        t_idx i_x = s.x/l_dx;
        t_idx i_y = s.y/l_dy;
        t_idx finalIndex = i_x + i_y * l_waveProp->getStride();

        // write data
        writeCSVLine(time, l_waveProp->getHeight()[finalIndex], l_waveProp->getMomentumX()[finalIndex], (l_waveProp->getMomentumY() == nullptr) ? 0.0 : l_waveProp->getMomentumY()[finalIndex], csv_file);
        csv_file << std::flush;
    }

    lastTime = time;
}