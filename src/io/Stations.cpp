/**
 * @author Magalena Schwarzkopf
 * @author Dominik Münch
 */

#include "Stations.h"
#include <fstream>
#include <iostream>



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




tsunami_lab::io::Stations::Stations(t_real i_x, t_real i_y, std::string i_name)
    : m_x(i_x), m_y(i_y), m_name(i_name) {
}



void tsunami_lab::io::Stations::writeToCSV(t_real time,
                          t_real water_height,
                          t_real momentum_x,
                          t_real momentum_y,
                          std::ostream & /*io_stream*/) {

    std::string filename = m_name + ".csv";

    // check if file exists / is empty
    std::ifstream check(filename);
    bool is_empty = !check.good() ||
                    check.peek() == std::ifstream::traits_type::eof();
    check.close();

    // open file 
    std::ofstream csv_file(filename, std::ios::app);

    // write header if needed
    if (is_empty) {
        csv_file << "time,water_height,momentum_x,momentum_y\n";
    }

    // write data
    writeCSVLine(time, water_height, momentum_x, momentum_y, csv_file);

    csv_file << std::flush;
}




