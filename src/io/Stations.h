#ifndef TSUNAMI_LAB_IO_STATIONS_H
#define TSUNAMI_LAB_IO_STATIONS_H

#include "../constants.h"
#include <string>
#include <ostream>


namespace tsunami_lab {
  namespace io {
    class Stations;
  }
}

class tsunami_lab::io::Stations {
    public:
        Stations(t_real i_x, t_real i_y, std::string i_name);

        void writeToCSV(t_real time,
                        t_real water_height,
                        t_real momentum_x,
                        t_real momentum_y,
                        std::ostream & io_stream);

    private:
        t_real m_x;
        t_real m_y;
        std::string m_name;


};



#endif