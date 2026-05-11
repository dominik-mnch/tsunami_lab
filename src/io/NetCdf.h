/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * IO-routines for writing a snapshot as Comma Separated Values (CSV).
 **/
#ifndef TSUNAMI_LAB_IO_NETCDF
#define TSUNAMI_LAB_IO_NETCDF

#include "../constants.h"
#include "../patches/WavePropagation.h"
#include <cstring>
#include <netcdf>

namespace tsunami_lab {
  namespace io {
    class NetCdf;
  }
}

class tsunami_lab::io::NetCdf {
  private:
    t_real               m_dx;
    t_real               m_dy;
    t_real               m_originX;
    t_real               m_originY;
    t_idx                m_nx;
    t_idx                m_ny;
    t_idx                m_stride;
    t_real               m_time;
    t_real       const * m_h;
    t_real       const * m_b;
    t_real       const * m_hu;
    t_real       const * m_hv;
    netCDF::NcFile m_file;

  public:

    /**
     * Initializes the netCDF writer, opens the file and writes the header.
     *
     * @param i_dx cell width in x-direction.
     * @param i_dy cell width in y-direction.
     * @param i_originX origin of the x-coordinate system.
     * @param i_originY origin of the y-coordinate system.
     * @param i_nx number of cells in x-direction.
     * @param i_ny number of cells in y-direction.
     * @param i_stride stride of the data arrays in y-direction (x is assumed to be stride-1).
     * @param i_time simulation time of the snapshot.
     * @param i_h water height of the cells; optional: use nullptr if not required.
     * @param i_b bathymetry of the cells; optional: use nullptr if not required.
     * @param i_hu momentum in x-direction of the cells; optional: use nullptr if not required.
     * @param i_hv momentum in y-direction of the cells; optional: use nullptr if not required.
     **/

    NetCdf( t_real               i_dx,
            t_real               i_dy,
            t_real               i_originX,
            t_real               i_originY,
            t_idx                i_nx,
            t_idx                i_ny,
            t_idx                i_stride,
            t_real               i_time,
            t_real       const * i_h,
            t_real       const * i_b,
            t_real       const * i_hu,
            t_real       const * i_hv);

    /** 
     * Writes the current time step to the netCDF file.
     **/
    void writeTimeStep(t_real simTime);
};



#endif