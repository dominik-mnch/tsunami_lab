/**
 * @author Magadalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * IO-routines for writing a snapshot as Comma Separated Values (CSV).
 **/
#ifndef TSUNAMI_LAB_IO_NETCDF
#define TSUNAMI_LAB_IO_NETCDF

#include "../constants.h"
#include "../patches/WavePropagation.h"
#include <cstring>
#include <string>
#include <netcdf.h>

namespace tsunami_lab {
  namespace io {
    class NetCdf;
  }
}

class tsunami_lab::io::NetCdf {
  private:
    //! Cell width in x-direction.
    t_real m_dx;
    //! Cell width in y-direction.
    t_real m_dy;
    //! Origin of the x-coordinate system.
    t_real m_originX;
    //! Origin of the y-coordinate system.
    t_real m_originY;
    //! Number of interior cells in x-direction.
    t_idx m_nx;
    //! Number of interior cells in y-direction.
    t_idx m_ny;
    //! Row stride of the field arrays.
    t_idx m_stride;
    //! Stored simulation time value.
    t_real m_time;
    //! Pointer to water height values.
    t_real const * m_h;
    //! Pointer to bathymetry values.
    t_real const * m_b;
    //! Pointer to x-momentum values.
    t_real const * m_hu;
    //! Pointer to y-momentum values.
    t_real const * m_hv;
    //! Current index in the unlimited time dimension.
    t_idx m_timeStep;
    //! netCDF file handle identifier.
    t_idx m_ncId;
    //! Variable identifier of the time variable.
    t_idx m_varTimeId;
    //! Variable identifier of the bathymetry variable.
    t_idx m_varBathyId;
    //! Variable identifier of the water height variable.
    t_idx m_varHeightId;
    //! Variable identifier of the x-momentum variable.
    t_idx m_varMomentumXId;
    //! Variable identifier of the y-momentum variable.
    t_idx m_varMomentumYId;

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
     * @param i_h pointer to water height of the cells
     * @param i_b pointer to bathymetry of the cells
     * @param i_hu pointer to momentum in x-direction of the cells
     * @param i_hv pointer to momentum in y-direction of the cells
    * @param i_filePath output path of the netCDF file.
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
            t_real       const * i_hv,
            std::string const & i_filePath = "solutions/solution.nc" );

        /**
         * Closes the netCDF file.
         **/
        ~NetCdf();

    /** 
     * Writes the current time step to the netCDF file.
     * 
     * @param simTime simulation time of the snapshot to write.
     **/
    void writeTimeStep(t_real simTime);
};



#endif