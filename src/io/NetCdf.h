/**
 * @author Magadalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * IO-routines for reading and writing netCDF files.
 **/
#ifndef TSUNAMI_LAB_IO_NETCDF
#define TSUNAMI_LAB_IO_NETCDF

#include "../constants.h"
#include "../patches/WavePropagation.h"
#include <cstring>
#include <string>
#include <vector>
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
    //! Number of cells that get averaged to reduce output (k x k block gets averaged to 1 cell in output).
    t_idx m_k;
    //! Row stride of the field arrays.
    t_idx m_stride;
    //! Stored simulation time value.
    t_real m_time;
    //! Stored simulation end time.
    t_real m_endTime;
    //! Solver mode written to checkpoint (1 for F-Wave, 0 for Roe).
    int m_solverMode;
    //! Propagation mode written to checkpoint.
    std::string m_propagation;
    //! Setup name written to checkpoint.
    std::string m_setup;
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
    //! netCDF file handle for the checkpoint file (c_invalidId if no checkpoint).
    t_idx m_checkpointNcId;
    //! Simulation time of the last written time step (used by overwriteCheckpointEndTime).
    t_real m_lastSimTime;


    /**
     * Helper function that writes a 2d field (height, momentum, bathymetry) to the netCDF file.
     * Averages k x k blocks to reduce output resolution and writes the data of the current time step.
     * 
     * @param i_data pointer to the field data to write (row-major, ny * stride).
     * @param i_varId variable identifier of the netCDF variable to write to.
     */
    void helperWritingData( t_real const * i_data, t_idx i_varId ) const;
public:

    /**
     * Holds grid data loaded from netCDF files and provides nearest-neighbor lookups.
     **/
    struct Data {
      //! Grid coordinates in x-direction (bathymetry grid).
      std::vector<t_real> xCoords;
      //! Grid coordinates in y-direction (bathymetry grid).
      std::vector<t_real> yCoords;
      //! Bathymetry values on the grid (row-major, y * gridNx + x).
      std::vector<t_real> bathymetryData;
      //! Displacement values on the displacement grid (row-major, y * dispNx + x).
      std::vector<t_real> displacementData;
      //! x-coordinates of the displacement grid (may differ from bathymetry).
      std::vector<t_real> dispXCoords;
      //! y-coordinates of the displacement grid (may differ from bathymetry).
      std::vector<t_real> dispYCoords;
      //! Number of grid points in x-direction (bathymetry grid).
      t_idx gridNx = 0;
      //! Number of grid points in y-direction (bathymetry grid).
      t_idx gridNy = 0;
      //! Number of grid points in x-direction (displacement grid).
      t_idx dispNx = 0;
      //! Number of grid points in y-direction (displacement grid).
      t_idx dispNy = 0;

      t_real getBathymetry( t_real i_x, t_real i_y ) const;
      t_real getDisplacement( t_real i_x, t_real i_y ) const;
    };

    /**
     * Reads bathymetry (and optionally displacement) from netCDF files.
     *
     * @param i_bathymetryFile path to the netCDF file containing bathymetry data.
     * @param i_displacementFile path to the netCDF file containing displacement data (optional).
     * @return Data struct with all loaded grid data and lookup helpers.
     **/
    static Data read( std::string const & i_bathymetryFile,
                      std::string const & i_displacementFile = "" );

    /**
     * Holds all data loaded from a checkpoint + solution file pair.
     **/
    struct CheckpointData {
      //! Number of cells in x-direction.
      t_idx nx = 0;
      //! Number of cells in y-direction.
      t_idx ny = 0;
      //! Cell width in x-direction.
      t_real dx = 0;
      //! Cell width in y-direction.
      t_real dy = 0;
      //! Origin of the x-coordinate system.
      t_real originX = 0;
      //! Origin of the y-coordinate system.
      t_real originY = 0;
      //! Simulation time of the checkpoint.
      t_real simTime = 0;
      //! Simulation end time.
      t_real endTime = 0;
      //! Number of cells that get averaged to reduce output.
      t_idx k = 1;
      //! Solver mode (1 for F-Wave, 0 for Roe).
      int solverMode = 1;
      //! Propagation mode.
      std::string propagation;
      //! Setup name.
      std::string setup;
      //! Water height at the last valid checkpoint step (row-major, ny * nx).
      std::vector<t_real> h;
      //! x-momentum at the last valid checkpoint step (row-major, ny * nx).
      std::vector<t_real> hu;
      //! y-momentum at the last valid checkpoint step (row-major, ny * nx). Empty for 1d.
      std::vector<t_real> hv;
      //! Bathymetry (row-major, ny * nx).
      std::vector<t_real> b;
    };

    /**
     * Reads simulation parameters from a checkpoint file and the last valid
     * state (h, hu, hv, b) from the sibling solution.nc.
     *
     * @param i_checkpointFile path to the checkpoint netCDF file.
     * @return CheckpointData struct with all data needed to resume.
     **/
    static CheckpointData readCheckpoint( std::string const & i_checkpointFile );

    /**
     * Initializes the netCDF writer, opens the file and writes the header.
     *
     * @param i_dx cell width in x-direction.
     * @param i_dy cell width in y-direction.
     * @param i_originX origin of the x-coordinate system.
     * @param i_originY origin of the y-coordinate system.
     * @param i_nx number of cells in x-direction.
     * @param i_ny number of cells in y-direction.
     * @param i_k number of cells that get averaged to reduce output 
     * @param i_stride stride of the data arrays in y-direction (x is assumed to be stride-1).
     * @param i_time simulation time of the snapshot.
     * @param i_endTime simulation end time (for checkpoint)
    * @param i_solverMode solver mode (1 for F-Wave, 0 for Roe)
    * @param i_propagation propagation mode
    * @param i_setup setup name
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
            t_idx                i_k,
            t_idx                i_stride,
            t_real               i_time,
            t_real               i_endTime,
            int                  i_solverMode,
            std::string const &  i_propagation,
            std::string const &  i_setup,
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
    void writeTimeStep( t_real simTime );

    /**
     * Creates a checkpoint file and writes all simulation parameters as global attributes.
     * The sim_time attribute is initialised to 0 and updated after every write cycle.
     *
     * @param i_filePath output path of the checkpoint file.
     **/
    void defineCheckpoint( std::string const & i_filePath );

    /**
     * Overwrites the current simulation time of the netCDF checkpoint file.
     **/
    void overwriteCheckpointSimTime();
};



#endif