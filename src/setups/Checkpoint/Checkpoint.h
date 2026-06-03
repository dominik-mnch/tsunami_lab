/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Checkpoint setup that loads a netCDF checkpoint file.
 **/
#ifndef TSUNAMI_LAB_SETUPS_CHECKPOINT_H
#define TSUNAMI_LAB_SETUPS_CHECKPOINT_H

#include "../Setup.h"
#include <string>
#include <vector>

namespace tsunami_lab {
  namespace setups {
    class Checkpoint;
  }
}

/**
 * Checkpoint setup that loads a netCDF checkpoint file.
 **/
class tsunami_lab::setups::Checkpoint: public Setup {

  private:
    //! Number of cells in x-direction.
    t_idx m_nx = 0;
    //! Number of cells in y-direction.
    t_idx m_ny = 0;
    //! Cell width in x-direction.
    t_real m_dx = 0;
    //! Cell width in y-direction.
    t_real m_dy = 0;
    //! Origin of the x-coordinate system.
    t_real m_originX = 0;
    //! Origin of the y-coordinate system.
    t_real m_originY = 0;
    //! Simulation end time read from checkpoint.
    t_real m_endTime = 0;
    //! Loaded water height values (row-major, ny * nx).
    std::vector<t_real> m_h;
    //! Loaded x-momentum values (row-major, ny * nx).
    std::vector<t_real> m_hu;
    //! Loaded y-momentum values (row-major, ny * nx).
    std::vector<t_real> m_hv;
    //! Loaded bathymetry values (row-major, ny * nx).
    std::vector<t_real> m_b;

    /**
     * Returns the flat array index for coordinates (i_x, i_y).
     **/
    t_idx getIndex( t_real i_x, t_real i_y ) const;

  public:

    /**
     * Constructor: reads checkpoint.nc for simulation parameters and
     * solution.nc (in the same directory) for the last valid state.
     *
     * @param i_checkpointFile path to the checkpoint netCDF file.
     */
    Checkpoint( std::string const & i_checkpointFile );

    /**
     * Gets the number of cells in x-direction. 
     */
    t_idx  getNx()      const { return m_nx; }

    /**
     * Gets the number of cells in y-direction.
     */
    t_idx  getNy()      const { return m_ny; }

    /**
     * Gets the cell width in x-direction.
     */
    t_real getDx()      const { return m_dx; }

    /**
     * Gets the cell width in y-direction.
     */
    t_real getDy()      const { return m_dy; }

    /**
     * Gets the origin of the x-coordinate system.
     */
    t_real getOriginX() const { return m_originX; }

    /**
     * Gets the origin of the y-coordinate system.
     */
    t_real getOriginY() const { return m_originY; }

    /**
     * Gets the simulation end time read from checkpoint.
     */
    t_real getEndTime() const { return m_endTime; }

    /**
     * Gets the water height at a given point.
     *
     * @param i_x x-coordinate.
     * @param i_y y-coordinate.
     * @return water height.
     **/
    t_real getHeight( t_real i_x,
                      t_real i_y ) const;

    /**
     * Gets the bathymetry at a given point.
     *
     * @param i_x x-coordinate.
     * @param i_y y-coordinate.
     * @return bathymetry.
     **/
    t_real getBathymetry( t_real i_x,
                          t_real i_y ) const;

    /**
     * Gets the momentum in x-direction at a given point.
     *
     * @param i_x x-coordinate.
     * @param i_y y-coordinate.
     * @return momentum in x-direction.
     **/
    t_real getMomentumX( t_real i_x,
                         t_real i_y ) const;

    /**
     * Gets the momentum in y-direction at a given point.
     *
     * @param i_x x-coordinate.
     * @param i_y y-coordinate.
     * @return momentum in y-direction.
     **/
    t_real getMomentumY( t_real i_x,
                         t_real i_y ) const;
};

#endif
