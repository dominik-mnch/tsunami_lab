/**
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Two-dimensional Tsunami Event setup with netCDF input.
 **/
#ifndef TSUNAMI_LAB_SETUPS_TSUNAMI_EVENT_2D_H
#define TSUNAMI_LAB_SETUPS_TSUNAMI_EVENT_2D_H

#include "../Setup.h"
#include "../../io/NetCdf.h"

namespace tsunami_lab {
  namespace setups {
    class TsunamiEvent2d;
  }
}

/**
 * 2D Tsunami Event setup with netCDF support.
 **/
class tsunami_lab::setups::TsunamiEvent2d: public Setup {

  private:
    //! Grid data loaded from the bathymetry and displacement netCDF files.
    tsunami_lab::io::NetCdf::Data m_data;
    //! Epsilon for shallow water threshold.
    static constexpr t_real m_delta = 20.0;

  public:

    /**
     * Constructor for TsunamiEvent2d setup.
     *
     * @param i_bathymetryFile path to the bathymetry netCDF file.
     * @param i_displacementFile path to the displacement netCDF file.
     **/
    TsunamiEvent2d( std::string const & i_bathymetryFile,
                    std::string const & i_displacementFile );

    /**
     * Gets the water height at a given point.
     *
     * @param i_x x-coordinate.
     * @param i_y y-coordinate.
     * @return water height at the given point.
     **/
    t_real getHeight( t_real i_x,
                      t_real i_y ) const;

    /**
     * Gets the bathymetry at a given point.
     *
     * @param i_x x-coordinate.
     * @param i_y y-coordinate.
     * @return bathymetry at the given point.
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
