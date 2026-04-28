/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * One-dimensional Rare Rare problem.
 **/
#ifndef TSUNAMI_LAB_SETUPS_RARE_RARE_1D_H
#define TSUNAMI_LAB_SETUPS_RARE_RARE_1D_H

#include "../Setup.h"

namespace tsunami_lab {
  namespace setups {
    class RareRare1d;
  }
}

/**
 * 1d Rare Rare setup.
 **/
class tsunami_lab::setups::RareRare1d: public Setup {
  private:
    //! inital height of the water streams
    t_real m_height = 0;

    //! momentum of the water streams.
    t_real m_momentum = 0;

    //! location of the discontinuity
    t_real m_locationDiscontinuity = 0;

    //! true if a parabolic bathymetry should be used
    bool m_useBathymetryParabola = false;

    //! constant term of the parabola
    t_real m_bathymetryOffset = 0;

    //! x-position of the parabola's vertex
    t_real m_bathymetryCenter = 0;

    //! quadratic coefficient of the parabola
    t_real m_bathymetryScale = 0;

  public:
    /**
     * Constructor.
     *
     * @param i_height initial water height of the water streams.
     * @param i_momentum momentum of the water streams.
     * @param i_locationDiscontinuity location (x-coordinate) of the discontinuity.
     * @param i_useBathymetryParabola enable parabolic bathymetry.
     * @param i_bathymetryOffset constant term of the bathymetry parabola.
     * @param i_bathymetryCenter x-position of the parabola's vertex.
     * @param i_bathymetryScale quadratic coefficient of the bathymetry parabola.
     **/
    RareRare1d( t_real i_height,
                t_real i_momentum,
                t_real i_locationDiscontinuity,
                bool i_useBathymetryParabola = false,
                t_real i_bathymetryOffset = 0,
                t_real i_bathymetryCenter = 0,
          t_real i_bathymetryScale = 0 );

    /**
     * Gets the water height at a given point.
     *
     * @param i_x x-coordinate of the queried point.
     * @return height at the given point.
     **/
    t_real getHeight( t_real i_x,
                      t_real      ) const;

    /**
     * Gets the bathymetry at a given point.
     *
     * @return bathymetry at the given point.
     **/
    t_real getBathymetry( t_real,
                          t_real   ) const;
                          
    /**
     * Gets the momentum in x-direction.
     *
     * @param i_x x-coordinate of the queried point.
     * @return momentum in x-direction.
     **/
    t_real getMomentumX( t_real i_x,
                         t_real ) const;

    /**
     * Gets the momentum in y-direction.
     *
     * @param i_x x-coordinate of the queried point.
     * @return momentum in y-direction.
     **/
    t_real getMomentumY( t_real i_x,
                         t_real ) const;

};

#endif