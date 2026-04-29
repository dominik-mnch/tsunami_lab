/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * One-dimensional supercritical problem.
 **/
#ifndef TSUNAMI_LAB_SETUPS_SUPERCRITICAL_1D_H
#define TSUNAMI_LAB_SETUPS_SUPERCRITICAL_1D_H

#include "../Setup.h"

namespace tsunami_lab {
  namespace setups {
    class Supercritical1d;
  }
}

/**
 * 1d supercritical setup.
 **/
class tsunami_lab::setups::Supercritical1d: public Setup {
  private:
    //! height on the left side 
    t_real m_heightLeft = 0;

    //! momentum on the left side
    t_real m_momentumLeft = 0;

    //! height on the right side
    t_real m_heightRight = 0;

    //! momentum on the right side
    t_real m_momentumRight = 0;

    //! location of the dam
    t_real m_locationDam = 0;

    //! bathymetry of the cells
    t_real * m_b = nullptr;

  public:
    /**
     * Constructor.
     **/
    Supercritical1d();
   

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
     * @param i_ix index of the cell
     * @return bathymetry at the given point.
     **/
    t_real getBathymetry( t_real i_ix,
                          t_real       ) const;

    /**
     * Gets the momentum in x-direction.
     *
     * @return momentum in x-direction.
     **/
    t_real getMomentumX( t_real i_x,
                         t_real ) const;

    /**
     * Gets the momentum in y-direction.
     *
     * @return momentum in y-direction.
     **/
    t_real getMomentumY( t_real,
                         t_real ) const;

};

#endif