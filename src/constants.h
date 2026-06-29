/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Constants / typedefs used throughout the code.
 **/
#ifndef TSUNAMI_LAB_CONSTANTS_H
#define TSUNAMI_LAB_CONSTANTS_H

#include <cstddef>

namespace tsunami_lab {
  //! integral type for cell-ids, pointer arithmetic, etc.
  typedef std::size_t t_idx;

  //! floating point type
  typedef float t_real;

  //! Wet/dry threshold: cells with h below this are considered dry and clamped to 0
  //! This prevents NaN generation from division by zero in velocity calculations (u=hu/h)
  static constexpr t_real WET_DRY_THRESHOLD = 1.0e-3f;
}

#endif