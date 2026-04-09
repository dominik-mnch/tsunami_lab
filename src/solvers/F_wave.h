/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 * 
 * @section DESCRIPTION
 * F-Wave solver for the one-dimensional shallow water equations.
 */

 #ifndef TSUNAMI_LAB_SOLVERS_F_WAVE
 #define TSUNAMI_LAB_SOLVERS_F_WAVE

 #include "../constants.h"

namespace tsunami_lab {
    namespace solvers {
        class F_wave;
    }
}

class tsunami_lab::solvers::F_wave {

    private:
    //! square root of gravity
    static t_real constexpr m_gSqrt = 3.131557121;

    /**
     * Computes the wave speeds.
     *
     * @param i_hL height of the left side.
     * @param i_hR height of the right side.
     * @param i_uL particle velocity of the leftside.
     * @param i_uR particles velocity of the right side.
     * @param o_waveSpeedL will be set to the speed of the wave propagating to the left.
     * @param o_waveSpeedR will be set to the speed of the wave propagating to the right.
     **/
    static void eigenvalues( t_real   i_hL,
                            t_real   i_hR,
                            t_real   i_uL,
                            t_real   i_uR,
                            t_real & o_waveSpeedL,
                            t_real & o_waveSpeedR );
};


#endif