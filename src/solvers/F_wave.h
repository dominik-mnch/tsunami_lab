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

    public:
        /**
         * Computes the net-updates.
         *
         * @param i_hL height of the left side.
         * @param i_hR height of the right side.
         * @param i_huL momentum of the left side.
         * @param i_huR momentum of the right side.
         * @param i_bL bathymetry of the left side.
         * @param i_bR bathymetry of the right side.
         * @param o_netUpdateL will be set to the net-updates for the left side; 0: height, 1: momentum.
         * @param o_netUpdateR will be set to the net-updates for the right side; 0: height, 1: momentum.
         **/
        static void netUpdates( t_real i_hL,
                                t_real i_hR,
                                t_real i_huL,
                                t_real i_huR,
                                t_real i_bL,
                                t_real i_bR,
                                t_real o_netUpdateL[2],
                                t_real o_netUpdateR[2] );
};

#endif