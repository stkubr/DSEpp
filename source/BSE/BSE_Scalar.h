//
// Created by stkubr on 01.05.15.
//

#ifndef DSEPP_BSE_SCALAR_H
#define DSEPP_BSE_SCALAR_H

#include "BSE_Matrix.h"

namespace BSE {

    class C_BSE_Scalar : public C_BSE_Matrix {

    public:
        C_BSE_Scalar() {
            SetNameID("BSE_Scalar", 1);
            num_amplitudes = 4; // number of amplitudes for Scalar
            Initialization();
            threadloc_Projectors.resize(omp_get_max_threads(), std::vector<t_cmplxDirac>(num_amplitudes));
            threadloc_WeightCoeff.resize(omp_get_max_threads(), t_cmplxArray1D(num_amplitudes));
        }

    protected:
        void setDiracStructures(t_cmplxVector _k, t_cmplxVector _P, std::vector<t_cmplxDirac> &DiracStructure);

        t_cmplxMatrix disentangleAmps(t_cmplxMatrix *pre_result);
    };

}
#endif //DSEPP_BSE_SCALAR_H
