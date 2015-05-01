//
// Created by stkubr on 22.04.15.
//

#ifndef DSEPP_BSE_PION_H
#define DSEPP_BSE_PION_H

#include "BSE_Matrix.h"

class C_BSE_PseudoScalar: public C_BSE_Matrix {

public:
    C_BSE_PseudoScalar(){
        SetNameID("BSE_PseudoScalar",1);
        num_amplitudes=4; // number of amplitudes for PseudoScalar
        Initialization();
        threadloc_Projectors.resize(omp_get_max_threads(), std::vector<t_cmplxDirac>(num_amplitudes));
        threadloc_WeightCoeff.resize(omp_get_max_threads(), t_cmplxArray1D(num_amplitudes));
    }

    void SetDiracStructures(t_cmplxVector _k, t_cmplxVector _P, std::vector<t_cmplxDirac> & DiracStructure){
        (DiracStructure)[0]=ii*Y5;
        (DiracStructure)[1]=Y5*(_P*Z);
        (DiracStructure)[2]=Y5*(_k*Z)*(_k*_P);
        (DiracStructure)[3]=Y5*((_k*SIG)*_P);
    }

    t_cmplxMatrix DisentangleAmps(t_cmplxMatrix * pre_result){
        t_cmplx p2 = threadloc_Momenta[omp_get_thread_num()].p2;
        t_cmplx N2_Factor = threadloc_Momenta[omp_get_thread_num()].N2_Factor;
        t_cmplx p_P = threadloc_Momenta[omp_get_thread_num()].p_P;
        t_cmplx P2 = threadloc_Momenta[omp_get_thread_num()].P2;
        t_cmplx weight = threadloc_WeightCoeff[omp_get_thread_num()][0];
		t_cmplxMatrix result(num_amplitudes,1);
		result(0,0)=1.0/weight*(*pre_result)(0,0);
		result(1,0)=1.0/weight*(p2/N2_Factor*(*pre_result)(1,0) - 1.0/N2_Factor*(*pre_result)(2,0));
		result(2,0)=1.0/weight*(-1.0/N2_Factor*(*pre_result)(1,0) + P2/p_P/p_P/N2_Factor*(*pre_result)(2,0));
		result(3,0)=1.0/weight*(-1.0/N2_Factor*(*pre_result)(3,0));
		return result;
    }
};


#endif //DSEPP_BSE_PION_H