//
// Created by stkubr on 22.04.15.
//

#ifndef DSEPP_BSE_PION_H
#define DSEPP_BSE_PION_H


#include "BSE_TwoBody.h"

class C_BSE_Pion: public C_BSE_TwoBody {

public:
    C_BSE_Pion(){
        //SetNameID("BSE_PseudoScalar",1);
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

    virtual C_BSE_Pion * MakeCopy()
    {
        return new C_BSE_Pion(*this);
    }

    t_cmplxMatrix DisentangleAmps(t_cmplxMatrix * pre_result){
        t_cmplx p2 = threadloc_Momenta[omp_get_thread_num()].p2;
        t_cmplx N2_Factor = threadloc_Momenta[omp_get_thread_num()].N2_Factor;
        t_cmplx p_P = threadloc_Momenta[omp_get_thread_num()].p_P;
        t_cmplx P2 = threadloc_Momenta[omp_get_thread_num()].P2;
        t_cmplx weight = threadloc_WeightCoeff[omp_get_thread_num()][0];
		t_cmplxMatrix result(4,1);
		result(0,0)=1.0/weight*(*pre_result)(0,0);
		result(1,0)=1.0/weight*(p2/N2_Factor*(*pre_result)(1,0) - 1.0/N2_Factor*(*pre_result)(2,0));
		result(2,0)=1.0/weight*(-1.0/N2_Factor*(*pre_result)(1,0) + P2/p_P/p_P/N2_Factor*(*pre_result)(2,0));
		result(3,0)=1.0/weight*(-1.0/N2_Factor*(*pre_result)(3,0));
		return result;
/*		#endif

#if BASIS_TYPE == 2

		t_cmplxMatrix result(num_amplitudes,1);
		for (int i = 0; i < num_amplitudes; i++) result(i,0)=(*pre_result)(i,0)/WeightCoeff[i];
		return result;

		#endif*/
    }

    void setPropagators(t_cmplxVector & K_plus, t_cmplxVector  & K_minus, t_cmplxDirac & S_p, t_cmplxDirac & S_m)
    {
        t_cmplxArray1D quark_temp_sigma;
        quark_temp_sigma= Parton_P->PropagatorAtPoint(K_plus * K_plus);
        S_p=(-1.0*ii*(K_plus*Z)*quark_temp_sigma[3] + I*quark_temp_sigma[4]);

        quark_temp_sigma= Parton_M->PropagatorAtPoint(K_minus * K_minus);
        S_m=(-1.0*ii*(K_minus*Z)*(quark_temp_sigma[3]) + I*(quark_temp_sigma[4]));
    }
};


#endif //DSEPP_BSE_PION_H
