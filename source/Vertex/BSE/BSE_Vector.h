//
// Created by stkubr on 01.05.15.
//

#ifndef DSEPP_BSE_VECTOR_H
#define DSEPP_BSE_VECTOR_H

#include "BSE_Matrix.h"

class C_BSE_Vector: public C_BSE_Matrix {

public:
    C_BSE_Vector(){
        SetNameID("BSE_Vector",1);
        num_amplitudes=8; // number of amplitudes for PseudoScalar
        Initialization();
        threadloc_Projectors.resize(omp_get_max_threads(), std::vector<t_cmplxDirac>(num_amplitudes));
        threadloc_WeightCoeff.resize(omp_get_max_threads(), t_cmplxArray1D(num_amplitudes));
    }

    void SetDiracStructures(t_cmplxVector _k, t_cmplxVector _P, std::vector<t_cmplxDirac> & DiracStructure){
        t_cmplxDirac Y_T=C_Kinematics_1loop::TransIn(Y,_P);
        t_cmplxVector _k_T = C_Kinematics_1loop::TransIn(_k,_P);
        t_cmplx _k2_T=_k_T*_k_T;
        t_cmplx _k2=_k*_k;
        t_cmplx _k_P=_k*_P;
        t_cmplx _P2=_P*_P;

		(DiracStructure)[0]=Y_T;
		(DiracStructure)[1]=(_k_T*(_k_T*Z) - 1.0/3.0*Y_T*(_k2_T))/_k2;
		(DiracStructure)[2]=_k_T*(_P*Z)*_k_P/(_k2*_P2);
		(DiracStructure)[3]=-1.0*(Y_T*( (_P*Z)*(_k*Z) - (_k*Z)*(_P*Z) ) + 2.0*_k_T*(_P*Z) )/_k2/2.0;
		(DiracStructure)[4]=ii*_k_T*I/_k2;
		(DiracStructure)[5]=ii*(Y_T*(_k_T*Z) - (_k_T*Z)*Y_T)*_k_P/_k2;
		(DiracStructure)[6]=ii*(Y_T*(_P*Z) - (_P*Z)*Y_T)*(1.0 - _k_P*_k_P/_k2/_P2) - ii*2.0*_k_T*(_k_T*Z)*(_P*Z)/_k2;
		(DiracStructure)[7]=ii*_k_T*(_k_T*Z)*(_P*Z)/_k2;
    }

    t_cmplxMatrix DisentangleAmps(t_cmplxMatrix * pre_result){
        t_cmplxMatrix result(num_amplitudes,1);
        for (int i = 0; i < num_amplitudes; i++) result(i,0)=(*pre_result)(i,0)/threadloc_WeightCoeff[omp_get_thread_num()][i];
        return result;
    }
};


#endif //DSEPP_BSE_VECTOR_H
