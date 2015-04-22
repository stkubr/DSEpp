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
        Amplitudes.resize(4);
        Projectors.resize(4);
        WaveFunctions.resize(4);
        WeightCoeff.resize(4);
        //SaveBSEPath="Data_files/SaveBSE_PseudoScalar.dat";
    }

    void SetDiracStructures(t_cmplxVector _k, t_cmplxVector _P, std::vector<t_cmplxDirac> * DiracStructure){
        //#if BASIS_TYPE == 1
        (*DiracStructure)[0]=ii*Y5;
        (*DiracStructure)[1]=Y5*(_P*Z);
        (*DiracStructure)[2]=Y5*(_k*Z)*(_k*_P);
        (*DiracStructure)[3]=Y5*((_k*SIG)*_P);

        //#endif

        /*#if BASIS_TYPE == 2

        Y_T=Momenta.TransIn(Y,_P);
        t_cmplxVector _k_T;
        _k_T=Momenta.TransIn(_k,_P);
        t_cmplx _k2_T,_k2, _k_P,_P2;
        _k2_T=_k_T*_k_T;
        _k2=_k*_k;
        _k_P=_k*_P;
        _P2=_P*_P;

        (*DiracStructure)[0]=Y5;
        (*DiracStructure)[1]=Y5*(_P*Z);
        (*DiracStructure)[2]=Y5*(_k_T*Z);
        (*DiracStructure)[3]=ii/2.0*Y5*((_k_T*Z)*(_P*Z) - (_P*Z)*(_k_T*Z));

        #endif*/
    }

    virtual C_BSE_Pion * MakeCopy()
    {
        return new C_BSE_Pion(*this);
    }

    t_cmplxMatrix DisentangleAmps(t_cmplxMatrix * pre_result){
//#if BASIS_TYPE == 1
		t_cmplxMatrix result(4,1);
		result(0,0)=1.0/WeightCoeff[0]*(*pre_result)(0,0);
		result(1,0)=1.0/WeightCoeff[0]*(Momenta.p2/Momenta.N2_Factor*(*pre_result)(1,0) - 1.0/Momenta.N2_Factor*(*pre_result)(2,0));
		result(2,0)=1.0/WeightCoeff[0]*(-1.0/Momenta.N2_Factor*(*pre_result)(1,0) + Momenta.P2/Momenta.p_P/Momenta.p_P/Momenta.N2_Factor*(*pre_result)(2,0));
		result(3,0)=1.0/WeightCoeff[0]*(-1.0/Momenta.N2_Factor*(*pre_result)(3,0));
		return result;
/*		#endif

#if BASIS_TYPE == 2

		t_cmplxMatrix result(num_amplitudes,1);
		for (int i = 0; i < num_amplitudes; i++) result(i,0)=(*pre_result)(i,0)/WeightCoeff[i];
		return result;

		#endif*/
    }

    void setPropagators(t_cmplxVector *K_plus, t_cmplxVector *K_minus)
    {
        t_cmplxArray1D quark_temp_sigma;
        quark_temp_sigma= Parton_P->PropagatorAtPoint((*K_plus) * (*K_plus));
        S_p=(-1.0*ii*((*K_plus)*Z)*quark_temp_sigma[3] + I*quark_temp_sigma[4]);

        quark_temp_sigma= Parton_M->PropagatorAtPoint((*K_minus) * (*K_minus));
        S_m=(-1.0*ii*((*K_minus)*Z)*(quark_temp_sigma[3]) + I*(quark_temp_sigma[4]));
    }
};


#endif //DSEPP_BSE_PION_H
