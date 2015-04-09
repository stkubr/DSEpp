/*
 * RLandPseudoScalar.hpp
 *
 *  Created on: Jan 28, 2014
 *      Author: stkubr
 */

#ifndef RLANDPSEUDOSCALAR_HPP_
#define RLANDPSEUDOSCALAR_HPP_


#include "TwoQuarkKernel.hpp"

class C_Kernel_RL_PS: public C_TwoQuarkKernel{
	friend class C_Kernel_Factory;
protected:
	t_cmplxTensor (C_Kernel_RL_PS::*SetPSMatrix) (t_cmplx&,t_cmplx&);
    t_cmplx Z2;
    t_cmplx PseudoMesonMass;
	double PionDecayConst;
	PS_type_ID Exchange_type_ID;

	C_Kernel_RL_PS(){
		SetNameID("Kernel RainbowLadder + PseudoScalar exchange",1);
		Memory->VertexDressings.resize(2, t_cmplxArray3D(1));
        Z2=1.0;
        PseudoMesonMass = 0.138;
		PionDecayConst = 0.093;
	}

    void setMesonExchangeMass(t_cmplx _M);

    void setMediators(t_cmplxVector& k, t_cmplxVector& p, t_cmplxVector& P, std::vector<t_cmplxTensor>& Mediators);

	t_cmplx ElementKmatrix(int t, int s, int r, int u,
							  std::vector<t_cmplxTensor>& Mediators);

    t_cmplx VertexDressingAt(int kernel_type, int num_P, int num_amp, t_cmplx coordin);

    t_cmplx SetInterpolation(t_cmplx vertex_momenta, t_cmplx prop_momenta);

	t_cmplxTensor SetPSMatrix_etta_quark(t_cmplx& vertex_momenta, t_cmplx& prop_momenta);

	t_cmplxTensor SetPSMatrix_etta_BSE(t_cmplx& vertex_momenta, t_cmplx& prop_momenta);

	t_cmplxTensor SetPSMatrix_pion_quark(t_cmplx& vertex_momenta, t_cmplx& prop_momenta);

	t_cmplxTensor SetPSMatrix_pion_BSE(t_cmplx& vertex_momenta, t_cmplx& prop_momenta);

public:

	void info() { std::cout << "Kernel RainbowLadder + PseudoScalar exchange" << std::endl; }

	void setExchangeID(PS_type_ID exchange_id){
		Exchange_type_ID=exchange_id;
	}

	void setConvolutionType(int type);
};

#endif /* RLANDPSEUDOSCALAR_HPP_ */
