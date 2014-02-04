/*
 * RainbowLadderKernel.cpp
 *
 *  Created on: Jan 28, 2014
 *      Author: stkubr
 */

#include "RainbowLadderKernel.hpp"

C_Kernel_RL::C_Kernel_RL(){
	SetNameID("Kernel RainbowLadder",1);
	Memory->VertexDressings.resize(2, t_cmplxArray3D(1));
}

void C_Kernel_RL::SetMediators(t_cmplxVector& k, t_cmplxVector& p, t_cmplxVector& P, std::vector<t_cmplxTensor>& Mediators){
	t_cmplx k2_product;
	t_cmplxTensor Gluon_Matrix(2);
	k2_product=(k*k);
	t_cmplx Gluon_factor=Z2*Z2*4.0/3.0*Gluon->GetGluonAt(k2_product);
	Gluon_Matrix=Gluon_factor*(g-((k)%(k))/(k2_product));

	(Mediators).resize(1);
	(Mediators)[0]=Gluon_Matrix;
}

t_cmplx C_Kernel_RL::getElementKmatrix(int t, int s, int r, int u, std::vector<t_cmplxTensor>& Mediators){
	t_cmplx Gluon_contrib=0.0;
	for (int mu = 0; mu < 4; mu++){
		for (int nu = 0; nu < 4; nu++){
			Gluon_contrib+=(Y.Element(t,s).Element(mu) * Mediators[0].Element(mu,nu) * Y.Element(r,u).Element(nu) );
		}
	}
	return Gluon_contrib;
}
