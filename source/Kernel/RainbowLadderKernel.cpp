/*
 * RainbowLadderKernel.cpp
 *
 *  Created on: Jan 28, 2014
 *      Author: stkubr
 */

#include "RainbowLadderKernel.hpp"

using namespace Kernels;

// Sets gluon propagating inside of Kernel
void C_Kernel_RL::setMediators(t_cmplxVector& k, t_cmplxVector& p, t_cmplxVector& P, std::vector<t_cmplxTensor>& Mediators){
	t_cmplx k2_product;
	t_cmplxTensor Gluon_Matrix(2);
    t_cmplx Z2 = Propagators[1]->Z2Factor();
	k2_product=(k*k);
	t_cmplx Gluon_factor=Z2*Z2*4.0/3.0* Propagators[0]->PropagatorAtPoint(k2_product)[0];
	Gluon_Matrix=Gluon_factor*(g-((k)%(k))/(k2_product));

	(Mediators).resize(1);
	(Mediators)[0]=Gluon_Matrix;
}

// Set a (t,s,r,u) element of the K_matrix for a given class
t_cmplx C_Kernel_RL::ElementKmatrix(int t, int s, int r, int u, std::vector<t_cmplxTensor>& Mediators){
	t_cmplx Gluon_contrib=0.0;
	for (int mu = 0; mu < 4; mu++){
		for (int nu = 0; nu < 4; nu++){
			Gluon_contrib+=(Y.Element(t,s).Element(mu) * Mediators[0].Element(mu,nu) * Y.Element(r,u).Element(nu) );
		}
	}
	return Gluon_contrib;
}
