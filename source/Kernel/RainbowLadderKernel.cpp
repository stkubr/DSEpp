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

void C_Kernel_RL::SetKmatrix(t_cmplxMatrix2D& K_matrix, std::vector<t_cmplxTensor>& Mediators){
	for (int t = 0; t < 4; t++){
		for (int s = 0; s < 4; s++){
			for (int r = 0; r < 4; r++){
				for (int u = 0; u < 4; u++){
					t_cmplx Gluon_contrib=0.0;
					for (int i = 0; i < 4; i++){
						for (int j = 0; j < 4; j++){
							Gluon_contrib+=(Y.Element(t,s).Element(i) * Mediators[0].Element(i,j) * Y.Element(r,u).Element(j) );
						}
					}
					K_matrix(t,s)(r,u)=Gluon_contrib;
				}
			}
		}
	}
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


