/*
 * RLandPseudoScalar.cpp
 *
 *  Created on: Jan 28, 2014
 *      Author: stkubr
 */

#include "RLandPseudoScalar.hpp"

    void C_Kernel_RL_PS::setMesonExchangeMass(t_cmplx _M){
        PseudoMesonMass=_M;
    }

	void C_Kernel_RL_PS::setMediators(t_cmplxVector& k, t_cmplxVector& p, t_cmplxVector& P, 
									  std::vector<t_cmplxTensor>& Mediators){
		t_cmplx k2_product;
		k2_product=(k*k);
		t_cmplx Gluon_factor=1.0;
		t_cmplx Z2 = Propagators[1]->DressingFactor();
		t_cmplxTensor Gluon_Matrix(2),PS_Matrix(0);
		Gluon_factor=Z2*Z2*4.0/3.0* Propagators[0]->PropagatorAtPoint(k2_product)[0];
		Gluon_Matrix=Gluon_factor*(g-((k)%(k))/(k2_product));
		t_cmplx PS_momenta=(P)*(P);

		PS_Matrix=(this->*SetPSMatrix)(PS_momenta,k2_product);

		(Mediators).resize(2);
		(Mediators)[0]=Gluon_Matrix;
		(Mediators)[1]=PS_Matrix;
	}

	t_cmplx C_Kernel_RL_PS::ElementKmatrix(int t, int s, int r, int u, std::vector<t_cmplxTensor>& Mediators){
		t_cmplx Gluon_contrib=0.0;
		t_cmplx Pion_contrib=0.0;
		for (int mu = 0; mu < 4; mu++){
			for (int nu = 0; nu < 4; nu++){
				Gluon_contrib+=(Y.Element(t,s).Element(mu) * Mediators[0].Element(mu,nu) * Y.Element(r,u).Element(nu) );
			}
		}
		Pion_contrib=(Y5.Element(t,s)*Y5.Element(r,u)) * Mediators[1];
		return Gluon_contrib + Pion_contrib;
	}

    t_cmplx C_Kernel_RL_PS::VertexDressingAt(int kernel_type, int num_P, int num_amp, t_cmplx coordin){
        t_cmplx result;
        t_cmplx F1,N,temp;
        t_cmplx z_i,dz_i;
        F1=t_cmplx(0.0,0.0);
        N=t_cmplx(0.0,0.0);
        for (int j=1;j<= Memory->VertexDressings[kernel_type][num_P][1].size();j++){
            z_i=Memory->VertexDressings[kernel_type][num_P][1][j-1];
            dz_i=Memory->VertexDressings[kernel_type][num_P][2][j-1];
            F1+=conj(z_i-coordin)/norm(z_i-coordin)*dz_i*Memory->VertexDressings[kernel_type][num_P][num_amp+3][j-1];
            N+=conj(z_i-coordin)/norm(z_i-coordin)*dz_i;
        }
        result=(F1/N);
        return result;
    }


    t_cmplx C_Kernel_RL_PS::SetInterpolation(t_cmplx vertex_momenta, t_cmplx prop_momenta){
		double real_prop_momenta=real(prop_momenta);
		t_cmplx P2,P1,C1,C2,L,res;
		int up,down,end;
		bool interpol_exit=false;
		end = Memory->VertexDressings[1].size()-1;

		if(real_prop_momenta < real(Memory->VertexDressings[1][0][0][0])){
			res = (VertexDressingAt(1,0,0,vertex_momenta));
		} else {
			int i = 0;
			if(real_prop_momenta > real(Memory->VertexDressings[1][end][0][0])){
				res =0.0;
			} else {
			while(!interpol_exit){
				if(real_prop_momenta < real(Memory->VertexDressings[1][i][0][0])) {
					P2=Memory->VertexDressings[1][i][0][0]; up=i;
					P1=Memory->VertexDressings[1][i-1][0][0]; down=i-1;
					L=P2-P1; C1=real_prop_momenta-P1; C2=L-C1;
					interpol_exit=true;
				}
				i++;
			}
			res = VertexDressingAt(1,up,0,vertex_momenta)*C1/L + VertexDressingAt(1,down,0,vertex_momenta)*C2/L;
			}
		}
		return res;
	}

	t_cmplxTensor C_Kernel_RL_PS::SetPSMatrix_etta_quark(t_cmplx& vertex_momenta, t_cmplx& prop_momenta){
		t_cmplx PS_factor;
		t_cmplxTensor PS_Matrix(0);
		if (real(vertex_momenta)<1000.0){
			PS_factor=SetInterpolation(vertex_momenta,prop_momenta);
		}
		else PS_factor=0.0;
		PS_Matrix=-1.0*PS_factor*Z2/(1.0 + (real(vertex_momenta) + 13.0)/100.0);
		return PS_Matrix;
	}

	t_cmplxTensor C_Kernel_RL_PS::SetPSMatrix_etta_BSE(t_cmplx& vertex_momenta, t_cmplx& prop_momenta){
		t_cmplx PS_factor;
		t_cmplxTensor PS_Matrix(0);
		if (real(vertex_momenta)<1000.0){
			PS_factor=2.0*real(SetInterpolation(vertex_momenta,prop_momenta));
		}
		else PS_factor=0.0;
		PS_Matrix=1.0*Z2*PS_factor/(1.0 + (real(prop_momenta))/100.0);
		return PS_Matrix;
	}

	t_cmplxTensor C_Kernel_RL_PS::SetPSMatrix_pion_quark(t_cmplx& vertex_momenta, t_cmplx& prop_momenta){
		t_cmplx PS_factor;
		double PionDecayConst = 0.093;
		t_cmplxTensor PS_Matrix(0);
		//PS_factor=GetDressingAt(0,0,0,vertex_momenta);
		PS_factor= Propagators[1]->PropagatorAtPoint(vertex_momenta)[1];
		PS_Matrix=-3.0*PS_factor*Z2/(prop_momenta + PseudoMesonMass*PseudoMesonMass)/PionDecayConst/(1.0 + (real(vertex_momenta))/100.0);
		return PS_Matrix;
	}

	t_cmplxTensor C_Kernel_RL_PS::SetPSMatrix_pion_BSE(t_cmplx& vertex_momenta, t_cmplx& prop_momenta){
		t_cmplx PS_factor;
		double PionDecayConst = 0.093;
		t_cmplxTensor PS_Matrix(0);
		//PS_factor=GetDressingAt(0,0,0,vertex_momenta);
		PS_factor= Propagators[1]->PropagatorAtPoint(vertex_momenta)[1];
		PS_Matrix=3.0*real(PS_factor)*Z2/(prop_momenta + PseudoMesonMass*PseudoMesonMass)/PionDecayConst/(1.0 + (real(vertex_momenta))/100.0);
		return PS_Matrix;
	}

	void C_Kernel_RL_PS::setConvolutionType(int type){
		if (Exchange_type_ID==Pion_exchange_ID) {
			if (type == 0) SetPSMatrix=&C_Kernel_RL_PS::SetPSMatrix_pion_quark;
			if (type == 1) SetPSMatrix=&C_Kernel_RL_PS::SetPSMatrix_pion_BSE;
		}

		if (Exchange_type_ID==Etta_exchange_ID) {
			if (type == 0) SetPSMatrix=&C_Kernel_RL_PS::SetPSMatrix_etta_quark;
			if (type == 1) SetPSMatrix=&C_Kernel_RL_PS::SetPSMatrix_etta_BSE;
		}
	}




