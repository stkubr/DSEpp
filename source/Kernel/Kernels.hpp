#pragma once

#include "AbstractKernel.hpp"

class C_Kernel_RL: public C_AbstractKernel{
	public:
	
	C_Kernel_RL(){
		SetNameID("Kernel RainbowLadder",1);
		Memory->VertexDressings.resize(2, t_cmplxArray3D(1));
		//Memory->VertexDressings.resize(1);
		//Memory->VertexDressings[0].resize(1);
	}
	
    void info() { std::cout << "Kernel RainbowLadder" << std::endl; } 
    
    void SetKmatrix(t_cmplxMatrix2D (*K_matrix), std::vector<t_cmplxTensor> (*Mediators)){
		for (int t = 0; t < 4; t++){
			for (int s = 0; s < 4; s++){
				for (int r = 0; r < 4; r++){
					for (int u = 0; u < 4; u++){
						t_cmplx Gluon_contrib=0.0;
						for (int i = 0; i < 4; i++){
							for (int j = 0; j < 4; j++){
								Gluon_contrib+=(Y.Element(t,s).Element(i) * (*Mediators)[0].Element(i,j) * Y.Element(r,u).Element(j) );
							}
						}
						(*K_matrix)(t,s)(r,u)=Gluon_contrib;
					}
				}
			}
		}
	}
	
	void SetMediators(t_cmplxVector *k, t_cmplxVector *p, t_cmplxVector *P, std::vector<t_cmplxTensor> (*Mediators)){
		t_cmplx k2_product;
		k2_product=(*k)*(*k);
		t_cmplx Gluon_factor=1.0;
		
		t_cmplxTensor Gluon_Matrix(2); 
		
		Gluon_factor=Z2*Z2*4.0/3.0*Gluon->GetGluonAt(&k2_product);
		Gluon_Matrix=Gluon_factor*(g-((*k)%(*k))/(k2_product));
		
		(*Mediators).resize(1);
		(*Mediators)[0]=Gluon_Matrix;
	}    
};

class C_Kernel_RL_PS: public C_AbstractKernel{
	protected:
	t_cmplxTensor (C_Kernel_RL_PS::*SetPSMatrix) (t_cmplx,t_cmplx);
	
	public:
	
	C_Kernel_RL_PS(){
		SetNameID("Kernel RainbowLadder + PseudoScalar exchange",1);
		Memory->VertexDressings.resize(2, t_cmplxArray3D(1));
		//Memory->VertexDressings[0].resize(1);
	}
	
    void info() { std::cout << "Kernel RainbowLadder + PseudoScalar exchange" << std::endl; }     
    
    void SetKmatrix(t_cmplxMatrix2D (*K_matrix), std::vector<t_cmplxTensor> (*Mediators)){
		for (int t = 0; t < 4; t++){
			for (int s = 0; s < 4; s++){
				for (int r = 0; r < 4; r++){
					for (int u = 0; u < 4; u++){
						t_cmplx Gluon_contrib=0.0;
						t_cmplx Pion_contrib=0.0;
						for (int i = 0; i < 4; i++){
							for (int j = 0; j < 4; j++){
								Gluon_contrib+=(Y.Element(t,s).Element(i) * (*Mediators)[0].Element(i,j) * Y.Element(r,u).Element(j) );
							}
						}
						Pion_contrib=(Y5.Element(t,s)*Y5.Element(r,u))*(*Mediators)[1];
						(*K_matrix)(t,s)(r,u)=Gluon_contrib + Pion_contrib;
					}
				}
			}
		}
	}
	
	void SetMediators(t_cmplxVector *k, t_cmplxVector *p, t_cmplxVector *P, std::vector<t_cmplxTensor> (*Mediators)){
		t_cmplx k2_product;
		k2_product=(*k)*(*k);
		t_cmplx Gluon_factor=1.0;
		t_cmplx Pion_factor=1.0;
		Pion_switcher=1.0;
		t_cmplxTensor Gluon_Matrix(2),PS_Matrix(0); 
		
		Gluon_factor=Z2*Z2*4.0/3.0*Gluon->GetGluonAt(&k2_product);
		Gluon_Matrix=Gluon_factor*(g-((*k)%(*k))/(k2_product));
		
		t_cmplx PS_momenta=(*P)*(*P);
		

		
		PS_Matrix=(this->*SetPSMatrix)(PS_momenta,k2_product);

		//Pion_factor=1.0/(k2_product + 10.0);//(pion_momenta+10.0);//*2.0*real(GetDressingAt(1,0,0,pion_momenta));
		
		//Pion_Matrix=10.0*Z2/(k2_product + 9.0)/(1.0 + real(pion_momenta)/50.0);//(k2_product + PseudoMesonMass*PseudoMesonMass);
		
		(*Mediators).resize(2);
		(*Mediators)[0]=Gluon_Matrix;
		(*Mediators)[1]=PS_Matrix;
	}
    
    t_cmplx SetInterpolation(t_cmplx vertex_momenta, t_cmplx prop_momenta){
		double re_prop_mom=real(prop_momenta);
		t_cmplx P2,P1,C1,C2,L,res;
		int up,down,end;
		bool interpol_exit=false;
		
		end = Memory->VertexDressings[1].size()-1;
		//std::cout << end << "  " << real(Memory->VertexDressings[1][end][0][0]) << std::endl;
		//cin.get();
				
		if(re_prop_mom < real(Memory->VertexDressings[1][0][0][0])){ res = (GetDressingAt(1,0,0,vertex_momenta));}
		else{	
			int i = 0;
			if(re_prop_mom > real(Memory->VertexDressings[1][end][0][0])){res =0.0;}// 2.0*real(GetDressingAt(1,end,0,vertex_momenta));}
			else {
			while(!interpol_exit){
				if(re_prop_mom < real(Memory->VertexDressings[1][i][0][0])) {
					P2=Memory->VertexDressings[1][i][0][0]; up=i;
					P1=Memory->VertexDressings[1][i-1][0][0]; down=i-1;
					L=P2-P1; C1=re_prop_mom-P1; C2=L-C1;
					interpol_exit=true;
				}
				i++;
			}
			res = GetDressingAt(1,up,0,vertex_momenta)*C1/L + GetDressingAt(1,down,0,vertex_momenta)*C2/L;
			}//std::cout << real(Memory->VertexDressings[1][up][0][0]) << "  " << C1/L << "  " << real(Memory->VertexDressings[1][down][0][0]) << "  " << C2/L << std::endl;
		}
		//std::cout << real(res) << "  " << real(prop_momenta) << std::endl;
		return res;
	}
	
	t_cmplxTensor SetPSMatrix_etta_quark(t_cmplx vertex_momenta, t_cmplx prop_momenta){
		t_cmplx PS_factor;

		t_cmplxTensor PS_Matrix(0);
		if (real(vertex_momenta)<1000.0){
			PS_factor=SetInterpolation(vertex_momenta,prop_momenta);
		}
		else PS_factor=0.0;
						
		PS_Matrix=-1.0*PS_factor*Z2/(1.0 + (real(vertex_momenta) + 13.0)/10.0);
		return PS_Matrix;
	}
	
	t_cmplxTensor SetPSMatrix_etta_BSE(t_cmplx vertex_momenta, t_cmplx prop_momenta){
		t_cmplx PS_factor;
		t_cmplxTensor PS_Matrix(0);
		if (real(vertex_momenta)<1000.0){
			PS_factor=2.0*real(SetInterpolation(vertex_momenta,prop_momenta));
		}
		else PS_factor=0.0;
						
		PS_Matrix=5.0*Z2*PS_factor/(1.0 + (real(prop_momenta))/100.0);;
		return PS_Matrix;
	}
	
	t_cmplxTensor SetPSMatrix_pion_quark(t_cmplx vertex_momenta, t_cmplx prop_momenta){
		t_cmplx PS_factor;
		t_cmplxTensor PS_Matrix(0);
		PS_factor=GetDressingAt(0,0,0,vertex_momenta);
						
		PS_Matrix=-3.0*PS_factor*Z2/(prop_momenta + PseudoMesonMass*PseudoMesonMass)/0.093/(1.0 + (real(vertex_momenta))/100.0);
		return PS_Matrix;
	}
	
	t_cmplxTensor SetPSMatrix_pion_BSE(t_cmplx vertex_momenta, t_cmplx prop_momenta){
		t_cmplx PS_factor;
		t_cmplxTensor PS_Matrix(0);
		PS_factor=GetDressingAt(0,0,0,vertex_momenta);
						
		PS_Matrix=3.0*real(PS_factor)*Z2/(prop_momenta + PseudoMesonMass*PseudoMesonMass)/0.093/(1.0 + (real(vertex_momenta))/100.0);
		return PS_Matrix;
	}
	
	void SetConvolutionType(int type){
		if (Exchange_type_ID==Pion_exchange_ID) {
			if (type == 0) SetPSMatrix=&C_Kernel_RL_PS::SetPSMatrix_pion_quark;
			if (type == 1) SetPSMatrix=&C_Kernel_RL_PS::SetPSMatrix_pion_BSE;
		}
		
		if (Exchange_type_ID==Etta_exchange_ID) {
			if (type == 0) SetPSMatrix=&C_Kernel_RL_PS::SetPSMatrix_etta_quark;
			if (type == 1) SetPSMatrix=&C_Kernel_RL_PS::SetPSMatrix_etta_BSE;
		} 
	}
};

C_AbstractKernel* C_AbstractKernel::createKernel(Kernel_ID * id){
	C_AbstractKernel * p;
    switch (*id)
    {
        case RL_ID:
            p = new C_Kernel_RL();
            p->Kernel_type_ID=(*id);
            break;  
        case RL_PS_ID:
            p = new C_Kernel_RL_PS(); 
            p->Kernel_type_ID=(*id);         
            break;              
        default:
            assert( false);
    }
    return p;
}
 
class C_Kernel_Factory{
	public:
	C_AbstractKernel* Create(Kernel_ID * _id) {
		return C_AbstractKernel::createKernel( _id );
    }
	~C_Kernel_Factory() {}
};

C_Kernel_Factory * KernelFactory = new C_Kernel_Factory;

