#pragma once
#include "BSE_Base.h"

enum Dirac_ID {PseudoScalar_ID=0, Scalar_ID, Vector_ID, AxialVector_ID, Tensor_ID, Dirac_ID_End};

class C_BSE_Hadron_Meson: public C_BSE_Hadron_Base {
	public:

	C_BSE_Hadron_Meson(){
		//Gluon=C_Gluon::getInstance();
	}
	
	t_cmplxMatrix setResultBSA(){
		t_cmplxMatrix result;
		
		SetWaveFunctions();
		result=(this->*GetBSA_ref)();
		return result;
	}
	
	virtual C_BSE_Hadron_Meson * MakeCopy()
	{
		return new C_BSE_Hadron_Meson(*this);
	}
	
	void SetWaveFunctions(){
		if(!flag_precalculation){
		SetDiracStructures(Momenta.k,Momenta.P,&Amplitudes);
		setPropagators(&Momenta.k_p, &Momenta.k_m);
		for (int i = 0; i < num_amplitudes; i++){ WaveFunctions[i]=S_p*Amplitudes[i]*S_m;}
		}
		else{
			for (int i = 0; i < num_amplitudes; i++){ WaveFunctions[i]=Memory->AmpStorage[i][Int_counter];}
		}
	}
	
	void SetFullWaveFunction(){
		FullWaveFunction.Zero();
		for (int i = 0; i < num_amplitudes; i++){ FullWaveFunction+=WaveFunctions[i]*U_amp[i]; }
	}
	
	void SetFullAmplitude(){
		FullAmplitude=Amplitudes[0];
		FullAmplitude.Zero();
		for (int i = 0; i < num_amplitudes; i++){ FullAmplitude+=Amplitudes[i]*U_amp[i]; }
	}
	
	t_cmplxMatrix GetBSA(){
		t_cmplxMatrix result,pre_result;
		t_cmplxVector k_p_P;
		k_p_P=(Momenta.k + Momenta.p - Momenta.P)/2.0;
		result.Resize(num_amplitudes,1);
		pre_result.Resize(num_amplitudes,1);
		SetFullWaveFunction();
		bool flag_reset_kernel=true;
		for (int i = 0; i < num_amplitudes; i++) {pre_result(i,0)=Kernel->TraceKernelWithoutStoring(&(Projectors[i]),&FullWaveFunction,&Momenta.q,&Momenta.k,&k_p_P,flag_reset_kernel); flag_reset_kernel=false;}
		result=DisentangleAmps(&pre_result);
		return result;
	}
	
	t_cmplxMatrix GetBSA_matrix(){
		t_cmplxMatrix result_M,result,pre_result;
		t_cmplxVector k_p_P;
		k_p_P=(Momenta.k + Momenta.p - Momenta.P)/2.0;
		result.Resize(num_amplitudes,1);
		pre_result.Resize(num_amplitudes,1);
		result_M.Resize(num_amplitudes,num_amplitudes);
		bool flag_reset_kernel=true;
		for (int i = 0; i < num_amplitudes; i++){
			for (int j = 0; j < num_amplitudes; j++){
				pre_result(j,0)=Kernel->TraceKernelWithoutStoring(&(Projectors[j]),&(WaveFunctions[i]),&Momenta.q,&Momenta.k,&k_p_P,flag_reset_kernel); 
				flag_reset_kernel=false;
			}
			result=DisentangleAmps(&pre_result);
			for (int j = 0; j < num_amplitudes; j++){
				result_M(i,j)=(result(j,0));
			}
		}
		return result_M;
	}
	
	t_cmplxMatrix GetBSA_norm(){
		double Norm_factor;
		t_cmplxMatrix result;
		result.Resize(num_amplitudes,1);
		
		SetDiracStructures(Momenta.k,Momenta.P,&Amplitudes);
		SetFullAmplitude();
		setPropagators(&Momenta.k_p, &Momenta.k_m);
		
		Norm_factor=pi/2.0;
		
		result(0,0)=Norm_factor*(FullAmplitude*S_p*FullAmplitude*S_m).Tr();
		return result;
	}
	
	virtual t_cmplxMatrix DisentangleAmps(t_cmplxMatrix * pre_result){t_cmplxMatrix dummy; std::cout << "Error virtual call" << std::endl; StopLine(); return dummy;}
	
	
	static C_BSE_Hadron_Meson * createMesonBSE(Dirac_ID id);
	
};

class C_Mesons_PseudoScalar: public C_BSE_Hadron_Meson {
	
	public:
	C_Mesons_PseudoScalar(){
		SetNameID("BSE_PseudoScalar",1);
		num_amplitudes=4; // number of amplitudes for PseudoScalar
		Initialization();
		Amplitudes.resize(4);
		Projectors.resize(4);
		WaveFunctions.resize(4);
		WeightCoeff.resize(4);
		SaveBSEPath="Data_files/SaveBSE_PseudoScalar.dat";
	}
	
	void SetDiracStructures(t_cmplxVector _k, t_cmplxVector _P, std::vector<t_cmplxDirac> * DiracStructure){
		#if BASIS_TYPE == 1
		(*DiracStructure)[0]=ii*Y5;
		(*DiracStructure)[1]=Y5*(_P*Z);
		(*DiracStructure)[2]=Y5*(_k*Z)*(_k*_P);
		(*DiracStructure)[3]=Y5*((_k*SIG)*_P);
		
		#endif
		
		#if BASIS_TYPE == 2
		
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
		
		#endif
	}
	
	virtual C_Mesons_PseudoScalar * MakeCopy()
	{
		return new C_Mesons_PseudoScalar(*this);
	}
	
	t_cmplxMatrix DisentangleAmps(t_cmplxMatrix * pre_result){
		#if BASIS_TYPE == 1
		t_cmplxMatrix result(4,1);
		result(0,0)=1.0/WeightCoeff[0]*(*pre_result)(0,0);
		result(1,0)=1.0/WeightCoeff[0]*(Momenta.p2/Momenta.N2_Factor*(*pre_result)(1,0) - 1.0/Momenta.N2_Factor*(*pre_result)(2,0));
		result(2,0)=1.0/WeightCoeff[0]*(-1.0/Momenta.N2_Factor*(*pre_result)(1,0) + Momenta.P2/Momenta.p_P/Momenta.p_P/Momenta.N2_Factor*(*pre_result)(2,0));
		result(3,0)=1.0/WeightCoeff[0]*(-1.0/Momenta.N2_Factor*(*pre_result)(3,0));
		return result;
		#endif
		
		#if BASIS_TYPE == 2
		
		t_cmplxMatrix result(num_amplitudes,1);
		for (int i = 0; i < num_amplitudes; i++) result(i,0)=(*pre_result)(i,0)/WeightCoeff[i];
		return result;

		#endif
	}
};

class C_Mesons_Scalar: public C_BSE_Hadron_Meson {
	
	public:
	C_Mesons_Scalar(){
		SetNameID("BSE_Scalar",1);
		num_amplitudes=4; // number of amplitudes for PseudoScalar
		Initialization();
		Amplitudes.resize(4);
		Projectors.resize(4);
		WaveFunctions.resize(4);
		WeightCoeff.resize(4);
		SaveBSEPath="Data_files/SaveBSE_Scalar.dat";
	}
	
	void SetDiracStructures(t_cmplxVector _k, t_cmplxVector _P, std::vector<t_cmplxDirac> * DiracStructure){
		#if BASIS_TYPE == 1
		(*DiracStructure)[0]=ii*I;
		(*DiracStructure)[1]=I*(_P*Z)*(_k*_P);
		(*DiracStructure)[2]=I*(_k*Z);
		(*DiracStructure)[3]=I*((_k*SIG)*_P);
		
		#else 
		
		Y_T=Momenta.TransIn(Y,_P);
		t_cmplxVector _k_T;
		_k_T=Momenta.TransIn(_k,_P);
		t_cmplx _k2_T,_k2, _k_P,_P2;
		_k2_T=_k_T*_k_T;
		_k2=_k*_k;
		_k_P=_k*_P;
		_P2=_P*_P;
		
		(*DiracStructure)[0]=I;
		(*DiracStructure)[1]=(_P*Z);
		(*DiracStructure)[2]=(_k_T*Z);
		(*DiracStructure)[3]=ii/2.0*((_k_T*Z)*(_P*Z) - (_P*Z)*(_k_T*Z));
		
		#endif
	}
	
	virtual C_Mesons_Scalar * MakeCopy()
	{
		return new C_Mesons_Scalar(*this);
	}
	
	t_cmplxMatrix DisentangleAmps(t_cmplxMatrix * pre_result){
		#if BASIS_TYPE == 1
		t_cmplxMatrix result(4,1);
		result(0,0)=1.0/WeightCoeff[0]*(*pre_result)(0,0);
		result(1,0)=-1.0/WeightCoeff[0]*(Momenta.p2/Momenta.p_P/Momenta.p_P/Momenta.N2_Factor*(*pre_result)(1,0) - 1.0/Momenta.N2_Factor*(*pre_result)(2,0));
		result(2,0)=-1.0/WeightCoeff[0]*(-1.0/Momenta.N2_Factor*(*pre_result)(1,0) + Momenta.P2/Momenta.N2_Factor*(*pre_result)(2,0));
		result(3,0)=-1.0/WeightCoeff[0]/Momenta.N2_Factor*(*pre_result)(3,0);
		return result;
		
		#endif
		
		#if BASIS_TYPE == 2
		
		t_cmplxMatrix result(num_amplitudes,1);
		for (int i = 0; i < num_amplitudes; i++) result(i,0)=(*pre_result)(i,0)/WeightCoeff[i];
		return result;

		#endif
	}
};

class C_Mesons_Vector: public C_BSE_Hadron_Meson {
	
	public:
	C_Mesons_Vector(){
		SetNameID("BSE_Vector",1);
		num_amplitudes=5; // number of amplitudes for PseudoScalar
		Initialization();
		Amplitudes.resize(8);
		Projectors.resize(8);
		WaveFunctions.resize(8);
		WeightCoeff.resize(8);
		SaveBSEPath="Data_files/SaveBSE_Vector.dat";
	}
	
	void SetDiracStructures(t_cmplxVector _k, t_cmplxVector _P, std::vector<t_cmplxDirac> * DiracStructure){
		Y_T=Momenta.TransIn(Y,_P);
		t_cmplxVector _k_T;
		_k_T=Momenta.TransIn(_k,_P);
		t_cmplx _k2_T,_k2, _k_P,_P2;
		_k2_T=_k_T*_k_T;
		_k2=_k*_k;
		_k_P=_k*_P;
		_P2=_P*_P;
		t_cmplxDirac Y_T_T;
		Y_T_T=Y_T - _k_T*(_k_T*Z)/_k2_T;

		t_cmplx k2_sqrt;
		k2_sqrt=(_k2);
		
		#if BASIS_TYPE == 1
		(*DiracStructure)[0]=Y_T;
		(*DiracStructure)[1]=(_k_T*(_k_T*Z) - 1.0/3.0*Y_T*(_k2_T))/_k2;
		(*DiracStructure)[2]=_k_T*(_P*Z)*_k_P/(_k2*_P2);
		(*DiracStructure)[3]=-1.0*(Y_T*( (_P*Z)*(_k*Z) - (_k*Z)*(_P*Z) ) + 2.0*_k_T*(_P*Z) )/k2_sqrt/2.0;
		(*DiracStructure)[4]=ii*_k_T*I/k2_sqrt;
		(*DiracStructure)[5]=ii*(Y_T*(_k_T*Z) - (_k_T*Z)*Y_T)*_k_P/_k2;
		(*DiracStructure)[6]=ii*(Y_T*(_P*Z) - (_P*Z)*Y_T)*(1.0 - _k_P*_k_P/_k2/_P2) - ii*2.0*_k_T*(_k_T*Z)*(_P*Z)/_k2;
		(*DiracStructure)[7]=ii*_k_T*(_k_T*Z)*(_P*Z)/_k2;
		
		#endif
		
		#if BASIS_TYPE == 2

		(*DiracStructure)[0]=_k_T*I;
		(*DiracStructure)[1]=(_P*Z)*_k_T;
		(*DiracStructure)[2]=(_k_T*Z)*_k_T;
		(*DiracStructure)[3]=ii/2.0*((_k_T*Z)*(_P*Z) - (_P*Z)*(_k_T*Z))*_k_T;
		(*DiracStructure)[4]=Y_T_T;
		(*DiracStructure)[5]=(_P*Z)*Y_T_T;
		(*DiracStructure)[6]=(_k_T*Z)*Y_T_T;
		(*DiracStructure)[7]=ii/2.0*((_k_T*Z)*(_P*Z) - (_P*Z)*(_k_T*Z))*Y_T_T;

		#endif
	}
	
	virtual C_Mesons_Vector * MakeCopy()
	{
		return new C_Mesons_Vector(*this);
	}
	
	t_cmplxMatrix DisentangleAmps(t_cmplxMatrix * pre_result){
		t_cmplxMatrix result(num_amplitudes,1);
		for (int i = 0; i < num_amplitudes; i++) result(i,0)=(*pre_result)(i,0)/WeightCoeff[i];
		return result;
	}
};


class C_Mesons_AxialVector: public C_BSE_Hadron_Meson {
	
	public:
	C_Mesons_AxialVector(){
		SetNameID("BSE_AxialVector",1);
		num_amplitudes=8; // number of amplitudes for PseudoScalar
		Initialization();
		Amplitudes.resize(8);
		Projectors.resize(8);
		WaveFunctions.resize(8);
		WeightCoeff.resize(8);
		SaveBSEPath="Data_files/SaveBSE_AxialVector.dat";
	}
	
	void SetDiracStructures(t_cmplxVector _k, t_cmplxVector _P, std::vector<t_cmplxDirac> * DiracStructure){
		Y_T=Momenta.TransIn(Y,_P);
		t_cmplxVector _k_T;
		_k_T=Momenta.TransIn(_k,_P);
		t_cmplx _k2_T,_k2, _k_P,_P2;
		_k2_T=_k_T*_k_T;
		_k2=_k*_k;
		_k_P=_k*_P;
		_P2=_P*_P;
		
		t_cmplxDirac Y_T_T;
		Y_T_T=Y_T - _k_T*(_k_T*Z)/_k2_T;

		t_cmplx k2_sqrt;
		k2_sqrt=(sqrt(_k2));
		
		#if BASIS_TYPE == 1
		(*DiracStructure)[0]=Y5*Y_T*_k_P;
		(*DiracStructure)[7]=Y5*ii*_k_T*I/k2_sqrt*_k_P;
		(*DiracStructure)[2]=Y5*(_k_T*(_k_T*Z) - 1.0/3.0*Y_T*(_k2_T))/_k2;
		(*DiracStructure)[3]=Y5*_k_T*(_P*Z)*_k_P/(_k2*_P2);
		(*DiracStructure)[4]=-1.0*Y5*(Y_T*( (_P*Z)*(_k*Z) - (_k*Z)*(_P*Z) ) + 2.0*_k_T*(_P*Z) )/k2_sqrt/2.0;
		(*DiracStructure)[5]=Y5*ii*(Y_T*(_k_T*Z) - (_k_T*Z)*Y_T)*_k_P/_k2;
		(*DiracStructure)[6]=Y5*(ii*(Y_T*(_P*Z) - (_P*Z)*Y_T)*(1.0 - _k_P*_k_P/_k2/_P2) - ii*2.0*_k_T*(_k_T*Z)*(_P*Z)/_k2);
		(*DiracStructure)[1]=Y5*ii*_k_T*(_k_T*Z)*(_P*Z)/_k2;
		
		#endif
		
		#if BASIS_TYPE == 2
		
		(*DiracStructure)[0]=Y5*_k_T;
		(*DiracStructure)[1]=Y5*(_P*Z)*_k_T;
		(*DiracStructure)[2]=Y5*(_k_T*Z)*_k_T;
		(*DiracStructure)[3]=ii/2.0*Y5*((_k_T*Z)*(_P*Z) - (_P*Z)*(_k_T*Z))*_k_T;
		(*DiracStructure)[4]=Y5*Y_T_T;
		(*DiracStructure)[5]=Y5*(_P*Z)*Y_T_T;
		(*DiracStructure)[6]=Y5*(_k_T*Z)*Y_T_T;
		(*DiracStructure)[7]=ii/2.0*Y5*((_k_T*Z)*(_P*Z) - (_P*Z)*(_k_T*Z))*Y_T_T;
		
		#endif
	}
	
	virtual C_Mesons_AxialVector * MakeCopy()
	{
		return new C_Mesons_AxialVector(*this);
	}
	
	t_cmplxMatrix DisentangleAmps(t_cmplxMatrix * pre_result){
		t_cmplxMatrix result(num_amplitudes,1);
		for (int i = 0; i < num_amplitudes; i++) result(i,0)=(*pre_result)(i,0)/WeightCoeff[i];
		return result;
	}
};


class C_Mesons_Tensor: public C_BSE_Hadron_Meson {
	
	public:
	C_Mesons_Tensor(){
		SetNameID("BSE_Tensor",1);
		num_amplitudes=8; // number of amplitudes for PseudoScalar
		Initialization();
		Amplitudes.resize(num_amplitudes);
		Projectors.resize(num_amplitudes);
		WaveFunctions.resize(num_amplitudes);
		WeightCoeff.resize(num_amplitudes);
		SaveBSEPath="Data_files/SaveBSE_Tensor.dat";
	}
	
	void SetDiracStructures(t_cmplxVector _k, t_cmplxVector _P, std::vector<t_cmplxDirac> * DiracStructure){
		Y_T=Momenta.TransIn(Y,_P);
		t_cmplxVector _k_T;
		_k_T=Momenta.TransIn(_k,_P);
		t_cmplx _k2_T,_k2, _k_P,_P2;
		_k2_T=_k_T*_k_T;
		_k2=_k*_k;
		_k_P=_k*_P;
		_P2=_P*_P;
		
		t_cmplxDirac Y_T_T;
		Y_T_T=Y_T - _k_T*(_k_T*Z)/_k2_T;
		
		t_cmplxDirac M,N,M_T,N_T;
		
		t_cmplxTensor g_T;
		
		M = (Y_T_T % _k_T) + (_k_T % Y_T_T);
 		N = (_k_T%_k_T)*I;
		g_T = g - (_P%_P)/_P2;
		
		t_cmplx temp1,temp2;
		//temp1=(M|g_T);
		//std::cout << (M|g_T) << std::endl;
		//temp2=(N|g_T);
		
		temp1=(g_T).Mag2();
		//std::cout << temp1 << std::endl;
		M_T=M - g_T*(M|g_T)/temp1;
		N_T=N - g_T*(M|g_T)/temp1;
		
		
		//DebugLine("after");
		(*DiracStructure)[0]=M_T;
		(*DiracStructure)[1]=(_P*Z)*M_T;
		(*DiracStructure)[2]=(_k_T*Z)*M_T;
		(*DiracStructure)[3]=ii/2.0*((_k_T*Z)*(_P*Z) - (_P*Z)*(_k_T*Z))*M_T;
		
		(*DiracStructure)[4]=N_T;
		(*DiracStructure)[5]=(_P*Z)*N_T;
		(*DiracStructure)[6]=(_k_T*Z)*N_T;
		(*DiracStructure)[7]=ii/2.0*((_k_T*Z)*(_P*Z) - (_P*Z)*(_k_T*Z))*N_T;
		
	}
	
	virtual C_Mesons_Tensor * MakeCopy()
	{
		return new C_Mesons_Tensor(*this);
	}
	
	t_cmplxMatrix DisentangleAmps(t_cmplxMatrix * pre_result){
		t_cmplxMatrix result(num_amplitudes,1);
		for (int i = 0; i < num_amplitudes; i++) result(i,0)=(*pre_result)(i,0)/WeightCoeff[i];
		return result;
	}
};

C_BSE_Hadron_Meson * C_BSE_Hadron_Meson::createMesonBSE(Dirac_ID id){
	C_BSE_Hadron_Meson * p;
	switch(id){
		case PseudoScalar_ID:
			p = new C_Mesons_PseudoScalar();          
			break;  
		case Scalar_ID:
			p = new C_Mesons_Scalar();          
			break; 
		case Vector_ID:
			p = new C_Mesons_Vector();          
			break; 
		case AxialVector_ID:
			p = new C_Mesons_AxialVector();          
			break; 
		case Tensor_ID:
			p = new C_Mesons_Tensor();          
			break;
        default:
            assert(false);
	}
	return p;
}

