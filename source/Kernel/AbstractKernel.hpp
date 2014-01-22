#pragma once

enum Kernel_ID {RL_ID=0, RL_PS_ID, Kernel_ID_End};

enum Gluon_ID {RL_MT_Light_ID=0, RL_MT_Heavy_ID, PS_Light_ID, RL_MT_Heavy_DD_ID, Gluon_ID_End};

enum PS_type_ID {Pion_exchange_ID=0, Etta_exchange_ID, PS_type_ID_End};

class C_AbstractKernel: public C_AbsDiagram{
	
	protected:
	t_cmplx Z2;
	C_Gluon * Gluon;
	t_cmplx PseudoMesonMass;
	std::vector<t_cmplxMatrix2D> K_Store;
	Kernel_ID Kernel_type_ID;
	PS_type_ID Exchange_type_ID;
	//static C_Kernel * p_instance;
	
	C_AbstractKernel(){
		SetNameID("Kernel",1);	
		Memory_abs=DedicMemFactory_Kernel->CreateMemory();
		Memory=(C_DedicMem_Kernel*)Memory_abs;
		Z2=1.0;
		Pion_switcher=1.0;
		K_Store.resize(omp_get_max_threads());
	}

	public:
	virtual void info()=0;
	double Pion_switcher;
	C_DedicMem_Kernel* Memory;
	
	void SpecifyGluon(Gluon_ID gluon_id){
		switch (gluon_id){
			case RL_MT_Light_ID:
				Gluon=C_Gluon::getInstance("../Parameters_files/Gluons/RL_MT_Light_List.txt");
				break;
			case PS_Light_ID:
				Gluon=C_Gluon::getInstance("../Parameters_files/Gluons/PS_Light_List.txt");
				break;
			case RL_MT_Heavy_ID:
				Gluon=C_Gluon::getInstance("../Parameters_files/Gluons/RL_MT_Heavy_List.txt");
				break;
			case RL_MT_Heavy_DD_ID:
				Gluon=C_Gluon::getInstance("../Parameters_files/Gluons/RL_MT_Heavy_DD_List.txt");
				break;
			default:
				assert( false);
		}
	}
		
	void SetMesonExchangeMass(t_cmplx _M){
		PseudoMesonMass=_M;
	}
	
	void SetExchangeID(PS_type_ID exchange_id){
		Exchange_type_ID=exchange_id;
	}
	
	Kernel_ID GetKernelID(){
		return Kernel_type_ID;
	}
	
	
	virtual void SetConvolutionType(int type){}
	
	static C_AbstractKernel* createKernel( Kernel_ID * id );

// Allocate Memory for K_Matrix_storage
	void AllocateMemFor(int _K_ctr){
		Memory->ResizeKstorage(_K_ctr);
	}
	
	void setZ2DressingFactor(t_cmplx _Z2){
		Z2=_Z2;
	}
	
	void SetDressingStorages(t_cmplxVector *k, t_cmplxVector *p, t_cmplxVector *P){
		
	}
	

	t_cmplx GetDressingAt(int kernel_type, int num_P, int num_amp, t_cmplx coordin){
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
	
// Set K storage at (k,p,q) 	
	void SetKstorage(int K_ctr,t_cmplxVector *k, t_cmplxVector *p, t_cmplxVector *q){
		t_cmplx k2_product;
		k2_product=(*k)*(*k);
		t_cmplxMatrix2D K_matrix;
		t_cmplx Gluon_factor;
		ResizeKmatrix(&K_matrix);
		t_cmplxTensor Gluon_Matrix(2);
		Gluon_factor=Gluon->GetGluonAt(&k2_product);
		Gluon_Matrix=Gluon_factor*(g-((*k)%(*k))/(k2_product));
		//SetKmatrix(&K_matrix,&Gluon_Matrix);
		Memory->SetKmatrixAt(K_ctr,&K_matrix);
	}
	
	t_cmplxMatrix2D * GetKmatrixAt(int i){
		return Memory->GetKmatrixAt(i);
	}
	
	t_cmplx TraceKernelWithoutStoring(t_cmplxDirac *Projector, t_cmplxDirac *WaveFunc, t_cmplxVector *k, t_cmplxVector *p, t_cmplxVector *P,bool flag_reset_kernel){
		t_cmplx k2_product,result,dummy;
		int rank=0;
		result=0.0;
		//k2_product=(*k)*(*k);
		//dcx Gluon_factor=1.0;
		//dcx Pion_factor=1.0;
		t_cmplxMatrix2D K_matrix;
		ResizeKmatrix(&K_matrix);
		
		if (flag_reset_kernel){
			std::vector<t_cmplxTensor> MediatorKernel;
			SetMediators(k,p,P,&MediatorKernel);
			SetKmatrix(&K_matrix,&MediatorKernel);
			//std::cout << flag_reset_kernel << "  " << omp_get_thread_num() << "  " << K_Store.size() << std::endl;
			K_Store[omp_get_thread_num()]=K_matrix;	
		}
		else {
			K_matrix=K_Store[omp_get_thread_num()];
		}

		rank=(*Projector)(0,0).Rank();
		for (int t = 0; t < 4; t++){
			for (int s = 0; s < 4; s++){
				for (int r = 0; r < 4; r++){
					for (int u = 0; u < 4; u++){
						//std::cout << t << "  " << s << "  " << r << "  " << u << std::endl;
						dummy=InnerTensorProduct(&((*Projector)(u,t)),&((*WaveFunc)(s,r)),rank)*K_matrix(t,s)(r,u);
						result+=dummy;
					}
				}
			}
		}	
		
		
		return result;
	}
	
	t_cmplx InnerTensorProduct(t_cmplxTensor *A, t_cmplxTensor *B, int rank){
		t_cmplx result=0.0;
		if(rank==2){
			for (int i = 0; i < 4; i++){
				for (int j = 0; j < 4; j++){
					result+=(*A)(i,j)*(*B)(j,i);} } }
		if(rank==1){
			for (int i = 0; i < 4; i++){
				result+=(*A)(i)*(*B)(i);} }
		if(rank==0){
			result+=(*A)(0)*(*B)(0);}
		return result;
	}
	
// K matrix resizing to (4,4,4,4)
	void ResizeKmatrix(t_cmplxMatrix2D (*K_matrix)){
		(*K_matrix).Resize(4,4);
		for (int i = 0; i < 4; i++){
			for (int j = 0; j < 4; j++){
				(*K_matrix)(i,j).Resize(4,4);
			}
		}
	}
	 
	virtual void SetMediators(t_cmplxVector *k, t_cmplxVector *p, t_cmplxVector *P, std::vector<t_cmplxTensor> (*Mediators))=0;
// Set K matrix
	virtual void SetKmatrix(t_cmplxMatrix2D (*K_matrix), std::vector<t_cmplxTensor> (*Mediators))=0;	
	
	t_cmplxMatrix Check(t_cmplxMatrix2D (*K_matrix)){
		t_cmplxMatrix temp(4,4);
		t_cmplxVector k_dummy;
		t_cmplx dummy;
		k_dummy.SetP4(1.0,1.0,1.0,1.0);
		//for (int i = 0; i < 4; i++){
		for (int t = 0; t < 4; t++){
			for (int s = 0; s < 4; s++){
				for (int r = 0; r < 4; r++){
						temp(t,s)+=(*K_matrix)(t,r)(r,s);
				}
			}
		}
		return temp;
	}
	
};


