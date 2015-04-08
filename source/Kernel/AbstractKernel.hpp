#pragma once

enum Kernel_ID {RL_ID=0, RL_PS_ID, Kernel_ID_End};

enum PS_type_ID {Pion_exchange_ID=0, Etta_exchange_ID, PS_type_ID_End};

#include "../Abs/AbsDiagram.hpp"
class C_Propagator;
#include "../DSE/Propagator.hpp"



class C_AbstractKernel: public C_AbsDiagram{

protected:
	std::vector<C_Propagator*> Propagators;
	/// TODO to include vector of general vertexes
	std::vector<t_cmplxMatrix2D> KMatrixThreadStorage;
	Kernel_ID Kernel_type_ID;
	PS_type_ID Exchange_type_ID;

	C_AbstractKernel(){
		SetNameID("Kernel",1);	
		Memory=static_pointer_cast<C_DedicMem_Kernel*>(DedicMemFactory_Kernel->CreateMemory());
		KMatrixThreadStorage.resize(omp_get_max_threads());
	}

public:
	C_DedicMem_Kernel * Memory;

	virtual void info()=0;

	void setExchangeID(PS_type_ID exchange_id){
		Exchange_type_ID=exchange_id;
	}

	Kernel_ID KernelID(){
		return Kernel_type_ID;
	}

	void setPropagators(std::vector<C_Propagator*> __Propagators){
		Propagators=__Propagators;
	}

	virtual void setConvolutionType(int type){}

	static C_AbstractKernel* createKernel( Kernel_ID id );

	// Take a trace of Tr(Projector * Kernel * WaveFunction) at provided momenta
	t_cmplx TraceKernelWithoutStoring(t_cmplxDirac& Projector,
									  t_cmplxDirac& WaveFunc,
									  t_cmplxVector& k,
									  t_cmplxVector& p,
									  t_cmplxVector& P, bool flag_reset_kernel){

		t_cmplx result, temp_storage;
		if (flag_reset_kernel){
			setKMatrixThreadStorage(k,p,P);
		}
		int rank=Projector(0,0).Rank();

		for (int t = 0; t < 4; t++){
			for (int s = 0; s < 4; s++){
				for (int r = 0; r < 4; r++){
					for (int u = 0; u < 4; u++){
                        temp_storage = takeInnerProduct((Projector)(u,t), (WaveFunc)(s,r), rank)
								    * KMatrixThreadStorage[omp_get_thread_num()](t,s)(r,u);
                        result+=temp_storage;
					}
				}
			}
		}
		return result;
	}

	// Sets K_matrix for each thread calling the trace
	void setKMatrixThreadStorage(t_cmplxVector& k, t_cmplxVector& p, t_cmplxVector& P){
		std::vector<t_cmplxTensor> MediatorKernel;
		resizeKmatrix(KMatrixThreadStorage[omp_get_thread_num()]);
		setMediators(k,p,P,MediatorKernel);
		setKmatrix(KMatrixThreadStorage[omp_get_thread_num()],MediatorKernel);
	}

	// Takes inner product of two provided tensors
	t_cmplx takeInnerProduct(t_cmplxTensor& A, t_cmplxTensor& B, int rank){
		t_cmplx result=0.0;
		if(rank==2){
			for (int i = 0; i < 4; i++){
				for (int j = 0; j < 4; j++){
					result+=(A)(i,j)*(B)(j,i);} } }
		if(rank==1){
			for (int i = 0; i < 4; i++){
				result+=(A)(i)*(B)(i);} }
		if(rank==0){
			result+=(A)(0)*(B)(0);}
		return result;
	}

	// K matrix resizing to (4,4,4,4)
	void resizeKmatrix(t_cmplxMatrix2D& K_matrix){
		K_matrix.Resize(4,4);
		for (int i = 0; i < 4; i++){
			for (int j = 0; j < 4; j++){
				K_matrix(i,j).Resize(4,4);
			}
		}
	}

	// Sets particles propagating inside of Kernel (specified by inherited kernels classes)
	virtual void setMediators(t_cmplxVector& k, t_cmplxVector& p, t_cmplxVector& P,
			std::vector<t_cmplxTensor>& Mediators)=0;

	// Set a (t,s,r,u) element of the K_matrix (specified by inherited kernels classes)
	virtual t_cmplx ElementKmatrix(int t, int s, int r, int u,
			std::vector<t_cmplxTensor>& Mediators)=0;

	// Composes K matrix for tracer
	void setKmatrix(t_cmplxMatrix2D& K_matrix, std::vector<t_cmplxTensor>& Mediators){
		for (int t = 0; t < 4; t++){
			for (int s = 0; s < 4; s++){
				for (int r = 0; r < 4; r++){
					for (int u = 0; u < 4; u++){
						K_matrix(t,s)(r,u)=ElementKmatrix(t, s, r, u, Mediators);
					}
				}
			}
		}
	}

};


