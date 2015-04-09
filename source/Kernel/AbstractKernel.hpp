/*
 * AbstractKernel.hpp
 *      Author: stkubr
 */

#ifndef ABSTRACTKERNEL_HPP_
#define ABSTRACTKERNEL_HPP_

/**
 * Enumeration of all available types of two-body scattering kernels
 * RL - Rainbow-ladder single gluon exchange
 * RL_PS - Rainbow-ladder single gluon exchange + PseudoScalar exchange (aka pion cloud)
*/
enum Kernel_ID {RL_ID=0, RL_PS_ID, Kernel_ID_End};

/**
 * Enumeration of all available types PseudoScalar exchange
 * Pion_exchange - the pion cloud effect
 * Etta_exchange - /f$ \eta_c /f$ exchange for charmonium
*/
enum PS_type_ID {Pion_exchange_ID=0, Etta_exchange_ID, PS_type_ID_End};

#include "../Abs/AbsDiagram.hpp"
class C_Propagator;
#include "../DSE/Propagator.hpp"
#include "../DedicMem/MemoryFactories.hpp"

class C_AbstractKernel: public C_AbsDiagram{
protected:

	/// vector of pointer to prapagator as two-body scattering mediators
	std::vector<C_Propagator*> Propagators;
	/// TODO to include vector of general vertexes

	std::vector<t_cmplxMatrix2D> threadloc_KMatrix;

	Kernel_ID Kernel_type_ID;
	PS_type_ID Exchange_type_ID;

	C_AbstractKernel(){
		SetNameID("Kernel",1);
		Memory=static_cast<C_DedicMem_Kernel*>(DedicMemFactory_Kernel->CreateMemory());
		threadloc_KMatrix.resize(omp_get_max_threads());
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
								    * threadloc_KMatrix[omp_get_thread_num()](t,s)(r,u);
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
		resizeKmatrix(threadloc_KMatrix[omp_get_thread_num()]);
		setMediators(k,p,P,MediatorKernel);
		setKmatrix(threadloc_KMatrix[omp_get_thread_num()],MediatorKernel);
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

#endif /* ABSTRACTKERNEL_HPP_ */
