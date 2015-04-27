//
// Created by stkubr on 09.04.15.
//

#ifndef _DSEPP_TWOQUARKKERNEL_HPP_
#define _DSEPP_TWOQUARKKERNEL_HPP_

#include "AbstractKernel.hpp"

class C_TwoQuarkKernel: public C_AbstractKernel{
protected:
    /// Kernel's Memory handler class
    C_DedicMem_Kernel * Memory;

    /// Thread-local KMatrix storage
    std::vector<t_cmplxMatrix2D> threadloc_KMatrix;


    /**
     * Pointers to Propagators that will be used in inside the Kernel;
     * It can be that during the calculation of quark propagator
     * it will be traced with the kernel that contains a pointer to the very same
     * quark propagator; this fact reflects the essence of self-consistency of
     * DSE integral equations
    */
    std::vector<C_Propagator*> Propagators;

    Kernel_ID Kernel_type_ID;

    /**
     * The constractor
    */
    C_TwoQuarkKernel(){
        Memory=static_cast<C_DedicMem_Kernel*>(DedicMemFactory_Kernel->CreateMemory());
        threadloc_KMatrix.resize(omp_get_max_threads());
    }

    /**
     * Sets K_matrix for each thread calling the trace
    */
    void setKMatrixThreadStorage(t_cmplxVector& k, t_cmplxVector& p, t_cmplxVector& P){
        std::vector<t_cmplxTensor> MediatorKernel;
        resizeKmatrix(threadloc_KMatrix[omp_get_thread_num()]);
        setMediators(k,p,P,MediatorKernel);
        setKmatrix(threadloc_KMatrix[omp_get_thread_num()],MediatorKernel);
    }

    /**
     * Takes inner product of two provided tensors of the given rank
    */
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

    /**
     * Resize the K matrix to (4,4,4,4)
    */
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

public:
    void setPropagators(std::vector<C_Propagator*> __Propagators){
        Propagators=__Propagators;
    }

    t_cmplx TraceKernelWithoutStoring(t_cmplxDirac &Projector,
            t_cmplxDirac &WaveFunc,
            t_cmplxVector &k,
            t_cmplxVector &p,
            t_cmplxVector &P,
            bool flag_reset_kernel) {

        t_cmplx result, temp_storage;
        if (flag_reset_kernel) {
            setKMatrixThreadStorage(k, p, P);
        }
        int rank = Projector(0, 0).Rank();

        for (int t = 0; t < 4; t++) {
            for (int s = 0; s < 4; s++) {
                for (int r = 0; r < 4; r++) {
                    for (int u = 0; u < 4; u++) {
                        temp_storage = takeInnerProduct((Projector)(u, t), (WaveFunc)(s, r), rank)
                                * threadloc_KMatrix[omp_get_thread_num()](t, s)(r, u);
                        result += temp_storage;
                    }
                }
            }
        }
        return result;
    }
};

#endif //_DSEPP_TWOQUARKKERNEL_HPP_
