//
// Created by stkubr on 21.04.15.
//

#ifndef DSEPP_BSE_H
#define DSEPP_BSE_H


#include <source/types.h>
#include <source/Abs/AbsDiagram.hpp>
#include <source/Kernel/AbstractKernel.hpp>

class C_BSE: public C_AbsDiagram {
public:
    //virtual t_cmplxArray2D calcEigenvalues(t_cmplxVector P)=0;

    virtual void dressBSE(t_cmplxVector P)=0;

   // virtual void setBSEonPath(t_cmplxArray2D &AmplitudesOnPath, t_cmplxArray1D &Path)=0;

    //virtual double checkSum()=0;

    virtual void linkToKernel(C_AbstractKernel * _Kernel)=0;

    virtual void linkToPartons(std::vector<C_Propagator*> _Partons)=0;

    virtual ~C_BSE(){}
};


#endif //DSEPP_BSE_H
