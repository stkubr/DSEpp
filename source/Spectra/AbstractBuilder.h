//
// Created by stkubr on 02.05.15.
//

#ifndef DSEPP_ABSTRACTBUILDER_H
#define DSEPP_ABSTRACTBUILDER_H

#include "Spectra.h"

class C_AbstractBuilder {
public:
    C_Spectra * Spectra;

    virtual ~C_AbstractBuilder() {}

    C_Spectra * GetPhysicalState() {return Spectra;}

    void createNewPhysicalState() {Spectra = new C_Spectra();}

    virtual void buildPropagators()=0;
    virtual void buildKernels()=0;
    virtual void buildBSEs()=0;
    virtual void linkThemAll()=0;
};

#endif //DSEPP_ABSTRACTBUILDER_H
