//
// Created by stkubr on 02.05.15.
//

#ifndef DSEPP_BINDER_H
#define DSEPP_BINDER_H

#include "AbstractBuilder.h"

class C_Binder {
private:
    C_AbstractBuilder * BSE_Builder;
public:
    void SetBSEBuilder(C_AbstractBuilder * b) {BSE_Builder=b;}

    C_Spectra * GetPhysicalState() {return BSE_Builder->GetPhysicalState();}

    void ConstructPhysState(){
        BSE_Builder->createNewPhysicalState();
        BSE_Builder->buildPropagators();
        BSE_Builder->buildKernels();
        BSE_Builder->buildBSEs();
        BSE_Builder->linkThemAll();
    }
};

extern C_Binder Binder;

#endif //DSEPP_BINDER_H
