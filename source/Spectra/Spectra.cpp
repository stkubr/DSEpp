//
// Created by stkubr on 02.05.15.
//

#include "Spectra.h"


void C_Spectra::linkAllPieces() {
    Kernels[0]->setPropagators({Propagators[0],Propagators[1]});
    Propagators[1]->linkToKernel(Kernels[0]);
    Propagators[2]->linkToKernel(Kernels[0]);
    for (int i = 0; i < BSEs.size(); i++){
        BSEs[i]->linkToKernel(Kernels[0]);
        BSEs[i]->linkToPartons({Propagators[1],Propagators[2]});
    }
}

void C_Spectra::checkAllPieces() {
    for (int i = 0; i < Propagators.size(); i++){
        Propagators[i]->GetNameID();
    }
    for (int i = 0; i < Kernels.size(); i++){
        Kernels[i]->GetNameID();
    }
    for (int i = 0; i < BSEs.size()-1; i++){
        BSEs[i]->GetNameID();
    }
}

void C_Spectra::deleteAllPieces() {
    for (int i = 0; i < Propagators.size(); i++){
        if(Propagators[i]) delete (Propagators[i]);
    }
    for (int i = 0; i < Kernels.size(); i++){
        if(Kernels[i]) delete (Kernels[i]);
    }
    for (int i = 0; i < BSEs.size(); i++){
        if(BSEs[i]) delete (BSEs[i]);
    }
}