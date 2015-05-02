//
// Created by stkubr on 02.05.15.
//

#ifndef DSEPP_SPECTRA_H
#define DSEPP_SPECTRA_H

#include "../DSEpp.hpp"

class C_Spectra {
public:
    std::vector<C_Propagator*> Propagators;
    std::vector<C_AbstractKernel*> Kernels;
    std::vector<C_BSE*>  BSEs;

public:
    ~C_Spectra(){
        deleteAllPieces();
    }

    void setPropagators(std::vector<C_Propagator *> _Propagators){
        Propagators=_Propagators;
    }

    void setKernels(std::vector<C_AbstractKernel *> _Kernels){
        Kernels=_Kernels;
    }

    void setBSEs(std::vector<C_BSE *> _BSEs){
        BSEs=_BSEs;
    }

    void linkAllPieces(){
        Kernels[0]->setPropagators({Propagators[0],Propagators[1]});
        Propagators[1]->linkToKernel(Kernels[0]);
        Propagators[2]->linkToKernel(Kernels[0]);
        for (int i = 0; i < BSEs.size(); i++){
            BSEs[i]->linkToKernel(Kernels[0]);
            BSEs[i]->linkToPartons({Propagators[1],Propagators[2]});
        }
    }

    void checkAllPieces(){
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

    void deleteAllPieces(){
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
};


#endif //DSEPP_SPECTRA_H
