//
// Created by stkubr on 02.05.15.
//

#ifndef DSEPP_SPECTRA_H
#define DSEPP_SPECTRA_H

#include "../DSEpp.hpp"

class C_Spectra {
protected:
    std::vector<Propagators::C_Propagator*> Propagators;
    std::vector<Kernels::C_AbstractKernel*> Kernels;
    std::vector<BSE::C_BSE*>  BSEs;

public:
    ~C_Spectra(){
        deleteAllPieces();
    }

    void setPropagators(std::vector<Propagators::C_Propagator *> _Propagators){
        Propagators=_Propagators;
    }

    void setKernels(std::vector<Kernels::C_AbstractKernel *> _Kernels){
        Kernels=_Kernels;
    }

    void setBSEs(std::vector<BSE::C_BSE *> _BSEs){
        BSEs=_BSEs;
    }

    Propagators::C_Propagator * getPropagators(int i) const {
        if(Propagators[i]) return Propagators[i];
    }

    Kernels::C_AbstractKernel * getKernels(int i) const {
        if(Kernels[i]) return Kernels[i];
    }

    BSE::C_BSE * getBSEs(int i) const {
        if(BSEs[i]) return BSEs[i];
    }

    void linkAllPieces();

    void checkAllPieces();

    void deleteAllPieces();
};


#endif //DSEPP_SPECTRA_H
