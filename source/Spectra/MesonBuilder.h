//
// Created by stkubr on 02.05.15.
//

#ifndef DSEPP_MESONBUILDER_H
#define DSEPP_MESONBUILDER_H


#include "AbstractBuilder.h"

class C_MesonBuilder: public C_AbstractBuilder {
public:
    Kernel_ID kernel_ID;
    Quark_ID quark_1_ID;
    Quark_ID quark_2_ID;
    Gluon_ID gluon_ID;

    std::vector<C_Propagator*> Propagators;
    std::vector<C_AbstractKernel*> Kernels;
    std::vector<C_BSE*>  BSEs;

    C_MesonBuilder(Quark_ID Q1_id, Quark_ID Q2_id, Kernel_ID K1_id, Gluon_ID G1_id){
        quark_1_ID=Q1_id;
        quark_2_ID=Q2_id;
        kernel_ID=K1_id;
        gluon_ID=G1_id;
        Propagators.resize(3);
        Kernels.resize(1);
        BSEs.resize((int)BSE_ID_End);
    }

    void SymmetricQuarksDetector(){
        C_Quark_Factory QuarkFactory = C_Quark_Factory::instance();
        if(quark_1_ID==quark_2_ID){
            Propagators[1]=QuarkFactory.Create(quark_1_ID);
            Propagators[2]=Propagators[1];
        }
        else {
            Propagators[1]=QuarkFactory.Create(quark_1_ID);
            Propagators[2]=QuarkFactory.Create(quark_2_ID);
        }
    }

    void buildPropagators(){
        C_Gluon_Factory GluonFactory = C_Gluon_Factory::instance();
        Propagators[0]=GluonFactory.Create(gluon_ID);
        SymmetricQuarksDetector();
        Spectra->setPropagators(Propagators);
    }

    void buildKernels(){
        C_Kernel_Factory KernelFactory = C_Kernel_Factory::instance();
        Kernels[0]=KernelFactory.Create(kernel_ID);
        Spectra->setKernels(Kernels);
    }

    void buildBSEs(){
        C_BSE_Factory BSEFactory = C_BSE_Factory::instance();
        for (int i = 0; i < BSE_ID_End; i++){
            BSEs[i]=BSEFactory.Create((BSE_ID)(i));
        }
        Spectra->setBSEs(BSEs);
    }

    void linkThemAll(){
        Spectra->linkAllPieces();
    }
};

extern C_MesonBuilder TestBuilder;//(Test_ID, Test_ID, RL_ID, Test_Gluon_ID);

#endif //DSEPP_MESONBUILDER_H
