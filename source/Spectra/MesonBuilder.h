//
// Created by stkubr on 02.05.15.
//

#ifndef DSEPP_MESONBUILDER_H
#define DSEPP_MESONBUILDER_H


#include "AbstractBuilder.h"

class C_MesonBuilder: public C_AbstractBuilder {
public:
    Kernels::Kernel_ID kernel_ID;
    Propagators::Quark_ID quark_1_ID;
    Propagators::Quark_ID quark_2_ID;
    Propagators::Gluon_ID gluon_ID;

    std::vector<Propagators::C_Propagator*> Propagators;
    std::vector<Kernels::C_AbstractKernel*> Kernels;
    std::vector<BSE::C_BSE*>  BSEs;

    C_MesonBuilder(Propagators::Quark_ID Q1_id, Propagators::Quark_ID Q2_id, Kernels::Kernel_ID K1_id, Propagators::Gluon_ID G1_id){
        quark_1_ID=Q1_id;
        quark_2_ID=Q2_id;
        kernel_ID=K1_id;
        gluon_ID=G1_id;
        Propagators.resize(3);
        Kernels.resize(1);
        BSEs.resize((int)BSE::BSE_ID_End);
    }

    void SymmetricQuarksDetector(){
        Propagators::C_Quark_Factory QuarkFactory = Propagators::C_Quark_Factory::instance();
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
        Propagators::C_Gluon_Factory GluonFactory = Propagators::C_Gluon_Factory::instance();
        Propagators[0]=GluonFactory.Create(gluon_ID);
        SymmetricQuarksDetector();
        Spectra->setPropagators(Propagators);
    }

    void buildKernels(){
        Kernels::C_Kernel_Factory KernelFactory = Kernels::C_Kernel_Factory::instance();
        Kernels[0]=KernelFactory.Create(kernel_ID);
        Spectra->setKernels(Kernels);
    }

    void buildBSEs(){
        BSE::C_BSE_Factory BSEFactory = BSE::C_BSE_Factory::instance();
        for (int i = 0; i < BSE::BSE_ID_End; i++){
            BSEs[i]=BSEFactory.Create((BSE::BSE_ID)(i));
        }
        Spectra->setBSEs(BSEs);
    }

    void linkThemAll(){
        Spectra->linkAllPieces();
    }
};

extern C_MesonBuilder TestBuilder;//(Test_ID, Test_ID, RL_ID, Test_Gluon_ID);

#endif //DSEPP_MESONBUILDER_H
