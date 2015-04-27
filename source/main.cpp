#include <time.h>
#include <iomanip>
#include "omp.h"

#define NUM_PRECISION 7
#define ORTH_CHECK 0
#define BASIS_TYPE 1

#include "DSEpp.hpp"

/// main()
int main(int __argc, char *__argv[]) {
    // Time initialization
    double tstart, tstop, ttime;
    tstart = (double) clock() / CLOCKS_PER_SEC;

     C_Propagator * up_quark;
     C_Propagator * gluon;
     C_AbstractKernel * kernel;

     Quark_ID quark_type;
     Kernel_ID kernel_type;
     Gluon_ID gluon_type;

     quark_type=Test_ID;
     kernel_type=RL_ID;
     gluon_type=Test_Gluon_ID;

     C_Kernel_Factory * KernelFactory = new C_Kernel_Factory;
     C_Quark_Factory * QuarkFactory = new C_Quark_Factory;
     C_Gluon_Factory * GluonFactory = new C_Gluon_Factory;

     up_quark=QuarkFactory->Create((int)quark_type);
     kernel=KernelFactory->Create(kernel_type);
     gluon=GluonFactory->Create((int)gluon_type);

     std::vector<C_Propagator *> Props;
     Props.push_back(gluon);
     Props.push_back(up_quark);

     kernel->setPropagators(Props);

     up_quark->linkToKernel(kernel);
     //up_quark->dressPropagator();

     std::cout << std::setprecision(16) << up_quark->checkSum() << endl;
     //61.49288903070202

/*
    C_BoundState_parameters params;
    std::string ParamPath = "Parameters_files/Mesons/Scalar_symmetric.txt";
    params.setParams(ParamPath);
    params.Print();
*/
    std::vector<C_Propagator *> Partons;
    Partons.push_back(up_quark);
    Partons.push_back(up_quark);
    C_BSE_Pion pion;
    pion.linkToPartons(Partons);
    pion.linkToKernel(kernel);
    t_cmplxVector P;
    P.SetP4(0,0,0,t_cmplx(0,0.1));
    pion.dressBSE(P);

    tstop = (double) clock() / CLOCKS_PER_SEC;
    ttime = (tstop - tstart) / omp_get_max_threads();
    cout << "\n\n" << "calculation time=" << "  " << ttime << "  " << "seconds" << endl;
    return 0;
}