#include <time.h>
#include "omp.h"

#define NUM_PRECISION 7
#define ORTH_CHECK 0
#define BASIS_TYPE 1

#include "../source/types.h"
#include "../source/DSE/Gluon.hpp"
#include "../source/Kernel/KernelFactory.hpp"
#include "../source/DSE/Quark_id.hpp"


/// main()
int main(int __argc,char *__argv[]){
    // Time initialization
    double tstart, tstop, ttime;
    tstart = (double)clock()/CLOCKS_PER_SEC;


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

    up_quark->LinkToKernel(kernel);
    up_quark->DressPropagator();

    t_dArray1D ref_value(2,0);
    ref_value[0]=134.092697854;
    ref_value[1]=38.5110544794;

    t_dArray1D value(2,0);
    value=up_quark->GetTotalSum();
//134.092691503 38.5110260692
    tstop = (double)clock()/CLOCKS_PER_SEC;
    ttime= (tstop-tstart)/omp_get_max_threads();
    cout << "\n\n" << "calculation time=" << "  " << ttime << "  " << "seconds" << endl;
    return 0;
}