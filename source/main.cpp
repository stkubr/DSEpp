#include <time.h>
#include <iomanip>

#include "omp.h"

#define NUM_PRECISION 7
#define ORTH_CHECK 0
#define BASIS_TYPE 1

#include "DSEpp.hpp"
#include "Spectra/Binder.h"
#include "Spectra/MesonBuilder.h"

/// main()
int main(int __argc, char *__argv[]) {
  // Time initialization
  double tstart, tstop, ttime;
  tstart = (double) clock() / CLOCKS_PER_SEC;

  /*  Quark_ID quark_type = Test_ID;
    Gluon_ID gluon_type = Test_Gluon_ID;
    Kernel_ID kernel_type = RL_ID;

    C_Kernel_Factory KernelFactory = C_Kernel_Factory::instance();
    C_Quark_Factory QuarkFactory = C_Quark_Factory::instance();
    C_Gluon_Factory GluonFactory = C_Gluon_Factory::instance();

    C_Propagator *up_quark = QuarkFactory.Create((int) quark_type);
    C_Propagator *gluon = GluonFactory.Create((int) gluon_type);
    C_AbstractKernel *kernel = KernelFactory.Create(kernel_type);

    std::vector<C_Propagator *> Props;
    Props.push_back(gluon);
    Props.push_back(up_quark);

    kernel->setPropagators(Props);

    up_quark->linkToKernel(kernel);
    //up_quark->dressPropagator();

    std::cout << std::setprecision(16) << up_quark->checkSum() << endl;

    std::vector<C_Propagator *> partons;
    partons.push_back(up_quark);
    partons.push_back(up_quark);

    C_BSE_Factory BSE_Factory = C_BSE_Factory::instance();

    C_BSE *pion = BSE_Factory.Create(PseudoScalar_ID);
    pion->linkToPartons(partons);
    pion->linkToKernel(kernel);

    t_cmplxVector P;
    P.SetP4(0, 0, 0, t_cmplx(0, 0.138));
    //pion.dressBSE(P);
    pion->calcEigenStates(P, 10);*/

  C_Spectra * spectra;

  Binder.SetBSEBuilder(&TestBuilder);
  Binder.ConstructPhysState();
  spectra=Binder.GetPhysicalState();
  spectra->checkAllPieces();

  t_cmplxVector P;
  P.SetP4(0, 0, 0, t_cmplx(0, 0.138));
  spectra->BSEs[0]->calcEigenStates(P, 10);

  P.SetP4(0, 0, 0, t_cmplx(0, 0.650));
  spectra->BSEs[1]->calcEigenStates(P, 10);

  P.SetP4(0, 0, 0, t_cmplx(0, 0.757));
  spectra->BSEs[2]->calcEigenStates(P, 10);



  tstop = (double) clock() / CLOCKS_PER_SEC;
  ttime = (tstop - tstart) / omp_get_max_threads();
  cout << "\n\n" << "calculation time=" << "  " << ttime << "  " << "seconds" << endl;
  return 0;
}