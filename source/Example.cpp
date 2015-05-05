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

  // create pointer to container for DSE objects (like: quarks and gluons) and BSE objects (like: mesons)
  C_Spectra * spectra;

  // set test builder
  // - will create two test quarks and Maris-Tandy gluon kernel
  // - link them
  // - will create a tower of BSE(mesons) for all possible quantum numbers (listed in enum BSE_ID)
  // - link them to aforementioned kernel and quarks
  // - after that is fully operational
  Binder.SetBSEBuilder(&TestBuilder);
  Binder.ConstructPhysState();
  spectra=Binder.GetPhysicalState();
  spectra->checkAllPieces();

  // create a 4-vector momenta with pion mass in the rest frame
  t_cmplxVector P;
  P.SetP4(0, 0, 0, t_cmplx(0, 0.138));

  // calculate eigenvalues using matrix for the pseudoscalar mesons
  spectra->getBSEs(0)->calcEigenStates(P, 10);

  P.SetP4(0, 0, 0, t_cmplx(0, 0.650));
  // calculate eigenvalues using matrix for the scalar mesons
  spectra->getBSEs(1)->calcEigenStates(P, 10);

  P.SetP4(0, 0, 0, t_cmplx(0, 0.757));
  // calculate eigenvalues using matrix for the vector mesons
  spectra->getBSEs(2)->calcEigenStates(P, 10);



  tstop = (double) clock() / CLOCKS_PER_SEC;
  ttime = (tstop - tstart) / omp_get_max_threads();
  cout << "\n\n" << "calculation time=" << "  " << ttime << "  " << "seconds" << endl;
  return 0;
}