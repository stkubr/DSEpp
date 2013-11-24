#include <math.h>
#include <complex>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <getopt.h>
#include <time.h>
#include <sstream>
#include <vector>
#include <memory>
#include "omp.h"

#include "../eigen/Eigen/Eigenvalues"

#include "../qft++/include/matrix.h"
#include "../qft++/include/tensor.h"
#include "../qft++/include/relativistic-quantum-mechanics.h"

using namespace std;

typedef complex<double> t_cmplx;
typedef Vector4<complex<double> > t_cmplxVector;
typedef Matrix<complex<double> > t_cmplxMatrix;
typedef Matrix<Tensor<complex<double> > > t_cmplxDirac;
typedef Matrix<double> t_dMatrix;
typedef Tensor<complex<double> > t_cmplxTensor;
typedef Tensor<double> t_dTensor;
typedef vector < vector<t_cmplxTensor> > t_cmplx2DVecTen;
typedef Matrix< t_cmplxMatrix > t_cmplxMatrix2D;

typedef vector<double> t_dArray1D;
typedef vector<vector<double> > t_dArray2D;
typedef vector<vector<vector<double> > > t_dArray3D;

typedef vector<t_cmplx> t_cmplxArray1D;
typedef vector<vector<t_cmplx> > t_cmplxArray2D;
typedef vector<vector<vector<t_cmplx> > > t_cmplxArray3D;
typedef vector<t_cmplxArray3D> t_cmplxArray4D;

#define NUM_PRECISION 7
#define ORTH_CHECK 0
#define BASIS_TYPE 1

#define SHIFT 1.006

#include "Abs/AbstractClass.hpp"
#include "NumLibs/Linear_interpolation.h"
#include "NumLibs/Integrator.h"
#include "NumLibs/Support_functions.h"
#include "DedicMem/DedicMem.h"
#include "Abs/AbsDiagram.h"
#include "Abs/Kinematics.hpp"
#include "Abs/Loop_Integrator.hpp"
#include "Kernel/Gluon.hpp"
#include "Kernel/AbstractKernel.hpp"
#include "Kernel/Kernels.hpp"
#include "DSE/DSE.h"
#include "DSE/Quark.h"
#include "DSE/Quark_id.h"
#include "Vertex/BSE/BSE_Mesons.hpp"
#include "Abs/1LoopDiagram.hpp"
#include "Abs/2LoopDiagram.hpp"
#include "Abs/FormFactor.hpp"
#include "Abs/FF_SelfCoupling.hpp"
#include "Abs/FF_Seagull.hpp"
#include "Vertex/BSE/BSE_Builder.hpp"
#include "Abs/Manipulator.hpp"
#include "Abs/Invoker.hpp"

//#include "test_proj/test_EV.hpp"

/// main()
int main(int __argc,char *__argv[]){
	// Time initialization
	double tstart, tstop, ttime;
	tstart = (double)clock()/CLOCKS_PER_SEC;
	/*
	C_Pion_Invoker invoker;
	invoker.execute();
		*/
	tstop = (double)clock()/CLOCKS_PER_SEC;
ttime= (tstop-tstart)/omp_get_max_threads();
cout << "\n\n" << "calculation time=" << "  " << ttime << "  " << "seconds" << endl;
 return 0;
}
