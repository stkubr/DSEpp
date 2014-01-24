#define BOOST_TEST_MODULE Test
#include <boost/test/included/unit_test.hpp>

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


//using namespace std;

typedef complex<double> t_cmplx;
typedef Vector4<complex<double> > t_cmplxVector;
typedef Matrix<complex<double> > t_cmplxMatrix;
typedef Matrix<Tensor<complex<double> > > t_cmplxDirac;
typedef Matrix<double> t_dMatrix;
typedef Tensor<complex<double> > t_cmplxTensor;
typedef Tensor<double> t_dTensor;
typedef std::vector < std::vector<t_cmplxTensor> > t_cmplx2DVecTen;
typedef Matrix< t_cmplxMatrix > t_cmplxMatrix2D;

typedef std::vector<double> t_dArray1D;
typedef std::vector<std::vector<double> > t_dArray2D;
typedef std::vector<std::vector<std::vector<double> > > t_dArray3D;

typedef std::vector<t_cmplx> t_cmplxArray1D;
typedef std::vector<std::vector<t_cmplx> > t_cmplxArray2D;
typedef std::vector<std::vector<std::vector<t_cmplx> > > t_cmplxArray3D;
typedef std::vector<t_cmplxArray3D> t_cmplxArray4D;

#define NUM_PRECISION 7
#define ORTH_CHECK 0
#define BASIS_TYPE 1
#define DEBUGLEVEL 0



/*
#include "../source/NumLibs/Linear_interpolation.h"
#include "../source/NumLibs/Integrator.h"
#include "../source/NumLibs/Support_functions.h"
#include "../source/DedicMem/DedicMem.h"
#include "../source/Abs/AbsDiagram.h"
#include "../source/Kernel/Gluon.hpp"
#include "../source/Kernel/AbstractKernel.hpp"*/
#include "../source/Kernel/Kernels.hpp"
#include "../source/DSE/Quark_id.h"
#include "../source/Mockups/IntegrationMockups.hpp"
#include "IntegrationUnitTest.hpp"
#include "GeometryTests/PathsUnitTest.hpp"



BOOST_AUTO_TEST_SUITE(GeometryTest)
BOOST_AUTO_TEST_CASE(SinglePointTest)
{
	C_PathsUnitTest ParabolaTest;
	t_cmplx ReturnValue=ParabolaTest.SinglePointTest(0.0,0);
	
	BOOST_CHECK_CLOSE(real(ReturnValue),-1.0,1.e-8);
	BOOST_CHECK_EQUAL(imag(ReturnValue),0.0);
	
	C_PathsUnitTest LineTest;
	ReturnValue=LineTest.SinglePointTest(0.0,1);

	BOOST_CHECK_CLOSE(real(ReturnValue),1.0,1.e-8);
	BOOST_CHECK_EQUAL(imag(ReturnValue),0.0);

}

BOOST_AUTO_TEST_CASE(PathTest)
{
	C_PathsUnitTest ParabolaTest;
	t_cmplxArray1D ReturnPoints;
	ReturnPoints = ParabolaTest.PathTest(0);
	
	BOOST_CHECK_CLOSE(real(ReturnPoints[0]),-1.0,1.e-8);
	BOOST_CHECK_EQUAL(imag(ReturnPoints[0]),0.0);
	
	BOOST_CHECK_CLOSE(real(ReturnPoints[1]),-0.75,1.e-8);
	BOOST_CHECK_CLOSE(imag(ReturnPoints[1]),1.0,1.e-8);
	
	BOOST_CHECK_CLOSE(real(ReturnPoints[2]),0,1.e-8);
	BOOST_CHECK_CLOSE(imag(ReturnPoints[2]),2,1.e-8);


	C_PathsUnitTest LineTest;
	ReturnPoints = LineTest.PathTest(1);

	BOOST_CHECK_CLOSE(real(ReturnPoints[0]),1.0,1.e-8);
	BOOST_CHECK_EQUAL(imag(ReturnPoints[0]),0.0);

	BOOST_CHECK_CLOSE(real(ReturnPoints[1]),1,1.e-8);
	BOOST_CHECK_CLOSE(imag(ReturnPoints[1]),0.5,1.e-8);

	BOOST_CHECK_CLOSE(real(ReturnPoints[2]),1,1.e-8);
	BOOST_CHECK_CLOSE(imag(ReturnPoints[2]),1,1.e-8);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(IntegrationTest)
BOOST_AUTO_TEST_CASE(int_test)
{
	double tstart, tstop, ttime;
	tstart = (double)clock()/CLOCKS_PER_SEC;
		
	C_Propagator * up_quark;
	C_AbstractKernel * kernel;
	
	Quark_ID quark_type;
	Kernel_ID kernel_type;

	quark_type=Up_ID;
	kernel_type=RL_ID;
	
	C_Kernel_Factory * KernelFactory = new C_Kernel_Factory;
	C_Quark_Factory * QuarkFactory = new C_Quark_Factory;
	
	up_quark=QuarkFactory->Create(&quark_type);
	kernel=KernelFactory->Create(&kernel_type);
	kernel->SpecifyGluon(RL_MT_Light_ID);
	up_quark->LinkToKernel(kernel);
	
	up_quark->DressPropagator();

	
	t_dArray1D ref_value(2,0);
	ref_value[0]=134.095477101;
	ref_value[1]=38.513658633;
	
	t_dArray1D value(2,0);
	value=up_quark->GetTotalSum();
	
	BOOST_CHECK_CLOSE(ref_value[0],value[0],1.e-4);
	BOOST_CHECK_CLOSE(ref_value[1],value[1],1.e-4);
	
		
	tstop = (double)clock()/CLOCKS_PER_SEC;
	ttime= (tstop-tstart)/omp_get_max_threads();
	std::cout << "\n\n" << "calculation time=" <<" "<< ttime <<" "<< "seconds" << std::endl;

}

BOOST_AUTO_TEST_SUITE_END()

//----------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE(GauLeg)
BOOST_AUTO_TEST_CASE(GauLegTest1)
{
	C_IntegrationUnitTest GauLegTest;
	double UnitTestResult=GauLegTest.IntegrationUnitTest(qgausleg_sym_ID,0);
	BOOST_CHECK_SMALL(UnitTestResult,1.e-10);
}

BOOST_AUTO_TEST_CASE(GauLegTest2)
{
	C_IntegrationUnitTest GauLegTest;
	double UnitTestResult=GauLegTest.IntegrationUnitTest(qgausleg_sym_ID,1);
	BOOST_CHECK_CLOSE(1.570796326794897,UnitTestResult,1.e-6);
}

BOOST_AUTO_TEST_CASE(GauLegTest3)
{
	C_IntegrationUnitTest GauLegTest;
	double UnitTestResult=GauLegTest.IntegrationUnitTest(qgausleg_sym_ID,2);
	BOOST_CHECK_CLOSE(2.434347529565719,UnitTestResult,1.e-6);
}

BOOST_AUTO_TEST_SUITE_END()

//----------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE(GauCheb)
BOOST_AUTO_TEST_CASE(GauChebTest1)
{
	C_IntegrationUnitTest GauLegTest;
	double UnitTestResult=GauLegTest.IntegrationUnitTest(qgauscheb_ID,0);
	BOOST_CHECK_SMALL(UnitTestResult,1.e-10);
}

BOOST_AUTO_TEST_CASE(GauChebTest2)
{
	C_IntegrationUnitTest GauLegTest;
	double UnitTestResult=GauLegTest.IntegrationUnitTest(qgauscheb_ID,3);
	BOOST_CHECK_CLOSE(1.570796326794897,UnitTestResult,1.e-6);
}

BOOST_AUTO_TEST_CASE(GauChebTest3)
{
	C_IntegrationUnitTest GauLegTest;
	double UnitTestResult=GauLegTest.IntegrationUnitTest(qgauscheb_ID,4);
	BOOST_CHECK_CLOSE(2.434347529565719,UnitTestResult,1.e-3);
}
BOOST_AUTO_TEST_SUITE_END()

//----------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE(GauLegLin)
BOOST_AUTO_TEST_CASE(GauLegLinTest1)
{
	C_IntegrationUnitTest GauLegTest;
	double UnitTestResult=GauLegTest.IntegrationUnitTest(qgausleg_lin_ID,0);
	BOOST_CHECK_SMALL(UnitTestResult,1.e-6);
}

BOOST_AUTO_TEST_CASE(GauLegLinTest2)
{
	C_IntegrationUnitTest GauLegTest;
	double UnitTestResult=GauLegTest.IntegrationUnitTest(qgausleg_lin_ID,1);
	BOOST_CHECK_CLOSE(1.570796326794897,UnitTestResult,1.e-6);
}

BOOST_AUTO_TEST_CASE(GauLegLinTest3)
{
	C_IntegrationUnitTest GauLegTest;
	double UnitTestResult=GauLegTest.IntegrationUnitTest(qgausleg_lin_ID,2);
	BOOST_CHECK_CLOSE(2.434347529565719,UnitTestResult,1.e-6);
}

BOOST_AUTO_TEST_SUITE_END()
