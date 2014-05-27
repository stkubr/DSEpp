#define BOOST_TEST_MODULE Test
#include <boost/test/included/unit_test.hpp>

#include <complex>
#include <time.h>
#include "omp.h"

#define NUM_PRECISION 7
#define ORTH_CHECK 0
#define BASIS_TYPE 1

#include "../source/types.h"
#include "../source/DSE/Gluon.hpp"
#include "../source/Kernel/KernelFactory.hpp"
#include "../source/DSE/Quark_id.h"

#include "../source/Mockups/IntegrationMockups.hpp"
#include "IntegrationTests/IntegrationUnitTest.hpp"
#include "GeometryTests/PathsUnitTest.hpp"

//#include "../source/Vertex/BSE/BSE_Mesons.hpp"
//#include "../source/Vertex/BSE/BSE_Builder.hpp"
//----------------------------------------------------------------------
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

//----------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE(IntegrationTest_RainbowLadder_RL)
BOOST_AUTO_TEST_CASE(int_test)
{
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

	BOOST_CHECK_CLOSE(ref_value[0],value[0],1.e-4);
	BOOST_CHECK_CLOSE(ref_value[1],value[1],1.e-4);


	tstop = (double)clock()/CLOCKS_PER_SEC;
	ttime= (tstop-tstart)/omp_get_max_threads();
	std::cout << "\n\n" << "calculation time=" <<" "<< ttime <<" "<< "seconds" << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()

//----------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE(IntegrationTest_PseudoScalar_PS)
BOOST_AUTO_TEST_CASE(int_test)
{
	double tstart, tstop, ttime;
	tstart = (double)clock()/CLOCKS_PER_SEC;
		
	C_Propagator * up_quark;
	C_Propagator * gluon;
	C_AbstractKernel * kernel;
	
	Quark_ID quark_type;
	Kernel_ID kernel_type;
	Gluon_ID gluon_type;

	quark_type=Test_ID;
	kernel_type=RL_PS_ID;
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


	kernel->setExchangeID(Pion_exchange_ID);
	kernel->setConvolutionType(0);
	kernel->setPropagators(Props);
	kernel->setMesonExchangeMass(0.0);
	
	up_quark->LinkToKernel(kernel);
	up_quark->DressPropagator();

	
	t_dArray1D ref_value(2,0);
	ref_value[0]=129.66153264449591;
	ref_value[1]= 34.55210436280955;
	
	t_dArray1D value(2,0);
	value=up_quark->GetTotalSum();
	
	BOOST_CHECK_CLOSE(ref_value[0],value[0],1.e-3);
	BOOST_CHECK_CLOSE(ref_value[1],value[1],1.e-3);
	
		
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
