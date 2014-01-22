#pragma once

class C_IntegrationUnitTest{
	public:
	C_IntegrationUnitTest(){
		
	}
	
	double IntegrationUnitTest(Integrator_ID integrator_id, int selectFunctionType=0){
		double dUnitTestResult;
		C_Integration_Mockup IntegrationUnitTester;
		IntegrationUnitTester.setFunctionType(selectFunctionType);
		C_Integrator_Line<t_dMatrix,C_Integration_Mockup,double> * IntegratorGaussLegendre;
		IntegratorGaussLegendre=C_Integrator_Line<t_dMatrix,C_Integration_Mockup,double>::createIntegrator(IntegrationUnitTester.getNumIntegrationPoints(),
																										  IntegrationUnitTester.getDownLimit(),
																							              IntegrationUnitTester.getUpLimit(), 1, integrator_id);
		dUnitTestResult = IntegratorGaussLegendre->getResult(&C_Integration_Mockup::IntegrandMockup, &IntegrationUnitTester)(0,0);
		std::cout.precision(16);
		//std::cout << dUnitTestResult << std::endl;
		//BOOST_CHECK_CLOSE(dRefResult,dUnitTestResult,Precision);
		delete IntegratorGaussLegendre;
		return dUnitTestResult;
	}
};

class C_PathesUnitTest{
	public:
	
	// PathTest
	t_cmplxArray1D PathTest(int PathType){
		t_cmplx TestValue;
		int NumSamplePoints = 3;
		t_cmplxArray1D SamplePoints(NumSamplePoints,0);
		t_cmplxArray1D ReturnPoints(NumSamplePoints,0);
		t_cmplx ParabolaApex = t_cmplx(0.0,1.0);
		t_cmplx AxisShiftY = t_cmplx(1.0,0.0);
		
		SamplePoints = SamplePointsForTest(NumSamplePoints);
		
		switch (PathType) {
			case 0:{
				Geometry::C_Parabola ParabolaInstance(ParabolaApex);
				ReturnPoints = getPathOnVector(SamplePoints, &ParabolaInstance);
			}break;
			case 1:{
				Geometry::C_Line LineInstance(ParabolaApex, AxisShiftY);
				ReturnPoints = getPathOnVector(SamplePoints, &LineInstance);
			}break;
			default:
				std::cout << "Error: PathTest got wrong PathType!" << std::endl;
				assert(false);
		}
		
		return ReturnPoints;
	}
	
	// SinglePointTest
	t_cmplx SinglePointTest(double t_paramtr, int PathType){
		t_cmplx ParabolaApex = t_cmplx(0.0,1.0);
		t_cmplx AxisShiftY = t_cmplx(1.0,0.0);
		t_cmplx TestValue=0.0;
		
		switch (PathType) {
			case 0:{
				Geometry::C_Parabola ParabolaInstance(ParabolaApex);
				TestValue = getPathAtPoint(t_paramtr, &ParabolaInstance);
			}break;
			case 1:{
				Geometry::C_Line LineInstance(ParabolaApex, AxisShiftY);
				TestValue = getPathAtPoint(t_paramtr, &LineInstance);
			}break;
			default:
				std::cout << "Error: PathTest got wrong PathType!" << std::endl;
				assert(false);
		}
		
		return TestValue;
	}
	
	t_cmplxArray1D getPathOnVector(t_cmplxArray1D SamplePoints, Geometry::C_Path * PathInstance){
		return PathInstance->getPathOnVector(SamplePoints);
	}

	t_cmplx getPathAtPoint(t_cmplx t_paramtr, Geometry::C_Path * PathInstance ){
		return PathInstance->getPathAt(t_paramtr);
	}
	
	t_cmplxArray1D SamplePointsForTest(int NumTestPoints){
		double DownLimit,UpLimit,x,dx;
		t_cmplxArray1D ReturnSamplePoints(NumTestPoints,0);
		
		DownLimit=0.0;
		UpLimit=1.0;
		dx=(UpLimit + DownLimit)/(NumTestPoints-1);
		x=DownLimit;
		
		for (int i = 0; i <= NumTestPoints-1; i++){
			ReturnSamplePoints[i]=x;
			x+=dx;
		}
		
		return ReturnSamplePoints;
	}
		
};

