/*
 * IntegrationUnitTest.cpp
 *
 *
 *  Created on: Jan 23, 2014
 *      Author: stkubr
 */

#include "../source/NumLibs/Integrator.h"
#include "../source/Mockups/IntegrationMockups.hpp"
#include "IntegrationUnitTest.hpp"

double C_IntegrationUnitTest::IntegrationUnitTest(Integrator_ID integrator_id, int selectFunctionType) {
	double dUnitTestResult;
	C_Integration_Mockup IntegrationUnitTester;
	IntegrationUnitTester.setFunctionType(selectFunctionType);
	C_Integrator_Line<t_dMatrix, C_Integration_Mockup, double> * IntegratorGaussLegendre;

	IntegratorGaussLegendre = C_Integrator_Line
			<t_dMatrix, C_Integration_Mockup,double>::createIntegrator(
			IntegrationUnitTester.getNumIntegrationPoints(),
			IntegrationUnitTester.getDownLimit(),
			IntegrationUnitTester.getUpLimit(), 1, integrator_id);

	dUnitTestResult = IntegratorGaussLegendre->getResult(
			&C_Integration_Mockup::IntegrandMockup, &IntegrationUnitTester)(0,
			0);
	std::cout.precision(16);

	delete IntegratorGaussLegendre;

	return dUnitTestResult;
}

