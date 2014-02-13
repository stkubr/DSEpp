/*
 * IntegrationUnitTest.hpp
 *
 *  Created on: Jan 23, 2014
 *      Author: stkubr
 */

#ifndef INTEGRATIONUNITTEST_HPP_
#define INTEGRATIONUNITTEST_HPP_

#include "../../source/NumLibs/Integrator.h"
#include "../../source/Mockups/IntegrationMockups.hpp"

class C_IntegrationUnitTest{
	public:
	double IntegrationUnitTest(Integrator_ID integrator_id, int selectFunctionType=0);
};

#endif /* INTEGRATIONUNITTEST_HPP_ */
