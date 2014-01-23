/*
 * PathesUnitTest.hpp
 *
 *  Created on: Jan 23, 2014
 *      Author: stkubr
 */

#ifndef PATHESUNITTEST_HPP_
#define PATHESUNITTEST_HPP_

#include "../../source/NumLibs/Paths.hpp"


class C_PathsUnitTest{

public:
	// PathTest
	t_cmplxArray1D PathTest(int PathType);

	// SinglePointTest
	t_cmplx SinglePointTest(double t_paramtr, int PathType);

private:
	t_cmplxArray1D getPathOnVector(t_cmplxArray1D SamplePoints,
			Geometry::C_Path * PathInstance);

	t_cmplx getPathAtPoint(t_cmplx t_paramtr,
			Geometry::C_Path * PathInstance);

	t_cmplxArray1D SamplePointsForTest(int NumTestPoints);

};

#endif /* PATHESUNITTEST_HPP_ */
