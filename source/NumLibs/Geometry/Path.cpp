/*
 * Path.cpp
 *
 *  Created on: Jan 23, 2014
 *      Author: stkubr
 */

#include "Path.hpp"

namespace Geometry {

t_cmplxArray1D C_Path::getPathOnVector(t_cmplxArray1D SamplePoints) {
	int NumSamplePoints;
	NumSamplePoints = SamplePoints.size();
	t_cmplxArray1D PathOnSamlePoints(NumSamplePoints, 0);

	for (int i = 0; i < NumSamplePoints; i++) {
		PathOnSamlePoints[i] = getPathAt(SamplePoints[i]);
	}

	return PathOnSamlePoints;
}

t_cmplx C_Path::getPathAt(t_cmplx t_paramtr) {
	assert(false);
	return 0.0;
}
t_cmplx C_Path::getDerivativePathAt(t_cmplx t_paramtr) {
	assert(false);
	return 0.0;
}

}


