/*
 * Path.hpp
 *
 *  Created on: Jan 23, 2014
 *      Author: stkubr
 */

#ifndef PATH_HPP_
#define PATH_HPP_

#include "../../types.h"

namespace Geometry {

class C_Path {
public:
	t_cmplxArray1D getPathOnVector(t_cmplxArray1D SamplePoints);

	virtual t_cmplx getPathAt(t_cmplx t_paramtr);

	virtual t_cmplx getDerivativePathAt(t_cmplx t_paramtr);
};

}

#endif /* PATH_HPP_ */
