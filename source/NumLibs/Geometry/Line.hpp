/*
 * Line.hpp
 *
 *  Created on: Jan 23, 2014
 *      Author: stkubr
 */

#ifndef LINE_HPP_
#define LINE_HPP_

#include "Path.hpp"

namespace Geometry {

class C_Line: public C_Path {
private:
	t_cmplx k_coeff, b_coeff;

public:
	C_Line(t_cmplx __k_coeff, t_cmplx __b_coeff);

	// f(t)=t*k_coeff + b_coeff
	t_cmplx getPathAt(t_cmplx t_paramtr);

	// f^{\prime}(t)=k_coeff
	t_cmplx getDerivativePathAt(t_cmplx t_paramtr);
};

}

#endif /* LINE_HPP_ */
