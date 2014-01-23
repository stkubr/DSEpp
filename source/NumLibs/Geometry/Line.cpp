/*
 * Line.cpp
 *
 *  Created on: Jan 23, 2014
 *      Author: stkubr
 */

#include "Line.hpp"

namespace Geometry {

C_Line::C_Line(t_cmplx __k_coeff, t_cmplx __b_coeff) {
	k_coeff = __k_coeff;
	b_coeff = __b_coeff;
}

// f(t)=t*k_coeff + b_coeff
t_cmplx C_Line::getPathAt(t_cmplx t_paramtr) {
	return t_paramtr * k_coeff + b_coeff;
}

// f^{\prime}(t)=k_coeff
t_cmplx C_Line::getDerivativePathAt(t_cmplx t_paramtr) {
	return k_coeff;
}

}

