/*
 * Parabola.cpp
 *
 *  Created on: Jan 23, 2014
 *      Author: stkubr
 */

#include "Parabola.hpp"

namespace Geometry {

C_Parabola::C_Parabola(t_cmplx __ParabolaApex) {
	ParabolaApex = __ParabolaApex;
}

// f(t)=t*t + 2.0*t*Apex + Apex*Apex
t_cmplx C_Parabola::getPathAt(t_cmplx t_paramtr) {
	return t_paramtr * t_paramtr
			+ 2.0 * t_paramtr * ParabolaApex
			+ ParabolaApex * ParabolaApex;
}

// f^{\prime}(t)=2*t + 2.0*Apex
t_cmplx C_Parabola::getDerivativePathAt(t_cmplx t_paramtr) {
	return 2.0 * t_paramtr + 2.0 * ParabolaApex;
}

}
