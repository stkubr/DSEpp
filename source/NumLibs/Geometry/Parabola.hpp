/*
 * Parabola.hpp
 *
 *  Created on: Jan 23, 2014
 *      Author: stkubr
 */

#ifndef PARABOLA_HPP_
#define PARABOLA_HPP_

#include "../../types.h"
#include "Path.hpp"

namespace Geometry {

class C_Parabola: public C_Path{
private:
	t_cmplx ParabolaApex;

public:
	C_Parabola(t_cmplx __ParabolaApex);

	// f(t)=t*t + 2.0*t*Apex + Apex*Apex
	t_cmplx getPathAt(t_cmplx t_paramtr);

	// f^{\prime}(t)=2*t + 2.0*Apex
	t_cmplx getDerivativePathAt(t_cmplx t_paramtr);
};

}

#endif /* PARABOLA_HPP_ */
