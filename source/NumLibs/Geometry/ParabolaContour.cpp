/*
 * ParabolaContour.cpp
 *
 *  Created on: Jan 23, 2014
 *      Author: stkubr
 */

#include "ParabolaContour.hpp"

namespace Geometry {

C_ParabolaContour::C_ParabolaContour(t_cmplx __ParabolaApex, t_cmplx __k_coeff, t_cmplx __b_coeff){
	ContourPath.resize(2);
	parabola = std::make_shared<C_Parabola> (__ParabolaApex);
	line = std::make_shared<C_Line> (__k_coeff, __b_coeff);
}

void C_ParabolaContour::setParabolaContour(const t_dArray1D& p_parabola,
										   const t_dArray1D& w_parabola,
										   const t_dArray1D& p_line,
										   const t_dArray1D& w_line){
	/*
	 * Normally Cauchy contour integration goes counter-clockwise,
	 * so the weights of some pieces (going clock-wise) of contour has to have
	 * "-" sign to contribute to total contour integral with a right sign.
	 *
	 * Also the contour is symmetric in respect to real line - so
	 * the lower parts are just complex conjugation of upper parts.
	 */

	// Upper parabola of contour (clock-wise)
	int i=0;
	for_each(p_parabola.begin(), p_parabola.end(), [&](t_cmplx t_param) {
		ContourPath[0].push_back(parabola->getPathAt(t_param));
		ContourPath[1].push_back(parabola->getDerivativePathAt(t_param)*-1.0*w_parabola[i]);
		i++;
	});

	// Upper line of contour
	i=0;
	for_each(p_line.begin(), p_line.end(), [&](t_cmplx t_param) {
		ContourPath[0].push_back(line->getPathAt(t_param));
		ContourPath[1].push_back(line->getDerivativePathAt(t_param)*w_line[i]);
		i++;
	});

	// Lower line of contour (clock-wise)
	i=0;
	for_each(p_line.begin(), p_line.end(), [&](t_cmplx t_param) {
		ContourPath[0].push_back(conj(line->getPathAt(t_param)));
		ContourPath[1].push_back(conj(line->getDerivativePathAt(t_param))*-1.0*w_line[i]);
		i++;
	});

	// Lower parabola of contour
	i=0;
	for_each(p_parabola.begin(), p_parabola.end(), [&](t_cmplx t_param) {
		ContourPath[0].push_back(conj(parabola->getPathAt(t_param)));
		ContourPath[1].push_back(conj(parabola->getDerivativePathAt(t_param))*w_parabola[i]);
		i++;
	});

}

}



