/*
 * ParabolaContour.hpp
 *
 *  Created on: Jan 23, 2014
 *      Author: stkubr
 */

#ifndef PARABOLACONTOUR_HPP_
#define PARABOLACONTOUR_HPP_

#include <vector>
#include <memory>
#include <algorithm>

#include "Path.hpp"
#include "Parabola.hpp"
#include "Line.hpp"


namespace Geometry {

class C_ParabolaContour {
private:
	std::shared_ptr<C_Path>  parabola;
	std::shared_ptr<C_Path>  line;

	t_cmplxArray2D ContourPath;

public:
	C_ParabolaContour(t_cmplx __ParabolaApex,
					  t_cmplx __k_coeff, t_cmplx __b_coeff);

	void setParabolaContour(const t_dArray1D& p_parabola,
							const t_dArray1D& w_parabola,
							const t_dArray1D& p_line,
							const t_dArray1D& w_line);

	t_cmplxArray2D getParabolaContour();

};

}

#endif /* PARABOLACONTOUR_HPP_ */
