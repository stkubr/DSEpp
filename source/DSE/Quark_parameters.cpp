/*
 * Quark_parameters.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: stkubr
 */

#include "Quark_parameters.hpp"

#include <iostream>

void C_Quark_parameters::Print() {
	 std::cout << "num_prop_steps" << " - " << num_prop_steps << std::endl <<
	 "num_angle_steps" << " - " << num_angle_steps <<std::endl <<
	 "m0" << " - " << m0 <<std::endl <<
	 "mu" << " - " << mu <<std::endl <<
	 "M2_contour" << " - " << M2_contour <<std::endl <<
	 "LimDk" << " - " << LimDk <<std::endl <<
	 "LimUk" << " - " << LimUk <<std::endl <<
	 "EffectiveCutoff" << " - " << EffectiveCutoff <<std::endl <<
	 "Accuracy" << " - " << Accuracy << std::endl <<
	 "ReCalcProp" << " - " << ReCalcProp << std::endl <<
	 "HeavyLight" << " - " << HeavyLight << std::endl;
}
