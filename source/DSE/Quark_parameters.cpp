/*
 * Quark_parameters.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: stkubr
 */

#include "Quark_parameters.hpp"
#include <fstream>
#include <iostream>
#include <string.h>

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

void C_Quark_parameters::ReadParameters(std::ifstream & _ParamList) {
	std::string line;
	if ((_ParamList).is_open()) {
		while ((_ParamList).good()) {
			(_ParamList) >> line >> num_prop_steps;
			(_ParamList) >> line >> num_angle_steps;
			(_ParamList) >> line >> m0;
			(_ParamList) >> line >> mu;
			(_ParamList) >> line >> M2_contour;
			(_ParamList) >> line >> LimDk;
			(_ParamList) >> line >> LimUk;
			(_ParamList) >> line >> EffectiveCutoff;
			(_ParamList) >> line >> Accuracy;
			(_ParamList) >> line >> ReCalcProp;
			(_ParamList) >> line >> HeavyLight;
		}
		Print();
	} else
		std::cout << "Cant open Quark parameters file!" << std::endl;
}
