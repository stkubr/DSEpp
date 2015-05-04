/*
 * Quark_parameters.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: stkubr
 */

#include "Quark_parameters.hpp"

using namespace Propagators;

void C_Quark_parameters::Print() {
	std::cout << "Points on Parabola" << " - " << num_prop_steps << std::endl <<
			"Points on Cutoff-line" << " - " << num_cutoff_steps << std::endl <<
			"num_angle_steps" << " - " << num_angle_steps <<std::endl <<
			"m0" << " - " << m0 <<std::endl <<
			"mu" << " - " << mu <<std::endl <<
			"M2_contour" << " - " << M2_contour <<std::endl <<
			"LimDk" << " - " << LimDk <<std::endl <<
			"LimUk" << " - " << LimUk <<std::endl <<
			"EffectiveCutoff" << " - " << EffectiveCutoff <<std::endl <<
			"Accuracy" << " - " << Accuracy << std::endl <<
			"flag_loadPropagator" << " - " << flag_loadPropagator << std::endl <<
			"flag_LightOrHeavyQuark" << " - " << flag_LightOrHeavyQuark << std::endl;
}

void C_Quark_parameters::ReadParameters(std::string& _ParamPath) {
	std::string line;
	std::ifstream _ParamList(_ParamPath);
	if ((_ParamList).is_open()) {
		while ((_ParamList).good()) {
			(_ParamList) >> line >> num_prop_steps;
			(_ParamList) >> line >> num_cutoff_steps;
			(_ParamList) >> line >> num_angle_steps;
			(_ParamList) >> line >> m0;
			(_ParamList) >> line >> mu;
			(_ParamList) >> line >> M2_contour;
			(_ParamList) >> line >> LimDk;
			(_ParamList) >> line >> LimUk;
			(_ParamList) >> line >> EffectiveCutoff;
			(_ParamList) >> line >> Accuracy;
			(_ParamList) >> line >> flag_loadPropagator;
			(_ParamList) >> line >> flag_LightOrHeavyQuark;
		}
		Print();
	} else {
		std::cout << "Cant open Quark parameters file!" << std::endl;
		assert(false);
	}
}

void C_Quark_parameters::setDefault(){
	std::cout << "Default Quark parameters!" << std::endl;
	num_prop_steps = 64;
	num_cutoff_steps = 10;
	num_angle_steps = 8;
	m0 = 0.0037;
	mu = 19.0;
	M2_contour = 1.0;
	LimDk = 0.0001;
	LimUk = 50.0;
	EffectiveCutoff = 0.8;
	Accuracy = 0.0001;
	flag_loadPropagator = 1;
	flag_LightOrHeavyQuark = 0;
}
