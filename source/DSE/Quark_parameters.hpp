/*
 * Quark_parameters.hpp
 *
 *  Created on: Jan 13, 2014
 *      Author: stkubr
 */

#ifndef QUARK_PARAMETERS_HPP_
#define QUARK_PARAMETERS_HPP_
#include <fstream>
#include <iostream>
#include <string>
#include <assert.h>

namespace Propagators {

/**
 * /brief Parameters for a quark
 *
 * Integration points, bare mass, etc..
*/
	class C_Quark_parameters {
	public:
		int num_prop_steps, num_cutoff_steps, num_angle_steps;
		double m0, mu, LimUk, LimDk, M2_contour, EffectiveCutoff, flag_LightOrHeavyQuark, Accuracy;
		bool flag_loadPropagator;

	public:
		void Print();

		void ReadParameters(std::string &_ParamPath);

		void setDefault();

	};

}
#endif /* QUARK_PARAMETERS_HPP_ */
