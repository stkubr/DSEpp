/*
 * Quark_parameters.hpp
 *
 *  Created on: Jan 13, 2014
 *      Author: stkubr
 */

#ifndef QUARK_PARAMETERS_HPP_
#define QUARK_PARAMETERS_HPP_
#include <fstream>

class C_Quark_parameters{
	public:
	int num_prop_steps,num_angle_steps;
	double m0,mu,LimUk,LimDk,M2_contour,EffectiveCutoff,HeavyLight,Accuracy;
	bool ReCalcProp;

	public:
	void Print();
	void ReadParameters(std::ifstream & _ParamList);

};

#endif /* QUARK_PARAMETERS_HPP_ */
