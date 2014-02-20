/*
 * OneLoopIntegrator.hpp
 *
 *  Created on: Feb 12, 2014
 *      Author: stkubr
 */

#ifndef ONELOOPINTEGRATOR_HPP_
#define ONELOOPINTEGRATOR_HPP_

class C_OneLoopIntegrator{
	public:
	t_cmplxArray1D integrand_args;
	int k_col,Int_counter;
	C_Kinematics_1loop Momenta;
	C_Integrator_Line<t_cmplxMatrix,C_OneLoopIntegrator,double> * Integ_radial;
	C_Integrator_Line<t_cmplxMatrix,C_OneLoopIntegrator,double> * Integ_angle_cheb;
	C_Integrator_Line<t_cmplxMatrix,C_OneLoopIntegrator,double> * Integ_angle_Y;

	bool flag_sigma;

	C_OneLoopIntegrator(){
		integrand_args.resize(3);
	}

	//______________________________________________________________________
	// Multidimensional integration routine
	t_cmplxMatrix quad3d(){
		//integrand=func;
		k_col=0;
		Int_counter=0;
		return Integ_radial->getResult(&C_OneLoopIntegrator::f1,this);
	}
	t_cmplxMatrix f1(double k){
		k_col++;
		integrand_args[0]=(k);
		return Integ_angle_cheb->getResult(&C_OneLoopIntegrator::f2,this);
	}
	t_cmplxMatrix f2(double z){
		integrand_args[1]=z;
		flag_sigma=true;
		return Integ_angle_Y->getResult(&C_OneLoopIntegrator::f3,this);
	}
	t_cmplxMatrix f3(double y){
		integrand_args[2]=y;
		return (this->Integrand)(integrand_args);
	}

	virtual t_cmplxMatrix Integrand(t_cmplxArray1D args) {t_cmplxMatrix dummy; std::cout << "Error virtual call" << std::endl; assert(false); return dummy;}
};



#endif /* ONELOOPINTEGRATOR_HPP_ */
