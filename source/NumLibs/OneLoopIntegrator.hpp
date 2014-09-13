/*
 * OneLoopIntegrator.hpp
 *
 *  Created on: Feb 12, 2014
 *      Author: stkubr
 */

#ifndef ONELOOPINTEGRATOR_HPP_
#define ONELOOPINTEGRATOR_HPP_

#include "../Abs/Kinematics.hpp"

template <typename T_out, typename T_in, typename T_arg> class C_OneLoopIntegrator{
	public:
    T_arg integrand_args;
	int k_col,Int_counter;
	C_Kinematics_1loop Momenta;
	C_Integrator_Line<T_out,C_OneLoopIntegrator,T_in> * Integ_radial;
	C_Integrator_Line<T_out,C_OneLoopIntegrator,T_in> * Integ_angle_Z;
	C_Integrator_Line<T_out,C_OneLoopIntegrator,T_in> * Integ_angle_Y;

    std::function<t_cmplxMatrix(t_cmplxArray1D)>  *integrand;

	bool flag_sigma;

	C_OneLoopIntegrator(){
		//integrand_args.resize(3);
	}

	//______________________________________________________________________
	// 3d-dimensional integration routine
	T_out MultiDimInt3D(std::function<t_cmplxMatrix(t_cmplxArray1D)>  *func_to_int){
		integrand=func_to_int;
		k_col=0;
		Int_counter=0;
		return Integ_radial->getResult(&C_OneLoopIntegrator::f1_3D,this);
	}
	T_out f1_3D(T_in modul_momenta){
		k_col++;
		integrand_args[0]=(modul_momenta);
		return Integ_angle_Z->getResult(&C_OneLoopIntegrator::f2_3D,this);
	}
	T_out f2_3D(T_in angle_z){
		integrand_args[1]= angle_z;
		flag_sigma=true;
		return Integ_angle_Y->getResult(&C_OneLoopIntegrator::f3_3D,this);
	}
	T_out f3_3D(T_in angle_y){
		integrand_args[2]= angle_y;
		return (*integrand)(integrand_args);
	}


    //______________________________________________________________________
    // 2d-dimensional integration routine
    T_out MultiDimInt2D(std::function<t_cmplxMatrix(t_cmplxArray1D)>  *func_to_int){
        integrand=func_to_int;
        k_col=0;
        Int_counter=0;
        return Integ_radial->getResult(&C_OneLoopIntegrator::f1_2D,this);
    }
    T_out f1_2D(T_in modul_momenta){
        k_col++;
        integrand_args[0]=(modul_momenta);
        return Integ_angle_Z->getResult(&C_OneLoopIntegrator::f2_2D,this);
    }
    T_out f2_2D(T_in angle_z){
        integrand_args[1]= angle_z;
        return (*integrand)(integrand_args);
    }
};



#endif /* ONELOOPINTEGRATOR_HPP_ */
