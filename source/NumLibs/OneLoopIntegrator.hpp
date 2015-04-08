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
	int num_IntegDimentions;
	C_Kinematics_1loop Momenta;
	C_Integrator_Line<T_out,T_in> * Integrator_momentum;
	C_Integrator_Line<T_out,T_in> * Integrator_angle_Z;
	C_Integrator_Line<T_out,T_in> * Integrator_angle_Y;

    std::function<T_out(T_arg)>  *integrand;

	/*C_OneLoopIntegrator(){
        integrand_args.resize(3);
	}*/

	//______________________________________________________________________
	// 3d-dimensional integration routine
/*	T_out MultiDimInt3D(std::function<T_out(T_arg)>  *func_to_int){
		integrand=func_to_int;
		return Integrator_momentum->getResult(&C_OneLoopIntegrator::f1_3D,this);
	}
	T_out f1_3D(T_in modul_momenta){
		integrand_args[0]=(modul_momenta);
		return Integrator_angle_Z->getResult(&C_OneLoopIntegrator::f2_3D,this);
	}
	T_out f2_3D(T_in angle_z){
		integrand_args[1]= angle_z;
		return Integrator_angle_Y->getResult(&C_OneLoopIntegrator::f3_3D,this);
	}
	T_out f3_3D(T_in angle_y){
		integrand_args[2]= angle_y;
		return (*integrand)(integrand_args);
	}
*/

    //______________________________________________________________________
    // 2d-dimensional integration routine
 /*   T_out MultiDimInt2D(std::function<T_out(T_arg)>  *func_to_int){
        integrand=func_to_int;
        return Integrator_momentum->getResult(&C_OneLoopIntegrator::f1_2D,this);
    }
    T_out f1_2D(T_in modul_momenta){
        integrand_args[0]=(modul_momenta);
        return Integrator_angle_Z->getResult(&C_OneLoopIntegrator::f2_2D,this);
    }
    T_out f2_2D(T_in angle_z){
        integrand_args[1]= angle_z;
        return (*integrand)(integrand_args);
    }
*/

    T_out MultiDimInt2D_wo_nested(std::function<T_out(T_arg)> *func_to_int, int numRows, int numCols){
        t_dArray1D x0,x1;
        t_dArray1D w0,w1;
        T_arg integrand_args_local;
        Integrator_momentum->getNodes(x0, w0);
        Integrator_angle_Z->getNodes(x1, w1);
        integrand=func_to_int;
        T_out result,sum;
        integrand_args_local.resize(num_IntegDimentions);
        sum.Resize(numRows,numCols);
        double w1_temp;
        for (int i = 1; i < x0.size(); ++i) {
            integrand_args_local[0]=x0[i];
            w1_temp = w0[i];
            for (int j = 1; j < x1.size(); ++j) {
                integrand_args_local[1]=x1[j];
                result = w1_temp*w1[j]*(*integrand)(integrand_args_local);
                sum += result;
            }
        }
        return sum;
    }


    //todo implement conversion from multidimensional to one dimension non-recursive summation
/*    T_out MultiDimND(std::function<T_out(T_arg)>  *func_to_int){
        t_dArray2D coordinates(integrand_args.size());
        t_dArray2D weights(integrand_args.size());
        for (int i = 0; i < integrand_args.size(); ++i) {
            Integrator[i]->getNotes(coordinates[i],weights[i]);
        }

    }
*/


};



#endif /* ONELOOPINTEGRATOR_HPP_ */
