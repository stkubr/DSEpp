/*
 * OneLoopIntegrator.hpp
 *
 *  Created on: Feb 12, 2014
 *      Author: stkubr
 */

#ifndef ONELOOPINTEGRATOR_HPP_
#define ONELOOPINTEGRATOR_HPP_

#include <functional>
#include <omp.h>
#include "../Abs/Kinematics.hpp"
#include "Integrator.hpp"

namespace Integration {

    template<typename T_out, typename T_in, typename T_arg>
    class C_OneLoopIntegrator {
    protected:
        int numIntegDimentions;
        C_Integrator_Line<T_out, T_in> *Integrator_momentum;
        C_Integrator_Line<T_out, T_in> *Integrator_angle_Z;
        C_Integrator_Line<T_out, T_in> *Integrator_angle_Y;
        std::vector<int> threadloc_Integ_ctr;
        std::vector<int> threadloc_momentum_inx;

        /// the variables over the loop integration is done, for example vector \f$ (k,z,y,\phi) \f$
        T_arg integrand_args;

        /// the momenta 4-vectors associated with singe loop integration
        std::vector<C_Kinematics_1loop> threadloc_Momenta;

        C_OneLoopIntegrator() {
            threadloc_Integ_ctr.resize(omp_get_max_threads(), 0);
            threadloc_momentum_inx.resize(omp_get_max_threads(), 0);
            threadloc_Momenta.resize(omp_get_max_threads());
        }

        virtual ~C_OneLoopIntegrator() { }

        T_out calcIntegral3D(std::function<T_out(T_arg)> *func_to_int, int numRows, int numCols) {
            t_dArray1D x, z, y;
            t_dArray1D w_x, w_z, w_y;
            T_arg integrand_args_local;
            Integrator_momentum->getNodes(x, w_x);
            Integrator_angle_Z->getNodes(z, w_z);
            Integrator_angle_Y->getNodes(y, w_y);
            T_out result, sum;
            integrand_args_local.resize(numIntegDimentions);
            sum.Resize(numRows, numCols);
            double w_x_temp, w_z_temp;
            threadloc_Integ_ctr[omp_get_thread_num()] = 0;
            for (int i = 1; i < x.size(); ++i) {
                threadloc_momentum_inx[omp_get_thread_num()] = i;
                integrand_args_local[0] = x[i];
                w_x_temp = w_x[i];
                for (int j = 1; j < z.size(); ++j) {
                    integrand_args_local[1] = z[j];
                    w_z_temp = w_z[j];
                    for (int k = 1; k < y.size(); ++k) {
                        integrand_args_local[2] = y[k];
                        result = w_x_temp * w_z_temp * w_y[k] * (*func_to_int)(integrand_args_local);
                        threadloc_Integ_ctr[omp_get_thread_num()]++;
                        sum += result;
                    }
                }
            }
            return sum;
        }

        T_out calcIntegral2D(std::function<T_out(T_arg)> *func_to_int, int numRows, int numCols) {
            t_dArray1D x0, x1;
            t_dArray1D w0, w1;
            T_arg integrand_args_local;
            Integrator_momentum->getNodes(x0, w0);
            Integrator_angle_Z->getNodes(x1, w1);
            T_out result, sum;
            integrand_args_local.resize(numIntegDimentions);
            sum.Resize(numRows, numCols);
            double w1_temp;
            for (int i = 1; i < x0.size(); ++i) {
                integrand_args_local[0] = x0[i];
                w1_temp = w0[i];
                for (int j = 1; j < x1.size(); ++j) {
                    integrand_args_local[1] = x1[j];
                    result = w1_temp * w1[j] * (*func_to_int)(integrand_args_local);
                    sum += result;
                }
            }
            return sum;
        }
    };

}
#endif /* ONELOOPINTEGRATOR_HPP_ */
