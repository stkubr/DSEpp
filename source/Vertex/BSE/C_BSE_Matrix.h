//
// Created by stkubr on 29.04.15.
//

#ifndef DSEPP_C_BSE_MATRIX_H
#define DSEPP_C_BSE_MATRIX_H


#include "BSE_TwoBody.h"

class C_BSE_Matrix: public C_BSE_TwoBody {
public:

    t_cmplxMatrix IntegAngleY(std::function<t_cmplxMatrix(double)> & bound_member_fn){
        return Integrator_angle_Y->getResult(&bound_member_fn, num_amplitudes, num_amplitudes);
    }

    t_cmplxMatrix func_warp(t_cmplxArray1D & args, double y){
        args[2]=y;
        return Integrand_matrix(args);
    }

    t_cmplxMatrix Integrand_matrix(t_cmplxArray1D args){
        t_cmplxMatrix result(num_amplitudes,1), pre_result(num_amplitudes,1);
        t_cmplx Kinematic_factor;
        t_cmplx x,y,z;
        x=sqrt(args[0]);
        z=args[1];
        y=args[2];
        C_Kinematics_1loop momenta = SetIntMomenta(x,y,z);
        t_cmplxArray1D U_amp = SetDressing_normal(z);
        Kinematic_factor=-1.0/(8.0*pi*pi*pi*pi)*(args[0]);
        std::vector<t_cmplxDirac> WaveFunctions = SetWaveFunctions(momenta);
        t_cmplxDirac FullWaveFunction = SetFullWaveFunction(WaveFunctions, U_amp);
        pre_result = traceWithKernel_Matrix(FullWaveFunction, momenta);
        result = Kinematic_factor*pre_result;
        return result;
    }

    t_cmplxMatrix traceWithKernel_Matrix(t_cmplxDirac & FullWaveFunction,C_Kinematics_1loop & momenta){
        t_cmplxMatrix result(num_amplitudes,1), pre_result(num_amplitudes,1), result_M(num_amplitudes,num_amplitudes);
        t_cmplxVector k_p_P;
        k_p_P = (momenta.k + momenta.p - momenta.P)/2.0;
        bool flag_reset_kernel=true;
        for (int i = 0; i < num_amplitudes; i++){
            for (int j = 0; j < num_amplitudes; j++){
            pre_result(i,0)=Kernel->TraceKernelWithoutStoring(threadloc_Projectors[omp_get_thread_num()][i],
                                                              FullWaveFunction,
                                                              momenta.q, momenta.k, k_p_P,
                                                              flag_reset_kernel);
            flag_reset_kernel=false;
            }
            result=DisentangleAmps(&pre_result);
            for (int j = 0; j < num_amplitudes; j++){
                result_M(i,j)=result(j,0);
            }
        }
        return result;
    }

    t_cmplxArray2D SetEVMatrix(t_cmplx _P){
        threadloc_Momenta[omp_get_thread_num()].SetVector_P(_P);
        PreCalculation();
        t_dArray1D zz_rad_temp,zz_cheb_temp;
        setInitialAMP();
        Eigen::ComplexEigenSolver<Eigen::MatrixXcf> ces;
        Memory->ResizeEVMatrix(params.NumRadial,params.NumCheb_nod1,num_amplitudes,1);
        //std::cout << "Memory->EVMatrix call" << std::endl;
//#pragma omp parallel
        {//start of pragma
            t_cmplxMatrix Temp_Matrix;
            double _p2,_k2,_z,_y,_zp;
            double _w_zp,_w_k2,_w_z;
            t_cmplxArray1D integrand_args_local;
            integrand_args_local.resize(numIntegDimentions);

            threadloc_momentum_inx[omp_get_thread_num()]=0;
            threadloc_Integ_ctr[omp_get_thread_num()]=0;
            std::function<t_cmplxMatrix(double)> bound_member_fn =
                    std::bind(&C_BSE_Matrix::func_warp, this, integrand_args_local, std::placeholders::_1);
//#pragma omp for
            for (int p_ctr = 1; p_ctr < zz_rad.size() ; p_ctr++){
                threadloc_momentum_inx[omp_get_thread_num()]=p_ctr;
                _p2=zz_rad[p_ctr];
                for (int zp_ctr = 1; zp_ctr < zz_cheb.size() ; zp_ctr++){
                    _zp=zz_cheb[zp_ctr];
                    threadloc_Momenta[omp_get_thread_num()].SetVectors_p(_zp,sqrt(_p2));
                    SetDiracStructures(threadloc_Momenta[omp_get_thread_num()].p,
                                       threadloc_Momenta[omp_get_thread_num()].P,
                                       threadloc_Projectors[omp_get_thread_num()]);
                    SetWeightCoeff();
                    for (int k_ctr = 1; k_ctr < zz_rad.size(); k_ctr++){
                        _k2=zz_rad[k_ctr];
                        _w_k2=w_rad[k_ctr];
                        integrand_args_local[0]=_k2;
                        for (int z_ctr = 1; z_ctr < zz_cheb.size() ; z_ctr++){
                            _z=zz_cheb[z_ctr];
                            _w_z=w_cheb[z_ctr];
                            integrand_args_local[1]=_z;
                            Temp_Matrix=IntegAngleY(bound_member_fn);
                            for (int p_amp_ctr = 0; p_amp_ctr < num_amplitudes ; p_amp_ctr++){
                                for (int k_amp_ctr = 0; k_amp_ctr < num_amplitudes ; k_amp_ctr++){

                                    //std::cout << k_ctr << "  " << z_ctr << "  " << z_ctr-1 + (k_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*k_amp_ctr << std::endl;
                                    //std::cout << k_ctr << "  " << z_ctr << "  " << zp_ctr + (p_ctr-1)*(NumRadialPoints-z_ctr*(NumRadialPoints-1))-1+ NumRadialPoints*(num_cheb_nod1)*p_amp_ctr << "  " << z_ctr + (k_ctr-1)*(NumRadialPoints-z_ctr*(NumRadialPoints-1)) + NumRadialPoints*(num_cheb_nod1)*k_amp_ctr-1 << std::endl;
                                    Memory->EVMatrix(zp_ctr-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*p_amp_ctr, z_ctr-1 + (k_ctr-1)*(params.NumCheb_nod1) + params.NumRadial*(params.NumCheb_nod1)*k_amp_ctr)=pi/2.0*_w_k2*_w_z*Temp_Matrix(p_amp_ctr,k_amp_ctr);
                                }
                            }
                            //std::cout << Temp_Matrix(0,0) << std::endl;
                            //cin.get();
                        }
                        threadloc_momentum_inx[omp_get_thread_num()]++;
                    }
                    threadloc_momentum_inx[omp_get_thread_num()]=0;
                    threadloc_Integ_ctr[omp_get_thread_num()]=0;
                }
            }
        }// end of pragma
        //std::cout << "EVMatrix is full. EigenValues computation engaged..." << std::endl;
        ces.compute(Memory->EVMatrix);
        flag_precalculation=false;
        //CalcEVMatrix(&ces);
        Eigen::VectorXcf EV=ces.eigenvalues();


        auto Parity = [&](int num_state) -> t_cmplx {
            t_cmplx parity=0.0;
            for (int p_ctr = 1; p_ctr < 2 ; p_ctr++){
                for (int p_amp_ctr = 0; p_amp_ctr < 1 ; p_amp_ctr++){
                    for (int zp_ctr = 1; zp_ctr < zz_cheb.size() ; zp_ctr++){
                        t_cmplx summ,diff;

                        summ = ces.eigenvectors().col(num_state)(zp_ctr-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*p_amp_ctr) +
                               ces.eigenvectors().col(num_state)(zz_cheb.size() - zp_ctr-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*p_amp_ctr);

                        diff = ces.eigenvectors().col(num_state)(zp_ctr-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*p_amp_ctr) -
                               ces.eigenvectors().col(num_state)(zz_cheb.size() - zp_ctr-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*p_amp_ctr);

                        if(fabs(real(summ))<=fabs(real(diff))) {parity = -1.0;}
                        else{ parity = 1.0; }
                    }
                }
            }
            return parity;
        };

        //std::cout << "EigenValues computation is done. The eigenvalues of EVMatrix are obtained." << std::endl;
        int i= EV.size()-1;
        t_cmplxArray2D Dominant_EV_and_parity(2);
        while( i > EV.size()-10)
        {
            //std::cout << EV[i] << "  " << i << std::endl;
            Dominant_EV_and_parity[0].push_back(EV[i]);
            Dominant_EV_and_parity[1].push_back(Parity(i));
            std::cout << i << "  " << EV[i] << "  " << Parity(i)<< std::endl;
            i--;
        }
        //std::cout << ces.eigenvectors().col(EV.size()-1) << std::endl;
        return Dominant_EV_and_parity;
    }

  /*  void DrawBSA_matrix(t_cmplx _P, int _state, int amp_num){
        int num_state;
        Momenta.SetVector_P(_P);
        PreCalculation();
        t_dArray1D zz_rad_temp,zz_cheb_temp;
        SetDressing_ref=&C_BSE_Hadron_Base::SetDressing_normal;
        GetBSA_ref=&C_BSE_Hadron_Base::GetBSA_matrix;

        flag_precalculation=false;
        Eigen::ComplexEigenSolver<Eigen::MatrixXcf> ces;
        //ces.compute(Memory->EVMatrix);
        CalcEVMatrix(&ces);
        Eigen::VectorXcf EV=ces.eigenvalues();
        num_state = EV.size()- _state;

        ofstream DrawBSA_matrix;
        DrawBSA_matrix.open ("Data_files/BSEs_matrix.dat");

        int i=EV.size()-1;

        while(i > EV.size()-6)
        {
            //std::cout << EV[i] << "  " << i << std::endl;
            i--;
        }

        t_cmplxArray1D result(zz_rad.size());
        for (int p_ctr = 1; p_ctr <= params.NumRadial ; p_ctr++){
            result[p_ctr-1]=0;
            for (int j=1;j<=params.NumCheb_nod1;j++)
            {
                t_cmplx z_v,element;
                element=ces.eigenvectors().col(num_state)(j-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*amp_num);
                z_v=zz_cheb[j];
                result[p_ctr-1] += w_cheb[j]*Cheb_polinoms(z_v,0)*real(element);
                //std::cout << w_cheb[j] << "  " << U_ex(z_v) << "  " << real(element) << std::endl;
            }
            DrawBSA_matrix << zz_rad[p_ctr] << "  " << real(result[p_ctr-1]) << std::endl;
            //std::cout << zz_rad[p_ctr] << "  " << result[p_ctr-1] << std::endl;
            //cin.get();
        }
        DrawBSA_matrix.close();
    }

    */

};


#endif //DSEPP_C_BSE_MATRIX_H
