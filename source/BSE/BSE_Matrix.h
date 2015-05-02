//
// Created by stkubr on 29.04.15.
//

#ifndef DSEPP_C_BSE_MATRIX_H
#define DSEPP_C_BSE_MATRIX_H

#include "BSE_TwoBody.h"

class C_BSE_Matrix: public C_BSE_TwoBody {
protected:
    virtual ~C_BSE_Matrix(){}

    t_cmplxMatrix IntegAngleY(std::function<t_cmplxMatrix(double)> & bound_member_fn){
        return Integrator_angle_Y->getResult(&bound_member_fn, num_amplitudes, num_amplitudes);
    }

    t_cmplxMatrix func_warp(t_cmplxArray1D * args, double y){
        (*args)[2]=y;
        return Integrand_matrix(args);
    }

    t_cmplxMatrix Integrand_matrix(t_cmplxArray1D * args){
        t_cmplxMatrix result(num_amplitudes,1), pre_result(num_amplitudes,1);
        t_cmplx Kinematic_factor;
        t_cmplx x,y,z;
        x=sqrt((*args)[0]);
        z=(*args)[1];
        y=(*args)[2];
        C_Kinematics_1loop momenta = SetIntMomenta(x,y,z);
        Kinematic_factor=-1.0/(8.0*pi*pi*pi*pi)*((*args)[0]);
        std::vector<t_cmplxDirac> WaveFunctions = SetWaveFunctions(momenta);
        pre_result = traceWithKernel_Matrix(WaveFunctions, momenta);
        result = Kinematic_factor*pre_result;
        threadloc_Integ_ctr[omp_get_thread_num()]++;
        return result;
    }

    t_cmplxMatrix traceWithKernel_Matrix(std::vector<t_cmplxDirac> & WaveFunctions,C_Kinematics_1loop & momenta){
        t_cmplxMatrix result(num_amplitudes,1), pre_result(num_amplitudes,1), result_M(num_amplitudes,num_amplitudes);
        t_cmplxVector k_p_P;
        k_p_P = (momenta.k + momenta.p - momenta.P)/2.0;
        bool flag_reset_kernel=true;
        for (int i = 0; i < num_amplitudes; i++){
            for (int j = 0; j < num_amplitudes; j++){
                pre_result(j,0)=Kernel->TraceKernelWithoutStoring(threadloc_Projectors[omp_get_thread_num()][j],
                                                                  WaveFunctions[i],
                                                                  momenta.q, momenta.k, k_p_P,
                                                                  flag_reset_kernel);
                flag_reset_kernel=false;
            }
            result=DisentangleAmps(&pre_result);
            for (int j = 0; j < num_amplitudes; j++){
                result_M(i,j)=result(j,0);
            }
        }
        return result_M;
    }

    void calcEVMatrix(t_cmplx _P){
        threadloc_Momenta[omp_get_thread_num()].SetVector_P(_P);
        PreCalculation();
        Memory->ResizeEVMatrix(params.NumRadial,params.NumCheb_nod1,num_amplitudes,1);
#pragma omp parallel
        {//start of pragma
            t_cmplxArray1D integrand_args_local;
            integrand_args_local.resize(numIntegDimentions);
            threadloc_Momenta[omp_get_thread_num()].SetVector_P(_P);
            threadloc_Integ_ctr[omp_get_thread_num()]=0;
            std::function<t_cmplxMatrix(double)> bound_member_fn =
                    std::bind(&C_BSE_Matrix::func_warp, this, &integrand_args_local, std::placeholders::_1);
#pragma omp for
            for (int p_ctr = 1; p_ctr < zz_rad.size() ; p_ctr++){
                double _p2=zz_rad[p_ctr];
                for (int zp_ctr = 1; zp_ctr < zz_cheb.size() ; zp_ctr++){
                    double _zp=zz_cheb[zp_ctr];
                    threadloc_Momenta[omp_get_thread_num()].SetVectors_p(_zp,sqrt(_p2));
                    SetDiracStructures(threadloc_Momenta[omp_get_thread_num()].p,
                                       threadloc_Momenta[omp_get_thread_num()].P,
                                       threadloc_Projectors[omp_get_thread_num()]);
                    SetWeightCoeff();
                    setInternalGrid(bound_member_fn, &integrand_args_local, p_ctr, zp_ctr);
                    threadloc_Integ_ctr[omp_get_thread_num()]=0;
                }
            }
        }// end of pragma
        flag_precalculation=false;
    }

    void setInternalGrid(std::function<t_cmplxMatrix(double)> & bound_member_fn,
                         t_cmplxArray1D * integrand_args_local,
                         int p_ctr, int zp_ctr){
        t_cmplxMatrix Temp_Matrix;
        for (int k_ctr = 1; k_ctr < zz_rad.size(); k_ctr++){
            (*integrand_args_local)[0]=zz_rad[k_ctr];
            double _w_k2=w_rad[k_ctr];
            for (int z_ctr = 1; z_ctr < zz_cheb.size() ; z_ctr++){
                (*integrand_args_local)[1]=zz_cheb[z_ctr];
                double _w_z=w_cheb[z_ctr];
                Temp_Matrix=_w_k2*_w_z*IntegAngleY(bound_member_fn);
                setEVMatrix(Temp_Matrix, p_ctr, zp_ctr, k_ctr, z_ctr);
            }
            threadloc_momentum_inx[omp_get_thread_num()]++;
        }
    }

    void setEVMatrix(t_cmplxMatrix & Temp_Matrix, int p_ctr, int zp_ctr, int k_ctr, int z_ctr){
        for (int p_amp_ctr = 0; p_amp_ctr < num_amplitudes ; p_amp_ctr++){
            for (int k_amp_ctr = 0; k_amp_ctr < num_amplitudes ; k_amp_ctr++){
                int inx_external = zp_ctr-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*p_amp_ctr;
                int inx_internal = z_ctr-1 + (k_ctr-1)*(params.NumCheb_nod1) + params.NumRadial*(params.NumCheb_nod1)*k_amp_ctr;
                Memory->EVMatrix(inx_external, inx_internal) = pi/2.0*Temp_Matrix(p_amp_ctr,k_amp_ctr);
            }
        }
    }

    t_cmplx detectSymmetricity(Eigen::MatrixXcf & eigenvectors, int num_state){
        t_cmplx symmetricity;
        int p_ctr = 1;
        int p_amp_ctr = 0;
        int inx_momenta_amp = (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*p_amp_ctr;
        for (int zp_ctr = 1; zp_ctr < zz_cheb.size() ; zp_ctr++){
            t_cmplx summ,diff;
            int inx_angle_plus =  zp_ctr-1 + inx_momenta_amp;
            int inx_angle_minus =  zz_cheb.size() - zp_ctr-1 + inx_momenta_amp;
            summ = eigenvectors.col(num_state)(inx_angle_plus) + eigenvectors.col(num_state)(inx_angle_minus);
            diff = eigenvectors.col(num_state)(inx_angle_plus) - eigenvectors.col(num_state)(inx_angle_minus);
            if(fabs(real(summ))<=fabs(real(diff))) {symmetricity = -1.0;}
            else{ symmetricity = 1.0; }
        }
        return symmetricity;
    }

public:
    t_cmplxArray2D calcEigenStates(t_cmplxVector P, int numberOfStates){
        preparePropagators();
        calcEVMatrix(P[3]);
        std::cout << "EVMatrix is full. EigenValues computation engaged..." << std::endl;
        Eigen::ComplexEigenSolver<Eigen::MatrixXcf> ces;
        ces.compute(Memory->EVMatrix);
        Eigen::VectorXcf eigenvalues = ces.eigenvalues();
        Eigen::MatrixXcf eigenvectors = (Eigen::MatrixXcf)ces.eigenvectors();
        std::cout << "EigenValues computation is done. The eigenvalues of EVMatrix are obtained." << std::endl;
        int i= eigenvalues.size()-1;
        t_cmplxArray2D Dominant_EV_and_parity(2);
        while( i > eigenvalues.size() - numberOfStates) {
            t_cmplx symmetricity = detectSymmetricity(eigenvectors,i);
            Dominant_EV_and_parity[0].push_back(eigenvalues[i]);
            Dominant_EV_and_parity[1].push_back(symmetricity);
            std::cout << i << "  " << eigenvalues[i] << "  " << symmetricity << std::endl;
            i--;
        }
        return Dominant_EV_and_parity;
    }

    double checkSum_EVMatrixNorm(){
        calcEVMatrix(t_cmplx(0,0.1));
        return Memory->EVMatrix.norm();
    }
};


#endif //DSEPP_C_BSE_MATRIX_H
