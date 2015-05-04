//
// Created by stkubr on 29.04.15.
//

#include "BSE_Matrix.h"

using namespace BSE;

t_cmplxMatrix C_BSE_Matrix::integrateYAngle(std::function<t_cmplxMatrix(double)> &bound_member_fn)  {
    return Integrator_angle_Y->getResult(&bound_member_fn, num_amplitudes, num_amplitudes);
}

t_cmplxMatrix C_BSE_Matrix::wrapIntegrand(t_cmplxArray1D *args, double y)  {
    (*args)[2] = y;
    return Integrand_matrix(args);
}

t_cmplxMatrix C_BSE_Matrix::Integrand_matrix(t_cmplxArray1D *args)  {
    t_cmplxMatrix result(num_amplitudes, 1), pre_result(num_amplitudes, 1);
    t_cmplx Kinematic_factor;
    t_cmplx x, y, z;
    x = sqrt((*args)[0]);
    z = (*args)[1];
    y = (*args)[2];
    C_Kinematics_1loop momenta = IntMomenta(x, y, z);
    Kinematic_factor = -1.0 / (8.0 * pi * pi * pi * pi) * ((*args)[0]);
    std::vector<t_cmplxDirac> WaveFunctions = getWaveFunctions(momenta);
    pre_result = traceWithKernel_Matrix(WaveFunctions, momenta);
    result = Kinematic_factor * pre_result;
    threadloc_Integ_ctr[omp_get_thread_num()]++;
    return result;
}

t_cmplxMatrix C_BSE_Matrix::traceWithKernel_Matrix(std::vector<t_cmplxDirac> &WaveFunctions,
                                                   C_Kinematics_1loop &momenta) {
    t_cmplxMatrix result(num_amplitudes, 1), pre_result(num_amplitudes, 1), result_M(num_amplitudes,
                                                                                     num_amplitudes);
    t_cmplxVector k_p_P;
    k_p_P = (momenta.k + momenta.p - momenta.P) / 2.0;
    bool flag_reset_kernel = true;
    for (int i = 0; i < num_amplitudes; i++) {
        for (int j = 0; j < num_amplitudes; j++) {
            pre_result(j, 0) = Kernel->TraceKernelWithoutStoring(threadloc_Projectors[omp_get_thread_num()][j],
                                                                 WaveFunctions[i],
                                                                 momenta.q, momenta.k, k_p_P,
                                                                 flag_reset_kernel);
            flag_reset_kernel = false;
        }
        result = disentangleAmps(&pre_result);
        for (int j = 0; j < num_amplitudes; j++) {
            result_M(i, j) = result(j, 0);
        }
    }
    return result_M;
}

void C_BSE_Matrix::calcEVMatrix(t_cmplx _P)  {
    threadloc_Momenta[omp_get_thread_num()].SetVector_P(_P);
    PreCalculation();
    Memory->ResizeEVMatrix(params.NumRadial, params.NumCheb_nod1, num_amplitudes, 1);
#pragma omp parallel
    {//start of pragma
        t_cmplxArray1D integrand_args_local;
        integrand_args_local.resize(numIntegDimentions);
        threadloc_Momenta[omp_get_thread_num()].SetVector_P(_P);
        threadloc_Integ_ctr[omp_get_thread_num()] = 0;
        std::function<t_cmplxMatrix(double)> bound_member_fn =
                std::bind(&C_BSE_Matrix::wrapIntegrand, this, &integrand_args_local, std::placeholders::_1);
#pragma omp for
        for (int p2_inx = 1; p2_inx < zz_rad.size(); p2_inx++) {
            double p2 = zz_rad[p2_inx];
            for (int z_p_inx = 1; z_p_inx < zz_cheb.size(); z_p_inx++) {
                double z_p = zz_cheb[z_p_inx];
                threadloc_Momenta[omp_get_thread_num()].SetVectors_p(z_p, sqrt(p2));
                setDiracStructures(threadloc_Momenta[omp_get_thread_num()].p,
                                   threadloc_Momenta[omp_get_thread_num()].P,
                                   threadloc_Projectors[omp_get_thread_num()]);
                setWeightCoeff();

                setInternalGrid(bound_member_fn, &integrand_args_local, p2_inx, z_p_inx);
                threadloc_Integ_ctr[omp_get_thread_num()] = 0;
            }
        }
    }// end of pragma
    flag_precalculation = false;
}

void C_BSE_Matrix::setInternalGrid(std::function<t_cmplxMatrix(double)> &bound_member_fn,
                                   t_cmplxArray1D *integrand_args_local, int p2_inx, int z_p_inx)  {
    t_cmplxMatrix Temp_Matrix;
    for (int k2_inx = 1; k2_inx < zz_rad.size(); k2_inx++) {
        (*integrand_args_local)[0] = zz_rad[k2_inx];
        double w_k2 = w_rad[k2_inx];
        for (int z_k_inx = 1; z_k_inx < zz_cheb.size(); z_k_inx++) {
            (*integrand_args_local)[1] = zz_cheb[z_k_inx];
            double w_z_k = w_cheb[z_k_inx];
            Temp_Matrix = w_k2 * w_z_k * integrateYAngle(bound_member_fn);
            setEVMatrix(Temp_Matrix, p2_inx, z_p_inx, k2_inx, z_k_inx);
        }
        threadloc_momentum_inx[omp_get_thread_num()]++;
    }
}

void C_BSE_Matrix::setEVMatrix(t_cmplxMatrix &Temp_Matrix,
                               int p2_inx, int z_p_inx, int k2_inx, int z_k_inx)  {
    for (int proj_amp_inx = 0; proj_amp_inx < num_amplitudes; proj_amp_inx++) {
        for (int wave_amp_inx = 0; wave_amp_inx < num_amplitudes; wave_amp_inx++) {
            int inx_external = z_p_inx - 1 + (p2_inx - 1) * (params.NumCheb_nod1) +
                               params.NumRadial * (params.NumCheb_nod1) * proj_amp_inx;
            int inx_internal = z_k_inx - 1 + (k2_inx - 1) * (params.NumCheb_nod1) +
                               params.NumRadial * (params.NumCheb_nod1) * wave_amp_inx;
            Memory->EVMatrix(inx_external, inx_internal) = pi / 2.0 * Temp_Matrix(proj_amp_inx, wave_amp_inx);
        }
    }
}

t_cmplx C_BSE_Matrix::detectSymmetricity(Eigen::MatrixXcf &eigenvectors, int num_state)  {
    t_cmplx symmetricity;
    int p_ctr = 1;
    int p_amp_ctr = 0;
    int inx_momenta_amp =
            (p_ctr - 1) * (params.NumCheb_nod1) + params.NumRadial * (params.NumCheb_nod1) * p_amp_ctr;
    for (int zp_ctr = 1; zp_ctr < zz_cheb.size(); zp_ctr++) {
        t_cmplx summ, diff;
        int inx_angle_plus = zp_ctr - 1 + inx_momenta_amp;
        int inx_angle_minus = zz_cheb.size() - zp_ctr - 1 + inx_momenta_amp;
        summ = eigenvectors.col(num_state)(inx_angle_plus) + eigenvectors.col(num_state)(inx_angle_minus);
        diff = eigenvectors.col(num_state)(inx_angle_plus) - eigenvectors.col(num_state)(inx_angle_minus);
        if (fabs(real(summ)) <= fabs(real(diff))) { symmetricity = -1.0; }
        else { symmetricity = 1.0; }
    }
    return symmetricity;
}

t_cmplxArray2D C_BSE_Matrix::calcEigenStates(t_cmplxVector P, int numberOfStates)  {
    preparePropagators();
    calcEVMatrix(P[3]);
    std::cout << "EVMatrix is full. EigenValues computation engaged..." << std::endl;
    Eigen::ComplexEigenSolver<Eigen::MatrixXcf> ces;
    ces.compute(Memory->EVMatrix);
    Eigen::VectorXcf eigenvalues = ces.eigenvalues();
    Eigen::MatrixXcf eigenvectors = (Eigen::MatrixXcf) ces.eigenvectors();
    std::cout << "EigenValues computation is done. The eigenvalues of EVMatrix are obtained." << std::endl;
    int i = eigenvalues.size() - 1;
    t_cmplxArray2D Dominant_EV_and_parity(2);
    while (i > eigenvalues.size() - numberOfStates) {
        t_cmplx symmetricity = detectSymmetricity(eigenvectors, i);
        Dominant_EV_and_parity[0].push_back(eigenvalues[i]);
        Dominant_EV_and_parity[1].push_back(symmetricity);
        std::cout << i << "  " << eigenvalues[i] << "  " << symmetricity << std::endl;
        i--;
    }
    return Dominant_EV_and_parity;
}

double C_BSE_Matrix::checkSum_EVMatrixNorm()  {
    calcEVMatrix(t_cmplx(0, 0.1));
    return Memory->EVMatrix.norm();
}
