//
// Created by stkubr on 01.05.15.
//

#include "BSE_Vector.h"

using namespace BSE;

void C_BSE_Vector::setDiracStructures(t_cmplxVector _k, t_cmplxVector _P, std::vector<t_cmplxDirac> &DiracBasis)  {
    t_cmplxDirac Y_T = C_Kinematics_1loop::TransIn(Y, _P);
    t_cmplxVector _k_T;
    _k_T = C_Kinematics_1loop::TransIn(_k, _P);
    t_cmplx _k2_T = _k_T * _k_T;
    t_cmplx _k2 = _k * _k;
    t_cmplx _k_P = _k * _P;
    t_cmplx _P2 = _P * _P;

    (DiracBasis)[0] = Y_T;
    (DiracBasis)[1] = (_k_T * (_k_T * Z) - 1.0 / 3.0 * Y_T * (_k2_T)) / _k2;
    (DiracBasis)[2] = _k_T * (_P * Z) * _k_P / (_k2 * _P2);
    (DiracBasis)[3] = -1.0 * (Y_T * ((_P * Z) * (_k * Z) - (_k * Z) * (_P * Z)) + 2.0 * _k_T * (_P * Z)) / _k2 / 2.0;
    (DiracBasis)[4] = ii * _k_T * I / _k2;
    (DiracBasis)[5] = ii * (Y_T * (_k_T * Z) - (_k_T * Z) * Y_T) * _k_P / _k2;
    (DiracBasis)[6] = ii * (Y_T * (_P * Z) - (_P * Z) * Y_T) * (1.0 - _k_P * _k_P / _k2 / _P2) -
                      ii * 2.0 * _k_T * (_k_T * Z) * (_P * Z) / _k2;
    (DiracBasis)[7] = ii * _k_T * (_k_T * Z) * (_P * Z) / _k2;
}

t_cmplxMatrix C_BSE_Vector::disentangleAmps(t_cmplxMatrix *pre_result)  {
    t_cmplxMatrix result(num_amplitudes, 1);
    for (int i = 0; i < num_amplitudes; i++)
        result(i, 0) = (*pre_result)(i, 0) / threadloc_WeightCoeff[omp_get_thread_num()][i];
    return result;
}