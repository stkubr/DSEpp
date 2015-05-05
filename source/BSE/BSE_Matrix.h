//
// Created by stkubr on 29.04.15.
//

#ifndef DSEPP_C_BSE_MATRIX_H
#define DSEPP_C_BSE_MATRIX_H

#include "BSE_TwoBody.h"

namespace BSE {

/**
 * \brief Two-quark Bethe-Salpeter Equations
 *
 * This class adds the lambda eigenvalue calculation by
 * matrix eigenvalue decomposition method. The method itself is provided by Eigen library.
 * So the steps are:
 * - set all momenta, Dirac basis projectors and partons(propagators)
 * - calculate WaveFunction and Projectors for every amplitude permutation
 * - trace them with scattering Kernel on every step of Y angle integration
 * - use the integrated out resulting matrix from only Y angle integration to populate the EVMatrix
 * - do for every point of external (p^2, z_p) grid and internal (k^2, z_k) grid
 * - once the EV_matrix is full use Eigen library to compute eigenvalues
 * - detect the eigenvector Z angle symmetricity
 * - output the eigenvalues and corresponding symmetricity
 */
    class C_BSE_Matrix : public C_BSE_TwoBody {
    protected:
        virtual ~C_BSE_Matrix() { }

        /// integrates over only Y angle
        t_cmplxMatrix integrateYAngle(std::function<t_cmplxMatrix(double)> &bound_member_fn);

        /// wraps the integrand for Y angle only integration
        t_cmplxMatrix wrapIntegrand(t_cmplxArray1D *args, double y);

        /// Loop integrand
        /// "matrix" denote the fact that it outputs the matrix of amplitudes projected out from WaveFunctions
        t_cmplxMatrix Integrand_matrix(t_cmplxArray1D *args);

        /// trace the WaveFunctions and threadloc_Projectors with scattering Kernel
        t_cmplxMatrix traceWithKernel_Matrix(std::vector<t_cmplxDirac> &WaveFunctions, C_Kinematics_1loop &momenta);

        /// calls setInternalGrid() on every point of external grid
        void calcEVMatrix(t_cmplx _P);

        /// calls setEVMatrix() on every point of internal grid
        void setInternalGrid(std::function<t_cmplxMatrix(double)> &bound_member_fn,
                             t_cmplxArray1D *integrand_args_local,
                             int p2_inx, int z_p_inx);

        /// populates EVMatrix from Temp_Matrix according to provided indexes
        void setEVMatrix(t_cmplxMatrix &Temp_Matrix, int p2_inx, int z_p_inx, int k2_inx, int z_k_inx);

        /// checks Z angle symmetricity of eigenvectors
        t_cmplx detectSymmetricity(Eigen::MatrixXcf &eigenvectors, int num_state);

    public:
        /// calculates the first numberOfStates biggest eigenvalues and their symmetricity
        t_cmplxArray2D calcEigenStates(t_cmplxVector P, int numberOfStates);

        /// checks the EVMatrix norm for tests
        double checkSum_EVMatrixNorm();
    };

}

#endif //DSEPP_C_BSE_MATRIX_H
