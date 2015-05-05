//
// Created by stkubr on 21.04.15.
//

#ifndef SHOW_PRECISION
#define SHOW_PRECISION 6
#endif

#ifndef DSEPP_BSE_TWOBODY_H
#define DSEPP_BSE_TWOBODY_H

#include <source/NumLibs/OneLoopIntegrator.hpp>
#include <source/NumLibs/Integrator.hpp>
#include <source/NumLibs/Integrator_Path.hpp>
#include "source/NumLibs/ChebychevPolinoms.h"
#include "source/NumLibs/Extra_functions.h"
#include <iomanip>
#include "BSE.h"
#include "BoundState_parameters.h"

namespace BSE {

/**
 * \brief Two-quark Bethe-Salpeter Equations
 *
 * In this class it uses the Power Method to obtain the solutions
 * overall this method can be described as:
 * - set all momenta, Dirac basis projectors and partons(propagators)
 * - use AmplitudesMatrix as initial guess for the BSAs (Bethe-Salpeter Amplutudes) to construct BSEVertex
 * - calculate FullWaveFunction and trace it and Projectors with scattering Kernel on every step of loop integration
 * - project out the resulting BSAs into buffer_AmplitudesMatrix
 * - calculate Lambda eigenvalue
 * - repeat for desired number of steps
 */
    class C_BSE_TwoBody : public C_BSE, public Integration::C_OneLoopIntegrator<t_cmplxMatrix, double, t_cmplxArray1D> {
    protected:
        /// Pointer to scattering kernel
        Kernels::C_AbstractKernel *Kernel;

        /// Pointer to parton with momenta \f$ k_+ \f$
        Propagators::C_Propagator *Parton_P;

        /// Pointer to parton with momenta \f$ k_- \f$
        Propagators::C_Propagator *Parton_M;

        /// Pointer to external dedicated memory storage
        C_DedicMem_BSE *Memory;

        /// BSE parameters
        C_BoundState_parameters params;

        /// Number of projected out amplitudes (ex:pion has four)
        int num_amplitudes;

        /// threadlocal storage of projectors
        std::vector<std::vector<t_cmplxDirac>> threadloc_Projectors;

        /// threadlocal storage of scalar products of two Projectors of the each amplitude
        t_cmplxArray2D threadloc_WeightCoeff;

        /// Matrix of all dressing functions of all amplitudes (for all Chebychev orders)
        /// Columns denote amplitude number and Chebychev order
        /// Rows denote momenta \f$ p^2 \f$ dependence
        t_cmplxMatrix AmplitudesMatrix;

        t_cmplxMatrix buffer_AmplitudesMatrix, buffer_BSA;

        /// "zz" are arrays of integration notes, "w" are corresponing weights
        t_dArray1D zz_rad, zz_cheb, zz_angleY, w_rad, w_cheb, w_angleY;

        /// flag if precalculation is done
        bool flag_precalculation;

        /// Constructor
        C_BSE_TwoBody();

        /// Destructor
        virtual ~C_BSE_TwoBody();

        /// sets threadloc_WeightCoeff
        void setWeightCoeff();

        /// check Dirac basic Orthogonality
        void checkOrthogonality();

        /// Initialize integrators and resize storages
        void Initialization();

        /// Resize storages and get integration nodes
        void resizeAllStorages();

        /// Sets all enties in AmplitudesMatrix to 1.0
        void setInitialAmplitudes();

        /// outputs the kinematic scheme of inner loop integration momenta
        C_Kinematics_1loop IntMomenta(t_cmplx x, t_cmplx y, t_cmplx z);

        /// precalculates and stores getWaveFunctions for loop integration
        void PreCalculation();

        /// outputs getWaveFunctions
        std::vector<t_cmplxDirac> getWaveFunctions(C_Kinematics_1loop &Momenta);

        /// outputs FullWaveFunction
        t_cmplxDirac getFullWaveFunction(std::vector<t_cmplxDirac> &WaveFunctions, t_cmplxArray1D &U_amp);

        /// disentangles the amplutudes in case the chosen Dirac basis is not orthogonal
        virtual t_cmplxMatrix disentangleAmps(t_cmplxMatrix *pre_result) {
            t_cmplxMatrix dummy;
            std::cout << "Error - virtual call" << std::endl;
            assert(false);
            return dummy;
        };

        /// sets the Dirac basis according to meson quantum numbers
        virtual void setDiracStructures(t_cmplxVector _k, t_cmplxVector _P, std::vector<t_cmplxDirac> &DiracStructure) {
            std::cout << "Error - virtual call" << std::endl;
            assert(false);
        };

        /// dress the partons
        void preparePropagators();

        /// sets the Propagators(Partons)
        void setPropagators(t_cmplxVector &K_plus, t_cmplxVector &K_minus,
                            t_cmplxDirac &S_p, t_cmplxDirac &S_m);

        /// outputs the sum of product all Chebychev order amplitude decompositions in AmplitudesMatrix with corresponding Chebychev polinoms
        t_cmplxArray1D BSEVertex(t_cmplx z);

        /// Loop integrand
        /// "vector" denote the fact that it outputs the vector of amplitudes projected out from FullWaveFunction
        t_cmplxMatrix Integrand_vector(t_cmplxArray1D args);

        /// trace the FullWaveFunction and threadloc_Projectors with scattering Kernel
        t_cmplxMatrix traceWithKernel_Vector(t_cmplxDirac &FullWaveFunction, C_Kinematics_1loop &momenta);

        /// Outputs the left hands side of the BSE, the result of integration, the amplitudes with angle dependance
        t_cmplxMatrix outputBSA(t_cmplx _p, t_cmplx _P, int proj,
                                std::function<t_cmplxMatrix(t_cmplxArray1D)> bound_member_fn);

        /// Projects the amplitudes without angle dependence for a given Chebychev order
        t_cmplxMatrix projectBSA(int cheb_order, t_cmplxMatrix &BSA_with_angles);

        /// saves the resulting amplitudes, the right hand side of BSE
        void setBufferIn(int i, t_cmplxMatrix &BSA_with_angles);

        /// loads the resulting amplitudes, the right hand side of BSE
        void setBufferOut(int i, t_cmplxMatrix &BSA_with_angles);

        /// does one full iteration of Power Method to solve BSE
        void doOneIterationPM(t_cmplx P);

        t_cmplx solveBSE(t_cmplx P, int steps);

        t_cmplx calcLambdaEV();

        /// normalizes BSA
        /// since the eigenvectors defined up to arbitrary multiplicative factor
        /// it is wise to normalize them in order to prevent precision overflow
        /// here it is just normalized to the first element
        void normBSA();

        void showAmplitudesMatrix();

    public:
        /// switcher for homogenious or inhomogenious equations (on-shell or off-shell)
        bool flag_off_shell;

        /// switcher for the detailed description of the amplitudes
        bool flag_amp_desciption;

        /// similar to the quarks the QCD bound states are also dressed by surrounding quark-gluon medium
        void dressBSE(t_cmplxVector P) {
            solveBSE(P[3], 10); // 4th component (rest frame)
        }

        /// checkSum for tests
        double checkSum_PowerMethod();

        virtual double checkSum_EVMatrixNorm() {
            std::cout << "Base class call\n";
            return 0;
        }

        /// exports the obtained BSAs for a further use
        void setBSAonPath(t_cmplxArray2D &AmplitudePath, t_cmplxArray1D &Path, t_cmplx P);

        /// links to externally provided Kernel
        void linkToKernel(Kernels::C_AbstractKernel *_Kernel) {
            Kernel = _Kernel;
        }

        /// links to externally provided propagators(partons)
        void linkToPartons(std::vector<Propagators::C_Propagator *> _partons) {
            Parton_P = _partons[0];
            Parton_M = _partons[1];
        }
    };

}
#endif //DSEPP_BSE_TWOBODY_H
