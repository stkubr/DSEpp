//
// Created by stkubr on 21.04.15.
//

#ifndef DSEPP_BSE_H
#define DSEPP_BSE_H


#include <source/types.h>
#include <source/Abs/AbsDiagram.hpp>
#include <source/Kernel/AbstractKernel.hpp>

namespace BSE {

/**
 * \brief The interface class for all BSEs: two-body (mesons) or more
 *
*/
    class C_BSE : public C_AbsDiagram {
    public:
        virtual t_cmplxArray2D calcEigenStates(t_cmplxVector P, int numberOfStates) = 0;

        virtual void dressBSE(t_cmplxVector P) = 0;

        virtual void setBSAonPath(t_cmplxArray2D &AmplitudePath, t_cmplxArray1D &Path, t_cmplx P) = 0;

        virtual double checkSum_PowerMethod() = 0;

        virtual double checkSum_EVMatrixNorm() = 0;

        virtual void linkToKernel(Kernels::C_AbstractKernel *_Kernel) = 0;

        virtual void linkToPartons(std::vector<Propagators::C_Propagator *> _Partons) = 0;

        virtual ~C_BSE() { }
    };

}
#endif //DSEPP_BSE_H
