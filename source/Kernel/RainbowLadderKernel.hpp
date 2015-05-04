/*
 * RainbowLadderKernel.hpp
 *
 *  Created on: Jan 28, 2014
 *      Author: stkubr
 */

#ifndef RAINBOWLADDERKERNEL_HPP_
#define RAINBOWLADDERKERNEL_HPP_

#include "TwoQuarkKernel.hpp"

namespace Kernels {

/**
 * \brief The Rainbow-Ladder scattering Kernel
 *
 * Pure single gluon exchange. setMediators() sets the gluon polarization tensor with effective gluon dressing function.
*/
    class C_Kernel_RL : public C_TwoQuarkKernel {
        friend class C_Kernel_Factory;

    protected:
        t_cmplx Z2;

        C_Kernel_RL() {
            SetNameID("Kernel RainbowLadder", 1);
            Memory->VertexDressings.resize(2, t_cmplxArray3D(1));
            Z2 = 1.0;
        }

        /// Sets gluon propagating inside of Kernel
        void setMediators(t_cmplxVector &k,
                          t_cmplxVector &p,
                          t_cmplxVector &P,
                          std::vector<t_cmplxTensor> &Mediators);

        /// Set a (t,s,r,u) element of the K_matrix for a given class
        t_cmplx ElementKmatrix(int t, int s, int r, int u,
                               std::vector<t_cmplxTensor> &Mediators);

    public:
        void info() { std::cout << "Kernel RainbowLadder" << std::endl; }
    };
}
#endif /* RAINBOWLADDERKERNEL_HPP_ */
