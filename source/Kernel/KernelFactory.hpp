/*
 * KernelFactory.hpp
 *
 *  Created on: Jan 28, 2014
 *      Author: stkubr
 */

#ifndef KERNELFACTORY_HPP_
#define KERNELFACTORY_HPP_

#include "RainbowLadderKernel.hpp"
#include "RLandPseudoScalar.hpp"

class C_Kernel_Factory{
private:
    C_Kernel_Factory(){}

public:
    static C_Kernel_Factory& instance(){
        static C_Kernel_Factory ins;
        return ins;
    }

	C_AbstractKernel* Create(Kernel_ID _id) {
        C_AbstractKernel * p;
        switch (_id)
        {
            case RL_ID:
                p = new C_Kernel_RL();
                break;
            case RL_PS_ID:
                p = new C_Kernel_RL_PS();
                break;
            default:
                assert(false);
        }
        return p;
    }
};

//C_Kernel_Factory * KernelFactory = new C_Kernel_Factory;

#endif /* KERNELFACTORY_HPP_ */
