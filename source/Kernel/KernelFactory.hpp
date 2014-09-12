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

C_AbstractKernel* C_AbstractKernel::createKernel(Kernel_ID id){
	C_AbstractKernel * p;
    switch (id)
    {
        case RL_ID:
            p = new C_Kernel_RL();
            p->Kernel_type_ID=(id);
            break;
        case RL_PS_ID:
            p = new C_Kernel_RL_PS();
            p->Kernel_type_ID=(id);
            break;
        default:
            assert(false);
    }
    return p;
}

class C_Kernel_Factory{
	public:
	C_AbstractKernel* Create(Kernel_ID _id) {
		return C_AbstractKernel::createKernel( _id );
    }
	~C_Kernel_Factory() {}
};

//C_Kernel_Factory * KernelFactory = new C_Kernel_Factory;

#endif /* KERNELFACTORY_HPP_ */
