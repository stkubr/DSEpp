/*
 * RainbowLadderKernel.hpp
 *
 *  Created on: Jan 28, 2014
 *      Author: stkubr
 */

#ifndef RAINBOWLADDERKERNEL_HPP_
#define RAINBOWLADDERKERNEL_HPP_

#include "AbstractKernel.hpp"

class C_Kernel_RL: public C_AbstractKernel{
	public:

	C_Kernel_RL();

    void info() { std::cout << "Kernel RainbowLadder" << std::endl; }

	void SetMediators(t_cmplxVector& k,
					  t_cmplxVector& p,
					  t_cmplxVector& P,
					  std::vector<t_cmplxTensor>& Mediators);

	t_cmplx getElementKmatrix(int t, int s, int r, int u,
							  std::vector<t_cmplxTensor>& Mediators);

	void SpecifyGluon(Gluon_ID gluon_id);
};



#endif /* RAINBOWLADDERKERNEL_HPP_ */
