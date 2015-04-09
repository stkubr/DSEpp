/*
 * AbstractKernel.hpp
 *      Author: stkubr
 */

#ifndef ABSTRACTKERNEL_HPP_
#define ABSTRACTKERNEL_HPP_

/**
 * Enumeration of all available types of two-body scattering kernels
 * RL - Rainbow-ladder single gluon exchange
 * RL_PS - Rainbow-ladder single gluon exchange + PseudoScalar exchange (aka pion cloud)
*/
enum Kernel_ID {RL_ID=0, RL_PS_ID, Kernel_ID_End};

/**
 * Enumeration of all available types PseudoScalar exchange
 * Pion_exchange - the pion cloud effect
 * Etta_exchange - /f$ \eta_c /f$ exchange for charmonium
*/
enum PS_type_ID {Pion_exchange_ID=0, Etta_exchange_ID, PS_type_ID_End};

#include "../Abs/AbsDiagram.hpp"
class C_Propagator;
class C_Kernel_Factory;
#include "../DSE/Propagator.hpp"
#include "../DedicMem/MemoryFactories.hpp"

class C_AbstractKernel: public C_AbsDiagram{
public:
	virtual void info()=0;

	virtual void setMesonExchangeMass(t_cmplx _M){}

	virtual void setExchangeID(PS_type_ID exchange_id){}

	virtual void setPropagators(std::vector<C_Propagator*> __Propagators)=0;

	virtual void setConvolutionType(int type){}

	// Take a trace of Tr(Projector * Kernel * WaveFunction) at provided momenta
	virtual t_cmplx TraceKernelWithoutStoring(t_cmplxDirac& Projector,
									  t_cmplxDirac& WaveFunc,
									  t_cmplxVector& k,
									  t_cmplxVector& p,
									  t_cmplxVector& P, bool flag_reset_kernel)=0;
};

#endif /* ABSTRACTKERNEL_HPP_ */
