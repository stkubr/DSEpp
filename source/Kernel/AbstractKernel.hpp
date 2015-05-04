/* \file AbstractKernel.hpp
 *      Author: stkubr
 */

#ifndef ABSTRACTKERNEL_HPP_
#define ABSTRACTKERNEL_HPP_

namespace Propagators{
	class C_Propagator;
}

#include "../Abs/AbsDiagram.hpp"
#include "../DSE/Propagator.hpp"
#include "../DedicMem/MemoryFactories.hpp"

namespace Kernels {
/**
 * Enumeration of all available types of two-body scattering kernels
 * RL - Rainbow-ladder single gluon exchange
 * RL_PS - Rainbow-ladder single gluon exchange + PseudoScalar exchange (aka pion cloud)
 */
	enum Kernel_ID {
		RL_ID = 0, RL_PS_ID, Kernel_ID_End
	};

/**
 * Enumeration of all available types PseudoScalar exchange
 * Pion_exchange - the pion cloud effect
 * Etta_exchange - /f$ \eta_c /f$ exchange for charmonium
*/
	enum PS_type_ID {
		Pion_exchange_ID = 0, Etta_exchange_ID, PS_type_ID_End
	};

/**
 * \brief The interface class for all possible scattering kernels
*/
	class C_AbstractKernel : public C_AbsDiagram {
	public:
		virtual void info() = 0;

		virtual void setMesonExchangeMass(t_cmplx _M) { }

		virtual void setExchangeID(PS_type_ID exchange_id) { }

		/// Sets the the propagators used inside the kernel as mediator
		/// for example: in case of single gluon exchange kernel the only mediating particle inside the kernel is gluon
		/// and Z_2 renormalization factor is given by the quark to be dressed by this kernel,
		/// so that the calculations are selfconsistent.
		virtual void setPropagators(std::vector<Propagators::C_Propagator *> _Propagators) = 0;

		/// In case of non-trivial kernel (like with pion exchange) the BSE and the DSE might should be dressed by kernel
		/// differently. This function handles that.
		virtual void setConvolutionType(int type) { }

		/// Take a trace of Tr(Projector * Kernel * WaveFunction) at provided momenta
		virtual t_cmplx TraceKernelWithoutStoring(t_cmplxDirac &Projector,
												  t_cmplxDirac &WaveFunc,
												  t_cmplxVector &k,
												  t_cmplxVector &p,
												  t_cmplxVector &P, bool flag_reset_kernel) = 0;
	};
}
#endif /* ABSTRACTKERNEL_HPP_ */
