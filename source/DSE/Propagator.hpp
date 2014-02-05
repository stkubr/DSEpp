/*
 * Propagator.hpp
 *
 *  Created on: Jan 21, 2014
 *      Author: stkubr
 */

#ifndef PROPAGATOR_HPP_
#define PROPAGATOR_HPP_

#include "Quark_parameters.hpp"
#include "../types.h"
#include "../Abs/AbsDiagram.h"
#include "../Kernel/AbstractKernel.hpp"


enum Quark_ID { Up_ID=0, Down_ID, Strange_ID, Charm_ID, Quark_ID_End };

class C_Propagator: public C_AbsDiagram{
	public:
	C_AbstractKernel * Kernel;
	int num_amplitudes;
	bool flag_dressed;

	public:
	C_Quark_parameters params;
	//virtual void info() = 0;
	virtual void DressPropagator()=0;
	virtual void InitialState()=0;
	virtual t_cmplxArray1D getPropAt(t_cmplx q);
	//virtual t_cmplxDirac getTensorExpression(t_cmplxVector& p);
	virtual void SetQuarkonPath(std::vector<t_cmplxMatrix> (*AmplitudePath),t_cmplxArray1D (*Path));
	virtual t_cmplx getDressingFactor();
	virtual void setContourApex(double M2);
	virtual void ExportForKernel( t_cmplxArray2D * dummy);
	virtual t_dArray1D GetTotalSum();
	virtual ~C_Propagator();

	void LinkToKernel(C_AbstractKernel * _K);
};

class C_Propagator_Factory{
	public:
	virtual C_Propagator * Create(int)=0;
	virtual ~C_Propagator_Factory() {}
};

#endif /* PROPAGATOR_HPP_ */
