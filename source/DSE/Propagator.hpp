/*
 * Propagator.hpp
 *
 *  Created on: Jan 21, 2014
 *      Author: stkubr
 */

#ifndef PROPAGATOR_HPP_
#define PROPAGATOR_HPP_

#include "../types.h"
#include "../Abs/AbsDiagram.hpp"
class C_AbstractKernel;
#include "../Kernel/AbstractKernel.hpp"


class C_Propagator: public C_AbsDiagram{
	public:

	virtual void DressPropagator()=0;

	virtual void InitialState()=0;

	virtual t_cmplxArray1D getPropAt(t_cmplx q);

	// TODO virtual t_cmplxDirac getTensorExpression(t_cmplxVector& p);

	virtual void SetPropagatorOnPath(std::vector<t_cmplxMatrix> & AmplitudesOnPath, t_cmplxArray1D & Path);

	virtual t_cmplx getDressingFactor();

	virtual void setContourApex(double M2);

	virtual double checkSum();

	virtual void LinkToKernel(C_AbstractKernel * _K);

	virtual ~C_Propagator() {}
};

class C_Propagator_Factory{
	public:
	virtual C_Propagator * Create(int)=0;
	virtual ~C_Propagator_Factory() {}
};

#endif /* PROPAGATOR_HPP_ */
