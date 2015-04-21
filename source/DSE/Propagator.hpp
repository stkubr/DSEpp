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

	/**
 	* Dress Propagator according to defined in derived class DSE scheme
 	* this is the function where all !!SCIENCE!! of DSE happens
 	*/
	virtual void dressPropagator()=0;

	/**
 	* Reset parameters to initial values
 	*/
	virtual void setToInitialState()=0;

	/**
 	* Return the value of all possible for this kind of propagator dressing functions at point
 	*/
	virtual t_cmplxArray1D PropagatorAtPoint(t_cmplx q);

	// TODO virtual t_cmplxDirac getTensorExpression(t_cmplxVector& p);

	/**
 	* Save propagator dressing functions, evaluated on provided "Path", in provided "AmplitudeStorage"
 	*/
	virtual void setPropagatorOnPath(t_cmplxArray2D &AmplitudesOnPath, t_cmplxArray1D &Path);

	virtual t_cmplx DressingFactor();

	virtual void setContourApex(double M2);

	virtual double checkSum();

	virtual void linkToKernel(C_AbstractKernel * _Kernel);

	virtual ~C_Propagator() {}
};

#endif /* PROPAGATOR_HPP_ */
