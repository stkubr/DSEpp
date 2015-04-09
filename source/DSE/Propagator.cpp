/*
 * Propagator.cpp
 *
 *      Author: stkubr
 */

#include "Propagator.hpp"

t_cmplxArray1D C_Propagator::PropagatorAtPoint(t_cmplx q) {
	t_cmplxArray1D dummy(1, 0);
	std::cout << "Pure virtual call:C_Propagator::PropagatorAtPoint" << std::endl;
	assert(false);
	return dummy;
}

void C_Propagator::SetPropagatorOnPath(std::vector<t_cmplxMatrix> & AmplitudesOnPath,
								 t_cmplxArray1D & Path) {
	std::cout << "Pure virtual call:C_Propagator::SetQuarkonPath" << std::endl;
	assert(false);
}

t_cmplx C_Propagator::DressingFactor() {
	t_cmplx dummy;
	std::cout << "Pure virtual call:C_Propagator::DressingFactor" << std::endl;
	assert(false);
	return dummy;
}

void C_Propagator::setContourApex(double M2) {
	std::cout << "Pure virtual call:C_Propagator::setContourApex" << std::endl;
	assert(false);
}

double C_Propagator::checkSum() {
	double dummy;
	std::cout << "Pure virtual call:C_Propagator::checkSum" << std::endl;
	assert(false);
	return dummy;
}

void C_Propagator::linkToKernel(C_AbstractKernel *_K) {
	std::cout << "Pure virtual call:C_Propagator::linkToKernel" << std::endl;
	assert(false);
}
