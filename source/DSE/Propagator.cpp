/*
 * Propagator.cpp
 *
 *      Author: stkubr
 */

#include "Propagator.hpp"

using namespace Propagators;

t_cmplxArray1D C_Propagator::PropagatorAtPoint(t_cmplx q) {
	t_cmplxArray1D dummy(1, 0);
	std::cout << "Pure virtual call:C_Propagator::PropagatorAtPoint" << std::endl;
	assert(false);
	return dummy;
}

void C_Propagator::setPropagatorOnPath(t_cmplxArray2D &AmplitudesOnPath,
									   t_cmplxArray1D &Path) {
	std::cout << "Pure virtual call:C_Propagator::SetQuarkonPath" << std::endl;
	assert(false);
}

t_cmplx C_Propagator::Z2Factor() {
	t_cmplx dummy;
	std::cout << "Pure virtual call:C_Propagator::Z2Factor" << std::endl;
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

void C_Propagator::linkToKernel(Kernels::C_AbstractKernel *_K) {
	std::cout << "Pure virtual call:C_Propagator::linkToKernel" << std::endl;
	assert(false);
}
