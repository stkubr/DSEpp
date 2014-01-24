/*
 * Propagator.cpp
 *
 *  Created on: Jan 21, 2014
 *      Author: stkubr
 */

#include "Propagator.hpp"

t_cmplxArray1D C_Propagator::getPropAt(t_cmplx q) {
	t_cmplxArray1D dummy(1, 0);
	std::cout << "Error virtual call" << std::endl;
	return dummy;
}
void C_Propagator::SetQuarkonPath(std::vector<t_cmplxMatrix> (*AmplitudePath),
								 t_cmplxArray1D (*Path)) {
	std::cout << "Error virtual call" << std::endl;
}
t_cmplx C_Propagator::getDressingFactor() {
	t_cmplx dummy;
	std::cout << "Error virtual call" << std::endl;
	return dummy;
}
void C_Propagator::setContourApex(double M2) {
	t_cmplx dummy;
	std::cout << "Error virtual call" << std::endl;
}
void C_Propagator::ExportForKernel(t_cmplxArray2D * dummy) {
	std::cout << "Error virtual call" << std::endl;
}
t_dArray1D C_Propagator::GetTotalSum() {
	t_dArray1D dummy;
	std::cout << "Error virtual call" << std::endl;
	return dummy;
}

void C_Propagator::LinkToKernel(C_AbstractKernel * _K) {
	Kernel = _K;
}

C_Propagator::~C_Propagator() {
}
