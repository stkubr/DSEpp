/*
 * DedicMem.cpp
 *
 *  Created on: Jan 24, 2014
 *      Author: stkubr
 */

#include "DedicMem.hpp"


void C_DedicMem_Quark::resizeGrid(int amps_num, int exter_num, int inter_num) {
    S_grid.resize(amps_num, t_cmplxArray2D(exter_num, t_cmplxArray1D(inter_num)));
}

void C_DedicMem_Quark::resizeContour(int amps_num, int cont_num) {
    S_cont.resize(amps_num, t_cmplxArray1D(cont_num));
}

void C_DedicMem_Quark::RemoveGrid() {
    S_grid.clear();
    std::cout << "Grid storage has been erased." << std::endl;
}

void C_DedicMem_Quark::RemoveContour() {
    S_cont.clear();
    std::cout << "Grid storage has been erased." << std::endl;
}

void C_DedicMem_Quark::info() {
    std::cout << "Dedicated Memory for Quark" << std::endl;
}


t_cmplxArray4D VertexDressings;

void C_DedicMem_Kernel::info() {
    std::cout << "Dedicated Memory for Kernel" << std::endl;
}

void C_DedicMem_Kernel::ResizeKstorage(int i) {
    K_Matrix_Storage.resize(i);
}

void C_DedicMem_Kernel::SetKmatrixAt(int i, t_cmplxMatrix2D *_K) {
    K_Matrix_Storage[i] = (*_K);
}

void C_DedicMem_Kernel::EraseKstorage() {
    K_Matrix_Storage.clear();
}

t_cmplxMatrix2D *C_DedicMem_Kernel::GetKmatrixAt(int i) {
    return &K_Matrix_Storage[i];
}


// AmpStorage Manipulation
//----------------------------------------------------------------------
void C_DedicMem_BSE::resizeAmpStorage(int amps, int points) {
    AmpStorage.resize(amps, std::vector<t_cmplxDirac>(points));
}

void C_DedicMem_BSE::setAmpStorage(int amp, int point, t_cmplxDirac *Amp) {
    AmpStorage[amp][point] = (*Amp);
}

t_cmplxDirac C_DedicMem_BSE::getAmpStorage(int amp, int point) {
    return AmpStorage[amp][point];
}

void C_DedicMem_BSE::clearAmpStorage() {
    AmpStorage.clear();
}
//----------------------------------------------------------------------


//EVMatrix Manipulation
//----------------------------------------------------------------------
void C_DedicMem_BSE::ResizeEVMatrix(int num_rad, int num_angle, int num_amplitudes, int num_chebs) {
    EVMatrix = Eigen::MatrixXcf::Ones(num_rad * num_angle * num_amplitudes * num_chebs, num_rad * num_angle * num_amplitudes * num_chebs);
}
//----------------------------------------------------------------------

//BSE_onCauchy Contour ang Grid Manipulation
//----------------------------------------------------------------------
void C_DedicMem_BSE::resizeBSEContour(int num_amplitudes, int num_points) {
    CauchyContour.resize(num_amplitudes, t_cmplxArray1D(num_points));
}

void C_DedicMem_BSE::resizeBSEGrid(int num_amplitudes, int num_ex, int num_in) {
    CauchyGrid.resize(num_amplitudes, t_cmplxArray2D(num_ex, t_cmplxArray1D(num_in)));
}

void C_DedicMem_BSE::clearCauchyGrid() {
    CauchyGrid.clear();
}

//----------------------------------------------------------------------

void C_DedicMem_BSE::info() {
    std::cout << "Dedicated Memory for BSE" << std::endl;
}


