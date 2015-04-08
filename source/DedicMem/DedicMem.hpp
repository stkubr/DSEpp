/*
 * DedicMem.hpp
 *
 *  Created on: Jan 21, 2014
 *      Author: stkubr
 */

#ifndef DEDICMEM_HPP_
#define DEDICMEM_HPP_

#include <vector>
#include <fstream>
#include "Eigen/Eigenvalues"
#include "../types.h"

class C_DedicMem_Abs{
	
	public:
	virtual void info() = 0;
	virtual ~C_DedicMem_Abs() {}
};

class C_DedicMem_Quark: public C_DedicMem_Abs {
	public:
	t_cmplxArray3D S_grid;
	t_cmplxArray2D S_cont;
	
	public:
	void resizeGrid(int amps_num, int exter_num, int inter_num);
	void resizeContour(int amps_num, int cont_num);
	
	void RemoveGrid();

	void info();
};

class C_DedicMem_Kernel: public C_DedicMem_Abs {
	
	private:
	std::vector<t_cmplxMatrix2D > K_Matrix_Storage;
	
	public:
	t_cmplxArray4D VertexDressings;
	void info();

	void ResizeKstorage(int i);
	
	void SetKmatrixAt(int i ,t_cmplxMatrix2D * _K);
	
	void EraseKstorage();
	
	t_cmplxMatrix2D * GetKmatrixAt(int i);
	
};

class C_DedicMem_BSE: public C_DedicMem_Abs {
	
	public:
	std::vector< std::vector <t_cmplxDirac> > AmpStorage;
	Eigen::MatrixXcf EVMatrix;
	t_cmplxArray2D CauchyContour;
	t_cmplxArray3D CauchyGrid;
	
	public:
	
	// AmpStorage Manipulation
//----------------------------------------------------------------------
	void resizeAmpStorage(int amps, int points);
	
	void setAmpStorage(int amp, int point, t_cmplxDirac * Amp);
	
	t_cmplxDirac getAmpStorage(int amp, int point);
	
	void clearAmpStorage();
//----------------------------------------------------------------------


	//EVMatrix Manipulation
//----------------------------------------------------------------------
	void ResizeEVMatrix(int num_rad, int num_angle, int num_amplitudes, int num_chebs);
//----------------------------------------------------------------------
	
	//BSE_onCauchy Contour ang Grid Manipulation
//----------------------------------------------------------------------
	void resizeBSEContour(int num_amplitudes, int num_points);
	
	void resizeBSEGrid(int num_amplitudes, int num_ex, int num_in);
	
	void clearCauchyGrid();
	
//----------------------------------------------------------------------
	
	void info();
};


#endif /* DEDICMEM_HPP_ */
