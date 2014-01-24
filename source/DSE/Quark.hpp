#pragma once

#include "Propagator.hpp"
#include "../NumLibs/Geometry/ParabolaContour.hpp"
#include "../NumLibs/Integrator.h"
#include "../Kernel/Gluon.hpp"
#include "../NumLibs/Support_functions.h"

class C_Quark: public C_Propagator {

protected:
	const char * SavePropPath;
	C_Gluon * Gluon;

	t_cmplxMatrix (C_Quark::*integrand)(t_cmplxArray1D);

	C_Integrator_Line<t_cmplxMatrix, C_Quark, double> * Integ_radial_leg;
	C_Integrator_Line<t_cmplxMatrix, C_Quark, double> * Integ_radial_short_leg;
	C_Integrator_Line<t_cmplxMatrix, C_Quark, double> * Integ_angle_cheb;
	C_Integrator_Cauchy<t_cmplxArray1D, t_cmplxArray3D, t_cmplx> * Integ_cauchy_long;
	int index_p, grid1_num;
	t_cmplx x;
	t_dArray1D zz_rad, w_rad, zz_angle, w_angle, z_circus, w_circus;
	t_cmplxArray1D integrand_args;
	double B_renorm, B_mu, A_renorm, Z2, check_res, eps;
	bool flag_normalized;

	t_cmplxVector k, p;

	C_Quark();

	void ReadParameters(ifstream & _ParamList);

	t_cmplx getDressingFactor();

	void setContourApex(double M2);

// Set parameters to initial values
//----------------------------------------------------------------------
	void InitialState();

// Resize all storages (internal and external), also create side objects like (Integrators, Kernels and etc.)
//----------------------------------------------------------------------
	void InitializateIntegrators();

	void ResizeMemory();

// Copier of "this" (used in parallel sections)
//----------------------------------------------------------------------	
	virtual C_Quark * MakeCopy();

// Set Cauchy contour
//----------------------------------------------------------------------
	void setContourAndGrid();

// Initial guesses for A and B
//----------------------------------------------------------------------	
	t_cmplx InitStepA(t_cmplx z);
	t_cmplx InitStepB(t_cmplx z);

// Multidimensional integration on complex plane
//----------------------------------------------------------------------
	t_cmplxMatrix Multi_INT_cx(t_cmplxMatrix (C_Quark::*func_to_int)(t_cmplxArray1D));
	t_cmplxMatrix f1(double y);
	t_cmplxMatrix f2(double z);

// Evaluate Cauchy integral on contour, at certain point
//----------------------------------------------------------------------
	t_cmplxArray1D getCauchyAt_embedded(t_cmplx coordin);

// Evaluate Cauchy integral on contour, obtain Propogator on grid
//----------------------------------------------------------------------
	void CalcPropGrid();

// Evaluate DSE integral on grid, obtain Propagator on contour
//----------------------------------------------------------------------
	void CalcPropCont();

// Analytic form of the Integrand (available only for RL or Pion Contribution)
//----------------------------------------------------------------------
	t_cmplxMatrix Integrand_analitic(t_cmplxArray1D values);

// Set k and p vectors for the Numerical Integrand
//----------------------------------------------------------------------	
	void setKinematic(t_cmplx x, t_cmplx y, t_cmplx z);

// Numerical form of the Integrand (avaible for general Kernel)
//----------------------------------------------------------------------
	t_cmplxMatrix Integrand_numerical(t_cmplxArray1D values);

// Calculation A,B,M,Sigma_V,Sigma_S at point in contour
//----------------------------------------------------------------------
	t_cmplxArray1D getPropAt(t_cmplx q);

// Calculate consequently Grid and Contour until converge
//----------------------------------------------------------------------	
	void PropSetAndCheck();

// Check convergence
//----------------------------------------------------------------------
	void PropCheck(int s);

// Initialization (Dressing) of the Propagator
//----------------------------------------------------------------------	
	void DressPropagator();

// Drow Propagator at real line
//----------------------------------------------------------------------
	void write_Prop_re(int s);

// Write Propagator to file
//----------------------------------------------------------------------
	void SavePropCountour();

// Read Propagator from file
//----------------------------------------------------------------------
	void LoadPropCountour();

	void ExportPropagator();

	void SetQuarkonPath(std::vector<t_cmplxMatrix> (*Amplitudes),
			t_cmplxArray1D (*Path));

	t_dArray1D GetTotalSum();

public:
	static C_Quark* createQuark(Quark_ID id);
	C_DedicMem_Quark * Memory;

};
