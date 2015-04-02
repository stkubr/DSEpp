#pragma once

#include "Propagator.hpp"
#include "../NumLibs/Geometry/ParabolaContour.hpp"
#include "../NumLibs/Integrator.h"
#include "../NumLibs/Support_functions.h"
#include "../NumLibs/OneLoopIntegrator.hpp"
#include "Quark_parameters.hpp"

//enum Quark_ID { Up_ID=0, Down_ID, Strange_ID, Charm_ID, Test_ID, Quark_ID_End };

class C_Quark: public C_Propagator, public C_OneLoopIntegrator<t_cmplxMatrix, double, t_cmplxArray1D> {

protected:

	C_AbstractKernel * Kernel;
	const char * SavePropPath;
	int num_amplitudes;
	bool flag_dressed;
	C_Quark_parameters params;

    std::function<t_cmplx(t_cmplx, t_cmplx)> CauchyIntegratonWeight_lambda;

	C_Integrator_Line<t_cmplxMatrix, C_Quark, double> * Integrator_momentum_short;
	C_Integrator_Path<t_cmplx, t_cmplxArray2D, t_cmplx> * Integrator_cauchy;
	t_dArray1D zz_rad, w_rad, zz_line, w_line, zz_angle, w_angle, z_circus, w_circus;
	double kinematicFactor, B_renorm, B_mu, A_renorm, Z2;
	bool flag_renormalization;
	vector<int> threadloc_p_momenta_inx, threadloc_integr_inx;


	// Constructor
	C_Quark();

	//t_cmplx getTensorExpression(t_cmplxVector& p);

	void LinkToKernel(C_AbstractKernel * _K);

	void ReadParameters(string & _ParamPath);

	t_cmplx getDressingFactor();

	void setContourApex(double M2);

	// Set parameters to initial values
	void InitialState();

	// Resize all storages (internal and external), also create side objects like (Integrators, Kernels and etc.)
	void InitializateIntegrators();

	void ResizeMemory();

	// Copier of "this" (used in parallel sections)
	virtual C_Quark * MakeCopy();

	// Set Cauchy contour
	void setContour();

	// Set Cauchy contour
	void setGrid();

	// Initial guesses for A and B
	void setInitialAandB();

	// Evaluate Cauchy integral on contour, at certain point
	t_cmplxArray1D PropagatorOnPoint(t_cmplx coordin);

	// Evaluate Cauchy integral on contour, obtain Propagator on grid
	void calcPropOnGrid();

	// Evaluate DSE integral on grid, obtain Propagator on contour
	void calcPropOnContour();

	// Set k and p vectors for the Numerical Integrand
	void setKinematic(t_cmplxVector& k, t_cmplxVector& p, t_cmplx x, t_cmplx y, t_cmplx z);

	// Numerical form of the Integrand (avaible for general Kernel)
	t_cmplxMatrix Integrand_numerical(t_cmplxArray1D values);

	// Calculation A,B,M,Sigma_V,Sigma_S at point in contour
	t_cmplxArray1D getPropAt(t_cmplx q);

	// Calculate consequently Grid and Contour until converge
	void calcPropagator();

	// Apply renormalization
	void renormalizeProp();

	// Check convergence
	double checkConvergence(double previous_checksum);

	// Initialization (Dressing) of the Propagator
	void DressPropagator();

	// Draw Propagator at real line
	void drawOnRealAxis(int s);

	// Write Propagator to file
	void savePropagator();

	// Read Propagator from file
	void loadPropagator();

	// Export Propagator to file
	// (exports all what is needed to perform Cauchy integration outside of this library: contour, weights, etc.)
	void exportPropagator();

	// Saves Quark's A and B function on provided "Path" in provided "AmplutudeStorage"
	void setPropagatorOnPath(std::vector<t_cmplxMatrix> & Amplitudes, t_cmplxArray1D & Path);

	// Gets sum A and B at 100 points. Used for unit tests.
	double checkSum();


public:
	C_DedicMem_Quark * Memory;

};

