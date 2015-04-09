/**
 * Quark.hpp
 * Author: stkubr
 */

#ifndef QUARK_HPP_
#define QUARK_HPP_

#include "Propagator.hpp"
#include "../NumLibs/Geometry/ParabolaContour.hpp"
#include "../NumLibs/Integrator.hpp"
#include "../NumLibs/Integrator_Path.hpp"
#include "../NumLibs/Support_functions.h"
#include "../NumLibs/OneLoopIntegrator.hpp"
#include "Quark_parameters.hpp"
#include "../DedicMem/MemoryFactories.hpp"

class C_Quark: public C_Propagator, public C_OneLoopIntegrator<t_cmplxMatrix, double, t_cmplxArray1D> {

protected:

	/// two-body scattering Kernel
	C_AbstractKernel * Kernel;

	/// quark's Memory handler class, contains the contour and the grid at which quark DSE evaluated
	C_DedicMem_Quark * Memory;

	/// path to file to store the calculated propagator
	string SavePropPath;

	/// number of projected out amplitudes (quark has two, A and B)
	int num_amplitudes;

	/// contains all params the quark DSE need to have
	C_Quark_parameters params;

	/// defines the Cauchy integral measure
    std::function<t_cmplx(t_cmplx, t_cmplx)> CauchyIntegratonWeight_lambda;

	/// Line Integrator for the cutoff part of parabolic contour
	C_Integrator_Line<t_cmplxMatrix, double> *Integrator_momentum_cutoff;

	/// Path Integrator for the complete contour
	C_Integrator_Path<t_cmplx, t_cmplxArray2D, t_cmplx> * Integrator_cauchy;

	/// zz are the integration nodes, w are the weights
	t_dArray1D zz_radial, w_radial, zz_cutoff, w_cutoff, zz_angle, w_angle;

	/// just a momentum volume factor
	double kinematicFactor;

	/// Renormalization constants
	double B_renorm, B_mu, A_renorm, Z2;

	/// thread-local storages of indexes
	vector<int> threadloc_p_momenta_inx, threadloc_integr_inx;

	bool flag_dressed;
	bool flag_renormalization;

	// the constructor
	C_Quark();

	// the destructor
	~C_Quark();

	//t_cmplx getTensorExpression(t_cmplxVector& p);

	void linkToKernel(C_AbstractKernel *_K);

	void readQuarkParameters(string &_ParamPath);

	t_cmplx DressingFactor();

	void setContourApex(double M2);

	// Set parameters to initial values
	void setToInitialState();

	// Create Integrators and get the integration points and weights out of them
	void initializateIntegrators();

	// Resize the dedicated memory storages for the contour and the grid according to loaded "params"
	void resizeMemory();

	// Set Contour in complex plane
	void setContour();

	// Set Grid in complex plane
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
	t_cmplxArray1D PropagatorAtPoint(t_cmplx q);

	// Calculate consequently Grid and Contour until converge
	void calcPropagator();

	// Apply renormalization
	void renormalizeProp();

	// Check convergence
	double checkConvergence(double previous_checksum);

	// Dress the Propagator
	void dressPropagator();

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

};

#endif /* QUARK_HPP_ */