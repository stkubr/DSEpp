#pragma once

#include "../NumLibs/Linear_interpolation.hpp"
#include "Propagator.hpp"
#include <string>

enum Gluon_ID {RL_MT_Light_ID=0, RL_MT_Heavy_ID, PS_Light_ID, RL_MT_Heavy_DD_ID, Arbitrary_Gluon_ID, Test_Gluon_ID, Gluon_ID_End};

class C_Gluon: public C_Propagator {
private:
	double D, w2;
	int LogTail;
	t_cmplx (C_Gluon::*Gluon_ref)(t_cmplx);
	std::string GluonParamPath;
	Interpolation::Linear<t_cmplx, t_cmplx> * FuncToInterpolate;

	// Constructor
	C_Gluon(std::string& _GluonParamPath);
	C_Gluon();

	// Initialization
	void setToInitialState();

	void setGluonDefaultParameters();

	// Load interpolation points for Gluon
	void setInterpolatorPoints(std::string& _InterpolationPointsPath);

	// Read parameters from file
	void ReadParameters();

	// Gluon check on the real line
	void GluonCheck();


public:
	// Parametrized Factory Method function
	static C_Gluon* createGluon(Gluon_ID id);
	static C_Gluon* createGluon(Gluon_ID id, std::string& _InterpolationPointsPath);

	void dressPropagator() {
		// So far the Gluon is given by function, so it is considered to be already dressed.
	}

	// Get value of Gluon at k
	t_cmplxArray1D PropagatorAtPoint(t_cmplx k);

	// Maris-Tandy gluon model
	t_cmplx GluonMT(t_cmplx k);
	
	// Set manually Maris-Tandy gluon model parameters (redefined by Lambda and Etta)
	void setMTParams(double Lambda, double Etta);

	// Model given by interpolation over provided points set
	t_cmplx GluonByInterpolation(t_cmplx k);

};


