#pragma once

#include "../NumLibs/Linear_interpolation.h"
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
	void InitialState();

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

	void DressPropagator() {
		// So far the Gluon is given by function, so it is considered to be already dressed.
	}

    t_cmplx getDressingFactor(){
        std::cout << "Maris-Tandy gluon cannot provide the dressing factor!" << std::endl;
        assert(false);
    }

	// Get value of Gluon at k
	t_cmplxArray1D getPropAt(t_cmplx k);

	// Maris-Tandy gluon model
	t_cmplx GluonMT(t_cmplx k);
	
	// Set manually Maris-Tandy gluon model parameters (redefined by Lambda and Etta)
	void setMTParams(double Lambda, double Etta);

	// Model given by interpolation over provided points set
	t_cmplx GluonByInterpolation(t_cmplx k);

};

class C_Gluon_Factory: public C_Propagator_Factory{
	public:
	C_Propagator* Create(int _id) {
		Gluon_ID id=(Gluon_ID)_id;
		return C_Gluon::createGluon( id );
    }

	C_Propagator* Create(int _id, std::string & _InterpolationPointsPath) {
		Gluon_ID id=(Gluon_ID)_id;
		return C_Gluon::createGluon( id, _InterpolationPointsPath );
	}
};
