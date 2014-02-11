#pragma once

#include "../NumLibs/Linear_interpolation.h"
#include "Propagator.hpp"

enum Gluon_ID {RL_MT_Light_ID=0, RL_MT_Heavy_ID, PS_Light_ID, RL_MT_Heavy_DD_ID, Gluon_ID_End};

class C_Gluon: public C_Propagator
{
	private:
	
	double D,w2;
	int LogTail,GluonType;
	bool Init_flag;
	t_cmplx (C_Gluon::*Gluon_ref)(t_cmplx);
	const char * GluonParamPath;
	Interpolation::Linear<t_cmplx,t_cmplx> * FuncToInterpolate;
	
// Constructor
	C_Gluon(const char * __GluonParamPath);
	
	public:

	void DressPropagator(){}

	static C_Gluon* createGluon( Gluon_ID id );

// Initialization
	void InitialState();
	
	void SetInterpolator_for_ChristianGluon();
	
	
// Read parameters from file
	void ReadParameters();
	
// Get value of Gluon at k
	t_cmplxArray1D getPropAt(t_cmplx k);
	
// Kernel check on real line
	void GluonCheck();

// Maris-Tandy gluon model	
	t_cmplx GluonMT(t_cmplx k);
	
	t_cmplx	GluonFischer(t_cmplx k);

};

class C_Gluon_Factory: public C_Propagator_Factory{
	public:
	C_Propagator* Create(int _id) {
		Gluon_ID id=(Gluon_ID)_id;
		return C_Gluon::createGluon( id );
    }
};
