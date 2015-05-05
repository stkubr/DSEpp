

#pragma once

#include "../NumLibs/LinearInterpolation.h"
#include "Propagator.hpp"
#include <string>

namespace Propagators {
/**
 * Enumeration of all available types of gluon sets
 * every ID denotes a specific parameters for Maris-Tandy model
 * except Arbitrary_Gluon_ID - in this case the gluon dressing function will be loaded from file
*/
	enum Gluon_ID {
		RL_MT_Light_ID = 0,
		RL_MT_Heavy_ID,
		PS_Light_ID,
		RL_MT_Heavy_DD_ID,
		Arbitrary_Gluon_ID,
		Test_Gluon_ID,
		Gluon_ID_End
	};

/**
 * \brief The Gluon model
 * Although the gluon itself can posses its own DSE here we use an approximation:
 * the effective gluon dressing function given by Maris-Tandy model. However any another
 * model can be leaded from file.
*/
	class C_Gluon : public C_Propagator {
	private:
		double D, w2;
		int LogTail;

		t_cmplx (C_Gluon::*Gluon_ref)(t_cmplx);

		std::string GluonParamPath;
		Interpolation::C_Linear<t_cmplx, t_cmplx> *FuncToInterpolate;

		/// Constructor
		C_Gluon(std::string &_GluonParamPath);

		C_Gluon();

		/// Initialization
		void setToInitialState();

		void setGluonDefaultParameters();

		/// Load interpolation points for Gluon
		void setInterpolatorPoints(std::string &_InterpolationPointsPath);

		/// Read parameters from file
		void readParameters();

		/// Gluon check on the real line
		void checkGluon();

		/// Maris-Tandy gluon model
		t_cmplx GluonByMarisTandy(t_cmplx k);

		/// Model given by interpolation over provided points set
		t_cmplx GluonByInterpolation(t_cmplx k);


	public:
		/// Parametrized Factory Method function
		static C_Gluon *createGluon(Gluon_ID id);

		static C_Gluon *createGluon(Gluon_ID id, std::string &_InterpolationPointsPath);

		void dressPropagator() {
			// So far the Gluon is given by function, so it is considered to be already dressed.
		}

		/// Get value of Gluon at /f$ k^2 /f$
		t_cmplxArray1D PropagatorAtPoint(t_cmplx k);

		/// Set manually Maris-Tandy gluon model parameters (redefined by Lambda and Etta)
		void setMTParams(double Lambda, double Etta);
	};

}