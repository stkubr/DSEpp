/*
 * Propagator.hpp
 *
 *  Created on: Jan 21, 2014
 *      Author: stkubr
 */

#ifndef PROPAGATOR_HPP_
#define PROPAGATOR_HPP_

namespace Kernels{
	class C_AbstractKernel;
}

#include "../types.h"
#include "../Abs/AbsDiagram.hpp"
#include "../Kernel/AbstractKernel.hpp"

namespace Propagators {

/**
 * \brief The interface class for all propagators: quarks, gluons or ghosts
 *
 * (yep, fancy name for a particle).
*/
	class C_Propagator : public C_AbsDiagram {
	public:
		/**
 	* Dress Propagator according to defined in derived class DSE scheme
 	* this is the function where all !!SCIENCE!! of DSE happens
 	*/
		virtual void dressPropagator() = 0;

		/**
 	* Reset parameters to initial values
 	*/
		virtual void setToInitialState() = 0;

		/**
 	* Return the value of all possible for this kind of propagator dressing functions at point
 	*/
		virtual t_cmplxArray1D PropagatorAtPoint(t_cmplx q);

		/**
 	* Save propagator dressing functions, evaluated on provided "Path", in provided "AmplitudeStorage"
 	*/
		virtual void setPropagatorOnPath(t_cmplxArray2D &AmplitudesOnPath, t_cmplxArray1D &Path);

		/// Outputs the /f$ Z_2(\mu^2) /f$ renormalization factor
		virtual t_cmplx Z2Factor();

		/// Sets the /f$ -\frac{M^2}{4} /f$ for the contour
		virtual void setContourApex(double M2);

		/// Gets sum A and B at 100 points. Used for unit tests.
		virtual double checkSum();

		/// Sets the class member Kernel pointer to one provided externally
		virtual void linkToKernel(Kernels::C_AbstractKernel *_Kernel);

		virtual ~C_Propagator() { }
	};

}
#endif /* PROPAGATOR_HPP_ */
