/*
 * AbsDiagram.hpp
 *      Author: stkubr
 */

#ifndef ABSDIAGRAM_HPP_
#define ABSDIAGRAM_HPP_

#include "../types.h"

class C_AbsDiagram{
	protected:
	string name;
	int ID;

	/// Dirac gamma matrices \f$ \gamma_\mu \f$.
	DiracGamma Z,Y;

	/// Dirac gamma matrix \f$ \gamma_5 \f$.
	DiracGamma5 _Y5;

	/// Dirac sigma matrix \f$ \sigma_{\mu\nu} \f$.
	DiracSigma SIG;

	/// Euclidean metric tensor \f$ g_{\mu\nu} \f$.
	MetricTensor g;

	/// Unit matrix (4,4) in Dirac space
	t_cmplxDirac I;

	t_cmplxDirac Y5;
	double pi;

	/// Imaginary unit
	t_cmplx ii;
	

	C_AbsDiagram(){
		ii=t_cmplx(0.0,1.0);
		pi=3.14159265358979;
		I=(1.0/4.0*(Z*Z));
		Y5=I*_Y5;
	}

	public:
	void SetNameID(string _name, int _ID){
		name=_name;
		ID=_ID;
		//std::cout << "Object -" << "  " << name << "  " << "with ID -" << "  " << ID << "  " << "has been created." << std::endl;
	}
	
	void GetNameID(){
		std::cout << "This is object -" <<"  "<< name <<"  "<< "with ID -" <<"  "<< ID  << std::endl;
	}

	virtual ~C_AbsDiagram(){}
};

#endif /* ABSDIAGRAM_HPP_ */