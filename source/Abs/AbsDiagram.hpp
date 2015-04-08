#pragma once

#include "../types.h"

class C_AbsDiagram{
	protected:
	string name;
	int ID;
	DiracGamma Z,Y;
	DiracGamma5 _Y5;
	DiracSigma SIG;
	MetricTensor g;
	t_cmplxDirac I,Y5;
	double pi;
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
};
