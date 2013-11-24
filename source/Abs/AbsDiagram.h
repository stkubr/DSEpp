#pragma once

class C_AbsDiagram: public C_AbstractClass{	
	protected:
	string name;
	int ID;
	C_DedicMem_Abs * Memory_abs;
	DiracGamma Z,Y;
	DiracGamma5 _Y5;
	DiracSigma SIG;
	MetricTensor g;
	t_cmplxDirac I,Y5;
	double pi;
	t_cmplx ii;
	
	public:
	C_AbsDiagram(){
		ii=t_cmplx(0.0,1.0);
		pi=3.14159265358979;
		I=(1.0/4.0*(Z*Z));
		Y5=I*_Y5;
	}
	
	void SetNameID(string _name, int _ID){
		name=_name;
		ID=_ID;
		//cout << "Object -" << "  " << name << "  " << "with ID -" << "  " << ID << "  " << "has been created." << endl;
	}
	
	void GetNameID(){
		cout << "This is object -" <<"  "<< name <<"  "<< "with ID -" <<"  "<< ID  << endl;
	}
};
