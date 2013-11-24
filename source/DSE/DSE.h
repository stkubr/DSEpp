#pragma once 

enum Quark_ID { Up_ID=0, Down_ID, Strange_ID, Charm_ID, Quark_ID_End };

class C_Quark_parameters{
	public:
	int num_prop_steps,num_angle_steps;
	double m0,mu,LimUk,LimDk,M2_contour,EffectiveCutoff,HeavyLight,Accuracy;
	bool ReCalcProp;
	
	void Print(std::ostream &__os = std::cout) const {
		 __os << "num_prop_steps" << " - " << num_prop_steps << endl <<
		 "num_angle_steps" << " - " << num_angle_steps <<endl <<
		 "m0" << " - " << m0 <<endl <<
		 "mu" << " - " << mu <<endl <<
		 "M2_contour" << " - " << M2_contour <<endl <<
		 "LimDk" << " - " << LimDk <<endl <<
		 "LimUk" << " - " << LimUk <<endl <<
		 "EffectiveCutoff" << " - " << EffectiveCutoff <<endl <<
		 "Accuracy" << " - " << Accuracy << endl <<
		 "ReCalcProp" << " - " << ReCalcProp << endl <<
		 "HeavyLight" << " - " << HeavyLight << endl; 
	}
};

class C_Propagator: public C_AbsDiagram{
	public:
	C_AbstractKernel * Kernel;
	int num_amplitudes;
	bool flag_dressed;
	
	public:
	C_Quark_parameters params;
	//virtual void info() = 0;
	virtual void DressPropagator() = 0;
	virtual void InitialState() = 0;
	virtual t_cmplxArray1D getPropAt(t_cmplx q) {t_cmplxArray1D dummy(1,0); cout << "Error virtual call" << endl; return dummy; }
	virtual void SetQuarkonPath(vector<t_cmplxMatrix> (*AmplitudePath),t_cmplxArray1D (*Path)){cout << "Error virtual call" << endl; StopLine();}
	virtual t_cmplx getDressingFactor() { t_cmplx dummy; cout << "Error virtual call" << endl; return dummy; }
	virtual void setContourApex(double M2) { t_cmplx dummy; cout << "Error virtual call" << endl; StopLine(); }
	virtual void ExportForKernel( t_cmplxArray2D * dummy) {cout << "Error virtual call" << endl; StopLine(); }
	virtual t_dArray1D GetResult(){t_dArray1D dummy; cout << "Error virtual call" << endl; return dummy;};
	
	void LinkToKernel(C_AbstractKernel * _K){
		Kernel=_K;
		//Kernel->GetNameID();
	}
};

class C_Propagator_Factory: public C_AbstractClass{
	public:
	virtual C_Propagator * Create(void *)=0;
	virtual ~C_Propagator_Factory() {}
};

ostream& operator<<(ostream &__os, const C_Quark_parameters &__params){
		__params.Print(__os);
		return __os;
	}
