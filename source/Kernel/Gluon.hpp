#pragma once

class C_Gluon: protected C_AbsDiagram
{
	private:
	
	double D,w2;
	int LogTail,GluonType;
	bool Init_flag;
	t_cmplx (C_Gluon::*Gluon_ref)(t_cmplx*);
	const char * GluonParamPath;
	Interpolation::Linear<t_cmplx,t_cmplx> * FuncToInterpolate;
	
// Constructor
	C_Gluon(const char * __GluonParamPath){
		SetNameID("Gluon", 1);
		Init_flag=false;
		GluonParamPath=__GluonParamPath;
		Initialization();
		GluonCheck();
	}
	
	public:
// Singleton creation
	static C_Gluon * getInstance(const char * GluonParamPath) {
        return new C_Gluon(GluonParamPath);
    }
	
	
// Initialization
	void Initialization(){
		if(!Init_flag){ 
			ReadParameters();
			switch(GluonType){
				case 1:{
					Gluon_ref=&C_Gluon::GluonMT;
				} break;
				case 2:{
					//cout << "Special case  " << GluonType << endl;
					SetInterpolator_for_ChristianGluon();
					Gluon_ref=&C_Gluon::GluonFischer;
				} break;
			}
		}
		
	}
	
	void SetInterpolator_for_ChristianGluon(){
		
		t_cmplxArray2D GlounTempStorage(2);
		ifstream GluonInterpStream;
			GluonInterpStream.open("../Data_files/FWC_Nov_1.dat"); 
			int i=0;
			t_cmplx p2, dressing;
			GluonInterpStream >> p2 >> dressing;
			while(GluonInterpStream.good())
			{
				GlounTempStorage[0].push_back(p2);
				GlounTempStorage[1].push_back(dressing);
				GluonInterpStream >> p2 >> dressing;
				//cout << GlounTempStorage[0][i] << GlounTempStorage[1][i] << endl;
				i++;
			}
			//else cout << "Cant open file!" << endl;
			GluonInterpStream.close();
		
		
		FuncToInterpolate = new Interpolation::Linear<t_cmplx,t_cmplx>(GlounTempStorage[0].size(), &GlounTempStorage[0], &GlounTempStorage[1]);
	}
	
	
// Read parameters from file
	void ReadParameters(){
		string line;
		//ifstream ParamList("../Parameters_files/GluonParamList.txt");
		ifstream ParamList(GluonParamPath);
		if (ParamList.is_open())
		{
			while ( ParamList.good() )
			{
				ParamList >> line >> D;
				ParamList >> line >> w2;
				ParamList >> line >> LogTail;
				ParamList >> line >> GluonType;
			}
			cout << "D -" <<"  "<< D <<"  "<< "w2 -" <<"  "<< w2 <<"  "<< "LogTail -" <<"  "<< LogTail <<"  "<< "GluonType -" <<"  "<< GluonType << endl;
		}
		else cout << "Cant open file!" << endl;
	}
	
// Get value of Gluon at k
	t_cmplx GetGluonAt(t_cmplx * k){
		return (this->*Gluon_ref)(k);
	}
	
// Kernel check on real line
	void GluonCheck(){
		double Pu,Pd,x,scale,dp;
		scale = 100;
		Pd=0.01;
		Pu=10.0;
		dp=pow(10,(log10(Pu/Pd)/scale));
		
		x=Pd;
		ofstream Gluon;
		Gluon.open ("../Data_files/Gluon.dat");
		
		for (int i = 0; i <= scale; i++)
		{		
			t_cmplx point = x*x; 
			Gluon << real(point) <<"  "<< real(GetGluonAt(&point)) << endl;
			
			x*=dp;
			//cout << "Computing Gluon of type ...  " << GluonType << "  " << 100.0/(scale+1)*(i+1) << "%\n";
		}
		
		Gluon.close();
	}

// Maris-Tandy gluon model	
	t_cmplx GluonMT(t_cmplx * k){
		const double gamma_m=12.0/(33.0 - 2.0*4.0);
		const double m_t=0.5;
		const double tau=7.389056 - 1.0;
		const double LambdaQCD=0.234;
		//const double D_w=0.372;
		const double D_w=D;
		const double pi=3.14159265358979;
		return 4.0*pi*pi*D_w*(*k)*exp(-1.0*(*k)/w2)/w2/w2/w2  + LogTail*4.0*pi*pi*(gamma_m*(1.0-1.0*exp(-(*k)/4.0/m_t/m_t))/(*k))/(0.5*log(tau+(1.0+(*k)/LambdaQCD/LambdaQCD)*(1.0+(*k)/LambdaQCD/LambdaQCD)));
	}
	
	t_cmplx	GluonFischer(t_cmplx * k){
		return 4.0*pi*FuncToInterpolate->getValue((*k)/5.0)/(*k)/1.0;
	}
	
};

//C_Gluon * C_Gluon::p_instance=0;
