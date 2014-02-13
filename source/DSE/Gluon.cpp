/*
 * Gluon.cpp
 *
 *  Created on: Feb 7, 2014
 *      Author: stkubr
 */

#include "Gluon.hpp"

// Constructor
//----------------------------------------------------------------------
C_Gluon::C_Gluon(std::string& __GluonParamPath){
	SetNameID("Gluon", 1);
	GluonParamPath=__GluonParamPath;
	Gluon_ref=&C_Gluon::GluonMT;
	ReadParameters();
	GluonCheck();
}

C_Gluon::C_Gluon(){
	SetNameID("Gluon", 1);
}


// Parameterized Factory Method function
//----------------------------------------------------------------------
C_Gluon* C_Gluon::createGluon(Gluon_ID id){
	C_Gluon * p;
	std::string __GluonParamPath;
	switch (id){
	case RL_MT_Light_ID:
		__GluonParamPath=("Parameters_files/Gluons/RL_MT_Light_List.txt");
		p = new C_Gluon(__GluonParamPath);
		break;
	case RL_MT_Heavy_ID:
		__GluonParamPath=("Parameters_files/Gluons/RL_MT_Heavy_List.txt");
		p = new C_Gluon(__GluonParamPath);
		break;
	case RL_MT_Heavy_DD_ID:
		__GluonParamPath=("Parameters_files/Gluons/RL_MT_Heavy_DD_List.txt");
		p = new C_Gluon(__GluonParamPath);
		break;
	case PS_Light_ID:
		__GluonParamPath=("Parameters_files/Gluons/RL_MT_Heavy_DD_List.txt");
		p = new C_Gluon(__GluonParamPath);
		break;
	case Test_Gluon_ID:
		p = new C_Gluon();
		p -> setGluonDefaultParameters();
		break;
	default:
		std::cout << "No such type of Maris-Tandy-like Gluon!" << std::endl;
		assert(false);
	}
	return p;
};

// Parameterized Factory Method function
//----------------------------------------------------------------------
C_Gluon* C_Gluon::createGluon( Gluon_ID id, std::string& _InterpolationPointsPath){
	C_Gluon * p;
	switch (id){
	case Arbitrary_Gluon_ID:
		p = new C_Gluon();
		p -> SetInterpolatorPoints(_InterpolationPointsPath);
		break;
	default:
		std::cout << "No such type of RainbowLadder Gluon!" << std::endl;
		assert(false);
	}
	return p;
};

// Initialization
//----------------------------------------------------------------------
void C_Gluon::InitialState(){
	ReadParameters();
	GluonCheck();
}

// Sets Maris-Tandy-like Gluon with
void C_Gluon::setGluonDefaultParameters(){
	Gluon_ref=&C_Gluon::GluonMT;
	D=0.93;
	w2=0.16;
    LogTail=1;
}

// Load interpolation points for gluon
//----------------------------------------------------------------------
void C_Gluon::SetInterpolatorPoints(std::string& _InterpolationPointsPath){
	Gluon_ref=&C_Gluon::GluonByInterpolation;
	t_cmplxArray2D GluonTempStorage(2);
	std::ifstream GluonInterpStream;
	t_cmplx coordinate, value;
	GluonInterpStream.open(_InterpolationPointsPath);
	if (GluonInterpStream.is_open()){
		GluonInterpStream >> coordinate >> value;
		while(GluonInterpStream.good()){
			GluonTempStorage[0].push_back(coordinate);
			GluonTempStorage[1].push_back(value);
			GluonInterpStream >> coordinate >> value;
		}
	}
	else {std::cout << "Cant open file!" << std::endl; assert(false);}
	GluonInterpStream.close();
	FuncToInterpolate = new Interpolation::Linear<t_cmplx,t_cmplx>(GluonTempStorage[0].size(), &GluonTempStorage[0], &GluonTempStorage[1]);
}

// Read parameters from file
//----------------------------------------------------------------------
void C_Gluon::ReadParameters(){
	std::string line;
	std::ifstream ParamList(GluonParamPath);
	if (ParamList.is_open()){
		while (ParamList.good()){
			ParamList >> line >> D;
			ParamList >> line >> w2;
			ParamList >> line >> LogTail;
		}
		std::cout << "D -" <<"  "<< D <<"  "<< "w2 -" <<"  "<< w2 <<"  "<< "LogTail -" <<"  "<< LogTail  << std::endl;
	}
	else {std::cout << "Cant open file!" << std::endl; assert(false);}
}

// Get value of Gluon at k
//----------------------------------------------------------------------
t_cmplxArray1D C_Gluon::getPropAt(t_cmplx k){
	t_cmplxArray1D return_vec(1,0);
	return_vec[0]=(this->*Gluon_ref)(k);
	return return_vec;
}

// Kernel check on real line
//----------------------------------------------------------------------
void C_Gluon::GluonCheck(){
	double Pu,Pd,x,scale,dp;
	scale = 100;
	Pd=0.01;
	Pu=10.0;
	dp=pow(10,(log10(Pu/Pd)/scale));

	x=Pd;
	ofstream Gluon;
	Gluon.open ("Data_files/Gluon.dat");
	for (int i = 0; i <= scale; i++){
		t_cmplx point = x*x;
		Gluon << real(point) <<"  "<< real(getPropAt(point)[0]) << std::endl;
		x*=dp;
	}
	Gluon.close();
}

// Maris-Tandy gluon model
//----------------------------------------------------------------------
t_cmplx C_Gluon::GluonMT(t_cmplx k){
	const double gamma_m=12.0/(33.0 - 2.0*4.0);
	const double m_t=0.5;
	const double tau=7.389056 - 1.0;
	const double LambdaQCD=0.234;
	const double pi=3.14159265358979;
	return 4.0*pi*pi*D*(k)*exp(-1.0*(k)/w2)/w2/w2/w2  + LogTail*4.0*pi*pi*(gamma_m*(1.0-1.0*exp(-(k)/4.0/m_t/m_t))/(k))/(0.5*log(tau+(1.0+(k)/LambdaQCD/LambdaQCD)*(1.0+(k)/LambdaQCD/LambdaQCD)));
}

// Model given by interpolation over provided points set
//----------------------------------------------------------------------
t_cmplx	C_Gluon::GluonByInterpolation(t_cmplx k){
	return 4.0*pi*FuncToInterpolate->getValue(k)/k;
}


