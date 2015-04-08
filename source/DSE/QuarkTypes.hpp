#pragma once
#include "Quark.hpp"
#include <string>

enum Quark_ID { Up_ID=0, Down_ID, Strange_ID, Charm_ID, Test_ID, Quark_ID_End };

class C_Up_Quark: public C_Quark{
public:
	C_Up_Quark(){
		SetNameID("UP_Quark",1);
		std::string ParamsPath=("Parameters_files/Propagators/UP_quark.txt");
		SavePropPath=("Data_files/SaveQuarkContour_Chiral.dat");
		ReadParameters(ParamsPath);
		InitialState();
		InitializateIntegrators();
	}

	void info() { std::cout << "UP_Quark initialization..." << std::endl; }
};

class C_Down_Quark: public C_Quark{
public:
	void info() { std::cout << "Down_Quark" << std::endl; }
	// TODO implement down quark
};

class C_Strange_Quark: public C_Quark{
public:
	void info() { std::cout << "Strange_Quark" << std::endl; }
	// TODO implement strange quark
};


class C_Charm_Quark: public C_Quark{
public:
	C_Charm_Quark(){
		SetNameID("CHARM_Quark",1);
		std::string ParamsPath=("Parameters_files/Propagators/CHARM_quark.txt");
		SavePropPath=("Data_files/SaveQuarkContour_CHARM.dat");
		ReadParameters(ParamsPath);
		InitialState();
		InitializateIntegrators();
	}

	void info() { std::cout << "Charm_Quark initialization..." << std::endl; }
};

class C_Test_Quark: public C_Quark{
public:
	C_Test_Quark(){
		SetNameID("TEST_Quark",1);
		SavePropPath=("Data_files/SaveQuarkContour_TEST.dat");
		params.setDefault();
		InitialState();
		InitializateIntegrators();
	}

	void info() { std::cout << "Test_Quark initialization..." << std::endl; }
};


