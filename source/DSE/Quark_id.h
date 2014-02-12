#pragma once
#include "Quark.hpp"
#include <string>

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

	C_Up_Quark * MakeCopy(){
		return new C_Up_Quark(*this);
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

	C_Charm_Quark * MakeCopy(){
		return new C_Charm_Quark(*this);
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

	C_Test_Quark * MakeCopy(){
		return new C_Test_Quark(*this);
	}

	void info() { std::cout << "Charm_Quark initialization..." << std::endl; }
};


// Parameterized Factory Method function
C_Quark* C_Quark::createQuark( Quark_ID id ){
	C_Quark * p;
	switch (id){
	case Up_ID:
		p = new C_Up_Quark();
		break;
	case Down_ID:
		p = new C_Down_Quark();
		break;
	case Strange_ID:
		p = new C_Strange_Quark();
		break;
	case Charm_ID:
		p = new C_Charm_Quark();
		break;
	case Test_ID:
		p = new C_Test_Quark();
		break;
	default:
		assert( false);
	}
	return p;
};


class C_Quark_Factory: public C_Propagator_Factory{
public:
	C_Propagator* Create(int _id) {
		Quark_ID id=(Quark_ID)_id;
		return C_Quark::createQuark( (id) );
	}
};

