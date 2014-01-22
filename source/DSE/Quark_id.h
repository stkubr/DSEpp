#pragma once
#include "Quark.h"
class C_Up_Quark: public C_Quark{
	
	
	public:
	
		C_Up_Quark(){
			SetNameID("UP_Quark",1);
			ifstream ParamsList("../Parameters_files/Propagators/UP_quark.txt");
			//const char * PathToSave;
			SavePropPath=("../Data_files/SaveQuarkContour_Chiral.dat");
			//SavePropPath=&_SaveProp;
			ReadParameters(ParamsList);
			InitialState();
			InitializateIntegrators();
			//DressPropagator();
		}
		
		//virtual ~C_Up_Quark();

		C_Up_Quark * MakeCopy(){
			return new C_Up_Quark(*this);
		}
		
		
		void info() { std::cout << "UP_Quark initialization..." << std::endl; }  
		
};
  
class C_Down_Quark: public C_Quark{
	public:
    void info() { std::cout << "Down_Quark" << std::endl; }

    //virtual ~C_Down_Quark();
};
  
class C_Strange_Quark: public C_Quark{
	public:   
    void info() { std::cout << "Strange_Quark" << std::endl; }    

    //virtual ~C_Strange_Quark();
};

class C_Charm_Quark: public C_Quark{
	public:   
	
	C_Charm_Quark(){
			SetNameID("CHARM_Quark",1);
			ifstream ParamsList("../Parameters_files/Propagators/CHARM_quark.txt");
			//const char * PathToSave;
			SavePropPath=("../Data_files/SaveQuarkContour_CHARM.dat");
			//SavePropPath=&_SaveProp;
			ReadParameters(ParamsList);
			InitialState();
			InitializateIntegrators();
			//DressPropagator();
		}

	//virtual ~C_Charm_Quark();
		
		C_Charm_Quark * MakeCopy(){
			return new C_Charm_Quark(*this);
		}
	
	
    void info() { std::cout << "Charm_Quark" << std::endl; }    
};
  
  
// Realization parametrized Factory Method function
C_Quark* C_Quark::createQuark( Quark_ID id ){
    C_Quark * p;
    switch (id)
    {
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
        default:
            assert( false);
    }
    return p;
};


class C_Quark_Factory: public C_Propagator_Factory{
	public:
	C_Propagator* Create(void * _id) {
		Quark_ID * id=(Quark_ID *)_id;
		return C_Quark::createQuark( (*id) );
    }
};

C_Quark_Factory * QuarkFactory = new C_Quark_Factory;

