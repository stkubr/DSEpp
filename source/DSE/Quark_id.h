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
			ReadParameters(&ParamsList);
			InitialState();
			InitializateIntegrators();
			//DressPropagator();
		}
		
		C_Up_Quark * MakeCopy(){
			return new C_Up_Quark(*this);
		}
		
		
		void info() { cout << "UP_Quark initialization..." << endl; }  
		
};
  
class C_Down_Quark: public C_Quark{
	public:
    void info() { cout << "Down_Quark" << endl; }    
};
  
class C_Strange_Quark: public C_Quark{
	public:   
    void info() { cout << "Strange_Quark" << endl; }    
};

class C_Charm_Quark: public C_Quark{
	public:   
	
	C_Charm_Quark(){
			SetNameID("CHARM_Quark",1);
			ifstream ParamsList("../Parameters_files/Propagators/CHARM_quark.txt");
			//const char * PathToSave;
			SavePropPath=("../Data_files/SaveQuarkContour_CHARM.dat");
			//SavePropPath=&_SaveProp;
			ReadParameters(&ParamsList);
			InitialState();
			InitializateIntegrators();
			//DressPropagator();
		}
		
		C_Charm_Quark * MakeCopy(){
			return new C_Charm_Quark(*this);
		}
	
	
    void info() { cout << "Charm_Quark" << endl; }    
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

