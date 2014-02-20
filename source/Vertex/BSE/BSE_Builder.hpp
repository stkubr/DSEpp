#pragma once

//#include "../../DSE/Propagator.hpp"



class C_Manipulator;

// Product 
class C_Physical_State{
	public:
	std::vector<C_Propagator*> Propagators;
	std::vector<C_AbstractKernel*> Kernels;
	std::vector<C_BSE_Hadron_Base*>  BSEs;
	
	public:
	C_Physical_State() {}
	~C_Physical_State() {}
	
	friend class C_Manipulator;
	
	void SetPropagators(std::vector<C_Propagator*> _Propagators){
		Propagators=_Propagators;
	}
	
	void SetKernels(std::vector<C_AbstractKernel*> _Kernels){
		Kernels=_Kernels;
	}
	
	void SetBSEs(std::vector<C_BSE_Hadron_Base*> _BSEs){
		BSEs=_BSEs;
	}
	
	void CheckPieces(){
		for (int i = 0; i < Propagators.size(); i++){
			Propagators[i]->GetNameID();
		}
		
		for (int i = 0; i < Kernels.size(); i++){
			Kernels[i]->GetNameID();
		}
		
		for (int i = 0; i < BSEs.size(); i++){
			BSEs[i]->GetNameID();
		}
	}	
};


// Abstract builder
class C_BSE_Abstract_Builder{
	public:
	C_Physical_State * PhysicalState;
	
	C_Quark_Factory * QuarkFactory;
	C_Gluon_Factory * GluonFactory;
	C_Kernel_Factory * KernelFactory;

	C_BSE_Abstract_Builder() {
		QuarkFactory = new C_Quark_Factory;
		GluonFactory = new C_Gluon_Factory;
		KernelFactory = new C_Kernel_Factory;
	}
	~C_BSE_Abstract_Builder() {}
	
	C_Physical_State * GetPhysicalState() {return PhysicalState;} 
	
	void createNewPhysicalState() {PhysicalState= new C_Physical_State();}
	
	virtual void buildPropagators()=0;
	virtual void buildKernels()=0;
	virtual void buildBSEs()=0;
	virtual void LinkThemAll()=0;
};


// Concrete builder
class C_Meson: public C_BSE_Abstract_Builder{
	public:
	Kernel_ID kernel_ID;
	Quark_ID quark_1_ID;
	Quark_ID quark_2_ID;
	Gluon_ID gluon_ID;
	PS_type_ID exchange_type_ID;
	
	std::vector<C_Propagator*> Propagators;
	std::vector<C_AbstractKernel*> Kernels;
	std::vector<C_BSE_Hadron_Base*>  BSEs;
	
	C_Meson(Quark_ID Q1_id, Quark_ID Q2_id, Kernel_ID K1_id, Gluon_ID G1_id, PS_type_ID Ex_id){
		quark_1_ID=Q1_id;
		quark_2_ID=Q2_id;
		kernel_ID=K1_id;
		gluon_ID=G1_id;
		exchange_type_ID=Ex_id;
		Propagators.resize(2);
		Kernels.resize(1);
		BSEs.resize(5);
	}
	~C_Meson(){ }
	
	void SymmetricQuarksDetector(){
		if(quark_1_ID==quark_2_ID){
			Propagators[0]=QuarkFactory->Create(quark_1_ID);
			Propagators[1]=Propagators[0];
		}
		else {
			Propagators[0]=QuarkFactory->Create(quark_1_ID);
			Propagators[1]=QuarkFactory->Create(quark_2_ID);
		}
	}
	
	void buildPropagators(){
		SymmetricQuarksDetector();
		PhysicalState->SetPropagators(Propagators);
	}
	
	void buildKernels(){
		//vector<C_AbstractKernel*> Kernels(1);
		Kernels[0]=KernelFactory->Create(kernel_ID);
		PhysicalState->SetKernels(Kernels);
	}
	
	void buildBSEs(){
		//vector<C_BSE_Hadron_Base*> BSEs(1);
		for (int i = 0; i < Dirac_ID_End; i++){
			BSEs[i]=C_BSE_Hadron_Meson::createMesonBSE((Dirac_ID)(i));
			//std::cout << i << std::endl;
		}
		PhysicalState->SetBSEs(BSEs);
	}	
	
	void LinkThemAll(){
		Propagators[0]->LinkToKernel(Kernels[0]);
		Propagators[1]->LinkToKernel(Kernels[0]);
		for (int i = 0; i < BSEs.size(); i++){
			BSEs[i]->LinkToKernel(Kernels[0]);
			BSEs[i]->LinkToPartons(Propagators[0],Propagators[1]);
		}
	}
};

// Director
class C_BSE_Binder{
	private:
	C_BSE_Abstract_Builder * BSE_Builder;
	public:
	C_BSE_Binder() : BSE_Builder(NULL){	}
	~C_BSE_Binder() {}
	
	void SetBSEBuilder(C_BSE_Abstract_Builder * b) {BSE_Builder=b;}
	
	C_Physical_State * GetPhysicalState() {return BSE_Builder->GetPhysicalState();}
	
	void ConstructPhysState(){
		BSE_Builder->createNewPhysicalState();
		BSE_Builder->buildPropagators();
		BSE_Builder->buildKernels();
		BSE_Builder->buildBSEs();
		BSE_Builder->LinkThemAll();
	}
};

C_BSE_Binder Binder;
C_Meson PionBuilder_RL(Up_ID, Up_ID, RL_ID, RL_MT_Light_ID, Pion_exchange_ID);

C_Meson PionBuilder_PS(Up_ID, Up_ID, RL_PS_ID, PS_Light_ID, Pion_exchange_ID);

C_Meson CharmBuilder_RL(Charm_ID, Charm_ID, RL_ID, RL_MT_Heavy_ID, Etta_exchange_ID);

C_Meson CharmUpBuilder_RL(Charm_ID, Up_ID, RL_ID, RL_MT_Heavy_DD_ID, Etta_exchange_ID);

C_Meson CharmBuilder_PS(Charm_ID, Charm_ID, RL_PS_ID, RL_MT_Heavy_ID, Etta_exchange_ID);










