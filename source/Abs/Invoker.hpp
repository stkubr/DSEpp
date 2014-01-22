class C_Invoker: public C_AbstractClass{
	public:
	//C_BSE_Binder * binder;
	//C_Manipulator * Control;
	
	
	/*C_Manipulator * Control;
	void SetControl(C_Manipulator * __Control){Control=__Control; }
	C_Manipulator * GetControl() {return Control;}*/
		
	virtual void execute()=0;
};

class C_Charm_Invoker: public C_Invoker{
	public:
	void execute(){
		C_Physical_State * charm, * charm_PS;
		Binder.SetBSEBuilder(&CharmBuilder_RL);
		Binder.ConstructPhysState();
		charm=Binder.GetPhysicalState();
		
		charm->CheckPieces();
		
		Binder.SetBSEBuilder(&CharmBuilder_PS);
		Binder.ConstructPhysState();
		charm_PS=Binder.GetPhysicalState();
		
		charm_PS->CheckPieces();
		
		C_Manipulator CharmManipulator;
		
		const char * path_draw;
		const char * path_store;
		const char * path_spectra;
		const char * path_formfactor;
		path_draw=("../Data_files/DrawBSA_OFFshell_etta.dat");
		path_store=("../Data_files/DrawBSA_OFFshell_etta_STORAGE.dat");
		path_spectra=("../Data_files/CharmSpectra_FWC_November_PS_2.dat");
		path_formfactor=("../Data_files/FormFactor_Charm_PS.dat");

		ifstream charmstream(path_store);
		
		//CharmManipulator.SetPhysicalState(charm);
		//CharmManipulator.DressPropagators();
		//CharmManipulator.CalcFormFactor_RL(20.0,path_formfactor);
		//CharmManipulator.CalcSpectra(path_spectra);
		//charm->BSEs[0]->DrawBSA_matrix(dcx(0.0,3.6),3,1);
		//CharmManipulator.CheckGroundStates();
		//CharmManipulator.DressBSEs(4,dcx(0.0,2.95),10);
		//charm->BSEs[0]->flag_off_shell=true;
		//CharmManipulator.DrawBSA_complex(16.5,path_draw);
		//CharmManipulator.ExportBSA_complex(path_store,path_draw);
		//MemoryManager->CopyMemoryFrom(&charmstream,charm->Kernels[0]->Memory);
		

		MemoryManager->CopyMemoryFrom(&charmstream,charm_PS->Kernels[0]->Memory);
		CharmManipulator.SetPhysicalState(charm_PS);
		charm_PS->Kernels[0]->SetMesonExchangeMass(2.98);
		CharmManipulator.DressPropagators();
		CharmManipulator.CalcSpectra(2.5,4.5,path_spectra);
		//CharmManipulator.CheckGroundStates();
	}
};

class C_Charm_UP_Invoker: public C_Invoker{
	public:
	void execute(){
		C_Physical_State * charm_up;
		Binder.SetBSEBuilder(&CharmUpBuilder_RL);
		Binder.ConstructPhysState();
		charm_up=Binder.GetPhysicalState();
		
		charm_up->CheckPieces();
		
		C_Manipulator CharmManipulator;
		
		const char * path_draw;
		const char * path_store;
		const char * path_spectra;
		//path_draw=("../Data_files/DrawBSA_OFFshell_etta.dat");
		//path_store=("../Data_files/DrawBSA_OFFshell_etta_STORAGE.dat");
		path_spectra=("../Data_files/Charm_Up_Spectra.dat");

		ifstream charmstream(path_store);
		
		CharmManipulator.SetPhysicalState(charm_up);
		CharmManipulator.DressPropagators();
		//CharmManipulator.CalcSpectra(path_spectra);
		
		t_cmplxArray2D temp;
		/*temp=charm_up->BSEs[0]->SetEVMatrix(dcx(0.0,2.0));
		std::cout << temp[0][0] << "  " << temp[1][0] << std::endl;*/
		CharmManipulator.DressBSEs(0,t_cmplx(0.0,0.7),10);
		//charm->BSEs[0]->flag_off_shell=true;
		//CharmManipulator.DrawBSA_complex(16.5,path_draw);
		//CharmManipulator.ExportBSA_complex(path_store,path_draw);
		//MemoryManager->CopyMemoryFrom(&charmstream,charm->Kernels[0]->Memory);
		
		/*MemoryManager->CopyMemoryFrom(&charmstream,charm_PS->Kernels[0]->Memory);
		CharmManipulator.SetPhysicalState(charm_PS);
		charm_PS->Propagators[0]->flag_dressed=false;
		charm_PS->Kernels[0]->SetConvolutionType(0);
		CharmManipulator.DressPropagators();
		charm_PS->Kernels[0]->SetConvolutionType(1);
		CharmManipulator.CalcSpectra(path_spectra);*/
	}
};

class C_Pion_Invoker: public C_Invoker{
	public:
	void execute(){
		C_Physical_State * up_up;
		Binder.SetBSEBuilder(&PionBuilder_PS);
		Binder.ConstructPhysState();
		up_up=Binder.GetPhysicalState();
		
		up_up->CheckPieces();
		
		C_Manipulator up_upManipulator;
		
		const char * path_draw;
		const char * path_store;
		const char * path_spectra;
		const char * path_formfactor;
		const char * path_quark_mass_dependence;
		const char * path_charge_radius_qmd;
		path_draw=("../Data_files/DrawBSA_OFFshell_pion_extra.dat");
		path_store=("../Data_files/DrawBSA_OFFshell_pion_STORAGE.dat");
		path_formfactor=("../Data_files/FormFactor_PION_Pion_Cloud_long.dat");
		path_quark_mass_dependence=("../Data_files/QuarkMassDependence_PionCloud.dat");
		path_charge_radius_qmd=("../Data_files/ChargeRadiusQMD_PionCloud_reduced.dat");
		//ifstream charmstream("../Data_files/Charm_E_amp_onCauchy.dat");
		
		
		
		up_upManipulator.SetPhysicalState(up_up);
		//up_upManipulator.CalcChargeRadiusQuarkMassDependence(0.00048724,0.15,path_charge_radius_qmd);
		//std::cout << up_upManipulator.getPionMass(0.0037) << std::endl;
		//up_upManipulator.DressPropagators();
		//up_upManipulator.CheckGroundStates(0.138,0.585,0.718,0.80);
		//up_upManipulator.CalcFormFactor_PionCloud(28.7,path_formfactor);
		
		//up_upManipulator.SetPhysicalState(up_up);
		//up_upManipulator.CalcQuarkMassDependence(0.00048724,0.15,path_quark_mass_dependence);
		up_upManipulator.CalcChargeRadiusQuarkMassDependence(0.00048724,0.15,path_charge_radius_qmd);
		//up_upManipulator.DressPropagators();
		//up_upManipulator.CheckGroundStates(0.138,0.64,0.76,0.85);
		//up_upManipulator.CalcChargeRadius_RL(20.3,0.138);
		//up_up->BSEs[0]->flag_amp_desciption=true;
		//up_upManipulator.DressBSEs(0,dcx(0.0,0.138),4);
		//up_upManipulator.NormalizeBSEs(0,dcx(0,0.138));
		/*dcxArray2D temp;
		temp=up_up->BSEs[2]->SetEVMatrix(dcx(0.0,0.677));
		std::cout << temp[0][0] << "  " << temp[1][0] << std::endl;
		std::cout << temp[0][1] << "  " << temp[1][1] << std::endl;
		std::cout << temp[0][2] << "  " << temp[1][2] << std::endl;*/
		//up_up->BSEs[0]->flag_off_shell=true;
		//up_upManipulator.ExportBSA_complex(path_store,path_draw);
		//up_up->BSEs[0]->flag_off_shell=false;
		//up_up->BSEs[0]->
		//up_upManipulator.DrawBSA_complex(20.5,path);
		//up_up->Kernels[0]->SetPionSign(1.0);
		//
		/*CharmManipulator.NormalizeBSEs(0,dcx(0.0,2.985));
		CharmManipulator.SetDressingFiles(0,dcx(0.0,2.985),path);
		
		MemoryManager->CopyMemoryFrom(&charmstream,charm->Kernels[0]->Memory);
		
		dcx z_i,E_amp;
		ofstream temp_continuation;
		temp_continuation.open ("../Data_files/E_amp_onCauchy.dat");
		for (int j=0;j<= charm->Kernels[0]->Memory->VertexDressings[1][0][0].size()/4;j++) 
		{
			z_i=real(charm->Kernels[0]->Memory->VertexDressings[1][0][0][j]);
			E_amp=real(charm->Kernels[0]->GetDressingAt(1,0,z_i));
			temp_continuation << z_i << '\t' << E_amp << std::endl;
			std::cout << z_i << '\t' << E_amp << std::endl;
		}
		temp_continuation.close();*/
	}
};
