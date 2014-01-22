#pragma once

class C_Manipulator: public C_AbstractClass{
	public:
	C_Physical_State * PhysState;
	
	t_dArray1D z_rad, w_rad;
	bool IfSpurious;
	
	C_Integrator_Line<t_cmplxMatrix,C_Manipulator,double> * Integ_radial_leg;
	
	C_Manipulator(){
		Integ_radial_leg=C_Integrator_Line<t_cmplxMatrix,C_Manipulator,double>::createIntegrator(40, 0.0001, 2000, 1, qgausleg_log_ID);
		Integ_radial_leg->getNodes(&z_rad,&w_rad);
	}
	
	
	void SetPhysicalState(C_Physical_State * __PhysState) { PhysState=__PhysState; }
	C_Physical_State * GetPhysicalState() {return PhysState;}
	
	void DressPropagators(){
		C_Quark * quark_pointer;
		quark_pointer=(C_Quark*)(PhysState->Propagators[0]);
		//MemoryManager->CopyMemoryFrom(quark_pointer->Memory,PhysState->Kernels[0]->Memory);
		PhysState->Kernels[0]->SetConvolutionType(0);
		
		PhysState->Propagators[0]->DressPropagator();
		PhysState->Propagators[1]->DressPropagator();
		
		PhysState->Kernels[0]->SetConvolutionType(1);
	}
	
	void ReDressPropagators(double m1, double m2, double M2_contour){
		PhysState->Propagators[0]->InitialState();
		PhysState->Propagators[1]->InitialState();
		
		PhysState->Propagators[0]->setContourApex(M2_contour);
		PhysState->Propagators[1]->setContourApex(M2_contour);
		
		PhysState->Propagators[0]->params.m0=m1;
		PhysState->Propagators[1]->params.m0=m2;
		DressPropagators();
	}
	
	void DressBSEs(int num_BSE,t_cmplx Q, int steps){
		PhysState->BSEs[num_BSE]->DressBSA(Q,steps);
	}
	
	t_cmplxArray2D CalcEVBSEs(int num_BSE,t_cmplx Q){
		return PhysState->BSEs[num_BSE]->SetEVMatrix(Q);
	}
	
	double NormalizeBSEs(int num_BSE,t_cmplx Q){
		return PhysState->BSEs[num_BSE]->NormalizeBSA(Q);
	}
	
	void DressBSEs_complex(int num_BSE,t_cmplx Q, int steps, double M_contour){
		PhysState->BSEs[num_BSE]->DressBSA_complex(Q,steps,M_contour);
	}
	
	void SetDressingFiles(int num_BSE,t_cmplx Q, const char * SaveDressingPath){
		t_cmplxArray1D E_amplitute;
		PhysState->BSEs[num_BSE]->SetBSAonPath(&E_amplitute,&PhysState->Kernels[0]->Memory->VertexDressings[0][0][0],Q);
		ofstream SaveDressingStream;
		SaveDressingStream.open(SaveDressingPath);
		for (int i = 0; i < E_amplitute.size(); i++){
			SaveDressingStream << PhysState->Kernels[0]->Memory->VertexDressings[0][0][0][i] << '\t' << PhysState->Kernels[0]->Memory->VertexDressings[0][0][1][i] << '\t' << E_amplitute[i] << std::endl;
		}
		SaveDressingStream.close();
		std::cout << "Dressing was saved to file " << (*SaveDressingPath) << std::endl;		
	}
	
	void LoadDressingFiles(){
		t_cmplx z_i,E_amp;
		ofstream temp_continuation;
		temp_continuation.open ("Data_files/E_amp_onCauchy.dat");
		/*for (int j=0;j<= PhysState->Kernels[0]->Memory->VertexDressings[1][0][0].size()/4;j++) 
		{
			z_i=real(PhysState->Kernels[0]->Memory->VertexDressings[1][0][0][j]);
			E_amp=real(PhysState->Kernels[0]->GetDressingAt(1,0,z_i));
			temp_continuation << z_i << '\t' << E_amp << std::endl;
			std::cout << z_i << '\t' << E_amp << std::endl;
		}*/
		temp_continuation.close();
	}
	
	void DrawBSA_complex(double Norm, const char * path_to_bse){
		double Pu,Pd,x,scale,dp;
		scale = 20;
		Pd=0.01;
		Pu=40.0;
		dp=pow(10,(log10(Pu/Pd)/scale));
		x=Pd;
		
		if(PhysState->BSEs[0]->flag_off_shell){
			ofstream DrawBSA;
			DrawBSA.open (path_to_bse);
			for (int i = 0; i <= scale; i++)
			{		
				t_cmplx point = x; 
				this->DressBSEs_complex(0,point,10, 1.0);
				for (int j = 1; j < z_rad.size()-1; j++){
					t_cmplxArray1D amps;
					
					amps=PhysState->BSEs[0]->getCauchyAt_embedded(z_rad[j]);
					
					DrawBSA << real(point*point) << "  " << z_rad[j];
					for (int k = 0; k < amps.size(); k++){
						if (k==0) DrawBSA << "  " << real(amps[k]) - 1.0;
						else DrawBSA << "  " << real(amps[k]);
					}
					DrawBSA << std::endl;
				}				
				x*=dp;
			}
			DrawBSA.close();
		}
		else{
			ofstream DrawBSA;
			DrawBSA.open ("Data_files/DrawBSA_ONshell.dat");
			for (int i = 0; i <= scale; i++)
			{		
				t_cmplx point = x; 
				this->DressBSEs_complex(0,point,13,1.0);
				DrawBSA << real(point*point) << "  " << Norm*real(PhysState->BSEs[0]->getCauchyAt_embedded(0.001)[0]) << std::endl;
				x*=dp;
			}
			DrawBSA.close();
		}
	
	}
	
	void ExportBSA_complex(const char * path_to_store, const char * path_to_draw){
		double Pu,Pd,x,scale,dp;
		scale = 30;
		Pd=0.01;
		Pu=45.0;
		dp=pow(10,(log10(Pu/Pd)/scale));
		x=Pd;
		PhysState->BSEs[0]->setContourAndGrid(1.0);
		if(PhysState->BSEs[0]->flag_off_shell){
			
			ofstream StoreBSA;
			StoreBSA.open (path_to_store);
			
			ofstream DrawBSA;
			DrawBSA.open (path_to_draw);
			
			StoreBSA << scale  << '	' <<  PhysState->BSEs[0]->Memory->CauchyContour[0].size() << std::endl;
			for (int i = 0; i <= scale; i++)
			{		
				t_cmplx point = x; 
				this->DressBSEs_complex(0,point,15,1.0);
				
				for (int j = 1; j < z_rad.size()-1; j++){
					t_cmplxArray1D amps;
					
					amps=PhysState->BSEs[0]->getCauchyAt_embedded(z_rad[j]);
					
					DrawBSA << real(point*point) << "  " << z_rad[j];
					for (int k = 0; k < amps.size(); k++){
						if (k==0) DrawBSA << "  " << real(amps[k]) - 1.0;
						else DrawBSA << "  " << real(amps[k]);
					}
					DrawBSA << std::endl;
				}
				
				for (int j = 0; j < PhysState->BSEs[0]->Memory->CauchyContour[0].size(); j++){					
					StoreBSA << real(point*point)  << '	' <<  PhysState->BSEs[0]->Memory->CauchyContour[0][j]  << '	' <<  PhysState->BSEs[0]->Memory->CauchyContour[1][j];
					for (int k = 0; k < PhysState->BSEs[0]->Memory->CauchyContour.size()-2; k++){
						if (k==0) StoreBSA  << '	' <<  (PhysState->BSEs[0]->Memory->CauchyContour[k+2][j]) - 1.0;
						else StoreBSA  << '	' <<  (PhysState->BSEs[0]->Memory->CauchyContour[k+2][j]);
					}
					StoreBSA << std::endl;
				}				
				x*=dp;
			}
			DrawBSA.close();
			StoreBSA.close();
		}	
	}
	
	void CalcSpectra(double down_limit, double up_limit, const char * path_to_spectra){
		t_cmplxArray3D Spectra(5,t_cmplxArray2D(2));
		
		ofstream CharmSpectra;
		CharmSpectra.open (path_to_spectra);
		for (int i = 0; i < 5; i++){
			int num_exited=0;
			IfSpurious=false;
			std::cout << "Now doing - " << "  " << i << std::endl;
			while (!IfSpurious){
				t_cmplxArray1D state = RootFind(i,num_exited,down_limit,up_limit);
				if (!IfSpurious){
					std::cout << num_exited << "  " << state[0] << "  " << state[1] << std::endl;
					Spectra[i][0].push_back(state[0]);
					Spectra[i][1].push_back(state[1]);
					num_exited++;
				}
			}
			
			for (int j = 0; j < Spectra[i][0].size(); j++){
				CharmSpectra << 2*i + (real(Spectra[i][1][j])+1.0)/2.0 << "  " << imag(Spectra[i][0][j])  << std::endl;
				std::cout << 2*i + (real(Spectra[i][1][j])+1.0)/2.0 << "  " << imag(Spectra[i][0][j]) << "  " << (real(Spectra[i][1][j]))<< std::endl;
			}
		}
		CharmSpectra.close();
	}
	
	t_cmplxArray1D RootFind(int num_bse, int num_state, double down_limit, double up_limit){
		t_cmplx start,end,mid,root,parity;
		t_cmplx s_P,e_P,m_P,r_P;
		double accuracy=0.00001;
		double eps=1.0;
		
		s_P=t_cmplx(0.0,down_limit);e_P=t_cmplx(0.0,up_limit);
		int ctr=0;
		while (eps>accuracy){
			t_cmplxArray2D temp_EV_parity;
			m_P=(e_P+s_P)/2.0;
			temp_EV_parity=CalcEVBSEs(num_bse,m_P);
			mid=temp_EV_parity[0][num_state];
			parity=temp_EV_parity[1][num_state];
			eps = fabs(real(mid) - 1.0);
			std::cout << mid << "  " << m_P << "  " << eps << std::endl;
			if(fabs(imag(mid))>0.1 || real(mid)<0.001 || ctr>30 || fabs(imag(t_cmplx(0.0,up_limit) - m_P)) < 0.1  ) {IfSpurious=true; break;}
			if(real(mid)>1.0) e_P=m_P; 
			if(real(mid)<1.0) s_P=m_P; 
			ctr++;
		}
		r_P=m_P;
		
		t_cmplxArray1D result_EV_parity(2);
		
		result_EV_parity[0]=r_P;
		result_EV_parity[1]=parity;
		
		return result_EV_parity;		
	}
	
	void CheckGroundStates(double M_PS, double M_S, double M_V, double M_AV){
		t_cmplxArray3D temp(4);
		std::cout << "PseudoScalar EV calculation" << std::endl;
		temp[0]=CalcEVBSEs(0,t_cmplx(0.0,M_PS));
		std::cout << "Scalar EV calculation" << std::endl;
		temp[1]=CalcEVBSEs(1,t_cmplx(0.0,M_S));
		std::cout << "Vector EV calculation" << std::endl;
		temp[2]=CalcEVBSEs(2,t_cmplx(0.0,M_V));
		std::cout << "AxialVector EV calculation" << std::endl;
		temp[3]=CalcEVBSEs(3,t_cmplx(0.0,M_AV));
		
		std::vector <string> names(4);
		names[0]="PS-";
		names[1]="S-";
		names[2]="V-";
		names[3]="AV-";
		std::cout << fixed;
		std::cout << setprecision (NUM_PRECISION);
		for (int i = 0; i < 5; i++){
			for (int j = 0; j < 4; j++){
				std::cout << names[j] << "  " << real(temp[j][0][i]) << "  " << real(temp[j][1][i]) << "||   ";
			}
			std::cout << std::endl;
		}
	}
	
	double getPionMass(double m_quark){
		t_dArray2D QuarkAndPionMasses(2);
		double m_q, M_pi;
		
		ifstream QuarkAndPionMassesStream;
		QuarkAndPionMassesStream.open("Data_files/QuarkMassDependence_RL.dat"); 
		QuarkAndPionMassesStream >> m_q >> M_pi;
		int i=0;
		while(QuarkAndPionMassesStream.good())
		{
			QuarkAndPionMasses[0].push_back(m_q);
			QuarkAndPionMasses[1].push_back(M_pi);
			QuarkAndPionMassesStream >> m_q >> M_pi;
			i++;
		}
		QuarkAndPionMassesStream.close();
		Interpolation::Linear<double,double> PionMassToInterpolate(QuarkAndPionMasses[0].size(), &QuarkAndPionMasses[0], &QuarkAndPionMasses[1]);
		return PionMassToInterpolate.getValue(m_quark);
	}
	
	
	/*dcx RootFind_new(C_BSE_Hadron_Base * obj, int num_state){
		dcx start,end,mid,root;
		dcx s_P,e_P,m_P,r_P;
		double accuracy=0.0001;
		double eps=1.0;
		
		//if(_ID==Pion_ID) {s_P=prev_PionMass;e_P=dcx(0.0,LimitBSMass);accuracy=0.0001;}
		//if(_ID==Pho_ID) {s_P=prev_RhoMass;e_P=dcx(0.0,LimitBSMass);accuracy=0.001;}
		s_P=dcx(0.0,2.5);e_P=dcx(0.0,3.7);accuracy=0.001;
		//start=obj->SetEVMatrix(s_P)[0];
		//end=obj->SetEVMatrix(e_P)[0];
		start=obj->SetEVMatrix(s_P)[0][num_state];
		std::cout << start << "  " << s_P << std::endl;
		end=obj->SetEVMatrix(e_P)[0][num_state];
		std::cout << end << "  " << e_P << std::endl;
		while (eps>accuracy){
			
			double gamma=(real(end)-real(start))/(imag(e_P) - imag(s_P));
			m_P=dcx(0.0,1.0)*((1.0-real(start))/gamma + imag(s_P));
			
			mid=obj->SetEVMatrix(m_P)[0][num_state];
			eps = fabs(real(mid) - 1.0);
			std::cout << mid << "  " << m_P << "  " << eps << std::endl;
			if(real(mid)>1.0) e_P=m_P; 
			if(real(mid)<1.0) s_P=m_P; 
		}
		
		r_P=m_P;
		return r_P;		
	}*/
	
	void CalcQuarkMassDependence(double initial_mass, double final_mass, const char * path_quark_mass_dependence){
		double m_quark=initial_mass;
		double M2_contour;
		t_cmplxArray1D MassAndParity(2);
		
		ofstream MassDependenceStream;
		MassDependenceStream.open(path_quark_mass_dependence);
		
		while (m_quark < final_mass ){
			double M2_contour=0.5 + sqrt(m_quark)*2.0;
			PhysState->Kernels[0]->SetMesonExchangeMass(getPionMass(m_quark));
			
			ReDressPropagators(m_quark,m_quark,M2_contour);
			PhysState->BSEs[0]->flag_amp_desciption=false;
			MassAndParity=RootFind(0,0,0.0,1.2);
			
			std::cout << m_quark  << '	' <<  imag(MassAndParity[0]) << std::endl;
			MassDependenceStream << m_quark  << '	' <<  imag(MassAndParity[0]) << std::endl;
			
			m_quark=1.5*m_quark;
		}
		MassDependenceStream.close();
	}
	
	void CalcChargeRadiusQuarkMassDependence(double initial_mass, double final_mass, const char * path_charge_qmd){
		double m_quark=initial_mass;
		double M2_contour,Norm_bs,Mass_bs,ChargeRad_bs;
		
		ofstream ChargeMassDepStream;
		ChargeMassDepStream.open(path_charge_qmd);
		
		while (m_quark < final_mass){
			double M2_contour=0.6 + sqrt(m_quark)*2.0;
			PhysState->Kernels[0]->SetMesonExchangeMass(getPionMass(m_quark));
			
			ReDressPropagators(m_quark,m_quark,M2_contour);
			PhysState->BSEs[0]->flag_amp_desciption=false;
			
			Mass_bs=imag(RootFind(0,0,0.0,1.2)[0]);
			Norm_bs=NormalizeBSEs(0,t_cmplx(0.0,Mass_bs));
			
			if (PhysState->Kernels[0]->GetKernelID() == RL_ID){
				ChargeRad_bs=CalcChargeRadius_RL(Norm_bs,Mass_bs);
			}
			if (PhysState->Kernels[0]->GetKernelID() == RL_PS_ID){
				ChargeRad_bs=CalcChargeRadius_PionCloud(Norm_bs,Mass_bs);
			}
		
			std::cout << m_quark  << '	' <<  Mass_bs  << '	' <<  Norm_bs  << '	' <<  ChargeRad_bs << std::endl;
			ChargeMassDepStream << m_quark  << '	' <<  Mass_bs  << '	' <<  Norm_bs  << '	' <<  ChargeRad_bs << std::endl;
			
			m_quark=1.5*m_quark;
		}
		ChargeMassDepStream.close();
	}
	
	
	void CalcFormFactor_RL(double _Norm,double Mass_bs, const char * path_to_formfactor){
		double Pu,Pd,x,scale,dp;
		scale = 20;
		Pd=0.01;
		Pu=2.45;
		dp=Pu/scale;
		x=Pd;
		
		C_FormFactor_RL FF_1;
		FF_1.SetParams(_Norm,Mass_bs);
		FF_1.IncomingBS = PhysState->BSEs[0];
		FF_1.OutgoingBS = PhysState->BSEs[0];
		PhysState->BSEs[2]->flag_off_shell=true;
		FF_1.TransmitingBS = PhysState->BSEs[2];
		FF_1.Quark = PhysState->Propagators[0];
				
		ofstream FormFactorStream;
		FormFactorStream.open (path_to_formfactor);			
		
		t_cmplx FF_1_res;
		for (int i = 0; i <= scale; i++)
		{		
			FF_1_res=FF_1.GetResultAtPoint(x);
			//std::cout << x*x << "  " << FF_1_res+FF_2_res-FF_3_res << "  " << x*x*(FF_1_res+FF_2_res-FF_3_res) << "  " << FF_1_res << "  " << FF_2_res << "  " << FF_3_res << std::endl;	
			//FormFactorStream << x*x << "  " << FF_1_res+FF_2_res+FF_3_res << "  " << x*x*(FF_1_res+FF_2_res+FF_3_res) << "  " << FF_1_res << "  " << FF_2_res << "  " << FF_3_res << std::endl;		
			std::cout << x*x << "  " << FF_1_res << std::endl;	
			FormFactorStream << x*x << "  " << real(FF_1_res) << std::endl;		
			x+=dp;
		}
		FormFactorStream.close();
	}
	
	
	double CalcChargeRadius_RL(double _Norm, double Mass_bs){
		double Q,h,Deriv_total;
		t_cmplxArray1D deriv_FF(4);
		
		C_FormFactor_RL FF_1;
		FF_1.SetParams(_Norm,Mass_bs);
		FF_1.IncomingBS = PhysState->BSEs[0];
		FF_1.OutgoingBS = PhysState->BSEs[0];
		PhysState->BSEs[2]->flag_off_shell=true;
		FF_1.TransmitingBS = PhysState->BSEs[2];
		FF_1.Quark = PhysState->Propagators[0];
				
		Q=0.0316;
		h=Q/5.0;
		
		for (int i = 1; i <= 2; i++)
		{
			double i_count;
			i_count=i;
			deriv_FF[(i-1)]=FF_1.GetResultAtPoint(Q+i_count*h);
			std::cout << "FormFactor at point " << "  " << Q+i_count*h << "  " << real(deriv_FF[(i-1)]) << std::endl;

			deriv_FF[(i-1)+2]=FF_1.GetResultAtPoint(Q-i_count*h);
			std::cout << "FormFactor at point " << "  " << Q-i_count*h << "  " << real(deriv_FF[(i-1)+2]) << std::endl;
		}
		
		
		Deriv_total=real(-deriv_FF[1] + 8.0*deriv_FF[0] - 8.0*deriv_FF[2] + deriv_FF[3])/(12.0*(h))/(2.0*Q)*0.197*0.197*6.0;
		
		std::cout << "Charge Radius is -" << "  " << Deriv_total << std::endl;
		return Deriv_total;	
	}
	
	void CalcFormFactor_PionCloud(double _Norm, double Mass_bs, const char * path_to_formfactor){
		double Pu,Pd,x,scale,dp;
		scale = 20;
		Pd=0.01;
		Pu=2.45;
		dp=Pu/scale;
		x=Pd;
		
		C_FormFactor_RL FF_1;
		FF_1.SetParams(_Norm,Mass_bs);
		FF_1.IncomingBS = PhysState->BSEs[0];
		FF_1.OutgoingBS = PhysState->BSEs[0];
		PhysState->BSEs[2]->flag_off_shell=true;
		FF_1.TransmitingBS = PhysState->BSEs[2];
		FF_1.Quark = PhysState->Propagators[0];
		
		C_FF_SelfCoupling FF_2;
		FF_2.SetParams(_Norm,Mass_bs);
		FF_2.Vertexes[0] = PhysState->BSEs[0];
		FF_2.Vertexes[1] = PhysState->BSEs[0];
		FF_2.Propagators[0] = PhysState->Propagators[0];
		
		C_FF_Seagull FF_3;
		FF_3.SetParams(_Norm,Mass_bs);
		FF_3.Vertexes[0] = PhysState->BSEs[0];
		FF_3.Vertexes[1] = PhysState->BSEs[0];
		FF_3.Propagators[0] = PhysState->Propagators[0];
		
		ofstream FormFactorStream;
		FormFactorStream.open (path_to_formfactor);			
		
		t_cmplx FF_1_res,FF_2_res,FF_3_res;
		for (int i = 0; i <= scale; i++)
		{		
			FF_1_res=FF_1.GetResultAtPoint(x);
			FF_2_res=FF_2.GetResultAtPoint(x);
			FF_3_res=FF_3.GetResultAtPoint(x);
			std::cout << x*x << "  " << FF_1_res << "  " << FF_2_res << "  " << FF_3_res << std::endl;	
			FormFactorStream << x*x << "  " << real(FF_1_res) << "  " << real(FF_2_res) << "  " << real(FF_3_res) << std::endl;			
			x+=dp;
		}
		FormFactorStream.close();
	}
	
	double CalcChargeRadius_PionCloud(double _Norm, double Mass_bs){
		double Q,h,Deriv_total;
		t_cmplxArray1D deriv_FF(4);
		
		C_FormFactor_RL FF_1;
		FF_1.SetParams(_Norm,Mass_bs);
		FF_1.IncomingBS = PhysState->BSEs[0];
		FF_1.OutgoingBS = PhysState->BSEs[0];
		PhysState->BSEs[2]->flag_off_shell=true;
		FF_1.TransmitingBS = PhysState->BSEs[2];
		FF_1.Quark = PhysState->Propagators[0];
		
		C_FF_SelfCoupling FF_2;
		FF_2.SetParams(_Norm,Mass_bs);
		FF_2.Vertexes[0] = PhysState->BSEs[0];
		FF_2.Vertexes[1] = PhysState->BSEs[0];
		FF_2.Propagators[0] = PhysState->Propagators[0];
		
		C_FF_Seagull FF_3;
		FF_3.SetParams(_Norm,Mass_bs);
		FF_3.Vertexes[0] = PhysState->BSEs[0];
		FF_3.Vertexes[1] = PhysState->BSEs[0];
		FF_3.Propagators[0] = PhysState->Propagators[0];
						
		Q=0.0316;
		h=Q/5.0;
		
		for (int i = 1; i <= 2; i++)
		{
			double i_count;
			i_count=i;
			
			t_cmplx FF_1_res,FF_2_res,FF_3_res;	
			
			FF_1_res=FF_1.GetResultAtPoint(Q+i_count*h);
			FF_2_res=FF_2.GetResultAtPoint(Q+i_count*h);
			FF_3_res=FF_3.GetResultAtPoint(Q+i_count*h);
			
			deriv_FF[(i-1)]=FF_1_res + FF_2_res + FF_3_res;
			std::cout << "FormFactor at point " << "  " << Q+i_count*h << "  " << real(deriv_FF[(i-1)]) << "  " << real(FF_1_res) << "  " << real(FF_2_res) << "  " << real(FF_3_res)  << std::endl;
			
			FF_1_res=FF_1.GetResultAtPoint(Q+i_count*h);
			FF_2_res=FF_2.GetResultAtPoint(Q+i_count*h);
			FF_3_res=FF_3.GetResultAtPoint(Q+i_count*h);
			
			deriv_FF[(i-1)+2]=FF_1_res + FF_2_res + FF_3_res;
			std::cout << "FormFactor at point " << "  " << Q-i_count*h << "  " << real(deriv_FF[(i-1)+2]) << "  " << real(FF_1_res) << "  " << real(FF_2_res) << "  " << real(FF_3_res) << std::endl;
		}
		
		Deriv_total=real(-deriv_FF[1] + 8.0*deriv_FF[0] - 8.0*deriv_FF[2] + deriv_FF[3])/(12.0*(h))/(2.0*Q)*0.197*0.197*6.0;
		
		std::cout << "Charge Radius is -" << "  " << Deriv_total << std::endl;
		return Deriv_total;	
	}
	
	void CalcFormFactorWTI(){
		t_cmplx FF_1_res,FF_2_res,FF_3_res;
		t_cmplx Q=t_cmplx(0.447,0.0);
		
		C_FormFactor_RL FF_1;
		FF_1.IncomingBS = PhysState->BSEs[0];
		FF_1.OutgoingBS = PhysState->BSEs[0];
		PhysState->BSEs[2]->flag_off_shell=true;
		FF_1.TransmitingBS = PhysState->BSEs[2];
		FF_1.Quark = PhysState->Propagators[0];
						
		FF_1_res=FF_1.GetResultAtPoint(Q);

		std::cout << FF_1_res << std::endl;	
	}
};
