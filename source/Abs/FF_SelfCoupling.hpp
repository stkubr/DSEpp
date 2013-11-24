#pragma once 

class C_FF_SelfCoupling: public C_2_LoopDiagram{
	public:
	
	t_cmplxVector Q,P,P_p,P_m;
	t_cmplxVector q_m_m,q_p_p,q_m_p,q_p_m;
	t_cmplxVector q_B_p,q_B_m;
	t_cmplx Q_value,M2_bs;
	t_cmplxDirac S_m_m,S_p_p,S_p_m,S_m_p,Pion_p,Pion_m;
	double Norm;
	
	C_FF_SelfCoupling(){
		num_amplitudes=1;
		LimDk=0.0001;
		LimUk=2000.0;
		
		Q_value=0.01;
		M2_bs=0.138*0.138;
		Norm=29.1;
		
		Initialization();
		SetNumVertexesAndPropagators(2,1);
		Memory->SetNumPathesAndStorages(7,6);
	}
	
	virtual C_FF_SelfCoupling * MakeCopy()
	{
		return new C_FF_SelfCoupling(*this);
	}
	
	void SetParams(double _Norm, double _Mass_bs){
		Norm=_Norm;
		M2_bs=_Mass_bs*_Mass_bs;
	}	
	
	void SetTotalKinematics(t_dArray1D args){
		double x1,z1,y1;
		double x2,z2,y2,phi2;
		x1=args[0];z1=args[1];y1=args[2];
		x2=args[3];z2=args[4];y2=args[5];phi2=args[6];
		
		SetExternalKinematics(Q_value,M2_bs,M2_bs);
		
		Set1LoopKinematics(x1,z1,y1);
		
		Set2LoopKinematics(x2,z2,y2,phi2);
		
		q_B_p=(k1+k2+P)/2.0;
		q_B_m=(k1+k2-P)/2.0;
	}
	
	void SetExternalKinematics(t_cmplx _Q, t_cmplx M2_in, t_cmplx M2_out){
		Q.SetP4(0.0,0.0,_Q,0.0);
		P.SetP4(0.0,0.0,0.0,ii*sqrt(M2_bs + _Q*_Q/4.0));
		P_p=P + Q/2.0;
		P_m=P - Q/2.0;
	}
	
	void Set1LoopKinematics(double x1, double z1, double y1){
		k1.SetP4(0.0, sqrt((1.0-z1*z1)*x1*(1.0-y1*y1)), y1*sqrt((1.0-z1*z1)*x1), sqrt(x1)*z1);
		q_m_m=k1 - Q/4.0 - P/2.0;
		q_p_p=k1 + Q/4.0 + P/2.0;
	}
	
	void Set2LoopKinematics(double x2, double z2, double y2, double phi2){
		k2.SetP4(sqrt((1.0-z2*z2)*x2*(1.0-y2*y2))*sin(phi2), sqrt((1.0-z2*z2)*x2*(1.0-y2*y2))*cos(phi2), y2*sqrt((1.0-z2*z2)*x2), sqrt(x2)*z2);
		q_p_m=k2 + Q/4.0 - P/2.0;
		q_m_p=k2 - Q/4.0 + P/2.0;
	}
	
	
	void SetTotalPathes(t_cmplx _Q){
		Q_value=_Q;
		Int_counter_total=0;
		Int_counter_1=0;
		vector<int> counter(7);
		for (counter[0] = 1; counter[0] < zz_vector[0].size(); counter[0]++){
			for (counter[1] = 1; counter[1] < zz_vector[1].size(); counter[1]++){
				for (counter[2] = 1; counter[2] < zz_vector[2].size(); counter[2]++){
					Int_counter_1++;
					Int_counter_2=0;
					for (counter[3] = 1; counter[3]< zz_vector[3].size(); counter[3]++){
						for (counter[4] = 1; counter[4] < zz_vector[4].size(); counter[4]++){
							for (counter[5] = 1; counter[5] < zz_vector[5].size(); counter[5]++){
								for (counter[6] = 1; counter[6] < zz_vector[6].size(); counter[6]++){
									SetTotalKinematics(Memory_Loop->PointsAndWieghts[0][Int_counter_total]);
									//Memory->Pathes[6].push_back((k1+k2+P)*(k1+k2+P)/4.0);
									//cout << Memory->Pathes[6][Int_counter_total] << endl;
									Int_counter_2++;
									Int_counter_total++;
								}
							}
						}
					}
				}
			}
		}
	}
	
	void Set1LoopPathes(t_cmplx _Q){
		Q_value=_Q;
		SetExternalKinematics(Q_value,M2_bs,M2_bs);
		Int_counter_1=0;
		vector<int> ctr(3);
		for (ctr[0] = 1; ctr[0] < zz_vector[0].size(); ctr[0]++){
			for (ctr[1] = 1; ctr[1] < zz_vector[1].size(); ctr[1]++){
				for (ctr[2] = 1; ctr[2] < zz_vector[2].size(); ctr[2]++){
					Set1LoopKinematics(zz_vector[0][ctr[0]],zz_vector[1][ctr[1]],zz_vector[2][ctr[2]]);
					Memory->Pathes[0].push_back(q_m_m*q_m_m);
					Memory->Pathes[1].push_back(q_p_p*q_p_p);
					Memory->Pathes[2].push_back(k1*k1);
					//cout << Memory->Pathes[0][Int_counter_1] << "  " << Memory->Pathes[1][Int_counter_1] << "  " << Memory->Pathes[2][Int_counter_1] << endl;
					Int_counter_1++;	
				}
			}
		}
	}
	
	void Set2LoopPathes(t_cmplx _Q){
		Q_value=_Q;
		SetExternalKinematics(Q_value,M2_bs,M2_bs);
		Int_counter_2=0;
		vector<int> ctr(7);
		for (ctr[3] = 1; ctr[3]< zz_vector[3].size(); ctr[3]++){
			for (ctr[4] = 1; ctr[4] < zz_vector[4].size(); ctr[4]++){
				for (ctr[5] = 1; ctr[5] < zz_vector[5].size(); ctr[5]++){
					for (ctr[6] = 1; ctr[6] < zz_vector[6].size(); ctr[6]++){
						Set2LoopKinematics(zz_vector[3][ctr[3]],zz_vector[4][ctr[4]],zz_vector[5][ctr[5]],zz_vector[6][ctr[6]]);
						Memory->Pathes[3].push_back(q_p_m*q_p_m);
						Memory->Pathes[4].push_back(q_m_p*q_m_p);
						Memory->Pathes[5].push_back(k2*k2);
						//cout << Memory->Pathes[3][Int_counter_2] << "  " << Memory->Pathes[4][Int_counter_2] << "  " << Memory->Pathes[5][Int_counter_2] << endl;
						Int_counter_2++;
					}
				}
			}
		}
	}
	
	void FillStorages(t_cmplx _Q, t_cmplx M_in, t_cmplx M_out)
	{
		Memory->ErasePathesAndStorages();
		Memory->SetNumPathesAndStorages(7,6);
		
		Set1LoopPathes(_Q);
		Set2LoopPathes(_Q);
		
		cout << "Pathes are set."<< endl;
		
		Vertexes[0]->flag_amp_desciption=true;
		//cout << Memory->Pathes[2].size() << endl;
		Vertexes[0]->SetBSAonPath(&Memory->Storages[2],&Memory->Pathes[2],M_in,true);
		Vertexes[1]->SetBSAonPath(&Memory->Storages[5],&Memory->Pathes[5],M_out,true);
		
		cout << "Quarks" << endl;
		
		Propagators[0]->SetQuarkonPath(&Memory->Storages[0],&Memory->Pathes[0]);
		Propagators[0]->SetQuarkonPath(&Memory->Storages[1],&Memory->Pathes[1]);
		Propagators[0]->SetQuarkonPath(&Memory->Storages[3],&Memory->Pathes[3]);
		Propagators[0]->SetQuarkonPath(&Memory->Storages[4],&Memory->Pathes[4]);
				
		cout << endl;
		cout << "Storages at " << _Q << " have been precalculated" << endl;
		cout << "We are ready for calculation of FormFactor at " << _Q << endl;
	}
	
	/*dcxMatrix Integrand(dArray1D args) {
		dcxMatrix result;
		result.Resize(num_amplitudes,1);
		vector<dcxDirac> Photon_DS(8);
		dcx Kinematic_factor;
		dcx x,y,z;
		x=(args[0]);
		z=args[1];
		y=args[2];
		SetKinematics(real(x),real(z),real(y),Q_value);
		Kinematic_factor=1.0/(16.0*pi*pi*pi)*(x)/(1.05+x/100.0);
		if(flag_sigma){ 
			TransmitingBS->SetDiracStructures(q_p,Q,&Photon_DS);
			S_m_m=-1.0*ii*(q_m*Z)*Memory->Quark_m_m_Stg[Sub_Int_counter](0,0) + I*Memory->Quark_m_m_Stg[Sub_Int_counter](1,0);
			Photon_m=Photon_DS[0];
			Photon_m.Zero();
			for (int i = 0; i < TransmitingBS->num_amplitudes; i++){
				Photon_m += Photon_DS[i]*Memory->Photon_Stg[Sub_Int_counter](i,0);
			}
			Sub_Int_counter++; flag_sigma=false;
		}
		S_p_p=-1.0*ii*(q_p_p*Z)*Memory->Quark_p_p_Stg[Int_counter](0,0) + I*Memory->Quark_p_p_Stg[Int_counter](1,0);
		S_p_m=-1.0*ii*(q_p_m*Z)*Memory->Quark_p_m_Stg[Int_counter](0,0) + I*Memory->Quark_p_m_Stg[Int_counter](1,0);
		Pion_p=ii*Y5*Memory->Pion_p_Stg[Int_counter](0,0) + Y5*(P_p*Z)*Memory->Pion_p_Stg[Int_counter](1,0) + Y5*(k_p*Z)*(k_p*P_p)*Memory->Pion_p_Stg[Int_counter](2,0) + Y5*((k_p*SIG)*P_p)*Memory->Pion_p_Stg[Int_counter](3,0);
		Pion_m=ii*Y5*Memory->Pion_m_Stg[Int_counter](0,0) + Y5*(P_m*Z)*Memory->Pion_m_Stg[Int_counter](1,0) + Y5*(k_m*Z)*(k_m*P_m)*Memory->Pion_m_Stg[Int_counter](2,0) + Y5*((k_m*SIG)*P_m)*Memory->Pion_m_Stg[Int_counter](3,0);

		result(0,0)= -0.5*Pion_norm*Pion_norm*Kinematic_factor*P/(P*P)*(Pion_p*S_p_m*ii*Photon_m*S_p_p*Pion_m*S_m_m).Tr();

		Int_counter++;
		return result;
	}*/
	
	t_cmplxMatrix Integrand(t_dArray1D args,int Int_counter_total) {
		t_cmplxMatrix result; result.Resize(1,1);
		t_cmplx Kinematic_factor,res;
		t_cmplx B_p,B_m;
		Kinematic_factor=args[0]/(16.0*pi*pi*pi)*args[3]/(16.0*pi*pi*pi*pi);
		SetTotalKinematics(args);
		Int_counter_1=Memory_Loop->Counters[0][Int_counter_total];
		Int_counter_2=Memory_Loop->Counters[1][Int_counter_total];
		//cout << Int_counter_total << "  " << Int_counter_1 << "  " << Int_counter_2 << endl;
		//cin.get();
		
		S_m_m=-1.0*ii*(q_m_m*Z)*Memory->Storages[0][Int_counter_1](0,0) + I*Memory->Storages[0][Int_counter_1](1,0);
		S_p_p=-1.0*ii*(q_p_p*Z)*Memory->Storages[1][Int_counter_1](0,0) + I*Memory->Storages[1][Int_counter_1](1,0);
		S_p_m=-1.0*ii*(q_p_m*Z)*Memory->Storages[3][Int_counter_2](0,0) + I*Memory->Storages[3][Int_counter_2](1,0);
		S_m_p=-1.0*ii*(q_m_p*Z)*Memory->Storages[4][Int_counter_2](0,0) + I*Memory->Storages[4][Int_counter_2](1,0);
		
		Pion_p=ii*Y5*Memory->Storages[2][Int_counter_1](0,0) + Y5*(P_p*Z)*Memory->Storages[2][Int_counter_1](1,0) + Y5*(k1*Z)*(k1*P_p)*Memory->Storages[2][Int_counter_1](2,0) + Y5*((k1*SIG)*P_p)*Memory->Storages[2][Int_counter_1](3,0);
		
		Pion_m=ii*Y5*Memory->Storages[5][Int_counter_2](0,0) + Y5*(P_m*Z)*Memory->Storages[5][Int_counter_2](1,0) + Y5*(k2*Z)*(k2*P_m)*Memory->Storages[5][Int_counter_2](2,0) + Y5*((k2*SIG)*P_m)*Memory->Storages[5][Int_counter_2](3,0);
		
		t_cmplx PiProp_p,PiProp_m;
		t_cmplx k2_p,k2_m;
		k2_p=(Q/2.0-k1+k2)*(Q/2.0-k1+k2);
		k2_m=(-1.0*Q/2.0-k1+k2)*(-1.0*Q/2.0-k1+k2);
		PiProp_p=1.0/(M2_bs + k2_p)/0.093;
		PiProp_m=1.0/(M2_bs + k2_m)/0.093;
		
		B_p=Propagators[0]->getPropAt(q_B_p*q_B_p)[1]*PiProp_p;
		B_m=Propagators[0]->getPropAt(q_B_m*q_B_m)[1]*PiProp_m;		
		
		result(0,0)=Norm*Norm*Kinematic_factor*((S_m_m*Pion_p*S_p_p) *Y5*(B_p)* (S_m_p*Pion_m*S_p_m) * Y5*(B_m)).Tr()/4.0;
		
		return result;
	}
	
	t_cmplx GetResultAtPoint(t_cmplx _Q){
		t_cmplx res;
		//res=quad7d()(0,0);
		//cout << res/10.0/10.0/2.0/pi/pi << "  " << Int_counter_total << "  " << Int_counter_1 << "  " << Int_counter_2 << endl;
		Q_value=_Q;
		SetIntegrationVectors(NumPointsTotal);
		
		FillStorages(Q_value,ii*sqrt(M2_bs),ii*sqrt(M2_bs));
		//SetTotalPathes(1.0);
		//SetIntegrationPathes();
		
		res=quad7d_vector()(0,0);
		return res;
	}
	
};
