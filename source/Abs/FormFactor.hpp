#pragma once

class C_FormFactor_RL: public C_1_LoopDiagram{
	public:
	double Norm;
	t_cmplx Q_value,M2_bs;
	t_cmplxDirac Photon_m,S_m_m,S_p_p,S_p_m,Pion_p,Pion_m;
	
	C_FormFactor_RL(){
		ReadParameters();
		
		num_amplitudes=1;
		LimDk=0.0001;
		LimUk=2000.0;
		M2_bs=0.138*0.138;
		Norm=20.3;
		
		Initialization();
	}
	
	void SetParams(double _Norm, double _Mass_bs){
		Norm=_Norm;
		M2_bs=_Mass_bs*_Mass_bs;
	}
	
	void SetKinematics(double x, double z, double y, t_cmplx _Q)
	{
		q.SetP4(0.0, sqrt((1.0-z*z)*x*(1.0-y*y)), y*sqrt((1.0-z*z)*x), sqrt(x)*z);
		Q.SetP4(0.0,0.0,_Q,0.0);
		P.SetP4(0.0,0.0,0.0,ii*sqrt(M2_bs + _Q*_Q/4.0));
		P_p=P + Q/2.0;
		P_m=P - Q/2.0;
		q_p=q + P/2.0;
		q_m=q - P/2.0;
		q_p_p=q_p + Q/2.0;
		q_p_m=q_p - Q/2.0;
		k_p=q + Q/4.0;
		k_m=q - Q/4.0;
	}
	
	void SetIntegrationPathes(t_cmplx _Q){
		Memory->ErasePathes();
		std::cout << "Setting Pathes...     ";
		for (unsigned int i = 1; i < zz_rad.size(); i++){
			for (unsigned int j = 1; j < zz_cheb.size() ; j++) {
				for (unsigned int k = 1; k < zz_angleY.size() ; k++) {
					SetKinematics(zz_rad[i],zz_cheb[j],zz_angleY[k],_Q);
					Memory->Path_Pion_p.push_back(k_p*k_p);
					Memory->Path_Pion_m.push_back(k_m*k_m);
					Memory->Path_Quark_p_p.push_back(q_p_p*q_p_p);
					Memory->Path_Quark_p_m.push_back(q_p_m*q_p_m);
				}
				Memory->Path_Photon.push_back(q_p*q_p);
				Memory->Path_Quark_m_m.push_back(q_m*q_m);
			}
		}
		std::cout << "Pathes have been set." << std::endl;
	}
	
	void FillStorages(t_cmplx _Q, t_cmplx M_in, t_cmplx M_out)
	{
		SetIntegrationPathes(_Q);
		IncomingBS->flag_amp_desciption=true;
		IncomingBS->SetBSAonPath(&Memory->Pion_p_Stg,&Memory->Path_Pion_p,M_in,true);
		OutgoingBS->SetBSAonPath(&Memory->Pion_m_Stg,&Memory->Path_Pion_m,M_out,true);
		TransmitingBS->flag_amp_desciption=true;
		TransmitingBS->SetBSAonPath_Cauchy(&Memory->Photon_Stg, &Memory->Path_Photon,_Q,real(sqrt(M2_bs+_Q*_Q/4.0)));
		//TransmitingBS->SetBSAonPath(&Memory->Photon_Stg, &Memory->Path_Photon,_Q,false);
		std::cout << "Quarks" << std::endl;
		
		Quark->SetQuarkonPath(&Memory->Quark_m_m_Stg,&Memory->Path_Quark_m_m);
		Quark->SetQuarkonPath(&Memory->Quark_p_p_Stg,&Memory->Path_Quark_p_p);
		Quark->SetQuarkonPath(&Memory->Quark_p_m_Stg,&Memory->Path_Quark_p_m);
		std::cout << std::endl;
		std::cout << "Storages at " << _Q << " have been precalculated" << std::endl;
		std::cout << "We are ready for calculation of FormFactor at " << _Q << std::endl;
	}
	
	t_cmplxMatrix Integrand(t_cmplxArray1D args) {
		t_cmplxMatrix result;
		result.Resize(num_amplitudes,1);
		std::vector<t_cmplxDirac> Photon_DS(8);
		t_cmplx Kinematic_factor;
		t_cmplx x,y,z;
		x=(args[0]);
		z=args[1];
		y=args[2];
		SetKinematics(real(x),real(z),real(y),Q_value);
		Kinematic_factor=1.0/(16.0*pi*pi*pi)*(x)/(1.0+x/10.0);
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

		result(0,0)= -0.5*Norm*Norm*Kinematic_factor*P/(P*P)*(Pion_p*S_p_p*ii*Photon_m*S_p_m*Pion_m*S_m_m).Tr();

		Int_counter++;
		return result;
	}
	
	t_cmplx GetResultAtPoint(t_cmplx _Q)
	{
		Q_value=_Q;
		t_cmplxMatrix result(num_amplitudes,1);
		FillStorages(_Q,ii*sqrt(M2_bs),ii*sqrt(M2_bs));
		Int_counter=0;
		Sub_Int_counter=0;
		std::cout << "Integration started..." << std::endl;
		result = quad3d();
		return result(0,0);
	}	
};
