#pragma once

class C_1_LoopDiagram: public C_AbsDiagram,public C_1_Loop_Int{
	public:
	C_DedicMem_1_LoopDiagram * Memory;
	
	C_BSE_Hadron_Base * IncomingBS;
	C_BSE_Hadron_Base * OutgoingBS;
	C_BSE_Hadron_Base * TransmitingBS;
	C_Propagator * Quark;
	
	int num_amplitudes,NumRadialPoints,NumChebPoints,NumAngleY,Sub_Int_counter;
	double LimDk,LimUk;
	
	/*dcxDirac S_m_m,S_p_m,S_p_p;
	dcxDirac IncomingBS_p,OutgoingBS_m;
	dcxDirac Photon_m;*/
	t_cmplxVector q,Q,P,P_p,P_m,q_p,q_m,q_p_m,q_p_p,k_p,k_m,k_T;
	t_dArray1D zz_rad,zz_cheb,zz_angleY,w_rad,w_cheb,w_angleY;
	//dcxArray1D Path_TransmitingBS,Path_IncomingBS,Path_OutgoingBS,Path_Quark_p_m,Path_Quark_p_p,Path_Quark_m_m;
	//dcxArray1D Path_Photon,Path_Pion_p,Path_Pion_m,Path_Quark_p_m,Path_Quark_p_p,Path_Quark_m_m;
	
	C_1_LoopDiagram(){
		Memory=(C_DedicMem_1_LoopDiagram*)DedicMemFactory_1_LoopDiagram->CreateMemory();
		
	}
	
	
	void ReadParameters()
	{
		std::cout << "Parameters initialization..."  << std::endl;
		string line;
		ifstream ParamList("Parameters_files/1LoopDiagram.txt");
		if (ParamList.is_open())
		{
			while ( ParamList.good() )
			{
				ParamList >> line >> NumRadialPoints;
				ParamList >> line >> NumChebPoints;
				ParamList >> line >> NumAngleY;
			}
			std::cout << NumRadialPoints << "  " << NumChebPoints << "  " << NumAngleY << std::endl;
		}
		else std::cout << "Cant open file!" << std::endl;
	}
	
	void Initialization(){
		
		Integ_radial=C_Integrator_Line<t_cmplxMatrix,C_1_Loop_Int,double>::createIntegrator(NumRadialPoints, LimDk, LimUk, num_amplitudes, qgausleg_log_ID);
		Integ_angle_cheb=C_Integrator_Line<t_cmplxMatrix,C_1_Loop_Int,double>::createIntegrator(NumChebPoints, LimDk, LimUk, num_amplitudes, qgauscheb_ID);
		Integ_angle_Y=C_Integrator_Line<t_cmplxMatrix,C_1_Loop_Int,double>::createIntegrator(NumAngleY, -1.0, 1.0, num_amplitudes, qgausleg_sym_ID);
		
		Integ_radial->getNodes(&zz_rad,&w_rad);
		Integ_angle_cheb->getNodes(&zz_cheb,&w_cheb);
		Integ_angle_Y->getNodes(&zz_angleY,&w_angleY);
	}
	
	/// This method does something
	/*void ErasePathes()
	{
		Path_Photon.clear();
		Path_Pion_p.clear();
		Path_Pion_m.clear();
		Path_Quark_p_m.clear();
		Path_Quark_p_p.clear();
		Path_Quark_m_m.clear();
	}*/
	
	virtual void SetKinematics(double x, double z, double y, t_cmplx _Q)=0;
	
	virtual void SetIntegrationPathes(t_cmplx _Q)=0;
	
	virtual void FillStorages(t_cmplx _Q, t_cmplx M_in, t_cmplx M_out)=0;
		
};
