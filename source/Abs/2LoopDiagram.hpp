#pragma once 

class C_2_LoopDiagram: public C_AbsDiagram,public C_2_Loop_Int{
	public:
	C_DedicMem_2_LoopDiagram * Memory;
	
	std::vector<C_BSE_Hadron_Base *> Vertexes;
	std::vector<C_Propagator *> Propagators;
	
	int num_amplitudes,NumRadialPoints_1,NumChebPoints_1,NumAngleY_1,Sub_Int_counter,NumPointsTotal;
	int NumRadialPoints_2,NumChebPoints_2,NumAngleY_2, NumAnglePhi_2;
	double LimDk,LimUk;
	
	t_cmplxVector k1,k2;

	
	C_2_LoopDiagram(){
		Memory=(C_DedicMem_2_LoopDiagram*)DedicMemFactory_2_LoopDiagram->CreateMemory();
		ReadParameters();
		zz_vector.resize(7);
		w_vector.resize(7);
	}
	
	virtual C_2_LoopDiagram * MakeCopy()
	{
		return new C_2_LoopDiagram(*this);
	}
	
	void ReadParameters()
	{
		std::cout << "Parameters initialization..."  << std::endl;
		string line;
		ifstream ParamList("Parameters_files/2LoopDiagram.txt");
		if (ParamList.is_open())
		{
			while ( ParamList.good() )
			{
				ParamList >> line >> NumRadialPoints_1;
				ParamList >> line >> NumChebPoints_1;
				ParamList >> line >> NumAngleY_1;
				
				ParamList >> line >> NumRadialPoints_2;
				ParamList >> line >> NumChebPoints_2;
				ParamList >> line >> NumAngleY_2;
				ParamList >> line >> NumAnglePhi_2;
			}
			std::cout << NumRadialPoints_1 << "  " << NumChebPoints_1 << "  " << NumAngleY_1 << "  " << NumRadialPoints_2 << "  " << NumChebPoints_2 << "  " << NumAngleY_2 << "  " << NumAnglePhi_2 << std::endl;
		}
		else std::cout << "Cant open file!" << std::endl;
	}
	
	void Initialization(){
		
		Integ_radial_1=C_Integrator_Line<t_cmplxMatrix,C_2_Loop_Int,double>::createIntegrator(NumRadialPoints_1, LimDk, LimUk, num_amplitudes, qgausleg_log_ID);
		Integ_angle_cheb_1=C_Integrator_Line<t_cmplxMatrix,C_2_Loop_Int,double>::createIntegrator(NumChebPoints_1, LimDk, LimUk, num_amplitudes, qgauscheb_ID);
		Integ_angle_Y_1=C_Integrator_Line<t_cmplxMatrix,C_2_Loop_Int,double>::createIntegrator(NumAngleY_1, -1.0, 1.0, num_amplitudes, qgausleg_sym_ID);
		
		Integ_radial_2=C_Integrator_Line<t_cmplxMatrix,C_2_Loop_Int,double>::createIntegrator(NumRadialPoints_2, LimDk, LimUk, num_amplitudes, qgausleg_log_ID);
		Integ_angle_cheb_2=C_Integrator_Line<t_cmplxMatrix,C_2_Loop_Int,double>::createIntegrator(NumChebPoints_2, LimDk, LimUk, num_amplitudes, qgauscheb_ID);
		Integ_angle_Y_2=C_Integrator_Line<t_cmplxMatrix,C_2_Loop_Int,double>::createIntegrator(NumAngleY_2, -1.0, 1.0, num_amplitudes, qgausleg_sym_ID);
		Integ_angle_Phi_2=C_Integrator_Line<t_cmplxMatrix,C_2_Loop_Int,double>::createIntegrator(NumAnglePhi_2, 0, 2.0*pi, num_amplitudes, qgausleg_sym_ID);
		
		Integ_radial_1->getNodes(&zz_vector[0],&w_vector[0]);
		Integ_angle_cheb_1->getNodes(&zz_vector[1],&w_vector[1]);
		Integ_angle_Y_1->getNodes(&zz_vector[2],&w_vector[2]);
		
		Integ_radial_2->getNodes(&zz_vector[3],&w_vector[3]);
		Integ_angle_cheb_2->getNodes(&zz_vector[4],&w_vector[4]);
		Integ_angle_Y_2->getNodes(&zz_vector[5],&w_vector[5]);
		Integ_angle_Phi_2->getNodes(&zz_vector[6],&w_vector[6]);
		
		NumPointsTotal=NumRadialPoints_1 * NumChebPoints_1 * NumAngleY_1 * NumRadialPoints_2 * NumChebPoints_2 * NumAngleY_2 * NumAnglePhi_2;
	}
	
	/*virtual C_2_LoopDiagram * MakeCopy()
	{
		return new C_2_LoopDiagram(*this);
	}*/
	
	/*dcxMatrix Integrand(dArray1D args) {
		dcxMatrix dummy; dummy.Resize(1,1);
		dcx Kinematic_factor;
		Kinematic_factor=args[0]/(16.0*pi*pi*pi)*args[3]/(32.0*pi*pi*pi*pi);
		//dcxVector k1,k2;
		k1.SetP4(0,0,0,sqrt(args[0]));
		k2.SetP4(0,0,0,sqrt(args[3]));
		
		//dummy(0,0)=(q1*q1)*(q2*q2);
		dummy(0,0)=(((k1*Z + I)*Y5*(k1*Z + I)) *Y5* ((k2*Z + I)*Y5*(k2*Z + I)) * Y5).Tr()/4.0;
		
		Int_counter_total++; 
		return dummy;
	}
	
	
	void GetResult(){
		dcx res;
		//res=quad7d()(0,0);
		//std::cout << res/10.0/10.0/2.0/pi/pi << "  " << Int_counter_total << "  " << Int_counter_1 << "  " << Int_counter_2 << std::endl;
		
		SetIntegrationVectors(NumPointsTotal);
		
		res=quad7d_vector()(0,0);
		std::cout << res << "  " << Int_counter_total << "  " << Int_counter_1 << "  " << Int_counter_2 << std::endl;
		//cin.get();
	}*/
	
	void SetNumVertexesAndPropagators(int num_vertexes, int num_propagators){
		Vertexes.resize(num_vertexes);
		Propagators.resize(num_propagators);
	}
	
	//virtual void SetKinematics(dArray1D int_args)=0;
	
	//virtual void SetIntegrationPathes()=0;
	
	//virtual void FillStorages(dcx _Q, dcx M_in, dcx M_out)=0;
	
};

