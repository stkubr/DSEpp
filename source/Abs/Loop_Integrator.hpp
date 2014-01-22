#pragma once

class C_Loop_Integrator: public C_AbstractClass{
	
};

class C_1_Loop_Int: public C_Loop_Integrator{
	public:
	t_cmplxArray1D integrand_args;
	int k_col,Int_counter;
	C_Kinematics_1loop Momenta;
	C_Integrator_Line<t_cmplxMatrix,C_1_Loop_Int,double> * Integ_radial;
	C_Integrator_Line<t_cmplxMatrix,C_1_Loop_Int,double> * Integ_angle_cheb;
	C_Integrator_Line<t_cmplxMatrix,C_1_Loop_Int,double> * Integ_angle_Y;
	
	bool flag_sigma;
	
	C_1_Loop_Int(){
		integrand_args.resize(3);
	}
	
	//______________________________________________________________________
	// Multidimensional integration routine
	t_cmplxMatrix quad3d(){
		//integrand=func;
		k_col=0;
		Int_counter=0;
		return Integ_radial->getResult(&C_1_Loop_Int::f1,this);
	}
	t_cmplxMatrix f1(double k){
		k_col++;
		integrand_args[0]=(k);
		return Integ_angle_cheb->getResult(&C_1_Loop_Int::f2,this);
	}
	t_cmplxMatrix f2(double z){
		integrand_args[1]=z;
		flag_sigma=true;
		return Integ_angle_Y->getResult(&C_1_Loop_Int::f3,this);
	}
	t_cmplxMatrix f3(double y){
		integrand_args[2]=y;
		return (this->Integrand)(integrand_args);
	}
	
	virtual t_cmplxMatrix Integrand(t_cmplxArray1D args) {t_cmplxMatrix dummy; std::cout << "Error virtual call" << std::endl; StopLine(); return dummy;}
};

class C_2_Loop_Int: public C_Loop_Integrator{
	public:
	C_DedicMem_2_Loop_Int * Memory_Loop;
	t_dArray1D integrand_args;
	//dArray3D PointsAndWieghts;
	t_dArray2D zz_vector,w_vector;
	int Int_counter_total;
	int k_col_1,Int_counter_1;
	int k_col_2,Int_counter_2;
	C_Kinematics_1loop Momenta;
	C_Integrator_Line<t_cmplxMatrix,C_2_Loop_Int,double> * Integ_radial_1;
	C_Integrator_Line<t_cmplxMatrix,C_2_Loop_Int,double> * Integ_angle_cheb_1;
	C_Integrator_Line<t_cmplxMatrix,C_2_Loop_Int,double> * Integ_angle_Y_1;
	
	C_Integrator_Line<t_cmplxMatrix,C_2_Loop_Int,double> * Integ_radial_2;
	C_Integrator_Line<t_cmplxMatrix,C_2_Loop_Int,double> * Integ_angle_cheb_2;
	C_Integrator_Line<t_cmplxMatrix,C_2_Loop_Int,double> * Integ_angle_Y_2;
	C_Integrator_Line<t_cmplxMatrix,C_2_Loop_Int,double> * Integ_angle_Phi_2;
	
	bool flag_sigma=true;
	
	C_2_Loop_Int(){
		integrand_args.resize(7);
		zz_vector.resize(7);
		w_vector.resize(7);
		Memory_Loop=(C_DedicMem_2_Loop_Int*)DedicMemFactory_2_Loop_Int->CreateMemory();
	}
		
	//______________________________________________________________________
	// Multidimensional integration routine
	t_cmplxMatrix quad7d(){
		//integrand=func;
		Int_counter_total=0;
		k_col_1=0;
		Int_counter_1=0;
		return Integ_radial_1->getResult(&C_2_Loop_Int::f1,this);
	}
	t_cmplxMatrix f1(double k){
		k_col_1++;
		integrand_args[0]=(k);
		return Integ_angle_cheb_1->getResult(&C_2_Loop_Int::f2,this);
	}
	t_cmplxMatrix f2(double z){
		integrand_args[1]=z;
		return Integ_angle_Y_1->getResult(&C_2_Loop_Int::f3,this);
	}
	t_cmplxMatrix f3(double y){
		integrand_args[2]=y;
		
		return quad4d();
	}
	
	//______________________________________________________________________
	// Multidimensional integration routine
	t_cmplxMatrix quad4d(){
		//integrand=func;
		k_col_2=0;
		Int_counter_2=0;
		Int_counter_1++;
		return Integ_radial_2->getResult(&C_2_Loop_Int::s1,this);
	}
	t_cmplxMatrix s1(double k){
		k_col_2++;
		integrand_args[3]=(k);
		return Integ_angle_cheb_2->getResult(&C_2_Loop_Int::s2,this);
	}
	t_cmplxMatrix s2(double z){
		integrand_args[4]=z;
		return Integ_angle_Y_2->getResult(&C_2_Loop_Int::s3,this);
	}
	t_cmplxMatrix s3(double y){
		integrand_args[5]=y;
		return Integ_angle_Phi_2->getResult(&C_2_Loop_Int::s4,this);
	}
	t_cmplxMatrix s4(double phi){
		integrand_args[6]=phi;
		Int_counter_2++;
		return (this->Integrand)(integrand_args,Int_counter_total);
	}
	
	void SetIntegrationVectors(int NumPointsTotal){
		std::cout << NumPointsTotal << std::endl;
		Memory_Loop->ResizePointsAndWieghts(NumPointsTotal,7);
		Int_counter_total=0;
		Int_counter_1=0;
		std::vector<int> counter(7);
		for (counter[0] = 1; counter[0] < zz_vector[0].size(); counter[0]++){
			for (counter[1] = 1; counter[1] < zz_vector[1].size(); counter[1]++){
				for (counter[2] = 1; counter[2] < zz_vector[2].size(); counter[2]++){
					Int_counter_2=0;
					for (counter[3] = 1; counter[3]< zz_vector[3].size(); counter[3]++){
						for (counter[4] = 1; counter[4] < zz_vector[4].size(); counter[4]++){
							for (counter[5] = 1; counter[5] < zz_vector[5].size(); counter[5]++){
								for (counter[6] = 1; counter[6] < zz_vector[6].size(); counter[6]++){
									for (int i = 0; i < counter.size(); i++){
										Memory_Loop->PointsAndWieghts[0][Int_counter_total][i]=zz_vector[i][counter[i]];
										Memory_Loop->PointsAndWieghts[1][Int_counter_total][i]=w_vector[i][counter[i]];
										Memory_Loop->Counters[0][Int_counter_total]=Int_counter_1;
										Memory_Loop->Counters[1][Int_counter_total]=Int_counter_2;
										Memory_Loop->Counters[2][Int_counter_total]=Int_counter_total;
										//std::cout << Int_counter_total << "  " << Int_counter_1 << "  " << Int_counter_2 << std::endl;
									}
									Int_counter_2++;
									Int_counter_total++;
								}
							}
						}
					}
					Int_counter_1++;
				}
			}
		}
		
		std::cout << Int_counter_total << "  " << Int_counter_1 << "  " << Int_counter_2 << std::endl;
		std::cout << Memory_Loop->PointsAndWieghts[0][0].size() << std::endl;
	}
	
	t_cmplxMatrix quad7d_vector(){
		t_cmplxMatrix result(1,1);
		/*result=(this->Integrand)(Memory_Loop->PointsAndWieghts[0][0]);
		result.Zero();*/
		Int_counter_1=0;
		Int_counter_2=0;
		Int_counter_total=0;
#pragma omp parallel //num_threads(1)
{//start parallel region
		C_2_Loop_Int * Parallel_Copy;
		Parallel_Copy=MakeCopy();
		#pragma omp for
		for (int i = 0; i < Memory_Loop->PointsAndWieghts[0].size(); i++){
			double Weights=1.0;
			for (int j = 0; j < Memory_Loop->PointsAndWieghts[1][0].size(); j++){
				double temp=(Memory_Loop->PointsAndWieghts[1][i][j]);
				Weights*=temp;
			}
			
			t_cmplxMatrix temp_matrix;
			temp_matrix=(Parallel_Copy->Integrand)(Memory_Loop->PointsAndWieghts[0][i],i)*Weights;
			
			#pragma omp critical
			{
			result+=temp_matrix;
			}
		}
		delete Parallel_Copy;
}//end parallel region
		return result;
	}
	
	virtual C_2_Loop_Int * MakeCopy()
	{
		return new C_2_Loop_Int(*this);
	}
	
	virtual t_cmplxMatrix Integrand(t_dArray1D args, int ctr) {t_cmplxMatrix dummy; std::cout << "Error virtual call" << std::endl; StopLine(); return dummy;}
};
