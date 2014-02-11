#pragma once
#include "../AbsVertex.h"

class C_BSE_Base: public C_AbsVertex{
	protected:
	// 4-momenta vectors

	
	// Basic Dirac structure
	t_cmplxDirac S_p,S_m,Y_T;
	
	//Side Objects
	C_Propagator * Parton_P;
	C_Propagator * Parton_M;
	
	// Constructor
	C_BSE_Base()
	{
		Memory=(C_DedicMem_BSA*)DedicMemFactory_BSA->CreateMemory();
	}
	
	public:
	C_DedicMem_BSA * Memory;
	void LinkToKernel(C_AbstractKernel * _K){
		Kernel=_K;
	}
};

class C_BSE_Hadron_parameters{
	public:
	int NumRadial,NumCheb_nod1,NumCheb_nod2,Cheb_order,NumAngleY;
	int NumRadial_Contour, NumCheb_Contour, NumAngleY_Contour; 
	int NumRadial_Matrix, NumCheb_Matrix, NumAngleY_Matrix;
	double LimUk,LimDk,zetta_part;
	bool OffShell;
	
	void Print(std::ostream &__os = std::cout) const {
		__os << 
		"LimDk" << " - " << LimDk << std::endl <<
		"LimUk" << " - " << LimUk << std::endl <<
		"zetta_part" << " - " << zetta_part << std::endl <<
		"OffShell" << " - " << OffShell << std::endl <<
		
		"NumRadial" << " - " << NumRadial << std::endl <<
		"Cheb_order" << " - " << Cheb_order << std::endl <<
		"NumCheb_nod1" << " - " << NumCheb_nod1 << std::endl <<
		"NumCheb_nod2" << " - " << NumCheb_nod2 << std::endl <<
		"NumAngleY" << " - " << NumAngleY << std::endl <<
		
		"NumRadial_Contour" << " - " << NumRadial_Contour << std::endl <<
		"NumCheb_Contour" << " - " << NumCheb_Contour << std::endl <<
		"NumAngleY_Contour" << " - " << NumAngleY_Contour << std::endl <<
		
		"NumRadial_Matrix" << " - " << NumRadial_Matrix << std::endl <<
		"NumCheb_Matrix" << " - " << NumCheb_Matrix << std::endl <<
		"NumAngleY_Matrix" << " - " << NumAngleY_Matrix << std::endl <<
		
		std::endl;	
	}
	
	void setParams(std::ifstream * _ParamList){
			string line;
			if ((*_ParamList).is_open())
			{
				while ( (*_ParamList).good() )
				{
					(*_ParamList) >> line >> LimDk;
					(*_ParamList) >> line >> LimUk;
					(*_ParamList) >> line >> zetta_part;
					(*_ParamList) >> line >> OffShell;
					
					(*_ParamList) >> line >> NumRadial;
					(*_ParamList) >> line >> Cheb_order;
					(*_ParamList) >> line >> NumCheb_nod1;
					(*_ParamList) >> line >> NumCheb_nod2;
					(*_ParamList) >> line >> NumAngleY;
					
					(*_ParamList) >> line >> NumRadial_Contour;
					(*_ParamList) >> line >> NumCheb_Contour;
					(*_ParamList) >> line >> NumAngleY_Contour;
					
					(*_ParamList) >> line >> NumRadial_Matrix;
					(*_ParamList) >> line >> NumCheb_Matrix;
					(*_ParamList) >> line >> NumAngleY_Matrix;
					
				}
			}
			else std::cout << "Cant open file!" << std::endl;
		}
};

ostream& operator<<(ostream &__os, const C_BSE_Hadron_parameters &__params){
		__params.Print(__os);
		return __os;
	}


class C_BSE_Hadron_Base: public C_BSE_Base, public C_1_Loop_Int{
	protected:
	int /*Int_counter,*/Complex_Int_counter,index_zp,index_p;
	
	std::vector<t_cmplxDirac> Amplitudes,Projectors,WaveFunctions;
	t_cmplxDirac FullWaveFunction,FullAmplitude;
	
	ifstream * ParamList;
	const char * SaveBSEPath;
	Quark_ID Parton_ID_P,Parton_ID_M;
	
	C_BSE_Hadron_parameters params;
	//C_Kinematics_1loop Momenta;
	t_cmplxArray1D /*integrand_args,*/U_amp,WeightCoeff;
	t_dArray1D zz_rad, zz_cheb, zz_angleY , zz_cauchy, w_rad, w_cheb, w_angleY, w_cauchy, proj_amp;

	/*C_Integrator_Line<dcxMatrix,C_BSE_Hadron_Base,double> * Integ_radial;
	C_Integrator_Line<dcxMatrix,C_BSE_Hadron_Base,double> * Integ_angle_cheb;
	C_Integrator_Line<dcxMatrix,C_BSE_Hadron_Base,double> * Integ_angle_Y;*/
	C_Integrator_Cauchy<t_cmplxArray1D,t_cmplxArray3D,t_cmplx> * Integ_cauchy_long;
	
	t_cmplxMatrix Zero,dataAmp,AMP,BUFFER_AMP,BUFFER_F_ex,BUFFER_dataAmp_ex,NextDir;
	
	t_cmplxMatrix (C_BSE_Hadron_Base::*integrand) (t_cmplxArray1D);
	void (C_BSE_Hadron_Base::*SetDressing_ref) (t_cmplx);
	t_cmplxMatrix (C_BSE_Hadron_Base::*GetBSA_ref) ();
	// counters
	//int k_col;
	
	//flags
	public:
	double NORM;
	int num_amplitudes;
	bool flag_off_shell,flag_amp_desciption,flag_precalculation;
	
	public:
	C_BSE_Hadron_Base(){
		ifstream ParamList("Parameters_files/Mesons/Scalar_symmetric.txt");
		params.setParams(&ParamList);
		std::cout << params << std::endl;
		flag_off_shell=false;
		flag_amp_desciption=false;
		flag_precalculation=false;
		NORM=1.0;
	}
	
	~C_BSE_Hadron_Base(){
	}
	
	
	void setPropagators(t_cmplxVector *K_plus, t_cmplxVector *K_minus)
	{
		t_cmplxArray1D quark_temp_sigma;
		quark_temp_sigma=Parton_P->getPropAt((*K_plus)*(*K_plus));
		S_p=(-1.0*ii*((*K_plus)*Z)*quark_temp_sigma[3] + SHIFT*I*quark_temp_sigma[4]);
		
		quark_temp_sigma=Parton_M->getPropAt((*K_minus)*(*K_minus));
		S_m=(-1.0*ii*((*K_minus)*Z)*(quark_temp_sigma[3]) + SHIFT*I*(quark_temp_sigma[4]));
	}
	
	//virtual void SetProjectors(dcxVector _p, dcxVector _P){std::cout << "Error - virtual call" << std::endl; StopLine();};
	//virtual void SetAmplitudes(dcxVector _k, dcxVector _P){std::cout << "Error - virtual call" << std::endl; StopLine();};
	virtual void SetDiracStructures(t_cmplxVector _k, t_cmplxVector _P, std::vector<t_cmplxDirac> * DiracStructure){std::cout << "Error - virtual call" << std::endl; StopLine();};
	virtual void SetWaveFunctions(){std::cout << "Error - virtual call" << std::endl; StopLine();};
	virtual t_cmplxMatrix GetBSA(){t_cmplxMatrix dummy; std::cout << "Error - virtual call" << std::endl; StopLine(); return dummy;};
	virtual t_cmplxMatrix GetBSA_matrix(){t_cmplxMatrix dummy; std::cout << "Error - virtual call" << std::endl; StopLine(); return dummy;};
	virtual t_cmplxMatrix GetBSA_norm(){t_cmplxMatrix dummy; std::cout << "Error - virtual call" << std::endl; StopLine(); return dummy;};
	
	void SetWeightCoeff(){
		for (int i = 0; i < num_amplitudes; i++) { WeightCoeff[i]=(Projectors[i]|Projectors[i]).Tr(); }
	}
	
	void OrthogonalityCheck(){
		for (int i = 0; i < num_amplitudes; i++) { 
			for (int j = 0; j < num_amplitudes; j++) { 
				double res;
				t_cmplx dummy;
				dummy=(Projectors[i]|Projectors[j]).Tr();
				res=real(dummy);
				std::cout << i+1 << "  " << j+1 << "  " << res << "  || ";
			}	
			std::cout << std::endl << std::endl;
		}
		cin.get();
	}
	
	void Initialization(){
		Int_counter=0;
		Integ_radial=C_Integrator_Line<t_cmplxMatrix,C_1_Loop_Int,double>::createIntegrator(params.NumRadial, params.LimDk, params.LimUk, num_amplitudes, qgausleg_log_ID);
		Integ_angle_cheb=C_Integrator_Line<t_cmplxMatrix,C_1_Loop_Int,double>::createIntegrator(params.NumCheb_nod1, params.LimDk, params.LimUk, num_amplitudes, qgauscheb_ID);
		Integ_angle_Y=C_Integrator_Line<t_cmplxMatrix,C_1_Loop_Int,double>::createIntegrator(params.NumAngleY, -1.0, 1.0, num_amplitudes, qgausleg_sym_ID);
		Integ_cauchy_long=C_Integrator_Cauchy<t_cmplxArray1D,t_cmplxArray3D,t_cmplx>::createIntegrator(params.NumRadial_Contour, sqrt(params.LimDk), sqrt(params.LimUk), num_amplitudes, qcauchyleg_lin_ID);
		integrand_args.resize(3);
		
		//Parton_P=QuarkFactory->Create(&Parton_ID_P);
		//Parton_M=QuarkFactory->Create(&Parton_ID_M);
		//Parton_M=Parton_P;
		
		//Kernel->setZ2DressingFactor(Parton_P->getDressingFactor());
		
		ResizeALL();
		
		setInitialAMP();
	}
	
	void LinkToPartons(C_Propagator * _P,C_Propagator * _M){
		Parton_P=_P;
		Parton_M=_M;
	}
	
	void ResizeALL(){

		Integ_radial->getNodes(&zz_rad,&w_rad);
		Integ_angle_cheb->getNodes(&zz_cheb,&w_cheb);
		Integ_angle_Y->getNodes(&zz_angleY,&w_angleY);
		Integ_cauchy_long->getNodes(&zz_cauchy,&w_cauchy);
		
		Zero.Resize(num_amplitudes,1);
		dataAmp.Resize(num_amplitudes,1);
		BUFFER_AMP.Resize(params.NumRadial+1,params.Cheb_order*num_amplitudes+1);
		BUFFER_F_ex.Resize(params.NumRadial+1,(params.NumCheb_nod2)*num_amplitudes+1);
		BUFFER_dataAmp_ex.Resize(num_amplitudes,params.NumCheb_nod2+1);
		AMP.Resize(params.NumRadial+1,params.Cheb_order*num_amplitudes+1);
		NextDir.Resize(params.NumRadial+1,params.Cheb_order*num_amplitudes+1);
		U_amp.resize(num_amplitudes);
		proj_amp.resize(params.Cheb_order+1);
	}

	
	void setInitialAMP(){
		for (int i = 1; i <= params.NumRadial; i++){
			for (int j = 1; j <=params.Cheb_order*num_amplitudes ; j++){
				AMP(i,j)=1.0;
			}	
		}
	}
	
	void writeBSA()
	{
		ofstream BSA_Stream;
		BSA_Stream.open (SaveBSEPath);
		for (int j = 1; j <= params.NumRadial ; j++)
		{
			BSA_Stream << AMP(j,0) << '\t' ;
			for (int k = 1; k <= num_amplitudes; k++)
			{
				for (int i = 1; i <= params.Cheb_order ; i++)
				{
					BSA_Stream << AMP(j,i + params.Cheb_order*(k-1)) << '\t' ;
				}
			}
			BSA_Stream << std::endl;
		}
		BSA_Stream.close();
	}	
	
	void set_BSA_on_grid()
	{
		ifstream BSA_Stream(SaveBSEPath);
		if (BSA_Stream.is_open())
		{
			for (int i=1;i <= params.NumRadial;i++)
			{
				for (int j = 0; j <= num_amplitudes ; j++)
				{
					BSA_Stream >> AMP(i,j);
				}
			}
			std::cout << "Initial AMP has been set!" << std::endl; 
		}
		else {std::cout << "Cant open file!" << std::endl; setInitialAMP();}
	}
	
	void setProj(int k){
			for (int i = 1; i <= params.Cheb_order; i++){
				proj_amp[i]=0.0;
			}
			proj_amp[k]=1.0;
	}
	
	t_cmplx U_ex(t_cmplx z_ex)
	{
		t_cmplx sum_ex; sum_ex=0;
		for (int i = 1; i <= params.Cheb_order ; i++)
		{
			sum_ex+=Cheb_polinoms(z_ex,(i-1)*2)*proj_amp[i];
		}
		return sum_ex;
	}
	
	virtual C_BSE_Hadron_Base * MakeCopy()
	{
		return new C_BSE_Hadron_Base(*this);
	}
		
	t_cmplxMatrix IntegAngleY(){
		return Integ_angle_Y->getResult(&C_BSE_Hadron_Base::f3,this);
	}
	
	void PreCalculation(){
		Memory->resizeAmpStorage(num_amplitudes,params.NumRadial*params.NumCheb_nod1*params.NumAngleY);
		int int_ctr=0;
		for (int i = 1; i < zz_rad.size(); i++){
			for (int j = 1; j < zz_cheb.size(); j++){
				for (int k = 1; k < zz_angleY.size(); k++){
					SetIntMomenta(sqrt(zz_rad[i]),zz_angleY[k],zz_cheb[j]);
					SetWaveFunctions();
					for (int amp_ctr = 0; amp_ctr < num_amplitudes; amp_ctr++){ Memory->setAmpStorage(amp_ctr, int_ctr, &WaveFunctions[amp_ctr]); }
					int_ctr++;	
				}
			}
		}
		flag_precalculation=true;
		//cin.get();
	}
	
	void SetIntMomenta(t_cmplx x, t_cmplx y, t_cmplx z){
		Momenta.SetVectors_k(params.zetta_part,x,y,z);
		Momenta.SetVectors_q();
		Momenta.SetVestors_k_for_S(params.zetta_part,Momenta.k);
	}
	
	void CalcEVMatrix(Eigen::ComplexEigenSolver<Eigen::MatrixXcf> * ces){
		Memory->ResizeEVMatrix(params.NumRadial,params.NumCheb_nod1,num_amplitudes,1);
		//std::cout << "Memory->EVMatrix call" << std::endl;
		#pragma omp parallel 
		{//start of pragma	
			C_BSE_Hadron_Base * bsa_copy_omp;
			bsa_copy_omp=MakeCopy();
			t_cmplxMatrix Temp_Matrix;
			double _p2,_k2,_z,_y,_zp;
			double _w_zp,_w_k2,_w_z;
			bsa_copy_omp->k_col=0;
			bsa_copy_omp->Int_counter=0;
		#pragma omp for 
			for (int p_ctr = 1; p_ctr < zz_rad.size() ; p_ctr++){
				bsa_copy_omp->index_zp=0;
				bsa_copy_omp->index_p=p_ctr;
				_p2=zz_rad[p_ctr];
				for (int zp_ctr = 1; zp_ctr < zz_cheb.size() ; zp_ctr++){
					bsa_copy_omp->index_zp=zp_ctr;
					_zp=zz_cheb[zp_ctr];
					bsa_copy_omp->Momenta.SetVectors_p(_zp,sqrt(_p2));
					bsa_copy_omp->SetDiracStructures(bsa_copy_omp->Momenta.p,bsa_copy_omp->Momenta.P,&bsa_copy_omp->Projectors);
					bsa_copy_omp->SetWeightCoeff();
					//FullWaveFunction=Projectors[0];
					//_w_zp=w10_ch_ex[zp_ctr];
					for (int k_ctr = 1; k_ctr < zz_rad.size(); k_ctr++){
						//k_ctr = 1;
						_k2=zz_rad[k_ctr];	
						_w_k2=w_rad[k_ctr];
						bsa_copy_omp->integrand_args[0]=_k2;
						bsa_copy_omp->flag_sigma=false;
						for (int z_ctr = 1; z_ctr < zz_cheb.size() ; z_ctr++){
							_z=zz_cheb[z_ctr];
							_w_z=w_cheb[z_ctr];
							bsa_copy_omp->integrand_args[1]=_z;
							//Temp_Matrix=bsa_copy_omp->SetMatrixAtPoint(_P,_p2,_zp,_k2,_z);
							Temp_Matrix=bsa_copy_omp->IntegAngleY();
							for (int p_amp_ctr = 0; p_amp_ctr < num_amplitudes ; p_amp_ctr++){
								for (int k_amp_ctr = 0; k_amp_ctr < num_amplitudes ; k_amp_ctr++){
									
									//std::cout << k_ctr << "  " << z_ctr << "  " << z_ctr-1 + (k_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*k_amp_ctr << std::endl;
									//std::cout << k_ctr << "  " << z_ctr << "  " << zp_ctr + (p_ctr-1)*(NumRadialPoints-z_ctr*(NumRadialPoints-1))-1+ NumRadialPoints*(num_cheb_nod1)*p_amp_ctr << "  " << z_ctr + (k_ctr-1)*(NumRadialPoints-z_ctr*(NumRadialPoints-1)) + NumRadialPoints*(num_cheb_nod1)*k_amp_ctr-1 << std::endl;
									Memory->EVMatrix(zp_ctr-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*p_amp_ctr, z_ctr-1 + (k_ctr-1)*(params.NumCheb_nod1) + params.NumRadial*(params.NumCheb_nod1)*k_amp_ctr)=pi/2.0*_w_k2*_w_z*Temp_Matrix(p_amp_ctr,k_amp_ctr);
								}
							}
							//std::cout << Temp_Matrix(0,0) << std::endl;
							//cin.get();
						}
						bsa_copy_omp->k_col++;
					}
					bsa_copy_omp->k_col=0;
					bsa_copy_omp->Int_counter=0;
				}
			}
			delete bsa_copy_omp;
		}// end of pragma
		//std::cout << "EVMatrix is full. EigenValues computation engaged..." << std::endl;
		ces->compute(Memory->EVMatrix);
	}
	
	t_cmplxArray2D SetEVMatrix(t_cmplx _P){
		Momenta.SetVector_P(_P);
		PreCalculation();
		t_dArray1D zz_rad_temp,zz_cheb_temp;
		SetDressing_ref=&C_BSE_Hadron_Base::SetDressing_normal;
		GetBSA_ref=&C_BSE_Hadron_Base::GetBSA_matrix;
		
		setInitialAMP();
		
		Eigen::ComplexEigenSolver<Eigen::MatrixXcf> ces;
		
		Memory->ResizeEVMatrix(params.NumRadial,params.NumCheb_nod1,num_amplitudes,1);
		//std::cout << "Memory->EVMatrix call" << std::endl;
		#pragma omp parallel 
		{//start of pragma	
			C_BSE_Hadron_Base * bsa_copy_omp;
			bsa_copy_omp=MakeCopy();
			t_cmplxMatrix Temp_Matrix;
			double _p2,_k2,_z,_y,_zp;
			double _w_zp,_w_k2,_w_z;
			bsa_copy_omp->k_col=0;
			bsa_copy_omp->Int_counter=0;
		#pragma omp for 
			for (int p_ctr = 1; p_ctr < zz_rad.size() ; p_ctr++){
				bsa_copy_omp->index_zp=0;
				bsa_copy_omp->index_p=p_ctr;
				_p2=zz_rad[p_ctr];
				for (int zp_ctr = 1; zp_ctr < zz_cheb.size() ; zp_ctr++){
					bsa_copy_omp->index_zp=zp_ctr;
					_zp=zz_cheb[zp_ctr];
					bsa_copy_omp->Momenta.SetVectors_p(_zp,sqrt(_p2));
					bsa_copy_omp->SetDiracStructures(bsa_copy_omp->Momenta.p,bsa_copy_omp->Momenta.P,&bsa_copy_omp->Projectors);
					bsa_copy_omp->SetWeightCoeff();
					//FullWaveFunction=Projectors[0];
					//_w_zp=w10_ch_ex[zp_ctr];
					for (int k_ctr = 1; k_ctr < zz_rad.size(); k_ctr++){
						//k_ctr = 1;
						_k2=zz_rad[k_ctr];	
						_w_k2=w_rad[k_ctr];
						bsa_copy_omp->integrand_args[0]=_k2;
						bsa_copy_omp->flag_sigma=false;
						for (int z_ctr = 1; z_ctr < zz_cheb.size() ; z_ctr++){
							_z=zz_cheb[z_ctr];
							_w_z=w_cheb[z_ctr];
							bsa_copy_omp->integrand_args[1]=_z;
							//Temp_Matrix=bsa_copy_omp->SetMatrixAtPoint(_P,_p2,_zp,_k2,_z);
							Temp_Matrix=bsa_copy_omp->IntegAngleY();
							for (int p_amp_ctr = 0; p_amp_ctr < num_amplitudes ; p_amp_ctr++){
								for (int k_amp_ctr = 0; k_amp_ctr < num_amplitudes ; k_amp_ctr++){
									
									//std::cout << k_ctr << "  " << z_ctr << "  " << z_ctr-1 + (k_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*k_amp_ctr << std::endl;
									//std::cout << k_ctr << "  " << z_ctr << "  " << zp_ctr + (p_ctr-1)*(NumRadialPoints-z_ctr*(NumRadialPoints-1))-1+ NumRadialPoints*(num_cheb_nod1)*p_amp_ctr << "  " << z_ctr + (k_ctr-1)*(NumRadialPoints-z_ctr*(NumRadialPoints-1)) + NumRadialPoints*(num_cheb_nod1)*k_amp_ctr-1 << std::endl;
									Memory->EVMatrix(zp_ctr-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*p_amp_ctr, z_ctr-1 + (k_ctr-1)*(params.NumCheb_nod1) + params.NumRadial*(params.NumCheb_nod1)*k_amp_ctr)=pi/2.0*_w_k2*_w_z*Temp_Matrix(p_amp_ctr,k_amp_ctr);
								}
							}
							//std::cout << Temp_Matrix(0,0) << std::endl;
							//cin.get();
						}
						bsa_copy_omp->k_col++;
					}
					bsa_copy_omp->k_col=0;
					bsa_copy_omp->Int_counter=0;
				}
			}
			delete bsa_copy_omp;
		}// end of pragma
		//std::cout << "EVMatrix is full. EigenValues computation engaged..." << std::endl;
		ces.compute(Memory->EVMatrix);
		flag_precalculation=false;
		//CalcEVMatrix(&ces);
		Eigen::VectorXcf EV=ces.eigenvalues();
	
	
		auto Parity = [&](int num_state) -> t_cmplx { 
			t_cmplx parity=0.0;
			for (int p_ctr = 1; p_ctr < 2 ; p_ctr++){
				for (int p_amp_ctr = 0; p_amp_ctr < 1 ; p_amp_ctr++){
					for (int zp_ctr = 1; zp_ctr < zz_cheb.size() ; zp_ctr++){
						t_cmplx summ,diff;
						
						summ = ces.eigenvectors().col(num_state)(zp_ctr-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*p_amp_ctr) +
						ces.eigenvectors().col(num_state)(zz_cheb.size() - zp_ctr-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*p_amp_ctr);
						
						diff = ces.eigenvectors().col(num_state)(zp_ctr-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*p_amp_ctr) -
						ces.eigenvectors().col(num_state)(zz_cheb.size() - zp_ctr-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*p_amp_ctr);

						if(fabs(real(summ))<=fabs(real(diff))) {parity = -1.0;}
						else{ parity = 1.0; }
					}
				}
			}
			return parity;
		};
		
		//std::cout << "EigenValues computation is done. The eigenvalues of EVMatrix are obtained." << std::endl;
		int i= EV.size()-1;
		t_cmplxArray2D Dominant_EV_and_parity(2);
		while( i > EV.size()-10)
		{
			//std::cout << EV[i] << "  " << i << std::endl;
			Dominant_EV_and_parity[0].push_back(EV[i]);
			Dominant_EV_and_parity[1].push_back(Parity(i));
			//std::cout << i << "  " << EV[i] << "  " << Parity(i)<< std::endl;
			i--;
		}
		//std::cout << ces.eigenvectors().col(EV.size()-1) << std::endl;
		return Dominant_EV_and_parity;
	}
	
	void DrawBSA_matrix(t_cmplx _P, int _state, int amp_num){
		int num_state;
		Momenta.SetVector_P(_P);
		PreCalculation();
		t_dArray1D zz_rad_temp,zz_cheb_temp;
		SetDressing_ref=&C_BSE_Hadron_Base::SetDressing_normal;
		GetBSA_ref=&C_BSE_Hadron_Base::GetBSA_matrix;
		
		flag_precalculation=false;
		Eigen::ComplexEigenSolver<Eigen::MatrixXcf> ces;
		//ces.compute(Memory->EVMatrix);
		CalcEVMatrix(&ces);
		Eigen::VectorXcf EV=ces.eigenvalues();
		num_state = EV.size()- _state;
		
		ofstream DrawBSA_matrix;
		DrawBSA_matrix.open ("Data_files/BSEs_matrix.dat");
		
		int i=EV.size()-1;
		
		while(i > EV.size()-6)
		{
			//std::cout << EV[i] << "  " << i << std::endl;
			i--;
		}
		
		t_cmplxArray1D result(zz_rad.size());
		for (int p_ctr = 1; p_ctr <= params.NumRadial ; p_ctr++){
			result[p_ctr-1]=0;
			for (int j=1;j<=params.NumCheb_nod1;j++) 
			{
				t_cmplx z_v,element;
				element=ces.eigenvectors().col(num_state)(j-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*amp_num);
				z_v=zz_cheb[j];	
				result[p_ctr-1] += w_cheb[j]*Cheb_polinoms(z_v,0)*real(element);
				//std::cout << w_cheb[j] << "  " << U_ex(z_v) << "  " << real(element) << std::endl;
			}
			DrawBSA_matrix << zz_rad[p_ctr] << "  " << real(result[p_ctr-1]) << std::endl;
			//std::cout << zz_rad[p_ctr] << "  " << result[p_ctr-1] << std::endl;
			//cin.get();	
		}
		DrawBSA_matrix.close();
	}
		
	virtual t_cmplxMatrix setResultBSA(){t_cmplxMatrix dummy; std::cout << "Error virtual call" << std::endl; StopLine(); return dummy;}

	void SetDressing_normal(t_cmplx z){
		if (flag_sigma == true){ 
			flag_sigma = false;
			for (int j = 0; j < num_amplitudes ; j++){	
				U_amp[j]=0.0;    
			}
			for (int i = 1; i <= params.Cheb_order ; i++){
				t_cmplx temp_in;
				temp_in=Cheb_polinoms(z,(i-1)*2);
				for (int j = 0; j < num_amplitudes ; j++){
					U_amp[j]+=temp_in*(AMP(k_col,i+j*params.Cheb_order));    
				}
			}
	    }
	}
	
	void SetDressing_shifted(t_cmplx z){
		Momenta.ShiftMomenta(params.zetta_part);
		for (int j = 0; j < num_amplitudes ; j++){	
			U_amp[j]=0.0;    
		}
		for (int i = 1; i <= params.Cheb_order ; i++){
			t_cmplx temp_in;
			temp_in=Cheb_polinoms(z,(i-1)*2);
			for (int j = 0; j < num_amplitudes ; j++){
				U_amp[j]+=temp_in*(Memory->CauchyGrid[j+1][index_p][Complex_Int_counter]); 
			}
		}
	}
	
	t_cmplxMatrix Integrand(t_cmplxArray1D args){
		t_cmplxMatrix result;
		t_cmplx Kinematic_factor;
		t_cmplx x,y,z;
		x=sqrt(args[0]);
		z=args[1];
		y=args[2];
		SetIntMomenta(x, y, z);
		
		/*#pragma omp master
		{
			std::cout << flag_sigma << std::endl;

		}
		cin.get();*/
		
		(this->*SetDressing_ref)(z);
		
		Kinematic_factor=-1.0/(8.0*pi*pi*pi*pi)*(args[0]);//(1.0 + args[0]/1000.0);//(1.0 + real(Momenta.P2)/1000.0)/(1.0 + real(Momenta.p2)/1000.0);
		result=Kinematic_factor*setResultBSA();
					
		Int_counter++;
		Complex_Int_counter++;
		return result;
	}
	
	t_cmplxMatrix CalcBSA(t_cmplx _p, t_cmplx _P, int proj){
		t_cmplx Bare_vertex;
		int z_ex_counter=0;
		t_cmplxMatrix result(num_amplitudes,1);
		Momenta.SetVector_P(_P);
		Complex_Int_counter=0;
		setProj(proj);
		if(proj==1){
			for (int j=1;j<=params.NumCheb_nod2;j++) 
			{
				t_cmplx z_v;
				z_v=zz_cheb[j];
				Momenta.SetVectors_p(z_v,_p);
				
				SetDiracStructures(Momenta.p,Momenta.P,&Projectors);
				SetWeightCoeff();
				FullWaveFunction=Projectors[0];
				
				#if ORTH_CHECK
				#pragma omp master
				{
					std::cout << "Orthogonality Check - " << std::endl;
					OrthogonalityCheck();
				}
				cin.get();
				#endif

				
				index_zp=z_ex_counter;
				z_ex_counter++;	
				
				//DebugLine("before");
				
				result=quad3d();

				for (int i = 0; i < num_amplitudes ; i++)
				{
					Bare_vertex=0.0;
					if (flag_off_shell && i==0) {Bare_vertex=2.0/pi*Parton_P->getDressingFactor();}
					BUFFER_dataAmp_ex(i,z_ex_counter)=result(i,0)+Bare_vertex;
				}
			}
		}
		
		for (int i = 0; i < num_amplitudes ; i++)
		{	
			result(i,0)=0.0;
			for (int j=1;j<=params.NumCheb_nod2;j++) 
			{
				t_cmplx z_v;
				z_v=zz_cheb[j];	
				result(i,0) += w_cheb[j]*BUFFER_dataAmp_ex(i,j)*U_ex(z_v);
				//std::cout << BUFFER_dataAmp_ex(0,j) << "  " << z_v << "  " << j << "  " << w_cheb[j] << "  " << U_ex(z_v) << std::endl;
			}	
		}
		return result;
	}
	
// 
//------------------------------------------------------------------
	void BSA_step(t_cmplx P)
	{
		for (int proj_cheb = 1; proj_cheb <=params.Cheb_order ; proj_cheb++)
		{	
#pragma omp parallel 
{//start of pragma		
				t_cmplx p2;
				C_BSE_Hadron_Base * bsa_copy_omp;
				bsa_copy_omp=MakeCopy();
				//bsa_copy_omp->ResetPointers();
				t_cmplxMatrix Temp_matrix(num_amplitudes,1);
				#pragma omp for	
				for (int i = 1; i <= params.NumRadial; i++)
				{
					if(proj_cheb>=2)
					{
						for (int ampl = 0; ampl < num_amplitudes ; ampl++)
						{
						//std::cout << "more Chebychevs !! ";
							for (int w = 1; w <=params.NumCheb_nod2; w++)
							{
							bsa_copy_omp->BUFFER_dataAmp_ex(ampl,w)=BUFFER_F_ex(i,w + params.NumCheb_nod2*(ampl))*1.0;
							//std::cout << bsa_copy_omp.F_buffer[w] << "  ";
							}
						//std::cout << std::endl;
						}
					}
							
					bsa_copy_omp->index_p=i;
					p2=zz_rad[i];
					BUFFER_AMP(i,0)=p2;
					//std::cout << "Hi, I`m Debug Line" << "   Before Integr!" << std::endl;
					
					Temp_matrix=bsa_copy_omp->CalcBSA(sqrt(p2),P,proj_cheb);
					
					//std::cout << "Hi, I`m Debug Line" << "   After Integr!" << std::endl;
					for (int ampl = 0; ampl < num_amplitudes ; ampl++)
					{
						
						BUFFER_AMP(i,proj_cheb+params.Cheb_order*(ampl))=Temp_matrix(ampl,0);
						//std::cout << BUFFER_AMP(i,proj_cheb+params.Cheb_order*(ampl)) << "  " << i << "  " << proj_cheb+params.Cheb_order*(ampl) << std::endl;
					}
					//std::cout << "Blah blah, I`m debug line!" << "  " << i << "  " << proj_cheb << "  " << BUFFER_AMP(i,0) << "  " << BUFFER_AMP(i,proj_cheb) << "  " << p2  << std::endl;
					if(proj_cheb==1)
					{
						//std::cout << " The first Chebychevs !! "<< std::endl;
						for (int ampl = 0; ampl < num_amplitudes ; ampl++)
						{
							for (int w = 1; w <= params.NumCheb_nod2 ; w++)
							{
							BUFFER_F_ex(i,w + params.NumCheb_nod2*(ampl))=bsa_copy_omp->BUFFER_dataAmp_ex(ampl,w);
							//std::cout << BUFFER_F_ex(i,w + num_cheb_nod2*(ampl)) << "  " << w << "  " << ampl << std::endl;
							}
						}
					}
				}
				delete bsa_copy_omp;
}//end of pragma
			}
		
		}
			
// 
//------------------------------------------------------------------
	t_cmplx DressBSA(t_cmplx P, int steps)
	{
		t_cmplx ff1,ff2,ff3,Lambda_EV;
		ff1=0;ff2=0;ff3=0;
		Momenta.SetVector_P(P);
		PreCalculation();
		SetDressing_ref=&C_BSE_Hadron_Base::SetDressing_normal;
		GetBSA_ref=&C_BSE_Hadron_Base::GetBSA;
		std::cout << "Dressing BSA... " << std::endl << std::endl;
		std::cout << fixed;
		std::cout << setprecision (NUM_PRECISION);
	
		//if(flag_load_AMP) 
		set_BSA_on_grid();
	
		for (int kk = 1; kk <= steps; kk++)
		{
			std::cout << "Step number - " << kk << std::endl;
			
			BSA_step(P);
			
			//std::cout << "Time now - " << (Get_Time() - Init_time)/NumThreadsUse << std::endl;
			
			for (int k = 1; k <=1; k++)
			{
				for (int i = 1; i <= params.NumRadial; i++)
				{
					ff1+=(BUFFER_AMP(i,k)*AMP(i,k));
					ff2+=(AMP(i,k)*AMP(i,k));
				}
			}
			
			Lambda_EV=(ff1/ff2);
			PrintLine('-');
			std::cout  << "Eigenvalue.v1 - " << ff1/ff2 << "  at P = " << "  " << P << std::endl;
			PrintLine('-');
			ff1=0;ff2=0;ff3=0;
			
			PreNormBSA();
			if(flag_amp_desciption) WriteAmplitudeDescyption();
			
			PrintLine('#');
		}
		//std::cout << "Time now - " << (Get_Time() - Init_time) << std::endl;
		//DrawBSA(P);

		writeBSA();

		flag_precalculation=false;
		return Lambda_EV;
	}
	
	void PreNormBSA()
	{
		for (int i = 1; i <= params.NumRadial; i++)
		{
			for (int j = 1; j <= params.Cheb_order*(num_amplitudes); j++)
			{
			AMP(i,0)=BUFFER_AMP(i,0);
			if(flag_off_shell)AMP(i,j)=BUFFER_AMP(i,j);//BUFFER_AMP(1,1);
			else AMP(i,j)=BUFFER_AMP(i,j)/(real(BUFFER_AMP(1,1)));
			
			}
		}
	}
	
	void WriteAmplitudeDescyption()
	{
		std::cout << fixed;
		std::cout << setprecision (NUM_PRECISION);
		std::cout << "Dirac Amplitudes" << std::endl;
		PrintSpaces((NUM_PRECISION+2)*2+3+2);
		for (int i = 0; i < num_amplitudes; i++)
		{
			std::cout << i;
			for (int j = 1; j <= params.Cheb_order; j++)
			{
				PrintSpaces((NUM_PRECISION+2)*2+3);
			}
			std::cout << "";
		}
		std::cout << std::endl;
		std::cout << "Chebyshev Projections" << std::endl << "    p^2";
		PrintSpaces((NUM_PRECISION+2)*2+3-5);
		for (int j = 0; j <num_amplitudes; j++)
		{
			for (int i = 1; i <=params.Cheb_order ; i++)
			{
			std::cout << (i-1)*2; 
			PrintSpaces((NUM_PRECISION+2)*2+3);
			}
		}
		std::cout << std::endl << AMP << std::endl;
	}
	
	void DrawBSA(t_cmplx _P){
		t_cmplxMatrix TempArray(num_amplitudes,1);
		ofstream temp_continuation;
		temp_continuation.open ("Data_files/BSEs.dat");
		for (int i = 1; i < zz_rad.size(); i++)
		{
			TempArray=CalcBSA(sqrt(zz_rad[i]),_P,1);
			temp_continuation << (zz_rad[i]);
			for (int amp = 0; amp < num_amplitudes; amp++) temp_continuation  << '	' <<  real(TempArray(amp,0));
			temp_continuation << std::endl;
		}
		temp_continuation.close();
	}
	
	
	double NormalizeBSA(t_cmplx _K){
		std::cout << std::endl;
		PrintLine('=');
		std::cout << std::endl;
		std::cout << "Calculation Norm Factor..." << std::endl;
		std::cout << std::endl;
		double BSE_Norm_factor;
		t_cmplx h=_K/5.0;
		t_cmplx res_int,N,Lambda_EV;
		flag_amp_desciption=false;
		std::vector <t_cmplx> deriv_step_norm(4);
		for (int i = 1; i <= 2; i++)
		{
			double i_count;
			i_count=i;
			std::cout << "Derivative at point " << _K+(i_count)*h << std::endl;
			deriv_step_norm[(i-1)]=DressBSA(_K+(i_count)*h,10);
			
			std::cout << "Derivative at point " << _K-(i_count)*h << std::endl;
			deriv_step_norm[(i-1)+2]=DressBSA(_K-(i_count)*h,10);
		}
		//flag_normalization=true;
		
		Lambda_EV=DressBSA(_K,10);
		
		SetDressing_ref=&C_BSE_Hadron_Base::SetDressing_normal;
		GetBSA_ref=&C_BSE_Hadron_Base::GetBSA_norm;
			
		res_int=quad3d()(0,0);
		//flag_normalization=false;
		N=(-deriv_step_norm[1] + 8.0*deriv_step_norm[0] - 8.0*deriv_step_norm[2] + deriv_step_norm[3])/(12.0*imag(h));
		std::cout << deriv_step_norm[0] << "  " << deriv_step_norm[1] << "  " << deriv_step_norm[2] << "  " << deriv_step_norm[3] << std::endl;
		std::cout << N << "  " << res_int << std::endl;
		//NORM=sqrt(1.0/fabs(real(N))/3.0/fabs(real(res_int))/Lambda_EV)/sqrt(2.0*imag(_K))*0.5;
		//NORM=sqrt(fabs(real(N))/fabs(real(res_int)))/sqrt(2.0*imag(K_v));
		
		BSE_Norm_factor=real(sqrt(Lambda_EV/fabs(real(N))/fabs(real(res_int)))*sqrt(2.0*imag(_K)));
		
		std::cout << "NormFactor = "<< BSE_Norm_factor << std::endl;
		
		return BSE_Norm_factor;
	}
	
	void SetBSAonPath(t_cmplxArray1D (*AmplitudePath),t_cmplxArray1D (*Path) ,t_cmplx Q)
	{
		std::cout << std::endl;
		std::cout << "On Path calculation..." << std::endl;
		flag_amp_desciption=true;
		DressBSA(Q,1);
		int num_points;
		double Time;
		num_points=(*Path).size();
		(*AmplitudePath).resize(num_points);
		t_cmplxMatrix Temp(num_amplitudes,1);
		std::cout << "Initialized" << std::endl;
		std::cout << "points on the path - " << num_points << std::endl;
		Time=Get_Time();
	#pragma omp parallel 
	{//start of pragma		
			C_BSE_Hadron_Base * bsa_copy_omp;
			bsa_copy_omp=MakeCopy();
			t_cmplxMatrix Temp_matrix(num_amplitudes,1);
			#pragma omp for 
			for (int i = 0; i < num_points; i++)
			{
				(*AmplitudePath)[i]=NORM*bsa_copy_omp->CalcBSA(sqrt((*Path)[i]),Q,1)(0,0);
			}
			delete bsa_copy_omp;
	}//end of pragma
		std::cout << " Time for continuation spent - " << (Get_Time() - Time)/8.0 << std::endl;
		std::cout << "On Path calculation finished." << std::endl;
	}
	
	
void SetBSAonPath(std::vector<t_cmplxMatrix> (*AmplitudePath),t_cmplxArray1D (*Path) ,t_cmplx Q, bool Interpolate_flag)
{
	std::cout << std::endl;
	std::cout << "On Path calculation..." << std::endl;
	flag_amp_desciption=false;
	if(!Interpolate_flag) DressBSA(Q, 12); else DressBSA(Q, 8);
	int num_points;
	double Time;
	num_points=(*Path).size();
	(*AmplitudePath).clear();
	(*AmplitudePath).resize(num_points);
	t_cmplxMatrix Temp(num_amplitudes+1,1);
	std::cout << "Initialized" << std::endl;
	std::cout << "points on the path - " << num_points << std::endl;
	if(!Interpolate_flag)
	{
		Time=Get_Time();
#pragma omp parallel 
{//start of pragma		
		C_BSE_Hadron_Base * bsa_copy_omp;
		bsa_copy_omp=MakeCopy();
		t_cmplxMatrix Temp_matrix(num_amplitudes,1);
		#pragma omp for 
		for (int i = 0; i < num_points; i++)
		{
			(*AmplitudePath)[i]=bsa_copy_omp->CalcBSA(sqrt((*Path)[i]),Q,1);
		}
		delete bsa_copy_omp;
}//end of pragma
	std::cout << " Time for continuation spent - " << (Get_Time() - Time)/8.0 << std::endl;
	} 
	else {
		t_cmplxArray2D Amplitude(num_amplitudes+1, t_cmplxArray1D(params.NumRadial));
		std::cout << "Enter Interpolation routine..." << std::endl;
		for (int i = 0; i <num_amplitudes+1 ; i++)
		{
			for (int j = 0; j < params.NumRadial; j++)
			{
				Amplitude[i][j]=AMP(j+1,i);
				//std::cout << i << "  " << j << "  " << Amplitude[i][j] << std::endl;
			}
			//cin.get();
		}
		std::cout << "Temp Amplitute matrix has been set." << std::endl;
		for (int j = 0; j < num_points; j++)
		{
			(*AmplitudePath)[j].Resize(num_amplitudes,1);
				//std::cout << i << "  " << j << "  " << Amplitude[i][j] << std::endl;
		}
		std::cout << "Passed storage has been resized." << std::endl;
		for (int i = 1; i < num_amplitudes+1; i++)
		{
			Interpolation::Linear<t_cmplx,t_cmplx> FuncToInterpolate(params.NumRadial, &Amplitude[0], &Amplitude[i]);
			for (int j = 0; j < num_points-1; j++)
			{
				if(real((*Path)[j])<real(Amplitude[0][params.NumRadial-1]))
				{
				(*AmplitudePath)[j](i-1,0)=FuncToInterpolate.getValue((*Path)[j]);
				//std::cout << i << "  " << j << "  " << (*Path)[j] << "  " << (*AmplitudePath)[j](i,0) << std::endl;
				}
			}
			//cin.get();
		}
	}
std::cout << "On Path calculation finished." << std::endl;
}

void SetBSAonPath_Cauchy(std::vector<t_cmplxMatrix> (*AmplitudePath),t_cmplxArray1D (*Path) ,t_cmplx Q, double M_contour)
{
	std::cout << std::endl;
	std::cout << "On Path calculation..." << std::endl;
	//flag_amp_desciption=false;
	DressBSA_complex(Q, 8, M_contour);
	int num_points;
	double Time;
	num_points=(*Path).size();
	(*AmplitudePath).clear();
	(*AmplitudePath).resize(num_points);
	t_cmplxMatrix Temp(num_amplitudes+1,1);
	std::cout << "Initialized" << std::endl;
	std::cout << "points on the path - " << num_points << std::endl;
		Time=Get_Time();
#pragma omp parallel 
{//start of pragma		
		C_BSE_Hadron_Base * bsa_copy_omp;
		bsa_copy_omp=MakeCopy();
		t_cmplxMatrix Temp_matrix(num_amplitudes,1);
		t_cmplxArray1D Temp_array;
		#pragma omp for 
		for (int i = 0; i < num_points; i++)
		{
			Temp_array=bsa_copy_omp->getCauchyAt_embedded((*Path)[i]);
			(*AmplitudePath)[i].Resize(num_amplitudes,1);
			for (int j = 0; j < num_amplitudes; j++){
				(*AmplitudePath)[i].Element(j,0)=Temp_array[j];
			}
		}
		delete bsa_copy_omp;
}//end of pragma
	std::cout << " Time for continuation spent - " << (Get_Time() - Time)/8.0 << std::endl;
	
std::cout << "On Path calculation finished." << std::endl;
}	
	void setContourAndGrid(double M_contour);
	t_cmplxArray1D getCauchyAt_embedded(t_cmplx coordin);
	void CalcVectorGrid();
	void CalcVectorCont(t_cmplx _P);
	void DressBSA_complex(t_cmplx _P, int steps,double M_contour);
	double CheckBSA_complex(t_cmplx _OldAmp);
	void PreNormBSA_complex();
	void DrawBSA_complex();
};

#include "BSE_Base_onCauchy.hpp"
