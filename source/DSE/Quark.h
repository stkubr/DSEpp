#pragma once
#include "DSE.h"
class C_Quark: public C_Propagator{
	
	protected:
	ifstream * ParamList;
	const char * SavePropPath;
	char * LoadPropPath;
	C_Gluon * Gluon;
	t_cmplxMatrix (C_Quark::*integrand) (t_cmplxArray1D);
	C_Integrator_Line<t_cmplxMatrix,C_Quark,double> * Integ_radial_leg;
	C_Integrator_Line<t_cmplxMatrix,C_Quark,double> * Integ_angle_cheb;
	C_Integrator_Cauchy<t_cmplxArray1D,t_cmplxArray3D,t_cmplx> * Integ_cauchy_long;
	//C_Integrator_Cauchy<dcxArray1D,dcxArray2D,dcx> * Integ_cauchy_short;
	int index_p,grid1_num;
	t_cmplx x;
	t_dArray1D zz_rad, w_rad, zz_angle, w_angle, z_circus, w_circus;
	t_cmplxArray1D integrand_args;
	//C_Quark_parameters params;
	double B_renorm,B_mu,A_renorm,Z2,check_res,eps;
	bool flag_normalized;
	
	t_cmplxVector k,p;
	
	C_Quark(){
		Memory_abs=DedicMemFactory_Quark->CreateMemory();
		Memory=(C_DedicMem_Quark*)Memory_abs;
		
		flag_dressed=false;
		flag_normalized=false;
	}
	
	void ReadParameters(ifstream * _ParamList){
			string line;
			if ((*_ParamList).is_open())
			{
				while ( (*_ParamList).good() )
				{
					(*_ParamList) >> line >> num_amplitudes;
					(*_ParamList) >> line >> params.num_prop_steps;
					(*_ParamList) >> line >> params.num_angle_steps;
					(*_ParamList) >> line >> params.m0;
					(*_ParamList) >> line >> params.mu;
					(*_ParamList) >> line >> params.M2_contour;
					(*_ParamList) >> line >> params.LimDk;
					(*_ParamList) >> line >> params.LimUk;
					(*_ParamList) >> line >> params.EffectiveCutoff;
					(*_ParamList) >> line >> params.Accuracy;
					(*_ParamList) >> line >> params.ReCalcProp;
					(*_ParamList) >> line >> params.HeavyLight;
				}
				cout << "num_amplitudes - " << num_amplitudes << endl;
				cout << params << endl;
			}
			else cout << "Cant open file!" << endl;
		}
		
	t_cmplx getDressingFactor(){
		return Z2;
	}
	
	void setContourApex(double M2){
		params.M2_contour=M2;
	}

// Set parameters to initial values
//----------------------------------------------------------------------
	void InitialState(){
		A_renorm=1.0;
		B_renorm=0.0;
		B_mu=0.0;
		Z2=1.0;
		eps=1.0;
		check_res=0.0;
		flag_dressed=false;
		flag_normalized=false;
		
		ResizeMemory();
	}

// Resize all storages (internal and external), also create side objects like (Integrators, Kernels and etc.)
//----------------------------------------------------------------------
	void InitializateIntegrators(){
		Integ_radial_leg=C_Integrator_Line<t_cmplxMatrix,C_Quark,double>::createIntegrator(params.num_prop_steps, params.LimDk, params.LimUk, num_amplitudes, qgausleg_lin_ID);
		Integ_angle_cheb=C_Integrator_Line<t_cmplxMatrix,C_Quark,double>::createIntegrator(params.num_angle_steps, params.LimDk, params.LimUk, num_amplitudes, qgauscheb_ID);
		Integ_cauchy_long=C_Integrator_Cauchy<t_cmplxArray1D,t_cmplxArray3D,t_cmplx>::createIntegrator(params.num_prop_steps, params.LimDk, params.LimUk, num_amplitudes, qcauchyleg_lin_ID);
		Integ_radial_leg->getNodes(&zz_rad,&w_rad);
		Integ_angle_cheb->getNodes(&zz_angle,&w_angle);
		integrand_args.resize(2);
	}


	void ResizeMemory(){
		Memory->resizeContour(4, num_amplitudes+2, params.num_prop_steps);
		Memory->resizeGrid(1,num_amplitudes+1, 4*params.num_prop_steps, params.num_prop_steps*params.num_angle_steps);
	}
	

// Copier of "this" (used im parallel sections)
//----------------------------------------------------------------------	
	virtual C_Quark * MakeCopy(){
		return new C_Quark(*this);
	}
	
	
// Set Cauchy contour
//----------------------------------------------------------------------
	void setContourAndGrid()
	{		
		double M2=params.M2_contour;
		int n,m;
		t_cmplx Ku,Kd,dk,temp,temp2,t,weight;
		
		for (int i = 0; i <= params.num_prop_steps-1; i++)		
		{
			//cout << grid_k_cx[i] << "   ";
			t=(zz_rad[i+1]);
			weight=-1.0*w_rad[i+1];
			temp=t*t - M2/4.0 + ii*t*sqrt(M2);
			temp2=(2.0*t + ii*sqrt(M2))*weight;
			//cout << t << "  " << temp << endl;
			//temp=t*t + ii*t*t;
			//k_cont_1.push_back(temp);
			Memory->S_cont[0][0][i]=(temp);
			Memory->S_cont[0][1][i]=(temp2);
			Memory->S_cont[0][2][i]=(InitStepA(temp));
			Memory->S_cont[0][3][i]=(InitStepB(temp));
		}
		Kd=t*t - M2/4.0;
		for (int i = 0; i <= params.num_prop_steps-1; i++)		
		{
			//cout << grid_k_cx[i] << "   ";
			t=(zz_rad[i+1]);
			weight=1.0*w_rad[i+1];
			temp=Kd + ii*t*sqrt(M2) + params.LimUk*params.LimDk;
			temp2=ii*sqrt(M2)*weight;
			//k_cont_1.push_back(temp);
			Memory->S_cont[1][0][i]=(temp);
			Memory->S_cont[1][1][i]=(temp2);
			Memory->S_cont[1][2][i]=(InitStepA(temp));
			Memory->S_cont[1][3][i]=(InitStepB(temp));
			//cout << k_cont_2[i] << endl;
		}
		for (int i = 0; i <= params.num_prop_steps-1; i++)		
		{
			//cout << grid_k_cx[i] << "   ";
			t=(zz_rad[i+1]);
			weight=1.0*w_rad[i+1];
			temp=t*t - M2/4.0 - ii*t*sqrt(M2);
			temp2=(2.0*t - ii*sqrt(M2))*weight;
			//temp=t*t - ii*t*t;
			//k_cont_1.push_back(temp);
			Memory->S_cont[2][0][i]=(temp);
			Memory->S_cont[2][1][i]=(temp2);
			Memory->S_cont[2][2][i]=(InitStepA(temp));
			Memory->S_cont[2][3][i]=(InitStepB(temp));
			//cout << k_cont_2[i] << endl;
		}
		for (int i = 0; i <= params.num_prop_steps-1; i++)		
		{
			//cout << grid_k_cx[i] << "   ";
			t=(zz_rad[i+1]);
			weight=-1.0*w_rad[i+1];
			temp=Kd - ii*t*sqrt(M2) + params.LimUk*params.LimDk;
			temp2=-ii*sqrt(M2)*weight;
			//k_cont_1.push_back(temp);
			Memory->S_cont[3][0][i]=(temp);
			Memory->S_cont[3][1][i]=(temp2);
			Memory->S_cont[3][2][i]=(InitStepA(temp));
			Memory->S_cont[3][3][i]=(InitStepB(temp));
			//cout << k_cont_2[i] << endl;
		}
		ofstream countur_data;
		countur_data.open ("../Data_files/countur_data.dat");
		for (int i = 0; i < Memory->S_cont.size(); i++)
		{
			for (int j = 0; j < Memory->S_cont[0][0].size(); j++){
				countur_data << real(Memory->S_cont[i][0][j]) <<"  "<< imag(Memory->S_cont[i][0][j]) << endl;
				//cout << real(Memory->S_cont[i][0][j]) << "  " << imag(Memory->S_cont[i][0][j]) << endl; 
			}
		}
		countur_data.close();
		cout << "Contour has been set! points - " << Memory->S_cont[0][0].size()*Memory->S_cont[0].size() << endl;
		//cout << Memory->S_cont[0][0].size() << "  " << Memory->S_cont[0].size() << "  " << Memory->S_cont.size() << endl;
		//cin.get();
// set (p - k)^2 grid
		for (int num_part = 0; num_part < Memory->S_cont.size(); num_part++)
		{
			for (int i = 0; i < Memory->S_cont[num_part][0].size(); i++)
			{
				//int num_rad_steps=Memory->S_cont[num_part][0].size();
				for (int j = 0; j < params.num_prop_steps; j++)
				{
					for (int k = 0; k < params.num_angle_steps; k++)
					{
						 Memory->S_grid[0][0][i+num_part*params.num_prop_steps][params.num_angle_steps*j + k]=(Memory->S_cont[num_part][0][i] + zz_rad[j+1]*zz_rad[j+1] - 2.0*sqrt(Memory->S_cont[num_part][0][i]*zz_rad[j+1]*zz_rad[j+1])*zz_angle[k+1] );
						 //cout << Memory->S_cont[num_part][0][i] << "  " << zz_rad[j+1]*zz_rad[j+1] << "  " << 2.0*sqrt(Memory->S_cont[num_part][0][i]*zz_rad[j+1]*zz_rad[j+1])*zz_angle[k+1] << endl;
						// cin.get();
					}
					//cout << num_part << "  " << i+num_part*params.num_prop_steps << "  " << j << endl;
				}
				//cout << num_part << "  " << i+num_part*params.num_prop_steps << "  " << Memory->S_cont[num_part][0][i] << "  " << Memory->S_grid[0][0][i+num_part*params.num_prop_steps][0] << endl;
			}
		}
// set (p - k/2)^2 grid 
	/*	for (int i = 0; i < Memory->S_cont[0][0].size(); i++)
		{
			for (int j = 0; j < params.num_prop_steps; j++)
			{
				for (int k = 0; k < params.num_angle_steps; k++)
				{
					 //S_grid2_storage[0][i][params.num_angle_steps*j + k]=(S_cont_storage[0][i] + zz_rad[j]*zz_rad[j]/4.0 - sqrt(S_cont_storage[0][i]*zz_rad[j]*zz_rad[j])*zz_angle[k] );
				}
			}
		}*/
	
//check set p-k grid
		ofstream grid_data;
		grid_data.open ("../Data_files/grid_data.dat");
			int i=0;
			for (int j = 0; j < Memory->S_grid[0][0][0].size(); j++)
			{
				grid_data << real(Memory->S_grid[0][0][i][j]) <<"  "<< imag(Memory->S_grid[0][0][i][j]) << endl;
				//cout << real(Memory->S_grid[0][0][i][j]) << "  " << imag(Memory->S_grid[0][0][i][j]) << endl;
			}
		grid_data.close();
		cout << "Grid has been set! points - " << Memory->S_grid[0][0][0].size()*Memory->S_grid[0][0].size() << endl << endl;
	}
	
// Initial guesses for A and B
//----------------------------------------------------------------------	
	t_cmplx InitStepA (t_cmplx z)
	{
		return  1.0 + exp(-z*z);
	}
	t_cmplx InitStepB (t_cmplx z)
	{
		return 1.0*exp(-z*z);
	}
	
// Multidimentional integration on complex plane
//----------------------------------------------------------------------
	t_cmplxMatrix Multi_INT_cx(t_cmplxMatrix (C_Quark::*func_to_int) (t_cmplxArray1D) )
	{
		integrand=func_to_int;
		return Integ_radial_leg->getResult(&C_Quark::f1,this);
	}
	t_cmplxMatrix f1 (double y)
	{
		integrand_args[0]=y;
		return Integ_angle_cheb->getResult(&C_Quark::f2,this);
	}
	t_cmplxMatrix f2 (double z)
	{
		integrand_args[1]=z;
		return (this->*integrand)(integrand_args);
	}
	
// Evaluate Cauchy integral on contour, at certain point 
//----------------------------------------------------------------------
	t_cmplxArray1D getCauchyAt(t_cmplx point){
		t_cmplx coordin;
		t_cmplxArray1D result;
		t_cmplxArray2D Temp_return;
		double sign;
		t_cmplx tempN=t_cmplx(0.0,0.0);
		result.resize(num_amplitudes);
		Temp_return.resize(4,t_cmplxArray1D(num_amplitudes+1,0.0));
		coordin=point;
		for (int step = 0; step < 4; step++){
			Temp_return[step]=Integ_cauchy_long->getResult(&(Memory->S_cont),step,&coordin);
		}
		
		for (int i = 0; i < 4; i++){
			if(i==0 || i == 3) sign=-1.0; else sign=1.0;
			tempN+=sign*Temp_return[i][0];
			//cout << tempN << endl;
		}
		
		for (int i = 0; i < num_amplitudes; i++){
			for (int j = 0; j < 4; j++){
				if(j==0 || j == 3) sign=-1.0; else sign=1.0;
				result[i]+=sign*Temp_return[j][i+1]/tempN;
				//cout << i << "  " << result[i] << endl;
			}
		}
		return result;
	}
	
// Cauchy integration routine
//----------------------------------------------------------------------
	t_cmplxArray1D getCauchyAt_embedded(t_cmplx coordin)
	{
		int j;//,sumj;sumj=0;
		double zr,zm,dz,t;
		t_cmplxArray1D result(2);
		t_cmplx F1,F2,N,temp,sumF1,sumF2,sumN;
		double aa=params.LimDk, bb=params.LimUk;
		t_cmplx z_i,dz_i;
		t_cmplx resF1,resF2,resN;
		
		sumF1=t_cmplx(0.0,0.0);
		sumF2=t_cmplx(0.0,0.0);
		sumN=t_cmplx(0.0,0.0);
		resF1=t_cmplx(0.0,0.0);
		resF2=t_cmplx(0.0,0.0);
		resN=t_cmplx(0.0,0.0);
        zm=aa;
		zr=(bb-aa);

		for (j=1;j<=params.num_prop_steps;j++) 
		{
			temp = 1.0;
			
			z_i=(Memory->S_cont)[0][0][j-1];
			dz_i=temp*(Memory->S_cont)[0][1][j-1];
			F1=conj(z_i-coordin)/norm(z_i-coordin)*dz_i*(Memory->S_cont)[0][2][j-1];
			F2=conj(z_i-coordin)/norm(z_i-coordin)*dz_i*(Memory->S_cont)[0][3][j-1];
			N=conj(z_i-coordin)/norm(z_i-coordin)*dz_i;
			
			z_i=(Memory->S_cont)[1][0][j-1];
			dz_i=temp*(Memory->S_cont)[1][1][j-1];
			F1+=conj(z_i-coordin)/norm(z_i-coordin)*dz_i*(Memory->S_cont)[1][2][j-1];
			F2+=conj(z_i-coordin)/norm(z_i-coordin)*dz_i*(Memory->S_cont)[1][3][j-1];
			N+=conj(z_i-coordin)/norm(z_i-coordin)*dz_i;
			
			z_i=(Memory->S_cont)[2][0][j-1];
			dz_i=temp*(Memory->S_cont)[2][1][j-1];
			F1+=conj(z_i-coordin)/norm(z_i-coordin)*dz_i*(Memory->S_cont)[2][2][j-1];
			F2+=conj(z_i-coordin)/norm(z_i-coordin)*dz_i*(Memory->S_cont)[2][3][j-1];
			N+=conj(z_i-coordin)/norm(z_i-coordin)*dz_i;
			
			z_i=(Memory->S_cont)[3][0][j-1];
			dz_i=temp*(Memory->S_cont)[3][1][j-1];
			F1+=conj(z_i-coordin)/norm(z_i-coordin)*dz_i*(Memory->S_cont)[3][2][j-1];
			F2+=conj(z_i-coordin)/norm(z_i-coordin)*dz_i*(Memory->S_cont)[3][3][j-1];
			N+=conj(z_i-coordin)/norm(z_i-coordin)*dz_i;
		
			sumF1 += F1;
			sumF2 += F2;
			sumN += N;
			
		}
		
		resN = sumN*zr;
		resF1 = sumF1*zr;
		resF2 = sumF2*zr;
		//cout << sumF1*zr << "  " << sumF2*zr << "  " << sumN *zr<< endl;
		result[0]=(resF1/resN);
		result[1]=(resF2/resN);
		//cout << resF1/resN << "  " << resF2/resN << endl;
		return result;
	}
	
	void ExportForKernel(t_cmplxArray2D * KernelQuarkStorage)
	{
		int j;//,sumj;sumj=0;
		double zr,zm,dz,t;
		t_cmplx F1,F2,N,temp,sumF1,sumF2,sumN;
		double aa=params.LimDk, bb=params.LimUk;
		t_cmplx z_i,dz_i;
		t_cmplx resF1,resF2,resN;
		
		sumF1=t_cmplx(0.0,0.0);
		sumF2=t_cmplx(0.0,0.0);
		sumN=t_cmplx(0.0,0.0);
		resF1=t_cmplx(0.0,0.0);
		resF2=t_cmplx(0.0,0.0);
		resN=t_cmplx(0.0,0.0);
        zm=aa;
		zr=(bb-aa);
		
		(*KernelQuarkStorage).resize(3,t_cmplxArray1D(4*params.num_prop_steps));
		
		for (j=1;j<=params.num_prop_steps;j++) 
		{
			(*KernelQuarkStorage)[0][j-1]=(Memory->S_cont)[0][0][j-1];
			(*KernelQuarkStorage)[1][j-1]=(Memory->S_cont)[0][1][j-1];
			(*KernelQuarkStorage)[2][j-1]=(Memory->S_cont)[0][3][j-1];
			
			(*KernelQuarkStorage)[0][j-1 + params.num_prop_steps]=(Memory->S_cont)[1][0][j-1];
			(*KernelQuarkStorage)[1][j-1 + params.num_prop_steps]=(Memory->S_cont)[1][1][j-1];
			(*KernelQuarkStorage)[2][j-1 + params.num_prop_steps]=(Memory->S_cont)[1][3][j-1];
			
			(*KernelQuarkStorage)[0][j-1 + 2*params.num_prop_steps]=(Memory->S_cont)[2][0][j-1];
			(*KernelQuarkStorage)[1][j-1 + 2*params.num_prop_steps]=(Memory->S_cont)[2][1][j-1];
			(*KernelQuarkStorage)[2][j-1 + 2*params.num_prop_steps]=(Memory->S_cont)[2][3][j-1];
			//cout << -1.0*(Memory->S_cont)[2][3][j-1] << "  " << temp*(Memory->S_cont)[2][1][j-1] << endl;
			
			(*KernelQuarkStorage)[0][j-1 + 3*params.num_prop_steps]=(Memory->S_cont)[3][0][j-1];
			(*KernelQuarkStorage)[1][j-1 + 3*params.num_prop_steps]=(Memory->S_cont)[3][1][j-1];
			(*KernelQuarkStorage)[2][j-1 + 3*params.num_prop_steps]=(Memory->S_cont)[3][3][j-1];
		}
	}
		
// Evaluate Cauchy integral on contour, obtain Propogator on grid
//----------------------------------------------------------------------
	void CalcPropGrid()
	{
		int num_grid;
		int num_threads=4;
		//double L=LimUk*LimUk*EffectiveCutoff;
		num_grid=Memory->S_grid[0][0][0].size();
		cout << "Grid" << endl;
		
		//Init_time=Get_Time();
		
#pragma omp parallel num_threads(4)
{//begin of pragma
		t_cmplx coordin;
		C_Quark * quark_copy;
		quark_copy=MakeCopy();
		//quark_copy->ResetIntegrators();
		//C_Integrator_Cauchy<dcxArray1D,dcxArray3D,dcx> * copy_Integ_Cauchy
		t_cmplxArray1D S_temp_storage(num_amplitudes);
		#pragma omp for// ordered
		for (int step = 0; step < 4 ; step++)
		{
		for (int i = 0; i < params.num_prop_steps; i++)
		{
			for (int j = 0; j < num_grid; j++)
			{
				coordin=Memory->S_grid[0][0][i+step*params.num_prop_steps][j];
				if (real(coordin)<params.LimUk*params.LimUk*params.EffectiveCutoff)//LimUk*LimUk*EffectiveCutoff
				{
					S_temp_storage=quark_copy->getCauchyAt_embedded(coordin);
					Memory->S_grid[0][1][i+step*params.num_prop_steps][j]=(S_temp_storage[0]);
					Memory->S_grid[0][2][i+step*params.num_prop_steps][j]=(S_temp_storage[1]);
					//cout <<  coordin << "  " << getCauchyAt(coordin)[0] << endl;
					//cin.get();
				}
				else
				{
					//cout << coordin << "  " << i+step*params.num_prop_steps << endl;
					Memory->S_grid[0][1][i+step*params.num_prop_steps][j]=Z2*(1.0*params.HeavyLight -(params.HeavyLight-1.0)*real(coordin)) ;
					Memory->S_grid[0][2][i+step*params.num_prop_steps][j]=Z2*(params.m0);
				}
			}
		}
		}
		delete quark_copy;
}//end of pragma

	if (flag_normalized==true)
	{
		A_renorm=real(getCauchyAt_embedded(params.mu*params.mu)[0])-Z2*1.0;
		Z2=1.0 - (A_renorm);	
		Kernel->setZ2DressingFactor(this->getDressingFactor());
	}
	//Kernel->SetPionSign(-1.0);
	//MemoryManager->CopyMemoryFrom(this->Memory,Kernel->Memory);
	
	//cout << "Z2 - " << "  " << Z2 << "  " << "m_renorm - " B_mu << endl;
		
	}
		
// Evaluate DSE integral on grid, obtain Propagator on contour
//----------------------------------------------------------------------
	void CalcPropCont()
	{
		t_cmplx temp;
		cout << "Contour" << endl;
		int n;
		n=Memory->S_cont[0][0].size();
#pragma omp parallel 
{
		C_Quark * quark_copy;
		quark_copy=MakeCopy();
		//quark_copy->ResetIntegrators();
		t_cmplxMatrix Temp_return(num_amplitudes,1),Temp_bare(num_amplitudes,1),Temp_return2(num_amplitudes,1);
		#pragma omp for
		for (int i = 0; i < n; i++)
		{
			quark_copy->index_p=i;
			quark_copy->x=Memory->S_cont[0][0][i];
			
			Temp_bare(0,0)=Z2;
			Temp_bare(1,0)=params.m0 - B_renorm + params.HeavyLight;
			
			quark_copy->grid1_num=0;
			
			if(real(quark_copy->x)<params.LimUk*params.LimUk*1.1){
			Temp_return=Temp_bare + quark_copy->Multi_INT_cx(&C_Quark::Integrand_numerical);quark_copy->grid1_num=0;
			Memory->S_cont[0][2][i]=Temp_return(0,0);
			Memory->S_cont[0][3][i]=Temp_return(1,0); 
			Memory->S_cont[1][2][i]=Z2;
			Memory->S_cont[1][3][i]=Z2*params.m0;
			Memory->S_cont[2][2][i]=conj(Temp_return(0,0));
			Memory->S_cont[2][3][i]=conj(Temp_return(1,0));
			Memory->S_cont[3][2][i]=Z2;
			Memory->S_cont[3][3][i]=Z2*params.m0;
			}
			else{
				Memory->S_cont[0][2][i]=Z2;
				Memory->S_cont[0][3][i]=Z2*params.m0; 
				Memory->S_cont[1][2][i]=Z2;
				Memory->S_cont[1][3][i]=Z2*params.m0;
				Memory->S_cont[2][2][i]=Z2;
				Memory->S_cont[2][3][i]=Z2*params.m0;
				Memory->S_cont[3][2][i]=Z2;
				Memory->S_cont[3][3][i]=Z2*params.m0;
			}
			//cout << Memory->S_cont[0][0][i] << "  " << Memory->S_cont[0][2][index_p] << "  " << Memory->S_cont[0][3][index_p] << endl;
			//cin.get();
		}
		delete quark_copy;
}
	}

// Analitic form of the Integrand (avaible only for RL or Pion Contribution)
//----------------------------------------------------------------------
	t_cmplxMatrix Integrand_analitic (t_cmplxArray1D values)
	{
		grid1_num++;
		t_cmplx y,z;
		y=values[0];
		z=values[1];
		t_cmplxMatrix result(num_amplitudes,1);
		t_cmplx _A,_B,_B2,y2,pk,epsilon;
		y2=y*y;
		pk=Memory->S_grid[0][0][index_p][grid1_num-1];
		_A=Memory->S_grid[0][1][index_p][grid1_num-1];
		_B=Memory->S_grid[0][2][index_p][grid1_num-1];
		//_B2=Memory->S_grid[0][2][index_p][grid1_num-1];
		epsilon=1.0/(pk*_A*_A+_B*_B);
		//if(abs(epsilon)>100.0) epsilon=0.0;
		result(0,0)=Z2*Z2*2.0/(8.0*pi*pi*pi)*(_A*y*y*y*epsilon)*4.0/3.0*Gluon->GetGluonAt(&y2)*(1.0 + 2.0*z*z - 3.0*sqrt(y*y/(x))*z);
		result(1,0)=Z2*Z2*2.0*3.0/(8.0*pi*pi*pi)*(_B*y*y*y*epsilon*4.0/3.0*Gluon->GetGluonAt(&y2));
		//cout << y2 << "  " << pk << "  " << (_A*y*y*y*epsilon)*Gluon->GetGluonAt(&y2)*(1.0 + 2.0*z*z - 3.0*sqrt(y*y/(x))*z) << endl;
		//cin.get();
		//- PionContrib*Flavor_factor*Z2*2.0/(8.0*pi*pi*pi)*(1.0 - sqrt(y*y/(x))*z)*(_A*y*y*y*epsilon)*_B2/f_pion*(1.0/(y*y  + M_pion*M_pion));
		return result;
	}
	
// Set k and p vectors for the Numerical Integrand
//----------------------------------------------------------------------	
	void setKinematic(t_cmplx x, t_cmplx y, t_cmplx z){
		k.SetP4(0.0,0.0,sqrt(1.0-z*z)*y,y*z);
		p.SetP4(0.0,0.0,0.0,sqrt(x));
	}

// Numerical form of the Integrand (avaible for general Kernel)
//----------------------------------------------------------------------
	t_cmplxMatrix Integrand_numerical (t_cmplxArray1D values)
	{
		grid1_num++;
		t_cmplx y,z,kinematic_factor,temp;
		y=values[0];
		z=values[1];
		setKinematic(x,y,z);
		t_cmplxMatrix result(num_amplitudes,1);
		t_cmplxDirac S(4,4),Proj1(4,4),Proj2(4,4);
		t_cmplxVector pion_momenta;
		pion_momenta=(p-k/2.0);
		t_cmplx _A,_B,_B2,y2,pk,epsilon;
		y2=y*y;
		pk=Memory->S_grid[0][0][index_p][grid1_num-1];
		
		_A=Memory->S_grid[0][1][index_p][grid1_num-1];
		_B=Memory->S_grid[0][2][index_p][grid1_num-1];
		//_B2=Memory->S_grid[0][2][index_p][grid1_num-1];
		epsilon=1.0/(pk*_A*_A+_B*_B);
		
		kinematic_factor=1.0/(16.0*pi*pi*pi);
		S=-ii*((p-k)*Y)*_A + I*_B;
		
		Proj1=-ii*(p*Y)/(p*p);
		Proj2=I;

		result(0,0)=kinematic_factor*(y*y*y*epsilon)*Kernel->TraceKernelWithoutStoring(&Proj1,&S,&k,&p,&pion_momenta,true);
		result(1,0)=kinematic_factor*(y*y*y*epsilon)*Kernel->TraceKernelWithoutStoring(&Proj2,&S,&k,&p,&pion_momenta,false);

		//temp=Z2*Z2*2.0*4.0/(3.0*8.0*pi*pi*pi)*(_A*y*y*y*epsilon)*Gluon->GetGluonAt(&y2)*(1.0 + 2.0*z*z - 3.0*sqrt(y*y/(x))*z)
		//temp=-3.0*Z2*2.0/(8.0*pi*pi*pi)*(1.0 - sqrt(y*y/(x))*z)*(_A*y*y*y*epsilon)*getCauchyAt_embedded(pk)[1]/0.093*(1.0/(y*y  + 0.138*0.138));
		
	/*	#pragma omp master
		{
		cout << result(0,0) << "  " << temp << "  " << getCauchyAt_embedded(1.0)[1] << "  " << pk << endl;
		cin.get();
		}
	*/	//- PionContrib*Flavor_factor*Z2*2.0/(8.0*pi*pi*pi)*(1.0 - sqrt(y*y/(x))*z)*(_A*y*y*y*epsilon)*_B2/f_pion*(1.0/(y*y  + M_pion*M_pion));
		return result;
	}
	
// Calculation A,B,M,Sigma_V,Sigma_S at point in contour
//----------------------------------------------------------------------
	t_cmplxArray1D getPropAt(t_cmplx q) 
	{
		t_cmplx temp,t;
		t_cmplxArray1D tempvector;
		t_cmplxArray1D storage(5);
		
		DressPropagator();
		
		//if (real(q)<params.LimUk*params.LimUk){
			tempvector=getCauchyAt_embedded(q);
		/*}
		else{
			dcxArray1D asymptotics(2);
			asymptotics[0]=Z2;
			asymptotics[0]=Z2*params.m0;
			tempvector=asymptotics;
		}*/
		
				
		// A
		storage[0]=tempvector[0];
		// B
		storage[1]=tempvector[1];
		// M
		storage[2]=storage[1]/storage[0];
		// Sigma_V
		storage[3]=(storage[0]/(q*storage[0]*storage[0]+storage[1]*storage[1]));
		// Sigma_S
		storage[4]=(storage[1]/(q*storage[0]*storage[0]+storage[1]*storage[1]));	
		
		return storage;
	}
	
// Calculate consequently Grid and Contour until converge
//----------------------------------------------------------------------	
	void PropSetAndCheck()
	{
		if(params.ReCalcProp)
		{
			for (int i = 0; i <2 ; i++)
			{
				CalcPropGrid();
				CalcPropCont();
				write_Prop_re(100);
			}
			//StopLine();
			flag_normalized=true;
			check_res=1.0;
			eps=1.0;
			while (eps>params.Accuracy)
			{
				CalcPropGrid();
				CalcPropCont();
				B_mu=real(getCauchyAt_embedded(params.mu*params.mu)[1]);
				if (flag_normalized==true) B_renorm+=real(getCauchyAt_embedded(params.mu*params.mu)[1])-params.m0;
				//Z4=2.0-B_renorm/m0;
				PropCheck(100);
				write_Prop_re(100);
			}
		}	
	}
	
// Check convergence
//----------------------------------------------------------------------
	void PropCheck(int s)
	{
		double Pu,Pd,x,scale,dp;
		t_cmplxArray1D storage;
		double res=0;
		scale = s;
		Pd=0.01;
		Pu=params.LimUk*0.8;
		dp=pow(10,(log10(Pu/Pd)/scale));
		vector<double> A_check(scale),B_check(scale);
			
		x=Pd;
		for (int i = 0; i < scale; i++)
		{
		storage=getCauchyAt_embedded(x*x);
		A_check[i]=real(storage[0]);
		B_check[i]=real(storage[1]);
		res+=A_check[i]*B_check[i];
		x*=dp;
		}
		
		eps=fabs(res - check_res)/fabs(res);
		cout << "Z2 - " <<"  "<< Z2 <<"  "<< "A_mu" <<"  "<< real(getCauchyAt_embedded(params.mu*params.mu)[0]) <<"  "<<
				"m_renorm - " <<"  "<< B_mu /*_ "Z4 - " << "  " << Z4*/ <<"  "<< "Accuracy - " <<"  "<< eps << endl;
		check_res=res;
	}
	
// Initialization (Dressing) of the Propagator
//----------------------------------------------------------------------	
	void DressPropagator(){
		//Kernel->GetNameID();
		if (flag_dressed==false)
		{
			PrintLine('-');
			cout << "Start Dressing for " << name << " with quark mass -" <<"  "<< params.m0 <<"  "<< "and contour -" <<"  "<< params.M2_contour << endl;
			PrintLine('-');
			setContourAndGrid();
			LoadPropCountour();
			MemoryManager->CopyMemoryFrom(this->Memory,Kernel->Memory);
			PropSetAndCheck();
			flag_dressed=true;
			Kernel->setZ2DressingFactor(this->getDressingFactor());
			//MemoryManager->CopyMemoryFrom(this->Memory,Kernel->Memory);
			//SavePropCountour();
			//Kernel->TakeObjectForDressing(this);
			//ExportPropagator();
			Memory->RemoveGrid();
			//CalcCondensate();
			//ExportPropagator();
			write_Prop_re(100);
		}
	}
	
// Drow Propagator at real line
//----------------------------------------------------------------------
	void write_Prop_re( int s)
	{
		double Pu,Pd,x,scale,dp;
		t_cmplxArray1D storage;
		scale = s;
		Pd=0.01;
		Pu=params.LimUk;
		dp=pow(10,(log10(Pu/Pd)/scale));
		//dp=(Pu-Pd)/scale;
			
		x=Pd;
		ofstream data_Prop_re;
		data_Prop_re.open ("../Data_files/data_Prop_re_new.dat");
		
		for (int i = 0; i <= scale; i++)
		{
			storage=getCauchyAt_embedded(x*x);
			data_Prop_re << x*x <<"  "<< real(storage[0]) <<"  "<< real(storage[1])  << endl;
			x*=dp;
			//cout << "computing real part...  " << 100.0/(scale+1)*(i+1) << "%\n";
		}
		data_Prop_re.close();
	}
	
// Write Propagator to file
//----------------------------------------------------------------------
	void SavePropCountour()
	{
		ofstream SavePropStream;
		SavePropStream.open(SavePropPath);
		(SavePropStream) << "Z2" << '\t' << Z2 << endl;
		for (int i = 0; i < 4; i++){
			for (int j=1;j<=params.num_prop_steps;j++) 
			{
				(SavePropStream) << Memory->S_cont[i][0][j-1] << '\t' << Memory->S_cont[i][2][j-1] << '\t' << Memory->S_cont[i][3][j-1]  << endl;
			}
		}
		SavePropStream.close();
		cout << "Propagator was saved." << endl;	
	}
	
// Read Propagator from file
//----------------------------------------------------------------------
	void LoadPropCountour()
	{
		string dummy;
		if (!params.ReCalcProp)
		{
			ifstream PropContourStream;
			PropContourStream.open(SavePropPath); 
			PropContourStream >> dummy >> Z2;
			if (PropContourStream.is_open())
			{
				for (int i = 0; i < 4; i++){
					for (int j=1;j<=params.num_prop_steps;j++) 
					{
						
						PropContourStream >> Memory->S_cont[i][0][j-1] >> Memory->S_cont[i][2][j-1] >> Memory->S_cont[i][3][j-1];
						//cout << i << "  " << j << "  " << Memory->S_cont[i][0][j-1] << "  " << Memory->S_cont[i][2][j-1] << "  " << Memory->S_cont[i][3][j-1] << endl;
					}
				}
				cout << "Propagator was loaded. And WILL NOT be recalculated" << endl;
			}
			else cout << "Cant open file!" << endl;
			PropContourStream.close();
		} else cout << "Propagator was NOT loaded. But WILL be recalculated." << endl;
	}
	
	void ExportPropagator()
	{
		int j;//,sumj;sumj=0;
		t_cmplx z_i,dz_i,temp;
		ofstream PropContourStream;
		PropContourStream.open ("../Data_files/PropContourExportData_new.dat");
		cout << "Writing to file - \"Data_files/PropContourExportData_new.dat\"  " << endl;
		PropContourStream << fixed;
		PropContourStream << setprecision (10);
		PropContourStream << 0.0001 << '\t' << 2000 << '\t' << Z2  << '\t' << endl;
		for (int i = 0; i < 4; i++){
			for (j=1;j<=params.num_prop_steps;j++) 
			{
				t_cmplx temp=1.0;
				z_i=(Memory->S_cont)[i][0][j-1];
				dz_i=temp*(Memory->S_cont)[i][1][j-1];
				PropContourStream << real(z_i) << '\t' << imag(z_i) << '\t' << real(dz_i) << '\t' << imag(dz_i)  << '\t' << real(temp) << '\t' << real((Memory->S_cont)[i][2][j-1]) << '\t' << imag((Memory->S_cont)[i][2][j-1]) << '\t' << real((Memory->S_cont)[i][3][j-1]) << '\t' << imag((Memory->S_cont)[i][3][j-1]) << endl;
			}
		}	
		PropContourStream.close();
		
		cout << "The quark propagator has beed exported" << endl;
	}
	
	void SetQuarkonPath(vector<t_cmplxMatrix> (*AmplitudePath),t_cmplxArray1D (*Path))
	{
		cout << endl;
		cout << "Quark on Path calculation..." << endl;
		int num_points;
		num_points=(*Path).size();
		(*AmplitudePath).resize(num_points);
		t_cmplxArray1D temp_quark(5);
		t_cmplxMatrix Temp(2,1);
		cout << "Initialized" << endl;
		cout << "points on the path - " << num_points << endl;
		for (int i = 0; i < num_points; i++)
		{
			temp_quark=getPropAt((*Path)[i]);
			Temp(0,0)=temp_quark[3];
			Temp(1,0)=temp_quark[4];
			//cout << i << "  " << (*Path)[i] << "  " << Temp(0,0) << "  " << Temp(1,0) << endl;
			(*AmplitudePath)[i]=Temp;
		}
		cout << "Quark on Path calculation finished." << endl;
	}
	
	t_dArray1D GetResult()
		{
			double Pu,Pd,x,scale,dp;
			t_cmplxArray1D storage;
			scale = 100;
			Pd=0.01;
			Pu=params.LimUk;
			dp=pow(10,(log10(Pu/Pd)/scale));

			t_dArray1D temp_result(2,0);

			x=Pd;
			for (int i = 0; i < scale; i++)
			{
				storage=getCauchyAt_embedded(x*x);
				temp_result[0]+=real(storage[0]);
				temp_result[1]+=real(storage[1]);
				x*=dp;
			}
			cout.precision(12);
			cout << temp_result[0] <<" "<< temp_result[1] <<endl;
			return temp_result;
		}

	public:
    static C_Quark* createQuark( Quark_ID id );
    C_DedicMem_Quark * Memory;
		
};

