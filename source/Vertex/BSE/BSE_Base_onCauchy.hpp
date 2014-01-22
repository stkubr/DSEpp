#pragma once

	void C_BSE_Hadron_Base::setContourAndGrid(double M_contour){
		t_cmplx M=t_cmplx(0.0,M_contour);
		int n,m;
		t_cmplx Ku,Kd,temp,dif_temp,t;
		Memory->resizeBSEContour(num_amplitudes+2,params.NumRadial_Contour*4);
		Memory->resizeBSEGrid(num_amplitudes+1, params.NumRadial_Contour*4, params.NumRadial*params.NumCheb_nod1*params.NumCheb_nod1*params.NumAngleY);
		
		
		
		for (int i = 1; i <= params.NumRadial_Contour; i++)		
		{
			t=(zz_cauchy[i]);
			temp=t*t + M*M/4.0 + t*M;
			dif_temp=-1.0*w_cauchy[i]*(2.0*t + M);
			Memory->CauchyContour[0][i-1]=temp;
			Memory->CauchyContour[1][i-1]=dif_temp;
			for (int j = 0; j < num_amplitudes ; j++){Memory->CauchyContour[j+2][i-1]=exp(-real(temp));}
		}
		Kd=t*t + M*M/4.0;
		for (int i = 1; i <= params.NumRadial_Contour; i++)		
		{
			t=(zz_cauchy[i]);
			temp=Kd + t*M + params.LimUk*params.LimDk;
			dif_temp=w_cauchy[i]*M;
			Memory->CauchyContour[0][i-1+params.NumRadial_Contour]=temp;
			Memory->CauchyContour[1][i-1+params.NumRadial_Contour]=dif_temp;
			for (int j = 0; j < num_amplitudes ; j++){Memory->CauchyContour[j+2][i-1+params.NumRadial_Contour]=exp(-real(temp));}
		}
		for (int i = 1; i <= params.NumRadial_Contour; i++)		
		{
			t=(zz_cauchy[i]);
			temp=t*t + M*M/4.0 + conj(t*M);
			dif_temp=w_cauchy[i]*(2.0*t + conj(M));
			Memory->CauchyContour[0][i-1+params.NumRadial_Contour*2]=temp;
			Memory->CauchyContour[1][i-1+params.NumRadial_Contour*2]=dif_temp;
			for (int j = 0; j < num_amplitudes ; j++){Memory->CauchyContour[j+2][i-1+params.NumRadial_Contour*2]=exp(-real(temp));}
		}
		for (int i = 1; i <= params.NumRadial_Contour; i++)		
		{
			t=(zz_cauchy[i]);
			temp=Kd + conj(t*M) + params.LimUk*params.LimDk;
			dif_temp=-1.0*w_cauchy[i]*conj(M);
			Memory->CauchyContour[0][i-1+params.NumRadial_Contour*3]=temp;
			Memory->CauchyContour[1][i-1+params.NumRadial_Contour*3]=dif_temp;
			for (int j = 0; j < num_amplitudes ; j++){Memory->CauchyContour[j+2][i-1+params.NumRadial_Contour*3]=exp(-real(temp));}
		}
		ofstream countur_data;
		countur_data.open ("../Data_files/Vector_Contour.dat");
		n=Memory->CauchyContour[0].size();
		for (int i = 0; i < n; i++)
		{
			countur_data << real(Memory->CauchyContour[0][i]) << "  " << imag(Memory->CauchyContour[0][i]) << std::endl; 
		}
		countur_data.close();
		std::cout << "Contour has been set! points - "<< n << std::endl;
		
		
		// set (p - k)^2 grid 
		t_cmplxVector k_temp, p_temp, P_temp;
		t_cmplx _p2;
		double _k2,_z,_y,_zp;
		for (int i = 0; i < Memory->CauchyContour[0].size(); i++)
		{
			_p2=Memory->CauchyContour[0][i];
			int z_and_Int_ctr=0;
			for (int zp_ctr = 1; zp_ctr <= params.NumCheb_nod1 ; zp_ctr++)
			{
				_zp=zz_cheb[zp_ctr];
				p_temp.SetP4(0.0, 0.0, sqrt(_p2*(1.0-_zp*_zp)), _zp*sqrt(_p2));
				for (int k_ctr = 1; k_ctr <= params.NumRadial; k_ctr++)
				{
					_k2=zz_rad[k_ctr];
					for (int z_ctr = 1; z_ctr <= params.NumCheb_nod1 ; z_ctr++)
					{
						_z=zz_cheb[z_ctr];
						for (int y_ctr = 1; y_ctr <= params.NumAngleY; y_ctr++)
						{
							_y=zz_angleY[y_ctr];
							k_temp.SetP4(0.0, sqrt((1.0-_z*_z)*_k2*(1.0-_y*_y)), _y*sqrt((1.0-_z*_z)*_k2), sqrt(_k2)*_z );
							Memory->CauchyGrid[0][i][z_and_Int_ctr]=((p_temp-k_temp)*(p_temp-k_temp));
							//std::cout << i << "  " << zp_ctr << "  " << k_ctr << "  " << z_ctr << "  " << y_ctr << "  " << z_and_Int_ctr << std::endl;
							z_and_Int_ctr++;
						}
					}
				}
			}
		}
		
		m=Memory->CauchyGrid[0][0].size();
		ofstream grid_data;
		grid_data.open ("../Data_files/Vector_Grid.dat");
			int i=0;
			for (int j = 0; j < m; j++)
			{
				grid_data << real(Memory->CauchyGrid[0][i][j]) << "  " << imag(Memory->CauchyGrid[0][i][j]) << std::endl;
			}
		grid_data.close();
		std::cout << "Grid has been set! points - " << m*n << std::endl << std::endl;
	}
	
	t_cmplxArray1D C_BSE_Hadron_Base::getCauchyAt_embedded(t_cmplx coordin){
		t_cmplxArray1D result(num_amplitudes);
		t_cmplxArray1D F(num_amplitudes);
		t_cmplx z_i,dz_i,N;

		for (int i = 0; i < num_amplitudes; i++){F[i]=0.0;}
		N=t_cmplx(0.0,0.0);
		
		//std::cout << Memory->VertexDressings[0][0][0].size() << std::endl;
		
		for (int j=1;j<= Memory->CauchyContour[0].size();j++) 
		{
			z_i=Memory->CauchyContour[0][j-1];
			dz_i=Memory->CauchyContour[1][j-1];
			for (int i = 0; i < num_amplitudes; i++) F[i]+=conj(z_i-coordin)/norm(z_i-coordin)*dz_i*Memory->CauchyContour[i+2][j-1];
			N+=conj(z_i-coordin)/norm(z_i-coordin)*dz_i;
		}
		for (int i = 0; i < num_amplitudes; i++) result[i]=(F[i]/N);		
		return result;
	}
	
	
	
	void C_BSE_Hadron_Base::CalcVectorGrid(){
		std::cout << "Vector Grid" << std::endl;
#pragma omp parallel //num_threads(1)
{//begin of pragma
		t_cmplx coordin;
		t_cmplxArray1D Vector_temp_storage(num_amplitudes);
		#pragma omp for
		for (int i = 0; i < Memory->CauchyContour[0].size(); i++)
		{
			for (int j = 0; j < Memory->CauchyGrid[0][0].size(); j++)
			{
				coordin=Memory->CauchyGrid[0][i][j];
				if (real(coordin)<params.LimUk)
				{
					Vector_temp_storage=getCauchyAt_embedded(coordin);
					for (int amp = 0; amp < num_amplitudes; amp++)Memory->CauchyGrid[amp+1][i][j]=Vector_temp_storage[amp];
				}
				else
				{
					for (int amp = 0; amp < num_amplitudes; amp++)Memory->CauchyGrid[amp+1][i][j]=0.0;//1.0/real(coordin);
					if(flag_off_shell) Memory->CauchyGrid[1][i][j]=1.0;
				}
			}
		}
}//end of pragma
	}
	
	void C_BSE_Hadron_Base::CalcVectorCont(t_cmplx _P){
		double temp_time;
		temp_time=Get_Time();
		std::cout << "Vector Contour" << std::endl;
#pragma omp parallel //num_threads(1)
{//start of pragma
		C_BSE_Hadron_Base * temp_copy;
		temp_copy=MakeCopy();
		t_cmplxMatrix Temp_matrix;
		#pragma omp for schedule(dynamic)
		for (int i = 0; i < Memory->CauchyContour[0].size()/2; i++)
		{
			temp_copy->index_p=i;
			if (real(Memory->CauchyContour[0][i])<params.LimUk*0.9)
			{
				Temp_matrix=temp_copy->CalcBSA(sqrt(Memory->CauchyContour[0][i]),_P,1); 
				for (int amp = 0; amp < num_amplitudes; amp++)Memory->CauchyContour[amp+2][i]=Temp_matrix(amp,0);
				for (int amp = 0; amp < num_amplitudes; amp++)Memory->CauchyContour[amp+2][i+Memory->CauchyContour[0].size()/2]=conj(Temp_matrix(amp,0));
			} 
			else
			{
				double x_cutoff=real(Memory->CauchyContour[0][i]);
				
				for (int amp = 0; amp < num_amplitudes; amp++)Memory->CauchyContour[amp+2][i]=1.0/x_cutoff/x_cutoff;
				for (int amp = 0; amp < num_amplitudes; amp++)Memory->CauchyContour[amp+2][i+Memory->CauchyContour[0].size()/2]=1.0/x_cutoff/x_cutoff;
				if(flag_off_shell) { 
					Memory->CauchyContour[2][i]=Parton_P->getDressingFactor()/(1+sqrt(x_cutoff)/2000.0); 
					Memory->CauchyContour[2][i+Memory->CauchyContour[0].size()/2]=Parton_P->getDressingFactor()/(1+sqrt(x_cutoff)/2000.0); 
				}
			}
		}
		delete temp_copy;
}//end of pragma

	std::cout << "Time spent - " << (Get_Time()-temp_time)/omp_get_max_threads() << std::endl;
	t_cmplxArray1D TempArray(num_amplitudes);
	ofstream temp_continuation;
	temp_continuation.open ("../Data_files/BSEs_onCauchy.dat");
	for (int i = 1; i < zz_rad.size(); i++)
		{
			TempArray=getCauchyAt_embedded(zz_rad[i]);
			
			temp_continuation << (zz_rad[i]);
			for (int amp = 0; amp < num_amplitudes; amp++) temp_continuation  << '	' <<  real(TempArray[amp]);
			temp_continuation << std::endl;
			
			std::cout << (zz_rad[i]);
			for (int amp = 0; amp < num_amplitudes; amp++) std::cout << "  " << (TempArray[amp]);
			std::cout << std::endl;
		}
	temp_continuation.close();
	}	
	
	void C_BSE_Hadron_Base::DressBSA_complex(t_cmplx _P, int steps, double M_contour){
		SetDressing_ref=&C_BSE_Hadron_Base::SetDressing_shifted;
		GetBSA_ref=&C_BSE_Hadron_Base::GetBSA;
		setContourAndGrid(M_contour);
		double eps=0.01;
		double Accuracy=1.0;
		t_cmplx OldAmp=0.0;
		int dummy_ctr=0;
		while((dummy_ctr < steps) || Accuracy>eps*0.2)
		{
			CalcVectorGrid();
			CalcVectorCont(_P);
			
			Accuracy=CheckBSA_complex(OldAmp);
			OldAmp=getCauchyAt_embedded(zz_rad[1])[0];
			
			PreNormBSA_complex();
			std::cout << "Interation Step - " << dummy_ctr << std::endl;
			dummy_ctr++;
		}
			CalcVectorGrid();
			CalcVectorCont(_P);
		t_cmplxArray1D TempArray(num_amplitudes);
		ofstream temp_continuation;
		temp_continuation.open ("../Data_files/BSEs_onCauchy.dat");
		for (int i = 1; i < zz_rad.size(); i++)
		{
			TempArray=getCauchyAt_embedded(zz_rad[i]);
			
			temp_continuation << (zz_rad[i]);
			for (int amp = 0; amp < num_amplitudes; amp++) temp_continuation  << '	' <<  real(TempArray[amp]);
			temp_continuation << std::endl;
			
			std::cout << (zz_rad[i]);
			for (int amp = 0; amp < num_amplitudes; amp++) std::cout << "  " << (TempArray[amp]);
			std::cout << std::endl;
		}
		temp_continuation.close();
		DrawBSA_complex();
	}
	
	void C_BSE_Hadron_Base::PreNormBSA_complex(){
		t_cmplx prenorm_factor;
		if(!flag_off_shell){
			prenorm_factor=getCauchyAt_embedded(zz_rad[1])[0];
			for (int j=1;j<= Memory->CauchyContour[0].size();j++) 
			{
				for (int i = 0; i < num_amplitudes; i++) Memory->CauchyContour[i+2][j-1]=Memory->CauchyContour[i+2][j-1]/prenorm_factor;
			}
		}
	}
	
	double C_BSE_Hadron_Base::CheckBSA_complex(t_cmplx _OldAmp){
		t_cmplx prenorm_factor;
		prenorm_factor=getCauchyAt_embedded(zz_rad[1])[0];
		return abs(real(_OldAmp) - real(prenorm_factor));
	}
	
	void C_BSE_Hadron_Base::DrawBSA_complex(){
		t_cmplx z_i,dz_i,E_amp;
		ofstream temp_continuation;
		temp_continuation.open ("../Data_files/E_amp_onCauchy.dat");
		for (int j=1;j<= Memory->CauchyContour[0].size();j++) 
		{
			z_i=Memory->CauchyContour[0][j-1];
			dz_i=Memory->CauchyContour[1][j-1];
			E_amp=Memory->CauchyContour[2][j-1];
			temp_continuation << z_i << '\t' << dz_i << '\t' << E_amp << std::endl;
		}
		temp_continuation.close();
	}
	
	
	
	
	
	
