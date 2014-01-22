#pragma once

enum Integrator_ID { qgausleg_log_ID=0, qgausleg_lin_ID, qgauscheb_ID, qgausleg_sym_ID, qcauchyleg_lin_ID };

class C_Integrator: public C_AbstractClass{
	protected:
	double pi,LimUp,LimDown;
	int NumPoints,NumAps;
	t_dArray1D zz,w,x;
	Integrator_ID id;
	C_Integrator(int _NumPoints, double _LimDown, double _LimUp, int _NumAps, Integrator_ID _id){
		pi=3.14159265358979;
		NumPoints=_NumPoints;
		LimUp=_LimUp;
		LimDown=_LimDown;
		NumAps=_NumAps;
		id=_id;
		zz.resize(NumPoints+1);
		x.resize(NumPoints+1);
		w.resize(NumPoints+1);
		setNodes(id);
		//std::cout << "Creating Integrator with ID - " << id << std::endl;
	}
	
	//Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
	//arrays x[1..n] and w[1..n] of length n, containing the abscissas and weights of the Gauss-
	//Legendre n-point quadrature formula.
	void gauleg(double x1, double x2, t_dArray1D *x, t_dArray1D *w, int n)
	{
		int m,j,i;
		double z1,z,xm,xl,pp,p3,p2,p1;
		//High precision is a good idea for this routine.
		m=(n+1)/2;
		//The roots are symmetric in the interval, so we only have to ﬁnd half of them.
		xm=0.5*(x2+x1);
		xl=0.5*(x2-x1);
		for (i=1;i<=m;i++) {
			//Loop over the desired roots.
			z=cos(3.141592654*(i-0.25)/(n+0.5));
			//Starting with the above approximation to the ith root, we enter the main loop of reﬁnement by Newton’s method.
			do {
				p1=1.0;
				p2=0.0;
				for (j=1;j<=n;j++) {
					//Loop up the recurrence relation to get the
					p3=p2;
					//Legendre polynomial evaluated at z.
					p2=p1;
					p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
				}
			//p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
			//by a standard relation involving also p2, the polynomial of one lower order.
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
			
			//Newton’s method.
			} while (fabs(z-z1) > 3.0e-11);
			(*x)[i]=xm-xl*z;
			//std::cout << (*x)[i] << std::endl;
			//Scale the root to the desired interval,
			(*x)[n+1-i]=xm+xl*z;
			//and put in its symmetric counterpart.
			(*w)[i]=2.0*xl/((1.0-z*z)*pp*pp);
			//Compute the weight
			(*w)[n+1-i]=(*w)[i];
			//and its symmetric counterpart.
		}
	}

	//This routine returns arrays x[1..n] and w[1..n] of length n,
	//containing the abscissas and weights of the Gauss-Chebyshev n-point quadrature formula.
	void gaucheb(double x1, double x2, t_dArray1D *x, t_dArray1D *w, int n)
	{
		double zr;
		zr=0.5*(x2-x1);
		for (int i = 0; i <= n; i++)
		{
			(*x)[i]=0.0;
			(*w)[i]=0.0;
		}
		
		for (int i = 1; i <= n ; i++)
		{
			(*x)[i]=zr*cos((i/(n+1.0))*pi);
			(*w)[i]=(pi/(n+1.0))*sin(i/(n+1.0)*pi)*sin(i/(n+1.0)*pi);
		}
	}	
		
	void setNodes(Integrator_ID _id){
		double aa,bb,zm,zr,dz;
		switch (_id)
	    {
	        case qgausleg_log_ID:
	            this->gauleg(0.0,1.0,&zz,&w,NumPoints);
	            aa=log(LimDown);
	            bb=log(LimUp);		
				zm=(aa);
				zr=(bb-aa);
				for (int j=1;j<=NumPoints;j++) 
				{
					dz=zr*zz[j];
					x[j]=exp(zm+dz);
					w[j]=w[j]*exp(zm+dz)*zr;
				}
	            break;  
	              
	        case qgausleg_lin_ID:
				this->gauleg(0.0,1.0,&zz,&w,NumPoints);
				aa=(LimDown);
				bb=(LimUp);		
				zm=(aa);
				zr=(bb-aa);
				for (int j=1;j<=NumPoints;j++) 
				{
					dz=zr*zz[j];
					x[j]=(zm+dz - zr*zz[1]);
					w[j]=w[j]*zr;
				}
	            break;
	            
	        case qgauscheb_ID:
	            this->gaucheb(-1.0,1.0,&zz,&w,NumPoints);
	            for (int j=1;j<=NumPoints;j++) 
				{
					x[j]=zz[j];
					//w[j]=w[j]; //weight stays same
				}
	            break;
	            
	        case qgausleg_sym_ID:
	            this->gauleg(-1.0,1.0,&zz,&w,NumPoints);
	            aa=(LimDown);
				bb=(LimUp);		
				zm=0.5*(bb+aa);
				zr=0.5*(bb-aa);
				for (int j=1;j<=NumPoints;j++) 
				{
					dz=zr*zz[j];
					x[j]=(zm+dz);
					w[j]=w[j]*zr;
					//std::cout <<x[j] << "  " <<  w[j] << "  " << x[j]*w[j] <<std::endl;
					//cin.get();
				}
	            break;
	            
	        case qcauchyleg_lin_ID:
	            this->gauleg(0.0,1.0,&zz,&w,NumPoints);
	            aa=(LimDown);
				bb=(LimUp);		
				zm=(aa);
				zr=(bb-aa);
				for (int j=1;j<=NumPoints;j++) 
				{
					dz=zr*zz[j];
					x[j]=(zm+dz - zr*zz[1]);
					w[j]=w[j]*zr;
				}
	            break;
	            
	        default:
	            assert( false);
	    }
	}
	
	public:
	void getNodes(t_dArray1D * _x, t_dArray1D * _w){
		(*_x)=x;
		(*_w)=w;
	}
	
};

template <typename T_out, typename F, typename T_in> class C_Integrator_Line: public C_Integrator{
	protected:		
	T_out (C_Integrator_Line::*Integrator_ref)(T_out (F::*func) (T_in), F * obj_ref);
	
	C_Integrator_Line(int _NumPoints, double _LimDown, double _LimUp, int _NumAps, Integrator_ID _id):C_Integrator(_NumPoints, _LimDown, _LimUp, _NumAps, _id){
	}
	
	public:
	static C_Integrator_Line * createIntegrator(int _NumPoints, double _LimDown, double _LimUp, int _NumAps, Integrator_ID _id ){
		C_Integrator_Line * p = 0; 
		p = new C_Integrator_Line(_NumPoints, _LimDown, _LimUp, _NumAps, _id); 
	    /*switch (_id)
	    {
	        case qgausleg_log_ID:
				p = new C_Integrator_Line(_NumPoints, _LimDown, _LimUp, _NumAps, _id); 
				//p->Integrator_ref=&C_Integrator_Line::qgaus_log;
	            break;  
	            
	        case qgausleg_lin_ID:
				p = new C_Integrator_Line(_NumPoints, _LimDown, _LimUp, _NumAps, _id); 
	            //p->Integrator_ref=&C_Integrator_Line::qgaus_lin;
	            break;     
	            
	        case qgauscheb_ID:
				p = new C_Integrator_Line(_NumPoints, _LimDown, _LimUp, _NumAps, _id); 
	            //p->Integrator_ref=&C_Integrator_Line::qgauscheb;
	            break;        
	        default:
				std::cout << "Integrator_ID Error" << std::endl;
	            assert( false);
	    }*/
	    p->Integrator_ref=&C_Integrator_Line::qgaus;
	    return p;
	}
	
	T_out getResult(T_out (F::*func) (T_in), F * obj_ref){
		return (this->*Integrator_ref)(func,obj_ref);
	}
	
	T_out qgaus(T_out (F::*func) (T_in), F * obj_ref)
	{
		T_out s, Rplus;
		int num_row,num_cols;
		for (int j=1;j<=NumPoints;j++) 
		{
			Rplus=(obj_ref->*func)(x[j]);
			num_row=Rplus.NumRows();
			num_cols=Rplus.NumCols();
			s.Resize(num_row,num_cols);
			s += w[j]*(Rplus);
		}
		return s;
	}

};

template <typename T_out,typename T_contour ,typename T_in> class C_Integrator_Cauchy: public C_Integrator{
	protected:		
	T_out (C_Integrator_Cauchy::*Integrator_ref)(T_contour * Contour,int num_part ,T_in * Point);
	
	C_Integrator_Cauchy(int _NumPoints, double _LimDown, double _LimUp, int _NumAps, Integrator_ID _id):C_Integrator(_NumPoints, _LimDown, _LimUp, _NumAps, _id){}
	
	public:
	static C_Integrator_Cauchy * createIntegrator(int _NumPoints, double _LimDown, double _LimUp, int _NumAps, Integrator_ID _id ){
		C_Integrator_Cauchy * p; 
	    switch (_id)
	    {
	        case qcauchyleg_lin_ID:
				p = new C_Integrator_Cauchy(_NumPoints, _LimDown, _LimUp, _NumAps, _id); 
	            p->Integrator_ref=&C_Integrator_Cauchy::qcauchyleg_lin;
	            break;       
	        default:
				std::cout << "Integrator_ID Error" << std::endl;
	            assert( false);
	    }
	    return p;
	}
	
	T_out getResult(T_contour * Contour, int num_part ,T_in * Point){
		return (this->*Integrator_ref)(Contour,num_part,Point);
	}
	
	T_out qcauchyleg_lin(T_contour * Contour, int num_part,T_in * Point)
	{
		T_out result(NumAps+1),sumF(NumAps);
		T_in sumN;
		for (int i = 0; i < NumAps ; i++){sumF[i] = T_in(0.0,0.0);}
		sumN = T_in(0.0,0.0);
//#pragma omp parallel 
{
		//#pragma omp for
		for (int j=1;j<=NumPoints;j++) 
		{
			T_out F(NumAps);
			T_in N,temp;
			T_in z_i,dz_i,coordin;
			temp = w[j];
			z_i=(*Contour)[num_part][0][j-1];
			dz_i=(*Contour)[num_part][1][j-1];
			coordin=(*Point);
			for (int i = 0; i < NumAps ; i++){F[i]=conj(z_i-coordin)/norm(z_i-coordin)*dz_i*(*Contour)[num_part][i+2][j-1];}
			N=conj(z_i-coordin)/norm(z_i-coordin)*dz_i;
			#pragma omp critical
			{
			for (int i = 0; i < NumAps ; i++){sumF[i] += temp*F[i];}
			sumN += temp*N;
			}
			//std::cout << coordin << "  " << z_i << "  " << dz_i << "  " << sumN << "  " << sumF[0] << "  " << (*Contour)[num_part][2][j-1]  << std::endl;
		}			
}		
		result[0]=(sumN);
		for (int i = 0; i < NumAps ; i++){	result[i+1]=(sumF[i]);  }
		return result;
	}
};

