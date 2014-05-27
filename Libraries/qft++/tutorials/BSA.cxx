#include <math.h>
#include <complex>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <getopt.h>
#include <time.h>
#include <sstream>

#include "omp.h"
#include "/home/nucleus/stkubr/Documents/code/qft++/include/matrix.h"
#include "/home/nucleus/stkubr/Documents/code/qft++/include/tensor.h"
#include "/home/nucleus/stkubr/Documents/code/qft++/include/relativistic-quantum-mechanics.h"
#include "/home/nucleus/stkubr/Documents/code/DS_eq/DSE_prop_cont.cpp"

#include "/home/nucleus/stkubr/Documents/code/Support_functions/Support_functions.cpp"

#include "/home/nucleus/stkubr/Documents/code/BSA/BSA_func_descript.cpp"

//#define F Y5*F0*Y*Y/4.0// - i*(Y5*(P*Z))*F0
#define _  <<"   "<<

#define D_VALUE2 0.93
#define w2_VALUE2 0.4*0.4
#define pi 3.14159265
#define m2pi (0.138*0.138)
#define NUM_STEPS 30

using namespace std;

typedef complex<double> dcx;

//_____________________________________________________________________________

class BSA
{
	//--------------------------------------------------------------------------	
	private: 
	//Private
	double ksav,zsav,ysav,m0,w2,zet,D,LimKd,LimKu;
    dcx i,T,N,pion;
    // Private Flags
	bool sigma_flag,flag_grids;
	int sigma_counter,e,NumIntergation;
	
    //--------------------------------------------------------------------------	
	public: 
	//Public
	double (BSA::*nrfunc)(double,double,double,double); //pointer to Integrand function
	double p_v, P_v, z_v, F0, lambda_EV, norm_factor, AMP, grid_k[640], grid_z[640], NORM, P_EV, NORM_EV, pion_AMP,zz[640],w[640],zz5[64],w5[64];
	//Public Flags
	bool flag_sigma_buffer,normalize,flag_calc_EV,flag_calc_BSA,flag_BSA_buffer;
	
	int n,m; //grid counters
	int k_col,z_row;
	dcx k_pv, k_mv, buffer_sigma_p_s[4100], buffer_sigma_p_v[4100], buffer_sigma_m_s[4100], buffer_sigma_m_v[4100];
	
	// Tensor Objects
	Matrix<complex<double> > S_p,S_m; //Propagators
	Matrix<complex<double> > F;
	Matrix<complex<double> > I;
	
    Vector4<complex<double> > p,P,k,q,k_p,k_m;
    
    Matrix<double> Buffer_AMP;
	Matrix<double> AMP_m;
	
	DiracGamma Z,Y;
	DiracGamma5 Y5;
	MetricTensor g;
	
	// Side Objects
	CQuark quark_p;
	CQuark quark_m;
	
	
//Constructor
//------------------------------------------------------------------
	BSA():S_p(4,4),S_m(4,4),F(4,4),I(4,4),Buffer_AMP(NUM_STEPS+1,2),AMP_m(NUM_STEPS+1,2),quark_p(m0),quark_m(m0)
	{
		m0=0.005;
		
		flag_sigma_buffer=false; //by default turn on calculation of quark propagators 
		flag_calc_EV=false;
		flag_BSA_buffer=false;
		normalize=false;
		
		w2=w2_VALUE;
		D=D_VALUE;
		zet=0.5;
		
		NumIntergation=NUM_STEPS;
		
		gauleg(0.0,1.0,zz,w,NUM_STEPS);
		
		gauleg(0.0,1.0,zz5,w5,5);
		
		F0=1.0;
		norm_factor=1.0;
		
		LimKd=0.0025;
		LimKu=25.0;
		
		I=(1.0/4.0*(Z*Z));
		
		i=dcx(0.0,1.0); // the imaginary unit
		
		//Counters to Zero!
		sigma_counter=0;
		n=0;
		m=0;
		
		k_col=0;
		z_row=-1;
	}
	


void gauleg(double x1, double x2, double x[], double w[], int n);/*
//Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
//arrays x[1..n] and w[1..n] of length n, containing the abscissas and weights of the Gauss-
//Legendre n-point quadrature formula.
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
	x[i]=xm-xl*z;
	//Scale the root to the desired interval,
	x[n+1-i]=xm+xl*z;
	//and put in its symmetric counterpart.
	w[i]=2.0*xl/((1.0-z*z)*pp*pp);
	//Compute the weight
	w[n+1-i]=w[i];
	//and its symmetric counterpart.
}
}
	*/
//_____________________________________________________________________________
// Integrator Gauss-Legendre (n points) LOG scale
	double qgaus_n(double (BSA::*func) (double), double a, double b)
	{
	int j;
	double zr,zm,dz,s;
	s=0;
	double Rplus;
	
	//double aa=a, bb=b;	                       
	//double aa=sqrt(1.0)*pow(a,1.0/1.8), bb=sqrt(1.0)*pow(b,1.0/1.8);	          
	double aa=(log((a))), bb=(log((b)));	
		
	zm=(aa);
	zr=(bb-aa);
		
	for (j=1;j<=NumIntergation;j++) 
	{
		dz=zr*zz[j];
			
		Rplus=(this->*func)(zm+dz);
		//Rminus=(this->*func)(zm-dz);
			
		s += w[j]*(Rplus);
	}	
	
	return s *= zr;
	}

//_____________________________________________________________________________
// Integrator Gauss-Legendre (5 points) MIRRORED
double qgaus5(double (BSA::*func) (double), double a, double b)
	{
	int j;
	double zr,zm,dz,s;
	s=0;
	double Rplus, Rminus;
	

    static float zz5[]={0.0,0.1488743389,0.4333953941,
                       0.6794095682,0.8650633666,0.9739065285};
    static float w5[]={0.0,0.2955242247,0.2692667193,
                       0.2190863625,0.1494513491,0.0666713443};
	
	
	double aa=a, bb=b;	                       
	zm=0.5*(bb+aa);
	zr=0.5*(bb-aa);
		
	for (j=1;j<=5;j++) 
	{
		dz=zr*zz5[j];
			
		Rplus=(this->*func)(zm+dz);
		Rminus=(this->*func)(zm-dz);
			
		s += w5[j]*(Rplus+Rminus);
	}
			
	return s *= zr;
	}

//_____________________________________________________________________________
// Integrator Gauss-Legendre (5 points)
double qgaus5_redused(double (BSA::*func) (double), double a, double b)
	{
	int j;
	double zr,zm,dz,s;
	s=0;
	double Rplus;
	/*
    static float zz5[]={0.0,0.1488743389,0.4333953941,
                       0.6794095682,0.8650633666,0.9739065285};
    static float w5[]={0.0,0.2955242247,0.2692667193,
                       0.2190863625,0.1494513491,0.0666713443};
	*/
	double aa=a, bb=b;	                       
	                      	
	zm=(aa);
	zr=(bb-aa);
		
	for (j=1;j<=5;j++) 
	{
		dz=zr*zz5[j];
			
		Rplus=(this->*func)(zm+dz);
		//Rminus=(this->*func)(zm-dz);
			
		s += w5[j]*(Rplus);
	}
			
	return s *= zr;
	}

//_____________________________________________________________________________
// Multidimensional integration routine
double quad3d(double (BSA::*func)(double, double, double,double), double x1, double x2)
{
nrfunc=func;

	return qgaus_n(&BSA::f1,x1,x2);
}

double f1(double k)
{
ksav=exp(k);

if(flag_grids==true) {
	m=0;
	grid_k[n+1]=exp(k);
	//cout << grid_k[n+1] _ n+1 << endl;
	n++;
}

if (flag_BSA_buffer==true) 
{
	k_col++;
	F0=AMP_m(k_col,1);
	//cout << ksav _ F0 _ k_col  << endl;
}

if(normalize==true) {
	z_row=-1;
	k_col++;	
}

	return ksav*qgaus5(&BSA::f2,-1.0,1.0);
}

double f2(double z)
{
zsav=z;
sigma_flag=true;

if(flag_grids==true) {
	grid_z[m]=zsav;
	m++;
}

if(normalize==true) {
	z_row++;
	//cout << k_col _ z_row _ Buffer_AMP(k_col,z_row) << endl;
}
	return qgaus5(&BSA::f3,-1.0,1.0);
}

double f3(double y)
{
ysav=y;

	return 2.0*qgaus5_redused(&BSA::f4, 0.0,1.0);
}

double f4(double z_p)
{
	return (this->*nrfunc)(ksav, zsav, ysav, z_p);
}

// Calculation routine for dressing functions
//------------------------------------------------------------------
void calc_sigma()
{
		
	if (sigma_flag == true){ 
		
		if(flag_sigma_buffer==false) {
			quark_p.calc_valueAB_cx(k_pv);
			buffer_sigma_p_s[sigma_counter]=quark_p.sigma_s;
			buffer_sigma_p_v[sigma_counter]=quark_p.sigma_v;
			//quark_p.sigma_s=buffer_sigma_p_s[sigma_counter];
			//quark_p.sigma_v=buffer_sigma_p_v[sigma_counter];
			//cout << quark_p.sigma_v _ (k_pv) _ sigma_counter << endl;
		}
		else {
			quark_p.sigma_s=buffer_sigma_p_s[sigma_counter];
			quark_p.sigma_v=buffer_sigma_p_v[sigma_counter];
			//cout << quark_p.sigma_v _ (k_pv) _ sigma_counter << endl;
		}
				
	}
	
	if (sigma_flag == true){ 
		
		if(flag_sigma_buffer==false) {
			quark_m.calc_valueAB_cx(k_mv); 
			buffer_sigma_m_s[sigma_counter]=quark_m.sigma_s;
			buffer_sigma_m_v[sigma_counter]=quark_m.sigma_v;
			sigma_counter++;
			//cout << quark_m.sigma_v _ (k_mv) _ sigma_counter << endl;
		}
		else {
			quark_m.sigma_s=buffer_sigma_m_s[sigma_counter];
			quark_m.sigma_v=buffer_sigma_m_v[sigma_counter];
			sigma_counter++;
			//cout << quark_m.sigma_v _ (k_mv) _ sigma_counter << endl;
		}
				
		sigma_flag = false;
		S_p=(-1.0*i*(k_p*Z)*quark_p.sigma_v + I*quark_p.sigma_s);
		S_m=(-1.0*i*(k_m*Z)*quark_m.sigma_v + I*quark_m.sigma_s);
	}
}

// Just fake function to get integration grid
//------------------------------------------------------------------
double fake_func(double x, double z, double y, double z_p)
{
	return 0;
}

//_____________________________________________________________________________
// Integrand BSA 
double Integrand_F (double x, double z, double y, double z_p)
{		
	double U_in,U_ex,T_factor;
	dcx Kernel_factor;
	
	p.SetP4(0.0, 0.0, i*sqrt(p_v)*sqrt(1.0-z_p*z_p), i*sqrt(p_v)*z_p); 
    k.SetP4(0.0, i*sqrt(x)*sqrt(1.0-z*z)*sqrt(1.0-y*y), i*sqrt(x)*y*sqrt(1.0-z*z), i*sqrt(x)*z );
    
    k_p=k+zet*P;
	k_pv=(k_p*k_p);
	k_m=k+(zet-1)*P;
	k_mv=(k_m*k_m);
    
    calc_sigma();
    
 	q=k-p;
 	
 	F=Y5*F0;
 	
 	U_in=(16.0*z*z*z*z - 12.0*z*z +1.0)             + (4*z*z - 1)      + (1.0);
 	U_ex=(16.0*z_p*z_p*z_p*z_p - 12.0*z_p*z_p +1.0)*0  + (4*z_p*z_p - 1)*0  + (1.0); 	
 	
 	T_factor=-D/(3.0*2.0*pi*pi*w2*w2*w2)*(x)*sqrt(1.0-z*z)*sqrt(1.0-z_p*z_p)*(U_in)*(U_ex);
 	
 	Kernel_factor= (q*q)*exp( -1.0*(q*q)/w2 ) + ((w2*w2*w2*0.48*(1.0-1.0*exp(-1.0*(q*q)))/(q*q))/(0.5*log(exp(2.0)-1.0+(1.0+(q*q)/0.234/0.234)*(1.0+(q*q)/0.234/0.234))));
 	
 	
	T=T_factor*(( ( (Y5*Y*(S_p)* F *(S_m)) % Y) | (g-1.0*(q%q)/(q*q)) )*Kernel_factor ).Tr();
	return real(T);
}

// Calculation BSA at point (p^2,P)
//------------------------------------------------------------------
double calc_BSA (double a, double b)
{
	p_v=a;
	P_v=b;
	
	// set P vectors
	P.SetP4(0.0,0.0,0.0,-1.0*sqrt(P_v*P_v) ); 

	if(flag_calc_EV==true){
		
		flag_sigma_buffer=false; // turn on calculation of quark propagators on a grid
		flag_calc_BSA=false;
		
			AMP= quad3d(&BSA::Integrand_F,LimKd,LimKu);
		
			// prenorm for eigenvalue	
			F0=AMP/norm_factor;
			norm_factor=F0;
		
			sigma_counter=0; 
			flag_sigma_buffer=true; // now let "sigma_buffer" turn on
			//cout << "i=" << i << "   " << F0  << "   " << AMP  << endl;
	} 
	
	if(flag_calc_BSA==true) {
			//if(flag_sigma_buffer==true){cout << "flag_sigma_buffer is switched ON" << endl;} 
			//else {cout << "flag_sigma_buffer is switched OFF" << endl;}
		
			flag_calc_EV=false;
				
			F0= quad3d(&BSA::Integrand_F,LimKd,LimKu);
				
			sigma_counter=0; 
			k_col=0;
			flag_sigma_buffer=true; // now let "sigma_buffer" turn on
		}
	return F0;
}


// Calculation BSA array at p^2 axis
//------------------------------------------------------------------
double BSA_curve(double P)
{
	double p2,ff1,ff2,ff3,ff4,temp;
	set_grids();
	ff1=0;ff2=0;ff3=0;
	
	cout << "Computating BSA array... " << endl;
		
	flag_calc_BSA=true;
	F0=1.0;
	cout << calc_BSA(0.0,P) << endl;
	
	
#pragma omp parallel num_threads(NUM_STEPS) private(p2)
{	
	BSA bsa_copy_omp (*this);
	bsa_copy_omp.flag_calc_BSA=true;
		
	#pragma omp for	
	for (int i = 1; i <= n; i++)
	{
		p2=grid_k[i];
		Buffer_AMP(i,0)=p2;
		Buffer_AMP(i,1)=bsa_copy_omp.calc_BSA(p2,P);
		//cout << "Blah blah, I`m debug line!" _ Buffer_AMP(i,0) _ Buffer_AMP(i,1) << endl;
	}
}//end of pragma

	cout << Buffer_AMP << endl;
	PrintLine('-');
		
	flag_BSA_buffer=true;
	
	for (int kk = 0; kk < 1 ; kk++)
	{
		for (int i = 1; i <= n; i++)
		{
			AMP_m(i,0)=Buffer_AMP(i,0);
			AMP_m(i,1)=Buffer_AMP(i,1);
		}
		
#pragma omp parallel num_threads(NUM_STEPS) private(p2)
{//start of pragma
	BSA bsa_copy_omp (*this);
	bsa_copy_omp.flag_calc_BSA=true;
	#pragma omp for	
		for (int j = 1; j <= n; j++)
		{
			p2=grid_k[j];
			Buffer_AMP(j,0)=p2;
			Buffer_AMP(j,1)=bsa_copy_omp.calc_BSA(p2,P);
			//cout << "Blah blah, I`m debug line!" << endl;
		}
}//end of pragma

		/*
		temp=Buffer_AMP(0,1);
		for (int i = 0; i < n; i++)
		{
			Buffer_AMP(i,0)=Buffer_AMP(i,0);
			Buffer_AMP(i,1)=Buffer_AMP(i,1)/temp;
		}*/
		
		PrintLine('*');
		cout << Buffer_AMP << endl;
		PrintLine('-');
		cout  << endl;
		
		PrintLine('-');
		cout << AMP_m << endl;
		PrintLine('*');
		cout << endl << endl;
		
		for (int i = 1; i <= n; i++)
		{
			ff1+=(Buffer_AMP(i,1)*AMP_m(i,1));
			ff2+=(AMP_m(i,1)*AMP_m(i,1));
		}
	
	cout  << ff1 _ ff2 _ ff1/ff2 _ (Buffer_AMP(1,1)*AMP_m(1,1))/(AMP_m(1,1)*AMP_m(1,1)) << endl;
		ff1=0;
		ff2=0;
		ff3=0;
		
		temp=Buffer_AMP(1,1);
		for (int i = 1; i <= n; i++)
		{
			Buffer_AMP(i,0)=Buffer_AMP(i,0);
			Buffer_AMP(i,1)=Buffer_AMP(i,1)/1.0;
		}
	}
	
	flag_BSA_buffer=false;
	
	ofstream UNnorm_BSA_grid;
	
	UNnorm_BSA_grid.open ("/home/nucleus/stkubr/Documents/code/qft++/tutorials/UNnorm_BSA_grid.dat");
		for (int i = 1; i <= n; i++)
		{
			UNnorm_BSA_grid << sqrt(Buffer_AMP(i,0)) << " " << Buffer_AMP(i,1) << " ;" << endl;
		}
		
	UNnorm_BSA_grid.close();
	
		
	PrintLine('-');
	cout << AMP_m << endl;
	PrintLine('-');
	cout << endl << endl;
	
	for (int i = 1; i <= n; i++)
		{
			ff1+=(Buffer_AMP(i,1)*AMP_m(i,1));
			ff2+=(AMP_m(i,1)*AMP_m(i,1));
			ff3+=Buffer_AMP(i,1)*(Buffer_AMP(i,0)-Buffer_AMP(i-1,0));
		}
		
	lambda_EV=ff1/ff2;
	cout  << ff1 _ ff2 _ ff3 _ lambda_EV _ (Buffer_AMP(1,1)*AMP_m(1,1))/(AMP_m(1,1)*AMP_m(1,1)) << endl;
	
	return lambda_EV;
}


// Eigenvalue curve calculation routine and writing into "EigenValue.dat"
//------------------------------------------------------------------
void set_dom_EV_curve(double scale)
{
	double dP,EV,Diff_EV,P2;
	double Pu=0.15;
	double Pd=0.12;
	dP=(Pu-Pd)/scale;
	P2=Pd;
	
	flag_calc_EV=true;
	
	ofstream EigenValue;
	EigenValue.open ("EigenValue.dat");
	for (int i = 0; i <= scale; i++)
	{
		EV=calc_BSA(0,P2);
		Diff_EV=0.5/P2*( log(calc_BSA(0,P2+dP))-log(calc_BSA(0,P2)) )/dP;
		EigenValue << P2 _ EV _ Diff_EV  << endl;
		cout << P2 _ EV _ Diff_EV  << endl;
		P2+=dP;
	}
	EigenValue.close();
}

// Eigenvalue curve root-finding and derivative at root point
//------------------------------------------------------------------
void set_dom_EV(double down_P, double up_P)
{
	double up,down,mid;
	double BSAup,BSAdown,BSAmid;
	double res;
	double EPS_root=0.0001;
	
	res=1.0;	
	up=up_P;
	down=down_P;
	mid=(up_P+down_P)/2.0;
	
	cout << endl << "Eigenvalue Calculation" << endl;
	
	while (abs(res)>EPS_root)
	{
		
		BSAup=BSA_curve(up)-1.0;
		cout << up << "   " << "UP point = " << BSAup << endl;
		BSAdown=BSA_curve(down)-1.0;
		cout << down << "   " << "DOWN point = " << BSAdown << endl;
		BSAmid=BSA_curve(mid)-1.0;
		cout << mid << "   " << "MID point = " << BSAmid << endl;
		
		res=BSAmid;
		P_EV=mid;
		cout << mid _ res+1.0 <<  endl;
		if(BSAup*BSAmid > 0) {down=down; up=mid; mid=(up+down)/2.0;}
		if(BSAdown*BSAmid > 0) {up=up; down=mid; mid=(up+down)/2.0;}
	}
	
	NORM_EV=deriv_EV(P_EV);
	cout << endl << P_EV _ res+1.0 _ NORM_EV << endl;
}

// Derivative at orbitary point
//------------------------------------------------------------------
double deriv_EV(double P_value)
{
	double h=0.01;
	double N_v;
	
	cout << "Derivative Calculation" << endl;
	
	N_v=0.5/P_value*(-log(BSA_curve(P_value+2.0*h)) + 8.0*log(BSA_curve(P_value+h)) - 8.0*log(BSA_curve(P_value-h)) + log(BSA_curve(P_value-2.0*h)) )/(12.0*h);
	return N_v;
}

// Get integration grids (k,z)
//------------------------------------------------------------------
void set_grids()
{
	flag_grids=true;
	quad3d(&BSA::fake_func,LimKd,LimKu);
	flag_grids=false;	
}

// Set BSA on grid to "Buffer_AMP" and write it to file "UNnorm_BSA_grid.dat" (code parallelized)
//------------------------------------------------------------------
void set_BSA_on_grid (double P)
{
	double k,z;
	set_grids();
	
	cout << "Computating BSA on grid... " << endl;
	
	BSA bsa_copy_grid;
	
#pragma omp parallel num_threads(10) private(k,z,bsa_copy_grid)
{
	bsa_copy_grid.flag_calc_BSA=true;
	#pragma omp for
	for (int i = 0; i < n ; i++)
	{
		k=grid_k[i];
		for (int j = 0; j < m ; j++)
		{
			z=grid_z[j];
			Buffer_AMP(i,j)=bsa_copy_grid.calc_BSA(k,P);
		}
	}
}//end of pragma
write_BSA();
}

// Write BSA to file "UNnorm_BSA_grid.dat"
//------------------------------------------------------------------
void write_BSA()
{
	cout << Buffer_AMP << endl << endl << endl;
	ofstream UNnorm_BSA_grid;
	
	UNnorm_BSA_grid.open ("/home/nucleus/stkubr/Documents/code/qft++/tutorials/UNnorm_BSA_grid.dat");
	for (int i = 0; i < n ; i++)
	{
		for (int j = 0; j < m ; j++)
		{
			 UNnorm_BSA_grid << Buffer_AMP(i,j) << " ";
		}
		UNnorm_BSA_grid <<";"<< endl;
	}
	UNnorm_BSA_grid.close();
}

// Read BSA from file "UNnorm_BSA_grid.dat" and set it to "Buffer_AMP"
//------------------------------------------------------------------
void get_BSA_on_grid()
{
	set_grids();
	string number;
	int i=0;
	//int j=0;
	double value;
	char c;
	
	ifstream UNnorm_BSA_grid("/home/nucleus/stkubr/Documents/code/qft++/tutorials/UNnorm_BSA_grid.dat");
	if (UNnorm_BSA_grid.is_open())
	{
		while ( UNnorm_BSA_grid.good() )
		{
			c=UNnorm_BSA_grid.get();
			if(c!=' ' && c!=';'){number += c;}
			if(c==' ')
			{
				value=atof( number.c_str() );
				AMP_m(i,1)=value; 
				cout << number _ i  << endl;
				number.clear();
			}
			if(c==';')
			{ 
				//cout << number _ i _ i << endl;
				number.clear();
				i++; 
			}
		}
		cout << AMP_m << endl;
	}
	else 
	{
		cout << "Cant open file!" << endl;
	}
	//cout << endl << Buffer_AMP << endl << endl;
}
/*
double Integrand_norm_F(double t, double z, double y)
{
	double x; 
	x=t;
	
	k.SetP4(0, i*sqrt(x)*sqrt(1.0-z*z)*sqrt(1.0-y*y), i*sqrt(x)*y*sqrt(1.0-z*z), i*sqrt(x)*z );
	
	k_p=k+zet*P;
	k_pv=(k_p*k_p);
	k_m=k+(zet-1)*P;
	k_mv=(k_m*k_m);

 	calc_sigma();
	
	N = 1.0/(64.0*1.23*pi*pi*pi)*x*(k*k)*sqrt(1-z*z)*( Y5*Buffer_AMP(k_col,z_row)*S_p*Y5*Buffer_AMP(k_col,z_row)*S_m ).Tr();
	
	return real(N);
}


// Normalization BSA Procedure
void calc_norm_BSA(double Plocal)
{
	P_v=Plocal;
	normalize=true;
	P.SetP4(0.0,0.0,0.0,-1.0*sqrt(P_v*P_v) );
	flag_sigma_buffer=false;
	k_col=-1;
	z_row=-1;
	NORM=quad3d(&BSA::Integrand_norm_F,LimKd,LimKu);
	
	cout << 1/sqrt(NORM) << endl;
}

double Integrand_pion(double t, double z, double y)
{
	double x; 
	x=t;
	
	k.SetP4(0, i*sqrt(x)*sqrt(1.0-z*z)*sqrt(1.0-y*y), i*sqrt(x)*y*sqrt(1.0-z*z), i*sqrt(x)*z );
	
	k_p=k+zet*P;
	k_pv=(k_p*k_p);
	k_m=k+(zet-1)*P;
	k_mv=(k_m*k_m);

 	calc_sigma();
	pion=-3.0/64.0/pi/pi/pi/P_v/P_v*x*(k*k)*sqrt(1-z*z)*( i*Y5*Buffer_F(k_col,z_row)*S_p*(Y*P)*Y5*S_m ).Tr();
	
	return real(pion);
}


// Pion Decay Constant Calculation
void calc_pion_decay(double Plocal)
{
	P_v=Plocal;
	normalize=true;
	k_col=-1;
	z_row=-1;
	P.SetP4(0,0,0,-1.0*sqrt(P_v*P_v) );
	flag_sigma_buffer=false;
	
	pion_AMP=quad3d(&BSA::Integrand_pion,LimKd,LimKu);
	
	cout << pion_AMP << endl;
}

*/

};
//_____________________________________________________________________________




//_____________________________________________________________________________
/// main()
int main(int __argc,char *__argv[]){
	// Time initialization
	double tstart, tstop, ttime;
	tstart = (double)clock()/CLOCKS_PER_SEC;
	
	BSA bsa;
	//cout << bsa.calc_BSA(0,0.139,0) << endl;
	//bsa.set_dom_EV(0.11, 0.14);
	bsa.BSA_curve(0.121);
	//bsa.BSA_curve(0.139);
	//bsa.set_dom_EV_curve(2);
	//bsa.get_BSA_on_grid();
	

	
	//BSA bsa_m;
	//bsa_m.set_BSA_on_grid(0.139);
	
	/*
	BSA bsa_1;
	bsa_1.get_BSA_on_grid();
	bsa_1.calc_norm_BSA(0.139);
	BSA bsa_2;
	bsa_2.get_BSA_on_grid();
	bsa_2.calc_pion_decay(0.139);
	cout << 1/sqrt(bsa_1.NORM)*bsa_2.pion_AMP << endl;
	*/
	
	//bsa.get_dom_EV_curve(30);
	//bsa.set_dom_EV(0.12, 0.14);
	//bsa.set_BSA_on_grid(0.1318);
	//bsa.get_BSA_on_grid();
	//bsa.omp_test();
	//bsa.calc_BSA(0,0.135,0);
	//bsa.dom_EV(30);

// Time of the calculation.
tstop = (double)clock()/CLOCKS_PER_SEC;
ttime= tstop-tstart;
cout << "\n\n\n" << "calculation time=" _ ttime/8.0 _ "seconds" << endl;

 return 0;
}
