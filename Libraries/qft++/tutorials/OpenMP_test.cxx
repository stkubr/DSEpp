#include <math.h>
#include <complex>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <omp.h>
#include "/home/nucleus/stkubr/Documents/code/qft++/include/matrix.h"
#include "/home/nucleus/stkubr/Documents/code/qft++/include/tensor.h"
#include "/home/nucleus/stkubr/Documents/code/qft++/include/relativistic-quantum-mechanics.h"


typedef complex<double> dcx;

using namespace std;

double f1(double a)
{
	return a;
}



class NiceClass 
{
	public:
	int n;
	double result[3],ksav,zsav;
	double grid1[3];
	double grid2[3];
	Matrix<dcx> M,S;
	MetricTensor g;
	Vector4<double> P;
	DiracGamma Z,Y;
	
	double (NiceClass::*nrfunc)(double,double,double);

//constructor
NiceClass(void):M(3,3),S(4,4)
{
	n=3;
}

void SetGrids()
{
	for(int i = 0; i < n; i++) {grid1[i] = 1.0; grid2[i] = 1.0*i;}
}


dcx NiceFunc (double a, double b)
{
	dcx res;
	P.SetP4(a,b,0,0);
	res=(P*P);
	return res;
}


// Integrator Gauss-Legendre (5 points)
double qgaus5(double (NiceClass::*func) (double), double a, double b)
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
// Multidimensional integration routine
double quad3d(double (NiceClass::*func)(double, double, double), double x1, double x2)
{
nrfunc=func;

return qgaus5(&NiceClass::f1,x1,x2);
}

double f1(double k)
{
ksav=k;

return qgaus5(&NiceClass::f2,-1,1);
}

double f2(double z)
{
zsav=z;

return qgaus5(&NiceClass::f3,-1,1);
}

double f3(double y)
{
	
return (this->*nrfunc)(ksav,zsav,y);
}



void obj_call()
{
	NiceClass Initial_object;
	for (int i = 0; i < 2; i++)
	{
		NiceClass one_of_obj = Initial_object;
		one_of_obj.calc_func();
	}
}



void calc_func()
{
	double g1,g2;
	int i,k,q;
	Vector4<double> P;
	SetGrids();
	
	
	NiceClass one_of_obj;
	
#pragma omp parallel num_threads(10) private(g1,g2,i,k,P,one_of_obj)
{
	one_of_obj.SetGrids();
	#pragma omp for
	for (i = 0; i < 3; i++)
	{
		q=omp_get_num_threads();
		g1=one_of_obj.grid1[i];
		for (k = 0; k < 3; k++)
		{
			g2=one_of_obj.grid2[k];
			M(i,k)=one_of_obj.NiceFunc(g1,g2);
		}
	}
}//end pragma 

cout << M << "   " << q << endl;
}

};

void func ()
{
	for (int i = 0; i < omp_get_num_threads(); i++)
	{
		
	}
	
}


int main(int __argc,char *__argv[]){

	NiceClass object;
	object.calc_func();
	cout << object.M << endl;

 return 0;
}
