#include <math.h>
#include <complex>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <getopt.h>
#include <time.h>
#include "/home/nucleus/stkubr/Documents/code/qft++/include/tensor.h"
#include "/home/nucleus/stkubr/Documents/code/qft++/include/relativistic-quantum-mechanics.h"
#include "/home/nucleus/stkubr/Documents/code/Support_functions/Support_functions.cpp"

#define I 1.0/4.0*(Z*Z)

using namespace std;
//_____________________________________________________________________________

/// Prints program usage to the screen
void PrintUsage();

//_____________________________________________________________________________

typedef complex<double> dcx;
typedef Matrix<Tensor<complex<double> > > dcxDirac;

int main(int __argc,char *__argv[]){

  /*__________________________Parse the Command Line_________________________*/
  int c;
  //extern char* optarg;
  //extern int optind;
  
  while((c = getopt(__argc,__argv,"h")) != -1){
    switch(c){
    case 'h': // help option
      PrintUsage();
      return EXIT_SUCCESS;
      break;
    default:
      break;
    }
  }  
	double Init_Time;
	Init_Time=Get_Time();
  /*_____________________________Creating Tensors____________________________*/
  // Tensor is a template class, here we'll work with double's and 
  // complex<double>'s...but you could use any data type that defines the 
  // proper operators. The constructor argument is the rank.
  Tensor<double> x1(1),x2(2),x3(3);
  Tensor<complex<double> > y1(1),y(2),B(2),A_1(3),A_2(2);
  Vector4<double> v;
  Matrix <complex<double> > S1(4,4),S2(4,4),F(4,4);
  Matrix<Tensor<complex<double> > > A[4];
  LeviCivitaTensor eps; // totally anti-symetric 4th rank tensor
  MetricTensor g; // g_{mu,nu} 
  DiracGamma Z,Y,Q;
  DiracGamma5 Y5;
  DiracSigma SIG;
  dcxDirac T;
  Vector4<complex<double> > P,k,k_t;
  /*_____________________________Setting Tensors_____________________________*/
  // The metric and Levi-Civita tensors are set when they're created. To set
  // elements of our other tensors we have several options. We can either set
  // them using other Tensor's (using the = operator), or by direct access to
  // the elements.

  PrintLine(':');
  cout << "We've set up the following tensors: " << endl;
  cout << "tensors storing double's:" << endl;
  // element access
  dcx i(0,1);
  dcx P2;
  for(int e = 0; e < 4; e++) {
    x1(e) = e;
    x2(e,e) = 1.0; // make it diagnol
  }
  cout << "->x1: " << x1 << endl;
  cout << "->x2: " << x2.Mag2() << endl;
 P.SetP4(0,0,0,1.0);
 k.SetP4(1.0,1.0,1.0,1.0);
 
  // just make sure these don't throw errors
  // x3.Symmetric();
  //x3.AntiSymmetric();

  cout << "->metric: g^{mu,nu}: " << g << endl;
  

	struct Amp
	{
		Matrix<Tensor<complex<double> > > W[5];
	} amp;
	

  

  cout << "tensors storing complex double's:" << endl;
  y1(0) = complex<double>(0.,1.);
  y1(3) = complex<double>(0.,1.);
  cout << "->y1: " << y1 << endl;

  /*_____________________________Basic Operations____________________________*/
  PrintLine(':');

  cout << "Tensor contractions:" << endl;
 

  // To contract the last index of x with the 1st index of y, just do x*y
  cout << "contraction of 1 index: " << endl;
  cout << "x1_{mu}x1^{mu}:\n->" << x1*x1 << endl;
  cout << "x1_{mu}x2^{mu,nu}:\n->" << x1 * x2 << endl;
  cout << "x1_{mu}x3^{mu,nu,rho}:\n->" << x1 * x3 << endl;
  cout << "x1_{rho}x3^{mu,nu,rho}:\n->" << x3 * x1 << endl;
  cout << (x3.Permute(2,3)) * x1 << endl;
  cout << "sqrt(v^{mu}v_{mu}):\n->" << P*P << endl;
  cout << "DIRAC TEST 1  " << endl  << endl;
	T=SIG;
	cout << T << endl;
	//cout << T.Size() << endl;

	cout << "Ololo" << endl;

  // To contract the last n indicies of x with the 1st n indicies of y, do
  // x.Contract(y,n). If n is the rank of either x or y, then you can just do
  // (x|y)...a tensor inner product.
  cout << "contraction of multiple indicies:" << endl;
  cout << "x2_{mu,nu}x3^{mu,nu,rho}:\n->" << (x2|x3) << endl;
  cout << "x1^{mu}x2^{nu,rho}x3_{mu,nu,rho}:\n->" << ((x1%x2)|x3) << endl;
  cout << "(x1^{mu}x2^{nu,rho} + 2*x3^{mu,nu,rho})x3_{mu,nu,rho}:\n->"
       << (((x1%x2) + 2*x3)|x3) << endl;
  cout << "epsilon^{mu,nu,rho,pi} x2_{mu,nu}:\n->" << (eps|x2) << endl;
  cout << "etc...as complicated as you want, it's still easy in the code."
       << endl;


  /*_________________________Special 4-Vector Methods________________________*/
  /*PrintLine(':');
  cout << "There are also a number of special methods defined for 4-vectors:"
       << endl;
  cout << "4-vector v:->" << v << endl;
  cout << "invariant mass:-> " << v.Mass() << endl;
  cout << "beta:-> " << v.Beta() << endl;
  cout << "cos(theta):-> " << v.CosTheta() << endl;
  cout << "see the Vector4 class documentation for a complete list." << endl;
  PrintLine(':');

  return EXIT_SUCCESS;*/
}
//_____________________________________________________________________________

void PrintUsage(){ 
  cout << "Usage: tensor " << endl;
  cout << "This executable provides a number of example usages of the tensor "
       << "package. Run\nthe executable to see what's being done, then look at"
       << " the source file to see \nhow it's done in the code."
       << endl;
}
//_____________________________________________________________________________
