//
// Created by stkubr on 21.04.15.
//

#ifndef KINEMATICS_1LOOP_H_
#define KINEMATICS_1LOOP_H_

#include "../types.h"

class C_Kinematics_1loop{
	private:
	DiracGamma Z;
	
	public:
	t_cmplxVector k,k_T, q,q_T, p,p_T, P, k_p,k_m,K;
	t_cmplx p_P, k_P, q_P, p2, k2, q2, P2, p2_T, k2_T, q2_T, N2_Factor;
	
	void SetVector_P(t_cmplx P_v)
	{
		 P.SetP4(0.0,0.0,0.0, (P_v) );
		 P2=P*P;
	}
	
	void SetVector_K(t_cmplx K_v)
	{
		 P.SetP4(0.0,0.0,0.0, (K_v) );
	}
	
	void SetVectors_p(t_cmplx z_ex, t_cmplx p_ex)
	{
		 p.SetP4(0.0, 0.0, sqrt(1.0-z_ex*z_ex)*p_ex, p_ex*z_ex); 
		 p_T=TransIn(p,P);
		 p_P=(p*P);
		 p2=p*p;
		 p2_T=p_T*p_T;
		 N2_Factor=(p2*P2 - p_P*p_P);
	}
	
	void SetVectors_k(double _zetta, t_cmplx x, t_cmplx y, t_cmplx z)
	{
		 k.SetP4(0.0, sqrt((1.0-z*z)*(1.0-y*y))*x, y*sqrt((1.0-z*z))*x, (x)*z );
		 k_T=TransIn(k,P);
		 k_P=(k*P);
		 k2=k*k;
		 k2_T=k_T*k_T;
	}
	
	void SetVectors_q()
	{
		 q=p-k;
		 q_T=TransIn(q,P);
		 q_P=(q*P);
		 q2=q*q;
		 q2_T=q_T*q_T;
	}
	
	void SetVestors_k_for_S(double _zetta, t_cmplxVector _k){
		k_p=_k+_zetta*P;
		k_m=_k+(_zetta-1.0)*P;
	}
	
	t_cmplxDirac TransIn(t_cmplxDirac _T, t_cmplxVector _P)
	{
		return _T - (_P*_T)*_P/(_P*_P);
	}
	
	t_cmplxTensor TransIn(t_cmplxVector _T, t_cmplxVector _P)
	{
		return _T - (_P*_T)*_P/(_P*_P);
	}
	
	void ShiftMomenta(double _zetta){
		q=k;
		
		q_T=TransIn(q,P);
		q_P=(q*P);
		q2=q*q;
		q2_T=q_T*q_T;
		
		k=p-k;
		
		k_T=TransIn(k,P);
		k_P=(k*P);
		k2=k*k;
		k2_T=k_T*k_T;
		
		k_p=k+_zetta*P;
		k_m=k+(_zetta-1.0)*P;	
	}
	
		
};

#endif /* KINEMATICS_1LOOP_H_ */