//
// Created by stkubr on 07.04.15.
//

#ifndef _DSEPP_INTEGRATOR_HPP_
#define _DSEPP_INTEGRATOR_HPP_

#include "../types.h"
#include "IntegrationNodes.hpp"

template <typename T_out, typename T_in> class C_Integrator_Line: public C_IntegrationNodes{
	protected:		

	std::function <T_out(T_in)> * Integrand;


	
	C_Integrator_Line(int _NumPoints, double _LimDown, double _LimUp, int _NumAps, Integrator_ID _id):
			C_IntegrationNodes(_NumPoints, _LimDown, _LimUp, _NumAps, _id) {
	}
	
	public:
	static C_Integrator_Line * createIntegrator(int _NumPoints, double _LimDown, double _LimUp, int _NumAps, Integrator_ID _id ){
		C_Integrator_Line * p = 0; 
		p = new C_Integrator_Line(_NumPoints, _LimDown, _LimUp, _NumAps, _id);
	    return p;
	}

	T_out getResult(std::function <T_out(T_in)> * Integrand)
	{
		T_out s, Rplus;
		int num_row,num_cols;

        Rplus=(*Integrand)(x[1]);
        num_row=Rplus.NumRows();
        num_cols=Rplus.NumCols();
        s.Resize(num_row,num_cols);
        s = w[1]*(Rplus);

		for (int j=2;j<=NumPoints;j++)
		{
			Rplus=(*Integrand)(x[j]);
			s += w[j]*(Rplus);
		}
		return s;
	}

};

#endif //_DSEPP_INTEGRATOR_HPP_

