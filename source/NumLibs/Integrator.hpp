//
// Created by stkubr on 07.04.15.
//

#ifndef _DSEPP_INTEGRATOR_HPP_
#define _DSEPP_INTEGRATOR_HPP_

#include "../types.h"
#include "IntegrationNodes.hpp"

template <typename T_out, typename T_in> class C_Integrator_Line: public C_IntegrationNodes{
	protected:
	C_Integrator_Line(int _NumPoints, double _LimDown, double _LimUp, int _NumAps, Integrator_ID _id):
			C_IntegrationNodes(_NumPoints, _LimDown, _LimUp, _NumAps, _id) {
	}
	
	public:
	static C_Integrator_Line * createIntegrator(int _NumPoints, double _LimDown, double _LimUp, int _NumAps, Integrator_ID _id ){
		C_Integrator_Line * p = 0; 
		p = new C_Integrator_Line(_NumPoints, _LimDown, _LimUp, _NumAps, _id);
	    return p;
	}

	T_out getResult(std::function <T_out(T_in)> * func_to_int,  int numRows, int numCols) {
		T_out sum, result;
		sum.Resize(numRows,numCols);
		result.Resize(numRows,numCols);
		for (int j=1;j<=NumPoints;j++) {
			result=(*func_to_int)(x[j]);
			sum += w[j]*(result);
		}
		return sum;
	}

};

#endif //_DSEPP_INTEGRATOR_HPP_

