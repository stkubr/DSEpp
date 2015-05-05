/*
 * Support_functions.h
 *
 *  Created on: Jan 21, 2014
 *      Author: stkubr
 */

#ifndef EXTRA_FUNCTIONS_H_
#define EXTRA_FUNCTIONS_H_

#include "../types.h"

/// root finder lambda(P)=1
inline t_cmplx findState(std::function<t_cmplx(t_cmplx)> Func, double down_limit, double up_limit){
	t_cmplx start,end,mid;
	t_cmplx s_P,e_P,m_P,r_P;
	double accuracy=0.001;
	double eps=1.0;

	s_P=down_limit;e_P=up_limit;
	int ctr=0;
	while (eps>accuracy){
		m_P=(e_P+s_P)/2.0;
		mid=Func(m_P);
		eps = fabs(real(mid) - 1.0);
		std::cout << mid << " " << m_P << " " << eps << std::endl;
		if(fabs(imag(mid))>0.1 || ctr>30 || fabs(real(up_limit - m_P)/up_limit) < 0.01  ) { break;}
		if(real(mid)<1.0) e_P=m_P;
		if(real(mid)>1.0) s_P=m_P;
		ctr++;
	}
	r_P=m_P;
	return r_P;
}

inline double derivative(std::function<double(double)> & Func, double coordin, double delta){
	return (Func(coordin + delta) - Func(coordin - delta))/(2.0*delta);
}

/// Prints a line accross the screen
inline void PrintLine(char __c) {
	for (int i = 0; i < 80; i++)
		std::cout << __c;
	std::cout << std::endl;
}

/// Prints "h" spaces
inline void PrintSpaces(int h) {
	for (int i = 0; i < h; i++)
		std::cout << " ";
}

#endif /* EXTRA_FUNCTIONS_H_ */
