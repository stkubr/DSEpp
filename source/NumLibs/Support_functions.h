/*
 * Support_functions.h
 *
 *  Created on: Jan 21, 2014
 *      Author: stkubr
 */

#ifndef SUPPORT_FUNCTIONS_H_
#define SUPPORT_FUNCTIONS_H_


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

inline t_cmplx Cheb_polinoms(t_cmplx x, int order) {
	t_cmplx U[order + 10];

	U[0] = 1.0;
	U[1] = 2.0 * x;
	U[2] = 4.0 * x * x - 1.0;

	if (order == 0) {
		return U[0];
	}

	if (order == 1) {
		return U[1];
	}

	if (order == 2) {
		return U[2];
	}

	if (order > 2) {
		for (int i = 3; i <= order; i++) {
			U[i] = 2.0 * x * U[i - 1] - U[i - 2];
		}
		return U[order];
	}

	std::cout << "Chebys Error!" << std::endl;
	return 0;
}

inline double Get_Time() {
	return (double) clock() / CLOCKS_PER_SEC;
}

inline void DebugLine(string _word) {
//#if DEBUG_MODE==1
	std::cout << " Debug Line at place - " << _word << std::endl;
//#endif
}

#endif /* SUPPORT_FUNCTIONS_H_ */
