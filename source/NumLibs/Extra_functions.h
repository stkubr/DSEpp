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

#endif /* SUPPORT_FUNCTIONS_H_ */
