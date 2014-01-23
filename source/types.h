/*
 * types.hpp
 *
 *  Created on: Jan 21, 2014
 *      Author: stkubr
 */

#ifndef TYPES_HPP_
#define TYPES_HPP_

#include <complex>
#include "../qft++/include/matrix.h"
#include "../qft++/include/tensor.h"
#include "../qft++/include/relativistic-quantum-mechanics.h"
#include <vector>


typedef std::complex<double> t_cmplx;

typedef Vector4<complex<double> > t_cmplxVector;

typedef Matrix<double> t_dMatrix;
typedef Matrix<complex<double> > t_cmplxMatrix;

typedef Tensor<double> t_dTensor;
typedef Tensor<complex<double> > t_cmplxTensor;

typedef Matrix<Tensor<complex<double> > > t_cmplxDirac;

typedef std::vector<std::vector<t_cmplxTensor> > t_cmplx2DVecTen;
typedef Matrix< t_cmplxMatrix > t_cmplxMatrix2D;

typedef std::vector<double> t_dArray1D;
typedef std::vector<std::vector<double> > t_dArray2D;
typedef std::vector<std::vector<std::vector<double> > > t_dArray3D;

typedef std::vector<t_cmplx> t_cmplxArray1D;
typedef std::vector<std::vector<t_cmplx> > t_cmplxArray2D;
typedef std::vector<std::vector<std::vector<t_cmplx> > > t_cmplxArray3D;
typedef std::vector<t_cmplxArray3D> t_cmplxArray4D;

#endif /* TYPES_HPP_ */
