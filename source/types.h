/*
 * types.hpp
 *
 *  Created on: Jan 21, 2014
 *      Author: stkubr
 */



#ifndef TYPES_HPP_
#define TYPES_HPP_

#include <complex>
#include "matrix.h"
#include "tensor.h"
#include "relativistic-quantum-mechanics.h"
#include <vector>

#ifndef _NUM_THREADS
#define _NUM_THREADS omp_get_max_threads()
#endif

typedef std::complex<double> t_cmplx;

/// 4-vector in Euclidean space
typedef Vector4<complex<double> > t_cmplxVector;

/// Matrix (if in Dirac space then (4,4) )
typedef Matrix<double> t_dMatrix;

/// Matrix (if in Dirac space then (4,4) )
typedef Matrix<complex<double> > t_cmplxMatrix;

/// Tensor in Euclidean space with arbitrary rank
typedef Tensor<double> t_dTensor;

/// Tensor in Euclidean space with arbitrary rank
typedef Tensor<complex<double> > t_cmplxTensor;

/// Composite object which is Dirac matrix and Lorentz tensor
typedef Matrix<Tensor<complex<double> > > t_cmplxDirac;

/// Two-quark scattering Kernel data format;
/// it has two incoming dirac indexes and two outgoing - /f$ K_{(ts)(ru)} /f$
/// so it is (4,4)(4,4) matrix of rank 4
typedef Matrix< t_cmplxMatrix > t_cmplxMatrix2D;

typedef std::vector<double> t_dArray1D;
typedef std::vector<t_dArray1D> t_dArray2D;
typedef std::vector<t_dArray2D> t_dArray3D;

typedef std::vector<t_cmplx> t_cmplxArray1D;
typedef std::vector<t_cmplxArray1D> t_cmplxArray2D;
typedef std::vector<t_cmplxArray2D> t_cmplxArray3D;
typedef std::vector<t_cmplxArray3D> t_cmplxArray4D;

#endif /* TYPES_HPP_ */
