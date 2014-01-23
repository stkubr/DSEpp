#pragma once

#include <memory>
#include <vector>

namespace Geometry {

class C_Path {
public:
	t_cmplxArray1D getPathOnVector(t_cmplxArray1D SamplePoints) {
		int NumSamplePoints;
		NumSamplePoints = SamplePoints.size();
		t_cmplxArray1D PathOnSamlePoints(NumSamplePoints, 0);

		for (int i = 0; i < NumSamplePoints; i++) {
			PathOnSamlePoints[i] = getPathAt(SamplePoints[i]);
		}

		return PathOnSamlePoints;
	}

	virtual t_cmplx getPathAt(t_cmplx t_paramtr) {
		assert(false);
		return 0.0;
	}
	virtual t_cmplx getDerivativePathAt(t_cmplx t_paramtr) {
		assert(false);
		return 0.0;
	}
};

class C_Parabola: public C_Path {
private:
	t_cmplx ParabolaApex;

public:
	C_Parabola(t_cmplx __ParabolaApex) {
		ParabolaApex = __ParabolaApex;
	}

	// f(t)=t*t + 2.0*t*Apex + Apex*Apex
	t_cmplx getPathAt(t_cmplx t_paramtr) {
		return t_paramtr * t_paramtr + 2.0 * t_paramtr * ParabolaApex
				+ ParabolaApex * ParabolaApex;
	}

	// f^{\prime}(t)=2*t + 2.0*Apex
	t_cmplx getDerivativePathAt(t_cmplx t_paramtr) {
		return 2.0 * t_paramtr + 2.0 * t_paramtr * ParabolaApex;
	}

};

class C_Line: public C_Path {
private:
	t_cmplx k_coeff, b_coeff;

public:
	C_Line(t_cmplx __k_coeff, t_cmplx __b_coeff) {
		k_coeff = __k_coeff;
		b_coeff = __b_coeff;
	}

	// f(t)=t*k_coeff + b_coeff
	t_cmplx getPathAt(t_cmplx t_paramtr) {
		return t_paramtr * k_coeff + b_coeff;
	}

	// f^{\prime}(t)=k_coeff
	t_cmplx getDerivativePathAt(t_cmplx t_paramtr) {
		return k_coeff;
	}
};

class C_ParabolaContour {
private:
	std::shared_ptr<C_Path>  parabola;
	std::shared_ptr<C_Path>  line;
	//
	t_cmplxArray2D ContourPath;

public:
	C_ParabolaContour(t_cmplx __ParabolaApex, t_cmplx __k_coeff, t_cmplx __b_coeff){
		ContourPath.resize(2);
		parabola = std::make_shared<C_Parabola> (__ParabolaApex);
		line = std::make_shared<C_Line> (__k_coeff, __b_coeff);
	}

	void setParabolaContour(const t_dArray1D& p_parabola,
							const t_dArray1D& w_parabola,
							const t_dArray1D& p_line,
							const t_dArray1D& w_line){

		/*
		 * Normally Cauchy contour integration goes counter-clockwise,
		 * so the weights of some pieces (going clock-wise) of contour has to have
		 * "-" sign to contribute to total contour integral with a right sign.
		 *
		 * Also the contour is symmetric - so the lower pars are just
		 * complex conjugation of upper parts.
		 */

		// Upper parabola of contour (clock-wise)
		int i=0;
		for_each(p_parabola.begin(), p_parabola.end(), [&](t_cmplx t_param) {
			ContourPath[0].push_back(parabola->getPathAt(t_param));
			ContourPath[1].push_back(parabola->getDerivativePathAt(t_param)*-1.0*w_parabola[i]);
			i++;
		});

		// Upper line of contour
		i=0;
		for_each(p_line.begin(), p_line.end(), [&](t_cmplx t_param) {
			ContourPath[0].push_back(line->getPathAt(t_param));
			ContourPath[1].push_back(line->getDerivativePathAt(t_param)*w_line[i]);
			i++;
		});

		// Lower line of contour (clock-wise)
		i=0;
		for_each(p_line.begin(), p_line.end(), [&](t_cmplx t_param) {
			ContourPath[0].push_back(conj(line->getPathAt(t_param)));
			ContourPath[1].push_back(conj(line->getDerivativePathAt(t_param))*-1.0*w_line[i]);
			i++;
		});

		// Lower parabola of contour
		i=0;
		for_each(p_parabola.begin(), p_parabola.end(), [&](t_cmplx t_param) {
			ContourPath[0].push_back(conj(parabola->getPathAt(t_param)));
			ContourPath[1].push_back(conj(parabola->getDerivativePathAt(t_param))*w_parabola[i]);
			i++;
		});

	}

};

}
