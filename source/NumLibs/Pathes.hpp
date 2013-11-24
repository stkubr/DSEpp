#pragma once

namespace Geometry{



class C_Path{
public:
	t_cmplxArray1D getPathOnVector(t_cmplxArray1D SamplePoints ) {
			int NumSamplePoints;
			NumSamplePoints = SamplePoints.size();
			t_cmplxArray1D PathOnSamlePoints(NumSamplePoints,0);

			for (int i = 0; i < NumSamplePoints ; i++){
				PathOnSamlePoints[i]=getPathAt(SamplePoints[i]);
			}

			return PathOnSamlePoints;
		}

	virtual	t_cmplx getPathAt(t_cmplx t_paramtr){assert(false); return 0.0;}
};

class C_Parabola: public C_Path{
private:
	t_cmplx ParabolaApex;
	
public:
	C_Parabola(t_cmplx __ParabolaApex){
		ParabolaApex=__ParabolaApex;
	}
		
	t_cmplx getPathAt(t_cmplx t_paramtr){
		// f(t)=t*t + i*2.0*t*Apex + Apex*Apex
		return t_paramtr*t_paramtr + 2.0*t_paramtr*ParabolaApex + ParabolaApex*ParabolaApex;
	}
	
};

class C_Line: public C_Path{
private:
	t_cmplx k_coeff, b_coeff;

public:
	C_Line(t_cmplx __k_coeff, t_cmplx __b_coeff){
		k_coeff=__k_coeff;
		b_coeff=__b_coeff;
	}

	t_cmplx getPathAt(t_cmplx t_paramtr){
		// f(t)=t*k_coeff + b_coeff
		return t_paramtr*k_coeff + b_coeff;
	}
};



}
