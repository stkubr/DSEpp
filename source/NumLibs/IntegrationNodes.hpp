//
// Created by stkubr on 07.04.15.
//

#ifndef _DSEPP_INTEGRATIONNODES_HPP_
#define _DSEPP_INTEGRATIONNODES_HPP_

#include "../types.h"

enum Integrator_ID { qgausleg_log_ID=0, qgausleg_lin_ID, qgauscheb_ID, qgausleg_sym_ID };

class C_IntegrationNodes{
protected:
    const double pi=3.14159265358979;
    double LimUp,LimDown;
    int NumPoints,NumAps;
    t_dArray1D zz,w,x;
    Integrator_ID id;

    C_IntegrationNodes(int _NumPoints, double _LimDown, double _LimUp, int _NumAps, Integrator_ID _id){
        NumPoints=_NumPoints;
        LimUp=_LimUp;
        LimDown=_LimDown;
        NumAps=_NumAps;
        id=_id;
        zz.resize(NumPoints+1);
        x.resize(NumPoints+1);
        w.resize(NumPoints+1);
        setNodes(id);
    }

    //Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
    //arrays x[1..n] and w[1..n] of length n, containing the abscissas and weights of the Gauss-
    //Legendre n-point quadrature formula.
    void gauleg(double x1, double x2, t_dArray1D & x, t_dArray1D & w, int n)
    {
        int m,j,i;
        double z1,z,xm,xl,pp,p3,p2,p1;
        m=(n+1)/2;

        xm=0.5*(x2+x1);
        xl=0.5*(x2-x1);
        for (i=1;i<=m;i++) {

            z=cos(3.141592654*(i-0.25)/(n+0.5));
            do {
                p1=1.0;
                p2=0.0;
                for (j=1;j<=n;j++) {
                    p3=p2;
                    p2=p1;

                    p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
                }
                pp=n*(z*p1-p2)/(z*z-1.0);
                z1=z;
                z=z1-p1/pp;


            } while (fabs(z-z1) > 3.0e-11);
            x[i]=xm-xl*z;

            x[n+1-i]=xm+xl*z;

            w[i]=2.0*xl/((1.0-z*z)*pp*pp);

            w[n+1-i]=w[i];
        }
    }

    //This routine returns arrays x[1..n] and w[1..n] of length n,
    //containing the abscissas and weights of the Gauss-Chebyshev n-point quadrature formula.
    void gaucheb(double x1, double x2, t_dArray1D & x, t_dArray1D & w, int n)
    {
        double zr;
        zr=0.5*(x2-x1);
        for (int i = 0; i <= n; i++){
            x[i]=0.0;
            w[i]=0.0;
        }

        for (int i = 1; i <= n ; i++){
            x[i]=zr*cos((i/(n+1.0))*pi);
            w[i]=(pi/(n+1.0))*sin(i/(n+1.0)*pi)*sin(i/(n+1.0)*pi);
        }
    }

    void setNodes(Integrator_ID _id){
        double aa,bb,zm,zr,dz;
        switch (_id){
            case qgausleg_log_ID:
                gauleg(0.0,1.0,zz,w,NumPoints);
                aa=log(LimDown);
                bb=log(LimUp);
                zm=(aa);
                zr=(bb-aa);
                for (int j=1;j<=NumPoints;j++)
                {
                    dz=zr*zz[j];
                    x[j]=exp(zm+dz);
                    w[j]=w[j]*exp(zm+dz)*zr;
                }
                break;

            case qgausleg_lin_ID:
                gauleg(0.0,1.0,zz,w,NumPoints);
                aa=(LimDown);
                bb=(LimUp);
                zm=(aa);
                zr=(bb-aa);
                for (int j=1;j<=NumPoints;j++)
                {
                    dz=zr*zz[j];
                    x[j]=(zm+dz - zr*zz[1]);
                    w[j]=w[j]*zr;
                }
                break;

            case qgauscheb_ID:
                gaucheb(-1.0,1.0,zz,w,NumPoints);
                for (int j=1;j<=NumPoints;j++)
                {
                    x[j]=zz[j];
                    //w[j]=w[j]; //weight stays the same
                }
                break;

            case qgausleg_sym_ID:
                gauleg(-1.0,1.0,zz,w,NumPoints);
                aa=(LimDown);
                bb=(LimUp);
                zm=0.5*(bb+aa);
                zr=0.5*(bb-aa);
                for (int j=1;j<=NumPoints;j++)
                {
                    dz=zr*zz[j];
                    x[j]=(zm+dz);
                    w[j]=w[j]*zr;
                }
                break;

            default:
                assert( false);
        }
    }

public:
    void getNodes(t_dArray1D & _x, t_dArray1D & _w) {
        _x = x;
        _w = w;
    }
};


#endif //_DSEPP_INTEGRATIONNODES_HPP_
