#ifndef DSEPP_BARYCENTRICINTERPOLATION_H
#define DSEPP_BARYCENTRICINTERPOLATION_H

#include "../types.h"

namespace Interpolation{

    class C_Base_interp{
    public:
        int n, mm, jsav, cor, dj;
        t_dArray1D xx, yy;
        C_Base_interp(t_dArray1D &x, t_dArray1D &y, int m)
                :n(x.size()), mm(m), jsav(0), cor(0) {

            xx=x;
            yy=y;
            dj = min(1,(int)pow((double)n,0.25));
        }

        double interp(double x) {
            int jlo = cor ? hunt(x) : locate(x);
            return rawinterp(jlo,x);
        }

        int locate(const double x);
        int hunt(const double x);

        double virtual rawinterp(int jlo, double x) = 0;
    };

    int C_Base_interp::locate(const double x){
        int ju,jm,jl;
        if (n < 2 || mm < 2 || mm > n) throw("locate size error");
        bool ascnd=(xx[n-1] >= xx[0]);
        jl=0;
        ju=n-1;
        while (ju-jl > 1) {
            jm = (ju+jl) >> 1;
            if (x >= xx[jm] == ascnd)
                jl=jm;
            else
                ju=jm;
        }
        cor = abs(jl-jsav) > dj ? 0 : 1;
        jsav = jl;
        return max(0,min(n-mm,jl-((mm-2)>>1)));
    }

    int C_Base_interp::hunt(const double x){
        int jl=jsav, jm, ju, inc=1;
        if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
        bool ascnd=(xx[n-1] >= xx[0]);
        if (jl < 0 || jl > n-1) {
            jl=0;
            ju=n-1;
        } else {
            if (x >= xx[jl] == ascnd) {
                for (;;) {
                    ju = jl + inc;
                    if (ju >= n-1) { ju = n-1; break;}
                    else if (x < xx[ju] == ascnd) break;
                    else {
                        jl = ju;
                        inc += inc;
                    }
                }
            } else {
                ju = jl;
                for (;;) {
                    jl = jl - inc;
                    if (jl <= 0) { jl = 0; break;}
                    else if (x >= xx[jl] == ascnd) break;
                    else {
                        ju = jl;
                        inc += inc;
                    }
                }
            }
        }

        while (ju-jl > 1) {
            jm = (ju+jl) >> 1;
            if (x >= xx[jm] == ascnd)
                jl=jm;
            else
                ju=jm;
        }
        cor = abs(jl-jsav) > dj ? 0 : 1;
        jsav = jl;
        return max(0,min(n-mm,jl-((mm-2)>>1)));
    }

    class C_BaryRat_interp : public C_Base_interp{
    public:
        t_dArray1D w;
        int d;
        C_BaryRat_interp(t_dArray1D &xv, t_dArray1D &yv, int dd);
        double rawinterp(int jl, double x);
        double interp(double x);
    };

    C_BaryRat_interp::C_BaryRat_interp(t_dArray1D &xv, t_dArray1D &yv, int dd)
            : C_Base_interp(xv,yv,xv.size()), w(n), d(dd){
        if (n<=d) throw("d too big");
        for (int k=0;k<n;k++) {
            int imin=max(k-d,0);
            int imax = k >= n-d ? n-d-1 : k;
            double temp = imin & 1 ? -1.0 : 1.0;
            double sum=0.0;
            for (int i=imin;i<=imax;i++) {
                int jmax=min(i+d,n-1);
                double term=1.0;
                for (int j=i;j<=jmax;j++) {
                    if (j==k) continue;
                    term *= (xx[k]-xx[j]);
                }
                term=temp/term;
                temp=-temp;
                sum += term;
            }
            w[k]=sum;
        }
    }

    double C_BaryRat_interp::rawinterp(int jl, double x){
        double num=0,den=0;
        for (int i=0;i<n;i++) {
            double h=x-xx[i];
            if (h == 0.0) {
                return yy[i];
            } else {
                double temp=w[i]/h;
                num += temp*yy[i];
                den += temp;
            }
        }
        return num/den;
    }

    double C_BaryRat_interp::interp(double x) {
        return rawinterp(1,x);
    }

}

#endif //DSEPP_BARYCENTRICINTERPOLATION_H
