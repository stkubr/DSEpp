#ifndef DSEPP_INTERPOLATION_H
#define DSEPP_INTERPOLATION_H

#include <vector>
#include "../types.h"

namespace Interpolation
{

    template<typename x_type, typename y_type> class C_Linear {
    public:
        C_Linear(int n, std::vector<x_type> & x, std::vector<y_type> & y) {
            m_x.clear();
            m_y.clear();
            m_x.resize(n);
            m_y.resize(n);

            for (int i = 0; i < n; ++i) {
                m_x[i] = x[i];
                m_y[i] = y[i];
            }
        }

        double inline Real(t_cmplx x) {
            return real(x);
        }
        double inline Real(double x) {
            return (x);
        }
        y_type getValue(x_type x) {
            int i = 0;
            while (Real(x) > Real(m_x[++i]));
            y_type a = (x - m_x[i - 1]) / (m_x[i] - m_x[i - 1]);
            return m_y[i - 1] + a * (m_y[i] - m_y[i - 1]);
        }

    private:
        std::vector<x_type> m_x;
        std::vector<y_type> m_y;
    };
}

#endif //DSEPP_INTERPOLATION_H
