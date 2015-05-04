#ifndef DSEPP_CHEBYCHEVPOLINOMS_H
#define DSEPP_CHEBYCHEVPOLINOMS_H

#include <source/types.h>

/// returns Chebychev polinom of the order at x
inline t_cmplx Cheb_polinoms(t_cmplx x, int order) {
    if (order == 0) {
        return 1.0;
    }

    if (order == 1) {
        return 2.0 * x;
    }

    if (order == 2) {
        return 4.0 * x * x - 1.0;
    }

    if (order > 2) {
        t_cmplx U[order + 1];
        U[0] = 1.0;
        U[1] = 2.0 * x;
        U[2] = 4.0 * x * x - 1.0;
        for (int i = 3; i <= order; i++) {
            U[i] = 2.0 * x * U[i - 1] - U[i - 2];
        }
        return U[order];
    }
    else {
        std::cout << "Chebys Error!" << std::endl;
        assert(false);
    }
    return 0;
}

#endif //DSEPP_CHEBYCHEVPOLINOMS_H
