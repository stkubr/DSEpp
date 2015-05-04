//
// Created by stkubr on 21.04.15.
//

#ifndef DSEPP_BOUNDSTATE_PARAMETERS_H
#define DSEPP_BOUNDSTATE_PARAMETERS_H

#include <iostream>
#include <string>
#include <fstream>
#include <assert.h>

namespace BSE {

    class C_BoundState_parameters {
    public:
        int NumRadial, NumCheb_nod1, NumCheb_nod2, Cheb_order, NumAngleY;
        int NumRadial_Contour, NumCheb_Contour, NumAngleY_Contour;
        double LimUk, LimDk, zetta_part;
        bool OffShell;

        void print();

        void setParams(std::string &_ParamPath);

        void setDefault();
    };

}
#endif //DSEPP_BOUNDSTATE_PARAMETERS_H
