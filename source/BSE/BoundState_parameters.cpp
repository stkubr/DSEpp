//
// Created by stkubr on 21.04.15.
//

#include "BoundState_parameters.h"

using namespace BSE;

void C_BoundState_parameters::print() {
    std::cout <<
    "LimDk" << " - " << LimDk << std::endl <<
    "LimUk" << " - " << LimUk << std::endl <<
    "zetta_part" << " - " << zetta_part << std::endl <<
    "OffShell" << " - " << OffShell << std::endl <<

    "NumRadial" << " - " << NumRadial << std::endl <<
    "Cheb_order" << " - " << Cheb_order << std::endl <<
    "NumCheb_nod1" << " - " << NumCheb_nod1 << std::endl <<
    "NumCheb_nod2" << " - " << NumCheb_nod2 << std::endl <<
    "NumAngleY" << " - " << NumAngleY << std::endl <<

    "NumRadial_Contour" << " - " << NumRadial_Contour << std::endl <<
    "NumCheb_Contour" << " - " << NumCheb_Contour << std::endl <<
    "NumAngleY_Contour" << " - " << NumAngleY_Contour << std::endl;
}

void C_BoundState_parameters::setParams(std::string &_ParamPath) {
    std::string line;
    std::ifstream _ParamList(_ParamPath);
    if ((_ParamList).is_open()) {
        while ( (_ParamList).good() ) {
            (_ParamList) >> line >> LimDk;
            (_ParamList) >> line >> LimUk;
            (_ParamList) >> line >> zetta_part;
            (_ParamList) >> line >> OffShell;

            (_ParamList) >> line >> NumRadial;
            (_ParamList) >> line >> Cheb_order;
            (_ParamList) >> line >> NumCheb_nod1;
            (_ParamList) >> line >> NumCheb_nod2;
            (_ParamList) >> line >> NumAngleY;

            (_ParamList) >> line >> NumRadial_Contour;
            (_ParamList) >> line >> NumCheb_Contour;
            (_ParamList) >> line >> NumAngleY_Contour;
        }
    }
    else {
        std::cout << "Cant open BoundState parameters file!" << std::endl;
        assert(false);
    }
}

void C_BoundState_parameters::setDefault() {
    std::cout << "Default BoundState parameters!" << std::endl;
    LimDk = 0.0001;
    LimUk = 2000;
    zetta_part = 0.5;
    OffShell = 0;

    NumRadial = 20;
    Cheb_order = 2;
    NumCheb_nod1 = 4;
    NumCheb_nod2 = 4;
    NumAngleY = 4;

    NumRadial_Contour = 48;
    NumCheb_Contour = 8;
    NumAngleY_Contour = 8;
}