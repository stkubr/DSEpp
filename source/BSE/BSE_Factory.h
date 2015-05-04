//
// Created by stkubr on 01.05.15.
//

#ifndef DSEPP_BSE_FACTORY_H
#define DSEPP_BSE_FACTORY_H

#include "BSE_PseudoScalar.h"
#include "BSE_Scalar.h"
#include "BSE_Vector.h"

namespace BSE {

    enum BSE_ID {
        PseudoScalar_ID = 0, Scalar_ID, Vector_ID, BSE_ID_End
    };

    class C_BSE_Factory {
    private:
        C_BSE_Factory() { }

    public:
        static C_BSE_Factory &instance() {
            static C_BSE_Factory ins;
            return ins;
        }

        C_BSE *Create(BSE_ID id) {
            C_BSE *ptr;
            switch (id) {
                case PseudoScalar_ID:
                    ptr = new C_BSE_PseudoScalar();
                    break;
                case Scalar_ID:
                    ptr = new C_BSE_Scalar();
                    break;
                case Vector_ID:
                    ptr = new C_BSE_Vector();
                    break;
                default:
                    assert(false);
            }
            return ptr;
        }
    };

}

#endif //DSEPP_BSE_FACTORY_H
