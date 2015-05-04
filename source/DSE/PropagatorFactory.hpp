//
// Created by stkubr.
//

#ifndef _DSEPP_PROPAGATORFACTORY_HPP_
#define _DSEPP_PROPAGATORFACTORY_HPP_

#include "Gluon.hpp"
#include "QuarkTypes.hpp"

namespace Propagators {

    class C_Gluon_Factory {
    private:
        C_Gluon_Factory() { }

    public:
        static C_Gluon_Factory &instance() {
            static C_Gluon_Factory ins;
            return ins;
        }

        C_Propagator *Create(int _id) {
            Gluon_ID id = (Gluon_ID) _id;
            return C_Gluon::createGluon(id);
        }

        C_Propagator *Create(int _id, std::string &_InterpolationPointsPath) {
            Gluon_ID id = (Gluon_ID) _id;
            return C_Gluon::createGluon(id, _InterpolationPointsPath);
        }
    };

    class C_Quark_Factory {
    private:
        C_Quark_Factory() { }

    public:
        static C_Quark_Factory &instance() {
            static C_Quark_Factory ins;
            return ins;
        }

        C_Propagator *Create(int _id) {
            Quark_ID id = (Quark_ID) _id;
            C_Quark *p;
            switch (id) {
                case Up_ID:
                    p = new C_Up_Quark();
                    break;
                case Down_ID:
                    p = new C_Down_Quark();
                    break;
                case Strange_ID:
                    p = new C_Strange_Quark();
                    break;
                case Charm_ID:
                    p = new C_Charm_Quark();
                    break;
                case Test_ID:
                    p = new C_Test_Quark();
                    break;
                default:
                    assert(false);
            }
            return p;
        }
    };

}
#endif //_DSEPP_PROPAGATORFACTORY_HPP_
