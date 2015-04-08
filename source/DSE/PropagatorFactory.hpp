//
// Created by stkubr.
//

#ifndef _DSEPP_PROPAGATORFACTORY_HPP_
#define _DSEPP_PROPAGATORFACTORY_HPP_

#include "Gluon.hpp"
#include "QuarkTypes.hpp"

class C_Propagator_Factory{
public:
    virtual C_Propagator * Create(int)=0;
    virtual ~C_Propagator_Factory() {}
};

class C_Gluon_Factory: public C_Propagator_Factory{
public:
    C_Propagator* Create(int _id) {
        Gluon_ID id=(Gluon_ID)_id;
        return C_Gluon::createGluon( id );
    }

    C_Propagator* Create(int _id, std::string & _InterpolationPointsPath) {
        Gluon_ID id=(Gluon_ID)_id;
        return C_Gluon::createGluon( id, _InterpolationPointsPath );
    }
};

class C_Quark_Factory: public C_Propagator_Factory{
public:

    C_Propagator* Create(int _id) {
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

#endif //_DSEPP_PROPAGATORFACTORY_HPP_
