//
// Created by stkubr on 08.04.15.
//

#ifndef _DSEPP_MEMORYFACTORIES_HPP_
#define _DSEPP_MEMORYFACTORIES_HPP_

#include "DedicMem.hpp"

class C_MemoryFactory_Quark{
public:
    C_DedicMem_Abs* CreateMemory() {
        return new C_DedicMem_Quark;
    }
};

class C_MemoryFactory_Kernel{
public:
    C_DedicMem_Abs* CreateMemory() {
        return new C_DedicMem_Kernel;
    }
};

class C_MemoryFactory_BSE{
public:
    C_DedicMem_Abs* CreateMemory() {
        return new C_DedicMem_BSE;
    }
};

extern C_MemoryFactory_Quark * DedicMemFactory_Quark;
extern C_MemoryFactory_Kernel * DedicMemFactory_Kernel;
extern C_MemoryFactory_BSE * DedicMemFactory_BSE;

#endif //_DSEPP_MEMORYFACTORIES_HPP_
