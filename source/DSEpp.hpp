//
// Created by stkubr on 01.04.15.
//

#ifndef _DSEPP_DSEPP_HPP_
#define _DSEPP_DSEPP_HPP_

#ifndef _NUM_THREADS
#define _NUM_THREADS omp_get_max_threads()
#endif

class C_AbstractKernel;
class C_Propagator;

#include "../source/types.h"
#include "../source/DSE/Gluon.hpp"
#include "../source/Kernel/KernelFactory.hpp"
#include "../source/DSE/Quark_id.hpp"


#endif //_DSEPP_DSEPP_HPP_
