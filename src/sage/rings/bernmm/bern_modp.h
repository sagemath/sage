/*
   bern_modp.cpp:  computing isolated Bernoulli numbers modulo p

   Copyright (C) 2008, 2009, David Harvey

   This file is part of the bernmm package (version 1.1).

   bernmm is released under a BSD-style license. See the README file in
   the source distribution for details.
*/

#ifndef BERNMM_BERN_MODP_H
#define BERNMM_BERN_MODP_H

#include <NTL/ZZ.h>

namespace bernmm {


/*
   Returns B_k mod p, in [0, p), or -1 if B_k is not p-integral.

   PRECONDITIONS:
      2 <= p < NTL_SP_BOUND, p prime
      k >= 0
*/
long bern_modp(long p, long k);


/*
   Exported for testing.
*/
long _bern_modp_powg(long p, NTL::mulmod_t pinv, long k);
long _bern_modp_pow2(long p, NTL::mulmod_t pinv, long k);


};


#endif

// end of file ================================================================
