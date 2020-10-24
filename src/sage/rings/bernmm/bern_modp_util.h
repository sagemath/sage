/*
   bern_modp_util.h:  number-theoretic utility functions

   Copyright (C) 2008, 2009, David Harvey

   This file is part of the bernmm package (version 1.1).

   bernmm is released under a BSD-style license. See the README file in
   the source distribution for details.
*/


#ifndef BERNMM_BERN_MODP_UTIL_H
#define BERNMM_BERN_MODP_UTIL_H


#include <vector>
#include <cassert>
#include <climits>

#include <NTL/ZZ.h>

#if ULONG_MAX == 4294967295U
#define ULONG_BITS 32
#elif ULONG_MAX == 18446744073709551615U
#define ULONG_BITS 64
#else
#error Oops! Unsigned long is neither 32 nor 64 bits.
#error You need to update bern_modp_util.h.
#endif


namespace bernmm {


/*
   Same as NTL's PowerMod, but also accepts an _ninv_ parameter, which is the
   same as the ninv parameter for NTL's MulMod routines, i.e. should have
   ninv = PrepMulMod(n).

   (Implementation is adapted from ZZ.c in NTL 5.4.1.)
*/
long PowerMod(long a, long ee, long n, NTL::mulmod_t ninv);


/*
   Represents the factorisation of an integer n into distinct prime factors.

   (Very naive implementation!)
*/
class Factorisation
{
protected:
   /*
      Finds distinct prime factors of m in the range k < p <= m.
      Assumes that m does not have any prime factors p <= k.
      Appends factors found to _factors_.
   */
   void helper(long k, long m);

public:
   // the integer
   long n;

   // the distinct factors (in increasing order)
   std::vector<long> factors;

   // initialises with given integer
   Factorisation(long n);
};



class PrimeTable
{
private:
   std::vector<long> data;   // bit-vector; 0 means prime, 1 means composite

   // read bit from index i
   inline bool get(long i) const
   {
      return (data[i / ULONG_BITS] >> (i % ULONG_BITS)) & 1;
   }

   // set bit at index i
   inline void set(long i)
   {
      data[i / ULONG_BITS] |= (1L << (i % ULONG_BITS));
   }


public:
   // initialise with primes up to given bound
   PrimeTable(long bound);

   // test whether n is prime by table lookup
   inline bool is_prime(long n) const
   {
      return !get(n);
   }

   // returns smallest prime p that is larger than n
   long next_prime(long n) const
   {
      for (n++; !is_prime(n); n++);
      return n;
   }
};



/*
   Returns 1 if n is prime.
*/
int is_prime(long n);


/*
   Returns smallest prime larger than p.
*/
long next_prime(long p);


/*
   Computes order of x mod p, given the factorisation F of p-1.
*/
long order(long x, long p, NTL::mulmod_t pinv, const Factorisation& F);


/*
   Finds the smallest primitive root mod p, given the factorisation F of p-1.
*/
long primitive_root(long p, NTL::mulmod_t pinv, const Factorisation& F);


};    // end namespace


#endif

// end of file ================================================================
