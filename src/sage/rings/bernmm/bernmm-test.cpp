/*
   bernmm-test.cpp:  test module

   Copyright (C) 2008, 2009, David Harvey

   This file is part of the bernmm package (version 1.1).

   bernmm is released under a BSD-style license. See the README file in
   the source distribution for details.
*/

#include <iostream>
#include <NTL/ZZ.h>
#include <gmp.h>
#include "bern_modp_util.h"
#include "bern_modp.h"
#include "bern_rat.h"


NTL_CLIENT;


using namespace bernmm;
using namespace std;


/*
   Computes B_0, B_1, ..., B_{n-1} using naive algorithm, writes them to res.
*/
void bern_naive(mpq_t* res, long n)
{
   mpq_t t, u;
   mpq_init(t);
   mpq_init(u);

   // compute res[j] = B_j / j! for 0 <= j < n
   if (n > 0)
      mpq_set_si(res[0], 1, 1);

   for (long j = 1; j < n; j++)
   {
      mpq_set_si(res[j], 0, 1);
      mpq_set_ui(t, 1, 1);
      for (long k = 0; k < j; k++)
      {
         mpz_mul_ui(mpq_denref(t), mpq_denref(t), k + 2);
         mpq_mul(u, res[j - 1 - k], t);
         mpq_sub(res[j], res[j], u);
      }
   }

   // multiply through by j! for 0 <= j < n
   mpq_set_ui(t, 1, 1);
   for (long j = 2; j < n; j++)
   {
      mpz_mul_ui(mpq_numref(t), mpq_numref(t), j);
      mpq_mul(res[j], res[j], t);
   }

   mpq_clear(u);
   mpq_clear(t);
}


/*
   Tests _bern_modp_powg() for a given p and k by comparing against the
   rational number B_k (must be supplied in b).

   Returns 1 on success.
*/
int testcase__bern_modp_powg(long p, long k, mpq_t b)
{
   mulmod_t pinv = PrepMulMod(p);

   // compute B_k mod p using _bern_modp_powg()
   long x = _bern_modp_powg(p, pinv, k);
   x = MulMod(x, k, p, pinv);

   // compute B_k mod p from rational B_k
   long y = mpz_fdiv_ui(mpq_numref(b), p);
   long z = mpz_fdiv_ui(mpq_denref(b), p);
   return y == MulMod(z, x, p, pinv);
}



/*
   Tests _bern_modp_powg() by comparing against naive computation of B_k
   (as a rational) for a range of small p and k.

   Returns 1 on success.
*/
int test__bern_modp_powg()
{
   int success = 1;

   const long MAX = 300;
   mpq_t bern[MAX];

   // compute B_k's as rational numbers using naive algorithm
   for (long i = 0; i < MAX; i++)
      mpq_init(bern[i]);
   bern_naive(bern, MAX);

   // try a range of k's
   for (long k = 2; k < MAX && success; k += 2)
   {
      // try a range of small p's
      for (long p = k + 3; p < 2*MAX && success; p += 2)
      {
         if (!ProbPrime(p))
            continue;
         success = success && testcase__bern_modp_powg(p, k, bern[k]);
      }

      // try a single larger p
      success = success && testcase__bern_modp_powg(1000003, k, bern[k]);
   }

   // if we're on a 32-bit machine, try a single example with p right near
   // NTL's boundary (this is infeasible on a 64-bit machine)
   if (NTL_SP_NBITS <= 32)
   {
      long p = NTL_SP_BOUND - 1;
      while (!ProbPrime(p))
         p--;

      long k = (MAX/2)*2 - 2;
      success = success && testcase__bern_modp_powg(p, k, bern[k]);
   }

   for (long i = 0; i < MAX; i++)
      mpq_clear(bern[i]);

   return success;
}



/*
   Tests _bern_modp_pow2() for a given p and k by comparing against result
   from _bern_modp_powg().

   Returns 1 on success.

   If 2^k = 1 mod p, then _bern_modp_pow2() won't work, so it just returns 1.
*/
int testcase__bern_modp_pow2(long p, long k)
{
   mulmod_t pinv = PrepMulMod(p);

   if (PowerMod(2, k, p, pinv) == 1)
      return 1;

   long x = _bern_modp_powg(p, pinv, k);
   long y = _bern_modp_pow2(p, pinv, k);

   return x == y;
}



/*
   Tests _bern_modp_pow2() by comparing against _bern_modp_powg() for
   a range of p and k.

   Returns 1 on success.
*/
int test__bern_modp_pow2()
{
   int success = 1;

   // exhaustive comparison over some small p and k
   for (long p = 5; p < 2000 && success; p += 2)
   {
      if (!ProbPrime(p))
         continue;

      for (long k = 2; k <= p - 3 && success; k += 2)
         success = success && testcase__bern_modp_pow2(p, k);
   }

   // a few larger values of p
   for (long p = 1000000; p < 1030000; p++)
   {
      if (!ProbPrime(p))
         continue;

      long k = 2 * (rand() % ((p-3)/2)) + 2;
      success = success && testcase__bern_modp_pow2(p, k);
   }

   // if we're on a 32-bit machine, try a single example with p right near
   // NTL's boundary (this is infeasible on a 64-bit machine)
   if (NTL_SP_NBITS <= 32)
   {
      long p = NTL_SP_BOUND - 1;
      while (!ProbPrime(p))
         p--;
      success = success & testcase__bern_modp_pow2(p, 10);
   }

   // try a few just below the REDC barrier
   if (ULONG_BITS == 32)
   {
      long boundary = 1L << (ULONG_BITS/2 - 1);
      for (long p = boundary - 1000; p < boundary && success; p++)
      {
         if (ProbPrime(p))
         {
            for (long trial = 0; trial < 1000 && success; trial++)
            {
               long k = 2 * (rand() % ((p-3)/2)) + 2;
               success = success && testcase__bern_modp_pow2(p, k);
            }
         }
      }
   }
   else
   {
      // on a 64-bit machine, only try one, since these are huge!
      long p = 1L << (ULONG_BITS/2 - 1);
      while (!ProbPrime(p))
         p--;
      success = success && testcase__bern_modp_pow2(p, 10);
   }

   return success;
}


/*
   Tests bern_rat() by comparing against the naive algorithm for several small
   k, and testing against bern_modp() for a couple of larger k.

   Returns 1 on success.
*/
int test_bern_rat()
{
   int success = 1;

   const long MAX = 300;
   mpq_t bern[MAX];

   // compute B_k's as rational numbers using naive algorithm
   for (long i = 0; i < MAX; i++)
      mpq_init(bern[i]);
   bern_naive(bern, MAX);

   mpq_t x;
   mpq_init(x);

   // exhaustive test for small k
   for (long k = 0; k < MAX && success; k++)
   {
      bern_rat(x, k, 4);    // try with 4 threads just for fun
      success = success && mpq_equal(x, bern[k]);
   }

   // try a few larger k
   for (long i = 0; i < 50 && success; i++)
   {
      long k = ((random() % 20000) / 2) * 2;
      bern_rat(x, k, 4);

      // compare with modular information
      long p = 1000003;
      long num = mpz_fdiv_ui(mpq_numref(x), p);
      long den = mpz_fdiv_ui(mpq_denref(x), p);
      success = success && (MulMod(bern_modp(p, k), den, p) == num);
   }

   mpq_clear(x);
   for (long i = 0; i < MAX; i++)
      mpq_clear(bern[i]);

   return success;
}


void report(int success)
{
   if (success)
      cout << "ok" << endl;
   else
   {
      cout << "failed!" << endl;
      abort();
   }
}


int main(int argc, char* argv[])
{
   if (argc == 1)
   {
      cout << "bernmm test module" << endl;
      cout << endl;
      cout << "   bernmm-test --test" << endl;
      cout << "        runs test suite" << endl;
      cout << "   bernmm-test --rational <k> <threads>" << endl;
      cout << "        computes B_k with <threads> threads" << endl;
      cout << "   bernmm-test --modular <p> <k>" << endl;
      cout << "        computes B_k mod p" << endl;
      return 0;
   }

   if (!strcmp(argv[1], "--test"))
   {
      cout << "testing _bern_modp_powg()... " << flush;
      report(test__bern_modp_powg());

      cout << "testing _bern_modp_pow2()... " << flush;
      report(test__bern_modp_pow2());

      cout << "testing bern_rat()... " << flush;
      report(test_bern_rat());
   }
   else if (!strcmp(argv[1], "--rational"))
   {
      if (argc <= 3)
      {
         cout << "not enough arguments" << endl;
         return 0;
      }
      long k = atol(argv[2]);
      long threads = atol(argv[3]);
      mpq_t r;
      mpq_init(r);
      bern_rat(r, k, threads);
      gmp_printf("%Zd/%Zd\n", mpq_numref(r), mpq_denref(r));
      mpq_clear(r);
   }
   else if (!strcmp(argv[1], "--modular"))
   {
      if (argc <= 3)
      {
         cout << "not enough arguments" << endl;
         return 0;
      }
      long p = atol(argv[2]);
      long k = atol(argv[3]);
      cout << bern_modp(p, k) << endl;
   }
   else
   {
      cout << "unknown command" << endl;
   }

   return 0;
}


// end of file ================================================================
