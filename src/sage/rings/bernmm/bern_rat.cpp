/*
   bern_rat.cpp:  multi-modular algorithm for computing Bernoulli numbers

   Copyright (C) 2008, 2009, David Harvey

   This file is part of the bernmm package (version 1.1).

   bernmm is released under a BSD-style license. See the README file in
   the source distribution for details.
*/

#include <gmp.h>
#include <NTL/ZZ.h>
#include <cmath>
#include <vector>
#include <set>
#include "bern_modp_util.h"
#include "bern_modp.h"
#include "bern_rat.h"

#ifdef USE_THREADS
#include <pthread.h>
#endif


using namespace std;
using namespace NTL;


namespace bernmm {


/*
   Computes the denominator of B_k using Clausen/von Staudt.
*/
void bern_den(mpz_t res, long k, const PrimeTable& table)
{
   mpz_set_ui(res, 1);

   // loop through factors of k
   for (long f = 1; f*f <= k; f++)
   {
      // if f divides k....
      if (k % f == 0)
      {
         // ... then both f + 1 and k/f + 1 are candidates for primes
         // dividing the denominator of B_k
         if (table.is_prime(f + 1))
            mpz_mul_ui(res, res, f + 1);

         if (f*f != k)
            if (table.is_prime(k/f + 1))
               mpz_mul_ui(res, res, k/f + 1);
      }
   }
}


// width of interval for each block
#define BLOCK_SIZE 1000


/*
   Represents that B_k is congruent to _residue_ modulo _modulus_.
*/
struct Item
{
   mpz_t modulus;
   mpz_t residue;

   Item()
   {
      mpz_init(modulus);
      mpz_init(residue);
   }

   ~Item()
   {
      mpz_clear(residue);
      mpz_clear(modulus);
   }
};


/*
   Items get sorted by modulus.
*/
struct Item_cmp
{
   bool operator()(const Item* x, const Item* y) const
   {
      return mpz_cmp(x->modulus, y->modulus) < 0;
   }
};


/*
   Returns new Item that combines information from op1 and op2 via CRT.
*/
Item* CRT(Item* op1, Item* op2)
{
   Item* res = new Item;

   // let n1, n2 be the moduli, and r1, r2 be the residues

   // res->modulus = t, where t = 0 mod n1, t = 1 mod n2
   mpz_invert(res->modulus, op1->modulus, op2->modulus);
   mpz_mul(res->modulus, res->modulus, op1->modulus);

   // res->residue = r2 - r1
   mpz_sub(res->residue, op2->residue, op1->residue);
   // res->residue = t * (r2 - r1)
   mpz_mul(res->residue, res->residue, res->modulus);
   // res->residue = r1 + t * (r2 - r1)
   mpz_add(res->residue, res->residue, op1->residue);
   // res->modulus = n1 * n2
   mpz_mul(res->modulus, op1->modulus, op2->modulus);
   // res->residue = r1 mod n1, r2 = mod n2
   mpz_mod(res->residue, res->residue, res->modulus);

   return res;
}


struct State
{
   long k;
   long bound;   // only use primes less than this bound
   const PrimeTable* table;

   // index of block that should be processed next
   long next;

   std::set<Item*, Item_cmp> items;
#ifdef USE_THREADS
   pthread_mutex_t lock;
#endif

   State(long k, long bound, const PrimeTable& table)
   {
      this->k = k;
      this->bound = bound;
      this->next = 0;
      this->table = &table;
#ifdef USE_THREADS
      pthread_mutex_init(&lock, NULL);
#endif
   }

   ~State()
   {
#ifdef USE_THREADS
      pthread_mutex_destroy(&lock);
#endif
   }
};


void* worker(void* arg)
{
   State& state = *((State*) arg);
   long k = state.k;

#ifdef USE_THREADS
   pthread_mutex_lock(&state.lock);
#endif

   while (1)
   {
      if (state.next * BLOCK_SIZE < state.bound)
      {
         // need to generate more modular data

         long next = state.next++;
#ifdef USE_THREADS
         pthread_mutex_unlock(&state.lock);
#endif

         Item* item = new Item;

         mpz_set_ui(item->modulus, 1);
         mpz_set_ui(item->residue, 0);

         for (long p = max(5, state.table->next_prime(next * BLOCK_SIZE));
              p < state.bound && p < (next+1) * BLOCK_SIZE;
              p = state.table->next_prime(p))
         {
            if (k % (p-1) == 0)
               continue;

            // compute B_k mod p
            long b = bern_modp(p, k);

            // CRT into running total
            long x = MulMod(SubMod(b, mpz_fdiv_ui(item->residue, p), p),
                            InvMod(mpz_fdiv_ui(item->modulus, p), p), p);
            mpz_addmul_ui(item->residue, item->modulus, x);
            mpz_mul_ui(item->modulus, item->modulus, p);
         }

#ifdef USE_THREADS
         pthread_mutex_lock(&state.lock);
#endif
         state.items.insert(item);
      }
      else
      {
         // all modular data has been generated

         if (state.items.size() <= 1)
         {
            // no more CRTs for this thread to perform
#ifdef USE_THREADS
            pthread_mutex_unlock(&state.lock);
#endif
            return NULL;
         }

         // CRT two smallest items together
         Item* item1 = *(state.items.begin());
         state.items.erase(state.items.begin());
         Item* item2 = *(state.items.begin());
         state.items.erase(state.items.begin());
#ifdef USE_THREADS
         pthread_mutex_unlock(&state.lock);
#endif

         Item* item3 = CRT(item1, item2);
         delete item1;
         delete item2;

#ifdef USE_THREADS
         pthread_mutex_lock(&state.lock);
#endif
         state.items.insert(item3);
      }
   }
}


void bern_rat(mpq_t res, long k, int num_threads)
{
   // special cases

   if (k == 0)
   {
      // B_0 = 1
      mpq_set_ui(res, 1, 1);
      return;
   }

   if (k == 1)
   {
      // B_1 = -1/2
      mpq_set_si(res, -1, 2);
      return;
   }

   if (k == 2)
   {
      // B_2 = 1/6
      mpq_set_si(res, 1, 6);
      return;
   }

   if (k & 1)
   {
      // B_k = 0 if k is odd
      mpq_set_ui(res, 0, 1);
      return;
   }

   if (num_threads <= 0)
      num_threads = 1;

   mpz_t num, den;
   mpz_init(num);
   mpz_init(den);

   const double log2 =    0.69314718055994528622676;
   const double invlog2 = 1.44269504088896340735992;   // = 1/log(2)

   // compute preliminary prime bound and build prime table
   long bound1 = (long) max(37.0, ceil((k + 0.5) * log(k) * invlog2));
   PrimeTable table(bound1);

   // compute denominator of B_k
   bern_den(den, k, table);

   // compute number of bits we need to resolve the numerator
   long bits = (long) ceil((k + 0.5) * log(k) * invlog2 - 4.094 * k + 2.470
                                      + log(mpz_get_d(den)) * invlog2);

   // compute tighter prime bound
   // (note: we can safely get away with double-precision here. It would
   // only start being insufficient around k = 10^13 or so, which is totally
   // impractical at present.)
   double prod = 1.0;
   long prod_bits = 0;
   long p;
   for (p = 5; prod_bits < bits + 1; p = table.next_prime(p))
   {
      if (p >= NTL_SP_BOUND)
         abort();   // !!!!! not sure what else we can do here...
      if (k % (p-1) != 0)
         prod *= (double) p;
      int exp;
      prod = frexp(prod, &exp);
      prod_bits += exp;
   }
   long bound2 = p;

   State state(k, bound2, table);

#ifdef USE_THREADS
   vector<pthread_t> threads(num_threads - 1);

   pthread_attr_t attr;
   pthread_attr_init(&attr);
#ifdef THREAD_STACK_SIZE
   pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE * 1024);
#endif

   // spawn worker threads to process blocks
   for (long i = 0; i < num_threads - 1; i++)
      pthread_create(&threads[i], &attr, worker, &state);
#endif

   worker(&state);    // make this thread a worker too

#ifdef USE_THREADS
   for (long i = 0; i < num_threads - 1; i++)
      pthread_join(threads[i], NULL);

   pthread_attr_destroy (&attr);
#endif

   // reconstruct B_k as a rational number
   Item* item = *(state.items.begin());
   mpz_mul(num, item->residue, den);
   mpz_mod(num, num, item->modulus);

   if (k % 4 == 0)
   {
      // B_k is negative
      mpz_sub(num, item->modulus, num);
      mpz_neg(num, num);
   }

   delete item;

   mpz_swap(num, mpq_numref(res));
   mpz_swap(den, mpq_denref(res));

   mpz_clear(num);
   mpz_clear(den);
}



};    // end namespace


// end of file ================================================================
