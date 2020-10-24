/*
   bern_modp.cpp:  computing isolated Bernoulli numbers modulo p

   Copyright (C) 2008, 2009, David Harvey

   This file is part of the bernmm package (version 1.1).

   bernmm is released under a BSD-style license. See the README file in
   the source distribution for details.
*/


#include <limits.h>
#include <signal.h>
#include <cstring>
#include <gmp.h>
#include <NTL/ZZ.h>
#include "bern_modp_util.h"
#include "bern_modp.h"


NTL_CLIENT;


using namespace std;


namespace bernmm {


/******************************************************************************

   Computing the main sum (general case)

******************************************************************************/

/*
   Returns (1 - g^k) B_k / 2k mod p.

   PRECONDITIONS:
      5 <= p < NTL_SP_BOUND, p prime
      2 <= k <= p-3, k even
      pinv = PrepMulMod(p)
      g = a multiplicative generator of GF(p), in [0, p)
*/
long bernsum_powg(long p, mulmod_t pinv, long k, long g)
{
   long half_gm1 = (g + ((g & 1) ? 0 : p) - 1) / 2;    // (g-1)/2 mod p
   long g_to_jm1 = 1;
   long g_to_km1 = PowerMod(g, k-1, p, pinv);
   long g_to_km1_to_j = g_to_km1;
   long sum = 0;
   muldivrem_t g_pinv = PrepMulDivRem(g, p);
   mulmod_precon_t g_to_km1_pinv = PrepMulModPrecon(g_to_km1, p, pinv);

   for (long j = 1; j <= (p-1)/2; j++)
   {
      // at this point,
      //    g_to_jm1 holds g^(j-1) mod p
      //    g_to_km1_to_j holds (g^(k-1))^j mod p

      // update g_to_jm1 and compute q = (g*(g^(j-1) mod p) - (g^j mod p)) / p
      long q;
      g_to_jm1 = MulDivRem(q, g_to_jm1, g, p, g_pinv);

      // compute h = -h_g(g^j) = q - (g-1)/2
      long h = SubMod(q, half_gm1, p);

      // add h_g(g^j) * (g^(k-1))^j to running total
      sum = SubMod(sum, MulMod(h, g_to_km1_to_j, p, pinv), p);

      // update g_to_km1_to_j
      g_to_km1_to_j = MulModPrecon(g_to_km1_to_j, g_to_km1, p, g_to_km1_pinv);
   }

   return sum;
}



/******************************************************************************

   Computing the main sum (c = 1/2 case)

******************************************************************************/


/*
   The Expander class stores precomputed information for a fixed integer p,
   that subsequently permits fast computation of the binary expansion of s/p
   for any 0 < s < p.

   The constructor takes p and max_words as input. Must have 1 <= max_words <=
   MAX_INV. It computes an approximation to 1/p.

   The function expand(word_t* res, long s, long n) computes n words of s/p.
   Must have 0 < s < p and 1 <= n <= max_words. The output is written to res.
   The first word of output is junk. The next n words are the digits of s/p,
   from least to most significant. The buffer must be at least n+2 words long
   (even though the first and last words are never used for output).

   A "word" is a word_t, and contains WORD_BITS bits. On most systems this
   will be an mp_limb_t.
*/

#define MAX_INV 256

#if (GMP_NAIL_BITS == 0) && (GMP_LIMB_BITS >= ULONG_BITS)
// fast mpn-based version

typedef mp_limb_t word_t;
#define WORD_BITS GMP_LIMB_BITS

class Expander
{
private:
   // Approximation to 1/p. We store (max_words + 1) limbs.
   mp_limb_t pinv[MAX_INV + 2];
   mp_limb_t p;
   int max_words;

public:
   Expander(long p, int max_words)
   {
      assert(max_words >= 1);
      assert(max_words <= MAX_INV);

      this->max_words = max_words;
      this->p = p;
      mp_limb_t one = 1;
      mpn_divrem_1(pinv, max_words + 1, &one, 1, p);
   }

   void expand(word_t* res, long s, int n)
   {
      assert(s > 0 && s < p);
      assert(n >= 1);
      assert(n <= max_words);

      if (s == 1)
      {
         // already have 1/p; just copy it
         for (int i = 1; i <= n; i++)
            res[i] = pinv[max_words - n + i];
      }
      else
      {
         mpn_mul_1(res, pinv + max_words - n, n + 1, (mp_limb_t) s);

         // If the first output limb is really close to 0xFFFF..., then there's
         // a possibility of overflow, so fall back on doing division directly.
         // This should happen extremely rarely --- essentially never on a
         // 64-bit system, and very occasionally on a 32-bit system.
         if (res[0] > -((mp_limb_t) s))
         {
            mp_limb_t ss = s;
            mpn_divrem_1(res, n + 1, &ss, 1, p);
         }
      }
   }
};


#else
// slow mpz-based version, since GMP is using nails, or mp_limb_t is
// absurdly narrow

typedef unsigned long word_t;
#define WORD_BITS ULONG_BITS

class Expander
{
private:
   mp_limb_t p;
   mpz_t temp;

public:
   Expander(long p, int max_words)
   {
      this->p = p;
      mpz_init(temp);
   }

   ~Expander()
   {
      mpz_clear(temp);
   }

   void expand(word_t* res, long s, int n)
   {
      assert(s > 0 && s < p);
      assert(n >= 1);

      mpz_set_ui(temp, s);
      mpz_mul_2exp(temp, temp, WORD_BITS * n);
      mpz_fdiv_q_ui(temp, temp, p);
      mpz_export(res + 1, NULL, -1, sizeof(word_t), 0, 0, temp);
   }
};

#endif



/*
   Returns (2^(-k) - 1) 2 B_k / k  mod p.

   (Note: this is useless if 2^k = 1 mod p.)

   PRECONDITIONS:
      5 <= p < NTL_SP_BOUND, p prime
      2 <= k <= p-3, k even
      pinv = PrepMulMod(p)
      g = a multiplicative generator of GF(p), in [0, p)
      n = multiplicative order of 2 in GF(p)
*/

#define TABLE_LG_SIZE 8
#define TABLE_SIZE (((word_t) 1) << TABLE_LG_SIZE)
#define TABLE_MASK (TABLE_SIZE - 1)
#define NUM_TABLES (WORD_BITS / TABLE_LG_SIZE)

#if WORD_BITS % TABLE_LG_SIZE != 0
#error Number of bits in a long must be divisible by TABLE_LG_SIZE
#endif

long bernsum_pow2(long p, mulmod_t pinv, long k, long g, long n)
{
   // In the main summation loop we accumulate data into the _tables_ array;
   // tables[y][z] contributes to the final answer with a weight of
   //
   // sum(-(-1)^z[t] * (2^(k-1))^(WORD_BITS - 1 - y * TABLE_LG_SIZE - t) :
   //                                                  0 <= t < TABLE_LG_SIZE),
   //
   // where z[t] denotes the t-th binary digit of z (LSB is t = 0).
   // The memory footprint for _tables_ is 4KB on a 32-bit machine, or 16KB
   // on a 64-bit machine, so should fit easily into L1 cache.
   long tables[NUM_TABLES][TABLE_SIZE];
   memset(tables, 0, sizeof(long) * NUM_TABLES * TABLE_SIZE);

   long m = (p-1) / n;

   // take advantage of symmetry (n' and m' from the paper)
   if (n & 1)
      m >>= 1;
   else
      n >>= 1;

   // g^(k-1)
   long g_to_km1 = PowerMod(g, k-1, p, pinv);
   // 2^(k-1)
   long two_to_km1 = PowerMod(2, k-1, p, pinv);
   // B^(k-1), where B = 2^WORD_BITS
   long B_to_km1 = PowerMod(two_to_km1, WORD_BITS, p, pinv);
   // B^(MAX_INV)
   long s_jump = PowerMod(2, MAX_INV * WORD_BITS, p, pinv);

   // help speed up modmuls
   mulmod_precon_t g_pinv = PrepMulModPrecon(g, p, pinv);
   mulmod_precon_t g_to_km1_pinv = PrepMulModPrecon(g_to_km1, p, pinv);
   mulmod_precon_t two_to_km1_pinv = PrepMulModPrecon(two_to_km1, p, pinv);
   mulmod_precon_t B_to_km1_pinv = PrepMulModPrecon(B_to_km1, p, pinv);
   mulmod_precon_t s_jump_pinv = PrepMulModPrecon(s_jump, p, pinv);

   long g_to_km1_to_i = 1;
   long g_to_i = 1;
   long sum = 0;

   // Precompute some of the binary expansion of 1/p; at most MAX_INV words,
   // or possibly less if n is sufficiently small
   Expander expander(p, (n >= MAX_INV * WORD_BITS)
                                       ? MAX_INV : ((n - 1) / WORD_BITS + 1));

   // =========== phase 1: main summation loop

   // loop over outer sum
   for (long i = 0; i < m; i++)
   {
      // s keeps track of g^i*2^j mod p
      long s = g_to_i;
      // x keeps track of (g^i*2^j)^(k-1) mod p
      long x = g_to_km1_to_i;

      // loop over inner sum; break it up into chunks of length at most
      // MAX_INV * WORD_BITS. If n is large, this allows us to do most of
      // the work with mpn_mul_1 instead of mpn_divrem_1, and also improves
      // memory locality.
      for (long nn = n; nn > 0; nn -= MAX_INV * WORD_BITS)
      {
         word_t s_over_p[MAX_INV + 2];
         long bits, words;

         if (nn >= MAX_INV * WORD_BITS)
         {
            // do one chunk of length exactly MAX_INV * WORD_BITS
            bits = MAX_INV * WORD_BITS;
            words = MAX_INV;
         }
         else
         {
            // last chunk of length less than MAX_INV * WORD_BITS
            bits = nn;
            words = (nn - 1) / WORD_BITS + 1;
         }

         // compute some bits of the binary expansion of s/p
         expander.expand(s_over_p, s, words);
         word_t* next = s_over_p + words;

         // loop over whole words
         for (; bits >= WORD_BITS; bits -= WORD_BITS, next--)
         {
            word_t y = *next;

#if NUM_TABLES != 8 && NUM_TABLES != 4
            // generic version
            for (long h = 0; h < NUM_TABLES; h++)
            {
               long& target = tables[h][y & TABLE_MASK];
               target = SubMod(target, x, p);
               y >>= TABLE_LG_SIZE;
            }
#else
            // unrolled versions for 32-bit/64-bit machines
            long& target0 = tables[0][y & TABLE_MASK];
            target0 = SubMod(target0, x, p);

            long& target1 = tables[1][(y >> TABLE_LG_SIZE) & TABLE_MASK];
            target1 = SubMod(target1, x, p);

            long& target2 = tables[2][(y >> (2*TABLE_LG_SIZE)) & TABLE_MASK];
            target2 = SubMod(target2, x, p);

            long& target3 = tables[3][(y >> (3*TABLE_LG_SIZE)) & TABLE_MASK];
            target3 = SubMod(target3, x, p);
#if NUM_TABLES == 8
            long& target4 = tables[4][(y >> (4*TABLE_LG_SIZE)) & TABLE_MASK];
            target4 = SubMod(target4, x, p);

            long& target5 = tables[5][(y >> (5*TABLE_LG_SIZE)) & TABLE_MASK];
            target5 = SubMod(target5, x, p);

            long& target6 = tables[6][(y >> (6*TABLE_LG_SIZE)) & TABLE_MASK];
            target6 = SubMod(target6, x, p);

            long& target7 = tables[7][(y >> (7*TABLE_LG_SIZE)) & TABLE_MASK];
            target7 = SubMod(target7, x, p);
#endif
#endif

            x = MulModPrecon(x, B_to_km1, p, B_to_km1_pinv);
         }

         // loop over remaining bits in the last word
         word_t y = *next;
         for (; bits > 0; bits--)
         {
            if (y & (((word_t) 1) << (WORD_BITS - 1)))
               sum = SubMod(sum, x, p);
            else
               sum = AddMod(sum, x, p);

            x = MulModPrecon(x, two_to_km1, p, two_to_km1_pinv);
            y <<= 1;
         }

         // update s
         s = MulModPrecon(s, s_jump, p, s_jump_pinv);
      }

      // update g^i and (g^(k-1))^i
      g_to_i = MulModPrecon(g_to_i, g, p, g_pinv);
      g_to_km1_to_i = MulModPrecon(g_to_km1_to_i, g_to_km1, p, g_to_km1_pinv);
   }

   // =========== phase 2: consolidate table data

   // compute weights[z] = sum((-1)^z[t] * (2^(k-1))^(TABLE_LG_SIZE - 1 - t) :
   //                                                  0 <= t < TABLE_LG_SIZE).

   long weights[TABLE_SIZE];
   weights[0] = 0;
   for (long h = 0, x = 1; h < TABLE_LG_SIZE;
        h++, x = MulModPrecon(x, two_to_km1, p, two_to_km1_pinv))
   {
      for (long i = (1L << h) - 1; i >= 0; i--)
      {
         weights[2*i+1] = SubMod(weights[i], x, p);
         weights[2*i]   = AddMod(weights[i], x, p);
      }
   }

   // combine table data with weights

   long x_jump = PowerMod(two_to_km1, TABLE_LG_SIZE, p, pinv);

   for (long h = NUM_TABLES - 1, x = 1; h >= 0; h--)
   {
      mulmod_precon_t x_pinv = PrepMulModPrecon(x, p, pinv);

      for (long i = 0; i < TABLE_SIZE; i++)
      {
         long y = MulMod(tables[h][i], weights[i], p, pinv);
         y = MulModPrecon(y, x, p, x_pinv);
         sum = SubMod(sum, y, p);
      }

      x = MulModPrecon(x_jump, x, p, x_pinv);
   }

   return sum;
}


/******************************************************************************

   Computing the main sum (c = 1/2 case, with REDC arithmetic)

   Throughout this section F denotes 2^(ULONG_BITS / 2).

******************************************************************************/


/*
   Returns x/F mod n. Output is in [0, 2n), i.e. *not* reduced completely
   into [0, n).

   PRECONDITIONS:
      3 <= n < F, n odd
      0 <= x < nF    (if n < F/2)
      0 <= x < nF/2  (if n > F/2)
      ninv2 = -1/n mod F
*/
#define LOW_MASK ((1L << (ULONG_BITS / 2)) - 1)
static inline long RedcFast(long x, long n, long ninv2)
{
   unsigned long y = (x * ninv2) & LOW_MASK;
   unsigned long z = x + (n * y);
   return z >> (ULONG_BITS / 2);
}


/*
   Same as RedcFast(), but reduces output into [0, n).
*/
static inline long Redc(long x, long n, long ninv2)
{
   long y = RedcFast(x, n, ninv2);
   if (y >= n)
      y -= n;
   return y;
}


/*
   Computes -1/n mod F, in [0, F).

   PRECONDITIONS:
      3 <= n < F, n odd
*/
long PrepRedc(long n)
{
   long ninv2 = -n;   // already correct mod 8

   // newton's method for 2-adic inversion
   for (long bits = 3; bits < ULONG_BITS/2; bits *= 2)
      ninv2 = 2*ninv2 + n * ninv2 * ninv2;

   return ninv2 & LOW_MASK;
}


/*
   Same as bernsum_pow2(), but uses REDC arithmetic, and various delayed
   reduction strategies.

   PRECONDITIONS:
      Same as bernsum_pow2(), and in addition:
      p < 2^(ULONG_BITS/2 - 1)

   (See bernsum_pow2() for code comments; we only add comments here where
   something is different from bernsum_pow2())
*/
long bernsum_pow2_redc(long p, mulmod_t pinv, long k, long g, long n)
{
   long pinv2 = PrepRedc(p);
   long F = (1L << (ULONG_BITS/2)) % p;

   long tables[NUM_TABLES][TABLE_SIZE];
   memset(tables, 0, sizeof(long) * NUM_TABLES * TABLE_SIZE);

   long m = (p-1) / n;

   if (n & 1)
      m >>= 1;
   else
      n >>= 1;

   long g_to_km1 = PowerMod(g, k-1, p, pinv);
   long two_to_km1 = PowerMod(2, k-1, p, pinv);
   long B_to_km1 = PowerMod(two_to_km1, WORD_BITS, p, pinv);
   long s_jump = PowerMod(2, MAX_INV * WORD_BITS, p, pinv);

   long g_redc = MulMod(g, F, p, pinv);
   long g_to_km1_redc = MulMod(g_to_km1, F, p, pinv);
   long two_to_km1_redc = MulMod(two_to_km1, F, p, pinv);
   long B_to_km1_redc = MulMod(B_to_km1, F, p, pinv);
   long s_jump_redc = MulMod(s_jump, F, p, pinv);

   long g_to_km1_to_i = 1;    // always in [0, 2p)
   long g_to_i = 1;           // always in [0, 2p)
   long sum = 0;

   Expander expander(p, (n >= MAX_INV * WORD_BITS)
                                       ? MAX_INV : ((n - 1) / WORD_BITS + 1));

   // =========== phase 1: main summation loop

   for (long i = 0; i < m; i++)
   {
      long s = g_to_i;           // always in [0, p)
      if (s >= p)
         s -= p;

      long x = g_to_km1_to_i;    // always in [0, 2p)

      for (long nn = n; nn > 0; nn -= MAX_INV * WORD_BITS)
      {
         word_t s_over_p[MAX_INV + 2];
         long bits, words;

         if (nn >= MAX_INV * WORD_BITS)
         {
            bits = MAX_INV * WORD_BITS;
            words = MAX_INV;
         }
         else
         {
            bits = nn;
            words = (nn - 1) / WORD_BITS + 1;
         }

         expander.expand(s_over_p, s, words);
         word_t* next = s_over_p + words;

         for (; bits >= WORD_BITS; bits -= WORD_BITS, next--)
         {
            word_t y = *next;

            // note: we add the values into tables *without* reduction mod p

#if NUM_TABLES != 8 && NUM_TABLES != 4
            // generic version
            for (long h = 0; h < NUM_TABLES; h++)
            {
               tables[h][y & TABLE_MASK] += x;
               y >>= TABLE_LG_SIZE;
            }
#else
            // unrolled versions for 32-bit/64-bit machines
            tables[0][ y                       & TABLE_MASK] += x;
            tables[1][(y >>    TABLE_LG_SIZE ) & TABLE_MASK] += x;
            tables[2][(y >> (2*TABLE_LG_SIZE)) & TABLE_MASK] += x;
            tables[3][(y >> (3*TABLE_LG_SIZE)) & TABLE_MASK] += x;
#if NUM_TABLES == 8
            tables[4][(y >> (4*TABLE_LG_SIZE)) & TABLE_MASK] += x;
            tables[5][(y >> (5*TABLE_LG_SIZE)) & TABLE_MASK] += x;
            tables[6][(y >> (6*TABLE_LG_SIZE)) & TABLE_MASK] += x;
            tables[7][(y >> (7*TABLE_LG_SIZE)) & TABLE_MASK] += x;
#endif
#endif

            x = RedcFast(x * B_to_km1_redc, p, pinv2);
         }

         // bring x into [0, p) for next loop
         if (x >= p)
            x -= p;

         word_t y = *next;
         for (; bits > 0; bits--)
         {
            if (y & (((word_t) 1) << (WORD_BITS - 1)))
               sum = SubMod(sum, x, p);
            else
               sum = AddMod(sum, x, p);

            x = Redc(x * two_to_km1_redc, p, pinv2);
            y <<= 1;
         }

         s = Redc(s * s_jump_redc, p, pinv2);
      }

      g_to_i = RedcFast(g_to_i * g_redc, p, pinv2);
      g_to_km1_to_i = RedcFast(g_to_km1_to_i * g_to_km1_redc, p, pinv2);
   }

   // At this point, each table entry is at most p^2 (since x was always
   // in [0, 2p), and the inner loop was called at most (p/2) / WORD_BITS
   // times, and 2p * p/2 / WORD_BITS * TABLE_LG_SIZE <= p^2).

   // =========== phase 2: consolidate table data

   long weights[TABLE_SIZE];
   weights[0] = 0;
   // we store the weights multiplied by a factor of 2^(3*ULONG_BITS/2) to
   // compensate for the three rounds of REDC reduction in the loop below
   for (long h = 0, x = PowerMod(2, 3*ULONG_BITS/2, p, pinv);
        h < TABLE_LG_SIZE; h++, x = Redc(x * two_to_km1_redc, p, pinv2))
   {
      for (long i = (1L << h) - 1; i >= 0; i--)
      {
         weights[2*i+1] = SubMod(weights[i], x, p);
         weights[2*i]   = AddMod(weights[i], x, p);
      }
   }

   long x_jump = PowerMod(two_to_km1, TABLE_LG_SIZE, p, pinv);
   long x_jump_redc = MulMod(x_jump, F, p, pinv);

   for (long h = NUM_TABLES - 1, x = 1; h >= 0; h--)
   {
      for (long i = 0; i < TABLE_SIZE; i++)
      {
         long y;
         y = RedcFast(tables[h][i], p, pinv2);
         y = RedcFast(y * weights[i], p, pinv2);
         y = RedcFast(y * x, p, pinv2);
         sum += y;
      }

      x = Redc(x * x_jump_redc, p, pinv2);
   }

   return sum % p;
}



/******************************************************************************

   Wrappers for bernsum_*

******************************************************************************/


/*
   Returns B_k/k mod p, in the range [0, p).

   PRECONDITIONS:
      5 <= p < NTL_SP_BOUND, p prime
      2 <= k <= p-3, k even
      pinv = PrepMulMod(p)

   Algorithm: uses bernsum_powg() to compute the main sum.
*/
long _bern_modp_powg(long p, mulmod_t pinv, long k)
{
   Factorisation F(p-1);
   long g = primitive_root(p, pinv, F);

   // compute main sum
   long x = bernsum_powg(p, pinv, k, g);

   // divide by (1 - g^k) and multiply by 2
   long g_to_k = PowerMod(g, k, p, pinv);
   long t = InvMod(p + 1 - g_to_k, p);
   x = MulMod(x, t, p, pinv);
   x = AddMod(x, x, p);

   return x;
}


/*
   Returns B_k/k mod p, in the range [0, p).

   PRECONDITIONS:
      5 <= p < NTL_SP_BOUND, p prime
      2 <= k <= p-3, k even
      pinv = PrepMulMod(p)
      2^k != 1 mod p

   Algorithm: uses bernsum_pow2() (or bernsum_pow2_redc() if p is small
   enough) to compute the main sum.
*/
long _bern_modp_pow2(long p, mulmod_t pinv, long k)
{
   Factorisation F(p-1);
   long g = primitive_root(p, pinv, F);
   long n = order(2, p, pinv, F);

   // compute main sum
   long x;
   if (p < (1L << (ULONG_BITS/2 - 1)))
      x = bernsum_pow2_redc(p, pinv, k, g, n);
   else
      x = bernsum_pow2(p, pinv, k, g, n);

   // divide by 2*(2^(-k) - 1)
   long t = PowerMod(2, -k, p, pinv) - 1;
   t = AddMod(t, t, p);
   t = InvMod(t, p);
   x = MulMod(x, t, p, pinv);

   return x;
}



/*
   Returns B_k/k mod p, in the range [0, p).

   PRECONDITIONS:
      5 <= p < NTL_SP_BOUND, p prime
      2 <= k <= p-3, k even
      pinv = PrepMulMod(p)
*/
long _bern_modp(long p, mulmod_t pinv, long k)
{
   if (PowerMod(2, k, p, pinv) != 1)
      // 2^k != 1 mod p, so we use the faster version
      return _bern_modp_pow2(p, pinv, k);
   else
      // forced to use slower version
      return _bern_modp_powg(p, pinv, k);
}



/******************************************************************************

   Main bern_modp() routine

******************************************************************************/

long bern_modp(long p, long k)
{
   assert(k >= 0);
   assert(2 <= p && p < NTL_SP_BOUND);

   // B_0 = 1
   if (k == 0)
      return 1;

   // B_1 = -1/2 mod p
   if (k == 1)
   {
      if (p == 2)
         return -1;
      return (p-1)/2;
   }

   // B_k = 0 for odd k >= 3
   if (k & 1)
      return 0;

   // denominator of B_k is always divisible by 6 for k >= 2
   if (p <= 3)
      return -1;

   // use Kummer's congruence (k = m mod p-1  =>  B_k/k = B_m/m mod p)
   long m = k % (p-1);
   if (m == 0)
      return -1;

   mulmod_t pinv = PrepMulMod(p);
   long x = _bern_modp(p, pinv, m);    // = B_m/m mod p
   return MulMod(x, k%p, p, pinv);
}


};    // end namespace



// end of file ================================================================
