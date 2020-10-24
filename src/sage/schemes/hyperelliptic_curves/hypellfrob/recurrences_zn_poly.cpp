/* ============================================================================

   recurrences_zn_poly.cpp:  recurrences solved via zn_poly arithmetic

   This file is part of hypellfrob (version 2.1.1).

   Copyright (C) 2007, 2008, David Harvey

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

============================================================================ */


#include "recurrences_zn_poly.h"


NTL_CLIENT


namespace hypellfrob {


Shifter::~Shifter()
{
   zn_array_mulmid_precomp1_clear(kernel_precomp);
   free(input_twist);
}


Shifter::Shifter(ulong d, ulong a, ulong b, const zn_mod_t mod)
{
   this->d = d;
   this->mod = mod;

   input_twist = (ulong*) malloc(sizeof(ulong) * (3*d + 3));
   output_twist = input_twist + d + 1;
   scratch = output_twist + d + 1;

   ZZ modulus;
   modulus = zn_mod_get (mod);

   // ------------------------ compute input_twist -------------------------

   // prod = (d!)^(-1)
   ulong prod = 1;
   for (ulong i = 2; i <= d; i++)
      prod = zn_mod_mul(prod, i, mod);
   prod = to_ulong(InvMod(to_ZZ(prod), modulus));

   // input_twist[i] = ((d-i)!)^(-1)
   input_twist[0] = prod;
   for (ulong i = 1; i <= d; i++)
      input_twist[i] = zn_mod_mul(input_twist[i-1], d - (i-1), mod);

   // input_twist[i] = ((d-i)!*i!)^(-1)
   for (ulong i = 0; i <= d/2; i++)
   {
      input_twist[i] = zn_mod_mul(input_twist[i], input_twist[d-i], mod);
      input_twist[d-i] = input_twist[i];
   }

   // input_twist[i] = \prod_{0 <= j <= d, j != i} (i-j)^(-1)
   //                = (-1)^(d-i) ((d-i)!*i!)^(-1)
   for (long i = d - 1; i >= 0; i -= 2)
      input_twist[i] = zn_mod_neg(input_twist[i], mod);

   // ----------------- compute output_twist and kernel --------------------

   // need some temp space:
   // c, accum, accum_inv each of length 2d+1
   ulong* c = (ulong*) malloc(sizeof(ulong) * (6*d+3));
   ulong* accum = c + 2*d + 1;
   ulong* accum_inv = accum + 2*d + 1;
   ulong* kernel = c;   // overwrites c

   // c[i] = c_i = a + (i-d)*b     for 0 <= i <= 2d
   c[0] = zn_mod_sub(a, zn_mod_mul(zn_mod_reduce(d, mod), b, mod), mod);
   for (ulong i = 1; i <= 2*d; i++)
      c[i] = zn_mod_add(c[i-1], b, mod);

   // accum[i] = c_0 * c_1 * ... * c_i    for 0 <= i <= 2d
   accum[0] = c[0];
   for (ulong i = 1; i <= 2*d; i++)
      accum[i] = zn_mod_mul(accum[i-1], c[i], mod);

   // accum_inv[i] = (c_0 * c_1 * ... * c_i)^(-1)    for 0 <= i <= 2d
   accum_inv[2*d] = to_ulong(InvMod(to_ZZ(accum[2*d]), modulus));

   for (long i = 2*d - 1; i >= 0; i--)
      accum_inv[i] = zn_mod_mul(accum_inv[i+1], c[i+1], mod);

   // output_twist[i] = b^{-d} * c_i * c_{i+1} * ... * c_{i+d}
   // for 0 <= i <= d
   ulong factor = to_long(PowerMod(to_ZZ(b), -((long)d), modulus));
   output_twist[0] = zn_mod_mul(factor, accum[d], mod);
   for (ulong i = 1; i <= d; i++)
      output_twist[i] = zn_mod_mul(zn_mod_mul(factor, accum[i+d], mod),
                                   accum_inv[i-1], mod);

   // kernel[i] = (c_i)^(-1)    for 0 <= i <= 2d
   kernel[0] = accum_inv[0];
   for (ulong i = 1; i <= 2*d; i++)
      kernel[i] = zn_mod_mul(accum_inv[i], accum[i-1], mod);

   // precompute FFT of kernel
   zn_array_mulmid_precomp1_init(kernel_precomp, kernel, 2*d+1, d+1, mod);

   free(c);
}


void Shifter::shift(ulong* output, const ulong* input)
{
   // multiply inputs pointwise by input_twist
   for (ulong i = 0; i <= d; i++)
      scratch[i] = zn_mod_mul(input[i], input_twist[i], mod);

   // do middle product
   zn_array_mulmid_precomp1_execute(output, scratch, kernel_precomp);

   // multiply outputs pointwise by output_twist
   for (ulong i = 0; i <= d; i++)
      output[i] = zn_mod_mul(output[i], output_twist[i], mod);
}


/*
   Checks whether large_evaluate() will succeed for these choices of k and u.
   Returns 1 if okay, otherwise 0.
*/
int check_params(ulong k, ulong u, const zn_mod_t mod)
{
   ulong n = zn_mod_get (mod);

   if (k >= n || u >= n)
      return 0;

   if (k <= 1)
      return 1;

   if (k == n - 1)
      return 0;

   ulong k2 = k / 2;

   // need the following elements to be invertible:
   // u
   // 1, 2, ..., k + 1
   // k2 + i*u for -k2 <= i <= k2
   ulong prod = u;
   for (ulong i = 2; i <= k; i++)
      prod = zn_mod_mul(prod, i, mod);
   ulong temp = zn_mod_mul(k2, zn_mod_sub(1, u, mod), mod);
   for (ulong i = 0; i <= 2*k2; i++)
   {
      prod = zn_mod_mul(prod, temp, mod);
      temp = zn_mod_add(temp, u, mod);
   }

   ZZ x, y;
   x = prod;
   y = n;
   if (GCD(x, y) != 1)
      return 0;

   // check recursively below
   return check_params(k2, u, mod);
}



/*
   Let M0 and M1 be square matrices of size r*r. Let M(x) = M0 + x*M1; this
   is a matrix of linear polys in x. The matrices M0 and M1 are passed in
   row-major order.

   Let P(x) = M(x+1) M(x+2) ... M(x+k); this is a matrix of polynomials of
   degree k.

   This class attempts to compute the matrices
       P(0), P(u), P(2u), ..., P(ku).

   Usage:

   * Call the constructor.

   * Call evaluate() with half == 0. This computes the first half of the
   outputs, specifically P(mu) for 0 <= m <= k2, where k2 = floor(k/2).
   The (i, j) entry of P(mu) is stored in output[i*r+j][m+offset]. Each array
   output[i*r+j] must be preallocated to length k2 + 3. (The extra two matrices
   at the end are used for scratch space.)

   * Call evaluate() with half == 1. This computes the second half, namely
   P(mu) for k2 + 1 <= m <= k. The (i, j) entry of P(mu) is stored in
   output[i*r+j][m-(k2+1)+offset]. Each array output[i*r+j] must be
   preallocated to length k2 + 3. It's okay for this second half to overwrite
   the first half generated earlier (in fact this property is used by
   zn_poly_interval_products to conserve memory).

   The computation may fail for certain bad choices of k and u. Let p = residue
   characteristic (i.e. assume mod is p^m for some m). Typically this function
   gets called for k ~ u and k*u ~ M*p, for some M much smaller than p. In this
   situation, I expect most choices of (k, u) are not bad. Failure is very
   bad: the program will crash. Luckily, you can call check_params() to test
   whether (k, u) is bad. If it's bad, you probably should just increment u
   until you find one that's not bad. (Sorry, I haven't done a proper analysis
   of the situation yet.) The main routines fall back on the zz_pX version if
   they can't find a good parameter choice.

   Must have 0 <= k < n and 0 < u < n.

*/
LargeEvaluator::LargeEvaluator(int r, ulong k, ulong u,
                               const vector<vector<ulong> >& M0,
                               const vector<vector<ulong> >& M1,
                               const zn_mod_t& mod) : M0(M0), M1(M1), mod(mod)
{
   assert(k < zn_mod_get(mod));
   assert(u < zn_mod_get(mod));
   assert(r >= 1);
   assert(k >= 1);

   this->r = r;
   this->k = k;
   this->k2 = k / 2;
   this->odd = k & 1;
   this->u = u;
   this->shifter = NULL;
}


LargeEvaluator::~LargeEvaluator()
{
   if (this->shifter)
      delete shifter;
}


void LargeEvaluator::evaluate(int half, vector<ulong_array>& output,
                              ulong offset)
{
   // base cases

   if (k == 0)
   {
      // identity matrix
      for (int x = 0; x < r; x++)
      for (int y = 0; y < r; y++)
         output[y*r + x].data[offset] = (x == y);
      return;
   }

   if (k == 1)
   {
      // evaluate M(1) (for half == 0) or M(u+1) (for half == 1)
      for (int x = 0; x < r; x++)
      for (int y = 0; y < r; y++)
      {
         ulong temp = zn_mod_add(M0[y][x], M1[y][x], mod);
         output[y*r + x].data[offset] = half ?
               zn_mod_add(temp, zn_mod_mul(u, M1[y][x], mod), mod) : temp;
      }
      return;
   }

   // recursive case

   // Let Q(x) = M(x+1) ... M(x+k2), where k2 = floor(k/2).
   // Then we have either
   //    P(x) = Q(x) Q(x+k2)           if k is even
   //    P(x) = Q(x) Q(x+k2) M(x+k)    if k is odd

   if (half == 0)
   {
      // This is the first time through this function.
      // Allocate scratch space.
      scratch.resize(r*r);
      for (int i = 0; i < r*r; i++)
         scratch[i].resize(k2 + 3);

      {
         // Recursively compute Q(0), Q(u), ..., Q(k2*u).
         // These are stored in scratch.
         LargeEvaluator recurse(r, k2, u, M0, M1, mod);
         recurse.evaluate_all(scratch);
      }

      // Precomputations for value-shifting
      shifter = new Shifter(k2, k2, u, mod);
   }
   else   // half == 1
   {
      // Shift original sequence by (k2+1)*u to obtain
      // Q((k2+1)*u), Q((k2+2)*u), ..., Q((2*k2+1)*u)
      Shifter big_shifter(k2, zn_mod_mul(k2+1, u, mod), u, mod);
      for (int i = 0; i < r*r; i++)
         big_shifter.shift(scratch[i].data, scratch[i].data);
   }

   // Let H = (k2+1)*u*half, so now scratch contains
   // Q(H), Q(H + u), ..., Q(H + k2*u).

   // Shift by k2 to obtain Q(H + k2), Q(H + u + k2), ..., Q(H + k2*u + k2).
   // Results are stored directly in output array. We put them one slot
   // to the right, to make room for inplace matrix multiplies later on.
   // (If k is odd, they're shifted right by yet another slot, to make room
   // for another round of multiplies.)
   for (int i = 0; i < r*r; i++)
      shifter->shift(output[i].data + offset + odd + 1, scratch[i].data);

   // If k is odd, right-multiply each Q(H + i*u + k2) by M(H + i*u + k)
   // (results are stored in output array, shifted one entry to the left)
   if (odd)
   {
      ulong_array cruft(r*r);    // for storing M(H + i*u + k)
      ulong point = k;           // evaluation point
      if (half)
         point = zn_mod_add(point, zn_mod_mul(k2 + 1, u, mod), mod);

      for (ulong i = 0; i <= k2; i++, point = zn_mod_add(point, u, mod))
      {
         // compute M(H + i*u + k) = M0 + M1*point
         for (int x = 0; x < r; x++)
         for (int y = 0; y < r; y++)
            cruft.data[y*r + x] = zn_mod_add(M0[y][x],
                                  zn_mod_mul(M1[y][x], point, mod), mod);

         // multiply
         for (int x = 0; x < r; x++)
         for (int y = 0; y < r; y++)
         {
            ulong accum = 0;
            for (int z = 0; z < r; z++)
               accum = zn_mod_add(accum,
                       zn_mod_mul(output[y*r + z].data[offset + i + 2],
                                  cruft.data[z*r + x], mod), mod);
            output[y*r + x].data[offset + i + 1] = accum;
         }
      }
   }

   ulong n = zn_mod_get(mod);

   // Multiply to obtain P(H), P(H + u), ..., P(H + k2*u)
   // (except for the last one, in the second half, if k is even)
   // Store results directly in output array.
   for (ulong i = 0; i + (half && !odd) <= k2; i++)
   for (int x = 0; x < r; x++)
   for (int y = 0; y < r; y++)
   {
      ulong sum_hi = 0;
      ulong sum_lo = 0;
      for (int z = 0; z < r; z++)
      {
         ulong hi, lo;
         ZNP_MUL_WIDE(hi, lo, scratch[y*r + z].data[i],
                              output[z*r + x].data[offset + i + 1]);
         ZNP_ADD_WIDE(sum_hi, sum_lo, sum_hi, sum_lo, hi, lo);
         if (sum_hi >= n)
            sum_hi -= n;
      }
      output[y*r + x].data[offset + i] =
                           zn_mod_reduce_wide(sum_hi, sum_lo, mod);
   }
}


/*
   Evaluates both halves, storing results in output array.
*/
void LargeEvaluator::evaluate_all(vector<ulong_array>& output)
{
   evaluate(0, output, 0);
   evaluate(1, output, k/2 + 1);
}



/*
See interval_products_wrapper().

NOTE 1:
   I haven't proved that the algorithm here always succeeds, although I expect
   that it almost always will, especially when p is very large. If it
   succeeds, it returns 1. If it fails, it returns 0 (practically instantly).
   In the latter case, at least the caller can fall back on
   ntl_interval_products().

NOTE 2:
   The algorithm here is similar to ntl_interval_products(). However, it
   doesn't do the "refining step" -- it just handles the smaller intervals
   in the naive fashion. Also, instead of breaking intervals into power-of-four
   lengths, it just does the whole thing in one chunk. The performance ends up
   being smoother, but it's harder to prove anything about avoiding
   non-invertible elements. Hence the caveat in Note 1.

*/
int zn_poly_interval_products(vector<vector<vector<ulong> > >& output,
                              const vector<vector<ulong> >& M0,
                              const vector<vector<ulong> >& M1,
                              const vector<ZZ>& target, const zn_mod_t& mod)
{
   output.resize(target.size() / 2);

   // select k such that k*(k+1) >= total length of intervals
   ulong k;
   {
      ZZ len = target.back() - target.front();
      ZZ kk = SqrRoot(len);
      k = to_ulong(kk);
      if (kk * (kk + 1) < len)
         k++;
   }

   int r = M0.size();

   // try to find good parameters u and k
   ulong u = k;
   for (int trial = 0; ; trial++, u++)
   {
      if (check_params(k, u, mod))
         break;    // found some good parameters
      if (trial == 5)
         return 0;    // too many failures, give up
   }

   ulong n = zn_mod_get(mod);

   // shift M0 over to account for starting index
   vector<vector<ulong> > M0_shifted = M0;
   ulong shift = target.front() % n;
   for (int x = 0; x < r; x++)
   for (int y = 0; y < r; y++)
      M0_shifted[y][x] = zn_mod_add(M0[y][x],
                            zn_mod_mul(shift, M1[y][x], mod), mod);

   // prepare for evaluating products over the big intervals
   // [0, k), [u, u+k), ..., [ku, ku+k)
   LargeEvaluator evaluator(r, k, u, M0_shifted, M1, mod);

   vector<ulong_array> big(r*r);
   for (int i = 0; i < r*r; i++)
      big[i].resize(k/2 + 3);    // space for half the products, plus two more

   // evaluate the first half of the products
   evaluator.evaluate(0, big, 0);
   // flag indicating which half we currently have stored in the "big" array
   int half = 0;

   vector<vector<ulong> > accum(r, vector<ulong>(r));
   vector<vector<ulong> > temp1(r, vector<ulong>(r));
   vector<vector<ulong> > temp2(r, vector<ulong>(r));

   // for each target interval....
   for (size_t i = 0; i < target.size()/2; i++)
   {
      // doing interval [s0, s1)
      ZZ s0 = target[2*i];
      ZZ s1 = target[2*i+1];

      // product accumulated so far is [t0, t1).
      ZZ t0 = s0;
      ZZ t1 = s0;
      for (int x = 0; x < r; x++)
      for (int y = 0; y < r; y++)
         accum[x][y] = (x == y);     // identity matrix

      while (t1 < s1)
      {
         // if we are exactly on the left end of a big interval, and the
         // big interval fits inside the target interval, then roll it in
         if ((t1 - target[0]) % u == 0  &&  t1 + k <= s1)
         {
            // compute which "big" interval we are rolling in
            size_t index = to_ulong((t1 - target[0]) / u);

            if (index >= k/2 + 1)
            {
               // if the "big" interval is in the second half, and we haven't
               // computed the second half of the intervals yet, then go and
               // compute them (overwriting the first half)
               if (half == 0)
               {
                  evaluator.evaluate(1, big, 0);
                  half = 1;
               }
               index -= (k/2 + 1);
            }

            for (int x = 0; x < r; x++)
            for (int y = 0; y < r; y++)
               temp1[y][x] = big[y*r + x].data[index];

            t1 += k;
         }
         else
         {
            // otherwise just multiply by the single matrix M(t1 + 1)
            ulong e = (t1 + 1) % n;

            for (int x = 0; x < r; x++)
            for (int y = 0; y < r; y++)
               temp1[y][x] = zn_mod_add(M0[y][x],
                             zn_mod_mul(M1[y][x], e, mod), mod);

            t1++;
         }

         // multiply by whichever matrix we picked up above
         for (int x = 0; x < r; x++)
         for (int y = 0; y < r; y++)
         {
            ulong sum_hi = 0;
            ulong sum_lo = 0;
            for (int z = 0; z < r; z++)
            {
               ulong hi, lo;
               ZNP_MUL_WIDE(hi, lo, accum[y][z], temp1[z][x]);
               ZNP_ADD_WIDE(sum_hi, sum_lo, sum_hi, sum_lo, hi, lo);
               if (sum_hi >= n)
                  sum_hi -= n;
            }
            temp2[y][x] = zn_mod_reduce_wide(sum_hi, sum_lo, mod);
         }

         accum.swap(temp2);
      }

      // store result in output array
      output[i] = accum;
   }

   return 1;
}


};  // namespace hypellfrob


// ----------------------- end of file
