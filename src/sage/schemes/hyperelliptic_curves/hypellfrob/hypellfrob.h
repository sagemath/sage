/* ============================================================================

   hypellfrob.h:  header for hypellfrob.cpp

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


#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <vector>


namespace hypellfrob {


/*

Let M(x) be the matrix M0 + x*M1; this is a matrix of linear polys in x.
Let M(a, b) = M(a + 1) M(a + 2) ... M(b). This function evaluates the products
M(a[i], b[i]) for some sequence of intervals
  a[0] < b[0] <= a[1] < b[1] <= ... <= a[n-1] < b[n-1].

The intervals are supplied in "target", simply as the list
  a[0], b[0], a[1], b[1], ...

There are three possible underlying implementations:
   * ntl_interval_products (ZZ_p version),
   * ntl_interval_products (zz_p version)
This function is a wrapper which takes ZZ_p input, calls one of the two
above implementations depending on the size of the current ZZ_p modulus, and
produces output in ZZ_p format.

PRECONDITIONS:
   Let d = b[n-1] - a[0]. Then 2, 3, ... 1 + floor(sqrt(d)) must all be
   invertible.

*/
void hypellfrob_interval_products_wrapper(NTL::mat_ZZ_p& output,
                               const NTL::mat_ZZ_p& M0, const NTL::mat_ZZ_p& M1,
                               const std::vector<NTL::ZZ>& target);

/*
Computes frobenius matrix for given p, to precision p^N, for the
hyperelliptic curve y^2 = Q(x), on the standard basis of cohomology.

PRECONDITIONS:
   p must be a prime > (2g+1)(2N-1).
   N >= 1.
   Degree of Q should be 2g+1 for some g >= 1.
   Q must be monic. The reduction of Q mod p must have no multiple roots.

RETURN VALUE:
   1 on success, in which case "output" holds the resulting 2g * 2g matrix.
   0 if any of the above conditions are not satisfied. (EXCEPTION: matrix()
       will not check that p is prime. That's up to you.)

*/
int matrix(NTL::mat_ZZ& output, const NTL::ZZ& p, int N, const NTL::ZZX& Q);


};

// ----------------------- end of file
