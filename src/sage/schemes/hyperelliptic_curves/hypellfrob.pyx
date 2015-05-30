r"""
Frobenius on Monsky-Washnitzer cohomology of a hyperelliptic curve over GF(p),
for largish p

This is a wrapper for the matrix() function in hypellfrob.cpp.

AUTHOR:
    -- David Harvey (2007-05)
    -- David Harvey (2007-12): rewrote for hypellfrob version 2.0
"""

#################################################################################
#       Copyright (C) 2007 David Harvey <dmharvey@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX
from sage.libs.ntl.ntl_mat_ZZ cimport ntl_mat_ZZ
from sage.libs.ntl.all import ZZ, ZZX
from sage.matrix.all import Matrix
from sage.rings.all import Qp, O as big_oh
from sage.rings.arith import is_prime

include "sage/libs/ntl/decl.pxi"
include "sage/ext/interrupt.pxi"


cdef extern from "hypellfrob.h":
     int hypellfrob_matrix "hypellfrob::matrix" (mat_ZZ_c output, ZZ_c p, int N, ZZX_c Q)


def hypellfrob(p, N, Q):
   r"""
   Computes the matrix of Frobenius acting on the Monsky-Washnitzer cohomology
   of a hyperelliptic curve $y^2 = Q(x)$, with respect to the basis $x^i dx/y$,
   $0 \leq i < 2g$.

   INPUT:
      p -- a prime
      Q -- a monic polynomial in $\Z[x]$ of odd degree. Must have no multiple
           roots mod p.
      N -- precision parameter; the output matrix will be correct modulo $p^N$.

   PRECONDITIONS:
      Must have $p > (2g+1)(2N-1)$, where $g = (\deg(Q)-1)/2$ is the genus
      of the curve.

   ALGORITHM:
      Described in ``Kedlaya's algorithm in larger characteristic'' by David
      Harvey. Running time is theoretically soft-$O(p^{1/2} N^{5/2} g^3)$.

   TODO:
      Remove the restriction on $p$. Probably by merging in Robert's code,
      which eventually needs a fast C++/NTL implementation.

   EXAMPLES:
      sage: from sage.schemes.hyperelliptic_curves.hypellfrob import hypellfrob
      sage: R.<x> = PolynomialRing(ZZ)
      sage: f = x^5 + 2*x^2 + x + 1; p = 101
      sage: M = hypellfrob(p, 4, f); M
       [ 91844754 + O(101^4)  38295665 + O(101^4)  44498269 + O(101^4)  11854028 + O(101^4)]
       [ 93514789 + O(101^4)  48987424 + O(101^4)  53287857 + O(101^4)  61431148 + O(101^4)]
       [ 77916046 + O(101^4)  60656459 + O(101^4) 101244586 + O(101^4)  56237448 + O(101^4)]
       [ 58643832 + O(101^4)  81727988 + O(101^4)  85294589 + O(101^4)  70104432 + O(101^4)]
      sage: -M.trace()
       7 + O(101^4)
      sage: sum([legendre_symbol(f(i), p) for i in range(p)])
       7
      sage: ZZ(M.det())
       10201
      sage: M = hypellfrob(p, 1, f); M
       [ 0 + O(101)  0 + O(101) 93 + O(101) 62 + O(101)]
       [ 0 + O(101)  0 + O(101) 55 + O(101) 19 + O(101)]
       [ 0 + O(101)  0 + O(101) 65 + O(101) 42 + O(101)]
       [ 0 + O(101)  0 + O(101) 89 + O(101) 29 + O(101)]

   AUTHORS:
      -- David Harvey (2007-05)
      -- David Harvey (2007-12): updated for hypellfrob version 2.0

   """

   # Sage objects that wrap the NTL objects
   cdef ntl_ZZ pp
   cdef ntl_ZZX QQ
   cdef ntl_mat_ZZ mm   # the result will go in mm
   cdef int i, j

   if N < 1:
      raise ValueError, "N must be an integer >= 1"

   Q = Q.list()
   if len(Q) < 4 or len(Q) % 2 or Q[-1] != 1:
      raise ValueError, "Q must be a monic polynomial of odd degree >= 3"
   QQ = ZZX(Q)

   bound = (len(Q) - 1) * (2*N - 1)
   if p <= bound:
      raise ValueError, "In the current implementation, p must be greater than (2g+1)(2N-1) = %s" % bound

   if not is_prime(p):
      raise ValueError, "p (= %s) must be prime" % p

   pp = ZZ(p)

   cdef int g    # the genus
   g = (len(Q) / 2) - 1

   # Note: the C++ code actually resets the size of the matrix, but this seems
   # to confuse the Sage NTL wrapper. So to be safe I'm setting it ahead of
   # time.
   mm = ntl_mat_ZZ(2*g, 2*g)

   cdef int result
   sig_on()
   cdef mat_ZZ_c *mm_x = &mm.x    # workaround for Cython misfeature
   result = hypellfrob_matrix(mm_x[0], pp.x, N, QQ.x)
   sig_off()

   if not result:
      raise ValueError, "Could not compute frobenius matrix, because the curve is singular at p."

   R = Qp(p, N, print_mode="terse")
   prec = big_oh(p**N)
   data = [[mm[j, i]._integer_() + prec for i from 0 <= i < 2*g] for j from 0 <= j < 2*g]
   return Matrix(R, data)


################ end of file
