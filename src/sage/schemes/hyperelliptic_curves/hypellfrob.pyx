# distutils: language = c++
# distutils: sources = sage/schemes/hyperelliptic_curves/hypellfrob/hypellfrob.cpp sage/schemes/hyperelliptic_curves/hypellfrob/recurrences_ntl.cpp
# distutils: depends = sage/schemes/hyperelliptic_curves/hypellfrob/hypellfrob.h sage/schemes/hyperelliptic_curves/hypellfrob/recurrences_ntl.h
# distutils: include_dirs = sage/libs/ntl/ sage/schemes/hyperelliptic_curves/hypellfrob/ NTL_INCDIR
# distutils: libraries = gmp NTL_LIBRARIES
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# sage.doctest: needs sage.libs.ntl sage.modules sage.rings.padics

r"""
Frobenius on Monsky-Washnitzer cohomology of a hyperelliptic curve

This module provides :func:`hypellfrob`, that is a wrapper for the ``matrix()``
function in ``hypellfrob.cpp``.

``hypellfrob.cpp`` is a C++ program for computing the zeta function of a
hyperelliptic curve over a largish prime finite field, based on the method
described in the paper [Harv2007]_. More precisely, it computes the matrix of
Frobenius on the Monsky-Washnitzer cohomology of the curve; the zeta function
can be recovered via the characteristic polynomial of the matrix.

AUTHORS:

- David Harvey (2007-05): initial version
- David Harvey (2007-12): rewrote for ``hypellfrob`` version 2.0
- Alex J. Best (2018-02): added wrapper
"""

# *****************************************************************************
#       Copyright (C) 2007 David Harvey <dmharvey@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from cysignals.signals cimport sig_on, sig_off
from libcpp.vector cimport vector

from sage.libs.ntl.ntl_ZZ_pContext import ZZ_pContext_factory
from sage.libs.ntl.all import ZZ, ZZX
from sage.matrix.constructor import Matrix
from sage.rings.padics.factory import Qp
from sage.rings.big_oh import O as big_oh
from sage.arith.misc import is_prime

from sage.libs.ntl.ntl_ZZ_p cimport ntl_ZZ_p
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX
from sage.libs.ntl.ntl_mat_ZZ cimport ntl_mat_ZZ
from sage.libs.ntl.ntl_ZZ_pContext cimport (ntl_ZZ_pContext_class,
                                            ntl_ZZ_pContext_factory)

from sage.libs.ntl.conversion cimport set_ntl_matrix_modn_dense

include "sage/libs/ntl/decl.pxi"


cdef extern from "hypellfrob.h":
    int hypellfrob_matrix "hypellfrob::matrix" (mat_ZZ_c output, ZZ_c p,
                                                int N, ZZX_c Q)
    void interval_products_wrapper \
        "hypellfrob::hypellfrob_interval_products_wrapper" \
        (mat_ZZ_p_c &output, const mat_ZZ_p_c &M0, const mat_ZZ_p_c &M1,
         const vector[ZZ_c] target)


def interval_products(M0, M1, target):
    r"""
    Given matrices `M(t)` with entries linear in `t` over `\ZZ/N\ZZ` and a list
    of integers `a_0 < b_0 \le a_1 < b_1 \le \cdots \le a_n < b_n`, compute the
    matrices `\prod_{t = a_i + 1}^{b_i} M(t)` for `i = 0` to `n`.

    INPUT:

    - ``M0``, ``M1`` -- matrices over `\ZZ/N\ZZ`, so that `M(t) = M_0 + M_1t`
    - ``target`` -- list of integers `a_0, b_0, \dots, a_n, b_n`

    ALGORITHM:

    Described in [Harv2007]_, Theorem 10. Based on the work of Bostan-Gaudry-Schost
    [BGS2007]_.

    EXAMPLES::

        sage: from sage.schemes.hyperelliptic_curves.hypellfrob import interval_products
        sage: interval_products(Matrix(Integers(9), 2,2, [1,0,1,0]),
        ....:   Matrix(Integers(9), 2, 2, [1, 1, 0, 2]),[0,2,2,4])
        [
        [7 8]  [5 4]
        [5 1], [2 7]
        ]
        sage: [prod(Matrix(Integers(9), 2, 2, [t + 1, t, 1, 2*t])
        ....:  for t in range(2*i + 1, 2*i + 1 + 2)) for i in range(2)]
        [
        [7 8]  [5 4]
        [5 1], [2 7]
        ]

    An example with larger modulus::

        sage: interval_products(Matrix(Integers(3^8), 1, 1, [1]),
        ....:   Matrix(Integers(3^8), 1, 1, [1]), [2,4])
        [[20]]
        sage: [prod(Matrix(Integers(3^8), 1, 1, [t + 1]) for t in range(3,5))]
        [[20]]

    An even larger modulus::

        sage: interval_products(Matrix(Integers(3^18), 1, 1, [1]),
        ....:   Matrix(Integers(3^18), 1, 1, [1]), [2,4])
        [[20]]
        sage: [prod(Matrix(Integers(3^18), 1, 1, [t + 1]) for t in range(3,5))]
        [[20]]
    """
    # Sage objects that wrap the NTL objects
    cdef mat_ZZ_p_c mm0, mm1
    cdef mat_ZZ_p_c out
    cdef vector[ZZ_c] targ
    cdef ntl_ZZ_pContext_class c = \
        (<ntl_ZZ_pContext_factory>ZZ_pContext_factory).make_c(
            ntl_ZZ(M0.base_ring().characteristic()))
    cdef long dim = M0.nrows()
    sig_on()
    c.restore_c()
    sig_off()
    set_ntl_matrix_modn_dense(mm0, c, M0)
    set_ntl_matrix_modn_dense(mm1, c, M1)
    sig_on()
    for t in target:
        targ.push_back(ntl_ZZ(t).x)
    numintervals = len(target)/2
    out.SetDims(dim, dim*numintervals)

    interval_products_wrapper(out, mm0, mm1, targ)
    sig_off()

    R = M0.matrix_space()
    mats = [R(0) for k in range(numintervals)]
    cdef ntl_ZZ_p tmp
    tmp = ntl_ZZ_p(modulus=c)
    for k in range(numintervals):
        for j in range(dim):
            for i in range(dim):
                sig_on()
                tmp.x = out.get(j, i + dim * k)
                sig_off()
                mats[k][j, i] = tmp._integer_()

    return mats


def hypellfrob(p, N, Q):
    r"""
    Compute the matrix of Frobenius acting on the Monsky-Washnitzer cohomology
    of a hyperelliptic curve `y^2 = Q(x)`, with respect to the basis
    `x^i dx/y`, `0 \leq i < 2g`.

    INPUT:

    - ``p`` -- a prime
    - ``Q`` -- a monic polynomial in `\ZZ[x]` of odd degree; must have no
      multiple roots mod `p`
    - ``N`` -- precision parameter; the output matrix will be correct modulo `p^N`

    The prime `p` should satisfy `p > (2g+1)(2N-1)`, where `g =
    \left(\deg Q - 1\right) / 2` is the genus of the curve.

    ALGORITHM:

    Described in [Harv2007]_, Section 7. Running time is theoretically
    `\widetilde{O}(p^{1/2} N^{5/2} g^3)`.

    EXAMPLES::

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
        sage: sum(legendre_symbol(f(i), p) for i in range(p))
        7
        sage: ZZ(M.det())
        10201
        sage: M = hypellfrob(p, 1, f); M
        [     O(101)      O(101) 93 + O(101) 62 + O(101)]
        [     O(101)      O(101) 55 + O(101) 19 + O(101)]
        [     O(101)      O(101) 65 + O(101) 42 + O(101)]
        [     O(101)      O(101) 89 + O(101) 29 + O(101)]

    .. TODO::

        Remove the restriction on `p`. Probably by merging in Robert's code,
        which eventually needs a fast C++/NTL implementation.
    """
    # Sage objects that wrap the NTL objects
    cdef ntl_ZZ pp
    cdef ntl_ZZX QQ
    cdef ntl_mat_ZZ mm   # the result will go in mm
    cdef int i, j

    if N < 1:
        raise ValueError("N must be an integer >= 1")

    Q = Q.list()
    if len(Q) < 4 or len(Q) % 2 or Q[-1] != 1:
        raise ValueError("Q must be a monic polynomial of odd degree >= 3")
    QQ = ZZX(Q)

    bound = (len(Q) - 1) * (2*N - 1)
    if p <= bound:
        raise ValueError("In the current implementation, p must be greater "
                         "than (2g+1)(2N-1) = %s" % bound)

    if not is_prime(p):
        raise ValueError("p (= %s) must be prime" % p)

    pp = ZZ(p)

    cdef int g    # the genus
    g = (len(Q) / 2) - 1

    # Note: the C++ code actually resets the size of the matrix, but this seems
    # to confuse the Sage NTL wrapper. So to be safe I'm setting it ahead of
    # time.
    mm = ntl_mat_ZZ(2 * g, 2 * g)

    cdef int result
    sig_on()
    cdef mat_ZZ_c *mm_x = &mm.x    # workaround for Cython misfeature
    result = hypellfrob_matrix(mm_x[0], pp.x, N, QQ.x)
    sig_off()

    if not result:
        raise ValueError("Could not compute frobenius matrix"
                         ", because the curve is singular at p.")

    R = Qp(p, N, print_mode='terse')
    prec = big_oh(p**N)
    data = [[mm[j, i]._integer_() + prec for i in range(2 * g)]
            for j in range(2 * g)]
    return Matrix(R, data)
