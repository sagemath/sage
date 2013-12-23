r"""
Bernoulli numbers modulo p

AUTHOR:
    - David Harvey (2006-07-26): initial version
    - William Stein (2006-07-28): some touch up.
    - David Harvey (2006-08-06): new, faster algorithm, also using faster NTL interface
    - David Harvey (2007-08-31): algorithm for a single Bernoulli number mod p
    - David Harvey (2008-06): added interface to bernmm, removed old code
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 David Harvey <dmharvey@math.harvard.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cimport sage.rings.fast_arith
import sage.rings.fast_arith
cdef sage.rings.fast_arith.arith_int arith_int
arith_int  = sage.rings.fast_arith.arith_int()

ctypedef long long llong

import sage.rings.arith

from sage.libs.ntl import all as ntl
from sage.libs.ntl.ntl_ZZ_pX cimport ntl_ZZ_pX
from sage.rings.finite_rings.integer_mod_ring import Integers
from sage.rings.bernmm import bernmm_bern_modp



def verify_bernoulli_mod_p(data):
    """
    Computes checksum for bernoulli numbers.

    It checks the identity
        $$ \sum_{n=0}^{(p-3)/2} 2^{2n} (2n+1) B_{2n}  \equiv  -2 \pmod p $$

    (see "Irregular Primes to One Million", Buhler et al)

    INPUT:
        data -- list, same format as output of bernoulli_mod_p function

    OUTPUT:
        bool -- True if checksum passed

    EXAMPLES:
        sage: from sage.rings.bernoulli_mod_p import verify_bernoulli_mod_p
        sage: verify_bernoulli_mod_p(bernoulli_mod_p(next_prime(3)))
        True
        sage: verify_bernoulli_mod_p(bernoulli_mod_p(next_prime(1000)))
        True
        sage: verify_bernoulli_mod_p([1, 2, 4, 5, 4])
        True
        sage: verify_bernoulli_mod_p([1, 2, 3, 4, 5])
        False

    This one should test that long longs are working:
        sage: verify_bernoulli_mod_p(bernoulli_mod_p(next_prime(20000)))
        True

    AUTHOR: David Harvey
    """

    cdef int N, p, i, product, sum, value, element
    N = len(data)
    p = N*2 + 1
    product = 1
    sum = 0
    for i from 0 <= i < N:
        element = data[i]
        value = <int> (((((<llong> product) * (2*i+1)) % p) * element) % p)
        sum = (sum + value) % p
        product = (4 * product) % p

    if (sum + 2) % p == 0:
        return True
    else:
        return False


def bernoulli_mod_p(int p):
    r"""
    Returns the bernoulli numbers $B_0, B_2, ... B_{p-3}$ modulo $p$.

    INPUT:
        p -- integer, a prime

    OUTPUT:
        list -- Bernoulli numbers modulo $p$ as a list
                of integers [B(0), B(2), ... B(p-3)].

    ALGORITHM:
        Described in accompanying latex file.

    PERFORMANCE:
        Should be complexity $O(p \log p)$.

    EXAMPLES:
    Check the results against PARI's C-library implementation (that
    computes exact rationals) for $p = 37$:

        sage: bernoulli_mod_p(37)
         [1, 31, 16, 15, 16, 4, 17, 32, 22, 31, 15, 15, 17, 12, 29, 2, 0, 2]
        sage: [bernoulli(n) % 37 for n in xrange(0, 36, 2)]
         [1, 31, 16, 15, 16, 4, 17, 32, 22, 31, 15, 15, 17, 12, 29, 2, 0, 2]

    Boundary case:
        sage: bernoulli_mod_p(3)
         [1]

    AUTHOR:
        -- David Harvey (2006-08-06)

    """

    if p <= 2:
        raise ValueError, "p (=%s) must be a prime >= 3"%p

    if not sage.rings.arith.is_prime(p):
        raise ValueError, "p (=%s) must be a prime"%p

    cdef int g, gSqr, gInv, gInvSqr, isOdd

    g = sage.rings.arith.primitive_root(p)
    gInv = arith_int.c_inverse_mod_int(g, p)
    gSqr = ((<llong> g) * g) % p
    gInvSqr = ((<llong> gInv) * gInv) % p
    isOdd = ((p-1)/2) % 2

    # STEP 1: compute the polynomials G(X) and J(X)

    # These hold g^{i-1} and g^{-i} at the beginning of each iteration
    cdef llong gPower, gPowerInv
    gPower = gInv
    gPowerInv = 1

    # "constant" is (g-1)/2 mod p
    cdef int constant
    if g % 2:
        constant = (g-1)/2
    else:
        constant = (g+p-1)/2

    # fudge holds g^{i^2}, fudgeInv holds g^{-i^2}
    cdef int fudge, fudgeInv
    fudge = fudgeInv = 1

    cdef ntl_ZZ_pX G, J
    G = ntl.ZZ_pX([], ntl.ZZ(p))
    J = ntl.ZZ_pX([], ntl.ZZ(p))
    G.preallocate_space((p-1)/2)
    J.preallocate_space((p-1)/2)

    cdef int i
    cdef llong temp, h
    for i from 0 <= i < (p-1)/2:
        # compute h = h(g^i)/g^i  (h(x) is as in latex notes)
        temp = g * gPower
        h = ((p + constant - (temp / p)) * gPowerInv) % p
        gPower = temp % p
        gPowerInv = (gPowerInv * gInv) % p

        # store the coefficient  g^{i^2} h(g^i)/g^i
        G.setitem_from_int(i, <int> ((h * fudge) % p))

        # store the coefficient  g^{-i^2}
        J.setitem_from_int(i, fudgeInv)

        # update fudge and fudgeInv
        fudge = (((fudge * gPower) % p) * ((gPower * g) % p)) % p
        fudgeInv = (((fudgeInv * gPowerInv) % p) * ((gPowerInv * g) % p)) % p

    J.setitem_from_int(0, 0)

    # STEP 2: multiply the polynomials

    cdef ntl_ZZ_pX product
    product = G * J

    # STEP 3: assemble the result

    cdef int gSqrPower, value
    output = [1]
    gSqrPower = gSqr
    fudge = g
    for i from 1 <= i < (p-1)/2:
        value = product.getitem_as_int(i + (p-1)/2)

        if isOdd:
            value = (G.getitem_as_int(i) + product.getitem_as_int(i) - value + p) % p
        else:
            value = (G.getitem_as_int(i) + product.getitem_as_int(i) + value) % p

        value = (((4 * i * (<llong> fudge)) % p) * (<llong> value)) % p
        value = ((<llong> value) * (arith_int.c_inverse_mod_int(1 - gSqrPower, p))) % p

        output.append(value)

        gSqrPower = ((<llong> gSqrPower) * g) % p
        fudge = ((<llong> fudge) * gSqrPower) % p
        gSqrPower = ((<llong> gSqrPower) * g) % p

    return output



def bernoulli_mod_p_single(long p, long k):
    r"""
    Returns the bernoulli number $B_k$ mod $p$.

    If $B_k$ is not $p$-integral, an ArithmeticError is raised.

    INPUT:
        p -- integer, a prime
        k -- non-negative integer

    OUTPUT:
        The $k$-th bernoulli number mod $p$.

    EXAMPLES:
        sage: bernoulli_mod_p_single(1009, 48)
        628
        sage: bernoulli(48) % 1009
        628

        sage: bernoulli_mod_p_single(1, 5)
        Traceback (most recent call last):
        ...
        ValueError: p (=1) must be a prime >= 3

        sage: bernoulli_mod_p_single(100, 4)
        Traceback (most recent call last):
        ...
        ValueError: p (=100) must be a prime

        sage: bernoulli_mod_p_single(19, 5)
        0

        sage: bernoulli_mod_p_single(19, 18)
        Traceback (most recent call last):
        ...
        ArithmeticError: B_k is not integral at p

        sage: bernoulli_mod_p_single(19, -4)
        Traceback (most recent call last):
        ...
        ValueError: k must be non-negative

    Check results against bernoulli_mod_p:

        sage: bernoulli_mod_p(37)
         [1, 31, 16, 15, 16, 4, 17, 32, 22, 31, 15, 15, 17, 12, 29, 2, 0, 2]
        sage: [bernoulli_mod_p_single(37, n) % 37 for n in xrange(0, 36, 2)]
         [1, 31, 16, 15, 16, 4, 17, 32, 22, 31, 15, 15, 17, 12, 29, 2, 0, 2]

        sage: bernoulli_mod_p(31)
         [1, 26, 1, 17, 1, 9, 11, 27, 14, 23, 13, 22, 14, 8, 14]
        sage: [bernoulli_mod_p_single(31, n) % 31 for n in xrange(0, 30, 2)]
         [1, 26, 1, 17, 1, 9, 11, 27, 14, 23, 13, 22, 14, 8, 14]

        sage: bernoulli_mod_p(3)
         [1]
        sage: [bernoulli_mod_p_single(3, n) % 3 for n in xrange(0, 2, 2)]
         [1]

        sage: bernoulli_mod_p(5)
         [1, 1]
        sage: [bernoulli_mod_p_single(5, n) % 5 for n in xrange(0, 4, 2)]
         [1, 1]

        sage: bernoulli_mod_p(7)
         [1, 6, 3]
        sage: [bernoulli_mod_p_single(7, n) % 7 for n in xrange(0, 6, 2)]
         [1, 6, 3]

    AUTHOR:
        -- David Harvey (2007-08-31)
        -- David Harvey (2008-06): rewrote to use bernmm library

    """
    if p <= 2:
        raise ValueError, "p (=%s) must be a prime >= 3"%p

    if not sage.rings.arith.is_prime(p):
        raise ValueError, "p (=%s) must be a prime"%p

    R = Integers(p)

    cdef long x = bernmm_bern_modp(p, k)
    if x == -1:
        raise ArithmeticError, "B_k is not integral at p"
    return x


# ============ end of file
