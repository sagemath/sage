r"""
Fast calculation of cyclotomic polynomials

This module provides a function :func:`cyclotomic_coeffs`, which calculates the
coefficients of cyclotomic polynomials. This is not intended to be invoked
directly by the user, but it is called by the method
:meth:`~sage.rings.polynomial.polynomial_ring.PolynomialRing_general.cyclotomic_polynomial`
method of univariate polynomial ring objects and the top-level
:func:`~sage.misc.functional.cyclotomic_polynomial` function.
"""

#*****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sys

include "sage/ext/stdsage.pxi"
include "sage/ext/interrupt.pxi"
include "sage/ext/cdefs.pxi"

from sage.rings.arith import factor
from sage.rings.infinity import infinity
from sage.misc.misc import prod, subsets
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.libs.pari.gen cimport gen
from sage.libs.pari.all import pari

def cyclotomic_coeffs(nn, sparse=None):
    u"""
    This calculates the coefficients of the n-th cyclotomic polynomial
    by using the formula

    .. math::

        \\Phi_n(x) = \\prod_{d|n} (1-x^{n/d})^{\\mu(d)}

    where `\\mu(d)` is the Moebius function that is 1 if d has an even
    number of distinct prime divisors, -1 if it has an odd number of
    distinct prime divisors, and 0 if d is not squarefree.

    Multiplications and divisions by polynomials of the
    form `1-x^n` can be done very quickly in a single pass.

    If sparse is True, the result is returned as a dictionary of the non-zero
    entries, otherwise the result is returned as a list of python ints.

    EXAMPLES::

        sage: from sage.rings.polynomial.cyclotomic import cyclotomic_coeffs
        sage: cyclotomic_coeffs(30)
        [1, 1, 0, -1, -1, -1, 0, 1, 1]
        sage: cyclotomic_coeffs(10^5)
        {0: 1, 10000: -1, 40000: 1, 30000: -1, 20000: 1}
        sage: R = QQ['x']
        sage: R(cyclotomic_coeffs(30))
        x^8 + x^7 - x^5 - x^4 - x^3 + x + 1

    Check that it has the right degree::

        sage: euler_phi(30)
        8
        sage: R(cyclotomic_coeffs(14)).factor()
        x^6 - x^5 + x^4 - x^3 + x^2 - x + 1

    The coefficients are not always +/-1::

        sage: cyclotomic_coeffs(105)
        [1, 1, 1, 0, 0, -1, -1, -2, -1, -1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, -1, -1, -2, -1, -1, 0, 0, 1, 1, 1]

    In fact the height is not bounded by any polynomial in n (Erdos),
    although takes a while just to exceed linear::

        sage: v = cyclotomic_coeffs(1181895)
        sage: max(v)
        14102773

    The polynomial is a palindrome for any n::

        sage: n = ZZ.random_element(50000)
        sage: factor(n)
        3 * 10009
        sage: v = cyclotomic_coeffs(n, sparse=False)
        sage: v == list(reversed(v))
        True

    AUTHORS:

    - Robert Bradshaw (2007-10-27): initial version (inspired by work of Andrew
      Arnold and Michael Monagan)
    """
    factors = factor(nn)
    if any([e != 1 for p, e in factors]):
        # If there are primes that occur in the factorization with multiplicity
        # greater than one we use the fact that Phi_ar(x) = Phi_r(x^a) when all
        # primes dividing a divide r.
        rad = prod([p for p, e in factors])
        rad_coeffs = cyclotomic_coeffs(rad, sparse=True)
        pow = int(nn // rad)
        if sparse is None or sparse:
            L = {}
        else:
            L = [0] * (1 + pow * prod([p-1 for p, e in factors]))
        for mon, c in rad_coeffs.items():
            L[mon*pow] = c
        return L

    elif len(factors) == 1 and not sparse:
        # \Phi_p is easy to calculate for p prime.
        return [1] * factors[0][0]

    # The following bounds are from Michael Monagan:
    #    For all n < 169,828,113, the height of Phi_n(x) is less than 60 bits.
    #    At n = 169828113, we get a height of 31484567640915734951 which is 65 bits
    #    For n=10163195, the height of Phi_n(x) is 1376877780831,  40.32 bits.
    #    For n<10163195, the height of Phi_n(x) is <= 74989473, 26.16 bits.
    cdef long fits_long_limit = 169828113 if sizeof(long) >= 8 else 10163195
    if nn >= fits_long_limit and bateman_bound(nn) > sys.maxint:
        # Do this to avoid overflow.
        print "Warning: using PARI (slow!)"
        from sage.interfaces.gp import pari
        return [int(a) for a in pari.polcyclo(nn).Vecrev()]

    cdef long d, max_deg = 0, n = nn
    primes = [int(p) for p, e in factors]
    prime_subsets = list(subsets(primes))
    if n > 5000:
        prime_subsets.sort(my_cmp)

    for s in prime_subsets:
        if len(s) % 2 == 0:
            d = prod(s)
            max_deg += n / d

    if (<object>max_deg)*sizeof(long) > sys.maxint:
        raise MemoryError, "Not enough memory to calculate cyclotomic polynomial of %s" % n
    cdef long* coeffs = <long*>sage_malloc(sizeof(long) * (max_deg+1))
    if coeffs == NULL:
        raise MemoryError, "Not enough memory to calculate cyclotomic polynomial of %s" % n
    memset(coeffs, 0, sizeof(long) * (max_deg+1))
    coeffs[0] = 1

    cdef long k, dd, offset = 0, deg = 0
    for s in prime_subsets:
        if len(s) % 2 == 0:
            d = prod(s)
            dd = n / d
#            f *= (1-x^dd)
            sig_on()
            for k from deg+dd >= k >= dd:
                coeffs[k] -= coeffs[k-dd]
            deg += dd
            sig_off()

    prime_subsets.reverse()
    for s in prime_subsets:
        if len(s) % 2 == 1:
            d = prod(s)
            dd = n / d
#            f /= (1-x^dd)
            sig_on()
            for k from deg >= k > deg-dd:
                coeffs[k] = -coeffs[k]
            for k from deg-dd >= k >= offset:
                coeffs[k] = coeffs[k+dd] - coeffs[k]
            offset += dd
            sig_off()

    cdef long non_zero = 0
    if sparse is None:
        for k from offset <= k <= deg:
            non_zero += coeffs[k] != 0
        sparse = non_zero < 0.25*(deg-offset)

    if sparse:
        L = {}
        for k from offset <= k <= deg:
            if coeffs[k]:
                L[k-offset] = coeffs[k]
    else:
        L = [coeffs[k] for k from offset <= k <= deg]

    sage_free(coeffs)
    return L

def cyclotomic_value(n, x):
    """
    Returns the value of the `n`-th cyclotomic polynomial evaulated at `x`.

    INPUT:

    - n -- an Integer, specifying which cyclotomic polynomial is to be
      evaluated.
    - x -- an element of a ring.

    OUTPUT:

    - the value of the cyclotomic polynomial `\Phi_n` at `x`.

    ALGORITHM:

    - Reduce to the case that n is squarefree: use the identity

    .. MATH::

        \Phi_n(x) = \Phi_q(x^{n/q})

    where `q` is the radical of `n`.

    - Use the identity

    .. MATH::

        \Phi_n(x) = \prod_{d | n} (x^d - 1)^{\mu(n / d)},

    where `\mu` is the Moebius function.

    - Handles the case that x^d = 1 for some d, but not the case that
      x^d - 1 is non-invertible: in this case polynomial evaluation is
      used instead.

    EXAMPLES::

        sage: cyclotomic_value(51, 3)
        1282860140677441
        sage: cyclotomic_polynomial(51)(3)
        1282860140677441

    It works for non-integral values as well::

        sage: cyclotomic_value(144, 4/3)
        79148745433504023621920372161/79766443076872509863361
        sage: cyclotomic_polynomial(144)(4/3)
        79148745433504023621920372161/79766443076872509863361

    TESTS::

        sage: K.<i> = NumberField(polygen(QQ)^2 + 1)
        sage: R.<x> = QQ[]
        sage: for y in [-1, 0, 1, 2, 1/2, mod(3, 8), GF(9,'a').gen(), Zp(3)(54), i, x^2+2]:
        ...       for n in range(1, 61):
        ...           val1 = cyclotomic_value(n, y)
        ...           val2 = cyclotomic_polynomial(n)(y)
        ...           assert val1 == val2 and val1.parent() is val2.parent()

        sage: cyclotomic_value(20, I)
        5
        sage: a = cyclotomic_value(10, mod(3, 11)); a
        6
        sage: a.parent()
        Ring of integers modulo 11
        sage: cyclotomic_value(30, -1.0)
        1.00000000000000
        sage: S.<t> = R.quotient(R.cyclotomic_polynomial(15))
        sage: cyclotomic_value(15, t)
        0
        sage: cyclotomic_value(30, t)
        2*t^7 - 2*t^5 - 2*t^3 + 2*t
        sage: S.<t> = R.quotient(x^10)
        sage: cyclotomic_value(2^128-1, t)
        -t^7 - t^6 - t^5 + t^2 + t + 1
        sage: cyclotomic_value(10,mod(3,4))
        1
    """
    n = int(n)
    if n == 1:
        return x - 1
    if n <= 0:
        raise ValueError, "n must be positive"
    try:
        return x.parent()(pari.polcyclo_eval(n, x._pari_()))
    except StandardError:
        pass
    # The following is modeled on the implementation in Pari
    factors = factor(n)
    cdef Py_ssize_t i, j, ti, L, root_of_unity = -1
    primes = [p for p, e in factors]
    L = len(primes)
    if any([e != 1 for p, e in factors]):
        # If there are primes that occur in the factorization with multiplicity
        # greater than one we use the fact that Phi_ar(x) = Phi_r(x^a) when all
        # primes dividing a divide r.
        rad = prod(primes)
        pow = n // rad
        x = x**pow
        n = rad
    if x == 1:
        # if n is prime, return n
        if L == 1:
            return n * x # in case the parent of x has nonzero characteristic
        else:
            return x
    xd = [x] # the x^d for d | n
    cdef char mu
    cdef char* md = <char*>sage_malloc(sizeof(char) * (1 << L)) # the mu(d) for d | n
    one = x.parent()(1)
    try:
        md[0] = 1
        if L & 1:
            mu = -1
            num = 1
            den = x - 1
        else:
            mu = 1
            num = x - 1
            den = 1
        for i in range(L):
            ti = 1 << i
            p = primes[i]
            for j in range(ti):
                xpow = xd[j]**p
                xd.append(xpow)
                md[ti+j] = -md[j]
                # if xpow = 1, we record such smallest index,
                # and deal with the corresponding factors at the end.
                if xpow == one:
                    if root_of_unity == -1:
                        root_of_unity = ti+j
                elif mu == md[ti+j]:
                    num *= xpow - one
                else:
                    den *= xpow - one
    finally:
        sage_free(md)
    try:
        ans = num / den
    except ZeroDivisionError:
        # We fall back on evaluation of the cyclotomic polynomial.
        # This case is triggered in cyclotomic_value(10, mod(3, 4)) for example.
        from sage.misc.functional import cyclotomic_polynomial
        return cyclotomic_polynomial(n)(x)
    if root_of_unity >= 0:
        # x is a root of unity.  If root_of_unity=2^L, x is a primitive
        # root of unity and the value is zero
        if root_of_unity == (1 << L) - 1:
            return x - x # preserves the parent, as well as precision for p-adic x
        # x is a primitive d-th root of unity, where d|n and d<n.
        # If root_of_unity = (1<<L) - (1<<(i-1)) - 1 for some i < L,
        # then n/d == primes[i] and we need to multiply by primes[i],
        # otherwise n/d is composite and nothing more needs to be done.
        for i in range(L):
            if root_of_unity + (1 << i) + 1 == 1 << L:
                ans *= primes[i]
                break
    return x.parent()(ans)

def bateman_bound(nn):
    _, n = nn.val_unit(2)
    primes = [p for p, _ in factor(n)]
    j = len(primes)
    return prod([primes[k]^(2^(j-k-2)-1) for k in range(j-2)])

def my_cmp(a, b):
    return int(prod(b) - prod(a))
