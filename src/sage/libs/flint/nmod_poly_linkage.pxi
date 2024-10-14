r"""
Linkage for arithmetic with FLINT's nmod_poly_t elements.

This file provides the backend for \class{Polynomial_zmod_flint} via
templating.

AUTHOR:
    -- Martin Albrecht (2009-01) another initial implementation
    -- Burcin Erocal (2008-11) initial implementation
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Burcin Erocal <burcin@erocal.org>
#       Copyright (C) 2009 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.signals cimport sig_on, sig_off
from cysignals.memory cimport sig_malloc, sig_free

from sage.libs.flint.types cimport *
from sage.libs.flint.nmod_poly cimport *
from sage.libs.flint.nmod_poly_factor cimport *
from sage.libs.flint.ulong_extras cimport *
from sage.structure.factorization import Factorization

cdef inline celement *celement_new(unsigned long n) noexcept:
    cdef celement *g = <celement *>sig_malloc(sizeof(nmod_poly_t))
    nmod_poly_init(g, n)
    return g

cdef inline int celement_delete(nmod_poly_t e, unsigned long n) noexcept:
    nmod_poly_clear(e)
    sig_free(e)

cdef inline int celement_construct(nmod_poly_t e, unsigned long n) noexcept:
    """
    EXAMPLES::

        sage: P.<x> = GF(32003)[]

        sage: Q.<x> = GF(7)[]
    """
    nmod_poly_init(e, n)

cdef inline int celement_destruct(nmod_poly_t e, unsigned long n) noexcept:
    """
    EXAMPLES::

        sage: P.<x> = GF(32003)[]
        sage: del x

        sage: Q.<x> = GF(7)[]
        sage: del x
    """
    nmod_poly_clear(e)

cdef inline int celement_gen(nmod_poly_t e, long i, unsigned long n) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(32003)[]

        sage: Q.<x> = GF(7)[]
    """
    nmod_poly_zero(e)
    nmod_poly_set_coeff_ui(e, 1, 1)

cdef object celement_repr(nmod_poly_t e, unsigned long n):
    raise NotImplementedError

cdef inline int celement_set(nmod_poly_t res, nmod_poly_t a, unsigned long n) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(32003)[]
        sage: y = copy(x)
        sage: y is x
        False
        sage: y == x
        True

        sage: Q.<x> = GF(7)[]
        sage: y = copy(x)
        sage: y is x
        False
        sage: y == x
        True

        sage: R.<x> = PolynomialRing(Integers(121))
        sage: S.<y> = PolynomialRing(Integers(11))
        sage: S(50*x)
        6*y
        sage: R(S(50*x))
        6*x
    """
    cdef unsigned long i
    if a.mod.n <= n:
        nmod_poly_set(res, a)
    else:
        nmod_poly_zero(res)
        for i from 0 <= i < a.length:
            nmod_poly_set_coeff_ui(res, i, nmod_poly_get_coeff_ui(a, i) % n)

cdef inline int celement_set_si(nmod_poly_t res, long i, unsigned long n) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(32003)[]
        sage: P(32003)
        0
        sage: P(1)
        1
        sage: P(32004)
        1

        sage: Q.<x> = GF(7)[]
        sage: Q(7)
        0
        sage: Q(1)
        1
        sage: Q(8)
        1
    """
    while i < 0:
        i += n
    nmod_poly_zero(res)
    if i:
        nmod_poly_set_coeff_ui(res, 0, <unsigned long>i)

cdef inline long celement_get_si(nmod_poly_t res, unsigned long n) except -2:
    raise NotImplementedError

cdef inline bint celement_is_zero(nmod_poly_t a, unsigned long n) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(32003)[]
        sage: P(1).is_zero()
        False
        sage: P(0).is_zero()
        True

        sage: Q.<x> = GF(7)[]
        sage: Q(1).is_zero()
        False
        sage: Q(0).is_zero()
        True
    """
    return nmod_poly_is_zero(a)

cdef inline bint celement_is_one(nmod_poly_t a, unsigned long n) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(32003)[]
        sage: P(1).is_one()
        True
        sage: P(0).is_one()
        False

        sage: Q.<x> = GF(7)[]
        sage: Q(1).is_one()
        True
        sage: Q(0).is_one()
        False
    """

    return nmod_poly_is_one(a)

cdef inline bint celement_equal(nmod_poly_t a, nmod_poly_t b, unsigned long n) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(32003)[]
        sage: (3*2)*x == 3*(2*x)
        True
        sage: (3*2)*x + 2 == 3*(2*x) + 1 + 1
        True
        sage: (3*2)*x + 7 == 3*(2*x) + 1 + 1
        False

        sage: Q.<x> = GF(7)[]
        sage: (3*2)*x == 3*(2*x)
        True
        sage: (3*2)*x + 2 == 3*(2*x) + 1 + 1
        True
        sage: (3*2)*x + 7 == 3*(2*x) + 1 + 1
        False
    """
    return nmod_poly_equal(a, b)

cdef inline int celement_cmp(nmod_poly_t l, nmod_poly_t r, unsigned long n) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(32003)[]
        sage: x > x
        False
        sage: x^2 > x
        True
        sage: 3*x > x
        True

        sage: Q.<x> = GF(7)[]
        sage: x > x
        False
        sage: x^2 > x
        True
        sage: 3*x > x
        True

        sage: f = x^64 + x^20 + 1
        sage: g = x^63 + x^20 + 1
        sage: f > g
        True

        sage: f = x^64 + x^10 + 1
        sage: g = x^64 + x^20 + 1
        sage: f < g
        True

        sage: f = x^64 + x^20
        sage: g = x^64 + x^20 + 1
        sage: f < g
        True
    """
    cdef int deg_right = nmod_poly_degree(r)
    cdef int degdiff = deg_right - nmod_poly_degree(l)
    cdef int i
    cdef unsigned long rcoeff, lcoeff
    if degdiff > 0:
        return -1
    elif degdiff < 0:
        return 1
    else:
        if nmod_poly_equal(l, r):
            return 0
        i = deg_right
        while i >= 0:
            rcoeff = nmod_poly_get_coeff_ui(r, i)
            lcoeff = nmod_poly_get_coeff_ui(l, i)
            if lcoeff < rcoeff:
                return -1
            if lcoeff > rcoeff:
                return 1
            i -= 1
        return 0

cdef long celement_len(nmod_poly_t a, unsigned long n) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(32003)[]
        sage: (x + 1).degree()
        1
        sage: (x).degree()
        1
        sage: P(0).degree()
        -1

        sage: Q.<x> = GF(7)[]
        sage: (x + 1).degree()
        1
        sage: (x).degree()
        1
        sage: P(0).degree()
        -1
    """
    return <long>nmod_poly_length(a)

cdef inline int celement_add(nmod_poly_t res, nmod_poly_t a, nmod_poly_t b, unsigned long n) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(32003)[]
        sage: x + 1
        x + 1

        sage: Q.<x> = GF(7)[]
        sage: x + 1
        x + 1
    """
    nmod_poly_add(res, a, b)

cdef inline int celement_sub(nmod_poly_t res, nmod_poly_t a, nmod_poly_t b, unsigned long n) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(32003)[]
        sage: x - 1
        x + 32002

        sage: Q.<x> = GF(7)[]
        sage: x - 1
        x + 6
    """
    nmod_poly_sub(res, a, b)

cdef inline int celement_neg(nmod_poly_t res, nmod_poly_t a, unsigned long n) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(32003)[]
        sage: -(x + 2)
        32002*x + 32001

        sage: Q.<x> = GF(7)[]
        sage: -(x + 2)
        6*x + 5
    """
    nmod_poly_neg(res, a)

cdef inline int celement_mul_scalar(nmod_poly_t res, nmod_poly_t p,
        object c, unsigned long n) except -2:
    """
    TESTS::

        sage: P.<x> = GF(32003)[]
        sage: p = P.random_element(degree=2)
        sage: (389*p).coefficients() == [389*x for x in p.coefficients()]
        True
        sage: p = P.random_element(degree=8)
        sage: (p*9836).coefficients() == [x*9836 for x in p.coefficients()]
        True
    """
    nmod_poly_scalar_mul_nmod(res, p, (<unsigned long>c)%n)

cdef inline int celement_mul(nmod_poly_t res, nmod_poly_t a, nmod_poly_t b, unsigned long n) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(32003)[]
        sage: (x + 1) * (x + 2)
        x^2 + 3*x + 2

        sage: Q.<x> = GF(7)[]
        sage: (x + 1) * (x + 2)
        x^2 + 3*x + 2
    """
    nmod_poly_mul(res, a, b)

cdef inline int celement_div(nmod_poly_t res, nmod_poly_t a, nmod_poly_t b, unsigned long n) except -2:
    raise NotImplementedError

cdef inline int celement_truncate(nmod_poly_t res, nmod_poly_t a, long len, unsigned long n) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(7)[]
        sage: p = 4*x^4 + 3*x^3 + 2*x^2 + x
        sage: p.truncate(3)
        2*x^2 + x

        sage: Q.<x> = GF(32003)[]
        sage: q = 1 + x + x^2 * Q.random_element()
        sage: q.truncate(2)
        x + 1
    """
    nmod_poly_set(res, a)
    nmod_poly_truncate(res, len)

cdef inline int celement_floordiv(nmod_poly_t res, nmod_poly_t a, nmod_poly_t b, unsigned long n) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(32003)[]
        sage: (x + 1) // (x + 2)
        1
        sage: (3*x^2 + 1) // (x + 2)
        3*x + 31997
        sage: (x^2 + 3*x + 2)//(x + 1)
        x + 2
        sage: (x^2 + 3*x + 2)//(x + 2)
        x + 1

        sage: Q.<x> = GF(7)[]
        sage: (x + 1) // (x + 2)
        1
        sage: (3*x^2 + 1) // (x + 2)
        3*x + 1
        sage: (x^2 + 3*x + 2)//(x + 1)
        x + 2
        sage: (x^2 + 3*x + 2)//(x + 2)
        x + 1
    """
    nmod_poly_div(res, a, b)

cdef inline int celement_mod(nmod_poly_t res, nmod_poly_t a, nmod_poly_t b, unsigned long n) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(32003)[]
        sage: f = 24998*x^2 + 29761*x + 2252
        sage: g = 20778*x^2 + 15346*x + 12697

        sage: f % g
        5815*x + 10280

        sage: f^5 % g
        7231*x + 17274

        sage: R.<x> = Integers(81)[]
        sage: f = x^7 + x + 1; g = x^3
        sage: r = f % g; r
        x + 1
        sage: g * x^4 + r
        x^7 + x + 1

        sage: f = x^3 + 1
        sage: f % 3
        Traceback (most recent call last):
        ...
        ValueError: Leading coefficient of a must be invertible.

        sage: f % 0
        Traceback (most recent call last):
        ...
        ZeroDivisionError
    """
    cdef nmod_poly_t q
    cdef unsigned long leadcoeff, modulus

    nmod_poly_init(q, n)
    leadcoeff = nmod_poly_get_coeff_ui(b, nmod_poly_degree(b))
    modulus = nmod_poly_modulus(b)
    if (leadcoeff > 1 and n_gcd(modulus,leadcoeff) != 1):
        raise ValueError("Leading coefficient of a must be invertible.")

    nmod_poly_divrem(q, res, a, b)
    nmod_poly_clear(q)

cdef inline int celement_quorem(nmod_poly_t q, nmod_poly_t r, nmod_poly_t a, nmod_poly_t b, unsigned long n) except -2:
    """
    EXAMPLES::

        sage: R.<x> = Integers(125)[]
        sage: f = x^5+1; g = (x+1)^2
        sage: q, r = f.quo_rem(g)
        sage: q
        x^3 + 123*x^2 + 3*x + 121
        sage: r
        5*x + 5
        sage: q*g + r
        x^5 + 1

        sage: x.quo_rem(5*x)
        Traceback (most recent call last):
        ...
        ValueError: Leading coefficient of a must be invertible.

        sage: x.quo_rem(0)
        Traceback (most recent call last):
        ...
        ZeroDivisionError
    """
    cdef unsigned long leadcoeff, modulus

    leadcoeff = nmod_poly_get_coeff_ui(b, nmod_poly_degree(b))
    modulus = nmod_poly_modulus(b)
    if (leadcoeff > 1 and n_gcd(modulus,leadcoeff) != 1):
        raise ValueError("Leading coefficient of a must be invertible.")

    nmod_poly_divrem(q, r, a, b)

cdef inline int celement_inv(nmod_poly_t res, nmod_poly_t a, unsigned long n) except -2:
    raise NotImplementedError

cdef inline int celement_pow(nmod_poly_t res, nmod_poly_t x, long e, nmod_poly_t modulus, unsigned long n) except -2:
    """
    Compute `x^e`, possibly modulo ``modulus``.

    INPUT:

    - ``x`` -- polynomial; the base

    - ``e`` -- integer; the exponent

    - ``modulus`` -- polynomial or NULL; if not NULL, then perform a modular exponentiation

    - ``n`` -- integer; not used, but all polynomials' coefficients are understood modulo ``n``

    EXAMPLES::

        sage: P.<x> = GF(32003)[]
        sage: f = 24998*x^2 + 29761*x + 2252

        sage: f*f
        9426*x^4 + 15477*x^3 + 6531*x^2 + 14980*x + 15030
        sage: f^2
        9426*x^4 + 15477*x^3 + 6531*x^2 + 14980*x + 15030

        sage: f*f*f
        25062*x^6 + 30670*x^5 + 15816*x^4 + 20746*x^3 + 9142*x^2 + 5697*x + 20389
        sage: f^3
        25062*x^6 + 30670*x^5 + 15816*x^4 + 20746*x^3 + 9142*x^2 + 5697*x + 20389

        sage: f*f*f*f*f
        20269*x^10 + 20535*x^9 + 7313*x^8 + 7311*x^7 + 16853*x^6 + 142*x^5 + 23853*x^4 + 12065*x^3 + 516*x^2 + 8473*x + 17945
        sage: f^5
        20269*x^10 + 20535*x^9 + 7313*x^8 + 7311*x^7 + 16853*x^6 + 142*x^5 + 23853*x^4 + 12065*x^3 + 516*x^2 + 8473*x + 17945

        sage: f^0
        1

        sage: f^1
        24998*x^2 + 29761*x + 2252

        sage: f^-1
        18649/(x^2 + 16863*x + 9612)

        sage: f^-5
        24620/(x^10 + 20309*x^9 + 29185*x^8 + 11948*x^7 + 1965*x^6 + 7713*x^5 + 5810*x^4 + 20457*x^3 + 30732*x^2 + 9706*x + 4485)

     Testing the modulus::

        sage: g = 20778*x^2 + 15346*x + 12697

        sage: pow(f, 2, g)
        15328*x + 6968
        sage: f^2 % g
        15328*x + 6968

        sage: pow(f, -2, g)
        16346/(x + 251)
        sage: (f^2 % g)^-1
        16346/(x + 251)

        sage: pow(f, 5, g)
        7231*x + 17274
        sage: f^5 % g
        7231*x + 17274

    Make sure that exponentiation can be interrupted, see :issue:`17470`::

        sage: n = 2^23
        sage: alarm(0.2); x^n; cancel_alarm()
        Traceback (most recent call last):
        ...
        AlarmInterrupt
    """
    if modulus != NULL:
        sig_on()
        nmod_poly_powmod_ui_binexp(res, x, e, modulus)
        sig_off()
    else:
        sig_on()
        nmod_poly_pow(res, x, e)
        sig_off()

cdef inline int celement_gcd(nmod_poly_t res, nmod_poly_t a, nmod_poly_t b, unsigned long n) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(32003)[]
        sage: f = P.random_element(degree=4)
        sage: g = P.random_element(degree=3)
        sage: h = P.random_element(degree=2)
        sage: F = f*g
        sage: G = f*h
        sage: d = (F).gcd(G)
        sage: (F//d)*d == F
        True
        sage: (G//d)*d == G
        True

        sage: Q.<x> = GF(7)[]
        sage: f = Q.random_element(degree=4)
        sage: g = Q.random_element(degree=3)
        sage: h = Q.random_element(degree=2)
        sage: F = f*g
        sage: G = f*h
        sage: d = (F).gcd(G)
        sage: (F//d)*d == F
        True
        sage: (G//d)*d == G
        True

    Check that we catch a case where FLINT fails.  These generate the unit
    ideal, so their GCD should be 1 (or another unit)::

        sage: R.<x> = Integers(121)[]
        sage: f = 11*x^2 + 1
        sage: f - 61*x*f.derivative()
        1
        sage: gcd(f, f.derivative())
        Traceback (most recent call last):
        ...
        RuntimeError: FLINT gcd calculation failed
    """
    if celement_is_zero(b, n):
        nmod_poly_set(res, a)
        return 0

    # FLINT provides no interface for detecting errors here
    try:
        sig_on()
        try:
            nmod_poly_gcd(res, a, b)
        finally:
            sig_off()
    except RuntimeError:
        raise RuntimeError("FLINT gcd calculation failed")

    cdef unsigned long leadcoeff = nmod_poly_get_coeff_ui(res, nmod_poly_degree(res))
    cdef unsigned long modulus = nmod_poly_modulus(res)
    if n_gcd(modulus, leadcoeff) == 1:
        nmod_poly_make_monic(res, res)

cdef inline int celement_xgcd(nmod_poly_t res, nmod_poly_t s, nmod_poly_t t, nmod_poly_t a, nmod_poly_t b, unsigned long n) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(32003)[]
        sage: f = P.random_element(degree=4)
        sage: g = P.random_element(degree=3)
        sage: h = P.random_element(degree=2)
        sage: F = f*g
        sage: G = f*h
        sage: d,s,t = (F).xgcd(G)
        sage: (F//d)*d == F
        True
        sage: (G//d)*d == G
        True

        sage: Q.<x> = GF(7)[]
        sage: f = Q.random_element(degree=4)
        sage: g = Q.random_element(degree=3)
        sage: h = Q.random_element(degree=2)
        sage: F = f*g
        sage: G = f*h
        sage: d,s,t = (F).xgcd(G)
        sage: (F//d)*d == F
        True
        sage: (G//d)*d == G
        True

    TESTS:

    Ensure that :issue:`38537` is fixed::

        sage: k = Zmod(2**16)
        sage: R.<x> = k[]
        sage: u = x + 10161
        sage: v = x + 10681
        sage: u.xgcd(v)
        Traceback (most recent call last):
        ...
        ValueError: non-invertible elements encountered during XGCD
    """
    try:
        sig_on()
        nmod_poly_xgcd(res, s, t, a, b)
        sig_off()
    except RuntimeError:
        raise ValueError("non-invertible elements encountered during XGCD")


cdef factor_helper(Polynomial_zmod_flint poly, bint squarefree=False):
    """
    EXAMPLES::

        sage: P.<x> = GF(1009)[]
        sage: factors = (prod(P.random_element(degree=2) for i in range(5))).factor()
        sage: all(factor.is_irreducible() for factor, mult in factors)
        True

        sage: def nonzero_random(P):
        ....:     r = P.random_element()
        ....:     while r == 0:
        ....:         r = P.random_element()
        ....:     return r

        sage: p = (prod(nonzero_random(P)^i for i in range(5)))
        sage: decomp = p.squarefree_decomposition()
        sage: any(factor.is_square() for factor, mult in decomp)
        False
        sage: start = 0
        sage: for factor, mult in decomp:
        ....:     assert mult > start
        ....:     start = mult
    """
    cdef nmod_poly_factor_t factors_c
    nmod_poly_factor_init(factors_c)

    sig_on()
    if squarefree:
        nmod_poly_factor_squarefree(factors_c, &poly.x)
    else:
        nmod_poly_factor(factors_c, &poly.x)
    sig_off()

    factor_list = []
    cdef Polynomial_zmod_flint t
    for i in range(factors_c.num):
        t = poly._new()
        nmod_poly_swap(&t.x, &factors_c.p[i])
        factor_list.append((t, factors_c.exp[i]))

    nmod_poly_factor_clear(factors_c)

    return Factorization(factor_list, unit=poly.leading_coefficient(),
            sort=(not squarefree))
