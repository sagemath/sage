r"""
Linkage for arithmetic with NTL's ZZ_pX elements.

This file provides the backend for \class{Polynomial_ZZ_pX} via
templating.
"""

from cysignals.signals cimport sig_on, sig_off

from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class
from sage.libs.ntl.ZZ_pX cimport *
from sage.libs.ntl.ntl_ZZ_p cimport ntl_ZZ_p
from sage.libs.ntl.types cimport ZZ_pX_c, ZZ_pX_c

cdef ZZ_pX_c *celement_new(cparent parent):
    """
    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**80)),implementation='NTL')
    """
    if parent != NULL:
        parent[0].restore()
    return new ZZ_pX_c()

cdef int celement_delete(ZZ_pX_c *e, cparent parent):
    """
    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**80)),implementation='NTL')
        sage: del x
    """
    if parent != NULL:
        parent[0].restore()
    del e

cdef int celement_construct(ZZ_pX_c *e, cparent parent):
    """
    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**80)),implementation='NTL')
    """
    if parent != NULL:
        parent[0].restore()

cdef int celement_destruct(ZZ_pX_c *e, cparent parent):
    """
    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**80)),implementation='NTL')
        sage: del x
    """
    # do not call restore here
    # 1) the NTL context might have already been destroyed when exiting Python
    # 2) you better not make any NTL calls after destruct, no need to set the context

cdef int celement_gen(ZZ_pX_c *e, long i, cparent parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**80)),implementation='NTL')
    """
    if parent != NULL:
        parent[0].restore()
    ZZ_pX_SetX(e[0])

cdef object celement_repr(ZZ_pX_c *e, cparent parent):
    """
    We ignore NTL's printing.

    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**80)),implementation='NTL')
        sage: x
        x
    """
    raise NotImplementedError

cdef inline int celement_set(ZZ_pX_c* res, ZZ_pX_c* a, cparent parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**80)),implementation='NTL')
        sage: y = x
        sage: y
        x
    """
    res[0] = a[0]

cdef inline int celement_set_si(ZZ_pX_c* res, long i, cparent parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**80)),implementation='NTL')
        sage: P(0)
        0
        sage: P(17)
        17
        sage: P(next_prime(2**80))
        0
    """
    if parent != NULL:
        parent[0].restore()
    ZZ_pX_SetCoeff_long(res[0], 0, i)

cdef inline long celement_get_si(ZZ_pX_c* res, cparent parent) except -2:
    raise NotImplementedError

cdef inline bint celement_is_zero(ZZ_pX_c* a, cparent parent) except -2:
    """
    TESTS::

        sage: R.<x> = PolynomialRing(Integers(12^29), implementation='NTL')
        sage: f = x^4 + 1
        sage: not f
        False
        sage: not (x-x)
        True

        sage: P.<x> = PolynomialRing(GF(next_prime(2**80)),implementation='NTL')
        sage: bool(x), x.is_zero()
        (True, False)
        sage: bool(P(0)), P(0).is_zero()
        (False, True)
    """
    return ZZ_pX_IsZero(a[0])

cdef inline bint celement_is_one(ZZ_pX_c *a, cparent parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**80)),implementation='NTL')
        sage: x.is_one()
        False
        sage: P(1).is_one()
        True
    """
    if parent != NULL:
        parent[0].restore()
    return ZZ_pX_IsOne(a[0])

cdef inline bint celement_equal(ZZ_pX_c *a, ZZ_pX_c *b, cparent parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
        sage: x == x
        True
        sage: y = x; x == y
        True
        sage: x^2 + 1 == x^2 + x
        False
    """
    if parent != NULL:
        parent[0].restore()
    return a[0] == b[0]

cdef inline int celement_cmp(ZZ_pX_c *a, ZZ_pX_c *b, cparent parent) except -2:
    """
    Not used.

    Comparison is implemented in
    ``sage/rings/polynomial/polynomial_zz_pex.pyx`` instead.
    """
    raise NotImplementedError

cdef long celement_len(ZZ_pX_c *a, cparent parent) except -2:
    """
    EXAMPLES::

        sage: R.<x> = PolynomialRing(Integers(14^34), implementation='NTL')
        sage: f = x^4 - x - 1
        sage: f.degree()
        4
        sage: f = 14^43*x + 1
        sage: f.degree()
        0
        sage: P.<x> = PolynomialRing(GF(next_prime(2**80)),implementation='NTL')
        sage: x.degree()
        1
        sage: (x+1).degree()
        1
    """
    if parent != NULL:
        parent[0].restore()
    return int(ZZ_pX_deg(a[0]))+1

cdef inline int celement_add(ZZ_pX_c *res, ZZ_pX_c *a, ZZ_pX_c *b, cparent parent) except -2:
    """
    TESTS::

        sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
        sage: (x+5) + (x^2 - 6)
        x^2 + x + 999999999999999999999999999999
    """
    if parent != NULL:
        parent[0].restore()
    ZZ_pX_add(res[0], a[0], b[0])

cdef inline int celement_sub(ZZ_pX_c* res, ZZ_pX_c* a, ZZ_pX_c* b, cparent parent) except -2:
    """
    TESTS::

        sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
        sage: (x+5) - (x^2 - 6)
        999999999999999999999999999999*x^2 + x + 11
    """
    if parent != NULL:
        parent[0].restore()
    ZZ_pX_sub(res[0], a[0], b[0])

cdef inline int celement_neg(ZZ_pX_c* res, ZZ_pX_c* a, cparent parent) except -2:
    """
    EXAMPLES::

        sage: K.<a> = GF(next_prime(2**60)**3)
        sage: P.<x> = PolynomialRing(K,implementation='NTL')
        sage: -x
        1152921504606847008*x
    """
    if parent != NULL:
        parent[0].restore()
    ZZ_pX_negate(res[0], a[0])

cdef inline int celement_mul_scalar(ZZ_pX_c* res, ZZ_pX_c* p, object c, cparent parent) except -1:
    raise NotImplementedError

cdef inline int celement_mul(ZZ_pX_c* res, ZZ_pX_c* a, ZZ_pX_c* b, cparent parent) except -2:
    """
    TESTS::

        sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
        sage: (x+5) * (x^2 - 1)
        x^3 + 5*x^2 + 999999999999999999999999999999*x + 999999999999999999999999999995
    """
    if parent != NULL:
        parent[0].restore()
    ZZ_pX_mul(res[0], a[0], b[0])

cdef inline int celement_truncate(ZZ_pX_c* res, ZZ_pX_c* a, long len, cparent parent) except -2:
    """
    Returns this polynomial mod `x^n`.

    EXAMPLES::

        sage: R.<x> = PolynomialRing(Integers(15^30), implementation='NTL')
        sage: f = sum(x^n for n in range(10)); f
        x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
        sage: f.truncate(6)
        x^5 + x^4 + x^3 + x^2 + x + 1
    """
    ZZ_pX_trunc(res[0], a[0], len)

cdef inline int celement_div(ZZ_pX_c* res, ZZ_pX_c* a, ZZ_pX_c* b, cparent parent) except -2:
    if parent != NULL:
        parent[0].restore()
    return ZZ_pX_divide(res[0], a[0], b[0])

cdef inline int celement_floordiv(ZZ_pX_c* res, ZZ_pX_c* a, ZZ_pX_c* b, cparent parent) except -2:
    """
    Returns the whole part of self/right, without remainder.

    For q = n // d, we have deg(n - q*d) < deg(d)

    EXAMPLES::

        sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
        sage: f = x^7 + 1; g = x^2 - 1
        sage: q = f // g; q
        x^5 + x^3 + x
        sage: f - q*g
        x + 1
    """
    if parent != NULL:
        parent[0].restore()
    ZZ_pX_div(res[0], a[0], b[0])

cdef inline int celement_mod(ZZ_pX_c* res, ZZ_pX_c* a, ZZ_pX_c* b, cparent parent) except -2:
    """
    EXAMPLES::

        sage: R.<x> = PolynomialRing(Integers(9^30), implementation='NTL')
        sage: f = x^7 + x + 1; g = x^3 - 1
        sage: r = f % g; r
        2*x + 1
        sage: g * (x^4 + x) + r
        x^7 + x + 1
    """
    if parent != NULL:
        parent[0].restore()
    ZZ_pX_rem(res[0], a[0], b[0])

cdef inline int celement_quorem(ZZ_pX_c* q, ZZ_pX_c* r, ZZ_pX_c* a, ZZ_pX_c* b, cparent parent) except -2:
    """
    Returns `q` and `r`, with the degree of `r` less than the degree of `right`,
    such that `q * right + r = self`.

    EXAMPLES::

        sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
        sage: f = x^5+1; g = (x+1)^2
        sage: q, r = f.quo_rem(g)
        sage: q
        x^3 + 999999999999999999999999999998*x^2 + 3*x + 999999999999999999999999999996
        sage: r
        5*x + 5
        sage: q*g + r
        x^5 + 1
    """
    if parent != NULL:
        parent[0].restore()
    ZZ_pX_DivRem(q[0], r[0], a[0], b[0])

cdef inline int celement_inv(ZZ_pX_c* res, ZZ_pX_c* a, cparent parent) except -2:
    raise NotImplementedError

cdef inline int celement_pow(ZZ_pX_c* res, ZZ_pX_c* x, long e, ZZ_pX_c *modulus, cparent parent) except -2:
    """
    TESTS::

        sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
        sage: (x+1)^5
        x^5 + 5*x^4 + 10*x^3 + 10*x^2 + 5*x + 1

    We define ``0^0`` to be unity, :trac:`13895`::

        sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
        sage: R(0)^0
        1

    The value returned from ``0^0`` should belong to our ring::

        sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
        sage: type(R(0)^0) == type(R(0))
        True

    """
    if parent != NULL:
        parent[0].restore()

    cdef ZZ_pX_Modulus_c mod
    cdef ZZ_pX_c y
    if modulus == NULL:
        if ZZ_pX_IsX(x[0]):
            sig_on()
            ZZ_pX_LeftShift(res[0], x[0], e - 1)
            sig_off()
        else:
            sig_on()
            ZZ_pX_power(res[0], x[0], e)
            sig_off()
    else:
        if ZZ_pX_deg(modulus[0]) == 1:
             ZZ_pX_rem(y, x[0], modulus[0])
             sig_on()
             ZZ_pX_power(res[0], y, e)
             sig_off()
             return 0
        ZZ_pX_Modulus_build(mod, modulus[0])
        if ZZ_pX_IsX(x[0]):
            sig_on()
            ZZ_pX_PowerXMod_long_pre(res[0], e, mod)
            sig_off()
        elif ZZ_pX_deg(x[0]) < ZZ_pX_deg(modulus[0]):
            sig_on()
            ZZ_pX_PowerMod_long_pre(res[0], x[0], e, mod)
            sig_off()
        else:
            ZZ_pX_rem_pre(y, x[0], mod)
            sig_on()
            ZZ_pX_PowerMod_long_pre(res[0], y, e, mod)
            sig_off()

cdef inline int celement_gcd(ZZ_pX_c* res, ZZ_pX_c* a, ZZ_pX_c *b, cparent parent) except -2:
    """
    Return the greatest common divisor of this polynomial and ``other``, as
    a monic polynomial.

    INPUT:

    - ``other`` -- a polynomial defined over the same ring as ``self``

    EXAMPLES::

        sage: R.<x> = PolynomialRing(GF(3),implementation="NTL")
        sage: f,g = x + 2, x^2 - 1
        sage: f.gcd(g)
        x + 2

    """
    if parent != NULL:
        parent[0].restore()
    ZZ_pX_GCD(res[0], a[0], b[0])

cdef inline int celement_xgcd(ZZ_pX_c* res, ZZ_pX_c* s, ZZ_pX_c *t, ZZ_pX_c* a, ZZ_pX_c *b, cparent parent) except -2:
    r"""
    Compute the extended gcd of this element and ``other``.

    INPUT:

    - ``other`` -- an element in the same polynomial ring

    OUTPUT:

    A tuple ``(r,s,t)`` of elements in the polynomial ring such
    that ``r = s*self + t*other``.

    EXAMPLES::

        sage: R.<x> = PolynomialRing(GF(3),implementation='NTL')
        sage: x.xgcd(x)
        (x, 0, 1)
        sage: (x^2 - 1).xgcd(x - 1)
        (x + 2, 0, 1)
        sage: R.zero().xgcd(R.one())
        (1, 0, 1)
        sage: (x^3 - 1).xgcd((x - 1)^2)
        (x^2 + x + 1, 0, 1)
        sage: ((x - 1)*(x + 1)).xgcd(x*(x - 1))
        (x + 2, 1, 2)

    """
    if parent != NULL:
        parent[0].restore()
    ZZ_pX_XGCD(res[0], s[0], t[0], a[0], b[0])
