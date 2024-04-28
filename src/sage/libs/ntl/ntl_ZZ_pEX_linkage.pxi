r"""
Linkage for arithmetic with NTL's ZZ_pEX elements.

This file provides the backend for \class{Polynomial_ZZ_pEX} via
templating.

AUTHOR:
    -- Yann Laigle-Chapuy (2010-01): initial version
"""

#*****************************************************************************
#       Copyright (C) 2010 Yann Laigle-Chapuy
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.signals cimport sig_on, sig_off

from sage.libs.ntl.ntl_ZZ_pEContext cimport ntl_ZZ_pEContext_class
from sage.libs.ntl.ZZ_pEX cimport *
from sage.libs.ntl.ntl_ZZ_pE cimport ntl_ZZ_pE
from sage.libs.ntl.types cimport ZZ_pX_c, ZZ_pEX_c

cdef ZZ_pEX_c *celement_new(cparent parent) noexcept:
    """
    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    return new ZZ_pEX_c()

cdef int celement_delete(ZZ_pEX_c *e, cparent parent) noexcept:
    """
    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
        sage: del x
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    del e

cdef int celement_construct(ZZ_pEX_c *e, cparent parent) noexcept:
    """
    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()

cdef int celement_destruct(ZZ_pEX_c *e, cparent parent) noexcept:
    """
    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
        sage: del x
    """
    # do not call restore here
    # 1) the NTL context might have already been destroyed when exiting Python
    # 2) you better not make any NTL calls after destruct, no need to set the context

cdef int celement_gen(ZZ_pEX_c *e, long i, cparent parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    ZZ_pEX_SetX(e[0])

cdef object celement_repr(ZZ_pEX_c *e, cparent parent):
    """
    We ignore NTL's printing.

    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
        sage: x
        x
    """
    raise NotImplementedError

cdef inline int celement_set(ZZ_pEX_c* res, ZZ_pEX_c* a, cparent parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
        sage: y = x
        sage: y
        x
    """
    res[0] = a[0]

cdef inline int celement_set_si(ZZ_pEX_c* res, long i, cparent parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
        sage: P(0)
        0
        sage: P(17)
        17
        sage: P(next_prime(2**60))
        0
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    ZZ_pEX_SetCoeff_long(res[0], 0, i)

cdef inline long celement_get_si(ZZ_pEX_c* res, cparent parent) except -2:
    raise NotImplementedError

cdef inline bint celement_is_zero(ZZ_pEX_c* a, cparent parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
        sage: bool(x), x.is_zero()
        (True, False)
        sage: bool(P(0)), P(0).is_zero()
        (False, True)
    """
    return ZZ_pEX_IsZero(a[0])

cdef inline bint celement_is_one(ZZ_pEX_c *a, cparent parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
        sage: x.is_one()
        False
        sage: P(1).is_one()
        True
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    return ZZ_pEX_IsOne(a[0])

cdef inline bint celement_equal(ZZ_pEX_c *a, ZZ_pEX_c *b, cparent parent) except -2:
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
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    return a[0] == b[0]

cdef inline int celement_cmp(ZZ_pEX_c *a, ZZ_pEX_c *b, cparent parent) except -2:
    """
    Not used.

    Comparison is implemented in
    ``sage/rings/polynomial/polynomial_zz_pex.pyx`` instead.
    """
    raise NotImplementedError

cdef long celement_len(ZZ_pEX_c *a, cparent parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
        sage: x.degree()
        1
        sage: (x+1).degree()
        1
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    return int(ZZ_pEX_deg(a[0]))+1

cdef inline int celement_add(ZZ_pEX_c *res, ZZ_pEX_c *a, ZZ_pEX_c *b, cparent parent) except -2:
    """
    EXAMPLES::

        sage: K.<a> = GF(next_prime(2**60)**3)
        sage: P.<x> = PolynomialRing(K,implementation='NTL')
        sage: (1+a+a^2)*x + (1+x+x^2)
        x^2 + (a^2 + a + 2)*x + 1
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    ZZ_pEX_add(res[0], a[0], b[0])

cdef inline int celement_sub(ZZ_pEX_c* res, ZZ_pEX_c* a, ZZ_pEX_c* b, cparent parent) except -2:
    """
    EXAMPLES::

        sage: K.<a> = GF(next_prime(2**60)**3)
        sage: P.<x> = PolynomialRing(K,implementation='NTL')
        sage: (1+a+a^2)*x - (1+x+x^2)
        1152921504606847008*x^2 + (a^2 + a)*x + 1152921504606847008
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    ZZ_pEX_sub(res[0], a[0], b[0])

cdef inline int celement_neg(ZZ_pEX_c* res, ZZ_pEX_c* a, cparent parent) except -2:
    """
    EXAMPLES::

        sage: K.<a> = GF(next_prime(2**60)**3)
        sage: P.<x> = PolynomialRing(K,implementation='NTL')
        sage: -x
        1152921504606847008*x
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    ZZ_pEX_negate(res[0], a[0])

cdef inline int celement_mul_scalar(ZZ_pEX_c* res, ZZ_pEX_c* p, object c, cparent parent) except -1:
    raise NotImplementedError

cdef inline int celement_mul(ZZ_pEX_c* res, ZZ_pEX_c* a, ZZ_pEX_c* b, cparent parent) except -2:
    """
    EXAMPLES::

        sage: K.<a> = GF(next_prime(2**60)**3)
        sage: P.<x> = PolynomialRing(K,implementation='NTL')
        sage: (1+a+a^2)*x * (1+x+x^2)
        (a^2 + a + 1)*x^3 + (a^2 + a + 1)*x^2 + (a^2 + a + 1)*x
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    ZZ_pEX_mul(res[0], a[0], b[0])

cdef inline int celement_truncate(ZZ_pEX_c* res, ZZ_pEX_c* a, long len, cparent parent) except -2:
    """
    EXAMPLES::

        sage: K.<a> = GF(next_prime(2**60)**3)
        sage: P.<x> = PolynomialRing(K,implementation='NTL')
        sage: p = (a^2 + 1)*x^3 + (a + 1)*x^2 + (a^2 + a + 1)*x + a
        sage: p.truncate(2)   # indirect doctest
        (a^2 + a + 1)*x + a
    """
    ZZ_pEX_trunc(res[0], a[0], len)

cdef inline int celement_div(ZZ_pEX_c* res, ZZ_pEX_c* a, ZZ_pEX_c* b, cparent parent) except -2:
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    return ZZ_pEX_divide(res[0], a[0], b[0])

cdef inline int celement_floordiv(ZZ_pEX_c* res, ZZ_pEX_c* a, ZZ_pEX_c* b, cparent parent) except -2:
    """
    EXAMPLES::

        sage: K.<a> = GF(next_prime(2**60)**3)
        sage: P.<x> = PolynomialRing(K,implementation='NTL')
        sage: (x^2+2*a*x+a^2)//(x+a)
        x + a
        sage: (x^2+2*a*x)//(x+a)
        x + a
        sage: x//(x+1)
        1
        sage: (x+1)//x
        1
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    ZZ_pEX_div_ZZ_pEX(res[0], a[0], b[0])

cdef inline int celement_mod(ZZ_pEX_c* res, ZZ_pEX_c* a, ZZ_pEX_c* b, cparent parent) except -2:
    """
    EXAMPLES::

        sage: K.<a> = GF(next_prime(2**60)**3)
        sage: P.<x> = PolynomialRing(K,implementation='NTL')
        sage: (x^2-2*a*x) % (x+a)
        3*a^2
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    ZZ_pEX_rem(res[0], a[0], b[0])

cdef inline int celement_quorem(ZZ_pEX_c* q, ZZ_pEX_c* r, ZZ_pEX_c* a, ZZ_pEX_c* b, cparent parent) except -2:
    """
    EXAMPLES::

        sage: K.<a> = GF(next_prime(2**60)**3)
        sage: P.<x> = PolynomialRing(K,implementation='NTL')
        sage: (x^2+2*a*x).quo_rem(x-a)
        (x + 3*a, 3*a^2)
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    ZZ_pEX_DivRem(q[0], r[0], a[0], b[0])

cdef inline int celement_inv(ZZ_pEX_c* res, ZZ_pEX_c* a, cparent parent) except -2:
    raise NotImplementedError

cdef inline int celement_pow(ZZ_pEX_c* res, ZZ_pEX_c* x, long e, ZZ_pEX_c *modulus, cparent parent) except -2:
    """
    EXAMPLES::

        sage: K.<a> = GF(next_prime(2**60)**3)
        sage: P.<x> = PolynomialRing(K,implementation='NTL')
        sage: x^1000
        x^1000
        sage: (x+1)^2
        x^2 + 2*x + 1
        sage: (x+1)^(-2)
        1/(x^2 + 2*x + 1)
        sage: f = x+(a+1)
        sage: f**50 == sum(binomial(50,i)*(a+1)**i*x**(50-i) for i in range(51))        # needs sage.symbolic
        True

    TESTS:

    Check that :issue:`15777` is fixed::

        sage: k.<t> = GF(5**5)
        sage: x = polygen(k)
        sage: pow(x+1,100,x)
        1
        sage: pow(x+2,3,x)
        3
        sage: pow(x**3+1,2,x**2+2)
        x + 3
        sage: pow(x**3+1,10**7,x**2+2)
        x + 2
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()

    cdef ZZ_pEX_Modulus_c mod
    cdef ZZ_pEX_c y
    if modulus == NULL:
        if ZZ_pEX_IsX(x[0]):
            sig_on()
            ZZ_pEX_LeftShift(res[0], x[0], e - 1)
            sig_off()
        else:
            sig_on()
            ZZ_pEX_power(res[0], x[0], e)
            sig_off()
    else:
        if ZZ_pEX_deg(modulus[0]) == 1:
            ZZ_pEX_rem(y, x[0], modulus[0])
            sig_on()
            ZZ_pEX_power(res[0], y, e)
            sig_off()
            return 0
        ZZ_pEX_Modulus_build(mod, modulus[0])
        if ZZ_pEX_deg(x[0]) < ZZ_pEX_deg(modulus[0]):
            sig_on()
            ZZ_pEX_PowerMod_pre(res[0], x[0], e, mod)
            sig_off()
        else:
            ZZ_pEX_rem_pre(y, x[0], mod)
            sig_on()
            ZZ_pEX_PowerMod_pre(res[0], y, e, mod)
            sig_off()

cdef inline int celement_gcd(ZZ_pEX_c* res, ZZ_pEX_c* a, ZZ_pEX_c *b, cparent parent) except -2:
    """
    EXAMPLES::

        sage: K.<a> = GF(next_prime(2**60)**3)
        sage: P.<x> = PolynomialRing(K,implementation='NTL')
        sage: f = (x+3)*(x^7+a*x^5+1)
        sage: f.gcd(x+3)
        x + 3
        sage: f.gcd(x+4)
        1
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    ZZ_pEX_GCD(res[0], a[0], b[0])

cdef inline int celement_xgcd(ZZ_pEX_c* res, ZZ_pEX_c* s, ZZ_pEX_c *t, ZZ_pEX_c* a, ZZ_pEX_c *b, cparent parent) except -2:
    """
    EXAMPLES::

        sage: K.<a> = GF(next_prime(2**60)**3)
        sage: P.<x> = PolynomialRing(K,implementation='NTL')
        sage: f = (x+3)*(x^7+a*x^5+1)
        sage: f.xgcd(x+3)
        (x + 3, 0, 1)
        sage: (a+1+x).xgcd(a+x)
        (1, 1, 1152921504606847008)
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    ZZ_pEX_XGCD(res[0], s[0], t[0], a[0], b[0])
