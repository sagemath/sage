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
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/interrupt.pxi"

from sage.libs.ntl.ntl_ZZ_pEContext cimport ntl_ZZ_pEContext_class
from sage.libs.ntl.ntl_ZZ_pEContext_decl cimport *
from sage.libs.ntl.ntl_ZZ_pEX_decl cimport *
from sage.libs.ntl.ntl_ZZ_pE_decl cimport ZZ_pE_from_str
from sage.libs.ntl.ntl_ZZ_pE cimport ntl_ZZ_pE

cdef ZZ_pEX_c *celement_new(cparent parent):
    """
    EXAMPLE:
        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    cdef ZZ_pEX_c *e = ZZ_pEX_new()
    return e

cdef int celement_delete(ZZ_pEX_c *e, cparent parent):
    """
    EXAMPLE:
        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
        sage: del x
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    ZZ_pEX_delete(e)

cdef int celement_construct(ZZ_pEX_c *e, cparent parent):
    """
    EXAMPLE:
        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    ZZ_pEX_construct(e)

cdef int celement_destruct(ZZ_pEX_c *e, cparent parent):
    """
    EXAMPLE:
        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
        sage: del x
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    ZZ_pEX_destruct(e)

cdef int celement_gen(ZZ_pEX_c *e, long i, cparent parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    ZZ_pEX_SetX(e[0])

cdef object celement_repr(ZZ_pEX_c *e, cparent parent):
    """
    We ignore NTL's printing.

    EXAMPLE:
        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
        sage: x
        x
    """
    raise NotImplementedError

cdef inline int celement_set(ZZ_pEX_c* res, ZZ_pEX_c* a, cparent parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
        sage: y = x
        sage: y
        x
    """
    res[0] = a[0]

cdef inline int celement_set_si(ZZ_pEX_c* res, long i, cparent parent) except -2:
    """
    EXAMPLE:
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
    EXAMPLE:
        sage: P.<x> = PolynomialRing(GF(next_prime(2**60)**3,'a'),implementation='NTL')
        sage: bool(x), x.is_zero()
        (True, False)
        sage: bool(P(0)), P(0).is_zero()
        (False, True)
    """
#    if parent != NULL:
#        parent[0].zzpc[0].restore()
#        parent[0].zzpec[0].restore()
    return ZZ_pEX_IsZero(a[0])

cdef inline bint celement_is_one(ZZ_pEX_c *a, cparent parent) except -2:
    """
    EXAMPLE:
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
    EXAMPLE:
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
    return ZZ_pEX_equal(a[0], b[0])

cdef inline int celement_cmp(ZZ_pEX_c *a, ZZ_pEX_c *b, cparent parent) except -2:
    """
    Not used.
    """
    return NotImplementedError

cdef long celement_len(ZZ_pEX_c *a, cparent parent) except -2:
    """
    EXAMPLE:
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
    EXAMPLE:
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
    EXAMPLE:
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
    EXAMPLE:
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
    EXAMPLE:
        sage: K.<a> = GF(next_prime(2**60)**3)
        sage: P.<x> = PolynomialRing(K,implementation='NTL')
        sage: (1+a+a^2)*x * (1+x+x^2)
        (a^2 + a + 1)*x^3 + (a^2 + a + 1)*x^2 + (a^2 + a + 1)*x
    """
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    ZZ_pEX_mul(res[0], a[0], b[0])

cdef inline int celement_div(ZZ_pEX_c* res, ZZ_pEX_c* a, ZZ_pEX_c* b, cparent parent) except -2:
    if parent != NULL:
        parent[0].zzpc[0].restore()
        parent[0].zzpec[0].restore()
    return ZZ_pEX_divide(res[0], a[0], b[0])

cdef inline int celement_floordiv(ZZ_pEX_c* res, ZZ_pEX_c* a, ZZ_pEX_c* b, cparent parent) except -2:
    """
    EXAMPLE:
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
    EXAMPLE:
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
    EXAMPLE:
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
    EXAMPLE::

        sage: K.<a> = GF(next_prime(2**60)**3)
        sage: P.<x> = PolynomialRing(K,implementation='NTL')
        sage: x^1000
        x^1000
        sage: (x+1)^2
        x^2 + 2*x + 1
        sage: (x+1)^(-2)
        1/(x^2 + 2*x + 1)
        sage: f = x+(a+1)
        sage: f**50 == sum(binomial(50,i)*(a+1)**i*x**(50-i) for i in range(51))
        True

    TESTS:

    Check that :trac:`15777` is fixed::

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
    EXAMPLE:
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
    EXAMPLE:
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
