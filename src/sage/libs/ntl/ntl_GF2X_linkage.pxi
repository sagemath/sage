r"""
Linkage for arithmetic with NTL's GF2X elements.

This file provides the backend for \class{Polynomial_GF2X} via
templating.

AUTHOR:
    -- Martin Albrecht (2008-10): initial version
"""

#*****************************************************************************
#       Copyright (C) 2008 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.signals cimport sig_on, sig_off

from sage.libs.ntl.GF2 cimport *
from sage.libs.ntl.GF2X cimport *


cdef GF2X_c *celement_new(long parent) noexcept:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
    """
    return new GF2X_c()

cdef int celement_delete(GF2X_c *e, long parent) noexcept:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: del x
    """
    del e

cdef int celement_construct(GF2X_c *e, long parent) noexcept:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
    """
    pass

cdef int celement_destruct(GF2X_c *e, long parent) noexcept:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: del x
    """
    pass

cdef int celement_gen(GF2X_c *e, long i, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
    """
    cdef unsigned char g = 2
    GF2XFromBytes(e[0], <unsigned char *>(&g), 1)

cdef object celement_repr(GF2X_c *e, long parent):
    """
    We ignore NTL's printing.

    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: x
        x
    """
    raise NotImplementedError

cdef inline int celement_set(GF2X_c* res, GF2X_c* a, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: y = x; y
        x
    """
    res[0] = a[0]

cdef inline int celement_set_si(GF2X_c* res, long i, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: P(0)
        0
        sage: P(2)
        0
        sage: P(1)
        1
    """
    GF2X_conv_long(res[0], i)

cdef inline long celement_get_si(GF2X_c* res, long parent) except -2:
    raise NotImplementedError

cdef inline bint celement_is_zero(GF2X_c* a, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: bool(x), x.is_zero()
        (True, False)
        sage: bool(P(0)), P(0).is_zero()
        (False, True)
    """
    return GF2X_IsZero(a[0])

cdef inline bint celement_is_one(GF2X_c *a, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: x.is_one()
        False
        sage: P(1).is_one()
        True
    """
    return GF2X_IsOne(a[0])

cdef inline bint celement_equal(GF2X_c *a, GF2X_c *b, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: x == x
        True
        sage: y = x; x == y
        True
        sage: x^2 + 1 == x^2 + x
        False
    """
    return a[0] == b[0]

cdef inline int celement_cmp(GF2X_c *a, GF2X_c *b, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: x != 1
        True
        sage: x < 1
        False
        sage: x > 1
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
    cdef bint t
    cdef long diff
    cdef long ca, cb
    diff = GF2X_NumBits(a[0]) - GF2X_NumBits(b[0])
    if diff > 0:
        return 1
    elif diff < 0:
        return -1
    else:
        for i in range(GF2X_NumBits(a[0])-1, -1, -1):
            ca = GF2_conv_to_long(GF2X_coeff(a[0], i))
            cb = GF2_conv_to_long(GF2X_coeff(b[0], i))
            if ca < cb:
                return -1
            elif ca > cb:
                return 1
    return 0

cdef long celement_len(GF2X_c *a, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: x.degree()
        1
        sage: (x+1).degree()
        1
    """
    return int(GF2X_NumBits(a[0]))

cdef inline int celement_add(GF2X_c *res, GF2X_c *a, GF2X_c *b, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: x + 1
        x + 1
    """
    GF2X_add(res[0], a[0], b[0])

cdef inline int celement_sub(GF2X_c* res, GF2X_c* a, GF2X_c* b, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: x - 1
        x + 1
    """
    GF2X_sub(res[0], a[0], b[0])

cdef inline int celement_neg(GF2X_c* res, GF2X_c* a, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: -x
        x
    """
    res[0] = a[0]

cdef inline int celement_mul_scalar(GF2X_c* res, GF2X_c* p, object c,
        long parent) except -1:
    """
    TESTS::

        sage: P.<x> = GF(2)[]
        sage: p = P.random_element()
        sage: 0*p
        0
        sage: 1*p == p
        True
        sage: (3^97)*p == p
        True
    """
    if int(c) == 0:
        GF2X_conv_long(res[0], 0)
    else:
        res[0] = p[0]

cdef inline int celement_mul(GF2X_c* res, GF2X_c* a, GF2X_c* b, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: x*(x+1)
        x^2 + x
    """
    GF2X_mul(res[0], a[0], b[0])

cdef inline int celement_truncate(GF2X_c* res, GF2X_c* a, long len, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: p = x^5 + x^3 + x^2 + x + 1
        sage: p.truncate(2)   # indirect doctest
        x + 1
    """
    GF2X_trunc(res[0], a[0], len)

cdef inline int celement_div(GF2X_c* res, GF2X_c* a, GF2X_c* b, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
    """
    return GF2X_divide(res[0], a[0], b[0])

cdef inline int celement_floordiv(GF2X_c* res, GF2X_c* a, GF2X_c* b, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: x//(x + 1)
        1
        sage: (x + 1)//x
        1
    """
    GF2X_div(res[0], a[0], b[0])

cdef inline int celement_mod(GF2X_c* res, GF2X_c* a, GF2X_c* b, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: (x^2 + 1) % x^2
        1
    """
    GF2X_rem(res[0], a[0], b[0])

cdef inline int celement_quorem(GF2X_c* q, GF2X_c* r, GF2X_c* a, GF2X_c* b, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: f = x^2 + x + 1
        sage: f.quo_rem(x + 1)
        (x, 1)
    """
    GF2X_DivRem(q[0], r[0], a[0], b[0])

cdef inline int celement_inv(GF2X_c* res, GF2X_c* a, long parent) except -2:
    """
    We ignore NTL here and use the fraction field constructor.

    EXAMPLES::

        sage: P.<x> = GF(2)[]
    """
    raise NotImplementedError

cdef inline int celement_pow(GF2X_c* res, GF2X_c* x, long e, GF2X_c *modulus, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: x^1000
        x^1000
        sage: (x+1)^2
        x^2 + 1
        sage: (x+1)^(-2)
        1/(x^2 + 1)
        sage: f = x^9 + x^7 + x^6 + x^5 + x^4 + x^2 + x
        sage: h = x^10 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + 1
        sage: (f^2) % h
        x^9 + x^8 + x^7 + x^5 + x^3
        sage: pow(f, 2, h)
        x^9 + x^8 + x^7 + x^5 + x^3
        sage: pow(x, 1000, h)
        x^8 + x^7 + x^4

    Check that deg x >= deg modulus works (:issue:`35324`)::

        sage: pow(x+1, 2, x^2+x+1)
        x
        sage: pow(x^2+1, 2, x^2+x+1)
        x + 1
    """
    cdef GF2XModulus_c mod
    cdef GF2X_c xmod

    if modulus == NULL:
        if GF2X_IsX(x[0]):
            GF2X_LeftShift(res[0], x[0], e - 1)
        else:
            do_sig = GF2X_deg(x[0]) > 1e5
            if do_sig:
                sig_on()
            GF2X_power(res[0], x[0], e)
            if do_sig:
                sig_off()
    else:
        GF2X_rem(xmod, x[0], modulus[0])
        GF2XModulus_build(mod, modulus[0])
        if GF2X_IsX(xmod):
            GF2X_PowerXMod_long_pre(res[0], e, mod)
        else:
            do_sig = GF2X_deg(x[0]) > 1e5
            if do_sig:
                sig_on()
            GF2X_PowerMod_long_pre(res[0], xmod, e, mod)
            if do_sig:
                sig_off()


cdef inline int celement_gcd(GF2X_c* res, GF2X_c* a, GF2X_c *b, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: f = x*(x+1)
        sage: f.gcd(x+1)
        x + 1
        sage: f.gcd(x^2)
        x
    """
    GF2X_GCD(res[0], a[0], b[0])

cdef inline int celement_xgcd(GF2X_c* res, GF2X_c* s, GF2X_c *t, GF2X_c* a, GF2X_c *b, long parent) except -2:
    """
    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: f = x*(x+1)
        sage: f.xgcd(x+1)
        (x + 1, 0, 1)
        sage: f.xgcd(x^2)
        (x, 1, 1)
    """
    GF2X_XGCD(res[0], s[0], t[0], a[0], b[0])
