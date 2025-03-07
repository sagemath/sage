"""
This module provides the core implementation of multiprecision
floating-point arithmetic. Operations are done in-place.

TESTS:

See if :issue:`15118` is fixed::

    sage: import mpmath
    sage: mpmath.mpf(0)^(-2)
    Traceback (most recent call last):
    ...
    ZeroDivisionError
    sage: mpmath.zeta(2r, -3r)
    Traceback (most recent call last):
    ...
    ZeroDivisionError
"""

#*****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.long cimport *
from cpython.float cimport *
from cpython.complex cimport *
from cpython.number cimport *

from libc.math cimport sqrt as fsqrt
from libc.math cimport frexp

from cysignals.signals cimport sig_check

from sage.ext.stdsage cimport PY_NEW
from sage.libs.gmp.all cimport *
from sage.libs.mpfr cimport *
from sage.rings.integer cimport Integer

from sage.libs.gmp.pylong cimport *

cdef mpz_set_integer(mpz_t v, x):
    if isinstance(x, int):
        mpz_set_pylong(v, x)
    elif isinstance(x, Integer):
        mpz_set(v, (<Integer>x).value)
    else:
        raise TypeError("cannot convert %s to an integer" % x)

cdef inline void mpz_add_si(mpz_t a, mpz_t b, long x) noexcept:
    if x >= 0:
        mpz_add_ui(a, b, x)
    else:
        # careful: overflow when negating INT_MIN
        mpz_sub_ui(a, b, <unsigned long>(-x))

cdef inline mpzi(mpz_t n):
    return mpz_get_pyintlong(n)

cdef inline mpzl(mpz_t n):
    return mpz_get_pylong(n)

# This should be done better
cdef int mpz_tstbit_abs(mpz_t z, unsigned long bit_index) noexcept:
    cdef int res
    if mpz_sgn(z) < 0:
        mpz_neg(z, z)
        res = mpz_tstbit(z, bit_index)
        mpz_neg(z, z)
    else:
        res = mpz_tstbit(z, bit_index)
    return res

cdef void mpz_set_fixed(mpz_t t, MPF *x, int prec, bint abs=False) noexcept:
    """
    Set t = x, or t = |x|, as a fixed-point number with prec bits.
    """
    cdef int offset
    offset = mpz_get_si(x.exp) + prec
    if offset >= 0:
        mpz_mul_2exp(t, x.man, offset)
    else:
        mpz_tdiv_q_2exp(t, x.man, -offset)
    if abs:
        mpz_abs(t, t)

cdef unsigned long mpz_bitcount(mpz_t z) noexcept:
    if mpz_sgn(z) == 0:
        return 0
    return mpz_sizeinbase(z, 2)

# The following limits allowed exponent shifts. We could use mpz_fits_slong_p,
# but then (-LONG_MIN) wraps around; we may also not be able to add large
# shifts safely. A higher limit could be used on 64-bit systems, but
# it is unlikely that anyone will run into this (adding numbers
# that differ by 2^(2^30), at precisions of 2^30 bits).

# Note: MPFR's emax is 1073741823
DEF MAX_SHIFT = 536870912      # 2^29

cdef int mpz_reasonable_shift(mpz_t z) noexcept:
    if mpz_sgn(z) > 0:
        return mpz_cmp_ui(z, MAX_SHIFT) < 0
    else:
        return mpz_cmp_si(z, -MAX_SHIFT) > 0

DEF ROUND_N = 0
DEF ROUND_F = 1
DEF ROUND_C = 2
DEF ROUND_D = 3
DEF ROUND_U = 4

DEF S_NORMAL = 0
DEF S_ZERO = 1
DEF S_NZERO = 2
DEF S_INF = 3
DEF S_NINF = 4
DEF S_NAN = 5

cdef inline str rndmode_to_python(int rnd):
    if rnd == ROUND_N: return 'n'
    if rnd == ROUND_F: return 'f'
    if rnd == ROUND_C: return 'c'
    if rnd == ROUND_D: return 'd'
    if rnd == ROUND_U: return 'u'

cdef inline rndmode_from_python(str rnd):
    if rnd == 'n': return ROUND_N
    if rnd == 'f': return ROUND_F
    if rnd == 'c': return ROUND_C
    if rnd == 'd': return ROUND_D
    if rnd == 'u': return ROUND_U

cdef inline mpfr_rnd_t rndmode_to_mpfr(int rnd) noexcept:
    if rnd == ROUND_N: return MPFR_RNDN
    if rnd == ROUND_F: return MPFR_RNDD
    if rnd == ROUND_C: return MPFR_RNDU
    if rnd == ROUND_D: return MPFR_RNDZ
    if rnd == ROUND_U: return MPFR_RNDA

cdef inline int reciprocal_rnd(int rnd) noexcept:
    if rnd == ROUND_N: return ROUND_N
    if rnd == ROUND_D: return ROUND_U
    if rnd == ROUND_U: return ROUND_D
    if rnd == ROUND_C: return ROUND_F
    if rnd == ROUND_F: return ROUND_C

cdef MPopts opts_exact
cdef MPopts opts_double_precision
cdef MPopts opts_mini_prec

opts_exact.prec = 0
opts_exact.rounding = ROUND_N
opts_double_precision.prec = 53
opts_double_precision.rounding = ROUND_N
opts_mini_prec.prec = 5
opts_mini_prec.rounding = ROUND_D

cdef double _double_inf = float("1e300") * float("1e300")
cdef double _double_ninf = -_double_inf
cdef double _double_nan = _double_inf - _double_inf

cdef inline void MPF_init(MPF *x) noexcept:
    """Allocate space and set value to zero.
    Must be called exactly once when creating a new MPF."""
    x.special = S_ZERO
    mpz_init(x.man)
    mpz_init(x.exp)

cdef inline void MPF_clear(MPF *x) noexcept:
    """Deallocate space. Must be called exactly once when finished with an MPF."""
    mpz_clear(x.man)
    mpz_clear(x.exp)

cdef inline void MPF_set(MPF *dest, MPF *src) noexcept:
    """Clone MPF value. Assumes source value is already normalized."""
    if src is dest:
        return
    dest.special = src.special
    mpz_set(dest.man, src.man)
    mpz_set(dest.exp, src.exp)

cdef inline void MPF_set_zero(MPF *x) noexcept:
    """Set value to 0."""
    x.special = S_ZERO

cdef inline void MPF_set_one(MPF *x) noexcept:
    """Set value to 1."""
    x.special = S_NORMAL
    mpz_set_ui(x.man, 1)
    mpz_set_ui(x.exp, 0)

cdef inline void MPF_set_nan(MPF *x) noexcept:
    """Set value to NaN (not a number)."""
    x.special = S_NAN

cdef inline void MPF_set_inf(MPF *x) noexcept:
    """Set value to +infinity."""
    x.special = S_INF

cdef inline void MPF_set_ninf(MPF *x) noexcept:
    """Set value to -infinity."""
    x.special = S_NINF

cdef MPF_set_si(MPF *x, long n):
    """Set value to that of a given C (long) integer."""
    if n:
        x.special = S_NORMAL
        mpz_set_si(x.man, n)
        mpz_set_ui(x.exp, 0)
        MPF_normalize(x, opts_exact)
    else:
        MPF_set_zero(x)

cdef MPF_set_int(MPF *x, n):
    """Set value to that of a given Python integer."""
    x.special = S_NORMAL
    mpz_set_integer(x.man, n)
    if mpz_sgn(x.man):
        mpz_set_ui(x.exp, 0)
        MPF_normalize(x, opts_exact)
    else:
        MPF_set_zero(x)

cdef MPF_set_man_exp(MPF *x, man, exp):
    """
    Set value to man*2^exp where man, exp may be of any appropriate
    Python integer types.
    """
    x.special = S_NORMAL
    mpz_set_integer(x.man, man)
    mpz_set_integer(x.exp, exp)
    MPF_normalize(x, opts_exact)


# Temporary variables. Note: not thread-safe.
# Used by MPF_add/MPF_sub/MPF_div
cdef mpz_t tmp_exponent
mpz_init(tmp_exponent)
cdef MPF tmp0
MPF_init(&tmp0)

# Used by MPF_hypot and MPF_cmp, which may call MPF_add/MPF_sub
cdef MPF tmp1
MPF_init(&tmp1)
cdef MPF tmp2
MPF_init(&tmp2)


# Constants needed in a few places
cdef MPF MPF_C_1
MPF_init(&MPF_C_1)
MPF_set_si(&MPF_C_1, 1)
cdef Integer MPZ_ZERO = Integer(0)
cdef tuple _mpf_fzero = (0, MPZ_ZERO, 0, 0)
cdef tuple _mpf_fnan = (0, MPZ_ZERO, -123, -1)
cdef tuple _mpf_finf = (0, MPZ_ZERO, -456, -2)
cdef tuple _mpf_fninf = (1, MPZ_ZERO, -789, -3)

cdef MPF_set_tuple(MPF *x, tuple value):
    """
    Set value of an MPF to that of a normalized (sign, man, exp, bc) tuple
    in the format used by mpmath.libmp.
    """
    #cdef int sign
    cdef Integer man
    sign, _man, exp, bc = value
    if isinstance(_man, Integer):
        man = <Integer>_man
    else:
        # This is actually very unlikely; it should never happen
        # in internal code that man isn't an Integer. Maybe the check
        # can be avoided by doing checks in e.g. MPF_set_any?
        man = Integer(_man)
    if mpz_sgn(man.value):
        MPF_set_man_exp(x, man, exp)
        if sign:
            mpz_neg(x.man, x.man)
        return
    if value == _mpf_fzero:
        MPF_set_zero(x)
    elif value == _mpf_finf:
        MPF_set_inf(x)
    elif value == _mpf_fninf:
        MPF_set_ninf(x)
    else:
        MPF_set_nan(x)

cdef MPF_to_tuple(MPF *x):
    """Convert MPF value to (sign, man, exp, bc) tuple."""
    cdef Integer man
    if x.special:
        if x.special == S_ZERO: return _mpf_fzero
        #if x.special == S_NZERO: return _mpf_fnzero
        if x.special == S_INF: return _mpf_finf
        if x.special == S_NINF: return _mpf_fninf
        return _mpf_fnan
    man = PY_NEW(Integer)
    if mpz_sgn(x.man) < 0:
        mpz_neg(man.value, x.man)
        sign = 1
    else:
        mpz_set(man.value, x.man)
        sign = 0
    exp = mpz_get_pyintlong(x.exp)
    bc = mpz_sizeinbase(x.man, 2)
    return (sign, man, exp, bc)

cdef MPF_set_double(MPF *r, double x):
    """
    Set r to the value of a C double x.
    """
    cdef int exp
    cdef double man
    if x != x:
        MPF_set_nan(r)
        return
    if x == _double_inf:
        MPF_set_inf(r)
        return
    if x == _double_ninf:
        MPF_set_ninf(r)
        return
    man = frexp(x, &exp)
    man *= 9007199254740992.0
    mpz_set_d(r.man, man)
    mpz_set_si(r.exp, exp-53)
    r.special = S_NORMAL
    MPF_normalize(r, opts_exact)

import math as pymath

# TODO: implement this function safely without using the Python math module
cdef double MPF_to_double(MPF *x, bint strict) noexcept:
    """Convert MPF value to a Python float."""
    if x.special == S_NORMAL:
        man = mpzi(x.man)
        exp = mpzi(x.exp)
        bc = mpz_sizeinbase(x.man, 2)
        try:
            if bc < 100:
                return pymath.ldexp(man, exp)
            # Try resizing the mantissa. Overflow may still happen here.
            n = bc - 53
            m = man >> n
            return pymath.ldexp(m, exp + n)
        except OverflowError:
            if strict:
                raise
            # Overflow to infinity
            if exp + bc > 0:
                if man < 0:
                    return _double_ninf
                else:
                    return _double_inf
            # Underflow to zero
            return 0.0
    if x.special == S_ZERO:
        return 0.0
    if x.special == S_INF:
        return _double_inf
    if x.special == S_NINF:
        return _double_ninf
    return _double_nan

cdef MPF_to_fixed(mpz_t r, MPF *x, long prec, bint truncate):
    """
    Set r = x, r being in the format of a fixed-point number with prec bits.
    Floor division is used unless truncate=True in which case
    truncating division is used.
    """
    cdef long shift
    if x.special:
        if x.special == S_ZERO or x.special == S_NZERO:
            mpz_set_ui(r, 0)
            return
        raise ValueError("cannot create fixed-point number from special value")
    if mpz_reasonable_shift(x.exp):
        # XXX: signed integer overflow
        shift = mpz_get_si(x.exp) + prec
        if shift >= 0:
            mpz_mul_2exp(r, x.man, shift)
        else:
            if truncate:
                mpz_tdiv_q_2exp(r, x.man, -shift)
            else:
                mpz_fdiv_q_2exp(r, x.man, -shift)
        return
    # Underflow
    if mpz_sgn(x.exp) < 0:
        mpz_set_ui(r, 0)
        return
    raise OverflowError("cannot convert huge number to fixed-point format")

cdef int MPF_sgn(MPF *x) noexcept:
    """
    Give the sign of an MPF (-1, 0, or 1).
    """
    if x.special:
        if x.special == S_INF:
            return 1
        if x.special == S_NINF:
            return -1
        return 0
    return mpz_sgn(x.man)

cdef void MPF_neg(MPF *r, MPF *s) noexcept:
    """
    Set r = -s. MPF_neg(x, x) negates in place.
    """
    if s.special:
        if s.special == S_ZERO:
            r.special = S_ZERO  # r.special = S_NZERO
        elif s.special == S_NZERO:
            r.special = S_ZERO
        elif s.special == S_INF:
            r.special = S_NINF
        elif s.special == S_NINF:
            r.special = S_INF
        else:
            r.special = s.special
        return
    r.special = s.special
    mpz_neg(r.man, s.man)
    if r is not s:
        mpz_set(r.exp, s.exp)

cdef void MPF_abs(MPF *r, MPF *s) noexcept:
    """
    Set r = abs(s). MPF_abs(r, r) sets the absolute value in place.
    """
    if s.special:
        if s.special == S_NINF:
            r.special = S_INF
        else:
            r.special = s.special
        return
    r.special = s.special
    mpz_abs(r.man, s.man)
    if r is not s:
        mpz_set(r.exp, s.exp)

cdef MPF_normalize(MPF *x, MPopts opts):
    """
    Normalize.

    With prec = 0, trailing zero bits are stripped but no rounding
    is performed.
    """
    cdef int sign
    cdef long trail, bc, shift
    if x.special != S_NORMAL:
        return
    sign = mpz_sgn(x.man)
    if sign == 0:
        x.special = S_ZERO
        mpz_set_ui(x.exp, 0)
        return
    bc = mpz_sizeinbase(x.man, 2)
    shift = bc - opts.prec
    # Ok if mantissa small and no trailing zero bits
    if (shift <= 0 or not opts.prec) and mpz_odd_p(x.man):
        return
    # Mantissa is too large, so divide by appropriate power of 2
    # Need to be careful about rounding
    if shift > 0 and opts.prec:
        if opts.rounding == ROUND_N:
            if mpz_tstbit_abs(x.man, shift-1):
                if mpz_tstbit_abs(x.man, shift) or mpz_scan1(x.man, 0) < (shift-1):
                    if sign < 0:
                        mpz_fdiv_q_2exp(x.man, x.man, shift)
                    else:
                        mpz_cdiv_q_2exp(x.man, x.man, shift)
                else:
                    mpz_tdiv_q_2exp(x.man, x.man, shift)
            else:
                mpz_tdiv_q_2exp(x.man, x.man, shift)
        elif opts.rounding == ROUND_D:
            mpz_tdiv_q_2exp(x.man, x.man, shift)
        elif opts.rounding == ROUND_F:
            mpz_fdiv_q_2exp(x.man, x.man, shift)
        elif opts.rounding == ROUND_C:
            mpz_cdiv_q_2exp(x.man, x.man, shift)
        elif opts.rounding == ROUND_U:
            if sign < 0:
                mpz_fdiv_q_2exp(x.man, x.man, shift)
            else:
                mpz_cdiv_q_2exp(x.man, x.man, shift)
        else:
            raise ValueError("bad rounding mode")
    else:
        shift = 0
    # Strip trailing bits
    trail = mpz_scan1(x.man, 0)
    if 0 < trail < bc:
        mpz_tdiv_q_2exp(x.man, x.man, trail)
        shift += trail
    mpz_add_si(x.exp, x.exp, shift)

cdef void MPF_pos(MPF *x, MPF *y, MPopts opts) noexcept:
    """
    Set x = +y (i.e. copy the value, and round if the
    working precision is smaller than the width
    of the mantissa of y).
    """
    MPF_set(x, y)
    MPF_normalize(x, opts)

cdef void _add_special(MPF *r, MPF *s, MPF *t) noexcept:
    if s.special == S_ZERO:
        # (+0) + (-0) = +0
        if t.special == S_NZERO:
            MPF_set(r, s)
        # (+0) + x = x
        else:
            MPF_set(r, t)
    elif t.special == S_ZERO:
        # (-0) + (+0) = +0
        if s.special == S_NZERO:
            MPF_set(r, t)
        # x + (+0) = x
        else:
            MPF_set(r, s)
    # (+/- 0) + x = x
    elif s.special == S_NZERO:
        MPF_set(r, t)
    elif t.special == S_NZERO:
        MPF_set(r, s)
    # (+/- inf) + (-/+ inf) = nan
    elif ((s.special == S_INF and t.special == S_NINF) or
       (s.special == S_NINF and t.special == S_INF)):
        MPF_set_nan(r)
    # nan or +/- inf trumps any finite number
    elif s.special == S_NAN or t.special == S_NAN:
        MPF_set_nan(r)
    elif s.special:
        MPF_set(r, s)
    else:
        MPF_set(r, t)
    return

cdef void _sub_special(MPF *r, MPF *s, MPF *t) noexcept:
    if s.special == S_ZERO:
        # (+0) - (+/-0) = (+0)
        if t.special == S_NZERO:
            MPF_set(r, s)
        else:
            # (+0) - x = (-x)
            MPF_neg(r, t)
    elif t.special == S_ZERO:
        # x - (+0) = x; also covers (-0) - (+0) = (-0)
        MPF_set(r, s)
    # (-0) - x = x
    elif s.special == S_NZERO:
        # (-0) - (-0) = (+0)
        if t.special == S_NZERO:
            MPF_set_zero(r)
        # (-0) - x = -x
        else:
            MPF_neg(r, t)
    elif t.special == S_NZERO:
        # x - (-0) = x
        MPF_set(r, s)
    # (+/- inf) - (+/- inf) = nan
    elif ((s.special == S_INF and t.special == S_INF) or
       (s.special == S_NINF and t.special == S_NINF)):
        MPF_set_nan(r)
    elif s.special == S_NAN or t.special == S_NAN:
        MPF_set_nan(r)
    # nan - x or (+/-inf) - x = l.h.s
    elif s.special:
        MPF_set(r, s)
    # x - nan or x - (+/-inf) = (- r.h.s)
    else:
        MPF_neg(r, t)

cdef void _mul_special(MPF *r, MPF *s, MPF *t) noexcept:
    if s.special == S_ZERO:
        if t.special == S_NORMAL or t.special == S_ZERO:
            MPF_set(r, s)
        elif t.special == S_NZERO:
            MPF_set(r, t)
        else:
            MPF_set_nan(r)
    elif t.special == S_ZERO:
        if s.special == S_NORMAL:
            MPF_set(r, t)
        elif s.special == S_NZERO:
            MPF_set(r, s)
        else:
            MPF_set_nan(r)
    elif s.special == S_NZERO:
        if t.special == S_NORMAL:
            if mpz_sgn(t.man) < 0:
                MPF_set_zero(r)
            else:
                MPF_set(r, s)
        else:
            MPF_set_nan(r)
    elif t.special == S_NZERO:
        if s.special == S_NORMAL:
            if mpz_sgn(s.man) < 0:
                MPF_set_zero(r)
            else:
                MPF_set(r, t)
        else:
            MPF_set_nan(r)
    elif s.special == S_NAN or t.special == S_NAN:
        MPF_set_nan(r)
    else:
        if MPF_sgn(s) == MPF_sgn(t):
            MPF_set_inf(r)
        else:
            MPF_set_ninf(r)

cdef _div_special(MPF *r, MPF *s, MPF *t):
    # TODO: handle signed zeros correctly
    if s.special == S_NAN or t.special == S_NAN:
        MPF_set_nan(r)
    elif t.special == S_ZERO or t.special == S_NZERO:
        raise ZeroDivisionError
    elif s.special == S_ZERO or s.special == S_NZERO:
        MPF_set_zero(r)
    elif s.special == S_NORMAL:
        MPF_set_zero(r)
    elif s.special == S_INF or s.special == S_NINF:
        if t.special == S_INF or t.special == S_NINF:
            MPF_set_nan(r)
        elif MPF_sgn(s) == MPF_sgn(t):
            MPF_set_inf(r)
        else:
            MPF_set_ninf(r)
    # else:
    elif t.special == S_INF or t.special == S_NINF:
        MPF_set_zero(r)

cdef _add_perturbation(MPF *r, MPF *s, int sign, MPopts opts):
    cdef long shift
    if opts.rounding == ROUND_N:
        MPF_set(r, s)
    else:
        shift = opts.prec - mpz_sizeinbase(s.man, 2) + 8
        if shift < 0:
            shift = 8
        mpz_mul_2exp(r.man, s.man, shift)
        mpz_add_si(r.man, r.man, sign)
        mpz_sub_ui(r.exp, s.exp, shift)
        MPF_normalize(r, opts)

cdef MPF_add(MPF *r, MPF *s, MPF *t, MPopts opts):
    """
    Set r = s + t, with exact rounding.

    With prec = 0, the addition is performed exactly. Note that this
    may cause overflow if the exponents are huge.
    """
    cdef long shift, sbc, tbc
    #assert (r is not s) and (r is not t)
    if s.special or t.special:
        _add_special(r, s, t)
        return
    r.special = S_NORMAL
    # Difference between exponents
    mpz_sub(tmp_exponent, s.exp, t.exp)
    if mpz_reasonable_shift(tmp_exponent):
        shift = mpz_get_si(tmp_exponent)
        if shift >= 0:
            # |s| >> |t|
            if shift > 2*opts.prec and opts.prec:
                sbc = mpz_sizeinbase(s.man, 2)
                tbc = mpz_sizeinbase(t.man, 2)
                if shift + sbc - tbc > opts.prec+8:
                    _add_perturbation(r, s, mpz_sgn(t.man), opts)
                    return
            # |s| > |t|
            mpz_mul_2exp(tmp0.man, s.man, shift)
            mpz_add(r.man, tmp0.man, t.man)
            mpz_set(r.exp, t.exp)
            MPF_normalize(r, opts)
        elif shift < 0:
            shift = -shift
            # |s| << |t|
            if shift > 2*opts.prec and opts.prec:
                sbc = mpz_sizeinbase(s.man, 2)
                tbc = mpz_sizeinbase(t.man, 2)
                if shift + tbc - sbc > opts.prec+8:
                    _add_perturbation(r, t, mpz_sgn(s.man), opts)
                    return
            # |s| < |t|
            mpz_mul_2exp(tmp0.man, t.man, shift)
            mpz_add(r.man, tmp0.man, s.man)
            mpz_set(r.exp, s.exp)
            MPF_normalize(r, opts)
    else:
        if not opts.prec:
            raise OverflowError("the exact result does not fit in memory")
        # |s| >>> |t|
        if mpz_sgn(tmp_exponent) > 0:
            _add_perturbation(r, s, mpz_sgn(t.man), opts)
        # |s| <<< |t|
        else:
            _add_perturbation(r, t, mpz_sgn(s.man), opts)

cdef MPF_sub(MPF *r, MPF *s, MPF *t, MPopts opts):
    """
    Set r = s - t, with exact rounding.

    With prec = 0, the addition is performed exactly. Note that this
    may cause overflow if the exponents are huge.
    """
    cdef long shift, sbc, tbc
    #assert (r is not s) and (r is not t)
    if s.special or t.special:
        _sub_special(r, s, t)
        return
    r.special = S_NORMAL
    # Difference between exponents
    mpz_sub(tmp_exponent, s.exp, t.exp)
    if mpz_reasonable_shift(tmp_exponent):
        shift = mpz_get_si(tmp_exponent)
        if shift >= 0:
            # |s| >> |t|
            if shift > 2*opts.prec and opts.prec:
                sbc = mpz_sizeinbase(s.man, 2)
                tbc = mpz_sizeinbase(t.man, 2)
                if shift + sbc - tbc > opts.prec+8:
                    _add_perturbation(r, s, -mpz_sgn(t.man), opts)
                    return
            # |s| > |t|
            mpz_mul_2exp(tmp0.man, s.man, shift)
            mpz_sub(r.man, tmp0.man, t.man)
            mpz_set(r.exp, t.exp)
            MPF_normalize(r, opts)
        elif shift < 0:
            shift = -shift
            # |s| << |t|
            if shift > 2*opts.prec and opts.prec:
                sbc = mpz_sizeinbase(s.man, 2)
                tbc = mpz_sizeinbase(t.man, 2)
                if shift + tbc - sbc > opts.prec+8:
                    _add_perturbation(r, t, -mpz_sgn(s.man), opts)
                    MPF_neg(r, r)
                    return
            # |s| < |t|
            mpz_mul_2exp(tmp0.man, t.man, shift)
            mpz_sub(r.man, s.man, tmp0.man)
            mpz_set(r.exp, s.exp)
            MPF_normalize(r, opts)
    else:
        if not opts.prec:
            raise OverflowError("the exact result does not fit in memory")
        # |s| >>> |t|
        if mpz_sgn(tmp_exponent) > 0:
            _add_perturbation(r, s, -mpz_sgn(t.man), opts)
        # |s| <<< |t|
        else:
            _add_perturbation(r, t, -mpz_sgn(s.man), opts)
            MPF_neg(r, r)

cdef bint MPF_eq(MPF *s, MPF *t) noexcept:
    """
    Evaluates s == t.
    """
    if s.special == S_NAN or t.special == S_NAN:
        return False
    if s.special == t.special:
        if s.special == S_NORMAL:
            return (mpz_cmp(s.man, t.man) == 0) and (mpz_cmp(s.exp, t.exp) == 0)
        else:
            return True
    return False

cdef bint MPF_ne(MPF *s, MPF *t) noexcept:
    """
    Evaluates s != t.
    """
    if s.special == S_NAN or t.special == S_NAN:
        return True
    if s.special == S_NORMAL and t.special == S_NORMAL:
        return (mpz_cmp(s.man, t.man) != 0) or (mpz_cmp(s.exp, t.exp) != 0)
    return s.special != t.special

cdef int MPF_cmp(MPF *s, MPF *t) noexcept:
    """
    Evaluates cmp(s,t). Conventions for nan follow those
    of the mpmath.libmp function.
    """
    cdef long sbc, tbc
    cdef int cm
    if MPF_eq(s, t):
        return 0
    if s.special != S_NORMAL or t.special != S_NORMAL:
        if s.special == S_ZERO: return -MPF_sgn(t)
        if t.special == S_ZERO: return MPF_sgn(s)
        if t.special == S_NAN: return 1
        if s.special == S_INF: return 1
        if t.special == S_NINF: return 1
        return -1
    if mpz_sgn(s.man) != mpz_sgn(t.man):
        if mpz_sgn(s.man) < 0:
            return -1
        else:
            return 1
    if not mpz_cmp(s.exp, t.exp):
        return mpz_cmp(s.man, t.man)
    mpz_add_ui(tmp1.exp, s.exp, mpz_sizeinbase(s.man, 2))
    mpz_add_ui(tmp2.exp, t.exp, mpz_sizeinbase(t.man, 2))
    cm = mpz_cmp(tmp1.exp, tmp2.exp)
    if mpz_sgn(s.man) < 0:
        if cm < 0: return 1
        if cm > 0: return -1
    else:
        if cm < 0: return -1
        if cm > 0: return 1
    MPF_sub(&tmp1, s, t, opts_mini_prec)
    return MPF_sgn(&tmp1)

cdef bint MPF_lt(MPF *s, MPF *t) noexcept:
    """
    Evaluates s < t.
    """
    if s.special == S_NAN or t.special == S_NAN:
        return False
    return MPF_cmp(s, t) < 0

cdef bint MPF_le(MPF *s, MPF *t) noexcept:
    """
    Evaluates s <= t.
    """
    if s.special == S_NAN or t.special == S_NAN:
        return False
    return MPF_cmp(s, t) <= 0

cdef bint MPF_gt(MPF *s, MPF *t) noexcept:
    """
    Evaluates s > t.
    """
    if s.special == S_NAN or t.special == S_NAN:
        return False
    return MPF_cmp(s, t) > 0

cdef bint MPF_ge(MPF *s, MPF *t) noexcept:
    """
    Evaluates s >= t.
    """
    if s.special == S_NAN or t.special == S_NAN:
        return False
    return MPF_cmp(s, t) >= 0

cdef MPF_mul(MPF *r, MPF *s, MPF *t, MPopts opts):
    """
    Set r = s * t, with correct rounding.

    With prec = 0, the multiplication is performed exactly,
    i.e. no rounding is performed.
    """
    if s.special or t.special:
        _mul_special(r, s, t)
    else:
        r.special = S_NORMAL
        mpz_mul(r.man, s.man, t.man)
        mpz_add(r.exp, s.exp, t.exp)
        if opts.prec:
            MPF_normalize(r, opts)

cdef MPF_div(MPF *r, MPF *s, MPF *t, MPopts opts):
    """
    Set r = s / t, with correct rounding.
    """
    cdef int sign
    cdef long sbc, tbc, extra
    cdef mpz_t rem
    #assert (r is not s) and (r is not t)
    if s.special or t.special:
        _div_special(r, s, t)
        return
    r.special = S_NORMAL
    # Division by a power of two <=> shift exponents
    if mpz_cmp_si(t.man, 1) == 0:
        MPF_set(&tmp0, s)
        mpz_sub(tmp0.exp, tmp0.exp, t.exp)
        MPF_normalize(&tmp0, opts)
        MPF_set(r, &tmp0)
        return
    elif mpz_cmp_si(t.man, -1) == 0:
        MPF_neg(&tmp0, s)
        mpz_sub(tmp0.exp, tmp0.exp, t.exp)
        MPF_normalize(&tmp0, opts)
        MPF_set(r, &tmp0)
        return
    sign = mpz_sgn(s.man) != mpz_sgn(t.man)
    # Same strategy as for addition: if there is a remainder, perturb
    # the result a few bits outside the precision range before rounding
    extra = opts.prec - mpz_sizeinbase(s.man,2) + mpz_sizeinbase(t.man,2) + 5
    if extra < 5:
        extra = 5
    mpz_init(rem)
    mpz_mul_2exp(tmp0.man, s.man, extra)
    mpz_tdiv_qr(r.man, rem, tmp0.man, t.man)
    if mpz_sgn(rem):
        mpz_mul_2exp(r.man, r.man, 1)
        if sign:
            mpz_sub_ui(r.man, r.man, 1)
        else:
            mpz_add_ui(r.man, r.man, 1)
        extra  += 1
    mpz_clear(rem)
    mpz_sub(r.exp, s.exp, t.exp)
    mpz_sub_ui(r.exp, r.exp, extra)
    MPF_normalize(r, opts)

cdef int MPF_sqrt(MPF *r, MPF *s, MPopts opts) noexcept:
    """
    Set r = sqrt(s), with correct rounding.
    """
    cdef long shift
    cdef mpz_t rem
    #assert r is not s
    if s.special:
        if s.special == S_ZERO or s.special == S_INF:
            MPF_set(r, s)
        else:
            MPF_set_nan(r)
        return 0
    if mpz_sgn(s.man) < 0:
        MPF_set_nan(r)
        return 1
    r.special = S_NORMAL
    if mpz_odd_p(s.exp):
        mpz_sub_ui(r.exp, s.exp, 1)
        mpz_mul_2exp(r.man, s.man, 1)
    elif mpz_cmp_ui(s.man, 1) == 0:
        # Square of a power of two
        mpz_set_ui(r.man, 1)
        mpz_tdiv_q_2exp(r.exp, s.exp, 1)
        MPF_normalize(r, opts)
        return 0
    else:
        mpz_set(r.man, s.man)
        mpz_set(r.exp, s.exp)
    shift = 2*opts.prec - mpz_sizeinbase(r.man,2) + 4
    if shift < 4:
        shift = 4
    shift += shift & 1
    mpz_mul_2exp(r.man, r.man, shift)
    if opts.rounding == ROUND_F or opts.rounding == ROUND_D:
        mpz_sqrt(r.man, r.man)
    else:
        mpz_init(rem)
        mpz_sqrtrem(r.man, rem, r.man)
        if mpz_sgn(rem):
            mpz_mul_2exp(r.man, r.man, 1)
            mpz_add_ui(r.man, r.man, 1)
            shift += 2
        mpz_clear(rem)
    mpz_add_si(r.exp, r.exp, -shift)
    mpz_tdiv_q_2exp(r.exp, r.exp, 1)
    MPF_normalize(r, opts)
    return 0

cdef MPF_hypot(MPF *r, MPF *a, MPF *b, MPopts opts):
    """
    Set r = sqrt(a^2 + b^2)
    """
    cdef MPopts tmp_opts
    if a.special == S_ZERO:
        MPF_abs(r, b)
        MPF_normalize(r, opts)
        return
    if b.special == S_ZERO:
        MPF_abs(r, a)
        MPF_normalize(r, opts)
        return
    tmp_opts = opts
    tmp_opts.prec += 30
    MPF_mul(&tmp1, a, a, opts_exact)
    MPF_mul(&tmp2, b, b, opts_exact)
    MPF_add(r, &tmp1, &tmp2, tmp_opts)
    MPF_sqrt(r, r, opts)

cdef MPF_pow_int(MPF *r, MPF *x, mpz_t n, MPopts opts):
    """
    Set r = x ** n. Currently falls back to mpmath.libmp
    unless n is tiny.
    """
    cdef long m, absm
    cdef unsigned long bc
    cdef int nsign
    if x.special != S_NORMAL:
        nsign = mpz_sgn(n)
        if x.special == S_ZERO:
            if nsign < 0:
                raise ZeroDivisionError
            elif nsign == 0:
                MPF_set(r, &MPF_C_1)
            else:
                MPF_set_zero(r)
        elif x.special == S_INF:
            if nsign > 0:
                MPF_set(r, x)
            elif nsign == 0:
                MPF_set_nan(r)
            else:
                MPF_set_zero(r)
        elif x.special == S_NINF:
            if nsign > 0:
                if mpz_odd_p(n):
                    MPF_set(r, x)
                else:
                    MPF_neg(r, x)
            elif nsign == 0:
                MPF_set_nan(r)
            else:
                MPF_set_zero(r)
        else:
            MPF_set_nan(r)
        return
    bc = mpz_sizeinbase(r.man,2)
    r.special = S_NORMAL
    if mpz_reasonable_shift(n):
        m = mpz_get_si(n)
        if m == 0:
            MPF_set(r, &MPF_C_1)
            return
        if m == 1:
            MPF_set(r, x)
            MPF_normalize(r, opts)
            return
        if m == 2:
            MPF_mul(r, x, x, opts)
            return
        if m == -1:
            MPF_div(r, &MPF_C_1, x, opts)
            return
        if m == -2:
            MPF_mul(r, x, x, opts_exact)
            MPF_div(r, &MPF_C_1, r, opts)
            return
        absm = abs(m)
        if bc * absm < 10000:
            mpz_pow_ui(r.man, x.man, absm)
            mpz_mul_ui(r.exp, x.exp, absm)
            if m < 0:
                MPF_div(r, &MPF_C_1, r, opts)
            else:
                MPF_normalize(r, opts)
            return
    r.special = S_NORMAL
    # (2^p)^n
    if mpz_cmp_si(x.man, 1) == 0:
        mpz_set(r.man, x.man)
        mpz_mul(r.exp, x.exp, n)
        return
    # (-2^p)^n
    if mpz_cmp_si(x.man, -1) == 0:
        if mpz_odd_p(n):
            mpz_set(r.man, x.man)
        else:
            mpz_neg(r.man, x.man)
        mpz_mul(r.exp, x.exp, n)
        return
    # TODO: implement efficiently here
    import mpmath.libmp
    MPF_set_tuple(r,
        mpmath.libmp.mpf_pow_int(MPF_to_tuple(x), mpzi(n),
        opts.prec, rndmode_to_python(opts.rounding)))

cdef mpz_t _pi_value
cdef int _pi_prec = -1

cdef mpz_t _ln2_value
cdef int _ln2_prec = -1

cdef mpz_set_pi(mpz_t x, int prec):
    """
    Set x = pi as a fixed-point number.
    """
    global _pi_value
    global _pi_prec
    if prec <= _pi_prec:
        mpz_tdiv_q_2exp(x, _pi_value, _pi_prec-prec)
    else:
        from mpmath.libmp import pi_fixed
        if _pi_prec < 0:
            mpz_init(_pi_value)
        mpz_set_integer(_pi_value, pi_fixed(prec))
        mpz_set(x, _pi_value)
        _pi_prec = prec

cdef mpz_set_ln2(mpz_t x, int prec):
    """
    Set x = ln(2) as a fixed-point number.
    """
    global _ln2_value
    global _ln2_prec
    if prec <= _ln2_prec:
        mpz_tdiv_q_2exp(x, _ln2_value, _ln2_prec-prec)
    else:
        from mpmath.libmp import ln2_fixed
        if _ln2_prec < 0:
            mpz_init(_ln2_value)
        mpz_set_integer(_ln2_value, ln2_fixed(prec))
        mpz_set(x, _ln2_value)
        _ln2_prec = prec

cdef void _cy_exp_mpfr(mpz_t y, mpz_t x, int prec) noexcept:
    """
    Compute y = exp(x) for fixed-point numbers y and x using MPFR,
    assuming that no overflow will occur.
    """
    cdef mpfr_t yf, xf
    mpfr_init2(xf, mpz_bitcount(x)+2)
    mpfr_init2(yf, prec+2)
    mpfr_set_z(xf, x, MPFR_RNDN)
    mpfr_div_2exp(xf, xf, prec, MPFR_RNDN)
    mpfr_exp(yf, xf, MPFR_RNDN)
    mpfr_mul_2exp(yf, yf, prec, MPFR_RNDN)
    mpfr_get_z(y, yf, MPFR_RNDN)
    mpfr_clear(yf)
    mpfr_clear(xf)

cdef cy_exp_basecase(mpz_t y, mpz_t x, int prec):
    """
    Compute y = exp(x) for fixed-point numbers y and x, assuming
    that x is small (|x| ~< 1). At small precisions, this function
    is equivalent to the exp_basecase function in
    mpmath.libmp.exp_fixed.
    """
    cdef int k, r, u
    cdef mpz_t s0, s1, x2, a
    # TODO: could use custom implementation here; for now switch to MPFR
    if prec > 2000:
        _cy_exp_mpfr(y, x, prec)
        return
    mpz_init(s0)
    mpz_init(s1)
    mpz_init(x2)
    mpz_init(a)
    r = <int>fsqrt(prec)
    prec += r
    mpz_set_ui(s0, 1)
    mpz_mul_2exp(s0, s0, prec)
    mpz_set(s1, s0)
    k = 2
    mpz_mul(x2, x, x)
    mpz_fdiv_q_2exp(x2, x2, prec)
    mpz_set(a, x2)
    while mpz_sgn(a):
        sig_check()
        mpz_fdiv_q_ui(a, a, k)
        mpz_add(s0, s0, a)
        k += 1
        mpz_fdiv_q_ui(a, a, k)
        mpz_add(s1, s1, a)
        k += 1
        mpz_mul(a, a, x2)
        mpz_fdiv_q_2exp(a, a, prec)
    mpz_mul(s1, s1, x)
    mpz_fdiv_q_2exp(s1, s1, prec)
    mpz_add(s0, s0, s1)
    u = r
    while r:
        sig_check()
        mpz_mul(s0, s0, s0)
        mpz_fdiv_q_2exp(s0, s0, prec)
        r -= 1
    mpz_fdiv_q_2exp(y, s0, u)
    mpz_clear(s0)
    mpz_clear(s1)
    mpz_clear(x2)
    mpz_clear(a)


cdef MPF_exp(MPF *y, MPF *x, MPopts opts):
    """
    Set y = exp(x).
    """
    cdef bint sign, is_int
    cdef long wp, wpmod, offset, mag
    cdef mpz_t t, u
    cdef tuple w
    if x.special:
        if x.special == S_ZERO:
            MPF_set_si(y, 1)
        elif x.special == S_NINF:
            MPF_set_zero(y)
        elif x.special == S_INF:
            MPF_set_inf(y)
        else:
            MPF_set_nan(y)
        return
    wp = opts.prec + 14
    sign = mpz_sgn(x.man) < 0
    is_int = mpz_sgn(x.exp) >= 0
    # note: bogus if not reasonable shift
    mag = mpz_bitcount(x.man) + mpz_get_si(x.exp)
    if (not mpz_reasonable_shift(x.exp)) or mag < -wp:
        if mpz_sgn(x.exp) <= 0:
            # perturb
            MPF_set_one(y)
            if opts.rounding != ROUND_N:
                mpz_mul_2exp(y.man, y.man, wp)
                if sign:
                    mpz_sub_ui(y.man, y.man, 1)
                else:
                    mpz_add_ui(y.man, y.man, 1)
                mpz_set_si(y.exp, -wp)
                MPF_normalize(y, opts)
            return
        else:
            raise OverflowError("exp of a huge number")
    #offset = mpz_get_si(x.exp) + wp
    mpz_init(t)
    if mag > 1:
        wpmod = wp + mag
        mpz_set_fixed(t, x, wpmod, False)
        mpz_init(u)
        mpz_set_ln2(u, wpmod)
        # y.exp, t = divmod(t, ln2)
        mpz_fdiv_qr(y.exp,t,t,u)
        mpz_clear(u)
        mpz_fdiv_q_2exp(t, t, mag)
    else:
        mpz_set_fixed(t, x, wp, False)
        mpz_set_ui(y.exp, 0)
    cy_exp_basecase(y.man, t, wp)
    mpz_add_si(y.exp, y.exp, -wp)
    y.special = S_NORMAL
    mpz_clear(t)
    MPF_normalize(y, opts)


cdef MPF_complex_sqrt(MPF *c, MPF *d, MPF *a, MPF *b, MPopts opts):
    """
    Set c+di = sqrt(a+bi).

    c, a and d, b may be the same objects.
    """
    cdef int apos, bneg
    cdef MPF t, u, v
    cdef MPopts wpopts
    if b.special == S_ZERO:
        if a.special == S_ZERO:
            MPF_set_zero(c)
            MPF_set_zero(d)
        # a+bi, a < 0, b = 0
        elif MPF_sgn(a) < 0:
            MPF_abs(d, a)
            MPF_sqrt(d, d, opts)
            MPF_set_zero(c)
        # a > 0
        else:
            MPF_sqrt(c, a, opts)
            MPF_set_zero(d)
        return
    wpopts.prec = opts.prec + 20
    wpopts.rounding = ROUND_D
    MPF_init(&t)
    MPF_init(&u)
    MPF_init(&v)
    apos = MPF_sgn(a) >= 0
    bneg = MPF_sgn(b) <= 0
    if apos:
        # real part
        MPF_hypot(&t, a, b, wpopts)  # t = abs(a+bi) + a
        MPF_add(&t, &t, a, wpopts)
        MPF_set(&u, &t)
        mpz_sub_ui(u.exp, u.exp, 1)  # u = t / 2
        MPF_sqrt(c, &u, opts)   # re = sqrt(u)
        # imag part
        mpz_add_ui(t.exp, t.exp, 1)  # t = 2*t
        MPF_sqrt(&u, &t, wpopts)   # u = sqrt(t)
        MPF_div(d, b, &u, opts) # im = b / u
    else:
        MPF_set(&v, b)
        MPF_hypot(&t, a, b, wpopts)  # t = abs(a+bi) - a
        MPF_sub(&t, &t, a, wpopts)
        MPF_set(&u, &t)
        mpz_sub_ui(u.exp, u.exp, 1)  # u = t / 2
        MPF_sqrt(d, &u, opts)  # im = sqrt(u)
        mpz_add_ui(t.exp, t.exp, 1)  # t = 2*t
        MPF_sqrt(&u, &t, wpopts)   # u = sqrt(t)
        MPF_div(c, &v, &u, opts)  # re = b / u
        if bneg:
            MPF_neg(c, c)
            MPF_neg(d, d)
    MPF_clear(&t)
    MPF_clear(&u)
    MPF_clear(&v)

cdef int MPF_get_mpfr_overflow(mpfr_t y, MPF *x) noexcept:
    """
    Store the mpmath number x exactly in the MPFR variable y. The precision
    of y will be adjusted if necessary. If the exponent overflows, only
    the mantissa is stored and 1 is returned; if no overflow occurs,
    the function returns 0.
    """
    cdef long prec, exp
    if x.special != S_NORMAL:
        if x.special == S_ZERO:
            mpfr_set_ui(y, 0, MPFR_RNDN)
        elif x.special == S_INF:
            mpfr_set_inf(y, 1)
        elif x.special == S_NINF:
            mpfr_set_inf(y, -1)
        else:
            mpfr_set_nan(y)
        return 0
    prec = mpz_bitcount(x.man)
    # Minimum precision for MPFR
    if prec < 2:
        prec = 2
    mpfr_set_prec(y, prec)
    mpfr_set_z(y, x.man, MPFR_RNDN)
    if mpz_reasonable_shift(x.exp):
        exp = mpz_get_si(x.exp)
        if exp >= 0:
            mpfr_mul_2exp(y, y, exp, MPFR_RNDN)
        else:
            mpfr_div_2exp(y, y, -exp, MPFR_RNDN)
        return 0
    else:
        return 1

cdef MPF_set_mpfr(MPF *y, mpfr_t x, MPopts opts):
    """
    Convert the MPFR number x to a normalized MPF y.
    inf/nan and zero are handled.
    """
    cdef long exp
    # TODO: use mpfr_regular_p with MPFR 3
    if mpfr_nan_p(x):
        MPF_set_nan(y)
        return
    if mpfr_inf_p(x):
        if mpfr_sgn(x) > 0:
            MPF_set_inf(y)
        else:
            MPF_set_ninf(y)
        return
    if mpfr_zero_p(x):
        MPF_set_zero(y)
        return
    exp = mpfr_get_z_exp(y.man, x)
    mpz_set_si(y.exp, exp)
    y.special = S_NORMAL
    MPF_normalize(y, opts)

cdef int MPF_log(MPF *y, MPF *x, MPopts opts) noexcept:
    """
    Set y = log(|x|). Returns 1 if x is negative.
    """
    cdef MPF t
    cdef bint negative, overflow
    cdef mpfr_rnd_t rndmode
    cdef mpfr_t yy, xx
    if x.special != S_NORMAL:
        if x.special == S_ZERO:
            MPF_set_ninf(y)
            return 0
        if x.special == S_INF:
            MPF_set_inf(y)
            return 0
        if x.special == S_NAN:
            MPF_set_nan(y)
            return 0
        if x.special == S_NINF:
            MPF_set_inf(y)
            return 1

    negative = MPF_sgn(x) < 0
    mpfr_init2(xx, opts.prec)
    mpfr_init2(yy, opts.prec)

    overflow = MPF_get_mpfr_overflow(xx, x)
    rndmode = rndmode_to_mpfr(opts.rounding)

    if overflow:
        MPF_init(&t)
        # Copy x exponent in case x and y are aliased
        mpz_set(t.exp, x.exp)

        # log(m * 2^e) = log(m) + e*log(2)
        mpfr_abs(xx, xx, MPFR_RNDN)
        mpfr_log(yy, xx, rndmode)
        MPF_set_mpfr(y, yy, opts)

        mpz_set_ln2(t.man, opts.prec+20)
        mpz_mul(t.man, t.man, t.exp)
        mpz_set_si(t.exp, -(opts.prec+20))
        t.special = S_NORMAL

        MPF_add(y, y, &t, opts)
        MPF_clear(&t)
    else:
        mpfr_abs(xx, xx, MPFR_RNDN)
        mpfr_log(yy, xx, rndmode)
        MPF_set_mpfr(y, yy, opts)

    mpfr_clear(xx)
    mpfr_clear(yy)
    return negative

cdef MPF_set_pi(MPF *x, MPopts opts):
    """
    Set x = pi.
    """
    x.special = S_NORMAL
    mpz_set_pi(x.man, (opts.prec+20))
    mpz_set_si(x.exp, -(opts.prec+20))
    MPF_normalize(x, opts)

cdef MPF_set_ln2(MPF *x, MPopts opts):
    """
    Set x = ln(2).
    """
    x.special = S_NORMAL
    mpz_set_ln2(x.man, (opts.prec+20))
    mpz_set_si(x.exp, -(opts.prec+20))
    MPF_normalize(x, opts)


def exp_fixed(Integer x, int prec, ln2=None):
    """
    Return a fixed-point approximation of exp(x) where x is a fixed-point
    number.

    EXAMPLES::

        sage: from sage.libs.mpmath.ext_impl import exp_fixed
        sage: y = exp_fixed(1<<53, 53)
        sage: float(y) / 2^53
        2.718281828459044
    """
    cdef Integer v
    cdef mpz_t n, t
    cdef long nn
    mpz_init(n)
    mpz_init(t)
    if ln2 is None:
        mpz_set_ln2(t, prec)
        mpz_fdiv_qr(n, t, x.value, t)
    else:
        mpz_fdiv_qr(n, t, x.value, (<Integer>ln2).value)
    nn = mpz_get_si(n)
    v = PY_NEW(Integer)
    cy_exp_basecase(v.value, t, prec)
    if nn >= 0:
        mpz_mul_2exp(v.value, v.value, nn)
    else:
        mpz_fdiv_q_2exp(v.value, v.value, -nn)
    mpz_clear(t)
    mpz_clear(n)
    return v


def cos_sin_fixed(Integer x, int prec, pi2=None):
    """
    Return fixed-point approximations of cos(x), sin(x) where
    x is a fixed-point number.

    EXAMPLES::

        sage: from sage.libs.mpmath.ext_impl import cos_sin_fixed
        sage: c, s = cos_sin_fixed(1<<53, 53)
        sage: float(c) / 2^53
        0.5403023058681398
        sage: float(s) / 2^53
        0.8414709848078965
    """
    cdef Integer cv, sv
    cdef mpfr_t t, cf, sf
    mpfr_init2(t, mpz_bitcount(x.value)+2)
    mpfr_init2(cf, prec)
    mpfr_init2(sf, prec)
    mpfr_set_z(t, x.value, MPFR_RNDN)
    mpfr_div_2exp(t, t, prec, MPFR_RNDN)
    mpfr_sin_cos(sf, cf, t, MPFR_RNDN)
    mpfr_mul_2exp(cf, cf, prec, MPFR_RNDN)
    mpfr_mul_2exp(sf, sf, prec, MPFR_RNDN)
    cv = PY_NEW(Integer)
    sv = PY_NEW(Integer)
    mpfr_get_z(cv.value, cf, MPFR_RNDN)
    mpfr_get_z(sv.value, sf, MPFR_RNDN)
    mpfr_clear(t)
    mpfr_clear(cf)
    mpfr_clear(sf)
    return cv, sv


DEF MAX_LOG_INT_CACHE = 2000

cdef mpz_t log_int_cache[MAX_LOG_INT_CACHE+1]
cdef long log_int_cache_prec[MAX_LOG_INT_CACHE+1]
cdef bint log_int_cache_initialized = 0

cdef mpz_log_int(mpz_t v, mpz_t n, int prec):
    """
    Set v = log(n) where n is an integer and v is a fixed-point number
    with the specified precision.
    """
    cdef mpfr_t f
    mpfr_init2(f, prec+15)
    mpfr_set_z(f, n, MPFR_RNDN)
    mpfr_log(f, f, MPFR_RNDN)
    mpfr_mul_2exp(f, f, prec, MPFR_RNDN)
    mpfr_get_z(v, f, MPFR_RNDN)
    mpfr_clear(f)


def log_int_fixed(n, long prec, ln2=None):
    """
    Return fixed-point approximation of log(n).

    EXAMPLES::

        sage: from sage.libs.mpmath.ext_impl import log_int_fixed
        sage: float(log_int_fixed(5, 53)) / 2^53
        1.6094379124341003
        sage: float(log_int_fixed(5, 53)) / 2^53   # exercise cache
        1.6094379124341003
    """
    global log_int_cache_initialized
    cdef Integer t
    cdef int i
    t = PY_NEW(Integer)
    mpz_set_integer(t.value, n)
    if mpz_sgn(t.value) <= 0:
        mpz_set_ui(t.value, 0)
    elif mpz_cmp_ui(t.value, MAX_LOG_INT_CACHE) <= 0:
        if not log_int_cache_initialized:
            for i in range(MAX_LOG_INT_CACHE+1):
                mpz_init(log_int_cache[i])
                log_int_cache_prec[i] = 0
            log_int_cache_initialized = 1
        i = mpz_get_si(t.value)
        if log_int_cache_prec[i] < prec:
            mpz_log_int(log_int_cache[i], t.value, prec+64)
            log_int_cache_prec[i] = prec+64
        mpz_tdiv_q_2exp(t.value, log_int_cache[i], log_int_cache_prec[i]-prec)

    else:
        mpz_log_int(t.value, t.value, prec)
    return t


cdef _MPF_cos_python(MPF *c, MPF *x, MPopts opts):
    """
    Compute c = cos(x) by calling the mpmath.libmp Python implementation.
    """
    from mpmath.libmp.libelefun import mpf_cos_sin
    ct = mpf_cos_sin(MPF_to_tuple(x), opts.prec,
            rndmode_to_python(opts.rounding), 1, False)
    MPF_set_tuple(c, ct)

cdef _MPF_sin_python(MPF *s, MPF *x, MPopts opts):
    """
    Compute s = sin(x) by calling the mpmath.libmp Python implementation.
    """
    from mpmath.libmp.libelefun import mpf_cos_sin
    st = mpf_cos_sin(MPF_to_tuple(x), opts.prec,
            rndmode_to_python(opts.rounding), 2, False)
    MPF_set_tuple(s, st)


cdef MPF_cos(MPF *c, MPF *x, MPopts opts):
    """
    Set c = cos(x)
    """
    cdef mpfr_t cf, xf
    cdef bint overflow
    if x.special != S_NORMAL:
        if x.special == S_ZERO:
            MPF_set_one(c)
        else:
            MPF_set_nan(c)
        return
    mpfr_init(xf)
    mpfr_init2(cf, opts.prec)
    overflow = MPF_get_mpfr_overflow(xf, x)
    if overflow or opts.rounding == ROUND_U:
        _MPF_cos_python(c, x, opts)
    else:
        mpfr_cos(cf, xf, rndmode_to_mpfr(opts.rounding))
        MPF_set_mpfr(c, cf, opts)
    mpfr_clear(xf)
    mpfr_clear(cf)

cdef MPF_sin(MPF *s, MPF *x, MPopts opts):
    """
    Set s = sin(x)
    """
    cdef mpfr_t sf, xf
    cdef bint overflow
    if x.special != S_NORMAL:
        if x.special == S_ZERO:
            MPF_set_zero(s)
        else:
            MPF_set_nan(s)
        return
    mpfr_init(xf)
    mpfr_init2(sf, opts.prec)
    overflow = MPF_get_mpfr_overflow(xf, x)
    if overflow or opts.rounding == ROUND_U:
        _MPF_sin_python(s, x, opts)
    else:
        mpfr_sin(sf, xf, rndmode_to_mpfr(opts.rounding))
        MPF_set_mpfr(s, sf, opts)
    mpfr_clear(xf)
    mpfr_clear(sf)

cdef MPF_cos_sin(MPF *c, MPF *s, MPF *x, MPopts opts):
    """
    Set c = cos(x), s = sin(x)
    """
    cdef mpfr_t cf, sf, xf
    cdef bint overflow
    if x.special != S_NORMAL:
        if x.special == S_ZERO:
            MPF_set_one(c)
            MPF_set_zero(s)
        else:
            MPF_set_nan(c)
            MPF_set_nan(s)
        return
    mpfr_init(xf)
    mpfr_init2(sf, opts.prec)
    mpfr_init2(cf, opts.prec)
    overflow = MPF_get_mpfr_overflow(xf, x)
    if overflow or opts.rounding == ROUND_U:
        _MPF_cos_python(c, x, opts)
        _MPF_sin_python(s, x, opts)
    else:
        mpfr_sin_cos(sf, cf, xf, rndmode_to_mpfr(opts.rounding))
        MPF_set_mpfr(s, sf, opts)
        MPF_set_mpfr(c, cf, opts)
    mpfr_clear(xf)
    mpfr_clear(cf)
    mpfr_clear(sf)


cdef MPF_complex_exp(MPF *re, MPF *im, MPF *a, MPF *b, MPopts opts):
    """
    Set re+im*i = exp(a+bi)
    """
    cdef MPF mag, c, s
    cdef MPopts wopts
    if a.special == S_ZERO:
        MPF_cos_sin(re, im, b, opts)
        return
    if b.special == S_ZERO:
        MPF_exp(re, a, opts)
        MPF_set_zero(im)
        return
    MPF_init(&mag)
    MPF_init(&c)
    MPF_init(&s)
    wopts = opts
    wopts.prec += 4
    MPF_exp(&mag, a, wopts)
    MPF_cos_sin(&c, &s, b, wopts)
    MPF_mul(re, &mag, &c, opts)
    MPF_mul(im, &mag, &s, opts)
    MPF_clear(&mag)
    MPF_clear(&c)
    MPF_clear(&s)

cdef int MPF_pow(MPF *z, MPF *x, MPF *y, MPopts opts) except -1:
    """
    Set z = x^y for real x and y and returns 0 if the result is real-valued.
    If the result is complex, does nothing and returns 1.
    """
    cdef MPopts wopts
    cdef mpz_t t
    cdef MPF w
    cdef int xsign, ysign
    cdef mpz_t tm

    # Integer exponentiation, if reasonable
    if y.special == S_NORMAL and mpz_sgn(y.exp) >= 0:
        mpz_init(tm)
        # check if size is reasonable
        mpz_add_ui(tm, y.exp, mpz_bitcount(y.man))
        mpz_abs(tm, tm)
        if mpz_cmp_ui(tm, 10000) < 0:
            # man * 2^exp
            mpz_mul_2exp(tm, y.man, mpz_get_ui(y.exp))
            MPF_pow_int(z, x, tm, opts)
            mpz_clear(tm)
            return 0
        mpz_clear(tm)

    # x ^ 0
    if y.special == S_ZERO:
        if x.special == S_NORMAL or x.special == S_ZERO:
            MPF_set_one(z)
        else:
            MPF_set_nan(z)
        return 0

    xsign = MPF_sgn(x)
    ysign = MPF_sgn(y)

    if xsign < 0:
        return 1

    # Square root or integer power thereof
    if y.special == S_NORMAL and mpz_cmp_si(y.exp, -1) == 0:
        # x^(1/2)
        if mpz_cmp_ui(y.man, 1) == 0:
            MPF_sqrt(z, x, opts)
            return 0
        # x^(-1/2)
        if mpz_cmp_si(y.man, -1) == 0:
            wopts = opts
            wopts.prec += 10
            wopts.rounding = reciprocal_rnd(wopts.rounding)
            MPF_sqrt(z, x, wopts)
            MPF_div(z, &MPF_C_1, z, opts)
            return 0
        # x^(n/2)
        wopts = opts
        wopts.prec += 10
        if mpz_sgn(y.man) < 0:
            wopts.rounding = reciprocal_rnd(wopts.rounding)
        mpz_init_set(t, y.man)
        MPF_sqrt(z, x, wopts)
        MPF_pow_int(z, z, t, opts)
        mpz_clear(t)
        return 0

    if x.special != S_NORMAL or y.special != S_NORMAL:
        if x.special == S_NAN or y.special == S_NAN:
            MPF_set_nan(z)
            return 0
        if y.special == S_ZERO:
            if x.special == S_NORMAL or x.special == S_ZERO:
                MPF_set_one(z)
            else:
                MPF_set_nan(z)
            return 0
        if x.special == S_ZERO and y.special == S_NORMAL:
            if mpz_sgn(y.man) > 0:
                MPF_set_zero(z)
            return 0

    wopts = opts
    wopts.prec += 10
    MPF_init(&w)
    MPF_log(&w, x, wopts)
    MPF_mul(&w, &w, y, opts_exact)
    MPF_exp(z, &w, opts)
    MPF_clear(&w)
    return 0

cdef MPF_complex_square(MPF *re, MPF *im, MPF *a, MPF *b, MPopts opts):
    """
    Set re+im*i = (a+bi)^2 = a^2-b^2, 2ab*i.
    """
    cdef MPF t, u
    MPF_init(&t)
    MPF_init(&u)
    MPF_mul(&t,a,a,opts_exact)
    MPF_mul(&u,b,b,opts_exact)
    MPF_sub(re, &t, &u, opts)
    MPF_mul(im, a, b, opts)
    if im.special == S_NORMAL:
        mpz_add_ui(im.exp, im.exp, 1)
    MPF_clear(&t)
    MPF_clear(&u)


cdef MPF_complex_reciprocal(MPF *re, MPF *im, MPF *a, MPF *b, MPopts opts):
    """
    Set re+im*i = 1/(a+bi), i.e. compute the reciprocal of
    a complex number.
    """
    cdef MPopts wopts
    cdef MPF t, u, m
    wopts = opts
    wopts.prec += 10
    MPF_init(&t)
    MPF_init(&u)
    MPF_init(&m)
    MPF_mul(&t, a, a,opts_exact)
    MPF_mul(&u, b, b,opts_exact)
    MPF_add(&m, &t, &u,wopts)
    MPF_div(&t, a, &m, opts)
    MPF_div(&u, b, &m, opts)
    MPF_set(re, &t)
    MPF_neg(im, &u)
    MPF_clear(&t)
    MPF_clear(&u)
    MPF_clear(&m)


cdef MPF_complex_pow_int(MPF *zre, MPF *zim, MPF *xre, MPF *xim, mpz_t n, MPopts opts):
    """
    Set zre+zim*i = (xre+xim) ^ n, i.e. raise a complex number to an integer power.
    """
    cdef MPopts wopts
    cdef long m

    if xim.special == S_ZERO:
        MPF_pow_int(zre, xre, n, opts)
        MPF_set_zero(zim)
        return

    if xre.special == S_ZERO:
        # n % 4
        m = mpz_get_si(n) % 4
        if m == 0:
            MPF_pow_int(zre, xim, n, opts)
            MPF_set_zero(zim)
            return
        if m == 1:
            MPF_set_zero(zre)
            MPF_pow_int(zim, xim, n, opts)
            return
        if m == 2:
            MPF_pow_int(zre, xim, n, opts)
            MPF_neg(zre, zre)
            MPF_set_zero(zim)
            return
        if m == 3:
            MPF_set_zero(zre)
            MPF_pow_int(zim, xim, n, opts)
            MPF_neg(zim, zim)
            return

    if mpz_reasonable_shift(n):
        m = mpz_get_si(n)
        if m == 0:
            MPF_set_one(zre)
            MPF_set_zero(zim)
            return
        if m == 1:
            MPF_pos(zre, xre, opts)
            MPF_pos(zim, xim, opts)
            return
        if m == 2:
            MPF_complex_square(zre, zim, xre, xim, opts)
            return
        if m == -1:
            MPF_complex_reciprocal(zre, zim, xre, xim, opts)
            return
        if m == -2:
            wopts = opts
            wopts.prec += 10
            MPF_complex_square(zre, zim, xre, xim, wopts)
            MPF_complex_reciprocal(zre, zim, zre, zim, opts)
            return

    xret = MPF_to_tuple(xre)
    ximt = MPF_to_tuple(xim)
    from mpmath.libmp import mpc_pow_int
    vr, vi = mpc_pow_int((xret, ximt), mpzi(n),
                         opts.prec, rndmode_to_python(opts.rounding))
    MPF_set_tuple(zre, vr)
    MPF_set_tuple(zim, vi)


cdef MPF_complex_pow_re(MPF *zre, MPF *zim, MPF *xre, MPF *xim, MPF *y, MPopts opts):
    """
    Set (zre+zim*i) = (xre+xim*i) ^ y, i.e. raise a complex number
    to a real power.
    """

    cdef mpz_t tm
    cdef MPopts wopts

    if y.special == S_ZERO:
        if xre.special == S_NORMAL and xim.special == S_NORMAL:
            # x ^ 0
            MPF_set_one(zre)
            MPF_set_zero(zim)
            return

    wopts = opts
    wopts.prec += 10

    if y.special == S_NORMAL:
        # Integer
        if mpz_cmp_ui(y.exp, 0) >= 0 and mpz_reasonable_shift(y.exp):
            mpz_init_set(tm, y.man)
            mpz_mul_2exp(tm, tm, mpz_get_ui(y.exp))
            MPF_complex_pow_int(zre, zim, xre, xim, tm, opts)
            mpz_clear(tm)
            return
        # x ^ (n/2)
        if mpz_cmp_si(y.exp, -1) == 0:
            mpz_init_set(tm, y.man)
            MPF_complex_sqrt(zre, zim, xre, xim, wopts)
            MPF_complex_pow_int(zre, zim, zre, zim, tm, opts)
            mpz_clear(tm)
            return

    xret = MPF_to_tuple(xre)
    ximt = MPF_to_tuple(xim)
    yret = MPF_to_tuple(y)
    from mpmath.libmp import mpc_pow_mpf
    vr, vi = mpc_pow_mpf((xret, ximt), yret,
                         opts.prec, rndmode_to_python(opts.rounding))
    MPF_set_tuple(zre, vr)
    MPF_set_tuple(zim, vi)


cdef MPF_complex_pow(MPF *zre, MPF *zim, MPF *xre, MPF *xim, MPF *yre, MPF *yim, MPopts opts):
    """
    Set (zre + zim*i) = (xre+xim*i) ^ (yre+yim*i).
    """
    if yim.special == S_ZERO:
        MPF_complex_pow_re(zre, zim, xre, xim, yre, opts)
        return
    xret = MPF_to_tuple(xre)
    ximt = MPF_to_tuple(xim)
    yret = MPF_to_tuple(yre)
    yimt = MPF_to_tuple(yim)
    from mpmath.libmp import mpc_pow
    vr, vi = mpc_pow((xret, ximt), (yret, yimt),
                     opts.prec, rndmode_to_python(opts.rounding))
    MPF_set_tuple(zre, vr)
    MPF_set_tuple(zim, vi)


cdef mpz_set_tuple_fixed(mpz_t x, tuple t, long prec):
    """
    Set the integer x to a fixed-point number with specified precision
    and the value of t = (sign,man,exp,bc). Truncating division is used
    if the value cannot be represented exactly.
    """
    cdef long offset
    sign, man, exp, bc = t
    mpz_set_integer(x, man)
    if sign:
        mpz_neg(x, x)
    offset = exp + prec
    if offset >= 0:
        mpz_mul_2exp(x, x, offset)
    else:
        mpz_tdiv_q_2exp(x, x, -offset)

cdef mpz_set_complex_tuple_fixed(mpz_t x, mpz_t y, tuple t, long prec):
    """
    Set the integers (x,y) to fixed-point numbers with the values of
    the mpf pair t = ((xsign,xman,xexp,xbc), (ysign,yman,yexp,ybc)).
    """
    mpz_set_tuple_fixed(x, t[0], prec)
    mpz_set_tuple_fixed(y, t[1], prec)

cdef MPF_set_fixed(MPF *x, mpz_t man, long wp, long prec, int rnd):
    """
    Set value of an MPF given a fixed-point mantissa of precision wp,
    rounding to the given precision and rounding mode.
    """
    cdef MPopts opts
    opts.prec = prec
    opts.rounding = rnd
    x.special = S_NORMAL
    mpz_set(x.man, man)
    mpz_set_si(x.exp, -wp)
    MPF_normalize(x, opts)

# TODO: we should allocate these dynamically
DEF MAX_PARAMS = 128
cdef mpz_t AINT[MAX_PARAMS]
cdef mpz_t BINT[MAX_PARAMS]
cdef mpz_t AP[MAX_PARAMS]
cdef mpz_t AQ[MAX_PARAMS]
cdef mpz_t BP[MAX_PARAMS]
cdef mpz_t BQ[MAX_PARAMS]
cdef mpz_t AREAL[MAX_PARAMS]
cdef mpz_t BREAL[MAX_PARAMS]
cdef mpz_t ACRE[MAX_PARAMS]
cdef mpz_t ACIM[MAX_PARAMS]
cdef mpz_t BCRE[MAX_PARAMS]
cdef mpz_t BCIM[MAX_PARAMS]


cdef MPF_hypsum(MPF *a, MPF *b, int p, int q, param_types, str ztype, coeffs, z,
    long prec, long wp, long epsshift, dict magnitude_check, kwargs):
    """
    Evaluates a+bi = pFq(..., z) by summing the hypergeometric
    series in fixed-point arithmetic.

    This basically a Cython version of
    mpmath.libmp.libhyper.make_hyp_summator(). It should produce identical
    results (the calculations are exactly the same, and the same rounding
    is used for divisions).

    This function is not intended to be called directly; it is wrapped
    by the hypsum_internal function in ext_main.pyx.
    """
    cdef long i, j, k, n, p_mag, cancellable_real, MAX, magn
    cdef int have_complex_param, have_complex_arg, have_complex

    cdef mpz_t SRE, SIM, PRE, PIM, ZRE, ZIM, TRE, TIM, URE, UIM, MUL, DIV, HIGH, LOW, one
    # Count number of parameters
    cdef int aint, bint, arat, brat, areal, breal, acomplex, bcomplex
    cdef int have_multiplier

    if p >= MAX_PARAMS or q >= MAX_PARAMS:
        raise NotImplementedError("too many hypergeometric function parameters")

    have_complex_param = 'C' in param_types
    have_complex_arg = ztype == 'C'
    have_complex = have_complex_param or have_complex_arg

    mpz_init(one)
    mpz_init(SRE)
    mpz_init(SIM)
    mpz_init(PRE)
    mpz_init(PIM)
    mpz_init(ZRE)
    mpz_init(ZIM)
    mpz_init(MUL)
    mpz_init(DIV)
    mpz_init(HIGH)
    mpz_init(LOW)
    mpz_init(TRE)
    mpz_init(TIM)
    mpz_init(URE)
    mpz_init(UIM)

    aint = bint = arat = brat = areal = breal = acomplex = bcomplex = 0

    MAX = kwargs.get('maxterms', wp*100)
    mpz_set_ui(HIGH, 1)
    mpz_mul_2exp(HIGH, HIGH, epsshift)
    mpz_neg(LOW, HIGH)

    mpz_set_ui(one, 1)
    mpz_mul_2exp(one, one, wp)
    mpz_set(SRE, one)
    mpz_set(PRE, one)

    # Copy input data to mpzs
    if have_complex_arg:
        mpz_set_complex_tuple_fixed(ZRE, ZIM, z, wp)
    else:
        mpz_set_tuple_fixed(ZRE, z, wp)
    for i in range(p):
        sig_check()
        if param_types[i] == 'Z':
            mpz_init(AINT[aint])
            mpz_set_integer(AINT[aint], coeffs[i])
            aint += 1
        elif param_types[i] == 'Q':
            mpz_init(AP[arat])
            mpz_init(AQ[arat])
            __p, __q = coeffs[i]._mpq_
            mpz_set_integer(AP[arat], __p)
            mpz_set_integer(AQ[arat], __q)
            arat += 1
        elif param_types[i] == 'R':
            mpz_init(AREAL[areal])
            mpz_set_tuple_fixed(AREAL[areal], coeffs[i]._mpf_, wp)
            areal += 1
        elif param_types[i] == 'C':
            mpz_init(ACRE[acomplex])
            mpz_init(ACIM[acomplex])
            mpz_set_complex_tuple_fixed(ACRE[acomplex], ACIM[acomplex], coeffs[i]._mpc_, wp)
            acomplex += 1
        else:
            raise ValueError
    for i in range(p,p+q):
        sig_check()
        if param_types[i] == 'Z':
            mpz_init(BINT[bint])
            mpz_set_integer(BINT[bint], coeffs[i])
            bint += 1
        elif param_types[i] == 'Q':
            mpz_init(BP[brat])
            mpz_init(BQ[brat])
            __p, __q = coeffs[i]._mpq_
            mpz_set_integer(BP[brat], __p)
            mpz_set_integer(BQ[brat], __q)
            brat += 1
        elif param_types[i] == 'R':
            mpz_init(BREAL[breal])
            mpz_set_tuple_fixed(BREAL[breal], coeffs[i]._mpf_, wp)
            breal += 1
        elif param_types[i] == 'C':
            mpz_init(BCRE[bcomplex])
            mpz_init(BCIM[bcomplex])
            mpz_set_complex_tuple_fixed(BCRE[bcomplex], BCIM[bcomplex], coeffs[i]._mpc_, wp)
            bcomplex += 1
        else:
            raise ValueError

    cancellable_real = min(areal, breal)

    # Main loop
    for n in range(1, 10**8):

        if n in magnitude_check:
            p_mag = mpz_bitcount(PRE)
            if have_complex:
                p_mag = max(p_mag, mpz_bitcount(PIM))
            magnitude_check[n] = wp - p_mag

        # Update rational part of product
        mpz_set_ui(MUL, 1)
        mpz_set_ui(DIV, n)

        for i in range(aint): mpz_mul(MUL, MUL, AINT[i])
        for i in range(arat): mpz_mul(MUL, MUL, AP[i])
        for i in range(brat): mpz_mul(MUL, MUL, BQ[i])
        for i in range(bint): mpz_mul(DIV, DIV, BINT[i])
        for i in range(brat): mpz_mul(DIV, DIV, BP[i])
        for i in range(arat): mpz_mul(DIV, DIV, AQ[i])

        # Check for singular terms
        if mpz_sgn(DIV) == 0:
            if mpz_sgn(MUL) == 0:
                break
            raise ZeroDivisionError

        # Multiply real factors
        for k in range(cancellable_real):
            sig_check()
            mpz_mul(PRE, PRE, AREAL[k])
            mpz_fdiv_q(PRE, PRE, BREAL[k])
        for k in range(cancellable_real, areal):
            sig_check()
            mpz_mul(PRE, PRE, AREAL[k])
            mpz_fdiv_q_2exp(PRE, PRE, wp)
        for k in range(cancellable_real, breal):
            sig_check()
            mpz_mul_2exp(PRE, PRE, wp)
            mpz_fdiv_q(PRE, PRE, BREAL[k])
        if have_complex:
            for k in range(cancellable_real):
                sig_check()
                mpz_mul(PIM, PIM, AREAL[k])
                mpz_fdiv_q(PIM, PIM, BREAL[k])
            for k in range(cancellable_real, areal):
                sig_check()
                mpz_mul(PIM, PIM, AREAL[k])
                mpz_fdiv_q_2exp(PIM, PIM, wp)
            for k in range(cancellable_real, breal):
                sig_check()
                mpz_mul_2exp(PIM, PIM, wp)
                mpz_fdiv_q(PIM, PIM, BREAL[k])

        # Update product
        if have_complex:
            if have_complex_arg:
                # PRE = ((mul*(PRE*ZRE-PIM*ZIM))//div)>>wp
                # PIM = ((mul*(PIM*ZRE+PRE*ZIM))//div)>>wp
                mpz_mul(TRE, PRE, ZRE)
                mpz_submul(TRE, PIM, ZIM)
                mpz_mul(TRE, TRE, MUL)

                mpz_mul(TIM, PIM, ZRE)
                mpz_addmul(TIM, PRE, ZIM)
                mpz_mul(TIM, TIM, MUL)

                mpz_fdiv_q(PRE, TRE, DIV)
                mpz_fdiv_q_2exp(PRE, PRE, wp)

                mpz_fdiv_q(PIM, TIM, DIV)
                mpz_fdiv_q_2exp(PIM, PIM, wp)
            else:
                mpz_mul(PRE, PRE, MUL)
                mpz_mul(PRE, PRE, ZRE)
                mpz_fdiv_q_2exp(PRE, PRE, wp)
                mpz_fdiv_q(PRE, PRE, DIV)

                mpz_mul(PIM, PIM, MUL)
                mpz_mul(PIM, PIM, ZRE)
                mpz_fdiv_q_2exp(PIM, PIM, wp)
                mpz_fdiv_q(PIM, PIM, DIV)

            for i in range(acomplex):
                sig_check()
                mpz_mul(TRE, PRE, ACRE[i])
                mpz_submul(TRE, PIM, ACIM[i])
                mpz_mul(TIM, PIM, ACRE[i])
                mpz_addmul(TIM, PRE, ACIM[i])
                mpz_fdiv_q_2exp(PRE, TRE, wp)
                mpz_fdiv_q_2exp(PIM, TIM, wp)

            for i in range(bcomplex):
                sig_check()
                mpz_mul(URE, BCRE[i], BCRE[i])
                mpz_addmul(URE, BCIM[i], BCIM[i])
                mpz_mul(TRE, PRE, BCRE[i])
                mpz_addmul(TRE, PIM, BCIM[i])
                mpz_mul(TIM, PIM, BCRE[i])
                mpz_submul(TIM, PRE, BCIM[i])
                mpz_mul_2exp(PRE, TRE, wp)
                mpz_fdiv_q(PRE, PRE, URE)
                mpz_mul_2exp(PIM, TIM, wp)
                mpz_fdiv_q(PIM, PIM, URE)
        else:
            mpz_mul(PRE, PRE, MUL)
            mpz_mul(PRE, PRE, ZRE)
            mpz_fdiv_q_2exp(PRE, PRE, wp)
            mpz_fdiv_q(PRE, PRE, DIV)

        # Add product to sum
        if have_complex:
            mpz_add(SRE, SRE, PRE)
            mpz_add(SIM, SIM, PIM)
            if mpz_cmpabs(PRE, HIGH) < 0 and mpz_cmpabs(PIM, HIGH) < 0:
                break
        else:
            mpz_add(SRE, SRE, PRE)
            if mpz_cmpabs(PRE, HIGH) < 0:
                break

        if n > MAX:
            from mpmath.libmp import NoConvergence
            raise NoConvergence('Hypergeometric series converges too slowly. Try increasing maxterms.')

        # +1 all parameters for next iteration
        for i in range(aint): mpz_add_ui(AINT[i], AINT[i], 1)
        for i in range(bint): mpz_add_ui(BINT[i], BINT[i], 1)
        for i in range(arat): mpz_add(AP[i], AP[i], AQ[i])
        for i in range(brat): mpz_add(BP[i], BP[i], BQ[i])
        for i in range(areal): mpz_add(AREAL[i], AREAL[i], one)
        for i in range(breal): mpz_add(BREAL[i], BREAL[i], one)
        for i in range(acomplex): mpz_add(ACRE[i], ACRE[i], one)
        for i in range(bcomplex): mpz_add(BCRE[i], BCRE[i], one)

    # Done
    if have_complex:
        MPF_set_fixed(a, SRE, wp, prec, ROUND_N)
        MPF_set_fixed(b, SIM, wp, prec, ROUND_N)
        if mpz_sgn(SRE):
            if mpz_sgn(SIM):
                magn = max(mpz_get_si(a.exp) + <long>mpz_bitcount(a.man),
                           mpz_get_si(b.exp) + <long>mpz_bitcount(b.man))
            else:
                magn = mpz_get_si(a.exp) + <long>mpz_bitcount(a.man)
        elif mpz_sgn(SIM):
            magn = mpz_get_si(b.exp) + <long>mpz_bitcount(b.man)
        else:
            magn = -wp
    else:
        MPF_set_fixed(a, SRE, wp, prec, ROUND_N)
        if mpz_sgn(SRE):
            magn = mpz_get_si(a.exp) + <long>mpz_bitcount(a.man)
        else:
            magn = -wp

    mpz_clear(one)
    mpz_clear(SRE)
    mpz_clear(SIM)
    mpz_clear(PRE)
    mpz_clear(PIM)
    mpz_clear(ZRE)
    mpz_clear(ZIM)
    mpz_clear(MUL)
    mpz_clear(DIV)
    mpz_clear(HIGH)
    mpz_clear(LOW)
    mpz_clear(TRE)
    mpz_clear(TIM)
    mpz_clear(URE)
    mpz_clear(UIM)

    for i in range(aint): mpz_clear(AINT[i])
    for i in range(bint): mpz_clear(BINT[i])
    for i in range(arat):
        mpz_clear(AP[i])
        mpz_clear(AQ[i])
    for i in range(brat):
        mpz_clear(BP[i])
        mpz_clear(BQ[i])
    for i in range(areal): mpz_clear(AREAL[i])
    for i in range(breal): mpz_clear(BREAL[i])
    for i in range(acomplex):
        mpz_clear(ACRE[i])
        mpz_clear(ACIM[i])
    for i in range(bcomplex):
        mpz_clear(BCRE[i])
        mpz_clear(BCIM[i])

    return have_complex, magn
