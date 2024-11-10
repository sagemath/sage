"""
mpmath floating-point numbers

Implements mpf and mpc types, with binary operations and support
for interaction with other types. Also implements the main
context class, and related utilities.
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

from cysignals.signals cimport sig_check

from sage.ext.stdsage cimport PY_NEW

from sage.libs.gmp.all cimport *
from sage.rings.integer cimport Integer

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

from sage.libs.mpmath.ext_impl cimport *

import mpmath.rational as rationallib
import mpmath.libmp as libmp
import mpmath.function_docs as function_docs
from mpmath.libmp import to_str
from mpmath.libmp import repr_dps, prec_to_dps, dps_to_prec

DEF OP_ADD = 0
DEF OP_SUB = 1
DEF OP_MUL = 2
DEF OP_DIV = 3
DEF OP_POW = 4
DEF OP_MOD = 5
DEF OP_RICHCMP = 6
DEF OP_EQ = (OP_RICHCMP + 2)
DEF OP_NE = (OP_RICHCMP + 3)
DEF OP_LT = (OP_RICHCMP + 0)
DEF OP_GT = (OP_RICHCMP + 4)
DEF OP_LE = (OP_RICHCMP + 1)
DEF OP_GE = (OP_RICHCMP + 5)

cdef MPopts opts_exact
cdef MPopts opts_double_precision
cdef MPopts opts_mini_prec

opts_exact.prec = 0
opts_exact.rounding = ROUND_N
opts_double_precision.prec = 53
opts_double_precision.rounding = ROUND_N
opts_mini_prec.prec = 5
opts_mini_prec.rounding = ROUND_D

cdef MPF MPF_C_0
cdef MPF MPF_C_1
cdef MPF MPF_C_2

MPF_init(&MPF_C_0)
MPF_set_zero(&MPF_C_0)
MPF_init(&MPF_C_1)
MPF_set_si(&MPF_C_1, 1)
MPF_init(&MPF_C_2)
MPF_set_si(&MPF_C_2, 2)

# Temporaries used for operands in binary operations
cdef mpz_t tmp_mpz
mpz_init(tmp_mpz)

cdef MPF tmp1
cdef MPF tmp2
cdef MPF tmp_opx_re
cdef MPF tmp_opx_im
cdef MPF tmp_opy_re
cdef MPF tmp_opy_im

MPF_init(&tmp1)
MPF_init(&tmp2)
MPF_init(&tmp_opx_re)
MPF_init(&tmp_opx_im)
MPF_init(&tmp_opy_re)
MPF_init(&tmp_opy_im)

cdef class Context
cdef class mpnumber
cdef class mpf_base
cdef class mpf
cdef class mpc
cdef class constant
cdef class wrapped_libmp_function
cdef class wrapped_specfun

cdef __isint(MPF *v):
    return v.special == S_ZERO or (v.special == S_NORMAL and mpz_sgn(v.exp) >= 0)

cdef int MPF_set_any(MPF *re, MPF *im, x, MPopts opts, bint str_tuple_ok) except -1:
    """
    Set re + im*i = x, where x is any Python number.

    Returns 0 if unable to coerce x; 1 if x is real and re was set;
    2 if x is complex and both re and im were set.

    If str_tuple_ok=True, strings and tuples are accepted and converted
    (useful for parsing arguments, but not for arithmetic operands).
    """
    if isinstance(x, mpf):
        MPF_set(re, &(<mpf>x).value)
        return 1
    if isinstance(x, mpc):
        MPF_set(re, &(<mpc>x).re)
        MPF_set(im, &(<mpc>x).im)
        return 2
    if isinstance(x, (int, Integer)):
        MPF_set_int(re, x)
        return 1
    if isinstance(x, float):
        MPF_set_double(re, x)
        return 1
    if isinstance(x, complex):
        MPF_set_double(re, x.real)
        MPF_set_double(im, x.imag)
        return 2
    if isinstance(x, constant):
        MPF_set_tuple(re, x.func(opts.prec, rndmode_to_python(opts.rounding)))
        return 1
    if hasattr(x, "_mpf_"):
        MPF_set_tuple(re, x._mpf_)
        return 1
    if hasattr(x, "_mpc_"):
        r, i = x._mpc_
        MPF_set_tuple(re, r)
        MPF_set_tuple(im, i)
        return 2
    if hasattr(x, "_mpmath_"):
        return MPF_set_any(re, im, x._mpmath_(opts.prec,
            rndmode_to_python(opts.rounding)), opts, False)
    if isinstance(x, rationallib.mpq):
        p, q = x._mpq_
        MPF_set_int(re, p)
        MPF_set_int(im, q)
        MPF_div(re, re, im, opts)
        #MPF_set_tuple(re, libmp.from_rational(p, q, opts.prec,
        #    rndmode_to_python(opts.rounding)))
        return 1
    if hasattr(x, '_mpi_'):
        a, b = x._mpi_
        if a == b:
            MPF_set_tuple(re, a)
            return 1
        raise ValueError("can only create mpf from zero-width interval")
    if str_tuple_ok:
        if isinstance(x, tuple):
            if len(x) == 2:
                MPF_set_man_exp(re, x[0], x[1])
                return 1
            elif len(x) == 4:
                MPF_set_tuple(re, x)
                return 1
        if isinstance(x, str):
            try:
                st = libmp.from_str(x, opts.prec,
                    rndmode_to_python(opts.rounding))
            except ValueError:
                return 0
            MPF_set_tuple(re, st)
            return 1
    return 0

cdef binop(int op, x, y, MPopts opts):
    cdef int typx
    cdef int typy
    cdef MPF xre, xim, yre, yim
    cdef mpf rr
    cdef mpc rc
    cdef MPopts altopts

    if isinstance(x, mpf):
        xre = (<mpf>x).value
        typx = 1
    elif isinstance(x, mpc):
        xre = (<mpc>x).re
        xim = (<mpc>x).im
        typx = 2
    else:
        typx = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, opts, False)
        if typx == 0:
            return NotImplemented
        xre = tmp_opx_re
        xim = tmp_opx_im

    if isinstance(y, mpf):
        yre = (<mpf>y).value
        typy = 1
    elif isinstance(y, mpc):
        yre = (<mpc>y).re
        yim = (<mpc>y).im
        typy = 2
    else:
        typy = MPF_set_any(&tmp_opy_re, &tmp_opy_im, y, opts, False)
        if typy == 0:
            return NotImplemented
        yre = tmp_opy_re
        yim = tmp_opy_im

    if op == OP_ADD:
        if typx == 1 and typy == 1:
            # Real result
            rr = mpf.__new__(mpf)
            MPF_add(&rr.value, &xre, &yre, opts)
            return rr
        else:
            # Complex result
            rc = mpc.__new__(mpc)
            MPF_add(&rc.re, &xre, &yre, opts)
            if typx == 1:
                MPF_set(&rc.im, &yim)
                #MPF_normalize(&rc.im, opts)
            elif typy == 1:
                MPF_set(&rc.im, &xim)
                #MPF_normalize(&rc.im, opts)
            else:
                MPF_add(&rc.im, &xim, &yim, opts)
            return rc

    elif op == OP_SUB:
        if typx == 1 and typy == 1:
            # Real result
            rr = mpf.__new__(mpf)
            MPF_sub(&rr.value, &xre, &yre, opts)
            return rr
        else:
            # Complex result
            rc = mpc.__new__(mpc)
            MPF_sub(&rc.re, &xre, &yre, opts)
            if typx == 1:
                MPF_neg(&rc.im, &yim)
                MPF_normalize(&rc.im, opts)
            elif typy == 1:
                MPF_set(&rc.im, &xim)
                MPF_normalize(&rc.im, opts)
            else:
                MPF_sub(&rc.im, &xim, &yim, opts)
            return rc

    elif op == OP_MUL:
        if typx == 1 and typy == 1:
            # Real result
            rr = mpf.__new__(mpf)
            MPF_mul(&rr.value, &xre, &yre, opts)
            return rr
        else:
            # Complex result
            rc = mpc.__new__(mpc)
            if typx == 1:
                MPF_mul(&rc.re, &yre, &xre, opts)
                MPF_mul(&rc.im, &yim, &xre, opts)
            elif typy == 1:
                MPF_mul(&rc.re, &xre, &yre, opts)
                MPF_mul(&rc.im, &xim, &yre, opts)
            else:
                # a*c - b*d
                MPF_mul(&rc.re, &xre,   &yre,  opts_exact)
                MPF_mul(&tmp1,  &xim,   &yim,  opts_exact)
                MPF_sub(&rc.re, &rc.re, &tmp1, opts)
                # a*d + b*c
                MPF_mul(&rc.im, &xre,   &yim,  opts_exact)
                MPF_mul(&tmp1,  &xim,   &yre,  opts_exact)
                MPF_add(&rc.im, &rc.im, &tmp1, opts)
            return rc

    elif op == OP_DIV:
        if typx == 1 and typy == 1:
            # Real result
            rr = mpf.__new__(mpf)
            MPF_div(&rr.value, &xre, &yre, opts)
            return rr
        else:
            rc = mpc.__new__(mpc)
            if typy == 1:
                MPF_div(&rc.re, &xre, &yre, opts)
                MPF_div(&rc.im, &xim, &yre, opts)
            else:
                if typx == 1:
                    xim = MPF_C_0
                altopts = opts
                altopts.prec += 10
                # m = c*c + d*d
                MPF_mul(&tmp1, &yre,  &yre, opts_exact)
                MPF_mul(&tmp2, &yim,  &yim, opts_exact)
                MPF_add(&tmp1, &tmp1, &tmp2, altopts)
                # (a*c+b*d)/m
                MPF_mul(&rc.re, &xre,   &yre, opts_exact)
                MPF_mul(&tmp2,  &xim,   &yim, opts_exact)
                MPF_add(&rc.re, &rc.re, &tmp2, altopts)
                MPF_div(&rc.re, &rc.re, &tmp1, opts)
                # (b*c-a*d)/m
                MPF_mul(&rc.im, &xim,   &yre, opts_exact)
                MPF_mul(&tmp2,  &xre,   &yim, opts_exact)
                MPF_sub(&rc.im, &rc.im, &tmp2, altopts)
                MPF_div(&rc.im, &rc.im, &tmp1, opts)
            return rc

    elif op == OP_POW:
        if typx == 1 and typy == 1:
            rr = mpf.__new__(mpf)
            if not MPF_pow(&rr.value, &xre, &yre, opts):
                return rr
        if typx == 1: xim = MPF_C_0
        if typy == 1: yim = MPF_C_0
        rc = mpc.__new__(mpc)
        MPF_complex_pow(&rc.re, &rc.im, &xre, &xim, &yre, &yim, opts)
        return rc

    elif op == OP_MOD:
        if typx != 1 or typx != 1:
            raise TypeError("mod for complex numbers")
        xret = MPF_to_tuple(&xre)
        yret = MPF_to_tuple(&yre)
        v = libmp.mpf_mod(xret, yret, opts.prec, rndmode_to_python(opts.rounding))
        rr = mpf.__new__(mpf)
        MPF_set_tuple(&rr.value, v)
        return rr

    elif op == OP_EQ:
        if typx == 1 and typy == 1:
            return MPF_eq(&xre, &yre)
        if typx == 1:
            return MPF_eq(&xre, &yre) and MPF_eq(&yim, &MPF_C_0)
        if typy == 1:
            return MPF_eq(&xre, &yre) and MPF_eq(&xim, &MPF_C_0)
        return MPF_eq(&xre, &yre) and MPF_eq(&xim, &yim)

    elif op == OP_NE:
        if typx == 1 and typy == 1:
            return MPF_ne(&xre, &yre)
        if typx == 1:
            return MPF_ne(&xre, &yre) or MPF_ne(&yim, &MPF_C_0)
        if typy == 1:
            return MPF_ne(&xre, &yre) or MPF_ne(&xim, &MPF_C_0)
        return MPF_ne(&xre, &yre) or MPF_ne(&xim, &yim)

    elif op == OP_LT:
        if typx != 1 or typy != 1:
            raise ValueError("cannot compare complex numbers")
        return MPF_lt(&xre, &yre)

    elif op == OP_GT:
        if typx != 1 or typy != 1:
            raise ValueError("cannot compare complex numbers")
        return MPF_gt(&xre, &yre)

    elif op == OP_LE:
        if typx != 1 or typy != 1:
            raise ValueError("cannot compare complex numbers")
        return MPF_le(&xre, &yre)

    elif op == OP_GE:
        if typx != 1 or typy != 1:
            raise ValueError("cannot compare complex numbers")
        return MPF_ge(&xre, &yre)

    return NotImplemented


cdef MPopts global_opts

global_context = None

cdef class Context:
    cdef public mpf, mpc, constant  # , def_mp_function
    cdef public trap_complex
    cdef public pretty

    def __cinit__(ctx):
        """
        At present, only a single global context should exist::

            sage: from mpmath import mp
            sage: type(mp)
            <class 'mpmath.ctx_mp.MPContext'>
        """
        global global_opts, global_context
        global_opts = opts_double_precision
        global_context = ctx
        ctx.mpf = mpf
        ctx.mpc = mpc
        ctx.constant = constant
        #ctx.def_mp_function = def_mp_function
        ctx._mpq = rationallib.mpq

    def default(ctx):
        """
        Set defaults.

        TESTS::

            sage: import mpmath
            sage: mpmath.mp.prec = 100
            sage: mpmath.mp.default()
            sage: mpmath.mp.prec
            53
        """
        global global_opts
        global_opts = opts_double_precision
        ctx.trap_complex = False
        ctx.pretty = False

    def _get_prec(ctx):
        """
        Controls the working precision in bits::

            sage: from mpmath import mp
            sage: mp.prec = 100
            sage: mp.prec
            100
            sage: mp.dps
            29
            sage: mp.prec = 53
        """
        return global_opts.prec

    def _set_prec(ctx, prec):
        """
        Controls the working precision in bits::

            sage: from mpmath import mp
            sage: mp.prec = 100
            sage: mp.prec
            100
            sage: mp.dps
            29
            sage: mp.prec = 53
        """
        global_opts.prec = prec

    def _set_dps(ctx, n):
        """
        Controls the working precision in decimal digits::

            sage: from mpmath import mp
            sage: mp.dps = 100
            sage: mp.prec
            336
            sage: mp.dps
            100
            sage: mp.prec = 53
        """
        global_opts.prec = dps_to_prec(int(n))

    def _get_dps(ctx):
        """
        Controls the working precision in decimal digits::

            sage: from mpmath import mp
            sage: mp.dps = 100
            sage: mp.prec
            336
            sage: mp.dps
            100
            sage: mp.prec = 53
        """
        return prec_to_dps(global_opts.prec)

    dps = property(_get_dps, _set_dps, doc=_get_dps.__doc__)
    prec = property(_get_prec, _set_prec, doc=_get_dps.__doc__)
    _dps = property(_get_dps, _set_dps, doc=_get_dps.__doc__)
    _prec = property(_get_prec, _set_prec, doc=_get_dps.__doc__)

    def _get_prec_rounding(ctx):
        """
        Return the precision and rounding mode::

            sage: from mpmath import mp
            sage: mp._get_prec_rounding()
            (53, 'n')
        """
        return global_opts.prec, rndmode_to_python(global_opts.rounding)

    _prec_rounding = property(_get_prec_rounding)

    cpdef mpf make_mpf(ctx, tuple v):
        """
        Create an mpf from tuple data::

            sage: import mpmath
            sage: float(mpmath.mp.make_mpf((0,1,-1,1)))
            0.5
        """
        cdef mpf x
        x = mpf.__new__(mpf)
        MPF_set_tuple(&x.value, v)
        return x

    cpdef mpc make_mpc(ctx, tuple v):
        """
        Create an mpc from tuple data::

            sage: import mpmath
            sage: complex(mpmath.mp.make_mpc(((0,1,-1,1), (1,1,-2,1))))
            (0.5-0.25j)
        """
        cdef mpc x
        x = mpc.__new__(mpc)
        MPF_set_tuple(&x.re, v[0])
        MPF_set_tuple(&x.im, v[1])
        return x

    def convert(ctx, x, strings=True):
        """
        Convert *x* to an ``mpf``, ``mpc`` or ``mpi``. If *x* is of type ``mpf``,
        ``mpc``, ``int``, ``float``, ``complex``, the conversion
        will be performed losslessly.

        If *x* is a string, the result will be rounded to the present
        working precision. Strings representing fractions or complex
        numbers are permitted.

        TESTS::

            sage: from mpmath import mp, convert
            sage: mp.dps = 15; mp.pretty = False
            sage: convert(3.5)
            mpf('3.5')
            sage: convert('2.1')
            mpf('2.1000000000000001')
            sage: convert('3/4')
            mpf('0.75')
            sage: convert('2+3j')
            mpc(real='2.0', imag='3.0')
        """
        cdef mpf rr
        cdef mpc rc
        if isinstance(x, mpnumber):
            return x
        typx = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, global_opts, strings)
        if typx == 1:
            rr = mpf.__new__(mpf)
            MPF_set(&rr.value, &tmp_opx_re)
            return rr
        if typx == 2:
            rc = mpc.__new__(mpc)
            MPF_set(&rc.re, &tmp_opx_re)
            MPF_set(&rc.im, &tmp_opx_im)
            return rc
        return ctx._convert_fallback(x, strings)

    def isnan(ctx, x):
        """
        For an ``mpf`` *x*, determines whether *x* is not-a-number (nan).

        TESTS::

            sage: from mpmath import isnan, nan
            sage: isnan(nan), isnan(3)
            (True, False)
        """
        cdef int s, t, typ
        if isinstance(x, mpf):
            return (<mpf>x).value.special == S_NAN
        if isinstance(x, mpc):
            s = (<mpc>x).re.special
            t = (<mpc>x).im.special
            return s == S_NAN or t == S_NAN
        if type(x) is int or type(x) is long or isinstance(x, Integer) \
            or isinstance(x, rationallib.mpq):
            return False
        typ = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, global_opts, 0)
        if typ == 1:
            s = tmp_opx_re.special
            return s == S_NAN
        if typ == 2:
            s = tmp_opx_re.special
            t = tmp_opx_im.special
            return s == S_NAN or t == S_NAN
        raise TypeError("isnan() needs a number as input")

    def isinf(ctx, x):
        """
        Return *True* if the absolute value of *x* is infinite;
        otherwise return *False*.

        TESTS::

            sage: from mpmath import isinf, inf, mpc
            sage: isinf(inf)
            True
            sage: isinf(-inf)
            True
            sage: isinf(3)
            False
            sage: isinf(3+4j)
            False
            sage: isinf(mpc(3,inf))
            True
            sage: isinf(mpc(inf,3))
            True
        """
        cdef int s, t, typ
        if isinstance(x, mpf):
            s = (<mpf>x).value.special
            return s == S_INF or s == S_NINF
        if isinstance(x, mpc):
            s = (<mpc>x).re.special
            t = (<mpc>x).im.special
            return s == S_INF or s == S_NINF or t == S_INF or t == S_NINF
        if type(x) is int or type(x) is long or isinstance(x, Integer) \
            or isinstance(x, rationallib.mpq):
            return False
        typ = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, global_opts, 0)
        if typ == 1:
            s = tmp_opx_re.special
            return s == S_INF or s == S_NINF
        if typ == 2:
            s = tmp_opx_re.special
            t = tmp_opx_im.special
            return s == S_INF or s == S_NINF or t == S_INF or t == S_NINF
        raise TypeError("isinf() needs a number as input")

    def isnormal(ctx, x):
        """
        Determine whether *x* is "normal" in the sense of floating-point
        representation; that is, return *False* if *x* is zero, an
        infinity or NaN; otherwise return *True*. By extension, a
        complex number *x* is considered "normal" if its magnitude is
        normal.

        TESTS::

            sage: from mpmath import isnormal, inf, nan, mpc
            sage: isnormal(3)
            True
            sage: isnormal(0)
            False
            sage: isnormal(inf); isnormal(-inf); isnormal(nan)
            False
            False
            False
            sage: isnormal(0+0j)
            False
            sage: isnormal(0+3j)
            True
            sage: isnormal(mpc(2,nan))
            False
        """
        # TODO: optimize this
        if hasattr(x, "_mpf_"):
            return bool(x._mpf_[1])
        if hasattr(x, "_mpc_"):
            re, im = x._mpc_
            re_normal = bool(re[1])
            im_normal = bool(im[1])
            if re == libmp.fzero: return im_normal
            if im == libmp.fzero: return re_normal
            return re_normal and im_normal
        if type(x) is int or type(x) is long or isinstance(x, Integer) \
            or isinstance(x, rationallib.mpq):
            return bool(x)
        x = ctx.convert(x)
        if hasattr(x, '_mpf_') or hasattr(x, '_mpc_'):
            return ctx.isnormal(x)
        raise TypeError("isnormal() needs a number as input")

    def isint(ctx, x, gaussian=False):
        """
        Return *True* if *x* is integer-valued; otherwise return
        *False*.

        TESTS::

            sage: from mpmath import isint, mpf, inf
            sage: isint(3)
            True
            sage: isint(mpf(3))
            True
            sage: isint(3.2)
            False
            sage: isint(inf)
            False

        Optionally, Gaussian integers can be checked for::

            sage: isint(3+0j)
            True
            sage: isint(3+2j)
            False
            sage: isint(3+2j, gaussian=True)
            True
        """
        cdef MPF v
        cdef MPF w
        cdef int typ
        if type(x) is int or type(x) is long or isinstance(x, Integer):
            return True
        if isinstance(x, mpf):
            v = (<mpf>x).value
            return __isint(&v)
        if isinstance(x, mpc):
            v = (<mpc>x).re
            w = (<mpc>x).im
            if gaussian:
                return __isint(&v) and __isint(&w)
            return (w.special == S_ZERO) and __isint(&v)
        if isinstance(x, rationallib.mpq):
            p, q = x._mpq_
            return not (p % q)
        typ = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, global_opts, 0)
        if typ == 1:
            return __isint(&tmp_opx_re)
        if typ == 2:
            v = tmp_opx_re
            w = tmp_opx_im
            if gaussian:
                return __isint(&v) and __isint(&w)
            return (w.special == S_ZERO) and __isint(&v)
        raise TypeError("isint() needs a number as input")

    def fsum(ctx, terms, bint absolute=False, bint squared=False):
        """
        Calculate a sum containing a finite number of terms (for infinite
        series, see :func:`nsum`). The terms will be converted to
        mpmath numbers. For len(terms) > 2, this function is generally
        faster and produces more accurate results than the builtin
        Python function :func:`sum`.

        With ``squared=True`` each term is squared, and with ``absolute=True``
        the absolute value of each term is used.

        TESTS::

            sage: from mpmath import mp, fsum
            sage: mp.dps = 15; mp.pretty = False
            sage: fsum([1, 2, 0.5, 7])
            mpf('10.5')

        Check that the regression from `mpmath/issues/723 <https://github.com/mpmath/mpmath/issues/723>`__
        has been fixed::

            sage: from mpmath import *
            sage: mp.dps=16
            sage: zeta(-0.01 + 1000j)
            mpc(real='-8.9714595...', imag='8.7321793...')
        """
        cdef MPF sre, sim, tre, tim, tmp
        cdef mpf rr
        cdef mpc rc
        cdef MPopts workopts
        cdef int styp, ttyp
        workopts = global_opts
        workopts.prec = workopts.prec * 2 + 50
        workopts.rounding = ROUND_D
        unknown = global_context.zero
        try:
            sig_check()
            MPF_init(&sre)
            MPF_init(&sim)
            MPF_init(&tre)
            MPF_init(&tim)
            MPF_init(&tmp)
            styp = 1
            for term in terms:
                ttyp = MPF_set_any(&tre, &tim, term, workopts, 0)
                if ttyp == 0:
                    if absolute: term = ctx.absmax(term)
                    if squared: term = term**2
                    unknown += term
                    continue
                if absolute:
                    if squared:
                        if ttyp == 1:
                            MPF_mul(&tre, &tre, &tre, opts_exact)
                            MPF_add(&sre, &sre, &tre, workopts)
                        elif ttyp == 2:
                            # |(a+bi)^2| = a^2+b^2
                            MPF_mul(&tre, &tre, &tre, opts_exact)
                            MPF_add(&sre, &sre, &tre, workopts)
                            MPF_mul(&tim, &tim, &tim, opts_exact)
                            MPF_add(&sre, &sre, &tim, workopts)
                    else:
                        if ttyp == 1:
                            MPF_abs(&tre, &tre)
                            MPF_add(&sre, &sre, &tre, workopts)
                        elif ttyp == 2:
                            # |a+bi| = sqrt(a^2+b^2)
                            MPF_mul(&tre, &tre, &tre, opts_exact)
                            MPF_mul(&tim, &tim, &tim, opts_exact)
                            MPF_add(&tre, &tre, &tim, workopts)
                            MPF_sqrt(&tre, &tre, workopts)
                            MPF_add(&sre, &sre, &tre, workopts)
                elif squared:
                    if ttyp == 1:
                        MPF_mul(&tre, &tre, &tre, opts_exact)
                        MPF_add(&sre, &sre, &tre, workopts)
                    elif ttyp == 2:
                        # (a+bi)^2 = a^2-b^2 + 2i*ab
                        MPF_mul(&tmp, &tre, &tim, opts_exact)
                        MPF_mul(&tmp, &tmp, &MPF_C_2, opts_exact)
                        MPF_add(&sim, &sim, &tmp, workopts)
                        MPF_mul(&tre, &tre, &tre, opts_exact)
                        MPF_add(&sre, &sre, &tre, workopts)
                        MPF_mul(&tim, &tim, &tim, opts_exact)
                        MPF_sub(&sre, &sre, &tim, workopts)
                        styp = 2
                else:
                    if ttyp == 1:
                        MPF_add(&sre, &sre, &tre, workopts)
                    elif ttyp == 2:
                        MPF_add(&sre, &sre, &tre, workopts)
                        MPF_add(&sim, &sim, &tim, workopts)
                        styp = 2
            MPF_clear(&tre)
            MPF_clear(&tim)
            if styp == 1:
                rr = mpf.__new__(mpf)
                MPF_set(&rr.value, &sre)
                MPF_clear(&sre)
                MPF_clear(&sim)
                MPF_normalize(&rr.value, global_opts)
                if unknown is not global_context.zero:
                    return ctx._stupid_add(rr, unknown)
                return rr
            elif styp == 2:
                rc = mpc.__new__(mpc)
                MPF_set(&rc.re, &sre)
                MPF_set(&rc.im, &sim)
                MPF_clear(&sre)
                MPF_clear(&sim)
                MPF_normalize(&rc.re, global_opts)
                MPF_normalize(&rc.im, global_opts)
                if unknown is not global_context.zero:
                    return ctx._stupid_add(rc, unknown)
                return rc
            else:
                MPF_clear(&sre)
                MPF_clear(&sim)
                return +unknown
        except KeyboardInterrupt:
            raise KeyboardInterrupt('Ctrl-C pressed while running fsum')

    def fdot(ctx, A, B=None, bint conjugate=False):
        r"""
        Compute the dot product of the iterables `A` and `B`.

        .. MATH::

            \sum_{k=0} A_k B_k.

        Alternatively, :func:`fdot` accepts a single iterable of pairs.
        In other words, ``fdot(A,B)`` and ``fdot(zip(A,B))`` are equivalent.

        The elements are automatically converted to mpmath numbers.

        TESTS::

            sage: from mpmath import mp, fdot
            sage: mp.dps = 15; mp.pretty = False
            sage: A = [2, 1.5r, 3]
            sage: B = [1, -1, 2]
            sage: fdot(A, B)
            mpf('6.5')
            sage: list(zip(A, B))
            [(2, 1), (1.5, -1), (3, 2)]
            sage: fdot(_)
            mpf('6.5')
        """
        if B:
            A = zip(A, B)
        cdef MPF sre, sim, tre, tim, ure, uim, tmp
        cdef mpf rr
        cdef mpc rc
        cdef MPopts workopts
        cdef int styp, ttyp, utyp
        MPF_init(&sre)
        MPF_init(&sim)
        MPF_init(&tre)
        MPF_init(&tim)
        MPF_init(&ure)
        MPF_init(&uim)
        MPF_init(&tmp)
        workopts = global_opts
        workopts.prec = workopts.prec * 2 + 50
        workopts.rounding = ROUND_D
        unknown = global_context.zero
        styp = 1
        for a, b in A:
            ttyp = MPF_set_any(&tre, &tim, a, workopts, 0)
            utyp = MPF_set_any(&ure, &uim, b, workopts, 0)
            if utyp == 2 and conjugate:
                MPF_neg(&uim, &uim)
            if ttyp == 0 or utyp == 0:
                if conjugate:
                    b = b.conj()
                unknown += a * b
                continue
            styp = max(styp, ttyp)
            styp = max(styp, utyp)
            if ttyp == 1:
                if utyp == 1:
                    MPF_mul(&tre, &tre, &ure, opts_exact)
                    MPF_add(&sre, &sre, &tre, workopts)
                elif utyp == 2:
                    MPF_mul(&ure, &ure, &tre, opts_exact)
                    MPF_mul(&uim, &uim, &tre, opts_exact)
                    MPF_add(&sre, &sre, &ure, workopts)
                    MPF_add(&sim, &sim, &uim, workopts)
                    styp = 2
            elif ttyp == 2:
                styp = 2
                if utyp == 1:
                    MPF_mul(&tre, &tre, &ure, opts_exact)
                    MPF_mul(&tim, &tim, &ure, opts_exact)
                    MPF_add(&sre, &sre, &tre, workopts)
                    MPF_add(&sim, &sim, &tim, workopts)
                elif utyp == 2:
                    MPF_mul(&tmp, &tre, &ure, opts_exact)
                    MPF_add(&sre, &sre, &tmp, workopts)
                    MPF_mul(&tmp, &tim, &uim, opts_exact)
                    MPF_sub(&sre, &sre, &tmp, workopts)
                    MPF_mul(&tmp, &tim, &ure, opts_exact)
                    MPF_add(&sim, &sim, &tmp, workopts)
                    MPF_mul(&tmp, &tre, &uim, opts_exact)
                    MPF_add(&sim, &sim, &tmp, workopts)
        MPF_clear(&tre)
        MPF_clear(&tim)
        MPF_clear(&ure)
        MPF_clear(&uim)
        MPF_clear(&tmp)
        if styp == 1:
            rr = mpf.__new__(mpf)
            MPF_set(&rr.value, &sre)
            MPF_clear(&sre)
            MPF_clear(&sim)
            MPF_normalize(&rr.value, global_opts)
            if unknown is not global_context.zero:
                return ctx._stupid_add(rr, unknown)
            return rr
        elif styp == 2:
            rc = mpc.__new__(mpc)
            MPF_set(&rc.re, &sre)
            MPF_set(&rc.im, &sim)
            MPF_clear(&sre)
            MPF_clear(&sim)
            MPF_normalize(&rc.re, global_opts)
            MPF_normalize(&rc.im, global_opts)
            if unknown is not global_context.zero:
                return ctx._stupid_add(rc, unknown)
            return rc
        else:
            MPF_clear(&sre)
            MPF_clear(&sim)
            return +unknown

    # Doing a+b directly doesn't work with mpi, presumably due to
    # Cython trying to be clever with the operation resolution
    cdef _stupid_add(ctx, a, b):
        return a + b

    def _convert_param(ctx, x):
        """
        Internal function for parsing a hypergeometric function parameter.
        Retrurns (T, x) where T = 'Z', 'Q', 'R', 'C' depending on the
        type of the parameter, and with x converted to the canonical
        mpmath type.

        TESTS::

            sage: from mpmath import mp
            sage: (x, T) = mp._convert_param(3)
            sage: (x, type(x).__name__, T)
            (3, 'int', 'Z')
            sage: (x, T) = mp._convert_param(2.5)
            sage: (x, type(x).__name__, T)
            (mpq(5,2), 'mpq', 'Q')
            sage: (x, T) = mp._convert_param(2.3)
            sage: (str(x), type(x).__name__, T)
            ('2.3', 'mpf', 'R')
            sage: (x, T) = mp._convert_param(2+3j)
            sage: (x, type(x).__name__, T)
            (mpc(real='2.0', imag='3.0'), 'mpc', 'C')
            sage: mp.pretty = False
        """
        cdef MPF v
        cdef bint ismpf, ismpc
        if isinstance(x, (int, Integer)):
            return int(x), 'Z'
        if isinstance(x, tuple):
            p, q = x
            p = int(p)
            q = int(q)
            if not p % q:
                return p // q, 'Z'
            return rationallib.mpq((p,q)), 'Q'
        if isinstance(x, str) and '/' in x:
            p, q = x.split('/')
            p = int(p)
            q = int(q)
            if not p % q:
                return p // q, 'Z'
            return rationallib.mpq((p,q)), 'Q'
        if isinstance(x, constant):
            return x, 'R'
        ismpf = isinstance(x, mpf)
        ismpc = isinstance(x, mpc)
        if not (ismpf or ismpc):
            x = global_context.convert(x)
            ismpf = isinstance(x, mpf)
            ismpc = isinstance(x, mpc)
            if not (ismpf or ismpc):
                return x, 'U'
        if ismpf:
            v = (<mpf>x).value
        elif ismpc:
            if (<mpc>x).im.special != S_ZERO:
                return x, 'C'
            x = (<mpc>x).real
            v = (<mpc>x).re
        # A real number
        if v.special == S_ZERO:
            return 0, 'Z'
        if v.special != S_NORMAL:
            return x, 'U'
        if mpz_sgn(v.exp) >= 0:
            return (mpzi(v.man) << mpzi(v.exp)), 'Z'
        if mpz_cmp_si(v.exp, -4) > 0:
            p = mpzi(v.man)
            q = 1 << (-mpzi(v.exp))
            return rationallib.mpq((p,q)), 'Q'
        return x, 'R'

    def mag(ctx, x):
        """
        Quick logarithmic magnitude estimate of a number. Returns an
        integer or infinity `m` such that `|x| <= 2^m`. It is not
        guaranteed that `m` is an optimal bound, but it will never
        be too large by more than 2 (and probably not more than 1).

        TESTS::

            sage: from mpmath import *
            sage: mag(10), mag(10.0), mag(mpf(10)), int(ceil(log(10,2)))
            (4, 4, 4, 4)
            sage: mag(10j), mag(10+10j)
            (4, 5)
            sage: mag(0.01), int(ceil(log(0.01,2)))
            (-6, -6)
            sage: mag(0), mag(inf), mag(-inf), mag(nan)
            (mpf('-inf'), mpf('+inf'), mpf('+inf'), mpf('nan'))

    ::

            sage: class MyInt(int):
            ....:     pass
            sage: class MyFloat(float):
            ....:     pass
            sage: mag(MyInt(10))
            4
        """
        cdef int typ
        if isinstance(x, (int, Integer)):
            mpz_set_integer(tmp_opx_re.man, x)
            if mpz_sgn(tmp_opx_re.man) == 0:
                return global_context.ninf
            else:
                return mpz_sizeinbase(tmp_opx_re.man,2)
        if isinstance(x, rationallib.mpq):
            p, q = x._mpq_
            mpz_set_integer(tmp_opx_re.man, int(p))
            if mpz_sgn(tmp_opx_re.man) == 0:
                return global_context.ninf
            mpz_set_integer(tmp_opx_re.exp, int(q))
            return 1 + mpz_sizeinbase(tmp_opx_re.man,2) + mpz_sizeinbase(tmp_opx_re.exp,2)
        typ = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, global_opts, False)
        if typ == 1:
            if tmp_opx_re.special == S_ZERO:
                return global_context.ninf
            if tmp_opx_re.special == S_INF or tmp_opx_re.special == S_NINF:
                return global_context.inf
            if tmp_opx_re.special != S_NORMAL:
                return global_context.nan
            mpz_add_ui(tmp_opx_re.exp, tmp_opx_re.exp, mpz_sizeinbase(tmp_opx_re.man, 2))
            return mpzi(tmp_opx_re.exp)
        if typ == 2:
            if tmp_opx_re.special == S_NAN or tmp_opx_im.special == S_NAN:
                return global_context.nan
            if tmp_opx_re.special == S_INF or tmp_opx_im.special == S_NINF or \
               tmp_opx_im.special == S_INF or tmp_opx_im.special == S_NINF:
                return global_context.inf
            if tmp_opx_re.special == S_ZERO:
                if tmp_opx_im.special == S_ZERO:
                    return global_context.ninf
                else:
                    mpz_add_ui(tmp_opx_im.exp, tmp_opx_im.exp, mpz_sizeinbase(tmp_opx_im.man, 2))
                    return mpzi(tmp_opx_im.exp)
            elif tmp_opx_im.special == S_ZERO:
                mpz_add_ui(tmp_opx_re.exp, tmp_opx_re.exp, mpz_sizeinbase(tmp_opx_re.man, 2))
                return mpzi(tmp_opx_re.exp)
            mpz_add_ui(tmp_opx_im.exp, tmp_opx_im.exp, mpz_sizeinbase(tmp_opx_im.man, 2))
            mpz_add_ui(tmp_opx_re.exp, tmp_opx_re.exp, mpz_sizeinbase(tmp_opx_re.man, 2))
            if mpz_cmp(tmp_opx_re.exp, tmp_opx_im.exp) >= 0:
                mpz_add_ui(tmp_opx_re.exp, tmp_opx_re.exp, 1)
                return mpzi(tmp_opx_re.exp)
            else:
                mpz_add_ui(tmp_opx_im.exp, tmp_opx_im.exp, 1)
                return mpzi(tmp_opx_im.exp)
        raise TypeError("requires an mpf/mpc")

    def _wrap_libmp_function(ctx, mpf_f, mpc_f=None, mpi_f=None, doc="<no doc>"):
        """
        Create a high-level mpmath function from base functions working
        on mpf, mpc and mpi tuple data.

        TESTS::

            sage: from mpmath import mp
            sage: mp.pretty = False
            sage: f = lambda x, prec, rnd: x
            sage: g = mp._wrap_libmp_function(f)
            sage: g(mp.mpf(2))
            mpf('2.0')
        """
        name = mpf_f.__name__[4:]
        doc = function_docs.__dict__.get(name, "Computes the %s of x" % doc)
        # workaround lack of closures in Cython
        f_cls = type(name, (wrapped_libmp_function,), {'__doc__':doc})
        f = f_cls(mpf_f, mpc_f, mpi_f, doc)
        return f

    @classmethod
    def _wrap_specfun(cls, name, f, wrap):
        """
        Add the given function as a method to the context object,
        optionally wrapping it to do automatic conversions and
        allocating a small number of guard bits to suppress
        typical roundoff error.

        TESTS::

            sage: from mpmath import mp
            sage: mp._wrap_specfun("foo", lambda ctx, x: ctx.prec + x, True)
            sage: mp.pretty = False; mp.prec = 53
            sage: mp.foo(5)       # 53 + 10 guard bits + 5
            mpf('68.0')
        """
        doc = function_docs.__dict__.get(name, getattr(f, '__doc__', '<no doc>'))
        if wrap:
            # workaround lack of closures in Cython
            f_wrapped_cls = type(name, (wrapped_specfun,), {'__doc__':doc})
            f_wrapped = f_wrapped_cls(name, f)
        else:
            f_wrapped = f
        f_wrapped.__doc__ = doc
        setattr(cls, name, f_wrapped)

    cdef MPopts _fun_get_opts(ctx, kwargs) noexcept:
        """
        Helper function that extracts precision and rounding information
        from kwargs, or returns the global working precision and rounding
        if no options are specified.
        """
        cdef MPopts opts
        opts.prec = global_opts.prec
        opts.rounding = global_opts.rounding
        if kwargs:
            if 'prec' in kwargs:
                opts.prec = int(kwargs['prec'])
            if 'dps'  in kwargs:
                opts.prec = libmp.dps_to_prec(int(kwargs['dps']))
            if 'rounding' in kwargs:
                opts.rounding = rndmode_from_python(kwargs['rounding'])
        return opts

    def _sage_sqrt(ctx, x, **kwargs):
        """
        Square root of an mpmath number x, using the Sage mpmath backend.

        TESTS::

            sage: from mpmath import mp
            sage: mp.dps = 15
            sage: print(mp.sqrt(2))   # indirect doctest
            1.4142135623731
            sage: print(mp.sqrt(-2))
            (0.0 + 1.4142135623731j)
            sage: print(mp.sqrt(2+2j))
            (1.55377397403004 + 0.643594252905583j)
        """
        cdef MPopts opts
        cdef int typx
        cdef mpf rr
        cdef mpc rc
        cdef tuple rev, imv
        opts = ctx._fun_get_opts(kwargs)
        typx = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, opts, 1)
        if typx == 1:
            rr = mpf.__new__(mpf)
            if MPF_sqrt(&rr.value, &tmp_opx_re, opts):
                rc = mpc.__new__(mpc)
                MPF_complex_sqrt(&rc.re, &rc.im, &tmp_opx_re, &MPF_C_0, opts)
                return rc
            return rr
        elif typx == 2:
            rc = mpc.__new__(mpc)
            MPF_complex_sqrt(&rc.re, &rc.im, &tmp_opx_re, &tmp_opx_im, opts)
            return rc
        else:
            raise NotImplementedError("unknown argument")

    def _sage_exp(ctx, x, **kwargs):
        """
        Exponential function of an mpmath number x, using the Sage
        mpmath backend.

        EXAMPLES::

            sage: from mpmath import mp
            sage: mp.dps = 15
            sage: print(mp.exp(2))   # indirect doctest
            7.38905609893065
            sage: print(mp.exp(2+2j))
            (-3.07493232063936 + 6.71884969742825j)
        """
        cdef MPopts opts
        cdef int typx
        cdef mpf rr
        cdef mpc rc
        cdef tuple rev, imv
        opts = ctx._fun_get_opts(kwargs)
        typx = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, opts, 1)
        if typx == 1:
            rr = mpf.__new__(mpf)
            MPF_exp(&rr.value, &tmp_opx_re, opts)
            return rr
        elif typx == 2:
            rc = mpc.__new__(mpc)
            MPF_complex_exp(&rc.re, &rc.im, &tmp_opx_re, &tmp_opx_im, opts)
            return rc
        else:
            raise NotImplementedError("unknown argument")

    def _sage_cos(ctx, x, **kwargs):
        """
        Cosine of an mpmath number x, using the Sage mpmath backend.

        EXAMPLES::

            sage: from mpmath import mp
            sage: mp.dps = 15
            sage: print(mp.cos(2))   # indirect doctest
            -0.416146836547142
            sage: print(mp.cos(2+2j))
            (-1.56562583531574 - 3.29789483631124j)
        """
        cdef MPopts opts
        cdef int typx
        cdef mpf rr
        cdef mpc rc
        cdef tuple rev, imv
        opts = ctx._fun_get_opts(kwargs)
        typx = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, global_opts, 1)
        if typx == 1:
            rr = mpf.__new__(mpf)
            MPF_cos(&rr.value, &tmp_opx_re, opts)
            return rr
        elif typx == 2:
            rev = MPF_to_tuple(&tmp_opx_re)
            imv = MPF_to_tuple(&tmp_opx_im)
            cxu = libmp.mpc_cos((rev, imv), opts.prec,
                rndmode_to_python(opts.rounding))
            rc = mpc.__new__(mpc)
            MPF_set_tuple(&rc.re, cxu[0])
            MPF_set_tuple(&rc.im, cxu[1])
            return rc
        else:
            raise NotImplementedError("unknown argument")

    def _sage_sin(ctx, x, **kwargs):
        """
        Sine of an mpmath number x, using the Sage mpmath backend.

        EXAMPLES::

            sage: from mpmath import mp
            sage: mp.dps = 15
            sage: print(mp.sin(2))   # indirect doctest
            0.909297426825682
            sage: print(mp.sin(2+2j))
            (3.42095486111701 - 1.50930648532362j)
        """
        cdef MPopts opts
        cdef int typx
        cdef mpf rr
        cdef mpc rc
        cdef tuple rev, imv
        opts = ctx._fun_get_opts(kwargs)
        typx = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, global_opts, 1)
        if typx == 1:
            rr = mpf.__new__(mpf)
            MPF_sin(&rr.value, &tmp_opx_re, opts)
            return rr
        elif typx == 2:
            rev = MPF_to_tuple(&tmp_opx_re)
            imv = MPF_to_tuple(&tmp_opx_im)
            cxu = libmp.mpc_sin((rev, imv), opts.prec,
                rndmode_to_python(opts.rounding))
            rc = mpc.__new__(mpc)
            MPF_set_tuple(&rc.re, cxu[0])
            MPF_set_tuple(&rc.im, cxu[1])
            return rc
        else:
            raise NotImplementedError("unknown argument")

    def _sage_ln(ctx, x, **kwargs):
        """
        Natural logarithm of an mpmath number x, using the Sage mpmath backend.

        EXAMPLES::

            sage: from mpmath import mp
            sage: print(mp.ln(2))   # indirect doctest
            0.693147180559945
            sage: print(mp.ln(-2))
            (0.693147180559945 + 3.14159265358979j)
            sage: print(mp.ln(2+2j))
            (1.03972077083992 + 0.785398163397448j)
        """
        cdef MPopts opts
        cdef int typx
        cdef mpf rr
        cdef mpc rc
        cdef tuple rev, imv
        typx = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, global_opts, 1)
        prec = global_opts.prec
        rounding = rndmode_to_python(global_opts.rounding)
        opts.prec = global_opts.prec
        opts.rounding = global_opts.rounding
        if kwargs:
            if 'prec' in kwargs: opts.prec = int(kwargs['prec'])
            if 'dps'  in kwargs: opts.prec = libmp.dps_to_prec(int(kwargs['dps']))
            if 'rounding' in kwargs: opts.rounding = rndmode_from_python(kwargs['rounding'])
        if typx == 1:
            if MPF_sgn(&tmp_opx_re) < 0:
                rc = mpc.__new__(mpc)
                MPF_log(&rc.re, &tmp_opx_re, opts)
                MPF_set_pi(&rc.im, opts)
                return rc
            else:
                rr = mpf.__new__(mpf)
                MPF_log(&rr.value, &tmp_opx_re, opts)
                return rr
        elif typx == 2:
            rev = MPF_to_tuple(&tmp_opx_re)
            imv = MPF_to_tuple(&tmp_opx_im)
            cxu = libmp.mpc_log((rev, imv), opts.prec, rndmode_to_python(opts.rounding))
            rc = mpc.__new__(mpc)
            MPF_set_tuple(&rc.re, cxu[0])
            MPF_set_tuple(&rc.im, cxu[1])
            return rc

            #rc = mpc.__new__(mpc)
            #MPF_complex_log(&rc.re, &rc.im, &tmp_opx_re, &tmp_opx_im, opts)
            #return rc
        else:
            raise NotImplementedError("unknown argument")


cdef class wrapped_libmp_function:

    cdef public mpf_f, mpc_f, mpi_f, name, __name__, __doc__

    def __init__(self, mpf_f, mpc_f=None, mpi_f=None, doc="<no doc>"):
        """
        Create a high-level mpmath function from base functions working
        on mpf, mpc and mpi tuple data.

        TESTS::

            sage: from sage.libs.mpmath.ext_main import wrapped_libmp_function
            sage: from mpmath import mp
            sage: from mpmath.libmp import mpf_exp, mpf_sqrt
            sage: f = lambda x, prec, rnd: mpf_exp(mpf_sqrt(x, prec, rnd), prec, rnd)
            sage: g = wrapped_libmp_function(f)
            sage: g(mp.mpf(3))
            mpf('5.6522336740340915')
            sage: mp.exp(mp.sqrt(3))
            mpf('5.6522336740340915')
            sage: g(mp.mpf(3), prec=10)
            mpf('5.65625')
            sage: g(mp.mpf(3), prec=10, rounding='u')
            mpf('5.65625')
            sage: g(mp.mpf(3), prec=10, rounding='d')
            mpf('5.640625')
        """
        self.mpf_f = mpf_f
        self.mpc_f = mpc_f
        self.mpi_f = mpi_f
        self.name = self.__name__ = mpf_f.__name__[4:]
        self.__doc__ = function_docs.__dict__.get(self.name, "Computes the %s of x" % doc)

    def __call__(self, x, **kwargs):
        """
        A wrapped mpmath library function performs automatic
        conversions and uses the default working precision
        unless overridden::

            sage: from mpmath import mp
            sage: mp.sinh(2)
            mpf('3.6268604078470186')
            sage: mp.sinh(2, prec=10)
            mpf('3.625')
            sage: mp.sinh(2, prec=10, rounding='d')
            mpf('3.625')
            sage: mp.sinh(2, prec=10, rounding='u')
            mpf('3.62890625')
        """
        cdef int typx
        cdef tuple rev, imv, reu, cxu
        cdef mpf rr
        cdef mpc rc
        prec = global_opts.prec
        rounding = rndmode_to_python(global_opts.rounding)
        if kwargs:
            if 'prec' in kwargs: prec = int(kwargs['prec'])
            if 'dps'  in kwargs: prec = libmp.dps_to_prec(int(kwargs['dps']))
            if 'rounding' in kwargs: rounding = kwargs['rounding']
        typx = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, global_opts, 1)
        if typx == 1:
            rev = MPF_to_tuple(&tmp_opx_re)
            try:
                reu = self.mpf_f(rev, prec, rounding)
                rr = mpf.__new__(mpf)
                MPF_set_tuple(&rr.value, reu)
                return rr
            except libmp.ComplexResult:
                if global_context.trap_complex:
                    raise
                cxu = self.mpc_f((rev, libmp.fzero), prec, rounding)
                rc = mpc.__new__(mpc)
                MPF_set_tuple(&rc.re, cxu[0])
                MPF_set_tuple(&rc.im, cxu[1])
                return rc
        if typx == 2:
            rev = MPF_to_tuple(&tmp_opx_re)
            imv = MPF_to_tuple(&tmp_opx_im)
            cxu = self.mpc_f((rev, imv), prec, rounding)
            rc = mpc.__new__(mpc)
            MPF_set_tuple(&rc.re, cxu[0])
            MPF_set_tuple(&rc.im, cxu[1])
            return rc
        x = global_context.convert(x)
        if hasattr(x, "_mpf_") or hasattr(x, "_mpc_"):
            return self(x, **kwargs)
        #if hasattr(x, "_mpi_"):
        #    if self.mpi_f:
        #        return global_context.make_mpi(self.mpi_f(x._mpi_, prec))
        raise NotImplementedError("%s of a %s" % (self.name, type(x)))


cdef class wrapped_specfun:
    cdef public f, name, __name__, __doc__

    def __init__(self, name, f):
        """
        Create an object holding a wrapped mpmath function
        along with metadata (name and documentation).

        TESTS::

            sage: import mpmath
            sage: from sage.libs.mpmath.ext_main import wrapped_specfun
            sage: f = wrapped_specfun("f", lambda ctx, x: x)
            sage: f.name
            'f'
        """
        self.name = self.__name__ = name
        self.f = f
        self.__doc__ = function_docs.__dict__.get(name, "<no doc>")

    def __call__(self, *args, **kwargs):
        """
        Call wrapped mpmath function. Arguments are automatically converted
        to mpmath number, and the internal working precision is increased
        by a few bits to suppress typical rounding errors::

            sage: from mpmath import mp
            sage: from sage.libs.mpmath.ext_main import wrapped_specfun
            sage: f = wrapped_specfun("f", lambda ctx, x: x + ctx.prec)
            sage: f("1")     # 53 + 10 guard bits + 1
            mpf('64.0')
        """
        cdef int origprec
        args = [global_context.convert(a) for a in args]
        origprec = global_opts.prec
        global_opts.prec += 10
        try:
            retval = self.f(global_context, *args, **kwargs)
        finally:
            global_opts.prec = origprec
        return +retval


cdef class mpnumber:

    def __richcmp__(self, other, int op):
        """
        Comparison of mpmath numbers. Compatible numerical types
        are automatically converted to mpmath numbers::

            sage: from mpmath import mpf, mpc
            sage: mpf(3) == mpc(3)
            True
            sage: mpf(3) == mpc(4)
            False
            sage: mpf(3) == float(3)
            True
            sage: mpf(3) < float(2.5)
            False
            sage: mpc(3) < mpc(4)
            Traceback (most recent call last):
                ...
            ValueError: cannot compare complex numbers
        """
        return binop(OP_RICHCMP+op, self, other, global_opts)

    def __add__(self, other):
        """
        Addition of mpmath numbers. Compatible numerical types
        are automatically converted to mpmath numbers::

            sage: from mpmath import mpf, mpc
            sage: mpf(3) + mpc(3)
            mpc(real='6.0', imag='0.0')
            sage: float(4) + mpf(3)
            mpf('7.0')
        """
        return binop(OP_ADD, self, other, global_opts)

    def __sub__(self, other):
        """
        Subtraction of mpmath numbers. Compatible numerical types
        are automatically converted to mpmath numbers::

            sage: from mpmath import mpf, mpc
            sage: mpf(5) - mpc(3)
            mpc(real='2.0', imag='0.0')
            sage: float(4) - mpf(3)
            mpf('1.0')
        """
        return binop(OP_SUB, self, other, global_opts)

    def __mul__(self, other):
        """
        Multiplication of mpmath numbers. Compatible numerical types
        are automatically converted to mpmath numbers::

            sage: from mpmath import mpf, mpc
            sage: mpf(5) * mpc(3)
            mpc(real='15.0', imag='0.0')
            sage: float(4) * mpf(3)
            mpf('12.0')
        """
        return binop(OP_MUL, self, other, global_opts)

    def __truediv__(self, other):
        """
        Division of mpmath numbers. Compatible numerical types
        are automatically converted to mpmath numbers::

            sage: from mpmath import mpf, mpc
            sage: mpf(10) / mpc(5)
            mpc(real='2.0', imag='0.0')
            sage: float(9) / mpf(3)
            mpf('3.0')
        """
        return binop(OP_DIV, self, other, global_opts)

    def __mod__(self, other):
        """
        Remainder of mpmath numbers. Compatible numerical types
        are automatically converted to mpmath numbers::

            sage: from mpmath import mpf
            sage: mpf(12) % float(7)
            mpf('5.0')
        """
        return binop(OP_MOD, self, other, global_opts)

    def __pow__(self, other, mod):
        """
        Exponentiation of mpmath numbers. Compatible numerical types
        are automatically converted to mpmath numbers::

            sage: from mpmath import mpf, mpc
            sage: mpf(10) ** mpc(3)
            mpc(real='1000.0', imag='0.0')
            sage: mpf(3) ** float(2)
            mpf('9.0')
        """
        if mod is not None:
            raise ValueError("three-argument pow not supported")
        return binop(OP_POW, self, other, global_opts)

    def ae(s, t, rel_eps=None, abs_eps=None):
        """
        Check if two numbers are approximately equal to within the specified
        tolerance (see mp.almosteq for documentation)::

            sage: from mpmath import mpf, mpc
            sage: mpf(3).ae(mpc(3,1e-10))
            False
            sage: mpf(3).ae(mpc(3,1e-10), rel_eps=1e-5)
            True
        """
        return global_context.almosteq(s, t, rel_eps, abs_eps)


cdef class mpf_base(mpnumber):

    # Shared methods for mpf, constant. However, somehow some methods
    # (hash?, __richcmp__?) are not inherited, so they have to
    # be defined multiple times. TODO: fix this.

    def __hash__(self):
        """
        Support hashing of derived classes::

            sage: from mpmath import mpf
            sage: from sage.libs.mpmath.ext_main import mpf_base
            sage: class X(mpf_base): _mpf_ = mpf(3.25)._mpf_
            sage: hash(X()) == hash(float(X()))
            True
        """
        return libmp.mpf_hash(self._mpf_)

    def __repr__(self):
        """
        Support repr() of derived classes::

            sage: from mpmath import mpf
            sage: from sage.libs.mpmath.ext_main import mpf_base
            sage: class X(mpf_base): _mpf_ = mpf(3.25)._mpf_
            sage: repr(X())
            "mpf('3.25')"
        """
        if global_context.pretty:
            return str(self)
        n = repr_dps(global_opts.prec)
        return "mpf('%s')" % to_str(self._mpf_, n)

    def __str__(self):
        """
        Support str() of derived classes::

            sage: from mpmath import mpf
            sage: from sage.libs.mpmath.ext_main import mpf_base
            sage: class X(mpf_base): _mpf_ = mpf(3.25)._mpf_
            sage: str(X())
            '3.25'
        """
        return to_str(self._mpf_, global_context._str_digits)

    @property
    def real(self):
        """
        Support real part of derived classes::

            sage: from mpmath import mpf
            sage: from sage.libs.mpmath.ext_main import mpf_base
            sage: class X(mpf_base): _mpf_ = mpf(3.25)._mpf_
            sage: X().real
            mpf('3.25')
        """
        return self

    @property
    def imag(self):
        """
        Support imaginary part of derived classes::

            sage: from mpmath import mpf
            sage: from sage.libs.mpmath.ext_main import mpf_base
            sage: class X(mpf_base): _mpf_ = mpf(3.25)._mpf_
            sage: X().imag
            mpf('0.0')
        """
        return global_context.zero

    def conjugate(self):
        """
        Support complex conjugate of derived classes::

            sage: from mpmath import mpf
            sage: from sage.libs.mpmath.ext_main import mpf_base
            sage: class X(mpf_base): _mpf_ = mpf(3.25)._mpf_
            sage: X().conjugate()
            mpf('3.25')
        """
        return self

    @property
    def man(self):
        """
        Support mantissa extraction of derived classes::

            sage: from mpmath import mpf
            sage: from sage.libs.mpmath.ext_main import mpf_base
            sage: class X(mpf_base): _mpf_ = mpf(3.25)._mpf_
            sage: X().man
            13
        """
        return self._mpf_[1]

    @property
    def exp(self):
        """
        Support exponent extraction of derived classes::

            sage: from mpmath import mpf
            sage: from sage.libs.mpmath.ext_main import mpf_base
            sage: class X(mpf_base): _mpf_ = mpf(3.25)._mpf_
            sage: X().exp
            -2
        """
        return self._mpf_[2]

    @property
    def bc(self):
        """
        Support bitcount extraction of derived classes::

            sage: from mpmath import mpf
            sage: from sage.libs.mpmath.ext_main import mpf_base
            sage: class X(mpf_base): _mpf_ = mpf(3.25)._mpf_
            sage: X().bc
            4
        """
        return self._mpf_[3]

    # XXX: optimize
    def __int__(self):
        """
        Support integer conversion for derived classes::

            sage: from mpmath import mpf
            sage: from sage.libs.mpmath.ext_main import mpf_base
            sage: class X(mpf_base): _mpf_ = mpf(3.25)._mpf_
            sage: int(X())
            3
        """
        return int(libmp.to_int(self._mpf_))

    def __float__(self):
        """
        Support float conversion for derived classes::

            sage: from mpmath import mpf
            sage: from sage.libs.mpmath.ext_main import mpf_base
            sage: class X(mpf_base): _mpf_ = mpf(3.25)._mpf_
            sage: float(X())
            3.25
        """
        return libmp.to_float(self._mpf_)

    def __complex__(self):
        """
        Support complex conversion for derived classes::

            sage: from mpmath import mpf
            sage: from sage.libs.mpmath.ext_main import mpf_base
            sage: class X(mpf_base): _mpf_ = mpf(3.25)._mpf_
            sage: complex(X())
            (3.25+0j)
        """
        return complex(float(self))

    def to_fixed(self, prec):
        """
        Support conversion to a fixed-point integer for derived classes::

            sage: from mpmath import mpf
            sage: from sage.libs.mpmath.ext_main import mpf_base
            sage: class X(mpf_base): _mpf_ = mpf(3.25)._mpf_
            sage: X().to_fixed(30)
            3489660928
        """
        return libmp.to_fixed(self._mpf_, prec)

    def __getstate__(self):
        return libmp.to_pickable(self._mpf_)

    def __setstate__(self, val):
        self._mpf_ = libmp.from_pickable(val)


cdef class mpf(mpf_base):
    """
    An mpf instance holds a real-valued floating-point number. mpf:s
    work analogously to Python floats, but support arbitrary-precision
    arithmetic.
    """

    cdef MPF value

    def __init__(self, x=0, **kwargs):
        """
        Create an mpf from a recognized type, optionally rounding
        to a precision and in a direction different from the default.

        TESTS::

            sage: from mpmath import mpf
            sage: mpf()
            mpf('0.0')
            sage: mpf(5)
            mpf('5.0')
            sage: mpf('inf')
            mpf('+inf')
            sage: mpf('nan')
            mpf('nan')
            sage: mpf(float(2.5))
            mpf('2.5')
            sage: mpf("2.5")
            mpf('2.5')
            sage: mpf("0.3")
            mpf('0.29999999999999999')
            sage: mpf("0.3", prec=10)
            mpf('0.2998046875')
            sage: mpf("0.3", prec=10, rounding='u')
            mpf('0.30029296875')
        """
        cdef MPopts opts
        opts = global_opts
        if kwargs:
            if 'prec' in kwargs: opts.prec = int(kwargs['prec'])
            if 'dps'  in kwargs: opts.prec = libmp.dps_to_prec(int(kwargs['dps']))
            if 'rounding' in kwargs: opts.rounding = rndmode_from_python(kwargs['rounding'])
        if MPF_set_any(&self.value, &self.value, x, opts, 1) != 1:
            raise TypeError

    def __reduce__(self):
        """
        Support pickling::

            sage: from mpmath import mpf
            sage: loads(dumps(mpf(0.5))) == mpf(0.5)
            True
        """
        return (mpf, (), self._mpf_)

    def _get_mpf(self):
        """
        Return internal representation of ``self`` as a tuple
        of (sign bit, mantissa, exponent, bitcount)::

            sage: from mpmath import mp
            sage: mp.mpf(-3)._mpf_
            (1, 3, 0, 2)
        """
        return MPF_to_tuple(&self.value)

    def _set_mpf(self, v):
        """
        Set tuple value of ``self`` (warning: unsafe)::

            sage: from mpmath import mp
            sage: x = mp.mpf(-3)
            sage: x._mpf_ = (1, 3, -1, 2)
            sage: x
            mpf('-1.5')
        """
        MPF_set_tuple(&self.value, v)

    _mpf_ = property(_get_mpf, _set_mpf, doc=_get_mpf.__doc__)

    def __bool__(self):
        """
        Return whether the number is nonzero::

            sage: from mpmath import mpf
            sage: bool(mpf(3.5))
            True
            sage: bool(mpf(0.0))
            False
        """
        return self.value.special != S_ZERO

    def __hash__(self):
        """
        Hash values are compatible with builtin Python floats
        when the precision is small enough::

            sage: from mpmath import mpf
            sage: hash(mpf(2.5)) == hash(float(2.5))
            True
            sage: hash(mpf('inf')) == hash(float(Infinity))
            True
        """
        return libmp.mpf_hash(self._mpf_)

    @property
    def real(self):
        """
        Real part, leaves ``self`` unchanged::

            sage: from mpmath import mpf
            sage: mpf(2.5).real
            mpf('2.5')
        """
        return self

    @property
    def imag(self):
        """
        Imaginary part, equal to zero::

            sage: from mpmath import mpf
            sage: mpf(2.5).imag
            mpf('0.0')
        """
        return global_context.zero

    def conjugate(self):
        """
        Complex conjugate, leaves ``self`` unchanged::

            sage: from mpmath import mpf
            sage: mpf(2.5).conjugate()
            mpf('2.5')
        """
        return self

    @property
    def man(self):
        """
        Return the binary mantissa of ``self``. The result is a Sage
        integer::

            sage: from mpmath import mpf
            sage: mpf(-500.5).man
            1001
            sage: type(_)
            <class 'sage.rings.integer.Integer'>
        """
        return self._mpf_[1]

    @property
    def exp(self):
        """
        Return the binary exponent of ``self``::

            sage: from mpmath import mpf
            sage: mpf(1/64.).exp
            -6
        """
        return self._mpf_[2]

    @property
    def bc(self):
        """
        Return the number of bits in the mantissa of ``self``::

            sage: from mpmath import mpf
            sage: mpf(-256).bc
            1
            sage: mpf(-255).bc
            8
        """
        return self._mpf_[3]

    def to_fixed(self, long prec):
        """
        Convert to a fixed-point integer of the given precision::

            sage: from mpmath import mpf
            sage: mpf(7.25).to_fixed(30)
            7784628224
            sage: ZZ(7.25 * 2**30)
            7784628224
        """
        # return libmp.to_fixed(self._mpf_, prec)
        MPF_to_fixed(tmp_mpz, &self.value, prec, False)
        cdef Integer r
        r = PY_NEW(Integer)
        mpz_set(r.value, tmp_mpz)
        return r

    def __int__(self):
        """
        Convert to a Python integer (truncating if necessary)::

            sage: from mpmath import mpf
            sage: int(mpf(2.5))
            2
            sage: type(_)
            <... 'int'>
        """
        MPF_to_fixed(tmp_mpz, &self.value, 0, True)
        return mpzi(tmp_mpz)

    def __float__(self):
        """
        Convert to a double-precision Python float::

            sage: from mpmath import mpf
            sage: float(mpf(2.5))
            2.5
            sage: type(_)
            <... 'float'>
        """
        return MPF_to_double(&self.value, False)

    def __getstate__(self):
        """
        Support pickling::

            sage: from mpmath import mpf
            sage: loads(dumps(mpf(3))) == mpf(3)
            True
        """
        return libmp.to_pickable(self._mpf_)

    def __setstate__(self, val):
        """
        Support pickling::

            sage: from mpmath import mpf
            sage: loads(dumps(mpf(3))) == mpf(3)
            True
        """
        self._mpf_ = libmp.from_pickable(val)

    def __cinit__(self):
        """
        Create a new mpf::

            sage: from mpmath import mpf
            sage: x = mpf()
        """
        MPF_init(&self.value)

    def __dealloc__(self):
        MPF_clear(&self.value)

    def __neg__(s):
        """
        Negate ``self``, rounded to the current working precision::

            sage: from mpmath import mpf
            sage: -mpf(2)
            mpf('-2.0')
        """
        cdef mpf r = mpf.__new__(mpf)
        MPF_neg(&r.value, &s.value)
        MPF_normalize(&r.value, global_opts)
        return r

    def __pos__(s):
        """
        Round the number to the current working precision::

            sage: from mpmath import mp, mpf
            sage: mp.prec = 200
            sage: x = mpf(1) / 3
            sage: x.man
            1071292029505993517027974728227441735014801995855195223534251
            sage: mp.prec = 53
            sage: (+x).man
            6004799503160661
            sage: print(+x)
            0.333333333333333
        """
        cdef mpf r = mpf.__new__(mpf)
        MPF_set(&r.value, &s.value)
        MPF_normalize(&r.value, global_opts)
        return r

    def __abs__(s):
        """
        Compute the absolute value, rounded to the current
        working precision::

            sage: from mpmath import mpf
            sage: abs(mpf(-2))
            mpf('2.0')
        """
        cdef mpf r = mpf.__new__(mpf)
        MPF_abs(&r.value, &s.value)
        MPF_normalize(&r.value, global_opts)
        return r

    def sqrt(s):
        """
        Compute the square root, rounded to the current
        working precision::

            sage: from mpmath import mpf
            sage: mpf(2).sqrt()
            mpf('1.4142135623730951')
        """
        cdef mpf r = mpf.__new__(mpf)
        MPF_sqrt(&r.value, &s.value, global_opts)
        return r

    def __round__(self, *args):
        return round(float(self), *args)

    def __richcmp__(self, other, int op):
        """
        Compare numbers::

            sage: from mpmath import mpf
            sage: mpf(3) > 2
            True
            sage: mpf(3) == 3
            True
            sage: mpf(3) == 4
            False
        """
        return binop(OP_RICHCMP+op, self, other, global_opts)


cdef class constant(mpf_base):
    """
    Represent a mathematical constant with dynamic precision.
    When printed or used in an arithmetic operation, a constant
    is converted to a regular mpf at the working precision. A
    regular mpf can also be obtained using the operation +x.
    """

    cdef public name, func, __doc__

    def __init__(self, func, name, docname=''):
        """
        Create a constant from a function computing an mpf
        tuple value::

            sage: from mpmath import mp, mpf
            sage: q = mp.constant(lambda prec, rnd: mpf(0.25)._mpf_, "quarter", "q")
            sage: q
            <quarter: 0.25~>
            sage: q + 1
            mpf('1.25')
        """
        self.name = name
        self.func = func
        self.__doc__ = getattr(function_docs, docname, '')

    def __call__(self, prec=None, dps=None, rounding=None):
        """
        Calling a constant is equivalent to rounding it. A
        custom precision and rounding direction can also be passed::

            sage: from mpmath import pi
            sage: print(pi(dps=5, rounding='d'))
            3.1415901184082
            sage: print(pi(dps=5, rounding='u'))
            3.14159393310547
        """
        prec2 = global_opts.prec
        rounding2 = rndmode_to_python(global_opts.rounding)
        if not prec: prec = prec2
        if not rounding: rounding = rounding2
        if dps: prec = dps_to_prec(dps)
        return global_context.make_mpf(self.func(prec, rounding))

    @property
    def _mpf_(self):
        """
        Return the tuple value of the constant as if rounded
        to an mpf at the present working precision::

            sage: from mpmath import pi
            sage: pi._mpf_
            (0, 884279719003555, -48, 50)
            sage: 884279719003555 / 2.0**48
            3.14159265358979
        """
        prec = global_opts.prec
        rounding = rndmode_to_python(global_opts.rounding)
        return self.func(prec, rounding)

    def __repr__(self):
        """
        Represent ``self`` as a string. With mp.pretty=False, the
        representation differs from that of an ordinary mpf::

            sage: from mpmath import mp
            sage: mp2 = mp.clone()
            sage: mp2.pretty = True
            sage: repr(mp2.pi)
            '3.14159265358979'
            sage: mp2.pretty = False
            sage: repr(mp2.pi)
            '<pi: 3.14159~>'
        """
        if global_context.pretty:
            return str(self)
        return "<%s: %s~>" % (self.name, global_context.nstr(self))

    def __bool__(self):
        """
        Return whether the constant is nonzero::

            sage: from mpmath import pi
            sage: bool(pi)
            True
        """
        return self._mpf_ != libmp.fzero

    def __neg__(self):
        """
        Negate the constant::

            sage: from mpmath import pi
            sage: -pi
            mpf('-3.1415926535897931')
        """
        return -mpf(self)

    def __pos__(self):
        """
        Instantiate the constant as an mpf::

            sage: from mpmath import pi
            sage: +pi
            mpf('3.1415926535897931')
        """
        return mpf(self)

    def __abs__(self):
        """
        Compute the absolute value of the constant::

            sage: from mpmath import pi
            sage: abs(pi)
            mpf('3.1415926535897931')
        """
        return abs(mpf(self))

    def sqrt(self):
        """
        Compute the square root of the constant::

            sage: from mpmath import pi
            sage: print(pi.sqrt())
            1.77245385090552
        """
        return mpf(self).sqrt()

    # XXX: optimize
    def to_fixed(self, prec):
        """
        Convert to a fixed-point integer::

            sage: from mpmath import pi
            sage: float(pi.to_fixed(10) / 2.0**10)
            3.140625
        """
        return libmp.to_fixed(self._mpf_, prec)

    def __getstate__(self):
        return libmp.to_pickable(self._mpf_)

    def __setstate__(self, val):
        self._mpf_ = libmp.from_pickable(val)

    # WHY is this method not inherited from the base class by Cython?
    def __hash__(self):
        """
        A constant hashes as if instantiated to a number::

            sage: from mpmath import pi
            sage: hash(pi) == hash(+pi)
            True
        """
        return libmp.mpf_hash(self._mpf_)

    def __richcmp__(self, other, int op):
        """
        A constant hashes as if instantiated to a number::

            sage: from mpmath import pi
            sage: pi == pi
            True
            sage: pi > 3.14
            True
            sage: pi < 3.14
            False
        """
        return binop(OP_RICHCMP+op, self, other, global_opts)


cdef class mpc(mpnumber):
    """
    An mpc represents a complex number using a pair of mpf:s (one
    for the real part and another for the imaginary part.) The mpc
    class behaves fairly similarly to Python's complex type.
    """

    cdef MPF re
    cdef MPF im

    def __init__(self, real=0, imag=0):
        """
        Create a new mpc::

            sage: from mpmath import mpc
            sage: mpc() == mpc(0,0) == mpc(1,0)-1 == 0
            True
        """
        cdef int typx, typy
        typx = MPF_set_any(&self.re, &self.im, real, global_opts, 1)
        if typx == 2:
            typy = 1
        else:
            typy = MPF_set_any(&self.im, &self.im, imag, global_opts, 1)
        if typx == 0 or typy != 1:
            raise TypeError

    def __cinit__(self):
        """
        Create a new mpc::

            sage: from mpmath import mpc
            sage: x = mpc()
        """
        MPF_init(&self.re)
        MPF_init(&self.im)

    def __dealloc__(self):
        MPF_clear(&self.re)
        MPF_clear(&self.im)

    def __reduce__(self):
        """
        Support pickling::

            sage: from mpmath import mpc
            sage: loads(dumps(mpc(1,3))) == mpc(1,3)
            True
        """
        return (mpc, (), self._mpc_)

    def __setstate__(self, val):
        """
        Support pickling::

            sage: from mpmath import mpc
            sage: loads(dumps(mpc(1,3))) == mpc(1,3)
            True
        """
        self._mpc_ = val[0], val[1]

    def __repr__(self):
        """
        TESTS::

            sage: from mpmath import mp
            sage: mp2 = mp.clone()
            sage: mp2.pretty = True
            sage: repr(mp2.mpc(2,3))
            '(2.0 + 3.0j)'
            sage: mp2.pretty = False
            sage: repr(mp2.mpc(2,3))
            "mpc(real='2.0', imag='3.0')"
        """
        if global_context.pretty:
            return str(self)
        re, im = self._mpc_
        n = repr_dps(global_opts.prec)
        return "mpc(real='%s', imag='%s')" % (to_str(re, n), to_str(im, n))

    def __str__(s):
        """
        TESTS::

            sage: from mpmath import mp
            sage: str(mp.mpc(2,3))
            '(2.0 + 3.0j)'
        """
        return "(%s)" % libmp.mpc_to_str(s._mpc_, global_context._str_digits)

    def __bool__(self):
        """
        TESTS::

            sage: from mpmath import mp
            sage: bool(mp.mpc(0,1))
            True
            sage: bool(mp.mpc(1,0))
            True
            sage: bool(mp.mpc(0,0))
            False
        """
        return self.re.special != S_ZERO or self.im.special != S_ZERO

    #def __complex__(self):
    #    a, b = self._mpc_
    #    return complex(libmp.to_float(a), libmp.to_float(b))

    def __complex__(self):
        """
        TESTS::

            sage: from mpmath import mp
            sage: complex(mp.mpc(1,2)) == complex(1,2)
            True
        """
        return complex(MPF_to_double(&self.re, False), MPF_to_double(&self.im, False))

    def _get_mpc(self):
        """
        Return tuple value of ``self``::

            sage: from mpmath import mp
            sage: mp.mpc(2,3)._mpc_
            ((0, 1, 1, 1), (0, 3, 0, 2))
        """
        return MPF_to_tuple(&self.re), MPF_to_tuple(&self.im)

    def _set_mpc(self, tuple v):
        """
        Set tuple value of ``self`` (warning: unsafe)::

            sage: from mpmath import mp
            sage: x = mp.mpc(2,3)
            sage: x._mpc_ = (x._mpc_[1], x._mpc_[0])
            sage: x
            mpc(real='3.0', imag='2.0')
        """
        MPF_set_tuple(&self.re, v[0])
        MPF_set_tuple(&self.im, v[1])

    _mpc_ = property(_get_mpc, _set_mpc, doc=_get_mpc.__doc__)

    @property
    def real(self):
        """
        Return the real part of ``self`` as an mpf::

            sage: from mpmath import mp
            sage: mp.mpc(1,2).real
            mpf('1.0')
        """
        cdef mpf r = mpf.__new__(mpf)
        MPF_set(&r.value, &self.re)
        return r

    @property
    def imag(self):
        """
        Return the imaginary part of ``self`` as an mpf::

            sage: from mpmath import mp
            sage: mp.mpc(1,2).imag
            mpf('2.0')
        """
        cdef mpf r = mpf.__new__(mpf)
        MPF_set(&r.value, &self.im)
        return r

    def __hash__(self):
        """
        Return the hash value of ``self``::

        EXAMPLES::

            sage: from mpmath import mp
            sage: hash(mp.mpc(2,3)) == hash(complex(2,3))
            True

        TESTS:

        Check that :issue:`31676` is fixed::

            sage: from mpmath import mpc
            sage: hash(mpc(1, -1)) == hash(mpc(-1, -1))  # should not return OverflowError: Python int too large to convert to C ssize_t
            False
        """
        return hash(libmp.mpc_hash(self._mpc_))

    def __neg__(s):
        """
        Negate the number::

            sage: from mpmath import mpc
            sage: -mpc(1,2)
            mpc(real='-1.0', imag='-2.0')
        """
        cdef mpc r = mpc.__new__(mpc)
        MPF_neg(&r.re, &s.re)
        MPF_neg(&r.im, &s.im)
        MPF_normalize(&r.re, global_opts)
        MPF_normalize(&r.im, global_opts)
        return r

    def conjugate(s):
        """
        Return the complex conjugate::

            sage: from mpmath import mpc
            sage: mpc(1,2).conjugate()
            mpc(real='1.0', imag='-2.0')
        """
        cdef mpc r = mpc.__new__(mpc)
        MPF_set(&r.re, &s.re)
        MPF_neg(&r.im, &s.im)
        MPF_normalize(&r.re, global_opts)
        MPF_normalize(&r.im, global_opts)
        return r

    def __pos__(s):
        """
        Round the number to the current working precision::

            sage: from mpmath import mp
            sage: mp.prec = 200
            sage: x = mp.mpc(1) / 3
            sage: x.real.man
            1071292029505993517027974728227441735014801995855195223534251
            sage: mp.prec = 53
            sage: +x
            mpc(real='0.33333333333333331', imag='0.0')
            sage: (+x).real.man
            6004799503160661
        """
        cdef mpc r = mpc.__new__(mpc)
        MPF_set(&r.re, &s.re)
        MPF_set(&r.im, &s.im)
        MPF_normalize(&r.re, global_opts)
        MPF_normalize(&r.im, global_opts)
        return r

    def __abs__(s):
        """
        Return the absolute value of ``self``::

            sage: from mpmath import mpc
            sage: abs(mpc(3,4))
            mpf('5.0')
        """
        cdef mpf r = mpf.__new__(mpf)
        MPF_hypot(&r.value, &s.re, &s.im, global_opts)
        return r

    def __richcmp__(self, other, int op):
        """
        Complex numbers can be compared for equality::

            sage: from mpmath import mpc
            sage: mpc(2,3) == complex(2,3)
            True
            sage: mpc(-2,3) == complex(2,3)
            False
            sage: mpc(-2,3) != complex(2,3)
            True
        """
        return binop(OP_RICHCMP+op, self, other, global_opts)


def hypsum_internal(int p, int q, param_types, str ztype, coeffs, z,
    long prec, long wp, long epsshift, dict magnitude_check, kwargs):
    """
    Internal summation routine for hypergeometric series (wraps
    extension function MPF_hypsum to handle mpf/mpc types).

    EXAMPLES::

        sage: from mpmath import mp  # indirect doctest
        sage: mp.dps = 15
        sage: print(mp.hyp1f1(1,2,3))
        6.36184564106256

    .. TODO::

        convert mpf/mpc parameters to fixed-point numbers here
        instead of converting to tuples within MPF_hypsum.
    """
    cdef mpf f
    cdef mpc c
    c = mpc.__new__(mpc)
    have_complex, magn = MPF_hypsum(&c.re, &c.im, p, q, param_types,
        ztype, coeffs, z, prec, wp, epsshift, magnitude_check, kwargs)
    if have_complex:
        v = c
    else:
        f = mpf.__new__(mpf)
        MPF_set(&f.value, &c.re)
        v = f
    return v, have_complex, magn
