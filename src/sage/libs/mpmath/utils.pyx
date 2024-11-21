"""
Utilities for Sage-mpmath interaction

Also patches some mpmath functions for speed
"""

from sage.ext.stdsage cimport PY_NEW

from sage.rings.integer cimport Integer
from sage.rings.real_mpfr cimport RealNumber
from sage.rings.complex_mpfr cimport ComplexNumber
from sage.structure.element cimport Element

from sage.libs.mpfr cimport *
from sage.libs.gmp.all cimport *

from sage.rings.real_mpfr cimport RealField

cpdef int bitcount(n) noexcept:
    """
    Bitcount of a Sage Integer or Python int/long.

    EXAMPLES::

        sage: from mpmath.libmp import bitcount
        sage: bitcount(0)
        0
        sage: bitcount(1)
        1
        sage: bitcount(100)
        7
        sage: bitcount(-100)
        7
        sage: bitcount(2r)
        2
        sage: bitcount(2L)
        2
    """
    cdef Integer m
    if isinstance(n, Integer):
        m = <Integer>n
    else:
        m = Integer(n)
    if mpz_sgn(m.value) == 0:
        return 0
    return mpz_sizeinbase(m.value, 2)

cpdef isqrt(n):
    """
    Square root (rounded to floor) of a Sage Integer or Python int/long.
    The result is a Sage Integer.

    EXAMPLES::

        sage: from mpmath.libmp import isqrt
        sage: isqrt(0)
        0
        sage: isqrt(100)
        10
        sage: isqrt(10)
        3
        sage: isqrt(10r)
        3
        sage: isqrt(10L)
        3
    """
    cdef Integer m, y
    if isinstance(n, Integer):
        m = <Integer>n
    else:
        m = Integer(n)
    if mpz_sgn(m.value) < 0:
        raise ValueError("square root of negative integer not defined.")
    y = PY_NEW(Integer)
    mpz_sqrt(y.value, m.value)
    return y

cpdef from_man_exp(man, exp, long prec=0, str rnd='d'):
    """
    Create normalized mpf value tuple from mantissa and exponent.

    With prec > 0, rounds the result in the desired direction
    if necessary.

    EXAMPLES::

        sage: from mpmath.libmp import from_man_exp
        sage: from_man_exp(-6, -1)
        (1, 3, 0, 2)
        sage: from_man_exp(-6, -1, 1, 'd')
        (1, 1, 1, 1)
        sage: from_man_exp(-6, -1, 1, 'u')
        (1, 1, 2, 1)
    """
    cdef Integer res
    cdef long bc
    res = Integer(man)
    bc = mpz_sizeinbase(res.value, 2)
    if not prec:
        prec = bc
    if mpz_sgn(res.value) < 0:
        mpz_neg(res.value, res.value)
        return normalize(1, res, exp, bc, prec, rnd)
    else:
        return normalize(0, res, exp, bc, prec, rnd)

cpdef normalize(long sign, Integer man, exp, long bc, long prec, str rnd):
    """
    Create normalized mpf value tuple from full list of components.

    EXAMPLES::

        sage: from mpmath.libmp import normalize
        sage: normalize(0, 4, 5, 3, 53, 'n')
        (0, 1, 7, 1)
    """
    cdef long shift
    cdef Integer res
    cdef unsigned long trail
    if mpz_sgn(man.value) == 0:
        from mpmath.libmp import fzero
        return fzero
    if bc <= prec and mpz_odd_p(man.value):
        return (sign, man, exp, bc)
    shift = bc - prec
    res = PY_NEW(Integer)
    if shift > 0:
        if rnd == 'n':
            if mpz_tstbit(man.value, shift-1) and (mpz_tstbit(man.value, shift)
                or (mpz_scan1(man.value, 0) < (shift-1))):
                mpz_cdiv_q_2exp(res.value, man.value, shift)
            else:
                mpz_fdiv_q_2exp(res.value, man.value, shift)
        elif rnd == 'd':
            mpz_fdiv_q_2exp(res.value, man.value, shift)
        elif rnd == 'f':
            if sign: mpz_cdiv_q_2exp(res.value, man.value, shift)
            else:    mpz_fdiv_q_2exp(res.value, man.value, shift)
        elif rnd == 'c':
            if sign: mpz_fdiv_q_2exp(res.value, man.value, shift)
            else:    mpz_cdiv_q_2exp(res.value, man.value, shift)
        elif rnd == 'u':
            mpz_cdiv_q_2exp(res.value, man.value, shift)
        exp += shift
    else:
        mpz_set(res.value, man.value)
    # Strip trailing bits
    trail = mpz_scan1(res.value, 0)
    if 0 < trail < bc:
        mpz_tdiv_q_2exp(res.value, res.value, trail)
        exp += trail
    bc = mpz_sizeinbase(res.value, 2)
    return (sign, res, int(exp), bc)

cdef mpfr_from_mpfval(mpfr_t res, tuple x):
    """
    Set value of an MPFR number (in place) to that of a given mpmath mpf
    data tuple.
    """
    cdef int sign
    cdef Integer man
    cdef long exp
    cdef long bc
    sign, man, exp, bc = x
    if man:
        mpfr_set_z(res, man.value, MPFR_RNDZ)
        if sign:
            mpfr_neg(res, res, MPFR_RNDZ)
        mpfr_mul_2si(res, res, exp, MPFR_RNDZ)
        return
    from mpmath.libmp import finf, fninf
    if exp == 0:
        mpfr_set_ui(res, 0, MPFR_RNDZ)
    elif x == finf:
        mpfr_set_inf(res, 1)
    elif x == fninf:
        mpfr_set_inf(res, -1)
    else:
        mpfr_set_nan(res)

cdef mpfr_to_mpfval(mpfr_t value):
    """
    Given an MPFR value, return an mpmath mpf data tuple representing
    the same number.
    """
    if mpfr_nan_p(value):
        from mpmath.libmp import fnan
        return fnan
    if mpfr_inf_p(value):
        from mpmath.libmp import finf, fninf
        if mpfr_sgn(value) > 0:
            return finf
        else:
            return fninf
    if mpfr_sgn(value) == 0:
        from mpmath.libmp import fzero
        return fzero
    sign = 0
    cdef Integer man = PY_NEW(Integer)
    exp = mpfr_get_z_exp(man.value, value)
    if mpz_sgn(man.value) < 0:
        mpz_neg(man.value, man.value)
        sign = 1
    cdef unsigned long trailing
    trailing = mpz_scan1(man.value, 0)
    if trailing:
        mpz_tdiv_q_2exp(man.value, man.value, trailing)
        exp += trailing
    bc = mpz_sizeinbase(man.value, 2)
    return (sign, man, int(exp), bc)


def mpmath_to_sage(x, prec):
    """
    Convert any mpmath number (mpf or mpc) to a Sage RealNumber or
    ComplexNumber of the given precision.

    EXAMPLES::

        sage: import sage.libs.mpmath.all as a
        sage: a.mpmath_to_sage(a.mpf('2.5'), 53)
        2.50000000000000
        sage: a.mpmath_to_sage(a.mpc('2.5','-3.5'), 53)
        2.50000000000000 - 3.50000000000000*I
        sage: a.mpmath_to_sage(a.mpf('inf'), 53)
        +infinity
        sage: a.mpmath_to_sage(a.mpf('-inf'), 53)
        -infinity
        sage: a.mpmath_to_sage(a.mpf('nan'), 53)
        NaN
        sage: a.mpmath_to_sage(a.mpf('0'), 53)
        0.000000000000000

    A real example::

        sage: RealField(100)(pi)
        3.1415926535897932384626433833
        sage: t = RealField(100)(pi)._mpmath_(); t
        mpf('3.1415926535897932')
        sage: a.mpmath_to_sage(t, 100)
        3.1415926535897932384626433833

    We can ask for more precision, but the result is undefined::

        sage: a.mpmath_to_sage(t, 140) # random
        3.1415926535897932384626433832793333156440
        sage: ComplexField(140)(pi)
        3.1415926535897932384626433832795028841972

    A complex example::

        sage: ComplexField(100)([0, pi])
        3.1415926535897932384626433833*I
        sage: t = ComplexField(100)([0, pi])._mpmath_(); t
        mpc(real='0.0', imag='3.1415926535897932')
        sage: sage.libs.mpmath.all.mpmath_to_sage(t, 100)
        3.1415926535897932384626433833*I

    Again, we can ask for more precision, but the result is undefined::

        sage: sage.libs.mpmath.all.mpmath_to_sage(t, 140) # random
        3.1415926535897932384626433832793333156440*I
        sage: ComplexField(140)([0, pi])
        3.1415926535897932384626433832795028841972*I
    """
    cdef RealNumber y
    cdef ComplexNumber z
    if hasattr(x, "_mpf_"):
        y = RealField(prec)()
        mpfr_from_mpfval(y.value, x._mpf_)
        return y
    elif hasattr(x, "_mpc_"):
        from sage.rings.complex_mpfr import ComplexField
        z = ComplexField(prec)(0)
        re, im = x._mpc_
        mpfr_from_mpfval(z.__re, re)
        mpfr_from_mpfval(z.__im, im)
        return z
    else:
        raise TypeError("cannot convert %r to Sage", x)


def sage_to_mpmath(x, prec):
    """
    Convert any Sage number that can be coerced into a RealNumber
    or ComplexNumber of the given precision into an mpmath mpf or mpc.
    Integers are currently converted to int.

    Lists, tuples and dicts passed as input are converted
    recursively.

    EXAMPLES::

        sage: import sage.libs.mpmath.all as a
        sage: a.mp.dps = 15
        sage: print(a.sage_to_mpmath(2/3, 53))
        0.666666666666667
        sage: print(a.sage_to_mpmath(2./3, 53))
        0.666666666666667
        sage: print(a.sage_to_mpmath(3+4*I, 53))
        (3.0 + 4.0j)
        sage: print(a.sage_to_mpmath(1+pi, 53))
        4.14159265358979
        sage: a.sage_to_mpmath(infinity, 53)
        mpf('+inf')
        sage: a.sage_to_mpmath(-infinity, 53)
        mpf('-inf')
        sage: a.sage_to_mpmath(NaN, 53)
        mpf('nan')
        sage: a.sage_to_mpmath(0, 53)
        0
        sage: a.sage_to_mpmath([0.5, 1.5], 53)
        [mpf('0.5'), mpf('1.5')]
        sage: a.sage_to_mpmath((0.5, 1.5), 53)
        (mpf('0.5'), mpf('1.5'))
        sage: a.sage_to_mpmath({'n':0.5}, 53)
        {'n': mpf('0.5')}
    """
    cdef RealNumber y
    if isinstance(x, Element):
        if isinstance(x, Integer):
            return int(<Integer>x)
        try:
            if isinstance(x, RealNumber):
                return x._mpmath_()
            else:
                x = RealField(prec)(x)
                return x._mpmath_()
        except TypeError:
            if isinstance(x, ComplexNumber):
                return x._mpmath_()
            else:
                from sage.rings.complex_mpfr import ComplexField
                x = ComplexField(prec)(x)
                return x._mpmath_()
    if isinstance(x, (tuple, list)):
        return type(x)([sage_to_mpmath(v, prec) for v in x])
    if isinstance(x, dict):
        return dict([(k, sage_to_mpmath(v, prec)) for (k, v) in x.items()])
    return x


def call(func, *args, **kwargs):
    """
    Call an mpmath function with Sage objects as inputs and
    convert the result back to a Sage real or complex number.

    By default, a RealNumber or ComplexNumber with the current
    working precision of mpmath (mpmath.mp.prec) will be returned.

    If prec=n is passed among the keyword arguments, the temporary
    working precision will be set to n and the result will also
    have this precision.

    If parent=P is passed, P.prec() will be used as working
    precision and the result will be coerced to P (or the
    corresponding complex field if necessary).

    Arguments should be Sage objects that can be coerced into RealField
    or ComplexField elements. Arguments may also be tuples, lists or
    dicts (which are converted recursively), or any type that mpmath
    understands natively (e.g. Python floats, strings for options).

    EXAMPLES::

        sage: import sage.libs.mpmath.all as a
        sage: a.mp.prec = 53
        sage: a.call(a.erf, 3+4*I)
        -120.186991395079 - 27.7503372936239*I
        sage: a.call(a.polylog, 2, 1/3+4/5*I)
        0.153548951541433 + 0.875114412499637*I
        sage: a.call(a.barnesg, 3+4*I)
        -0.000676375932234244 - 0.0000442236140124728*I
        sage: a.call(a.barnesg, -4)
        0.000000000000000
        sage: a.call(a.hyper, [2,3], [4,5], 1/3)
        1.10703578162508
        sage: a.call(a.hyper, [2,3], [4,(2,3)], 1/3)
        1.95762943509305
        sage: a.call(a.quad, a.erf, [0,1])
        0.486064958112256
        sage: a.call(a.gammainc, 3+4*I, 2/3, 1-pi*I, prec=100)
        -274.18871130777160922270612331 + 101.59521032382593402947725236*I
        sage: x = (3+4*I).n(100)
        sage: y = (2/3).n(100)
        sage: z = (1-pi*I).n(100)
        sage: a.call(a.gammainc, x, y, z, prec=100)
        -274.18871130777160922270612331 + 101.59521032382593402947725236*I
        sage: a.call(a.erf, infinity)
        1.00000000000000
        sage: a.call(a.erf, -infinity)
        -1.00000000000000
        sage: a.call(a.gamma, infinity)
        +infinity
        sage: a.call(a.polylog, 2, 1/2, parent=RR)
        0.582240526465012
        sage: a.call(a.polylog, 2, 2, parent=RR)
        2.46740110027234 - 2.17758609030360*I
        sage: a.call(a.polylog, 2, 1/2, parent=RealField(100))
        0.58224052646501250590265632016
        sage: a.call(a.polylog, 2, 2, parent=RealField(100))
        2.4674011002723396547086227500 - 2.1775860903036021305006888982*I
        sage: a.call(a.polylog, 2, 1/2, parent=CC)
        0.582240526465012
        sage: type(_)
        <class 'sage.rings.complex_mpfr.ComplexNumber'>
        sage: a.call(a.polylog, 2, 1/2, parent=RDF)
        0.5822405264650125
        sage: type(_)
        <class 'sage.rings.real_double...RealDoubleElement...'>

    Check that :issue:`11885` is fixed::

        sage: a.call(a.ei, 1.0r, parent=float)
        1.8951178163559366

    Check that :issue:`14984` is fixed::

        sage: a.call(a.log, -1.0r, parent=float)
        3.141592653589793j
    """
    from mpmath import mp
    orig = mp.prec
    prec = kwargs.pop('prec', orig)
    parent = kwargs.pop('parent', None)
    if parent is not None:
        try:
            prec = parent.prec()
        except AttributeError:
            pass
    prec2 = prec + 20
    args = sage_to_mpmath(args, prec2)
    kwargs = sage_to_mpmath(kwargs, prec2)
    try:
        mp.prec = prec
        y = func(*args, **kwargs)
    finally:
        mp.prec = orig
    y = mpmath_to_sage(y, prec)
    if parent is None:
        return y
    try:
        return parent(y)
    except TypeError as error:
        try:
            return parent.complex_field()(y)
        except AttributeError:
            if parent is float:
                return complex(y)
            else:
                raise TypeError(error)
