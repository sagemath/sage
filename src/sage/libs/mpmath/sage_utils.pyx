# sage_setup: distribution = sagemath-mpmath
r"""
Utility functions called by the Sage library
"""

from sage.ext.stdsage cimport PY_NEW
from sage.libs.mpfr cimport *
from sage.libs.gmp.all cimport *
from sage.libs.mpmath._vendor.mpmath import mp
from sage.libs.mpmath._vendor.mpmath.libmp import finf, fnan, fninf, fzero
from sage.rings.complex_mpfr cimport ComplexNumber
from sage.rings.complex_mpfr import ComplexField
from sage.rings.integer cimport Integer
from sage.rings.real_mpfr cimport RealField, RealNumber
from sage.structure.element cimport Element


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
        return fnan
    if mpfr_inf_p(value):
        if mpfr_sgn(value) > 0:
            return finf
        else:
            return fninf
    if mpfr_sgn(value) == 0:
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

        sage: RealField(100)(pi)                                                        # needs sage.symbolic
        3.1415926535897932384626433833
        sage: t = RealField(100)(pi)._mpmath_(); t                                      # needs sage.symbolic
        mpf('3.1415926535897932')
        sage: a.mpmath_to_sage(t, 100)                                                  # needs sage.symbolic
        3.1415926535897932384626433833

    We can ask for more precision, but the result is undefined::

        sage: a.mpmath_to_sage(t, 140)  # random                                        # needs sage.symbolic
        3.1415926535897932384626433832793333156440
        sage: ComplexField(140)(pi)                                                     # needs sage.symbolic
        3.1415926535897932384626433832795028841972

    A complex example::

        sage: ComplexField(100)([0, pi])                                                # needs sage.symbolic
        3.1415926535897932384626433833*I
        sage: t = ComplexField(100)([0, pi])._mpmath_(); t                              # needs sage.symbolic
        mpc(real='0.0', imag='3.1415926535897932')
        sage: sage.libs.mpmath.all.mpmath_to_sage(t, 100)                               # needs sage.symbolic
        3.1415926535897932384626433833*I

    Again, we can ask for more precision, but the result is undefined::

        sage: sage.libs.mpmath.all.mpmath_to_sage(t, 140)  # random                     # needs sage.symbolic
        3.1415926535897932384626433832793333156440*I
        sage: ComplexField(140)([0, pi])                                                # needs sage.symbolic
        3.1415926535897932384626433832795028841972*I
    """
    cdef RealNumber y
    cdef ComplexNumber z
    if hasattr(x, "_mpf_"):
        y = RealField(prec)()
        mpfr_from_mpfval(y.value, x._mpf_)
        return y
    elif hasattr(x, "_mpc_"):
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
        sage: print(a.sage_to_mpmath(3+4*I, 53))                                        # needs sage.symbolic
        (3.0 + 4.0j)
        sage: print(a.sage_to_mpmath(1+pi, 53))                                         # needs sage.symbolic
        4.14159265358979
        sage: a.sage_to_mpmath(infinity, 53)
        mpf('+inf')
        sage: a.sage_to_mpmath(-infinity, 53)
        mpf('-inf')
        sage: a.sage_to_mpmath(NaN, 53)                                                 # needs sage.symbolic
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
        sage: a.call(a.erf, 3+4*I)                                                      # needs sage.symbolic
        -120.186991395079 - 27.7503372936239*I
        sage: a.call(a.polylog, 2, 1/3+4/5*I)                                           # needs sage.symbolic
        0.153548951541433 + 0.875114412499637*I
        sage: a.call(a.barnesg, 3+4*I)                                                  # needs sage.symbolic
        -0.000676375932234244 - 0.0000442236140124728*I
        sage: a.call(a.barnesg, -4)
        0.000000000000000
        sage: a.call(a.hyper, [2,3], [4,5], 1/3)
        1.10703578162508
        sage: a.call(a.hyper, [2,3], [4,(2,3)], 1/3)
        1.95762943509305
        sage: a.call(a.quad, a.erf, [0,1])
        0.486064958112256
        sage: a.call(a.gammainc, 3+4*I, 2/3, 1-pi*I, prec=100)                          # needs sage.symbolic
        -274.18871130777160922270612331 + 101.59521032382593402947725236*I
        sage: x = (3+4*I).n(100)                                                        # needs sage.symbolic
        sage: y = (2/3).n(100)                                                          # needs sage.symbolic
        sage: z = (1-pi*I).n(100)                                                       # needs sage.symbolic
        sage: a.call(a.gammainc, x, y, z, prec=100)                                     # needs sage.symbolic
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
