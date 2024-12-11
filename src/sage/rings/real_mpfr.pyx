r"""
Arbitrary precision floating point real numbers using GNU MPFR

AUTHORS:

- Kyle Schalm (2005-09)

- William Stein: bug fixes, examples, maintenance

- Didier Deshommes (2006-03-19): examples

- David Harvey (2006-09-20): compatibility with Element._parent

- William Stein (2006-10): default printing truncates to avoid base-2
  rounding confusing (fix suggested by Bill Hart)

- Didier Deshommes: special constructor for QD numbers

- Paul Zimmermann (2008-01): added new functions from mpfr-2.3.0,
  replaced some, e.g., sech = 1/cosh, by their original mpfr version.

- Carl Witty (2008-02): define floating-point rank and associated
  functions; add some documentation

- Robert Bradshaw (2009-09): decimal literals, optimizations

- Jeroen Demeyer (2012-05-27): set the MPFR exponent range to the
  maximal possible value (:issue:`13033`)

- Travis Scrimshaw (2012-11-02): Added doctests for full coverage

- Eviatar Bach (2013-06): Fixing numerical evaluation of log_gamma

- Vincent Klein (2017-06): RealNumber constructor support gmpy2.mpfr
  , gmpy2.mpq or gmpy2.mpz parameter.
  Add __mpfr__ to class RealNumber.

This is a binding for the MPFR arbitrary-precision floating point
library.

We define a class :class:`RealField`, where each instance of
:class:`RealField` specifies a field of floating-point
numbers with a specified precision and rounding mode. Individual
floating-point numbers are of :class:`RealNumber`.

In Sage (as in MPFR), floating-point numbers of precision
`p` are of the form `s m 2^{e-p}`, where
`s \in \{-1, 1\}`, `2^{p-1} \leq m < 2^p`, and
`-2^B + 1 \leq e \leq 2^B - 1` where `B = 30` on 32-bit systems
and `B = 62` on 64-bit systems;
additionally, there are the special values ``+0``, ``-0``,
``+infinity``, ``-infinity`` and ``NaN`` (which stands for Not-a-Number).

Operations in this module which are direct wrappers of MPFR
functions are "correctly rounded"; we briefly describe what this
means. Assume that you could perform the operation exactly, on real
numbers, to get a result `r`. If this result can be
represented as a floating-point number, then we return that
number.

Otherwise, the result `r` is between two floating-point
numbers. For the directed rounding modes (round to plus infinity,
round to minus infinity, round to zero), we return the
floating-point number in the indicated direction from `r`.
For round to nearest, we return the floating-point number which is
nearest to `r`.

This leaves one case unspecified: in round to nearest mode, what
happens if `r` is exactly halfway between the two nearest
floating-point numbers? In that case, we round to the number with
an even mantissa (the mantissa is the number `m` in the
representation above).

Consider the ordered set of floating-point numbers of precision
`p`. (Here we identify ``+0`` and
``-0``, and ignore ``NaN``.) We can give a
bijection between these floating-point numbers and a segment of the
integers, where 0 maps to 0 and adjacent floating-point numbers map
to adjacent integers. We call the integer corresponding to a given
floating-point number the "floating-point rank" of the number.
(This is not standard terminology; I just made it up.)

EXAMPLES:

A difficult conversion::

    sage: RR(sys.maxsize)
    9.22337203685478e18      # 64-bit
    2.14748364700000e9       # 32-bit

TESTS::

    sage: -1e30
    -1.00000000000000e30
    sage: (-1. + 2^-52).hex()
    '-0xf.ffffffffffffp-4'

Make sure we don't have a new field for every new literal::

    sage: parent(2.0) is parent(2.0)
    True
    sage: RealField(100, rnd='RNDZ') is RealField(100, rnd='RNDD')
    False
    sage: RealField(100, rnd='RNDZ') is RealField(100, rnd='RNDZ')
    True
    sage: RealField(100, rnd='RNDZ') is RealField(100, rnd=1)
    True
"""

# *****************************************************************************
#       Copyright (C) 2005-2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

import math  # for log
import operator
import re
import sys

from cpython.object cimport Py_NE, Py_EQ
from cysignals.signals cimport sig_on, sig_off

import sage.arith.misc
import sage.misc.weak_dict
import sage.rings.abc
import sage.rings.infinity
import sage.rings.rational_field

from sage.arith.constants cimport M_LN2_LN10
from sage.arith.long cimport is_small_python_int
from sage.arith.numerical_approx cimport digits_to_bits
from sage.categories.map cimport Map
from sage.cpython.string cimport char_to_str, str_to_bytes
from sage.ext.stdsage cimport PY_NEW
from sage.libs.gmp.mpz cimport *
from sage.libs.gmp.pylong cimport mpz_set_pylong
from sage.libs.mpfr cimport *
from sage.libs.mpmath.utils cimport mpfr_to_mpfval
from sage.misc.randstate cimport randstate, current_randstate

from sage.structure.element cimport Element
from sage.structure.parent cimport Parent
from sage.structure.element cimport have_same_parent
from sage.structure.richcmp cimport rich_to_bool_sgn
cdef bin_op
from sage.structure.element import bin_op

from sage.libs.mpmath.utils cimport mpfr_to_mpfval

from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.real_double cimport RealDoubleElement

try:
    from cypari2 import Gen
    from sage.libs.pari.convert_sage_real_mpfr import new_gen_from_real_mpfr_element, set_real_mpfr_element_from_gen
except ImportError:
    Gen = ()

cimport gmpy2
gmpy2.import_gmpy2()


# ****************************************************************************
#
#       Implementation
#
# ****************************************************************************

_re_skip_zeroes = re.compile(r'^(.+?)0*$')

cdef object numpy_double_interface = {'typestr': '=f8'}
cdef object numpy_object_interface = {'typestr': '|O'}

# Avoid signal handling for cheap operations when the
# precision is below this threshold.
cdef enum:
    SIG_PREC_THRESHOLD = 1000


# ****************************************************************************
#
#       External Python access to constants
#
# ****************************************************************************

def mpfr_prec_min():
    """
    Return the mpfr variable ``MPFR_PREC_MIN``.

    EXAMPLES::

        sage: from sage.rings.real_mpfr import mpfr_prec_min
        sage: mpfr_prec_min()
        1
        sage: R = RealField(2)
        sage: R(2) + R(1)
        3.0
        sage: R(4) + R(1)
        4.0

        sage: R = RealField(0)
        Traceback (most recent call last):
        ...
        ValueError: prec (=0) must be >= 1 and <= ...
    """
    return MPFR_PREC_MIN


# see Issue #11666 for the origin of this magical constant
def mpfr_prec_max():
    """
    Return the mpfr variable ``MPFR_PREC_MAX``.

    EXAMPLES::

        sage: from sage.rings.real_mpfr import mpfr_prec_max
        sage: mpfr_prec_max()
        2147483391              # 32-bit
        9223372036854775551     # 64-bit

        sage: R = RealField(2^31-257); R
        Real Field with 2147483391 bits of precision

        sage: R = RealField(2^31-256)
        Traceback (most recent call last):                     # 32-bit
        ...                                                    # 32-bit
        ValueError: prec (=...) must be >= 1 and <= ...        # 32-bit
    """
    return MPFR_PREC_MAX


def mpfr_get_exp_min():
    """
    Return the current minimal exponent for MPFR numbers.

    EXAMPLES::

        sage: from sage.rings.real_mpfr import mpfr_get_exp_min
        sage: mpfr_get_exp_min()
        -1073741823            # 32-bit
        -4611686018427387903   # 64-bit
        sage: 0.5 >> (-mpfr_get_exp_min())
        2.38256490488795e-323228497            # 32-bit
        8.50969131174084e-1388255822130839284  # 64-bit
        sage: 0.5 >> (-mpfr_get_exp_min()+1)
        0.000000000000000
    """
    return mpfr_get_emin()


def mpfr_get_exp_max():
    """
    Return the current maximal exponent for MPFR numbers.

    EXAMPLES::

        sage: from sage.rings.real_mpfr import mpfr_get_exp_max
        sage: mpfr_get_exp_max()
        1073741823            # 32-bit
        4611686018427387903   # 64-bit
        sage: 0.5 << mpfr_get_exp_max()
        1.04928935823369e323228496            # 32-bit
        2.93782689455579e1388255822130839282  # 64-bit
        sage: 0.5 << (mpfr_get_exp_max()+1)
        +infinity
    """
    return mpfr_get_emax()


def mpfr_set_exp_min(mp_exp_t e):
    """
    Set the minimal exponent for MPFR numbers.

    EXAMPLES::

        sage: from sage.rings.real_mpfr import mpfr_get_exp_min, mpfr_set_exp_min
        sage: old = mpfr_get_exp_min()
        sage: mpfr_set_exp_min(-1000)
        sage: 0.5 >> 1000
        4.66631809251609e-302
        sage: 0.5 >> 1001
        0.000000000000000
        sage: mpfr_set_exp_min(old)
        sage: 0.5 >> 1001
        2.33315904625805e-302
    """
    if mpfr_set_emin(e) != 0:
        raise OverflowError("bad value for mpfr_set_exp_min()")


def mpfr_set_exp_max(mp_exp_t e):
    """
    Set the maximal exponent for MPFR numbers.

    EXAMPLES::

        sage: from sage.rings.real_mpfr import mpfr_get_exp_max, mpfr_set_exp_max
        sage: old = mpfr_get_exp_max()
        sage: mpfr_set_exp_max(1000)
        sage: 0.5 << 1000
        5.35754303593134e300
        sage: 0.5 << 1001
        +infinity
        sage: mpfr_set_exp_max(old)
        sage: 0.5 << 1001
        1.07150860718627e301
    """
    if mpfr_set_emax(e) != 0:
        raise OverflowError("bad value for mpfr_set_exp_max()")


def mpfr_get_exp_min_min():
    """
    Get the minimal value allowed for :func:`mpfr_set_exp_min`.

    EXAMPLES::

        sage: from sage.rings.real_mpfr import mpfr_get_exp_min_min, mpfr_set_exp_min
        sage: mpfr_get_exp_min_min()
        -1073741823            # 32-bit
        -4611686018427387903   # 64-bit

    This is really the minimal value allowed::

        sage: mpfr_set_exp_min(mpfr_get_exp_min_min() - 1)
        Traceback (most recent call last):
        ...
        OverflowError: bad value for mpfr_set_exp_min()
    """
    return mpfr_get_emin_min()


def mpfr_get_exp_max_max():
    """
    Get the maximal value allowed for :func:`mpfr_set_exp_max`.

    EXAMPLES::

        sage: from sage.rings.real_mpfr import mpfr_get_exp_max_max, mpfr_set_exp_max
        sage: mpfr_get_exp_max_max()
        1073741823            # 32-bit
        4611686018427387903   # 64-bit

    This is really the maximal value allowed::

        sage: mpfr_set_exp_max(mpfr_get_exp_max_max() + 1)
        Traceback (most recent call last):
        ...
        OverflowError: bad value for mpfr_set_exp_max()
    """
    return mpfr_get_emax_max()


# On Sage startup, set the exponent range to the maximum allowed
mpfr_set_exp_min(mpfr_get_emin_min())
mpfr_set_exp_max(mpfr_get_emax_max())

# *****************************************************************************
#
#       Real Field
#
# *****************************************************************************
# The real field is in Cython, so mpfr elements will have access to
# their parent via direct C calls, which will be faster.

from sage.arith.long cimport (pyobject_to_long, integer_check_long_py,
                              ERR_OVERFLOW)
cdef dict rounding_modes = dict(RNDN=MPFR_RNDN, RNDZ=MPFR_RNDZ,
                                RNDD=MPFR_RNDD, RNDU=MPFR_RNDU,
                                RNDA=MPFR_RNDA, RNDF=MPFR_RNDF)

cdef double LOG_TEN_TWO_PLUS_EPSILON = 3.321928094887363
# a small overestimate of log(10,2)

cdef object RealField_cache = sage.misc.weak_dict.WeakValueDictionary()

cpdef RealField(mpfr_prec_t prec=53, int sci_not=0, rnd=MPFR_RNDN):
    """
    RealField(prec, sci_not, rnd):

    INPUT:

    - ``prec`` -- integer (default: 53); precision ``prec`` is
      the number of bits used to represent the mantissa of a
      floating-point number. The precision can be any integer between
      :func:`mpfr_prec_min()` and :func:`mpfr_prec_max()`. In the current
      implementation, :func:`mpfr_prec_min()` is equal to 2.

    - ``sci_not`` -- boolean (default: ``False``); if ``True``, always display
      using scientific notation. If ``False``, display using scientific
      notation only for very large or very small numbers.

    - ``rnd`` -- string; the rounding mode:

      - ``'RNDN'`` -- (default) round to nearest (ties go to the even
        number): Knuth says this is the best choice to prevent "floating
        point drift"
      - ``'RNDD'`` -- round towards minus infinity
      - ``'RNDZ'`` -- round towards zero
      - ``'RNDU'`` -- round towards plus infinity
      - ``'RNDA'`` -- round away from zero
      - ``'RNDF'`` -- faithful rounding (currently experimental; not
        guaranteed correct for every operation)
      - for specialized applications, the rounding mode can also be
        given as an integer value of type ``mpfr_rnd_t``. However, the
        exact values are unspecified.

    EXAMPLES::

        sage: RealField(10)
        Real Field with 10 bits of precision
        sage: RealField()
        Real Field with 53 bits of precision
        sage: RealField(100000)
        Real Field with 100000 bits of precision

    Here we show the effect of rounding::

        sage: R17d = RealField(17,rnd='RNDD')
        sage: a = R17d(1)/R17d(3); a.exact_rational()
        87381/262144
        sage: R17u = RealField(17,rnd='RNDU')
        sage: a = R17u(1)/R17u(3); a.exact_rational()
        43691/131072

    .. NOTE::

       The default precision is 53, since according to the MPFR
       manual: 'mpfr should be able to exactly reproduce all
       computations with double-precision machine floating-point
       numbers (double type in C), except the default exponent range
       is much wider and subnormal numbers are not implemented.'

    .. SEEALSO::

        - :mod:`sage.rings.real_mpfr`
        - :class:`sage.rings.real_arb.RealBallField` (real numbers with rigorous
          error bounds)
    """
    # We allow specifying the rounding mode as string or integer.
    # But we pass an integer to __init__
    cdef long r
    try:
        r = pyobject_to_long(rnd)
    except TypeError:
        try:
            r = rounding_modes[rnd]
        except KeyError:
            msg = "rounding mode (={!r}) must be one of {}"
            raise ValueError(msg.format(rnd, sorted(rounding_modes)))

    try:
        return RealField_cache[prec, sci_not, r]
    except KeyError:
        R = RealField_class(prec=prec, sci_not=sci_not, rnd=r)
        RealField_cache[prec, sci_not, r] = R
        return R


cdef class RealField_class(sage.rings.abc.RealField):
    """
    An approximation to the field of real numbers using floating point
    numbers with any specified precision. Answers derived from
    calculations in this approximation may differ from what they would
    be if those calculations were performed in the true field of real
    numbers. This is due to the rounding errors inherent to finite
    precision calculations.

    .. SEEALSO::

        - :mod:`sage.rings.real_mpfr`
        - :class:`sage.rings.real_arb.RealBallField` (real numbers with rigorous
          error bounds)
        - :mod:`sage.rings.complex_mpfr`
    """
    def __init__(self, mpfr_prec_t prec=53, int sci_not=0, long rnd=MPFR_RNDN):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: RealField()
            Real Field with 53 bits of precision
            sage: RealField(100000)
            Real Field with 100000 bits of precision
            sage: RealField(17,rnd='RNDD')
            Real Field with 17 bits of precision and rounding RNDD

        TESTS:

        Test the various rounding modes::

            sage: RealField(100, rnd="RNDN")
            Real Field with 100 bits of precision
            sage: RealField(100, rnd="RNDZ")
            Real Field with 100 bits of precision and rounding RNDZ
            sage: RealField(100, rnd="RNDU")
            Real Field with 100 bits of precision and rounding RNDU
            sage: RealField(100, rnd="RNDD")
            Real Field with 100 bits of precision and rounding RNDD
            sage: RealField(100, rnd="RNDA")
            Real Field with 100 bits of precision and rounding RNDA
            sage: RealField(100, rnd="RNDF")
            Real Field with 100 bits of precision and rounding RNDF
            sage: RealField(100, rnd=0)
            Real Field with 100 bits of precision
            sage: RealField(100, rnd=1)
            Real Field with 100 bits of precision and rounding RNDZ
            sage: RealField(100, rnd=2)
            Real Field with 100 bits of precision and rounding RNDU
            sage: RealField(100, rnd=3)
            Real Field with 100 bits of precision and rounding RNDD
            sage: RealField(100, rnd=4)
            Real Field with 100 bits of precision and rounding RNDA
            sage: RealField(100, rnd=5)
            Real Field with 100 bits of precision and rounding RNDF
            sage: RealField(100, rnd=3.14)
            Traceback (most recent call last):
            ...
            ValueError: rounding mode (=3.14000000000000) must be one of ['RNDA', 'RNDD', 'RNDF', 'RNDN', 'RNDU', 'RNDZ']
            sage: RealField(100, rnd=6)
            Traceback (most recent call last):
            ...
            ValueError: unknown rounding mode 6
            sage: RealField(100, rnd=10^100)
            Traceback (most recent call last):
            ...
            OverflowError: Sage Integer too large to convert to C long

        Check methods inherited from categories::

            sage: RealField(10).is_finite()
            False
        """
        if prec < MPFR_PREC_MIN or prec > MPFR_PREC_MAX:
            raise ValueError("prec (=%s) must be >= %s and <= %s" % (
                prec, MPFR_PREC_MIN, MPFR_PREC_MAX))
        self._prec = prec
        self.sci_not = sci_not

        self.rnd = <mpfr_rnd_t>rnd
        cdef const char* rnd_str = mpfr_print_rnd_mode(self.rnd)
        if rnd_str is NULL:
            raise ValueError("unknown rounding mode {}".format(rnd))
        self.rnd_str = char_to_str(rnd_str + 5)  # Strip "MPFR_"

        from sage.categories.fields import Fields
        Parent.__init__(self, self, names=tuple(), normalize=False,
                        category=Fields().Infinite().Metric().Complete())

        # Initialize zero and one
        cdef RealNumber rn
        rn = self._new()
        mpfr_set_zero(rn.value, 1)
        self._zero_element = rn

        rn = self._new()
        mpfr_set_ui(rn.value, 1, MPFR_RNDZ)
        self._one_element = rn

        self._populate_coercion_lists_(convert_method_name='_mpfr_')

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: RealField() # indirect doctest
            Real Field with 53 bits of precision
            sage: RealField(100000) # indirect doctest
            Real Field with 100000 bits of precision
            sage: RealField(17,rnd='RNDD') # indirect doctest
            Real Field with 17 bits of precision and rounding RNDD
        """
        s = "Real Field with %s bits of precision" % self._prec
        if self.rnd != MPFR_RNDN:
            s = s + " and rounding %s" % (self.rnd_str)
        return s

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(RealField()) # indirect doctest
            \Bold{R}
        """
        return "\\Bold{R}"

    def _sage_input_(self, sib, coerce):
        r"""
        Produce an expression which will reproduce this value when
        evaluated.

        EXAMPLES::

            sage: sage_input(RR, verify=True)
            # Verified
            RR
            sage: sage_input(RealField(25, rnd='RNDZ'), verify=True)
            # Verified
            RealField(25, rnd='RNDZ')
            sage: k = (RR, RealField(75, rnd='RNDU'), RealField(13))
            sage: sage_input(k, verify=True)
            # Verified
            (RR, RealField(75, rnd='RNDU'), RealField(13))
            sage: sage_input((k, k), verify=True)
            # Verified
            RR75u = RealField(75, rnd='RNDU')
            RR13 = RealField(13)
            ((RR, RR75u, RR13), (RR, RR75u, RR13))
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: RealField(99, rnd='RNDD')._sage_input_(SageInputBuilder(), False)
            {call: {atomic:RealField}({atomic:99}, rnd={atomic:'RNDD'})}
        """
        if self.rnd == MPFR_RNDN and self.prec() == 53:
            return sib.name('RR')

        if self.rnd != MPFR_RNDN:
            rnd_abbrev = self.rnd_str[-1:].lower()
            v = sib.name('RealField')(sib.int(self.prec()), rnd=self.rnd_str)
        else:
            rnd_abbrev = ''
            v = sib.name('RealField')(sib.int(self.prec()))

        name = 'RR%d%s' % (self.prec(), rnd_abbrev)
        sib.cache(self, v, name)
        return v

    cpdef bint is_exact(self) except -2:
        """
        Return ``False``, since a real field (represented using finite
        precision) is not exact.

        EXAMPLES::

            sage: RR.is_exact()
            False
            sage: RealField(100).is_exact()
            False
        """
        return False

    def _element_constructor_(self, x, base=10):
        """
        Coerce ``x`` into this real field.

        EXAMPLES::

            sage: R = RealField(20)
            sage: R('1.234')
            1.2340
            sage: R('2', base=2)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert '2' to a real number
            sage: a = R('1.1001', base=2); a
            1.5625
            sage: a.str(2)
            '1.1001000000000000000'
            sage: R(oo)
            +infinity
            sage: R(unsigned_infinity)
            Traceback (most recent call last):
            ...
            ValueError: can only convert signed infinity to RR
            sage: R(CIF(NaN))                                                           # needs sage.symbolic
            NaN
            sage: R(complex(1.7))
            1.7000
        """
        if hasattr(x, '_mpfr_'):
            return x._mpfr_(self)
        cdef RealNumber z
        z = self._new()
        z._set(x, base)
        return z

    cpdef _coerce_map_from_(self, S):
        """
        Canonical coercion of x to this MPFR real field.

        The rings that canonically coerce to this MPFR real field are:

        - Any MPFR real field with precision that is as large as this one

        - int, integer, and rational rings.

        - the field of algebraic reals

        - floats and RDF if self.prec = 53

        EXAMPLES::

            sage: RR.has_coerce_map_from(ZZ) # indirect doctest
            True
            sage: RR.has_coerce_map_from(float)
            True
            sage: RealField(100).has_coerce_map_from(float)
            False
            sage: RR.has_coerce_map_from(RealField(200))
            True
            sage: RR.has_coerce_map_from(RealField(20))
            False
            sage: RR.has_coerce_map_from(RDF)
            True
            sage: RR.coerce_map_from(ZZ)(2)
            2.00000000000000
            sage: RR.coerce(3.4r)
            3.40000000000000
            sage: RR.coerce(3.4)
            3.40000000000000
            sage: RR.coerce(3.4r)
            3.40000000000000
            sage: RR.coerce(3.400000000000000000000000000000000000000000)
            3.40000000000000
            sage: RealField(100).coerce(3.4)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Real Field with 53 bits of precision to Real Field with 100 bits of precision
            sage: RR.coerce(17/5)
            3.40000000000000
            sage: RR.coerce(2^4000)
            1.31820409343094e1204
            sage: RR.coerce_map_from(float)
            Generic map:
              From: Set of Python objects of class 'float'
              To:   Real Field with 53 bits of precision

        TESTS::

            sage: 1.0 - ZZ(1) - int(1) - 1 - QQ(1) - RealField(100)(1) - AA(1) - RLF(1)             # needs sage.rings.number_field
            -6.00000000000000
            sage: R = RR['x']   # Hold reference to avoid garbage collection, see Issue #24709
            sage: R.get_action(ZZ)
            Right scalar multiplication by Integer Ring on Univariate Polynomial Ring in x over Real Field with 53 bits of precision
        """
        from sage.rings.integer_ring import ZZ
        from sage.rings.rational_field import QQ
        from sage.rings.real_double import RDF

        if S is ZZ:
            return ZZtoRR(ZZ, self)
        elif S is QQ:
            return QQtoRR(QQ, self)
        elif (S is RDF or S is float) and self._prec <= 53:
            return double_toRR(S, self)
        elif S is long:
            return int_toRR(long, self)
        elif S is int:
            return int_toRR(int, self)
        elif isinstance(S, RealField_class) and S.prec() >= self._prec:
            return RRtoRR(S, self)
        elif QQ.has_coerce_map_from(S):
            return QQtoRR(QQ, self) * QQ._internal_coerce_map_from(S)
        from sage.rings.real_lazy import RLF
        if isinstance(S, sage.rings.abc.AlgebraicRealField) or S is RLF:
            return self.convert_method_map(S, "_mpfr_")
        return self._coerce_map_via([RLF], S)

    def __richcmp__(RealField_class self, other, int op):
        """
        Compare two real fields, returning ``True`` if they are equivalent
        and ``False`` if they are not.

        EXAMPLES::

            sage: RealField(10) == RealField(11)
            False
            sage: RealField(10) == RealField(10)
            True
            sage: RealField(10,rnd='RNDN') == RealField(10,rnd='RNDZ')
            False

        Scientific notation affects only printing, not mathematically how
        the field works, so this does not affect equality testing::

            sage: RealField(10,sci_not=True) == RealField(10,sci_not=False)
            True
            sage: RealField(10) == IntegerRing()
            False

        ::

            sage: RS = RealField(sci_not=True)
            sage: RR == RS
            True
            sage: RS.scientific_notation(False)
            sage: RR == RS
            True
        """
        if op != Py_EQ and op != Py_NE:
            return NotImplemented
        if not isinstance(other, RealField_class):
            return NotImplemented

        _other = <RealField_class>other  # to access C structure
        return (self._prec == _other._prec and
                self.rnd == _other.rnd) == (op == Py_EQ)

    def __reduce__(self):
        """
        Return the arguments sufficient for pickling.

        EXAMPLES::

            sage: R = RealField(sci_not=1, prec=200, rnd='RNDU')
            sage: loads(dumps(R)) == R
            True
        """
        return __create__RealField_version0, (self._prec, self.sci_not, self.rnd_str)

    def construction(self):
        r"""
        Return the functorial construction of ``self``, namely,
        completion of the rational numbers with respect to the prime
        at `\infty`.

        Also preserves other information that makes this field unique (e.g.
        precision, rounding, print mode).

        EXAMPLES::

            sage: R = RealField(100, rnd='RNDU')
            sage: c, S = R.construction(); S
            Rational Field
            sage: R == c(S)
            True
        """
        from sage.categories.pushout import CompletionFunctor
        return (CompletionFunctor(sage.rings.infinity.Infinity,
                                  self.prec(),
                                  {'type': 'MPFR',
                                   'sci_not': self.scientific_notation(),
                                   'rnd': self.rnd}),
                sage.rings.rational_field.QQ)

    def gen(self, i=0):
        """
        Return the ``i``-th generator of ``self``.

        EXAMPLES::

            sage: R=RealField(100)
            sage: R.gen(0)
            1.0000000000000000000000000000
            sage: R.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: self has only one generator
        """
        if i == 0:
            return self(1)
        else:
            raise IndexError("self has only one generator")

    def complex_field(self):
        """
        Return complex field of the same precision.

        EXAMPLES::

            sage: RR.complex_field()
            Complex Field with 53 bits of precision
            sage: RR.complex_field() is CC
            True
            sage: RealField(100,rnd='RNDD').complex_field()
            Complex Field with 100 bits of precision
            sage: RealField(100).complex_field()
            Complex Field with 100 bits of precision
        """
        from sage.rings.complex_mpfr import ComplexField
        return ComplexField(self.prec())

    def algebraic_closure(self):
        """
        Return the algebraic closure of ``self``, i.e., the complex field with
        the same precision.

        EXAMPLES::

            sage: RR.algebraic_closure()
            Complex Field with 53 bits of precision
            sage: RR.algebraic_closure() is CC
            True
            sage: RealField(100,rnd='RNDD').algebraic_closure()
            Complex Field with 100 bits of precision
            sage: RealField(100).algebraic_closure()
            Complex Field with 100 bits of precision
        """
        return self.complex_field()

    def ngens(self):
        """
        Return the number of generators.

        EXAMPLES::

            sage: RR.ngens()
            1
        """
        return 1

    def gens(self):
        """
        Return a list of generators.

        EXAMPLES::

            sage: RR.gens()
            [1.00000000000000]
        """
        return [self.gen()]

    def _is_valid_homomorphism_(self, codomain, im_gens, base_map=None):
        """
        Return ``True`` if the map from ``self`` to ``codomain`` sending
        ``self(1)`` to the unique element of ``im_gens`` is a valid field
        homomorphism. Otherwise, return ``False``.

        EXAMPLES::

            sage: RR._is_valid_homomorphism_(RDF,[RDF(1)])
            True
            sage: RR._is_valid_homomorphism_(CDF,[CDF(1)])
            True
            sage: RR._is_valid_homomorphism_(CDF,[CDF(-1)])
            False
            sage: R=RealField(100)
            sage: RR._is_valid_homomorphism_(R,[R(1)])
            False
            sage: RR._is_valid_homomorphism_(CC,[CC(1)])
            True
            sage: RR._is_valid_homomorphism_(GF(2),GF(2)(1))
            False
        """

        try:
            s = codomain.coerce(self(1))
        except TypeError:
            return False
        return s == im_gens[0]

    def _repr_option(self, key):
        """
        Metadata about the :meth:`_repr_` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: RealField(10)._repr_option('element_is_atomic')
            True
        """
        if key == 'element_is_atomic':
            return True
        return super()._repr_option(key)

    def characteristic(self):
        """
        Return 0, since the field of real numbers has characteristic 0.

        EXAMPLES::

            sage: RealField(10).characteristic()
            0
        """
        return Integer(0)

    def name(self):
        """
        Return the name of ``self``, which encodes the precision and
        rounding convention.

        EXAMPLES::

            sage: RR.name()
            'RealField53_0'
            sage: RealField(100,rnd='RNDU').name()
            'RealField100_2'
        """
        return "RealField%s_%s" % (self._prec, self.rnd)

    def __hash__(self):
        """
        Return a hash function of the field, which takes into account
        the precision and rounding convention.

        EXAMPLES::

            sage: hash(RealField(100,rnd='RNDU')) == hash(RealField(100,rnd='RNDU'))
            True
            sage: hash(RR) == hash(RealField(53))
            True
        """
        return hash(self.name())

    def precision(self):
        """
        Return the precision of ``self``.

        EXAMPLES::

            sage: RR.precision()
            53
            sage: RealField(20).precision()
            20
        """
        return Integer(self._prec)

    prec = precision  # an alias

    def _magma_init_(self, magma):
        r"""
        Return a string representation of ``self`` in the Magma language.

        .. WARNING::

            This ignores the rounding convention of ``self``.

        EXAMPLES::

            sage: magma(RealField(70)) # optional - magma # indirect doctest
            Real field of precision 21
            sage: 10^21 < 2^70 < 10^22
            True
            sage: s = magma(RealField(70)).sage(); s # optional - magma # indirect doctest
            Real Field with 70 bits of precision
        """
        return "RealField(%s : Bits := true)" % self.prec()

    def to_prec(self, prec):
        """
        Return the real field that is identical to ``self``, except that
        it has the specified precision.

        EXAMPLES::

            sage: RR.to_prec(212)
            Real Field with 212 bits of precision
            sage: R = RealField(30, rnd="RNDZ")
            sage: R.to_prec(300)
            Real Field with 300 bits of precision and rounding RNDZ
        """
        if prec == self._prec:
            return self
        else:
            return RealField(prec, self.sci_not, self.rnd)

    def pi(self):
        r"""
        Return `\pi` to the precision of this field.

        EXAMPLES::

            sage: R = RealField(100)
            sage: R.pi()
            3.1415926535897932384626433833
            sage: R.pi().sqrt()/2
            0.88622692545275801364908374167
            sage: R = RealField(150)
            sage: R.pi().sqrt()/2
            0.88622692545275801364908374167057259139877473
        """
        cdef RealNumber x = self._new()
        if self._prec > SIG_PREC_THRESHOLD:
            sig_on()
        # The docs for mpfr_free_cache say "Free the cache used by
        # the functions computing constants if needed (currently
        # mpfr_const_log2, mpfr_const_pi and mpfr_const_euler)", so
        # this isn't a seriously bad thing to do.  This prevents trac
        # #5689.  This is needed for all constants, despite what the docs say.
        # NOTE: The MPFR docs at this time leave off several mpfr_const
        # functions, but this free is needed for them too!
        mpfr_free_cache()
        mpfr_const_pi(x.value, self.rnd)
        if self._prec > SIG_PREC_THRESHOLD:
            sig_off()
        return x

    def euler_constant(self):
        """
        Return Euler's gamma constant to the precision of this field.

        EXAMPLES::

            sage: RealField(100).euler_constant()
            0.57721566490153286060651209008
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_free_cache()
        mpfr_const_euler(x.value, self.rnd)
        sig_off()
        return x

    def catalan_constant(self):
        """
        Return Catalan's constant to the precision of this field.

        EXAMPLES::

            sage: RealField(100).catalan_constant()
            0.91596559417721901505460351493
        """
        cdef RealNumber x = self._new()
        if self._prec > SIG_PREC_THRESHOLD:
            sig_on()
        mpfr_free_cache()
        mpfr_const_catalan(x.value, self.rnd)
        if self._prec > SIG_PREC_THRESHOLD:
            sig_off()
        return x

    def log2(self):
        r"""
        Return `\log(2)` (i.e., the natural log of 2) to the precision
        of this field.

        EXAMPLES::

            sage: R=RealField(100)
            sage: R.log2()
            0.69314718055994530941723212146
            sage: R(2).log()
            0.69314718055994530941723212146
        """
        cdef RealNumber x = self._new()
        if self._prec > SIG_PREC_THRESHOLD:
            sig_on()
        mpfr_free_cache()
        mpfr_const_log2(x.value, self.rnd)
        if self._prec > SIG_PREC_THRESHOLD:
            sig_off()
        return x

    def random_element(self, min=-1, max=1, distribution=None):
        """
        Return a uniformly distributed random number between ``min`` and
        ``max`` (default -1 to 1).

        .. WARNING::

            The argument ``distribution`` is ignored---the random number
            is from the uniform distribution.

        EXAMPLES::

            sage: r = RealField(100).random_element(-5, 10)
            sage: r.parent() is RealField(100)
            True
            sage: -5 <= r <= 10
            True

        TESTS::

            sage: RealField(31).random_element().parent()
            Real Field with 31 bits of precision
            sage: RealField(32).random_element().parent()
            Real Field with 32 bits of precision
            sage: RealField(33).random_element().parent()
            Real Field with 33 bits of precision
            sage: RealField(63).random_element().parent()
            Real Field with 63 bits of precision
            sage: RealField(64).random_element().parent()
            Real Field with 64 bits of precision
            sage: RealField(65).random_element().parent()
            Real Field with 65 bits of precision
            sage: RealField(10).random_element().parent()
            Real Field with 10 bits of precision
            sage: RR.random_element().parent()
            Real Field with 53 bits of precision
        """
        cdef RealNumber x = self._new()
        cdef randstate rstate = current_randstate()
        mpfr_urandomb(x.value, rstate.gmp_state)
        if min == 0 and max == 1:
            return x
        else:
            return (max-min)*x + min

    def factorial(self, int n):
        """
        Return the factorial of the integer ``n`` as a real number.

        EXAMPLES::

            sage: RR.factorial(0)
            1.00000000000000
            sage: RR.factorial(1000000)
            8.26393168833124e5565708
            sage: RR.factorial(-1)
            Traceback (most recent call last):
            ...
            ArithmeticError: n must be nonnegative
        """
        cdef RealNumber x
        if n < 0:
            raise ArithmeticError("n must be nonnegative")
        x = self._new()
        if self._prec > SIG_PREC_THRESHOLD and n < SIG_PREC_THRESHOLD:
            sig_on()
        mpfr_fac_ui(x.value, n, self.rnd)
        if self._prec > SIG_PREC_THRESHOLD and n < SIG_PREC_THRESHOLD:
            sig_off()
        return x

    def rounding_mode(self):
        """
        Return the rounding mode.

        EXAMPLES::

            sage: RR.rounding_mode()
            'RNDN'
            sage: RealField(20,rnd='RNDZ').rounding_mode()
            'RNDZ'
            sage: RealField(20,rnd='RNDU').rounding_mode()
            'RNDU'
            sage: RealField(20,rnd='RNDD').rounding_mode()
            'RNDD'
        """
        return self.rnd_str

    def scientific_notation(self, status=None):
        """
        Set or return the scientific notation printing flag. If this flag
        is ``True`` then real numbers with this space as parent print using
        scientific notation.

        INPUT:

        - ``status`` -- boolean optional flag

        EXAMPLES::

            sage: RR.scientific_notation()
            False
            sage: elt = RR(0.2512); elt
            0.251200000000000
            sage: RR.scientific_notation(True)
            sage: elt
            2.51200000000000e-1
            sage: RR.scientific_notation()
            True
            sage: RR.scientific_notation(False)
            sage: elt
            0.251200000000000
            sage: R = RealField(20, sci_not=1)
            sage: R.scientific_notation()
            True
            sage: R(0.2512)
            2.5120e-1
        """
        if status is None:
            return self.sci_not
        else:
            self.sci_not = status

    def zeta(self, n=2):
        """
        Return an `n`-th root of unity in the real field, if one
        exists, or raise a :exc:`ValueError` otherwise.

        EXAMPLES::

            sage: R = RealField()
            sage: R.zeta()
            -1.00000000000000
            sage: R.zeta(1)
            1.00000000000000
            sage: R.zeta(5)
            Traceback (most recent call last):
            ...
            ValueError: No 5th root of unity in self
        """
        if n == 1:
            return self(1)
        elif n == 2:
            return self(-1)
        raise ValueError("No %sth root of unity in self" % n)

    def _factor_univariate_polynomial(self, f):
        """
        Factor the univariate polynomial ``f``.

        INPUT:

        - ``f`` -- a univariate polynomial defined over the real numbers

        OUTPUT:

        - A factorization of ``f`` over the real numbers into a unit and monic
          irreducible factors

        .. NOTE::

            This is a helper method for
            :meth:`sage.rings.polynomial.polynomial_element.Polynomial.factor`.

            This method calls PARI to compute the factorization.

        TESTS::

            sage: # needs sage.libs.pari
            sage: k = RealField(100)
            sage: R.<x> = k[]
            sage: k._factor_univariate_polynomial( x )
            x
            sage: k._factor_univariate_polynomial( 2*x )
            (2.0000000000000000000000000000) * x
            sage: k._factor_univariate_polynomial( x^2 )
            x^2
            sage: k._factor_univariate_polynomial( x^2 + 1 )
            x^2 + 1.0000000000000000000000000000
            sage: k._factor_univariate_polynomial( x^2 - 1 )
            (x - 1.0000000000000000000000000000) * (x + 1.0000000000000000000000000000)
            sage: k._factor_univariate_polynomial( (x - 1)^3 )
            (x - 1.0000000000000000000000000000)^3
            sage: k._factor_univariate_polynomial( x^2 - 3 )
            (x - 1.7320508075688772935274463415) * (x + 1.7320508075688772935274463415)
        """
        R = f.parent()
        F = list(f._pari_with_name().factor())

        from sage.structure.factorization import Factorization
        return Factorization([(R(g).monic(), e) for g, e in zip(*F)],
                             f.leading_coefficient())

# ****************************************************************************
#
#     RealNumber -- element of Real Field
#
# ****************************************************************************

cdef class RealLiteral(RealNumber)

cdef class RealNumber(sage.structure.element.RingElement):
    """
    A floating point approximation to a real number using any specified
    precision. Answers derived from calculations with such
    approximations may differ from what they would be if those
    calculations were performed with true real numbers. This is due to
    the rounding errors inherent to finite precision calculations.

    The approximation is printed to slightly fewer digits than its
    internal precision, in order to avoid confusing roundoff issues
    that occur because numbers are stored internally in binary.
    """
    def __cinit__(self, parent, x=None, base=None):
        """
        Initialize the parent of this element and allocate memory.

        TESTS::

            sage: from sage.rings.real_mpfr import RealNumber
            sage: RealNumber.__new__(RealNumber, None)
            Traceback (most recent call last):
            ...
            TypeError: Cannot convert NoneType to sage.rings.real_mpfr.RealField_class
            sage: RealNumber.__new__(RealNumber, ZZ)
            Traceback (most recent call last):
            ...
            TypeError: Cannot convert sage.rings.integer_ring.IntegerRing_class to sage.rings.real_mpfr.RealField_class
            sage: RealNumber.__new__(RealNumber, RR)
            NaN
        """
        cdef RealField_class p = <RealField_class?>parent
        mpfr_init2(self.value, p._prec)
        self._parent = p

    def __init__(self, parent, x=0, int base=10):
        """
        Create a real number. Should be called by first creating a
        RealField, as illustrated in the examples.

        EXAMPLES::

            sage: R = RealField()
            sage: R('1.2456')
            1.24560000000000
            sage: R = RealField(3)
            sage: R('1.2456')
            1.2

        EXAMPLES: Rounding Modes

        ::

            sage: w = RealField(3)(5/2)
            sage: RealField(2, rnd="RNDZ")(w).str(2)
            '10.'
            sage: RealField(2, rnd="RNDD")(w).str(2)
            '10.'
            sage: RealField(2, rnd="RNDU")(w).str(2)
            '11.'
            sage: RealField(2, rnd="RNDN")(w).str(2)
            '10.'

        Conversion from gmpy2 numbers::

            sage: from gmpy2 import *
            sage: RR(mpz(5))
            5.00000000000000
            sage: RR(mpq(1/2))
            0.500000000000000
            sage: RR(mpfr('42.1'))
            42.1000000000000

        .. NOTE::

           A real number is an arbitrary precision mantissa with a
           limited precision exponent. A real number can have three
           special values: Not-a-Number (NaN) or plus or minus
           Infinity. NaN represents an uninitialized object, the
           result of an invalid operation (like 0 divided by 0), or a
           value that cannot be determined (like +Infinity minus
           +Infinity). Moreover, like in the IEEE 754-1985 standard,
           zero is signed, i.e. there are both +0 and -0; the behavior
           is the same as in the IEEE 754-1985 standard and it is
           generalized to the other functions supported by MPFR.

        TESTS::

            sage: TestSuite(R).run()

        Test underscores as digit separators (PEP 515,
        https://www.python.org/dev/peps/pep-0515/)::

            sage: RealNumber('1_3.1e-32_45')
            1.31000000000000e-3244

        Test conversion from base different from `10`::

            sage: RR('0xabc')
            Traceback (most recent call last):
            ...
            TypeError: unable to convert '0xabc' to a real number
            sage: RR("0x123.e1", base=0)  # rel tol 1e-12
            291.878906250000
            sage: RR("0x123.@1", base=0)  # rel tol 1e-12
            4656.00000000000
            sage: RR("1Xx", base=36)  # rel tol 1e-12
            2517.00000000000
            sage: RR("-1Xx@-10", base=62)  # rel tol 1e-12
            -7.08805492048139e-15
            sage: RR("1", base=1)
            Traceback (most recent call last):
            ...
            ValueError: base (=1) must be 0 or between 2 and 62
            sage: RR("1", base=-1)
            Traceback (most recent call last):
            ...
            ValueError: base (=-1) must be 0 or between 2 and 62
            sage: RR("1", base=63)
            Traceback (most recent call last):
            ...
            ValueError: base (=63) must be 0 or between 2 and 62
        """
        if x is not None:
            self._set(x, base)

    def _magma_init_(self, magma):
        r"""
        Return a string representation of ``self`` in the Magma language.

        EXAMPLES::

            sage: magma(RR(10.5)) # indirect doctest, optional - magma
            10.5000000000000
            sage: magma(RealField(200)(1/3)) # indirect, optional - magma
            0.333333333333333333333333333333333333333333333333333333333333
            sage: magma(RealField(1000)(1/3)) # indirect, optional - magma
            0.3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
        """
        real_string = self.str()
        digit_precision_upper_bound = len(real_string)
        return "%s!%sp%s" % (self.parent()._magma_init_(magma),
                             real_string, digit_precision_upper_bound)

    @property
    def __array_interface__(self):
        """
        Used for NumPy conversion.

        EXAMPLES::

            sage: # needs numpy
            sage: import numpy
            sage: numpy.arange(10.0)
            array([0., 1., 2., 3., 4., 5., 6., 7., 8., 9.])
            sage: numpy.array([1.0, 1.1, 1.2]).dtype
            dtype('float64')
            sage: numpy.array([1.000000000000000000000000000000000000]).dtype
            dtype('O')
        """
        if (<RealField_class>self._parent)._prec <= 53:
            return numpy_double_interface
        else:
            return numpy_object_interface

    cdef _set(self, x, int base):
        # This should not be called except when the number is being created.
        # Real Numbers are supposed to be immutable.
        cdef RealField_class parent
        parent = self._parent
        if base != 0 and (base < 2 or base > 62):
            raise ValueError(f"base (={base}) must be 0 or between 2 and 62")
        if isinstance(x, RealNumber):
            if isinstance(x, RealLiteral):
                s = (<RealLiteral>x).literal
                base = (<RealLiteral>x).base
                if mpfr_set_str(self.value, str_to_bytes(s), base, parent.rnd):
                    self._set(s, base)
            else:
                mpfr_set(self.value, (<RealNumber>x).value, parent.rnd)
        elif isinstance(x, Integer):
            mpfr_set_z(self.value, (<Integer>x).value, parent.rnd)
        elif isinstance(x, Rational):
            mpfr_set_q(self.value, (<Rational>x).value, parent.rnd)
        elif isinstance(x, Gen) and set_real_mpfr_element_from_gen(self, x):
            pass
        elif isinstance(x, int):
            x = Integer(x)
            mpfr_set_z(self.value, (<Integer>x).value, parent.rnd)
        elif isinstance(x, float):
            mpfr_set_d(self.value, x, parent.rnd)
        elif isinstance(x, complex) and x.imag == 0:
            mpfr_set_d(self.value, x.real, parent.rnd)
        elif isinstance(x, RealDoubleElement):
            mpfr_set_d(self.value, (<RealDoubleElement>x)._value, parent.rnd)
        elif type(x) is gmpy2.mpfr:
            mpfr_set(self.value, (<gmpy2.mpfr>x).f, parent.rnd)
        elif type(x) is gmpy2.mpq:
            mpfr_set_q(self.value, (<gmpy2.mpq>x).q, parent.rnd)
        elif type(x) is gmpy2.mpz:
            mpfr_set_z(self.value, (<gmpy2.mpz>x).z, parent.rnd)
        else:
            s = str(x).replace(' ', '').replace('_', '')
            s_lower = s.lower()
            if s_lower == 'infinity':
                raise ValueError('can only convert signed infinity to RR')
            elif mpfr_set_str(self.value, str_to_bytes(s), base, parent.rnd) == 0:
                pass
            elif s == 'NaN' or s == '@NaN@' or s == '[..NaN..]' or s == 'NaN+NaN*I':
                mpfr_set_nan(self.value)
            elif s_lower == '+infinity':
                mpfr_set_inf(self.value, 1)
            elif s_lower == '-infinity':
                mpfr_set_inf(self.value, -1)
            else:
                raise TypeError("unable to convert {!r} to a real number".format(s))

    def __reduce__(self):
        """
        EXAMPLES::

            sage: R = RealField(sci_not=1, prec=200, rnd='RNDU')
            sage: b = R('393.39203845902384098234098230948209384028340')
            sage: loads(dumps(b)) == b
            True
            sage: b = R(1)/R(0); b
            +infinity
            sage: loads(dumps(b)) == b
            True
            sage: b = R(-1)/R(0); b
            -infinity
            sage: loads(dumps(b)) == b
            True
            sage: b = R(-1).sqrt(); b
            1.0000000000000000000000000000000000000000000000000000000000*I
            sage: loads(dumps(b)) == b
            True
        """
        s = self.str(32, no_sci=False, e='@')
        return (__create__RealNumber_version0, (self._parent, s, 32))

    def __dealloc__(self):
        if self._parent is not None:
            mpfr_clear(self.value)

    def __repr__(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: RR(2.1) # indirect doctest
            2.10000000000000
        """
        return self.str(truncate=True)

    def __format__(self, format_spec):
        """
        Return a formatted string representation of this real number.

        EXAMPLES::

            sage: format(RR(32/3), '.4f')
            '10.6667'
            sage: '{:.4e}'.format(RR(2/3))
            '6.6667e-1'
            sage: format(RealField(240)(1/7), '.60f')
            '0.142857142857142857142857142857142857142857142857142857142857'

        TESTS::

            sage: format(RR(oo), '.4')
            'Infinity'
            sage: format(RR(-oo), '.4')
            '-Infinity'
            sage: format(RR(NaN), '.4')                                                 # needs sage.symbolic
            'NaN'
            sage: '{}'.format(RR(oo))
            '+infinity'
        """
        if not format_spec:
            # Since there are small formatting differences between RR and
            # Decimal, we avoid converting to Decimal when spec is empty, which
            # commonly occurs in f-strings. This could be improved by fully
            # parsing format_spec in order to avoid using Decimal altogether.
            return repr(self)
        from decimal import Decimal
        return format(Decimal(repr(self)), format_spec)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(RR(2.1)) # indirect doctest
            2.10000000000000
            sage: latex(RR(2e100)) # indirect doctest
            2.00000000000000 \times 10^{100}
        """
        s = repr(self)
        parts = s.split('e')
        if len(parts) > 1:
            # scientific notation
            if parts[1][0] == '+':
                parts[1] = parts[1][1:]
            s = "%s \\times 10^{%s}" % (parts[0], parts[1])
        return s

    def _interface_init_(self, I=None):
        """
        Return string representation of ``self`` in base 10, avoiding
        scientific notation except for very large or very small numbers.

        This is most likely to make sense in other computer algebra systems
        (this function is the default for exporting to other computer
        algebra systems).

        EXAMPLES::

            sage: n = 1.3939494594
            sage: n._interface_init_()
            '1.3939494593999999'
            sage: s1 = RR(sin(1)); s1                                                   # needs sage.symbolic
            0.841470984807897
            sage: s1._interface_init_()                                                 # needs sage.symbolic
            '0.84147098480789650'
            sage: s1 == RR(gp(s1))                                                      # needs sage.symbolic
            True
        """
        return self.str(10, no_sci=True)

    def _mathematica_init_(self):
        """
        TESTS:

        Check that :issue:`28814` is fixed::

            sage: mathematica(3.5e-15)           # optional - mathematica
            3.5*^-15
            sage: mathematica.Log(3.5e6).sage()  # optional - mathematica
            15.06827352645964
        """
        return self.str(10, e='*^')

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: values = [-infinity, -20, 0, 1, 2^500, -2^4000, -2^-500, 2^-4000]
            sage: values += [NaN, -e]                                                   # needs sage.symbolic
            sage: for prec in (2, 53, 200):
            ....:     for rnd_dir in ('RNDN', 'RNDD', 'RNDU', 'RNDZ'):
            ....:         fld = RealField(prec, rnd=rnd_dir)
            ....:         var = polygen(fld)
            ....:         for v in [NaN, -infinity, -20, -e, 0, 1, 2^500, -2^4000, -2^-500, 2^-4000] + [fld.random_element() for _ in range(5)]:
            ....:             for preparse in (True, False, None):
            ....:                 _ = sage_input(fld(v), verify=True, preparse=preparse)
            ....:                 _ = sage_input(fld(v) * var, verify=True, preparse=preparse)
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sib_np = SageInputBuilder(preparse=False)
            sage: RR60 = RealField(60)
            sage: RR(-infinity)._sage_input_(sib, True)
            {unop:- {call: {atomic:RR}({atomic:Infinity})}}
            sage: RR(NaN)._sage_input_(sib, True)                                       # needs sage.symbolic
            {call: {atomic:RR}({atomic:NaN})}
            sage: RR(12345)._sage_input_(sib, True)
            {atomic:12345}
            sage: RR(-12345)._sage_input_(sib, False)
            {unop:- {call: {atomic:RR}({atomic:12345})}}
            sage: RR(1.579)._sage_input_(sib, True)
            {atomic:1.579}
            sage: RR(1.579)._sage_input_(sib_np, True)
            {atomic:1.579}
            sage: RR60(1.579)._sage_input_(sib, True)
            {atomic:1.5790000000000000008}
            sage: RR60(1.579)._sage_input_(sib_np, True)
            {call: {call: {atomic:RealField}({atomic:60})}({atomic:'1.5790000000000000008'})}
            sage: RR(1.579)._sage_input_(sib_np, False)
            {call: {atomic:RR}({atomic:1.579})}
            sage: RR(1.579)._sage_input_(sib_np, 2)
            {atomic:1.579}
            sage: RealField(150)(pi)._sage_input_(sib, True)                            # needs sage.symbolic
            {atomic:3.1415926535897932384626433832795028841971694008}
            sage: RealField(150)(pi)._sage_input_(sib_np, True)                         # needs sage.symbolic
            {call: {call: {atomic:RealField}({atomic:150})}({atomic:'3.1415926535897932384626433832795028841971694008'})}
        """
        # We have a bewildering array of conditions to deal with:
        # 1) number, or NaN or infinity
        # 2) rounding direction: up, down, or nearest
        # 3) preparser enabled: yes, no, or maybe
        # 4) is this number equal to some Python float: yes or no
        # 5) coerced

        # First, handle NaN and infinity
        if not mpfr_number_p(self.value):
            if mpfr_inf_p(self.value):
                v = sib.name('Infinity')
            else:
                v = sib.name('NaN')
            v = sib(self._parent)(v)
            if mpfr_sgn(self.value) < 0:
                v = -v
            return v

        from sage.rings.integer_ring import ZZ

        cdef mpfr_rnd_t rnd = (<RealField_class>self._parent).rnd

        cdef bint negative = mpfr_sgn(self.value) < 0
        if negative:
            self = -self

        # There are five possibilities for printing this floating-point
        # number, ordered from prettiest to ugliest (IMHO).
        # 1) An integer: 42
        # 2) A simple literal: 3.14159
        # 3) A coerced integer: RR(42)
        # 4) A coerced literal: RR(3.14159)
        # 5) A coerced string literal: RR('3.14159')

        # To use choice 1 or choice 3, this number must be an integer.
        cdef bint can_use_int_literal = \
            self.abs() < (Integer(1) << self.prec()) and self in ZZ

        # If "not coerced", then we will introduce a conversion
        # ourselves.  If coerced==2, then there will be an external conversion.
        # On the other hand, if coerced==1 (or True), then we only
        # have a coercion, not a conversion; which means we need to read
        # in a number with at least the number of bits we need.
        will_convert = (coerced == 2 or not coerced)

        self_str = self.str(skip_zeroes=(will_convert or self.prec() <= 53))

        # To use choice 2 or choice 4, we must be able to read
        # numbers of this precision as a literal.  We support this
        # only for the default rounding mode; "pretty" output for
        # other rounding modes is a lot of work for very little gain
        # (since other rounding modes are very rarely used).
        # (One problem is that we don't know whether it will be the
        # positive or negative value that will be coerced into the
        # desired parent; for example, this differs between "x^2 - 1.3*x"
        # and "-1.3*x".)
        cdef bint can_use_float_literal = \
            rnd == MPFR_RNDN and (sib.preparse() or
                                  ((will_convert or self.prec() <= 53) and
                                   self._parent(float(self_str)) == self))

        if can_use_int_literal or can_use_float_literal:
            if can_use_int_literal:
                v = sib.int(self._integer_())
            else:
                v = sib.float_str(self_str)
            if not coerced and (can_use_int_literal or not sib.preparse() or create_RealNumber(self_str).prec() != self.prec()):
                v = sib(self.parent())(v)
        else:
            if rnd == MPFR_RNDN:
                s = self_str
            else:
                # This is tricky.  str() uses mpfr_get_str() with
                # digits=0; this guarantees to give enough digits
                # to recreate the input, if we print and read with
                # round-to-nearest.  However, we are not going to
                # read with round-to-nearest, so we might need more digits.
                # If we're going to read with round-down, then we need
                # to print with round-up; and we'll use one more bit
                # to make sure we have enough digits.
                # Since we always read nonnegative numbers, reading with
                # RNDZ is the same as reading with RNDD.
                if rnd == MPFR_RNDD or rnd == MPFR_RNDZ:
                    fld = RealField(self.prec() + 1, rnd='RNDU')
                else:
                    fld = RealField(self.prec() + 1, rnd='RNDD')
                s = fld(self).str()
            v = sib(self.parent())(sib.float_str(repr(s)))

        if negative:
            v = -v

        return v

    def __hash__(self):
        """
        Return the hash of self, which coincides with the Python float
        (and often int) type.

        This has the drawback that two very close high precision numbers
        will have the same hash, but allows them to play nicely with other
        real types.

        EXAMPLES::

            sage: hash(RR(1.2)) == hash(1.2r)
            True
        """
        return hash(float(self))

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        Return the image of ``self`` under the homomorphism from the rational
        field to ``codomain``.

        This always just returns ``self`` coerced into the ``codomain``.

        EXAMPLES::

            sage: RR(2.1)._im_gens_(RDF, [RDF(1)])
            2.1
            sage: R = RealField(20)
            sage: RR(2.1)._im_gens_(R, [R(1)])
            2.1000
        """
        return codomain(self)  # since 1 |--> 1

    def real(self):
        """
        Return the real part of ``self``.

        (Since ``self`` is a real number, this simply returns ``self``.)

        EXAMPLES::

            sage: RR(2).real()
            2.00000000000000
            sage: RealField(200)(-4.5).real()
            -4.5000000000000000000000000000000000000000000000000000000000
        """
        return self

    def imag(self):
        """
        Return the imaginary part of ``self``.

        (Since ``self`` is a real number, this simply returns exactly 0.)

        EXAMPLES::

            sage: RR.pi().imag()
            0
            sage: RealField(100)(2).imag()
            0
        """
        from sage.rings.integer_ring import ZZ

        return ZZ(0)

    def str(self, int base=10, size_t digits=0, *, no_sci=None,
            e=None, bint truncate=False, bint skip_zeroes=False):
        """
        Return a string representation of ``self``.

        INPUT:

        - ``base`` -- (default: 10) base for output

        - ``digits`` -- (default: 0) number of digits to display; when
          ``digits`` is zero, choose this automatically

        - ``no_sci`` -- if 2, never print using scientific notation; if
          ``True``, use scientific notation only for large or small
          numbers; if ``False`` always print with scientific notation;
          if ``None`` (the default), print how the parent prints.

        - ``e`` -- symbol used in scientific notation; defaults to 'e' for
          base=10, and '@' otherwise

        - ``truncate`` -- boolean (default: ``False``); if ``True``, round off the
          last digits in base-10 printing to lessen confusing base-2
          roundoff issues. This flag may not be used in other bases or
          when ``digits`` is given.

        - ``skip_zeroes`` -- boolean (default: ``False``); if ``True``, skip
          trailing zeroes in mantissa

        EXAMPLES::

            sage: a = 61/3.0; a
            20.3333333333333
            sage: a.str()
            '20.333333333333332'
            sage: a.str(truncate=True)
            '20.3333333333333'
            sage: a.str(2)
            '10100.010101010101010101010101010101010101010101010101'
            sage: a.str(no_sci=False)
            '2.0333333333333332e1'
            sage: a.str(16, no_sci=False)
            '1.4555555555555@1'
            sage: a.str(digits=5)
            '20.333'
            sage: a.str(2, digits=5)
            '10100.'

            sage: b = 2.0^99
            sage: b.str()
            '6.3382530011411470e29'
            sage: b.str(no_sci=False)
            '6.3382530011411470e29'
            sage: b.str(no_sci=True)
            '6.3382530011411470e29'
            sage: c = 2.0^100
            sage: c.str()
            '1.2676506002282294e30'
            sage: c.str(no_sci=False)
            '1.2676506002282294e30'
            sage: c.str(no_sci=True)
            '1.2676506002282294e30'
            sage: c.str(no_sci=2)
            '1267650600228229400000000000000.'
            sage: 0.5^53
            1.11022302462516e-16
            sage: 0.5^54
            5.55111512312578e-17
            sage: (0.01).str()
            '0.010000000000000000'
            sage: (0.01).str(skip_zeroes=True)
            '0.01'
            sage: (-10.042).str()
            '-10.042000000000000'
            sage: (-10.042).str(skip_zeroes=True)
            '-10.042'
            sage: (389.0).str(skip_zeroes=True)
            '389.'

        Test various bases::

            sage: print((65536.0).str(base=2))
            1.0000000000000000000000000000000000000000000000000000e16
            sage: print((65536.0).str(base=36))
            1ekg.00000000
            sage: print((65536.0).str(base=62))
            H32.0000000

        String conversion respects rounding::

            sage: x = -RR.pi()
            sage: x.str(digits=1)
            '-3.'
            sage: y = RealField(53, rnd="RNDD")(x)
            sage: y.str(digits=1)
            '-4.'
            sage: y = RealField(53, rnd="RNDU")(x)
            sage: y.str(digits=1)
            '-3.'
            sage: y = RealField(53, rnd="RNDZ")(x)
            sage: y.str(digits=1)
            '-3.'
            sage: y = RealField(53, rnd="RNDA")(x)
            sage: y.str(digits=1)
            '-4.'

        Zero has the correct number of digits::

            sage: zero = RR.zero()
            sage: print(zero.str(digits=3))
            0.00
            sage: print(zero.str(digits=3, no_sci=False))
            0.00e0
            sage: print(zero.str(digits=3, skip_zeroes=True))
            0.

        The output always contains a decimal point, except when using
        scientific notation with exactly one digit::

            sage: print((1e1).str(digits=1))
            10.
            sage: print((1e10).str(digits=1))
            1e10
            sage: print((1e-1).str(digits=1))
            0.1
            sage: print((1e-10).str(digits=1))
            1e-10
            sage: print((-1e1).str(digits=1))
            -10.
            sage: print((-1e10).str(digits=1))
            -1e10
            sage: print((-1e-1).str(digits=1))
            -0.1
            sage: print((-1e-10).str(digits=1))
            -1e-10

        TESTS::

            sage: x = RR.pi()
            sage: x.str(base=1)
            Traceback (most recent call last):
            ...
            ValueError: base (=1) must be an integer between 2 and 62
            sage: x.str(base=63)
            Traceback (most recent call last):
            ...
            ValueError: base (=63) must be an integer between 2 and 62
            sage: x.str(digits=-10)
            Traceback (most recent call last):
            ...
            OverflowError: can...t convert negative value to size_t
            sage: x.str(base=16, truncate=True)
            Traceback (most recent call last):
            ...
            ValueError: truncate is only supported in base 10
            sage: x.str(digits=10, truncate=True)
            Traceback (most recent call last):
            ...
            ValueError: cannot truncate when digits is given
        """
        if base < 2 or base > 62:
            raise ValueError("base (=%s) must be an integer between 2 and 62" % base)
        if mpfr_nan_p(self.value):
            if base >= 24:
                return "@NaN@"
            else:
                return "NaN"
        elif mpfr_inf_p(self.value):
            if mpfr_sgn(self.value) > 0:
                return "+infinity"
            else:
                return "-infinity"

        if e is None:
            if base > 10:
                e = '@'
            else:
                e = 'e'

        if truncate:
            if base != 10:
                raise ValueError("truncate is only supported in base 10")
            if digits:
                raise ValueError("cannot truncate when digits is given")
            # This computes digits = floor(log_{10}(2^(b-1))),
            # which is the number of *decimal* digits that are
            # "right", given that the last binary bit of the binary
            # number can be off.  That is, if this real is within a
            # relative error of 2^(-b) of an exact decimal with
            # `digits` digits, that decimal will be returned.
            # This is equivalent to saying that exact decimals with
            # `digits` digits differ by at least 2*2^(-b) (relative).

            # (Depending on the precision and the exact number involved,
            # adjacent exact decimals can differ by far more than 2*2^(-b)
            # (relative).)

            # This avoids the confusion a lot of people have with the last
            # 1-2 binary digits being wrong due to rounding coming from
            # representing numbers in binary.
            digits = <size_t>(((<RealField_class>self._parent)._prec - 1) * M_LN2_LN10)
            if digits < 2:
                digits = 2

            # For backwards compatibility, add one extra digit for 0.0
            if mpfr_zero_p(self.value):
                digits += 1

        sig_on()
        cdef char *s
        cdef mp_exp_t exponent
        s = mpfr_get_str(NULL, &exponent, base, digits,
                         self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        if s is NULL:
            raise RuntimeError("unable to convert an mpfr number to a string")
        # t contains just digits (no sign, decimal point or exponent)
        t = char_to_str(s)
        if t[0] == '-':
            sgn = "-"
            t = char_to_str(s + 1)
        else:
            sgn = ""
        mpfr_free_str(s)

        if skip_zeroes:
            t = _re_skip_zeroes.match(t).group(1)

        # Treat 0.0 as having exponent 1, this gives better results
        # (effectively treating it as 0. instead of .0)
        if mpfr_zero_p(self.value):
            exponent = 1

        if no_sci is None:
            use_sci = (<RealField_class>self._parent).sci_not or abs(exponent-1) >= 6
        elif no_sci is True:
            use_sci = abs(exponent-1) >= 6
        else:
            use_sci = not no_sci

        if use_sci:
            if len(t) > 1:
                t = t[0] + "." + t[1:]
            return "%s%s%s%s" % (sgn, t, e, exponent-1)

        if exponent <= 0:
            # Only digits after decimal point
            return sgn + '0.' + '0' * (-exponent) + t
        elif exponent >= len(t):
            # Only digits before decimal point
            return sgn + t + '0' * (exponent - len(t)) + "."
        else:
            # Digits before and after decimal point
            return sgn + t[:exponent] + "." + t[exponent:]

    def hex(self):
        """
        Return a hexadecimal floating-point representation of ``self``, in the
        style of C99 hexadecimal floating-point constants.

        EXAMPLES::

            sage: RR(-1/3).hex()
            '-0x5.5555555555554p-4'
            sage: Reals(100)(123.456e789).hex()
            '0xf.721008e90630c8da88f44dd2p+2624'
            sage: (-0.).hex()
            '-0x0p+0'

        ::

            sage: [(a.hex(), float(a).hex()) for a in [.5, 1., 2., 16.]]
            [('0x8p-4', '0x1.0000000000000p-1'),
            ('0x1p+0', '0x1.0000000000000p+0'),
            ('0x2p+0', '0x1.0000000000000p+1'),
            ('0x1p+4', '0x1.0000000000000p+4')]

        Special values::

            sage: [RR(s).hex() for s in ['+inf', '-inf', 'nan']]
            ['inf', '-inf', 'nan']
        """
        cdef char *s
        cdef int r
        sig_on()
        r = mpfr_asprintf(&s, "%Ra", self.value)
        sig_off()
        if r < 0:  # MPFR free()s its buffer itself in this case
            raise RuntimeError("unable to convert an mpfr number to a string")
        t = char_to_str(s)
        mpfr_free_str(s)
        return t

    def __copy__(self):
        """
        Return copy of ``self`` - since ``self`` is immutable, we just return
        ``self`` again.

        EXAMPLES::

            sage: a = 3.5
            sage: copy(a) is a
            True
        """
        return self    # since object is immutable.

    def __deepcopy__(self, memo):
        """
        EXAMPLES::

            sage: a = 3.5
            sage: deepcopy(a) is a
            True
        """
        return self    # since object is immutable.

    def _integer_(self, Z=None):
        """
        If this floating-point number is actually an integer, return that
        integer. Otherwise, raise an exception.

        EXAMPLES::

            sage: ZZ(237.0) # indirect doctest
            237
            sage: ZZ(0.0/0.0)
            Traceback (most recent call last):
            ...
            TypeError: Attempt to coerce non-integral RealNumber to Integer
            sage: ZZ(1.0/0.0)
            Traceback (most recent call last):
            ...
            TypeError: Attempt to coerce non-integral RealNumber to Integer
            sage: ZZ(-123456789.0)
            -123456789
            sage: ZZ(RealField(300)(2.0)^290)
            1989292945639146568621528992587283360401824603189390869761855907572637988050133502132224
            sage: ZZ(-2345.67)
            Traceback (most recent call last):
            ...
            TypeError: Attempt to coerce non-integral RealNumber to Integer
        """
        cdef Integer n

        if mpfr_integer_p(self.value):
            n = PY_NEW(Integer)
            mpfr_get_z(n.value, self.value, MPFR_RNDN)
            return n

        raise TypeError("Attempt to coerce non-integral RealNumber to Integer")

    def integer_part(self):
        """
        If in decimal this number is written ``n.defg``, returns ``n``.

        OUTPUT: a Sage Integer

        EXAMPLES::

            sage: a = 119.41212
            sage: a.integer_part()
            119
            sage: a = -123.4567
            sage: a.integer_part()
            -123

        A big number with no decimal point::

            sage: a = RR(10^17); a
            1.00000000000000e17
            sage: a.integer_part()
            100000000000000000
        """
        if not mpfr_number_p(self.value):
            raise ValueError('cannot convert infinity or NaN to Sage Integer')

        cdef Integer z = Integer()
        mpfr_get_z(z.value, self.value, MPFR_RNDZ)
        return z

    def fp_rank(self):
        r"""
        Return the floating-point rank of this number. That is, if you
        list the floating-point numbers of this precision in order, and
        number them starting with `0.0 \rightarrow 0` and extending
        the list to positive and negative infinity, returns the number
        corresponding to this floating-point number.

        EXAMPLES::

            sage: RR(0).fp_rank()
            0
            sage: RR(0).nextabove().fp_rank()
            1
            sage: RR(0).nextbelow().nextbelow().fp_rank()
            -2
            sage: RR(1).fp_rank()
            4835703278458516698824705            # 32-bit
            20769187434139310514121985316880385  # 64-bit
            sage: RR(-1).fp_rank()
            -4835703278458516698824705            # 32-bit
            -20769187434139310514121985316880385  # 64-bit
            sage: RR(1).fp_rank() - RR(1).nextbelow().fp_rank()
            1
            sage: RR(-infinity).fp_rank()
            -9671406552413433770278913            # 32-bit
            -41538374868278621023740371006390273  # 64-bit
            sage: RR(-infinity).fp_rank() - RR(-infinity).nextabove().fp_rank()
            -1
        """
        if mpfr_nan_p(self.value):
            raise ValueError("Cannot compute fp_rank of NaN")

        cdef Integer z = PY_NEW(Integer)

        cdef mp_exp_t EXP_MIN = mpfr_get_exp_min()
        cdef mp_exp_t EXP_MAX = mpfr_get_exp_max()
        # fp_rank(0.0) = 0
        # fp_rank(m*2^e-p) = (m-2^{p-1})+(e-EXP_MIN)*2^{p-1}+1
        #                  = m+(e-EXP_MIN-1)*2^{p-1}+1
        # fp_rank(infinity) = (EXP_MAX+1-EXP_MIN)*2^{p-1}+1
        # fp_rank(-x) = -fp_rank(x)

        cdef int sgn = mpfr_sgn(self.value)

        if sgn == 0:
            return z

        cdef mpfr_prec_t prec = (<RealField_class>self._parent)._prec

        if mpfr_inf_p(self.value):
            mpz_set_ui(z.value, EXP_MAX+1-EXP_MIN)
            mpz_mul_2exp(z.value, z.value, prec-1)
            mpz_add_ui(z.value, z.value, 1)
            if sgn < 0:
                mpz_neg(z.value, z.value)
            return z

        cdef mpz_t mantissa
        mpz_init(mantissa)
        cdef mp_exp_t exponent = mpfr_get_z_exp(mantissa, self.value)
        mpz_set_si(z.value, exponent+prec-EXP_MIN-1)
        mpz_mul_2exp(z.value, z.value, prec-1)
        mpz_add_ui(z.value, z.value, 1)
        if sgn > 0:
            mpz_add(z.value, z.value, mantissa)
        else:
            mpz_sub(z.value, z.value, mantissa)
            mpz_neg(z.value, z.value)
        mpz_clear(mantissa)
        return z

    def fp_rank_delta(self, RealNumber other):
        r"""
        Return the floating-point rank delta between ``self``
        and ``other``. That is, if the return value is
        positive, this is the number of times you have to call
        ``.nextabove()`` to get from ``self`` to ``other``.

        EXAMPLES::

            sage: [x.fp_rank_delta(x.nextabove()) for x in                              # needs sage.symbolic
            ....:    (RR(-infinity), -1.0, 0.0, 1.0, RR(pi), RR(infinity))]
            [1, 1, 1, 1, 1, 0]

        In the 2-bit floating-point field, one subsegment of the
        floating-point numbers is: 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32

        ::

            sage: R2 = RealField(2)
            sage: R2(1).fp_rank_delta(R2(2))
            2
            sage: R2(2).fp_rank_delta(R2(1))
            -2
            sage: R2(1).fp_rank_delta(R2(1048576))
            40
            sage: R2(24).fp_rank_delta(R2(4))
            -5
            sage: R2(-4).fp_rank_delta(R2(-24))
            -5

        There are lots of floating-point numbers around 0::

            sage: R2(-1).fp_rank_delta(R2(1))
            4294967298            # 32-bit
            18446744073709551618  # 64-bit
        """
        # We create the API for forward compatibility, because it can have
        # a (somewhat) more efficient implementation than this; but for now,
        # we just go with the stupid implementation.
        return other.fp_rank() - self.fp_rank()

    ########################
    #   Basic Arithmetic
    ########################

    def __add__(left, right):
        r"""
        TESTS::

            sage: RR(1) + RIF(1)
            doctest:...:
            DeprecationWarning: automatic conversions from floating-point numbers to intervals are deprecated
            See https://github.com/sagemath/sage/issues/15114 for details.
            2
            sage: import warnings; warnings.resetwarnings()
        """
        if have_same_parent(left, right):
            return (<RealNumber> left)._add_(right)
        try:
            from sage.rings.real_mpfi import RealIntervalFieldElement
        except ImportError:
            RealIntervalFieldElement = None
        if type(right) is RealIntervalFieldElement:
            return right.__add__(left)
        elif isinstance(left, RealNumber):
            return Element.__add__(left, right)
        else:
            return Element.__add__(right, left)

    def __sub__(left, right):
        r"""
        TESTS::

            sage: RR(2) - RIF(1)
            doctest:...:
            DeprecationWarning: automatic conversions from floating-point numbers to intervals are deprecated
            See https://github.com/sagemath/sage/issues/15114 for details.
            1
            sage: import warnings; warnings.resetwarnings()
        """
        if have_same_parent(left, right):
            return (<RealNumber> left)._sub_(right)
        try:
            from sage.rings.real_mpfi import RealIntervalFieldElement
        except ImportError:
            RealIntervalFieldElement = None
        if type(right) is RealIntervalFieldElement:
            return (-right).__add__(left)
        elif isinstance(left, RealNumber):
            return Element.__sub__(left, right)
        else:
            return Element.__rsub__(right, left)

    def __mul__(left, right):
        r"""
        TESTS::

            sage: RR(1) * RIF(1)
            doctest:...:
            DeprecationWarning: automatic conversions from floating-point numbers to intervals are deprecated
            See https://github.com/sagemath/sage/issues/15114 for details.
            1
            sage: import warnings; warnings.resetwarnings()
        """
        if have_same_parent(left, right):
            return (<RealNumber> left)._mul_(right)
        try:
            from sage.rings.real_mpfi import RealIntervalFieldElement
        except ImportError:
            RealIntervalFieldElement = None
        if type(right) is RealIntervalFieldElement:
            return right.__mul__(left)
        elif isinstance(left, RealNumber):
            return Element.__mul__(left, right)
        else:
            return Element.__rmul__(right, left)

    def __truediv__(left, right):
        r"""
        TESTS::

            sage: RR(1) / RIF(1/2)
            doctest:...:
            DeprecationWarning: automatic conversions from floating-point numbers to intervals are deprecated
            See https://github.com/sagemath/sage/issues/15114 for details.
            2
            sage: import warnings; warnings.resetwarnings()
        """
        if have_same_parent(left, right):
            return (<RealNumber> left)._div_(right)
        try:
            from sage.rings.real_mpfi import RealIntervalFieldElement
        except ImportError:
            RealIntervalFieldElement = None
        if type(right) is RealIntervalFieldElement:
            return right.__rtruediv__(left)
        elif isinstance(left, RealNumber):
            return Element.__truediv__(left, right)
        else:
            return Element.__rtruediv__(right, left)

    cpdef _add_(self, other):
        """
        Add two real numbers with the same parent.

        EXAMPLES::

            sage: R = RealField()
            sage: R(-1.5) + R(2.5) # indirect doctest
            1.00000000000000
        """
        cdef RealNumber x = self._new()
        mpfr_add(x.value, self.value, (<RealNumber>other).value, (<RealField_class>self._parent).rnd)
        return x

    def __invert__(self):
        """
        Return the reciprocal of ``self``.

        EXAMPLES::

            sage: ~RR(5/2)
            0.400000000000000
            sage: ~RR(1)
            1.00000000000000
            sage: ~RR(0)
            +infinity
        """
        return self._parent(1) / self

    cpdef _sub_(self, right):
        """
        Subtract two real numbers with the same parent.

        EXAMPLES::

            sage: R = RealField()
            sage: R(-1.5) - R(2.5) # indirect doctest
            -4.00000000000000
        """
        cdef RealNumber x = self._new()
        mpfr_sub(x.value, self.value, (<RealNumber>right).value, (<RealField_class> self._parent).rnd)
        return x

    def _sympy_(self):
        """
        Return a sympy object of ``self``.

        EXAMPLES::

            sage: RealField(100)(1/7)._sympy_()                                         # needs sympy
            0.14285714285714285714285714286
            sage: type(_)                                                               # needs sympy
            <class 'sympy.core.numbers.Float'>

        TESTS:

        An indirect doctest to check this (see :issue:`14915`)::

            sage: x,y = var('x, y')                                                     # needs sage.symbolic
            sage: integrate(y, y, 0.5, 8*log(x), algorithm='sympy')                     # needs sympy sage.symbolic
            32*log(x)^2 - 0.125000000000000

        Check that :issue:`28903` is fixed::

            sage: (10.0^400)._sympy_()                                                  # needs sympy
            1.00000000000000e+400
        """
        import sympy
        return sympy.Float(self, precision=self._parent.precision())

    cpdef _mul_(self, right):
        """
        Multiply two real numbers with the same parent.

        EXAMPLES::

            sage: R = RealField()
            sage: R(-1.5) * R(2.5) # indirect doctest
            -3.75000000000000

        If two elements have different precision, arithmetic operations are
        performed after coercing to the lower precision::

            sage: R20 = RealField(20)
            sage: R100 = RealField(100)
            sage: a = R20('393.3902834028345')
            sage: b = R100('393.3902834028345')
            sage: a
            393.39
            sage: b
            393.39028340283450000000000000
            sage: a*b
            154760.
            sage: b*a
            154760.
            sage: parent(b*a)
            Real Field with 20 bits of precision
        """
        cdef RealNumber x = self._new()
        mpfr_mul(x.value, self.value, (<RealNumber>right).value, (<RealField_class>self._parent).rnd)
        return x

    cpdef _div_(self, right):
        """
        Divide ``self`` by ``other``, where both are real numbers with the same
        parent.

        EXAMPLES::

            sage: RR(1)/RR(3) # indirect doctest
            0.333333333333333
            sage: RR(1)/RR(0)
            +infinity
            sage: R = RealField()
            sage: R(-1.5) / R(2.5)
            -0.600000000000000
        """
        cdef RealNumber x = self._new()
        mpfr_div((<RealNumber>x).value, self.value,
                 (<RealNumber>right).value, (<RealField_class>self._parent).rnd)
        return x

    cpdef _neg_(self):
        """
        Return the negative of ``self``.

        EXAMPLES::

            sage: RR(1)._neg_()
            -1.00000000000000
            sage: RR('inf')._neg_()
            -infinity
            sage: RR('nan')._neg_()
            NaN
        """
        cdef RealNumber x = self._new()
        mpfr_neg(x.value, self.value, (<RealField_class>self._parent).rnd)
        return x

    def __abs__(self):
        """
        Return the absolute value of ``self``.

        EXAMPLES::

            sage: RR(-1).__abs__()
            1.00000000000000
            sage: RR('-inf').__abs__()
            +infinity
            sage: RR('inf').__abs__()
            +infinity
            sage: RR('nan').__abs__()
            NaN
        """
        return self.abs()

    cdef RealNumber abs(RealNumber self):
        """
        Return the absolute value of ``self``.

        EXAMPLES::

            sage: RR('-1').abs()
            1.00000000000000
            sage: RR('-inf').abs()
            +infinity
            sage: RR('inf').abs()
            +infinity
            sage: RR('nan').abs()
            NaN
        """
        cdef RealNumber x = self._new()
        mpfr_abs(x.value, self.value, (<RealField_class>self._parent).rnd)
        return x

    # Bit shifting
    def _lshift_(RealNumber self, n):
        """
        Return ``self * (2^n)`` for an integer ``n``.

        EXAMPLES::

            sage: RR(1.0)._lshift_(32)
            4.29496729600000e9
            sage: RR(1.5)._lshift_(2)
            6.00000000000000
        """
        cdef RealNumber x
        if n > sys.maxsize:
            raise OverflowError("n (=%s) must be <= %s" % (n, sys.maxsize))
        x = self._new()
        mpfr_mul_2ui(x.value, self.value, n, (<RealField_class>self._parent).rnd)
        return x

    def __lshift__(x, y):
        """
        Return ``self * (2^n)`` for an integer ``n``.

        EXAMPLES::

            sage: 1.0 << 32
            4.29496729600000e9
            sage: 1.5 << 2.5
            Traceback (most recent call last):
            ...
            TypeError: unsupported operands for <<
        """
        if not isinstance(x, RealNumber):
            raise TypeError("unsupported operands for <<")
        try:
            return x._lshift_(Integer(y))
        except TypeError:
            raise TypeError("unsupported operands for <<")

    def _rshift_(RealNumber self, n):
        """
        Return ``self / (2^n)`` for an integer ``n``.

        EXAMPLES::

            sage: RR(1.0)._rshift_(32)
            2.32830643653870e-10
            sage: RR(1.5)._rshift_(2)
            0.375000000000000
        """
        if n > sys.maxsize:
            raise OverflowError("n (=%s) must be <= %s" % (n, sys.maxsize))
        cdef RealNumber x = self._new()
        mpfr_div_2exp(x.value, self.value, n, (<RealField_class>self._parent).rnd)
        return x

    def __rshift__(x, y):
        """
        Return ``self / (2^n)`` for an integer ``n``.

        EXAMPLES::

            sage: 1024.0 >> 7
            8.00000000000000
            sage: 1.5 >> 2.5
            Traceback (most recent call last):
            ...
            TypeError: unsupported operands for >>
        """
        if not isinstance(x, RealNumber):
            raise TypeError("unsupported operands for >>")
        try:
            return x._rshift_(Integer(y))
        except TypeError:
            raise TypeError("unsupported operands for >>")

    def multiplicative_order(self):
        """
        Return the multiplicative order of ``self``.

        EXAMPLES::

            sage: RR(1).multiplicative_order()
            1
            sage: RR(-1).multiplicative_order()
            2
            sage: RR(3).multiplicative_order()
            +Infinity
        """

        if self == 1:
            return 1
        elif self == -1:
            return 2
        return sage.rings.infinity.infinity

    def sign(self):
        """
        Return ``+1`` if ``self`` is positive, ``-1`` if ``self`` is negative,
        and ``0`` if ``self`` is zero.

        EXAMPLES::

            sage: R=RealField(100)
            sage: R(-2.4).sign()
            -1
            sage: R(2.1).sign()
            1
            sage: R(0).sign()
            0
        """
        return mpfr_sgn(self.value)

    def precision(self):
        """
        Return the precision of ``self``.

        EXAMPLES::

            sage: RR(1.0).precision()
            53
            sage: RealField(101)(-1).precision()
            101
        """
        return Integer((<RealField_class>self._parent)._prec)

    prec = precision  # alias

    def conjugate(self):
        """
        Return the complex conjugate of this real number, which is the
        number itself.

        EXAMPLES::

            sage: x = RealField(100)(1.238)
            sage: x.conjugate()
            1.2380000000000000000000000000
        """
        return self

    def ulp(self, field=None):
        """
        Return the unit of least precision of ``self``, which is the
        weight of the least significant bit of ``self``. This is always
        a strictly positive number. It is also the gap between this
        number and the closest number with larger absolute value that
        can be represented.

        INPUT:

        - ``field`` -- :class:`RealField` used as parent of the result;
          if not specified, use ``parent(self)``

        .. NOTE::

            The ulp of zero is defined as the smallest representable
            positive number. For extremely small numbers, underflow
            occurs and the output is also the smallest representable
            positive number (the rounding mode is ignored, this
            computation is done by rounding towards +infinity).

        .. SEEALSO::

            :meth:`epsilon` for a scale-invariant version of this.

        EXAMPLES::

            sage: a = 1.0
            sage: a.ulp()
            2.22044604925031e-16
            sage: (-1.5).ulp()
            2.22044604925031e-16
            sage: a + a.ulp() == a
            False
            sage: a + a.ulp()/2 == a
            True

            sage: a = RealField(500).pi()
            sage: b = a + a.ulp()
            sage: (a+b)/2 in [a,b]
            True

        The ulp of zero is the smallest nonzero number::

            sage: a = RR(0).ulp()
            sage: a
            2.38256490488795e-323228497            # 32-bit
            8.50969131174084e-1388255822130839284  # 64-bit
            sage: a.fp_rank()
            1

        The ulp of very small numbers results in underflow, so the
        smallest nonzero number is returned instead::

            sage: a.ulp() == a
            True

        We use a different field::

            sage: a = RealField(256).pi()
            sage: a.ulp()
            3.454467422037777850154540745120159828446400145774512554009481388067436721265e-77
            sage: e = a.ulp(RealField(64))
            sage: e
            3.45446742203777785e-77
            sage: parent(e)
            Real Field with 64 bits of precision
            sage: e = a.ulp(QQ)
            Traceback (most recent call last):
            ...
            TypeError: field argument must be a RealField

        For infinity and NaN, we get back positive infinity and NaN::

            sage: a = RR(infinity)
            sage: a.ulp()
            +infinity
            sage: (-a).ulp()
            +infinity
            sage: a = RR('nan')
            sage: a.ulp()
            NaN
            sage: parent(RR('nan').ulp(RealField(42)))
            Real Field with 42 bits of precision
        """
        cdef RealField_class _parent
        if field is None:
            _parent = self._parent
        else:
            try:
                _parent = field
            except TypeError:
                raise TypeError("field argument must be a RealField")

        cdef RealNumber x = _parent._new()
        cdef mp_exp_t e
        if mpfr_regular_p(self.value):
            # Non-zero number
            e = mpfr_get_exp(self.value) - mpfr_get_prec(self.value)
            # Round up in case of underflow
            mpfr_set_ui_2exp(x.value, 1, e, MPFR_RNDU)
        else:
            # Special cases: zero, infinity, NaN
            if mpfr_zero_p(self.value):
                mpfr_set_zero(x.value, 1)
                mpfr_nextabove(x.value)
            elif mpfr_inf_p(self.value):
                mpfr_set_inf(x.value, 1)
            else:  # NaN
                mpfr_set_nan(x.value)
        return x

    def epsilon(self, field=None):
        """
        Return ``abs(self)`` divided by `2^b` where `b` is the
        precision in bits of ``self``. Equivalently, return
        ``abs(self)`` multiplied by the :meth:`ulp` of 1.

        This is a scale-invariant version of :meth:`ulp` and it lies
        in `[u/2, u)` where `u` is ``self.ulp()`` (except in the case
        of zero or underflow).

        INPUT:

        - ``field`` -- :class:`RealField` used as parent of the result
          If not specified, use ``parent(self)``

        OUTPUT:

        ``field(self.abs() / 2^self.precision())``

        EXAMPLES::

            sage: RR(2^53).epsilon()
            1.00000000000000
            sage: RR(0).epsilon()
            0.000000000000000
            sage: a = RR.pi()
            sage: a.epsilon()
            3.48786849800863e-16
            sage: a.ulp()/2, a.ulp()
            (2.22044604925031e-16, 4.44089209850063e-16)
            sage: a / 2^a.precision()
            3.48786849800863e-16
            sage: (-a).epsilon()
            3.48786849800863e-16

        We use a different field::

            sage: a = RealField(256).pi()
            sage: a.epsilon()
            2.713132368784788677624750042896586252980746500631892201656843478528498954308e-77
            sage: e = a.epsilon(RealField(64))
            sage: e
            2.71313236878478868e-77
            sage: parent(e)
            Real Field with 64 bits of precision
            sage: e = a.epsilon(QQ)
            Traceback (most recent call last):
            ...
            TypeError: field argument must be a RealField

        Special values::

            sage: RR('nan').epsilon()
            NaN
            sage: parent(RR('nan').epsilon(RealField(42)))
            Real Field with 42 bits of precision
            sage: RR('+Inf').epsilon()
            +infinity
            sage: RR('-Inf').epsilon()
            +infinity
        """
        cdef RealField_class _parent
        if field is None:
            _parent = self._parent
        else:
            try:
                _parent = field
            except TypeError:
                raise TypeError("field argument must be a RealField")

        cdef RealNumber x = _parent._new()
        mpfr_div_2exp(x.value, self.value, mpfr_get_prec(self.value), _parent.rnd)
        mpfr_abs(x.value, x.value, _parent.rnd)
        return x

    ###################
    # Rounding etc.
    ###################

    cpdef _mod_(left, right):
        """
        Return the value of ``left - n*right``, rounded according to the
        rounding mode of the parent, where ``n`` is the integer quotient of
        ``x`` divided by ``y``. The integer ``n`` is rounded toward the
        nearest integer (ties rounded to even).

        EXAMPLES::

            sage: 10.0 % 2r
            0.000000000000000
            sage: 20r % .5
            0.000000000000000

            sage 1.1 % 0.25
            0.100000000000000
        """
        cdef RealNumber x
        x = (<RealNumber>left)._new()
        mpfr_remainder (x.value, (<RealNumber>left).value,
                        (<RealNumber>right).value,
                        (<RealField_class>(<RealNumber>left)._parent).rnd)
        return x

    def round(self):
        """
         Round ``self`` to the nearest representable integer, rounding halfway
         cases away from zero.

         .. NOTE::

             The rounding mode of the parent field does not affect the result.

         EXAMPLES::

             sage: RR(0.49).round()
             0
             sage: RR(0.5).round()
             1
             sage: RR(-0.49).round()
             0
             sage: RR(-0.5).round()
             -1
         """
        cdef RealNumber x = self._new()
        mpfr_round(x.value, self.value)
        return x.integer_part()

    def floor(self):
        """
        Return the floor of ``self``.

        EXAMPLES::

            sage: R = RealField()
            sage: (2.99).floor()
            2
            sage: (2.00).floor()
            2
            sage: floor(RR(-5/2))
            -3
            sage: floor(RR(+infinity))
            Traceback (most recent call last):
            ...
            ValueError: Calling floor() on infinity or NaN
        """
        cdef RealNumber x
        if not mpfr_number_p(self.value):
            raise ValueError('Calling floor() on infinity or NaN')
        x = self._new()
        mpfr_floor(x.value, self.value)
        return x.integer_part()

    def ceil(self):
        """
        Return the ceiling of ``self``.

        EXAMPLES::

            sage: (2.99).ceil()
            3
            sage: (2.00).ceil()
            2
            sage: (2.01).ceil()
            3

        ::

            sage: ceil(10^16 * 1.0)
            10000000000000000
            sage: ceil(10^17 * 1.0)
            100000000000000000
            sage: ceil(RR(+infinity))
            Traceback (most recent call last):
            ...
            ValueError: Calling ceil() on infinity or NaN
        """
        cdef RealNumber x
        if not mpfr_number_p(self.value):
            raise ValueError('Calling ceil() on infinity or NaN')
        x = self._new()
        mpfr_ceil(x.value, self.value)
        return x.integer_part()

    ceiling = ceil  # alias

    def trunc(self):
        """
        Truncate ``self``.

        EXAMPLES::

            sage: (2.99).trunc()
            2
            sage: (-0.00).trunc()
            0
            sage: (0.00).trunc()
            0
        """
        cdef RealNumber x = self._new()
        mpfr_trunc(x.value, self.value)
        return x.integer_part()

    def frac(self):
        """
        Return a real number such that
        ``self = self.trunc() + self.frac()``.  The return value will also
        satisfy ``-1 < self.frac() < 1``.

        EXAMPLES::

            sage: (2.99).frac()
            0.990000000000000
            sage: (2.50).frac()
            0.500000000000000
            sage: (-2.79).frac()
            -0.790000000000000
            sage: (-2.79).trunc() + (-2.79).frac()
            -2.79000000000000
        """
        cdef RealNumber x
        x = self._new()
        mpfr_frac(x.value, self.value, (<RealField_class>self._parent).rnd)
        return x

    def nexttoward(self, other):
        """
        Return the floating-point number adjacent to ``self`` which is closer
        to ``other``. If ``self`` or other is ``NaN``, returns ``NaN``; if
        ``self`` equals ``other``, returns ``self``.

        EXAMPLES::

            sage: (1.0).nexttoward(2).str()
            '1.0000000000000002'
            sage: (1.0).nexttoward(RR('-infinity')).str()
            '0.99999999999999989'
            sage: RR(infinity).nexttoward(0)
            2.09857871646739e323228496            # 32-bit
            5.87565378911159e1388255822130839282  # 64-bit
            sage: RR(pi).str()                                                          # needs sage.symbolic
            '3.1415926535897931'
            sage: RR(pi).nexttoward(22/7).str()                                         # needs sage.symbolic
            '3.1415926535897936'
            sage: RR(pi).nexttoward(21/7).str()                                         # needs sage.symbolic
            '3.1415926535897927'
        """
        cdef RealNumber other_rn
        if isinstance(other, RealNumber):
            other_rn = other
        else:
            other_rn = self._parent(other)

        cdef RealNumber x = self._new()

        mpfr_set(x.value, self.value, MPFR_RNDN)
        mpfr_nexttoward(x.value, other_rn.value)

        return x

    def nextabove(self):
        """
        Return the next floating-point number larger than ``self``.

        EXAMPLES::

            sage: RR('-infinity').nextabove()
            -2.09857871646739e323228496            # 32-bit
            -5.87565378911159e1388255822130839282  # 64-bit
            sage: RR(0).nextabove()
            2.38256490488795e-323228497            # 32-bit
            8.50969131174084e-1388255822130839284  # 64-bit
            sage: RR('+infinity').nextabove()
            +infinity
            sage: RR(-sqrt(2)).str()                                                    # needs sage.symbolic
            '-1.4142135623730951'
            sage: RR(-sqrt(2)).nextabove().str()                                        # needs sage.symbolic
            '-1.4142135623730949'
        """

        cdef RealNumber x = self._new()
        mpfr_set(x.value, self.value, MPFR_RNDN)
        mpfr_nextabove(x.value)

        return x

    def nextbelow(self):
        """
        Return the next floating-point number smaller than ``self``.

        EXAMPLES::

            sage: RR('-infinity').nextbelow()
            -infinity
            sage: RR(0).nextbelow()
            -2.38256490488795e-323228497            # 32-bit
            -8.50969131174084e-1388255822130839284  # 64-bit
            sage: RR('+infinity').nextbelow()
            2.09857871646739e323228496              # 32-bit
            5.87565378911159e1388255822130839282    # 64-bit
            sage: RR(-sqrt(2)).str()                                                    # needs sage.symbolic
            '-1.4142135623730951'
            sage: RR(-sqrt(2)).nextbelow().str()                                        # needs sage.symbolic
            '-1.4142135623730954'
        """

        cdef RealNumber x = self._new()
        mpfr_set(x.value, self.value, MPFR_RNDN)
        mpfr_nextbelow(x.value)

        return x

    ###########################################
    # Conversions
    ###########################################

    def __float__(self):
        """
        Return a Python float approximating ``self``.

        EXAMPLES::

            sage: RR(pi).__float__()                                                    # needs sage.symbolic
            3.141592653589793
            sage: type(RR(pi).__float__())                                              # needs sage.symbolic
            <... 'float'>
        """
        return mpfr_get_d(self.value, (<RealField_class>self._parent).rnd)

    def _rpy_(self):
        """
        Return ``self.__float__()`` for rpy to convert into the
        appropriate R object.

        EXAMPLES::

            sage: n = RealNumber(2.0)
            sage: n._rpy_()
            2.0
            sage: type(n._rpy_())
            <... 'float'>
        """
        return self.__float__()

    def __int__(self):
        """
        Return the Python integer truncation of ``self``.

        EXAMPLES::

            sage: RR(pi).__int__()                                                      # needs sage.symbolic
            3
            sage: type(RR(pi).__int__())                                                # needs sage.symbolic
            <... 'int'>
        """
        if not mpfr_number_p(self.value):
            raise ValueError('cannot convert infinity or NaN to Python int')

        cdef Integer z = Integer()
        mpfr_get_z(z.value, self.value, MPFR_RNDZ)
        return z.__int__()

    def __complex__(self):
        """
        Return a Python complex number equal to ``self``.

        EXAMPLES::

            sage: RR(pi).__complex__()                                                  # needs sage.symbolic
            (3.141592653589793+0j)
            sage: type(RR(pi).__complex__())                                            # needs sage.symbolic
            <... 'complex'>
        """

        return complex(float(self))

    def _complex_number_(self):
        """
        Return a Sage complex number equal to ``self``.

        EXAMPLES::

            sage: RR(pi)._complex_number_()                                             # needs sage.symbolic
            3.14159265358979
            sage: parent(RR(pi)._complex_number_())                                     # needs sage.symbolic
            Complex Field with 53 bits of precision
        """
        from sage.rings.complex_mpfr import ComplexField
        return ComplexField(self.prec())(self)

    def _axiom_(self, axiom):
        """
        Return ``self`` as a floating point number in Axiom.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: R = RealField(100)
            sage: R(pi)
            3.1415926535897932384626433833
            sage: axiom(R(pi))  # indirect doctest      # optional - axiom
            3.1415926535 8979323846 26433833
            sage: fricas(R(pi))                         # optional - fricas
            3.1415926535_8979323846_26433833
        """
        prec = self.parent().prec()

        # Set the precision in Axiom
        old_prec = axiom('precision(%s)$Float' % prec)
        res = axiom('%s :: Float' % self.exact_rational())
        axiom.eval('precision(%s)$Float' % old_prec)

        return res

    _fricas_ = _axiom_

    def __pari__(self):
        """
        Return ``self`` as a Pari floating-point number.

        EXAMPLES::

            sage: RR(2.0).__pari__()                                                    # needs sage.libs.pari
            2.00000000000000

        The current Pari precision affects the printing of this number, but
        Pari does maintain the same 250-bit number on both 32-bit and
        64-bit platforms::

            sage: # needs sage.libs.pari
            sage: RealField(250).pi().__pari__()
            3.14159265358979
            sage: RR(0.0).__pari__()
            0.E-19
            sage: RR(-1.234567).__pari__()
            -1.23456700000000
            sage: RR(2.0).sqrt().__pari__()
            1.41421356237310
            sage: RR(2.0).sqrt().__pari__().sage()
            1.41421356237309515
            sage: RR(2.0).sqrt().__pari__().sage().prec()
            64
            sage: RealField(70)(pi).__pari__().sage().prec()                            # needs sage.symbolic
            96                                         # 32-bit
            128                                        # 64-bit
            sage: for i in range(100, 200):
            ....:     assert(RR(i).sqrt() == RR(i).sqrt().__pari__().sage())

        TESTS:

        Check that we create real zeros without mantissa::

            sage: RDF(0).__pari__().sizeword()                                          # needs sage.libs.pari
            2
            sage: RealField(100)(0.0).__pari__().sizeword()                             # needs sage.libs.pari
            2

        Check that the largest and smallest exponents representable by
        PARI convert correctly::

            sage: # needs sage.libs.pari
            sage: a = pari(0.5) << (sys.maxsize+1)/4
            sage: RR(a) >> (sys.maxsize+1)/4
            0.500000000000000
            sage: a = pari(0.5) >> (sys.maxsize-3)/4
            sage: RR(a) << (sys.maxsize-3)/4
            0.500000000000000
        """
        return new_gen_from_real_mpfr_element(self)

    def _mpmath_(self, prec=None, rounding=None):
        """
        Return an mpmath version of this :class:`RealNumber`.

        .. NOTE::

           Currently the rounding mode is ignored.

        EXAMPLES::

            sage: RR(-1.5)._mpmath_()
            mpf('-1.5')
        """
        if prec is not None:
            return RealField(prec)(self)._mpmath_()
        from sage.libs.mpmath.all import make_mpf
        return make_mpf(mpfr_to_mpfval(self.value))

    def sign_mantissa_exponent(self):
        r"""
        Return the sign, mantissa, and exponent of ``self``.

        In Sage (as in MPFR), floating-point numbers of precision `p`
        are of the form `s m 2^{e-p}`, where `s \in \{-1, 1\}`,
        `2^{p-1} \leq m < 2^p`, and `-2^{30} + 1 \leq e \leq 2^{30} -
        1`; plus the special values ``+0``, ``-0``, ``+infinity``,
        ``-infinity``, and ``NaN`` (which stands for Not-a-Number).

        This function returns `s`, `m`, and `e-p`.  For the special values:

        - ``+0`` returns ``(1, 0, 0)`` (analogous to IEEE-754;
          note that MPFR actually stores the exponent as "smallest exponent
          possible")
        - ``-0`` returns ``(-1, 0, 0)`` (analogous to IEEE-754;
          note that MPFR actually stores the exponent as "smallest exponent
          possible")
        - the return values for ``+infinity``, ``-infinity``, and ``NaN`` are
          not specified.

        EXAMPLES::

            sage: R = RealField(53)
            sage: a = R(exp(1.0)); a
            2.71828182845905
            sage: sign, mantissa, exponent = R(exp(1.0)).sign_mantissa_exponent()
            sage: sign, mantissa, exponent
            (1, 6121026514868073, -51)
            sage: sign*mantissa*(2**exponent) == a
            True

        The mantissa is always a nonnegative number (see :issue:`14448`)::

            sage: RR(-1).sign_mantissa_exponent()
            (-1, 4503599627370496, -52)

        We can also calculate this also using `p`-adic valuations::

            sage: a = R(exp(1.0))
            sage: b = a.exact_rational()
            sage: valuation, unit = b.val_unit(2)
            sage: (b/abs(b), unit, valuation)
            (1, 6121026514868073, -51)
            sage: a.sign_mantissa_exponent()
            (1, 6121026514868073, -51)

        TESTS::

            sage: R('+0').sign_mantissa_exponent()
            (1, 0, 0)
            sage: R('-0').sign_mantissa_exponent()
            (-1, 0, 0)
        """
        cdef Integer mantissa
        cdef Integer sign
        cdef mp_exp_t exponent

        if mpfr_signbit(self.value) == 0:
            sign = Integer(1)
        else:
            sign = Integer(-1)

        if mpfr_zero_p(self.value):
            mantissa = Integer(0)
            exponent = Integer(0)
        else:
            mantissa = Integer()
            exponent = mpfr_get_z_exp(mantissa.value, self.value)

        return sign, sign*mantissa, Integer(exponent)

    def exact_rational(self):
        """
        Return the exact rational representation of this floating-point
        number.

        EXAMPLES::

            sage: RR(0).exact_rational()
            0
            sage: RR(1/3).exact_rational()
            6004799503160661/18014398509481984
            sage: RR(37/16).exact_rational()
            37/16
            sage: RR(3^60).exact_rational()
            42391158275216203520420085760
            sage: RR(3^60).exact_rational() - 3^60
            6125652559
            sage: RealField(5)(-pi).exact_rational()                                    # needs sage.symbolic
            -25/8

        TESTS::

            sage: RR('nan').exact_rational()
            Traceback (most recent call last):
            ...
            ValueError: unable to convert NaN to a rational number
            sage: RR('-infinity').exact_rational()
            Traceback (most recent call last):
            ...
            ValueError: unable to convert -infinity to a rational number
        """
        if not mpfr_number_p(self.value):
            raise ValueError(f"unable to convert {self} to a rational number")

        if mpfr_sgn(self.value) == 0:
            return Rational(0)

        mantissa = <Integer>Integer.__new__(Integer)
        exponent = mpfr_get_z_exp(mantissa.value, self.value)

        return Rational(mantissa) << exponent

    def as_integer_ratio(self):
        """
        Return a coprime pair of integers ``(a, b)`` such that ``self``
        equals ``a / b`` exactly.

        EXAMPLES::

            sage: RR(0).as_integer_ratio()
            (0, 1)
            sage: RR(1/3).as_integer_ratio()
            (6004799503160661, 18014398509481984)
            sage: RR(37/16).as_integer_ratio()
            (37, 16)
            sage: RR(3^60).as_integer_ratio()
            (42391158275216203520420085760, 1)
            sage: RR('nan').as_integer_ratio()
            Traceback (most recent call last):
            ...
            ValueError: unable to convert NaN to a rational number

        This coincides with Python floats::

            sage: pi = RR.pi()
            sage: pi.as_integer_ratio()
            (884279719003555, 281474976710656)
            sage: float(pi).as_integer_ratio() == pi.as_integer_ratio()
            True
        """
        return self.exact_rational().as_integer_ratio()

    def simplest_rational(self):
        """
        Return the simplest rational which is equal to ``self`` (in the Sage
        sense). Recall that Sage defines the equality operator by coercing
        both sides to a single type and then comparing; thus, this finds
        the simplest rational which (when coerced to this RealField) is
        equal to ``self``.

        Given rationals `a / b` and `c / d` (both in lowest terms), the former
        is simpler if `b < d` or if `b = d` and `|a| < |c|`.

        The effect of rounding modes is slightly counter-intuitive.
        Consider the case of round-toward-minus-infinity. This rounding is
        performed when coercing a rational to a floating-point number; so
        the :meth:`simplest_rational()` of a round-to-minus-infinity number
        will be either exactly equal to or slightly larger than the number.

        EXAMPLES::

            sage: RRd = RealField(53, rnd='RNDD')
            sage: RRz = RealField(53, rnd='RNDZ')
            sage: RRu = RealField(53, rnd='RNDU')
            sage: RRa = RealField(53, rnd='RNDA')
            sage: def check(x):
            ....:     rx = x.simplest_rational()
            ....:     assert x == rx
            ....:     return rx
            sage: RRd(1/3) < RRu(1/3)
            True
            sage: check(RRd(1/3))
            1/3
            sage: check(RRu(1/3))
            1/3
            sage: check(RRz(1/3))
            1/3
            sage: check(RRa(1/3))
            1/3
            sage: check(RR(1/3))
            1/3
            sage: check(RRd(-1/3))
            -1/3
            sage: check(RRu(-1/3))
            -1/3
            sage: check(RRz(-1/3))
            -1/3
            sage: check(RRa(-1/3))
            -1/3
            sage: check(RR(-1/3))
            -1/3
            sage: check(RealField(20)(pi))                                              # needs sage.symbolic
            355/113
            sage: check(RR(pi))                                                         # needs sage.symbolic
            245850922/78256779
            sage: check(RR(2).sqrt())
            131836323/93222358
            sage: check(RR(1/2^210))
            1/1645504557321205859467264516194506011931735427766374553794641921
            sage: check(RR(2^210))
            1645504557321205950811116849375918117252433820865891134852825088
            sage: (RR(17).sqrt()).simplest_rational()^2 - 17
            -1/348729667233025
            sage: (RR(23).cube_root()).simplest_rational()^3 - 23
            -1404915133/264743395842039084891584
            sage: RRd5 = RealField(5, rnd='RNDD')
            sage: RRu5 = RealField(5, rnd='RNDU')
            sage: RR5 = RealField(5)
            sage: below1 = RR5(1).nextbelow()
            sage: check(RRd5(below1))
            31/32
            sage: check(RRu5(below1))
            16/17
            sage: check(below1)
            21/22
            sage: below1.exact_rational()
            31/32
            sage: above1 = RR5(1).nextabove()
            sage: check(RRd5(above1))
            10/9
            sage: check(RRu5(above1))
            17/16
            sage: check(above1)
            12/11
            sage: above1.exact_rational()
            17/16
            sage: check(RR(1234))
            1234
            sage: check(RR5(1234))
            1185
            sage: check(RR5(1184))
            1120
            sage: RRd2 = RealField(2, rnd='RNDD')
            sage: RRu2 = RealField(2, rnd='RNDU')
            sage: RR2 = RealField(2)
            sage: check(RR2(8))
            7
            sage: check(RRd2(8))
            8
            sage: check(RRu2(8))
            7
            sage: check(RR2(13))
            11
            sage: check(RRd2(13))
            12
            sage: check(RRu2(13))
            13
            sage: check(RR2(16))
            14
            sage: check(RRd2(16))
            16
            sage: check(RRu2(16))
            13
            sage: check(RR2(24))
            21
            sage: check(RRu2(24))
            17
            sage: check(RR2(-24))
            -21
            sage: check(RRu2(-24))
            -24

        TESTS::

            sage: RR('nan').simplest_rational()
            Traceback (most recent call last):
            ...
            ValueError: cannot convert NaN or infinity to rational number
            sage: RR('-infinity').simplest_rational()
            Traceback (most recent call last):
            ...
            ValueError: cannot convert NaN or infinity to rational number
        """
        if not mpfr_number_p(self.value):
            raise ValueError('cannot convert NaN or infinity to rational number')

        if mpfr_zero_p(self.value):
            return Rational(0)

        from sage.rings.real_mpfi import RealIntervalField

        cdef mpfr_rnd_t rnd = (<RealField_class>self._parent).rnd
        cdef mpfr_prec_t prec = (<RealField_class>self._parent)._prec

        cdef RealNumber low, high
        cdef int odd

        if rnd == MPFR_RNDN:
            # hp == "high precision"
            hp_field = RealField(prec + 1)
            hp_val = hp_field(self)
            hp_intv_field = RealIntervalField(prec + 1)
            low = hp_val.nextbelow()
            high = hp_val.nextabove()
            hp_intv = hp_intv_field(low, high)
            # In MPFR_RNDN mode, we round to nearest, preferring even mantissas
            # if we are exactly halfway between representable floats.
            # Thus, the values low and high will round to self iff the
            # mantissa of self is even.  (Note that this only matters
            # if low or high is an integer; if they are not integers,
            # then self is simpler than either low or high.)
            # Is there a better (faster) way to check this?
            odd = self._parent(low) != self
            return hp_intv.simplest_rational(low_open=odd, high_open=odd)

        if rnd == MPFR_RNDZ:
            if mpfr_sgn(self.value) > 0:
                rnd = MPFR_RNDD
            else:
                rnd = MPFR_RNDU
        elif rnd == MPFR_RNDA:
            if mpfr_sgn(self.value) > 0:
                rnd = MPFR_RNDU
            else:
                rnd = MPFR_RNDD

        intv_field = RealIntervalField(prec)

        if rnd == MPFR_RNDD:
            intv = intv_field(self, self.nextabove())
            return intv.simplest_rational(high_open = True)
        if rnd == MPFR_RNDU:
            intv = intv_field(self.nextbelow(), self)
            return intv.simplest_rational(low_open = True)

    def nearby_rational(self, max_error=None, max_denominator=None):
        """
        Find a rational near to ``self``. Exactly one of ``max_error`` or
        ``max_denominator`` must be specified.

        If ``max_error`` is specified, then this returns the simplest rational
        in the range ``[self-max_error .. self+max_error]``. If
        ``max_denominator`` is specified, then this returns the rational
        closest to ``self`` with denominator at most ``max_denominator``.
        (In case of ties, we pick the simpler rational.)

        EXAMPLES::

            sage: (0.333).nearby_rational(max_error=0.001)
            1/3
            sage: (0.333).nearby_rational(max_error=1)
            0
            sage: (-0.333).nearby_rational(max_error=0.0001)
            -257/772

        ::

            sage: (0.333).nearby_rational(max_denominator=100)
            1/3
            sage: RR(1/3 + 1/1000000).nearby_rational(max_denominator=2999999)
            777780/2333333
            sage: RR(1/3 + 1/1000000).nearby_rational(max_denominator=3000000)
            1000003/3000000
            sage: (-0.333).nearby_rational(max_denominator=1000)
            -333/1000
            sage: RR(3/4).nearby_rational(max_denominator=2)
            1

            sage: # needs sage.symbolic
            sage: RR(pi).nearby_rational(max_denominator=120)
            355/113
            sage: RR(pi).nearby_rational(max_denominator=10000)
            355/113
            sage: RR(pi).nearby_rational(max_denominator=100000)
            312689/99532
            sage: RR(pi).nearby_rational(max_denominator=1)
            3

            sage: RR(-3.5).nearby_rational(max_denominator=1)
            -3

        TESTS::

            sage: RR('nan').nearby_rational(max_denominator=1000)
            Traceback (most recent call last):
            ...
            ValueError: cannot convert NaN or infinity to rational number
            sage: RR('nan').nearby_rational(max_error=0.01)
            Traceback (most recent call last):
            ...
            ValueError: cannot convert NaN or infinity to rational number
            sage: RR(oo).nearby_rational(max_denominator=1000)
            Traceback (most recent call last):
            ...
            ValueError: cannot convert NaN or infinity to rational number
            sage: RR(oo).nearby_rational(max_error=0.01)
            Traceback (most recent call last):
            ...
            ValueError: cannot convert NaN or infinity to rational number
        """

        if not mpfr_number_p(self.value):
            raise ValueError('cannot convert NaN or infinity to rational number')

        if ((max_error is None and max_denominator is None) or
                (max_error is not None and max_denominator is not None)):
            raise ValueError('Must specify exactly one of max_error or max_denominator in nearby_rational()')

        if max_error is not None:
            from sage.rings.real_mpfi import RealIntervalField

            intv_field = RealIntervalField(self.prec())
            intv = intv_field(self - max_error, self + max_error)

            return intv.simplest_rational()

        cdef int sgn = mpfr_sgn(self.value)

        if sgn == 0:
            return Rational(0)

        cdef Rational self_r = self.exact_rational()

        cdef Integer self_d = self_r.denominator()

        if self_d <= max_denominator:
            return self_r

        if sgn < 0:
            self_r = -self_r

        cdef Integer fl = self_r.floor()
        cdef Rational target = self_r - fl

        cdef int low_done = 0
        cdef int high_done = 0

        # We use the Stern-Brocot tree to find the nearest neighbors of
        # self with denominator at most max_denominator.  However,
        # navigating the Stern-Brocot tree in the straightforward way
        # can be very slow; for instance, to get to 1/1000000 takes a
        # million steps.  Instead, we perform many steps at once;
        # this probably slows down the average case, but it drastically
        # speeds up the worst case.

        # Suppose we have a/b < c/d < e/f, where a/b and e/f are
        # neighbors in the Stern-Brocot tree and c/d is the target.
        # We alternate between moving the low and the high end toward
        # the target as much as possible.  Suppose that there are
        # k consecutive rightward moves in the Stern-Brocot tree
        # traversal; then we end up with (a+k*e)/(b+k*f).  We have
        # two constraints on k.  First, the result must be <= c/d;
        # this gives us the following:
        # (a+k*e)/(b+k*f) <= c/d
        # d*a + k*(d*e) <= c*b + k*(c*f)
        # k*(d*e) - k*(c*f) <= c*b - d*a
        # k <= (c*b - d*a)/(d*e - c*f)
        # when moving the high side, we get
        # (k*a+e)/(k*b+f) >= c/d
        # k*(d*a) + d*e >= k*(c*b) + c*f
        # d*e - c*f >= k*(c*b - d*a)
        # k <= (d*e - c*f)/(c*b - d*a)

        # We also need the denominator to be <= max_denominator; this
        # gives (b+k*f) <= max_denominator or
        # k <= (max_denominator - b)/f
        # or
        # k <= (max_denominator - f)/b

        # We use variables with the same names as in the math above.

        cdef Integer a = Integer(0)
        cdef Integer b = Integer(1)
        cdef Integer c = target.numerator()
        cdef Integer d = target.denominator()
        cdef Integer e = Integer(1)
        cdef Integer f = Integer(1)

        cdef Integer k

        while not low_done or not high_done:
            # Move the low side
            k = (c*b - d*a) // (d*e - c*f)

            if b+k*f > max_denominator:
                k = (max_denominator - b) // f
                low_done = True

            if k == 0:
                low_done = True

            a = a + k*e
            b = b + k*f

            # Move the high side
            k = (d*e - c*f) // (c*b - d*a)

            if k*b + f >= max_denominator:
                k = (max_denominator - f) // b
                high_done = True

            if k == 0:
                high_done = True

            e = k*a + e
            f = k*b + f

        # Now a/b and e/f are rationals surrounding c/d.  We know that
        # neither is equal to c/d, since d > max_denominator and
        # b and f are both <= max_denominator.  (We know that
        # d > max_denominator because we return early (before we
        # get here) if d <= max_denominator.)

        low = a / b
        high = e / f

        D0 = target - low
        D1 = high - target
        if D1 < D0:
            result = high
        elif D0 < D1:
            result = low
        else:
            if f < b:
                result = high
            elif b < f:
                result = low
            else:
                if e < a:
                    result = high
                else:
                    result = low

        result += fl

        if sgn < 0:
            return -result
        return result

    def __mpfr__(self):
        """
        Convert Sage ``RealNumber`` to gmpy2 ``mpfr``.

        EXAMPLES::

            sage: r = RR(4.12)
            sage: r.__mpfr__()
            mpfr('4.1200000000000001')
            sage: from gmpy2 import mpfr
            sage: mpfr(RR(4.5))
            mpfr('4.5')
            sage: R = RealField(127)
            sage: mpfr(R.pi()).precision
            127
            sage: R = RealField(42)
            sage: mpfr(R.pi()).precision
            42
            sage: R = RealField(256)
            sage: x = mpfr(R.pi())
            sage: x.precision
            256
            sage: y = R(x)
            sage: mpfr(y) == x
            True
            sage: x = mpfr('2.567e42', precision=128)
            sage: y = RealField(128)(x)
            sage: mpfr(y) == x
            True
        """
        return gmpy2.GMPy_MPFR_From_mpfr(self.value)

    ###########################################
    # Comparisons: ==, !=, <, <=, >, >=
    ###########################################

    def is_NaN(self):
        """
        Return ``True`` if ``self`` is Not-a-Number ``NaN``.

        EXAMPLES::

            sage: a = RR(0) / RR(0); a
            NaN
            sage: a.is_NaN()
            True
        """
        return mpfr_nan_p(self.value)

    def is_positive_infinity(self):
        r"""
        Return ``True`` if ``self`` is `+\infty`.

        EXAMPLES::

            sage: a = RR('1.494') / RR(0); a
            +infinity
            sage: a.is_positive_infinity()
            True
            sage: a = -RR('1.494') / RR(0); a
            -infinity
            sage: RR(1.5).is_positive_infinity()
            False
            sage: a.is_positive_infinity()
            False
        """
        return mpfr_inf_p(self.value) and mpfr_sgn(self.value) > 0

    def is_negative_infinity(self):
        r"""
        Return ``True`` if ``self`` is `-\infty`.

        EXAMPLES::

            sage: a = RR('1.494') / RR(0); a
            +infinity
            sage: a.is_negative_infinity()
            False
            sage: a = -RR('1.494') / RR(0); a
            -infinity
            sage: RR(1.5).is_negative_infinity()
            False
            sage: a.is_negative_infinity()
            True
        """
        return mpfr_inf_p(self.value) and mpfr_sgn(self.value) < 0

    def is_infinity(self):
        r"""
        Return ``True`` if ``self`` is `\infty` and ``False`` otherwise.

        EXAMPLES::

            sage: a = RR('1.494') / RR(0); a
            +infinity
            sage: a.is_infinity()
            True
            sage: a = -RR('1.494') / RR(0); a
            -infinity
            sage: a.is_infinity()
            True
            sage: RR(1.5).is_infinity()
            False
            sage: RR('nan').is_infinity()
            False
        """
        return mpfr_inf_p(self.value)

    def is_unit(self):
        """
        Return ``True`` if ``self`` is a unit (has a multiplicative inverse)
        and ``False`` otherwise.

        EXAMPLES::

            sage: RR(1).is_unit()
            True
            sage: RR('0').is_unit()
            False
            sage: RR('-0').is_unit()
            False
            sage: RR('nan').is_unit()
            False
            sage: RR('inf').is_unit()
            False
            sage: RR('-inf').is_unit()
            False
        """
        return mpfr_sgn(self.value) != 0 and not mpfr_inf_p(self.value)

    def is_real(self):
        """
        Return ``True`` if ``self`` is real (of course, this always returns
        ``True`` for a finite element of a real field).

        EXAMPLES::

            sage: RR(1).is_real()
            True
            sage: RR('-100').is_real()
            True
            sage: RR(NaN).is_real()                                                     # needs sage.symbolic
            False
        """
        return not mpfr_nan_p(self.value)

    def is_integer(self):
        """
        Return ``True`` if this number is a integer.

        EXAMPLES::

            sage: RR(1).is_integer()
            True
            sage: RR(0.1).is_integer()
            False
        """
        return mpfr_integer_p(self.value) != 0

    def __bool__(self):
        """
        Return ``True`` if ``self`` is nonzero.

        EXAMPLES::

            sage: bool(RR(1))
            True
            sage: bool(RR(0))
            False
            sage: bool(RR('inf'))
            True

        TESTS:

        Check that :issue:`20502` is fixed::

            sage: bool(RR('nan'))
            True
            sage: RR('nan').is_zero()
            False
        """
        return not mpfr_zero_p(self.value)

    cpdef _richcmp_(self, other, int op):
        """
        Compare ``self`` and ``other`` according to the rich
        comparison operator ``op``.

        EXAMPLES::

            sage: RR('-inf') < RR('inf')
            True
            sage: RR('-inf') < RR(-10000)
            True
            sage: RR(100000000) < RR('inf')
            True
            sage: RR(100000000) > RR('inf')
            False
            sage: RR(100000000) < RR('inf')
            True
            sage: RR(-1000) < RR(1000)
            True
            sage: RR(1) < RR(1).nextabove()
            True
            sage: RR(1) <= RR(1).nextabove()
            True
            sage: RR(1) <= RR(1)
            True
            sage: RR(1) < RR(1)
            False
            sage: RR(1) > RR(1)
            False
            sage: RR(1) >= RR(1)
            True
            sage: RR('inf') == RR('inf')
            True
            sage: RR('inf') == RR('-inf')
            False

        A ``NaN`` is not equal to anything, including itself::

            sage: RR('nan') == RR('nan')
            False
            sage: RR('nan') != RR('nan')
            True
            sage: RR('nan') < RR('nan')
            False
            sage: RR('nan') > RR('nan')
            False
            sage: RR('nan') <= RR('nan')
            False
            sage: RR('nan') >= RR('nan')
            False
            sage: RR('nan') == RR(0)
            False
            sage: RR('nan') != RR(0)
            True
            sage: RR('nan') < RR(0)
            False
            sage: RR('nan') > RR(0)
            False
            sage: RR('nan') <= RR(0)
            False
            sage: RR('nan') >= RR(0)
            False
        """
        cdef RealNumber y = <RealNumber>other
        if mpfr_nan_p(self.value) or mpfr_nan_p(y.value):
            return op == Py_NE
        cdef int c = mpfr_cmp(self.value, y.value)
        return rich_to_bool_sgn(op, c)

    ############################
    # Special Functions
    ############################

    def sqrt(self, extend=True, all=False):
        r"""
        The square root function.

        INPUT:

        - ``extend`` -- boolean (default: ``True``); if ``True``, return a
          square root in a complex field if necessary if ``self`` is negative.
          Otherwise raise a :exc:`ValueError`.

        - ``all`` -- boolean (default: ``False``); if ``True``, return a
          list of all square roots

        EXAMPLES::

            sage: r = -2.0
            sage: r.sqrt()
            1.41421356237310*I

        ::

            sage: r = 4.0
            sage: r.sqrt()
            2.00000000000000
            sage: r.sqrt()^2 == r
            True

        ::

            sage: r = 4344
            sage: r.sqrt()                                                              # needs sage.symbolic
            2*sqrt(1086)

        ::

            sage: r = 4344.0
            sage: r.sqrt()^2 == r
            True
            sage: r.sqrt()^2 - r
            0.000000000000000

        ::

            sage: r = -2.0
            sage: r.sqrt()
            1.41421356237310*I
        """
        cdef RealNumber x
        if mpfr_cmp_ui(self.value, 0) >= 0:
            x = self._new()
            if (<RealField_class>self._parent)._prec > 10*SIG_PREC_THRESHOLD:
                sig_on()
            mpfr_sqrt(x.value, self.value, (<RealField_class>self._parent).rnd)
            if (<RealField_class>self._parent)._prec > 10*SIG_PREC_THRESHOLD:
                sig_off()
            if all:
                if x.is_zero():
                    return [x]
                else:
                    return [x, -x]
            return x
        if not extend:
            raise ValueError("negative number %s does not have a square root in the real field" % self)
        return self._complex_number_().sqrt(all=all)

    def is_square(self):
        """
        Return whether or not this number is a square in this field. For
        the real numbers, this is ``True`` if and only if ``self`` is
        nonnegative.

        EXAMPLES::

            sage: r = 3.5
            sage: r.is_square()
            True
            sage: r = 0.0
            sage: r.is_square()
            True
            sage: r = -4.0
            sage: r.is_square()
            False
        """
        return mpfr_sgn(self.value) >= 0

    def cube_root(self):
        """
        Return the cubic root (defined over the real numbers) of ``self``.

        EXAMPLES::

            sage: r = 125.0; r.cube_root()
            5.00000000000000
            sage: r = -119.0
            sage: r.cube_root()^3 - r       # illustrates precision loss
            -1.42108547152020e-14
        """
        cdef RealNumber x = self._new()
        if (<RealField_class>self._parent)._prec > 10*SIG_PREC_THRESHOLD:
            sig_on()
        mpfr_cbrt(x.value, self.value, (<RealField_class>self._parent).rnd)
        if (<RealField_class>self._parent)._prec > 10*SIG_PREC_THRESHOLD:
            sig_off()
        return x

    def _pow(self, RealNumber exponent):
        """
        Compute ``self`` raised to the power of exponent, rounded in the
        direction specified by the parent of ``self``.

        EXAMPLES::

            sage: R = RealField(30)
            sage: a = R('1.23456')
            sage: a._pow(20.0)
            67.646297
        """
        cdef RealNumber x
        x = self._new()
        sig_on()
        mpfr_pow(x.value, self.value, exponent.value, (<RealField_class>self._parent).rnd)
        sig_off()
        if mpfr_nan_p(x.value):
            return self._complex_number_()**exponent._complex_number_()
        return x

    def __pow__(self, exponent, modulus):
        """
        Compute ``self`` raised to the power of exponent, rounded in the
        direction specified by the parent of ``self``.

        If the result is not a real number, ``self`` and the exponent are both
        coerced to complex numbers (with sufficient precision), then the
        exponentiation is computed in the complex numbers. Thus this
        function can return either a real or complex number.

        EXAMPLES::

            sage: R = RealField(30)
            sage: a = R('1.23456')
            sage: a^20
            67.646297
            sage: a^a
            1.2971115
            sage: b = R(-1)
            sage: b^(1/2)
            -8.7055157e-10 + 1.0000000*I

        We raise a real number to a symbolic object::

            sage: x, y = var('x,y')                                                     # needs sage.symbolic
            sage: 1.5^x                                                                 # needs sage.symbolic
            1.50000000000000^x
            sage: -2.3^(x+y^3+sin(x))                                                   # needs sage.symbolic
            -2.30000000000000^(y^3 + x + sin(x))

        TESTS:

        We see that :issue:`10736` is fixed::

            sage: 16^0.5
            4.00000000000000
            sage: int(16)^0.5
            4.00000000000000
            sage: (1/2)^2.0
            0.250000000000000
            sage: [n^(1.5) for n in range(10)]
            [0.000000000000000, 1.00000000000000, 2.82842712474619, 5.19615242270663,
             8.00000000000000, 11.1803398874989, 14.6969384566991, 18.5202591774521,
             22.6274169979695, 27.0000000000000]
            sage: int(-2)^(0.333333)
            0.629961522017056 + 1.09112272417509*I
            sage: int(0)^(1.0)
            0.000000000000000
            sage: int(0)^(0.0)
            1.00000000000000
        """
        cdef RealNumber base, x
        cdef mpfr_rnd_t rounding_mode
        if isinstance(self, RealNumber):
            base = <RealNumber>self
            rounding_mode = (<RealField_class>base._parent).rnd
            x = base._new()
            sig_on()
            if is_small_python_int(exponent):
                mpfr_pow_si(x.value, base.value, exponent, rounding_mode)
            elif isinstance(exponent, Integer):
                mpfr_pow_z(x.value, base.value, (<Integer>exponent).value, rounding_mode)
            elif isinstance(exponent, RealNumber):
                mpfr_pow(x.value, base.value, (<RealNumber>exponent).value, rounding_mode)
            else:
                sig_off()
                return bin_op(self, exponent, operator.pow)
            sig_off()
            if mpfr_nan_p(x.value):
                return base._complex_number_() ** exponent
            return x
        else:
            return bin_op(self, exponent, operator.pow)

    def log(self, base=None):
        """
        Return the logarithm of ``self`` to the ``base``.

        EXAMPLES::

            sage: R = RealField()
            sage: R(2).log()
            0.693147180559945
            sage: log(RR(2))
            0.693147180559945
            sage: log(RR(2), "e")
            0.693147180559945
            sage: log(RR(2), e)                                                         # needs sage.symbolic
            0.693147180559945

        ::

            sage: r = R(-1); r.log()
            3.14159265358979*I
            sage: log(RR(-1), e)                                                        # needs sage.symbolic
            3.14159265358979*I
            sage: r.log(2)
            4.53236014182719*I

        For the error value NaN (Not A Number), log will return NaN::

            sage: r = R(NaN); r.log()                                                   # needs sage.symbolic
            NaN
        """
        if mpfr_nan_p(self.value):
            return self

        cdef RealNumber x
        if self < 0:
            if base is None or base == 'e':
                return self._complex_number_().log()
            else:
                return self._complex_number_().log(base)
        if base is None or base == 'e':
            x = self._new()
            if (<RealField_class>self._parent)._prec > SIG_PREC_THRESHOLD:
                sig_on()
            mpfr_log(x.value, self.value, (<RealField_class>self._parent).rnd)
            if (<RealField_class>self._parent)._prec > SIG_PREC_THRESHOLD:
                sig_off()
            return x
        elif base == 10:
            return self.log10()
        elif base == 2:
            return self.log2()
        else:
            return self.log() / (self.parent()(base)).log()

    def log2(self):
        """
        Return log to the base 2 of ``self``.

        EXAMPLES::

            sage: r = 16.0
            sage: r.log2()
            4.00000000000000

        ::

            sage: r = 31.9; r.log2()
            4.99548451887751

        ::

            sage: r = 0.0
            sage: r.log2()
            -infinity

        ::

            sage: r = -3.0; r.log2()
            1.58496250072116 + 4.53236014182719*I
        """
        cdef RealNumber x
        if self < 0:
            return self._complex_number_().log(2)
        x = self._new()
        if (<RealField_class>self._parent)._prec > SIG_PREC_THRESHOLD:
            sig_on()
        mpfr_log2(x.value, self.value, (<RealField_class>self._parent).rnd)
        if (<RealField_class>self._parent)._prec > SIG_PREC_THRESHOLD:
            sig_off()
        return x

    def log10(self):
        """
        Return log to the base 10 of ``self``.

        EXAMPLES::

            sage: r = 16.0; r.log10()
            1.20411998265592
            sage: r.log() / log(10.0)
            1.20411998265592

        ::

            sage: r = 39.9; r.log10()
            1.60097289568675

        ::

            sage: r = 0.0
            sage: r.log10()
            -infinity

        ::

            sage: r = -1.0
            sage: r.log10()
            1.36437635384184*I
        """
        cdef RealNumber x
        if self < 0:
            return self._complex_number_().log(10)
        x = self._new()
        if (<RealField_class>self._parent)._prec > SIG_PREC_THRESHOLD:
            sig_on()
        mpfr_log10(x.value, self.value, (<RealField_class>self._parent).rnd)
        if (<RealField_class>self._parent)._prec > SIG_PREC_THRESHOLD:
            sig_off()
        return x

    def log1p(self):
        """
        Return log base `e` of ``1 + self``.

        EXAMPLES::

            sage: r = 15.0; r.log1p()
            2.77258872223978
            sage: (r+1).log()
            2.77258872223978

        For small values, this is more accurate than computing
        ``log(1 + self)`` directly, as it avoids cancellation issues::

            sage: r = 3e-10
            sage: r.log1p()
            2.99999999955000e-10
            sage: (1+r).log()
            3.00000024777111e-10
            sage: r100 = RealField(100)(r)
            sage: (1+r100).log()
            2.9999999995500000000978021372e-10

        ::

            sage: r = 38.9; r.log1p()
            3.68637632389582

        ::

            sage: r = -1.0
            sage: r.log1p()
            -infinity

        ::

            sage: r = -2.0
            sage: r.log1p()
            3.14159265358979*I
        """
        cdef RealNumber x
        if self < -1:
            return (self+1.0)._complex_number_().log()
        x = self._new()
        if (<RealField_class>self._parent)._prec > SIG_PREC_THRESHOLD:
            sig_on()
        mpfr_log1p(x.value, self.value, (<RealField_class>self._parent).rnd)
        if (<RealField_class>self._parent)._prec > SIG_PREC_THRESHOLD:
            sig_off()
        return x

    def exp(self):
        r"""
        Return `e^\mathtt{self}`.

        EXAMPLES::

            sage: r = 0.0
            sage: r.exp()
            1.00000000000000

        ::

            sage: r = 32.3
            sage: a = r.exp(); a
            1.06588847274864e14
            sage: a.log()
            32.3000000000000

        ::

            sage: r = -32.3
            sage: r.exp()
            9.38184458849869e-15
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_exp(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def exp2(self):
        r"""
        Return `2^\mathtt{self}`.

        EXAMPLES::

            sage: r = 0.0
            sage: r.exp2()
            1.00000000000000

        ::

            sage: r = 32.0
            sage: r.exp2()
            4.29496729600000e9

        ::

            sage: r = -32.3
            sage: r.exp2()
            1.89117248253021e-10
        """
        cdef RealNumber x = self._new()
        if (<RealField_class>self._parent)._prec > SIG_PREC_THRESHOLD:
            sig_on()
        mpfr_exp2(x.value, self.value, (<RealField_class>self._parent).rnd)
        if (<RealField_class>self._parent)._prec > SIG_PREC_THRESHOLD:
            sig_off()
        return x

    def exp10(self):
        r"""
        Return `10^\mathtt{self}`.

        EXAMPLES::

            sage: r = 0.0
            sage: r.exp10()
            1.00000000000000

        ::

            sage: r = 32.0
            sage: r.exp10()
            1.00000000000000e32

        ::

            sage: r = -32.3
            sage: r.exp10()
            5.01187233627276e-33
        """
        cdef RealNumber x = self._new()
        if (<RealField_class>self._parent)._prec > SIG_PREC_THRESHOLD:
            sig_on()
        mpfr_exp10(x.value, self.value, (<RealField_class>self._parent).rnd)
        if (<RealField_class>self._parent)._prec > SIG_PREC_THRESHOLD:
            sig_off()
        return x

    def expm1(self):
        r"""
        Return `e^\mathtt{self}-1`, avoiding cancellation near 0.

        EXAMPLES::

            sage: r = 1.0
            sage: r.expm1()
            1.71828182845905

        ::

            sage: r = 1e-16
            sage: exp(r)-1
            0.000000000000000
            sage: r.expm1()
            1.00000000000000e-16
        """
        cdef RealNumber x = self._new()
        if (<RealField_class>self._parent)._prec > SIG_PREC_THRESHOLD:
            sig_on()
        mpfr_expm1(x.value, self.value, (<RealField_class>self._parent).rnd)
        if (<RealField_class>self._parent)._prec > SIG_PREC_THRESHOLD:
            sig_off()
        return x

    def eint(self):
        """
        Return the exponential integral of this number.

        EXAMPLES::

            sage: r = 1.0
            sage: r.eint()
            1.89511781635594

        ::

            sage: r = -1.0
            sage: r.eint()
            -0.219383934395520
        """
        cdef RealNumber x = self._new()
        if (<RealField_class>self._parent)._prec > SIG_PREC_THRESHOLD:
            sig_on()
        mpfr_eint(x.value, self.value, (<RealField_class>self._parent).rnd)
        if (<RealField_class>self._parent)._prec > SIG_PREC_THRESHOLD:
            sig_on()
        return x

    def cos(self):
        """
        Return the cosine of ``self``.

        EXAMPLES::

            sage: t=RR.pi()/2
            sage: t.cos()
            6.12323399573677e-17
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_cos(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def sin(self):
        """
        Return the sine of ``self``.

        EXAMPLES::

            sage: R = RealField(100)
            sage: R(2).sin()
            0.90929742682568169539601986591
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_sin(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def tan(self):
        """
        Return the tangent of ``self``.

        EXAMPLES::

            sage: q = RR.pi()/3
            sage: q.tan()
            1.73205080756888
            sage: q = RR.pi()/6
            sage: q.tan()
            0.577350269189626
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_tan(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def sincos(self):
        """
        Return a pair consisting of the sine and cosine of ``self``.

        EXAMPLES::

            sage: R = RealField()
            sage: t = R.pi()/6
            sage: t.sincos()
            (0.500000000000000, 0.866025403784439)
        """
        cdef RealNumber x, y
        x = self._new()
        y = self._new()
        sig_on()
        mpfr_sin_cos(x.value, y.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x, y

    def arccos(self):
        """
        Return the inverse cosine of ``self``.

        EXAMPLES::

            sage: q = RR.pi()/3
            sage: i = q.cos()
            sage: i.arccos() == q
            True
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_acos(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def arcsin(self):
        """
        Return the inverse sine of ``self``.

        EXAMPLES::

            sage: q = RR.pi()/5
            sage: i = q.sin()
            sage: i.arcsin() == q
            True
            sage: i.arcsin() - q
            0.000000000000000
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_asin(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def arctan(self):
        """
        Return the inverse tangent of ``self``.

        EXAMPLES::

            sage: q = RR.pi()/5
            sage: i = q.tan()
            sage: i.arctan() == q
            True
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_atan(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def cosh(self):
        """
        Return the hyperbolic cosine of ``self``.

        EXAMPLES::

            sage: q = RR.pi()/12
            sage: q.cosh()
            1.03446564009551
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_cosh(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def sinh(self):
        """
        Return the hyperbolic sine of ``self``.

        EXAMPLES::

            sage: q = RR.pi()/12
            sage: q.sinh()
            0.264800227602271
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_sinh(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def tanh(self):
        """
        Return the hyperbolic tangent of ``self``.

        EXAMPLES::

            sage: q = RR.pi()/11
            sage: q.tanh()
            0.278079429295850
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_tanh(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def coth(self):
        """
        Return the hyperbolic cotangent of ``self``.

        EXAMPLES::

            sage: RealField(100)(2).coth()
            1.0373147207275480958778097648
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_coth(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def arccoth(self):
        """
        Return the inverse hyperbolic cotangent of ``self``.

        EXAMPLES::

            sage: q = RR.pi()/5
            sage: i = q.coth()
            sage: i.arccoth() == q
            True
        """
        return (~self).arctanh()

    def cot(self):
        """
        Return the cotangent of ``self``.

        EXAMPLES::

            sage: RealField(100)(2).cot()
            -0.45765755436028576375027741043
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_cot(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def csch(self):
        """
        Return the hyperbolic cosecant of ``self``.

        EXAMPLES::

            sage: RealField(100)(2).csch()
            0.27572056477178320775835148216
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_csch(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def arccsch(self):
        """
        Return the inverse hyperbolic cosecant of ``self``.

        EXAMPLES::

            sage: i = RR.pi()/5
            sage: q = i.csch()
            sage: q.arccsch() == i
            True
        """
        return (~self).arcsinh()

    def csc(self):
        """
        Return the cosecant of ``self``.

        EXAMPLES::

            sage: RealField(100)(2).csc()
            1.0997501702946164667566973970
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_csc(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def sech(self):
        """
        Return the hyperbolic secant of ``self``.

        EXAMPLES::

            sage: RealField(100)(2).sech()
            0.26580222883407969212086273982
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_sech(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def arcsech(self):
        """
        Return the inverse hyperbolic secant of ``self``.

        EXAMPLES::

            sage: i = RR.pi()/3
            sage: q = i.sech()
            sage: q.arcsech() == i
            True
        """
        return (~self).arccosh()

    def sec(self):
        """
        Return the secant of this number.

        EXAMPLES::

            sage: RealField(100)(2).sec()
            -2.4029979617223809897546004014
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_sec(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def arccosh(self):
        """
        Return the hyperbolic inverse cosine of ``self``.

        EXAMPLES::

            sage: q = RR.pi()/2
            sage: i = q.cosh() ; i
            2.50917847865806
            sage: q == i.arccosh()
            True
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_acosh(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def arcsinh(self):
        """
        Return the hyperbolic inverse sine of ``self``.

        EXAMPLES::

            sage: q = RR.pi()/7
            sage: i = q.sinh() ; i
            0.464017630492991
            sage: i.arcsinh() - q
            0.000000000000000
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_asinh(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def arctanh(self):
        """
        Return the hyperbolic inverse tangent of ``self``.

        EXAMPLES::

            sage: q = RR.pi()/7
            sage: i = q.tanh() ; i
            0.420911241048535
            sage: i.arctanh() - q
            0.000000000000000
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_atanh(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def agm(self, other):
        r"""
        Return the arithmetic-geometric mean of ``self`` and ``other``.

        The arithmetic-geometric mean is the common limit of the sequences
        `u_n` and `v_n`, where `u_0` is ``self``,
        `v_0` is other, `u_{n+1}` is the arithmetic mean
        of `u_n` and `v_n`, and `v_{n+1}` is the
        geometric mean of `u_n` and `v_n`. If any operand is negative, the
        return value is ``NaN``.

        INPUT:

        - ``right`` -- another real number

        OUTPUT: the AGM of ``self`` and ``other``

        EXAMPLES::

            sage: a = 1.5
            sage: b = 2.5
            sage: a.agm(b)
            1.96811775182478
            sage: RealField(200)(a).agm(b)
            1.9681177518247777389894630877503739489139488203685819712291
            sage: a.agm(100)
            28.1189391225320

        The AGM always lies between the geometric and arithmetic mean::

            sage: sqrt(a*b) < a.agm(b) < (a+b)/2
            True

        It is, of course, symmetric::

            sage: b.agm(a)
            1.96811775182478

        and satisfies the relation `AGM(ra, rb) = r AGM(a, b)`::

            sage: (2*a).agm(2*b) / 2
            1.96811775182478
            sage: (3*a).agm(3*b) / 3
            1.96811775182478

        It is also related to the elliptic integral

        .. MATH::

            \int_0^{\pi/2} \frac{d\theta}{\sqrt{1-m\sin^2\theta}}.

        ::

            sage: m = (a-b)^2/(a+b)^2
            sage: E = numerical_integral(1/sqrt(1-m*sin(x)^2), 0, RR.pi()/2)[0]         # needs sage.symbolic
            sage: RR.pi()/4 * (a+b)/E                                                   # needs sage.symbolic
            1.96811775182478

        TESTS::

            sage: 1.5.agm(0)
            0.000000000000000
        """
        cdef RealNumber x, _other
        if isinstance(other, RealNumber) and ((<Element>other)._parent is self._parent):
            _other = <RealNumber>other
        else:
            _other = self._parent(other)

        x = self._new()
        if (<RealField_class>self._parent)._prec > 10000:
            sig_on()
        mpfr_agm(x.value, self.value, _other.value, (<RealField_class>self._parent).rnd)
        if (<RealField_class>self._parent)._prec > 10000:
            sig_off()
        return x

    def erf(self):
        """
        Return the value of the error function on ``self``.

        EXAMPLES::

            sage: R = RealField(53)
            sage: R(2).erf()
            0.995322265018953
            sage: R(6).erf()
            1.00000000000000
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_erf(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def erfc(self):
        r"""
        Return the value of the complementary error function on ``self``,
        i.e., `1-\mathtt{erf}(\mathtt{self})`.

        EXAMPLES::

            sage: R = RealField(53)
            sage: R(2).erfc()
            0.00467773498104727
            sage: R(6).erfc()
            2.15197367124989e-17
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_erfc(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def j0(self):
        """
        Return the value of the Bessel `J` function of order 0 at ``self``.

        EXAMPLES::

            sage: R = RealField(53)
            sage: R(2).j0()
            0.223890779141236
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_j0(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def j1(self):
        """
        Return the value of the Bessel `J` function of order 1 at ``self``.

        EXAMPLES::

            sage: R = RealField(53)
            sage: R(2).j1()
            0.576724807756873
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_j1(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def jn(self, long n):
        """
        Return the value of the Bessel `J` function of order `n` at ``self``.

        EXAMPLES::

            sage: R = RealField(53)
            sage: R(2).jn(3)
            0.128943249474402
            sage: R(2).jn(-17)
            -2.65930780516787e-15
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_jn(x.value, n, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def y0(self):
        """
        Return the value of the Bessel `Y` function of order 0 at ``self``.

        EXAMPLES::

            sage: R = RealField(53)
            sage: R(2).y0()
            0.510375672649745
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_y0(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def y1(self):
        """
        Return the value of the Bessel `Y` function of order 1 at ``self``.

        EXAMPLES::

            sage: R = RealField(53)
            sage: R(2).y1()
            -0.107032431540938
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_y1(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def yn(self, long n):
        """
        Return the value of the Bessel `Y` function of order `n` at ``self``.

        EXAMPLES::

            sage: R = RealField(53)
            sage: R(2).yn(3)
            -1.12778377684043
            sage: R(2).yn(-17)
            7.09038821729481e12
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_yn(x.value, n, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def gamma(self):
        """
        Return the value of the Euler gamma function on ``self``.

        EXAMPLES::

            sage: R = RealField()
            sage: R(6).gamma()
            120.000000000000
            sage: R(1.5).gamma()
            0.886226925452758
        """
        cdef RealNumber x = self._new()
        if (<RealField_class>self._parent)._prec > SIG_PREC_THRESHOLD:
            sig_on()
        mpfr_gamma(x.value, self.value, (<RealField_class>self._parent).rnd)
        if (<RealField_class>self._parent)._prec > SIG_PREC_THRESHOLD:
            sig_off()
        return x

    def log_gamma(self):
        """
        Return the principal branch of the log gamma of ``self``.

        Note that this is not in general equal to log(gamma(``self``))
        for negative input.

        EXAMPLES::

            sage: R = RealField(53)
            sage: R(6).log_gamma()
            4.78749174278205
            sage: R(1e10).log_gamma()
            2.20258509288811e11
            sage: log_gamma(-2.1)
            1.53171380819509 - 9.42477796076938*I
            sage: log(gamma(-1.1)) == log_gamma(-1.1)
            False
        """
        cdef RealNumber x = self._new()
        parent = (<RealField_class>self._parent)
        if not mpfr_sgn(self.value) < 0:
            if parent._prec > SIG_PREC_THRESHOLD:
                sig_on()
            mpfr_lngamma(x.value, self.value, parent.rnd)
            if parent._prec > SIG_PREC_THRESHOLD:
                sig_off()
            return x
        from sage.libs.mpmath.utils import call
        from mpmath import loggamma
        return call(loggamma, mpfr_to_mpfval(self.value), parent=parent)

    def zeta(self):
        r"""
        Return the Riemann zeta function evaluated at this real number.

        .. NOTE::

           PARI is vastly more efficient at computing the Riemann zeta
           function. See the example below for how to use it.

        EXAMPLES::

            sage: R = RealField()
            sage: R(2).zeta()
            1.64493406684823
            sage: R.pi()^2/6
            1.64493406684823
            sage: R(-2).zeta()
            0.000000000000000
            sage: R(1).zeta()
            +infinity

        Computing zeta using PARI is much more efficient in difficult
        cases. Here's how to compute zeta with at least a given precision::

            sage: z = pari(2).zeta(precision=53); z                                     # needs sage.libs.pari
            1.64493406684823
            sage: pari(2).zeta(precision=128).sage().prec()                             # needs sage.libs.pari
            128
            sage: pari(2).zeta(precision=65).sage().prec()                              # needs sage.libs.pari
            128                                                # 64-bit
            96                                                 # 32-bit

        Note that the number of bits of precision in the constructor only
        effects the internal precision of the pari number, which is rounded
        up to the nearest multiple of 32 or 64. To increase the number of
        digits that gets displayed you must use
        ``pari.set_real_precision``.

        ::

            sage: type(z)                                                               # needs sage.libs.pari
            <class 'cypari2.gen.Gen'>
            sage: R(z)                                                                  # needs sage.libs.pari
            1.64493406684823
        """
        cdef RealNumber x = self._new()
        sig_on()
        mpfr_zeta(x.value, self.value, (<RealField_class>self._parent).rnd)
        sig_off()
        return x

    def algebraic_dependency(self, n):
        """
        Return a polynomial of degree at most `n` which is
        approximately satisfied by this number.

        .. NOTE::

            The resulting polynomial need not be irreducible, and indeed
            usually won't be if this number is a good approximation to an
            algebraic number of degree less than `n`.

        ALGORITHM:

        Uses the PARI C-library :pari:`algdep` command.

        EXAMPLES::

            sage: r = sqrt(2.0); r
            1.41421356237310
            sage: r.algebraic_dependency(5)
            x^2 - 2
        """
        return sage.arith.misc.algdep(self, n)

    algdep = algebraic_dependency

    def nth_root(self, int n, int algorithm=0):
        r"""
        Return an `n`-th root of ``self``.

        INPUT:

        - ``n`` -- a positive number, rounded down to the nearest integer;
          note that `n` should be less than ``sys.maxsize``

        - ``algorithm`` -- set this to 1 to call mpfr directly, set this to 2
          to use interval arithmetic and logarithms, or leave it at the default
          of 0 to choose the algorithm which is estimated to be faster

        AUTHORS:

        - Carl Witty (2007-10)

        EXAMPLES::

            sage: R = RealField()
            sage: R(8).nth_root(3)
            2.00000000000000
            sage: R(8).nth_root(3.7)    # illustrate rounding down
            2.00000000000000
            sage: R(-8).nth_root(3)
            -2.00000000000000
            sage: R(0).nth_root(3)
            0.000000000000000
            sage: R(32).nth_root(-1)
            Traceback (most recent call last):
            ...
            ValueError: n must be positive
            sage: R(32).nth_root(1.0)
            32.0000000000000
            sage: R(4).nth_root(4)
            1.41421356237310
            sage: R(4).nth_root(40)
            1.03526492384138
            sage: R(4).nth_root(400)
            1.00347174850950
            sage: R(4).nth_root(4000)
            1.00034663365385
            sage: R(4).nth_root(4000000)
            1.00000034657365
            sage: R(-27).nth_root(3)
            -3.00000000000000
            sage: R(-4).nth_root(3999999)
            -1.00000034657374

        Note that for negative numbers, any even root throws an exception::

            sage: R(-2).nth_root(6)
            Traceback (most recent call last):
            ...
            ValueError: taking an even root of a negative number

        The `n`-th root of 0 is defined to be 0, for any `n`::

            sage: R(0).nth_root(6)
            0.000000000000000
            sage: R(0).nth_root(7)
            0.000000000000000

        TESTS:

        The old and new algorithms should give exactly the same results in
        all cases::

            sage: def check(x, n):
            ....:     answers = []
            ....:     for sign in (1, -1):
            ....:         if is_even(n) and sign == -1:
            ....:             continue
            ....:         for rounding in ('RNDN', 'RNDD', 'RNDU', 'RNDZ'):
            ....:             fld = RealField(x.prec(), rnd=rounding)
            ....:             fx = fld(sign * x)
            ....:             alg_mpfr = fx.nth_root(n, algorithm=1)
            ....:             alg_mpfi = fx.nth_root(n, algorithm=2)
            ....:             assert(alg_mpfr == alg_mpfi)
            ....:             if sign == 1: answers.append(alg_mpfr)
            ....:     return answers

        Check some perfect powers (and nearby numbers)::

            sage: check(16.0, 4)
            [2.00000000000000, 2.00000000000000, 2.00000000000000, 2.00000000000000]
            sage: check((16.0).nextabove(), 4)
            [2.00000000000000, 2.00000000000000, 2.00000000000001, 2.00000000000000]
            sage: check((16.0).nextbelow(), 4)
            [2.00000000000000, 1.99999999999999, 2.00000000000000, 1.99999999999999]
            sage: check(((9.0 * 256)^7), 7)
            [2304.00000000000, 2304.00000000000, 2304.00000000000, 2304.00000000000]
            sage: check(((9.0 * 256)^7).nextabove(), 7)
            [2304.00000000000, 2304.00000000000, 2304.00000000001, 2304.00000000000]
            sage: check(((9.0 * 256)^7).nextbelow(), 7)
            [2304.00000000000, 2303.99999999999, 2304.00000000000, 2303.99999999999]
            sage: check(((5.0 / 512)^17), 17)
            [0.00976562500000000, 0.00976562500000000, 0.00976562500000000, 0.00976562500000000]
            sage: check(((5.0 / 512)^17).nextabove(), 17)
            [0.00976562500000000, 0.00976562500000000, 0.00976562500000001, 0.00976562500000000]
            sage: check(((5.0 / 512)^17).nextbelow(), 17)
            [0.00976562500000000, 0.00976562499999999, 0.00976562500000000, 0.00976562499999999]

        And check some non-perfect powers::

            sage: check(2.0, 3)
            [1.25992104989487, 1.25992104989487, 1.25992104989488, 1.25992104989487]
            sage: check(2.0, 4)
            [1.18920711500272, 1.18920711500272, 1.18920711500273, 1.18920711500272]
            sage: check(2.0, 5)
            [1.14869835499704, 1.14869835499703, 1.14869835499704, 1.14869835499703]

        And some different precisions::

            sage: check(RealField(20)(22/7), 19)
            [1.0621, 1.0621, 1.0622, 1.0621]
            sage: check(RealField(200)(e), 4)                                           # needs sage.symbolic
            [1.2840254166877414840734205680624364583362808652814630892175,
             1.2840254166877414840734205680624364583362808652814630892175,
             1.2840254166877414840734205680624364583362808652814630892176,
             1.2840254166877414840734205680624364583362808652814630892175]

        Check that :issue:`12105` is fixed::

            sage: RealField(53)(0.05).nth_root(7 * 10^8)
            0.999999995720382
        """
        if n <= 0:
            raise ValueError("n must be positive")

        cdef int odd = (n & 1)

        cdef int sgn = mpfr_sgn(self.value)

        if sgn < 0 and not odd:
            raise ValueError("taking an even root of a negative number")

        if sgn == 0 or n == 1 or not mpfr_number_p(self.value):
            return self

        cdef RealField_class fld = <RealField_class>self._parent

        if algorithm == 0 and n <= 10000 / fld._prec:
            # This is a rough estimate for when it is probably
            # faster to call mpfr directly.  (This is a pretty
            # good estimate on one particular machine, a
            # Core 2 Duo in 32-bit mode, but has not been tested
            # on other machines.)
            algorithm = 1

        cdef RealNumber x

        if algorithm == 1:
            x = self._new()
            sig_on()
            mpfr_rootn_ui(x.value, self.value, n, (<RealField_class>self._parent).rnd)
            sig_off()
            return x

        cdef mpfr_rnd_t rnd = (<RealField_class>self._parent).rnd

        cdef Integer mantissa
        cdef mp_exp_t exponent
        cdef int pow2
        cdef int exact

        if rnd != MPFR_RNDN:
            # We are going to implement nth_root using interval
            # arithmetic.  To guarantee correct rounding, we will
            # increase the precision of the interval arithmetic until
            # the resulting interval is small enough...until the
            # interval is entirely within the interval represented
            # by a single floating-point number.

            # This always works, unless the correct answer is exactly
            # on the boundary point between the intervals of two
            # floating-point numbers.  In round-to-nearest mode, this
            # is impossible; the boundary points are the
            # numbers which can be exactly represented as precision-{k+1}
            # floating-point numbers, but not precision-{k} numbers.
            # A precision-{k} floating-point number cannot be a perfect
            # n'th power (n >= 2) of such a number.

            # However, in the directed rounding modes, the boundary points
            # are the floating-point numbers themselves.  So in a
            # directed rounding mode, we need to check whether this
            # floating-point number is a perfect n'th power.

            # Suppose this number is (a * 2^k)^n, for odd integer a
            # and arbitrary integer k.  Then this number is
            # (a^n) * 2^(k*n), where a^n is odd.

            # We start by extracting the mantissa and exponent of (the
            # absolute value of) this number.

            mantissa = Integer()
            sig_on()
            exponent = mpfr_get_z_exp(mantissa.value, self.value)
            sig_off()
            mpz_abs(mantissa.value, mantissa.value)

            # Now, we want to divide out any powers of two in mantissa,
            # leaving it as an odd number.

            sig_on()
            pow2 = mpz_scan1(mantissa.value, 0)
            sig_off()

            if pow2 > 0:
                exponent = exponent + pow2
                sig_on()
                mpz_fdiv_q_2exp(mantissa.value, mantissa.value, pow2)
                sig_off()

            # Our floating-point number is equal to mantissa * 2^exponent,
            # and we know that mantissa is odd.

            if exponent % n == 0:
                # The exponent is a multiple of n, so it's possible that
                # we have a perfect power.  Now we need to check the
                # mantissa.

                sig_on()
                exact = mpz_root(mantissa.value, mantissa.value, n)
                sig_off()

                if exact:
                    # Yes, we are a perfect power.  We've replaced mantissa
                    # with its n'th root, so we can just build
                    # the resulting floating-point number.

                    x = self._new()

                    sig_on()
                    mpfr_set_z(x.value, mantissa.value, MPFR_RNDN)
                    sig_off()
                    mpfr_mul_2si(x.value, x.value, exponent / n, MPFR_RNDN)
                    if sgn < 0:
                        mpfr_neg(x.value, x.value, MPFR_RNDN)

                    return x

        # If we got here, then we're not a perfect power of a boundary
        # point, so it's safe to use the interval arithmetic technique.

        from sage.rings.real_mpfi import RealIntervalField

        cdef mpfr_prec_t prec = fld._prec + 10

        cdef RealNumber lower
        cdef RealNumber upper

        while True:
            ifld = RealIntervalField(prec)
            intv = ifld(self)
            if sgn < 0:
                intv = -intv
            intv = (intv.log() / n).exp()
            if sgn < 0:
                intv = -intv
            lower = fld(intv.lower())
            upper = fld(intv.upper())

            if mpfr_equal_p(lower.value, upper.value):
                # Yes, we found the answer
                return lower

            prec = prec + 20

cdef class RealLiteral(RealNumber):
    """
    Real literals are created in preparsing and provide a way to allow
    casting into higher precision rings.
    """

    cdef readonly literal
    cdef readonly int base

    def __init__(self, RealField_class parent, x=0, int base=10):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: RealField(200)(float(1.3))
            1.3000000000000000444089209850062616169452667236328125000000
            sage: RealField(200)(1.3)
            1.3000000000000000000000000000000000000000000000000000000000
            sage: 1.3 + 1.2
            2.50000000000000
            sage: RR(1_0000.000000000000000000000000000000000000)
            10000.0000000000
        """
        RealNumber.__init__(self, parent, x, base)
        if isinstance(x, str):
            self.base = base
            self.literal = x.replace('_', '')

    def __neg__(self):
        """
        Return the negative of ``self``.

        EXAMPLES::

            sage: RealField(300)(-1.2)
            -1.20000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
            sage: RealField(300)(-(-1.2))
            1.20000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        """
        if self.literal is not None and self.literal[0] == '-':
            return RealLiteral(self._parent, self.literal[1:], self.base)
        else:
            return RealLiteral(self._parent, '-'+self.literal, self.base)

    def numerical_approx(self, prec=None, digits=None, algorithm=None):
        """
        Change the precision of ``self`` to ``prec`` bits
        or ``digits`` decimal digits.

        INPUT:

        - ``prec`` -- precision in bits

        - ``digits`` -- precision in decimal digits (only used if
          ``prec`` is not given)

        - ``algorithm`` -- ignored for real numbers

        If neither ``prec`` nor ``digits`` is given, the default
        precision is 53 bits (roughly 16 digits).

        OUTPUT: a ``RealNumber`` with the given precision

        EXAMPLES::

            sage: (1.3).numerical_approx()
            1.30000000000000
            sage: n(1.3, 120)
            1.3000000000000000000000000000000000

        Compare with::

            sage: RealField(120)(RR(13/10))
            1.3000000000000000444089209850062616
            sage: n(RR(13/10), 120)
            Traceback (most recent call last):
            ...
            TypeError: cannot approximate to a precision of 120 bits, use at most 53 bits

        The result is a non-literal::

            sage: type(1.3)
            <class 'sage.rings.real_mpfr.RealLiteral'>
            sage: type(n(1.3))
            <class 'sage.rings.real_mpfr.RealNumber'>
        """
        if prec is None:
            prec = digits_to_bits(digits)
        return RealField(prec)(self.literal)


RR = RealField()

RR_min_prec = RealField(MPFR_PREC_MIN)


def create_RealNumber(s, int base=10, int pad=0, rnd='RNDN', int min_prec=53):
    r"""
    Return the real number defined by the string ``s`` as an element of
    ``RealField(prec=n)``, where ``n`` potentially has slightly
    more (controlled by pad) bits than given by ``s``.

    INPUT:

    - ``s`` -- string that defines a real number (or
      something whose string representation defines a number)

    - ``base`` -- integer between 2 and 62

    - ``pad`` -- nonnegative integer

    - ``rnd`` -- rounding mode:

      - ``'RNDN'`` -- round to nearest
      - ``'RNDZ'`` -- round toward zero
      - ``'RNDD'`` -- round down
      - ``'RNDU'`` -- round up

    - ``min_prec`` -- number will have at least this many
      bits of precision, no matter what

    EXAMPLES::

        sage: RealNumber('2.3') # indirect doctest
        2.30000000000000
        sage: RealNumber(10)
        10.0000000000000
        sage: RealNumber('1.0000000000000000000000000000000000')
        1.000000000000000000000000000000000
        sage: RealField(200)(1.2)
        1.2000000000000000000000000000000000000000000000000000000000
        sage: (1.2).parent() is RR
        True

    We can use various bases::

        sage: RealNumber("10101e2",base=2)
        84.0000000000000
        sage: RealNumber("deadbeef", base=16)
        3.73592855900000e9
        sage: RealNumber("deadbeefxxx", base=16)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert 'deadbeefxxx' to a real number
        sage: RealNumber("z", base=36)
        35.0000000000000
        sage: RealNumber("AAA", base=37)
        14070.0000000000
        sage: RealNumber("aaa", base=37)
        50652.0000000000
        sage: RealNumber("3.4", base='foo')
        Traceback (most recent call last):
        ...
        TypeError: an integer is required
        sage: RealNumber("3.4", base=63)
        Traceback (most recent call last):
        ...
        ValueError: base (=63) must be an integer between 2 and 62

    The rounding mode is respected in all cases::

        sage: RealNumber("1.5", rnd="RNDU").parent()
        Real Field with 53 bits of precision and rounding RNDU
        sage: RealNumber("1.50000000000000000000000000000000000000", rnd="RNDU").parent()
        Real Field with 130 bits of precision and rounding RNDU

    TESTS::

        sage: RealNumber('.000000000000000000000000000000001').prec()
        53
        sage: RealNumber('-.000000000000000000000000000000001').prec()
        53

        sage: RealNumber('-.123456789123456789').prec()
        60
        sage: RealNumber('.123456789123456789').prec()
        60
        sage: RealNumber('0.123456789123456789').prec()
        60
        sage: RealNumber('00.123456789123456789').prec()
        60
        sage: RealNumber('123456789.123456789').prec()
        60

    Make sure we've rounded up ``log(10,2)`` enough to guarantee
    sufficient precision (:issue:`10164`)::

        sage: ks = 5*10**5, 10**6
        sage: all(RealNumber("1." + "0"*k +"1")-1 > 0 for k in ks)
        True
    """
    if not isinstance(s, str):
        s = str(s)

    # Check for a valid base
    if base < 2 or base > 62:
        raise ValueError(f"base (={base}) must be an integer between 2 and 62")

    if base == 10 and min_prec == 53 and len(s) <= 15 and rnd == "RNDN":
        R = RR
    else:
        # For bases 15 and up, treat 'e' as digit
        if base <= 14 and ('e' in s or 'E' in s):
            # Figure out the exponent
            index = max(s.find('e'), s.find('E'))
            mantissa = s[:index]
        else:
            mantissa = s

        # Find the first nonzero entry in rest
        sigfig_mantissa = mantissa.lstrip('-0.')
        sigfigs = len(sigfig_mantissa) - ('.' in sigfig_mantissa)

        if base == 10:
            # hard-code the common case
            bits = int(LOG_TEN_TWO_PLUS_EPSILON * sigfigs) + 1
        else:
            bits = int(math.log(base, 2) * 1.00001 * sigfigs) + 1

        R = RealField(prec=max(bits + pad, min_prec), rnd=rnd)

    return RealLiteral(R, s, base)


def is_RealNumber(x):
    """
    Return ``True`` if ``x`` is of type :class:`RealNumber`, meaning that it
    is an element of the MPFR real field with some precision.

    EXAMPLES::

        sage: from sage.rings.real_mpfr import is_RealNumber
        sage: is_RealNumber(2.5)
        doctest:warning...
        DeprecationWarning: The function is_RealNumber is deprecated;
        use 'isinstance(..., RealNumber)' instead.
        See https://github.com/sagemath/sage/issues/38128 for details.
        True
        sage: is_RealNumber(float(2.3))
        False
        sage: is_RealNumber(RDF(2))
        False
        sage: is_RealNumber(pi)                                                         # needs sage.symbolic
        False
    """
    from sage.misc.superseded import deprecation_cython
    deprecation_cython(38128,
                       "The function is_RealNumber is deprecated; "
                       "use 'isinstance(..., RealNumber)' instead.")
    return isinstance(x, RealNumber)


def __create__RealField_version0(prec, sci_not, rnd):
    """
    Create a :class:`RealField_class` by calling :func:`RealField()`.

    EXAMPLES::

        sage: sage.rings.real_mpfr.__create__RealField_version0(53, 0, 'RNDN')
        Real Field with 53 bits of precision
    """
    return RealField(prec, sci_not, rnd)


def __create__RealNumber_version0(parent, x, base=10):
    """
    Create a :class:`RealNumber`.

    EXAMPLES::

        sage: sage.rings.real_mpfr.__create__RealNumber_version0(RR, 22,base=3)
        22.0000000000000
    """
    return RealNumber(parent, x, base=base)


cdef class RRtoRR(Map):
    cpdef Element _call_(self, x):
        """
        EXAMPLES::

            sage: from sage.rings.real_mpfr import RRtoRR
            sage: R10 = RealField(10)
            sage: R100 = RealField(100)
            sage: f = RRtoRR(R100, R10)
            sage: a = R100(1.2)
            sage: f(a)
            1.2
            sage: g = f.section()
            sage: g
            Generic map:
              From: Real Field with 10 bits of precision
              To:   Real Field with 100 bits of precision
            sage: g(f(a)) # indirect doctest
            1.1992187500000000000000000000
            sage: b = R10(2).sqrt()
            sage: f(g(b))
            1.4
            sage: f(g(b)) == b
            True
        """
        cdef RealField_class parent = <RealField_class>self._codomain
        cdef RealNumber y = parent._new()
        if type(x) is RealLiteral:
            mpfr_set_str(y.value, str_to_bytes((<RealLiteral>x).literal),
                         (<RealLiteral>x).base, parent.rnd)
        else:
            mpfr_set(y.value, (<RealNumber>x).value, parent.rnd)
        return y

    def section(self):
        """
        EXAMPLES::

            sage: from sage.rings.real_mpfr import RRtoRR
            sage: R10 = RealField(10)
            sage: R100 = RealField(100)
            sage: f = RRtoRR(R100, R10)
            sage: f.section()
            Generic map:
              From: Real Field with 10 bits of precision
              To:   Real Field with 100 bits of precision
        """
        return RRtoRR(self._codomain, self.domain())

cdef class ZZtoRR(Map):
    cpdef Element _call_(self, x):
        """
        EXAMPLES::

            sage: from sage.rings.real_mpfr import ZZtoRR
            sage: f = ZZtoRR(ZZ, RealField(20))
            sage: f(123456789) # indirect doctest
            1.2346e8
        """
        cdef RealField_class parent = <RealField_class>self._codomain
        cdef RealNumber y = parent._new()
        mpfr_set_z(y.value, (<Integer>x).value, parent.rnd)
        return y

cdef class QQtoRR(Map):
    cpdef Element _call_(self, x):
        """
        EXAMPLES::

            sage: from sage.rings.real_mpfr import QQtoRR
            sage: f = QQtoRR(QQ, RealField(200))
            sage: f(-1/3) # indirect doctest
            -0.33333333333333333333333333333333333333333333333333333333333
        """
        cdef RealField_class parent = <RealField_class>self._codomain
        cdef RealNumber y = parent._new()
        mpfr_set_q(y.value, (<Rational>x).value, parent.rnd)
        return y

cdef class double_toRR(Map):
    cpdef Element _call_(self, x):
        """
        Take anything that can be converted to a double.

        EXAMPLES::

            sage: from sage.rings.real_mpfr import double_toRR
            sage: f = double_toRR(RDF, RealField(22))
            sage: f(RDF.pi()) # indirect doctest
            3.14159
            sage: f = double_toRR(RDF, RealField(200))
            sage: f(RDF.pi())
            3.1415926535897931159979634685441851615905761718750000000000
        """
        cdef RealField_class parent = <RealField_class>self._codomain
        cdef RealNumber y = parent._new()
        mpfr_set_d(y.value, x, parent.rnd)
        return y

cdef class int_toRR(Map):
    cpdef Element _call_(self, x):
        """
        Take Python int/long instances.

        EXAMPLES::

            sage: from sage.rings.real_mpfr import int_toRR
            sage: f = int_toRR(int, RR)
            sage: f(-10r) # indirect doctest
            -10.0000000000000
            sage: f(2^75)
            3.77789318629572e22

        Also accepts objects that can be converted to int/long::

            sage: R.<x> = ZZ[]
            sage: f = int_toRR(R, RR)
            sage: f(x-x+1)
            1.00000000000000
        """
        cdef RealField_class parent = <RealField_class>self._codomain
        cdef RealNumber y = parent._new()
        cdef int err = 0
        cdef long x_long
        cdef mpz_t x_mpz

        if not isinstance(x, int):
            x = int(x)

        integer_check_long_py(x, &x_long, &err)

        if not err:
            mpfr_set_si(y.value, x_long, parent.rnd)
        elif err == ERR_OVERFLOW:
            mpz_init(x_mpz)
            mpz_set_pylong(x_mpz, x)
            mpfr_set_z(y.value, x_mpz, parent.rnd)
            mpz_clear(x_mpz)
        else:
            # This should never happen
            raise TypeError("argument cannot be converted to a Python int/long")

        return y


# Hook into Rational
from sage.rings.rational import _register_real_number_class
_register_real_number_class(RealNumber)


# Support Python's numbers abstract base class
import numbers
numbers.Real.register(RealNumber)
