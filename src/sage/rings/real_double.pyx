r"""
Double precision floating point real numbers

EXAMPLES:

We create the real double vector space of dimension `3`::

    sage: V = RDF^3; V                                                                  # needs sage.modules
    Vector space of dimension 3 over Real Double Field

Notice that this space is unique::

    sage: V is RDF^3                                                                    # needs sage.modules
    True
    sage: V is FreeModule(RDF, 3)                                                       # needs sage.modules
    True
    sage: V is VectorSpace(RDF, 3)                                                      # needs sage.modules
    True

Also, you can instantly create a space of large dimension::

    sage: V = RDF^10000                                                                 # needs sage.modules

TESTS:

Test NumPy conversions::

    sage: RDF(1).__array_interface__
    {'typestr': '=f8'}
    sage: import numpy                                                                  # needs numpy
    sage: numpy.array([RDF.pi()]).dtype                                                 # needs numpy
    dtype('float64')
"""

# ****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cimport libc.math
from libc.string cimport memcpy
from cpython.object cimport *
from cpython.float cimport *

from sage.ext.stdsage cimport PY_NEW
from sage.cpython.python_debug cimport if_Py_TRACE_REFS_then_PyObject_INIT

import math

import sage.arith.misc
import sage.rings.integer
import sage.rings.rational

from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ

from sage.categories.morphism cimport Morphism
from sage.structure.coerce cimport is_numpy_type
from sage.misc.randstate cimport randstate, current_randstate
from sage.structure.richcmp cimport rich_to_bool
from sage.arith.constants cimport *

cimport gmpy2


new_gen_from_real_double_element = None


cdef class RealDoubleField_class(sage.rings.abc.RealDoubleField):
    """
    An approximation to the field of real numbers using double
    precision floating point numbers. Answers derived from calculations
    in this approximation may differ from what they would be if those
    calculations were performed in the true field of real numbers. This
    is due to the rounding errors inherent to finite precision
    calculations.

    EXAMPLES::

        sage: RR == RDF                                                                 # needs sage.rings.real_mpfr
        False
        sage: RDF == RealDoubleField()    # RDF is the shorthand
        True

    ::

        sage: RDF(1)
        1.0
        sage: RDF(2/3)
        0.6666666666666666

    A :exc:`TypeError` is raised if the coercion doesn't make sense::

        sage: RDF(QQ['x'].0)
        Traceback (most recent call last):
        ...
        TypeError: cannot convert nonconstant polynomial
        sage: RDF(QQ['x'](3))
        3.0

    One can convert back and forth between double precision real
    numbers and higher-precision ones, though of course there may be
    loss of precision::

        sage: # needs sage.rings.real_mpfr
        sage: a = RealField(200)(2).sqrt(); a
        1.4142135623730950488016887242096980785696718753769480731767
        sage: b = RDF(a); b
        1.4142135623730951
        sage: a.parent()(b)
        1.4142135623730951454746218587388284504413604736328125000000
        sage: a.parent()(b) == b
        True
        sage: b == RR(a)
        True

    TESTS::

        sage: RDF.is_finite()
        False
    """
    def __init__(self):
        """
        Initialize ``self``.

        TESTS::

            sage: R = RealDoubleField()
            sage: TestSuite(R).run()
        """
        from sage.categories.fields import Fields
        Field.__init__(self, self,
                       category=Fields().Infinite().Metric().Complete())
        self._populate_coercion_lists_(init_no_parent=True,
                                       convert_method_name='_real_double_')

    _element_constructor_ = RealDoubleElement

    def __reduce__(self):
        """
        For pickling.

        EXAMPLES::

            sage: loads(dumps(RDF)) is RDF
            True
        """
        return RealDoubleField, ()

    cpdef bint is_exact(self) except -2:
        """
        Return ``False``, because doubles are not exact.

        EXAMPLES::

            sage: RDF.is_exact()
            False
        """
        return False

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(RDF)  # indirect doctest
            \Bold{R}
        """
        return "\\Bold{R}"

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: sage_input(RDF, verify=True)
            # Verified
            RDF
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: RDF._sage_input_(SageInputBuilder(), False)
            {atomic:RDF}
        """
        return sib.name('RDF')

    def __repr__(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: RealDoubleField()  # indirect doctest
            Real Double Field
            sage: RDF
            Real Double Field
        """
        return "Real Double Field"

    def _repr_option(self, key):
        """
        Metadata about the :meth:`_repr_` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: RDF._repr_option('element_is_atomic')
            True
        """
        if key == 'element_is_atomic':
            return True
        return super()._repr_option(key)

    def __richcmp__(self, x, op):
        """
        Compare ``self`` to ``x``.

        EXAMPLES::

            sage: RDF == 5
            False
            sage: loads(dumps(RDF)) == RDF
            True
        """
        if isinstance(x, RealDoubleField_class):
            return rich_to_bool(op, 0)
        if op == Py_NE:
            return True
        return NotImplemented

    def construction(self):
        r"""
        Return the functorial construction of ``self``, namely, completion of
        the rational numbers with respect to the prime at `\infty`.

        Also preserves other information that makes this field unique (i.e.
        the Real Double Field).

        EXAMPLES::

            sage: c, S = RDF.construction(); S
            Rational Field
            sage: RDF == c(S)
            True
        """
        from sage.categories.pushout import CompletionFunctor
        return (CompletionFunctor(sage.rings.infinity.Infinity,
                                  53,
                                  {'type': 'RDF'}),
                sage.rings.rational_field.QQ)

    def complex_field(self):
        """
        Return the complex field with the same precision as ``self``, i.e.,
        the complex double field.

        EXAMPLES::

            sage: RDF.complex_field()                                                   # needs sage.rings.complex_double
            Complex Double Field
        """
        from sage.rings.complex_double import CDF
        return CDF

    def algebraic_closure(self):
        """
        Return the algebraic closure of ``self``, i.e., the complex double
        field.

        EXAMPLES::

            sage: RDF.algebraic_closure()                                               # needs sage.rings.complex_double
            Complex Double Field
        """
        from sage.rings.complex_double import CDF
        return CDF

    cpdef _coerce_map_from_(self, S):
        """
        Canonical coercion of ``S`` to the real double field.

        The rings that canonically coerce to the real double field are:

        - the real double field itself
        - int, long, integer, and rational rings
        - numpy integers and floatings
        - the real lazy field
        - the MPFR real field with at least 53 bits of precision

        EXAMPLES::

            sage: RDF.coerce(5)  # indirect doctest
            5.0
            sage: RDF.coerce(9499294r)
            9499294.0
            sage: RDF.coerce(61/3)
            20.333333333333332
            sage: parent(RDF(3) + CDF(5))                                               # needs sage.rings.complex_double
            Complex Double Field
            sage: parent(CDF(5) + RDF(3))                                               # needs sage.rings.complex_double
            Complex Double Field
            sage: CDF.gen(0) + 5.0                                                      # needs sage.rings.complex_double
            5.0 + 1.0*I
            sage: RLF(2/3) + RDF(1)
            1.6666666666666665

            sage: import numpy                                                          # needs numpy
            sage: RDF.coerce(numpy.int8('1'))                                           # needs numpy
            1.0
            sage: RDF.coerce(numpy.float64('1'))                                        # needs numpy
            1.0

            sage: RDF.coerce(pi)                                                        # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Symbolic Ring to Real Double Field

        Test that :issue:`15695` is fixed (see also :issue:`18076`)::

            sage: 1j + numpy.float64(2)                                                 # needs numpy
            2.00000000000000 + 1.00000000000000*I
            sage: parent(_)                                                             # needs numpy
            Complex Field with 53 bits of precision
        """
        if S is int or S is float:
            return ToRDF(S)

        from sage.rings.rational_field import QQ
        try:
            from sage.rings.real_lazy import RLF
        except ImportError:
            RLF = None

        if S is ZZ or S is QQ or S is RLF:
            return ToRDF(S)

        if isinstance(S, sage.rings.abc.RealField):
            if S.prec() >= 53:
                return ToRDF(S)
            else:
                return None
        elif is_numpy_type(S):
            import numpy
            if issubclass(S, numpy.integer) or issubclass(S, numpy.floating):
                return ToRDF(S)
            else:
                return None

        try:
            from sage.rings.real_mpfr import RR
        except ImportError:
            pass
        else:
            connecting = RR._internal_coerce_map_from(S)
            if connecting is not None:
                return ToRDF(RR) * connecting

    def _magma_init_(self, magma):
        r"""
        Return a string representation of ``self`` in the Magma language.

        EXAMPLES:

        Magma handles precision in decimal digits, so we lose a bit::

            sage: magma(RDF)        # indirect doctest  # optional - magma
            Real field of precision 15
            sage: 10^15 < 2^53 < 10^16
            True

        When we convert back from Magma, we convert to a generic real field
        that has 53 bits of precision::

            sage: magma(RDF).sage()                     # optional - magma
            Real Field with 53 bits of precision
        """
        return "RealField(%s : Bits := true)" % self.prec()

    def _fricas_init_(self):
        r"""
        Return the FriCAS representation of the real double field.

        EXAMPLES::

            sage: fricas(RDF)       # indirect doctest  # optional - fricas
            DoubleFloat
        """
        return "DoubleFloat"

    def _polymake_init_(self):
        r"""
        Return the polymake representation of the real double field.

        EXAMPLES::

            sage: polymake(RDF)     # indirect doctest  # optional - jupymake
            Float
        """
        return '"Float"'

    def precision(self):
        """
        Return the precision of this real double field in bits.

        Always returns 53.

        EXAMPLES::

            sage: RDF.precision()
            53
        """
        return 53

    prec = precision

    def to_prec(self, prec):
        """
        Return the real field to the specified precision. As doubles have
        fixed precision, this will only return a real double field if ``prec``
        is exactly 53.

        EXAMPLES::

            sage: RDF.to_prec(52)                                                       # needs sage.rings.real_mpfr
            Real Field with 52 bits of precision
            sage: RDF.to_prec(53)
            Real Double Field
        """
        if prec == 53:
            return self
        from sage.rings.real_mpfr import RealField
        return RealField(prec)

    def gen(self, n=0):
        """
        Return the generator of the real double field.

        EXAMPLES::

            sage: RDF.0
            1.0
            sage: RDF.gens()
            (1.0,)
        """
        if n != 0:
            raise ValueError("only 1 generator")
        return RealDoubleElement(1)

    def ngens(self):
        """
        Return the number of generators which is always 1.

        EXAMPLES::

            sage: RDF.ngens()
            1
        """
        return 1

    def characteristic(self):
        """
        Return 0, since the field of real numbers has characteristic 0.

        EXAMPLES::

            sage: RDF.characteristic()
            0
        """
        return Integer(0)

    cdef _new_c(self, double value):
        cdef RealDoubleElement x
        x = PY_NEW(RealDoubleElement)
        x._value = value
        return x

    def random_element(self, double min=-1, double max=1):
        """
        Return a random element of this real double field in the interval
        ``[min, max]``.

        EXAMPLES::

            sage: RDF.random_element().parent() is RDF
            True
            sage: -1 <= RDF.random_element() <= 1
            True
            sage: 100 <= RDF.random_element(min=100, max=110) <= 110
            True
        """
        cdef randstate rstate = current_randstate()

        return self._new_c((max-min)*rstate.c_rand_double() + min)

    def name(self):
        """
        The name of ``self``.

        EXAMPLES::

            sage: RDF.name()
            'RealDoubleField'
        """
        return "RealDoubleField"

    def __hash__(self):
        """
        Return the hash value of ``self``.

        This class is intended for use as a singleton so any instance
        of it should be equivalent from a hashing perspective.

        TESTS::

            sage: from sage.rings.real_double import RealDoubleField_class
            sage: hash(RDF) == hash(RealDoubleField_class())
            True
        """
        return 1157042230

    def pi(self):
        r"""
        Return `\pi` to double-precision.

        EXAMPLES::

            sage: RDF.pi()
            3.141592653589793
            sage: RDF.pi().sqrt()/2
            0.8862269254527579
        """
        return self(M_PI)

    def euler_constant(self):
        """
        Return Euler's gamma constant to double precision.

        EXAMPLES::

            sage: RDF.euler_constant()
            0.5772156649015329
        """
        return self(M_EULER)

    def log2(self):
        r"""
        Return `\log(2)` to the precision of this field.

        EXAMPLES::

            sage: RDF.log2()
            0.6931471805599453
            sage: RDF(2).log()
            0.6931471805599453
        """
        return self(M_LN2)

    def factorial(self, int n):
        """
        Return the factorial of the integer `n` as a real number.

        EXAMPLES::

            sage: RDF.factorial(100)
            9.332621544394415e+157
        """
        return global_dummy_element._factorial(n)

    def zeta(self, n=2):
        """
        Return an `n`-th root of unity in the real field, if one
        exists, or raise a :exc:`ValueError` otherwise.

        EXAMPLES::

            sage: RDF.zeta()
            -1.0
            sage: RDF.zeta(1)
            1.0
            sage: RDF.zeta(5)
            Traceback (most recent call last):
            ...
            ValueError: No 5th root of unity in self
        """
        if n == 1:
            return self(1)
        elif n == 2:
            return self(-1)
        raise ValueError("No %sth root of unity in self" % n)

    def NaN(self):
        """
        Return Not-a-Number ``NaN``.

        EXAMPLES::

            sage: RDF.NaN()
            NaN
        """
        return self(0)/self(0)

    nan = NaN

    def _factor_univariate_polynomial(self, f):
        """
        Factor the univariate polynomial ``f``.

        INPUT:

        - ``f`` -- a univariate polynomial defined over the double precision
          real numbers

        OUTPUT:

        - A factorization of ``f`` over the double precision real numbers
          into a unit and monic irreducible factors

        .. NOTE::

            This is a helper method for
            :meth:`sage.rings.polynomial.polynomial_element.Polynomial.factor`.

        TESTS::

            sage: # needs numpy
            sage: R.<x> = RDF[]
            sage: RDF._factor_univariate_polynomial(x)
            x
            sage: RDF._factor_univariate_polynomial(2*x)
            (2.0) * x
            sage: RDF._factor_univariate_polynomial(x^2)
            x^2
            sage: RDF._factor_univariate_polynomial(x^2 + 1)
            x^2 + 1.0
            sage: RDF._factor_univariate_polynomial(x^2 - 1)
            (x - 1.0) * (x + 1.0)

        The implementation relies on the ``roots()`` method which often reports
        roots not to be real even though they are::

            sage: f = (x-1)^3                                                           # needs numpy
            sage: f.roots(ring=CDF)  # abs tol 2e-5                                     # needs numpy
            [(1.0000065719436413, 1),
             (0.9999967140281792 - 5.691454546815028e-06*I, 1),
             (0.9999967140281792 + 5.691454546815028e-06*I, 1)]

        This leads to the following incorrect factorization::

            sage: f.factor()  # abs tol 2e-5                                            # needs numpy
            (x - 1.0000065719436413) * (x^2 - 1.9999934280563585*x + 0.9999934280995487)
        """
        from sage.rings.complex_double import CDF
        roots = f.roots(CDF)

        # collect real roots and conjugate pairs of non-real roots
        real_roots = [(r, e) for r, e in roots if r.imag().is_zero()]
        non_real_roots = {r: e for r, e in roots if not r.imag().is_zero()}
        assert all(non_real_roots[r.conj()] == e for r, e in non_real_roots.items()), "Bug in root finding code over RDF - roots must always come in conjugate pairs"
        non_real_roots = [(r, e) for r, e in non_real_roots.items() if r.imag() > 0]

        # turn the roots into irreducible factors
        x = f.parent().gen()
        real_factors = [(x - r.real(), e) for r, e in real_roots]
        non_real_factors = [(x**2 - (r + r.conj()).real()*x + (r*r.conj()).real(), e) for r, e in non_real_roots]

        # make the factors monic
        from sage.structure.factorization import Factorization
        return Factorization([(g.monic(), e) for g, e in real_factors + non_real_factors], f.leading_coefficient())


cdef class RealDoubleElement(FieldElement):
    """
    An approximation to a real number using double precision floating
    point numbers. Answers derived from calculations with such
    approximations may differ from what they would be if those
    calculations were performed with true real numbers. This is due to
    the rounding errors inherent to finite precision calculations.
    """

    __array_interface__ = {'typestr': '=f8'}

    def __cinit__(self):
        """
        Initialize ``self`` for cython.

        EXAMPLES::

            sage: RDF(2.3)  # indirect doctest
            2.3
        """
        (<Element>self)._parent = _RDF

    def __init__(self, x):
        """
        Create a new ``RealDoubleElement`` with value ``x``.

        EXAMPLES::

            sage: RDF(10^100)
            1e+100

        TESTS::

            sage: from gmpy2 import *
            sage: RDF(mpz(42))
            42.0
            sage: RDF(mpq(3/4))
            0.75
            sage: RDF(mpq('4.1'))
            4.1
        """
        self._value = float(x)

    def _magma_init_(self, magma):
        r"""
        Return a string representation of ``self`` in the Magma language.

        EXAMPLES::

            sage: RDF(10.5)
            10.5
            sage: magma(RDF(10.5))  # indirect doctest  # optional - magma
            10.5000000000000
        """
        return "%s!%s" % (self.parent()._magma_init_(magma), self)

    def __reduce__(self):
        """
        For pickling.

        EXAMPLES::

            sage: a = RDF(-2.7)
            sage: loads(dumps(a)) == a
            True
        """
        return RealDoubleElement, (self._value, )

    cdef _new_c(self, double value):
        cdef RealDoubleElement x
        x = PY_NEW(RealDoubleElement)
        x._value = value
        return x

    def prec(self):
        """
        Return the precision of this number in bits.

        Always returns 53.

        EXAMPLES::

            sage: RDF(0).prec()
            53
        """
        return 53

    def ulp(self):
        """
        Return the unit of least precision of ``self``, which is the
        weight of the least significant bit of ``self``. This is always
        a strictly positive number. It is also the gap between this
        number and the closest number with larger absolute value that
        can be represented.

        EXAMPLES::

            sage: a = RDF(pi)                                                           # needs sage.symbolic
            sage: a.ulp()                                                               # needs sage.symbolic
            4.440892098500626e-16
            sage: b = a + a.ulp()                                                       # needs sage.symbolic

        Adding or subtracting an ulp always gives a different number::

            sage: # needs sage.symbolic
            sage: a + a.ulp() == a
            False
            sage: a - a.ulp() == a
            False
            sage: b + b.ulp() == b
            False
            sage: b - b.ulp() == b
            False

        Since the default rounding mode is round-to-nearest, adding or
        subtracting something less than half an ulp always gives the
        same number, unless the result has a smaller ulp. The latter
        can only happen if the input number is (up to sign) exactly a
        power of 2::

            sage: # needs sage.symbolic
            sage: a - a.ulp()/3 == a
            True
            sage: a + a.ulp()/3 == a
            True
            sage: b - b.ulp()/3 == b
            True
            sage: b + b.ulp()/3 == b
            True

            sage: c = RDF(1)
            sage: c - c.ulp()/3 == c
            False
            sage: c.ulp()
            2.220446049250313e-16
            sage: (c - c.ulp()).ulp()
            1.1102230246251565e-16

        The ulp is always positive::

            sage: RDF(-1).ulp()
            2.220446049250313e-16

        The ulp of zero is the smallest positive number in RDF::

            sage: RDF(0).ulp()
            5e-324
            sage: RDF(0).ulp()/2
            0.0

        Some special values::

            sage: a = RDF(1)/RDF(0); a
            +infinity
            sage: a.ulp()
            +infinity
            sage: (-a).ulp()
            +infinity
            sage: a = RDF('nan')
            sage: a.ulp() is a
            True

        The ulp method works correctly with small numbers::

            sage: u = RDF(0).ulp()
            sage: u.ulp() == u
            True
            sage: x = u * (2^52-1)  # largest denormal number
            sage: x.ulp() == u
            True
            sage: x = u * 2^52  # smallest normal number
            sage: x.ulp() == u
            True
        """
        # First, check special values
        if self._value == 0:
            return RealDoubleElement(libc.math.ldexp(1.0, -1074))
        if libc.math.isnan(self._value):
            return self
        if libc.math.isinf(self._value):
            return self.abs()

        # Normal case
        cdef int e
        libc.math.frexp(self._value, &e)
        e -= 53
        # Correction for denormals
        if e < -1074:
            e = -1074
        return RealDoubleElement(libc.math.ldexp(1.0, e))

    def real(self):
        """
        Return ``self`` - we are already real.

        EXAMPLES::

            sage: a = RDF(3)
            sage: a.real()
            3.0
        """
        return self

    def imag(self):
        """
        Return the imaginary part of this number, which is zero.

        EXAMPLES::

            sage: a = RDF(3)
            sage: a.imag()
            0.0
        """
        return RealDoubleElement(0)

    def __complex__(self):
        """
        Return ``self`` as a python complex number.

        EXAMPLES::

            sage: a = 2303
            sage: RDF(a)
            2303.0
            sage: complex(RDF(a))
            (2303+0j)
        """
        return complex(self._value, 0)

    def _integer_(self, ZZ=None):
        """
        If this floating-point number is actually an integer, return
        that integer.  Otherwise, raise an exception.

        EXAMPLES::

            sage: ZZ(RDF(237.0))  # indirect doctest
            237
            sage: ZZ(RDF(0.0/0.0))
            Traceback (most recent call last):
            ...
            ValueError: cannot convert float NaN to integer
            sage: ZZ(RDF(1.0/0.0))
            Traceback (most recent call last):
            ...
            OverflowError: cannot convert float infinity to integer
            sage: ZZ(RDF(-123456789.0))
            -123456789
            sage: ZZ(RDF((2.0))^290)
            1989292945639146568621528992587283360401824603189390869761855907572637988050133502132224
            sage: ZZ(RDF(-2345.67))
            Traceback (most recent call last):
            ...
            TypeError: cannot convert non-integral float to integer
        """
        return Integer(self._value)

    def __mpfr__(self):
        """
        Convert Sage ``RealDoubleElement`` to gmpy2 ``mpfr``.

        EXAMPLES::

            sage: RDF(42.2).__mpfr__()
            mpfr('42.200000000000003')
            sage: from gmpy2 import mpfr
            sage: mpfr(RDF(5.1))
            mpfr('5.0999999999999996')

        TESTS::

            sage: RDF().__mpfr__(); raise NotImplementedError("gmpy2 is not installed")
            Traceback (most recent call last):
            ...
            NotImplementedError: gmpy2 is not installed
        """
        return gmpy2.mpfr(self._value)

    def _interface_init_(self, I=None):
        """
        Return ``self`` formatted as a string, suitable as input to another
        computer algebra system. (This is the default function used for
        exporting to other computer algebra systems.)

        EXAMPLES::

            sage: s1 = RDF(sin(1)); s1                                                  # needs sage.symbolic
            0.8414709848078965
            sage: s1._interface_init_()                                                 # needs sage.symbolic
            '0.8414709848078965'
            sage: s1 == RDF(gp(s1))                                                     # needs sage.libs.pari sage.symbolic
            True
        """
        return repr(self._value)

    def _mathematica_init_(self):
        """
        TESTS:

        Check that :issue:`28814` is fixed::

            sage: mathematica(RDF(1e25))   # optional - mathematica
            1.*^25
            sage: mathematica(RDF(1e-25))  # optional - mathematica
            1.*^-25
        """
        from sage.rings.real_mpfr import RR
        return RR(self._value)._mathematica_init_()

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: sage_input(RDF(NaN))                                                  # needs sage.symbolic
            RDF(NaN)
            sage: sage_input(RDF(-infinity), verify=True)
            # Verified
            -RDF(infinity)
            sage: sage_input(RDF(-infinity)*polygen(RDF))
            R.<x> = RDF[]
            -RDF(infinity)*x + RDF(NaN)
            sage: sage_input(RDF(pi), verify=True)                                      # needs sage.symbolic
            # Verified
            RDF(3.1415926535897931)
            sage: sage_input(RDF(-e), verify=True, preparse=False)                      # needs sage.symbolic
            # Verified
            -RDF(2.718281828459045...)
            sage: sage_input(RDF(pi)*polygen(RDF), verify=True, preparse=None)          # needs sage.symbolic
            # Verified
            R = RDF['x']
            x = R.gen()
            3.1415926535897931*x
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: RDF(22/7)._sage_input_(sib, True)                                     # needs sage.sage.rings.real_mpfr
            {atomic:3.1428571428571428}
            sage: RDF(22/7)._sage_input_(sib, False)                                    # needs sage.sage.rings.real_mpfr
            {call: {atomic:RDF}({atomic:3.1428571428571428})}
        """
        cdef bint isinf = libc.math.isinf(self._value)
        cdef bint isnan = libc.math.isnan(self._value)
        if isinf or isnan:
            if isnan:
                v = sib.name('NaN')
            else:
                v = sib.name('infinity')
            v = sib(self.parent())(v)
            if self._value < 0:
                v = -v
            return v

        from sage.rings.integer_ring import ZZ
        from sage.rings.real_mpfr import RR

        cdef bint negative = self._value < 0
        if negative:
            self = -self

        # There are five possibilities for printing this floating-point
        # number, ordered from prettiest to ugliest (IMHO).
        # 1) An integer: 42
        # 2) A simple literal: 3.14159
        # 3) A coerced integer: RDF(42)
        # 4) A coerced literal: RDF(3.14159)
        # 5) A coerced RR value: RDF(RR('3.14159'))

        # str(self) works via libc, which we don't necessarily trust
        # to produce the best possible answer.  So this function prints
        # via RR/MPFR.  Without the preparser, input works via libc as
        # well, but we don't have a choice about that.

        # To use choice 1 or choice 3, this number must be an integer.
        cdef bint can_use_int_literal = \
            self.abs() < (Integer(1) << self.prec()) and self in ZZ

        self_str = RR(self._value).str(truncate=False, skip_zeroes=True)

        # To use choice 2 or choice 4, we must be able to read
        # numbers of this precision as a literal.
        cdef bint can_use_float_literal = \
            (sib.preparse() or float(self_str) == self)

        if can_use_int_literal or can_use_float_literal:
            if can_use_int_literal:
                v = sib.int(self._integer_())
            else:
                v = sib.float_str(self_str)
        else:
            v = sib(RR(self))
        if not coerced:
            v = sib(self.parent())(v)

        if negative:
            v = -v

        return v

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: RDF(-2/3)
            -0.6666666666666666
            sage: a = RDF(2); a
            2.0
            sage: a^2
            4.0
            sage: RDF("nan")
            NaN
            sage: RDF(1)/RDF(0)
            +infinity
            sage: RDF(-1)/RDF(0)
            -infinity
        """
        return double_repr(self._value)

    def __format__(self, format_spec):
        """
        Return a formatted string representation of this real number.

        EXAMPLES::

            sage: format(RDF(32/3), '.4f')
            '10.6667'
            sage: '{:.4e}'.format(RDF(2/3))
            '6.6667e-01'
        """
        return format(float(self), format_spec)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(RDF(3.4))       # indirect doctest
            3.4
            sage: latex(RDF(2e-100))    # indirect doctest
            2 \times 10^{-100}
        """
        s = self.str()
        parts = s.split('e')
        if len(parts) > 1:
            # scientific notation
            if parts[1][0] == '+':
                parts[1] = parts[1][1:]
            s = "%s \\times 10^{%s}" % (parts[0], parts[1])
        return s

    def __hash__(self):
        """
        Return the hash of ``self``, which coincides with the python float
        (and often int) type.

        EXAMPLES::

            sage: hash(RDF(1.2)) == hash(1.2r)
            True
        """
        return hash(self._value)

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        Return the image of ``self`` under the homomorphism from the rational
        field to ``codomain``.

        This always just returns ``self`` coerced into the ``codomain``.

        EXAMPLES::

            sage: RDF(2.1)._im_gens_(RR, [RR(1)])
            2.10000000000000
            sage: R = RealField(20)                                                     # needs sage.rings.real_mpfr
            sage: RDF(2.1)._im_gens_(R, [R(1)])                                         # needs sage.rings.real_mpfr
            2.1000
        """
        return codomain(self)  # since 1 |--> 1

    def str(self):
        """
        Return the informal string representation of ``self``.

        EXAMPLES::

            sage: a = RDF('4.5'); a.str()
            '4.5'
            sage: a = RDF('49203480923840.2923904823048'); a.str()
            '49203480923840.29'
            sage: a = RDF(1)/RDF(0); a.str()
            '+infinity'
            sage: a = -RDF(1)/RDF(0); a.str()
            '-infinity'
            sage: a = RDF(0)/RDF(0); a.str()
            'NaN'

        We verify consistency with ``RR`` (mpfr reals)::

            sage: str(RR(RDF(1)/RDF(0))) == str(RDF(1)/RDF(0))
            True
            sage: str(RR(-RDF(1)/RDF(0))) == str(-RDF(1)/RDF(0))
            True
            sage: str(RR(RDF(0)/RDF(0))) == str(RDF(0)/RDF(0))
            True
        """
        return double_repr(self._value)

    def __copy__(self):
        """
        Return copy of ``self``, which since ``self`` is immutable, is just
        ``self``.

        EXAMPLES::

            sage: r = RDF('-1.6')
            sage: r.__copy__() is r
            True
        """
        return self

    def __deepcopy__(self, memo):
        """
        EXAMPLES::

            sage: r = RDF('-1.6')
            sage: deepcopy(r) is r
            True
        """
        return self

    def integer_part(self):
        """
        If in decimal this number is written ``n.defg``, returns ``n``.

        EXAMPLES::

            sage: r = RDF('-1.6')
            sage: a = r.integer_part(); a
            -1
            sage: type(a)
            <class 'sage.rings.integer.Integer'>
            sage: r = RDF(0.0/0.0)
            sage: a = r.integer_part()
            Traceback (most recent call last):
            ...
            TypeError: Attempt to get integer part of NaN
        """
        if libc.math.isnan(self._value):
            raise TypeError("Attempt to get integer part of NaN")
        else:
            return Integer(int(self._value))

    def sign_mantissa_exponent(self):
        r"""
        Return the sign, mantissa, and exponent of ``self``.

        In Sage (as in MPFR), floating-point numbers of precision `p`
        are of the form `s m 2^{e-p}`, where `s \in \{-1, 1\}`,
        `2^{p-1} \leq m < 2^p`, and `-2^{30} + 1 \leq e \leq 2^{30} -
        1`; plus the special values ``+0``, ``-0``, ``+infinity``,
        ``-infinity``, and ``NaN`` (which stands for Not-a-Number).

        This function returns `s`, `m`, and `e-p`.  For the special values:

        - ``+0`` returns ``(1, 0, 0)``
        - ``-0`` returns ``(-1, 0, 0)``
        - the return values for ``+infinity``, ``-infinity``, and ``NaN`` are
          not specified.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: a = RDF(exp(1.0)); a
            2.718281828459045
            sage: sign, mantissa, exponent = RDF(exp(1.0)).sign_mantissa_exponent()
            sage: sign, mantissa, exponent
            (1, 6121026514868073, -51)
            sage: sign*mantissa*(2**exponent) == a
            True

        The mantissa is always a nonnegative number::

            sage: RDF(-1).sign_mantissa_exponent()                                      # needs sage.rings.real_mpfr
            (-1, 4503599627370496, -52)

        TESTS::

            sage: RDF('+0').sign_mantissa_exponent()                                    # needs sage.rings.real_mpfr
            (1, 0, 0)
            sage: RDF('-0').sign_mantissa_exponent()                                    # needs sage.rings.real_mpfr
            (-1, 0, 0)
        """
        from sage.rings.real_mpfr import RR
        return RR(self._value).sign_mantissa_exponent()

    def as_integer_ratio(self):
        """
        Return a coprime pair of integers ``(a, b)`` such that ``self``
        equals ``a / b`` exactly.

        EXAMPLES::

            sage: RDF(0).as_integer_ratio()
            (0, 1)
            sage: RDF(1/3).as_integer_ratio()
            (6004799503160661, 18014398509481984)
            sage: RDF(37/16).as_integer_ratio()
            (37, 16)
            sage: RDF(3^60).as_integer_ratio()
            (42391158275216203520420085760, 1)
        """
        nd = float.as_integer_ratio(self._value)
        return (Integer(nd[0]), Integer(nd[1]))

    ########################
    #   Basic Arithmetic
    ########################
    def __invert__(self):
        """
        Compute the multiplicative inverse of ``self``.

        EXAMPLES::

            sage: a = RDF(-1.5)*RDF(2.5)
            sage: a.__invert__()
            -0.26666666666666666
            sage: ~a
            -0.26666666666666666
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = 1.0 / self._value
        return x

    cpdef _add_(self, right):
        """
        Add two real numbers with the same parent.

        EXAMPLES::

            sage: RDF('-1.5') + RDF('2.5')  # indirect doctest
            1.0
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = self._value + (<RealDoubleElement>right)._value
        return x

    cpdef _sub_(self, right):
        """
        Subtract two real numbers with the same parent.

        EXAMPLES::

            sage: RDF('-1.5') - RDF('2.5')  # indirect doctest
            -4.0
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = self._value - (<RealDoubleElement>right)._value
        return x

    cpdef _mul_(self, right):
        """
        Multiply two real numbers with the same parent.

        EXAMPLES::

            sage: RDF('-1.5') * RDF('2.5')  # indirect doctest
            -3.75
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = self._value * (<RealDoubleElement>right)._value
        return x

    cpdef _div_(self, right):
        """
        Divide ``self`` by ``right``.

        EXAMPLES::

            sage: RDF('-1.5') / RDF('2.5')  # indirect doctest
            -0.6
            sage: RDF(1)/RDF(0)
            +infinity
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = self._value / (<RealDoubleElement>right)._value
        return x

    def __neg__(self):
        """
        Negate ``self``.

        EXAMPLES::

            sage: -RDF('-1.5')
            1.5
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = -self._value
        return x

    def conjugate(self):
        r"""
        Return the complex conjugate of this real number, which is
        the real number itself.

        EXAMPLES::

            sage: RDF(4).conjugate()
            4.0
        """
        return self

    def __abs__(self):
        """
        Return the absolute value of ``self``.

        EXAMPLES::

            sage: abs(RDF(1.5))
            1.5
            sage: abs(RDF(-1.5))
            1.5
            sage: abs(RDF(0.0))
            0.0
            sage: abs(RDF(-0.0))
            0.0
        """
        # Use signbit instead of >= to handle -0.0 correctly
        if not libc.math.signbit(self._value):
            return self
        else:
            return self._new_c(-self._value)

    cpdef RealDoubleElement abs(RealDoubleElement self):
        """
        Return the absolute value of ``self``.

        EXAMPLES::

            sage: RDF(1e10).abs()
            10000000000.0
            sage: RDF(-1e10).abs()
            10000000000.0
        """
        if self._value >= 0:
            return self
        else:
            return self._new_c(-self._value)

    def __lshift__(x, y):
        """
        LShifting a double is not supported; nor is lshifting a
        :class:`RealDoubleElement`.

        TESTS::

            sage: RDF(2) << 3
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for <<
        """
        raise TypeError("unsupported operand type(s) for <<")

    def __rshift__(x, y):
        """
        RShifting a double is not supported; nor is rshifting a
        :class:`RealDoubleElement`.

        TESTS::

            sage: RDF(2) >> 3
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for >>
        """
        raise TypeError("unsupported operand type(s) for >>")

    def multiplicative_order(self):
        r"""
        Return `n` such that ``self^n == 1``.

        Only `\pm 1` have finite multiplicative order.

        EXAMPLES::

            sage: RDF(1).multiplicative_order()
            1
            sage: RDF(-1).multiplicative_order()
            2
            sage: RDF(3).multiplicative_order()
            +Infinity
        """
        if self._value == 1:
            return 1
        elif self._value == -1:
            return 2
        return sage.rings.infinity.infinity

    def sign(self):
        """
        Return -1, 0, or 1 if ``self`` is negative, zero, or positive;
        respectively.

        EXAMPLES::

            sage: RDF(-1.5).sign()
            -1
            sage: RDF(0).sign()
            0
            sage: RDF(2.5).sign()
            1
        """
        if not self._value:
            return 0
        if self._value > 0:
            return 1
        return -1

    ###################
    # Rounding etc
    ###################

    def round(self):
        """
        Round ``self`` to the nearest integer.

        This uses the convention of rounding half to even
        (i.e., if the fractional part of ``self`` is `0.5`, then it
        is rounded to the nearest even integer).

        EXAMPLES::

            sage: RDF(0.49).round()
            0
            sage: a=RDF(0.51).round(); a
            1
            sage: RDF(0.5).round()
            0
            sage: RDF(1.5).round()
            2
        """
        return Integer(round(self._value))

    def floor(self):
        """
        Return the floor of ``self``.

        EXAMPLES::

            sage: RDF(2.99).floor()
            2
            sage: RDF(2.00).floor()
            2
            sage: RDF(-5/2).floor()
            -3
        """
        return Integer(math.floor(self._value))

    def ceil(self):
        """
        Return the ceiling of ``self``.

        EXAMPLES::

            sage: RDF(2.99).ceil()
            3
            sage: RDF(2.00).ceil()
            2
            sage: RDF(-5/2).ceil()
            -2
        """
        return Integer(math.ceil(self._value))

    ceiling = ceil

    def trunc(self):
        """
        Truncates this number (returns integer part).

        EXAMPLES::

            sage: RDF(2.99).trunc()
            2
            sage: RDF(-2.00).trunc()
            -2
            sage: RDF(0.00).trunc()
            0
        """
        return Integer(int(self._value))

    def frac(self):
        """
        Return a real number in `(-1, 1)`. It satisfies the relation:
        ``x = x.trunc() + x.frac()``

        EXAMPLES::

            sage: RDF(2.99).frac()
            0.9900000000000002
            sage: RDF(2.50).frac()
            0.5
            sage: RDF(-2.79).frac()
            -0.79
        """
        return self._new_c(self._value - int(self._value))

    ###########################################
    # Conversions
    ###########################################

    def __float__(self):
        """
        Return ``self`` as a python float.

        EXAMPLES::

            sage: float(RDF(1.5))
            1.5
            sage: type(float(RDF(1.5)))
            <... 'float'>
        """
        return self._value

    def _rpy_(self):
        """
        Return ``self.__float__()`` for rpy to convert into the
        appropriate R object.

        EXAMPLES::

            sage: n = RDF(2.0)
            sage: n._rpy_()
            2.0
            sage: type(n._rpy_())
            <... 'float'>
        """
        return self.__float__()

    def __int__(self):
        """
        Return integer truncation of this real number.

        EXAMPLES::

            sage: int(RDF(2.99))
            2
            sage: int(RDF(-2.99))
            -2
        """
        return int(self._value)

    def _complex_mpfr_field_(self, CC):
        """
        EXAMPLES::

            sage: a = RDF(1/3)
            sage: CC(a)                                                                 # needs sage.rings.real_mpfr
            0.333333333333333
            sage: a._complex_mpfr_field_(CC)                                            # needs sage.rings.real_mpfr
            0.333333333333333

        If we coerce to a higher-precision field the extra bits appear
        random; they are actually 0s in base 2.

        ::

            sage: a._complex_mpfr_field_(ComplexField(100))                             # needs sage.rings.real_mpfr
            0.33333333333333331482961625625
            sage: a._complex_mpfr_field_(ComplexField(100)).str(2)                      # needs sage.rings.real_mpfr
            '0.01010101010101010101010101010101010101010101010101010100000000000000000000000000000000000000000000000'
        """
        return CC(self._value)

    def _complex_double_(self, CDF):
        """
        Return ``self`` as a complex double.

        EXAMPLES::

            sage: CDF(RDF(1/3))  # indirect doctest                                     # needs sage.rings.complex_double
            0.3333333333333333
        """
        return CDF(self._value)

    def __pari__(self):
        """
        Return a PARI representation of ``self``.

        EXAMPLES::

            sage: RDF(1.5).__pari__()                                                   # needs sage.libs.pari
            1.50000000000000
        """
        global new_gen_from_real_double_element
        if new_gen_from_real_double_element is None:
            from sage.libs.pari.convert_sage_real_double import new_gen_from_real_double_element
        return new_gen_from_real_double_element(self)

    ###########################################
    # Comparisons: ==, !=, <, <=, >, >=
    ###########################################

    def is_NaN(self):
        """
        Check if ``self`` is ``NaN``.

        EXAMPLES::

            sage: RDF(1).is_NaN()
            False
            sage: a = RDF(0)/RDF(0)
            sage: a.is_NaN()
            True
        """
        return bool(libc.math.isnan(self._value))

    def is_positive_infinity(self):
        r"""
        Check if ``self`` is `+\infty`.

        EXAMPLES::

            sage: a = RDF(1)/RDF(0)
            sage: a.is_positive_infinity()
            True
            sage: a = RDF(-1)/RDF(0)
            sage: a.is_positive_infinity()
            False
        """
        if not libc.math.isinf(self._value):
            return False
        return self._value > 0

    def is_negative_infinity(self):
        r"""
        Check if ``self`` is `-\infty`.

        EXAMPLES::

            sage: a = RDF(2)/RDF(0)
            sage: a.is_negative_infinity()
            False
            sage: a = RDF(-3)/RDF(0)
            sage: a.is_negative_infinity()
            True
        """
        if not libc.math.isinf(self._value):
            return False
        return self._value < 0

    def is_infinity(self):
        r"""
        Check if ``self`` is `\infty`.

        EXAMPLES::

            sage: a = RDF(2); b = RDF(0)
            sage: (a/b).is_infinity()
            True
            sage: (b/a).is_infinity()
            False
        """
        return bool(libc.math.isinf(self._value))

    cpdef _richcmp_(left, right, int op):
        """
        Rich comparison of ``left`` and ``right``.

        EXAMPLES::

            sage: RDF(2) < RDF(0)
            False
            sage: RDF(2) == RDF(4/2)
            True
            sage: RDF(-2) > RDF(-4)
            True

        TESTS:

        Check comparisons with ``NaN`` (:issue:`16515`)::

            sage: n = RDF('NaN')
            sage: n == n
            False
            sage: n == RDF(1)
            False
        """
        # We really need to use the correct operators, to deal
        # correctly with NaNs.
        cdef double x = (<RealDoubleElement>left)._value
        cdef double y = (<RealDoubleElement>right)._value
        if op == Py_LT:
            return x < y
        elif op == Py_LE:
            return x <= y
        elif op == Py_EQ:
            return x == y
        elif op == Py_NE:
            return x != y
        elif op == Py_GT:
            return x > y
        else:
            return x >= y

    ############################
    # Special Functions
    ############################

    def NaN(self):
        """
        Return Not-a-Number ``NaN``.

        EXAMPLES::

            sage: RDF.NaN()
            NaN
        """
        return self(0)/self(0)

    nan = NaN

    def sqrt(self, extend=True, all=False):
        """
        The square root function.

        INPUT:

        - ``extend`` -- boolean (default: ``True``); if ``True``, return a
          square root in a complex field if necessary if ``self`` is negative.
          Otherwise raise a :exc:`ValueError`.

        - ``all`` -- boolean (default: ``False``); if ``True``, return a
          list of all square roots

        EXAMPLES::

            sage: r = RDF(4.0)
            sage: r.sqrt()
            2.0
            sage: r.sqrt()^2 == r
            True

        ::

            sage: r = RDF(4344)
            sage: r.sqrt()
            65.90902821313632
            sage: r.sqrt()^2 - r
            0.0

        ::

            sage: r = RDF(-2.0)
            sage: r.sqrt()                                                              # needs sage.rings.complex_double
            1.4142135623730951*I

        ::

            sage: RDF(2).sqrt(all=True)
            [1.4142135623730951, -1.4142135623730951]
            sage: RDF(0).sqrt(all=True)
            [0.0]
            sage: RDF(-2).sqrt(all=True)                                                # needs sage.rings.complex_double
            [1.4142135623730951*I, -1.4142135623730951*I]
        """
        if self._value >= 0:
            x = self._new_c(libc.math.sqrt(self._value))
            if all:
                if x.is_zero():
                    return [x]
                else:
                    return [x, -x]
            else:
                return x
        if not extend:
            raise ValueError("negative number %s does not have a square root in the real field" % self)
        import sage.rings.complex_double
        return self._complex_double_(sage.rings.complex_double.CDF).sqrt(all=all)

    def is_square(self):
        """
        Return whether or not this number is a square in this field. For
        the real numbers, this is ``True`` if and only if ``self`` is
        nonnegative.

        EXAMPLES::

            sage: RDF(3.5).is_square()
            True
            sage: RDF(0).is_square()
            True
            sage: RDF(-4).is_square()
            False
        """
        return self._value >= 0

    def is_integer(self):
        """
        Return ``True`` if this number is a integer.

        EXAMPLES::

            sage: RDF(3.5).is_integer()
            False
            sage: RDF(3).is_integer()
            True
        """
        return self._value in ZZ

    def cube_root(self):
        """
        Return the cubic root (defined over the real numbers) of ``self``.

        EXAMPLES::

            sage: r = RDF(125.0); r.cube_root()
            5.000000000000001
            sage: r = RDF(-119.0)
            sage: r.cube_root()^3 - r  # rel tol 1
            -1.4210854715202004e-14
        """
        return self.nth_root(3)

    def agm(self, other):
        r"""
        Return the arithmetic-geometric mean of ``self`` and ``other``. The
        arithmetic-geometric mean is the common limit of the sequences
        `u_n` and `v_n`, where `u_0` is ``self``,
        `v_0` is other, `u_{n+1}` is the arithmetic mean
        of `u_n` and `v_n`, and `v_{n+1}` is the
        geometric mean of `u_n` and `v_n`. If any operand is negative, the
        return value is ``NaN``.

        EXAMPLES::

            sage: a = RDF(1.5)
            sage: b = RDF(2.3)
            sage: a.agm(b)
            1.8786484558146697

        The arithmetic-geometric mean always lies between the geometric and
        arithmetic mean::

            sage: sqrt(a*b) < a.agm(b) < (a+b)/2
            True
        """
        cdef double a = self._value
        cdef double b = other
        cdef double eps = 2.0**-51
        if a < 0 or b < 0:
            return self._parent.nan()
        while True:
            a1 = (a+b)/2
            b1 = libc.math.sqrt(a*b)
            if abs((b1/a1)-1) < eps:
                return self._new_c(a1)
            a, b = a1, b1

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

            sage: r = sqrt(RDF(2)); r
            1.4142135623730951
            sage: r.algebraic_dependency(5)                                             # needs sage.libs.pari
            x^2 - 2
        """
        return sage.arith.misc.algdep(self, n)

    algdep = algebraic_dependency

cdef class ToRDF(Morphism):
    def __init__(self, R):
        """
        Fast morphism from anything with a ``__float__`` method to an ``RDF``
        element.

        EXAMPLES::

            sage: f = RDF.coerce_map_from(ZZ); f
            Native morphism:
              From: Integer Ring
              To:   Real Double Field
            sage: f(4)
            4.0
            sage: f = RDF.coerce_map_from(QQ); f
            Native morphism:
              From: Rational Field
              To:   Real Double Field
            sage: f(1/2)
            0.5
            sage: f = RDF.coerce_map_from(int); f
            Native morphism:
              From: Set of Python objects of class 'int'
              To:   Real Double Field
            sage: f(3r)
            3.0
            sage: f = RDF.coerce_map_from(float); f
            Native morphism:
              From: Set of Python objects of class 'float'
              To:   Real Double Field
            sage: f(3.5)
            3.5
        """
        from sage.categories.homset import Hom
        if isinstance(R, type):
            from sage.sets.pythonclass import Set_PythonType
            R = Set_PythonType(R)
        Morphism.__init__(self, Hom(R, RDF))

    cpdef Element _call_(self, x):
        """
        Send ``x`` to the image under this map.

        EXAMPLES::

            sage: f = RDF.coerce_map_from(float)
            sage: f(3.5)  # indirect doctest
            3.5
        """
        cdef RealDoubleElement r = <RealDoubleElement>PY_NEW(RealDoubleElement)
        r._value = PyFloat_AsDouble(x)
        return r

    def _repr_type(self):
        """
        Return the representation type of ``self``.

        EXAMPLES::

            sage: RDF.coerce_map_from(float)._repr_type()
            'Native'
        """
        return "Native"


#####################################################
# unique objects
#####################################################
cdef RealDoubleField_class _RDF
_RDF = RealDoubleField_class()

RDF = _RDF   # external interface


def RealDoubleField():
    """
    Return the unique instance of the
    :class:`real double field<RealDoubleField_class>`.

    EXAMPLES::

        sage: RealDoubleField() is RealDoubleField()
        True
    """
    global _RDF
    return _RDF


def is_RealDoubleElement(x):
    """
    Check if ``x`` is an element of the real double field.

    EXAMPLES::

        sage: from sage.rings.real_double import is_RealDoubleElement
        sage: is_RealDoubleElement(RDF(3))
        doctest:warning...
        DeprecationWarning: The function is_RealDoubleElement is deprecated;
        use 'isinstance(..., RealDoubleElement)' instead.
        See https://github.com/sagemath/sage/issues/38128 for details.
        True
        sage: is_RealDoubleElement(RIF(3))                                              # needs sage.rings.real_interval_field
        False
    """
    from sage.misc.superseded import deprecation_cython
    deprecation_cython(38128,
                       "The function is_RealDoubleElement is deprecated; "
                       "use 'isinstance(..., RealDoubleElement)' instead.")
    return isinstance(x, RealDoubleElement)


# ################ FAST CREATION CODE ######################
#            Based on fast integer creation code
#         There is nothing to see here, move along

# We use a global element to steal all the references
# from.  DO NOT INITIALIZE IT AGAIN and DO NOT REFERENCE IT!
cdef RealDoubleElement global_dummy_element

try:
    from sage.rings.real_double_element_gsl import RealDoubleElement_gsl
except ImportError:
    global_dummy_element = RealDoubleElement(0)
else:
    global_dummy_element = RealDoubleElement_gsl(0)

# A global pool for performance when elements are rapidly created and destroyed.
# It operates on the following principles:
#
# - The pool starts out empty.
# - When a new element is needed, one from the pool is returned
#   if available, otherwise a new RealDoubleElement object is created
# - When an element is collected, it will add it to the pool
#   if there is room, otherwise it will be deallocated.
DEF element_pool_size = 50

cdef PyObject* element_pool[element_pool_size]
cdef int element_pool_count = 0

# used for profiling the pool
cdef int total_alloc = 0
cdef int use_pool = 0


cdef PyObject* fast_tp_new(type t, args, kwds) noexcept:
    global element_pool, element_pool_count, total_alloc, use_pool

    cdef PyObject* new

    # for profiling pool usage
    # total_alloc += 1

    # If there is a ready real double in the pool, we will
    # decrement the counter and return that.

    if element_pool_count > 0:

        # for profiling pool usage
        # use_pool += 1

        element_pool_count -= 1
        new = <PyObject *> element_pool[element_pool_count]

    # Otherwise, we have to create one.

    else:

        # allocate enough room for the RealDoubleElement,
        # sizeof_RealDoubleElement is sizeof(RealDoubleElement).
        # The use of PyObject_Malloc directly assumes
        # that RealDoubleElements are not garbage collected, i.e.
        # they do not possess references to other Python
        # objects (As indicated by the Py_TPFLAGS_HAVE_GC flag).
        # See below for a more detailed description.

        new = <PyObject*>PyObject_Malloc(sizeof(RealDoubleElement))

        # Now set every member as set in z, the global dummy RealDoubleElement
        # created before this tp_new started to operate.

        memcpy(new, (<void*>global_dummy_element), sizeof(RealDoubleElement))

    # This line is only needed if Python is compiled in debugging mode
    # './configure --with-pydebug' or SAGE_DEBUG=yes. If that is the
    # case a Python object has a bunch of debugging fields which are
    # initialized with this macro.
    if_Py_TRACE_REFS_then_PyObject_INIT(
            new, Py_TYPE(global_dummy_element))

    # The global_dummy_element may have a reference count larger than
    # one, but it is expected that newly created objects have a
    # reference count of one. This is potentially unneeded if
    # everybody plays nice, because the global_dummy_element has only
    # one reference in that case.

    # Objects from the pool have reference count zero, so this
    # needs to be set in this case.

    new.ob_refcnt = 1

    return new

cdef void fast_tp_dealloc(PyObject* o) noexcept:

    # If there is room in the pool for a used integer object,
    # then put it in rather than deallocating it.

    global element_pool, element_pool_count

    if element_pool_count < element_pool_size:

        # And add it to the pool.
        element_pool[element_pool_count] = o
        element_pool_count += 1
        return

    # Free the object. This assumes that Py_TPFLAGS_HAVE_GC is not
    # set. If it was set another free function would need to be
    # called.

    PyObject_Free(o)


from sage.misc.allocator cimport hook_tp_functions, hook_tp_functions_type
hook_tp_functions(global_dummy_element, <newfunc>(&fast_tp_new), <destructor>(&fast_tp_dealloc), False)
try:
    from sage.rings.real_double_element_gsl import RealDoubleElement_gsl
except Exception:
    pass
else:
    # global_dummy_element is of type RealDoubleElement_gsl,
    # so hook the base class now.
    hook_tp_functions_type(RealDoubleElement, <newfunc>(&fast_tp_new), <destructor>(&fast_tp_dealloc), False)
    # From here on, calling PY_NEW(RealDoubleElement) actually creates an instance of RealDoubleElement_gsl


cdef double_repr(double x):
    """
    Convert a double to a string with maximum precision.
    """
    if libc.math.isfinite(x):
        return repr(x)
    if libc.math.isinf(x):
        if x > 0:
            return "+infinity"
        if x < 0:
            return "-infinity"
    return "NaN"


# Support Python's numbers abstract base class
import numbers
numbers.Real.register(RealDoubleElement)
