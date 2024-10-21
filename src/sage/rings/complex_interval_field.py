# sage_setup: distribution = sagemath-flint
r"""
Arbitrary precision complex intervals (parent class)

AUTHORS:

- William Stein wrote complex_field.py.

- William Stein (2006-01-26): complete rewrite

Then ``complex_field.py`` was copied to ``complex_interval_field.py`` and
heavily modified:

- Carl Witty (2007-10-24): rewrite for intervals

- Niles Johnson (2010-08): :issue:`3893`: ``random_element()``
  should pass on ``*args`` and ``**kwds``.

- Travis Scrimshaw (2012-10-18): Added documentation to get full coverage.

.. NOTE::

    The :class:`ComplexIntervalField` differs from :class:`ComplexField` in
    that :class:`ComplexIntervalField` only gives the digits with exact
    precision, then a ``?`` signifying that the last digit can have an error of
    ``+/-1``.
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


import weakref

import sage.rings.abc
from sage.misc.cachefunc import cached_method
from sage.rings import complex_interval, integer
from sage.rings.complex_mpfr import ComplexField
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_mpfi import RealIntervalField, RealIntervalField_class
from sage.rings.ring import Field
from sage.structure.parent import Parent

cache = {}
def ComplexIntervalField(prec=53, names=None):
    """
    Return the complex interval field with real and imaginary parts having
    ``prec`` *bits* of precision.

    EXAMPLES::

        sage: ComplexIntervalField()
        Complex Interval Field with 53 bits of precision
        sage: ComplexIntervalField(100)
        Complex Interval Field with 100 bits of precision
        sage: ComplexIntervalField(100).base_ring()
        Real Interval Field with 100 bits of precision
        sage: i = ComplexIntervalField(200).gen()
        sage: i^2
        -1
        sage: i^i
        0.207879576350761908546955619834978770033877841631769608075136?
    """
    global cache
    if prec in cache:
        X = cache[prec]
        C = X()
        if C is not None:
            return C
    C = ComplexIntervalField_class(prec)
    cache[prec] = weakref.ref(C)
    return C


class ComplexIntervalField_class(sage.rings.abc.ComplexIntervalField):
    """
    The field of complex (interval) numbers.

    EXAMPLES::

        sage: C = ComplexIntervalField(); C
        Complex Interval Field with 53 bits of precision
        sage: Q = RationalField()
        sage: C(1/3)
        0.3333333333333334?
        sage: C(1/3, 2)
        0.3333333333333334? + 2*I

    We can also coerce rational numbers and integers into ``C``, but
    coercing a polynomial will raise an exception::

        sage: Q = RationalField()
        sage: C(1/3)
        0.3333333333333334?
        sage: S.<x> = PolynomialRing(Q)
        sage: C(x)
        Traceback (most recent call last):
        ...
        TypeError: cannot convert nonconstant polynomial

    This illustrates precision::

        sage: CIF = ComplexIntervalField(10); CIF(1/3, 2/3)
        0.334? + 0.667?*I
        sage: CIF
        Complex Interval Field with 10 bits of precision
        sage: CIF = ComplexIntervalField(100); CIF
        Complex Interval Field with 100 bits of precision
        sage: z = CIF(1/3, 2/3); z
        0.333333333333333333333333333334? + 0.666666666666666666666666666667?*I

    We can load and save complex numbers and the complex interval field::

        sage: saved_z = loads(z.dumps())
        sage: saved_z.endpoints() == z.endpoints()
        True
        sage: loads(CIF.dumps()) == CIF
        True
        sage: k = ComplexIntervalField(100)
        sage: loads(dumps(k)) == k
        True

    This illustrates basic properties of a complex (interval) field::

        sage: CIF = ComplexIntervalField(200)
        sage: CIF.is_field()
        True
        sage: CIF.characteristic()
        0
        sage: CIF.precision()
        200
        sage: CIF.variable_name()
        'I'
        sage: CIF == ComplexIntervalField(200)
        True
        sage: CIF == ComplexIntervalField(53)
        False
        sage: CIF == 1.1
        False
        sage: CIF = ComplexIntervalField(53)

        sage: CIF.category()
        Category of infinite fields
        sage: TestSuite(CIF).run(skip='_test_gcd_vs_xgcd')

    TESTS:

    This checks that :issue:`15355` is fixed::

        sage: x + CIF(RIF(-2,2), 0)
        x + 0.?e1
        sage: x + CIF(RIF(-2,2), RIF(-2,2))
        x + 0.?e1 + 0.?e1*I
        sage: x + RIF(-2,2)
        x + 0.?e1
        sage: x + CIF(RIF(3.14,3.15), RIF(3.14, 3.15))
        x + 3.15? + 3.15?*I
        sage: CIF(RIF(-2,2), RIF(-2,2))
        0.?e1 + 0.?e1*I
        sage: x + CIF(RIF(3.14,3.15), 0)
        x + 3.15?

    Methods inherited from categories::

        sage: CIF.is_finite()
        False

    .. SEEALSO::

        - :mod:`sage.rings.real_mpfi`
        - :class:`sage.rings.complex_arb.ComplexBallField` (alternative
          implementation of complex intervals, with more features)
    """
    Element = complex_interval.ComplexIntervalFieldElement

    def __init__(self, prec=53):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: ComplexIntervalField()
            Complex Interval Field with 53 bits of precision
            sage: ComplexIntervalField(200)
            Complex Interval Field with 200 bits of precision
        """
        self._prec = int(prec)
        from sage.categories.fields import Fields
        Field.__init__(self, self.real_field(), ("I",), False,
                category=Fields().Infinite())
        self._populate_coercion_lists_(convert_method_name='_complex_mpfi_')

    def __reduce__(self):
        """
        Used for pickling.

        TESTS::

            sage: loads(dumps(CIF)) == CIF
            True
        """
        return ComplexIntervalField, (self._prec, )

    def construction(self):
        """
        Return the functorial construction of this complex interval field,
        namely as the algebraic closure of the real interval field with
        the same precision.

        EXAMPLES::

            sage: c, S = CIF.construction(); c, S
            (AlgebraicClosureFunctor,
             Real Interval Field with 53 bits of precision)
            sage: CIF == c(S)
            True

        TESTS:

        Test that :issue:`19922` is fixed::

            sage: c = ComplexIntervalField(128).an_element()
            sage: r = RealIntervalField(64).an_element()
            sage: c + r
            1 + 1*I
            sage: r + c
            1 + 1*I
            sage: parent(c+r)
            Complex Interval Field with 64 bits of precision
            sage: R = ComplexIntervalField(128)['x']
            sage: (R.gen() * RIF.one()).parent()
            Univariate Polynomial Ring in x over Complex Interval Field with 53 bits of precision
        """
        from sage.categories.pushout import AlgebraicClosureFunctor
        return (AlgebraicClosureFunctor(), self.real_field())

    def is_exact(self):
        """
        The complex interval field is not exact.

        EXAMPLES::

            sage: CIF.is_exact()
            False
        """
        return False

    def prec(self):
        """
        Return the precision of ``self`` (in bits).

        EXAMPLES::

            sage: CIF.prec()
            53
            sage: ComplexIntervalField(200).prec()
            200
        """
        return self._prec

    def to_prec(self, prec):
        """
        Return a complex interval field with the given precision.

        EXAMPLES::

            sage: CIF.to_prec(150)
            Complex Interval Field with 150 bits of precision
            sage: CIF.to_prec(15)
            Complex Interval Field with 15 bits of precision
            sage: CIF.to_prec(53) is CIF
            True
        """
        return ComplexIntervalField(prec)

    def _magma_init_(self, magma):
        r"""
        Return a string representation of ``self`` in the Magma language.

        EXAMPLES::

            sage: magma(ComplexIntervalField(100)) # optional - magma # indirect doctest
            Complex field of precision 30
            sage: floor(RR(log(2**100, 10)))                                            # needs sage.symbolic
            30
        """
        return "ComplexField(%s : Bits := true)" % self.prec()

    def _sage_input_(self, sib, coerce):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: sage_input(CIF, verify=True)
            # Verified
            CIF
            sage: sage_input(ComplexIntervalField(25), verify=True)
            # Verified
            ComplexIntervalField(25)
            sage: k = (CIF, ComplexIntervalField(37), ComplexIntervalField(1024))
            sage: sage_input(k, verify=True)
            # Verified
            (CIF, ComplexIntervalField(37), ComplexIntervalField(1024))
            sage: sage_input((k, k), verify=True)
            # Verified
            CIF37 = ComplexIntervalField(37)
            CIF1024 = ComplexIntervalField(1024)
            ((CIF, CIF37, CIF1024), (CIF, CIF37, CIF1024))
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: ComplexIntervalField(2)._sage_input_(SageInputBuilder(), False)
            {call: {atomic:ComplexIntervalField}({atomic:2})}
        """
        if self.prec() == 53:
            return sib.name('CIF')

        v = sib.name('ComplexIntervalField')(sib.int(self.prec()))
        name = 'CIF%d' % self.prec()
        sib.cache(self, v, name)
        return v

    precision = prec

    @cached_method
    def real_field(self):
        """
        Return the underlying :class:`RealIntervalField`.

        EXAMPLES::

            sage: R = CIF.real_field(); R
            Real Interval Field with 53 bits of precision
            sage: ComplexIntervalField(200).real_field()
            Real Interval Field with 200 bits of precision
            sage: CIF.real_field() is R
            True
        """
        return RealIntervalField(self._prec)

    # For compatibility with other complex number implementations
    # such as CC.
    _real_field = real_field

    @cached_method
    def middle_field(self):
        """
        Return the corresponding :class:`ComplexField` with the same precision
        as ``self``.

        EXAMPLES::

            sage: CIF.middle_field()
            Complex Field with 53 bits of precision
            sage: ComplexIntervalField(200).middle_field()
            Complex Field with 200 bits of precision
        """
        return ComplexField(self._prec)

    def __eq__(self, other):
        """
        Test whether ``self`` is equal to ``other``.

        If ``other`` is not a :class:`ComplexIntervalField_class`,
        return ``False``.  Otherwise, return ``True`` if ``self`` and
        ``other`` have the same precision.

        EXAMPLES::

            sage: CIF == ComplexIntervalField(200)
            False
            sage: CIF == CC
            False
            sage: CIF == CIF
            True
        """
        if not isinstance(other, ComplexIntervalField_class):
            return False
        return self._prec == other._prec

    def __hash__(self):
        """
        Return the hash.

        EXAMPLES::

            sage: C = ComplexIntervalField(200)
            sage: from sage.rings.complex_interval_field import ComplexIntervalField_class
            sage: D = ComplexIntervalField_class(200)
            sage: hash(C) == hash(D)
            True
        """
        return hash((self.__class__, self._prec))

    def __ne__(self, other):
        """
        Test whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: CIF != ComplexIntervalField(200)
            True
            sage: CIF != CC
            True
            sage: CIF != CIF
            False
        """
        return not (self == other)

    def __call__(self, x=None, im=None, **kwds):
        """
        Construct an element.

        EXAMPLES::

            sage: CIF(2) # indirect doctest
            2
            sage: CIF(CIF.0)
            1*I
            sage: CIF('1+I')
            Traceback (most recent call last):
            ...
            TypeError: unable to convert '1+I' to real interval
            sage: CIF(2,3)
            2 + 3*I
            sage: CIF(pi, e)                                                            # needs sage.symbolic
            3.141592653589794? + 2.718281828459046?*I
            sage: ComplexIntervalField(100)(CIF(RIF(2,3)))
            3.?

            sage: QQi.<i> = QuadraticField(-1)
            sage: CIF(i)
            1*I
            sage: QQi.<i> = QuadraticField(-1, embedding=CC(0,-1))
            sage: CIF(i)
            -1*I
            sage: QQi.<i> = QuadraticField(-1, embedding=None)
            sage: CIF(i)
            1*I

        ::

            sage: R.<x> = CIF[]
            sage: a = R(CIF(0,1)); a
            I
            sage: CIF(a)
            1*I
        """
        # Note: we override Parent.__call__ because we want to support
        # CIF(a, b) and that is hard to do using coerce maps.
        if im is not None or kwds:
            return self.element_class(self, x, im, **kwds)
        return Parent.__call__(self, x)

    def _coerce_map_from_(self, S):
        """
        Canonical coercion from ``S`` to this complex interval field.

        The rings that canonically coerce to the MPFI complex field are:

        - this MPFI complex field, or any other of higher precision

        - anything that canonically coerces to the real interval field
          with this precision

        - some exact or lazy parents representing subsets of the complex
          numbers, such as ``QQbar`` and ``CLF``.

        EXAMPLES::

            sage: CIF((2,1)) + 2 + I # indirect doctest
            4 + 2*I
            sage: CIF((2,1)) + CC.pi()
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for +: 'Complex Interval Field with 53 bits of precision' and 'Complex Field with 53 bits of precision'

            sage: CIF.coerce_map_from(QQ)
            Coercion map:
              From: Rational Field
              To:   Complex Interval Field with 53 bits of precision
            sage: CIF.coerce_map_from(int)
            Coercion map:
              From: Set of Python objects of class 'int'
              To:   Complex Interval Field with 53 bits of precision
            sage: CIF.coerce_map_from(GaussianIntegers())                               # needs sage.rings.number_field
            Conversion via _complex_mpfi_ method map:
              From: Gaussian Integers generated by I in Number Field in I
               with defining polynomial x^2 + 1 with I = 1*I
              To:   Complex Interval Field with 53 bits of precision
            sage: CIF.coerce_map_from(QQbar)                                            # needs sage.rings.number_field
            Conversion via _complex_mpfi_ method map:
              From: Algebraic Field
              To:   Complex Interval Field with 53 bits of precision
            sage: CIF.coerce_map_from(AA)                                               # needs sage.rings.number_field
            Conversion via _complex_mpfi_ method map:
              From: Algebraic Real Field
              To:   Complex Interval Field with 53 bits of precision
            sage: CIF.coerce_map_from(UniversalCyclotomicField())                       # needs sage.libs.gap sage.rings.number_field
            Conversion via _complex_mpfi_ method map:
              From: Universal Cyclotomic Field
              To:   Complex Interval Field with 53 bits of precision

        TESTS::

            sage: CIF.has_coerce_map_from(RR)
            False
            sage: CIF.has_coerce_map_from(RDF)
            False
            sage: CIF.has_coerce_map_from(float)
            False
        """
        # Direct and efficient conversions
        if S is ZZ or S is QQ or S is int:
            return True
        if isinstance(S, (ComplexIntervalField_class,
                          RealIntervalField_class)):
            return S.precision() >= self._prec

        # If coercion to CC is possible and there is a _complex_mpfi_
        # method, assume that it defines a coercion to CIF
        if self.middle_field().has_coerce_map_from(S):
            f = self._convert_method_map(S)
            if f is not None:
                return f

        return self._coerce_map_via( (self.real_field(),), S)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: ComplexIntervalField()  # indirect doctest
            Complex Interval Field with 53 bits of precision
            sage: ComplexIntervalField(100)  # indirect doctest
            Complex Interval Field with 100 bits of precision
        """
        return "Complex Interval Field with %s bits of precision" % self._prec

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(ComplexIntervalField()) # indirect doctest
            \Bold{C}
        """
        return "\\Bold{C}"

    def characteristic(self):
        """
        Return the characteristic of the complex (interval) field, which is 0.

        EXAMPLES::

            sage: CIF.characteristic()
            0
        """
        return integer.Integer(0)

    def gen(self, n=0):
        """
        Return the generator of the complex (interval) field.

        EXAMPLES::

            sage: CIF.0
            1*I
            sage: CIF.gen(0)
            1*I
        """
        if n != 0:
            raise IndexError("n must be 0")
        return self.element_class(self, 0, 1)

    def random_element(self, *args, **kwds):
        """
        Create a random element of ``self``.

        This simply chooses the real and imaginary part randomly, passing
        arguments and keywords to the underlying real interval field.

        EXAMPLES::

            sage: CIF.random_element().parent() is CIF
            True
            sage: re, im = CIF.random_element(10, 20)
            sage: 10 <= re <= 20
            True
            sage: 10 <= im <= 20
            True

        Passes extra positional or keyword arguments through::

            sage: re, im = CIF.random_element(max=0, min=-5)
            sage: -5 <= re <= 0
            True
            sage: -5 <= im <= 0
            True
        """
        rand = self.real_field().random_element
        re = rand(*args, **kwds)
        im = rand(*args, **kwds)
        return self.element_class(self, re, im)

    def is_field(self, proof=True):
        """
        Return ``True``, since the complex numbers are a field.

        EXAMPLES::

            sage: CIF.is_field()
            True
        """
        return True

    def pi(self):
        r"""
        Return `\pi` as an element in the complex (interval) field.

        EXAMPLES::

            sage: ComplexIntervalField(100).pi()
            3.14159265358979323846264338328?
        """
        return self.element_class(self, self.real_field().pi())

    def ngens(self):
        r"""
        The number of generators of this complex (interval) field as an
        `\RR`-algebra.

        There is one generator, namely ``sqrt(-1)``.

        EXAMPLES::

            sage: CIF.ngens()
            1
        """
        return 1

    def zeta(self, n=2):
        r"""
        Return a primitive `n`-th root of unity.

        .. TODO::

            Implement :class:`ComplexIntervalFieldElement` multiplicative order
            and set this output to have multiplicative order ``n``.

        INPUT:

        - ``n`` -- integer (default: 2)

        OUTPUT: a complex `n`-th root of unity

        EXAMPLES::

            sage: CIF.zeta(2)
            -1
            sage: CIF.zeta(5)
            0.309016994374948? + 0.9510565162951536?*I
        """
        from .integer import Integer
        n = Integer(n)
        if n == 1:
            x = self.element_class(self, 1)
        elif n == 2:
            x = self.element_class(self, -1)
        elif n >= 3:
            # Use De Moivre
            # e^(2*pi*i/n) = cos(2pi/n) + i *sin(2pi/n)
            RIF = self.real_field()
            pi = RIF.pi()
            z = 2*pi/n
            x = self.element_class(self, z.cos(), z.sin())
        # Uncomment after implemented
        #x._set_multiplicative_order( n )
        return x

    def scientific_notation(self, status=None):
        """
        Set or return the scientific notation printing flag.

        If this flag is ``True`` then complex numbers with this space as parent
        print using scientific notation.

        EXAMPLES::

            sage: CIF((0.025, 2))
            0.025000000000000002? + 2*I
            sage: CIF.scientific_notation(True)
            sage: CIF((0.025, 2))
            2.5000000000000002?e-2 + 2*I
            sage: CIF.scientific_notation(False)
            sage: CIF((0.025, 2))
            0.025000000000000002? + 2*I
        """
        return self.real_field().scientific_notation(status)
