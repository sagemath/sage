r"""
Arbitrary precision complex balls

This is an incomplete interface to the `acb module of FLINT <https://flintlib.org/doc/acb.html>`_;
it may be useful to refer to its documentation for more details.

Parts of the documentation for this module are copied or adapted from Arb's
(now FLINT's) own documentation, licenced at the time under the GNU General
Public License version 2, or later.

.. SEEALSO::

    - :mod:`Real balls <sage.rings.real_arb>`
    - :mod:`Complex interval field (using MPFI) <sage.rings.complex_interval_field>`
    - :mod:`Complex intervals (using MPFI) <sage.rings.complex_interval>`

Data Structure
==============

A :class:`ComplexBall` represents a complex number with error bounds. It wraps
an object of type ``acb_t``, which  consists of a pair of real number balls
representing the real and imaginary part with separate error bounds. (See the
documentation of :mod:`sage.rings.real_arb` for more information.)

A :class:`ComplexBall` thus represents a rectangle `[m_1-r_1, m_1+r_1] +
[m_2-r_2, m_2+r_2] i` in the complex plane. This is used instead of a
disk or square representation (consisting of a complex floating-point midpoint
with a single radius), since it allows implementing many operations more
conveniently by splitting into ball operations on the real and imaginary parts.
It also allows tracking when complex numbers have an exact (for example exactly
zero) real part and an inexact imaginary part, or vice versa.

The parents of complex balls are instances of :class:`ComplexBallField`.
The name ``CBF`` is bound to the complex ball field with the default precision
of 53 bits::

    sage: CBF is ComplexBallField() is ComplexBallField(53)
    True

Comparison
==========

.. WARNING::

    In accordance with the semantics of FLINT/Arb, identical :class:`ComplexBall`
    objects are understood to give permission for algebraic simplification.
    This assumption is made to improve performance. For example, setting ``z =
    x*x`` sets `z` to a ball enclosing the set `\{t^2 : t \in x\}` and not the
    (generally larger) set `\{tu : t \in x, u \in x\}`.

Two elements are equal if and only if they are exact and equal (in spite of the
above warning, inexact balls are not considered equal to themselves)::

    sage: a = CBF(1, 2)
    sage: b = CBF(1, 2)
    sage: a is b
    False
    sage: a == a
    True
    sage: a == b
    True

::

    sage: a = CBF(1/3, 1/5)
    sage: b = CBF(1/3, 1/5)
    sage: a.is_exact()
    False
    sage: b.is_exact()
    False
    sage: a is b
    False
    sage: a == a
    False
    sage: a == b
    False

A ball is nonzero in the sense of usual comparison if and only if it does not
contain zero::

    sage: a = CBF(RIF(-0.5, 0.5))
    sage: a != 0
    False
    sage: b = CBF(1/3, 1/5)
    sage: b != 0
    True

However, ``bool(b)`` returns ``False`` for a ball ``b`` only if ``b`` is exactly
zero::

    sage: bool(a)
    True
    sage: bool(b)
    True
    sage: bool(CBF.zero())
    False

Coercion
========

Automatic coercions work as expected::

    sage: # needs sage.symbolic
    sage: bpol = 1/3*CBF(i) + AA(sqrt(2))
    sage: bpol += polygen(RealBallField(20), 'x') + QQbar(i)
    sage: bpol
    x + [1.41421 +/- ...e-6] + [1.33333 +/- ...e-6]*I
    sage: bpol.parent()
    Univariate Polynomial Ring in x over Complex ball field with 20 bits of precision
    sage: bpol/3
    ([0.333333 +/- ...e-7])*x + [0.47140 +/- ...e-6] + [0.44444 +/- ...e-6]*I

TESTS::

    sage: polygen(CBF, 'x')^3
    x^3

::

    sage: SR.coerce(CBF(0.42 + 3.33*I))                                                 # needs sage.symbolic
    [0.4200000000000000 +/- ...e-17] + [3.330000000000000 +/- ...e-17]*I

Check that :issue:`19839` is fixed::

    sage: log(SR(CBF(0.42))).pyobject().parent()                                        # needs sage.symbolic
    Complex ball field with 53 bits of precision

:issue:`24621`::

    sage: CBF(NumberField(polygen(QQ, 'y')^3 + 20, 'a', embedding=CC(1.35,2.35)).gen())
    [1.35720880829745...] + [2.35075461245119...]*I

Classes and Methods
===================
"""
#*****************************************************************************
# Copyright (C) 2014 Clemens Heuberger <clemens.heuberger@aau.at>
#               2017 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                https://www.gnu.org/licenses/
#*****************************************************************************

import operator
import sys
import warnings
from cysignals.signals cimport sig_on, sig_str, sig_off, sig_error, sig_block, sig_unblock

import sage.categories.fields

cimport sage.rings.abc
cimport sage.rings.rational

from cpython.long cimport PyLong_AsLong
from cpython.object cimport Py_EQ, Py_NE
from cpython.complex cimport PyComplex_FromDoubles

from sage.ext.stdsage cimport PY_NEW

from sage.libs.mpfr cimport MPFR_RNDU, MPFR_RNDD, MPFR_PREC_MIN, mpfr_get_d_2exp
from sage.libs.flint.types cimport ARF_RND_NEAR, arf_t, mag_t
from sage.libs.flint.arb cimport *
from sage.libs.flint.acb cimport *
from sage.libs.flint.acb_calc cimport *
from sage.libs.flint.acb_hypgeom cimport *
from sage.libs.flint.acb_elliptic cimport *
from sage.libs.flint.acb_modular cimport *
from sage.libs.flint.acb_poly cimport *
from sage.libs.flint.arf cimport arf_init, arf_get_d, arf_get_mpfr, arf_clear, arf_set, arf_is_nan
from sage.libs.flint.mag cimport (mag_init, mag_clear, mag_set_d,
        MAG_BITS, mag_zero, mag_set_ui_2exp_si,
        mag_mul_2exp_si)
from sage.libs.flint.fmpz cimport fmpz_t, fmpz_init, fmpz_get_mpz, fmpz_set_mpz, fmpz_clear
from sage.libs.flint.fmpq cimport fmpq_t, fmpq_init, fmpq_set_mpq, fmpq_clear
from sage.libs.gmp.mpz cimport mpz_fits_slong_p, mpz_get_si
from sage.libs.gsl.complex cimport gsl_complex_rect
from sage.rings.real_double cimport RealDoubleElement
from sage.rings.complex_double cimport ComplexDoubleElement
from sage.rings.integer cimport Integer
from sage.rings.rational_field import QQ
from sage.rings.number_field.number_field_base import NumberField
from sage.rings.polynomial.polynomial_complex_arb cimport Polynomial_complex_arb
from sage.rings.real_arb cimport mpfi_to_arb, arb_to_mpfi
from sage.rings.real_arb import RealBallField
from sage.rings.real_mpfi cimport RealIntervalField_class
from sage.rings.real_mpfr cimport RealField_class, RealField, RealNumber
from sage.rings.ring import Field
from sage.structure.element cimport Element
from sage.structure.unique_representation import UniqueRepresentation
from sage.arith.long cimport is_small_python_int

from sage.misc.lazy_string import lazy_string
from sage.rings.complex_interval_field import ComplexIntervalField, ComplexIntervalField_class
from sage.rings.integer_ring import ZZ


cdef void ComplexIntervalFieldElement_to_acb(
    acb_t target,
    ComplexIntervalFieldElement source) noexcept:
    """
    Convert a :class:`ComplexIntervalFieldElement` to an ``acb``.

    INPUT:

    - ``target`` -- an ``acb_t``

    - ``source`` -- a :class:`ComplexIntervalFieldElement`

    OUTPUT:

    None.
    """
    cdef long precision
    precision = source.parent().precision()
    mpfi_to_arb(acb_realref(target), source.__re, precision)
    mpfi_to_arb(acb_imagref(target), source.__im, precision)

cdef int acb_to_ComplexIntervalFieldElement(
    ComplexIntervalFieldElement target,
    const acb_t source) except -1:
    """
    Convert an ``acb`` to a :class:`ComplexIntervalFieldElement`.

    INPUT:

    - ``target`` -- a :class:`ComplexIntervalFieldElement`

    - ``source`` -- an ``acb_t``

    OUTPUT:

    A :class:`ComplexIntervalFieldElement`.
    """
    cdef long precision = target._prec

    arb_to_mpfi(target.__re, acb_realref(source), precision)
    arb_to_mpfi(target.__im, acb_imagref(source), precision)
    return 0

cdef class IntegrationContext:
    r"""
    Used to wrap the integrand and hold some context information during
    numerical integration.
    """
    cdef object f
    cdef object parent
    cdef object exn_type
    cdef object exn_obj
    cdef object exn_tb

cdef int acb_calc_func_callback(acb_ptr out, const acb_t inp, void * param,
        long order, long prec) noexcept:
    r"""
    Callback used for numerical integration

    TESTS::

        sage: CBF.integral(lambda x, flag: 24, 0, 2)
        48.00000000000000

        sage: CBF.integral(lambda x, flag: "a", 0, 1)
        Traceback (most recent call last):
        ...
        TypeError: no canonical coercion ... to Complex ball field with 53 bits
        of precision

        sage: def foo(*args):
        ....:     raise RuntimeError
        sage: CBF.integral(foo, 0, 2)
        Traceback (most recent call last):
        ...
        RuntimeError

        sage: points = []
        sage: def foo(x, flag):
        ....:     points.append(x)
        ....:     return x
        sage: CBF.integral(foo, 0, 1)
        [0.50000000000000...]
        sage: points
        [[+/- 1.01], ..., [0.788...], [0.211...]]
    """
    cdef IntegrationContext ctx
    cdef ComplexBall x
    sig_block()
    try:
        ctx = <IntegrationContext>param
        if ctx.exn_type is not None or order >= 2:
            acb_indeterminate(out)
            return 0
        try:
            x = ComplexBall.__new__(ComplexBall)
            assert prec == ctx.parent._prec
            x._parent = ctx.parent
            acb_set(x.value, inp)
            y = ctx.f(x, (order == 1))
            if not isinstance(y, ComplexBall):
                y = ctx.parent.coerce(y)
            acb_set(out, (<ComplexBall> y).value)
        except BaseException:
            ctx.exn_type, ctx.exn_obj, ctx.exn_tb = sys.exc_info()
            acb_indeterminate(out)
        return 0
    finally:
        sig_unblock()


class ComplexBallField(UniqueRepresentation, sage.rings.abc.ComplexBallField):
    r"""
    An approximation of the field of complex numbers using pairs of mid-rad
    intervals.

    INPUT:

    - ``precision`` -- integer `\ge 2`

    EXAMPLES::

        sage: CBF(1)
        1.000000000000000

    TESTS::

        sage: ComplexBallField(0)
        Traceback (most recent call last):
        ...
        ValueError: precision must be at least 2
        sage: ComplexBallField(1)
        Traceback (most recent call last):
        ...
        ValueError: precision must be at least 2

        sage: ComplexBallField().is_finite()
        False

        sage: loads(dumps(ComplexBallField(60))) is ComplexBallField(60)
        True

    .. SEEALSO::

        - :mod:`sage.rings.complex_arb`
        - :mod:`sage.rings.complex_mpfr`
        - :mod:`sage.rings.complex_mpfi`
        - :mod:`sage.rings.real_arb`
    """
    Element = ComplexBall

    @staticmethod
    def __classcall__(cls, long precision=53):
        r"""
        Normalize the arguments for caching.

        TESTS::

            sage: ComplexBallField(53) is ComplexBallField()
            True
        """
        return super().__classcall__(cls, precision)

    def __init__(self, long precision=53):
        r"""
        Initialize the complex ball field.

        INPUT:

        - ``precision`` -- integer `\ge 2`

        EXAMPLES::

            sage: CBF(1)
            1.000000000000000

        TESTS::

            sage: CBF.base()
            Real ball field with 53 bits of precision
            sage: CBF.base_ring()
            Real ball field with 53 bits of precision

        There are direct coercions from ZZ and QQ (for which FLINT provides
        construction functions)::

            sage: CBF.coerce_map_from(ZZ)
            Coercion map:
              From: Integer Ring
              To:   Complex ball field with 53 bits of precision
            sage: CBF.coerce_map_from(QQ)
            Coercion map:
              From: Rational Field
              To:   Complex ball field with 53 bits of precision

        Various other coercions are available through real ball fields or CLF::

            sage: CBF.has_coerce_map_from(RLF)
            True
            sage: CBF.has_coerce_map_from(AA)
            True
            sage: CBF.has_coerce_map_from(QuadraticField(-1))
            True
            sage: CBF.has_coerce_map_from(QQbar)
            True
            sage: CBF.has_coerce_map_from(CLF)
            True
        """
        if precision < 2:
            raise ValueError("precision must be at least 2")
        self._prec = precision
        real_field = RealBallField(self._prec)
        Field.__init__(self,
                base_ring=real_field,
                category=sage.categories.fields.Fields().Infinite())
        from sage.rings.rational_field import QQ
        self._populate_coercion_lists_([ZZ, QQ], convert_method_name='_acb_')

    def _real_field(self):
        """
        TESTS::

            sage: CBF._real_field() is RBF
            True
        """
        return self._base

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: ComplexBallField()
            Complex ball field with 53 bits of precision
            sage: ComplexBallField(106)
            Complex ball field with 106 bits of precision
        """
        return "Complex ball field with {} bits of precision".format(self._prec)

    def construction(self):
        """
        Return the construction of a complex ball field as the algebraic
        closure of the real ball field with the same precision.

        EXAMPLES::

            sage: functor, base = CBF.construction()
            sage: functor, base
            (AlgebraicClosureFunctor, Real ball field with 53 bits of precision)
            sage: functor(base) is CBF
            True
        """
        from sage.categories.pushout import AlgebraicClosureFunctor
        return (AlgebraicClosureFunctor(), self._base)

    def complex_field(self):
        """
        Return the complex ball field with the same precision, i.e. ``self``.

        EXAMPLES::

            sage: CBF.complex_field() is CBF
            True
        """
        return ComplexBallField(self._prec)

    def ngens(self):
        r"""
        Return 1 as the only generator is the imaginary unit.

        EXAMPLES::

            sage: CBF.ngens()
            1
        """
        return 1

    def gen(self, i):
        r"""
        For i = 0, return the imaginary unit in this complex ball field.

        EXAMPLES::

            sage: CBF.0
            1.000000000000000*I
            sage: CBF.gen(1)
            Traceback (most recent call last):
            ...
            ValueError: only one generator
        """
        if i == 0:
            return self(0, 1)
        else:
            raise ValueError("only one generator")

    def gens(self):
        r"""
        Return the tuple of generators of this complex ball field, i.e.
        ``(i,)``.

        EXAMPLES::

            sage: CBF.gens()
            (1.000000000000000*I,)
            sage: CBF.gens_dict()
            {'1.000000000000000*I': 1.000000000000000*I}
        """
        return (self(0, 1),)

    def _coerce_map_from_(self, other):
        r"""
        Parents that canonically coerce into complex ball fields include:

        - anything that coerces into the corresponding real ball field;

        - real and complex ball fields with a larger precision;

        - various exact or lazy parents representing subsets of the complex
          numbers, such as ``QQbar``, ``CLF``, and number fields equipped
          with complex embeddings.

        TESTS::

            sage: CBF.coerce_map_from(CBF)
            Identity endomorphism of Complex ball field with 53 bits of precision
            sage: CBF.coerce_map_from(ComplexBallField(100))
            Coercion map:
              From: Complex ball field with 100 bits of precision
              To:   Complex ball field with 53 bits of precision
            sage: CBF.has_coerce_map_from(ComplexBallField(42))
            False
            sage: CBF.has_coerce_map_from(RealBallField(54))
            True
            sage: CBF.has_coerce_map_from(RealBallField(52))
            False
            sage: CBF.has_coerce_map_from(QuadraticField(-2))
            True
            sage: CBF.has_coerce_map_from(QuadraticField(2, embedding=None))
            False

        Check that there are no coercions from interval or floating-point parents::

            sage: CBF.has_coerce_map_from(RIF)
            False
            sage: CBF.has_coerce_map_from(CIF)
            False
            sage: CBF.has_coerce_map_from(RR)
            False
            sage: CBF.has_coerce_map_from(CC)
            False

        Check that the map goes through the ``_acb_`` method::

            sage: CBF.coerce_map_from(QuadraticField(-2, embedding=AA(-2).sqrt()))
            Conversion via _acb_ method map:
            ...
            sage: CBF.convert_map_from(QuadraticField(-2))
            Conversion via _acb_ method map:
            ...
            sage: x = polygen(ZZ, 'x')
            sage: CBF.coerce_map_from(NumberField(x^7 + 2, 'a',
            ....:                                 embedding=QQbar(-2)^(1/7)))
            Conversion via _acb_ method map:
            ...
            sage: CBF.coerce_map_from(QQbar)
            Conversion via _acb_ method map:
            ...
        """
        if isinstance(other, RealBallField):
            return other._prec >= self._prec
        elif isinstance(other, ComplexBallField):
            return other._prec >= self._prec

        if other is not QQ and isinstance(other, NumberField):
            emb = other.coerce_embedding()
            return emb is not None and self.has_coerce_map_from(emb.codomain())

        from sage.rings.qqbar import AA, QQbar
        from sage.rings.real_lazy import RLF, CLF

        if other in [AA, QQbar, RLF, CLF]:
            return True

    def _element_constructor_(self, x=None, y=None):
        r"""
        Convert (x, y) to an element of this complex ball field, perhaps
        non-canonically.

        INPUT:

        - ``x``, ``y`` -- (optional) either a complex number, interval or ball,
          or two real ones (see examples below for more information on accepted
          number types).

        EXAMPLES::

            sage: CBF()
            0
            sage: CBF(1) # indirect doctest
            1.000000000000000
            sage: CBF(1, 1)
            1.000000000000000 + 1.000000000000000*I
            sage: CBF(pi, sqrt(2))                                                      # needs sage.symbolic
            [3.141592653589793 +/- ...e-16] + [1.414213562373095 +/- ...e-16]*I
            sage: CBF(I)
            1.000000000000000*I
            sage: CBF(pi + I/3)                                                         # needs sage.symbolic
            [3.141592653589793 +/- ...e-16] + [0.3333333333333333 +/- ...e-17]*I
            sage: CBF(QQbar(i/7)) # abs tol 1e-16
            [0.1428571428571429 +/- 4.29e-17]*I
            sage: CBF(AA(sqrt(2)))                                                      # needs sage.symbolic
            [1.414213562373095 +/- ...e-16]
            sage: CBF(CIF(0, 1))
            1.000000000000000*I
            sage: CBF(RBF(1/3))
            [0.3333333333333333 +/- ...e-17]
            sage: CBF(RBF(1/3), RBF(1/6))
            [0.3333333333333333 +/- ...e-17] + [0.1666666666666667 +/- ...e-17]*I
            sage: CBF(1/3)
            [0.3333333333333333 +/- ...e-17]
            sage: CBF(1/3, 1/6)
            [0.3333333333333333 +/- ...e-17] + [0.1666666666666667 +/- ...e-17]*I
            sage: ComplexBallField(106)(1/3, 1/6)
            [0.33333333333333333333333333333333 +/- ...e-33] + [0.16666666666666666666666666666666 +/- ...e-33]*I
            sage: NF.<a> = QuadraticField(-2)
            sage: CBF(1/5 + a/2)
            [0.2000000000000000 +/- ...e-17] + [0.707106781186547 +/- ...e-16]*I
            sage: CBF(infinity, NaN)                                                    # needs sage.symbolic
            [+/- inf] + nan*I
            sage: CBF(x)                                                                # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: unable to convert x to a ComplexBall

        .. SEEALSO::

            :meth:`sage.rings.real_arb.RealBallField._element_constructor_`

        TESTS::

            sage: CBF(1+I, 2)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert ... to a RealBall

        The following conversions used to yield incorrect enclosures::

            sage: # needs sage.symbolic
            sage: a = CBF(airy_ai(1)); a
            [0.1352924163128814 +/- 6.95e-17]
            sage: a.overlaps(ComplexBallField(100).one().airy_ai())
            True
            sage: v = CBF(zetaderiv(1, 3/2)); v
            [-3.932239737431101 +/- 5.58e-16]
            sage: v.overlaps(ComplexBallField(100)(3/2).zetaderiv(1))
            True
        """
        try:
            return self.element_class(self, x, y)
        except TypeError:
            pass

        if y is None:
            try:
                _x = self._base(x)
                return self.element_class(self, _x)
            except (TypeError, ValueError):
                pass
            # Handle symbolic expressions in a special way in order to avoid
            # unsafe conversions as much as possible. Unlike the real case,
            # this is not implemented via an _acb_() method, because such a
            # conversion method would also be called by things like
            # CBF(re_expr, im_expr).
            from sage.structure.element import Expression
            if isinstance(x, Expression):
                # Parse the expression. Despite the method name, the result
                # will be a complex ball.
                _x = x._arb_(self)
                return self.element_class(self, _x)
            try:
                _y = self._base(x.imag())
                _x = self._base(x.real())
                return self.element_class(self, _x, _y)
            except (AttributeError, TypeError):
                pass
            try:
                _x = ComplexIntervalField(self._prec)(x)
                return self.element_class(self, _x)
            except TypeError:
                pass
            raise TypeError(lazy_string("unable to convert %s to a ComplexBall", x))
        else:
            _x = self._base(x)
            _y = self._base(y)
            return self.element_class(self, _x, _y)

    def _an_element_(self):
        r"""
        Construct an element.

        EXAMPLES::

            sage: CBF.an_element() # indirect doctest
            [0.3333333333333333 +/- ...e-17] - [0.1666666666666667 +/- ...e-17]*I
        """
        return self(1.0/3, -1.0/6)

    def precision(self):
        """
        Return the bit precision used for operations on elements of this field.

        EXAMPLES::

            sage: ComplexBallField().precision()
            53
        """
        return self._prec

    prec = precision

    def is_exact(self):
        """
        Complex ball fields are not exact.

        EXAMPLES::

            sage: ComplexBallField().is_exact()
            False
        """
        return False

    def characteristic(self):
        """
        Complex ball fields have characteristic zero.

        EXAMPLES::

            sage: ComplexBallField().characteristic()
            0
        """
        return 0

    def some_elements(self):
        """
        Complex ball fields contain elements with exact, inexact, infinite, or
        undefined real and imaginary parts.

        EXAMPLES::

            sage: CBF.some_elements()
            [1.000000000000000,
             -0.5000000000000000*I,
             1.000000000000000 + [0.3333333333333333 +/- ...e-17]*I,
             [-0.3333333333333333 +/- ...e-17] + 0.2500000000000000*I,
             [-2.175556475109056e+181961467118333366510562 +/- ...e+181961467118333366510545],
             [+/- inf],
             [0.3333333333333333 +/- ...e-17] + [+/- inf]*I,
             [+/- inf] + [+/- inf]*I,
             nan,
             nan + nan*I,
             [+/- inf] + nan*I]
        """
        return [self(1), self(0, -1./2), self(1, 1./3), self(-1./3, 1./4),
                -self(1, 1)**(Integer(2)**80),
                self('inf'), self(1./3, 'inf'), self('inf', 'inf'),
                self('nan'), self('nan', 'nan'), self('inf', 'nan')]

    def _roots_univariate_polynomial(self, pol, ring, multiplicities,
                                     algorithm, proof=True):
        r"""
        Compute the roots of ``pol``.

        This method is used internally by the
        :meth:`sage.rings.polynomial.polynomial_element.Polynomial.roots`
        method of polynomials with complex ball coefficients. See its
        documentation for details.

        EXAMPLES::

            sage: import warnings
            sage: warnings.simplefilter("always")

            sage: Pol.<x> = CBF[]
            sage: i = CBF.gen(0)

            sage: (x^4 - 1/3).roots()
            Traceback (most recent call last):
            ...
            ValueError: polynomial with interval coefficients, use multiplicities=False

            sage: (x^4 - 1/3).roots(multiplicities=False) # indirect doctest
            [[-0.759835685651593 +/- ...e-16] + [+/- ...e-16]*I,
             [0.759835685651593 +/- ...e-16] + [+/- ...e-16]*I,
             [+/- ...e-16] + [0.759835685651593 +/- ...e-16]*I,
             [+/- ...e-16] + [-0.759835685651593 +/- ...e-16]*I]

            sage: (x^4 - 1/3).roots(RBF, multiplicities=False)
            [[-0.759835685651593 +/- ...e-16], [0.759835685651593 +/- ...e-16]]

            sage: (x^4 - 3).roots(RealBallField(100), multiplicities=False)
            [[-1.316074012952492460819218901797 +/- ...e-34],
             [1.316074012952492460819218901797 +/- ...e-34]]

            sage: (x^4 - 3).roots(ComplexIntervalField(100), multiplicities=False)
            [-1.31607401295249246081921890180? + 0.?e-37*I,
             1.31607401295249246081921890180? + 0.?e-37*I,
             0.?e-37 + 1.31607401295249246081921890180?*I,
             0.?e-37 - 1.31607401295249246081921890180?*I]

            sage: (x^2 - i/3).roots(ComplexBallField(2), multiplicities=False)
            [[+/- 0.409] + [+/- 0.409]*I, [+/- 0.409] + [+/- 0.409]*I]

            sage: ((x - 1)^2).roots(multiplicities=False)
            Traceback (most recent call last):
            ...
            ValueError: unable to isolate the roots (try using proof=False or
            increasing the precision)
            sage: ((x - 1)^2).roots(multiplicities=False, proof=False)
            doctest:...
            UserWarning: roots may have been lost...
            [[1.00000000000 +/- ...e-12] + [+/- ...e-11]*I,
             [1.0000000000 +/- ...e-12] + [+/- ...e-12]*I]

            sage: pol = x^7 - 2*(1000*x - 1)^2 # Mignotte polynomial
            sage: pol.roots(multiplicities=False)
            Traceback (most recent call last):
            ...
            ValueError: unable to isolate the roots (try using proof=False or
            increasing the precision)
            sage: sorted(pol.roots(multiplicities=False, proof=False), key=str)
            doctest:...
            UserWarning: roots may have been lost...
            [[-14.72907378354557 +/- ...e-15] + [-10.70100790294238 +/- ...e-15]*I,
             [-14.72907378354557 +/- ...e-15] + [10.70100790294238 +/- ...e-15]*I,
             [0.00100000 +/- ...e-10] + [+/- ...e-10]*I,
             [0.001000000 +/- ...e-10] + [+/- ...e-10]*I,
             [18.20524201487994 +/- ...e-15] + [+/- ...e-37]*I,
             [5.625452776105595 +/- ...e-16] + [-17.31459450084417 +/- ...e-15]*I,
             [5.625452776105595 +/- ...e-16] + [17.31459450084417 +/- ...e-15]*I]
            sage: pol.roots(ComplexBallField(100), multiplicities=False)
            [[0.00099999999997763932022675...] + [+/- ...]*I,
             ...]

            sage: ((x - 1)^2 + 2^(-70)*i/3).roots(RBF, multiplicities=False)
            Traceback (most recent call last):
            ...
            ValueError: unable to determine which roots are real

        TESTS::

            sage: CBF._roots_univariate_polynomial(CBF['x'].zero(), CBF, False, None)
            Traceback (most recent call last):
            ...
            ArithmeticError: taking the roots of the zero polynomial
        """
        if algorithm is not None:
            raise NotImplementedError
        if multiplicities:
            raise ValueError("polynomial with interval coefficients, "
                             "use multiplicities=False")

        cdef bint real = False
        if ring is None:
            ring = self
        elif isinstance(ring, (ComplexBallField, ComplexIntervalField_class)):
            pass
        elif isinstance(ring, (RealBallField, RealIntervalField_class)):
            real = True
        elif ring.has_coerce_map_from(self):
            pass
        elif (ring.has_coerce_map_from(self._base)
                and RealField(MPFR_PREC_MIN).has_coerce_map_from(ring)):
            real = True
        else:
            raise NotImplementedError

        cdef Polynomial_complex_arb poly = <Polynomial_complex_arb?> pol
        cdef acb_poly_t rounded_poly
        cdef long tgtprec = ring.precision()
        cdef long maxprec = 3*max(poly._parent._base._prec, tgtprec)
        cdef long initial_prec = min(32, maxprec)
        cdef long prec = initial_prec
        cdef long isolated = 0
        cdef RealBall rb
        cdef ComplexBall cb
        acb_poly_init(rounded_poly)
        cdef long deg = acb_poly_degree(poly._poly)
        if deg < 0:
            raise ArithmeticError("taking the roots of the zero polynomial")
        cdef acb_ptr roots = _acb_vec_init(deg)
        try:
            sig_on()
            while ((isolated < deg or any(acb_rel_accuracy_bits(&roots[i]) < tgtprec
                                        for i in range(deg)))
                and prec < maxprec):
                acb_poly_set_round(rounded_poly, poly._poly, prec)
                maxiter = min(max(deg, 32), prec)
                if (prec == initial_prec):
                    isolated = acb_poly_find_roots(roots, rounded_poly, NULL, maxiter, prec)
                else:
                    isolated = acb_poly_find_roots(roots, rounded_poly, roots, maxiter, prec)
                prec *= 2
            sig_off()

            if isolated < deg:
                if proof:
                    raise ValueError("unable to isolate the roots (try using "
                            "proof=False or increasing the precision)")
                else:
                    warnings.warn("roots may have been lost")

            _acb_vec_sort_pretty(roots, deg)

            res = []
            if real:
                if not acb_poly_validate_real_roots(roots, rounded_poly, prec):
                    raise ValueError("unable to determine which roots are real")
                for i in range(deg):
                    if arb_contains_zero(acb_imagref(&roots[i])):
                        rb = RealBall.__new__(RealBall)
                        rb._parent = self._base
                        arb_set(rb.value, acb_realref(&roots[i]))
                        res.append(ring(rb))
            else:
                for i in range(deg):
                    cb = ComplexBall.__new__(ComplexBall)
                    cb._parent = self
                    acb_set(cb.value, &roots[i])
                    res.append(ring(cb))
        finally:
            _acb_vec_clear(roots, deg)
            acb_poly_clear(rounded_poly)
        return res

    def _sum_of_products(self, terms):
        r"""
        Compute a sum of product of complex balls without creating temporary
        Python objects

        The input objects should be complex balls, but need not belong to this
        parent. The computation is performed at the precision of this parent.

        EXAMPLES::

            sage: Pol.<x> = ComplexBallField(1000)[]
            sage: pol = (x + 1/3)^100
            sage: CBF._sum_of_products((c, c) for c in pol)
            [6.3308767660842e+23 +/- ...e+9]

        TESTS::

            sage: CBF._sum_of_products([])
            0
            sage: CBF._sum_of_products([[]])
            1.000000000000000
            sage: CBF._sum_of_products([["a"]])
            Traceback (most recent call last):
            ...
            TypeError: Cannot convert str to sage.rings.complex_arb.ComplexBall
        """
        cdef ComplexBall res = ComplexBall.__new__(ComplexBall)
        cdef ComplexBall factor
        cdef acb_t tmp
        res._parent = self
        acb_zero(res.value)
        acb_init(tmp)
        try:
            for term in terms:
                acb_one(tmp)
                for factor in term:
                    acb_mul(tmp, tmp, factor.value, self._prec)
                acb_add(res.value, res.value, tmp, self._prec)
        finally:
            acb_clear(tmp)
        return res

    # Constants

    def pi(self):
        r"""
        Return a ball enclosing `\pi`.

        EXAMPLES::

            sage: CBF.pi()
            [3.141592653589793 +/- ...e-16]
            sage: ComplexBallField(128).pi()
            [3.1415926535897932384626433832795028842 +/- ...e-38]

            sage: CBF.pi().parent()
            Complex ball field with 53 bits of precision
        """
        cdef ComplexBall res = ComplexBall.__new__(ComplexBall)
        res._parent = self
        if _do_sig(self._prec): sig_on()
        arb_const_pi(acb_realref(res.value), self._prec)
        arb_zero(acb_imagref(res.value))
        if _do_sig(self._prec): sig_off()
        return res

    def integral(self, func, a, b, params=None,
            rel_tol=None, abs_tol=None,
            deg_limit=None, eval_limit=None, depth_limit=None,
            use_heap=None, verbose=None):
        r"""
        Compute a rigorous enclosure of the integral of ``func`` on the
        interval [``a``, ``b``].

        INPUT:

        - ``func`` -- a callable object accepting two parameters, a complex
          ball ``x`` and a boolean flag ``analytic``, and returning an element
          of this ball field (or some value that coerces into this ball field),
          such that:

          - ``func(x, False)`` evaluates the integrand `f` on the ball ``x``.
            There are no restrictions on the behavior of `f` on ``x``; in
            particular, it can be discontinuous.

          - ``func(x, True)`` evaluates `f(x)` if  `f` is analytic on the
            whole ``x``, and returns some non-finite ball (e.g., ``self(NaN)``)
            otherwise.

          (The ``analytic`` flag only needs to be checked for integrands that
          are non-analytic but bounded in some regions, typically complex
          functions with branch cuts, like `\sqrt{z}`. In particular, it can be
          ignored for meromorphic functions.)

        - ``a``, ``b`` -- integration bounds. The bounds can be real or complex
          balls, or elements of any parent that coerces into this ball field,
          e.g. rational or algebraic numbers.

        - ``rel_tol`` -- relative accuracy goal (default: `2^{-p}` where `p` is
          the precision of the ball field)

        - ``abs_tol`` -- absolute accuracy goal (default: `2^{-p}` where `p` is
          the precision of the ball field)

        Additionally, the following optional parameters can be used to control
        the integration algorithm. See the `FLINT documentation <https://flintlib.org/doc/acb_calc.html>`_
        for more information.

        - ``deg_limit`` -- maximum quadrature degree for each
          subinterval

        - ``eval_limit`` -- maximum number of function
          evaluations

        - ``depth_limit`` -- maximum search depth for
          adaptive subdivision

        - ``use_heap`` -- boolean (default: ``False``); if ``True``, use a
          priority queue instead of a stack to manage subintervals. This
          sometimes gives better results for integrals with slow convergence but
          may require more memory and increasing ``depth_limit``.

        - ``verbose`` -- integer (default: 0); if set to 1, some information
          about the overall integration process is printed to standard
          output. If set to 2, information about each subinterval is printed.

        EXAMPLES:

        Some analytic integrands::

            sage: CBF.integral(lambda x, _: x, 0, 1)
            [0.500000000000000 +/- ...e-16]

            sage: CBF.integral(lambda x, _: x.gamma(), 1 - CBF(i), 1 + CBF(i)) # abs tol 1e-13
            [+/- 1.39e-15] + [1.57239266949806 +/- 8.33e-15]*I

            sage: C = ComplexBallField(100)
            sage: C.integral(lambda x, _: x.cos() * x.sin(), 0, 1)
            [0.35403670913678559674939205737 +/- ...e-30]

            sage: CBF.integral(lambda x, _: (x + x.exp()).sin(), 0, 8)
            [0.34740017266 +/- ...e-12]

            sage: C = ComplexBallField(2000)
            sage: C.integral(lambda x, _: (x + x.exp()).sin(), 0, 8) # long time
            [0.34740017...55347713 +/- ...e-598]

        Here the integration path crosses the branch cut of the square root::

            sage: def my_sqrt(z, analytic):
            ....:     if (analytic and not z.real() > 0
            ....:                  and z.imag().contains_zero()):
            ....:         return CBF(NaN)
            ....:     else:
            ....:         return z.sqrt()
            sage: CBF.integral(my_sqrt, -1 + CBF(i), -1 - CBF(i))                       # needs sage.symbolic
            [+/- ...e-14] + [-0.4752076627926 +/- 5...e-14]*I

        Note, though, that proper handling of the ``analytic`` flag is required
        even when the path does not touch the branch cut::

            sage: correct = CBF.integral(my_sqrt, 1, 2); correct
            [1.21895141649746 +/- ...e-15]
            sage: RBF(integral(sqrt(x), x, 1, 2))       # long time                     # needs sage.symbolic
            [1.21895141649746 +/- ...e-15]
            sage: wrong = CBF.integral(lambda z, _: z.sqrt(), 1, 2) # WRONG!
            sage: correct - wrong
            [-5.640636259e-5 +/- ...e-15]

        We can integrate the real absolute value function by defining a
        piecewise holomorphic extension::

            sage: def real_abs(z, analytic):
            ....:     if z.real().contains_zero():
            ....:         if analytic:
            ....:             return z.parent()(NaN)
            ....:         else:
            ....:             return z.union(-z)
            ....:     elif z.real() > 0:
            ....:         return z
            ....:     else:
            ....:         return -z
            sage: CBF.integral(real_abs, -1, 1)
            [1.00000000000...]
            sage: CBF.integral(lambda z, analytic: real_abs(z.sin(), analytic), 0, 2*CBF.pi())
            [4.00000000000...]

        Some methods of complex balls natively support the ``analytic`` flag::

            sage: CBF.integral(lambda z, analytic: z.log(analytic=analytic),
            ....:              -1-CBF(i), -1+CBF(i))
            [+/- ...e-14] + [0.26394350735484 +/- ...e-15]*I
            sage: from sage.rings.complex_arb import ComplexBall
            sage: CBF.integral(ComplexBall.sqrt, -1+CBF(i), -1-CBF(i))
            [+/- ...e-14] + [-0.4752076627926 +/- 5...e-14]*I

        Here the integrand has a pole on or very close to the integration path,
        but there is no need to explicitly handle the ``analytic`` flag since
        the integrand is unbounded::

            sage: CBF.integral(lambda x, _: 1/x, -1, 1)
            nan + nan*I
            sage: CBF.integral(lambda x, _: 1/x, 10^-1000, 1)
            nan + nan*I
            sage: CBF.integral(lambda x, _: 1/x, 10^-1000, 1, abs_tol=1e-10)
            [2302.5850930 +/- ...e-8]

        Tolerances::

            sage: CBF.integral(lambda x, _: x.exp(), -1020, -1010)
            [+/- ...e-438]
            sage: CBF.integral(lambda x, _: x.exp(), -1020, -1010, abs_tol=1e-450)
            [2.304377150950e-439 +/- ...e-452]
            sage: CBF.integral(lambda x, _: x.exp(), -1020, -1010, abs_tol=0)
            [2.304377150950e-439 +/- 7...e-452]
            sage: CBF.integral(lambda x, _: x.exp(), -1020, -1010, rel_tol=1e-2, abs_tol=0)
            [2.3044e-439 +/- ...e-444]

            sage: epsi = CBF(1e-10)
            sage: CBF.integral(lambda x, _: x*(1/x).sin(), epsi, 1)
            [0.38 +/- ...e-3]
            sage: CBF.integral(lambda x, _: x*(1/x).sin(), epsi, 1, use_heap=True)
            [0.37853002 +/- ...e-9]

        ALGORITHM:

        Uses the `acb_calc <https://flintlib.org/doc/acb_calc.html>`_ module of
        the FLINT library.

        TESTS::

            sage: CBF.integral(lambda x, _: x, 0, 10, rel_tol=1e-10,
            ....:     abs_tol=1e-10, deg_limit=1, eval_limit=20, depth_limit=4,
            ....:     use_heap=False)
            [50.00000000 +/- ...e-9]

            sage: i = QuadraticField(-1).gen()
            sage: CBF.integral(lambda x, _: (1 + i*x).gamma(), -1, 1) # abs tol 1e-13
            [1.57239266949806 +/- 8.33e-15] + [+/- 1.39e-15]*I

            sage: ComplexBallField(10000).integral(lambda x, _: x.sin(), 0, 1, rel_tol=1e-300)
            [0.459... +/- ...e-3...]
            sage: CBF.integral(lambda x, _: x.sin(), 0, 100, rel_tol=10)
            [0.1377 +/- ...e-5]

            sage: ComplexBallField(10000).integral(lambda x, _: x.sin(), 0, 1, abs_tol=1e-400)
            [0.459697... +/- ...e-4...]
            sage: CBF.integral(lambda x, _: x.sin(), 0, 1, abs_tol=10)
            [+/- 0.842]

            sage: ComplexBallField(100).integral(lambda x, _: sin(x), RBF(0), RBF(1))
            [0.4596976941318602825990633926 +/- ...e-29]

            sage: from cysignals.alarm import alarm
            sage: alarm(0.1r)
            sage: C = ComplexBallField(1000000)
            sage: C.integral(lambda x, _: x.cos() * x.sin(), 0, 1)
            Traceback (most recent call last):
            ...
            AlarmInterrupt
        """
        cdef IntegrationContext ctx = IntegrationContext()
        cdef acb_calc_integrate_opt_t arb_opts
        cdef long cgoal, expo
        cdef mag_t ctol
        cdef RealNumber tmp
        cdef ComplexBall ca, cb

        if isinstance(a, (RealBall, ComplexBall)):
            ca = <ComplexBall> self(a)
        else:
            ca = <ComplexBall> self.coerce(a)
        if isinstance(b, (RealBall, ComplexBall)):
            cb = <ComplexBall> self(b)
        else:
            cb = <ComplexBall> self.coerce(b)

        mag_init(ctol)

        ctx.f = func
        ctx.parent = self
        ctx.exn_type = None

        acb_calc_integrate_opt_init(arb_opts)
        if deg_limit is not None:
            arb_opts.deg_limit = deg_limit
        if eval_limit is not None:
            arb_opts.eval_limit = eval_limit
        if depth_limit is not None:
            arb_opts.depth_limit = depth_limit
        if use_heap is not None:
            arb_opts.use_heap = use_heap
        if verbose is not None:
            arb_opts.verbose = verbose

        RR = RealField()
        if rel_tol is None:
            cgoal = self._prec
        else:
            tmp = <RealNumber> RR(rel_tol)
            mpfr_get_d_2exp(&cgoal, tmp.value, MPFR_RNDD)
            cgoal = -cgoal

        if abs_tol is None:
            mag_set_ui_2exp_si(ctol, 1, -self._prec)
        else:
            tmp = <RealNumber> RR(abs_tol)
            mag_set_d(ctol, mpfr_get_d_2exp(&expo, tmp.value, MPFR_RNDD))
            mag_mul_2exp_si(ctol, ctol, expo)

        cdef ComplexBall res = ComplexBall.__new__(ComplexBall)
        res._parent = self

        try:
            sig_on()
            acb_calc_integrate(
                    res.value,
                    <acb_calc_func_t> acb_calc_func_callback,
                    <void *> ctx,
                    ca.value, cb.value,
                    cgoal, ctol, arb_opts, self._prec)
            sig_off()
        finally:
            mag_clear(ctol)

        if ctx.exn_type is not None:
            raise ctx.exn_type, ctx.exn_obj, ctx.exn_tb

        return res


cdef inline bint _do_sig(long prec) noexcept:
    """
    Whether signal handlers should be installed for calls to FLINT.
    """
    return (prec > 1000)

cdef inline long prec(ComplexBall ball) noexcept:
    return ball._parent._prec

cdef bint arb_contained_unit_interval(arb_t b) noexcept:
    r"""
    Test if a real ball is contained in [-1,1]. Useful for dealing with branch
    cuts of inverse trigonometric functions.
    """
    cdef arb_t u
    arb_init(u)
    try:
        arb_one(u)
        if not arb_lt(b, u):
            return False
        arb_neg(u, u)
        if not arb_gt(b, u):
            return False
        return True
    finally:
        arb_clear(u)

cdef bint arb_gt_neg_one(arb_t b) noexcept:
    r"""
    Test if a real ball is contained in [-1,∞). Useful for dealing with branch
    cuts.
    """
    cdef arb_t neg_one
    arb_init(neg_one)
    arb_set_si(neg_one, -1)
    cdef bint res = arb_gt(b, neg_one)
    arb_clear(neg_one)
    return res

cdef inline real_ball_field(ComplexBall ball):
    return ball._parent._base

cdef class ComplexBall(RingElement):
    """
    Hold one ``acb_t`` of the `FLINT library <https://flintlib.org>`_.

    EXAMPLES::

        sage: a = ComplexBallField()(1, 1)
        sage: a
        1.000000000000000 + 1.000000000000000*I
    """
    def __cinit__(self):
        """
        Allocate memory for the encapsulated value.

        EXAMPLES::

            sage: ComplexBallField(2)(0) # indirect doctest
            0
        """
        acb_init(self.value)

    def __dealloc__(self):
        """
        Deallocate memory of the encapsulated value.

        EXAMPLES::

            sage: a = ComplexBallField(2)(0) # indirect doctest
            sage: del a
        """
        acb_clear(self.value)

    def __init__(self, parent, x=None, y=None):
        """
        Initialize the :class:`ComplexBall`.

        INPUT:

        - ``parent`` -- a :class:`ComplexBallField`

        - ``x``, ``y`` -- (optional) either a complex number, interval or ball,
          or two real ones

        .. SEEALSO:: :meth:`ComplexBallField._element_constructor_`

        TESTS::

            sage: from sage.rings.complex_arb import ComplexBall
            sage: CBF53, CBF100 = ComplexBallField(53), ComplexBallField(100)
            sage: ComplexBall(CBF100)
            0
            sage: ComplexBall(CBF100, ComplexBall(CBF53, ComplexBall(CBF100, 1/3)))
            [0.333333333333333333333333333333 +/- ...e-31]
            sage: ComplexBall(CBF100, RBF(pi))                                          # needs sage.symbolic
            [3.141592653589793 +/- ...e-16]

            sage: ComplexBall(CBF100, -3r)
            Traceback (most recent call last):
            ...
            TypeError: unsupported initializer
            sage: CBF100(-3r)
            -3.000000000000000000000000000000

            sage: ComplexBall(CBF100, CIF(1, 2))
            1.000000000000000000000000000000 + 2.000000000000000000000000000000*I
            sage: ComplexBall(CBF100, RBF(1/3), RBF(1))
            [0.3333333333333333 +/- ...e-17] + 1.000000000000000000000000000000*I
            sage: ComplexBall(CBF100, 10^100)
            [1.000000000000000000000000000000e+100 +/- ...]

            sage: NF.<a> = QuadraticField(-1, embedding=CC(0, -1))
            sage: CBF(a)
            -1.000000000000000*I

            sage: NF.<a> = QuadraticField(-1, embedding=None)
            sage: CBF(a)
            1.000000000000000*I
            sage: CBF.coerce(a)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion ...

            sage: NF.<a> = QuadraticField(-2)
            sage: CBF(1/3 + a).real()
            [0.3333333333333333 +/- ...e-17]

            sage: ComplexBall(CBF, 1, 1/2)
            1.000000000000000 + 0.5000000000000000*I
            sage: ComplexBall(CBF, 1, 1)
            1.000000000000000 + 1.000000000000000*I
            sage: ComplexBall(CBF, 1, 1/2)
            1.000000000000000 + 0.5000000000000000*I
            sage: ComplexBall(CBF, 1/2, 1)
            0.5000000000000000 + 1.000000000000000*I
            sage: ComplexBall(CBF, 1/2, 1/2)
            0.5000000000000000 + 0.5000000000000000*I
            sage: ComplexBall(CBF, 1/2, 'a')
            Traceback (most recent call last):
            ...
            TypeError: unsupported initializer
            sage: ComplexBall(CBF, 'a')
            Traceback (most recent call last):
            ...
            TypeError: unsupported initializer
        """
        cdef fmpz_t tmpz
        cdef fmpq_t tmpq
        cdef bint cplx = False

        Element.__init__(self, parent)

        if x is None:
            return
        elif isinstance(x, ComplexBall):
            acb_set(self.value, (<ComplexBall> x).value)
            cplx = True
        elif isinstance(x, RealBall):
            arb_set(acb_realref(self.value), (<RealBall> x).value)
        elif isinstance(x, sage.rings.integer.Integer):
            if _do_sig(prec(self)): sig_on()
            fmpz_init(tmpz)
            fmpz_set_mpz(tmpz, (<sage.rings.integer.Integer> x).value)
            arb_set_fmpz(acb_realref(self.value), tmpz)
            fmpz_clear(tmpz)
            if _do_sig(prec(self)): sig_off()
        elif isinstance(x, sage.rings.rational.Rational):
            if _do_sig(prec(self)): sig_on()
            fmpq_init(tmpq)
            fmpq_set_mpq(tmpq, (<sage.rings.rational.Rational> x).value)
            arb_set_fmpq(acb_realref(self.value), tmpq, prec(self))
            fmpq_clear(tmpq)
            if _do_sig(prec(self)): sig_off()
        elif isinstance(x, ComplexIntervalFieldElement):
            ComplexIntervalFieldElement_to_acb(self.value,
                                               <ComplexIntervalFieldElement> x)
            cplx = True
        else:
            raise TypeError("unsupported initializer")

        if y is None:
            return
        elif cplx:
            raise TypeError("unsupported initializer")
        elif isinstance(y, RealBall):
            arb_set(acb_imagref(self.value), (<RealBall> y).value)
        elif isinstance(y, sage.rings.integer.Integer):
            if _do_sig(prec(self)): sig_on()
            fmpz_init(tmpz)
            fmpz_set_mpz(tmpz, (<sage.rings.integer.Integer> y).value)
            arb_set_fmpz(acb_imagref(self.value), tmpz)
            fmpz_clear(tmpz)
            if _do_sig(prec(self)): sig_off()
        elif isinstance(y, sage.rings.rational.Rational):
            if _do_sig(prec(self)): sig_on()
            fmpq_init(tmpq)
            fmpq_set_mpq(tmpq, (<sage.rings.rational.Rational> y).value)
            arb_set_fmpq(acb_imagref(self.value), tmpq, prec(self))
            fmpq_clear(tmpq)
            if _do_sig(prec(self)): sig_off()
        else:
            raise TypeError("unsupported initializer")

    def __hash__(self):
        """
        TESTS::

            sage: hash(CBF(1/3)) == hash(RBF(1/3))
            True
            sage: hash(CBF(1/3 + 2*i)) != hash(CBF(1/3 + i))
            True
        """
        if self.is_real():
            return hash(self.real())
        else:
            return (hash(self.real()) // 3) ^ hash(self.imag())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        OUTPUT: string

        EXAMPLES::

           sage: CBF(1/3)
           [0.3333333333333333 +/- ...e-17]
           sage: CBF(0, 1/3)
           [0.3333333333333333 +/- ...e-17]*I
           sage: CBF(1/3, 1/6)
           [0.3333333333333333 +/- ...e-17] + [0.1666666666666667 +/- ...e-17]*I

        TESTS::

           sage: CBF(1-I/2)
           1.000000000000000 - 0.5000000000000000*I
        """
        cdef arb_t real = acb_realref(self.value)
        cdef arb_t imag = acb_imagref(self.value)
        if arb_is_zero(imag):
            return self.real()._repr_()
        elif arb_is_zero(real):
            return "{}*I".format(self.imag()._repr_())
        elif arb_is_exact(imag) and arb_is_negative(imag):
            return "{} - {}*I".format(self.real()._repr_(),
                                        (-self.imag())._repr_())
        else:
            return "{} + {}*I".format(self.real()._repr_(),
                                        self.imag()._repr_())

    def __reduce__(self):
        r"""
        Serialize a ComplexBall.

        TESTS::

            sage: [loads(dumps(b)).identical(b) for b in                                # needs sage.symbolic
            ....:     [ComplexBallField(60)(1/3 + i*pi), CBF(NaN)]]
            [True, True]
        """
        return self.__class__, (self._parent, self.real(), self.imag())

    def _is_atomic(self):
        r"""
        Declare that complex balls print atomically in some cases.

        TESTS::

            sage: CBF(-1/3)._is_atomic()
            True

        This method should in principle ensure that ``CBF['x']([1, -1/3])``
        is printed as::

            sage: CBF['x']([1, -1/3]) # todo - not tested
            [-0.3333333333333333 +/- ...e-17]*x + 1.000000000000000

        However, this facility is not really used in Sage at this point, and we
        still get::

            sage: CBF['x']([1, -1/3])
            ([-0.3333333333333333 +/- ...e-17])*x + 1.000000000000000
        """
        return self.is_real() or self.real().is_zero()

    # Conversions

    cpdef ComplexIntervalFieldElement _complex_mpfi_(self, parent):
        """
        Return :class:`ComplexIntervalFieldElement` of the same value.

        EXAMPLES::

            sage: CIF(CBF(1/3, 1/3)) # indirect doctest
            0.3333333333333333? + 0.3333333333333333?*I
        """
        cdef ComplexIntervalFieldElement res = parent.zero()
        res = res._new() # FIXME after modernizing CIF
        acb_to_ComplexIntervalFieldElement(res, self.value)
        return res

    def _integer_(self, _):
        """
        Check that this ball contains a single integer and return that integer.

        EXAMPLES::

            sage: ZZ(CBF(-42, RBF(.1, rad=.2))) # indirect doctest
            -42
            sage: ZZ(CBF(i))
            Traceback (most recent call last):
            ...
            ValueError: 1.000000000000000*I does not contain a unique integer
        """
        cdef Integer res
        cdef fmpz_t tmp
        fmpz_init(tmp)
        try:
            if acb_get_unique_fmpz(tmp, self.value):
                res = PY_NEW(Integer)
                fmpz_get_mpz(res.value, tmp)
            else:
                raise ValueError("{} does not contain a unique integer".format(self))
        finally:
            fmpz_clear(tmp)
        return res

    def _rational_(self):
        """
        Check that this ball contains a single rational number and return that
        number.

        EXAMPLES::

            sage: QQ(CBF(12345/2^5))
            12345/32
            sage: QQ(CBF(i))
            Traceback (most recent call last):
            ...
            ValueError: 1.000000000000000*I does not contain a unique rational number
        """
        if acb_is_real(self.value) and acb_is_exact(self.value):
            return self.real().mid().exact_rational()
        else:
            raise ValueError("{} does not contain a unique rational number".format(self))

    def _complex_mpfr_field_(self, parent):
        r"""
        Convert this complex ball to a complex number.

        INPUT:

        - ``parent`` -- :class:`~sage.rings.complex_mpfr.ComplexField_class`,
          target parent

        EXAMPLES::

            sage: CC(CBF(1/3, 1/3))
            0.333333333333333 + 0.333333333333333*I
            sage: ComplexField(100)(CBF(1/3, 1/3))
            0.33333333333333331482961625625 + 0.33333333333333331482961625625*I

        Check nan and inf::

            sage: CC(CBF('nan', 1/3))
            NaN + 0.333333333333333*I
            sage: CC(CBF('+inf'))
            +infinity
        """
        real_field = parent._base
        return parent(real_field(self.real()), real_field(self.imag()))

    def _real_mpfi_(self, parent):
        r"""
        Try to convert this complex ball to a real interval.

        Fail if the imaginary part is not exactly zero.

        INPUT:

        - ``parent`` -- :class:`~sage.rings.real_mpfi.RealIntervalField_class`,
          target parent

        EXAMPLES::

            sage: RIF(CBF(RBF(1/3, rad=1e-5)))
            0.3334?
            sage: RIF(CBF(RBF(1/3, rad=1e-5), 1e-10))
            Traceback (most recent call last):
            ...
            ValueError: nonzero imaginary part
        """
        if acb_is_real(self.value):
            return parent(self.real())
        else:
            raise ValueError("nonzero imaginary part")

    def _mpfr_(self, parent):
        r"""
        Try to convert this complex ball to a real number.

        Fail if the imaginary part is not exactly zero.

        INPUT:

        - ``parent`` -- :class:`~sage.rings.real_mpfr.RealField_class`,
          target parent

        EXAMPLES::

            sage: RR(CBF(1/3))
            0.333333333333333
            sage: RR(CBF(1, 1/3) - CBF(0, 1/3))
            Traceback (most recent call last):
            ...
            ValueError: nonzero imaginary part
        """
        if acb_is_real(self.value):
            return parent(self.real())
        else:
            raise ValueError("nonzero imaginary part")

    def __float__(self):
        """
        Convert ``self`` to a ``float``.

        EXAMPLES::

            sage: float(CBF(1))
            1.0
            sage: float(CBF(1/3))
            0.3333333333333333
            sage: float(CBF(1,1))
            Traceback (most recent call last):
            ...
            TypeError: can...t convert complex ball to float
        """
        if not arb_is_zero(acb_imagref(self.value)):
            raise TypeError("can't convert complex ball to float")
        return arf_get_d(arb_midref(acb_realref(self.value)), ARF_RND_NEAR)

    def __complex__(self):
        """
        Convert ``self`` to a ``complex``.

        EXAMPLES::

            sage: complex(CBF(1))
            (1+0j)
            sage: complex(CBF(1,1))
            (1+1j)

        Check nan and inf::

            sage: complex(CBF(0, 'nan'))
            nanj
            sage: complex(CBF('+inf', '-inf'))
            (inf-infj)
        """
        return PyComplex_FromDoubles(
                arf_get_d(arb_midref(acb_realref(self.value)), ARF_RND_NEAR),
                arf_get_d(arb_midref(acb_imagref(self.value)), ARF_RND_NEAR))

    def _real_double_(self, parent):
        r"""
        Convert this complex ball to a real double.

        EXAMPLES::

            sage: RDF(CBF(3))
            3.0
            sage: RDF(CBF(3/7))
            0.42857142857142855
            sage: RDF(CBF(1 + I))
            Traceback (most recent call last):
            ...
            TypeError: can...t convert complex ball to float

            sage: RDF(CBF(3)).parent()
            Real Double Field

        Check nan and inf::

            sage: RDF(CBF('nan'))
            NaN
            sage: RDF(CBF('nan')).parent()
            Real Double Field
            sage: RDF(CBF('+inf'))
            +infinity
            sage: RDF(CBF('+inf')).parent()
            Real Double Field

        TESTS:

        Check that conversions go through this method::

            sage: RDF.convert_map_from(CBF)
            Conversion via _real_double_ method map:
              From: Complex ball field with 53 bits of precision
              To:   Real Double Field
        """
        if not arb_is_zero(acb_imagref(self.value)):
            raise TypeError("can't convert complex ball to float")
        cdef RealDoubleElement x
        x = PY_NEW(RealDoubleElement)
        x._value = arf_get_d(arb_midref(acb_realref(self.value)), ARF_RND_NEAR)
        return x

    def _complex_double_(self, parent):
        r"""
        Convert this complex ball to a complex double.

        EXAMPLES::

            sage: CDF(CBF(1+I))
            1.0 + 1.0*I
            sage: CDF(CBF(3/7, -21/13))
            0.42857142857142855 - 1.6153846153846152*I

        Check nan and inf::

            sage: CDF(CBF('nan', 'nan'))
            NaN + NaN*I
            sage: CDF(CBF('+inf', 'nan'))
            +infinity + NaN*I
            sage: CDF(CBF('+inf', '-inf'))
            +infinity - +infinity*I

        TESTS:

        The conversion map should go through this method. However, it is
        manually called in the constructor of complex double elements::

            sage: CDF.convert_map_from(CBF)  # bug
            Conversion map:
              From: Complex ball field with 53 bits of precision
              To:   Complex Double Field
        """
        cdef ComplexDoubleElement x = ComplexDoubleElement.__new__(ComplexDoubleElement)
        x._complex = gsl_complex_rect(
                arf_get_d(arb_midref(acb_realref(self.value)), ARF_RND_NEAR),
                arf_get_d(arb_midref(acb_imagref(self.value)), ARF_RND_NEAR))
        return x

    # Real and imaginary part, midpoint, radius

    cpdef RealBall real(self):
        """
        Return the real part of this ball.

        OUTPUT: a :class:`~sage.rings.real_arb.RealBall`

        EXAMPLES::

           sage: a = CBF(1/3, 1/5)
           sage: a.real()
           [0.3333333333333333 +/- ...e-17]
           sage: a.real().parent()
           Real ball field with 53 bits of precision
        """
        cdef RealBall r = RealBall.__new__(RealBall)
        r._parent = real_ball_field(self)
        arb_set(r.value, acb_realref(self.value))
        return r

    cpdef RealBall imag(self):
        """
        Return the imaginary part of this ball.

        OUTPUT: a :class:`~sage.rings.real_arb.RealBall`

        EXAMPLES::

           sage: a = CBF(1/3, 1/5)
           sage: a.imag()
           [0.2000000000000000 +/- ...e-17]
           sage: a.imag().parent()
           Real ball field with 53 bits of precision
        """
        cdef RealBall r = RealBall.__new__(RealBall)
        r._parent = real_ball_field(self)
        arb_set(r.value, acb_imagref(self.value))
        return r

    def __abs__(self):
        """
        Return the absolute value of this complex ball.

        EXAMPLES::

            sage: CBF(1 + i).abs() # indirect doctest
            [1.414213562373095 +/- ...e-16]
            sage: abs(CBF(i))
            1.000000000000000

            sage: CBF(1 + i).abs().parent()
            Real ball field with 53 bits of precision
        """
        cdef RealBall r = RealBall(real_ball_field(self))
        acb_abs(r.value, self.value, prec(self))
        return r

    def below_abs(self, test_zero=False):
        """
        Return a lower bound for the absolute value of this complex ball.

        INPUT:

        - ``test_zero`` -- boolean (default: ``False``); if ``True``,
          make sure that the returned lower bound is positive, raising
          an error if the ball contains zero.

        OUTPUT: a ball with zero radius

        EXAMPLES::

            sage: b = ComplexBallField(8)(1+i).below_abs()
            sage: b
            [1.4 +/- 0.0141]
            sage: b.is_exact()
            True
            sage: QQ(b)*128
            181
            sage: (CBF(1/3) - 1/3).below_abs()
            0
            sage: (CBF(1/3) - 1/3).below_abs(test_zero=True)
            Traceback (most recent call last):
            ...
            ValueError: ball contains zero

        .. SEEALSO:: :meth:`above_abs`
        """
        cdef RealBall res = RealBall(real_ball_field(self))
        acb_get_abs_lbound_arf(arb_midref(res.value), self.value, prec(self))
        if test_zero and arb_contains_zero(res.value):
            assert acb_contains_zero(self.value)
            raise ValueError("ball contains zero")
        return res

    def above_abs(self):
        """
        Return an upper bound for the absolute value of this complex ball.

        OUTPUT: a ball with zero radius

        EXAMPLES::

            sage: b = ComplexBallField(8)(1+i).above_abs()
            sage: b
            [1.4 +/- 0.0219]
            sage: b.is_exact()
            True
            sage: QQ(b)*128
            182

        .. SEEALSO:: :meth:`below_abs`
        """
        cdef RealBall res = RealBall(real_ball_field(self))
        acb_get_abs_ubound_arf(arb_midref(res.value), self.value, prec(self))
        return res

    def arg(self):
        """
        Return the argument of this complex ball.

        EXAMPLES::

            sage: CBF(1 + i).arg()
            [0.7853981633974483 +/- ...e-17]
            sage: CBF(-1).arg()
            [3.141592653589793 +/- ...e-16]
            sage: CBF(-1).arg().parent()
            Real ball field with 53 bits of precision
        """
        cdef RealBall r = RealBall(real_ball_field(self))
        acb_arg(r.value, self.value, prec(self))
        return r

    def mid(self):
        """
        Return the midpoint of this ball.

        OUTPUT:

        :class:`~sage.rings.complex_mpfr.ComplexNumber`, floating-point
        complex number formed by the centers of the real and imaginary parts of
        this ball.

        EXAMPLES::

            sage: CBF(1/3, 1).mid()
            0.333333333333333 + 1.00000000000000*I
            sage: CBF(1/3, 1).mid().parent()
            Complex Field with 53 bits of precision
            sage: CBF('inf', 'nan').mid()
            +infinity + NaN*I
            sage: CBF('nan', 'inf').mid()
            NaN + +infinity*I
            sage: CBF('nan').mid()
            NaN
            sage: CBF('inf').mid()
            +infinity
            sage: CBF(0, 'inf').mid()
            +infinity*I

        .. SEEALSO:: :meth:`squash`
        """
        re, im = self.real().mid(), self.imag().mid()
        field = sage.rings.complex_mpfr.ComplexField(
                max(prec(self), re.prec(), im.prec()))
        return field(re, im)

    def squash(self):
        """
        Return an exact ball with the same midpoint as this ball.

        OUTPUT: a :class:`ComplexBall`

        EXAMPLES::

            sage: mid = CBF(1/3, 1/10).squash()
            sage: mid
            [0.3333333333333333 +/- ...e-17] + [0.09999999999999999 +/- ...e-18]*I
            sage: mid.parent()
            Complex ball field with 53 bits of precision
            sage: mid.is_exact()
            True

        .. SEEALSO:: :meth:`mid`
        """
        cdef ComplexBall res = self._new()
        arf_set(arb_midref(acb_realref(res.value)), arb_midref(acb_realref(self.value)))
        arf_set(arb_midref(acb_imagref(res.value)), arb_midref(acb_imagref(self.value)))
        mag_zero(arb_radref(acb_realref(res.value)))
        mag_zero(arb_radref(acb_imagref(res.value)))
        return res

    def rad(self):
        """
        Return an upper bound for the error radius of this ball.

        OUTPUT:

        A :class:`~sage.rings.real_mpfr.RealNumber` of the same precision as
        the radii of real balls.

        .. WARNING::

            Unlike a :class:`~sage.rings.real_arb.RealBall`,
            a :class:`ComplexBall` is *not* defined
            by its midpoint and radius. (Instances of :class:`ComplexBall` are
            actually rectangles, not balls.)

        EXAMPLES::

            sage: CBF(1 + i).rad()
            0.00000000
            sage: CBF(i/3).rad()
            1.1102230e-16
            sage: CBF(i/3).rad().parent()
            Real Field with 30 bits of precision

        .. SEEALSO:: :meth:`diameter`, :meth:`mid`

        TESTS::

            sage: (CBF(0, 1/3) << (2^64)).rad()
            Traceback (most recent call last):
            ...
            RuntimeError: unable to convert the radius to MPFR (exponent out of range?)
        """
        # Should we return a real number with rounding towards +∞ (or away from
        # zero if/when implemented)?
        cdef RealField_class rad_field = RealField(MAG_BITS)
        cdef RealNumber rad = RealNumber(rad_field, None)
        cdef arf_t tmp
        arf_init(tmp)
        acb_get_rad_ubound_arf(tmp, self.value, MAG_BITS)
        sig_str("unable to convert the radius to MPFR (exponent out of range?)")
        if arf_get_mpfr(rad.value, tmp, MPFR_RNDU):
            sig_error()
        sig_off()
        arf_clear(tmp)
        return rad

    # Should we implement rad_as_ball? If we do, should it return an enclosure
    # of the radius (which radius?), or an upper bound?

    def diameter(self):
        r"""
        Return the diameter of this ball.

        EXAMPLES::

            sage: CBF(1 + i).diameter()
            0.00000000
            sage: CBF(i/3).diameter()
            2.2204460e-16
            sage: CBF(i/3).diameter().parent()
            Real Field with 30 bits of precision
            sage: CBF(CIF(RIF(1.02, 1.04), RIF(2.1, 2.2))).diameter()
            0.20000000

        .. SEEALSO:: :meth:`rad`, :meth:`mid`
        """
        return 2 * self.rad()

    def union(self, other):
        r"""
        Return a ball containing the convex hull of ``self`` and ``other``.

        EXAMPLES::

            sage: b = CBF(1 + i).union(0)
            sage: b.real().endpoints()
            (-9.31322574615479e-10, 1.00000000093133)
        """
        cdef ComplexBall my_other = self._parent.coerce(other)
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_union(res.value, self.value, my_other.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    # Precision and accuracy

    def nbits(self):
        r"""
        Return the minimum precision sufficient to represent this ball exactly.

        More precisely, the output is the number of bits needed to represent
        the absolute value of the mantissa of both the real and the imaginary
        part of the midpoint.

        EXAMPLES::

            sage: CBF(17, 1023).nbits()
            10
            sage: CBF(1/3, NaN).nbits()                                                 # needs sage.symbolic
            53
            sage: CBF(NaN).nbits()                                                      # needs sage.symbolic
            0
        """
        return acb_bits(self.value)

    def round(self):
        """
        Return a copy of this ball rounded to the precision of the parent.

        EXAMPLES:

        It is possible to create balls whose midpoint is more precise that
        their parent's nominal precision (see :mod:`~sage.rings.real_arb` for
        more information)::

            sage: b = CBF(exp(I*pi/3).n(100))                                           # needs sage.symbolic
            sage: b.mid()                                                               # needs sage.symbolic
            0.50000000000000000000000000000 + 0.86602540378443864676372317075*I

        The ``round()`` method rounds such a ball to its parent's precision::

            sage: b.round().mid()                                                       # needs sage.symbolic
            0.500000000000000 + 0.866025403784439*I

        .. SEEALSO:: :meth:`trim`
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_set_round(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def accuracy(self):
        """
        Return the effective relative accuracy of this ball measured in bits.

        This is computed as if calling
        :meth:`~sage.rings.real_arb.RealBall.accuracy()`
        on the real ball whose midpoint is the larger out of the real and
        imaginary midpoints of this complex ball, and whose radius is the
        larger out of the real and imaginary radii of this complex ball.

        EXAMPLES::

            sage: CBF(exp(I*pi/3)).accuracy()                                           # needs sage.symbolic
            51
            sage: CBF(I/2).accuracy() == CBF.base().maximal_accuracy()
            True
            sage: CBF('nan', 'inf').accuracy() == -CBF.base().maximal_accuracy()
            True

        .. SEEALSO::

            :meth:`~sage.rings.real_arb.RealBallField.maximal_accuracy`
        """
        return acb_rel_accuracy_bits(self.value)

    def trim(self):
        """
        Return a trimmed copy of this ball.

        Return a copy of this ball with both the real and imaginary parts
        trimmed (see :meth:`~sage.rings.real_arb.RealBall.trim()`).

        EXAMPLES::

            sage: b = CBF(1/3, RBF(1/3, rad=.01))
            sage: b.mid()
            0.333333333333333 + 0.333333333333333*I
            sage: b.trim().mid()
            0.333333333333333 + 0.333333015441895*I

        .. SEEALSO:: :meth:`round`
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_trim(res.value, self.value)
        if _do_sig(prec(self)): sig_off()
        return res

    def add_error(self, ampl):
        """
        Increase the radii of the real and imaginary parts by (an upper bound
        on) ``ampl``.

        If ``ampl`` is negative, the radii remain unchanged.

        INPUT:

        - ``ampl`` -- a **real** ball (or an object that can be coerced to a
          real ball)

        OUTPUT: a new complex ball

        EXAMPLES::

            sage: CBF(1+i).add_error(10^-16)
            [1.000000000000000 +/- ...e-16] + [1.000000000000000 +/- ...e-16]*I
        """
        return ComplexBall(self._parent, self.real().add_error(ampl), self.imag().add_error(ampl))

    # Comparisons and predicates

    def is_NaN(self):
        """
        Return ``True`` iff either the real or the imaginary part
        is not-a-number.

        EXAMPLES::

            sage: CBF(NaN).is_NaN()                                                     # needs sage.symbolic
            True
            sage: CBF(-5).gamma().is_NaN()
            True
            sage: CBF(oo).is_NaN()
            False
            sage: CBF(42+I).is_NaN()
            False
        """
        return (arf_is_nan(arb_midref(acb_realref(self.value)))
                or arf_is_nan(arb_midref(acb_imagref(self.value))))

    def is_zero(self):
        """
        Return ``True`` iff the midpoint and radius of this ball are both zero.

        EXAMPLES::

            sage: CBF(0).is_zero()
            True
            sage: CBF(RIF(-0.5, 0.5)).is_zero()
            False

        .. SEEALSO:: :meth:`is_nonzero`
        """
        return acb_is_zero(self.value)

    def is_nonzero(self):
        """
        Return ``True`` iff zero is not contained in the interval represented
        by this ball.

        .. NOTE::

            This method is not the negation of :meth:`is_zero`: it only
            returns ``True`` if zero is known not to be contained in the ball.

            Use ``bool(b)`` (or, equivalently, ``not b.is_zero()``) to check if
            a ball ``b`` **may** represent a nonzero number (for instance, to
            determine the “degree” of a polynomial with ball coefficients).

        EXAMPLES::

            sage: CBF(pi, 1/3).is_nonzero()                                             # needs sage.symbolic
            True
            sage: CBF(RIF(-0.5, 0.5), 1/3).is_nonzero()
            True
            sage: CBF(1/3, RIF(-0.5, 0.5)).is_nonzero()
            True
            sage: CBF(RIF(-0.5, 0.5), RIF(-0.5, 0.5)).is_nonzero()
            False

        .. SEEALSO:: :meth:`is_zero`
        """
        return (arb_is_nonzero(acb_realref(self.value))
                or arb_is_nonzero(acb_imagref(self.value)))

    def __bool__(self):
        """
        Return ``True`` iff this complex ball is not the zero ball, i.e. if the
        midpoint and radius of its real and imaginary parts are not all zero.

        This is the preferred way, for instance, to determine the “degree” of a
        polynomial with ball coefficients.

        .. WARNING::

            A “nonzero” ball in the sense of this method may represent the
            value zero. Use :meth:`is_nonzero` to check that a complex number
            represented by a ``ComplexBall`` object is known to be nonzero.

        EXAMPLES::

            sage: bool(CBF(0)) # indirect doctest
            False
            sage: bool(CBF(i))
            True
            sage: bool(CBF(RIF(-0.5, 0.5)))
            True
        """
        return not acb_is_zero(self.value)

    def is_exact(self):
        """
        Return ``True`` iff the radius of this ball is zero.

        EXAMPLES::

            sage: CBF(1).is_exact()
            True
            sage: CBF(1/3, 1/3).is_exact()
            False
        """
        return acb_is_exact(self.value)

    def is_real(self):
        """
        Return ``True`` iff the imaginary part of this ball is exactly zero.

        EXAMPLES::

            sage: CBF(1/3, 0).is_real()
            True
            sage: (CBF(i/3) - CBF(1, 1/3)).is_real()
            False
            sage: CBF('inf').is_real()
            True
        """
        return acb_is_real(self.value)

    cpdef _richcmp_(left, right, int op):
        """
        Compare ``left`` and ``right``.

        For more information, see :mod:`sage.rings.complex_arb`.

        EXAMPLES::

            sage: a = CBF(1)
            sage: b = CBF(1)
            sage: a is b
            False
            sage: a == b
            True
            sage: a = CBF(1/3)
            sage: a.is_exact()
            False
            sage: b = CBF(1/3)
            sage: b.is_exact()
            False
            sage: a == b
            False
            sage: a = CBF(1, 2)
            sage: b = CBF(1, 2)
            sage: a is b
            False
            sage: a == b
            True

        TESTS:

        Balls whose intersection consists of one point::

            sage: a = CBF(RIF(1, 2), RIF(1, 2))
            sage: b = CBF(RIF(2, 4), RIF(2, 4))
            sage: a < b
            Traceback (most recent call last):
            ...
            TypeError: No order is defined for ComplexBalls.
            sage: a > b
            Traceback (most recent call last):
            ...
            TypeError: No order is defined for ComplexBalls.
            sage: a <= b
            Traceback (most recent call last):
            ...
            TypeError: No order is defined for ComplexBalls.
            sage: a >= b
            Traceback (most recent call last):
            ...
            TypeError: No order is defined for ComplexBalls.
            sage: a == b
            False
            sage: a != b
            False

        Balls with non-trivial intersection::

            sage: a = CBF(RIF(1, 4), RIF(1, 4))
            sage: a = CBF(RIF(2, 5), RIF(2, 5))
            sage: a == b
            False
            sage: a != b
            False

        One ball contained in another::

            sage: a = CBF(RIF(1, 4), RIF(1, 4))
            sage: b = CBF(RIF(2, 3), RIF(2, 3))
            sage: a == b
            False
            sage: a != b
            False

        Disjoint balls::

            sage: a = CBF(1/3, 1/3)
            sage: b = CBF(1/5, 1/5)
            sage: a == b
            False
            sage: a != b
            True

        Exact elements::

            sage: a = CBF(2, 2)
            sage: b = CBF(2, 2)
            sage: a.is_exact()
            True
            sage: b.is_exact()
            True
            sage: a == b
            True
            sage: a != b
            False
        """
        cdef ComplexBall lt, rt

        lt = left
        rt = right

        if op == Py_EQ:
            return acb_eq(lt.value, rt.value)
        if op == Py_NE:
            return acb_ne(lt.value, rt.value)
        raise TypeError("No order is defined for ComplexBalls.")

    def identical(self, ComplexBall other):
        """
        Return whether ``self`` and ``other`` represent the same ball.

        INPUT:

        - ``other`` -- a :class:`ComplexBall`

        OUTPUT:

        Return ``True`` iff ``self`` and ``other`` are equal as sets, i.e. if their
        real and imaginary parts each have the same midpoint and radius.

        Note that this is not the same thing as testing whether both ``self``
        and ``other`` certainly represent the complex real number, unless
        either ``self`` or ``other`` is exact (and neither contains NaN). To
        test whether both operands might represent the same mathematical
        quantity, use :meth:`overlaps` or ``in``, depending on the
        circumstance.

        EXAMPLES::

            sage: CBF(1, 1/3).identical(1 + CBF(0, 1)/3)
            True
            sage: CBF(1, 1).identical(1 + CBF(0, 1/3)*3)
            False
        """
        return acb_equal(self.value, other.value)

    def overlaps(self, ComplexBall other):
        """
        Return ``True`` iff ``self`` and ``other`` have some point in common.

        INPUT:

        - ``other`` -- a :class:`ComplexBall`

        EXAMPLES::

            sage: CBF(1, 1).overlaps(1 + CBF(0, 1/3)*3)
            True
            sage: CBF(1, 1).overlaps(CBF(1, 'nan'))
            True
            sage: CBF(1, 1).overlaps(CBF(0, 'nan'))
            False
        """
        return acb_overlaps(self.value, other.value)

    def contains_exact(self, other):
        """
        Return ``True`` *iff* ``other`` is contained in ``self``.

        Use ``other in self`` for a test that works for a wider range of inputs
        but may return false negatives.

        INPUT:

        - ``other`` -- :class:`ComplexBall`,
          :class:`~sage.rings.integer.Integer`,
          or :class:`~sage.rings.rational.Rational`

        EXAMPLES::

            sage: CBF(RealBallField(100)(1/3), 0).contains_exact(1/3)
            True
            sage: CBF(1).contains_exact(1)
            True
            sage: CBF(1).contains_exact(CBF(1))
            True

            sage: CBF(sqrt(2)).contains_exact(sqrt(2))                                  # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: unsupported type: <class 'sage.symbolic.expression.Expression'>
        """
        cdef fmpz_t tmpz
        cdef fmpq_t tmpq
        if _do_sig(prec(self)): sig_on()
        try:
            if isinstance(other, ComplexBall):
                res = acb_contains(self.value, (<ComplexBall> other).value)
            elif isinstance(other, Integer):
                fmpz_init(tmpz)
                fmpz_set_mpz(tmpz, (<Integer> other).value)
                res = acb_contains_fmpz(self.value, tmpz)
                fmpz_clear(tmpz)
            elif isinstance(other, sage.rings.rational.Rational):
                fmpq_init(tmpq)
                fmpq_set_mpq(tmpq, (<sage.rings.rational.Rational> other).value)
                res = acb_contains_fmpq(self.value, tmpq)
                fmpq_clear(tmpq)
            else:
                raise TypeError("unsupported type: " + str(type(other)))
        finally:
            if _do_sig(prec(self)): sig_off()
        return res

    def __contains__(self, other):
        """
        Return ``True`` if ``other`` can be verified to be contained in ``self``.

        Depending on the type of ``other``, the test may use interval
        arithmetic with a precision determined by the parent of ``self`` and
        may return false negatives.

        EXAMPLES::

            sage: 1/3*i in CBF(0, 1/3)
            True

        A false negative::

            sage: RLF(1/3) in CBF(RealBallField(100)(1/3), 0)
            False

        .. SEEALSO:: :meth:`contains_exact`
        """
        if not isinstance(other, (
                ComplexBall,
                Integer,
                sage.rings.rational.Rational)):
            other = self._parent(other)
        return self.contains_exact(other)

    def contains_zero(self):
        """
        Return ``True`` iff this ball contains zero.

        EXAMPLES::

            sage: CBF(0).contains_zero()
            True
            sage: CBF(RIF(-1,1)).contains_zero()
            True
            sage: CBF(i).contains_zero()
            False
        """
        return acb_contains_zero(self.value)

    def contains_integer(self):
        """
        Return ``True`` iff this ball contains any integer.

        EXAMPLES::

            sage: CBF(3, RBF(0.1)).contains_integer()
            False
            sage: CBF(3, RBF(0.1,0.1)).contains_integer()
            True
        """
        return acb_contains_int(self.value)

    # Arithmetic

    def __neg__(self):
        """
        Return the opposite of this ball.

        EXAMPLES::

            sage: -CBF(1/3 + I)
            [-0.3333333333333333 +/- ...e-17] - 1.000000000000000*I
        """
        cdef ComplexBall res = self._new()
        acb_neg(res.value, self.value)
        return res

    def conjugate(self):
        """
        Return the complex conjugate of this ball.

        EXAMPLES::

            sage: CBF(-2 + I/3).conjugate()
            -2.000000000000000 + [-0.3333333333333333 +/- ...e-17]*I
        """
        cdef ComplexBall res = self._new()
        acb_conj(res.value, self.value)
        return res

    cpdef _add_(self, other):
        """
        Return the sum of two balls, rounded to the ambient field's precision.

        The resulting ball is guaranteed to contain the sums of any two points
        of the respective input balls.

        EXAMPLES::

            sage: CBF(1) + CBF(I)
            1.000000000000000 + 1.000000000000000*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_add(res.value, self.value, (<ComplexBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cpdef _sub_(self, other):
        """
        Return the difference of two balls, rounded to the ambient field's
        precision.

        The resulting ball is guaranteed to contain the differences of any two
        points of the respective input balls.

        EXAMPLES::

            sage: CBF(1) - CBF(I)
            1.000000000000000 - 1.000000000000000*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_sub(res.value, self.value, (<ComplexBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def __invert__(self):
        """
        Return the inverse of this ball.

        The result is guaranteed to contain the inverse of any point of the
        input ball.

        EXAMPLES::

            sage: ~CBF(i/3)
            [-3.00000000000000 +/- ...e-16]*I
            sage: ~CBF(0)
            nan
            sage: ~CBF(RIF(10,11))
            [0.1 +/- ...e-3]
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_inv(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cpdef _mul_(self, other):
        """
        Return the product of two balls, rounded to the ambient field's
        precision.

        The resulting ball is guaranteed to contain the products of any two
        points of the respective input balls.

        EXAMPLES::

            sage: CBF(-2, 1)*CBF(1, 1/3)
            [-2.333333333333333 +/- ...e-16] + [0.333333333333333 +/- ...e-16]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_mul(res.value, self.value, (<ComplexBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def __lshift__(val, shift):
        r"""
        If ``val`` is a ``ComplexBall`` and ``shift`` is an integer, return the
        ball obtained by shifting the center and radius of ``val`` to the left
        by ``shift`` bits.

        INPUT:

        - ``shift`` -- integer; may be negative

        EXAMPLES::

            sage: CBF(i/3) << 2
            [1.333333333333333 +/- ...e-16]*I
            sage: CBF(i) << -2
            0.2500000000000000*I

        TESTS::

            sage: CBF(i) << (2^65)
            [3.636549880934858e+11106046577046714264 +/- ...e+11106046577046714248]*I
            sage: 'a' << CBF(1/3)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for <<: 'str' and
            'ComplexBall'
            sage: CBF(1) << 1/2
            Traceback (most recent call last):
            ...
            TypeError: shift should be an integer
        """
        cdef fmpz_t tmpz
        # the ComplexBall might be shift, not val
        if not isinstance(val, ComplexBall):
            raise TypeError("unsupported operand type(s) for <<: '{}' and '{}'"
                            .format(type(val).__name__, type(shift).__name__))
        cdef ComplexBall self = val
        cdef ComplexBall res = self._new()
        if is_small_python_int(shift):
            acb_mul_2exp_si(res.value, self.value, PyLong_AsLong(shift))
        elif isinstance(shift, Integer):
            sig_on()
            fmpz_init(tmpz)
            fmpz_set_mpz(tmpz, (<Integer> shift).value)
            acb_mul_2exp_fmpz(res.value, self.value, tmpz)
            fmpz_clear(tmpz)
            sig_off()
        else:
            raise TypeError("shift should be an integer")
        return res

    def __rshift__(val, shift):
        r"""
        If ``val`` is a ``ComplexBall`` and ``shift`` is an integer, return the
        ball obtained by shifting the center and radius of ``val`` to the right
        by ``shift`` bits.

        INPUT:

        - ``shift`` -- integer; may be negative

        EXAMPLES::

            sage: CBF(1+I) >> 2
            0.2500000000000000 + 0.2500000000000000*I

        TESTS::

            sage: 'a' >> CBF(1/3)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for >>: 'str' and
            'ComplexBall'
        """
        # the ComplexBall might be shift, not val
        if isinstance(val, ComplexBall):
            return val << (-shift)
        else:
            raise TypeError("unsupported operand type(s) for >>: '{}' and '{}'"
                            .format(type(val).__name__, type(shift).__name__))

    cpdef _div_(self, other):
        """
        Return the quotient of two balls, rounded to the ambient field's
        precision.

        The resulting ball is guaranteed to contain the quotients of any two
        points of the respective input balls.

        EXAMPLES::

            sage: CBF(-2, 1)/CBF(1, 1/3)
            [-1.500000000000000 +/- ...e-16] + [1.500000000000000 +/- ...e-16]*I
            sage: CBF(2+I)/CBF(0)
            nan + nan*I
            sage: CBF(1)/CBF(0)
            nan
            sage: CBF(1)/CBF(RBF(0, 1.))
            nan
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_div(res.value, self.value, (<ComplexBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def __pow__(base, expo, _):
        """
        EXAMPLES::

            sage: CBF(-1)**(1/2)
            1.000000000000000*I
            sage: CBF(e)**CBF(i*pi)                                                     # needs sage.symbolic
            [-1.00000000000000 +/- ...e-16] + [+/- ...e-15]*I
            sage: CBF(0, 1)**AA(2)**(1/2)
            [-0.60569986707881 +/- ...e-15] + [0.79569320156748 +/- ...e-15]*I

            sage: CBF(i)**RBF(2**1000)
            [+/- 1.01] + [+/- 1.01]*I
            sage: CBF(i)**(2**1000)
            1.000000000000000

            sage: CBF(0)^(1/3)
            0
            sage: CBF(0)^(-1)
            nan
            sage: CBF(0)^(-2)
            nan + nan*I

        TESTS::

            sage: (CBF(e)**CBF(i))**RBF(pi)                                             # needs sage.symbolic
            [-1.0000000000000 +/- ...e-15] + [+/- ...e-15]*I
            sage: CBF(2*i)**10r
            -1024.000000000000
            sage: CBF(1,1) ^ -1r
            0.5000000000000000 - 0.5000000000000000*I
            sage: CBF(2)**SR.var('x')                                                   # needs sage.symbolic
            2.000000000000000^x
        """
        if (isinstance(base, ComplexBall)
                # explicit whitelist due to difference in semantics:
                # ball**non_ball may need to coerce both its arguments
                and isinstance(expo, (int, Integer, RealBall, ComplexBall))):
            return (<ComplexBall> base).pow(expo)
        else:
            return sage.structure.element.bin_op(base, expo, operator.pow)

    cpdef pow(self, expo, analytic=False):
        r"""
        Raise this ball to the power of ``expo``.

        INPUT:

        - ``analytic`` -- boolean (default: ``False``); if ``True``, return an
          indeterminate (not-a-number) value when the exponent is not an
          integer and the base ball touches the branch cut of the logarithm

        EXAMPLES::

            sage: CBF(-1).pow(CBF(i))
            [0.0432139182637723 +/- ...e-17]
            sage: CBF(-1).pow(CBF(i), analytic=True)
            nan + nan*I
            sage: CBF(-10).pow(-2)
            [0.0100000000000000 +/- ...e-18]
            sage: CBF(-10).pow(-2, analytic=True)
            [0.0100000000000000 +/- ...e-18]

        TESTS::

            sage: CBF(2).pow(SR.var('x'))                                               # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Symbolic Ring to Complex ball
            field with 53 bits of precision
        """
        cdef fmpz_t tmpz
        cdef ComplexBall res = self._new()
        if is_small_python_int(expo):
            if _do_sig(prec(self)): sig_on()
            acb_pow_si(res.value, self.value, PyLong_AsLong(expo), prec(self))
            if _do_sig(prec(self)): sig_off()
        elif isinstance(expo, Integer):
            if _do_sig(prec(self)): sig_on()
            fmpz_init(tmpz)
            fmpz_set_mpz(tmpz, (<Integer> expo).value)
            acb_pow_fmpz(res.value, self.value, tmpz, prec(self))
            fmpz_clear(tmpz)
            if _do_sig(prec(self)): sig_off()
        elif isinstance(expo, RealBall):
            if (analytic and not arb_is_int((<RealBall> expo).value)
                    and arb_contains_zero(acb_imagref(self.value))
                    and arb_contains_nonpositive(acb_realref(self.value))):
                acb_indeterminate(res.value)
            else:
                if _do_sig(prec(self)): sig_on()
                acb_pow_arb(res.value, self.value, (<RealBall> expo).value, prec(self))
                if _do_sig(prec(self)): sig_off()
        else:
            if not isinstance(expo, ComplexBall):
                expo = self._parent.coerce(expo)
            if (analytic and not acb_is_int((<ComplexBall> expo).value)
                    and arb_contains_zero(acb_imagref(self.value))
                    and arb_contains_nonpositive(acb_realref(self.value))):
                acb_indeterminate(res.value)
            else:
                if _do_sig(prec(self)): sig_on()
                acb_pow(res.value, self.value, (<ComplexBall> expo).value, prec(self))
                if _do_sig(prec(self)): sig_off()
        return res

    def sqrt(self, analytic=False):
        """
        Return the square root of this ball.

        If either the real or imaginary part is exactly zero, only a single
        real square root is needed.

        INPUT:

        - ``analytic`` -- boolean (default: ``False``); if ``True``, return an
          indeterminate (not-a-number) value when the input ball touches
          the branch cut

        EXAMPLES::

            sage: CBF(-2).sqrt()
            [1.414213562373095 +/- ...e-16]*I
            sage: CBF(-2).sqrt(analytic=True)
            nan + nan*I
        """
        cdef ComplexBall res = self._new()
        if (analytic and arb_contains_zero(acb_imagref(self.value))
                     and arb_contains_nonpositive(acb_realref(self.value))):
            acb_indeterminate(res.value)
        else:
            if _do_sig(prec(self)): sig_on()
            acb_sqrt(res.value, self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        return res

    def rsqrt(self, analytic=False):
        """
        Return the reciprocal square root of ``self``.

        If either the real or imaginary part is exactly zero, only a single
        real reciprocal square root is needed.

        INPUT:

        - ``analytic`` -- boolean (default: ``False``); if ``True``, return an
          indeterminate (not-a-number) value when the input ball touches
          the branch cut

        EXAMPLES::

            sage: CBF(-2).rsqrt()
            [-0.707106781186547 +/- ...e-16]*I
            sage: CBF(-2).rsqrt(analytic=True)
            nan + nan*I
            sage: CBF(0, 1/2).rsqrt()
            1.000000000000000 - 1.000000000000000*I
            sage: CBF(0).rsqrt()
            nan + nan*I
        """
        cdef ComplexBall res = self._new()
        if (analytic and arb_contains_zero(acb_imagref(self.value))
                     and arb_contains_nonpositive(acb_realref(self.value))):
            acb_indeterminate(res.value)
        else:
            if _do_sig(prec(self)): sig_on()
            acb_rsqrt(res.value, self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        return res

    def cube(self):
        """
        Return the cube of this ball.

        The result is computed efficiently using two real squarings, two real
        multiplications, and scalar operations.

        EXAMPLES::

            sage: CBF(1, 1).cube()
            -2.000000000000000 + 2.000000000000000*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_cube(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def rising_factorial(self, n):
        r"""
        Return the ``n``-th rising factorial of this ball.

        The `n`-th rising factorial of `x` is equal to `x (x+1) \cdots (x+n-1)`.

        For complex `n`, it is a quotient of gamma functions.

        EXAMPLES::

            sage: CBF(1).rising_factorial(5)
            120.0000000000000
            sage: CBF(1/3, 1/2).rising_factorial(300)
            [-3.87949484514e+612 +/- 5...e+600] + [-3.52042209763e+612 +/- 5...e+600]*I

            sage: CBF(1).rising_factorial(-1)
            nan
            sage: CBF(1).rising_factorial(2**64)
            [+/- ...e+347382171326740403407]
            sage: ComplexBallField(128)(1).rising_factorial(2**64)
            [2.34369112679686134...e+347382171305201285713 +/- ...]
            sage: CBF(1/2).rising_factorial(CBF(2,3)) # abs tol 1e-15
            [-0.123060451458124 +/- 3.06e-16] + [0.0406412631676552 +/- 7.57e-17]*I
        """
        cdef ComplexBall result = self._new()
        cdef ComplexBall my_n = self._parent.coerce(n)
        if _do_sig(prec(self)): sig_on()
        acb_rising(result.value, self.value, my_n.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    # Elementary functions

    def log(self, base=None, analytic=False):
        """
        General logarithm (principal branch).

        INPUT:

        - ``base`` -- (optional) complex ball or number; if ``None``, return
          the principal branch of the natural logarithm ``ln(self)``,
          otherwise, return the general logarithm ``ln(self)/ln(base)``

        - ``analytic`` -- boolean (default: ``False``); if ``True``, return an
          indeterminate (not-a-number) value when the input ball touches
          the branch cut (with respect to ``self``)

        EXAMPLES::

            sage: CBF(2*i).log()
            [0.693147180559945 +/- ...e-16] + [1.570796326794897 +/- ...e-16]*I
            sage: CBF(-1).log()
            [3.141592653589793 +/- ...e-16]*I

            sage: CBF(2*i).log(2)
            [1.000000000000000 +/- ...e-16] + [2.26618007091360 +/- ...e-15]*I
            sage: CBF(2*i).log(CBF(i))
            [1.000000000000000 +/- ...e-16] + [-0.441271200305303 +/- ...e-16]*I

            sage: CBF('inf').log()
            [+/- inf]
            sage: CBF(2).log(0)
            nan + nan*I

            sage: CBF(-1).log(2)
            [4.53236014182719 +/- ...e-15]*I
            sage: CBF(-1).log(2, analytic=True)
            nan + nan*I
            sage: CBF(-1, RBF(0, rad=.1r)).log(analytic=False)
            [+/- ...e-3] + [+/- 3.15]*I
        """
        cdef ComplexBall cst
        cdef ComplexBall res = self._new()
        if (analytic and arb_contains_zero(acb_imagref(self.value))
                     and arb_contains_nonpositive(acb_realref(self.value))):
            acb_indeterminate(res.value)
        else:
            if _do_sig(prec(self)): sig_on()
            acb_log(res.value, self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
            if base is not None:
                cst = self._parent.coerce(base).log()
                if _do_sig(prec(self)): sig_on()
                acb_div(res.value, res.value, cst.value, prec(self))
                if _do_sig(prec(self)): sig_off()
        return res

    def log1p(self, analytic=False):
        """
        Return ``log(1 + self)``, computed accurately when ``self`` is close to
        zero.

        INPUT:

        - ``analytic`` -- boolean (default: ``False``); if ``True``, return an
          indeterminate (not-a-number) value when the input ball touches
          the branch cut

        EXAMPLES::

            sage: eps = RBF(1e-50)
            sage: CBF(1+eps, eps).log()
            [+/- ...e-16] + [1.000000000000000e-50 +/- ...e-66]*I
            sage: CBF(eps, eps).log1p()
            [1.000000000000000e-50 +/- ...e-68] + [1.00000000000000e-50 +/- ...e-66]*I
            sage: CBF(-3/2).log1p(analytic=True)
            nan + nan*I

        TESTS::

            sage: CBF(-1/2).log1p(analytic=True)
            [-0.693147180559945 +/- ...e-16]
        """
        cdef ComplexBall res = self._new()
        if (analytic and arb_contains_zero(acb_imagref(self.value))
                     and not arb_gt_neg_one(acb_realref(self.value))):
            acb_indeterminate(res.value)
        else:
            if _do_sig(prec(self)): sig_on()
            acb_log1p(res.value, self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        return res

    def exp(self):
        """
        Return the exponential of this ball.

        .. SEEALSO:: :meth:`exppii`

        EXAMPLES::

            sage: CBF(i*pi).exp()                                                       # needs sage.symbolic
            [-1.00000000000000 +/- ...e-16] + [+/- ...e-16]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_exp(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def exppii(self):
        """
        Return ``exp(pi*i*self)``.

        EXAMPLES::

            sage: CBF(1/2).exppii()
            1.000000000000000*I
            sage: CBF(0, -1/pi).exppii()                                                # needs sage.symbolic
            [2.71828182845904 +/- ...e-15]
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_exp_pi_i(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def sin(self):
        """
        Return the sine of this ball.

        EXAMPLES::

            sage: CBF(i*pi).sin()                                                       # needs sage.symbolic
            [11.54873935725775 +/- ...e-15]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_sin(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def cos(self):
        """
        Return the cosine of this ball.

        EXAMPLES::

            sage: CBF(i*pi).cos()                                                       # needs sage.symbolic
            [11.59195327552152 +/- ...e-15]
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_cos(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def tan(self):
        """
        Return the tangent of this ball.

        EXAMPLES::

            sage: CBF(pi/2, 1/10).tan()                                                 # needs sage.symbolic
            [+/- ...e-14] + [10.03331113225399 +/- ...e-15]*I
            sage: CBF(pi/2).tan()                                                       # needs sage.symbolic
            nan
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_tan(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def cot(self):
        """
        Return the cotangent of this ball.

        EXAMPLES::

            sage: CBF(pi, 1/10).cot()                                                   # needs sage.symbolic
            [+/- ...e-14] + [-10.03331113225399 +/- ...e-15]*I
            sage: CBF(pi).cot()                                                         # needs sage.symbolic
            nan
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_cot(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def sec(self):
        """
        Return the secant of this ball.

        EXAMPLES::

            sage: CBF(1, 1).sec()
            [0.498337030555187 +/- ...e-16] + [0.591083841721045 +/- ...e-16]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_sec(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def csc(self):
        """
        Return the cosecant of this ball.

        EXAMPLES::

            sage: CBF(1, 1).csc()
            [0.621518017170428 +/- ...e-16] + [-0.303931001628426 +/- ...e-16]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_csc(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def sinh(self):
        """
        Return the hyperbolic sine of this ball.

        EXAMPLES::

            sage: CBF(1, 1).sinh()
            [0.634963914784736 +/- ...e-16] + [1.298457581415977 +/- ...e-16]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_sinh(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def cosh(self):
        """
        Return the hyperbolic cosine of this ball.

        EXAMPLES::

            sage: CBF(1, 1).cosh()
            [0.833730025131149 +/- ...e-16] + [0.988897705762865 +/- ...e-16]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_cosh(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def tanh(self):
        """
        Return the hyperbolic tangent of this ball.

        EXAMPLES::

            sage: CBF(1, 1).tanh()
            [1.083923327338694 +/- ...e-16] + [0.2717525853195117 +/- ...e-17]*I
            sage: CBF(0, pi/2).tanh()                                                   # needs sage.symbolic
            nan*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_tanh(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def coth(self):
        """
        Return the hyperbolic cotangent of this ball.

        EXAMPLES::

            sage: CBF(1, 1).coth()
            [0.868014142895925 +/- ...e-16] + [-0.2176215618544027 +/- ...e-17]*I
            sage: CBF(0, pi).coth()                                                     # needs sage.symbolic
            nan*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_coth(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def sech(self):
        """
        Return the hyperbolic secant of this ball.

        EXAMPLES::

            sage: CBF(pi/2, 1/10).sech()                                                # needs sage.symbolic
            [0.397174529918189 +/- ...e-16] + [-0.0365488656274242 +/- ...e-17]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_sech(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def csch(self):
        """
        Return the hyperbolic cosecant of this ball.

        EXAMPLES::

            sage: CBF(1, 1).csch()
            [0.303931001628426 +/- ...e-16] + [-0.621518017170428 +/- ...e-16]*I
            sage: CBF(i*pi).csch()                                                      # needs sage.symbolic
            nan*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_csch(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arcsin(self, analytic=False):
        """
        Return the arcsine of this ball.

        INPUT:

        - ``analytic`` -- boolean (default: ``False``); if ``True``, return an
          indeterminate (not-a-number) value when the input ball touches
          the branch cut

        EXAMPLES::

            sage: CBF(1+i).arcsin()
            [0.66623943249252 +/- ...e-15] + [1.06127506190504 +/- ...e-15]*I
            sage: CBF(1, RIF(0,1/1000)).arcsin()
            [1.6 +/- 0.0619] + [+/- 0.0322]*I
            sage: CBF(1, RIF(0,1/1000)).arcsin(analytic=True)
            nan + nan*I
        """
        cdef ComplexBall res = self._new()
        if (analytic and arb_contains_zero(acb_imagref(self.value))
                     and not arb_contained_unit_interval(acb_realref(self.value))):
            acb_indeterminate(res.value)
        else:
            if _do_sig(prec(self)): sig_on()
            acb_asin(res.value, self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        return res

    def arccos(self, analytic=False):
        """
        Return the arccosine of this ball.

        INPUT:

        - ``analytic`` -- boolean (default: ``False``); if ``True``, return an
          indeterminate (not-a-number) value when the input ball touches
          the branch cut

        EXAMPLES::

            sage: CBF(1+i).arccos()
            [0.90455689430238 +/- ...e-15] + [-1.06127506190504 +/- ...e-15]*I
            sage: CBF(-1).arccos()
            [3.141592653589793 +/- ...e-16]
            sage: CBF(-1).arccos(analytic=True)
            nan + nan*I
        """
        cdef ComplexBall res = self._new()
        if (analytic and arb_contains_zero(acb_imagref(self.value))
                     and not arb_contained_unit_interval(acb_realref(self.value))):
            acb_indeterminate(res.value)
        else:
            if _do_sig(prec(self)): sig_on()
            acb_acos(res.value, self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        return res

    def arctan(self, analytic=False):
        """
        Return the arctangent of this ball.

        INPUT:

        - ``analytic`` -- boolean (default: ``False``); if ``True``, return an
          indeterminate (not-a-number) value when the input ball touches
          the branch cut

        EXAMPLES::

            sage: CBF(1+i).arctan()
            [1.017221967897851 +/- ...e-16] + [0.4023594781085251 +/- ...e-17]*I
            sage: CBF(i).arctan()
            nan + nan*I
            sage: CBF(2*i).arctan()
            [1.570796326794897 +/- ...e-16] + [0.549306144334055 +/- ...e-16]*I
            sage: CBF(2*i).arctan(analytic=True)
            nan + nan*I
        """
        cdef ComplexBall res = self._new()
        if (analytic and arb_contains_zero(acb_realref(self.value))
                     and not arb_contained_unit_interval(acb_imagref(self.value))):
            acb_indeterminate(res.value)
        else:
            if _do_sig(prec(self)): sig_on()
            acb_atan(res.value, self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        return res

    def arcsinh(self, analytic=False):
        """
        Return the hyperbolic arcsine of this ball.

        INPUT:

        - ``analytic`` -- boolean (default: ``False``); if ``True``, return an
          indeterminate (not-a-number) value when the input ball touches
          the branch cut

        EXAMPLES::

            sage: CBF(1+i).arcsinh()
            [1.06127506190504 +/- ...e-15] + [0.66623943249252 +/- ...e-15]*I
            sage: CBF(2*i).arcsinh()
            [1.31695789692482 +/- ...e-15] + [1.570796326794897 +/- ...e-16]*I
            sage: CBF(2*i).arcsinh(analytic=True)
            nan + nan*I
        """
        cdef ComplexBall res = self._new()
        if (analytic and arb_contains_zero(acb_realref(self.value))
                     and not arb_contained_unit_interval(acb_imagref(self.value))):
            acb_indeterminate(res.value)
        else:
            if _do_sig(prec(self)): sig_on()
            acb_asinh(res.value, self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        return res

    def arccosh(self, analytic=False):
        """
        Return the hyperbolic arccosine of this ball.

        INPUT:

        - ``analytic`` -- boolean (default: ``False``); if ``True``, return an
          indeterminate (not-a-number) value when the input ball touches
          the branch cut

        EXAMPLES::

            sage: CBF(1+i).arccosh()
            [1.061275061905035 +/- ...e-16] + [0.904556894302381 +/- ...e-16]*I
            sage: CBF(-2).arccosh()
            [1.316957896924817 +/- ...e-16] + [3.141592653589793 +/- ...e-16]*I
            sage: CBF(-2).arccosh(analytic=True)
            nan + nan*I
        """
        cdef ComplexBall res = self._new()
        if (analytic and arb_contains_zero(acb_imagref(self.value))
                     and not arb_gt_neg_one(acb_realref(self.value))):
            acb_indeterminate(res.value)
        else:
            if _do_sig(prec(self)): sig_on()
            acb_acosh(res.value, self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        return res

    def arctanh(self, analytic=False):
        """
        Return the hyperbolic arctangent of this ball.

        INPUT:

        - ``analytic`` -- boolean (default: ``False``); if ``True``, return an
          indeterminate (not-a-number) value when the input ball touches
          the branch cut

        EXAMPLES::

            sage: CBF(1+i).arctanh()
            [0.4023594781085251 +/- ...e-17] + [1.017221967897851 +/- ...e-16]*I
            sage: CBF(-2).arctanh()
            [-0.549306144334055 +/- ...e-16] + [1.570796326794897 +/- ...e-16]*I
            sage: CBF(-2).arctanh(analytic=True)
            nan + nan*I
        """
        cdef ComplexBall res = self._new()
        if (analytic and arb_contains_zero(acb_imagref(self.value))
                     and not arb_contained_unit_interval(acb_realref(self.value))):
            acb_indeterminate(res.value)
        else:
            if _do_sig(prec(self)): sig_on()
            acb_atanh(res.value, self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        return res

    # Special functions

    def gamma(self, z=None):
        """
        Return the image of this ball by the Euler Gamma function (if
        ``z = None``) or the incomplete Gamma function (otherwise).

        EXAMPLES::

            sage: CBF(1, 1).gamma() # abs tol 1e-15
            [0.498015668118356 +/- 1.26e-16] + [-0.1549498283018107 +/- 8.43e-17]*I
            sage: CBF(-1).gamma()
            nan
            sage: CBF(1, 1).gamma(0) # abs tol 1e-15
            [0.498015668118356 +/- 1.26e-16] + [-0.1549498283018107 +/- 8.43e-17]*I
            sage: CBF(1, 1).gamma(100)
            [-3.6143867454139e-45 +/- ...e-59] + [-3.7022961377791e-44 +/- ...e-58]*I
            sage: CBF(1, 1).gamma(CLF(i)) # abs tol 1e-14
            [0.328866841935004 +/- 7.07e-16] + [-0.189749450456210 +/- 9.05e-16]*I
        """
        cdef ComplexBall my_z
        cdef ComplexBall res = self._new()
        if z is None:
            if _do_sig(prec(self)): sig_on()
            acb_gamma(res.value, self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        else:
            my_z = self._parent.coerce(z)
            if _do_sig(prec(self)): sig_on()
            acb_hypgeom_gamma_upper(res.value, self.value, my_z.value, 0, prec(self))
            if _do_sig(prec(self)): sig_off()
        return res

    gamma_inc = gamma

    def log_gamma(self, analytic=False):
        r"""
        Return the image of this ball by the logarithmic Gamma function.

        The branch cut of the logarithmic gamma function is placed on the
        negative half-axis, which means that
        ``log_gamma(z) + log z = log_gamma(z+1)`` holds for all `z`,
        whereas ``log_gamma(z) != log(gamma(z))`` in general.

        INPUT:

        - ``analytic`` -- boolean (default: ``False``); if ``True``, return an
          indeterminate (not-a-number) value when the input ball touches
          the branch cut

        EXAMPLES::

            sage: CBF(1000, 1000).log_gamma()
            [5466.22252162990 +/- ...e-12] + [7039.33429191119 +/- ...e-12]*I
            sage: CBF(-1/2).log_gamma()
            [1.265512123484645 +/- ...e-16] + [-3.141592653589793 +/- ...e-16]*I
            sage: CBF(-1).log_gamma()
            nan + ...*I
            sage: CBF(-3/2).log_gamma() # abs tol 1e-14
            [0.860047015376481 +/- 3.82e-16] + [-6.283185307179586 +/- 6.77e-16]*I
            sage: CBF(-3/2).log_gamma(analytic=True)
            nan + nan*I
        """
        cdef ComplexBall res = self._new()
        if (analytic and arb_contains_zero(acb_imagref(self.value))
                     and arb_contains_nonpositive(acb_realref(self.value))):
            acb_indeterminate(res.value)
        else:
            if _do_sig(prec(self)): sig_on()
            acb_lgamma(res.value, self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        return res

    def rgamma(self):
        """
        Compute the reciprocal gamma function with argument ``self``.

        EXAMPLES::

            sage: CBF(6).rgamma()
            [0.00833333333333333 +/- ...e-18]
            sage: CBF(-1).rgamma()
            0
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_rgamma(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def psi(self, n=None):
        """
        Compute the digamma function with argument ``self``.

        If ``n`` is provided, compute the polygamma function of order ``n``
        and argument ``self``.

        EXAMPLES::

            sage: CBF(1, 1).psi()
            [0.0946503206224770 +/- ...e-17] + [1.076674047468581 +/- ...e-16]*I
            sage: CBF(-1).psi()
            nan
            sage: CBF(1,1).psi(10)
            [56514.8269344249 +/- ...e-11] + [56215.1218005823 +/- ...e-11]*I
        """
        cdef ComplexBall my_n
        cdef ComplexBall result = self._new()
        if n is None:
            if _do_sig(prec(self)): sig_on()
            acb_digamma(result.value, self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        else:
            my_n = self._parent.coerce(n)
            if _do_sig(prec(self)): sig_on()
            acb_polygamma(result.value, my_n.value, self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        return result

    def zeta(self, a=None):
        """
        Return the image of this ball by the Hurwitz zeta function.

        For ``a = None``, this computes the Riemann zeta function.

        EXAMPLES::

            sage: CBF(1, 1).zeta()
            [0.5821580597520036 +/- ...e-17] + [-0.9268485643308071 +/- ...e-17]*I
            sage: CBF(1, 1).zeta(1)
            [0.5821580597520036 +/- ...e-17] + [-0.9268485643308071 +/- ...e-17]*I
            sage: CBF(1, 1).zeta(1/2)
            [1.497919876084167 +/- ...e-16] + [0.2448655353684164 +/- ...e-17]*I
            sage: CBF(1, 1).zeta(CBF(1, 1))
            [-0.3593983122202835 +/- ...e-17] + [-2.875283329756940 +/- ...e-16]*I
            sage: CBF(1, 1).zeta(-1)
            nan + nan*I
        """
        cdef ComplexBall a_ball
        cdef ComplexBall res = self._new()
        if a is None:
            if _do_sig(prec(self)): sig_on()
            acb_zeta(res.value, self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        else:
            a_ball = self._parent.coerce(a)
            if _do_sig(prec(self)): sig_on()
            acb_hurwitz_zeta(res.value, self.value, a_ball.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        return res

    def zetaderiv(self, k):
        r"""
        Return the image of this ball by the `k`-th derivative of the Riemann
        zeta function.

        For a more flexible interface, see the low-level method
        ``_zeta_series`` of polynomials with complex ball coefficients.

        EXAMPLES::

            sage: CBF(1/2, 3).zetaderiv(1)
            [0.191759884092721...] + [-0.073135728865928...]*I
            sage: CBF(2).zetaderiv(3)
            [-6.0001458028430...]
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        Pol = PolynomialRing(self._parent, 'x')
        ser = Pol([self, 1])._zeta_series(k + 1)
        return ser[k]*ZZ.coerce(k).factorial()

    def lambert_w(self, branch=0):
        r"""
        Return the image of this ball by the specified branch of the Lambert W
        function.

        EXAMPLES::

            sage: CBF(1 + I).lambert_w()
            [0.6569660692304...] + [0.3254503394134...]*I
            sage: CBF(1 + I).lambert_w(2)
            [-2.1208839379437...] + [11.600137110774...]*I
            sage: CBF(1 + I).lambert_w(2^100)
            [-70.806021532123...] + [7.9648836259913...]*I
        """
        cdef fmpz_t _branch
        fmpz_init(_branch)
        fmpz_set_mpz(_branch, (<Integer> Integer(branch)).value)
        cdef ComplexBall res = self._new()
        sig_on()
        acb_lambertw(res.value, self.value, _branch, 0, prec(self))
        sig_off()
        fmpz_clear(_branch)
        return res

    def polylog(self, s):
        r"""
        Return the polylogarithm `\operatorname{Li}_s(\mathrm{self})`.

        EXAMPLES::

            sage: CBF(2).polylog(1)
            [+/- ...e-15] + [-3.14159265358979 +/- ...e-15]*I
            sage: CBF(1, 1).polylog(CBF(1, 1))
            [0.3708160030469 +/- ...e-14] + [2.7238016577979 +/- ...e-14]*I

        TESTS::

            sage: CBF(2).polylog(1r)
            [+/- ...e-15] + [-3.14159265358979 +/- ...e-15]*I
        """
        cdef ComplexBall s_as_ball
        cdef Integer s_as_Integer
        cdef ComplexBall res = self._new()
        try:
            s_as_Integer = ZZ.coerce(s)
            if mpz_fits_slong_p(s_as_Integer.value):
                if _do_sig(prec(self)): sig_on()
                acb_polylog_si(res.value, mpz_get_si(s_as_Integer.value), self.value, prec(self))
                if _do_sig(prec(self)): sig_off()
                return res
        except TypeError:
            pass
        s_as_ball = self._parent.coerce(s)
        if _do_sig(prec(self)): sig_on()
        acb_polylog(res.value, s_as_ball.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def barnes_g(self):
        """
        Return the Barnes G-function of ``self``.

        EXAMPLES::

            sage: CBF(-4).barnes_g()
            0
            sage: CBF(8).barnes_g()
            24883200.00000000
            sage: CBF(500,10).barnes_g()
            [4.54078781e+254873 +/- ...e+254864] + [8.65835455e+254873 +/- ...e+254864]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_barnes_g(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def log_barnes_g(self):
        """
        Return the logarithmic Barnes G-function of ``self``.

        EXAMPLES::

            sage: CBF(10^100).log_barnes_g()
            [1.14379254649702e+202 +/- ...e+187]
            sage: CBF(0,1000).log_barnes_g()
            [-2702305.04929258 +/- ...e-9] + [-790386.325561423 +/- ...e-10]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_log_barnes_g(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def agm1(self):
        """
        Return the arithmetic-geometric mean of 1 and ``self``.

        The arithmetic-geometric mean is defined such that the function is
        continuous in the complex plane except for a branch cut along the
        negative half axis (where it is continuous from above). This
        corresponds to always choosing an "optimal" branch for the square root
        in the arithmetic-geometric mean iteration.

        EXAMPLES::

            sage: CBF(0, -1).agm1()
            [0.599070117367796 +/- 3.9...e-16] + [-0.599070117367796 +/- 5.5...e-16]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_agm1(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def hypergeometric(self, a, b, bint regularized=False):
        r"""
        Return the generalized hypergeometric function of ``self``.

        INPUT:

        - ``a`` -- upper parameters; list of complex numbers that coerce into
          this ball's parent

        - ``b`` -- lower parameters; list of complex numbers that coerce into
          this ball's parent

        - ``regularized`` -- if ``True``, the regularized generalized
          hypergeometric function is computed

        OUTPUT: the generalized hypergeometric function defined by

        .. MATH::

            {}_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;z)
            = \sum_{k=0}^\infty \frac{(a_1)_k\dots(a_p)_k}{(b_1)_k\dots(b_q)_k} \frac {z^k} {k!}

        extended using analytic continuation or regularization when the sum
        does not converge.

        The regularized generalized hypergeometric function

        .. MATH::

            {}_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;z)
            = \sum_{k=0}^\infty \frac{(a_1)_k\dots(a_p)_k}{\Gamma(b_1+k)\dots\Gamma(b_q+k)} \frac {z^k} {k!}

        is well-defined even when the lower parameters are nonpositive
        integers. Currently, this is only supported for some `p` and `q`.

        EXAMPLES::

            sage: CBF(1, pi/2).hypergeometric([], [])                                   # needs sage.symbolic
            [+/- ...e-16] + [2.71828182845904 +/- ...e-15]*I

            sage: CBF(1, pi).hypergeometric([1/4], [1/4])                               # needs sage.symbolic
            [-2.7182818284590 +/- ...e-14] + [+/- ...e-14]*I

            sage: CBF(1000, 1000).hypergeometric([10], [AA(sqrt(2))])                   # needs sage.symbolic
            [9.79300951360e+454 +/- ...e+442] + [5.522579106816e+455 +/- ...e+442]*I
            sage: CBF(1000, 1000).hypergeometric([100], [AA(sqrt(2))])                  # needs sage.symbolic
            [1.27967355557e+590 +/- ...e+578] + [-9.32333491987e+590 +/- ...e+578]*I

            sage: CBF(0, 1).hypergeometric([], [1/2, 1/3, 1/4])
            [-3.7991962344383 +/- ...e-14] + [23.878097177805 +/- ...e-13]*I

            sage: CBF(0).hypergeometric([1], [])
            1.000000000000000
            sage: CBF(1, 1).hypergeometric([1], [])
            1.000000000000000*I

            sage: CBF(2+3*I).hypergeometric([1/4,1/3],[1/2]) # abs tol 1e-14
            [0.7871684267473 +/- 6.79e-14] + [0.2749254173721 +/- 8.82e-14]*I
            sage: CBF(2+3*I).hypergeometric([1/4,1/3],[1/2],regularized=True)
            [0.4441122268685 +/- 3...e-14] + [0.1551100567338 +/- 5...e-14]*I

            sage: CBF(5).hypergeometric([2,3], [-5])
            nan + nan*I
            sage: CBF(5).hypergeometric([2,3], [-5], regularized=True)
            [5106.925964355 +/- ...e-10]

            sage: CBF(2016).hypergeometric([], [2/3]) # abs tol 1e+26
            [2.0256426923278e+38 +/- 9.59e+24]
            sage: CBF(-2016).hypergeometric([], [2/3], regularized=True)
            [-0.0005428550847 +/- ...e-14]

            sage: CBF(-7).hypergeometric([4], [])
            0.0002441406250000000

            sage: CBF(0, 3).hypergeometric([CBF(1,1)], [-4], regularized=True)
            [239.514000752841 +/- ...e-13] + [105.175157349015 +/- ...e-13]*I

        TESTS::

            sage: CBF(0, 1).hypergeometric([QQbar(sqrt(2)), RLF(pi)], [1r, 1/2])        # needs sage.symbolic
            [-8.7029449215408 +/- ...e-14] + [-0.8499070546106 +/- ...e-14]*I
        """
        cdef ComplexBall tmp, my_a, my_b, my_c
        cdef ComplexBall res = self._new()
        cdef long p = len(a)
        cdef long q = len(b)
        if p == q == 0:
            return self.exp()
        if p == 1 and q == 0:
            my_a = self._parent.coerce(a[0])
            return (1-self)**(-my_a)
        if p == 0 and q == 1:
            my_b = self._parent.coerce(b[0])
            if _do_sig(prec(self)): sig_on()
            acb_hypgeom_0f1(res.value, my_b.value, self.value,
                regularized, prec(self))
            if _do_sig(prec(self)): sig_off()
            return res
        if p == q == 1:
            my_a = self._parent.coerce(a[0])
            my_b = self._parent.coerce(b[0])
            if _do_sig(prec(self)): sig_on()
            acb_hypgeom_m(res.value, my_a.value, my_b.value, self.value,
                regularized, prec(self))
            if _do_sig(prec(self)): sig_off()
            return res
        if p == 2 and q == 1:
            my_a = self._parent.coerce(a[0])
            my_b = self._parent.coerce(a[1])
            my_c = self._parent.coerce(b[0])
            if _do_sig(prec(self)): sig_on()
            acb_hypgeom_2f1(res.value, my_a.value, my_b.value, my_c.value,
                self.value, regularized, prec(self))
            if _do_sig(prec(self)): sig_off()
            return res
        if regularized:
            raise NotImplementedError("regularized=True not yet supported in this case")
        cdef long i1 = -1
        cdef long s
        try:
            i1 = a.index(1)
            s = 1
        except ValueError:
            s = 0
        cdef acb_ptr vec_a = _acb_vec_init(p - s)
        cdef acb_ptr vec_b = _acb_vec_init(q + 1 - s)
        cdef long j = 0
        for i in range(p):
            if i != i1:
                tmp = self._parent.coerce(a[i])
                acb_set(&(vec_a[j]), tmp.value)
                j += 1
        for i in range(q):
            tmp = self._parent.coerce(b[i])
            acb_set(&(vec_b[i]), tmp.value)
        if s == 0:
            acb_one(&(vec_b[q]))
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_pfq_direct(res.value, vec_a, p - s, vec_b, q + 1 - s,
                               self.value, -1, prec(self))
        if _do_sig(prec(self)): sig_off()
        _acb_vec_clear(vec_b, q + 1 - s)
        _acb_vec_clear(vec_a, p - s)
        return res

    def hypergeometric_U(self, a, b):
        """
        Return the Tricomi confluent hypergeometric function U(a, b, self) of
        this ball.

        EXAMPLES::

            sage: CBF(1000, 1000).hypergeometric_U(RLF(pi), -100)                       # needs sage.symbolic
            [-7.261605907166e-11 +/- ...e-24] + [-7.928136216391e-11 +/- ...e-24]*I
            sage: CBF(1000, 1000).hypergeometric_U(0, -100)
            1.000000000000000
        """
        cdef ComplexBall res = self._new()
        cdef ComplexBall my_a = self._parent.coerce(a)
        cdef ComplexBall my_b = self._parent.coerce(b)
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_u(res.value, my_a.value, my_b.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def erf(self):
        """
        Return the error function with argument ``self``.

        EXAMPLES::

            sage: CBF(1, 1).erf()
            [1.316151281697947 +/- ...e-16] + [0.1904534692378347 +/- ...e-17]*I
        """

        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_erf(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def erfc(self):
        """
        Compute the complementary error function with argument ``self``.

        EXAMPLES::

            sage: CBF(20).erfc() # abs tol 1e-190
            [5.39586561160790e-176 +/- 6.73e-191]
            sage: CBF(100, 100).erfc()
            [0.00065234366376858 +/- ...e-18] + [-0.00393572636292141 +/- ...e-18]*I
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_erfc(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def airy(self):
        """
        Return the Airy functions Ai, Ai', Bi, Bi' with argument ``self``,
        evaluated simultaneously.

        EXAMPLES::

            sage: CBF(10*pi).airy()                                                     # needs sage.symbolic
            ([1.2408955946101e-52 +/- ...e-66],
             [-6.965048886977e-52 +/- ...e-65],
             [2.2882956833435e+50 +/- ...e+36],
             [1.2807602335816e+51 +/- ...e+37])
            sage: ai, aip, bi, bip = CBF(1,2).airy()
            sage: (ai * bip - bi * aip) * CBF(pi)                                       # needs sage.symbolic
            [1.0000000000000 +/- ...e-15] + [+/- ...e-16]*I
        """
        cdef ComplexBall ai = self._new()
        cdef ComplexBall aip = self._new()
        cdef ComplexBall bi = self._new()
        cdef ComplexBall bip = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_airy(ai.value, aip.value, bi.value, bip.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return (ai, aip, bi, bip)

    def airy_ai(self):
        """
        Return the Airy function Ai with argument ``self``.

        EXAMPLES::

            sage: CBF(1,2).airy_ai()
            [-0.2193862549814276 +/- ...e-17] + [-0.1753859114081094 +/- ...e-17]*I
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_airy(result.value, NULL, NULL, NULL, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def airy_ai_prime(self):
        """
        Return the Airy function derivative Ai' with argument ``self``.

        EXAMPLES::

            sage: CBF(1,2).airy_ai_prime()
            [0.1704449781789148 +/- ...e-17] + [0.387622439413295 +/- ...e-16]*I
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_airy(NULL, result.value, NULL, NULL, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def airy_bi(self):
        """
        Return the Airy function Bi with argument ``self``.

        EXAMPLES::

            sage: CBF(1,2).airy_bi()
            [0.0488220324530612 +/- ...e-17] + [0.1332740579917484 +/- ...e-17]*I
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_airy(NULL, NULL, result.value, NULL, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def airy_bi_prime(self):
        """
        Return the Airy function derivative Bi' with argument ``self``.

        EXAMPLES::

            sage: CBF(1,2).airy_bi_prime()
            [-0.857239258605362 +/- ...e-16] + [0.4955063363095674 +/- ...e-17]*I
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_airy(NULL, NULL, NULL, result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def bessel_J(self, nu):
        """
        Return the Bessel function of the first kind with argument ``self``
        and index ``nu``.

        EXAMPLES::

            sage: CBF(1, 1).bessel_J(1)
            [0.614160334922903 +/- ...e-16] + [0.365028028827088 +/- ...e-16]*I
            sage: CBF(100, -100).bessel_J(1/3)
            [1.108431870251e+41 +/- ...e+28] + [-8.952577603125e+41 +/- ...e+28]*I
        """
        cdef ComplexBall result = self._new()
        cdef ComplexBall my_nu = self._parent.coerce(nu)
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_bessel_j(result.value, my_nu.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def bessel_J_Y(self, nu):
        """
        Return the Bessel function of the first and second kind with argument
        ``self`` and index ``nu``, computed simultaneously.

        EXAMPLES::

            sage: J, Y = CBF(1, 1).bessel_J_Y(1)
            sage: J - CBF(1, 1).bessel_J(1)
            [+/- ...e-16] + [+/- ...e-16]*I
            sage: Y - CBF(1, 1).bessel_Y(1)
            [+/- ...e-14] + [+/- ...e-14]*I
        """
        cdef ComplexBall result1 = self._new()
        cdef ComplexBall result2 = self._new()
        cdef ComplexBall my_nu = self._parent.coerce(nu)
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_bessel_jy(result1.value, result2.value,
            my_nu.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result1, result2

    def bessel_Y(self, nu):
        """
        Return the Bessel function of the second kind with argument ``self``
        and index ``nu``.

        EXAMPLES::

            sage: CBF(1, 1).bessel_Y(1)
            [-0.6576945355913 +/- ...e-14] + [0.6298010039929 +/- ...e-14]*I
            sage: CBF(100, -100).bessel_Y(1/3)
            [-8.952577603125e+41 +/- ...e+28] + [-1.108431870251e+41 +/- ...e+28]*I
        """
        cdef ComplexBall result = self._new()
        cdef ComplexBall my_nu = self._parent.coerce(nu)
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_bessel_y(result.value, my_nu.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def bessel_I(self, nu):
        """
        Return the modified Bessel function of the first kind with argument ``self``
        and index ``nu``.

        EXAMPLES::

            sage: CBF(1, 1).bessel_I(1)
            [0.365028028827088 +/- ...e-16] + [0.614160334922903 +/- ...e-16]*I
            sage: CBF(100, -100).bessel_I(1/3)
            [5.4362189595644e+41 +/- ...e+27] + [7.1989436985321e+41 +/- ...e+27]*I
        """
        cdef ComplexBall result = self._new()
        cdef ComplexBall my_nu = self._parent.coerce(nu)
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_bessel_i(result.value, my_nu.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def bessel_K(self, nu):
        """
        Return the modified Bessel function of the second kind with argument
        ``self`` and index ``nu``.

        EXAMPLES::

            sage: CBF(1, 1).bessel_K(0)
            [0.08019772694652 +/- ...e-15] + [-0.357277459285330 +/- ...e-16]*I
            sage: CBF(1, 1).bessel_K(1)
            [0.02456830552374 +/- ...e-15] + [-0.45971947380119 +/- ...e-15]*I
            sage: CBF(100, 100).bessel_K(QQbar(i))
            [3.8693896656383e-45 +/- ...e-59] + [5.507100423418e-46 +/- ...e-59]*I
        """
        cdef ComplexBall result = self._new()
        cdef ComplexBall my_nu = self._parent.coerce(nu)
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_bessel_k(result.value, my_nu.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def exp_integral_e(self, s):
        """
        Return the image of this ball by the generalized exponential integral
        with index ``s``.

        EXAMPLES::

            sage: CBF(1+i).exp_integral_e(1)
            [0.00028162445198 +/- ...e-15] + [-0.17932453503936 +/- ...e-15]*I
            sage: CBF(1+i).exp_integral_e(QQbar(i))
            [-0.10396361883964 +/- ...e-15] + [-0.16268401277783 +/- ...e-15]*I
        """
        cdef ComplexBall res = self._new()
        cdef ComplexBall my_s = self._parent.coerce(s)
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_expint(res.value, my_s.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def Ei(self):
        """
        Return the exponential integral with argument ``self``.

        EXAMPLES::

            sage: CBF(1, 1).Ei()
            [1.76462598556385 +/- ...e-15] + [2.38776985151052 +/- ...e-15]*I
            sage: CBF(0).Ei()
            nan

        TESTS:

            sage: CBF(Ei(I))  # abs tol 1e-16                                           # needs sage.symbolic
            [0.337403922900968 +/- 3.76e-16] + [2.51687939716208 +/- 2.01e-15]*I
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_ei(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def Si(self):
        """
        Return the sine integral with argument ``self``.

        EXAMPLES::

            sage: CBF(1, 1).Si()
            [1.10422265823558 +/- ...e-15] + [0.88245380500792 +/- ...e-15]*I
            sage: CBF(0).Si()
            0

        TESTS:

            sage: CBF(Si(I))                                                            # needs sage.symbolic
            [1.05725087537573 +/- 2.77e-15]*I
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_si(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    sin_integral = Si # as for the symbolic function

    def Ci(self):
        """
        Return the cosine integral with argument ``self``.

        EXAMPLES::

            sage: CBF(1, 1).Ci()
            [0.882172180555936 +/- ...e-16] + [0.287249133519956 +/- ...e-16]*I
            sage: CBF(0).Ci()
            nan + nan*I

        TESTS:

            sage: CBF(Ci(I))  # abs tol 1e-17                                           # needs sage.symbolic
            [0.837866940980208 +/- 4.72e-16] + [1.570796326794897 +/- 5.54e-16]*I
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_ci(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    cos_integral = Ci # as for the symbolic function

    def Shi(self):
        """
        Return the hyperbolic sine integral with argument ``self``.

        EXAMPLES::

            sage: CBF(1, 1).Shi()
            [0.88245380500792 +/- ...e-15] + [1.10422265823558 +/- ...e-15]*I
            sage: CBF(0).Shi()
            0

        TESTS:

            sage: CBF(Shi(I))  # abs tol 1e-15                                          # needs sage.symbolic
            [0.946083070367183 +/- 9.22e-16]*I
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_shi(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    sinh_integral = Shi

    def Chi(self):
        """
        Return the hyperbolic cosine integral with argument ``self``.

        EXAMPLES::

            sage: CBF(1, 1).Chi()
            [0.882172180555936 +/- ...e-16] + [1.28354719327494 +/- ...e-15]*I
            sage: CBF(0).Chi()
            nan + nan*I

        TESTS:

            sage: CBF(Chi(I))  # abs tol 1e-16                                          # needs sage.symbolic
            [0.337403922900968 +/- 3.25e-16] + [1.570796326794897 +/- 5.54e-16]*I
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_chi(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    cosh_integral = Chi

    def li(self, bint offset=False):
        """
        Return the logarithmic integral with argument ``self``.

        If ``offset`` is True, return the offset logarithmic integral.

        EXAMPLES::

            sage: CBF(1, 1).li()
            [0.61391166922120 +/- ...e-15] + [2.05958421419258 +/- ...e-15]*I
            sage: CBF(0).li()
            0
            sage: CBF(0).li(offset=True)
            [-1.045163780117493 +/- ...e-16]
            sage: li(0).n()                                                             # needs sage.symbolic
            0.000000000000000
            sage: Li(0).n()                                                             # needs sage.symbolic
            -1.04516378011749

        TESTS::

            sage: CBF(li(0))                                                            # needs sage.symbolic
            0
            sage: CBF(Li(0))                                                            # needs sage.symbolic
            [-1.04516378011749...]
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_li(result.value, self.value, offset, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    log_integral = li

    def Li(self):
        """
        Offset logarithmic integral.

        EXAMPLES::

            sage: CBF(0).Li()
            [-1.045163780117493 +/- ...e-16]
            sage: li(0).n()                                                             # needs sage.symbolic
            0.000000000000000
            sage: Li(0).n()                                                             # needs sage.symbolic
            -1.04516378011749
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_li(result.value, self.value, True, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    log_integral_offset = Li

    def jacobi_theta(self, tau):
        r"""
        Return the four Jacobi theta functions evaluated at the argument
        ``self`` (representing `z`) and the parameter ``tau`` which should lie
        in the upper half plane.

        The following definitions are used:

        .. MATH::

            \theta_1(z,\tau) = 2 q_{1/4} \sum_{n=0}^{\infty} (-1)^n q^{n(n+1)} \sin((2n+1) \pi z)

            \theta_2(z,\tau) = 2 q_{1/4} \sum_{n=0}^{\infty} q^{n(n+1)} \cos((2n+1) \pi z)

            \theta_3(z,\tau) = 1 + 2 \sum_{n=1}^{\infty} q^{n^2} \cos(2n \pi z)

            \theta_4(z,\tau) = 1 + 2 \sum_{n=1}^{\infty} (-1)^n q^{n^2} \cos(2n \pi z)

        where `q = \exp(\pi i \tau)` and `q_{1/4} = \exp(\pi i \tau / 4)`.
        Note that `z` is multiplied by `\pi`; some authors omit this factor.

        EXAMPLES::

            sage: CBF(3,-1/2).jacobi_theta(CBF(1/4,2))
            ([-0.186580562274757 +/- ...e-16] + [0.93841744788594 +/- ...e-15]*I,
             [-1.02315311037951 +/- ...e-15] + [-0.203600094532010 +/- ...e-16]*I,
             [1.030613911309632 +/- ...e-16] + [0.030613917822067 +/- ...e-16]*I,
             [0.969386075665498 +/- ...e-16] + [-0.030613917822067 +/- ...e-16]*I)

            sage: CBF(3,-1/2).jacobi_theta(CBF(1/4,-2))
            (nan + nan*I, nan + nan*I, nan + nan*I, nan + nan*I)

            sage: CBF(0).jacobi_theta(CBF(0,1))
            (0,
             [0.913579138156117 +/- ...e-16],
             [1.086434811213308 +/- ...e-16],
             [0.913579138156117 +/- ...e-16])
        """

        cdef ComplexBall res1 = self._new()
        cdef ComplexBall res2 = self._new()
        cdef ComplexBall res3 = self._new()
        cdef ComplexBall res4 = self._new()
        cdef ComplexBall my_tau = self._parent.coerce(tau)
        if _do_sig(prec(self)): sig_on()
        acb_modular_theta(res1.value, res2.value, res3.value, res4.value,
            self.value, my_tau.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res1, res2, res3, res4

    def modular_j(self):
        """
        Return the modular j-invariant with *tau* given by ``self``.

        EXAMPLES::

            sage: CBF(0,1).modular_j()
            [1728.0000000000 +/- ...e-11]
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_modular_j(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def modular_eta(self):
        """
        Return the Dedekind eta function with *tau* given by ``self``.

        EXAMPLES::

            sage: CBF(0,1).modular_eta()
            [0.768225422326057 +/- ...e-16]
            sage: CBF(12,1).modular_eta()
            [-0.768225422326057 +/- ...e-16]
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_modular_eta(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def modular_lambda(self):
        """
        Return the modular lambda function with *tau* given by ``self``.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: tau = CBF(sqrt(2),pi)
            sage: tau.modular_lambda()
            [-0.00022005123884157 +/- ...e-18] + [-0.00079787346459944 +/- ...e-18]*I
            sage: (tau + 2).modular_lambda()
            [-0.00022005123884157 +/- ...e-18] + [-0.00079787346459944 +/- ...e-18]*I
            sage: (tau / (1 - 2*tau)).modular_lambda()
            [-0.00022005123884 +/- ...e-15] + [-0.00079787346460 +/- ...e-15]*I
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_modular_lambda(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def modular_delta(self):
        """
        Return the modular discriminant with *tau* given by ``self``.

        EXAMPLES::

            sage: CBF(0,1).modular_delta()
            [0.0017853698506421 +/- ...e-17]
            sage: a, b, c, d = 2, 5, 1, 3
            sage: tau = CBF(1,3)
            sage: ((a*tau+b)/(c*tau+d)).modular_delta()
            [0.20921376655 +/- ...e-12] + [1.57611925523 +/- ...e-12]*I
            sage: (c*tau+d)^12 * tau.modular_delta()
            [0.20921376654986 +/- ...e-15] + [1.5761192552253 +/- ...e-14]*I
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_modular_delta(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def eisenstein(self, long n):
        r"""
        Return the first ``n`` entries in the sequence of Eisenstein series
        `G_4(\tau), G_6(\tau), G_8(\tau), \ldots` where *tau* is given
        by ``self``. The output is a list.

        EXAMPLES::

            sage: a, b, c, d = 2, 5, 1, 3
            sage: tau = CBF(1,3)
            sage: tau.eisenstein(4)
            [[2.1646498507193 +/- ...e-14],
             [2.0346794456073 +/- ...e-14],
             [2.0081609898081 +/- ...e-14],
             [2.0019857082706 +/- ...e-14]]
            sage: ((a*tau+b)/(c*tau+d)).eisenstein(3)[2]
            [331011.2004330 +/- ...e-8] + [-711178.1655746 +/- ...e-8]*I
            sage: (c*tau+d)^8 * tau.eisenstein(3)[2]
            [331011.20043304 +/- ...e-9] + [-711178.1655746 +/- ...e-8]*I
        """
        if n < 0:
            raise ValueError("n must be nonnegative")
        cdef acb_ptr vec_r = _acb_vec_init(n)
        if _do_sig(prec(self)): sig_on()
        acb_modular_eisenstein(vec_r, self.value, n, prec(self))
        if _do_sig(prec(self)): sig_off()
        result = [self._new() for i in range(n)]
        for i in range(n):
            acb_swap((<ComplexBall>(result[i])).value, &(vec_r[i]))
        _acb_vec_clear(vec_r, n)
        return result

    def elliptic_p(self, tau, n=None):
        r"""
        Return the Weierstrass elliptic function with lattice parameter ``tau``,
        evaluated at ``self``. The function is doubly periodic in ``self``
        with periods 1 and ``tau``, which should lie in the upper half plane.

        If ``n`` is given, return a list containing the first ``n``
        terms in the Taylor expansion at ``self``. In particular, with
        ``n`` = 2, compute the Weierstrass elliptic function together
        with its derivative, which generate the field of elliptic
        functions with periods 1 and ``tau``.

        EXAMPLES::

            sage: tau = CBF(1,4)
            sage: z = CBF(sqrt(2), sqrt(3))                                             # needs sage.symbolic
            sage: z.elliptic_p(tau)                                                     # needs sage.symbolic
            [-3.28920996772709 +/- ...e-15] + [-0.0003673767302933 +/- ...e-17]*I
            sage: (z + tau).elliptic_p(tau)                                             # needs sage.symbolic
            [-3.28920996772709 +/- ...e-15] + [-0.000367376730293 +/- ...e-16]*I
            sage: (z + 1).elliptic_p(tau)                                               # needs sage.symbolic
            [-3.28920996772709 +/- ...e-15] + [-0.0003673767302933 +/- ...e-17]*I

            sage: z.elliptic_p(tau, 3)                                                  # needs sage.symbolic
            [[-3.28920996772709 +/- ...e-15] + [-0.0003673767302933 +/- ...e-17]*I,
             [0.002473055794309 +/- ...e-16] + [0.003859554040267 +/- ...e-16]*I,
             [-0.01299087561709 +/- ...e-15] + [0.00725027521915 +/- ...e-15]*I]
            sage: (z + 3 + 4*tau).elliptic_p(tau, 3)                                    # needs sage.symbolic
            [[-3.28920996772709 +/- ...e-15] + [-0.00036737673029 +/- ...e-15]*I,
             [0.0024730557943 +/- ...e-14] + [0.0038595540403 +/- ...e-14]*I,
             [-0.01299087562 +/- ...e-12] + [0.00725027522 +/- ...e-12]*I]
        """
        cdef ComplexBall my_tau = self._parent.coerce(tau)
        cdef ComplexBall result
        cdef long nn
        cdef acb_ptr vec_r
        if n is None:
            result = self._new()
            if _do_sig(prec(self)): sig_on()
            acb_elliptic_p(result.value, self.value,
                my_tau.value, prec(self))
            if _do_sig(prec(self)): sig_off()
            return result
        else:
            nn = n
            if nn < 0:
                raise ValueError("n must be nonnegative")
            vec_r = _acb_vec_init(nn)
            if _do_sig(prec(self)): sig_on()
            acb_modular_elliptic_p_zpx(vec_r, self.value, my_tau.value, nn, prec(self))
            if _do_sig(prec(self)): sig_off()
            result_list = [self._new() for i in range(nn)]
            for i in range(nn):
                acb_swap((<ComplexBall>(result_list[i])).value, &(vec_r[i]))
            _acb_vec_clear(vec_r, nn)
            return result_list

    def elliptic_invariants(self):
        r"""
        Return the lattice invariants ``(g2, g3)``.

        EXAMPLES::

            sage: CBF(0,1).elliptic_invariants()
            ([189.07272012923 +/- ...e-12], [+/- ...e-12])
            sage: CBF(sqrt(2)/2, sqrt(2)/2).elliptic_invariants()                       # needs sage.symbolic
            ([+/- ...e-12] + [-332.5338031465...]*I,
             [1254.46842157...] + [1254.46842157...]*I)
        """
        cdef ComplexBall g2 = self._new()
        cdef ComplexBall g3 = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_elliptic_invariants(g2.value, g3.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return (g2, g3)

    def elliptic_roots(self):
        r"""
        Return the lattice roots ``(e1, e2, e3)`` of `4 z^3 - g_2 z - g_3`.

        EXAMPLES::

            sage: e1, e2, e3 = CBF(0,1).elliptic_roots()
            sage: e1, e2, e3
            ([6.8751858180204 +/- ...e-14],
             [+/- ...e-14],
             [-6.8751858180204 +/- ...e-14])
            sage: g2, g3 = CBF(0,1).elliptic_invariants()
            sage: 4 * e1^3 - g2 * e1 - g3
            [+/- ...e-11]
        """
        cdef ComplexBall e1 = self._new()
        cdef ComplexBall e2 = self._new()
        cdef ComplexBall e3 = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_elliptic_roots(e1.value, e2.value, e3.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return (e1, e2, e3)

    def elliptic_k(self):
        """
        Return the complete elliptic integral of the first kind evaluated
        at *m* given by ``self``.

        EXAMPLES::

            sage: CBF(2,3).elliptic_k()
            [1.04291329192852 +/- ...e-15] + [0.62968247230864 +/- ...e-15]*I
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_elliptic_k(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def elliptic_e(self):
        """
        Return the complete elliptic integral of the second kind evaluated
        at *m* given by ``self``.

        EXAMPLES::

            sage: CBF(2,3).elliptic_e()
            [1.472797144959 +/- ...e-13] + [-1.231604783936 +/- ...e-14]*I
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_elliptic_e(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def elliptic_pi(self, m):
        """
        Return the complete elliptic integral of the third kind evaluated
        at *m* given by ``self``.

        EXAMPLES::

            sage: CBF(2,3).elliptic_pi(CBF(1,1))
            [0.2702999736198...] + [0.715676058329...]*I
        """
        cdef ComplexBall result = self._new()
        cdef ComplexBall my_m = self._parent.coerce(m)
        if _do_sig(prec(self)): sig_on()
        acb_elliptic_pi(result.value, self.value, my_m.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def elliptic_f(self, m):
        r"""
        Return the incomplete elliptic integral of the first kind evaluated
        at *m*.

        See :meth:`elliptic_k` for the corresponding complete integral

        INPUT:

        - ``m`` -- complex ball

        EXAMPLES::

            sage: CBF(1,2).elliptic_f(CBF(0,1))
            [0.6821522911854 +/- ...e-14] + [1.2482780628143 +/- ...e-14]*I

        At parameter `\pi/2` it is a complete integral::

            sage: phi = CBF(1,1)
            sage: (CBF.pi()/2).elliptic_f(phi)
            [1.5092369540513 +/- ...e-14] + [0.6251464152027 +/- ...e-15]*I
            sage: phi.elliptic_k()
            [1.50923695405127 +/- ...e-15] + [0.62514641520270 +/- ...e-15]*I

            sage: phi = CBF(2, 3/7)
            sage: (CBF.pi()/2).elliptic_f(phi)
            [1.3393589639094 +/- ...e-14] + [1.1104369690719 +/- ...e-14]*I
            sage: phi.elliptic_k()
            [1.33935896390938 +/- ...e-15] + [1.11043696907194 +/- ...e-15]*I
        """
        cdef ComplexBall result = self._new()
        cdef ComplexBall my_m = self._parent.coerce(m)
        if _do_sig(prec(self)): sig_on()
        acb_elliptic_f(result.value, self.value, my_m.value, 0, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def elliptic_e_inc(self, m):
        r"""
        Return the incomplete elliptic integral of the second kind evaluated
        at *m*.

        See :meth:`elliptic_e` for the corresponding complete integral

        INPUT:

        - ``m`` -- complex ball

        EXAMPLES::

            sage: CBF(1,2).elliptic_e_inc(CBF(0,1))
            [1.906576998914 +/- ...e-13] + [3.6896645289411 +/- ...e-14]*I

        At parameter `\pi/2` it is a complete integral::

            sage: phi = CBF(1,1)
            sage: (CBF.pi()/2).elliptic_e_inc(phi)
            [1.2838409578982 +/- ...e-14] + [-0.5317843366915 +/- ...e-14]*I
            sage: phi.elliptic_e()
            [1.2838409578982 +/- 5...e-14] + [-0.5317843366915 +/- 3...e-14]*I

            sage: phi = CBF(2, 3/7)
            sage: (CBF.pi()/2).elliptic_e_inc(phi)
            [0.787564350925 +/- ...e-13] + [-0.686896129145 +/- ...e-13]*I
            sage: phi.elliptic_e()
            [0.7875643509254 +/- ...e-14] + [-0.686896129145 +/- ...e-13]*I
        """
        cdef ComplexBall result = self._new()
        cdef ComplexBall my_m = self._parent.coerce(m)
        if _do_sig(prec(self)): sig_on()
        acb_elliptic_e_inc(result.value, self.value, my_m.value, 0, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def elliptic_pi_inc(self, phi, m):
        r"""
        Return the Legendre incomplete elliptic integral of the third kind.

        See: :meth:`elliptic_pi` for the complete integral.

        INPUT:

        - ``phi`` -- complex ball

        - ``m`` -- complex ball

        EXAMPLES::

            sage: CBF(1,2).elliptic_pi_inc(CBF(0,1), CBF(2,-3))
            [0.05738864021418 +/- ...e-15] + [0.55557494549951 +/- ...e-15]*I

        At parameter `\pi/2` it is a complete integral::

            sage: n = CBF(1,1)
            sage: m = CBF(-2/3, 3/5)
            sage: n.elliptic_pi_inc(CBF.pi()/2, m) # this is a regression, see :issue:28623
            nan + nan*I
            sage: n.elliptic_pi(m)
            [0.8934793755173...] + [0.957078687107...]*I

            sage: n = CBF(2, 3/7)
            sage: m = CBF(-1/3, 2/9)
            sage: n.elliptic_pi_inc(CBF.pi()/2, m) # this is a regression, see :issue:28623
            nan + nan*I
            sage: n.elliptic_pi(m)
            [0.296958874641...] + [1.318879533273...]*I
        """
        cdef ComplexBall result = self._new()
        cdef ComplexBall my_phi = self._parent.coerce(phi)
        cdef ComplexBall my_m = self._parent.coerce(m)
        if _do_sig(prec(self)): sig_on()
        acb_elliptic_pi_inc(result.value, self.value, my_phi.value, my_m.value, 0, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def elliptic_rf(self, y, z):
        r"""
        Return the Carlson symmetric elliptic integral of the first kind evaluated
        at ``(self, y, z)``.

        INPUT:

        - ``y`` -- complex ball

        - ``z`` -- complex ball

        EXAMPLES::

            sage: CBF(0,1).elliptic_rf(CBF(-1/2,1), CBF(-1,-1))
            [1.469800396738515 +/- ...e-16] + [-0.2358791199824196 +/- ...e-17]*I
        """
        cdef ComplexBall result = self._new()
        cdef ComplexBall my_y = self._parent.coerce(y)
        cdef ComplexBall my_z = self._parent.coerce(z)
        if _do_sig(prec(self)): sig_on()
        acb_elliptic_rf(result.value, self.value, my_y.value, my_z.value, 0, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def elliptic_rg(self, y, z):
        r"""
        Return the Carlson symmetric elliptic integral of the second kind evaluated
        at ``(self, y, z)``.

        INPUT:

        - ``y`` -- complex ball

        - ``z`` -- complex ball

        EXAMPLES::

            sage: CBF(0,1).elliptic_rg(CBF(-1/2,1), CBF(-1,-1))
            [0.1586786770922370 +/- ...e-17] + [0.2239733128130531 +/- ...e-17]*I
        """
        cdef ComplexBall result = self._new()
        cdef ComplexBall my_y = self._parent.coerce(y)
        cdef ComplexBall my_z = self._parent.coerce(z)
        if _do_sig(prec(self)): sig_on()
        acb_elliptic_rg(result.value, self.value, my_y.value, my_z.value, 0, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def elliptic_rj(self, y, z, p):
        r"""
        Return the Carlson symmetric elliptic integral of the third kind evaluated
        at ``(self, y, z)``.

        INPUT:

        - ``y`` -- complex ball

        - ``z`` -- complex ball

        - ``p`` -- complex bamm

        EXAMPLES::

            sage: CBF(0,1).elliptic_rj(CBF(-1/2,1), CBF(-1,-1), CBF(2))
            [1.00438675628573...] + [-0.24516268343916...]*I
        """
        cdef ComplexBall result = self._new()
        cdef ComplexBall my_y = self._parent.coerce(y)
        cdef ComplexBall my_z = self._parent.coerce(z)
        cdef ComplexBall my_p = self._parent.coerce(p)
        if _do_sig(prec(self)): sig_on()
        acb_elliptic_rj(result.value, self.value, my_y.value, my_z.value, my_p.value, 0, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def elliptic_zeta(self, tau):
        r"""
        Return the value of the Weierstrass zeta function at ``(self, tau)``.

        EXAMPLES::

        - ``tau`` -- a complex ball with positive imaginary part

        EXAMPLES::

            sage: CBF(1,1).elliptic_zeta(CBF(1,3))
            [3.2898676194970 +/- ...e-14] + [0.1365414361782 +/- ...e-14]*I
        """
        cdef ComplexBall result = self._new()
        cdef ComplexBall my_tau = self._parent.coerce(tau)
        if _do_sig(prec(self)): sig_on()
        acb_elliptic_zeta(result.value, self.value, my_tau.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def elliptic_sigma(self, tau):
        r"""
        Return the value of the Weierstrass sigma function at ``(self, tau)``.

        EXAMPLES::

        - ``tau`` -- a complex ball with positive imaginary part

        EXAMPLES::

            sage: CBF(1,1).elliptic_sigma(CBF(1,3))
            [-0.543073363596 +/- ...e-13] + [3.6357291186244 +/- ...e-14]*I
        """
        cdef ComplexBall result = self._new()
        cdef ComplexBall my_tau = self._parent.coerce(tau)
        if _do_sig(prec(self)): sig_on()
        acb_elliptic_sigma(result.value, self.value, my_tau.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def chebyshev_T(self, n):
        """
        Return the Chebyshev function of the first kind of order ``n``
        evaluated at ``self``.

        EXAMPLES::

            sage: CBF(1/3).chebyshev_T(20)
            [0.8710045668809 +/- ...e-14]
            sage: CBF(1/3).chebyshev_T(CBF(5,1))
            [1.84296854518763 +/- ...e-15] + [0.20053614301799 +/- ...e-15]*I
        """
        cdef ComplexBall my_n = self._parent.coerce(n)
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_chebyshev_t(result.value, my_n.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def chebyshev_U(self, n):
        """
        Return the Chebyshev function of the second kind of order ``n``
        evaluated at ``self``.

        EXAMPLES::

            sage: CBF(1/3).chebyshev_U(20)
            [0.6973126541184 +/- ...e-14]
            sage: CBF(1/3).chebyshev_U(CBF(5,1))
            [1.75884964893425 +/- ...e-15] + [0.7497317165104 +/- ...e-14]*I
        """
        cdef ComplexBall my_n = self._parent.coerce(n)
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_chebyshev_u(result.value, my_n.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def jacobi_P(self, n, a, b):
        r"""
        Return the Jacobi polynomial (or function) `P_n^{(a,b)}(z)`
        evaluated at ``self``.

        EXAMPLES::

            sage: CBF(5,-6).jacobi_P(8, CBF(1,2), CBF(2,3))
            [-920983000.45982 +/- ...e-6] + [6069919969.92857 +/- ...e-6]*I
        """
        cdef ComplexBall my_n = self._parent.coerce(n)
        cdef ComplexBall my_a = self._parent.coerce(a)
        cdef ComplexBall my_b = self._parent.coerce(b)
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_jacobi_p(result.value, my_n.value,
            my_a.value, my_b.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def gegenbauer_C(self, n, m):
        r"""
        Return the Gegenbauer polynomial (or function) `C_n^m(z)`
        evaluated at ``self``.

        EXAMPLES::

            sage: CBF(-10).gegenbauer_C(7, 1/2)
            [-263813415.6250000 +/- ...e-8]
        """
        cdef ComplexBall my_n = self._parent.coerce(n)
        cdef ComplexBall my_m = self._parent.coerce(m)
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_gegenbauer_c(result.value, my_n.value,
            my_m.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def laguerre_L(self, n, m=0):
        r"""
        Return the Laguerre polynomial (or function) `L_n^m(z)`
        evaluated at ``self``.

        EXAMPLES::

            sage: CBF(10).laguerre_L(3)
            [-45.6666666666666 +/- ...e-14]
            sage: CBF(10).laguerre_L(3, 2)
            [-6.666666666667 +/- ...e-13]
            sage: CBF(5,7).laguerre_L(CBF(2,3), CBF(1,-2)) # abs tol 1e-9
            [5515.3150302713 +/- 5.02e-11] + [-12386.9428452714 +/- 6.21e-11]*I
        """
        cdef ComplexBall my_n = self._parent.coerce(n)
        cdef ComplexBall my_m = self._parent.coerce(m)
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_laguerre_l(result.value, my_n.value,
            my_m.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def hermite_H(self, n):
        """
        Return the Hermite function (or polynomial) of order ``n``
        evaluated at ``self``.

        EXAMPLES::

            sage: CBF(10).hermite_H(1)
            20.00000000000000
            sage: CBF(10).hermite_H(30)
            [8.0574670961707e+37 +/- ...e+23]
        """
        cdef ComplexBall my_n = self._parent.coerce(n)
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_hermite_h(result.value, my_n.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def legendre_P(self, n, m=0, type=2):
        r"""
        Return the Legendre function of the first kind `P_n^m(z)`
        evaluated at ``self``.

        The ``type`` parameter can be either 2 or 3. This selects between
        different branch cut conventions. The definitions of the "type 2"
        and "type 3" functions are the same as those used by *Mathematica*
        and *mpmath*.

        EXAMPLES::

            sage: CBF(1/2).legendre_P(5)
            [0.0898437500000000 +/- 7...e-17]
            sage: CBF(1,2).legendre_P(CBF(2,3), CBF(0,1))
            [0.10996180744364 +/- ...e-15] + [0.14312767804055 +/- ...e-15]*I
            sage: CBF(-10).legendre_P(5, 325/100)
            [-22104403.487377 +/- ...e-7] + [53364750.687392 +/- ...e-7]*I
            sage: CBF(-10).legendre_P(5, 325/100, type=3)
            [-57761589.914581 +/- ...e-7] + [+/- ...e-7]*I
        """
        cdef ComplexBall my_n = self._parent.coerce(n)
        cdef ComplexBall my_m = self._parent.coerce(m)
        cdef ComplexBall result = self._new()
        cdef int my_type = type
        if my_type != 2 and my_type != 3:
            raise ValueError("expected type 2 or 3")
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_legendre_p(result.value, my_n.value,
            my_m.value, self.value, my_type - 2, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def legendre_Q(self, n, m=0, type=2):
        r"""
        Return the Legendre function of the second kind `Q_n^m(z)`
        evaluated at ``self``.

        The ``type`` parameter can be either 2 or 3. This selects between
        different branch cut conventions. The definitions of the "type 2"
        and "type 3" functions are the same as those used by *Mathematica*
        and *mpmath*.

        EXAMPLES::

            sage: CBF(1/2).legendre_Q(5)
            [0.55508089057168 +/- ...e-15]
            sage: CBF(1,2).legendre_Q(CBF(2,3), CBF(0,1))
            [0.167678710 +/- ...e-10] + [-0.161558598 +/- ...e-10]*I
            sage: CBF(-10).legendre_Q(5, 325/100)
            [-83825154.36008 +/- ...e-6] + [-34721515.80396 +/- ...e-6]*I
            sage: CBF(-10).legendre_Q(5, 325/100, type=3)
            [-4.797306921692e-6 +/- ...e-19] + [-4.797306921692e-6 +/- ...e-19]*I
        """
        cdef ComplexBall my_n = self._parent.coerce(n)
        cdef ComplexBall my_m = self._parent.coerce(m)
        cdef ComplexBall result = self._new()
        cdef int my_type = type
        if my_type != 2 and my_type != 3:
            raise ValueError("expected type 2 or 3")
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_legendre_q(result.value, my_n.value,
            my_m.value, self.value, my_type - 2, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def spherical_harmonic(self, phi, n, m):
        r"""
        Return the spherical harmonic `Y_n^m(\theta,\phi)`
        evaluated at `\theta` given by ``self``.
        In the current implementation, ``n`` and ``m`` must be small integers.

        EXAMPLES::

            sage: CBF(1+I).spherical_harmonic(1/2, -3, -2)
            [0.80370071745224 +/- ...e-15] + [-0.07282031864711 +/- ...e-15]*I
        """
        cdef ComplexBall my_phi = self._parent.coerce(phi)
        cdef long my_n = n
        cdef long my_m = m
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_spherical_y(result.value, my_n, my_m,
            self.value, my_phi.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result


CBF = ComplexBallField()
