# -*- coding: utf-8
r"""
Arbitrary precision complex balls using Arb

This is a binding to the `Arb library <http://fredrikj.net/arb/>`_; it
may be useful to refer to its documentation for more details.

Parts of the documentation for this module are copied or adapted from
Arb's own documentation, licenced under the GNU General Public License
version 2, or later.

.. SEEALSO::

    - :mod:`Real balls using Arb <sage.rings.real_arb>`
    - :mod:`Complex interval field (using MPFI) <sage.rings.complex_interval_field>`
    - :mod:`Complex intervals (using MPFI) <sage.rings.complex_interval>`

Data Structure
==============

A :class:`ComplexBall` represents a complex number with error bounds. It wraps
an Arb object of type ``acb_t``, which  consists of a pair of real number balls
representing the real and imaginary part with separate error bounds. (See the
documentation of :mod:`sage.rings.real_arb` for more information.)

A :class:`ComplexBall` thus represents a rectangle `[m_1-r_1, m_1+r_1] +
[m_2-r_2, m_2+r_2] i` in the complex plane. This is used in Arb instead of a
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

    In accordance with the semantics of Arb, identical :class:`ComplexBall`
    objects are understood to give permission for algebraic simplification.
    This assumption is made to improve performance. For example, setting ``z =
    x*x`` sets `z` to a ball enclosing the set `\{t^2 : t \in x\}` and not the
    (generally larger) set `\{tu : t \in x, u \in x\}`.

Two elements are equal if and only if they are the same object
or if both are exact and equal::

    sage: a = CBF(1, 2)
    sage: b = CBF(1, 2)
    sage: a is b
    False
    sage: a == b
    True
    sage: a = CBF(1/3, 1/5)
    sage: b = CBF(1/3, 1/5)
    sage: a.is_exact()
    False
    sage: b.is_exact()
    False
    sage: a is b
    False
    sage: a == b
    False

A ball is non-zero in the sense of usual comparison if and only if it does not
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

    sage: bpol = 1/3*CBF(i) + AA(sqrt(2)) + (polygen(RealBallField(20), 'x') + QQbar(i))
    sage: bpol
    x + [1.41421 +/- 5.09e-6] + [1.33333 +/- 3.97e-6]*I
    sage: bpol.parent()
    Univariate Polynomial Ring in x over Complex ball field with 20 bits precision
    sage: bpol/3
    ([0.333333 +/- 4.93e-7])*x + [0.47140 +/- 5.39e-6] + [0.44444 +/- 4.98e-6]*I

TESTS::

    sage: polygen(CBF, x)^3
    x^3

::

    sage: SR.coerce(CBF(0.42 + 3.33*I))
    [0.4200000000000000 +/- 1.56e-17] + [3.330000000000000 +/- 7.11e-17]*I

Check that :trac:`19839` is fixed::

    sage: log(SR(CBF(0.42))).pyobject().parent()
    Complex ball field with 53 bits precision

Classes and Methods
===================
"""
#*****************************************************************************
# Copyright (C) 2014 Clemens Heuberger <clemens.heuberger@aau.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                http://www.gnu.org/licenses/
#*****************************************************************************
include "cysignals/signals.pxi"

import operator

import sage.categories.fields

cimport sage.rings.integer
cimport sage.rings.rational

import sage.rings.number_field.number_field as number_field

from cpython.float cimport PyFloat_AS_DOUBLE
from cpython.int cimport PyInt_AS_LONG
from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_NE, Py_GT, Py_GE

from sage.libs.mpfr cimport MPFR_RNDU
from sage.libs.arb.arb cimport *
from sage.libs.arb.acb cimport *
from sage.libs.arb.acb_hypgeom cimport *
from sage.libs.arb.acb_modular cimport *
from sage.libs.arb.arf cimport arf_init, arf_get_mpfr, arf_set_mpfr, arf_clear, arf_set_mag, arf_set
from sage.libs.arb.mag cimport mag_init, mag_clear, mag_add, mag_set_d, MAG_BITS, mag_is_inf, mag_is_finite, mag_zero
from sage.libs.flint.fmpz cimport fmpz_t, fmpz_init, fmpz_get_mpz, fmpz_set_mpz, fmpz_clear, fmpz_abs
from sage.libs.flint.fmpq cimport fmpq_t, fmpq_init, fmpq_set_mpq, fmpq_clear
from sage.libs.gmp.mpz cimport mpz_fits_ulong_p, mpz_fits_slong_p, mpz_get_ui, mpz_get_si, mpz_sgn
from sage.rings.complex_field import ComplexField
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field_element_quadratic cimport NumberFieldElement_quadratic
from sage.rings.real_arb cimport mpfi_to_arb, arb_to_mpfi, real_part_of_quadratic_element_to_arb
from sage.rings.real_arb import RealBallField
from sage.rings.real_mpfr cimport RealField_class, RealField, RealNumber
from sage.rings.ring import Field
from sage.structure.element cimport Element, ModuleElement
from sage.structure.parent cimport Parent
from sage.structure.unique_representation import UniqueRepresentation

cdef void ComplexIntervalFieldElement_to_acb(
    acb_t target,
    ComplexIntervalFieldElement source):
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


class ComplexBallField(UniqueRepresentation, Field):
    r"""
    An approximation of the field of complex numbers using pairs of mid-rad
    intervals.

    INPUT:

    - ``precision`` -- an integer `\ge 2`.

    EXAMPLES::

        sage: CBF(1)
        1.000000000000000

    TESTS::

        sage: ComplexBallField(0)
        Traceback (most recent call last):
        ...
        ValueError: Precision must be at least 2.
        sage: ComplexBallField(1)
        Traceback (most recent call last):
        ...
        ValueError: Precision must be at least 2.
    """
    Element = ComplexBall

    @staticmethod
    def __classcall__(cls, long precision=53, category=None):
        r"""
        Normalize the arguments for caching.

        TESTS::

            sage: ComplexBallField(53) is ComplexBallField()
            True
        """
        return super(ComplexBallField, cls).__classcall__(cls, precision, category)

    def __init__(self, precision, category):
        r"""
        Initialize the complex ball field.

        INPUT:

        - ``precision`` -- an integer `\ge 2`.

        EXAMPLES::

            sage: CBF(1)
            1.000000000000000

        TESTS::

            sage: CBF.base()
            Real ball field with 53 bits precision
            sage: CBF.base_ring()
            Real ball field with 53 bits precision

        There are direct coercions from ZZ and QQ (for which arb provides
        construction functions)::

            sage: CBF.coerce_map_from(ZZ)
            Conversion map:
            From: Integer Ring
            To:   Complex ball field with 53 bits precision
            sage: CBF.coerce_map_from(QQ)
            Conversion map:
            From: Rational Field
            To:   Complex ball field with 53 bits precision

        Various other coercions are available through real ball fields or CLF::

            sage: CBF.coerce_map_from(RLF)
            Composite map:
            From: Real Lazy Field
            To:   Complex ball field with 53 bits precision
            Defn:   Conversion map:
                    From: Real Lazy Field
                    To:   Real ball field with 53 bits precision
                    then
                    Conversion map:
                    From: Real ball field with 53 bits precision
                    To:   Complex ball field with 53 bits precision
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
            raise ValueError("Precision must be at least 2.")
        real_field = RealBallField(precision)
        super(ComplexBallField, self).__init__(
                base_ring=real_field,
                category=category or sage.categories.fields.Fields().Infinite())
        self._prec = precision
        from sage.rings.rational_field import QQ
        from sage.rings.real_lazy import CLF
        self._populate_coercion_lists_([ZZ, QQ, real_field, CLF])

    def _real_field(self):
        """
        TESTS::

            sage: CBF._real_field()
            Real ball field with 53 bits precision
        """
        return self._base

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: ComplexBallField()
            Complex ball field with 53 bits precision
            sage: ComplexBallField(106)
            Complex ball field with 106 bits precision
        """
        return "Complex ball field with {} bits precision".format(self._prec)

    def construction(self):
        """
        Return the construction of a complex ball field as the algebraic
        closure of the real ball field with the same precision.

        EXAMPLES::

            sage: functor, base = CBF.construction()
            sage: functor, base
            (AlgebraicClosureFunctor, Real ball field with 53 bits precision)
            sage: functor(base) is CBF
            True
        """
        from sage.categories.pushout import AlgebraicClosureFunctor
        return (AlgebraicClosureFunctor(), self._base)

    def complex_field(self):
        """
        Return the complex ball field with the same precision, i.e. ``self``

        EXAMPLES::

            sage: CBF.complex_field() is CBF
            True
        """
        return ComplexBallField(self._prec)

    def ngens(self):
        r"""
        Return 1 as the only generator is the imaginary unit.

        EXAMPLE::

            sage: CBF.ngens()
            1
        """
        return 1

    def gen(self, i):
        r"""
        For i = 0, return the imaginary unit in this complex ball field.

        EXAMPLE::

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

        EXAMPLE::

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
            Identity endomorphism of Complex ball field with 53 bits precision
            sage: CBF.coerce_map_from(ComplexBallField(100))
            Conversion map:
            From: Complex ball field with 100 bits precision
            To:   Complex ball field with 53 bits precision
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
        """
        if isinstance(other, (RealBallField, ComplexBallField)):
            return (other._prec >= self._prec)
        elif isinstance(other, number_field.NumberField_quadratic):
            emb = other.coerce_embedding()
            if emb is not None:
                return self.has_coerce_map_from(emb.codomain())

    def _element_constructor_(self, x=None, y=None):
        r"""
        Convert (x, y) to an element of this complex ball field, perhaps
        non-canonically.

        INPUT:

        - ``x``, ``y`` (optional) -- either a complex number, interval or ball,
          or two real ones (see examples below for more information on accepted
          number types).

        EXAMPLES::

            sage: CBF()
            0
            sage: CBF(1) # indirect doctest
            1.000000000000000
            sage: CBF(1, 1)
            1.000000000000000 + 1.000000000000000*I
            sage: CBF(pi, sqrt(2))
            [3.141592653589793 +/- 5.61e-16] + [1.414213562373095 +/- 4.10e-16]*I
            sage: CBF(I)
            1.000000000000000*I
            sage: CBF(pi+I/3)
            [3.141592653589793 +/- 5.61e-16] + [0.3333333333333333 +/- 7.04e-17]*I
            sage: CBF(QQbar(i/7))
            [0.1428571428571428 +/- 9.09e-17]*I
            sage: CBF(AA(sqrt(2)))
            [1.414213562373095 +/- 4.10e-16]
            sage: CBF(CIF(0, 1))
            1.000000000000000*I
            sage: CBF(RBF(1/3))
            [0.3333333333333333 +/- 7.04e-17]
            sage: CBF(RBF(1/3), RBF(1/6))
            [0.3333333333333333 +/- 7.04e-17] + [0.1666666666666667 +/- 7.04e-17]*I
            sage: CBF(1/3)
            [0.3333333333333333 +/- 7.04e-17]
            sage: CBF(1/3, 1/6)
            [0.3333333333333333 +/- 7.04e-17] + [0.1666666666666667 +/- 7.04e-17]*I
            sage: ComplexBallField(106)(1/3, 1/6)
            [0.33333333333333333333333333333333 +/- 6.94e-33] + [0.16666666666666666666666666666666 +/- 7.70e-33]*I
            sage: NF.<a> = QuadraticField(-2)
            sage: CBF(1/5 + a/2)
            [0.2000000000000000 +/- 4.45e-17] + [0.707106781186547 +/- 5.86e-16]*I
            sage: CBF(infinity, NaN)
            [+/- inf] + nan*I
            sage: CBF(x)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert x to a ComplexBall

        .. SEEALSO::

            :meth:`sage.rings.real_arb.RealBallField._element_constructor_`

        TESTS::

            sage: CBF(1+I, 2)
            Traceback (most recent call last):
            ...
            ValueError: nonzero imaginary part
        """
        try:
            return self.element_class(self, x, y)
        except TypeError:
            pass

        if y is None:
            try:
                x = self._base(x)
                return self.element_class(self, x)
            except (TypeError, ValueError):
                pass
            try:
                y = self._base(x.imag())
                x = self._base(x.real())
                return self.element_class(self, x, y)
            except (AttributeError, TypeError):
                pass
            try:
                x = ComplexIntervalField(self._prec)(x)
                return self.element_class(self, x)
            except TypeError:
                pass
            raise TypeError("unable to convert {!r} to a ComplexBall".format(x))
        else:
            x = self._base(x)
            y = self._base(y)
            return self.element_class(self, x, y)

    def _an_element_(self):
        r"""
        Construct an element.

        EXAMPLES::

            sage: CBF.an_element() # indirect doctest
            [0.3333333333333333 +/- 1.49e-17] - [0.1666666666666667 +/- 4.26e-17]*I
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

    def is_exact(self):
        """
        Complex ball fields are not exact.

        EXAMPLES::

            sage: ComplexBallField().is_exact()
            False
        """
        return False

    def is_finite(self):
        """
        Complex ball fields are infinite.

        They already specify it via their category, but we currently need to
        re-implement this method due to the legacy implementation in
        :class:`sage.rings.ring.Ring`.

        EXAMPLES::

            sage: ComplexBallField().is_finite()
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
             1.000000000000000 + [0.3333333333333333 +/- 1.49e-17]*I,
             [-0.3333333333333333 +/- 1.49e-17] + 0.2500000000000000*I,
             [-2.175556475109056e+181961467118333366510562 +/- 1.29e+181961467118333366510545],
             [+/- inf],
             [0.3333333333333333 +/- 1.49e-17] + [+/- inf]*I,
             [+/- inf] + [+/- inf]*I,
             nan,
             nan + nan*I,
             [+/- inf] + nan*I]
        """
        return [self(1), self(0, -1./2), self(1, 1./3), self(-1./3, 1./4),
                -self(1, 1)**(sage.rings.integer.Integer(2)**80),
                self('inf'), self(1./3, 'inf'), self('inf', 'inf'),
                self('nan'), self('nan', 'nan'), self('inf', 'nan')]

cdef inline bint _do_sig(long prec):
    """
    Whether signal handlers should be installed for calls to arb.
    """
    return (prec > 1000)

cdef inline long prec(ComplexBall ball):
    return ball._parent._prec

cdef inline Parent real_ball_field(ComplexBall ball):
    return ball._parent._base

cdef class ComplexBall(RingElement):
    """
    Hold one ``acb_t`` of the `Arb library
    <http://fredrikj.net/arb/>`_

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

        - ``parent`` -- a :class:`ComplexBallField`.

        - ``x``, ``y`` (optional) -- either a complex number, interval or ball,
          or two real ones.

        .. SEEALSO:: :meth:`ComplexBallField._element_constructor_`

        TESTS::

            sage: from sage.rings.complex_arb import ComplexBall
            sage: CBF53, CBF100 = ComplexBallField(53), ComplexBallField(100)
            sage: ComplexBall(CBF100)
            0
            sage: ComplexBall(CBF100, ComplexBall(CBF53, ComplexBall(CBF100, 1/3)))
            [0.333333333333333333333333333333 +/- 4.65e-31]
            sage: ComplexBall(CBF100, RBF(pi))
            [3.141592653589793 +/- 5.61e-16]
            sage: ComplexBall(CBF100, -3r)
            -3.000000000000000000000000000000
            sage: ComplexBall(CBF100, 10^100)
            1.000000000000000000000000000000e+100
            sage: ComplexBall(CBF100, CIF(1, 2))
            1.000000000000000000000000000000 + 2.000000000000000000000000000000*I
            sage: ComplexBall(CBF100, RBF(1/3), RBF(1))
            [0.3333333333333333 +/- 7.04e-17] + 1.000000000000000000000000000000*I
            sage: ComplexBall(CBF100, 1, 2)
            Traceback (most recent call last):
            ...
            TypeError: unsupported initializer
            sage: NF.<a> = QuadraticField(-1, embedding=CC(0, -1))
            sage: CBF(a)
            -1.000000000000000*I
            sage: NF.<a> = QuadraticField(-1, embedding=None)
            sage: CBF(a)
            Traceback (most recent call last):
            ...
            ValueError: need an embedding
        """
        cdef fmpz_t tmpz
        cdef fmpq_t tmpq
        cdef NumberFieldElement_quadratic x_as_qe
        cdef long myprec

        RingElement.__init__(self, parent)

        if x is None:
            return
        elif y is None:
            if isinstance(x, ComplexBall):
                acb_set(self.value, (<ComplexBall> x).value)
            elif isinstance(x, RealBall):
                acb_set_arb(self.value, (<RealBall> x).value)
            elif isinstance(x, int):
                acb_set_si(self.value, PyInt_AS_LONG(x))
            elif isinstance(x, sage.rings.integer.Integer):
                if _do_sig(prec(self)): sig_on()
                fmpz_init(tmpz)
                fmpz_set_mpz(tmpz, (<sage.rings.integer.Integer> x).value)
                acb_set_fmpz(self.value, tmpz)
                fmpz_clear(tmpz)
                if _do_sig(prec(self)): sig_off()
            elif isinstance(x, sage.rings.rational.Rational):
                if _do_sig(prec(self)): sig_on()
                fmpq_init(tmpq)
                fmpq_set_mpq(tmpq, (<sage.rings.rational.Rational> x).value)
                acb_set_fmpq(self.value, tmpq, prec(self))
                fmpq_clear(tmpq)
                if _do_sig(prec(self)): sig_off()
            elif isinstance(x, ComplexIntervalFieldElement):
                ComplexIntervalFieldElement_to_acb(self.value,
                                                   <ComplexIntervalFieldElement> x)
            elif isinstance(x, NumberFieldElement_quadratic):
                x_as_qe = <NumberFieldElement_quadratic> x
                real_part_of_quadratic_element_to_arb(acb_realref(self.value),
                        x_as_qe, prec(self))
                myprec = prec(self) + 4
                if mpz_sgn(x_as_qe.D.value) < 0:
                    if x_as_qe._parent._embedding is None:
                        raise ValueError("need an embedding")
                    fmpz_init(tmpz)
                    fmpz_set_mpz(tmpz, x_as_qe.D.value)
                    fmpz_abs(tmpz, tmpz)
                    arb_sqrt_fmpz(acb_imagref(self.value), tmpz, myprec)
                    fmpz_set_mpz(tmpz, x_as_qe.b)
                    arb_mul_fmpz(acb_imagref(self.value), acb_imagref(self.value), tmpz, myprec)
                    fmpz_set_mpz(tmpz, x_as_qe.denom)
                    arb_div_fmpz(acb_imagref(self.value), acb_imagref(self.value), tmpz, prec(self))
                    fmpz_clear(tmpz)
                    if not x_as_qe.standard_embedding:
                        acb_conj(self.value, self.value)
            else:
                raise TypeError("unsupported initializer")
        elif isinstance(x, RealBall) and isinstance(y, RealBall):
            arb_set(acb_realref(self.value), (<RealBall> x).value)
            arb_set(acb_imagref(self.value), (<RealBall> y).value)
        else:
            raise TypeError("unsupported initializer")

    cdef ComplexBall _new(self):
        """
        Return a new complex ball element with the same parent as ``self``.
        """
        cdef ComplexBall x
        x = ComplexBall.__new__(ComplexBall)
        x._parent = self._parent
        return x

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

        OUTPUT:

        A string.

        EXAMPLES::

           sage: CBF(1/3)
           [0.3333333333333333 +/- 7.04e-17]
           sage: CBF(0, 1/3)
           [0.3333333333333333 +/- 7.04e-17]*I
           sage: CBF(1/3, 1/6)
           [0.3333333333333333 +/- 7.04e-17] + [0.1666666666666667 +/- 7.04e-17]*I

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

    def _is_atomic(self):
        r"""
        Declare that complex balls print atomically in some cases.

        TESTS::

            sage: CBF(-1/3)._is_atomic()
            True

        This method should in principle ensure that ``CBF['x']([1, -1/3])``
        is printed as::

            sage: CBF['x']([1, -1/3]) # todo - not tested
            [-0.3333333333333333 +/- 7.04e-17]*x + 1.000000000000000

        However, this facility is not really used in Sage at this point, and we
        still get::

            sage: CBF['x']([1, -1/3])
            ([-0.3333333333333333 +/- 7.04e-17])*x + 1.000000000000000
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
        cdef sage.rings.integer.Integer res
        cdef fmpz_t tmp
        fmpz_init(tmp)
        try:
            if acb_get_unique_fmpz(tmp, self.value):
                res = sage.rings.integer.Integer.__new__(sage.rings.integer.Integer)
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

        - ``parent`` - :class:`~sage.rings.complex_field.ComplexField_class`,
          target parent.

        EXAMPLES::

            sage: CC(CBF(1/3, 1/3))
            0.333333333333333 + 0.333333333333333*I
            sage: ComplexField(100)(CBF(1/3, 1/3))
            0.33333333333333331482961625625 + 0.33333333333333331482961625625*I
        """
        real_field = parent._base
        return parent(real_field(self.real()), real_field(self.imag()))

    def _real_mpfi_(self, parent):
        r"""
        Try to convert this complex ball to a real interval.

        Fail if the imaginary part is not exactly zero.

        INPUT:

        - ``parent`` - :class:`~sage.rings.real_mpfi.RealIntervalField_class`,
          target parent.

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

        - ``parent`` - :class:`~sage.rings.real_mpfr.RealField_class`,
          target parent.

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

    # Real and imaginary part, midpoint, radius

    cpdef RealBall real(self):
        """
        Return the real part of this ball.

        OUTPUT:

        A :class:`~sage.rings.real_arb.RealBall`.

        EXAMPLES::

           sage: a = CBF(1/3, 1/5)
           sage: a.real()
           [0.3333333333333333 +/- 7.04e-17]
        """
        cdef RealBall r = RealBall(real_ball_field(self))
        arb_set(r.value, acb_realref(self.value))
        return r

    cpdef RealBall imag(self):
        """
        Return the imaginary part of this ball.

        OUTPUT:

        A :class:`~sage.rings.real_arb.RealBall`.

        EXAMPLES::

           sage: a = CBF(1/3, 1/5)
           sage: a.imag()
           [0.2000000000000000 +/- 4.45e-17]
        """
        cdef RealBall r = RealBall(real_ball_field(self))
        arb_set(r.value, acb_imagref(self.value))
        return r

    def __abs__(self):
        """
        Return the absolute value of this complex ball.

        EXAMPLES::

            sage: CBF(1 + i).abs() # indirect doctest
            [1.414213562373095 +/- 2.99e-16]
            sage: abs(CBF(i))
            1.000000000000000

            sage: CBF(1 + i).abs().parent()
            Real ball field with 53 bits precision
        """
        cdef RealBall r = RealBall(real_ball_field(self))
        acb_abs(r.value, self.value, prec(self))
        return r

    def below_abs(self, test_zero=False):
        """
        Return a lower bound for the absolute value of this complex ball.

        INPUT:

        - ``test_zero`` (boolean, default ``False``) -- if ``True``,
          make sure that the returned lower bound is positive, raising
          an error if the ball contains zero.

        OUTPUT:

        A ball with zero radius

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

        OUTPUT:

        A ball with zero radius

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
            [0.785398163397448 +/- 3.91e-16]
            sage: CBF(-1).arg()
            [3.141592653589793 +/- 5.61e-16]
            sage: CBF(-1).arg().parent()
            Real ball field with 53 bits precision
        """
        cdef RealBall r = RealBall(real_ball_field(self))
        acb_arg(r.value, self.value, prec(self))
        return r

    def mid(self):
        """
        Return the midpoint of this ball.

        OUTPUT:

        :class:`~sage.rings.complex_number.ComplexNumber`, floating-point
        complex number formed by the centers of the real and imaginary parts of
        this ball.

        EXAMPLES::

            sage: CBF(1/3, 1).mid()
            0.333333333333333 + 1.00000000000000*I
            sage: CBF(1/3, 1).mid().parent()
            Complex Field with 53 bits of precision
            sage: CBF('inf', 'nan').mid()
            +infinity - NaN*I
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
        field = sage.rings.complex_field.ComplexField(
                max(prec(self), re.prec(), im.prec()))
        return field(re, im)

    def squash(self):
        """
        Return an exact ball with the same midpoint as this ball.

        OUTPUT:

        A :class:`ComplexBall`.

        EXAMPLES::

            sage: mid = CBF(1/3, 1/10).squash()
            sage: mid
            [0.3333333333333333 +/- 1.49e-17] + [0.09999999999999999 +/- 1.68e-18]*I
            sage: mid.parent()
            Complex ball field with 53 bits precision
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

    # Precision and accuracy

    def round(self):
        """
        Return a copy of this ball rounded to the precision of the parent.

        EXAMPLES:

        It is possible to create balls whose midpoint is more precise that
        their parent's nominal precision (see :mod:`~sage.rings.real_arb` for
        more information)::

            sage: b = CBF(exp(I*pi/3).n(100))
            sage: b.mid()
            0.50000000000000000000000000000 + 0.86602540378443864676372317075*I

        The ``round()`` method rounds such a ball to its parent's precision::

            sage: b.round().mid()
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

            sage: CBF(exp(I*pi/3)).accuracy()
            52
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

        - ``ampl`` - A **real** ball (or an object that can be coerced to a
          real ball).

        OUTPUT:

        A new complex ball.

        EXAMPLES::

            sage: CBF(1+i).add_error(10^-16)
            [1.000000000000000 +/- 1.01e-16] + [1.000000000000000 +/- 1.01e-16]*I
        """
        return ComplexBall(self._parent, self.real().add_error(ampl), self.imag().add_error(ampl))

    # Comparisons and predicates

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

            sage: CBF(pi, 1/3).is_nonzero()
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

    def __nonzero__(self):
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
        cdef acb_t difference

        lt = left
        rt = right

        if op == Py_EQ:
            return lt is rt or (
                acb_is_exact(lt.value) and acb_is_exact(rt.value)
                and acb_equal(lt.value, rt.value))

        if op == Py_NE:
            return not acb_overlaps(lt.value, rt.value)

        elif op == Py_GT or op == Py_GE or op == Py_LT or op == Py_LE:
            raise TypeError("No order is defined for ComplexBalls.")

    def identical(self, ComplexBall other):
        """
        Return whether ``self`` and ``other`` represent the same ball.

        INPUT:

        - ``other`` -- a :class:`ComplexBall`.

        OUTPUT:

        Return True iff ``self`` and ``other`` are equal as sets, i.e. if their
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
        Return True iff ``self`` and ``other`` have some point in common.

        INPUT:

        - ``other`` -- a :class:`ComplexBall`.

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

            sage: CBF(sqrt(2)).contains_exact(sqrt(2))
            Traceback (most recent call last):
            ...
            TypeError: unsupported type: <type 'sage.symbolic.expression.Expression'>
        """
        cdef fmpz_t tmpz
        cdef fmpq_t tmpq
        if _do_sig(prec(self)): sig_on()
        try:
            if isinstance(other, ComplexBall):
                res = acb_contains(self.value, (<ComplexBall> other).value)
            elif isinstance(other, sage.rings.integer.Integer):
                fmpz_init(tmpz)
                fmpz_set_mpz(tmpz, (<sage.rings.integer.Integer> other).value)
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
        Return True if ``other`` can be verified to be contained in ``self``.

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
                sage.rings.integer.Integer,
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
            [-0.3333333333333333 +/- 7.04e-17] - 1.000000000000000*I
        """
        cdef ComplexBall res = self._new()
        acb_neg(res.value, self.value)
        return res

    def conjugate(self):
        """
        Return the complex conjugate of this ball.

        EXAMPLES::

            sage: CBF(-2 + I/3).conjugate()
            -2.000000000000000 + [-0.3333333333333333 +/- 7.04e-17]*I
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
            [-3.00000000000000 +/- 9.44e-16]*I
            sage: ~CBF(0)
            [+/- inf]
            sage: ~CBF(RIF(10,11))
            [0.1 +/- 9.53e-3]
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
            [-2.333333333333333 +/- 5.37e-16] + [0.333333333333333 +/- 4.82e-16]*I
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

        - ``shift`` -- integer, may be negative.

        EXAMPLES::

            sage: CBF(i/3) << 2
            [1.333333333333333 +/- 4.82e-16]*I
            sage: CBF(i) << -2
            0.2500000000000000*I

        TESTS::

            sage: CBF(i) << (2^65)
            [3.636549880934858e+11106046577046714264 +/- 1.91e+11106046577046714248]*I
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
        if isinstance(shift, int):
             acb_mul_2exp_si(res.value, self.value, PyInt_AS_LONG(shift))
        elif isinstance(shift, sage.rings.integer.Integer):
            sig_on()
            fmpz_init(tmpz)
            fmpz_set_mpz(tmpz, (<sage.rings.integer.Integer> shift).value)
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

        - ``shift`` -- integer, may be negative.

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
            [-1.50000000000000 +/- 1.27e-15] + [1.500000000000000 +/- 8.94e-16]*I
            sage: CBF(2+I)/CBF(0)
            [+/- inf] + [+/- inf]*I
            sage: CBF(1)/CBF(0)
            [+/- inf]
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
            sage: CBF(e)**CBF(i*pi)
            [-1.0000000000000 +/- 1.98e-15] + [+/- 2.32e-15]*I
            sage: CBF(0, 1)**AA(2)**(1/2)
            [-0.60569986707881 +/- 4.36e-15] + [0.79569320156748 +/- 2.53e-15]*I

            sage: CBF(i)**RBF(2**1000)
            [+/- 2.51] + [+/- 2.87]*I
            sage: CBF(i)**(2**1000)
            1.000000000000000

            sage: CBF(0)^(1/3)
            0
            sage: CBF(0)^(-1)
            [+/- inf]
            sage: CBF(0)^(-2)
            [+/- inf] + [+/- inf]*I

        TESTS::

            sage: (CBF(e)**CBF(i))**RBF(pi)
            [-1.0000000000000 +/- 5.48e-15] + [+/- 4.14e-15]*I
            sage: CBF(2*i)**10r
            -1024.000000000000
            sage: CBF(1,1) ^ -1r
            0.5000000000000000 - 0.5000000000000000*I

        """
        cdef fmpz_t tmpz
        if not isinstance(base, ComplexBall):
            return sage.structure.element.bin_op(base, expo, operator.pow)
        cdef ComplexBall self = base
        cdef ComplexBall res = self._new()
        if isinstance(expo, int):
            if _do_sig(prec(self)): sig_on()
            acb_pow_si(res.value, self.value, PyInt_AS_LONG(expo), prec(self))
            if _do_sig(prec(self)): sig_off()
        elif isinstance(expo, sage.rings.integer.Integer):
            if _do_sig(prec(self)): sig_on()
            fmpz_init(tmpz)
            fmpz_set_mpz(tmpz, (<sage.rings.integer.Integer> expo).value)
            acb_pow_fmpz(res.value, self.value, tmpz, prec(self))
            fmpz_clear(tmpz)
            if _do_sig(prec(self)): sig_off()
        elif isinstance(expo, ComplexBall):
            if _do_sig(prec(self)): sig_on()
            acb_pow(res.value, self.value, (<ComplexBall> expo).value, prec(self))
            if _do_sig(prec(self)): sig_off()
        elif isinstance(expo, RealBall):
            if _do_sig(prec(self)): sig_on()
            acb_pow_arb(res.value, self.value, (<RealBall> expo).value, prec(self))
            if _do_sig(prec(self)): sig_off()
        else:
            return sage.structure.element.bin_op(base, expo, operator.pow)
        return res

    def sqrt(self):
        """
        Return the square root of this ball.

        If either the real or imaginary part is exactly zero, only a single
        real square root is needed.

        EXAMPLES::

            sage: CBF(-2).sqrt()
            [1.414213562373095 +/- 2.99e-16]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_sqrt(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def rsqrt(self):
        """
        Return the reciprocal square root of ``self``.

        If either the real or imaginary part is exactly zero, only a single
        real reciprocal square root is needed.

        EXAMPLES::

            sage: CBF(-2).rsqrt()
            [-0.707106781186547 +/- 5.73e-16]*I
            sage: CBF(0, 1/2).rsqrt()
            1.000000000000000 - 1.000000000000000*I
            sage: CBF(0).rsqrt()
            nan
        """
        cdef ComplexBall res = self._new()
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
        """
        Return the ``n``-th rising factorial of this ball.

        The `n`-th rising factorial of `x` is equal to `x (x+1) \cdots (x+n-1)`.

        For complex `n`, it is a quotient of gamma functions.

        EXAMPLES::

            sage: CBF(1).rising_factorial(5)
            120.0000000000000
            sage: CBF(1/3, 1/2).rising_factorial(300)
            [-3.87949484514e+612 +/- 5.24e+600] + [-3.52042209763e+612 +/- 5.56e+600]*I

            sage: CBF(1).rising_factorial(-1)
            nan
            sage: CBF(1).rising_factorial(2**64)
            [+/- 2.30e+347382171305201370464]
            sage: ComplexBallField(128)(1).rising_factorial(2**64)
            [2.343691126796861348e+347382171305201285713 +/- 4.71e+347382171305201285694]
            sage: CBF(1/2).rising_factorial(CBF(2,3))
            [-0.123060451458124 +/- 4.46e-16] + [0.040641263167655 +/- 3.75e-16]*I

        """
        cdef ComplexBall result = self._new()
        cdef ComplexBall my_n = self._parent.coerce(n)
        if _do_sig(prec(self)): sig_on()
        acb_rising(result.value, self.value, my_n.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    # Elementary functions

    def log(self, base=None):
        """
        General logarithm (principal branch).

        INPUT:

        - ``base`` (optional, complex ball or number) -- if ``None``, return
          the principal branch of the natural logarithm ``ln(self)``,
          otherwise, return the general logarithm ``ln(self)/ln(base)``

        EXAMPLES::

            sage: CBF(2*i).log()
            [0.6931471805599453 +/- 4.16e-17] + [1.570796326794897 +/- 6.65e-16]*I
            sage: CBF(-1).log()
            [3.141592653589793 +/- 5.61e-16]*I

            sage: CBF(2*i).log(2)
            [1.000000000000000 +/- 8.01e-17] + [2.26618007091360 +/- 4.23e-15]*I
            sage: CBF(2*i).log(CBF(i))
            [1.000000000000000 +/- 2.83e-16] + [-0.441271200305303 +/- 2.82e-16]*I

            sage: CBF('inf').log()
            nan + nan*I
            sage: CBF(2).log(0)
            nan + nan*I
        """
        cdef ComplexBall cst
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_log(res.value, self.value, prec(self))
        if base is not None:
            cst = self._parent.coerce(base).log()
            if _do_sig(prec(self)): sig_on()
            acb_div(res.value, res.value, cst.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        if _do_sig(prec(self)): sig_off()
        return res

    def log1p(self):
        """
        Return ``log(1 + self)``, computed accurately when ``self`` is close to
        zero.

        EXAMPLES::

            sage: eps = RBF(1e-50)
            sage: CBF(1+eps, eps).log()
            [+/- 2.23e-16] + [1.000000000000000e-50 +/- 2.30e-66]*I
            sage: CBF(eps, eps).log1p()
            [1.000000000000000e-50 +/- 7.63e-68] + [1.00000000000000e-50 +/- 2.30e-66]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_log1p(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def exp(self):
        """
        Return the exponential of this ball.

        .. SEEALSO:: :meth:`exppii`

        EXAMPLES::

            sage: CBF(i*pi).exp()
            [-1.00000000000000 +/- 6.67e-16] + [+/- 5.68e-16]*I
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
            sage: CBF(0, -1/pi).exppii()
            [2.71828182845904 +/- 6.20e-15]
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

            sage: CBF(i*pi).sin()
            [11.5487393572577 +/- 5.34e-14]*I
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

            sage: CBF(i*pi).cos()
            [11.59195327552152 +/- 8.38e-15]
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

            sage: CBF(pi/2, 1/10).tan()
            [+/- 2.87e-14] + [10.0333111322540 +/- 3.16e-14]*I
            sage: CBF(pi/2).tan()
            [+/- inf]
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

            sage: CBF(pi, 1/10).cot()
            [+/- 5.74e-14] + [-10.0333111322540 +/- 4.05e-14]*I
            sage: CBF(pi).cot()
            [+/- inf]
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_cot(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arcsin(self):
        """
        Return the arcsine of this ball.

        EXAMPLES::

            sage: CBF(1+i).arcsin()
            [0.66623943249252 +/- 5.40e-15] + [1.06127506190504 +/- 5.04e-15]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_asin(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arccos(self):
        """
        Return the arccosine of this ball.

        EXAMPLES::

            sage: CBF(1+i).arccos()
            [0.90455689430238 +/- 2.18e-15] + [-1.06127506190504 +/- 5.04e-15]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_acos(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arctan(self):
        """
        Return the arctangent of this ball.

        EXAMPLES::

            sage: CBF(1+i).arctan()
            [1.017221967897851 +/- 4.93e-16] + [0.4023594781085251 +/- 8.52e-17]*I
            sage: CBF(i).arctan()
            nan + nan*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_atan(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arcsinh(self):
        """
        Return the hyperbolic arcsine of this ball.

        EXAMPLES::

            sage: CBF(1+i).arcsinh()
            [1.06127506190504 +/- 5.04e-15] + [0.66623943249252 +/- 5.40e-15]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_asinh(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arccosh(self):
        """
        Return the hyperbolic arccosine of this ball.

        EXAMPLES::

            sage: CBF(1+i).arccosh()
            [1.061275061905035 +/- 8.44e-16] + [0.904556894302381 +/- 8.22e-16]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_acosh(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arctanh(self):
        """
        Return the hyperbolic arctangent of this ball.

        EXAMPLES::

            sage: CBF(1+i).arctanh()
            [0.4023594781085251 +/- 8.52e-17] + [1.017221967897851 +/- 4.93e-16]*I
        """
        cdef ComplexBall res = self._new()
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

            sage: CBF(1, 1).gamma()
            [0.49801566811836 +/- 4.98e-15] + [-0.154949828301811 +/- 7.67e-16]*I
            sage: CBF(-1).gamma()
            nan
            sage: CBF(1, 1).gamma(0)
            [0.49801566811836 +/- 4.98e-15] + [-0.154949828301811 +/- 7.67e-16]*I
            sage: CBF(1, 1).gamma(100)
            [-3.6143867454139e-45 +/- 7.26e-59] + [-3.7022961377791e-44 +/- 4.71e-58]*I
            sage: CBF(1, 1).gamma(CLF(i))
            [0.32886684193500 +/- 5.49e-15] + [-0.18974945045621 +/- 1.49e-15]*I
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

    def log_gamma(self):
        r"""
        Return the image of this ball by the logarithmic Gamma function.

        The branch cut of the logarithmic gamma function is placed on the
        negative half-axis, which means that
        ``log_gamma(z) + log z = log_gamma(z+1)`` holds for all `z`,
        whereas ``log_gamma(z) != log(gamma(z))`` in general.

        EXAMPLES::

            sage: CBF(1000, 1000).log_gamma()
            [5466.22252162990 +/- 3.05e-12] + [7039.33429191119 +/- 3.81e-12]*I
            sage: CBF(-1/2).log_gamma()
            [1.265512123484645 +/- 8.82e-16] + [-3.141592653589793 +/- 5.68e-16]*I
            sage: CBF(-1).log_gamma()
            nan + [-3.141592653589793 +/- 5.68e-16]*I
        """
        cdef ComplexBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_lgamma(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def rgamma(self):
        """
        Compute the reciprocal gamma function with argument ``self``.

        EXAMPLES::

            sage: CBF(6).rgamma()
            [0.00833333333333333 +/- 4.96e-18]
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
            [0.0946503206224770 +/- 7.34e-17] + [1.076674047468581 +/- 2.63e-16]*I
            sage: CBF(-1).psi()
            nan
            sage: CBF(1,1).psi(10)
            [56514.8269344249 +/- 4.70e-11] + [56215.1218005823 +/- 5.70e-11]*I

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
            [0.5821580597520036 +/- 5.27e-17] + [-0.9268485643308071 +/- 2.81e-17]*I
            sage: CBF(1, 1).zeta(1)
            [0.5821580597520036 +/- 5.27e-17] + [-0.9268485643308071 +/- 2.81e-17]*I
            sage: CBF(1, 1).zeta(1/2)
            [1.497919876084167 +/- 2.91e-16] + [0.2448655353684164 +/- 4.22e-17]*I
            sage: CBF(1, 1).zeta(CBF(1, 1))
            [-0.3593983122202835 +/- 3.01e-17] + [-2.875283329756940 +/- 4.52e-16]*I
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

    def polylog(self, s):
        """
        Return the polylogarithm `\operatorname{Li}_s(\mathrm{self})`.

        EXAMPLES::

            sage: CBF(2).polylog(1)
            [+/- 4.65e-15] + [-3.14159265358979 +/- 8.15e-15]*I
            sage: CBF(1, 1).polylog(CBF(1, 1))
            [0.3708160030469 +/- 2.38e-14] + [2.7238016577979 +/- 4.22e-14]*I

        TESTS::

            sage: CBF(2).polylog(1r)
            [+/- 4.65e-15] + [-3.14159265358979 +/- 8.15e-15]*I
        """
        cdef ComplexBall s_as_ball
        cdef sage.rings.integer.Integer s_as_Integer
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
            [4.54078781e+254873 +/- 5.43e+254864] + [8.65835455e+254873 +/- 7.27e+254864]*I
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
            [1.14379254649702e+202 +/- 4.09e+187]
            sage: CBF(0,1000).log_barnes_g()
            [-2702305.04929258 +/- 2.60e-9] + [-790386.325561423 +/- 9.72e-10]*I
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
            [0.5990701173678 +/- 1.14e-14] + [-0.5990701173678 +/- 1.22e-14]*I
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

        - ``a`` -- upper parameters, list of complex numbers that coerce into
          this ball's parent;

        - ``b`` -- lower parameters, list of complex numbers that coerce into
          this ball's parent.

        - ``regularized`` -- if True, the regularized generalized hypergeometric
          function is computed.

        OUTPUT:

        The generalized hypergeometric function defined by

        .. math::

            {}_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;z)
            = \sum_{k=0}^\infty \frac{(a_1)_k\dots(a_p)_k}{(b_1)_k\dots(b_q)_k} \frac {z^k} {k!}

        extended using analytic continuation or regularization when the sum
        does not converge.

        The regularized generalized hypergeometric function

        .. math::

            {}_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;z)
            = \sum_{k=0}^\infty \frac{(a_1)_k\dots(a_p)_k}{\Gamma(b_1+k)\dots\Gamma(b_q+k)} \frac {z^k} {k!}

        is well-defined even when the lower parameters are nonpositive
        integers. Currently, this is only supported for some `p` and `q`.

        EXAMPLES::

            sage: CBF(1, pi/2).hypergeometric([], [])
            [+/- 7.72e-16] + [2.71828182845904 +/- 6.45e-15]*I

            sage: CBF(1, pi).hypergeometric([1/4], [1/4])
            [-2.7182818284590 +/- 8.63e-14] + [+/- 3.69e-14]*I

            sage: CBF(1000, 1000).hypergeometric([10], [AA(sqrt(2))])
            [9.79300951360e+454 +/- 4.83e+442] + [5.522579106816e+455 +/- 3.39e+442]*I
            sage: CBF(1000, 1000).hypergeometric([100], [AA(sqrt(2))])
            [+/- 8.88e+596] + [+/- 8.88e+596]*I

            sage: CBF(0, 1).hypergeometric([], [1/2, 1/3, 1/4])
            [-3.7991962344383 +/- 4.98e-14] + [23.8780971778049 +/- 5.40e-14]*I

            sage: CBF(0).hypergeometric([1], [])
            1.000000000000000
            sage: CBF(1, 1).hypergeometric([1], [])
            1.000000000000000*I

            sage: CBF(2+3*I).hypergeometric([1/4,1/3],[1/2])
            [0.7871684267473 +/- 3.57e-14] + [0.2749254173721 +/- 6.45e-14]*I
            sage: CBF(2+3*I).hypergeometric([1/4,1/3],[1/2],regularized=True)
            [0.4441122268685 +/- 1.83e-14] + [0.1551100567338 +/- 4.07e-14]*I

            sage: CBF(5).hypergeometric([2,3], [-5])
            nan + nan*I
            sage: CBF(5).hypergeometric([2,3], [-5], regularized=True)
            [5106.925964355 +/- 5.41e-10]

            sage: CBF(2016).hypergeometric([], [2/3])
            [2.0256426923278e+38 +/- 8.90e+24]
            sage: CBF(-2016).hypergeometric([], [2/3], regularized=True)
            [-0.0005428550847 +/- 5.00e-14]

            sage: CBF(-7).hypergeometric([4], [])
            0.0002441406250000000

            sage: CBF(0, 3).hypergeometric([CBF(1,1)], [-4], regularized=True)
            [239.5140007528 +/- 5.01e-11] + [105.1751573490 +/- 2.36e-11]*I

        TESTS::

            sage: CBF(0, 1).hypergeometric([QQbar(sqrt(2)), RLF(pi)], [1r, 1/2])
            [-8.7029449215408 +/- 6.89e-14] + [-0.8499070546106 +/- 4.98e-14]*I

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
        for i in xrange(p):
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

            sage: CBF(1000, 1000).hypergeometric_U(RLF(pi), -100)
            [-7.261605907166e-11 +/- 4.89e-24] + [-7.928136216391e-11 +/- 5.36e-24]*I
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
            [1.31615128169795 +/- 8.80e-15] + [0.19045346923783 +/- 9.19e-15]*I
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

            sage: CBF(20).erfc()
            [5.3958656116079e-176 +/- 1.08e-190]
            sage: CBF(100, 100).erfc()
            [0.00065234366376858 +/- 6.52e-18] + [-0.00393572636292141 +/- 5.16e-18]*I
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

            sage: CBF(10*pi).airy()
            ([1.2408955946101e-52 +/- 5.50e-66],
             [-6.965048886977e-52 +/- 5.23e-65],
             [2.288295683344e+50 +/- 5.10e+37],
             [1.2807602335816e+51 +/- 4.97e+37])
            sage: ai, aip, bi, bip = CBF(1,2).airy()
            sage: (ai * bip - bi * aip) * CBF(pi)
            [1.0000000000000 +/- 1.25e-15] + [+/- 3.27e-16]*I

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
            [-0.2193862549814276 +/- 7.47e-17] + [-0.1753859114081094 +/- 4.06e-17]*I
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
            [0.1704449781789148 +/- 3.12e-17] + [0.387622439413295 +/- 1.06e-16]*I
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
            [0.0488220324530612 +/- 1.30e-17] + [0.1332740579917484 +/- 6.25e-17]*I
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
            [-0.857239258605362 +/- 3.47e-16] + [0.4955063363095674 +/- 9.22e-17]*I
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
            [0.614160334922903 +/- 8.48e-16] + [0.365028028827088 +/- 6.62e-16]*I
            sage: CBF(100, -100).bessel_J(1/3)
            [1.108431870251e+41 +/- 5.53e+28] + [-8.952577603125e+41 +/- 2.91e+28]*I
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
            [+/- 7.95e-16] + [+/- 6.84e-16]*I
            sage: Y - CBF(1, 1).bessel_Y(1)
            [+/- 2.31e-14] + [+/- 2.26e-14]*I

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
            [-0.6576945355913 +/- 5.62e-14] + [0.6298010039929 +/- 2.77e-14]*I
            sage: CBF(100, -100).bessel_Y(1/3)
            [-8.952577603125e+41 +/- 4.65e+28] + [-1.108431870251e+41 +/- 6.29e+28]*I
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
            [0.365028028827088 +/- 6.62e-16] + [0.614160334922903 +/- 8.48e-16]*I
            sage: CBF(100, -100).bessel_I(1/3)
            [5.4362189595644e+41 +/- 6.48e+27] + [7.1989436985321e+41 +/- 2.69e+27]*I
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
            [0.08019772694652 +/- 3.19e-15] + [-0.35727745928533 +/- 1.08e-15]*I
            sage: CBF(1, 1).bessel_K(1)
            [0.02456830552374 +/- 6.22e-15] + [-0.45971947380119 +/- 6.74e-15]*I
            sage: CBF(100, 100).bessel_K(QQbar(i))
            [3.8693896656383e-45 +/- 2.38e-59] + [5.5071004234177e-46 +/- 5.86e-60]*I
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
            [0.00028162445198 +/- 2.78e-15] + [-0.17932453503936 +/- 2.56e-15]*I
            sage: CBF(1+i).exp_integral_e(QQbar(i))
            [-0.10396361883964 +/- 4.92e-15] + [-0.16268401277783 +/- 4.78e-15]*I
        """
        cdef ComplexBall res = self._new()
        cdef ComplexBall my_s = self._parent.coerce(s)
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_expint(res.value, my_s.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def ei(self):
        """
        Return the exponential integral with argument ``self``.

        EXAMPLES::

            sage: CBF(1, 1).ei()
            [1.76462598556385 +/- 6.65e-15] + [2.38776985151052 +/- 4.34e-15]*I
            sage: CBF(0).ei()
            nan
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_ei(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def si(self):
        """
        Return the sine integral with argument ``self``.

        EXAMPLES::

            sage: CBF(1, 1).si()
            [1.10422265823558 +/- 2.16e-15] + [0.88245380500792 +/- 3.15e-15]*I
            sage: CBF(0).si()
            0
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_si(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def ci(self):
        """
        Return the cosine integral with argument ``self``.

        EXAMPLES::

            sage: CBF(1, 1).ci()
            [0.882172180555936 +/- 4.85e-16] + [0.287249133519956 +/- 3.47e-16]*I
            sage: CBF(0).ci()
            nan + nan*I
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_ci(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def shi(self):
        """
        Return the hyperbolic sine integral with argument ``self``.

        EXAMPLES::

            sage: CBF(1, 1).shi()
            [0.88245380500792 +/- 3.15e-15] + [1.10422265823558 +/- 2.16e-15]*I
            sage: CBF(0).shi()
            0
        """

        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_shi(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def chi(self):
        """
        Return the hyperbolic cosine integral with argument ``self``.

        EXAMPLES::

            sage: CBF(1, 1).chi()
            [0.882172180555936 +/- 4.85e-16] + [1.28354719327494 +/- 1.05e-15]*I
            sage: CBF(0).chi()
            nan + nan*I
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_chi(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def li(self, bint offset=False):
        """
        Return the logarithmic integral with argument ``self``.

        If ``offset`` is True, return the offset logarithmic integral.

        EXAMPLES::

            sage: CBF(1, 1).li()
            [0.61391166922119 +/- 7.03e-15] + [2.05958421419258 +/- 8.25e-15]*I
            sage: CBF(0).li()
            0
            sage: CBF(0).li(offset=True)
            [-1.045163780117493 +/- 5.54e-16]
            sage: li(0).n()
            0.000000000000000
            sage: Li(0).n()
            -1.04516378011749
        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_hypgeom_li(result.value, self.value, offset, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def jacobi_theta(self, tau):
        r"""
        Return the four Jacobi theta functions evaluated at the argument
        ``self`` (representing `z`) and the parameter ``tau`` which should lie
        in the upper half plane.

        The following definitions are used:

        .. math ::

            \theta_1(z,\tau) = 2 q_{1/4} \sum_{n=0}^{\infty} (-1)^n q^{n(n+1)} \sin((2n+1) \pi z)

            \theta_2(z,\tau) = 2 q_{1/4} \sum_{n=0}^{\infty} q^{n(n+1)} \cos((2n+1) \pi z)

            \theta_3(z,\tau) = 1 + 2 \sum_{n=1}^{\infty} q^{n^2} \cos(2n \pi z)

            \theta_4(z,\tau) = 1 + 2 \sum_{n=1}^{\infty} (-1)^n q^{n^2} \cos(2n \pi z)

        where `q = \exp(\pi i \tau)` and `q_{1/4} = \exp(\pi i \tau / 4)`.
        Note that `z` is multiplied by `\pi`; some authors omit this factor.

        EXAMPLES::

            sage: CBF(3,-1/2).jacobi_theta(CBF(1/4,2))
            ([-0.186580562274757 +/- 5.52e-16] + [0.93841744788594 +/- 2.48e-15]*I,
             [-1.02315311037951 +/- 4.10e-15] + [-0.203600094532010 +/- 7.33e-16]*I,
             [1.030613911309632 +/- 4.25e-16] + [0.030613917822067 +/- 1.89e-16]*I,
             [0.969386075665498 +/- 4.65e-16] + [-0.030613917822067 +/- 1.89e-16]*I)

            sage: CBF(3,-1/2).jacobi_theta(CBF(1/4,-2))
            ([+/- inf] + [+/- inf]*I,
             [+/- inf] + [+/- inf]*I,
             [+/- inf] + [+/- inf]*I,
             [+/- inf] + [+/- inf]*I)

            sage: CBF(0).jacobi_theta(CBF(0,1))
            (0,
             [0.91357913815612 +/- 3.96e-15],
             [1.086434811213308 +/- 8.16e-16],
             [0.913579138156117 +/- 8.89e-16])

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
            [1728.0000000000 +/- 5.33e-11]

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
            [0.768225422326057 +/- 9.18e-16]
            sage: CBF(12,1).modular_eta()
            [-0.768225422326057 +/- 9.18e-16]

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

            sage: tau = CBF(sqrt(2),pi)
            sage: tau.modular_lambda()
            [-0.00022005123884157 +/- 6.64e-18] + [-0.0007978734645994 +/- 5.23e-17]*I
            sage: (tau + 2).modular_lambda()
            [-0.00022005123884157 +/- 6.64e-18] + [-0.0007978734645994 +/- 5.23e-17]*I
            sage: (tau / (1 - 2*tau)).modular_lambda()
            [-0.00022005123884 +/- 2.75e-15] + [-0.00079787346460 +/- 3.34e-15]*I

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
            [0.0017853698506421 +/- 6.15e-17]
            sage: a, b, c, d = 2, 5, 1, 3
            sage: tau = CBF(1,3)
            sage: ((a*tau+b)/(c*tau+d)).modular_delta()
            [0.20921376655 +/- 6.94e-12] + [1.5761192552 +/- 3.49e-11]*I
            sage: (c*tau+d)^12 * tau.modular_delta()
            [0.20921376654986 +/- 4.89e-15] + [1.5761192552253 +/- 4.45e-14]*I

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
            [[2.1646498507193 +/- 6.30e-14],
             [2.0346794456073 +/- 2.44e-14],
             [2.0081609898081 +/- 3.67e-14],
             [2.0019857082706 +/- 4.60e-14]]
            sage: ((a*tau+b)/(c*tau+d)).eisenstein(3)[2]
            [331011.200433 +/- 1.36e-7] + [-711178.165575 +/- 5.18e-7]*I
            sage: (c*tau+d)^8 * tau.eisenstein(3)[2]
            [331011.20043304 +/- 7.62e-9] + [-711178.1655746 +/- 1.34e-8]*I

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
            sage: z = CBF(sqrt(2), sqrt(3))
            sage: z.elliptic_p(tau)
            [-3.28920996772709 +/- 7.68e-15] + [-0.000367376730293 +/- 3.58e-16]*I
            sage: (z + tau).elliptic_p(tau)
            [-3.28920996772709 +/- 8.87e-15] + [-0.00036737673029 +/- 4.37e-15]*I
            sage: (z + 1).elliptic_p(tau)
            [-3.28920996772709 +/- 7.68e-15] + [-0.000367376730293 +/- 3.58e-16]*I

            sage: z.elliptic_p(tau, 3)
            [[-3.28920996772709 +/- 7.66e-15] + [-0.000367376730293 +/- 3.45e-16]*I,
             [0.00247305579431 +/- 2.06e-15] + [0.00385955404027 +/- 4.07e-15]*I,
             [-0.0129908756171 +/- 2.12e-14] + [0.0072502752191 +/- 6.05e-14]*I]
            sage: (z + 3 + 4*tau).elliptic_p(tau, 3)
            [[-3.2892099677271 +/- 3.40e-14] + [-0.0003673767303 +/- 2.49e-14]*I,
             [0.00247305579 +/- 6.04e-12] + [0.00385955404 +/- 1.86e-12]*I,
             [-0.012990876 +/- 4.77e-10] + [0.007250275 +/- 3.20e-10]*I]

        """
        cdef ComplexBall my_tau = self._parent.coerce(tau)
        cdef ComplexBall result
        cdef long nn
        cdef acb_ptr vec_r
        if n is None:
            result = self._new()
            if _do_sig(prec(self)): sig_on()
            acb_modular_elliptic_p(result.value, self.value,
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

    def elliptic_k(self):
        """
        Return the complete elliptic integral of the first kind evaluated
        at *m* given by ``self``.

        EXAMPLES::

            sage: CBF(2,3).elliptic_k()
            [1.0429132919285 +/- 5.77e-14] + [0.6296824723086 +/- 7.16e-14]*I

        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_modular_elliptic_k(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def elliptic_e(self):
        """
        Return the complete elliptic integral of the second kind evaluated
        at *m* given by ``self``.

        EXAMPLES::

            sage: CBF(2,3).elliptic_e()
            [1.472797144959 +/- 5.13e-13] + [-1.231604783936 +/- 1.61e-13]*I

        """
        cdef ComplexBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        acb_modular_elliptic_e(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def chebyshev_T(self, n):
        """
        Return the Chebyshev function of the first kind of order ``n``
        evaluated at ``self``.

        EXAMPLES::

            sage: CBF(1/3).chebyshev_T(20)
            [0.8710045668809 +/- 6.15e-14]
            sage: CBF(1/3).chebyshev_T(CBF(5,1))
            [1.8429685451876 +/- 3.57e-14] + [0.20053614301799 +/- 7.05e-15]*I

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
            [0.6973126541184 +/- 2.83e-14]
            sage: CBF(1/3).chebyshev_U(CBF(5,1))
            [1.7588496489342 +/- 5.99e-14] + [0.7497317165104 +/- 4.35e-14]*I

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
            [-920983000.460 +/- 7.60e-4] + [6069919969.93 +/- 2.03e-3]*I

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
            [-263813415.6250000 +/- 9.57e-8]

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
            [-45.6666666666666 +/- 9.28e-14]
            sage: CBF(10).laguerre_L(3, 2)
            [-6.666666666667 +/- 4.15e-13]
            sage: CBF(5,7).laguerre_L(CBF(2,3), CBF(1,-2))
            [5515.31503027 +/- 7.93e-9] + [-12386.94284527 +/- 8.00e-9]*I

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
            [8.05746709617e+37 +/- 1.32e+25]

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
            0.08984375000000000
            sage: CBF(1,2).legendre_P(CBF(2,3), CBF(0,1))
            [0.10996180744 +/- 4.71e-12] + [0.14312767804 +/- 1.62e-12]*I
            sage: CBF(-10).legendre_P(5, 325/100)
            [-22104403.487377 +/- 6.81e-7] + [53364750.687392 +/- 7.25e-7]*I
            sage: CBF(-10).legendre_P(5, 325/100, type=3)
            [-57761589.914581 +/- 6.99e-7] + [+/- 5.14e-7]*I

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
            [0.55508089057168 +/- 5.32e-15]
            sage: CBF(1,2).legendre_Q(CBF(2,3), CBF(0,1))
            [0.167678710 +/- 5.89e-10] + [-0.161558598 +/- 8.76e-10]*I
            sage: CBF(-10).legendre_Q(5, 325/100)
            [-83825154.36008 +/- 4.95e-6] + [-34721515.80396 +/- 5.40e-6]*I
            sage: CBF(-10).legendre_Q(5, 325/100, type=3)
            [-4.797306921692e-6 +/- 6.92e-19] + [-4.797306921692e-6 +/- 6.68e-19]*I

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
            [0.80370071745224 +/- 4.02e-15] + [-0.07282031864711 +/- 4.69e-15]*I
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
