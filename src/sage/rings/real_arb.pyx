# -*- coding: utf-8
r"""
Arbitrary precision real balls using Arb

This is a binding to the `Arb library <http://fredrikj.net/arb/>`_ for ball
arithmetic. It may be useful to refer to its documentation for more details.

Parts of the documentation for this module are copied or adapted from
Arb's own documentation, licenced under the GNU General Public License
version 2, or later.

.. SEEALSO::

    - :mod:`Complex balls using Arb <sage.rings.complex_arb>`
    - :mod:`Real intervals using MPFI <sage.rings.real_mpfi>`

Data Structure
==============

Ball arithmetic, also known as mid-rad interval arithmetic, is an extension of
floating-point arithmetic in which an error bound is attached to each variable.
This allows doing rigorous computations over the real numbers, while avoiding
the overhead of traditional (inf-sup) interval arithmetic at high precision,
and eliminating much of the need for time-consuming and bug-prone manual error
analysis associated with standard floating-point arithmetic.

Sage :class:`RealBall` objects wrap Arb objects of type ``arb_t``. A real
ball represents a ball over the real numbers, that is, an interval `[m-r,m+r]`
where the midpoint `m` and the radius `r` are (extended) real numbers::

    sage: RBF(pi)
    [3.141592653589793 +/- 5.61e-16]
    sage: RBF(pi).mid(), RBF(pi).rad()
    (3.14159265358979, 4.4408921e-16)

The midpoint is represented as an arbitrary-precision floating-point number
with arbitrary-precision exponent. The radius is a floating-point number with
fixed-precision mantissa and arbitrary-precision exponent. ::

    sage: RBF(2)^(2^100)
    [2.285367694229514e+381600854690147056244358827360 +/- 2.98e+381600854690147056244358827344]

:class:`RealBallField` objects (the parents of real balls) model the field of
real numbers represented by balls on which computations are carried out with a
certain precision::

    sage: RBF
    Real ball field with 53 bits precision

It is possible to construct a ball whose parent is the real ball field with
precision `p` but whose midpoint does not fit on `p` bits. However, the results
of operations involving such a ball will (usually) be rounded to its parent's
precision::

    sage: RBF(factorial(50)).mid(), RBF(factorial(50)).rad()
    (3.0414093201713378043612608166064768844377641568961e64, 0.00000000)
    sage: (RBF(factorial(50)) + 0).mid()
    3.04140932017134e64

Comparison
==========

.. WARNING::

    In accordance with the semantics of Arb, identical :class:`RealBall`
    objects are understood to give permission for algebraic simplification.
    This assumption is made to improve performance.  For example, setting ``z =
    x*x`` may set `z` to a ball enclosing the set `\{t^2 : t \in x\}` and not
    the (generally larger) set `\{tu : t \in x, u \in x\}`.

Two elements are equal if and only if they are the same object
or if both are exact and equal::

    sage: a = RBF(1)
    sage: b = RBF(1)
    sage: a is b
    False
    sage: a == b
    True
    sage: a = RBF(1/3)
    sage: b = RBF(1/3)
    sage: a.is_exact()
    False
    sage: b.is_exact()
    False
    sage: a is b
    False
    sage: a == b
    False

A ball is non-zero in the sense of comparison if and only if it does not
contain zero. ::

    sage: a = RBF(RIF(-0.5, 0.5))
    sage: a != 0
    False
    sage: b = RBF(1/3)
    sage: b != 0
    True

However, ``bool(b)`` returns ``False`` for a ball ``b`` only if ``b`` is exactly
zero::

    sage: bool(a)
    True
    sage: bool(b)
    True
    sage: bool(RBF.zero())
    False

A ball ``left`` is less than a ball ``right`` if all elements of
``left`` are less than all elements of ``right``. ::

    sage: a = RBF(RIF(1, 2))
    sage: b = RBF(RIF(3, 4))
    sage: a < b
    True
    sage: a <= b
    True
    sage: a > b
    False
    sage: a >= b
    False
    sage: a = RBF(RIF(1, 3))
    sage: b = RBF(RIF(2, 4))
    sage: a < b
    False
    sage: a <= b
    False
    sage: a > b
    False
    sage: a >= b
    False

Comparisons with Sage symbolic infinities work with some limitations::

    sage: -infinity < RBF(1) < +infinity
    True
    sage: -infinity < RBF(infinity)
    True
    sage: RBF(infinity) < infinity
    False
    sage: RBF(NaN) < infinity
    Traceback (most recent call last):
    ...
    ValueError: infinite but not with +/- phase
    sage: 1/RBF(0) <= infinity
    Traceback (most recent call last):
    ...
    ValueError: infinite but not with +/- phase

Comparisons between elements of real ball fields, however, support special
values and should be preferred::

    sage: RBF(NaN) < RBF(infinity)
    False
    sage: 1/RBF(0) <= RBF(infinity)
    True

TESTS::

    sage: (RBF(pi) * identity_matrix(QQ, 3)).parent()
    Full MatrixSpace of 3 by 3 dense matrices over Real ball field
    with 53 bits precision

    sage: polygen(RBF, x)^3
    x^3

::

    sage: SR.coerce(RBF(0.42))
    [0.4200000000000000 +/- 1.56e-17]
    sage: RBF(0.42) + SR(1)
    [1.420000000000000 +/- 2.94e-16]
    sage: _.parent()
    Symbolic Ring

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

from cpython.float cimport PyFloat_AS_DOUBLE
from cpython.int cimport PyInt_AS_LONG
from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_NE, Py_GT, Py_GE
from libc.stdlib cimport abort

from sage.libs.arb.arb cimport *
from sage.libs.arb.arf cimport (
        arf_init, arf_get_mpfr, arf_set_mpfr, arf_clear, arf_set_mag,
        arf_set, arf_get_d, arf_get_fmpz_2exp, arf_abs_bound_lt_2exp_si,
        ARF_RND_UP, ARF_PREC_EXACT
)
from sage.libs.arb.arf cimport arf_equal, arf_is_nan, arf_is_neg_inf, arf_is_pos_inf, arf_get_mag
from sage.libs.arb.mag cimport mag_init, mag_clear, mag_add, mag_set_d, MAG_BITS, mag_is_inf, mag_is_finite, mag_zero
from sage.libs.flint.flint cimport flint_free
from sage.libs.flint.fmpz cimport (
        fmpz_t, fmpz_init, fmpz_get_mpz, fmpz_set_mpz, fmpz_clear, fmpz_fdiv_ui
)
from sage.libs.flint.fmpq cimport fmpq_t, fmpq_init, fmpq_set_mpq, fmpq_clear
from sage.libs.gmp.mpz cimport mpz_fits_ulong_p, mpz_fits_slong_p, mpz_get_ui, mpz_get_si, mpz_sgn
from sage.libs.mpfi cimport mpfi_get_left, mpfi_get_right, mpfi_interv_fr
from sage.libs.mpfr cimport mpfr_t, mpfr_init2, mpfr_clear, mpfr_sgn, MPFR_PREC_MIN, mpfr_equal_p
from sage.libs.mpfr cimport MPFR_RNDN, MPFR_RNDU, MPFR_RNDD, MPFR_RNDZ
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.real_double cimport RealDoubleElement
from sage.rings.real_mpfr cimport RealField_class, RealField, RealNumber
from sage.rings.ring import Field
from sage.structure.element cimport Element, ModuleElement, RingElement

import operator

import sage.categories.fields
import sage.rings.number_field.number_field as number_field

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_mpfi import RealIntervalField, RealIntervalField_class
from sage.structure.unique_representation import UniqueRepresentation

cdef void mpfi_to_arb(arb_t target, const mpfi_t source, const long precision):
    """
    Convert an MPFI interval to an Arb ball.

    INPUT:

    - ``target`` -- an ``arb_t``.

    - ``source`` -- an ``mpfi_t``.

    - ``precision`` -- an integer `\ge 2`.

    TESTS::

        sage: RBF(RIF(infinity)).endpoints()
        (+infinity, +infinity)
        sage: RBF(RIF(-infinity)).endpoints()
        (-infinity, -infinity)
        sage: RIF(RBF(infinity)).endpoints()
        (+infinity, +infinity)
        sage: RIF(RBF(-infinity)).endpoints()
        (-infinity, -infinity)
    """
    cdef mpfr_t left
    cdef mpfr_t right

    mpfr_init2(left, precision)
    mpfr_init2(right, precision)

    if _do_sig(precision): sig_on()
    mpfi_get_left(left, source)
    mpfi_get_right(right, source)
    arb_set_interval_mpfr(target, left, right, precision)
    # Work around weakness of arb_set_interval_mpfr(tgt, inf, inf)
    if mpfr_equal_p(left, right):
        mag_zero(arb_radref(target))
    if _do_sig(precision): sig_off()

    mpfr_clear(left)
    mpfr_clear(right)

cdef int arb_to_mpfi(mpfi_t target, arb_t source, const long precision) except -1:
    """
    Convert an Arb ball to an MPFI interval.

    INPUT:

    - ``target`` -- an ``mpfi_t``.

    - ``source`` -- an ``arb_t``.

    - ``precision`` -- an integer `\ge 2`.

    EXAMPLES::

        sage: RIF(RBF(2)**(2**100)) # indirect doctest
        Traceback (most recent call last):
        ...
        ArithmeticError: Error converting arb to mpfi. Overflow?
    """
    cdef mpfr_t left
    cdef mpfr_t right

    mpfr_init2(left, precision)
    mpfr_init2(right, precision)

    try:
        sig_on()
        arb_get_interval_mpfr(left, right, source)
        mpfi_interv_fr(target, left, right)
        sig_off()
    except RuntimeError:
        raise ArithmeticError("Error converting arb to mpfi. Overflow?")
    finally:
        mpfr_clear(left)
        mpfr_clear(right)

cdef int real_part_of_quadratic_element_to_arb(arb_t res,
        nfeq.NumberFieldElement_quadratic x, const long prec) except -1:
    r"""
    Convert the real part of a quadratic element to an arb object of type
    ``arb_t``.

    This function does *not* check that the parent has a real or complex
    embedding.

    TESTS::

        sage: NF.<sqrt2> = QuadraticField(2)
        sage: a = (sqrt2 - 1)^1000
        sage: RBF(a)
        [1.676156872756536e-383 +/- 4.39e-399]

        sage: NF.<a> = QuadraticField(-2)
        sage: CBF(1/3 + a).real()
        [0.3333333333333333 +/- 7.04e-17]
    """
    cdef fmpz_t tmpz
    cdef arb_t rootD
    cdef long myprec = prec + 6
    fmpz_init(tmpz)
    arb_init(rootD)
    while True: # a-b√D might cancel
        fmpz_set_mpz(tmpz, x.a)
        arb_set_fmpz(res, tmpz)
        if mpz_sgn(x.D.value) > 0:
            if _do_sig(myprec): sig_on()
            fmpz_set_mpz(tmpz, x.D.value)
            arb_sqrt_fmpz(rootD, tmpz, myprec)
            fmpz_set_mpz(tmpz, x.b)
            if x.standard_embedding:
                arb_addmul_fmpz(res, rootD, tmpz, myprec)
            else:
                arb_submul_fmpz(res, rootD, tmpz, myprec)
            if _do_sig(myprec): sig_off()
            if arb_rel_accuracy_bits(res) < prec - 4:
                myprec *= 2
                continue
        break
    if _do_sig(myprec): sig_on()
    fmpz_set_mpz(tmpz, x.denom)
    arb_div_fmpz(res, res, tmpz, prec)
    arb_clear(rootD)
    fmpz_clear(tmpz)
    if _do_sig(myprec): sig_off()
    return 0

class RealBallField(UniqueRepresentation, Field):
    r"""
    An approximation of the field of real numbers using mid-rad intervals, also
    known as balls.

    INPUT:

    - ``precision`` -- an integer `\ge 2`.

    EXAMPLES::

        sage: RBF = RealBallField() # indirect doctest
        sage: RBF(1)
        1.000000000000000

    ::

        sage: (1/2*RBF(1)) + AA(sqrt(2)) - 1 + polygen(QQ, x)
        x + [0.914213562373095 +/- 4.10e-16]

    TESTS::

        sage: RBF.bracket(RBF(1/2), RBF(1/3))
        [+/- 5.56e-17]
        sage: RBF.cardinality()
        +Infinity
        sage: RBF.cartesian_product(QQ).an_element()**2
        ([1.440000000000000 +/- 4.98e-16], 1/4)
        sage: RBF.coerce_embedding() is None
        True
        sage: loads(dumps(RBF)) is RBF
        True
        sage: RBF['x'].gens_dict_recursive()
        {'x': x}
        sage: RBF.is_finite()
        False
        sage: RBF.is_zero()
        False
        sage: RBF.one()
        1.000000000000000
        sage: RBF.zero()
        0
    """
    Element = RealBall

    @staticmethod
    def __classcall__(cls, long precision=53, category=None):
        r"""
        Normalize the arguments for caching.

        TESTS::

            sage: RealBallField(53) is RealBallField()
            True
        """
        return super(RealBallField, cls).__classcall__(cls, precision, category)

    def __init__(self, precision, category):
        r"""
        Initialize the real ball field.

        INPUT:

        - ``precision`` -- an integer `\ge 2`.

        EXAMPLES::

            sage: RBF = RealBallField()
            sage: RBF(1)
            1.000000000000000
            sage: RealBallField(0)
            Traceback (most recent call last):
            ...
            ValueError: Precision must be at least 2.
            sage: RealBallField(1)
            Traceback (most recent call last):
            ...
            ValueError: Precision must be at least 2.

        TESTS::

            sage: RBF.base()
            Real ball field with 53 bits precision
            sage: RBF.base_ring()
            Real ball field with 53 bits precision

        """
        if precision < 2:
            raise ValueError("Precision must be at least 2.")
        super(RealBallField, self).__init__(
                base_ring=self,
                category=category or sage.categories.fields.Fields().Infinite())
        self._prec = precision
        from sage.rings.qqbar import AA
        from sage.rings.real_lazy import RLF
        self._populate_coercion_lists_([ZZ, QQ, AA, RLF])

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: RealBallField()
            Real ball field with 53 bits precision
            sage: RealBallField(106)
            Real ball field with 106 bits precision
        """
        return "Real ball field with {} bits precision".format(self._prec)

    def _coerce_map_from_(self, other):
        r"""
        Parents that canonically coerce into real ball fields include:

        - some exact or lazy parents representing subsets of the reals, such as
          ``ZZ``, ``QQ``, ``AA``, and ``RLF``;

        - real ball fields with a larger precision.

        TESTS::

            sage: RealBallField().has_coerce_map_from(RealBallField(54))
            True
            sage: RealBallField().has_coerce_map_from(RealBallField(52))
            False
            sage: RealBallField().has_coerce_map_from(RIF)
            False
            sage: RealBallField().has_coerce_map_from(SR)
            False
            sage: RealBallField().has_coerce_map_from(RR)
            False
            sage: RBF.has_coerce_map_from(QuadraticField(2))
            True
            sage: RBF.has_coerce_map_from(QuadraticField(2, embedding=None))
            False
            sage: RBF.has_coerce_map_from(QuadraticField(-2))
            False
        """
        if isinstance(other, RealBallField):
            return (other._prec >= self._prec)
        elif isinstance(other, number_field.NumberField_quadratic):
            emb = other.coerce_embedding()
            if emb is not None:
                return self.has_coerce_map_from(emb.codomain())

    def _element_constructor_(self, mid=None, rad=None):
        """
        Convert ``mid`` to an element of this real ball field, perhaps
        non-canonically.

        In addition to the inputs supported by :meth:`RealBall.__init__`,
        anything that is convertible to a real interval can also be used to
        construct a real ball::

            sage: RBF(RIF(0, 1))                  # indirect doctest
            [+/- 1.01]
            sage: RBF(1)
            1.000000000000000
            sage: RBF(x)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert x to a RealBall

        Various symbolic constants can be converted without going through real
        intervals. (This is faster and yields tighter error bounds.) ::

            sage: RBF(e)
            [2.718281828459045 +/- 5.35e-16]
            sage: RBF(pi)
            [3.141592653589793 +/- 5.61e-16]

        Symbolic expressions are parsed ::

            sage: RBF(4*zeta(3))
            [4.808227612638377 +/- 9.42e-16]
        """
        try:
            return self.element_class(self, mid, rad)
        except TypeError:
            pass
        try:
            return self.element_class(self, mid.pyobject(), rad)
        except (AttributeError, TypeError):
            pass
        try:
            return mid.operator()(*[self(operand) for operand in mid.operands()])
        except (AttributeError, TypeError):
            pass
        try:
            mid = RealIntervalField(self._prec)(mid)
            return self.element_class(self, mid, rad)
        except TypeError:
            pass
        raise TypeError("unable to convert {!r} to a RealBall".format(mid))

    def _repr_option(self, key):
        """
        Declare that real balls print atomically.

        TESTS::

            sage: RBF._repr_option('element_is_atomic')
            True
            sage: RBF['x']([-2,-2,-2/3])
            [-0.666666666666667 +/- 4.82e-16]*x^2 - 2.000000000000000*x
            - 2.000000000000000
            sage: RBF._repr_option('element_is_atomic_typo')
            Traceback (most recent call last):
            ...
            KeyError: 'element_is_atomic_typo'
        """
        if key == 'element_is_atomic':
            return True

        return super(RealBallField, self)._repr_option(key)

    def gens(self):
        r"""
        EXAMPLE::

            sage: RBF.gens()
            (1.000000000000000,)
            sage: RBF.gens_dict()
            {'1.000000000000000': 1.000000000000000}
        """
        return (self.one(),)

    def construction(self):
        """
        Return the construction of a real ball field as a completion of the
        rationals.

        EXAMPLES::

            sage: RBF = RealBallField(42)
            sage: functor, base = RBF.construction()
            sage: functor, base
            (Completion[+Infinity], Rational Field)
            sage: functor(base) is RBF
            True
        """
        from sage.categories.pushout import CompletionFunctor
        functor = CompletionFunctor(sage.rings.infinity.Infinity,
                                    self._prec,
                                    {'type': 'Ball'})
        return functor, QQ

    def complex_field(self):
        """
        Return the complex ball field with the same precision.

        EXAMPLES::

            sage: from sage.rings.complex_arb import ComplexBallField
            sage: RBF.complex_field()
            Complex ball field with 53 bits precision
            sage: RealBallField(3).algebraic_closure()
            Complex ball field with 3 bits precision
        """
        from sage.rings.complex_arb import ComplexBallField
        return ComplexBallField(self._prec)

    algebraic_closure = complex_field

    def precision(self):
        """
        Return the bit precision used for operations on elements of this field.

        EXAMPLES::

            sage: RealBallField().precision()
            53
        """
        return self._prec

    def is_exact(self):
        """
        Real ball fields are not exact.

        EXAMPLES::

            sage: RealBallField().is_exact()
            False
        """
        return False

    def is_finite(self):
        """
        Real ball fields are infinite.

        They already specify it via their category, but we currently need to
        re-implement this method due to the legacy implementation in
        :class:`sage.rings.ring.Ring`.

        EXAMPLES::

            sage: RealBallField().is_finite()
            False
        """
        return False

    def characteristic(self):
        """
        Real ball fields have characteristic zero.

        EXAMPLES::

            sage: RealBallField().characteristic()
            0
        """
        return 0

    def some_elements(self):
        """
        Real ball fields contain exact balls, inexact balls, infinities, and
        more.

        EXAMPLES::

            sage: RBF.some_elements()
            [1.000000000000000,
            [0.3333333333333333 +/- 7.04e-17],
            [-4.733045976388941e+363922934236666733021124 +/- 3.46e+363922934236666733021108],
            [+/- inf],
            [+/- inf],
            nan]
        """
        import sage.symbolic.constants
        return [self(1), self(1)/3,
                -self(2)**(Integer(2)**80),
                self(sage.rings.infinity.Infinity), ~self(0),
                self.element_class(self, sage.symbolic.constants.NotANumber())]

    # Ball functions of non-ball arguments

    def sinpi(self, x):
        """
        Return a ball enclosing sin(πx).

        This works even if ``x`` itself is not a ball, and may be faster or
        more accurate where ``x`` is a rational number.

        EXAMPLES::

            sage: RBF.sinpi(1)
            0
            sage: RBF.sinpi(1/3)
            [0.866025403784439 +/- 5.15e-16]
            sage: RBF.sinpi(1 + 2^(-100))
            [-2.478279624546525e-30 +/- 5.90e-46]

        .. SEEALSO:: :meth:`~sage.rings.real_arb.RealBall.sin`

        TESTS::

            sage: RBF.sinpi(RLF(sqrt(2)))
            [-0.96390253284988 +/- 4.11e-15]
        """
        cdef RealBall res, x_as_ball
        cdef Rational x_as_Rational
        cdef fmpq_t tmpq
        res = self.element_class(self)
        try:
            x_as_Rational = QQ.coerce(x)
            try:
                if _do_sig(self._prec): sig_on()
                fmpq_init(tmpq)
                fmpq_set_mpq(tmpq, x_as_Rational.value)
                arb_sin_pi_fmpq(res.value, tmpq, self._prec)
                if _do_sig(self._prec): sig_off()
            finally:
                fmpq_clear(tmpq)
            return res
        except TypeError:
            pass
        x_as_ball = self.coerce(x)
        if _do_sig(self._prec): sig_on()
        arb_sin_pi(res.value, x_as_ball.value, self._prec)
        if _do_sig(self._prec): sig_off()
        return res

    def cospi(self, x):
        """
        Return a ball enclosing cos(πx).

        This works even if ``x`` itself is not a ball, and may be faster or
        more accurate where ``x`` is a rational number.

        EXAMPLES::

            sage: RBF.cospi(1)
            -1.000000000000000
            sage: RBF.cospi(1/3)
            0.5000000000000000

        .. SEEALSO:: :meth:`~sage.rings.real_arb.RealBall.cos`

        TESTS::

            sage: RBF.cospi(RLF(sqrt(2)))
            [-0.26625534204142 +/- 5.38e-15]
        """
        cdef RealBall res, x_as_ball
        cdef Rational x_as_Rational
        cdef fmpq_t tmpq
        res = self.element_class(self)
        try:
            x_as_Rational = QQ.coerce(x)
            try:
                if _do_sig(self._prec): sig_on()
                fmpq_init(tmpq)
                fmpq_set_mpq(tmpq, x_as_Rational.value)
                arb_cos_pi_fmpq(res.value, tmpq, self._prec)
                if _do_sig(self._prec): sig_off()
            finally:
                fmpq_clear(tmpq)
            return res
        except TypeError:
            pass
        x_as_ball = self.coerce(x)
        if _do_sig(self._prec): sig_on()
        arb_cos_pi(res.value, x_as_ball.value, self._prec)
        if _do_sig(self._prec): sig_off()
        return res

    def gamma(self, x):
        """
        Return a ball enclosing the gamma function of ``x``.

        This works even if ``x`` itself is not a ball, and may be more
        efficient in the case where ``x`` is an integer or a rational number.

        EXAMPLES::

            sage: RBF.gamma(5)
            24.00000000000000
            sage: RBF.gamma(10**20)
            [+/- 5.92e+1956570551809674821757]
            sage: RBF.gamma(1/3)
            [2.678938534707747 +/- 8.99e-16]
            sage: RBF.gamma(-5)
            nan

        .. SEEALSO:: :meth:`~sage.rings.real_arb.RealBall.gamma`

        TESTS::

            sage: RBF.gamma(RLF(pi))
            [2.2880377953400 +/- 4.29e-14]
        """
        cdef RealBall res
        cdef Integer x_as_Integer
        cdef Rational x_as_Rational
        cdef fmpz_t tmpz
        cdef fmpq_t tmpq
        res = self.element_class(self)
        try:
            x_as_Integer = ZZ.coerce(x)
            try:
                if _do_sig(self._prec): sig_on()
                fmpz_init(tmpz)
                fmpz_set_mpz(tmpz, x_as_Integer.value)
                arb_gamma_fmpz(res.value, tmpz, self._prec)
                if _do_sig(self._prec): sig_off()
            finally:
                fmpz_clear(tmpz)
            return res
        except TypeError:
            pass
        try:
            x_as_Rational = QQ.coerce(x)
            try:
                if _do_sig(self._prec): sig_on()
                fmpq_init(tmpq)
                fmpq_set_mpq(tmpq, x_as_Rational.value)
                arb_gamma_fmpq(res.value, tmpq, self._prec)
                if _do_sig(self._prec): sig_off()
            finally:
                fmpq_clear(tmpq)
            return res
        except TypeError:
            pass
        return self.coerce(x).gamma()

    def zeta(self, s):
        """
        Return a ball enclosing the Riemann zeta function of ``s``.

        This works even if ``s`` itself is not a ball, and may be more
        efficient in the case where ``s`` is an integer.

        EXAMPLES::

            sage: RBF.zeta(3) # abs tol 5e-16
            [1.202056903159594 +/- 2.87e-16]
            sage: RBF.zeta(1)
            nan
            sage: RBF.zeta(1/2)
            [-1.460354508809587 +/- 1.94e-16]

        .. SEEALSO:: :meth:`~sage.rings.real_arb.RealBall.zeta`
        """
        cdef RealBall res
        cdef Integer s_as_Integer
        try:
            s_as_Integer = ZZ.coerce(s)
            if mpz_fits_ulong_p(s_as_Integer.value):
                res = self.element_class(self)
                if _do_sig(self._prec): sig_on()
                arb_zeta_ui(res.value, mpz_get_ui(s_as_Integer.value), self._prec)
                if _do_sig(self._prec): sig_off()
                return res
        except TypeError:
            pass
        return self.coerce(s).zeta()

    def bernoulli(self, n):
        """
        Return a ball enclosing the ``n``-th Bernoulli number.

        EXAMPLES::

            sage: [RBF.bernoulli(n) for n in range(4)]
            [1.000000000000000, -0.5000000000000000, [0.1666666666666667 +/- 7.04e-17], 0]
            sage: RBF.bernoulli(2**20)
            [-1.823002872104961e+5020717 +/- 7.16e+5020701]
            sage: RBF.bernoulli(2**1000)
            Traceback (most recent call last):
            ...
            ValueError: argument too large

        TESTS::

            sage: RBF.bernoulli(2r)
            [0.1666666666666667 +/- 7.04e-17]
            sage: RBF.bernoulli(2/3)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Rational Field to Integer Ring
            sage: RBF.bernoulli(-1)
            Traceback (most recent call last):
            ...
            ValueError: expected a nonnegative index
        """
        cdef RealBall res
        cdef Integer n_as_Integer = ZZ.coerce(n)
        if mpz_fits_ulong_p(n_as_Integer.value):
            res = self.element_class(self)
            if _do_sig(self._prec): sig_on()
            arb_bernoulli_ui(res.value, mpz_get_ui(n_as_Integer.value), self._prec)
            if _do_sig(self._prec): sig_off()
            return res
        elif n_as_Integer < 0:
            raise ValueError("expected a nonnegative index")
        else:
            # TODO: Fall back to a Sage implementation in this case?
            raise ValueError("argument too large")

    def fibonacci(self, n):
        """
        Return a ball enclosing the ``n``-th Fibonacci number.

        EXAMPLES::

            sage: [RBF.fibonacci(n) for n in xrange(7)]
            [0,
            1.000000000000000,
            1.000000000000000,
            2.000000000000000,
            3.000000000000000,
            5.000000000000000,
            8.000000000000000]
            sage: RBF.fibonacci(-2)
            -1.000000000000000
            sage: RBF.fibonacci(10**20)
            [3.78202087472056e+20898764024997873376 +/- 4.01e+20898764024997873361]
        """
        cdef fmpz_t tmpz
        cdef RealBall res = self.element_class(self)
        cdef Integer n_as_Integer = ZZ.coerce(n)
        try:
            if _do_sig(self._prec): sig_on()
            fmpz_init(tmpz)
            fmpz_set_mpz(tmpz, n_as_Integer.value)
            arb_fib_fmpz(res.value, tmpz, self._prec)
            if _do_sig(self._prec): sig_off()
        finally:
            fmpz_clear(tmpz)
        return res

    def bell_number(self, n):
        """
        Return a ball enclosing the ``n``-th Bell number.

        EXAMPLES::

            sage: [RBF.bell_number(n) for n in xrange(7)]
            [1.000000000000000,
             1.000000000000000,
             2.000000000000000,
             5.000000000000000,
             15.00000000000000,
             52.00000000000000,
             203.0000000000000]
            sage: RBF.bell_number(-1)
            Traceback (most recent call last):
            ...
            ValueError: expected a nonnegative index
            sage: RBF.bell_number(10**20)
            [5.38270113176282e+1794956117137290721328 +/- 5.44e+1794956117137290721313]
        """
        cdef fmpz_t tmpz
        cdef RealBall res = self.element_class(self)
        cdef Integer n_as_Integer = ZZ.coerce(n)
        if n_as_Integer < 0:
            raise ValueError("expected a nonnegative index")
        try:
            if _do_sig(self._prec): sig_on()
            fmpz_init(tmpz)
            fmpz_set_mpz(tmpz, n_as_Integer.value)
            arb_bell_fmpz(res.value, tmpz, self._prec)
            if _do_sig(self._prec): sig_off()
        finally:
            fmpz_clear(tmpz)
        return res

    def double_factorial(self, n):
        """
        Return a ball enclosing the ``n``-th double factorial.

        EXAMPLES::

            sage: [RBF.double_factorial(n) for n in range(7)]
            [1.000000000000000,
             1.000000000000000,
             2.000000000000000,
             3.000000000000000,
             8.000000000000000,
             15.00000000000000,
             48.00000000000000]
            sage: RBF.double_factorial(2**20)
            [1.4483729903e+2928836 +/- 8.96e+2928825]
            sage: RBF.double_factorial(2**1000)
            Traceback (most recent call last):
            ...
            ValueError: argument too large
            sage: RBF.double_factorial(-1)
            Traceback (most recent call last):
            ...
            ValueError: expected a nonnegative index

        """
        cdef RealBall res
        cdef Integer n_as_Integer = ZZ.coerce(n)
        if mpz_fits_ulong_p(n_as_Integer.value):
            res = self.element_class(self)
            if _do_sig(self._prec): sig_on()
            arb_doublefac_ui(res.value, mpz_get_ui(n_as_Integer.value), self._prec)
            if _do_sig(self._prec): sig_off()
            return res
        elif n_as_Integer < 0:
            raise ValueError("expected a nonnegative index")
        else:
            # TODO: Fall back to a Sage implementation in this case?
            raise ValueError("argument too large")

    def maximal_accuracy(self):
        r"""
        Return the relative accuracy of exact elements measured in bits.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: RBF.maximal_accuracy()
            9223372036854775807 # 64-bit
            2147483647          # 32-bit

        .. SEEALSO:: :meth:`RealBall.accuracy`
        """
        return ARF_PREC_EXACT

cdef inline bint _do_sig(long prec):
    """
    Whether signal handlers should be installed for calls to arb.

    TESTS::

        sage: _ = RealBallField()(1).psi() # indirect doctest
        sage: _ = RealBallField(1500)(1).psi()
    """
    return (prec > 1000)

cdef inline long prec(RealBall ball):
    return ball._parent._prec

cdef class RealBall(RingElement):
    """
    Hold one ``arb_t`` of the `Arb library
    <http://fredrikj.net/arb/>`_

    EXAMPLES::

        sage: a = RealBallField()(RIF(1))                     # indirect doctest
        sage: b = a.psi()
        sage: b
        [-0.577215664901533 +/- 3.85e-16]
        sage: RIF(b)
        -0.577215664901533?
    """

    def __cinit__(self):
        """
        Allocate memory for the encapsulated value.

        EXAMPLES::

            sage: RealBallField()(RIF(1)) # indirect doctest
            1.000000000000000
        """
        arb_init(self.value)

    def __dealloc__(self):
        """
        Deallocate memory of the encapsulated value.

        EXAMPLES::

            sage: a = RealBallField()(RIF(1)) # indirect doctest
            sage: del a
        """
        arb_clear(self.value)

    def __init__(self, parent, mid=None, rad=None):
        """
        Initialize the :class:`RealBall`.

        INPUT:

        - ``parent`` -- a :class:`RealBallField`.

        - ``mid`` (optional) --  ball midpoint, see examples below. If omitted,
          initialize the ball to zero, ignoring the ``rad`` argument.

        - ``rad`` (optional) -- a :class:`RealNumber` or a Python float, ball
          radius. If the midpoint is not exactly representable in
          floating-point, the radius is adjusted to account for the roundoff
          error.

        EXAMPLES::

            sage: RBF = RealBallField()
            sage: RBF()
            0

        One can create exact real balls using elements of various exact parents,
        or using floating-point numbers::

            sage: RBF(3)
            3.000000000000000
            sage: RBF(3r)
            3.000000000000000
            sage: RBF(1/3)
            [0.3333333333333333 +/- 7.04e-17]
            sage: RBF(3.14)
            [3.140000000000000 +/- 1.25e-16]

        ::

            sage: RBF(3, 0.125)
            [3e+0 +/- 0.126]
            sage: RBF(pi, 0.125r)
            [3e+0 +/- 0.267]

        ::

            sage: NF.<sqrt2> = QuadraticField(2)
            sage: RBF(1/5 + sqrt2/2)
            [0.907106781186547 +/- 5.33e-16]

        Note that integers and floating-point numbers are ''not'' rounded to
        the parent's precision::

            sage: b = RBF(11111111111111111111111111111111111111111111111); b
            [1.111111111111111e+46 +/- 1.12e+30]
            sage: b.mid().exact_rational()
            11111111111111111111111111111111111111111111111

        Similarly, converting a real ball from one real ball field to another
        (with a different precision) only changes the way it is displayed and
        the precision of operations involving it, not the actual representation
        of its center::

            sage: RBF100 = RealBallField(100)
            sage: b100 = RBF100(1/3); b100
            [0.333333333333333333333333333333 +/- 4.65e-31]
            sage: b53 = RBF(b100); b53
            [0.3333333333333333 +/- 3.34e-17]
            sage: RBF100(b53)
            [0.333333333333333333333333333333 +/- 4.65e-31]

        Special values are supported::

            sage: RBF(oo).mid(), RBF(-oo).mid(), RBF(unsigned_infinity).mid()
            (+infinity, -infinity, 0.000000000000000)
            sage: RBF(NaN)
            nan

        Strings can be given as input. Strings must contain decimal
        floating-point literals. A valid string must consist of a midpoint,
        a midpoint and a radius separated by "+/-", or just a
        radius prefixed by "+/-". Optionally, the whole string can be enclosed
        in square brackets. In general, the string representation of a
        real ball as returned by ``str()`` can be parsed back (the result
        will be larger than the original ball if rounding occurs).
        A few examples::

            sage: RBF("1.1")
            [1.100000000000000 +/- 3.56e-16]
            sage: RBF(str(RBF("1.1")))
            [1.100000000000000 +/- 7.12e-16]
            sage: RBF("3.25")
            3.250000000000000
            sage: RBF("-3.1 +/- 1e-10")
            [-3.100000000 +/- 1.01e-10]
            sage: RBF("[+/-1]")
            [+/- 1.01]
            sage: RBF("inf +/- inf")
            [+/- inf]

        .. SEEALSO:: :meth:`RealBallField._element_constructor_`

        TESTS::

            sage: from sage.rings.real_arb import RealBall
            sage: RealBall(RBF, sage.symbolic.constants.Pi()) # abs tol 1e-16
            [3.141592653589793 +/- 5.62e-16]
            sage: RealBall(RBF, sage.symbolic.constants.Log2()) # abs tol 1e-16
            [0.693147180559945 +/- 4.06e-16]
            sage: RealBall(RBF, sage.symbolic.constants.Catalan())
            [0.915965594177219 +/- 1.23e-16]
            sage: RealBall(RBF, sage.symbolic.constants.Khinchin())
            [2.685452001065306 +/- 6.82e-16]
            sage: RealBall(RBF, sage.symbolic.constants.Glaisher())
            [1.282427129100623 +/- 6.02e-16]
            sage: RealBall(RBF, sage.symbolic.constants.e)
            [2.718281828459045 +/- 5.35e-16]
            sage: RealBall(RBF, sage.symbolic.constants.EulerGamma())
            [0.577215664901533 +/- 3.57e-16]
            sage: RBF("1 +/- 0.001")
            [1.00 +/- 1.01e-3]
            sage: RBF("2.3e10000000000000000000000 +/- 0.00005e10000000000000000000000")
            [2.3000e+10000000000000000000000 +/- 5.01e+9999999999999999999995]
            sage: RBF("0.3 +/- 0.2 +/- 0.1")
            Traceback (most recent call last):
            ...
            ValueError: unsupported string format

            sage: NF.<a> = QuadraticField(2)
            sage: RBF.coerce(a)
            [1.414213562373095 +/- 3.03e-16]
            sage: NF.<a> = QuadraticField(2, embedding=-1.4)
            sage: RBF(a)
            [-1.414213562373095 +/- 3.03e-16]
            sage: NF.<a> = QuadraticField(2, embedding=None)
            sage: RBF(a)
            Traceback (most recent call last):
            ...
            ValueError: need an embedding
            sage: RBF.coerce(a)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion...
            sage: QQi.<i> = QuadraticField(-1)
            sage: RBF(QQi(3))
            3.000000000000000
            sage: RBF(i)
            Traceback (most recent call last):
            ...
            ValueError: nonzero imaginary part
            sage: RBF.coerce(QQi(3))
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion...
        """
        import sage.symbolic.constants
        cdef fmpz_t tmpz
        cdef fmpq_t tmpq
        cdef arf_t  tmpr
        cdef mag_t  tmpm
        cdef nfeq.NumberFieldElement_quadratic mid_as_qe

        Element.__init__(self, parent)

        if mid is None:
            return

        elif isinstance(mid, RealBall):
            arb_set(self.value, (<RealBall> mid).value) # no rounding!
        elif isinstance(mid, int):
            arb_set_si(self.value, PyInt_AS_LONG(mid)) # no rounding!
        elif isinstance(mid, Integer):
            if _do_sig(prec(self)): sig_on()
            fmpz_init(tmpz)
            fmpz_set_mpz(tmpz, (<Integer> mid).value)
            arb_set_fmpz(self.value, tmpz) # no rounding!
            fmpz_clear(tmpz)
            if _do_sig(prec(self)): sig_off()
        elif isinstance(mid, Rational):
            if _do_sig(prec(self)): sig_on()
            fmpq_init(tmpq)
            fmpq_set_mpq(tmpq, (<Rational> mid).value)
            arb_set_fmpq(self.value, tmpq, prec(self))
            fmpq_clear(tmpq)
            if _do_sig(prec(self)): sig_off()
        elif isinstance(mid, RealNumber):
            if _do_sig(prec(self)): sig_on()
            arf_init(tmpr)
            arf_set_mpfr(tmpr, (<RealNumber> mid).value)
            arb_set_arf(self.value, tmpr) # no rounding!
            arf_clear(tmpr)
            if _do_sig(prec(self)): sig_off()
        elif isinstance(mid, nfeq.NumberFieldElement_quadratic):
            mid_as_qe = <nfeq.NumberFieldElement_quadratic> mid
            if mpz_sgn(mid_as_qe.b) != 0:
                if mpz_sgn(mid_as_qe.D.value) < 0:
                    raise ValueError("nonzero imaginary part")
                elif mid_as_qe._parent._embedding is None:
                    raise ValueError("need an embedding")
            real_part_of_quadratic_element_to_arb(self.value, mid_as_qe, prec(self))
        elif isinstance(mid, sage.rings.infinity.AnInfinity):
            if isinstance(mid, sage.rings.infinity.PlusInfinity):
                arb_pos_inf(self.value)
            elif isinstance(mid, sage.rings.infinity.MinusInfinity):
                arb_neg_inf(self.value)
            else:
                arb_zero_pm_inf(self.value)
        elif isinstance(mid, sage.symbolic.constants.Constant):
            if _do_sig(prec(self)): sig_on()
            try:
                if isinstance(mid, sage.symbolic.constants.NotANumber):
                    arb_indeterminate(self.value)
                elif isinstance(mid, sage.symbolic.constants.Pi):
                    arb_const_pi(self.value, prec(self))
                elif isinstance(mid, sage.symbolic.constants.Log2):
                    arb_const_log2(self.value, prec(self))
                elif isinstance(mid, sage.symbolic.constants.Catalan):
                    arb_const_catalan(self.value, prec(self))
                elif isinstance(mid, sage.symbolic.constants.Khinchin):
                    arb_const_khinchin(self.value, prec(self))
                elif isinstance(mid, sage.symbolic.constants.Glaisher):
                    arb_const_glaisher(self.value, prec(self))
                elif isinstance(mid, sage.symbolic.constants.EulerGamma):
                    arb_const_euler(self.value, prec(self))
                else:
                    raise TypeError("unsupported constant")
            finally:
                if _do_sig(prec(self)): sig_off()
        elif isinstance(mid, sage.symbolic.constants_c.E):
            if _do_sig(prec(self)): sig_on()
            arb_const_e(self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        elif isinstance(mid, RealIntervalFieldElement):
            mpfi_to_arb(self.value,
                (<RealIntervalFieldElement> mid).value,
                prec(self))
        elif isinstance(mid, str):
            if arb_set_str(self.value, mid, prec(self)) != 0:
                raise ValueError("unsupported string format")
        else:
            raise TypeError("unsupported midpoint type")

        if rad is not None:
            mag_init(tmpm)
            if isinstance(rad, RealNumber):
                arf_init(tmpr)
                arf_set_mpfr(tmpr, (<RealNumber> rad).value)
                arf_get_mag(tmpm, tmpr)
                arf_clear(tmpr)
            elif isinstance(rad, float):
                mag_set_d(tmpm, PyFloat_AS_DOUBLE(rad))
            else:
                raise TypeError("rad should be a RealNumber or a Python float")
            mag_add(arb_radref(self.value), arb_radref(self.value), tmpm)
            mag_clear(tmpm)

    cdef RealBall _new(self):
        """
        Return a new real ball element with the same parent as ``self``.

        TESTS::

            sage: RealBallField()(2)**2 # indirect doctest
            4.000000000000000

        """
        cdef RealBall x
        x = RealBall.__new__(RealBall)
        x._parent = self._parent
        return x

    def __hash__(self):
        """
        TESTS::

            sage: hash(RealBallField(10)(1)) == hash(RealBallField(20)(1))
            True
            sage: hash(RBF(1/3)) == hash(RBF(1/3, rad=.1))
            False
            sage: vals = [0, 1, 3/4, 5/8, 7/8, infinity, 'nan', 2^1000 - 1]
            sage: len({hash(RBF(v)) for v in vals}) == len(vals)
            True
        """
        cdef arf_t mid = arb_midref(self.value)
        cdef fmpz_t mant, expo
        fmpz_init(mant)
        fmpz_init(expo)
        arf_get_fmpz_2exp(mant, expo, mid)
        cdef long h = (
                fmpz_fdiv_ui(mant, 1073741789)
                ^ fmpz_fdiv_ui(expo, 2**30)
                ^ (arf_abs_bound_lt_2exp_si(mid) << 10)
                ^ arb_rel_error_bits(self.value) << 20)
        fmpz_clear(expo)
        fmpz_clear(mant)
        return h

    def _repr_(self):
        """
        Return a string representation of ``self``.

        OUTPUT:

        A string.

        EXAMPLES::

           sage: RealBallField()(RIF(1.9, 2))
           [2e+0 +/- 0.101]
        """
        cdef char* c_result
        cdef bytes py_string

        c_result = arb_get_str(self.value, (prec(self) * 31) // 100, 0)
        try:
            py_string = c_result
        finally:
            flint_free(c_result)

        return py_string

    # Conversions

    cpdef RealIntervalFieldElement _real_mpfi_(self, RealIntervalField_class parent):
        """
        Return a :mod:`real interval <sage.rings.real_mpfi>` containing this ball.

        OUTPUT:

        A :class:`~sage.rings.real_mpfi.RealIntervalFieldElement`.

        EXAMPLES::

            sage: a = RealBallField()(RIF(2))
            sage: RIF(a)                                        # indirect doctest
            2
        """
        cdef RealIntervalFieldElement result
        result = parent(0)
        arb_to_mpfi(result.value, self.value, prec(self))
        return result

    def _integer_(self, _):
        """
        Check that this ball contains a single integer and return that integer.

        EXAMPLES::

            sage: ZZ(RBF(1, rad=0.1r))
            1
            sage: ZZ(RBF(1, rad=1.0r))
            Traceback (most recent call last):
            ...
            ValueError: [+/- 2.01] does not contain a unique integer
            sage: ZZ(RBF(pi))
            Traceback (most recent call last):
            ...
            ValueError: [3.141592653589793 +/- 5.61e-16] does not contain a unique integer

        """
        cdef Integer res
        cdef fmpz_t tmp
        fmpz_init(tmp)
        try:
            if arb_get_unique_fmpz(tmp, self.value):
                res = Integer.__new__(Integer)
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

            sage: QQ(RBF(123456/2^12))
            1929/64
            sage: QQ(RBF(1/3))
            Traceback (most recent call last):
            ...
            ValueError: [0.3333333333333333 +/- 7.04e-17] does not contain a unique rational number
        """
        if arb_is_exact(self.value):
            return self.mid().exact_rational()
        else:
            raise ValueError("{} does not contain a unique rational number".format(self))

    def _mpfr_(self, RealField_class field):
        """
        Convert this real ball to a real number.

        This attempts to do something sensible for all rounding modes, as
        illustrated below.

        EXAMPLES::

            sage: mypi = RBF(pi)
            sage: RR(mypi)
            3.14159265358979
            sage: Reals(rnd='RNDU')(mypi)
            3.14159265358980
            sage: Reals(rnd='RNDD')(mypi)
            3.14159265358979
            sage: Reals(rnd='RNDZ')(mypi)
            3.14159265358979
            sage: Reals(rnd='RNDZ')(-mypi)
            -3.14159265358979
            sage: Reals(rnd='RNDU')(-mypi)
            -3.14159265358979

        ::

            sage: b = RBF(RIF(-1/2, 1))
            sage: RR(b)
            0.250000000000000
            sage: Reals(rnd='RNDU')(b)
            1.00000000093133
            sage: Reals(rnd='RNDD')(b)
            -0.500000000931323
            sage: Reals(rnd='RNDZ')(b)
            0.250000000000000
        """
        cdef RealNumber left, mid, right
        cdef long prec = field.precision()
        cdef int sl, sr
        if (field.rnd == MPFR_RNDN or
                field.rnd == MPFR_RNDZ and arb_contains_zero(self.value)):
            mid = RealNumber(field, None)
            sig_str("unable to convert to MPFR (exponent out of range?)")
            arf_get_mpfr(mid.value, arb_midref(self.value), field.rnd)
            sig_off()
            return mid
        else:
            left = RealNumber(field, None)
            right = RealNumber(field, None)
            sig_str("unable to convert to MPFR (exponent out of range?)")
            arb_get_interval_mpfr(left.value, right.value, self.value)
            sig_off()
            if field.rnd == MPFR_RNDD:
                return left
            elif field.rnd == MPFR_RNDU:
                return right
            elif field.rnd == MPFR_RNDZ:
                sl, sr = mpfr_sgn(left.value), mpfr_sgn(left.value)
                if sr > 0 and sl > 0:
                    return left
                elif sr < 0 and sl < 0:
                    return right
                else:
                    return field(0)
        raise ValueError("unknown rounding mode")

    # Center and radius, absolute value, endpoints

    def mid(self):
        """
        Return the center of this ball.

        EXAMPLES::

            sage: RealBallField(16)(1/3).mid()
            0.3333
            sage: RealBallField(16)(1/3).mid().parent()
            Real Field with 16 bits of precision
            sage: RealBallField(16)(RBF(1/3)).mid().parent()
            Real Field with 53 bits of precision
            sage: RBF('inf').mid()
            +infinity

        ::

            sage: b = RBF(2)^(2^1000)
            sage: b.mid()
            Traceback (most recent call last):
            ...
            RuntimeError: unable to convert to MPFR (exponent out of range?)

        .. SEEALSO:: :meth:`rad`, :meth:`squash`
        """
        cdef long mid_prec = max(arb_bits(self.value), prec(self))
        if mid_prec < MPFR_PREC_MIN:
            mid_prec = MPFR_PREC_MIN
        cdef RealField_class mid_field = RealField(mid_prec)
        return self._mpfr_(mid_field)

    center = mid

    def rad(self):
        """
        Return the radius of this ball.

        EXAMPLES::

            sage: RBF(1/3).rad()
            5.5511151e-17
            sage: RBF(1/3).rad().parent()
            Real Field with 30 bits of precision

        .. SEEALSO:: :meth:`mid`, :meth:`rad_as_ball`, :meth:`diameter`

        TESTS::

            sage: (RBF(1, rad=.1) << (2^64)).rad()
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
        sig_str("unable to convert the radius to MPFR (exponent out of range?)")
        arf_set_mag(tmp, arb_radref(self.value))
        if arf_get_mpfr(rad.value, tmp, MPFR_RNDN):
            abort()
        sig_off()
        arf_clear(tmp)
        return rad

    def diameter(self):
        r"""
        Return the diameter of this ball.

        EXAMPLES::

            sage: RBF(1/3).diameter()
            1.1102230e-16
            sage: RBF(1/3).diameter().parent()
            Real Field with 30 bits of precision
            sage: RBF(RIF(1.02, 1.04)).diameter()
            0.020000000

        .. SEEALSO:: :meth:`rad`, :meth:`rad_as_ball`, :meth:`mid`
        """
        return 2 * self.rad()

    def squash(self):
        """
        Return an exact ball with the same center as this ball.

        EXAMPLES::

            sage: mid = RealBallField(16)(1/3).squash()
            sage: mid
            [0.3333 +/- 2.83e-5]
            sage: mid.is_exact()
            True
            sage: mid.parent()
            Real ball field with 16 bits precision

        .. SEEALSO:: :meth:`mid`, :meth:`rad_as_ball`
        """
        cdef RealBall res = self._new()
        arf_set(arb_midref(res.value), arb_midref(self.value))
        mag_zero(arb_radref(res.value))
        return res

    def rad_as_ball(self):
        """
        Return an exact ball with center equal to the radius of this ball.

        EXAMPLES::

            sage: rad = RBF(1/3).rad_as_ball()
            sage: rad
            [5.55111512e-17 +/- 3.13e-26]
            sage: rad.is_exact()
            True
            sage: rad.parent()
            Real ball field with 30 bits precision

        .. SEEALSO:: :meth:`squash`, :meth:`rad`
        """
        cdef RealBall res = self._parent.element_class(RealBallField(MAG_BITS))
        arf_set_mag(arb_midref(res.value), arb_radref(self.value))
        mag_zero(arb_radref(res.value))
        return res

    def __abs__(self):
        """
        Return the absolute value of this ball.

        EXAMPLES::

            sage: RBF(-1/3).abs() # indirect doctest
            [0.3333333333333333 +/- 7.04e-17]
            sage: abs(RBF(-1))
            1.000000000000000
        """
        cdef RealBall r = self._new()
        arb_abs(r.value, self.value)
        return r

    def below_abs(self, test_zero=False):
        """
        Return a lower bound for the absolute value of this ball.

        INPUT:

        - ``test_zero`` (boolean, default ``False``) -- if ``True``,
          make sure that the returned lower bound is positive, raising
          an error if the ball contains zero.

        OUTPUT:

        A ball with zero radius

        EXAMPLES::

            sage: RealBallField(8)(1/3).below_abs()
            [0.33 +/- 7.82e-5]
            sage: b = RealBallField(8)(1/3).below_abs()
            sage: b
            [0.33 +/- 7.82e-5]
            sage: b.is_exact()
            True
            sage: QQ(b)
            169/512

            sage: RBF(0).below_abs()
            0
            sage: RBF(0).below_abs(test_zero=True)
            Traceback (most recent call last):
            ...
            ValueError: ball contains zero

        .. SEEALSO:: :meth:`above_abs`
        """
        cdef RealBall res = self._new()
        arb_get_abs_lbound_arf(arb_midref(res.value), self.value, prec(self))
        if test_zero and arb_contains_zero(res.value):
            assert arb_contains_zero(self.value)
            raise ValueError("ball contains zero")
        return res

    def above_abs(self):
        """
        Return an upper bound for the absolute value of this ball.

        OUTPUT:

        A ball with zero radius

        EXAMPLES::

            sage: b = RealBallField(8)(1/3).above_abs()
            sage: b
            [0.33 +/- 3.99e-3]
            sage: b.is_exact()
            True
            sage: QQ(b)
            171/512

        .. SEEALSO:: :meth:`below_abs`
        """
        cdef RealBall res = self._new()
        arb_get_abs_ubound_arf(arb_midref(res.value), self.value, prec(self))
        return res

    def upper(self, rnd=None):
        """
        Return the right endpoint of this ball, rounded upwards.

        INPUT:

        - ``rnd`` (string) -- rounding mode for the parent of the result (does
          not affect its value!), see
          :meth:`sage.rings.real_mpfi.RealIntervalFieldElement.upper`

        OUTPUT:

        A real number.

        EXAMPLES::

            sage: RBF(-1/3).upper()
            -0.333333333333333
            sage: RBF(-1/3).upper().parent()
            Real Field with 53 bits of precision and rounding RNDU

        .. SEEALSO::

           :meth:`lower`, :meth:`endpoints`
        """
        # naive and slow
        return self._real_mpfi_(RealIntervalField(prec(self))).upper(rnd)

    def lower(self, rnd=None):
        """
        Return the right endpoint of this ball, rounded downwards.

        INPUT:

        - ``rnd`` (string) -- rounding mode for the parent of the result (does
          not affect its value!), see
          :meth:`sage.rings.real_mpfi.RealIntervalFieldElement.lower`

        OUTPUT:

        A real number.

        EXAMPLES::

            sage: RBF(-1/3).lower()
            -0.333333333333334
            sage: RBF(-1/3).lower().parent()
            Real Field with 53 bits of precision and rounding RNDD

        .. SEEALSO:: :meth:`upper`, :meth:`endpoints`
        """
        # naive and slow
        return self._real_mpfi_(RealIntervalField(prec(self))).lower(rnd)

    def endpoints(self, rnd=None):
        """
        Return the endpoints of this ball, rounded outwards.

        INPUT:

        - ``rnd`` (string) -- rounding mode for the parent of the resulting
          floating-point numbers (does not affect their values!), see
          :meth:`sage.rings.real_mpfi.RealIntervalFieldElement.upper`

        OUTPUT:

        A pair of real numbers.

        EXAMPLES::

            sage: RBF(-1/3).endpoints()
            (-0.333333333333334, -0.333333333333333)

        .. SEEALSO:: :meth:`lower`, :meth:`upper`
        """
        # naive and slow
        return self._real_mpfi_(RealIntervalField(prec(self))).endpoints(rnd)

    def union(self, other):
        r"""
        Return a ball containing the convex hull of ``self`` and ``other``.

        EXAMPLES::

            sage: RBF(0).union(1).endpoints()
            (0.000000000000000, 1.00000000000000)
        """
        cdef RealBall my_other = self._parent.coerce(other)
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_union(res.value, self.value, my_other.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    # Precision and accuracy

    def round(self):
        """
        Return a copy of this ball with center rounded to the precision of the
        parent.

        EXAMPLES:

        It is possible to create balls whose midpoint is more precise that
        their parent's nominal precision (see :mod:`~sage.rings.real_arb` for
        more information)::

            sage: b = RBF(pi.n(100))
            sage: b.mid()
            3.141592653589793238462643383

        The ``round()`` method rounds such a ball to its parent's precision::

            sage: b.round().mid()
            3.14159265358979

        .. SEEALSO:: :meth:`trim`
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_set_round(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def accuracy(self):
        """
        Return the effective relative accuracy of this ball measured in bits.

        The accuracy is defined as the difference between the position of the
        top bit in the midpoint and the top bit in the radius, minus one.
        The result is clamped between plus/minus
        :meth:`~RealBallField.maximal_accuracy`.

        EXAMPLES::

            sage: RBF(pi).accuracy()
            51
            sage: RBF(1).accuracy() == RBF.maximal_accuracy()
            True
            sage: RBF(NaN).accuracy() == -RBF.maximal_accuracy()
            True

        .. SEEALSO:: :meth:`~RealBallField.maximal_accuracy`
        """
        return arb_rel_accuracy_bits(self.value)

    def trim(self):
        """
        Return a trimmed copy of this ball.

        Round ``self`` to a number of bits equal to the :meth:`accuracy` of
        ``self`` (as indicated by its radius), plus a few guard bits. The
        resulting ball is guaranteed to contain ``self``, but is more economical
        if ``self`` has less than full accuracy.

        EXAMPLES::

            sage: b = RBF(0.11111111111111, rad=.001)
            sage: b.mid()
            0.111111111111110
            sage: b.trim().mid()
            0.111111104488373

        .. SEEALSO:: :meth:`round`
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_trim(res.value, self.value)
        if _do_sig(prec(self)): sig_off()
        return res

    def add_error(self, ampl):
        """
        Increase the radius of this ball by (an upper bound on) ``ampl``.

        If ``ampl`` is negative, the radius is unchanged.

        INPUT:

        - ``ampl`` -- A real ball (or an object that can be coerced to a real
          ball).

        OUTPUT:

        A new real ball.

        EXAMPLES::

            sage: err = RBF(10^-16)
            sage: RBF(1).add_error(err)
            [1.000000000000000 +/- 1.01e-16]

        TESTS::

            sage: RBF(1).add_error(-1)
            1.000000000000000
            sage: RBF(0).add_error(RBF(1, rad=2.)).endpoints()
            (-3.00000000745059, 3.00000000745059)
        """
        cdef RealBall res = self._new()
        cdef RealBall my_ampl = self._parent(ampl)
        if my_ampl < 0:
            my_ampl = self._parent.zero()
        arb_set(res.value, self.value)
        arb_add_error(res.value, my_ampl.value)
        return res

    # Comparisons and predicates

    def is_zero(self):
        """
        Return ``True`` iff the midpoint and radius of this ball are both zero.

        EXAMPLES::

            sage: RBF = RealBallField()
            sage: RBF(0).is_zero()
            True
            sage: RBF(RIF(-0.5, 0.5)).is_zero()
            False

        .. SEEALSO:: :meth:`is_nonzero`
        """
        return arb_is_zero(self.value)

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

            sage: RBF = RealBallField()
            sage: RBF(pi).is_nonzero()
            True
            sage: RBF(RIF(-0.5, 0.5)).is_nonzero()
            False

        .. SEEALSO:: :meth:`is_zero`
        """
        return arb_is_nonzero(self.value)

    def __nonzero__(self):
        """
        Return ``True`` iff this ball is not the zero ball, i.e. if it its
        midpoint and radius are not both zero.

        This is the preferred way, for instance, to determine the “degree” of a
        polynomial with ball coefficients.

        .. WARNING::

            A “nonzero” ball in the sense of this method may represent the
            value zero. Use :meth:`is_nonzero` to check that a real number
            represented by a ``RealBall`` object is known to be nonzero.

        EXAMPLES::

            sage: bool(RBF(0)) # indirect doctest
            False
            sage: bool(RBF(1/3))
            True
            sage: bool(RBF(RIF(-0.5, 0.5)))
            True
        """
        return not arb_is_zero(self.value)

    def is_exact(self):
        """
        Return ``True`` iff the radius of this ball is zero.

        EXAMPLES::

            sage: RBF = RealBallField()
            sage: RBF(1).is_exact()
            True
            sage: RBF(RIF(0.1, 0.2)).is_exact()
            False
        """
        return arb_is_exact(self.value)

    cpdef _richcmp_(left, right, int op):
        """
        Compare ``left`` and ``right``.

        For more information, see :mod:`sage.rings.real_arb`.

        EXAMPLES::

            sage: RBF = RealBallField()
            sage: a = RBF(1)
            sage: b = RBF(1)
            sage: a is b
            False
            sage: a == b
            True
            sage: a = RBF(1/3)
            sage: a.is_exact()
            False
            sage: b = RBF(1/3)
            sage: b.is_exact()
            False
            sage: a == b
            False

        TESTS:

        Balls whose intersection consists of one point::

            sage: a = RBF(RIF(1, 2))
            sage: b = RBF(RIF(2, 4))
            sage: a < b
            False
            sage: a > b
            False
            sage: a <= b
            False
            sage: a >= b
            False
            sage: a == b
            False
            sage: a != b
            False

        Balls with non-trivial intersection::

            sage: a = RBF(RIF(1, 4))
            sage: a = RBF(RIF(2, 5))
            sage: a < b
            False
            sage: a <= b
            False
            sage: a > b
            False
            sage: a >= b
            False
            sage: a == b
            False
            sage: a != b
            False

        One ball contained in another::

            sage: a = RBF(RIF(1, 4))
            sage: b = RBF(RIF(2, 3))
            sage: a < b
            False
            sage: a <= b
            False
            sage: a > b
            False
            sage: a >= b
            False
            sage: a == b
            False
            sage: a != b
            False

        Disjoint balls::

            sage: a = RBF(1/3)
            sage: b = RBF(1/2)
            sage: a < b
            True
            sage: a <= b
            True
            sage: a > b
            False
            sage: a >= b
            False
            sage: a == b
            False
            sage: a != b
            True

        Exact elements::

            sage: a = RBF(2)
            sage: b = RBF(2)
            sage: a.is_exact()
            True
            sage: b.is_exact()
            True
            sage: a < b
            False
            sage: a <= b
            True
            sage: a > b
            False
            sage: a >= b
            True
            sage: a == b
            True
            sage: a != b
            False

        Special values::

            sage: inf = RBF(+infinity)
            sage: other_inf = RBF(+infinity, 42.r)
            sage: neg_inf = RBF(-infinity)
            sage: extended_line = 1/RBF(0)
            sage: exact_nan = inf - inf
            sage: exact_nan.mid(), exact_nan.rad()
            (NaN, 0.00000000)
            sage: other_exact_nan = inf - inf

        ::

            sage: exact_nan == exact_nan, exact_nan <= exact_nan, exact_nan >= exact_nan
            (False, False, False)
            sage: exact_nan != exact_nan, exact_nan < exact_nan, exact_nan > exact_nan
            (False, False, False)
            sage: from operator import eq, ne, le, lt, ge, gt
            sage: ops = [eq, ne, le, lt, ge, gt]
            sage: any(op(exact_nan, other_exact_nan) for op in ops)
            False
            sage: any(op(exact_nan, b) for op in ops for b in [RBF(1), extended_line, inf, neg_inf])
            False

        ::

            sage: neg_inf < a < inf and inf > a > neg_inf
            True
            sage: neg_inf <= b <= inf and inf >= b >= neg_inf
            True
            sage: neg_inf <= extended_line <= inf and inf >= extended_line >= neg_inf
            True
            sage: neg_inf < extended_line or extended_line < inf
            False
            sage: inf > extended_line or extended_line > neg_inf
            False

        ::

            sage: all(b <= b == b >= b and not (b < b or b != b or b > b)
            ....:     for b in [inf, neg_inf, other_inf])
            True
            sage: any(b1 == b2 for b1 in [inf, neg_inf, a, extended_line]
            ....:              for b2 in [inf, neg_inf, a, extended_line]
            ....:              if not b1 is b2)
            False
            sage: all(b1 != b2 and not b1 == b2
            ....:     for b1 in [inf, neg_inf, a]
            ....:     for b2 in [inf, neg_inf, a]
            ....:     if not b1 is b2)
            True
            sage: neg_inf <= -other_inf == neg_inf == -other_inf < other_inf == inf <= other_inf
            True
            sage: any(inf < b or b > inf
            ....:     for b in [inf, other_inf,  a, extended_line])
            False
            sage: any(inf <= b or b >= inf for b in [a, extended_line])
            False
        """
        cdef RealBall lt, rt
        cdef arb_t difference

        lt = left
        rt = right

        if arb_is_finite(lt.value) or arb_is_finite(rt.value):
            if lt is rt:
                return op == Py_EQ or op == Py_GE or op == Py_LE
            if op == Py_EQ:
                return arb_is_exact(lt.value) and arb_equal(lt.value, rt.value)
            arb_init(difference)
            arb_sub(difference, lt.value, rt.value, prec(lt))
            if op == Py_NE:
                result = arb_is_nonzero(difference)
            elif op == Py_GT:
                result = arb_is_positive(difference)
            elif op == Py_GE:
                result = arb_is_nonnegative(difference)
            elif op == Py_LT:
                result = arb_is_negative(difference)
            elif op == Py_LE:
                result = arb_is_nonpositive(difference)
            arb_clear(difference)
            return result
        elif arf_is_nan(arb_midref(lt.value)) or arf_is_nan(arb_midref(rt.value)):
            return False
        elif mag_is_inf(arb_radref(lt.value)):
            # left is the whole extended real line
            if op == Py_GE:
                return arf_is_neg_inf(arb_midref(rt.value)) and mag_is_finite(arb_radref(rt.value))
            elif op == Py_LE:
                return arf_is_pos_inf(arb_midref(rt.value)) and mag_is_finite(arb_radref(rt.value))
            else:
                return False
        elif mag_is_inf(arb_radref(rt.value)):
            # right is the whole extended real line
            if op == Py_GE:
                return arf_is_pos_inf(arb_midref(lt.value)) and mag_is_finite(arb_radref(lt.value))
            elif op == Py_LE:
                return arf_is_neg_inf(arb_midref(lt.value)) and mag_is_finite(arb_radref(lt.value))
            else:
                return False
        else:
            # both left and right are special, neither is nan, and neither is
            # [-∞,∞], so they are both points at infinity
            if op == Py_EQ:
                return arf_equal(arb_midref(lt.value), arb_midref(rt.value))
            elif op == Py_NE:
                return not arf_equal(arb_midref(lt.value), arb_midref(rt.value))
            elif op == Py_GT:
                return (arf_is_pos_inf(arb_midref(lt.value))
                        and arf_is_neg_inf(arb_midref(rt.value)))
            elif op == Py_LT:
                return (arf_is_neg_inf(arb_midref(lt.value))
                        and arf_is_pos_inf(arb_midref(rt.value)))
            elif op == Py_GE or op == Py_LE:
                return True
        assert False, "not reached"

    def min(self, *others):
        """
        Return a ball containing the minimum of this ball and the
        remaining arguments.

        EXAMPLES::

            sage: RBF(1, rad=.5).min(0)
            0

            sage: RBF(0, rad=2.).min(RBF(0, rad=1.)).endpoints()
            (-2.00000000651926, 1.00000000465662)

            sage: RBF(infinity).min(3, 1/3)
            [0.3333333333333333 +/- 7.04e-17]

        Note that calls involving NaNs try to return a number when possible.
        This is consistent with IEEE-754-2008 but may be surprising. ::

            sage: RBF('nan').min(0)
            0
            sage: RBF('nan').min(RBF('nan'))
            nan

        .. SEEALSO:: :meth:`max`

        TESTS::

            sage: RBF(0).min()
            0
            sage: RBF(infinity).min().rad()
            0.00000000
        """
        iv = self._real_mpfi_(RealIntervalField(prec(self)))
        my_others = [self._parent.coerce(x) for x in others]
        return self._parent(iv.min(*my_others))

    def max(self, *others):
        """
        Return a ball containing the maximum of this ball and the
        remaining arguments.

        EXAMPLES::

            sage: RBF(-1, rad=.5).max(0)
            0

            sage: RBF(0, rad=2.).max(RBF(0, rad=1.)).endpoints()
            (-1.00000000465662, 2.00000000651926)

            sage: RBF(-infinity).max(-3, 1/3)
            [0.3333333333333333 +/- 7.04e-17]

        Note that calls involving NaNs try to return a number when possible.
        This is consistent with IEEE-754-2008 but may be surprising. ::

            sage: RBF('nan').max(0)
            0
            sage: RBF('nan').max(RBF('nan'))
            nan

        .. SEEALSO:: :meth:`min`

        TESTS::

            sage: RBF(0).max()
            0
        """
        iv = self._real_mpfi_(RealIntervalField(prec(self)))
        my_others = [self._parent.coerce(x) for x in others]
        return self._parent(iv.max(*my_others))

    def is_finite(self):
        """
        Return True iff the midpoint and radius of this ball are both
        finite floating-point numbers, i.e. not infinities or NaN.

        EXAMPLES::

            sage: (RBF(2)^(2^1000)).is_finite()
            True
            sage: RBF(oo).is_finite()
            False
        """
        return arb_is_finite(self.value)

    def identical(self, RealBall other):
        """
        Return True iff ``self`` and ``other`` are equal as balls, i.e.
        have both the same midpoint and radius.

        Note that this is not the same thing as testing whether both ``self``
        and ``other`` certainly represent the same real number, unless either
        ``self`` or ``other`` is exact (and neither contains NaN). To test
        whether both operands might represent the same mathematical quantity,
        use :meth:`overlaps` or :meth:`contains`, depending on the
        circumstance.

        EXAMPLES::

            sage: RBF(1).identical(RBF(3)-RBF(2))
            True
            sage: RBF(1, rad=0.25r).identical(RBF(1, rad=0.25r))
            True
            sage: RBF(1).identical(RBF(1, rad=0.25r))
            False
        """
        return arb_equal(self.value, other.value)

    def overlaps(self, RealBall other):
        """
        Return True iff ``self`` and ``other`` have some point in common.

        If either ``self`` or ``other`` contains NaN, this method always
        returns nonzero (as a NaN could be anything, it could in particular
        contain any number that is included in the other operand).

        EXAMPLES::

            sage: RBF(pi).overlaps(RBF(pi) + 2**(-100))
            True
            sage: RBF(pi).overlaps(RBF(3))
            False
        """
        return arb_overlaps(self.value, other.value)

    def contains_exact(self, other):
        """
        Return ``True`` *iff* the given number (or ball) ``other`` is contained
        in the interval represented by ``self``.

        If ``self`` contains NaN, this function always returns ``True`` (as
        it could represent anything, and in particular could represent all the
        points included in ``other``). If ``other`` contains NaN and ``self``
        does not, it always returns ``False``.

        Use ``other in self`` for a test that works for a wider range of inputs
        but may return false negatives.

        EXAMPLES::

            sage: b = RBF(1)
            sage: b.contains_exact(1)
            True
            sage: b.contains_exact(QQ(1))
            True
            sage: b.contains_exact(1.)
            True
            sage: b.contains_exact(b)
            True

        ::

            sage: RBF(1/3).contains_exact(1/3)
            True
            sage: RBF(sqrt(2)).contains_exact(sqrt(2))
            Traceback (most recent call last):
            ...
            TypeError: unsupported type: <type 'sage.symbolic.expression.Expression'>

        TESTS::

            sage: b.contains_exact(1r)
            True

        """
        cdef fmpz_t tmpz
        cdef fmpq_t tmpq
        if _do_sig(prec(self)): sig_on()
        try:
            if isinstance(other, RealBall):
                res = arb_contains(self.value, (<RealBall> other).value)
            elif isinstance(other, int):
                res = arb_contains_si(self.value, PyInt_AS_LONG(other))
            elif isinstance(other, Integer):
                fmpz_init(tmpz)
                fmpz_set_mpz(tmpz, (<Integer> other).value)
                res = arb_contains_fmpz(self.value, tmpz)
                fmpz_clear(tmpz)
            elif isinstance(other, Rational):
                fmpq_init(tmpq)
                fmpq_set_mpq(tmpq, (<Rational> other).value)
                res = arb_contains_fmpq(self.value, tmpq)
                fmpq_clear(tmpq)
            elif isinstance(other, RealNumber):
                res = arb_contains_mpfr(self.value, (<RealNumber> other).value)
            else:
                raise TypeError("unsupported type: " + str(type(other)))
        finally:
            if _do_sig(prec(self)): sig_off()
        return res

    def __contains__(self, other):
        """
        Return True if ``other`` can be verified to be contained in ``self``.

        The test is done using interval arithmetic with a precision determined
        by the parent of ``self`` and may return false negatives.

        EXAMPLES::

            sage: sqrt(2) in RBF(sqrt(2))
            True

        A false negative::

            sage: sqrt(2) in RBF(RealBallField(100)(sqrt(2)))
            False

        .. SEEALSO:: :meth:`contains_exact`
        """
        return self.contains_exact(self._parent(other))

    def contains_zero(self):
        """
        Return ``True`` iff this ball contains zero.

        EXAMPLES::

            sage: RBF(0).contains_zero()
            True
            sage: RBF(RIF(-1, 1)).contains_zero()
            True
            sage: RBF(1/3).contains_zero()
            False
        """
        return arb_contains_zero(self.value)

    def contains_integer(self):
        """
        Return ``True`` iff this ball contains any integer.

        EXAMPLES::

            sage: RBF(3.1, 0.1).contains_integer()
            True
            sage: RBF(3.1, 0.05).contains_integer()
            False
        """
        return arb_contains_int(self.value)

    def is_negative_infinity(self):
        """
        Return ``True`` if this ball is the point -∞.

        EXAMPLES::

            sage: RBF(-infinity).is_negative_infinity()
            True
        """
        return (arf_is_neg_inf(arb_midref(self.value))
                and mag_is_finite(arb_radref(self.value)))

    def is_positive_infinity(self):
        """
        Return ``True`` if this ball is the point +∞.

        EXAMPLES::

            sage: RBF(infinity).is_positive_infinity()
            True
        """
        return (arf_is_pos_inf(arb_midref(self.value))
                and mag_is_finite(arb_radref(self.value)))

    def is_infinity(self):
        """
        Return ``True`` if this ball contains or may represent a point at
        infinity.

        This is the exact negation of :meth:`is_finite`, used in comparisons
        with Sage symbolic infinities.

        .. WARNING::

            Contrary to the usual convention, a return value of True does
            not imply that all points of the ball satisfy the predicate.
            This is due to the way comparisons with symbolic infinities work in
            sage.

        EXAMPLES::

            sage: RBF(infinity).is_infinity()
            True
            sage: RBF(-infinity).is_infinity()
            True
            sage: RBF(NaN).is_infinity()
            True
            sage: (~RBF(0)).is_infinity()
            True
            sage: RBF(42, rad=1.r).is_infinity()
            False
        """
        return not self.is_finite()

    # Arithmetic

    def __neg__(self):
        """
        Return the opposite of this ball.

        EXAMPLES::

            sage: -RBF(1/3)
            [-0.3333333333333333 +/- 7.04e-17]
        """
        cdef RealBall res = self._new()
        arb_neg(res.value, self.value)
        return res

    def __invert__(self):
        """
        Return the inverse of this ball.

        The result is guaranteed to contain the inverse of any point of the
        input ball.

        EXAMPLES::

            sage: ~RBF(5)
            [0.2000000000000000 +/- 4.45e-17]
            sage: ~RBF(0)
            [+/- inf]
            sage: RBF(RIF(-0.1,0.1))
            [+/- 0.101]

        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_inv(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cpdef _add_(self, other):
        """
        Return the sum of two balls, rounded to the ambient field's precision.

        The resulting ball is guaranteed to contain the sums of any two points
        of the respective input balls.

        EXAMPLES::

            sage: RBF(1) + RBF(1/3)
            [1.333333333333333 +/- 5.37e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_add(res.value, self.value, (<RealBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cpdef _sub_(self, other):
        """
        Return the difference of two balls, rounded to the ambient field's
        precision.

        The resulting ball is guaranteed to contain the differences of any two
        points of the respective input balls.

        EXAMPLES::

            sage: RBF(1) - RBF(1/3)
            [0.666666666666667 +/- 5.37e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_sub(res.value, self.value, (<RealBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cpdef _mul_(self, other):
        """
        Return the product of two balls, rounded to the ambient field's
        precision.

        The resulting ball is guaranteed to contain the products of any two
        points of the respective input balls.

        EXAMPLES::

            sage: RBF(-2) * RBF(1/3)
            [-0.666666666666667 +/- 4.82e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_mul(res.value, self.value, (<RealBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cpdef _div_(self, other):
        """
        Return the quotient of two balls, rounded to the ambient field's
        precision.

        The resulting ball is guaranteed to contain the quotients of any two
        points of the respective input balls.

        EXAMPLES::

            sage: RBF(pi)/RBF(e)
            [1.155727349790922 +/- 8.43e-16]
            sage: RBF(2)/RBF(0)
            [+/- inf]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_div(res.value, self.value, (<RealBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def __pow__(base, expo, _):
        """
        EXAMPLES::

            sage: RBF(e)^17
            [24154952.7535753 +/- 9.30e-8]
            sage: RBF(e)^(-1)
            [0.367879441171442 +/- 4.50e-16]
            sage: RBF(e)^(1/2)
            [1.648721270700128 +/- 4.96e-16]
            sage: RBF(e)^RBF(pi)
            [23.1406926327793 +/- 9.16e-14]

        ::

            sage: RBF(-1)^(1/3)
            nan
            sage: RBF(0)^(-1)
            [+/- inf]
            sage: RBF(-e)**RBF(pi)
            nan

        TESTS::

            sage: RBF(e)**(2r)
            [7.38905609893065 +/- 4.68e-15]
            sage: RBF(e)**(-1r)
            [0.367879441171442 +/- 4.50e-16]
        """
        cdef fmpz_t tmpz
        if not isinstance(base, RealBall):
            return sage.structure.element.bin_op(base, expo, operator.pow)
        cdef RealBall self = base
        cdef RealBall res = self._new()
        if isinstance(expo, int) and expo > 0:
            if _do_sig(prec(self)): sig_on()
            arb_pow_ui(res.value, self.value, PyInt_AS_LONG(expo), prec(self))
            if _do_sig(prec(self)): sig_off()
        elif isinstance(expo, Integer):
            if _do_sig(prec(self)): sig_on()
            fmpz_init(tmpz)
            fmpz_set_mpz(tmpz, (<Integer> expo).value)
            arb_pow_fmpz(res.value, self.value, tmpz, prec(self))
            fmpz_clear(tmpz)
            if _do_sig(prec(self)): sig_off()
        elif isinstance(expo, RealBall):
            if _do_sig(prec(self)): sig_on()
            arb_pow(res.value, self.value, (<RealBall> expo).value, prec(self))
            if _do_sig(prec(self)): sig_off()
        else:
            return sage.structure.element.bin_op(base, expo, operator.pow)
        return res

    def sqrt(self):
        """
        Return the square root of this ball.

        EXAMPLES::

            sage: RBF(2).sqrt()
            [1.414213562373095 +/- 2.99e-16]
            sage: RBF(-1/3).sqrt()
            nan
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_sqrt(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def sqrtpos(self):
        """
        Return the square root of this ball, assuming that it represents a
        nonnegative number.

        Any negative numbers in the input interval are discarded.

        EXAMPLES::

            sage: RBF(2).sqrtpos()
            [1.414213562373095 +/- 2.99e-16]
            sage: RBF(-1/3).sqrtpos()
            0
            sage: RBF(0, rad=2.r).sqrtpos()
            [+/- 1.42]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_sqrtpos(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def rsqrt(self):
        """
        Return the reciprocal square root of ``self``.

        At high precision, this is faster than computing a square root.

        EXAMPLES::

            sage: RBF(2).rsqrt()
            [0.707106781186547 +/- 5.73e-16]
            sage: RBF(0).rsqrt()
            nan
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_rsqrt(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def sqrt1pm1(self):
        """
        Return `\sqrt{1+\mathrm{self}}-1`, computed accurately when ``self`` is
        close to zero.

        EXAMPLES::

            sage: eps = RBF(10^(-20))
            sage: (1 + eps).sqrt() - 1
            [+/- 1.12e-16]
            sage: eps.sqrt1pm1()
            [5.00000000000000e-21 +/- 2.54e-36]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_sqrt1pm1(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    # Floor, ceil, etc.

    def floor(self):
        """
        Return the floor of this ball.

        EXAMPLES::

            sage: RBF(1000+1/3, rad=1.r).floor()
            [1.00e+3 +/- 1.01]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_floor(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def ceil(self):
        """
        Return the ceil of this ball.

        EXAMPLES::

            sage: RBF(1000+1/3, rad=1.r).ceil()
            [1.00e+3 +/- 2.01]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_ceil(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def __lshift__(val, shift):
        r"""
        If ``val`` is a ``RealBall`` and ``shift`` is an integer, return the
        ball obtained by shifting the center and radius of ``val`` to the left
        by ``shift`` bits.

        INPUT:

        - ``shift`` -- integer, may be negative.

        EXAMPLES::

            sage: RBF(1/3) << 2 # indirect doctest
            [1.333333333333333 +/- 4.82e-16]
            sage: RBF(1) << -1
            0.5000000000000000

        TESTS::

            sage: RBF(1) << (2^100)
            [2.285367694229514e+381600854690147056244358827360 +/- 2.98e+381600854690147056244358827344]
            sage: RBF(1) << (-2^100)
            [4.375663498372584e-381600854690147056244358827361 +/- 9.57e-381600854690147056244358827378]

            sage: "a" << RBF(1/3)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for <<: 'str' and 'RealBall'
            sage: RBF(1) << RBF(1/3)
            Traceback (most recent call last):
            ...
            TypeError: shift should be an integer
        """
        cdef fmpz_t tmpz
        # the RealBall might be shift, not val
        if not isinstance(val, RealBall):
            raise TypeError("unsupported operand type(s) for <<: '{}' and '{}'"
                            .format(type(val).__name__, type(shift).__name__))
        cdef RealBall self = val
        cdef RealBall res = self._new()
        if isinstance(shift, int):
            arb_mul_2exp_si(res.value, self.value, PyInt_AS_LONG(shift))
        elif isinstance(shift, Integer):
            sig_on()
            fmpz_init(tmpz)
            fmpz_set_mpz(tmpz, (<Integer> shift).value)
            arb_mul_2exp_fmpz(res.value, self.value, tmpz)
            fmpz_clear(tmpz)
            sig_off()
        else:
            raise TypeError("shift should be an integer")
        return res

    def __rshift__(val, shift):
        r"""
        If ``val`` is a ``RealBall`` and ``shift`` is an integer, return the
        ball obtained by shifting the center and radius of ``val`` to the right
        by ``shift`` bits.

        INPUT:

        - ``shift`` -- integer, may be negative.

        EXAMPLES::

            sage: RBF(4) >> 2
            1.000000000000000
            sage: RBF(1/3) >> -2
            [1.333333333333333 +/- 4.82e-16]

        TESTS::

            sage: "a" >> RBF(1/3)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for >>: 'str' and 'RealBall'
        """
        # the RealBall might be shift, not val
        if isinstance(val, RealBall):
            return val << (-shift)
        else:
            raise TypeError("unsupported operand type(s) for >>: '{}' and '{}'"
                            .format(type(val).__name__, type(shift).__name__))

    # Elementary functions

    def log(self, base=None):
        """
        Return the logarithm of this ball.

        INPUT:

        - ``base`` (optional, positive real ball or number) -- if ``None``,
          return the natural logarithm ``ln(self)``, otherwise, return the
          general logarithm ``ln(self)/ln(base)``

        EXAMPLES::

            sage: RBF(3).log()
            [1.098612288668110 +/- 6.63e-16]
            sage: RBF(3).log(2)
            [1.584962500721156 +/- 7.53e-16]

            sage: RBF(-1/3).log()
            nan
            sage: RBF(3).log(-1)
            nan
            sage: RBF(2).log(0)
            nan
        """
        cdef RealBall cst
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_log(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        if base is not None:
            cst = self._parent.coerce(base).log()
            if _do_sig(prec(self)): sig_on()
            arb_div(res.value, res.value, cst.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        return res

    def log1p(self):
        """
        Return ``log(1 + self)``, computed accurately when ``self`` is close to
        zero.

        EXAMPLES::

            sage: eps = RBF(1e-30)
            sage: (1 + eps).log()
            [+/- 2.23e-16]
            sage: eps.log1p()
            [1.00000000000000e-30 +/- 2.68e-46]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_log1p(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def exp(self):
        """
        Return the exponential of this ball.

        EXAMPLES::

            sage: RBF(1).exp()
            [2.718281828459045 +/- 5.41e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_exp(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def expm1(self):
        """
        Return ``exp(self) - 1``, computed accurately when ``self`` is close to
        zero.

        EXAMPLES::

            sage: eps = RBF(1e-30)
            sage: exp(eps) - 1
            [+/- 3.16e-30]
            sage: eps.expm1()
            [1.000000000000000e-30 +/- 8.34e-47]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_expm1(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def sin(self):
        """
        Return the sine of this ball.

        EXAMPLES::

            sage: RBF(pi).sin() # abs tol 1e-16
            [+/- 5.69e-16]

        .. SEEALSO:: :meth:`~sage.rings.real_arb.RealBallField.sinpi`
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_sin(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def cos(self):
        """
        Return the cosine of this ball.

        EXAMPLES::

            sage: RBF(pi).cos() # abs tol 1e-16
            [-1.00000000000000 +/- 6.69e-16]

        .. SEEALSO:: :meth:`~sage.rings.real_arb.RealBallField.cospi`
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_cos(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def tan(self):
        """
        Return the tangent of this ball.

        EXAMPLES::

            sage: RBF(1).tan()
            [1.557407724654902 +/- 3.26e-16]
            sage: RBF(pi/2).tan()
            [+/- inf]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_tan(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def cot(self):
        """
        Return the cotangent of this ball.

        EXAMPLES::

            sage: RBF(1).cot()
            [0.642092615934331 +/- 4.79e-16]
            sage: RBF(pi).cot()
            [+/- inf]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_cot(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arcsin(self):
        """
        Return the arcsine of this ball.

        EXAMPLES::

            sage: RBF(1).arcsin()
            [1.570796326794897 +/- 6.65e-16]
            sage: RBF(1, rad=.125r).arcsin()
            nan
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_asin(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arccos(self):
        """
        Return the arccosine of this ball.

        EXAMPLES::

            sage: RBF(1).arccos()
            0
            sage: RBF(1, rad=.125r).arccos()
            nan
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_acos(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arctan(self):
        """
        Return the arctangent of this ball.

        EXAMPLES::

            sage: RBF(1).arctan()
            [0.785398163397448 +/- 3.91e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_atan(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def sinh(self):
        """
        Return the hyperbolic sine of this ball.

        EXAMPLES::

            sage: RBF(1).sinh()
            [1.175201193643801 +/- 6.18e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_sinh(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def cosh(self):
        """
        Return the hyperbolic cosine of this ball.

        EXAMPLES::

            sage: RBF(1).cosh()
            [1.543080634815244 +/- 5.28e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_cosh(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def tanh(self):
        """
        Return the hyperbolic tangent of this ball.

        EXAMPLES::

            sage: RBF(1).tanh()
            [0.761594155955765 +/- 2.81e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_tanh(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def coth(self):
        """
        Return the hyperbolic cotangent of this ball.

        EXAMPLES::

            sage: RBF(1).coth()
            [1.313035285499331 +/- 4.97e-16]
            sage: RBF(0).coth()
            [+/- inf]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_coth(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arcsinh(self):
        """
        Return the inverse hyperbolic sine of this ball.

        EXAMPLES::

            sage: RBF(1).arcsinh()
            [0.881373587019543 +/- 1.87e-16]
            sage: RBF(0).arcsinh()
            0
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_asinh(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arccosh(self):
        """
        Return the inverse hyperbolic cosine of this ball.

        EXAMPLES::

            sage: RBF(2).arccosh()
            [1.316957896924817 +/- 6.61e-16]
            sage: RBF(1).arccosh()
            0
            sage: RBF(0).arccosh()
            nan
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_acosh(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arctanh(self):
        """
        Return the inverse hyperbolic tangent of this ball.

        EXAMPLES::

            sage: RBF(0).arctanh()
            0
            sage: RBF(1/2).arctanh()
            [0.549306144334055 +/- 3.32e-16]
            sage: RBF(1).arctanh()
            nan
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_atanh(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    # Special functions

    def gamma(self):
        """
        Return the image of this ball by the Euler Gamma function.

        For integer and rational arguments,
        :meth:`~sage.rings.real_arb.RealBallField.gamma` may be faster.

        EXAMPLES::

            sage: RBF(1/2).gamma()
            [1.772453850905516 +/- 3.41e-16]

        .. SEEALSO::
            :meth:`~sage.rings.real_arb.RealBallField.gamma`
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_gamma(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def log_gamma(self):
        """
        Return the image of this ball by the logarithmic Gamma function.

        The complex branch structure is assumed, so if ``self`` <= 0, the result
        is an indeterminate interval.

        EXAMPLES::

            sage: RBF(1/2).log_gamma()
            [0.572364942924700 +/- 2.67e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_lgamma(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def rgamma(self):
        """
        Return the image of this ball by the function 1/Γ, avoiding division by
        zero at the poles of the gamma function.

        EXAMPLES::

            sage: RBF(-1).rgamma()
            0
            sage: RBF(3).rgamma()
            0.5000000000000000
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_rgamma(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def rising_factorial(self, n):
        """
        Return the ``n``-th rising factorial of this ball.

        The `n`-th rising factorial of `x` is equal to `x (x+1) \cdots (x+n-1)`.

        For real `n`, it is a quotient of gamma functions.

        EXAMPLES::

            sage: RBF(1).rising_factorial(5)
            120.0000000000000
            sage: RBF(1/2).rising_factorial(1/3)
            [0.63684988431797 +/- 5.71e-15]
        """
        cdef RealBall result = self._new()
        cdef RealBall my_n = self._parent.coerce(n)
        if _do_sig(prec(self)): sig_on()
        arb_rising(result.value, self.value, my_n.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    cpdef RealBall psi(self):
        """
        Compute the digamma function with argument self.

        EXAMPLES::

            sage: RBF(1).psi()
            [-0.577215664901533 +/- 3.85e-16]
        """

        cdef RealBall result = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_digamma(result.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    def zeta(self, a=None):
        """
        Return the image of this ball by the Hurwitz zeta function.

        For ``a = 1`` (or ``a = None``), this computes the Riemann zeta function.

        Use :meth:`RealBallField.zeta` to compute the Riemann zeta function of
        a small integer without first converting it to a real ball.

        EXAMPLES::

            sage: RBF(-1).zeta()
            [-0.0833333333333333 +/- 4.26e-17]
            sage: RBF(-1).zeta(1)
            [-0.0833333333333333 +/- 4.26e-17]
            sage: RBF(-1).zeta(2)
            [-1.083333333333333 +/- 4.08e-16]
        """
        cdef RealBall a_ball
        cdef RealBall res = self._new()
        if a is None:
            if _do_sig(prec(self)): sig_on()
            arb_zeta(res.value, self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        else:
            a_ball = self._parent.coerce(a)
            if _do_sig(prec(self)): sig_on()
            arb_hurwitz_zeta(res.value, self.value, a_ball.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        return res

    def polylog(self, s):
        """
        Return the polylogarithm `\operatorname{Li}_s(\mathrm{self})`.

        EXAMPLES::

            sage: polylog(0, -1)
            -1/2
            sage: RBF(-1).polylog(0)
            [-0.50000000000000 +/- 1.29e-15]
            sage: polylog(1, 1/2)
            -log(1/2)
            sage: RBF(1/2).polylog(1)
            [0.6931471805599 +/- 5.02e-14]
            sage: RBF(1/3).polylog(1/2)
            [0.44210883528067 +/- 6.75e-15]
            sage: RBF(1/3).polylog(RLF(pi))
            [0.34728895057225 +/- 5.51e-15]

        TESTS::

            sage: RBF(1/3).polylog(2r)
            [0.36621322997706 +/- 4.62e-15]
        """
        cdef RealBall s_as_ball
        cdef Integer s_as_Integer
        cdef RealBall res = self._new()
        try:
            s_as_Integer = ZZ.coerce(s)
            if mpz_fits_slong_p(s_as_Integer.value):
                if _do_sig(prec(self)): sig_on()
                arb_polylog_si(res.value, mpz_get_si(s_as_Integer.value), self.value, prec(self))
                if _do_sig(prec(self)): sig_off()
                return res
        except TypeError:
            pass
        s_as_ball = self._parent.coerce(s)
        if _do_sig(prec(self)): sig_on()
        arb_polylog(res.value, s_as_ball.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def chebyshev_T(self, n):
        """
        Evaluate the Chebyshev polynomial of the first kind ``T_n`` at this
        ball.

        EXAMPLES::

            sage: RBF(pi).chebyshev_T(0)
            1.000000000000000
            sage: RBF(pi).chebyshev_T(1) # abs tol 1e-16
            [3.141592653589793 +/- 5.62e-16]
            sage: RBF(pi).chebyshev_T(10**20)
            Traceback (most recent call last):
            ...
            ValueError: index too large
            sage: RBF(pi).chebyshev_T(-1)
            Traceback (most recent call last):
            ...
            ValueError: expected a nonnegative index
        """
        cdef RealBall res = self._new()
        cdef Integer n_as_Integer = ZZ.coerce(n)
        if mpz_fits_ulong_p(n_as_Integer.value):
            if _do_sig(prec(self)): sig_on()
            arb_chebyshev_t_ui(res.value, mpz_get_ui(n_as_Integer.value), self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
            return res
        elif n_as_Integer < 0:
            raise ValueError("expected a nonnegative index")
        else:
            raise ValueError("index too large")

    def chebyshev_U(self, n):
        """
        Evaluate the Chebyshev polynomial of the second kind ``U_n`` at this
        ball.

        EXAMPLES::

            sage: RBF(pi).chebyshev_U(0)
            1.000000000000000
            sage: RBF(pi).chebyshev_U(1)
            [6.28318530717959 +/- 4.66e-15]
            sage: RBF(pi).chebyshev_U(10**20)
            Traceback (most recent call last):
            ...
            ValueError: index too large
            sage: RBF(pi).chebyshev_U(-1)
            Traceback (most recent call last):
            ...
            ValueError: expected a nonnegative index
        """
        cdef RealBall res = self._new()
        cdef Integer n_as_Integer = ZZ.coerce(n)
        if mpz_fits_ulong_p(n_as_Integer.value):
            if _do_sig(prec(self)): sig_on()
            arb_chebyshev_u_ui(res.value, mpz_get_ui(n_as_Integer.value), self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
            return res
        elif n_as_Integer < 0:
            raise ValueError("expected a nonnegative index")
        else:
            raise ValueError("index too large")

    def agm(self, other):
        """
        Return the arithmetic-geometric mean of ``self`` and ``other``.

        EXAMPLES::

            sage: RBF(1).agm(1)
            1.000000000000000
            sage: RBF(sqrt(2)).agm(1)^(-1)
            [0.83462684167407 +/- 4.31e-15]
        """
        cdef RealBall other_as_ball
        cdef RealBall res = self._new()
        other_as_ball = self._parent.coerce(other)
        if _do_sig(prec(self)): sig_on()
        arb_agm(res.value, self.value, other_as_ball.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

RBF = RealBallField()
