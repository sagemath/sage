r"""
Arbitrary precision real balls

This is a binding to the `arb module of FLINT <https://flintlib.org/doc/arb.html>`_.
It may be useful to refer to its documentation for more details.

Parts of the documentation for this module are copied or adapted from Arb's
(now FLINT's) own documentation, licenced (at the time) under the GNU General
Public License version 2, or later.

.. SEEALSO::

    - :mod:`Complex balls <sage.rings.complex_arb>`
    - :mod:`Real intervals using MPFI <sage.rings.real_mpfi>`

Data Structure
==============

Ball arithmetic, also known as mid-rad interval arithmetic, is an extension of
floating-point arithmetic in which an error bound is attached to each variable.
This allows doing rigorous computations over the real numbers, while avoiding
the overhead of traditional (inf-sup) interval arithmetic at high precision,
and eliminating much of the need for time-consuming and bug-prone manual error
analysis associated with standard floating-point arithmetic.

Sage :class:`RealBall` objects wrap FLINT objects of type ``arb_t``. A real
ball represents a ball over the real numbers, that is, an interval `[m-r,m+r]`
where the midpoint `m` and the radius `r` are (extended) real numbers::

    sage: RBF(pi)                                                                       # needs sage.symbolic
    [3.141592653589793 +/- ...e-16]
    sage: RBF(pi).mid(), RBF(pi).rad()                                                  # needs sage.symbolic
    (3.14159265358979, ...e-16)

The midpoint is represented as an arbitrary-precision floating-point number
with arbitrary-precision exponent. The radius is a floating-point number with
fixed-precision mantissa and arbitrary-precision exponent. ::

    sage: RBF(2)^(2^100)
    [2.285367694229514e+381600854690147056244358827360 +/- ...e+381600854690147056244358827344]

:class:`RealBallField` objects (the parents of real balls) model the field of
real numbers represented by balls on which computations are carried out with a
certain precision::

    sage: RBF
    Real ball field with 53 bits of precision

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

    In accordance with the semantics of FLINT/Arb, identical :class:`RealBall`
    objects are understood to give permission for algebraic simplification.
    This assumption is made to improve performance.  For example, setting ``z =
    x*x`` may set `z` to a ball enclosing the set `\{t^2 : t \in x\}` and not
    the (generally larger) set `\{tu : t \in x, u \in x\}`.

Two elements are equal if and only if they are exact and equal (in spite of the
above warning, inexact balls are not considered equal to themselves)::

    sage: a = RBF(1)
    sage: b = RBF(1)
    sage: a is b
    False
    sage: a == a
    True
    sage: a == b
    True

::

    sage: a = RBF(1/3)
    sage: b = RBF(1/3)
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

A ball is nonzero in the sense of comparison if and only if it does not
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
    sage: RBF(NaN) < infinity                                                           # needs sage.symbolic
    Traceback (most recent call last):
    ...
    ValueError: infinite but not with +/- phase
    sage: 1/RBF(0) <= infinity
    Traceback (most recent call last):
    ...
    ValueError: infinite but not with +/- phase

Comparisons between elements of real ball fields, however, support special
values and should be preferred::

    sage: RBF(NaN) < RBF(infinity)                                                      # needs sage.symbolic
    False
    sage: RBF(0).add_error(infinity) <= RBF(infinity)
    True

TESTS::

    sage: (RBF(pi) * identity_matrix(QQ, 3)).parent()                                   # needs sage.symbolic
    Full MatrixSpace of 3 by 3 dense matrices over Real ball field
    with 53 bits of precision

    sage: polygen(RBF, 'x')^3
    x^3

::

    sage: SR.coerce(RBF(0.42))                                                          # needs sage.symbolic
    [0.4200000000000000 +/- ...e-17]
    sage: RBF(0.42) + SR(1)                                                             # needs sage.symbolic
    [1.420000000000000 +/- ...e-16]
    sage: _.parent()                                                                    # needs sage.symbolic
    Symbolic Ring

Classes and Methods
===================
"""

# ****************************************************************************
#       Copyright (C) 2014 Clemens Heuberger <clemens.heuberger@aau.at>
#                     2017 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_str, sig_off

from cpython.float cimport PyFloat_AS_DOUBLE
from cpython.long cimport PyLong_AsLong
from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_NE, Py_GT, Py_GE
from libc.stdlib cimport abort

from sage.libs.flint.arb cimport *
from sage.libs.flint.arb_hypgeom cimport *
from sage.libs.flint.arf cimport *
from sage.libs.flint.arf cimport *
from sage.libs.flint.mag cimport *
from sage.libs.flint.flint cimport flint_free
from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpq cimport *
from sage.libs.gmp.mpz cimport *
from sage.libs.mpfi cimport *
from sage.libs.mpfr cimport *
from sage.libs.mpfr cimport MPFR_RNDN, MPFR_RNDU, MPFR_RNDD, MPFR_RNDZ

from sage.structure.element cimport Element, RingElement
from sage.rings.ring cimport Field
import sage.rings.abc
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.real_mpfr cimport RealField_class, RealField, RealNumber
from sage.arith.long cimport is_small_python_int

import operator

import sage.categories.fields

from sage.misc.lazy_string import lazy_string
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_mpfi import RealIntervalField, RealIntervalField_class
from sage.structure.unique_representation import UniqueRepresentation
from sage.cpython.string cimport char_to_str, str_to_bytes

cdef void mpfi_to_arb(arb_t target, const mpfi_t source, const long precision) noexcept:
    r"""
    Convert an MPFI interval to an Arb ball.

    INPUT:

    - ``target`` -- an ``arb_t``

    - ``source`` -- an ``mpfi_t``

    - ``precision`` -- integer `\ge 2`

    TESTS::

        sage: RBF(RIF(infinity)).endpoints()
        (+infinity, +infinity)
        sage: RBF(RIF(-infinity)).endpoints()
        (-infinity, -infinity)
        sage: RBF(RIF(-infinity, infinity)).endpoints()
        (-infinity, +infinity)
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
    if mpfr_inf_p(left) and mpfr_inf_p(right) and mpfr_sgn(left) < 0 < mpfr_sgn(right):
        # Work around a weakness of arb_set_interval_mpfr(tgt, -inf, inf)
        arb_zero_pm_inf(target)
    else:
        arb_set_interval_mpfr(target, left, right, precision)
    if _do_sig(precision): sig_off()

    mpfr_clear(left)
    mpfr_clear(right)

cdef int arb_to_mpfi(mpfi_t target, arb_t source, const long precision) except -1:
    r"""
    Convert an Arb ball to an MPFI interval.

    INPUT:

    - ``target`` -- an ``mpfi_t``

    - ``source`` -- an ``arb_t``

    - ``precision`` -- integer `\ge 2`

    EXAMPLES::

        sage: RIF(RBF(2)**(2**100))
        [5.8756537891115869e1388255822130839282 .. +infinity] # 64-bit
        [2.098... .. +infinity]                               # 32-bit
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


class RealBallField(UniqueRepresentation, sage.rings.abc.RealBallField):
    r"""
    An approximation of the field of real numbers using mid-rad intervals, also
    known as balls.

    INPUT:

    - ``precision`` -- integer `\ge 2`

    EXAMPLES::

        sage: RBF = RealBallField() # indirect doctest
        sage: RBF(1)
        1.000000000000000

    ::

        sage: (1/2*RBF(1)) + AA(sqrt(2)) - 1 + polygen(QQ, 'x')                         # needs sage.symbolic
        x + [0.914213562373095 +/- ...e-16]

    TESTS::

        sage: RBF.bracket(RBF(1/2), RBF(1/3))
        [+/- ...e-17]
        sage: RBF.cardinality()
        +Infinity
        sage: RBF.cartesian_product(QQ).an_element()**2
        ([1.440000000000000 +/- ...e-16], 1/4)
        sage: RBF.coerce_embedding() is None
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

        sage: NF.<sqrt2> = QuadraticField(2, embedding=AA(2).sqrt())                    # needs sage.rings.number_field
        sage: a = (sqrt2 - 1)^1000                                                      # needs sage.rings.number_field
        sage: RBF(a)                                                                    # needs sage.rings.number_field
        [1.676156872756536e-383 +/- ...e-399]

        sage: RealBallField().is_finite()
        False

        sage: loads(dumps(RealBallField(60))) is RealBallField(60)
        True

    .. SEEALSO::

        - :mod:`sage.rings.real_arb`
        - :mod:`sage.rings.real_mpfr`
        - :mod:`sage.rings.real_mpfi` (real intervals represented by their
          endpoints)
        - :mod:`sage.rings.complex_arb`
    """
    Element = RealBall

    @staticmethod
    def __classcall__(cls, long precision=53):
        r"""
        Normalize the arguments for caching.

        TESTS::

            sage: RealBallField(53) is RealBallField() is RBF
            True
        """
        return super().__classcall__(cls, precision)

    def __init__(self, long precision=53):
        r"""
        Initialize the real ball field.

        INPUT:

        - ``precision`` -- integer `\ge 2`

        EXAMPLES::

            sage: RBF = RealBallField()
            sage: RBF(1)
            1.000000000000000
            sage: RealBallField(0)
            Traceback (most recent call last):
            ...
            ValueError: precision must be at least 2
            sage: RealBallField(1)
            Traceback (most recent call last):
            ...
            ValueError: precision must be at least 2

        TESTS::

            sage: RBF.base()
            Real ball field with 53 bits of precision
            sage: RBF.base_ring()
            Real ball field with 53 bits of precision
        """
        if precision < 2:
            raise ValueError("precision must be at least 2")
        Field.__init__(self,
                base_ring=self,
                category=sage.categories.fields.Fields().Infinite())
        self._prec = precision
        from sage.rings.real_lazy import RLF
        self._populate_coercion_lists_(coerce_list=[ZZ, QQ], convert_method_name='_arb_')

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: RealBallField()
            Real ball field with 53 bits of precision
            sage: RealBallField(106)
            Real ball field with 106 bits of precision
        """
        return "Real ball field with {} bits of precision".format(self._prec)

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
            sage: RealBallField().has_coerce_map_from(SR)                               # needs sage.symbolic
            False
            sage: RealBallField().has_coerce_map_from(RR)
            False
            sage: K = QuadraticField(2, embedding=AA(2).sqrt())
            sage: RBF.has_coerce_map_from(K)
            True
            sage: RBF.has_coerce_map_from(QuadraticField(2, embedding=None))
            False
            sage: RBF.has_coerce_map_from(QuadraticField(-2))
            False

        Check that the map goes through the ``_arb_`` method::

            sage: RBF.coerce_map_from(QuadraticField(2, embedding=AA(2).sqrt()))
            Conversion via _arb_ method map:
            ...
            sage: RBF.convert_map_from(QuadraticField(2))
            Conversion via _arb_ method map:
            ...
            sage: RBF.coerce_map_from(AA)
            Conversion via _arb_ method map:
            ...
            sage: RBF.convert_map_from(QQbar)
            Conversion via _arb_ method map:
            ...
        """
        if isinstance(other, RealBallField):
            return other._prec >= self._prec

        from sage.rings.qqbar import AA
        from sage.rings.real_lazy import RLF
        if other in [AA, RLF]:
            return True

        from sage.rings.number_field.number_field_base import NumberField
        if isinstance(other, NumberField):
            emb = other.coerce_embedding()
            return emb is not None and self.has_coerce_map_from(emb.codomain())

    def _element_constructor_(self, mid=None, rad=None):
        """
        Convert ``mid`` to an element of this real ball field, perhaps
        non-canonically.

        In addition to the inputs supported by :meth:`RealBall.__init__`,
        elements that can be coerced to real intervals can also be used to
        construct a real ball::

            sage: RBF(RIF(0, 1))                  # indirect doctest
            [+/- 1.01]
            sage: RBF(1)
            1.000000000000000
            sage: RBF(x)                                                                # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: unable to convert x to a RealBall

        Various symbolic constants are supported::

            sage: RBF(e)                                                                # needs sage.symbolic
            [2.718281828459045 +/- ...e-16]
            sage: RBF(pi)                                                               # needs sage.symbolic
            [3.141592653589793 +/- ...e-16]

        Symbolic expressions are parsed ::

            sage: RBF(4*zeta(3))                                                        # needs sage.symbolic
            [4.8082276126383...]
            sage: RBF(exp(1), 0.01)                                                     # needs sage.symbolic
            [2.7 +/- ...]

        TESTS:

        The following conversions used to yield incorrect results::

            sage: RBF(airy_ai(1))                                                       # needs sage.symbolic
            [0.135292416312881...]
            sage: v = RBF(zetaderiv(1, 3/2)); v                                         # needs sage.symbolic
            [-3.932239737431101 +/- 5.58e-16]
            sage: v.overlaps(RealBallField(100)(3/2).zetaderiv(1))                      # needs sage.symbolic
            True
        """
        # Symbolic expressions are handled in a special way, see
        # Expression._arb_(). A call like RBF(expr, rad) converts expr to a
        # ball using its _arb_() method and sends us the result to adjust the
        # radius thanks to the general mechanism implemented in
        # NamedConvertMap.
        try:
            return self.element_class(self, mid, rad)
        except TypeError:
            pass
        try:
            _mid = RealIntervalField(self._prec)(mid)
            return self.element_class(self, _mid, rad)
        except (TypeError, ValueError):
            pass
        raise TypeError(lazy_string("unable to convert %s to a RealBall", mid))

    def _repr_option(self, key):
        """
        Declare that real balls print atomically.

        TESTS::

            sage: RBF._repr_option('element_is_atomic')
            True
            sage: RBF['x']([-2,-2,-2/3])
            [-0.666666666666667 +/- ...e-16]*x^2 - 2.000000000000000*x
            - 2.000000000000000
            sage: RBF._repr_option('element_is_atomic_typo')
            Traceback (most recent call last):
            ...
            KeyError: 'element_is_atomic_typo'
        """
        if key == 'element_is_atomic':
            return True

        return super()._repr_option(key)

    def gens(self) -> tuple:
        r"""
        EXAMPLES::

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
            (Completion[+Infinity, prec=42], Rational Field)
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
            Complex ball field with 53 bits of precision
            sage: RealBallField(3).algebraic_closure()
            Complex ball field with 3 bits of precision
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

    prec = precision

    def is_exact(self):
        """
        Real ball fields are not exact.

        EXAMPLES::

            sage: RealBallField().is_exact()
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
            [0, 1.000000000000000, [0.3333333333333333 +/- ...e-17],
            [-4.733045976388941e+363922934236666733021124 +/- ...e+363922934236666733021108],
            [+/- inf], [+/- inf], [+/- inf], nan]
        """
        inf = self(sage.rings.infinity.Infinity)
        elements = [self(0), self(1), self(1)/3, -self(2)**(Integer(2)**80),
                    inf, -inf, self.zero().add_error(inf)]
        try:
            from sage.symbolic.constants import NotANumber
        except ImportError:
            pass
        else:
            elements.append(self.element_class(self, NotANumber()))
        return elements

    def _sum_of_products(self, terms):
        r"""
        Compute a sum of product of real balls without creating temporary
        Python objects

        The input objects should be real balls, but need not belong to this
        parent. The computation is performed at the precision of this parent.

        EXAMPLES::

            sage: Pol.<x> = RealBallField(1000)[]
            sage: pol = (x + 1/3)^100
            sage: RBF._sum_of_products((c, c) for c in pol)
            [6.3308767660842e+23 +/- ...e+9]

        TESTS::

            sage: RBF._sum_of_products([])
            0
            sage: RBF._sum_of_products([[]])
            1.000000000000000
            sage: RBF._sum_of_products([["a"]])
            Traceback (most recent call last):
            ...
            TypeError: Cannot convert str to sage.rings.real_arb.RealBall
        """
        cdef RealBall res = RealBall.__new__(RealBall)
        cdef RealBall factor
        cdef arb_t tmp
        res._parent = self
        arb_zero(res.value)
        arb_init(tmp)
        try:
            for term in terms:
                arb_one(tmp)
                for factor in term:
                    arb_mul(tmp, tmp, factor.value, self._prec)
                arb_add(res.value, res.value, tmp, self._prec)
        finally:
            arb_clear(tmp)
        return res
    # Constants

    def pi(self):
        r"""
        Return a ball enclosing `\pi`.

        EXAMPLES::

            sage: RBF.pi()
            [3.141592653589793 +/- ...e-16]
            sage: RealBallField(128).pi()
            [3.1415926535897932384626433832795028842 +/- ...e-38]
        """
        cdef RealBall res = RealBall.__new__(RealBall)
        res._parent = self
        if _do_sig(self._prec): sig_on()
        arb_const_pi(res.value, self._prec)
        if _do_sig(self._prec): sig_off()
        return res

    def log2(self):
        r"""
        Return a ball enclosing `\log(2)`.

        EXAMPLES::

            sage: RBF.log2()
            [0.6931471805599453 +/- ...e-17]
            sage: RealBallField(128).log2()
            [0.69314718055994530941723212145817656807 +/- ...e-39]
        """
        cdef RealBall res = RealBall.__new__(RealBall)
        res._parent = self
        if _do_sig(self._prec): sig_on()
        arb_const_log2(res.value, self._prec)
        if _do_sig(self._prec): sig_off()
        return res

    def euler_constant(self):
        r"""
        Return a ball enclosing the Euler constant.

        EXAMPLES::

            sage: RBF.euler_constant() # abs tol 1e-15
            [0.5772156649015329 +/- 9.00e-17]
            sage: RealBallField(128).euler_constant()
            [0.57721566490153286060651209008240243104 +/- ...e-39]
        """
        cdef RealBall res = RealBall.__new__(RealBall)
        res._parent = self
        if _do_sig(self._prec): sig_on()
        arb_const_euler(res.value, self._prec)
        if _do_sig(self._prec): sig_off()
        return res

    def catalan_constant(self):
        r"""
        Return a ball enclosing the Catalan constant.

        EXAMPLES::

            sage: RBF.catalan_constant()
            [0.915965594177219 +/- ...e-16]
            sage: RealBallField(128).catalan_constant()
            [0.91596559417721901505460351493238411077 +/- ...e-39]
        """
        cdef RealBall res = RealBall.__new__(RealBall)
        res._parent = self
        if _do_sig(self._prec): sig_on()
        arb_const_catalan(res.value, self._prec)
        if _do_sig(self._prec): sig_off()
        return res

    # Ball functions of non-ball arguments

    def sinpi(self, x):
        r"""
        Return a ball enclosing `\sin(\pi x)`.

        This works even if ``x`` itself is not a ball, and may be faster or
        more accurate where ``x`` is a rational number.

        EXAMPLES::

            sage: RBF.sinpi(1)
            0
            sage: RBF.sinpi(1/3)
            [0.866025403784439 +/- ...e-16]
            sage: RBF.sinpi(1 + 2^(-100))
            [-2.478279624546525e-30 +/- ...e-46]

        .. SEEALSO:: :meth:`~sage.rings.real_arb.RealBall.sin`

        TESTS::

            sage: RBF.sinpi(RLF(sqrt(2)))                                               # needs sage.symbolic
            [-0.963902532849877 +/- ...e-16]
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
        r"""
        Return a ball enclosing `\cos(\pi x)`.

        This works even if ``x`` itself is not a ball, and may be faster or
        more accurate where ``x`` is a rational number.

        EXAMPLES::

            sage: RBF.cospi(1)
            -1.000000000000000
            sage: RBF.cospi(1/3)
            0.5000000000000000

        .. SEEALSO:: :meth:`~sage.rings.real_arb.RealBall.cos`

        TESTS::

            sage: RBF.cospi(RLF(sqrt(2)))                                               # needs sage.symbolic
            [-0.26625534204142 +/- ...e-15]
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
            [1.932849514310098...+1956570551809674817225 +/- ...]
            sage: RBF.gamma(1/3)
            [2.678938534707747 +/- ...e-16]
            sage: RBF.gamma(-5)
            nan

        .. SEEALSO:: :meth:`~sage.rings.real_arb.RealBall.gamma`

        TESTS::

            sage: RBF.gamma(RLF(pi))  # abs tol 1e-13                                   # needs sage.symbolic
            [2.28803779534003 +/- 4.12e-15]
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

            sage: RBF.zeta(3)
            [1.202056903159594 +/- ...e-16]
            sage: RBF.zeta(1)
            nan
            sage: RBF.zeta(1/2)
            [-1.460354508809587 +/- ...e-16]

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
            [1.000000000000000, -0.5000000000000000, [0.1666666666666667 +/- ...e-17], 0]
            sage: RBF.bernoulli(2**20)
            [-1.823002872104961e+5020717 +/- ...e+5020701]
            sage: RBF.bernoulli(2**1000)
            Traceback (most recent call last):
            ...
            ValueError: argument too large

        TESTS::

            sage: RBF.bernoulli(2r)
            [0.1666666666666667 +/- ...e-17]
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

            sage: [RBF.fibonacci(n) for n in range(7)]
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
            [3.78202087472056e+20898764024997873376 +/- ...e+20898764024997873361]
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

            sage: [RBF.bell_number(n) for n in range(7)]
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
            [5.38270113176282e+1794956117137290721328 +/- ...e+1794956117137290721313]
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
            [1.448372990...e+2928836 +/- ...]
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

        OUTPUT: integer

        EXAMPLES::

            sage: RBF.maximal_accuracy()
            9223372036854775807 # 64-bit
            2147483647          # 32-bit

        .. SEEALSO:: :meth:`RealBall.accuracy`
        """
        return ARF_PREC_EXACT


cdef inline bint _do_sig(long prec) noexcept:
    """
    Whether signal handlers should be installed for calls to FLINT.

    TESTS::

        sage: _ = RealBallField()(1).psi() # indirect doctest
        sage: _ = RealBallField(1500)(1).psi()
    """
    return (prec > 1000)

cdef inline long prec(RealBall ball) noexcept:
    return ball._parent._prec


def create_RealBall(parent, serialized):
    r"""
    Create a RealBall from a serialized representation.

    TESTS::

        sage: from sage.rings.real_arb import create_RealBall
        sage: create_RealBall(RBF, b'15555555555555 -36 1 -36')
        [0.3333333333333333 +/- 7.04e-17]
        sage: create_RealBall(RBF, b'foo')
        Traceback (most recent call last):
        ...
        ValueError: incorrect format
    """
    cdef RealBall res = RealBall.__new__(RealBall)
    res._parent = parent
    sig_on()
    cdef bint error = arb_load_str(res.value, serialized)
    sig_off()
    if error:
        raise ValueError("incorrect format")
    else:
        return res


cdef class RealBall(RingElement):
    """
    Hold one ``arb_t``.

    EXAMPLES::

        sage: a = RealBallField()(RIF(1))                     # indirect doctest
        sage: b = a.psi()
        sage: b # abs tol 1e-15
        [-0.5772156649015329 +/- 4.84e-17]
        sage: RIF(b)
        -0.577215664901533?
    """

    def __cinit__(self):
        """
        Initialize the parent and allocate memory.

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

        - ``parent`` -- a :class:`RealBallField`

        - ``mid`` -- (optional) ball midpoint, see examples below. If omitted,
          initialize the ball to zero, ignoring the ``rad`` argument.

        - ``rad`` -- (optional) a :class:`RealNumber` or a Python float, ball
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
            [0.3333333333333333 +/- ...e-17]
            sage: RBF(3.14)
            [3.140000000000000 +/- ...e-16]

        ::

            sage: RBF(3, 0.125)
            [3e+0 +/- 0.126]
            sage: RBF(pi, 0.125r)                                                       # needs sage.symbolic
            [3e+0 +/- 0.267]
            sage: RBF(3, 1/8)
            [3e+0 +/- 0.126]
            sage: RBF(13, 1)
            [1e+1 +/- 4.01]

        ::

            sage: NF.<sqrt2> = QuadraticField(2)                                        # needs sage.rings.number_field
            sage: RBF(1/5 + sqrt2/2)                                                    # needs sage.rings.number_field
            [0.907106781186547 +/- ...e-16]

        Note that integers and floating-point numbers are ''not'' rounded to
        the parent's precision::

            sage: b = RBF(11111111111111111111111111111111111111111111111); b
            [1.111111111111111e+46 +/- ...e+30]
            sage: b.mid().exact_rational()
            11111111111111111111111111111111111111111111111

        Similarly, converting a real ball from one real ball field to another
        (with a different precision) only changes the way it is displayed and
        the precision of operations involving it, not the actual representation
        of its center::

            sage: RBF100 = RealBallField(100)
            sage: b100 = RBF100(1/3); b100
            [0.333333333333333333333333333333 +/- ...e-31]
            sage: b53 = RBF(b100); b53
            [0.3333333333333333 +/- ...e-17]
            sage: RBF100(b53)
            [0.333333333333333333333333333333 +/- ...e-31]

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
            [1.100000000000000 +/- ...e-16]
            sage: RBF(str(RBF("1.1")))
            [1.100000000000000 +/- ...e-16]
            sage: RBF("3.25")
            3.250000000000000
            sage: RBF("-3.1 +/- 1e-10")
            [-3.100000000 +/- ...e-10]
            sage: RBF("[+/-1]")
            [+/- 1.01]
            sage: RBF("inf +/- inf")
            [+/- inf]

        .. SEEALSO:: :meth:`RealBallField._element_constructor_`

        TESTS::

            sage: from sage.rings.real_arb import RealBall
            sage: RealBall(RBF, sage.symbolic.constants.Pi())
            [3.141592653589793 +/- ...e-16]
            sage: RealBall(RBF, sage.symbolic.constants.Log2())
            [0.6931471805599453 +/- ...e-17]
            sage: RealBall(RBF, sage.symbolic.constants.Catalan())
            [0.915965594177219 +/- ...e-16]
            sage: RealBall(RBF, sage.symbolic.constants.Khinchin())
            [2.685452001065306 +/- ...e-16]
            sage: RealBall(RBF, sage.symbolic.constants.Glaisher())
            [1.282427129100623 +/- ...e-16]
            sage: RealBall(RBF, sage.symbolic.constants.e)
            [2.718281828459045 +/- ...e-16]
            sage: RealBall(RBF, sage.symbolic.constants.EulerGamma()) # abs tol 1e-15
            [0.5772156649015329 +/- 9.00e-17]
            sage: RBF("1 +/- 0.001")
            [1.00 +/- ...e-3]
            sage: RBF("2.3e10000000000000000000000 +/- 0.00005e10000000000000000000000")
            [2.3000e+10000000000000000000000 +/- ...e+9999999999999999999995]
            sage: RBF("0.3 +/- 0.2 +/- 0.1")
            Traceback (most recent call last):
            ...
            ValueError: unsupported string format

            sage: NF.<a> = QuadraticField(2, embedding=AA(2).sqrt())
            sage: RBF.coerce(a)
            [1.414213562373095 +/- ...e-16]
            sage: NF.<a> = QuadraticField(2, embedding=-AA(2).sqrt())
            sage: RBF.coerce(a)
            [-1.414213562373095 +/- ...e-16]
            sage: NF.<a> = QuadraticField(2, embedding=None)
            sage: RBF.coerce(a)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion ...
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
        cdef fmpz_t tmpz
        cdef fmpq_t tmpq
        cdef arf_t  tmpr
        cdef mag_t  tmpm

        Element.__init__(self, parent)

        if mid is None:
            return

        elif isinstance(mid, RealBall):
            arb_set(self.value, (<RealBall> mid).value) # no rounding!
        elif is_small_python_int(mid):
            arb_set_si(self.value, PyLong_AsLong(mid)) # no rounding!
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
        elif isinstance(mid, RealIntervalFieldElement):
            mpfi_to_arb(self.value,
                (<RealIntervalFieldElement> mid).value,
                prec(self))
        elif isinstance(mid, str):
            if arb_set_str(self.value, str_to_bytes(mid), prec(self)) != 0:
                raise ValueError("unsupported string format")
        else:
            # the initializers that trigger imports
            import sage.rings.infinity
            if isinstance(mid, sage.rings.infinity.AnInfinity):
                if isinstance(mid, sage.rings.infinity.PlusInfinity):
                    arb_pos_inf(self.value)
                elif isinstance(mid, sage.rings.infinity.MinusInfinity):
                    arb_neg_inf(self.value)
                else:
                    arb_zero_pm_inf(self.value)
            else:
                try:
                    import sage.symbolic.constants
                    import sage.symbolic.expression
                except ImportError:
                    raise TypeError("unsupported midpoint type")
                if isinstance(mid, sage.symbolic.constants.Constant):
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
                elif isinstance(mid, sage.symbolic.expression.E):
                    if _do_sig(prec(self)): sig_on()
                    arb_const_e(self.value, prec(self))
                    if _do_sig(prec(self)): sig_off()
                else:
                    raise TypeError("unsupported midpoint type")

        if rad is not None:
            mag_init(tmpm)
            if isinstance(rad, RealNumber):
                arf_init(tmpr)
                arf_set_mpfr(tmpr, (<RealNumber> rad).value)
                arf_get_mag(tmpm, tmpr)
                arf_clear(tmpr)
            elif isinstance(rad, Integer):
                arf_init(tmpr)
                arf_set_mpz(tmpr, (<Integer> rad).value)
                arf_get_mag(tmpm, tmpr)
                arf_clear(tmpr)
            elif isinstance(rad, Rational):
                arf_init(tmpr)
                arf_set_mpz(tmpr, (<Integer> rad.numerator()).value)
                fmpz_init(tmpz)
                fmpz_set_mpz(tmpz, (<Integer> rad.denominator()).value)
                arf_div_fmpz(tmpr, tmpr, tmpz, prec(self), ARF_RND_UP)
                arf_get_mag(tmpm, tmpr)
                arf_clear(tmpr)
            elif isinstance(rad, float):
                mag_set_d(tmpm, PyFloat_AS_DOUBLE(rad))
            else:
                raise TypeError("rad should be a RealNumber or a Python float")
            mag_add(arb_radref(self.value), arb_radref(self.value), tmpm)
            mag_clear(tmpm)

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

        OUTPUT: string

        EXAMPLES::

           sage: RealBallField()(RIF(1.9, 2))
           [2e+0 +/- 0.101]
        """
        cdef char* c_result

        c_result = arb_get_str(self.value, (prec(self) * 31) // 100, 0)
        try:
            py_string = char_to_str(c_result)
        finally:
            flint_free(c_result)

        return py_string

    def __reduce__(self):
        r"""
        Serialize a RealBall.

        TESTS::

            sage: [loads(dumps(b)).identical(b)
            ....:  for b in [RealBallField(60).pi(), RBF(infinity)]]
            [True, True]
            sage: b = RBF(NaN); loads(dumps(b)).identical(b)                            # needs sage.symbolic
            True
        """
        cdef bytes py_val
        sig_on()
        cdef char* c_val = arb_dump_str(self.value)
        sig_off()
        try:
            py_val = <bytes> c_val
        finally:
            flint_free(c_val)
        return create_RealBall, (self._parent, py_val)

    # Conversions

    cpdef RealIntervalFieldElement _real_mpfi_(self, RealIntervalField_class parent):
        """
        Return a :mod:`real interval <sage.rings.real_mpfi>` containing this ball.

        OUTPUT: a :class:`~sage.rings.real_mpfi.RealIntervalFieldElement`

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
            sage: ZZ(RBF(pi))                                                           # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: [3.141592653589793 +/- ...e-16] does not contain a unique integer
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
            ValueError: [0.3333333333333333 +/- ...e-17] does not contain a unique rational number
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

            sage: # needs sage.symbolic
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

    def __float__(self):
        """
        Convert ``self`` to a ``float``.

        EXAMPLES::

            sage: float(RBF(1))
            1.0
        """
        return float(self.n(prec(self)))

    def __complex__(self):
        """
        Convert ``self`` to a ``complex``.

        EXAMPLES::

            sage: complex(RBF(1))
            (1+0j)
        """
        return complex(self.n(prec(self)))

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
            +infinity

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
            [0.3333 +/- ...e-5]
            sage: mid.is_exact()
            True
            sage: mid.parent()
            Real ball field with 16 bits of precision

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
            [5.55111512e-17 +/- ...e-26]
            sage: rad.is_exact()
            True
            sage: rad.parent()
            Real ball field with 30 bits of precision

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
            [0.3333333333333333 +/- ...e-17]
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

        - ``test_zero`` -- boolean (default: ``False``); if ``True``,
          make sure that the returned lower bound is positive, raising
          an error if the ball contains zero.

        OUTPUT: a ball with zero radius

        EXAMPLES::

            sage: RealBallField(8)(1/3).below_abs()
            [0.33 +/- ...e-5]
            sage: b = RealBallField(8)(1/3).below_abs()
            sage: b
            [0.33 +/- ...e-5]
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

        OUTPUT: a ball with zero radius

        EXAMPLES::

            sage: b = RealBallField(8)(1/3).above_abs()
            sage: b
            [0.33 +/- ...e-3]
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

        - ``rnd`` -- string; rounding mode for the parent of the result (does
          not affect its value!), see
          :meth:`sage.rings.real_mpfi.RealIntervalFieldElement.upper`

        OUTPUT: a real number

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

        - ``rnd`` -- string; rounding mode for the parent of the result (does
          not affect its value!), see
          :meth:`sage.rings.real_mpfi.RealIntervalFieldElement.lower`

        OUTPUT: a real number

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

        - ``rnd`` -- string; rounding mode for the parent of the resulting
          floating-point numbers (does not affect their values!), see
          :meth:`sage.rings.real_mpfi.RealIntervalFieldElement.upper`

        OUTPUT: a pair of real numbers

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
            (-9.31322574615479e-10, 1.00000000093133)
        """
        cdef RealBall my_other = self._parent.coerce(other)
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_union(res.value, self.value, my_other.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def real(self):
        r"""
        Return the real part of this ball.

        EXAMPLES::

            sage: RBF(1/3).real()
            [0.3333333333333333 +/- 7.04e-17]
        """
        return self

    def imag(self):
        r"""
        Return the imaginary part of this ball.

        EXAMPLES::

            sage: RBF(1/3).imag()
            0
        """
        return self._parent.zero()

    # Precision and accuracy

    def nbits(self):
        r"""
        Return the minimum precision sufficient to represent this ball exactly.

        In other words, return the number of bits needed to represent the
        absolute value of the mantissa of the midpoint of this ball. The result
        is 0 if the midpoint is a special value.

        EXAMPLES::

            sage: RBF(1/3).nbits()
            53
            sage: RBF(1023, .1).nbits()
            10
            sage: RBF(1024, .1).nbits()
            1
            sage: RBF(0).nbits()
            0
            sage: RBF(infinity).nbits()
            0
        """
        return arb_bits(self.value)

    def round(self):
        """
        Return a copy of this ball with center rounded to the precision of the
        parent.

        EXAMPLES:

        It is possible to create balls whose midpoint is more precise than
        their parent's nominal precision (see :mod:`~sage.rings.real_arb` for
        more information)::

            sage: b = RBF(pi.n(100))                                                    # needs sage.symbolic
            sage: b.mid()                                                               # needs sage.symbolic
            3.141592653589793238462643383

        The ``round()`` method rounds such a ball to its parent's precision::

            sage: b.round().mid()                                                       # needs sage.symbolic
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

            sage: RBF(pi).accuracy()                                                    # needs sage.symbolic
            52
            sage: RBF(1).accuracy() == RBF.maximal_accuracy()
            True
            sage: RBF(NaN).accuracy() == -RBF.maximal_accuracy()                        # needs sage.symbolic
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

        - ``ampl`` -- a real ball (or an object that can be coerced to a real
          ball)

        OUTPUT: a new real ball

        EXAMPLES::

            sage: err = RBF(10^-16)
            sage: RBF(1).add_error(err)
            [1.000000000000000 +/- ...e-16]

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
            sage: RBF(pi).is_nonzero()                                                  # needs sage.symbolic
            True
            sage: RBF(RIF(-0.5, 0.5)).is_nonzero()
            False

        .. SEEALSO:: :meth:`is_zero`
        """
        return arb_is_nonzero(self.value)

    def __bool__(self):
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
            sage: extended_line = RBF(0).add_error(infinity)
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

        lt = left
        rt = right

        if op == Py_EQ:
            return arb_eq(lt.value, rt.value)
        if op == Py_NE:
            return arb_ne(lt.value, rt.value)
        if op == Py_GT:
            return arb_gt(lt.value, rt.value)
        if op == Py_LT:
            return arb_lt(lt.value, rt.value)
        if op == Py_GE:
            return arb_ge(lt.value, rt.value)
        if op == Py_LE:
            return arb_le(lt.value, rt.value)

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
            [0.3333333333333333 +/- ...e-17]

            sage: RBF('nan').min(0)
            nan

        .. SEEALSO:: :meth:`max`

        TESTS::

            sage: RBF(0).min()
            0
            sage: RBF(infinity).min().rad()
            0.00000000
        """
        cdef RealBall res = self._new()
        cdef long p = prec(self)
        arb_set(res.value, self.value)
        for b in others:
            if not isinstance(b, RealBall):
                b = self._parent.coerce(b)
            arb_min(res.value, res.value, (<RealBall> b).value, p)
        return res

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
            [0.3333333333333333 +/- ...e-17]

            sage: RBF('nan').max(0)
            nan

        .. SEEALSO:: :meth:`min`

        TESTS::

            sage: RBF(0).max()
            0
        """
        cdef RealBall res = self._new()
        cdef long p = prec(self)
        arb_set(res.value, self.value)
        for b in others:
            if not isinstance(b, RealBall):
                b = self._parent.coerce(b)
            arb_max(res.value, res.value, (<RealBall> b).value, p)
        return res

    def is_finite(self):
        """
        Return ``True`` iff the midpoint and radius of this ball are both
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
        Return ``True`` iff ``self`` and ``other`` are equal as balls, i.e.
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
        Return ``True`` iff ``self`` and ``other`` have some point in common.

        If either ``self`` or ``other`` contains NaN, this method always
        returns nonzero (as a NaN could be anything, it could in particular
        contain any number that is included in the other operand).

        EXAMPLES::

            sage: RBF(pi).overlaps(RBF(pi) + 2**(-100))                                 # needs sage.symbolic
            True
            sage: RBF(pi).overlaps(RBF(3))                                              # needs sage.symbolic
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
            sage: RBF(sqrt(2)).contains_exact(sqrt(2))                                  # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: unsupported type: <class 'sage.symbolic.expression.Expression'>

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
            elif is_small_python_int(other):
                res = arb_contains_si(self.value, PyLong_AsLong(other))
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
        Return ``True`` if ``other`` can be verified to be contained in ``self``.

        The test is done using interval arithmetic with a precision determined
        by the parent of ``self`` and may return false negatives.

        EXAMPLES::

            sage: sqrt(2) in RBF(sqrt(2))                                               # needs sage.symbolic
            True

        A false negative::

            sage: sqrt(2) in RBF(RealBallField(100)(sqrt(2)))                           # needs sage.symbolic
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

            Contrary to the usual convention, a return value of ``True`` does
            not imply that all points of the ball satisfy the predicate.
            This is due to the way comparisons with symbolic infinities work in
            sage.

        EXAMPLES::

            sage: RBF(infinity).is_infinity()
            True
            sage: RBF(-infinity).is_infinity()
            True
            sage: RBF(NaN).is_infinity()                                                # needs sage.symbolic
            True
            sage: (~RBF(0)).is_infinity()
            True
            sage: RBF(42, rad=1.r).is_infinity()
            False
        """
        return not self.is_finite()

    def is_NaN(self):
        """
        Return ``True`` if this ball is not-a-number.

        EXAMPLES::

            sage: RBF(NaN).is_NaN()                                                     # needs sage.symbolic
            True
            sage: RBF(-5).gamma().is_NaN()
            True
            sage: RBF(infinity).is_NaN()
            False
            sage: RBF(42, rad=1.r).is_NaN()
            False
        """
        return arf_is_nan(arb_midref(self.value))

    # Arithmetic

    def __neg__(self):
        """
        Return the opposite of this ball.

        EXAMPLES::

            sage: -RBF(1/3)
            [-0.3333333333333333 +/- ...e-17]
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
            [0.2000000000000000 +/- ...e-17]
            sage: ~RBF(0)
            nan
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
            [1.333333333333333 +/- ...e-16]
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
            [0.666666666666667 +/- ...e-16]
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
            [-0.666666666666667 +/- ...e-16]
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

            sage: RBF(pi)/RBF(e)                                                        # needs sage.symbolic
            [1.155727349790922 +/- ...e-16]
            sage: RBF(2)/RBF(0)
            nan
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_div(res.value, self.value, (<RealBall> other).value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def __pow__(base, expo, _):
        """
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: RBF(e)^17
            [24154952.7535753 +/- ...e-8]
            sage: RBF(e)^(-1)
            [0.367879441171442 +/- ...e-16]
            sage: RBF(e)^(1/2)
            [1.648721270700128 +/- ...e-16]
            sage: RBF(e)^RBF(pi)
            [23.1406926327793 +/- ...e-14]

        ::

            sage: RBF(-1)^(1/3)
            nan
            sage: RBF(0)^(-1)
            nan
            sage: RBF(-e)**RBF(pi)                                                      # needs sage.symbolic
            nan

        TESTS::

            sage: RBF(e)**(2r)                                                          # needs sage.symbolic
            [7.38905609893065 +/- ...e-15]
            sage: RBF(e)**(-1r)                                                         # needs sage.symbolic
            [0.367879441171442 +/- ...e-16]
        """
        cdef fmpz_t tmpz
        if not isinstance(base, RealBall):
            return sage.structure.element.bin_op(base, expo, operator.pow)
        cdef RealBall self = base
        cdef RealBall res = self._new()
        if is_small_python_int(expo) and expo > 0:
            if _do_sig(prec(self)): sig_on()
            arb_pow_ui(res.value, self.value, PyLong_AsLong(expo), prec(self))
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
            [1.414213562373095 +/- ...e-16]
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
            [1.414213562373095 +/- ...e-16]
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
            [0.707106781186547 +/- ...e-16]
            sage: RBF(0).rsqrt()
            nan
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_rsqrt(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def sqrt1pm1(self):
        r"""
        Return `\sqrt{1+\mathrm{self}}-1`, computed accurately when ``self`` is
        close to zero.

        EXAMPLES::

            sage: eps = RBF(10^(-20))
            sage: (1 + eps).sqrt() - 1
            [+/- ...e-16]
            sage: eps.sqrt1pm1()
            [5.00000000000000e-21 +/- ...e-36]
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

        - ``shift`` -- integer; may be negative

        EXAMPLES::

            sage: RBF(1/3) << 2 # indirect doctest
            [1.333333333333333 +/- ...e-16]
            sage: RBF(1) << -1
            0.5000000000000000

        TESTS::

            sage: RBF(1) << (2^100)
            [2.285367694229514e+381600854690147056244358827360 +/- ...e+381600854690147056244358827344]
            sage: RBF(1) << (-2^100)
            [4.375663498372584e-381600854690147056244358827361 +/- ...e-381600854690147056244358827378]

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
        if is_small_python_int(shift):
            arb_mul_2exp_si(res.value, self.value, PyLong_AsLong(shift))
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

        - ``shift`` -- integer; may be negative

        EXAMPLES::

            sage: RBF(4) >> 2
            1.000000000000000
            sage: RBF(1/3) >> -2
            [1.333333333333333 +/- ...e-16]

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

        - ``base`` -- (optional) positive real ball or number; if ``None``,
          return the natural logarithm ``ln(self)``, otherwise, return the
          general logarithm ``ln(self)/ln(base)``

        EXAMPLES::

            sage: RBF(3).log()
            [1.098612288668110 +/- ...e-16]
            sage: RBF(3).log(2)
            [1.58496250072116 +/- ...e-15]
            sage: log(RBF(5), 2)
            [2.32192809488736 +/- ...e-15]

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
            [+/- ...e-16]
            sage: eps.log1p()
            [1.00000000000000e-30 +/- ...e-46]
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
            [2.718281828459045 +/- ...e-16]
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
            [+/- ...e-30]
            sage: eps.expm1()
            [1.000000000000000e-30 +/- ...e-47]
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

            sage: RBF(pi).sin()                                                         # needs sage.symbolic
            [+/- ...e-16]

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

            sage: RBF(pi).cos()                                                         # needs sage.symbolic
            [-1.00000000000000 +/- ...e-16]

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
            [1.557407724654902 +/- ...e-16]
            sage: RBF(pi/2).tan()                                                       # needs sage.symbolic
            nan
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
            [0.642092615934331 +/- ...e-16]
            sage: RBF(pi).cot()                                                         # needs sage.symbolic
            nan
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_cot(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def sec(self):
        """
        Return the secant of this ball.

        EXAMPLES::

            sage: RBF(1).sec()
            [1.850815717680925 +/- ...e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_sec(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def csc(self):
        """
        Return the cosecant of this ball.

        EXAMPLES::

            sage: RBF(1).csc()
            [1.188395105778121 +/- ...e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_csc(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arcsin(self):
        """
        Return the arcsine of this ball.

        EXAMPLES::

            sage: RBF(1).arcsin()
            [1.570796326794897 +/- ...e-16]
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
            [0.7853981633974483 +/- ...e-17]
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
            [1.175201193643801 +/- ...e-16]
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
            [1.543080634815244 +/- ...e-16]
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
            [0.761594155955765 +/- ...e-16]
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
            [1.313035285499331 +/- ...e-16]
            sage: RBF(0).coth()
            nan
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_coth(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def sech(self):
        """
        Return the hyperbolic secant of this ball.

        EXAMPLES::

            sage: RBF(1).sech()
            [0.648054273663885 +/- ...e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_sech(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def csch(self):
        """
        Return the hyperbolic cosecant of this ball.

        EXAMPLES::

            sage: RBF(1).csch()
            [0.850918128239321 +/- ...e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_csch(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def arcsinh(self):
        """
        Return the inverse hyperbolic sine of this ball.

        EXAMPLES::

            sage: RBF(1).arcsinh()
            [0.881373587019543 +/- ...e-16]
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
            [1.316957896924817 +/- ...e-16]
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
            [0.549306144334055 +/- ...e-16]
            sage: RBF(1).arctanh()
            nan
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_atanh(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    # Special functions

    def erf(self):
        """
        Error function.

        EXAMPLES::

            sage: RBF(1/2).erf() # abs tol 1e-16
            [0.520499877813047 +/- 6.10e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_hypgeom_erf(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def erfi(self):
        """
        Imaginary error function.

        EXAMPLES::

            sage: RBF(1/2).erfi()
            [0.614952094696511 +/- 2.22e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_hypgeom_erfi(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def Ei(self):
        """
        Exponential integral.

        EXAMPLES::

            sage: RBF(1).Ei()  # abs tol 5e-16
            [1.89511781635594 +/- 4.94e-15]

        TESTS::

            sage: RBF(Ei(1))  # abs tol 5e-16                                           # needs sage.symbolic
            [1.89511781635594 +/- 4.94e-15]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_hypgeom_ei(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def Si(self):
        """
        Sine integral.

        EXAMPLES::

            sage: RBF(1).Si() # abs tol 1e-15
            [0.946083070367183 +/- 9.22e-16]

        TESTS::

            sage: RBF(Si(1))  # abs tol 1e-15                                           # needs sage.symbolic
            [0.946083070367183 +/- 9.22e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_hypgeom_si(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    sin_integral = Si

    def Ci(self):
        """
        Cosine integral.

        EXAMPLES::

            sage: RBF(1).Ci()  # abs tol 5e-16
            [0.337403922900968 +/- 3.25e-16]

        TESTS::

            sage: RBF(Ci(1))  # abs tol 5e-16                                           # needs sage.symbolic
            [0.337403922900968 +/- 3.25e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_hypgeom_ci(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cos_integral = Ci

    def Shi(self):
        """
        Hyperbolic sine integral.

        EXAMPLES::

            sage: RBF(1).Shi()
            [1.05725087537573 +/- 2.77e-15]

        TESTS::

            sage: RBF(Shi(1))                                                           # needs sage.symbolic
            [1.05725087537573 +/- 2.77e-15]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_hypgeom_shi(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    sinh_integral = Shi

    def Chi(self):
        """
        Hyperbolic cosine integral.

        EXAMPLES::

            sage: RBF(1).Chi()  # abs tol 1e-17
            [0.837866940980208 +/- 4.72e-16]

        TESTS::

            sage: RBF(Chi(1))  # abs tol 1e-17                                          # needs sage.symbolic
            [0.837866940980208 +/- 4.72e-16]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_hypgeom_chi(res.value, self.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    cosh_integral = Chi

    def li(self):
        """
        Logarithmic integral.

        EXAMPLES::

            sage: RBF(3).li()  # abs tol 1e-15
            [2.16358859466719 +/- 4.72e-15]

        TESTS::

            sage: RBF(li(0))                                                            # needs sage.symbolic
            0
            sage: RBF(Li(0))                                                            # needs sage.symbolic
            [-1.04516378011749 +/- 4.23e-15]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_hypgeom_li(res.value, self.value, False, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    log_integral = li

    def Li(self):
        """
        Offset logarithmic integral.

        EXAMPLES::

            sage: RBF(3).Li()  # abs tol 1e-15
            [1.11842481454970 +/- 7.61e-15]
        """
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        arb_hypgeom_li(res.value, self.value, True, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    log_integral_offset = Li

    def beta(self, a, z=1):
        """
        (Incomplete) beta function.

        INPUT:

        - ``a``, ``z`` -- (optional) real balls

        OUTPUT:

        The lower incomplete beta function `B(self, a, z)`.

        With the default value of ``z``, the complete beta function `B(self, a)`.

        EXAMPLES::

            sage: RBF(sin(3)).beta(RBF(2/3).sqrt())  # abs tol 1e-13                    # needs sage.symbolic
            [7.407661629415 +/- 1.07e-13]
            sage: RealBallField(100)(7/2).beta(1)  # abs tol 1e-30
            [0.28571428571428571428571428571 +/- 5.23e-30]
            sage: RealBallField(100)(7/2).beta(1, 1/2)
            [0.025253813613805268728601584361 +/- 2.53e-31]

        .. TODO::

            At the moment RBF(beta(a,b)) does not work, one needs
            RBF(a).beta(b) for this to work. See :issue:`32851`
            and :issue:`24641`.
        """
        cdef RealBall a_ball, z_ball
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        a_ball = self._parent.coerce(a)
        z_ball = self._parent.coerce(z)
        arb_hypgeom_beta_lower(res.value, self.value, a_ball.value,
                               z_ball.value, False, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def gamma(self, a=None):
        r"""
        Image of this ball by the (upper incomplete) Euler Gamma function.

        For `a` real, return the upper incomplete Gamma function
        `\Gamma(self,a)`.

        For integer and rational arguments,
        :meth:`~sage.rings.real_arb.RealBallField.gamma` may be faster.

        EXAMPLES::

            sage: RBF(1/2).gamma()
            [1.772453850905516 +/- ...e-16]
            sage: RBF(gamma(3/2, RBF(2).sqrt()))  # abs tol 2e-17
            [0.37118875695353 +/- 3.00e-15]
            sage: RBF(3/2).gamma_inc(RBF(2).sqrt())  # abs tol 2e-17
            [0.37118875695353 +/- 3.00e-15]

        .. SEEALSO::
            :meth:`~sage.rings.real_arb.RealBallField.gamma`

        TESTS::

            sage: RealBallField(100).gamma(1/2)
            [1.77245385090551602729816748334 +/- 1.90e-30]
        """
        cdef RealBall a_ball
        cdef RealBall res = self._new()
        if a is None:
            if _do_sig(prec(self)): sig_on()
            arb_gamma(res.value, self.value, prec(self))
            if _do_sig(prec(self)): sig_off()
        else:
            if _do_sig(prec(self)): sig_on()
            a_ball = self._parent.coerce(a)
            arb_hypgeom_gamma_upper(res.value, self.value, a_ball.value, 0, prec(self))
            if _do_sig(prec(self)): sig_off()
        return res

    gamma_inc = gamma

    def gamma_inc_lower(self, a):
        r"""
        Image of this ball by the lower incomplete Euler Gamma function.

        For `a` real, return the lower incomplete Gamma function
        of `\Gamma(self,a)`.

        EXAMPLES::

            sage: RBF(gamma_inc_lower(1/2, RBF(2).sqrt()))
            [1.608308637729248 +/- 8.14e-16]
            sage: RealBallField(100)(7/2).gamma_inc_lower(5)
            [2.6966551541863035516887949614 +/- 8.91e-29]
        """
        cdef RealBall a_ball
        cdef RealBall res = self._new()
        if _do_sig(prec(self)): sig_on()
        a_ball = RBF(a)
        arb_hypgeom_gamma_lower(res.value, self.value, a_ball.value, 0, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

    def log_gamma(self):
        """
        Return the image of this ball by the logarithmic Gamma function.

        The complex branch structure is assumed, so if ``self <= 0``, the result
        is an indeterminate interval.

        EXAMPLES::

            sage: RBF(1/2).log_gamma()
            [0.572364942924700 +/- ...e-16]
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
        r"""
        Return the ``n``-th rising factorial of this ball.

        The `n`-th rising factorial of `x` is equal to `x (x+1) \cdots (x+n-1)`.

        For real `n`, it is a quotient of gamma functions.

        EXAMPLES::

            sage: RBF(1).rising_factorial(5)
            120.0000000000000
            sage: RBF(1/2).rising_factorial(1/3) # abs tol 1e-14
            [0.636849884317974 +/- 8.98e-16]
        """
        cdef RealBall result = self._new()
        cdef RealBall my_n = self._parent.coerce(n)
        if _do_sig(prec(self)): sig_on()
        arb_rising(result.value, self.value, my_n.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return result

    cpdef RealBall psi(self):
        """
        Compute the digamma function with argument ``self``.

        EXAMPLES::

            sage: RBF(1).psi() # abs tol 1e-15
            [-0.5772156649015329 +/- 4.84e-17]
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

        Otherwise, it computes the Hurwitz zeta function.

        Use :meth:`RealBallField.zeta` to compute the Riemann zeta function of
        a small integer without first converting it to a real ball.

        EXAMPLES::

            sage: RBF(-1).zeta()
            [-0.0833333333333333 +/- ...e-17]
            sage: RBF(-1).zeta(1)
            [-0.0833333333333333 +/- ...e-17]
            sage: RBF(-1).zeta(2)
            [-1.083333333333333 +/- ...e-16]
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

    def zetaderiv(self, k):
        r"""
        Return the image of this ball by the `k`-th derivative of the Riemann
        zeta function.

        For a more flexible interface, see the low-level method
        ``_zeta_series`` of polynomials with complex ball coefficients.

        EXAMPLES::

            sage: RBF(1/2).zetaderiv(1)
            [-3.92264613920915...]
            sage: RBF(2).zetaderiv(3)
            [-6.0001458028430...]
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        Pol = PolynomialRing(self._parent.complex_field(), 'x')
        ser = Pol([self, 1])._zeta_series(k + 1)
        return ser[k].real()*ZZ.coerce(k).factorial()

    def lambert_w(self):
        r"""
        Return the image of this ball by the Lambert W function.

        EXAMPLES::

            sage: RBF(1).lambert_w()
            [0.5671432904097...]
        """
        cdef RealBall res = self._new()
        sig_on()
        arb_lambertw(res.value, self.value, 0, prec(self))
        sig_off()
        return res

    def polylog(self, s):
        r"""
        Return the polylogarithm `\operatorname{Li}_s(\mathrm{self})`.

        EXAMPLES::

            sage: polylog(0, -1)                                                        # needs sage.symbolic
            -1/2
            sage: RBF(-1).polylog(0)
            [-0.50000000000000 +/- ...e-16]
            sage: polylog(1, 1/2)                                                       # needs sage.symbolic
            -log(1/2)
            sage: RBF(1/2).polylog(1)
            [0.69314718055995 +/- ...e-15]
            sage: RBF(1/3).polylog(1/2)
            [0.44210883528067 +/- 6.7...e-15]
            sage: RBF(1/3).polylog(RLF(pi))                                             # needs sage.symbolic
            [0.34728895057225 +/- ...e-15]

        TESTS::

            sage: RBF(1/3).polylog(2r)
            [0.366213229977063 +/- ...e-16]
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

            sage: # needs sage.symbolic
            sage: RBF(pi).chebyshev_T(0)
            1.000000000000000
            sage: RBF(pi).chebyshev_T(1)
            [3.141592653589793 +/- ...e-16]
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

            sage: # needs sage.symbolic
            sage: RBF(pi).chebyshev_U(0)
            1.000000000000000
            sage: RBF(pi).chebyshev_U(1)
            [6.283185307179586 +/- ...e-16]
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
            sage: RBF(sqrt(2)).agm(1)^(-1)                                              # needs sage.symbolic
            [0.8346268416740...]
        """
        cdef RealBall other_as_ball
        cdef RealBall res = self._new()
        other_as_ball = self._parent.coerce(other)
        if _do_sig(prec(self)): sig_on()
        arb_agm(res.value, self.value, other_as_ball.value, prec(self))
        if _do_sig(prec(self)): sig_off()
        return res

RBF = RealBallField()
