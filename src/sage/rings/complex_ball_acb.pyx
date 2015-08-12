r"""
Arbitrary Precision Complex Intervals using Arb

AUTHORS:

- Clemens Heuberger (2014-10-25): Initial version.

This is a rudimentary binding to the optional `Arb library
<http://fredrikj.net/arb/>`_; it may be useful to refer to its
documentation for more details.

You may have to run ``sage -i arb`` to use the arb library.

.. WARNING::

    Identical :class:`ComplexBall` objects are understood to give
    permission for algebraic simplification. This assumption is made
    to improve performance. For example, setting ``z = x*x`` sets `z`
    to a ball enclosing the set `\{t^2 : t \in x\}` and not the
    (generally larger) set `\{tu : t \in x, u \in x\}`.

Comparison
==========

Two elements are equal if and only if they are the same object
or if both are exact and equal::

    sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
    sage: CBF = ComplexBallField() # optional - arb
    doctest:...: FutureWarning: This class/method/function is marked as experimental.
    It, its functionality or its interface might change without a formal deprecation.
    See http://trac.sagemath.org/17218 for details.
    sage: a = CBF(1, 2) # optional - arb
    sage: b = CBF(1, 2) # optional - arb
    sage: a is b # optional - arb
    False
    sage: a == b # optional - arb
    True
    sage: a = CBF(1/3, 1/5) # optional - arb
    sage: b = CBF(1/3, 1/5) # optional - arb
    sage: a.is_exact() # optional - arb
    False
    sage: b.is_exact() # optional - arb
    False
    sage: a is b # optional - arb
    False
    sage: a == b # optional - arb
    False

A ball is non-zero if and only if it does not contain zero. ::

    sage: a = CBF(RIF(-0.5, 0.5)) # optional - arb
    sage: bool(a) # optional - arb
    False
    sage: a != 0 # optional - arb
    False
    sage: b = CBF(1/3, 1/5) # optional - arb
    sage: bool(b) # optional - arb
    True
    sage: b != 0 # optional - arb
    True

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
include 'sage/ext/interrupt.pxi'
include "sage/ext/python.pxi"
include "sage/ext/stdsage.pxi"

import sage.categories.fields
from sage.libs.arb.arb cimport *
from sage.libs.arb.acb cimport *
from sage.misc.superseded import experimental
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.real_arb cimport mpfi_to_arb, arb_to_mpfi
from sage.rings.real_arb import RealBallField
from sage.structure.parent cimport Parent
from sage.structure.unique_representation import UniqueRepresentation

cdef inline bint acb_is_nonzero(const acb_t z):
    return arb_is_nonzero(&z.real) or arb_is_nonzero(&z.imag)

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
    mpfi_to_arb(&target.real, source.__re, precision)
    mpfi_to_arb(&target.imag, source.__im, precision)

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

    arb_to_mpfi(target.__re, &source.real, precision)
    arb_to_mpfi(target.__im, &source.imag, precision)
    return 0

class ComplexBallField(UniqueRepresentation, Parent):
    r"""
    An approximation of the field of complex numbers using mid-rad
    intervals, also known as balls.

    INPUT:

    - ``precision`` -- an integer `\ge 2`.

    EXAMPLES::

        sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
        sage: CBF = ComplexBallField() # optional - arb; indirect doctest
        sage: CBF(1) # optional - arb
        1.000000000000000

    TESTS::

        sage: ComplexBallField(0) # optional - arb
        Traceback (most recent call last):
        ...
        ValueError: Precision must be at least 2.
        sage: ComplexBallField(1) # optional - arb
        Traceback (most recent call last):
        ...
        ValueError: Precision must be at least 2.
    """
    Element = ComplexBall

    @staticmethod
    def __classcall__(cls, long precision=53):
        r"""
        Normalize the arguments for caching.

        TESTS::

            sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
            sage: ComplexBallField(53) is ComplexBallField() # optional - arb
            True
        """
        return super(ComplexBallField, cls).__classcall__(cls, precision)

    @experimental(17218)
    def __init__(self, precision):
        r"""
        Initialize the complex ball field.

        INPUT:

        - ``precision`` -- an integer `\ge 2`.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
            sage: CBF = ComplexBallField() # optional - arb
            sage: CBF(1) # optional - arb
            1.000000000000000
        """
        if precision < 2:
            raise ValueError("Precision must be at least 2.")
        super(ComplexBallField, self).__init__(category=[sage.categories.fields.Fields()])
        self._prec = precision
        self.RealBallField = RealBallField(precision)

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
            sage: ComplexBallField() # optional - arb
            Complex ball field with 53 bits precision
            sage: ComplexBallField(106) # optional - arb
            Complex ball field with 106 bits precision
        """
        return "Complex ball field with {} bits precision".format(self._prec)

    def _coerce_map_from_(self, S):
        r"""
        Currently, there is no coercion.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
            sage: ComplexBallField()._coerce_map_from_(CIF) # optional - arb
            False
            sage: ComplexBallField()._coerce_map_from_(SR) # optional - arb
            False
        """
        return False

    def _element_constructor_(self, *args, **kwds):
        r"""
        Construct a :class:`ComplexBall`.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField # optional: arb
            sage: CBF = ComplexBallField() # optional: arb
            sage: CBF(1) # optional: arb; indirect doctest
            1.000000000000000
        """
        return self.element_class(self, *args, **kwds)

    def _an_element_(self):
        r"""
        Construct an element.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField # optional: arb
            sage: CBF = ComplexBallField() # optional: arb
            sage: CBF._an_element_() # optional: arb; indirect doctest
            [0.3333333333333333 +/- 1.49e-17] + [0.1666666666666667 +/- 4.26e-17]*I
        """
        return self(1.0/3, 1.0/6)

    def precision(self):
        """
        Return the bit precision used for operations on elements of this field.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
            sage: ComplexBallField().precision() # optional - arb
            53
        """
        return self._prec

    def is_exact(self):
        """
        Complex ball fields are not exact.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
            sage: ComplexBallField().is_exact() # optional - arb
            False
        """
        return False

    def is_finite(self):
        """
        Complex ball fields are infinite.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
            sage: ComplexBallField().is_finite() # optional - arb
            False
        """
        return False

    def characteristic(self):
        """
        Complex ball fields have characteristic zero.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
            sage: ComplexBallField().characteristic() # optional - arb
            0
        """
        return 0


cdef inline bint _do_sig(long prec):
    """
    Whether signal handlers should be installed for calls to arb.

    TESTS::

        sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
    """
    return (prec > 1000)

cdef inline long prec(ComplexBall ball):
    return ball._parent._prec

cdef inline Parent real_ball_field(ComplexBall ball):
    return ball._parent.RealBallField

cdef class ComplexBall(Element):
    """
    Hold one ``acb_t`` of the `Arb library
    <http://fredrikj.net/arb/>`_

    EXAMPLES::

        sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
        sage: a = ComplexBallField()(1, 1) # optional - arb
        sage: a # optional - arb
        1.000000000000000 + 1.000000000000000*I
        sage: a._interval() # optional - arb
        1 + 1*I
    """
    def __cinit__(self):
        """
        Allocate memory for the encapsulated value.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
            sage: ComplexBallField(2)(0) # optional - arb; indirect doctest
            0
        """
        acb_init(self.value)

    def __dealloc__(self):
        """
        Deallocate memory of the encapsulated value.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
            sage: a = ComplexBallField(2)(0) # optional - arb; indirect doctest
            sage: del a # optional - arb
        """
        acb_clear(self.value)

    def __init__(self, parent, x=None, y=None):
        """
        Initialize the :class:`ComplexBall` using `x` and `y`.

        INPUT:

        - ``parent`` -- a :class:`ComplexBallField`.

        - ``x`` -- (default: ``None``) ``None`` or a
          :class:`~sage.rings.complex_interval.ComplexIntervalFieldElement` or
          a :class:`sage.rings.real_arb.RealBall`.

        - ``y`` -- (default: ``None``) ``None`` or a
          :class:`sage.rings.real_arb.RealBall`.

        OUTPUT:

        None.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
            sage: from sage.rings.real_arb import RealBallField # optional - arb
            sage: CBF = ComplexBallField() # optional - arb
            sage: CBF(CIF(0, 1)) # optional - arb
            1.000000000000000*I
            sage: CBF(1) # optional - arb
            1.000000000000000
            sage: CBF(1, 1) # optional - arb
            1.000000000000000 + 1.000000000000000*I
            sage: CBF(x) # optional - arb
            Traceback (most recent call last):
            ...
            TypeError: unable to convert to a ComplexIntervalFieldElement
            sage: RBF = RealBallField() # optional - arb
            sage: CBF(RBF(1/3)) # optional - arb
            [0.3333333333333333 +/- 7.04e-17]
            sage: CBF(RBF(1/3), RBF(1/6)) # optional - arb
            [0.3333333333333333 +/- 7.04e-17] + [0.1666666666666667 +/- 7.04e-17]*I
            sage: CBF(1/3) # optional - arb
            [0.333333333333333 +/- 3.99e-16]
            sage: CBF(1/3, 1/6) # optional - arb
            [0.3333333333333333 +/- 7.04e-17] + [0.1666666666666667 +/- 7.04e-17]*I
            sage: ComplexBallField(106)(1/3, 1/6) # optional - arb
            [0.33333333333333333333333333333333 +/- 6.94e-33] + [0.16666666666666666666666666666666 +/- 7.70e-33]*I
        """
        Element.__init__(self, parent)

        if x is None:
            return

        if y is None:
            # we assume x to be a complex number
            if isinstance(x, RealBall):
                arb_set(&self.value.real, (<RealBall> x).value)
                arb_set_ui(&self.value.imag, 0)
            else:
                if not isinstance(x, ComplexIntervalFieldElement):
                    try:
                        x = ComplexIntervalField(prec(self))(x)
                    except TypeError:
                        raise TypeError("unable to convert to a "
                                        "ComplexIntervalFieldElement")
                ComplexIntervalFieldElement_to_acb(self.value,
                                                   <ComplexIntervalFieldElement> x)
        else:
            if not isinstance(x, RealBall):
                try:
                    x = real_ball_field(self)(x)
                except TypeError:
                    raise TypeError("unable to convert to a "
                                    "RealBall")
            if not isinstance(y, RealBall):
                try:
                    y = real_ball_field(self)(y)
                except TypeError:
                    raise TypeError("unable to convert to a "
                                    "RealBall")
            arb_set(&self.value.real, (<RealBall> x).value)
            arb_set(&self.value.imag, (<RealBall> y).value)

    cdef ComplexBall _new(self):
        """
        Return a new complex ball element with the same parent as ``self``.
        """
        cdef ComplexBall x
        x = ComplexBall.__new__(ComplexBall)
        x._parent = self._parent
        return x

    cpdef RealBall real(self):
        """
        Return the real part of this ball.

        OUTPUT:

        A :class:`RealBall`.

        EXAMPLES::

           sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
           sage: CBF = ComplexBallField() # optional - arb
           sage: a = CBF(1/3, 1/5) # optional - arb
           sage: a.real() # optional - arb
           [0.3333333333333333 +/- 7.04e-17]
        """
        cdef RealBall r
        r = real_ball_field(self)(0)
        arb_set(r.value, &self.value.real)
        return r

    cpdef RealBall imag(self):
        """
        Return the imaginary part of this ball.

        OUTPUT:

        A :class:`RealBall`.

        EXAMPLES::

           sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
           sage: CBF = ComplexBallField() # optional - arb
           sage: a = CBF(1/3, 1/5) # optional - arb
           sage: a.imag() # optional - arb
           [0.2000000000000000 +/- 4.45e-17]
        """
        cdef RealBall r
        r = real_ball_field(self)(0)
        arb_set(r.value, &self.value.imag)
        return r

    def _repr_(self):
        """
        Return a string representation of ``self``.

        OUTPUT:

        A string.

        EXAMPLES::

           sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
           sage: CBF = ComplexBallField() # optional - arb
           sage: CBF(1/3) # optional - arb
           [0.333333333333333 +/- 3.99e-16]
           sage: CBF(0, 1/3) # optional - arb
           [0.3333333333333333 +/- 7.04e-17]*I
           sage: ComplexBallField()(1/3, 1/6)  # optional - arb
           [0.3333333333333333 +/- 7.04e-17] + [0.1666666666666667 +/- 7.04e-17]*I
        """
        if arb_is_zero(&self.value.imag):
            return self.real()._repr_()
        elif arb_is_zero(&self.value.real):
            return "{}*I".format(self.imag()._repr_())
        else:
            return "{} + {}*I".format(self.real()._repr_(),
                                        self.imag()._repr_())

    cpdef ComplexIntervalFieldElement _interval(self):
        """
        Return :class:`ComplexIntervalFieldElement` of the same value.

        OUTPUT:

        A :class:`ComplexIntervalFieldElement`.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
            sage: CBF = ComplexBallField() # optional - arb
            sage: a = CBF(CIF(2, 2))                # optional - arb
            sage: a._interval()   # optional - arb
            2 + 2*I
        """
        cdef ComplexIntervalFieldElement target = ComplexIntervalField(prec(self))(0)
        acb_to_ComplexIntervalFieldElement(target, self.value)
        return target

    # Comparisons and predicates

    def is_zero(self):
        """
        Return ``True`` iff the midpoint and radius of this ball are both zero.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
            sage: CBF = ComplexBallField() # optional - arb
            sage: CBF(0).is_zero() # optional - arb
            True
            sage: CBF(RIF(-0.5, 0.5)).is_zero() # optional - arb
            False
        """
        return acb_is_zero(self.value)

    def __nonzero__(self):
        """
        Return ``True`` iff zero is not contained in the interval represented
        by this ball.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
            sage: CBF = ComplexBallField() # optional - arb
            sage: bool(CBF(pi, 1/3)) # optional - arb
            True
            sage: bool(CBF(RIF(-0.5, 0.5), 1/3)) # optional - arb
            True
            sage: bool(CBF(1/3, RIF(-0.5, 0.5))) #optional - arb
            True
            sage: bool(CBF(RIF(-0.5, 0.5), RIF(-0.5, 0.5))) #optional - arb
            False
        """
        return acb_is_nonzero(self.value)

    def is_exact(self):
        """
        Return ``True`` iff the radius of this ball is zero.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
            sage: CBF = ComplexBallField() # optional - arb
            sage: CBF(1).is_exact() # optional - arb
            True
            sage: CBF(1/3, 1/3).is_exact() # optional - arb
            False
        """
        return acb_is_exact(self.value)

    def __richcmp__(left, right, int op):
        """
        Compare ``left`` and ``right``.

        For more information, see :mod:`sage.rings.complex_ball_acb`.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
            sage: CBF = ComplexBallField() # optional - arb
            sage: a = CBF(1) # optional - arb
            sage: b = CBF(1) # optional - arb
            sage: a is b # optional - arb
            False
            sage: a == b # optional - arb
            True
            sage: a = CBF(1/3) # optional - arb
            sage: a.is_exact() # optional - arb
            False
            sage: b = CBF(1/3) # optional - arb
            sage: b.is_exact() # optional - arb
            False
            sage: a == b # optional - arb
            False
        """
        return (<Element>left)._richcmp(right, op)

    cpdef _richcmp_(left, Element right, int op):
        """
        Compare ``left`` and ``right``.

        For more information, see :mod:`sage.rings.complex_ball_acb`.

        EXAMPLES::

            sage: from sage.rings.complex_ball_acb import ComplexBallField # optional - arb
            sage: CBF = ComplexBallField() # optional - arb
            sage: a = CBF(1, 2) # optional - arb
            sage: b = CBF(1, 2) # optional - arb
            sage: a is b # optional - arb
            False
            sage: a == b # optional - arb
            True

        TESTS:

            Balls whose intersection consists of one point::

                sage: a = CBF(RIF(1, 2), RIF(1, 2)) # optional - arb
                sage: b = CBF(RIF(2, 4), RIF(2, 4)) # optional - arb
                sage: a < b # optional - arb
                Traceback (most recent call last):
                ...
                TypeError: No order is defined for ComplexBalls.
                sage: a > b # optional - arb
                Traceback (most recent call last):
                ...
                TypeError: No order is defined for ComplexBalls.
                sage: a <= b # optional - arb
                Traceback (most recent call last):
                ...
                TypeError: No order is defined for ComplexBalls.
                sage: a >= b # optional - arb
                Traceback (most recent call last):
                ...
                TypeError: No order is defined for ComplexBalls.
                sage: a == b # optional - arb
                False
                sage: a != b # optional - arb
                False

            Balls with non-trivial intersection::

                sage: a = CBF(RIF(1, 4), RIF(1, 4)) # optional - arb
                sage: a = CBF(RIF(2, 5), RIF(2, 5)) # optional - arb
                sage: a == b # optional - arb
                False
                sage: a != b # optional - arb
                False

            One ball contained in another::

                sage: a = CBF(RIF(1, 4), RIF(1, 4)) # optional - arb
                sage: b = CBF(RIF(2, 3), RIF(2, 3)) # optional - arb
                sage: a == b # optional - arb
                False
                sage: a != b # optional - arb
                False

            Disjoint balls::

                sage: a = CBF(1/3, 1/3) # optional - arb
                sage: b = CBF(1/5, 1/5) # optional - arb
                sage: a == b # optional - arb
                False
                sage: a != b # optional - arb
                True

            Exact elements::

                sage: a = CBF(2, 2) # optional - arb
                sage: b = CBF(2, 2) # optional - arb
                sage: a.is_exact() # optional - arb
                True
                sage: b.is_exact() # optional - arb
                True
                sage: a == b # optional - arb
                True
                sage: a != b # optional - arb
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
