# sage_setup: distribution = sagemath-flint
# sage.doctest: # needs sage.symbolic
"""
Real intervals with a fixed absolute precision
"""

from sage.ext.stdsage cimport PY_NEW

from sage.libs.gmp.mpz cimport *

from sage.structure.factory import UniqueFactory
from sage.structure.element cimport Element, FieldElement
from sage.rings.ring cimport Field
from sage.rings.integer cimport Integer
import sage.rings.abc

from sage.structure.parent cimport Parent
from sage.structure.element cimport parent as parent_of

from sage.rings.real_mpfr import RR_min_prec
from sage.rings.real_mpfi import RealIntervalField, RealIntervalFieldElement
from sage.rings.rational_field import QQ

cdef Integer zero = Integer(0)
cdef Integer one = Integer(1)

cpdef inline Integer shift_floor(Integer x, long shift):
    r"""
    Return `x / 2^s` where `s` is the value of ``shift``, rounded towards
    `-\infty`. For internal use.

    EXAMPLES::

        sage: from sage.rings.real_interval_absolute import shift_floor
        sage: shift_floor(15, 2)
        3
        sage: shift_floor(-15, 2)
        -4
    """
    cdef Integer z = PY_NEW(Integer)
    mpz_fdiv_q_2exp(z.value, x.value, shift)
    return z

cpdef inline Integer shift_ceil(Integer x, long shift):
    r"""
    Return `x / 2^s` where `s` is the value of ``shift``, rounded towards
    `+\infty`. For internal use.

    EXAMPLES::

        sage: from sage.rings.real_interval_absolute import shift_ceil
        sage: shift_ceil(15, 2)
        4
        sage: shift_ceil(-15, 2)
        -3
        sage: shift_ceil(32, 2)
        8
        sage: shift_ceil(-32, 2)
        -8
    """
    cdef Integer z = PY_NEW(Integer)
    mpz_cdiv_q_2exp(z.value, x.value, shift)
    return z


class Factory(UniqueFactory):
    def create_key(self, prec):
        """
        The only piece of data is the precision.

        TESTS::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: RealIntervalAbsoluteField.create_key(1000)
            1000
        """
        return prec

    def create_object(self, version, prec):
        """
        Ensures uniqueness.

        TESTS::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: RealIntervalAbsoluteField(23) is RealIntervalAbsoluteField(23) # indirect doctest
            True
        """
        return RealIntervalAbsoluteField_class(prec)


RealIntervalAbsoluteField = Factory('sage.rings.real_interval_absolute.RealIntervalAbsoluteField')
RealIntervalAbsoluteField.__doc__ = RealIntervalAbsoluteField_class.__doc__


cdef class RealIntervalAbsoluteField_class(Field):
    """
    This field is similar to the :class:`RealIntervalField` except instead of
    truncating everything to a fixed relative precision, it maintains a
    fixed absolute precision.

    Note that unlike the standard real interval field, elements in this
    field can have different size and experience coefficient blowup. On
    the other hand, it avoids precision loss on addition and subtraction.
    This is useful for, e.g., series computations for special functions.

    EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(10); R
            Real Interval Field with absolute precision 2^-10
            sage: R(3/10)
            0.300?
            sage: R(1000003/10)
            100000.300?
            sage: R(1e100) + R(1) - R(1e100)
            1
    """

    cdef long _absprec

    def __init__(self, absprec):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: RealIntervalAbsoluteField(100)
            Real Interval Field with absolute precision 2^-100
            sage: RealIntervalAbsoluteField(-100)
            Traceback (most recent call last):
              File "<ipython console>", line 1, in <module>
              File "real_interval_absolute.pyx", line 81, in sage.rings.real_interval_absolute.RealIntervalAbsoluteField.__init__ (sage/rings/real_interval_absolute.c:2463)
            ValueError: Absolute precision must be positive.
        """
        if absprec < 0:
            raise ValueError("Absolute precision must be positive.")
        Field.__init__(self, self)
        self._absprec = absprec

    def __reduce__(self):
        """
        Used for pickling.

        TESTS::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(100)
            sage: loads(dumps(R))
            Real Interval Field with absolute precision 2^-100
        """
        return RealIntervalAbsoluteField, (self._absprec,)

    def _element_constructor_(self, x):
        """
        Construct an element with ``self`` as the parent.

        TESTS::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(100)
            sage: R(1/2) # indirect doctest
            0.50000000000000000000000000000000?
        """
        return RealIntervalAbsoluteElement(self, x)

    cpdef _coerce_map_from_(self, R):
        """
        Anything that coerces into the reals coerces into this ring.

        TESTS::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(100)
            sage: R.has_coerce_map_from(RR) # indirect doctest
            True
            sage: R.has_coerce_map_from(QQ)
            True
            sage: R.has_coerce_map_from(Qp(5))                                          # needs sage.rings.padics
            False

            sage: R(1/2) + 100
            100.5000000000000000000000000000000?
            sage: R(1/2) + 1.75
            2.2500000000000000000000000000000?

            sage: R10 = RealIntervalAbsoluteField(10)
            sage: R10(1/4) + R(1/4)
            0.50000?
        """
        if isinstance(R, RealIntervalAbsoluteField_class):
            return self._absprec < (<RealIntervalAbsoluteField_class>R)._absprec
        elif isinstance(R, sage.rings.abc.RealIntervalField):
            return True
        else:
            return RR_min_prec.has_coerce_map_from(R)

    def _repr_(self):
        """
        Return the string representation of ``self``.

        TESTS::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(100)
            sage: print(R)
            Real Interval Field with absolute precision 2^-100
            sage: R._repr_()
            'Real Interval Field with absolute precision 2^-100'
        """
        return "Real Interval Field with absolute precision 2^-%s" % self._absprec

    def absprec(self):
        """
        Return the absolute precision of ``self``.

        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(100)
            sage: R.absprec()
            100
            sage: RealIntervalAbsoluteField(5).absprec()
            5
        """
        return self._absprec


cdef inline shift_left(value, shift):
    """
    Utility function for operands that don't support the ``<<`` operator.
    """
    try:
        return value << shift
    except TypeError:
        if isinstance(value, (str, list, tuple)):
            # Better than the OverflowError we would get from trying to multiply.
            raise
        else:
            return value * (one << shift)

cdef class RealIntervalAbsoluteElement(FieldElement):

    # This could be optimized by letting these be raw mpz_t.
    cdef Integer _mantissa  # left endpoint
    cdef Integer _diameter

    def __init__(self, RealIntervalAbsoluteField_class parent, value):
        """
        Create a :class:`RealIntervalAbsoluteElement`.

        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(50)
            sage: R(1)
            1
            sage: R(1/3)
            0.333333333333334?
            sage: R(1.3)
            1.300000000000000?
            sage: R(pi)
            3.141592653589794?
            sage: R((11, 12))
            12.?
            sage: R((11, 11.00001))
            11.00001?

            sage: R100 = RealIntervalAbsoluteField(100)
            sage: R(R100((5,6)))
            6.?
            sage: R100(R((5,6)))
            6.?
            sage: RIF(CIF(NaN))
            [.. NaN ..]
        """
        Element.__init__(self, parent)

        if isinstance(value, RealIntervalAbsoluteElement):
            prec_diff = (<RealIntervalAbsoluteField_class>(<Element>value)._parent)._absprec - parent._absprec
            if prec_diff > 0:
                self._mantissa = shift_floor((<RealIntervalAbsoluteElement>value)._mantissa, prec_diff)
                self._diameter = shift_ceil((<RealIntervalAbsoluteElement>value)._diameter, prec_diff)
            else:
                self._mantissa = (<RealIntervalAbsoluteElement>value)._mantissa << -prec_diff
                self._diameter = (<RealIntervalAbsoluteElement>value)._diameter << -prec_diff
            return

        if isinstance(value, tuple):
            value, upper = value
        elif isinstance(value, RealIntervalFieldElement):
            value, upper = value.lower(), value.upper()
        else:
            upper = None
        value = shift_left(value, parent._absprec)
        if upper is not None:
            upper = shift_left(upper, parent._absprec)
        else:
            upper = value
        from sage.functions.other import floor, ceil
        try:
            self._mantissa = floor(value)
            self._diameter = ceil(upper) - self._mantissa
        except OverflowError:
            raise TypeError(type(value))

    def __reduce__(self):
        """
        Used for pickling.

        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(50)
            sage: loads(dumps(R(1/16)))
            0.06250000000000000?
            sage: R = RealIntervalAbsoluteField(100)
            sage: loads(dumps(R(1/3)))
            0.333333333333333333333333333334?
            sage: loads(dumps(R(pi))).endpoints() == R(pi).endpoints()
            True
        """
        return RealIntervalAbsoluteElement, (self._parent, self.endpoints())

    cdef _new_c(self, Integer _mantissa, Integer _diameter):
        cdef RealIntervalAbsoluteElement x
        x = <RealIntervalAbsoluteElement>RealIntervalAbsoluteElement.__new__(RealIntervalAbsoluteElement)
        x._parent = self._parent
        x._mantissa = _mantissa
        x._diameter = _diameter
        return x

    cpdef lower(self):
        """
        Return the lower bound of ``self``.

        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(50)
            sage: R(1/4).lower()
            1/4
        """
        return QQ(self._mantissa) >> (<RealIntervalAbsoluteField_class>self._parent)._absprec

    cpdef midpoint(self):
        """
        Return the midpoint of ``self``.

        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(100)
            sage: R(1/4).midpoint()
            1/4
            sage: R(pi).midpoint()
            7964883625991394727376702227905/2535301200456458802993406410752
            sage: R(pi).midpoint().n()
            3.14159265358979
        """
        return (self._mantissa + self._diameter / 2) >> (<RealIntervalAbsoluteField_class>self._parent)._absprec

    cpdef upper(self):
        """
        Return the upper bound of ``self``.

        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(50)
            sage: R(1/4).upper()
            1/4
        """
        return QQ(self._mantissa + self._diameter) >> (<RealIntervalAbsoluteField_class>self._parent)._absprec

    cpdef absolute_diameter(self):
        """
        Return the diameter ``self``.

        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(10)
            sage: R(1/4).absolute_diameter()
            0
            sage: a = R(pi)
            sage: a.absolute_diameter()
            1/1024
            sage: a.upper() - a.lower()
            1/1024
        """
        return QQ(self._diameter) >> (<RealIntervalAbsoluteField_class>self._parent)._absprec

    diameter = absolute_diameter

    cpdef endpoints(self):
        """
        Return the left and right endpoints of ``self``, as a tuple.

        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(10)
            sage: R(1/4).endpoints()
            (1/4, 1/4)
            sage: R((1,2)).endpoints()
            (1, 2)
        """
        return self.lower(), self.upper()

    def _real_mpfi_(self, R):
        """
        Create a (relative) real interval out of this absolute real interval.

        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(10)
            sage: R(1/2)._real_mpfi_(RIF)
            0.50000000000000000?

            sage: a = RealIntervalAbsoluteField(100)(1/3)
            sage: RIF(a)
            0.3333333333333334?
        """
        return R(self._mantissa, self._mantissa + self._diameter) >> (<RealIntervalAbsoluteField_class>self._parent)._absprec

    cpdef long mpfi_prec(self) noexcept:
        """
        Return the precision needed to represent this value as an mpfi interval.

        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(10)
            sage: R(10).mpfi_prec()
            14
            sage: R(1000).mpfi_prec()
            20
        """
        return max(mpz_sizeinbase(self._mantissa.value, 2), mpz_sizeinbase(self._diameter.value, 2))

    def _repr_(self):
        """
        Leverage real interval printing.

        TESTS::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(10)
            sage: R(1/3) # indirect doctest
            0.334?
            sage: R(10^50/3)
            3.3333333333333333333333333333333333333333333333333334?e49
            sage: R(0)
            0
        """
        prec = max(self.mpfi_prec(), 5)
        return repr(self._real_mpfi_(RealIntervalField(prec)))

    def __hash__(self):
        """
        Hash to the midpoint of the interval.

        TESTS::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(10)
            sage: hash(R(10))
            10
            sage: hash(R((11,13)))
            12
            sage: hash(R(1/4)) == hash(1/4)
            True
            sage: hash(R(pi))
            891658780           # 32-bit
            532995478001132060  # 64-bit
        """
        return hash(self.midpoint())

    def __contains__(self, x):
        """
        Return whether the given value lies in this interval.

        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(50)
            sage: 1 in R((1,2))
            True
            sage: 2 in R((1,2))
            True
            sage: 3 in R((1,2))
            False
            sage: 1.75 in R((1,2))
            True
        """
        x *= (one << self._parent.absprec())
        return self._mantissa <= x <= self._mantissa + self._diameter

    cpdef bint is_positive(self) noexcept:
        """
        Return whether ``self`` is definitely positive.

        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(10)
            sage: R(10).is_positive()
            True
            sage: R((10,11)).is_positive()
            True
            sage: R((0,11)).is_positive()
            False
            sage: R((-10,11)).is_positive()
            False
            sage: R((-10,-1)).is_positive()
            False
            sage: R(pi).is_positive()
            True
        """
        return mpz_sgn(self._mantissa.value) == 1

    cpdef bint contains_zero(self) noexcept:
        """
        Return whether ``self`` contains zero.

        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(10)
            sage: R(10).contains_zero()
            False
            sage: R((10,11)).contains_zero()
            False
            sage: R((0,11)).contains_zero()
            True
            sage: R((-10,11)).contains_zero()
            True
            sage: R((-10,-1)).contains_zero()
            False
            sage: R((-10,0)).contains_zero()
            True
            sage: R(pi).contains_zero()
            False
        """
        return (mpz_sgn(self._mantissa.value) == 0
                or (mpz_sgn(self._mantissa.value) == -1 and mpz_cmpabs(self._mantissa.value, self._diameter.value) <= 0))

    cpdef bint is_negative(self) noexcept:
        """
        Return whether ``self`` is definitely negative.

        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(100)
            sage: R(10).is_negative()
            False
            sage: R((10,11)).is_negative()
            False
            sage: R((0,11)).is_negative()
            False
            sage: R((-10,11)).is_negative()
            False
            sage: R((-10,-1)).is_negative()
            True
            sage: R(pi).is_negative()
            False
        """
        return (mpz_sgn(self._mantissa.value) == -1
                and mpz_cmpabs(self._mantissa.value, self._diameter.value) > 0)

    cdef bint is_exact(self) noexcept:
        return not self._diameter

    def __bool__(self):
        """
        Return ``True`` for anything except exact zero.

        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(10)
            sage: bool(R(1))
            True
            sage: bool(R(0))
            False
            sage: bool(R((0,1)))
            True
        """
        return not not self._mantissa or not not self._diameter

    def __neg__(self):
        """
        TESTS::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(100)
            sage: -R(1/2)
            -0.50000000000000000000000000000000?
            sage: -R((101,102))
            -102.?
        """
        return self._new_c(-self._mantissa - self._diameter, self._diameter)

    def __abs__(self):
        """
        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(100)
            sage: abs(-R(1/4))
            0.2500000000000000000000000000000?
        """
        return self.abs()

    cpdef abs(self):
        """
        Return the absolute value of ``self``.

        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(100)
            sage: R(1/3).abs()
            0.333333333333333333333333333334?
            sage: R(-1/3).abs()
            0.333333333333333333333333333334?
            sage: R((-1/3, 1/2)).abs()
            1.?
            sage: R((-1/3, 1/2)).abs().endpoints()
            (0, 1/2)
            sage: R((-3/2, 1/2)).abs().endpoints()
            (0, 3/2)
        """
        if self.is_positive():
            return self
        elif self.is_negative():
            return -self
        else:
            return self._new_c(zero, max(-self._mantissa, self._mantissa + self._diameter))

    cpdef _add_(self, _other):
        """
        TESTS::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(10)
            sage: R(1) + R(2) # indirect doctest
            3
            sage: R(1e100) + R(0.1) + R(-1e100)
            0.100?
            sage: (R((1,2)) + 1).endpoints()
            (2, 3)
            sage: (R((1,2)) + R((-10,0))).endpoints()
            (-9, 2)
        """
        cdef RealIntervalAbsoluteElement other = <RealIntervalAbsoluteElement>_other
        return self._new_c(self._mantissa + other._mantissa, self._diameter + other._diameter)

    cpdef _sub_(self, _other):
        """
        TESTS::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(10)
            sage: R(1) - R(2) # indirect doctest
            -1
            sage: R(1e100) - R(0.1) - R(1e100)
            -0.100?
            sage: (R((1,2)) - 1).endpoints()
            (0, 1)
            sage: (R((1,2)) - R((10,100))).endpoints()
            (-99, -8)
            sage: R(pi) - R(pi)
            0.000?
        """
        cdef RealIntervalAbsoluteElement other = <RealIntervalAbsoluteElement>_other
        return self._new_c(self._mantissa - other._mantissa - other._diameter, self._diameter + other._diameter)

    cpdef _mul_(self, _other):
        """
        TESTS::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(10)
            sage: R(2) * R(3) # indirect doctest
            6
            sage: R(2) * R(-3)
            -6
            sage: elts = [R((left, right)) for left in [-2..2] for right in [left+1..2]]
            sage: elts
            [-2.?, -1.?, 0.?e1, 0.?e1, -1.?, 0.?, 0.?e1, 1.?, 1.?, 2.?]
            sage: for a in elts:
            ....:     for b in elts:
            ....:         if (a*b).lower() != (a._real_mpfi_(RIF)*b._real_mpfi_(RIF)).lower():
            ....:             print(a, b)
            ....:         if (a*b).upper() != (a._real_mpfi_(RIF)*b._real_mpfi_(RIF)).upper():
            ....:             print(a, b)
            sage: R(pi) * R(pi) - R(pi^2)
            0.00?
        """
        cdef bint negate = False
        cdef RealIntervalAbsoluteElement other = <RealIntervalAbsoluteElement>_other
        cdef long absprec = (<RealIntervalAbsoluteField_class>self._parent)._absprec

        # Break symmetry.
        if self.is_negative():
            negate = not negate
            self = -self
        if other.is_negative():
            negate = not negate
            other = -other
        elif other.contains_zero():
            self, other = other, self

        if self.is_positive():
            if other.is_positive():
                res = self._new_c(shift_floor(self._mantissa * other._mantissa, absprec),
                                  shift_ceil(self._mantissa * other._diameter + (other._mantissa + other._diameter) * self._diameter, absprec))
            else:
                res = self._new_c(shift_floor((self._mantissa + self._diameter) * other._mantissa, absprec),
                                  shift_ceil((self._mantissa + self._diameter) * other._diameter, absprec))
        else:
            # They both contain zero.
            self_lower, self_upper = self._mantissa, self._mantissa + self._diameter
            other_lower, other_upper = other._mantissa, other._mantissa + other._diameter
            lower = min(self_lower * other_upper, self_upper * other_lower)
            upper = max(self_lower * other_lower, self_upper * other_upper)
            res = self._new_c(shift_floor(lower, absprec),
                              shift_ceil(upper - lower, absprec))
        if negate:
            res = -res
        return res

    cpdef _acted_upon_(self, x, bint self_on_left):
        """
        ``Absprec * relprec -> absprec`` works better than coercing both
        operands to absolute precision first.

        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(4)
            sage: 1/4 * R(100)
            25
            sage: 1/100 * R(100)
            1
            sage: R(1/100) * R(100)
            1.?e1
            sage: RIF(1/100) * R(100)
            1.0?

            sage: R(1.5)._acted_upon_(3, True)
            4.500?
        """
        if x < 0:
            neg = self._acted_upon_(-x, self_on_left)
            return None if neg is None else -neg
        if type(x) in (int, Integer):
            return self._new_c(self._mantissa * x, self._diameter * x)
        P = parent_of(x)
        if isinstance(P, Parent):
            if P.is_exact():
                left = (self._mantissa * x).floor()
                right = ((self._mantissa + self._diameter) * x).ceil()
                return self._new_c(left, right - left)
            elif isinstance(x, RealIntervalFieldElement):
                if x.contains_zero() or self.contains_zero():
                    return self * RealIntervalAbsoluteElement(self._parent, (x.lower(), x.upper()))
                # Remember, x > 0
                if self.is_positive():
                    left = (self._mantissa * x.lower()).floor()
                    right = ((self._mantissa + self._diameter) * x.upper()).ceil()
                else:
                    left = (self._mantissa * x.upper()).floor()
                    right = ((self._mantissa + self._diameter) * x.lower()).ceil()
                return self._new_c(left, right - left)

    def __invert__(self):
        """
        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(10)
            sage: ~R(2)
            0.50000?
            sage: ~~R(2)
            2
            sage: ~R(3)
            0.334?
            sage: ~~R(3)
            3.00?

            sage: R = RealIntervalAbsoluteField(200)
            sage: ~R(1e10)
            1.00000000000000000000000000000000000000000000000000?e-10
            sage: ~~R(1e10)
            1.00000000000000000000000000000000000000000000000000?e10
            sage: (~R((1,2))).endpoints()
            (1/2, 1)
            sage: (~R((1/4,8))).endpoints()
            (1/8, 4)
            sage: R(1/pi) - 1/R(pi)
            0.?e-60
        """
        if self.contains_zero():
            raise ZeroDivisionError("Inversion of an interval containing zero.")
        cdef long absprec = (<RealIntervalAbsoluteField_class>self._parent)._absprec
        cdef bint negate
        if self.is_negative():
            self = -self
            negate = True
        else:
            negate = False

        # Let our (positive) interval be [2^-B a, 2^-B b].
        # Then its inverse is [2^-B 2^(2B)/b , 2^-B 2^(2B)/a].

        cdef Integer mantissa = <Integer>PY_NEW(Integer)
        cdef Integer diameter = <Integer>PY_NEW(Integer)
        cdef mpz_t scaling_factor
        mpz_init_set_ui(scaling_factor, 1)
        try:
            mpz_set_ui(scaling_factor, 1)
            mpz_mul_2exp(scaling_factor, scaling_factor, 2 * absprec)
            # Use diameter as temp value for right endpoint...
            mpz_add(diameter.value, self._mantissa.value, self._diameter.value)
            mpz_fdiv_q(mantissa.value, scaling_factor, diameter.value)
            # Divide a second time to get the rounding correct...
            mpz_cdiv_q(diameter.value, scaling_factor, self._mantissa.value)
            mpz_sub(diameter.value, diameter.value, mantissa.value)
        finally:
            mpz_clear(scaling_factor)

        res = self._new_c(mantissa, diameter)
        if negate:
            res = -res
        return res

    cpdef _div_(self, _other):
        """
        TESTS::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(10)
            sage: R(5)/R(2) # indirect doctest
            2.5000?
            sage: a = R((5,6))/R((2,4))
            sage: a.endpoints()
            (5/4, 3)
            sage: a = R((-5,6))/R((2,4))
            sage: a.endpoints()
            (-5/2, 3)
            sage: R(1e100) / R(2e100)
            0.50000?
        """
        cdef RealIntervalAbsoluteElement other = <RealIntervalAbsoluteElement>_other
        if other.contains_zero():
            raise ZeroDivisionError("Division by an interval containing zero.")

        cdef Integer mantissa = <Integer>PY_NEW(Integer)
        cdef Integer diameter = <Integer>PY_NEW(Integer)
        cdef long absprec = (<RealIntervalAbsoluteField_class>self._parent)._absprec
        cdef bint negate = False
        cdef mpz_t temp
        mpz_init(temp)

        try:

            if self.is_negative():
                negate = True
                self = -self

            mpz_mul_2exp(mantissa.value, self._mantissa.value, absprec)
            if self.contains_zero():
                mpz_fdiv_q(mantissa.value, mantissa.value, other._mantissa.value)
            else:
                mpz_add(temp, other._mantissa.value, other._diameter.value)
                mpz_fdiv_q(mantissa.value, mantissa.value, temp)

            mpz_add(temp, self._mantissa.value, self._diameter.value)
            mpz_mul_2exp(diameter.value, temp, absprec)
            mpz_cdiv_q(diameter.value, diameter.value, other._mantissa.value)
            mpz_sub(diameter.value, diameter.value, mantissa.value)

        finally:
            mpz_clear(temp)

        res = self._new_c(mantissa, diameter)
        if negate:
            res = -res
        return res

    def __lshift__(RealIntervalAbsoluteElement self, long n):
        """
        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(10)
            sage: R(1) << 2
            4
            sage: R(3) << -2
            0.75000?
            sage: (R((1/2, 5)) << 10).endpoints()
            (512, 5120)
        """
        return self.shift(n)

    def __rshift__(RealIntervalAbsoluteElement self, long n):
        """
        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(10)
            sage: R(1) >> 2
            0.2500?
            sage: R(3) >> -2
            12
            sage: (R((1/2, 5)) >> 10).endpoints()
            (0, 5/1024)
        """
        return self.shift(-n)

    cdef shift(self, long n):
        if n >= 0:
            return self._new_c(self._mantissa << n, self._diameter << n)
        else:
            return self._new_c(shift_floor(self._mantissa, -n), shift_ceil(self._diameter, -n))

    def __pow__(self, exponent, dummy):
        """
        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(10)
            sage: R(10)^10
            10000000000
            sage: R(10)^100
            10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
            sage: R(10)^-2
            0.010?
            sage: R(10)^-10
            0.001?
            sage: R(100)^(1/2)
            10.00?
            sage: (R((2,3))^2).endpoints()
            (4, 9)
            sage: (R((-5,3))^3).endpoints()
            (-125, 75)
            sage: (R((-5,3))^4).endpoints()
            (0, 625)
        """
        cdef RealIntervalAbsoluteElement base
        try:
            base = self
        except TypeError:
            base = exponent.parent()(self)
        cdef RealIntervalAbsoluteField_class parent = self.parent()
        if not base:
            if not exponent:
                raise ZeroDivisionError
            else:
                return base
        import math
        cdef double height
        base_height = max(abs(base.lower()), abs(base.upper()))
        height = math.log(base_height) / math.log(2)
        height *= <double>(exponent.midpoint() if isinstance(exponent, RealIntervalAbsoluteElement) else exponent)
        relprec = max(<long>height + parent._absprec, 10)
        RIF = RealIntervalField(relprec)
        if isinstance(exponent, RealIntervalAbsoluteElement):
            exponent = exponent._real_mpfi_(RIF)
        return parent(base._real_mpfi_(RIF)**exponent)

    def __getattr__(self, name):
        """
        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(100)
            sage: R(1).sin()
            0.841470984807896506652502321631?
            sage: R(2).log()
            0.693147180559945309417232121458?
            sage: R(1).exp().log()
            1.00000000000000000000000000000?

            sage: R((0,10)).sin().endpoints()
            (-1, 1)
            sage: R(1).sin().parent()
            Real Interval Field with absolute precision 2^-100
        """
        if name[0] != '_' and hasattr(RealIntervalFieldElement, name):
            return MpfrOp(self, name)
        else:
            raise AttributeError(name)

    def sqrt(self):
        """
        Return the square root of ``self``.

        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
            sage: R = RealIntervalAbsoluteField(100)
            sage: R(2).sqrt()
            1.414213562373095048801688724210?
            sage: R((4,9)).sqrt().endpoints()
            (2, 3)
        """
        return self._parent(self._real_mpfi_(RealIntervalField(self.mpfi_prec())).sqrt())

cdef class MpfrOp:
    """
    This class is used to endow absolute real interval field elements with
    all the methods of (relative) real interval field elements.

    EXAMPLES::

        sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField
        sage: R = RealIntervalAbsoluteField(100)
        sage: R(1).sin()
        0.841470984807896506652502321631?
    """
    cdef object name
    cdef RealIntervalAbsoluteElement value

    def __init__(self, value, name):
        """
        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField, MpfrOp
            sage: R = RealIntervalAbsoluteField(100)
            sage: MpfrOp(R(1), 'tan')()
            1.557407724654902230506974807459?
        """
        self.name = name
        self.value = value

    def __call__(self, *args):
        """
        EXAMPLES::

            sage: from sage.rings.real_interval_absolute import RealIntervalAbsoluteField, MpfrOp
            sage: R = RealIntervalAbsoluteField(100)
            sage: curried_log = MpfrOp(R(2), 'log')
            sage: curried_log()
            0.693147180559945309417232121458?
            sage: curried_log(2)
            1
        """
        cdef RealIntervalAbsoluteField_class parent = self.value._parent
        cdef long absprec = parent._absprec
        cdef long relprec = self.value.mpfi_prec()
        for a in args:
            if isinstance(a, RealIntervalAbsoluteElement):
                absprec = min(absprec, (<RealIntervalAbsoluteField_class>(<RealIntervalAbsoluteElement>a)._parent)._absprec)
                relprec = max(absprec, (<RealIntervalAbsoluteElement>a).mpfi_prec())
        if parent._absprec > absprec:
            parent = RealIntervalAbsoluteField(absprec)
        if relprec < 53:
            relprec = 53
        R = RealIntervalField(relprec)
        new_args = []
        for a in args:
            if isinstance(a, RealIntervalAbsoluteElement):
                new_args.append(a._real_mpfi_(R))
            else:
                new_args.append(a)
        return parent(getattr(self.value._real_mpfi_(R), self.name)(*new_args))
