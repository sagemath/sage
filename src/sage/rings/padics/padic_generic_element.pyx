"""
`p`-adic Generic Element

Elements of `p`-adic Rings and Fields

AUTHORS:

- David Roe

- Genya Zaytman: documentation

- David Harvey: doctests

- Julian Rueth: fixes for exp() and log(), implemented gcd, xgcd
"""
# ****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed@math.harvard.edu>
#                     2007      William Stein <wstein@gmail.com>
#                     2013-2014 Julian Rueth <julian.rueth@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.ext.stdsage cimport PY_NEW

cimport sage.rings.padics.local_generic_element
from sage.libs.gmp.mpz cimport mpz_set_si
from sage.arith.srange import srange
from sage.rings.padics.local_generic_element cimport LocalGenericElement
from sage.rings.padics.precision_error import PrecisionError
from sage.rings.rational cimport Rational
from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import infinity
from sage.structure.element import coerce_binop
from sage.structure.richcmp cimport rich_to_bool

cdef long maxordp = (1L << (sizeof(long) * 8 - 2)) - 1

cdef class pAdicGenericElement(LocalGenericElement):
    cpdef _richcmp_(left, right, int op):
        r"""
        First compare valuations, then compare normalized
        residue of unit part.

        EXAMPLES::

            sage: R = Zp(19, 5,'capped-rel','series'); K = Qp(19, 5, 'capped-rel', 'series')
            sage: a = R(2); a
            2 + O(19^5)
            sage: b = R(3); b
            3 + O(19^5)
            sage: a < b
            True
            sage: a = K(2); a
            2 + O(19^5)
            sage: b = K(3); b
            3 + O(19^5)
            sage: a < b
            True

        ::

            sage: R = Zp(5); a = R(5, 6); b = R(5 + 5^6, 8)
            sage: a == b  # indirect doctest
            True

        ::

            sage: R = Zp(5)
            sage: a = R(17)
            sage: b = R(21)
            sage: a == b
            False
            sage: a < b
            True

        ::

            sage: R = ZpCA(5)
            sage: a = R(17)
            sage: b = R(21)
            sage: a == b
            False
            sage: a < b
            True

        ::

            sage: R = ZpFM(5)
            sage: a = R(17)
            sage: b = R(21)
            sage: a == b
            False
            sage: a < b
            True
        """
        # handle exact zeros
        left_zero = left._is_exact_zero()
        right_zero = right._is_exact_zero()
        if left_zero and right_zero:
            return rich_to_bool(op, 0)
        elif left_zero:
            if right.is_zero():
                return rich_to_bool(op, 0)
            return rich_to_bool(op, 1)
        elif right_zero:
            if left.is_zero():
                return rich_to_bool(op, 0)
            return rich_to_bool(op, -1)

        m = min(left.precision_absolute(), right.precision_absolute())
        x_ordp = left.valuation()
        left_zero =  bool(x_ordp >= m)
        y_ordp = right.valuation()
        right_zero = bool(y_ordp >= m)
        # handle approximate zeros
        if left_zero and right_zero:
            return rich_to_bool(op, 0)
        elif left_zero:
            return rich_to_bool(op, 1)
        elif right_zero:
            return rich_to_bool(op, -1)

        # no more zeros
        if x_ordp < y_ordp:
            return rich_to_bool(op, -1)
        elif x_ordp > y_ordp:
            return rich_to_bool(op, 1)
        else:  # equal ordp
            test = (<pAdicGenericElement>left.unit_part())._cmp_units(right.unit_part())
            return rich_to_bool(op, test)

    cdef int _cmp_units(left, pAdicGenericElement right) except -2:
        raise NotImplementedError

    cdef int _set_from_Integer(self, Integer x, absprec, relprec) except -1:
        raise NotImplementedError
    cdef int _set_from_mpz(self, mpz_t x) except -1:
        raise NotImplementedError
    cdef int _set_from_mpz_rel(self, mpz_t x, long relprec) except -1:
        raise NotImplementedError
    cdef int _set_from_mpz_abs(self, mpz_t value, long absprec) except -1:
        raise NotImplementedError
    cdef int _set_from_mpz_both(self, mpz_t x, long absprec, long relprec) except -1:
        raise NotImplementedError

    cdef int _set_from_Rational(self, Rational x, absprec, relprec) except -1:
        raise NotImplementedError
    cdef int _set_from_mpq(self, mpq_t x) except -1:
        raise NotImplementedError
    cdef int _set_from_mpq_rel(self, mpq_t x, long relprec) except -1:
        raise NotImplementedError
    cdef int _set_from_mpq_abs(self, mpq_t value, long absprec) except -1:
        raise NotImplementedError
    cdef int _set_from_mpq_both(self, mpq_t x, long absprec, long relprec) except -1:
        raise NotImplementedError

    cdef int _pshift_self(self, long shift) except -1:
        raise NotImplementedError

    cdef int _set_inexact_zero(self, long absprec) except -1:
        raise NotImplementedError
    cdef int _set_exact_zero(self) except -1:
        raise TypeError("this type of p-adic does not support exact zeros")

    cpdef bint _is_exact_zero(self) except -1:
        """
        Return ``True`` if ``self`` is exactly zero.  Since
        non-capped-relative elements cannot be exact, this function
        always returns False.

        EXAMPLES::

            sage: ZpCA(17)(0)._is_exact_zero()
            False
        """
        return False

    cpdef bint _is_inexact_zero(self) except -1:
        """
        Return ``True`` if ``self`` is indistinguishable from zero, but not
        exactly zero.

        EXAMPLES::

            sage: Zp(5)(0,5)._is_inexact_zero()
            True
        """
        raise NotImplementedError

    cpdef bint _is_zero_rep(self) except -1:
        """
        Return ``True`` if ``self`` is indistinguishable from zero.

        EXAMPLES::

            sage: ZpCA(17)(0,15)._is_zero_rep()
            True
        """
        return self._is_inexact_zero() or self._is_exact_zero()

    cdef bint _set_prec_abs(self, long absprec) except -1:
        self._set_prec_both(absprec, (<PowComputer_class>self.parent().prime_pow).prec_cap)

    cdef bint _set_prec_rel(self, long relprec) except -1:
        self._set_prec_both((<PowComputer_class>self.parent().prime_pow).prec_cap, relprec)

    cdef bint _set_prec_both(self, long absprec, long relprec) except -1:
        return 0

    def _quo_rem(self, right):
        """
        Quotient with remainder.

        We choose the remainder to have the same `p`-adic expansion
        as the numerator, but truncated at the valuation of the denominator.

        EXAMPLES::

            sage: R = Zp(3, 5)
            sage: R(12).quo_rem(R(2))
            (2*3 + O(3^6), 0)
            sage: R(2).quo_rem(R(12))
            (O(3^4), 2 + O(3^5))
            sage: q, r = R(4).quo_rem(R(12)); q, r
            (1 + 2*3 + 2*3^3 + O(3^4), 1 + O(3^5))
            sage: 12*q + r == 4
            True

        In general, the remainder is returned with maximal precision.
        However, it is not the case when the valuation of the divisor
        is greater than the absolute precision on the numerator::

            sage: R(1,2).quo_rem(R(81))
            (O(3^0), 1 + O(3^2))

        For fields the normal quotient always has remainder 0:

            sage: K = Qp(3, 5)
            sage: K(12).quo_rem(K(2))
            (2*3 + O(3^6), 0)
            sage: q, r = K(4).quo_rem(K(12)); q, r
            (3^-1 + O(3^4), 0)
            sage: 12*q + r == 4
            True

        You can get the same behavior for fields as for rings
        by using this underscored method::

            sage: K(12)._quo_rem(K(2))
            (2*3 + O(3^6), 0)
            sage: K(2)._quo_rem(K(12))
            (O(3^4), 2 + O(3^5))
        """
        if right._is_exact_zero():
            raise ZeroDivisionError("cannot divide by zero")
        if right.is_zero():
            raise PrecisionError("cannot divide by something indistinguishable from zero")
        K = self.parent()
        R = K.integer_ring()
        sval = self.valuation()
        diff = sval - right.valuation()
        srelprec = self.precision_relative()
        if diff + srelprec < 0:
            return K(0,0), self
        if srelprec == 0:
            return K(0,diff), K(0)
        if diff >= 0:  # remainder is 0
            return K(self/right), K(0)
        unit = R(self.unit_part())
        high = (unit << diff) >> (diff - sval)
        return K(high/right), K(self-high).lift_to_precision()

    def __floordiv__(self, right):
        """
        Divide ``self`` by ``right`` and throws away the nonintegral part if
        ``self.parent()`` is not a field.

        There are a number of reasonable definitions for floor
        division.  Any definition should satisfy the following
        identity:

        (1) a = (a // b) * b + a % b

        If a and b lie in a field, then setting a % b = 0 and a // b =
        a / b provides a canonical way of satisfying this equation.

        However, for elements of integer rings, there are many choices
        of definitions for a // b and a % b that satisfy this
        equation.  Since `p`-adic rings in Sage come equipped with a
        uniformizer pi, we can use the choice of uniformizer in our
        definitions.  Here are some other criteria we might ask for:

        (2) If b = pi^k, the series expansion (in terms of pi) of a //
        b is just the series expansion of a, shifted over by k terms.

        (2') The series expansion of a % pi^k has no terms above
        pi^(k-1).

        The conditions (2) and (2') are equivalent.  But when we
        generalize these conditions to arbitrary b they diverge.

        (3) For general b, the series expansion of a // b is just the
        series expansion of a / b, truncating terms with negative
        exponents of pi.

        (4) For general b, the series expansion of a % b has no terms
        above b.valuation() - 1.

        In order to satisfy (3), one defines

        a // b = (a / b.unit_part()) >> b.valuation()
        a % b = a - (a // b) * b

        In order to satisfy (4), one defines

        a % b = a.lift() % pi.lift()^b.valuation()
        a // b = ((a - a % b) >> b.valuation()) / b.unit_part()


        In Sage we choose option (4) since it has better precision behavior.

        EXAMPLES::

            sage: R = ZpCA(5); a = R(129378); b = R(2398125)
            sage: a // b  # indirect doctest
            1 + 2*5 + 2*5^3 + 4*5^4 + 5^6 + 5^7 + 5^8 + 4*5^9 + 2*5^10 + 4*5^11 + 4*5^12 + 2*5^13 + 3*5^14 + O(5^16)
            sage: a / b
            4*5^-4 + 3*5^-3 + 2*5^-2 + 5^-1 + 3 + 3*5 + 4*5^2 + 2*5^4 + 2*5^6 + 4*5^7 + 5^9 + 5^10 + 5^11 + O(5^12)
            sage: a % b
            3 + O(5^20)
            sage: a
            3 + 2*5^4 + 5^5 + 3*5^6 + 5^7 + O(5^20)
            sage: (a // b) * b + a % b
            3 + 2*5^4 + 5^5 + 3*5^6 + 5^7 + O(5^20)

        The alternative definition::

            sage: c = (a // b.unit_part()) >> b.valuation(); c
            3 + 3*5 + 4*5^2 + 2*5^4 + 2*5^6 + 4*5^7 + 5^9 + 5^10 + 5^11 + O(5^12)
            sage: othermod = a - c*b; othermod
            3 + 5^4 + 3*5^5 + 2*5^6 + 4*5^7 + 5^8 + O(5^16)
        """
        right = self.parent()(right)
        if right._is_inexact_zero():
            raise PrecisionError("cannot divide by something indistinguishable from zero")
        elif right._is_exact_zero():
            raise ZeroDivisionError("cannot divide by zero")
        return self._floordiv_(right)

    cpdef _floordiv_(self, right):
        """
        Implement floor division.

        EXAMPLES::

            sage: R = Zp(5, 5); a = R(77)
            sage: a // 15  # indirect doctest
            5 + O(5^4)
        """
        return self.quo_rem(right, integral=True)[0]

    def __getitem__(self, n):
        r"""
        Return the coefficient of `p^n` in the series expansion of this
        element, as an integer in the range `0` to `p-1`.

        EXAMPLES::

            sage: R = Zp(7,4,'capped-rel','series'); a = R(1/3); a
            5 + 4*7 + 4*7^2 + 4*7^3 + O(7^4)
            sage: a[0]  # indirect doctest
            doctest:warning
            ...
            DeprecationWarning: __getitem__ is changing to match the behavior of number fields. Please use expansion instead.
            See https://github.com/sagemath/sage/issues/14825 for details.
            5
            sage: a[1]
            4

        Negative indices do not have the special meaning they have for regular
        python lists. In the following example, ``a[-1]`` is simply the
        coefficient of `7^{-1}`::

            sage: K = Qp(7,4,'capped-rel')
            sage: b = K(1/7 + 7); b
            7^-1 + 7 + O(7^3)
            sage: b[-2]
            0
            sage: b[-1]
            1
            sage: b[0]
            0
            sage: b[1]
            1
            sage: b[2]
            0

        It is an error to access coefficients which are beyond the precision
        bound::

            sage: b[3]
            Traceback (most recent call last):
            ...
            PrecisionError
            sage: b[-2]
            0

        Slices also work::

            sage: a[0:2]
            5 + 4*7 + O(7^2)
            sage: a[-1:3:2]
            5 + 4*7^2 + O(7^3)
            sage: b[0:2]
            7 + O(7^2)
            sage: b[-1:3:2]
            7^-1 + 7 + O(7^3)

        If the slice includes coefficients which are beyond the precision
        bound, they are ignored. This is similar to the behaviour of slices of
        python lists::

            sage: a[3:7]
            4*7^3 + O(7^4)
            sage: b[3:7]
            O(7^3)

        For extension elements, "zeros" match the behavior of
        ``list``::

            sage: # needs sage.libs.ntl
            sage: S.<a> = Qq(125)
            sage: a[-2]                                                                 # needs sage.rings.padics
            []

        .. SEEALSO::

            :meth:`sage.rings.padics.local_generic_element.LocalGenericElement.slice`
        """
        from sage.misc.superseded import deprecation_cython as deprecation
        deprecation(14825, "__getitem__ is changing to match the behavior of number fields. Please use expansion instead.")
        return self.expansion(n)

    def __invert__(self):
        r"""
        Return the multiplicative inverse of ``self``.

        EXAMPLES::

            sage: R = Zp(7,4,'capped-rel','series'); a = R(3); a
            3 + O(7^4)
            sage: ~a  # indirect doctest
            5 + 4*7 + 4*7^2 + 4*7^3 + O(7^4)

        .. NOTE::

            The element returned is an element of the fraction field.
        """
        return ~self.parent().fraction_field()(self, relprec = self.precision_relative())

    cpdef _mod_(self, right):
        """
        If ``self`` is in a field, returns 0.  If in a ring, returns a
        `p`-adic integer such that

        (1) a = (a // b) * b + a % b

        holds.

        WARNING: The series expansion of a % b continues above the
        valuation of b.

        The definitions of a // b and a % b are intertwined by
        equation (1).  If a and b lie in a field, then setting a % b =
        0 and a // b = a / b provides a canonical way of satisfying
        this equation.

        However, for elements of integer rings, there are many choices
        of definitions for a // b and a % b that satisfy this
        equation.  Since `p`-adic rings in Sage come equipped with a
        uniformizer pi, we can use the choice of uniformizer in our
        definitions.  Here are some other criteria we might ask for:

        (2) If b = pi^k, the series expansion (in terms of pi) of a //
        b is just the series expansion of a, shifted over by k terms.

        (2') The series expansion of a % pi^k has no terms above
        pi^(k-1).

        The conditions (2) and (2') are equivalent.  But when we
        generalize these conditions to arbitrary b they diverge.

        (3) For general b, the series expansion of a // b is just the
        series expansion of a / b, truncating terms with negative
        exponents of pi.

        (4) For general b, the series expansion of a % b has no terms
        above b.valuation() - 1.

        In order to satisfy (3), one defines

        a // b = (a / b.unit_part()) >> b.valuation()
        a % b = a - (a // b) * b

        In order to satisfy (4), one defines

        a % b = a.lift() % pi.lift()^b.valuation()
        a // b = ((a - a % b) >> b.valuation()) / b.unit_part()


        In Sage we choose option (4) because it has better precision behavior.

        EXAMPLES::

            sage: R = ZpCA(5); a = R(129378); b = R(2398125)
            sage: a % b
            3 + O(5^20)
        """
        return self._quo_rem(right)[1]

    #def _is_exact_zero(self):
    #    return False

    #def _is_inexact_zero(self):
    #    return self.is_zero() and not self._is_exact_zero()

    def str(self, mode=None):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Zp(5,5,print_mode='bars')(1/3).str()[3:]
            '1|3|1|3|2'
        """
        return self._repr_(mode=mode)

    def _repr_(self, mode=None, do_latex=False):
        r"""
        Return a string representation of this element.

        INPUT:

        - ``mode`` -- allows one to override the default print mode of
          the parent (default: ``None``)

        - ``do_latex`` -- whether to return a latex representation or
          a normal one

        EXAMPLES::

            sage: Zp(5,5)(1/3)  # indirect doctest
            2 + 3*5 + 5^2 + 3*5^3 + 5^4 + O(5^5)

        We check that :issue:`26479` is fixed::

            sage: # needs sage.libs.ntl
            sage: x = polygen(ZZ, 'x')
            sage: K.<pi> = Qp(2).extension(x^3 - 2)
            sage: latex(pi)                                                             # needs sage.symbolic
            \pi + O(\pi^{61})
        """
        return self.parent()._printer.repr_gen(self, do_latex, mode=mode)

    def additive_order(self, prec=None):
        r"""
        Return the additive order of this element truncated
        at precision ``prec``.

        INPUT:

        - ``prec`` -- integer or ``None`` (default: ``None``)

        OUTPUT: the additive order of this element

        EXAMPLES::

            sage: R = Zp(7, 4, 'capped-rel', 'series'); a = R(7^3); a.additive_order(3)
            1
            sage: a.additive_order(4)
            +Infinity
            sage: R = Zp(7, 4, 'fixed-mod', 'series'); a = R(7^5); a.additive_order(6)
            1
        """
        if self.is_zero(prec):
            return Integer(1)
        else:
            return infinity

    def artin_hasse_exp(self, prec=None, algorithm=None):
        r"""
        Return the Artin-Hasse exponential of this element.

        INPUT:

        - ``prec`` -- integer or ``None`` (default: ``None``);
          the desired precision on the result; if ``None``, the
          precision is derived from the precision on the input

        - ``algorithm`` -- ``'direct'``, ``'series'``, ``'newton'`` or
          ``None`` (default)

          The direct algorithm computes the Artin-Hasse exponential
          of `x`, namely ``AH(x)`` as

          .. MATH::

              AH(x) = \exp(x + \frac{x^p}{p} + \frac{x^{p^2}}{p^2} + \dots

          It runs roughly as fast as the computation of the exponential
          (since the computation of the argument is not that costly).

          The series algorithm computes the series defining the
          Artin-Hasse exponential and evaluates it.

          The ``'newton'`` algorithm solves the equation

          .. MATH::

              \log(AH(x)) = x + \frac{x^p}{p} + \frac{x^{p^2}}{p^2} + \dots

          using a Newton scheme. It runs roughly as fast as the computation
          of the logarithm.

          By default, we use the ``'direct'`` algorithm if a fast algorithm for
          computing the exponential is available.
          If not, we use the ``'newton'`` algorithm if a fast algorithm for
          computing the logarithm is available.
          Otherwise we switch to the ``'series'`` algorithm.

        OUTPUT: the Artin-Hasse exponential of this element

        See :wikipedia:`Artin-Hasse_exponential` for more information.

        EXAMPLES::

            sage: x = Zp(5)(45/7)
            sage: y = x.artin_hasse_exp(); y
            1 + 2*5 + 4*5^2 + 3*5^3 + 5^7 + 2*5^8 + 3*5^10 + 2*5^11 + 2*5^12 +
            2*5^13 + 5^14 + 3*5^17 + 2*5^18 + 2*5^19 + O(5^20)

            sage: y * (-x).artin_hasse_exp()
            1 + O(5^20)

        The function respects your precision::

            sage: x = Zp(3,30)(45/7)
            sage: x.artin_hasse_exp()
            1 + 2*3^2 + 3^4 + 2*3^5 + 3^6 + 2*3^7 + 2*3^8 + 3^9 + 2*3^10 + 3^11 +
            3^13 + 2*3^15 + 2*3^16 + 2*3^17 + 3^19 + 3^20 + 2*3^21 + 3^23 + 3^24 +
            3^26 + 3^27 + 2*3^28 + O(3^30)

        Unless you tell it not to::

            sage: x = Zp(3,30)(45/7)
            sage: x.artin_hasse_exp()
            1 + 2*3^2 + 3^4 + 2*3^5 + 3^6 + 2*3^7 + 2*3^8 + 3^9 + 2*3^10 + 3^11 +
            3^13 + 2*3^15 + 2*3^16 + 2*3^17 + 3^19 + 3^20 + 2*3^21 + 3^23 + 3^24 +
            3^26 + 3^27 + 2*3^28 + O(3^30)
            sage: x.artin_hasse_exp(10)
            1 + 2*3^2 + 3^4 + 2*3^5 + 3^6 + 2*3^7 + 2*3^8 + 3^9 + O(3^10)

        For precision 1 the function just returns 1 since the
        exponential is always a 1-unit::

            sage: x = Zp(3).random_element()
            sage: while x.dist(0) >= 1:
            ....:     x = Zp(3).random_element()
            sage: x.artin_hasse_exp(1)
            1 + O(3)

        TESTS:

        Using Theorem 2.5 of [Conr]_::

            sage: x1 = 5*Zp(5).random_element()
            sage: x2 = 5*Zp(5).random_element()
            sage: y1 = x1.artin_hasse_exp()
            sage: y2 = x2.artin_hasse_exp()
            sage: (y1 - y2).abs() == (x1 - x2).abs()
            True

        Comparing with the formal power series definition::

            sage: x = PowerSeriesRing(QQ, 'x', default_prec=82).gen()
            sage: AH = sum(x**(3**i)/(3**i) for i in range(5)).O(82).exp()
            sage: z = Zp(3)(33/7)
            sage: ahz = AH(z); ahz                                                      # needs sage.libs.ntl
            1 + 2*3 + 3^2 + 3^3 + 2*3^5 + 3^6 + 2*3^7 + 3^9 + 3^11 + 3^12 +
            3^13 + 3^14 + 2*3^15 + 3^16 + 2*3^18 + 2*3^19 + O(3^20)
            sage: ahz - z.artin_hasse_exp()                                             # needs sage.libs.ntl
            O(3^20)

        Out of convergence domain::

            sage: Zp(5)(1).artin_hasse_exp()
            Traceback (most recent call last):
            ...
            ValueError: Artin-Hasse exponential does not converge on this input

        AUTHORS:

        - Mitchell Owen, Sebastian Pancrantz (2012-02): initial version.

        - Xavier Caruso (2018-08): extend to any p-adic rings and fields
          and implement several algorithms.
        """
        if self.valuation() < 1:
            raise ValueError("Artin-Hasse exponential does not converge on this input")
        R = self.parent()
        if prec is None:
            prec = min(self.precision_absolute(), R.precision_cap())
        else:
            prec = min(prec, self.precision_absolute(), R.precision_cap())

        if algorithm is None:
            try:
                R(0).exp(1, algorithm='binary_splitting')  # we check that binary splitting is available
                ans = self._AHE_direct(prec, exp_algorithm='binary_splitting')
            except NotImplementedError:
                try:
                    R(1).log(1, algorithm='binary_splitting')  # we check that binary splitting is available
                    ans = self._AHE_newton(prec, log_algorithm='binary_splitting')
                except NotImplementedError:
                    ans = self._AHE_series(prec)
        elif algorithm == 'direct':
            ans = self._AHE_direct(prec)
        elif algorithm == 'series':
            ans = self._AHE_series(prec)
        elif algorithm == 'newton':
            ans = self._AHE_newton(prec)
        else:
            raise ValueError("Algorithm must be 'direct', 'series', 'newton' or None")
        return ans

    def _AHE_direct(self, prec, exp_algorithm=None):
        r"""
        Return the Artin-Hasse exponential of this element.

        If `x` denotes the input element, its Artin-Hasse
        exponential is computed by taking the exponential of

        .. MATH::

            x + \frac{x^p}{p} + \frac{x^{p^2}}{p^2} + \cdots.

        INPUT:

        - ``prec`` -- integer; the precision at which the
          result should be computed

        - ``exp_algorithm`` -- string; the algorithm called
          for computing the exponential

        EXAMPLES::

            sage: W = Zp(3,10)
            sage: W(123456).artin_hasse_exp(algorithm='direct')  # indirect doctest
            1 + 3 + 2*3^3 + 2*3^4 + 3^5 + 2*3^6 + 2*3^7 + 3^8 + O(3^10)

        When `x^{p^i}/p^i` is not in the domain of convergence of the
        exponential for some nonnegative integer `i`, an error is raised::

            sage: # needs sage.libs.ntl
            sage: S.<x> = W[]
            sage: R.<pi> = W.extension(x^2 + 3)
            sage: pi.artin_hasse_exp(algorithm='direct')  # indirect doctest            # needs sage.symbolic
            Traceback (most recent call last):
            ...
            NotImplementedError: one factor of the Artin-Hasse exponential does not converge

        There is however an important exception.
        When we are working over `\ZZ_2` or `\QQ_2` and `x` is congruent to `2`
        modulo `4`, then `x` and `x^2/2` are not in the domain of convergence of
        the exponential. However, `\exp(x + x^2/2)` does converge.
        In this case, the Artin-Hasse exponential of `x`, denoted by `AH(x)`, is

        .. MATH::

            AH(x) = - \exp(x + \frac{x^2}{2} + \frac{x^4}{4} + \dots)

        with a negative sign. This method knows about this fact and handles
        the computation correctly::

            sage: W = Zp(2,8)
            sage: x = W(1234); x
            2 + 2^4 + 2^6 + 2^7 + O(2^9)
            sage: y1 = x.artin_hasse_exp(algorithm='direct'); y1
            1 + 2 + 2^2 + 2^6 + O(2^8)
            sage: y2 = exp(x + x^2/2 + x^4/4 + x^8/8); y2
            1 + 2^3 + 2^4 + 2^5 + 2^7 + O(2^8)
            sage: y1 == -y2
            True
            sage: y1 == x.artin_hasse_exp(algorithm='series')                           # needs sage.symbolic
            True

        .. SEEALSO::

            :meth:`artin_hasse_exp`, :meth:`_AHE_series`, :meth:`_AHE_newton`
        """
        R = self.parent()
        p = R.prime()
        pow = self.add_bigoh(prec)
        arg = pow
        denom = 1
        trunc = prec
        if R.absolute_degree() == 1:
            # Special code for Zp and Qp
            while pow != 0:
                trunc += 1
                pow = (pow**p).add_bigoh(trunc)
                denom *= p
                arg += pow/denom
            AH = arg.exp(algorithm=exp_algorithm)
            if p == 2 and self.add_bigoh(2) == 2:
                AH = -AH
        else:
            e = R.absolute_e()
            ep = e // (p-1)
            while pow != 0:
                trunc += e
                pow = (pow**p).add_bigoh(trunc)
                denom *= p
                s = pow/denom
                if s.valuation() <= ep:
                    raise NotImplementedError("one factor of the Artin-Hasse exponential does not converge")
                arg += s
            AH = arg.exp(algorithm=exp_algorithm)
        return AH

    def _AHE_series(self, prec):
        r"""
        Return the Artin-Hasse exponential of this element.

        This method first evaluates the Artin-Hasse series

        .. MATH::

            AH(x) = \exp(x + \frac{x^p}{p} + \frac{x^{p^2}}{p^2} + \dots)

        at enough precision and then plugs in the input element.

        INPUT:

        - ``prec`` -- integer; the precision at which the
          result should be computed

        EXAMPLES::

            sage: W = Zp(3,10)
            sage: W(123456).artin_hasse_exp(algorithm='series')  # indirect doctest     # needs sage.symbolic
            1 + 3 + 2*3^3 + 2*3^4 + 3^5 + 2*3^6 + 2*3^7 + 3^8 + O(3^10)

            sage: # needs sage.libs.ntl
            sage: S.<x> = W[]
            sage: R.<pi> = W.extension(x^2 + 3)
            sage: pi.artin_hasse_exp(algorithm='series')  # indirect doctest            # needs sage.symbolic
            1 + pi + 2*pi^2 + 2*pi^3 + 2*pi^4 + 2*pi^10 + 2*pi^11 + pi^13 + pi^18 + pi^19 + O(pi^20)

        .. SEEALSO::

            :meth:`artin_hasse_exp`, :meth:`_AHE_direct`, :meth:`_AHE_newton`
        """
        R = self.parent()
        p = R.prime()
        e = R.absolute_e()

        # We compute the Artin-Hasse series at the requested precision
        L = _AHE_coefficients(p, prec, 1 + (prec-1)//e)
        # We evaluate it using Horner algorithm
        y = R(0)
        x = self.add_bigoh(prec)
        for i in range(prec-1, -1, -1):
            y = y*x + R(L[i])

        return y

    def _AHE_newton(self, prec, log_algorithm=None):
        r"""
        Return the Artin-Hasse exponential of this element.

        If `x` denotes the input element, its Artin-Hasse exponential
        is computed by solving the following equation in `y`

        .. MATH::

            \log(y) = x + \frac{x^p}{p} + \frac{x^{p^2}}{p^2} + \dots

        using a Newton scheme.

        The first approximation used for initializing the Newton iteration
        is computed using the ``series`` algorithm (see :meth:`_AHE_series`).

        INPUT:

        - ``prec`` -- integer; the precision at which the
          result should be computed

        EXAMPLES::

            sage: W = Zp(3,10)
            sage: W(123456).artin_hasse_exp(algorithm='newton')  # indirect doctest
            1 + 3 + 2*3^3 + 2*3^4 + 3^5 + 2*3^6 + 2*3^7 + 3^8 + O(3^10)

            sage: # needs sage.libs.ntl
            sage: S.<x> = W[]
            sage: R.<pi> = W.extension(x^2 + 3)
            sage: pi.artin_hasse_exp(algorithm='newton')  # indirect doctest            # needs sage.symbolic
            1 + pi + 2*pi^2 + 2*pi^3 + 2*pi^4 + 2*pi^10 + 2*pi^11 + pi^13 + pi^18 + pi^19 + O(pi^20)

        .. SEEALSO::

            :meth:`artin_hasse_exp`, :meth:`_AHE_direct`, :meth:`_AHE_series`
        """
        R = self.parent()
        p = R.prime()
        e = R.absolute_e()

        # Step 1:
        # We compute a sufficiently good approximation of the result
        # in order to bootstrap the Newton iteration

        # We compute the Artin-Hasse series at the requested precision
        ep = e // (p-1)
        startprec = min(prec, ep+1)
        L = _AHE_coefficients(p, startprec, 1)
        # We evaluate it using Horner algorithm
        y = R(0)
        x = self.add_bigoh(startprec)
        for i in range(startprec-1, -1, -1):
            y = y*x + R(L[i])

        # Step 2:
        # We use Newton iteration to solve the equation
        # log(AH(x)) = x + x^p/p + x^(p^2)/p^2 + ...

        # We compute b = 1 + x + x^p/p + x^(p^2)/p^2 + ...
        pow = self.add_bigoh(prec)
        b = 1 + pow
        denom = 1
        trunc = prec
        while pow != 0:
            trunc += e
            pow = (pow**p).add_bigoh(trunc)
            denom *= p
            b += pow/denom
        # We iterate the Newton scheme: y_(n+1) = y_n * (b - log(y_n))
        curprec = startprec
        while curprec < prec:
            if p == 2:
                curprec = 2*curprec - e
            else:
                curprec = 2*curprec
            y = y.lift_to_precision(min(prec,curprec))
            y *= b - y.log(algorithm=log_algorithm)

        return R(y)

    def minimal_polynomial(self, name='x', base=None):
        """
        Return the minimal polynomial of this element over ``base``.

        INPUT:

        - ``name`` -- string (default: ``'x'``); the name of the variable

        - ``base`` -- a ring (default: the base ring of the parent);
          the base ring over which the minimal polynomial is computed

        EXAMPLES::

            sage: Zp(5,5)(1/3).minimal_polynomial('x')                                  # needs sage.libs.ntl
            (1 + O(5^5))*x + 3 + 5 + 3*5^2 + 5^3 + 3*5^4 + O(5^5)

            sage: Zp(5,5)(1/3).minimal_polynomial('foo')                                # needs sage.libs.ntl
            (1 + O(5^5))*foo + 3 + 5 + 3*5^2 + 5^3 + 3*5^4 + O(5^5)

        ::

            sage: # needs sage.libs.ntl
            sage: K.<a> = QqCR(2^3,5)
            sage: S.<x> = K[]
            sage: L.<pi> = K.extension(x^4 - 2*a)
            sage: pi.minimal_polynomial()                                               # needs sage.symbolic
            (1 + O(2^5))*x^4 + a*2 + a*2^2 + a*2^3 + a*2^4 + a*2^5 + O(2^6)
            sage: (pi^2).minimal_polynomial()                                           # needs sage.symbolic
            (1 + O(2^5))*x^2 + a*2 + a*2^2 + a*2^3 + a*2^4 + a*2^5 + O(2^6)
            sage: (1/pi).minimal_polynomial()                                           # needs sage.symbolic
            (1 + O(2^5))*x^4 + (a^2 + 1)*2^-1 + O(2^4)
            sage: elt = L.random_element()
            sage: P = elt.minimal_polynomial()  # not tested, known bug (see :issue:`32111`)
            sage: P(elt) == 0  # not tested
            True
        """
        parent = self.parent()
        R = parent.base_ring()
        if base is None:
            base = R
        polring = base[name]
        if base is parent:
            return polring([-self,1])
        elif base is R:
            from sage.modules.free_module import VectorSpace
            L = parent.fraction_field()
            K = base.fraction_field()
            deg = L.relative_degree()
            V = VectorSpace(K, deg)
            vector = [K(1)] + (deg-1)*[K(0)]
            vectors = [vector]
            W = V.span(vectors)
            elt = self
            while True:
                poly = elt.polynomial()
                vector = V([ poly[i] for i in range(deg) ])
                if vector in W: break
                vectors.append(vector)
                W += V.span([vector])
                elt *= self
            W = V.span_of_basis(vectors)
            coeffs = [ -c for c in W.coordinate_vector(vector) ] + [K(1)]
            return polring(coeffs)
        else:
            raise NotImplementedError

    def norm(self, base=None):
        """
        Return the norm of this `p`-adic element over ``base``.

        .. WARNING::

            This is not the `p`-adic absolute value.  This is a field
            theoretic norm down to a base ring.  If you want the
            `p`-adic absolute value, use the method :meth:`abs`
            instead.

        INPUT:

        - ``base`` -- a subring of the parent (default: base ring)

        OUTPUT: the norm of this `p`-adic element over the given base

        EXAMPLES::

            sage: Zp(5)(5).norm()                                                       # needs sage.libs.ntl
            5 + O(5^21)

        ::

            sage: # needs sage.libs.ntl
            sage: K.<a> = QqCR(2^3,5)
            sage: S.<x> = K[]
            sage: L.<pi> = K.extension(x^4 - 2*a)
            sage: pi.norm()  # norm over K                                              # needs sage.symbolic
            a*2 + a*2^2 + a*2^3 + a*2^4 + a*2^5 + O(2^6)
            sage: (pi^2).norm()                                                         # needs sage.symbolic
            a^2*2^2 + O(2^7)
            sage: pi.norm()^2                                                           # needs sage.symbolic
            a^2*2^2 + O(2^7)

        TESTS::

            sage: x = L.random_element()                                                # needs sage.libs.ntl
            sage: y = L.random_element()                                                # needs sage.libs.ntl
            sage: (x*y).norm() == x.norm() * y.norm()  # not tested, known bug (see :issue:`32085`)
            True
        """
        parent = self.parent()
        if base is None:
            base = parent.base_ring()
        poly = self.minimal_polynomial(base=base)
        polydeg = poly.degree()
        extdeg = parent.absolute_degree() // (base.absolute_degree() * polydeg)
        return ((-1)**polydeg * poly[0]) ** extdeg

    def trace(self, base=None):
        """
        Return the trace of this `p`-adic element over the base ring.

        INPUT:

        - ``base`` -- a subring of the parent (default: base ring)

        OUTPUT: the trace of this `p`-adic element over the given base

        EXAMPLES::

            sage: Zp(5,5)(5).trace()                                                    # needs sage.libs.ntl
            5 + O(5^6)

            sage: # needs sage.libs.ntl
            sage: K.<a> = QqCR(2^3,7)
            sage: S.<x> = K[]
            sage: L.<pi> = K.extension(x^4 - 4*a*x^3 + 2*a)
            sage: pi.trace()  # trace over K                                            # needs sage.symbolic
            a*2^2 + O(2^8)
            sage: (pi+1).trace()                                                        # needs sage.symbolic
            (a + 1)*2^2 + O(2^7)

        TESTS::

            sage: x = L.random_element()                                                # needs sage.libs.ntl
            sage: y = L.random_element()                                                # needs sage.libs.ntl
            sage: (x+y).trace() == x.trace() + y.trace()  # not tested, known bug (see :issue:`32085`)
            True
        """
        parent = self.parent()
        if base is None:
            base = parent.base_ring()
        poly = self.minimal_polynomial(base=base)
        polydeg = poly.degree()
        extdeg = parent.absolute_degree() // (base.absolute_degree() * polydeg)
        return -extdeg * poly[polydeg-1]

    def algdep(self, n):
        """
        Return a polynomial of degree at most `n` which is approximately
        satisfied by this number. Note that the returned polynomial need not be
        irreducible, and indeed usually won't be if this number is a good
        approximation to an algebraic number of degree less than `n`.

        ALGORITHM: Uses the PARI C-library :pari:`algdep` command.

        INPUT:

        - ``self`` -- a `p`-adic element
        - ``n`` -- integer

        OUTPUT: polynomial; degree `n` polynomial approximately satisfied by ``self``

        EXAMPLES::

            sage: K = Qp(3,20,'capped-rel','series'); R = Zp(3,20,'capped-rel','series')
            sage: a = K(7/19); a
            1 + 2*3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^8 + 2*3^9 + 3^11 + 3^12
              + 2*3^15 + 2*3^16 + 3^17 + 2*3^19 + O(3^20)
            sage: a.algdep(1)
            19*x - 7
            sage: K2 = Qp(7,20,'capped-rel')
            sage: b = K2.zeta(); b.algdep(2)
            x^2 - x + 1
            sage: K2 = Qp(11,20,'capped-rel')
            sage: b = K2.zeta(); b.algdep(4)
            x^4 - x^3 + x^2 - x + 1
            sage: a = R(7/19); a
            1 + 2*3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^8 + 2*3^9 + 3^11 + 3^12
              + 2*3^15 + 2*3^16 + 3^17 + 2*3^19 + O(3^20)
            sage: a.algdep(1)
            19*x - 7
            sage: R2 = Zp(7,20,'capped-rel')
            sage: b = R2.zeta(); b.algdep(2)
            x^2 - x + 1
            sage: R2 = Zp(11,20,'capped-rel')
            sage: b = R2.zeta(); b.algdep(4)
            x^4 - x^3 + x^2 - x + 1
        """
        # TODO: figure out if this works for extension rings.  If not, move this to padic_base_generic_element.
        from sage.arith.misc import algdep
        return algdep(self, n)

    def algebraic_dependency(self, n):
        """
        Return a polynomial of degree at most `n` which is approximately
        satisfied by this number.  Note that the returned polynomial need not
        be irreducible, and indeed usually won't be if this number is a good
        approximation to an algebraic number of degree less than `n`.

        ALGORITHM: Uses the PARI C-library :pari:`algdep` command.

        INPUT:

        - ``self`` -- a `p`-adic element
        - ``n`` -- integer

        OUTPUT: polynomial; degree `n` polynomial approximately satisfied by ``self``

        EXAMPLES::

            sage: K = Qp(3,20,'capped-rel','series'); R = Zp(3,20,'capped-rel','series')
            sage: a = K(7/19); a
            1 + 2*3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^8 + 2*3^9 + 3^11 + 3^12
              + 2*3^15 + 2*3^16 + 3^17 + 2*3^19 + O(3^20)
            sage: a.algebraic_dependency(1)
            19*x - 7
            sage: K2 = Qp(7,20,'capped-rel')
            sage: b = K2.zeta(); b.algebraic_dependency(2)
            x^2 - x + 1
            sage: K2 = Qp(11,20,'capped-rel')
            sage: b = K2.zeta(); b.algebraic_dependency(4)
            x^4 - x^3 + x^2 - x + 1
            sage: a = R(7/19); a
            1 + 2*3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^8 + 2*3^9 + 3^11 + 3^12
              + 2*3^15 + 2*3^16 + 3^17 + 2*3^19 + O(3^20)
            sage: a.algebraic_dependency(1)
            19*x - 7
            sage: R2 = Zp(7,20,'capped-rel')
            sage: b = R2.zeta(); b.algebraic_dependency(2)
            x^2 - x + 1
            sage: R2 = Zp(11,20,'capped-rel')
            sage: b = R2.zeta(); b.algebraic_dependency(4)
            x^4 - x^3 + x^2 - x + 1
        """
        return self.algdep(n)

    #def exp_artin_hasse(self):
    #    """
    #    Returns the Artin-Hasse exponential of self.

    #    This is defined by: E_p(x) = exp(x + x^p/p + x^(p^2)/p^2 + ...)
    #    """
    #    raise NotImplementedError

    def dwork_expansion(self, bd=20, a=0):
        r"""
        Return the value of a function defined by Dwork.

        Used to compute the `p`-adic Gamma function, see :meth:`gamma`.

        INPUT:

        - ``bd`` -- integer (default: 20); precision bound
        - ``a`` -- integer (default: 0); offset parameter

        OUTPUT: a `p`-adic integer

        .. NOTE::

            This is based on GP code written by Fernando Rodriguez
            Villegas (http://www.ma.utexas.edu/cnt/cnt-frames.html).
            William Stein sped it up for GP
            (http://sage.math.washington.edu/home/wstein/www/home/wbhart/pari-2.4.2.alpha/src/basemath/trans2.c).
            The output is a `p`-adic integer from Dwork's expansion,
            used to compute the `p`-adic gamma function as in [RV2007]_
            section 6.2.
            The coefficients of the expansion are now cached to speed up
            multiple evaluation, as in the trace formula for hypergeometric
            motives.

        EXAMPLES::

            sage: R = Zp(17)
            sage: x = R(5+3*17+13*17^2+6*17^3+12*17^5+10*17^(14)+5*17^(17)+O(17^(19)))
            sage: x.dwork_expansion(18)
            16 + 7*17 + 11*17^2 + 4*17^3 + 8*17^4 + 10*17^5 + 11*17^6 + 6*17^7
            + 17^8 + 8*17^10 + 13*17^11 + 9*17^12 + 15*17^13  + 2*17^14 + 6*17^15
            + 7*17^16 + 6*17^17 + O(17^18)

            sage: R = Zp(5)
            sage: x = R(3*5^2+4*5^3+1*5^4+2*5^5+1*5^(10)+O(5^(20)))
            sage: x.dwork_expansion()
            4 + 4*5 + 4*5^2 + 4*5^3 + 2*5^4 + 4*5^5 + 5^7 + 3*5^9 + 4*5^10 + 3*5^11
            + 5^13 + 4*5^14 + 2*5^15 + 2*5^16 + 2*5^17 + 3*5^18 + O(5^20)

        TESTS:

        This test was added in :issue:`24433`::

            sage: F = Qp(7)
            sage: F(4).gamma()
            6 + O(7^20)
            sage: -F(1).dwork_expansion(a=3)
            6 + 4*7^19 + O(7^20)
        """
        R = self.parent()
        p = R.prime()

        # If p == 2, must work in Qp rather than Zp.
        if p == 2 and not R.is_field():
            S = R.fraction_field()
            return R(S(self).dwork_expansion(bd, a))
        try:
            v = R.dwork_coeffs
            if len(v) < p*bd:
                raise AttributeError
        except AttributeError:
            v = dwork_mahler_coeffs(R, bd)
        return evaluate_dwork_mahler(v, self, p, bd, a)

    def gamma(self, algorithm='pari'):
        r"""
        Return the value of the `p`-adic Gamma function.

        INPUT:

        - ``algorithm`` -- string. Can be set to ``'pari'`` to call
          the PARI function, or ``'sage'`` to call the function
          implemented in Sage. The default is ``'pari'`` since
          PARI is about 10 times faster than Sage.

        OUTPUT: a `p`-adic integer

        .. NOTE::

            This is based on GP code written by Fernando Rodriguez
            Villegas (http://www.ma.utexas.edu/cnt/cnt-frames.html).
            William Stein sped it up for GP
            (http://sage.math.washington.edu/home/wstein/www/home/wbhart/pari-2.4.2.alpha/src/basemath/trans2.c).
            The ``'sage'`` version uses dwork_expansion() to compute the
            `p`-adic gamma function of ``self`` as in [RV2007]_ section 6.2.

        EXAMPLES:

        This example illustrates ``x.gamma()`` for `x` a `p`-adic unit::

            sage: R = Zp(7)
            sage: x = R(2+3*7^2+4*7^3+O(7^20))
            sage: x.gamma('pari')
            1 + 2*7^2 + 4*7^3 + 5*7^4 + 3*7^5 + 7^8 + 7^9 + 4*7^10 + 3*7^12
            + 7^13 + 5*7^14 + 3*7^15 + 2*7^16 + 2*7^17 + 5*7^18 + 4*7^19 + O(7^20)
            sage: x.gamma('sage')
            1 + 2*7^2 + 4*7^3 + 5*7^4 + 3*7^5 + 7^8 + 7^9 + 4*7^10 + 3*7^12
            + 7^13 + 5*7^14 + 3*7^15 + 2*7^16 + 2*7^17 + 5*7^18 + 4*7^19 + O(7^20)
            sage: x.gamma('pari') == x.gamma('sage')
            True

        Now ``x.gamma()`` for `x` a `p`-adic integer but not a unit::

            sage: R = Zp(17)
            sage: x = R(17+17^2+3*17^3+12*17^8+O(17^13))
            sage: x.gamma('pari')
            1 + 12*17 + 13*17^2 + 13*17^3 + 10*17^4 + 7*17^5 + 16*17^7
            + 13*17^9 + 4*17^10 + 9*17^11 + 17^12 + O(17^13)
            sage: x.gamma('sage')
            1 + 12*17 + 13*17^2 + 13*17^3 + 10*17^4 + 7*17^5 + 16*17^7
            + 13*17^9 + 4*17^10 + 9*17^11 + 17^12 + O(17^13)
            sage: x.gamma('pari') == x.gamma('sage')
            True

        Finally, this function is not defined if `x` is not a `p`-adic integer::

            sage: K = Qp(7)
            sage: x = K(7^-5 + 2*7^-4 + 5*7^-3 + 2*7^-2 + 3*7^-1 + 3 + 3*7
            ....:       + 7^3 + 4*7^4 + 5*7^5 + 6*7^8 + 3*7^9 + 6*7^10 + 5*7^11 + 6*7^12
            ....:       + 3*7^13 + 5*7^14 + O(7^15))
            sage: x.gamma()
            Traceback (most recent call last):
            ...
            ValueError: The p-adic gamma function only works on elements of Zp

        TESTS:

        We check that :issue:`23784` is resolved::

            sage: Zp(5)(0).gamma()
            1 + O(5^20)

        Check the cached version of ``dwork_expansion`` from :issue:`24433`::

            sage: p = next_prime(200)
            sage: F = Qp(p)
            sage: l1 = [F(a/(p-1)).gamma(algorithm='pari') for a in range(p-1)]     # long time
            sage: l2 = [F(a/(p-1)).gamma(algorithm='sage') for a in range(p-1)]     # long time
            sage: all(l1[i] == l2[i] for i in range(p-1))  # long time
            True

        The `p`-adic Gamma function has anomalous behavior for the prime 2::

            sage: F = Qp(2)
            sage: x = F(-1) + O(2^2)
            sage: x.gamma(algorithm='pari')
            1 + O(2)
            sage: x.gamma(algorithm='sage')
            1 + O(2)
            sage: x = F(-1) + O(2^3)
            sage: x.gamma(algorithm='pari')
            1 + O(2^3)
            sage: x.gamma(algorithm='sage')
            1 + O(2^3)
        """
        if self.parent().absolute_degree() > 1 or self.valuation() < 0:
            raise ValueError('The p-adic gamma function only works '
                             'on elements of Zp')
        parent = self.parent()
        n = self.precision_absolute()
        if n is infinity:
            # Have to deal with exact zeros separately
            return parent(1)
        if algorithm == 'pari':
            return parent(self.__pari__().gamma())
        elif algorithm == 'sage':
            p = parent.prime()
            bd = -((-n*p) // (p-1))
            k = (-self) % p
            x = (self+k) >> 1
            return -x.dwork_expansion(bd, k.lift())

    @coerce_binop
    def gcd(self, other):
        r"""
        Return a greatest common divisor of ``self`` and ``other``.

        INPUT:

        - ``other`` -- an element in the same ring as ``self``

        AUTHORS:

        - Julian Rueth (2012-10-19): initial version

        .. NOTE::

            Since the elements are only given with finite precision,
            their greatest common divisor is in general not unique (not even up
            to units). For example `O(3)` is a representative for the elements
            0 and 3 in the 3-adic ring `\ZZ_3`. The greatest common
            divisor of `O(3)` and `O(3)` could be (among others) 3 or 0 which
            have different valuation. The algorithm implemented here, will
            return an element of minimal valuation among the possible greatest
            common divisors.

        EXAMPLES:

        The greatest common divisor is either zero or a power of the
        uniformizing parameter::

            sage: R = Zp(3)
            sage: R.zero().gcd(R.zero())
            0
            sage: R(3).gcd(9)
            3 + O(3^21)

        A nonzero result is always lifted to the maximal precision possible in
        the ring::

            sage: a = R(3,2); a
            3 + O(3^2)
            sage: b = R(9,3); b
            3^2 + O(3^3)
            sage: a.gcd(b)
            3 + O(3^21)
            sage: a.gcd(0)
            3 + O(3^21)

        If both elements are zero, then the result is zero with the precision
        set to the smallest of their precisions::

            sage: a = R.zero(); a
            0
            sage: b = R(0,2); b
            O(3^2)
            sage: a.gcd(b)
            O(3^2)

        One could argue that it is mathematically correct to return `9 +
        O(3^{22})` instead. However, this would lead to some confusing
        behaviour::

            sage: alternative_gcd = R(9,22); alternative_gcd
            3^2 + O(3^22)
            sage: a.is_zero()
            True
            sage: b.is_zero()
            True
            sage: alternative_gcd.is_zero()
            False

        If exactly one element is zero, then the result depends on the
        valuation of the other element::

            sage: R(0,3).gcd(3^4)
            O(3^3)
            sage: R(0,4).gcd(3^4)
            O(3^4)
            sage: R(0,5).gcd(3^4)
            3^4 + O(3^24)

        Over a field, the greatest common divisor is either zero (possibly with
        finite precision) or one::

            sage: K = Qp(3)
            sage: K(3).gcd(0)
            1 + O(3^20)
            sage: K.zero().gcd(0)
            0
            sage: K.zero().gcd(K(0,2))
            O(3^2)
            sage: K(3).gcd(4)
            1 + O(3^20)

        TESTS:

        The implementation also works over extensions::

            sage: # needs sage.libs.ntl
            sage: K = Qp(3)
            sage: R.<a> = K[]
            sage: L.<a> = K.extension(a^3-3)
            sage: (a+3).gcd(3)
            1 + O(a^60)
            sage: R = Zp(3)
            sage: S.<a> = R[]
            sage: S.<a> = R.extension(a^3-3)
            sage: (a+3).gcd(3)
            a + O(a^61)
            sage: K = Qp(3)
            sage: R.<a> = K[]
            sage: L.<a> = K.extension(a^2-2)
            sage: (a+3).gcd(3)
            1 + O(3^20)
            sage: R = Zp(3)
            sage: S.<a> = R[]
            sage: S.<a> = R.extension(a^2-2)
            sage: (a+3).gcd(3)
            1 + O(3^20)

        For elements with a fixed modulus::

            sage: R = ZpFM(3)
            sage: R(3).gcd(9)
            3

        And elements with a capped absolute precision::

            sage: R = ZpCA(3)
            sage: R(3).gcd(9)
            3 + O(3^20)
        """
        if self.is_zero() and other.is_zero():
            if self.valuation() < other.valuation():
                return self
            else:
                return other

        if self.parent().is_field():
            return self.parent().one()

        if min(self.valuation(),other.valuation()) >= min(self.precision_absolute(),other.precision_absolute()):
            return self.parent().zero().add_bigoh(min(self.precision_absolute(),other.precision_absolute()))

        return self.parent().uniformiser_pow( min(self.valuation(),other.valuation()) )

    @coerce_binop
    def xgcd(self, other):
        r"""
        Compute the extended gcd of this element and ``other``.

        INPUT:

        - ``other`` -- an element in the same ring

        OUTPUT:

        A tuple ``r``, ``s``, ``t`` such that ``r`` is a greatest common
        divisor of this element and ``other`` and ``r = s*self + t*other``.

        AUTHORS:

        - Julian Rueth (2012-10-19): initial version

        .. NOTE::

            Since the elements are only given with finite precision, their
            greatest common divisor is in general not unique (not even up to
            units). For example `O(3)` is a representative for the elements 0
            and 3 in the 3-adic ring `\ZZ_3`. The greatest common
            divisor of `O(3)` and `O(3)` could be (among others) 3 or 0 which
            have different valuation. The algorithm implemented here, will
            return an element of minimal valuation among the possible greatest
            common divisors.

        EXAMPLES:

        The greatest common divisor is either zero or a power of the
        uniformizing parameter::

            sage: R = Zp(3)
            sage: R.zero().xgcd(R.zero())
            (0, 1 + O(3^20), 0)
            sage: R(3).xgcd(9)
            (3 + O(3^21), 1 + O(3^20), 0)

        Unlike for :meth:`gcd`, the result is not lifted to the maximal
        precision possible in the ring; it is such that ``r = s*self +
        t*other`` holds true::

            sage: a = R(3,2); a
            3 + O(3^2)
            sage: b = R(9,3); b
            3^2 + O(3^3)
            sage: a.xgcd(b)
            (3 + O(3^2), 1 + O(3), 0)
            sage: a.xgcd(0)
            (3 + O(3^2), 1 + O(3), 0)

        If both elements are zero, then the result is zero with
        the precision set to the smallest of their precisions::

            sage: a = R.zero(); a
            0
            sage: b = R(0,2); b
            O(3^2)
            sage: a.xgcd(b)
            (O(3^2), 0, 1 + O(3^20))

        If only one element is zero, then the result depends on its precision::

            sage: # needs sage.rings.padics
            sage: R(9).xgcd(R(0,1))
            (O(3), 0, 1 + O(3^20))
            sage: R(9).xgcd(R(0,2))
            (O(3^2), 0, 1 + O(3^20))
            sage: R(9).xgcd(R(0,3))
            (3^2 + O(3^22), 1 + O(3^20), 0)
            sage: R(9).xgcd(R(0,4))
            (3^2 + O(3^22), 1 + O(3^20), 0)

        Over a field, the greatest common divisor is either zero (possibly with
        finite precision) or one::

            sage: K = Qp(3)
            sage: K(3).xgcd(0)
            (1 + O(3^20), 3^-1 + O(3^19), 0)
            sage: K.zero().xgcd(0)
            (0, 1 + O(3^20), 0)
            sage: K.zero().xgcd(K(0,2))
            (O(3^2), 0, 1 + O(3^20))
            sage: K(3).xgcd(4)
            (1 + O(3^20), 3^-1 + O(3^19), 0)

        TESTS:

        The implementation also works over extensions::

            sage: # needs sage.libs.ntl
            sage: K = Qp(3)
            sage: R.<a> = K[]
            sage: L.<a> = K.extension(a^3-3)
            sage: (a+3).xgcd(3)
            (1 + O(a^60),
             a^-1 + 2*a + a^3 + 2*a^4 + 2*a^5 + 2*a^8 + 2*a^9
              + 2*a^12 + 2*a^13 + 2*a^16 + 2*a^17 + 2*a^20 + 2*a^21 + 2*a^24
              + 2*a^25 + 2*a^28 + 2*a^29 + 2*a^32 + 2*a^33 + 2*a^36 + 2*a^37
              + 2*a^40 + 2*a^41 + 2*a^44 + 2*a^45 + 2*a^48 + 2*a^49 + 2*a^52
              + 2*a^53 + 2*a^56 + 2*a^57 + O(a^59),
             0)
            sage: R = Zp(3)
            sage: S.<a> = R[]
            sage: S.<a> = R.extension(a^3-3)
            sage: (a+3).xgcd(3)
            (a + O(a^61),
             1 + 2*a^2 + a^4 + 2*a^5 + 2*a^6 + 2*a^9 + 2*a^10
              + 2*a^13 + 2*a^14 + 2*a^17 + 2*a^18 + 2*a^21 + 2*a^22 + 2*a^25
              + 2*a^26 + 2*a^29 + 2*a^30 + 2*a^33 + 2*a^34 + 2*a^37 + 2*a^38
              + 2*a^41 + 2*a^42 + 2*a^45 + 2*a^46 + 2*a^49 + 2*a^50 + 2*a^53
              + 2*a^54 + 2*a^57 + 2*a^58 + O(a^60),
             0)
            sage: K = Qp(3)
            sage: R.<a> = K[]
            sage: L.<a> = K.extension(a^2-2)
            sage: (a+3).xgcd(3)
            (1 + O(3^20),
             2*a + (a + 1)*3 + (2*a + 1)*3^2 + (a + 2)*3^4 + 3^5
              + (2*a + 2)*3^6 + a*3^7 + (2*a + 1)*3^8 + (a + 2)*3^10 + 3^11
              + (2*a + 2)*3^12 + a*3^13 + (2*a + 1)*3^14 + (a + 2)*3^16
              + 3^17 + (2*a + 2)*3^18 + a*3^19 + O(3^20),
             0)
            sage: R = Zp(3)
            sage: S.<a> = R[]
            sage: S.<a> = R.extension(a^2-2)
            sage: (a+3).xgcd(3)
            (1 + O(3^20),
             2*a + (a + 1)*3 + (2*a + 1)*3^2 + (a + 2)*3^4 + 3^5
              + (2*a + 2)*3^6 + a*3^7 + (2*a + 1)*3^8 + (a + 2)*3^10 + 3^11
              + (2*a + 2)*3^12 + a*3^13 + (2*a + 1)*3^14 + (a + 2)*3^16 + 3^17
              + (2*a + 2)*3^18 + a*3^19 + O(3^20),
             0)

        For elements with a fixed modulus::

            sage: R = ZpFM(3)
            sage: R(3).xgcd(9)
            (3, 1, 0)

        And elements with a capped absolute precision::

            sage: R = ZpCA(3)
            sage: R(3).xgcd(9)
            (3 + O(3^20), 1 + O(3^19), O(3^20))
        """
        s,t = self.parent().zero(), self.parent().zero()
        if self.is_zero() and other.is_zero():
            if self.valuation() <= other.valuation():
                s = self.parent().one()
            else:
                t = self.parent().one()
        elif self.parent().is_field():
            if not self.is_zero():
                s = ~self
            else:
                t = ~other
        elif self.valuation() < other.valuation():
            if self.is_zero():
                s = self.parent().one()
            else:
                s = self.unit_part().inverse_of_unit()
        else:
            if other.is_zero():
                t = self.parent().one()
            else:
                t = other.unit_part().inverse_of_unit()

        return s*self+t*other,s,t

    def is_square(self):
        """
        Return whether this element is a square.

        INPUT:

        - ``self`` -- a `p`-adic element

        EXAMPLES::

            sage: R = Zp(3,20,'capped-rel')
            sage: R(0).is_square()
            True
            sage: R(1).is_square()
            True
            sage: R(2).is_square()
            False

        TESTS::

            sage: R(3).is_square()
            False
            sage: R(4).is_square()
            True
            sage: R(6).is_square()
            False
            sage: R(9).is_square()
            True

            sage: R2 = Zp(2,20,'capped-rel')
            sage: R2(0).is_square()
            True
            sage: R2(1).is_square()
            True
            sage: R2(2).is_square()
            False
            sage: R2(3).is_square()
            False
            sage: R2(4).is_square()
            True
            sage: R2(5).is_square()
            False
            sage: R2(6).is_square()
            False
            sage: R2(7).is_square()
            False
            sage: R2(8).is_square()
            False
            sage: R2(9).is_square()
            True

            sage: K = Qp(3,20,'capped-rel')
            sage: K(0).is_square()
            True
            sage: K(1).is_square()
            True
            sage: K(2).is_square()
            False
            sage: K(3).is_square()
            False
            sage: K(4).is_square()
            True
            sage: K(6).is_square()
            False
            sage: K(9).is_square()
            True
            sage: K(1/3).is_square()
            False
            sage: K(1/9).is_square()
            True

            sage: K2 = Qp(2,20,'capped-rel')
            sage: K2(0).is_square()
            True
            sage: K2(1).is_square()
            True
            sage: K2(2).is_square()
            False
            sage: K2(3).is_square()
            False
            sage: K2(4).is_square()
            True
            sage: K2(5).is_square()
            False
            sage: K2(6).is_square()
            False
            sage: K2(7).is_square()
            False
            sage: K2(8).is_square()
            False
            sage: K2(9).is_square()
            True
            sage: K2(1/2).is_square()
            False
            sage: K2(1/4).is_square()
            True
        """
        if self._is_exact_zero():
            return True
        parent = self.parent()
        if parent.prime() != 2:
            if self.is_zero():
                raise PrecisionError("not enough precision to be sure that this element has a square root")
            return (self.valuation() % 2 == 0) and (self.unit_part().residue(1).is_square())
        else:
            e = parent.absolute_e()
            try:
                self.add_bigoh(self.valuation() + 2*e + 1).nth_root(2)
            except ValueError:
                return False
            except PrecisionError:
                raise PrecisionError("not enough precision to be sure that this element has a square root")
            return True

    def is_squarefree(self):
        r"""
        Return whether this element is squarefree, i.e., whether there exists
        no non-unit `g` such that `g^2` divides this element.

        EXAMPLES:

        The zero element is never squarefree::

            sage: K = Qp(2)
            sage: K.zero().is_squarefree()
            False

        In `p`-adic rings, only elements of valuation at most 1 are
        squarefree::

            sage: R = Zp(2)
            sage: R(1).is_squarefree()
            True
            sage: R(2).is_squarefree()
            True
            sage: R(4).is_squarefree()
            False

        This works only if the precision is known sufficiently well::

            sage: R(0,1).is_squarefree()
            Traceback (most recent call last):
            ...
            PrecisionError: element not known to sufficient precision to decide squarefreeness
            sage: R(0,2).is_squarefree()
            False
            sage: R(1,1).is_squarefree()
            True

        For fields we are not so strict about the precision and treat inexact
        zeros as the zero element::

            K(0,0).is_squarefree()
            False
        """
        if self.parent().is_field():
            if self.is_zero():
                return False
            return True
        else:
            v = self.valuation()
            if v >= 2:
                return False
            elif self.is_zero():
                raise PrecisionError("element not known to sufficient precision to decide squarefreeness")
            else:
                return True

    #def log_artin_hasse(self):
    #    raise NotImplementedError

    def multiplicative_order(self, prec=None):
        r"""
        Return the multiplicative order of ``self``, where ``self`` is
        considered to be one if it is one modulo `p^{\mbox{prec}}`.

        INPUT:

        - ``self`` -- a `p`-adic element
        - ``prec`` -- integer

        OUTPUT: integer; the multiplicative order of ``self``

        EXAMPLES::

            sage: K = Qp(5,20,'capped-rel')
            sage: K(-1).multiplicative_order(20)
            2
            sage: K(1).multiplicative_order(20)
            1
            sage: K(2).multiplicative_order(20)
            +Infinity
            sage: K(5).multiplicative_order(20)
            +Infinity
            sage: K(1/5).multiplicative_order(20)
            +Infinity
            sage: K.zeta().multiplicative_order(20)
            4

        Over unramified extensions::

            sage: # needs sage.libs.ntl
            sage: L1.<a> = Qq(5^3)
            sage: c = L1.teichmuller(a)
            sage: c.multiplicative_order()
            124
            sage: c^124
            1 + O(5^20)

        Over totally ramified extensions::

            sage: # needs sage.libs.ntl
            sage: x = polygen(ZZ, 'x')
            sage: L2.<pi> = Qp(5).extension(x^4 + 5*x^3 + 10*x^2 + 10*x + 5)
            sage: u = 1 + pi
            sage: u.multiplicative_order()
            5
            sage: v = L2.teichmuller(2)
            sage: v.multiplicative_order()
            4
            sage: (u*v).multiplicative_order()
            20

        TESTS::

            sage: R = Zp(5,20,'capped-rel')
            sage: R(-1).multiplicative_order(20)
            2
            sage: R(1).multiplicative_order(20)
            1
            sage: R(2).multiplicative_order(20)
            +Infinity
            sage: R(3).multiplicative_order(20)
            +Infinity
            sage: R(4).multiplicative_order(20)
            +Infinity
            sage: R(5).multiplicative_order(20)
            +Infinity
            sage: R(25).multiplicative_order(20)
            +Infinity
            sage: R.zeta().multiplicative_order(20)
            4
        """
        if prec is not None:
            self = self.add_bigoh(prec)
        if self == 0 or self.valuation() != 0:
            return infinity
        parent = self.parent()
        p = parent.prime()

        # Compute the multiplicative order outside p
        res = self.residue()
        order = res.multiplicative_order()
        self /= parent.teichmuller(self)
        if self == 1:
            return order

        # Compute multiplicative order at p
        e = parent.absolute_e()
        if not (p-1).divides(e):
            return infinity
        n = e.valuation(p)
        for _ in range(n+1):
            order *= p
            self = self**p
            if self == 1:
                return order
        return infinity

    def valuation(self, p=None):
        r"""
        Return the valuation of this element.

        INPUT:

        - ``self`` -- a `p`-adic element
        - ``p`` -- a prime (default: ``None``); if specified, will make sure
          that ``p == self.parent().prime()``

        .. NOTE::

            The optional argument `p` is used for consistency with the valuation
            methods on integers and rationals.

        OUTPUT: integer; the valuation of ``self``

        EXAMPLES::

            sage: R = Zp(17, 4,'capped-rel')
            sage: a = R(2*17^2)
            sage: a.valuation()
            2
            sage: R = Zp(5, 4,'capped-rel')
            sage: R(0).valuation()
            +Infinity

        TESTS::

            sage: R(1).valuation()
            0
            sage: R(2).valuation()
            0
            sage: R(5).valuation()
            1
            sage: R(10).valuation()
            1
            sage: R(25).valuation()
            2
            sage: R(50).valuation()
            2
            sage: R = Qp(17, 4)
            sage: a = R(2*17^2)
            sage: a.valuation()
            2
            sage: R = Qp(5, 4)
            sage: R(0).valuation()
            +Infinity
            sage: R(1).valuation()
            0
            sage: R(2).valuation()
            0
            sage: R(5).valuation()
            1
            sage: R(10).valuation()
            1
            sage: R(25).valuation()
            2
            sage: R(50).valuation()
            2
            sage: R(1/2).valuation()
            0
            sage: R(1/5).valuation()
            -1
            sage: R(1/10).valuation()
            -1
            sage: R(1/25).valuation()
            -2
            sage: R(1/50).valuation()
            -2

            sage: K.<a> = Qq(25)                                                        # needs sage.libs.ntl
            sage: K(0).valuation()                                                      # needs sage.libs.ntl
            +Infinity

            sage: R(1/50).valuation(5)
            -2
            sage: R(1/50).valuation(3)
            Traceback (most recent call last):
            ...
            ValueError: Ring (5-adic Field with capped relative precision 4) residue field of the wrong characteristic.
        """
        if p is not None and p != self.parent().prime():
            raise ValueError('Ring (%s) residue field of the wrong characteristic.' % self.parent())
        cdef long v = self.valuation_c()
        if v == maxordp:
            return infinity
        if v == -maxordp:
            return -infinity
        cdef Integer ans = PY_NEW(Integer)
        mpz_set_si(ans.value, v)
        return ans

    cdef long valuation_c(self) noexcept:
        """
        This function is overridden in subclasses to provide an
        actual implementation of valuation.

        For exact zeros, ``maxordp`` is returned, rather than infinity.

        EXAMPLES:

        For example, the valuation function on pAdicCappedRelativeElements
        uses an overridden version of this function.

        ::

            sage: Zp(5)(5).valuation()  # indirect doctest
            1
        """
        raise NotImplementedError

    cpdef val_unit(self):
        r"""
        Return ``(self.valuation(), self.unit_part())``. To be overridden in
        derived classes.

        EXAMPLES::

            sage: Zp(5,5)(5).val_unit()
            (1, 1 + O(5^5))
        """
        raise NotImplementedError

    def ordp(self, p=None):
        r"""
        Return the valuation of ``self``, normalized so that the valuation of `p` is 1.

        INPUT:

        - ``self`` -- a `p`-adic element
        - ``p`` -- a prime (default: ``None``); if specified, will make sure
          that ``p == self.parent().prime()``

        .. NOTE::

            The optional argument `p` is used for consistency with the valuation
            methods on integers and rationals.

        OUTPUT:

        integer -- the valuation of ``self``, normalized so that the valuation of `p` is 1

        EXAMPLES::

            sage: R = Zp(5,20,'capped-rel')
            sage: R(0).ordp()
            +Infinity
            sage: R(1).ordp()
            0
            sage: R(2).ordp()
            0
            sage: R(5).ordp()
            1
            sage: R(10).ordp()
            1
            sage: R(25).ordp()
            2
            sage: R(50).ordp()
            2
            sage: R(1/2).ordp()
            0
        """
        return self.valuation(p) / self.parent().absolute_e()

    def is_prime(self):
        """
        Return whether this element is prime in its parent.

        EXAMPLES::

            sage: A = Zp(2)
            sage: A(1).is_prime()
            False
            sage: A(2).is_prime()
            True

            sage: K = A.fraction_field()
            sage: K(2).is_prime()
            False

        ::

            sage: # needs sage.libs.ntl
            sage: x = polygen(ZZ, 'x')
            sage: B.<pi> = A.extension(x^5 - 2)
            sage: pi.is_prime()                                                         # needs sage.symbolic
            True
            sage: B(2).is_prime()
            False
        """
        if self.is_zero():
            return True
        if self.parent().is_field():
            return False
        return self.valuation() == 1

    def rational_reconstruction(self):
        r"""
        Return a rational approximation to this `p`-adic number.

        This will raise an :exc:`ArithmeticError` if there are no valid
        approximations to the unit part with numerator and
        denominator bounded by ``sqrt(p^absprec / 2)``.

        .. SEEALSO::

            :meth:`_rational_`

        OUTPUT: rational; an approximation to ``self``

        EXAMPLES::

            sage: R = Zp(5,20,'capped-rel')
            sage: for i in range(11):
            ....:     for j in range(1,10):
            ....:         if j == 5:
            ....:             continue
            ....:         assert i/j == R(i/j).rational_reconstruction()
        """
        if self.is_zero(self.precision_absolute()):
            return Rational(0)
        p = self.parent().prime()
        alpha = self.unit_part().lift()
        m = Integer(p**self.precision_relative())
        from sage.arith.misc import rational_reconstruction
        r = rational_reconstruction(alpha, m)
        return (Rational(p)**self.valuation())*r

    def _rational_(self):
        r"""
        Return a rational approximation to this `p`-adic number.

        If there is no good rational approximation to the unit part,
        will just return the integer approximation.

        EXAMPLES::

            sage: R = Zp(7,5)
            sage: QQ(R(125))  # indirect doctest
            125
        """
        try:
            return self.rational_reconstruction()
        except ArithmeticError:
            p = self.parent().prime()
            return Rational(p**self.valuation() * self.unit_part().lift())

    def _number_field_(self, K):
        r"""
        Return an element of `K` approximating this `p`-adic number.

        INPUT:

        - ``K`` -- a number field

        EXAMPLES::

            sage: R.<a> = Zq(125)                                                       # needs sage.libs.ntl
            sage: K = R.exact_field()                                                   # needs sage.libs.ntl
            sage: a._number_field_(K)                                                   # needs sage.libs.ntl
            a
        """
        Kbase = K.base_ring()
        if K.defining_polynomial() != self.parent().defining_polynomial(exact=True):
            # Might convert to K's base ring.
            return Kbase(self)
        L = [Kbase(c) for c in self.polynomial().list()]
        if len(L) < K.relative_degree():
            L += [Kbase(0)] * (K.relative_degree() - len(L))
        return K(L)

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        Return the image of this element under the morphism defined by
        ``im_gens`` in ``codomain``, where elements of the
        base ring are mapped by ``base_map``.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<x> = ZZ[]
            sage: K.<a> = Qq(25, modulus=x^2-2)
            sage: L.<b> = Qq(625, modulus=x^4-2)
            sage: phi = K.hom([b^2]); phi(a+1)
            (b^2 + 1) + O(5^20)
            sage: z = L(-1).sqrt()
            sage: psi = L.hom([z*b]); psi(phi(a) + 5*b) == psi(phi(a)) + 5*psi(b)
            True
            sage: z = (1+5*b).log()
            sage: w = (5 - 5*b).exp()
            sage: psi(z*w) == psi(z) * psi(w)
            True

            sage: # needs sage.libs.ntl
            sage: P.<pi> = K.extension(x^2 - 5)
            sage: cc = K.hom([-a])
            sage: alpha = P.hom([pi], base_map=cc); alpha(a) + a
            O(pi^40)
            sage: zz = (1 + a*pi).log()
            sage: ww = pi.exp()                                                         # needs sage.symbolic
            sage: beta = P.hom([-pi], base_map=cc)
            sage: beta(ww*zz) == beta(ww)*beta(zz)                                      # needs sage.symbolic
            True
        """
        L = self.parent()
        K = L.base_ring()
        if L is K:
            # Qp or Zp, so there is a unique map
            if base_map is None:
                return codomain.coerce(self)
            else:
                return base_map(self)
        f = self.polynomial()
        if base_map is not None:
            f = f.change_ring(base_map)
        return f(im_gens[0])

    def _log_generic(self, aprec, mina=0):
        r"""
        Return ``\log(self)`` for ``self`` equal to 1 in the residue field.

        This is a helper method for :meth:`log`.

        INPUT:

        - ``aprec`` -- integer; the precision to which the result is
          correct. ``aprec`` must not exceed the precision cap of the ring over
          which this element is defined.

        - ``mina`` -- integer (default: 0); the series will check `n` up to
          this valuation (and beyond) to see if they can contribute to the
          series

        ALGORITHM:

        The computation uses the series

        .. MATH::

            -\log(1-x)=\sum_{n=1}^\infty \frac{x^n}{n}.

        For the result to be correct to precision ``aprec``, we sum all terms
        for which the valuation of `x^n/n` is strictly smaller than ``aprec``.

        EXAMPLES::

            sage: r = Qp(5,prec=4)(6)
            sage: r._log_generic(2)
            5 + O(5^2)
            sage: r._log_generic(4)
            5 + 2*5^2 + 4*5^3 + O(5^4)
            sage: r._log_generic(100)
            5 + 2*5^2 + 4*5^3 + O(5^4)

            sage: r = Zp(5,prec=4,type='fixed-mod')(6)
            sage: r._log_generic(5)
            5 + 2*5^2 + 4*5^3

        Only implemented for elements congruent to 1 modulo the maximal ideal::

            sage: r = Zp(5,prec=4,type='fixed-mod')(2)
            sage: r._log_generic(5)
            Traceback (most recent call last):
            ...
            ValueError: Input value (=2) must be 1 in the residue field
        """
        x = 1-self
        R = self.parent()
        # to get the precision right over capped-absolute rings, we have to
        # work over the capped-relative fraction field
        if R.is_capped_absolute():
            R = R.fraction_field()
            x = R(x)

        alpha=x.valuation()
        if alpha<=0:
            raise ValueError('Input value (=%s) must be 1 in the residue field' % self)

        e=R.absolute_e()
        p=R.prime()

        # we sum all terms of the power series of log into total
        total=R.zero()

        # pre-compute x^p/p into x2p_p
        if mina == 0 and alpha*p - e > aprec:
            # The value of x^p/p is not needed in that case
            x2p_p = R(0)
        elif R.is_capped_relative():
            if p*alpha >= e:
                x2p_p = x**p/p
            else:
                # x^p/p has negative valuation, so we need to be much
                # more careful about precision.
                x = x.lift_to_precision()
                x2p_p = x**p/p
        else:
            xu=x.unit_part()
            pu=R(p).unit_part()
            x2p_p=((xu**p)*pu.inverse_of_unit())*R.uniformizer_pow(p*alpha-e)

        # To get result right to precision aprec, we need all terms for which
        # the valuation of x^n/n is strictly smaller than aprec.
        # If we rewrite n=u*p^a with u a p-adic unit, then these are the terms
        # for which u<(aprec+a*v(p))/(v(x)*p^a).
        # Two sum over these terms, we run two nested loops, the outer one
        # iterates over the possible values for a, the inner one iterates over
        # the possible values for u.
        upper_u = (aprec/alpha).floor()
        if mina > 0 or upper_u > 0:
            a=0
            p2a=1       # p^a
            x2pa = x    # x^(p^a)

            # In the unramified case, we can stop summing terms as soon as
            # there are no u for a given a to sum over. In the ramified case,
            # it can happen that for some initial a there are no such u but
            # later in the series there are such u again. mina can be set to
            # take care of this by summing at least to a=mina-1
            while True:
                # we compute the sum for the possible values for u using Horner's method
                inner_sum = R.zero()
                for u in range(upper_u,0,-1):
                    # We want u to be a p-adic unit
                    if u%p==0:
                        new_term = R.zero()
                    else:
                        new_term = ~R(u)

                    # This hack is to deal with rings that don't lift to fields
                    if u > 1 or x2p_p.is_zero():
                        inner_sum = (inner_sum+new_term)*x2pa
                    else:
                        inner_sum = (inner_sum+new_term)*(x2p_p**a)*(x**(p2a-a*p))

                total -= inner_sum

                # Now increase a and check if a new iteration of the loop is needed
                a += 1
                p2a = p2a*p
                upper_u = ((aprec+a*e)/(alpha*p2a)).floor()
                if a >= mina and upper_u <= 0: break

                # We perform this last operation after the test
                # because it is costly and may raise OverflowError
                x2pa = x2pa**p

        return total.add_bigoh(aprec)

    def _log_binary_splitting(self, aprec, mina=0):
        r"""
        Return ``\log(self)`` for ``self`` equal to 1 in the residue field.

        This is a helper method for :meth:`log`.
        It uses a fast binary splitting algorithm.

        INPUT:

        - ``aprec`` -- integer; the precision to which the result is
          correct. ``aprec`` must not exceed the precision cap of the ring over
          which this element is defined.

        - ``mina`` -- integer (default: 0); the series will check `n` up to
          this valuation (and beyond) to see if they can contribute to the
          series

        .. NOTE::

            The function does not check that its argument ``self`` is
            1 in the residue field. If this assumption is not fulfilled
            the behaviour of the function is not specified.

        ALGORITHM:

        1. Raise `u` to the power `p^v` for a suitable `v` in order
           to make it closer to 1. (`v` is chosen such that `p^v` is
           close to the precision.)

        2. Write

        .. MATH::

            u^{p-1} = \prod_{i=1}^\infty (1 - a_i p^{(v+1)*2^i})

        with `0 \leq a_i < p^{(v+1)*2^i}` and compute
        `\log(1 - a_i p^{(v+1)*2^i})` using the standard Taylor expansion

        .. MATH::

            \log(1 - x) = -x - 1/2 x^2 - 1/3 x^3 - 1/4 x^4 - 1/5 x^5 - \cdots

        together with a binary splitting method.

        3. Divide the result by `p^v`

        The complexity of this algorithm is quasi-linear.

        EXAMPLES::

            sage: r = Qp(5,prec=4)(6)
            sage: r._log_binary_splitting(2)
            5 + O(5^2)
            sage: r._log_binary_splitting(4)
            5 + 2*5^2 + 4*5^3 + O(5^4)
            sage: r._log_binary_splitting(100)
            5 + 2*5^2 + 4*5^3 + O(5^4)

            sage: r = Zp(5,prec=4,type='fixed-mod')(6)
            sage: r._log_binary_splitting(5)
            5 + 2*5^2 + 4*5^3
        """
        raise NotImplementedError

    def log(self, p_branch=None, pi_branch=None, aprec=None, change_frac=False, algorithm=None):
        r"""
        Compute the `p`-adic logarithm of this element.

        The usual power series for the logarithm with values in the additive
        group of a `p`-adic ring only converges for 1-units (units congruent to
        1 modulo `p`).  However, there is a unique extension of the logarithm
        to a homomorphism defined on all the units: If `u = a \cdot v` is a
        unit with `v \equiv 1 \pmod{p}` and `a` a Teichmuller representative,
        then we define `log(u) = log(v)`.  This is the correct extension
        because the units `U` split as a product `U = V \times \langle w
        \rangle`, where `V` is the subgroup of 1-units and `w` is a fundamental
        root of unity.  The `\langle w \rangle` factor is torsion, so must go
        to 0 under any homomorphism to the fraction field, which is a torsion
        free group.

        INPUT:

        - ``p_branch`` -- an element in the base ring or its fraction
          field; the implementation will choose the branch of the
          logarithm which sends `p` to ``branch``

        - ``pi_branch`` -- an element in the base ring or its fraction
          field; the implementation will choose the branch of the
          logarithm which sends the uniformizer to ``branch``; you
          may specify at most one of ``p_branch`` and ``pi_branch``,
          and must specify one of them if this element is not a unit

        - ``aprec`` -- integer or ``None`` (default: ``None``); if not
          ``None``, then the result will only be correct to precision
          ``aprec``

        - ``change_frac`` -- in general the codomain of the logarithm should be
          in the `p`-adic field, however, for most neighborhoods of 1, it lies
          in the ring of integers. This flag decides if the codomain should be
          the same as the input (default) or if it should change to the
          fraction field of the input.

        - ``algorithm`` -- ``'generic'``, ``'binary_splitting'`` or ``None`` (default)
          The generic algorithm evaluates naively the series defining the log,
          namely

          .. MATH::

              \log(1-x) = -x - 1/2 x^2 - 1/3 x^3 - 1/4 x^4 - 1/5 x^5 - \cdots.

          Its binary complexity is quadratic with respect to the precision.

          The ``'binary_splitting'`` algorithm is faster, it has a quasi-linear
          complexity.
          By default, we use ``'binary_splitting'`` if it is available. Otherwise
          we switch to the ``'generic'`` algorithm.

        .. NOTE::

            What some other systems do:

            - PARI: Seems to define the logarithm for units not congruent
              to 1 as we do.

            - MAGMA: Only implements logarithm for 1-units (version 2.19-2)

        .. TODO::

            There is a soft-linear time algorithm for logarithm described
            by Dan Bernstein at
            http://cr.yp.to/lineartime/multapps-20041007.pdf

        EXAMPLES::

            sage: Z13 = Zp(13, 10)
            sage: a = Z13(14); a
            1 + 13 + O(13^10)
            sage: a.log()
            13 + 6*13^2 + 2*13^3 + 5*13^4 + 10*13^6 + 13^7 + 11*13^8 + 8*13^9 + O(13^10)

            sage: Q13 = Qp(13, 10)
            sage: a = Q13(14); a
            1 + 13 + O(13^10)
            sage: a.log()
            13 + 6*13^2 + 2*13^3 + 5*13^4 + 10*13^6 + 13^7 + 11*13^8 + 8*13^9 + O(13^10)

        Note that the relative precision decreases when we take log.
        Precisely the absolute precision on ``log(a)`` agrees with the relative
        precision on ``a`` thanks to the relation `d\log(a) = da/a`.

        The call ``log(a)`` works as well::

            sage: log(a)
            13 + 6*13^2 + 2*13^3 + 5*13^4 + 10*13^6 + 13^7 + 11*13^8 + 8*13^9 + O(13^10)
            sage: log(a) == a.log()
            True

        The logarithm is not only defined for 1-units::

            sage: R = Zp(5,10)
            sage: a = R(2)
            sage: a.log()
            2*5 + 3*5^2 + 2*5^3 + 4*5^4 + 2*5^6 + 2*5^7 + 4*5^8 + 2*5^9 + O(5^10)

        If you want to take the logarithm of a non-unit you must specify either
        ``p_branch`` or ``pi_branch``::

            sage: b = R(5)
            sage: b.log()
            Traceback (most recent call last):
            ...
            ValueError: you must specify a branch of the logarithm for non-units
            sage: b.log(p_branch=4)
            4 + O(5^10)
            sage: c = R(10)
            sage: c.log(p_branch=4)
            4 + 2*5 + 3*5^2 + 2*5^3 + 4*5^4 + 2*5^6 + 2*5^7 + 4*5^8 + 2*5^9 + O(5^10)

        The branch parameters are only relevant for elements of nonzero
        valuation::

            sage: a.log(p_branch=0)
            2*5 + 3*5^2 + 2*5^3 + 4*5^4 + 2*5^6 + 2*5^7 + 4*5^8 + 2*5^9 + O(5^10)
            sage: a.log(p_branch=1)
            2*5 + 3*5^2 + 2*5^3 + 4*5^4 + 2*5^6 + 2*5^7 + 4*5^8 + 2*5^9 + O(5^10)

        Logarithms can also be computed in extension fields. First, in an
        Eisenstein extension::

            sage: R = Zp(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^4 + 15*x^2 + 625*x - 5
            sage: W.<w> = R.ext(f)                                                      # needs sage.libs.ntl
            sage: z = 1 + w^2 + 4*w^7; z                                                # needs sage.libs.ntl
            1 + w^2 + 4*w^7 + O(w^20)
            sage: z.log()                                                               # needs sage.libs.ntl
            w^2 + 2*w^4 + 3*w^6 + 4*w^7 + w^9 + 4*w^10 + 4*w^11 + 4*w^12
             + 3*w^14 + w^15 + w^17 + 3*w^18 + 3*w^19 + O(w^20)

        In an extension, there will usually be a difference between
        specifying ``p_branch`` and ``pi_branch``::

            sage: # needs sage.libs.ntl
            sage: b = W(5)
            sage: b.log()
            Traceback (most recent call last):
            ...
            ValueError: you must specify a branch of the logarithm for non-units
            sage: b.log(p_branch=0)
            O(w^20)
            sage: b.log(p_branch=w)
            w + O(w^20)
            sage: b.log(pi_branch=0)
            3*w^2 + 2*w^4 + 2*w^6 + 3*w^8 + 4*w^10 + w^13 + w^14 + 2*w^15
             + 2*w^16 + w^18 + 4*w^19 + O(w^20)
            sage: b.unit_part().log()
            3*w^2 + 2*w^4 + 2*w^6 + 3*w^8 + 4*w^10 + w^13 + w^14 + 2*w^15
             + 2*w^16 + w^18 + 4*w^19 + O(w^20)
            sage: y = w^2 * 4*w^7; y
            4*w^9 + O(w^29)
            sage: y.log(p_branch=0)
            2*w^2 + 2*w^4 + 2*w^6 + 2*w^8 + w^10 + w^12 + 4*w^13 + 4*w^14 + 3*w^15
             + 4*w^16 + 4*w^17 + w^18 + 4*w^19 + O(w^20)
            sage: y.log(p_branch=w)
            w + 2*w^2 + 2*w^4 + 4*w^5 + 2*w^6 + 2*w^7 + 2*w^8 + 4*w^9 + w^10
             + 3*w^11 + w^12 + 4*w^14 + 4*w^16 + 2*w^17 + w^19 + O(w^20)

        Check that log is multiplicative::

            sage: y.log(p_branch=0) + z.log() - (y*z).log(p_branch=0)                   # needs sage.libs.ntl
            O(w^20)

        Now an unramified example::

            sage: # needs sage.libs.ntl
            sage: g = x^3 + 3*x + 3
            sage: A.<a> = R.ext(g)
            sage: b = 1 + 5*(1 + a^2) + 5^3*(3 + 2*a)
            sage: b.log()
            (a^2 + 1)*5 + (3*a^2 + 4*a + 2)*5^2 + (3*a^2 + 2*a)*5^3
             + (3*a^2 + 2*a + 2)*5^4 + O(5^5)

        Check that log is multiplicative::

            sage: # needs sage.libs.ntl
            sage: c = 3 + 5^2*(2 + 4*a)
            sage: b.log() + c.log() - (b*c).log()
            O(5^5)

        We illustrate the effect of the ``precision`` argument::

            sage: R = ZpCA(7,10)
            sage: x = R(41152263); x
            5 + 3*7^2 + 4*7^3 + 3*7^4 + 5*7^5 + 6*7^6 + 7^9 + O(7^10)
            sage: x.log(aprec = 5)
            7 + 3*7^2 + 4*7^3 + 3*7^4 + O(7^5)
            sage: x.log(aprec = 7)
            7 + 3*7^2 + 4*7^3 + 3*7^4 + 7^5 + 3*7^6 + O(7^7)
            sage: x.log()
            7 + 3*7^2 + 4*7^3 + 3*7^4 + 7^5 + 3*7^6 + 7^7 + 3*7^8 + 4*7^9 + O(7^10)

        The logarithm is not defined for zero::

            sage: R.zero().log()
            Traceback (most recent call last):
            ...
            ValueError: logarithm is not defined at zero

        For elements in a `p`-adic ring, the logarithm will be returned in the
        same ring::

            sage: x = R(2)
            sage: x.log().parent()
            7-adic Ring with capped absolute precision 10
            sage: x = R(14)
            sage: x.log(p_branch=0).parent()
            7-adic Ring with capped absolute precision 10

        This is not possible if the logarithm has negative valuation::

            sage: R = ZpCA(3,10)
            sage: S.<x> = R[]
            sage: f = x^3 - 3
            sage: W.<w> = R.ext(f)                                                      # needs sage.libs.ntl sage.rings.padics
            sage: w.log(p_branch=2)                                                     # needs sage.libs.ntl sage.rings.padics
            Traceback (most recent call last):
            ...
            ValueError: logarithm is not integral, use change_frac=True
            to obtain a result in the fraction field
            sage: w.log(p_branch=2, change_frac=True)                                   # needs sage.libs.ntl sage.rings.padics
            2*w^-3 + O(w^24)

        TESTS:

        Check that the generic algorithm and the binary splitting algorithm
        returns the same answers::

            sage: p = 17
            sage: R = Zp(p)
            sage: a = 1 + p*R.random_element()
            sage: l1 = a.log(algorithm='generic')
            sage: l2 = a.log(algorithm='binary_splitting')
            sage: l1 == l2
            True
            sage: l1.precision_absolute() == l2.precision_absolute()
            True

        Check multiplicativity::

            sage: p = 11
            sage: R = Zp(p, prec=1000)

            sage: x = 1 + p*R.random_element()
            sage: y = 1 + p*R.random_element()
            sage: log(x*y) == log(x) + log(y)
            True

            sage: x = y = 0
            sage: while x == 0:
            ....:     x = R.random_element()
            sage: while y == 0:
            ....:     y = R.random_element()
            sage: branch = R.random_element()
            sage: (x*y).log(p_branch=branch) == x.log(p_branch=branch) + y.log(p_branch=branch)
            True

        Note that multiplicativity may fail in the fixed modulus setting
        due to rounding errors::

            sage: R = ZpFM(2, prec=5)
            sage: R(180).log(p_branch=0) == R(30).log(p_branch=0) + R(6).log(p_branch=0)
            False

        Check that log is the inverse of exp::

            sage: p = 11
            sage: R = Zp(p, prec=1000)
            sage: a = 1 + p*R.random_element()
            sage: exp(log(a)) == a
            True

            sage: a = p*R.random_element()
            sage: log(exp(a)) == a
            True

        Check that results are consistent over a range of precision::

            sage: max_prec = 40
            sage: p = 3
            sage: K = Zp(p, max_prec)
            sage: full_log = (K(1 + p)).log()
            sage: for prec in range(2, max_prec):
            ....:     ll1 = (K(1+p).add_bigoh(prec)).log()
            ....:     ll2 = K(1+p).log(prec)
            ....:     assert ll1 == full_log
            ....:     assert ll2 == full_log
            ....:     assert ll1.precision_absolute() == prec

        Check that ``aprec`` works for fixed-mod elements::

            sage: R = ZpFM(7,10)
            sage: x = R(41152263); x
            5 + 3*7^2 + 4*7^3 + 3*7^4 + 5*7^5 + 6*7^6 + 7^9
            sage: x.log(aprec = 5)
            7 + 3*7^2 + 4*7^3 + 3*7^4
            sage: x.log(aprec = 7)
            7 + 3*7^2 + 4*7^3 + 3*7^4 + 7^5 + 3*7^6
            sage: x.log()                                                               # needs sage.symbolic
            7 + 3*7^2 + 4*7^3 + 3*7^4 + 7^5 + 3*7^6 + 7^7 + 3*7^8 + 4*7^9

        Check that precision is computed correctly in highly ramified
        extensions::

            sage: S.<x> = ZZ[]
            sage: K = Qp(5,5)
            sage: f = x^625 - 5*x - 5

            sage: # needs sage.libs.ntl sage.rings.padics
            sage: W.<w> = K.extension(f)
            sage: z = 1 - w^2 + O(w^11)
            sage: x = 1 - z
            sage: z.log().precision_absolute()
            -975
            sage: (x^5/5).precision_absolute()
            -570
            sage: (x^25/25).precision_absolute()
            -975
            sage: (x^125/125).precision_absolute()
            -775

            sage: # needs sage.libs.ntl
            sage: z = 1 - w + O(w^2)
            sage: x = 1 - z
            sage: z.log().precision_absolute()                                          # needs sage.rings.padics
            -1625
            sage: (x^5/5).precision_absolute()
            -615
            sage: (x^25/25).precision_absolute()
            -1200
            sage: (x^125/125).precision_absolute()
            -1625
            sage: (x^625/625).precision_absolute()
            -1250

            sage: z.log().precision_relative()                                          # needs sage.libs.ntl sage.rings.padics
            250

        Performances::

            sage: R = Zp(17, prec=10^5)
            sage: a = R.random_element()
            sage: b = a.log(p_branch=0)   # should be rather fast

        AUTHORS:

        - William Stein: initial version

        - David Harvey (2006-09-13): corrected subtle precision bug (need to
          take denominators into account! -- see :issue:`53`)

        - Genya Zaytman (2007-02-14): adapted to new `p`-adic class

        - Amnon Besser, Marc Masdeu (2012-02-21): complete rewrite, valid for
          generic `p`-adic rings.

        - Soroosh Yazdani (2013-02-1): Fixed a precision issue in
          :meth:`_log_generic`.  This should really fix the issue with
          divisions.

        - Julian Rueth (2013-02-14): Added doctests, some changes for
          capped-absolute implementations.

        - Xavier Caruso (2017-06): Added binary splitting type algorithms
          over Qp
        """
        if self.is_zero():
            raise ValueError('logarithm is not defined at zero')
        if p_branch is not None and pi_branch is not None:
            raise ValueError("you may only specify a branch of the logarithm in one way")
        R = self.parent()
        p = R.prime()
        q = p**R.absolute_f()
        e = R.absolute_e()

        if self.is_padic_unit():
            total = R.zero()
        else:
            if pi_branch is None:
                if p_branch is None:
                    raise ValueError("you must specify a branch of the logarithm for non-units")
                pi_branch = (p_branch - R._log_unit_part_p()) / e
                # Be careful: in ramified extensions, R._log_unit_part_p() is theoretically known at higher precision than the cap
                # In some cases, this may result in a loss of precision on pi_branch, and then on the final result
            total = self.valuation() * pi_branch
        y = self.unit_part()
        x = 1 - y

        if x.valuation()>0:
            denom=Integer(1)
        else:
            y=y**(q-1) # Is it better to multiply it by Teichmuller element?
            denom=Integer(q-1)
            x = 1 - y

        minaprec = y.precision_absolute()
        minn = 0
        if e != 1:
            xval = x.valuation()
            lamb = minaprec - xval
            if lamb > 0 and lamb*(p-1) <= e:
                # This is the precision region where the absolute
                # precision of the answer might be less than the
                # absolute precision of the input

                # kink is the number of times we multiply the relative
                # precision by p before starting to add e instead.
                kink = (e // (lamb * (p-1))).exact_log(p) + 1

                # deriv0 is within 1 of the n yielding the minimal
                # absolute precision
                tmp = (e / (minaprec * p.log(prec=53))).floor()
                if tmp > 0:
                    deriv0 = tmp.exact_log(p)
                else:
                    deriv0 = 0

                # These are the absolute precisions of x^(p^n) at potential minimum points
                L = [(minaprec * p**n - n * e, n) for n in [0, kink, deriv0, deriv0+1]]
                L.sort()
                minaprec = L[0][0]
                minn = L[0][1]

        if aprec is None or aprec > minaprec:
            aprec=minaprec

        if algorithm is None:
            try:
                # The binary splitting algorithm is supposed to be faster
                log_unit = y._log_binary_splitting(aprec, minn)
            except NotImplementedError:
                log_unit = y._log_generic(aprec, minn)
        elif algorithm == "generic":
            log_unit = y._log_generic(aprec, minn)
        elif algorithm == "binary_splitting":
            log_unit = y._log_binary_splitting(aprec, minn)
        else:
            raise ValueError("algorithm must be either 'generic', 'binary_splitting' or None")

        retval = total + log_unit*R(denom).inverse_of_unit()
        if not change_frac:
            if retval.valuation() < 0 and not R.is_field():
                raise ValueError("logarithm is not integral, use change_frac=True to obtain a result in the fraction field")
            retval = R(retval)
        return retval.add_bigoh(aprec)

    def _exp_generic(self, aprec):
        r"""
        Compute the exponential power series of this element, using Horner's
        evaluation and only one division.

        This is a helper method for :meth:`exp`.

        INPUT:

        - ``aprec`` -- integer; the precision to which to compute the
          exponential

        EXAMPLES::

            sage: R.<w> = Zq(7^2,5)                                                     # needs sage.libs.ntl
            sage: x = R(7*w)                                                            # needs sage.libs.ntl
            sage: x.exp(algorithm='generic')   # indirect doctest                       # needs sage.libs.ntl
            1 + w*7 + (4*w + 2)*7^2 + (w + 6)*7^3 + 5*7^4 + O(7^5)

        AUTHORS:

        - Genya Zaytman (2007-02-15)

        - Amnon Besser, Marc Masdeu (2012-02-23): Complete rewrite

        - Soroosh Yazdani (2013-02-01): Added the code for capped relative

        - Julian Rueth (2013-02-14): Rewrite to solve some precision problems
          in the capped-absolute case
        """
        R=self.parent()
        p=self.parent().prime()
        e=self.parent().absolute_e()
        x_unit=self.unit_part()
        p_unit=R(p).unit_part().lift_to_precision()
        x_val=self.valuation()

        # the valuation of n! is bounded by e*n/(p-1), therefore the valuation
        # of self^n/n! is bigger or equal to n*x_val - e*n/(p-1). So, we only
        # have to sum terms for which n does not exceed N
        N = (aprec // (x_val - e/(p-1))).floor()

        # We evaluate the exponential series:
        # First, we compute the value of x^N+N*x^(N-1)+...+x*N!+N! using
        # Horner's method. Then, we divide by N!.
        # This would only work for capped relative elements since for other
        # elements, we would lose too much precision in the multiplications
        # with natural numbers. Therefore, we emulate the behaviour of
        # capped-relative elements and keep track of the unit part and the
        # valuation separately.

        # the value of x^N + N*x^(N-1) + ... + (N-1)!*x + N!
        series_unit,series_val = R.one(), 0

        # we compute the value of N! as we go through the loop
        nfactorial_unit,nfactorial_val = R.one(),0

        nmodp = N%p
        for n in range(N,0,-1):
            # multiply everything by x
            series_val += x_val
            series_unit *= x_unit

            # compute the new value of N*(N-1)*...
            if nmodp == 0:
                n_pval, n_punit = Integer(n).val_unit(p)
                nfactorial_unit *= R(n_punit) * p_unit**n_pval
                nfactorial_val += n_pval*e
                nmodp = p
            else:
                nfactorial_unit *= n
            nmodp -= 1

            # now add N*(N-1)*...
            common_val = min(nfactorial_val, series_val)
            series_unit = (series_unit<<(series_val-common_val)) + (nfactorial_unit<<(nfactorial_val-common_val))
            series_val = common_val

        # multiply the result by N!
        return series_unit*nfactorial_unit.inverse_of_unit()<<(series_val-nfactorial_val)

    def _exp_binary_splitting(self, aprec):
        r"""
        Compute the exponential power series of this element.

        This is a helper method for :meth:`exp`.

        INPUT:

        - ``aprec`` -- integer; the precision to which to compute the
          exponential

        .. NOTE::

            The function does not check that its argument ``self`` is
            the disk of convergence of ``exp``. If this assumption is not
            fulfilled the behaviour of the function is not specified.

        ALGORITHM:

        Write

        .. MATH::

            self = \sum_{i=1}^\infty a_i p^{2^i}

        with `0 \leq a_i < p^{2^i}` and compute
        `\exp(a_i p^{2^i})` using the standard Taylor expansion

        .. MATH::

            \exp(x) = 1 + x + x^2/2 + x^3/6 + x^4/24 + \cdots

        together with a binary splitting method.

        The binary complexity of this algorithm is quasi-linear.

        EXAMPLES::

            sage: R = Zp(7,5)
            sage: x = R(7)
            sage: x.exp(algorithm='binary_splitting')   # indirect doctest
            1 + 7 + 4*7^2 + 2*7^3 + O(7^5)
        """
        raise NotImplementedError("the binary splitting algorithm is not implemented for the parent: %s" % self.parent())

    def _exp_newton(self, aprec, log_algorithm=None):
        r"""
        Compute the exponential power series of this element.

        This is a helper method for :meth:`exp`.

        INPUT:

        - ``aprec`` -- integer; the precision to which to compute the
          exponential

        - ``log_algorithm`` -- (default: ``None``) the algorithm used for
          computing the logarithm. This attribute is passed to the log
          method. See :meth:`log` for more details about the possible
          algorithms.

        .. NOTE::

            The function does not check that its argument ``self`` is
            the disk of convergence of ``exp``. If this assumption is not
            fulfilled the behaviour of the function is not specified.

        ALGORITHM:

        Solve the equation `\log(x) = self` using the Newton scheme::

        .. MATH::

            x_{i+1} = x_i \cdot (1 + self - \log(x_i))

        The binary complexity of this algorithm is roughly the same
        than that of the computation of the logarithm.

        EXAMPLES::

            sage: R.<w> = Zq(7^2,5)                                                     # needs sage.libs.ntl
            sage: x = R(7*w)                                                            # needs sage.libs.ntl
            sage: x.exp(algorithm='newton')   # indirect doctest                        # needs sage.libs.ntl
            1 + w*7 + (4*w + 2)*7^2 + (w + 6)*7^3 + 5*7^4 + O(7^5)
        """
        R = self.parent()
        e = R.absolute_e()
        a = R(1,aprec)
        l = R(0,aprec)
        if R.prime() == 2:
            trunc = e + 1
            while trunc < aprec:
                trunc = 2*trunc - e
                b = (self-l).add_bigoh(trunc).lift_to_precision(aprec)
                a *= 1+b
                l += (1+b).log(aprec, algorithm=log_algorithm)
        else:
            trunc = 1
            while trunc < aprec:
                trunc = 2*trunc
                b = (self-l).add_bigoh(trunc).lift_to_precision(aprec)
                a *= 1+b
                l += (1+b).log(aprec, algorithm=log_algorithm)
        return a

    def exp(self, aprec=None, algorithm=None):
        r"""
        Compute the `p`-adic exponential of this element if the exponential
        series converges.

        INPUT:

        - ``aprec`` -- integer or ``None`` (default: ``None``); if
          specified, computes only up to the indicated precision

        - ``algorithm`` -- ``'generic'``, ``'binary_splitting'``, ``'newton'``
          or ``None`` (default)

          The ``'generic'`` algorithm evaluates naively the series defining the
          exponential, namely

          .. MATH::

              \exp(x) = 1 + x + x^2/2 + x^3/6 + x^4/24 + \cdots.

          Its binary complexity is quadratic with respect to the precision.

          The ``'binary_splitting'`` algorithm is faster, it has a quasi-linear
          complexity.

          The ``'newton'`` algorithms solve the equation `\log(x) =` ``self``
          using a Newton scheme. It runs roughly as fast as the computation
          of the logarithm.

          By default, we use the ``'binary_splitting'`` if it is available.
          If it is not, we use the ``'newton'`` algorithm if a fast algorithm for
          computing the logarithm is available.
          Otherwise we switch to the ``'generic'`` algorithm.

        EXAMPLES:

        :meth:`log` and :meth:`exp` are inverse to each other::

            sage: Z13 = Zp(13, 10)
            sage: a = Z13(14); a
            1 + 13 + O(13^10)
            sage: a.log().exp()
            1 + 13 + O(13^10)

        An error occurs if this is called with an element for which the
        exponential series does not converge::

            sage: Z13.one().exp()
            Traceback (most recent call last):
            ...
            ValueError: Exponential does not converge for that input.

        The next few examples illustrate precision when computing `p`-adic
        exponentials::

            sage: R = Zp(5,10)
            sage: e = R(2*5 + 2*5**2 + 4*5**3 + 3*5**4
            ....:        + 5**5 + 3*5**7 + 2*5**8 + 4*5**9).add_bigoh(10); e
            2*5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^5 + 3*5^7 + 2*5^8 + 4*5^9 + O(5^10)
            sage: e.exp()*R.teichmuller(4)
            4 + 2*5 + 3*5^3 + O(5^10)

        ::

            sage: K = Qp(5,10)
            sage: e = K(2*5 + 2*5**2 + 4*5**3 + 3*5**4
            ....:        + 5**5 + 3*5**7 + 2*5**8 + 4*5**9).add_bigoh(10); e
            2*5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^5 + 3*5^7 + 2*5^8 + 4*5^9 + O(5^10)
            sage: e.exp()*K.teichmuller(4)
            4 + 2*5 + 3*5^3 + O(5^10)

        Logarithms and exponentials in extension fields. First, in an
        Eisenstein extension::

            sage: # needs sage.libs.ntl
            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^4 + 15*x^2 + 625*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = 1 + w^2 + 4*w^7; z
            1 + w^2 + 4*w^7 + O(w^20)
            sage: z.log().exp()
            1 + w^2 + 4*w^7 + O(w^20)

        Now an unramified example::

            sage: # needs sage.libs.ntl
            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: g = x^3 + 3*x + 3
            sage: A.<a> = R.ext(g)
            sage: b = 1 + 5*(1 + a^2) + 5^3*(3 + 2*a); b
            1 + (a^2 + 1)*5 + (2*a + 3)*5^3 + O(5^5)
            sage: b.log().exp()
            1 + (a^2 + 1)*5 + (2*a + 3)*5^3 + O(5^5)

        TESTS:

        Check that results are consistent over a range of precision::

            sage: max_prec = 40
            sage: p = 3
            sage: K = Zp(p, max_prec)
            sage: full_exp = (K(p)).exp()
            sage: for prec in range(2, max_prec):
            ....:     ll = (K(p).add_bigoh(prec)).exp()
            ....:     assert ll == full_exp
            ....:     assert ll.precision_absolute() == prec
            sage: K = Qp(p, max_prec)
            sage: full_exp = (K(p)).exp()
            sage: for prec in range(2, max_prec):
            ....:     ll = (K(p).add_bigoh(prec)).exp()
            ....:     assert ll == full_exp
            ....:     assert ll.precision_absolute() == prec

        Check that this also works for capped-absolute implementations::

            sage: Z13 = ZpCA(13, 10)
            sage: a = Z13(14); a
            1 + 13 + O(13^10)
            sage: a.log().exp()
            1 + 13 + O(13^10)

            sage: # needs sage.libs.ntl
            sage: R = ZpCA(5,5)
            sage: S.<x> = R[]
            sage: f = x^4 + 15*x^2 + 625*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = 1 + w^2 + 4*w^7; z
            1 + w^2 + 4*w^7 + O(w^20)
            sage: z.log().exp()
            1 + w^2 + 4*w^7 + O(w^20)

        Check that this also works for fixed-mod implementations::

            sage: Z13 = ZpFM(13, 10)
            sage: a = Z13(14); a
            1 + 13
            sage: a.log().exp()
            1 + 13

            sage: # needs sage.libs.ntl
            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^4 + 15*x^2 + 625*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = 1 + w^2 + 4*w^7; z
            1 + w^2 + 4*w^7
            sage: z.log().exp()
            1 + w^2 + 4*w^7

        Some corner cases::

            sage: Z2 = Zp(2, 5)
            sage: Z2(2).exp()
            Traceback (most recent call last):
            ...
            ValueError: Exponential does not converge for that input.

            sage: # needs sage.libs.ntl
            sage: S.<x> = Z2[]
            sage: W.<w> = Z2.ext(x^3-2)
            sage: (w^2).exp()
            Traceback (most recent call last):
            ...
            ValueError: Exponential does not converge for that input.
            sage: (w^3).exp()
            Traceback (most recent call last):
            ...
            ValueError: Exponential does not converge for that input.
            sage: (w^4).exp()
            1 + w^4 + w^5 + w^7 + w^9 + w^10 + w^14 + O(w^15)

        Check that all algorithms output the same result::

            sage: R = Zp(5,50)
            sage: a = 5 * R.random_element()
            sage: bg = a.exp(algorithm='generic')
            sage: bbs = a.exp(algorithm='binary_splitting')
            sage: bn = a.exp(algorithm='newton')
            sage: bg == bbs
            True
            sage: bg == bn
            True

        Performances::

            sage: R = Zp(17,10^5)
            sage: a = 17 * R.random_element()
            sage: b = a.exp()    # should be rather fast

        AUTHORS:

        - Genya Zaytman (2007-02-15)

        - Amnon Besser, Marc Masdeu (2012-02-23): Complete rewrite

        - Julian Rueth (2013-02-14): Added doctests, fixed some corner cases

        - Xavier Caruso (2017-06): Added binary splitting and Newton algorithms
        """
        p = self.parent().prime()

        if (p-1)*self.valuation() <= self.parent().absolute_e():
            raise ValueError('Exponential does not converge for that input.')

        # The optimal absolute precision on exp(self)
        # is the absolution precision on self
        maxprec = min(self.precision_absolute(), self.parent().precision_cap())
        if aprec is None or aprec > maxprec:
            aprec = maxprec

        if algorithm is None:
            try:
                ans = self._exp_binary_splitting(aprec)
            except NotImplementedError:
                try:
                    ans = self._exp_newton(aprec, log_algorithm='binary_splitting')
                except NotImplementedError:
                    ans = self._exp_generic(aprec)
        elif algorithm == 'generic':
            ans = self._exp_generic(aprec)
        elif algorithm == 'binary_splitting':
            ans = self._exp_binary_splitting(aprec)
        elif algorithm == 'newton':
            ans = self._exp_newton(aprec)
        else:
            raise ValueError("algorithm must be 'generic', 'binary_splitting', 'newton' or None")
        return ans.add_bigoh(aprec)

    def square_root(self, extend=True, all=False, algorithm=None):
        r"""
        Return the square root of this `p`-adic number.

        INPUT:

        - ``self`` -- a `p`-adic element

        - ``extend`` -- boolean (default: ``True``); if ``True``, return a
          square root in an extension if necessary; if ``False`` and no root
          exists in the given ring or field, raise a :exc:`ValueError`.

        - ``all`` -- boolean (default: ``False``); if ``True``, return a
          list of all square roots

        - ``algorithm`` -- ``'pari'``, ``'sage'`` or ``None`` (default:
          ``None``); Sage provides an implementation for any extension of
          `\QQ_p`, whereas only square roots over `\QQ_p` are implemented in
          PARI. The default is ``'pari'`` if the ground field is `\QQ_p`,
          ``'sage'`` otherwise.

        OUTPUT:

        The square root or the list of all square roots of this `p`-adic
        number.

        NOTE:

        The square root is chosen (resp. the square roots are ordered) in
        a deterministic way.

        EXAMPLES::

            sage: R = Zp(3, 20)
            sage: R(0).square_root()
            0

            sage: R(1).square_root()
            1 + O(3^20)

            sage: R(2).square_root(extend=False)
            Traceback (most recent call last):
            ...
            ValueError: element is not a square

            sage: -R(4).square_root()
            2 + O(3^20)

            sage: R(9).square_root()
            3 + O(3^21)

        When `p = 2`, the precision of the square root is less
        than the input::

            sage: R2 = Zp(2, 20)
            sage: R2(1).square_root()
            1 + O(2^19)
            sage: R2(4).square_root()
            2 + O(2^20)

            sage: R.<t> = Zq(2^10, 10)                                                  # needs sage.libs.ntl
            sage: u = 1 + 8*t                                                           # needs sage.libs.ntl
            sage: u.square_root()                                                       # needs sage.libs.ntl
            1 + t*2^2 + t^2*2^3 + t^2*2^4 + (t^4 + t^3 + t^2)*2^5 + (t^4 + t^2)*2^6
             + (t^5 + t^2)*2^7 + (t^6 + t^5 + t^4 + t^2)*2^8 + O(2^9)

            sage: # needs sage.libs.ntl
            sage: x = polygen(ZZ, 'x')
            sage: R.<a> = Zp(2).extension(x^3 - 2)
            sage: u = R(1 + a^4 + a^5 + a^7 + a^8, 10); u
            1 + a^4 + a^5 + a^7 + a^8 + O(a^10)
            sage: v = u.square_root(); v
            1 + a^2 + a^4 + a^6 + O(a^7)

        However, observe that the precision increases to its original value
        when we recompute the square of the square root::

            sage: v^2                                                                   # needs sage.libs.ntl
            1 + a^4 + a^5 + a^7 + a^8 + O(a^10)

        If the input does not have enough precision in order to determine if
        the given element has a square root in the ground field, an error is
        raised::

            sage: # needs sage.libs.ntl
            sage: R(1, 6).square_root()
            Traceback (most recent call last):
            ...
            PrecisionError: not enough precision to be sure that this element has a square root
            sage: R(1, 7).square_root()
            1 + O(a^4)
            sage: R(1+a^6, 7).square_root(extend=False)
            Traceback (most recent call last):
            ...
            ValueError: element is not a square

        In particular, an error is raised when we try to compute the square
        root of an inexact zero.

        TESTS::

            sage: R = Qp(5, 100)
            sage: c = R.random_element()
            sage: s = (c^2).square_root()
            sage: s == c or s == -c
            True

            sage: Q2 = Qp(2,20,'capped-rel')
            sage: Q2(1).square_root()
            1 + O(2^19)
            sage: Q2(4).square_root()
            2 + O(2^20)

            sage: Q5 = Qp(5,20,'capped-rel')
            sage: Q5(1).square_root()
            1 + O(5^20)
            sage: Q5(-1).square_root() == Q5.teichmuller(2) or Q5(-1).square_root() == Q5.teichmuller(3)
            True

            sage: Z3 = Zp(3,20,'capped-abs')
            sage: Z3(1).square_root()
            1 + O(3^20)
            sage: Z3(4).square_root() in [ Z3(2), -Z3(2) ]
            True
            sage: Z3(9).square_root()
            3 + O(3^19)

            sage: Z2 = Zp(2,20,'capped-abs')
            sage: Z2(1).square_root()
            1 + O(2^19)
            sage: Z2(4).square_root()
            2 + O(2^18)
            sage: Z2(9).square_root() in [ Z2(3), -Z2(3) ]
            True
            sage: Z2(17).square_root()
            1 + 2^3 + 2^5 + 2^6 + 2^7 + 2^9 + 2^10 + 2^13 + 2^16 + 2^17 + O(2^19)

            sage: Z5 = Zp(5,20,'capped-abs')
            sage: Z5(1).square_root()
            1 + O(5^20)
            sage: Z5(-1).square_root() == Z5.teichmuller(2) or Z5(-1).square_root() == Z5.teichmuller(3)
            True
        """
        # We first check trivial cases and precision
        if self._is_exact_zero():
            return self
        parent = self.parent()
        if self.is_zero() or (parent.prime() == 2 and self.precision_relative() < 1 + 2*parent.absolute_e()):
            raise PrecisionError("not enough precision to be sure that this element has a square root")

        if algorithm is None:
            if parent.absolute_degree() == 1:
                algorithm = "pari"
            else:
                algorithm = "sage"

        ans = None
        if algorithm == "pari":
            from cypari2.handle_error import PariError
            try:
                # use pari
                ans = parent(self.__pari__().sqrt())
            except PariError:
                # todo: should eventually change to return an element of
                # an extension field
                pass
        elif algorithm == "sage":
            try:
                ans = self.nth_root(2)
            except ValueError:
                pass
        if ans is not None:
            ans2 = -ans
            E1 = ans.expansion()
            E2 = ans2.expansion()
            if ans.parent().is_field():
                i = 0
            else:
                i = ans.valuation()
            while True:
                try:
                    d1 = E1[i]
                    d2 = E2[i]
                except (PrecisionError, IndexError):
                    break
                if d1 > d2:
                    ans, ans2 = ans2, ans
                    break
                if d1 < d2:
                    break
                i += 1
            if all:
                return [ans, ans2]
            else:
                return ans
        if extend:
            raise NotImplementedError("extending using the sqrt function not yet implemented")
        elif all:
            return []
        else:
            raise ValueError("element is not a square")

    def nth_root(self, n, all=False):
        """
        Return the `n`-th root of this element.

        INPUT:

        - ``n`` -- integer

        - ``all`` -- boolean (default: ``False``); if ``True``,
          return all `n`-th roots of this element, instead of just one

        EXAMPLES::

            sage: A = Zp(5,10)
            sage: x = A(61376); x
            1 + 5^3 + 3*5^4 + 4*5^5 + 3*5^6 + O(5^10)
            sage: y = x.nth_root(4); y
            2 + 5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^6 + O(5^10)
            sage: y^4 == x
            True

            sage: x.nth_root(4, all=True)
            [2 + 5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^6 + O(5^10),
             4 + 4*5 + 4*5^2 + 4*5^4 + 3*5^5 + 5^6 + 3*5^7 + 5^8 + 5^9 + O(5^10),
             3 + 3*5 + 2*5^2 + 5^4 + 4*5^5 + 3*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + O(5^10),
             1 + 4*5^3 + 5^5 + 3*5^6 + 5^7 + 3*5^8 + 3*5^9 + O(5^10)]

        When `n` is divisible by the underlying prime `p`, we
        are losing precision (which is consistent with the fact
        that raising to the `p`-th power increases precision)::

            sage: z = x.nth_root(5); z
            1 + 5^2 + 3*5^3 + 2*5^4 + 5^5 + 3*5^7 + 2*5^8 + O(5^9)
            sage: z^5
            1 + 5^3 + 3*5^4 + 4*5^5 + 3*5^6 + O(5^10)

        Everything works over extensions as well::

            sage: # needs sage.libs.ntl
            sage: W.<a> = Zq(5^3)
            sage: S.<x> = W[]
            sage: R.<pi> = W.extension(x^7 - 5)
            sage: R(5).nth_root(7)
            pi + O(pi^141)
            sage: R(5).nth_root(7, all=True)
            [pi + O(pi^141)]

        An error is raised if the given element is not an `n`-th power
        in the ring::

            sage: R(5).nth_root(11)                                                     # needs sage.libs.ntl
            Traceback (most recent call last):
            ...
            ValueError: this element is not a nth power

        Similarly, when precision on the input is too small, an error
        is raised::

            sage: # needs sage.libs.ntl
            sage: x = R(1,6); x
            1 + O(pi^6)
            sage: x.nth_root(5)
            Traceback (most recent call last):
            ...
            PrecisionError: not enough precision to be sure that this element is a nth power

        Check that :issue:`30314` is fixed::

            sage: # needs sage.libs.ntl
            sage: K = Qp(29)
            sage: x = polygen(K)
            sage: L.<a> = K.extension(x^2 - 29)
            sage: L(4).nth_root(2)
            2 + O(a^40)

        TESTS:

        We check that it works over different fields::

            sage: # needs sage.libs.ntl
            sage: K.<a> = Qq(2^3)
            sage: S.<x> = K[]
            sage: L.<pi> = K.extension(x^2 + 2*x + 2)
            sage: elt = L.random_element()
            sage: elt in (elt^8).nth_root(8, all=True)
            True
            sage: elt = L.random_element()
            sage: elt in (elt^16).nth_root(16, all=True)
            True
            sage: elt = L.random_element()
            sage: elt in (elt^56).nth_root(56, all=True)
            True

            sage: # needs sage.libs.ntl
            sage: K.<a> = Qq(3^2)
            sage: S.<x> = K[]
            sage: Z = (1+x)^3
            sage: E = Z^2 + Z + 1
            sage: L.<pi> = K.extension(E)
            sage: elt = L.random_element()
            sage: elt in (elt^9).nth_root(9, all=True)
            True
            sage: elt = L.random_element()
            sage: elt in (elt^27).nth_root(27, all=True)
            True
            sage: elt = L.random_element()
            sage: elt in (elt^108).nth_root(108, all=True)
            True

            sage: # needs sage.libs.flint sage.libs.ntl
            sage: K.<a> = ZqCA(3^2)
            sage: S.<x> = K[]
            sage: Z = (1+x)^3 + 3*x^2
            sage: E = Z^2 + Z + 1
            sage: L.<pi> = K.extension(E)
            sage: elt = L.random_element()
            sage: elt in (elt^9).nth_root(9, all=True)
            True
            sage: elt = L.random_element()
            sage: try:                                                                  # needs sage.rings.real_double
            ....:     assert elt in (elt^27).nth_root(27, all=True)
            ....: except sage.rings.padics.precision_error.PrecisionError:
            ....:     pass
            sage: elt = L.random_element()
            sage: try:                                                                  # needs sage.rings.real_double
            ....:     assert elt in (elt^108).nth_root(108, all=True)
            ....: except sage.rings.padics.precision_error.PrecisionError:
            ....:     pass

            sage: # needs sage.libs.ntl
            sage: K.<a> = Qq(3^2)
            sage: S.<x> = K[]
            sage: Z = (1+x)^3 + 3*x^3
            sage: E = (Z^2 + Z + 1)(a*x).monic()
            sage: L.<pi> = K.extension(E)
            sage: elt = L.random_element()
            sage: elt in (elt^9).nth_root(9, all=True)
            True
            sage: elt = L.random_element()
            sage: elt in (elt^27).nth_root(27, all=True)
            True
            sage: elt = L.random_element()
            sage: elt in (elt^108).nth_root(108, all=True)
            True
        """
        n = ZZ(n)
        if n == 0:
            raise ValueError("n must be a nonzero integer")
        elif n == 1:
            return self
        elif n < 0:
            return (~self).nth_root(-n)
        parent = self.parent()
        K = parent.fraction_field()  # due to conversion issues
        p = parent.prime()
        e = parent.absolute_e()
        ep = e // (p-1)

        # We first check trivial cases
        if self._is_exact_zero():
            return self
        if self.is_zero():
            raise PrecisionError("not enough precision to be sure that this element is a nth power")

        v = n.valuation(p)
        m = n // (p**v)

        # We check the valuation
        val = self.valuation()
        if val % n != 0:
            raise ValueError("this element is not a nth power")
        # and the residue
        a = K(self) >> val
        abar = a.residue()
        try:
            xbar = abar.nth_root(m)
        except ValueError:
            raise ValueError("this element is not a nth power")

        # We take the inverse mth root at small precision
        prec = a.precision_absolute()
        minprec = v*e + ep + 1
        if m == 1:
            parity = 0
            root = a.add_bigoh(minprec)
        else:
            parity = 1
            root = K(~xbar)
            invm = K(1/m)
            curprec = 1
            while curprec < min(minprec,prec):
                curprec *= 2
                root = root.lift_to_precision(min(minprec,prec,curprec))
                root += invm * root * (1 - a*(root**m))

        # We now extract the `(p^v)`-th root
        zeta, s, nextzeta = K._primitive_qth_root_of_unity(v)
        if v:
            nextzeta = (parent(nextzeta[0]), nextzeta[1])  # nextzeta[0] may have a wrong parent (with more precision)
        for i in range(v):
            if s > 0 and i >= s:
                root, accuracy = root._inverse_pth_root(twist=zeta, hint=nextzeta)
            else:
                root, accuracy = root._inverse_pth_root()
            if accuracy is not infinity and accuracy is not None:
                raise ValueError("this element is not a nth power")

        # We check the precision
        if v > 0 and prec < minprec:
            raise PrecisionError("not enough precision to be sure that this element is a nth power")

        # We lift the root using Newton iteration
        if v % 2 == parity:
            root = ~root
        invn = K(1/n)
        curprec = minprec
        while curprec < prec:
            curprec -= v*e
            curprec = min(2*curprec + v*e, p*curprec + (v-1)*e)
            root = root.lift_to_precision(min(prec,curprec))
            root += invn * root * (1 - a*(root**n))
        root = (~root) << (val // n)

        if all:
            return [ parent(root*zeta) for zeta in K.roots_of_unity(n) ]
        else:
            return parent(root)

    def _inverse_pth_root(self, twist=None, hint=None):
        r"""
        In its simplest form, computes the inverse of
        ``p``-th root of this element.

        This is a helper function used in :meth:`nth_root`
        and :meth:`primitive_root_of_unity`.

        INPUT:

        - ``twist`` -- an element in the same parent or ``None``
          (default: ``None``)

        - ``hint`` -- tuple or ``None`` (default: ``None``); if not
          ``None``, it has to be the output of ``twist._inverse_pth_root()``

        OUTPUT:

        When ``twist`` is ``None``, the output is a couple
        ``(invroot, accuracy)`` where:

        - ``accuracy`` is the highest valuation of an element of
          the form ``self * x^p - 1`` for `x` varying in the
          parent of this element, and

        - ``invroot`` is an element realizing this maximum.

        If the precision on the element is not enough to determine
        ``accuracy``, the value ``None`` is returned.

        When ``twist`` is not ``None``, the maximum is taken over
        all elements of the form ``self * x^p * twist^i - 1`` for
        for `x` varying in the parent of this element and `i`
        varying in the range `\{0, 1, \ldots, p-1\}`

        .. NOTE::

            This function assumes that the input element and ``twist``
            (if given) are units in the integer ring.

        TESTS::

            sage: R = Zp(11)
            sage: [ R.teichmuller(x).nth_root(11) == R.teichmuller(x) for x in range(1,11) ]  # indirect doctest
            [True, True, True, True, True, True, True, True, True, True]

            sage: # needs sage.libs.ntl
            sage: W.<a> = Zq(5^3)
            sage: S.<x> = W[]
            sage: R.<pi> = W.extension(x^8 + 15*a*x - 5)
            sage: y = R.random_element()
            sage: for n in [5, 10, 15]:                                                       # indirect doctest
            ....:     z = y**n
            ....:     assert z.nth_root(n)**n == z
        """
        ring = self.parent()
        p = ring.prime()
        e = ring.absolute_e()
        ep = e // (p-1)

        if twist is None:
            accuracy = None
        else:
            if hint is None:
                invroottwist, accuracy = twist._inverse_pth_root()
            else:
                invroottwist, accuracy = hint
            if accuracy is None:
                raise NotImplementedError("try to increase the precision cap of the parent...")

        a = self
        prec = a.precision_absolute()

        # We will need 1/a at higher precision
        ainv = ~(a.add_bigoh(e+ep+1))

        # We lift modulo pi^(e // (p-1))
        k = ring.residue_field()
        x = ring(0)
        curprec = 0  # curprec is the valuation of (ainv - x^p)
        while curprec < min(prec, e+ep):
            # recomputing x^p is not necessary:
            # we can alternatively update it after each update of x
            # (which is theoretically a bit faster)
            b = ainv - x**p
            if b == 0: break
            curprec = b.valuation()
            bexp = iter(b.unit_part().expansion())
            maxprec = prec
            while curprec < maxprec:
                try:
                    coeff = k(next(bexp))
                except StopIteration:
                    coeff = k(0)
                if coeff != 0:
                    if curprec % p == 0:
                        cbar = coeff.nth_root(p)
                        c = ring(cbar).lift_to_precision()
                        exponent = curprec // p
                        x += c << exponent
                        maxprec = min(maxprec, exponent + e)
                    elif accuracy == curprec:
                        alpha = (twist * invroottwist.add_bigoh(1 + curprec // p)**p - 1) >> curprec
                        exponent = coeff / (ainv.residue() * alpha.residue())
                        try:
                            exponent = ZZ(exponent)
                        except TypeError:
                            return x, curprec
                        else:
                            ainv //= twist**exponent
                            a *= twist**exponent
                            x *= invroottwist**exponent
                            break
                    else:
                        return x, curprec
                curprec += 1

        # We check if the precision was enough
        # We didn't do it before because we could have proved
        # that there is no pth root in the previous step.
        if prec < e + ep + 1:
            x = x.add_bigoh((prec+p-1) // p)
            return x, None

        # We lift one step further
        curprec = e + ep
        if e % (p-1) == 0:
            rho = k(ring(p).expansion(e))
            b = ainv - x**p
            b >>= curprec
            coeff = -a.residue()*b.residue()
            if accuracy == curprec:
                sigma = rho * (-rho).nth_root(p-1)  # should never fail
                alpha = (twist * invroottwist.add_bigoh(ep+1)**p - 1) >> curprec
                alpha = alpha.residue()
                tr = (alpha/sigma).trace()
                if tr != 0:
                    exponent = ZZ(-(coeff/sigma).trace() / tr)
                    coeff += exponent*alpha
                    ainv //= twist**exponent
                    a *= twist**exponent
                    x *= invroottwist**exponent
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            S = PolynomialRing(k, name='x')
            AS = S([ coeff, rho ] + (p-2)*[0] + [1])
            roots = AS.roots()
            if len(roots) == 0:
                return x, curprec
            x += ring(roots[0][0] * x.residue()) << ep

        # We perform Newton iteration
        curprec += 1
        while curprec < prec:
            curprec -= e
            curprec = min(2*curprec + e, p*curprec)
            x = x.lift_to_precision(min(prec,curprec))
            x += x * (1 - a*x**p) / p

        return x, infinity

    def __abs__(self):
        """
        Return the `p`-adic absolute value of ``self``.

        This is normalized so that the absolute value of `p` is `1/p`.

        EXAMPLES::

            sage: abs(Qp(5)(15))
            1/5
            sage: abs(Qp(7)(0))
            0

        An unramified extension::

            sage: # needs sage.libs.ntl
            sage: R = Zp(5,5)
            sage: P.<x> = PolynomialRing(R)
            sage: Z25.<u> = R.ext(x^2 - 3)
            sage: abs(u)
            1
            sage: abs(u^24-1)
            1/5

        A ramified extension::

            sage: # needs sage.libs.ntl
            sage: W.<w> = R.ext(x^5 + 75*x^3 - 15*x^2 + 125*x - 5)
            sage: abs(w)
            0.724779663677696
            sage: abs(W(0))
            0.000000000000000
        """
        return self.abs()

    cpdef abs(self, prec=None):
        """
        Return the `p`-adic absolute value of ``self``.

        This is normalized so that the absolute value of `p` is `1/p`.

        INPUT:

        - ``prec`` -- integer; the precision of the real field in which
          the answer is returned.  If ``None``, returns a rational for
          absolutely unramified fields, or a real with 53 bits of
          precision for ramified fields.

        EXAMPLES::

            sage: a = Qp(5)(15); a.abs()
            1/5
            sage: a.abs(53)                                                             # needs sage.rings.real_mpfr
            0.200000000000000
            sage: Qp(7)(0).abs()
            0
            sage: Qp(7)(0).abs(prec=20)                                                 # needs sage.rings.real_mpfr
            0.00000

        An unramified extension::

            sage: # needs sage.libs.ntl
            sage: R = Zp(5,5)
            sage: P.<x> = PolynomialRing(R)
            sage: Z25.<u> = R.ext(x^2 - 3)
            sage: u.abs()
            1
            sage: (u^24-1).abs()
            1/5

        A ramified extension::

            sage: # needs sage.libs.ntl
            sage: W.<w> = R.ext(x^5 + 75*x^3 - 15*x^2 + 125*x - 5)
            sage: w.abs()
            0.724779663677696
            sage: W(0).abs()
            0.000000000000000
        """
        K = self.parent()
        if not prec and K.absolute_e() > 1:
            prec = 53
        if prec:
            from sage.rings.real_mpfr import RealField
            if self.is_zero():
                return RealField(prec).zero()
            return RealField(prec)(K.prime())**(-self.ordp())
        else:
            if self.is_zero():
                return Rational(0)
            return Rational(K.prime())**(-self.valuation())

    cpdef bint _is_base_elt(self, p) except -1:
        r"""
        Return ``True`` if this element is an element of `\ZZ_p` or `\QQ_p` (rather than
        an extension).

        INPUT:

        - ``p`` -- a prime, which is compared with the parent of this element

        EXAMPLES::

            sage: a = Zp(5)(3); a._is_base_elt(5)
            True
            sage: a._is_base_elt(17)
            False
        """
        raise NotImplementedError

    def _polylog_res_1(self, n, p_branch=0):
        r"""
        Return `Li_n(`self`)`, the `n`-th `p`-adic polylogarithm of ``self``,
        assuming that ``self`` is congruent to `1 \pmod p`.

        This is an internal function, used by :meth:`polylog`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: `Li_n(`self`)`

        EXAMPLES::

            sage: Qp(2)(-1)._polylog_res_1(6) == 0                                      # needs sage.symbolic
            True

            sage: Qp(5)(1)._polylog_res_1(1)
            Traceback (most recent call last):
            ...
            ValueError: Polylogarithm is not defined for 1.

        Only polylogarithms for `n` at least two are defined by this function ::

            sage: Qp(11)(2)._polylog_res_1(1)
            Traceback (most recent call last):
            ...
            ValueError: Polylogarithm only implemented for n at least 2.
        """
        from sage.rings.power_series_ring import PowerSeriesRing
        from sage.arith.misc import integer_ceil as ceil, integer_floor as floor
        from sage.rings.padics.factory import Qp
        from sage.misc.verbose import verbose

        if self == 1:
            raise ValueError('Polylogarithm is not defined for 1.')
        if n <= 1:
            raise ValueError('Polylogarithm only implemented for n at least 2.')

        p = self.parent().prime()
        prec = self.precision_absolute()

        K = self.parent().fraction_field()
        z = K(self)

        hsl = max(prec / ((z - 1).valuation()) + 1, prec*(p == 2), 2)
        N = floor(prec - n*(hsl - 1).log(p).n())

        verbose(hsl, level=3)

        def bound(m):
            return prec - m + Integer(1-2**(m-1)).valuation(p) - m*(hsl - 1).log(p).n()

        gsl = max(_findprec(1/(p-1), 1, _polylog_c(m,p) + bound(m), p) for m in range(2,n+1))
        gsl = max(gsl, 2)
        verbose(gsl, level=3)
        g = _compute_g(p, n, max(bound(m) + m*floor((gsl-1).log(p).n()) for m in range(2, n+1)), gsl)
        verbose(g, level=3)
        S = PowerSeriesRing(K, default_prec = ceil(hsl), names='t')
        t = S.gen()
        log1plust = (1+t).log()
        log1plusti = 1

        G = (n+1)*[0]
        for i in range(n+1):
            G[i] = (log1plusti)/Integer(i).factorial()
            log1plusti *= log1plust

        verbose(G, level=3)

        H = (n+1)*[0]
        H[2] = -sum([((-t)**i)/i**2 for i in srange(1,hsl+2)])
        for i in range(2, n):
            H[i+1] = (H[i]/(1+t) + G[i]/t).integral()
            if (i + 1) % 2 == 1:
                if p != 2:
                    H[i+1] += (2**i*p**(i+1)*g[i+1](1/K(2)))/((1-2**i)*(p**(i+1) - 1))
                else:
                    H[i+1] += (2**i*H[i+1](K(-2)))/(1 - 2**(i+1))

        verbose(H, level=3)
        return (H[n](z - 1) - ((z.log(p_branch))**(n-1)*(1 - z).log(p_branch))/Integer(n-1).factorial()).add_bigoh(N)

    def polylog(self, n, p_branch=0):
        r"""
        Return `Li_n(\mathrm{self})`, the `n`-th `p`-adic polylogarithm
        of this element.

        INPUT:

        - ``n`` -- nonnegative integer
        - ``p_branch`` -- an element in the base ring or its fraction
          field; the implementation will choose the branch of the
          logarithm which sends `p` to ``p_branch``

        EXAMPLES:

        The `n`-th polylogarithm of `-1` is `0` for even `n`::

            sage: Qp(13)(-1).polylog(6) == 0                                            # needs sage.rings.real_mpfr sage.symbolic
            True

        We can check some identities, for example those mentioned in [DCW2016]_::

            sage: x = Qp(7, prec=30)(1/3)
            sage: (x^2).polylog(4) - 8*x.polylog(4) - 8*(-x).polylog(4) == 0            # needs sage.symbolic
            True

        ::

            sage: x = Qp(5, prec=30)(4)
            sage: x.polylog(2) + (1/x).polylog(2) + x.log(0)**2/2 == 0                  # needs sage.symbolic
            True

        ::

            sage: x = Qp(11, prec=30)(2)
            sage: x.polylog(2) + (1-x).polylog(2) + x.log(0)**2*(1-x).log(0) == 0       # needs sage.symbolic
            True

        `Li_1(z) = -\log(1-z)` for `|z| < 1`::

            sage: Qp(5)(10).polylog(1) == -Qp(5)(1-10).log(0)
            True

        The dilogarithm of 1 is zero::

            sage: Qp(5)(1).polylog(2)                                                   # needs sage.rings.real_mpfr sage.symbolic
            O(5^20)

        The cubing relation holds for the trilogarithm at 1::

            sage: K = Qp(7)
            sage: z = K.zeta(3)
            sage: -8*K(1).polylog(3) == 9*(K(z).polylog(3) + K(z^2).polylog(3))         # needs sage.rings.padics sage.rings.real_mpfr sage.symbolic
            True

        The polylogarithm of 0 is 0::

            sage: Qp(11)(0).polylog(7)
            0

        Only polylogarithms for positive `n` are defined::

            sage: Qp(11)(2).polylog(-1)
            Traceback (most recent call last):
            ...
            ValueError: polylogarithm only implemented for n at least 0

        Check that :issue:`29222` is fixed::

            sage: K = Qp(7)
            sage: print(K(1 + 7^11).polylog(4))                                         # needs sage.symbolic
            6*7^14 + 3*7^15 + 7^16 + 7^17 + O(7^18)

        ALGORITHM:

        The algorithm of Besser-de Jeu, as described in [BdJ2008]_ is used.

        AUTHORS:

        - Jennifer Balakrishnan - Initial implementation
        - Alex J. Best (2017-07-21) - Extended to other residue disks

        .. TODO::

            - Implement for extensions.
            - Use the change method to create K from ``self.parent()``.
        """
        from sage.rings.power_series_ring import PowerSeriesRing
        from sage.rings.padics.factory import Qp
        from sage.misc.verbose import verbose
        from sage.arith.misc import integer_ceil as ceil, integer_floor as floor
        from sage.rings.infinity import PlusInfinity

        if self.parent().absolute_degree() != 1:
            raise NotImplementedError("polylogarithms are not currently implemented for elements of extensions")
            # TODO implement this (possibly after the change method for padic generic elements is added).
        if n == 0:
            return self/(1-self)
        if n == 1:
            return -(1-self).log(p_branch)
        if n < 0:
            raise ValueError('polylogarithm only implemented for n at least 0')

        prec = self.precision_absolute()

        p = self.parent().prime()
        K = self.parent().fraction_field()

        z = K(self)
        n = Integer(n)

        if z.valuation() < 0:
            verbose("residue oo, using functional equation for reciprocal. %d %s"%(n,str(self)), level=2)
            return (-1)**(n+1)*(1/z).polylog(n)-(z.log(p_branch)**n)/K(n.factorial())

        zeta = K.teichmuller(z)

        # Which residue disk are we in?
        if zeta == 0:
            if z.precision_absolute() == PlusInfinity():
                return K(0)
            verbose("residue 0, using series. %d %s"%(n,str(self)), level=2)
            M = ceil((prec/z.valuation()).log(p).n())
            N = prec - n*M
            ret = K(0)
            for m in range(M + 1):
                zpm = z**(p**m)
                ret += p**(-m*n)*sum(zpm**k/Integer(k)**n for k in
                        range(1, _findprec(p**m*z.valuation(), 0, N + n*m, p)) if k % p != 0)

            return ret.add_bigoh(N)

        if zeta == 1:
            if z == 1:
                return Integer(2)**(n-1)*K(-1).polylog(n, p_branch=p_branch)/(1-Integer(2)**(n-1))
            verbose("residue 1, using _polylog_res_1. %d %s"%(n,str(self)), level=2)
            return self._polylog_res_1(n, p_branch)

        # Set up precision bounds
        tsl = prec / (z - zeta).valuation() + 1
        N = floor(prec - n*(tsl - 1).log(p).n())
        gsl = max(_findprec(1/(p-1), 1, prec - m + _polylog_c(m,p) - m*(tsl - 1).log(p).n(), p) for m in range(1,n+1))
        gsl = max(gsl,2)

        gtr = _compute_g(p, n, prec + n*(gsl - 1).log(p).n(), gsl)

        K = Qp(p, prec)

        # Residue disk around zeta
        verbose("general case. %d %s"%(n, str(self)), level=2)
        Li_i_zeta = [0] + [p**i/(p**i-1)*gtr[i](1/(1-zeta)) for i in range(1,n+1)]

        T = PowerSeriesRing(K, default_prec=ceil(tsl), names='t')
        t = T.gen()
        F = (n+1)*[0]
        F[0] = (zeta+t)/(1-zeta-t)
        for i in range(n):
            F[i+1] = Li_i_zeta[i+1] + (F[i]/(zeta + t)).integral()

        return (F[n](z - zeta)).add_bigoh(N)


# Artin-Hasse exponential
_AHE_coefficients_cache = {}


def _AHE_coefficients(p, N, prec):
    r"""
    Compute the first ``N`` coefficients of the ``p``-adic
    Artin-Hasse exponential series at precision ``prec``.

    The output is a list of coefficients. The common parent
    of these coefficients is the ring of ``p``-adic integers
    with fixed modulus (with some internal precision which
    could be strictly higher than ``prec``).

    The result is cached.

    EXAMPLES::

        sage: from sage.rings.padics.padic_generic_element import _AHE_coefficients

        sage: L = _AHE_coefficients(101, 10, 3); L
        [1,
         1,
         51 + 50*101 + 50*101^2,
         17 + 84*101 + 16*101^2,
         80 + 96*101 + 79*101^2,
         16 + 100*101 + 15*101^2,
         70 + 16*101 + 53*101^2,
         10 + 60*101 + 7*101^2,
         77 + 32*101 + 89*101^2,
         31 + 37*101 + 32*101^2]
        sage: L == [ 1/factorial(i) for i in range(10) ]
        True

    We check the parent::

        sage: [ elt.parent() for elt in L ]
        [101-adic Ring of fixed modulus 101^3,
         101-adic Ring of fixed modulus 101^3,
         101-adic Ring of fixed modulus 101^3,
         101-adic Ring of fixed modulus 101^3,
         101-adic Ring of fixed modulus 101^3,
         101-adic Ring of fixed modulus 101^3,
         101-adic Ring of fixed modulus 101^3,
         101-adic Ring of fixed modulus 101^3,
         101-adic Ring of fixed modulus 101^3,
         101-adic Ring of fixed modulus 101^3]

    Sometimes the precision on the result seems to be higher
    that the requested precision.
    However, the result is *not* guaranteed to be correct
    beyond the requested precision::

        sage: L = _AHE_coefficients(2, 513, 1); L                                       # needs sage.symbolic
        [1,
         1,
         1,
         2 + 2^2 + 2^4 + 2^6 + 2^8,
         ...
         1 + 2 + 2^2 + 2^5 + 2^8,
         2^2 + 2^6 + 2^9,
         1]

    We check that the result is correct modulo `2^1`::

        sage: S.<x> = PowerSeriesRing(QQ, 513)
        sage: AH = exp(sum(x^(2^i) / 2^i for i in range(10)))
        sage: R = ZpFM(2, 1)
        sage: [ R(c) for c in L ] == [ R(c) for c in AH.list() ]                        # needs sage.rings.padics
        True

    But it is not modulo `2^{10}`::

        sage: R = ZpFM(2, 10)
        sage: [ R(c) for c in L ] == [ R(c) for c in AH.list() ]                        # needs sage.rings.padics
        False
    """
    from sage.rings.padics.factory import ZpFM
    from sage.arith.misc import integer_floor as floor

    if N < p:
        internal_prec = prec
    else:
        internal_prec = prec + floor((N-1).log()/p.log())
    if p in _AHE_coefficients_cache:
        cache_internal_prec, values = _AHE_coefficients_cache[p]
    else:
        cache_internal_prec = 0
    if cache_internal_prec < internal_prec:
        parent = ZpFM(p, internal_prec)
        values = [ parent(1) ]
    for i in range(len(values), N):
        c = 0
        dec = 1
        while dec <= i:
            c += values[i-dec]
            dec *= p
        values.append(c // i)
    _AHE_coefficients_cache[p] = (internal_prec, values)
    return values


# Module functions used by polylog
def _polylog_c(n, p):
    """
    Return c(n, p) = p/(p-1) - (n-1)/log(p) + (n-1)*log(n*(p-1)/log(p),p) + log(2*p*(p-1)*n/log(p), p)
    as defined in Prop 6.1 of [BdJ2008]_ which is used as a precision bound.
    This is an internal function, used by :meth:`polylog`.

    EXAMPLES::

        sage: sage.rings.padics.padic_generic_element._polylog_c(1, 2)                  # needs sage.symbolic
        4.52876637294490
    """
    return p/(p-1) - (n-1)/p.log().n() + (n-1)*(n*(p-1)/p.log().n()).log(p).n() + (2*p*(p-1)*n/p.log().n()).log(p).n()


def _findprec(c_1, c_2, c_3, p):
    r"""
    Return an integer k such that c_1*k - c_2*log_p(k) > c_3.
    This is an internal function, used by :meth:`polylog`.

    INPUT:

    - `c_1`, `c_2`, `c_3` -- positive integers
    - ``p`` -- prime

    OUTPUT:

    ``sl`` such that `kc_1 - c_2 \log_p(k) > c_3` for all `k \ge sl`

    EXAMPLES::

        sage: sage.rings.padics.padic_generic_element._findprec(1, 1, 2, 2)             # needs sage.rings.real_double sage.rings.real_mpfr sage.symbolic
        5
        sage: 5*1 - 5*log(1, 2) > 2
        True

    See Remark 7.11 of [BdJ2008]_.
    """
    from sage.arith.misc import integer_ceil as ceil

    k = Integer(max(ceil(c_2/c_1), 2))
    while True:
        if c_1*k - c_2*k.log(p).n() > c_3:
            return k
        k += 1


def _compute_g(p, n, prec, terms):
    r"""
    Return the list of power series `g_i = \int(-g_{i-1}/(v-v^2))` used in the computation of polylogarithms.

    This is an internal function, used by :meth:`polylog`.

    EXAMPLES::

        sage: sage.rings.padics.padic_generic_element._compute_g(7, 3, 3, 3)[0]         # needs sage.libs.ntl sage.rings.real_double
        O(7^3)*v^2 + (1 + O(7^3))*v + O(7^3)
    """
    from sage.rings.power_series_ring import PowerSeriesRing
    from sage.arith.misc import integer_ceil as ceil
    from sage.rings.padics.factory import Qp

    # Compute the sequence of power series g
    R = PowerSeriesRing(Qp(p, ceil(prec)), default_prec=terms, names='v')
    v = R.gen()
    g = (n+1)*[0]
    g[0] = v - 1 - ((v-1)**p)/(v**p-(v-1)**p)
    for i in range(n):
        g[i+1] = -(g[i]/(v-v**2)).integral()
    return [x.truncate(terms) for x in g]


cpdef dwork_mahler_coeffs(R, int bd=20):
    r"""
    Compute Dwork's formula for Mahler coefficients of `p`-adic Gamma.

    This is called internally when one computes Gamma for a `p`-adic
    integer. Normally there is no need to call it directly.

    INPUT:

    - ``R`` -- `p`-adic ring in which to compute
    - ``bd`` -- integer; number of terms in the expansion to use

    OUTPUT: a list of `p`-adic integers

    EXAMPLES::

        sage: from sage.rings.padics.padic_generic_element import dwork_mahler_coeffs, evaluate_dwork_mahler
        sage: R = Zp(3)
        sage: v = dwork_mahler_coeffs(R)
        sage: x = R(1/7)
        sage: evaluate_dwork_mahler(v, x, 3, 20, 1)
        2 + 2*3 + 3^2 + 3^3 + 3^4 + 3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^11
         + 2*3^12 + 3^13 + 3^14 + 2*3^16 + 3^17 + 3^19 + O(3^20)
        sage: x.dwork_expansion(a=1)  # Same result
        2 + 2*3 + 3^2 + 3^3 + 3^4 + 3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^11
         + 2*3^12 + 3^13 + 3^14 + 2*3^16 + 3^17 + 3^19 + O(3^20)
    """
    from sage.rings.padics.factory import Qp
    cdef int i
    cdef long k, p

    v = [R.one()]
    p = R.prime()
    for k in range(1, p):
        v.append(v[-1] / R(k))
    if bd > 1:
        R1 = Qp(p, prec=bd) # Need divisions in this calculation
        u = [R1(x) for x in v]
        for i in range(1, bd):
            u[0] = ((u[-1] + u[0]) / i) >> 1
            for j in range(1, p):
                u[j] = (u[j-1] + u[j]) / (j + i * p)
            for x in u:
                v.append(R(x << i))
    return v

cpdef evaluate_dwork_mahler(v, x, long long p, int bd, long long a):
    """
    Evaluate Dwork's Mahler series for `p`-adic Gamma.

    EXAMPLES::

        sage: from sage.rings.padics.padic_generic_element import dwork_mahler_coeffs, evaluate_dwork_mahler
        sage: R = Zp(3)
        sage: v = dwork_mahler_coeffs(R)
        sage: x = R(1/7)
        sage: evaluate_dwork_mahler(v, x, 3, 20, 1)
        2 + 2*3 + 3^2 + 3^3 + 3^4 + 3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^11
         + 2*3^12 + 3^13 + 3^14 + 2*3^16 + 3^17 + 3^19 + O(3^20)
        sage: x.dwork_expansion(a=1)  # Same result
        2 + 2*3 + 3^2 + 3^3 + 3^4 + 3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^11
         + 2*3^12 + 3^13 + 3^14 + 2*3^16 + 3^17 + 3^19 + O(3^20)
    """
    cdef int k
    bd -= 1
    a1 = a + bd*p
    s = v[a1]
    u = x + bd
    one = x.parent().one()
    for k in range(bd):
        a1 -= p
        u -= one
        s = s*u + v[a1]
    return -s

cdef long long evaluate_dwork_mahler_long(array.array v, long long x, long long p, int bd,
                                     long long a, long long q) noexcept:
    cdef int k
    cdef long long a1, s, u
    bd -= 1
    a1 = a + bd*p
    s = v[a1]
    u = x + bd
    for k in range(bd):
        a1 -= p
        u -= 1
        s = s*u + v[a1] # force cast to long long
        s = s % q
    return -s

cpdef gauss_table(long long p, int f, int prec, bint use_longs):
    r"""
    Compute a table of Gauss sums using the Gross-Koblitz formula.

    This is used in the computation of `L`-functions of hypergeometric motives.
    The Gross-Koblitz formula is used as in `sage.rings.padics.misc.gauss_sum`,
    but further unpacked for efficiency.

    INPUT:

    - ``p`` -- prime
    - ``f``, ``prec`` -- positive integers
    - ``use_longs`` -- boolean; if ``True``, computations are done in C ``long long``
        integers rather than Sage `p`-adics, and the results are returned
        as a Python array rather than a list.

    OUTPUT:

    A list of length `q-1=p^f-1`. The entries are `p`-adic units created with
    absolute precision ``prec``.

    EXAMPLES::

        sage: from sage.rings.padics.padic_generic_element import gauss_table
        sage: gauss_table(2,2,4,False)
        [1 + 2 + 2^2 + 2^3, 1 + 2 + 2^2 + 2^3, 1 + 2 + 2^2 + 2^3]
        sage: gauss_table(3,2,4,False)[3]
        2 + 3 + 2*3^2
    """
    from sage.rings.padics.factory import Zp, Qp

    cdef int i, j, bd
    cdef long long q, q1, q3, r, r1, r2, s1, s2, k
    cdef array.array vv, ans1

    if (f == 1 and prec == 1): # Shortcut for this key special case
        ans1 = array.array('l', [0]) * p
        ans1[0] = p-1
        for r in range(1, p-1):
            k = ans1[r-1]
            ans1[r] = k * r % p
        return ans1

    q = p ** f
    q1 = q - 1
    bd = (p*prec+p-2) // (p-1) - 1
    R = Zp(p, prec, 'fixed-mod')
    if p == 2: # Dwork expansion has denominators when p = 2
        R1 = Qp(p, prec)
        use_longs = False
    else:
        R1 = R
    d = ~R1(q1)
    v = dwork_mahler_coeffs(R1, bd)
    if use_longs:
        q3 = p ** prec
        r2 = d.lift() % q3
        vv = array.array('l', [0]) * len(v)
        for k in range(len(v)):
            vv[k] = v[k].lift() % q3
        ans1 = array.array('l', [0]) * q1
        ans1[0] = -1
        ans = ans1
    else:
        u = R1.one()
        ans = [0 for r in range(q1)]
        ans[0] = -u
    for r in range(1, q1):
        if ans[r]: continue
        if use_longs:
            s1 = 1
        else:
            s = u
        r1 = r
        for j in range(1, f+1):
            k = r1 % p
            r1 = (r1 + k * q1) // p
            if use_longs: # Use Dwork expansion to compute p-adic Gamma
                s1 *= -evaluate_dwork_mahler_long(vv, r1*r2%q3, p, bd, k, q3)
                s1 %= q3
            else:
                s *= -evaluate_dwork_mahler(v, R1(r1)*d, p, bd, k)
            if r1 == r:
                break
        if use_longs:
            if j < f:
                s2 = s1
                for i in range(f//j-1):
                    s1 = s1 * s2 % q3
            ans1[r] = -s1
        else:
            if j < f:
                s **= f // j
            ans[r] = -s
        for i in range(j-1):
            r1 = r1 * p % q1 # Initially r1 == r
            ans[r1] = ans[r]
    if p != 2: return ans
    return [R(x) for x in ans]
