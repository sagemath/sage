# sage.doctest: needs sage.rings.padics
"""
Local Generic Element

This file contains a common superclass for `p`-adic elements and power
series elements.

AUTHORS:

- David Roe: initial version

- Julian Rueth (2012-10-15, 2014-06-25, 2017-08-04): added inverse_of_unit(); improved
  add_bigoh(); added _test_expansion()
"""
# ****************************************************************************
#       Copyright (C) 2007-2017 David Roe <roed@math.harvard.edu>
#                     2012-2017 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.infinity import infinity
from sage.structure.element cimport CommutativeRingElement
from sage.structure.element import coerce_binop
from itertools import islice


cdef class LocalGenericElement(CommutativeRingElement):
    #cpdef _add_(self, right):
    #    raise NotImplementedError

    cpdef _div_(self, right):
        r"""
        Returns the quotient of ``self`` by ``right``.

        INPUT:

        - ``self`` -- a `p`-adic element.

        - ``right`` -- a `p`-adic element distinguishable from zero.
          In a fixed-modulus ring, this element must be a unit.

        EXAMPLES::

            sage: R = Zp(7, 4, 'capped-rel', 'series'); R(3)/R(5)
            2 + 4*7 + 5*7^2 + 2*7^3 + O(7^4)
            sage: R(2/3) / R(1/3) #indirect doctest
            2 + O(7^4)
            sage: R(49) / R(7)
            7 + O(7^5)
            sage: R = Zp(7, 4, 'capped-abs', 'series'); 1/R(7)
            7^-1 + O(7^2)
            sage: R = Zp(7, 4, 'fixed-mod'); 1/R(7)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
        """
        # this doctest doesn't actually test the function, since it's overridden.
        return self * ~right

    def inverse_of_unit(self):
        r"""
        Returns the inverse of ``self`` if ``self`` is a unit.

        OUTPUT:

            - an element in the same ring as ``self``

        EXAMPLES::

            sage: R = ZpCA(3,5)
            sage: a = R(2); a
            2 + O(3^5)
            sage: b = a.inverse_of_unit(); b
            2 + 3 + 3^2 + 3^3 + 3^4 + O(3^5)

        A ``ZeroDivisionError`` is raised if an element has no inverse in the
        ring::

            sage: R(3).inverse_of_unit()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: inverse of 3 + O(3^5) does not exist

        Unlike the usual inverse of an element, the result is in the same ring
        as ``self`` and not just in its fraction field::

            sage: c = ~a; c
            2 + 3 + 3^2 + 3^3 + 3^4 + O(3^5)
            sage: a.parent()
            3-adic Ring with capped absolute precision 5
            sage: b.parent()
            3-adic Ring with capped absolute precision 5
            sage: c.parent()
            3-adic Field with capped relative precision 5

        For fields this does of course not make any difference::

            sage: R = QpCR(3,5)
            sage: a = R(2)
            sage: b = a.inverse_of_unit()
            sage: c = ~a
            sage: a.parent()
            3-adic Field with capped relative precision 5
            sage: b.parent()
            3-adic Field with capped relative precision 5
            sage: c.parent()
            3-adic Field with capped relative precision 5

        TESTS:

        Test that this works for all kinds of p-adic base elements::

            sage: ZpCA(3,5)(2).inverse_of_unit()
            2 + 3 + 3^2 + 3^3 + 3^4 + O(3^5)
            sage: ZpCR(3,5)(2).inverse_of_unit()
            2 + 3 + 3^2 + 3^3 + 3^4 + O(3^5)
            sage: ZpFM(3,5)(2).inverse_of_unit()
            2 + 3 + 3^2 + 3^3 + 3^4
            sage: ZpFP(3,5)(2).inverse_of_unit()
            2 + 3 + 3^2 + 3^3 + 3^4
            sage: QpCR(3,5)(2).inverse_of_unit()
            2 + 3 + 3^2 + 3^3 + 3^4 + O(3^5)

        Over unramified extensions::

            sage: R = ZpCA(3,5); S.<t> = R[]; W.<t> = R.extension( t^2 + 1 )
            sage: t.inverse_of_unit()
            2*t + 2*t*3 + 2*t*3^2 + 2*t*3^3 + 2*t*3^4 + O(3^5)

            sage: R = ZpCR(3,5); S.<t> = R[]; W.<t> = R.extension( t^2 + 1 )
            sage: t.inverse_of_unit()
            2*t + 2*t*3 + 2*t*3^2 + 2*t*3^3 + 2*t*3^4 + O(3^5)

            sage: R = ZpFM(3,5); S.<t> = R[]; W.<t> = R.extension( t^2 + 1 )
            sage: t.inverse_of_unit()
            2*t + 2*t*3 + 2*t*3^2 + 2*t*3^3 + 2*t*3^4

            sage: R = ZpFP(3,5); S.<t> = R[]; W.<t> = R.extension( t^2 + 1 )
            sage: t.inverse_of_unit()
            2*t + 2*t*3 + 2*t*3^2 + 2*t*3^3 + 2*t*3^4

            sage: R = QpCR(3,5); S.<t> = R[]; W.<t> = R.extension( t^2 + 1 )
            sage: t.inverse_of_unit()
            2*t + 2*t*3 + 2*t*3^2 + 2*t*3^3 + 2*t*3^4 + O(3^5)

        Over Eisenstein extensions::

            sage: R = ZpCA(3,5); S.<t> = R[]; W.<t> = R.extension( t^2 - 3 )
            sage: (t - 1).inverse_of_unit()
            2 + 2*t + t^2 + t^3 + t^4 + t^5 + t^6 + t^7 + t^8 + t^9 + O(t^10)

            sage: R = ZpCR(3,5); S.<t> = R[]; W.<t> = R.extension( t^2 - 3 )
            sage: (t - 1).inverse_of_unit()
            2 + 2*t + t^2 + t^3 + t^4 + t^5 + t^6 + t^7 + t^8 + t^9 + O(t^10)

            sage: R = ZpFM(3,5); S.<t> = R[]; W.<t> = R.extension( t^2 - 3 )
            sage: (t - 1).inverse_of_unit()
            2 + 2*t + t^2 + t^3 + t^4 + t^5 + t^6 + t^7 + t^8 + t^9

            sage: R = QpCR(3,5); S.<t> = R[]; W.<t> = R.extension( t^2 - 3 )
            sage: (t - 1).inverse_of_unit()
            2 + 2*t + t^2 + t^3 + t^4 + t^5 + t^6 + t^7 + t^8 + t^9 + O(t^10)

        """
        if not self.is_unit():
            raise ZeroDivisionError(f"inverse of {self} does not exist")
        return self.parent()(~self)

    def __iter__(self):
        """
        Local elements should not be iterable, so this method correspondingly
        raises a ``TypeError``.

        .. NOTE::

            Typically, local elements provide a implementation for
            ``__getitem__``. If they do not provide a method ``__iter__``, then
            iterating over them is realized by calling ``__getitem__``,
            starting from index 0. However, there are several issues with this.
            For example, terms with negative valuation would be excluded from
            the iteration, and an exact value of zero would lead to an infinite
            iterable.

            There doesn't seem to be an obvious behaviour that iteration over
            such elements should produce, so it is disabled; see :trac:`13592`.

        TESTS::

            sage: x = Qp(3).zero()
            sage: for v in x: pass
            Traceback (most recent call last):
            ...
            TypeError: this local element is not iterable

        """
        raise TypeError("this local element is not iterable")

    def slice(self, i, j, k = 1, lift_mode='simple'):
        r"""
        Returns the sum of the `pi^{i + l \cdot k}` terms of the series
        expansion of this element, where pi is the uniformizer,
        for `i + l \cdot k` between ``i`` and ``j-1`` inclusive, and
        nonnegative integers `l`. Behaves analogously to the slice
        function for lists.

        INPUT:

        - ``i`` -- an integer; if set to ``None``, the sum will start with the
          first non-zero term of the series.

        - ``j`` -- an integer; if set to ``None`` or `\infty`, this method
          behaves as if it was set to the absolute precision of this element.

        - ``k`` -- (default: 1) a positive integer

        EXAMPLES::

            sage: R = Zp(5, 6, 'capped-rel')
            sage: a = R(1/2); a
            3 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 2*5^5 + O(5^6)
            sage: a.slice(2, 4)
            2*5^2 + 2*5^3 + O(5^4)
            sage: a.slice(1, 6, 2)
            2*5 + 2*5^3 + 2*5^5 + O(5^6)

        The step size ``k`` has to be positive::

            sage: a.slice(0, 3, 0)
            Traceback (most recent call last):
            ...
            ValueError: slice step must be positive
            sage: a.slice(0, 3, -1)
            Traceback (most recent call last):
            ...
            ValueError: slice step must be positive

        If ``i`` exceeds ``j``, then the result will be zero, with the
        precision given by ``j``::

            sage: a.slice(5, 4)
            O(5^4)
            sage: a.slice(6, 5)
            O(5^5)

        However, the precision cannot exceed the precision of the element::

            sage: a.slice(101,100)
            O(5^6)
            sage: a.slice(0,5,2)
            3 + 2*5^2 + 2*5^4 + O(5^5)
            sage: a.slice(0,6,2)
            3 + 2*5^2 + 2*5^4 + O(5^6)
            sage: a.slice(0,7,2)
            3 + 2*5^2 + 2*5^4 + O(5^6)

        If start is left blank, it is set to the valuation::

            sage: K = Qp(5, 6)
            sage: x = K(1/25 + 5); x
            5^-2 + 5 + O(5^4)
            sage: x.slice(None, 3)
            5^-2 + 5 + O(5^3)
            sage: x[:3]
            doctest:warning
            ...
            DeprecationWarning: __getitem__ is changing to match the behavior of number fields. Please use expansion instead.
            See https://github.com/sagemath/sage/issues/14825 for details.
            5^-2 + 5 + O(5^3)

        TESTS:

        Test that slices also work over fields::

            sage: a = K(1/25); a
            5^-2 + O(5^4)
            sage: b = K(25); b
            5^2 + O(5^8)

            sage: a.slice(2, 4)
            O(5^4)
            sage: b.slice(2, 4)
            5^2 + O(5^4)
            sage: a.slice(-3, -1)
            5^-2 + O(5^-1)
            sage: b.slice(-1, 1)
            O(5)
            sage: b.slice(-3, -1)
            O(5^-1)
            sage: b.slice(101, 100)
            O(5^8)
            sage: b.slice(0,7,2)
            5^2 + O(5^7)
            sage: b.slice(0,8,2)
            5^2 + O(5^8)
            sage: b.slice(0,9,2)
            5^2 + O(5^8)

        Test that slices also work over eisenstein extensions::

            sage: F = Qp(5)
            sage: H.<x> = F[]
            sage: T.<t> = F.extension(x^2-5)
            sage: a = T(3*t^-2 + 1 + 4*t + 2*t^2)
            sage: a.slice(0, 1)
            1 + O(t)
            sage: a.slice(-3, 4)
            3*t^-2 + 1 + 4*t + 2*t^2 + O(t^4)
            sage: a.slice(-2, 6, 3)
            3*t^-2 + 4*t + O(t^6)

        Test that slices also work over unramified extensions::

            sage: F = Qp(5)
            sage: H.<x> = F[]
            sage: T.<t> = F.extension(x^2-2)
            sage: a = T(3*5^-1 + 1 + (3*t + 4)*5^2)
            sage: a.slice(0, 1)
            1 + O(5)
            sage: a.slice(-3, 4)
            3*5^-1 + 1 + (3*t + 4)*5^2 + O(5^4)
            sage: a.slice(-1, 6, 3)
            3*5^-1 + (3*t + 4)*5^2 + O(5^6)

        Test that slices also work over 2-step extensions (unramified followed by eisenstein)::

            sage: F = Qp(5)
            sage: H.<x> = F[]
            sage: T.<t> = F.extension(x^2-3)
            sage: D.<y> = T[]
            sage: W.<w> = T.extension((4*5^-2 + 2*5^-1 + 4 + (2*t + 2)*5 + 3*t*5^3 + 4*5^4 + 3*5^5 + (2*t + 2)*5^8 + (4*t + 3)*5^9 + 2*t*5^10 + (3*t + 3)*5^11 + (3*t + 1)*5^12 + (3*t + 2)*5^13 + 4*5^14 + (2*t + 4)*5^15 + (4*t + 1)*5^16 + (t + 1)*5^17 + O(5^18))*y^2 + (t + 2*t*5 + t*5^2 + 4*t*5^3 + (2*t + 4)*5^4 + (3*t + 4)*5^5 + (t + 1)*5^6 + t*5^7 + (2*t + 4)*5^8 + 3*5^9 + 2*5^10 + 5^12 + (4*t + 2)*5^13 + 5^14 + 5^15 + 3*t*5^16 + (t + 2)*5^17 + 4*5^18 + (3*t + 1)*5^19 + O(5^20))*y + (2*t + 2)*5^-1 + 3 + 5 + t*5^2 + (4*t + 2)*5^3 + (4*t + 1)*5^4 + (3*t + 4)*5^5 + (4*t + 4)*5^6 + (3*t + 2)*5^7 + (4*t + 4)*5^8 + 3*5^9 + (t + 3)*5^10 + (4*t + 3)*5^11 + 5^12 + (2*t + 2)*5^14 + 4*t*5^15 + (2*t + 2)*5^16 + (4*t + 4)*5^17 + O(5^18))
            sage: a = W(3*w^-36 + (2*t + 2)*w^-23)
            sage: a.slice(-25,2)
            (2*t + 2)*w^-23 + O(w^2)
            sage: a.slice(0, 1)
            O(w)

        Verify that :trac:`14106` has been fixed::

            sage: R = Zp(5,7)
            sage: a = R(300)
            sage: a
            2*5^2 + 2*5^3 + O(5^9)
            sage: a[:5]
            2*5^2 + 2*5^3 + O(5^5)
            sage: a.slice(None, 5, None)
            2*5^2 + 2*5^3 + O(5^5)

        Verify that :trac:`30695` has been fixed::

            sage: F=Qp(3)
            sage: a=F(0)
            sage: a.slice(0,None)
            0

        """
        if k is None:
            k = 1
        if k <= 0:
            raise ValueError("slice step must be positive")
        if i is None:
            i = self.valuation()
        if j is None or j is infinity:
            j = self.precision_absolute()
            if j is infinity:
                return self.parent()(0)

        start = i
        stop = j

        # for fields, self.list() contains only the coefficients starting from
        # self.valuation(), so we have to shift the indices around to make up
        # for this
        if self.parent().is_field():
            start -= self.valuation()
            stop -= self.valuation()

        # make sure that start and stop are non-negative
        if start<0:
            i += -start # fix the value of ppow below
            start = 0
        stop = max(stop, 0)

        # the increase of the pi-power in every step
        pk = self.parent().uniformizer_pow(k)
        # the pi-power of the first term
        ppow = self.parent().uniformizer_pow(i)

        # construct the return value
        ans = self.parent().zero()
        unramified_generator = self.parent()(self.parent().residue_field().gen()).lift_to_precision()
        for c in islice(self.expansion(lift_mode=lift_mode), int(start), int(stop), int(k)):
            genpow = 1
            if not isinstance(c, list): c = [c] # relevant for the case of base-rings, or one-step
                                                # eisenstein extensions
            for d in c:
                ans += d * genpow * ppow
                genpow *= unramified_generator
            ppow *= pk

        # fix the precision of the return value
        if j < ans.precision_absolute() or self.precision_absolute() < ans.precision_absolute():
            ans = ans.add_bigoh(min(j, self.precision_absolute()))

        return ans

    def _latex_(self):
        r"""
        Return a latex representation of self.

        EXAMPLES::

            sage: R = Zp(5); a = R(17)
            sage: latex(a) #indirect doctest
            2 + 3 \cdot 5 + O(5^{20})
        """
        # TODO: add a bunch more documentation of latexing elements
        return self._repr_(do_latex = True)

    #def __mod__(self, right):
    #    raise NotImplementedError

    #cpdef _mul_(self, right):
    #    raise NotImplementedError

    #cdef _neg_(self):
    #    raise NotImplementedError

    #def __pow__(self, right):
    #    raise NotImplementedError

    cpdef _sub_(self, right):
        r"""
        Returns the difference between ``self`` and ``right``.

        EXAMPLES::

            sage: R = Zp(7, 4, 'capped-rel', 'series'); a = R(12); b = R(5); a - b
            7 + O(7^4)
            sage: R(4/3) - R(1/3) #indirect doctest
            1 + O(7^4)
        """
        # this doctest doesn't actually test this function, since _sub_ is overridden.
        return self + (-right)

    def add_bigoh(self, absprec):
        """
        Return a copy of this element with absolute precision decreased to
        ``absprec``.

        INPUT:

        - ``absprec`` -- an integer or positive infinity

        EXAMPLES::

            sage: K = QpCR(3,4)
            sage: o = K(1); o
            1 + O(3^4)
            sage: o.add_bigoh(2)
            1 + O(3^2)
            sage: o.add_bigoh(-5)
            O(3^-5)

        One cannot use ``add_bigoh`` to lift to a higher precision; this
        can be accomplished with :meth:`lift_to_precision`::

            sage: o.add_bigoh(5)
            1 + O(3^4)

        Negative values of ``absprec`` return an element in the fraction field
        of the element's parent::

            sage: R = ZpCA(3,4)
            sage: R(3).add_bigoh(-5)
            O(3^-5)

        For fixed-mod elements this method truncates the element::

            sage: R = ZpFM(3,4)
            sage: R(3).add_bigoh(1)
            0

        If ``absprec`` exceeds the precision of the element, then this method
        has no effect::

            sage: R(3).add_bigoh(5)
            3

        A negative value for ``absprec`` returns an element in the fraction field::

            sage: R(3).add_bigoh(-1).parent()
            3-adic Field with floating precision 4

        TESTS:

        Test that this also works for infinity::

            sage: R = ZpCR(3,4)
            sage: R(3).add_bigoh(infinity)
            3 + O(3^5)
            sage: R(0).add_bigoh(infinity)
            0

        Check that :trac:`23464` has been resolved::

            sage: R.<pi> = Qp(7).extension(x^3 - 7)
            sage: (pi^93).add_bigoh(-10)
            O(pi^-10)

        """
        parent = self.parent()
        if absprec >= self.precision_absolute():
            return self
        if absprec < 0:
            parent = parent.fraction_field()
        return parent(self, absprec=absprec)

    #def copy(self):
    #    raise NotImplementedError

    #def exp(self):
    #    raise NotImplementedError

    def is_integral(self):
        """
        Returns whether self is an integral element.

        INPUT:

        - ``self`` -- a local ring element

        OUTPUT:

        - boolean -- whether ``self`` is an integral element.

        EXAMPLES::

            sage: R = Qp(3,20)
            sage: a = R(7/3); a.is_integral()
            False
            sage: b = R(7/5); b.is_integral()
            True
        """
        return self.valuation() >= 0

    #def is_square(self):
    #    raise NotImplementedError

    def is_padic_unit(self):
        """
        Returns whether self is a `p`-adic unit. That is, whether it has zero valuation.

        INPUT:

        - ``self`` -- a local ring element

        OUTPUT:

        - boolean -- whether ``self`` is a unit

        EXAMPLES::

            sage: R = Zp(3,20,'capped-rel'); K = Qp(3,20,'capped-rel')
            sage: R(0).is_padic_unit()
            False
            sage: R(1).is_padic_unit()
            True
            sage: R(2).is_padic_unit()
            True
            sage: R(3).is_padic_unit()
            False
            sage: Qp(5,5)(5).is_padic_unit()
            False

        TESTS::

            sage: R(4).is_padic_unit()
            True
            sage: R(6).is_padic_unit()
            False
            sage: R(9).is_padic_unit()
            False
            sage: K(0).is_padic_unit()
            False
            sage: K(1).is_padic_unit()
            True
            sage: K(2).is_padic_unit()
            True
            sage: K(3).is_padic_unit()
            False
            sage: K(4).is_padic_unit()
            True
            sage: K(6).is_padic_unit()
            False
            sage: K(9).is_padic_unit()
            False
            sage: K(1/3).is_padic_unit()
            False
            sage: K(1/9).is_padic_unit()
            False
            sage: Qq(3^2,5,names='a')(3).is_padic_unit()
            False
        """
        return self.valuation() == 0

    def is_unit(self):
        """
        Returns whether self is a unit

        INPUT:

        - ``self`` -- a local ring element

        OUTPUT:

        - boolean -- whether ``self`` is a unit

        .. NOTE::

            For fields all nonzero elements are units. For DVR's, only
            those elements of valuation 0 are. An older implementation
            ignored the case of fields, and returned always the
            negation of self.valuation()==0. This behavior is now
            supported with self.is_padic_unit().

        EXAMPLES::

            sage: R = Zp(3,20,'capped-rel'); K = Qp(3,20,'capped-rel')
            sage: R(0).is_unit()
            False
            sage: R(1).is_unit()
            True
            sage: R(2).is_unit()
            True
            sage: R(3).is_unit()
            False
            sage: Qp(5,5)(5).is_unit() # Note that 5 is invertible in `QQ_5`, even if it has positive valuation!
            True
            sage: Qp(5,5)(5).is_padic_unit()
            False

        TESTS::

            sage: R(4).is_unit()
            True
            sage: R(6).is_unit()
            False
            sage: R(9).is_unit()
            False
            sage: K(0).is_unit()
            False
            sage: K(1).is_unit()
            True
            sage: K(2).is_unit()
            True
            sage: K(3).is_unit()
            True
            sage: K(4).is_unit()
            True
            sage: K(6).is_unit()
            True
            sage: K(9).is_unit()
            True
            sage: K(1/3).is_unit()
            True
            sage: K(1/9).is_unit()
            True
            sage: Qq(3^2,5,names='a')(3).is_unit()
            True
            sage: R(0,0).is_unit()
            False
            sage: K(0,0).is_unit()
            False
        """
        if self.is_zero():
            return False
        if self.parent().is_field():
            return True
        return self.valuation() == 0

    #def is_zero(self, prec):
    #    raise NotImplementedError

    #def is_equal_to(self, right, prec):
    #    raise NotImplementedError

    #def lift(self):
    #    raise NotImplementedError

    #def list(self):
    #    raise NotImplementedError

    #def log(self):
    #    raise NotImplementedError

    #def multiplicative_order(self, prec):
    #    raise NotImplementedError

    #def padded_list(self):
    #    raise NotImplementedError

    #def precision_absolute(self):
    #    raise NotImplementedError

    #def precision_relative(self):
    #    raise NotImplementedError

    #def residue(self, prec):
    #    raise NotImplementedError

    def sqrt(self, extend=True, all=False, algorithm=None):
        r"""
        Return the square root of this element.

        INPUT:

        - ``self`` -- a `p`-adic element.

        - ``extend`` -- a boolean (default: ``True``); if ``True``, return a
          square root in an extension if necessary; if ``False`` and no root
          exists in the given ring or field, raise a ValueError.

        - ``all`` -- a boolean (default: ``False``); if ``True``, return a
          list of all square roots.

        - ``algorithm`` -- ``"pari"``, ``"sage"`` or ``None`` (default:
          ``None``); Sage provides an implementation for any extension of
          `Q_p` whereas only square roots over `Q_p` is implemented in Pari;
          the default is ``"pari"`` if the ground field is `Q_p`, ``"sage"``
          otherwise.

        OUTPUT:

        The square root or the list of all square roots of this element.

        .. NOTE::

            The square root is chosen (resp. the square roots are ordered) in
            a deterministic way, which is compatible with change of precision.

        EXAMPLES::

            sage: R = Zp(3, 20)
            sage: sqrt(R(0))
            0

            sage: sqrt(R(1))
            1 + O(3^20)

            sage: R(2).sqrt(extend=False)
            Traceback (most recent call last):
            ...
            ValueError: element is not a square

            sage: s = sqrt(R(4)); -s
            2 + O(3^20)

            sage: s = sqrt(R(9)); s
            3 + O(3^21)

        Over the `2`-adics, the precision of the square root is less
        than the input::

            sage: R2 = Zp(2, 20)
            sage: sqrt(R2(1))
            1 + O(2^19)
            sage: sqrt(R2(4))
            2 + O(2^20)

            sage: R.<t> = Zq(2^10, 10)
            sage: u = 1 + 8*t
            sage: sqrt(u)
            1 + t*2^2 + t^2*2^3 + t^2*2^4 + (t^4 + t^3 + t^2)*2^5 + (t^4 + t^2)*2^6 + (t^5 + t^2)*2^7 + (t^6 + t^5 + t^4 + t^2)*2^8 + O(2^9)

            sage: R.<a> = Zp(2).extension(x^3 - 2)
            sage: u = R(1 + a^4 + a^5 + a^7 + a^8, 10); u
            1 + a^4 + a^5 + a^7 + a^8 + O(a^10)
            sage: v = sqrt(u); v
            1 + a^2 + a^4 + a^6 + O(a^7)

        However, observe that the precision increases to its original value
        when we recompute the square of the square root::

            sage: v^2
            1 + a^4 + a^5 + a^7 + a^8 + O(a^10)

        If the input does not have enough precision in order to determine if
        the given element has a square root in the ground field, an error is
        raised::

            sage: R(1, 6).sqrt()
            Traceback (most recent call last):
            ...
            PrecisionError: not enough precision to be sure that this element has a square root

            sage: R(1, 7).sqrt()
            1 + O(a^4)

            sage: R(1+a^6, 7).sqrt(extend=False)
            Traceback (most recent call last):
            ...
            ValueError: element is not a square

        In particular, an error is raised when we try to compute the square
        root of an inexact

        TESTS::

            sage: R = Qp(5, 100)
            sage: c = R.random_element()
            sage: s = sqrt(c^2)
            sage: s == c or s == -c
            True

            sage: c2 = c^2
            sage: c2 = c2.add_bigoh(c2.valuation() + 50)
            sage: s == sqrt(c2)
            True
        """
        return self.square_root(extend, all)

    #def square_root(self, extend = True, all = False):
    #    raise NotImplementedError

    #def unit_part(self):
    #    raise NotImplementedError

    #def valuation(self):
    #    raise NotImplementedError

    def normalized_valuation(self):
        r"""
        Returns the normalized valuation of this local ring element,
        i.e., the valuation divided by the absolute ramification index.

        INPUT:

        ``self`` -- a local ring element.

        OUTPUT:

        rational -- the normalized valuation of ``self``.

        EXAMPLES::

            sage: Q7 = Qp(7)
            sage: R.<x> = Q7[]
            sage: F.<z> = Q7.ext(x^3+7*x+7)
            sage: z.normalized_valuation()
            1/3
        """
        F = self.parent()
        return self.valuation()/F.absolute_e()

    def _min_valuation(self):
        r"""
        Returns the valuation of this local ring element.

        This function only differs from valuation for relaxed elements.

        INPUT:

        - ``self`` -- a local ring element.

        OUTPUT:

        - integer -- the valuation of ``self``.

        EXAMPLES::

            sage: R = Qp(7, 4, 'capped-rel', 'series')
            sage: R(7)._min_valuation()
            1
            sage: R(1/7)._min_valuation()
            -1
        """
        return self.valuation()

    def euclidean_degree(self):
        r"""
        Return the degree of this element as an element of an Euclidean domain.

        EXAMPLES:

        For a field, this is always zero except for the zero element::

            sage: K = Qp(2)
            sage: K.one().euclidean_degree()
            0
            sage: K.gen().euclidean_degree()
            0
            sage: K.zero().euclidean_degree()
            Traceback (most recent call last):
            ...
            ValueError: euclidean degree not defined for the zero element

        For a ring which is not a field, this is the valuation of the element::

            sage: R = Zp(2)
            sage: R.one().euclidean_degree()
            0
            sage: R.gen().euclidean_degree()
            1
            sage: R.zero().euclidean_degree()
            Traceback (most recent call last):
            ...
            ValueError: euclidean degree not defined for the zero element
        """
        if self.is_zero():
            raise ValueError("euclidean degree not defined for the zero element")

        from sage.categories.fields import Fields
        if self.parent() in Fields():
            from sage.rings.integer import Integer
            return Integer(0)
        return self.valuation()

    @coerce_binop
    def quo_rem(self, other, integral=False):
        r"""
        Return the quotient with remainder of the division of this element by
        ``other``.

        INPUT:

        - ``other`` -- an element in the same ring
        - ``integral`` -- if True, use integral-style remainders even when the parent is a field.
          Namely, the remainder will have no terms in its p-adic expansion above
          the valuation of ``other``.

        EXAMPLES::

            sage: R = Zp(3, 5)
            sage: R(12).quo_rem(R(2))
            (2*3 + O(3^6), 0)
            sage: R(2).quo_rem(R(12))
            (O(3^4), 2 + O(3^5))

            sage: K = Qp(3, 5)
            sage: K(12).quo_rem(K(2))
            (2*3 + O(3^6), 0)
            sage: K(2).quo_rem(K(12))
            (2*3^-1 + 1 + 3 + 3^2 + 3^3 + O(3^4), 0)

        You can get the same behavior for fields as for rings
        by using integral=True::

            sage: K(12).quo_rem(K(2), integral=True)
            (2*3 + O(3^6), 0)
            sage: K(2).quo_rem(K(12), integral=True)
            (O(3^4), 2 + O(3^5))
        """
        if other.is_zero():
            raise ZeroDivisionError

        from sage.categories.fields import Fields
        if not integral and self.parent() in Fields():
            return (self / other, self.parent().zero())
        else:
            return self._quo_rem(other)

    def _test_trivial_powers(self, **options):
        r"""
        Check that taking trivial powers of elements works as expected.

        EXAMPLES::

            sage: x = Zp(3, 5).zero()
            sage: x._test_trivial_powers()

        """
        tester = self._tester(**options)

        x = self**1
        tester.assertEqual(x, self)
        tester.assertEqual(x.precision_absolute(), self.precision_absolute())

        z = self**0
        one = self.parent().one()
        tester.assertEqual(z, one)
        tester.assertEqual(z.precision_absolute(), one.precision_absolute())

    def _test_expansion(self, **options):
        r"""
        Check that ``expansion`` works as expected.

        EXAMPLES::

            sage: x = Zp(3, 5).zero()
            sage: x._test_expansion()

        """
        tester = self._tester(**options)

        shift = self.parent().one()
        v = 0
        # so that this test doesn't take too long for large precision cap
        prec_cutoff = int(min((10000 / (1 + self.precision_relative())).ceil(), 100))

        from sage.categories.fields import Fields
        if self.parent() in Fields():
            v = self.valuation()
            from sage.rings.infinity import infinity
            if self.valuation() is not infinity:
                shift = shift << v

        if self.parent().is_lattice_prec() or self.parent().is_relaxed():
            modes = ['simple']
        else:
            modes = ['simple', 'smallest', 'teichmuller']
        for mode in modes:
            expansion = self.expansion(lift_mode=mode)
            expansion_sum = sum(self.parent().maximal_unramified_subextension()(c) *
                                (self.parent().one()<<i)
                                for i, c in enumerate(islice(expansion, prec_cutoff))) * shift

            tester.assertEqual(self.add_bigoh(prec_cutoff), expansion_sum.add_bigoh(prec_cutoff))

            for i, c in enumerate(islice(expansion, prec_cutoff)):
                tester.assertEqual(c, self.expansion(lift_mode=mode, n=i+v))

            if mode == 'teichmuller':
                q = self.parent().residue_field().cardinality()
                for c in islice(expansion, prec_cutoff):
                    tester.assertEqual(c, c**q)
