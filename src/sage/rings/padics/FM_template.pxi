"""
Fixed modulus template for complete discrete valuation rings

In order to use this template you need to write a linkage file and
gluing file.  For an example see ``mpz_linkage.pxi`` (linkage file) and
``padic_fixed_modulus_element.pyx`` (gluing file).

The linkage file implements a common API that is then used in the
class :class:`FMElement` defined here.  See ``sage/libs/linkages/padics/API.pxi``
for the functions needed.

The gluing file does the following:

- ``ctypedef``'s ``celement`` to be the appropriate type (e.g. ``mpz_t``)
- includes the linkage file
- includes this template
- defines a concrete class inheriting from :class:`FMElement`, and implements
  any desired extra methods

AUTHORS:

- David Roe (2012-03-01) -- initial version
"""

#*****************************************************************************
#       Copyright (C) 2007-2012 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "padic_template_element.pxi"
from cpython.long cimport *

from sage.structure.element cimport Element
from sage.rings.padics.common_conversion cimport comb_prec, _process_args_and_kwds
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.categories.sets_cat import Sets
from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
from sage.categories.homset import Hom

cdef class FMElement(pAdicTemplateElement):
    cdef int _set(self, x, long val, long xprec, absprec, relprec) except -1:
        """
        Set the value of this element from given defining data.

        This function is intended for use in conversion, and should
        not be called on an element created with :meth:`_new_c`.

        INPUT:

        - ``x`` -- data defining a `p`-adic element: int, long,
          Integer, Rational, other `p`-adic element...

        - ``val`` -- the valuation of the resulting element (unused;
          for compatibility with other `p`-adic precision modes)

        - ``xprec -- an inherent precision of ``x`` (unused; for
          compatibility with other `p`-adic precision modes)

        - ``absprec`` -- an absolute precision cap for this element
          (unused; for compatibility with other `p`-adic precision
          modes)

        - ``relprec`` -- a relative precision cap for this element
          (unused; for compatibility with other `p`-adic precision
          modes)

        TESTS::

            sage: R = ZpFM(5)
            sage: a = R(17,5); a  # indirect doctest
            2 + 3*5
            sage: R = ZpFM(5,5)
            sage: a = R(25/9); a  # indirect doctest
            4*5^2 + 2*5^3
        """
        IF CELEMENT_IS_PY_OBJECT:
            polyt = type(self.prime_pow.modulus)
            self.value = <celement>polyt.__new__(polyt)
        cconstruct(self.value, self.prime_pow)
        if isinstance(x, FMElement) and x.parent() is self.parent():
            cshift_notrunc(self.value, (<FMElement>x).value, 0, 0, self.prime_pow, False)
        else:
            cconv(self.value, x, self.prime_pow.ram_prec_cap, 0, self.prime_pow)

    cdef FMElement _new_c(self):
        """
        Create a new element with the same basic info.

        TESTS::

            sage: R = ZpFM(5); R(6) * R(7)  # indirect doctest
            2 + 3*5 + 5^2
        """
        cdef type t = type(self)
        cdef FMElement ans = t.__new__(t)
        ans._parent = self._parent
        ans.prime_pow = self.prime_pow
        IF CELEMENT_IS_PY_OBJECT:
            polyt = type(self.prime_pow.modulus)
            ans.value = <celement>polyt.__new__(polyt)
        cconstruct(ans.value, ans.prime_pow)
        return ans

    cdef pAdicTemplateElement _new_with_value(self, celement value, long absprec):
        """
        Create a new element with a given value and absolute precision.

        Used by code that doesn't know the precision type.
        """
        cdef FMElement ans = self._new_c()
        creduce(ans.value, value, ans.prime_pow.prec_cap, ans.prime_pow)
        return ans

    cdef int _get_unit(self, celement value) except -1:
        """
        Set ``value`` to the unit of this `p`-adic element.
        """
        cremove(value, self.value, self.prime_pow.ram_prec_cap, self.prime_pow)

    cdef int check_preccap(self) except -1:
        """
        Check that the precision of this element does not exceed the
        precision cap. Does nothing for fixed mod elements.

        TESTS::

            sage: ZpFM(5)(1).lift_to_precision(30)  # indirect doctest
            1
        """
        pass

    def __copy__(self):
        """
        Return a copy of this element.

        EXAMPLES::

            sage: a = ZpFM(5,6)(17); b = copy(a)
            sage: a == b
            True
            sage: a is b
            False
        """
        cdef FMElement ans = self._new_c()
        ccopy(ans.value, self.value, ans.prime_pow)
        return ans

    def __dealloc__(self):
        """
        Deallocate the underlying data structure.

        TESTS::

            sage: R = ZpFM(5)
            sage: a = R(17)
            sage: del(a)
        """
        cdestruct(self.value, self.prime_pow)

    def __reduce__(self):
        """
        Return a tuple of a function and data that can be used to unpickle this
        element.

        EXAMPLES::

            sage: a = ZpFM(5)(-3)
            sage: type(a)
            <class 'sage.rings.padics.padic_fixed_mod_element.pAdicFixedModElement'>
            sage: loads(dumps(a)) == a
            True
        """
        return unpickle_fme_v2, (self.__class__, self.parent(), cpickle(self.value, self.prime_pow))

    cpdef _neg_(self):
        r"""
        Return the additive inverse of this element.

        EXAMPLES::

            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: -R(7)  # indirect doctest
            6*7 + 6*7^2 + 6*7^3
        """
        cdef FMElement ans = self._new_c()
        cneg(ans.value, self.value, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        creduce_small(ans.value, ans.value, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        return ans

    cpdef _add_(self, _right):
        r"""
        Return the sum of this element and ``_right``.

        EXAMPLES::

            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: x = R(1721); x
            6 + 5*7^3
            sage: y = R(1373); y
            1 + 4*7^3
            sage: x + y  # indirect doctest
            7 + 2*7^3
        """
        cdef FMElement right = _right
        cdef FMElement ans = self._new_c()
        cadd(ans.value, self.value, right.value, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        creduce_small(ans.value, ans.value, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        return ans

    cpdef _sub_(self, _right):
        r"""
        Return the difference of this element and ``_right``.

        EXAMPLES::

            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: x = R(1721); x
            6 + 5*7^3
            sage: y = R(1373); y
            1 + 4*7^3
            sage: x - y  # indirect doctest
            5 + 7^3
        """
        cdef FMElement right = _right
        cdef FMElement ans = self._new_c()
        csub(ans.value, self.value, right.value, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        creduce_small(ans.value, ans.value, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        return ans

    def __invert__(self):
        r"""
        Return multiplicative inverse of this element. The valuation
        of ``self`` must be zero.

        EXAMPLES::

            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: ~R(2)
            4 + 3*7 + 3*7^2 + 3*7^3
            sage: ~R(0)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
            sage: ~R(7)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
        """
        if not cisunit(self.value, self.prime_pow):
            raise ValueError("cannot invert non-unit")
        cdef FMElement ans = self._new_c()
        cinvert(ans.value, self.value, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        return ans

    cpdef _mul_(self, _right):
        r"""
        Return the product of this element and ``_right``.

        EXAMPLES::

            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: R(3) * R(2)  # indirect doctest
            6
            sage: R(1/2) * R(2)
            1
        """
        cdef FMElement right = _right
        cdef FMElement ans = self._new_c()
        cmul(ans.value, self.value, right.value, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        creduce(ans.value, ans.value, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        return ans

    cpdef _div_(self, _right):
        r"""
        Return the quotient of this element and ``right``. ``right`` must have
        valuation zero.

        EXAMPLES::

            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: R(3) / R(2)  # indirect doctest
            5 + 3*7 + 3*7^2 + 3*7^3
            sage: R(5) / R(0)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
            sage: R(7) / R(49)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
        """
        cdef FMElement right = _right
        cdef FMElement ans = self._new_c()
        if not cisunit(right.value, self.prime_pow):
            raise ValueError("cannot invert non-unit")
        cdivunit(ans.value, self.value, right.value, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        creduce(ans.value, ans.value, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        return ans

    def _quo_rem(self, _right):
        """
        Quotient with remainder.

        EXAMPLES::

            sage: R = ZpFM(3, 5)
            sage: R(12).quo_rem(R(2))  # indirect doctest
            (2*3, 0)
            sage: R(2).quo_rem(R(12))
            (0, 2)
            sage: q, r = R(4).quo_rem(R(12)); q, r
            (1 + 2*3 + 2*3^3, 1)
            sage: 12*q + r == 4
            True
        """
        cdef FMElement right = _right
        if ciszero(right.value, right.prime_pow):
            raise ZeroDivisionError
        cdef FMElement q = self._new_c()
        cdef FMElement r = self._new_c()
        cdef long sval, rval, diff, pcap = self.prime_pow.ram_prec_cap
        sval = self.valuation_c()
        rval = right.valuation_c()
        diff = sval - rval
        if ciszero(self.value, self.prime_pow):
            csetzero(q.value, q.prime_pow)
            csetzero(r.value, r.prime_pow)
        elif diff >= 0:
            # shift right and self by the same power of the uniformizer
            cshift_notrunc(r.value, right.value, -rval, pcap, r.prime_pow, False)
            cshift_notrunc(q.value, self.value, -rval, pcap, q.prime_pow, False)
            # divide
            cdivunit(q.value, q.value, r.value, pcap, q.prime_pow)
            csetzero(r.value, r.prime_pow)
        else:
            cshift(q.value, r.value, self.value, -rval, pcap, q.prime_pow, False)
            cshift_notrunc(q.prime_pow.shift_rem, right.value, -rval, pcap, q.prime_pow, False)
            cdivunit(q.value, q.value, q.prime_pow.shift_rem, pcap, q.prime_pow)
        creduce(q.value, q.value, pcap, q.prime_pow)
        return q, r

    def __pow__(FMElement self, _right, dummy): # NOTE: dummy ignored, always use self.prime_pow.ram_prec_cap
        """
        Exponentiation by an integer

        EXAMPLES::

            sage: R = ZpFM(11, 5)
            sage: R(1/2)^5
            10 + 7*11 + 11^2 + 5*11^3 + 4*11^4
            sage: R(1/32)
            10 + 7*11 + 11^2 + 5*11^3 + 4*11^4
            sage: R(1/2)^5 == R(1/32)
            True
            sage: R(3)^1000  # indirect doctest
            1 + 4*11^2 + 3*11^3 + 7*11^4

        TESTS:

        We check that :issue:`15640` is resolved::

            sage: R(11)^-1
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
        """
        cdef FMElement ans = self._new_c()
        cdef Integer right = Integer(_right)
        if right < 0:
            self = ~self
            mpz_neg(right.value, right.value)
        cpow(ans.value, self.value, right.value, self.prime_pow.ram_prec_cap, self.prime_pow)
        return ans

    cdef pAdicTemplateElement _lshift_c(self, long shift):
        r"""
        Multiply ``self`` by `\pi^{shift}`.

        If shift < -self.valuation(), digits will be truncated.  See
        :meth:`__rshift__` for details.

        EXAMPLES:

        We create a fixed modulus ring::

            sage: R = ZpFM(5, 20); a = R(1000); a
            3*5^3 + 5^4

        Shifting to the right is the same as dividing by a power of
        the uniformizer `\pi` of the `p`-adic ring.::

            sage: a >> 1
            3*5^2 + 5^3

        Shifting to the left is the same as multiplying by a power of
        `\pi`::

            sage: a << 2
            3*5^5 + 5^6
            sage: a*5^2
            3*5^5 + 5^6

        Shifting by a negative integer to the left is the same as
        right shifting by the absolute value::

            sage: a << -3
            3 + 5
            sage: a >> 3
            3 + 5
        """
        if shift < 0:
            return self._rshift_c(-shift)
        elif shift == 0:
            return self
        cdef FMElement ans = self._new_c()
        if shift >= self.prime_pow.ram_prec_cap:
            csetzero(ans.value, ans.prime_pow)
        else:
            cshift_notrunc(ans.value, self.value, shift, ans.prime_pow.ram_prec_cap, ans.prime_pow, True)
        return ans

    cdef pAdicTemplateElement _rshift_c(self, long shift):
        r"""
        Divide by `\pi^{shift}`, and truncate.

        Note that this operation will insert arbitrary digits (in
        practice, currently all zero) in the least significant digits.

        EXAMPLES::

            sage: R = ZpFM(997, 7); a = R(123456878908); a
            964*997 + 572*997^2 + 124*997^3

        Shifting to the right divides by a power of `\pi`, but
        dropping terms with negative valuation::

            sage: a >> 3
            124

        A negative shift multiplies by that power of `\pi`::

            sage: a >> -3
            964*997^4 + 572*997^5 + 124*997^6
        """
        if shift < 0:
            return self._lshift_c(-shift)
        elif shift == 0:
            return self
        cdef FMElement ans = self._new_c()
        if shift >= self.prime_pow.ram_prec_cap:
            csetzero(ans.value, ans.prime_pow)
        else:
            cshift(ans.value, ans.prime_pow.shift_rem, self.value, -shift, ans.prime_pow.ram_prec_cap, ans.prime_pow, True)
        return ans

    def add_bigoh(self, absprec):
        r"""
        Return a new element truncated modulo `\pi^{\mbox{absprec}}`.

        INPUT:

        - ``absprec`` -- integer or infinity

        OUTPUT: a new element truncated modulo `\pi^{\mbox{absprec}}`

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod','series'); a = R(8); a.add_bigoh(1)
            1

        TESTS:

        We handle very large and very small values for ``absprec`` correctly::

            sage: a = R(7)
            sage: a.add_bigoh(2^1000)
            7
            sage: a.add_bigoh(-2^1000)
            0
        """
        cdef long aprec
        if absprec is infinity:
            return self
        if isinstance(absprec, int):
            aprec = absprec
        else:
            if not isinstance(absprec, Integer):
                absprec = Integer(absprec)
            if mpz_sgn((<Integer>absprec).value) == -1:
                return self.parent().fraction_field()(0)
            elif mpz_fits_slong_p((<Integer>absprec).value) == 0:
                return self
            else:
                aprec = mpz_get_si((<Integer>absprec).value)
        if aprec < 0:
            return self.parent().fraction_field()(self, absprec)
        elif aprec >= self.prime_pow.prec_cap:
            return self
        cdef FMElement ans = self._new_c()
        creduce(ans.value, self.value, aprec, ans.prime_pow)
        return ans

    cpdef bint _is_exact_zero(self) except -1:
        """
        Test whether this element is an exact zero, which is always
        False for fixed modulus elements.

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod','series'); a = R(8); a._is_exact_zero()
            False
        """
        return False

    cpdef bint _is_inexact_zero(self) except -1:
        """
        Return ``True`` if ``self`` is indistinguishable from zero.

        EXAMPLES::

            sage: R = ZpFM(7, 5)
            sage: R(14)._is_inexact_zero()
            False
            sage: R(0)._is_inexact_zero()
            True
        """
        return ciszero(self.value, self.prime_pow)

    def is_zero(self, absprec=None):
        r"""
        Return whether ``self`` is zero modulo `\pi^{\mbox{absprec}}`.

        INPUT:

        - ``absprec`` -- integer

        EXAMPLES::

            sage: R = ZpFM(17, 6)
            sage: R(0).is_zero()
            True
            sage: R(17^6).is_zero()
            True
            sage: R(17^2).is_zero(absprec=2)
            True
        """
        cdef bint iszero = ciszero(self.value, self.prime_pow)
        if absprec is None:
            return iszero
        if not isinstance(absprec, Integer):
            absprec = Integer(absprec)
        if mpz_cmp_si((<Integer>absprec).value, self.prime_pow.ram_prec_cap) >= 0:
            return iszero
        cdef long val = self.valuation_c()
        return mpz_cmp_si((<Integer>absprec).value, val) <= 0

    def __bool__(self):
        """
        Return ``True`` if this element is distinguishable from zero.

        For most applications, explicitly specifying the power of p
        modulo which the element is supposed to be nonzero is preferable.

        EXAMPLES::

            sage: R = ZpFM(5); a = R(0); b = R(75)
            sage: bool(a), bool(b)  # indirect doctest
            (False, True)
        """
        return not ciszero(self.value, self.prime_pow)

    def is_equal_to(self, _right, absprec=None):
        r"""
        Return whether this element is equal to ``right`` modulo `p^{\mbox{absprec}}`.

        If ``absprec`` is ``None``, returns if ``self == 0``.

        INPUT:

        - ``right`` -- a `p`-adic element with the same parent
        - ``absprec`` -- positive integer or ``None`` (default: ``None``)

        EXAMPLES::

            sage: R = ZpFM(2, 6)
            sage: R(13).is_equal_to(R(13))
            True
            sage: R(13).is_equal_to(R(13+2^10))
            True
            sage: R(13).is_equal_to(R(17), 2)
            True
            sage: R(13).is_equal_to(R(17), 5)
            False
        """
        cdef FMElement right
        cdef long aprec, rprec, sval, rval
        if self.parent() is _right.parent():
            right = _right
        else:
            right = self.parent()(_right)
        if absprec is None:
            # The default absolute precision is given by the precision cap
            aprec = self.prime_pow.ram_prec_cap
        else:
            if not isinstance(absprec, Integer):
                absprec = Integer(absprec)
            # If absprec is not positive, then self and right are always
            # equal.
            if mpz_sgn((<Integer>absprec).value) < 0:
                return True
            # If absprec is bigger than the precision cap, we use it
            # instead.
            if mpz_cmp_si((<Integer>absprec).value, self.prime_pow.ram_prec_cap) >= 0:
                aprec = self.prime_pow.ram_prec_cap
            else:
                aprec = mpz_get_si((<Integer>absprec).value)
        return ccmp(self.value,
                    right.value,
                    aprec,
                    aprec < self.prime_pow.ram_prec_cap,
                    aprec < right.prime_pow.ram_prec_cap,
                    self.prime_pow) == 0

    cdef int _cmp_units(self, pAdicGenericElement _right) except -2:
        """
        Comparison of units, used in equality testing.

        EXAMPLES::

            sage: R = ZpFM(5)
            sage: a = R(17); b = R(0,3); c = R(85,7); d = R(2, 1)
            sage: any([a == b, a == c, b == c, b == d, c == d, a == d])  # indirect doctest
            False
            sage: all([a == a, b == b, c == c, d == d])
            True
        """
        cdef FMElement right = _right
        return ccmp(self.value, right.value, self.prime_pow.ram_prec_cap, False, False, self.prime_pow)

    cdef pAdicTemplateElement lift_to_precision_c(self, long absprec):
        """
        Lifts this element to another with precision at least absprec.

        Since fixed modulus elements don't track precision, this
        function just returns the same element.

        EXAMPLES::

            sage: R = ZpFM(5)
            sage: a = R(77, 2); a
            2 + 3*5^2
            sage: a.lift_to_precision(17)  # indirect doctest
            2 + 3*5^2
        """
        return self

    def _teichmuller_set_unsafe(self):
        """
        Set this element to the Teichmuller representative with the
        same residue.

        .. WARNING::

            This function modifies the element, which is not safe.
            Elements are supposed to be immutable.

        EXAMPLES::

            sage: R = ZpFM(17,5); a = R(11)
            sage: a
            11
            sage: a._teichmuller_set_unsafe(); a
            11 + 14*17 + 2*17^2 + 12*17^3 + 15*17^4
            sage: E = a.expansion(lift_mode='teichmuller'); E
            17-adic expansion of 11 + 14*17 + 2*17^2 + 12*17^3 + 15*17^4 (teichmuller)
            sage: list(E)
            [11 + 14*17 + 2*17^2 + 12*17^3 + 15*17^4, 0, 0, 0, 0]

        Note that if you set an element which is congruent to 0 you
        get 0 to maximum precision::

            sage: b = R(17*5); b
            5*17
            sage: b._teichmuller_set_unsafe(); b
            0
        """
        if cisunit(self.value, self.prime_pow):
            cteichmuller(self.value, self.value, self.prime_pow.ram_prec_cap, self.prime_pow)
        else:
            csetzero(self.value, self.prime_pow)

    def _polynomial_list(self, pad=False):
        """
        Return the coefficient list for a polynomial over the base ring
        yielding this element.

        INPUT:

        - ``pad`` -- whether to pad the result with zeros of the appropriate precision

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<x> = ZZ[]
            sage: K.<a> = ZqFM(25)
            sage: W.<w> = K.extension(x^3 - 5)
            sage: (1 + w)._polynomial_list()
            [1, 1]
            sage: (1 + w)._polynomial_list(pad=True)
            [1, 1, 0]
        """
        R = self.base_ring()
        e = self.parent().relative_e()
        L = ccoefficients(self.value, 0, self.prime_pow.ram_prec_cap, self.prime_pow)
        if pad:
            n = self.parent().relative_degree()
            L.extend([R.zero()] * (n - len(L)))
        return L

    def polynomial(self, var='x'):
        """
        Return a polynomial over the base ring that yields this element
        when evaluated at the generator of the parent.

        INPUT:

        - ``var`` -- string, the variable name for the polynomial

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFM(5^3)
            sage: a.polynomial()
            x
            sage: a.polynomial(var='y')
            y
            sage: (5*a^2 + 25).polynomial()
            5*x^2 + 5^2
        """
        R = self.base_ring()
        S = R[var]
        return S(self._polynomial_list())

    def precision_absolute(self):
        """
        The absolute precision of this element.

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod'); a = R(7); a.precision_absolute()
            4
        """
        cdef Integer ans = Integer.__new__(Integer)
        mpz_set_si(ans.value, self.prime_pow.ram_prec_cap)
        return ans

    def precision_relative(self):
        r"""
        The relative precision of this element.

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod'); a = R(7); a.precision_relative()
            3
            sage: a = R(0); a.precision_relative()
            0
        """
        cdef Integer ans = Integer.__new__(Integer)
        mpz_set_si(ans.value, self.prime_pow.ram_prec_cap - self.valuation_c())
        return ans

    cpdef pAdicTemplateElement unit_part(FMElement self):
        r"""
        Return the unit part of ``self``.

        If the valuation of ``self`` is positive, then the high digits of the
        result will be zero.

        EXAMPLES::

            sage: R = Zp(17, 4, 'fixed-mod')
            sage: R(5).unit_part()
            5
            sage: R(18*17).unit_part()
            1 + 17
            sage: R(0).unit_part()
            0
            sage: type(R(5).unit_part())
            <class 'sage.rings.padics.padic_fixed_mod_element.pAdicFixedModElement'>
            sage: R = ZpFM(5, 5); a = R(75); a.unit_part()
            3
        """
        cdef FMElement ans = (<FMElement>self)._new_c()
        cremove(ans.value, (<FMElement>self).value, (<FMElement>self).prime_pow.ram_prec_cap, (<FMElement>self).prime_pow)
        return ans

    cdef long valuation_c(self) noexcept:
        """
        Return the valuation of this element.

        TESTS::

            sage: R = ZpFM(5, 5); R(0).valuation()  # indirect doctest
            5
            sage: R = Zp(17, 4,'fixed-mod')
            sage: a = R(2*17^2)
            sage: a.valuation()
            2
            sage: R = Zp(5, 4,'fixed-mod')
            sage: R(0).valuation()
            4
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
        """
        # for backward compatibility
        return cvaluation(self.value, self.prime_pow.ram_prec_cap, self.prime_pow)

    cpdef val_unit(self):
        """
        Return a 2-tuple, the first element set to the valuation of
        ``self``, and the second to the unit part of ``self``.

        If ``self == 0``, then the unit part is ``O(p^self.parent().precision_cap())``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: a = R(75); b = a - a
            sage: a.val_unit()
            (2, 3)
            sage: b.val_unit()
            (5, 0)
        """
        cdef FMElement unit = self._new_c()
        cdef Integer valuation = Integer.__new__(Integer)
        mpz_set_si(valuation.value, cremove(unit.value, self.value, self.prime_pow.ram_prec_cap, self.prime_pow))
        return valuation, unit

    def __hash__(self):
        """
        Hashing.

        EXAMPLES::

            sage: R = ZpFM(11, 5)
            sage: hash(R(3)) == hash(3)
            True
        """
        return chash(self.value, 0, self.prime_pow.ram_prec_cap, self.prime_pow)

cdef class pAdicCoercion_ZZ_FM(RingHomomorphism):
    r"""
    The canonical inclusion from `\ZZ` to a fixed modulus ring.

    EXAMPLES::

        sage: f = ZpFM(5).coerce_map_from(ZZ); f
        Ring morphism:
          From: Integer Ring
          To:   5-adic Ring of fixed modulus 5^20

    TESTS::

        sage: TestSuite(f).run()
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = ZpFM(5).coerce_map_from(ZZ); type(f)
            <class 'sage.rings.padics.padic_fixed_mod_element.pAdicCoercion_ZZ_FM'>
        """
        RingHomomorphism.__init__(self, ZZ.Hom(R))
        self._zero = R.element_class(R, 0)
        self._section = pAdicConvert_FM_ZZ(R)

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpFM(5).coerce_map_from(ZZ)
            sage: g = copy(f)  # indirect doctest
            sage: g == f
            True
            sage: g(6)
            1 + 5
            sage: g(6) == f(6)
            True
        """
        _slots = RingHomomorphism._extra_slots(self)
        _slots['_zero'] = self._zero
        _slots['_section'] = self.section() # use method since it copies coercion-internal sections.
        return _slots

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpFM(5).coerce_map_from(ZZ)
            sage: g = copy(f)  # indirect doctest
            sage: g == f
            True
            sage: g(6)
            1 + 5
            sage: g(6) == f(6)
            True
        """
        self._zero = _slots['_zero']
        self._section = _slots['_section']
        RingHomomorphism._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = ZpFM(5).coerce_map_from(ZZ)
            sage: f(0).parent()
            5-adic Ring of fixed modulus 5^20
            sage: f(5)
            5
        """
        if mpz_sgn((<Integer>x).value) == 0:
            return self._zero
        cdef FMElement ans = self._zero._new_c()
        cconv_mpz_t(ans.value, (<Integer>x).value, ans.prime_pow.ram_prec_cap, True, ans.prime_pow)
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in (relative or absolute or both).

        INPUT:

        - ``x`` -- an Integer

        - ``absprec``, or the first positional argument -- the maximum
          absolute precision (unused for fixed modulus elements)

        - ``relprec``, or the second positional argument -- the
          maximum relative precision (unused for fixed modulus
          elements)

        EXAMPLES::

            sage: R = ZpFM(5,4)
            sage: type(R(10,2))
            <class 'sage.rings.padics.padic_fixed_mod_element.pAdicFixedModElement'>
            sage: R(30,2)
            5 + 5^2
            sage: R(30,3,1)
            5 + 5^2
            sage: R(30,absprec=2)
            5 + 5^2
            sage: R(30,relprec=2)
            5 + 5^2
            sage: R(30,absprec=1)
            5 + 5^2
            sage: R(30,empty=True)
            5 + 5^2
        """
        if mpz_sgn((<Integer>x).value) == 0:
            return self._zero
        cdef FMElement ans = self._zero._new_c()
        cconv_mpz_t(ans.value, (<Integer>x).value, ans.prime_pow.ram_prec_cap, True, ans.prime_pow)
        return ans

    def section(self):
        r"""
        Return a map back to `\ZZ` that approximates an element of this
        `p`-adic ring by an integer.

        EXAMPLES::

            sage: f = ZpFM(5).coerce_map_from(ZZ).section()
            sage: f(ZpFM(5)(-1)) - 5^20
            -1
        """
        from sage.misc.constant_function import ConstantFunction
        if not isinstance(self._section.domain, ConstantFunction):
            import copy
            self._section = copy.copy(self._section)
        return self._section

cdef class pAdicConvert_FM_ZZ(RingMap):
    r"""
    The map from a fixed modulus ring back to `\ZZ` that returns the smallest
    nonnegative integer approximation to its input which is accurate up to the precision.

    If the input is not in the closure of the image of `\ZZ`, raises a :exc:`ValueError`.

    EXAMPLES::

        sage: f = ZpFM(5).coerce_map_from(ZZ).section(); f
        Set-theoretic ring morphism:
          From: 5-adic Ring of fixed modulus 5^20
          To:   Integer Ring
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = ZpFM(5).coerce_map_from(ZZ).section(); type(f)
            <class 'sage.rings.padics.padic_fixed_mod_element.pAdicConvert_FM_ZZ'>
            sage: f.category()
            Category of homsets of sets
        """
        if R.absolute_degree() > 1 or R.characteristic() != 0 or R.residue_characteristic() == 0:
            RingMap.__init__(self, Hom(R, ZZ, SetsWithPartialMaps()))
        else:
            RingMap.__init__(self, Hom(R, ZZ, Sets()))

    cpdef Element _call_(self, _x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = ZpFM(5).coerce_map_from(ZZ).section()
            sage: f(ZpFM(5)(-1)) - 5^20
            -1
            sage: f(ZpFM(5)(0))
            0
        """
        cdef Integer ans = Integer.__new__(Integer)
        cdef FMElement x = _x
        cconv_mpz_t_out(ans.value, x.value, 0, x.prime_pow.ram_prec_cap, x.prime_pow)
        return ans

cdef class pAdicConvert_QQ_FM(Morphism):
    r"""
    The inclusion map from `\QQ` to a fixed modulus ring that is defined
    on all elements with nonnegative `p`-adic valuation.

    EXAMPLES::

        sage: f = ZpFM(5).convert_map_from(QQ); f
        Generic morphism:
          From: Rational Field
          To:   5-adic Ring of fixed modulus 5^20
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = ZpFM(5).convert_map_from(QQ); type(f)
            <class 'sage.rings.padics.padic_fixed_mod_element.pAdicConvert_QQ_FM'>
        """
        Morphism.__init__(self, Hom(QQ, R, SetsWithPartialMaps()))
        self._zero = R.element_class(R, 0)

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpFM(5).convert_map_from(QQ)
            sage: g = copy(f)  # indirect doctest
            sage: g == f # todo: comparison not implemented
            True
            sage: g(1/6)
            1 + 4*5 + 4*5^3 + 4*5^5 + 4*5^7 + 4*5^9 + 4*5^11 + 4*5^13 + 4*5^15 + 4*5^17 + 4*5^19
            sage: g(1/6) == f(1/6)
            True
        """
        _slots = Morphism._extra_slots(self)
        _slots['_zero'] = self._zero
        return _slots

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpFM(5).convert_map_from(QQ)
            sage: g = copy(f)  # indirect doctest
            sage: g == f # todo: comparison not implemented
            True
            sage: g(1/6)
            1 + 4*5 + 4*5^3 + 4*5^5 + 4*5^7 + 4*5^9 + 4*5^11 + 4*5^13 + 4*5^15 + 4*5^17 + 4*5^19
            sage: g(1/6) == f(1/6)
            True
        """
        self._zero = _slots['_zero']
        Morphism._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = ZpFM(5,4).convert_map_from(QQ)
            sage: f(1/7)
            3 + 3*5 + 2*5^3
            sage: f(0)
            0
        """
        if mpq_sgn((<Rational>x).value) == 0:
            return self._zero
        cdef FMElement ans = self._zero._new_c()
        cconv_mpq_t(ans.value, (<Rational>x).value, ans.prime_pow.ram_prec_cap, True, ans.prime_pow)
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in (relative or absolute or both).

        INPUT:

        - ``x`` -- a Rational

        - ``absprec``, or the first positional argument -- the maximum
          absolute precision (unused for fixed modulus elements)

        - ``relprec``, or the second positional argument -- the
          maximum relative precision (unused for fixed modulus
          elements)

        EXAMPLES::

            sage: R = ZpFM(5,4)
            sage: type(R(1/7,2))
            <class 'sage.rings.padics.padic_fixed_mod_element.pAdicFixedModElement'>
            sage: R(1/7,2)
            3 + 3*5 + 2*5^3
            sage: R(1/7,3,1)
            3 + 3*5 + 2*5^3
            sage: R(1/7,absprec=2)
            3 + 3*5 + 2*5^3
            sage: R(1/7,relprec=2)
            3 + 3*5 + 2*5^3
            sage: R(1/7,absprec=1)
            3 + 3*5 + 2*5^3
            sage: R(1/7,empty=True)
            3 + 3*5 + 2*5^3
        """
        if mpq_sgn((<Rational>x).value) == 0:
            return self._zero
        cdef FMElement ans = self._zero._new_c()
        cconv_mpq_t(ans.value, (<Rational>x).value, ans.prime_pow.ram_prec_cap, True, ans.prime_pow)
        return ans

cdef class pAdicCoercion_FM_frac_field(RingHomomorphism):
    r"""
    The canonical inclusion of `\ZZ_q` into its fraction field.

    EXAMPLES::

        sage: # needs sage.libs.flint
        sage: R.<a> = ZqFM(27, implementation='FLINT')
        sage: K = R.fraction_field()
        sage: f = K.coerce_map_from(R); f
        Ring morphism:
          From: 3-adic Unramified Extension Ring in a defined by x^3 + 2*x + 1
          To:   3-adic Unramified Extension Field in a defined by x^3 + 2*x + 1

    TESTS::

        sage: TestSuite(f).run()                                                        # needs sage.libs.flint
    """
    def __init__(self, R, K):
        """
        Initialization.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFM(27)
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R); type(f)
            <class 'sage.rings.padics.qadic_flint_FM.pAdicCoercion_FM_frac_field'>
        """
        RingHomomorphism.__init__(self, R.Hom(K))
        self._zero = K(0)
        self._section = pAdicConvert_FM_frac_field(K, R)

    cpdef Element _call_(self, _x):
        """
        Evaluation.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFM(27)
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f(a)
            a
        """
        cdef FMElement x = _x
        if ciszero(x.value, x.prime_pow):
            return self._zero
        cdef FPElement ans = self._zero._new_c()
        ans.ordp = cremove(ans.unit, x.value, x.prime_pow.ram_prec_cap, x.prime_pow)
        IF CELEMENT_IS_PY_OBJECT:
            # The base ring is wrong, so we fix it.
            K = ans.unit.base_ring()
            ans.unit._coeffs = [K(c) for c in ans.unit._coeffs]
        return ans

    cpdef Element _call_with_args(self, _x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in
        (relative or absolute or both).

        See the documentation for
        :meth:`pAdicCappedAbsoluteElement.__init__` for more details.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFM(27)
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f(a, 3)
            a
            sage: b = 117*a
            sage: f(b, 3)
            a*3^2
            sage: f(b, 4, 1)
            a*3^2
            sage: f(b, 4, 3)
            a*3^2 + a*3^3
            sage: f(b, absprec=4)
            a*3^2 + a*3^3
            sage: f(b, relprec=3)
            a*3^2 + a*3^3 + a*3^4
            sage: f(b, absprec=1)
            0
            sage: f(R(0))
            0
        """
        cdef long aprec, rprec
        cdef FMElement x = _x
        if ciszero(x.value, x.prime_pow):
            return self._zero
        cdef FPElement ans = self._zero._new_c()
        cdef bint reduce = False
        _process_args_and_kwds(&aprec, &rprec, args, kwds, False, x.prime_pow)
        ans.ordp = cremove(ans.unit, x.value, aprec, x.prime_pow)
        if aprec < ans.ordp + rprec:
            rprec = aprec - ans.ordp
        if rprec <= 0:
            return self._zero
        creduce(ans.unit, ans.unit, rprec, x.prime_pow)
        IF CELEMENT_IS_PY_OBJECT:
            # The base ring is wrong, so we fix it.
            K = ans.unit.base_ring()
            ans.unit._coeffs = [K(c) for c in ans.unit._coeffs]
        return ans

    def section(self):
        """
        Return a map back to the ring that converts elements of
        nonnegative valuation.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFM(27)
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f.section()(K.gen())
            a
        """
        from sage.misc.constant_function import ConstantFunction
        if not isinstance(self._section.domain, ConstantFunction):
            import copy
            self._section = copy.copy(self._section)
        return self._section

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        TESTS::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFM(27)
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Ring morphism:
              From: 3-adic Unramified Extension Ring in a defined by x^3 + 2*x + 1
              To:   3-adic Unramified Extension Field in a defined by x^3 + 2*x + 1
            sage: g == f
            True
            sage: g is f
            False
            sage: g(a)
            a
            sage: g(a) == f(a)
            True
        """
        _slots = RingHomomorphism._extra_slots(self)
        _slots['_zero'] = self._zero
        _slots['_section'] = self.section() # use method since it copies coercion-internal sections.
        return _slots

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        TESTS::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFM(9)
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Ring morphism:
              From: 3-adic Unramified Extension Ring in a defined by x^2 + 2*x + 2
              To:   3-adic Unramified Extension Field in a defined by x^2 + 2*x + 2
            sage: g == f
            True
            sage: g is f
            False
            sage: g(a)
            a
            sage: g(a) == f(a)
            True
        """
        self._zero = _slots['_zero']
        self._section = _slots['_section']
        RingHomomorphism._update_slots(self, _slots)

    def is_injective(self):
        r"""
        Return whether this map is injective.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFM(9)
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f.is_injective()
            True
        """
        return True

    def is_surjective(self):
        r"""
        Return whether this map is surjective.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFM(9)
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f.is_surjective()
            False
        """
        return False


cdef class pAdicConvert_FM_frac_field(Morphism):
    r"""
    The section of the inclusion from `\ZZ_q` to its fraction field.

    EXAMPLES::

        sage: # needs sage.libs.flint
        sage: R.<a> = ZqFM(27)
        sage: K = R.fraction_field()
        sage: f = R.convert_map_from(K); f
        Generic morphism:
          From: 3-adic Unramified Extension Field in a defined by x^3 + 2*x + 1
          To:   3-adic Unramified Extension Ring in a defined by x^3 + 2*x + 1
    """
    def __init__(self, K, R):
        """
        Initialization.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFM(27)
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K); type(f)
            <class 'sage.rings.padics.qadic_flint_FM.pAdicConvert_FM_frac_field'>
        """
        Morphism.__init__(self, Hom(K, R, SetsWithPartialMaps()))
        self._zero = R(0)

    cpdef Element _call_(self, _x):
        """
        Evaluation.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFM(27)
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K)
            sage: f(K.gen())
            a
        """
        cdef FPElement x = _x
        if x.ordp < 0:
            raise ValueError("negative valuation")
        if x.ordp >= self._zero.prime_pow.ram_prec_cap:
            return self._zero
        cdef FMElement ans = self._zero._new_c()
        cshift_notrunc(ans.value, x.unit, x.ordp, ans.prime_pow.ram_prec_cap, ans.prime_pow, x.ordp > 0)
        IF CELEMENT_IS_PY_OBJECT:
            # The base ring is wrong, so we fix it.
            R = ans.value.base_ring()
            ans.value._coeffs = [R(c) for c in ans.value._coeffs]
        return ans

    cpdef Element _call_with_args(self, _x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in
        (relative or absolute or both).

        See the documentation for
        :meth:`pAdicCappedAbsoluteElement.__init__` for more details.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFM(27)
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K); a = K(a)
            sage: f(a, 3)
            a
            sage: b = 117*a
            sage: f(b, 3)
            a*3^2
            sage: f(b, 4, 1)
            a*3^2
            sage: f(b, 4, 3)
            a*3^2 + a*3^3
            sage: f(b, absprec=4)
            a*3^2 + a*3^3
            sage: f(b, relprec=3)
            a*3^2 + a*3^3 + a*3^4
            sage: f(b, absprec=1)
            0
            sage: f(K(0))
            0
        """
        cdef long aprec, rprec
        cdef FPElement x = _x
        if x.ordp < 0:
            raise ValueError("negative valuation")
        if x.ordp >= self._zero.prime_pow.ram_prec_cap:
            return self._zero
        cdef FMElement ans = self._zero._new_c()
        _process_args_and_kwds(&aprec, &rprec, args, kwds, True, ans.prime_pow)
        if rprec < aprec - x.ordp:
            aprec = x.ordp + rprec
        cshift_notrunc(ans.value, x.unit, x.ordp, aprec, ans.prime_pow, x.ordp > 0)
        IF CELEMENT_IS_PY_OBJECT:
            # The base ring is wrong, so we fix it.
            R = ans.value.base_ring()
            ans.value._coeffs = [R(c) for c in ans.value._coeffs]
        return ans

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        TESTS::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFM(27)
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K)
            sage: a = K(a)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Generic morphism:
             From: 3-adic Unramified Extension Field in a defined by x^3 + 2*x + 1
             To:   3-adic Unramified Extension Ring in a defined by x^3 + 2*x + 1
            sage: g == f
            True
            sage: g is f
            False
            sage: g(a)
            a
            sage: g(a) == f(a)
            True
        """
        _slots = Morphism._extra_slots(self)
        _slots['_zero'] = self._zero
        return _slots

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        TESTS::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFM(9)
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K)
            sage: a = f(a)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Generic morphism:
              From: 3-adic Unramified Extension Field in a defined by x^2 + 2*x + 2
              To:   3-adic Unramified Extension Ring in a defined by x^2 + 2*x + 2
            sage: g == f
            True
            sage: g is f
            False
            sage: g(a)
            a
            sage: g(a) == f(a)
            True
        """
        self._zero = _slots['_zero']
        Morphism._update_slots(self, _slots)


def unpickle_fme_v2(cls, parent, value):
    """
    Unpickles a fixed-mod element.

    EXAMPLES::

        sage: from sage.rings.padics.padic_fixed_mod_element import pAdicFixedModElement, unpickle_fme_v2
        sage: R = ZpFM(5)
        sage: a = unpickle_fme_v2(pAdicFixedModElement, R, 17*25); a
        2*5^2 + 3*5^3
        sage: a.parent() is R
        True
    """
    cdef FMElement ans = cls.__new__(cls)
    ans._parent = parent
    ans.prime_pow = <PowComputer_?>parent.prime_pow
    IF CELEMENT_IS_PY_OBJECT:
        polyt = type(ans.prime_pow.modulus)
        ans.value = <celement>polyt.__new__(polyt)
    cconstruct(ans.value, ans.prime_pow)
    cunpickle(ans.value, value, ans.prime_pow)
    return ans
