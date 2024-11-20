r"""
Capped absolute template for complete discrete valuation rings

In order to use this template you need to write a linkage file and gluing file.
For an example see ``mpz_linkage.pxi`` (linkage file) and
``padic_capped_absolute_element.pyx`` (gluing file).

The linkage file implements a common API that is then used in the class
:class:`CAElement` defined here.
See the documentation of ``mpz_linkage.pxi`` for the functions needed.

The gluing file does the following:

- ctypedef's celement to be the appropriate type (e.g. ``mpz_t``)
- includes the linkage file
- includes this template
- defines a concrete class inheriting from :class:`CAElement`, and implements
  any desired extra methods

AUTHORS:

- David Roe (2012-3-1) -- initial version
"""

#*****************************************************************************
#       Copyright (C) 2012 David Roe <roed.math@gmail.com>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# This file implements common functionality among template elements
include "padic_template_element.pxi"

from collections.abc import Iterable
from sage.structure.element cimport Element
from sage.rings.padics.common_conversion cimport comb_prec, _process_args_and_kwds
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.categories.sets_cat import Sets
from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
from sage.categories.homset import Hom

cdef class CAElement(pAdicTemplateElement):
    cdef int _set(self, x, long val, long xprec, absprec, relprec) except -1:
        """
        Set the value of this element from given defining data.

        This function is intended for use in conversion, and should
        not be called on an element created with :meth:`_new_c`.

        INPUT:

        - ``x`` -- data defining a `p`-adic element: int,
          Integer, Rational, other `p`-adic element...

        - ``val`` -- the valuation of the resulting element

        - ``xprec -- an inherent precision of ``x``

        - ``absprec`` -- an absolute precision cap for this element

        - ``relprec`` -- a relative precision cap for this element

        TESTS::

            sage: R = ZpCA(5)
            sage: a = R(17,5); a  # indirect doctest
            2 + 3*5 + O(5^5)
            sage: a = R(75, absprec = 5, relprec = 4); a  # indirect doctest
            3*5^2 + O(5^5)
            sage: a = R(25/9, absprec = 5); a  # indirect doctest
            4*5^2 + 2*5^3 + O(5^5)
            sage: a = R(25/9, absprec = 5, relprec = 4); a  # indirect doctest
            4*5^2 + 2*5^3 + O(5^5)
        """
        IF CELEMENT_IS_PY_OBJECT:
            polyt = type(self.prime_pow.modulus)
            self.value = <celement>polyt.__new__(polyt)
        cconstruct(self.value, self.prime_pow)
        cdef long rprec = comb_prec(relprec, self.prime_pow.ram_prec_cap)
        cdef long aprec = comb_prec(absprec, min(self.prime_pow.ram_prec_cap, xprec))
        if aprec <= val:
            csetzero(self.value, self.prime_pow)
            self.absprec = aprec
        else:
            self.absprec = min(aprec, val + rprec)
            if isinstance(x, CAElement) and x.parent() is self.parent():
                cshift_notrunc(self.value, (<CAElement>x).value, 0, self.absprec, self.prime_pow, True)
            else:
                cconv(self.value, x, self.absprec, 0, self.prime_pow)

    cdef CAElement _new_c(self):
        """
        Create a new element with the same basic info.

        TESTS::

            sage: R = ZpCA(5); R(6,5) * R(7,8)  # indirect doctest
            2 + 3*5 + 5^2 + O(5^5)

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqCA(25)
            sage: S.<x> = ZZ[]
            sage: W.<w> = R.ext(x^2 - 5)
            sage: w * (w+1)  # indirect doctest
            w + w^2 + O(w^40)
        """
        cdef type t = type(self)
        cdef CAElement ans = t.__new__(t)
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
        cdef CAElement ans = self._new_c()
        ans.absprec = absprec
        self.check_preccap()
        creduce(ans.value, value, absprec, ans.prime_pow)
        return ans

    cdef int _get_unit(self, celement value) except -1:
        """
        Set ``value`` to the unit of this `p`-adic element.
        """
        cremove(value, self.value, self.absprec, self.prime_pow, True)

    cdef int check_preccap(self) except -1:
        """
        Check that this element doesn't have precision higher than
        allowed by the precision cap.

        TESTS::

            sage: ZpCA(5)(1).lift_to_precision(30)  # indirect doctest
            Traceback (most recent call last):
            ...
            PrecisionError: precision higher than allowed by the precision cap
        """
        if self.absprec > self.prime_pow.ram_prec_cap:
            raise PrecisionError("precision higher than allowed by the precision cap")

    def __copy__(self):
        """
        Return a copy of this element.

        EXAMPLES::

            sage: a = ZpCA(5,6)(17); b = copy(a)
            sage: a == b
            True
            sage: a is b
            False
        """
        cdef CAElement ans = self._new_c()
        ans.absprec = self.absprec
        ccopy(ans.value, self.value, ans.prime_pow)
        return ans

    def __dealloc__(self):
        """
        Deallocate the underlying data structure.

        TESTS::

            sage: R = ZpCA(5)
            sage: a = R(17)
            sage: del(a)
        """
        cdestruct(self.value, self.prime_pow)

    def __reduce__(self):
        """
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: a = ZpCA(5)(-3)
            sage: type(a)
            <class 'sage.rings.padics.padic_capped_absolute_element.pAdicCappedAbsoluteElement'>
            sage: loads(dumps(a)) == a
            True
        """
        return unpickle_cae_v2, (self.__class__, self.parent(), cpickle(self.value, self.prime_pow), self.absprec)

    cpdef _neg_(self):
        """
        Return the additive inverse of this element.

        EXAMPLES::

            sage: R = Zp(5, prec=10, type='capped-abs')
            sage: a = R(1)
            sage: -a  # indirect doctest
            4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + O(5^10)
        """
        cdef CAElement ans = self._new_c()
        ans.absprec = self.absprec
        cneg(ans.value, self.value, ans.absprec, ans.prime_pow)
        creduce_small(ans.value, ans.value, ans.absprec, ans.prime_pow)
        return ans

    cpdef _add_(self, _right):
        """
        Return the sum of this element and ``_right``.

        EXAMPLES::

            sage: R = ZpCA(13, 4)
            sage: R(2) + R(3)  # indirect doctest
            5 + O(13^4)
            sage: R(12) + R(1)
            13 + O(13^4)

        Check that :issue:`20245` is resolved::

            sage: R(1,1) + R(169,3)
            1 + O(13)
        """
        cdef CAElement right = _right
        cdef CAElement ans = self._new_c()
        ans.absprec = min(self.absprec, right.absprec)
        cadd(ans.value, self.value, right.value, ans.absprec, ans.prime_pow)
        creduce(ans.value, ans.value, ans.absprec, ans.prime_pow)
        return ans

    cpdef _sub_(self, _right):
        """
        Return the difference of this element and ``_right``.

        EXAMPLES::

            sage: R = ZpCA(13, 4)
            sage: R(10) - R(10)  # indirect doctest
            O(13^4)
            sage: R(10) - R(11)
            12 + 12*13 + 12*13^2 + 12*13^3 + O(13^4)
        """
        cdef CAElement right = _right
        cdef CAElement ans = self._new_c()
        ans.absprec = min(self.absprec, right.absprec)
        csub(ans.value, self.value, right.value, ans.absprec, ans.prime_pow)
        creduce(ans.value, ans.value, ans.absprec, ans.prime_pow)
        return ans

    def __invert__(self):
        """
        Return the multiplicative inverse of this element.

        .. NOTE::

            The result always lives in the fraction field, even if this element
            is a unit.

        EXAMPLES::

            sage: R = ZpCA(17)
            sage: ~R(-1) == R(-1)
            True
            sage: ~R(5) * 5
            1 + O(17^20)
            sage: ~R(5)
            7 + 3*17 + 10*17^2 + 13*17^3 + 6*17^4 + 3*17^5 + 10*17^6 + 13*17^7
             + 6*17^8 + 3*17^9 + 10*17^10 + 13*17^11 + 6*17^12 + 3*17^13
             + 10*17^14 + 13*17^15 + 6*17^16 + 3*17^17 + 10*17^18 + 13*17^19 + O(17^20)
            sage: ~R(-1) == R(-1)  # indirect doctest
            True
        """
        return ~self.parent().fraction_field()(self)

    cpdef _mul_(self, _right):
        """
        Return the product of this element and ``_right``.

        EXAMPLES::

            sage: R = ZpCA(5)
            sage: a = R(20,5); b = R(75, 4); a * b  # indirect doctest
            2*5^3 + 2*5^4 + O(5^5)
        """
        cdef CAElement right = _right
        cdef CAElement ans = self._new_c()
        cdef long vals, valr
        if self.absprec == self.prime_pow.ram_prec_cap and right.absprec == self.prime_pow.ram_prec_cap:
            ans.absprec = self.absprec
        else:
            vals = self.valuation_c()
            valr = right.valuation_c()
            ans.absprec = min(vals + valr + min(self.absprec - vals, right.absprec - valr), self.prime_pow.ram_prec_cap)
        cmul(ans.value, self.value, right.value, ans.absprec, ans.prime_pow)
        creduce(ans.value, ans.value, ans.absprec, ans.prime_pow)
        return ans

    cpdef _div_(self, right):
        """
        Return the quotient of this element and ``right``.

        .. NOTE::

            The result always lives in the fraction field, even if ``right`` is
            a unit.

        EXAMPLES::

            sage: R = ZpCA(13, 4)
            sage: R(2) / R(3)  # indirect doctest
            5 + 4*13 + 4*13^2 + 4*13^3 + O(13^4)
            sage: a = R(169 * 2) / R(13); a
            2*13 + O(13^3)
            sage: R(13) / R(169 * 2)
            7*13^-1 + 6 + O(13)
            sage: ~a
            7*13^-1 + 6 + O(13)
            sage: 1 / a
            7*13^-1 + 6 + O(13)
        """
        K = self.parent().fraction_field()
        return K(self) / K(right)

    def _quo_rem(self, _right):
        """
        Quotient with remainder.

        EXAMPLES::

            sage: R = ZpCA(3, 5)
            sage: R(12).quo_rem(R(2))  # indirect doctest
            (2*3 + O(3^5), O(3^5))
            sage: R(2).quo_rem(R(12))  # indirect doctest
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
        """
        cdef CAElement right = _right
        if right._is_inexact_zero():
            raise ZeroDivisionError
        cdef CAElement q = self._new_c()
        cdef CAElement r = self._new_c()
        cdef long sval, srprec, rval, rrprec, diff, qrprec
        sval = self.valuation_c()
        srprec = self.absprec - sval
        rval = right.valuation_c()
        rrprec = right.absprec - rval
        diff = sval - rval
        rprec = min(srprec, rrprec)
        r.absprec = r.prime_pow.ram_prec_cap
        qrprec = diff + srprec
        if qrprec < 0:
            csetzero(q.value, q.prime_pow)
            q.absprec = 0
            r = self
        elif qrprec == 0:
            q._set_inexact_zero(0)
            ccopy(r.value, self.value, r.prime_pow)
        elif ciszero(self.value, self.prime_pow):
            q.absprec = diff + rprec
            csetzero(q.value, q.prime_pow)
            csetzero(r.value, r.prime_pow)
        elif diff >= 0:
            q.absprec = diff + rprec
            # shift right and self by the same power of the uniformizer
            cshift_notrunc(r.value, right.value, -rval, q.absprec, r.prime_pow, False)
            # We use shift_rem as a temp variable
            cshift_notrunc(self.prime_pow.shift_rem, self.value, -rval, q.absprec, q.prime_pow, False)
            # divide
            cdivunit(q.value, self.prime_pow.shift_rem, r.value, q.absprec, q.prime_pow)
            csetzero(r.value, r.prime_pow)
        else:
            q.absprec = min(qrprec, rrprec)
            cshift(q.value, r.value, self.value, -rval, q.absprec, q.prime_pow, False)
            cshift_notrunc(q.prime_pow.shift_rem, right.value, -rval, q.absprec, q.prime_pow, False)
            cdivunit(q.value, q.value, q.prime_pow.shift_rem, q.absprec, q.prime_pow)
        creduce(q.value, q.value, q.absprec, q.prime_pow)
        return q, r

    def __pow__(CAElement self, _right, dummy):
        r"""
        Exponentiation.

        When ``right`` is divisible by `p` then one can get more
        precision than expected.  See the documentation in
        :mod:`sage.rings.padics.CR_template.pxi` for more details.

        For `p`-adic exponents, `a^b` is defined as `\exp(b \log(a))`.
        Since the `p`-adic logarithm is defined for `a` a unit, the
        same is true of exponentiation.

        INPUT:

        - ``_right`` -- currently integers and `p`-adic exponents are supported

        - ``dummy`` -- not used (Python's ``__pow__`` signature includes it)

        EXAMPLES::

            sage: R = ZpCA(11, 5)
            sage: R(1/2)^5
            10 + 7*11 + 11^2 + 5*11^3 + 4*11^4 + O(11^5)
            sage: R(1/32)
            10 + 7*11 + 11^2 + 5*11^3 + 4*11^4 + O(11^5)
            sage: R(1/2)^5 == R(1/32)
            True
            sage: R(3)^1000
            1 + 4*11^2 + 3*11^3 + 7*11^4 + O(11^5)

        `p`-adic exponents are supported::

            sage: R = ZpCA(11, 5, print_mode='terse')
            sage: a = R(3/14, 3); b = R(8/9); c = R(11,2)
            sage: a
            1046 + O(11^3)
            sage: b
            35790 + O(11^5)
            sage: a^b
            177 + O(11^3)
            sage: a^35790
            177 + O(11^3)
            sage: a^c
            848 + O(11^3)
            sage: (a.log()*c).exp()
            848 + O(11^3)

            sage: R = ZpCA(19, 5, print_mode='series')
            sage: a = R(8/5,4); a
            13 + 7*19 + 11*19^2 + 7*19^3 + O(19^4)
            sage: a^(R(19/7))
            1 + 14*19^2 + 11*19^3 + 13*19^4 + O(19^5)
            sage: (a // R.teichmuller(13))^(R(19/7))
            1 + 14*19^2 + 11*19^3 + 13*19^4 + O(19^5)
            sage: (a.log() * 19/7).exp()
            1 + 14*19^2 + 11*19^3 + 13*19^4 + O(19^5)

        Check that :issue:`31875` is fixed::

            sage: R(1)^R(0)
            1 + O(19^5)

            sage: # needs sage.libs.flint
            sage: S.<a> = ZqCA(4)
            sage: S(1)^S(0)
            1 + O(2^20)
        """
        cdef long relprec, val, rval
        cdef mpz_t tmp
        cdef Integer right
        cdef CAElement pright, ans
        cdef bint exact_exp
        if isinstance(_right, (Integer, int, Rational)):
            if _right < 0:
                base = ~self
                return base.__pow__(-_right, dummy)
            exact_exp = True
        elif self.parent() is _right.parent():
            # For extension elements, we need to switch to the
            # fraction field sometimes in highly ramified extensions.
            exact_exp = (<CAElement>_right)._is_exact_zero()
            pright = _right
        else:
            self, _right = canonical_coercion(self, _right)
            return self.__pow__(_right, dummy)
        ans = self._new_c()
        if exact_exp and _right == 0:
            # return 1 to maximum precision
            ans.absprec = self.prime_pow.ram_prec_cap
            csetone(ans.value, ans.prime_pow)
        elif ciszero(self.value, self.prime_pow):
            # We may assume from above that right > 0 if exact.
            # So we return a zero of precision right * self.ordp.
            if isinstance(_right, int):
                _right = Integer(_right)
            if isinstance(_right, Integer):
                right = <Integer>_right
                if self.absprec == 0:
                    ans.absprec = 0
                else:
                    mpz_init(tmp)
                    mpz_mul_si(tmp, right.value, self.absprec)
                    if mpz_cmp_si(tmp, self.prime_pow.ram_prec_cap) >= 0:
                        ans.absprec = self.prime_pow.ram_prec_cap
                    else:
                        ans.absprec = mpz_get_si(tmp)
                    mpz_clear(tmp)
                csetzero(ans.value, ans.prime_pow)
            else:
                if not exact_exp and self.absprec > 0:
                    raise ValueError("in order to raise to a p-adic exponent, base must be a unit")
                raise PrecisionError("need more precision")
        else:
            val = self.valuation_c()
            if exact_exp:
                # exact_pow_helper is defined in padic_template_element.pxi
                right = exact_pow_helper(&relprec, self.absprec - val, _right, self.prime_pow)
                mpz_init(tmp)
                mpz_mul_si(tmp, right.value, val)
                if mpz_cmp_si(tmp, self.prime_pow.ram_prec_cap) >= 0:
                    ans.absprec = self.prime_pow.ram_prec_cap
                    csetzero(ans.value, ans.prime_pow)
                else:
                    ans.absprec = min(mpz_get_si(tmp) + relprec, self.prime_pow.ram_prec_cap)
                    cpow(ans.value, self.value, right.value, ans.absprec, ans.prime_pow)
                mpz_clear(tmp)
            else:
                rval = pright.valuation_c()
                if rval != 0:
                    pright = pright.unit_part()
                # We may assume that val = 0 since the following will quickly raise an error otherwise.
                # padic_pow_helper is defined in padic_template_element.pxi
                ans.absprec = padic_pow_helper(ans.value, self.value, val, self.absprec,
                                               pright.value, rval, pright.absprec, self.prime_pow)
        return ans

    cdef pAdicTemplateElement _lshift_c(self, long shift):
        r"""
        Multiply by `\pi^{\mbox{shift}}`.

        Negative shifts may truncate the result.

        TESTS::

            sage: R = ZpCA(5); a = R(17); a << 2
            2*5^2 + 3*5^3 + O(5^20)
            sage: a << -1
            3 + O(5^19)
            sage: a << 0 == a
            True
            sage: a << 400
            O(5^20)
            sage: a << -400
            O(5^0)
        """
        if shift < 0:
            return self._rshift_c(-shift)
        elif shift == 0:
            return self
        cdef CAElement ans = self._new_c()
        if shift >= self.prime_pow.ram_prec_cap:
            csetzero(ans.value, ans.prime_pow)
            ans.absprec = self.prime_pow.ram_prec_cap
        else:
            ans.absprec = min(self.absprec + shift, self.prime_pow.ram_prec_cap)
            cshift_notrunc(ans.value, self.value, shift, ans.absprec, ans.prime_pow, self.prime_pow.e > 1)
        return ans

    cdef pAdicTemplateElement _rshift_c(self, long shift):
        r"""
        Divide by ``π^{\mbox{shift}}``.

        Positive shifts may truncate the result.

        TESTS::

            sage: R = ZpCA(5); a = R(77); a >> 1
            3*5 + O(5^19)
            sage: a >> -1
            2*5 + 3*5^3 + O(5^20)
            sage: a >> 0 == a
            True
            sage: a >> 400
            O(5^0)
            sage: a >> -400
            O(5^20)
        """
        if shift < 0:
            return self._lshift_c(-shift)
        elif shift == 0:
            return self
        cdef CAElement ans = self._new_c()
        if shift >= self.absprec:
            csetzero(ans.value, ans.prime_pow)
            ans.absprec = 0
        else:
            ans.absprec = self.absprec - shift
            cshift(ans.value, ans.prime_pow.shift_rem, self.value, -shift, ans.absprec, ans.prime_pow, self.prime_pow.e > 1)
        return ans

    def add_bigoh(self, absprec):
        """
        Return a new element with absolute precision decreased to
        ``absprec``.  The precision never increases.

        INPUT:

        - ``absprec`` -- integer or infinity

        OUTPUT: ``self`` with precision set to the minimum of ``self``'s
        precision and ``prec``

        EXAMPLES::

            sage: R = Zp(7,4,'capped-abs','series'); a = R(8); a.add_bigoh(1)
            1 + O(7)

            sage: k = ZpCA(3,5)
            sage: a = k(41); a
            2 + 3 + 3^2 + 3^3 + O(3^5)
            sage: a.add_bigoh(7)
            2 + 3 + 3^2 + 3^3 + O(3^5)
            sage: a.add_bigoh(3)
            2 + 3 + 3^2 + O(3^3)

        TESTS:

        Verify that :issue:`13591` has been resolved::

            sage: k(3).add_bigoh(-1)
            O(3^-1)
        """
        cdef long aprec, newprec
        if absprec is infinity:
            return self
        if isinstance(absprec, int):
            aprec = absprec
        else:
            if not isinstance(absprec, Integer):
                absprec = Integer(absprec)
            if mpz_fits_slong_p((<Integer>absprec).value) == 0:
                if mpz_sgn((<Integer>absprec).value) == -1:
                    raise ValueError("absprec must fit into a signed long")
                else:
                    aprec = self.prime_pow.ram_prec_cap
            else:
                aprec = mpz_get_si((<Integer>absprec).value)
        if aprec >= self.absprec:
            return self
        if aprec < 0:
            return self.parent().fraction_field()(self).add_bigoh(absprec)
        cdef CAElement ans = self._new_c()
        ans.absprec = aprec
        creduce(ans.value, self.value, ans.absprec, ans.prime_pow)
        return ans

    cpdef bint _is_exact_zero(self) except -1:
        """
        Test whether this element is an exact zero, which is always
        ``False`` for capped absolute elements.

        This function exists for compatibility with capped relative
        elements.

        EXAMPLES::

            sage: ZpCA(5)(0)._is_exact_zero()
            False
        """
        return False

    cpdef bint _is_inexact_zero(self) except -1:
        """
        Determine whether this element is indistinguishable from
        zero.

        EXAMPLES::

            sage: R = ZpCA(7, 5)
            sage: R(7^5)._is_inexact_zero()
            True
            sage: R(0,4)._is_inexact_zero()
            True
            sage: R(0)._is_inexact_zero()
            True
        """
        return ciszero(self.value, self.prime_pow)

    def is_zero(self, absprec=None):
        r"""
        Determine whether this element is zero modulo
        `\pi^{\mbox{absprec}}`.

        If ``absprec is None``, returns ``True`` if this element is
        indistinguishable from zero.

        INPUT:

        - ``absprec`` -- integer, infinity, or ``None``

        EXAMPLES::

            sage: R = ZpCA(17, 6)
            sage: R(0).is_zero()
            True
            sage: R(17^6).is_zero()
            True
            sage: R(17^2).is_zero(absprec=2)
            True
            sage: R(17^6).is_zero(absprec=10)
            Traceback (most recent call last):
            ...
            PrecisionError: not enough precision to determine if element is zero
        """
        if absprec is infinity:
            raise PrecisionError("not enough precision to determine if element is zero")
        cdef bint iszero = ciszero(self.value, self.prime_pow)
        if absprec is None:
            return iszero
        cdef long val = self.valuation_c()
        if isinstance(absprec, int):
            if iszero and absprec > self.absprec:
                raise PrecisionError("not enough precision to determine if element is zero")
            return val >= absprec
        if not isinstance(absprec, Integer):
            absprec = Integer(absprec)
        if iszero:
            if mpz_cmp_si((<Integer>absprec).value, val) > 0:
                raise PrecisionError("not enough precision to determine if element is zero")
            else:
                return True
        return mpz_cmp_si((<Integer>absprec).value, val) <= 0

    def __bool__(self):
        """
        Whether this element should be considered true in a boolean context.

        For most applications, explicitly specifying the power of p
        modulo which the element is supposed to be nonzero is
        preferable.

        EXAMPLES::

            sage: R = ZpCA(5); a = R(0); b = R(0,5); c = R(75)
            sage: bool(a), bool(b), bool(c)
            (False, False, True)
        """
        return not ciszero(self.value, self.prime_pow)

    def is_equal_to(self, _right, absprec=None):
        r"""
        Determine whether the inputs are equal modulo
        `\pi^{\mbox{absprec}}`.

        INPUT:

        - ``right`` -- a `p`-adic element with the same parent

        - ``absprec`` -- integer, infinity, or ``None``

        EXAMPLES::

            sage: R = ZpCA(2, 6)
            sage: R(13).is_equal_to(R(13))
            True
            sage: R(13).is_equal_to(R(13+2^10))
            True
            sage: R(13).is_equal_to(R(17), 2)
            True
            sage: R(13).is_equal_to(R(17), 5)
            False
            sage: R(13).is_equal_to(R(13+2^10),absprec=10)
            Traceback (most recent call last):
            ...
            PrecisionError: elements not known to enough precision
        """
        if absprec is infinity:
            raise PrecisionError("elements not known to enough precision")
        cdef CAElement right
        cdef long aprec, rprec, sval, rval
        if self.parent() is _right.parent():
            right = _right
        else:
            right = self.parent()(_right)
        if absprec is None:
            aprec = min(self.absprec, right.absprec)
        else:
            if not isinstance(absprec, Integer):
                absprec = Integer(absprec)
            if mpz_fits_slong_p((<Integer>absprec).value) == 0:
                if mpz_sgn((<Integer>absprec).value) < 0:
                    return True
                else:
                    raise PrecisionError("elements not known to enough precision")
            aprec = mpz_get_si((<Integer>absprec).value)
            if aprec > self.absprec or aprec > right.absprec:
                raise PrecisionError("elements not known to enough precision")
        return ccmp(self.value, right.value, aprec, aprec < self.absprec, aprec < right.absprec, self.prime_pow) == 0

    cdef int _cmp_units(self, pAdicGenericElement _right) except -2:
        """
        This function is used in comparing `p`-adic elements.

        EXAMPLES::

            sage: R = ZpCA(37)
            sage: R(17) == R(17+37^6)  # indirect doctest
            False
        """
        cdef CAElement right = _right
        cdef long aprec = min(self.absprec, right.absprec)
        if aprec == 0:
            return 0
        return ccmp(self.value, right.value, aprec, aprec < self.absprec, aprec < right.absprec, self.prime_pow)

    cdef pAdicTemplateElement lift_to_precision_c(self, long absprec):
        """
        Return an arbitrary lift of this element to higher precision.

        If ``absprec`` is less than the absolute precision of this
        element this function will return the input element.

        INPUT:

        - ``absprec`` -- integer; at most the precision cap of the parent

        EXAMPLES::

            sage: R = ZpCA(19)
            sage: a = R(19, 7); a
            19 + O(19^7)
            sage: a.lift_to_precision(12)  # indirect doctest
            19 + O(19^12)
            sage: a.lift_to_precision(4) is a
            True
        """
        cdef CAElement ans
        if absprec == maxordp:
            absprec = self.prime_pow.ram_prec_cap
        if absprec <= self.absprec:
            return self
        ans = self._new_c()
        ccopy(ans.value, self.value, ans.prime_pow)
        ans.absprec = absprec
        return ans

    def _cache_key(self):
        r"""
        Return a hashable key which identifies this element for caching.

        TESTS::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqCA(9)
            sage: (9*a)._cache_key()
            (..., ((), (), (0, 1)), 20)

        .. SEEALSO::

            :meth:`sage.misc.cachefunc._cache_key`
        """
        def tuple_recursive(l):
            return tuple(tuple_recursive(x) for x in l) if isinstance(l, Iterable) else l

        return (self.parent(), tuple_recursive(trim_zeros(list(self.expansion()))), self.precision_absolute())

    def _teichmuller_set_unsafe(self):
        """
        Set this element to the Teichmuller representative with the
        same residue.

        .. WARNING::

            This function modifies the element, which is not safe.
            Elements are supposed to be immutable.

        EXAMPLES::

            sage: R = ZpCA(17,5); a = R(11)
            sage: a
            11 + O(17^5)
            sage: a._teichmuller_set_unsafe(); a
            11 + 14*17 + 2*17^2 + 12*17^3 + 15*17^4 + O(17^5)
            sage: E = a.expansion(lift_mode='teichmuller'); E
            17-adic expansion of 11 + 14*17 + 2*17^2 + 12*17^3 + 15*17^4 + O(17^5) (teichmuller)
            sage: list(E)
            [11 + 14*17 + 2*17^2 + 12*17^3 + 15*17^4 + O(17^5), O(17^5), O(17^5), O(17^5), O(17^5)]

        Note that if you set an element which is congruent to 0 you
        get 0 to maximum precision::

            sage: b = R(17*5); b
            5*17 + O(17^5)
            sage: b._teichmuller_set_unsafe(); b
            O(17^5)
        """
        if self.valuation_c() > 0:
            csetzero(self.value, self.prime_pow)
            self.absprec = self.prime_pow.ram_prec_cap
        elif self.absprec == 0:
            raise ValueError("not enough precision")
        else:
            cteichmuller(self.value, self.value, self.absprec, self.prime_pow)

    def _polynomial_list(self, pad=False):
        """
        Return the coefficient list for a polynomial over the base ring
        yielding this element.

        INPUT:

        - ``pad`` -- whether to pad the result with zeros of the appropriate precision

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<x> = ZZ[]
            sage: K.<a> = ZqCA(25)
            sage: W.<w> = K.extension(x^3 - 5)
            sage: (1 + w + O(w^11))._polynomial_list()
            [1 + O(5^4), 1 + O(5^4)]
            sage: (1 + w + O(w^11))._polynomial_list(pad=True)
            [1 + O(5^4), 1 + O(5^4), O(5^3)]
            sage: W(0)._polynomial_list()
            []
            sage: W(0)._polynomial_list(pad=True)
            [O(5^20), O(5^20), O(5^20)]
            sage: W(O(w^7))._polynomial_list()
            []
            sage: W(O(w^7))._polynomial_list(pad=True)
            [O(5^3), O(5^2), O(5^2)]
        """
        R = self.base_ring()
        prec = self.precision_absolute()
        e = self.parent().relative_e()
        L = ccoefficients(self.value, 0, self.absprec, self.prime_pow)
        if pad:
            n = self.parent().relative_degree()
            L.extend([R.zero()] * (n - len(L)))
        if e == 1:
            return [R(c, prec) for c in L]
        else:
            return [R(c, (prec - i - 1) // e + 1) for i, c in enumerate(L)]

    def polynomial(self, var='x'):
        """
        Return a polynomial over the base ring that yields this element
        when evaluated at the generator of the parent.

        INPUT:

        - ``var`` -- string; the variable name for the polynomial

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqCA(5^3)
            sage: a.polynomial()
            (1 + O(5^20))*x + O(5^20)
            sage: a.polynomial(var='y')
            (1 + O(5^20))*y + O(5^20)
            sage: (5*a^2 + R(25, 4)).polynomial()
            (5 + O(5^4))*x^2 + O(5^4)*x + 5^2 + O(5^4)
        """
        R = self.base_ring()
        S = R[var]
        return S(self._polynomial_list())

    def precision_absolute(self):
        """
        The absolute precision of this element.

        This is the power of the maximal ideal modulo which this
        element is defined.

        EXAMPLES::

            sage: R = Zp(7,4,'capped-abs'); a = R(7); a.precision_absolute()
            4
        """
        cdef Integer ans = Integer.__new__(Integer)
        mpz_set_si(ans.value, self.absprec)
        return ans

    def precision_relative(self):
        """
        The relative precision of this element.

        This is the power of the maximal ideal modulo which the unit
        part of this element is defined.

        EXAMPLES::

            sage: R = Zp(7,4,'capped-abs'); a = R(7); a.precision_relative()
            3
        """
        cdef Integer ans = Integer.__new__(Integer)
        mpz_set_si(ans.value, self.absprec - self.valuation_c())
        return ans

    cpdef pAdicTemplateElement unit_part(CAElement self):
        r"""
        Return the unit part of this element.

        EXAMPLES::

            sage: R = Zp(17,4,'capped-abs', 'val-unit')
            sage: a = R(18*17)
            sage: a.unit_part()
            18 + O(17^3)
            sage: type(a)
            <class 'sage.rings.padics.padic_capped_absolute_element.pAdicCappedAbsoluteElement'>
            sage: R(0).unit_part()
            O(17^0)
        """
        cdef CAElement ans = (<CAElement>self)._new_c()
        cdef long val = cremove(ans.value, (<CAElement>self).value, (<CAElement>self).absprec, (<CAElement>self).prime_pow, True)
        ans.absprec = (<CAElement>self).absprec - val
        return ans

    cdef long valuation_c(self) noexcept:
        """
        Return the valuation of this element.

        TESTS::

            sage: R = ZpCA(5)
            sage: R(5^5*1827).valuation()
            5
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
            sage: R(0).valuation()
            20
            sage: R(0,6).valuation()
            6
        """
        return cvaluation(self.value, self.absprec, self.prime_pow)

    cpdef val_unit(self):
        r"""
        Return a 2-tuple, the first element set to the valuation of this
        element, and the second to the unit part of this element.

        For a zero element, the unit part is ``O(p^0)``.

        EXAMPLES::

            sage: R = ZpCA(5)
            sage: a = R(75, 6); b = a - a
            sage: a.val_unit()
            (2, 3 + O(5^4))
            sage: b.val_unit()
            (6, O(5^0))
        """
        cdef CAElement unit = self._new_c()
        cdef Integer valuation = Integer.__new__(Integer)
        cdef long val = cremove(unit.value, self.value, self.absprec, self.prime_pow, True)
        mpz_set_si(valuation.value, val)
        unit.absprec = self.absprec - val
        return valuation, unit

    def __hash__(self):
        """
        Hashing.

        .. WARNING::

            Hashing of `p`-adic elements will likely be deprecated soon.  See :issue:`11895`.

        EXAMPLES::

            sage: R = ZpCA(11, 5)
            sage: hash(R(3)) == hash(3)
            True
        """
        return chash(self.value, 0, self.absprec, self.prime_pow)

cdef class pAdicCoercion_ZZ_CA(RingHomomorphism):
    """
    The canonical inclusion from the ring of integers to a capped absolute
    ring.

    EXAMPLES::

        sage: f = ZpCA(5).coerce_map_from(ZZ); f
        Ring morphism:
          From: Integer Ring
          To:   5-adic Ring with capped absolute precision 20

    TESTS::

        sage: TestSuite(f).run()
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = ZpCA(5).coerce_map_from(ZZ); type(f)
            <class 'sage.rings.padics.padic_capped_absolute_element.pAdicCoercion_ZZ_CA'>
        """
        RingHomomorphism.__init__(self, ZZ.Hom(R))
        self._zero = R.element_class(R, 0)
        self._section = pAdicConvert_CA_ZZ(R)

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpCA(5).coerce_map_from(ZZ)
            sage: g = copy(f)  # indirect doctest
            sage: g == f
            True
            sage: g(6)
            1 + 5 + O(5^20)
            sage: f(6) == g(6)
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

            sage: f = ZpCA(5).coerce_map_from(ZZ)
            sage: g = copy(f)  # indirect doctest
            sage: g == f
            True
            sage: g(6)
            1 + 5 + O(5^20)
            sage: f(6) == g(6)
            True
        """
        self._zero = _slots['_zero']
        self._section = _slots['_section']
        RingHomomorphism._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = ZpCA(5).coerce_map_from(ZZ)
            sage: f(0).parent()
            5-adic Ring with capped absolute precision 20
            sage: f(5)
            5 + O(5^20)
        """
        if mpz_sgn((<Integer>x).value) == 0:
            return self._zero
        cdef CAElement ans = self._zero._new_c()
        ans.absprec = ans.prime_pow.ram_prec_cap
        cconv_mpz_t(ans.value, (<Integer>x).value, ans.absprec, True, ans.prime_pow)
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in
        (relative or absolute or both).

        See the documentation for
        :meth:`pAdicCappedAbsoluteElement.__init__` for more details.

        EXAMPLES::

            sage: R = ZpCA(5,4)
            sage: type(R(10,2))
            <class 'sage.rings.padics.padic_capped_absolute_element.pAdicCappedAbsoluteElement'>
            sage: R(10,2)  # indirect doctest
            2*5 + O(5^2)
            sage: R(10,3,1)
            2*5 + O(5^2)
            sage: R(10,absprec=2)
            2*5 + O(5^2)
            sage: R(10,relprec=2)
            2*5 + O(5^3)
            sage: R(10,absprec=1)
            O(5)
            sage: R(10,empty=True)
            O(5^0)
        """
        cdef long val, aprec, rprec
        cdef CAElement ans
        _process_args_and_kwds(&aprec, &rprec, args, kwds, True, self._zero.prime_pow)
        if mpz_sgn((<Integer>x).value) == 0:
            if aprec >= self._zero.prime_pow.ram_prec_cap:
                return self._zero
            ans = self._zero._new_c()
            csetzero(ans.value, ans.prime_pow)
            ans.absprec = aprec
        else:
            val = get_ordp(x, self._zero.prime_pow)
            ans = self._zero._new_c()
            if aprec <= val:
                csetzero(ans.value, ans.prime_pow)
                ans.absprec = aprec
            else:
                ans.absprec = min(aprec, val + rprec)
                cconv_mpz_t(ans.value, (<Integer>x).value, ans.absprec, True, self._zero.prime_pow)
        return ans

    def section(self):
        """
        Return a map back to the ring of integers that approximates an element
        by an integer.

        EXAMPLES::

            sage: f = ZpCA(5).coerce_map_from(ZZ).section()
            sage: f(ZpCA(5)(-1)) - 5^20
            -1
        """
        from sage.misc.constant_function import ConstantFunction
        if not isinstance(self._section.domain, ConstantFunction):
            import copy
            self._section = copy.copy(self._section)
        return self._section

cdef class pAdicConvert_CA_ZZ(RingMap):
    """
    The map from a capped absolute ring back to the ring of integers that
    returns the smallest nonnegative integer approximation to its input
    which is accurate up to the precision.

    Raises a :exc:`ValueError` if the input is not in the closure of the image
    of the ring of integers.

    EXAMPLES::

        sage: f = ZpCA(5).coerce_map_from(ZZ).section(); f
        Set-theoretic ring morphism:
          From: 5-adic Ring with capped absolute precision 20
          To:   Integer Ring
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = ZpCA(5).coerce_map_from(ZZ).section(); type(f)
            <class 'sage.rings.padics.padic_capped_absolute_element.pAdicConvert_CA_ZZ'>
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

            sage: f = ZpCA(5).coerce_map_from(ZZ).section()
            sage: f(ZpCA(5)(-1)) - 5^20
            -1
            sage: f(ZpCA(5)(0))
            0
        """
        cdef Integer ans = Integer.__new__(Integer)
        cdef CAElement x = _x
        cconv_mpz_t_out(ans.value, x.value, 0, x.absprec, x.prime_pow)
        return ans

cdef class pAdicConvert_QQ_CA(Morphism):
    """
    The inclusion map from the rationals to a capped absolute ring that is
    defined on all elements with nonnegative `p`-adic valuation.

    EXAMPLES::

        sage: f = ZpCA(5).convert_map_from(QQ); f
        Generic morphism:
          From: Rational Field
          To:   5-adic Ring with capped absolute precision 20
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = ZpCA(5).convert_map_from(QQ); type(f)
            <class 'sage.rings.padics.padic_capped_absolute_element.pAdicConvert_QQ_CA'>
        """
        Morphism.__init__(self, Hom(QQ, R, SetsWithPartialMaps()))
        self._zero = R.element_class(R, 0)

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpCA(5).convert_map_from(QQ)
            sage: g = copy(f)  # indirect doctest
            sage: g == f # todo: comparison not implemented
            True
            sage: g(1/6)
            1 + 4*5 + 4*5^3 + 4*5^5 + 4*5^7 + 4*5^9 + 4*5^11 + 4*5^13 + 4*5^15 + 4*5^17 + 4*5^19 + O(5^20)
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

            sage: f = ZpCA(5).convert_map_from(QQ)
            sage: g = copy(f)  # indirect doctest
            sage: g == f # todo: comparison not implemented
            True
            sage: g(1/6)
            1 + 4*5 + 4*5^3 + 4*5^5 + 4*5^7 + 4*5^9 + 4*5^11 + 4*5^13 + 4*5^15 + 4*5^17 + 4*5^19 + O(5^20)
            sage: g(1/6) == f(1/6)
            True
        """
        self._zero = _slots['_zero']
        Morphism._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = ZpCA(5,4).convert_map_from(QQ)
            sage: f(1/7)
            3 + 3*5 + 2*5^3 + O(5^4)
            sage: f(0)
            O(5^4)
        """
        if mpq_sgn((<Rational>x).value) == 0:
            return self._zero
        cdef CAElement ans = self._zero._new_c()
        cconv_mpq_t(ans.value, (<Rational>x).value, ans.prime_pow.ram_prec_cap, True, ans.prime_pow)
        ans.absprec = ans.prime_pow.ram_prec_cap
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in (relative or absolute or both).

        See the documentation for :meth:`pAdicCappedAbsoluteElement.__init__` for more details.

        EXAMPLES::

            sage: R = ZpCA(5,4)
            sage: type(R(10/3,2))
            <class 'sage.rings.padics.padic_capped_absolute_element.pAdicCappedAbsoluteElement'>
            sage: R(10/3,2)  # indirect doctest
            4*5 + O(5^2)
            sage: R(10/3,3,1)
            4*5 + O(5^2)
            sage: R(10/3,absprec=2)
            4*5 + O(5^2)
            sage: R(10/3,relprec=2)
            4*5 + 5^2 + O(5^3)
            sage: R(10/3,absprec=1)
            O(5)
            sage: R(10/3,empty=True)
            O(5^0)
            sage: R(3/100,relprec=3)
            Traceback (most recent call last):
            ...
            ValueError: p divides denominator
        """
        cdef long val, aprec, rprec
        cdef CAElement ans
        _process_args_and_kwds(&aprec, &rprec, args, kwds, True, self._zero.prime_pow)
        if mpq_sgn((<Rational>x).value) == 0:
            if aprec >= self._zero.prime_pow.ram_prec_cap:
                return self._zero
            ans = self._zero._new_c()
            csetzero(ans.value, ans.prime_pow)
            ans.absprec = aprec
        else:
            val = get_ordp(x, self._zero.prime_pow)
            ans = self._zero._new_c()
            if aprec <= val:
                csetzero(ans.value, ans.prime_pow)
                ans.absprec = aprec
            else:
                ans.absprec = min(aprec, val + rprec)
                cconv_mpq_t(ans.value, (<Rational>x).value, ans.absprec, True, self._zero.prime_pow)
        return ans

cdef class pAdicCoercion_CA_frac_field(RingHomomorphism):
    """
    The canonical inclusion of Zq into its fraction field.

    EXAMPLES::

        sage: # needs sage.libs.flint
        sage: R.<a> = ZqCA(27, implementation='FLINT')
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
            sage: R.<a> = ZqCA(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R); type(f)
            <class 'sage.rings.padics.qadic_flint_CA.pAdicCoercion_CA_frac_field'>
        """
        RingHomomorphism.__init__(self, R.Hom(K))
        self._zero = K(0)
        self._section = pAdicConvert_CA_frac_field(K, R)

    cpdef Element _call_(self, _x):
        """
        Evaluation.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqCA(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f(a)
            a + O(3^20)
        """
        cdef CAElement x = _x
        cdef CRElement ans = self._zero._new_c()
        ans.ordp = 0
        ans.relprec = x.absprec
        ccopy(ans.unit, x.value, x.prime_pow)
        IF CELEMENT_IS_PY_OBJECT:
            # The base ring is wrong, so we fix it.
            K = ans.unit.base_ring()
            ans.unit._coeffs = [K(c) for c in ans.unit._coeffs]
        ans._normalize()
        return ans

    cpdef Element _call_with_args(self, _x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in
        (relative or absolute or both).

        See the documentation for
        :meth:`pAdicCappedAbsoluteElement.__init__` for more details.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqCA(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f(a, 3)  # indirect doctest
            a + O(3^3)
            sage: b = 9*a
            sage: f(b, 3)
            a*3^2 + O(3^3)
            sage: f(b, 4, 1)
            a*3^2 + O(3^3)
            sage: f(b, 4, 3)
            a*3^2 + O(3^4)
            sage: f(b, absprec=4)
            a*3^2 + O(3^4)
            sage: f(b, relprec=3)
            a*3^2 + O(3^5)
            sage: f(b, absprec=1)
            O(3)
            sage: f(R(0))
            O(3^20)
        """
        cdef long aprec, rprec
        cdef CAElement x = _x
        cdef CRElement ans = self._zero._new_c()
        _process_args_and_kwds(&aprec, &rprec, args, kwds, False, ans.prime_pow)
        if x.absprec < aprec:
            aprec = x.absprec
        ans.ordp = cremove(ans.unit, x.value, aprec, x.prime_pow, True)
        ans.relprec = aprec - ans.ordp
        if rprec < ans.relprec:
            ans.relprec = rprec
        if ans.relprec < 0:
            ans.relprec = 0
            ans.ordp = aprec
            csetzero(ans.unit, x.prime_pow)
        else:
            IF CELEMENT_IS_PY_OBJECT:
                # The base ring is wrong, so we fix it.
                K = ans.unit.base_ring()
                ans.unit._coeffs = [K(c) for c in ans.unit._coeffs]
            pass
        return ans

    def section(self):
        """
        Return a map back to the ring that converts elements of
        nonnegative valuation.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqCA(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f(K.gen())
            a + O(3^20)
            sage: f.section()
            Generic morphism:
              From: 3-adic Unramified Extension Field in a defined by x^3 + 2*x + 1
              To:   3-adic Unramified Extension Ring in a defined by x^3 + 2*x + 1
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
            sage: R.<a> = ZqCA(27, implementation='FLINT')
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
            a + O(3^20)
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
            sage: R.<a> = ZqCA(9, implementation='FLINT')
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
            a + O(3^20)
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
            sage: R.<a> = ZqCA(9, implementation='FLINT')
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
            sage: R.<a> = ZqCA(9, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f.is_surjective()
            False
        """
        return False


cdef class pAdicConvert_CA_frac_field(Morphism):
    r"""
    The section of the inclusion from `\ZZ_q` to its fraction field.

    EXAMPLES::

        sage: # needs sage.libs.flint
        sage: R.<a> = ZqCA(27, implementation='FLINT')
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
            sage: R.<a> = ZqCA(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K); type(f)
            <class 'sage.rings.padics.qadic_flint_CA.pAdicConvert_CA_frac_field'>
        """
        Morphism.__init__(self, Hom(K, R, SetsWithPartialMaps()))
        self._zero = R(0)

    cpdef Element _call_(self, _x):
        """
        Evaluation.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqCA(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K)
            sage: f(K.gen())  # indirect doctest
            a + O(3^20)
        """
        cdef CRElement x = _x
        if x.ordp < 0:
            raise ValueError("negative valuation")
        cdef CAElement ans = self._zero._new_c()
        cdef bint reduce = (x.prime_pow.e > 1)
        ans.absprec = x.relprec + x.ordp
        if ans.absprec > ans.prime_pow.ram_prec_cap:
            ans.absprec = ans.prime_pow.ram_prec_cap
            reduce = True
        if x.ordp >= ans.absprec:
            csetzero(ans.value, ans.prime_pow)
        else:
            cshift_notrunc(ans.value, x.unit, x.ordp, ans.absprec, ans.prime_pow, reduce)
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
            sage: R.<a> = ZqCA(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K); a = K(a)
            sage: f(a, 3)  # indirect doctest
            a + O(3^3)
            sage: b = 9*a
            sage: f(b, 3)
            a*3^2 + O(3^3)
            sage: f(b, 4, 1)
            a*3^2 + O(3^3)
            sage: f(b, 4, 3)
            a*3^2 + O(3^4)
            sage: f(b, absprec=4)
            a*3^2 + O(3^4)
            sage: f(b, relprec=3)
            a*3^2 + O(3^5)
            sage: f(b, absprec=1)
            O(3)
            sage: f(K(0))
            O(3^20)
        """
        cdef long aprec, rprec
        cdef CRElement x = _x
        if x.ordp < 0:
            raise ValueError("negative valuation")
        cdef CAElement ans = self._zero._new_c()
        cdef bint reduce = False
        _process_args_and_kwds(&aprec, &rprec, args, kwds, True, ans.prime_pow)
        if x.relprec < rprec:
            rprec = x.relprec
            reduce = True
        ans.absprec = rprec + x.ordp
        if aprec < ans.absprec:
            ans.absprec = aprec
            reduce = True
        if x.ordp >= ans.absprec:
            csetzero(ans.value, ans.prime_pow)
        else:
            cshift_notrunc(ans.value, x.unit, x.ordp, ans.absprec, ans.prime_pow, reduce)
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
            sage: R.<a> = ZqCA(27, implementation='FLINT')
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
            a + O(3^20)
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
            sage: R.<a> = ZqCA(9, implementation='FLINT')
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
            a + O(3^20)
            sage: g(a) == f(a)
            True
        """
        self._zero = _slots['_zero']
        Morphism._update_slots(self, _slots)


def unpickle_cae_v2(cls, parent, value, absprec):
    r"""
    Unpickle capped absolute elements.

    INPUT:

    - ``cls`` -- the class of the capped absolute element

    - ``parent`` -- a `p`-adic ring

    - ``value`` -- a Python object wrapping a celement, of the kind
      accepted by the cunpickle function

    - ``absprec`` -- a Python int or Sage integer

    EXAMPLES::

        sage: from sage.rings.padics.padic_capped_absolute_element import unpickle_cae_v2, pAdicCappedAbsoluteElement
        sage: R = ZpCA(5,8)
        sage: a = unpickle_cae_v2(pAdicCappedAbsoluteElement, R, 42, int(6)); a
        2 + 3*5 + 5^2 + O(5^6)
        sage: a.parent() is R
        True
    """
    cdef CAElement ans = cls.__new__(cls)
    ans._parent = parent
    ans.prime_pow = <PowComputer_?>parent.prime_pow
    IF CELEMENT_IS_PY_OBJECT:
        polyt = type(ans.prime_pow.modulus)
        ans.value = <celement>polyt.__new__(polyt)
    cconstruct(ans.value, ans.prime_pow)
    cunpickle(ans.value, value, ans.prime_pow)
    ans.absprec = absprec
    return ans
