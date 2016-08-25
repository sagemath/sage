"""
Capped absolute template for complete discrete valuation rings

In order to use this template you need to write a linkage file and gluing file.
For an example see mpz_linkage.pxi (linkage file) and padic_capped_absolute_element.pyx (gluing file).

The linkage file implements a common API that is then used in the class CAElement defined here.
See the documentation of mpz_linkage.pxi for the functions needed.

The gluing file does the following:

- ctypedef's celement to be the appropriate type (e.g. mpz_t)
- includes the linkage file
- includes this template
- defines a concrete class inheriting from CAElement, and implements
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
        Sets the value of this element from given defining data.

        This function is intended for use in conversion, and should
        not be called on an element created with :meth:`_new_c`.

        INPUT:

        - ``x`` -- data defining a `p`-adic element: int, long,
          Integer, Rational, other `p`-adic element...

        - ``val`` -- the valuation of the resulting element

        - ``xprec -- an inherent precision of ``x``

        - ``absprec`` -- an absolute precision cap for this element

        - ``relprec`` -- a relative precision cap for this element

        TESTS::

            sage: R = ZpCA(5)
            sage: a = R(17,5); a #indirect doctest
            2 + 3*5 + O(5^5)
            sage: a = R(75, absprec = 5, relprec = 4); a #indirect doctest
            3*5^2 + O(5^5)
            sage: a = R(25/9, absprec = 5); a #indirect doctest
            4*5^2 + 2*5^3 + O(5^5)
            sage: a = R(25/9, absprec = 5, relprec = 4); a #indirect doctest
            4*5^2 + 2*5^3 + O(5^5)
        """
        cconstruct(self.value, self.prime_pow)
        cdef long rprec = comb_prec(relprec, self.prime_pow.prec_cap)
        cdef long aprec = comb_prec(absprec, min(self.prime_pow.prec_cap, xprec))
        if aprec <= val:
            csetzero(self.value, self.prime_pow)
            self.absprec = aprec
        else:
            self.absprec = min(aprec, val + rprec)
            if isinstance(x,CAElement) and x.parent() is self.parent():
                cshift(self.value, (<CAElement>x).value, 0, self.absprec, self.prime_pow, True)
            else:
                cconv(self.value, x, self.absprec, 0, self.prime_pow)

    cdef CAElement _new_c(self):
        """
        Creates a new element with the same basic info.

        TESTS::

            sage: R = ZpCA(5); R(6,5) * R(7,8) #indirect doctest
            2 + 3*5 + 5^2 + O(5^5)
        """
        cdef type t = type(self)
        cdef CAElement ans = t.__new__(t)
        ans._parent = self._parent
        ans.prime_pow = self.prime_pow
        cconstruct(ans.value, ans.prime_pow)
        return ans

    cdef int check_preccap(self) except -1:
        """
        Checks that this element doesn't have precision higher than
        allowed by the precision cap.

        TESTS::

            sage: ZpCA(5)(1).lift_to_precision(30) # indirect doctest
            Traceback (most recent call last):
            ...
            PrecisionError: Precision higher than allowed by the precision cap.
        """
        if self.absprec > self.prime_pow.prec_cap:
            raise PrecisionError("Precision higher than allowed by the precision cap.")

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
            <type 'sage.rings.padics.padic_capped_absolute_element.pAdicCappedAbsoluteElement'>
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
            sage: -a #indirect doctest
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
            sage: R(2) + R(3) #indirect doctest
            5 + O(13^4)
            sage: R(12) + R(1)
            13 + O(13^4)

        Check that :trac:`20245` is resolved::

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
            sage: R(10) - R(10) #indirect doctest
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
            7 + 3*17 + 10*17^2 + 13*17^3 + 6*17^4 + 3*17^5 + 10*17^6 + 13*17^7 + 6*17^8 + 3*17^9 + 10*17^10 + 13*17^11 + 6*17^12 + 3*17^13 + 10*17^14 + 13*17^15 + 6*17^16 + 3*17^17 + 10*17^18 + 13*17^19 + O(17^20)
            sage: ~R(-1) == R(-1) #indirect doctest
            True
        """
        return ~self.parent().fraction_field()(self)

    cpdef _mul_(self, _right):
        """
        Return the product of this element and ``_right``.

        EXAMPLES::

            sage: R = ZpCA(5)
            sage: a = R(20,5); b = R(75, 4); a * b #indirect doctest
            2*5^3 + 2*5^4 + O(5^5)
        """
        cdef CAElement right = _right
        cdef CAElement ans = self._new_c()
        cdef long vals, valr
        if self.absprec == self.prime_pow.prec_cap and right.absprec == self.prime_pow.prec_cap:
            ans.absprec = self.absprec
        else:
            vals = self.valuation_c()
            valr = right.valuation_c()
            ans.absprec = min(vals + valr + min(self.absprec - vals, right.absprec - valr), self.prime_pow.prec_cap)
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
            sage: R(2) / R(3) # indirect doctest
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

    def __pow__(CAElement self, _right, dummy):
        """
        Exponentiation.

        When ``right`` is divisible by `p` then one can get more
        precision than expected.  See the documentation in
        :mod:`sage.rings.padics.CR_template.pxi` for more details.

        For `p`-adic exponents, `a^b` is defined as `\exp(b \log(a))`.
        Since the `p`-adic logarithm is defined for `a` a unit, the
        same is true of exponentiation.

        INPUT:

        - ``_right`` -- currently integers and `p`-adic exponents are
          supported.

        - ``dummy`` -- not used (Python's ``__pow__`` signature
          includes it)

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
        """
        cdef long relprec, val, rval
        cdef mpz_t tmp
        cdef Integer right
        cdef CAElement pright, ans
        cdef bint exact_exp
        if isinstance(_right, Integer) or isinstance(_right, (int, long)) \
                                          or isinstance(_right, Rational):
            if _right < 0:
                base = ~self
                return base.__pow__(-_right, dummy)
            exact_exp = True
        elif self.parent() is _right.parent():
            ## For extension elements, we need to switch to the
            ## fraction field sometimes in highly ramified extensions.
            exact_exp = False
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
            if isinstance(_right, (int, long)):
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
                raise PrecisionError("Need more precision")
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
        """
        Multiplies by `\pi^{\mbox{shift}}`.

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
        if shift >= self.prime_pow.prec_cap:
            csetzero(ans.value, ans.prime_pow)
            ans.absprec = self.prime_pow.prec_cap
        else:
            ans.absprec = min(self.absprec + shift, self.prime_pow.prec_cap)
            cshift(ans.value, self.value, shift, ans.absprec, ans.prime_pow, False)
        return ans

    cdef pAdicTemplateElement _rshift_c(self, long shift):
        """
        Divides by ``\pi^{\mbox{shift}}``.

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
            cshift(ans.value, self.value, -shift, ans.absprec, ans.prime_pow, False)
        return ans

    def add_bigoh(self, absprec):
        """
        Returns a new element with absolute precision decreased to
        ``absprec``.  The precision never increases.

        INPUT:

        - ``absprec`` -- an integer

        OUTPUT:

        ``self`` with precision set to the minimum of ``self's`` precision and ``prec``

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
        """
        cdef long aprec, newprec
        if isinstance(absprec, int):
            aprec = absprec
        else:
            if not isinstance(absprec, Integer):
                absprec = Integer(absprec)
            aprec = mpz_get_si((<Integer>absprec).value)
        if aprec >= self.absprec:
            return self
        cdef CAElement ans = self._new_c()
        ans.absprec = aprec
        creduce(ans.value, self.value, ans.absprec, ans.prime_pow)
        return ans

    cpdef bint _is_exact_zero(self) except -1:
        """
        Tests whether this element is an exact zero, which is always
        False for capped absolute elements.

        This function exists for compatibility with capped relative
        elements.

        EXAMPLES::

            sage: ZpCA(5)(0)._is_exact_zero()
            False
        """
        return False

    cpdef bint _is_inexact_zero(self) except -1:
        """
        Determines whether this element is indistinguishable from
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

    def is_zero(self, absprec = None):
        r"""
        Determines whether this element is zero modulo
        `\pi^{\mbox{absprec}}`.

        If ``absprec is None``, returns ``True`` if this element is
        indistinguishable from zero.

        INPUT:

        - ``absprec`` -- an integer, infinity, or ``None``

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
            PrecisionError: Not enough precision to determine if element is zero
        """
        if absprec is infinity:
            raise PrecisionError("Not enough precision to determine if element is zero")
        cdef bint iszero = ciszero(self.value, self.prime_pow)
        if absprec is None:
            return iszero
        cdef long val = self.valuation_c()
        if isinstance(absprec, int):
            if iszero and absprec > self.absprec:
                raise PrecisionError("Not enough precision to determine if element is zero")
            return val >= absprec
        if not isinstance(absprec, Integer):
            absprec = Integer(absprec)
        if iszero:
            if mpz_cmp_si((<Integer>absprec).value, val) > 0:
                raise PrecisionError("Not enough precision to determine if element is zero")
            else:
                return True
        return mpz_cmp_si((<Integer>absprec).value, val) <= 0

    def __nonzero__(self):
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
        Determines whether the inputs are equal modulo
        `\pi^{\mbox{absprec}}`.

        INPUT:

        - ``right`` -- a `p`-adic element with the same parent

        - ``absprec`` -- an integer, infinity, or ``None``

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
            PrecisionError: Elements not known to enough precision
        """
        if absprec is infinity:
            raise PrecisionError("Elements not known to enough precision")
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
                    raise PrecisionError("Elements not known to enough precision")
            aprec = mpz_get_si((<Integer>absprec).value)
            if aprec > self.absprec or aprec > right.absprec:
                raise PrecisionError("Elements not known to enough precision")
        return ccmp(self.value, right.value, aprec, aprec < self.absprec, aprec < right.absprec, self.prime_pow) == 0

    cdef int _cmp_units(self, pAdicGenericElement _right) except -2:
        """
        This function is used in comparing `p`-adic elements.

        EXAMPLES::

            sage: R = ZpCA(37)
            sage: R(17) == R(17+37^6) # indirect doctest
            False
        """
        cdef CAElement right = _right
        cdef long aprec = min(self.absprec, right.absprec)
        if aprec == 0:
            return 0
        return ccmp(self.value, right.value, aprec, aprec < self.absprec, aprec < right.absprec, self.prime_pow)

    cdef pAdicTemplateElement lift_to_precision_c(self, long absprec):
        """
        Returns an arbitrary lift of this element to higher precision.

        If ``absprec`` is less than the absolute precision of this
        element this function will return the input element.

        INPUT:

        - ``absprec`` -- an integer, at most the precision cap of the
          parent.

        EXAMPLES::

            sage: R = ZpCA(19)
            sage: a = R(19, 7); a
            19 + O(19^7)
            sage: a.lift_to_precision(12) # indirect doctest
            19 + O(19^12)
            sage: a.lift_to_precision(4) is a
            True
        """
        cdef CAElement ans
        if absprec == maxordp:
            absprec = self.prime_pow.prec_cap
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

            sage: R.<a> = ZqCA(9)
            sage: (9*a)._cache_key()
            (..., ((), (), (0, 1)), 20)

        .. SEEALSO::

            :meth:`sage.misc.cachefunc._cache_key`
        """
        tuple_recursive = lambda l: tuple(tuple_recursive(x) for x in l) if isinstance(l, list) else l
        return (self.parent(), tuple_recursive(self.list()), self.precision_absolute())

    def list(self, lift_mode = 'simple', start_val = None):
        """
        Returns a list of coefficients of `p` starting with `p^0`.

        For each lift mode, this function returns a list of `a_i` so
        that this element can be expressed as

        .. MATH::

            \pi^v \cdot \sum_{i=0}^\infty a_i \pi^i

        where `v` is the valuation of this element when the parent is
        a field, and `v = 0` otherwise.

        Different lift modes affect the choice of `a_i`.  When
        ``lift_mode`` is ``'simple'``, the resulting `a_i` will be
        non-negative: if the residue field is `\mathbb{F}_p` then they
        will be integers with `0 \le a_i < p`; otherwise they will be
        a list of integers in the same range giving the coefficients
        of a polynomial in the indeterminant representing the maximal
        unramified subextension.

        Choosing ``lift_mode`` as ``'smallest'`` is similar to
        ``'simple'``, but uses a balanced representation `-p/2 < a_i
        \le p/2`.

        Finally, setting ``lift_mode = 'teichmuller'`` will yield
        Teichmuller representatives for the `a_i`: `a_i^q = a_i`.  In
        this case the `a_i` will also be `p`-adic elements.

        INPUT:

        - ``lift_mode`` -- ``'simple'``, ``'smallest'`` or
          ``'teichmuller'`` (default ``'simple'``)

        - ``start_val`` -- start at this valuation rather than the
          default (`0` or the valuation of this element).  If
          ``start_val`` is larger than the valuation of this element
          a ``ValueError`` is raised.

        .. NOTE::

            Use slice operators to get a particular range.

        EXAMPLES::

            sage: R = ZpCA(7,6); a = R(12837162817); a
            3 + 4*7 + 4*7^2 + 4*7^4 + O(7^6)
            sage: L = a.list(); L
            [3, 4, 4, 0, 4]
            sage: sum([L[i] * 7^i for i in range(len(L))]) == a
            True
            sage: L = a.list('smallest'); L
            [3, -3, -2, 1, -3, 1]
            sage: sum([L[i] * 7^i for i in range(len(L))]) == a
            True
            sage: L = a.list('teichmuller'); L
            [3 + 4*7 + 6*7^2 + 3*7^3 + 2*7^5 + O(7^6),
            O(7^5),
            5 + 2*7 + 3*7^3 + O(7^4),
            1 + O(7^3),
            3 + 4*7 + O(7^2),
            5 + O(7)]
            sage: sum([L[i] * 7^i for i in range(len(L))])
            3 + 4*7 + 4*7^2 + 4*7^4 + O(7^6)

        If the element has positive valuation then the list will start
        with some zeros::

            sage: a = R(7^3 * 17)
            sage: a.list()
            [0, 0, 0, 3, 2]
        """
        if ciszero(self.value, self.prime_pow):
            return []
        if lift_mode == 'teichmuller':
            vlist = self.teichmuller_list()
        elif lift_mode == 'simple':
            vlist = clist(self.value, self.absprec, True, self.prime_pow)
        elif lift_mode == 'smallest':
            vlist = clist(self.value, self.absprec, False, self.prime_pow)
        else:
            raise ValueError("unknown lift_mode")
        if start_val is not None:
            if start_val > 0:
                if start_val > self.valuation_c():
                    raise ValueError("starting valuation must be smaller than the element's valuation.  See slice()")
                vlist = vlist[start_val:]
            elif start_val < 0:
                if lift_mode == 'teichmuller':
                    zero = self.parent()(0)
                else:
                    # needs to be defined in the linkage file.
                    zero = _list_zero
                vlist = [zero] * (-start_val) + vlist
        return vlist

    def teichmuller_list(self):
        r"""
        Returns a list `[a_0, a_1,\ldots, a_n]` such that

        - `a_i^q = a_i`, where `q` is the cardinality of the residue field,

        - ``self`` equals `\sum_{i = 0}^n a_i \pi^i`, and

        - if `a_i \ne 0`, the absolute precision of `a_i` is
          ``self.precision_relative() - i``

        EXAMPLES::

            sage: R = ZpCA(5,5); R(14).list('teichmuller') #indirect doctest
            [4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5),
            3 + 3*5 + 2*5^2 + 3*5^3 + O(5^4),
            2 + 5 + 2*5^2 + O(5^3),
            1 + O(5^2),
            4 + O(5)]
        """
        ans = PyList_New(0)
        if ciszero(self.value, self.prime_pow):
            return ans
        cdef long curpower = self.absprec
        cdef CAElement list_elt
        cdef CAElement tmp = self._new_c()
        ccopy(tmp.value, self.value, self.prime_pow)
        while not ciszero(tmp.value, tmp.prime_pow) and curpower > 0:
            list_elt = self._new_c()
            cteichmuller(list_elt.value, tmp.value, curpower, self.prime_pow)
            if ciszero(list_elt.value, self.prime_pow):
                cshift_notrunc(tmp.value, tmp.value, -1, curpower-1, self.prime_pow)
            else:
                csub(tmp.value, tmp.value, list_elt.value, curpower, self.prime_pow)
                cshift_notrunc(tmp.value, tmp.value, -1, curpower-1, self.prime_pow)
                creduce(tmp.value, tmp.value, curpower-1, self.prime_pow)
            list_elt.absprec = curpower
            curpower -= 1
            PyList_Append(ans, list_elt)
        return ans

    def _teichmuller_set_unsafe(self):
        """
        Sets this element to the Teichmuller representative with the
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
            sage: a.list('teichmuller')
            [11 + 14*17 + 2*17^2 + 12*17^3 + 15*17^4 + O(17^5)]

        Note that if you set an element which is congruent to 0 you
        get 0 to maximum precision::

            sage: b = R(17*5); b
            5*17 + O(17^5)
            sage: b._teichmuller_set_unsafe(); b
            O(17^5)
        """
        if self.valuation_c() > 0:
            csetzero(self.value, self.prime_pow)
            self.absprec = self.prime_pow.prec_cap
        elif self.absprec == 0:
            raise ValueError("not enough precision")
        else:
            cteichmuller(self.value, self.value, self.absprec, self.prime_pow)

    def precision_absolute(self):
        """
        The absolute precision of this element.

        This is the power of the maximal ideal modulo which this
        element is defined.

        EXAMPLES::

            sage: R = Zp(7,4,'capped-abs'); a = R(7); a.precision_absolute()
            4
        """
        cdef Integer ans = PY_NEW(Integer)
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
        cdef Integer ans = PY_NEW(Integer)
        mpz_set_si(ans.value, self.absprec - self.valuation_c())
        return ans

    cpdef pAdicTemplateElement unit_part(CAElement self):
        r"""
        Returns the unit part of this element.

        EXAMPLES::

            sage: R = Zp(17,4,'capped-abs', 'val-unit')
            sage: a = R(18*17)
            sage: a.unit_part()
            18 + O(17^3)
            sage: type(a)
            <type 'sage.rings.padics.padic_capped_absolute_element.pAdicCappedAbsoluteElement'>
            sage: R(0).unit_part()
            O(17^0)
        """
        cdef CAElement ans = (<CAElement>self)._new_c()
        cdef long val = cremove(ans.value, (<CAElement>self).value, (<CAElement>self).absprec, (<CAElement>self).prime_pow)
        ans.absprec = (<CAElement>self).absprec - val
        return ans

    cdef long valuation_c(self):
        """
        Returns the valuation of this element.

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
        """
        Returns a 2-tuple, the first element set to the valuation of this
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
        cdef Integer valuation = PY_NEW(Integer)
        cdef long val = cremove(unit.value, self.value, self.absprec, self.prime_pow)
        mpz_set_si(valuation.value, val)
        unit.absprec = self.absprec - val
        return valuation, unit

    def __hash__(self):
        """
        Hashing.

        .. WARNING::

            Hashing of `p`-adic elements will likely be deprecated soon.  See :trac:`11895`.

        EXAMPLES::

            sage: R = ZpCA(11, 5)
            sage: hash(R(3)) == hash(3)
            True
        """
        return chash(self.value, 0, self.absprec, self.prime_pow)

cdef class pAdicCoercion_ZZ_CA(RingHomomorphism_coercion):
    """
    The canonical inclusion from the ring of integers to a capped absolute
    ring.

    EXAMPLES::

        sage: f = ZpCA(5).coerce_map_from(ZZ); f
        Ring Coercion morphism:
          From: Integer Ring
          To:   5-adic Ring with capped absolute precision 20
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = ZpCA(5).coerce_map_from(ZZ); type(f)
            <type 'sage.rings.padics.padic_capped_absolute_element.pAdicCoercion_ZZ_CA'>
        """
        RingHomomorphism_coercion.__init__(self, ZZ.Hom(R), check=False)
        self._zero = R._element_constructor(R, 0)
        self._section = pAdicConvert_CA_ZZ(R)

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpCA(5).coerce_map_from(ZZ)
            sage: g = copy(f) # indirect doctest
            sage: g == f
            True
            sage: g(6)
            1 + 5 + O(5^20)
            sage: f(6) == g(6)
            True
        """
        _slots['_zero'] = self._zero
        _slots['_section'] = self._section
        return RingHomomorphism_coercion._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpCA(5).coerce_map_from(ZZ)
            sage: g = copy(f) # indirect doctest
            sage: g == f
            True
            sage: g(6)
            1 + 5 + O(5^20)
            sage: f(6) == g(6)
            True
        """
        self._zero = _slots['_zero']
        self._section = _slots['_section']
        RingHomomorphism_coercion._update_slots(self, _slots)

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
        ans.absprec = ans.prime_pow.prec_cap
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
            <type 'sage.rings.padics.padic_capped_absolute_element.pAdicCappedAbsoluteElement'>
            sage: R(10,2)
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
            if aprec >= self._zero.prime_pow.prec_cap:
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
        Returns a map back to the ring of integers that approximates an element
        by an integer.

        EXAMPLES::

            sage: f = ZpCA(5).coerce_map_from(ZZ).section()
            sage: f(ZpCA(5)(-1)) - 5^20
            -1
        """
        return self._section

cdef class pAdicConvert_CA_ZZ(RingMap):
    """
    The map from a capped absolute ring back to the ring of integers that
    returns the the smallest non-negative integer approximation to its input
    which is accurate up to the precision.

    Raises a ``ValueError`` if the input is not in the closure of the image of
    the ring of integers.

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
            <type 'sage.rings.padics.padic_capped_absolute_element.pAdicConvert_CA_ZZ'>
            sage: f.category()
            Category of homsets of sets
        """
        if R.degree() > 1 or R.characteristic() != 0 or R.residue_characteristic() == 0:
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
        cdef Integer ans = PY_NEW(Integer)
        cdef CAElement x = _x
        cconv_mpz_t_out(ans.value, x.value, 0, x.absprec, x.prime_pow)
        return ans

cdef class pAdicConvert_QQ_CA(Morphism):
    """
    The inclusion map from the rationals to a capped absolute ring that is
    defined on all elements with non-negative `p`-adic valuation.

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
            <type 'sage.rings.padics.padic_capped_absolute_element.pAdicConvert_QQ_CA'>
        """
        Morphism.__init__(self, Hom(QQ, R, SetsWithPartialMaps()))
        self._zero = R._element_constructor(R, 0)

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpCA(5).convert_map_from(QQ)
            sage: g = copy(f) # indirect doctest
            sage: g == f # todo: comparison not implemented
            True
            sage: g(1/6)
            1 + 4*5 + 4*5^3 + 4*5^5 + 4*5^7 + 4*5^9 + 4*5^11 + 4*5^13 + 4*5^15 + 4*5^17 + 4*5^19 + O(5^20)
            sage: g(1/6) == f(1/6)
            True
        """
        _slots['_zero'] = self._zero
        return Morphism._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpCA(5).convert_map_from(QQ)
            sage: g = copy(f) # indirect doctest
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
        cconv_mpq_t(ans.value, (<Rational>x).value, ans.prime_pow.prec_cap, True, ans.prime_pow)
        ans.absprec = ans.prime_pow.prec_cap
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in (relative or absolute or both).

        See the documentation for :meth:`pAdicCappedAbsoluteElement.__init__` for more details.

        EXAMPLES::

            sage: R = ZpCA(5,4)
            sage: type(R(10/3,2))
            <type 'sage.rings.padics.padic_capped_absolute_element.pAdicCappedAbsoluteElement'>
            sage: R(10/3,2)
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
            if aprec >= self._zero.prime_pow.prec_cap:
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

cdef class pAdicCoercion_CA_frac_field(RingHomomorphism_coercion):
    """
    The canonical inclusion of Zq into its fraction field.

    EXAMPLES::

        sage: R.<a> = ZqCA(27, implementation='FLINT')
        sage: K = R.fraction_field()
        sage: K.coerce_map_from(R)
        Ring Coercion morphism:
          From: Unramified Extension of 3-adic Ring with capped absolute precision 20 in a defined by (1 + O(3^20))*x^3 + (O(3^20))*x^2 + (2 + O(3^20))*x + (1 + O(3^20))
          To:   Unramified Extension of 3-adic Field with capped relative precision 20 in a defined by (1 + O(3^20))*x^3 + (O(3^20))*x^2 + (2 + O(3^20))*x + (1 + O(3^20))
    """
    def __init__(self, R, K):
        """
        Initialization.

        EXAMPLES::

            sage: R.<a> = ZqCA(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R); type(f)
            <type 'sage.rings.padics.qadic_flint_CA.pAdicCoercion_CA_frac_field'>
        """
        RingHomomorphism_coercion.__init__(self, R.Hom(K), check=False)
        self._zero = K(0)
        self._section = pAdicConvert_CA_frac_field(K, R)

    cpdef Element _call_(self, _x):
        """
        Evaluation.

        EXAMPLES::

            sage: R.<a> = ZqCA(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f(a)
            a + O(3^20)
        """
        cdef CAElement x = _x
        cdef CRElement ans = self._zero._new_c()
        ans.ordp = cremove(ans.unit, x.value, x.absprec, x.prime_pow)
        ans.relprec = x.absprec - ans.ordp
        return ans

    cpdef Element _call_with_args(self, _x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in
        (relative or absolute or both).

        See the documentation for
        :meth:`pAdicCappedAbsoluteElement.__init__` for more details.

        EXAMPLES::

            sage: R.<a> = ZqCA(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f(a, 3)
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
        cdef bint reduce = False
        _process_args_and_kwds(&aprec, &rprec, args, kwds, False, ans.prime_pow)
        if x.absprec < aprec:
            aprec = x.absprec
            reduce = True
        ans.ordp = cremove(ans.unit, x.value, aprec, x.prime_pow)
        ans.relprec = aprec - ans.ordp
        if rprec < ans.relprec:
            ans.relprec = rprec
            reduce = True
        if ans.relprec < 0:
            ans.relprec = 0
            ans.ordp = aprec
            csetzero(ans.unit, x.prime_pow)
        elif reduce:
            creduce(ans.unit, ans.unit, ans.relprec, x.prime_pow)
        return ans

    def section(self):
        """
        Returns a map back to the ring that converts elements of
        non-negative valuation.

        EXAMPLES::

            sage: R.<a> = ZqCA(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f(K.gen())
            a + O(3^20)
        """
        return self._section

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        TESTS::

            sage: R.<a> = ZqCA(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Ring Coercion morphism:
              From: Unramified Extension of 3-adic Ring with capped absolute precision 20 in a defined by (1 + O(3^20))*x^3 + (O(3^20))*x^2 + (2 + O(3^20))*x + (1 + O(3^20))
              To:   Unramified Extension of 3-adic Field with capped relative precision 20 in a defined by (1 + O(3^20))*x^3 + (O(3^20))*x^2 + (2 + O(3^20))*x + (1 + O(3^20))
            sage: g == f
            True
            sage: g is f
            False
            sage: g(a)
            a + O(3^20)
            sage: g(a) == f(a)
            True

        """
        _slots['_zero'] = self._zero
        _slots['_section'] = self._section
        return RingHomomorphism_coercion._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        TESTS::

            sage: R.<a> = ZqCA(9, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Ring Coercion morphism:
              From: Unramified Extension of 3-adic Ring with capped absolute precision 20 in a defined by (1 + O(3^20))*x^2 + (2 + O(3^20))*x + (2 + O(3^20))
              To:   Unramified Extension of 3-adic Field with capped relative precision 20 in a defined by (1 + O(3^20))*x^2 + (2 + O(3^20))*x + (2 + O(3^20))
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
        RingHomomorphism_coercion._update_slots(self, _slots)

cdef class pAdicConvert_CA_frac_field(Morphism):
    """
    The section of the inclusion from `\ZZ_q`` to its fraction field.

    EXAMPLES::

        sage: R.<a> = ZqCA(27, implementation='FLINT')
        sage: K = R.fraction_field()
        sage: f = R.convert_map_from(K); f
        Generic morphism:
          From: Unramified Extension of 3-adic Field with capped relative precision 20 in a defined by (1 + O(3^20))*x^3 + (O(3^20))*x^2 + (2 + O(3^20))*x + (1 + O(3^20))
          To:   Unramified Extension of 3-adic Ring with capped absolute precision 20 in a defined by (1 + O(3^20))*x^3 + (O(3^20))*x^2 + (2 + O(3^20))*x + (1 + O(3^20))
    """
    def __init__(self, K, R):
        """
        Initialization.

        EXAMPLES::

            sage: R.<a> = ZqCA(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K); type(f)
            <type 'sage.rings.padics.qadic_flint_CA.pAdicConvert_CA_frac_field'>
        """
        Morphism.__init__(self, Hom(K, R, SetsWithPartialMaps()))
        self._zero = R(0)

    cpdef Element _call_(self, _x):
        """
        Evaluation.

        EXAMPLES::

            sage: R.<a> = ZqCA(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K)
            sage: f(K.gen())
            a + O(3^20)
        """
        cdef CRElement x = _x
        if x.ordp < 0: raise ValueError("negative valuation")
        cdef CAElement ans = self._zero._new_c()
        cdef bint reduce = False
        ans.absprec = x.relprec + x.ordp
        if ans.absprec > ans.prime_pow.ram_prec_cap:
            ans.absprec = ans.prime_pow.ram_prec_cap
            reduce = True
        if x.ordp >= ans.absprec:
            csetzero(ans.value, ans.prime_pow)
        else:
            cshift(ans.value, x.unit, x.ordp, ans.absprec, ans.prime_pow, reduce)
        return ans

    cpdef Element _call_with_args(self, _x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in
        (relative or absolute or both).

        See the documentation for
        :meth:`pAdicCappedAbsoluteElement.__init__` for more details.

        EXAMPLES::

            sage: R.<a> = ZqCA(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K); a = K(a)
            sage: f(a, 3)
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
        if x.ordp < 0: raise ValueError("negative valuation")
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
            sig_on()
            cshift(ans.value, x.unit, x.ordp, ans.absprec, ans.prime_pow, reduce)
            sig_off()
        return ans

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        TESTS::

            sage: R.<a> = ZqCA(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K)
            sage: a = K(a)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Generic morphism:
              From: Unramified Extension of 3-adic Field with capped relative precision 20 in a defined by (1 + O(3^20))*x^3 + (O(3^20))*x^2 + (2 + O(3^20))*x + (1 + O(3^20))
              To:   Unramified Extension of 3-adic Ring with capped absolute precision 20 in a defined by (1 + O(3^20))*x^3 + (O(3^20))*x^2 + (2 + O(3^20))*x + (1 + O(3^20))
            sage: g == f
            True
            sage: g is f
            False
            sage: g(a)
            a + O(3^20)
            sage: g(a) == f(a)
            True

        """
        _slots['_zero'] = self._zero
        return Morphism._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        TESTS::

            sage: R.<a> = ZqCA(9, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K)
            sage: a = f(a)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Generic morphism:
              From: Unramified Extension of 3-adic Field with capped relative precision 20 in a defined by (1 + O(3^20))*x^2 + (2 + O(3^20))*x + (2 + O(3^20))
              To:   Unramified Extension of 3-adic Ring with capped absolute precision 20 in a defined by (1 + O(3^20))*x^2 + (2 + O(3^20))*x + (2 + O(3^20))
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
    """
    Unpickle capped absolute elements.

    INPUT:

    - ``cls`` -- the class of the capped absolute element.

    - ``parent`` -- the parent, a `p`-adic ring

    - ``value`` -- a Python object wrapping a celement, of the kind
      accepted by the cunpickle function.

    - ``absprec`` -- a Python int or Sage integer.

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
    cconstruct(ans.value, ans.prime_pow)
    cunpickle(ans.value, value, ans.prime_pow)
    ans.absprec = absprec
    return ans
