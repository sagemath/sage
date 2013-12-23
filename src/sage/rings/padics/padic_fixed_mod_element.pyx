"""
p-Adic Fixed-Mod Element

Elements of p-Adic Rings with Fixed Modulus

AUTHORS:

- David Roe
- Genya Zaytman: documentation
- David Harvey: doctests
"""

#*****************************************************************************
#       Copyright (C) 2007 David Roe <roed@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/libs/ntl/decl.pxi"
include "sage/ext/gmp.pxi"
include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"

cimport sage.rings.padics.padic_generic_element
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement

cimport sage.rings.integer
cimport sage.rings.rational
cimport sage.rings.padics.pow_computer
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.padics.pow_computer cimport PowComputer_base
from sage.rings.padics.padic_printing cimport pAdicPrinter_class

#import sage.rings.padics.padic_ring_generic_element
#import sage.rings.padics.padic_field_generic_element
#import sage.rings.padics.padic_lazy_element
import sage.rings.finite_rings.integer_mod
import sage.rings.integer
import sage.rings.rational

from sage.rings.infinity import infinity
from sage.rings.finite_rings.integer_mod import Mod
from sage.rings.padics.precision_error import PrecisionError

#pAdicLazyElement = sage.rings.padics.padic_lazy_element.pAdicLazyElement
#pAdicGenericElement = sage.rings.padics.padic_generic_element.pAdicGenericElement

from sage.libs.pari.gen cimport gen as pari_gen
from sage.interfaces.gp import GpElement

cdef class pAdicFixedModElement(pAdicBaseGenericElement):
    def __init__(pAdicFixedModElement self, parent, x, absprec = None, relprec = None, empty = False):
        r"""
        INPUT:
            parent -- a pAdicRingFixedMod object.

        Types currently supported:
            Integers
            Rationals -- denominator must be relatively prime to p
            FixedMod p-adics

        Types that should be supported:
            Finite precision p-adics
            Lazy p-adics
            Elements of local extensions of THIS p-adic ring that actually lie in Zp
            Elements of IntegerModRing(p^k) for k less than or equal to the modulus

        EXAMPLES:
            sage: R = Zp(5, 20, 'fixed-mod', 'terse')

        Construct from integers:
            sage: R(3)
            3 + O(5^20)
            sage: R(75)
            75 + O(5^20)
            sage: R(0)
            0 + O(5^20)

            sage: R(-1)
            95367431640624 + O(5^20)
            sage: R(-5)
            95367431640620 + O(5^20)

        Construct from rationals:
            sage: R(1/2)
            47683715820313 + O(5^20)
            sage: R(-7875/874)
            9493096742250 + O(5^20)
            sage: R(15/425)
            Traceback (most recent call last):
            ...
            ValueError: p divides denominator

        # todo: the above error message does not agree with the error message
        # in the corresponding capped-relative constructor

        Construct from IntegerMod:
            sage: R(Integers(125)(3))
            3 + O(5^20)
            sage: R(Integers(5)(3))
            3 + O(5^20)
            sage: R(Integers(5^30)(3))
            3 + O(5^20)
            sage: R(Integers(5^30)(1+5^23))
            1 + O(5^20)
            sage: R(Integers(49)(3))
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce from the given integer mod ring (not a power of the same prime)

            sage: R(Integers(48)(3))
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce from the given integer mod ring (not a power of the same prime)

        Some other conversions:
            sage: R(R(5))
            5 + O(5^20)

        # todo: doctests for converting from other types of p-adic rings

        """
        mpz_init(self.value)
        pAdicBaseGenericElement.__init__(self,parent)
        if empty:
            return
        cdef Integer tmp
        if PY_TYPE_CHECK(x, pAdicGenericElement):
            if x.valuation() < 0:
                raise ValueError, "element has negative valuation"
            if parent.prime() != x.parent().prime():
                raise TypeError, "Cannot coerce between p-adic parents with different primes."
        #if PY_TYPE_CHECK(x, pAdicLazyElement):
        #    try:
        #        x.set_precision_absolute(absprec)
        #    except PrecisionError:
        #        pass
        #    if mpz_cmp((<pAdicLazyElement>x).value, self.prime_pow.top_power) >= 0:
        #        mpz_mod(self.value, (<pAdicLazyElement>x).value, self.prime_pow.top_power)
        #    else:
        #        mpz_set(self.value, (<pAdicLazyElement>x).value)
        #    return
        if PY_TYPE_CHECK(x, pAdicBaseGenericElement):
            tmp = <Integer> x.lift()
            if mpz_cmp(tmp.value, self.prime_pow.pow_mpz_t_top()[0]) >= 0:
                mpz_mod(self.value, tmp.value, self.prime_pow.pow_mpz_t_top()[0])
            else:
                mpz_set(self.value, tmp.value)
            return

        if isinstance(x, pari_gen) or isinstance(x, GpElement):
            if isinstance(x, GpElement):
                x = x._pari_()
            if x.type() == "t_PADIC":
                if x.variable() != self.prime_pow.prime:
                    raise TypeError, "Cannot coerce a pari p-adic with the wrong prime."
                x = x.lift()
            if x.type() == 't_INT':
                x = Integer(x)
            elif x.type() == 't_FRAC':
                x = Rational(x)
            else:
                raise TypeError, "unsupported coercion from pari: only p-adics, integers and rationals allowed"

        if sage.rings.finite_rings.integer_mod.is_IntegerMod(x):
            if (<Integer>x.modulus())._is_power_of(<Integer>parent.prime()):
                x = x.lift()
            else:
                raise TypeError, "cannot coerce from the given integer mod ring (not a power of the same prime)"

        #if sage.rings.finite_field_element.is_FiniteFieldElement(x):
        #    if x.parent().order() != parent.prime():
        #        raise TypeError, "can only create p-adic element out of finite field when order of field is p"
        #    x = x.lift()

        #Now use the code below to convert from integer or rational, so don't make the next line elif

        if PY_TYPE_CHECK(x, Integer):
            self._set_from_mpz((<Integer>x).value)
        elif isinstance(x, Rational):
            if self.prime_pow.prime.divides(x.denominator()):
                raise ValueError, "p divides the denominator"
            else:
                tmp = <Integer> x % parent.prime_pow(parent.precision_cap())
                self._set_from_mpz(tmp.value)
        elif isinstance(x, (int, long)):
            tmp = <Integer> Integer(x)
            self._set_from_mpz(tmp.value)
        else:
            self._set_from_mpq((<Rational>Rational(x)).value)

    def __dealloc__(self):
        """
        Deallocation.

        TESTS::

            sage: R = ZpFM(5)
            sage: a = R(17)
            sage: del(a)
        """
        mpz_clear(self.value)

    def __reduce__(self):
        """
        Pickling.

        EXAMPLES::

            sage: a = ZpFM(5)(-3)
            sage: type(a)
            <type 'sage.rings.padics.padic_fixed_mod_element.pAdicFixedModElement'>
            sage: loads(dumps(a)) == a
            True
        """
        return make_pAdicFixedModElement, (self.parent(), self.lift())

    cdef int _set_from_mpz(pAdicFixedModElement self, mpz_t value) except -1:
        """
        Sets self from an mpz_t.

        TESTS::

            sage: R = ZpFM(5)
            sage: a = R(17,5); a #indirect doctest
            2 + 3*5 + O(5^20)
        """
        if mpz_sgn(value) == -1 or mpz_cmp(value, self.prime_pow.pow_mpz_t_top()[0]) >= 0:
            sig_on()
            mpz_mod(self.value, value, self.prime_pow.pow_mpz_t_top()[0])
            sig_off()
        else:
            mpz_set(self.value, value)
        return 0

    cdef int _set_from_mpq(pAdicFixedModElement self, mpq_t x) except -1:
        """
        Sets self from an mpq_t.

        TESTS::

            sage: R = ZpFM(5,5)
            sage: a = R(25/9); a #indirect doctest
            4*5^2 + 2*5^3 + O(5^5)
        """
        if mpz_divisible_p(mpq_denref(x), self.prime_pow.prime.value):
            raise ValueError, "p divides denominator"
        sig_on()
        mpz_invert(self.value, mpq_denref(x), self.prime_pow.pow_mpz_t_top()[0])
        mpz_mul(self.value, self.value, mpq_numref(x))
        mpz_mod(self.value, self.value, self.prime_pow.pow_mpz_t_top()[0])
        sig_off()
        return 0

    cdef int _set_mpz_into(pAdicFixedModElement self, mpz_t dest) except -1:
        """
        Sets dest to a lift of self.

        TESTS::

            sage: R = ZpCA(5); S.<a> = ZqFM(25)
            sage: S(R(17))
            2 + 3*5 + O(5^20)
        """
        mpz_set(dest, self.value)
        return 0

    cdef int _set_mpq_into(pAdicFixedModElement self, mpq_t dest) except -1:
        """
        Sets dest to a lift of self.

        Not currently used internally.
        """
        mpq_set_z(dest, self.value)
        return 0

    cdef pAdicFixedModElement _new_c(self):
        """
        Creates a new element with the same basic info.

        TESTS::

            sage: R = ZpFM(5); R(6) * R(7) #indirect doctest
            2 + 3*5 + 5^2 + O(5^20)
        """
        cdef pAdicFixedModElement x
        x = PY_NEW(pAdicFixedModElement)
        x._parent = self._parent
        mpz_init(x.value)
        x.prime_pow = self.prime_pow
        return x

    cpdef bint _is_inexact_zero(self) except -1:
        """
        Returns True if self is indistinguishable from zero.

        EXAMPLES:
            sage: R = ZpFM(7, 5)
            sage: R(14)._is_inexact_zero()
            False
            sage: R(0)._is_inexact_zero()
            True
        """
        return mpz_sgn(self.value) == 0

    def __richcmp__(left, right, op):
        """
        Comparison.

        TESTS::

            sage: R = ZpFM(5)
            sage: a = R(17)
            sage: b = R(21)
            sage: a == b
            False
            sage: a < b
            True
        """
        return (<Element>left)._richcmp(right, op)

    def __invert__(self):
        r"""
        Returns multiplicative inverse of this element. The valuation
        of self must be zero.

        EXAMPLES::

            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: ~R(2)
            4 + 3*7 + 3*7^2 + 3*7^3 + O(7^4)
            sage: ~R(0)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
            sage: ~R(7)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
        """
        return self._invert_c_impl()

    cpdef RingElement _invert_c_impl(self):
        """
        Returns multiplicative inverse of this element. The valuation
        of self must be zero.

        EXAMPLES::

            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: ~R(2) # indirect doctest
            4 + 3*7 + 3*7^2 + 3*7^3 + O(7^4)
            sage: ~R(0)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
            sage: ~R(7)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
        """
        cdef pAdicFixedModElement ans
        if mpz_divisible_p(self.value, self.prime_pow.prime.value) != 0:
            raise ValueError, "cannot invert non-unit"
        else:
            ans = self._new_c()
            sig_on()
            mpz_invert(ans.value, self.value, self.prime_pow.pow_mpz_t_top()[0])
            sig_off()
            return ans

    cdef pAdicFixedModElement _lshift_c(pAdicFixedModElement self, long shift):
        """
        Multiplies self by p^shift.

        If shift < -self.ordp(), digits will be truncated.  See __rshift__ for details.

        EXAMPLES::

            sage: R = ZpFM(5); a = R(17); a << 2 #indirect doctest
            2*5^2 + 3*5^3 + O(5^20)
        """
        cdef pAdicFixedModElement ans
        cdef unsigned long prec_cap
        if shift < 0:
            return self._rshift_c(-shift)
        prec_cap = self.prime_pow.prec_cap
        if shift >= prec_cap:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            return ans
        elif shift > 0:
            ans = self._new_c()
            if mpz_cmp(self.value, self.prime_pow.pow_mpz_t_tmp(prec_cap - shift)[0]) >= 0:
                mpz_mod(ans.value, self.value, self.prime_pow.pow_mpz_t_tmp(prec_cap - shift)[0])
            else:
                mpz_set(ans.value, self.value)
            mpz_mul(ans.value, ans.value, self.prime_pow.pow_mpz_t_tmp(shift)[0])
            return ans
        else:
            return self

    def __lshift__(pAdicFixedModElement self, shift):
        """
        Multiplies self by p^shift.

        If shift < -self.ordp(), digits will be truncated.  See __rshift__ for details.

        EXAMPLES::

            We create a fixed modulus ring:
            sage: R = ZpFM(5, 20); a = R(1000); a
            3*5^3 + 5^4 + O(5^20)

            Shifting to the right is the same as dividing by a power of
            the uniformizer $p$ of the $p$-adic ring.
            sage: a >> 1
            3*5^2 + 5^3 + O(5^20)

            Shifting to the left is the same as multiplying by a power of $p$:
            sage: a << 2
            3*5^5 + 5^6 + O(5^20)
            sage: a*5^2
            3*5^5 + 5^6 + O(5^20)

            Shifting by a negative integer to the left is the same as right shifting
            by the absolute value:
            sage: a << -3
            3 + 5 + O(5^20)
            sage: a >> 3
            3 + 5 + O(5^20)
        """
        cdef pAdicFixedModElement ans
        if not PY_TYPE_CHECK(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_slong_p((<Integer>shift).value) == 0:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            return ans
        return self._lshift_c(mpz_get_si((<Integer>shift).value))

    cdef pAdicFixedModElement _rshift_c(pAdicFixedModElement self, long shift):
        """
        Divides by p^shift and truncates.

        Note that this operation loses precision if shift > 0.

        EXAMPLES::

            sage: R = ZpFM(5); a = R(77); a >> 1 #indirect doctest
            3*5 + O(5^20)
        """
        cdef pAdicFixedModElement ans
        cdef unsigned long prec_cap
        if shift < 0:
            return self._lshift_c(-shift)
        prec_cap = self.prime_pow.prec_cap
        if shift >= prec_cap:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            return ans
        elif shift > 0:
            ans = self._new_c()
            mpz_fdiv_q(ans.value, self.value, self.prime_pow.pow_mpz_t_tmp(shift)[0])
            return ans
        else:
            return self

    def __rshift__(pAdicFixedModElement self, shift):
        """
        Divides by p^shift, and truncates.

        Note that this operation will insert arbitrary digits (in practice, currently all zero) in the least significant digits.

        EXAMPLES::

            sage: R = ZpFM(997, 7); a = R(123456878908); a
            964*997 + 572*997^2 + 124*997^3 + O(997^7)

            Shifting to the right divides by a power of p, but dropping terms with
            negative valuation:
            sage: a >> 3
            124 + O(997^7)

            A negative shift multiplies by that power of p.
            sage: a >> -3
            964*997^4 + 572*997^5 + 124*997^6 + O(997^7)
        """
        cdef pAdicFixedModElement ans
        if not PY_TYPE_CHECK(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_slong_p((<Integer>shift).value) == 0:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            return ans
        return self._rshift_c(mpz_get_si((<Integer>shift).value))

    cpdef ModuleElement _neg_(self):
        r"""
        Returns negative of self.

        EXAMPLES::

            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: -R(7) #indirect doctest
            6*7 + 6*7^2 + 6*7^3 + O(7^4)
        """
        if mpz_sgn(self.value) == 0:
            return self
        cdef pAdicFixedModElement ans
        ans = self._new_c()
        mpz_sub(ans.value, self.prime_pow.pow_mpz_t_top()[0], self.value)
        return ans

    def __pow__(pAdicFixedModElement self, right, m): # NOTE: m ignored, always use self.prime_pow.pow_mpz_t_top()[0]
        """
        Exponentiation.

        EXAMPLES::

            sage: R = ZpFM(11, 5)
            sage: R(1/2)^5
            10 + 7*11 + 11^2 + 5*11^3 + 4*11^4 + O(11^5)
            sage: R(1/32)
            10 + 7*11 + 11^2 + 5*11^3 + 4*11^4 + O(11^5)
            sage: R(1/2)^5 == R(1/32)
            True
            sage: R(3)^1000 #indirect doctest
            1 + 4*11^2 + 3*11^3 + 7*11^4 + O(11^5)

        TESTS:

        We define ``0^0`` to be unity, :trac:`13786`::

            sage: R = ZpFM(11, 5)
            sage: R(0)^0
            1 + O(11^5)
            sage: R(0)^0 == R(1)
            True

        The value returned from ``0^0`` should belong to our ring::

            sage: R = ZpFM(11, 5)
            sage: type(R(0)^0) == type(R(0))
            True

        Test that #3865 is fixed::

            sage: R(gp('5 + O(11^2)'))
            5 + O(11^5)

        """
        if not PY_TYPE_CHECK(right, Integer):
            right = Integer(right) #Need to make sure that this works for p-adic exponents
        if right == 0 and self == 0:
            return self.parent(1)
        cdef pAdicFixedModElement ans
        ans = self._new_c()
        sig_on()
        mpz_powm(ans.value, self.value, (<Integer>right).value, self.prime_pow.pow_mpz_t_top()[0])
        sig_off()
        return ans

    cpdef ModuleElement _add_(self, ModuleElement right):
        r"""
        Returns sum of self and right.

        EXAMPLES::

            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: x = R(1721); x
            6 + 5*7^3 + O(7^4)
            sage: y = R(1373); y
            1 + 4*7^3 + O(7^4)
            sage: x + y #indirect doctest
            7 + 2*7^3 + O(7^4)
        """
        cdef pAdicFixedModElement ans
        ans = self._new_c()
        mpz_add(ans.value, self.value, (<pAdicFixedModElement>right).value)
        if mpz_cmp(ans.value, self.prime_pow.pow_mpz_t_top()[0]) >= 0:
            mpz_sub(ans.value, ans.value, self.prime_pow.pow_mpz_t_top()[0])
        return ans

    cpdef RingElement _mul_(self, RingElement right):
        r"""
        Returns product of self and right.

        EXAMPLES::

            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: R(3) * R(2) #indirect doctest
            6 + O(7^4)
            sage: R(1/2) * R(2)
            1 + O(7^4)
        """
        cdef pAdicFixedModElement ans
        ans = self._new_c()
        mpz_mul(ans.value, self.value, (<pAdicFixedModElement>right).value)
        mpz_fdiv_r(ans.value, ans.value, self.prime_pow.pow_mpz_t_top()[0])
        return ans

    cpdef ModuleElement _sub_(self, ModuleElement right):
        r"""
        Returns difference of self and right.

        EXAMPLES:
            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: x = R(1721); x
            6 + 5*7^3 + O(7^4)
            sage: y = R(1373); y
            1 + 4*7^3 + O(7^4)
            sage: x - y #indirect doctest
            5 + 7^3 + O(7^4)
        """
        cdef pAdicFixedModElement ans
        ans = self._new_c()
        mpz_sub(ans.value, self.value, (<pAdicFixedModElement>right).value)
        if mpz_sgn(ans.value) == -1:
            mpz_add(ans.value, ans.value, self.prime_pow.pow_mpz_t_top()[0])
        return ans

    cpdef RingElement _div_(self, RingElement right):
        r"""
        Returns quotient of self and right. The latter must have
        valuation zero.

        EXAMPLES:
            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: R(3) / R(2) #indirect doctest
            5 + 3*7 + 3*7^2 + 3*7^3 + O(7^4)
            sage: R(5) / R(0)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
            sage: R(7) / R(49)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
        """
        cdef int t
        cdef pAdicFixedModElement ans
        if mpz_divisible_p((<pAdicFixedModElement>right).value, self.prime_pow.prime.value) != 0:
            raise ValueError, "cannot invert non-unit"
        else:
            ans = self._new_c()
            sig_on()
            mpz_invert(ans.value, (<pAdicFixedModElement>right).value, self.prime_pow.pow_mpz_t_top()[0])
            mpz_mul(ans.value, ans.value, self.value)
            mpz_fdiv_r(ans.value, ans.value, self.prime_pow.pow_mpz_t_top()[0])
            sig_off()
            return ans

    def add_bigoh(self, absprec):
        """
        Returns a new element truncated modulo p^absprec.

        INPUT::

            - self -- a p-adic element
            - absprec -- an integer

        OUTPUT::

            - element -- a new element truncated modulo p^absprec.

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod','series'); a = R(8); a.add_bigoh(1)
            1 + O(7^4)
        """
        cdef pAdicFixedModElement ans
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        if mpz_cmp_ui((<Integer>absprec).value, self.prime_pow.prec_cap) >= 0:
            return self
        ans = self._new_c()
        mpz_mod(ans.value, self.value, self.prime_pow.pow_mpz_t_tmp(mpz_get_ui((<Integer>absprec).value))[0])
        return ans

    def __copy__(self):
        """
        Returns a copy of self.

        EXAMPLES::

            sage: a = ZpFM(5,6)(17); b = copy(a)
            sage: a == b
            True
            sage: a is b
            False
        """
        cdef pAdicFixedModElement ans
        ans = self._new_c()
        mpz_set(ans.value, self.value)
        return ans

    def is_zero(self, absprec = None):
        r"""
        Returns whether self is zero modulo $p^{\mbox{absprec}}$.

        INPUT::

            - self -- a p-adic element
            - absprec -- an integer

        OUTPUT::

            boolean -- whether self is zero

        EXAMPLES::

            sage: R = ZpFM(17, 6)
            sage: R(0).is_zero()
            True
            sage: R(17^6).is_zero()
            True
            sage: R(17^2).is_zero(absprec=2)
            True
        """
        if absprec is None:
            return mpz_sgn(self.value) == 0
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        cdef unsigned long aprec
        aprec = mpz_get_ui((<Integer>absprec).value)
        if aprec >= self.prime_pow.prec_cap:
            return mpz_sgn(self.value) == 0
        cdef mpz_t tmp
        mpz_init(tmp)
        mpz_mod(tmp, self.value, self.prime_pow.pow_mpz_t_tmp(aprec)[0])
        if mpz_sgn(tmp) == 0:
            mpz_clear(tmp)
            return True
        else:
            mpz_clear(tmp)
            return False

    def is_equal_to(self, right, absprec = None): #assumes they have the same parent
        r"""
        Returns whether self is equal to right modulo $p^{\mbox{absprec}}$.

        If absprec is None, returns if self == 0.

        INPUT::

            - self -- a p-adic element
            - right -- a p-addic element with the same parent
            - absprec -- a positive integer (or None)

        OUTPUT::

            boolean -- whether self is equal to right

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
        if absprec is None:
            return mpz_cmp(self.value, (<pAdicFixedModElement>right).value) == 0
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        if absprec < 0:
            return True
        cdef unsigned long aprec
        aprec = mpz_get_ui((<Integer>absprec).value)
        if aprec >= self.prime_pow.prec_cap:
            return mpz_cmp(self.value, (<pAdicFixedModElement>right).value) == 0
        cdef mpz_t tmp1, tmp2
        mpz_init(tmp1)
        mpz_init(tmp2)
        mpz_mod(tmp1, self.value, self.prime_pow.pow_mpz_t_tmp(aprec)[0])
        mpz_mod(tmp2, (<pAdicFixedModElement>right).value, self.prime_pow.pow_mpz_t_tmp(aprec)[0])
        if mpz_cmp(tmp1, tmp2) == 0:
            mpz_clear(tmp1)
            mpz_clear(tmp2)
            return True
        else:
            mpz_clear(tmp1)
            mpz_clear(tmp2)
            return False

    def lift(self):
        r"""
        Return an integer congruent to self modulo self's precision.

        INPUT::

            - self -- a p-adic element

        OUTPUT::

            - integer -- a integer congruent to self mod $p^{\mbox{prec}}$

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod'); a = R(8); a.lift()
            8
            sage: type(a.lift())
            <type 'sage.rings.integer.Integer'>
        """
        return self.lift_c()

    cdef Integer lift_c(pAdicFixedModElement self):
        """
        Returns an integer congruent to self modulo self's precision.

        EXAMPLES::

            sage: R = ZpFM(7,4); a = R(8); a.lift() # indirect doctest
            8
        """
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set(ans.value, self.value)
        return ans

    def lift_to_precision(self, absprec=None):
        """
        Returns self.

        For compatibility with other p-adic types.

        EXAMPLES::

            sage: R = ZpFM(5); a = R(5); a.lift_to_precision(7)
            5 + O(5^20)
        """
        return self

    def list(self, lift_mode = 'simple'):
        r"""
        Returns a list of coefficients of p starting with $p^0$.

        INPUT::

            - self -- a p-adic element
            - lift_mode -- 'simple', 'smallest' or 'teichmuller' (default 'simple')

        OUTPUT::

            - list -- the list of coefficients of self

        NOTES::

            Returns a list [a_0, a_1, \ldots, a_n] so that each a_i is an integer
            and \sum_{i = 0}^n a_i * p^i = self, modulo the precision cap.
            If lift_mode = 'simple', 0 <= a_i < p.
            If lift_mode = 'smallest', -p/2 < a_i <= p/2.
            If lift_mode = 'teichmuller', a_i^p = a_i, modulo the precision cap.

        EXAMPLES::

            sage: R = ZpFM(7,6); a = R(12837162817); a
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
            O(7^6),
            5 + 2*7 + 3*7^3 + 6*7^4 + 4*7^5 + O(7^6),
            1 + O(7^6),
            3 + 4*7 + 6*7^2 + 3*7^3 + 2*7^5 + O(7^6),
            5 + 2*7 + 3*7^3 + 6*7^4 + 4*7^5 + O(7^6)]
            sage: sum([L[i] * 7^i for i in range(len(L))])
            3 + 4*7 + 4*7^2 + 4*7^4 + O(7^6)
        """
        if lift_mode == 'teichmuller':
            return self.teichmuller_list()
        if lift_mode == 'simple':
            return (<pAdicPrinter_class>self.parent()._printer).base_p_list(self.value, True)
        elif lift_mode == 'smallest':
            return (<pAdicPrinter_class>self.parent()._printer).base_p_list(self.value, False)
        else:
            raise ValueError

    cdef object teichmuller_list(pAdicFixedModElement self):
        r"""
        Returns a list [$a_0$, $a_1$,..., $a_n$] such that
            - $a_i^p = a_i$
            - self.unit_part() = $\sum_{i = 0}^n a_i p^i$

        EXAMPLES::

            sage: R = ZpCA(5,5); R(14).list('teichmuller') #indirect doctest
            [4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5),
            3 + 3*5 + 2*5^2 + 3*5^3 + O(5^4),
            2 + 5 + 2*5^2 + O(5^3),
            1 + O(5^2),
            4 + O(5)]
        """
        # May eventually want to add a dict to store teichmuller lifts already seen, if p small enough
        cdef unsigned long curpower, preccap
        cdef mpz_t tmp, tmp2
        cdef pAdicFixedModElement list_elt
        ans = PyList_New(0)
        preccap = self.prime_pow.prec_cap
        curpower = preccap
        mpz_init_set(tmp, self.value)
        mpz_init(tmp2)
        while mpz_sgn(tmp) != 0:
            curpower -= 1
            list_elt = self._new_c()
            mpz_mod(list_elt.value, tmp, self.prime_pow.prime.value)
            mpz_set(tmp2, self.prime_pow.pow_mpz_t_tmp(preccap)[0])
            self.teichmuller_set_c(list_elt.value, tmp2)
            mpz_sub(tmp, tmp, list_elt.value)
            mpz_divexact(tmp, tmp, self.prime_pow.prime.value)
            mpz_mod(tmp, tmp, self.prime_pow.pow_mpz_t_tmp(curpower)[0])
            PyList_Append(ans, list_elt)
        mpz_clear(tmp)
        mpz_clear(tmp2)
        return ans

    def _teichmuller_set(self):
        """
        Sets self to be the Teichmuller representative with the same residue as self.

        WARNING: This function modifies self, which is not safe.  Elements are supposed to be immutable.

        EXAMPLES::

            sage: R = ZpFM(17,5); a = R(11)
            sage: a
            11 + O(17^5)
            sage: a._teichmuller_set(); a
            11 + 14*17 + 2*17^2 + 12*17^3 + 15*17^4 + O(17^5)
            sage: a.list('teichmuller')
            [11 + 14*17 + 2*17^2 + 12*17^3 + 15*17^4 + O(17^5)]
        """
        cdef mpz_t tmp
        mpz_init_set(tmp, self.prime_pow.pow_mpz_t_top()[0])
        self.teichmuller_set_c(self.value, tmp)
        mpz_clear(tmp)


    def multiplicative_order(self):
        r"""
        Returns the minimum possible multiplicative order of self.

        INPUT::

            - self -- a p-adic element

        OUTPUT::

            - integer -- the multiplicative order of self.  This is
              the minimum multiplicative order of all elements of Z_p
              lifting self to infinite precision.

        EXAMPLES::

            sage: R = ZpFM(7, 6)
            sage: R(1/3)
            5 + 4*7 + 4*7^2 + 4*7^3 + 4*7^4 + 4*7^5 + O(7^6)
            sage: R(1/3).multiplicative_order()
            +Infinity
            sage: R(7).multiplicative_order()
            +Infinity
            sage: R(1).multiplicative_order()
            1
            sage: R(-1).multiplicative_order()
            2
            sage: R.teichmuller(3).multiplicative_order()
            6
        """
        cdef mpz_t tmp
        cdef Integer ans
        if mpz_divisible_p(self.value, self.prime_pow.prime.value):
            return infinity
        if mpz_cmp_ui(self.value, 1) == 0:
            ans = PY_NEW(Integer)
            mpz_set_ui(ans.value, 1)
            return ans
        mpz_init(tmp)
        mpz_sub_ui(tmp, self.prime_pow.pow_mpz_t_top()[0], 1)
        if mpz_cmp(self.value, tmp) == 0:
            ans = PY_NEW(Integer)
            mpz_set_ui(ans.value, 2)
            return ans
        # check if self is an approximation to a teichmuller lift:
        mpz_powm(tmp, self.value, self.prime_pow.prime.value, self.prime_pow.pow_mpz_t_top()[0])
        if mpz_cmp(tmp, self.value) == 0:
            mpz_clear(tmp)
            return self.residue(1).multiplicative_order()
        else:
            mpz_clear(tmp)
            return infinity

    def padded_list(self, n, list_mode = 'simple'):
        """
        Returns a list of coefficients of p starting with $p^0$ up to
        $p^n$ exclusive (padded with zeros if needed)

        INPUT::

            - self -- a p-adic element
            - n -- an integer

        OUTPUT::

            - list -- the list of coefficients of self

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod'); a = R(2*7+7**2); a.padded_list(5)
            [0, 2, 1, 0, 0]

        NOTE::

            For elements with positive valuation, this function will
            return a list with leading 0s, unlike for field elements.

            The slice operators throw an error if asked for a slice
            above the precision, while this function works
        """
        if list_mode == 'simple' or list_mode == 'smallest':
            zero = Integer(0)
        else:
            zero = self.parent()(0)
        L = self.list()
        return L[:n] + [zero] * (n - len(L))

    def precision_absolute(self):
        """
        Returns the absolute precision of self.

        INPUT::

            - self -- a p-adic element

        OUTPUT::

            - integer -- the absolute precision of self

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod'); a = R(7); a.precision_absolute()
            4
        """
        return self.parent().precision_cap()

    def precision_relative(self):
        r"""
        Returns the relative precision of self

        INPUT::

            - self -- a p-adic element

        OUTPUT::

            - integer -- the relative precision of self

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod'); a = R(7); a.precision_relative()
            3
            sage: a = R(0); a.precision_relative()
            0
        """
        cdef unsigned long diff
        cdef Integer ans
        ans = PY_NEW(Integer)
        diff = self.prime_pow.prec_cap - self.valuation_c()
        mpz_set_si(ans.value, diff)
        return ans

    def residue(self, absprec=1):
        r"""
        Reduces this mod $p^{\mbox{prec}}$

        INPUT::

            - self -- a p-adic element
            - absprec - an integer (default 1)

        OUTPUT::

            - element of Z/(p^prec Z) -- self reduced mod p^prec

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod'); a = R(8); a.residue(1)
            1
        """
        cdef Integer selfvalue, modulus
        selfvalue = PY_NEW(Integer)
        modulus = PY_NEW(Integer)
        mpz_set(selfvalue.value, self.value)
        cdef unsigned long aprec
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        if mpz_fits_ulong_p((<Integer>absprec).value) == 0:
            raise ValueError, "When calling residue, use the exponent of p, not the integer p^exp."
        else:
            aprec = mpz_get_ui((<Integer>absprec).value)
        if aprec > self.prime_pow.prec_cap:
            sig_on()
            mpz_pow_ui(modulus.value, self.prime_pow.prime.value, aprec)
            sig_off()
        else:
            mpz_set(modulus.value, self.prime_pow.pow_mpz_t_tmp(aprec)[0])
        return Mod(selfvalue, modulus)

    #def square_root(self):
    #    r"""
    #    Returns the square root of this p-adic number

    #    INPUT:
    #        self -- a p-adic element
    #    OUTPUT:
    #        p-adic element -- the square root of this p-adic number

    #        The square root chosen is the one whose reduction mod p is in
    #        the range [0, p/2).

    #        Note that because this is a fixed modulus ring, garbage digits
    #        may be introduced, if either
    #        (a) the valuation of the input is positive, or
    #        (b) p = 2.

    #        If no square root exists, a ValueError is raised.
    #        (This may be changed later to return an element of an extension
    #        field.)

    #    EXAMPLES:
    #        sage: R = Zp(3,20,'fixed-mod')
    #        sage: R(0).square_root()
    #            O(3^20)
    #        sage: R(1).square_root()
    #            1 + O(3^20)
    #        sage: R(2).square_root()
    #        Traceback (most recent call last):
    #        ...
    #        ValueError: element is not a square
    #        sage: R(4).square_root() == R(-2)
    #            True
    #        sage: R(9).square_root()
    #            3 + O(3^20)
    #        sage: R2 = Zp(2,20,'fixed-mod')
    #        sage: R2(0).square_root()
    #            O(2^20)
    #        sage: R2(1).square_root()
    #            1 + O(2^20)
    #        sage: R2(4).square_root()
    #            2 + O(2^20)
    #        sage: R2(9).square_root() == R2(3) or R2(9).square_root() == R2(-3)
    #            True
    #        sage: R2(17).square_root()
    #            1 + 2^3 + 2^5 + 2^6 + 2^7 + 2^9 + 2^10 + 2^13 + 2^16 + 2^17 + O(2^20)
    #        sage: R3 = Zp(5,20,'fixed-mod', 'terse')
    #        sage: R3(0).square_root()
    #            0 + O(5^20)
    #        sage: R3(1).square_root()
    #            1 + O(5^20)
    #        sage: R3(-1).square_root() == R3.teichmuller(2) or R3(-1).square_root() == R3.teichmuller(3)
    #            True
    #    """
    #    #todo: make more efficient
    #    try:
    #        # use pari
    #        return self.parent()(pari(self).sqrt())
    #    except PariError:
    #        # todo: should eventually change to return an element of
    #        # an extension field
    #        raise ValueError, "element is not a square"

    cpdef pAdicFixedModElement unit_part(pAdicFixedModElement self):
        r"""
        Returns the unit part of self.

        If the valuation of self is positive, then the high digits of the
        result will be zero.

        INPUT::

            - self -- a p-adic element

        OUTPUT::

            - p-adic element -- the unit part of self

        EXAMPLES::

            sage: R = Zp(17, 4, 'fixed-mod')
            sage: R(5).unit_part()
            5 + O(17^4)
            sage: R(18*17).unit_part()
            1 + 17 + O(17^4)
            sage: R(0).unit_part()
            O(17^4)
            sage: type(R(5).unit_part())
            <type 'sage.rings.padics.padic_fixed_mod_element.pAdicFixedModElement'>
            sage: R = ZpFM(5, 5); a = R(75); a.unit_part()
            3 + O(5^5)
        """
        cdef pAdicFixedModElement ans
        if mpz_sgn(self.value) == 0:
            return self
        elif mpz_divisible_p(self.value, self.prime_pow.prime.value):
            ans = self._new_c()
            mpz_remove(ans.value, self.value, self.prime_pow.prime.value)
            return ans
        else:
            return self

    def valuation(self):
        """
        Returns the valuation of self.

        If self is zero, the valuation returned is the precision of the ring.

        INPUT::

            - self -- a p-adic element

        OUTPUT::

            - integer -- the valuation of self.

        EXAMPLES::

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
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set_si(ans.value, self.valuation_c())
        return ans

    cdef long valuation_c(self):
        """
        Returns the valuation of self.

        EXAMPLES::

            sage: R = ZpFM(5, 5); R(0).valuation() #indirect doctest
            5
        """
        if mpz_sgn(self.value) == 0:
            return self.prime_pow.prec_cap
        cdef mpz_t tmp
        cdef long ans
        mpz_init(tmp)
        ans = mpz_remove(tmp, self.value, self.prime_pow.prime.value)
        mpz_clear(tmp)
        return ans

    cpdef val_unit(self):
        """
        Returns a 2-tuple, the first element set to the valuation of
        self, and the second to the unit part of self.

        If self == 0, then the unit part is O(p^self.parent().precision_cap()).

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: a = R(75); b = a - a
            sage: a.val_unit()
            (2, 3 + O(5^5))
            sage: b.val_unit()
            (5, O(5^5))
        """
        cdef Integer val
        cdef pAdicFixedModElement unit
        if mpz_sgn(self.value) == 0:
            return (self.parent().precision_cap(), self)
        val = PY_NEW(Integer)
        unit = self._new_c()
        mpz_set_ui(val.value, mpz_remove(unit.value, self.value, self.prime_pow.prime.value))
        return (val, unit)

    def __hash__(self):
        """
        Hashing.

        EXAMPLES::

            sage: R = ZpCA(11, 5)
            sage: hash(R(3)) == hash(3)
            True
        """
        return hash(self.lift_c())

def make_pAdicFixedModElement(parent, value):
    """
    Unpickles a capped relative element.

    EXAMPLES::

        sage: from sage.rings.padics.padic_fixed_mod_element import make_pAdicFixedModElement
        sage: R = ZpFM(5)
        sage: a = make_pAdicFixedModElement(R, 17*25); a
        2*5^2 + 3*5^3 + O(5^20)
    """
    return parent(value)

