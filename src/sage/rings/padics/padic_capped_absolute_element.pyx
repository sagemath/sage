"""
`p`-Adic Capped Absolute Element

Elements of `p`-Adic Rings with Absolute Precision Cap

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

cimport sage.rings.rational
cimport sage.rings.padics.pow_computer
from sage.rings.padics.pow_computer cimport PowComputer_base
from sage.rings.padics.padic_printing cimport pAdicPrinter_class
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement

import sage.rings.padics.padic_generic_element
import sage.rings.finite_rings.integer_mod
import sage.rings.integer
import sage.rings.rational

from sage.rings.infinity import infinity
from sage.rings.finite_rings.integer_mod import Mod
from sage.rings.padics.precision_error import PrecisionError

from sage.libs.pari.gen cimport gen as pari_gen
gp_element = sage.interfaces.gp.GpElement

cdef class pAdicCappedAbsoluteElement(pAdicBaseGenericElement):
    def __init__(pAdicCappedAbsoluteElement self, parent, x, absprec=infinity, relprec = infinity, empty=False):
        """
        EXAMPLES::

            sage: R = ZpCA(3, 5)
            sage: R(2)
            2 + O(3^5)
            sage: R(2, absprec=2)
            2 + O(3^2)
            sage: R(3, relprec=2)
            3 + O(3^3)
            sage: R(Qp(3)(10))
            1 + 3^2 + O(3^5)
            sage: R(pari(6))
            2*3 + O(3^5)
            sage: R(pari(1/2))
            2 + 3 + 3^2 + 3^3 + 3^4 + O(3^5)
            sage: R(1/2)
            2 + 3 + 3^2 + 3^3 + 3^4 + O(3^5)
            sage: R(mod(-1, 3^7))
            2 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + O(3^5)
            sage: R(mod(-1, 3^2))
            2 + 2*3 + O(3^2)
            sage: R(3 + O(3^2))
            3 + O(3^2)

        Test that #3865 is fixed::

            sage: ZpCA(7, 10)(gp('7 + O(7^2)'))
            7 + O(7^2)
        """
        mpz_init(self.value)
        pAdicBaseGenericElement.__init__(self,parent)
        if empty:
            return
        if absprec is infinity or absprec > parent.precision_cap():
            absprec = parent.precision_cap()
        elif not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        if relprec is not infinity and not PY_TYPE_CHECK(relprec, Integer):
            relprec = Integer(relprec)
        if PY_TYPE_CHECK(x, pAdicGenericElement):
            if x.valuation() < 0:
                raise ValueError, "element valuation cannot be negative."
            if parent.prime() != x.parent().prime():
                raise TypeError, "Cannot coerce between p-adic parents with different primes"
        cdef mpz_t modulus
        cdef Integer tmp
        cdef unsigned long k
        if PY_TYPE_CHECK(x, pAdicBaseGenericElement):
            self._set_from_Integer(x._integer_(), min(absprec, x.precision_absolute()), relprec)
            return
        elif isinstance(x, (int, long)):
            x = Integer(x)
        elif isinstance(x, gp_element) or isinstance(x, pari_gen):
            if isinstance(x, gp_element): x = x._pari_()
            if x.type() == "t_PADIC":
                if x.variable() != self.prime_pow.prime:
                    raise TypeError, "Cannot coerce a pari p-adic with the wrong prime."
                absprec = min(Integer(x.padicprec(parent.prime())), absprec)
                x = x.lift()
            if x.type() == "t_INT":
                x = Integer(x)
            elif x.type() == "t_FRAC":
                x = Rational(x)
            else:
                raise TypeError, "unsupported coercion from pari: only p-adics, integers and rationals allowed"

        elif sage.rings.finite_rings.integer_mod.is_IntegerMod(x):
            mpz_init_set(modulus, (<Integer>x.modulus()).value)
            k = mpz_remove(modulus, modulus, self.prime_pow.prime.value)
            if mpz_cmp_ui(modulus, 1) == 0:
                if mpz_cmp_ui((<Integer>absprec).value, k) > 0:
                    absprec = PY_NEW(Integer)
                    mpz_set_ui((<Integer>absprec).value, k)
                x = x.lift()
                mpz_clear(modulus)
            else:
                mpz_clear(modulus)
                raise TypeError, "cannot coerce from the given integer mod ring (not a power of the same prime)"
        if PY_TYPE_CHECK(x, Integer):
            self._set_from_Integer(x, absprec, relprec)
            return
        if PY_TYPE_CHECK(x, Rational):
            self._set_from_Rational(x, absprec, relprec)
            return
        self._set_from_Rational(Rational(x), absprec, relprec)

    def __dealloc__(self):
        """
        Deallocation.

        TESTS::

            sage: R = ZpCA(5)
            sage: a = R(17)
            sage: del(a)
        """
        mpz_clear(self.value)

    def __reduce__(self):
        """
        Pickling.

        TESTS::

            sage: a = ZpCA(5)(-3)
            sage: type(a)
            <type 'sage.rings.padics.padic_capped_absolute_element.pAdicCappedAbsoluteElement'>
            sage: loads(dumps(a)) == a
            True
        """
        return make_pAdicCappedAbsoluteElement, (self.parent(), self.lift(), self.absprec)

    cdef bint _set_prec_abs(self, long absprec) except -1:
        """
        Sets ``self.absprec``.

        TESTS::

            sage: R = ZpCA(5)
            sage: a = R(17,5); a #indirect doctest
            2 + 3*5 + O(5^5)
        """
        self.absprec = absprec

    cdef bint _set_prec_both(self, long absprec, long relprec) except -1:
        """
        Sets ``self.absprec``.

        TESTS::

            sage: R = ZpCA(5)
            sage: a = R(75, absprec = 5, relprec = 4); a #indirect doctest
            3*5^2 + O(5^5)
        """
        self.absprec = absprec
        cdef long ordp
        if relprec < absprec:
            ordp = self.valuation_c()
            if ordp + relprec < absprec:
                self.absprec = ordp + relprec

    cdef int _set_from_Integer(pAdicCappedAbsoluteElement self, Integer x, absprec, relprec) except -1:
        """
        Set ``self.value`` and ``self.absprec``.

        ``absprec`` must be positive.

        TESTS::

            sage: R = ZpCA(5)
            sage: a = R(75, absprec = 5, relprec = 4); a #indirect doctest
            3*5^2 + O(5^5)
        """
        if relprec is not infinity and mpz_sgn((<Integer>relprec).value) == -1:
            raise ValueError, "relprec must be positive"
        if absprec is not infinity and mpz_sgn((<Integer>absprec).value) == -1:
            raise ValueError, "absprec must be positive"
        if relprec is infinity or mpz_fits_ulong_p((<Integer>relprec).value) == 0:
            if absprec is infinity or mpz_cmp_ui((<Integer>absprec).value, self.prime_pow.prec_cap) >= 0:
                return self._set_from_mpz_abs(x.value, self.prime_pow.prec_cap)
            else:
                return self._set_from_mpz_abs(x.value, mpz_get_si((<Integer>absprec).value))
        else:
            if absprec is infinity or mpz_cmp_ui((<Integer>absprec).value, self.prime_pow.prec_cap) >= 0:
                return self._set_from_mpz_both(x.value, self.prime_pow.prec_cap, mpz_get_si((<Integer>relprec).value))
            else:
                return self._set_from_mpz_both(x.value, mpz_get_si((<Integer>absprec).value), mpz_get_si((<Integer>relprec).value))

    cdef int _set_from_Rational(pAdicCappedAbsoluteElement self, Rational x, absprec, relprec) except -1:
        """
        Set ``self.value`` and ``self.absprec``.
        ``absprec`` must be positive.

        TESTS::

            sage: R = ZpCA(5)
            sage: a = R(25/9, absprec = 5, relprec = 4); a #indirect doctest
            4*5^2 + 2*5^3 + O(5^5)
        """
        if relprec is not infinity and mpz_sgn((<Integer>relprec).value) == -1:
            raise ValueError, "relprec must be positive"
        if absprec is not infinity and mpz_sgn((<Integer>absprec).value) == -1:
            raise ValueError, "absprec must be positive"
        if relprec is infinity or mpz_fits_ulong_p((<Integer>relprec).value) == 0:
            if absprec is infinity or mpz_cmp_ui((<Integer>absprec).value, self.prime_pow.prec_cap) >= 0:
                return self._set_from_mpq_abs(x.value, self.prime_pow.prec_cap)
            else:
                return self._set_from_mpq_abs(x.value, mpz_get_si((<Integer>absprec).value))
        else:
            if absprec is infinity or mpz_cmp_ui((<Integer>absprec).value, self.prime_pow.prec_cap) >= 0:
                return self._set_from_mpq_both(x.value, self.prime_pow.prec_cap, mpz_get_si((<Integer>relprec).value))
            else:
                return self._set_from_mpq_both(x.value, mpz_get_si((<Integer>absprec).value), mpz_get_si((<Integer>relprec).value))

    cdef int _set_from_mpz_abs(self, mpz_t x, long absprec) except -1:
        """
        ``self.prime_pow`` must already be set.

        TESTS::

            sage: R = ZpCA(5)
            sage: a = R(17,5); a #indirect doctest
            2 + 3*5 + O(5^5)
        """
        self._set_prec_abs(absprec)
        if mpz_sgn(x) == -1 or mpz_cmp(x, self.prime_pow.pow_mpz_t_tmp(self.absprec)[0]) >= 0:
            sig_on()
            mpz_mod(self.value, x, self.prime_pow.pow_mpz_t_tmp(self.absprec)[0])
            sig_off()
        else:
            mpz_set(self.value, x)
        return 0

    cdef int _set_from_mpz_both(self, mpz_t x, long absprec, long relprec) except -1:
        """
        ``self.prime_pow`` must already be set

        TESTS::

            sage: R = ZpCA(5)
            sage: a = R(75, absprec = 5, relprec = 4); a #indirect doctest
            3*5^2 + O(5^5)
        """
        if mpz_sgn(x) == 0:
            mpz_set(self.value, x)
            return 0
        mpz_set(self.value, x)
        self._set_prec_both(absprec, relprec)
        if mpz_sgn(x) == -1 or mpz_cmp(x, self.prime_pow.pow_mpz_t_tmp(self.absprec)[0]) >= 0:
            sig_on()
            mpz_mod(self.value, x, self.prime_pow.pow_mpz_t_tmp(self.absprec)[0])
            sig_off()
        return 0

    cdef int _set_from_mpq_abs(self, mpq_t x, long absprec) except -1:
        """
        ``self.prime_pow`` must already be set.

        TESTS::

            sage: R = ZpCA(5)
            sage: a = R(25/9, absprec = 5); a #indirect doctest
            4*5^2 + 2*5^3 + O(5^5)
        """
        self._set_prec_abs(absprec)
        if mpz_divisible_p(mpq_denref(x), self.prime_pow.prime.value):
            raise ValueError, "p divides denominator"
        sig_on()
        mpz_invert(self.value, mpq_denref(x), self.prime_pow.pow_mpz_t_tmp(absprec)[0])
        mpz_mul(self.value, self.value, mpq_numref(x))
        mpz_mod(self.value, self.value, self.prime_pow.pow_mpz_t_tmp(absprec)[0])
        sig_off()
        return 0

    cdef int _set_from_mpq_both(self, mpq_t x, long absprec, long relprec) except -1:
        """
        ``self.prime_pow`` must already be set

        TESTS::

            sage: R = ZpCA(5)
            sage: a = R(25/9, absprec = 5, relprec = 4); a #indirect doctest
            4*5^2 + 2*5^3 + O(5^5)
        """
        cdef long k
        cdef mpz_t tmp
        if mpq_sgn(x) == 0:
            mpz_set_ui(self.value, 0)
            return 0
        mpz_init(tmp)
        sig_on()
        k = mpz_remove(tmp, mpq_numref(x), self.prime_pow.prime.value)
        sig_off()
        mpz_clear(tmp)
        self.absprec = k + relprec
        if self.absprec > absprec:
            self.absprec = absprec
        return self._set_from_mpq_abs(x, self.absprec)

    cdef int _set_mpz_into(pAdicCappedAbsoluteElement self, mpz_t dest) except -1:
        """
        Sets ``dest`` to a lift of ``self``.

        TESTS::

            sage: R = ZpCA(5); S.<a> = ZqCA(25)
            sage: S(R(17))
            2 + 3*5 + O(5^20)
        """
        mpz_set(dest, self.value)
        return 0

    cdef int _set_mpq_into(pAdicCappedAbsoluteElement self, mpq_t dest) except -1:
        """
        Sets ``dest`` to a lift of ``self``.

        Not currently used internally.
        """
        mpq_set_z(dest, self.value)
        return 0

    cdef pAdicCappedAbsoluteElement _new_c(self):
        """
        Creates a new element with the same basic info.

        TESTS::

            sage: R = ZpCA(5); R(6,5) * R(7,8) #indirect doctest
            2 + 3*5 + 5^2 + O(5^5)
        """
        cdef pAdicCappedAbsoluteElement x
        x = PY_NEW(pAdicCappedAbsoluteElement)
        x._parent = self._parent
        x.prime_pow = self.prime_pow
        mpz_init(x.value)
        return x

    cpdef bint _is_inexact_zero(self) except -1:
        """
        Returns ``True`` if ``self`` is indistinguishable from zero.

        EXAMPLES::

            sage: R = ZpCA(7, 5)
            sage: R(7^5)._is_inexact_zero()
            True
            sage: R(0,4)._is_inexact_zero()
            True
            sage: R(0)._is_inexact_zero()
            True
        """
        return mpz_sgn(self.value) == 0

    def __richcmp__(left, right, op):
        """
        Comparison.

        TESTS::

            sage: R = ZpCA(5)
            sage: a = R(17)
            sage: b = R(21)
            sage: a == b
            False
            sage: a < b
            True
        """
        return (<Element>left)._richcmp(right, op)

    def __invert__(self):
        """
        Returns the multiplicative inverse of ``self``.

        EXAMPLES::

            sage: R = ZpCA(17)
            sage: ~R(-1) == R(-1)
            True
            sage: ~R(5) * 5
            1 + O(17^20)
            sage: ~R(5)
            7 + 3*17 + 10*17^2 + 13*17^3 + 6*17^4 + 3*17^5 + 10*17^6 + 13*17^7 + 6*17^8 + 3*17^9 + 10*17^10 + 13*17^11 + 6*17^12 + 3*17^13 + 10*17^14 + 13*17^15 + 6*17^16 + 3*17^17 + 10*17^18 + 13*17^19 + O(17^20)
        """
        return self._invert_c_impl()

    cpdef RingElement _invert_c_impl(self):
        """
        Returns the multiplicative inverse of ``self``.

        EXAMPLES::

            sage: R = ZpCA(17)
            sage: ~R(-1) == R(-1) #indirect doctest
            True
        """
        return self.parent().fraction_field()(self).__invert__()

    cpdef ModuleElement _neg_(self):
        """
        Returns the negation of self.

        EXAMPLES::

            sage: R = Zp(5, prec=10, type='capped-abs')
            sage: a = R(1)
            sage: -a #indirect doctest
            4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + O(5^10)
        """
        if mpz_sgn(self.value) == 0:
            return self
        cdef pAdicCappedAbsoluteElement ans
        ans = self._new_c()
        ans.absprec = self.absprec
        mpz_sub(ans.value, self.prime_pow.pow_mpz_t_tmp(self.absprec)[0], self.value)
        return ans

    def __pow__(pAdicCappedAbsoluteElement self, right, dummy):
        """
        Exponentiation.

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

        TESTS:

        We define ``0^0`` to be unity, :trac:`13941`::

            sage: R = ZpCA(11, 5)
            sage: R(0)^0
            1 + O(11^5)
            sage: R(0)^0 == R(1)
            True

        The value returned from ``0^0`` should belong to our ring::

            sage: R = ZpCA(11, 5)
            sage: type(R(0)^0) == type(R(1))
            True

        """
        # p-adic exponents and correct precisions!!!
        if (self == 0) and (right == 0):
            return self.parent(1)
        cdef Integer new, absprec
        cdef long val
        new = Integer(right) #Need to make sure that this works for p-adic exponents
        val = self.valuation_c()
        if (val > 0) and isinstance(right, pAdicBaseGenericElement):
            raise ValueError, "Can only have p-adic exponent if base is a unit"
        if (new < 0):
            return (~self).__pow__(-new)
        cdef pAdicCappedAbsoluteElement ans
        cdef Integer preccap
        preccap = <Integer>self.parent().precision_cap()
        ans = self._new_c()
        if val > 0:
            if new * val >= preccap:
                ans._set_from_Integer(Integer(0), preccap, infinity)
            else:
                absprec = <Integer> min(self.precision_relative() + new * val, preccap)
                ans._set_prec_abs(mpz_get_ui(absprec.value))
                sig_on()
                mpz_powm(ans.value, self.value, new.value, self.prime_pow.pow_mpz_t_tmp(ans.absprec)[0])
                sig_off()
        else:
            ans.absprec = self.absprec
            sig_on()
            mpz_powm(ans.value, self.value, new.value, self.prime_pow.pow_mpz_t_tmp(self.absprec)[0])
            sig_off()
        return ans

    cpdef ModuleElement _add_(self, ModuleElement _right):
        """
        Addition.

        EXAMPLES::

            sage: R = ZpCA(13, 4)
            sage: R(2) + R(3) #indirect doctest
            5 + O(13^4)
            sage: R(12) + R(1)
            13 + O(13^4)
        """
        cdef pAdicCappedAbsoluteElement ans
        cdef pAdicCappedAbsoluteElement right = <pAdicCappedAbsoluteElement> _right
        ans = self._new_c()
        if self.absprec < right.absprec:
            ans.absprec = self.absprec
        else:
            ans.absprec = right.absprec
        mpz_add(ans.value, self.value, right.value)
        if mpz_cmp(ans.value, self.prime_pow.pow_mpz_t_tmp(ans.absprec)[0]) >= 0:
            mpz_mod(ans.value, ans.value, self.prime_pow.pow_mpz_t_tmp(ans.absprec)[0])
        return ans

    cpdef RingElement _div_(self, RingElement right):
        """
        Division.

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
        return self * (~right)

    def __lshift__(pAdicCappedAbsoluteElement self, shift):
        """
        Multiplies ``self`` by ``p^shift``.

        If ``shift < -self.ordp()``, digits will be truncated.  See
        ``__rshift__`` for details.

        EXAMPLES:

        We create a capped absolute ring::

            sage: R = Zp(5, 20, 'capped-abs'); a = R(1000); a
            3*5^3 + 5^4 + O(5^20)

        Shifting to the right is the same as dividing by a power of
        the uniformizer `p` of the `p`-adic ring.::

            sage: a >> 1
            3*5^2 + 5^3 + O(5^19)

        Shifting to the left is the same as multiplying by a power of
        `p`::

            sage: a << 2
            3*5^5 + 5^6 + O(5^20)
            sage: a*5^2
            3*5^5 + 5^6 + O(5^20)

        Shifting by a negative integer to the left is the same as
        right shifting by the absolute value::

            sage: a << -3
            3 + 5 + O(5^17)
            sage: a >> 3
            3 + 5 + O(5^17)
        """
        cdef pAdicCappedAbsoluteElement ans
        if not PY_TYPE_CHECK(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_slong_p((<Integer>shift).value) == 0:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            if mpz_sgn((<Integer>shift).value) == 1:
                ans._set_prec_abs(self.prime_pow.prec_cap)
            else:
                ans._set_prec_abs(0)
            return ans
        return self._lshift_c(mpz_get_si((<Integer>shift).value))

    cdef pAdicCappedAbsoluteElement _lshift_c(pAdicCappedAbsoluteElement self, long shift):
        """
        Multiplies ``self`` by ``p^shift``.

        If ``shift < -self.ordp()``, digits will be truncated.  See
        ``__rshift__`` for details.

        EXAMPLES::

            sage: R = ZpCA(5); a = R(17); a << 2
            2*5^2 + 3*5^3 + O(5^20)
        """
        cdef unsigned long prec_cap, ansprec
        cdef pAdicCappedAbsoluteElement ans
        cdef PowComputer_class powerer
        if shift < 0:
            return self._rshift_c(-shift)
        prec_cap = self.prime_pow.prec_cap
        if shift >= prec_cap: # if we ever start caching valuation, this check should change
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            ans._set_prec_abs(prec_cap)
            return ans
        elif shift > 0:
            ans = self._new_c()
            ansprec = shift + self.absprec
            if ansprec > prec_cap:
                ansprec = prec_cap
            ans._set_prec_abs(ansprec)
            mpz_mul(ans.value, self.value, self.prime_pow.pow_mpz_t_tmp(shift)[0])
            if mpz_cmp(ans.value, ans.prime_pow.pow_mpz_t_top()[0]) >= 0:
                mpz_mod(ans.value, ans.value, ans.prime_pow.pow_mpz_t_top()[0])
            return ans
        else:
            return self

    def __rshift__(pAdicCappedAbsoluteElement self, shift):
        """
        Divides by ``p^shift``, and truncates.

        EXAMPLES::

            sage: R = Zp(997, 7, 'capped-abs'); a = R(123456878908); a
            964*997 + 572*997^2 + 124*997^3 + O(997^7)

        Shifting to the right divides by a power of `p`, but dropping
        terms with negative valuation::

            sage: a >> 3
            124 + O(997^4)

        A negative shift multiplies by that power of `p`.::

            sage: a >> -3
            964*997^4 + 572*997^5 + 124*997^6 + O(997^7)
        """
        cdef pAdicCappedAbsoluteElement ans
        if not PY_TYPE_CHECK(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_slong_p((<Integer>shift).value) == 0:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            if mpz_sgn((<Integer>shift).value) == -1:
                ans._set_prec_abs(self.prime_pow.prec_cap)
            else:
                ans._set_prec_abs(0)
            return ans
        return self._rshift_c(mpz_get_si((<Integer>shift).value))

    cdef pAdicCappedAbsoluteElement _rshift_c(pAdicCappedAbsoluteElement self, long shift):
        """
        Divides by ``p^shift`` and truncates.

        EXAMPLES::

            sage: R = ZpCA(5); a = R(77); a >> 1
            3*5 + O(5^19)
        """
        cdef unsigned long prec_cap, ansprec
        cdef pAdicCappedAbsoluteElement ans
        if shift < 0:
            return self._lshift_c(-shift)
        prec_cap = self.prime_pow.prec_cap
        if shift >= self.absprec:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            ans._set_prec_abs(0)
            return ans
        elif shift > 0:
            ans = self._new_c()
            ansprec = self.absprec - shift
            ans._set_prec_abs(ansprec)
            mpz_fdiv_q(ans.value, self.value, self.prime_pow.pow_mpz_t_tmp(shift)[0])
            return ans
        else:
            return self

    cpdef RingElement _mul_(self, RingElement _right):
        """
        Multiplication.

        EXAMPLES::

            sage: R = ZpCA(5)
            sage: a = R(20,5); b = R(75, 4); a * b #indirect doctest
            2*5^3 + 2*5^4 + O(5^5)
        """
        cdef pAdicCappedAbsoluteElement ans
        cdef pAdicCappedAbsoluteElement right = <pAdicCappedAbsoluteElement> _right
        cdef unsigned long sval, rval, prec1, prec2, prec_cap
        ans = self._new_c()
        prec_cap = self.prime_pow.prec_cap
        sval = self.valuation_c()
        rval = right.valuation_c()
        prec1 = sval + right.absprec
        prec2 = rval + self.absprec
        if prec1 < prec2:
            if prec_cap < prec1:
                ans._set_prec_abs(prec_cap)
            else:
                ans._set_prec_abs(prec1)
        else:
            if prec_cap < prec2:
                ans._set_prec_abs(prec_cap)
            else:
                ans._set_prec_abs(prec2)
        mpz_mul(ans.value, self.value, right.value)
        mpz_mod(ans.value, ans.value, self.prime_pow.pow_mpz_t_tmp(ans.absprec)[0])
        return ans

    cpdef ModuleElement _sub_(self, ModuleElement _right):
        """
        Subtraction.

        EXAMPLES::

            sage: R = ZpCA(13, 4)
            sage: R(10) - R(10) #indirect doctest
            O(13^4)
            sage: R(10) - R(11)
            12 + 12*13 + 12*13^2 + 12*13^3 + O(13^4)
        """
        cdef pAdicCappedAbsoluteElement ans
        cdef pAdicCappedAbsoluteElement right = <pAdicCappedAbsoluteElement> _right
        ans = self._new_c()
        if self.absprec < right.absprec:
            ans.absprec = self.absprec
        else:
            ans.absprec = right.absprec
        mpz_sub(ans.value, self.value, right.value)
        if mpz_sgn(ans.value) == -1 or mpz_cmp(ans.value, self.prime_pow.pow_mpz_t_tmp(ans.absprec)[0]) >= 0:
            mpz_mod(ans.value, ans.value, self.prime_pow.pow_mpz_t_tmp(ans.absprec)[0])
        return ans

    def add_bigoh(pAdicCappedAbsoluteElement self, absprec):
        """
        Returns a new element with absolute precision decreased to
        ``prec``.  The precision never increases.

        INPUT:

        - ``self`` -- a `p`-adic element

        - ``prec`` -- an integer

        OUTPUT:

        - ``element`` -- ``self`` with precision set to the minimum of ``self's`` precision and ``prec``

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
        cdef pAdicCappedAbsoluteElement ans
        cdef unsigned long newprec
        cdef Integer _absprec
        if not PY_TYPE_CHECK(absprec, Integer):
            _absprec = Integer(absprec)
        else:
            _absprec = <Integer>absprec
        if mpz_fits_ulong_p(_absprec.value) == 0:
            if mpz_sgn((<Integer>absprec).value) == -1:
                ans = self._new_c()
                ans._set_prec_abs(0)
                mpz_set_ui(ans.value, 0)
                return ans
            else:
                return self
        newprec = mpz_get_ui(_absprec.value)
        if newprec >= self.absprec:
            return self
        elif newprec <= 0:
            ans = self._new_c()
            ans._set_prec_abs(0)
            mpz_set_ui(ans.value, 0)
            return ans
        else:
            ans = self._new_c()
            ans._set_prec_abs(newprec)
            mpz_set(ans.value, self.value)
            if mpz_cmp(ans.value, self.prime_pow.pow_mpz_t_tmp(ans.absprec)[0]) >= 0:
                mpz_mod(ans.value, ans.value, self.prime_pow.pow_mpz_t_tmp(ans.absprec)[0])
            return ans

    def __copy__(pAdicCappedAbsoluteElement self):
        """
        Returns a copy of ``self``.

        EXAMPLES::

            sage: a = ZpCA(5,6)(17); b = copy(a)
            sage: a == b
            True
            sage: a is b
            False
        """
        cdef pAdicCappedAbsoluteElement ans
        ans = self._new_c()
        ans.absprec = self.absprec
        mpz_set(ans.value, self.value)
        return ans

    def is_zero(self, absprec = None):
        r"""
        Returns whether ``self`` is zero modulo ``p^absprec``.

        INPUT:

        - ``self`` -- a `p`-adic element

        - ``prec`` -- an integer

        OUTPUT:

        - ``boolean`` -- whether self is zero

        EXAMPLES::

            sage: R = ZpCA(17, 6)
            sage: R(0).is_zero()
            True
            sage: R(17^6).is_zero()
            True
            sage: R(17^2).is_zero(absprec=2)
            True
        """
        if absprec is None:
            return mpz_sgn(self.value) == 0
        cdef Integer _absprec
        if not PY_TYPE_CHECK(absprec, Integer):
            _absprec = Integer(absprec)
        else:
            _absprec = <Integer>absprec
        if mpz_cmp_ui(_absprec.value, self.absprec) > 0:
            raise PrecisionError, "Not enough precision to determine if element is zero"
        cdef mpz_t tmp
        mpz_init(tmp)
        cdef unsigned long aprec
        aprec = mpz_get_ui(_absprec.value)
        mpz_mod(tmp, self.value, self.prime_pow.pow_mpz_t_tmp(aprec)[0])
        if mpz_sgn(tmp) == 0:
            mpz_clear(tmp)
            return True
        else:
            mpz_clear(tmp)
            return False

    def is_equal_to(pAdicCappedAbsoluteElement self, pAdicCappedAbsoluteElement right, absprec = None):
        r"""
        Returns whether ``self`` is equal to ``right`` modulo ``p^absprec``.

        INPUT:

        - ``self`` -- a `p`-adic element

        - ``right`` -- a `p`-adic element with the same parent

        - ``absprec`` -- an integer

        OUTPUT:

        - ``boolean`` -- whether ``self`` is equal to ``right``

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
        """
        cdef unsigned long aprec
        if absprec is None:
            if self.absprec < right.absprec:
                aprec = self.absprec
            else:
                aprec = right.absprec
        else:
            if not PY_TYPE_CHECK(absprec, Integer):
                absprec = Integer(absprec)
            if mpz_fits_slong_p((<Integer>absprec).value) == 0:
                if mpz_sgn((<Integer>absprec).value) == 1:
                    raise PrecisionError, "Elements are not known to specified precision"
                else:
                    return True
            aprec = mpz_get_ui((<Integer>absprec).value)
            if aprec > self.absprec or aprec > right.absprec:
                raise PrecisionError, "Elements are not known to specified precision"
        cdef mpz_t tmp1, tmp2
        mpz_init(tmp1)
        mpz_init(tmp2)
        mpz_mod(tmp1, self.value, self.prime_pow.pow_mpz_t_tmp(aprec)[0])
        mpz_mod(tmp2, (right).value, self.prime_pow.pow_mpz_t_tmp(aprec)[0])
        if mpz_cmp(tmp1, tmp2) == 0:
            mpz_clear(tmp1)
            mpz_clear(tmp2)
            return True
        else:
            mpz_clear(tmp1)
            mpz_clear(tmp2)
            return False

    cpdef Integer lift(self):
        """
        Returns an integer congruent to this `p`-adic element modulo
        ``p^self.absprec()``.

        EXAMPLES::

            sage: R = ZpCA(3)
            sage: R(10).lift()
            10
            sage: R(-1).lift()
            3486784400
        """
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set(ans.value, self.value)
        return ans

    def lift_to_precision(self, absprec = None):
        """
        Returns another element of the same parent, with absolute
        precision at least absprec, congruent to this one modulo the
        known precision.

        INPUT:

        - ``absprec`` -- (default ``None``) the absolute precision of
          the result.  If ``None``, lifts to the maximum precision
          allowed.

        .. NOTE::

            If setting ``absprec`` that high would violate the
            precision cap, raises a precision error.

        EXAMPLES::

            sage: R = ZpCA(17)
            sage: R(-1,2).lift_to_precision(10)
            16 + 16*17 + O(17^10)
            sage: R(1,15).lift_to_precision(10)
            1 + O(17^15)
            sage: R(1,15).lift_to_precision(30)
            Traceback (most recent call last):
            ...
            PrecisionError: Precision higher than allowed by the precision cap.
            sage: R(-1,2).lift_to_precision().precision_absolute() == R.precision_cap()
            True
        """
        cdef pAdicCappedAbsoluteElement ans
        cdef unsigned long prec_cap
        cdef unsigned long dest_prec
        prec_cap = self.prime_pow.prec_cap
        if absprec is not None and not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        if absprec is None:
            dest_prec = prec_cap
        elif mpz_fits_slong_p((<Integer>absprec).value) == 0:
            if mpz_sgn((<Integer>absprec).value) < 0:
                return self
            else:
                raise PrecisionError("Precision higher than allowed by the precision cap.")
        else:
            dest_prec = mpz_get_ui((<Integer>absprec).value)
            if prec_cap < dest_prec:
                raise PrecisionError("Precision higher than allowed by the precision cap.")
        if dest_prec <= self.absprec:
            return self
        else:
            ans = self._new_c()
            ans._set_prec_abs(dest_prec)
            mpz_set(ans.value, self.value)
        return ans

    def list(pAdicCappedAbsoluteElement self, lift_mode = 'simple'):
        """
        Returns a list of coefficients of `p` starting with `p^0`

        INPUT:

        - ``self`` -- a `p`-adic element

        - ``lift_mode`` -- ``'simple'``, ``'smallest'`` or
          ``'teichmuller'`` (default ``'simple'``)

        OUTPUT:

        - ``list`` -- the list of coefficients of ``self``

        NOTES:

        - Returns a list `[a_0, a_1, \ldots, a_n]` so that:

          + If ``lift_mode = 'simple'``, `a_i` is an integer with `0
            \le a_i < p`.

          + If ``lift_mode = 'smallest'``, `a_i` is an integer with
            `-p/2 < a_i \le p/2`.

          + If ``lift_mode = 'teichmuller'``, `a_i` has the same
            parent as `self` and `a_i^p \equiv a_i` modulo
            ``p^(self.precision_absolute() - i)``

        - `\sum_{i = 0}^n a_i \cdot p^i =` ``self``, modulo the
          precision of ``self``.


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
        """
        if lift_mode == 'teichmuller':
            return self.teichmuller_list()
        if lift_mode == 'simple':
            return (<pAdicPrinter_class>self.parent()._printer).base_p_list(self.value, True)
        elif lift_mode == 'smallest':
            return (<pAdicPrinter_class>self.parent()._printer).base_p_list(self.value, False)
        else:
            raise ValueError

    cdef object teichmuller_list(pAdicCappedAbsoluteElement self):
        r"""
        Returns a list `[a_0, a_1,\ldots, a_n]` such that

        - `a_i^p = a_i`

        - ``self`` equals `\sum_{i = 0}^n a_i p^i`

        - if `a_i \ne 0`, the absolute precision of `a_i` is
          ``self.precision_relative() - i``

        EXAMPLES::

            sage: R = ZpFM(5,5); R(14).list('teichmuller') #indirect doctest
            [4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5),
            3 + 3*5 + 2*5^2 + 3*5^3 + 5^4 + O(5^5),
            2 + 5 + 2*5^2 + 5^3 + 3*5^4 + O(5^5),
            1 + O(5^5),
            4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5)]
        """
        # May eventually want to add a dict to store teichmuller lifts already seen, if p small enough
        cdef unsigned long curpower
        cdef mpz_t tmp, tmp2
        cdef pAdicCappedAbsoluteElement list_elt
        ans = PyList_New(0)
        curpower = self.absprec
        mpz_init_set(tmp, self.value)
        mpz_init(tmp2)
        while mpz_sgn(tmp) != 0:
            list_elt = self._new_c()
            list_elt._set_prec_abs(curpower)
            mpz_mod(list_elt.value, tmp, self.prime_pow.prime.value)
            mpz_set(tmp2, self.prime_pow.pow_mpz_t_tmp(curpower)[0])
            self.teichmuller_set_c(list_elt.value, tmp2)
            mpz_sub(tmp, tmp, list_elt.value)
            mpz_divexact(tmp, tmp, self.prime_pow.prime.value)
            curpower -= 1
            mpz_mod(tmp, tmp, self.prime_pow.pow_mpz_t_tmp(curpower)[0])
            PyList_Append(ans, list_elt)
        mpz_clear(tmp2)
        mpz_clear(tmp)
        return ans

    def _teichmuller_set(self):
        """
        Sets ``self`` to be the Teichmuller representative with the
        same residue as ``self``.

        WARNING:

        This function modifies ``self``, which is not safe.  Elements
        are supposed to be immutable.

        EXAMPLES::

            sage: R = ZpCA(17,5); a = R(11)
            sage: a
            11 + O(17^5)
            sage: a._teichmuller_set(); a
            11 + 14*17 + 2*17^2 + 12*17^3 + 15*17^4 + O(17^5)
            sage: a.list('teichmuller')
            [11 + 14*17 + 2*17^2 + 12*17^3 + 15*17^4 + O(17^5)]
        """
        cdef mpz_t tmp
        mpz_init_set(tmp, self.prime_pow.pow_mpz_t_tmp(self.absprec)[0])
        self.teichmuller_set_c(self.value, tmp)
        mpz_clear(tmp)

    def multiplicative_order(self):
        r"""
        Returns the minimum possible multiplicative order of ``self``.

        INPUT:

        - ``self`` -- a `p`-adic element

        OUTPUT:

        - the multiplicative order of self.  This is the minimum
          multiplicative order of all elements of `\mathbb{Z}_p`
          lifting ``self`` to infinite precision.

        EXAMPLES::

            sage: R = ZpCA(7, 6)
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
        cdef mpz_t ppow_minus_one
        cdef Integer ans
        if mpz_divisible_p(self.value, self.prime_pow.prime.value):
            return infinity
        if mpz_cmp_ui(self.value, 1) == 0:
            ans = PY_NEW(Integer)
            mpz_set_ui(ans.value, 1)
            return ans
        mpz_init(ppow_minus_one)
        mpz_sub_ui(ppow_minus_one, self.prime_pow.pow_mpz_t_tmp(self.absprec)[0], 1)
        if mpz_cmp(self.value, ppow_minus_one) == 0:
            ans = PY_NEW(Integer)
            mpz_set_ui(ans.value, 2)
            mpz_clear(ppow_minus_one)
            return ans
        # check if self is an approximation to a teichmuller lift:
        mpz_powm(ppow_minus_one, self.value, self.prime_pow.prime.value, self.prime_pow.pow_mpz_t_tmp(self.absprec)[0])
        if mpz_cmp(ppow_minus_one, self.value) == 0:
            mpz_clear(ppow_minus_one)
            return self.residue(1).multiplicative_order()
        else:
            mpz_clear(ppow_minus_one)
            return infinity

    def padded_list(self, n, list_mode = 'simple'):
        """
        Returns a list of coefficients of `p` starting with `p^0` up
        to `p^n` exclusive (padded with zeros if needed)

        INPUT:

        - ``self`` -- a `p`-adic element

        - ``n`` - an integer

        OUTPUT:

        - ``list`` -- the list of coefficients of ``self``

        EXAMPLES::

            sage: R = Zp(7,4,'capped-abs'); a = R(2*7+7**2); a.padded_list(5)
            [0, 2, 1, 0, 0]

        NOTE:

        this differs from the padded_list method of padic_field_element

        the slice operators throw an error if asked for a slice above
        the precision, while this function works
        """
        if list_mode == 'simple' or list_mode == 'smallest':
            zero = Integer(0)
        else:
            zero = self.parent()(0, 0)
        L = self.list()
        return L[:n] + [zero] * (n - len(L))

    def precision_absolute(self):
        """
        Returns the absolute precision of ``self``.

        This is the power of the maximal ideal modulo which this
        element is defined.

        INPUT:

        - ``self`` -- a `p`-adic element

        OUTPUT:

        - ``integer`` -- the absolute precision of ``self``

        EXAMPLES::

            sage: R = Zp(7,4,'capped-abs'); a = R(7); a.precision_absolute()
            4
       """
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set_si(ans.value, self.absprec)
        return ans

    def precision_relative(self):
        """
        Returns the relative precision of ``self``.

        This is the power of the maximal ideal modulo which the unit
        part of ``self`` is defined.

        INPUT:

        - ``self`` -- a `p`-adic element

        OUTPUT:

        - ``integer`` -- the relative precision of ``self``

        EXAMPLES::

            sage: R = Zp(7,4,'capped-abs'); a = R(7); a.precision_relative()
            3
       """
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set_si(ans.value, self.absprec - self.valuation_c())
        return ans

    def residue(self, absprec = 1):
        r"""
        Reduces ``self`` modulo ``p^absprec``

        INPUT:

        - ``self`` -- a `p`-adic element

        - ``absprec`` - an integer

        OUTPUT:

        - element of `\mathbb{Z}/p^{\mbox{absprec}} \mathbb{Z}` --
          ``self`` reduced modulo ``p^absprec``.

        EXAMPLES::

            sage: R = Zp(7,4,'capped-abs'); a = R(8); a.residue(1)
            1
        """
        if absprec > self.precision_absolute():
            raise PrecisionError, "Not enough precision known in order to compute residue."
        elif absprec < 0:
            raise ValueError, "cannot reduce modulo a negative power of p"
        cdef Integer selfvalue
        selfvalue = PY_NEW(Integer)
        mpz_set(selfvalue.value, self.value)
        return Mod(selfvalue, self.parent().prime_pow(absprec))


    #def square_root(self):
    #    r"""
    #    Returns the square root of this p-adic number

    #    INPUT:
    #        self -- a p-adic element
    #    OUTPUT:
    #        p-adic element -- the square root of this p-adic number
    #    EXAMPLES:
    #        sage: R = Zp(3,20,'capped-abs')
    #        sage: R(0).square_root()
    #            O(3^10)
    #        sage: R(1).square_root()
    #            1 + O(3^20)
    #        sage: R(4).square_root() == R(-2)
    #            True
    #        sage: R(9).square_root()
    #            3 + O(3^19)
    #        sage: R2 = Zp(2,20,'capped-abs')
    #        sage: R2(0).square_root()
    #            O(2^10)
    #        sage: R2(1).square_root()
    #            1 + O(2^19)
    #        sage: R2(4).square_root()
    #            2 + O(2^18)
    #        sage: R2(9).square_root() == R2(3) or R2(9).square_root() == R2(-3)
    #            True
    #        sage: R2(17).square_root()
    #            1 + 2^3 + 2^5 + 2^6 + 2^7 + 2^9 + 2^10 + 2^13 + 2^16 + 2^17 + O(2^19)
    #        sage: R3 = Zp(5,20,'capped-abs')
    #        sage: R3(0).square_root()
    #            O(5^10)
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
    #        raise ValueError, "element is not a square" # should eventually change to return an element of an extension field

    cpdef pAdicCappedAbsoluteElement unit_part(self):
        r"""
        Returns the unit part of ``self``.

        INPUT:

        - ``self`` -- a `p`-adic element

        OUTPUT:

        - `p`-adic element -- the unit part of ``self``

        EXAMPLES::

            sage: R = Zp(17,4,'capped-abs', 'val-unit')
            sage: a = R(18*17)
            sage: a.unit_part()
            18 + O(17^3)
            sage: type(a)
            <type 'sage.rings.padics.padic_capped_absolute_element.pAdicCappedAbsoluteElement'>
        """
        cdef pAdicCappedAbsoluteElement ans
        cdef unsigned long v
        if mpz_sgn(self.value) == 0:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            ans._set_prec_abs(0)
            return ans
        elif mpz_divisible_p(self.value, self.prime_pow.prime.value):
            ans = self._new_c()
            sig_on()
            v = mpz_remove(ans.value, self.value, self.prime_pow.prime.value)
            sig_off()
            ans._set_prec_abs(self.absprec - v)
            return ans
        else:
            return self

    def valuation(self):
        """
        Returns the valuation of ``self``, ie the largest power of `p`
        dividing ``self``.

        EXAMPLES::

            sage: R = ZpCA(5)
            sage: R(5^5*1827).valuation()
            5

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
            sage: R(0).valuation()
            20
        """
        # We override this, rather than using the valuation in
        # padic_generic_element, for speed reasons.

        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set_ui(ans.value, self.valuation_c())
        return ans

    cdef long valuation_c(self):
        """
        Returns the valuation of ``self``, ie the largest power of `p`
        dividing ``self``.

        EXAMPLES::

            sage: R = ZpCA(5)
            sage: a = R(0,6)
            sage: a.valuation() #indirect doctest
            6
        """
        if mpz_sgn(self.value) == 0:
            return self.absprec
        cdef mpz_t tmp
        cdef long ans
        mpz_init(tmp)
        ans = mpz_remove(tmp, self.value, self.prime_pow.prime.value)
        mpz_clear(tmp)
        return ans

    cpdef val_unit(self):
        """
        Returns a 2-tuple, the first element set to the valuation of
        ``self``, and the second to the unit part of ``self``.

        If ``self = 0``, then the unit part is ``O(p^0)``.

        EXAMPLES::

            sage: R = ZpCA(5)
            sage: a = R(75, 6); b = a - a
            sage: a.val_unit()
            (2, 3 + O(5^4))
            sage: b.val_unit()
            (6, O(5^0))
        """
        cdef pAdicCappedAbsoluteElement unit
        cdef Integer val
        cdef unsigned long v
        val = PY_NEW(Integer)
        if mpz_sgn(self.value) == 0:
            unit = self._new_c()
            mpz_set_ui(unit.value, 0)
            unit._set_prec_abs(0)
            mpz_set_ui(val.value, self.absprec)
            return (val, unit)
        elif mpz_divisible_p(self.value, self.prime_pow.prime.value):
            unit = self._new_c()
            v = mpz_remove(unit.value, self.value, self.prime_pow.prime.value)
            unit._set_prec_abs(self.absprec - v)
            mpz_set_ui(val.value, v)
            return (val, unit)
        else:
            mpz_set_ui(val.value, 0)
            return (val, self)

    def __hash__(self):
        """
        Hashing.

        EXAMPLES::

            sage: R = ZpCA(11, 5)
            sage: hash(R(3)) == hash(3)
            True
        """
        return hash(self.lift())

def make_pAdicCappedAbsoluteElement(parent, x, absprec):
    """
    Unpickles a capped absolute element.

    EXAMPLES::

        sage: from sage.rings.padics.padic_capped_absolute_element import make_pAdicCappedAbsoluteElement
        sage: R = ZpCA(5)
        sage: a = make_pAdicCappedAbsoluteElement(R, 17*25, 5); a
        2*5^2 + 3*5^3 + O(5^5)
    """
    return parent(x, absprec=absprec)
