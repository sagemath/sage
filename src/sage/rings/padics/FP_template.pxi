"""
Floating point template for complete discrete valuation rings

In order to use this template you need to write a linkage file and
gluing file.  For an example see ``mpz_linkage.pxi`` (linkage file) and
``padic_floating_point_element.pyx`` (gluing file).

The linkage file implements a common API that is then used in the
class :class:`FPElement` defined here.  See ``sage/libs/linkages/padics/API.pxi``
for the functions needed.

The gluing file does the following:

- ``ctypedef``'s ``celement`` to be the appropriate type (e.g. ``mpz_t``)
- includes the linkage file
- includes this template
- defines a concrete class inheriting from :class:`FPElement`, and implements
  any desired extra methods

AUTHORS:

- David Roe (2016-03-21) -- initial version
"""

#*****************************************************************************
#       Copyright (C) 2007-2016 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.ext.stdsage cimport PY_NEW
include "padic_template_element.pxi"
from cpython.long cimport *

from sage.structure.element cimport Element
from sage.rings.padics.common_conversion cimport comb_prec, _process_args_and_kwds
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.categories.sets_cat import Sets
from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
from sage.categories.homset import Hom

cdef inline bint overunderflow(long* ordp, celement unit, PowComputer_ prime_pow) noexcept:
    """
    Check for over and underflow.  If detected, sets ordp and unit
    appropriately, and returns True.  If not, returns False.
    """
    if ordp[0] >= maxordp:
        ordp[0] = maxordp
        csetzero(unit, prime_pow)
    elif ordp[0] <= minusmaxordp:
        ordp[0] = minusmaxordp
        csetone(unit, prime_pow)
    else:
        return False
    return True

cdef inline bint overunderflow_mpz(long* ordp, mpz_t ordp_mpz, celement unit, PowComputer_ prime_pow) noexcept:
    """
    Check for over and underflow with an mpz_t ordp.  If detected, sets ordp and unit
    appropriately, and returns True.  If not, returns False.
    """
    if mpz_fits_slong_p(ordp_mpz) == 0 or mpz_cmp_si(ordp_mpz, maxordp) >= 0 or mpz_cmp_si(ordp_mpz, minusmaxordp) <= 0:
        if mpz_sgn(ordp_mpz) > 0:
            ordp[0] = maxordp
            csetzero(unit, prime_pow)
        else:
            ordp[0] = minusmaxordp
            csetone(unit, prime_pow)
        return True
    return False

cdef inline bint very_pos_val(long ordp) noexcept:
    return ordp >= maxordp

cdef inline bint very_neg_val(long ordp) noexcept:
    return ordp <= minusmaxordp

cdef inline bint huge_val(long ordp) noexcept:
    return very_pos_val(ordp) or very_neg_val(ordp)

cdef class FPElement(pAdicTemplateElement):
    cdef int _set(self, x, long val, long xprec, absprec, relprec) except -1:
        """
        Set the value of this element from given defining data.

        This function is intended for use in conversion, and should
        not be called on an element created with :meth:`_new_c`.

        INPUT:

        - ``x`` -- data defining a `p`-adic element: int, long,
          Integer, Rational, other `p`-adic element...

        - ``val`` -- the valuation of the resulting element

        - ``xprec -- an inherent precision of ``x``, if ``val``
          is larger then the result will be zero

        - ``absprec`` -- an absolute precision cap for this element,
          if ``val`` is larger then the result will be zero

        - ``relprec`` -- a relative precision cap for this element
          (unused; for compatibility with other `p`-adic precision
          modes)

        TESTS::

            sage: R = ZpFP(5)
            sage: a = R(17,5); a  # indirect doctest
            2 + 3*5
            sage: R(15)  # indirect doctest
            3*5

            sage: R = ZpFP(5,5)
            sage: a = R(25/9); a  # indirect doctest
            4*5^2 + 2*5^3 + 5^5 + 2*5^6
            sage: R(ZpCR(5)(25/9)) == a
            True
            sage: R(5) - R(5)
            0

        We check that :issue:`23966` is resolved::

            sage: R = ZpFM(2)
            sage: K = R.fraction_field()
            sage: K(R.zero())
            0
        """
        IF CELEMENT_IS_PY_OBJECT:
            polyt = type(self.prime_pow.modulus)
            self.unit = <celement>polyt.__new__(polyt)
        cconstruct(self.unit, self.prime_pow)
        if val >= xprec or val >= absprec:
            self._set_exact_zero()
        elif very_neg_val(val):
            self._set_infinity()
        else:
            self.ordp = val
            if isinstance(x, FPElement) and x.parent() is self.parent():
                ccopy(self.unit, (<FPElement>x).unit, self.prime_pow)
            else:
                cconv(self.unit, x, self.prime_pow.ram_prec_cap, val, self.prime_pow)

    cdef int _set_exact_zero(self) except -1:
        """
        Set this element to zero.

        TESTS::

            sage: R = Zp(5); R(0)  # indirect doctest
            0
        """
        csetzero(self.unit, self.prime_pow)
        self.ordp = maxordp

    cdef int _set_infinity(self) except -1:
        """
        Set this element to zero.

        TESTS::

            sage: R = Zp(5); R(0)  # indirect doctest
            0
        """
        csetone(self.unit, self.prime_pow)
        self.ordp = minusmaxordp

    cdef FPElement _new_c(self):
        """
        Create a new element with the same basic info.

        TESTS::

            sage: R = ZpFP(5); R(6) * R(7)  # indirect doctest
            2 + 3*5 + 5^2

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFP(25)
            sage: S.<x> = ZZ[]
            sage: W.<w> = R.ext(x^2 - 5)
            sage: w * (w+1)  # indirect doctest
            w + w^2
        """
        cdef type t = type(self)
        cdef FPElement ans = t.__new__(t)
        ans._parent = self._parent
        ans.prime_pow = self.prime_pow
        IF CELEMENT_IS_PY_OBJECT:
            polyt = type(self.prime_pow.modulus)
            ans.unit = <celement>polyt.__new__(polyt)
        cconstruct(ans.unit, ans.prime_pow)
        return ans

    cdef pAdicTemplateElement _new_with_value(self, celement value, long absprec):
        """
        Create a new element with a given value and absolute precision.

        Used by code that doesn't know the precision type.
        """
        cdef FPElement ans = self._new_c()
        ans.ordp = 0
        ccopy(ans.unit, value, ans.prime_pow)
        ans._normalize()
        return ans

    cdef int _get_unit(self, celement value) except -1:
        """
        Set ``value`` to the unit of this `p`-adic element.
        """
        ccopy(value, self.unit, self.prime_pow)

    cdef int check_preccap(self) except -1:
        """
        Check that the precision of this element does not exceed the
        precision cap. Does nothing for floating point elements.

        TESTS::

            sage: ZpFP(5)(1).lift_to_precision(30)  # indirect doctest
            1
        """
        pass

    def __copy__(self):
        """
        Return a copy of this element.

        EXAMPLES::

            sage: a = ZpFP(5,6)(17); b = copy(a)
            sage: a == b
            True
            sage: a is b
            False
        """
        cdef FPElement ans = self._new_c()
        ans.ordp = self.ordp
        ccopy(ans.unit, self.unit, ans.prime_pow)
        return ans

    cdef int _normalize(self) except -1:
        """
        Normalize this element, so that ``self.ordp`` is correct.

        TESTS::

            sage: R = ZpFP(5)
            sage: R(6) + R(4)  # indirect doctest
            2*5
        """
        cdef long diff
        cdef bint is_zero
        if very_pos_val(self.ordp):
            self._set_exact_zero()
        elif very_neg_val(self.ordp):
            self._set_infinity()
        else:
            is_zero = creduce(self.unit, self.unit, self.prime_pow.ram_prec_cap, self.prime_pow)
            if is_zero:
                self.ordp = maxordp
            else:
                diff = cremove(self.unit, self.unit, self.prime_pow.ram_prec_cap, self.prime_pow)
                self.ordp += diff
                if very_pos_val(self.ordp):
                    self._set_exact_zero()

    def __dealloc__(self):
        """
        Deallocate the underlying data structure.

        TESTS::

            sage: R = ZpFP(5)
            sage: a = R(17)
            sage: del(a)
        """
        cdestruct(self.unit, self.prime_pow)

    def __reduce__(self):
        """
        Return a tuple of a function and data that can be used to unpickle this
        element.

        EXAMPLES::

            sage: a = ZpFP(5)(-3)
            sage: type(a)
            <class 'sage.rings.padics.padic_floating_point_element.pAdicFloatingPointElement'>
            sage: loads(dumps(a)) == a
            True
        """
        return unpickle_fpe_v2, (self.__class__, self.parent(), cpickle(self.unit, self.prime_pow), self.ordp)

#    def __richcmp__(self, right, int op):
#        """
#        Compare this element to ``right`` using the comparison operator ``op``.
#
#        TESTS::
#
#            sage: R = ZpFP(5)
#            sage: a = R(17)
#            sage: b = R(21)
#            sage: a == b
#            False
#            sage: a < b
#            True
#        """
#        return (<Element>self)._richcmp(right, op)

    cpdef _neg_(self):
        r"""
        Return the additive inverse of this element.

        EXAMPLES::

            sage: R = Zp(7, 4, 'floating-point', 'series')
            sage: -R(7)  # indirect doctest
            6*7 + 6*7^2 + 6*7^3 + 6*7^4
        """
        cdef FPElement ans = self._new_c()
        ans.ordp = self.ordp
        if huge_val(self.ordp): # zero or infinity
            ccopy(ans.unit, self.unit, ans.prime_pow)
        else:
            cneg(ans.unit, self.unit, ans.prime_pow.ram_prec_cap, ans.prime_pow)
            creduce_small(ans.unit, ans.unit, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        return ans

    cpdef _add_(self, _right):
        r"""
        Return the sum of this element and ``_right``.

        EXAMPLES::

            sage: R = Zp(7, 4, 'floating-point', 'series')
            sage: x = R(1721); x
            6 + 5*7^3
            sage: y = R(1373); y
            1 + 4*7^3
            sage: x + y  # indirect doctest
            7 + 2*7^3
        """
        cdef FPElement ans
        cdef FPElement right = _right
        cdef long tmpL
        if self.ordp == right.ordp:
            ans = self._new_c()
            ans.ordp = self.ordp
            if huge_val(ans.ordp):
                ccopy(ans.unit, self.unit, ans.prime_pow)
            else:
                cadd(ans.unit, self.unit, right.unit, ans.prime_pow.ram_prec_cap, ans.prime_pow)
                ans._normalize() # safer than trying to leave unnormalized
        else:
            if self.ordp > right.ordp:
                # Addition is commutative, swap so self.ordp < right.ordp
                ans = right
                right = self
                self = ans
            tmpL = right.ordp - self.ordp
            if tmpL > self.prime_pow.ram_prec_cap:
                return self
            ans = self._new_c()
            ans.ordp = self.ordp
            if huge_val(ans.ordp):
                ccopy(ans.unit, self.unit, ans.prime_pow)
            else:
                cshift_notrunc(ans.unit, right.unit, tmpL, ans.prime_pow.ram_prec_cap, ans.prime_pow, False)
                cadd(ans.unit, ans.unit, self.unit, ans.prime_pow.ram_prec_cap, ans.prime_pow)
                creduce(ans.unit, ans.unit, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        return ans

    cpdef _sub_(self, _right):
        r"""
        Return the difference of this element and ``_right``.

        EXAMPLES::

            sage: R = Zp(7, 4, 'floating-point', 'series')
            sage: x = R(1721); x
            6 + 5*7^3
            sage: y = R(1373); y
            1 + 4*7^3
            sage: x - y  # indirect doctest
            5 + 7^3
        """
        cdef FPElement ans
        cdef FPElement right = _right
        cdef long tmpL
        if self.ordp == right.ordp:
            ans = self._new_c()
            ans.ordp = self.ordp
            if huge_val(ans.ordp):
                ccopy(ans.unit, self.unit, ans.prime_pow)
            else:
                csub(ans.unit, self.unit, right.unit, ans.prime_pow.ram_prec_cap, ans.prime_pow)
                ans._normalize() # safer than trying to leave unnormalized
        elif self.ordp < right.ordp:
            tmpL = right.ordp - self.ordp
            if tmpL > self.prime_pow.ram_prec_cap:
                return self
            ans = self._new_c()
            ans.ordp = self.ordp
            if huge_val(ans.ordp):
                ccopy(ans.unit, self.unit, ans.prime_pow)
            else:
                cshift_notrunc(ans.unit, right.unit, tmpL, ans.prime_pow.ram_prec_cap, ans.prime_pow, False)
                csub(ans.unit, self.unit, ans.unit, ans.prime_pow.ram_prec_cap, ans.prime_pow)
                creduce(ans.unit, ans.unit, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        else:
            tmpL = self.ordp - right.ordp
            if tmpL > self.prime_pow.ram_prec_cap:
                return right._neg_()
            ans = self._new_c()
            ans.ordp = right.ordp
            if huge_val(ans.ordp):
                ccopy(ans.unit, self.unit, ans.prime_pow)
            else:
                cshift_notrunc(ans.unit, self.unit, tmpL, ans.prime_pow.ram_prec_cap, ans.prime_pow, False)
                csub(ans.unit, ans.unit, right.unit, ans.prime_pow.ram_prec_cap, ans.prime_pow)
                creduce(ans.unit, ans.unit, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        return ans

    def __invert__(self):
        r"""
        Return multiplicative inverse of this element.

        EXAMPLES::

            sage: R = Zp(7, 4, 'floating-point', 'series')
            sage: ~R(2)
            4 + 3*7 + 3*7^2 + 3*7^3
            sage: ~R(0)
            infinity
            sage: ~R(7)
            7^-1
        """
        # Input should be normalized!
        cdef FPElement ans = self._new_c()
        if ans.prime_pow.in_field == 0:
            ans._parent = self._parent.fraction_field()
            ans.prime_pow = ans._parent.prime_pow
        ans.ordp = -self.ordp
        if very_pos_val(ans.ordp):
            csetone(ans.unit, ans.prime_pow)
        elif very_neg_val(ans.ordp):
            csetzero(ans.unit, ans.prime_pow)
        else:
            cinvert(ans.unit, self.unit, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        return ans

    cpdef _mul_(self, _right):
        r"""
        Return the product of this element and ``_right``.

        EXAMPLES::

            sage: R = Zp(7, 4, 'floating-point', 'series')
            sage: R(3) * R(2)  # indirect doctest
            6
            sage: R(1/2) * R(2)
            1
        """
        cdef FPElement right = _right
        if very_pos_val(self.ordp):
            if very_neg_val(right.ordp):
                raise ZeroDivisionError("Cannot multiply 0 by infinity")
            return self
        elif very_pos_val(right.ordp):
            if very_neg_val(self.ordp):
                raise ZeroDivisionError("Cannot multiply 0 by infinity")
            return right
        elif very_neg_val(self.ordp):
            return self
        elif very_neg_val(right.ordp):
            return right
        cdef FPElement ans = self._new_c()
        ans.ordp = self.ordp + right.ordp
        if overunderflow(&ans.ordp, ans.unit, ans.prime_pow):
            return ans
        cmul(ans.unit, self.unit, right.unit, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        creduce(ans.unit, ans.unit, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        return ans

    cpdef _div_(self, _right):
        r"""
        Return the quotient of this element and ``right``.

        EXAMPLES::

            sage: R = Zp(7, 4, 'floating-point', 'series')
            sage: R(3) / R(2)  # indirect doctest
            5 + 3*7 + 3*7^2 + 3*7^3
            sage: R(5) / R(0)
            infinity
            sage: R(7) / R(49)
            7^-1
        """
        # Input should be normalized!
        cdef FPElement right = _right
        cdef FPElement ans = self._new_c()
        if ans.prime_pow.in_field == 0:
            ans._parent = self._parent.fraction_field()
            ans.prime_pow = ans._parent.prime_pow
        if very_pos_val(self.ordp):
            if very_pos_val(right.ordp):
                raise ZeroDivisionError("Cannot divide 0 by 0")
            ans._set_exact_zero()
        elif very_neg_val(right.ordp):
            if very_neg_val(self.ordp):
                raise ZeroDivisionError("Cannot divide infinity by infinity")
            ans._set_exact_zero()
        elif very_neg_val(self.ordp) or very_pos_val(right.ordp):
            ans._set_infinity()
        else:
            ans.ordp = self.ordp - right.ordp
            if overunderflow(&ans.ordp, ans.unit, ans.prime_pow):
                return ans
            cdivunit(ans.unit, self.unit, right.unit, ans.prime_pow.ram_prec_cap, ans.prime_pow)
            creduce(ans.unit, ans.unit, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        return ans

    def _quo_rem(self, _right):
        """
        Quotient with remainder.

        We choose the remainder to have the same `p`-adic expansion
        as the numerator, but truncated at the valuation of the denominator.

        EXAMPLES::

            sage: R = ZpFP(3, 5)
            sage: R(12).quo_rem(R(2))  # indirect doctest
            (2*3, 0)
            sage: R(2).quo_rem(R(12))
            (0, 2)
            sage: q, r = R(4).quo_rem(R(12)); q, r
            (1 + 2*3 + 2*3^3, 1)
            sage: 12*q + r == 4
            True

        For fields the normal quotient always has remainder 0:

            sage: K = QpFP(3, 5)
            sage: K(12).quo_rem(K(2))
            (2*3, 0)
            sage: q, r = K(4).quo_rem(K(12)); q, r
            (3^-1, 0)
            sage: 12*q + r == 4
            True

        You can get the same behavior for fields as for rings
        by using integral=True::

            sage: K(12).quo_rem(K(2), integral=True)
            (2*3, 0)
            sage: K(2).quo_rem(K(12), integral=True)
            (0, 2)
        """
        cdef FPElement right = _right
        if very_pos_val(right.ordp):
            raise ZeroDivisionError("Cannot find quo_rem by 0")
        elif very_neg_val(right.ordp):
            raise ZeroDivisionError("Cannot find quo_rem by infinity")
        if huge_val(self.ordp):
            return self, self
        cdef FPElement q = self._new_c()
        cdef FPElement r = self._new_c()
        cdef long diff = self.ordp - right.ordp
        if diff >= 0:
            q.ordp = diff
            cdivunit(q.unit, self.unit, right.unit, q.prime_pow.ram_prec_cap, q.prime_pow)
            r._set_exact_zero()
        else:
            r.ordp = self.ordp
            q.ordp = 0
            cshift(q.prime_pow.shift_rem, r.unit, self.unit, diff, q.prime_pow.ram_prec_cap, q.prime_pow, False)
            cdivunit(q.unit, q.prime_pow.shift_rem, right.unit, q.prime_pow.ram_prec_cap, q.prime_pow)
        q._normalize()
        return q, r

    def __pow__(FPElement self, _right, dummy): # NOTE: dummy ignored, always use self.prime_pow.ram_prec_cap
        """
        Exponentiation by an integer

        EXAMPLES::

            sage: R = ZpFP(11, 5)
            sage: R(1/2)^5
            10 + 7*11 + 11^2 + 5*11^3 + 4*11^4
            sage: R(1/32)
            10 + 7*11 + 11^2 + 5*11^3 + 4*11^4
            sage: R(1/2)^5 == R(1/32)
            True
            sage: R(3)^1000  # indirect doctest
            1 + 4*11^2 + 3*11^3 + 7*11^4
            sage: R(11)^-1
            11^-1

        TESTS:

        Check that :issue:`31875` is fixed::

            sage: R(1)^R(0)
            1

            sage: S.<a> = ZqFP(4)                                                       # needs sage.libs.flint
            sage: S(1)^S(0)                                                             # needs sage.libs.flint
            1
        """
        cdef long dummyL
        cdef mpz_t tmp
        cdef Integer right
        cdef FPElement base, pright, ans
        cdef bint exact_exp
        if isinstance(_right, (Integer, int, Rational)):
            if _right < 0:
                self = ~self
                _right = -_right
            exact_exp = True
        elif self.parent() is _right.parent():
            # For extension elements, we need to switch to the
            # fraction field sometimes in highly ramified extensions.
            exact_exp = (<FPElement>_right)._is_exact_zero()
            pright = _right
        else:
            self, _right = canonical_coercion(self, _right)
            return self.__pow__(_right, dummy)
        if exact_exp and _right == 0:
            ans = self._new_c()
            ans.ordp = 0
            csetone(ans.unit, ans.prime_pow)
            return ans
        if huge_val(self.ordp):
            if exact_exp:
                # We may assume from above that right > 0
                return self
            else:
                # log(0) and log(infinity) not defined
                raise ValueError("0^x and inf^x not defined for p-adic x")
        ans = self._new_c()
        if exact_exp:
            # exact_pow_helper is defined in padic_template_element.pxi
            right = exact_pow_helper(&dummyL, self.prime_pow.ram_prec_cap, _right, self.prime_pow)
            mpz_init(tmp)
            try:
                mpz_mul_si(tmp, right.value, self.ordp)
                if overunderflow_mpz(&ans.ordp, tmp, ans.unit, ans.prime_pow):
                    return ans
                else:
                    ans.ordp = mpz_get_si(tmp)
            finally:
                mpz_clear(tmp)
            cpow(ans.unit, self.unit, right.value, ans.prime_pow.ram_prec_cap, ans.prime_pow)
        else:
            # padic_pow_helper is defined in padic_template_element.pxi
            dummyL = padic_pow_helper(ans.unit, self.unit, self.ordp, self.prime_pow.ram_prec_cap,
                                      pright.unit, pright.ordp, pright.prime_pow.ram_prec_cap, self.prime_pow)
            ans.ordp = 0
        return ans

    cdef pAdicTemplateElement _lshift_c(self, long shift):
        r"""
        Multiply ``self`` by `\pi^{shift}`.

        Negative shifts may truncate the result if the parent is not a
        field.

        EXAMPLES:

        We create a floating point ring::

            sage: R = ZpFP(5, 20); a = R(1000); a
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
        cdef FPElement ans = self._new_c()
        # check both in case of overflow in sum; this case also includes self.ordp = maxordp
        if very_pos_val(shift) or very_pos_val(self.ordp + shift):
            # need to check that we're not shifting infinity
            if very_neg_val(self.ordp):
                raise ZeroDivisionError("Cannot multiply zero by infinity")
            ans.ordp = maxordp
            csetzero(ans.unit, ans.prime_pow)
        else:
            ans.ordp = self.ordp + shift
            ccopy(ans.unit, self.unit, ans.prime_pow)
        return ans

    cdef pAdicTemplateElement _rshift_c(self, long shift):
        r"""
        Divide by `\pi^{shift}`.

        Positive shifts may truncate the result if the parent is not a
        field.

        EXAMPLES::

            sage: R = ZpFP(997, 7); a = R(123456878908); a
            964*997 + 572*997^2 + 124*997^3

        Shifting to the right divides by a power of `\pi`, but
        dropping terms with negative valuation::

            sage: a >> 3
            124

        A negative shift multiplies by that power of `\pi`::

            sage: a >> -3
            964*997^4 + 572*997^5 + 124*997^6
        """
        if shift == 0:
            return self
        elif very_pos_val(self.ordp):
            if very_pos_val(shift):
                raise ZeroDivisionError("Cannot divide zero by zero")
            return self
        elif very_neg_val(self.ordp):
            if very_neg_val(shift):
                raise ZeroDivisionError("Cannot divide infinity by infinity")
            return self
        cdef FPElement ans = self._new_c()
        cdef long diff
        if self.prime_pow.in_field == 1 or shift <= self.ordp:
            if very_pos_val(shift):
                ans._set_infinity()
            elif very_neg_val(shift):
                ans._set_exact_zero()
            else:
                ans.ordp = self.ordp - shift
                ccopy(ans.unit, self.unit, ans.prime_pow)
        else:
            diff = shift - self.ordp
            if diff >= self.prime_pow.ram_prec_cap:
                ans._set_exact_zero()
            else:
                ans.ordp = 0
                cshift(ans.unit, ans.prime_pow.shift_rem, self.unit, -diff, ans.prime_pow.ram_prec_cap, ans.prime_pow, False)
                ans._normalize()
        return ans

    def _repr_(self, mode=None, do_latex=False):
        """
        Return a string representation of this element.

        INPUT:

        - ``mode`` -- allows one to override the default print mode of
          the parent (default: ``None``)

        - ``do_latex`` -- whether to return a latex representation or
          a normal one

        EXAMPLES::

            sage: ZpFP(5,5)(1/3)  # indirect doctest
            2 + 3*5 + 5^2 + 3*5^3 + 5^4
            sage: ~QpFP(5,5)(0)
            infinity
        """
        if very_neg_val(self.ordp):
            return "infinity"
        return self.parent()._printer.repr_gen(self, do_latex, mode=mode)

    def add_bigoh(self, absprec):
        r"""
        Return a new element truncated modulo `\pi^{\mbox{absprec}}`.

        INPUT:

        - ``absprec`` -- integer or infinity

        OUTPUT: a new element truncated modulo `\pi^{\mbox{absprec}}`

        EXAMPLES::

            sage: R = Zp(7,4,'floating-point','series'); a = R(8); a.add_bigoh(1)
            1
        """
        cdef long aprec, newprec
        if absprec is infinity or very_neg_val(self.ordp):
            return self
        elif isinstance(absprec, int):
            aprec = absprec
        else:
            if not isinstance(absprec, Integer):
                absprec = Integer(absprec)
            if mpz_fits_slong_p((<Integer>absprec).value) == 0:
                if mpz_sgn((<Integer>absprec).value) > 0:
                    return self
                aprec = minusmaxordp
            else:
                aprec = mpz_get_si((<Integer>absprec).value)
        if aprec >= self.ordp + self.prime_pow.ram_prec_cap:
            return self
        cdef FPElement ans = self._new_c()
        if aprec <= self.ordp:
            ans._set_exact_zero()
            if aprec < 0:
                return self.parent().fraction_field()(0)
        else:
            ans.ordp = self.ordp
            creduce(ans.unit, self.unit, aprec - self.ordp, ans.prime_pow)
        return ans

    cpdef bint _is_exact_zero(self) except -1:
        """
        Test whether this element is exactly zero.

        EXAMPLES::

            sage: R = Zp(7,4,'floating-point','series'); a = R(8); a._is_exact_zero()
            False
            sage: b = R(0); b._is_exact_zero()
            True
        """
        return very_pos_val(self.ordp)

    cpdef bint _is_inexact_zero(self) except -1:
        """
        Return ``True`` if ``self`` is indistinguishable from zero.

        EXAMPLES::

            sage: R = ZpFP(7, 5)
            sage: R(14)._is_inexact_zero()
            False
            sage: R(0)._is_inexact_zero()
            True
        """
        return very_pos_val(self.ordp)

    def is_zero(self, absprec = None):
        r"""
        Return whether ``self`` is zero modulo `\pi^{\mbox{absprec}}`.

        INPUT:

        - ``absprec`` -- integer

        EXAMPLES::

            sage: R = ZpFP(17, 6)
            sage: R(0).is_zero()
            True
            sage: R(17^6).is_zero()
            False
            sage: R(17^2).is_zero(absprec=2)
            True
        """
        if absprec is None:
            return very_pos_val(self.ordp)
        if very_pos_val(self.ordp):
            return True
        if absprec is infinity:
            return False
        if isinstance(absprec, int):
            return self.ordp >= absprec
        if not isinstance(absprec, Integer):
            absprec = Integer(absprec)
        return mpz_cmp_si((<Integer>absprec).value, self.ordp) <= 0

    def __bool__(self):
        """
        Return ``True`` if this element is distinguishable from zero.

        For most applications, explicitly specifying the power of `p`
        modulo which the element is supposed to be nonzero is preferable.

        EXAMPLES::

            sage: R = ZpFP(5); a = R(0); b = R(75)
            sage: bool(a), bool(b)  # indirect doctest
            (False, True)
        """
        return not very_pos_val(self.ordp)

    def is_equal_to(self, _right, absprec=None):
        r"""
        Return whether this element is equal to ``right`` modulo `p^{\mbox{absprec}}`.

        If ``absprec`` is ``None``, determines whether ``self`` and ``right``
        have the same value.

        INPUT:

        - ``right`` -- a `p`-adic element with the same parent
        - ``absprec`` -- positive integer or ``None`` (default: ``None``)

        EXAMPLES::

            sage: R = ZpFP(2, 6)
            sage: R(13).is_equal_to(R(13))
            True
            sage: R(13).is_equal_to(R(13+2^10))
            True
            sage: R(13).is_equal_to(R(17), 2)
            True
            sage: R(13).is_equal_to(R(17), 5)
            False
        """
        cdef FPElement right
        cdef long aprec, rprec
        if self.parent() is _right.parent():
            right = _right
        else:
            right = self.parent().coerce(_right)
        if very_neg_val(self.ordp):
            if very_neg_val(right.ordp):
                return True
            return False
        elif very_neg_val(right.ordp):
            return False
        if absprec is None or absprec is infinity:
            return ((self.ordp == right.ordp) and
                    (ccmp(self.unit, right.unit, self.prime_pow.ram_prec_cap, False, False, self.prime_pow) == 0))
        if not isinstance(absprec, Integer):
            absprec = Integer(absprec)
        if mpz_cmp_si((<Integer>absprec).value, self.ordp) <= 0:
            if mpz_cmp_si((<Integer>absprec).value, right.ordp) <= 0:
                return True
            return False
        elif mpz_cmp_si((<Integer>absprec).value, right.ordp) <= 0:
            return False
        if self.ordp != right.ordp:
            return False
        if mpz_cmp_si((<Integer>absprec).value, maxordp) >= 0:
            return ccmp(self.unit, right.unit, self.prime_pow.ram_prec_cap, False, False, self.prime_pow) == 0
        aprec = mpz_get_si((<Integer>absprec).value)
        rprec = aprec - self.ordp
        if rprec > self.prime_pow.ram_prec_cap:
            rprec = self.prime_pow.ram_prec_cap
        return ccmp(self.unit,
                    right.unit,
                    rprec,
                    rprec < self.prime_pow.ram_prec_cap,
                    rprec < right.prime_pow.ram_prec_cap,
                    self.prime_pow) == 0

    cdef int _cmp_units(self, pAdicGenericElement _right) except -2:
        """
        Comparison of units, used in equality testing.

        EXAMPLES::

            sage: R = ZpFP(5)
            sage: a = R(17); b = R(0,3); c = R(85,7); d = R(2, 1)
            sage: any([a == b, a == c, b == c, b == d, c == d, a == d])  # indirect doctest
            False
            sage: all([a == a, b == b, c == c, d == d])
            True
        """
        cdef FPElement right = _right
        return ccmp(self.unit, right.unit, self.prime_pow.ram_prec_cap, False, False, self.prime_pow)

    cdef pAdicTemplateElement lift_to_precision_c(self, long absprec):
        """
        Lift this element to another with precision at least absprec.

        Since floating point elements don't track precision, this
        function just returns the same element.

        EXAMPLES::

            sage: R = ZpFP(5)
            sage: a = R(77, 2); a
            2
            sage: a.lift_to_precision(17)  # indirect doctest
            2
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

            sage: R = ZpFP(17,5); a = R(11)
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
        if self.ordp > 0:
            self._set_exact_zero()
        elif self.ordp < 0:
            raise ValueError("cannot set negative valuation element to Teichmuller representative.")
        else:
            cteichmuller(self.unit, self.unit, self.prime_pow.ram_prec_cap, self.prime_pow)

    def _polynomial_list(self, pad=False):
        """
        Return the coefficient list for a polynomial over the base ring
        yielding this element.

        INPUT:

        - ``pad`` -- whether to pad the result with zeros of the appropriate precision

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<x> = ZZ[]
            sage: K.<a> = QqFP(5^3)
            sage: W.<w> = K.extension(x^3 - 5)
            sage: (1 + w)._polynomial_list()
            [1, 1]
            sage: (1 + w)._polynomial_list(pad=True)
            [1, 1, 0]
        """
        R = self.base_ring()
        if very_pos_val(self.ordp):
            L = []
        elif very_neg_val(self.ordp):
            L = [~R(0)]
        else:
            L = ccoefficients(self.unit, self.ordp, self.prime_pow.ram_prec_cap, self.prime_pow)
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
            sage: K.<a> = QqFP(5^3)
            sage: a.polynomial()
            x
            sage: a.polynomial(var='y')
            y
            sage: (5*a^2 + K(25, 4)).polynomial()
            5*x^2 + 5^2
        """
        R = self.base_ring()
        S = R[var]
        return S(self._polynomial_list())

    def precision_absolute(self):
        """
        The absolute precision of this element.

        EXAMPLES::

            sage: R = Zp(7,4,'floating-point'); a = R(7); a.precision_absolute()
            5
            sage: R(0).precision_absolute()
            +Infinity
            sage: (~R(0)).precision_absolute()
            -Infinity
        """
        if very_pos_val(self.ordp):
            return infinity
        elif very_neg_val(self.ordp):
            return -infinity
        cdef Integer ans = PY_NEW(Integer)
        mpz_set_si(ans.value, self.ordp + self.prime_pow.ram_prec_cap)
        return ans

    def precision_relative(self):
        r"""
        The relative precision of this element.

        EXAMPLES::

            sage: R = Zp(7,4,'floating-point'); a = R(7); a.precision_relative()
            4
            sage: R(0).precision_relative()
            0
            sage: (~R(0)).precision_relative()
            0
        """
        cdef Integer ans = PY_NEW(Integer)
        if huge_val(self.ordp):
            mpz_set_si(ans.value, 0)
        else:
            mpz_set_si(ans.value, self.prime_pow.ram_prec_cap)
        return ans

    cpdef pAdicTemplateElement unit_part(FPElement self):
        r"""
        Return the unit part of this element.

        If the valuation of this element is positive, then the high
        digits of the result will be zero.

        EXAMPLES::

            sage: R = Zp(17, 4, 'floating-point')
            sage: R(5).unit_part()
            5
            sage: R(18*17).unit_part()
            1 + 17
            sage: R(0).unit_part()
            Traceback (most recent call last):
            ...
            ValueError: unit part of 0 and infinity not defined
            sage: type(R(5).unit_part())
            <class 'sage.rings.padics.padic_floating_point_element.pAdicFloatingPointElement'>
            sage: R = ZpFP(5, 5); a = R(75); a.unit_part()
            3
        """
        if huge_val(self.ordp):
            raise ValueError("unit part of 0 and infinity not defined")
        cdef FPElement ans = (<FPElement>self)._new_c()
        ans.ordp = 0
        ccopy(ans.unit, (<FPElement>self).unit, ans.prime_pow)
        return ans

    cdef long valuation_c(self) noexcept:
        """
        Return the valuation of this element.

        If this element is an exact zero, returns ``maxordp``, which is defined as
        ``(1L << (sizeof(long) * 8 - 2))-1``.

        If this element is infinity, returns ``-maxordp``.

        TESTS::

            sage: R = ZpFP(5, 5); R(1).valuation()  # indirect doctest
            0
            sage: R = Zp(17, 4,'floating-point')
            sage: a = R(2*17^2)
            sage: a.valuation()
            2
            sage: R = Zp(5, 4,'floating-point')
            sage: R(0).valuation()
            +Infinity
            sage: (~R(0)).valuation()
            -Infinity
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
        return self.ordp

    cpdef val_unit(self, p=None):
        """
        Return a 2-tuple, the first element set to the valuation of
        this element, and the second to the unit part.

        If this element is either zero or infinity, raises an error.

        EXAMPLES::

            sage: R = ZpFP(5,5)
            sage: a = R(75); b = a - a
            sage: a.val_unit()
            (2, 3)
            sage: b.val_unit()
            Traceback (most recent call last):
            ...
            ValueError: unit part of 0 and infinity not defined
        """
        if p is not None and p != self.parent().prime():
            raise ValueError('Ring (%s) residue field of the wrong characteristic.' % self.parent())
        if huge_val(self.ordp):
            raise ValueError("unit part of 0 and infinity not defined")
        cdef Integer valuation = PY_NEW(Integer)
        mpz_set_si(valuation.value, self.ordp)
        cdef FPElement unit = self._new_c()
        unit.ordp = 0
        ccopy(unit.unit, self.unit, unit.prime_pow)
        return valuation, unit

    def __hash__(self):
        """
        Hashing.

        EXAMPLES::

            sage: R = ZpFP(11, 5)
            sage: hash(R(3)) == hash(3)
            True
        """
        if very_pos_val(self.ordp):
            return 0
        if very_neg_val(self.ordp):
            return 314159
        return chash(self.unit, self.ordp, self.prime_pow.ram_prec_cap, self.prime_pow) ^ self.ordp

cdef class pAdicCoercion_ZZ_FP(RingHomomorphism):
    """
    The canonical inclusion from the integer ring to a floating point ring.

    EXAMPLES::

        sage: f = ZpFP(5).coerce_map_from(ZZ); f
        Ring morphism:
          From: Integer Ring
          To:   5-adic Ring with floating precision 20

    TESTS::

        sage: TestSuite(f).run()
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = ZpFP(5).coerce_map_from(ZZ); type(f)
            <class 'sage.rings.padics.padic_floating_point_element.pAdicCoercion_ZZ_FP'>
        """
        RingHomomorphism.__init__(self, ZZ.Hom(R))
        self._zero = R.element_class(R, 0)
        self._section = pAdicConvert_FP_ZZ(R)

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpFP(5).coerce_map_from(ZZ)
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

            sage: f = ZpFP(5).coerce_map_from(ZZ)
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

            sage: f = ZpFP(5).coerce_map_from(ZZ)
            sage: f(0).parent()
            5-adic Ring with floating precision 20
            sage: f(5)
            5
        """
        if mpz_sgn((<Integer>x).value) == 0:
            return self._zero
        cdef FPElement ans = self._zero._new_c()
        ans.ordp = cconv_mpz_t(ans.unit, (<Integer>x).value, ans.prime_pow.ram_prec_cap, False, ans.prime_pow)
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in (relative or absolute or both).

        INPUT:

        - ``x`` -- an Integer

        - ``absprec``, or the first positional argument -- the maximum
          absolute precision (unused for floating point elements)

        - ``relprec``, or the second positional argument -- the
          maximum relative precision (unused for floating point
          elements)

        EXAMPLES::

            sage: R = ZpFP(5,4)
            sage: type(R(10,2))
            <class 'sage.rings.padics.padic_floating_point_element.pAdicFloatingPointElement'>
            sage: R(30,2)
            5
            sage: R(30,3,2)
            5 + 5^2
            sage: R(30,absprec=2)
            5
            sage: R(30,relprec=2)
            5 + 5^2
            sage: R(30,absprec=1)
            0
            sage: R(30,empty=True)
            0
        """
        cdef long val, aprec, rprec
        if mpz_sgn((<Integer>x).value) == 0:
            return self._zero
        _process_args_and_kwds(&aprec, &rprec, args, kwds, False, self._zero.prime_pow)
        val = get_ordp(x, self._zero.prime_pow)
        if aprec - val < rprec:
            rprec = aprec - val
        if rprec <= 0:
            return self._zero
        cdef FPElement ans = self._zero._new_c()
        ans.ordp = cconv_mpz_t(ans.unit, (<Integer>x).value, rprec, False, ans.prime_pow)
        return ans

    def section(self):
        r"""
        Return a map back to `\ZZ` that approximates an element of this
        `p`-adic ring by an integer.

        EXAMPLES::

            sage: f = ZpFP(5).coerce_map_from(ZZ).section()
            sage: f(ZpFP(5)(-1)) - 5^20
            -1
        """
        from sage.misc.constant_function import ConstantFunction
        if not isinstance(self._section.domain, ConstantFunction):
            import copy
            self._section = copy.copy(self._section)
        return self._section


cdef class pAdicConvert_FP_ZZ(RingMap):
    r"""
    The map from a floating point ring back to `\ZZ` that returns the smallest
    nonnegative integer approximation to its input which is accurate up to the precision.

    If the input is not in the closure of the image of `\ZZ`, raises a :exc:`ValueError`.

    EXAMPLES::

        sage: f = ZpFP(5).coerce_map_from(ZZ).section(); f
        Set-theoretic ring morphism:
          From: 5-adic Ring with floating precision 20
          To:   Integer Ring
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = ZpFP(5).coerce_map_from(ZZ).section(); type(f)
            <class 'sage.rings.padics.padic_floating_point_element.pAdicConvert_FP_ZZ'>
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

            sage: f = QpFP(5).coerce_map_from(ZZ).section()
            sage: f(QpFP(5)(-1)) - 5^20
            -1
            sage: f(QpFP(5)(5))
            5
            sage: f(QpFP(5)(0))
            0
            sage: f(~QpFP(5)(5))
            Traceback (most recent call last):
            ...
            ValueError: negative valuation
            sage: f(~QpFP(5)(0))
            Traceback (most recent call last):
            ...
            ValueError: Infinity cannot be converted to a rational
        """
        cdef Integer ans = PY_NEW(Integer)
        cdef FPElement x = _x
        if very_pos_val(x.ordp):
            mpz_set_ui(ans.value, 0)
        elif very_neg_val(x.ordp):
            raise ValueError("Infinity cannot be converted to a rational")
        else:
            cconv_mpz_t_out(ans.value, x.unit, x.ordp, x.prime_pow.ram_prec_cap, x.prime_pow)
        return ans

cdef class pAdicCoercion_QQ_FP(RingHomomorphism):
    """
    The canonical inclusion from the rationals to a floating point field.

    EXAMPLES::

        sage: f = QpFP(5).coerce_map_from(QQ); f
        Ring morphism:
          From: Rational Field
          To:   5-adic Field with floating precision 20

    TESTS::

        sage: TestSuite(f).run()
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = QpFP(5).coerce_map_from(QQ); type(f)
            <class 'sage.rings.padics.padic_floating_point_element.pAdicCoercion_QQ_FP'>
        """
        RingHomomorphism.__init__(self, QQ.Hom(R))
        self._zero = R.element_class(R, 0)
        self._section = pAdicConvert_FP_QQ(R)

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = QpFP(5).coerce_map_from(QQ)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Ring morphism:
              From: Rational Field
              To:   5-adic Field with floating precision 20
            sage: g == f
            True
            sage: g is f
            False
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

            sage: f = QpFP(5).coerce_map_from(QQ)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Ring morphism:
              From: Rational Field
              To:   5-adic Field with floating precision 20
            sage: g == f
            True
            sage: g is f
            False
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

            sage: f = QpFP(5).coerce_map_from(QQ)
            sage: f(0).parent()
            5-adic Field with floating precision 20
            sage: f(1/5)
            5^-1
            sage: f(1/4)
            4 + 3*5 + 3*5^2 + 3*5^3 + 3*5^4 + 3*5^5 + 3*5^6 + 3*5^7 + 3*5^8 + 3*5^9 + 3*5^10 + 3*5^11 + 3*5^12 + 3*5^13 + 3*5^14 + 3*5^15 + 3*5^16 + 3*5^17 + 3*5^18 + 3*5^19
            sage: f(1/4, 5)
            4 + 3*5 + 3*5^2 + 3*5^3 + 3*5^4
        """
        if mpq_sgn((<Rational>x).value) == 0:
            return self._zero
        cdef FPElement ans = self._zero._new_c()
        ans.ordp = cconv_mpq_t(ans.unit, (<Rational>x).value, ans.prime_pow.ram_prec_cap, False, self._zero.prime_pow)
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in
        (relative or absolute or both).

        See the documentation for
        :meth:`pAdicCappedRelativeElement.__init__` for more details.

        EXAMPLES::

            sage: R = QpFP(5,4)
            sage: type(R(10/3,2))
            <class 'sage.rings.padics.padic_floating_point_element.pAdicFloatingPointElement'>
            sage: R(10/3,2)
            4*5
            sage: R(10/3,3,1)
            4*5
            sage: R(10/3,absprec=2)
            4*5
            sage: R(10/3,relprec=2)
            4*5 + 5^2
            sage: R(10/3,absprec=1)
            0
            sage: R(3/100,absprec=-1)
            2*5^-2
        """
        cdef long val, aprec, rprec
        cdef FPElement ans
        _process_args_and_kwds(&aprec, &rprec, args, kwds, False, self._zero.prime_pow)
        if mpq_sgn((<Rational>x).value) == 0:
            return self._zero
        val = get_ordp(x, self._zero.prime_pow)
        if aprec <= val:
            return self._zero
        ans = self._zero._new_c()
        rprec = min(rprec, aprec - val)
        ans.ordp = cconv_mpq_t(ans.unit, (<Rational>x).value, rprec, False, self._zero.prime_pow)
        return ans

    def section(self):
        """
        Return a map back to the rationals that approximates an element by
        a rational number.

        EXAMPLES::

            sage: f = QpFP(5).coerce_map_from(QQ).section()
            sage: f(QpFP(5)(1/4))
            1/4
            sage: f(QpFP(5)(1/5))
            1/5
        """
        from sage.misc.constant_function import ConstantFunction
        if not isinstance(self._section.domain, ConstantFunction):
            import copy
            self._section = copy.copy(self._section)
        return self._section

cdef class pAdicConvert_FP_QQ(RingMap):
    """
    The map from the floating point ring back to the rationals that returns a
    rational approximation of its input.

    EXAMPLES::

        sage: f = QpFP(5).coerce_map_from(QQ).section(); f
        Set-theoretic ring morphism:
          From: 5-adic Field with floating precision 20
          To:   Rational Field
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = QpFP(5).coerce_map_from(QQ).section(); type(f)
            <class 'sage.rings.padics.padic_floating_point_element.pAdicConvert_FP_QQ'>
            sage: f.category()
            Category of homsets of sets with partial maps
        """
        RingMap.__init__(self, Hom(R, QQ, SetsWithPartialMaps()))

    cpdef Element _call_(self, _x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = QpFP(5).coerce_map_from(QQ).section()
            sage: f(QpFP(5)(-1))
            -1
            sage: f(QpFP(5)(0))
            0
            sage: f(QpFP(5)(1/5))
            1/5
        """
        cdef Rational ans = Rational.__new__(Rational)
        cdef FPElement x = _x
        if very_pos_val(x.ordp):
            mpq_set_ui(ans.value, 0, 1)
        elif very_neg_val(x.ordp):
            raise ValueError("Infinity cannot be converted to a rational")
        else:
            cconv_mpq_t_out(ans.value, x.unit, x.ordp, x.prime_pow.ram_prec_cap, x.prime_pow)
        return ans

cdef class pAdicConvert_QQ_FP(Morphism):
    r"""
    The inclusion map from `\QQ` to a floating point ring.

    EXAMPLES::

        sage: f = ZpFP(5).convert_map_from(QQ); f
        Generic morphism:
          From: Rational Field
          To:   5-adic Ring with floating precision 20
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = ZpFP(5).convert_map_from(QQ); type(f)
            <class 'sage.rings.padics.padic_floating_point_element.pAdicConvert_QQ_FP'>
        """
        Morphism.__init__(self, Hom(QQ, R, SetsWithPartialMaps()))
        self._zero = R.element_class(R, 0)

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpFP(5).convert_map_from(QQ)
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

            sage: f = ZpFP(5).convert_map_from(QQ)
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

            sage: f = ZpFP(5,4).convert_map_from(QQ)
            sage: f(1/7)
            3 + 3*5 + 2*5^3
            sage: f(0/1)
            0
        """
        if mpq_sgn((<Rational>x).value) == 0:
            return self._zero
        cdef FPElement ans = self._zero._new_c()
        ans.ordp = cconv_mpq_t(ans.unit, (<Rational>x).value, ans.prime_pow.ram_prec_cap, False, ans.prime_pow)
        if ans.ordp < 0:
            raise ValueError("p divides the denominator")
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in (relative or absolute or both).

        INPUT:

        - ``x`` -- a Rational

        - ``absprec``, or the first positional argument -- the maximum
          absolute precision (unused for floating point elements)

        - ``relprec``, or the second positional argument -- the
          maximum relative precision (unused for floating point
          elements)

        EXAMPLES::

            sage: R = ZpFP(5,4)
            sage: type(R(1/7,2))
            <class 'sage.rings.padics.padic_floating_point_element.pAdicFloatingPointElement'>
            sage: R(1/7,2)
            3 + 3*5
            sage: R(1/7,3,2)
            3 + 3*5
            sage: R(1/7,absprec=2)
            3 + 3*5
            sage: R(5/7,relprec=2)
            3*5 + 3*5^2
            sage: R(1/7,absprec=1)
            3
            sage: R(1/7,empty=True)
            0
        """
        cdef long val, aprec, rprec
        if mpq_sgn((<Rational>x).value) == 0:
            return self._zero
        _process_args_and_kwds(&aprec, &rprec, args, kwds, False, self._zero.prime_pow)
        val = get_ordp(x, self._zero.prime_pow)
        rprec = min(rprec, aprec - val)
        if rprec <= 0:
            return self._zero
        cdef FPElement ans = self._zero._new_c()
        ans.ordp = cconv_mpq_t(ans.unit, (<Rational>x).value, rprec, False, ans.prime_pow)
        if ans.ordp < 0:
            raise ValueError("p divides the denominator")
        return ans

cdef class pAdicCoercion_FP_frac_field(RingHomomorphism):
    r"""
    The canonical inclusion of `\ZZ_q` into its fraction field.

    EXAMPLES::

        sage: # needs sage.libs.flint
        sage: R.<a> = ZqFP(27, implementation='FLINT')
        sage: K = R.fraction_field()
        sage: f = K.coerce_map_from(R); f
        Ring morphism:
          From: 3-adic Unramified Extension Ring in a defined by x^3 + 2*x + 1
          To:   3-adic Unramified Extension Field in a defined by x^3 + 2*x + 1

    TESTS::

        sage: TestSuite(f).run()                                                        # needs sage.libs.flint
    """
    def __init__(self, R, K):
        r"""
        Initialization.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFP(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R); type(f)
            <class 'sage.rings.padics.qadic_flint_FP.pAdicCoercion_FP_frac_field'>
        """
        RingHomomorphism.__init__(self, R.Hom(K))
        self._zero = K(0)
        self._section = pAdicConvert_FP_frac_field(K, R)

    cpdef Element _call_(self, _x):
        r"""
        Evaluation.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFP(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f(a)
            a
            sage: f(R(0))
            0
        """
        cdef FPElement x = _x
        cdef FPElement ans = self._zero._new_c()
        ans.ordp = x.ordp
        cshift_notrunc(ans.unit, x.unit, 0, ans.prime_pow.ram_prec_cap, x.prime_pow, False)
        IF CELEMENT_IS_PY_OBJECT:
            # The base ring is wrong, so we fix it.
            K = ans.unit.base_ring()
            ans.unit._coeffs = [K(c) for c in ans.unit._coeffs]
        return ans

    cpdef Element _call_with_args(self, _x, args=(), kwds={}):
        r"""
        This function is used when some precision cap is passed in
        (relative or absolute or both).

        See the documentation for
        :meth:`pAdicCappedAbsoluteElement.__init__` for more details.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFP(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f(a, 3)
            a
            sage: b = 9*a + 27
            sage: f(b, 3)
            a*3^2
            sage: f(b, 4, 1)
            a*3^2
            sage: f(b, 4, 3)
            a*3^2 + 3^3
            sage: f(b, absprec=4)
            a*3^2 + 3^3
            sage: f(b, relprec=3)
            a*3^2 + 3^3
            sage: f(b, absprec=1)
            0
            sage: f(R(0))
            0
        """
        cdef long aprec, rprec
        cdef FPElement x = _x
        cdef FPElement ans = self._zero._new_c()
        cdef bint reduce = False
        _process_args_and_kwds(&aprec, &rprec, args, kwds, False, ans.prime_pow)
        if aprec <= x.ordp:
            ans._set_exact_zero()
        else:
            if rprec < ans.prime_pow.ram_prec_cap:
                reduce = True
            else:
                rprec = ans.prime_pow.ram_prec_cap
            if aprec < rprec + x.ordp:
                rprec = aprec - x.ordp
                reduce = True
            ans.ordp = x.ordp
            cshift_notrunc(ans.unit, x.unit, 0, rprec, x.prime_pow, reduce)
            IF CELEMENT_IS_PY_OBJECT:
                # The base ring is wrong, so we fix it.
                K = ans.unit.base_ring()
                ans.unit._coeffs = [K(c) for c in ans.unit._coeffs]
        return ans

    def section(self):
        r"""
        Return a map back to the ring that converts elements of
        nonnegative valuation.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFP(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f(K.gen())
            a
        """
        return self._section

    cdef dict _extra_slots(self):
        r"""
        Helper for copying and pickling.

        TESTS::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFP(27, implementation='FLINT')
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
        r"""
        Helper for copying and pickling.

        TESTS::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFP(9, implementation='FLINT')
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

cdef class pAdicConvert_FP_frac_field(Morphism):
    r"""
    The section of the inclusion from `\ZZ_q` to its fraction field.

    EXAMPLES::

        sage: # needs sage.libs.flint
        sage: R.<a> = ZqFP(27, implementation='FLINT')
        sage: K = R.fraction_field()
        sage: f = R.convert_map_from(K); f
        Generic morphism:
          From: 3-adic Unramified Extension Field in a defined by x^3 + 2*x + 1
          To:   3-adic Unramified Extension Ring in a defined by x^3 + 2*x + 1
    """
    def __init__(self, K, R):
        r"""
        Initialization.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFP(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K); type(f)
            <class 'sage.rings.padics.qadic_flint_FP.pAdicConvert_FP_frac_field'>
        """
        Morphism.__init__(self, Hom(K, R, SetsWithPartialMaps()))
        self._zero = R(0)

    cpdef Element _call_(self, _x):
        r"""
        Evaluation.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFP(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K)
            sage: f(K.gen())
            a
        """
        cdef FPElement x = _x
        if x.ordp < 0:
            raise ValueError("negative valuation")
        cdef FPElement ans = self._zero._new_c()
        ans.ordp = x.ordp
        cshift_notrunc(ans.unit, x.unit, 0, ans.prime_pow.ram_prec_cap, ans.prime_pow, False)
        IF CELEMENT_IS_PY_OBJECT:
            # The base ring is wrong, so we fix it.
            K = ans.unit.base_ring()
            ans.unit._coeffs = [K(c) for c in ans.unit._coeffs]
        return ans

    cpdef Element _call_with_args(self, _x, args=(), kwds={}):
        r"""
        This function is used when some precision cap is passed in
        (relative or absolute or both).

        See the documentation for
        :meth:`pAdicCappedAbsoluteElement.__init__` for more details.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFP(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K); a = K(a)
            sage: f(a, 3)
            a
            sage: b = 9*a + 27
            sage: f(b, 3)
            a*3^2
            sage: f(b, 4, 1)
            a*3^2
            sage: f(b, 4, 3)
            a*3^2 + 3^3
            sage: f(b, absprec=4)
            a*3^2 + 3^3
            sage: f(b, relprec=3)
            a*3^2 + 3^3
            sage: f(b, absprec=1)
            0
            sage: f(K(0))
            0
        """
        cdef long aprec, rprec
        cdef FPElement x = _x
        if x.ordp < 0:
            raise ValueError("negative valuation")
        cdef FPElement ans = self._zero._new_c()
        cdef bint reduce = False
        _process_args_and_kwds(&aprec, &rprec, args, kwds, False, ans.prime_pow)
        if aprec <= x.ordp:
            ans._set_exact_zero()
        else:
            if rprec < ans.prime_pow.ram_prec_cap:
                reduce = True
            else:
                rprec = ans.prime_pow.ram_prec_cap
            if aprec < rprec + x.ordp:
                rprec = aprec - x.ordp
                reduce = True
            ans.ordp = x.ordp
            cshift_notrunc(ans.unit, x.unit, 0, rprec, x.prime_pow, reduce)
            IF CELEMENT_IS_PY_OBJECT:
                # The base ring is wrong, so we fix it.
                K = ans.unit.base_ring()
                ans.unit._coeffs = [K(c) for c in ans.unit._coeffs]
        return ans

    cdef dict _extra_slots(self):
        r"""
        Helper for copying and pickling.

        TESTS::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFP(27, implementation='FLINT')
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
        r"""
        Helper for copying and pickling.

        TESTS::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFP(9, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K)
            sage: a = K(a)
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


def unpickle_fpe_v2(cls, parent, unit, ordp):
    """
    Unpickles a floating point element.

    EXAMPLES::

        sage: from sage.rings.padics.padic_floating_point_element import pAdicFloatingPointElement, unpickle_fpe_v2
        sage: R = ZpFP(5)
        sage: a = unpickle_fpe_v2(pAdicFloatingPointElement, R, 17, 2); a
        2*5^2 + 3*5^3
        sage: a.parent() is R
        True
    """
    cdef FPElement ans = cls.__new__(cls)
    ans._parent = parent
    ans.prime_pow = <PowComputer_?>parent.prime_pow
    IF CELEMENT_IS_PY_OBJECT:
        polyt = type(ans.prime_pow.modulus)
        ans.unit = <celement>polyt.__new__(polyt)
    cconstruct(ans.unit, ans.prime_pow)
    cunpickle(ans.unit, unit, ans.prime_pow)
    ans.ordp = ordp
    return ans
