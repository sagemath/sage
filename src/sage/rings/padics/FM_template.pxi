"""
Fixed modulus template for complete discrete valuation rings.

In order to use this template you need to write a linkage file and
gluing file.  For an example see mpz_linkage.pxi (linkage file) and
padic_fixed_modulus_element.pyx (gluing file).

The linkage file implements a common API that is then used in the
class FMElement defined here.  See sage/libs/linkages/padics/API.pxi
for the functions needed.

The gluing file does the following:

- ctypedef's celement to be the appropriate type (e.g. mpz_t)
- includes the linkage file
- includes this template
- defines a concrete class inheriting from FMElement, and implements
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

from sage.ext.stdsage cimport PY_NEW
include "padic_template_element.pxi"
from cpython.int cimport *

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
        Sets the value of this element from given defining data.

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
            sage: a = R(17,5); a #indirect doctest
            2 + 3*5 + O(5^20)
            sage: R = ZpFM(5,5)
            sage: a = R(25/9); a #indirect doctest
            4*5^2 + 2*5^3 + O(5^5)
        """
        cconstruct(self.value, self.prime_pow)
        if isinstance(x,FMElement) and x.parent() is self.parent():
            cshift(self.value, (<FMElement>x).value, 0, 0, self.prime_pow, False)
        else:
            cconv(self.value, x, self.prime_pow.prec_cap, 0, self.prime_pow)

    cdef FMElement _new_c(self):
        """
        Creates a new element with the same basic info.

        TESTS::

            sage: R = ZpFM(5); R(6) * R(7) #indirect doctest
            2 + 3*5 + 5^2 + O(5^20)
        """
        cdef type t = type(self)
        cdef FMElement ans = t.__new__(t)
        ans._parent = self._parent
        ans.prime_pow = self.prime_pow
        cconstruct(ans.value, ans.prime_pow)
        return ans

    cdef int check_preccap(self) except -1:
        """
        Check that the precision of this element does not exceed the
        precision cap. Does nothing for fixed mod elements.

        TESTS::

            sage: ZpFM(5)(1).lift_to_precision(30) # indirect doctest
            1 + O(5^20)
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
            <type 'sage.rings.padics.padic_fixed_mod_element.pAdicFixedModElement'>
            sage: loads(dumps(a)) == a
            True
        """
        return unpickle_fme_v2, (self.__class__, self.parent(), cpickle(self.value, self.prime_pow))

    cpdef _neg_(self):
        r"""
        Return the additive inverse of this element.

        EXAMPLES::

            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: -R(7) #indirect doctest
            6*7 + 6*7^2 + 6*7^3 + O(7^4)
        """
        cdef FMElement ans = self._new_c()
        cneg(ans.value, self.value, ans.prime_pow.prec_cap, ans.prime_pow)
        creduce_small(ans.value, ans.value, ans.prime_pow.prec_cap, ans.prime_pow)
        return ans

    cpdef _add_(self, _right):
        r"""
        Return the sum of this element and ``_right``.

        EXAMPLES::

            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: x = R(1721); x
            6 + 5*7^3 + O(7^4)
            sage: y = R(1373); y
            1 + 4*7^3 + O(7^4)
            sage: x + y #indirect doctest
            7 + 2*7^3 + O(7^4)
        """
        cdef FMElement right = _right
        cdef FMElement ans = self._new_c()
        cadd(ans.value, self.value, right.value, ans.prime_pow.prec_cap, ans.prime_pow)
        creduce_small(ans.value, ans.value, ans.prime_pow.prec_cap, ans.prime_pow)
        return ans

    cpdef _sub_(self, _right):
        r"""
        Return the difference of this element and ``_right``.

        EXAMPLES::

            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: x = R(1721); x
            6 + 5*7^3 + O(7^4)
            sage: y = R(1373); y
            1 + 4*7^3 + O(7^4)
            sage: x - y #indirect doctest
            5 + 7^3 + O(7^4)
        """
        cdef FMElement right = _right
        cdef FMElement ans = self._new_c()
        csub(ans.value, self.value, right.value, ans.prime_pow.prec_cap, ans.prime_pow)
        creduce_small(ans.value, ans.value, ans.prime_pow.prec_cap, ans.prime_pow)
        return ans

    def __invert__(self):
        r"""
        Returns multiplicative inverse of this element. The valuation
        of ``self`` must be zero.

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
        if not cisunit(self.value, self.prime_pow):
            raise ValueError("cannot invert non-unit")
        cdef FMElement ans = self._new_c()
        cinvert(ans.value, self.value, ans.prime_pow.prec_cap, ans.prime_pow)
        return ans

    cpdef _mul_(self, _right):
        r"""
        Return the product of this element and ``_right``.

        EXAMPLES::

            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: R(3) * R(2) #indirect doctest
            6 + O(7^4)
            sage: R(1/2) * R(2)
            1 + O(7^4)
        """
        cdef FMElement right = _right
        cdef FMElement ans = self._new_c()
        cmul(ans.value, self.value, right.value, ans.prime_pow.prec_cap, ans.prime_pow)
        creduce(ans.value, ans.value, ans.prime_pow.prec_cap, ans.prime_pow)
        return ans

    cpdef _div_(self, _right):
        r"""
        Return the quotient of this element and ``right``. ``right`` must have
        valuation zero.

        EXAMPLES::

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
        cdef FMElement right = _right
        cdef FMElement ans = self._new_c()
        if not cisunit(right.value, self.prime_pow):
            raise ValueError("cannot invert non-unit")
        cdivunit(ans.value, self.value, right.value, ans.prime_pow.prec_cap, ans.prime_pow)
        creduce(ans.value, ans.value, ans.prime_pow.prec_cap, ans.prime_pow)
        return ans

    def __pow__(FMElement self, _right, dummy): # NOTE: dummy ignored, always use self.prime_pow.prec_cap
        """
        Exponentiation by an integer

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

        We check that :trac:`15640` is resolved::

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
        cpow(ans.value, self.value, right.value, self.prime_pow.prec_cap, self.prime_pow)
        return ans

    cdef pAdicTemplateElement _lshift_c(self, long shift):
        """
        Multiplies self by `\pi^{shift}`.

        If shift < -self.valuation(), digits will be truncated.  See
        :meth:`__rshift__` for details.

        EXAMPLES:

        We create a fixed modulus ring::

            sage: R = ZpFM(5, 20); a = R(1000); a
            3*5^3 + 5^4 + O(5^20)

        Shifting to the right is the same as dividing by a power of
        the uniformizer `\pi` of the `p`-adic ring.::

            sage: a >> 1
            3*5^2 + 5^3 + O(5^20)

        Shifting to the left is the same as multiplying by a power of
        `\pi`::

            sage: a << 2
            3*5^5 + 5^6 + O(5^20)
            sage: a*5^2
            3*5^5 + 5^6 + O(5^20)

        Shifting by a negative integer to the left is the same as
        right shifting by the absolute value::

            sage: a << -3
            3 + 5 + O(5^20)
            sage: a >> 3
            3 + 5 + O(5^20)
        """
        if shift < 0:
            return self._rshift_c(-shift)
        elif shift == 0:
            return self
        cdef FMElement ans = self._new_c()
        if shift >= self.prime_pow.prec_cap:
            csetzero(ans.value, ans.prime_pow)
        else:
            cshift(ans.value, self.value, shift, ans.prime_pow.prec_cap, ans.prime_pow, False)
        return ans

    cdef pAdicTemplateElement _rshift_c(self, long shift):
        """
        Divides by `\pi^{shift}`, and truncates.

        Note that this operation will insert arbitrary digits (in
        practice, currently all zero) in the least significant digits.

        EXAMPLES::

            sage: R = ZpFM(997, 7); a = R(123456878908); a
            964*997 + 572*997^2 + 124*997^3 + O(997^7)

        Shifting to the right divides by a power of `\pi`, but
        dropping terms with negative valuation::

            sage: a >> 3
            124 + O(997^7)

        A negative shift multiplies by that power of `\pi`::

            sage: a >> -3
            964*997^4 + 572*997^5 + 124*997^6 + O(997^7)
        """
        if shift < 0:
            return self._lshift_c(-shift)
        elif shift == 0:
            return self
        cdef FMElement ans = self._new_c()
        if shift >= self.prime_pow.prec_cap:
            csetzero(ans.value, ans.prime_pow)
        else:
            cshift(ans.value, self.value, -shift, ans.prime_pow.prec_cap, ans.prime_pow, False)
        return ans

    def add_bigoh(self, absprec):
        """
        Returns a new element truncated modulo `\pi^{\mbox{absprec}}`.

        INPUT:

        - ``absprec`` -- an integer

        OUTPUT:

            - a new element truncated modulo `\pi^{\mbox{absprec}}`.

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod','series'); a = R(8); a.add_bigoh(1)
            1 + O(7^4)
        """
        cdef long aprec, newprec
        if isinstance(absprec, int):
            aprec = absprec
        else:
            if not isinstance(absprec, Integer):
                absprec = Integer(absprec)
            aprec = mpz_get_si((<Integer>absprec).value)
        if aprec >= self.prime_pow.prec_cap:
            return self
        cdef FMElement ans = self._new_c()
        creduce(ans.value, self.value, aprec, ans.prime_pow)
        return ans

    cpdef bint _is_exact_zero(self) except -1:
        """
        Tests whether this element is an exact zero, which is always
        False for fixed modulus elements.

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod','series'); a = R(8); a._is_exact_zero()
            False
        """
        return False

    cpdef bint _is_inexact_zero(self) except -1:
        """
        Returns True if self is indistinguishable from zero.

        EXAMPLES::

            sage: R = ZpFM(7, 5)
            sage: R(14)._is_inexact_zero()
            False
            sage: R(0)._is_inexact_zero()
            True
        """
        return ciszero(self.value, self.prime_pow)

    def is_zero(self, absprec = None):
        r"""
        Returns whether self is zero modulo `\pi^{\mbox{absprec}}`.

        INPUT:

        - ``absprec`` -- an integer

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
        if mpz_cmp_si((<Integer>absprec).value, self.prime_pow.prec_cap) >= 0:
            return iszero
        cdef long val = self.valuation_c()
        return mpz_cmp_si((<Integer>absprec).value, val) <= 0

    def __nonzero__(self):
        """
        Returns True if this element is distinguishable from zero.

        For most applications, explicitly specifying the power of p
        modulo which the element is supposed to be nonzero is
        preferrable.

        EXAMPLES::

            sage: R = ZpFM(5); a = R(0); b = R(75)
            sage: bool(a), bool(b) # indirect doctest
            (False, True)
        """
        return not ciszero(self.value, self.prime_pow)

    def is_equal_to(self, _right, absprec=None):
        r"""
        Returns whether this element is equal to ``right`` modulo `p^{\mbox{absprec}}`.

        If ``absprec`` is ``None``, returns if ``self == 0``.

        INPUT:

        - ``right`` -- a p-adic element with the same parent
        - ``absprec`` -- a positive integer or ``None`` (default: ``None``)

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
            aprec = self.prime_pow.prec_cap
        else:
            if not isinstance(absprec, Integer):
                absprec = Integer(absprec)
            # If absprec is not positive, then self and right are always
            # equal.
            if mpz_sgn((<Integer>absprec).value) < 0:
                return True
            # If absprec is bigger than the precision cap, we use it
            # instead.
            if mpz_cmp_si((<Integer>absprec).value, self.prime_pow.prec_cap) >= 0:
                aprec = self.prime_pow.prec_cap
            else:
                aprec = mpz_get_si((<Integer>absprec).value)
        return ccmp(self.value,
                    right.value,
                    aprec,
                    aprec < self.prime_pow.prec_cap,
                    aprec < right.prime_pow.prec_cap,
                    self.prime_pow) == 0

    cdef int _cmp_units(self, pAdicGenericElement _right) except -2:
        """
        Comparison of units, used in equality testing.

        EXAMPLES::

            sage: R = ZpFM(5)
            sage: a = R(17); b = R(0,3); c = R(85,7); d = R(2, 1)
            sage: any([a == b, a == c, b == c, b == d, c == d, a == d]) # indirect doctest
            False
            sage: all([a == a, b == b, c == c, d == d])
            True
        """
        cdef FMElement right = _right
        return ccmp(self.value, right.value, self.prime_pow.prec_cap, False, False, self.prime_pow)

    cdef pAdicTemplateElement lift_to_precision_c(self, long absprec):
        """
        Lifts this element to another with precision at least absprec.

        Since fixed modulus elements don't track precision, this
        function just returns the same element.

        EXAMPLES::

            sage: R = ZpFM(5);
            sage: a = R(77, 2); a
            2 + 3*5^2 + O(5^20)
            sage: a.lift_to_precision(17) # indirect doctest
            2 + 3*5^2 + O(5^20)
        """
        return self

    def list(self, lift_mode = 'simple'):
        r"""
        Returns a list of coefficients of `\pi^i` starting with `\pi^0`.

        INPUT:

        - ``lift_mode`` -- ``'simple'``, ``'smallest'`` or ``'teichmuller'``
          (default: ``'simple'``:)

        OUTPUT:

        The list of coefficients of this element.

        .. NOTE::

            - Returns a list `[a_0, a_1, \ldots, a_n]` so that each `a_i`
              is an integer and `\sum_{i = 0}^n a_i \cdot p^i` is equal to
              this element modulo the precision cap.

            - If ``lift_mode`` is ``'simple'``, `0 \leq a_i < p`.

            - If ``lift_mode`` is ``'smallest'``, `-p/2 < a_i \leq p/2`.

            - If ``lift_mode`` is ``'teichmuller'``, `a_i^q = a_i`, modulo
              the precision cap.

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
        if ciszero(self.value, self.prime_pow):
            return []
        if lift_mode == 'teichmuller':
            return self.teichmuller_list()
        elif lift_mode == 'simple':
            return clist(self.value, self.prime_pow.prec_cap, True, self.prime_pow)
        elif lift_mode == 'smallest':
            return clist(self.value, self.prime_pow.prec_cap, False, self.prime_pow)
        else:
            raise ValueError("unknown lift_mode")

    def teichmuller_list(self):
        r"""
        Returns a list [`a_0`, `a_1`,..., `a_n`] such that

        - `a_i^q = a_i`

        - self.unit_part() = `\sum_{i = 0}^n a_i \pi^i`

        EXAMPLES::

            sage: R = ZpFM(5,5); R(14).list('teichmuller') #indirect doctest
            [4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5),
            3 + 3*5 + 2*5^2 + 3*5^3 + 5^4 + O(5^5),
            2 + 5 + 2*5^2 + 5^3 + 3*5^4 + O(5^5),
            1 + O(5^5),
            4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5)]
        """
        ans = PyList_New(0)
        if ciszero(self.value, self.prime_pow):
            return ans
        cdef long curpower = self.prime_pow.prec_cap
        cdef long prec_cap = self.prime_pow.prec_cap
        cdef FMElement list_elt
        cdef FMElement tmp = self._new_c()
        ccopy(tmp.value, self.value, self.prime_pow)
        while not ciszero(tmp.value, tmp.prime_pow) and curpower > 0:
            list_elt = self._new_c()
            cteichmuller(list_elt.value, tmp.value, prec_cap, self.prime_pow)
            if ciszero(list_elt.value, self.prime_pow):
                cshift_notrunc(tmp.value, tmp.value, -1, prec_cap, self.prime_pow)
            else:
                csub(tmp.value, tmp.value, list_elt.value, prec_cap, self.prime_pow)
                cshift_notrunc(tmp.value, tmp.value, -1, prec_cap, self.prime_pow)
                creduce(tmp.value, tmp.value, prec_cap, self.prime_pow)
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

            sage: R = ZpFM(17,5); a = R(11)
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
        if cisunit(self.value, self.prime_pow):
            cteichmuller(self.value, self.value, self.prime_pow.prec_cap, self.prime_pow)
        else:
            csetzero(self.value, self.prime_pow)

    def precision_absolute(self):
        """
        The absolute precision of this element.

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod'); a = R(7); a.precision_absolute()
            4
        """
        cdef Integer ans = PY_NEW(Integer)
        mpz_set_si(ans.value, self.prime_pow.prec_cap)
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
        cdef Integer ans = PY_NEW(Integer)
        mpz_set_si(ans.value, self.prime_pow.prec_cap - self.valuation_c())
        return ans

    cpdef pAdicTemplateElement unit_part(FMElement self):
        r"""
        Returns the unit part of self.

        If the valuation of self is positive, then the high digits of the
        result will be zero.

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
        cdef FMElement ans = (<FMElement>self)._new_c()
        cremove(ans.value, (<FMElement>self).value, (<FMElement>self).prime_pow.prec_cap, (<FMElement>self).prime_pow)
        return ans

    cdef long valuation_c(self):
        """
        Returns the valuation of this element.

        TESTS::

            sage: R = ZpFM(5, 5); R(0).valuation() #indirect doctest
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
        return cvaluation(self.value, self.prime_pow.prec_cap, self.prime_pow)

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
        cdef FMElement unit = self._new_c()
        cdef Integer valuation = PY_NEW(Integer)
        mpz_set_si(valuation.value, cremove(unit.value, self.value, self.prime_pow.prec_cap, self.prime_pow))
        return valuation, unit

    def __hash__(self):
        """
        Hashing.

        EXAMPLES::

            sage: R = ZpCA(11, 5)
            sage: hash(R(3)) == hash(3)
            True
        """
        return chash(self.value, 0, self.prime_pow.prec_cap, self.prime_pow)

cdef class pAdicCoercion_ZZ_FM(RingHomomorphism_coercion):
    """
    The canonical inclusion from ZZ to a fixed modulus ring.

    EXAMPLES::

        sage: f = ZpFM(5).coerce_map_from(ZZ); f
        Ring Coercion morphism:
          From: Integer Ring
          To:   5-adic Ring of fixed modulus 5^20
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = ZpFM(5).coerce_map_from(ZZ); type(f)
            <type 'sage.rings.padics.padic_fixed_mod_element.pAdicCoercion_ZZ_FM'>
        """
        RingHomomorphism_coercion.__init__(self, ZZ.Hom(R), check=False)
        self._zero = R._element_constructor(R, 0)
        self._section = pAdicConvert_FM_ZZ(R)

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpFM(5).coerce_map_from(ZZ)
            sage: g = copy(f) # indirect doctest
            sage: g == f
            True
            sage: g(6)
            1 + 5 + O(5^20)
            sage: g(6) == f(6)
            True
        """
        _slots['_zero'] = self._zero
        _slots['_section'] = self._section
        return RingHomomorphism_coercion._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpFM(5).coerce_map_from(ZZ)
            sage: g = copy(f) # indirect doctest
            sage: g == f
            True
            sage: g(6)
            1 + 5 + O(5^20)
            sage: g(6) == f(6)
            True
        """
        self._zero = _slots['_zero']
        self._section = _slots['_section']
        RingHomomorphism_coercion._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = ZpFM(5).coerce_map_from(ZZ)
            sage: f(0).parent()
            5-adic Ring of fixed modulus 5^20
            sage: f(5)
            5 + O(5^20)
        """
        if mpz_sgn((<Integer>x).value) == 0:
            return self._zero
        cdef FMElement ans = self._zero._new_c()
        cconv_mpz_t(ans.value, (<Integer>x).value, ans.prime_pow.prec_cap, True, ans.prime_pow)
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in (relative or absolute or both).

        INPUT:

        - ``x`` -- an Integer

        - ``absprec``, or the first positional argument -- the maximum
          absolute precision (unused for fixed modulus elements).

        - ``relprec``, or the second positional argument -- the
          maximum relative precision (unused for fixed modulus
          elements)

        EXAMPLES::

            sage: R = ZpFM(5,4)
            sage: type(R(10,2))
            <type 'sage.rings.padics.padic_fixed_mod_element.pAdicFixedModElement'>
            sage: R(30,2)
            5 + 5^2 + O(5^4)
            sage: R(30,3,1)
            5 + 5^2 + O(5^4)
            sage: R(30,absprec=2)
            5 + 5^2 + O(5^4)
            sage: R(30,relprec=2)
            5 + 5^2 + O(5^4)
            sage: R(30,absprec=1)
            5 + 5^2 + O(5^4)
            sage: R(30,empty=True)
            5 + 5^2 + O(5^4)
        """
        if mpz_sgn((<Integer>x).value) == 0:
            return self._zero
        cdef FMElement ans = self._zero._new_c()
        cconv_mpz_t(ans.value, (<Integer>x).value, ans.prime_pow.prec_cap, True, ans.prime_pow)
        return ans

    def section(self):
        """
        Returns a map back to ZZ that approximates an element of this
        `p`-adic ring by an integer.

        EXAMPLES::

            sage: f = ZpFM(5).coerce_map_from(ZZ).section()
            sage: f(ZpFM(5)(-1)) - 5^20
            -1
        """
        return self._section

cdef class pAdicConvert_FM_ZZ(RingMap):
    """
    The map from a fixed modulus ring back to ZZ that returns the the smallest
    non-negative integer approximation to its input which is accurate up to the precision.

    If the input is not in the closure of the image of ZZ, raises a ValueError.

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
            <type 'sage.rings.padics.padic_fixed_mod_element.pAdicConvert_FM_ZZ'>
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

            sage: f = ZpFM(5).coerce_map_from(ZZ).section()
            sage: f(ZpFM(5)(-1)) - 5^20
            -1
            sage: f(ZpFM(5)(0))
            0
        """
        cdef Integer ans = PY_NEW(Integer)
        cdef FMElement x = _x
        cconv_mpz_t_out(ans.value, x.value, 0, x.prime_pow.prec_cap, x.prime_pow)
        return ans

cdef class pAdicConvert_QQ_FM(Morphism):
    """
    The inclusion map from QQ to a fixed modulus ring that is defined
    on all elements with non-negative p-adic valuation.

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
            <type 'sage.rings.padics.padic_fixed_mod_element.pAdicConvert_QQ_FM'>
        """
        Morphism.__init__(self, Hom(QQ, R, SetsWithPartialMaps()))
        self._zero = R._element_constructor(R, 0)

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = ZpFM(5).convert_map_from(QQ)
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

            sage: f = ZpFM(5).convert_map_from(QQ)
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

            sage: f = ZpFM(5,4).convert_map_from(QQ)
            sage: f(1/7)
            3 + 3*5 + 2*5^3 + O(5^4)
            sage: f(0)
            O(5^4)
        """
        if mpq_sgn((<Rational>x).value) == 0:
            return self._zero
        cdef FMElement ans = self._zero._new_c()
        cconv_mpq_t(ans.value, (<Rational>x).value, ans.prime_pow.prec_cap, True, ans.prime_pow)
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in (relative or absolute or both).

        INPUT:

        - ``x`` -- a Rational

        - ``absprec``, or the first positional argument -- the maximum
          absolute precision (unused for fixed modulus elements).

        - ``relprec``, or the second positional argument -- the
          maximum relative precision (unused for fixed modulus
          elements)

        EXAMPLES::

            sage: R = ZpFM(5,4)
            sage: type(R(1/7,2))
            <type 'sage.rings.padics.padic_fixed_mod_element.pAdicFixedModElement'>
            sage: R(1/7,2)
            3 + 3*5 + 2*5^3 + O(5^4)
            sage: R(1/7,3,1)
            3 + 3*5 + 2*5^3 + O(5^4)
            sage: R(1/7,absprec=2)
            3 + 3*5 + 2*5^3 + O(5^4)
            sage: R(1/7,relprec=2)
            3 + 3*5 + 2*5^3 + O(5^4)
            sage: R(1/7,absprec=1)
            3 + 3*5 + 2*5^3 + O(5^4)
            sage: R(1/7,empty=True)
            3 + 3*5 + 2*5^3 + O(5^4)
        """
        if mpq_sgn((<Rational>x).value) == 0:
            return self._zero
        cdef FMElement ans = self._zero._new_c()
        cconv_mpq_t(ans.value, (<Rational>x).value, ans.prime_pow.prec_cap, True, ans.prime_pow)
        return ans

def unpickle_fme_v2(cls, parent, value):
    """
    Unpickles a capped relative element.

    EXAMPLES::

        sage: from sage.rings.padics.padic_fixed_mod_element import pAdicFixedModElement, unpickle_fme_v2
        sage: R = ZpFM(5)
        sage: a = unpickle_fme_v2(pAdicFixedModElement, R, 17*25); a
        2*5^2 + 3*5^3 + O(5^20)
        sage: a.parent() is R
        True
    """
    cdef FMElement ans = cls.__new__(cls)
    ans._parent = parent
    ans.prime_pow = <PowComputer_?>parent.prime_pow
    cconstruct(ans.value, ans.prime_pow)
    cunpickle(ans.value, value, ans.prime_pow)
    return ans
