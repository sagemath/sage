"""
p-Adic Generic Element

Elements of `p`-Adic Rings and Fields

AUTHORS:

- David Roe

- Genya Zaytman: documentation

- David Harvey: doctests

- Julian Rueth (2013-02-14): fixes for exp() and log()

"""

#*****************************************************************************
#       Copyright (C) 2007 David Roe <roed@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/gmp.pxi"
include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"

import sys

cimport sage.rings.padics.local_generic_element
from sage.rings.padics.local_generic_element cimport LocalGenericElement
from sage.rings.rational cimport Rational
from sage.rings.integer cimport Integer
from sage.rings.infinity import infinity
import sage.rings.rational_field

cdef long maxordp = (1L << (sizeof(long) * 8 - 2)) - 1

cdef class pAdicGenericElement(LocalGenericElement):
    def __richcmp__(left, right, int op):
        """
        Comparison.

        EXAMPLES::

            sage: R = Zp(5); a = R(5, 6); b = R(5 + 5^6, 8)
            sage: a == b #indirect doctest
            True
        """
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
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
        """
        m = min(left.precision_absolute(), right.precision_absolute())
        x_ordp = left.valuation()
        if x_ordp >= m :
            x_ordp = infinity
        y_ordp = right.valuation()
        if y_ordp >= m :
            y_ordp = infinity
        if x_ordp < y_ordp:
            return -1
        elif x_ordp > y_ordp:
            return 1
        else:  # equal ordp
            if x_ordp is infinity:
                return 0 # since both are zero
            else:
                return (<pAdicGenericElement>left.unit_part())._cmp_units(right.unit_part())

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
        raise TypeError, "this type of p-adic does not support exact zeros"

    cpdef bint _is_exact_zero(self) except -1:
        """
        Returns True if self is exactly zero.  Since
        non-capped-relative elements cannot be exact, this function
        always returns False.

        EXAMPLES::

            sage: ZpCA(17)(0)._is_exact_zero()
            False
        """
        return False

    cpdef bint _is_inexact_zero(self) except -1:
        """
        Returns True if self is indistinguishable from zero, but not
        exactly zero.

        EXAMPLES::

            sage: Zp(5)(0,5)._is_inexact_zero()
            True
        """
        raise NotImplementedError

    cpdef bint _is_zero_rep(self) except -1:
        """
        Returns True is self is indistinguishable from zero.

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

    #def _pari_(self):
    #    """
    #    Returns a pari version of this element.

    #    EXAMPLES::

    #        sage: R = Zp(5)
    #        sage: pari(R(1777))
    #        2 + 5^2 + 4*5^3 + 2*5^4 + O(5^20)
    #    """
    #    return pari(self._pari_init_())

    def __floordiv__(self, right):
        """
        Divides self by right and throws away the nonintegral part if
        self.parent() is not a field.

        There are a number of reasonable definitions for floor
        division.  Any definition should satisfy the following
        identity:

        (1) a = (a // b) * b + a % b

        If a and b lie in a field, then setting a % b = 0 and a // b =
        a / b provides a canonical way of satisfying this equation.

        However, for elements of integer rings, there are many choices
        of definitions for a // b and a % b that satisfy this
        equation.  Since p-adic rings in Sage come equipped with a
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


        In Sage we choose option (3), mainly because it is more easily
        defined in terms of shifting and thus generalizes more easily
        to extension rings.

        EXAMPLES::

            sage: R = ZpCA(5); a = R(129378); b = R(2398125)
            sage: a // b #indirect doctest
            3 + 3*5 + 4*5^2 + 2*5^4 + 2*5^6 + 4*5^7 + 5^9 + 5^10 + 5^11 + O(5^12)
            sage: a / b
            4*5^-4 + 3*5^-3 + 2*5^-2 + 5^-1 + 3 + 3*5 + 4*5^2 + 2*5^4 + 2*5^6 + 4*5^7 + 5^9 + 5^10 + 5^11 + O(5^12)
            sage: a % b
            3 + 5^4 + 3*5^5 + 2*5^6 + 4*5^7 + 5^8 + O(5^16)
            sage: (a // b) * b + a % b
            3 + 2*5^4 + 5^5 + 3*5^6 + 5^7 + O(5^16)

            The alternative definition:

            sage: a
            3 + 2*5^4 + 5^5 + 3*5^6 + 5^7 + O(5^20)
            sage: c = ((a - 3)>>4)/b.unit_part(); c
            1 + 2*5 + 2*5^3 + 4*5^4 + 5^6 + 5^7 + 5^8 + 4*5^9 + 2*5^10 + 4*5^11 + 4*5^12 + 2*5^13 + 3*5^14 + O(5^16)
            sage: c*b + 3
            3 + 2*5^4 + 5^5 + 3*5^6 + 5^7 + O(5^20)
        """
        if right == 0:
            raise ZeroDivisionError
        P = self.parent()
        if P.is_field():
            return self / right
        else:
            right = P(right)
            v, u = right.val_unit()
            return P(self / u) >> v


    def __getitem__(self, n):
        r"""
        Returns the coefficient of `p^n` in the series expansion of this
        element, as an integer in the range `0` to `p-1`.

        EXAMPLES::

            sage: R = Zp(7,4,'capped-rel','series'); a = R(1/3); a
            5 + 4*7 + 4*7^2 + 4*7^3 + O(7^4)
            sage: a[0] #indirect doctest
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
            IndexError: list index out of range
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

            sage: S.<a> = Qq(125)
            sage: a[-2]
            []

        .. SEEALSO::

            :meth:`sage.rings.padics.local_generic_element.LocalGenericElement.slice`
        """
        if isinstance(n, slice):
            return self.slice(n.start, n.stop, n.step)
        if self.parent().f() == 1:
            zero = Integer(0)
        else:
            zero = []
        if n < self.valuation():
            return zero
        if n >= self.precision_absolute():
            raise IndexError("list index out of range")

        if self.parent().is_field():
            n -= self.valuation()

        # trailing coefficients which are zero are not stored in self.list() -
        # we catch an IndexError to check for this.
        try:
            return self.list()[n]
        except IndexError:
            return zero

    def __invert__(self):
        r"""
        Returns the multiplicative inverse of self.

        EXAMPLE::

            sage: R = Zp(7,4,'capped-rel','series'); a = R(3); a
            3 + O(7^4)
            sage: ~a #indirect doctest
            5 + 4*7 + 4*7^2 + 4*7^3 + O(7^4)

        .. NOTE::

            The element returned is an element of the fraction field.
        """
        return self.parent().fraction_field()(self, relprec = self.precision_relative()).__invert__()

    def __mod__(self, right):
        """
        If self is in a field, returns 0.  If in a ring, returns a
        p-adic integer such that

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
        equation.  Since p-adic rings in Sage come equipped with a
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


        In Sage we choose option (3), mainly because it is more easily
        defined in terms of shifting and thus generalizes more easily
        to extension rings.

        EXAMPLES::

            sage: R = ZpCA(5); a = R(129378); b = R(2398125)
            sage: a // b
            3 + 3*5 + 4*5^2 + 2*5^4 + 2*5^6 + 4*5^7 + 5^9 + 5^10 + 5^11 + O(5^12)
            sage: a / b
            4*5^-4 + 3*5^-3 + 2*5^-2 + 5^-1 + 3 + 3*5 + 4*5^2 + 2*5^4 + 2*5^6 + 4*5^7 + 5^9 + 5^10 + 5^11 + O(5^12)
            sage: a % b #indirect doctest
            3 + 5^4 + 3*5^5 + 2*5^6 + 4*5^7 + 5^8 + O(5^16)

            The alternative definition:

            sage: a
            3 + 2*5^4 + 5^5 + 3*5^6 + 5^7 + O(5^20)
            sage: c = ((a - 3)>>4)/b.unit_part(); c
            1 + 2*5 + 2*5^3 + 4*5^4 + 5^6 + 5^7 + 5^8 + 4*5^9 + 2*5^10 + 4*5^11 + 4*5^12 + 2*5^13 + 3*5^14 + O(5^16)
            sage: c*b + 3
            3 + 2*5^4 + 5^5 + 3*5^6 + 5^7 + O(5^20)
        """
        if right == 0:
            raise ZeroDivisionError
        if self.parent().is_field():
            return self.parent()(0)
        else:
            return self - (self // right) * right

    #def _is_exact_zero(self):
    #    return False

    #def _is_inexact_zero(self):
    #    return self.is_zero() and not self._is_exact_zero()

    def str(self, mode=None):
        """
        Returns a string representation of self.

        EXAMPLES::

            sage: Zp(5,5,print_mode='bars')(1/3).str()[3:]
            '1|3|1|3|2'
        """
        return self._repr_(mode=mode)

    def _repr_(self, mode=None, do_latex=False):
        """
        Returns a string representation of self.

        INPUTS:

        - ``mode`` -- allows one to override the default print mode of
          the parent (default: ``None``).

        - ``do_latex`` -- whether to return a latex representation or
          a normal one.

        EXAMPLES::

            sage: Zp(5,5)(1/3) # indirect doctest
            2 + 3*5 + 5^2 + 3*5^3 + 5^4 + O(5^5)
        """
        return self.parent()._printer.repr_gen(self, do_latex, mode=mode)

    def additive_order(self, prec):
        r"""
        Returns the additive order of self, where self is considered
        to be zero if it is zero modulo `p^{\mbox{prec}}`.

        INPUT:

        - ``self`` -- a p-adic element
        - ``prec`` -- an integer

        OUTPUT:

        integer -- the additive order of self

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

    def algdep(self, n):
        """
        Returns a polynomial of degree at most `n` which is approximately
        satisfied by this number. Note that the returned polynomial need not be
        irreducible, and indeed usually won't be if this number is a good
        approximation to an algebraic number of degree less than `n`.

        ALGORITHM: Uses the PARI C-library ``algdep`` command.

        INPUT:

        - ``self`` -- a p-adic element
        - ``n`` -- an integer

        OUTPUT:

        polynomial -- degree n polynomial approximately satisfied by self

        EXAMPLES::

            sage: K = Qp(3,20,'capped-rel','series'); R = Zp(3,20,'capped-rel','series')
            sage: a = K(7/19); a
            1 + 2*3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^8 + 2*3^9 + 3^11 + 3^12 + 2*3^15 + 2*3^16 + 3^17 + 2*3^19 + O(3^20)
            sage: a.algdep(1)
            19*x - 7
            sage: K2 = Qp(7,20,'capped-rel')
            sage: b = K2.zeta(); b.algdep(2)
            x^2 - x + 1
            sage: K2 = Qp(11,20,'capped-rel')
            sage: b = K2.zeta(); b.algdep(4)
            x^4 - x^3 + x^2 - x + 1
            sage: a = R(7/19); a
            1 + 2*3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^8 + 2*3^9 + 3^11 + 3^12 + 2*3^15 + 2*3^16 + 3^17 + 2*3^19 + O(3^20)
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
        return sage.rings.arith.algdep(self, n)

    def algebraic_dependency(self, n):
        """
        Returns a polynomial of degree at most `n` which is approximately
        satisfied by this number.  Note that the returned polynomial need not
        be irreducible, and indeed usually won't be if this number is a good
        approximation to an algebraic number of degree less than `n`.

        ALGORITHM: Uses the PARI C-library algdep command.

        INPUT:

        - ``self`` -- a p-adic element
        - ``n`` -- an integer

        OUTPUT:

        polynomial -- degree n polynomial approximately satisfied by self

        EXAMPLES::

            sage: K = Qp(3,20,'capped-rel','series'); R = Zp(3,20,'capped-rel','series')
            sage: a = K(7/19); a
            1 + 2*3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^8 + 2*3^9 + 3^11 + 3^12 + 2*3^15 + 2*3^16 + 3^17 + 2*3^19 + O(3^20)
            sage: a.algebraic_dependency(1)
            19*x - 7
            sage: K2 = Qp(7,20,'capped-rel')
            sage: b = K2.zeta(); b.algebraic_dependency(2)
            x^2 - x + 1
            sage: K2 = Qp(11,20,'capped-rel')
            sage: b = K2.zeta(); b.algebraic_dependency(4)
            x^4 - x^3 + x^2 - x + 1
            sage: a = R(7/19); a
            1 + 2*3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^8 + 2*3^9 + 3^11 + 3^12 + 2*3^15 + 2*3^16 + 3^17 + 2*3^19 + O(3^20)
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

    #def gamma(self):
    #    raise NotImplementedError

    def is_square(self): #should be overridden for lazy elements
        """
        Returns whether self is a square

        INPUT:

        - ``self`` -- a p-adic element

        OUTPUT:

        boolean -- whether self is a square

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
        if self._is_exact_zero() or self._is_inexact_zero():
            return True
        elif self.parent().prime() != 2:
            return (self.valuation() % 2 == 0) and (self.unit_part().residue(1).is_square())
        else:
            #won't work for general extensions...
            return (self.valuation() % 2 == 0) and (self.unit_part().residue(3) == 1)

    #def log_artin_hasse(self):
    #    raise NotImplementedError

    def multiplicative_order(self, prec = None): #needs to be rewritten for lazy elements
        r"""
        Returns the multiplicative order of self, where self is
        considered to be one if it is one modulo `p^{\mbox{prec}}`.

        INPUT:

        - ``self`` -- a p-adic element
        - ``prec`` -- an integer

        OUTPUT:

        - integer -- the multiplicative order of self

        EXAMPLES::

            sage: K = Qp(5,20,'capped-rel')
            sage: K(-1).multiplicative_order(20)
            2
            sage: K(1).multiplicative_order(20)
            1
            sage: K(2).multiplicative_order(20)
            +Infinity
            sage: K(3).multiplicative_order(20)
            +Infinity
            sage: K(4).multiplicative_order(20)
            +Infinity
            sage: K(5).multiplicative_order(20)
            +Infinity
            sage: K(25).multiplicative_order(20)
            +Infinity
            sage: K(1/5).multiplicative_order(20)
            +Infinity
            sage: K(1/25).multiplicative_order(20)
            +Infinity
            sage: K.zeta().multiplicative_order(20)
            4

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
        if self.valuation() != 0:
            return infinity
        res = self.residue(1)
        if self.is_equal_to(self.parent().teichmuller(res.lift()),prec): #should this be made more efficient?
            return res.multiplicative_order()
        else:
            return infinity

    def valuation(self, p = None):
        """
        Returns the valuation of self.

        INPUT:

        - ``self`` -- a p-adic element
        - ``p`` -- a prime (default: None). If specified, will make sure that p==self.parent().prime()

        NOTE: The optional argument p is used for consistency with the valuation methods on integer and rational.

        OUTPUT:

        integer -- the valuation of self

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

            sage: K.<a> = Qq(25)
            sage: K(0).valuation()
            +Infinity

            sage: R(1/50).valuation(5)
            -2
            sage: R(1/50).valuation(3)
            Traceback (most recent call last):
            ...
            ValueError: Ring (5-adic Field with capped relative precision 4) residue field of the wrong characteristic.
        """
        if not p is None and p != self.parent().prime():
            raise ValueError, 'Ring (%s) residue field of the wrong characteristic.'%self.parent()
        cdef long v = self.valuation_c()
        if v == maxordp:
            return infinity
        cdef Integer ans = PY_NEW(Integer)
        mpz_set_si(ans.value, v)
        return ans

    cdef long valuation_c(self):
        """
        This function is overridden in subclasses to provide an
        actual implementation of valuation.

        For exact zeros, ``maxordp`` is returned, rather than infinity.

        EXAMPLES:

        For example, the valuation function on pAdicCappedRelativeElements
        uses an overridden version of this function.

        ::

            sage: Zp(5)(5).valuation() #indirect doctest
            1
        """
        raise NotImplementedError

    cpdef val_unit(self):
        """
        Return ``(self.valuation(), self.unit_part())``. To be overridden in
        derived classes.

        EXAMPLES::

            sage: Zp(5,5)(5).val_unit()
            (1, 1 + O(5^5))
        """
        raise NotImplementedError

    def ordp(self, p = None):
        r"""
        Returns the valuation of self, normalized so that the valuation of `p` is 1

        INPUT:

        - ``self`` -- a p-adic element
        - ``p`` -- a prime (default: ``None``). If specified, will make sure that ``p == self.parent().prime()``

        NOTE: The optional argument p is used for consistency with the valuation methods on integer and rational.


        OUTPUT:

        integer -- the valuation of self, normalized so that the valuation of `p` is 1

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
        return self.valuation(p) / self.parent().ramification_index()

    def rational_reconstruction(self):
        r"""
        Returns a rational approximation to this p-adic number

        INPUT:

        - ``self`` -- a p-adic element

        OUTPUT:

        rational -- an approximation to self

        EXAMPLES::

            sage: R = Zp(5,20,'capped-rel')
            sage: for i in range(11):
            ...       for j in range(1,10):
            ...           if j == 5:
            ...               continue
            ...           assert i/j == R(i/j).rational_reconstruction()
        """
        if self.is_zero(self.precision_absolute()):
            return Rational(0)
        p = self.parent().prime()
        alpha = self.unit_part().lift()
        m = Integer(p**self.precision_relative())
        r = sage.rings.arith.rational_reconstruction(alpha, m)
        return (Rational(p)**self.valuation())*r

    def _shifted_log(self, aprec, mina=0):
        r"""
        Return ``-\log(1-self)`` for elements of positive valuation.

        This is a helper method for :meth:`log`.

        INPUT:

        - ``aprec`` -- an integer, the precision to which the result is
          correct. ``aprec`` must not exceed the precision cap of the ring over
          which this element is defined.

        - ``mina`` -- an integer (default: 0), the series will check `n` up to
          this valuation (and beyond) to see if they can contribute to the
          series.

        ALGORITHM:

        The computation uses the series

        .. MATH::

            -\log(1-x)=\sum_{n=1}^\infty \frac{x^n}{n}.

        For the result to be correct to precision ``aprec``, we sum all terms
        for which the valuation of `x^n/n` is stricly smaller than ``aprec``.

        EXAMPLES::

            sage: r = Qp(5,prec=4)(5)
            sage: r._shifted_log(2)
            5 + O(5^2)
            sage: r._shifted_log(4)
            5 + 3*5^2 + 4*5^3 + O(5^4)
            sage: r._shifted_log(100)
            5 + 3*5^2 + 4*5^3 + O(5^5)

            sage: r = Zp(5,prec=4,type='fixed-mod')(5)
            sage: r._shifted_log(5)
            5 + 3*5^2 + 4*5^3 + O(5^4)

        Only implemented for elements of positive valuation::

            sage: r = Zp(5,prec=4,type='fixed-mod')(1)
            sage: r._shifted_log(5)
            Traceback (most recent call last):
            ...
            ValueError: Input value (=1 + O(5^4)) must have strictly positive valuation

        """
        x = self
        R = self.parent()
        # to get the precision right over capped-absolute rings, we have to
        # work over the capped-relative fraction field
        if R.is_capped_absolute():
            R = R.fraction_field()
            x = R(x)

        alpha=x.valuation()
        if alpha<=0:
            raise ValueError('Input value (=%s) must have strictly positive valuation' % self)

        e=R.ramification_index()
        p=R.prime()

        # we sum all terms of the power series of log into total
        total=R.zero()

        # pre-compute x^p/p into x2p_p
        if R.is_capped_relative():
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
        a=0
        p2a=1       # p^a
        x2pa = x    # x^(p^a)

        while True:
            upper_u = ((aprec+a*e)/(alpha*p2a)).floor()
            # In the unramified case, we can stop summing terms as soon as
            # there are no u for a given a to sum over. In the ramified case,
            # it can happen that for some initial a there are no such u but
            # later in the series there are such u again. mina can be set to
            # take care of this by summing at least to a=mina-1
            if a >= mina and upper_u<=0:
                break
            # we compute the sum for the possible values for u using Horner's method
            inner_sum = R.zero()
            for u in xrange(upper_u,0,-1):
                # We want u to be a p-adic unit
                if u%p==0:
                    new_term = R.zero()
                else:
                    new_term = ~R(u)

                # This hack is to deal with rings that don't lift to fields
                if u>1 or x2p_p.is_zero():
                    inner_sum = (inner_sum+new_term)*x2pa
                else:
                    inner_sum = (inner_sum+new_term)*(x2p_p**a)*(x**(p2a-a*p))

            total += inner_sum

            # Now increase the power of p
            a += 1
            p2a = p2a*p
            x2pa = x2pa**p

        return total.add_bigoh(aprec)

    def log(self, p_branch=None, pi_branch=None, branch=None, aprec=None, change_frac=False):
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

        INPUTS:

        - ``p_branch`` -- an element in the base ring or its fraction
          field; the implementation will choose the branch of the
          logarithm which sends `p` to ``branch``.

        - ``pi_branch`` -- an element in the base ring or its fraction
          field; the implementation will choose the branch of the
          logarithm which sends the uniformizer to ``branch``.  You
          may specify at most one of ``p_branch`` and ``pi_branch``,
          and must specify one of them if this element is not a unit.

        - ``aprec`` -- an integer or ``None`` (default: ``None``) if not
          ``None``, then the result will only be correct to precision
          ``aprec``.

        - ``change_frac`` -- In general the codomain of the logarithm should be
          in the `p`-adic field, however, for most neighborhoods of 1, it lies
          in the ring of integers. This flag decides if the codomain should be
          the same as the input (default) or if it should change to the
          fraction field of the input.

        NOTES:

        What some other systems do:

        - PARI: Seems to define the logarithm for units not congruent
          to 1 as we do.

        - MAGMA: Only implements logarithm for 1-units (as of version 2.19-2)

        .. TODO::

        There is a soft-linear time algorith for logarithm described
        by Dan Berstein at
        http://cr.yp.to/lineartime/multapps-20041007.pdf

        ALGORITHM:

        1. Take the unit part `u` of the input.

        2. Raise `u` to `q-1` where `q` is the inertia degree of the ring
        extension, to obtain a 1-unit.

        3. Use the series expansion

        .. MATH::

            \log(1-x) = -x - 1/2 x^2 - 1/3 x^3 - 1/4 x^4 - 1/5 x^5 - \cdots

        to compute the logarithm `\log(u)`.

        4. Divide the result by ``q-1`` and multiply by ``self.valuation()*log(pi)``

        EXAMPLES::

            sage: Z13 = Zp(13, 10)
            sage: a = Z13(14); a
            1 + 13 + O(13^10)

        Note that the relative precision decreases when we take log -- it is
        the absolute precision that is preserved::

            sage: a.log()
            13 + 6*13^2 + 2*13^3 + 5*13^4 + 10*13^6 + 13^7 + 11*13^8 + 8*13^9 + O(13^10)
            sage: Q13 = Qp(13, 10)
            sage: a = Q13(14); a
            1 + 13 + O(13^10)
            sage: a.log()
            13 + 6*13^2 + 2*13^3 + 5*13^4 + 10*13^6 + 13^7 + 11*13^8 + 8*13^9 + O(13^10)

        The next few examples illustrate precision when computing `p`-adic
        logarithms::

            sage: R = Zp(5,10)
            sage: e = R(389); e
            4 + 2*5 + 3*5^3 + O(5^10)
            sage: e.log()
            2*5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^5 + 3*5^7 + 2*5^8 + 4*5^9 + O(5^10)
            sage: K = Qp(5,10)
            sage: e = K(389); e
            4 + 2*5 + 3*5^3 + O(5^10)
            sage: e.log()
            2*5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^5 + 3*5^7 + 2*5^8 + 4*5^9 + O(5^10)

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
            ValueError: You must specify a branch of the logarithm for non-units
            sage: b.log(p_branch=4)
            4 + O(5^10)
            sage: c = R(10)
            sage: c.log(p_branch=4)
            4 + 2*5 + 3*5^2 + 2*5^3 + 4*5^4 + 2*5^6 + 2*5^7 + 4*5^8 + 2*5^9 + O(5^10)

        The branch parameters are only relevant for elements of non-zero
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
            sage: W.<w> = R.ext(f)
            sage: z = 1 + w^2 + 4*w^7; z
            1 + w^2 + 4*w^7 + O(w^20)
            sage: z.log()
            w^2 + 2*w^4 + 3*w^6 + 4*w^7 + w^9 + 4*w^10 + 4*w^11 + 4*w^12 + 3*w^14 + w^15 + w^17 + 3*w^18 + 3*w^19 + O(w^20)

        In an extension, there will usually be a difference between
        specifying ``p_branch`` and ``pi_branch``::

            sage: b = W(5)
            sage: b.log()
            Traceback (most recent call last):
            ...
            ValueError: You must specify a branch of the logarithm for non-units
            sage: b.log(p_branch=0)
            O(w^20)
            sage: b.log(p_branch=w)
            w + O(w^20)
            sage: b.log(pi_branch=0)
            3*w^2 + 2*w^4 + 2*w^6 + 3*w^8 + 4*w^10 + w^13 + w^14 + 2*w^15 + 2*w^16 + w^18 + 4*w^19 + O(w^20)
            sage: b.unit_part().log()
            3*w^2 + 2*w^4 + 2*w^6 + 3*w^8 + 4*w^10 + w^13 + w^14 + 2*w^15 + 2*w^16 + w^18 + 4*w^19 + O(w^20)
            sage: y = w^2 * 4*w^7; y
            4*w^9 + O(w^29)
            sage: y.log(p_branch=0)
            2*w^2 + 2*w^4 + 2*w^6 + 2*w^8 + w^10 + w^12 + 4*w^13 + 4*w^14 + 3*w^15 + 4*w^16 + 4*w^17 + w^18 + 4*w^19 + O(w^20)
            sage: y.log(p_branch=w)
            w + 2*w^2 + 2*w^4 + 4*w^5 + 2*w^6 + 2*w^7 + 2*w^8 + 4*w^9 + w^10 + 3*w^11 + w^12 + 4*w^14 + 4*w^16 + 2*w^17 + w^19 + O(w^20)

        Check that log is multiplicative::

            sage: y.log(p_branch=0) + z.log() - (y*z).log(p_branch=0)
            O(w^20)

        Now an unramified example::

            sage: g = x^3 + 3*x + 3
            sage: A.<a> = R.ext(g)
            sage: b = 1 + 5*(1 + a^2) + 5^3*(3 + 2*a)
            sage: b.log()
            (a^2 + 1)*5 + (3*a^2 + 4*a + 2)*5^2 + (3*a^2 + 2*a)*5^3 + (3*a^2 + 2*a + 2)*5^4 + O(5^5)

        Check that log is multiplicative::

            sage: c = 3 + 5^2*(2 + 4*a)
            sage: b.log() + c.log() - (b*c).log()
            O(5^5)

        We illustrate the effect of the precision argument::

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
            sage: W.<w> = R.ext(f)
            sage: w.log(branch=2)
            Traceback (most recent call last):
            ...
            ValueError: logarithm is not integral, use change_frac=True to obtain a result in the fraction field
            sage: w.log(branch=2, change_frac=True)
            2*w^-3 + O(w^21)

        TESTS:

        Check that results are consistent over a range of precision::

            sage: max_prec = 40
            sage: p = 3
            sage: K = Zp(p, max_prec)
            sage: full_log = (K(1 + p)).log()
            sage: for prec in range(2, max_prec):
            ...       ll1 = (K(1+p).add_bigoh(prec)).log()
            ...       ll2 = K(1+p).log(prec)
            ...       assert ll1 == full_log
            ...       assert ll2 == full_log
            ...       assert ll1.precision_absolute() == prec

        Check that ``aprec`` works for fixed-mod elements::

            sage: R = ZpFM(7,10)
            sage: x = R(41152263); x
            5 + 3*7^2 + 4*7^3 + 3*7^4 + 5*7^5 + 6*7^6 + 7^9 + O(7^10)
            sage: x.log(aprec = 5)
            7 + 3*7^2 + 4*7^3 + 3*7^4 + O(7^10)
            sage: x.log(aprec = 7)
            7 + 3*7^2 + 4*7^3 + 3*7^4 + 7^5 + 3*7^6 + O(7^10)
            sage: x.log()
            7 + 3*7^2 + 4*7^3 + 3*7^4 + 7^5 + 3*7^6 + 7^7 + 3*7^8 + 4*7^9 + O(7^10)

        Check that precision is computed correctly in highly ramified
        extensions::

            sage: S.<x> = ZZ[]
            sage: K = Qp(5,5)
            sage: f = x^625 - 5*x - 5
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

            sage: z = 1 - w + O(w^2)
            sage: x = 1 - z
            sage: z.log().precision_absolute()
            -1625
            sage: (x^5/5).precision_absolute()
            -615
            sage: (x^25/25).precision_absolute()
            -1200
            sage: (x^125/125).precision_absolute()
            -1625
            sage: (x^625/625).precision_absolute()
            -1250

            sage: z.log().precision_relative()
            250

        AUTHORS:

        - William Stein: initial version

        - David Harvey (2006-09-13): corrected subtle precision bug (need to
          take denominators into account! -- see :trac:`53`)

        - Genya Zaytman (2007-02-14): adapted to new `p`-adic class

        - Amnon Besser, Marc Masdeu (2012-02-21): complete rewrite, valid for
          generic `p`-adic rings.

        - Soroosh Yazdani (2013-02-1): Fixed a precision issue in
          :meth:`_shifted_log`.  This should really fix the issue with
          divisions.

        - Julian Rueth (2013-02-14): Added doctests, some changes for
          capped-absolute implementations.

        """
        if self.is_zero():
            raise ValueError('logarithm is not defined at zero')
        if branch is not None:
            from sage.misc.superseded import deprecation
            deprecation(12575, "The keyword branch is deprecated.  Please use p_branch or pi_branch instead")
            p_branch = branch
        if p_branch is not None and pi_branch is not None:
            raise ValueError("You may only specify a branch of the logarithm in one way")
        R = self.parent()
        p = R.prime()
        q = p**R.f()

        if self.is_padic_unit():
            total = R.zero()
        else:
            if pi_branch is None:
                if p_branch is None:
                    raise ValueError("You must specify a branch of the logarithm for non-units")
                pi_branch = (p_branch - R._log_unit_part_p()) / R.e()
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
        e = R.e()
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
                deriv0 = (e / (minaprec * p.log(prec=53))).floor().exact_log(p)

                # These are the absolute precisions of x^(p^n) at potential minimum points
                L = [(minaprec * p**n - n * e, n) for n in [0, kink, deriv0, deriv0+1]]
                L.sort()
                minaprec = L[0][0]
                minn = L[0][1]

        if aprec is None or aprec > minaprec:
            aprec=minaprec

        retval = total - x._shifted_log(aprec, minn)*R(denom).inverse_of_unit()
        if not change_frac:
            if retval.valuation() < 0 and not R.is_field():
                raise ValueError("logarithm is not integral, use change_frac=True to obtain a result in the fraction field")
            retval=R(retval)
        return retval.add_bigoh(aprec)

    def exp(self, aprec = None):
        r"""
        Compute the `p`-adic exponential of this element if the exponential
        series converges.

        INPUT:

        - ``aprec`` -- an integer or ``None`` (default: ``None``); if
          specified, computes only up to the indicated precision.

        ALGORITHM: If self has a ``lift`` method (which should happen for
        elements of `\QQ_p` and `\ZZ_p`), then one uses the rule:
        `\exp(x)=\exp(p)^{x/p}` modulo the precision. The value of `\exp(p)` is
        precomputed. Otherwise, use the power series expansion of `\exp`,
        evaluating a certain number of terms which does about `O(\mbox{prec})`
        multiplications.

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
            sage: e = R(2*5 + 2*5**2 + 4*5**3 + 3*5**4 + 5**5 + 3*5**7 + 2*5**8 + 4*5**9).add_bigoh(10); e
            2*5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^5 + 3*5^7 + 2*5^8 + 4*5^9 + O(5^10)
            sage: e.exp()*R.teichmuller(4)
            4 + 2*5 + 3*5^3 + O(5^10)

        ::

            sage: K = Qp(5,10)
            sage: e = K(2*5 + 2*5**2 + 4*5**3 + 3*5**4 + 5**5 + 3*5**7 + 2*5**8 + 4*5**9).add_bigoh(10); e
            2*5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^5 + 3*5^7 + 2*5^8 + 4*5^9 + O(5^10)
            sage: e.exp()*K.teichmuller(4)
            4 + 2*5 + 3*5^3 + O(5^10)

        Logarithms and exponentials in extension fields. First, in an
        Eisenstein extension::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^4 + 15*x^2 + 625*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = 1 + w^2 + 4*w^7; z
            1 + w^2 + 4*w^7 + O(w^20)
            sage: z.log().exp()
            1 + w^2 + 4*w^7 + O(w^20)

        Now an unramified example::

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
            ...       ll = (K(p).add_bigoh(prec)).exp()
            ...       assert ll == full_exp
            ...       assert ll.precision_absolute() == prec
            sage: K = Qp(p, max_prec)
            sage: full_exp = (K(p)).exp()
            sage: for prec in range(2, max_prec):
            ...       ll = (K(p).add_bigoh(prec)).exp()
            ...       assert ll == full_exp
            ...       assert ll.precision_absolute() == prec

        Check that this also works for capped-absolute implementations::

            sage: Z13 = ZpCA(13, 10)
            sage: a = Z13(14); a
            1 + 13 + O(13^10)
            sage: a.log().exp()
            1 + 13 + O(13^10)

            sage: R = ZpCA(5,5)
            sage: S.<x> = R[]
            sage: f = x^4 + 15*x^2 + 625*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = 1 + w^2 + 4*w^7; z
            1 + w^2 + 4*w^7 + O(w^16)
            sage: z.log().exp()
            1 + w^2 + 4*w^7 + O(w^16)

        Check that this also works for fixed-mod implementations::

            sage: Z13 = ZpFM(13, 10)
            sage: a = Z13(14); a
            1 + 13 + O(13^10)
            sage: a.log().exp()
            1 + 13 + O(13^10)

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^4 + 15*x^2 + 625*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = 1 + w^2 + 4*w^7; z
            1 + w^2 + 4*w^7 + O(w^20)
            sage: z.log().exp()
            1 + w^2 + 4*w^7 + O(w^20)

        Some corner cases::

            sage: Z2 = Zp(2, 5)
            sage: Z2(2).exp()
            Traceback (most recent call last):
            ...
            ValueError: Exponential does not converge for that input.

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

        AUTHORS:

        - Genya Zaytman (2007-02-15)

        - Amnon Besser, Marc Masdeu (2012-02-23): Complete rewrite

        - Julian Rueth (2013-02-14): Added doctests, fixed some corner cases

        """
        p = self.parent().prime()

        if (p-1)*self.valuation() <= self.parent().ramification_index():
            raise ValueError('Exponential does not converge for that input.')

        if aprec is None or aprec > self.parent().precision_cap():
            aprec=self.parent().precision_cap()

        if hasattr(self,'lift'):
            y=Rational(self.lift())
            if p == 2:
                # in Z_2, the element has at least valuation 2, so we can
                # divide it by 4 (and use the exponential of 4 since the
                # exponential of 2 does not exist)
                p = 4
            return (self.parent()._exp_p()**Integer(y/p)).add_bigoh(min(aprec,self.precision_absolute()))
        else:
            return self._exp(aprec)

    def _exp(self, aprec):
        r"""
        Compute the exponential power series of this element, using Horner's
        evaluation and only one division.

        This is a helper method for :meth:`exp`.

        INPUT:

        - ``aprec`` -- an integer, the precison to which to compute the
          exponential

        EXAMPLES::

            sage: R.<w> = Zq(7^2,5)
            sage: x = R(7*w)
            sage: x._exp(5)
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
        e=self.parent().ramification_index()
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

        # the value of x^N+N*x^(N-1)+...+x*N!+N!
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

    def square_root(self, extend = True, all = False):
        r"""
        Returns the square root of this p-adic number

        INPUT:

        - ``self`` -- a p-adic element
        - ``extend`` -- bool (default: True); if True, return a square root in
          an extension if necessary; if False and no root exists in the given
          ring or field, raise a ValueError
        - ``all`` -- bool (default: False); if True, return a list of all
          square roots

        OUTPUT:

        p-adic element -- the square root of this p-adic number

        If ``all=False``, the square root chosen is the one whose
        reduction mod `p` is in the range `[0, p/2)`.

        EXAMPLES::

            sage: R = Zp(3,20,'capped-rel', 'val-unit')
            sage: R(0).square_root()
            0
            sage: R(1).square_root()
            1 + O(3^20)
            sage: R(2).square_root(extend = False)
            Traceback (most recent call last):
            ...
            ValueError: element is not a square
            sage: R(4).square_root() == R(-2)
            True
            sage: R(9).square_root()
            3 * 1 + O(3^21)

        When p = 2, the precision of the square root is one less
        than the input::

            sage: R2 = Zp(2,20,'capped-rel')
            sage: R2(0).square_root()
            0
            sage: R2(1).square_root()
            1 + O(2^19)
            sage: R2(4).square_root()
            2 + O(2^20)

            sage: R2(9).square_root() == R2(3, 19) or R2(9).square_root() == R2(-3, 19)
            True

            sage: R2(17).square_root()
            1 + 2^3 + 2^5 + 2^6 + 2^7 + 2^9 + 2^10 + 2^13 + 2^16 + 2^17 + O(2^19)

            sage: R3 = Zp(5,20,'capped-rel')
            sage: R3(0).square_root()
            0
            sage: R3(1).square_root()
            1 + O(5^20)
            sage: R3(-1).square_root() == R3.teichmuller(2) or R3(-1).square_root() == R3.teichmuller(3)
            True


        TESTS::

            sage: R = Qp(3,20,'capped-rel')
            sage: R(0).square_root()
            0
            sage: R(1).square_root()
            1 + O(3^20)
            sage: R(4).square_root() == R(-2)
            True
            sage: R(9).square_root()
            3 + O(3^21)
            sage: R(1/9).square_root()
            3^-1 + O(3^19)

            sage: R2 = Qp(2,20,'capped-rel')
            sage: R2(0).square_root()
            0
            sage: R2(1).square_root()
            1 + O(2^19)
            sage: R2(4).square_root()
            2 + O(2^20)
            sage: R2(9).square_root() == R2(3,19) or R2(9).square_root() == R2(-3,19)
            True
            sage: R2(17).square_root()
            1 + 2^3 + 2^5 + 2^6 + 2^7 + 2^9 + 2^10 + 2^13 + 2^16 + 2^17 + O(2^19)

            sage: R3 = Qp(5,20,'capped-rel')
            sage: R3(0).square_root()
            0
            sage: R3(1).square_root()
            1 + O(5^20)
            sage: R3(-1).square_root() == R3.teichmuller(2) or R3(-1).square_root() == R3.teichmuller(3)
            True
        """
        # need special case for zero since pari(self) is the *integer* zero
        # whose square root is a real number....!
        if self.valuation() is infinity:
            return self

        from sage.libs.pari.all import PariError
        try:
            # use pari
            ans = self.parent()(self._pari_().sqrt())
            if all:
                return [ans, -ans]
            else:
                return ans
        except PariError:
            # todo: should eventually change to return an element of
            # an extension field
            if extend:
                raise NotImplementedError, "extending using the sqrt function not yet implemented"
            elif all:
                return []
            else:
                raise ValueError, "element is not a square"

    #def _unit_part(self):
    #    raise NotImplementedError

    cpdef abs(self, prec=None):
        """
        Return the `p`-adic absolute value of ``self``.

        This is normalized so that the absolute value of `p` is `1/p`.

        INPUT:

        - ``prec`` -- Integer.  The precision of the real field in which
          the answer is returned.  If ``None``, returns a rational for
          absolutely unramified fields, or a real with 53 bits of
          precision for ramified fields.

        EXAMPLES::

            sage: a = Qp(5)(15); a.abs()
            1/5
            sage: a.abs(53)
            0.200000000000000
            sage: Qp(7)(0).abs()
            0
            sage: Qp(7)(0).abs(prec=20)
            0.00000

        An unramified extension::

            sage: R = Zp(5,5)
            sage: P.<x> = PolynomialRing(R)
            sage: Z25.<u> = R.ext(x^2 - 3)
            sage: u.abs()
            1
            sage: (u^24-1).abs()
            1/5

        A ramified extension::

            sage: W.<w> = R.ext(x^5 + 75*x^3 - 15*x^2 + 125*x - 5)
            sage: w.abs()
            0.724779663677696
            sage: W(0).abs()
            0.000000000000000
        """
        K = self.parent()
        if not prec and K.e() > 1:
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

