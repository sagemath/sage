"""
This file gives a class from which all the `p`-adic templates inherit.

This file is included in the various precision template files (such as
sage/rings/padics/CR_template.pxi).  While this choice means that routines here
do have access to the inlined functions defined in the linkage files, the
downside is that there will be multiple ``pAdicTemplateElement``
classes, one for each `p`-adic template.

AUTHORS:

- David Roe -- initial version (2012-3-1)
"""
# ****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython.long cimport *

from sage.libs.gmp.all cimport *
import sage.rings.finite_rings.integer_mod
from cypari2.types cimport *
from cypari2.gen cimport Gen as pari_gen
from sage.libs.pari.convert_gmp cimport INT_to_mpz
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.padics.common_conversion cimport get_ordp, get_preccap
from sage.rings.integer cimport Integer
from sage.rings.infinity import infinity
from sage.rings.rational import Rational
from sage.rings.padics.precision_error import PrecisionError
from sage.rings.padics.misc import trim_zeros
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.structure.element import canonical_coercion

import itertools

cdef long maxordp = (1L << (sizeof(long) * 8 - 2)) - 1
cdef long minusmaxordp = -maxordp

cdef inline int check_ordp(long ordp) except -1:
    """
    Check for overflow after addition or subtraction of valuations.

    There is another variant, :meth:`check_ordp_mpz`, for ``mpz_t`` input.

    If overflow is detected, raises a :exc:`OverflowError`.
    """
    if ordp >= maxordp or ordp <= minusmaxordp:
        raise OverflowError("valuation overflow")

cdef class pAdicTemplateElement(pAdicGenericElement):
    """
    A class for common functionality among the `p`-adic template classes.

    INPUT:

    - ``parent`` -- a local ring or field

    - ``x`` -- data defining this element.  Various types are supported,
      including ints, Integers, Rationals, PARI `p`-adics, integers mod `p^k`
      and other Sage `p`-adics.

    - ``absprec`` -- a cap on the absolute precision of this element

    - ``relprec`` -- a cap on the relative precision of this element

    EXAMPLES::

        sage: Zp(17)(17^3, 8, 4)
        17^3 + O(17^7)
    """
    def __init__(self, parent, x, absprec=infinity, relprec=infinity):
        """
        Initialization.

        .. NOTE::

            This initialization function is not called for Integers
            and Rationals since a conversion morphism has been
            implemented.  It is, however, used for python ints and longs.

        EXAMPLES::

            sage: a = Zp(5)(1/2,3); a
            3 + 2*5 + 2*5^2 + O(5^3)
            sage: type(a)
            <class 'sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement'>
            sage: TestSuite(a).run()

        TESTS::

            sage: # needs sage.libs.ntl
            sage: QQq.<zz> = Qq(25,4)
            sage: FFp = Zp(5,5).residue_field()
            sage: QQq(FFp.zero())
            O(5)
            sage: QQq(FFp.one())
            1 + O(5)
            sage: QQq(IntegerModRing(25)(15))
            3*5 + O(5^2)
            sage: QQq(IntegerModRing(9)(0))
            Traceback (most recent call last):
            ...
            TypeError: p does not divide modulus 9

        ::

            sage: Zp(2)(Zp(5)(1))
            Traceback (most recent call last):
            ...
            TypeError: no conversion between padics when prime numbers differ
        """
        self.prime_pow = <PowComputer_?>parent.prime_pow
        pAdicGenericElement.__init__(self, parent)
        cdef long val, xprec
        cdef GEN pari_tmp
        if isinstance(x, int):
            x = Integer(x)
        elif isinstance(x, pari_gen):
            pari_tmp = (<pari_gen>x).g
            if typ(pari_tmp) == t_INT:
                x = PY_NEW(Integer)
                INT_to_mpz((<Integer>x).value, pari_tmp)
            elif typ(pari_tmp) == t_FRAC:
                x = Rational(x)
        elif isinstance(x, pAdicGenericElement):
            if self.prime_pow.prime != x.parent().prime():
                raise TypeError("no conversion between padics when prime numbers differ")
            if not ((<pAdicGenericElement>x)._is_base_elt(self.prime_pow.prime) or x.parent() is self.parent()):
                if x.parent().modulus().change_ring(self.base_ring()) == self.parent().modulus():
                    x = x.polynomial().change_ring(self.base_ring()).list()
                else:
                    x = self.base_ring()(x)
                    if x.is_zero():
                        absprec = min(absprec, x.precision_absolute()*self.prime_pow.e)
                        x = []
                    else:
                        x = [x]
        elif isinstance(x, sage.rings.finite_rings.integer_mod.IntegerMod_abstract):
            if not Integer(self.prime_pow.prime).divides(x.parent().order()):
                raise TypeError("p does not divide modulus %s" % x.parent().order())
        elif isinstance(x, Element) and isinstance(x.parent(), FiniteField):
            k = self.parent().residue_field()
            if not k.has_coerce_map_from(x.parent()):
                raise NotImplementedError("conversion from finite fields which do not embed into the residue field not implemented")

            x = k(x)
            if not k.is_prime_field():
                x = [k.prime_subfield()(c) for c in x.polynomial().list()]
                x = x + [k.prime_subfield().zero()] * (k.degree() - len(x))
        elif isinstance(x, (Integer, Rational, list, tuple)):
            pass
        elif isinstance(x, Polynomial) and x.variable_name() == self.parent().variable_name():
            x = x.list()
        else:
            x = Rational(x)
        val = get_ordp(x, self.prime_pow)
        if val < 0 and self.prime_pow.in_field == 0:
            raise ValueError("negative valuation")
        xprec = get_preccap(x, self.prime_pow)
        self._set(x, val, xprec, absprec, relprec)

    cdef int _set(self, x, long val, long xprec, absprec, relprec) except -1:
        """
        Set this element from given data computed in :meth:`__init__`.

        INPUT:

        - ``x`` -- an int, Integer, Rational, PARI `p`-adic, integer mod `p^k` or
          Sage `p`-adic

        - ``val`` -- a long; the valuation of ``x``

        - ``xprec`` -- a long; the cap on the absolute precision imposed by ``x``

        - ``absprec`` -- integer or infinity; a bound on the absolute precision

        - ``relprec`` -- integer or infinity; a bound on the relative precision
        """
        raise NotImplementedError

    cdef pAdicTemplateElement _new_with_value(self, celement value, long absprec):
        """
        Create a new element with a given value and absolute precision.

        Used by code that doesn't know the precision type.
        """
        raise NotImplementedError

    cdef int _get_unit(self, celement value) except -1:
        """
        Set ``value`` to the unit of this `p`-adic element.
        """
        raise NotImplementedError

    def __lshift__(pAdicTemplateElement self, shift):
        """
        Multiply ``self`` by ``pi^shift``.

        If ``shift`` is negative and this element does not lie in a field,
        digits may be truncated.  See ``__rshift__`` for details.

        EXAMPLES:

        We create some `p`-adic rings::

            sage: R = Zp(5, 20, 'capped-abs'); a = R(1000); a
            3*5^3 + 5^4 + O(5^20)
            sage: S = Zp(5); b = S(1000); b
            3*5^3 + 5^4 + O(5^23)

        Shifting to the right is the same as dividing by a power of
        the uniformizer `p` of the `p`-adic ring.::

            sage: a >> 1
            3*5^2 + 5^3 + O(5^19)
            sage: b >> 1
            3*5^2 + 5^3 + O(5^22)

        Shifting to the left is the same as multiplying by a power of
        `p`::

            sage: a << 2
            3*5^5 + 5^6 + O(5^20)
            sage: a*5^2
            3*5^5 + 5^6 + O(5^20)
            sage: b << 2
            3*5^5 + 5^6 + O(5^25)
            sage: b*5^2
            3*5^5 + 5^6 + O(5^25)

        Shifting by a negative integer to the left is the same as
        right shifting by the absolute value::

            sage: a << -3
            3 + 5 + O(5^17)
            sage: a >> 3
            3 + 5 + O(5^17)
            sage: b << -3
            3 + 5 + O(5^20)
            sage: b >> 3
            3 + 5 + O(5^20)
        """
        # TODO: move this up the hierarchy, perhaps this should go all the way to element?
        # The "verify that shift is an integer" part could be shared
        cdef long s
        if isinstance(shift, int):
            s = PyLong_AsLong(shift)
        else:
            if not isinstance(shift, Integer):
                shift = Integer(shift)
            if mpz_fits_slong_p((<Integer>shift).value) == 0:
                raise ValueError("valuation overflow")
            s = mpz_get_si((<Integer>shift).value)
        check_ordp(s)
        return self._lshift_c(s)

    cdef pAdicTemplateElement _lshift_c(self, long shift):
        raise NotImplementedError

    def __rshift__(pAdicTemplateElement self, shift):
        """
        Divide by ``p^shift``, and truncate (if the parent is not a field).

        EXAMPLES::

            sage: R = Zp(997, 7, 'capped-abs'); a = R(123456878908); a
            964*997 + 572*997^2 + 124*997^3 + O(997^7)
            sage: S = Zp(5); K = Qp(5); b = S(17); b
            2 + 3*5 + O(5^20)

        Shifting to the right divides by a power of `p`, but drops
        terms with negative valuation::

            sage: a >> 3
            124 + O(997^4)
            sage: b >> 1
            3 + O(5^19)
            sage: b >> 40
            O(5^0)

        If the parent is a field no truncation is performed::

            sage: K(17) >> 1
            2*5^-1 + 3 + O(5^19)

        A negative shift multiplies by that power of `p`::

            sage: a >> -3
            964*997^4 + 572*997^5 + 124*997^6 + O(997^7)
            sage: K(17) >> -5
            2*5^5 + 3*5^6 + O(5^25)
        """
        cdef long s
        if isinstance(shift, int):
            s = PyLong_AsLong(shift)
        else:
            if not isinstance(shift, Integer):
                shift = Integer(shift)
            if mpz_fits_slong_p((<Integer>shift).value) == 0:
                raise ValueError("valuation overflow")
            s = mpz_get_si((<Integer>shift).value)
        check_ordp(s)
        return self._rshift_c(s)

    cdef pAdicTemplateElement _rshift_c(self, long shift):
        """
        Divide by ``p^shift`` and truncate (if the parent is not a field).
        """
        raise NotImplementedError

    cdef int check_preccap(self) except -1:
        """
        Check that this element doesn't have precision higher than allowed by
        the precision cap.
        """
        raise NotImplementedError

    def lift_to_precision(self, absprec=None):
        r"""
        Return another element of the same parent with absolute precision at
        least ``absprec``, congruent to this `p`-adic element modulo the
        precision of this element.

        INPUT:

        - ``absprec`` -- integer or ``None`` (default: ``None``); the
          absolute precision of the result. If ``None``, lifts to the maximum
          precision allowed.

        .. NOTE::

            If setting ``absprec`` that high would violate the precision cap,
            raises a precision error.  Note that the new digits will not
            necessarily be zero.

        EXAMPLES::

            sage: R = ZpCA(17)
            sage: R(-1,2).lift_to_precision(10)
            16 + 16*17 + O(17^10)
            sage: R(1,15).lift_to_precision(10)
            1 + O(17^15)
            sage: R(1,15).lift_to_precision(30)
            Traceback (most recent call last):
            ...
            PrecisionError: precision higher than allowed by the precision cap
            sage: R(-1,2).lift_to_precision().precision_absolute() == R.precision_cap()
            True

            sage: R = Zp(5); c = R(17,3); c.lift_to_precision(8)
            2 + 3*5 + O(5^8)
            sage: c.lift_to_precision().precision_relative() == R.precision_cap()
            True

        Fixed modulus elements don't raise errors::

            sage: R = ZpFM(5); a = R(5); a.lift_to_precision(7)
            5
            sage: a.lift_to_precision(10000)
            5
        """
        if absprec is None:
            absprec = maxordp
        if not isinstance(absprec, Integer):
            absprec = Integer(absprec)
        if mpz_fits_slong_p((<Integer>absprec).value) == 0:
            raise PrecisionError("precision higher than allowed by the precision cap")
        ans = self.lift_to_precision_c(mpz_get_si((<Integer>absprec).value))
        ans.check_preccap()
        return ans

    cdef pAdicTemplateElement lift_to_precision_c(self, long absprec):
        """
        Lift this element to another with precision at least ``absprec``.
        """
        raise NotImplementedError

    def expansion(self, n = None, lift_mode = 'simple', start_val = None):
        r"""
        Return the coefficients in a `\pi`-adic expansion.
        If this is a field element, start at
        `\pi^{\mbox{valuation}}`, if a ring element at `\pi^0`.

        For each lift mode, this function returns a list of `a_i` so
        that this element can be expressed as

        .. MATH::

            \pi^v \cdot \sum_{i=0}^\infty a_i \pi^i,

        where `v` is the valuation of this element when the parent is
        a field, and `v = 0` otherwise.

        Different lift modes affect the choice of `a_i`.  When
        ``lift_mode`` is ``'simple'``, the resulting `a_i` will be
        nonnegative: if the residue field is `\GF{p}` then they
        will be integers with `0 \le a_i < p`; otherwise they will be
        a list of integers in the same range giving the coefficients
        of a polynomial in the indeterminant representing the maximal
        unramified subextension.

        Choosing ``lift_mode`` as ``'smallest'`` is similar to
        ``'simple'``, but uses a balanced representation `-p/2 < a_i
        \le p/2`.

        Finally, setting ``lift_mode = 'teichmuller'`` will yield
        Teichmuller representatives for the `a_i`: `a_i^q = a_i`.  In
        this case the `a_i` will lie in the ring of integers of the
        maximal unramified subextension of the parent of this element.

        INPUT:

        - ``n`` -- integer (default: ``None``); if given, returns the
          corresponding entry in the expansion. Can also accept a slice (see
          :meth:`slice`).

        - ``lift_mode`` -- ``'simple'``, ``'smallest'`` or
          ``'teichmuller'`` (default: ``'simple'``)

        - ``start_val`` -- start at this valuation rather than the
          default (`0` or the valuation of this element)

        OUTPUT:

        - If ``n`` is ``None``, an iterable giving a `\pi`-adic expansion of this
          element.  For base elements the contents will be integers if
          ``lift_mode`` is ``'simple'`` or ``'smallest'``, and
          elements of ``self.parent()`` if ``lift_mode`` is
          ``'teichmuller'``.

        - If ``n`` is an integer, the coefficient of `\pi^n` in the
          `\pi`-adic expansion of this element.

        .. NOTE::

            Use slice operators to get a particular range.

        EXAMPLES::

            sage: R = Zp(7,6); a = R(12837162817); a
            3 + 4*7 + 4*7^2 + 4*7^4 + O(7^6)
            sage: E = a.expansion(); E
            7-adic expansion of 3 + 4*7 + 4*7^2 + 4*7^4 + O(7^6)
            sage: list(E)
            [3, 4, 4, 0, 4, 0]
            sage: sum([c * 7^i for i, c in enumerate(E)]) == a
            True
            sage: E = a.expansion(lift_mode='smallest'); E
            7-adic expansion of 3 + 4*7 + 4*7^2 + 4*7^4 + O(7^6) (balanced)
            sage: list(E)
            [3, -3, -2, 1, -3, 1]
            sage: sum([c * 7^i for i, c in enumerate(E)]) == a
            True
            sage: E = a.expansion(lift_mode='teichmuller'); E
            7-adic expansion of 3 + 4*7 + 4*7^2 + 4*7^4 + O(7^6) (teichmuller)
            sage: list(E)
            [3 + 4*7 + 6*7^2 + 3*7^3 + 2*7^5 + O(7^6),
            0,
            5 + 2*7 + 3*7^3 + O(7^4),
            1 + O(7^3),
            3 + 4*7 + O(7^2),
            5 + O(7)]
            sage: sum(c * 7^i for i, c in enumerate(E))
            3 + 4*7 + 4*7^2 + 4*7^4 + O(7^6)

        If the element has positive valuation then the list will start
        with some zeros::

            sage: a = R(7^3 * 17)
            sage: E = a.expansion(); E
            7-adic expansion of 3*7^3 + 2*7^4 + O(7^9)
            sage: list(E)
            [0, 0, 0, 3, 2, 0, 0, 0, 0]

        The expansion of 0 is truncated::

            sage: E = R(0, 7).expansion(); E
            7-adic expansion of O(7^7)
            sage: len(E)
            0
            sage: list(E)
            []

        In fields, on the other hand, the expansion starts at the valuation::

            sage: R = Qp(7,4); a = R(6*7+7**2); E = a.expansion(); E
            7-adic expansion of 6*7 + 7^2 + O(7^5)
            sage: list(E)
            [6, 1, 0, 0]
            sage: list(a.expansion(lift_mode='smallest'))
            [-1, 2, 0, 0]
            sage: list(a.expansion(lift_mode='teichmuller'))
            [6 + 6*7 + 6*7^2 + 6*7^3 + O(7^4),
            2 + 4*7 + 6*7^2 + O(7^3),
            3 + 4*7 + O(7^2),
            3 + O(7)]

        You can ask for a specific entry in the expansion::

            sage: a.expansion(1)
            6
            sage: a.expansion(1, lift_mode='smallest')
            -1
            sage: a.expansion(2, lift_mode='teichmuller')
            2 + 4*7 + 6*7^2 + O(7^3)

        TESTS:

        Check to see that :issue:`10292` is resolved::

            sage: # needs sage.schemes
            sage: E = EllipticCurve('37a')
            sage: R = E.padic_regulator(7)
            sage: len(R.expansion())
            19
        """
        if isinstance(n, slice):
            return self.slice(n.start, n.stop, n.step, lift_mode=lift_mode)

        cdef long shift, prec, val
        cdef expansion_mode mode
        prec = self.precision_relative()
        val = self.valuation_c()
        if prec == 0:
            shift = 0
        elif start_val is not None:
            if n is not None:
                raise ValueError("n and start_val are incompatible options")
            shift = val - start_val
        elif self.prime_pow.in_field:
            shift = 0
        else:
            shift = val

        if lift_mode == 'simple':
            mode = simple_mode
        elif lift_mode == 'smallest':
            mode = smallest_mode
        elif lift_mode == 'teichmuller':
            mode = teichmuller_mode
        else:
            raise ValueError("unknown lift_mode")

        cdef ExpansionIterable expansion = ExpansionIterable(self, prec, shift, mode)
        if n is None:
            return expansion
        else:
            if n < val:
                return _zero(mode, expansion.teich_ring)
            elif self.prime_pow.in_field:
                return expansion[n - val]
            else:
                return expansion[n]

    def teichmuller_expansion(self, n = None):
        r"""
        Return an iterator over coefficients `a_0, a_1, \dots, a_n` such that

        - `a_i^q = a_i`, where `q` is the cardinality of the residue field,

        - this element can be expressed as

        .. MATH::

            \pi^v \cdot \sum_{i=0}^\infty a_i \pi^i

        where `v` is the valuation of this element when the parent is
        a field, and `v = 0` otherwise.

        - if `a_i \ne 0`, the precision of `a_i` is `i` less
          than the precision of this element (relative in the case that
          the parent is a field, absolute otherwise)

        .. NOTE::

            The coefficients will lie in the ring of integers of the
            maximal unramified subextension.

        INPUT:

        - ``n`` -- integer (default: ``None``); if given, returns the
          coefficient of `\pi^n` in the expansion

        EXAMPLES:

        For fields, the expansion starts at the valuation::

            sage: R = Qp(5,5); list(R(70).teichmuller_expansion())
            [4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5),
            3 + 3*5 + 2*5^2 + 3*5^3 + O(5^4),
            2 + 5 + 2*5^2 + O(5^3),
            1 + O(5^2),
            4 + O(5)]

        But if you specify ``n``, you get the coefficient of `\pi^n`::

            sage: R(70).teichmuller_expansion(2)
            3 + 3*5 + 2*5^2 + 3*5^3 + O(5^4)
        """
        return self.expansion(n, lift_mode='teichmuller')

    def _ext_p_list(self, pos):
        """
        Return the `p`-adic expansion of the unit part.  Used in printing.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<a> = Qq(125)
            sage: b = a^2 + 5*a + 1
            sage: b._ext_p_list(True)
            [[1, 0, 1], [0, 1]]
        """
        if pos:
            return trim_zeros(list(self.unit_part().expansion(lift_mode='simple')))
        else:
            return trim_zeros(list(self.unit_part().expansion(lift_mode='smallest')))

    cpdef pAdicTemplateElement unit_part(self):
        r"""
        Return the unit part of this element.

        This is the `p`-adic element `u` in the same ring so that this
        element is `\pi^v u`, where `\pi` is a uniformizer and `v` is
        the valuation of this element.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<a> = Zq(125)
            sage: (5*a).unit_part()
            a + O(5^20)
        """
        raise NotImplementedError

    cpdef bint _is_base_elt(self, p) except -1:
        """
        Return ``True`` if this element is an element of Zp or Qp (rather than
        an extension).

        INPUT:

        - ``p`` -- a prime, which is compared with the parent of this element

        EXAMPLES::

            sage: a = Zp(5)(3); a._is_base_elt(5)
            True
            sage: a._is_base_elt(17)
            False
        """
        return self.prime_pow.prime == p and self.prime_pow.deg == 1

    def _prime_pow(self):
        """
        Provides access to this element's ``prime_pow``.

        EXAMPLES::

            sage: R = ZpCR(5,5)
            sage: R(1)._prime_pow()
            PowComputer for 5
        """
        return self.prime_pow

    def residue(self, absprec=1, field=None, check_prec=True):
        r"""
        Reduce this element modulo `p^\mathrm{absprec}`.

        INPUT:

        - ``absprec`` -- ``0`` or ``1``

        - ``field`` -- boolean (default: ``None``); for precision 1, whether to return
          an element of the residue field or a residue ring.  Currently unused.

        - ``check_prec`` -- boolean (default: ``True``); whether to raise an error if this
          element has insufficient precision to determine the reduction.  Errors are never
          raised for fixed-mod or floating-point types.

        OUTPUT:

        This element reduced modulo `p^\mathrm{absprec}` as an element of the
        residue field or the null ring.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<a> = Zq(27, 4)
            sage: (3 + 3*a).residue()
            0
            sage: (a + 1).residue()
            a0 + 1

        TESTS::

            sage: # needs sage.libs.ntl
            sage: a.residue(0)
            0
            sage: a.residue(2)
            Traceback (most recent call last):
            ...
            NotImplementedError: reduction modulo p^n with n>1
            sage: a.residue(10)
            Traceback (most recent call last):
            ...
            PrecisionError: insufficient precision to reduce modulo p^10
            sage: a.residue(10, check_prec=False)
            Traceback (most recent call last):
            ...
            NotImplementedError: reduction modulo p^n with n>1

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqCA(27, 4)
            sage: (3 + 3*a).residue()
            0
            sage: (a + 1).residue()
            a0 + 1

            sage: # needs sage.libs.ntl
            sage: R.<a> = Qq(27, 4)
            sage: (3 + 3*a).residue()
            0
            sage: (a + 1).residue()
            a0 + 1
            sage: (a/3).residue()
            Traceback (most recent call last):
            ...
            ValueError: element must have nonnegative valuation in order to compute residue
        """
        if absprec < 0:
            raise ValueError("cannot reduce modulo a negative power of the uniformizer")
        if self.valuation() < 0:
            raise ValueError("element must have nonnegative valuation in order to compute residue")
        R = self.parent()
        if check_prec and (R.is_fixed_mod() or R.is_floating_point()):
            check_prec = False
        if check_prec and absprec > self.precision_absolute():
            raise PrecisionError("insufficient precision to reduce modulo p^%s" % absprec)
        if field and absprec != 1:
            raise ValueError("field keyword may only be set at precision 1")
        if absprec == 0:
            from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
            return IntegerModRing(1).zero()
        elif absprec == 1:
            parent = R.residue_field()
            if self.valuation() > 0:
                return parent.zero()
            return parent(self.expansion(0))
        else:
            raise NotImplementedError("reduction modulo p^n with n>1")

cdef Integer exact_pow_helper(long *ansrelprec, long relprec, _right, PowComputer_ prime_pow):
    """
    This function is used by exponentiation in both ``CR_template.pxi``
    and ``CA_template.pxi`` to determine the extra precision gained from
    an exponent of positive valuation.  See ``__pow__`` there and in
    ``padic_ZZ_pX_CR_element.pyx`` for more details on this phenomenon.

    INPUT:

    - ``ansrelprec`` -- (return value) the relative precision of the answer

    - ``relprec`` -- positive integer: the relative precision of the base

    - ``_right`` -- the exponent, nonzero

    - ``prime_pow`` -- the Powcomputer for the ring

    OUTPUT: an Integer congruent to the given exponent
    """
    cdef Integer right, p = prime_pow.prime
    cdef long exp_val
    cdef bint isbase
    if isinstance(_right, int):
        _right = Integer(_right)
    if isinstance(_right, Integer):
        right = <Integer> _right
        # Be careful: prime_pow.e is assumed to be the absolute index of ramification!
        exp_val = mpz_get_si((<Integer>(right.valuation(p) * prime_pow.e)).value)
    elif isinstance(_right, Rational):
        raise NotImplementedError
    ansrelprec[0] = relprec + exp_val
    # Over Z_2 or Q_2, the square of an odd number is congruent to 1 mod 8
    if exp_val > 0 and prime_pow.deg == 1 and mpz_cmp_ui(p.value, 2) == 0 and relprec == 1:
        ansrelprec[0] += 1

    return right

cdef long padic_pow_helper(celement result, celement base, long base_val, long base_relprec,
                           celement right_unit, long right_val, long right_relprec, PowComputer_ prime_pow) except -1:
    """
    INPUT:

    - ``result`` -- the result of exponentiation

    - ``base`` -- a celement; the base of the exponentiation

    - ``base_val`` -- a long; used to check that the base is a unit

    - ``base_relprec`` -- positive integer; the relative precision
      of the base

    - ``right_unit`` -- the unit part of the exponent

    - ``right_val`` -- the valuation of the exponent

    - ``right_relprec`` -- the relative precision of the exponent

    - ``prime_pow`` -- the Powcomputer for the ring

    OUTPUT: the precision of the result

    EXAMPLES::

        sage: R = Zp(17,print_mode='digits')
        sage: a = R(9283732, 6); b = R(17^3*237, 7)
        sage: str(a)
        '...692AAF'
        sage: str(a^b)  # indirect doctest
        '...55GA0001'
        sage: str((a // R.teichmuller(15))^b)
        '...55GA0001'
        sage: str((a.log()*b).exp())
        '...55GA0001'
    """
    if base_val != 0:
        raise ValueError("in order to raise to a p-adic exponent, base must be a unit")
    # TODO NOTE:  this function needs to be updated for extension elements.
    cdef long loga_val, loga_aprec, bloga_val, bloga_aprec
    cdef Integer expcheck, right
    cteichmuller(prime_pow.powhelper_oneunit, base, base_relprec, prime_pow)
    cdivunit(prime_pow.powhelper_oneunit, base, prime_pow.powhelper_oneunit, base_relprec, prime_pow)
    csetone(prime_pow.powhelper_teichdiff, prime_pow)
    csub(prime_pow.powhelper_teichdiff, prime_pow.powhelper_oneunit, prime_pow.powhelper_teichdiff, base_relprec, prime_pow)
    # For extension elements in ramified extensions, the computation of the
    # valuation and precision of log(a) is more complicated)
    loga_val = cvaluation(prime_pow.powhelper_teichdiff, base_relprec, prime_pow)
    loga_aprec = base_relprec
    # valuation of b*log(a)
    bloga_val = loga_val + right_val
    bloga_aprec = bloga_val + min(right_relprec, loga_aprec - loga_val)
    if bloga_aprec > prime_pow.ram_prec_cap:
        bloga_aprec = prime_pow.ram_prec_cap
    expcheck = PY_NEW(Integer)
    mpz_sub_ui(expcheck.value, prime_pow.prime.value, 1)
    mpz_mul_si(expcheck.value, expcheck.value, bloga_val)
    if mpz_cmp_ui(expcheck.value, prime_pow.e) <= 0:
        raise ValueError("exponential does not converge")
    right = PY_NEW(Integer)
    try:
        cconv_mpz_t_out(right.value, right_unit, right_val, right_relprec, prime_pow)
    except ValueError:
        # Here we need to use the exp(b log(a)) definition,
        # since we can't convert the exponent to an integer
        raise NotImplementedError("exponents with negative valuation not yet supported")
    # For extension elements in ramified extensions
    # the following precision might need to be changed.
    cpow(result, prime_pow.powhelper_oneunit, right.value, bloga_aprec, prime_pow)
    return bloga_aprec

cdef _zero(expansion_mode mode, teich_ring):
    """
    Return an appropriate zero for a given expansion mode.

    INPUT:

    - ``mode`` -- either ``simple_mode`` or ``smallest_mode`` or ``teichmuller_mode``
    - ``teich_ring`` -- the integer ring of the maximal unramified subextension
      of the parent.  Only used in ``teichmuller_mode``.
    """
    if mode == teichmuller_mode:
        return teich_ring(0)
    else:
        return _expansion_zero

cdef class ExpansionIter():
    """
    An iterator over a `p`-adic expansion.

    This class should not be instantiated directly, but instead
    using :meth:`expansion`.

    INPUT:

    - ``elt`` -- the `p`-adic element
    - ``prec`` -- the number of terms to be emitted
    - ``mode`` -- either ``simple_mode``, ``smallest_mode`` or ``teichmuller_mode``

    EXAMPLES::

        sage: E = Zp(5,4)(373).expansion()
        sage: I = iter(E)  # indirect doctest
        sage: type(I)
        <class 'sage.rings.padics.padic_capped_relative_element.ExpansionIter'>
    """
    cdef pAdicTemplateElement elt
    cdef celement tmp
    cdef celement curvalue
    cdef long curpower
    cdef bint tracks_prec
    cdef expansion_mode mode
    cdef object teich_ring

    def __cinit__(self, pAdicTemplateElement elt, long prec, expansion_mode mode):
        """
        Allocate memory for the iterator.

        TESTS::

            sage: E = Zp(5,4)(373).expansion()
            sage: I = iter(E)
            sage: next(I)
            3
        """
        self.elt = elt
        self.curpower = prec
        self.mode = mode
        if mode == teichmuller_mode:
            R = elt.parent()
            self.tracks_prec = R.is_capped_relative() or R.is_capped_absolute()
            self.teich_ring = R.maximal_unramified_subextension().integer_ring()
        IF CELEMENT_IS_PY_OBJECT:
            polyt = type(elt.prime_pow.modulus)
            self.tmp = <celement>polyt.__new__(polyt)
            self.curvalue = <celement>polyt.__new__(polyt)
        cconstruct(self.tmp, elt.prime_pow)
        cconstruct(self.curvalue, elt.prime_pow)
        elt._get_unit(self.curvalue)

    def __dealloc__(self):
        """
        Deallocate memory for the iterator.

        TESTS::

            sage: E = Zp(5,4)(373).expansion()
            sage: I = iter(E)
            sage: del I
        """
        cdestruct(self.tmp, self.elt.prime_pow)
        cdestruct(self.curvalue, self.elt.prime_pow)

    def __iter__(self):
        """
        Characteristic property of an iterator: ``__iter__`` returns itself.

        TESTS::

            sage: E = Zp(5,4)(373).expansion()
            sage: I = iter(E)
            sage: I is iter(I)
            True
        """
        return self

    def __len__(self):
        """
        Return the number of terms that will be emitted.

        TESTS::

            sage: E = Zp(5,4)(373).expansion()
            sage: I = iter(E)
            sage: len(I)
            4
            sage: c = next(I); len(I)
            3
        """
        return self.curpower

    def __next__(self):
        """
        Provide the next coefficient in the `p`-adic expansion.

        EXAMPLES::

            sage: E = Zp(5,4)(373).expansion()
            sage: I = iter(E)
            sage: next(I)
            3
            sage: next(I), next(I), next(I)
            (4, 4, 2)
        """
        if self.curpower <= 0:
            raise StopIteration
        self.curpower -= 1
        cdef pAdicTemplateElement ans
        cdef PowComputer_ pp = self.elt.prime_pow
        cdef long prec
        if ciszero(self.curvalue, pp):
            return _zero(self.mode, self.teich_ring)
        if self.mode == teichmuller_mode:
            prec = self.curpower+1 if self.tracks_prec else pp.ram_prec_cap
            cteichmuller(self.tmp, self.curvalue, prec, pp)
            if ciszero(self.tmp, pp):
                cshift_notrunc(self.curvalue, self.curvalue, -1, prec-1, pp, True)
                return _zero(teichmuller_mode, self.teich_ring)
            else:
                csub(self.curvalue, self.curvalue, self.tmp, prec, pp)
                cshift_notrunc(self.curvalue, self.curvalue, -1, prec-1, pp, True)
                return self.teich_ring(self.elt._new_with_value(self.tmp, prec))
        else:
            return cexpansion_next(self.curvalue, self.mode, self.curpower, pp)

cdef class ExpansionIterable():
    r"""
    An iterable storing a `p`-adic expansion of an element.

    This class should not be instantiated directly, but instead using :meth:`expansion`.

    INPUT:

    - ``elt`` -- the `p`-adic element
    - ``prec`` -- the number of terms to be emitted
    - ``val_shift`` -- how many zeros to add at the beginning of the expansion,
      or the number of initial terms to truncate (if negative)
    - ``mode`` -- one of the following:

      * ``'simple_mode'``
      * ``'smallest_mode'``
      * ``'teichmuller_mode'``

    EXAMPLES::

        sage: E = Zp(5,4)(373).expansion()  # indirect doctest
        sage: type(E)
        <class 'sage.rings.padics.padic_capped_relative_element.ExpansionIterable'>
    """
    cdef pAdicTemplateElement elt
    cdef celement tmp
    cdef long prec
    cdef long val_shift
    cdef expansion_mode mode
    cdef object teich_ring

    def __cinit__(self, pAdicTemplateElement elt, long prec, long val_shift, expansion_mode mode):
        """
        Allocate memory for the iteratable.

        TESTS::

            sage: Zp(5,4)(373).expansion()
            5-adic expansion of 3 + 4*5 + 4*5^2 + 2*5^3 + O(5^4)
        """
        self.elt = elt
        IF CELEMENT_IS_PY_OBJECT:
            polyt = type(elt.prime_pow.modulus)
            self.tmp = <celement>polyt.__new__(polyt)
        cconstruct(self.tmp, elt.prime_pow)
        self.prec = prec
        self.val_shift = val_shift
        self.mode = mode
        if mode == teichmuller_mode:
            self.teich_ring = elt.parent().maximal_unramified_subextension().integer_ring()

    def __dealloc__(self):
        """
        Deallocate memory for the iteratable.

        TESTS::

            sage: E = Zp(5,4)(373).expansion()
            sage: del E
        """
        cdestruct(self.tmp, self.elt.prime_pow)

    def __iter__(self):
        """
        Return an iterator, based on a corresponding :class:`ExpansionIter`.

        If ``val_shift`` is positive, will first emit that many zeros
        (of the appropriate type: ``[]`` instead when the inertia degree
        is larger than one.

        If ``val_shift`` is negative, will truncate that many terms at
        the start of the expansion.

        EXAMPLES::

            sage: E = Zp(5,4)(373).expansion()
            sage: type(iter(E))
            <class 'sage.rings.padics.padic_capped_relative_element.ExpansionIter'>
            sage: E = Zp(5,4)(373).expansion(start_val=-1)
            sage: type(iter(E))
            <class 'itertools.chain'>
            sage: E = Zp(5,4)(373).expansion(start_val=1)
            sage: type(iter(E))
            <class 'itertools.islice'>
        """
        cdef ExpansionIter expansion = ExpansionIter(self.elt, self.prec, self.mode)
        if self.val_shift == 0:
            return expansion
        elif self.val_shift < 0:
            return itertools.islice(expansion, -self.val_shift, None)
        else:
            return itertools.chain(itertools.repeat(_zero(self.mode, self.teich_ring), self.val_shift), expansion)

    def __len__(self):
        """
        Return the number of terms that will be emitted.

        TESTS::

            sage: len(Zp(5,4)(373).expansion())
            4
            sage: len(Zp(5,4)(373).expansion(start_val=-1))
            5
            sage: len(Zp(5,4)(373).expansion(start_val=1))
            3
            sage: len(Zp(5,4)(0).expansion())
            0
        """
        return self.prec + self.val_shift

    def __getitem__(self, n):
        """
        Return the ``n``-th entry in the expansion.

        Negative indices are not allowed.

        EXAMPLES::

            sage: E = Zp(5,4)(373).expansion()
            sage: E[0]
            3
            sage: E[3]
            2
            sage: list(E[::2])
            [3, 4]
            sage: a = E[-1]
            Traceback (most recent call last):
            ...
            ValueError: negative indices not supported
            sage: Zp(5,4)(373).expansion(lift_mode='smallest')[3]
            -2
        """
        if isinstance(n, slice):
            start = int(n.start) if n.start is not None else None
            stop = int(n.stop) if n.stop is not None else None
            step = int(n.step) if n.step is not None else None
            return itertools.islice(iter(self), start, stop, step)
        cdef long m = n - self.val_shift
        cdef celement value
        if n < 0:
            raise ValueError("negative indices not supported")
        elif m < 0:
            return _zero(self.mode, self.teich_ring)
        elif m >= self.prec:
            raise PrecisionError
        elif self.mode == simple_mode:
            self.elt._get_unit(self.tmp)
            return cexpansion_getitem(self.tmp, m, self.elt.prime_pow)
        else:
            expansion = ExpansionIter(self.elt, self.prec, self.mode)
            # We do this in a naive way, though it should be feasible
            # to implement something better for smallest_mode
            return next(itertools.islice(expansion, m, m+1))

    def __repr__(self):
        """
        String representation.

        EXAMPLES::

            sage: Zp(5,4)(373).expansion()
            5-adic expansion of 3 + 4*5 + 4*5^2 + 2*5^3 + O(5^4)
            sage: Zp(5,4)(373).expansion(lift_mode='smallest')
            5-adic expansion of 3 + 4*5 + 4*5^2 + 2*5^3 + O(5^4) (balanced)
            sage: Zp(5,4)(373).expansion(lift_mode='teichmuller')
            5-adic expansion of 3 + 4*5 + 4*5^2 + 2*5^3 + O(5^4) (teichmuller)
        """
        if self.mode == simple_mode:
            modestr = ""
        elif self.mode == smallest_mode:
            modestr = " (balanced)"
        else:
            modestr = " (teichmuller)"
        p = self.elt.prime_pow.prime
        return "%s-adic expansion of %s%s" % (p, self.elt, modestr)
