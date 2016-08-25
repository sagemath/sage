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

#*****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.int cimport *

from sage.libs.gmp.all cimport *
import sage.rings.finite_rings.integer_mod
from sage.libs.pari.types cimport *
from sage.libs.pari.gen cimport gen as pari_gen
from sage.libs.pari.pari_instance cimport INT_to_mpz
from sage.rings.padics.common_conversion cimport get_ordp, get_preccap
from sage.rings.integer cimport Integer
from sage.rings.infinity import infinity
from sage.rings.rational import Rational
from sage.rings.padics.precision_error import PrecisionError
from sage.structure.element import canonical_coercion

cdef long maxordp = (1L << (sizeof(long) * 8 - 2)) - 1
cdef long minusmaxordp = -maxordp

cdef inline int check_ordp(long ordp) except -1:
    """
    Checks for overflow after addition or subtraction of valuations.

    There is another variant, :meth:`check_ordp_mpz`, for ``mpz_t`` input.

    If overflow is detected, raises a OverflowError.
    """
    if ordp >= maxordp or ordp <= minusmaxordp:
        raise OverflowError("valuation overflow")

cdef class pAdicTemplateElement(pAdicGenericElement):
    """
    A class for common functionality among the `p`-adic template classes.

    INPUT:

    - ``parent`` -- a local ring or field

    - ``x`` -- data defining this element.  Various types are supported,
      including ints, Integers, Rationals, PARI p-adics, integers mod `p^k`
      and other Sage p-adics.

    - ``absprec`` -- a cap on the absolute precision of this element

    - ``relprec`` -- a cap on the relative precision of this element

    EXAMPLES::

        sage: Zp(17)(17^3, 8, 4)
        17^3 + O(17^7)
    """
    def __init__(self, parent, x, absprec=infinity, relprec=infinity):
        """
        Initialization.

        .. NOTE:

            This initialization function is not called for Integers
            and Rationals since a conversion morphism has been
            implemented.  It is, however, used for python ints and longs.

        EXAMPLES::

            sage: a = Zp(5)(1/2,3); a
            3 + 2*5 + 2*5^2 + O(5^3)
            sage: type(a)
            <type 'sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement'>
            sage: TestSuite(a).run()

        TESTS::

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

        """
        self.prime_pow = <PowComputer_?>parent.prime_pow
        pAdicGenericElement.__init__(self, parent)
        cdef long val, xprec
        cdef GEN pari_tmp
        if isinstance(x, (int, long)):
            x = Integer(x)
        elif isinstance(x, pari_gen):
            pari_tmp = (<pari_gen>x).g
            if typ(pari_tmp) == t_INT:
                x = PY_NEW(Integer)
                INT_to_mpz((<Integer>x).value, pari_tmp)
            elif typ(pari_tmp) == t_FRAC:
                x = Rational(x)
        elif isinstance(x, pAdicGenericElement):
            if not ((<pAdicGenericElement>x)._is_base_elt(self.prime_pow.prime) or x.parent() is self.parent()):
                raise NotImplementedError("conversion between padic extensions not implemented")
        elif sage.rings.finite_rings.integer_mod.is_IntegerMod(x):
            if not Integer(self.prime_pow.prime).divides(x.parent().order()):
                raise TypeError("p does not divide modulus %s"%x.parent().order())
        elif sage.rings.finite_rings.element_base.is_FiniteFieldElement(x):
            k = self.parent().residue_field()
            if not k.has_coerce_map_from(x.parent()):
                raise NotImplementedError("conversion from finite fields which do not embed into the residue field not implemented.")

            x = k(x)
            if not k.is_prime_field():
                x = [k.prime_subfield()(c) for c in x.polynomial().list()]
                x = x + [k.prime_subfield().zero()] * (k.degree() - len(x))
        elif isinstance(x, (Integer, Rational, list, tuple)):
            pass
        elif sage.rings.polynomial.polynomial_element.is_Polynomial(x) and x.variable_name() == self.parent().variable_name():
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
        Sets this element from given data computed in :meth:`__init__`.

        INPUT:

        - ``x`` -- an int, Integer, Rational, PARI p-adic, integer mod `p^k` or Sage p-adic

        - ``val`` -- a long, the valuation of ``x``

        - ``xprec`` -- a long, the cap on the absolute precision imposed by ``x``

        - ``absprec`` -- an integer or infinity, a bound on the absolute precision

        - ``relprec`` -- an integer or infinity, a bound on the relative precision

        """
        raise NotImplementedError

    def __lshift__(pAdicTemplateElement self, shift):
        """
        Multiplies ``self`` by ``pi^shift``.

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
            s = PyInt_AS_LONG(shift)
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
        Divides by ``p^shift``, and truncates (if the parent is not a field).

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
            s = PyInt_AS_LONG(shift)
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
        Divides by ``p^shift`` and truncates (if the parent is not a field).
        """
        raise NotImplementedError

    cdef int check_preccap(self) except -1:
        """
        Checks that this element doesn't have precision higher than allowed by
        the precision cap.
        """
        raise NotImplementedError

    def lift_to_precision(self, absprec=None):
        """
        Returns another element of the same parent with absolute precision at
        least ``absprec``, congruent to this `p`-adic element modulo the
        precision of this element.

        INPUT:

        - ``absprec`` -- an integer or ``None`` (default: ``None``), the
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
            PrecisionError: Precision higher than allowed by the precision cap.
            sage: R(-1,2).lift_to_precision().precision_absolute() == R.precision_cap()
            True

            sage: R = Zp(5); c = R(17,3); c.lift_to_precision(8)
            2 + 3*5 + O(5^8)
            sage: c.lift_to_precision().precision_relative() == R.precision_cap()
            True

        Fixed modulus elements don't raise errors::

            sage: R = ZpFM(5); a = R(5); a.lift_to_precision(7)
            5 + O(5^20)
            sage: a.lift_to_precision(10000)
            5 + O(5^20)

        """
        if absprec is None:
            absprec = maxordp
        if not isinstance(absprec, Integer):
            absprec = Integer(absprec)
        if mpz_fits_slong_p((<Integer>absprec).value) == 0:
            raise PrecisionError("Precision higher than allowed by the precision cap")
        ans = self.lift_to_precision_c(mpz_get_si((<Integer>absprec).value))
        ans.check_preccap()
        return ans

    cdef pAdicTemplateElement lift_to_precision_c(self, long absprec):
        """
        Lifts this element to another with precision at least absprec.
        """
        raise NotImplementedError

    def padded_list(self, n, lift_mode = 'simple'):
        """
        Returns a list of coefficients of the uniformizer `\pi`
        starting with `\pi^0` up to `\pi^n` exclusive (padded with
        zeros if needed).

        For a field element of valuation `v`, starts at `\pi^v`
        instead.

        INPUT:

        - ``n`` - an integer

        - ``lift_mode`` - 'simple', 'smallest' or 'teichmuller'

        EXAMPLES::

            sage: R = Zp(7,4,'capped-abs'); a = R(2*7+7**2); a.padded_list(5)
            [0, 2, 1, 0, 0]
            sage: R = Zp(7,4,'fixed-mod'); a = R(2*7+7**2); a.padded_list(5)
            [0, 2, 1, 0, 0]

        For elements with positive valuation, this function will
        return a list with leading 0s if the parent is not a field::

            sage: R = Zp(7,3,'capped-rel'); a = R(2*7+7**2); a.padded_list(5)
            [0, 2, 1, 0, 0]
            sage: R = Qp(7,3); a = R(2*7+7**2); a.padded_list(5)
            [2, 1, 0, 0]
            sage: a.padded_list(3)
            [2, 1]
        """
        if lift_mode == 'simple' or lift_mode == 'smallest':
            # needs to be defined in the linkage file.
            zero = _list_zero
        elif lift_mode == 'teichmuller':
            zero = self.parent()(0,0)
        else:
            raise ValueError("%s not a recognized lift mode"%lift_mode)
        L = self.list(lift_mode)
        if self.prime_pow.in_field == 1:
            if self._is_exact_zero():
                n = 0
            else:
                n -= self.valuation()
        return L[:n] + [zero] * (n - len(L))

    def _ext_p_list(self, pos):
        if pos:
            return self.unit_part().list('simple')
        else:
            return self.unit_part().list('smallest')

    cpdef pAdicTemplateElement unit_part(self):
        """
        Returns the unit part of this element.

        This is the `p`-adic element `u` in the same ring so that this
        element is `\pi^v u`, where `\pi` is a uniformizer and `v` is
        the valuation of this element.
        """
        raise NotImplementedError

    cpdef bint _is_base_elt(self, p) except -1:
        """
        Return ``True`` if this element is an element of Zp or Qp (rather than
        an extension).

        INPUT:

        - ``p`` -- a prime, which is compared with the parent of this element.

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

    def residue(self, absprec=1):
        r"""
        Reduce this element modulo `p^\mathrm{absprec}`.

        INPUT:

        - ``absprec`` -- ``0`` or ``1``.

        OUTPUT:

        This element reduced modulo `p^\mathrm{absprec}` as an element of the
        residue field or the null ring.

        EXAMPLES::

            sage: R.<a> = ZqFM(27, 4)
            sage: (3 + 3*a).residue()
            0
            sage: (a + 1).residue()
            a0 + 1

        TESTS::

            sage: a.residue(0)
            0
            sage: a.residue(2)
            Traceback (most recent call last):
            ...
            NotImplementedError: reduction modulo p^n with n>1.
            sage: a.residue(10)
            Traceback (most recent call last):
            ...
            PrecisionError: insufficient precision to reduce modulo p^10.

            sage: R.<a> = ZqCA(27, 4)
            sage: (3 + 3*a).residue()
            0
            sage: (a + 1).residue()
            a0 + 1

            sage: R.<a> = Qq(27, 4)
            sage: (3 + 3*a).residue()
            0
            sage: (a + 1).residue()
            a0 + 1
            sage: (a/3).residue()
            Traceback (most recent call last):
            ...
            ValueError: element must have non-negative valuation in order to compute residue.
            
        """
        if absprec < 0:
            raise ValueError("cannot reduce modulo a negative power of the uniformizer.")
        if self.valuation() < 0:
            raise ValueError("element must have non-negative valuation in order to compute residue.")
        if absprec > self.precision_absolute():
            raise PrecisionError("insufficient precision to reduce modulo p^%s."%absprec)
        if absprec == 0:
            from sage.rings.all import IntegerModRing
            return IntegerModRing(1).zero()
        elif absprec == 1:
            parent = self.parent().residue_field()
            if self.valuation() > 0:
                return parent.zero()
            digits = self.padded_list(1)
            return parent(digits[0])
        else:
            raise NotImplementedError("reduction modulo p^n with n>1.")

cdef Integer exact_pow_helper(long *ansrelprec, long relprec, _right, PowComputer_ prime_pow):
    """
    This function is used by exponentiation in both CR_template.pxi
    and CA_template.pxi to determine the extra precision gained from
    an exponent of positive valuation.  See __pow__ there and in
    padic_ZZ_pX_CR_element.pyx for more details on this phenomenon.

    INPUT:

    - ``ansrelprec`` -- (return value) the relative precision of the answer

    - ``relprec`` -- a positive integer: the relative precision of the base

    - ``_right`` -- the exponent, nonzero

    - ``prime_pow`` -- the Powcomputer for the ring.

    OUTPUT:

    an Integer congruent to the given exponent

    """
    ####### NOTE:  this function needs to be updated for extension elements. #######
    cdef Integer right, p = prime_pow.prime
    cdef long exp_val
    cdef bint isbase
    if isinstance(_right, (int, long)):
        _right = Integer(_right)
    if isinstance(_right, Integer):
        right = <Integer> _right
        exp_val = mpz_get_si((<Integer>right.valuation(p)).value)
    elif isinstance(_right, Rational):
        raise NotImplementedError
    ansrelprec[0] = relprec + exp_val
    if exp_val > 0 and mpz_cmp_ui(p.value, 2) == 0 and relprec == 1:
        ansrelprec[0] += 1

    return right

cdef long padic_pow_helper(celement result, celement base, long base_val, long base_relprec,
                           celement right_unit, long right_val, long right_relprec, PowComputer_ prime_pow) except -1:
    """
    INPUT:

    - ``result`` -- the result of exponentiation.

    - ``base`` -- a celement, the base of the exponentiation.

    - ``base_val`` -- a long, used to check that the base is a unit

    - ``base_relprec`` -- a positive integer: the relative precision
      of the base.

    - ``right_unit`` -- the unit part of the exponent

    - ``right_val`` -- the valuation of the exponent

    - ``right_relprec`` -- the relative precision of the exponent

    - ``prime_pow`` -- the Powcomputer for the ring.

    OUTPUT:

    the precision of the result

    EXAMPLES::

        sage: R = Zp(17,print_mode='digits')
        sage: a = R(9283732, 6); b = R(17^3*237, 7)
        sage: str(a)
        '...692AAF'
        sage: str(a^b) # indirect doctest
        '...55GA0001'
        sage: str((a // R.teichmuller(15))^b)
        '...55GA0001'
        sage: str((a.log()*b).exp())
        '...55GA0001'
    """
    if base_val != 0:
        raise ValueError("in order to raise to a p-adic exponent, base must be a unit")
    ####### NOTE:  this function needs to be updated for extension elements. #######
    cdef celement oneunit, teichdiff
    cdef long loga_val, loga_aprec, bloga_val, bloga_aprec
    cdef Integer expcheck, right
    try:
        cconstruct(oneunit, prime_pow)
        cconstruct(teichdiff, prime_pow)
        cteichmuller(oneunit, base, base_relprec, prime_pow)
        cdivunit(oneunit, base, oneunit, base_relprec, prime_pow)
        csetone(teichdiff, prime_pow)
        csub(teichdiff, oneunit, teichdiff, base_relprec, prime_pow)
        ## For extension elements in ramified extensions, the computation of the
        ## valuation and precision of log(a) is more complicated)
        loga_val = cvaluation(teichdiff, base_relprec, prime_pow)
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
        ## For extension elements in ramified extensions
        ## the following precision might need to be changed.
        cpow(result, oneunit, right.value, bloga_aprec, prime_pow)
    finally:
        cdestruct(oneunit, prime_pow)
        cdestruct(teichdiff, prime_pow)
    return bloga_aprec
