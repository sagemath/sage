"""
`p`-Adic Capped Relative Elements

Elements of `p`-Adic Rings with Capped Relative Precision

AUTHORS:

- David Roe: initial version, rewriting to use templates (2012-3-1)
- Genya Zaytman: documentation
- David Harvey: doctests
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

include "sage/libs/linkages/padics/mpz.pxi"
include "CR_template.pxi"

from sage.libs.pari.pari_instance cimport PariInstance
cdef PariInstance P = sage.libs.pari.pari_instance.pari
from sage.rings.finite_rings.integer_mod import Mod
from sage.rings.padics.pow_computer cimport PowComputer_class

cdef class PowComputer_(PowComputer_base):
    """
    A PowComputer for a capped-relative padic ring or field.
    """
    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field):
        """
        Initialization.

        EXAMPLES::

            sage: R = ZpCR(5)
            sage: type(R.prime_pow)
            <type 'sage.rings.padics.padic_capped_relative_element.PowComputer_'>
            sage: R.prime_pow._prec_type
            'capped-rel'
        """
        self._prec_type = 'capped-rel'
        PowComputer_base.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field)

cdef class pAdicCappedRelativeElement(CRElement):
    """
    Constructs new element with given parent and value.

    INPUT:

    - ``x`` -- value to coerce into a capped relative ring or field

    - ``absprec`` -- maximum number of digits of absolute precision

    - ``relprec`` -- maximum number of digits of relative precision

    EXAMPLES::

        sage: R = Zp(5, 10, 'capped-rel')

    Construct from integers::

        sage: R(3)
        3 + O(5^10)
        sage: R(75)
        3*5^2 + O(5^12)
        sage: R(0)
        0
        sage: R(-1)
        4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + O(5^10)
        sage: R(-5)
        4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + O(5^11)
        sage: R(-7*25)
        3*5^2 + 3*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + O(5^12)

    Construct from rationals::

        sage: R(1/2)
        3 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 2*5^5 + 2*5^6 + 2*5^7 + 2*5^8 + 2*5^9 + O(5^10)
        sage: R(-7875/874)
        3*5^3 + 2*5^4 + 2*5^5 + 5^6 + 3*5^7 + 2*5^8 + 3*5^10 + 3*5^11 + 3*5^12 + O(5^13)
        sage: R(15/425)
        Traceback (most recent call last):
        ...
        ValueError: p divides the denominator

    Construct from IntegerMod::

        sage: R(Integers(125)(3))
        3 + O(5^3)
        sage: R(Integers(5)(3))
        3 + O(5)
        sage: R(Integers(5^30)(3))
        3 + O(5^10)
        sage: R(Integers(5^30)(1+5^23))
        1 + O(5^10)
        sage: R(Integers(49)(3))
        Traceback (most recent call last):
        ...
        TypeError: p does not divide modulus 49

    ::

        sage: R(Integers(48)(3))
        Traceback (most recent call last):
        ...
        TypeError: p does not divide modulus 48

    Some other conversions::

        sage: R(R(5))
        5 + O(5^11)

    Construct from Pari objects::

        sage: R = Zp(5)
        sage: x = pari(123123) ; R(x)
        3 + 4*5 + 4*5^2 + 4*5^3 + 5^4 + 4*5^5 + 2*5^6 + 5^7 + O(5^20)
        sage: R(pari(R(5252)))
        2 + 2*5^3 + 3*5^4 + 5^5 + O(5^20)
        sage: R = Zp(5,prec=5)
        sage: R(pari(-1))
        4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5)
        sage: pari(R(-1))
        4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5)
        sage: pari(R(0))
        0
        sage: R(pari(R(0,5)))
        O(5^5)

    # todo: doctests for converting from other types of p-adic rings

    """
    def lift(self):
        """
        Return an integer or rational congruent to ``self`` modulo ``self``'s
        precision.  If a rational is returned, its denominator will equal
        ``p^ordp(self)``.

        EXAMPLES::

            sage: R = Zp(7,4,'capped-rel'); a = R(8); a.lift()
            8
            sage: R = Qp(7,4); a = R(8); a.lift()
            8
            sage: R = Qp(7,4); a = R(8/7); a.lift()
            8/7
        """
        return self.lift_c()

    cdef lift_c(self):
        """
        Implementation of lift.

        TESTS::

            sage: O(5^5).lift() #indirect doctest
            0
            sage: R = Qp(5); R(0).lift()
            0
            sage: R(5/9).lift()
            264909532335070
            sage: R(9/5).lift()
            9/5
        """
        cdef Integer ans
        cdef Rational ansr
        if self.ordp >= 0:
            ans = PY_NEW(Integer)
            if self.relprec == 0:
                mpz_set_ui(ans.value, 0)
            else:
                mpz_set(ans.value, self.unit)
                mpz_mul(ans.value, ans.value, self.prime_pow.pow_mpz_t_tmp(self.ordp))
            return ans
        else:
            ansr = Rational.__new__(Rational)
            if self.relprec == 0:
                mpq_set_si(ansr.value, 0, 1)
                return self
            else:
                mpz_set(mpq_numref(ansr.value), self.unit)
                mpz_set(mpq_denref(ansr.value), self.prime_pow.pow_mpz_t_tmp(-self.ordp))
            return ansr

    def _pari_(self):
        """
        Converts this element to an equivalent pari element.

        EXAMPLES::

            sage: R = Zp(17, 10); a = ~R(14); pari(a) #indirect doctest
            11 + 3*17 + 17^2 + 6*17^3 + 13*17^4 + 15*17^5 + 10*17^6 + 3*17^7 + 17^8 + 6*17^9 + O(17^10)
            sage: pari(R(0))
            0
            sage: pari(R(0,5))
            O(17^5)
        """
        return self._to_gen()

    cdef pari_gen _to_gen(self):
        """
        Converts this element to an equivalent pari element.

        EXAMPLES::

            sage: R = Zp(5, 10); a = R(17); pari(a) #indirect doctest
            2 + 3*5 + O(5^10)
            sage: pari(R(0))
            0
            sage: pari(R(0,5))
            O(5^5)
            sage: pari(R(0,5)).debug()
            [&=...] PADIC(lg=5):... (precp=0,valp=5):... ... ... ...
                p : [&=...] INT(lg=3):... (+,lgefint=3):... ... 
              p^l : [&=...] INT(lg=3):... (+,lgefint=3):... ... 
                I : [&=...] INT(lg=2):... (0,lgefint=2):... 
        """
        if exactzero(self.ordp):
            return P.new_gen_from_int(0)
        else:
            return P.new_gen_from_padic(self.ordp, self.relprec,
                                        self.prime_pow.prime.value,
                                        self.prime_pow.pow_mpz_t_tmp(self.relprec),
                                        self.unit)
    def _integer_(self, Z=None):
        """
        Returns an integer congruent to this element modulo
        ``p^self.absolute_precision()``.

         EXAMPLES::

            sage: R = Zp(5); a = R(-1); a._integer_()
            95367431640624
         """
        if self.ordp < 0:
            raise ValueError("Cannot form an integer out of a p-adic field element with negative valuation")
        return self.lift_c()

    def residue(self, absprec=1):
        """
        Reduces this element modulo `p^{\mathrm{absprec}}`.

        INPUT:

        - ``absprec`` - a non-negative integer (default: ``1``)

        OUTPUT:

        This element reduced modulo `p^\mathrm{absprec}` as an element of
        `\ZZ/p^\mathrm{absprec}\ZZ`

        EXAMPLES::

            sage: R = Zp(7,4)
            sage: a = R(8)
            sage: a.residue(1)
            1
            sage: a.residue(2)
            8

            sage: K = Qp(7,4)
            sage: a = K(8)
            sage: a.residue(1)
            1
            sage: a.residue(2)
            8
            sage: b = K(1/7)
            sage: b.residue()
            Traceback (most recent call last):
            ...
            ValueError: element must have non-negative valuation in order to compute residue.

        TESTS::

            sage: R = Zp(7,4)
            sage: a = R(8)
            sage: a.residue(0)
            0
            sage: a.residue(-1)
            Traceback (most recent call last):
            ...
            ValueError: cannot reduce modulo a negative power of p.
            sage: a.residue(5)
            Traceback (most recent call last):
            ...
            PrecisionError: not enough precision known in order to compute residue.

        """
        cdef Integer selfvalue, modulus
        cdef long aprec
        if not isinstance(absprec, Integer):
            absprec = Integer(absprec)
        if absprec > self.precision_absolute():
            raise PrecisionError("not enough precision known in order to compute residue.")
        elif absprec < 0:
            raise ValueError("cannot reduce modulo a negative power of p.")
        aprec = mpz_get_ui((<Integer>absprec).value)
        if self.ordp < 0:
            raise ValueError("element must have non-negative valuation in order to compute residue.")
        modulus = PY_NEW(Integer)
        mpz_set(modulus.value, self.prime_pow.pow_mpz_t_tmp(aprec))
        selfvalue = PY_NEW(Integer)
        if self.relprec == 0:
            mpz_set_ui(selfvalue.value, 0)
        else:
            # Need to do this better.
            mpz_mul(selfvalue.value, self.prime_pow.pow_mpz_t_tmp(self.ordp), self.unit)
        return Mod(selfvalue, modulus)

def unpickle_pcre_v1(R, unit, ordp, relprec):
    """
    Unpickles a capped relative element.

    EXAMPLES::

        sage: from sage.rings.padics.padic_capped_relative_element import unpickle_pcre_v1
        sage: R = Zp(5)
        sage: a = unpickle_pcre_v1(R, 17, 2, 5); a
        2*5^2 + 3*5^3 + O(5^7)
    """
    return unpickle_cre_v2(pAdicCappedRelativeElement, R, unit, ordp, relprec)

def base_p_list(Integer n, bint pos, PowComputer_class prime_pow):
    """
    Returns a base-`p` list of digits of ``n``.

    INPUT:

    - ``n`` -- a positive Integer.

    - ``pos`` -- a boolean.  If True, then returns the standard base `p` expansion.
                 Otherwise, the digits lie in the range `-p/2` to `p/2`.

    - ``prime_pow`` -- A PowComputer giving the prime.

    EXAMPLES::

        sage: from sage.rings.padics.padic_capped_relative_element import base_p_list
        sage: base_p_list(192837, True, Zp(5).prime_pow)
        [2, 2, 3, 2, 3, 1, 2, 2]
        sage: 2 + 2*5 + 3*5^2 + 2*5^3 + 3*5^4 + 5^5 + 2*5^6 + 2*5^7
        192837
        sage: base_p_list(192837, False, Zp(5).prime_pow)
        [2, 2, -2, -2, -1, 2, 2, 2]
        sage: 2 + 2*5 - 2*5^2 - 2*5^3 - 5^4 + 2*5^5 + 2*5^6 + 2*5^7
        192837
    """
    if mpz_sgn(n.value) < 0:
        raise ValueError("n must be nonnegative")
    return clist(n.value, prime_pow.prec_cap, pos, prime_pow)
