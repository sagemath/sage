"""
p-Adic Fixed-Mod Element

Elements of p-Adic Rings with Fixed Modulus

AUTHORS:

- David Roe
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
include "FM_template.pxi"

from sage.libs.pari.pari_instance cimport PariInstance
cdef PariInstance P = sage.libs.pari.pari_instance.pari
from sage.rings.finite_rings.integer_mod import Mod

cdef class PowComputer_(PowComputer_base):
    """
    A PowComputer for a fixed-modulus padic ring.
    """
    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field):
        """
        Initialization.

        EXAMPLES::

            sage: R = ZpFM(5)
            sage: type(R.prime_pow)
            <type 'sage.rings.padics.padic_fixed_mod_element.PowComputer_'>
            sage: R.prime_pow._prec_type
            'fixed-mod'
        """
        self._prec_type = 'fixed-mod'
        PowComputer_base.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field)

cdef class pAdicFixedModElement(FMElement):
    r"""
    INPUT:

    - ``parent`` -- a ``pAdicRingFixedMod`` object.

    - ``x`` -- input data to be converted into the parent.

    - ``absprec`` -- ignored; for compatibility with other `p`-adic rings

    - ``relprec`` -- ignored; for compatibility with other `p`-adic rings

    .. NOTE::

        The following types are currently supported for x:

        - Integers
        - Rationals -- denominator must be relatively prime to `p`
        - FixedMod `p`-adics
        - Elements of ``IntegerModRing(p^k)`` for ``k`` less than or equal
          to the modulus

        The following types should be supported eventually:

        - Finite precision `p`-adics
        - Lazy `p`-adics
        - Elements of local extensions of THIS `p`-adic ring that actually
          lie in `\ZZ_p`

    EXAMPLES::

        sage: R = Zp(5, 20, 'fixed-mod', 'terse')

    Construct from integers::

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

    Construct from rationals::

        sage: R(1/2)
        47683715820313 + O(5^20)
        sage: R(-7875/874)
        9493096742250 + O(5^20)
        sage: R(15/425)
        Traceback (most recent call last):
        ...
        ValueError: p divides denominator

    Construct from IntegerMod::

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
        TypeError: p does not divide modulus 49

        sage: R(Integers(48)(3))
        Traceback (most recent call last):
        ...
        TypeError: p does not divide modulus 48

    Some other conversions::

        sage: R(R(5))
        5 + O(5^20)

    .. TODO:: doctests for converting from other types of `p`-adic rings
    """
    def lift(self):
        r"""
        Return an integer congruent to ``self`` modulo the precision.

        .. WARNING::

            Since fixed modulus elements don't track their precision,
            the result may not be correct modulo
            `i^{\mathrm{prec_cap}}` if the element was defined by
            constructions that lost precision.

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod'); a = R(8); a.lift()
            8
            sage: type(a.lift())
            <type 'sage.rings.integer.Integer'>
        """
        return self.lift_c()

    cdef lift_c(self):
        """
        Returns an integer congruent to this element modulo the precision.

        .. WARNING::

            Since fixed modulus elements don't track their precision,
            the result may not be correct modulo
            `i^{\mbox{prec_cap}}` if the element was defined by
            constructions that lost precision.

        EXAMPLES::

            sage: R = ZpFM(7,4); a = R(8); a.lift() # indirect doctest
            8
        """
        cdef Integer ans = PY_NEW(Integer)
        mpz_set(ans.value, self.value)
        return ans

    def _pari_(self):
        """
        Conversion to PARI.

        EXAMPLES::

            sage: R = ZpCA(5)
            sage: pari(R(1777)) #indirect doctest
            2 + 5^2 + 4*5^3 + 2*5^4 + O(5^20)
        """
        return self._to_gen()

    cdef pari_gen _to_gen(self):
        """
        Convert ``self`` to an equivalent pari element.

        EXAMPLES::

            sage: R = ZpFM(5, 10); a = R(17); pari(a) # indirect doctest
            2 + 3*5 + O(5^10)
            sage: pari(R(0))
            O(5^10)
            sage: pari(R(0,5))
            O(5^10)
            sage: pari(R(0)).debug()
            [&=...] PADIC(lg=5):... (precp=0,valp=10):... ... ... ...
                p : [&=...] INT(lg=3):... (+,lgefint=3):... ... 
              p^l : [&=...] INT(lg=3):... (+,lgefint=3):... ... 
                I : [&=...] INT(lg=2):... (0,lgefint=2):... 

        This checks that :trac:`15653` is fixed::

            sage: x = polygen(ZpFM(3,10))
            sage: (x^3 + x + 1)._pari_().poldisc()
            2 + 3 + 2*3^2 + 3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + O(3^10)
        """
        cdef long val
        # Let val be the valuation of self, holder (defined in the
        # linkage file) be the unit part.
        if mpz_sgn(self.value) == 0:
            # Special case for zero: maximal valuation and 0 unit part
            val = self.prime_pow.prec_cap
            mpz_set_ui(holder.value, 0)
        else:
            val = mpz_remove(holder.value, self.value, self.prime_pow.prime.value)
        return P.new_gen_from_padic(val, self.prime_pow.prec_cap - val,
                                    self.prime_pow.prime.value,
                                    self.prime_pow.pow_mpz_t_tmp(self.prime_pow.prec_cap - val),
                                    holder.value)

    def _integer_(self, Z=None):
        """
        Return an integer congruent to ``self`` modulo the precision.

        .. WARNING::

            Since fixed modulus elements don't track their precision,
            the result may not be correct modulo
            `p^{\mathrm{prec_cap}}` if the element was defined by
            constructions that lost precision.

        EXAMPLES::

            sage: R = ZpFM(5); R(-1)._integer_()
            95367431640624
        """
        return self.lift_c()

    def residue(self, absprec=1):
        r"""
        Reduce ``self`` modulo `p^\mathrm{absprec}`.

        INPUT:

        - ``absprec`` -- an integer (default: ``1``)

        OUTPUT:

        This element reduced modulo `p^\mathrm{absprec}` as an element of
        `\ZZ/p^\mathrm{absprec}\ZZ`.

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod')
            sage: a = R(8)
            sage: a.residue(1)
            1
            sage: a.residue(2)
            8

        TESTS::

            sage: R = Zp(7,4,'fixed-mod')
            sage: a = R(8)
            sage: a.residue(0)
            0
            sage: a.residue(-1)
            Traceback (most recent call last):
            ...
            ValueError: Cannot reduce modulo a negative power of p.
            sage: a.residue(5)
            Traceback (most recent call last):
            ...
            PrecisionError: Not enough precision known in order to compute residue.

        """
        cdef Integer selfvalue, modulus
        if not isinstance(absprec, Integer):
            absprec = Integer(absprec)
        if absprec > self.precision_absolute():
            raise PrecisionError("Not enough precision known in order to compute residue.")
        elif absprec < 0:
            raise ValueError("Cannot reduce modulo a negative power of p.")
        cdef long aprec = mpz_get_ui((<Integer>absprec).value)
        modulus = PY_NEW(Integer)
        mpz_set(modulus.value, self.prime_pow.pow_mpz_t_tmp(aprec))
        selfvalue = PY_NEW(Integer)
        mpz_set(selfvalue.value, self.value)
        return Mod(selfvalue, modulus)

    def multiplicative_order(self):
        r"""
        Return the minimum possible multiplicative order of ``self``.

        OUTPUT:

        an integer -- the multiplicative order of this element.  This is the
        minimum multiplicative order of all elements of `\ZZ_p` lifting this
        element to infinite precision.

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
        mpz_sub_ui(tmp, self.prime_pow.pow_mpz_t_top(), 1)
        if mpz_cmp(self.value, tmp) == 0:
            ans = PY_NEW(Integer)
            mpz_set_ui(ans.value, 2)
            return ans
        # check if self is an approximation to a teichmuller lift:
        mpz_powm(tmp, self.value, self.prime_pow.prime.value, self.prime_pow.pow_mpz_t_top())
        if mpz_cmp(tmp, self.value) == 0:
            mpz_clear(tmp)
            return self.residue(1).multiplicative_order()
        else:
            mpz_clear(tmp)
            return infinity

def make_pAdicFixedModElement(parent, value):
    """
    Unpickles a fixed modulus element.

    EXAMPLES::

        sage: from sage.rings.padics.padic_fixed_mod_element import make_pAdicFixedModElement
        sage: R = ZpFM(5)
        sage: a = make_pAdicFixedModElement(R, 17*25); a
        2*5^2 + 3*5^3 + O(5^20)
    """
    return unpickle_fme_v2(pAdicFixedModElement, parent, value)

