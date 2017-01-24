"""
`p`-Adic Capped Absolute Elements

Elements of `p`-Adic Rings with Absolute Precision Cap

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
include "CA_template.pxi"

from sage.libs.pari.convert_gmp cimport new_gen_from_padic
from sage.rings.finite_rings.integer_mod import Mod

cdef class PowComputer_(PowComputer_base):
    """
    A PowComputer for a capped-absolute padic ring.
    """
    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field):
        """
        Initialization.

        EXAMPLES::

            sage: R = ZpCA(5)
            sage: type(R.prime_pow)
            <type 'sage.rings.padics.padic_capped_absolute_element.PowComputer_'>
            sage: R.prime_pow._prec_type
            'capped-abs'
        """
        self._prec_type = 'capped-abs'
        PowComputer_base.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field)

cdef class pAdicCappedAbsoluteElement(CAElement):
    """
    Constructs new element with given parent and value.

    INPUT:

    - ``x`` -- value to coerce into a capped absolute ring

    - ``absprec`` -- maximum number of digits of absolute precision

    - ``relprec`` -- maximum number of digits of relative precision

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
    """
    def lift(self):
        """
            sage: R = ZpCA(3)
            sage: R(10).lift()
            10
            sage: R(-1).lift()
            3486784400
        """
        return self.lift_c()

    cdef lift_c(self):
        """
        Implementation of lift.

        TESTS::

            sage: ZpCA(3,3)(1/4).lift() # indirect doctest
            7
        """
        cdef Integer ans = PY_NEW(Integer)
        mpz_set(ans.value, self.value)
        return ans

    def _pari_(self):
        """
        Conversion to pari.

        EXAMPLES::

            sage: R = ZpCA(5)
            sage: pari(R(1777)) #indirect doctest
            2 + 5^2 + 4*5^3 + 2*5^4 + O(5^20)
            sage: pari(R(0,0))
            O(5^0)
        """
        return self._to_gen()

    cdef pari_gen _to_gen(self):
        """
        Converts this element to an equivalent pari element.

        EXAMPLES::

            sage: R = ZpCA(5, 10); a = R(17); pari(a) #indirect doctest
            2 + 3*5 + O(5^10)
            sage: pari(R(0,5))
            O(5^5)
            sage: pari(R(0,5)).debug()
            [&=...] PADIC(lg=5):... (precp=0,valp=5):... ... ... ...
                p : [&=...] INT(lg=3):... (+,lgefint=3):... ... 
              p^l : [&=...] INT(lg=3):... (+,lgefint=3):... ... 
                I : [&=...] INT(lg=2):... (0,lgefint=2):... 
        """
        cdef long val
        # Let val be the valuation of self, holder (defined in the
        # linkage file) be the unit part.
        if mpz_sgn(self.value) == 0:
            # Special case for zero: maximal valuation and 0 unit part
            val = self.absprec
            mpz_set_ui(holder.value, 0)
        else:
            val = mpz_remove(holder.value, self.value, self.prime_pow.prime.value)
        return new_gen_from_padic(val, self.absprec - val,
                                  self.prime_pow.prime.value,
                                  self.prime_pow.pow_mpz_t_tmp(self.absprec - val),
                                  holder.value)

    def _integer_(self, Z=None):
        r"""
        Converts this element to an integer.

        TESTS::

            sage: R = ZpCA(5, prec = 4); a = R(642); a
            2 + 3*5 + O(5^4)
            sage: a._integer_()
            17
        """
        return self.lift_c()

    def residue(self, absprec=1):
        r"""
        Reduces ``self`` modulo `p^\mathrm{absprec}`.

        INPUT:

        - ``absprec`` - a non-negative integer (default: 1)

        OUTPUT:

        This element reduced modulo `p^\mathrm{absprec}` as an element of
        `\ZZ/p^\mathrm{absprec}\ZZ`

         EXAMPLES::

            sage: R = Zp(7,10,'capped-abs')
            sage: a = R(8)
            sage: a.residue(1)
            1

        This is different from applying ``% p^n`` which returns an element in
        the same ring::

            sage: b = a.residue(2); b
            8
            sage: b.parent()
            Ring of integers modulo 49
            sage: c = a % 7^2; c
            1 + 7 + O(7^8)
            sage: c.parent()
            7-adic Ring with capped absolute precision 10

        Note that reduction of ``c`` dropped to the precision of the unit part
        of ``7^2``, see :meth:`_mod_`::

            sage: R(7^2).unit_part()
            1 + O(7^8)

        TESTS::

            sage: a.residue(0)
            0
            sage: a.residue(-1)
            Traceback (most recent call last):
            ...
            ValueError: cannot reduce modulo a negative power of p.
            sage: a.residue(11)
            Traceback (most recent call last):
            ...
            PrecisionError: not enough precision known in order to compute residue.

        .. SEEALSO::

            :meth:`_mod_`

        """
        cdef Integer selfvalue, modulus
        if not isinstance(absprec, Integer):
            absprec = Integer(absprec)
        if mpz_cmp_si((<Integer>absprec).value, self.absprec) > 0:
            raise PrecisionError("not enough precision known in order to compute residue.")
        elif mpz_sgn((<Integer>absprec).value) < 0:
            raise ValueError("cannot reduce modulo a negative power of p.")
        cdef long aprec = mpz_get_ui((<Integer>absprec).value)
        modulus = PY_NEW(Integer)
        mpz_set(modulus.value, self.prime_pow.pow_mpz_t_tmp(aprec))
        selfvalue = PY_NEW(Integer)
        mpz_set(selfvalue.value, self.value)
        return Mod(selfvalue, modulus)

    def multiplicative_order(self):
        r"""
        Returns the minimum possible multiplicative order of this element.

        OUTPUT:
        the multiplicative order of self.  This is the minimum multiplicative
        order of all elements of `\ZZ_p` lifting ``self`` to infinite
        precision.

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
        mpz_sub_ui(ppow_minus_one, self.prime_pow.pow_mpz_t_tmp(self.absprec), 1)
        if mpz_cmp(self.value, ppow_minus_one) == 0:
            ans = PY_NEW(Integer)
            mpz_set_ui(ans.value, 2)
            mpz_clear(ppow_minus_one)
            return ans
        # check if self is an approximation to a teichmuller lift:
        mpz_powm(ppow_minus_one, self.value, self.prime_pow.prime.value, self.prime_pow.pow_mpz_t_tmp(self.absprec))
        if mpz_cmp(ppow_minus_one, self.value) == 0:
            mpz_clear(ppow_minus_one)
            return self.residue(1).multiplicative_order()
        else:
            mpz_clear(ppow_minus_one)
            return infinity

def make_pAdicCappedAbsoluteElement(parent, x, absprec):
    """
    Unpickles a capped absolute element.

    EXAMPLES::

        sage: from sage.rings.padics.padic_capped_absolute_element import make_pAdicCappedAbsoluteElement
        sage: R = ZpCA(5)
        sage: a = make_pAdicCappedAbsoluteElement(R, 17*25, 5); a
        2*5^2 + 3*5^3 + O(5^5)
    """
    return unpickle_cae_v2(pAdicCappedAbsoluteElement, parent, x, absprec)
