"""
`p`-adic Capped Absolute Elements

Elements of `p`-adic Rings with Absolute Precision Cap

AUTHORS:

- David Roe
- Genya Zaytman: documentation
- David Harvey: doctests
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
include "sage/libs/linkages/padics/mpz.pxi"
include "CA_template.pxi"

from sage.libs.pari.convert_gmp cimport new_gen_from_padic
from sage.rings.finite_rings.integer_mod import Mod

cdef extern from "transcendantal.c":
    cdef void padiclog(mpz_t ans, const mpz_t a, unsigned long p, unsigned long prec, const mpz_t modulo)
    cdef void padicexp(mpz_t ans, const mpz_t a, unsigned long p, unsigned long prec, const mpz_t modulo)
    cdef void padicexp_Newton(mpz_t ans, const mpz_t a, unsigned long p, unsigned long prec, unsigned long precinit, const mpz_t modulo)

cdef class PowComputer_(PowComputer_base):
    """
    A PowComputer for a capped-absolute `p`-adic ring.
    """
    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field):
        """
        Initialization.

        EXAMPLES::

            sage: R = ZpCA(5)
            sage: type(R.prime_pow)
            <class 'sage.rings.padics.padic_capped_absolute_element.PowComputer_'>
            sage: R.prime_pow._prec_type
            'capped-abs'
        """
        self._prec_type = 'capped-abs'
        PowComputer_base.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field)

cdef class pAdicCappedAbsoluteElement(CAElement):
    """
    Construct new element with given parent and value.

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
        EXAMPLES::

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

            sage: ZpCA(3,3)(1/4).lift()  # indirect doctest
            7
        """
        cdef Integer ans = Integer.__new__(Integer)
        mpz_set(ans.value, self.value)
        return ans

    def __pari__(self):
        """
        Conversion to pari.

        EXAMPLES::

            sage: R = ZpCA(5)
            sage: pari(R(1777))  # indirect doctest                                     # needs sage.libs.pari
            2 + 5^2 + 4*5^3 + 2*5^4 + O(5^20)
            sage: pari(R(0,0))                                                          # needs sage.libs.pari
            O(5^0)
        """
        return self._to_gen()

    cdef pari_gen _to_gen(self):
        """
        Convert this element to an equivalent pari element.

        EXAMPLES::

            sage: R = ZpCA(5, 10); a = R(17); pari(a)  # indirect doctest
            2 + 3*5 + O(5^10)
            sage: pari(R(0,5))                                                          # needs sage.libs.pari
            O(5^5)
            sage: pari(R(0,5)).debug()                                                  # needs sage.libs.pari
            [&=...] PADIC(lg=5):... (precp=0,valp=5):... ... ... ...
                p : [&=...] INT(lg=3):... (+,lgefint=3):... ...
              p^l : [&=...] INT(lg=3):... (+,lgefint=3):... ...
                I : gen_0
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
        Convert this element to an integer.

        TESTS::

            sage: R = ZpCA(5, prec = 4); a = R(642); a
            2 + 3*5 + O(5^4)
            sage: a._integer_()
            17
        """
        return self.lift_c()

    def residue(self, absprec=1, field=None, check_prec=True):
        r"""
        Reduces ``self`` modulo `p^\mathrm{absprec}`.

        INPUT:

        - ``absprec`` -- nonnegative integer (default: 1)

        - ``field`` -- boolean (default: ``None``); whether to return an
          element of GF(p) or Zmod(p)

        - ``check_prec`` -- boolean (default: ``True``); whether to raise an
          error if this element has insufficient precision to determine the
          reduction

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
            1 + 7 + O(7^10)
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
            ValueError: cannot reduce modulo a negative power of p
            sage: a.residue(11)
            Traceback (most recent call last):
            ...
            PrecisionError: not enough precision known in order to compute residue
            sage: a.residue(5, check_prec=False)
            8

            sage: a.residue(field=True).parent()
            Finite Field of size 7

        .. SEEALSO::

            :meth:`_mod_`
        """
        if not isinstance(absprec, Integer):
            absprec = Integer(absprec)
        if check_prec and mpz_cmp_si((<Integer>absprec).value, self.absprec) > 0:
            raise PrecisionError("not enough precision known in order to compute residue")
        elif mpz_sgn((<Integer>absprec).value) < 0:
            raise ValueError("cannot reduce modulo a negative power of p")
        if field is None:
            field = (absprec == 1)
        elif field and absprec != 1:
            raise ValueError("field keyword may only be set at precision 1")
        cdef long aprec = mpz_get_ui((<Integer>absprec).value)
        cdef Integer modulus = Integer.__new__(Integer)
        mpz_set(modulus.value, self.prime_pow.pow_mpz_t_tmp(aprec))
        cdef Integer selfvalue = Integer.__new__(Integer)
        mpz_set(selfvalue.value, self.value)
        if field:
            from sage.rings.finite_rings.finite_field_constructor import GF
            return GF(self.parent().prime())(selfvalue)
        else:
            return Mod(selfvalue, modulus)

    def multiplicative_order(self):
        r"""
        Return the minimum possible multiplicative order of this element.

        OUTPUT:

        The multiplicative order of ``self``.  This is the minimum multiplicative
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
            ans = Integer.__new__(Integer)
            mpz_set_ui(ans.value, 1)
            return ans
        mpz_init(ppow_minus_one)
        mpz_sub_ui(ppow_minus_one, self.prime_pow.pow_mpz_t_tmp(self.absprec), 1)
        if mpz_cmp(self.value, ppow_minus_one) == 0:
            ans = Integer.__new__(Integer)
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

    def _log_binary_splitting(self, aprec, mina=0):
        r"""
        Return ``\log(self)`` for ``self`` equal to 1 in the residue field.

        This is a helper method for :meth:`log`.
        It uses a fast binary splitting algorithm.

        INPUT:

        - ``aprec`` -- integer; the precision to which the result is
          correct. ``aprec`` must not exceed the precision cap of the ring over
          which this element is defined.
        - ``mina`` -- integer (default: 0); the series will check `n` up to
          this valuation (and beyond) to see if they can contribute to the
          series

        .. NOTE::

            The function does not check that its argument ``self`` is
            1 in the residue field. If this assumption is not fulfilled
            the behaviour of the function is not specified.

        ALGORITHM:

        1. Raise `u` to the power `p^v` for a suitable `v` in order
           to make it closer to 1. (`v` is chosen such that `p^v` is
           close to the precision.)

        2. Write

        .. MATH::

            u^{p-1} = \prod_{i=1}^\infty (1 - a_i p^{(v+1)*2^i})

        with `0 \leq a_i < p^{(v+1)*2^i}` and compute
        `\log(1 - a_i p^{(v+1)*2^i})` using the standard Taylor expansion

        .. MATH::

            \log(1 - x) = -x - 1/2 x^2 - 1/3 x^3 - 1/4 x^4 - 1/5 x^5 - \cdots

        together with a binary splitting method.

        3. Divide the result by `p^v`

        The complexity of this algorithm is quasi-linear.

        EXAMPLES::

            sage: r = Qp(5,prec=4)(6)
            sage: r._log_binary_splitting(2)
            5 + O(5^2)
            sage: r._log_binary_splitting(4)
            5 + 2*5^2 + 4*5^3 + O(5^4)
            sage: r._log_binary_splitting(100)
            5 + 2*5^2 + 4*5^3 + O(5^4)

            sage: r = Zp(5,prec=4,type='fixed-mod')(6)
            sage: r._log_binary_splitting(5)
            5 + 2*5^2 + 4*5^3
        """
        cdef unsigned long p
        cdef unsigned long prec = min(aprec, self.absprec)
        cdef pAdicCappedAbsoluteElement ans, unit

        if mpz_fits_slong_p(self.prime_pow.prime.value) == 0:
            raise NotImplementedError("the prime %s does not fit in a long" % self.prime_pow.prime)
        p = self.prime_pow.prime

        ans = self._new_c()
        ans.absprec = prec
        unit = self.unit_part()
        sig_on()
        padiclog(ans.value, unit.value, p, prec, self.prime_pow.pow_mpz_t_tmp(prec))
        sig_off()

        return ans

    def _exp_binary_splitting(self, aprec):
        r"""
        Compute the exponential power series of this element.

        This is a helper method for :meth:`exp`.

        INPUT:

        - ``aprec`` -- integer; the precision to which to compute the
          exponential

        .. NOTE::

            The function does not check that its argument ``self`` is
            the disk of convergence of ``exp``. If this assumption is not
            fulfilled the behaviour of the function is not specified.

        ALGORITHM:

        Write

        .. MATH::

            self = \sum_{i=1}^\infty a_i p^{2^i}

        with `0 \leq a_i < p^{2^i}` and compute
        `\exp(a_i p^{2^i})` using the standard Taylor expansion

        .. MATH::

            \exp(x) = 1 + x + x^2/2 + x^3/6 + x^4/24 + \cdots

        together with a binary splitting method.

        The binary complexity of this algorithm is quasi-linear.

        EXAMPLES::

            sage: R = Zp(7,5)
            sage: x = R(7)
            sage: x.exp(algorithm='binary_splitting')   # indirect doctest
            1 + 7 + 4*7^2 + 2*7^3 + O(7^5)
        """
        cdef unsigned long p
        cdef unsigned long prec = aprec
        cdef pAdicCappedAbsoluteElement ans

        if mpz_fits_slong_p(self.prime_pow.prime.value) == 0:
            raise NotImplementedError("the prime %s does not fit in a long" % self.prime_pow.prime)
        p = self.prime_pow.prime

        ans = self._new_c()
        ans.absprec = prec
        sig_on()
        padicexp(ans.value, self.value, p, prec, self.prime_pow.pow_mpz_t_tmp(prec))
        sig_off()

        return ans

    def _exp_newton(self, aprec, log_algorithm=None):
        r"""
        Compute the exponential power series of this element.

        This is a helper method for :meth:`exp`.

        INPUT:

        - ``aprec`` -- integer; the precision to which to compute the
          exponential

        - ``log_algorithm`` -- (default: ``None``) the algorithm used for
          computing the logarithm. This attribute is passed to the log
          method. See :meth:`log` for more details about the possible
          algorithms.

        .. NOTE::

            The function does not check that its argument ``self`` is
            the disk of convergence of ``exp``. If this assumption is not
            fulfilled the behaviour of the function is not specified.

        ALGORITHM:

        Solve the equation `\log(x) = self` using the Newton scheme::

        .. MATH::

            x_{i+1} = x_i \cdot (1 + self - \log(x_i))

        The binary complexity of this algorithm is roughly the same
        than that of the computation of the logarithm.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<w> = Zq(7^2,5)
            sage: x = R(7*w)
            sage: x.exp(algorithm='newton')   # indirect doctest
            1 + w*7 + (4*w + 2)*7^2 + (w + 6)*7^3 + 5*7^4 + O(7^5)
        """
        cdef unsigned long p
        cdef unsigned long prec = aprec
        cdef pAdicCappedAbsoluteElement ans

        if mpz_fits_slong_p(self.prime_pow.prime.value) == 0:
            raise NotImplementedError("the prime %s does not fit in a long" % self.prime_pow.prime)
        p = self.prime_pow.prime

        ans = self._new_c()
        ans.absprec = prec
        mpz_set_ui(ans.value, 1)
        sig_on()
        if p == 2:
            padicexp_Newton(ans.value, self.value, p, prec, 2, self.prime_pow.pow_mpz_t_tmp(prec))
        else:
            padicexp_Newton(ans.value, self.value, p, prec, 1, self.prime_pow.pow_mpz_t_tmp(prec))
        sig_off()

        return ans


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
