# sage_setup: distribution = sagemath-pari
# sage.doctest: needs sage.rings.padics
"""
`p`-adic Capped Relative Elements

Elements of `p`-adic Rings with Capped Relative Precision

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

from sage.libs.pari import pari
from sage.libs.pari.convert_gmp cimport new_gen_from_padic
from sage.rings.finite_rings.integer_mod import Mod
from sage.rings.padics.pow_computer cimport PowComputer_class

cdef extern from "transcendantal.c":
    cdef void padiclog(mpz_t ans, const mpz_t a, unsigned long p, unsigned long prec, const mpz_t modulo)
    cdef void padicexp(mpz_t ans, const mpz_t a, unsigned long p, unsigned long prec, const mpz_t modulo)
    cdef void padicexp_Newton(mpz_t ans, const mpz_t a, unsigned long p, unsigned long prec, unsigned long precinit, const mpz_t modulo)


cdef class PowComputer_(PowComputer_base):
    """
    A PowComputer for a capped-relative `p`-adic ring or field.
    """
    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field):
        """
        Initialization.

        EXAMPLES::

            sage: R = ZpCR(5)
            sage: type(R.prime_pow)
            <class 'sage.rings.padics.padic_capped_relative_element.PowComputer_'>
            sage: R.prime_pow._prec_type
            'capped-rel'
        """
        self._prec_type = 'capped-rel'
        PowComputer_base.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field)

cdef class pAdicCappedRelativeElement(CRElement):
    """
    Construct new element with given parent and value.

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
        sage: x = pari(123123) ; R(x)                                                   # needs sage.libs.pari
        3 + 4*5 + 4*5^2 + 4*5^3 + 5^4 + 4*5^5 + 2*5^6 + 5^7 + O(5^20)
        sage: R(pari(R(5252)))
        2 + 2*5^3 + 3*5^4 + 5^5 + O(5^20)
        sage: R = Zp(5,prec=5)
        sage: R(pari(-1))
        4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5)
        sage: pari(R(-1))                                                               # needs sage.libs.pari
        4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5)
        sage: pari(R(0))                                                                # needs sage.libs.pari
        0
        sage: R(pari(R(0,5)))
        O(5^5)

    .. TODO:: doctests for converting from other types of `p`-adic rings
    """
    def lift(self):
        r"""
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

            sage: O(5^5).lift()  # indirect doctest
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

    def __pari__(self):
        """
        Convert this element to an equivalent pari element.

        EXAMPLES::

            sage: R = Zp(17, 10); a = ~R(14); pari(a)  # indirect doctest
            11 + 3*17 + 17^2 + 6*17^3 + 13*17^4 + 15*17^5 + 10*17^6 + 3*17^7 + 17^8 + 6*17^9 + O(17^10)
            sage: pari(R(0))                                                            # needs sage.libs.pari
            0
            sage: pari(R(0,5))                                                          # needs sage.libs.pari
            O(17^5)
        """
        return self._to_gen()

    cdef pari_gen _to_gen(self):
        """
        Convert this element to an equivalent pari element.

        EXAMPLES::

            sage: R = Zp(5, 10); a = R(17); pari(a)  # indirect doctest
            2 + 3*5 + O(5^10)
            sage: pari(R(0))                                                            # needs sage.libs.pari
            0
            sage: pari(R(0,5))                                                          # needs sage.libs.pari
            O(5^5)
            sage: pari(R(0,5)).debug()                                                  # needs sage.libs.pari
            [&=...] PADIC(lg=5):... (precp=0,valp=5):... ... ... ...
                p : [&=...] INT(lg=3):... (+,lgefint=3):... ...
              p^l : [&=...] INT(lg=3):... (+,lgefint=3):... ...
                I : gen_0
        """
        if exactzero(self.ordp):
            return pari.zero()
        else:
            return new_gen_from_padic(self.ordp, self.relprec,
                                      self.prime_pow.prime.value,
                                      self.prime_pow.pow_mpz_t_tmp(self.relprec),
                                      self.unit)

    def _integer_(self, Z=None):
        r"""
        Return an integer congruent to this element modulo
        ``p^self.absolute_precision()``.

         EXAMPLES::

            sage: R = Zp(5); a = R(-1); a._integer_()
            95367431640624
         """
        if self.ordp < 0:
            raise ValueError("cannot form an integer out of a p-adic field element with negative valuation")
        return self.lift_c()

    def residue(self, absprec=1, field=None, check_prec=True):
        r"""
        Reduce this element modulo `p^{\mathrm{absprec}}`.

        INPUT:

        - ``absprec`` -- nonnegative integer (default: 1)

        - ``field`` -- boolean (default: ``None``); whether to return an element
          of `\GF{p}` or `\ZZ / p\ZZ`

        - ``check_prec`` -- boolean (default: ``True``); whether to raise
          an error if this element has insufficient precision to determine
          the reduction

        OUTPUT:

        This element reduced modulo `p^\mathrm{absprec}` as an element of
        `\ZZ/p^\mathrm{absprec}\ZZ`.

        EXAMPLES::

            sage: R = Zp(7,4)
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
            1 + 7 + O(7^4)
            sage: c.parent()
            7-adic Ring with capped relative precision 4

        For elements in a field, application of ``% p^n`` always returns
        zero, the remainder of the division by ``p^n``::

            sage: K = Qp(7,4)
            sage: a = K(8)
            sage: a.residue(2)
            8
            sage: a % 7^2
            1 + 7 + O(7^4)

            sage: b = K(1/7)
            sage: b.residue()
            Traceback (most recent call last):
            ...
            ValueError: element must have nonnegative valuation in order to compute residue

        TESTS::

            sage: R = Zp(7,4)
            sage: a = R(8)
            sage: a.residue(0)
            0
            sage: a.residue(-1)
            Traceback (most recent call last):
            ...
            ValueError: cannot reduce modulo a negative power of p
            sage: a.residue(5)
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
        cdef Integer selfvalue, modulus
        cdef long aprec
        if not isinstance(absprec, Integer):
            absprec = Integer(absprec)
        if check_prec and absprec > self.precision_absolute():
            raise PrecisionError("not enough precision known in order to compute residue")
        elif absprec < 0:
            raise ValueError("cannot reduce modulo a negative power of p")
        if self.ordp < 0:
            raise ValueError("element must have nonnegative valuation in order to compute residue")
        if field is None:
            field = (absprec == 1)
        elif field and absprec != 1:
            raise ValueError("field keyword may only be set at precision 1")
        modulus = PY_NEW(Integer)
        aprec = mpz_get_ui((<Integer>absprec).value)
        mpz_set(modulus.value, self.prime_pow.pow_mpz_t_tmp(aprec))
        selfvalue = PY_NEW(Integer)
        if self.relprec == 0:
            mpz_set_ui(selfvalue.value, 0)
        else:
            # Need to do this better.
            mpz_mul(selfvalue.value, self.prime_pow.pow_mpz_t_tmp(self.ordp), self.unit)
        if field:
            from sage.rings.finite_rings.finite_field_constructor import GF
            return GF(self.parent().prime())(selfvalue)
        else:
            return Mod(selfvalue, modulus)

    def _log_binary_splitting(self, aprec, mina=0):
        r"""
        Return ``\log(self)`` for ``self`` equal to 1 in the residue field.

        This is a helper method for :meth:`log`.
        It uses a fast binary splitting algorithm.

        INPUT:

        - ``aprec`` -- integer; the precision to which the result is
          correct; ``aprec`` must not exceed the precision cap of the ring over
          which this element is defined
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
        cdef unsigned long prec = min(aprec, self.relprec)
        cdef pAdicCappedRelativeElement ans

        if mpz_fits_slong_p(self.prime_pow.prime.value) == 0:
            raise NotImplementedError("the prime %s does not fit in a long" % self.prime_pow.prime)
        p = self.prime_pow.prime

        ans = self._new_c()
        ans.ordp = 0
        ans.relprec = prec
        sig_on()
        padiclog(ans.unit, self.unit, p, prec, self.prime_pow.pow_mpz_t_tmp(prec))
        sig_off()
        ans._normalize()

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
        cdef pAdicCappedRelativeElement ans
        cdef Integer selfint = self.lift_c()

        if mpz_fits_slong_p(self.prime_pow.prime.value) == 0:
            raise NotImplementedError("the prime %s does not fit in a long" % self.prime_pow.prime)
        p = self.prime_pow.prime

        ans = self._new_c()
        ans.ordp = 0
        ans.relprec = prec
        sig_on()
        padicexp(ans.unit, selfint.value, p, prec, self.prime_pow.pow_mpz_t_tmp(prec))
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
          computing the logarithm; this attribute is passed to the :meth:`log`
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
        cdef pAdicCappedRelativeElement ans
        cdef Integer selfint = self.lift_c()

        if mpz_fits_slong_p(self.prime_pow.prime.value) == 0:
            raise NotImplementedError("the prime %s does not fit in a long" % self.prime_pow.prime)
        p = self.prime_pow.prime

        ans = self._new_c()
        ans.ordp = 0
        ans.relprec = prec
        mpz_set_ui(ans.unit, 1)
        sig_on()
        if p == 2:
            padicexp_Newton(ans.unit, selfint.value, p, prec, 2, self.prime_pow.pow_mpz_t_tmp(prec))
        else:
            padicexp_Newton(ans.unit, selfint.value, p, prec, 1, self.prime_pow.pow_mpz_t_tmp(prec))
        sig_off()

        return ans


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
    r"""
    Return a base-`p` list of digits of ``n``.

    INPUT:

    - ``n`` -- a positive :class:`Integer`

    - ``pos`` -- boolean; if ``True``, then returns the standard base `p`
      expansion, otherwise the digits lie in the range `-p/2` to `p/2`

    - ``prime_pow`` -- a :class:`PowComputer` giving the prime

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
    cdef expansion_mode mode = simple_mode if pos else smallest_mode
    # We need a p-adic element to feed to ExpansionIter before resetting its curvalue
    from sage.rings.padics.factory import Zp
    p = prime_pow.prime
    dummy = Zp(p)(0)
    cdef ExpansionIter expansion = ExpansionIter(dummy, n.exact_log(p) + 2, mode)
    mpz_set(expansion.curvalue, n.value)
    return trim_zeros(list(expansion))
