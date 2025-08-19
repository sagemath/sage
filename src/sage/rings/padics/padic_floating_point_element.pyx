r"""
`p`-adic Floating Point Elements

Elements of `p`-adic Rings with Floating Point Precision

AUTHORS:

- David Roe: initial version (2016-12-6)
"""

#*****************************************************************************
#       Copyright (C) 2016 David Roe <roed.math@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/libs/linkages/padics/mpz.pxi"
include "FP_template.pxi"

from sage.libs.pari import pari
from sage.libs.pari.convert_gmp cimport new_gen_from_padic
from sage.rings.finite_rings.integer_mod import Mod

cdef extern from "transcendantal.c":
    cdef void padicexp(mpz_t ans, const mpz_t a, unsigned long p, unsigned long prec, const mpz_t modulo)
    cdef void padicexp_Newton(mpz_t ans, const mpz_t a, unsigned long p, unsigned long prec, unsigned long precinit, const mpz_t modulo)


cdef class PowComputer_(PowComputer_base):
    """
    A PowComputer for a floating-point `p`-adic ring or field.
    """
    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field):
        """
        Initialization.

        EXAMPLES::

            sage: R = ZpFP(5)
            sage: type(R.prime_pow)
            <class 'sage.rings.padics.padic_floating_point_element.PowComputer_'>
            sage: R.prime_pow._prec_type
            'floating-point'
        """
        self._prec_type = 'floating-point'
        PowComputer_base.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field)

cdef class pAdicFloatingPointElement(FPElement):
    """
    Construct new element with given parent and value.

    INPUT:

    - ``x`` -- value to coerce into a floating point ring or field

    - ``absprec`` -- maximum number of digits of absolute precision

    - ``relprec`` -- maximum number of digits of relative precision

    EXAMPLES::

        sage: R = Zp(5, 10, 'floating-point')

    Construct from integers::

        sage: R(3)
        3
        sage: R(75)
        3*5^2
        sage: R(0)
        0
        sage: R(-1)
        4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9
        sage: R(-5)
        4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10
        sage: R(-7*25)
        3*5^2 + 3*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11

    Construct from rationals::

        sage: R(1/2)
        3 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 2*5^5 + 2*5^6 + 2*5^7 + 2*5^8 + 2*5^9
        sage: R(-7875/874)
        3*5^3 + 2*5^4 + 2*5^5 + 5^6 + 3*5^7 + 2*5^8 + 3*5^10 + 3*5^11 + 3*5^12
        sage: R(15/425)
        Traceback (most recent call last):
        ...
        ValueError: p divides the denominator

    Construct from IntegerMod::

        sage: R(Integers(125)(3))
        3
        sage: R(Integers(5)(3))
        3
        sage: R(Integers(5^30)(3))
        3
        sage: R(Integers(5^30)(1+5^23))
        1
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
        5

    Construct from Pari objects::

        sage: R = ZpFP(5)
        sage: x = pari(123123) ; R(x)                                                   # needs sage.libs.pari
        3 + 4*5 + 4*5^2 + 4*5^3 + 5^4 + 4*5^5 + 2*5^6 + 5^7
        sage: R(pari(R(5252)))
        2 + 2*5^3 + 3*5^4 + 5^5
        sage: R = ZpFP(5,prec=5)
        sage: R(pari(-1))
        4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4
        sage: pari(R(-1))                                                               # needs sage.libs.pari
        4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5)
        sage: pari(R(0))                                                                # needs sage.libs.pari
        0
        sage: R(pari(R(0,5)))
        0

    # todo: doctests for converting from other types of p-adic rings
    """
    def lift(self):
        r"""
        Return an integer or rational congruent to ``self`` modulo ``self``'s
        precision.  If a rational is returned, its denominator will equal
        ``p^ordp(self)``.

        This method will raise a :exc:`ValueError` when this element
        is infinity.

        EXAMPLES::

            sage: R = Zp(7,4,'floating-point'); a = R(8); a.lift()
            8
            sage: R = QpFP(7,4); a = R(8); a.lift()
            8
            sage: R = QpFP(7,4); a = R(8/7); a.lift()
            8/7
        """
        return self.lift_c()

    cdef lift_c(self):
        r"""
        Implementation of lift.

        TESTS::

            sage: ZpFP(5)(0).lift()  # indirect doctest
            0
            sage: R = QpFP(5); R(0).lift()
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
            if very_pos_val(self.ordp):
                mpz_set_ui(ans.value, 0)
            else:
                mpz_set(ans.value, self.unit)
                mpz_mul(ans.value, ans.value, self.prime_pow.pow_mpz_t_tmp(self.ordp))
            return ans
        else:
            ansr = Rational.__new__(Rational)
            if very_neg_val(self.ordp):
                raise ValueError("infinity cannot be lifted to an integer or rational")
            mpz_set(mpq_numref(ansr.value), self.unit)
            mpz_set(mpq_denref(ansr.value), self.prime_pow.pow_mpz_t_tmp(-self.ordp))
            return ansr

    def __pari__(self):
        """
        Convert this element to an equivalent pari element.

        EXAMPLES::

            sage: R = ZpFP(17, 10); a = ~R(14); pari(a)  # indirect doctest
            11 + 3*17 + 17^2 + 6*17^3 + 13*17^4 + 15*17^5 + 10*17^6 + 3*17^7 + 17^8 + 6*17^9 + O(17^10)
            sage: pari(R(0))                                                            # needs sage.libs.pari
            0
        """
        return self._to_gen()

    cdef pari_gen _to_gen(self):
        """
        Convert this element to an equivalent pari element.

        EXAMPLES::

            sage: R = ZpFP(5, 10); a = R(17); pari(a)  # indirect doctest
            2 + 3*5 + O(5^10)
            sage: pari(R(0))                                                            # needs sage.libs.pari
            0
        """
        if very_pos_val(self.ordp):
            return pari.zero()
        elif very_neg_val(self.ordp):
            raise ValueError("no analogue of p-adic infinity in pari")
        else:
            return new_gen_from_padic(self.ordp, self.prime_pow.prec_cap,
                                      self.prime_pow.prime.value,
                                      self.prime_pow.pow_mpz_t_top(),
                                      self.unit)

    def _integer_(self, Z=None):
        r"""
        Return an integer congruent to this element modulo
        ``p^self.absolute_precision()``.

         EXAMPLES::

            sage: R = ZpFP(5); a = R(-1); a._integer_()
            95367431640624
         """
        if self.ordp < 0:
            raise ValueError("cannot form an integer out of a p-adic field element with negative valuation")
        return self.lift_c()

    def residue(self, absprec=1, field=None, check_prec=False):
        r"""
        Reduce this element modulo `p^{\mathrm{absprec}}`.

        INPUT:

        - ``absprec`` -- nonnegative integer (default: 1)

        - ``field`` -- boolean (default: ``None``); whether to return an
          element of `\GF{p}` or `\ZZ / p\ZZ`

        - ``check_prec`` -- ignored (for compatibility with other types)

        OUTPUT:

        This element reduced modulo `p^{\mathrm{absprec}}` as an element of
        `\ZZ/p^{\mathrm{absprec}}\ZZ`.

        EXAMPLES::

            sage: R = ZpFP(7,4)
            sage: a = R(8)
            sage: a.residue(1)
            1
            sage: a.residue(2)
            8

            sage: K = QpFP(7,4)
            sage: a = K(8)
            sage: a.residue(1)
            1
            sage: a.residue(2)
            8
            sage: b = K(1/7)
            sage: b.residue()
            Traceback (most recent call last):
            ...
            ValueError: element must have nonnegative valuation in order to compute residue

        TESTS::

            sage: R = ZpFP(7,4)
            sage: a = R(8)
            sage: a.residue(0)
            0
            sage: a.residue(-1)
            Traceback (most recent call last):
            ...
            ValueError: cannot reduce modulo a negative power of p
            sage: a.residue(5)
            8

            sage: a.residue(field=True).parent()
            Finite Field of size 7
        """
        cdef Integer selfvalue, modulus
        cdef long aprec
        if not isinstance(absprec, Integer):
            absprec = Integer(absprec)
        if mpz_sgn((<Integer>absprec).value) < 0:
            raise ValueError("cannot reduce modulo a negative power of p")
        if self.ordp < 0:
            raise ValueError("element must have nonnegative valuation in order to compute residue")
        if field is None:
            field = (absprec == 1)
        elif field and absprec != 1:
            raise ValueError("field keyword may only be set at precision 1")
        if mpz_fits_slong_p((<Integer>absprec).value) == 0:
            raise ValueError("absolute precision does not fit in a long")
        aprec = mpz_get_si((<Integer>absprec).value)
        modulus = PY_NEW(Integer)
        mpz_set(modulus.value, self.prime_pow.pow_mpz_t_tmp(aprec))
        selfvalue = PY_NEW(Integer)
        if very_pos_val(self.ordp):
            mpz_set_ui(selfvalue.value, 0)
        else:
            # Need to do this better.
            mpz_mul(selfvalue.value, self.prime_pow.pow_mpz_t_tmp(self.ordp), self.unit)
        if field:
            from sage.rings.finite_rings.finite_field_constructor import GF
            return GF(self.parent().prime())(selfvalue)
        else:
            return Mod(selfvalue, modulus)

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

            \mathrm{self} = \sum_{i=1}^\infty a_i p^{2^i}

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
        cdef pAdicFloatingPointElement ans
        cdef Integer selfint = self.lift_c()

        if mpz_fits_slong_p(self.prime_pow.prime.value) == 0:
            raise NotImplementedError("the prime %s does not fit in a long" % self.prime_pow.prime)
        p = self.prime_pow.prime

        ans = self._new_c()
        ans.ordp = 0
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
          computing the logarithm. This attribute is passed to the :meth:`log`
          method; see :meth:`log` for more details about the possible
          algorithms

        .. NOTE::

            The function does not check that its argument ``self`` is
            the disk of convergence of ``exp``. If this assumption is not
            fulfilled the behaviour of the function is not specified.

        ALGORITHM:

        Solve the equation `\log(x) = \mathrm{self}` using the Newton scheme::

        .. MATH::

            x_{i+1} = x_i \cdot (1 + \mathrm{self} - \log(x_i)).

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
        cdef pAdicFloatingPointElement ans
        cdef Integer selfint = self.lift_c()

        if mpz_fits_slong_p(self.prime_pow.prime.value) == 0:
            raise NotImplementedError("the prime %s does not fit in a long" % self.prime_pow.prime)
        p = self.prime_pow.prime

        ans = self._new_c()
        ans.ordp = 0
        mpz_set_ui(ans.unit, 1)
        sig_on()
        if p == 2:
            padicexp_Newton(ans.unit, selfint.value, p, prec, 2, self.prime_pow.pow_mpz_t_tmp(prec))
        else:
            padicexp_Newton(ans.unit, selfint.value, p, prec, 1, self.prime_pow.pow_mpz_t_tmp(prec))
        sig_off()

        return ans
