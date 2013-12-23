"""
Base class for generic `p`-adic polynomials

This provides common functionality for all `p`-adic polynomials, such
as printing and factoring.

AUTHORS:

- Jeroen Demeyer (2013-11-22): initial version, split off from other
  files, made Polynomial_padic the common base class for all p-adic
  polynomials.

"""

#*****************************************************************************
#       Copyright (C) 2007 David Roe <roed.math@gmail.com>
#       Copyright (C) 2013 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.padics.precision_error import PrecisionError
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.infinity import infinity

class Polynomial_padic(Polynomial):
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen, construct)

    def _repr(self, name=None):
        r"""
        EXAMPLES::

            sage: R.<w> = PolynomialRing(Zp(5, prec=5, type = 'capped-abs', print_mode = 'val-unit'))
            sage: f = 24 + R(4/3)*w + w^4
            sage: f._repr()
            '(1 + O(5^5))*w^4 + (O(5^5))*w^3 + (O(5^5))*w^2 + (1043 + O(5^5))*w + (24 + O(5^5))'
            sage: f._repr(name='z')
            '(1 + O(5^5))*z^4 + (O(5^5))*z^3 + (O(5^5))*z^2 + (1043 + O(5^5))*z + (24 + O(5^5))'

        TESTS::

            sage: k = Qp(5,10)
            sage: R.<x> = k[]
            sage: f = R([k(0,-3), 0, k(0,-1)]); f
            (O(5^-1))*x^2 + (O(5^-3))
            sage: f + f
            (O(5^-1))*x^2 + (O(5^-3))

        AUTHOR:

        - David Roe (2007-03-03), based on Polynomial_generic_dense._repr()
        """
        s = ""
        coeffs = self.list()
        if name is None:
            name = self.parent().variable_name()
        for n in reversed(xrange(len(coeffs))):
            x = coeffs[n]
            if x.valuation() != infinity:
                if s:
                    s += " + "
                x = "(%s)"%repr(x)
                if n > 1:
                    var = "*%s^%s"%(name,n)
                elif n==1:
                    var = "*%s"%name
                else:
                    var = ""
                s += (x + var)
        return s or "0"

    def factor(self):
        """
        Return the factorization of this polynomial.

        EXAMPLES::

            sage: R.<t> = PolynomialRing(Qp(3,3,print_mode='terse',print_pos=False))
            sage: pol = t^8 - 1
            sage: for p,e in pol.factor(): print e, p
            1 (1 + O(3^3))*t + (1 + O(3^3))
            1 (1 + O(3^3))*t + (-1 + O(3^3))
            1 (1 + O(3^3))*t^2 + (5 + O(3^3))*t + (-1 + O(3^3))
            1 (1 + O(3^3))*t^2 + (-5 + O(3^3))*t + (-1 + O(3^3))
            1 (1 + O(3^3))*t^2 + (1 + O(3^3))
            sage: R.<t> = PolynomialRing(Qp(5,6,print_mode='terse',print_pos=False))
            sage: pol = 100 * (5*t - 1) * (t - 5)
            sage: pol
            (500 + O(5^9))*t^2 + (-2600 + O(5^8))*t + (500 + O(5^9))
            sage: pol.factor()
            (500 + O(5^9)) * ((1 + O(5^5))*t + (-1/5 + O(5^5))) * ((1 + O(5^6))*t + (-5 + O(5^6)))
            sage: pol.factor().value()
            (500 + O(5^8))*t^2 + (-2600 + O(5^8))*t + (500 + O(5^8))

        The same factorization over `\ZZ_p`. In this case, the "unit"
        part is a `p`-adic unit and the power of `p` is considered to be
        a factor::

            sage: R.<t> = PolynomialRing(Zp(5,6,print_mode='terse',print_pos=False))
            sage: pol = 100 * (5*t - 1) * (t - 5)
            sage: pol
            (500 + O(5^9))*t^2 + (-2600 + O(5^8))*t + (500 + O(5^9))
            sage: pol.factor()
            (4 + O(5^6)) * ((5 + O(5^7)))^2 * ((1 + O(5^6))*t + (-5 + O(5^6))) * ((5 + O(5^6))*t + (-1 + O(5^6)))
            sage: pol.factor().value()
            (500 + O(5^8))*t^2 + (-2600 + O(5^8))*t + (500 + O(5^8))

        In the following example, the discriminant is zero, so the `p`-adic
        factorization is not well defined::

            sage: factor(t^2)
            Traceback (most recent call last):
            ...
            PrecisionError: p-adic factorization not well-defined since the discriminant is zero up to the requestion p-adic precision

        More examples over `\ZZ_p`::

            sage: R.<w> = PolynomialRing(Zp(5, prec=6, type = 'capped-abs', print_mode = 'val-unit'))
            sage: f = w^5-1
            sage: f.factor()
            ((1 + O(5^6))*w + (3124 + O(5^6))) * ((1 + O(5^6))*w^4 + (12501 + O(5^6))*w^3 + (9376 + O(5^6))*w^2 + (6251 + O(5^6))*w + (3126 + O(5^6)))

        See :trac:`4038`::

            sage: E = EllipticCurve('37a1')
            sage: K =Qp(7,10)
            sage: EK = E.base_extend(K)
            sage: E = EllipticCurve('37a1')
            sage: K = Qp(7,10)
            sage: EK = E.base_extend(K)
            sage: g = EK.division_polynomial_0(3)
            sage: g.factor()
            (3 + O(7^10)) * ((1 + O(7^10))*x + (1 + 2*7 + 4*7^2 + 2*7^3 + 5*7^4 + 7^5 + 5*7^6 + 3*7^7 + 5*7^8 + 3*7^9 + O(7^10))) * ((1 + O(7^10))*x^3 + (6 + 4*7 + 2*7^2 + 4*7^3 + 7^4 + 5*7^5 + 7^6 + 3*7^7 + 7^8 + 3*7^9 + O(7^10))*x^2 + (6 + 3*7 + 5*7^2 + 2*7^4 + 7^5 + 7^6 + 2*7^8 + 3*7^9 + O(7^10))*x + (2 + 5*7 + 4*7^2 + 2*7^3 + 6*7^4 + 3*7^5 + 7^6 + 4*7^7 + O(7^10)))

        TESTS:

        Check that :trac:`13293` is fixed::

            sage: R.<T> = Qp(3)[]
            sage: f = 1926*T^2 + 312*T + 387
            sage: f.factor()
            (3^2 + 2*3^3 + 2*3^4 + 3^5 + 2*3^6 + O(3^22)) * ((1 + O(3^19))*T + (2*3^-1 + 3 + 3^2 + 2*3^5 + 2*3^6 + 2*3^7 + 3^8 + 3^9 + 2*3^11 + 3^15 + 3^17 + O(3^19))) * ((1 + O(3^20))*T + (2*3 + 3^2 + 3^3 + 3^5 + 2*3^6 + 2*3^7 + 3^8 + 3^10 + 3^11 + 2*3^12 + 2*3^14 + 2*3^15 + 2*3^17 + 2*3^18 + O(3^20)))
        """
        if self == 0:
            raise ArithmeticError("factorization of 0 not defined")
        # Scale self such that 0 is the lowest valuation
        # amongst the coefficients
        try:
            val = self.valuation(val_of_var=0)
        except TypeError:
            val = min([c.valuation() for c in self])
        self_normal = self / self.base_ring().uniformizer_pow(val)

        absprec = min([x.precision_absolute() for x in self_normal])
        if self_normal.discriminant().valuation() >= absprec:
            raise PrecisionError(
                "p-adic factorization not well-defined since the discriminant is zero up to the requestion p-adic precision")
        G = self_normal._pari_().factorpadic(self.base_ring().prime(), absprec)
        return _pari_padic_factorization_to_sage(G, self.parent(), self.leading_coefficient())

def _pari_padic_factorization_to_sage(G, R, leading_coeff):
    """
    Given a PARI factorization matrix `G` representing a factorization
    of some polynomial in the `p`-adic polynomial ring `R`,
    return the corresponding Sage factorization. All factors in `G`
    are assumed to have content 1 (this is how PARI returns its
    factorizations).

    INPUT:

    - ``G`` -- PARI factorization matrix, returned by ``factorpadic``.

    - ``R`` -- polynomial ring to be used as parent ring of the factors

    - ``leading_coeff`` -- leading coefficient of the polynomial which
      was factored. This can belong to any ring which can be coerced
      into ``R.base_ring()``.

    OUTPUT:

    - A Sage :class:`Factorization`.

    """
    B = R.base_ring()
    p = B.prime()
    leading_coeff = B(leading_coeff)

    pols = [R(f, absprec=f.padicprec(p)) for f in G[0]]
    exps = [int(e) for e in G[1]]

    # Determine unit part (which is discarded by PARI)
    if B.is_field():
        # When the base ring is a field, we normalize
        # the irreducible factors so they have leading
        # coefficient 1.
        for i in xrange(len(pols)):
            lc = pols[i].leading_coefficient()
            lc = lc.lift_to_precision()  # Ensure we don't lose precision
            pols[i] *= ~lc
    else:
        # When the base ring is not a field, we normalize
        # the irreducible factors so that the leading term
        # is a power of p.
        c, leading_coeff = leading_coeff.val_unit()
        for i in xrange(len(pols)):
            v, upart = pols[i].leading_coefficient().val_unit()
            upart = upart.lift_to_precision()  # Ensure we don't lose precision
            pols[i] *= ~upart
            c -= exps[i] * v
        if c:
            # Add factor p^c
            pols.append(R(p))
            exps.append(c)

    from sage.structure.factorization import Factorization
    return Factorization(zip(pols, exps), leading_coeff)
