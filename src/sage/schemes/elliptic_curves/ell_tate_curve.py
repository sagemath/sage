# sage_setup: distribution = sagemath-schemes
r"""
Tate's parametrisation of `p`-adic curves with multiplicative reduction

Let `E` be an elliptic curve defined over the `p`-adic numbers `\QQ_p`.
Suppose that `E` has multiplicative reduction, i.e. that the `j`-invariant
of `E` has negative valuation, say `n`. Then there exists a parameter
`q` in `\ZZ_p` of valuation `n` such that the points of `E` defined over
the algebraic closure `\bar{\QQ}_p` are in bijection with
`\bar{\QQ}_p^{\times}\,/\, q^{\ZZ}`. More precisely there exists
the series `s_4(q)` and `s_6(q)` such that the
`y^2+x y = x^3 + s_4(q) x+s_6(q)` curve is isomorphic to `E` over
`\bar{\QQ}_p` (or over `\QQ_p` if the reduction is *split* multiplicative).
There is a `p`-adic analytic map from
`\bar{\QQ}^{\times}_p` to this curve with kernel `q^{\ZZ}`.
Points of good reduction correspond to points of valuation
`0` in `\bar{\QQ}^{\times}_p`.

See chapter V of [Sil1994]_ for more details.

AUTHORS:

- Chris Wuthrich (23/05/2007): first version

- William Stein (2007-05-29): added some examples; editing.

- Chris Wuthrich (04/09): reformatted docstrings.
"""
######################################################################
#       Copyright (C) 2007 chris wuthrich
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
######################################################################
from sage.arith.functions import lcm as LCM
from sage.matrix.constructor import matrix
from sage.misc.functional import denominator, log
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.modular.modform.constructor import EisensteinForms, CuspForms
from sage.rings.integer_ring import ZZ
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.structure.richcmp import richcmp, richcmp_method
from sage.structure.sage_object import SageObject

lazy_import('sage.rings.padics.factory', 'Qp')


@richcmp_method
class TateCurve(SageObject):
    r"""
    Tate's `p`-adic uniformisation of an elliptic curve with
    multiplicative reduction.

    .. NOTE::

        Some of the methods of this Tate curve only work when the
        reduction is split multiplicative over `\QQ_p`.

    EXAMPLES::

        sage: e = EllipticCurve('130a1')
        sage: eq = e.tate_curve(5); eq
        5-adic Tate curve associated to the Elliptic Curve
         defined by y^2 + x*y + y = x^3 - 33*x + 68 over Rational Field
        sage: eq == loads(dumps(eq))
        True

    REFERENCES: [Sil1994]_
    """
    def __init__(self, E, p):
        r"""
        INPUT:

        - ``E`` -- an elliptic curve over the rational numbers

        - ``p`` -- a prime where `E` has multiplicative reduction,
          i.e., such that `j(E)` has negative valuation

        EXAMPLES::

            sage: e = EllipticCurve('130a1')
            sage: eq = e.tate_curve(2); eq
            2-adic Tate curve associated to the Elliptic Curve
             defined by y^2 + x*y + y = x^3 - 33*x + 68 over Rational Field
        """
        if not p.is_prime():
            raise ValueError("p (=%s) must be a prime" % p)
        if E.j_invariant().valuation(p) >= 0:
            raise ValueError("the elliptic curve must have multiplicative reduction at %s" % p)
        self._p = ZZ(p)
        self._E = E
        self._q = self.parameter()

    def __richcmp__(self, other, op):
        r"""
        Compare ``self`` and ``other``.

        TESTS::

            sage: E = EllipticCurve('35a')
            sage: eq5 = E.tate_curve(5)
            sage: eq7 = E.tate_curve(7)
            sage: eq7 == eq7
            True
            sage: eq7 == eq5
            False
        """
        if type(self) is not type(other):
            return NotImplemented

        return richcmp((self._E, self._p), (other._E, other._p), op)

    def _repr_(self):
        r"""
        Return print representation.

        EXAMPLES::

            sage: e = EllipticCurve('130a1')
            sage: eq = e.tate_curve(2)
            sage: eq._repr_()
            '2-adic Tate curve associated to the Elliptic Curve defined by y^2 + x*y + y = x^3 - 33*x + 68 over Rational Field'
        """
        return "%s-adic Tate curve associated to the %s" % (self._p, self._E)

    def original_curve(self):
        r"""
        Return the elliptic curve the Tate curve was constructed from.

        EXAMPLES::

            sage: eq = EllipticCurve('130a1').tate_curve(5)
            sage: eq.original_curve()
            Elliptic Curve defined by y^2 + x*y + y = x^3 - 33*x + 68
             over Rational Field
        """
        return self._E

    def prime(self):
        r"""
        Return the residual characteristic `p`.

        EXAMPLES::

            sage: eq = EllipticCurve('130a1').tate_curve(5)
            sage: eq.original_curve()
            Elliptic Curve defined by y^2 + x*y + y = x^3 - 33*x + 68
             over Rational Field
            sage: eq.prime()
            5
        """
        return self._p

    def parameter(self, prec=20):
        r"""
        Return the Tate parameter `q` such that the curve is isomorphic
        over the algebraic closure of `\QQ_p` to the curve
        `\QQ_p^{\times}/q^{\ZZ}`.

        INPUT:

        - ``prec`` -- the `p`-adic precision (default: 20)

        EXAMPLES::

            sage: eq = EllipticCurve('130a1').tate_curve(5)
            sage: eq.parameter(prec=5)
            3*5^3 + 3*5^4 + 2*5^5 + 2*5^6 + 3*5^7 + O(5^8)
        """
        qE = getattr(self, "_q", None)
        if qE and qE.precision_relative() >= prec:
            return Qp(self._p, prec=prec)(qE)

        E4 = EisensteinForms(weight=4).basis()[0]
        Delta = CuspForms(weight=12).basis()[0]
        j = (E4.q_expansion(prec + 3)) ** 3 / Delta.q_expansion(prec + 3)
        jinv = (1 / j).power_series()
        q_in_terms_of_jinv = jinv.reverse()
        R = Qp(self._p, prec=prec)
        qE = q_in_terms_of_jinv(R(1 / self._E.j_invariant()))
        self._q = qE
        return qE

    def __sk(self, k, prec):
        q = self.parameter(prec=prec)
        return sum(n ** k * q ** n / (1 - q ** n)
                   for n in range(1, prec + 1))

    def __delta(self, prec):
        q = self.parameter(prec=prec)
        return q * prod([(1 - q**n)**24
                        for n in range(1, prec + 1)])

    def curve(self, prec=20):
        r"""
        Return the `p`-adic elliptic curve of the form
        `y^2+x y = x^3 + s_4 x+s_6`.

        This curve with split multiplicative reduction is isomorphic
        to the given curve over the algebraic closure of `\QQ_p`.

        INPUT:

        - ``prec`` -- the `p`-adic precision (default: 20)

        EXAMPLES::

            sage: eq = EllipticCurve('130a1').tate_curve(5)
            sage: eq.curve(prec=5)
            Elliptic Curve defined by y^2 + (1+O(5^5))*x*y =
             x^3 + (2*5^4+5^5+2*5^6+5^7+3*5^8+O(5^9))*x + (2*5^3+5^4+2*5^5+5^7+O(5^8))
             over 5-adic Field with capped relative precision 5
        """
        Eq = getattr(self, "__curve", None)
        if Eq and Eq.a6().precision_relative() >= prec:
            return Eq.change_ring(Qp(self._p, prec))

        qE = self.parameter(prec=prec)
        precp = prec + 2
        tate_a4 = -5 * self.__sk(3, precp)
        tate_a6 = (tate_a4 - 7 * self.__sk(5, precp)) / 12
        R = qE.parent()
        Eq = EllipticCurve([R.one(), R.zero(), R.zero(), R(tate_a4), R(tate_a6)])
        self.__curve = Eq
        return Eq

    def _Csquare(self, prec=20):
        r"""
        Return the square of the constant `C` such that the canonical
        Néron differential `\omega` and the canonical differential
        `\frac{du}{u}` on `\QQ^{\times}/q^{\ZZ}` are linked by `\omega
        = C \frac{du}{u}`.

        This constant is only a square in `\QQ_p` if the curve has split
        multiplicative reduction.

        INPUT:

        - ``prec`` -- the `p`-adic precision (default: 20)

        EXAMPLES::

            sage: eq = EllipticCurve('130a1').tate_curve(5)
            sage: eq._Csquare(prec=5)
            4 + 2*5^2 + 2*5^4 + O(5^5)
        """
        Csq = getattr(self, "__csquare", None)
        if Csq and Csq.precision_relative() >= prec:
            return Csq

        Eq = self.curve(prec=prec)
        tateCsquare = Eq.c6() * self._E.c4() / Eq.c4() / self._E.c6()
        self.__Csquare = tateCsquare
        return tateCsquare

    def E2(self, prec=20):
        r"""
        Return the value of the `p`-adic Eisenstein series of weight 2
        evaluated on the elliptic curve having split multiplicative
        reduction.

        INPUT:

        - ``prec`` -- the `p`-adic precision (default: 20)

        EXAMPLES::

            sage: eq = EllipticCurve('130a1').tate_curve(5)
            sage: eq.E2(prec=10)
            4 + 2*5^2 + 2*5^3 + 5^4 + 2*5^5 + 5^7 + 5^8 + 2*5^9 + O(5^10)

            sage: T = EllipticCurve('14').tate_curve(7)
            sage: T.E2(30)
            2 + 4*7 + 7^2 + 3*7^3 + 6*7^4 + 5*7^5 + 2*7^6 + 7^7 + 5*7^8 + 6*7^9 + 5*7^10 + 2*7^11 + 6*7^12 + 4*7^13 + 3*7^15 + 5*7^16 + 4*7^17 + 4*7^18 + 2*7^20 + 7^21 + 5*7^22 + 4*7^23 + 4*7^24 + 3*7^25 + 6*7^26 + 3*7^27 + 6*7^28 + O(7^30)
        """
        p = self._p
        Csq = self._Csquare(prec=prec)
        qE = self.parameter(prec=prec)
        n = qE.valuation()
        R = Qp(p, prec)
        e2 = Csq * (1 - 24 * sum(qE**i / (1 - qE**i)**2
                                 for i in range(1, prec // n + 5)))
        return R(e2)

    def is_split(self) -> bool:
        r"""
        Return ``True`` if the given elliptic curve has split multiplicative reduction.

        EXAMPLES::

            sage: eq = EllipticCurve('130a1').tate_curve(5)
            sage: eq.is_split()
            True

            sage: eq = EllipticCurve('37a1').tate_curve(37)
            sage: eq.is_split()
            False
        """
        return self._Csquare().is_square()

    def parametrisation_onto_tate_curve(self, u, prec=None):
        r"""
        Given an element `u` in `\QQ_p^{\times}`, this computes its image on the Tate curve
        under the `p`-adic uniformisation of `E`.

        INPUT:

        - ``u`` -- a nonzero `p`-adic number

        - ``prec`` -- the `p`-adic precision; default is the relative precision
          of ``u``, otherwise 20

        EXAMPLES::

            sage: eq = EllipticCurve('130a1').tate_curve(5)
            sage: eq.parametrisation_onto_tate_curve(1+5+5^2+O(5^10), prec=10)
            (5^-2 + 4*5^-1 + 1 + 2*5 + 3*5^2 + 2*5^5 + 3*5^6 + O(5^7)
             : 4*5^-3 + 2*5^-1 + 4 + 2*5 + 3*5^4 + 2*5^5 + O(5^6) : 1 + O(5^10))
            sage: eq.parametrisation_onto_tate_curve(1+5+5^2+O(5^10))
            (5^-2 + 4*5^-1 + 1 + 2*5 + 3*5^2 + 2*5^5 + 3*5^6 + O(5^7)
             : 4*5^-3 + 2*5^-1 + 4 + 2*5 + 3*5^4 + 2*5^5 + O(5^6) : 1 + O(5^10))
            sage: eq.parametrisation_onto_tate_curve(1+5+5^2+O(5^10), prec=20)
            Traceback (most recent call last):
            ...
            ValueError: requested more precision than the precision of u
        """
        if prec is None:
            prec = getattr(u, "precision_relative", lambda: 20)()
        u = Qp(self._p, prec)(u)
        if prec > u.precision_relative():
            raise ValueError("requested more precision than the precision of u")
        if u == 1:
            return self.curve(prec=prec)(0)

        q = self.parameter(prec=prec)
        un = u * q ** (-(u.valuation() // q.valuation()))

        precn = (prec // q.valuation()) + 4

        # formulas in Silverman II (Advanced Topics in the Arithmetic
        # of Elliptic curves, p. 425)

        powers_of_q = [(n, q**n) for n in range(1, precn)]
        xx = un / (1 - un)**2 + sum(qn * un / (1 - qn * un)**2 +
                                    qn / un / (1 - qn / un)**2 -
                                    2 * qn / (1 - qn)**2
                                    for n, qn in powers_of_q)

        yy = un**2 / (1 - un)**3 + sum(qn**2 * un**2 / (1 - qn * un)**3 -
                                       qn / un / (1 - qn / un)**3 +
                                       qn / (1 - qn)**2
                                       for n, qn in powers_of_q)

        return self.curve(prec=prec)([xx, yy])

    # From here on all functions need that the curve has split
    # multiplicative reduction.

    def L_invariant(self, prec=20):
        r"""
        Return the *mysterious* `\mathcal{L}`-invariant associated
        to an elliptic curve with split multiplicative reduction.

        One instance where this constant appears is in the exceptional
        case of the `p`-adic Birch and Swinnerton-Dyer conjecture as
        formulated in [MTT1986]_. See [Col2004]_ for a detailed discussion.

        INPUT:

        - ``prec`` -- the `p`-adic precision (default: 20)

        EXAMPLES::

            sage: eq = EllipticCurve('130a1').tate_curve(5)
            sage: eq.L_invariant(prec=10)
            5^3 + 4*5^4 + 2*5^5 + 2*5^6 + 2*5^7 + 3*5^8 + 5^9 + O(5^10)
        """
        if not self.is_split():
            raise RuntimeError("the curve must have split multiplicative "
                               "reduction")
        qE = self.parameter(prec=prec)
        n = qE.valuation()
        u = qE / self._p ** n
        # the p-adic logarithm of Iwasawa normalised by log(p) = 0
        return log(u) / n

    def _isomorphism(self, prec=20):
        r"""
        Return the isomorphism between ``self.curve()`` and the given
        curve in the form of a list ``[u,r,s,t]`` of `p`-adic numbers.

        For this to exist the given curve has to have split
        multiplicative reduction over `\QQ_p`.

        More precisely, if `E` has coordinates `x` and `y` and the Tate
        curve has coordinates `X`, `Y` with `Y^2 + XY = X^3 + s_4 X +s_6`
        then `X = u^2 x +r` and `Y = u^3 y +s u^2 x +t`.

        INPUT:

        - ``prec`` -- the `p`-adic precision (default: 20)

        EXAMPLES::

            sage: eq = EllipticCurve('130a1').tate_curve(5)
            sage: eq._isomorphism(prec=5)
            [2 + 3*5^2 + 2*5^3 + 4*5^4 + O(5^5),
             4 + 3*5 + 4*5^2 + 2*5^3 + O(5^5),
             3 + 2*5 + 5^2 + 5^3 + 2*5^4 + O(5^5),
             2 + 5 + 3*5^2 + 5^3 + 5^4 + O(5^5)]
        """
        if not self.is_split():
            raise RuntimeError("the curve must have split multiplicative "
                               "reduction")
        C = self._Csquare(prec=prec + 4).sqrt()
        R = Qp(self._p, prec)
        C = R(C)
        s = (C * R(self._E.a1()) - R.one()) / R(2)
        r = (C ** 2 * R(self._E.a2()) + s + s ** 2) / R(3)
        t = (C ** 3 * R(self._E.a3()) - r) / R(2)
        return [C, r, s, t]

    def _inverse_isomorphism(self, prec=20):
        r"""
        Return the isomorphism between the given curve and
        ``self.curve()`` in the form of a list ``[u,r,s,t]`` of
        `p`-adic numbers.

        For this to exist the given curve has to have split
        multiplicative reduction over `\QQ_p`.

        More precisely, if `E` has coordinates `x` and `y` and the Tate
        curve has coordinates `X`, `Y` with `Y^2 + XY = X^3 + s_4 X +s_6`
        then `x = u^2 X +r` and `y = u^3 Y +s u^2 X +t`.

        INPUT:

        - ``prec`` -- the `p`-adic precision (default: 20)

        EXAMPLES::

            sage: eq = EllipticCurve('130a1').tate_curve(5)
            sage: eq._inverse_isomorphism(prec=5)
            [3 + 2*5 + 3*5^3 + O(5^5), 4 + 2*5 + 4*5^3 + 3*5^4 + O(5^5),
             1 + 5 + 4*5^3 + 2*5^4 + O(5^5), 5 + 2*5^2 + 3*5^4 + O(5^5)]
        """
        if not self.is_split():
            raise RuntimeError("the curve must have split multiplicative "
                               "reduction")
        u, r, s, t = self._isomorphism(prec=prec)
        return [1 / u, -r / u ** 2, -s / u, (r * s - t) / u ** 3]

    def lift(self, P, prec=20):
        r"""
        Given a point `P` in the formal group of the elliptic curve `E` with split multiplicative reduction,
        this produces an element `u` in `\QQ_p^{\times}` mapped to the point `P` by the Tate parametrisation.
        The algorithm return the unique such element in `1+p\ZZ_p`.

        INPUT:

        - ``P`` -- a point on the elliptic curve

        - ``prec`` -- the `p`-adic precision (default: 20)

        EXAMPLES::

            sage: e = EllipticCurve('130a1')
            sage: eq = e.tate_curve(5)
            sage: P = e([-6,10])
            sage: l = eq.lift(12*P, prec=10); l
            1 + 4*5 + 5^3 + 5^4 + 4*5^5 + 5^6 + 5^7 + 4*5^8 + 5^9 + O(5^10)

        Now we map the lift l back and check that it is indeed right.::

            sage: eq.parametrisation_onto_original_curve(l)
            (4*5^-2 + 2*5^-1 + 4*5 + 3*5^3 + 5^4 + 2*5^5 + 4*5^6 + O(5^7)
             : 2*5^-3 + 5^-1 + 4 + 4*5 + 5^2 + 3*5^3 + 4*5^4 + O(5^6) : 1 + O(5^10))
            sage: e5 = e.change_ring(Qp(5,9))
            sage: e5(12*P)
            (4*5^-2 + 2*5^-1 + 4*5 + 3*5^3 + 5^4 + 2*5^5 + 4*5^6 + O(5^7)
             : 2*5^-3 + 5^-1 + 4 + 4*5 + 5^2 + 3*5^3 + 4*5^4 + O(5^6) : 1 + O(5^9))
        """
        p = self._p
        R = Qp(self._p, prec)
        if not self._E == P.curve():
            raise ValueError("the point must lie on the original curve.")
        if not self.is_split():
            raise ValueError("the curve must have split multiplicative reduction.")
        if P.is_zero():
            return R.one()
        if P[0].valuation(p) >= 0:
            raise ValueError("the point must lie in the formal group.")

        Eq = self.curve(prec=prec)
        C, r, s, t = self._isomorphism(prec=prec)
        xx = r + C ** 2 * P[0]
        yy = t + s * C ** 2 * P[0] + C ** 3 * P[1]
        assert Eq.defining_polynomial()(xx, yy, 1) == 0, f"bug: point ({xx}, {yy}) does not lie on the curve {Eq}"
        tt = -xx / yy
        eqhat = Eq.formal()
        eqlog = eqhat.log(prec + 3)
        z = eqlog(tt)
        u = ZZ.one()
        fac = ZZ.one()
        for i in range(1, 2 * prec + 1):
            fac *= i
            u += z ** i / fac
        return u

    def parametrisation_onto_original_curve(self, u, prec=None):
        r"""
        Given an element `u` in `\QQ_p^{\times}`, this computes its image on the original curve
        under the `p`-adic uniformisation of `E`.

        INPUT:

        - ``u`` -- a nonzero `p`-adic number

        - ``prec`` -- the `p`-adic precision; default is the relative precision
          of ``u``, otherwise 20

        EXAMPLES::

            sage: eq = EllipticCurve('130a1').tate_curve(5)
            sage: eq.parametrisation_onto_original_curve(1+5+5^2+O(5^10))
            (4*5^-2 + 4*5^-1 + 4 + 2*5^3 + 3*5^4 + 2*5^6 + O(5^7) :
             3*5^-3 + 5^-2 + 4*5^-1 + 1 + 4*5 + 5^2 + 3*5^5 + O(5^6) :
             1 + O(5^10))
            sage: eq.parametrisation_onto_original_curve(1+5+5^2+O(5^10), prec=20)
            Traceback (most recent call last):
            ...
            ValueError: requested more precision than the precision of u

        Here is how one gets a 4-torsion point on `E` over `\QQ_5`::

            sage: R = Qp(5,30)
            sage: i = R(-1).sqrt()
            sage: T = eq.parametrisation_onto_original_curve(i, prec=30); T
            (2 + 3*5 + 4*5^2 + 2*5^3 + 5^4 + 4*5^5 + 2*5^7 + 5^8 + 5^9 + 5^12 + 3*5^13 + 3*5^14 + 5^15 + 4*5^17 + 5^18 + 3*5^19 + 2*5^20 + 4*5^21 + 5^22 + 3*5^23 + 3*5^24 + 4*5^25 + 3*5^26 + 3*5^27 + 3*5^28 + 3*5^29 + O(5^30) : 3*5 + 5^2 + 5^4 + 3*5^5 + 3*5^7 + 2*5^8 + 4*5^9 + 5^10 + 2*5^11 + 4*5^13 + 2*5^14 + 4*5^15 + 4*5^16 + 3*5^17 + 2*5^18 + 4*5^20 + 2*5^21 + 2*5^22 + 4*5^23 + 4*5^24 + 4*5^25 + 5^26 + 3*5^27 + 2*5^28 + O(5^30) : 1 + O(5^30))
            sage: 4*T
            (0 : 1 + O(5^30) : 0)
        """
        if not self.is_split():
            raise ValueError("the curve must have split multiplicative "
                             "reduction.")
        if prec is None:
            prec = getattr(u, "precision_relative", lambda: 20)()

        P = self.parametrisation_onto_tate_curve(u, prec=prec)
        C, r, s, t = self._inverse_isomorphism(prec=prec)
        xx = r + C ** 2 * P[0]
        yy = t + s * C ** 2 * P[0] + C ** 3 * P[1]
        R = Qp(self._p, prec)
        E_over_Qp = self._E.base_extend(R)
        return E_over_Qp([xx, yy])

    def __padic_sigma_square(self, u, prec):
        q = self.parameter(prec=prec)
        return (u - 1)**2 / u * prod([((1 - q**n * u) * (1 - q**n / u) /
                                       (1 - q**n)**2)**2
                                      for n in range(1, prec + 1)])

    # the following functions are rather functions of the global curve
    # than the local curve
    # we use the same names as for elliptic curves over rationals.

    def padic_height(self, prec=20):
        r"""
        Return the canonical `p`-adic height function on the original curve.

        INPUT:

        - ``prec`` -- the `p`-adic precision (default: 20)

        OUTPUT: a function that can be evaluated on rational points of `E`

        EXAMPLES::

            sage: e = EllipticCurve('130a1')
            sage: eq = e.tate_curve(5)
            sage: h = eq.padic_height(prec=10)
            sage: P = e.gens()[0]
            sage: h(P)
            2*5^-1 + 1 + 2*5 + 2*5^2 + 3*5^3 + 3*5^6 + 5^7 + O(5^9)

        Check that it is a quadratic function::

            sage: h(3*P)-3^2*h(P)
            O(5^9)
        """
        if not self.is_split():
            raise NotImplementedError("the p-adic height is not implemented "
                                      "for non-split multiplicative reduction.")

        p = self._p

        # we will have to do it properly with David Harvey's _multiply_point(E, R, Q)
        n = LCM(self._E.tamagawa_numbers()) * (p - 1)

        # this function is a closure, I don't see how to doctest it (PZ)
        def _height(P, check=True):
            if check:
                assert P.curve() == self._E, "the point P must lie on the curve from which the height function was created"
            Q = n * P
            cQ = denominator(Q[0])
            q = self.parameter(prec=prec)
            nn = q.valuation()
            precp = prec + nn + 2
            uQ = self.lift(Q, prec=precp)
            si = self.__padic_sigma_square(uQ, prec=precp)
            q = self.parameter(prec=precp)
            nn = q.valuation()
            qEu = q / p ** nn
            res = -(log(si * self._Csquare(prec=precp) / cQ) + log(uQ)**2 / log(qEu)) / n**2
            R = Qp(self._p, prec)
            return R(res)

        return _height

    def padic_regulator(self, prec=20):
        r"""
        Compute the canonical `p`-adic regulator on the extended
        Mordell-Weil group as in [MTT1986]_
        (with the correction of [Wer1998]_ and sign convention in [SW2013]_.)

        The `p`-adic Birch and Swinnerton-Dyer conjecture predicts
        that this value appears in the formula for the leading term of
        the `p`-adic `L`-function.

        INPUT:

        - ``prec`` -- the `p`-adic precision (default: 20)

        EXAMPLES::

            sage: eq = EllipticCurve('130a1').tate_curve(5)
            sage: eq.padic_regulator()
            2*5^-1 + 1 + 2*5 + 2*5^2 + 3*5^3 + 3*5^6 + 5^7 + 3*5^9 + 3*5^10 + 3*5^12 + 4*5^13 + 3*5^15 + 2*5^16 + 3*5^18 + 4*5^19 +  4*5^20 + 3*5^21 + 4*5^22 + O(5^23)
        """
        prec = prec + 4

        K = Qp(self._p, prec=prec)
        rank = self._E.rank()
        if rank == 0:
            return K.one()

        if not self.is_split():
            raise NotImplementedError("the p-adic regulator is not implemented "
                                      "for non-split multiplicative reduction.")

        basis = self._E.gens()
        M = matrix(K, rank, rank, 0)

        height = self.padic_height(prec=prec)
        point_height = [height(P) for P in basis]
        for i in range(rank):
            for j in range(i + 1, rank):
                M[i, j] = M[j, i] = (- point_height[i] - point_height[j] + height(basis[i] + basis[j])) / 2
        for i in range(rank):
            M[i, i] = point_height[i]

        return M.determinant()
