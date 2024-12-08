r"""
Asymptotics of Multivariate Generating Series

Let `F(x) = \sum_{\nu \in \NN^d} F_{\nu} x^\nu` be a multivariate power series
with complex coefficients that converges in a neighborhood of the origin.
Assume that `F = G/H` for some functions `G` and `H` holomorphic in a
neighborhood of the origin. Assume also that `H` is a polynomial.

This computes asymptotics for the coefficients `F_{r \alpha}` as `r \to \infty`
with `r \alpha \in \NN^d` for `\alpha` in a permissible subset of `d`-tuples of
positive reals. More specifically, it computes arbitrary terms of the
asymptotic expansion for `F_{r \alpha}` when the asymptotics are controlled by
a strictly minimal multiple point of the algebraic variety `H = 0`.

The algorithms and formulas implemented here come from [RW2008]_
and [RW2012]_. For a general reference take a look in the book [PW2013]_.


Introductory Examples
=====================

::

    sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing

A univariate smooth point example::

    sage: R.<x> = PolynomialRing(QQ)
    sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
    sage: H = (x - 1/2)^3
    sage: Hfac = H.factor()
    sage: G = -1/(x + 3)/Hfac.unit()
    sage: F = FFPD(G, Hfac)
    sage: F
    (-1/(x + 3), [(x - 1/2, 3)])
    sage: alpha = [1]
    sage: decomp = F.asymptotic_decomposition(alpha)
    sage: decomp
    (0, []) +
    (-1/2*r^2*(x^2/(x^5 + 9*x^4 + 27*x^3 + 27*x^2)
     + 6*x/(x^5 + 9*x^4 + 27*x^3 + 27*x^2)
     + 9/(x^5 + 9*x^4 + 27*x^3 + 27*x^2))
     - 1/2*r*(5*x^2/(x^5 + 9*x^4 + 27*x^3 + 27*x^2)
     + 24*x/(x^5 + 9*x^4 + 27*x^3 + 27*x^2)
     + 27/(x^5 + 9*x^4 + 27*x^3 + 27*x^2))
     - 3*x^2/(x^5 + 9*x^4 + 27*x^3 + 27*x^2)
     - 9*x/(x^5 + 9*x^4 + 27*x^3 + 27*x^2)
     - 9/(x^5 + 9*x^4 + 27*x^3 + 27*x^2),
     [(x - 1/2, 1)])
    sage: F1 = decomp[1]
    sage: p = {x: 1/2}
    sage: asy = F1.asymptotics(p, alpha, 3)
    sage: asy
    (8/343*(49*r^2 + 161*r + 114)*2^r, 2, 8/7*r^2 + 184/49*r + 912/343)
    sage: F.relative_error(asy[0], alpha, [1, 2, 4, 8, 16], asy[1])
    [((1,), 7.555555556, [7.556851312], [-0.0001714971672]),
     ((2,), 14.74074074, [14.74052478], [0.00001465051901]),
     ((4,), 35.96502058, [35.96501458], [1.667911934e-7]),
     ((8,), 105.8425656, [105.8425656], [4.399565380e-11]),
     ((16,), 355.3119534, [355.3119534], [0.0000000000])]

Another smooth point example (Example 5.4 of [RW2008]_)::

    sage: R.<x,y> = PolynomialRing(QQ)
    sage: FFPD = FractionWithFactoredDenominatorRing(R)
    sage: q = 1/2
    sage: qq = q.denominator()
    sage: H = 1 - q*x + q*x*y - x^2*y
    sage: Hfac = H.factor()
    sage: G = (1 - q*x)/Hfac.unit()
    sage: F = FFPD(G, Hfac)
    sage: alpha = list(qq*vector([2, 1 - q]))
    sage: alpha
    [4, 1]
    sage: I = F.smooth_critical_ideal(alpha)
    sage: I
    Ideal (y^2 - 2*y + 1, x + 1/4*y - 5/4) of
     Multivariate Polynomial Ring in x, y over Rational Field
    sage: s = solve([SR(z) for z in I.gens()],
    ....:           [SR(z) for z in R.gens()], solution_dict=true)
    sage: s == [{SR(x): 1, SR(y): 1}]
    True
    sage: p = s[0]
    sage: asy = F.asymptotics(p, alpha, 1, verbose=True)
    Creating auxiliary functions...
    Computing derivatives of auxiliary functions...
    Computing derivatives of more auxiliary functions...
    Computing second order differential operator actions...
    sage: asy
    (1/24*2^(2/3)*(sqrt(3) + 4/(sqrt(3) + I) + I)*gamma(1/3)/(pi*r^(1/3)),
     1,
     1/24*2^(2/3)*(sqrt(3) + 4/(sqrt(3) + I) + I)*gamma(1/3)/(pi*r^(1/3)))
    sage: r = SR('r')
    sage: tuple((a*r^(1/3)).full_simplify() / r^(1/3) for a in asy)  # make nicer coefficients
    (1/12*sqrt(3)*2^(2/3)*gamma(1/3)/(pi*r^(1/3)),
     1,
     1/12*sqrt(3)*2^(2/3)*gamma(1/3)/(pi*r^(1/3)))
    sage: F.relative_error(asy[0], alpha, [1, 2, 4, 8, 16], asy[1])
    [((4, 1), 0.1875000000, [0.1953794675...], [-0.042023826...]),
     ((8, 2), 0.1523437500, [0.1550727862...], [-0.017913673...]),
     ((16, 4), 0.1221771240, [0.1230813519...], [-0.0074009592...]),
     ((32, 8), 0.09739671811, [0.09768973377...], [-0.0030084757...]),
     ((64, 16), 0.07744253816, [0.07753639308...], [-0.0012119297...])]

A multiple point example (Example 6.5 of [RW2012]_)::

    sage: R.<x,y> = PolynomialRing(QQ)
    sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
    sage: H = (1 - 2*x - y)**2 * (1 - x - 2*y)**2
    sage: Hfac = H.factor()
    sage: G = 1/Hfac.unit()
    sage: F = FFPD(G, Hfac)
    sage: F
    (1, [(x + 2*y - 1, 2), (2*x + y - 1, 2)])
    sage: I = F.singular_ideal()
    sage: I
    Ideal (x - 1/3, y - 1/3) of
     Multivariate Polynomial Ring in x, y over Rational Field
    sage: p = {x: 1/3, y: 1/3}
    sage: F.is_convenient_multiple_point(p)
    (True, 'convenient in variables [x, y]')
    sage: alpha = (var('a'), var('b'))
    sage: decomp =  F.asymptotic_decomposition(alpha); decomp
    (0, []) +
    (-1/9*r^2*(2*a^2/x^2 + 2*b^2/y^2 - 5*a*b/(x*y))
     - 1/9*r*(6*a/x^2 + 6*b/y^2 - 5*a/(x*y) - 5*b/(x*y))
     - 4/9/x^2 - 4/9/y^2 + 5/9/(x*y),
     [(x + 2*y - 1, 1), (2*x + y - 1, 1)])
    sage: F1 = decomp[1]
    sage: F1.asymptotics(p, alpha, 2)
    (-3*((2*a^2 - 5*a*b + 2*b^2)*r^2 + (a + b)*r + 3)*(1/((1/3)^a*(1/3)^b))^r,
     1/((1/3)^a*(1/3)^b), -3*(2*a^2 - 5*a*b + 2*b^2)*r^2 - 3*(a + b)*r - 9)
    sage: alpha = [4, 3]
    sage: decomp =  F.asymptotic_decomposition(alpha)
    sage: F1 = decomp[1]
    sage: asy = F1.asymptotics(p, alpha, 2)
    sage: asy
    (3*(10*r^2 - 7*r - 3)*2187^r, 2187, 30*r^2 - 21*r - 9)
    sage: F.relative_error(asy[0], alpha, [1, 2, 4, 8], asy[1])
    [((4, 3), 30.72702332, [0.0000000000], [1.000000000]),
     ((8, 6), 111.9315678, [69.00000000], [0.3835519207]),
     ((16, 12), 442.7813138, [387.0000000], [0.1259793763]),
     ((32, 24), 1799.879232, [1743.000000], [0.03160169385])]

TESTS::

    sage: R.<x,y> = PolynomialRing(QQ)
    sage: FFPD = FractionWithFactoredDenominatorRing(R)
    sage: H = (1 - 2*x - y) * (1 - x - 2*y)
    sage: G = 1
    sage: Hfac = H.factor()
    sage: G = G / Hfac.unit()
    sage: F = FFPD(G, Hfac); F
    (1, [(x + 2*y - 1, 1), (2*x + y - 1, 1)])
    sage: p = {x: 1, y: 1}
    sage: alpha = [1, 1]
    sage: F.asymptotics(p, alpha, 1)
    (1/3, 1, 1/3)

::

    sage: R.<x,y,t> = PolynomialRing(QQ)
    sage: FFPD = FractionWithFactoredDenominatorRing(R)
    sage: H = (1 - y) * (1 + x^2) * (1 - t*(1 + x^2 + x*y^2))
    sage: G = (1 + x) * (1 + x^2 - x*y^2)
    sage: Hfac = H.factor()
    sage: G = G / Hfac.unit()
    sage: F = FFPD(G, Hfac); F
    (-x^2*y^2 + x^3 - x*y^2 + x^2 + x + 1,
     [(y - 1, 1), (x^2 + 1, 1), (x*y^2*t + x^2*t + t - 1, 1)])
    sage: p = {x: 1, y: 1, t: 1/3}
    sage: alpha = [1, 1, 1]
    sage: F.asymptotics_multiple(p, alpha, 1, var('r'))  # not tested - see #19989


Various
=======

AUTHORS:

- Alexander Raichev (2008)
- Daniel Krenn (2014, 2016)


Classes and Methods
===================
"""
# ****************************************************************************
#       Copyright (C) 2008 Alexander Raichev <tortoise.said@gmail.com>
#       Copyright (C) 2014, 2016 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from functools import total_ordering
from itertools import combinations_with_replacement

from sage.categories.rings import Rings
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.element import RingElement
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp_by_eq_and_lt
from sage.structure.unique_representation import UniqueRepresentation

lazy_import("sage.calculus.var", "var")
lazy_import("sage.calculus.functional", "diff")
lazy_import("sage.symbolic.ring", "SR")


@total_ordering
class FractionWithFactoredDenominator(RingElement):
    r"""
    This element represents a fraction with a factored polynomial
    denominator. See also its parent
    :class:`FractionWithFactoredDenominatorRing` for details.

    Represents a fraction with factored polynomial denominator (FFPD)
    `p/(q_1^{e_1} \cdots q_n^{e_n})` by storing the parts `p` and
    `[(q_1, e_1), \ldots, (q_n, e_n)]`.
    Here `q_1, \ldots, q_n` are elements of a 0- or multi-variate factorial
    polynomial ring `R` , `q_1, \ldots, q_n` are distinct irreducible elements
    of `R` , `e_1, \ldots, e_n` are positive integers, and `p` is a function
    of the indeterminates of `R` (e.g., a Sage symbolic expression). An
    element `r` with no polynomial denominator is represented as ``(r, [])``.

    INPUT:

    - ``numerator`` -- an element `p`; this can be of any ring from which
      parent's base has coercion in
    - ``denominator_factored`` -- list of the form
      `[(q_1, e_1), \ldots, (q_n, e_n)]`, where the `q_1, \ldots, q_n` are
      distinct irreducible elements of `R` and the `e_i` are positive
      integers
    - ``reduce`` -- (optional) if ``True``, then represent
      `p/(q_1^{e_1} \cdots q_n^{e_n})` in lowest terms, otherwise
      this won't attempt to divide `p` by any of the `q_i`

    OUTPUT:

    An element representing the rational expression
    `p/(q_1^{e_1} \cdots q_n^{e_n})`.

    EXAMPLES::

        sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
        sage: R.<x,y> = PolynomialRing(QQ)
        sage: FFPD = FractionWithFactoredDenominatorRing(R)
        sage: df = [x, 1], [y, 1], [x*y+1, 1]
        sage: f = FFPD(x, df)
        sage: f
        (1, [(y, 1), (x*y + 1, 1)])
        sage: ff = FFPD(x, df, reduce=False)
        sage: ff
        (x, [(y, 1), (x, 1), (x*y + 1, 1)])

        sage: f = FFPD(x + y, [(x + y, 1)])
        sage: f
        (1, [])

    ::

        sage: R.<x> = PolynomialRing(QQ)
        sage: FFPD = FractionWithFactoredDenominatorRing(R)
        sage: f = 5*x^3 + 1/x + 1/(x-1) + 1/(3*x^2 + 1)
        sage: FFPD(f)
        (5*x^7 - 5*x^6 + 5/3*x^5 - 5/3*x^4 + 2*x^3 - 2/3*x^2 + 1/3*x - 1/3,
        [(x - 1, 1), (x, 1), (x^2 + 1/3, 1)])

    ::

        sage: R.<x,y> = PolynomialRing(QQ)
        sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
        sage: f = 2*y/(5*(x^3 - 1)*(y + 1))
        sage: FFPD(f)
        (2/5*y, [(y + 1, 1), (x - 1, 1), (x^2 + x + 1, 1)])

        sage: p = 1/x^2
        sage: q = 3*x**2*y
        sage: qs = q.factor()
        sage: f = FFPD(p/qs.unit(), qs)
        sage: f
        (1/3/x^2, [(y, 1), (x, 2)])

        sage: f = FFPD(cos(x)*x*y^2, [(x, 2), (y, 1)])
        sage: f
        (x*y^2*cos(x), [(y, 1), (x, 2)])

        sage: G = exp(x + y)
        sage: H = (1 - 2*x - y) * (1 - x - 2*y)
        sage: a = FFPD(G/H)
        sage: a
        (e^(x + y), [(x + 2*y - 1, 1), (2*x + y - 1, 1)])
        sage: a.denominator_ring
        Multivariate Polynomial Ring in x, y over Rational Field
        sage: b = FFPD(G, H.factor())
        sage: b
        (e^(x + y), [(x + 2*y - 1, 1), (2*x + y - 1, 1)])
        sage: b.denominator_ring
        Multivariate Polynomial Ring in x, y over Rational Field

    Singular throws a 'not implemented' error when trying to factor in
    a multivariate polynomial ring over an inexact field::

        sage: R.<x,y> = PolynomialRing(CC)
        sage: FFPD = FractionWithFactoredDenominatorRing(R)
        sage: f = (x + 1)/(x*y*(x*y + 1)^2)
        sage: FFPD(f)
        Traceback (most recent call last):
        ...
        TypeError: Singular error:
           ? not implemented
           ? error occurred in or before STDIN line ...:
           `def sage...=factorize(sage...);`

    AUTHORS:

    - Alexander Raichev (2012-07-26)
    - Daniel Krenn (2014-12-01)
    """

    def __init__(self, parent, numerator, denominator_factored, reduce=True):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: df = [x, 1], [y, 1], [x*y+1, 1]
            sage: f = FFPD(x, df)
            sage: TestSuite(f).run()
        """
        super().__init__(parent)

        from sage.rings.semirings.non_negative_integer_semiring import NN
        self._numerator = parent._numerator_ring(numerator)
        self._denominator_factored = list(
            (parent._denominator_ring(d), NN(n))
            for d, n in denominator_factored)

        R = self.denominator_ring
        if numerator in R and reduce:
            # Reduce fraction if possible.
            numer = R(self._numerator)
            df = self._denominator_factored
            new_df = []
            for (q, e) in df:
                ee = e
                quo, rem = numer.quo_rem(q)
                while rem == 0 and ee > 0:
                    ee -= 1
                    numer = quo
                    quo, rem = numer.quo_rem(q)
                if ee > 0:
                    new_df.append((q, ee))
            self._numerator = numer
            self._denominator_factored = new_df

    def numerator(self):
        r"""
        Return the numerator of ``self``.

        OUTPUT: the numerator

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: H = (1 - x - y - x*y)**2*(1-x)
            sage: Hfac = H.factor()
            sage: G = exp(y)/Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: F.numerator()
            -e^y
        """
        return self._numerator

    def denominator(self):
        r"""
        Return the denominator of ``self``.

        OUTPUT:

        The denominator (i.e., the product of the factored denominator).

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: H = (1 - x - y - x*y)**2*(1-x)
            sage: Hfac = H.factor()
            sage: G = exp(y)/Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: F.denominator()
            x^3*y^2 + 2*x^3*y + x^2*y^2 + x^3 - 2*x^2*y - x*y^2 - 3*x^2 - 2*x*y
            - y^2 + 3*x + 2*y - 1
        """
        return prod(q ** e for q, e in self.denominator_factored())

    def denominator_factored(self):
        r"""
        Return the factorization in ``self.denominator_ring`` of the denominator of
        ``self`` but without the unit part.

        OUTPUT:

        The factored denominator as a list of tuple ``(f, m)``, where `f` is
        a factor and `m` its multiplicity.

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: H = (1 - x - y - x*y)**2*(1-x)
            sage: Hfac = H.factor()
            sage: G = exp(y)/Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: F.denominator_factored()
            [(x - 1, 1), (x*y + x + y - 1, 2)]
        """
        return self._denominator_factored

    @property
    def denominator_ring(self):
        r"""
        Return the ring of the denominator.

        OUTPUT: a ring

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: H = (1 - x - y - x*y)**2*(1-x)
            sage: Hfac = H.factor()
            sage: G = exp(y)/Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: F.denominator_ring
            Multivariate Polynomial Ring in x, y over Rational Field
            sage: F = FFPD(G/H)
            sage: F
            (e^y, [(x - 1, 1), (x*y + x + y - 1, 2)])
            sage: F.denominator_ring
            Multivariate Polynomial Ring in x, y over Rational Field
        """
        return self.parent()._denominator_ring

    @property
    def numerator_ring(self):
        r"""
        Return the ring of the numerator.

        OUTPUT: a ring

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: H = (1 - x - y - x*y)**2*(1-x)
            sage: Hfac = H.factor()
            sage: G = exp(y)/Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: F.numerator_ring
            Symbolic Ring
            sage: F = FFPD(G/H)
            sage: F
            (e^y, [(x - 1, 1), (x*y + x + y - 1, 2)])
            sage: F.numerator_ring
            Symbolic Ring
        """
        return self.parent()._numerator_ring

    def dimension(self):
        r"""
        Return the number of indeterminates of ``self.denominator_ring``.

        OUTPUT: integer

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: H = (1 - x - y - x*y)**2*(1-x)
            sage: Hfac = H.factor()
            sage: G = exp(y)/Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: F.dimension()
            2
        """
        from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
        from sage.rings.polynomial.multi_polynomial_ring_base import MPolynomialRing_base
        R = self.denominator_ring
        if isinstance(R, PolynomialRing_general) or isinstance(R, MPolynomialRing_base):
            return R.ngens()
        raise NotImplementedError('only polynomial rings are supported as base')

    def quotient(self):
        r"""
        Convert ``self`` into a quotient.

        OUTPUT: an element

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: H = (1 - x - y - x*y)**2*(1-x)
            sage: Hfac = H.factor()
            sage: G = exp(y)/Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: F
            (-e^y, [(x - 1, 1), (x*y + x + y - 1, 2)])
            sage: F.quotient()
            -e^y/(x^3*y^2 + 2*x^3*y + x^2*y^2 + x^3 - 2*x^2*y - x*y^2 - 3*x^2 -
            2*x*y - y^2 + 3*x + 2*y - 1)
        """
        return self.numerator() / self.denominator()

    def _repr_(self) -> str:
        r"""
        Return a string representation of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: f = FFPD(x + y, [(y, 1), (x, 1)])
            sage: f
            (x + y, [(y, 1), (x, 1)])
        """
        return repr((self.numerator(), self.denominator_factored()))

    _richcmp_ = richcmp_by_eq_and_lt("_eq_", "_lt_")

    def _eq_(self, other):
        r"""
        Return whether the FFPD instance ``other`` is equal to
        this FFPD instance.

        Two FFPD instances are equal iff they represent the same
        fraction.

        INPUT:

        - ``other`` -- an instance of :class:`FractionWithFactoredDenominator`

        OUTPUT: boolean

        It can be assumed that ``self`` and ``other`` have the same parent.

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: df = [x, 1], [y, 1], [x*y+1, 1]
            sage: f = FFPD(x, df)
            sage: ff = FFPD(x, df, reduce=False)
            sage: f == ff
            True
            sage: g = FFPD(y, df)
            sage: g == f
            False

        ::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: G = exp(x + y)
            sage: H = (1 - 2*x - y) * (1 - x - 2*y)
            sage: a = FFPD(G/H)
            sage: b = FFPD(G, H.factor())
            sage: bool(a == b)
            True

        ::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: df = [x, 1], [y, 1], [x*y+1, 1]
            sage: f = FFPD(x, df)
            sage: ff = FFPD(x, df, reduce=False)
            sage: f != ff
            False
            sage: g = FFPD(y, df)
            sage: g != f
            True

        TESTS::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: f = FFPD(x, [])
            sage: f == x
            True

        ::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: f = FFPD(x*y, [(x-1, 1), (y-2, 2)])
            sage: g = FFPD(x, [(x-1, 1), (y-2, 2)])
            sage: f == f
            True
            sage: f == g
            False

            sage: f < g
            Traceback (most recent call last):
            ...
            AttributeError:
            'FractionWithFactoredDenominatorRing_with_category.element_class'
            object has no attribute '_lt_'
        """
        return (self.numerator() * other.denominator() ==
                other.numerator() * self.denominator())

    def _total_order_key_(self):
        r"""
        Return a key that can be used for sorting.

        FFPD ``A`` is less than FFPD ``B`` iff
        (the denominator factorization of ``A`` is shorter than that of ``B``)
        of (the denominator factorization lengths are equal and
        the denominator of ``A`` is less than that of ``B`` in their ring) or
        (the denominator factorization lengths are equal and the
        denominators are equal and the numerator of ``A`` is less than that
        of ``B`` in their ring).

        OUTPUT: a tuple

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: df = [x, 1], [y, 1], [x*y+1, 1]
            sage: f = FFPD(x, df); f
            (1, [(y, 1), (x*y + 1, 1)])
            sage: ff = FFPD(x, df, reduce=False); ff
            (x, [(y, 1), (x, 1), (x*y + 1, 1)])
            sage: g = FFPD(y, df)
            sage: h = FFPD(exp(x), df)
            sage: i = FFPD(sin(x + 2), df)
            sage: f._total_order_key_() < ff._total_order_key_()
            True
            sage: f._total_order_key_() < g._total_order_key_()
            True
            sage: g._total_order_key_() < h._total_order_key_()
            True
            sage: bool(h._total_order_key_() < i._total_order_key_())
            False
        """
        return (len(self.denominator_factored()),
                self.denominator(),
                self.numerator())

    def univariate_decomposition(self):
        r"""
        Return the usual univariate partial fraction decomposition
        of ``self``.

        Assume that the numerator of ``self`` lies in the same univariate
        factorial polynomial ring as the factors of the denominator.

        Let `f = p/q` be a rational expression where `p` and `q` lie in a
        univariate factorial polynomial ring `R`.
        Let `q_1^{e_1} \cdots q_n^{e_n}` be the
        unique factorization of `q` in `R` into irreducible factors.
        Then `f` can be written uniquely as:

        .. MATH::

            (*) \quad p_0 + \sum_{i=1}^{m} \frac{p_i}{q_i^{e_i}},

        for some `p_j \in R`.
        We call `(*)` the *usual partial fraction decomposition* of `f`.

        .. NOTE::

            This partial fraction decomposition can be computed using
            :meth:`~sage.symbolic.expression.Expression.partial_fraction` or
            :meth:`~sage.categories.quotient_fields.QuotientFields.ElementMethods.partial_fraction_decomposition`
            as well. However, here we use the already obtained/cached
            factorization of the denominator. This gives a speed up for
            non-small instances.

        OUTPUT: an instance of :class:`FractionWithFactoredDenominatorSum`

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing

        One variable::

            sage: R.<x> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: f = 5*x^3 + 1/x + 1/(x-1) + 1/(3*x^2 + 1)
            sage: f
            (5*x^7 - 5*x^6 + 5/3*x^5 - 5/3*x^4 + 2*x^3 - 2/3*x^2 + 1/3*x - 1/3)/(x^4 - x^3 + 1/3*x^2 - 1/3*x)
            sage: decomp = FFPD(f).univariate_decomposition()
            sage: decomp
            (5*x^3, []) +
            (1, [(x - 1, 1)]) +
            (1, [(x, 1)]) +
            (1/3, [(x^2 + 1/3, 1)])
            sage: decomp.sum().quotient() == f
            True

        One variable with numerator in symbolic ring::

            sage: R.<x> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: f = 5*x^3 + 1/x + 1/(x-1) + exp(x)/(3*x^2 + 1)
            sage: f
            (5*x^5 - 5*x^4 + 2*x - 1)/(x^2 - x) + e^x/(3*x^2 + 1)
            sage: decomp = FFPD(f).univariate_decomposition()
            sage: decomp
            (0, []) +
            (15/4*x^7 - 15/4*x^6 + 5/4*x^5 - 5/4*x^4 + 3/2*x^3 + 1/4*x^2*e^x -
             3/4*x^2 - 1/4*x*e^x + 1/2*x - 1/4, [(x - 1, 1)]) +
            (-15*x^7 + 15*x^6 - 5*x^5 + 5*x^4 - 6*x^3 -
             x^2*e^x + 3*x^2 + x*e^x - 2*x + 1, [(x, 1)]) +
            (1/4*(15*x^7 - 15*x^6 + 5*x^5 - 5*x^4 + 6*x^3 + x^2*e^x -
                  3*x^2 - x*e^x + 2*x - 1)*(3*x - 1), [(x^2 + 1/3, 1)])

        One variable over a finite field::

            sage: R.<x> = PolynomialRing(GF(2))
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: f = 5*x^3 + 1/x + 1/(x-1) + 1/(3*x^2 + 1)
            sage: f
            (x^6 + x^4 + 1)/(x^3 + x)
            sage: decomp = FFPD(f).univariate_decomposition()
            sage: decomp
            (x^3, []) + (1, [(x, 1)]) + (x, [(x + 1, 2)])
            sage: decomp.sum().quotient() == f
            True

        One variable over an inexact field::

            sage: R.<x> = PolynomialRing(CC)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: f = 5*x^3 + 1/x + 1/(x-1) + 1/(3*x^2 + 1)
            sage: f
            (5.00000000000000*x^7 - 5.00000000000000*x^6 + 1.66666666666667*x^5 - 1.66666666666667*x^4 + 2.00000000000000*x^3 - 0.666666666666667*x^2 + 0.333333333333333*x - 0.333333333333333)/(x^4 - x^3 + 0.333333333333333*x^2 - 0.333333333333333*x)
            sage: decomp = FFPD(f).univariate_decomposition()
            sage: decomp
            (5.00000000000000*x^3, []) +
            (1.00000000000000, [(x - 1.00000000000000, 1)]) +
            (-0.288675134594813*I, [(x - 0.577350269189626*I, 1)]) +
            (1.00000000000000, [(x, 1)]) +
            (0.288675134594813*I, [(x + 0.577350269189626*I, 1)])
            sage: decomp.sum().quotient() == f # Rounding error coming
            False

        TESTS::

            sage: R.<x> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: f = exp(x) / (x^2-x)
            sage: f
            e^x/(x^2 - x)
            sage: FFPD(f).univariate_decomposition()
            (0, []) + (e^x, [(x - 1, 1)]) + (-e^x, [(x, 1)])

        AUTHORS:

        - Robert Bradshaw (2007-05-31)
        - Alexander Raichev (2012-06-25)
        - Daniel Krenn (2014-12-01)
        """
        if self.dimension() > 1:
            return FractionWithFactoredDenominatorSum([self])

        R = self.denominator_ring
        p = self.numerator()
        q = self.denominator()
        try:
            whole, p = R(p).quo_rem(q)
            mn = R.one()
        except (TypeError, ValueError):
            whole = R(0)
            mn = p
            p = R.one()
        df = self.denominator_factored()
        decomp = [self.parent()(whole, [])]
        denominator = prod(b**n for b, n in df)
        for a, m in df:
            am = a**m
            q, r = denominator.quo_rem(am)
            assert r == 0
            numer = p * q.inverse_mod(am) % am
            # The inverse exists because the product and a**m
            # are relatively prime.
            decomp.append(self.parent()(mn * numer, [(a, m)]))
        return FractionWithFactoredDenominatorSum(decomp)

    def nullstellensatz_certificate(self):
        r"""
        Return a Nullstellensatz certificate of ``self`` if it exists.

        Let `[(q_1, e_1), \ldots, (q_n, e_n)]` be the denominator
        factorization of ``self``. The Nullstellensatz certificate is
        a list of polynomials `h_1, \ldots, h_m` in ``self.denominator_ring``
        that satisfies `h_1 q_1 + \cdots + h_m q_n = 1` if it exists.

        .. NOTE::

            Only works for multivariate base rings.

        OUTPUT:

        A list of polynomials or ``None`` if no Nullstellensatz
        certificate exists.

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: G = sin(x)
            sage: H = x^2 * (x*y + 1)
            sage: f = FFPD(G, H.factor())
            sage: L = f.nullstellensatz_certificate()
            sage: L
            [y^2, -x*y + 1]
            sage: df = f.denominator_factored()
            sage: sum(L[i]*df[i][0]**df[i][1] for i in range(len(df))) == 1
            True

        ::

            sage: f = 1/(x*y)
            sage: L = FFPD(f).nullstellensatz_certificate()
            sage: L is None
            True
        """
        R = self.denominator_ring
        df = self.denominator_factored()
        J = R.ideal([q ** e for q, e in df])
        if R.one() in J:
            return R.one().lift(J)
        return None

    def nullstellensatz_decomposition(self):
        r"""
        Return a Nullstellensatz decomposition of ``self``.

        Let `f = p/q` where `q` lies in a `d` -variate polynomial ring
        `K[X]` for some field `K` and `d \geq 1`.
        Let `q_1^{e_1} \cdots q_n^{e_n}` be the
        unique factorization of `q` in `K[X]` into irreducible factors and
        let `V_i` be the algebraic variety `\{x \in L^d \mid q_i(x) = 0\}`
        of `q_i` over the algebraic closure `L` of `K`.
        By [Rai2012]_, `f` can be written as

        .. MATH::

            (*) \quad \sum_A \frac{p_A}{\prod_{i \in A} q_i^{e_i}},

        where the `p_A` are products of `p` and elements in `K[X]` and
        the sum is taken over all subsets
        `A \subseteq \{1, \ldots, m\}` such that
        `\bigcap_{i\in A} T_i \neq \emptyset`.

        We call `(*)` a *Nullstellensatz decomposition* of `f`.
        Nullstellensatz decompositions are not unique.

        The algorithm used comes from [Rai2012]_.

        .. NOTE::

            Recursive. Only works for multivariate ``self``.

        OUTPUT: an instance of :class:`FractionWithFactoredDenominatorSum`

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import *
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: f = 1/(x*(x*y + 1))
            sage: decomp = FFPD(f).nullstellensatz_decomposition()
            sage: decomp
            (0, []) + (1, [(x, 1)]) + (-y, [(x*y + 1, 1)])
            sage: decomp.sum().quotient() == f
            True
            sage: [r.nullstellensatz_certificate() is None for r in decomp]
            [True, True, True]

        ::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: G = sin(y)
            sage: H = x*(x*y + 1)
            sage: f = FFPD(G, H.factor())
            sage: decomp = f.nullstellensatz_decomposition()
            sage: decomp
            (0, []) + (sin(y), [(x, 1)]) + (-y*sin(y), [(x*y + 1, 1)])
            sage: bool(decomp.sum().quotient() == G/H)
            True
            sage: [r.nullstellensatz_certificate() is None for r in decomp]
            [True, True, True]
        """
        L = self.nullstellensatz_certificate()
        if L is None:
            # No decomposing possible.
            return FractionWithFactoredDenominatorSum([self])

        # Otherwise decompose recursively.
        decomp = FractionWithFactoredDenominatorSum()
        p = self.numerator()
        df = self.denominator_factored()
        m = len(df)
        iteration1 = FractionWithFactoredDenominatorSum(
            [self.parent()(p * L[i], [df[j] for j in range(m) if j != i])
             for i in range(m) if L[i] != 0])

        # Now decompose each FFPD of iteration1.
        for r in iteration1:
            decomp.extend(r.nullstellensatz_decomposition())

        # Simplify and return result.
        return decomp._combine_like_terms_().whole_and_parts()

    def algebraic_dependence_certificate(self):
        r"""
        Return the algebraic dependence certificate of ``self``.

        The algebraic dependence certificate is the ideal `J` of
        annihilating polynomials for the set of polynomials
        ``[q^e for (q, e) in self.denominator_factored()]``,
        which could be the zero ideal.
        The ideal `J` lies in a polynomial ring over the field
        ``self.denominator_ring.base_ring()`` that has
        ``m = len(self.denominator_factored())`` indeterminates.

        OUTPUT: an ideal

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: f = 1/(x^2 * (x*y + 1) * y^3)
            sage: ff = FFPD(f)
            sage: J = ff.algebraic_dependence_certificate(); J
            Ideal (1 - 6*T2 + 15*T2^2 - 20*T2^3 + 15*T2^4 - T0^2*T1^3 -
             6*T2^5  + T2^6) of Multivariate Polynomial Ring in
             T0, T1, T2 over Rational Field
            sage: g = J.gens()[0]
            sage: df = ff.denominator_factored()
            sage: g(*(q**e for q, e in df)) == 0
            True

        ::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: G = exp(x + y)
            sage: H = x^2 * (x*y + 1) * y^3
            sage: ff = FFPD(G, H.factor())
            sage: J = ff.algebraic_dependence_certificate(); J
            Ideal (1 - 6*T2 + 15*T2^2 - 20*T2^3 + 15*T2^4 - T0^2*T1^3 -
            6*T2^5 + T2^6) of Multivariate Polynomial Ring in
            T0, T1, T2 over Rational Field
            sage: g = J.gens()[0]
            sage: df = ff.denominator_factored()
            sage: g(*(q**e for q, e in df)) == 0
            True

        ::

            sage: f = 1/(x^3 * y^2)
            sage: J = FFPD(f).algebraic_dependence_certificate()
            sage: J
            Ideal (0) of Multivariate Polynomial Ring in T0, T1 over Rational Field

        ::

            sage: f = sin(1)/(x^3 * y^2)
            sage: J = FFPD(f).algebraic_dependence_certificate()
            sage: J
            Ideal (0) of Multivariate Polynomial Ring in T0, T1 over Rational Field
        """
        R = self.denominator_ring
        df = self.denominator_factored()
        if not df:
            return R.ideal()    # The zero ideal.
        m = len(df)
        F = R.base_ring()
        Xs = list(R.gens())
        d = len(Xs)

        # Expand R by 2 * m new variables.
        S = 'S'
        while S in [str(x) for x in Xs]:
            S = S + 'S'
        Ss = [S + str(i) for i in range(m)]
        T = 'T'
        while T in [str(x) for x in Xs]:
            T = T + 'T'
        Ts = [T + str(i) for i in range(m)]

        Vs = [str(x) for x in Xs] + Ss + Ts
        RR = PolynomialRing(F, Vs)
        Xs = RR.gens()[:d]
        Ss = RR.gens()[d: d + m]
        Ts = RR.gens()[d + m: d + 2 * m]

        # Compute the appropriate elimination ideal.
        J = RR.ideal([Ss[j] - RR(df[j][0]) for j in range(m)] +
                     [Ss[j] ** df[j][1] - Ts[j] for j in range(m)])
        J = J.elimination_ideal(Xs + Ss)

        # Coerce J into the polynomial ring in the indeterminates Ts[m:].
        # I choose the negdeglex order because i find it useful in my work.
        RRR = PolynomialRing(F, [str(t) for t in Ts], order='negdeglex')
        return RRR.ideal(J)

    def algebraic_dependence_decomposition(self, whole_and_parts=True):
        r"""
        Return an algebraic dependence decomposition of ``self``.

        Let `f = p/q` where `q` lies in a `d`-variate polynomial ring
        `K[X]` for some field `K`.
        Let `q_1^{e_1} \cdots q_n^{e_n}` be the
        unique factorization of `q` in `K[X]` into irreducible factors and
        let `V_i` be the algebraic variety `\{x \in L^d \mid q_i(x) = 0\}`
        of `q_i` over the algebraic closure `L` of `K`.
        By [Rai2012]_, `f` can be written as

        .. MATH::

            (*) \quad \sum_A \frac{p_A}{\prod_{i \in A} q_i^{b_i}},

        where the `b_i` are positive integers, each `p_A` is a products
        of `p` and an element in `K[X]`,
        and the sum is taken over all subsets
        `A \subseteq \{1, \ldots, m\}` such that `|A| \leq d` and
        `\{q_i \mid i \in A\}` is algebraically independent.

        We call `(*)` an *algebraic dependence decomposition* of `f`.
        Algebraic dependence decompositions are not unique.

        The algorithm used comes from [Rai2012]_.

        OUTPUT: an instance of :class:`FractionWithFactoredDenominatorSum`

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: f = 1/(x^2 * (x*y + 1) * y^3)
            sage: ff = FFPD(f)
            sage: decomp = ff.algebraic_dependence_decomposition()
            sage: decomp
            (0, []) + (-x, [(x*y + 1, 1)]) +
            (x^2*y^2 - x*y + 1, [(y, 3), (x, 2)])
            sage: decomp.sum().quotient() == f
            True
            sage: for r in decomp:
            ....:     J = r.algebraic_dependence_certificate()
            ....:     J is None or J == J.ring().ideal()  # The zero ideal
            True
            True
            True

        ::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: G = sin(x)
            sage: H = x^2 * (x*y + 1) * y^3
            sage: f = FFPD(G, H.factor())
            sage: decomp = f.algebraic_dependence_decomposition()
            sage: decomp
            (0, []) + (x^4*y^3*sin(x), [(x*y + 1, 1)]) +
            (-(x^5*y^5 - x^4*y^4 + x^3*y^3 - x^2*y^2 + x*y - 1)*sin(x),
             [(y, 3), (x, 2)])
            sage: bool(decomp.sum().quotient() == G/H)
            True
            sage: for r in decomp:
            ....:     J = r.algebraic_dependence_certificate()
            ....:     J is None or J == J.ring().ideal()
            True
            True
            True
        """
        J = self.algebraic_dependence_certificate()
        if not J:
            # No decomposing possible.
            return FractionWithFactoredDenominatorSum([self])

        # Otherwise decompose recursively.
        decomp = FractionWithFactoredDenominatorSum()
        p = self.numerator()
        df = self.denominator_factored()
        m = len(df)
        g = J.gens()[0]     # An annihilating polynomial for df.
        new_vars = J.ring().gens()
        # Note that each new_vars[j] corresponds to df[j] such that
        # g([q**e for q, e in df]) = 0.
        # Assuming here that g.parent() has negdeglex term order
        # so that g.lt() is indeed the monomial we want below.
        # Use g to rewrite r into a sum of FFPDs,
        # each with < m distinct denominator factors.
        gg = (g.lt() - g) / (g.lc())
        numers = map(prod, zip(gg.coefficients(), gg.monomials()))
        e = list(g.lt().exponents())[0: m]
        denoms = [(new_vars[j], e[0][j] + 1) for j in range(m)]
        # Write r in terms of new_vars,
        # cancel factors in the denominator, and combine like terms.
        FFPD = FractionWithFactoredDenominatorRing(J.ring())
        iteration1_temp = FractionWithFactoredDenominatorSum(
            [FFPD(a, denoms) for a in numers])._combine_like_terms_()
        # Substitute in df.
        qpowsub = {new_vars[j]: df[j][0] ** df[j][1] for j in range(m)}
        iteration1 = FractionWithFactoredDenominatorSum()
        for r in iteration1_temp:
            num1 = p * J.ring()(r.numerator()).subs(qpowsub)
            denoms1 = []
            for q, e in r.denominator_factored():
                j = new_vars.index(q)
                denoms1.append((df[j][0], df[j][1] * e))
            iteration1.append(self.parent()(num1, denoms1))
        # Now decompose each FFPD of iteration1.
        for r in iteration1:
            decomp.extend(r.algebraic_dependence_decomposition())

        # Simplify and return result.
        return decomp._combine_like_terms_().whole_and_parts()

    def leinartas_decomposition(self):
        r"""
        Return a Leinartas decomposition of ``self``.

        Let `f = p/q` where `q` lies in a `d` -variate polynomial
        ring `K[X]` for some field `K`.
        Let `q_1^{e_1} \cdots q_n^{e_n}` be the
        unique factorization of `q` in `K[X]` into irreducible factors and
        let `V_i` be the algebraic variety
        `\{x\in L^d \mid q_i(x) = 0\}` of `q_i` over the algebraic closure
        `L` of `K`. By [Rai2012]_, `f` can be written as

        .. MATH::

            (*) \quad \sum_A \frac{p_A}{\prod_{i \in A} q_i^{b_i}},

        where the `b_i` are positive integers, each `p_A` is a product of
        `p` and an element of `K[X]`, and the sum is taken over all
        subsets `A \subseteq \{1, \ldots, m\}` such that

        1. `|A| \le d`,
        2. `\bigcap_{i\in A} T_i \neq \emptyset`, and
        3. `\{q_i \mid i\in A\}` is algebraically independent.

        In particular, any rational expression in `d` variables
        can be represented as a sum of rational expressions
        whose denominators each contain at most `d` distinct irreducible
        factors.

        We call `(*)` a *Leinartas decomposition* of `f`.
        Leinartas decompositions are not unique.

        The algorithm used comes from [Rai2012]_.

        OUTPUT: an instance of :class:`FractionWithFactoredDenominatorSum`

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: f = (x^2 + 1)/((x + 2)*(x - 1)*(x^2 + x + 1))
            sage: decomp = FFPD(f).leinartas_decomposition()
            sage: decomp
            (0, []) + (2/9, [(x - 1, 1)]) +
            (-5/9, [(x + 2, 1)]) + (1/3*x, [(x^2 + x + 1, 1)])
            sage: decomp.sum().quotient() == f
            True

        ::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: f = 1/x + 1/y + 1/(x*y + 1)
            sage: decomp = FFPD(f).leinartas_decomposition()
            sage: decomp
            (0, []) + (1, [(x*y + 1, 1)]) + (x + y, [(y, 1), (x, 1)])
            sage: decomp.sum().quotient() == f
            True
            sage: def check_decomp(r):
            ....:     L = r.nullstellensatz_certificate()
            ....:     J = r.algebraic_dependence_certificate()
            ....:     return L is None and (J is None or J == J.ring().ideal())
            sage: all(check_decomp(r) for r in decomp)
            True

        ::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: f = sin(x)/x + 1/y + 1/(x*y + 1)
            sage: G = f.numerator()
            sage: H = R(f.denominator())
            sage: ff = FFPD(G, H.factor())
            sage: decomp = ff.leinartas_decomposition()
            sage: decomp  # random - non canonical depends on singular version
            (0, []) +
            (-(x*y^2*sin(x) + x^2*y + x*y + y*sin(x) + x)*y, [(y, 1)]) +
            ((x*y^2*sin(x) + x^2*y + x*y + y*sin(x) + x)*x*y, [(x*y + 1, 1)]) +
            (x*y^2*sin(x) + x^2*y + x*y + y*sin(x) + x, [(y, 1), (x, 1)])
            sage: bool(decomp.sum().quotient() == f)
            True
            sage: all(check_decomp(r) for r in decomp)
            True

        ::

            sage: R.<x,y,z>= PolynomialRing(GF(2, 'a'))
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: f = 1/(x * y * z * (x*y + z))
            sage: decomp = FFPD(f).leinartas_decomposition()
            sage: decomp
            (0, []) + (1, [(z, 2), (x*y + z, 1)]) +
            (1, [(z, 2), (y, 1), (x, 1)])
            sage: decomp.sum().quotient() == f
            True
        """
        if self.dimension() == 1:
            # Sage's lift() function doesn't work in
            # univariate polynomial rings.
            # So nullstellensatz_decomposition() won't work.
            # Can use algebraic_dependence_decomposition(),
            # which is sufficient.
            # temp = FractionWithFactoredDenominatorSum([self])
            # Alternatively can use univariate_decomposition(),
            # which is more efficient.
            return self.univariate_decomposition()
        temp = self.nullstellensatz_decomposition()
        decomp = FractionWithFactoredDenominatorSum()
        for r in temp:
            decomp.extend(r.algebraic_dependence_decomposition())

        # Simplify and return result.
        return decomp._combine_like_terms_().whole_and_parts()

    def cohomology_decomposition(self):
        r"""
        Return the cohomology decomposition of ``self``.

        Let `p / (q_1^{e_1} \cdots q_n^{e_n})` be the fraction represented
        by ``self`` and let `K[x_1, \ldots, x_d]` be the polynomial ring
        in which the `q_i` lie.
        Assume that `n \leq d` and that the gradients of the `q_i` are linearly
        independent at all points in the intersection
        `V_1 \cap \ldots \cap V_n` of the algebraic varieties
        `V_i = \{x \in L^d \mid q_i(x) = 0 \}`, where `L` is the algebraic
        closure of the field `K`.
        Return a :class:`FractionWithFactoredDenominatorSum`
        `f` such that the differential form
        `f dx_1 \wedge \cdots \wedge dx_d` is de Rham cohomologous to the
        differential form
        `p / (q_1^{e_1} \cdots q_n^{e_n}) dx_1 \wedge \cdots \wedge dx_d`
        and such that the denominator of each summand of `f` contains
        no repeated irreducible factors.

        The algorithm used here comes from the proof of Theorem 17.4 of
        [AY1983]_.

        OUTPUT: an instance of :class:`FractionWithFactoredDenominatorSum`

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: f = 1/(x^2 + x + 1)^3
            sage: decomp = FFPD(f).cohomology_decomposition()
            sage: decomp
            (0, []) + (2/3, [(x^2 + x + 1, 1)])

        ::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: FFPD(1, [(x, 1), (y, 2)]).cohomology_decomposition()
            (0, [])

        The following example was fixed in :issue:`29465`::

            sage: p = 1
            sage: qs = [(x*y - 1, 1), (x**2 + y**2 - 1, 2)]
            sage: f = FFPD(p, qs)
            sage: f.cohomology_decomposition()
            (0, []) + (-4/3*x*y, [(x^2 + y^2 - 1, 1)]) +
            (1/3, [(x*y - 1, 1), (x^2 + y^2 - 1, 1)])
        """
        from sage.calculus.functions import jacobian
        from sage.arith.misc import XGCD as xgcd
        from sage.sets.set import Set

        R = self.denominator_ring
        df = self.denominator_factored()
        n = len(df)
        if sum(e for (q, e) in df) <= n:
            # No decomposing possible.
            return FractionWithFactoredDenominatorSum([self])

        # Otherwise decompose recursively.
        decomp = FractionWithFactoredDenominatorSum()
        p = self.numerator()
        qs = [q for (q, e) in df]
        # sort according to the term order of R
        X = sorted(R.gens())
        var_sets_n = sorted(sorted(s) for s in Set(X).subsets(n))
        Par = self.parent()

        # Compute Jacobian determinants for qs.
        dets = [R(jacobian(qs, x).determinant())
                for x in var_sets_n]

        # Get a Nullstellensatz certificate for qs and dets.
        if self.dimension() == 1:
            # Sage's lift() function doesn't work in
            # univariate polynomial rings.
            # So use xgcd(), which does the same thing in this case.
            # Note that by assumption qs and dets have length 1.
            L = xgcd(qs[0], dets[0])[1:]
        else:
            L = R.one().lift(R.ideal(qs + dets))

        # Do first iteration of decomposition.
        iteration1 = FractionWithFactoredDenominatorSum()
        # Contributions from qs.
        for i in range(n):
            if L[i] == 0:
                continue
            # Cancel one df[i] from denominator.
            new_df = [list(t) for t in df]
            if new_df[i][1] > 1:
                new_df[i][1] -= 1
            else:
                del new_df[i]
            iteration1.append(Par(p * L[i], new_df))

        # Contributions from dets.
        # Compute each contribution's cohomologous form using
        # the least index j such that new_df[j][1] > 1.
        # Know such an index exists by first 'if' statement at
        # the top.
        for j in range(n):
            if df[j][1] > 1:
                J = j
                break
        new_df = [list(t) for t in df]
        new_df[J][1] -= 1
        for k, x in enumerate(var_sets_n):
            if L[n + k] == 0:
                continue
            # Compute Jacobian in the Symbolic Ring.
            jac = jacobian([SR(p * L[n + k])] +
                           [SR(qs[j]) for j in range(n) if j != J],
                           [SR(xx) for xx in x])
            det = jac.determinant()
            # The parity epsilon from [AY1983, eq. (17.11)] does not
            # enter this computation, since we do not order the
            # coordinates x to the front of X in the representation of
            # this differential form (:issue:`29465`).
            iteration1.append(Par((-1) ** J * det / new_df[J][1], new_df))

        # Now decompose each FFPD of iteration1.
        for r in iteration1:
            decomp.extend(r.cohomology_decomposition())

        # Simplify and return result.
        return decomp._combine_like_terms_().whole_and_parts()

    def asymptotic_decomposition(self, alpha, asy_var=None):
        r"""
        Return the asymptotic decomposition of ``self``.

        The asymptotic decomposition of `F` is a sum that has the
        same asymptotic expansion as `f` in the direction ``alpha``
        but each summand has a denominator factorization of the form
        `[(q_1, 1), \ldots, (q_n, 1)]`, where `n` is at most the
        :meth:`dimension` of `F`.

        INPUT:

        - ``alpha`` -- a `d`-tuple of positive integers or symbolic variables
        - ``asy_var`` -- (default: ``None``) a symbolic variable with
          respect to which to compute asymptotics;
          if ``None`` is given, we set ``asy_var = var('r')``

        OUTPUT: an instance of :class:`FractionWithFactoredDenominatorSum`

        The output results from a Leinartas decomposition followed by a
        cohomology decomposition.

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: f = (x^2 + 1)/((x - 1)^3*(x + 2))
            sage: F = FFPD(f)
            sage: alpha = [var('a')]
            sage: F.asymptotic_decomposition(alpha)
            (0, []) +
            (1/54*(5*a^2 + 2*a^2/x + 11*a^2/x^2)*r^2
             - 1/54*(5*a - 2*a/x - 33*a/x^2)*r + 11/27/x^2,
            [(x - 1, 1)]) + (-5/27, [(x + 2, 1)])

        ::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: H = (1 - 2*x -y)*(1 - x -2*y)**2
            sage: Hfac = H.factor()
            sage: G = 1/Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: alpha = var('a, b')
            sage: F.asymptotic_decomposition(alpha)
            (0, []) +
            (-1/3*r*(a/x - 2*b/y) - 1/3/x + 2/3/y,
             [(x + 2*y - 1, 1), (2*x + y - 1, 1)])
        """
        R = self.denominator_ring
        d = self.dimension()
        n = len(self.denominator_factored())
        X = [SR(x) for x in R.gens()]

        # Reduce number of distinct factors in denominator of self
        # down to at most d.
        decomp1 = FractionWithFactoredDenominatorSum([self])
        if n > d:
            decomp1 = decomp1[0].leinartas_decomposition()

        # Reduce to no repeated factors in denominator of each element
        # of decomp1.
        # Compute the cohomology decomposition for each
        # Cauchy differential form generated by each element of decomp.
        if asy_var is None:
            asy_var = var('r')
        cauchy_stuff = prod([X[j] ** (-alpha[j] * asy_var - 1)
                             for j in range(d)])
        decomp2 = FractionWithFactoredDenominatorSum()
        for f in decomp1:
            ff = self.parent()(f.numerator() * cauchy_stuff,
                               f.denominator_factored())
            decomp2.extend(ff.cohomology_decomposition())
        decomp2 = decomp2._combine_like_terms_()

        # Divide out cauchy_stuff from integrands.
        decomp3 = FractionWithFactoredDenominatorSum()
        for f in decomp2:
            ff = self.parent()((f.numerator() /
                                cauchy_stuff).simplify_full().collect(asy_var),
                               f.denominator_factored())
            decomp3.append(ff)

        return decomp3

    def asymptotics(self, p, alpha, N, asy_var=None, numerical=0,
                    verbose=False):
        r"""
        Return the asymptotics in the given direction.

        This function returns the first `N` terms (some of which could be
        zero) of the asymptotic expansion of the Maclaurin ray coefficients
        `F_{r \alpha}` of the function `F` represented by ``self`` as
        `r \to \infty`, where `r` is ``asy_var`` and ``alpha`` is a tuple
        of positive integers of length `d` which is ``self.dimension()``.
        Assume that

        - `F` is holomorphic in a neighborhood of the origin;
        - the unique factorization of the denominator `H` of `F` in the local
          algebraic ring at `p` equals its unique factorization in the local
          analytic ring at `p`;
        - the unique factorization of `H` in the local algebraic ring at `p`
          has at most ``d`` irreducible factors, none of which are repeated
          (one can reduce to this case via :meth:`asymptotic_decomposition()`);
        - `p` is a convenient strictly minimal smooth or multiple point
          with all nonzero coordinates that is critical and nondegenerate
          for ``alpha``.

        The algorithms used here come from [RW2008]_ and [RW2012]_.

        INPUT:

        - ``p`` -- dictionary with keys that can be coerced to equal
          ``self.denominator_ring.gens()``
        - ``alpha`` -- tuple of length ``self.dimension()`` of positive
          integers or, if `p` is a smooth point, possibly of symbolic variables
        - ``N`` -- positive integer
        - ``asy_var`` -- (default: ``None``) a symbolic variable for the
          asymptotic expansion; if ``none`` is given, then ``var('r')`` will be
          assigned
        - ``numerical`` -- (default: 0) a natural number; if ``numerical`` is
          greater than 0, then return a numerical approximation of
          `F_{r \alpha}` with ``numerical`` digits of precision; otherwise
          return exact values
        - ``verbose`` -- boolean (default: ``False``); print the current state of
          the algorithm

        OUTPUT:

        The tuple ``(asy, exp_scale, subexp_part)``.
        Here ``asy`` is the sum of the first `N` terms (some of which might
        be 0) of the asymptotic expansion of `F_{r\alpha}` as `r \to \infty`;
        ``exp_scale**r`` is the exponential factor of ``asy``;
        ``subexp_part`` is the subexponential factor of ``asy``.

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing

        A smooth point example::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: H = (1 - x - y - x*y)**2
            sage: Hfac = H.factor()
            sage: G = 1/Hfac.unit()
            sage: F = FFPD(G, Hfac); print(F)
            (1, [(x*y + x + y - 1, 2)])
            sage: alpha = [4, 3]
            sage: decomp = F.asymptotic_decomposition(alpha); decomp
            (0, []) + (... - 1/2, [(x*y + x + y - 1, 1)])
            sage: F1 = decomp[1]
            sage: p = {y: 1/3, x: 1/2}
            sage: asy = F1.asymptotics(p, alpha, 2, verbose=True)
            Creating auxiliary functions...
            Computing derivatives of auxiliary functions...
            Computing derivatives of more auxiliary functions...
            Computing second order differential operator actions...
            sage: asy
            (1/6000*(3600*sqrt(5)*sqrt(3)*sqrt(2)*sqrt(r)/sqrt(pi)
              + 463*sqrt(5)*sqrt(3)*sqrt(2)/(sqrt(pi)*sqrt(r)))*432^r,
             432,
             3/5*sqrt(5)*sqrt(3)*sqrt(2)*sqrt(r)/sqrt(pi)
              + 463/6000*sqrt(5)*sqrt(3)*sqrt(2)/(sqrt(pi)*sqrt(r)))
            sage: F.relative_error(asy[0], alpha, [1, 2, 4, 8, 16], asy[1])  # abs tol 1e-10  # long time
            [((4, 3), 2.083333333, [2.092576110], [-0.004436533009]),
             ((8, 6), 2.787374614, [2.790732875], [-0.001204811281]),
             ((16, 12), 3.826259447, [3.827462310], [-0.0003143703383]),
             ((32, 24), 5.328112821, [5.328540787], [-0.00008032230388]),
             ((64, 48), 7.475927885, [7.476079664], [-0.00002030232879])]

        A multiple point example::

            sage: R.<x,y,z>= PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: H = (4 - 2*x - y - z)**2*(4 - x - 2*y - z)
            sage: Hfac = H.factor()
            sage: G = 16/Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: F
            (-16, [(x + 2*y + z - 4, 1), (2*x + y + z - 4, 2)])
            sage: alpha = [3, 3, 2]
            sage: decomp = F.asymptotic_decomposition(alpha); decomp
            (0, []) + (..., [(x + 2*y + z - 4, 1), (2*x + y + z - 4, 1)])
            sage: F1 = decomp[1]
            sage: p = {x: 1, y: 1, z: 1}
            sage: asy = F1.asymptotics(p, alpha, 2, verbose=True) # long time
            Creating auxiliary functions...
            Computing derivatives of auxiliary functions...
            Computing derivatives of more auxiliary functions...
            Computing second-order differential operator actions...
            sage: asy # long time
            (4/3*sqrt(3)*sqrt(r)/sqrt(pi) + 47/216*sqrt(3)/(sqrt(pi)*sqrt(r)),
             1, 4/3*sqrt(3)*sqrt(r)/sqrt(pi) + 47/216*sqrt(3)/(sqrt(pi)*sqrt(r)))
            sage: F.relative_error(asy[0], alpha, [1, 2, 4, 8], asy[1]) # long time
            [((3, 3, 2), 0.9812164307, [1.515572606], [-0.54458543...]),
             ((6, 6, 4), 1.576181132, [1.992989399], [-0.26444185...]),
             ((12, 12, 8), 2.485286378, [2.712196351], [-0.091301338...]),
             ((24, 24, 16), 3.700576827, [3.760447895], [-0.016178847...])]
        """
        R = self.denominator_ring

        # Coerce keys of p into R.
        p = coerce_point(R, p)

        if asy_var is None:
            asy_var = var('r')
        d = self.dimension()
        X = list(R.gens())
        alpha = list(alpha)
        df = self.denominator_factored()
        n = len(df)     # Number of smooth factors

        # Find greatest i such that X[i] is a convenient coordinate,
        # that is, such that for all (h, e) in df, we have
        # (X[i]*diff(h, X[i])).subs(p) != 0.
        # Assuming such an i exists.
        i = d - 1
        while 0 in [(X[i] * diff(h, X[i])).subs(p) for (h, e) in df]:
            i -= 1
        coordinate = i

        if n == 1:
            # Smooth point.
            return self.asymptotics_smooth(p, alpha, N, asy_var, coordinate,
                                           numerical, verbose=verbose)

        # Multiple point.
        return self.asymptotics_multiple(p, alpha, N, asy_var, coordinate,
                                         numerical, verbose=verbose)

    def asymptotics_smooth(self, p, alpha, N, asy_var, coordinate=None,
                           numerical=0, verbose=False):
        r"""
        Return the asymptotics in the given direction of a smooth point.

        This is the same as :meth:`asymptotics()`, but only in the
        case of a convenient smooth point.

        The formulas used for computing the asymptotic expansions are
        Theorems 3.2 and 3.3 [RW2008]_ with the exponent of `H`
        equal to 1. Theorem 3.2 is a specialization of Theorem 3.4
        of [RW2012]_ with `n = 1`.

        INPUT:

        - ``p`` -- dictionary with keys that can be coerced to equal
          ``self.denominator_ring.gens()``
        - ``alpha`` -- tuple of length ``d = self.dimension()`` of positive
          integers or, if `p` is a smooth point, possibly of symbolic variables
        - ``N`` -- positive integer
        - ``asy_var`` -- (default: ``None``) a symbolic variable; the variable
          of the asymptotic expansion, if none is given, ``var('r')`` will be
          assigned
        - ``coordinate`` -- (default: ``None``) an integer in
          `\{0, \ldots, d-1\}` indicating a convenient coordinate to base the
          asymptotic calculations on; if ``None`` is assigned, then choose
          ``coordinate=d-1``
        - ``numerical`` -- (default: 0) a natural number; if numerical is
          greater than 0, then return a numerical approximation of the
          Maclaurin ray coefficients of ``self`` with ``numerical`` digits of
          precision; otherwise return exact values
        - ``verbose`` -- boolean (default: ``False``); print the current state
          of the algorithm

        OUTPUT: the asymptotic expansion

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: H = 2 - 3*x
            sage: Hfac = H.factor()
            sage: G = 1/Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: F
            (-1/3, [(x - 2/3, 1)])
            sage: alpha = [2]
            sage: p = {x: 2/3}
            sage: asy = F.asymptotics_smooth(p, alpha, 3, asy_var=var('r'))
            sage: asy
            (1/2*(9/4)^r, 9/4, 1/2)

        ::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: H = 1-x-y-x*y
            sage: Hfac = H.factor()
            sage: G = 1/Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: alpha = [3, 2]
            sage: p = {y: 1/2*sqrt(13) - 3/2, x: 1/3*sqrt(13) - 2/3}
            sage: F.asymptotics_smooth(p, alpha, 2, var('r'), numerical=3, verbose=True)
            Creating auxiliary functions...
            Computing derivatives of auxiliary functions...
            Computing derivatives of more auxiliary functions...
            Computing second order differential operator actions...
            (71.2^r*(0.369/sqrt(r) - 0.018.../r^(3/2)), 71.2, 0.369/sqrt(r) - 0.018.../r^(3/2))

            sage: q = 1/2
            sage: qq = q.denominator()
            sage: H = 1 - q*x + q*x*y - x^2*y
            sage: Hfac = H.factor()
            sage: G = (1 - q*x)/Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: alpha = list(qq*vector([2, 1 - q]))
            sage: alpha
            [4, 1]
            sage: p = {x: 1, y: 1}
            sage: F.asymptotics_smooth(p, alpha, 5, var('r'), verbose=True) # not tested (140 seconds)
            Creating auxiliary functions...
            Computing derivatives of auxiliary functions...
            Computing derivatives of more auxiliary functions...
            Computing second order differential operator actions...
            (1/12*sqrt(3)*2^(2/3)*gamma(1/3)/(pi*r^(1/3))
              - 1/96*sqrt(3)*2^(1/3)*gamma(2/3)/(pi*r^(5/3)),
             1,
             1/12*sqrt(3)*2^(2/3)*gamma(1/3)/(pi*r^(1/3))
              - 1/96*sqrt(3)*2^(1/3)*gamma(2/3)/(pi*r^(5/3)))
        """
        from sage.calculus.functions import jacobian
        from sage.calculus.var import function
        from sage.functions.other import factorial
        from sage.misc.functional import sqrt
        from sage.functions.gamma import gamma
        from sage.functions.log import exp, log
        from sage.matrix.constructor import matrix
        from sage.modules.free_module_element import vector
        from sage.symbolic.constants import pi
        from sage.symbolic.relation import solve
        from sage.rings.cc import CC
        from sage.rings.rational_field import QQ

        R = self.denominator_ring
        d = self.dimension()
        I = sqrt(-ZZ.one())
        # Coerce everything into the Symbolic Ring.
        X = [SR(x) for x in R.gens()]
        G = SR(self.numerator())
        H = SR(self.denominator())
        p = {SR(x): p[x] for x in R.gens()}
        alpha = [SR(a) for a in alpha]

        # Put given convenient coordinate at end of variable list.
        if coordinate is not None:
            x = X.pop(coordinate)
            X.append(x)
            a = alpha.pop(coordinate)
            alpha.append(a)

        # Deal with the simple univariate case first.
        # Same as the multiple point case with n == d.
        # but with a negative sign.
        # I'll just past the code from the multiple point case.
        if d == 1:
            det = jacobian(H, X).subs(p).determinant().abs()
            exp_scale = prod([(p[X[i]] ** (-alpha[i])).subs(p)
                              for i in range(d)])
            subexp_part = -G.subs(p) / (det * prod(p.values()))
            if numerical:
                exp_scale = exp_scale.n(digits=numerical)
                subexp_part = subexp_part.n(digits=numerical)
            return (exp_scale ** asy_var * subexp_part, exp_scale, subexp_part)

        # If p is a tuple of rationals, then compute with it directly.
        # Otherwise, compute symbolically and plug in p at the end.
        if vector(p.values()) in QQ ** d:
            P = p
        else:
            sP = [var('p' + str(j)) for j in range(d)]
            P = {X[j]: sP[j] for j in range(d)}
            p = {sP[j]: p[X[j]] for j in range(d)}

        # Setup.
        if verbose:
            print("Creating auxiliary functions...")
        # Implicit functions.
        h = function('h')(*tuple(X[:d - 1]))
        U = function('U')(*tuple(X))
        # All other functions are defined in terms of h, U, and
        # explicit functions.
        Gcheck = -G / U * (h / X[d - 1])
        A = Gcheck.subs({X[d - 1]: ZZ.one() / h}) / h
        t = 't'
        L = [str(elt) for elt in X]
        while t in L:
            t = t + 't'
        T = [var(t + str(i)) for i in range(d - 1)]
        e = {X[i]: P[X[i]] * exp(I * T[i]) for i in range(d - 1)}
        ht = h.subs(e)
        At = A.subs(e)
        Phit = (-log(P[X[d - 1]] * ht) +
                I * sum([alpha[i] / alpha[d - 1] * T[i]
                         for i in range(d - 1)]))
        Tstar = {t: ZZ.zero() for t in T}
        # Store h and U and all their derivatives evaluated at P.
        atP = P.copy()
        atP.update({h.subs(P): ZZ.one() / P[X[d - 1]]})

        # Compute the derivatives of h up to order 2 * N, evaluate at P,
        # and store in atP.
        # Keep a copy of unevaluated h derivatives for use in the case
        # d = 2 and v > 2 below.
        hderivs1 = {}   # First derivatives of h.
        for i in range(d - 1):
            s = solve(diff(H.subs({X[d - 1]: ZZ.one() / h}), X[i]),
                      diff(h, X[i]))[0].rhs().simplify()
            hderivs1.update({diff(h, X[i]): s})
            atP.update({diff(h, X[i]).subs(P): s.subs(P).subs(atP)})
        hderivs = diff_all(h, X[0: d - 1], 2 * N, sub=hderivs1, rekey=h)
        for k in hderivs:
            atP.update({k.subs(P): hderivs[k].subs(atP)})

        # Compute the derivatives of U up to order 2 * N and evaluate at P.
        # To do this, differentiate H = U*Hcheck over and over, evaluate at P,
        # and solve for the derivatives of U at P.
        # Need the derivatives of H with short keys to pass on
        # to diff_prod later.
        Hderivs = diff_all(H, X, 2 * N, ending=[X[d - 1]], sub_final=P)
        if verbose:
            print("Computing derivatives of auxiliary functions...")
        # For convenience in checking if all the nontrivial derivatives of U
        # at p are zero a few line below, store the value of U(p) in atP
        # instead of in Uderivs.
        Uderivs = {}
        atP.update({U.subs(P): diff(H, X[d - 1]).subs(P)})
        end = [X[d - 1]]
        Hcheck = X[d - 1] - ZZ.one() / h
        k = H.polynomial(CC).degree() - 1
        if k == 0:
            # Then we can conclude that all higher derivatives of U are zero.
            for l in range(1, 2 * N + 1):
                for s in combinations_with_replacement(X, l):
                    Uderivs[diff(U, list(s)).subs(P)] = ZZ.zero()
        elif k > 0 and k < 2 * N:
            all_zero = True
            Uderivs = diff_prod(Hderivs, U, Hcheck, X,
                                range(1, k + 1), end, Uderivs, atP)
            # Check for a nonzero U derivative.
            if any(Uderivs.values()):
                all_zero = False
            if all_zero:
                # Then, using a proposition at the end of [RW2012], we can
                # conclude that all higher derivatives of U are zero.
                for l in range(k + 1, 2 * N + 1):
                    for s in combinations_with_replacement(X, l):
                        Uderivs.update({diff(U, list(s)).subs(P): ZZ.zero()})
            else:
                # Have to compute the rest of the derivatives.
                Uderivs = diff_prod(Hderivs, U, Hcheck, X,
                                    range(k + 1, 2 * N + 1), end, Uderivs, atP)
        else:
            Uderivs = diff_prod(Hderivs, U, Hcheck, X,
                                range(1, 2 * N + 1), end, Uderivs, atP)
        atP.update(Uderivs)

        # In general, this algorithm is not designed to handle the case of a
        # singular Phit''(Tstar).
        # However, when d = 2 the algorithm can cope.
        if d == 2:
            # Compute v, the order of vanishing at Tstar of Phit.
            # It is at least 2.
            v = Integer(2)
            Phitderiv = diff(Phit, T[0], 2)
            splat = Phitderiv.subs(Tstar).subs(atP).subs(p).simplify()
            while splat == 0:
                v += 1
                if v > 2 * N:
                    # Then need to compute more derivatives of h for atP.
                    hderivs.update({diff(h, X[0], v):
                                    diff(hderivs[diff(h, X[0], v - 1)],
                                    X[0]).subs(hderivs1)})
                    atP.update({diff(h, X[0], v).subs(P):
                                hderivs[diff(h, X[0], v)].subs(atP)})
                Phitderiv = diff(Phitderiv, T[0])
                splat = Phitderiv.subs(Tstar).subs(atP).subs(p).simplify()

        if d == 2 and v > 2:
            t = T[0]  # Simplify variable names.
            a = splat / factorial(v)
            Phitu = Phit - a * t ** v

            # Compute all partial derivatives of At and Phitu
            # up to orders 2*(N - 1) and 2*(N - 1) + v, respectively,
            # in case v is even.
            # Otherwise, compute up to orders N - 1 and N - 1 + v,
            # respectively.
            # To speed up later computations,
            # create symbolic functions AA and BB
            # to stand in for the expressions At and Phitu, respectively.
            if verbose:
                print("Computing derivatives of more auxiliary functions...")
            AA = function('AA')(t)
            BB = function('BB')(t)
            if v.mod(2) == 0:
                At_derivs = diff_all(At, T, 2 * N - 2, sub=hderivs1,
                                     sub_final=[Tstar, atP], rekey=AA)
                Phitu_derivs = diff_all(Phitu, T, 2 * N - 2 + v,
                                        sub=hderivs1, sub_final=[Tstar, atP],
                                        zero_order=v + 1, rekey=BB)
            else:
                At_derivs = diff_all(At, T, N - 1, sub=hderivs1,
                                     sub_final=[Tstar, atP], rekey=AA)
                Phitu_derivs = diff_all(Phitu, T, N - 1 + v,
                                        sub=hderivs1, sub_final=[Tstar, atP],
                                        zero_order=v + 1, rekey=BB)
            AABB_derivs = At_derivs
            AABB_derivs.update(Phitu_derivs)
            AABB_derivs[AA] = At.subs(Tstar).subs(atP)
            AABB_derivs[BB] = Phitu.subs(Tstar).subs(atP)
            if verbose:
                print("Computing second order differential operator actions...")
            DD = diff_op_simple(AA, BB, AABB_derivs, t, v, a, N)

            # Plug above into asymptotic formula.
            L = []
            if v.mod(2) == 0:
                for k in range(N):
                    L.append(sum([(-1) ** l * gamma((2 * k + v * l + 1) / v) /
                                  (factorial(l) * factorial(2 * k + v * l)) *
                                  DD[(k, l)] for l in range(0, 2 * k + 1)]))
                chunk = (a ** (-1 / v) / (pi * v) *
                         sum([alpha[d - 1] ** (-(2 * k + 1) / v) *
                              L[k] * asy_var ** (-(2 * k + 1) / v)
                              for k in range(N)]))
            else:
                zeta = exp(I * pi / (2 * v))
                for k in range(N):
                    L.append(sum([(-1) ** l * gamma((k + v * l + 1) / v) /
                                  (factorial(l) * factorial(k + v * l)) *
                                  (zeta ** (k + v * l + 1) +
                                   (-1) ** (k + v * l) *
                                   zeta ** (-(k + v * l + 1))) *
                                  DD[(k, l)] for l in range(0, k + 1)]))
                chunk = (abs(a) ** (-1 / v) / (2 * pi * v) *
                         sum([alpha[d - 1] ** (-(k + 1) / v) *
                              L[k] * asy_var ** (-(k + 1) / v)
                              for k in range(N)]))

        # Asymptotics for d >= 2 case.
        # A singular Phit''(Tstar) will cause a crash in this case.
        else:
            Phit1 = jacobian(Phit, T).subs(hderivs1)
            a = jacobian(Phit1, T).subs(hderivs1).subs(Tstar).subs(atP)
            a_inv = a.inverse()
            Phitu = (Phit - (1 / QQ(2)) * matrix([T]) *
                     a * matrix([T]).transpose())
            Phitu = Phitu[0][0]
            # Compute all partial derivatives of At and Phitu up to
            # orders 2 * N-2 and 2 * N, respectively.
            # Take advantage of the fact that At and Phitu
            # are sufficiently differentiable functions so that mixed partials
            # are equal.  Thus only need to compute representative partials.
            # Choose nondecreasing sequences as representative differentiation-
            # order sequences.
            # To speed up later computations,
            # create symbolic functions AA and BB
            # to stand in for the expressions At and Phitu, respectively.
            if verbose:
                print("Computing derivatives of more auxiliary functions...")
            AA = function('AA')(*tuple(T))
            At_derivs = diff_all(At, T, 2 * N - 2, sub=hderivs1,
                                 sub_final=[Tstar, atP], rekey=AA)
            BB = function('BB')(*tuple(T))
            Phitu_derivs = diff_all(Phitu, T, 2 * N, sub=hderivs1,
                                    sub_final=[Tstar, atP], rekey=BB, zero_order=3)
            AABB_derivs = At_derivs
            AABB_derivs.update(Phitu_derivs)
            AABB_derivs[AA] = At.subs(Tstar).subs(atP)
            AABB_derivs[BB] = Phitu.subs(Tstar).subs(atP)
            if verbose:
                print("Computing second order differential operator actions...")
            DD = diff_op(AA, BB, AABB_derivs, T, a_inv, 1, N)

            # Plug above into asymptotic formula.
            L = []
            for k in range(N):
                L.append(sum([DD[(0, k, l)] / ((-1) ** k * 2 ** (l + k) *
                                               factorial(l) * factorial(l + k))
                              for l in range(0, 2 * k + 1)]))
            chunk = sum([(2 * pi) ** ((1 - d) / Integer(2)) *
                         a.determinant() ** (-ZZ.one() / Integer(2)) *
                         alpha[d - 1] ** ((ZZ.one() - d) / Integer(2) - k) *
                         L[k] *
                         asy_var ** ((ZZ.one() - d) / Integer(2) - k)
                         for k in range(N)])

        chunk = chunk.subs(p).simplify()
        coeffs = chunk.coefficients(asy_var)
        coeffs.reverse()
        coeffs = coeffs[:N]
        if numerical:
            subexp_part = sum([co[0].subs(p).n(digits=numerical) *
                               asy_var ** co[1] for co in coeffs])
            exp_scale = prod([(P[X[i]] ** (-alpha[i])).subs(p)
                              for i in range(d)]).n(digits=numerical)
        else:
            subexp_part = sum([co[0].subs(p) * asy_var ** co[1]
                               for co in coeffs])
            exp_scale = prod([(P[X[i]] ** (-alpha[i])).subs(p)
                              for i in range(d)])
        return (exp_scale ** asy_var * subexp_part, exp_scale, subexp_part)

    def asymptotics_multiple(self, p, alpha, N, asy_var, coordinate=None,
                             numerical=0, verbose=False):
        r"""
        Return the asymptotics in the given direction of a multiple
        point nondegenerate for ``alpha``.

        This is the same as :meth:`asymptotics`, but only in the case
        of a convenient multiple point nondegenerate for ``alpha``.
        Assume also that ``self.dimension >= 2`` and that the
        ``p.values()`` are not symbolic variables.

        The formulas used for computing the asymptotic expansion are
        Theorem 3.4 and Theorem 3.7 of [RW2012]_.

        INPUT:

        - ``p`` -- dictionary with keys that can be coerced to equal
          ``self.denominator_ring.gens()``
        - ``alpha`` -- tuple of length ``d = self.dimension()`` of positive
          integers or, if `p` is a smooth point, possibly of symbolic variables
        - ``N`` -- positive integer
        - ``asy_var`` -- (default: ``None``) a symbolic variable; the variable
          of the asymptotic expansion, if none is given, ``var('r')`` will be
          assigned
        - ``coordinate`` -- (default: ``None``) an integer in
          `\{0, \ldots, d-1\}` indicating a convenient coordinate to base the
          asymptotic calculations on; if ``None`` is assigned, then choose
          ``coordinate=d-1``
        - ``numerical`` -- (default: 0) a natural number; if numerical is
          greater than 0, then return a numerical approximation of the
          Maclaurin ray coefficients of ``self`` with ``numerical`` digits of
          precision. Otherwise return exact values.
        - ``verbose`` -- boolean (default: ``False``); print the current state of
          the algorithm

        OUTPUT: the asymptotic expansion

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y,z>= PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: H = (4 - 2*x - y - z)*(4 - x -2*y - z)
            sage: Hfac = H.factor()
            sage: G = 16/Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: F
            (16, [(x + 2*y + z - 4, 1), (2*x + y + z - 4, 1)])
            sage: p = {x: 1, y: 1, z: 1}
            sage: alpha = [3, 3, 2]
            sage: F.asymptotics_multiple(p, alpha, 2, var('r'), verbose=True) # long time
            Creating auxiliary functions...
            Computing derivatives of auxiliary functions...
            Computing derivatives of more auxiliary functions...
            Computing second-order differential operator actions...
            (4/3*sqrt(3)/(sqrt(pi)*sqrt(r)) - 25/216*sqrt(3)/(sqrt(pi)*r^(3/2)),
             1,
             4/3*sqrt(3)/(sqrt(pi)*sqrt(r)) - 25/216*sqrt(3)/(sqrt(pi)*r^(3/2)))

            sage: H = (1 - x*(1 + y))*(1 - z*x**2*(1 + 2*y))
            sage: Hfac = H.factor()
            sage: G = 1/Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: F
            (1, [(x*y + x - 1, 1), (2*x^2*y*z + x^2*z - 1, 1)])
            sage: p = {x: 1/2, z: 4/3, y: 1}
            sage: alpha = [8, 3, 3]
            sage: F.asymptotics_multiple(p, alpha, 2, var('r'), coordinate=1, verbose=True) # long time
            Creating auxiliary functions...
            Computing derivatives of auxiliary functions...
            Computing derivatives of more auxiliary functions...
            Computing second-order differential operator actions...
            (1/172872*108^r*(24696*sqrt(7)*sqrt(3)/(sqrt(pi)*sqrt(r))
              - 1231*sqrt(7)*sqrt(3)/(sqrt(pi)*r^(3/2))),
             108,
             1/7*sqrt(7)*sqrt(3)/(sqrt(pi)*sqrt(r))
              - 1231/172872*sqrt(7)*sqrt(3)/(sqrt(pi)*r^(3/2)))

        ::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: H = (1 - 2*x - y) * (1 - x - 2*y)
            sage: Hfac = H.factor()
            sage: G = exp(x + y)/Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: F
            (e^(x + y), [(x + 2*y - 1, 1), (2*x + y - 1, 1)])
            sage: p = {x: 1/3, y: 1/3}
            sage: alpha = (var('a'), var('b'))
            sage: F.asymptotics_multiple(p, alpha, 2, var('r')) # long time
            (3*(1/((1/3)^a*(1/3)^b))^r*e^(2/3), 1/((1/3)^a*(1/3)^b), 3*e^(2/3))
        """
        from itertools import product
        from sage.calculus.functions import jacobian
        from sage.calculus.var import function
        from sage.combinat.combinat import stirling_number1
        from sage.functions.log import exp, log
        from sage.functions.other import factorial
        from sage.misc.functional import sqrt
        from sage.matrix.constructor import matrix
        from sage.misc.mrange import xmrange
        from sage.modules.free_module_element import vector
        from sage.rings.cc import CC
        from sage.arith.misc import binomial
        from sage.rings.rational_field import QQ
        from sage.symbolic.constants import pi
        from sage.symbolic.relation import solve

        R = self.denominator_ring

        # Coerce keys of p into R.
        p = coerce_point(R, p)

        d = self.dimension()
        I = sqrt(-ZZ.one())
        # Coerce everything into the Symbolic Ring.
        X = [SR(x) for x in R.gens()]
        G = SR(self.numerator())
        H = [SR(h) for (h, e) in self.denominator_factored()]
        Hprod = prod(H)
        n = len(H)
        P = {SR(x): p[x] for x in R.gens()}
        Sstar = self._crit_cone_combo(p, alpha, coordinate)

        # Put the given convenient variable at end of variable list.
        if coordinate is not None:
            x = X.pop(coordinate)
            X.append(x)
            a = alpha.pop(coordinate)
            alpha.append(a)

        # Case n = d.
        if n == d:
            det = jacobian(H, X).subs(P).determinant().abs()
            exp_scale = prod([(P[X[i]] ** (-alpha[i])).subs(P)
                              for i in range(d)])
            subexp_part = G.subs(P) / (det * prod(P.values()))
            if numerical:
                exp_scale = exp_scale.n(digits=numerical)
                subexp_part = subexp_part.n(digits=numerical)
            return (exp_scale ** asy_var * subexp_part, exp_scale, subexp_part)

        # Case n < d.
        # If P is a tuple of rationals, then compute with it directly.
        # Otherwise, compute symbolically and plug in P at the end.
        if vector(P.values()) not in QQ ** d:
            sP = [var('p' + str(j)) for j in range(d)]
            P = {X[j]: sP[j] for j in range(d)}
            p = {sP[j]: p[X[j]] for j in range(d)}

        # Setup.
        if verbose:
            print("Creating auxiliary functions...")
        # Create T and S variables.
        t = 't'
        L = [str(elt) for elt in X]
        while t in L:
            t = t + 't'
        T = [var(t + str(i)) for i in range(d - 1)]
        s = 's'
        while s in L:
            s = s + 't'
        S = [var(s + str(i)) for i in range(n - 1)]
        Sstar = {S[j]: Sstar[j] for j in range(n - 1)}
        thetastar = {t: ZZ.zero() for t in T}
        thetastar.update(Sstar)
        # Create implicit functions.
        h = [function('h' + str(j))(*tuple(X[:d - 1])) for j in range(n)]
        U = function('U')(*tuple(X))
        # All other functions are defined in terms of h, U, and
        # explicit functions.
        Hcheck = prod([X[d - 1] - ZZ.one() / h[j] for j in range(n)])
        Gcheck = -G / U * prod([-h[j] / X[d - 1] for j in range(n)])
        A = [(-1) ** (n - 1) * X[d - 1] ** (-n + j) *
             diff(Gcheck.subs({X[d - 1]: ZZ.one() / X[d - 1]}), X[d - 1], j)
             for j in range(n)]
        e = {X[i]: P[X[i]] * exp(I * T[i]) for i in range(d - 1)}
        ht = [hh.subs(e) for hh in h]
        hsumt = (sum([S[j] * ht[j] for j in range(n - 1)]) +
                 (ZZ.one() - sum(S)) * ht[n - 1])
        At = [AA.subs(e).subs({X[d - 1]: hsumt}) for AA in A]
        Phit = (-log(P[X[d - 1]] * hsumt) +
                I * sum([alpha[i] / alpha[d - 1] * T[i]
                         for i in range(d - 1)]))
        # atP Stores h and U and all their derivatives evaluated at C.
        atP = P.copy()
        atP.update({hh.subs(P): ZZ.one() / P[X[d - 1]] for hh in h})

        # Compute the derivatives of h up to order 2 * N and evaluate at P.
        hderivs1 = {}   # First derivatives of h.
        for (i, j) in xmrange([d - 1, n], tuple):
            s = solve(diff(H[j].subs({X[d - 1]: ZZ.one() / h[j]}), X[i]),
                      diff(h[j], X[i]))[0].rhs().simplify()
            hderivs1.update({diff(h[j], X[i]): s})
            atP.update({diff(h[j], X[i]).subs(P): s.subs(P).subs(atP)})
        hderivs = diff_all(h, X[0:d - 1], 2 * N, sub=hderivs1, rekey=h)
        for k in hderivs:
            atP.update({k.subs(P): hderivs[k].subs(atP)})

        # Compute the derivatives of U up to order 2 * N - 2 + min{n, N} - 1 and
        # evaluate at P.
        # To do this, differentiate H = U*Hcheck over and over, evaluate at P,
        # and solve for the derivatives of U at P.
        # Need the derivatives of H with short keys to pass on to
        # diff_prod later.
        if verbose:
            print("Computing derivatives of auxiliary functions...")
        m = min(n, N)
        end = [X[d - 1] for j in range(n)]
        Hprodderivs = diff_all(Hprod, X, 2 * N - 2 + n, ending=end, sub_final=P)
        atP.update({U.subs(P): diff(Hprod, X[d - 1], n).subs(P) / factorial(n)})
        Uderivs = {}
        k = Hprod.polynomial(CC).degree() - n
        if k == 0:
            # Then we can conclude that all higher derivatives of U are zero.
            for l in range(1, 2 * N - 2 + m):
                for s in combinations_with_replacement(X, l):
                    Uderivs[diff(U, list(s)).subs(P)] = ZZ.zero()
        elif k > 0 and k < 2 * N - 2 + m - 1:
            all_zero = True
            Uderivs = diff_prod(Hprodderivs, U, Hcheck, X,
                                range(1, k + 1), end, Uderivs, atP)
            # Check for a nonzero U derivative.
            if any(Uderivs.values()):
                all_zero = False
            if all_zero:
                # Then all higher derivatives of U are zero.
                for l in range(k + 1, 2 * N - 2 + m):
                    for s in combinations_with_replacement(X, l):
                        Uderivs.update({diff(U, list(s)).subs(P): ZZ.zero()})
            else:
                # Have to compute the rest of the derivatives.
                Uderivs = diff_prod(Hprodderivs, U, Hcheck, X,
                                    range(k + 1, 2 * N - 2 + m), end,
                                    Uderivs, atP)
        else:
            Uderivs = diff_prod(Hprodderivs, U, Hcheck, X,
                                range(1, 2 * N - 2 + m), end, Uderivs, atP)
        atP.update(Uderivs)
        Phit1 = jacobian(Phit, T + S).subs(hderivs1)
        a = jacobian(Phit1, T + S).subs(hderivs1).subs(thetastar).subs(atP)
        a_inv = a.inverse()
        Phitu = (Phit - (1 / Integer(2)) * matrix([T + S]) * a *
                 matrix([T + S]).transpose())
        Phitu = Phitu[0][0]

        # Compute all partial derivatives of At and Phitu up to orders 2 * N - 2
        # and 2 * N, respectively.  Take advantage of the fact that At and Phitu
        # are sufficiently differentiable functions so that mixed partials
        # are equal.  Thus only need to compute representative partials.
        # Choose nondecreasing sequences as representative differentiation-
        # order sequences.
        # To speed up later computations, create symbolic functions AA and BB
        # to stand in for the expressions At and Phitu respectively.
        if verbose:
            print("Computing derivatives of more auxiliary functions...")
        AA = [function('A' + str(j))(*tuple(T + S)) for j in range(n)]
        At_derivs = diff_all(At, T + S, 2 * N - 2, sub=hderivs1,
                             sub_final=[thetastar, atP], rekey=AA)
        BB = function('BB')(*tuple(T + S))
        Phitu_derivs = diff_all(Phitu, T + S, 2 * N, sub=hderivs1,
                                sub_final=[thetastar, atP], rekey=BB, zero_order=3)
        AABB_derivs = At_derivs
        AABB_derivs.update(Phitu_derivs)
        for j in range(n):
            AABB_derivs[AA[j]] = At[j].subs(thetastar).subs(atP)
        AABB_derivs[BB] = Phitu.subs(thetastar).subs(atP)

        if verbose:
            print("Computing second-order differential operator actions...")
        DD = diff_op(AA, BB, AABB_derivs, T + S, a_inv, n, N)
        L = {}
        for (j, k) in product(range(min(n, N)), range(max(0, N - 1 - n), N)):
            if j + k <= N - 1:
                L[(j, k)] = sum([DD[(j, k, l)] / ((-1) ** k * 2 ** (k + l) *
                                                  factorial(l) *
                                                  factorial(k + l))
                                 for l in range(2 * k + 1)])
                det = (a.determinant() ** (-1 / Integer(2)) *
                       (2 * pi) ** ((n - d) / Integer(2)))
        chunk = det * sum([(alpha[d - 1] * asy_var) ** ((n - d) /
                                                        Integer(2) - q) *
                           sum([L[(j, k)] * binomial(n - 1, j) *
                                stirling_number1(n - j, n + k - q) *
                                (-1) ** (q - j - k)
                                for (j, k) in product(range(min(n - 1, q) + 1),
                                                      range(max(0, q - n),
                                                            q + 1))
                                if j + k <= q])
                           for q in range(N)])
        chunk = chunk.subs(P).simplify()
        coeffs = chunk.coefficients(asy_var)
        coeffs.reverse()
        coeffs = coeffs[:N]
        if numerical:
            subexp_part = sum([co[0].subs(p).n(digits=numerical) *
                               asy_var ** co[1]
                               for co in coeffs])
            exp_scale = prod([(P[X[i]] ** (-alpha[i])).subs(p)
                              for i in range(d)]).n(digits=numerical)
        else:
            subexp_part = sum([co[0].subs(p) * asy_var ** co[1]
                               for co in coeffs])
            exp_scale = prod([(P[X[i]] ** (-alpha[i])).subs(p)
                              for i in range(d)])
        return (exp_scale ** asy_var * subexp_part, exp_scale, subexp_part)

    def _crit_cone_combo(self, p, alpha, coordinate=None):
        r"""
        Return an auxiliary point associated to the multiple
        point ``p`` of the factors ``self``.

        INPUT:

        - ``p`` -- dictionary with keys that can be coerced to equal
          ``self.denominator_ring.gens()``
        - ``alpha`` -- list of rationals

        OUTPUT:

        A solution of the matrix equation `y \Gamma = \alpha^{\prime}` for `y`,
        where `\Gamma` is the matrix given by
        ``[direction(v) for v in self.log_grads(p)]`` and
        `\alpha^{\prime}` is ``direction(alpha)``.

        .. SEEALSO::

            :func:`direction`

        .. NOTE::

            For internal use by
            :meth:`FractionWithFactoredDenominator.asymptotics_multiple()`.

        .. NOTE::

            Use this function only when `\Gamma` is well-defined and
            there is a unique solution to the matrix equation
            `y \Gamma = \alpha'`. Fails otherwise.

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: p = exp(x)
            sage: df = [(1 - 2*x - y, 1), (1 - x - 2*y, 1)]
            sage: f = FFPD(p, df)
            sage: p = {x: 1/3, y: 1/3}
            sage: alpha = (var('a'), var('b'))
            sage: f._crit_cone_combo(p, alpha)
            [1/3*(2*a - b)/b, -2/3*(a - 2*b)/b]
        """
        from sage.matrix.constructor import matrix
        from sage.symbolic.relation import solve

        # Assuming here that each log_grads(f) has nonzero final component.
        # Then 'direction' will not throw a division by zero error.
        R = self.denominator_ring

        # Coerce keys of p into R.
        p = coerce_point(R, p)

        d = self.dimension()
        n = len(self.denominator_factored())
        Gamma = matrix([direction(v, coordinate) for v in self.log_grads(p)])
        beta = direction(alpha, coordinate)
        # solve_left() fails when working in SR :-(.
        # So use solve() instead.
        # Gamma.solve_left(vector(beta))
        V = [var('sss' + str(i)) for i in range(n)]
        M = matrix(V) * Gamma
        eqns = [M[0][j] == beta[j] for j in range(d)]
        s = solve(eqns, V, solution_dict=True)[0]  # Assume a unique solution.
        return [s[v] for v in V]

    def grads(self, p):
        r"""
        Return a list of the gradients of the polynomials
        ``[q for (q, e) in self.denominator_factored()]`` evaluated at ``p``.

        INPUT:

        - ``p`` -- (default: ``None``) a dictionary whose keys are
          the generators of ``self.denominator_ring``

        OUTPUT: list

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: p = exp(x)
            sage: df = [(x^3 + 3*y^2, 5), (x*y, 2), (y, 1)]
            sage: f = FFPD(p, df)
            sage: f
            (e^x, [(y, 1), (x*y, 2), (x^3 + 3*y^2, 5)])
            sage: R.gens()
            (x, y)
            sage: p = None
            sage: f.grads(p)
            [(0, 1), (y, x), (3*x^2, 6*y)]

            sage: p = {x: sqrt(2), y: var('a')}
            sage: f.grads(p)
            [(0, 1), (a, sqrt(2)), (6, 6*a)]
        """
        R = self.denominator_ring

        # Coerce keys of p into R.
        p = coerce_point(R, p)

        X = R.gens()
        d = self.dimension()
        H = [h for (h, e) in self.denominator_factored()]
        n = len(H)
        return [tuple([diff(H[i], X[j]).subs(p) for j in range(d)])
                for i in range(n)]

    def log_grads(self, p):
        r"""
        Return a list of the logarithmic gradients of the polynomials
        ``[q for (q, e) in self.denominator_factored()]`` evaluated at ``p``.

        The logarithmic gradient of a function `f` at point `p` is the
        vector `(x_1 \partial_1 f(x), \ldots, x_d \partial_d f(x) )`
        evaluated at `p`.

        INPUT:

        - ``p`` -- (default: ``None``) a dictionary whose keys
          are the generators of ``self.denominator_ring``

        OUTPUT: list

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: p = exp(x)
            sage: df = [(x^3 + 3*y^2, 5), (x*y, 2), (y, 1)]
            sage: f = FFPD(p, df)
            sage: f
            (e^x, [(y, 1), (x*y, 2), (x^3 + 3*y^2, 5)])
            sage: R.gens()
            (x, y)
            sage: p = None
            sage: f.log_grads(p)
            [(0, y), (x*y, x*y), (3*x^3, 6*y^2)]

            sage: p = {x: sqrt(2), y: var('a')}
            sage: f.log_grads(p)
            [(0, a), (sqrt(2)*a, sqrt(2)*a), (6*sqrt(2), 6*a^2)]
        """
        R = self.denominator_ring

        # Coerce keys of p into R.
        p = coerce_point(R, p)

        X = R.gens()
        d = self.dimension()
        H = [h for (h, e) in self.denominator_factored()]
        n = len(H)
        return [tuple([(X[j] * diff(H[i], X[j])).subs(p) for j in range(d)])
                for i in range(n)]

    def critical_cone(self, p, coordinate=None):
        r"""
        Return the critical cone of the convenient multiple point ``p``.

        INPUT:

        - ``p`` -- dictionary with keys that can be coerced to equal
          ``self.denominator_ring.gens()`` and values in a field
        - ``coordinate`` -- (default: ``None``) a natural number

        OUTPUT: list of vectors

        This list of vectors generate the critical cone of ``p`` and
        the cone itself, which is ``None`` if the values of ``p`` don't lie in
        `\QQ`. Divide logarithmic gradients by their component ``coordinate``
        entries. If ``coordinate = None``, then search from `d-1` down to 0
        for the first index ``j`` such that for all ``i`` we have
        ``self.log_grads()[i][j] != 0`` and set ``coordinate = j``.

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y,z> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: G = 1
            sage: H = (1 - x*(1 + y)) * (1 - z*x**2*(1 + 2*y))
            sage: Hfac = H.factor()
            sage: G = 1/Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: p = {x: 1/2, y: 1, z: 4/3}
            sage: F.critical_cone(p)
            ([(2, 1, 0), (3, 1, 3/2)], 2-d cone in 3-d lattice N)
        """
        from sage.geometry.cone import Cone

        R = self.denominator_ring

        # Coerce keys of p into R.
        p = coerce_point(R, p)

        d = self.dimension()
        lg = self.log_grads(p)
        n = len(lg)
        if coordinate not in range(d):
            # Search from d-1 down to 0 for a coordinate j such that
            # for all i we have lg[i][j] != 0.
            # One is guaranteed to exist in the case of a convenient multiple
            # point.
            for j in reversed(range(d)):
                if 0 not in [lg[i][j] for i in range(n)]:
                    coordinate = j
                    break
        Gamma = [direction(v, coordinate) for v in lg]
        try:
            cone = Cone(Gamma)
        except TypeError:
            cone = None
        return (Gamma, cone)

    def is_convenient_multiple_point(self, p):
        r"""
        Test if ``p`` is a convenient multiple point of ``self``.

        In case ``p`` is a convenient multiple point, ``verdict = True`` and
        ``comment`` is a string stating which variables it's convenient to use.
        In case ``p`` is not, ``verdict = False`` and ``comment`` is a string
        explaining why ``p`` fails to be a convenient multiple point.

        See [RW2012]_ for more details.

        INPUT:

        - ``p`` -- dictionary with keys that can be coerced to equal
          ``self.denominator_ring.gens()``

        OUTPUT:

        A pair ``(verdict, comment)``.

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y,z> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: H = (1 - x*(1 + y)) * (1 - z*x**2*(1 + 2*y))
            sage: df = H.factor()
            sage: G = 1 / df.unit()
            sage: F = FFPD(G, df)
            sage: p1 = {x: 1/2, y: 1, z: 4/3}
            sage: p2 = {x: 1, y: 2, z: 1/2}
            sage: F.is_convenient_multiple_point(p1)
            (True, 'convenient in variables [x, y]')
            sage: F.is_convenient_multiple_point(p2)
            (False, 'not a singular point')
        """
        from sage.combinat.subset import Subsets
        from sage.matrix.constructor import matrix

        R = self.denominator_ring

        # Coerce keys of p into R.
        p = coerce_point(R, p)

        H = [h for (h, e) in self.denominator_factored()]
        n = len(H)
        d = self.dimension()

        # Test 1: Are the factors in H zero at p?
        if [h.subs(p) for h in H] != [0 for h in H]:
            # Failed test 1.  Move on to next point.
            return (False, 'not a singular point')

        # Test 2: Are the factors in H smooth at p?
        grads = self.grads(p)
        for v in grads:
            if v == [0 for i in range(d)]:
                return (False, 'not smooth point of factors')

        # Test 3: Do the factors in H intersect transversely at p?
        if n <= d:
            M = matrix(grads)
            if M.rank() != n:
                return (False, 'not a transverse intersection')
        else:
            # Check all sub-multisets of grads of size d.
            for S in Subsets(grads, d, submultiset=True):
                M = matrix(S)
                if M.rank() != d:
                    return (False, 'not a transverse intersection')

        # Test 4: Is p convenient?
        M = matrix(self.log_grads(p))
        cols = M.columns()
        convenient_coordinates = [j for j, c in enumerate(cols) if 0 not in c]
        if not convenient_coordinates:
            return (False, 'multiple point but not convenient')

        # Tests all passed
        X = R.gens()
        convenientX = [X[i] for i in convenient_coordinates]
        return (True, 'convenient in variables {}'.format(convenientX))

    def singular_ideal(self):
        r"""
        Return the singular ideal of ``self``.

        Let `R` be the ring of ``self`` and `H` its denominator.
        Let `H_{red}` be the reduction (square-free part) of `H`.
        Return the ideal in `R` generated by `H_{red}` and
        its partial derivatives.
        If the coefficient field of `R` is algebraically closed,
        then the output is the ideal of the singular locus (which
        is a variety) of the variety of `H`.

        OUTPUT: an ideal

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y,z> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: H = (1 - x*(1 + y))^3 * (1 - z*x**2*(1 + 2*y))
            sage: df = H.factor()
            sage: G = 1 / df.unit()
            sage: F = FFPD(G, df)
            sage: F.singular_ideal()
            Ideal (x*y + x - 1, y^2 - 2*y*z + 2*y - z + 1, x*z + y - 2*z + 1) of
             Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        R = self.denominator_ring

        Hred = prod([h for (h, e) in self.denominator_factored()])
        J = R.ideal([Hred] + Hred.gradient())
        return R.ideal(J.groebner_basis())

    def smooth_critical_ideal(self, alpha):
        r"""
        Return the smooth critical ideal of ``self``.

        Let `R` be the ring of ``self`` and `H` its denominator.
        Return the ideal in `R` of smooth critical points of the variety
        of `H` for the direction ``alpha``.
        If the variety `V` of `H` has no smooth points, then return the ideal
        in `R` of `V`.

        See [RW2012]_ for more details.

        INPUT:

        - ``alpha`` -- tuple of positive integers and/or symbolic entries
          of length ``self.denominator_ring.ngens()``

        OUTPUT: an ideal

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: H = (1 - x - y - x*y)^2
            sage: Hfac = H.factor()
            sage: G = 1/Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: alpha = var('a1, a2')
            sage: F.smooth_critical_ideal(alpha)
            Ideal (y^2 + (2*a1)/a2*y - 1, x + (-a2)/a1*y + (-a1 + a2)/a1) of
             Multivariate Polynomial Ring in x, y over Fraction Field of
             Multivariate Polynomial Ring in a1, a2 over Rational Field

            sage: H = (1-x-y-x*y)^2
            sage: Hfac = H.factor()
            sage: G = 1/Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: alpha = [7/3, var('a')]
            sage: F.smooth_critical_ideal(alpha)
            Ideal (y^2 + 14/(3*a)*y - 1, x + (-3*a)/7*y + (3*a - 7)/7) of Multivariate Polynomial Ring in x, y over Fraction Field of Univariate Polynomial Ring in a over Rational Field
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        R = self.denominator_ring
        Hred = prod([h for (h, e) in self.denominator_factored()])
        K = R.base_ring()
        d = self.dimension()

        # Expand K by the variables of alpha if there are any.
        indets = [a for a in alpha if a not in K and a in SR]

        indets = sorted(set(indets), key=str)   # Delete duplicates in indets.
        if indets:
            L = PolynomialRing(K, indets).fraction_field()
            S = R.change_ring(L)
            # Coerce alpha into L.
            alpha = [L(a) for a in alpha]
        else:
            S = R

        # Find smooth, critical points for alpha.
        X = S.gens()
        Hred = S(Hred)
        J = S.ideal([Hred] +
                    [alpha[d - 1] * X[i] * diff(Hred, X[i]) -
                     alpha[i] * X[d - 1] * diff(Hred, X[d - 1])
                     for i in range(d - 1)])
        return S.ideal(J.groebner_basis())

    def maclaurin_coefficients(self, multi_indices, numerical=0):
        r"""
        Return the Maclaurin coefficients of ``self`` with given
        ``multi_indices``.

        INPUT:

        - ``multi_indices`` -- list of tuples of positive integers, where
          each tuple has length ``self.dimension()``
        - ``numerical`` -- (default: 0) a natural number; if positive, return
          numerical approximations of coefficients with ``numerical`` digits of
          accuracy

        OUTPUT:

        A dictionary whose value of the key ``nu`` are the Maclaurin
        coefficient of index ``nu`` of ``self``.

        .. NOTE::

            Uses iterated univariate Maclaurin expansions. Slow.

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: H = 2 - 3*x
            sage: Hfac = H.factor()
            sage: G = 1 / Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: F
            (-1/3, [(x - 2/3, 1)])
            sage: F.maclaurin_coefficients([(2*k,) for k in range(6)])
            {(0,): 1/2,
             (2,): 9/8,
             (4,): 81/32,
             (6,): 729/128,
             (8,): 6561/512,
             (10,): 59049/2048}

        ::

            sage: R.<x,y,z> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: H = (4 - 2*x - y - z) * (4 - x - 2*y - z)
            sage: Hfac = H.factor()
            sage: G = 16 / Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: alpha = vector([3, 3, 2])
            sage: interval = [1, 2, 4]
            sage: S = [r*alpha for r in interval]
            sage: F.maclaurin_coefficients(S, numerical=10)  # long time
            {(3, 3, 2): 0.7849731445,
             (6, 6, 4): 0.7005249476,
             (12, 12, 8): 0.5847732654}
        """
        R = self.denominator_ring
        d = self.dimension()
        coeffs = {}

        # Deal with the simple univariate case first.
        if d == 1:
            f = SR(self.quotient())
            x = SR(R.gens()[0])
            m = max(multi_indices)[0]
            f = f.taylor(x, 0, m)
            F = R(f)
            tmp = F.coefficients()
            for nu in multi_indices:
                val = tmp[nu[0]]
                if numerical:
                    val = val.n(digits=numerical)
                coeffs[tuple(nu)] = val
            return coeffs

        # Create biggest multi-index needed.
        alpha = [max(nu[i] for nu in multi_indices) for i in range(d)]

        # Compute Maclaurin expansion of self up to index alpha.
        # Use iterated univariate expansions.
        # Slow!
        f = SR(self.quotient())
        X = [SR(g) for g in R.gens()]
        for i in range(d):
            f = f.taylor(X[i], 0, alpha[i])
        F = R(f)

        # Collect coefficients.
        X = R.gens()
        for nu in multi_indices:
            monomial = prod(X[i] ** nu[i] for i in range(d))
            val = F.monomial_coefficient(monomial)
            if numerical:
                val = val.n(digits=numerical)
            coeffs[tuple(nu)] = val
        return coeffs

    def relative_error(self, approx, alpha, interval, exp_scale=Integer(1),
                       digits=10):
        r"""
        Return the relative error between the values of the Maclaurin
        coefficients of ``self`` with multi-indices ``r alpha`` for ``r`` in
        ``interval`` and the values of the functions (of the variable ``r``)
        in ``approx``.

        INPUT:

        - ``approx`` -- an individual or list of symbolic expressions in
          one variable
        - ``alpha`` -- list of positive integers of length
          ``self.denominator_ring.ngens()``
        - ``interval`` -- list of positive integers
        - ``exp_scale`` -- (default: 1) a number

        OUTPUT: list of tuples with properties described below

        This outputs a list whose entries are a tuple
        ``(r*alpha, a_r, b_r, err_r)`` for ``r`` in ``interval``.
        Here ``r*alpha`` is a tuple; ``a_r`` is the ``r*alpha`` (multi-index)
        coefficient of the Maclaurin series for ``self`` divided by
        ``exp_scale**r``;
        ``b_r`` is a list of the values of the functions in ``approx``
        evaluated at ``r`` and divided by ``exp_scale**m``;
        ``err_r`` is the list of relative errors
        ``(a_r - f)/a_r`` for ``f`` in ``b_r``.
        All outputs are decimal approximations.

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: H = 1 - x - y - x*y
            sage: Hfac = H.factor()
            sage: G = 1 / Hfac.unit()
            sage: F = FFPD(G, Hfac)
            sage: alpha = [1, 1]
            sage: r = var('r')
            sage: a1 = (0.573/sqrt(r))*5.83^r
            sage: a2 = (0.573/sqrt(r) - 0.0674/r^(3/2))*5.83^r
            sage: es = 5.83
            sage: F.relative_error([a1, a2], alpha, [1, 2, 4, 8], es) # long time
            [((1, 1), 0.5145797599,
              [0.5730000000, 0.5056000000], [-0.1135300000, 0.01745066667]),
             ((2, 2), 0.3824778089,
              [0.4051721856, 0.3813426871], [-0.05933514614, 0.002967810973]),
             ((4, 4), 0.2778630595,
              [0.2865000000, 0.2780750000], [-0.03108344267, -0.0007627515584]),
             ((8, 8), 0.1991088276,
              [0.2025860928, 0.1996074055], [-0.01746414394, -0.002504047242])]
        """
        from sage.modules.free_module_element import vector

        if not isinstance(approx, (list, tuple)):
            approx = [approx]
        if approx[0].variables():
            av = approx[0].variables()[0]
        else:
            av = ZZ.one()

        # Get Maclaurin coefficients of self.
        alpha = vector(alpha)
        multi_indices = [r * alpha for r in interval]
        mac = self.maclaurin_coefficients(multi_indices, numerical=digits)
        # mac = self.old_maclaurin_coefficients(alpha, max(interval))
        mac_approx = {}
        stats = []
        for r in interval:
            exp_s_r = exp_scale ** r
            beta = tuple(r * alpha)
            mac[beta] = (mac[beta] / exp_s_r).n(digits=digits)
            mac_approx[beta] = [(f.subs({av: r}) / exp_s_r).n(digits=digits)
                                for f in approx]
            stats_row = [beta, mac[beta], mac_approx[beta]]
            if mac[beta] == 0:
                stats_row.extend([None for a in mac_approx[beta]])
            else:
                stats_row.append([(mac[beta] - a) / mac[beta]
                                  for a in mac_approx[beta]])
            stats.append(tuple(stats_row))
        return stats

    def _add_(left, right):
        r"""
        Return the sum of ``left`` with ``right``.

        INPUT:

        - ``left`` -- the left summand (i.e. ``self``)

        - ``right`` -- the right summand

        OUTPUT: the sum as a new element

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: df = (x, 1), (y, 1), (x*y + 1, 1)
            sage: f = FFPD(2, df)
            sage: g = FFPD(2*x*y, df)
            sage: f + g
            (2, [(y, 1), (x, 1)])
        """
        return FractionWithFactoredDenominatorSum([left, right]).sum()

    def _mul_(left, right):
        r"""
        Return the product of ``left`` with ``right``.

        INPUT:

        - ``left`` -- the left factor (i.e. ``self``)

        - ``right`` -- the right factor

        OUTPUT: the product as a new element

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: f = FFPD(2, [(x, 1), (x*y + 1, 1), (x*y^2 + 1, 1)])
            sage: g = FFPD(2*x*y, [(y, 1), (x*y + 1, 1), (x^2*y + 1, 1)])
            sage: f * g
            (4, [(x*y + 1, 1), (x*y + 1, 1), (x*y^2 + 1, 1), (x^2*y + 1, 1)])
        """
        numer = left.numerator() * right.numerator()
        df = left.denominator_factored() + right.denominator_factored()
        return left.parent()(numer, df)


class FractionWithFactoredDenominatorRing(UniqueRepresentation, Parent):
    r"""
    This is the ring of fractions with factored denominator.

    INPUT:

    - ``denominator_ring`` -- the base ring (a polynomial ring)

    - ``numerator_ring`` -- (optional) the numerator ring; the default is
      the ``denominator_ring``

    - ``category`` -- (default: :class:`Rings`) the category

    .. SEEALSO::

        :class:`FractionWithFactoredDenominator`,
        :mod:`~sage.rings.asymptotic.asymptotics_multivariate_generating_functions`

    EXAMPLES::

        sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
        sage: R.<x,y> = PolynomialRing(QQ)
        sage: FFPD = FractionWithFactoredDenominatorRing(R)
        sage: df = [x, 1], [y, 1], [x*y+1, 1]
        sage: f = FFPD(x, df)  # indirect doctest
        sage: f
        (1, [(y, 1), (x*y + 1, 1)])

    AUTHORS:

    - Daniel Krenn (2014-12-01)
    """

    @staticmethod
    def __classcall_private__(cls, denominator_ring, numerator_ring=None, category=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD1 = FractionWithFactoredDenominatorRing(R)
            sage: cat = Rings().Commutative()
            sage: FFPD2 = FractionWithFactoredDenominatorRing(R, R, cat)
            sage: FFPD1 is FFPD2
            True
        """
        if numerator_ring is None:
            numerator_ring = denominator_ring
        if not numerator_ring.has_coerce_map_from(denominator_ring):
            raise ValueError('numerator ring {} has no coercion map from the '
                             'denominator ring {}'.format(
                                 numerator_ring, denominator_ring))
        category = Rings().Commutative().or_subcategory(category)
        return super().__classcall__(cls, denominator_ring,
                                     numerator_ring, category)

    def __init__(self, denominator_ring, numerator_ring=None, category=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: P.<X, Y> = ZZ[]
            sage: FractionWithFactoredDenominatorRing(P)
            Ring of fractions with factored denominator
            over Multivariate Polynomial Ring in X, Y over Integer Ring
        """
        self._numerator_ring = numerator_ring
        self._denominator_ring = denominator_ring
        Parent.__init__(self, denominator_ring, category=category)

    def _repr_(self) -> str:
        r"""
        Return a representation.

        OUTPUT: string

        TESTS::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: P.<X, Y> = ZZ[]
            sage: FractionWithFactoredDenominatorRing(P)  # indirect doctest
            Ring of fractions with factored denominator
            over Multivariate Polynomial Ring in X, Y over Integer Ring
        """
        return ("Ring of fractions with factored denominator "
                "over {!r}".format(self.base()))

    def base_ring(self):
        r"""
        Return the base ring.

        OUTPUT: a ring

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: P.<X, Y> = ZZ[]
            sage: F = FractionWithFactoredDenominatorRing(P); F
            Ring of fractions with factored denominator
            over Multivariate Polynomial Ring in X, Y over Integer Ring
            sage: F.base_ring()
            Integer Ring
            sage: F.base()
            Multivariate Polynomial Ring in X, Y over Integer Ring
        """
        return self.base().base_ring()

    def _element_constructor_(self, *args, **kwargs):
        r"""
        Return an element of this ring.

        See :class:`FractionWithFactoredDenominator` for details.

        TESTS::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: df = [x, 1], [y, 1], [x*y+1, 1]
            sage: f = FFPD(x, df)  # indirect doctest
            sage: f
            (1, [(y, 1), (x*y + 1, 1)])
        """
        R = self.base()
        Q = R.fraction_field()

        # process keyword arguments
        reduce = kwargs.pop('reduce', None)

        if kwargs:
            raise ValueError('Unknown keyword arguments '
                             '%s given' % (kwargs,))

        # process arguments
        if len(args) > 2:
            raise ValueError('too many arguments given')

        elif not args:
            raise ValueError('No argument given. '
                             'We are in serious troubles...')

        # At this point we have one or two input arguments.

        x = args[0]
        try:
            P = x.parent()
        except AttributeError:
            P = None

        denominator_factored = None  # init
        reduce_default = True

        if len(args) == 2:
            numerator, denominator_factored = args
            if numerator is None:
                numerator = R(0)
            if denominator_factored is None:
                denominator_factored = []

            from sage.rings.semirings.non_negative_integer_semiring import NN
            try:
                denominator_factored = sorted(
                    (R(d[0]), NN(d[1])) for d in denominator_factored)
            except TypeError:
                raise TypeError('factored denominator is not well-formed '
                                'or of wrong type')

        # From now on we only have one input argument;
        # it's called x and has parent P.

        elif isinstance(P, FractionWithFactoredDenominatorRing):
            numerator = x._numerator
            denominator_factored = self._denominator_factored
            reduce_default = False

        elif P == SR:
            numerator = x.numerator()
            denominator = x.denominator()
            reduce_default = False

        elif x in R:
            numerator = R(x)
            denominator_factored = []

        elif x in Q:
            quotient = Q(x)
            numerator = quotient.numerator()
            denominator = quotient.denominator()
            reduce_default = False

        elif hasattr(x, 'numerator') and hasattr(x, 'denominator'):
            numerator = x.numerator()
            denominator = x.denominator()
            reduce_default = False

        else:
            raise TypeError('element {} is not contained in {}'.format(x, self))

        if reduce is None:
            reduce = reduce_default

        if denominator_factored is None:
            if denominator not in R:
                raise TypeError('extracted denominator {} is not in {}'.format(denominator, self))
            p = numerator
            q = R(denominator)

            from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
            from sage.rings.polynomial.multi_polynomial_ring_base import MPolynomialRing_base
            if isinstance(R, PolynomialRing_general) or isinstance(R, MPolynomialRing_base):
                if not R(q).is_unit():
                    # Factor denominator
                    try:
                        df = q.factor()
                    except NotImplementedError:
                        # Singular's factor() needs 'proof=False'.
                        df = q.factor(proof=False)
                    numerator = p / df.unit()
                    df = sorted(tuple(t) for t in df)  # sort for consistency
                    denominator_factored = df
            else:
                # At this point, denominator could not be factored.
                numerator = p / q
                denominator_factored = []

        return self.element_class(self,
                                  numerator=numerator,
                                  denominator_factored=denominator_factored,
                                  reduce=reduce)

    def _coerce_map_from_(self, P):
        r"""
        Check if there is a coercion from the given parent.

        INPUT:

        - ``P`` -- a parent

        OUTPUT: ``True`` if there is a coercion, otherwise ``False`` or ``None``

        TESTS::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: Q = QQ['x,y']
            sage: FFPD_QQ = FractionWithFactoredDenominatorRing(Q)
            sage: FFPD_QQ.has_coerce_map_from(Q)
            True
            sage: FFPD_QQ.has_coerce_map_from(QQ)
            True
            sage: FFPD_QQ.has_coerce_map_from(ZZ)
            True
            sage: FFPD_QQ.has_coerce_map_from(Q.fraction_field())
            True
            sage: Z = ZZ['x,y']
            sage: FFPD_ZZ = FractionWithFactoredDenominatorRing(Z)
            sage: FFPD_ZZ.has_coerce_map_from(FFPD_QQ)
            False
            sage: FFPD_QQ.has_coerce_map_from(FFPD_ZZ)
            True
            sage: FFPD_ZZ.has_coerce_map_from(QQ)
            False
            sage: FFPD_ZZ.has_coerce_map_from(Z.fraction_field())
            True
            sage: FFPD_ZZ.has_coerce_map_from(Q.fraction_field())
            False
            sage: FFPD_QQ.has_coerce_map_from(Z.fraction_field())
            True
        """
        if isinstance(P, FractionWithFactoredDenominatorRing):
            if self.base().has_coerce_map_from(P.base()):
                return True

        from sage.rings.fraction_field import FractionField_generic
        if isinstance(P, FractionField_generic):
            B = P.base()
            from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
            from sage.rings.polynomial.multi_polynomial_ring_base import MPolynomialRing_base
            if isinstance(B, PolynomialRing_general) or isinstance(B, MPolynomialRing_base):
                if self.base().has_coerce_map_from(B):
                    return True

        if self.base().has_coerce_map_from(P):
            return True

    def _an_element_(self):
        r"""
        Return an element.

        TESTS::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: FFPD.an_element()  # indirect doctest
            (42, [(x, 3)])
        """
        from sage.rings.semirings.non_negative_integer_semiring import NN
        return self(NN.an_element(), [(self.base().an_element(), NN(3))])

    Element = FractionWithFactoredDenominator


class FractionWithFactoredDenominatorSum(list):
    r"""
    A list representing the sum of :class:`FractionWithFactoredDenominator`
    objects with distinct denominator factorizations.

    AUTHORS:

    - Alexander Raichev (2012-06-25)

    - Daniel Krenn (2014-12-01)
    """

    def __repr__(self) -> str:
        r"""
        Return a string representation of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing, FractionWithFactoredDenominatorSum
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: f = FFPD(x + y, [(y, 1), (x, 1)])
            sage: g = FFPD(x**2 + y, [(y, 1), (x, 2)])
            sage: FractionWithFactoredDenominatorSum([f, g])
            (x + y, [(y, 1), (x, 1)]) + (x^2 + y, [(y, 1), (x, 2)])
        """
        return ' + '.join(repr(r) for r in self)

    def __eq__(self, other) -> bool:
        r"""
        Return ``True`` if ``self`` is equal to ``other``.

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing, FractionWithFactoredDenominatorSum
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: f = FFPD(x + y, [(y, 1), (x, 1)])
            sage: g = FFPD(x*(x + y), [(y, 1), (x, 2)])
            sage: s = FractionWithFactoredDenominatorSum([f]); s
            (x + y, [(y, 1), (x, 1)])
            sage: t = FractionWithFactoredDenominatorSum([g]); t
            (x + y, [(y, 1), (x, 1)])
            sage: s == t
            True
        """
        from operator import methodcaller
        return (sorted(self, key=methodcaller('_total_order_key_')) ==
                sorted(other, key=methodcaller('_total_order_key_')))

    def __ne__(self, other) -> bool:
        r"""
        Return ``True`` if ``self`` is not equal to ``other``.

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing, FractionWithFactoredDenominatorSum
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: f = FFPD(x + y, [(y, 1), (x, 1)])
            sage: g = FFPD(x + y, [(y, 1), (x, 2)])
            sage: s = FractionWithFactoredDenominatorSum([f]); s
            (x + y, [(y, 1), (x, 1)])
            sage: t = FractionWithFactoredDenominatorSum([g]); t
            (x + y, [(y, 1), (x, 2)])
            sage: s != t
            True
        """
        return not self.__eq__(other)

    @property
    def denominator_ring(self):
        r"""
        Return the polynomial ring of the denominators of ``self``.

        OUTPUT: a ring or ``None`` if the list is empty

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing, FractionWithFactoredDenominatorSum
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R)
            sage: f = FFPD(x + y, [(y, 1), (x, 1)])
            sage: s = FractionWithFactoredDenominatorSum([f])
            sage: s.denominator_ring
            Multivariate Polynomial Ring in x, y over Rational Field
            sage: g = FFPD(x + y, [])
            sage: t = FractionWithFactoredDenominatorSum([g])
            sage: t.denominator_ring
            Multivariate Polynomial Ring in x, y over Rational Field
        """
        for r in self:
            return r.denominator_ring
        return None

    def whole_and_parts(self):
        r"""
        Rewrite ``self`` as a sum of a (possibly zero) polynomial
        followed by reduced rational expressions.

        OUTPUT: an instance of :class:`FractionWithFactoredDenominatorSum`

        Only useful for multivariate decompositions.

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing, FractionWithFactoredDenominatorSum
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: f = x**2 + 3*y + 1/x + 1/y
            sage: f = FFPD(f); f
            (x^3*y + 3*x*y^2 + x + y, [(y, 1), (x, 1)])
            sage: FractionWithFactoredDenominatorSum([f]).whole_and_parts()
            (x^2 + 3*y, []) + (x + y, [(y, 1), (x, 1)])

            sage: f = cos(x)**2 + 3*y + 1/x + 1/y; f
            cos(x)^2 + 3*y + 1/x + 1/y
            sage: G = f.numerator()
            sage: H = R(f.denominator())
            sage: f = FFPD(G, H.factor()); f
            (x*y*cos(x)^2 + 3*x*y^2 + x + y, [(y, 1), (x, 1)])
            sage: FractionWithFactoredDenominatorSum([f]).whole_and_parts()
            (0, []) + (x*y*cos(x)^2 + 3*x*y^2 + x + y, [(y, 1), (x, 1)])
        """
        whole = 0
        parts = []
        R = self.denominator_ring
        for r in self:
            # Since r has already passed through FFPD.__init__()'s reducing
            # procedure, r is already in lowest terms.
            # Check if can write r as a mixed fraction: whole + fraction.
            p = r.numerator()
            q = r.denominator()
            if q == 1:
                # r is already whole
                whole += p
            else:
                try:
                    # Coerce p into R and divide p by q
                    p = R(p)
                    a, b = p.quo_rem(q)
                except TypeError:
                    # p is not in R and so can't divide p by q
                    a = 0
                    b = p
                whole += a
                parts.append(r.parent()(b, r.denominator_factored(), reduce=False))
        return FractionWithFactoredDenominatorSum(
            [r.parent()(whole, ())] + parts)  # r.parent() is not the nicest here

    def _combine_like_terms_(self):
        r"""
        Combine terms in ``self`` with the same denominator.
        Only useful for multivariate decompositions.

        OUTPUT: an instance of :class:`FractionWithFactoredDenominatorSum`

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing, FractionWithFactoredDenominatorSum
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: f = FFPD(1/(x * y * (x*y + 1)))
            sage: g = FFPD(x/(x * y * (x*y + 1)))
            sage: s = FractionWithFactoredDenominatorSum([f, g, f])
            sage: t = s._combine_like_terms_()
            sage: s
            (1, [(y, 1), (x, 1), (x*y + 1, 1)]) +
            (1, [(y, 1), (x*y + 1, 1)]) +
            (1, [(y, 1), (x, 1), (x*y + 1, 1)])
            sage: t
            (1, [(y, 1), (x*y + 1, 1)]) + (2, [(y, 1), (x, 1), (x*y + 1, 1)])

            sage: H = x * y * (x*y + 1)
            sage: f = FFPD(1, H.factor())
            sage: g = FFPD(exp(x + y), H.factor())
            sage: s = FractionWithFactoredDenominatorSum([f, g])
            sage: s
            (1, [(y, 1), (x, 1), (x*y + 1, 1)]) +
            (e^(x + y), [(y, 1), (x, 1), (x*y + 1, 1)])
            sage: t = s._combine_like_terms_()
            sage: t
            (e^(x + y) + 1, [(y, 1), (x, 1), (x*y + 1, 1)])
        """
        if not self:
            return self

        from operator import methodcaller
        # Combine like terms.
        FFPDs = sorted(self, key=methodcaller('_total_order_key_'))
        new_FFPDs = []
        temp = FFPDs[0]
        for f in FFPDs[1:]:
            if temp.denominator_factored() == f.denominator_factored():
                # Add f to temp.
                num = temp.numerator() + f.numerator()
                temp = f.parent()(num, temp.denominator_factored())
            else:
                # Append temp to new_FFPDs and update temp.
                new_FFPDs.append(temp)
                temp = f
        new_FFPDs.append(temp)
        return FractionWithFactoredDenominatorSum(new_FFPDs)

    def sum(self):
        r"""
        Return the sum of the elements in ``self``.

        OUTPUT: an instance of :class:`FractionWithFactoredDenominator`

        EXAMPLES::

            sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing, FractionWithFactoredDenominatorSum
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: FFPD = FractionWithFactoredDenominatorRing(R, SR)
            sage: df = (x, 1), (y, 1), (x*y + 1, 1)
            sage: f = FFPD(2, df)
            sage: g = FFPD(2*x*y, df)
            sage: FractionWithFactoredDenominatorSum([f, g])
            (2, [(y, 1), (x, 1), (x*y + 1, 1)]) + (2, [(x*y + 1, 1)])
            sage: FractionWithFactoredDenominatorSum([f, g]).sum()
            (2, [(y, 1), (x, 1)])

            sage: f = FFPD(cos(x), [(x, 2)])
            sage: g = FFPD(cos(y), [(x, 1), (y, 2)])
            sage: FractionWithFactoredDenominatorSum([f, g])
            (cos(x), [(x, 2)]) + (cos(y), [(y, 2), (x, 1)])
            sage: FractionWithFactoredDenominatorSum([f, g]).sum()
            (y^2*cos(x) + x*cos(y), [(y, 2), (x, 2)])
        """
        if not self:
            return self

        # Compute the sum's numerator and denominator.
        R = self.denominator_ring
        summy = sum(f.quotient() for f in self)
        numer = summy.numerator()
        denom = R(summy.denominator())

        # Compute the sum's denominator factorization.
        # Could use the factor() command, but it's probably faster to use
        # the irreducible factors of the denominators of self.
        df = []  # The denominator factorization for the sum.
        if denom == 1:
            # Done
            return FractionWithFactoredDenominatorRing(numer.parent())(numer, df, reduce=False)

        factors = []
        for f in self:
            factors.extend([q for (q, e) in f.denominator_factored()])

        # Eliminate repeats from factors and sort.
        factors = sorted(set(factors))

        # The irreducible factors of denom lie in factors.
        # Use this fact to build df.
        for q in factors:
            e = 0
            quo, rem = denom.quo_rem(q)
            while rem == 0:
                e += 1
                denom = quo
                quo, rem = denom.quo_rem(q)
            if e > 0:
                df.append((q, e))
        return FractionWithFactoredDenominatorRing(numer.parent())(numer, df, reduce=False)


#####################################################################
#  Helper functions


def diff_prod(f_derivs, u, g, X, interval, end, uderivs, atc):
    r"""
    Take various derivatives of the equation `f = ug`,
    evaluate them at a point `c`, and solve for the derivatives of `u`.

    INPUT:

    - ``f_derivs`` -- dictionary whose keys are all tuples of the form
      ``s + end``, where ``s`` is a sequence of variables from ``X`` whose
      length lies in ``interval``, and whose values are the derivatives
      of a function `f` evaluated at `c`
    - ``u`` -- a callable symbolic function
    - ``g`` -- an expression or callable symbolic function
    - ``X`` -- list of symbolic variables
    - ``interval`` -- list of positive integers
      Call the first and last values `n` and `nn`, respectively
    - ``end`` -- a possibly empty list of repetitions of the
      variable ``z``, where ``z`` is the last element of ``X``
    - ``uderivs`` -- dictionary whose keys are the symbolic
      derivatives of order 0 to order `n-1` of ``u`` evaluated at `c`
      and whose values are the corresponding derivatives evaluated at `c`
    - ``atc`` -- dictionary whose keys are the keys of `c` and all
      the symbolic derivatives of order 0 to order `nn` of ``g``
      evaluated `c` and whose values are the corresponding
      derivatives evaluated at `c`

    OUTPUT:

    A dictionary whose keys are the derivatives of ``u`` up to order
    `nn` and whose values are those derivatives evaluated at `c`.

    This function works by differentiating the equation `f = ug`
    with respect to the variable sequence ``s + end``,
    for all tuples ``s`` of ``X`` of lengths in ``interval``,
    evaluating at the point `c` ,
    and solving for the remaining derivatives of ``u``.
    This function assumes that ``u`` never appears in the
    differentiations of `f = ug` after evaluating at `c`.

    .. NOTE::

        For internal use by
        :meth:`FractionWithFactoredDenominator.asymptotics_multiple()`.

    EXAMPLES::

        sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import diff_prod
        sage: u = function('u')(x)
        sage: g = function('g')(x)
        sage: fd = {(x,):1,(x, x):1}
        sage: ud = {u(x=2): 1}
        sage: atc = {x: 2, g(x=2): 3, diff(g, x)(x=2): 5}
        sage: atc[diff(g, x, x)(x=2)] = 7
        sage: dd = diff_prod(fd, u, g, [x], [1, 2], [], ud, atc)
        sage: dd[diff(u, x, 2)(x=2)]
        22/9
    """
    from sage.symbolic.relation import solve

    for l in interval:
        D = {}
        rhs = []
        lhs = []
        new_vars = []
        for t in combinations_with_replacement(X, l):
            t = list(t)
            s = t + end
            lhs.append(f_derivs[tuple(s)])
            rhs.append(diff(u * g, s).subs(atc).subs(uderivs))
            # Since Sage's solve command can't take derivatives as variable
            # names, make new variables based on t to stand in for
            # diff(u, t) and store them in D.
            new_var = SR.temp_var()
            new_vars.append(new_var)
            D[diff(u, t).subs(atc)] = new_var
        eqns = [lhs[i] == rhs[i].subs(uderivs).subs(D)
                for i in range(len(lhs))]
        variables = D.values()
        sol = solve(eqns, *variables, solution_dict=True)
        uderivs.update(subs_all(D, sol[ZZ.zero()]))
        SR.cleanup_var(new_vars)
    return uderivs


def permutation_sign(s, u):
    r"""
    This function returns the sign of the permutation on
    ``1, ..., len(u)`` that is induced by the sublist ``s`` of ``u``.

    .. NOTE::

        This function was intended for internal use and is deprecated now
        (:issue:`29465`).

    INPUT:

    - ``s`` -- a sublist of ``u``
    - ``u`` -- list

    OUTPUT:

    The sign of the permutation obtained by taking indices
    within ``u`` of the list ``s + sc``, where ``sc`` is ``u``
    with the elements of ``s`` removed.

    EXAMPLES::

        sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import permutation_sign
        sage: u = ['a', 'b', 'c', 'd', 'e']
        sage: s = ['b', 'd']
        sage: permutation_sign(s, u)
        doctest:...: DeprecationWarning: the function permutation_sign is deprecated
        See https://github.com/sagemath/sage/issues/29465 for details.
        -1
        sage: s = ['d', 'b']
        sage: permutation_sign(s, u)
        1
    """
    from sage.misc.superseded import deprecation
    deprecation(29465, 'the function permutation_sign is deprecated')
    from sage.combinat.permutation import Permutation

    # Convert lists to lists of numbers in {1,..., len(u)}
    A = [i + 1 for i in range(len(u))]
    B = [u.index(x) + 1 for x in s]

    C = sorted(set(A).difference(set(B)))
    P = Permutation(B + C)
    return P.signature()


def subs_all(f, sub, simplify=False):
    r"""
    Return the items of `f` substituted by the dictionaries
    of ``sub`` in order of their appearance in ``sub``.

    INPUT:

    - ``f`` -- an individual or list of symbolic expressions
      or dictionaries
    - ``sub`` -- an individual or list of dictionaries
    - ``simplify`` -- boolean (default: ``False``); set to ``True`` to
      simplify the result

    OUTPUT:

    The items of ``f`` substituted by the dictionaries of ``sub`` in order
    of their appearance in ``sub``. The ``subs()`` command is used. If
    simplify is ``True``, then ``simplify()`` is used after substitution.

    EXAMPLES::

        sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import subs_all
        sage: var('x, y, z')
        (x, y, z)
        sage: a = {x:1}
        sage: b = {y:2}
        sage: c = {z:3}
        sage: subs_all(x + y + z, a)
        y + z + 1
        sage: subs_all(x + y + z, [c, a])
        y + 4
        sage: subs_all([x + y + z, y^2], b)
        [x + z + 2, 4]
        sage: subs_all([x + y + z, y^2], [b, c])
        [x + 5, 4]

    ::

        sage: var('x, y')
        (x, y)
        sage: a = {'foo': x**2 + y**2, 'bar': x - y}
        sage: b = {x: 1, y: 2}
        sage: subs_all(a, b)
        {'bar': -1, 'foo': 5}
    """
    singleton = False
    if not isinstance(f, (list, tuple)):
        f = [f]
        singleton = True
    if not isinstance(sub, (list, tuple)):
        sub = [sub]
    g = []
    for ff in f:
        for D in sub:
            if isinstance(ff, dict):
                ff = {k: ff[k].subs(D) for k in ff}
            else:
                ff = ff.subs(D)
        g.append(ff)

    if singleton and simplify:
        if isinstance(g[0], dict):
            return g[0]
        return g[0].simplify()

    if singleton and not simplify:
        return g[0]

    if not singleton and simplify:
        G = []
        for gg in g:
            if isinstance(gg, dict):
                G.append(gg)
            else:
                G.append(gg.simplify())
        return G

    return g


def diff_all(f, V, n, ending=[], sub=None, sub_final=None,
             zero_order=0, rekey=None):
    r"""
    Return a dictionary of representative mixed partial
    derivatives of `f` from order 1 up to order `n` with respect to the
    variables in `V`.

    The default is to key the dictionary by all nondecreasing sequences
    in `V` of length 1 up to length `n`.

    INPUT:

    - ``f`` -- an individual or list of `\mathcal{C}^{n+1}` functions
    - ``V`` -- list of variables occurring in `f`
    - ``n`` -- a natural number
    - ``ending`` -- list of variables in `V`
    - ``sub`` -- an individual or list of dictionaries
    - ``sub_final`` -- an individual or list of dictionaries
    - ``rekey`` -- a callable symbolic function in `V` or list thereof
    - ``zero_order`` -- a natural number

    OUTPUT: the dictionary ``{s_1:deriv_1, ..., sr:deriv_r}``

    Here ``s_1, ..., s_r`` is a listing of
    all nondecreasing sequences of length 1 up to length `n` over the
    alphabet `V`, where `w > v` in `X` if and only if ``str(w) > str(v)``,
    and ``deriv_j`` is the derivative of `f` with respect to the derivative
    sequence ``s_j`` and simplified with respect to the substitutions in
    ``sub`` and evaluated at ``sub_final``.
    Moreover, all derivatives with respect to sequences of length less than
    ``zero_order`` (derivatives of order less than ``zero_order`` )
    will be made zero.

    If ``rekey`` is nonempty, then ``s_1, ..., s_r`` will be replaced
    by the symbolic derivatives of the functions in ``rekey``.

    If ``ending`` is nonempty, then every derivative sequence ``s_j``
    will be suffixed by ``ending``.

    EXAMPLES::

        sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import diff_all
        sage: f = function('f')(x)
        sage: dd = diff_all(f, [x], 3)
        sage: dd[(x, x, x)]
        diff(f(x), x, x, x)

        sage: d1 = {diff(f, x): 4*x^3}
        sage: dd = diff_all(f, [x], 3, sub=d1)
        sage: dd[(x, x, x)]
        24*x

        sage: dd = diff_all(f, [x], 3, sub=d1, rekey=f)
        sage: dd[diff(f, x, 3)]
        24*x

        sage: a = {x:1}
        sage: dd = diff_all(f, [x], 3, sub=d1, rekey=f, sub_final=a)
        sage: dd[diff(f, x, 3)]
        24

    ::

        sage: X = var('x, y, z')
        sage: f = function('f')(*X)
        sage: dd = diff_all(f, X, 2, ending=[y, y, y])
        sage: dd[(z, y, y, y)]
        diff(f(x, y, z), y, y, y, z)

    ::

        sage: g = function('g')(*X)
        sage: dd = diff_all([f, g], X, 2)
        sage: dd[(0, y, z)]
        diff(f(x, y, z), y, z)

        sage: dd[(1, z, z)]
        diff(g(x, y, z), z, z)

        sage: f = exp(x*y*z)
        sage: ff = function('ff')(*X)
        sage: dd = diff_all(f, X, 2, rekey=ff)
        sage: dd[diff(ff, x, z)]
        x*y^2*z*e^(x*y*z) + y*e^(x*y*z)
    """
    singleton = False
    if not isinstance(f, list):
        f = [f]
        singleton = True

    # Build the dictionary of derivatives iteratively from a list
    # of nondecreasing derivative-order sequences.
    derivs = {}
    r = len(f)
    if ending:
        seeds = [ending]
        start = ZZ.one()
    else:
        seeds = [[v] for v in V]
        start = Integer(2)
    if singleton:
        for s in seeds:
            derivs[tuple(s)] = subs_all(diff(f[0], s), sub)
        for l in range(start, n + 1):
            for t in combinations_with_replacement(V, l):
                s = t + tuple(ending)
                derivs[s] = subs_all(diff(derivs[s[1:]], s[0]), sub)
    else:
        # Make the dictionary keys of the form (j, sequence of variables),
        # where j in range(r).
        for s in seeds:
            value = subs_all([diff(f[j], s) for j in range(r)], sub)
            derivs.update({tuple([j] + s): value[j] for j in range(r)})
        for l in range(start, n + 1):
            for t in combinations_with_replacement(V, l):
                s = t + tuple(ending)
                value = subs_all([diff(derivs[(j,) + s[1:]], s[0]) for j in range(r)], sub)
                derivs.update({(j,) + s: value[j] for j in range(r)})
    if zero_order:
        # Zero out all the derivatives of order < zero_order
        if singleton:
            for k in derivs:
                if len(k) < zero_order:
                    derivs[k] = ZZ.zero()
        else:
            # Ignore the first of element of k, which is an index.
            for k in derivs:
                if len(k) - 1 < zero_order:
                    derivs[k] = ZZ.zero()
    if sub_final:
        # Substitute sub_final into the values of derivs.
        for k in derivs:
            derivs[k] = subs_all(derivs[k], sub_final)
    if rekey:
        # Rekey the derivs dictionary by the value of rekey.
        F = rekey
        if singleton:
            # F must be a singleton.
            derivs = {diff(F, list(k)): derivs[k] for k in derivs}
        else:
            # F must be a list.
            derivs = {diff(F[k[0]], list(k)[1:]): derivs[k] for k in derivs}
    return derivs


def diff_op(A, B, AB_derivs, V, M, r, N):
    r"""
    Return the derivatives `DD^{(l+k)}(A[j] B^l)` evaluated at a point
    `p` for various natural numbers `j, k, l` which depend on `r` and `N`.

    Here `DD` is a specific second-order linear differential operator
    that depends on `M` , `A` is a list of symbolic functions,
    `B` is symbolic function, and ``AB_derivs`` contains all the derivatives
    of `A` and `B` evaluated at `p` that are necessary for the computation.

    INPUT:

    - ``A`` -- a single or length ``r`` list of symbolic functions in the
      variables ``V``
    - ``B`` -- a symbolic function in the variables ``V``
    - ``AB_derivs`` -- dictionary whose keys are the (symbolic)
      derivatives of ``A[0], ..., A[r-1]`` up to order ``2 * N-2`` and
      the (symbolic) derivatives of ``B`` up to order ``2 * N``;
      the values of the dictionary are complex numbers that are
      the keys evaluated at a common point `p`
    - ``V`` -- the variables of the ``A[j]`` and ``B``
    - ``M`` -- a symmetric `l \times l` matrix, where `l` is the
      length of ``V``
    - ``r``, ``N`` -- natural numbers

    OUTPUT: a dictionary

    The output is
    a dictionary whose keys are natural number tuples of the form
    `(j, k, l)`, where `l \leq 2k`, `j \leq r-1`, and `j+k \leq N-1`,
    and whose values are `DD^(l+k)(A[j] B^l)` evaluated at a point
    `p`, where `DD` is the linear second-order differential operator
    `-\sum_{i=0}^{l-1} \sum_{j=0}^{l-1} M[i][j]
    \partial^2 /(\partial V[j] \partial V[i])`.

    .. NOTE::

        For internal use by
        :meth:`FractionWithFactoredDenominator.asymptotics_smooth()` and
        :meth:`FractionWithFactoredDenominator.asymptotics_multiple()`.

    EXAMPLES::

        sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import diff_op
        sage: T = var('x, y')
        sage: A = function('A')(*tuple(T))
        sage: B = function('B')(*tuple(T))
        sage: AB_derivs = {}
        sage: M = matrix([[1, 2],[2, 1]])
        sage: DD = diff_op(A, B, AB_derivs, T, M, 1, 2)  # long time (see :issue:`35207`)
        sage: sorted(DD)                                 # long time
        [(0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 1, 2)]
        sage: DD[(0, 1, 2)].number_of_operands()         # long time
        246
    """
    from itertools import product
    from sage.misc.mrange import xmrange

    if not isinstance(A, list):
        A = [A]

    # First, compute the necessary product derivatives of A and B.
    product_derivs = {}
    for j, k in xmrange([r, N], tuple):
        if j + k < N:
            for l in range(2 * k + 1):
                for s in combinations_with_replacement(V, 2 * (k + l)):
                    DF = diff(A[j] * B ** l, list(s)).subs(AB_derivs)
                    product_derivs[(j, k, l) + s] = DF

    # Second, compute DD^(k+l)(A[j]*B^l)(p) and store values in dictionary.
    DD = {}
    rows = M.nrows()
    for j, k in xmrange([r, N], tuple):
        if j + k < N:
            for l in range(2 * k + 1):
                # Take advantage of the symmetry of M by ignoring
                # the upper-diagonal entries of M and multiplying by
                # appropriate powers of 2.
                if k + l == 0:
                    DD[(j, k, l)] = product_derivs[(j, k, l)]
                    continue
                S = [(a, b) for a, b in xmrange([rows, rows], tuple) if b <= a]
                P = product(S, repeat=k + l)
                diffo = ZZ.zero()
                for t in P:
                    idx = (j, k, l) + diff_seq(V, t)
                    if product_derivs[idx] != ZZ.zero():
                        MM = ZZ.one()
                        for (a, b) in t:
                            MM *= M[a][b]
                            if a != b:
                                MM *= Integer(2)
                        diffo += MM * product_derivs[idx]
                DD[(j, k, l)] = (-ZZ.one()) ** (k + l) * diffo
    return DD


def diff_seq(V, s):
    r"""
    Given a list ``s`` of tuples of natural numbers, return the
    list of elements of ``V`` with indices the elements of the elements
    of ``s``.

    INPUT:

    - ``V`` -- list
    - ``s`` -- list of tuples of natural numbers in the interval
      ``range(len(V))``

    OUTPUT:

    The tuple ``tuple([V[tt] for tt in sorted(t)])``, where ``t`` is the
    list of elements of the elements of ``s``.

    EXAMPLES::

        sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import diff_seq
        sage: V = list(var('x, t, z'))
        sage: diff_seq(V,([0, 1],[0, 2, 1],[0, 0]))
        (x, x, x, x, t, t, z)

    .. NOTE::

        This function is for internal use by :func:`diff_op()`.
    """
    t = []
    for ss in s:
        t.extend(ss)
    return tuple([V[tt] for tt in sorted(t)])


def diff_op_simple(A, B, AB_derivs, x, v, a, N):
    r"""
    Return `DD^(e k + v l)(A B^l)` evaluated at a point `p` for
    various natural numbers `e, k, l` that depend on `v` and `N`.

    Here `DD` is a specific linear differential operator that depends
    on `a` and `v` , `A` and `B` are symbolic functions, and ``AB_derivs``
    contains all the derivatives of `A` and `B` evaluated at `p` that are
    necessary for the computation.

    .. NOTE::

        For internal use by the function
        :meth:`FractionWithFactoredDenominator.asymptotics_smooth()`.

    INPUT:

    - ``A``, ``B`` -- symbolic functions in the variable ``x``
    - ``AB_derivs`` -- dictionary whose keys are the (symbolic)
      derivatives of ``A`` up to order ``2 * N`` if ``v`` is even or
      ``N`` if ``v`` is odd and the (symbolic) derivatives of ``B``
      up to order ``2 * N + v`` if ``v`` is even or ``N + v``
      if ``v`` is odd; the values of the dictionary are complex numbers
      that are the keys evaluated at a common point `p`
    - ``x`` -- a symbolic variable
    - ``a`` -- a complex number
    - ``v``, ``N`` -- natural numbers

    OUTPUT: a dictionary

    The output is
    a dictionary whose keys are natural number pairs of the form `(k, l)`,
    where `k < N` and `l \leq 2k` and whose values are
    `DD^(e k + v l)(A B^l)` evaluated at a point `p`.
    Here `e=2` if `v` is even, `e=1` if `v` is odd, and `DD` is the
    linear differential operator
    `(a^{-1/v} d/dt)` if `v` is even and
    `(|a|^{-1/v} i \text{sgn}(a) d/dt)` if `v` is odd.

    EXAMPLES::

        sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import diff_op_simple
        sage: A = function('A')(x)
        sage: B = function('B')(x)
        sage: AB_derivs = {}
        sage: sorted(diff_op_simple(A, B, AB_derivs, x, 3, 2, 2).items())
        [((0, 0), A(x)),
         ((1, 0), 1/2*I*2^(2/3)*diff(A(x), x)),
         ((1, 1),
          1/4*2^(2/3)*(B(x)*diff(A(x), x, x, x, x) + 4*diff(A(x), x, x, x)*diff(B(x), x) + 6*diff(A(x), x, x)*diff(B(x), x, x) + 4*diff(A(x), x)*diff(B(x), x, x, x) + A(x)*diff(B(x), x, x, x, x)))]
    """
    from sage.misc.functional import sqrt

    I = sqrt(-ZZ.one())
    DD = {}
    if v.mod(Integer(2)) == ZZ.zero():
        for k in range(N):
            for l in range(2 * k + 1):
                DD[(k, l)] = ((a ** (-ZZ.one() / v)) ** (2 * k + v * l) *
                              diff(A * B ** l, x,
                                   2 * k + v * l).subs(AB_derivs))
    else:
        for k in range(N):
            for l in range(k + 1):
                DD[(k, l)] = ((abs(a) ** (-ZZ.one() / v) * I *
                               a / abs(a)) ** (k + v * l) *
                              diff(A * B ** l, x,
                                   k + v * l).subs(AB_derivs))
    return DD


def direction(v, coordinate=None):
    r"""
    Return ``[vv/v[coordinate] for vv in v]`` where
    ``coordinate`` is the last index of ``v`` if not specified otherwise.

    INPUT:

    - ``v`` -- a vector
    - ``coordinate`` -- (default: ``None``) an index for ``v``

    EXAMPLES::

        sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import direction
        sage: direction([2, 3, 5])
        (2/5, 3/5, 1)
        sage: direction([2, 3, 5], 0)
        (1, 3/2, 5/2)
    """
    if coordinate is None:
        coordinate = len(v) - 1
    return tuple([vv / v[coordinate] for vv in v])


def coerce_point(R, p):
    r"""
    Coerce the keys of the dictionary ``p`` into the ring ``R``.

    .. WARNING::

        This method assumes that it is possible.

    EXAMPLES::

        sage: from sage.rings.asymptotic.asymptotics_multivariate_generating_functions import FractionWithFactoredDenominatorRing, coerce_point
        sage: R.<x,y> = PolynomialRing(QQ)
        sage: FFPD = FractionWithFactoredDenominatorRing(R)
        sage: f = FFPD()
        sage: p = {SR(x): 1, SR(y): 7/8}
        sage: for k in sorted(p, key=str):
        ....:     print("{} {} {}".format(k, k.parent(), p[k]))
        x Symbolic Ring 1
        y Symbolic Ring 7/8
        sage: q = coerce_point(R, p)
        sage: for k in sorted(q, key=str):
        ....:     print("{} {} {}".format(k, k.parent(), q[k]))
        x Multivariate Polynomial Ring in x, y over Rational Field 1
        y Multivariate Polynomial Ring in x, y over Rational Field 7/8
    """
    if p and list(p)[0].parent() != R:
        try:
            return {x: p[SR(x)] for x in R.gens()}
        except TypeError:
            pass
    return p
