r"""
Quotient fields
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_singleton import Category_singleton
from sage.misc.abstract_method import abstract_method
from sage.categories.fields import Fields

from sage.structure.element import coerce_binop

class QuotientFields(Category_singleton):
    """
    The category of quotient fields over an integral domain

    EXAMPLES::

        sage: QuotientFields()
        Category of quotient fields
        sage: QuotientFields().super_categories()
        [Category of fields]

    TESTS::

        sage: TestSuite(QuotientFields()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: QuotientFields().super_categories()
            [Category of fields]
        """
        return [Fields()]

    class ParentMethods:
        pass

    class ElementMethods:

        @abstract_method
        def numerator(self):
            pass

        @abstract_method
        def denominator(self):
            pass

        @coerce_binop
        def gcd(self, other):
            """
            Greatest common divisor

            .. NOTE::

                In a field, the greatest common divisor is not very informative,
                as it is only determined up to a unit. But in the fraction field
                of an integral domain that provides both gcd and lcm, it is
                possible to be a bit more specific and define the gcd uniquely
                up to a unit of the base ring (rather than in the fraction
                field).

            AUTHOR:

            - Simon King (2011-02): See :trac:`10771`

            EXAMPLES::

                sage: R.<x> = QQ['x']
                sage: p = (1+x)^3*(1+2*x^2)/(1-x^5)
                sage: q = (1+x)^2*(1+3*x^2)/(1-x^4)
                sage: factor(p)
                (-2) * (x - 1)^-1 * (x + 1)^3 * (x^2 + 1/2) * (x^4 + x^3 + x^2 + x + 1)^-1
                sage: factor(q)
                (-3) * (x - 1)^-1 * (x + 1) * (x^2 + 1)^-1 * (x^2 + 1/3)
                sage: gcd(p,q)
                (x + 1)/(x^7 + x^5 - x^2 - 1)
                sage: factor(gcd(p,q))
                (x - 1)^-1 * (x + 1) * (x^2 + 1)^-1 * (x^4 + x^3 + x^2 + x + 1)^-1
                sage: factor(gcd(p,1+x))
                (x - 1)^-1 * (x + 1) * (x^4 + x^3 + x^2 + x + 1)^-1
                sage: factor(gcd(1+x,q))
                (x - 1)^-1 * (x + 1) * (x^2 + 1)^-1

            TESTS:

            The following tests that the fraction field returns a correct gcd
            even if the base ring does not provide lcm and gcd::

                sage: R = ZZ.extension(x^2+1, names='i')
                sage: i = R.1
                sage: gcd(5, 3 + 4*i)
                -i - 2
                sage: P.<t> = R[]
                sage: gcd(t, i)
                Traceback (most recent call last):
                ...
                NotImplementedError: Gaussian Integers in Number Field in i with defining polynomial x^2 + 1 does not provide a gcd implementation for univariate polynomials
                sage: q = t/(t+1); q.parent()
                Fraction Field of Univariate Polynomial Ring in t over Gaussian Integers in Number Field in i with defining polynomial x^2 + 1
                sage: gcd(q, q)
                1
                sage: q.gcd(0)
                1
                sage: (q*0).gcd(0)
                0
            """
            P = self.parent()
            try:
                selfN = self.numerator()
                selfD = self.denominator()
                selfGCD = selfN.gcd(selfD)
                otherN = other.numerator()
                otherD = other.denominator()
                otherGCD = otherN.gcd(otherD)
                selfN = selfN // selfGCD
                selfD = selfD // selfGCD
                otherN = otherN // otherGCD
                otherD = otherD // otherGCD
                tmp = P(selfN.gcd(otherN))/P(selfD.lcm(otherD))
                return tmp
            except (AttributeError, NotImplementedError, TypeError, ValueError):
                zero = P.zero()
                if self == zero and other == zero:
                    return zero
                return P.one()

        @coerce_binop
        def lcm(self, other):
            """
            Least common multiple

            In a field, the least common multiple is not very informative, as it
            is only determined up to a unit. But in the fraction field of an
            integral domain that provides both gcd and lcm, it is reasonable to
            be a bit more specific and to define the least common multiple so
            that it restricts to the usual least common multiple in the base
            ring and is unique up to a unit of the base ring (rather than up to
            a unit of the fraction field).

            The least common multiple is easily described in terms of the
            prime decomposition. A rational number can be written as a product
            of primes with integer (positive or negative) powers in a unique
            way. The least common multiple of two rational numbers `x` and `y`
            can then be defined by specifying that the exponent of every prime
            `p` in `lcm(x,y)` is the supremum of the exponents of `p` in `x`,
            and the exponent of `p` in `y` (where the primes that does not
            appear in the decomposition of `x` or `y` are considered to have
            exponent zero).


            AUTHOR:

            - Simon King (2011-02): See :trac:`10771`

            EXAMPLES::

                sage: lcm(2/3, 1/5)
                2

            Indeed `2/3 = 2^1 3^{-1} 5^0` and `1/5 = 2^0 3^0
            5^{-1}`, so `lcm(2/3,1/5)= 2^1 3^0 5^0 = 2`.

                sage: lcm(1/3, 1/5)
                1
                sage: lcm(1/3, 1/6)
                1/3

            Some more involved examples::

                sage: R.<x> = QQ[]
                sage: p = (1+x)^3*(1+2*x^2)/(1-x^5)
                sage: q = (1+x)^2*(1+3*x^2)/(1-x^4)
                sage: factor(p)
                (-2) * (x - 1)^-1 * (x + 1)^3 * (x^2 + 1/2) * (x^4 + x^3 + x^2 + x + 1)^-1
                sage: factor(q)
                (-3) * (x - 1)^-1 * (x + 1) * (x^2 + 1)^-1 * (x^2 + 1/3)
                sage: factor(lcm(p,q))
                (x - 1)^-1 * (x + 1)^3 * (x^2 + 1/3) * (x^2 + 1/2)
                sage: factor(lcm(p,1+x))
                (x + 1)^3 * (x^2 + 1/2)
                sage: factor(lcm(1+x,q))
                (x + 1) * (x^2 + 1/3)

            TESTS:

            The following tests that the fraction field returns a correct lcm
            even if the base ring does not provide lcm and gcd::

                sage: R = ZZ.extension(x^2+1, names='i')
                sage: i = R.1
                sage: P.<t> = R[]
                sage: lcm(t, i)
                Traceback (most recent call last):
                ...
                NotImplementedError: Gaussian Integers in Number Field in i with defining polynomial x^2 + 1 does not provide a gcd implementation for univariate polynomials
                sage: q = t/(t+1); q.parent()
                Fraction Field of Univariate Polynomial Ring in t over Gaussian Integers in Number Field in i with defining polynomial x^2 + 1
                sage: lcm(q, q)
                1
                sage: q.lcm(0)
                0
                sage: (q*0).lcm(0)
                0

            Check that it is possible to take lcm of a rational and an integer
            (:trac:`17852`)::

                sage: (1/2).lcm(2)
                2
                sage: type((1/2).lcm(2))
                <class 'sage.rings.rational.Rational'>
            """
            P = self.parent()
            try:
                selfN = self.numerator()
                selfD = self.denominator()
                selfGCD = selfN.gcd(selfD)
                otherN = other.numerator()
                otherD = other.denominator()
                otherGCD = otherN.gcd(otherD)
                selfN = selfN // selfGCD
                selfD = selfD // selfGCD
                otherN = otherN // otherGCD
                otherD = otherD // otherGCD
                return P(selfN.lcm(otherN))/P(selfD.gcd(otherD))
            except (AttributeError, NotImplementedError, TypeError, ValueError):
                zero = P.zero()
                if self == zero or other == zero:
                    return zero
                return P.one()

        @coerce_binop
        def xgcd(self, other):
            """
            Return a triple ``(g,s,t)`` of elements of that field such that
            ``g`` is the greatest common divisor of ``self`` and ``other`` and
            ``g = s*self + t*other``.

            .. NOTE::

                In a field, the greatest common divisor is not very informative,
                as it is only determined up to a unit. But in the fraction field
                of an integral domain that provides both xgcd and lcm, it is
                possible to be a bit more specific and define the gcd uniquely
                up to a unit of the base ring (rather than in the fraction
                field).

            EXAMPLES::

                sage: QQ(3).xgcd(QQ(2))
                (1, 1, -1)
                sage: QQ(3).xgcd(QQ(1/2))
                (1/2, 0, 1)
                sage: QQ(1/3).xgcd(QQ(2))
                (1/3, 1, 0)
                sage: QQ(3/2).xgcd(QQ(5/2))
                (1/2, 2, -1)

                sage: R.<x> = QQ['x']
                sage: p = (1+x)^3*(1+2*x^2)/(1-x^5)
                sage: q = (1+x)^2*(1+3*x^2)/(1-x^4)
                sage: factor(p)
                (-2) * (x - 1)^-1 * (x + 1)^3 * (x^2 + 1/2) * (x^4 + x^3 + x^2 + x + 1)^-1
                sage: factor(q)
                (-3) * (x - 1)^-1 * (x + 1) * (x^2 + 1)^-1 * (x^2 + 1/3)
                sage: g,s,t = xgcd(p,q)
                sage: g
                (x + 1)/(x^7 + x^5 - x^2 - 1)
                sage: g == s*p + t*q
                True

            An example without a well defined gcd or xgcd on its base ring::

                sage: K = QuadraticField(5)
                sage: O = K.maximal_order()
                sage: R = PolynomialRing(O, 'x')
                sage: F = R.fraction_field()
                sage: x = F.gen(0)
                sage: x.gcd(x+1)
                1
                sage: x.xgcd(x+1)
                (1, 1/x, 0)
                sage: zero = F.zero()
                sage: zero.gcd(x)
                1
                sage: zero.xgcd(x)
                (1, 0, 1/x)
                sage: zero.xgcd(zero)
                (0, 0, 0)
            """
            P = self.parent()
            try:
                selfN = self.numerator()
                selfD = self.denominator()
                selfGCD = selfN.gcd(selfD)

                otherN = other.numerator()
                otherD = other.denominator()
                otherGCD = otherN.gcd(otherD)

                selfN = selfN // selfGCD
                selfD = selfD // selfGCD
                otherN = otherN // otherGCD
                otherD = otherD // otherGCD

                lcmD = selfD.lcm(otherD)
                g,s,t = selfN.xgcd(otherN)
                return (P(g)/P(lcmD), P(s*selfD)/P(lcmD),P(t*otherD)/P(lcmD))
            except (AttributeError, NotImplementedError, TypeError, ValueError):
                zero = self.parent().zero()
                one  = self.parent().one()
                if self != zero:
                    return (one, ~self, zero)
                elif other != zero:
                    return (one, zero, ~other)
                else:
                    return (zero, zero, zero)

        def factor(self, *args, **kwds):
            """
            Return the factorization of ``self`` over the base ring.

            INPUT:

            - ``*args`` - Arbitrary arguments suitable over the base ring
            - ``**kwds`` - Arbitrary keyword arguments suitable over the base ring

            OUTPUT:

            - Factorization of ``self`` over the base ring

            EXAMPLES::

                sage: K.<x> = QQ[]
                sage: f = (x^3+x)/(x-3)
                sage: f.factor()
                (x - 3)^-1 * x * (x^2 + 1)

            Here is an example to show that :trac:`7868` has been resolved::

                sage: R.<x,y> = GF(2)[]
                sage: f = x*y/(x+y)
                sage: f.factor()
                (x + y)^-1 * y * x
            """
            return (self.numerator().factor(*args, **kwds) /
                    self.denominator().factor(*args, **kwds))

        def partial_fraction_decomposition(self, decompose_powers=True):
            """
            Decomposes fraction field element into a whole part and a list of
            fraction field elements over prime power denominators.

            The sum will be equal to the original fraction.

            INPUT:

            - decompose_powers -- whether to decompose prime power
                                 denominators as opposed to having a single
                                 term for each irreducible factor of the
                                 denominator (default: True)

            OUTPUT:

            - Partial fraction decomposition of self over the base ring.

            AUTHORS:

            - Robert Bradshaw (2007-05-31)

            EXAMPLES::

                sage: S.<t> = QQ[]
                sage: q = 1/(t+1) + 2/(t+2) + 3/(t-3); q
                (6*t^2 + 4*t - 6)/(t^3 - 7*t - 6)
                sage: whole, parts = q.partial_fraction_decomposition(); parts
                [3/(t - 3), 1/(t + 1), 2/(t + 2)]
                sage: sum(parts) == q
                True
                sage: q = 1/(t^3+1) + 2/(t^2+2) + 3/(t-3)^5
                sage: whole, parts = q.partial_fraction_decomposition(); parts
                [1/3/(t + 1), 3/(t^5 - 15*t^4 + 90*t^3 - 270*t^2 + 405*t - 243), (-1/3*t + 2/3)/(t^2 - t + 1), 2/(t^2 + 2)]
                sage: sum(parts) == q
                True
                sage: q = 2*t / (t + 3)^2
                sage: q.partial_fraction_decomposition()
                (0, [2/(t + 3), -6/(t^2 + 6*t + 9)])
                sage: for p in q.partial_fraction_decomposition()[1]: print(p.factor())
                (2) * (t + 3)^-1
                (-6) * (t + 3)^-2
                sage: q.partial_fraction_decomposition(decompose_powers=False)
                (0, [2*t/(t^2 + 6*t + 9)])

            We can decompose over a given algebraic extension::

                sage: R.<x> = QQ[sqrt(2)][]
                sage: r =  1/(x^4+1)
                sage: r.partial_fraction_decomposition()
                (0,
                 [(-1/4*sqrt2*x + 1/2)/(x^2 - sqrt2*x + 1),
                  (1/4*sqrt2*x + 1/2)/(x^2 + sqrt2*x + 1)])

                sage: R.<x> = QQ[I][]  # of QQ[sqrt(-1)]
                sage: r =  1/(x^4+1)
                sage: r.partial_fraction_decomposition()
                (0, [(-1/2*I)/(x^2 - I), 1/2*I/(x^2 + I)])

            We can also ask Sage to find the least extension where the
            denominator factors in linear terms::

                sage: R.<x> = QQ[]
                sage: r = 1/(x^4+2)
                sage: N = r.denominator().splitting_field('a')
                sage: N
                Number Field in a with defining polynomial x^8 - 8*x^6 + 28*x^4 + 16*x^2 + 36
                sage: R1.<x1>=N[]
                sage: r1 = 1/(x1^4+2)
                sage: r1.partial_fraction_decomposition()
                (0,
                 [(-1/224*a^6 + 13/448*a^4 - 5/56*a^2 - 25/224)/(x1 - 1/28*a^6 + 13/56*a^4 - 5/7*a^2 - 25/28),
                  (1/224*a^6 - 13/448*a^4 + 5/56*a^2 + 25/224)/(x1 + 1/28*a^6 - 13/56*a^4 + 5/7*a^2 + 25/28),
                  (-5/1344*a^7 + 43/1344*a^5 - 85/672*a^3 - 31/672*a)/(x1 - 5/168*a^7 + 43/168*a^5 - 85/84*a^3 - 31/84*a),
                  (5/1344*a^7 - 43/1344*a^5 + 85/672*a^3 + 31/672*a)/(x1 + 5/168*a^7 - 43/168*a^5 + 85/84*a^3 + 31/84*a)])

            Or we may work directly over an algebraically closed field::

                sage: R.<x> = QQbar[]
                sage: r =  1/(x^4+1)
                sage: r.partial_fraction_decomposition()
                (0,
                 [(-0.1767766952966369? - 0.1767766952966369?*I)/(x - 0.7071067811865475? - 0.7071067811865475?*I),
                  (-0.1767766952966369? + 0.1767766952966369?*I)/(x - 0.7071067811865475? + 0.7071067811865475?*I),
                  (0.1767766952966369? - 0.1767766952966369?*I)/(x + 0.7071067811865475? - 0.7071067811865475?*I),
                  (0.1767766952966369? + 0.1767766952966369?*I)/(x + 0.7071067811865475? + 0.7071067811865475?*I)])

            We do the best we can over inexact fields::

                sage: R.<x> = RealField(20)[]
                sage: q = 1/(x^2 + x + 2)^2 + 1/(x-1); q
                (x^4 + 2.0000*x^3 + 5.0000*x^2 + 5.0000*x + 3.0000)/(x^5 + x^4 + 3.0000*x^3 - x^2 - 4.0000)
                sage: whole, parts = q.partial_fraction_decomposition(); parts
                [1.0000/(x - 1.0000), 1.0000/(x^4 + 2.0000*x^3 + 5.0000*x^2 + 4.0000*x + 4.0000)]
                sage: sum(parts)
                (x^4 + 2.0000*x^3 + 5.0000*x^2 + 5.0000*x + 3.0000)/(x^5 + x^4 + 3.0000*x^3 - x^2 - 4.0000)

            TESTS:

            We test partial fraction for irreducible denominators::

                sage: R.<x> = ZZ[]
                sage: q = x^2/(x-1)
                sage: q.partial_fraction_decomposition()
                (x + 1, [1/(x - 1)])
                sage: q = x^10/(x-1)^5
                sage: whole, parts = q.partial_fraction_decomposition()
                sage: whole + sum(parts) == q
                True

            And also over finite fields (see :trac:`6052`, :trac:`9945`)::

                sage: R.<x> = GF(2)[]
                sage: q = (x+1)/(x^3+x+1)
                sage: q.partial_fraction_decomposition()
                (0, [(x + 1)/(x^3 + x + 1)])

                sage: R.<x> = GF(11)[]
                sage: q = x + 1 + 1/(x+1) + x^2/(x^3 + 2*x + 9)
                sage: q.partial_fraction_decomposition()
                (x + 1, [1/(x + 1), x^2/(x^3 + 2*x + 9)])

            And even the rationals::

                sage: (26/15).partial_fraction_decomposition()
                (1, [1/3, 2/5])
                sage: (26/75).partial_fraction_decomposition()
                (-1, [2/3, 3/5, 2/25])

            A larger example::

                sage: S.<t> = QQ[]
                sage: r = t / (t^3+1)^5
                sage: r.partial_fraction_decomposition()
                (0,
                 [-35/729/(t + 1),
                  -35/729/(t^2 + 2*t + 1),
                  -25/729/(t^3 + 3*t^2 + 3*t + 1),
                  -4/243/(t^4 + 4*t^3 + 6*t^2 + 4*t + 1),
                  -1/243/(t^5 + 5*t^4 + 10*t^3 + 10*t^2 + 5*t + 1),
                  (35/729*t - 35/729)/(t^2 - t + 1),
                  (25/729*t - 8/729)/(t^4 - 2*t^3 + 3*t^2 - 2*t + 1),
                  (-1/81*t + 5/81)/(t^6 - 3*t^5 + 6*t^4 - 7*t^3 + 6*t^2 - 3*t + 1),
                  (-2/27*t + 1/9)/(t^8 - 4*t^7 + 10*t^6 - 16*t^5 + 19*t^4 - 16*t^3 + 10*t^2 - 4*t + 1),
                  (-2/27*t + 1/27)/(t^10 - 5*t^9 + 15*t^8 - 30*t^7 + 45*t^6 - 51*t^5 + 45*t^4 - 30*t^3 + 15*t^2 - 5*t + 1)])
                sage: sum(r.partial_fraction_decomposition()[1]) == r
                True

            Some special cases::

                sage: R = Frac(QQ['x']); x = R.gen()
                sage: x.partial_fraction_decomposition()
                (x, [])
                sage: R(0).partial_fraction_decomposition()
                (0, [])
                sage: R(1).partial_fraction_decomposition()
                (1, [])
                sage: (1/x).partial_fraction_decomposition()
                (0, [1/x])
                sage: (1/x+1/x^3).partial_fraction_decomposition()
                (0, [1/x, 1/x^3])

            This was fixed in :trac:`16240`::

                sage: R.<x> = QQ['x']
                sage: p=1/(-x + 1)
                sage: whole,parts = p.partial_fraction_decomposition()
                sage: p == sum(parts)
                True
                sage: p=3/(-x^4 + 1)
                sage: whole,parts = p.partial_fraction_decomposition()
                sage: p == sum(parts)
                True
                sage: p=(6*x^2 - 9*x + 5)/(-x^3 + 3*x^2 - 3*x + 1)
                sage: whole,parts = p.partial_fraction_decomposition()
                sage: p == sum(parts)
                True
            """
            denom = self.denominator()
            whole, numer = self.numerator().quo_rem(denom)
            factors = denom.factor()
            if not self.parent().is_exact():
                # factors not grouped in this case
                all = {}
                for r in factors:
                    all[r[0]] = 0
                for r in factors:
                    all[r[0]] += r[1]
                factors = sorted(all.items())

            # TODO(robertwb): Should there be a category of univariate polynomials?
            from sage.rings.fraction_field_element import FractionFieldElement_1poly_field
            is_polynomial_over_field = isinstance(self, FractionFieldElement_1poly_field)

            running_total = 0
            parts = []
            for r, e in factors:
                powers = [1]
                for ee in range(e):
                    powers.append(powers[-1] * r)
                d = powers[e]
                denom_div_d = denom // d
                # We know the inverse exists as the two are relatively prime.
                n = ((numer % d) * denom_div_d.inverse_mod(d)) % d
                if not is_polynomial_over_field:
                    running_total += n * denom_div_d
                # If the multiplicity is not one, further reduce.
                if decompose_powers:
                    r_parts = []
                    for ee in range(e, 0, -1):
                        n, n_part = n.quo_rem(r)
                        if n_part:
                            r_parts.append(n_part/powers[ee])
                    parts.extend(reversed(r_parts))
                else:
                    parts.append(n/powers[e])

            if not is_polynomial_over_field:
                # remainders not unique, need to re-compute whole to take into
                # account this freedom
                whole = (self.numerator() - running_total) // denom
            return whole, parts

        def derivative(self, *args):
            r"""
            The derivative of this rational function, with respect to variables
            supplied in args.

            Multiple variables and iteration counts may be supplied; see
            documentation for the global derivative() function for more
            details.

            .. SEEALSO::

               :meth:`_derivative`

            EXAMPLES::

                sage: F.<x> = Frac(QQ['x'])
                sage: (1/x).derivative()
                -1/x^2

            ::

                sage: (x+1/x).derivative(x, 2)
                2/x^3

            ::

                sage: F.<x,y> = Frac(QQ['x,y'])
                sage: (1/(x+y)).derivative(x,y)
                2/(x^3 + 3*x^2*y + 3*x*y^2 + y^3)
            """
            from sage.misc.derivative import multi_derivative
            return multi_derivative(self, args)

        def _derivative(self, var=None):
            r"""
            Returns the derivative of this rational function with respect to the
            variable ``var``.

            Over an ring with a working gcd implementation, the derivative of a
            fraction `f/g`, supposed to be given in lowest terms, is computed as
            `(f'(g/d) - f(g'/d))/(g(g'/d))`, where `d` is a greatest common
            divisor of `f` and `g`.

            INPUT:

            - ``var`` - Variable with respect to which the derivative is computed

            OUTPUT:

            - Derivative of ``self`` with respect to ``var``

            .. SEEALSO::

               :meth:`derivative`

            EXAMPLES::

                sage: F.<x> = Frac(QQ['x'])
                sage: t = 1/x^2
                sage: t._derivative(x)
                -2/x^3
                sage: t.derivative()
                -2/x^3

            ::

                sage: F.<x,y> = Frac(QQ['x,y'])
                sage: t = (x*y/(x+y))
                sage: t._derivative(x)
                y^2/(x^2 + 2*x*y + y^2)
                sage: t._derivative(y)
                x^2/(x^2 + 2*x*y + y^2)

            TESTS::

                sage: F.<t> = Frac(ZZ['t'])
                sage: F(0).derivative()
                0
                sage: F(2).derivative()
                0
                sage: t.derivative()
                1
                sage: (1+t^2).derivative()
                2*t
                sage: (1/t).derivative()
                -1/t^2
                sage: ((t+2)/(t-1)).derivative()
                -3/(t^2 - 2*t + 1)
                sage: (t/(1+2*t+t^2)).derivative()
                (-t + 1)/(t^3 + 3*t^2 + 3*t + 1)
            """
            R = self.parent()
            if var in R.gens():
                var = R.ring()(var)

            num = self.numerator()
            den = self.denominator()

            if (num.is_zero()):
                return R.zero()

            if R.is_exact():
                try:
                    numder = num._derivative(var)
                    dender = den._derivative(var)
                    d      = den.gcd(dender)
                    den    = den // d
                    dender = dender // d
                    tnum   = numder * den - num * dender
                    tden   = self.denominator() * den
                    if not tden.is_one() and tden.is_unit():
                        try:
                            tnum = tnum * tden.inverse_of_unit()
                            tden = R.ring().one()
                        except AttributeError:
                            pass
                        except NotImplementedError:
                            pass
                    return self.__class__(R, tnum, tden,
                        coerce=False, reduce=False)
                except AttributeError:
                    pass
                except NotImplementedError:
                    pass
                except TypeError:
                    pass
                num = self.numerator()
                den = self.denominator()

            num = num._derivative(var) * den - num * den._derivative(var)
            den = den**2

            return self.__class__(R, num, den,
                coerce=False, reduce=False)

