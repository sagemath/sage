r"""
Random variables and probability spaces

This introduces a class of random variables, with the focus on
discrete random variables (i.e. on a discrete probability space).
This avoids the problem of defining a measure space and measurable
functions.
"""

# ****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sage.rings.abc
from sage.structure.parent import Parent
from sage.functions.log import log
from sage.misc.functional import sqrt
from sage.rings.rational_field import RationalField
from sage.sets.set import Set
from pprint import pformat

################################################################################
################################################################################


def is_ProbabilitySpace(S):
    from sage.misc.superseded import deprecation
    deprecation(38184,
                "The function is_ProbabilitySpace is deprecated; "
                "use 'isinstance(..., ProbabilitySpace_generic)' instead.")
    return isinstance(S, ProbabilitySpace_generic)


def is_DiscreteProbabilitySpace(S):
    from sage.misc.superseded import deprecation
    deprecation(38184,
                "The function is_DiscreteProbabilitySpace is deprecated; "
                "use 'isinstance(..., DiscreteProbabilitySpace)' instead.")
    return isinstance(S, DiscreteProbabilitySpace)


def is_RandomVariable(X):
    from sage.misc.superseded import deprecation
    deprecation(38184,
                "The function is_RandomVariable is deprecated; "
                "use 'isinstance(..., RandomVariable_generic)' instead.")
    return isinstance(X, RandomVariable_generic)


def is_DiscreteRandomVariable(X):
    from sage.misc.superseded import deprecation
    deprecation(38184,
                "The function is_DiscreteRandomVariable is deprecated; "
                "use 'isinstance(..., DiscreteRandomVariable)' instead.")
    return isinstance(X, DiscreteRandomVariable)

################################################################################
################################################################################

# We could inherit from a functions class here but use Parent


class RandomVariable_generic(Parent):
    """
    A random variable.
    """
    def __init__(self, X, RR):
        if not isinstance(X, ProbabilitySpace_generic):
            raise TypeError("Argument X (= %s) must be a probability space" % X)
        Parent.__init__(self, X)
        self._codomain = RR

    def probability_space(self):
        return self.base()

    def domain(self):
        return self.base()

    def codomain(self):
        return self._codomain

    def field(self):
        return self._codomain


class DiscreteRandomVariable(RandomVariable_generic):
    """
    A random variable on a discrete probability space.
    """
    def __init__(self, X, f, codomain=None, check=False):
        r"""
        Create free binary string monoid on `n` generators.

        INPUT:

        - ``X`` -- a probability space
        - ``f`` -- dictionary such that X[x] = value for x in X
          is the discrete function on X
        """
        if not isinstance(X, DiscreteProbabilitySpace):
            raise TypeError("Argument X (= %s) must be a discrete probability space" % X)
        if check:
            raise NotImplementedError("Not implemented")
        if codomain is None:
            from sage.rings.real_mpfr import RealField
            RR = RealField()
        else:
            RR = codomain
        RandomVariable_generic.__init__(self, X, RR)
        self._function = f

    def __call__(self, x):
        """
        Return the value of the random variable at x.
        """
        RR = self.field()
        try:
            return RR(self._function[x])
        except KeyError:
            # Need some condition for x being a valid domain element:
            #    raise IndexError("Argument x (= %s) is not a valid domain element." % x)
            return RR(0)

    def __repr__(self):
        F = pformat(self._function)
        return "Discrete random variable defined by %s" % F

    def function(self):
        """
        The function defining the random variable.
        """
        return self._function

    def expectation(self):
        r"""
        The expectation of the discrete random variable, namely
        `\sum_{x \in S} p(x) X[x]`, where `X` = ``self`` and
        `S` is the probability space of `X`.
        """
        E = 0
        Omega = self.probability_space()
        for x in self._function:
            E += Omega(x) * self(x)
        return E

    def translation_expectation(self, map):
        r"""
        The expectation of the discrete random variable, namely
        `\sum_{x \in S} p(x) X[e(x)]`, where `X` = self,
        `S` is the probability space of `X`, and
        `e` = map.
        """
        E = 0
        Omega = self.probability_space()
        for x in Omega._function:
            E += Omega(x) * self(map(x))
        return E

    def variance(self):
        r"""
        The variance of the discrete random variable.

        Let `S` be the probability space of `X` = self,
        with probability function `p`, and `E(X)` be the
        expectation of `X`. Then the variance of `X` is:

        .. MATH::

           \mathrm{var}(X) = E((X-E(x))^2) = \sum_{x \in S} p(x) (X(x) - E(x))^2
        """
        Omega = self.probability_space()
        mu = self.expectation()
        var = 0
        for x in self._function:
            var += Omega(x) * (self(x) - mu)**2
        return var

    def translation_variance(self, map):
        r"""
        The variance of the discrete random variable `X \circ e`,
        where `X` = self, and `e` = map.

        Let `S` be the probability space of `X` = self,
        with probability function `p`, and `E(X)` be the
        expectation of `X`. Then the variance of `X` is:

        .. MATH::

           \mathrm{var}(X) = E((X-E(x))^2) = \sum_{x \in S} p(x) (X(x) - E(x))^2
        """
        Omega = self.probability_space()
        mu = self.translation_expectation(map)
        var = 0
        for x in Omega._function:
            var += Omega(x) * (self(map(x)) - mu)**2
        return var

    def covariance(self, other):
        r"""
        The covariance of the discrete random variable X = ``self`` with Y =
        ``other``.

        Let `S` be the probability space of `X` = self,
        with probability function `p`, and `E(X)` be the
        expectation of `X`. Then the variance of `X` is:

        .. MATH::

                     \text{cov}(X,Y) = E((X-E(X)\cdot (Y-E(Y)) = \sum_{x \in S} p(x) (X(x) - E(X))(Y(x) - E(Y))
        """
        Omega = self.probability_space()
        if Omega != other.probability_space():
            raise ValueError("Argument other (= %s) must be defined on the same probability space." % other)
        muX = self.expectation()
        muY = other.expectation()
        cov = 0
        for x in self._function:
            cov += Omega(x)*(self(x) - muX)*(other(x) - muY)
        return cov

    def translation_covariance(self, other, map):
        r"""
        The covariance of the probability space X = ``self`` with image of Y =
        ``other`` under the given map of the probability space.

        Let `S` be the probability space of `X` = self,
        with probability function `p`, and `E(X)` be the
        expectation of `X`. Then the variance of `X` is:

        .. MATH::

                     \text{cov}(X,Y) = E((X-E(X)\cdot (Y-E(Y)) = \sum_{x \in S} p(x) (X(x) - E(X))(Y(x) - E(Y))
        """
        Omega = self.probability_space()
        if Omega != other.probability_space():
            raise ValueError("Argument other (= %s) must be defined on the same probability space." % other)
        muX = self.expectation()
        muY = other.translation_expectation(map)
        cov = 0
        for x in Omega._function:
            cov += Omega(x)*(self(x) - muX)*(other(map(x)) - muY)
        return cov

    def standard_deviation(self):
        r"""
        The standard deviation of the discrete random variable.

        Let `S` be the probability space of `X` = self,
        with probability function `p`, and `E(X)` be the
        expectation of `X`. Then the standard deviation of
        `X` is defined to be

        .. MATH::

                     \sigma(X) = \sqrt{ \sum_{x \in S} p(x) (X(x) - E(x))^2}
        """
        return sqrt(self.variance())

    def translation_standard_deviation(self, map):
        r"""
        The standard deviation of the translated discrete random variable
        `X \circ e`, where `X` = ``self`` and `e` =
        map.

        Let `S` be the probability space of `X` = ``self``,
        with probability function `p`, and `E(X)` be the
        expectation of `X`. Then the standard deviation of
        `X` is defined to be

        .. MATH::

            \sigma(X) = \sqrt{ \sum_{x \in S} p(x) (X(x) - E(x))^2}
        """
        return sqrt(self.translation_variance(map))

    def correlation(self, other):
        """
        The correlation of the probability space X = ``self`` with Y = ``other``.
        """
        cov = self.covariance(other)
        sigX = self.standard_deviation()
        sigY = other.standard_deviation()
        if sigX == 0 or sigY == 0:
            raise ValueError("Correlation not defined if standard deviations are not both nonzero.")
        return cov/(sigX*sigY)

    def translation_correlation(self, other, map):
        """
        The correlation of the probability space X = ``self`` with image of Y =
        ``other`` under map.
        """
        cov = self.translation_covariance(other, map)
        sigX = self.standard_deviation()
        sigY = other.translation_standard_deviation(map)
        if sigX == 0 or sigY == 0:
            raise ValueError("Correlation not defined if standard deviations are not both nonzero.")
        return cov/(sigX*sigY)

################################################################################
################################################################################


class ProbabilitySpace_generic(RandomVariable_generic):
    r"""
    A probability space.
    """
    def __init__(self, domain, RR):
        """
        A generic probability space on given domain space and codomain
        ring.
        """
        if isinstance(domain, list):
            domain = tuple(domain)
        if not isinstance(domain, tuple):
            raise TypeError("Argument domain (= %s) must be a list, tuple, or set containing." % domain)
        self._domain = domain
        RandomVariable_generic.__init__(self, self, RR)

    def domain(self):
        return self._domain


class DiscreteProbabilitySpace(ProbabilitySpace_generic,DiscreteRandomVariable):
    r"""
    The discrete probability space
    """
    def __init__(self, X, P, codomain=None, check=False):
        r"""
        Create the discrete probability space with probabilities on the
        space X given by the dictionary P with values in the field
        real_field.

        EXAMPLES::

            sage: S = [ i for i in range(16) ]
            sage: P = {}
            sage: for i in range(15): P[i] = 2^(-i-1)
            sage: P[15] = 2^-15
            sage: X = DiscreteProbabilitySpace(S,P)
            sage: sum(X.function().values())
            1
            sage: X.domain()
            (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
            sage: X.set()
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}
            sage: X.entropy().n()                                                       # needs sage.libs.pari
            1.99993896484375

        A probability space can be defined on any list of elements::

            sage: AZ = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
            sage: S = [ AZ[i] for i in range(26) ]
            sage: P = { 'A':1/2, 'B':1/4, 'C':1/4 }
            sage: X = DiscreteProbabilitySpace(S,P)
            sage: X
            Discrete probability space defined by {'A': 1/2, 'B': 1/4, 'C': 1/4}
            sage: X.entropy().n()                                                       # needs sage.libs.pari
            1.50000000000000
        """
        if codomain is None:
            from sage.rings.real_mpfr import RealField
            codomain = RealField()
        if not isinstance(codomain, sage.rings.abc.RealField) and not isinstance(codomain, RationalField):
            raise TypeError("Argument codomain (= %s) must be the reals or rationals" % codomain)
        if check:
            one = sum(P.values())
            if isinstance(codomain, RationalField):
                if not one == 1:
                    raise TypeError("Argument P (= %s) does not define a probability function")
            else:
                if not abs(one - 1) < 2 ** (-codomain.precision() + 1):
                    raise TypeError("Argument P (= %s) does not define a probability function")
        ProbabilitySpace_generic.__init__(self, X, codomain)
        DiscreteRandomVariable.__init__(self, self, P, codomain, check)

    def __repr__(self):
        """
        TESTS::

            sage: S = list(range(4))
            sage: P = {i: 2^(-i-1) for i in range(3)}
            sage: P[4] = 2^-3
            sage: DiscreteProbabilitySpace(S,P)
            Discrete probability space defined by {0: 1/2, 1: 1/4, 2: 1/8, 4: 1/8}
        """
        F = pformat(self.function())
        return "Discrete probability space defined by %s" % F

    def set(self):
        r"""
        The set of values of the probability space taking possibly nonzero
        probability (a subset of the domain).
        """
        return Set(self.function())

    def entropy(self):
        """
        The entropy of the probability space.
        """
        def neg_xlog2x(p):
            if p == 0:
                return 0
            else:
                return -p*log(p,2)
        p = self.function()
        return sum([neg_xlog2x(p[x]) for x in p])
