# -*- coding: utf-8 -*-
r"""
Valuations on polynomial rings based on `\phi`-adic expansions

This file implements a base class for discrete valuations on polynomial rings,
defined by a `\phi`-adic expansion.

AUTHORS:

- Julian Rüth (2013-04-15): initial version

EXAMPLES:

The :mod:`Gauss valuation <sage.rings.valuation.gauss_valuation>` is a simple example of a valuation that relies on
`\phi`-adic expansions::

    sage: R.<x> = QQ[]
    sage: v = GaussValuation(R, QQ.valuation(2))

In this case, `\phi = x`, so the expansion simply lists the coefficients of the
polynomial::

    sage: f = x^2 + 2*x + 2
    sage: list(v.coefficients(f))
    [2, 2, 1]

Often only the first few coefficients are necessary in computations, so for
performance reasons, coefficients are computed lazily::

    sage: v.coefficients(f)
    <generator object ...coefficients at 0x...>

Another example of a :class:`DevelopingValuation` is an :mod:`augmented
valuation <sage.rings.valuation.augmented_valuation>`::

    sage: w = v.augmentation(x^2 + x + 1, 3)

Here, the expansion lists the remainders of repeated division by `x^2 + x + 1`::

    sage: list(w.coefficients(f))
    [x + 1, 1]

"""
#*****************************************************************************
#       Copyright (C) 2013-2017 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .valuation import DiscretePseudoValuation
from sage.misc.abstract_method import abstract_method

from sage.misc.cachefunc import cached_method


class DevelopingValuation(DiscretePseudoValuation):
    r"""
    Abstract base class for a discrete valuation of polynomials defined over
    the polynomial ring ``domain`` by the `\phi`-adic development.

    EXAMPLES::

        sage: R.<x> = QQ[]
        sage: v = GaussValuation(R, QQ.valuation(7))

    TESTS::

        sage: TestSuite(v).run() # long time

    """
    def __init__(self, parent, phi):
        r"""
        TESTS::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, QQ.valuation(7))
            sage: from sage.rings.valuation.developing_valuation import DevelopingValuation
            sage: isinstance(v, DevelopingValuation)
            True

        """
        DiscretePseudoValuation.__init__(self, parent)

        domain = parent.domain()
        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        if not is_PolynomialRing(domain) or not domain.ngens() == 1:
            raise TypeError("domain must be a univariate polynomial ring but %r is not"%(domain,))

        phi = domain.coerce(phi)
        if phi.is_constant() or not phi.is_monic():
            raise ValueError("phi must be a monic non-constant polynomial but %r is not"%(phi,))

        self._phi = phi

    def phi(self):
        r"""
        Return the polynomial `\phi`, the key polynomial of this valuation.

        EXAMPLES::

            sage: R = Zp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.phi()
            (1 + O(2^5))*x

        """
        return self._phi

    def effective_degree(self, f, valuations=None):
        r"""
        Return the effective degree of ``f`` with respect to this valuation.

        The effective degree of `f` is the largest `i` such that the valuation
        of `f` and the valuation of `f_i\phi^i` in the development `f=\sum_j
        f_j\phi^j` coincide (see [Mac1936II]_ p.497.)

        INPUT:

        - ``f`` -- a non-zero polynomial in the domain of this valuation

        EXAMPLES::

            sage: R = Zp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.effective_degree(x)
            1
            sage: v.effective_degree(2*x + 1)
            0

        """
        f = self.domain().coerce(f)

        if f.is_zero():
            raise ValueError("the effective degree is only defined for non-zero polynomials")

        if valuations is None:
            valuations = list(self.valuations(f))
        v = min(valuations)
        return [i for i,w in enumerate(valuations) if w == v][-1]

    @cached_method
    def _pow(self, f, e, error, effective_degree):
        r"""
        Return `f^e`.

        This method does not compute the exact value of `f^e` but only an
        element that differs from the correct result by an error with valuation
        at least ``error``. The output is assumed to have at most
        ``effective_degree``. If the effective degree is higher than
        ``effective_degree``, then the result may not be correct.

        EXAMPLES::

            sage: R = Zp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v._pow(2*x + 1, 10, effective_degree=0, error=5)
            1 + O(2^5)

        """
        if e == 0:
            return self.domain().one()
        if e == 1:
            return self.simplify(f, error=error)
        if e % 2 == 0:
            return self._pow(self.simplify(f*f, error=error*2/e, effective_degree=effective_degree*2/e),
                             e//2, error=error, effective_degree=effective_degree)
        else:
            return self.simplify(f*self._pow(f, e-1, error=error*(e-1)/e, effective_degree=effective_degree*(e-1)/e),
                                 error=error, effective_degree=effective_degree)

    def coefficients(self, f):
        r"""
        Return the `\phi`-adic expansion of ``f``.

        INPUT:

        - ``f`` -- a monic polynomial in the domain of this valuation

        OUTPUT:

        An iterator `f_0, f_1, \dots, f_n` of polynomials in the domain of this
        valuation such that `f=\sum_i f_i\phi^i`

        EXAMPLES::

            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: f = x^2 + 2*x + 3
            sage: list(v.coefficients(f)) # note that these constants are in the polynomial ring
            [1 + 2 + O(2^5), 2 + O(2^6), 1 + O(2^5)]
            sage: v = v.augmentation( x^2 + x + 1, 1)
            sage: list(v.coefficients(f))
            [(1 + O(2^5))*x + 2 + O(2^5), 1 + O(2^5)]

        """
        domain = self.domain()
        f = domain.coerce(f)

        if f.degree() < self.phi().degree():
            yield f
        elif self.phi().degree() == 1:
            if self.phi() != domain.gen() or not domain.is_exact():
                f = f(domain.gen() - self.phi()[0])
            for c in f.coefficients(sparse=False):
                yield domain(c)
        else:
            while f.degree() >= 0:
                f,r = self._quo_rem(f)
                yield r

    def _quo_rem(self, f):
        r"""
        Return the quotient and remainder of ``f`` divided by the key
        polynomial :meth:`phi`.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: v = GaussValuation(S, QQ.valuation(2))
            sage: v._quo_rem(x^2 + 1)
            (x, 1)

        """
        return f.quo_rem(self.phi())

    def newton_polygon(self, f, valuations=None):
        r"""
        Return the newton polygon of the `\phi`-adic development of ``f``.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        EXAMPLES::

            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: f = x^2 + 2*x + 3
            sage: v.newton_polygon(f)
            Finite Newton polygon with 2 vertices: (0, 0), (2, 0)

            sage: v = v.augmentation( x^2 + x + 1, 1)
            sage: v.newton_polygon(f)
            Finite Newton polygon with 2 vertices: (0, 0), (1, 1)
            sage: v.newton_polygon( f * v.phi()^3 )
            Finite Newton polygon with 2 vertices: (3, 3), (4, 4)

        """
        f = self.domain().coerce(f)

        from sage.geometry.newton_polygon import NewtonPolygon
        if valuations is None:
            valuations = self.valuations(f)
        return NewtonPolygon(list(enumerate(valuations)))

    def _call_(self, f):
        r"""
        Evaluate this valuation at ``f``.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        EXAMPLES::

            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: f = x^2 + 2*x + 3
            sage: v(f)
            0

            sage: v = v.augmentation( x^2 + x + 1, 1)
            sage: v(f)
            0
            sage: v(f * v.phi()^3 )
            3
            sage: v(S.zero())
            +Infinity

        """
        f = self.domain().coerce(f)

        from sage.rings.infinity import infinity
        if f.is_zero():
            return infinity

        ret = infinity
        for v in self.valuations(f, call_error=True):
            if ret is infinity or (v is not infinity and v < ret):
                # "ret is infinity" is redundant but much faster than < when ret is infinite
                ret = v
        return ret

    @abstract_method
    def valuations(self, f):
        r"""
        Return the valuations of the `f_i\phi^i` in the expansion `f=\sum f_i\phi^i`.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        OUTPUT:

        A list, each entry a rational numbers or infinity, the valuations of
        `f_0, f_1\phi, \dots`

        EXAMPLES::

            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S, R.valuation())
            sage: f = x^2 + 2*x + 16
            sage: list(v.valuations(f))
            [4, 1, 0]

        """

    def _test_effective_degree(self, **options):
        r"""
        Test the correctness of :meth:`effective_degree`.

        EXAMPLES::

            sage: R = Zp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v._test_effective_degree()
        
        """
        tester = self._tester(**options)
        S = tester.some_elements(self.domain().base_ring().some_elements())
        for x in S:
            if x == 0:
                continue
            tester.assertEqual(self.effective_degree(x), 0)
