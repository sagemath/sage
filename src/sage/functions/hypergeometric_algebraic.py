r"""
Hypergeometric functions over arbitrary rings

When the given variable `x` is not symbolic but lies in a polynomial
ring or a power series ring, the function
:func:`~sage.functions.hypergeometric.hypergeometric`
returns an instance of the class :class:`HypergeometricAlgebraic`::

    sage: S.<x> = QQ[]
    sage: f = hypergeometric([1/9, 4/9, 5/9], [1/3, 1], x)
    sage: f.parent()
    Hypergeometric functions in x over Rational Field

Below, we illustrate the main features provided by this class.
We introduce two additional hypergeometric series which will serve
as running examples::

    sage: g = hypergeometric([1/2, 5/6, 1], [5/3, 2], x)
    sage: h = hypergeometric([1/5, 1/5, 1/5, 1/5], [1/3, 27/5 - 1], x)

.. RUBRIC:: Hypergeometric functions over `\QQ`

A series `s(x)` is said globally bounded when it has positive radius
of convergence and there exist integers `a` and `b` such that `a \cdot
s(bx)` has integral coefficients.
The method :meth:`~HypergeometricAlgebraic_QQ.is_globally_bounded` checks
when this property is satisfied::

    sage: f.is_globally_bounded()
    True
    sage: g.is_globally_bounded()
    True
    sage: h.is_globally_bounded()
    False

More generally, the method
:meth:`HypergeometricAlgebraic_QQ.good_reduction_primes` returns the
set of primes modulo which the hypergeometric function can be reduced::

    sage: f.good_reduction_primes()
    Set of all prime numbers with 3 excluded: 2, 5, 7, 11, ...
    sage: g.good_reduction_primes()
    Set of all prime numbers with 2 excluded: 3, 5, 7, 11, ...
    sage: h.good_reduction_primes()
    Set of prime numbers congruent to 1, 8, 11 modulo 15 with 3, 17 included and 11 excluded: 3, 17, 23, 31, ...

On a different note, the method :meth:`~HypergeometricAlgebraic_QQ.is_algebraic`
checks whether an hypergeometric series defines an algebraic function
over `\QQ(x)`::

    sage: f.is_algebraic()
    False
    sage: g.is_algebraic()
    True
    sage: h.is_algebraic()
    False

.. RUBRIC:: Hypergeometric functions over finite fields

When `p` is a prime of good reduction of an hypergeometric function, we
can reduce the latter modulo `p` using the mod operator (``%``)::

    sage: f19 = f % 19
    sage: f19
    hypergeometric((1/9, 4/9, 5/9), (1/3, 1), x)
    sage: f19.base_ring()
    Finite Field of size 19

A remarkable feature of hypergeometric functions over finite fields is
that they are always algebraic!
The method :meth:`~HypergeometricAlgebraic_GFp.annihilating_ore_polynomial`
returns an annihilating polynomial (in the Frobenius)::

    sage: f19.annihilating_ore_polynomial()
    (18*x^76 + 13*x^57 + 6*x^38 + 17*x^19 + 12)*Frob^2 +
    (12*x^38 + 11*x^32 + 10*x^31 + ... + 18*x^12 + 7)*Frob +
    x^30 + 16*x^29 + 9*x^28 + ... + 6*x^13 + x^12

One subtlety is positive characteristic is that different set of
parameters may lead to the same series::

    sage: T.<y> = GF(13)[]
    sage: h1 = hypergeometric([1/12, 1/4], [1/2], y)
    sage: h2 = hypergeometric([1/12, 1/6], [1/3], y)
    sage: h1.power_series(500)
    1 + 6*y + 6*y^13 + 10*y^14 + 6*y^169 + 10*y^170 + 10*y^182 + 8*y^183 + O(y^500)
    sage: h2.power_series(500)
    1 + 6*y + 6*y^13 + 10*y^14 + 6*y^169 + 10*y^170 + 10*y^182 + 8*y^183 + O(y^500)

The method :meth:`~HypergeometricAlgebraic_GFp.is_equal_as_series` checks
when this happens::

    sage: h1.is_equal_as_series(h2)
    True

.. RUBRIC:: Hypergeometric functions over `p`-adic fields

Some methods related to `p`-adic properties of hypergeometric series
are also available,. This includes the computation of the `p`-adic
valuation::

    sage: hp3 = h.change_ring(Qp(3))
    sage: hp3.valuation()
    0

We can also compute the `p`-adic radius of convergence::

    sage: hp3.log_radius_of_convergence()
    2

Here, the log radius of convergence refers to the exponent on `p`
of the actual radius of convergence; in our example, the `p`-adic
radius of convergence of `h` is then `p^2`.

Evaluation of hypergeometric series at `p`-adic arguments also
works::

    sage: hp3(1/3)
    3 + 3^4 + 2*3^5 + 2*3^7 + 3^8 + 2*3^9 + 2*3^10 + 3^11 + 3^12 + 3^13 + 2*3^14 + 2*3^15 + 3^16 + 3^17 + 3^19 + O(3^20)
    sage: hp3(1/9)
    Traceback (most recent call last):
    ...
    ValueError: outside the domain of convergence

AUTHORS:

- Xavier Caruso, Florian Fürnsinn (2026-02): initial version
"""

# ***************************************************************************
#    Copyright (C) 2026 Xavier Caruso <xavier.caruso@normalesup.org>
#                       Florian Fürnsinn <florian.fuernsinn@univie.ac.at>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

import operator

from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex
from sage.misc.latex import latex_variable_name

from sage.misc.misc_c import prod
from sage.functions.other import floor
from sage.misc.functional import log
from sage.functions.hypergeometric import hypergeometric
from sage.arith.misc import gcd
from sage.matrix.constructor import matrix

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.element import coerce_binop
from sage.structure.category_object import normalize_names

from sage.categories.action import Action
from sage.categories.pushout import pushout
from sage.categories.map import Map
from sage.categories.sets_cat import Sets
from sage.categories.finite_fields import FiniteFields

from sage.matrix.special import companion_matrix
from sage.matrix.special import identity_matrix
from sage.combinat.subset import Subsets
from sage.geometry.newton_polygon import NewtonPolygon

from sage.rings.infinity import infinity
from sage.symbolic.ring import SR
from sage.sets.primes import Primes
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.padics.padic_generic import pAdicGeneric
from sage.rings.padics.factory import Qp
from sage.rings.number_field.number_field import CyclotomicField

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.lazy_series_ring import LazyPowerSeriesRing
from sage.rings.tate_algebra import TateAlgebra
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing

from sage.functions.hypergeometric_parameters import HypergeometricParameters


# Do we want to implement polynomial linear combination
# of hypergeometric functions?
# Advantages:
#  . reductions mod p of hypergeometric functions have this form in general
#  . many methods can be extended to this context
# Difficulty:
#  . not sure we can handle easily simplifications!

class HypergeometricAlgebraic(Element):
    r"""
    Class for (scalar multiples of) hypergeometric functions over arbitrary base rings.
    """
    def __init__(self, parent, arg1, arg2=None, scalar=None, check=True):
        r"""
        Initialize this hypergeometric function.

        INPUT:

        - ``parent`` -- the parent of this function

        - ``arg1``, ``arg2`` -- arguments defining this hypergeometric
          function, they can be:
          - the top and bottom parameters
          - a hypergeometric function and ``None``
          - an instance of the class :class:`HypergeometricParameters` and ``None``

        - ``scalar`` -- an element in the base ring, the scalar by
          which the hypergeometric function is multiplied

        TESTS::

            sage: S.<x> = QQ[]
            sage: h = hypergeometric((1/2, 1/3), (1,), x)
            sage: type(h)
            <class 'sage.functions.hypergeometric_algebraic.HypergeometricFunctions_with_category.element_class'>
            sage: TestSuite(h).run()

        ::

            sage: hypergeometric([-1], [-2], x)
            hypergeometric((-1,), (-2,), x)
            sage: hypergeometric([-2], [-1], x)
            Traceback (most recent call last):
            ...
            ValueError: the parameters (-2,) and (-1,) do not define a hypergeometric function
        """
        Element.__init__(self, parent)
        base = parent.base_ring()
        if scalar is None:
            scalar = base.one()
        else:
            scalar = base(scalar)
        if scalar == 0:
            parameters = None
        elif isinstance(arg1, HypergeometricAlgebraic):
            parameters = arg1._parameters
            scalar *= base(arg1._scalar)
        elif isinstance(arg1, HypergeometricParameters):
            parameters = arg1
        else:
            parameters = HypergeometricParameters(arg1, arg2)
        char = self.parent()._char
        if check and scalar:
            try:
                _ = parameters.degree()
            except ValueError:
                raise ValueError("the parameters %s and %s do not define a hypergeometric function"
                              % (parameters.top, parameters.bottom[:-1]))
            if char > 0:
                val, _, _ = parameters.valuation_position(char)
                if val < 0:
                    raise ValueError("the parameters %s and %s do not define a hypergeometric function in characteristic %s"
                                  % (parameters.top, parameters.bottom[:-1], char))
        self._scalar = scalar
        self._parameters = parameters
        self._coeffs = [scalar]
        self._coeffs_enriched = [[QQ(1), 0]]
        self._char = char

    def __hash__(self):
        r"""
        Return a hash of this hypergeometric function.

        TESTS::

            sage: S.<x> = QQ[]
            sage: h = hypergeometric((1/2, 1/3), (1,), x)
            sage: hash(h)  # random
            -5906731172464693436
        """
        return hash((self.base_ring(), self._parameters, self._scalar))

    @coerce_binop
    def is_equal_symbolically(self, other):
        r"""
        Return whether if the parameters defining the hypergeometric
        series ``self`` and ``other`` are the same.

        INPUT:

        - ``other`` -- an hypergeometric function

        EXAMPLES:

        The order of the parameters is not relevant::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/12, 1/6], [1/3], x)
            sage: g = hypergeometric([1/6, 1/12], [1/3], x)
            sage: f.is_equal_symbolically(g)
            True

        ::

            sage: h = hypergeometric([1/12, 1/4], [1/2], x)
            sage: g.is_equal_symbolically(h)
            False

        We emphasize that two hypergeometric functions are considered
        as different as soon as they have different parameters even if
        they define the same series::

            sage: Fq = GF(13)
            sage: g13 = g % 13
            sage: h13 = h % 13
            sage: g13 == h13
            False
            sage: g13.power_series(500)
            1 + 6*x + 6*x^13 + 10*x^14 + 6*x^169 + 10*x^170 + 10*x^182 + 8*x^183 + O(x^500)
            sage: h13.power_series(500)
            1 + 6*x + 6*x^13 + 10*x^14 + 6*x^169 + 10*x^170 + 10*x^182 + 8*x^183 + O(x^500)

        .. SEEALSO::

            :meth:`is_equal_as_series`
        """
        return self._parameters == other._parameters and self._scalar == other._scalar

    @coerce_binop
    def is_equal_as_series(self, other):
        r"""
        Return whether ``self`` and ``other`` define the same series.

        INPUT:

        - ``other`` -- an hypergeometric function over the same base

        EXAMPLES::

            sage: S.<x> = GF(13)[]
            sage: f = hypergeometric([1/12, 1/6], [1/3], x)
            sage: g = hypergeometric([1/12, 1/4], [1/2], x)
            sage: f.is_equal_as_series(g)
            True

        Note that this method is not implemented over all bases::

            sage: S.<x> = Integers(169)[]
            sage: f = hypergeometric([1/12, 1/6], [1/3], x)
            sage: g = hypergeometric([1/12, 1/4], [1/2], x)
            sage: f.is_equal_as_series(g)
            Traceback (most recent call last):
            ...
            NotImplementedError: equality as series is not implemented over Ring of integers modulo 169

        .. SEEALSO::

            :meth:`is_equal_symbolically`
        """
        if self._scalar != other._scalar:
            return False
        char = self.parent()._char
        if char == 0:
            return self.is_equal_symbolically(other)
        if char.is_prime():
            H = self.parent().change_ring(FiniteField(char))
            hs = H(self._parameters)
            ho = H(other._parameters)
            return hs.is_equal_as_series(ho)
        raise NotImplementedError("equality as series is not implemented over %s" % self.base_ring())

    def __eq__(self, other):
        r"""
        Return whether ``self`` is equal to ``other`` according to the
        equality convention defined in the parent.

        INPUT:

        - ``other`` -- an hypergeometric function

        TESTS::

            sage: A.<y> = GF(13)[]
            sage: S.<x> = A[]
            sage: f = y * hypergeometric([1/12, 1/6], [1/3], x)
            sage: g = y * hypergeometric([1/12, 1/4], [1/2], x)
            sage: f == g  # symbolic equality
            False

        ::

            sage: ff = y * hypergeometric([1/12, 1/6], [1/3], x, symbolic_equality=False)
            sage: gg = y * hypergeometric([1/12, 1/4], [1/2], x, symbolic_equality=False)
            sage: ff == gg  # equality as series
            True
        """
        if not isinstance(other, HypergeometricAlgebraic):
            return False
        if self.parent()._symbolic_equality and other.parent()._symbolic_equality:
            return self.is_equal_symbolically(other)
        return self.is_equal_as_series(other)

    def _repr_(self):
        r"""
        Return a string representation of this hypergeometric function.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f    # indirect doctest
            hypergeometric((1/3, 2/3), (1/2,), x)
            sage: 2*f  # indirect doctest
            2*hypergeometric((1/3, 2/3), (1/2,), x)
            sage: 0*f  # indirect doctest
            0
        """
        if self._parameters is None:
            return "0"
        scalar = self._scalar
        if scalar == 1:
            s = ""
        elif scalar._is_atomic():
            scalar = str(scalar)
            if scalar == "-1":
                s = "-"
            else:
                s = scalar + "*"
        else:
            s = "(%s)*" % scalar
        s += "hypergeometric(%s, %s, %s)" % (self.top(), self.bottom(), self.parent().variable_name())
        return s

    def _latex_(self):
        r"""
        Return a LaTex representation of this hypergeometric function.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f._latex_()
            '\\,_{2} F_{1} \\left(\\begin{matrix} \\frac{1}{3},\\frac{2}{3}\\\\\\frac{1}{2}\\end{matrix}; x \\right)'
        """
        if self._parameters is None:
            return "0"
        scalar = self._scalar
        if scalar == 1:
            s = ""
        elif scalar._is_atomic():
            scalar = latex(scalar)
            if scalar == "-1":
                s = "-"
            else:
                s = scalar
        else:
            s = r"\left(%s\right)" % scalar
        top = self.top()
        bottom = self.bottom()
        s += r"\,_{%s} F_{%s} " % (len(top), len(bottom))
        s += r"\left(\begin{matrix} "
        s += ",".join(latex(a) for a in top)
        s += r"\\"
        s += ",".join(latex(b) for b in bottom)
        s += r"\end{matrix}; %s \right)" % self.parent().latex_variable_name()
        return s

    def base_ring(self):
        r"""
        Return the ring over which this hypergeometric function is defined.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.base_ring()
            Rational Field

        ::

            sage: T.<y> = Qp(5)[]
            sage: g = hypergeometric([1/3, 2/3], [1/2], y)
            sage: g.base_ring()
            5-adic Field with capped relative precision 20

        ::

            sage: U.<z> = GF(5)[]
            sage: h = hypergeometric([1/3, 2/3], [1/2], z)
            sage: h.base_ring()
            Finite Field of size 5

        ::

            sage: V.<w> = CC[]
            sage: k = hypergeometric([1/3, 2/3], [1/2], w)
            sage: k.base_ring()
            Complex Field with 53 bits of precision
        """
        return self.parent().base_ring()

    def top(self):
        r"""
        Return the top parameters of this hypergeometric function.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.top()
            (1/3, 2/3)
        """
        return self._parameters.top

    def bottom(self):
        r"""
        Return the bottom parameters of this hypergeometric function (excluding
        the extra ``1``).

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.bottom()
            (1/2,)
        """
        return self._parameters.bottom[:-1]

    def scalar(self):
        r"""
        Return the scalar of this hypergeometric function.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.scalar()
            1
            sage: g = 4*f
            sage: g.scalar()
            4
        """
        return self._scalar

    def change_ring(self, R):
        r"""
        Return this hypergeometric function with changed base ring.

        INPUT:

        - ``R`` -- a commutative ring

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.base_ring()
            Rational Field
            sage: g = f.change_ring(Qp(5))
            sage: g.base_ring()
            5-adic Field with capped relative precision 20
        """
        H = self.parent().change_ring(R)
        return H(self._parameters, None, self._scalar)

    def change_variable_name(self, name):
        r"""
        Return this hypergeometric function with changed variable name

        INPUT:

        - ``name`` -- a string, the new variable name

        EXAMPLES::

            sage: S.<x> = Qp(5)[]
            sage: T.<y> = Qp(5)[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f
            hypergeometric((1/3, 2/3), (1/2,), x)
            sage: g = f.change_variable_name('y')
            sage: g
            hypergeometric((1/3, 2/3), (1/2,), y)
        """
        H = self.parent().change_variable_name(name)
        return H(self._parameters, None, self._scalar)

    def _add_(self, other):
        r"""
        Return the (formal) sum of the hypergeometric function
        and ``other``.

        INPUT:

        - ``other`` -- a hypergeometric function

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: g = 1/2 * hypergeometric([1/3, 2/3], [1/2], x)
            sage: h = hypergeometric([1/5, 2/5], [3/5], x)
            sage: f + g
            3/2*hypergeometric((1/3, 2/3), (1/2,), x)
            sage: f + h
            hypergeometric((1/3, 2/3), (1/2,), x) + hypergeometric((1/5, 2/5), (3/5,), x)

        ::

            sage: f + cos(x)
            cos(x) + hypergeometric((1/3, 2/3), (1/2,), x)
        """
        if self._parameters is None:
            return other
        if isinstance(other, HypergeometricAlgebraic):
            if other._parameters is None:
                return self
            if self._parameters == other._parameters:
                scalar = self._scalar + other._scalar
                return self.parent()(self._parameters, scalar=scalar)
        return SR(self) + SR(other)

    def _neg_(self):
        r"""
        Return the negative of this hypergeometric function.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = 2*hypergeometric([1/3, 2/3], [1/2], x)
            sage: -f
            -2*hypergeometric((1/3, 2/3), (1/2,), x)
        """
        if self._parameters is None:
            return self
        return self.parent()(self._parameters, scalar=-self._scalar)

    def _sub_(self, other):
        r"""
        Return the (formal) difference of the hypergeometric function
        with ``other``.

        INPUT:

        - ``other`` -- a hypergeometric function or a formal expression

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: g = 1/2 * hypergeometric([1/3, 2/3], [1/2], x)
            sage: h = hypergeometric([1/5, 2/5], [3/5], x)
            sage: f - g
            1/2*hypergeometric((1/3, 2/3), (1/2,), x)
            sage: f - h
            hypergeometric((1/3, 2/3), (1/2,), x) - hypergeometric((1/5, 2/5), (3/5,), x)

        ::

            sage: f - sin(x)
            hypergeometric((1/3, 2/3), (1/2,), x) - sin(x)
        """
        if self._parameters is None:
            return other
        if isinstance(other, HypergeometricAlgebraic):
            if other._parameters is None:
                return self
            if self._parameters == other._parameters:
                scalar = self._scalar - other._scalar
                return self.parent()(self._parameters, scalar=scalar)
        return SR(self) - SR(other)

    def _mul_(self, other):
        r"""
        Return the (formal) product of the hypergeometric function
        and ``other``

        INPUT:

        - ``other`` -- a hypergeometric function or a formal expression

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: g = 1/2 * hypergeometric([1/3, 2/3], [1/2], x)
            sage: h = hypergeometric([1/5, 2/5], [3/5], x)
            sage: f*g
            1/2*hypergeometric((1/3, 2/3), (1/2,), x)^2
            sage: f*h
            hypergeometric((1/3, 2/3), (1/2,), x)*hypergeometric((1/5, 2/5), (3/5,), x)

        ::

            sage: sin(x)*f + x
            hypergeometric((1/3, 2/3), (1/2,), x)*sin(x) + x
        """
        return SR(self) * SR(other)

    def __call__(self, x):
        r"""
        Return the value of this hypergeometric function at ``x``.

        INPUT:

        - ``x`` -- an element

        EXAMPLES::

            sage: S.<x> = RR[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f(0.5)
            1.36602540378444

        ::

            sage: g = 2*f
            sage: g(0.2)
            2.20941633798502
        """
        scalar = self._scalar
        if scalar == 0:
            return self.base_ring().zero()
        X = SR('X')
        h = hypergeometric(self.top(), self.bottom(), X)
        if scalar != 1:
            h *= scalar
        return h(X=x)

    def _compute_coeffs(self, prec):
        r"""
        Compute the coefficients of the series representation of this
        hypergeometric function up to a given precision, and store
        them in ``self._coeffs``.

        INPUT:

        - ``prec`` -- a positive integer

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f._coeffs
            [1]
            sage: f._compute_coeffs(3)
            sage: f._coeffs
            [1, 4/9, 80/243]

        ::

            sage: R.<y> = GF(13)[]
            sage: S.<x> = R[]
            sage: f = y * hypergeometric([1/12, 1/6], [1/3], x)
            sage: f.power_series(3)
            y + 6*y*x + O(x^3)
        """
        start = len(self._coeffs) - 1
        c, z = self._coeffs_enriched[-1]
        R = self.base_ring()
        scalar = self._scalar
        for i in range(start, prec - 1):
            for a in self._parameters.top:
                if a + i == 0:
                    z += 1
                else:
                    c *= a + i
            for b in self._parameters.bottom:
                if b + i == 0:
                    z -= 1
                else:
                    c /= b + i
            if z < 0:
                raise RuntimeError
            elif z > 0:
                self._coeffs.append(R.zero())
            else:
                self._coeffs.append(scalar * R(c))
            self._coeffs_enriched.append([c, z])

    def __getitem__(self, n):
        r"""
        Return the ``n``-th coefficient of the series representation of this
        hypergeometric function.

        INPUT:

        - ``n`` -- a non-negative integer

        EXAMPLES:

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f[9]
            409541017600/2541865828329
            sage: g = f % 5
            sage: g[9]
            0
        """
        self._compute_coeffs(n+1)
        S = self.base_ring()
        return S(self._coeffs[n])

    def coefficient(self, n):
        r"""
        Return the ``n``-th coefficient of the series representation of this
        hypergeometric function.

        INPUT:

        - ``n`` -- a non-negative integer

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.coefficient(9)
            409541017600/2541865828329
            sage: g = f % 5
            sage: g.coefficient(9)
            0
        """
        return self[n]

    def degree(self):
        r"""
        Return the degree of this hypergeometric function if it is
        a polynomial, ``+Infinity`` otherwise.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, -3], [1/2], x)
            sage: f.degree()
            3

        ::

            sage: g = hypergeometric([1/3, 2/3], [1/2], x)
            sage: g.degree()
            +Infinity

        Currently, this method is only implemented in characteristic
        zero::

            sage: T.<y> = GF(5)[]
            sage: h = hypergeometric([1/3, 2/3], [1/2], y)
            sage: h.degree()
            Traceback (most recent call last):
            ...
            NotImplementedError: degree is not implemented in positive characteristic
        """
        if not self._scalar:
            return ZZ(-1)
        if self._char:
            raise NotImplementedError("degree is not implemented in positive characteristic")
        return self._parameters.degree()

    def is_polynomial(self):
        r"""
        Return whether this hypergeometric series is a polynomial.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, -3], [1/2], x)
            sage: f.is_polynomial()
            True

        ::

            sage: g = hypergeometric([1/3, 2/3], [1/2], x)
            sage: g.is_polynomial()
            False
        """
        return self.degree() is not infinity

    def polynomial(self):
        r"""
        Return a polynomial representing a hypergeometric function,
        or raise an error if this hypergeometric function is not
        polynomial.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, -3], [1/2], x)
            sage: f.polynomial()
            -224/405*x^3 + 16/9*x^2 - 2*x + 1

        ::

            sage: g = hypergeometric([1/3, 2/3], [1/2], x)
            sage: g.polynomial()
            Traceback (most recent call last):
            ...
            ValueError: this hypergeometric series is not a polynomial
        """
        deg = self.degree()
        if deg is infinity:
            raise ValueError("this hypergeometric series is not a polynomial")
        S = self.parent().polynomial_ring()
        self._compute_coeffs(deg + 1)
        return S(self._coeffs)

    def power_series(self, prec=20):
        r"""
        Return the power series representation of this hypergeometric
        function up to a given precision.

        INPUT:

        - ``prec`` -- a positive integer (default: ``20``)

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.power_series(3)
            1 + 4/9*x + 80/243*x^2 + O(x^3)
        """
        if prec is infinity:
            S = self.parent().power_series_ring(infinity)
            return S(lambda n: self[n])
        S = self.parent().power_series_ring()
        self._compute_coeffs(prec)
        return S(self._coeffs, prec=prec)

    series = power_series

    def shift(self, s):
        r"""
        Return this hypergeometric function, where each parameter
        (including the additional ``1`` as a bottom parameter) is
        increased by ``s``.

        INPUT:

        - ``s`` -- a rational number

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: g = f.shift(3/2)
            sage: g
            hypergeometric((1, 11/6, 13/6), (2, 5/2), x)
        """
        return self.parent()(self._parameters.shift(s), scalar=self._scalar)

    @coerce_binop
    def hadamard_product(self, other):
        r"""
        Return the Hadamard product of this hypergeometric function
        and ``other``.

        INPUT:

        - ``other`` -- a hypergeometric function

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: h = 1/2*hypergeometric([1/5, 2/5], [3/5], x)
            sage: f.hadamard_product(h)
            1/2*hypergeometric((1/5, 1/3, 2/5, 2/3), (1/2, 3/5, 1), x)
        """
        if self._scalar == 0:
            return self
        if other._scalar == 0:
            return other
        top = self.top() + other.top()
        bottom = self._parameters.bottom + other.bottom()
        scalar = self._scalar * other._scalar
        return self.parent()(top, bottom, scalar=scalar)

    def _div_(self, other):
        r"""
        Return the (formal) quotient of the hypergeometric function
        and ``other``.

        INPUT:

        - ``other`` -- a hypergeometric function or a formal expression

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: g = 1/2 * hypergeometric([1/3, 2/3], [1/2], x)
            sage: h = hypergeometric([1/5, 2/5], [3/5], x)
            sage: f/g
            2
            sage: f/h
            hypergeometric((1/3, 2/3), (1/2,), x)/hypergeometric((1/5, 2/5), (3/5,), x)

        ::

            sage: f/sin(x) + x
            x + hypergeometric((1/3, 2/3), (1/2,), x)/sin(x)
        """
        return SR(self) / SR(other)

    def differential_operator(self, var='d'):
        # Differential equation might not be defined in positive characteristic
        # sage: f = hypergeometric([1/5, 1/5, 1/5], [1/3, 3/5], x)
        # sage: g = f % 3
        # sage: g.differential_operator()
        # Gives error message
        r"""
        Return the hypergeometric differential operator that annihilates
        this hypergeometric function as an Ore polynomial in the variable
        ``var``.

        INPUT:

        - ``var`` -- a string (default: ``d``), the variable name of
          the derivation

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.differential_operator(var='D')
            (-x^2 + x)*D^2 + (-2*x + 1/2)*D - 2/9

        Note that this does not necessarily give the minimal differential
        operator annihilating this hypergeometric function: in the example
        below, this method returns an operator of order `3` where `g` is
        solution of a differential equation of order `2`::

            sage: g = hypergeometric([1/3, 2/3, 6/5], [1/5, 1/2], x)
            sage: L = g.differential_operator()
            sage: L.degree()
            3
            sage: gs = g.power_series(100)
            sage: (72*x^3 - 234*x^2 + 162*x)*gs.derivative(2) + (144*x^2 - 450*x + 81)*gs.derivative() + (16*x - 216)*gs
            O(x^99)
        """
        S = self.parent().polynomial_ring()
        x = S.gen()
        D = OrePolynomialRing(S, S.derivation(), names=var)
        if self._scalar == 0:
            return D.one()
        t = x * D.gen()
        A = D.one()
        for a in self._parameters.top:
            A *= t + S(a)
        B = D.one()
        for b in self._parameters.bottom:
            B *= t + S(b-1)
        L = B - x*A
        return D([c//x for c in L.list()])

    def derivative(self):
        r"""
        Return the derivative of this hypergeometric function.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.derivative()
            4/9*hypergeometric((4/3, 5/3), (3/2,), x)
        """
        top = [a+1 for a in self.top()]
        bottom = [b+1 for b in self.bottom()]
        scalar = prod(self._parameters.top) / prod(self._parameters.bottom)
        scalar = self.base_ring()(scalar) * self._scalar
        return self.parent()(top, bottom, scalar)


# Over the rationals

class HypergeometricAlgebraic_QQ(HypergeometricAlgebraic):
    r"""
    Class for hypergeometric functions over `\QQ`.
    """
    def __mod__(self, p):
        r"""
        Return the reduction of the hypergeometric function modulo ``p``.

        INPUT:

        - ``p`` -- a prime number.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: g = f % 5
            sage: g
            hypergeometric((1/3, 2/3), (1/2,), x)
            sage: g.base_ring()
            Finite Field of size 5

        If the hypergeometric function does not have good reduction at `p`,
        an error is raised::

            sage: f % 3
            Traceback (most recent call last):
            ValueError: the parameters (1/3, 2/3) and (1/2,) do not
            define a hypergeometric function in characteristic 3
        """
        k = FiniteField(p)
        val = self._scalar.valuation(p)
        if val == 0:
            return self.change_ring(k)
        h = self.change_ring(Qp(p, 1))
        return h.residue()

    def valuation(self, p, position=False):
        r"""
        Return the `p`-adic valuation of this hypergeometric function, i.e., the
        maximal `s`, such that `p^{-s}` times this hypergeometric function has
        p-integral coefficients.

        INPUT:

        - ``p`` -- a prime number

        - ``position`` -- a boolean (default: ``False``); if ``True``, return
          also the first index in the series expansion at which the valuation
          is attained.

        ALGORITHM:

        See [CF2026]_, Section 2.2

        EXAMPLES::

            sage: S.<x> = QQ[x]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.valuation(5)
            0
            sage: g = 5*f
            sage: g.valuation(5)
            1

        An example where we ask for the position::

            sage: h = hypergeometric([1/5, 1/5, 1/5], [1/3, 9/5], x)
            sage: h.valuation(3, position=True)
            (-1, 1)

        We can check that the coefficient in `x` in the series expansion
        has indeed valuation `-1`::

            sage: s = h.power_series()
            sage: s
            1 + 1/75*x + 27/8750*x^2 + ... + O(x^20)
            sage: s[1].valuation(3)
            -1

        TESTS::

            sage: g.valuation(9)
            Traceback (most recent call last):
            ...
            ValueError: p must be a prime number
        """
        if not p.is_prime():
            raise ValueError("p must be a prime number")
        val, pos, _ = self._parameters.valuation_position(p)
        val += self._scalar.valuation(p)
        if position:
            return val, pos
        else:
            return val

    def has_good_reduction(self, p):
        r"""
        Return whether the `p`-adic valuation of this hypergeometric
        function is nonnegative, i.e., if its reduction modulo ``p``
        is well-defined.

        INPUT:

        - ``p`` -- a prime number

        EXAMPLES::

            sage: S.<x> = QQ[x]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.valuation(5)
            0
            sage: f.has_good_reduction(5)
            True
            sage: g = 1/5*f
            sage: g.has_good_reduction(5)
            False
        """
        return self.valuation(p) >= 0

    def good_reduction_primes(self):
        r"""
        Return the set of prime numbers modulo which this hypergeometric
        function can be reduced, i.e., the p-adic valuation is nonnegative.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.good_reduction_primes()
            Set of all prime numbers with 3 excluded: 2, 5, 7, 11, ...

        ALGORITHM:

        We implement the algorithm of [CF2026]_, Subsection 3.1

        EXAMPLES::

            sage: f = hypergeometric([1/5, 2/5, 3/5], [1/2, 1/7, 1/11], x)
            sage: f.good_reduction_primes()
            Finite set of prime numbers: 2, 7, 11

        ::

            sage: g = hypergeometric([1/4, 1/2, 3/4], [1/8], x)
            sage: g.good_reduction_primes()
            Set of prime numbers congruent to 3, 5, 7 modulo 8: 3, 5, 7, 11, ...
            sage: (73*g).good_reduction_primes()
            Set of prime numbers congruent to 3, 5, 7 modulo 8 with 73 included: 3, 5, 7, 11, ...
        """
        scalar = self._scalar
        if scalar == 0:
            return Primes()
        params = self._parameters
        d = params.d
        bound = params.bound

        exceptions = {}
        for p in Primes():
            if p > bound:
                break
            val, _, _ = params.valuation_position(p)
            exceptions[p] = (val + scalar.valuation(p) >= 0)

        classes = []
        F = None
        for c in range(bound, bound + d):
            if d.gcd(c) > 1:
                continue
            val, _, _ = params.valuation_position(c)
            if val >= 0:
                classes.append(c % d)
            if val is not -infinity:
                if F is None:
                    F = scalar.factor()
                for p, mult in F:
                    if p > bound and (p-c) % d == 0:
                        exceptions[p] = (val + mult >= 0)

        return Primes(modulus=d, classes=classes, exceptions=exceptions)

    def is_algebraic(self):
        r"""
        Return ``True`` if this hypergeometric function is algebraic over
        the rational functions, return ``False`` otherwise.

        ALGORITHM:

        We rely on the (Christol-)Beukers-Heckmann interlacing criterion
        (see [Chr1986]_, p.15, Cor.; [BeukersHeckman]_, Thm. 4.5). For
        integer differences between parameters we follow the flowchart in
        [FY2024]_, Fig. 1.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.is_algebraic()
            True
            sage: g = hypergeometric([1/3, 2/3, 1/4], [5/4, 1/2], x)
            sage: g.is_algebraic()
            False

        Using `fricas`, we can further compute minimal polynomials::

            sage: fricas.guessAlg(f.power_series(20).list())  # optional - fricas
            [
                 n                    3
              [[x ]f(x): (4 x - 4)f(x)  + 3 f(x) + 1 = 0,
                                    2         3
                          4 x   80 x    1792 x       4
               f(x) = 1 + --- + ----- + ------- + O(x )]
                           9     243      6561
              ]
        """
        parameters = self._parameters.remove_positive_integer_differences()
        if any(a in ZZ and a <= 0 for a in parameters.top):
            return True
        if not parameters.is_balanced():
            return False
        if parameters.has_negative_integer_differences():
            return False
        d = parameters.d
        return all(parameters.interlacing_criterion(c)
                   for c in range(d) if d.gcd(c) == 1)

    def is_globally_bounded(self, include_infinity=True):
        r"""
        Return whether this hypergeometric function is globally bounded
        (if ``include_infinity`` is ``False`` it is not checked whether
        the radius of convergence is finite).

        INPUT:

        - ``include_infinity`` -- a boolean (default: ``True``)

        ALGORITHM:

        We rely on Christol's classification of globally bounded
        hypergeometric functions (see [Chr1986]_, Prop. 1).

        EXAMPLES:

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/9, 4/9, 5/9], [1/3, 1], x)
            sage: f.is_globally_bounded()
            True
            sage: g = hypergeometric([1/9, 4/9, 5/9], [1/3], x)
            sage: g.is_globally_bounded()
            False
            sage: g.is_globally_bounded(include_infinity=False)
            True
        """
        if self.is_polynomial():
            return True
        if include_infinity and len(self.top()) > len(self.bottom()) + 1:
            return False
        d = self._parameters.d
        for c in range(d):
            if d.gcd(c) == 1:
                if not self._parameters.parenthesis_criterion(c):
                    return False
        return True

    def p_curvature_coranks(self):
        r"""
        Return a dictionary, where the integers from `1` to the number of
        parameters of this hypergeometric function are assigned the set of
        prime numbers for which the `p`-curvature has this given corank.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: g = hypergeometric([1/8, 3/8, 1/2], [1/4, 5/8], x)
            sage: g.p_curvature_coranks()
            {1: Set of prime numbers congruent to 3, 5 modulo 8: 3, 5, 11, 13, ...,
             2: Set of prime numbers congruent to 1, 7 modulo 8: 7, 17, 23, 31, ...,
             3: Empty set of prime numbers}
        """
        # Do we have an example with exceptional primes?
        if not self._parameters.is_balanced():
            raise NotImplementedError("Only implemented for nFn-1")
        d = self._parameters.d
        classes = dict.fromkeys(range(1, len(self.top())+1), Primes(modulus=0))
        for c in range(d):
            if gcd(c, d) == 1:
                Delta = QQ(1/c) % d
                j = self._parameters.interlacing_number(Delta)
                classes[j] = classes[j].union(Primes(modulus=d, classes=[c]))
        for p in Primes():
            # I am sure one can avoid computing the interlacing number again for
            # all primes here.
            if p > self._parameters.bound:
                break
            if gcd(p, d) > 1:
                # Do we exclude too many primes here? For which p is the
                # hypergeometric differential equation defined?
                continue
            qinterlacing = self._parameters.q_interlacing_number(p)
            cinterlacing = self._parameters.interlacing_number(QQ(1/p) % d)
            if qinterlacing != cinterlacing:
                classes[qinterlacing].include(p)
                classes[cinterlacing].exclude(p)
        return classes

    def monodromy(self, x=0, var='z'):
        r"""
        Return a local monodromy matrix of the hypergeometric differential
        equation associated to this hypergeometric function at the point
        ``x``.

        INPUT:

        - ``x`` -- a complex number (default: ``0``)

        - ``var`` -- a string (default: ``z``), the name of the variable
          representing a `d`-th root of unity for `d` being the least
          common multiple of the parameters.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.monodromy()
            [0 1]
            [1 0]

        The bases of the solution space are chosen in a compatible way
        across the three singularities of the differential equation::

            sage: g = hypergeometric([1/9, 4/9, 5/9], [1/3, 1], x)
            sage: g.monodromy(var='a')
            [ -a^3 + 1         1         0]
            [2*a^3 + 1         0         1]
            [ -a^3 - 1         0         0]
            sage: g.monodromy(x=Infinity) * g.monodromy(x=1) * g.monodromy()
            [1 0 0]
            [0 1 0]
            [0 0 1]

        ALGORITHM:

        We use the explicit formulas for the monodromy matrices presented in
        [BeukersHeckman]_, Thm. 3.5, attributed to Levelt.
        """
        params = self._parameters
        if not params.is_balanced():
            raise ValueError("hypergeometric equation is not Fuchsian")
        d = params.d
        K = CyclotomicField(d, names=var)
        z = K.gen()
        S = PolynomialRing(K, names='X')
        X = S.gen()
        if x == 0:
            B = prod(X - z**(b*d) for b in params.bottom)
            return companion_matrix(B, format='right').inverse()
        elif x == 1:
            A = prod(X - z**(a*d) for a in params.top)
            B = prod(X - z**(b*d) for b in params.bottom)
            return companion_matrix(A, format='right').inverse() * companion_matrix(B, format='right')
        elif x is infinity:
            A = prod(X - z**(a*d) for a in params.top)
            return companion_matrix(A, format='right')
        else:
            n = len(params.top)
            return identity_matrix(QQ, n)

    def is_maximum_unipotent_monodromy(self):
        r"""
        Return whether the hypergeometric differential operator associated
        to this hypergeometric function has maximal unipotent monodromy (MUM).

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.is_maximum_unipotent_monodromy()
            False
            sage: g = hypergeometric([1/9, 4/9, 5/9], [1, 2], x)
            sage: g.is_maximum_unipotent_monodromy()
            True
        """
        return all(b in ZZ for b in self.bottom())

    is_mum = is_maximum_unipotent_monodromy


# Over the p-adics

class HypergeometricAlgebraic_padic(HypergeometricAlgebraic):
    r"""
    Class for hypergeometric functions over `p`-adic fields.
    """
    def __init__(self, parent, arg1, arg2=None, scalar=None, check=True):
        r"""
        Initialize this hypergeometric function.

        INPUT:

        - ``parent`` -- the parent of this function, which has to be
          defined over the p-adics

        - ``arg1``, ``arg2`` -- arguments defining this hypergeometric
          function, they can be:

          - the top and bottom parameters

          - a hypergeometric function and ``None``

          - an instance of the class
            :class:`sage.functions.hypergeometric_parameters.HypergeometricParameters`
            and ``None``

        - ``scalar`` -- an element in the base ring, the scalar by
          which the hypergeometric function is multiplied

        TESTS::

            sage: S.<x> = Qp(5, 3)[]
            sage: h = hypergeometric((1/2, 1/3), (1,), x)
            sage: type(h)
            <class 'sage.functions.hypergeometric_algebraic.HypergeometricFunctions_with_category.element_class'>
            sage: TestSuite(h).run()
        """
        HypergeometricAlgebraic.__init__(self, parent, arg1, arg2, scalar, check)
        K = self.base_ring()
        self._p = K.prime()
        self._e = K.e()

    def residue(self):
        r"""
        Return the reduction of this hypergeometric function in the residue
        field of the p-adics over which this hypergeometric function is
        defined.

        EXAMPLES::

            sage: S.<x> = Qp(5)[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.parent()
            Hypergeometric functions in x over 5-adic Field with capped relative precision 20
            sage: g = f.residue()
            sage: g.parent()
            Hypergeometric functions in x over Finite Field of size 5
        """
        k = self.base_ring().residue_field()
        valscalar = self._scalar.valuation()
        if valscalar == 0:
            return self.change_ring(k)
        val, pos, _ = self._parameters.valuation_position(self._p)
        val += valscalar
        if val < 0:
            raise ValueError("bad reduction")
        if val > 0:
            H = self.parent().change_ring(k)
            return H.zero()
        raise NotImplementedError("the reduction is not a hypergeometric function")
        # In fact, it is x^s * h[s] * h, with
        # . s is pos
        # . h = self.shift(s)

    def dwork_image(self):
        r"""
        Return the hypergeometric function obtained from this one
        by applying the Dwork map to each of its parameters.

        EXAMPLES::

            sage: S.<x> = Qp(7)[]
            sage: f = hypergeometric([1/4, 1/3, 1/2], [2/5, 3/5, 1], x)
            sage: f.dwork_image()
            hypergeometric((1/3, 1/2, 3/4), (1/5, 4/5, 1), x)
        """
        parameters = self._parameters.dwork_image(self._p)
        return self.parent()(parameters, scalar=self._scalar)

    def _log_radius_of_convergence(self):
        r"""
        Helper function for :meth:`log_radius_of_convergence` and
        :meth:`_truncation_bound`.

        TESTS::

            sage: S.<x> = Qp(5)[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f._log_radius_of_convergence()
            0

        ::

            sage: g = hypergeometric([-3], [1/3], x)
            sage: g._log_radius_of_convergence()
            -1/4
            sage: g.log_radius_of_convergence()
            +Infinity
        """
        p = self._p
        step = self._e / (p - 1)
        log_radius = 0
        parameters = self._parameters
        for a in parameters.top:
            v = a.valuation(p)
            if v < 0:
                log_radius += v
            else:
                log_radius += step
        for b in parameters.bottom:
            v = b.valuation(p)
            if v < 0:
                log_radius -= v
            else:
                log_radius -= step
        return log_radius

    def log_radius_of_convergence(self):
        r"""
        Return the logarithmic `p`-adic radius of convergence of this
        hypergeometric function, that is the exponent on `p` on the
        `p`-adic radius of convergence.

        EXAMPLES::

            sage: S.<x> = Qp(5)[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.log_radius_of_convergence()
            0

        Here the `p`-adic radius of convergence is `p^0 = 1`,
        whereas, in the example below, it is `p^{5/4}`::

            sage: g = hypergeometric([1/3, 2/3], [1/5], x)
            sage: g.log_radius_of_convergence()
            5/4
        """
        if self.is_polynomial():
            return infinity
        return self._log_radius_of_convergence()

    def valuation(self, log_radius=0, position=False):
        r"""
        Return the p-adic valuation of this hypergeometric function on the
        disk of logarithmic radius ``log_radius``, and, if ``position`` is
        ``True`` the index of the first coefficient of the series that
        attains this valuation.

        INPUT:

        - ``log_radius`` -- a rational number

        - ``position`` -- a boolean (default: ``False``), if ``True`` the
          index of the first coefficient attaining the valuation is also
          returned

        ALGORITHM:

        See [CF2026]_, Section 2.2

        EXAMPLES::

            sage: S.<x> = Qp(5)[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.valuation()
            0

        ::

            sage: S.<x> = Qp(5)[]
            sage: g = 1/5 * hypergeometric([1/3, 2/3], [5^3/3], x)
            sage: g.valuation(-1, position=True)
            (-3, 1)
        """
        drift = -log_radius / self._e
        val, pos, _ = self._parameters.valuation_position(self._p, drift)
        val += self._scalar.valuation()
        if position:
            return val, pos
        return val

    def newton_polygon(self, log_radius=None):
        r"""
        Return the Newton polygon of this hypergeometric series.

        INPUT:

        - ``log_radius`` -- a rational number (default: ``None``);
          the last slope of the Newton polygon; if ``None``, the
          logarithmic `p`-adic radius of convergence of this
          hypergeometric function is used.

        ALGORITHM:

        See [CF2026]_, Section 2.3

        EXAMPLES::

            sage: S.<x> = Qp(19)[]
            sage: h = hypergeometric([1/5, 2/5, 3/5, 1/11], [1/2, 1/7], x)
            sage: h.newton_polygon()
            Traceback (most recent call last):
            ...
            ValueError: infinite Newton polygon; try to truncate it by giving a log radius less than 1/18

        Here the Newton polygon has an infinite number of vertices, so it
        cannot be computed entirely.
        As suggested by the error message, we can obtain a result by passing
        in a log radius or last slope: all the segments with slope less than
        this number will be discarded, resulting then in a finite number of
        vertices. When the given log radius gets closer to the actual log
        radius of convergence, the result gets more and more accurate::

            sage: h.newton_polygon(1/18 - 1/10)
            Infinite Newton polygon with 2 vertices: (0, 0), (10, -1) ending by an infinite line of slope -2/45
            sage: h.newton_polygon(1/18 - 1/1000)
            Infinite Newton polygon with 4 vertices: (0, 0), (10, -1), (11, -1), (144, 6) ending by an infinite line of slope 491/9000
        """
        scalar = self._scalar
        if scalar == 0:
            raise ValueError
        convergence = self.log_radius_of_convergence()
        if log_radius is None:
            log_radius = convergence
        start = -log_radius / self._e
        try:
            vertices = self._parameters.newton_polygon(self._p, start)
        except ValueError:
            raise ValueError("infinite Newton polygon; try to truncate it by giving a log radius less than %s" % convergence)
        valscalar = self._scalar.valuation()
        vertices = [[k, v + valscalar] for k, v in vertices]
        return NewtonPolygon(vertices, last_slope=log_radius)

    def _truncation_bound(self, log_radius, prec):
        r"""
        Return a bound of the number of terms needed to evaluate
        this hypergeometric function at precision ``prec`` at a
        `p`-adic argument of valuation ``-log_radius``.

        This is an helper method for :meth:`__call__` and
        :meth:`tate_series`.

        INPUT::

        - ``log_radius`` -- a rational number

        - ``prec`` -- an integer

        TESTS::

            sage: S.<x> = Qp(19)[]
            sage: h = hypergeometric([1/5, 2/5, 3/5, 1/11], [1/2, 1/7], x)
            sage: h._truncation_bound(0, 10)
            232
            sage: h._truncation_bound(0, 20)
            410
            sage: h._truncation_bound(1, 20)
            Traceback (most recent call last):
            ...
            ValueError: outside the domain of convergence
        """
        degree = self.degree()
        convergence = self._log_radius_of_convergence()
        margin = convergence - log_radius
        if margin <= 0:
            if degree is infinity:
                raise ValueError("outside the domain of convergence")
            else:
                return 1 + degree
        val = self.valuation(convergence)
        if val is not -infinity:
            lr = convergence
        else:
            # We choose an intermediate log_radius
            # It can be anything between convergence and log_radius
            # but it seems that the following works well (in the sense
            # that it gives good bounds at the end).
            lr = convergence - margin / max(prec, 2)
            val = self.valuation(lr)
        # Now, we know that
        #   val(h_k) >= -lr*k + val
        # and we want to find k such that
        #   val(h_k) >= -log_radius*k + prec
        # So we just solve the equation.
        k = (prec - val) / (lr - log_radius)
        return 1 + min(degree, max(0, floor(k)))

    def tate_series(self, log_radius, prec=None):
        r"""
        Return this hypergeometric series viewed in the Tate
        algebra with the given log radius.

        INPUT:

        - ``log_radius`` -- a rational number

        - ``prec`` -- a positive integer (default: ``None``);
          if ``None``, use the default precision of the base ring

        EXAMPLES::

            sage: K = Qp(7, prec=5, print_mode='digits')
            sage: S.<x> = K[]
            sage: h = hypergeometric([1/5, 2/5, 3/5, 1/11], [1/2, 1/7], x)
            sage: h.tate_series(0)
            ...00001 + ...40040*x + ...44000*x^2 + ...20000*x^4 + ...30000*x^3 + O(7^5 * <x>)
            sage: h.tate_series(1)
            ...562320000*x^4 + ...140040*x + ...00001 + ...5131000000*x^5 + ... + O(7^5 * <7*x>)

        The given log radius needs to be less than the `p`-adic logarithmic
        radius of convergence of the hypergeometric series. Otherwise, the
        hypergeometric series does not define an element in the corresponding
        Tate algebra and an error is raised::

            sage: h.log_radius_of_convergence()
            4/3
            sage: h.tate_series(2)
            Traceback (most recent call last):
            ...
            ValueError: outside the domain of convergence

        .. SEEALSO::

            :mod:`sage.rings.tate_algebra`
        """
        K = self.base_ring()
        name = self.parent().variable_name()
        S = TateAlgebra(K, log_radii=[log_radius], names=name)
        scalar = self._scalar
        if scalar == 0:
            return S.zero()
        if prec is None:
            prec = self.base_ring().precision_cap()
        trunc = self._truncation_bound(log_radius, prec - scalar.valuation())
        self._compute_coeffs(trunc)
        coeffs = {(i,): self._coeffs[i] for i in range(trunc)}
        return scalar * S(coeffs, prec)

    def __call__(self, x):
        r"""
        Return this hypergeometric function evaluated at ``x``.

        INPUT:

        - ``x`` -- a `p`-adic number

        EXAMPLES::

            sage: K = Qp(7, prec=5)
            sage: S.<x> = K[]
            sage: h = hypergeometric([1/5, 2/5, 3/5, 1/11], [1/2, 1/7], x)
            sage: h(1)
            1 + 4*7 + 4*7^3 + 6*7^4 + O(7^5)
            sage: h(1/7)
            5*7 + 2*7^2 + 3*7^3 + 7^4 + O(7^5)
            sage: h(1/49)
            Traceback (most recent call last):
            ...
            ValueError: outside the domain of convergence
        """
        K = self.base_ring()
        scalar = self._scalar
        if scalar == 0:
            return K.zero()
        x = K(x)
        val = min(x.valuation(), x.precision_absolute())
        if val is infinity:
            return K.one()
        w = self.valuation(-val)
        prec = w + K.precision_cap()
        trunc = self._truncation_bound(-val, prec - scalar.valuation())
        self._compute_coeffs(trunc)
        ans = sum(self._coeffs[i] * x**i for i in range(trunc))
        ans = ans.add_bigoh(prec)
        return scalar * ans


# Over prime finite fields

class HypergeometricAlgebraic_GFp(HypergeometricAlgebraic):
    r"""
    Class for hypergeometric functions over prime finite fields.
    """
    def __init__(self, parent, arg1, arg2=None, scalar=None, check=True):
        r"""
        Initialize this hypergeometric function.

        INPUT:

        - ``parent`` -- the parent of this function, which has to be
          defined over a finite field

        - ``arg1``, ``arg2`` -- arguments defining this hypergeometric
          function, they can be:

          - the top and bottom parameters

          - a hypergeometric function and ``None``

          - an instance of the class
            :class:`sage.functions.hypergeometric_parameters.HypergeometricParameters`
            and ``None``

        - ``scalar`` -- an element in the base ring, the scalar by
          which the hypergeometric function is multiplied

        TESTS::

            sage: S.<x> = GF(5)[]
            sage: h = hypergeometric((1/2, 1/3), (1,), x)
            sage: type(h)
            <class 'sage.functions.hypergeometric_algebraic.HypergeometricFunctions_with_category.element_class'>
            sage: TestSuite(h).run()

        ::

            sage: S.<x> = GF(5)[]
            sage: h = hypergeometric((1/2, 1/3), (1/7,), x)
            Traceback (most recent call last):
            ...
            ValueError: the parameters (1/3, 1/2) and (1/7,) do not define a hypergeometric function in characteristic 5
        """
        HypergeometricAlgebraic.__init__(self, parent, arg1, arg2, scalar, check)
        self._p = p = self.base_ring().cardinality()
        self._coeffs_enriched = [(Qp(p, 1).one(), 0)]

    # def __call__(self, x):
    #     return self.polynomial()(x)

    def __getitem__(self, n):
        r"""
        Return the ``n``-th coefficient of the series representation of this
        hypergeoimetric function.

        INPUT:

        - ``n`` -- a non-negative integer

        EXAMPLES:

            sage: S.<x> = GF(5)[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f[13]
            3

        Over finite fields, an algorithm with complexity `O(\log n)` is implemented.
        It is then safe to call this method with very large integers ``n``::

            sage: f[22204460492503130808472633361816408]  # very fast
            1
        """
        n = ZZ(n)
        p = self._p
        K = Qp(p, 1)
        scalar = self._scalar
        if not scalar:
            return scalar
        H = self.parent().change_ring(K)
        if n < len(self._coeffs) + p*log(n, p):
            self._compute_coeffs(n+1)
            return self._coeffs[n]
        parameters = self._parameters
        ans = K(scalar)
        while n > 0:
            n, r = n.quo_rem(p)
            h = H(parameters)
            ans *= h[r]
            parameters = parameters.shift(r).dwork_image(p)
        return self.base_ring()(ans)

    @coerce_binop
    def is_equal_as_series(self, other):
        r"""
        Return whether ``self`` and ``other`` define the same series.

        INPUT:

        - ``other`` -- an hypergeometric function over the same base

        EXAMPLES::

            sage: S.<x> = GF(13)[]
            sage: f = hypergeometric([1/12, 1/6], [1/3], x)
            sage: g = hypergeometric([1/12, 1/4], [1/2], x)
            sage: f.is_equal_as_series(g)
            True

        ::

            sage: f.power_series(1000)
            1 + 6*x + 6*x^13 + 10*x^14 + 6*x^169 + 10*x^170 + 10*x^182 + 8*x^183 + O(x^1000)
            sage: g.power_series(1000)
            1 + 6*x + 6*x^13 + 10*x^14 + 6*x^169 + 10*x^170 + 10*x^182 + 8*x^183 + O(x^1000)

        We emphasize that, although they define the same series,
        `f` and `g` are not considered as equal::

            sage: f == g
            False
        """
        if self.is_equal_symbolically(other):
            return True
        if self._scalar != other._scalar:
            return False
        H = self.parent()
        p = self._p
        queued = [(self._parameters, other._parameters)]
        checked = {}
        index = 0
        while index < len(queued):
            left, right = queued[index]
            index += 1
            if (left, right) in checked or (right, left) in checked:
                continue
            checked[(left, right)] = True
            criticals = [(1 - pa) % p
                         for pa in left.top + left.bottom + right.top + right.bottom
                         if pa.denominator() % p]
            criticals.sort()
            criticals.append(p)
            for i in range(len(criticals) - 1):
                ei = criticals[i]
                ej = criticals[i+1]
                if ei == ej:
                    continue
                ld = left.shift(ei).dwork_image(p).reduce(p)
                _, lpos, _ = ld.valuation_position(p)
                rd = right.shift(ei).dwork_image(p).reduce(p)
                _, rpos, _ = rd.valuation_position(p)
                if lpos is None or rpos is None:
                    if lpos != rpos:
                        return False
                    continue
                lh = H(left)
                rh = H(right)
                if lh[ei + lpos*p] == 0 and rh[ei + rpos*p] == 0:
                    continue
                if lpos != rpos or any(lh[r + lpos*p] != rh[r + rpos*p] for r in range(ei, ej)):
                    return False
                queued.append((ld.shift(lpos), rd.shift(rpos)))
        return True

    def is_algebraic(self):
        # I am convinced that this is true, but strictly speaking we only have
        # a statement for almost all primes in the literature.
        r"""
        Return whether this hypergeometric function is algebraic.

        This method always returns ``True`` since every hypergeometric
        function in characteristic `p` is algebraic.

        EXAMPLES::

            sage: S.<x> = GF(13)[]
            sage: f = hypergeometric([1/5, 2/5, 3/5, 1/11], [1/2, 1/7], x)
            sage: f.is_algebraic()
            True
        """
        return True

    def p_curvature(self):
        r"""
        Return the matrix of the `p`-curvature of the associated differential
        operator, in the standard basis.

        EXAMPLES::

            sage: S.<x> = GF(5)[]
            sage: f = hypergeometric ([1/9, 4/9, 5/9], [1/3, 1], x)
            sage: f.p_curvature()
            [              0 2/(x^5 + 4*x^4) 1/(x^4 + 4*x^3)]
            [              0               0               0]
            [              0               0               0]

        The following example defines an algebraic function over ``QQ``, thus
        its p-curvature vanishes for almost all of its reductions.::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.is_algebraic()
            True
            sage: g = f % 5
            sage: g.p_curvature()
            [0 0]
            [0 0]
        """
        L = self.differential_operator()
        K = L.base_ring().fraction_field()
        S = OrePolynomialRing(K, L.parent().twisting_derivation().extend_to_fraction_field(), names='d')
        L = S(L.list())
        d = S.gen()
        p = self._char
        rows = []
        n = L.degree()
        for i in range(p, p + n):
            Li = d**i % L
            rows.append([Li[j] for j in range(n)])
        return matrix(rows)

    def p_curvature_corank(self):  # maybe p_curvature_rank is preferable?
        # TODO: check if it is also correct when the parameters are not balanced
        r"""
        Return the corank of the ``p``-curvature matrix.

        ALGORITHM:

        We use [CFV2025]_, Thm. 3.1.17 and the fact that the corank of the
        p-curvature agrees with the number of solutions of the hypergeometric
        differential equation.

        EXAMPLES::

            sage: S.<x> = GF(5)[]
            sage: f = hypergeometric([1/9, 4/9, 5/9], [1/3, 1], x)
            sage: f.p_curvature_corank()
            2
        """
        return self._parameters.q_interlacing_number(self._char)

    def section(self, r):
        r"""
        Return the `r`-th section of this hypergeometric series:
        if this series reads `\sum_n a_n x^n`, it is by definition

        .. MATH:

            \sum_{n=0}^\infty a_{r+pn} x^n

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([7/8, 9/8, 11/8], [3/2, 7/4], x)
            sage: g = f % 5
            sage: g.section(0)
            hypergeometric((3/8, 5/8, 7/8), (1/2, 3/4), x)
            sage: g.section(1)
            2*hypergeometric((3/8, 5/8, 7/8), (1/2, 3/4), x)
            sage: g.section(2)
            hypergeometric((5/8, 7/8, 11/8), (3/4, 3/2), x)

        In certain rare cases, the section is not a scalar multiple of
        an hypergeometric function, by a monomial times a hypergeometric
        function.
        Since there is no support for such functions in SageMath at the
        time being, an error is raised in this case::

            sage: g = f % 3
            sage: g.section(1)
            Traceback (most recent call last):
            ...
            NotImplementedError: the reduction is not a hypergeometric function

        TESTS::

            sage: (0*g).section(1)
            0
            sage: g.section(2)
            0
        """
        scalar = self._scalar
        if scalar == 0:
            return self
        H = self.parent()
        p = self._p
        self._compute_coeffs(r+1)
        hr, z = scalar * self._coeffs_enriched[r]
        if z > 0:
            return H.zero()
        parameters = self._parameters.shift(r).dwork_image(p)
        val, _, _ = parameters.valuation_position(p)
        if val + hr.valuation() > 0:
            return H.zero()
        if val < 0:
            raise NotImplementedError("the reduction is not a hypergeometric function")
        scalar = self._scalar * H.base_ring()(hr)
        return H(parameters, scalar=scalar)

    def dwork_relation(self):
        r"""
        Return a list `(P_1, h_1), ..., (P_s, h_s)` where the
        `P_i` are polynomials and the `h_i` are hypergeometric
        functions such that `P_1 h_1^p + \cdots + P_s h_s^p` is
        equal to ``self``.

        .. NOTE::

            This method is used as a main ingrediant in the
            computation of an annihilating polynomial of ``self``
            (see :meth:`annihilating_ore_polynomial`).

        ALGORITHM:

        See [CF2026]_, Subsection 3.2

        EXAMPLES::

            sage: S.<x> = GF(3)[]
            sage: f = hypergeometric([7/8, 9/8, 11/8], [3/2, 7/4], x)
            sage: f.dwork_relation()
            {hypergeometric((1, 21/8, 25/8, 27/8), (3, 13/4, 7/2), x): 2*x^7,
             hypergeometric((3/8, 5/8, 9/8), (1/2, 5/4), x): 1}
        """
        if self._scalar == 0:
            return {}
        parameters = self._parameters
        p = self._char
        H = self.parent()
        S = H.polynomial_ring()
        criticals = [(1 - pa) % p
                     for pa in parameters.top + parameters.bottom
                     if pa.denominator() % p]
        criticals.sort()
        criticals.append(p)
        Ps = {}
        for i in range(len(criticals) - 1):
            ci = criticals[i]
            cj = criticals[i+1]
            if cj == ci:
                continue
            params = parameters.shift(ci).dwork_image(p)
            _, s, _ = params.valuation_position(p)
            if s is None:
                continue
            ci += s*p
            cj += s*p
            h = H(params.shift(s), check=False)
            self._compute_coeffs(cj + 1)
            P = S(self._coeffs[ci:cj])
            if P:
                Ps[h] = Ps.get(h, 0) + (P << ci)
        return Ps

    def annihilating_ore_polynomial(self, var='Frob'):
        r"""
        Return an Ore polynomaial in the Frobenius morphism, that
        annihilates this hypergeometric function.

        ALGORITHM:

        See [CF2026]_, Subsection 3.3

        INPUT:

        - ``var`` -- a string (default: ``Frob``), name of the variable
          representing the Frobenius morphism.

        EXAMPLES::

            sage: S.<x> = GF(5)[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.annihilating_ore_polynomial()
            (4*x^10 + 2*x^5 + 4)*Frob^2 + (4*x^3 + 4*x^2 + 1)*Frob + x^2
            sage: s = f.power_series(1000)
            sage: (4*x^10 + 2*x^5 + 4)*s^(5^2) + (4*x^3 + 4*x^2 + 1)*s^5 + x^2*s
            O(x^1000)

        There is no guarantee that the returned Ore polynomial is minimal.
        As an illustration, in the next example, the method outputs a Ore
        polynomial of degree `2` while `f` is already solution of a Frobenius
        equation of degree `1`::

            sage: S.<x> = GF(11)[]
            sage: f = hypergeometric([1/10, 5/24], [5/12], x)
            sage: f.annihilating_ore_polynomial()
            (8*x^12 + 6*x^11 + 6*x + 10)*Frob^2 + 1
            sage: s = f.power_series(1000)
            sage: s == (1 + 5*x)*s^11
            True
        """
        parameters = self._parameters
        if not parameters.is_balanced():
            raise NotImplementedError("the hypergeometric function is not a pFq with q = p-1")

        p = self._char
        S = self.parent().polynomial_ring()
        zero = S.zero()
        Frob = S.frobenius_endomorphism()
        Ore = OrePolynomialRing(S, Frob, names=var)

        # We remove the scalar
        if self._scalar == 0:
            return Ore.one()
        self = self.parent()(parameters)

        order = parameters.frobenius_order(p)
        bound = self.p_curvature_corank()

        rows = [{self: S.one()}]
        # If row is the i-th item of rows, we have:
        #   self = sum_g row[g] * g**(p**i)
        q = 1
        while True:
            row = {}
            previous_row = rows[-1]
            for _ in range(order):
                row = {}
                for g, P in previous_row.items():
                    for h, Q in g.dwork_relation().items():
                        # here g = sum(Q * h^p)
                        if h in row:
                            row[h] += P * insert_zeroes(Q, q)
                        else:
                            row[h] = P * insert_zeroes(Q, q)
                previous_row = row
                q *= p  # q = p**i
            rows.append(row)

            i = len(rows)
            Mrows = []
            Mqo = 1
            columns = {}
            for j in range(i-1, max(-1, i-2-bound), -1):
                for col in rows[j]:
                    columns[col] = None
            for j in range(i-1, max(-1, i-2-bound), -1):
                Mrow = []
                for col in columns:
                    Mrow.append(insert_zeroes(rows[j].get(col, zero), Mqo))
                Mrows.append(Mrow)
                Mqo *= p ** order
            M = matrix(S, Mrows)

            ker = kernel(M)
            if ker is not None:
                return insert_zeroes(Ore(ker), order)

    def is_lucas(self):
        r"""
        Return whether this hypergeometric function has the ``p``-Lucas
        property.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/5, 4/5], [1], x)
            sage: g = f % 19
            sage: g.is_lucas()
            True
            sage: h = f % 17
            sage: h.is_lucas()
            False

        ::

            sage: S.<x> = GF(11)[]
            sage: h = hypergeometric([1/10, 5/24], [5/12], x)
            sage: h.is_lucas()
            True
        """
        return all(P.degree() < self._p and self.is_equal_as_series(h)
                   for h, P in self.dwork_relation().items())


# Parent
########

class HypergeometricToSR(Map):
    r"""
    Map from hypergeometric series to symbolic ring
    """
    def _call_(self, h):
        r"""
        Return the symbolic expression representing ``h``.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: h = hypergeometric([1/5, 4/5], [1], x)
            sage: SR(h)  # indirect doctest
            hypergeometric((1/5, 4/5), (1,), x)
        """
        return h.scalar() * hypergeometric(h.top(), h.bottom(), SR.var(h.parent().variable_name()))


class ScalarMultiplication(Action):
    r"""
    Action on hypergeometric series by left multiplication
    by scalars.
    """
    def _act_(self, scalar, h):
        r"""
        Return the product ``scalar * h``.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: h = hypergeometric([1/5, 4/5], [1], x)
            sage: 2*h  # indirect doctest
            2*hypergeometric((1/5, 4/5), (1,), x)
        """
        return h.parent()(h, scalar=scalar)


class HypergeometricFunctions(Parent, UniqueRepresentation):
    r"""
    Hypergeometric functions over a base ring.
    """
    def __classcall__(cls, base, name, symbolic_equality=True):
        r"""
        Normalize parameters and call the init function.

        TESTS::

            sage: S.<x> = QQ[]
            sage: H1 = hypergeometric([], [], x, symbolic_equality=False).parent()
            sage: H2 = hypergeometric([], [], x, symbolic_equality=None).parent()
            sage: H1 is H2
            True
        """
        symbolic_equality = bool(symbolic_equality)
        name = normalize_names(1, name)[0]
        return super().__classcall__(cls, base, name, symbolic_equality)

    def __init__(self, base, name, symbolic_equality, category=None):
        r"""
        Initialize this set of hypergeometric functions.

        INPUT:

        - ``base`` -- the base ring

        - ``name`` -- a string, the name of the variable

        - ``symbolic_equality`` -- a boolean (default: ``True``); if
          ``True``, equality of elements in this parents are checked
          by comparing parameters; if ``False``, it is checked by
          comparing series

        - ``category`` -- a category (default: ``None``)

        .. NOTE::

            The option ``symbolic_equality=False`` is much slower
            and not implemented over all bases.

        TESTS::

            sage: S.<x> = QQ[]
            sage: H = hypergeometric([], [], x).parent()
            sage: TestSuite(H).run()
        """
        self._name = name
        self._latex_name = latex_variable_name(name)
        self._char = char = base.characteristic()
        self._symbolic_equality = symbolic_equality
        if char == 0:
            base = pushout(base, QQ)
        if base in FiniteFields() and base.is_prime_field():
            self.Element = HypergeometricAlgebraic_GFp
        elif base is QQ:
            self.Element = HypergeometricAlgebraic_QQ
        elif isinstance(base, pAdicGeneric):
            self.Element = HypergeometricAlgebraic_padic
        else:
            self.Element = HypergeometricAlgebraic
        if category is None:
            category = Sets()
        Parent.__init__(self, base, category=category)
        self.register_action(ScalarMultiplication(base, self, False, operator.mul))
        self.register_action(ScalarMultiplication(base, self, True, operator.mul))
        if char == 0:
            SR.register_coercion(HypergeometricToSR(self.Hom(SR)))

    def _repr_(self):
        r"""
        Return a string representation of this parent.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: H = hypergeometric([], [], x).parent()
            sage: H  # indirect doctest
            Hypergeometric functions in x over Rational Field
        """
        return "Hypergeometric functions in %s over %s" % (self._name, self._base)

    def _an_element_(self):
        r"""
        Return an element in this parent.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: H = hypergeometric([], [], x).parent()
            sage: H.an_element()  # indirect doctest
            hypergeometric((1,), (), x)
        """
        return self([1], [])

    @cached_method
    def zero(self):
        r"""
        Return the zero function in this parent.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: H = hypergeometric([], [], x).parent()
            sage: H.zero()
            0
        """
        return self(None, scalar=0)

    def _coerce_map_from_(self, other):
        r"""
        Return whether there is a coercion map from ``other``
        to ``self``.

        TESTS::

            sage: S.<x> = QQ[]
            sage: HS = hypergeometric([], [], x).parent()
            sage: T.<x> = RR[]
            sage: HT = hypergeometric([], [], x).parent()
            sage: HS.has_coerce_map_from(HT)  # indirect doctest
            False
            sage: HT.has_coerce_map_from(HS)  # indirect doctest
            True

        ::

            sage: HS2 = hypergeometric([], [], x, symbolic_equality=False).parent()
            sage: HS.has_coerce_map_from(HS2)
            False
            sage: HS2.has_coerce_map_from(HS)
            True
        """
        if (isinstance(other, HypergeometricFunctions)
        and self.base_ring().has_coerce_map_from(other.base_ring())):
            if self._symbolic_equality:
                return True
            else:
                return other._symbolic_equality

    def _pushout_(self, other):
        r"""
        Return a parent in which ``self`` and ``other`` both coerce.

        TESTS::

            sage: from sage.categories.pushout import pushout
            sage: S.<x> = QQ[]
            sage: HS = hypergeometric([], [], x).parent()
            sage: pushout(S, HS)
            Symbolic Ring

        ::

            sage: T.<x> = RR[]
            sage: HT = hypergeometric([], [], x).parent()
            sage: pushout(HT, HS) is HT
            True

        ::

            sage: HS2 = hypergeometric([], [], x, symbolic_equality=False).parent()
            sage: pushout(HS, HS2) is HS2
            True
        """
        if isinstance(other, HypergeometricFunctions) and self._name == other._name:
            base = pushout(self.base_ring(), other.base_ring())
            if base is not None:
                symbolic_equality = self._symbolic_equality and other._symbolic_equality
                return HypergeometricFunctions(base, self._name, symbolic_equality)
        if SR.has_coerce_map_from(other):
            return SR

    def base_ring(self):
        r"""
        Return the base ring over which the hypergeometric functions
        in this parent are defined.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: H = hypergeometric([], [], x).parent()
            sage: H.base_ring()
            Rational Field
        """
        return self._base

    def variable_name(self):
        r"""
        Return the variable name of the hypergeometric functions
        in this parent.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: H = hypergeometric([], [], x).parent()
            sage: H.variable_name()
            'x'
        """
        return self._name

    def latex_variable_name(self):
        r"""
        Return the LaTeX variable name of the hypergeometric functions
        in this parent.

        EXAMPLES::

            sage: S.<xi> = QQ[]
            sage: H = hypergeometric([], [], xi).parent()
            sage: H.latex_variable_name()
            '\\xi'
        """
        return self._latex_name

    def symbolic_equality(self):
        r"""
        Return whether or not the equality in the parent is checked
        symbolically.

        EXAMPLES::

            sage: S.<x> = GF(5)[]
            sage: f = hypergeometric([1/2, 1/3], [1], x)
            sage: f.parent().symbolic_equality()
            True

        ::

            sage: g = hypergeometric([1/2, 1/3], [1], x, symbolic_equality=False)
            sage: g.parent().symbolic_equality()
            False
        """
        return self._symbolic_equality

    def change_ring(self, R):
        r"""
        Return the parent for hypergeometric functions in the same
        variable over the ring ``R``.

        INPUT:

        - ``R`` -- a commutative ring

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: H = hypergeometric([], [], x).parent()
            sage: H
            Hypergeometric functions in x over Rational Field
            sage: H.change_ring(GF(5))
            Hypergeometric functions in x over Finite Field of size 5
        """
        return HypergeometricFunctions(R, self._name, self._symbolic_equality)

    def change_variable_name(self, name):
        r"""
        Return the parent for hypergeometric functions over the same
        ring with variable name ``name``.

        INPUT:

        - ``name`` -- a string

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: H = hypergeometric([], [], x).parent()
            sage: H
            Hypergeometric functions in x over Rational Field
            sage: H.change_variable_name('y')
            Hypergeometric functions in y over Rational Field
        """
        return HypergeometricFunctions(self._base, name, self._symbolic_equality)

    def polynomial_ring(self):
        r"""
        Return the polynomial ring with same variable name
        and same base field as ``self``.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: H = hypergeometric([], [], x).parent()
            sage: H.polynomial_ring() is S
            True
        """
        return PolynomialRing(self.base_ring(), self._name)

    def power_series_ring(self, default_prec=None):
        r"""
        Return the power series ring with same variable name
        and same base field as ``self``.

        INPUT:

        - ``default_prec`` -- a positive integer or ``Infinity``
          (default: ``20``)

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: H = hypergeometric([], [], x).parent()
            sage: H.power_series_ring()
            Power Series Ring in x over Rational Field

        When ``default_prec`` is set to ``Infinity``, a lazy
        power series ring is returned::

            sage: H.power_series_ring(infinity)
            Lazy Taylor Series Ring in x over Rational Field
        """
        if default_prec is infinity:
            return LazyPowerSeriesRing(self.base_ring(), self._name)
        else:
            return PowerSeriesRing(self.base_ring(), self._name, default_prec=default_prec)

# Helper functions
##################

def insert_zeroes(P, n):
    r"""
    Return `P(x^n)`.

    INPUT:

    - ``P`` -- a polynomial in `x`

    - ``n`` -- a positive integer

    EXAMPLES::

        sage: from sage.functions.hypergeometric_algebraic import insert_zeroes
        sage: S.<x> = QQ[]
        sage: insert_zeroes(x + 1, 5)
        x^5 + 1
    """
    cs = P.list()
    coeffs = n * len(cs) * [0]
    for i in range(len(cs)):
        coeffs[n*i] = cs[i]
    return P.parent()(coeffs)


def kernel(M, repeat=2):
    r"""
    Return a generator of the left kernel of the polynomial matrix
    `M`, assuming that the latter has rank at most `1`.

    INPUT:

    - ``repeat`` -- a positive integer (default: ``2``); the number
      of evaluation points we pick to check that the kernel is nonzero

    .. NOTE::

        The implementation is based on Cramer determinants.
        It is however currently faster than
        :meth:`sage.matrix.matrix_polynomial_dense.Matrix_polynomial_dense.minimal_kernel_basis`
        for matrices of small sizes with entries of large degrees.

    EXAMPLES::

        sage: from sage.functions.hypergeometric_algebraic import kernel
        sage: S.<x> = GF(5)[]

    When the kernel is zero, the function returns nothing::

        sage: M = matrix(2, 2, [x, x+1, x+2, x+3])
        sage: kernel(M)

    Otherwise, it returns the smallest generator as a list of polynomials::

        sage: M = matrix(3, 2, [x, x+1, x+2, x^2, x^2+2, x^2+4])
        sage: kernel(M)
        [x^4 + 4*x^3 + x + 2, 4*x^2 + 2*x + 3, 4*x^3 + x^2 + 3*x + 2]
    """
    n = M.nrows()
    m = M.ncols()
    if n > m + 1:
        raise RuntimeError
    if n <= m:
        K = M.base_ring().base_ring()
        for _ in range(repeat):
            a = K.random_element()
            Me = matrix(n, m, [f(a) for f in M.list()])
            if Me.rank() == n:
                return
    for J in Subsets(range(m), n-1):
        MJ = M.matrix_from_columns(J)
        minor = MJ.delete_rows([0]).determinant()
        if minor.is_zero():
            continue
        ker = [minor]
        for i in range(1, n):
            minor = MJ.delete_rows([i]).determinant()
            ker.append((-1)**i * minor)
        Z = matrix(ker) * M
        if not Z.is_zero():
            return
        g = ker[0].leading_coefficient() * gcd(ker)
        ker = [c//g for c in ker]
        return ker
