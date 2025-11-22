r"""
Hypergeometric functions over arbitrary rings

[Tutorial]

AUTHORS:

- Xavier Caruso, Florian Fürnsinn (2025-10): initial version
"""

# ***************************************************************************
#    Copyright (C) 2025 Xavier Caruso <xavier.caruso@normalesup.org>
#                       Florian Fürnsinn <florian.fuernsinn@univie.ac.at>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

import operator

from sage.misc.latex import latex
from sage.misc.latex import latex_variable_name

from sage.misc.misc_c import prod
from sage.misc.functional import log
from sage.functions.other import floor, ceil
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
from sage.rings.tate_algebra import TateAlgebra
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing

from sage.functions.hypergeometric_parameters import HypergeometricParameters


# Do we want to implement polynomial linear combinaison
# of hypergeometric functions?
# Advantages:
#  . reductions mod p of hypergeometric functions have this form in general
#  . many methods can be extended to this context
# Difficulty:
#  . not sure we can handle easily simplifications!

class HypergeometricAlgebraic(Element):
    r"""
    Class for hypergeometric functions over arbitrary base rings.
    """
    def __init__(self, parent, arg1, arg2=None, scalar=None):
        r"""
        Initialize this hypergeometric function.

        INPUT:

        - ``parent`` -- the parent of this function

        - ``arg1``, ``arg2`` -- arguments defining this hypergeometric
          function, they can be:
          - the top and bottom paramters
          - a hypergeometric function and ``None``
          - an instance of the class :class:`HypergeometricParameters` and ``None``

        - ``scalar`` -- an element in the base ring, the scalar by
          which the hypergeometric function is multiplied

        TESTS::

            sage: S.<x> = QQ[]
            sage: h = hypergeometric((1/2, 1/3), (1,), x)
            sage: type(h)
            <class 'sage.functions.hypergeometric_algebraic.HypergeometricFunctions.element_class'>
            sage: TestSuite(h).run()
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
        if scalar:
            if any(b in ZZ and b < 0 for b in parameters.bottom):
                raise ValueError("the parameters %s do not define a hypergeometric function" % parameters)
            if char > 0:
                val, _, _ = parameters.valuation_position(char)
                if val < 0:
                    raise ValueError("the parameters %s do not define a hypergeometric function in characteristic %s" % (parameters, char))
        self._scalar = scalar
        self._parameters = parameters
        self._coeffs = [scalar]
        self._char = char

    def __hash__(self):
        return hash((self.base_ring(), self._parameters, self._scalar))

    def __eq__(self, other):
        return (isinstance(other, HypergeometricAlgebraic)
            and self.base_ring() is other.base_ring()
            and self._parameters == other._parameters
            and self._scalar == other._scalar)

    def _repr_(self):
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
        """
        coeffs = self._coeffs
        start = len(coeffs) - 1
        c = coeffs[-1]
        for i in range(start, prec - 1):
            for a in self._parameters.top:
                c *= a + i
            for b in self._parameters.bottom:
                c /= b + i
            coeffs.append(c)

    def coefficient(self, n):
        r"""
        Return the ``n``-th coefficient of the series representation of this
        hypergeoimetric function.

        INPUT:

        - ``n`` -- a non-negative integer

        EXAMPLES:

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.coefficient(9)
            409541017600/2541865828329
            sage: g = f % 5
            sage: g.coefficient(9)
            0
        """
        self._compute_coeffs(n+1)
        S = self.base_ring()
        return S(self._coeffs[n])

    def power_series(self, prec=20):
        r"""
        Return the power series representation of this hypergeometric
        function up to a given precision.

        INPUT:

        - ``prec`` -- a positive integer

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.power_series(3)
            1 + 4/9*x + 80/243*x^2 + O(x^3)
        """
        S = self.parent().power_series_ring()
        self._compute_coeffs(prec)
        return S(self._coeffs, prec=prec)

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
        Return the hadamard product of this hypergeometric function
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

    def denominator(self):
        r"""
        Return the smallest common denominator of the parameters.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.denominator()
            6
        """
        return self._parameters.d

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

        We rely on Christol's criterion ([Chr1986]_, Prop. 1) for globally
        bounded hypergeometric function, from which a criterion can be deduced
        modulo which primes a hypergeometric function can be reduced
        ([CFV2025]_, Thm. 3.1.3). For small primes `p`, we compute the `p`-adic
        valuation of the hypergeometric function individually.
        """
        params = self._parameters
        d = params.d

        # We check the parenthesis criterion for c=1
        if not params.parenthesis_criterion(1):
            return Primes(modulus=0)

        # We check the parenthesis criterion for other c
        # and derive congruence classes with good reduction
        cs = [c for c in range(d) if d.gcd(c) == 1]
        goods = {c: None for c in cs}
        goods[1] = True
        for c in cs:
            if goods[c] is not None:
                continue
            cc = c
            goods[c] = True
            while cc != 1:
                if goods[cc] is False or not params.parenthesis_criterion(cc):
                    goods[c] = False
                    break
                cc = (cc * c) % d
            if goods[c]:
                cc = c
                while cc != 1:
                    goods[cc] = True
                    cc = (cc * c) % d

        # We treat exceptional primes
        bound = params.bound
        exceptions = {}
        for p in Primes():
            if p > bound:
                break
            if d % p == 0 and self.valuation(p) >= 0:
                exceptions[p] = True
            if d % p == 0 or not goods[p % d]:
                continue
            if self.valuation(p) < 0:
                exceptions[p] = False

        goods = [c for c, v in goods.items() if v]
        return Primes(modulus=d, classes=goods, exceptions=exceptions)

    def is_algebraic(self):
        r"""
        Return ``True`` if this hypergeometric function is algebraic over
        the rational functions, return ``False`` otherwise.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.is_algebraic()
            True
            sage: g = hypergeometric([1/3, 2/3, 1/4], [5/4, 1/2], x)
            sage: g.is_algebraic()
            False

        ALGORITHM:

        We rely on the (Christol-)Beukers-Heckmann interlacing criterion
        (see [Chr1986]_, p.15, Cor.; [BeukersHeckman]_, Thm. 4.5). For integer
        differences between parameters we follow the flowchart in
        [FY2024]_, Fig. 1.
        """
        if any(a in ZZ and a <= 0 for a in self.top()):
            return True
        if not self._parameters.is_balanced():
            return False
        simplified_parameters = self._parameters.remove_positive_integer_differences()
        if simplified_parameters.has_negative_integer_differences():
            return False
        d = simplified_parameters.d
        return all(simplified_parameters.interlacing_criterion(c)
                   for c in range(d) if d.gcd(c) == 1)

    def is_globally_bounded(self, include_infinity=True):
        r"""
        Return whether this hypergeometric function is globally bounded
        (if ``include_infinity`` is ``False`` it is not checked whether
        the radius of convergence is finite).

        INPUT:

        - ``include_infinity`` -- Boolean (default: ``True``)

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

        ALGORITHM:

        We rely on Christol's classification of globally bounded hypergeometric
        functions (see [Chr1986]_, Prop. 1).
        """
        if include_infinity and len(self.top()) > len(self.bottom()) + 1:
            return False
        d = self.denominator()
        for c in range(d):
            if d.gcd(c) == 1:
                if not self._parameters.parenthesis_criterion(c):
                    return False
        return True

    def p_curvature_coranks(self):
        r"""
        Return a dictonary, where the integers from ``1`` to the number of
        parameters of this hypergeometric function, are assigned the set of
        prime numbers for which the ``p``-curvature has this given corank.

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
        d = ZZ(self.denominator())
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
    def __init__(self, parent, arg1, arg2=None, scalar=None):
        r"""
        Initialize this hypergeometric function.

        INPUT:

        - ``parent`` -- the parent of this function, which has to be defined
        over the p-adics

        - ``arg1``, ``arg2`` -- arguments defining this hypergeometric
          function, they can be:

          - the top and bottom paramters

          - a hypergeometric function and ``None``

          - an instance of the class :class:`HypergeometricParameters` and ``None``

        - ``scalar`` -- an element in the base ring, the scalar by
          which the hypergeometric function is multiplied

        TESTS::

            sage: S.<x> = Qp(5, 3)[]
            sage: h = hypergeometric((1/2, 1/3), (1,), x)
            sage: type(h)
            <class 'sage.functions.hypergeometric_algebraic.HypergeometricFunctions.element_class'>
            sage: TestSuite(h).run()
        """
        HypergeometricAlgebraic.__init__(self, parent, arg1, arg2, scalar)
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
        if self._scalar.valuation() == 0:
            return self.change_ring(k)
        val, pos, _ = self._parameters.valuation_position(self._p)
        if val < 0:
            raise ValueError("bad reduction")
        if val > 0:
            H = self.parent().change_ring(k)
            return H(self._parameters, scalar=0)
        raise NotImplementedError("the reduction is not a hypergeometric function")
        # In fact, it is x^s * h[s] * h, with
        # . s = pos
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

    def log_radius_of_convergence(self):
        r"""
        Return the logarithmic p-adic radius of convergence of this
        hypergeometric function.

        EXAMPLES::

            sage: S.<x> = Qp(5)[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.log_radius_of_convergence()
            0
            sage: g = hypergeometric([1/3, 2/3], [1/5], x)
            sage: g.log_radius_of_convergence()
            5/4
        """
        p = self._p
        step = self._e / (p - 1)
        log_radius = 0
        for a in self._parameters.top:
            if a in ZZ and a <= 0:
                return infinity
            v = a.valuation(p)
            if v < 0:
                log_radius += v
            else:
                log_radius += step
        for b in self._parameters.bottom:
            v = b.valuation(p)
            if v < 0:
                log_radius -= v
            else:
                log_radius -= step
        return log_radius

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

        EXAMPLES::

            sage: S.<x> = Qp(5)[]
            sage: f = hypergeometric([1/3, 2/3], [1/2], x)
            sage: f.valuation()
            0

        ::

           sage: S.<x> = Qp(5)[]
           sage: g = 1/5 * hypergeometric([1/3, 2/3], [5^3/3], x)
           sage: g.valuation(-1, position=True)
           (-2, 1)
        """
        drift = -log_radius / self._e
        val, pos, _ = self._parameters.valuation_position(self._p, drift)
        if position:
            return val, pos
        else:
            return val

    def newton_polygon(self, log_radius=None):
        convergence = self.log_radius_of_convergence()
        if log_radius is None:
            log_radius = convergence
        start = -log_radius / self._e
        try:
            val = self._parameters.valuation_function(self._p, start)
        except ValueError:
            raise ValueError("infinite Newton polygon; try to truncate it by giving a log radius less than %s" % convergence)
        return NewtonPolygon(val, last_slope=log_radius)

    def _truncation_bound(self, log_radius, prec):
        convergence = self.log_radius_of_convergence()
        margin = convergence - log_radius
        if margin <= 0:
            raise ValueError("outside the domain of convergence")
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
        return 1 + max(0, floor(k))

    def tate_series(self, log_radius, prec=None):
        K = self.base_ring()
        name = self.parent().variable_name()
        S = TateAlgebra(K, log_radii=[log_radius], names=name)
        if prec is None:
            prec = self.base_ring().precision_cap()
        trunc = self._truncation_bound(log_radius, prec)
        self._compute_coeffs(trunc)
        coeffs = {(i,): self._coeffs[i] for i in range(trunc)}
        return self._scalar * S(coeffs, prec)

    def __call__(self, x):
        K = self.base_ring()
        x = K(x)
        val = min(x.valuation(), x.precision_absolute())
        if val is infinity:
            return K.one()
        w = self.valuation(-val)
        prec = w + K.precision_cap()
        trunc = self._truncation_bound(-val, prec)
        self._compute_coeffs(trunc)
        ans = sum(self._coeffs[i] * x**i for i in range(trunc))
        ans = ans.add_bigoh(prec)
        return self._scalar * ans


# Over prime finite fields

class HypergeometricAlgebraic_GFp(HypergeometricAlgebraic):
    def __init__(self, parent, arg1, arg2=None, scalar=None):
        # TODO: do we want to simplify automatically if the
        # hypergeometric series is a polynomial?
        HypergeometricAlgebraic.__init__(self, parent, arg1, arg2, scalar)
        self._p = p = self.base_ring().cardinality()
        self._coeffs = [Qp(p, 1)(self._scalar)]

    def power_series(self, prec=20):
        S = self.parent().power_series_ring()
        self._compute_coeffs(prec)
        try:
            f = S(self._coeffs, prec=prec)
        except ValueError:
            raise ValueError("denominator appears in the series at the required precision")
        return f

    def is_almost_defined(self):
        p = self._char
        d = self.denominator()
        if d.gcd(p) > 1:
            return False
        u = 1
        if not self._parameters.parenthesis_criterion(u):
            return False
        u = p % d
        while u != 1:
            if not self._parameters.parenthesis_criterion(u):
                return False
            u = p*u % d
        return True

    def is_defined(self):
        p = self._char
        if not self.is_almost_defined():
            return False
        bound = self._parameters.bound
        if bound < p:
            return True
        prec = 1 + p ** ceil(log(self._parameters.bound, p))
        try:
            self.series(prec)
        except ValueError:
            return False
        return True

    def is_defined_conjectural(self):
        p = self._char
        if not self.is_almost_defined():
            return False
        bound = self._parameters.bound
        q = p
        while q <= bound:
            if not self._parameters.q_parenthesis_criterion(q):
                return False
            q *= p
        return True

    def __call__(self, x):
        return self.polynomial()(x)

    def is_polynomial(self):
        raise NotImplementedError

    def degree(self):
        raise NotImplementedError

    def polynomial(self):
        raise NotImplementedError

    def is_algebraic(self):
        return True

    def p_curvature(self):
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
        return self._parameters.q_interlacing_number(self._char)

    def dwork_relation(self):
        r"""
        Return (P1, h1), ..., (Ps, hs) such that

            self = P1*h1^p + ... + Ps*hs^p
        """
        parameters = self._parameters
        if not parameters.is_balanced():
            raise ValueError("the hypergeometric function must be a pFq with q = p-1")
        p = self._char
        H = self.parent()
        F = H.base_ring()
        x = H.polynomial_ring().gen()
        coeffs = self._coeffs
        Ps = {}
        for r in range(p):
            params = parameters.shift(r).dwork_image(p)
            _, s, _ = params.valuation_position(p)
            h = H(params.shift(s))
            e = s*p + r
            if e >= len(coeffs):
                self._compute_coeffs(e + 1)
            c = F(coeffs[e])
            if c:
                if h in Ps:
                    Ps[h] += c * x**e
                else:
                    Ps[h] = c * x**e
        return Ps

    def annihilating_ore_polynomial(self, var='Frob'):
        # QUESTION: does this method actually return the
        # minimal Ore polynomial annihilating self?
        # Probably not :-(
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
        p = self._char
        if self._parameters.frobenius_order(p) > 1:
            # TODO: check this
            return False
        S = self.parent().polynomial_ring()
        K = S.fraction_field()
        Ore = OrePolynomialRing(K, K.frobenius_endomorphism(), names='F')
        Z = Ore(self.annihilating_ore_polynomial())
        Ap = self.series(p).polynomial()
        F = Ap * Ore.gen() - 1
        return (Z % F).is_zero()


# Parent
########

class HypergeometricToSR(Map):
    def _call_(self, h):
        return h.scalar() * hypergeometric(h.top(), h.bottom(), SR.var(h.parent().variable_name()))


class ScalarMultiplication(Action):
    def _act_(self, scalar, h):
        return h.parent()(h, scalar=scalar)


class HypergeometricFunctions(Parent, UniqueRepresentation):
    def __init__(self, base, name, category=None):
        self._name = normalize_names(1, name)[0]
        self._latex_name = latex_variable_name(self._name)
        self._char = char = base.characteristic()
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
        Parent.__init__(self, base, category=category)
        self.register_action(ScalarMultiplication(base, self, False, operator.mul))
        self.register_action(ScalarMultiplication(base, self, True, operator.mul))
        if char == 0:
            SR.register_coercion(HypergeometricToSR(self.Hom(SR)))

    def _repr_(self):
        return "Hypergeometric functions in %s over %s" % (self._name, self._base)

    def _element_constructor_(self, *args, **kwds):
        return self.element_class(self, *args, **kwds)

    def _coerce_map_from_(self, other):
        if (isinstance(other, HypergeometricFunctions)
        and other.has_coerce_map_from(self)):
            return True

    def _pushout_(self, other):
        if isinstance(other, HypergeometricFunctions) and self._name == other._name:
            base = pushout(self.base_ring(), other.base_ring())
            if base is not None:
                return HypergeometricFunctions(base, self._name)
        if SR.has_coerce_map_from(other):
            return SR

    def base_ring(self):
        return self._base

    def variable_name(self):
        return self._name

    def latex_variable_name(self):
        return self._latex_name

    def change_ring(self, R):
        return HypergeometricFunctions(R, self._name)

    def change_variable_name(self, name):
        return HypergeometricFunctions(self._base, name)

    def polynomial_ring(self):
        return PolynomialRing(self.base_ring(), self._name)

    def power_series_ring(self, default_prec=None):
        return PowerSeriesRing(self.base_ring(), self._name, default_prec=default_prec)


# Helper functions
##################

def insert_zeroes(P, n):
    cs = P.list()
    coeffs = n * len(cs) * [0]
    for i in range(len(cs)):
        coeffs[n*i] = cs[i]
    return P.parent()(coeffs)


def kernel(M, repeat=2):
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
