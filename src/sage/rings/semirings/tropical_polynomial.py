r"""
Univariate Tropical Polynomials

AUTHORS:

- Verrel Rievaldo Wijaya (2024-06): initial version

EXAMPLES::

    sage: T = TropicalSemiring(QQ, use_min=False)
    sage: R.<x> = PolynomialRing(T)
    sage: x.parent()
    Univariate Tropical Polynomial Semiring in x over Rational Field
    sage: (x + R(3)*x) * (x^2 + x)
    3*x^3 + 3*x^2
    sage: (x^2 + R(1)*x + R(-1))^2
    0*x^4 + 1*x^3 + 2*x^2 + 0*x + (-2)

REFERENCES:

- [Bru2014]_
- [Fil2017]_
"""

# ****************************************************************************
#       Copyright (C) 2024 Verrel Rievaldo Wijaya <verrelrievaldo@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.rings.polynomial.polynomial_element_generic import Polynomial_generic_sparse


class TropicalPolynomial(Polynomial_generic_sparse):
    r"""
    A univariate tropical polynomial.

    A tropical polynomial is a polynomial whose coefficients come from
    a tropical semiring. Tropical polynomial induces a function which is
    piecewise linear, where each segment of the function has a non-negative
    integer slope. Tropical roots (zeros) of polynomial `P(x)` is defined
    as all points `x_0` for which the graph of `P(x)` change its slope.
    The difference in the slopes of the two pieces adjacent to this root
    gives the order of the root.

    The tropical polynomials are implemented with a sparse format by using
    a ``dict`` whose keys are the exponent and values the corresponding
    coefficients. Each coefficient is a tropical number.

    EXAMPLES:

    First, we construct a tropical polynomial semiring by defining a base
    tropical semiring and then inputting it to :class:`PolynomialRing`::

        sage: T = TropicalSemiring(QQ, use_min=False)
        sage: R.<x> = PolynomialRing(T); R
        Univariate Tropical Polynomial Semiring in x over Rational Field
        sage: R.0
        0*x

    One way to construct an element is to provide a list or tuple of
    coefficients. Another way to define an element is to write a polynomial
    equation with each coefficient converted to the semiring::

        sage: p1 = R([1,4,None,0]); p1
        0*x^3 + 4*x + 1
        sage: p2 = R(1)*x^2 + R(2)*x + R(3); p2
        1*x^2 + 2*x + 3

    We can do some basic arithmetic operations for these tropical polynomials.
    Remember that numbers given have to be in the tropical semiring.
    If not, then it will raise an error::

        sage: p1 + p2
        0*x^3 + 1*x^2 + 4*x + 3
        sage: p1 * p2
        1*x^5 + 2*x^4 + 5*x^3 + 6*x^2 + 7*x + 4
        sage: 2 * p1
        Traceback (most recent call last):
        ...
        ArithmeticError: cannot negate any non-infinite element
        sage: T(2) * p1
        2*x^3 + 6*x + 3
        sage: p1(3)
        Traceback (most recent call last):
        ...
        TypeError: no common canonical parent for objects with parents:
        'Tropical semiring over Rational Field' and 'Integer Ring'
        sage: p1(T(3))
        9

    We can find all tropical roots of a tropical polynomial, counted
    with their multiplicities::

        sage: p1.roots()
        [-3, 2, 2]
        sage: p2.roots()
        [1, 1]

    Even though every tropical polynomials have tropical roots, this does not
    neccessarily means it can be factored into its linear factors::

        sage: p1.factor()
        (0) * (0*x^3 + 4*x + 1)
        sage: p2.factor()
        (1) * (0*x + 1)^2

    Every tropical polynomial `p(x)` have a corresponding unique tropical
    polynomial `\bar{p}(x)` with the same roots that can be factored. We
    call `\bar{p}(x)` the tropical polynomial split form of `p(x)`::

        sage: p1.split_form()
        0*x^3 + 2*x^2 + 4*x + 1
        sage: p2.split_form()
        1*x^2 + 2*x + 3

    Every tropical polynomial induce a piecewise linear function that can be
    invoked in the following way::

        sage: p1.piecewise_function()
        piecewise(x|-->1 on (-oo, -3], x|-->x + 4 on (-3, 2), x|-->3*x on [2, +oo); x)
        sage: p1.plot()
        Graphics object consisting of 1 graphics primitive

    .. PLOT::
        :width: 300 px

        T = TropicalSemiring(QQ, use_min=False)
        R = PolynomialRing(T, 'x')
        p1 = R([1,4,None,0])
        sphinx_plot(p1.plot())

    ::

        sage: p2.piecewise_function()
        piecewise(x|-->3 on (-oo, 1], x|-->2*x + 1 on (1, +oo); x)
        sage: p2.plot()
        Graphics object consisting of 1 graphics primitive

    .. PLOT::
        :width: 300 px

        T = TropicalSemiring(QQ, use_min=False)
        R = PolynomialRing(T, 'x')
        x = R.gen()
        p2 = R(1)*x**2 + R(2)*x + R(3)
        sphinx_plot(plot(p2, xmin=-1, xmax=3))

    TESTS:

    There is no subtraction for tropical polynomials because element in
    tropical semiring doesn't necessarily have additive inverse::

        sage: -p1
        Traceback (most recent call last):
        ...
        ArithmeticError: cannot negate any non-infinite element
    """
    def roots(self):
        r"""
        Return the list of all tropical roots of ``self``, counted with
        multiplicity.

        OUTPUT: a list of tropical numbers

        ALGORITHM:

        For each pair of monomials in the polynomial, we find the point
        where their values are equal.  This is the same as solving the
        equation `c_1 + a_1 \cdot x = c_2 + a_2 \cdot x` for `x`, where
        `(c_1, a_1)` and `(c_2, a_2)` are the coefficients and exponents
        of monomials.

        The solution to this equation is `x = (c_2-c_1)/(a_1-a_2)`. We
        substitute this `x` to each monomials in polynomial and check if
        the maximum (minimum) is achieved by the previous two monomials.
        If it is, then `x` is the root of tropical polynomial. In this
        case, the order of the root at `x` is the maximum of `|i-j|` for
        all possible pairs `i,j` which realise this maximum (minimum).

        EXAMPLES::

            sage: T = TropicalSemiring(QQ, use_min=True)
            sage: R.<x> = PolynomialRing(T)
            sage: p1 = R([5,4,1,0,2,4,3]); p1
            3*x^6 + 4*x^5 + 2*x^4 + 0*x^3 + 1*x^2 + 4*x + 5
            sage: p1.roots()
            [-1, -1, -1, 1, 2, 2]

        There will be no tropical root for constant polynomial. Additionaly,
        for a monomial, the tropical root is assumed to be the additive
        identity of its base tropical semiring::

            sage: p2 = R(2)
            sage: p2.roots()
            []
            sage: p3 = x^3
            sage: p3.roots()
            [+infinity, +infinity, +infinity]
        """
        from itertools import combinations
        tropical_roots = []
        data = self.monomial_coefficients()
        R = self.parent().base()
        if len(data) == 1:
            exponent = next(iter(data))
            if not exponent:
                return tropical_roots
            return [R.zero()] * exponent
        dict_root = {}
        dict_coeff = {i: c.lift() for i, c in data.items()}
        for comb in combinations(dict_coeff, 2):
            index1, index2 = comb[0], comb[1]
            root = (dict_coeff[index1]-dict_coeff[index2]) / (index2-index1)
            val_root = dict_coeff[index1] + index1*root
            check_maks = True
            for key in dict_coeff:
                if key not in comb:
                    val = dict_coeff[key] + key*root
                    if R._use_min:
                        if val < val_root:
                            check_maks = False
                            break
                    elif val > val_root:
                        check_maks = False
                        break
            if check_maks:
                order = abs(index1-index2)
                if root not in dict_root:
                    dict_root[root] = order
                elif order > dict_root[root]:
                    dict_root[root] = order
        for root in dict_root:
            tropical_roots += [root] * dict_root[root]
        return sorted(tropical_roots)

    def split_form(self):
        r"""
        Return the tropical polynomial which has the same roots as ``self``
        but which can be reduced to its linear factors.

        If a tropical polynomial has roots at `x_1, x_2, \ldots, x_n`, then
        its split form is the tropical product of linear terms of the form
        `(x + x_i)` for all `i=1,2,\ldots,n`.

        OUTPUT: :class:`TropicalPolynomial`

        EXAMPLES::

            sage: T = TropicalSemiring(QQ, use_min=True)
            sage: R.<x> = PolynomialRing(T)
            sage: p1 = R([5,4,1,0,2,4,3]); p1
            3*x^6 + 4*x^5 + 2*x^4 + 0*x^3 + 1*x^2 + 4*x + 5
            sage: p1.split_form()
            3*x^6 + 2*x^5 + 1*x^4 + 0*x^3 + 1*x^2 + 3*x + 5

        ::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R.<x> = PolynomialRing(T)
            sage: p1 = R([5,4,1,0,2,4,3])
            sage: p1.split_form()
            3*x^6 + 4*x^5 + 21/5*x^4 + 22/5*x^3 + 23/5*x^2 + 24/5*x + 5

        TESTS:

        We verify that the roots of tropical polynomial and its split form
        is really the same::

            sage: p1.roots() == p1.split_form().roots()
            True
        """
        roots = self.roots()
        R = self.parent()
        poly = R(self.monomial_coefficients()[self.degree()].lift())
        for root in roots:
            linear = R([root, 0])
            poly *= linear
        return poly

    def factor(self):
        r"""
        Return the factorization of ``self`` into its tropical linear factors.

        Note that the factor `x - x_0` in classical algebra gets transformed
        to the factor `x + x_0`, since the root of the tropical polynomial
        `x + x_0` is `x_0` and not `-x_0`. However, not every tropical
        polynomial can be factored.

        OUTPUT: a :class:'Factorization'

        EXAMPLES::

            sage: T = TropicalSemiring(QQ, use_min=True)
            sage: R.<x> = PolynomialRing(T)
            sage: p1 = R([6,3,1,0]); p1
            0*x^3 + 1*x^2 + 3*x + 6
            sage: factor(p1)
            (0) * (0*x + 1) * (0*x + 2) * (0*x + 3)

        Such factorization is not always possible::

            sage: p2 = R([4,4,2]); p2
            2*x^2 + 4*x + 4
            sage: p2.factor()
            (2) * (0*x^2 + 2*x + 2)

        TESTS:

        The factorization for a constant::

            sage: p3 = R(3)
            sage: p3.factor()
            (3) * 0
        """
        from sage.structure.factorization import Factorization
        unit = self.monomial_coefficients()[self.degree()]
        if self != self.split_form() or not self.roots():
            factor = [(self * self.parent(-unit.lift()), 1)]
            return Factorization(factor, unit=unit)

        R = self.parent()
        roots_order = {}
        for root in self.roots():
            if root in roots_order:
                roots_order[root] += 1
            else:
                roots_order[root] = 1
        factors = [(R([root, 0]), roots_order[root])
                   for root in roots_order]
        return Factorization(factors, unit=unit)

    def piecewise_function(self):
        r"""
        Return the tropical polynomial function of ``self``.

        The function is a piecewise linear function with the domains are
        divided by the roots. First we convert each term of polynomial to
        its corresponding linear function. Next, we must determine which
        term achieves the minimum (maximum) at each interval.

        OUTPUT: A piecewise function

        EXAMPLES::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R.<x> = PolynomialRing(T)
            sage: p1 = R([4,2,1,3]); p1
            3*x^3 + 1*x^2 + 2*x + 4
            sage: p1.piecewise_function()
            piecewise(x|-->4 on (-oo, 1/3], x|-->3*x + 3 on (1/3, +oo); x)

        A constant tropical polynomial will result in a constant function::

            sage: p2 = R(3)
            sage: p2.piecewise_function()
            3

        A monomial will resulted in a linear function::

            sage: p3 = R(1)*x^3
            sage: p3.piecewise_function()
            3*x + 1
        """
        from sage.symbolic.ring import SR
        from sage.functions.piecewise import piecewise
        from sage.sets.real_set import RealSet

        x = SR.var('x')
        data = self.monomial_coefficients()
        R = self.parent().base()
        if not self.roots():
            f = data[0].lift()
            return f

        if len(data) == 1:
            gradient = list(data)[0]
            intercept = data[gradient].lift()
            f = intercept + gradient*x
            return f

        unique_root = sorted(set(self.roots()))
        pieces = []
        domain = []
        for i in range(len(unique_root) + 1):
            if i == 0:
                test_number = R(unique_root[i] - 1)
            elif i == len(unique_root):
                test_number = R(unique_root[i-1] + 1)
            else:
                test_number = R((unique_root[i]+unique_root[i-1]) / 2)
            terms = {i: c * test_number**i for i, c in data.items()}
            if R._use_min:
                critical = min(terms.values())
            else:
                critical = max(terms.values())
            found_key = None
            for key, value in terms.items():
                if value == critical:
                    found_key = key
                    break
            gradient = found_key
            intercept = data[found_key].lift()

            # To make sure all roots is included in the domain
            if i == 0:
                interval = RealSet.unbounded_below_closed(unique_root[i])
                piecewise_linear = (interval, intercept + gradient*x)
                domain.append(interval)
            elif i == len(unique_root):
                if domain[i-1][0].upper_closed():
                    interval = RealSet.unbounded_above_open(unique_root[i-1])
                else:
                    interval = RealSet.unbounded_above_closed(unique_root[i-1])
                piecewise_linear = (interval, intercept + gradient*x)
                domain.append(interval)
            else:
                if domain[i-1][0].upper_closed():
                    interval = RealSet((unique_root[i-1], unique_root[i]))
                else:
                    interval = RealSet([unique_root[i-1], unique_root[i]])
                piecewise_linear = (interval, intercept + gradient*x)
                domain.append(interval)
            pieces.append(piecewise_linear)

        f = piecewise(pieces)
        return f

    def plot(self, xmin=None, xmax=None):
        r"""
        Return the plot of ``self``, which is the tropical polynomial
        function we get from ``self.piecewise_function()``.

        INPUT:

        - ``xmin`` -- (optional) real number
        - ``xmax`` -- (optional) real number

        OUTPUT:

        If ``xmin`` and ``xmax`` is given, then it return a plot of
        piecewise linear function of ``self`` with the axes start from
        ``xmin`` to ``xmax``. Otherwise, the domain will start from the
        the minimum root of ``self`` minus 1 to maximum root of ``self``
        plus 1. If the function of ``self`` is constant or linear, then
        the default domain will be `[-1,1]`.

        EXAMPLES:

        If the tropical semiring use a max-plus algebra, then the graph
        will be of piecewise linear convex function::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R.<x> = PolynomialRing(T)
            sage: p1 = R([4,2,1,3]); p1
            3*x^3 + 1*x^2 + 2*x + 4
            sage: p1.roots()
            [1/3, 1/3, 1/3]
            sage: p1.plot()
            Graphics object consisting of 1 graphics primitive

        .. PLOT::
            :width: 300 px

            T = TropicalSemiring(QQ, use_min=False)
            R = PolynomialRing(T, 'x')
            p1 = p1 = R([4,2,1,3])
            sphinx_plot(p1.plot())

        A different result will be obtained if the tropical semiring employs
        a min-plus algebra. Rather, a graph of the piecewise linear concave
        function will be obtained::

            sage: T = TropicalSemiring(QQ, use_min=True)
            sage: R.<x> = PolynomialRing(T)
            sage: p1 = R([4,2,1,3])
            sage: p1.roots()
            [-2, 1, 2]
            sage: p1.plot()
            Graphics object consisting of 1 graphics primitive

        .. PLOT::
            :width: 300 px

            T = TropicalSemiring(QQ, use_min=True)
            R = PolynomialRing(T, 'x')
            p1 = R([4,2,1,3])
            sphinx_plot(plot(p1, xmin=-4, xmax=4))

        TESTS:

        If ``xmin`` or ``xmax`` is given as an input, then the others also
        have to be given. Otherwise it will raise an error::

            sage: plot(p1, 5)
            Traceback (most recent call last):
            ...
            ValueError: expected 2 inputs for xmin and xmax, but got 1

        Error also occured when ``xmin`` is greater or equal than``xmax``::

            sage: plot(p1, 5, 3)
            Traceback (most recent call last):
            ...
            ValueError: xmin = 5 should be less than xmax = 3
        """
        from sage.plot.plot import plot
        f = self.piecewise_function()
        if (xmin is None) and (xmax is None):
            roots = sorted(self.roots())
            if (not roots) or (self.parent().base().zero() in roots):
                return plot(f, xmin=-1, xmax=1)
            else:
                return plot(f, xmin=roots[0]-1, xmax=roots[-1]+1)
        elif xmin is None or xmax is None:
            raise ValueError("expected 2 inputs for xmin and xmax, but got 1")
        elif xmin >= xmax:
            raise ValueError(f"xmin = {xmin} should be less than xmax = {xmax}")
        return plot(f, xmin=xmin, xmax=xmax)

    def _repr_(self):
        r"""
        Return a string represemtation of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x> = PolynomialRing(T)
            sage: R([-3,-1,2,-1])
            (-1)*x^3 + 2*x^2 + (-1)*x + (-3)
        """
        import re
        if not self.monomial_coefficients():
            return str(self.parent().base().zero())

        def replace_negatives(expr):
            modified_expr = re.sub(r'(-\d+/\d+|-\d+)', r'(\1)', expr)
            return modified_expr

        s = super()._repr()
        v = self.parent().variable_name()
        if s[0] == v:
            s = "1*" + s
        s = s.replace(" - ", " + -")
        s = s.replace(" + "+v, " + 1*"+v)
        s = s.replace("-"+v, "-1*"+v)
        s = replace_negatives(s)
        return s

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R = PolynomialRing(T, 'x')
            sage: p1 = R([-1,2,None,-3])
            sage: latex(p1)
            \left(-3\right) x^{3} + 2 x + \left(-1\right)
        """
        s = " "
        coeffs = self.list(copy=False)
        m = len(coeffs)
        name = self.parent().latex_variable_names()[0]
        for n in reversed(range(m)):
            x = coeffs[n]
            x = x._latex_()
            if x != self.parent().base().zero()._latex_():
                if n != m-1:
                    s += " + "
                if x.find("-") == 0:
                    x = "\\left(" + x + "\\right)"
                if n > 1:
                    v = f"|{name}^{{{n}}}"
                elif n == 1:
                    v = f"|{name}"
                else:
                    v = ""
                s += f"{x} {v}"
        s = s.replace("|", "")
        if s == " ":
            return self.parent().base().zero()._latex_()
        return s[1:].lstrip().rstrip()


class TropicalPolynomialSemiring(UniqueRepresentation, Parent):
    r"""
    The semiring of univariate tropical polynomials.

    This is the commutative semiring consisting of finite linear
    combinations of tropical monomials under (tropical)
    addition and multiplication with the identity element
    being `+\infty` (or `-\infty` depending on whether the
    base tropical semiring is using min-plus or max-plus algebra).

    EXAMPLES::

        sage: T = TropicalSemiring(QQ)
        sage: R.<x> = PolynomialRing(T)
        sage: f = T(1)*x
        sage: g = T(2)*x
        sage: f + g
        1*x
        sage: f * g
        3*x^2
        sage: f + R.zero() == f
        True
        sage: f * R.zero() == R.zero()
        True
        sage: f * R.one() == f
        True
    """
    @staticmethod
    def __classcall_private__(cls, base_semiring, names):
        """
        Ensures the names parameter is a tuple.

        EXAMPLES::

            sage: from sage.rings.semirings.tropical_polynomial import TropicalPolynomialSemiring
            sage: T = TropicalSemiring(ZZ)
            sage: TPS = TropicalPolynomialSemiring
            sage: TPS(T, 'x') == TPS(T, ('x'))
            True
        """
        if isinstance(names, str):
            names = (names,)
        return super().__classcall__(cls, base_semiring, tuple(names))

    def __init__(self, base_semiring, names):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R = PolynomialRing(T, 'x')
            sage: TestSuite(R).run()

        TESTS::

            sage: from sage.rings.semirings.tropical_polynomial import TropicalPolynomialSemiring
            sage: TropicalPolynomialSemiring(ZZ, names='x')
            Traceback (most recent call last):
            ...
            ValueError: Integer Ring is not a tropical semiring
        """
        from sage.categories.semirings import Semirings
        from sage.rings.semirings.tropical_semiring import TropicalSemiring
        if not isinstance(base_semiring, TropicalSemiring):
            raise ValueError(f"{base_semiring} is not a tropical semiring")
        Parent.__init__(self, base=base_semiring, names=names, category=Semirings())

    Element = TropicalPolynomial

    def _element_constructor_(self, x=None, check=True):
        """
        Convert ``x`` into ``self``, possibly non-canonically.

        INPUT:

        - ``x`` -- a list or tuple of coefficients, a polynomial, or a
          dictionary

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R = PolynomialRing(T, 'x')
            sage: R([1,2,3])
            3*x^2 + 2*x + 1
            sage: S.<x> = PolynomialRing(QQ)
            sage: R(x^2 - x + 1)
            1*x^2 + (-1)*x + 1

        If ``x`` is a tropical polynomial from different semiring, then it
        will converted to constant::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x> = PolynomialRing(T)
            sage: S.<y> = PolynomialRing(T)
            sage: p1 = R(y); p1
            0*y
            sage: p1.parent()
            Univariate Tropical Polynomial Semiring in x over Rational Field
            sage: p1.degree()
            0
        """
        if isinstance(x, (list, tuple)):
            for i, coeff in enumerate(x):
                if coeff == 0:
                    x[i] = self.base().one()
        elif isinstance(x, TropicalPolynomial):
            if x.parent() is not self:
                x = {0: x}
                check = False
        return self.element_class(self, x, check=check)

    @cached_method
    def one(self):
        """
        Return the multiplicative identity of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R = PolynomialRing(T, 'x')
            sage: R.one()
            0
        """
        return self.element_class(self, [self.base().one()])

    @cached_method
    def zero(self):
        """
        Return the additive identity of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R = PolynomialRing(T, 'x')
            sage: R.zero()
            +infinity
        """
        return self.element_class(self, [self.base().zero()])

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(ZZ)
            sage: R.<abc> = PolynomialRing(T); R
            Univariate Tropical Polynomial Semiring in abc over Integer Ring
        """
        return (f"Univariate Tropical Polynomial Semiring in {self.variable_name()}"
                f" over {self.base_ring().base_ring()}")

    def gen(self, n=0):
        """
        Return the indeterminate generator of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<abc> = PolynomialRing(T)
            sage: R.gen()
            0*abc

        TESTS::

            sage: R.gen(2)
            Traceback (most recent call last):
            ...
            IndexError: generator n not defined
        """
        if n != 0:
            raise IndexError("generator n not defined")
        return self.gens()[n]

    @cached_method
    def gens(self):
        """
        Return a tuple whose entries are the generators for ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<abc> = PolynomialRing(T)
            sage: R.gens()
            (0*abc,)
        """
        R = self.base()
        return (self.element_class(self, [R.zero(), R.one()]),)

    def ngens(self):
        """
        Return the number of generators of ``self``, which is 1
        since it is a univariate polynomial ring.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R = PolynomialRing(T, 'x')
            sage: R.ngens()
            1
        """
        from sage.rings.integer_ring import ZZ
        return ZZ.one()

    def random_element(self, degree=(-1, 2), monic=False, *args, **kwds):
        r"""
        Return a random tropical polynomial of given degrees (bounds).

        OUTPUT: a :class:`TropicalPolynomial`

        .. SEEALSO::

            :meth:`sage.rings.polynomial.polynomial_ring.PolynomialRing_general.random_element`

        EXAMPLES:

        Tropical polynomial over an integer with each coefficient bounded
        between ``x`` and ``y``::

            sage: T = TropicalSemiring(ZZ)
            sage: R = PolynomialRing(T, 'x')
            sage: f = R.random_element(8, x=3, y=10)
            sage: f.degree()
            8
            sage: f.parent() is R
            True
            sage: all(a >= T(3) for a in f.coefficients())
            True
            sage: all(a < T(10) for a in f.coefficients())
            True

        If a tuple of two integers is provided for the ``degree`` argument,
        a polynomial is selected with degrees within that range::

            sage: all(R.random_element(degree=(1, 5)).degree() in range(1, 6) for _ in range(10^3))
            True

        Note that the zero polynomial (`\pm \infty`) has degree `-1`.
        To include it, set the minimum degree to `-1`::

            sage: while R.random_element(degree=(-1,2), x=-1, y=1) != R.zero():
            ....:     pass

        Monic polynomials are chosen among all monic polynomials with degree
        between the given ``degree`` argument::

            sage: all(R.random_element(degree=(-1, 2), monic=True).is_monic() for _ in range(10^3))
            True
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        R = PolynomialRing(self.base().base_ring(), self.variable_names())
        f = R.random_element(degree=degree, monic=monic, *args, **kwds)
        new_dict = f.monomial_coefficients()
        if monic:
            new_dict[f.degree()] = 0
        return self.element_class(self, new_dict)

    def is_sparse(self):
        """
        Return ``True`` to indicate that the objects are sparse polynomials.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R = PolynomialRing(T, 'x')
            sage: R.is_sparse()
            True
        """
        return True

    def interpolation(self, points):
        """
        Return a tropical polynomial with its function is a linear
        interpolation of point in ``points`` if possible.

        If there is only one point, then it will give a constant polynomial.
        Because we are using linear interpolation, each point is actually
        a root of the resulted tropical polynomial.

        INPUT:

        - ``points`` -- a list of tuples ``(x, y)``

        OUTPUT: a :class:`TropicalPolynomial`

        EXAMPLES::

            sage: T = TropicalSemiring(QQ, use_min=True)
            sage: R = PolynomialRing(T, 'x')
            sage: points = [(-2,-3),(1,3),(2,4)]
            sage: p1 = R.interpolation(points); p1
            1*x^2 + 2*x + 4
            sage: p1.plot()
            Graphics object consisting of 1 graphics primitive

        .. PLOT::
            :width: 300 px

            T = TropicalSemiring(QQ, use_min=True)
            R = PolynomialRing(T, 'x')
            points = [(-2,-3),(1,3),(2,4)]
            p1 = R.interpolation(points)
            sphinx_plot(p1.plot())

        ::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R = PolynomialRing(T,'x')
            sage: points = [(0,0),(1,1),(2,4)]
            sage: p1 = R.interpolation(points); p1
            (-2)*x^3 + (-1)*x^2 + 0*x + 0
            sage: p1.plot()
            Graphics object consisting of 1 graphics primitive

        .. PLOT::
            :width: 300 px

            T = TropicalSemiring(QQ, use_min=False)
            R = PolynomialRing(T, 'x')
            points = [(0,0),(1,1),(2,4)]
            p1 = R.interpolation(points)
            sphinx_plot(p1.plot())

        TESTS:

        Every piecewise linear component of tropical polynomial function
        has to have an integer slope::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R = PolynomialRing(T,'x')
            sage: points = [(0,0),(2,3)]
            sage: R.interpolation(points)
            Traceback (most recent call last):
            ...
            ValueError: the slope is not an integer

        For max-plus algebra, the slope of the componenets has to be
        increasing as we move from left to right. Conversely for min-plus
        algebra, the slope of the componenets has to be decreasing from
        left to right::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R = PolynomialRing(T,'x')
            sage: points = [(-2,-3),(1,3),(2,4)]
            sage: R.interpolation(points)
            Traceback (most recent call last):
            ...
            ValueError: can not interpolate these points
        """
        points = sorted(points, key=lambda point: point[0])
        all_slope = [0]
        roots = {}
        R = self.base()
        if R._use_min:
            point_order = range(len(points)-1, 0, -1)
        else:
            point_order = range(len(points)-1)
        for i in point_order:
            if R._use_min:
                slope = (points[i-1][1]-points[i][1]) / (points[i-1][0]-points[i][0])
            else:
                slope = (points[i+1][1]-points[i][1]) / (points[i+1][0]-points[i][0])
            if not slope.is_integer():
                raise ValueError("the slope is not an integer")
            if slope < all_slope[-1]:
                raise ValueError("can not interpolate these points")
            elif slope > all_slope[-1]:
                order = slope - all_slope[-1]
                all_slope.append(slope)
                roots[points[i][0]] = order

        if len(all_slope) == 1:  # constant polynomial
            return self(points[0][1])

        result = self.one()
        for root, order in roots.items():
            result *= self([root, 0])**order
        test_value = result(R(points[0][0]))
        unit = R(points[0][1] - test_value.lift())
        result *= unit
        return result

    @classmethod
    def _implementation_names(cls, implementation, base_ring, sparse):
        """
        Return the only implementation currently available for this semiring,
        which is ``[None]``.

        EXAMPLES::

            sage: from sage.rings.semirings.tropical_polynomial import TropicalPolynomialSemiring
            sage: T = TropicalSemiring(QQ)
            sage: TropicalPolynomialSemiring._implementation_names(None, T, True)
            [None]
        """
        names = [None]
        assert isinstance(names, list)
        assert implementation in names
        return names
