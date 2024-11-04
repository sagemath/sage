# sage.doctest: needs sage.symbolic
r"""
Multivariate Tropical Polynomials

AUTHORS:

- Verrel Rievaldo Wijaya (2024-06): initial version

EXAMPLES::

    sage: T = TropicalSemiring(QQ, use_min=True)
    sage: R.<x,y,z> = PolynomialRing(T)
    sage: z.parent()
    Multivariate Tropical Polynomial Semiring in x, y, z over Rational Field
    sage: R(2)*x + R(-1)*x + R(5)*y + R(-3)
    (-1)*x + 5*y + (-3)
    sage: (x+y+z)^2
    0*x^2 + 0*x*y + 0*y^2 + 0*x*z + 0*y*z + 0*z^2

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
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.polynomial.multi_polynomial_element import MPolynomial_polydict


class TropicalMPolynomial(MPolynomial_polydict):
    r"""
    A multivariate tropical polynomial.

    Let `x_1, x_2, \ldots, x_n` be indeterminants. A tropical monomial is
    any product of these variables, possibly including repetitions:
    `x_1^{i_1}\cdots x_n^{i_n}` where `i_j \in \{0,1,\ldots\}`, for all
    `j\in \{1,\ldots,n\}`. A multivariate tropical polynomial is a finite
    linear combination of tropical monomials,
    `p(x_1, \ldots, x_n) = \sum_{i=1}^n c_i x_1^{i_1}\cdots x_n^{i_n}`.

    In classical arithmetic, we can rewrite the general form of a tropical
    monomial: `x_1^{i_1}\cdots x_n^{i_n} \mapsto i_1 x_1 + \cdots + i_n x_n`.
    Thus, the tropical polynomial can be viewed as the minimum (maximum) of
    a finite collection of linear functions.

    EXAMPLES:

    Construct a multivariate tropical polynomial semiring in two variables::

        sage: T = TropicalSemiring(QQ, use_min=False)
        sage: R.<a,b> = PolynomialRing(T); R
        Multivariate Tropical Polynomial Semiring in a, b over Rational Field

    Define some multivariate tropical polynomials::

        sage: p1 = R(3)*a*b + a + R(-1)*b; p1
        3*a*b + 0*a + (-1)*b
        sage: p2 = R(1)*a + R(1)*b + R(1)*a*b; p2
        1*a*b + 1*a + 1*b

    Some basic arithmetic operations for multivariate tropical polynomials::

        sage: p1 + p2
        3*a*b + 1*a + 1*b
        sage: p1 * p2
        4*a^2*b^2 + 4*a^2*b + 4*a*b^2 + 1*a^2 + 1*a*b + 0*b^2
        sage: T(2) * p1
        5*a*b + 2*a + 1*b
        sage: p1(T(1),T(2))
        6

    Let us look at the different result for tropical curve and 3d graph
    of tropical polynomial in two variables when different algebra is used.
    First for the min-plus algebra::

        sage: T = TropicalSemiring(QQ, use_min=True)
        sage: R.<a,b> = PolynomialRing(T)
        sage: p1 = R(3)*a*b + a + R(-1)*b
        sage: p1.tropical_variety()
        Tropical curve of 3*a*b + 0*a + (-1)*b
        sage: p1.tropical_variety().plot()
        Graphics object consisting of 3 graphics primitives

    .. PLOT::
        :width: 300 px

        T = TropicalSemiring(QQ, use_min=True)
        R = PolynomialRing(T, ('a,b'))
        a, b = R.gen(), R.gen(1)
        p1 = R(3)*a*b + a + R(-1)*b
        tv1 = p1.tropical_variety()
        sphinx_plot(tv1.plot())

    Tropical polynomial in two variables will induce a function in three
    dimension that consists of a number of surfaces::

        sage: p1.plot3d()
        Graphics3d Object

    .. PLOT::
        :width: 300 px

        T = TropicalSemiring(QQ, use_min=True)
        R = PolynomialRing(T, ('a,b'))
        a, b = R.gen(), R.gen(1)
        p1 = R(3)*a*b + a + R(-1)*b
        sphinx_plot(p1.plot3d())

    If we use a max-plus algebra, we will get a slightly different result::

        sage: T = TropicalSemiring(QQ, use_min=False)
        sage: R.<a,b> = PolynomialRing(T)
        sage: p1 = R(3)*a*b + a + R(-1)*b
        sage: p1.tropical_variety()
        Tropical curve of 3*a*b + 0*a + (-1)*b
        sage: p1.tropical_variety().plot()
        Graphics object consisting of 3 graphics primitives

    .. PLOT::
        :width: 300 px

        T = TropicalSemiring(QQ, use_min=False)
        R = PolynomialRing(T, ('a,b'))
        a, b = R.gen(), R.gen(1)
        p1 = R(3)*a*b + a + R(-1)*b
        tv1 = p1.tropical_variety()
        sphinx_plot(tv1.plot())

    ::

        sage: p1.plot3d()
        Graphics3d Object

    .. PLOT::
        :width: 300 px

        T = TropicalSemiring(QQ, use_min=False)
        R = PolynomialRing(T, ('a,b'))
        a, b = R.gen(), R.gen(1)
        p1 = R(3)*a*b + a + R(-1)*b
        sphinx_plot(p1.plot3d())

    TESTS:

    There is no subtraction defined for tropical polynomials::

        sage: T = TropicalSemiring(QQ)
        sage: R.<a,b> = PolynomialRing(T)
        sage: a - b
        Traceback (most recent call last):
        ...
        ArithmeticError: cannot negate any non-infinite element
    """
    def subs(self, fixed=None, **kwds):
        r"""
        Fix some given variables in ``self`` and return the changed
        tropical multivariate polynomials.

        .. SEEALSO::

            :meth:`sage.rings.polynomial.multi_polynomial_element.MPolynomial_polydict.subs`

        EXAMPLES::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = x^2 + y + R(3)
            sage: p1((R(4),y))
            0*y + 8
            sage: p1.subs({x: 4})
            0*y + 8
        """
        variables = list(self.parent().gens())
        for i in range(len(variables)):
            if str(variables[i]) in kwds:
                variables[i] = kwds[str(variables[i])]
            elif fixed:
                if variables[i] in fixed:
                    variables[i] = fixed[variables[i]]
                elif i in fixed:
                    variables[i] = fixed[i]
        if len(kwds) < len(variables):
            for i, v in enumerate(variables):
                variables[i] = self.parent()(v)
        return self(tuple(variables))

    def plot3d(self, color='random'):
        """
        Return the 3d plot of ``self``.

        Only implemented for tropical polynomial in two variables.
        The `x`-`y` axes for this 3d plot is the same as the `x`-`y`
        axes of the corresponding tropical curve.

        OUTPUT: Graphics3d Object

        EXAMPLES:

        A simple tropical polynomial that consist of only one surface::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = x^2
            sage: p1.plot3d()
            Graphics3d Object

        .. PLOT::
            :width: 300 px

            T = TropicalSemiring(QQ, use_min=False)
            R = PolynomialRing(T, ('x,y'))
            x, y = R.gen(), R.gen(1)
            p1 = x**2
            sphinx_plot(p1.plot3d())

        Tropical polynomials often have graphs that represent a combination
        of multiple surfaces::

            sage: p2 = R(3) + R(2)*x + R(2)*y + R(3)*x*y
            sage: p2.plot3d()
            Graphics3d Object

        .. PLOT::
            :width: 300 px

            T = TropicalSemiring(QQ, use_min=False)
            R = PolynomialRing(T, ('x,y'))
            x, y = R.gen(), R.gen(1)
            p2 = R(3) + R(2)*x + R(2)*y + R(3)*x*y
            sphinx_plot(p2.plot3d())

        ::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p3 = R(2)*x^2 + x*y + R(2)*y^2 + x + R(-1)*y + R(3)
            sage: p3.plot3d()
            Graphics3d Object

        .. PLOT::
            :width: 300 px

            T = TropicalSemiring(QQ)
            R = PolynomialRing(T, ('x,y'))
            x, y = R.gen(), R.gen(1)
            p3 = R(2)*x**2 + x*y + R(2)*y**2 + x + R(-1)*y + R(3)
            sphinx_plot(p3.plot3d())

        TESTS::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y,z> = PolynomialRing(T)
            sage: p1 = x*y*z + x
            sage: p1.plot3d()
            Traceback (most recent call last):
            ...
            NotImplementedError: can only plot the graph of tropical
            multivariate polynomial in two variables
        """
        from random import random
        from sage.plot.graphics import Graphics
        from sage.geometry.polyhedron.constructor import Polyhedron
        from sage.sets.real_set import RealSet
        from sage.symbolic.relation import solve

        if len(self.parent().variable_names()) != 2:
            raise NotImplementedError("can only plot the graph of tropical "
                                      "multivariate polynomial in two variables")
        tv = self.tropical_variety()
        axes = tv._axes()
        edge = set()
        if tv.components():
            v = tv._vars[0]
        T = self.parent().base()
        R = self.base_ring().base_ring()

        # Finding the point of curve that touch the edge of the axes
        for comp in tv.components():
            if len(comp[1]) == 1:
                valid_int = RealSet(comp[1][0])
            else:
                valid_int = RealSet(comp[1][0]).intersection(RealSet(comp[1][1]))
            for i, eqn in enumerate(comp[0]):
                j = (i+1) % 2
                if not eqn.is_numeric():
                    for k in range(2):
                        sol = solve(eqn == axes[i][k], v)
                        if sol[0].rhs() in valid_int:
                            valid_point = [R(eq.subs(**{str(v): sol[0].rhs()})) for eq in comp[0]]
                            if valid_point[j] in RealSet(axes[j]):
                                edge.add(tuple(valid_point))

        # Combine the edge, vertices, and corner point
        vertices = self.tropical_variety().vertices()
        corner = set()
        for i in axes[0]:
            for j in axes[1]:
                corner.add((i, j))
        marks = corner | vertices | edge

        # Calculate the value of polynomial at each marked point
        variables = self.parent().gens()
        terms = [a*variables[0]**b[0] * variables[1]**b[1] for a, b in zip(self.coefficients(), self.exponents())]
        point_terms = {}
        for mark in marks:
            mark_terms = []
            value = self(T(mark[0]), T(mark[1]))
            value_terms = [term(T(mark[0]), T(mark[1])) for term in terms]
            mark_terms.extend(terms[i] for i in range(len(terms))
                              if value_terms[i] == value)
            point_terms[(R(mark[0]), R(mark[1]), value.lift())] = mark_terms

        # Plot the points that attained its value at one term only
        combined_plot = Graphics()
        for elms in point_terms.values():
            if len(elms) == 1:
                poly_vert = []
                term = elms[0]
                for p, t in point_terms.items():
                    if term in t:
                        poly_vert.append(p)
                        t.remove(term)
                if color == 'random':
                    rand_color = (random(), random(), random())
                plot = Polyhedron(vertices=poly_vert).plot(color=rand_color)
                combined_plot += plot

        # Plot the remaining points
        for remain in point_terms.values():
            for term in remain:
                poly_vert = []
                for p, t in point_terms.items():
                    if term in t:
                        poly_vert.append(p)
                        t.remove(term)
                if color == 'random':
                    rand_color = (random(), random(), random())
                plot = Polyhedron(vertices=poly_vert).plot(color=rand_color)
                combined_plot += plot
        return combined_plot

    def tropical_variety(self):
        r"""
        Return tropical roots of ``self``.

        In the multivariate case, the roots can be represented by a
        tropical variety. In two dimensions, this is known as a tropical
        curve. For dimensions higher than two, it is referred to as a
        tropical hypersurface.

        OUTPUT: a :class:`sage.rings.semirings.tropical_variety.TropicalVariety`

        EXAMPLES:

        Tropical curve for tropical polynomials in two variables::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = x + y + R(0); p1
            0*x + 0*y + 0
            sage: p1.tropical_variety()
            Tropical curve of 0*x + 0*y + 0

        Tropical hypersurface for tropical polynomials in more than two
        variables::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y,z> = PolynomialRing(T)
            sage: p1 = R(1)*x*y + R(-1/2)*x*z + R(4)*z^2; p1
            1*x*y + (-1/2)*x*z + 4*z^2
            sage: p1.tropical_variety()
            Tropical surface of 1*x*y + (-1/2)*x*z + 4*z^2
        """
        from sage.rings.semirings.tropical_variety import TropicalCurve, TropicalSurface, TropicalVariety

        if self.parent().ngens() == 2:
            return TropicalCurve(self)
        if self.parent().ngens() == 3:
            return TropicalSurface(self)
        return TropicalVariety(self)

    def _repr_(self):
        r"""
        Return string representation of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: x + R(-1)*y + R(-3)
            0*x + (-1)*y + (-3)
        """
        if not self.monomial_coefficients():
            return str(self.parent().base().zero())
        s = super()._repr_()
        if self.monomials()[-1].is_constant():
            if self.monomial_coefficient(self.parent()(0)) < 0:
                s = s.replace(" - ", " + -")
                const = str(self.monomial_coefficient(self.parent(0)))
                s = s.replace(f" {const}", f" ({const})")
        return s

    def _latex_(self):
        r"""
        Return the latex representation of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = x^2 + R(-1)*x*y + R(-1)
            sage: latex(p1)
            0 x^{2} + \left(-1\right) x y + \left(-1\right)
            sage: latex(R.zero())
            \infty
        """
        if not self.monomial_coefficients():
            return self.parent().base().zero()._latex_()
        s = super()._latex_()
        if self.monomials()[-1].is_constant():
            if self.monomial_coefficient(self.parent()(0)) < 0:
                s = s.replace(" - ", " + -")
                const = str(self.monomial_coefficient(self.parent(0)))
                s = s.replace(f" {const}", f" \\left({const}\\right)")
        return s


class TropicalMPolynomialSemiring(UniqueRepresentation, Parent):
    r"""
    The semiring of tropical polynomials in multiple variables.

    This is the commutative semiring consisting of all finite linear
    combinations of tropical monomials under (tropical) addition
    and multiplication with coefficients in a tropical semiring.

    EXAMPLES::

        sage: T = TropicalSemiring(QQ)
        sage: R.<x,y> = PolynomialRing(T)
        sage: f = T(1)*x + T(-1)*y
        sage: g = T(2)*x + T(-2)*y
        sage: f + g
        1*x + (-2)*y
        sage: f * g
        3*x^2 + (-1)*x*y + (-3)*y^2
        sage: f + R.zero() == f
        True
        sage: f * R.zero() == R.zero()
        True
        sage: f * R.one() == f
        True
    """
    def __init__(self, base_semiring, n, names, order):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R = PolynomialRing(T, 5, 'x')
            sage: TestSuite(R).run()
        """
        from sage.rings.semirings.tropical_semiring import TropicalSemiring
        from sage.categories.semirings import Semirings
        if not isinstance(base_semiring, TropicalSemiring):
            raise ValueError(f"{base_semiring} is not a tropical semiring")
        Parent.__init__(self, base=base_semiring, names=names, category=Semirings())
        self._ngens = n
        self._term_order = order

    def term_order(self):
        """
        Return the defined term order of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y,z> = PolynomialRing(T)
            sage: R.term_order()
            Degree reverse lexicographic term order
        """
        return self._term_order

    Element = TropicalMPolynomial

    def _element_constructor_(self, x):
        r""""
        Convert ``x`` into ``self``.

        INPUT:

        - ``x`` -- ``dict``, constant, or :class:`MPolynomial`

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R = PolynomialRing(T, 'x,y')
            sage: dict1 = {(1,0):0, (0,1):-1, (1,1):3}
            sage: p1 = R(dict1); p1
            3*x*y + 0*x + (-1)*y
            sage: S.<x,y> = PolynomialRing(QQ)
            sage: f = -x*y + 1
            sage: R(f)
            (-1)*x*y + 1

        TESTS::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: S.<a,b> = PolynomialRing(T)
            sage: R(a + b)
            Traceback (most recent call last):
            ...
            ValueError: can not convert 0*a + 0*b to Multivariate Tropical
            Polynomial Semiring in x, y over Rational Field
        """
        from sage.rings.polynomial.multi_polynomial import MPolynomial
        if isinstance(x, TropicalMPolynomial):
            if x.parent() is not self:
                raise ValueError(f"can not convert {x} to {self}")
        if isinstance(x, MPolynomial):
            if x.parent().variable_names() == self.variable_names():
                x = x.monomial_coefficients()
            else:
                raise ValueError(f"can not convert {x} to {self}")
        elif (x in self.base().base_ring()) or (x in self.base()):
            term = [0] * self.ngens()
            x = {tuple(term): x}

        if isinstance(x, dict):
            for key, value in x.items():
                x[key] = self.base()(value)
        return self.element_class(self, x)

    @cached_method
    def one(self):
        r"""
        Return the multiplicative identity of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R = PolynomialRing(T, 'x,y')
            sage: R.one()
            0
        """
        exponent = [0] * self.ngens()
        return self.element_class(self, {tuple(exponent): self.base().one()})

    @cached_method
    def zero(self):
        r"""
        Return the additive identity of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R = PolynomialRing(T, 'x,y')
            sage: R.zero()
            +infinity
        """
        exponent = [0] * self.ngens()
        return self.element_class(self, {tuple(exponent): self.base().zero()})

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(RR)
            sage: R.<u,v,w> = PolynomialRing(T); R
            Multivariate Tropical Polynomial Semiring in u, v, w over Real
            Field with 53 bits of precision
        """
        if self._ngens == 0:
            return (f"Multivariate Tropical Polynomial Semiring in no variables"
                    f" over {self.base_ring().base_ring()}")
        return (f"Multivariate Tropical Polynomial Semiring in {', '.join(self.variable_names())}"
                f" over {self.base_ring().base_ring()}")

    def random_element(self, degree=2, terms=None, choose_degree=False,
                       *args, **kwargs):
        r"""
        Return a random multivariate tropical polynomial from ``self``.

        OUTPUT: a :class:`TropicalMPolynomial`

        .. SEEALSO::

            :meth:`sage.rings.polynomial.multi_polynomial_ring_base.MPolynomialRing_base.random_element`

        EXAMPLES:

        A random polynomial of at most degree `d` and at most `t` terms::

            sage: T = TropicalSemiring(QQ)
            sage: R.<a,b,c> = PolynomialRing(T)
            sage: f = R.random_element(2, 5)
            sage: f.degree() <= 2
            True
            sage: f.parent() is R
            True
            sage: len(list(f)) <= 5
            True

        Choose degrees of monomials randomly first rather than monomials
        uniformly random::

            sage: f = R.random_element(3, 6, choose_degree=True)
            sage: f.degree() <= 3
            True
            sage: f.parent() is R
            True
            sage: len(list(f)) <= 6
            True
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        R = PolynomialRing(self.base().base_ring(), self.variable_names())
        f = R.random_element(degree=degree, terms=terms, choose_degree=choose_degree,
                             *args, **kwargs)
        new_dict = {key: self.base()(value)
                    for key, value in f.monomial_coefficients().items()}
        return self.element_class(self, new_dict)

    def gen(self, n=0):
        r"""
        Return the ``n``-th generator of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<a,b,c> = PolynomialRing(T)
            sage: R.gen()
            0*a
            sage: R.gen(2)
            0*c

        TESTS::

            sage: R.gen(3)
            Traceback (most recent call last):
            ...
            IndexError: tuple index out of range
        """
        return self.gens()[n]

    @cached_method
    def gens(self):
        r"""
        Return the generators of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R = PolynomialRing(T, 5, 'x')
            sage: R.gens()
            (0*x0, 0*x1, 0*x2, 0*x3, 0*x4)
        """
        gens = []
        for i in range(self.ngens()):
            exponent = [0] * self.ngens()
            exponent[i] = 1
            element = self.element_class(self, {tuple(exponent): self.base().one()})
            gens.append(element)
        return tuple(gens)

    def ngens(self):
        r"""
        Return the number of generators of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R = PolynomialRing(T, 10, 'z')
            sage: R.ngens()
            10
        """
        from sage.rings.integer_ring import ZZ
        return ZZ(self._ngens)
