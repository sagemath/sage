r"""
Multivariate Tropical Polynomials

AUTHORS:

- Verrel Rievaldo Wijaya (2024-06): initial version

EXAMPLES::

    sage: T = TropicalSemiring(QQ, use_min=True)
    sage: R.<x,y,z> = PolynomialRing(T)
    sage: R(-1)*x + R(-1)*x + R(-10)*y + R(-3)
    (-1)*x + (-10)*y + (-3)
    sage: (x+y+z)^2
    0*x^2 + 0*x*y + 0*y^2 + 0*x*z + 0*y*z + 0*z^2

REFERENCES:

- [Bru2014]_
- [Fil2017]_
- [Hun2021]_
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
from sage.plot.plot3d.list_plot3d import list_plot3d
from sage.categories.semirings import Semirings

from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.multi_polynomial_element import MPolynomial_polydict
from sage.rings.semirings.tropical_semiring import TropicalSemiring

class TropicalMPolynomial(MPolynomial_polydict):
    r"""
    A multivariate tropical polynomial.

    Let `x_1, x_2, \ldots, x_n` be indeterminants. A tropical monomial is
    any product of these variables, possibly including repetitions:
    `x_1^{i_1}\dots x_n^{i_n}` where `i_j \in \{0,1,\ldots\}`, for all
    `j\in \{1,\ldots,n\}`. A multivariate tropical polynomial is a finite
    linear combination of tropical monomials,
    `p(x_1, \dots, x_n) = \sum_{i=1}^n c_i x_1^{i_1}\dots x_n^{i_n}`.

    In classical arithmetic, we can rewrite the general form of a tropical
    monomial: `x_1^{i_1}\dots x_n^{i_n} = i_1 x_1 + \dots + i_n x_n`. Thus,
    the tropical polynomial can be viewed as the minimum (maximum) of a
    finite collection of linear functions.

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

    Some basic arithmetic operations for tropical polynomials::

        sage: p1 + p2
        3*a*b + 1*a + 1*b
        sage: p1 * p2
        4*a^2*b^2 + 4*a^2*b + 4*a*b^2 + 1*a^2 + 1*a*b + 0*b^2
        sage: T(2) * p1
        5*a*b + 2*a + 1*b
        sage: p1(T(1),T(2))
        6

    Let us look at the different result for tropical curve and 3d graph
    of tropical polynomial in two variables when the min-plus or max-plus
    algebra is used::

        sage: T = TropicalSemiring(QQ, use_min=True)
        sage: R.<a,b> = PolynomialRing(T)
        sage: p1 = R(3)*a*b + a + R(-1)*b
        sage: p1.tropical_variety()
        Tropical curve of 3*a*b + 0*a + (-1)*b

    .. PLOT::
        :width: 300 px

        T = TropicalSemiring(QQ, use_min=True)
        R.<a,b> = PolynomialRing(T)
        p1 = R(3)*a*b + a + R(-1)*b
        tv1 = p1.tropical_variety()
        tv1.plot()
        
    .. PLOT::
        :width: 300 px

        T = TropicalSemiring(QQ, use_min=True)
        R.<a,b> = PolynomialRing(T)
        p1 = R(3)*a*b + a + R(-1)*b
        p1.plot3d()
    
    ::

        sage: T = TropicalSemiring(QQ, use_min=False)
        sage: R.<a,b> = PolynomialRing(T)
        sage: p1 = R(3)*a*b + a + R(-1)*b
        sage: p1.tropical_variety()
        Tropical curve of 3*a*b + 0*a + (-1)*b
    
    .. PLOT::
        :width: 300 px

        T = TropicalSemiring(QQ, use_min=False)
        R.<a,b> = PolynomialRing(T)
        p1 = R(3)*a*b + a + R(-1)*b
        tv1 = p1.tropical_variety()
        tv1.plot()

    .. PLOT::
        :width: 300 px

        T = TropicalSemiring(QQ, use_min=False)
        R.<a,b> = PolynomialRing(T)
        p1 = R(3)*a*b + a + R(-1)*b
        p1.plot3d()

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
        """
        Fix some given variables in ``self`` and return the changed
        tropical multivariate polynomials.

        INPUT:

        -  ``fixed`` -- (optional) dictionary of inputs

        -  ``**kwds`` -- named parameters

        OUTPUT: new :class:`TropicalMPolynomial`

        .. SEEALSO::
        
            :meth:`sage.rings.polynomial.multi_polynomial_element.MPolynomial_polydict.subs

        EXAMPLES::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = x^2 + y + R(3)
            sage: p1((R(4),y))
            0*y + 8
            sage: p1.subs({x: 4})
            0*y + 8
        """
        def check_type(dict, var):
            try:
                variables[i] = T(dict[var])
            except TypeError:
                variables[i] = dict[var]
        
        variables = list(self.parent().gens())
        T = self.parent().base()
        for i in range(len(variables)):
            if str(variables[i]) in kwds:
                check_type(kwds, str(variables[i]))
            elif fixed:
                if variables[i] in fixed:
                    check_type(fixed, variables[i])
                elif i in fixed:
                    check_type(fixed, i)
        if len(kwds) < len(variables):
            for i, var in enumerate(variables):
                if var.parent() is T:
                    variables[i] = self.parent()(var.lift())
                else:
                    variables[i] = self.parent()(var)
        return self(tuple(variables))

    def plot3d(self):
        """
        Return the 3d plot of ``self``.

        Only implemented for tropical polynomial in two variables.
        The `x`-`y` axes for this 3d plot is the same as the `x-y`
        axes of the corresponding tropical curve.

        OUTPUT: Graphics3d Object

        EXAMPLES::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = x^2
        
        .. PLOT::
        :width: 300 px

            T = TropicalSemiring(QQ, use_min=False)
            R.<x,y> = PolynomialRing(T)
            p1 = x^2
            p1.plot3d()

        ::

            sage: p2 = R(3)+R(2)*x+R(2)*y+R(3)*x*y; p2
            3*x*y + 2*x + 2*y + 3

        .. PLOT::
        :width: 300 px

            T = TropicalSemiring(QQ, use_min=False)
            R.<x,y> = PolynomialRing(T)
            p2 = R(3)+R(2)*x+R(2)*y+R(3)*x*y
            p2.plot3d()

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
        from sage.arith.srange import srange

        if len(self.parent().variable_names()) != 2:
            raise NotImplementedError("can only plot the graph of tropical "
                                      "multivariate polynomial in two variables")
        axes = self.tropical_variety()._axes()
        xmin, xmax = axes[0][0], axes[0][1]
        ymin, ymax = axes[1][0], axes[1][1]
        step = 0.5
        x_point = srange(xmin, xmax+step, step)
        y_point = srange(ymin, ymax+step, step)
        res = []
        T = self.parent().base()
        for x in x_point:
            for y in y_point:
                val = self(T(x),T(y)).lift()
                res.append([x,y,val])
        return list_plot3d(res, point_list=True)    

    def tropical_variety(self):
        r"""
        Return tropical roots of ``self``.
        
        In multivariate case, the roots can be represented by a tropical
        variety. For 2 dimensions, it is also called a tropical curve.

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
        Return a nice tropical polynomial string representation.
        
        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: x + R(-1)*y + R(-3)
            0*x + (-1)*y + (-3)
        """
        if not self.dict():
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
        Return the latex representation of this tropical polynomial.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = x^2 + R(-1)*x*y + R(-1)
            sage: latex(p1)
            0 x^{2} + \left(-1\right) x y + \left(-1\right)
            sage: latex(R.zero())
            \infty
        """
        if not self.dict():
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
    Semiring of tropical polynomials in multiple variables.

    Similar to the single-variable case, the set of multivariate tropical
    polynomials `R` also form a semiring because `(R,+)` is a commutative
    monoid and `(R,\cdot)` is a monoid. Additionally, R satisfy the
    distributive and annihilation property. Therefore, this semiring
    extends the concepts to polynomials with multiple variables.
    """

    def __init__(self, base_semiring, n, names, order='degrevlex'):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R = PolynomialRing(T, 5, 'x'); R
            Multivariate Tropical Polynomial Semiring in x0, x1, x2, x3, x4 
            over Rational Field
            sage: TestSuite(R).run()
        """
        from sage.rings.polynomial.term_order import TermOrder

        if not isinstance(base_semiring, TropicalSemiring):
            raise ValueError(f"{base_semiring} is not a tropical semiring")
        Parent.__init__(self, base=base_semiring, names=names, category=Semirings())
        self._ngens = n
        order = TermOrder(order, n)
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

        - ``x`` -- ``dict``, constant, or :class:`MPolynomial` #edit

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
        new_dict = {}
        if isinstance(x, TropicalMPolynomial):
            if x.parent() is not self:
                raise ValueError(f"can not convert {x} to {self}")
        if isinstance(x, MPolynomial):
            if x.parent().variable_names() == self.variable_names():
                x = x.dict()
            else:
                raise ValueError(f"can not convert {x} to {self}")
        elif x in self.base().base_ring(): # constant
            term = [0]*self.ngens()
            x = {tuple(term): x}
        if isinstance(x, dict):
            for key, value in x.items(): # convert coefficient to tropical
                new_dict[key] = self.base()(value)
        return self.element_class(self, new_dict)
    
    def one(self):
        r"""
        Return the multiplicative identity of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R = PolynomialRing(T, 'x')
            sage: R.one()
            0
        """
        return self(0)

    def zero(self):
        r"""
        Return the additive identity of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R = PolynomialRing(T, 'x')
            sage: R.zero()
            +infinity
        """
        return self(self.base().zero())

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

        SEEALSO::

            :meth:`sage.rings.polynomial.multi_polynomial_ring_base.MPolynomialRing_base.random_element`

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<a,b,c> = PolynomialRing(T)
            sage: f = R.random_element()
            sage: f.parent() is R
            True
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        R = PolynomialRing(self.base().base_ring(), self.variable_names())
        return self(R.random_element(degree=degree, terms=terms,
                                     choose_degree=choose_degree, 
                                     *args, **kwargs))
    
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
            gens.append(self({tuple(exponent):self.base()(0)}))
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
        return self._ngens
        