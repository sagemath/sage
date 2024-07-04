r"""
Multivariate Tropical Polynomial Semirings.

This module provides the implementation of parent and element class for 
multivariate tropical polynomials. When working with multivariate case, the
tropical roots is no longer a point. Instead it become a curve in 2d, a
surface in 3d, and a hypersurface in higher dimension.

AUTHORS:

- Verrel Rievaldo Wijaya

EXAMPLES:

Construct a multivariate tropical polynomial semiring in two variables::

    sage: T = TropicalSemiring(QQ, use_min=False)
    sage: R.<a,b> = PolynomialRing(T); R
    Multivariate Tropical Polynomial Semiring in a, b over Rational Field
    
Define some multivariate tropical polynomials::

    sage: p1 = R(3)*a*b + a + R(-1)*b
    3*a*b + 0*a + (-1)*b
    sage: p2 = R(1)*a + R(1)*b + R(1)*a*b; p2
    1*a*b + 1*a + 1*b

Some basic arithmetic operations for tropical polynomials::

    sage: p1 + p2
    3*a*b + 1*a + 1*b
    sage: p1 * p2
    4*a^2*b^2 + 4*a^2*b + 1*a^2 + 4*a*b^2 + 1*a*b + 0*b^2
    sage: p2^2
    2*a^2*b^2 + 2*a^2*b + 2*a*b^2 + 2*a^2 + 2*a*b + 2*b^2
    sage: T(2) * p1
    5*a*b + 2*a + 1*b
    sage: p1(T(1),T(2))
    6

Let's look at the different result for tropical curve and 3d graph of tropical
polynomial in two variables when the min-plus or max-plus algebra is used::

    sage: T = TropicalSemiring(QQ, use_min=True)
    sage: R.<a,b> = PolynomialRing(T)
    sage: p1 = R(3)*a*b + a + R(-1)*b
    sage: p1.tropical_variety()
    Tropical curve of 3*a*b + 0*a + (-1)*b are 
    [[(t1, -3), [t1 <= -4], 1]
    [(-4, t1), [t1 <= -3], 1]
    [(t1 - 1, t1), [t1 >= -3], 1]]
    sage: plot(p1.tropical_variety())
    sage: p1.plot3d()

    sage: T = TropicalSemiring(QQ, use_min=False)
    sage: R.<a,b> = PolynomialRing(T)
    sage: p1 = R(3)*a*b + a + R(-1)*b
    sage: p1.tropical_variety()
    Tropical curve of 3*a*b + 0*a + (-1)*b are 
    [[(t1, -3), [t1 >= -4], 1]
    [(-4, t1), [t1 >= -3], 1]
    [(t1 - 1, t1), [t1 <= -3], 1]]
    sage: plot(p1.tropical_variety())
    sage: p1.plot3d()

TESTS:

There is no subtraction defined for tropical polynomials::

    sage: -p1
    Traceback (most recent call last):
    ...
    ArithmeticError: cannot negate any non-infinite element

REFERENCES:

    - [Bru2013]_
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
from sage.plot.plot3d.list_plot3d import list_plot3d
from sage.symbolic.ring import SR
from sage.categories.sets_cat import Sets

from sage.rings.polynomial.term_order import TermOrder
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.multi_polynomial_element import MPolynomial_polydict
from sage.rings.semirings.tropical_semiring import TropicalSemiring
from sage.rings.semirings.tropical_variety import TropicalCurve, TropicalVariety

class TropicalMPolynomial(MPolynomial_polydict):
    r"""
    Generic multivariate tropical polynomial.

    """

    def plot3d(self):
        """
        Return the 3d plot of ``self``.

        OUTPUT: A Graphics3d Object.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = R(3)+R(2)*x+R(2)*y+R(3)*x*y; p1
            3*x*y + 2*x + 2*y + 3
            sage: p1.plot3d()

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
            raise NotImplementedError("can only plot the graph of tropical" \
                                " multivariate polynomial in two variables")
        axes = self.tropical_variety()._axes()
        xmin, xmax = axes[0][0], axes[0][1]
        ymin, ymax = axes[1][0], axes[1][1]
        step = 0.5
        x_point = srange(xmin-1, xmax+1+step, step)
        y_point = srange(ymin-1, ymax+1+step, step)
        res = []
        T = self.parent().base()
        for x in x_point:
            for y in y_point:
                val = self(T(x),T(y)).lift()
                res.append([x,y,val])
        return list_plot3d(res, point_list=True)    

    def tropical_variety(self):
        r"""
        Return tropical roots of ``self``. In multivariate case, the roots
        can be represented by a tropical variety. For 2 dimensions, it is 
        also called a tropical curve.

        OUTPUT: ``TropicalVariety`` object.
        
        EXAMPLES:

        Some examples of tropical curve for tropical polynomials in two 
        variables::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = x + y + R(0); p1
            0*x + 0*y + 0
            sage: p1.tropical_variety()
            Tropical curve of 0*x + 0*y + 0 are 
            [[(t1, t1), [t1 >= 0], 1]
            [(0, t1), [t1 <= 0], 1]
            [(t1, 0), [t1 <= 0], 1]]

        ::

            sage: T = TropicalSemiring(QQ, use_min=True)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p2 = R(-1)*x^2 + x + y^2 + R(0); p2
            (-1)*x^2 + 0*y^2 + 0*x + 0
            sage: p2.tropical_variety()
            Tropical curve of (-1)*x^2 + 0*x + 0*y^2 + 0 are 
            Tropical curve of (-1)*x^2 + 0*y^2 + 0*x + 0 are 
            [[(t1 + 1/2, t1), [t1 <= 0], 2]
            [(1/2, t1), [t1 >= 0], 2]
            [(t1, 0), [(1/2) <= t1], 2]]

        We can also find tropical hypersurface for any tropical polynomials 
        in `n\geq 2` variables::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y,z> = PolynomialRing(T)
            sage: p1 = R(1)*x*y + R(-1/2)*x*z + R(4)*z^2; p1
            1*x*y + (-1/2)*x*z + 4*z^2
            sage: p1.tropical_variety()
            Tropical hypersurface of 1*x*y + (-1/2)*x*z + 4*z^2 are 
            [[(t1, t2 - 3/2, t2), [t1 <= t2 + 9/2], 1]
            [(2*t1 - t2 + 3, t2, t1), [t2 + 3/2 <= t1], 1]
            [(t1 + 9/2, t2, t1), [t1 <= t2 + 3/2], 1]]
        """
        if self.parent().ngens() == 2:
            return TropicalCurve(self)
        else:
            return TropicalVariety(self)
    
    def _latex_(self):
        r"""
        Return a nice topical polynomial latex representation.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y,z> = PolynomialRing(T)
            sage: f = x^2 +R(-3)* y + R(1)*z; f
            0*x^2 + (-3)*y + 1*z
            sage: latex(f)
            0 x^{2} \oplus \left(-3\right) y \oplus 1 z
        """
        s = super()._latex_()
        s = s.replace("+", r'\oplus')
        return s

class TropicalMPolynomialSemiring(UniqueRepresentation, Parent):
    """
    Semiring structure of tropical polynomials in multiple variables.
    """

    def __init__(self, base_semiring, n, names, order='degrevlex'):
        """
        EXAMPLES::

            sage: T = TropicalSemiring(QQ, use_min=True)
            sage: R.<x,y,z> = PolynomialRing(T)
            sage: R.gens()
            (0*x, 0*y, 0*z)
            sage: R(1)*x*y*z + x
            1*x*y*z + 0*x
            sage: (x+y+z)^2
            0*x^2 + 0*x*y + 0*y^2 + 0*x*z + 0*y*z + 0*z^2
        """
        if not isinstance(base_semiring, TropicalSemiring):
            raise ValueError(f"{base_semiring} is not a tropical semiring")
        Parent.__init__(self, base=base_semiring, names=names, category=Sets())
        self._ngens = n
        order = TermOrder(order, n)
        self._term_order = order
        
    Element = TropicalMPolynomial

    def _element_constructor_(self, x):
        """"
        Convert ``x`` into this tropical multivariate polynomial semiring

        INPUT:

        - ``x`` -- dict or MPolynomial

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
        """
        C = self.element_class
        new_dict = {}
        if isinstance(x, MPolynomial):
            x = x.dict()
        elif x in self.base().base_ring(): # constant
            term = [0 for _ in range(self.ngens())]
            x = {tuple(term):x}
        
        if isinstance(x, dict):
            for key, value in x.items(): # convert each coefficient to tropical
                new_dict[key] = self.base()(value)

        return C(self, new_dict)
    
    def _repr_(self):
        if self._ngens == 0:
            return (f"Multivariate Tropical Polynomial Semiring in no variables"
            f" over {self.base_ring().base_ring()}")
        return (f"Multivariate Tropical Polynomial Semiring in {', '.join(self.variable_names())}"
            f" over {self.base_ring().base_ring()}")
    
    def term_order(self):
        return self._term_order
    
    def random_element(self, degree=2, terms=None, choose_degree=False,
                       *args, **kwargs):
        """
        Return a random multivariate tropical polynomial.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<a,b> = PolynomialRing(T)
            sage: f = R.random_element(); f
            1/9*a^2 + 1/13*a*b + 1/107*b^2 + 1*a
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        R = PolynomialRing(self.base().base_ring(), self.variable_names())
        return self(R.random_element(degree=degree, terms=terms, choose_degree=choose_degree, 
                                     *args, **kwargs))
    
    def gen(self, n=0):
        """
        Return the indeterminate generator of this polynomial semiring.
        """
        return self.gens()[n]
    
    @cached_method
    def gens(self):
        """
        Return a tuple whose entries are the generators for this object,
        in order.
        """
        gens = []
        for i in range(self.ngens()):
            exponent = [0] * self.ngens()
            exponent[i] = 1
            gens.append(self({tuple(exponent):self.base()(0)}))
        return tuple(gens)
    
    def ngens(self):
        """
        Return the number of generators of this polynomial semiring.
        """
        return self._ngens
        