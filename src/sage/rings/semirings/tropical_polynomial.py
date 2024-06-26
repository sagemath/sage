r"""
Univariate Tropical Polynomial Semirings

Tropical polynomial is a polynomial with coefficients from tropical semiring.
Tropical polynomial induces a function which is piecewise-linear and each 
piece has an integer slope. Tropical roots (zeros) of polynomial `P(x)` is 
defined as all points ``x_0`` for which the graph of ``P(x)`` change its slope.
The difference in the slopes of the two pieces adjacent to this root gives 
the order of the corresponding root. This module provides the implementation
of parent and element class for sparse tropical polynomials in one variable.

AUTHORS:

- Verrel Rievaldo Wijaya

EXAMPLES:

Construct a tropical polynomial semiring by first defining a base 
semiring and then inputting it to ``PolynomialRing`` constructor::

    sage: T = TropicalSemiring(QQ, use_min=False)
    sage: R = PolynomialRing(T,'y')
    sage: R
    Tropical Polynomial Semiring in y over Rational Field

We can define the element by giving a list or tuple of coefficients
that begins with constant. This is also the way to create a tropical
polynomial with `0` as coefficient::

    sage: p1 = R([1,4,None,0]); p1
    0*y^3 + 4*y + 1

Create an element by converting from classical polynomial::

    sage: S.<y> = PolynomialRing(QQ)
    sage: p2 = R(y^2+2*y+3); p2
    1*y^2 + 2*y + 3

We can do the addition, multiplication, and evaluation for tropical 
polynomials. When doing evaluation, make sure the input number is tropical.
If not, then it will raise an error::

    sage: p1 + p2
    0*y^3 + 1*y^2 + 4*y + 3
    sage: p1 * p2
    1*y^5 + 2*y^4 + 5*y^3 + 6*y^2 + 7*y + 4
    sage: p1(3)
    Traceback (most recent call last):
    ...
    TypeError: no common canonical parent for objects with parents: 
    'Tropical semiring over Rational Field' and 'Integer Ring'
    sage: p1(T(3))
    9

Beware that when multiplying tropical polynomial with a scalar, it
will raise an error if the scalar is not tropical number::

    sage: 2 * p1
    Traceback (most recent call last):
    ...
    ArithmeticError: cannot negate any non-infinite element
    sage: T(2) * p1
    2*y^3 + 6*y + 3 

We can also find all the tropical roots of tropical polynomial counted
with multiplicity. There will be no tropical root for constant polynomial. 
For a monomial, the tropical root is the additive identity of its base 
tropical semiring::

    sage: p1.roots()
    [-3, 2, 2]
    sage: p2.roots()
    [1, 1]
    sage: p3 = R(1)
    sage: p3.roots()
    []
    sage: p4 = R(y^2)
    sage: p4.roots()
    [-infinity, -infinity]

Even though some tropical polynomials have tropical roots, this does not
neccessarily means it can be factored into its linear factors::

    sage: p1.factor()
    (0) * (0*y^3 + 4*y + 1)
    sage: p2.factor()
    (1) * (0*y + 1)^2

Every tropical polynomial `p(x)` have a corresponding unique tropical 
polynomial `\bar{p}(x)` with the same roots which can be factored. Therefore
this two polynomial determine the same function. We call `\bar{p}(x)` the
conjugate of `p(x)`. The conjugate can be a convex or cancave function::

    sage: p1.convex_conjugate()
    0*y^3 + 2*y^2 + 4*y + 1
    sage: p1.concave_conjugate()
    0*y^3 + -3*y^2 + -1*y + 1
    sage: p1.convex_conjugate().factor()
    (0) * (0*y + -3) * (0*y + 2)^2
    
Because we are using max-plus algebra, then we can check that the induced 
tropical polynomial function of `p(x)` and its convex conjugate is equal::

    sage: p1.piecewise_function()
    piecewise(x|-->1 on (-oo, -3], x|-->x + 4 on (-3, 2), x|-->3*x on 
    [2, +oo); x)
    sage: p1.convex_conjugate().piecewise_function()
    piecewise(x|-->1 on (-oo, -3], x|-->x + 4 on (-3, 2), x|-->3*x on 
    (2, +oo); x)

Plot the graph of some tropical polynomials::
    sage: p1.plot()
    sage: plot(p2, xmin=0, xmax=2)

TESTS:

There is no subtraction for tropical polynomials because element in tropical 
semiring doesn't necessarily have additive inverse::

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

from sage.sets.real_set import RealSet
from sage.symbolic.ring import SR
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.rings.polynomial.polynomial_element_generic import \
Polynomial_generic_sparse
from sage.categories.semirings import Semirings
from sage.rings.semirings.tropical_semiring import TropicalSemiring
from itertools import combinations

class TropicalPolynomial(Polynomial_generic_sparse):
    """
    A generic sparse tropical polynomial.

    The `TropicalPolynomial`` class defines functionality for sparse
    polynomials over any tropical semiring. A sparse polynomial is 
    represented using a dictionary which maps each exponent to the
    corresponding coefficient. The coefficients is a tropical number.
    """                                                                            
    
    def roots(self):
        r"""
        Return the list of all tropical roots of ``self``

        OUTPUT:

        - ``tropical_roots`` -- list; Contains tropical roots of ``self``
         after being sorted counted with multiplicity

        EXAMPLES:

            sage: T = TropicalSemiring(QQ, use_min=True)
            sage: R = PolynomialRing(T, x)
            sage: p1 = R([5,4,1,0,2,4,3]); p1
            3*x^6 + 4*x^5 + 2*x^4 + 0*x^3 + 1*x^2 + 4*x + 5
            sage: p1.roots()
            [-1, -1, -1, 1, 2, 2]
            sage: p2 = R([0, None, 0]); p2
            0*x^2 + 0
            sage: p2.roots()
            [0, 0]

        """
        tropical_roots = []
        if len(self.dict()) == 1:
            exponent = list(self.dict())[0]
            if exponent == 0:
                return tropical_roots
            else:
                return [self.parent().base_ring().zero()]*exponent
        
        R = self.parent().base().base_ring()
        dict_root = {}
        dict_coeff = {i:c.lift() for i,c in self.dict().items()}
        for comb in combinations(dict_coeff, 2):
            index1, index2 = comb[0], comb[1]
            root = (dict_coeff[index1]-dict_coeff[index2])/(index2 - index1)
            val_root = dict_coeff[index1] + index1*root
            check_maks = True
            for key in dict_coeff:
                if key not in comb:
                    val = dict_coeff[key] + key*root
                    if self.base_ring()._use_min:
                        if val < val_root:
                            check_maks = False
                            break
                    else:
                        if val > val_root:
                            check_maks = False
                            break
            if check_maks:
                order = abs(index1-index2)
                if root not in  dict_root:
                    dict_root[root] = order
                else:
                    if order > dict_root[root]:
                        dict_root[root] = order
        
        for root in dict_root:
            tropical_roots += [root] * dict_root[root]
            
        return sorted(tropical_roots)
    
    def _conjugate(self, use_min):
        r"""
        Return the tropical polynomial which has the same roots as ``self``

        INPUT:

        - ``use_min`` -- bool. If `True` then the base tropical semiring
        will use a min-plus algebra. Otherwise it use a max-plus algebra
        
        OUTPUT: TropicalPolynomial

        EXAMPLES:

            sage: T = TropicalSemiring(QQ, use_min=True)
            sage: R = PolynomialRing(T, x)
            sage: p1 = R([5,4,1,0,2,4,3]); p1
            sage: p1._conjugate(True)
            3*x^6 + 2*x^5 + 1*x^4 + 0*x^3 + 1*x^2 + 3*x + 5
            sage: p1._conjugate(False)
            3*x^6 + 5*x^5 + 7*x^4 + 8*x^3 + 7*x^2 + 6*x + 5   

        """
        from sage.rings.polynomial.polynomial_ring_constructor import \
            PolynomialRing
        
        roots = self.roots()
        T = TropicalSemiring(self.parent().base().base_ring(), use_min=use_min)
        R = PolynomialRing(T, self.parent().variable_name())
        poly = R(self.dict()[self.degree()].lift())
        for root in roots:
            linear = R([root, 0])
            poly *= linear
        return poly

    def convex_conjugate(self):
        r"""
        Return the tropical polynomial which has the same roots as ``self``, 
        which can be factored and the function is convex

        OUTPUT: TropicalPolynomial

        EXAMPLES:

            sage: T = TropicalSemiring(QQ, use_min=True)
            sage: R = PolynomialRing(T, x)
            sage: p1 = R([5,3,0])
            sage: p1.factor()
            (0) * (0*x^2 + 3*x + 5)
            sage: p1.is_convex()
            False
            sage: p1.convex_conjugate()
            0*x^2 + 5/2*x + 5
            sage: p1.convex_conjugate().factor()
            (0) * (0*x + 5/2)^2
            sage: p1.convex_conjugate().is_convex()
            True

        """
        return self._conjugate(False)

    def concave_conjugate(self):
        r"""
        Return the tropical polynomial which has the same roots as ``self``, 
        which can be factored and the function is concave

        OUTPUT: TropicalPolynomial

        EXAMPLES:

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R = PolynomialRing(T, x)
            sage: p1 = R([4,0,2])
            sage: p1.factor()
            (2) * (0*x^2 + -2*x + 2)
            sage: p1.is_concave()
            False
            sage: p1.concave_conjugate()
            2*x^2 + 3*x + 4
            sage: p1.concave_conjugate().factor()
            (2) * (0*x + 1)^2
            sage: p1.concave_conjugate().is_concave()
            True

        """
        return self._conjugate(True)     

    def is_convex(self):
        """
        Return "True" if the induced function of the tropical polynomial
        is convex
        """
        if len(self.dict()) == 1:
            return True
        
        return not self.parent().base()._use_min
    
    def is_concave(self):
        """
        Return "True" if the induced function of the tropical polynomial
        is concave
        """
        if len(self.dict()) == 1:
            return True
        
        return self.parent().base()._use_min
    
    def factor(self):
        r"""
        Return the factorization of ``self`` into its tropical linear factors

        OUTPUT:

        A Factorization object of ``self``

        EXAMPLES:

            sage: T = TropicalSemiring(QQ, use_min=True)
            sage: R = PolynomialRing(T, x)
            sage: p1 = R([6,3,1,0]); p1
            0*x^3 + 1*x^2 + 3*x + 6
            sage: factor(p1)
            (0) * (0*x + 1) * (0*x + 2) * (0*x + 3)

        Such factorization is not always possible::
        
            sage: p2 = R([4,4,2]); p2
            2*x^2 + 4*x + 4
            sage: p2.factor()
            (2) * (0*x^2 + 2*x + 2)

        """
        from sage.structure.factorization import Factorization

        if self.parent().base()._use_min:
            conjugate = self.concave_conjugate()
        else:
            conjugate = self.convex_conjugate()

        unit = self.dict()[self.degree()]
        if self != conjugate or self.roots() == []:
            factor = [(self*self.parent(-unit.lift()), 1)]
            return Factorization(factor, unit=unit)

        R = self.parent()
        roots_order = {}
        for root in self.roots():
            if root in roots_order:
                roots_order[root] += 1
            else:
                roots_order[root] = 1 
        factors = []
        for root in roots_order:
            factors.append((R([root, 0]), roots_order[root]))
        
        return Factorization(factors, unit=unit)

    def piecewise_function(self):
        r"""
        Return the tropical polynomial function of ``self`` which is a 
        piecewise linear function with the domains are split by roots

        OUTPUT:

        - ``f`` -- a piecewise function

        EXAMPLES:

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R = PolynomialRing(T, x)
            sage: p1 = R([4,2,1,3]); p1
            3*x^3 + 1*x^2 + 2*x + 4
            sage: p1.piecewise_function()
            piecewise(x|-->4 on (-oo, 1/3], x|-->3*x + 3 on (1/3, +oo); x)

        A constant tropical polynomial will result in a constant function::

            sage: p2 = R(3)
            sage: p2.piecewise_function()
            3

        A monomial will result in a linear function::

            sage: S.<x> = PolynomialRing(QQ)
            sage: p3 = R(x^3)
            sage: p3.piecewise_function()
            3*x + 1
            
        """
        from sage.symbolic.ring import SR
        from sage.rings.infinity import infinity
        from sage.functions.piecewise import piecewise

        x = SR.var('x')
        R = self.parent().base().base_ring()
        if self.roots() == []:
            # f = R(str(self.dict()[0]))
            f = self.dict()[0].lift()
            return f
        
        if len(self.dict()) == 1:
            gradient = list(self.dict())[0]
            intercept = self.dict()[gradient].lift()
            f = intercept+gradient*x
            return f

        unique_root = sorted(list(set(self.roots())))
        pieces = []
        domain = []
        for i in range(len(unique_root)+1):
            if i == 0:
                test_number = self.base_ring()(unique_root[i]-1)
            elif i == len(unique_root):
                test_number = self.base_ring()(unique_root[i-1]+1)
            else:
                test_number = self.base_ring()((unique_root[i] + \
                                                unique_root[i-1])/2)

            terms = {i:c*(test_number**i) for i, c in self.dict().items()}
            if self.base_ring()._use_min:
                maximum = min(terms.values())
            else:
                maximum = max(terms.values())
            found_key = None
            for key, value in terms.items():
                if value == maximum:
                    found_key = key
                    break
            gradient = found_key
            intercept = self.dict()[found_key].lift()

            if i == 0:
                interval = RealSet.unbounded_below_closed(unique_root[i])
                piecewise_linear = (interval, intercept+gradient*x)
                domain.append(interval)
            elif i == len(unique_root):
                if domain[0][0].upper_closed():
                    interval = RealSet.unbounded_above_open(unique_root[i-1])
                else:
                    interval = RealSet.unbounded_above_closed(unique_root[i-1])
                piecewise_linear = (interval, intercept+gradient*x)
                domain.append(interval)
            else:
                if domain[i-1][0].upper_closed():
                    interval = RealSet((unique_root[i-1], unique_root[i]))
                else:
                    interval = RealSet([unique_root[i-1], unique_root[i]])
                piecewise_linear = (interval, intercept+gradient*x)
                domain.append(interval)
            pieces.append(piecewise_linear)

        f = piecewise(pieces)
        return f
    
    def plot(self, xmin=None, xmax=None):
        r"""
        Return the plot of tropical polynomial function of ``self``

        INPUT:

        - ``xmin`` -- (default: ``None``)
        - ``xmax`` -- (default: ``None``)

        OUTPUT:

        If ``xmin`` and ``xmax`` is given, then it will return a plot
        of piecewise linear function of ``self`` with domain start from
        ``xmin`` to ``xmax``. Otherwise, the domain will start from the
        the minimum root of ``self`` minus 1 to maximum root of ``self``
        plus 1. If the function of ``self`` is constant or linear, then 
        the default domain will be [-1,1].

        EXAMPLES:

        If the tropical semiring use a max-plus algebra, then the graph 
        will be of piecewise linear convex function::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R = PolynomialRing(T, x)
            sage: p1 = R([4,2,1,3]); p1
            3*x^3 + 1*x^2 + 2*x + 4
            sage: p1.roots()
            [1/3, 1/3, 1/3]
            sage: p1.plot()

        A different result will be obtained if the tropical semiring employs 
        a min-plus algebra. Rather, a graph of the piecewise linear concave 
        function will be obtained::

            sage: T = TropicalSemiring(QQ, use_min=True)
            sage: R = PolynomialRing(T, x)
            sage: p1 = R([4,2,1,3])
            sage: p1.roots()
            [-2, 1, 2]
            sage: plot(p1, xmin=-4, xmax=4)
        
        TESTS:

        If ``xmin`` or ``xmax`` is given as an input, then the others also
        have to be given. Otherwise it will raise an error::

            sage: plot(p1, 5)
            Traceback (most recent call last):
            ...
            ValueError: Expected 2 inputs for xmin and xmax, but got 1
        
        The error also occured when ``xmin`` is greater or equal than ``xmax``::

            sage: plot(p1, 5, 3)
            Traceback (most recent call last):
            ...
            ValueError: xmin = 5 should be less than xmax = 3

        """
        from sage.plot.plot import plot
        f = self.piecewise_function()
        if xmin is None and xmax is None:
            roots = sorted(self.roots())
            if roots==[] or self.parent().base().zero() in roots:
                return plot(f, xmin=-1, xmax=1)
            else:
                return plot(f, xmin=roots[0]-1, xmax=roots[-1]+1)
        elif xmin is None or xmax is None:
            raise ValueError(f"Expected 2 inputs for xmin and xmax, but got 1")
        elif (xmin>=xmax):
            raise ValueError(f"xmin = {xmin} should be less than xmax = {xmax}")
        else:
            return plot(f, xmin=xmin, xmax=xmax)


class TropicalPolynomialSemiring(UniqueRepresentation, Parent):
    """
    Semiring structure of tropical polynomials in one variable
    """

    @staticmethod
    def __classcall_private__(cls, base_semiring, names=None):
        if names is None:
           names = 'x'
        if isinstance(names, str):
            names = (names,)
        return super().__classcall__(cls, base_semiring, tuple(names))

    def __init__(self, base_semiring, names):
        Parent.__init__(self, base=base_semiring, names=names, category=Semirings())

    Element = TropicalPolynomial

    def _element_constructor_(self, x, check=True):
        """
        Convert ``x`` into this tropical polynomial semiring
        """
        C = self.element_class
        if isinstance(x, (list, tuple)):
            for i, coeff in enumerate(x):
                if coeff == 0:
                    x[i] = self.base()(0)
        return C(self, x, check=check)

    def _repr_(self):
        return (f"Univariate Tropical Polynomial Semiring in {self.variable_name()}"
            f" over {self.base_ring().base_ring()}")

    def gens(self):
        gens = []
        for v in self.variable_names():
            gen = SR.var(v)
            gens.append(gen)
        return tuple(gens)