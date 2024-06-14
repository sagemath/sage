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

- Travis Scrimshaw

EXAMPLES:

    Construct a tropical polynomial semiring by first defining a base 
    semiring and then inputting it to ``PolynomialRing`` constructor::

        sage: T = TropicalSemiring(QQ, use_min=False)
        sage: R = PolynomialRing(T,'y')
        sage: R
        Tropical Polynomial Semiring in y over Rational Field

    We can define the element by giving a list of coefficients that
    begins with constant. This is also the way to create a tropical
    polynomial with `0` as coefficient::

        sage: p1 = R([1,4,None,0]); p1
        0*y^3 + 4*y + 1

    Create an element by converting from classical polynomial::

        sage: S.<y> = PolynomialRing(QQ)
        sage: p2 = R(y^2+2*y+3); p2
        y^2 + 2*y + 3

    We can do the addition, multiplication, and evaluation for tropical 
    polynomials::

        sage: p1 + p2
        0*y^3 + y^2 + 4*y + 3
        sage: p1 * p2
        y^5 + 2*y^4 + 5*y^3 + 6*y^2 + 7*y + 4
        sage: p1(3)
        9

    Beware that when multiplying tropical polynomial with a scalar, it
    will give an error if the scalar is not tropical number::

        sage: 2 * p1
        Traceback (most recent call last):
        ...
        ArithmeticError: cannot negate any non-infinite element
        sage: T(2) * p1
        2*y^3 + 6*y + 3 

    We can find all the tropical roots of tropical polynomial counted
    with multiplicity. There will be no tropical root for constant
    polynomial. For a monomial, the tropical root is the additive identity
    of its base tropical semiring::

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

    The factorization of tropical polynomial to its linear factors::

        sage: p1.factorization()
        '(y - 3)*(y + 2)^2'
        sage: p2.factorization()
        '1*(y + 1)^2'
    
    To show the induced tropical polynomial function which is a piecewise 
    linear function::

        sage: p1.piecewise_function()
        piecewise(x|-->1 on (-oo, -3), x|-->x + 4 on [-3, 2], x|-->3*x on 
        (2, +oo); x)

    Draw the graph of tropical polynomial::
        sage: p1.plot()
        sage: plot(p2)

TESTS:

    There is no subtraction for tropical polynomials because element
    in tropical semiring doesn't necessarily have additive inverse::

        sage: -p1
        Traceback (most recent call last):
        ...
        ArithmeticError: cannot negate any non-infinite element

REFERENCES:

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

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.rings.polynomial.polynomial_element_generic import \
Polynomial_generic_sparse
from itertools import combinations

class TropicalPolynomial(Polynomial_generic_sparse):
    """
    A generic sparse tropical polynomial.

    The `TropicalPolynomial`` class defines functionality for sparse
    polynomials over any tropical semiring. A sparse polynomial is 
    represented using a dictionary which maps each exponent to the
    corresponding coefficient. The coefficients is a tropical number.
    """

    def __call__(self, val):
        r"""
        Return the value of ``self`` evaluated at ``val``.

        INPUT:

        - ``val`` -- a number from the base ring of tropical semiring

        OUTPUT:

        A single tropical number

        EXAMPLES:

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R = PolynomialRing(T, x)
            sage: p1 = R([1,4,None,T(0)])
            sage: p1(1)
            5

        """

        if val.parent() is not self.base_ring():
            val = self.base_ring()(val)
        terms = [c*(val**i) for i, c in self.dict().items()]
        if self.base_ring()._use_min:
            return min(terms)
        else:
            return max(terms)
    
    def roots(self):
        r"""
        Return the list of all tropical roots of ``self``

        OUTPUT:

        - ``tropical_roots`` -- A list containing tropical roots of ``self``
        counted with multiplicity

        EXAMPLES:

            sage: T = TropicalSemiring(QQ, use_min=True)
            sage: R = PolynomialRing(T, x)
            sage: p1 = R([5,4,1,0,2,4,3]); p1
            3*x^6 + 4*x^5 + 2*x^4 + 0*x^3 + x^2 + 4*x + 5
            sage: p1.roots()
            [2, 2, 1, -1, -1, -1]
            sage: p2 = R([0, None, 0]); p2
            0*x^2 + 0
            sage: p2.roots()
            [0, 0]

        """

        tropical_roots = []
        if len(self.dict()) == 1:
            exponent = list(self.dict().keys())[0]
            if exponent == 0:
                return tropical_roots
            else:
                return [self.parent().base_ring().zero()]*exponent
        
        R = self.parent().base().base_ring()
        dict_root = {}
        dict_coeff = {i:R(str(c)) for i,c in self.dict().items()}
        for comb in combinations(dict_coeff.keys(), 2):
            index1, index2 = comb[0], comb[1]
            root = (dict_coeff[index1]-dict_coeff[index2])/(index2 - index1)
            val_root = dict_coeff[index1] + index1*root
            check_maks = True
            for key in dict_coeff.keys():
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
            
        return tropical_roots
    
    def conjugate(self):
        r"""
        Return the tropical polynomial which has the same graph as ``self``
        but which can be factored

        OUTPUT:

        A tropical polynomial element

        EXAMPLES:

        There are tropical polynomials which is its own conjugate::

            sage: T = TropicalSemiring(QQ, use_min=True)
            sage: R = PolynomialRing(T, x)
            sage: p1 = R([5,3,1]); p1
            x^2 + 3*x + 5
            sage: p1.conjugate()
            x^2 + 3*x + 5
        
        Checking the plot of tropical polynomial and its conjugate::

            sage: p2 = R([4,4,2]); p2
            2*x^2 + 4*x + 4
            sage: p3 = p2.conjugate(); p3
            2*x^2 + 3*x + 4
            sage: p2.plot()
            sage: p3.plot()

        """

        roots = self.roots()
        poly = self.dict()[self.degree()]
        for root in roots:
            linear = self.parent()([root, 0])
            poly *= linear
        return poly
    
    def factorization(self):
        r"""
        Return the factorization of ``self`` to its linear factors

        OUTPUT:

        If ``self`` can be factored, then return a string that represents 
        its factorization. Otherwise, return the string 'Can not be factored'

        EXAMPLES:

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R = PolynomialRing(T, x)
            sage: S.<x> = PolynomialRing(QQ)
            sage: p1 = R(x^3+4*x^2+5*x-1)
            sage: p1.roots()
            [-6, 1, 3]
            sage: p1.factorization()
            '1*(x + (-6))*(x + 1)*(x + 3)'
        
        A tropical polynomial can not be factored if it is not the same with
        its conjugate::

            sage: p2 = R([2,1,-1]); p2
            -x^2 + x + 2
            sage: p2.factorization()
            'Can not be factored'
            sage: p2.conjugate()
            -x^2 + 1/2*x + 2

        The addition in tropical polynomial is defined in the tropical sense.
        Therefore the monomial will be factored to::

            sage: p3 = R(3*x^4)
            sage: p3.factorization()
            '3*(x + +infinity)^4'

        """

        if self != self.conjugate():
            return 'Can not be factored'

        if self.roots() == []:
            return str(self)
        
        leading = str(self.dict()[self.degree()])
        variable = self.variable_name()
        terms = []
        if leading != '0':
            terms.append(leading)
        roots_order = {}
        for root in self.roots():
            if root in roots_order:
                roots_order[root] += 1
            else:
                roots_order[root] = 1
        for root, order in roots_order.items():
            if root >= 0:
                if order == 1:
                    term = f"({variable} + {root})"
                else:
                    term = f"({variable} + {root})^{order}"
            else:
                if order == 1:
                    term = f"({variable} + ({root}))"
                else:
                    term = f"({variable} + ({root}))^{order}"
            terms.append(term)
        
        return '*'.join(terms)

    def piecewise_function(self):
        r"""
        Return the tropical polynomial function of ``self`` which is a 
        piecewise linear function with the domains are split by roots

        OUTPUT:

        - ``f`` -- a piecewise function in single variable

        EXAMPLES:

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R = PolynomialRing(T, x)
            sage: p1 = R([4,1,0]); p1
            3*x^3 + x^2 + 2*x + 4
            sage: p1.piecewise_function()
            piecewise(x|-->2*x on (-oo, 1), x|-->x + 1 on (1, 3), x|-->4 on 
            (3, +oo); x)

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

        from sage.all import var, piecewise, infinity

        x = var('x')
        R = self.parent().base().base_ring()
        if self.roots() == []:
            f = R(str(self.dict()[0]))
            return f
        
        if len(self.dict()) == 1:
            gradient = list(self.dict().keys())[0]
            intercept = R(str(self.dict()[gradient]))
            f = intercept+gradient*x
            return f

        unique_root = sorted(list(set(self.roots())))
        pieces = []
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
            intercept = R(str(self.dict()[found_key]))

            if i == 0:
                piecewise_linear = ((-infinity, unique_root[i]), \
                                    intercept+gradient*x)
            elif i == len(unique_root):
                piecewise_linear = ((unique_root[i-1], infinity), \
                                    intercept+gradient*x)
            else:
                piecewise_linear = ((unique_root[i-1], unique_root[i]), \
                                    intercept+gradient*x)

            pieces.append(piecewise_linear)

        f = piecewise(pieces)
        
        return f
    
    def plot(self):
        r"""
        Return the plot of tropical polynomial function of ``self``

        OUTPUT:

        A plot of piecewise linear function with its `xmin` and `xmax` 
        modified so the graph looks clear and well-presented

        EXAMPLES:

        If the tropical semiring use a max-plus algebra, then the graph 
        will be of piecewise linear convex function::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R = PolynomialRing(T, x)
            sage: p1 = R([4,2,1,3]); p1
            3*x^3 + x^2 + 2*x + 4
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
            [2, 1, -2]
            sage: plot(p1)
            
        """

        from sage.all import plot
        f = self.piecewise_function()
        roots = sorted(self.roots())
        return plot(f, (x, roots[0]-1, roots[-1]+1))


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
        Parent.__init__(self, base=base_semiring, names=names)

    Element = TropicalPolynomial

    def _element_constructor_(self, x, check=True):
        C = self.element_class
        if isinstance(x, list):
            for i, coeff in enumerate(x):
                if coeff == 0:
                    x[i] = self.base()(0)
        return C(self, x, check=check)

    def _repr_(self):
        return (f"Tropical Polynomial Semiring in {self.variable_name()}"
            f" over {self.base_ring().base_ring()}")