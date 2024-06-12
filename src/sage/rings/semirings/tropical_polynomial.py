r"""
Univariate Tropical Polynomial Semirings

Implements the parent class and element class for sparse polynomials over tropical semirings. We can do addition and multiplication for tropical polynomials. Added a method to find the tropical roots of tropical polynomials. Overriden the plot method so it specifically graph the tropical polynomial in cartesian coordinates.

AUTHORS:

- Verrel Rievaldo Wijaya

- Travis Scrimshaw

EXAMPLES::

    We can create a tropical polynomial semiring by first defining a tropical semiring and then inputting it to ``PolynomialRing`` constructor::

        sage: T = TropicalSemiring(QQ, use_min=False)
        sage: R = PolynomialRing(T,'y')
        sage: R
        Tropical Polynomial Semiring in y over Rational Field

    We can define the element by giving a list of coefficients that starts from constant::

        sage: p1 = R([1,4,None,T(0)]); p1
        0*y^3 + 4*y + 1

    Or by converting from classical polynomial::

        sage: S.<y> = PolynomialRing(QQ)
        sage: p2 = R(y^2+2*y+3); p2
        y^2 + 2*y + 3

    We can do the addition for two tropical polynomials::

        sage: p1 + p2
        0*y^3 + y^2 + 4*y + 3

    We can do the multiplication for two tropical polynomials::

        sage: p1 * p2
        y^5 + 2*y^4 + 5*y^3 + 6*y^2 + 7*y + 4

    Beware that when multiplying tropical polynomial with scalar, it will give an error if the scalar is not tropical number::

        sage: 2 * p1
        Traceback (most recent call last):
        ...
        ArithmeticError: cannot negate any non-infinite element
        sage: T(2) * p1
        2*y^3 + 6*y + 3

    We can do the evaluation of tropical polynomial at some value::
        
        sage: p1(3)
        9

    We can find all the tropical roots of tropical polynomial counted with multiplicity. There will no tropical root for constant polynomial. For a monomial, the tropical root is `-infinity`::

        sage: p1.roots()
        [-3, 2, 2]
        sage: p3 = R(1)
        sage: p3.roots()
        []
        sage: p4 = R(y^2)
        sage: p4.roots()
        [-infinity, -infinity]

    To show the plot of tropical polynomial we use::
        sage: p1.plot()

TESTS::

    There is no subtraction for tropical polynomials because element in tropical semiring doesn't necessarily have additive inverse::

        sage: -p1
        Traceback (most recent call last):
        ...
        ArithmeticError: cannot negate any non-infinite element

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

from sage.rings.semirings.tropical_polynomial_element import TropicalPolynomial
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.rings.polynomial.polynomial_element_generic import Polynomial_generic_sparse
from itertools import combinations

class TropicalPolynomial(Polynomial_generic_sparse):
    """
    A generic sparse tropical polynomial.

    The `TropicalPolynomial`` class defines functionality for sparse
    polynomials over any tropical semiring. A sparse polynomial is represented using a
    dictionary which maps each exponent to the corresponding coefficient. The
    coefficients in this case is a tropical number.
    """

    def __call__(self, val):
        r"""
        Return the value of ``self`` evaluated at ``val``.

        INPUT:

        - ``val`` -- a number that should be inside the base ring of tropical semiring

        OUTPUT:

        A single tropical number

        EXAMPLES::
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
        Return the list of all tropical roots of ``self``.

        OUTPUT:

        A list containing tropical roots of ``self`` counted with multiplicity 

        EXAMPLES::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R = PolynomialRing(T, x)
            sage: S.<x> = PolynomialRing(QQ)
            sage: p1 = R(x^3-2*x^2+3*x-1)
            sage: p1.roots()
            [-4, 1, 1]
            sage: p2 = TPS([T(0), None, T(0)]); p2
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
        dict_root = {i:R(str(c)) for i,c in self.dict().items()}
        for comb in combinations(dict_root.keys(), 2):
            index1, index2 = comb[0], comb[1]
            root = (dict_root[index1]-dict_root[index2])/(index2 - index1)
            val_root = dict_root[index1] + index1*root
            check_maks = True
            for key in dict_root.keys():
                if key not in comb:
                    val = dict_root[key] + key*root
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
                tropical_roots += [root] * order

                if len(tropical_roots) == self.degree():
                    break
            
        return tropical_roots
    
    def plot(self):
        r"""
        Return the graphs of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R = PolynomialRing(T, x)
            sage: p1 = R([T(0),1,-1]); p1
            -x^2 + x + 0
            sage: p1.roots()
            [-1, 2]
            sage: p1.plot()
            
        """
        from sage.all import var, plot, parametric_plot

        t = var('t')
        R = self.parent().base().base_ring()
        if self.roots() == []:
            return plot(R(str(self.dict()[0])))
        
        if len(self.dict()) == 1:
            gradient = list(self.dict().keys())[0]
            intercept = R(str(self.dict()[gradient]))
            return plot(intercept+gradient*t)

        unique_root = sorted(list(set(self.roots())))
        all_plot = []
        for i in range(len(unique_root)+1):
            if i == 0:
                test_number = self.base_ring()(unique_root[i]-1)
            elif i == len(unique_root):
                test_number = self.base_ring()(unique_root[i-1]+1)
            else:
                test_number = self.base_ring()((unique_root[i]+unique_root[i-1])/2)

            terms = {i:c*(test_number**i) for i, c in self.dict().items()}
            maximum = max(terms.values())
            found_key = None
            for key, value in terms.items():
                if value == maximum:
                    found_key = key
                    break
            gradient = found_key
            intercept = R(str(self.dict()[found_key]))

            if i == 0:
                piecewise_linear = parametric_plot((t, intercept+gradient*t), (t, unique_root[i]-1, unique_root[i]))
            elif i == len(unique_root):
                piecewise_linear = parametric_plot((t, intercept+gradient*t), (t, unique_root[i-1], unique_root[i-1]+1))
            else:
                piecewise_linear = parametric_plot((t, intercept+gradient*t), (t, unique_root[i-1], unique_root[i]))

            all_plot.append(piecewise_linear)
        
        return sum(all_plot)

class TropicalPolynomialSemiring(UniqueRepresentation, Parent):
    """
    Semiring structure of tropical polynomials in one variable.    
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
        return C(self, x, check=check)

    def _repr_(self):
        return f"Tropical Polynomial Semiring in {self.variable_name()} over {self.base_ring().base_ring()}"