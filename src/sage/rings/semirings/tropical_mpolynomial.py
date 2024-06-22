r"""
Multivariate Tropical Polynomial Semirings

<Description>

AUTHORS:

- Verrel Rievaldo Wijaya

EXAMPLES:

Construct multivariate tropical polynomial semirings::

    sage: T = TropicalSemiring(QQ, use_min=False)
    sage: R = PolynomialRing(T, 'a,b')
    sage: R
    Multivarite Tropical Polynomial Semiring in a, b over Rational Field
    
Create an element by inputting a dictionary::

    sage: dict1 = {(1,0):0, (0,1):-1, (1,1):3}
    sage: p1 = R(dict1); p1
    3*a*b + 0*a + (-1)*b

We can also create an element by converting from classical polynomial::

    sage: S.<a,b> = QQ[]
    sage: f = a + b + a*b
    sage: p2 = R(f); p2
    1*a*b + 1*a + 1*b

Some basic arithmetic operations::

    sage: p1 + p2
    3*a*b + 1*a + 1*b
    sage: p1 * p2
    4*a^2*b^2 + 4*a^2*b + 1*a^2 + 4*a*b^2 + 1*a*b + 0*b^2
    sage: T(2) * p1
    5*a*b + 2*a + 1*b

TESTS:

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

from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.multi_polynomial_element import MPolynomial_polydict
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.polynomial.polydict import PolyDict, ETuple

class TropicalCurve(Parent):
    def __init__(self, *args):
        pass

    def plot(self):
        pass


class TropicalMPolynomial(MPolynomial_polydict):
    r"""
    Generic multivariate tropical polynomial.

    """

    def add_zero_coefficient(self, *args):
        """
        Return a new tropical polynomial that include additional/modified
        terms with coefficient equals 0
        
        INPUT:

        *args: A variable number of arguments. Each argument should be of 
        type tuple that represent the exponent of terms

        EXAMPLES:

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R = PolynomialRing(T, 'x,y')
            sage: S.<x,y> = QQ[]
            sage: c1 = 5+3*x+2*y+x*y
            sage: p1 = R(c1); p1
            1*x*y + 3*x + 2*y + 5
            sage: p1.add_zero_coefficient((0,0),(3,0))
            0*x^3 + 1*x*y + 3*x + 2*y + 0

        TESTS:

            sage: p1.add_zero_coefficient((1,2,1))
            Traceback (most recent call last):
            ...
            ValueError: Lenght of tuple should not exceed 2

        """
        new_dict = self.dict()
        R = self.parent()
        for exponent in args:
            if len(exponent) > len(R.variable_names()):
                raise ValueError(f"Lenght of tuple should not exceed " \
                                 f"{len(R.variable_names())}")
            etuple = ETuple(exponent)
            new_dict[etuple] = 0

        return R(new_dict)
    
    def tropical_curve(self):
        r"""
        Return tropical roots

        OUTPUT:

        - tropical_roots -- a list of lists, where the inner list is of the
        form [point, condition for parameter, order]
        
        EXAMPLES:

        Some examples for tropical polynomials in two variables::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R = PolynomialRing(T, 'x,y')
            sage: dict1 = {(0,0):0, (1,0):0, (0,1):0}
            sage: p1 = R(dict1); p1
            0*x + 0*y + 0
            sage: p1.tropical_curve()
            [[(0, r1), [r1 < 0], 'order = 1'],
            [(r2, 0), [r2 < 0], 'order = 1'],
            [(r3, r3), [r3 > 0], 'order = 1']]

        ::

            sage: S.<x,y> = QQ[]
            sage: c1 = 3+2*x+2*y+3*x*y
            sage: dict2 = {(2,0):0, (0,2):0}
            sage: p2 = R(c1) + R(dict2); p2
            0*x^2 + 3*x*y + 2*x + 0*y^2 + 2*y + 3
            sage: p2.tropical_curve()
            [[(r4, -1), [1 < r4, r4 < 2], 'order = 1'],
            [(-1, r5), [1 < r5, r5 < 2], 'order = 1'],
            [(-r6, r6), [-1 < r6, r6 < 1], 'order = 1'],
            [(r7 + 3, r7), [-1 < r7], 'order = 1'],
            [(r8 - 3, r8), [2 < r8], 'order = 1'],
            [(1, r10), [r10 < -1], 'order = 1'],
            [(2, r11), [r11 < -1], 'order = 1'],
            [(r13, 1), [r13 < -1], 'order = 1'],
            [(r15, 2), [r15 < -1], 'order = 1']]

        ::

            sage: c2 = -1*x^2
            sage: dict3 = {(0,0):0, (1,0):0, (0,2):0}
            sage: p3 = R(c2) + R(dict3); p3
            (-1)*x^2 + 0*x + 0*y^2 + 0
            sage: p3.tropical_curve()
            [[(0, r34), [r34 < 0], 'order = 1'],
            [(r35, 0), [r35 < 0], 'order = 2'],
            [(2*r37, r37), [0 < r37, r37 < (1/2)], 'order = 1'],
            [(1, r38), [r38 < (1/2)], 'order = 1'],
            [(r39 + 1/2, r39), [(1/2) < r39], 'order = 2']]

        We can find tropical curve for any tropical polynomials in 
        `n\geq 2` variables:

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R = PolynomialRing(T, 'x,y,z')
            sage: S.<x,y,z> = QQ[]
            sage: p1 = R(x*y + (-1/2)*x*z + 4*z^2); p1
            1*x*y + (-1/2)*x*z + 4*z^2
            sage: p1.tropical_curve()
            [[(r32, r31 - 3/2, r31), [r31 + 9/2 < r32], 'order = 1'],
            [(2*r33 - r34 + 3, r34, r33), [r33 < r34 + 3/2], 'order = 1'],
            [(r35 + 9/2, r36, r35), [r36 + 3/2 < r35], 'order = 1']]

        """
        from sage.symbolic.ring import SR
        from sage.symbolic.relation import solve
        from itertools import combinations
        from sage.arith.misc import gcd

        tropical_roots = []
        variables = []
        for name in self.parent().variable_names():
            variables.append(SR.var(name))

        # convert each term to its linear function
        linear_eq = {}
        for key in self.dict():
            eq = 0
            for i,e in enumerate(key):
                eq += variables[i]*e
            eq += self.dict()[key].lift()
            linear_eq[key] = eq

        # checking for all possible combinations of two terms
        for keys in combinations(self.dict(), 2):
            sol = solve(linear_eq[keys[0]]==linear_eq[keys[1]], variables)
            
            # parametric solution of these two terms
            final_sol = []
            for s in sol[0]:
                final_sol.append(s.right())
            xy_interval = []
            xy_interval.append(tuple(final_sol))
            
            # comparing with other terms
            min_max = linear_eq[keys[0]]
            for i,v in enumerate(variables):
                min_max = min_max.subs(v==final_sol[i])
            
            all_sol_compare = []
            no_solution = False
            for compare in self.dict():
                if compare not in keys:
                    temp_compare = linear_eq[compare]
                    for i, v in enumerate(variables):
                        temp_compare = temp_compare.subs(v==final_sol[i])
                    if self.parent().base()._use_min:
                        sol_compare = solve(min_max < temp_compare, variables)
                    else:
                        sol_compare = solve(min_max > temp_compare, variables)
                    if sol_compare: # if there is solution
                        if isinstance(sol_compare[0], list):
                            if sol_compare[0]:
                                all_sol_compare.append(sol_compare[0][0])
                        else: # solution is unbounded on one side
                            all_sol_compare.append(sol_compare[0])
                    else:
                        no_solution = True
                        break

            # solve the condition for parameter
            if not no_solution:
                parameter = set()
                for sol in all_sol_compare:
                    parameter = parameter.union(set(sol.variables()))
                parameter_solution = solve(all_sol_compare, list(parameter))
                
                if parameter_solution:
                    xy_interval.append(parameter_solution[0])
                    # calculate order
                    index_diff = []
                    for i in range(len(keys[0])):
                        index_diff.append(abs(keys[0][i]-keys[1][i]))
                    order = gcd(index_diff)
                    xy_interval.append(f"order = {order}")

                    tropical_roots.append(xy_interval)

        return tropical_roots   


class TropicalMPolynomialSemiring(UniqueRepresentation, Parent):
    """
    Semiring structure of tropical polynomials in multiple variables
    """

    def __init__(self, base_semiring, names):
        Parent.__init__(self, base=base_semiring, names=names)

    Element = TropicalMPolynomial

    def _element_constructor_(self, x):
        """"
        Construct the element of this tropical multivariate polynomial
        semiring

        INPUT:

        - ``x`` -- dict or MPolynomial

        """

        C = self.element_class
        new_dict = {}

        if isinstance(x, MPolynomial):
            x = x.dict()

        for key, value in x.items(): # convert each coefficient to tropical
            new_dict[key] = self.base()(value)

        return C(self, new_dict)
    
    def _repr_(self):
        return (f"Multivarite Tropical Polynomial Semiring in {', '.join(self.variable_names())}"
            f" over {self.base_ring().base_ring()}")
    
    def random_element(self):
        """
        Return a random element from this semiring

        """

        from sage.rings.polynomial.polynomial_ring_constructor import \
            PolynomialRing
        R = PolynomialRing(self.base().base_ring(), self.variable_names())

        return self(R.random_element())
        