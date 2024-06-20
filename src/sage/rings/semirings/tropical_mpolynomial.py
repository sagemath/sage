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

class TropicalMPolynomial(MPolynomial_polydict):
    
    def roots(self):
        """

        OUTPUT:

        - tropical_roots -- a list of list, where the inner list is of the
        form [[x0,y0], [x1,y1], gradient, order] with [x0, y0] and [x1,y1]
        is the coordinates of point defining the line segment of tropical
        roots (two variables only)
        
        EXAMPLES:

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R = PolynomialRing(T, 'x,y')
            sage: dict1 = {(0,0):0, (1,0):0, (0,1):0}
            sage: p1 = R(dict1); p1
            0*x + 0*y + 0
            sage: p1.roots()
            [[[0, -Infinity], [0, 0], 'm = +Infinity', 'order = 1'],
            [[-Infinity, 0], [0, 0], 'm = 0', 'order = 1'],
            [[0, 0], [+Infinity, +Infinity], 'm = 1', 'order = 1']]

        ::

            sage: S.<x,y> = QQ[]
            sage: c1 = 3+2*x+2*y+3*x*y
            sage: dict2 = {(2,0):0, (0,2):0}
            sage: p2 = R(c1) + R(dict2); p2
            0*x^2 + 3*x*y + 2*x + 0*y^2 + 2*y + 3
            sage: p2.roots()
            [[[1, -1], [2, -1], 'm = 0', 'order = 1'],
            [[-1, 1], [-1, 2], 'm = +Infinity', 'order = 1'],
            [[1, -1], [-1, 1], 'm = -1', 'order = 1'],
            [[2, -1], [+Infinity, +Infinity], 'm = 1', 'order = 1'],
            [[-1, 2], [+Infinity, +Infinity], 'm = 1', 'order = 1'],
            [[1, -Infinity], [1, -1], 'm = +Infinity', 'order = 1'],
            [[2, -Infinity], [2, -1], 'm = +Infinity', 'order = 1'],
            [[-Infinity, 1], [-1, 1], 'm = 0', 'order = 1'],
            [[-Infinity, 2], [-1, 2], 'm = 0', 'order = 1']]

        ::

            sage: c2 = -1*x^2
            sage: dict3 = {(0,0):0, (1,0):0, (0,2):0}
            sage: p3 = R(c2) + R(dict3); p3
            (-1)*x^2 + 0*x + 0*y^2 + 0
            sage: p3.roots()
            [[[0, -Infinity], [0, 0], 'm = +Infinity', 'order = 1'],
            [[-Infinity, 0], [0, 0], 'm = 0', 'order = 2'],
            [[0, 0], [1, 1/2], 'm = 1/2', 'order = 1'],
            [[1, -Infinity], [1, 1/2], 'm = +Infinity', 'order = 1'],
            [[1, 1/2], [+Infinity, +Infinity], 'm = 1', 'order = 2']]

        """
        from sage.symbolic.ring import SR
        from sage.symbolic.relation import solve
        from sage.sets.real_set import RealSet
        from sage.rings.infinity import infinity
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

        # checking for all combinations of two terms
        for keys in combinations(self.dict(), 2):
            sol = solve(linear_eq[keys[0]]==linear_eq[keys[1]], variables)
            
            # parametric solution of these two terms
            final_sol = []
            for s in sol[0]:
                final_sol.append(s.right())
            
            # cheking if it is maximum with other terms
            temp_max = linear_eq[keys[0]]
            for i,v in enumerate(variables):
                temp_max = temp_max.subs(v==final_sol[i])
            sol_interval = RealSet().real_line()
            for compare in self.dict():
                if compare not in keys:
                    temp_compare = linear_eq[compare]
                    for i,v in enumerate(variables):
                        temp_compare = temp_compare.subs(v==final_sol[i])
                    sol_compare = solve(temp_max>=temp_compare, variables)
                    if len(sol_compare)==0: # no solution
                        compare_interval = RealSet()
                    elif isinstance(sol_compare[0], list):
                        if not sol_compare[0]: # all of real line
                            compare_interval = RealSet().real_line()
                        else:
                            compare_interval = RealSet(sol_compare[0][0])
                    else: # solution is unbounded in one side
                        compare_interval = RealSet(sol_compare[0])
                    sol_interval = sol_interval.intersection(compare_interval)

            # if it is really the maximum, then find two points that define 
            # the line segment
            if sol_interval:
                lowerpoint = sol_interval[0].lower()
                upperpoint = sol_interval[0].upper()
                xy_interval = [[], []]
                check_numeric = False
                for i, s in enumerate(final_sol):
                    if not s.is_numeric():
                        variable = s.variables()[0] # take the parameter
                        f_interval = []
                        f_lower = s.subs(variable==lowerpoint)
                        f_interval.append(f_lower)
                        f_upper = s.subs(variable==upperpoint)
                        f_interval.append(f_upper)
                    else:
                        check_numeric = True
                        index_numeric = i
                        f_interval = [s,s]
                    for j, point in enumerate(f_interval):
                        xy_interval[j].append(point)

                # add gradient
                if check_numeric:
                    if index_numeric==0:
                        gradient = infinity
                    else:
                        gradient = 0
                else:
                    variable = final_sol[0].variables()[0]
                    gradient = 1/final_sol[0].coefficient(variable, 1)
                xy_interval.append(f"m = {gradient}")
                
                # calculate order
                order = gcd(abs(keys[0][0]-keys[1][0]), abs(keys[0][1]-keys[1][1]))
                xy_interval.append(f"order = {order}")

                tropical_roots.append(xy_interval)

        return tropical_roots   

class TropicalMPolynomialSemiring(UniqueRepresentation, Parent):
    def __init__(self, base_semiring, names):
        Parent.__init__(self, base=base_semiring, names=names)

    Element = TropicalMPolynomial

    def _element_constructor_(self, x):
        C = self.element_class
        new_dict = {}

        if isinstance(x, MPolynomial):
            x = x.dict()

        for key, value in x.items():
            new_dict[key] = self.base()(value)

        return C(self, new_dict)
    
    def _repr_(self):
        return (f"Multivarite Tropical Polynomial Semiring in {', '.join(self.variable_names())}"
            f" over {self.base_ring().base_ring()}")