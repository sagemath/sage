r"""
Base class for tropical variety objects

The TropicalVariety represents the solution set of a system of tropical 
polynomial equations. In tropical geometry, a tropical variety is defined by 
taking the minimum (maximum) of polynomials and considering the loci where 
these minima (maxima) are attained multiple times. This module provides 
functionality to compute, visualize, and analyze tropical varieties.

AUTHORS:

- Verrel Rievaldo Wijaya

EXAMPLES:

    sage: T = TropicalSemiring(QQ, use_min=False)
    sage: R = PolynomialRing(T, 'x,y')
    sage: dict1 = {(0,0):0, (1,0):0, (0,1):0}
    sage: p1 = R(dict1)
    sage: th1 = p1.tropical_hypersurface(); th1
    Tropical curve of 0*x + 0*y + 0 are 
    [[(0, t1), [t1 <= 0], 1]
    [(t1, 0), [t1 <= 0], 1]
    [(t1, t1), [t1 >= 0], 1]]
    sage: th1.dimension()
    2
    sage: th1.number_of_components()
    3
    sage: plot(th1)

::

    sage: S.<x,y> = QQ[]
    sage: c2 = -1*x^2
    sage: dict2 = {(0,0):0, (1,0):0, (0,2):0}
    sage: p2 = R(c2) + R(dict2)
    sage: th2 = p2.tropical_hypersurface(); th2
    Tropical curve of (-1)*x^2 + 0*x + 0*y^2 + 0 are 
    [[(0, t1), [t1 <= 0], 1]
    [(t1, 0), [t1 <= 0], 2]
    [(2*t1, t1), [0 <= t1, t1 <= (1/2)], 1]
    [(1, t1), [t1 <= (1/2)], 1]
    [(t1 + 1/2, t1), [(1/2) <= t1], 2]]
    sage: th2.plot()

::

    sage: T = TropicalSemiring(QQ, use_min=True)
    sage: R = PolynomialRing(T, 'a,x,y,z')
    sage: S.<a,x,y,z> = QQ[]
    sage: p1 = R(x*y + (-1/2)*x*z + 4*z^2 + a*x)
    sage: p1.tropical_hypersurface()
    Tropical hypersurface of 1*a*x + 1*x*y + (-1/2)*x*z + 4*z^2 are 
    [[(t1, t2, t1, t3), [max(t1 + 3/2, 1/2*t1 + 1/2*t2 - 3/2) <= t3], 1]
    [(t1 - 3/2, t2, t3, t1), [t2 - 9/2 <= t1, t1 <= t3 + 3/2, t2 - 6 <= t3], 1]
    [(2*t1 - t2 + 3, t2, t3, t1), [t1 <= min(1/2*t2 + 1/2*t3 - 3/2, t2 - 9/2)], 1]
    [(t1, t2, t3 - 3/2, t3), [t2 - 9/2 <= t3, t3 <= t1 + 3/2, t2 - 6 <= t1], 1]
    [(t1, 2*t2 - t3 + 3, t3, t2), [t3 + 3/2 <= t2, t3 <= t1], 1]
    [(t1, t2 + 9/2, t3, t2), [t2 <= min(t3 + 3/2, t1 + 3/2)], 1]]

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

from sage.structure.sage_object import SageObject
from sage.rings.rational_field import QQ
from sage.plot.graphics import Graphics
from sage.plot.plot import parametric_plot
from sage.rings.infinity import infinity
from sage.symbolic.ring import SR
import operator

class TropicalVariety(SageObject):
    
    def __init__(self, poly):
        from sage.symbolic.relation import solve
        from itertools import combinations
        from sage.arith.misc import gcd

        tropical_roots = []
        variables = []
        for name in poly.parent().variable_names():
            variables.append(SR.var(name))

        # convert each term to its linear function
        linear_eq = {}
        for key in poly.dict():
            eq = 0
            for i,e in enumerate(key):
                eq += variables[i]*e
            eq += poly.dict()[key].lift()
            linear_eq[key] = eq

        # checking for all possible combinations of two terms
        for keys in combinations(poly.dict(), 2):
            sol = solve(linear_eq[keys[0]]==linear_eq[keys[1]], variables)
            
            # parametric solution of the chosen two terms
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
            for compare in poly.dict():
                if compare not in keys:
                    temp_compare = linear_eq[compare]
                    for i, v in enumerate(variables):
                        temp_compare = temp_compare.subs(v==final_sol[i])
                    if poly.parent().base()._use_min:
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
                    xy_interval.append(order)
                    tropical_roots.append(xy_interval)

        self.poly = poly
        self._hypersurface = []
        dim_param = len(tropical_roots[0][0]) - 1
        vars = [SR.var('t{}'.format(i)) for i in range(1, dim_param+1)]
        for arg in tropical_roots:
            if not isinstance(arg, list):
                raise ValueError("Input must be a list")
            else:
                subs_dict = {}
                index_vars = 0
                new_eq = []
                for eq in arg[0]:
                    var_eq = eq.variables()
                    for var in var_eq:
                        if var not in subs_dict:
                            subs_dict[var] = vars[index_vars]
                            index_vars += 1
                    new_eq.append(eq.subs(subs_dict))
                new_eq = tuple(new_eq)
                arg.remove(arg[0])
                arg.insert(0, new_eq)

                params = arg[1]
                arg.remove(params)
                new_param = []
                for param in params:
                    lhs = param.lhs().subs(subs_dict)
                    rhs = param.rhs().subs(subs_dict)
                    if param.operator() == operator.gt:
                        expr = lhs >= rhs
                    else:
                        expr = lhs <= rhs
                    new_param.append(expr)
                arg.insert(1, new_param)

            self._hypersurface.append(arg)   

    def dimension(self):
        return len(self._hypersurface[0][0])
    
    def number_of_components(self):
        return len(self._hypersurface)
    
    def _repr_(self):
        components = "\n".join([f"{row}" for row in self._hypersurface])
        return (f"Tropical hypersurface of {self.poly} are \n[{components}]")

class TropicalCurve(TropicalVariety):
    r""""
    Represent a tropical curve in `\mathbb{R}^2`. The representation is in 
    the form of list of lists, where the inner list represent each line 
    segments of tropical roots.

    """

    def __init__(self, poly):
        TropicalVariety.__init__(self, poly)        
    
    def _axes(self):
        """
        Set the default axes for ``self``

        OUTPUT: A list of two lists, where the first inner list represent
        value of x-axis and the second inner list represent value of y-axis

        """
        xmin = xmax = list(self.vertex())[0][0]
        for vertice in self.vertex():
            if vertice[0] < xmin:
                xmin = vertice[0]
            elif vertice[0] > xmax:
                xmax = vertice[0]
        
        ymin = ymax = list(self.vertex())[0][1]
        for vertice in self.vertex():
            if vertice[1] < ymin:
                ymin = vertice[1]
            elif vertice[1] > ymax:
                ymax = vertice[1]
        
        return [[xmin, xmax], [ymin, ymax]]
    
    def vertices(self):
        r"""
        Return all vertex of ``self``, which is the point where three or 
        more line segments intersect

        OUTPUT: A set of `(x,y)` points

        """
        vertex = set()
        for i, component in enumerate(self._hypersurface):
            parametric_function = component[0]
            var = component[1][0].variables()[0]
            interval = self._parameter_intervals()[i]
            lower = interval[0].lower()
            upper = interval[0].upper()
            if lower != -infinity:
                x = parametric_function[0].subs(var==lower)
                y = parametric_function[1].subs(var==lower)
                vertex.add((x,y))
            if upper != infinity:
                x = parametric_function[0].subs(var==upper)
                y = parametric_function[1].subs(var==upper)
                vertex.add((x,y))
        return vertex
    
    def _parameter_intervals(self):
        r"""
        Return the intervals of each parameter of ``self``

        OUTPUT: A list of ``RealSet``

        """
        from sage.sets.real_set import RealSet

        intervals = []
        for component in self._hypersurface:
            if len(component[1]) == 1:
                interval = RealSet(component[1][0])
            else:
                lower = QQ(component[1][0].left())
                upper = QQ(component[1][1].right())
                interval = RealSet([lower,upper])
            intervals.append(interval)
        
        return intervals
                

    def plot(self):
        """
        Return the plot of ``self``

        OUTPUT: A Graphics object

        """
        from sage.plot.text import text

        combined_plot = Graphics()
        large_int = 1000
        intervals = self._parameter_intervals()
        for i, component in enumerate(self._hypersurface):
            var = component[1][0].variables()[0]
            parametric_function = component[0]
            order = component[2]
            interval = intervals[i]
            if interval[0].lower() == -infinity:
                lower = interval[0].upper() - large_int
                upper = interval[0].upper()
                midpoint = upper - 0.5
            elif interval[0].upper() == infinity:
                lower = interval[0].lower()
                upper = interval[0].lower() + large_int
                midpoint = lower + 0.5
            else:
                lower = interval[0].lower()
                upper = interval[0].upper()
                midpoint = (lower+upper)/2
                
            plot = parametric_plot(parametric_function, (var, lower, upper),
                                    color='red')

            if component[2] > 1: # add order if >= 2
                point = []
                for eq in component[0]:
                    value = eq.subs(var==midpoint)
                    point.append(value)
                text_order = text(str(order), (point[0], point[1]), 
                                  fontsize=16, color='black')
                combined_plot += plot + text_order
            else:
                combined_plot += plot

        # set default axes
        axes = self._axes()
        xmin, xmax = axes[0][0], axes[0][1]
        ymin, ymax = axes[1][0], axes[1][1]
        combined_plot.set_axes_range(xmin=xmin-1, xmax=xmax+1, 
                                     ymin=ymin-1, ymax=ymax+1)
        return combined_plot

    def _repr_(self):
        components = "\n".join([f"{row}" for row in self._hypersurface])
        return (f"Tropical curve of {self.poly} are \n[{components}]")
