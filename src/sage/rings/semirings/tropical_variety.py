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
    Tropical Hypersurface in 2 dimensions: 
    [[(0, r1), [r1 <= 0], 1]
    [(r2, 0), [r2 <= 0], 1]
    [(r3, r3), [r3 >= 0], 1]]
    sage: th1.dimension()
    2
    sage: th1.components()
    3
    sage: plot(th1)

::

    sage: S.<x,y> = QQ[]
    sage: c2 = -1*x^2
    sage: dict2 = {(0,0):0, (1,0):0, (0,2):0}
    sage: p2 = R(c2) + R(dict2)
    sage: th2 = p2.tropical_hypersurface(); th2
    Tropical Hypersurface in 2 dimensions:
    [[(0, r4), [r4 <= 0], 1]
    [(r5, 0), [r5 <= 0], 2]
    [(2*r7, r7), [0 <= r7, r7 <= (1/2)], 1]
    [(1, r8), [r8 <= (1/2)], 1]
    [(r9 + 1/2, r9), [(1/2) <= r9], 2]]
    sage: th2.plot()

::

    sage: c3 = 7+4*x+4*x*y+3*y^2+(-3)*x^2
    sage: dict3 = {(0,1):0}
    sage: p3 = R(c3) + R(dict3)
    sage: th3 = p3.tropical_hypersurface()
    sage: th3.plot()

TESTS:

    sage: T = TropicalSemiring(QQ, use_min=False)
    sage: R = PolynomialRing(T, 'a,x,y,z')
    sage: S.<a,x,y,z> = QQ[]
    sage: p1 = R(x*y + (-1/2)*x*z + 4*z^2 + a*x)
    sage: p1.tropical_hypersurface()
    Tropical Hypersurface in 4 dimensions: 
    [[(r17, r18, r17, r16), [r16 < min(r17 + 3/2, 1/2*r17 + 1/2*r18 - 3/2)], 1]
    [(r19 - 3/2, r21, r20, r19), [r20 + 3/2 < r19, r19 < r21 - 9/2, r20 < r21 - 6], 1]
    [(2*r22 - r24 + 3, r24, r23, r22), [max(1/2*r23 + 1/2*r24 - 3/2, r24 - 9/2) < r22], 1]
    [(r26, r27, r25 - 3/2, r25), [r26 + 3/2 < r25, r25 < r27 - 9/2, r26 < r27 - 6], 1]
    [(r30, 2*r28 - r29 + 3, r29, r28), [r28 < r29 + 3/2, r30 < r29], 1]
    [(r33, r31 + 9/2, r32, r31), [max(r32 + 3/2, r33 + 3/2) < r31], 1]]
    sage: p1.tropical_hypersurface().plot()
    Traceback (most recent call last):
    ...
    NotImplementedError: Can't visualize tropical hypersuface in 4 dimensions
    

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
    
    def __init__(self, lhs, rhs):
        """
        INPUT:

        - lhs -- a list of tropical polynomials
        - rhs -- a list of tropical polynomials 
        """
        from sage.rings.semirings.tropical_polynomial import TropicalPolynomial
        from sage.rings.semirings.tropical_mpolynomial import TropicalMPolynomial
        
        SageObject.__init__(self)
        if len(lhs) != len(rhs):
            raise ValueError("The length of inputs has to be equal")
        
        for i in range(len(lhs)):
            if not isinstance(lhs[i], (TropicalPolynomial, TropicalMPolynomial)):
                raise ValueError("Each input has to be a tropical polynomial")
            if not isinstance(rhs[i], (TropicalPolynomial, TropicalMPolynomial)):
                raise ValueError("Each input has to be a tropical polynomial") 
        self.lhs = lhs
        self.rhs = rhs

    def _repr_(self):
        equations = []
        for i in range(len(self.lhs)):
            equation = " = ".join([f"{self.lhs[i]}", f"{self.rhs[i]}"])
            equations.append(equation)
        all_eq = "\n".join([f"{row}" for row in equations])
        return (f"Tropical variety defined by \n{all_eq} ")

class TropicalCurve(TropicalVariety):
    r""""
    Represents a tropical hypersurface in `n=2` dimensions. The representation
    is in the form of list of lists, where the inner list represent each
    line segments of tropical roots
    """

    def __init__(self, *args):
        SageObject.__init__(self)
        self._hypersurface = []
        dim_param = len(args[0][0]) - 1
        vars = [SR.var('t{}'.format(i)) for i in range(1, dim_param+1)]
        for arg in args:
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
    
    def vertex(self):
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

        if self.dimension() > 3:
            raise NotImplementedError(f"Can't visualize tropical hypersuface"
                                      f" in {self.dimension()} dimensions")
        
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
        return (f"Tropical Hypersurface in {self.dimension()} dimensions: \n"
                f"[{components}]")
