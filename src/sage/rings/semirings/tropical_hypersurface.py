r"""
Base class for tropical hypersurface objects

<Description>

AUTHORS:

- Verrel Rievaldo Wijaya

EXAMPLES:

    sage: T = TropicalSemiring(QQ, use_min=False)
    sage: R = PolynomialRing(T, 'x,y')
    sage: dict1 = {(0,0):0, (1,0):0, (0,1):0}
    sage: p1 = R(dict1)
    sage: th1 = p1.tropical_hypersurface(); th1
    Tropical Hypersurface in 2 dimensions: 
    [[(0, r1), [r1 < 0], 1]
    [(r2, 0), [r2 < 0], 1]
    [(r3, r3), [r3 > 0], 1]]
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
    [[(0, r4), [r4 < 0], 1]
    [(r5, 0), [r5 < 0], 2]
    [(2*r7, r7), [0 < r7, r7 < (1/2)], 1]
    [(1, r8), [r8 < (1/2)], 1]
    [(r9 + 1/2, r9), [(1/2) < r9], 2]]
    sage: th2.plot()

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

class TropicalHypersurface(SageObject):
    r""""
    Represents a tropical hypersurface in `n\geq 2` dimensions. The representation
    is in the form of list of lists, where the inner list represent each
    component of the surfaces

    """

    def __init__(self, *args):
        SageObject.__init__(self)
        self._hypersurface = []
        for arg in args:
            if not isinstance(arg, list):
                raise ValueError("Input must be a list")
            else:
                self._hypersurface.append(arg)

    def dimension(self):
        return len(self._hypersurface[0][0])
    
    def components(self):
        return len(self._hypersurface)

    def plot(self):
        """
        Return the plot of the tropical hypersurface ``self``


        """
        from sage.plot.text import text

        if self.dimension() > 3:
            raise NotImplementedError(f"Can't visualize tropical hypersuface"
                                      f" in {self.dimension()} dimensions")
        
        combined_plot = Graphics()
        for component in self._hypersurface:
            var = component[1][0].variables()[0]
            parametric_function = component[0]
            order = component[2]
            if len(component[1]) == 2:
                lower = QQ(component[1][0].left())
                upper = QQ(component[1][1].right())
            else:
                inequality = component[1][0]
                symbol = str(inequality).split()[1]
                if inequality.left().is_numeric():
                    number = QQ(inequality.left())
                    if symbol == '<':
                        pos_number = 'left'
                    else:
                        pos_number = 'right'
                else:
                    number = QQ(inequality.right())
                    if symbol == '<':
                        pos_number = 'right'
                    else:
                        pos_number = 'left'
                
                if pos_number == 'left':
                    lower = number
                    upper = number + 1
                else:
                    lower = number - 1
                    upper = number
            
            plot = parametric_plot(parametric_function, (var, lower, upper),
                                    color='red')

            if component[2] > 1:
                midpoint = (lower + upper)/2
                point = []
                for eq in component[0]:
                    value = eq.subs(var==midpoint)
                    point.append(value)
                text_order = text(str(order), (point[0], point[1]), 
                                  fontsize=16, color='black')
                combined_plot += plot + text_order
            else:
                combined_plot += plot
                    
        return combined_plot   

    def _repr_(self):
        components = "\n".join([f"{row}" for row in self._hypersurface])
        return (f"Tropical Hypersurface in {self.dimension()} dimensions: \n"
                f"[{components}]")
