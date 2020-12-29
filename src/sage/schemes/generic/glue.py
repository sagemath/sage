"""
Scheme obtained by gluing two other schemes
"""

#*******************************************************************************
#  Copyright (C) 2006 William Stein
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*******************************************************************************

from . import morphism
from . import scheme

class GluedScheme(scheme.Scheme):
    r"""
    INPUT:


    -  ``f`` - open immersion from a scheme U to a scheme
       X

    -  ``g`` - open immersion from U to a scheme Y


    OUTPUT: The scheme obtained by gluing X and Y along the open set
    U.

    .. note::

       Checking that `f` and `g` are open
       immersions is not implemented.
    """
    def __init__(self, f, g, check=True):
        if check:
            if not morphism.is_SchemeMorphism(f):
                raise TypeError("f (=%s) must be a scheme morphism"%f)
            if not morphism.is_SchemeMorphism(g):
                raise TypeError("g (=%s) must be a scheme morphism"%g)
            if f.domain() != g.domain():
                raise ValueError("f (=%s) and g (=%s) must have the same domain"%(f,g))
        self.__f = f
        self.__g = g

    def gluing_maps(self):
        return self.__f, self.__g

    def _repr_(self):
        return "Scheme obtained by gluing X and Y along U, where\n  X: %s\n  Y: %s\n  U: %s"%(
            self.__f.codomain(), self.__g.codomain(), self.__f.domain())

