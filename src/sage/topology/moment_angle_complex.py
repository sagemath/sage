r"""
Moment angle complexes

<Paragraph description>

EXAMPLES::

<Lots and lots of examples>

AUTHORS:

- Ognjen Petrov (2023-06-25): initial version

"""

# ****************************************************************************
#       Copyright (C) 2013 Ognjen Petrov <ognjenpetrov@yahoo.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.homology.homology_group import HomologyGroup
from sage.structure.sage_object import SageObject
from .simplicial_complex import SimplicialComplex
from .simplicial_complex_examples import Sphere, Simplex
from itertools import combinations

#TODO's:
# - Documentation
# - and a lot more ...
# - latex

class MomentAngleComplex(SageObject): # should this inherit SimplicialComplex
    def __init__(self, simplicial_complex, construct=False):
        if not isinstance(simplicial_complex, SimplicialComplex):
            raise ValueError("simplicial_complex must be a simplicial complex")

        self._simplicial_complex = simplicial_complex
        self._moment_angle_complex = None

        vertices = simplicial_complex.vertices()
        k = len(vertices)
        self._components = []
        self._symbolic_components = []

        for face in simplicial_complex.faces().values():
            for subcomplex in face:
                Y = []
                Ys = []
                for j in vertices:
                    if j in subcomplex:
                        Y.append(Simplex(2))
                        Ys.append("D_2")
                    else:
                        Y.append(Sphere(1))
                        Ys.append("S_1")

                self._components.append(Y)
                self._symbolic_components.append(Ys)

        #print('\n'.join(str(comps) for comps in self._components))
        #print('\n'.join(str(comps) for comps in self._symbolic_components))

        self._constructed = False
        if construct:
            self.construct()

    def __eq__(self, other):
        #we consider them equal if they have the same corresponding simplicial complexes
        return isinstance(other, MomentAngleComplex) and self._simplicial_complex.__eq__(other._simplicial_complex)

    def __ne__(self, right):
        # if not equal
        return not self.__eq__(right)

    def construct(self):
        if self._constructed:
            return

        self._moment_angle_complex = SimplicialComplex()
        for component in self._components:
            x = component[0]
            for j in range(1, len(component)-1):
                x = x.product(component[j])

            self._moment_angle_complex = self._moment_angle_complex.disjoint_union(x)

        self._constructed = True

    def vertices(self):
        return self._moment_angle_complex.vertices()

    def homology(self):
        return self._moment_angle_complex.homology()
