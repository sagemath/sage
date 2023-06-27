from sage.homology.homology_group import HomologyGroup
from sage.structure.sage_object import SageObject
from .simplicial_complex import SimplicialComplex
from .simplicial_complex_examples import Sphere, Simplex
from itertools import combinations

#TODO's : 
# - Documentation
# - __eq__, __ne__
# - and a lot more ...

class MomentAngleComplex(SimplicialComplex):
    def __init__(self, X, construct=False):
        # check whether X is an instance of SimplicialComplex
        self._simplicial_complex = X
        self._constructed = construct

        n = X.dimension()
        vertices = X.vertices()
        k = len(vertices);
        self._components = []
        self._symbolic_components = []

        for face in X.faces().values():
            for subcomplex in face:
                Y = []
                Ys = []
                for j in vertices:
                    if j in subcomplex:
                        Y.append(SimplicialComplex(Simplex(2)))
                        Ys.append("D2")
                    else:
                        Y.append(SimplicialComplex(Sphere(1)))
                        Ys.append("S1")

                self._components.append(Y)
                self._symbolic_components.append(Ys)

        print(self._components)
        print(self._symbolic_components)

        if construct:
            self.construct()

        Z = SimplicialComplex()

    def construct(self):
        for component in self._components:
            x = component[0]
            for y in component:
                # this is terribly slow, look for solutions
                # x = x.product(y, rename_vertices=False)
                #print(x)

        self._constructed = True

    def vertices(self):
        return self.z_delta.vertices()

    def simplices(self):
        return self.z_delta.simplices()

    def homology(self):
        return self.z_delta_topological.homology()
