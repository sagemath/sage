from sage.homology.homology_group import HomologyGroup
from sage.structure.sage_object import SageObject
from .simplicial_complex import SimplicialComplex
from .simplicial_complex_examples import Sphere, Simplex
from itertools import combinations

#TODO's:
# - Documentation
# - __eq__, __ne__
# - and a lot more ...
# - latex

class MomentAngleComplex(SimplicialComplex):
    def __init__(self, X, construct=False):
        # check whether X is an instance of SimplicialComplex
        self._simplicial_complex = X
        self._constructed = construct
        self._moment_angle_complex = None

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
                        Y.append(Simplex(2))
                        Ys.append("D_2")
                    else:
                        Y.append(Sphere(1))
                        Ys.append("S_1")

                self._components.append(Y)
                self._symbolic_components.append(Ys)

        #print('\n'.join(str(comps) for comps in self._components))
        #print('\n'.join(str(comps) for comps in self._symbolic_components))

        if construct:
            self.construct()

    def construct(self):
        #check whether it was constructed first
        self._moment_angle_complex = SimplicialComplex()

        for component in self._components:
            x = component[0]
            for j in range(1, len(component)-1):
                x = x.product(component[j])

            self._moment_angle_complex = self._moment_angle_complex.disjoint_union(x)
            #print(self._moment_angle_complex)

        self._constructed = True

    def vertices(self):
        return self.z_delta.vertices()

    def homology(self):
        return self._moment_angle_complex.homology()
