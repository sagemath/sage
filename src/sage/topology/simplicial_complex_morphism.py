r"""
Morphisms of simplicial complexes

AUTHORS:

- Benjamin Antieau <d.ben.antieau@gmail.com> (2009.06)

- Travis Scrimshaw (2012-08-18): Made all simplicial complexes immutable to
  work with the homset cache.

This module implements morphisms of simplicial complexes. The input is given
by a dictionary on the vertex set of a simplicial complex. The initialization
checks that faces are sent to faces.

There is also the capability to create the fiber product of two morphisms with
the same codomain.

EXAMPLES::

    sage: S = SimplicialComplex([[0,2],[1,5],[3,4]], is_mutable=False)
    sage: H = Hom(S,S.product(S, is_mutable=False))
    sage: H.diagonal_morphism()
    Simplicial complex morphism:
      From: Simplicial complex with vertex set (0, 1, 2, 3, 4, 5) and facets {(0, 2), (1, 5), (3, 4)}
      To: Simplicial complex with 36 vertices and 18 facets
      Defn: [0, 1, 2, 3, 4, 5] --> ['L0R0', 'L1R1', 'L2R2', 'L3R3', 'L4R4', 'L5R5']

    sage: S = SimplicialComplex([[0,2],[1,5],[3,4]], is_mutable=False)
    sage: T = SimplicialComplex([[0,2],[1,3]], is_mutable=False)
    sage: f = {0:0,1:1,2:2,3:1,4:3,5:3}
    sage: H = Hom(S,T)
    sage: x = H(f)
    sage: x.image()
    Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2), (1, 3)}
    sage: x.is_surjective()
    True
    sage: x.is_injective()
    False
    sage: x.is_identity()
    False

    sage: S = simplicial_complexes.Sphere(2)
    sage: H = Hom(S,S)
    sage: i = H.identity()
    sage: i.image()
    Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)}
    sage: i.is_surjective()
    True
    sage: i.is_injective()
    True
    sage: i.is_identity()
    True

    sage: S = simplicial_complexes.Sphere(2)
    sage: H = Hom(S,S)
    sage: i = H.identity()
    sage: j = i.fiber_product(i)
    sage: j
    Simplicial complex morphism:
      From: Simplicial complex with 4 vertices and 4 facets
      To:   Minimal triangulation of the 2-sphere
      Defn: L0R0 |--> 0
            L1R1 |--> 1
            L2R2 |--> 2
            L3R3 |--> 3
    sage: S = simplicial_complexes.Sphere(2)
    sage: T = S.product(SimplicialComplex([[0,1]]), rename_vertices = False, is_mutable=False)
    sage: H = Hom(T,S)
    sage: T
    Simplicial complex with 8 vertices and 12 facets
    sage: sorted(T.vertices())
    [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (3, 1)]
    sage: f = {(0, 0): 0, (0, 1): 0, (1, 0): 1, (1, 1): 1, (2, 0): 2, (2, 1): 2, (3, 0): 3, (3, 1): 3}
    sage: x = H(f)
    sage: U = simplicial_complexes.Sphere(1)
    sage: G = Hom(U,S)
    sage: U
    Minimal triangulation of the 1-sphere
    sage: g = {0:0,1:1,2:2}
    sage: y = G(g)
    sage: z = y.fiber_product(x)
    sage: z                                     # this is the mapping path space
    Simplicial complex morphism:
      From: Simplicial complex with 6 vertices and ... facets
      To:   Minimal triangulation of the 2-sphere
      Defn: ['L0R(0, 0)', 'L0R(0, 1)', 'L1R(1, 0)', 'L1R(1, 1)', 'L2R(2, 0)', 'L2R(2, 1)'] --> [0, 0, 1, 1, 2, 2]
"""

#*****************************************************************************
# Copyright (C) 2009 D. Benjamin Antieau <d.ben.antieau@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#
#*****************************************************************************

from .simplicial_complex import Simplex, SimplicialComplex
from sage.matrix.constructor import matrix, zero_matrix
from sage.rings.integer_ring import ZZ
from sage.homology.chain_complex_morphism import ChainComplexMorphism
from sage.combinat.permutation import Permutation
from sage.algebras.steenrod.steenrod_algebra_misc import convert_perm
from sage.categories.morphism import Morphism
from sage.categories.homset import Hom
from sage.categories.simplicial_complexes import SimplicialComplexes


def is_SimplicialComplexMorphism(x):
    """
    Return ``True`` if and only if ``x`` is a morphism of simplicial complexes.

    EXAMPLES::

        sage: from sage.topology.simplicial_complex_morphism import is_SimplicialComplexMorphism
        sage: S = SimplicialComplex([[0,1],[3,4]], is_mutable=False)
        sage: H = Hom(S,S)
        sage: f = {0:0,1:1,3:3,4:4}
        sage: x = H(f)
        sage: is_SimplicialComplexMorphism(x)
        True

    """
    return isinstance(x, SimplicialComplexMorphism)


class SimplicialComplexMorphism(Morphism):
    """
    An element of this class is a morphism of simplicial complexes.
    """
    def __init__(self,f,X,Y):
        """
        Input is a dictionary ``f``, the domain ``X``, and the codomain ``Y``.

        One can define the dictionary on the vertices of `X`.

        EXAMPLES::

            sage: S = SimplicialComplex([[0,1],[2],[3,4],[5]], is_mutable=False)
            sage: H = Hom(S,S)
            sage: f = {0:0,1:1,2:2,3:3,4:4,5:5}
            sage: g = {0:0,1:1,2:0,3:3,4:4,5:0}
            sage: x = H(f)
            sage: y = H(g)
            sage: x == y
            False
            sage: x.image()
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5) and facets {(2,), (5,), (0, 1), (3, 4)}
            sage: y.image()
            Simplicial complex with vertex set (0, 1, 3, 4) and facets {(0, 1), (3, 4)}
            sage: x.image() == y.image()
            False
        """
        if not isinstance(X,SimplicialComplex) or not isinstance(Y,SimplicialComplex):
            raise ValueError("X and Y must be SimplicialComplexes")
        if not set(f.keys()) == set(X.vertices()):
            raise ValueError("f must be a dictionary from the vertex set of X to single values in the vertex set of Y")
        dim = X.dimension()
        Y_faces = Y.faces()
        for k in range(dim+1):
            for i in X.faces()[k]:
                tup = i.tuple()
                fi = []
                for j in tup:
                    fi.append(f[j])
                v = Simplex(set(fi))
            if v not in Y_faces[v.dimension()]:
                raise ValueError("f must be a dictionary from the vertices of X to the vertices of Y")
        self._vertex_dictionary = f
        Morphism.__init__(self, Hom(X,Y,SimplicialComplexes()))

    def __eq__(self,x):
        """
        Return ``True`` if and only if ``self == x``.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(2)
            sage: H = Hom(S,S)
            sage: i = H.identity()
            sage: i
            Simplicial complex endomorphism of Minimal triangulation of the 2-sphere
              Defn: 0 |--> 0
                    1 |--> 1
                    2 |--> 2
                    3 |--> 3
            sage: f = {0:0,1:1,2:2,3:2}
            sage: j = H(f)
            sage: i==j
            False

            sage: T = SimplicialComplex([[1,2]], is_mutable=False)
            sage: T
            Simplicial complex with vertex set (1, 2) and facets {(1, 2)}
            sage: G = Hom(T,T)
            sage: k = G.identity()
            sage: g = {1:1,2:2}
            sage: l = G(g)
            sage: k == l
            True
        """
        if not isinstance(x,SimplicialComplexMorphism) or self.codomain() != x.codomain() or self.domain() != x.domain() or self._vertex_dictionary != x._vertex_dictionary:
            return False
        else:
            return True

    def __call__(self,x,orientation=False):
        """
        Input is a simplex of the domain. Output is the image simplex.

        If the optional argument ``orientation`` is ``True``, then this
        returns a pair ``(image simplex, oriented)`` where ``oriented``
        is 1 or `-1` depending on whether the map preserves or reverses
        the orientation of the image simplex.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(2)
            sage: T = simplicial_complexes.Sphere(3)
            sage: S
            Minimal triangulation of the 2-sphere
            sage: T
            Minimal triangulation of the 3-sphere
            sage: f = {0:0,1:1,2:2,3:3}
            sage: H = Hom(S,T)
            sage: x = H(f)
            sage: from sage.topology.simplicial_complex import Simplex
            sage: x(Simplex([0,2,3]))
            (0, 2, 3)

        An orientation-reversing example::

            sage: X = SimplicialComplex([[0,1]], is_mutable=False)
            sage: g = Hom(X,X)({0:1, 1:0})
            sage: g(Simplex([0,1]))
            (0, 1)
            sage: g(Simplex([0,1]), orientation=True)
            ((0, 1), -1)
        """
        dim = self.domain().dimension()
        if not isinstance(x, Simplex) or x.dimension() > dim or x not in self.domain().faces()[x.dimension()]:
            raise ValueError("x must be a simplex of the source of f")
        tup = x.tuple()
        fx = []
        for j in tup:
            fx.append(self._vertex_dictionary[j])
        if orientation:
            if len(set(fx)) == len(tup):
                oriented = Permutation(convert_perm(fx)).signature()
            else:
                oriented = 1
            return (Simplex(set(fx)), oriented)
        else:
            return Simplex(set(fx))

    def _repr_type(self):
        """
        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(1)
            sage: T = simplicial_complexes.Sphere(2)
            sage: H = Hom(S,T)
            sage: f = {0:0,1:1,2:2}
            sage: H(f)._repr_type()
            'Simplicial complex'
        """
        return "Simplicial complex"

    def _repr_defn(self):
        """
        If there are fewer than 5 vertices, print the image of each vertex
        on a separate line. Otherwise, print the map as a single line.

        EXAMPLES::

            sage: S = simplicial_complexes.Simplex(1)
            sage: print(Hom(S,S).identity()._repr_defn())
            0 |--> 0
            1 |--> 1
            sage: T = simplicial_complexes.Torus()
            sage: print(Hom(T,T).identity()._repr_defn())
            [0, 1, 2, 3, 4, 5, 6] --> [0, 1, 2, 3, 4, 5, 6]
        """
        vd = self._vertex_dictionary
        try:
            keys = sorted(vd.keys())
        except TypeError:
            keys = sorted(vd.keys(), key=str)
        if len(vd) < 5:
            return '\n'.join("{} |--> {}".format(v, vd[v]) for v in keys)
        domain = list(vd.keys())
        try:
            domain = sorted(domain)
        except TypeError:
            domain = sorted(domain, key=str)
        codomain = [vd[v] for v in domain]
        return "{} --> {}".format(domain, codomain)

    def associated_chain_complex_morphism(self,base_ring=ZZ,augmented=False,cochain=False):
        """
        Return the associated chain complex morphism of ``self``.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(1)
            sage: T = simplicial_complexes.Sphere(2)
            sage: H = Hom(S,T)
            sage: f = {0:0,1:1,2:2}
            sage: x = H(f)
            sage: x
            Simplicial complex morphism:
              From: Minimal triangulation of the 1-sphere
              To:   Minimal triangulation of the 2-sphere
              Defn: 0 |--> 0
                    1 |--> 1
                    2 |--> 2
            sage: a = x.associated_chain_complex_morphism()
            sage: a
            Chain complex morphism:
              From: Chain complex with at most 2 nonzero terms over Integer Ring
              To:   Chain complex with at most 3 nonzero terms over Integer Ring
            sage: a._matrix_dictionary
            {0: [1 0 0]
             [0 1 0]
             [0 0 1]
             [0 0 0], 1: [1 0 0]
             [0 1 0]
             [0 0 0]
             [0 0 1]
             [0 0 0]
             [0 0 0], 2: []}
            sage: x.associated_chain_complex_morphism(augmented=True)
            Chain complex morphism:
              From: Chain complex with at most 3 nonzero terms over Integer Ring
              To:   Chain complex with at most 4 nonzero terms over Integer Ring
            sage: x.associated_chain_complex_morphism(cochain=True)
            Chain complex morphism:
              From: Chain complex with at most 3 nonzero terms over Integer Ring
              To:   Chain complex with at most 2 nonzero terms over Integer Ring
            sage: x.associated_chain_complex_morphism(augmented=True,cochain=True)
            Chain complex morphism:
              From: Chain complex with at most 4 nonzero terms over Integer Ring
              To:   Chain complex with at most 3 nonzero terms over Integer Ring
            sage: x.associated_chain_complex_morphism(base_ring=GF(11))
            Chain complex morphism:
              From: Chain complex with at most 2 nonzero terms over Finite Field of size 11
              To:   Chain complex with at most 3 nonzero terms over Finite Field of size 11

        Some simplicial maps which reverse the orientation of a few simplices::

            sage: g = {0:1, 1:2, 2:0}
            sage: H(g).associated_chain_complex_morphism()._matrix_dictionary
            {0: [0 0 1]
             [1 0 0]
             [0 1 0]
             [0 0 0], 1: [ 0 -1  0]
             [ 0  0 -1]
             [ 0  0  0]
             [ 1  0  0]
             [ 0  0  0]
             [ 0  0  0], 2: []}
            sage: X = SimplicialComplex([[0, 1]], is_mutable=False)
            sage: Hom(X,X)({0:1, 1:0}).associated_chain_complex_morphism()._matrix_dictionary
            {0: [0 1]
             [1 0], 1: [-1]}
        """
        max_dim = max(self.domain().dimension(),self.codomain().dimension())
        min_dim = min(self.domain().dimension(),self.codomain().dimension())
        matrices = {}
        if augmented is True:
            m = matrix(base_ring,1,1,1)
            if not cochain:
                matrices[-1] = m
            else:
                matrices[-1] = m.transpose()
        for dim in range(min_dim+1):
            X_faces = self.domain()._n_cells_sorted(dim)
            Y_faces = self.codomain()._n_cells_sorted(dim)
            num_faces_X = len(X_faces)
            num_faces_Y = len(Y_faces)
            mval = [0 for i in range(num_faces_X*num_faces_Y)]
            for i in X_faces:
                y, oriented = self(i, orientation=True)
                if y.dimension() < dim:
                    pass
                else:
                    mval[X_faces.index(i)+(Y_faces.index(y)*num_faces_X)] = oriented
            m = matrix(base_ring,num_faces_Y,num_faces_X,mval,sparse=True)
            if not cochain:
                matrices[dim] = m
            else:
                matrices[dim] = m.transpose()
        for dim in range(min_dim+1,max_dim+1):
            try:
                l1 = len(self.codomain().n_cells(dim))
            except KeyError:
                l1 = 0
            try:
                l2 = len(self.domain().n_cells(dim))
            except KeyError:
                l2 = 0
            m = zero_matrix(base_ring,l1,l2,sparse=True)
            if not cochain:
                matrices[dim] = m
            else:
                matrices[dim] = m.transpose()
        if not cochain:
            return ChainComplexMorphism(matrices,
                    self.domain().chain_complex(base_ring=base_ring,augmented=augmented,cochain=cochain),
                    self.codomain().chain_complex(base_ring=base_ring,augmented=augmented,cochain=cochain))
        return ChainComplexMorphism(matrices,
                self.codomain().chain_complex(base_ring=base_ring,augmented=augmented,cochain=cochain),
                self.domain().chain_complex(base_ring=base_ring,augmented=augmented,cochain=cochain))

    def image(self):
        """
        Computes the image simplicial complex of `f`.

        EXAMPLES::

            sage: S = SimplicialComplex([[0,1],[2,3]], is_mutable=False)
            sage: T = SimplicialComplex([[0,1]], is_mutable=False)
            sage: f = {0:0,1:1,2:0,3:1}
            sage: H = Hom(S,T)
            sage: x = H(f)
            sage: x.image()
            Simplicial complex with vertex set (0, 1) and facets {(0, 1)}

            sage: S = SimplicialComplex(is_mutable=False)
            sage: H = Hom(S,S)
            sage: i = H.identity()
            sage: i.image()
            Simplicial complex with vertex set () and facets {()}
            sage: i.is_surjective()
            True
            sage: S = SimplicialComplex([[0,1]], is_mutable=False)
            sage: T = SimplicialComplex([[0,1], [0,2]], is_mutable=False)
            sage: f = {0:0,1:1}
            sage: g = {0:0,1:1}
            sage: k = {0:0,1:2}
            sage: H = Hom(S,T)
            sage: x = H(f)
            sage: y = H(g)
            sage: z = H(k)
            sage: x == y
            True
            sage: x == z
            False
            sage: x.image()
            Simplicial complex with vertex set (0, 1) and facets {(0, 1)}
            sage: y.image()
            Simplicial complex with vertex set (0, 1) and facets {(0, 1)}
            sage: z.image()
            Simplicial complex with vertex set (0, 2) and facets {(0, 2)}

        """
        fa = [self(i) for i in self.domain().facets()]
        return SimplicialComplex(fa, maximality_check=True)

    def is_surjective(self):
        """
        Return ``True`` if and only if ``self`` is surjective.

        EXAMPLES::

            sage: S = SimplicialComplex([(0,1,2)], is_mutable=False)
            sage: S
            Simplicial complex with vertex set (0, 1, 2) and facets {(0, 1, 2)}
            sage: T = SimplicialComplex([(0,1)], is_mutable=False)
            sage: T
            Simplicial complex with vertex set (0, 1) and facets {(0, 1)}
            sage: H = Hom(S,T)
            sage: x = H({0:0,1:1,2:1})
            sage: x.is_surjective()
            True

            sage: S = SimplicialComplex([[0,1],[2,3]], is_mutable=False)
            sage: T = SimplicialComplex([[0,1]], is_mutable=False)
            sage: f = {0:0,1:1,2:0,3:1}
            sage: H = Hom(S,T)
            sage: x = H(f)
            sage: x.is_surjective()
            True
        """
        return self.codomain() == self.image()

    def is_injective(self):
        """
        Return ``True`` if and only if ``self`` is injective.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(1)
            sage: T = simplicial_complexes.Sphere(2)
            sage: U = simplicial_complexes.Sphere(3)
            sage: H = Hom(T,S)
            sage: G = Hom(T,U)
            sage: f = {0:0,1:1,2:0,3:1}
            sage: x = H(f)
            sage: g = {0:0,1:1,2:2,3:3}
            sage: y = G(g)
            sage: x.is_injective()
            False
            sage: y.is_injective()
            True

        """
        v = [self._vertex_dictionary[i[0]] for i in self.domain().faces()[0]]
        for i in v:
            if v.count(i) > 1:
                return False
        return True

    def is_identity(self):
        """
        If ``self`` is an identity morphism, returns ``True``.
        Otherwise, ``False``.

        EXAMPLES::

            sage: T = simplicial_complexes.Sphere(1)
            sage: G = Hom(T,T)
            sage: T
            Minimal triangulation of the 1-sphere
            sage: j = G({0:0,1:1,2:2})
            sage: j.is_identity()
            True

            sage: S = simplicial_complexes.Sphere(2)
            sage: T = simplicial_complexes.Sphere(3)
            sage: H = Hom(S,T)
            sage: f = {0:0,1:1,2:2,3:3}
            sage: x = H(f)
            sage: x
            Simplicial complex morphism:
              From: Minimal triangulation of the 2-sphere
              To:   Minimal triangulation of the 3-sphere
              Defn: 0 |--> 0
                    1 |--> 1
                    2 |--> 2
                    3 |--> 3
            sage: x.is_identity()
            False
        """
        if self.domain() != self.codomain():
            return False
        else:
            f = dict()
            for i in self.domain().vertices():
                f[i] = i
            if self._vertex_dictionary != f:
                return False
            else:
                return True

    def fiber_product(self, other, rename_vertices = True):
        """
        Fiber product of ``self`` and ``other``. Both morphisms should have
        the same codomain. The method returns a morphism of simplicial
        complexes, which is the morphism from the space of the fiber product
        to the codomain.

        EXAMPLES::

            sage: S = SimplicialComplex([[0,1],[1,2]], is_mutable=False)
            sage: T = SimplicialComplex([[0,2],[1]], is_mutable=False)
            sage: U = SimplicialComplex([[0,1],[2]], is_mutable=False)
            sage: H = Hom(S,U)
            sage: G = Hom(T,U)
            sage: f = {0:0,1:1,2:0}
            sage: g = {0:0,1:1,2:1}
            sage: x = H(f)
            sage: y = G(g)
            sage: z = x.fiber_product(y)
            sage: z
            Simplicial complex morphism:
              From: Simplicial complex with 4 vertices and facets {...}
              To:   Simplicial complex with vertex set (0, 1, 2) and facets {(2,), (0, 1)}
              Defn: L0R0 |--> 0
                    L1R1 |--> 1
                    L1R2 |--> 1
                    L2R0 |--> 0
        """
        if self.codomain() != other.codomain():
            raise ValueError("self and other must have the same codomain.")
        X = self.domain().product(other.domain(),rename_vertices = rename_vertices)
        v = []
        f = dict()
        eff1 = self.domain().vertices()
        eff2 = other.domain().vertices()
        for i in eff1:
            for j in eff2:
                if self(Simplex([i])) == other(Simplex([j])):
                    if rename_vertices:
                        v.append("L"+str(i)+"R"+str(j))
                        f["L"+str(i)+"R"+str(j)] = self._vertex_dictionary[i]
                    else:
                        v.append((i,j))
                        f[(i,j)] = self._vertex_dictionary[i]
        return SimplicialComplexMorphism(f, X.generated_subcomplex(v), self.codomain())

    def mapping_torus(self):
        r"""
        The mapping torus of a simplicial complex endomorphism

        The mapping torus is the simplicial complex formed by taking
        the product of the domain of ``self`` with a `4` point
        interval `[I_0, I_1, I_2, I_3]` and identifying vertices of
        the form `(I_0, v)` with `(I_3, w)` where `w` is the image of
        `v` under the given morphism.

        See :wikipedia:`Mapping torus`

        EXAMPLES::

            sage: C = simplicial_complexes.Sphere(1)            # Circle
            sage: T = Hom(C,C).identity().mapping_torus() ; T   # Torus
            Simplicial complex with 9 vertices and 18 facets
            sage: T.homology() == simplicial_complexes.Torus().homology()
            True

            sage: f = Hom(C,C)({0:0,1:2,2:1})
            sage: K = f.mapping_torus() ; K  # Klein Bottle
            Simplicial complex with 9 vertices and 18 facets
            sage: K.homology() == simplicial_complexes.KleinBottle().homology()
            True

        TESTS::

            sage: g = Hom(simplicial_complexes.Simplex([1]),C)({1:0})
            sage: g.mapping_torus()
            Traceback (most recent call last):
            ...
            ValueError: self must have the same domain and codomain.
        """
        if self.domain() != self.codomain():
            raise ValueError("self must have the same domain and codomain.")
        map_dict = self._vertex_dictionary
        interval = SimplicialComplex([["I0","I1"],["I1","I2"]])
        product = interval.product(self.domain(),False)
        facets = list(product.maximal_faces())
        for facet in self.domain()._facets:
            left = [ ("I0",v) for v in facet ]
            right = [ ("I2",map_dict[v]) for v in facet ]
            for i in range(facet.dimension()+1):
                facets.append(tuple(left[:i+1]+right[i:]))
        return SimplicialComplex(facets)

    def induced_homology_morphism(self, base_ring=None, cohomology=False):
        """
        The map in (co)homology induced by this map

        INPUT:

        - ``base_ring`` -- must be a field (optional, default ``QQ``)

        - ``cohomology`` -- boolean (optional, default ``False``). If
          ``True``, the map induced in cohomology rather than homology.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(1)
            sage: T = S.product(S, is_mutable=False)
            sage: H = Hom(S,T)
            sage: diag = H.diagonal_morphism()
            sage: h = diag.induced_homology_morphism(QQ)
            sage: h
            Graded vector space morphism:
              From: Homology module of Minimal triangulation of the 1-sphere over Rational Field
              To:   Homology module of Simplicial complex with 9 vertices and 18 facets over Rational Field
              Defn: induced by:
                Simplicial complex morphism:
                  From: Minimal triangulation of the 1-sphere
                  To:   Simplicial complex with 9 vertices and 18 facets
                  Defn: 0 |--> L0R0
                        1 |--> L1R1
                        2 |--> L2R2

        We can view the matrix form for the homomorphism::

            sage: h.to_matrix(0) # in degree 0
            [1]
            sage: h.to_matrix(1) # in degree 1
            [1]
            [1]
            sage: h.to_matrix()  # the entire homomorphism
            [1|0]
            [-+-]
            [0|1]
            [0|1]
            [-+-]
            [0|0]

        The map on cohomology should be dual to the map on homology::

            sage: coh = diag.induced_homology_morphism(QQ, cohomology=True)
            sage: coh.to_matrix(1)
            [1 1]
            sage: h.to_matrix() == coh.to_matrix().transpose()
            True

        We can evaluate the map on (co)homology classes::

            sage: x,y = list(T.cohomology_ring(QQ).basis(1))
            sage: coh(x)
            h^{1,0}
            sage: coh(2*x+3*y)
            5*h^{1,0}

        Note that the complexes must be immutable for this to
        work. Many, but not all, complexes are immutable when
        constructed::

            sage: S.is_immutable()
            True
            sage: S.barycentric_subdivision().is_immutable()
            False
            sage: S2 = S.suspension()
            sage: S2.is_immutable()
            False
            sage: h = Hom(S,S2)({0: 0, 1:1, 2:2}).induced_homology_morphism()
            Traceback (most recent call last):
            ...
            ValueError: the domain and codomain complexes must be immutable
            sage: S2.set_immutable(); S2.is_immutable()
            True
            sage: h = Hom(S,S2)({0: 0, 1:1, 2:2}).induced_homology_morphism()
        """
        from sage.homology.homology_morphism import InducedHomologyMorphism
        return InducedHomologyMorphism(self, base_ring, cohomology)

    def is_contiguous_to(self, other):
        r"""
        Return ``True`` if ``self`` is contiguous to ``other``.

        Two morphisms `f_0, f_1: K \to L` are *contiguous* if for any
        simplex `\sigma \in K`, the union `f_0(\sigma) \cup
        f_1(\sigma)` is a simplex in `L`. This is not a transitive
        relation, but it induces an equivalence relation on simplicial
        maps: `f` is equivalent to `g` if there is a finite sequence
        `f_0 = f`, `f_1`, ..., `f_n = g` such that `f_i` and `f_{i+1}`
        are contiguous for each `i`.

        This is related to maps being homotopic: if they are
        contiguous, then they induce homotopic maps on the geometric
        realizations. Given two homotopic maps on the geometric
        realizations, then after barycentrically subdividing `n` times
        for some `n`, the maps have simplicial approximations which
        are in the same contiguity class. (This last fact is only true
        if the domain is a *finite* simplicial complex, by the way.)

        See Section 3.5 of Spanier [Spa1966]_ for details.

        ALGORITHM:

        It is enough to check when `\sigma` ranges over the facets.

        INPUT:

        - ``other`` -- a simplicial complex morphism with the same
          domain and codomain as ``self``

        EXAMPLES::

            sage: K = simplicial_complexes.Simplex(1)
            sage: L = simplicial_complexes.Sphere(1)
            sage: H = Hom(K, L)
            sage: f = H({0: 0, 1: 1})
            sage: g = H({0: 0, 1: 0})
            sage: f.is_contiguous_to(f)
            True
            sage: f.is_contiguous_to(g)
            True
            sage: h = H({0: 1, 1: 2})
            sage: f.is_contiguous_to(h)
            False

        TESTS::

            sage: one = Hom(K,K).identity()
            sage: one.is_contiguous_to(f)
            False
            sage: one.is_contiguous_to(3) # nonsensical input
            False
        """
        if not isinstance(other, SimplicialComplexMorphism):
            return False
        if self.codomain() != other.codomain() or self.domain() != other.domain():
            return False
        domain = self.domain()
        codomain = self.codomain()
        return all(Simplex(self(sigma).set().union(other(sigma))) in codomain
                   for sigma in domain.facets())

