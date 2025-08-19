# sage.doctest: needs sage.graphs
r"""
Homsets between simplicial complexes

AUTHORS:

- Travis Scrimshaw (2012-08-18): Made all simplicial complexes immutable to
  work with the homset cache.

EXAMPLES::

    sage: S = simplicial_complexes.Sphere(1)
    sage: T = simplicial_complexes.Sphere(2)
    sage: H = Hom(S,T)
    sage: f = {0:0,1:1,2:3}
    sage: x = H(f)
    sage: x
    Simplicial complex morphism:
      From: Minimal triangulation of the 1-sphere
      To: Minimal triangulation of the 2-sphere
      Defn: 0 |--> 0
            1 |--> 1
            2 |--> 3
    sage: x.is_injective()
    True
    sage: x.is_surjective()
    False
    sage: x.image()
    Simplicial complex with vertex set (0, 1, 3) and facets {(0, 1), (0, 3), (1, 3)}
    sage: from sage.topology.simplicial_complex import Simplex
    sage: s = Simplex([1,2])
    sage: x(s)
    (1, 3)

TESTS::

    sage: S = simplicial_complexes.Sphere(1)
    sage: T = simplicial_complexes.Sphere(2)
    sage: H = Hom(S,T)
    sage: loads(dumps(H)) == H
    True
"""

# ****************************************************************************
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
#                  https://www.gnu.org/licenses/
#
# ****************************************************************************

import sage.categories.homset
from .simplicial_complex_morphism import SimplicialComplexMorphism


def is_SimplicialComplexHomset(x) -> bool:
    """
    Return ``True`` if and only if ``x`` is a simplicial complex homspace.

    EXAMPLES::

        sage: S = SimplicialComplex(is_mutable=False)
        sage: T = SimplicialComplex(is_mutable=False)
        sage: H = Hom(S, T)
        sage: H
        Set of Morphisms from Simplicial complex with vertex set () and facets {()}
         to Simplicial complex with vertex set () and facets {()}
         in Category of finite simplicial complexes
        sage: from sage.topology.simplicial_complex_homset import is_SimplicialComplexHomset
        sage: is_SimplicialComplexHomset(H)
        doctest:warning...
        DeprecationWarning: the function is_SimplicialComplexHomset is deprecated;
        use 'isinstance(..., SimplicialComplexHomset)' instead
        See https://github.com/sagemath/sage/issues/37922 for details.
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(37922, "the function is_SimplicialComplexHomset is deprecated; use 'isinstance(..., SimplicialComplexHomset)' instead")
    return isinstance(x, SimplicialComplexHomset)


class SimplicialComplexHomset(sage.categories.homset.Homset):
    def __call__(self, f):
        """
        INPUT:

        - ``f`` -- dictionary with keys exactly the vertices of the domain
          and values vertices of the codomain

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(3)
            sage: T = simplicial_complexes.Sphere(2)
            sage: f = {0:0,1:1,2:2,3:2,4:2}
            sage: H = Hom(S,T)
            sage: x = H(f)
            sage: x
            Simplicial complex morphism:
              From: Minimal triangulation of the 3-sphere
              To: Minimal triangulation of the 2-sphere
              Defn: [0, 1, 2, 3, 4] --> [0, 1, 2, 2, 2]
        """
        return SimplicialComplexMorphism(f, self.domain(), self.codomain())

    def diagonal_morphism(self, rename_vertices=True):
        r"""
        Return the diagonal morphism in `Hom(S, S \times S)`.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(2)
            sage: H = Hom(S,S.product(S, is_mutable=False))
            sage: d = H.diagonal_morphism(); d
            Simplicial complex morphism:
              From: Minimal triangulation of the 2-sphere
              To:   Simplicial complex with 16 vertices and 96 facets
              Defn: 0 |--> L0R0
                    1 |--> L1R1
                    2 |--> L2R2
                    3 |--> L3R3

            sage: T = SimplicialComplex([[0], [1]], is_mutable=False)
            sage: U = T.product(T, rename_vertices=False, is_mutable=False)
            sage: G = Hom(T, U)
            sage: e = G.diagonal_morphism(rename_vertices=False); e
            Simplicial complex morphism:
              From: Simplicial complex with vertex set (0, 1) and facets {(0,), (1,)}
              To:   Simplicial complex with 4 vertices and
                    facets {((0, 0),), ((0, 1),), ((1, 0),), ((1, 1),)}
              Defn: 0 |--> (0, 0)
                    1 |--> (1, 1)
        """
        # Preserve whether the codomain is mutable when renaming the vertices.
        mutable = self._codomain.is_mutable()
        X = self._domain.product(self._domain, rename_vertices=rename_vertices,
                                 is_mutable=mutable)
        if self._codomain != X:
            raise TypeError("diagonal morphism is only defined for Hom(X,XxX)")
        f = {}
        if rename_vertices:
            f = {i: "L{0}R{0}".format(i) for i in self._domain.vertices()}
        else:
            f = {i: (i, i) for i in self._domain.vertices()}
        return SimplicialComplexMorphism(f, self._domain, X)

    def identity(self):
        """
        Return the identity morphism of `Hom(S,S)`.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(2)
            sage: H = Hom(S, S)
            sage: i = H.identity()
            sage: i.is_identity()
            True

            sage: T = SimplicialComplex([[0,1]], is_mutable=False)
            sage: G = Hom(T, T)
            sage: G.identity()
            Simplicial complex endomorphism of
             Simplicial complex with vertex set (0, 1) and facets {(0, 1)}
              Defn: 0 |--> 0
                    1 |--> 1
        """
        if not self.is_endomorphism_set():
            raise TypeError("identity map is only defined for endomorphism sets")
        f = {i: i for i in self._domain.vertices()}
        return SimplicialComplexMorphism(f, self._domain, self._codomain)

    def an_element(self):
        """
        Return a (non-random) element of ``self``.

        EXAMPLES::

            sage: S = simplicial_complexes.KleinBottle()
            sage: T = simplicial_complexes.Sphere(5)
            sage: H = Hom(S,T)
            sage: x = H.an_element()
            sage: x
            Simplicial complex morphism:
              From: Minimal triangulation of the Klein bottle
              To:   Minimal triangulation of the 5-sphere
              Defn: [0, 1, 2, 3, 4, 5, 6, 7] --> [0, 0, 0, 0, 0, 0, 0, 0]
        """
        X_vertices = self._domain.vertices()
        try:
            i = next(iter(self._codomain.vertices()))
        except StopIteration:
            if not X_vertices:
                return {}
            else:
                raise TypeError("there are no morphisms from a non-empty simplicial complex to an empty simplicial complex")
        f = {x: i for x in X_vertices}
        return SimplicialComplexMorphism(f, self._domain, self._codomain)
