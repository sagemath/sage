r"""
Induced morphisms on homology

This module implements morphisms on homology induced by morphisms of
simplicial complexes. It requires working with field coefficients.

See :class:`InducedHomologyMorphism` for documentation.

AUTHORS:

- John H. Palmieri (2015.09)
"""

########################################################################
#       Copyright (C) 2015 John H. Palmieri <palmieri@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
########################################################################

# To do: implement morphisms of cubical complexes, with methods
#   - domain
#   - codomain
#   - associated_chain_complex_morphism
# Once this is done, the code here ought to work without modification.

from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.graded_modules_with_basis import GradedModulesWithBasis
from sage.categories.morphism import Morphism
from sage.categories.homset import Hom
from sage.rings.rational_field import QQ
from sage.topology.simplicial_complex import SimplicialComplex

class InducedHomologyMorphism(Morphism):
    r"""
    An element of this class is a morphism of (co)homology groups
    induced by a map of simplicial complexes. It requires working
    with field coefficients.

    INPUT:

    - ``map`` -- the map of simplicial complexes
    - ``base_ring`` -- a field (optional, default ``QQ``)
    - ``cohomology`` -- boolean (optional, default ``False``). If
      ``True``, return the induced map in cohomology rather than
      homology.

    .. note::

        This is not intended to be used directly by the user, but instead
        via the method
        :meth:`~sage.topology.simplicial_complex_morphism.SimplicialComplexMorphism.induced_homology_morphism`.

    EXAMPLES::

        sage: S1 = simplicial_complexes.Sphere(1)
        sage: H = Hom(S1, S1)
        sage: f = H({0:0, 1:2, 2:1})  # f switches two vertices
        sage: f_star = f.induced_homology_morphism(QQ, cohomology=True)
        sage: f_star
        Graded algebra endomorphism of Cohomology ring of Minimal triangulation of the 1-sphere over Rational Field
          Defn: induced by:
            Simplicial complex endomorphism of Minimal triangulation of the 1-sphere
              Defn: 0 |--> 0
                    1 |--> 2
                    2 |--> 1
        sage: f_star.to_matrix(1)
        [-1]
        sage: f_star.to_matrix()
        [ 1| 0]
        [--+--]
        [ 0|-1]

        sage: T = simplicial_complexes.Torus()
        sage: y = T.homology_with_basis(QQ).basis()[(1,1)]
        sage: y.to_cycle()
        (0, 5) - (0, 6) + (5, 6)

    Since `(0,2) - (0,5) + (2,5)` is a cycle representing a homology
    class in the torus, we can define a map `S^1 \to T` inducing an
    inclusion on `H_1`::

        sage: Hom(S1, T)({0:0, 1:2, 2:5})
        Simplicial complex morphism:
          From: Minimal triangulation of the 1-sphere
          To: Minimal triangulation of the torus
          Defn: 0 |--> 0
                1 |--> 2
                2 |--> 5
        sage: g = Hom(S1, T)({0:0, 1:2, 2: 5})
        sage: g_star = g.induced_homology_morphism(QQ)
        sage: g_star.to_matrix(0)
        [1]
        sage: g_star.to_matrix(1)
        [-1]
        [ 0]
        sage: g_star.to_matrix()
        [ 1| 0]
        [--+--]
        [ 0|-1]
        [ 0| 0]
        [--+--]
        [ 0| 0]

    We can evaluate such a map on (co)homology classes::

        sage: H = S1.homology_with_basis(QQ)
        sage: a = H.basis()[(1,0)]
        sage: g_star(a)
        -h_{1,0}

        sage: T = S1.product(S1, is_mutable=False)
        sage: diag = Hom(S1,T).diagonal_morphism()
        sage: b,c = list(T.cohomology_ring().basis(1))
        sage: diag_c = diag.induced_homology_morphism(cohomology=True)
        sage: diag_c(b)
        h^{1,0}
        sage: diag_c(c)
        h^{1,0}
    """
    def __init__(self, map, base_ring=None, cohomology=False):
        """
        INPUT:

        - ``map`` -- the map of simplicial complexes
        - ``base_ring`` -- a field (optional, default ``QQ``)
        - ``cohomology`` -- boolean (optional, default ``False``). If
          ``True``, return the induced map in cohomology rather than
          homology.

        EXAMPLES::

            sage: from sage.homology.homology_morphism import InducedHomologyMorphism
            sage: K = simplicial_complexes.RandomComplex(8, 3)
            sage: H = Hom(K,K)
            sage: id = H.identity()
            sage: f = InducedHomologyMorphism(id, QQ)
            sage: f.to_matrix(0) == 1  and  f.to_matrix(1) == 1  and  f.to_matrix(2) == 1
            True
            sage: f = InducedHomologyMorphism(id, ZZ)
            Traceback (most recent call last):
            ...
            ValueError: the coefficient ring must be a field
            sage: S1 = simplicial_complexes.Sphere(1).barycentric_subdivision()
            sage: S1.is_mutable()
            True
            sage: g = Hom(S1, S1).identity()
            sage: h = g.induced_homology_morphism(QQ)
            Traceback (most recent call last):
            ...
            ValueError: the domain and codomain complexes must be immutable
            sage: S1.set_immutable()
            sage: g = Hom(S1, S1).identity()
            sage: h = g.induced_homology_morphism(QQ)
        """
        if (isinstance(map.domain(), SimplicialComplex)
            and (map.domain().is_mutable() or map.codomain().is_mutable())):
                raise ValueError('the domain and codomain complexes must be immutable')
        if base_ring is None:
            base_ring = QQ
        if not base_ring.is_field():
            raise ValueError('the coefficient ring must be a field')

        self._cohomology = cohomology
        self._map = map
        self._base_ring = base_ring
        if cohomology:
            domain = map.codomain().cohomology_ring(base_ring=base_ring)
            codomain = map.domain().cohomology_ring(base_ring=base_ring)
            Morphism.__init__(self, Hom(domain, codomain,
                                        category=GradedAlgebrasWithBasis(base_ring)))
        else:
            domain = map.domain().homology_with_basis(base_ring=base_ring, cohomology=cohomology)
            codomain = map.codomain().homology_with_basis(base_ring=base_ring, cohomology=cohomology)
            Morphism.__init__(self, Hom(domain, codomain,
                                        category=GradedModulesWithBasis(base_ring)))

    def base_ring(self):
        """
        The base ring for this map

        EXAMPLES::

            sage: K = simplicial_complexes.Simplex(2)
            sage: H = Hom(K,K)
            sage: id = H.identity()
            sage: id.induced_homology_morphism(QQ).base_ring()
            Rational Field
            sage: id.induced_homology_morphism(GF(13)).base_ring()
            Finite Field of size 13
        """
        return self._base_ring

    def to_matrix(self, deg=None):
        """
        The matrix for this map.

        If degree ``deg`` is specified, return the matrix just in that
        degree; otherwise, return the block matrix representing the
        entire map.

        INPUT:

        - ``deg`` -- (optional, default ``None``) the degree

        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: S1_b = S1.barycentric_subdivision()
            sage: S1_b.set_immutable()
            sage: d = {(0,): 0, (0,1): 1, (1,): 2, (1,2): 0, (2,): 1, (0,2): 2}
            sage: f = Hom(S1_b, S1)(d)
            sage: h = f.induced_homology_morphism(QQ)
            sage: h.to_matrix(1)
            [2]
            sage: h.to_matrix()
            [1|0]
            [-+-]
            [0|2]
        """
        base_ring = self.base_ring()
        # Compute homology case first.
        domain = self._map.domain()
        codomain = self._map.codomain()
        phi_codomain, H_codomain = codomain.algebraic_topological_model(base_ring)
        phi_domain, H_domain = domain.algebraic_topological_model(base_ring)
        mat = (phi_codomain.pi().to_matrix(deg)
               * self._map.associated_chain_complex_morphism(self.base_ring()).to_matrix(deg)
               * phi_domain.iota().to_matrix(deg))
        if self._cohomology:
            mat = mat.transpose()
            H_domain, H_codomain = H_codomain, H_domain
        if deg is None:
            import numpy as np
            betti_domain = [H_domain.free_module_rank(n)
                            for n in range(domain.dimension()+1)]
            betti_codomain = [H_codomain.free_module_rank(n)
                              for n in range(codomain.dimension()+1)]
            # Compute cumulative sums of Betti numbers to get subdivisions:
            row_subdivs = list(np.cumsum(betti_codomain[:-1]))
            col_subdivs = list(np.cumsum(betti_domain[:-1]))
            mat.subdivide(row_subdivs, col_subdivs)
        return mat

    def __call__(self, elt):
        """
        Evaluate this map on ``elt``, an element of (co)homology.

        INPUT:

        - ``elt`` -- informally, an element of the domain of this
          map. More formally, an element of
          :class:`homology_vector_space_with_basis.HomologyVectorSpaceWithBasis`.

        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: f = {0:0, 1:2, 2:1}
            sage: H = Hom(S1,S1)
            sage: g = H(f)
            sage: h = g.induced_homology_morphism(QQ)
            sage: x = S1.homology_with_basis().basis()[(1,0)]
            sage: x
            h_{1,0}
            sage: h(x)  # indirect doctest
            -h_{1,0}
        """
        base_ring = self.base_ring()
        if self._cohomology:
            codomain = self._map.domain().homology_with_basis(base_ring, cohomology=True)
            if elt.parent().complex() != self._map.codomain():
                raise ValueError('element is not a cohomology class for the correct complex')
        else:
            codomain = self._map.codomain().homology_with_basis(base_ring)
            if elt.parent().complex() != self._map.domain():
                raise ValueError('element is not a homology class for the correct complex')

        return codomain.from_vector(self.to_matrix() * elt.to_vector())

    def __eq__(self, other):
        """
        Return ``True`` if and only if this map agrees with ``other``.

        INPUT:

        - ``other`` -- another induced homology morphism

        This automatically returns ``False`` if the morphisms have
        different domains, codomains, base rings, or values for their
        cohomology flags

        Otherwise, determine this by computing the matrices for this
        map and ``other`` using the (same) basis for the homology
        vector spaces.

        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: K = simplicial_complexes.Simplex(2)
            sage: f = Hom(S1, K)({0: 0, 1:1, 2:2})
            sage: g = Hom(S1, K)({0: 0, 1:0, 2:0})
            sage: f.induced_homology_morphism(QQ) == g.induced_homology_morphism(QQ)
            True
            sage: f.induced_homology_morphism(QQ) == g.induced_homology_morphism(GF(2))
            False
            sage: id = Hom(K, K).identity()   # different domain
            sage: f.induced_homology_morphism(QQ) == id.induced_homology_morphism(QQ)
            False
        """
        if (self._map.domain() != other._map.domain()
            or self._map.codomain() != other._map.codomain()
            or self.base_ring() != other.base_ring()
            or self._cohomology != other._cohomology):
            return False
        dim = min(self._map.domain().dimension(), self._map.codomain().dimension())
        return all(self.to_matrix(d) == other.to_matrix(d) for d in range(dim+1))

    def is_identity(self):
        """
        True if this is the identity map on (co)homology.

        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: H = Hom(S1, S1)
            sage: flip = H({0:0, 1:2, 2:1})
            sage: flip.induced_homology_morphism(QQ).is_identity()
            False
            sage: flip.induced_homology_morphism(GF(2)).is_identity()
            True
            sage: rotate = H({0:1, 1:2, 2:0})
            sage: rotate.induced_homology_morphism(QQ).is_identity()
            True
        """
        return self.to_matrix().is_one()

    def is_surjective(self):
        """
        True if this map is surjective on (co)homology.

        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: K = simplicial_complexes.Simplex(2)
            sage: H = Hom(S1, K)
            sage: f = H({0:0, 1:1, 2:2})
            sage: f.induced_homology_morphism().is_surjective()
            True
            sage: f.induced_homology_morphism(cohomology=True).is_surjective()
            False
        """
        m = self.to_matrix()
        return m.rank() == m.nrows()

    def is_injective(self):
        """
        True if this map is injective on (co)homology.

        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: K = simplicial_complexes.Simplex(2)
            sage: H = Hom(S1, K)
            sage: f = H({0:0, 1:1, 2:2})
            sage: f.induced_homology_morphism().is_injective()
            False
            sage: f.induced_homology_morphism(cohomology=True).is_injective()
            True

            sage: T = simplicial_complexes.Torus()
            sage: g = Hom(S1, T)({0:0, 1:3, 2: 6})
            sage: g_star = g.induced_homology_morphism(QQ)
            sage: g.is_injective()
            True
        """
        return self.to_matrix().right_nullity() == 0

    def _repr_type(self):
        """
        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: K = simplicial_complexes.Simplex(2)
            sage: f = Hom(S1, K)({0: 0, 1:1, 2:2})
            sage: f.induced_homology_morphism()._repr_type()
            'Graded vector space'
            sage: f.induced_homology_morphism(cohomology=True)._repr_type()
            'Graded algebra'
        """
        return "Graded vector space" if not self._cohomology else "Graded algebra"

    def _repr_defn(self):
        """
        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: K = simplicial_complexes.Simplex(2)
            sage: f = Hom(S1, K)({0: 0, 1:1, 2:2})
            sage: print(f.induced_homology_morphism()._repr_defn())
            induced by:
              Simplicial complex morphism:
                From: Minimal triangulation of the 1-sphere
                To:   The 2-simplex
                Defn: 0 |--> 0
                      1 |--> 1
                      2 |--> 2
        """
        s = "induced by:"
        s += '\n  {}'.format('\n  '.join(self._map._repr_().split('\n')))
        return s
