# sage.doctest: needs sage.combinat sage.graphs sage.modules
r"""
Subword complex

Fix a Coxeter system `(W,S)`. The subword complex `\mathcal{SC}(Q,w)`
associated to a word `Q \in S^*` and an element `w \in W` is the
simplicial complex whose ground set is the set of positions in `Q` and whose
facets are complements of sets of positions defining a reduced
expression for `w`.

A subword complex is a shellable sphere if and only if the Demazure
product of `Q` equals `w`, otherwise it is a shellable ball.

The code is optimized to be used with :class:`ReflectionGroup`, it works as well
with :class:`CoxeterGroup`, but many methods fail for :class:`WeylGroup`.

EXAMPLES::

    sage: W = ReflectionGroup(['A',3]); I = list(W.index_set())         # optional - gap3
    sage: Q = I + W.w0.coxeter_sorting_word(I); Q                       # optional - gap3
    [1, 2, 3, 1, 2, 3, 1, 2, 1]

    sage: S = SubwordComplex(Q,W.w0)                                    # optional - gap3
    sage: for F in S: print("{} {}".format(F, F.root_configuration()))                 # optional - gap3
    (0, 1, 2) [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    (0, 1, 8) [(1, 0, 0), (0, 1, 0), (0, 0, -1)]
    (0, 2, 6) [(1, 0, 0), (0, 1, 1), (0, -1, 0)]
    (0, 6, 7) [(1, 0, 0), (0, 0, 1), (0, -1, -1)]
    (0, 7, 8) [(1, 0, 0), (0, -1, 0), (0, 0, -1)]
    (1, 2, 3) [(1, 1, 0), (0, 0, 1), (-1, 0, 0)]
    (1, 3, 8) [(1, 1, 0), (-1, 0, 0), (0, 0, -1)]
    (2, 3, 4) [(1, 1, 1), (0, 1, 0), (-1, -1, 0)]
    (2, 4, 6) [(1, 1, 1), (-1, 0, 0), (0, -1, 0)]
    (3, 4, 5) [(0, 1, 0), (0, 0, 1), (-1, -1, -1)]
    (3, 5, 8) [(0, 1, 0), (-1, -1, 0), (0, 0, -1)]
    (4, 5, 6) [(0, 1, 1), (-1, -1, -1), (0, -1, 0)]
    (5, 6, 7) [(-1, 0, 0), (0, 0, 1), (0, -1, -1)]
    (5, 7, 8) [(-1, 0, 0), (0, -1, 0), (0, 0, -1)]

Testing that the implementation also works with CoxeterGroup::

    sage: W = CoxeterGroup(['A',3]); I = list(W.index_set())
    sage: Q = I + W.w0.coxeter_sorting_word(I); Q
    [1, 2, 3, 1, 2, 3, 1, 2, 1]
    sage: S = SubwordComplex(Q,W.w0); S
    Subword complex of type ['A', 3] for Q = (1, 2, 3, 1, 2, 3, 1, 2, 1) and pi = [1, 2, 3, 1, 2, 1]
    sage: P = S.increasing_flip_poset(); P; len(P.cover_relations())
    Finite poset containing 14 elements
    21

The root configuration works::

    sage: for F in S: print("{} {}".format(F, F.root_configuration()))
    (0, 1, 2) [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    (0, 1, 8) [(1, 0, 0), (0, 1, 0), (0, 0, -1)]
    (0, 2, 6) [(1, 0, 0), (0, 1, 1), (0, -1, 0)]
    (0, 6, 7) [(1, 0, 0), (0, 0, 1), (0, -1, -1)]
    (0, 7, 8) [(1, 0, 0), (0, -1, 0), (0, 0, -1)]
    (1, 2, 3) [(1, 1, 0), (0, 0, 1), (-1, 0, 0)]
    (1, 3, 8) [(1, 1, 0), (-1, 0, 0), (0, 0, -1)]
    (2, 3, 4) [(1, 1, 1), (0, 1, 0), (-1, -1, 0)]
    (2, 4, 6) [(1, 1, 1), (-1, 0, 0), (0, -1, 0)]
    (3, 4, 5) [(0, 1, 0), (0, 0, 1), (-1, -1, -1)]
    (3, 5, 8) [(0, 1, 0), (-1, -1, 0), (0, 0, -1)]
    (4, 5, 6) [(0, 1, 1), (-1, -1, -1), (0, -1, 0)]
    (5, 6, 7) [(-1, 0, 0), (0, 0, 1), (0, -1, -1)]
    (5, 7, 8) [(-1, 0, 0), (0, -1, 0), (0, 0, -1)]

And the weight configuration also works::

    sage: W = CoxeterGroup(['A',2])
    sage: w = W.from_reduced_word([1,2,1])
    sage: SC = SubwordComplex([1,2,1,2,1],w)
    sage: F = SC([1,2])
    sage: F.extended_weight_configuration()
    [(4/3, 2/3), (2/3, 4/3), (-2/3, 2/3), (2/3, 4/3), (-2/3, 2/3)]
    sage: F.extended_weight_configuration(coefficients=(1,2))
    [(4/3, 2/3), (4/3, 8/3), (-2/3, 2/3), (4/3, 8/3), (-2/3, 2/3)]

One finally can compute the brick polytope, using all functionality
on weight configurations, though it does not realize to live in
real space::

    sage: W = CoxeterGroup(['A',3]); I = list(W.index_set())
    sage: Q = I + W.w0.coxeter_sorting_word(I)
    sage: S = SubwordComplex(Q,W.w0)
    sage: S.brick_polytope()                                                            # needs sage.geometry.polyhedron
    A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 14 vertices

    sage: W = CoxeterGroup(['H',3]); I = list(W.index_set())
    sage: Q = I + W.w0.coxeter_sorting_word(I)
    sage: S = SubwordComplex(Q,W.w0)
    sage: S.brick_polytope()                                                            # needs sage.geometry.polyhedron
    doctest:...: RuntimeWarning: the polytope is built with rational vertices
    A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 32 vertices

AUTHORS:

- Christian Stump: initial version
- Vincent Pilaud: greedy flip algorithm, minor improvements, documentation

REFERENCES:

.. [KnuMil] Knutson and Miller. *Subword complexes in Coxeter groups*. Adv. Math., 184(1):161-176, 2004.
.. [PilStu] Pilaud and Stump. *Brick polytopes of spherical subword complexes and generalized associahedra*. Adv. Math. 276:1-61, 2015.
"""
# ****************************************************************************
#       Copyright (C) 2015      Christian Stump <christian.stump@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from itertools import repeat
from copy import copy
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.structure.element import Element
from sage.structure.unique_representation import UniqueRepresentation
from sage.topology.simplicial_complex import SimplicialComplex, Simplex
from sage.categories.simplicial_complexes import SimplicialComplexes
from sage.combinat.subword_complex_c import _flip_c, _construct_facets_c

lazy_import('sage.geometry.polyhedron.constructor', 'Polyhedron')
lazy_import('sage.geometry.cone', 'Cone')


class SubwordComplexFacet(Simplex, Element):
    r"""
    A facet of a subword complex.

    Facets of the subword complex `\mathcal{SC}(Q,w)` are complements
    of sets of positions in `Q` defining a reduced expression for `w`.

    EXAMPLES::

        sage: # optional - gap3
        sage: W = ReflectionGroup(['A',2])
        sage: w = W.from_reduced_word([1,2,1])
        sage: SC = SubwordComplex([1,2,1,2,1], w)
        sage: F = SC[0]; F
        (0, 1)

        sage: W = CoxeterGroup(['A',2])
        sage: w = W.from_reduced_word([1,2,1])
        sage: SC = SubwordComplex([1,2,1,2,1], w)
        sage: F = SC[0]; F
        (0, 1)

    TESTS::

        sage: type(F)                                                   # optional - gap3
        <class 'sage.combinat.subword_complex.SubwordComplex_with_category.element_class'>
    """

    # standard functions

    def __init__(self, parent, positions, facet_test=True):
        r"""
        Initialize a facet of the subword complex ``parent``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])                          # optional - gap3
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)                # optional - gap3
            sage: F = SC([1,2]); F                                      # optional - gap3
            (1, 2)

            sage: W = CoxeterGroup(['A',2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)
            sage: F = SC([1,2]); F
            (1, 2)

        TESTS::

            sage: W = ReflectionGroup(['A',2])                          # optional - gap3
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)                # optional - gap3
            sage: SC([1,3])                                             # optional - gap3
            Traceback (most recent call last):
            ...
            ValueError: the given iterable [1, 3] is not a facet of the Subword complex of type ['A', 2] for Q = (1, 2, 1, 2, 1) and pi = [1, 2, 1]

            sage: W = ReflectionGroup(['A',2])                          # optional - gap3
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)                # optional - gap3
            sage: TestSuite(SC).run()                                   # optional - gap3
        """
        if facet_test and positions not in parent:
            raise ValueError("the given iterable %s is not a facet of the %s" % (positions, parent))
        Simplex.__init__(self, sorted(positions))
        Element.__init__(self, parent)
        self._extended_root_conf_indices = None
        self._extended_weight_conf = None

    # roots

    def _extended_root_configuration_indices(self):
        r"""
        Return the indices of the roots in ``self.group().roots()`` of
        the extended root configuration of ``self``.

        Let `Q = q_1 \dots q_m \in S^*` and `w \in W`. The extended
        root configuration of a facet `I` of `\mathcal{SC}(Q,w)` is
        the sequence `\mathsf{r}(I, 1), \dots, \mathsf{r}(I, m)` of
        roots defined by `\mathsf{r}(I, k) = \Pi Q_{[k-1]
        \smallsetminus I} (\alpha_{q_k})`, where `\Pi Q_{[k-1]
        \smallsetminus I}` is the product of the simple reflections
        `q_i` for `i \in [k-1] \smallsetminus I` in this order.

        .. SEEALSO::

            :meth:`extended_root_configuration`

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F._extended_root_configuration_indices()
            [0, 2, 3, 2, 1]

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F._extended_root_configuration_indices()
            [0, 1, 3, 1, 2]
        """
        if self._extended_root_conf_indices is None:
            self._extended_root_conf_indices = _extended_root_configuration_indices(self.parent().group(), self.parent().word(), self)
        return self._extended_root_conf_indices

    def _root_configuration_indices(self):
        r"""
        Return the indices of the roots in ``self.group().roots()`` of
        the root configuration of ``self``.

        Let `Q = q_1 \dots q_m \in S^*` and `w \in W`. The root
        configuration of a facet `I = [i_1, \dots, i_n]` of
        `\mathcal{SC}(Q,w)` is the sequence `\mathsf{r}(I, i_1),
        \dots, \mathsf{r}(I, i_n)` of roots defined by `\mathsf{r}(I,
        k) = \Pi Q_{[k-1] \smallsetminus I} (\alpha_{q_k})`, where
        `\Pi Q_{[k-1] \smallsetminus I}` is the product of the simple
        reflections `q_i` for `i \in [k-1] \smallsetminus I` in this
        order.

        .. SEEALSO::

            :meth:`root_configuration`

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F._root_configuration_indices()
            [2, 3]

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F._root_configuration_indices()                       # optional - gap3
            [1, 3]
        """
        indices = self._extended_root_configuration_indices()
        return [indices[i] for i in self]

    def extended_root_configuration(self):
        r"""
        Return the extended root configuration of ``self``.

        Let `Q = q_1 \dots q_m \in S^*` and `w \in W`. The extended
        root configuration of a facet `I` of `\mathcal{SC}(Q,w)` is
        the sequence `\mathsf{r}(I, 1), \dots, \mathsf{r}(I, m)` of
        roots defined by `\mathsf{r}(I, k) = \Pi Q_{[k-1]
        \smallsetminus I} (\alpha_{q_k})`, where `\Pi Q_{[k-1]
        \smallsetminus I}` is the product of the simple reflections
        `q_i` for `i \in [k-1] \smallsetminus I` in this order.

        The extended root configuration is used to perform flips efficiently.

        .. SEEALSO::

            :meth:`flip`

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.extended_root_configuration()
            [(1, 0), (1, 1), (-1, 0), (1, 1), (0, 1)]

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.extended_root_configuration()
            [(1, 0), (1, 1), (-1, 0), (1, 1), (0, 1)]
        """
        Phi = self.parent().group().roots()
        return [Phi[i] for i in self._extended_root_configuration_indices()]

    def root_configuration(self):
        r"""
        Return the root configuration of ``self``.

        Let `Q = q_1 \dots q_m \in S^*` and `w \in W`. The root
        configuration of a facet `I = [i_1, \dots, i_n]` of
        `\mathcal{SC}(Q,w)` is the sequence `\mathsf{r}(I, i_1),
        \dots, \mathsf{r}(I, i_n)` of roots defined by `\mathsf{r}(I,
        k) = \Pi Q_{[k-1] \smallsetminus I} (\alpha_{q_k})`, where
        `\Pi Q_{[k-1] \smallsetminus I}` is the product of the simple
        reflections `q_i` for `i \in [k-1] \smallsetminus I` in this
        order.

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.root_configuration()
            [(1, 1), (-1, 0)]

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.root_configuration()                                # optional - gap3
            [(1, 1), (-1, 0)]
        """
        Phi = self.parent().group().roots()
        return [Phi[i] for i in self._root_configuration_indices()]

    def kappa_preimage(self):
        r"""
        Return the fiber of ``self`` under the `\kappa` map.

        The `\kappa` map sends an element `w \in W` to the unique
        facet of `I \in \mathcal{SC}(Q,w)` such that the root
        configuration of `I` is contained in `w(\Phi^+)`.
        In other words, `w` is in the preimage of ``self`` under
        `\kappa` if and only if `w^{-1}` sends every root in the
        root configuration to a positive root.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])                          # optional - gap3
            sage: w = W.from_reduced_word([1,2,1])                      # optional - gap3
            sage: SC = SubwordComplex([1,2,1,2,1],w)                    # optional - gap3

            sage: F = SC([1,2]); F                                      # optional - gap3
            (1, 2)
            sage: F.kappa_preimage()                                    # optional - gap3
            [(1,4)(2,3)(5,6)]

            sage: F = SC([0,4]); F                                      # optional - gap3
            (0, 4)
            sage: F.kappa_preimage()                                    # optional - gap3
            [(1,3)(2,5)(4,6), (1,2,6)(3,4,5)]

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)

            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.kappa_preimage()
            [
            [-1  1]
            [ 0  1]
            ]

            sage: F = SC([0,4]); F
            (0, 4)
            sage: F.kappa_preimage()
            [
            [ 1  0]  [-1  1]
            [ 1 -1], [-1  0]
            ]
        """
        W = self.parent().group()
        N = len(W.long_element(as_word=True))
        root_conf = self._root_configuration_indices()
        return [~w for w in W
                if all(w.action_on_root_indices(i, side='left') < N
                       for i in root_conf)]

    def is_vertex(self):
        r"""
        Return ``True`` if ``self`` is a vertex of the brick polytope
        of ``self.parent``.

        A facet is a vertex of the brick polytope if its root cone is
        pointed. Note that this property is always satisfied for
        root-independent subword complexes.

        .. SEEALSO::

            :meth:`root_cone`

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',1])
            sage: w = W.from_reduced_word([1])
            sage: SC = SubwordComplex([1,1,1],w)
            sage: F = SC([0,1]); F.is_vertex()
            True
            sage: F = SC([0,2]); F.is_vertex()
            False

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1,2,1],w)
            sage: F = SC([0,1,2,3]); F.is_vertex()
            True
            sage: F = SC([0,1,2,6]); F.is_vertex()
            False

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1,2,1],w)
            sage: F = SC([0,1,2,3]); F.is_vertex()
            True
            sage: F = SC([0,1,2,6]); F.is_vertex()
            False
        """
        S = self.parent()
        if S.is_root_independent():
            return True
        return self.root_cone().is_strictly_convex()

    @cached_method
    def root_cone(self):
        r"""
        Return the polyhedral cone generated by the root configuration
        of ``self``.

        .. SEEALSO::

            :meth:`root_configuration`

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',1])
            sage: w = W.from_reduced_word([1])
            sage: SC = SubwordComplex([1,1,1],w)
            sage: F = SC([0,2]); F.root_cone()
            1-d cone in 1-d lattice N

            sage: W = CoxeterGroup(['A',1])
            sage: w = W.from_reduced_word([1])
            sage: SC = SubwordComplex([1,1,1],w)
            sage: F = SC([0,2]); F.root_cone()
            1-d cone in 1-d lattice N
        """
        return Cone(self.root_configuration())

    def upper_root_configuration(self):
        r"""
        Return the positive roots of the root configuration of ``self``.

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.root_configuration()
            [(1, 1), (-1, 0)]
            sage: F.upper_root_configuration()
            [(1, 0)]

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.upper_root_configuration()
            [(1, 0)]
        """
        conf = self._root_configuration_indices()
        W = self.parent().group()
        Phi = W.roots()
        N = len(Phi) // 2
        return [Phi[i - N] for i in conf if i >= N]

    # weights

    def extended_weight_configuration(self, coefficients=None):
        r"""
        Return the extended weight configuration of ``self``.

        Let `Q = q_1 \dots q_m \in S^*` and `w \in W`. The extended
        weight configuration of a facet `I` of `\mathcal{SC}(Q,w)` is
        the sequence `\mathsf{w}(I, 1), \dots, \mathsf{w}(I, m)` of
        weights defined by `\mathsf{w}(I, k) = \Pi Q_{[k-1]
        \smallsetminus I} (\omega_{q_k})`, where `\Pi Q_{[k-1]
        \smallsetminus I}` is the product of the simple reflections
        `q_i` for `i \in [k-1] \smallsetminus I` in this order.

        The extended weight configuration is used to compute the brick vector.

        INPUT:

        - ``coefficients`` -- (optional) a list of coefficients used to
          scale the fundamental weights

        .. SEEALSO::

            :meth:`brick_vector`

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2])
            sage: F.extended_weight_configuration()
            [(2/3, 1/3), (1/3, 2/3), (-1/3, 1/3), (1/3, 2/3), (-1/3, 1/3)]
            sage: F.extended_weight_configuration(coefficients=(1,2))
            [(2/3, 1/3), (2/3, 4/3), (-1/3, 1/3), (2/3, 4/3), (-1/3, 1/3)]

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2])
            sage: F.extended_weight_configuration()
            [(4/3, 2/3), (2/3, 4/3), (-2/3, 2/3), (2/3, 4/3), (-2/3, 2/3)]
            sage: F.extended_weight_configuration(coefficients=(1,2))
            [(4/3, 2/3), (4/3, 8/3), (-2/3, 2/3), (4/3, 8/3), (-2/3, 2/3)]
        """
        if coefficients is not None or self._extended_weight_conf is None:
            W = self.parent().group()
            I = W.index_set()
            Lambda = W.fundamental_weights()
            if coefficients is not None:
                coeff = {I[i]: coefficients[i]
                         for i in range(len(coefficients))}
                Lambda = {li: coeff[li] * Lambda[li] for li in Lambda.keys()}
            Q = self.parent().word()
            V_weights = []
            pi = W.one()
            for i, wi in enumerate(Q):
                fund_weight = Lambda[wi]
                V_weights.append(pi * fund_weight)
                if i not in self:
                    pi = pi.apply_simple_reflection_right(wi)
            if self._extended_weight_conf is None:
                self._extended_weight_conf = V_weights
            return V_weights
        else:
            return self._extended_weight_conf

    def weight_configuration(self):
        r"""
        Return the weight configuration of ``self``.

        Let `Q = q_1 \dots q_m \in S^*` and `w \in W`. The weight
        configuration of a facet `I = [i_1, \dots, i_n]` of
        `\mathcal{SC}(Q,w)` is the sequence `\mathsf{w}(I, i_1),
        \dots, \mathsf{w}(I, i_n)` of weights defined by
        `\mathsf{w}(I, k) = \Pi Q_{[k-1] \smallsetminus I}
        (\omega_{q_k})`, where `\Pi Q_{[k-1] \smallsetminus I}` is the
        product of the simple reflections `q_i` for `i \in [k-1]
        \smallsetminus I` in this order.

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.weight_configuration()
            [(1/3, 2/3), (-1/3, 1/3)]

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.weight_configuration()
            [(2/3, 4/3), (-2/3, 2/3)]
        """
        extended_configuration = self.extended_weight_configuration()
        return [extended_configuration[i] for i in self]

    @cached_method
    def weight_cone(self):
        r"""
        Return the polyhedral cone generated by the weight
        configuration of ``self``.

        .. SEEALSO::

            :meth:`weight_configuration`

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: WC = F.weight_cone(); WC
            2-d cone in 2-d lattice N
            sage: WC.rays()
            N( 1, 2),
            N(-1, 1)
            in 2-d lattice N

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: WC = F.weight_cone(); WC
            2-d cone in 2-d lattice N
        """
        return Cone(self.weight_configuration())

    def brick_vector(self, coefficients=None):
        r"""
        Return the brick vector of ``self``.

        This is the sum of the weight vectors in the extended weight
        configuration.

        INPUT:

        - ``coefficients`` -- (optional) a list of coefficients used to
          scale the fundamental weights

        .. SEEALSO::

            :meth:`extended_weight_configuration`

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.extended_weight_configuration()
            [(2/3, 1/3), (1/3, 2/3), (-1/3, 1/3), (1/3, 2/3), (-1/3, 1/3)]
            sage: F.brick_vector()
            (2/3, 7/3)
            sage: F.brick_vector(coefficients=[1,2])
            (4/3, 11/3)

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2])
            sage: F.brick_vector()
            (4/3, 14/3)
            sage: F.brick_vector(coefficients=[1,2])
            (8/3, 22/3)
        """
        return sum(self.extended_weight_configuration(coefficients=coefficients))

    # flip

    def flip(self, i, return_position=False):
        r"""
        Return the facet obtained after flipping position ``i`` in ``self``.

        INPUT:

        - ``i`` -- integer; position in the word `Q`
        - ``return_position`` -- boolean (default: ``False``); tells
          whether the new position should be returned as well

        OUTPUT:

        - The new subword complex facet.
        - The new position if ``return_position`` is ``True``.

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.flip(1)
            (2, 3)
            sage: F.flip(1, return_position=True)
            ((2, 3), 3)

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F
            (1, 2)
            sage: F.flip(1)
            (2, 3)
            sage: F.flip(1, return_position=True)
            ((2, 3), 3)
        """
        S = self.parent()
        F = set(self)
        R = list(self._extended_root_configuration_indices())
        j = _flip_c(self.parent().group(), F, R, i)  # F and R are changed here
        new_facet = S.element_class(self.parent(), F)
        new_facet._extended_root_conf_indices = tuple(R)
        if return_position:
            return new_facet, j
        else:
            return new_facet

    # plot and show

    def plot(self, list_colors=None, labels=[], thickness=3, fontsize=14,
             shift=(0, 0), compact=False, roots=True, **args):
        r"""
        In type `A` or `B`, plot a pseudoline arrangement representing
        the facet ``self``.

        Pseudoline arrangements are graphical representations of
        facets of types A or B subword complexes.

        INPUT:

        - ``list_colors`` -- list (default: ``[]``); to change the colors
          of the pseudolines
        - ``labels`` -- list (default: ``[]``); to change the labels
          of the pseudolines
        - ``thickness`` -- integer (default: 3); for the thickness
          of the pseudolines
        - ``fontsize`` -- integer (default: 14); for the size
          of the font used for labels
        - ``shift`` -- couple of coordinates (default: ``(0,0)``)
          to change the origin
        - ``compact`` -- boolean (default: ``False``); to require
          a more compact representation
        - ``roots`` -- boolean (default: ``True``); whether to print
          the extended root configuration

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F.plot()                                               # needs sage.plot
            Graphics object consisting of 26 graphics primitives

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F.plot()                                               # needs sage.plot
            Graphics object consisting of 26 graphics primitives

            sage: # optional - gap3
            sage: W = ReflectionGroup(['B',3])
            sage: c = W.from_reduced_word([1,2,3])
            sage: Q = c.reduced_word()*2 + W.w0.coxeter_sorting_word(c)
            sage: SC = SubwordComplex(Q, W.w0)
            sage: F = SC[15]; F.plot()                                                  # needs sage.plot
            Graphics object consisting of 53 graphics primitives

        TESTS::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['D',4])
            sage: c = W.from_reduced_word([1,2,3,4])
            sage: Q = c.reduced_word() + W.w0.coxeter_sorting_word(c)
            sage: SC = SubwordComplex(Q, W.w0)
            sage: F = SC[1]; F.plot()                                                   # needs sage.plot
            Traceback (most recent call last):
            ...
            ValueError: plotting is currently only implemented for irreducibles types A, B, and C.

            sage: W = CoxeterGroup(CoxeterMatrix((['A',2],['A',2])))
            sage: c = W.from_reduced_word([1,2,3,4])
            sage: Q = c.reduced_word() + W.w0.coxeter_sorting_word(c)
            sage: SC = SubwordComplex(Q, W.w0)
            sage: F = SC[1]; F.plot()                                                   # needs sage.plot
            Traceback (most recent call last):
            ...
            ValueError: plotting is currently only implemented for irreducibles types A, B, and C.

        REFERENCES: [PilStu]_
        """
        # check that the type is A or B
        # TODO in a better way
        S = self.parent()
        Q = S.word()
        W = S.group()
        n = W.rank()

        error_msg = "plotting is currently only implemented for irreducibles types A, B, and C."
        if S._cartan_type is not None:
            cartan_type = S._cartan_type
            type = cartan_type.type()
            G = cartan_type.coxeter_matrix().coxeter_graph()
        else:
            type = None

        if type not in ['A', 'B', 'C'] or not G.is_connected():
            raise ValueError(error_msg)

        # organization of the indexing
        # TODO: this might be better done in CoxeterType directly.
        index_set = None
        for a in G.vertex_iterator():
            if G.degree(a) == 1:
                b = G.neighbors(a)[0]
                if type == "A" or G.edge_label(a, b) == 4:
                    index_set = [a, b]
                    break
        assert index_set is not None, "Bug in the plot method"
        while G.degree(b) == 2:
            for c in G.neighbors(b):
                # picking the other neighbors of b
                if c != a:
                    index_set.append(c)
                    a = b
                    b = c
                    break

        # import plot facilities
        from sage.plot.line import line
        from sage.plot.text import text
        from sage.plot.colors import colors
        from sage.combinat.permutation import Permutation

        # get properties
        if type == 'A':
            last = n
        else:
            last = n - 1
        permutation = Permutation(range(1, last + 2))
        x_max = .5

        # list the pseudolines to be drawn
        pseudolines = [[(shift[0], shift[1] + i), .5] for i in range(last + 1)]
        pseudolines_type_B = [[] for _ in repeat(None, last + 1)]
        contact_points = []
        root_labels = []
        pseudoline_labels = []
        if labels is not False:
            pseudoline_labels += [(pseudoline,
                                   (shift[0] - .1, shift[1] + pseudoline),
                                   "center") for pseudoline in range(last + 1)]
        if roots:
            extended_root_conf = self.extended_root_configuration()
        for position in range(len(Q)):
            y = index_set.index(Q[position])
            if type in ['B', 'C'] and y == 0:
                pseudoline = permutation(1) - 1
                x = pseudolines[pseudoline].pop()
                if compact:
                    x_max = max(x + 1, x_max)
                else:
                    x = x_max
                    x_max += 1
                if position in self:
                    pseudolines[pseudoline] += [(shift[0] + x + 1,
                                                 shift[1]), x + 1]
                    contact_points += [[(shift[0] + x + .5, shift[1] - .2),
                                        (shift[0] + x + .5, shift[1])]]
                else:
                    pseudolines_type_B[pseudoline] = pseudolines[pseudoline] + [(shift[0] + x + .5, shift[1]), (shift[0] + x + .5, shift[1] - .2)]
                    pseudolines[pseudoline] = [(shift[0] + x + .6, shift[1] - .2), (shift[0] + x + .6, shift[1]), .5]
                if roots:
                    root_labels.append((extended_root_conf[position],
                                        (shift[0] + x + .25, shift[1] - .2)))
            else:
                if type in ['B', 'C']:
                    y -= 1
                pseudoline1 = permutation(y + 1) - 1
                pseudoline2 = permutation(y + 2) - 1
                x = max(pseudolines[pseudoline1].pop(),
                        pseudolines[pseudoline2].pop())
                if compact:
                    x_max = max(x + 1, x_max)
                else:
                    x = x_max
                    x_max += 1
                if position in self:
                    pseudolines[pseudoline1] += [(shift[0] + x + 1,
                                                  shift[1] + y), x + 1]
                    pseudolines[pseudoline2] += [(shift[0] + x + 1,
                                                  shift[1] + y + 1), x + 1]
                    contact_points += [[(shift[0] + x + .5, shift[1] + y),
                                        (shift[0] + x + .5, shift[1] + y + 1)]]
                else:
                    pseudolines[pseudoline1] += [(shift[0] + x + .6,
                                                  shift[1] + y),
                                                 (shift[0] + x + .6,
                                                  shift[1] + y + 1), x + 1]
                    pseudolines[pseudoline2] += [(shift[0] + x + .5,
                                                  shift[1] + y + 1),
                                                 (shift[0] + x + .5,
                                                  shift[1] + y), x + 1]
                    permutation = permutation._left_to_right_multiply_on_left(Permutation((y + 1, y + 2)))
                if roots:
                    root_labels.append((extended_root_conf[position],
                                        (shift[0] + x + .35,
                                         shift[1] + y + .5)))
                if labels is not False:
                    pseudoline_labels += [(pseudoline1, (shift[0] + x + .35,
                                                         shift[1] + y + .05),
                                           "bottom"),
                                          (pseudoline2, (shift[0] + x + .35,
                                                         shift[1] + y + .95),
                                           "top")]

        # transform list to real lines
        if list_colors is None:
            list_colors = []
        list_colors += ['red', 'blue', 'green', 'orange', 'yellow', 'purple']
        list_colors += list(colors)
        thickness = max(thickness, 2)
        L = line([(1, 1)])
        for contact_point in contact_points:
            L += line(contact_point, rgbcolor=[0, 0, 0],
                      thickness=thickness - 1)
        for pseudoline in range(last + 1):
            pseudolines[pseudoline].pop()
            pseudolines[pseudoline].append((shift[0] + x_max,
                                            shift[1] + permutation.inverse()(pseudoline + 1) - 1))
            L += line(pseudolines[pseudoline], color=list_colors[pseudoline],
                      thickness=thickness)
            if type in ['B', 'C']:
                L += line(pseudolines_type_B[pseudoline],
                          color=list_colors[pseudoline],
                          thickness=thickness, linestyle='--')
        for root_label in root_labels:
            L += text(root_label[0], root_label[1], rgbcolor=[0, 0, 0],
                      fontsize=fontsize, vertical_alignment='center',
                      horizontal_alignment='right')
        if len(labels) < last + 1:
            labels = list(range(1, last + 2))
        for pseudoline_label in pseudoline_labels:
            L += text(labels[pseudoline_label[0]], pseudoline_label[1],
                      color=list_colors[pseudoline_label[0]],
                      fontsize=fontsize,
                      vertical_alignment=pseudoline_label[2],
                      horizontal_alignment='right')
        if labels is not False:
            for pseudoline in range(last):
                L += text(labels[pseudoline],
                          (shift[0] + x_max + .1,
                           shift[1] + permutation.inverse()(pseudoline + 1) - 1),
                          color=list_colors[pseudoline], fontsize=fontsize,
                          vertical_alignment='center',
                          horizontal_alignment='left')
        L.axes(False)
        return L

    def show(self, *kwds, **args):
        """
        Show the facet ``self``.

        .. SEEALSO::

            :meth:`plot`

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1],w)
            sage: F = SC([1,2]); F.show()
            <BLANKLINE>
        """
        return self.plot().show(*kwds, **args)


class SubwordComplex(UniqueRepresentation, SimplicialComplex):
    r"""
    Fix a Coxeter system `(W,S)`. The subword complex
    `\mathcal{SC}(Q,w)` associated to a word `Q \in S^*` and an
    element `w \in W` is the simplicial complex whose ground set is the set of
    positions in `Q` and whose facets are complements of sets of
    positions defining a reduced expression for `w`.

    A subword complex is a shellable sphere if and only if the
    Demazure product of `Q` equals `w`, otherwise it is a shellable
    ball.

    .. WARNING::

        This implementation only works for groups build using ``CoxeterGroup``,
        and does not work with groups build using ``WeylGroup``.

    EXAMPLES:

    As an example, dual associahedra are subword complexes in type
    `A_{n-1}` given by the word `[1, \dots, n, 1, \dots, n, 1, \dots,
    n-1, \dots, 1, 2, 1]` and the permutation `w_0`.

    ::

        sage: # optional - gap3
        sage: W = ReflectionGroup(['A',2])
        sage: w = W.from_reduced_word([1,2,1])
        sage: SC = SubwordComplex([1,2,1,2,1], w); SC
        Subword complex of type ['A', 2] for Q = (1, 2, 1, 2, 1) and pi = [1, 2, 1]
        sage: SC.facets()
        [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]

        sage: W = CoxeterGroup(['A',2])
        sage: w = W.from_reduced_word([1,2,1])
        sage: SC = SubwordComplex([1,2,1,2,1], w); SC
        Subword complex of type ['A', 2] for Q = (1, 2, 1, 2, 1) and pi = [1, 2, 1]
        sage: SC.facets()
        [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]

    REFERENCES: [KnuMil]_, [PilStu]_

    TESTS::

        sage: # optional - gap3
        sage: W = ReflectionGroup(['A',2])
        sage: w = W.from_reduced_word([1,2,1])
        sage: SC1 = SubwordComplex([1,2,1,2,1], w)
        sage: SC2 = SubwordComplex([1,2,1,2,1], w)
        sage: SC1 == SC2
        True

        sage: W = CoxeterGroup(['A',2])
        sage: w = W.from_reduced_word([1,2,1])
        sage: SC1 = SubwordComplex([1,2,1,2,1], w)
        sage: SC2 = SubwordComplex([1,2,1,2,1], w)
        sage: SC1 == SC2
        True
    """

    # standard functions

    @staticmethod
    def __classcall__(cls, Q, w, algorithm='inductive'):
        r"""
        Making the input hashable.

        TESTS::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['B',2])
            sage: S = SubwordComplex((1,2)*3,W.w0)
            sage: T = SubwordComplex([1,2]*3,W.w0)
            sage: S is T
            True

            sage: W = CoxeterGroup(['B',2])
            sage: S = SubwordComplex((1,2)*3,W.w0)
            sage: T = SubwordComplex([1,2]*3,W.w0)
            sage: S is T
            True
        """
        Q = tuple(Q)
        return super().__classcall__(cls, Q, w, algorithm=algorithm)

    def __init__(self, Q, w, algorithm='inductive'):
        r"""
        Initialize the subword complex `\mathcal{SC}(Q,w)`.

        INPUT:

        - ``Q`` -- word on the simple generators of the Coxeter group
        - ``w`` -- element of the Coxeter group
        - ``algorithm`` -- (default: ``'inductive'``) choice of the
          algorithm to generate the subword complex. Options are
          ``'inductive'`` or ``'greedy'``. The second option is
          recommended when `|Q|` is closed to `\ell(w) + \mathrm{rank}(W)`.

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',3])
            sage: w = W.from_reduced_word([1,2,3,1,2,1])
            sage: SC = SubwordComplex([1,2,3,1,2,3,1,2,1], w); SC
            Subword complex of type ['A', 3] for Q = (1, 2, 3, 1, 2, 3, 1, 2, 1) and pi = [1, 2, 1, 3, 2, 1]
            sage: len(SC)
            14

            sage: W = CoxeterGroup(['A',3])
            sage: w = W.from_reduced_word([1,2,3,1,2,1])
            sage: SC = SubwordComplex([1,2,3,1,2,3,1,2,1], w); SC
            Subword complex of type ['A', 3] for Q = (1, 2, 3, 1, 2, 3, 1, 2, 1) and pi = [1, 2, 3, 1, 2, 1]
            sage: len(SC)
            14

        TESTS:

        Check for methods from the enumerated sets category::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: list(SC)
            [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: list(SC)
            [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,1,1])
            sage: SC = SubwordComplex([1,2,2,2,1], w)
            sage: len(SC)
            2
        """
        W = w.parent()
        I = W.index_set()
        if not all(i in I for i in Q):
            raise ValueError("all elements in Q = %s must be contained in the index set %s" % (Q, I))
        self._Q = Q
        self._pi = w
        if algorithm == "inductive":
            Fs = _construct_facets_c(Q, w)
        elif algorithm == "greedy":
            Fs, Rs = _greedy_flip_algorithm(Q, w)
        else:
            raise ValueError("the optional argument algorithm can be "
                             "either inductive or greedy")
        if not Fs:
            raise ValueError("the word %s does not contain a reduced expression for %s" % (Q, w.reduced_word()))
        cat = SimplicialComplexes().Finite().Enumerated()
        SimplicialComplex.__init__(self, maximal_faces=Fs,
                                   maximality_check=False,
                                   category=cat)
        self._W = W
        try:
            T = W.coxeter_matrix().coxeter_type()
            self._cartan_type = T.cartan_type()
        except AttributeError:
            self._cartan_type = None
        self._facets_dict = None
        if algorithm == "greedy":
            _facets_dict = {}
            for i in range(len(Fs)):
                X = self(Fs[i], facet_test=False)
                X._extended_root_conf_indices = Rs[i]
                _facets_dict[tuple(sorted(Fs[i]))] = X
            self._facets_dict = _facets_dict
        else:
            self._facets_dict = {}

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])                          # optional - gap3
            sage: w = W.from_reduced_word([1,2,1])                      # optional - gap3
            sage: SubwordComplex([1,2,1,2,1], w)                        # optional - gap3
            Subword complex of type ['A', 2] for Q = (1, 2, 1, 2, 1) and pi = [1, 2, 1]

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SubwordComplex([1,2,1,2,1], w)
            Subword complex of type ['A', 2] for Q = (1, 2, 1, 2, 1) and pi = [1, 2, 1]
        """
        if self._cartan_type is None:
            return "Subword complex of unknown type for Q = {} and pi = {}".format(self._Q, self._pi.reduced_word())
        else:
            return 'Subword complex of type {} for Q = {} and pi = {}'.format(self.cartan_type(), self._Q, self._pi.reduced_word())

    def __call__(self, F, facet_test=True):
        r"""
        Create a facet of ``self``.

        INPUT:

        - ``F`` -- an iterable of positions
        - ``facet_test`` -- boolean (default: ``True``); whether or
          not the facet ``F`` should be tested before creation

        OUTPUT: the facet of ``self`` at positions given by ``F``

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])                          # optional - gap3
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)                # optional - gap3
            sage: F = SC([1,2]); F                                      # optional - gap3
            (1, 2)

            sage: W = CoxeterGroup(['A',2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)
            sage: F = SC([1,2]); F
            (1, 2)
        """
        if hasattr(F, "parent") and F.parent() is self:
            return F
        return self.element_class(self, F, facet_test=facet_test)

    Element = SubwordComplexFacet

    def __contains__(self, F):
        r"""
        Test if ``self`` contains a given iterable ``F``.

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w  = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.facets()
            [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]
            sage: [0,1] in SC
            True
            sage: [0,2] in SC
            False
            sage: [0,1,5] in SC
            False
            sage: [0] in SC
            False
            sage: ['a','b'] in SC
            False

            sage: W = CoxeterGroup(['A',2])
            sage: w  = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.facets()
            [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]
            sage: [0,1] in SC
            True
            sage: [0,2] in SC
            False
            sage: [0,1,5] in SC
            False
            sage: [0] in SC
            False
            sage: ['a','b'] in SC
            False
        """
        W = self.group()
        Q = self.word()
        r = range(len(Q))
        if not all(i in r for i in F):
            return False
        return W.from_reduced_word(Qi for i, Qi in enumerate(Q) if i not in F) == self.pi()

    # getting the stored properties

    def group(self):
        r"""
        Return the group associated to ``self``.

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.group()
            Irreducible real reflection group of rank 2 and type A2

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.group()
            Finite Coxeter group over Integer Ring with Coxeter matrix:
            [1 3]
            [3 1]
        """
        return self._W

    def cartan_type(self):
        r"""
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.cartan_type()
            ['A', 2]

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.cartan_type()
            ['A', 2]
        """
        if self._cartan_type is None:
            raise ValueError("no Cartan type defined for {}".format(self._W))
        else:
            return self._cartan_type

    def word(self):
        r"""
        Return the word in the simple generators associated to ``self``.

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.word()
            (1, 2, 1, 2, 1)

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.word()
            (1, 2, 1, 2, 1)
        """
        return copy(self._Q)

    def pi(self):
        r"""
        Return the element in the Coxeter group associated to ``self``.

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.pi().reduced_word()
            [1, 2, 1]

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.pi().reduced_word()
            [1, 2, 1]
        """
        return self._pi

    def facets(self):
        r"""
        Return all facets of ``self``.

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.facets()
            [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.facets()
            [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]
        """
        if self._facets_dict:
            return [self._facets_dict[tuple(F)] for F in self._facets]
        else:
            return [self(F, facet_test=False) for F in self._facets]

    def __iter__(self):
        r"""
        Return an iterator on the facets of ``self``.

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: for I in SC: print(I)
            (0, 1)
            (0, 4)
            (1, 2)
            (2, 3)
            (3, 4)

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: for I in SC: print(I)
            (0, 1)
            (0, 4)
            (1, 2)
            (2, 3)
            (3, 4)
        """
        return iter(self.facets())

    def greedy_facet(self, side='positive'):
        r"""
        Return the negative (or positive) greedy facet of ``self``.

        This is the lexicographically last (or first) facet of ``self``.

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.greedy_facet(side='positive')
            (0, 1)
            sage: SC.greedy_facet(side='negative')
            (3, 4)

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.greedy_facet(side='positive')
            (0, 1)
            sage: SC.greedy_facet(side='negative')
            (3, 4)
        """
        return self.element_class(self, _greedy_facet(self.word(),
                                                      self.pi(), side=side))

    # topological properties

    def is_sphere(self):
        r"""
        Return ``True`` if the subword complex ``self`` is a sphere.

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',3])
            sage: w = W.from_reduced_word([2,3,2])
            sage: SC = SubwordComplex([3,2,3,2,3], w)
            sage: SC.is_sphere()
            True

            sage: SC = SubwordComplex([3,2,1,3,2,3], w)                 # optional - gap3
            sage: SC.is_sphere()                                        # optional - gap3
            False

            sage: W = CoxeterGroup(['A',3])
            sage: w = W.from_reduced_word([2,3,2])
            sage: SC = SubwordComplex([3,2,3,2,3], w)
            sage: SC.is_sphere()
            True
        """
        W = self._pi.parent()
        w = W.demazure_product(self._Q)
        return w == self._pi

    def is_ball(self):
        r"""
        Return ``True`` if the subword complex ``self`` is a ball.

        This is the case if and only if it is not a sphere.

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',3])
            sage: w = W.from_reduced_word([2,3,2])
            sage: SC = SubwordComplex([3,2,3,2,3], w)
            sage: SC.is_ball()
            False

            sage: SC = SubwordComplex([3,2,1,3,2,3], w)                 # optional - gap3
            sage: SC.is_ball()                                          # optional - gap3
            True

            sage: W = CoxeterGroup(['A',3])
            sage: w = W.from_reduced_word([2,3,2])
            sage: SC = SubwordComplex([3,2,3,2,3], w)
            sage: SC.is_ball()
            False
        """
        return not self.is_sphere()

    def is_pure(self):
        r"""
        Return ``True`` since all subword complexes are pure.

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',3])
            sage: w = W.from_reduced_word([2,3,2])
            sage: SC = SubwordComplex([3,2,3,2,3], w)
            sage: SC.is_pure()
            True

            sage: W = CoxeterGroup(['A',3])
            sage: w = W.from_reduced_word([2,3,2])
            sage: SC = SubwordComplex([3,2,3,2,3], w)
            sage: SC.is_pure()
            True
        """
        return True

    def dimension(self):
        r"""
        Return the dimension of ``self``.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])                          # optional - gap3
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)                # optional - gap3
            sage: SC.dimension()                                        # optional - gap3
            1

            sage: W = CoxeterGroup(['A',2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)
            sage: SC.dimension()
            1
        """
        return self._facets[0].dimension()

    # root and weight

    @cached_method
    def is_root_independent(self):
        r"""
        Return ``True`` if ``self`` is root-independent.

        This means that the root configuration
        of any (or equivalently all) facets is linearly independent.

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])                          # optional - gap3
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)                # optional - gap3
            sage: SC.is_root_independent()                              # optional - gap3
            True

            sage: SC = SubwordComplex([1,2,1,2,1,2], W.w0)              # optional - gap3
            sage: SC.is_root_independent()                              # optional - gap3
            False

            sage: W = CoxeterGroup(['A',2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)
            sage: SC.is_root_independent()
            True
        """
        from sage.matrix.constructor import matrix
        M = matrix(self.greedy_facet(side='negative').root_configuration())
        return M.rank() == max(M.ncols(), M.nrows())

    @cached_method
    def is_double_root_free(self):
        r"""
        Return ``True`` if ``self`` is double-root-free.

        This means that the root configurations
        of all facets do not contain a root twice.

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.is_double_root_free()
            True

            sage: SC = SubwordComplex([1,1,2,2,1,1], w)                 # optional - gap3
            sage: SC.is_double_root_free()                              # optional - gap3
            True

            sage: SC = SubwordComplex([1,2,1,2,1,2], w)                 # optional - gap3
            sage: SC.is_double_root_free()                              # optional - gap3
            False

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.is_double_root_free()
            True
        """
        if not self.is_root_independent():
            size = self.dimension() + 1
            for F in self:
                conf = F._root_configuration_indices()
                if len(set(conf)) < size:
                    return False
        return True

    def kappa_preimages(self):
        """
        Return a dictionary containing facets of ``self`` as keys,
        and list of elements of ``self.group()`` as values.

        .. SEEALSO::

            :meth:`kappa_preimage <sage.combinat.subword_complex.SubwordComplexFacet.kappa_preimage>`

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: kappa = SC.kappa_preimages()
            sage: for F in SC: print("{} {}".format(F, [w.reduced_word() for w in kappa[F]]))
            (0, 1) [[]]
            (0, 4) [[2], [2, 1]]
            (1, 2) [[1]]
            (2, 3) [[1, 2]]
            (3, 4) [[1, 2, 1]]

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: kappa = SC.kappa_preimages()
            sage: for F in SC: print("{} {}".format(F, [w.reduced_word() for w in kappa[F]]))
            (0, 1) [[]]
            (0, 4) [[2], [2, 1]]
            (1, 2) [[1]]
            (2, 3) [[1, 2]]
            (3, 4) [[1, 2, 1]]
        """
        return {F: F.kappa_preimage() for F in self}

    def brick_fan(self):
        r"""
        Return the brick fan of ``self``.

        It is the normal fan of the brick polytope of ``self``. It is
        formed by the cones generated by the weight configurations of
        the facets of ``self``.

        .. SEEALSO::

            :func:`weight_cone <sage.combinat.subword_complex.SubwordComplexFacet.weight_cone>`

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.brick_fan()
            Rational polyhedral fan in 2-d lattice N

            sage: W = CoxeterGroup(['A',2])
            sage: w = W.from_reduced_word([1,2,1])
            sage: SC = SubwordComplex([1,2,1,2,1], w)
            sage: SC.brick_fan()
            Rational polyhedral fan in 2-d lattice N
        """
        from sage.geometry.fan import Fan
        return Fan([F.weight_cone() for F in self])

    # brick polytope

    def brick_vectors(self, coefficients=None):
        r"""
        Return the list of all brick vectors of facets of ``self``.

        INPUT:

        - ``coefficients`` -- (optional) a list of coefficients used to
          scale the fundamental weights

        .. SEEALSO::

            :func:`brick_vector <sage.combinat.subword_complex.SubwordComplexFacet.brick_vector>`

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)
            sage: SC.brick_vectors()
            [(5/3, 7/3), (5/3, 1/3), (2/3, 7/3), (-1/3, 4/3), (-1/3, 1/3)]
            sage: SC.brick_vectors(coefficients=(1,2))
            [(7/3, 11/3), (7/3, 2/3), (4/3, 11/3), (-2/3, 5/3), (-2/3, 2/3)]

            sage: W = CoxeterGroup(['A',2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)
            sage: SC.brick_vectors()
            [(10/3, 14/3), (10/3, 2/3), (4/3, 14/3), (-2/3, 8/3), (-2/3, 2/3)]
            sage: SC.brick_vectors(coefficients=(1,2))
            [(14/3, 22/3), (14/3, 4/3), (8/3, 22/3), (-4/3, 10/3), (-4/3, 4/3)]
        """
        return [F.brick_vector(coefficients=coefficients) for F in self]

    def minkowski_summand(self, i):
        r"""
        Return the `i` th Minkowski summand of ``self``.

        INPUT:

        - ``i`` -- an integer defining a position in the word `Q`

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])                          # optional - gap3
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)                # optional - gap3
            sage: SC.minkowski_summand(1)                               # optional - gap3
            A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex

            sage: W = CoxeterGroup(['A',2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)
            sage: SC.minkowski_summand(1)
            A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex
        """
        G = self.group()
        from sage.rings.rational_field import QQ
        if G.coxeter_matrix().is_crystallographic():
            min_sum = [[QQ(v) for v in F.extended_weight_configuration()[i]] for F in self]
        else:
            from sage.rings.cc import CC
            from warnings import warn
            warn("the polytope is built with rational vertices", RuntimeWarning)
            min_sum = [[QQ(CC(v)) for v in F.extended_weight_configuration()[i]] for F in self]
        return Polyhedron(min_sum)

    def brick_polytope(self, coefficients=None):
        r"""
        Return the brick polytope of ``self``.

        This polytope is the convex hull of the brick vectors of ``self``.

        INPUT:

        - ``coefficients`` -- (optional) a list of coefficients used to
          scale the fundamental weights

        .. SEEALSO::

            :meth:`brick_vectors`

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])                          # optional - gap3
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)                # optional - gap3
            sage: X = SC.brick_polytope(); X                            # optional - gap3
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 5 vertices

            sage: Y = SC.brick_polytope(coefficients=[1,2]); Y          # optional - gap3
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 5 vertices

            sage: X == Y                                                # optional - gap3
            False

            sage: W = CoxeterGroup(['A',2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)
            sage: X = SC.brick_polytope(); X
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 5 vertices

            sage: # optional - gap3
            sage: W = ReflectionGroup(['H',3])
            sage: c = W.index_set(); Q = c + tuple(W.w0.coxeter_sorting_word(c))
            sage: SC = SubwordComplex(Q,W.w0)
            sage: SC.brick_polytope()
            doctest:...:
            RuntimeWarning: the polytope is built with rational vertices
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 32 vertices
        """
        BV = self.brick_vectors(coefficients=coefficients)
        G = self.group()
        from sage.rings.rational_field import QQ
        if G.coxeter_matrix().is_crystallographic():
            BV = [[QQ(v) for v in V] for V in BV]
        else:
            from sage.rings.cc import CC
            from warnings import warn
            warn("the polytope is built with rational vertices", RuntimeWarning)
            BV = [[QQ(CC(v).real()) for v in V] for V in BV]
        return Polyhedron(BV)

    def barycenter(self):
        """
        Return the barycenter of the brick polytope of ``self``.

        .. SEEALSO::

            :meth:`brick_polytope`

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])                          # optional - gap3
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)                # optional - gap3
            sage: SC.barycenter()                                       # optional - gap3
            (2/3, 4/3)

            sage: W = CoxeterGroup(['A',2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)
            sage: SC.barycenter()
            (4/3, 8/3)
        """
        facets = self.facets()
        if not self.is_root_independent():
            facets = [F for F in facets if F.is_vertex()]
        return sum(F.brick_vector() for F in facets) / len(facets)

    # cambrian constructions

    def cover_relations(self, label=False):
        """
        Return the set of cover relations in the associated poset.

        INPUT:

        - ``label`` -- boolean (default: ``False``); whether or not to label
          the cover relations by the position of flip

        OUTPUT: list of pairs of facets

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])                          # optional - gap3
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)                # optional - gap3
            sage: sorted(SC.cover_relations())                          # optional - gap3
            [((0, 1), (0, 4)),
             ((0, 1), (1, 2)),
             ((0, 4), (3, 4)),
             ((1, 2), (2, 3)),
             ((2, 3), (3, 4))]

            sage: W = CoxeterGroup(['A',2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)
            sage: sorted(SC.cover_relations())
            [((0, 1), (0, 4)),
             ((0, 1), (1, 2)),
             ((0, 4), (3, 4)),
             ((1, 2), (2, 3)),
             ((2, 3), (3, 4))]
        """
        N = len(self.group().long_element(as_word=True))
        F = self.greedy_facet(side='positive')
        Fs = {F}
        seen = {F}
        covers = []
        while Fs:
            F = Fs.pop()
            seen.add(F)
            conf = F._extended_root_configuration_indices()
            for i in F:
                if conf[i] < N:
                    G = F.flip(i)
                    if label:
                        covers.append((F, G, i))
                    else:
                        covers.append((F, G))
                    if G not in seen:
                        Fs.add(G)
        return covers

    def increasing_flip_graph(self, label=True):
        """
        Return the increasing flip graph of the subword complex.

        OUTPUT: a directed graph

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])                          # optional - gap3
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)                # optional - gap3
            sage: SC.increasing_flip_graph()                            # optional - gap3
            Digraph on 5 vertices

            sage: W = CoxeterGroup(['A',2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)
            sage: SC.increasing_flip_graph()
            Digraph on 5 vertices
        """
        from sage.graphs.digraph import DiGraph
        return DiGraph(self.cover_relations(label=label))

    def interval(self, I, J) -> set:
        """
        Return the interval [I,J] in the increasing flip graph subword complex.

        INPUT:

        - I, J -- two facets

        OUTPUT: a set of facets

        EXAMPLES::

            sage: # optional - gap3
            sage: W = ReflectionGroup(['A',2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)
            sage: F = SC([1,2])
            sage: SC.interval(F, F)
            {(1, 2)}

            sage: W = CoxeterGroup(['A',2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)
            sage: F = SC([1,2])
            sage: SC.interval(F, F)
            {(1, 2)}
        """
        G = self.increasing_flip_graph()
        paths = G.all_paths(I, J)
        return {K for path in paths for K in path}

    def increasing_flip_poset(self):
        """
        Return the increasing flip poset of the subword complex.

        OUTPUT: a poset

        EXAMPLES::

            sage: W = ReflectionGroup(['A',2])                          # optional - gap3
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)                # optional - gap3
            sage: SC.increasing_flip_poset()                            # optional - gap3
            Finite poset containing 5 elements

            sage: W = CoxeterGroup(['A',2])
            sage: SC = SubwordComplex([1,2,1,2,1], W.w0)
            sage: SC.increasing_flip_poset()
            Finite poset containing 5 elements
        """
        from sage.combinat.posets.posets import Poset
        cov = self.cover_relations()
        if not self.is_root_independent():
            Fs = [F for F in self if F.is_vertex()]
            cov = [(a, b) for a, b in cov if a in Fs and b in Fs]
        return Poset(((), cov), facade=True)


def _greedy_facet(Q, w, side='negative', n=None, pos=0, l=None, elems=[]):
    r"""
    Return the (positive or negative) *greedy facet* of the subword
    complex `SC(Q, w)`.

    INPUT:

    - ``Q`` -- a word
    - ``w`` -- an element in the Coxeter group
    - ``side`` -- string; either ``'negative'`` (default) or ``'positive'``
    - ``n`` -- integer (default: the length of `Q`)
    - ``pos`` -- integer (default: 0)
    - ``l`` -- integer (default: the length of `w`)
    - ``elems`` -- list (optional)

    OUTPUT: a set

    EXAMPLES::

        sage: from sage.combinat.subword_complex import _greedy_facet

        sage: # optional - gap3
        sage: W = ReflectionGroup(['A',2])
        sage: Q = [1,2,1,2,1]
        sage: w = W.from_reduced_word([1, 2, 1])
        sage: _greedy_facet(Q, w)
        {3, 4}

        sage: W = CoxeterGroup(['A',2])
        sage: Q = [1,2,1,2,1]
        sage: w = W.from_reduced_word([1, 2, 1])
        sage: _greedy_facet(Q, w)
        {3, 4}
    """
    if side == "negative":
        pass
    elif side == "positive":
        Q = Q[::-1]
        w = w.inverse()
    else:
        raise ValueError("the optional argument side is not positive "
                         "or negative")

    if n is None:
        n = len(Q)
    if l is None:
        l = w.length()

    if l == 0:
        return elems + list(range(pos, n))
    elif n < l:
        return []

    s = Q[pos]

    if w.has_left_descent(s):
        X = _greedy_facet(Q, w.apply_simple_reflection_left(s),
                          n=n, pos=pos + 1, l=l - 1, elems=elems)
    else:
        X = []

    if not X:
        X = _greedy_facet(Q, w, n=n, pos=pos + 1, l=l, elems=elems + [pos])

    if side == "positive":
        X = [n - 1 - i for i in X]

    return set(X)


def _extended_root_configuration_indices(W, Q, F):
    """
    Return the extended root configuration indices of the facet `F`.

    INPUT:

    - ``W`` -- a Coxeter group
    - ``Q`` -- a word representing an element of `W`
    - ``F`` -- a facet of the subword complex

    OUTPUT: list of root indices

    EXAMPLES::

        sage: from sage.combinat.subword_complex import _extended_root_configuration_indices

        sage: # optional - gap3
        sage: W = ReflectionGroup(['A',2])
        sage: w = W.from_reduced_word([1,2,1])
        sage: Q = [1,2,1,2,1]
        sage: SC = SubwordComplex(Q, w)
        sage: F = SC([1,2])
        sage: _extended_root_configuration_indices(W, Q, F)
        [0, 2, 3, 2, 1]

        sage: W = CoxeterGroup(['A',2])
        sage: w = W.from_reduced_word([1,2,1])
        sage: Q = [1,2,1,2,1]
        sage: SC = SubwordComplex(Q, w)
        sage: F = SC([1,2])
        sage: _extended_root_configuration_indices(W, Q, F)
        [0, 1, 3, 1, 2]
    """
    V_roots = []
    pi = W.one()
    for i, wi in enumerate(Q):
        V_roots.append(pi.action_on_root_indices(W.simple_root_index(wi),
                                                 side='left'))
        if i not in F:
            pi = pi.apply_simple_reflection_right(wi)
    return V_roots


def _greedy_flip_algorithm(Q, w):
    """
    INPUT:

    - ``Q`` -- a word in a Coxeter group `W`
    - ``w`` -- an element of `W`

    OUTPUT: a pair: the list of facets and the list of extended root conf. indices

    EXAMPLES::

        sage: from sage.combinat.subword_complex import _greedy_flip_algorithm

        sage: # optional - gap3
        sage: W = ReflectionGroup(['A',2])
        sage: Q = [1,2,1,2,1]
        sage: w = W.from_reduced_word([1,2,1])
        sage: _greedy_flip_algorithm(Q, w)
        ([{0, 1}, [1, 2], [2, 3], [3, 4], [0, 4]],
         [[0, 1, 0, 2, 1],
          [0, 2, 3, 2, 1],
          [0, 2, 1, 5, 1],
          [0, 2, 1, 3, 4],
          [0, 1, 2, 0, 4]])

        sage: W = CoxeterGroup(['A',2])
        sage: Q = [1,2,1,2,1]
        sage: w = W.from_reduced_word([1,2,1])
        sage: _greedy_flip_algorithm(Q, w)
        ([{0, 1}, [1, 2], [2, 3], [3, 4], [0, 4]],
         [[0, 2, 0, 1, 2],
          [0, 1, 3, 1, 2],
          [0, 1, 2, 4, 2],
          [0, 1, 2, 3, 5],
          [0, 2, 1, 0, 5]])
    """
    W = w.parent()
    F = _greedy_facet(Q, w, side='positive')
    R = _extended_root_configuration_indices(W, Q, F)
    facet_list = [F]
    extended_root_conf_indices_list = [R]
    flip_to_ancestors = [-1]
    next_index = 0
    while flip_to_ancestors:
        has_new_child = False
        for i in sorted(F):
            if (not has_new_child) and (i >= next_index):
                j = _flip_c(W, F, R, i, side='positive')
                if j != i:
                    flip_to_ancestors.append(j)
                    next_index = i + 1
                    has_new_child = True
                    facet_list.append(list(F))
                    extended_root_conf_indices_list.append(list(R))
        if not has_new_child:
            i = flip_to_ancestors.pop()
            if i != -1:
                j = _flip_c(W, F, R, i, side='negative')
                next_index = j + 1
    return facet_list, extended_root_conf_indices_list
