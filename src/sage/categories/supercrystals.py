# sage_setup: distribution = sagemath-categories
# sage.doctest: needs sage.graphs sage.combinat
r"""
Supercrystals
"""

#*****************************************************************************
#       Copyright (C) 2017 Franco Saliola <saliola@gmail.com>
#                     2017 Anne Schilling <anne at math.ucdavis.edu>
#                     2019 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.crystals import Crystals
from sage.categories.tensor import TensorProductsCategory


class SuperCrystals(Category_singleton):
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: from sage.categories.supercrystals import SuperCrystals
            sage: C = SuperCrystals()
            sage: C.super_categories()
            [Category of crystals]
        """
        return [Crystals()]

    class ParentMethods:
        def tensor(self, *crystals, **options):
            """
            Return the tensor product of ``self`` with the crystals ``B``.

            EXAMPLES::

                sage: B = crystals.Letters(['A',[1,2]])
                sage: C = crystals.Tableaux(['A',[1,2]], shape = [2,1])
                sage: T = C.tensor(B); T
                Full tensor product of the crystals
                 [Crystal of BKK tableaux of shape [2, 1] of gl(2|3),
                  The crystal of letters for type ['A', [1, 2]]]
                sage: S = B.tensor(C); S
                Full tensor product of the crystals
                 [The crystal of letters for type ['A', [1, 2]],
                  Crystal of BKK tableaux of shape [2, 1] of gl(2|3)]
                sage: G = T.digraph()
                sage: H = S.digraph()
                sage: G.is_isomorphic(H, edge_labels=True)
                True
            """
            cartan_type = self.cartan_type()
            if any(c.cartan_type() != cartan_type for c in crystals):
                raise ValueError("all crystals must be of the same Cartan type")

            if cartan_type.letter == 'Q':
                from sage.combinat.crystals.tensor_product import FullTensorProductOfQueerSuperCrystals
                return FullTensorProductOfQueerSuperCrystals((self,) + tuple(crystals), **options)
            else:
                from sage.combinat.crystals.tensor_product import FullTensorProductOfSuperCrystals
                return FullTensorProductOfSuperCrystals((self,) + tuple(crystals), **options)

    class Finite(CategoryWithAxiom):
        class ParentMethods:
            @cached_method(key=lambda s,i: tuple(i) if i is not None else s.index_set())
            def digraph(self, index_set=None):
                r"""
                Return the :class:`DiGraph` associated to ``self``.

                EXAMPLES::

                    sage: B = crystals.Letters(['A', [1,3]])
                    sage: G = B.digraph(); G
                    Multi-digraph on 6 vertices
                    sage: Q = crystals.Letters(['Q',3])
                    sage: G = Q.digraph(); G
                    Multi-digraph on 3 vertices
                    sage: G.edges(sort=True)
                    [(1, 2, -1), (1, 2, 1), (2, 3, -2), (2, 3, 2)]

                The edges of the crystal graph are by default colored using
                blue for edge 1, red for edge 2, green for edge 3, and dashed with
                the corresponding color for barred edges. Edge 0 is dotted black::

                    sage: view(G)  # optional - dot2tex graphviz, not tested (opens external window)
                """
                from sage.graphs.digraph import DiGraph
                from sage.misc.latex import LatexExpr
                from sage.combinat.root_system.cartan_type import CartanType

                if index_set is None:
                    index_set = self.index_set()

                G = DiGraph(multiedges=True)
                G.add_vertices(self)
                for i in index_set:
                    for x in G:
                        y = x.f(i)
                        if y is not None:
                            G.add_edge(x, y, i)

                def edge_options(data):
                    u, v, l = data
                    edge_opts = { 'edge_string': '->', 'color': 'black' }
                    if l > 0:
                        edge_opts['color'] = CartanType._colors.get(l, 'black')
                        edge_opts['label'] = LatexExpr(str(l))
                    elif l < 0:
                        edge_opts['color'] = "dashed," + CartanType._colors.get(-l, 'black')
                        edge_opts['label'] = LatexExpr("\\overline{%s}" % str(-l))
                    else:
                        edge_opts['color'] = "dotted," + CartanType._colors.get(l, 'black')
                        edge_opts['label'] = LatexExpr(str(l))
                    return edge_opts

                G.set_latex_options(format='dot2tex', edge_labels=True, edge_options=edge_options)
                return G

            def genuine_highest_weight_vectors(self):
                r"""
                Return the tuple of genuine highest weight elements of ``self``.

                EXAMPLES::

                    sage: B = crystals.Letters(['A', [1,2]])
                    sage: B.genuine_highest_weight_vectors()
                    (-2,)

                    sage: T = B.tensor(B)
                    sage: T.genuine_highest_weight_vectors()
                    ([-2, -1], [-2, -2])
                    sage: s1, s2 = T.connected_components()
                    sage: s = s1 + s2
                    sage: s.genuine_highest_weight_vectors()
                    ([-2, -1], [-2, -2])
                """
                return tuple([x[0] for x in self._genuine_highest_lowest_weight_vectors()])

            connected_components_generators = genuine_highest_weight_vectors

            def connected_components(self):
                r"""
                Return the connected components of ``self`` as subcrystals.

                EXAMPLES::

                    sage: B = crystals.Letters(['A', [1,2]])
                    sage: B.connected_components()
                    [Subcrystal of The crystal of letters for type ['A', [1, 2]]]

                    sage: T = B.tensor(B)
                    sage: T.connected_components()
                    [Subcrystal of Full tensor product of the crystals
                      [The crystal of letters for type ['A', [1, 2]],
                       The crystal of letters for type ['A', [1, 2]]],
                     Subcrystal of Full tensor product of the crystals
                      [The crystal of letters for type ['A', [1, 2]],
                       The crystal of letters for type ['A', [1, 2]]]]
                """
                category = SuperCrystals()
                from sage.categories.regular_supercrystals import RegularSuperCrystals
                if self in RegularSuperCrystals():
                    category = RegularSuperCrystals()
                index_set = self.index_set()
                cartan_type = self.cartan_type()
                CCs = []

                for mg in self.connected_components_generators():
                    if not isinstance(mg, tuple):
                        mg = (mg,)
                    subcrystal = self.subcrystal(generators=mg,
                                                 index_set=index_set,
                                                 cartan_type=cartan_type,
                                                 category=category)
                    CCs.append(subcrystal)

                return CCs

            def genuine_lowest_weight_vectors(self):
                r"""
                Return the tuple of genuine lowest weight elements of ``self``.

                EXAMPLES::

                    sage: B = crystals.Letters(['A', [1,2]])
                    sage: B.genuine_lowest_weight_vectors()
                    (3,)

                    sage: T = B.tensor(B)
                    sage: T.genuine_lowest_weight_vectors()
                    ([3, 3], [3, 2])
                    sage: s1, s2 = T.connected_components()
                    sage: s = s1 + s2
                    sage: s.genuine_lowest_weight_vectors()
                    ([3, 3], [3, 2])
                """
                return tuple([x[1] for x in self._genuine_highest_lowest_weight_vectors()])

            @cached_method
            def _genuine_highest_lowest_weight_vectors(self):
                r"""
                Return the genuine lowest and highest weight elements of ``self``.

                EXAMPLES::

                    sage: B = crystals.Letters(['A', [1,2]])
                    sage: B._genuine_highest_lowest_weight_vectors()
                    ((-2, 3),)

                    sage: T = B.tensor(B)
                    sage: T._genuine_highest_lowest_weight_vectors()
                    (([-2, -1], [3, 3]), ([-2, -2], [3, 2]))
                    sage: s1, s2 = T.connected_components()
                    sage: s = s1 + s2
                    sage: s._genuine_highest_lowest_weight_vectors()
                    (([-2, -1], [3, 3]), ([-2, -2], [3, 2]))

                An example with fake highest/lowest weight elements
                from [BKK2000]_::

                    sage: B = crystals.Tableaux(['A', [1,1]], shape=[3,2,1])
                    sage: B._genuine_highest_lowest_weight_vectors()
                    (([[-2, -2, -2], [-1, -1], [1]], [[-1, 1, 2], [1, 2], [2]]),)
                """
                X = []
                for G in self.digraph().connected_components_subgraphs():
                    src = G.sources()
                    sinks = G.sinks()
                    max_dist = -1
                    pair = None
                    for s in src:
                        for t in sinks:
                            d = G.distance(s, t)
                            if d < float('inf') and d > max_dist:
                                pair = (s, t)
                                max_dist = d
                    X.append(pair)
                return tuple(X)

            def character(self):
                """
                Return the character of ``self``.

                .. TODO::

                    Once the `WeylCharacterRing` is implemented, make this
                    consistent with the implementation in
                    :meth:`sage.categories.classical_crystals.ClassicalCrystals.ParentMethods.character`.

                EXAMPLES::

                    sage: B = crystals.Letters(['A',[1,2]])
                    sage: B.character()
                    B[(1, 0, 0, 0, 0)] + B[(0, 1, 0, 0, 0)] + B[(0, 0, 1, 0, 0)]
                     + B[(0, 0, 0, 1, 0)] + B[(0, 0, 0, 0, 1)]
                """
                from sage.rings.integer_ring import ZZ
                A = self.weight_lattice_realization().algebra(ZZ)
                return A.sum(A(x.weight()) for x in self)

            @cached_method
            def highest_weight_vectors(self):
                """
                Return the highest weight vectors of ``self``.

                EXAMPLES::

                    sage: B = crystals.Letters(['A', [1,2]])
                    sage: B.highest_weight_vectors()
                    (-2,)

                    sage: T = B.tensor(B)
                    sage: T.highest_weight_vectors()
                    ([-2, -2], [-2, -1])

                We give an example from [BKK2000]_ that has fake
                highest weight vectors::

                    sage: B = crystals.Tableaux(['A', [1,1]], shape=[3,2,1])
                    sage: B.highest_weight_vectors()
                    ([[-2, -2, -2], [-1, -1], [1]],
                     [[-2, -2, -2], [-1, 2], [1]],
                     [[-2, -2, 2], [-1, -1], [1]])
                    sage: B.genuine_highest_weight_vectors()
                    ([[-2, -2, -2], [-1, -1], [1]],)
                """
                return tuple(self.digraph().sources())

            @cached_method
            def lowest_weight_vectors(self):
                """
                Return the lowest weight vectors of ``self``.

                EXAMPLES::

                    sage: B = crystals.Letters(['A', [1,2]])
                    sage: B.lowest_weight_vectors()
                    (3,)

                    sage: T = B.tensor(B)
                    sage: sorted(T.lowest_weight_vectors())
                    [[3, 2], [3, 3]]

                We give an example from [BKK2000]_ that has fake
                lowest weight vectors::

                    sage: B = crystals.Tableaux(['A', [1,1]], shape=[3,2,1])
                    sage: sorted(B.lowest_weight_vectors())
                    [[[-2, 1, 2], [-1, 2], [1]],
                     [[-2, 1, 2], [-1, 2], [2]],
                     [[-1, 1, 2], [1, 2], [2]]]
                    sage: B.genuine_lowest_weight_vectors()
                    ([[-1, 1, 2], [1, 2], [2]],)
                """
                return tuple(self.digraph().sinks())

        class ElementMethods:
            def is_genuine_highest_weight(self, index_set=None):
                """
                Return whether ``self`` is a genuine highest weight element.

                INPUT:

                - ``index_set`` -- (optional) the index set of the (sub)crystal
                  on which to check

                EXAMPLES::

                    sage: B = crystals.Tableaux(['A', [1,1]], shape=[3,2,1])
                    sage: for b in B.highest_weight_vectors():
                    ....:     print("{} {}".format(b, b.is_genuine_highest_weight()))
                    [[-2, -2, -2], [-1, -1], [1]] True
                    [[-2, -2, -2], [-1, 2], [1]] False
                    [[-2, -2, 2], [-1, -1], [1]] False
                    sage: [b for b in B if b.is_genuine_highest_weight([-1,0])]
                    [[[-2, -2, -2], [-1, -1], [1]],
                     [[-2, -2, -2], [-1, -1], [2]],
                     [[-2, -2, -2], [-1, 2], [2]],
                     [[-2, -2, 2], [-1, -1], [2]],
                     [[-2, -2, 2], [-1, 2], [2]],
                     [[-2, -2, -2], [-1, 2], [1]],
                     [[-2, -2, 2], [-1, -1], [1]],
                     [[-2, -2, 2], [-1, 2], [1]]]
                """
                P = self.parent()
                if index_set is None or set(index_set) == set(P.index_set()):
                    return self in P.genuine_highest_weight_vectors()
                S = P.subcrystal(generators=P, index_set=index_set, category=P.category())
                return any(self == x.value for x in S.genuine_highest_weight_vectors())

            def is_genuine_lowest_weight(self, index_set=None):
                """
                Return whether ``self`` is a genuine lowest weight element.

                INPUT:

                - ``index_set`` -- (optional) the index set of the (sub)crystal
                  on which to check

                EXAMPLES::

                    sage: B = crystals.Tableaux(['A', [1,1]], shape=[3,2,1])
                    sage: for b in sorted(B.lowest_weight_vectors()):
                    ....:     print("{} {}".format(b, b.is_genuine_lowest_weight()))
                    [[-2, 1, 2], [-1, 2], [1]] False
                    [[-2, 1, 2], [-1, 2], [2]] False
                    [[-1, 1, 2], [1, 2], [2]] True
                    sage: [b for b in B if b.is_genuine_lowest_weight([-1,0])]
                    [[[-2, -1, 1], [-1, 1], [1]],
                     [[-2, -1, 1], [-1, 1], [2]],
                     [[-2, 1, 2], [-1, 1], [2]],
                     [[-2, 1, 2], [-1, 1], [1]],
                     [[-1, -1, 1], [1, 2], [2]],
                     [[-1, -1, 1], [1, 2], [1]],
                     [[-1, 1, 2], [1, 2], [2]],
                     [[-1, 1, 2], [1, 2], [1]]]
                """
                P = self.parent()
                if index_set is None or set(index_set) == set(P.index_set()):
                    return self in P.genuine_lowest_weight_vectors()
                S = P.subcrystal(generators=P, index_set=index_set, category=P.category())
                return any(self == x.value for x in S.genuine_lowest_weight_vectors())

    class TensorProducts(TensorProductsCategory):
        """
        The category of regular crystals constructed by tensor
        product of regular crystals.
        """
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: from sage.categories.supercrystals import SuperCrystals
                sage: SuperCrystals().TensorProducts().extra_super_categories()
                [Category of super crystals]
            """
            return [self.base_category()]
