r"""
Auslander-Reiten Quivers
"""

# ****************************************************************************
#  Copyright (C) 2024 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.richcmp import richcmp
from sage.structure.global_options import GlobalOptions
from sage.categories.sets_cat import Sets
from sage.sets.family import Family
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.free_module import CombinatorialFreeModule
from sage.rings.integer_ring import ZZ
from sage.misc.cachefunc import cached_method
from sage.graphs.digraph import DiGraph


class AuslanderReitenQuiver(UniqueRepresentation, Parent):
    r"""
    The Auslander-Reiten quiver.

    Let `Q = (Q_0, Q_1)` be a finite acyclic quiver. The
    *Auslander-Reiten quiver* (AR quiver) `\Gamma_Q` is the quiver
    whose vertices correspond to the indecompositible modules of `Q`
    (equivalently its path algebra over an algebraically closed field)
    and edges are irreducible morphisms.

    In this implementation, we denote the vertices of `\Gamma_Q` as
    certain pairs `\langle v, k \rangle`, where `v \in Q_0` and
    `k \in \ZZ \setminus \{0\}` is called the *level*. When `k > 0`
    (resp. `k < 0`), then it corresponds to a preprojective (resp.
    postinjective) module. When the quiver is a finite type Dynkin
    quiver, we consider all modules to be preprojectives and denoted
    by a positive level.

    .. NOTE::

        We use the terminology *postinjective* instead of *preinjective*
        given that they follow from injectives by AR translation.

    ALGORITHM:

    We compute the dimension vectors of a projective `\langle v, 1 \rangle`
    by counting the number of (directed) paths `u \to v` in `Q`. We then
    proceed inductively to compute all of the dimension vectors of level
    `k` by using the translation equation

    .. MATH::

        dim \langle v, k-1 \rangle + \dim \langle v, k \rangle
        = \sum_{u,k'} \dim \langle u, k' \rangle,

    where the sum is over all paths from `\langle v, k-1 \rangle` to
    `\langle v, k \rangle` in `\Gamma_Q`. More specifically, for each edge
    `(u, v, \ell) \in Q_1` (resp. `(v, u, \ell) \in Q_1`), we have
    `\langle u, k-1 \rangle` (resp. `\langle u, k \rangle`) in the sum
    (assuming the node is in the AR quiver).

    The algorithm for postinjectives is dual to the above.

    .. TODO::

        This only is implemented for the preprojectives and postinjectives
        when the quiver is not a finite type Dynkin quiver.

    .. TODO::

        Implement this for general Artinian algebras.

    EXAMPLES:

    We create the AR quivers for finite type `A_3` Dynkin quivers::

        sage: DA = DiGraph([[1, 2], [2, 3]])
        sage: AR = DA.auslander_reiten_quiver()
        sage: AR.digraph().edges(labels=False)
        [(<1, 1>, <2, 2>), (<2, 1>, <1, 1>), (<2, 1>, <3, 2>), (<3, 1>, <2, 1>),
         (<2, 2>, <3, 3>), (<3, 2>, <2, 2>)]

        sage: DA = DiGraph([[1, 2], [3, 2]])
        sage: AR = DA.auslander_reiten_quiver()
        sage: AR.digraph().edges(labels=False)
        [(<1, 1>, <2, 2>), (<2, 1>, <1, 1>), (<2, 1>, <3, 1>), (<3, 1>, <2, 2>),
         (<2, 2>, <1, 2>), (<2, 2>, <3, 2>)]

        sage: DA = DiGraph([[2, 1], [2, 3]])
        sage: AR = DA.auslander_reiten_quiver()
        sage: AR.digraph().edges(labels=False)
        [(<1, 1>, <2, 1>), (<2, 1>, <1, 2>), (<2, 1>, <3, 2>), (<3, 1>, <2, 1>),
         (<1, 2>, <2, 2>), (<3, 2>, <2, 2>)]

        sage: DA = DiGraph([[2, 1], [3, 2]])
        sage: AR = DA.auslander_reiten_quiver()
        sage: AR.digraph().edges(labels=False)
        [(<1, 1>, <2, 1>), (<2, 1>, <3, 1>), (<2, 1>, <1, 2>), (<3, 1>, <2, 2>),
         (<1, 2>, <2, 2>), (<2, 2>, <1, 3>)]

    An example for the type `D_5` Dynkin quiver::

        sage: DD = DiGraph([[5,3], [4,3], [3,2], [2,1]])
        sage: AR = DD.auslander_reiten_quiver()
        sage: AR
        Auslander-Reiten quiver of a ['D', 5] Dynkin quiver
        sage: len(list(DD))
        5

    An `E_8` Dynkin quiver::

        sage: DE = DiGraph([[8,7], [7,6], [5,6], [5,3], [3,4], [3,2], [2,1]])
        sage: AR = DE.auslander_reiten_quiver()
        sage: AR
        Auslander-Reiten quiver of a ['E', 8] Dynkin quiver
        sage: len(list(AR))
        120
        sage: len(list(RootSystem(['E', 8]).root_lattice().positive_roots()))
        120

    The Kronecker quiver::

        sage: D = DiGraph([[1,2,'a'], [1,2,'b']], multiedges=True)
        sage: AR = D.auslander_reiten_quiver()
        sage: for i in range(1, 5):
        ....:     for v in D.vertices():
        ....:         pp = AR(v, i)
        ....:         pi = AR(v, -i)
        ....:         print(pp, pp.dimension_vector(), "  ", pi, pi.dimension_vector())
        <1, 1> v1 + 2*v2      <1, -1> v1
        <2, 1> v2             <2, -1> 2*v1 + v2
        <1, 2> 3*v1 + 4*v2    <1, -2> 3*v1 + 2*v2
        <2, 2> 2*v1 + 3*v2    <2, -2> 4*v1 + 3*v2
        <1, 3> 5*v1 + 6*v2    <1, -3> 5*v1 + 4*v2
        <2, 3> 4*v1 + 5*v2    <2, -3> 6*v1 + 5*v2
        <1, 4> 7*v1 + 8*v2    <1, -4> 7*v1 + 6*v2
        <2, 4> 6*v1 + 7*v2    <2, -4> 8*v1 + 7*v2
    """
    @staticmethod
    def __classcall_private__(cls, quiver):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: D = DiGraph([[1,2], [2,3], [3,1]])
            sage: D.auslander_reiten_quiver()
            Traceback (most recent call last):
            ...
            ValueError: the quiver must not have cycles
        """
        if quiver.has_loops() or not quiver.is_directed_acyclic():
            raise ValueError("the quiver must not have cycles")
        quiver = quiver.copy(immutable=True)
        return super().__classcall__(cls, quiver)

    def __init__(self, quiver):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: DE = DiGraph([[7,6], [6,5], [5,3], [3,4], [2,3], [1,2]])
            sage: AR = DE.auslander_reiten_quiver()
            sage: TestSuite(AR).run()

            sage: D = DiGraph([[1,2], [3,4]])
            sage: AR = D.auslander_reiten_quiver()
            sage: TestSuite(AR).run()

            sage: D = DiGraph([[1,2,'a'], [1,2,'b']], multiedges=True)
            sage: AR = D.auslander_reiten_quiver()
            sage: TestSuite(AR).run()
        """
        self._quiver = quiver
        self._top_sort = quiver.topological_sort()
        self._dim_vec_space = CombinatorialFreeModule(ZZ, quiver.vertices(), prefix='v', bracket=False)
        self._max_level = float('inf')

        dynkin_type = detect_dynkin_quiver(quiver)
        if dynkin_type is not None:
            self._cartan_type = dynkin_type
        self._is_finite = dynkin_type is not None
        cat = Sets().Enumerated().Finite() if self._is_finite else Sets().Infinite()
        super().__init__(self, category=cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: DE = DiGraph([[5,6], [5,3], [3,4], [2,3], [1,2]])
            sage: DE.auslander_reiten_quiver()
            Auslander-Reiten quiver of a ['E', 6] Dynkin quiver

            sage: D = DiGraph([[1,2,'a'], [1,2,'b']], multiedges=True)
            sage: D.auslander_reiten_quiver()
            Auslander-Reiten quiver of Multi-digraph on 2 vertices
        """
        if self._is_finite:
            return "Auslander-Reiten quiver of a {} Dynkin quiver".format(self._cartan_type)
        return "Auslander-Reiten quiver of {}".format(self._quiver)

    # add options to class
    class options(GlobalOptions):
        r"""
        Set and display the global options for Auslander-Reiten quivers.
        If no parameters are set, then the function returns a copy of the
        options dictionary.

        The ``options`` to partitions can be accessed as the method
        :obj:`AuslanderReitenQuiver.options` of
        :class:`~sage.quivers.ar_quiver.AuslanderReitenQuiver`.

        @OPTIONS@

        EXAMPLES::

            sage: D = DiGraph([[1,2,'a'], [1,2,'b']], multiedges=True)
            sage: AR = D.auslander_reiten_quiver()
            sage: node = AR(2, 2)
            sage: latex(node)
            \left\langle 2, 2 \right\rangle
            sage: AR.options.latex = "dimension_vector"
            sage: latex(node)
            2 v_{1} + 3 v_{2}
            sage: AR.options.latex = "both"
            sage: latex(node)
            \begin{gathered} \left\langle 2, 2 \right\rangle \\ 2 v_{1} + 3 v_{2} \end{gathered}
            sage: AR.options._reset()
        """
        NAME = 'AuslanderReitenQuiver'
        module = 'sage.quivers.ar_quiver'
        latex = dict(default='node',
                     description='Specifies how nodes of the AR quiver should be latexed',
                     values=dict(node='latex as the node description',
                                 dimension_vector='latex as the dimension vector',
                                 both='latex as both'),
                     case_sensitive=False)

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: DE = DiGraph([[7,8], [7,6], [5,6], [3,5], [4,3], [2,3], [1,2]])
            sage: AR = DE.auslander_reiten_quiver()
            sage: AR._an_element_()
            <1, 1>

            sage: D = DiGraph([[1,2,'a'], [1,2,'b']], multiedges=True)
            sage: AR = D.auslander_reiten_quiver()
            sage: AR._an_element_()
            <1, 1>
        """
        return next(iter(self.projectives()))

    def quiver(self):
        r"""
        Return the quiver defining ``self``.

        EXAMPLES::

            sage: DE = DiGraph([[7,8], [7,6], [5,6], [3,5], [4,3], [2,3], [1,2]])
            sage: AR = DE.auslander_reiten_quiver()
            sage: AR.quiver() == DE
            True
        """
        return self._quiver

    def projectives(self):
        r"""
        Return the projectives of ``self``.

        EXAMPLES::

            sage: D = DiGraph([[1,2,'a'], [1,2,'b']], multiedges=True)
            sage: AR = D.auslander_reiten_quiver()
            sage: AR.projectives()
            Finite family {1: <1, 1>, 2: <2, 1>}
        """
        return Family({v: self.element_class(self, v, 1) for v in self._quiver.vertex_iterator()})

    @cached_method
    def simples(self):
        r"""
        Return the simples of ``self``.

        EXAMPLES::

            sage: DE = DiGraph([[7,8], [7,6], [5,6], [3,5], [4,3], [2,3], [1,2]])
            sage: AR = DE.auslander_reiten_quiver()
            sage: AR.simples()
            Finite family {1: <1, 15>,  2: <1, 14>,  3: <8, 4>,  4: <4, 15>,
                           5: <8, 3>,  6: <6, 1>,  7: <7, 15>,  8: <8, 1>}
        """
        ret = {}
        for elt in self:
            supp = elt.dimension_vector().support()
            if len(supp) != 1:
                continue
            ret[next(iter(supp))] = elt
        return Family(ret)

    def injectives(self):
        r"""
        Return the injectives of ``self``.

        EXAMPLES::

            sage: DE = DiGraph([[7,6], [6,5], [5,3], [4,3], [2,3], [1,2]])
            sage: AR = DE.auslander_reiten_quiver()
            sage: AR.injectives()
            Finite family {1: <1, 9>, 2: <2, 9>, 3: <3, 9>, 4: <4, 9>,
                           5: <5, 9>, 6: <6, 9>, 7: <7, 9>}
        """
        if self._is_finite:
            self.digraph()  # sets self._injective attribute
            return self._injectives
        return Family({v: self(v, -1) for v in self._quiver.vertex_iterator()})

    def _digraph_set_latex_options(self, G):
        """
        Set the latex options of the digraph ``G``.

        EXAMPLES::

            sage: D = DiGraph([[1,2,'a'], [1,2,'b']], multiedges=True)
            sage: AR = D.auslander_reiten_quiver()
            sage: G = AR.digraph_preprojectives(2)
            sage: G = AR._digraph_set_latex_options(G)
            sage: G.latex_options().get_option('edge_labels')
            True
        """
        G.set_latex_options(edge_labels=True)

        from sage.graphs.dot2tex_utils import have_dot2tex
        if have_dot2tex():
            from sage.misc.latex import LatexExpr

            def edge_options(data):
                u, v, l = data
                edge_opts = {}
                if l == 'ART':
                    edge_opts["color"] = "dashed,blue"
                    edge_opts["label"] = LatexExpr(r"\tau")
                return edge_opts

            G.set_latex_options(format='dot2tex', edge_options=edge_options)
        return G

    def digraph_preprojectives(self, max_depth, with_translations=False):
        r"""
        Return the diagraph of preprojectives of ``self`` up to ``max_depth``.

        EXAMPLES::

            sage: D = DiGraph([[1,2,'a'], [1,2,'b']], multiedges=True)
            sage: AR = D.auslander_reiten_quiver()
            sage: G = AR.digraph_preprojectives(3)
            sage: [node.dimension_vector() for node in G]
            [v1 + 2*v2, v2, 3*v1 + 4*v2, 2*v1 + 3*v2, 5*v1 + 6*v2, 4*v1 + 5*v2]
            sage: AR.digraph_preprojectives(0)
            Digraph on 0 vertices
        """
        if max_depth < 1:
            return self._digraph_set_latex_options(DiGraph())

        k = 2
        prev = dict(self.projectives())
        verts = list(prev.values())
        edges = [(prev[v], prev[u], l) for u, v, l in self._quiver.edge_iterator()]
        cur = self._dim_vecs_level(k)
        while k <= max_depth:
            # convert cur to the appropriate data
            cur = {v: self.element_class(self, v, k) for v in cur}
            verts.extend(cur.values())
            edges.extend((cur[v], cur[u], l)
                         for u in cur for _, v, l in self._quiver.outgoing_edge_iterator(u) if v in cur)
            edges.extend((prev[u], cur[v], l) for v in cur
                         for u, _, l in self._quiver.incoming_edge_iterator(v) if u in prev)
            if with_translations:
                edges.extend((cur[v], prev[v], 'ART') for v in cur if v in prev)
            k += 1
            prev = cur
            cur = self._dim_vecs_level(k)

        G = DiGraph([verts, edges], format='vertices_and_edges', multiedges=True, immutable=True)
        return self._digraph_set_latex_options(G)

    def digraph_postinjectives(self, max_depth, with_translations=False):
        """
        Return the diagraph of postinjectives of ``self`` up to ``max_depth``.

        EXAMPLES::

            sage: D = DiGraph([[1,2,'a'], [1,2,'b']], multiedges=True)
            sage: AR = D.auslander_reiten_quiver()
            sage: G = AR.digraph_postinjectives(3)
            sage: [node.dimension_vector() for node in G]
            [5*v1 + 4*v2, 6*v1 + 5*v2, 3*v1 + 2*v2, 4*v1 + 3*v2, v1, 2*v1 + v2]
            sage: AR.digraph_postinjectives(0)
            Digraph on 0 vertices
        """
        if max_depth < 1:
            return self._digraph_set_latex_options(DiGraph())

        k = 2
        prev = dict(self.injectives())
        verts = list(prev.values())
        edges = [(prev[u], prev[v], l) for u, v, l in self._quiver.edge_iterator()]
        cur = self._dim_vecs_level(-k)
        while k <= max_depth:
            # convert cur to the appropriate data
            cur = {v: self.element_class(self, v, -k) for v in cur}
            verts.extend(cur.values())
            edges.extend((cur[u], cur[v], l)
                         for u in cur for _, v, l in self._quiver.outgoing_edge_iterator(u) if v in cur)
            edges.extend((cur[v], prev[u], l) for v in cur
                         for u, _, l in self._quiver.incoming_edge_iterator(v) if u in prev)
            if with_translations:
                edges.extend((prev[v], cur[v], 'ART') for v in cur if v in prev)
            k += 1
            prev = cur
            cur = self._dim_vecs_level(-k)

        G = DiGraph([verts, edges], format='vertices_and_edges', multiedges=True, immutable=True)
        return self._digraph_set_latex_options(G)

    @cached_method
    def digraph(self, with_translations=False):
        r"""
        Return the diagraph of ``self``.

        INPUT:

        - ``with_translations`` -- boolean (default: ``False``); if ``True``, then
          include the arrows corresponding to the translations

        EXAMPLES::

            sage: DA = DiGraph([[1,2]])
            sage: AR = DA.auslander_reiten_quiver()
            sage: G = AR.digraph(); G
            Digraph on 3 vertices
            sage: G.edges()
            [(<1, 1>, <2, 2>, None), (<2, 1>, <1, 1>, None)]
            sage: GT = AR.digraph(with_translations=True)
            sage: GT.edges()
            [(<1, 1>, <2, 2>, None), (<2, 1>, <1, 1>, None), (<2, 2>, <2, 1>, 'ART')]
        """
        if not self._is_finite:
            raise TypeError("the AR quiver is not finite")

        if with_translations:
            G = self.digraph().copy(immutable=False)
            for v in G.vertex_iterator():
                u = v.translation()
                if u is not None:
                    G.add_edge(v, u, 'ART')
            G = G.copy(immutable=True)

        else:
            k = 2
            prev = dict(self.projectives())
            injectives = dict(prev)  # make a shallow copy since we will mutate it
            verts = list(prev.values())
            edges = [(prev[v], prev[u], l) for u, v, l in self._quiver.edge_iterator()]
            cur = self._dim_vecs_level(k)
            while cur:
                # convert cur to the appropriate data
                cur = {v: self.element_class(self, v, k) for v in cur}
                injectives.update(cur)
                verts.extend(cur.values())
                edges.extend((cur[v], cur[u], l)
                             for u in cur for _, v, l in self._quiver.outgoing_edge_iterator(u) if v in cur)
                edges.extend((prev[u], cur[v], l) for v in cur
                             for u, _, l in self._quiver.incoming_edge_iterator(v) if u in prev)
                k += 1
                prev = cur
                cur = self._dim_vecs_level(k)

            self._injectives = Family(injectives)
            G = DiGraph([verts, edges], format='vertices_and_edges', immutable=True)

        return self._digraph_set_latex_options(G)

    def __iter__(self):
        r"""
        Iterate over ``self`` when possible.

        EXAMPLES::

            sage: DD = DiGraph([[3,2], [4,2], [2,1]])
            sage: AR = DD.auslander_reiten_quiver()
            sage: list(AR)
            [<1, 1>, <2, 1>, <3, 1>, <4, 1>, <1, 2>, <2, 2>, <3, 2>, <4, 2>,
             <1, 3>, <2, 3>, <3, 3>, <4, 3>]
        """
        return iter(self.digraph())

    def _element_constructor_(self, vertex, level=None):
        r"""
        Construct an element of ``self``.

        EXAMPLES::

            sage: DA = DiGraph([[4,3], [3,2], [2,1]])
            sage: AR = DA.auslander_reiten_quiver()
            sage: AR(2, 2)
            <2, 2>
            sage: AR((2, 3))
            <2, 3>
            sage: AR(2, 4)
            Traceback (most recent call last):
            ...
            ValueError: no 2 at level 4
            sage: AR(2, 1/2)
            Traceback (most recent call last):
            ...
            ValueError: the level 1/2 must be an integer
            sage: AR(10, 1)
            Traceback (most recent call last):
            ...
            ValueError: 10 is not a vertex of the quiver
        """
        if level is None:
            if len(vertex) == 2:
                vertex, level = vertex
        if vertex not in self._quiver:
            raise ValueError(f"{vertex} is not a vertex of the quiver")
        if level == 1:
            return self.element_class(self, vertex, level)
        if level not in ZZ:
            raise ValueError(f"the level {level} must be an integer")

        if not self._is_finite:
            if level == -1 or vertex in self._dim_vecs_level(level):
                return self.element_class(self, vertex, level)
            # This is likely never true
            raise ValueError(f"no {vertex} at level {level}")

        # otherwise the AR quiver is finite
        if level < 0:
            self.digraph()  # computes the max level
            level = self._max_level - level
        if level > 1:
            if vertex in self._dim_vecs_level(level):
                return self.element_class(self, vertex, level)
            raise ValueError(f"no {vertex} at level {level}")

    @cached_method
    def _dim_vecs_level(self, k):
        r"""
        Return a ``dict`` of dimension vectors of level ``k``.

        .. WARNING::

            This is only meant to be used internally as the output is
            mutable but cached. Thus, the output should not be changed.

        EXAMPLES::

            sage: DA = DiGraph([[4,3], [3,2], [2,1]])
            sage: AR = DA.auslander_reiten_quiver()
            sage: AR._dim_vecs_level(1)
            {1: v1, 2: v1 + v2, 3: v1 + v2 + v3, 4: v1 + v2 + v3 + v4}
            sage: AR._dim_vecs_level(2)
            {1: v2, 2: v2 + v3, 3: v2 + v3 + v4}
            sage: AR._dim_vecs_level(3)
            {1: v3, 2: v3 + v4}
            sage: AR._dim_vecs_level(4)
            {1: v4}
            sage: AR._dim_vecs_level(-1)
            {1: v1 + v2 + v3 + v4, 2: v2 + v3 + v4, 3: v3 + v4, 4: v4}
            sage: AR._dim_vecs_level(-2)
            {2: v1 + v2 + v3, 3: v2 + v3, 4: v3}
            sage: AR._dim_vecs_level(-3)
            {3: v1 + v2, 4: v2}
            sage: AR._dim_vecs_level(-4)
            {4: v1}
        """
        if k == 0:
            raise ValueError("k must not be 0")
        M = self._dim_vec_space
        Q = self._quiver
        if k == 1:
            ret = {v: M._from_dict({u: ZZ(len(Q.all_paths(v, u, use_multiedges=True))) for u in Q.vertex_iterator()})
                   for v in Q.vertex_iterator()}
        elif k > 1:
            if k > self._max_level:
                return {}
            prev = self._dim_vecs_level(k-1)
            if k > self._max_level:  # this might get set on the recursive call
                return {}
            ret = {}
            for v in reversed(self._top_sort):
                if v not in prev:  # assumption: this vertex will never reappear
                    continue
                temp = -prev[v]
                for u, _, _ in Q.incoming_edge_iterator(v):
                    if u in prev:
                        temp += prev[u]
                for _, u, _ in Q.outgoing_edge_iterator(v):
                    if u in ret:
                        temp += ret[u]
                if all(coeff > 0 for key, coeff in temp):
                    ret[v] = temp
            if not ret:
                self._max_level = k

        elif k == -1:
            ret = {v: M._from_dict({u: ZZ(len(Q.all_paths(u, v, use_multiedges=True))) for u in Q.vertex_iterator()})
                   for v in Q.vertex_iterator()}

        elif k < -1:
            prev = self._dim_vecs_level(k+1)
            ret = {}
            for v in self._top_sort:
                if v not in prev:  # assumption: this vertex will never reappear
                    continue
                temp = -prev[v]
                for _, u, _ in Q.outgoing_edge_iterator(v):
                    if u in prev:
                        temp += prev[u]
                for u, _, _ in Q.incoming_edge_iterator(v):
                    if u in ret:
                        temp += ret[u]
                if all(coeff > 0 for key, coeff in temp):
                    ret[v] = temp

        return ret

    def dimension_vectors_of_level(self, k):
        r"""
        Return a :class:`Family` of dimension vectors of level ``k``.

        EXAMPLES::

            sage: DA = DiGraph([[4,3], [2,3], [2,1]])
            sage: AR = DA.auslander_reiten_quiver()
            sage: AR.dimension_vectors_of_level(1)
            {1: v1, 2: v1 + v2 + v3, 3: v3, 4: v3 + v4}
            sage: AR.dimension_vectors_of_level(3)
            {1: v4, 3: v2}
            sage: AR.dimension_vectors_of_level(10)
            {}
            sage: AR.dimension_vectors_of_level(-1)
             {1: v1 + v2, 2: v2, 3: v2 + v3 + v4, 4: v4}
            sage: AR.dimension_vectors_of_level(-2)
            {1: v3 + v4, 2: v1 + v2 + v3 + v4, 3: v1 + v2 + v3, 4: v2 + v3}

            sage: D = DiGraph([[1,2,'a'], [1,2,'b']], multiedges=True)
            sage: AR = D.auslander_reiten_quiver()
            sage: AR.dimension_vectors_of_level(1)
            {1: v1 + 2*v2, 2: v2}
            sage: AR.dimension_vectors_of_level(3)
            {1: 5*v1 + 6*v2, 2: 4*v1 + 5*v2}
            sage: AR.dimension_vectors_of_level(-1)
            {1: v1, 2: 2*v1 + v2}
            sage: AR.dimension_vectors_of_level(-3)
            {1: 5*v1 + 4*v2, 2: 6*v1 + 5*v2}
        """
        ret = self._dim_vecs_level(k)
        return dict(ret)  # make a (shallow) copy to allow a user to mutate it

    class Element(Element):
        r"""
        A node in the AR quiver.
        """
        def __init__(self, parent, vertex, level):
            r"""
            Initialize ``self``.

            EXAMPLES::

                sage: DA = DiGraph([[4,3], [3,2], [2,1]])
                sage: AR = DA.auslander_reiten_quiver()
                sage: TestSuite(AR(1, 3)).run()

                sage: D = DiGraph([[1,2,'a'], [1,2,'b']], multiedges=True)
                sage: AR = D.auslander_reiten_quiver()
                sage: TestSuite(AR(2, 3)).run()
                sage: TestSuite(AR(1, -4)).run()
            """
            self._vertex = vertex
            self._level = ZZ(level)
            Element.__init__(self, parent)

        def _repr_(self):
            r"""
            Return a string representation of ``self``.

            EXAMPLES::

                sage: DA = DiGraph([[4,3], [3,2], [2,1]])
                sage: AR = DA.auslander_reiten_quiver()
                sage: AR(1, 3)
                <1, 3>
            """
            return f"<{self._vertex}, {self._level}>"

        def _latex_(self):
            r"""
            Return a latex representation of ``self``.

            EXAMPLES::

                sage: DA = DiGraph([[4,3], [3,2], [2,1]])
                sage: AR = DA.auslander_reiten_quiver()
                sage: node = AR(2, 2)
                sage: latex(node)
                \left\langle 2, 2 \right\rangle
                sage: AR.options.latex = "dimension_vector"
                sage: latex(node)
                v_{2} + v_{3}
                sage: AR.options.latex = "both"
                sage: latex(node)
                \begin{gathered} \left\langle 2, 2 \right\rangle \\ v_{2} + v_{3} \end{gathered}
                sage: AR.options._reset()
            """
            from sage.misc.latex import latex
            node = r"\left\langle {}, {} \right\rangle".format(latex(self._vertex), self._level)
            latex_option = self.parent().options.latex
            if latex_option == "node":
                return node
            dim_vec = latex(self.dimension_vector())
            if latex_option == "dimension_vector":
                return dim_vec
            return r"\begin{{gathered}} {} \\ {} \end{{gathered}}".format(node, dim_vec)

        def _richcmp_(self, other, op):
            r"""
            Rich comparison of ``self`` to ``other`` by ``op``.

            EXAMPLES::

                sage: DA = DiGraph([[2,3], [2,1]])
                sage: AR = DA.auslander_reiten_quiver()
                sage: sorted(AR)
                [<1, 1>, <2, 1>, <3, 1>, <1, 2>, <2, 2>, <3, 2>]
            """
            return richcmp((self._level, self._vertex), (other._level, other._vertex), op)

        def __hash__(self):
            r"""
            Return the hash of ``self``.

            EXAMPLES::

                sage: DA = DiGraph([[2,3], [2,1]])
                sage: AR = DA.auslander_reiten_quiver()
                sage: node = AR(1, 2)
                sage: hash(node) == hash((2, 1))
                True
            """
            return hash((self._level, self._vertex))

        def vertex(self):
            r"""
            Return the vertex of the quiver corresponding to ``self``.

            EXAMPLES::

                sage: DA = DiGraph([[2,3], [2,1]])
                sage: AR = DA.auslander_reiten_quiver()
                sage: node = AR(1, 2)
                sage: node.vertex()
                1
            """
            return self._vertex

        def level(self):
            r"""
            Return the level of ``self``.

            EXAMPLES::

                sage: DA = DiGraph([[2,3], [2,1]])
                sage: AR = DA.auslander_reiten_quiver()
                sage: node = AR(1, 2)
                sage: node.level()
                2
            """
            return self._level

        def translation(self):
            r"""
            Return the AR translation of ``self``.

            EXAMPLES::

                sage: DA = DiGraph([[4,3], [3,2], [2,1]])
                sage: AR = DA.auslander_reiten_quiver()
                sage: node = AR(1, 1)
                sage: node.translation() is None
                True
                sage: node = AR(1, 2)
                sage: node.translation()
                <1, 1>
            """
            if self._level == 1:
                return None
            dim_vecs = self.parent()._dim_vecs_level(self._level - 1)
            if self._vertex not in dim_vecs:  # this likely never happens
                return None
            return type(self)(self.parent(), self._vertex, self._level - 1)

        def inverse_translation(self):
            r"""
            Return the inverse AR translation of ``self``.

            EXAMPLES::

                sage: DA = DiGraph([[2,3], [2,1]])
                sage: AR = DA.auslander_reiten_quiver()
                sage: node = AR(1, 1)
                sage: node.inverse_translation()
                <1, 2>
                sage: node = AR(1, 2)
                sage: node.inverse_translation() is None
                True

                sage: D = DiGraph([[1,2,'a'], [1,2,'b']], multiedges=True)
                sage: AR = D.auslander_reiten_quiver()
                sage: AR(2, -1).inverse_translation() is None
                True
            """
            if self._level == -1:
                return None
            dim_vecs = self.parent()._dim_vecs_level(self._level + 1)
            if self._vertex not in dim_vecs:
                return None
            return type(self)(self.parent(), self._vertex, self._level + 1)

        def dimension_vector(self):
            r"""
            Return the dimension vector of ``self``.

            EXAMPLES::

                sage: D = DiGraph([[1,2,'a'], [1,2,'b']], multiedges=True)
                sage: AR = D.auslander_reiten_quiver()
                sage: node = AR(2, -4)
                sage: node.dimension_vector()
                8*v1 + 7*v2
            """
            return self.parent()._dim_vecs_level(self._level)[self._vertex]


def detect_dynkin_quiver(quiver):
    """
    Determine if ``quiver`` is a finite type Dynkin quiver.

    EXAMPLES::

        sage: from sage.quivers.ar_quiver import detect_dynkin_quiver
        sage: D = DiGraph([[1,2], [2,3], [3, 4], [4,0], ['a','b'], ['b','c'], ['c','d'], ['c','e']])
        sage: detect_dynkin_quiver(D)
        D5xA5

        sage: D = DiGraph([[1,2,'a'], [1,2,'b']], multiedges=True)
        sage: detect_dynkin_quiver(D) is None
        True
        sage: D = DiGraph([[1, 2], [2, 3], [1, 3]])
        sage: detect_dynkin_quiver(D) is None
        True
        sage: D = DiGraph([[1,2], [1,3], [1,4], [1,5]])
        sage: detect_dynkin_quiver(D) is None
        True
        sage: D = DiGraph([[1,2], [2,3], [2,4], [4,5], [6,4]])
        sage: detect_dynkin_quiver(D) is None
        True
        sage: D = DiGraph([[1,2], [2,3], [3,4], [4,5], [5,6], [6,7], [7,8], [8,9], [0,3]])
        sage: detect_dynkin_quiver(D) is None
        True
        sage: D = DiGraph([[1,2], [2,3], [3,4], [4,5], [5,6], [6,7], [0,4]])
        sage: detect_dynkin_quiver(D) is None
        True
    """
    dynkin_type = []
    for Q in quiver.connected_components_subgraphs():
        if Q.has_multiple_edges():
            return None
        G = Q.to_undirected()
        if G.is_path():
            dynkin_type.append(['A', Q.num_verts()])
            continue
        degthree = G.vertices(degree=3)
        if len(degthree) != 1:
            return None
        G = G.copy(immutable=False)
        G.delete_vertex(degthree[0])
        path_lengths = sorted(G.connected_components_sizes())
        if len(path_lengths) != 3:
            return None
        if path_lengths[:2] == [1, 1]:
            dynkin_type.append(['D', G.num_verts() + 1])
        elif path_lengths[:2] == [1, 2] and path_lengths[2] in [2, 3, 4]:
            dynkin_type.append(['E', G.num_verts() + 1])
        else:
            return None
    if len(dynkin_type) == 1:
        return CartanType(dynkin_type[0])
    return CartanType(dynkin_type)
