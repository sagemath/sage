# sage.doctest: optional - sage.graphs
r"""
Graphic Matroids

Let `G = (V,E)` be a graph and let `C` be the collection of the edge sets
of cycles in `G`. The corresponding graphic matroid `M(G)` has groundset `E`
and circuits `C`.

Construction
============

The recommended way to create a graphic matroid is by using the
:func:`Matroid() <sage.matroids.constructor.Matroid>` function, with a
graph `G` as input. This function can accept many different kinds of input
to get a graphic matroid if the ``graph`` keyword is used, similar to the
:func:`Graph() <sage.graphs.graph.Graph>` constructor. However,
invoking the class directly is possible too. To get access to it, type::

    sage: from sage.matroids.advanced import *

See also :mod:`sage.matroids.advanced`.

Graphic matroids do not have a representation matrix or any of the
functionality of regular matroids. It is possible to get an instance of the
:class:`~sage.matroids.linear_matroid.RegularMatroid` class by using the
``regular`` keyword when constructing the matroid. It is also possible to cast
a :class:`GraphicMatroid` as a :class:`RegularMatroid` with the
:meth:`~sage.matroids.graphic_matroids.GraphicMatroid.regular_matroid`
method::

    sage: M1 = Matroid(graphs.DiamondGraph(), regular=True)
    sage: M2 = Matroid(graphs.DiamondGraph())
    sage: M3 = M2.regular_matroid()

Below are some examples of constructing a graphic matroid.

::

    sage: from sage.matroids.advanced import *
    sage: edgelist = [(0, 1, 'a'), (0, 2, 'b'), (1, 2, 'c')]
    sage: G = Graph(edgelist)
    sage: M1 = Matroid(G)
    sage: M2 = Matroid(graph=edgelist)
    sage: M3 = Matroid(graphs.CycleGraph(3))
    sage: M1 == M3
    False
    sage: M1.is_isomorphic(M3)
    True
    sage: M1.equals(M2)
    True
    sage: M1 == M2
    True
    sage: isinstance(M1, GraphicMatroid)
    True
    sage: isinstance(M1, RegularMatroid)
    False

Note that if there is not a complete set of unique edge labels, and
there are no parallel edges, then vertex tuples will be used for the
groundset. The user may wish to override this by specifying the
groundset, as the vertex tuples will not be updated if the matroid is
modified::

    sage: G = graphs.DiamondGraph()
    sage: M1 = Matroid(G)
    sage: N1 = M1.contract((0,1))
    sage: N1.graph().edges_incident(0, sort=True)
    [(0, 2, (0, 2)), (0, 2, (1, 2)), (0, 3, (1, 3))]
    sage: M2 = Matroid(range(G.num_edges()), G)
    sage: N2 = M2.contract(0)
    sage: N1.is_isomorphic(N2)
    True

AUTHORS:

- Zachary Gershkoff (2017-07-07): initial version
"""

# ****************************************************************************
#       Copyright (C) 2017 Zachary Gershkoff <zgersh2@lsu.edu>
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matroids.matroid cimport Matroid
from copy import copy
from sage.matroids.utilities import newlabel, split_vertex, sanitize_contractions_deletions
from itertools import combinations
from sage.rings.integer import Integer
from sage.sets.disjoint_set cimport DisjointSet_of_hashables

cdef class GraphicMatroid(Matroid):
    r"""
    The graphic matroid class.

    INPUT:

    - ``G`` -- :class:`Graph`
    - ``groundset`` -- list (optional); in 1-1 correspondence with ``G.edge_iterator()``

    OUTPUT: :class:`GraphicMatroid` where the groundset elements are the edges of `G`

    .. NOTE::

        If a disconnected graph is given as input, the instance of :class:`GraphicMatroid`
        will connect the graph components and store this as its graph.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = GraphicMatroid(graphs.BullGraph()); M
        Graphic matroid of rank 4 on 5 elements
        sage: N = GraphicMatroid(graphs.CompleteBipartiteGraph(3,3)); N
        Graphic matroid of rank 5 on 9 elements

    A disconnected input will get converted to a connected graph internally::

        sage: G1 = graphs.CycleGraph(3); G2 = graphs.DiamondGraph()
        sage: G = G1.disjoint_union(G2)
        sage: len(G)
        7
        sage: G.is_connected()
        False
        sage: M = GraphicMatroid(G)
        sage: M
        Graphic matroid of rank 5 on 8 elements
        sage: H = M.graph()
        sage: H
        Looped multi-graph on 6 vertices
        sage: H.is_connected()
        True
        sage: M.is_connected()
        False

    You can still locate an edge using the vertices of the input graph::

        sage: G1 = graphs.CycleGraph(3); G2 = graphs.DiamondGraph()
        sage: G = G1.disjoint_union(G2)
        sage: M = Matroid(G)
        sage: H = M.graph()
        sage: vm = M.vertex_map()
        sage: (u, v, l) = G.random_edge()
        sage: H.has_edge(vm[u], vm[v])
        True
    """

    # necessary (__init__, groundset, _rank)

    def __init__(self, G, groundset=None):
        """
        See class definition for full documentation.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: G1 = graphs.CycleGraph(3); G2 = graphs.DiamondGraph()
            sage: G = G1.disjoint_union(G2)
            sage: M = GraphicMatroid(G); M
            Graphic matroid of rank 5 on 8 elements
            sage: M.graph()
            Looped multi-graph on 6 vertices
            sage: M.graph().is_connected()
            True
            sage: M.is_connected()
            False

        TESTS::

            sage: TestSuite(M).run(verbose=True)
            running ._test_category() . . . pass
            running ._test_new() . . . pass
            running ._test_not_implemented_methods() . . . pass
            running ._test_pickling() . . . pass
        """
        from sage.graphs.graph import Graph

        if groundset is None:
            # Try to construct a groundset based on the edge labels.
            # If that fails, use range() to come up with a groundset.
            groundset = G.edge_labels()

        groundset_set = frozenset(groundset)

        # if the provided groundset is incomplete, it gets overwritten
        # invalidate ``None`` as label
        if None in groundset_set or len(groundset_set) != G.num_edges():
            groundset = range(G.num_edges())
            groundset_set = frozenset(groundset)

        self._groundset = groundset_set

        # Map vertices on input graph to vertices in self._G
        self._vertex_map = {v: v for v in G.vertex_iterator()}
        comps = G.connected_components(sort=False)
        while len(comps) > 1:
            comp = comps.pop()
            v1 = comps[-1][-1]
            v2 = comp[0]
            self._vertex_map[v2] = v1
            comps[-1].extend(comp)

        # Construct a graph and assign edge labels corresponding to the groundset
        edge_list = []
        for i, e in enumerate(G.edge_iterator()):
            # the ordering from edge_labels() respects edge_iterator() and not edges()
            edge_list.append((self._vertex_map[e[0]],
                              self._vertex_map[e[1]], groundset[i]))
        # If the matroid is empty, have the internal graph be a single vertex
        if edge_list:
            self._G = <GenericGraph_pyx?> Graph(edge_list, loops=True, multiedges=True,
                                                weighted=True, data_structure='static_sparse')
        else:
            self._G = <GenericGraph_pyx?> Graph(1, loops=True, multiedges=True,
                                                weighted=True, data_structure='static_sparse')
        # Map groundset elements to graph edges:
        # The edge labels should already be the elements.
        self._groundset_edge_map = ({l: (u, v) for (u, v, l) in self._G.edge_iterator()})

    cpdef frozenset groundset(self):
        """
        Return the groundset of the matroid as a frozenset.

        EXAMPLES::

            sage: M = Matroid(graphs.DiamondGraph())
            sage: sorted(M.groundset())
            [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]
            sage: G = graphs.CompleteGraph(3).disjoint_union(graphs.CompleteGraph(4))
            sage: M = Matroid(range(G.num_edges()), G); sorted(M.groundset())
            [0, 1, 2, 3, 4, 5, 6, 7, 8]
            sage: M = Matroid(Graph([(0, 1, 'a'), (0, 2, 'b'), (0, 3, 'c')]))
            sage: sorted(M.groundset())
            ['a', 'b', 'c']
        """
        return self._groundset

    cpdef int _rank(self, frozenset X):
        """
        Return the rank of a set ``X``.

        This method does no checking on ``X``, and ``X`` may be assumed to have
        the same interface as ``frozenset``.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface

        OUTPUT: the rank of `X` in the matroid

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: edgelist = [(0,0,0), (0,1,1), (0,2,2), (0,3,3), (1,2,4), (1,3,5)]
            sage: M = GraphicMatroid(Graph(edgelist, loops=True, multiedges=True))
            sage: M.rank([0])
            0
            sage: M.rank([1,2])
            2
            sage: M.rank([1,2,4])
            2
            sage: M.rank(M.groundset())
            3
            sage: edgelist = [(0,0,0), (1,2,1), (1,2,2), (2,3,3)]
            sage: M = GraphicMatroid(Graph(edgelist, loops=True, multiedges=True))
            sage: M.rank(M.groundset())
            2
            sage: M.rank([0,3])
            1
        """
        cdef DisjointSet_of_hashables DS_vertices
        cdef list edges = self.groundset_to_edges(X)
        cdef set vertices = set([u for (u, v, l) in edges]).union(
                                [v for (u, v, l) in edges])
        # This counts components:
        DS_vertices = DisjointSet_of_hashables(vertices)
        for (u, v, l) in edges:
            DS_vertices.union(u, v)
        return (len(vertices) - DS_vertices.number_of_subsets())

    # representation:

    def _repr_(self):
        """
        Return a string representation of the matroid.

        EXAMPLES::

            sage: M = Matroid(graphs.CompleteGraph(5)); M
            Graphic matroid of rank 4 on 10 elements
            sage: G = Graph([(0, 0), (0, 1), (0, 2), (1, 1), (2, 2)], loops=True)
            sage: M = Matroid(G); M
            Graphic matroid of rank 2 on 5 elements
        """
        r = self._rank(self._groundset)
        n = len(self._groundset)
        return f'Graphic matroid of rank {r} on {n} elements'

    # comparison:

    cpdef _vertex_stars(self):
        """
        Compute the set of edge labels around each vertex.

        Internal method for hashing purposes.

        OUTPUT: a ``frozenset`` of frozensets containing the edge labels
        around each vertex

        EXAMPLES::

            sage: M = Matroid(range(5), graphs.DiamondGraph())
            sage: sorted(M._vertex_stars(), key=str)
            [frozenset({0, 1}),
             frozenset({0, 2, 3}),
             frozenset({1, 2, 4}),
             frozenset({3, 4})]

            sage: N = Matroid(range(5), graphs.BullGraph())
            sage: sorted(N._vertex_stars(), key=str)
            [frozenset({0, 1}),
             frozenset({0, 2, 3}),
             frozenset({1, 2, 4}),
             frozenset({3}),
             frozenset({4})]
        """
        star_list = []
        for v in self._G.vertices(sort=False):
            star = [l for (_, _, l) in self._G.edges_incident(v)]
            star_list.append(frozenset(star))
        return frozenset(star_list)

    def __hash__(self):
        r"""
        Return an invariant of the matroid.

        This function is called when matroids are added to a set. It is very
        desirable to override it so it can distinguish matroids on the same
        groundset, which is a very typical use case!

        .. WARNING::

            This method is linked to __richcmp__ (in Cython) and __cmp__ or
            __eq__/__ne__ (in Python). If you override one, you should (and in
            Cython: MUST) override the other!

        EXAMPLES::

            sage: M = Matroid(graphs.CompleteGraph(3))
            sage: N = Matroid(graphs.CycleGraph(3))
            sage: O = Matroid(graphs.ButterflyGraph())
            sage: hash(M) == hash(N)
            True
            sage: hash(O) == hash(N)
            False
            sage: P = Matroid(Graph([(0, 1, 'a'), (0, 2, 'b'), (1, 2, 'c')]))
            sage: hash(P) == hash(M)
            False
        """
        return hash(self._vertex_stars())

    def __eq__(self, other):
        """
        Compare two matroids.

        For two graphic matroids to be equal, all attributes of the underlying
        graphs must be equal.

        INPUT:

        - ``other`` -- matroid

        OUTPUT: ``True`` if ``self`` and ``other`` have the same graph;
        ``False`` otherwise

        EXAMPLES::

            sage: M = Matroid(graphs.CompleteGraph(3))
            sage: N = Matroid(graphs.CycleGraph(3))
            sage: O = Matroid(graphs.ButterflyGraph())
            sage: P = Matroid(Graph([(0, 1, 'a'), (0, 2, 'b'), (1, 2, 'c')]))
            sage: M == N
            True
            sage: M == O
            False
            sage: M == P
            False

        A more subtle example where the vertex labels differ::

            sage: G1 = Graph([(0,1,0),(0,2,1),(1,2,2)])
            sage: G2 = Graph([(3,4,3),(3,5,4),(4,5,5),(4,6,6),(5,6,7)])
            sage: G = G1.disjoint_union(G2)
            sage: H = G2.disjoint_union(G1)
            sage: Matroid(G) == Matroid(H)
            False
            sage: Matroid(G).equals(Matroid(H))
            True

        Same except for vertex labels::

            sage: G1 = Graph([(0,1,0),(1,2,1),(2,0,2)])
            sage: G2 = Graph([(3,4,0),(4,5,1),(5,3,2)])
            sage: Matroid(G1) == Matroid(G2)
            False
        """
        # Graph.__eq__() will ignore edge labels unless we turn on weighted()
        # This will be done in __init__()
        if not isinstance(other, GraphicMatroid):
            return False
        return (self._G == other._G)

    def __ne__(self, other):
        """
        Compare two matroids.

        INPUT:

        - ``other`` -- matroid

        OUTPUT: ``False`` if ``self`` and ``other`` have the same graph;
        ``True`` otherwise

        EXAMPLES::

            sage: M = Matroid(range(4), graphs.CycleGraph(4))
            sage: N = Matroid(range(4), graphs.CompleteBipartiteGraph(2,2))
            sage: O = Matroid(graphs.PetersenGraph())
            sage: M != N
            True
            sage: M.equals(N)
            True
            sage: M != O
            True
        """
        return (not self == other)

    # copying, loading, saving

    def __reduce__(self):
        """
        Save the matroid for later reloading.

        EXAMPLES::

            sage: M = Matroid(graphs.PetersenGraph())
            sage: M == loads(dumps(M))
            True
            sage: loads(dumps(M))
            Graphic matroid of rank 9 on 15 elements
        """
        from sage.matroids.unpickling import unpickle_graphic_matroid
        data = (self._G, self.get_custom_name())
        version = 0
        return unpickle_graphic_matroid, (version, data)

    # overrides

    cpdef _minor(self, contractions, deletions):
        """
        Return a minor.

        INPUT:

        - ``contractions`` -- frozenset; subset of ``self.groundset()`` to be contracted
        - ``deletions`` -- frozenset; subset of ``self.groundset()`` to be deleted

        Assumptions: contractions are independent, deletions are coindependent,
        contractions and deletions are disjoint.

        OUTPUT: :class:`GraphicMatroid`

        EXAMPLES::

            sage: M = matroids.CompleteGraphic(5)
            sage: M._minor(deletions=frozenset([0,1,2]), contractions=frozenset([]))
            Graphic matroid of rank 4 on 7 elements
            sage: M._minor(deletions=frozenset([]), contractions=frozenset([0,1,2]))
            Graphic matroid of rank 1 on 7 elements
            sage: M = Matroid(range(15), graphs.PetersenGraph())
            sage: N = M._minor(deletions=frozenset([0, 3, 5, 9]),
            ....:              contractions=frozenset([1, 2, 11]))
            sage: N
            Graphic matroid of rank 6 on 8 elements
        """
        cdef GenericGraph_pyx g = self.graph()
        cdef list cont_edges = self._groundset_to_edges(contractions)
        cdef list del_edges = self._groundset_to_edges(deletions)
        # deletions first so contractions don't mess up the vertices
        g.delete_edges(del_edges)
        g.contract_edges(cont_edges)
        return GraphicMatroid(g)

    cpdef _has_minor(self, N, bint certificate=False):
        """
        Check if the matroid has a minor isomorphic to `M(H)`.

        INPUT:

        - ``N`` -- matroid
        - ``certificate`` -- boolean (default: ``False``); if ``True``, returns
          the certificate isomorphism from the minor of ``self`` to ``N``

        OUTPUT:

        boolean, or tuple if the ``certificate`` option is used. If
        ``certificate`` is ``True``, then the output will either be
        ``False, None`` or ``True, (X, Y, dic) where ``N`` is isomorphic to
        ``self.minor(X, Y)``, and ``dic`` is an isomorphism between ``N`` and
        ``self.minor(X, Y)``.

        EXAMPLES::

            sage: M = Matroid(range(9), graphs.CompleteBipartiteGraph(3, 3))
            sage: N = Matroid(range(3), graphs.CycleGraph(3))
            sage: N1 = Matroid(range(3), graph=graphs.CycleGraph(3),
            ....: regular=True)
            sage: _, cert = M._has_minor(N1, certificate=True)
            sage: Mp = M.minor(cert[0], cert[1])
            sage: N.is_isomorphism(Mp, cert[2])
            True
            sage: M._has_minor(N)
            True
            sage: M._has_minor(N1)
            True
            sage: _, cert = M._has_minor(N, certificate=True)
            sage: Mp = M.minor(cert[0], cert[1])
            sage: N.is_isomorphism(Mp, cert[2])
            True

        ::

            sage: M = matroids.CompleteGraphic(6)
            sage: N = Matroid(range(8), graphs.WheelGraph(5))
            sage: M._has_minor(N)
            True
            sage: _, cert = M._has_minor(N, certificate=True)
            sage: Mp = M.minor(cert[0], cert[1])
            sage: N.is_isomorphism(Mp, cert[2])
            True
            sage: N.has_minor(M)
            False
            sage: N.has_minor(M, certificate=True)
            (False, None)

        If the matroids are not 3-connected, then the default matroid
        algorithms are used::

            sage: M = matroids.CompleteGraphic(6)
            sage: N = Matroid(graphs.CycleGraph(4))
            sage: M.has_minor(N)
            True
            sage: N.has_minor(M)
            False
        """
        # The graph minor algorithm is faster but it doesn't make sense
        # to use it if M(H) is not 3-connected, because of all the possible
        # Whitney switches or 1-sums that will give the same matroid.
        if isinstance(N, GraphicMatroid) and N.is_3connected():
            # Graph.minor() does not work with multigraphs
            G = self.graph()
            G.allow_loops(False)
            G.allow_multiple_edges(False)
            H = N.graph()
            H.allow_loops(False)
            H.allow_multiple_edges(False)

            try:
                # Graph.minor() returns a certificate if there is one
                # and a ValueError if there isn't.
                cert = G.minor(H)
            except ValueError:
                if certificate:
                    return (False, None)
                else:
                    return False

            if certificate:
                # This is where it gets complicated.
                # The Graph.minor() method gives a dictionary of vertices
                # as its certificate. There is currently no easy way to
                # determine the edges.
                # From the dictionary, we can get an idea of what the
                # contractions are, and what vertices are not used.
                # So we'll merge the appropriate vertices, delete the
                # unused vertices, and pass to Matroid._has_minor().

                # Determine contractions:
                vertices_for_minor = cert.values()
                contractions = []
                for vertex_list in vertices_for_minor:
                    S = G.subgraph(vertex_list)
                    X = S.edge_labels()
                    contractions.extend(self.max_independent(X))

                # determine deletions:
                from itertools import chain
                deletions = []
                big_vertex_list = list(chain.from_iterable(vertices_for_minor))
                for v in G.vertices(sort=False):
                    if v not in big_vertex_list:
                        deletions.extend([l for (u0, v0, l) in G.edges_incident(v)])

                # take contractions and deletions with what we have so far
                # then use method from abstract matroid class
                conset, delset = sanitize_contractions_deletions(self, contractions, deletions)
                M = self._minor(contractions=conset, deletions=delset)
                should_be_true, elements = Matroid._has_minor(M, N, certificate=True)

                # elements is a tuple (contractions, deletions, dict)
                # There should be no more contractions
                delset = set(delset)
                delset.update(elements[1])
                return (True, (conset, frozenset(delset), elements[2]))

            else:
                return True
        else:
            # otherwise send it to regular matroids
            M = self.regular_matroid()
            if isinstance(N, GraphicMatroid):
                N = N.regular_matroid()
            return M._has_minor(N, certificate=certificate)

    cpdef int _corank(self, frozenset X):
        """
        Return the corank of the set `X` in the matroid.

        Internal version that does no input checking.

        INPUT:

        - ``X`` -- an iterable container of groundset elements

        OUTPUT: integer

        EXAMPLES::

            sage: M = Matroid(range(9), graphs.CompleteBipartiteGraph(3,3))
            sage: M._corank(frozenset([0,1,2]))
            2
            sage: M._corank(frozenset([1,2,3]))
            3
        """
        cdef DisjointSet_of_hashables DS_vertices
        cdef list all_vertices = self._G.vertices(sort=False)
        cdef list not_our_edges = self.groundset_to_edges(self._groundset.difference(X))
        DS_vertices = DisjointSet_of_hashables(all_vertices)
        for u, v, l in not_our_edges:
            DS_vertices.union(u, v)
        return len(X) - (DS_vertices.number_of_subsets() - Integer(1))

    cpdef bint _is_circuit(self, frozenset X):
        """
        Test if input is a circuit.

        INPUT:

        - ``X`` -- an iterable container of groundset elements

        OUTPUT: boolean

        EXAMPLES::

            sage: M = Matroid(range(5), graphs.DiamondGraph())
            sage: M._is_circuit(frozenset([0,1,2]))
            True
            sage: M._is_circuit(frozenset([0,1,2,3]))
            False
            sage: M._is_circuit(frozenset([0,1,3]))
            False
        """
        cdef GenericGraph_pyx g = self._subgraph_from_set(X)
        return g.is_cycle()

    cpdef frozenset _closure(self, frozenset X):
        """
        Return the closure of a set.

        INPUT:

        - ``X`` -- an iterable container of groundset elements

        OUTPUT: a subset of the groundset as a :class:`frozenset`

        EXAMPLES::

            sage: M = Matroid(range(5), graphs.DiamondGraph())
            sage: sorted(M._closure(frozenset([0])))
            [0]
            sage: sorted(M._closure(frozenset([0,1])))
            [0, 1, 2]
            sage: sorted(M._closure(M.groundset()))
            [0, 1, 2, 3, 4]

        TESTS:

        Make sure the closure gets loops::

            sage: edgelist = [(0, 0), (0, 1), (0, 2), (0, 3), (1, 2), (1, 2)]
            sage: M = Matroid(range(6), Graph(edgelist, loops=True, multiedges=True))
            sage: M.graph().edges(sort=True)
            [(0, 0, 0), (0, 1, 1), (0, 2, 2), (0, 3, 3), (1, 2, 4), (1, 2, 5)]
            sage: sorted(M._closure(frozenset([4])))
            [0, 4, 5]
        """
        cdef set XX = set(X)
        cdef frozenset Y = self.groundset().difference(XX)
        cdef list edgelist = self._groundset_to_edges(Y)
        cdef GenericGraph_pyx g = self._subgraph_from_set(XX)
        cdef list V = g.vertices(sort=False)
        cdef int components = g.connected_components_number()
        for e in edgelist:
            # a non-loop edge is in the closure iff both its vertices are
            # in the induced subgraph, and the edge doesn't connect components
            if e[0] in V and e[1] in V:
                g.add_edge(e)
                if g.connected_components_number() >= components:
                    XX.add(e[2])
                else:
                    g.delete_edge(e)
        # add all loops
        XX.update(set([l for (u, v, l) in self._G.loops()]))
        return frozenset(XX)

    cpdef frozenset _max_independent(self, frozenset X):
        """
        Compute a maximal independent subset.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: a subset of the groundset as a :class:`frozenset`

        EXAMPLES::

            sage: M = Matroid(range(5), graphs.DiamondGraph())
            sage: sorted(M._max_independent(M.groundset()))
            [0, 1, 3]
            sage: sorted(M._max_independent(frozenset([0,1,2])))
            [0, 1]
            sage: sorted(M._max_independent(frozenset([3,4])))
            [3, 4]
            sage: sorted(M._max_independent(frozenset([3])))
            [3]
            sage: N = M.graphic_extension(0, element='a')
            sage: sorted(N._max_independent(frozenset(['a'])))
            []
        """
        cdef DisjointSet_of_hashables DS_vertices
        cdef list edges = self.groundset_to_edges(X)
        cdef set vertices = set([u for (u, v, l) in edges])
        vertices.update([v for (u, v, l) in edges])

        cdef set our_set = set()
        DS_vertices = DisjointSet_of_hashables(vertices)
        for (u, v, l) in edges:
            if DS_vertices.find(u) != DS_vertices.find(v):
                DS_vertices.union(u, v)
                our_set.add(l)
        return frozenset(our_set)

    cpdef frozenset _max_coindependent(self, frozenset X):
        """
        Compute a maximal coindependent subset.

        INPUT:

        - ``X`` -- an iterable container of groundset elements

        OUTPUT: a subset of the groundset as a :class:`frozenset`

        EXAMPLES::

            sage: M = Matroid(range(5), graphs.DiamondGraph())
            sage: sorted(M._max_coindependent(M.groundset()))
            [2, 4]
            sage: sorted(M._max_coindependent(frozenset([2,3,4])))
            [2, 4]
            sage: N = M.graphic_extension(0, element=5)
            sage: sorted(N.max_coindependent(frozenset([0,1,2,5])))
            [1, 2, 5]
        """
        cdef DisjointSet_of_hashables DS_vertices
        cdef list edges = self.groundset_to_edges(X)
        cdef list all_vertices = self._G.vertices(sort=False)
        cdef list not_our_edges = self.groundset_to_edges(self._groundset.difference(X))

        cdef set our_set = set()
        DS_vertices = DisjointSet_of_hashables(all_vertices)
        for (u, v, l) in not_our_edges:
            DS_vertices.union(u, v)

        for (u, v, l) in edges:
            if DS_vertices.find(u) == DS_vertices.find(v):
                our_set.add(l)
            else:
                DS_vertices.union(u, v)
        return frozenset(our_set)

    cpdef frozenset _circuit(self, frozenset X):
        """
        Return a minimal dependent subset.

        INPUT:

        - ``X`` -- an iterable container of groundset elements

        OUTPUT: ``frozenset`` instance containing a subset of ``X``;
        a :exc:`ValueError` is raised if the set contains no circuit

        EXAMPLES::

            sage: M = Matroid(range(5), graphs.DiamondGraph())
            sage: sorted(M._circuit(M.groundset()))
            [0, 1, 2]
            sage: N = Matroid(range(9), graphs.CompleteBipartiteGraph(3,3))
            sage: sorted(N._circuit(frozenset([0, 1, 2, 6, 7, 8])))
            [0, 1, 6, 7]
            sage: N._circuit(frozenset([0, 1, 2]))
            Traceback (most recent call last):
            ...
            ValueError: no circuit in independent set

        TESTS:

        With two disjoint cycles in the graph::

            sage: edgelist = [(5,6), (0,1), (3,4), (1,2), (4,5), (2,0), (5,3)]
            sage: M = Matroid(range(7), Graph(edgelist))
            sage: M
            Graphic matroid of rank 5 on 7 elements
            sage: sorted(M._circuit(M.groundset()))
            [0, 1, 2]

        Giving it a long path before it finds a cycle::

            sage: edgelist = [(0,1), (1,2), (2,3), (3,4), (4,5), (4,5)]
            sage: M = Matroid(Graph(edgelist, multiedges=True))
            sage: M.graph().edges(sort=True)
            [(0, 1, 0), (1, 2, 1), (2, 3, 2), (3, 4, 3), (4, 5, 4), (4, 5, 5)]
            sage: sorted(M._circuit(M.groundset()))
            [4, 5]
        """
        cdef list edges = self.groundset_to_edges(X)
        cdef set vertices = set()
        cdef list vertex_list = []
        cdef list leaves
        cdef tuple leaf
        cdef set edge_set = set()
        cdef DisjointSet_of_hashables DS_vertices

        for (u, v, l) in edges:
            vertices.add(u)
            vertices.add(v)
        DS_vertices = DisjointSet_of_hashables(vertices)
        for (u, v, l) in edges:
            edge_set.add((u, v, l))
            if DS_vertices.find(u) != DS_vertices.find(v):
                DS_vertices.union(u, v)
            else:
                break
        else:
            raise ValueError("no circuit in independent set")

        for (u, v, l) in edge_set:
            vertex_list.extend([u, v])
        leaves = [(u, v, l) for (u, v, l) in edge_set
                  if vertex_list.count(u) == 1 or vertex_list.count(v) == 1]
        while leaves:
            for leaf in leaves:
                edge_set.remove(leaf)
                vertex_list.remove(leaf[0])
                vertex_list.remove(leaf[1])
            leaves = [(u, v, l) for (u, v, l) in edge_set
                      if vertex_list.count(u) == 1 or vertex_list.count(v) == 1]

        return frozenset([l for (u, v, l) in edge_set])

    cpdef frozenset _coclosure(self, frozenset X):
        """
        Return the coclosure of a set.

        INPUT:

        - ``X`` -- an iterable container of groundset elements

        OUTPUT: a subset of the groundset as a :class:`frozenset`

        EXAMPLES::

            sage: M = Matroid(range(5), graphs.DiamondGraph())
            sage: sorted(M._coclosure(frozenset([0])))
            [0, 1]
            sage: sorted(M._coclosure(frozenset([0,1])))
            [0, 1]
            sage: N = M.graphic_extension(0, element=5)
            sage: sorted(N._coclosure(frozenset([3])))
            [3, 4]
            sage: N = M.graphic_coextension(0, element=5)
            sage: sorted(N._coclosure(frozenset([3])))
            [3, 4, 5]
        """
        cdef GenericGraph_pyx g = self.graph()
        g.delete_edges(self._groundset_to_edges(X))
        cdef int components = g.connected_components_number()
        cdef set XX = set(X)
        cdef frozenset Y = self.groundset().difference(XX)
        for e in self._groundset_to_edges(Y):
            g.delete_edge(e)
            if g.connected_components_number() > components:
                XX.add(e[2])
            g.add_edge(e)
        return frozenset(XX)

    cpdef bint _is_closed(self, frozenset X):
        """
        Test if input is a closed set.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: boolean

        EXAMPLES::

            sage: edgelist = [(0, 0), (0, 1), (0, 2), (0, 3), (1, 2), (1, 2)]
            sage: M = Matroid(range(len(edgelist)), Graph(edgelist, loops=True,
            ....: multiedges=True))
            sage: M._is_closed(frozenset([0,4,5]))
            True
            sage: M._is_closed(frozenset([0,4]))
            False
            sage: M._is_closed(frozenset([1, 2, 3, 4 ,5]))
            False
        """
        # Take the set of vertices of the edges corresponding to the elements,
        # and check if there are other edges incident with two of those vertices.
        # Also, there must not be loops outside of X.
        cdef set XX = set(X)
        cdef set loop_labels = set([l for (u, v, l) in self._G.loops()])
        if not loop_labels.issubset(XX):
            return False

        # Remove loops from input since we don't want to count them as
        # components
        XX.difference_update(loop_labels)
        cdef list edge_list = self._groundset_to_edges(XX)

        cdef set vertex_set = set()
        cdef frozenset Y = self.groundset().difference(XX)
        cdef list edge_list2 = self._groundset_to_edges(Y)
        for e in edge_list:
            vertex_set.add(e[0])
            vertex_set.add(e[1])
        for e in edge_list2:
            if e[0] in vertex_set and e[1] in vertex_set:
                return False
        return True

    cpdef _is_isomorphic(self, other, certificate=False):
        """
        Test if ``self`` is isomorphic to ``other``.

        INPUT:

        - ``other`` -- matroid
        - ``certificate`` -- boolean

        OUTPUT:

        - If ``certificate`` is ``False``, boolean
        - If ``certificate`` is ``True``, a tuple containing a boolean and a dictionary
          giving the isomorphism or None

        EXAMPLES::

            sage: M = Matroid(range(5), graphs.DiamondGraph())
            sage: N = Matroid(graph=graphs.DiamondGraph(), regular=True)
            sage: M._is_isomorphic(N, certificate=True)
            (True, {0: (0, 1), 1: (0, 2), 2: (1, 2), 3: (1, 3), 4: (2, 3)})
            sage: O = Matroid(graphs.WheelGraph(5))
            sage: M._is_isomorphic(O, certificate=True)
            (False, None)

        ::

            sage: M1 = Matroid(range(4), graphs.CycleGraph(4))
            sage: M2 = Matroid(range(4), graphs.CompleteBipartiteGraph(2,2))
            sage: M3 = matroids.Uniform(3,4)
            sage: M1._is_isomorphic(M2)
            True
            sage: M1._is_isomorphic(M3)
            True

        ::

            sage: edgelist = [(0, 1, 'a'), (0, 2, 'b'), (0, 3, 'c'),
            ....:             (1, 2, 'd'), (1, 3, 'e'), (2, 3, 'f')]
            sage: M = Matroid(Graph(edgelist))
            sage: N = Matroid(range(6), graphs.WheelGraph(4))
            sage: M._is_isomorphic(N, certificate=True)
            (True, {'a': 2, 'b': 4, 'c': 5, 'd': 0, 'e': 1, 'f': 3})
            sage: N._is_isomorphic(M, certificate=True)
            (True, {0: 'd', 1: 'e', 2: 'a', 3: 'f', 4: 'b', 5: 'c'})
            sage: O = Matroid(range(6), graphs.CycleGraph(6))
            sage: M._is_isomorphic(O)
            False
        """
        # Check for 3-connectivity so we don't have to worry about Whitney twists
        if isinstance(other, GraphicMatroid) and other.is_3connected():
            G = self.graph()
            H = other.graph()
            G.allow_loops(False)
            G.allow_multiple_edges(False)
            H.allow_loops(False)
            H.allow_multiple_edges(False)

            result = G.is_isomorphic(H, certificate=certificate)
            if not certificate or result[0] is False:
                return result
            # If they are isomorphic and the user wants a certificate,
            # result[1] is a dictionary of vertices.
            # We need to translate this to edge labels.
            vertex_certif = result[1]
            elt_certif = {}
            for (u, v, l) in G.edge_iterator():
                l_maps_to = H.edge_label(vertex_certif[u], vertex_certif[v])
                elt_certif[l] = l_maps_to
            return (True, elt_certif)

        else:
            M = self.regular_matroid()
            if isinstance(other, GraphicMatroid):
                other = other.regular_matroid()
            if certificate:
                # iso0: isomorphism from M and self -- in this order,
                # to prevent an infinite recursion.
                iso0 = M._is_isomorphic(self, certificate=certificate)[1]
                # Now invert iso0 to get iso1, an isomorphism from self to M.
                iso1 = {iso0[e]: e for e in iso0}
                # iso2: isomorphism from M and other.
                isomorphic, iso2 = M._is_isomorphic(other, certificate=certificate)
                if not isomorphic:
                    return (False, None)
                # Compose iso1 and iso2, to go from self to other.
                return (True, {e: iso2[iso1[e]] for e in iso1})
            return M._is_isomorphic(other)

    cpdef _isomorphism(self, other):
        """
        Return isomorphism from ``self`` to ``other``, if such an isomorphism
        exists.

        Internal version that performs no checks on input.

        INPUT:

        - ``other`` -- matroid

        OUTPUT: dictionary or ``None``

        EXAMPLES::

            sage: M1 = Matroid(range(4), graphs.CycleGraph(4))
            sage: M2 = Matroid(range(4), graphs.CompleteBipartiteGraph(2,2))
            sage: M1._isomorphism(matroids.catalog.BetsyRoss())
            sage: M1._isomorphism(M2)
            {0: 0, 1: 1, 2: 2, 3: 3}
            sage: M3 = matroids.Uniform(3,4)
            sage: M1._isomorphism(M3)
            {0: 0, 1: 1, 2: 2, 3: 3}

        ::

            sage: edgelist = [(0, 1, 'a'), (0, 2, 'b'), (0, 3, 'c'),
            ....:             (1, 2, 'd'), (1, 3, 'e'), (2, 3, 'f')]
            sage: M = Matroid(Graph(edgelist))
            sage: N = Matroid(range(6), graphs.WheelGraph(4))
            sage: M._isomorphism(N)
            {'a': 2, 'b': 4, 'c': 5, 'd': 0, 'e': 1, 'f': 3}
            sage: O = Matroid(Graph(edgelist), regular=True)
            sage: iso = M._isomorphism(O)
            sage: M.is_isomorphism(O, iso)
            True
        """
        return self.is_isomorphic(other, certificate=True)[1]

    cpdef bint is_valid(self):
        """
        Test if the data obey the matroid axioms.

        Since a graph is used for the data, this is always the case.

        OUTPUT: ``True``

        EXAMPLES::

            sage: M = matroids.CompleteGraphic(4); M
            M(K4): Graphic matroid of rank 3 on 6 elements
            sage: M.is_valid()
            True
        """
        return True

    cpdef bint is_graphic(self):
        r"""
        Return if ``self`` is graphic.

        This is trivially ``True`` for a :class:`GraphicMatroid`.

        EXAMPLES::

            sage: M = Matroid(graphs.PetersenGraph())
            sage: M.is_graphic()
            True
        """
        return True

    cpdef bint is_regular(self):
        r"""
        Return if ``self`` is regular.

        This is always ``True`` for a :class:`GraphicMatroid`.

        EXAMPLES::

            sage: M = Matroid(graphs.DesarguesGraph())
            sage: M.is_regular()
            True
        """
        return True

    # graphic methods

    cpdef graph(self):
        """
        Return the graph that represents the matroid.

        The graph will always have loops and multiedges enabled.

        OUTPUT: graph

        EXAMPLES::

            sage: M = Matroid(Graph([(0, 1, 'a'), (0, 2, 'b'), (0, 3, 'c')]))
            sage: M.graph().edges(sort=True)
            [(0, 1, 'a'), (0, 2, 'b'), (0, 3, 'c')]
            sage: M = Matroid(graphs.CompleteGraph(5))
            sage: M.graph()
            Looped multi-graph on 5 vertices
        """
        # Return a mutable graph
        return self._G.copy(data_structure='sparse')

    cpdef vertex_map(self):
        """
        Return a dictionary mapping the input vertices to the current vertices.

        The graph for the matroid is always connected. If the constructor is
        given a graph with multiple components, it will connect them. The
        Python dictionary given by this method has the vertices from the
        input graph as keys, and the corresponding vertex label after any
        merging as values.

        OUTPUT: dictionary

        EXAMPLES::

            sage: G = Graph([(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5),
            ....: (6, 7), (6, 8), (7, 8), (8, 8), (7, 8)], multiedges=True, loops=True)
            sage: M = Matroid(range(G.num_edges()), G)
            sage: M.graph().edges(sort=True)
            [(0, 1, 0),
             (0, 2, 1),
             (1, 2, 2),
             (1, 4, 3),
             (1, 5, 4),
             (4, 5, 5),
             (4, 7, 6),
             (4, 8, 7),
             (7, 8, 8),
             (7, 8, 9),
             (8, 8, 10)]
            sage: M.vertex_map()
            {0: 0, 1: 1, 2: 2, 3: 1, 4: 4, 5: 5, 6: 4, 7: 7, 8: 8}
        """
        return copy(self._vertex_map)

    cpdef list groundset_to_edges(self, X):
        """
        Return a list of edges corresponding to a set of groundset elements.

        INPUT:

        - ``X`` -- subset of the groundset

        OUTPUT: list of graph edges

        EXAMPLES::

            sage: M = Matroid(range(5), graphs.DiamondGraph())
            sage: M.groundset_to_edges([2, 3, 4])
            [(1, 2, 2), (1, 3, 3), (2, 3, 4)]
            sage: M.groundset_to_edges([2, 3, 4, 5])
            Traceback (most recent call last):
            ...
            ValueError: input must be a subset of the groundset
        """
        for x in X:
            if x not in self._groundset:
                raise ValueError("input must be a subset of the groundset")
        return self._groundset_to_edges(X)

    cpdef _groundset_to_edges(self, X):
        """
        Return a list of edges corresponding to a set of groundset elements.

        INPUT:

        - ``X`` -- subset of the groundset

        OUTPUT: list of graph edges

        EXAMPLES::

            sage: M = Matroid(range(5), graphs.DiamondGraph())
            sage: M._groundset_to_edges([2, 3, 4])
            [(1, 2, 2), (1, 3, 3), (2, 3, 4)]
        """
        return [(self._groundset_edge_map[x][0], self._groundset_edge_map[x][1], x) for x in X]

    cpdef subgraph_from_set(self, X):
        """
        Return the subgraph corresponding to the matroid restricted to `X`.

        INPUT:

        - ``X`` -- subset of the groundset

        OUTPUT: graph

        EXAMPLES::

            sage: M = Matroid(range(5), graphs.DiamondGraph())
            sage: M.subgraph_from_set([0,1,2])
            Looped multi-graph on 3 vertices
            sage: M.subgraph_from_set([3,4,5])
            Traceback (most recent call last):
            ...
            ValueError: input must be a subset of the groundset
        """
        for x in X:
            if x not in self._groundset:
                raise ValueError("input must be a subset of the groundset")
        return self._subgraph_from_set(X)

    cpdef _subgraph_from_set(self, X):
        """
        Return the subgraph corresponding to `M` restricted to `X`.

        INPUT:

        - ``X`` -- subset of the groundset

        OUTPUT: graph

        EXAMPLES::

            sage: M = Matroid(range(5), graphs.DiamondGraph())
            sage: M._subgraph_from_set([0, 1, 2])
            Looped multi-graph on 3 vertices
        """
        from sage.graphs.graph import Graph
        edge_list = self._groundset_to_edges(X)
        return Graph(edge_list, loops=True, multiedges=True)

    cpdef graphic_extension(self, u, v=None, element=None):
        """
        Return a graphic matroid extended by a new element.

        A new edge will be added between ``u`` and ``v``. If ``v`` is not
        specified, then a loop is added on ``u``.

        INPUT:

        - ``u`` -- vertex in the matroid's graph
        - ``v`` -- (optional) another vertex
        - ``element`` -- (optional) the label of the new element

        OUTPUT:

        A :class:`GraphicMatroid` with the specified element added. Note that if
        ``v`` is not specified or if ``v`` is ``u``, then the new element will
        be a loop. If the new element's label is not specified, it will be
        generated automatically.

        EXAMPLES::

            sage: M = matroids.CompleteGraphic(4)
            sage: M1 = M.graphic_extension(0,1,'a'); M1
            Graphic matroid of rank 3 on 7 elements
            sage: list(M1.graph().edge_iterator())
            [(0, 1, 'a'), (0, 1, 0), (0, 2, 1), (0, 3, 2), (1, 2, 3), (1, 3, 4), (2, 3, 5)]
            sage: M2 = M1.graphic_extension(3); M2
            Graphic matroid of rank 3 on 8 elements

        ::

            sage: M = Matroid(range(10), graphs.PetersenGraph())
            sage: sorted(M.graphic_extension(0, 'b', 'c').graph().vertex_iterator(), key=str)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 'b']
            sage: M.graphic_extension('a', 'b', 'c').graph().vertices(sort=False)
            Traceback (most recent call last):
            ...
            ValueError: u must be an existing vertex

        TESTS::

            sage: M = Matroid(graphs.EmptyGraph())
            sage: M.graphic_extension(0)
            Graphic matroid of rank 0 on 1 elements
            sage: M.graphic_extension(0, 1, 'a')
            Graphic matroid of rank 1 on 1 elements
        """
        # This will possibly make a coloop if v is a new vertex
        if element is None:
            element = newlabel(self.groundset())
        elif element in self.groundset():
            raise ValueError("cannot extend by element already in groundset")
        if v is None:
            v = u
        G = self.graph()
        if u not in G:
            raise ValueError("u must be an existing vertex")
        G.add_edge(u, v, element)
        return GraphicMatroid(G)

    def graphic_extensions(self, element=None, vertices=None, simple=False):
        """
        Return an iterable containing the graphic extensions.

        This method iterates over the vertices in the input. If
        ``simple == False``, it first extends by a loop. It will then add an
        edge between every pair of vertices in the input, skipping pairs of
        vertices with an edge already between them if ``simple == True``.

        This method only considers the current graph presentation, and
        does not take 2-isomorphism into account. Use
        :meth:`twist <sage.matroids.graphic_matroid.GraphicMatroid.twist>` or
        :meth:`one_sum <sage.matroids.graphic_matroid.GraphicMatroid.one_sum>`
        if you wish to change the graph presentation.

        INPUT:

        - ``element`` -- (optional) the name of the newly added element in
          each extension
        - ``vertices`` -- (optional) a set of vertices over which the extension
          may be taken
        - ``simple`` -- boolean (default: ``False``); if ``True``, extensions
          by loops and parallel elements are not taken

        OUTPUT:

        An iterable containing instances of :class:`GraphicMatroid`. If
        ``vertices`` is not specified, every vertex is used.

        .. NOTE::

            The extension by a loop will always occur unless
            ``simple == True``. The extension by a coloop will never occur.

        EXAMPLES::

            sage: M = Matroid(range(5), graphs.DiamondGraph())
            sage: I = M.graphic_extensions('a')
            sage: for N in I:
            ....:     list(N.graph().edge_iterator())
            [(0, 0, 'a'), (0, 1, 0), (0, 2, 1), (1, 2, 2), (1, 3, 3), (2, 3, 4)]
            [(0, 1, 'a'), (0, 1, 0), (0, 2, 1), (1, 2, 2), (1, 3, 3), (2, 3, 4)]
            [(0, 1, 0), (0, 2, 'a'), (0, 2, 1), (1, 2, 2), (1, 3, 3), (2, 3, 4)]
            [(0, 1, 0), (0, 2, 1), (0, 3, 'a'), (1, 2, 2), (1, 3, 3), (2, 3, 4)]
            [(0, 1, 0), (0, 2, 1), (1, 2, 'a'), (1, 2, 2), (1, 3, 3), (2, 3, 4)]
            [(0, 1, 0), (0, 2, 1), (1, 2, 2), (1, 3, 'a'), (1, 3, 3), (2, 3, 4)]
            [(0, 1, 0), (0, 2, 1), (1, 2, 2), (1, 3, 3), (2, 3, 'a'), (2, 3, 4)]

        ::

            sage: M = Matroid(graphs.CompleteBipartiteGraph(3,3))
            sage: I = M.graphic_extensions(simple=True)
            sage: sum (1 for i in I)
            6
            sage: I = M.graphic_extensions(vertices=[0,1,2])
            sage: sum (1 for i in I)
            4
        """
        G = self.graph()
        if element is None:
            element = newlabel(self.groundset())
        elif element in self.groundset():
            raise ValueError("cannot extend by element already in groundset")
        if vertices is None:
            vertices = self._G.vertices(sort=False)
        elif not set(vertices).issubset(self._G.vertices(sort=False)):
            raise ValueError("vertices are not all in the graph")

        # First extend by a loop, then consider every pair of vertices.
        # Put the loop on the first vertex.
        if not simple:
            G.add_edge(vertices[0], vertices[0], element)
            yield GraphicMatroid(G)
            G.delete_edge(vertices[0], vertices[0], element)

        pairs = combinations(vertices, 2)
        for p in pairs:
            if not simple or not G.has_edge(p[0], p[1]):
                G.add_edge(p[0], p[1], element)
                yield GraphicMatroid(G)
                G.delete_edge(p[0], p[1], element)

    cpdef graphic_coextension(self, u, v=None, X=None, element=None):
        """
        Return a matroid coextended by a new element.

        A coextension in a graphic matroid is the opposite of contracting an
        edge; that is, a vertex is split, and a new edge is added between the
        resulting vertices. This method will create a new vertex `v` adjacent
        to `u`, and move the edges indicated by `X` from `u` to `v`.

        INPUT:

        - ``u`` -- the vertex to be split
        - ``v`` -- (optional) the name of the new vertex after splitting
        - ``X`` -- (optional) a list of the matroid elements corresponding to
          edges incident to ``u`` that move to the new vertex after splitting
        - ``element`` -- (optional) the name of the newly added element

        OUTPUT:

        An instance of :class:`GraphicMatroid` coextended by the new element.
        If ``X`` is not specified, the new element will be a coloop.

        .. NOTE::

            A loop on ``u`` will stay a loop unless it is in ``X``.

        EXAMPLES::

            sage: G = Graph([(0, 1, 0), (0, 2, 1), (0, 3, 2), (0, 4, 3),
            ....:            (1, 2, 4), (1, 4, 5), (2, 3, 6), (3, 4, 7)])
            sage: M = Matroid(G)
            sage: M1 = M.graphic_coextension(0, X=[1,2], element='a')
            sage: M1.graph().edges(sort=True)
            [(0, 1, 0),
             (0, 4, 3),
             (0, 5, 'a'),
             (1, 2, 4),
             (1, 4, 5),
             (2, 3, 6),
             (2, 5, 1),
             (3, 4, 7),
             (3, 5, 2)]

        TESTS::

            sage: M = Matroid(range(3), graphs.CycleGraph(3))
            sage: M = M.graphic_extension(0, element='a')
            sage: M.graph().edges(sort=True)
            [(0, 0, 'a'), (0, 1, 0), (0, 2, 1), (1, 2, 2)]
            sage: M1 = M.graphic_coextension(0, X=[1], element='b')
            sage: M1.graph().edges(sort=True)
            [(0, 0, 'a'), (0, 1, 0), (0, 3, 'b'), (1, 2, 2), (2, 3, 1)]
            sage: M2 = M.graphic_coextension(0, X=[1, 'a'], element='b')
            sage: M2.graph().edges(sort=True)
            [(0, 1, 0), (0, 3, 'a'), (0, 3, 'b'), (1, 2, 2), (2, 3, 1)]

        ::

            sage: M = Matroid(graphs.CycleGraph(3))
            sage: M = M.graphic_coextension(u=2, element='a')
            sage: M.graph()
            Looped multi-graph on 4 vertices
            sage: M.graph().loops()
            []
            sage: M = M.graphic_coextension(u=2, element='a')
            Traceback (most recent call last):
            ...
            ValueError: cannot extend by element already in groundset
            sage: M = M.graphic_coextension(u=4)
            Traceback (most recent call last):
            ...
            ValueError: u must be an existing vertex

        TESTS::

            sage: M = Matroid(graphs.EmptyGraph())
            sage: M.graphic_coextension(u=0)
            Graphic matroid of rank 1 on 1 elements

            sage: M = Matroid(graphs.DiamondGraph())
            sage: N = M.graphic_coextension(0,'q')
            sage: list(N.graph().vertex_iterator())
            ['q', 0, 1, 2, 3]

        ::

            sage: M = Matroid(range(5), graphs.DiamondGraph())
            sage: N = M.graphic_coextension(u=3, v=5, element='a')
            sage: N.graph().edges(sort=True)
            [(0, 1, 0), (0, 2, 1), (1, 2, 2), (1, 3, 3), (2, 3, 4), (3, 5, 'a')]
            sage: N = M.graphic_coextension(u=3, element='a')
            sage: N.graph().edges(sort=True)
            [(0, 1, 0), (0, 2, 1), (1, 2, 2), (1, 3, 3), (2, 3, 4), (3, 4, 'a')]
            sage: N = M.graphic_coextension(u=3, v=3, element='a')
            Traceback (most recent call last):
            ...
            ValueError: u and v must be distinct
        """
        if element is None:
            element = newlabel(self.groundset())
        else:
            if element in self.groundset():
                raise ValueError("cannot extend by element already in groundset")

        if u not in self._G.vertices(sort=False):
            raise ValueError("u must be an existing vertex")
        if v == u:
            raise ValueError("u and v must be distinct")
        # To prevent an error for iterating over None:
        if X is None:
            X = []

        G = self.graph()
        vertices = G.vertices(sort=False)
        if v is None:
            v = G.add_vertex()

        elif v in G:
            raise ValueError("vertex is already in the graph")
        if u not in vertices:
            G.add_edge(u, v, element)
            return GraphicMatroid(G)

        edgelist = self.groundset_to_edges(X)

        split_vertex(G, u, v, edgelist)
        G.add_edge(u, v, element)

        return GraphicMatroid(G)

    def graphic_coextensions(self, vertices=None, v=None, element=None, cosimple=False):
        """
        Return an iterator of graphic coextensions.

        This method iterates over the vertices in the input. If
        ``cosimple == False``, it first coextends by a coloop and series edge
        for every edge incident with the vertices. For vertices of degree four
        or higher, it will consider the ways to partition the vertex into two
        sets of cardinality at least two, and these will be the edges incident
        with the vertices after splitting.

        At most one series coextension will be taken for each series class.

        INPUT:

        - ``vertices`` -- (optional) the vertices to be split
        - ``v`` -- (optional) the name of the new vertex
        - ``element`` -- (optional) the name of the new element
        - ``cosimple`` -- boolean (default: ``False``); if ``True``, coextensions
          by a coloop or series elements will not be taken

        OUTPUT:

        An iterable containing instances of :class:`GraphicMatroid`. If
        ``vertices`` is not specified, the method iterates over all vertices.

        EXAMPLES::

            sage: G = Graph([(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 4), (2, 3), (3, 4)])
            sage: M = Matroid(range(8), G)
            sage: I = M.graphic_coextensions(vertices=[0], element='a')
            sage: sorted([N.graph().edges_incident(0, sort=True) for N in I], key=str)
            [[(0, 1, 0), (0, 2, 1), (0, 3, 2), (0, 4, 3), (0, 5, 'a')],
             [(0, 1, 0), (0, 2, 1), (0, 3, 2), (0, 5, 'a')],
             [(0, 1, 0), (0, 2, 1), (0, 4, 3), (0, 5, 'a')],
             [(0, 1, 0), (0, 2, 1), (0, 5, 'a')],
             [(0, 1, 0), (0, 3, 2), (0, 4, 3), (0, 5, 'a')],
             [(0, 1, 0), (0, 3, 2), (0, 5, 'a')],
             [(0, 2, 1), (0, 3, 2), (0, 4, 3), (0, 5, 'a')],
             [(0, 2, 1), (0, 3, 2), (0, 5, 'a')]]

        ::

            sage: N = Matroid(range(4), graphs.CycleGraph(4))
            sage: I = N.graphic_coextensions(element='a')
            sage: for N1 in I:  # random
            ....:     N1.graph().edges(sort=True)
            [(0, 1, 0), (0, 3, 1), (0, 4, 'a'), (1, 2, 2), (2, 3, 3)]
            [(0, 1, 0), (0, 3, 1), (1, 4, 2), (2, 3, 3), (2, 4, 'a')]
            sage: sum(1 for n in N.graphic_coextensions(cosimple=True))
            0

        TESTS::

            sage: M = Matroid(graphs.EmptyGraph())
            sage: M.graphic_coextension(0)
            Graphic matroid of rank 1 on 1 elements
            sage: I = M.graphic_coextensions(element='a')
            sage: for m in I:
            ....:     m.graph().edges(sort=True)
            [(0, 1, 'a')]
            sage: N = Matroid(graphs.CycleGraph(4))
            sage: I = N.graphic_coextensions(vertices=[3, 4], element='a')
            sage: next(I)
            Traceback (most recent call last):
            ...
            ValueError: vertices are not all in the graph

        We expect 136 graphic coextensions of an 8-spoked wheel: 128 extensions
        from the center vertex because there are (256/2) ways to put the 8
        center edges into 2 sets, and then 8 more for series extensions of the
        rims::

            sage: M = Matroid(graphs.WheelGraph(9))
            sage: I = M.graphic_coextensions()
            sage: sum(1 for N in I)
            136
            sage: I = M.graphic_coextensions(cosimple=True)
            sage: sum(1 for N in I)
            119
            sage: sum(1 for N in Matroid(graphs.WheelGraph(8)).graphic_coextensions())
            71

        This graph has max degree 3, so the only series extensions should be
        non-cosimple, ie. a coloop and one for every coseries class.
        12 total::

            sage: edgedict = {0: [1, 2, 3], 1: [2, 4], 2: [3], 3: [6],
            ....:             4: [5, 7], 5: [6, 7], 6: [7]}
            sage: M = Matroid(range(12), Graph(edgedict))
            sage: sorted(M.coclosure([4]))
            [4, 6]
            sage: sum(1 for N in M.graphic_coextensions())
            12
            sage: sum(1 for N in M.graphic_coextensions(cosimple=True))
            0
        """
        G = self.graph()
        if element is None:
            element = newlabel(self.groundset())
        elif element in self.groundset():
            raise ValueError("cannot extend by element already in groundset")
        if vertices is None:
            vertices = self._G.vertices(sort=False)
        elif not set(vertices).issubset(self._G.vertices(sort=False)):
            raise ValueError("vertices are not all in the graph")

        if v is None:
            # we just need to know what the vertex's name will be
            v = G.add_vertex()
            G.delete_vertex(v)
        elif v in G:
            raise ValueError("vertex is already in the graph")

        if not cosimple:
            # First extend by a coloop on the first vertex.
            G.add_edge(vertices[0], v, element)
            yield GraphicMatroid(G)
            G.delete_vertex(v)

            # Next add an edge in series, for every series class in the input
            edges = set(G.edges_incident(vertices))
            while edges:
                u0, v0, l = edges.pop()
                G.delete_edge(u0, v0, l)
                # place the new element on v0 if v0 is a vertex from input
                if v0 in vertices:
                    G.add_edge(u0, v, l)
                    G.add_edge(v, v0, element)
                else:
                    G.add_edge(u0, v, element)
                    G.add_edge(v, v0, l)
                yield GraphicMatroid(G)
                G.delete_vertex(v)
                G.add_edge(u0, v0, l)

                edges.difference_update(self.groundset_to_edges(self.coclosure([l])))

        # If a vertex has degree 1, or 2, or 3, we already handled it.
        for u in vertices:
            if G.degree(u) > 3:
                elts_incident = [l for (_, _, l) in G.edges_incident(u)]
                x = elts_incident.pop()
                for i in range(1, (len(elts_incident) - Integer(1))):
                    groups = combinations(elts_incident, i)
                    for g in groups:
                        g = list(g)
                        g.append(x)
                        yield self.graphic_coextension(
                            X=g, u=u, v=v, element=element)

    cpdef twist(self, X):
        """
        Perform a Whitney twist on the graph.

        `X` must be part of a 2-separation.
        The connectivity of `X` must be 1, and the subgraph induced by `X` must
        intersect the subgraph induced by the rest of the elements on exactly
        two vertices.

        INPUT:

        - ``X`` -- the set of elements to be twisted with respect
          to the rest of the matroid

        OUTPUT: :class:`GraphicMatroid` isomorphic to this matroid but
        with a graph that is not necessarily isomorphic

        EXAMPLES::

            sage: edgelist = [(0, 1, 0), (1, 2, 1), (1, 2, 2), (2, 3, 3),
            ....:             (2, 3, 4), (2, 3, 5), (3, 0, 6)]
            sage: M = Matroid(Graph(edgelist, multiedges=True))
            sage: M1 = M.twist([0, 1, 2]); M1.graph().edges(sort=True)
            [(0, 1, 1), (0, 1, 2), (0, 3, 6), (1, 2, 0), (2, 3, 3), (2, 3, 4), (2, 3, 5)]
            sage: M2 = M.twist([0, 1, 3])
            Traceback (most recent call last):
            ...
            ValueError: the input must display a 2-separation that is not a 1-separation

        TESTS::

            sage: edgedict = {0: [1, 2], 1: [2, 3], 2: [3], 3: [4, 5], 4: [5]}
            sage: M = Matroid(range(8), Graph(edgedict))
            sage: M.graph().edges(sort=True)
            [(0, 1, 0),
             (0, 2, 1),
             (1, 2, 2),
             (1, 3, 3),
             (2, 3, 4),
             (3, 4, 5),
             (3, 5, 6),
             (4, 5, 7)]
            sage: M1 = M.twist([0, 1]); M1.graph().edges(sort=True)
            [(0, 1, 1),
             (0, 2, 0),
             (1, 2, 2),
             (1, 3, 3),
             (2, 3, 4),
             (3, 4, 5),
             (3, 5, 6),
             (4, 5, 7)]
            sage: M2 = M1.twist([0, 1, 2]); M2.graph().edges(sort=True)
            [(0, 1, 0),
             (0, 2, 1),
             (1, 2, 2),
             (1, 3, 3),
             (2, 3, 4),
             (3, 4, 5),
             (3, 5, 6),
             (4, 5, 7)]
            sage: M1 == M
            False
            sage: M2 == M
            True
            sage: M2.twist([3, 4])
            Traceback (most recent call last):
            ...
            ValueError: too many vertices in the intersection
        """
        # We require two things:
        # (1) The connectivity of X is 1,
        # (2) X intersects the rest of the graph on 2 vertices
        if not set(X).issubset(self.groundset()):
            raise ValueError("X must be a subset of the groundset")
        connectivity = self.connectivity(X)
        if connectivity != 1:
            raise ValueError("the input must display a 2-separation "
                             "that is not a 1-separation")

        # Determine the vertices
        X_edges = self.groundset_to_edges(X)
        X_vertices = set([e[0] for e in X_edges]).union(
            [e[1] for e in X_edges])
        Y_edges = self.groundset_to_edges(self.groundset().difference(set(X)))
        Y_vertices = set([e[0] for e in Y_edges]).union(
            [e[1] for e in Y_edges])
        vertices = X_vertices.intersection(Y_vertices)
        if len(vertices) != 2:
            raise ValueError("too many vertices in the intersection")
        a = list(vertices)[0]
        b = list(vertices)[1]

        edges = [(u, v, l) for (u, v, l) in X_edges
                 if u in vertices or v in vertices]
        G = self.graph()
        for (u, v, l) in edges:
            G.delete_edge(u, v, l)
            if u == a:
                u = b
            elif u == b:
                u = a
            if v == a:
                v = b
            elif v == b:
                v = a
            G.add_edge(u, v, l)
        return GraphicMatroid(G)

    cpdef one_sum(self, X, u, v):
        """
        Arrange matroid components in the graph.

        The matroid's graph must be connected even if the matroid is not
        connected, but if there are multiple matroid components, the user may
        choose how they are arranged in the graph. This method will take the
        block of the graph that represents `X` and attach it by vertex `u` to
        another vertex `v` in the graph.

        INPUT:

        - ``X`` -- subset of the groundset
        - ``u`` -- vertex spanned by the edges of the elements in ``X``
        - ``v`` -- vertex spanned by the edges of the elements not in ``X``

        OUTPUT: :class:`GraphicMatroid` isomorphic to this matroid but
        with a graph that is not necessarily isomorphic

        EXAMPLES::

            sage: edgedict = {0:[1, 2], 1:[2, 3], 2:[3], 3:[4, 5], 6:[4, 5]}
            sage: M = Matroid(range(9), Graph(edgedict))
            sage: M.graph().edges(sort=True)
            [(0, 1, 0),
             (0, 2, 1),
             (1, 2, 2),
             (1, 3, 3),
             (2, 3, 4),
             (3, 4, 5),
             (3, 5, 6),
             (4, 6, 7),
             (5, 6, 8)]
            sage: M1 = M.one_sum(u=3, v=1, X=[5, 6, 7, 8])
            sage: M1.graph().edges(sort=True)
            [(0, 1, 0),
             (0, 2, 1),
             (1, 2, 2),
             (1, 3, 3),
             (1, 4, 5),
             (1, 5, 6),
             (2, 3, 4),
             (4, 6, 7),
             (5, 6, 8)]
            sage: M2 = M.one_sum(u=4, v=3, X=[5, 6, 7, 8])
            sage: M2.graph().edges(sort=True)
            [(0, 1, 0),
             (0, 2, 1),
             (1, 2, 2),
             (1, 3, 3),
             (2, 3, 4),
             (3, 6, 7),
             (3, 7, 5),
             (5, 6, 8),
             (5, 7, 6)]

        TESTS::

            sage: M = matroids.CompleteGraphic(4)
            sage: M.one_sum(u=1, v=2, X=[0,1])
            Traceback (most recent call last):
            ...
            ValueError: the input must display a 1-separation

        ::

            sage: M = Matroid(range(5), graphs.BullGraph())
            sage: M.graph().edges(sort=True)
            [(0, 1, 0), (0, 2, 1), (1, 2, 2), (1, 3, 3), (2, 4, 4)]
            sage: M1 = M.one_sum(u=3, v=0, X=[3,4])
            Traceback (most recent call last):
            ...
            ValueError: too many vertices in the intersection

            sage: M1 = M.one_sum(u=3, v=2, X=[3])
            sage: M1.graph().edges(sort=True)
            [(0, 1, 0), (0, 2, 1), (1, 2, 2), (2, 4, 4), (2, 5, 3)]

            sage: M2 = M1.one_sum(u=5, v=0, X=[3,4])
            sage: M2.graph().edges(sort=True)
            [(0, 1, 0), (0, 2, 1), (0, 3, 3), (1, 2, 2), (3, 4, 4)]

            sage: M = Matroid(range(5), graphs.BullGraph())
            sage: M.one_sum(u=0, v=1, X=[3])
            Traceback (most recent call last):
            ...
            ValueError: first vertex must be spanned by the input

            sage: M.one_sum(u=1, v=3, X=[3])
            Traceback (most recent call last):
            ...
            ValueError: second vertex must be spanned by the rest of the graph
        """
        # We require two things:
        # (1) The connectivity of X is 0,
        # (2) X intersects the rest of the graph on 1 vertex
        if not set(X).issubset(self.groundset()):
            raise ValueError("X must be a subset of the groundset")
        connectivity = self.connectivity(X)
        if connectivity != 0:
            raise ValueError("the input must display a 1-separation")
        G = self.graph()
        if u not in G or v not in G:
            raise ValueError("the vertices must already be in the graph")

        # Determine the vertex
        X_edges = self.groundset_to_edges(X)
        X_vertices = set([e[0] for e in X_edges]).union(
            [e[1] for e in X_edges])
        if u not in X_vertices:
            raise ValueError("first vertex must be spanned by the input")
        Y_edges = self.groundset_to_edges(self.groundset().difference(set(X)))
        Y_vertices = set([e[0] for e in Y_edges]).union(
            [e[1] for e in Y_edges])
        if v not in Y_vertices:
            raise ValueError("second vertex must be spanned by the rest of the graph")
        vertices = X_vertices.intersection(Y_vertices)
        if len(vertices) != 1:
            raise ValueError("too many vertices in the intersection")
        a = vertices.pop()
        b = G.add_vertex()

        edgeset = set(X_edges).intersection(set(G.edges_incident(a)))
        split_vertex(G, a, b, edgeset)
        # If u was the cut vertex, u is now detached from our component
        # so we merge the new vertex. Otherwise we can merge u
        if u == a:
            G.merge_vertices([v, b])
        else:
            G.merge_vertices([v, u])

        return GraphicMatroid(G)

    cpdef regular_matroid(self):
        """
        Return an instance of :class:`RegularMatroid` isomorphic to this
        :class:`GraphicMatroid`.

        EXAMPLES::

            sage: M = matroids.CompleteGraphic(5); M
            M(K5): Graphic matroid of rank 4 on 10 elements
            sage: N = M.regular_matroid(); N
            Regular matroid of rank 4 on 10 elements with 125 bases
            sage: M.equals(N)
            True
            sage: M == N
            False

        TESTS:

        Check that :issue:`28482` is fixed::

            sage: G = Graph([[3, 4], [4, 1], [1, 2], [2, 3], [3, 5], [5, 6], [6, 3]])
            sage: M = Matroid(G)
            sage: R = M.regular_matroid()
            sage: set(M.circuits()) == set(R.circuits())
            True
        """
        from sage.matroids.constructor import Matroid as ConstructorMatroid
        X = [l for u, v, l in self._G.edge_iterator()]
        return ConstructorMatroid(groundset=X, graph=self._G, regular=True)

    cpdef relabel(self, mapping):
        r"""
        Return an isomorphic matroid with relabeled groundset.

        The output is obtained by relabeling each element `e` by
        ``mapping[e]``, where ``mapping`` is a given injective map. If
        ``mapping[e]`` is not defined, then the identity map is assumed.

        INPUT:

        - ``mapping`` -- a Python object such that ``mapping[e]`` is the new
          label of `e`

        OUTPUT: matroid

        EXAMPLES::

            sage: M = matroids.CompleteGraphic(4)
            sage: sorted(M.groundset())
            [0, 1, 2, 3, 4, 5]
            sage: N = M.relabel({0: 6, 5: 'e'})
            sage: sorted(N.groundset(), key=str)
            [1, 2, 3, 4, 6, 'e']
            sage: N.is_isomorphic(M)
            True

        TESTS::

            sage: M = matroids.CompleteGraphic(4)
            sage: f = {0: 'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e', 5: 'f', 6: 'g'}
            sage: N = M.relabel(f)
            sage: for S in powerset(M.groundset()):
            ....:     assert M.rank(S) == N.rank([f[x] for x in S])
        """
        d = self._relabel_map(mapping)
        E = [d[x] for x in self.groundset()]
        M = GraphicMatroid(self.graph(), groundset=E)
        return M
