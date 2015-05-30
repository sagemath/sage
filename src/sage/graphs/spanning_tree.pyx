r"""
Spanning trees

This module is a collection of algorithms on spanning trees. Also included in
the collection are algorithms for minimum spanning trees. See the book
[JoynerNguyenCohen2010]_ for descriptions of spanning tree algorithms,
including minimum spanning trees.

.. SEEALSO::

   * :meth:`GenericGraph.min_spanning_tree
     <sage.graphs.generic_graph.GenericGraph.min_spanning_tree>`.

**Todo**

* Rewrite :func:`kruskal` to use priority queues. Once Cython has support
  for generators and the ``yield`` statement, rewrite :func:`kruskal` to use
  ``yield``.
* Prim's algorithm.
* Boruvka's algorithm.
* Parallel version of Boruvka's algorithm.
* Randomized spanning tree construction.

REFERENCES:

.. [Aldous90] D. Aldous, 'The random walk construction of
  uniform spanning trees', SIAM J Discrete Math 3 (1990),
  450-465.

.. [Broder89] A. Broder, 'Generating random spanning trees',
  Proceedings of the 30th IEEE Symposium on Foundations of
  Computer Science, 1989, pp. 442-447. :doi:`10.1109/SFCS.1989.63516`,
  <http://www.cs.cmu.edu/~15859n/RelatedWork/Broder-GenRanSpanningTrees.pdf>_

.. [CormenEtAl2001] Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest,
  and Clifford Stein. *Introduction to Algorithms*. 2nd edition, The MIT Press,
  2001.

.. [GoodrichTamassia2001] Michael T. Goodrich and Roberto Tamassia.
  *Data Structures and Algorithms in Java*. 2nd edition, John Wiley & Sons,
  2001.

.. [JoynerNguyenCohen2010] David Joyner, Minh Van Nguyen, and Nathann Cohen.
  *Algorithmic Graph Theory*. 2010,
  http://code.google.com/p/graph-theory-algorithms-book/

.. [Sahni2000] Sartaj Sahni. *Data Structures, Algorithms, and Applications
  in Java*. McGraw-Hill, 2000.


Methods
-------
"""

###########################################################################
# Copyright (c) 2007 Jason Grout <jason-sage@creativetrax.com>
# Copyright (c) 2009 Mike Hansen <mhansen@gmail.com>
# Copyright (c) 2010 Gregory McWhirter <gmcwhirt@uci.edu>
# Copyright (c) 2010 Minh Van Nguyen <nguyenminh2@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# http://www.gnu.org/licenses/
###########################################################################

include "sage/ext/cdefs.pxi"
include "sage/ext/interrupt.pxi"

cpdef kruskal(G, wfunction=None, bint check=False):
    r"""
    Minimum spanning tree using Kruskal's algorithm.

    This function assumes that we can only compute minimum spanning trees for
    undirected simple graphs. Such graphs can be weighted or unweighted.

    INPUT:

    - ``G`` -- A graph. This can be an undirected graph, a digraph, a
      multigraph, or a multidigraph. Note the following behaviours:

      - If ``G`` is unweighted, then consider the simple version of ``G``
        with all self-loops and multiple edges removed.

      - If ``G`` is directed, then we only consider its undirected version.

      - If ``G`` is weighted, we ignore all of its self-loops. Note that a
        weighted graph should only have numeric weights. You cannot assign
        numeric weights to some edges of ``G``, but have ``None`` as a
        weight for some other edge. If your input graph is weighted, you are
        responsible for assign numeric weight to each of its edges.
        Furthermore, we remove multiple edges as follows. First we convert
        ``G`` to be undirected. Suppose there are multiple edges from `u` to
        `v`. Among all such multiple edges, we choose one with minimum weight.

    - ``wfunction`` -- A weight function: a function that takes an edge and
      returns a numeric weight. Default: ``None``. The default is to
      assign each edge a weight of 1.

    - ``check`` -- Whether to first perform sanity checks on the input
      graph ``G``. Default: ``check=False``. If we toggle ``check=True``, the
      following sanity checks are first performed on ``G`` prior to running
      Kruskal's algorithm on that input graph:

      - Is ``G`` the null graph?
      - Is ``G`` disconnected?
      - Is ``G`` a tree?
      - Is ``G`` directed?
      - Does ``G`` have self-loops?
      - Does ``G`` have multiple edges?
      - Is ``G`` weighted?

      By default, we turn off the sanity checks for performance reasons. This
      means that by default the function assumes that its input graph is
      simple, connected, is not a tree, and has at least one vertex.
      If the input graph does not satisfy all of the latter conditions, you
      should set ``check=True`` to perform some sanity checks and
      preprocessing on the input graph. To further improve the runtime of this
      function, you should call it directly instead of using it indirectly
      via :meth:`sage.graphs.generic_graph.GenericGraph.min_spanning_tree`.

    OUTPUT:

    The edges of a minimum spanning tree of ``G``, if one exists, otherwise
    returns the empty list.

    - If ``G`` is a tree, return the edges of ``G`` regardless of whether
      ``G`` is weighted or unweighted, directed or undirected.

    - If ``G`` is unweighted, default to using unit weight for each edge of
      ``G``. The default behaviour is to use the already assigned weights of
      ``G`` provided that ``G`` is weighted.

    - If ``G`` is weighted and a weight function is also supplied, then use
      the already assigned weights of ``G``, not the weight function. If you
      really want to use a weight function for ``G`` even if ``G`` is
      weighted, first convert ``G`` to be unweighted and pass in the weight
      function.

    .. seealso::

        - :meth:`sage.graphs.generic_graph.GenericGraph.min_spanning_tree`

    EXAMPLES:

    An example from pages 727--728 in [Sahni2000]_. ::

        sage: from sage.graphs.spanning_tree import kruskal
        sage: G = Graph({1:{2:28, 6:10}, 2:{3:16, 7:14}, 3:{4:12}, 4:{5:22, 7:18}, 5:{6:25, 7:24}})
        sage: G.weighted(True)
        sage: E = kruskal(G, check=True); E
        [(1, 6, 10), (3, 4, 12), (2, 7, 14), (2, 3, 16), (4, 5, 22), (5, 6, 25)]

    Variants of the previous example. ::

        sage: H = Graph(G.edges(labels=False))
        sage: kruskal(H, check=True)
        [(1, 2, None), (1, 6, None), (2, 3, None), (2, 7, None), (3, 4, None), (4, 5, None)]
        sage: H = DiGraph(G.edges(labels=False))
        sage: kruskal(H, check=True)
        [(1, 2, None), (1, 6, None), (2, 3, None), (2, 7, None), (3, 4, None), (4, 5, None)]
        sage: G.allow_loops(True)
        sage: G.allow_multiple_edges(True)
        sage: G
        Looped multi-graph on 7 vertices
        sage: for i in range(20):
        ...       u = randint(1, 7)
        ...       v = randint(1, 7)
        ...       w = randint(0, 20)
        ...       G.add_edge(u, v, w)
        sage: H = copy(G)
        sage: H
        Looped multi-graph on 7 vertices
        sage: def sanitize(G):
        ...       G.allow_loops(False)
        ...       E = {}
        ...       for u, v, _ in G.multiple_edges():
        ...           E.setdefault(u, v)
        ...       for u in E:
        ...           W = sorted(G.edge_label(u, E[u]))
        ...           for w in W[1:]:
        ...               G.delete_edge(u, E[u], w)
        ...       G.allow_multiple_edges(False)
        sage: sanitize(H)
        sage: H
        Graph on 7 vertices
        sage: kruskal(G, check=True) == kruskal(H, check=True)
        True

    Note that we only consider an undirected version of the input graph. Thus
    if ``G`` is a weighted multidigraph and ``H`` is an undirected version of
    ``G``, then this function should return the same minimum spanning tree
    for both ``G`` and ``H``. ::

        sage: from sage.graphs.spanning_tree import kruskal
        sage: G = DiGraph({1:{2:[1,14,28], 6:[10]}, 2:{3:[16], 1:[15], 7:[14], 5:[20,21]}, 3:{4:[12,11]}, 4:{3:[13,3], 5:[22], 7:[18]}, 5:{6:[25], 7:[24], 2:[1,3]}}, multiedges=True)
        sage: G.multiple_edges(to_undirected=False)
        [(1, 2, 1), (1, 2, 14), (1, 2, 28), (5, 2, 1), (5, 2, 3), (4, 3, 3), (4, 3, 13), (3, 4, 11), (3, 4, 12), (2, 5, 20), (2, 5, 21)]
        sage: H = G.to_undirected()
        sage: H.multiple_edges(to_undirected=True)
        [(1, 2, 1), (1, 2, 14), (1, 2, 15), (1, 2, 28), (2, 5, 1), (2, 5, 3), (2, 5, 20), (2, 5, 21), (3, 4, 3), (3, 4, 11), (3, 4, 12), (3, 4, 13)]
        sage: kruskal(G, check=True)
        [(1, 2, 1), (1, 6, 10), (2, 3, 16), (2, 5, 1), (2, 7, 14), (3, 4, 3)]
        sage: kruskal(G, check=True) == kruskal(H, check=True)
        True
        sage: G.weighted(True)
        sage: H.weighted(True)
        sage: kruskal(G, check=True)
        [(1, 2, 1), (2, 5, 1), (3, 4, 3), (1, 6, 10), (2, 7, 14), (2, 3, 16)]
        sage: kruskal(G, check=True) == kruskal(H, check=True)
        True

    An example from pages 599--601 in [GoodrichTamassia2001]_. ::

        sage: G = Graph({"SFO":{"BOS":2704, "ORD":1846, "DFW":1464, "LAX":337},
        ...   "BOS":{"ORD":867, "JFK":187, "MIA":1258},
        ...   "ORD":{"PVD":849, "JFK":740, "BWI":621, "DFW":802},
        ...   "DFW":{"JFK":1391, "MIA":1121, "LAX":1235},
        ...   "LAX":{"MIA":2342},
        ...   "PVD":{"JFK":144},
        ...   "JFK":{"MIA":1090, "BWI":184},
        ...   "BWI":{"MIA":946}})
        sage: G.weighted(True)
        sage: kruskal(G, check=True)
        [('JFK', 'PVD', 144), ('BWI', 'JFK', 184), ('BOS', 'JFK', 187), ('LAX', 'SFO', 337), ('BWI', 'ORD', 621), ('DFW', 'ORD', 802), ('BWI', 'MIA', 946), ('DFW', 'LAX', 1235)]

    An example from pages 568--569 in [CormenEtAl2001]_. ::

        sage: G = Graph({"a":{"b":4, "h":8}, "b":{"c":8, "h":11},
        ...   "c":{"d":7, "f":4, "i":2}, "d":{"e":9, "f":14},
        ...   "e":{"f":10}, "f":{"g":2}, "g":{"h":1, "i":6}, "h":{"i":7}})
        sage: G.weighted(True)
        sage: kruskal(G, check=True)
        [('g', 'h', 1), ('c', 'i', 2), ('f', 'g', 2), ('a', 'b', 4), ('c', 'f', 4), ('c', 'd', 7), ('a', 'h', 8), ('d', 'e', 9)]

    TESTS:

    The input graph must not be empty. ::

        sage: from sage.graphs.spanning_tree import kruskal
        sage: kruskal(graphs.EmptyGraph(), check=True)
        []
        sage: kruskal(Graph(), check=True)
        []
        sage: kruskal(Graph(multiedges=True), check=True)
        []
        sage: kruskal(Graph(loops=True), check=True)
        []
        sage: kruskal(Graph(multiedges=True, loops=True), check=True)
        []
        sage: kruskal(DiGraph(), check=True)
        []
        sage: kruskal(DiGraph(multiedges=True), check=True)
        []
        sage: kruskal(DiGraph(loops=True), check=True)
        []
        sage: kruskal(DiGraph(multiedges=True, loops=True), check=True)
        []

    The input graph must be connected. ::

        sage: def my_disconnected_graph(n, ntries, directed=False, multiedges=False, loops=False):
        ...       G = Graph()
        ...       k = randint(1, n)
        ...       G.add_vertices(range(k))
        ...       if directed:
        ...           G = G.to_directed()
        ...       if multiedges:
        ...           G.allow_multiple_edges(True)
        ...       if loops:
        ...           G.allow_loops(True)
        ...       for i in range(ntries):
        ...           u = randint(0, k-1)
        ...           v = randint(0, k-1)
        ...           G.add_edge(u, v)
        ...       while G.is_connected():
        ...           u = randint(0, k-1)
        ...           v = randint(0, k-1)
        ...           G.delete_edge(u, v)
        ...       return G
        sage: G = my_disconnected_graph(100, 50, directed=False, multiedges=False, loops=False)  # long time
        sage: kruskal(G, check=True)  # long time
        []
        sage: G = my_disconnected_graph(100, 50, directed=False, multiedges=True, loops=False)  # long time
        sage: kruskal(G, check=True)  # long time
        []
        sage: G = my_disconnected_graph(100, 50, directed=False, multiedges=True, loops=True)  # long time
        sage: kruskal(G, check=True)  # long time
        []
        sage: G = my_disconnected_graph(100, 50, directed=True, multiedges=False, loops=False)  # long time
        sage: kruskal(G, check=True)  # long time
        []
        sage: G = my_disconnected_graph(100, 50, directed=True, multiedges=True, loops=False)  # long time
        sage: kruskal(G, check=True)  # long time
        []
        sage: G = my_disconnected_graph(100, 50, directed=True, multiedges=True, loops=True)  # long time
        sage: kruskal(G, check=True)  # long time
        []

    If the input graph is a tree, then return its edges. ::

        sage: T = graphs.RandomTree(randint(1, 50))  # long time
        sage: T.edges() == kruskal(T, check=True)  # long time
        True
    """
    g = G
    sortedE_iter = None
    # sanity checks
    if check:
        if G.order() == 0:
            return []
        if not G.is_connected():
            return []
        # G is now assumed to be a nonempty connected graph
        if G.num_verts() == G.num_edges() + 1:
            # G is a tree
            return G.edges()
        g = G.to_undirected()
        g.allow_loops(False)
        if g.weighted():
            # If there are multiple edges from u to v, retain the edge of
            # minimum weight among all such edges.
            if g.allows_multiple_edges():
                # By this stage, g is assumed to be an undirected, weighted
                # multigraph. Thus when we talk about a weighted multiedge
                # (u, v, w) of g, we mean that (u, v, w) and (v, u, w) are
                # one and the same undirected multiedge having the same weight
                # w.
                # If there are multiple edges from u to v, retain only the
                # start and end vertices of such edges. Let a and b be the
                # start and end vertices, respectively, of a weighted edge
                # (a, b, w) having weight w. Then there are multiple weighted
                # edges from a to b if and only if the set uniqE has the
                # tuple (a, b) as an element.
                uniqE = set()
                for u, v, _ in iter(g.multiple_edges(to_undirected=True)):
                    uniqE.add((u, v))
                # Let (u, v) be an element in uniqE. Then there are multiple
                # weighted edges from u to v. Let W be a list of all edge
                # weights of multiple edges from u to v, sorted in
                # nondecreasing order. If w is the first element in W, then
                # (u, v, w) is an edge of minimum weight (there may be
                # several edges of minimum weight) among all weighted edges
                # from u to v. If i >= 2 is the i-th element in W, delete the
                # multiple weighted edge (u, v, i).
                for u, v in uniqE:
                    W = sorted(g.edge_label(u, v))
                    for w in W[1:]:
                        g.delete_edge(u, v, w)
                # all multiple edges should now be removed; check this!
                assert g.multiple_edges() == []
                g.allow_multiple_edges(False)
            # sort edges by weights
            from operator import itemgetter
            sortedE_iter = iter(sorted(g.edges(), key=itemgetter(2)))
        else:
            g = g.to_simple()
            if wfunction is None:
                sortedE_iter = iter(sorted(g.edges()))
            else:
                sortedE_iter = iter(sorted(g.edges(), key=wfunction))
    # G is assumed to be simple, undirected, and unweighted
    else:
        if wfunction is None:
            sortedE_iter = iter(sorted(g.edges()))
        else:
            sortedE_iter = iter(sorted(g.edges(), key=wfunction))
    # Kruskal's algorithm
    T = []
    cdef int n = g.order()
    cdef int m = n - 1
    cdef int i = 0  # count the number of edges added so far
    union_find = dict()
    while i < m:
        e = sortedE_iter.next()
        components = []
        # acyclic test via union-find
        for startv in iter(e[0:2]):
            v = startv
            children = []
            # find the component a vertex lives in
            while v in union_find:
                children.append(v)
                v = union_find[v]
            # compress the paths as much as we can for efficiency reasons
            for c in children:
                union_find[c] = v
            components.append(v)
        if components[0] != components[1]:
            i += 1
            # NOTE: Once Cython supports generator and the yield statement,
            # we should replace the following line with a yield statement.
            # That way, we could access the edge of a minimum spanning tree
            # immediately after it is found, instead of waiting for all the
            # edges to be found and return the edges as a list.
            T.append(e)
            # union the components by making one the parent of the other
            union_find[components[0]] = components[1]
    return T


def random_spanning_tree(self, output_as_graph=False):
    r"""
    Return a random spanning tree of the graph.

    This uses the Aldous-Broder algorithm ([Broder89]_, [Aldous90]_)
    to generate a random spanning tree with the uniform distribution,
    as follows.

    Start from any vertex. Perform a random walk by choosing at every
    step one neighbor uniformly at random. Every time a new vertex `j`
    is met, add the edge `(i, j)` to the spanning tree, where `i` is
    the previous vertex in the random walk.

    INPUT:

    - ``output_as_graph`` -- boolean (default: ``False``) whether to return a
      list of edges or a graph.

    .. SEEALSO::

        :meth:`~sage.graphs.generic_graph.GenericGraph.spanning_trees_count`
        and :meth:`~sage.graphs.graph.Graph.spanning_trees`

    EXAMPLES::

        sage: G = graphs.TietzeGraph()
        sage: G.random_spanning_tree(output_as_graph=True)
        Graph on 12 vertices
        sage: rg = G.random_spanning_tree(); rg # random
        [(0, 9),
        (9, 11),
        (0, 8),
        (8, 7),
        (7, 6),
        (7, 2),
        (2, 1),
        (1, 5),
        (9, 10),
        (5, 4),
        (2, 3)]
        sage: Graph(rg).is_tree()
        True

    A visual example for the grid graph::

        sage: G = graphs.Grid2dGraph(6, 6)
        sage: pos = G.get_pos()
        sage: T = G.random_spanning_tree(True)
        sage: T.set_pos(pos)
        sage: T.show(vertex_labels=False)

    TESTS::

        sage: G = Graph()
        sage: G.random_spanning_tree()
        Traceback (most recent call last):
        ...
        ValueError: works only for non-empty connected graphs

        sage: G = graphs.CompleteGraph(3).complement()
        sage: G.random_spanning_tree()
        Traceback (most recent call last):
        ...
        ValueError: works only for non-empty connected graphs
    """
    from sage.misc.prandom import randint
    from sage.graphs.graph import Graph

    cdef int N = self.order()

    if N == 0 or not self.is_connected():
        raise ValueError('works only for non-empty connected graphs')

    s = self.vertex_iterator().next()
    found = set([s])
    cdef int found_nr = 1
    tree_edges = []
    while found_nr < N:
        neighbours = self.neighbors(s)
        new_s = neighbours[randint(0, len(neighbours) - 1)]
        if not(new_s in found):
            found.add(new_s)
            found_nr += 1
            tree_edges += [(s, new_s)]
        s = new_s

    if not output_as_graph:
        return tree_edges
    return Graph(tree_edges)
