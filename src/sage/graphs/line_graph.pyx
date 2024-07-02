# cython: binding=True
r"""
Line graphs

This module gather everything which is related to line graphs. Right now, this
amounts to the following functions :

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`line_graph` | Return the line graph of a given graph
    :meth:`is_line_graph` | Check whether a graph is a line graph
    :meth:`root_graph` | Return the root graph corresponding to the given graph

Author:

- Nathann Cohen (01-2013), :meth:`root_graph` method and module documentation.
  Written while listening to Nina Simone *"I wish I knew how it would feel to be
  free"*. Crazy good song. And *"Prendre ta douleur"*, too.

- David Coudert (10-2018), use maximal cliques iterator in :meth:`root_graph`,
  and use :meth:`root_graph` instead of forbidden subgraph search in
  :meth:`is_line_graph` (:issue:`26444`).

Definition
-----------

Given a graph `G`, the *line graph* `L(G)` of `G` is the graph such that

.. MATH::

    V(L(G)) =& E(G)\\
    E(L(G)) =& \{(e,e'):\text{ and }e,e'\text{ have a common endpoint in }G\}\\

The definition is extended to directed graphs. In this situation, there is an
arc `(e,e')` in `L(G)` if the destination of `e` is the origin of `e'`.

For more information, see the :wikipedia:`Line_graph`.

Root graph
----------

A graph whose line graph is `LG` is called the *root graph* of `LG`. The root
graph of a (connected) graph is unique ([Whi1932]_, [Har1969]_), except when
`LG=K_3`, as both `L(K_3)` and `L(K_{1,3})` are equal to `K_3`.

Here is how we can *"see"* `G` by staring (very intently) at `LG` :

  A graph `LG` is the line graph of `G` if there exists a collection
  `(S_v)_{v\in G}` of subsets of `V(LG)` such that :

  * Every `S_v` is a complete subgraph of `LG`.

  * Every `v\in LG` belongs to exactly two sets of the family `(S_v)_{v\in G}`.

  * Any two sets of `(S_v)_{v\in G}` have at most one common elements

  * For any edge `(u,v)\in LG` there exists a set of `(S_v)_{v\in G}` containing
    both `u` and `v`.

  In this family, each set `S_v` represent a vertex of `G`, and contains "the
  set of edges incident to `v` in `G`". Two elements `S_v,S_{v'}` have a
  nonempty intersection whenever `vv'` is an edge of `G`.

  Hence, finding the root graph of `LG` is the job of finding this collection of
  sets.

In particular, what we know for sure is that a maximal clique `S` of size `2` or
`\geq 4` in `LG` corresponds to a vertex of degree `|S|` in `G`, whose incident
edges are the elements of `S` itself.

The main problem lies with maximal cliques of size 3, i.e. triangles. Those we
have to split into two categories, *even* and *odd* triangles :

  A triangle `\{e_1,e_2,e_3\}\subseteq V(LG)` is said to be an *odd* triangle if
  there exists a vertex `e\in V(G)` incident to exactly *one* or *all* of
  `\{e_1,e_2,e_3\}`, and it is said to be *even* otherwise.

The very good point of this definition is that an inclusionwise maximal clique
which is an odd triangle will always correspond to a vertex of degree 3 in `G`,
while an even triangle could result from either a vertex of degree 3 in `G` or a
triangle in `G`. And in order to build the root graph we obviously have to
decide *which*.

Beineke proves in [Bei1970]_ that the collection of sets we are looking for
can be easily found. Indeed it turns out that it is the union of :

#. The family of all maximal cliques of `LG` of size 2 or `\geq 4`, as well as
   all odd triangles.

#. The family of all pairs of adjacent vertices which appear in exactly *one*
   maximal clique which is an even triangle.

There are actually four special cases to which the decomposition above does not
apply, i.e. graphs containing an edge which belongs to exactly two even
triangles. We deal with those independently.

* The :meth:`Complete graph
  <sage.graphs.graph_generators.GraphGenerators.CompleteGraph>` `K_3`.

* The :meth:`Diamond graph
  <sage.graphs.graph_generators.GraphGenerators.DiamondGraph>` -- the line graph
  of `K_{1,3}` plus an edge.

* The :meth:`Wheel graph
  <sage.graphs.graph_generators.GraphGenerators.WheelGraph>` on `4+1` vertices
  -- the line graph of the :meth:`Diamond graph
  <sage.graphs.graph_generators.GraphGenerators.DiamondGraph>`

* The :meth:`Octahedron
  <sage.graphs.graph_generators.GraphGenerators.OctahedralGraph>` -- the line
  graph of `K_4`.

This decomposition turns out to be very easy to implement :-)

.. WARNING::

    Even though the root graph is *NOT UNIQUE* for the triangle, this method
    returns `K_{1,3}` (and not `K_3`) in this case. Pay *very close* attention
    to that, for this answer is not theoretically correct : there is no unique
    answer in this case, and we deal with it by returning one of the two
    possible answers.

Functions
---------
"""


def is_line_graph(g, certificate=False):
    r"""
    Check whether the graph `g` is a line graph.

    INPUT:

    - ``certificate`` -- boolean; whether to return a certificate along with
      the boolean result. Here is what happens when ``certificate = True``:

      - If the graph is not a line graph, the method returns a pair ``(b,
        subgraph)`` where ``b`` is ``False`` and ``subgraph`` is a subgraph
        isomorphic to one of the 9 forbidden induced subgraphs of a line graph.

      - If the graph is a line graph, the method returns a triple ``(b,R,isom)``
        where ``b`` is ``True``, ``R`` is a graph whose line graph is the graph
        given as input, and ``isom`` is a map associating an edge of ``R`` to
        each vertex of the graph.

    .. NOTE::

        This method wastes a bit of time when the input graph is not connected.
        If you have performance in mind, it is probably better to only feed it
        with connected graphs only.

    .. SEEALSO::

        - The :mod:`line_graph <sage.graphs.line_graph>` module.

        - :meth:`~sage.graphs.graph_generators.GraphGenerators.line_graph_forbidden_subgraphs`
          -- the forbidden subgraphs of a line graph.

        - :meth:`~sage.graphs.generic_graph.GenericGraph.line_graph`

    EXAMPLES:

    A complete graph is always the line graph of a star::

        sage: graphs.CompleteGraph(5).is_line_graph()
        True

    The Petersen Graph not being claw-free, it is not a line
    graph::

        sage: graphs.PetersenGraph().is_line_graph()
        False

    This is indeed the subgraph returned::

        sage: C = graphs.PetersenGraph().is_line_graph(certificate=True)[1]             # needs sage.modules
        sage: C.is_isomorphic(graphs.ClawGraph())                                       # needs sage.modules
        True

    The house graph is a line graph::

        sage: g = graphs.HouseGraph()
        sage: g.is_line_graph()
        True

    But what is the graph whose line graph is the house ?::

        sage: # needs sage.modules
        sage: is_line, R, isom = g.is_line_graph(certificate=True)
        sage: R.sparse6_string()
        ':DaHI~'
        sage: R.show()                                                                  # needs sage.plot
        sage: isom
        {0: (0, 1), 1: (0, 2), 2: (1, 3), 3: (2, 3), 4: (3, 4)}

    TESTS:

    Disconnected graphs::

        sage: g = 2 * graphs.CycleGraph(3)
        sage: gl = g.line_graph().relabel(inplace=False)
        sage: new_g = gl.is_line_graph(certificate=True)[1]                             # needs sage.modules
        sage: g.line_graph().is_isomorphic(gl)                                          # needs sage.modules
        True

    Verify that :issue:`29740` is fixed::

        sage: g = Graph('O{e[{}^~z`MDZBZBkXzE^')
        sage: g.is_line_graph()
        False
    """
    g._scream_if_not_simple()

    def get_certificate(gg):
        """
        Assuming that gg is not a line graph, return a forbidden induced
        subgraph.
        """
        from sage.graphs.graph_generators import graphs
        for fg in graphs.line_graph_forbidden_subgraphs():
            h = gg.subgraph_search(fg, induced=True)
            if h is not None:
                return h

    if g.is_connected():
        try:
            # To check whether g is a line graph, we try to construct its root
            # graph
            R, isom = root_graph(g)
            if certificate:
                return True, R, isom
            return True
        except ValueError as VE:
            if str(VE) == "this graph is not a line graph !":
                # g is not a line graph
                if certificate:
                    return False, get_certificate(g)
                return False
            raise VE

    # g is not connected, so we apply the above procedure to each connected
    # component
    from sage.graphs.graph import Graph
    R = Graph()
    for gg in g.connected_components_subgraphs():
        try:
            RR, isom = root_graph(gg)
            R += RR
        except ValueError as VE:
            if str(VE) == "this graph is not a line graph !":
                # gg is not a line graph
                if certificate:
                    return False, get_certificate(gg)
                return False
            raise VE

    if certificate:
        _, isom = g.is_isomorphic(R.line_graph(labels=False), certificate=True)
        return True, R, isom
    return True


def line_graph(g, labels=True):
    """
    Return the line graph of the (di)graph ``g``.

    INPUT:

    - ``labels`` -- boolean (default: ``True``); whether edge labels should be
      taken in consideration. If ``labels=True``, the vertices of the line graph
      will be triples ``(u,v,label)``, and pairs of vertices otherwise.

    The line graph of an undirected graph G is an undirected graph H such that
    the vertices of H are the edges of G and two vertices e and f of H are
    adjacent if e and f share a common vertex in G. In other words, an edge in H
    represents a path of length 2 in G.

    The line graph of a directed graph G is a directed graph H such that the
    vertices of H are the edges of G and two vertices e and f of H are adjacent
    if e and f share a common vertex in G and the terminal vertex of e is the
    initial vertex of f. In other words, an edge in H represents a (directed)
    path of length 2 in G.

    .. NOTE::

        As a :class:`Graph` object only accepts hashable objects as vertices
        (and as the vertices of the line graph are the edges of the graph), this
        code will fail if edge labels are not hashable. You can also set the
        argument ``labels=False`` to ignore labels.

    .. SEEALSO::

        - The :mod:`line_graph <sage.graphs.line_graph>` module

        - :meth:`~sage.graphs.graph_generators.GraphGenerators.line_graph_forbidden_subgraphs`
          -- the forbidden subgraphs of a line graph

        - :meth:`~Graph.is_line_graph` -- tests whether a graph is a line graph

    EXAMPLES::

        sage: g = graphs.CompleteGraph(4)
        sage: h = g.line_graph()
        sage: h.vertices(sort=True)
        [(0, 1, None),
         (0, 2, None),
         (0, 3, None),
         (1, 2, None),
         (1, 3, None),
         (2, 3, None)]
        sage: h.am()                                                                    # needs sage.modules
        [0 1 1 1 1 0]
        [1 0 1 1 0 1]
        [1 1 0 0 1 1]
        [1 1 0 0 1 1]
        [1 0 1 1 0 1]
        [0 1 1 1 1 0]
        sage: h2 = g.line_graph(labels=False)
        sage: h2.vertices(sort=True)
        [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        sage: h2.am() == h.am()                                                         # needs sage.modules
        True
        sage: g = DiGraph([[1..4], lambda i,j: i < j])
        sage: h = g.line_graph()
        sage: h.vertices(sort=True)
        [(1, 2, None),
         (1, 3, None),
         (1, 4, None),
         (2, 3, None),
         (2, 4, None),
         (3, 4, None)]
        sage: h.edges(sort=True)
        [((1, 2, None), (2, 3, None), None),
         ((1, 2, None), (2, 4, None), None),
         ((1, 3, None), (3, 4, None), None),
         ((2, 3, None), (3, 4, None), None)]

    TESTS:

    :issue:`13787`::

        sage: g = graphs.KneserGraph(7,1)
        sage: C = graphs.CompleteGraph(7)
        sage: C.is_isomorphic(g)
        True
        sage: C.line_graph().is_isomorphic(g.line_graph())
        True
    """
    cdef dict conflicts = {}
    cdef list elist = []

    g._scream_if_not_simple()
    if g._directed:
        from sage.graphs.digraph import DiGraph
        G = DiGraph()
        G.add_vertices(g.edge_iterator(labels=labels))
        for v in g:
            # Connect appropriate incident edges of the vertex v
            G.add_edges((e, f) for e in g.incoming_edge_iterator(v, labels=labels)
                        for f in g.outgoing_edge_iterator(v, labels=labels))
        return G

    from sage.graphs.graph import Graph
    G = Graph()

    # We must sort the edges' endpoints so that (1,2,None) is seen as the
    # same edge as (2,1,None).
    #
    # We do so by comparing hashes, just in case all the natural order (<)
    # on vertices would not be a total order (for instance when vertices are
    # sets). If two adjacent vertices have the same hash, then we store the
    # pair in the dictionary of conflicts

    # 1) List of vertices in the line graph
    for e in g.edge_iterator(labels=labels):
        if hash(e[0]) < hash(e[1]):
            elist.append(e)
        elif hash(e[0]) > hash(e[1]):
            elist.append((e[1], e[0]) + e[2:])
        else:
            # Settle the conflict arbitrarily
            conflicts[e] = e
            conflicts[(e[1], e[0]) + e[2:]] = e
            elist.append(e)

    G.add_vertices(elist)

    # 2) adjacencies in the line graph
    for v in g:
        elist = []

        # Add the edge to the list, according to hashes, as previously
        for e in g.edge_iterator(v, labels=labels):
            if hash(e[0]) < hash(e[1]):
                elist.append(e)
            elif hash(e[0]) > hash(e[1]):
                elist.append((e[1], e[0]) + e[2:])
            else:
                elist.append(conflicts[e])

        # All pairs of elements in elist are edges of the line graph
        while elist:
            x = elist.pop()
            for y in elist:
                G.add_edge(x, y)

    return G


def root_graph(g, verbose=False):
    r"""
    Return the root graph corresponding to the given graph ``g``.

    See the documentation of :mod:`sage.graphs.line_graph` to know how it works.

    INPUT:

    - ``g`` -- a graph

    - ``verbose`` -- boolean (default: ``False``); display some information
      about what is happening inside of the algorithm

    .. WARNING::

        This code assumes that `g` is a line graph, and is a connected,
        undirected graph without multiple edges.

    TESTS:

    All connected graphs on 6 vertices::

        sage: from sage.graphs.line_graph import root_graph
        sage: def test(g):
        ....:    gl = g.line_graph(labels=False)
        ....:    d = root_graph(gl)
        sage: for i,g in enumerate(graphs(6)): # long time
        ....:   if not g.is_connected():
        ....:     continue
        ....:   test(g)

    Non line-graphs::

        sage: root_graph(graphs.PetersenGraph())
        Traceback (most recent call last):
        ...
        ValueError: this graph is not a line graph !

        sage: root_graph(Graph('O{e[{}^~z`MDZBZBkXzE^'))
        Traceback (most recent call last):
        ...
        ValueError: this graph is not a line graph !

    Small corner-cases::

        sage: from sage.graphs.line_graph import root_graph
        sage: root_graph(graphs.CompleteGraph(3))
        (Complete bipartite graph of order 1+3: Graph on 4 vertices,
         {0: (0, 1), 1: (0, 2), 2: (0, 3)})
        sage: G, D = root_graph(graphs.OctahedralGraph()); G
        Complete graph: Graph on 4 vertices
        sage: G, D = root_graph(graphs.DiamondGraph()); G
        Graph on 4 vertices
        sage: G, D = root_graph(graphs.WheelGraph(5)); G
        Diamond Graph: Graph on 4 vertices
    """
    from sage.graphs.digraph import DiGraph

    if isinstance(g, DiGraph):
        raise ValueError("g cannot be a DiGraph !")
    if g.has_multiple_edges():
        raise ValueError("g cannot have multiple edges !")
    if not g.is_connected():
        raise ValueError("g is not connected !")
    # is_line_graph expects a particular error message when g is not a line graph
    not_line_graph = "this graph is not a line graph !"

    # Complete Graph ?
    if g.is_clique():
        from sage.graphs.generators.basic import CompleteBipartiteGraph
        return (CompleteBipartiteGraph(1, g.order()),
                {v: (0, 1 + i) for i, v in enumerate(g)})

    # Diamond Graph ?
    elif g.order() == 4 and g.size() == 5:
        from sage.graphs.graph import Graph
        root = Graph([(0, 1), (1, 2), (2, 0), (0, 3)])
        return (root,
                g.is_isomorphic(root.line_graph(labels=False), certificate=True)[1])

    # Wheel on 5 vertices ?
    elif g.order() == 5 and g.size() == 8 and min(g.degree()) == 3:
        from sage.graphs.generators.basic import DiamondGraph
        root = DiamondGraph()
        return (root,
                g.is_isomorphic(root.line_graph(labels=False), certificate=True)[1])

    # Octahedron ?
    elif g.order() == 6 and g.size() == 12 and g.is_regular(k=4):
        from sage.graphs.generators.platonic_solids import OctahedralGraph
        if g.is_isomorphic(OctahedralGraph()):
            from sage.graphs.generators.basic import CompleteGraph
            root = CompleteGraph(4)
            return (root,
                    g.is_isomorphic(root.line_graph(labels=False), certificate=True)[1])

    # From now on we can assume (thanks to Beineke) that no edge belongs to two
    # even triangles at once.

    # Better to work on integers... Everything takes more time otherwise.
    G = g.relabel(inplace=False)

    # Dictionary of (pairs of) cliques, i.e. the two cliques associated with
    # each vertex.
    cdef dict v_cliques = {v: [] for v in G}
    # All the even triangles we meet
    cdef list even_triangles = []

    # We iterate over all maximal cliques of the graph.

    # As method G.cliques_maximal() returns the list of all maximal cliques
    # (which takes an exponential time on general graphs, while it is obviously
    # polynomial on line graphs), we instead use IndependentSet to iterate over
    # all the maximal cliques. This way, we can say early that the graph is not
    # a line graph.

    from sage.graphs.independent_sets import IndependentSets
    for S in IndependentSets(G, maximal=True, complement=True):

        # Triangles... even or odd ?
        if len(S) == 3:

            # If a vertex of G has an odd number of neighbors among the vertices
            # of S, then the triangle is odd. We compute the list of such
            # vertices by taking the symmetric difference of the neighborhood of
            # our three vertices.
            #
            # Note that the elements of S do not appear in this set as they are
            # all seen exactly twice.

            odd_neighbors = set(G.neighbor_iterator(S[0]))
            odd_neighbors.symmetric_difference_update(G.neighbor_iterator(S[1]))
            odd_neighbors.symmetric_difference_update(G.neighbor_iterator(S[2]))

            # Even triangles
            if not odd_neighbors:
                even_triangles.append(tuple(S))
                continue

            # We manage odd triangles the same way we manage other cliques ...

        # We now associate the clique to all the vertices it contains.
        for v in S:
            if len(v_cliques[v]) == 2:
                raise ValueError(not_line_graph)
            v_cliques[v].append(tuple(S))

        if verbose:
            print("Added clique", S)

    # Deal with even triangles
    for u, v, w in even_triangles:

        # According to Beineke, we must go through all even triangles, and for
        # each triangle uvw consider its three pairs of adjacent vertices uv,
        # vw, wu. For all pairs xy among those such that xy do not appear
        # together in any clique we have found so far, we add xy to the list of
        # cliques describing our covering.

        for x, y in [(u, v), (v, w), (w, u)]:

            # If edge xy does not appear in any of the cliques associated with y
            if all(x not in C for C in v_cliques[y]):
                if len(v_cliques[y]) >= 2 or len(v_cliques[x]) >= 2:
                    raise ValueError("this graph is not a line graph !")

                v_cliques[x].append((x, y))
                v_cliques[y].append((x, y))

                if verbose:
                    print("Adding pair", (x, y),
                          "appearing in the even triangle", (u, v, w))

    # Deal with vertices contained in only one clique. All edges must be defined
    # by TWO endpoints, so we add a fake clique.
    for x, clique_list in v_cliques.items():
        if len(clique_list) == 1:
            clique_list.append((x,))

    # We now have all our cliques. Let's build the root graph to check that it
    # all fits !
    from sage.graphs.graph import Graph
    R = Graph()

    # Associates an integer to each clique
    cdef dict relabel = {}

    # Associates to each vertex of G its pair of coordinates in R
    cdef dict vertex_to_map = {}

    for v, L in v_cliques.items():

        # Add cliques to relabel dictionary
        for S in L:
            if S not in relabel:
                relabel[S] = len(relabel)

        # The coordinates of edge v
        vertex_to_map[v] = relabel[L[0]], relabel[L[1]]

    if verbose:
        print("Final associations :")
        for v, L in v_cliques.items():
            print(v, L)

    # We now build R
    R.add_edges(vertex_to_map.values())

    # If g is a line graph, then it is isomorphic to the line graph of the graph
    # R that we have constructed, so we return R (and the isomorphism).
    is_isom, isom = g.is_isomorphic(R.line_graph(labels=False), certificate=True)
    if is_isom:
        return R, isom
    raise ValueError(not_line_graph)
