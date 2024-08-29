# distutils: language = c++
# distutils: extra_compile_args = -std=c++11
# distutils: libraries = bliss
# sage_setup: distribution = sagemath-bliss

r"""
Interface with bliss: graph (iso/auto)morphism

Implemented functions:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`automorphism_group` | Return the automorphism group of the given (di)graph
    :meth:`canonical_form` | Return a canonical label for the given (di)graph

AUTHORS:

    - Jernej Azarija
"""

# ****************************************************************************
#       Copyright (C) 2015 Jernej Azarija
#       Copyright (C) 2015 Nathann Cohen <nathann.cohen@gail.com>
#       Copyright (C) 2018 Christian Stump <christian.stump@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from libc.limits cimport LONG_MAX

from cysignals.memory cimport check_calloc, sig_free

cdef extern from "bliss/graph.hh" namespace "bliss":

    cdef cppclass Stats:
        pass

    cdef cppclass AbstractGraph:
        pass

    cdef cppclass Graph(AbstractGraph):
        Graph(const unsigned int)
        void add_edge(const unsigned int, const unsigned int)
        void change_color(const unsigned int, const unsigned int)
        const unsigned int* canonical_form(Stats&)

cdef extern from "bliss/digraph.hh" namespace "bliss":

    cdef cppclass Digraph(AbstractGraph):
        Digraph(const unsigned int)
        void add_edge(const unsigned int, const unsigned int)
        void change_color(const unsigned int, const unsigned int)
        const unsigned int* canonical_form(Stats&)
        unsigned int get_hash()

cdef extern from "bliss_cpp/bliss_find_automorphisms.h":

    void bliss_find_automorphisms(Graph*, void (*)(void*, unsigned int, const unsigned int*), void*, Stats&)
    void bliss_find_automorphisms(Digraph*, void (*)(void*, unsigned int, const unsigned int*), void*, Stats&)

cdef int encoding_numbits(int n) noexcept:
    r"""
    Return the number of bits needed to encode the `n` numbers from `1` to
    `n`. In other words, the last bit set in `n`.
    """
    if n <= 0:
        return 0
    cdef int i = 0
    while n:
        n >>= 1
        i += 1
    return i


cdef void add_gen(void *user_param, unsigned int n, const unsigned int *aut) noexcept:
    r"""
    Function called each time a new generator of the automorphism group is
    found.

    This function is used to append the new generators to a Python list. Its
    main job is to translate a permutation into disjoint cycles.

    INPUT:

    - ``user_param`` -- ``void *``; in the current implementation, points toward
      a Python object which is a pair ``(list_of_current_generators,
      vert_to_integer_labelling)``.

    - ``n`` -- integer; number of points in the graph

    - ``aut`` -- ``int *``; an automorphism of the graph
    """
    cdef int N
    cdef int tmp = 0
    cdef int cur = 0
    cdef list perm = []
    cdef bint* done = <bint*> check_calloc(n, sizeof(bint))

    gens, int_to_vertex, N = <object>user_param

    while cur < N:
        if not done[cur]:
            tmp = cur
            cycle = [int_to_vertex[cur]]
            done[cur] = True

            while aut[tmp] != cur:
                tmp = aut[tmp]
                done[tmp] = True
                cycle.append(int_to_vertex[tmp])

            perm.append(tuple(cycle))

        cur += 1

    gens.append(perm)

    sig_free(done)

#####################################################
# constructing bliss graphs from edge lists
#####################################################

cdef Graph *bliss_graph_from_labelled_edges(int Vnr, int Lnr, Vout, Vin, labels, partition) noexcept:
    r"""
    Return a bliss graph from the input data

    For edge labelled graphs, the bliss graph is constructed using `Vnr *
    \log(Lnr)` many vertices as described in Sec. 14 of the `nauty reference
    manual <http://pallini.di.uniroma1.it/Guide.html>`_.

    More precisely, let `V` the vertices of the original graph. We
    construct the new graph on the vertex set
    `V \times \{0, \ldots, log(Lnr)\}`. The integers in the second factor of
    `V \times \{0, \ldots, \log(Lnr)` encode the coloring of the edges. Then

    - for each vertex `v` in `G` and each `i` in `0`, ..., `log(Lnr)-1`, we
      add an edge from `(v, i)` to `(v, i+1)`

    - for each edge `e` from `u` to `v` with label `lab` in `G` we add edges
      from `(u, i)` to `(v, i)` for each `i` so that the `i`-th bit of `lab+1`
      is set (recall that the labels range from `0` to `Lnr-1` and hence the
      binary encoding of `1` to `Lnr` is so that at least one bit is set and
      their binary encoding has length at most `log(Lnr)`)

    .. WARNING::

        the input is not checked for correctness, any wrong input will result in
        a segfault

    INPUT:

    - ``Vnr`` -- integer; number of vertices, such that the vertices are `0,
      \ldots, Vnr-1`

    - ``Lnr`` -- integer; number of labels, such that the labels are `0, \ldots,
      Lnr-1`

    - ``Vout`` -- list; the list of vertices of outgoing edges

    - ``Vin`` -- list; the list of vertices of ingoing edges

    - ``labels`` -- list; the list of edge labels

    - ``partition`` -- an ordered partition of the vertex set
    """
    cdef Graph * g
    cdef int i, j, x, y, lab, Pnr, Enr, logLnr = 1

    if Lnr <= 1:
        g = new Graph(Vnr)
    else:
        logLnr = encoding_numbits(Lnr)
        g = new Graph(Vnr * logLnr)
    if not g:
        raise MemoryError("allocation failed")

    Enr = len(Vout)

    if Lnr <= 1:
        for i in range(Enr):
            g.add_edge(Vout[i], Vin[i])

    else:
        # arrows going up in layers
        for i in range(Vnr * (logLnr - 1)):
            g.add_edge(i, i + Vnr)

        # arrows inside layers shadowing the original graph
        for i in range(Enr):
            x = Vout[i]
            y = Vin[i]
            lab = labels[i] + 1
            j = 0
            while lab:
                if lab & 1:
                    g.add_edge(j * Vnr + x, j * Vnr + y)
                j += 1
                lab >>= 1

    # vertex partition gives colors
    if partition:
        Pnr = len(partition)
        for i in range(Pnr):
            for v in partition[i]:
                for j in range(logLnr):
                    g.change_color(j * Vnr + v, j * Pnr + i)
    else:
        for j in range(logLnr):
            for v in range(Vnr):
                g.change_color(j * Vnr + v, j)

    return g

cdef Digraph *bliss_digraph_from_labelled_edges(int Vnr, int Lnr, Vout, Vin, labels, partition) noexcept:
    r"""
    Return a bliss digraph from the input data

    For edge labelled graphs, the bliss graph is constructed using `Vnr *
    \log(Lnr)` many vertices as described in Sec 14 in the `nauty reference
    manual <http://pallini.di.uniroma1.it/Guide.html>`_.

    .. WARNING::

        the input is not checked for correctness, any wrong input will result in
        a segfault

    INPUT:

    - ``Vnr`` -- integer; number of vertices, such that the vertices are `0,
      \ldots, Vnr-1`

    - ``Lnr`` -- integer; number of labels, such that the labels are `0, \ldots,
      Lnr-1`

    - ``Vout`` -- list; the list of vertices of outgoing edges

    - ``Vin`` -- list; the list of vertices of ingoing edges

    - ``labels`` -- list; the list of edge labels

    - ``partition`` -- a partition of the vertex set
    """
    cdef Digraph *g
    cdef int i, j, x, y, lab, Pnr, Enr, logLnr = 1

    if Lnr <= 1:
        g = new Digraph(Vnr)
    else:
        logLnr = encoding_numbits(Lnr)
        g = new Digraph(Vnr * logLnr)
    if not g:
        raise MemoryError("allocation failed")

    Enr = len(Vout)

    if Lnr <= 1:
        for i in range(Enr):
            g.add_edge(Vout[i], Vin[i])
    else:
        # arrows going up in layers
        for i in range(Vnr * (logLnr - 1)):
            g.add_edge(i, i + Vnr)

        # arrows inside layers shadowing the original graph
        for i in range(Enr):
            x = Vout[i]
            y = Vin[i]
            lab = labels[i] + 1
            j = 0
            while lab:
                if lab & 1:
                    g.add_edge(j * Vnr + x, j * Vnr + y)
                j += 1
                lab >>= 1

    # vertex partition gives color
    if partition:
        Pnr = len(partition)
        for i in range(Pnr):
            for v in partition[i]:
                for j in range(logLnr):
                    g.change_color(j * Vnr + v, j * Pnr + i)
    else:
        for j in range(logLnr):
            for v in range(Vnr):
                g.change_color(j * Vnr + v, j)

    return g

#####################################################
# canonical form from graph or edge list
#####################################################

cdef canonical_form_from_edge_list(int Vnr, list Vout, list Vin, int Lnr=1, list labels=[],
                                   list partition=None, bint directed=False, bint certificate=False) noexcept:
    r"""
    Return an unsorted list of labelled edges of a canonical form.

    INPUT:

    - ``Vnr`` -- integer; number of vertices, such that the vertices are `0,
      \ldots, Vnr-1`

    - ``Vout`` -- list; the list of vertices of outgoing edges

    - ``Vin`` -- list; the list of vertices of ingoing edges

    - ``Lnr`` -- integer (default: 1); number of labels, such that the labels
      are `0, \ldots, Lnr-1`

    - ``labels`` -- list (default: ``[]``); the list of edge labels

    - ``partition`` -- list (default: ``None``); a partition of the vertex
      set

    - ``directed`` -- boolean (default: ``False``); whether the edges are
      directed or not

    - ``certificate`` -- boolean 'default: ``False``); whether to return the
      isomorphism to obtain the canonical labelling
    """
    # We need this to convert the numbers from <unsigned int> to
    # <long>. This assertion should be true simply for memory reasons.
    assert <unsigned long>(Vnr) <= <unsigned long>LONG_MAX

    cdef const unsigned int* aut
    cdef Graph* g
    cdef Digraph* d
    cdef Stats s
    cdef dict relabel

    cdef list new_edges = []
    cdef long e, f

    if directed:
        d = bliss_digraph_from_labelled_edges(Vnr, Lnr, Vout, Vin, labels, partition)
        aut = d.canonical_form(s)
    else:
        g = bliss_graph_from_labelled_edges(Vnr, Lnr, Vout, Vin, labels, partition)
        aut = g.canonical_form(s)

    for i in range(len(Vout)):
        x = Vout[i]
        y = Vin[i]
        e = aut[x]
        f = aut[y]
        if Lnr == 1:
            if not bool(labels):
                lab = 0
            else:
                lab = labels[0]
        else:
            lab = labels[i]
        if directed:
            new_edges.append((e, f, lab))
        else:
            new_edges.append((e, f, lab) if e > f else (f, e, lab))

    if certificate:
        relabel = {v: <long>aut[v] for v in range(Vnr)}

    if directed:
        del d
    else:
        del g

    if certificate:
        return new_edges, relabel
    return new_edges


cpdef canonical_form(G, partition=None, return_graph=False, use_edge_labels=True, certificate=False) noexcept:
    r"""
    Return a canonical label for the given (di)graph.

    A canonical label ``canonical_form(G)`` of ``G`` is a (di)graph defined on
    `\{0,...,n-1\}` such that ``G`` is isomorphic to ``H`` if and only if
    ``canonical_form(G)`` is equal to ``canonical_form(H)``.

    INPUT:

    - ``G`` -- a Sage (Di)Graph

    - ``partition`` -- list (default: ``None``); a partition of the vertices
      of ``G`` into color classes

    - ``return_graph`` -- boolean (default: ``False``); whether to return the
      canonical graph of ``G`` or its set of edges

    - ``use_edge_labels`` -- boolean (default: ``True``); whether to consider
      edge labels. The edge labels are assumed to be hashable and
      sortable. If this is not the case (ie a :exc:`TypeError` is
      raised), the algorithm will consider the string representations
      of the labels instead of the labels.

    - ``certificate`` -- boolean (default: ``False``); when set to ``True``,
      returns the labeling of G into a canonical graph

    TESTS::

        sage: from sage.graphs.bliss import canonical_form                  # optional - bliss
        sage: G = graphs.PetersenGraph()                                    # optional - bliss
        sage: canonical_form(G)                                             # optional - bliss
        [(2, 0, None),
         (2, 1, None),
         (3, 0, None),
         (4, 1, None),
         (5, 3, None),
         (5, 4, None),
         (6, 0, None),
         (6, 4, None),
         (7, 1, None),
         (7, 3, None),
         (8, 2, None),
         (8, 5, None),
         (9, 6, None),
         (9, 7, None),
         (9, 8, None)]

        sage: P = graphs.GeneralizedPetersenGraph(5, 2)                     # optional - bliss
        sage: Q = graphs.PetersenGraph()                                    # optional - bliss
        sage: canonical_form(P) == canonical_form(Q)                        # optional - bliss
        True

        sage: canonical_form(Graph(15), return_graph=True)                  # optional - bliss
        Graph on 15 vertices
        sage: g = digraphs.RandomTournament(40)                             # optional - bliss
        sage: g.is_isomorphic(canonical_form(g, return_graph=True))         # optional - bliss
        True

        sage: # optional - bliss
        sage: g1 = graphs.RandomGNP(100, .4)
        sage: r = Permutations(range(100)).random_element()
        sage: g2 = Graph([(r[u],r[v]) for u,v in g1.edges(sort=True, labels=False)])
        sage: g1 = canonical_form(g1, return_graph=True)
        sage: g2 = canonical_form(g2, return_graph=True)
        sage: g2 == g2
        True

        sage: g = Graph({1: [2]})
        sage: g_ = canonical_form(g, return_graph=True, certificate=True)   # optional - bliss
        sage: 0 in g_[0]                                                    # optional - bliss
        True

    Check that parameter ``use_edge_labels`` can be used (:issue:`27571`)::

        sage: g = Graph({1: {2: 'a'}})
        sage: canonical_form(g, use_edge_labels=True)                       # optional - bliss
        [(1, 0, 'a')]
        sage: canonical_form(g, use_edge_labels=False)                      # optional - bliss
        [(1, 0, None)]

    Check that :issue:`28531` is fixed::

        sage: from itertools import product, permutations
        sage: edges_list = [[(0,1), (1,2)],
        ....:               [(0,1),(1,2),(2,3)],
        ....:               [(0,1),(1,2),(2,3),(3,0)]]
        sage: for edges in edges_list:                                      # optional - bliss
        ....:     for labels in product([0,1], repeat=len(edges)):
        ....:         g = Graph([(u,v,l) for ((u,v),l) in zip(edges, labels)])
        ....:         gcan = canonical_form(g, use_edge_labels=True)
        ....:         for p in permutations(range(g.num_verts())):
        ....:             h = Graph([(p[u], p[v], lab) for u,v,lab in g.edges(sort=True)])
        ....:             hcan = canonical_form(h, use_edge_labels=True)
        ....:             if gcan != hcan: print(edges, labels, p)

    Check that it works with non hashable non sortable edge labels (relying
    on string representations of the labels)::

        sage: # needs sage.modules
        sage: g1 = Graph([(0, 1, matrix(ZZ, 2)), (0, 2, RDF.pi()), (1, 2, 'a')])
        sage: g2 = Graph([(1, 2, matrix(ZZ, 2)), (2, 0, RDF.pi()), (0, 1, 'a')])
        sage: g1can = canonical_form(g1, use_edge_labels=True)              # optional - bliss
        sage: g2can = canonical_form(g2, use_edge_labels=True)              # optional - bliss
        sage: g1can == g2can                                                # optional - bliss
        True

    Check that :issue:`32395` is fixed::

        sage: g = Graph([[0, 2]])  # 1 is not a vertex!
        sage: g.canonical_label(partition=[[0], [1], [2]], algorithm='bliss')  # optional - bliss
        Traceback (most recent call last):
        ...
        ValueError: vertex 1 of the partition is not a vertex of the graph
        sage: g.canonical_label(partition=[[0], [0, 2]], algorithm='bliss')  # optional - bliss
        Traceback (most recent call last):
        ...
        ValueError: vertex 0 can appear only once in the partition
        sage: g.canonical_label(partition=[[0, 0], [2]], algorithm='bliss')  # optional - bliss
        Traceback (most recent call last):
        ...
        ValueError: vertex 0 can appear only once in the partition
        sage: g.canonical_label(partition=[[0]], algorithm='bliss')  # optional - bliss
        Traceback (most recent call last):
        ...
        ValueError: some vertices of the graph are not in the partition
    """
    # We need this to convert the numbers from <unsigned int> to <long>.
    # This assertion should be true simply for memory reasons.
    cdef unsigned long Vnr = G.order()
    assert Vnr <= <unsigned long> LONG_MAX

    cdef bint directed = G.is_directed()

    cdef list Vout = []
    cdef list Vin = []
    cdef list labels = []

    cdef list int2vert
    cdef dict vert2int
    cdef dict lab_to_index
    cdef list edge_labels = [] if use_edge_labels else [None]
    cdef int Lnr = 1

    if partition:
        from itertools import chain
        int2vert = list(chain(*partition))
        # We check that the partition contains only vertices of the graph
        # and that it is actually a partition
        seen = set()
        for u in int2vert:
            if u not in G:
                raise ValueError("vertex {} of the partition is not a vertex of the graph".format(u))
            if u in seen:
                raise ValueError("vertex {} can appear only once in the partition".format(u))
            seen.add(u)
        if len(seen) != G.order():
            raise ValueError("some vertices of the graph are not in the partition")
    else:
        int2vert = list(G)
    vert2int = {v: i for i, v in enumerate(int2vert)}
    if partition:
        partition = [[vert2int[i] for i in part] for part in partition]

    # Create 3 lists to represent edges
    # - Vout[i] : source of the ith edge
    # - Vin[i] : destination of the ith edge
    # - labels[i] : label of the ith edge if use_edge_labels is True
    # On the way, assign a unique integer to each distinct label
    if use_edge_labels:
        try:
            edge_labels = sorted(set(G.edge_labels()))
        except TypeError:
            # NOTE: use edge labels might not be hashable or sortable...
            # rely loosely on string representation
            edge_labels = sorted(set(map(str, G.edge_labels())))
            lab_to_index = {lab: i for i, lab in enumerate(edge_labels)}
            for x, y, lab in G.edge_iterator(labels=True):
                Vout.append(vert2int[x])
                Vin.append(vert2int[y])
                labels.append(lab_to_index[str(lab)])

        else:
            lab_to_index = {lab: i for i, lab in enumerate(edge_labels)}
            for x, y, lab in G.edge_iterator(labels=True):
                Vout.append(vert2int[x])
                Vin.append(vert2int[y])
                labels.append(lab_to_index[lab])

        Lnr = len(lab_to_index)

    else:
        for x, y in G.edge_iterator(labels=False):
            Vout.append(vert2int[x])
            Vin.append(vert2int[y])

    new_edges, relabel = canonical_form_from_edge_list(Vnr, Vout, Vin, Lnr, labels, partition, directed, certificate=True)

    new_edges = [(x, y, edge_labels[lab]) for x, y, lab in new_edges]
    relabel = {int2vert[i]: j for i, j in relabel.items()}

    if return_graph:
        if directed:
            from sage.graphs.digraph import DiGraph
            H = DiGraph(new_edges, loops=G.allows_loops(), multiedges=G.allows_multiple_edges())
        else:
            from sage.graphs.graph import Graph
            H = Graph(new_edges, loops=G.allows_loops(), multiedges=G.allows_multiple_edges())

        H.add_vertices(range(G.order()))
        return (H, relabel) if certificate else H

    # Warning: this may break badly in Python 3 if the graph is not simple
    return (sorted(new_edges), relabel) if certificate else sorted(new_edges)


#####################################################
# automorphism group from graphs
#####################################################

cdef automorphism_group_gens_from_edge_list(int Vnr, Vout, Vin, int Lnr=1, labels=[],
                                            int2vert=[], partition=None, bint directed=False) noexcept:
    r"""
    Return an unsorted list of labelled edges of a canonical form.

    INPUT:

    - ``Vnr`` -- integer; number of vertices, such that the vertices are `0,
      \ldots, Vnr-1`

    - ``Vout`` -- list; the list of vertices of outgoing edges

    - ``Vin`` -- list; the list of vertices of ingoing edges

    - ``Lnr`` -- integer (default: 1); number of labels, such that the labels
      are `0, \ldots, Lnr-1`

    - ``labels`` -- list (default: ``[]``); the list of edge labels

    - ``int2vert`` -- list (default: ``[]``); ordering of the vertices

    - ``partition`` -- list (default: ``None``); a partition of the vertex
      set

    - ``directed`` -- boolean (default: ``False``); whether the edges are
      directed or not
    """
    # We need this to convert the numbers from <unsigned int> to
    # <long>. This assertion should be true simply for memory reasons.
    assert <unsigned long>(Vnr) <= <unsigned long>LONG_MAX

    cdef Graph* g
    cdef Digraph* d
    cdef Stats s

    if not int2vert:
        int2vert = list(range(Vnr))

    cdef list gens = []
    cdef tuple data = (gens, int2vert, Vnr)

    if directed:
        d = bliss_digraph_from_labelled_edges(Vnr, Lnr, Vout, Vin, labels, partition)
        bliss_find_automorphisms(d, add_gen, <void*>data, s)
        del d
    else:
        g = bliss_graph_from_labelled_edges(Vnr, Lnr, Vout, Vin, labels, partition)
        bliss_find_automorphisms(g, add_gen, <void*>data, s)
        del g

    return [[cyc for cyc in gen if cyc[0] is not None] for gen in gens]

cpdef automorphism_group(G, partition=None, use_edge_labels=True) noexcept:
    """
    Return the automorphism group of the given (di)graph.

    Compute the automorphism group of ``G`` subject to the vertex coloring
    ``partition``, if given.  The graph ``G`` can be a directed or undirected
    graph with or without edge labellings.

    Observe the neither the vertex colorings nor the edge colorings are
    interchangeable.

    INPUT:

    - ``G`` -- a Sage graph

    - ``partition`` -- list (default: ``None``); a partition of the vertices
      of ``G`` into color classes. Defaults to ``None``, which is equivalent to
      a partition of size 1.

    - ``use_edge_labels`` -- boolean (default: ``True``); whether to consider
      edge labels

    EXAMPLES::

        sage: from sage.graphs.bliss import automorphism_group                  # optional - bliss

    Computing the automorphism group of a graph or digraph::

        sage: # optional - bliss
        sage: G = graphs.CompleteMultipartiteGraph([1, 1, 1, 2])
        sage: automorphism_group(G).cardinality()
        12
        sage: D = DiGraph(G.edges(sort=True))
        sage: automorphism_group(D).cardinality()
        2

    Observe that the order 12 is given by permuting the first three vertices, or
    the last two in the case of a graph, while only the latter two are possible
    in the case of a directed graph.

    Partitioning the vertices into classes::

        sage: # optional - bliss
        sage: G = graphs.CompleteMultipartiteGraph([3, 2])
        sage: automorphism_group(G).cardinality()
        12
        sage: automorphism_group(G,partition=[[0],[1],[2],[3,4]]).cardinality()
        2
        sage: automorphism_group(G,partition=[[0],[1,2],[3,4]]).cardinality()
        4
        sage: automorphism_group(G,partition=[[1,2],[0,3],[4]]).cardinality()
        2

    Partitioning the edges into classes::

        sage: # optional - bliss
        sage: G = Graph(graphs.CompleteMultipartiteGraph([8, 2]), sparse=True)
        sage: for i,j in G.edges(labels=False, sort=False):
        ....:     if 0 <= i < 3:
        ....:         G.set_edge_label(i, j, "A")
        ....:     if 3 <= i < 6:
        ....:         G.set_edge_label(i, j, "B")
        ....:     if 6 <= i < 8:
        ....:         G.set_edge_label(i, j, "C")
        sage: factor(automorphism_group(G).cardinality())
        2^4 * 3^2
        sage: automorphism_group(G,[[0],[1],[2,3],[4,5],[6,7],[8],[9]]).cardinality()
        4

    TESTS::

        sage: from sage.graphs.bliss import automorphism_group                  # optional - bliss
        sage: G = graphs.PetersenGraph()                                        # optional - bliss
        sage: automorphism_group(G).is_isomorphic(G.automorphism_group())       # optional - bliss
        True

        sage: # optional - bliss
        sage: G = graphs.HeawoodGraph()
        sage: p = G.bipartite_sets()
        sage: A = G.automorphism_group(partition=[list(p[0]), list(p[1])])
        sage: automorphism_group(G, partition=p).is_isomorphic(A)
        True

        sage: G = graphs.CompleteMultipartiteGraph([5, 7, 11])
        sage: B = automorphism_group(G)                                         # optional - bliss
        sage: B.cardinality() == prod(factorial(n) for n in [5, 7, 11])         # optional - bliss
        True

        sage: # optional - bliss
        sage: G = Graph(graphs.CompleteMultipartiteGraph([8,8,8,5]),sparse=True)
        sage: for i,j in G.edges(labels=False, sort=False):
        ....:     if 0 <= i < 3:
        ....:         G.set_edge_label(i, j, "A")
        ....:     if 3 <= i < 6:
        ....:         G.set_edge_label(i, j, "B")
        ....:     if 6 <= i < 8:
        ....:         G.set_edge_label(i, j, "C")
        sage: card = automorphism_group(G).cardinality()
        sage: card == prod(factorial(n) for n in [3, 3, 2, 8, 8, 5, 2])
        True
        sage: card = automorphism_group(G, use_edge_labels=False).cardinality()
        sage: card == prod(factorial(n) for n in [8, 8, 8, 5, 3])
        True
        sage: card = automorphism_group(G, [[0 .. 7], [8 .. 11] ,[12 .. 28]]).cardinality()
        sage: card == prod(factorial(n) for n in [3, 3, 2, 4, 4, 8, 5])
        True

        sage: # optional - bliss
        sage: G = Graph()
        sage: G.add_edges((i,j,"A") for i in range(0, 2) for j in range(14,20))
        sage: G.add_edges((i,j,"B") for i in range(2, 5) for j in range(14,20))
        sage: G.add_edges((i,j,"C") for i in range(5, 9) for j in range(14,20))
        sage: G.add_edges((i,j,"D") for i in range(9,14) for j in range(14,20))
        sage: A = automorphism_group(G)
        sage: print(A.gens())               # random
        ((12,13), (11,12), (10,11), (9,10), (7,8), (6,7), (5,6), (3,4),
         (2,3), (1,0), (18,19), (17,18), (16,17), (15,16), (14,15))
        sage: A.cardinality() == prod(factorial(n) for n in [2,3,4,5,6])
        True

        sage: # optional - bliss
        sage: G = Graph()
        sage: alpha = "abcdefghijklmnopqrstuvwxyz"
        sage: G.add_edges((alpha[i],alpha[j],"A") for i in range(0, 2) for j in range(14,20))
        sage: G.add_edges((alpha[i],alpha[j],"B") for i in range(2, 5) for j in range(14,20))
        sage: G.add_edges((alpha[i],alpha[j],"C") for i in range(5, 9) for j in range(14,20))
        sage: G.add_edges((alpha[i],alpha[j],"D") for i in range(9,14) for j in range(14,20))
        sage: A = automorphism_group(G)
        sage: print(A.gens())
        (('m','n'), ('l','m'), ('k','l'), ('j','k'), ('h','i'),
         ('g','h'), ('f','g'), ('d','e'), ('c','d'), ('s','t'),
         ('r','s'), ('q','r'), ('p','q'), ('o','p'), ('a','b'))
        sage: A.cardinality() == prod(factorial(n) for n in [2, 3, 4, 5, 6])
        True

        sage: # optional - bliss
        sage: gg = graphs.CompleteGraph(5)
        sage: gg.allow_loops(True)
        sage: gg.add_edge(0, 0)
        sage: gg.add_edge(1, 1)
        sage: automorphism_group(gg).cardinality()
        12
        sage: automorphism_group(gg, [[0], [1, 2, 3, 4]]).cardinality()
        6
    """
    # We need this to convert the numbers from <unsigned int> to
    # <long>. This assertion should be true simply for memory reasons.
    cdef unsigned long Vnr = G.order()
    assert Vnr <= <unsigned long>LONG_MAX

    cdef bint directed = G.is_directed()

    cdef int labInd
    cdef list Vout = []
    cdef list Vin = []
    cdef list labels = []

    cdef list int2vert
    cdef dict vert2int
    cdef list edge_labels = []
    cdef int Lnr = 0 if use_edge_labels else 1

    if partition:
        from itertools import chain
        int2vert = list(chain(*partition))
    else:
        int2vert = list(G)
    vert2int = {v: i for i, v in enumerate(int2vert)}
    if partition:
        partition = [[vert2int[i] for i in part] for part in partition]

    # Create 3 lists to represent edges
    # - Vout[i] : source of the ith edge
    # - Vin[i] : destination of the ith edge
    # - labels[i] : label of the ith edge if use_edge_labels is True
    # On the way, assign a unique integer to each distinct label
    for x, y, lab in G.edge_iterator(labels=True):
        Vout.append(vert2int[x])
        Vin.append(vert2int[y])
        if use_edge_labels:
            try:
                labInd = edge_labels.index(lab)
            except ValueError:
                labInd = Lnr
                Lnr += 1
                edge_labels.append(lab)
            labels.append(labInd)

    gens = automorphism_group_gens_from_edge_list(Vnr, Vout, Vin, Lnr, labels, int2vert, partition, directed)

    from sage.groups.perm_gps.permgroup import PermutationGroup
    return PermutationGroup(gens, domain=int2vert[:G.order()])


#####################################################
# old direct interactions graphs <-> bliss graphs
#####################################################

cdef Graph *bliss_graph(G, partition, vert2int, int2vert) noexcept:
    r"""
    Return a bliss copy of a graph G

    INPUT:

    - ``G`` -- a Sage Graph

    - ``partition`` -- list; a partition of the vertex set

    - ``vert2int, int2vert`` -- a empty ``dict`` and a empty list; the
      entries of the dictionary are later set to record the labeling of our
      graph. They are taken as arguments to avoid technicalities of returning
      Python objects in Cython functions.
    """
    cdef Graph *g = new Graph(G.order())

    if not g:
        raise MemoryError("allocation failed")

    for i, v in enumerate(G):
        vert2int[v] = i
        int2vert[i] = v

    for x, y in G.edge_iterator(labels=False):
        g.add_edge(vert2int[x], vert2int[y])

    if partition:
        for i in range(1, len(partition)):
            for v in partition[i]:
                g.change_color(vert2int[v], i)
    return g


cdef Digraph *bliss_digraph(G, partition, vert2int, int2vert) noexcept:
    r"""
    Return a bliss copy of a digraph G

    INPUT:

    - ``G`` -- a Sage DiGraph

    - ``partition`` -- list; a partition of the vertex set

    - ``vert2int, int2vert`` -- a empty ``dict`` and a empty list; the
      entries of the dictionary are later set to record the labeling of our
      graph. They are taken as arguments to avoid technicalities of returning
      Python objects in Cython functions.
    """
    cdef Digraph *g = new Digraph(G.order())

    if not g:
        raise MemoryError("allocation failed")

    for i, v in enumerate(G):
        vert2int[v] = i
        int2vert[i] = v

    for x, y in G.edge_iterator(labels=False):
        g.add_edge(vert2int[x], vert2int[y])

    if partition:
        for i in range(1, len(partition)):
            for v in partition[i]:
                g.change_color(vert2int[v], i)
    return g
