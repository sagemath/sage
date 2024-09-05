# cython: binding=True
# distutils: language = c++
# distutils: extra_compile_args = -std=c++11
r"""
Slice decomposition

This module implements an extended lexBFS algorithm for computing the slice
decomposition of undirected graphs and the class :class:`~SliceDecomposition` to
represent such decompositions.

A formal definition of slice decompositions can be found in Section 3.2 of
[TCHP2008]_ and a description of the extended lexBFS algorithm is given in
Section 3.3 of [TCHP2008]_.

AUTHORS:

- Cyril Bouvier (2024-06-25): initial version
"""
# ****************************************************************************
# Copyright (C) 2024 Cyril Bouvier <cyril.bouvier@lirmm.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from libcpp.algorithm cimport swap
from cython.operator cimport dereference as deref

from sage.graphs.base.c_graph cimport CGraphBackend
from sage.data_structures.bitset_base cimport bitset_in

cdef void extended_lex_BFS(
        CGraph cg, vector[int] &sigma, vector[int] *sigma_inv,
        int initial_v_int, vector[int] *pred, vector[size_t] *xslice_len,
        vector[vector[int]] *lex_label) except *:
    r"""
    Perform a extended lexicographic breadth first search (LexBFS) on the
    undirected graph `G`.

    In addition to computing a LexBFS ordering, the extended LexBFS algorithm
    can be used to compute the slice decomposition of the graph.

    This function implements the `O(n+m)` time algorithm proposed in [HMPV2000]_
    and [TCHP2008]_.

    INPUT:

    - ``cg`` -- a ``CGraph``. This function ignores loops and multiple edges and
      assumes that the graph is undirected.

    - ``sigma`` -- vector of ``int`` to store the ordering of the vertices
      resulting from the LexBFS traversal. At the end, the vector will have size
      `n` (the number of vertices of the graph).

    - ``sigma_inv`` -- a pointer to a vector to store the inverse of the
      permutation ``sigma``. ``sigma_inv`` can be ``NULL`` if the caller does
      not need it (but, note that, the inverse of ``sigma`` is still needed by
      the algorithm, so it does not save time nor memory to have ``sigma_inv``
      equal to ``NULL``). At the end, if ``sigma_inv`` is not NULL, the vector
      pointer by it will have size `n` (the number of vertices of the graph)
      and will satisfy:
        * sigma[deref(sigma_inv)[v_int]] = v_int
        * deref(sigma_inv)[sigma[i]] = i

    - ``initial_v_int`` -- the first vertex to consider. It can be `-1`; in this
      case the first active vertex (corresponding to the first bit set in
      ``cg.active_vertices``) will be taken as first vertex.

    - ``pred`` -- a pointer to a vector of int to store the predecessor of a
      vertex in the LexBFS traversal. ``pred`` can be ``NULL`` if the caller
      does not need it (and the information will not be computed by the
      algorithm). At the end, if ``pred`` is not NULL, the vector pointer by it
      will have size `n` (the number of vertices of the graph) and pred[i] will
      be either -1 (if sigma[i] as no predecessor) or a positive value less than
      n such that the predecessor of sigma[i] is sigma[pred[i]].

    - ``xslice_len`` -- a pointer to a vector of size_t to store the length of
      the x-slices associated with the lexBFS traversal. ``xslice_len`` can be
      ``NULL`` if the caller does not need it (and the information will not be
      computed by the algorithm). At the end, if ``xslice_len`` is not NULL, the
      vector pointer by it will have size `n` (the number of vertices of the
      graph) and the length of the x-slice starting at sigma[i] will be
      xslice_len[i].

    - ``lex_label`` -- a pointer to a vector of vector[int] to store the
      lexicographic labels associated with the lexBFS traversal. ``lex_label``
      can be ``NULL`` if the caller does not need it (and the information will
      not be computed by the algorithm). At the end, if ``lex_label`` is not
      NULL, the vector pointer by it will have size `n` (the number of
      vertices of the graph) and the lexicographic label of sigma[i]
      will given by lex_label[i].

    ALGORITHM:

    This algorithm uses the notion of *partition refinement* to determine the
    exact position of the vertices in the ordering.

    Consider an ordering `\sigma` of the vertices. For a vertex `v`, we define
    `N_i(v) = \{u | u \in N(v) \text{ and } \sigma(u) < i\}`, that is the subset
    of neighbors of `v` appearing before the `i`-th vertex in the ordering
    `\sigma`. Now, a part of an ordering `\sigma` is a set of consecutive
    vertices, `S = \{u | i \leq \sigma(u) \leq j\}`, such that for any `u \in
    S`, we have `N_i(u) = N_i(\sigma^{-1}(i))` and for any `v` such that `j <
    \sigma(v)`, `N_i(v) \neq N_i(\sigma^{-1}(i))`. The *head* of a part is the
    first position of its vertices.

    The algorithm starts with a single part containing all vertices. Then, when
    the position of the `i`-th vertex `v` is fixed, it explores the neighbors of
    `v` that have not yet been ordered. Consider a part `S` such that `N(x)\cap
    S \neq \emptyset`. The algorithm will rearrange the ordering of the vertices
    in `S` so that the first vertices are the neighbors of `v`. The subpart
    containing the neighbors of `v` is assigned a new name, and the head of `S`
    is set to the position of the first vertex of `S \setminus N(v)` in the
    ordering `\sigma`.

    Observe that each arc of the graph can induce the subdivision of a part.
    Hence, the algorithm can use up to `m + 1` different parts.

    The time complexity of this algorithm is in `O(n + m)`, and our
    implementation follows that complexity ``SparseGraph``. For ``DenseGraph``,
    the complexity is `O(n^2)`. See [HMPV2000]_ and [TCHP2008]_ for more
    details.

    This implementation of extended LexBFS offers some guarantee on the order in
    which the vertices appear in the computed ordering: in case of a tie between
    lexicographic labels during the computation, this function will "choose" the
    vertices in the order in which they appear during the enumeration of the
    neighbors of their last common neighbor.
    For example, if `(u_0, ..., u_k)` is the beginning of the ordering being
    computed and the vertices `v` and `w` currently have the same lexicographic
    label (it means that they have the same neighbors in `(u_0, ..., u_k)`).
    Let call `u_j` their last neighbor in the current ordering (*i.e.*, for all
    `i > j`, `u_i` is not a neighbor of `v` and `w`). This implementation
    will choose `v` for the next vertex of the ordering if and only if `v`
    appeared before `w` when the neighbors of `u_j` where enumerated.

    One possible use of this guarantee is that the caller can reorder the
    adjacency list of vertices (by using, for example, a static sparse graph)
    to force the computed LexBFS order to respect a previous one.

    EXAMPLES::

    To see how it can be used, see the code of the lex_BFS method (in
    traversals.pyx) or of the class SliceDecomposition in this module.

    TESTS:

    Indirect doctests::

        sage: G = graphs.HouseGraph()
        sage: G.slice_decomposition()
        [0[1[2]] [3] [4]]
        sage: G.lex_BFS(algorithm="fast")
        [0, 1, 2, 3, 4]
    """
    cdef int n = <int> cg.num_verts
    # Variables for the partition refinement algorithm
    cdef size_t max_nparts = cg.num_arcs // 2 + 1
    cdef bint need_to_delete_sigma_inv = sigma_inv == NULL
    if sigma_inv == NULL:
        sigma_inv = new vector[int]()
    cdef vector[size_t] part_of = vector[size_t](n, 0)
    cdef vector[size_t] part_len  # only used if xslice_len != NULL (see below)
    cdef vector[size_t] part_head = vector[size_t](max_nparts)
    cdef vector[size_t] subpart = vector[size_t](max_nparts)
    cdef size_t p, part_of_i, nparts, old_nparts
    # Temporary variables
    cdef int max_degree = 0
    cdef int i, j, k, l, u_int, v_int, t_int

    # Resize vectors
    sigma.resize(n)
    deref(sigma_inv).resize(cg.active_vertices.size)
    if pred != NULL:
        deref(pred).clear()
        deref(pred).resize(n, -1)  # initialize pred[i] to -1 for 0 <= i < n
    if xslice_len != NULL:
        deref(xslice_len).resize(n)
        part_len.resize(max_nparts)
    if lex_label != NULL:
        deref(lex_label).resize(n)

    # Initialize the position of vertices in sigma (and compute max_degree)
    if initial_v_int >= 0:
        sigma[0] = initial_v_int
        deref(sigma_inv)[initial_v_int] = 0
        i = 1
    else:
        i = 0
    for v_int in range(<int> cg.active_vertices.size):
        if bitset_in(cg.active_vertices, v_int):
            if v_int != initial_v_int:
                sigma[i] = v_int
                deref(sigma_inv)[v_int] = i
                i = i + 1
            max_degree = max(max_degree, cg.out_degrees[v_int])

    # Variables needed to iterate over neighbors of a vertex
    cdef int nneighbors
    cdef vector[int] neighbors = vector[int](max_degree)

    # Initialize partition: one part containing all the vertices
    nparts = 1
    # all element of part_of are already initialized to 0
    part_head[0] = 0
    subpart[0] = 0
    if xslice_len != NULL:
        part_len[0] = n

    # Main loop
    for i in range(n):
        old_nparts = nparts

        part_of_i = part_of[i]

        # put i out of its part (updating part_len if needed)
        part_head[part_of_i] += 1
        if xslice_len != NULL:
            deref(xslice_len)[i] = part_len[part_of_i]
            part_len[part_of_i] -= 1

        v_int = sigma[i]

        # Iterate over the neighbors of v
        nneighbors = cg.out_neighbors_unsafe (v_int, neighbors.data(), max_degree)
        for k in range(nneighbors):
            u_int = neighbors[k]
            j = deref(sigma_inv)[u_int]  # get the position of u
            if j <= i:
                continue  # already taken care of

            if lex_label != NULL:
                deref(lex_label)[j].push_back (v_int)

            p = part_of[j]  # get the part of u
            l = part_head[p]  # get the beginning of the part containing u

            # if not last and next elem belongs in the same part (ie #part >= 2)
            if l < n - 1 and part_of[l + 1] == p:
                if l != j:  # not already first elem of the part
                    # Place u at the position of the head of the part
                    t_int = sigma[l]
                    deref(sigma_inv)[t_int], deref(sigma_inv)[u_int] = j, l
                    sigma[j], sigma[l] = t_int, u_int
                    if lex_label != NULL:
                        swap[vector[int]](deref(lex_label)[j],
                                          deref(lex_label)[l])
                    j = l
                part_head[p] += 1  # move the head of the part to next elem

            # if part p was not already cut in two during this iteration, we
            # create a new part using subpart
            if subpart[p] < old_nparts:
                subpart[p] = nparts
                part_head[nparts] = j
                if xslice_len != NULL:
                    part_len[nparts] = 0
                subpart[nparts] = 0
                nparts += 1

            # Finally, we update the name of the part for position j and set v
            # as predecessor of u
            part_of[j] = subpart[p]
            if xslice_len != NULL:
                part_len[p] -= 1
                part_len[subpart[p]] += 1
            if pred != NULL:
                deref(pred)[j] = i

    if need_to_delete_sigma_inv:
        del sigma_inv


def slice_decomposition(G, initial_vertex=None):
    r"""
    Compute a slice decomposition of the simple undirected graph

    INPUT:

    - ``G`` -- a Sage graph.

    - ``initial_vertex`` -- (default: ``None``); the first vertex to consider.

    OUTPUT:

    An object of type :class:`~sage.graphs.graph_decompositions.slice_decomposition.SliceDecomposition`
    that represents a slice decomposition of ``G``

    .. NOTE::

        Loops and multiple edges are ignored during the computation of the slice
        decomposition.

    ALGORITHM:

    The method use the algorithm based on "partition refinement" described in
    [HMPV2000]_ and [TCHP2008]_.
    The time complexity of this algorithm is in `O(n + m)`, and our
    implementation follows that complexity for ``SparseGraph``. For
    ``DenseGraph``, the complexity is `O(n^2)`.

    EXAMPLES:

    Slice decomposition of the Petersen Graph::

        sage: G = graphs.PetersenGraph()
        sage: SD = G.slice_decomposition(); SD
        [0[1[4[5]]] [2[6]] [3] [9] [7] [8]]

    The graph can have loops or multiple edges but they are ignored::

        sage: H = Graph(G,loops=True,multiedges=True)
        sage: H.add_edges([(4, 4), (2, 2), (1, 6)])
        sage: SD2 = H.slice_decomposition()
        sage: SD2 == SD
        True
        sage: SD2.underlying_graph() == G.to_simple(immutable=True)
        True

    The tree corresponding to the slice decomposition can be displayed using
    ``view``::

        sage: from sage.graphs.graph_latex import check_tkz_graph
        sage: check_tkz_graph()  # random - depends on Tex installation
        sage: view(G)  # not tested
        sage: latex(G)  # to obtain the corresponding LaTeX code
        \begin{tikzpicture}
        ...
        \end{tikzpicture}

    Slice decompositions are only defined for undirected graphs::

        sage: from sage.graphs.graph_decompositions.slice_decomposition import slice_decomposition
        sage: slice_decomposition(DiGraph())
        Traceback (most recent call last):
        ...
        ValueError: parameter G must be an undirected graph
    """
    return SliceDecomposition(G, initial_vertex=initial_vertex)


cdef class SliceDecomposition(SageObject):

    def __init__(self, G, initial_vertex=None):
        r"""
        Represents a slice decomposition of a simple directed graph.

        INPUT:

        - ``G`` -- a Sage graph.

        - ``initial_vertex`` -- (default: ``None``); the first vertex to
          consider.

        .. SEEALSO::

            * :meth:`~slice_decomposition` -- compute a slice decomposition of
              the simple undirected graph
            * Section 3.2 of [TCHP2008]_ for a formal definition.

        EXAMPLES:

        The constructor of the :class:`~SliceDecomposition` class is called by
        the :meth:`~slice_decomposition` method of undirected graphs::

            sage: from sage.graphs.graph_decompositions.slice_decomposition import SliceDecomposition
            sage: G = graphs.PetersenGraph()
            sage: SliceDecomposition(G) == G.slice_decomposition()
            True

        The vertex appearing first in the slice decomposition can be specified::

            sage: from sage.graphs.graph_decompositions.slice_decomposition import SliceDecomposition
            sage: SliceDecomposition(graphs.PetersenGraph(), initial_vertex=3)
            [3[2[4[8]]] [1[7]] [0] [9] [6] [5]]

        Slice decompositions are not defined for directed graphs::

            sage: from sage.graphs.graph_decompositions.slice_decomposition import SliceDecomposition
            sage: SliceDecomposition(DiGraph())
            Traceback (most recent call last):
            ...
            ValueError: parameter G must be an undirected graph

        .. automethod:: __getitem__
        """
        if G.is_directed():
            raise ValueError("parameter G must be an undirected graph")

        if initial_vertex is not None and initial_vertex not in G:
            raise LookupError(f"vertex ({initial_vertex}) is not a vertex of the graph")

        cdef CGraphBackend Gbackend = <CGraphBackend> G._backend
        cdef CGraph cg = Gbackend.cg()

        self._graph_class = type(G)

        cdef int initial_v_int
        if initial_vertex is not None:
            # we already checked that initial_vertex is in G
            initial_v_int = Gbackend.get_vertex(initial_vertex)
        else:
            initial_v_int = -1

        cdef vector[int] sigma
        cdef vector[vector[int]] lex_label

        # Compute the slice decomposition using the extended lexBFS algorithm
        extended_lex_BFS(cg, sigma, NULL, initial_v_int, NULL,
                         &(self.xslice_len), &lex_label)

        # Translate the results with the actual vertices of the graph
        self.sigma = tuple(Gbackend.vertex_label(v_int) for v_int in sigma)
        self.sigma_inv = {v: i  for i, v in enumerate(self.sigma)}
        self.lex_label = {i: tuple(Gbackend.vertex_label(v_int) for v_int in lli)
                                    for i, lli in enumerate(lex_label)}

    def __eq__(self, other):
        """
        Return whether ``self`` and ``other`` are equal.

        TESTS::

            sage: G = graphs.PetersenGraph()
            sage: SD = G.slice_decomposition()
            sage: SD == SD
            True
            sage: SD == G.slice_decomposition()
            True

            sage: P3 = graphs.PathGraph(3)
            sage: SD1 = P3.slice_decomposition(initial_vertex=0)
            sage: SD2 = P3.slice_decomposition(initial_vertex=2)
            sage: SD1 == SD2
            False
            sage: SD3 = graphs.CompleteGraph(3).slice_decomposition()
            sage: SD1 == SD3  # same lexBFS but different slice for 1
            False
            sage: SD4 = Graph([(0,1), (0,2)]).slice_decomposition()
            sage: SD3 == SD4  # same lexBFS and slices but different active edges
            False
        """
        if not isinstance(other, type(self)):
            return False

        cdef SliceDecomposition sd = <SliceDecomposition>other

        return self.sigma_inv == sd.sigma_inv \
                and self.lex_label == sd.lex_label \
                and self.xslice_len == sd.xslice_len

    def __hash__(self):
        r"""
        Compute a hash of a ``SliceDecomposition`` object.

        TESTS::

            sage: P3 = graphs.PathGraph(3)
            sage: hash(P3.slice_decomposition(initial_vertex=0))
            -7313201005658437102
            sage: hash(P3.slice_decomposition(initial_vertex=2))
            1181676064626878036
            sage: hash(graphs.CompleteGraph(3).slice_decomposition())
            6162668211142297415
            sage: hash(Graph([(0,1), (0,2)]).slice_decomposition())
            2898184589667302557
        """
        return hash((tuple(self.sigma_inv.items()),
                     tuple(self.lex_label.items()),
                     tuple(self.xslice_len)))

    def __getitem__(self, v):
        r"""
        Return the data about the x-slice of the vertex `v`.

        INPUT:

        - ``v`` -- a vertex of the graph corresponding to the slice
          decomposition.

        OUTPUT:

        A dictionnary with the keys:

        * ``"pivot"`` -- the vertex `v` given as parameter

        * ``"slice"`` -- the slice of `v` (see :meth:`~slice`)

        * ``"active_edges"`` -- the actives edges of `v` (see
          :meth:`~active_edges`)

        * ``"lexicographic_label"`` -- the lexicographic label of `v` (see
          :meth:`~lexicographic_label`)

        * ``"sequence"`` -- the x-slice sequence of `v` (see
          :meth:`~xslice_sequence`)

        This method can also be called via :meth:`xslice_data`.

        EXAMPLES:

        ::

            sage: G = Graph('L~mpn~Nrv{^o~_').relabel('abcdefguvwxyz',inplace=False)
            sage: SD = G.slice_decomposition(initial_vertex='x')
            sage: SD.xslice_data('a')
            {'active_edges': [('a', 'b'),
              ('a', 'c'),
              ('a', 'd'),
              ('a', 'e'),
              ('a', 'f'),
              ('c', 'g'),
              ('d', 'g'),
              ('f', 'g')],
             'lexicographic_label': ['x'],
             'pivot': 'a',
             'sequence': [['a'], ['b', 'c', 'd', 'e', 'f'], ['g']],
             'slice': ['a', 'b', 'c', 'd', 'e', 'f', 'g']}
            sage: SD.xslice_data('u')
            {'active_edges': [],
             'lexicographic_label': ['a', 'b', 'c', 'd', 'e', 'f', 'g'],
             'pivot': 'u',
             'sequence': [['u'], ['y', 'z']],
             'slice': ['u', 'y', 'z']}

        Some values of the returned dictionnary can be obtained via other
        methods (:meth:`~slice`, :meth:`~xslice_sequence`,
        :meth:`~active_edges`, :meth:`~lexicographic_label`)::

            sage: SD.slice('a')
            ['a', 'b', 'c', 'd', 'e', 'f', 'g']
            sage: SD.xslice_data('a')['slice']
            ['a', 'b', 'c', 'd', 'e', 'f', 'g']

            sage: SD.xslice_sequence('a')
            [['a'], ['b', 'c', 'd', 'e', 'f'], ['g']]
            sage: SD.xslice_data('a')['sequence']
            [['a'], ['b', 'c', 'd', 'e', 'f'], ['g']]

            sage: SD.active_edges('b') == SD.xslice_data('b')['active_edges']
            True

            sage: SD.lexicographic_label('u')
            ['a', 'b', 'c', 'd', 'e', 'f', 'g']
            sage: SD.xslice_data('u')['lexicographic_label']
            ['a', 'b', 'c', 'd', 'e', 'f', 'g']

        TESTS::

            sage: G = graphs.RandomGNP(15, 0.3)
            sage: SD = G.slice_decomposition()
            sage: all(SD[v]['slice'] == SD.slice(v) for v in G)
            True
            sage: all(SD[v]['sequence'] == SD.xslice_sequence(v) for v in G)
            True
            sage: all(SD[v]['active_edges'] == SD.active_edges(v) for v in G)
            True
            sage: all(SD[v]['lexicographic_label'] == SD.lexicographic_label(v) for v in G)
            True

            sage: SD = graphs.PetersenGraph().slice_decomposition()
            sage: SD['John']
            Traceback (most recent call last):
            ...
            LookupError: vertex (John) does not appear in the slice decomposition
        """
        if v not in self.sigma_inv:
            raise LookupError(f"vertex ({v}) does not appear in the slice "
                              "decomposition")
        cdef size_t i = self.sigma_inv[v]
        return {'pivot': v,
                'slice': self._slice(i),
                'sequence': self._xslice_sequence(i),
                'lexicographic_label': self._xslice_lex_label(i),
                'active_edges': self._xslice_active_edges(i),
                }

    def lexBFS_order(self):
        r"""
        Return the lexBFS order corresponding to the slice decomposition.

        EXAMPLES::

            sage: from sage.graphs.traversals import _is_valid_lex_BFS_order
            sage: G = graphs.PetersenGraph(); SD = G.slice_decomposition()
            sage: SD.lexBFS_order()
            [0, 1, 4, 5, 2, 6, 3, 9, 7, 8]
            sage: _is_valid_lex_BFS_order(G, SD.lexBFS_order())
            True

        TESTS::

            sage: from sage.graphs.traversals import _is_valid_lex_BFS_order
            sage: for _ in range(5):
            ....:   G = graphs.RandomGNP(15, 0.3)
            ....:   SD = G.slice_decomposition()
            ....:   _is_valid_lex_BFS_order(G, SD.lexBFS_order())
            True
            True
            True
            True
            True
        """
        return list(self.sigma)

    def xslice_data(self, v):
        r"""
        Return the data about the x-slice of the vertex `v`.

        This method is a wrapper around :meth:`SliceDecomposition.__getitem__`

        TESTS::

            sage: G = graphs.RandomGNP(15, 0.3)
            sage: SD = G.slice_decomposition()
            sage: all(SD[v] == SD.xslice_data(v) for v in G)
            True
        """
        return self[v]

    def slice(self, v):
        r"""
        Return the slice of the vertex `v`.

        The slice of `v` is the list of vertices `u` such that the neighbors of
        `u` that are before `v` in the lexBFS order are that same that the
        neighbors of `v` that are before `v` in the lexBFS order (*i.e.*, the
        lexicographic label of `v`). It can be shown that it is a factor of the
        lexBFS order.

        INPUT:

        - ``v`` -- a vertex of the graph corresponding to the slice
          decomposition.

        OUTPUT:

        A list of vertices

        EXAMPLES:

        ::

            sage: G = Graph('L~mpn~Nrv{^o~_').relabel('abcdefguvwxyz',inplace=False)
            sage: SD = G.slice_decomposition(initial_vertex='x')
            sage: SD.slice('a')
            ['a', 'b', 'c', 'd', 'e', 'f', 'g']

        The vertices of the slice have the same neighborhood "on the left"::

            sage: pos = lambda v: SD.lexBFS_order().index(v)
            sage: lla = set(SD.lexicographic_label('a'))
            sage: all(lla == {u for u in G.neighbors(v) if pos(u) < pos('a')} \
            ....:       for v in SD.slice('a'))
            True

        The slice is a factor of the lexBFS order::

            sage: ''.join(SD.slice('a')) in ''.join(SD.lexBFS_order())
            True

        The slice of the initial vertex is the whole graph::

            sage: SD.slice('x') == SD.lexBFS_order()
            True

        TESTS::

            sage: SD.slice('Michael')
            Traceback (most recent call last):
            ...
            LookupError: vertex (Michael) does not appear in the slice decomposition
        """
        if v not in self.sigma_inv:
            raise LookupError(f"vertex ({v}) does not appear in the slice "
                              "decomposition")
        cdef size_t i = self.sigma_inv[v]
        return self._slice(i)

    def xslice_sequence(self, v):
        r"""
        Return the x-slice sequence of the vertex `v`.

        INPUT:

        - ``v`` -- a vertex of the graph corresponding to the slice
          decomposition.

        OUTPUT:

        A list of list corresponding to the x-slice sequence of ``v``.

        EXAMPLES:

        ::

            sage: G = Graph('L~mpn~Nrv{^o~_').relabel('abcdefguvwxyz',inplace=False)
            sage: SD = G.slice_decomposition(initial_vertex='x')
            sage: SD.xslice_sequence('x')
            [['x'], ['a', 'b', 'c', 'd', 'e', 'f', 'g'], ['u', 'y', 'z'], ['v', 'w']]
            sage: SD.xslice_sequence('a')
            [['a'], ['b', 'c', 'd', 'e', 'f'], ['g']]

        The flatten x-slice sequence of a vertex corresponds to the slice of the
        same vertex::

            sage: from itertools import chain
            sage: all(list(chain(*SD.xslice_sequence(v))) == SD.slice(v) \
            ....:       for v in G)
            True

        The first list of the sequence is always a singleton containing the
        input vertex::

            sage: all(SD.xslice_sequence(v)[0] == [v] for v in G)
            True

        If the length of the slice if more than 1, the second list of the
        sequence is either, all the remaining vertices of the slice of `v`, if
        `v` is isolated in the subgraph induced by the slice of `v`, or the
        neighbors of `v` in the subgraph induced by the slice of `v`::

            sage: all(SD.xslice_sequence(v)[1] == SD.slice(v)[1:] for v in G \
            ....:           if G.subgraph(SD.slice(v)).degree(v) == 0 \
            ....:               and len(SD.slice(v)) > 1)
            True
            sage: for v in G:
            ....:     if len(SD.slice(v)) > 1:
            ....:         xslice_seq = SD.xslice_sequence(v)
            ....:         S = G.subgraph(SD.slice(v))
            ....:         if S.degree(v) > 0:
            ....:             set(xslice_seq[1]) == set(S.neighbor_iterator(v))
            True
            True
            True
            True

        TESTS::

            sage: SD = graphs.PetersenGraph().slice_decomposition()
            sage: SD.xslice_sequence('Terry')
            Traceback (most recent call last):
            ...
            LookupError: vertex (Terry) does not appear in the slice decomposition
        """
        if v not in self.sigma_inv:
            raise LookupError(f"vertex ({v}) does not appear in the slice "
                              "decomposition")
        cdef size_t i = self.sigma_inv[v]
        return self._xslice_sequence(i)

    def lexicographic_label(self, v):
        r"""
        Return the lexicographic label of the vertex `v`.

        The lexicographic label of a vertex `v` is the list of all the
        neighbors of `v` that appear before `v` in the lexBFS ordering
        corresponding to the slice decomposition.

        INPUT:

        - ``v`` -- a vertex of the graph corresponding to the slice
          decomposition.

        OUTPUT:

        A list of vertices.

        EXAMPLES::

            sage: G = Graph('L~mpn~Nrv{^o~_').relabel('abcdefguvwxyz',inplace=False)
            sage: SD = G.slice_decomposition(initial_vertex='x')
            sage: SD.lexicographic_label('f')
            ['x', 'a', 'c', 'd']
            sage: pos = lambda v: SD.lexBFS_order().index(v)
            sage: set(SD.lexicographic_label('f')) \
            ....:   == {v for v in G.neighbors('f') if pos(v) < pos('f')}
            True

        TESTS::

            sage: SD = graphs.PetersenGraph().slice_decomposition()
            sage: SD.lexicographic_label('Eric')
            Traceback (most recent call last):
            ...
            LookupError: vertex (Eric) does not appear in the slice decomposition
        """
        if v not in self.sigma_inv:
            raise LookupError(f"vertex ({v}) does not appear in the slice "
                              "decomposition")
        cdef size_t i = self.sigma_inv[v]
        return self._xslice_lex_label(i)

    def active_edges(self, v):
        r"""
        Return the active edges of the vertex `v`.

        An edge `(u, w)` is said to be active for `v` if `u` and `w` belongs
        to two differents slices of the x-slice sequence of `v`. Note that it
        defines a partition of the edges of the underlying graph.

        INPUT:

        - ``v`` -- a vertex of the graph corresponding to the slice
          decomposition.

        OUTPUT:

        A list of edges

        EXAMPLES:

        ::

            sage: G = Graph('L~mpn~Nrv{^o~_').relabel('abcdefguvwxyz',inplace=False)
            sage: SD = G.slice_decomposition(initial_vertex='x')
            sage: SD.xslice_sequence('a')
            [['a'], ['b', 'c', 'd', 'e', 'f'], ['g']]
            sage: ('c', 'g') in SD.active_edges('a')
            True
            sage: ('a', 'c') in SD.active_edges('a')
            True
            sage: ('c', 'd') in SD.active_edges('a')  # c and d in same slice
            False
            sage: ('a', 'u') in SD.active_edges('a')  # u not in x-slice of a
            False

        The set of active edges of every vertex is a partition of the edges::

            sage: from itertools import chain
            sage: E = list(chain(*(SD.active_edges(v) for v in G)))
            sage: G.size() == len(E) == len(set(E)) \
            ....:   and all(G.has_edge(u, w) for v in G for u, w in SD.active_edges(v))
            True

        TESTS::

            sage: SD = graphs.PetersenGraph().slice_decomposition()
            sage: SD.active_edges('Graham')
            Traceback (most recent call last):
            ...
            LookupError: vertex (Graham) does not appear in the slice decomposition
        """
        if v not in self.sigma_inv:
            raise LookupError(f"vertex ({v}) does not appear in the slice "
                              "decomposition")
        cdef size_t i = self.sigma_inv[v]
        return self._xslice_active_edges(i)

    def _slice(self, size_t idx):
        r"""
        This method is for internal use only

        TESTS:

        Indirect doctests::

            sage: SD = graphs.HouseGraph().slice_decomposition()
            sage: SD.slice(1)
            [1, 2]
        """
        return list(self.sigma[idx:idx+self.xslice_len[idx]])

    def _xslice_sequence(self, size_t idx):
        r"""
        This method is for internal use only

        TESTS:

        Indirect doctests::

            sage: SD = graphs.HouseGraph().slice_decomposition()
            sage: SD.xslice_sequence(0)
            [[0], [1, 2], [3], [4]]
        """
        cdef size_t l = self.xslice_len[idx]
        cdef size_t j = idx + 1
        cdef size_t lj

        S = [ [self.sigma[idx]] ]
        while j < idx + l:
            lj = self.xslice_len[j]
            S.append(list(self.sigma[j:j+lj]))
            j += lj
        assert j == idx + l, "slice decomposition is ill-formed"
        return S

    def _xslice_lex_label(self, size_t idx):
        r"""
        This method is for internal use only

        TESTS:

        Indirect doctests::

            sage: SD = graphs.HouseGraph().slice_decomposition()
            sage: SD.lexicographic_label(3)
            [1, 2]
        """
        return list(self.lex_label[idx])

    def _xslice_active_edges(self, size_t idx):
        r"""
        This method is for internal use only

        TESTS:

        Indirect doctests::

            sage: SD = graphs.HouseGraph().slice_decomposition()
            sage: SD.active_edges(0)
            [(0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (3, 4)]
        """
        cdef size_t l = self.xslice_len[idx]
        cdef size_t llv_prefix = len(self.lex_label[idx])
        cdef size_t j = idx + 1
        cdef size_t lj

        A = []
        while j < idx + l:
            lj = self.xslice_len[j]
            llj = self.lex_label[j]
            for u in self.sigma[j:j+lj]:
                for w in llj[llv_prefix:]:
                    A.append((w, u))
            j += lj
        assert j == idx + l, "slice decomposition is ill-formed"
        return A

    def underlying_graph(self):
        r"""
        Return the underlying graph corresponding to the slice decomposition.

        If `G` was the graph given as parameter to compute the slice
        decomposition, the underlying graph corresponds to ``G.to_simple()``
        where labels are ignored, *i.e.*, it is the input graph without loops,
        multiple edges and labels.

        .. NOTE::

            This method is mostly defined to test the computation of
            lexicographic labels and actives edges.

        EXAMPLES:

        ::

            sage: G = Graph('L~mpn~Nrv{^o~_').relabel('abcdefguvwxyz',inplace=False)
            sage: SD = G.slice_decomposition(initial_vertex='x')
            sage: SD.underlying_graph() == G
            True

        The graph can have loops or multiple edges but they are ignored::

            sage: G = graphs.CubeConnectedCycle(2)  # multiple edges
            sage: SD = G.slice_decomposition()
            sage: SD.underlying_graph() == G.to_simple(immutable=True)
            True

            sage: G = graphs.CubeConnectedCycle(1)  # loops
            sage: SD = G.slice_decomposition()
            sage: SD.underlying_graph() == G.to_simple(immutable=True)
            True

        TESTS::

            sage: for _ in range(5):
            ....:   G = graphs.RandomGNP(15, 0.3)
            ....:   SD = G.slice_decomposition()
            ....:   SD.underlying_graph() == G
            True
            True
            True
            True
            True
        """
        if not hasattr(self, '_underlying_graph'):
            vertices = self.sigma
            edges = [(u, v) for i, v in enumerate(self.sigma)
                            for u in self.lex_label[i]]
            data = [vertices, edges]
            Gclass = self._graph_class
            self._underlying_graph = Gclass(data, format='vertices_and_edges',
                                            immutable=True)
        return self._underlying_graph

    def _repr_(self):
        r"""
        Return a string representation of a ``SliceDecomposition`` object.

        TESTS::

            sage: G = graphs.PetersenGraph(); SD = G.slice_decomposition()
            sage: repr(SD)
            '[0[1[4[5]]] [2[6]] [3] [9] [7] [8]]'
            sage: G = Graph('L~mpn~Nrv{^o~_').relabel('abcdefguvwxyz',inplace=False)
            sage: SD = G.slice_decomposition(initial_vertex='x'); repr(SD)
            '[x[a[b[c[d]] [e[f]]] [g]] [u[y[z]]] [v[w]]]'
        """
        def inner_repr(idx):
            l = self.xslice_len[idx]
            S = []
            if l > 1:
                j = idx + 1
                while j < idx + l:
                    lj = self.xslice_len[j]
                    S.append(inner_repr(j))
                    j += lj
                assert j == idx + l, "slice decomposition is ill-formed"
            return f'{self.sigma[idx]}' + ' '.join(f'[{s}]' for s in S)
        return f'[{inner_repr(0)}]'

    def _latex_(self):
        r"""
        Return a string to render, using `\LaTeX`, the slice decomposition as a
        tree.

        TESTS::

            sage: from sage.graphs.graph_latex import check_tkz_graph
            sage: check_tkz_graph()  # random - depends on Tex installation
            sage: G = graphs.PetersenGraph(); SD = G.slice_decomposition()
            sage: latex(SD)
            \begin{tikzpicture}
            ...
              v0 -- {l0, v1, v4, v6, v7, v8, v9};
              v1 -- {l1, v2};
              v2 -- {l2, v3};
              v3 -- {l3};
              v4 -- {l4, v5};
              v5 -- {l5};
              v6 -- {l6};
              v7 -- {l7};
              v8 -- {l8};
              v9 -- {l9};
            ...
            \end{tikzpicture}

        """
        from sage.misc.latex import latex

        latex.add_package_to_preamble_if_available("tikz")
        latex.add_to_preamble(r"\usetikzlibrary{arrows,shapes,fit}")
        latex.add_to_preamble(r"\usetikzlibrary{graphs,graphdrawing}")
        latex.add_to_preamble(r"\usegdlibrary{trees}")

        # Call latex() on all vertices
        sigma_latex = [ latex(v) for v in self.sigma ]
        slices = [[] for _ in self.sigma]

        lines = [ r"\begin{tikzpicture}" ]
        lines.append(r"\graph [tree layout,level distance=0,level sep=1em,"
                     r"sibling distance=0,sibling sep=0.6em,"
                     r"tail anchor=center,head anchor=north,"
                     r"nodes={draw,rectangle,inner xsep=0.2em},edges={thick}]")
        lines.append("{")
        bo, bc = "{", "}"  # to write { and } in f-strings
        # Create the nodes and leaves of the slice decomposition tree
        for i in range(len(self.sigma)):
            l = self.xslice_len[i]
            label = r"\ ".join(sigma_latex[i:i+l])
            lines.append(f"  v{i}[as={bo}{label}{bc}];")
            lines.append(f"  l{i}[draw=none,as={bo}{sigma_latex[i]}{bc}];")
            j = i + 1
            slices[i].append(f"l{i}")
            while j < i + l:
                slices[i].append(f"v{j}")
                j += self.xslice_len[j]
        # Create the edges of the slice decomposition tree
        for i, S in enumerate(slices):
            lines.append(f"  v{i} -- " + "{" + ", ".join(S) + "};")
        lines.append("};")
        # Add dahsed red boxes around xslices
        for i, S in enumerate(slices):
            fit=" ".join(f"({s})" for s in S)
            lines.append(rf"\node (s{i}) [rectangle,inner xsep=0.2em,draw=red,"
                         f"densely dashed,fit={fit}]{bo}{bc};")

        lines.append(r"\end{tikzpicture}")
        return "\n".join(lines)
