# sage_setup: distribution = sagemath-graphs
# cython: binding=True
# distutils: libraries = rw
r"""
Rank Decompositions of graphs

This module wraps C code from Philipp Klaus Krause computing an optimal
rank-decomposition.

**Definitions :**

Given a graph `G` and a subset `S\subseteq V(G)` of vertices, the *rank-width*
of `S` in `G`, denoted `rw_G(S)`, is equal to the rank in `GF(2)` of the `|S|
\times (|V|-|S|)` matrix of the adjacencies between the vertices of `S` and
`V\backslash S`. By definition, `rw_G(S)` is equal to `rw_G(\overline S)` where
`\overline S` is the complement of `S` in `V(G)`.

A *rank-decomposition* of `G` is a tree whose `n` leaves are the elements of
`V(G)`, and whose internal nodes have degree 3. In a tree, any edge naturally
corresponds to a bipartition of the vertex set : indeed, the removal of any edge
splits the tree into two connected components, thus splitting the set of leaves
(i.e. vertices of `G`) into two sets. Hence we can define for any edge `e\in
E(G)` a width equal to the value `rw_G(S)` or `rw_G(\overline S)`, where
`S,\overline S` is the bipartition obtained from `e`. The *rank-width*
associated to the whole decomposition is then set to the maximum of the width of
all the edges it contains.

A *rank-decomposition* is said to be optimal for `G` if it is the decomposition
achieving the minimal *rank-width*.

**RW -- The original source code :**

RW is a program that calculates rank-width and
rank-decompositions. It is based on ideas from:

    * "Computing rank-width exactly" by Sang-il Oum [Oum2009]_
    * "Sopra una formula numerica" by Ernesto Pascal
    * "Generation of a Vector from the Lexicographical Index" by B.P. Buckles
      and M. Lybanon [BL1977]_
    * "Fast additions on masked integers" by Michael D. Adams and David S. Wise
      [AW2006]_

**OUTPUT:**

The rank decomposition is returned as a tree whose vertices are subsets of
`V(G)`. Its leaves, corresponding to the vertices of `G` are sets of 1 elements,
i.e. singletons.

::

    sage: g = graphs.PetersenGraph()
    sage: rw, tree = g.rank_decomposition()
    sage: all(len(v)==1 for v in tree if tree.degree(v) == 1)
    True

The internal nodes are sets of the decomposition. This way, it is easy to deduce
the bipartition associated to an edge from the tree. Indeed, two adjacent
vertices of the tree are comparable sets : they yield the bipartition obtained
from the smaller of the two and its complement.

::

    sage: g = graphs.PetersenGraph()
    sage: rw, tree = g.rank_decomposition()
    sage: u = Set([8, 9, 3, 7])
    sage: v = Set([8, 9])
    sage: tree.has_edge(u,v)
    True
    sage: m = min(u,v)
    sage: bipartition = (m, Set(g.vertices(sort=False)) - m)
    sage: bipartition
    ({8, 9}, {0, 1, 2, 3, 4, 5, 6, 7})

.. WARNING::

    * The current implementation cannot handle graphs of `\geq 32` vertices.

    * A bug that has been reported upstream make the code crash immediately on
      instances of size `30`. If you experience this kind of bug please report
      it to us, what we need is some information on the hardware you run to know
      where it comes from !

EXAMPLES::

        sage: g = graphs.PetersenGraph()
        sage: g.rank_decomposition()
        (3, Graph on 19 vertices)

AUTHORS:

- Philipp Klaus Krause : Implementation of the C algorithm
- Nathann Cohen : Interface with Sage and documentation

Methods
-------
"""

# ****************************************************************************
#       Copyright (C) 2011 Nathann Cohen <nathann.cohen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.memory cimport check_allocarray, sig_free
from cysignals.signals cimport *

from libc.string cimport memset

cdef list id_to_vertices
cdef dict vertices_to_id


def rank_decomposition(G, verbose=False):
    r"""
    Compute an optimal rank-decomposition of the given graph.

    This function is available as a method of the :class:`Graph
    <sage.graphs.graph>` class. See :meth:`rank_decomposition
    <sage.graphs.graph.Graph.rank_decomposition>`.

    INPUT:

    - ``verbose`` -- boolean (default: ``False``); whether to display progress
      information while computing the decomposition

    OUTPUT:

    A pair ``(rankwidth, decomposition_tree)``, where ``rankwidth`` is a
    numerical value and ``decomposition_tree`` is a ternary tree describing the
    decomposition (cf. the module's documentation).

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.rankwidth import rank_decomposition
        sage: g = graphs.PetersenGraph()
        sage: rank_decomposition(g)
        (3, Graph on 19 vertices)

    On more than 32 vertices::

        sage: g = graphs.RandomGNP(40, .5)
        sage: rank_decomposition(g)
        Traceback (most recent call last):
        ...
        RuntimeError: the rank decomposition cannot be computed on graphs of >= 32 vertices

    The empty graph::

        sage: g = Graph()
        sage: rank_decomposition(g)
        (0, Graph on 0 vertices)
    """
    cdef int n = G.order()

    if n >= 32:
        raise RuntimeError("the rank decomposition cannot be computed "
                           "on graphs of >= 32 vertices")

    elif not n:
        from sage.graphs.graph import Graph
        return (0, Graph())

    cdef int i

    if sage_graph_to_matrix(G):
        raise RuntimeError("there has been a mistake while converting the Sage "
                           "graph to a C structure, the memory is probably "
                           "insufficient (2^(n+1) is a *LOT*)")

    for i in range(n + 1):

        if verbose:
            print("Calculating for subsets of size ", i, "/", n + 1)

        # We want to properly deal with exceptions, in particular
        # KeyboardInterrupt. Whatever happens, when this code fails the memory
        # *HAS* to be freed, as it is particularly greedy (a graph of 29
        # vertices already requires more than 1GB of memory).

        if not sig_on_no_except():
            destroy_rw()
        cython_check_exception()

        # Actual computation
        calculate_level(i)
        sig_off()

    cdef int rank_width = <int> get_rw()

    # Original way of displaying the decomposition
    # print_rank_dec(0x7ffffffful >> (31 - num_vertices), 0)
    g = mkgraph(n)

    # Free the memory
    destroy_rw()

    return (rank_width, g)


cdef int sage_graph_to_matrix(G) noexcept:
    r"""
    Convert the given Sage graph as an adjacency matrix.
    """
    global id_to_vertices
    global vertices_to_id
    global adjacency_matrix
    global cslots

    cdef int num_vertices = G.order()

    # Prepares the C structure for the computation
    if init_rw_dec(num_vertices):
        # The original code does not completely frees the memory in case of
        # error
        destroy_rw()
        return 1

    memset(adjacency_matrix, 0, sizeof(subset_t) * num_vertices)

    # Initializing the lists of vertices
    cdef int i
    id_to_vertices = list(G)
    vertices_to_id = {v: i for i, v in enumerate(id_to_vertices)}

    # Filling the matrix
    for u, v in G.edge_iterator(labels=False):
        if u == v:
            continue
        set_am(vertices_to_id[u], vertices_to_id[v], 1)

    # All is fine.
    return 0


cdef uint_fast32_t bitmask(int i) noexcept:
    return (1ul << i)


cdef void set_am(int i, int j, int val) noexcept:
    r"""
    Set/Unset an arc between vertices i and j

    (this function is a copy of what can be found in rankwidth/rw.c)
    """
    global adjacency_matrix

    adjacency_matrix[i] &= ~bitmask(j)
    adjacency_matrix[j] &= ~bitmask(i)

    if val:
        adjacency_matrix[i] |= bitmask(j)
        adjacency_matrix[j] |= bitmask(i)


cdef void print_rank_dec(subset_t s, int l) noexcept:
    r"""
    Print the current rank decomposition as a text

    This function is a copy of the C routine printing the rank-decomposition is
    the original source code. It s not used at the moment, but can still prove
    useful when debugging the code.
    """
    global cslots

    print('\t' * l, end="")

    print("cslot: ", <unsigned int> s)
    if not cslots[s]:
        return
    print_rank_dec(cslots[s], l + 1)
    print_rank_dec(s & ~cslots[s], l + 1)


def mkgraph(int num_vertices):
    r"""
    Return the graph corresponding to the current rank-decomposition.

    (This function is for internal use)

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.rankwidth import rank_decomposition
        sage: g = graphs.PetersenGraph()
        sage: rank_decomposition(g)
        (3, Graph on 19 vertices)
    """
    global cslots
    global id_to_vertices

    cdef subset_t s

    from sage.graphs.graph import Graph
    g = Graph()

    cdef subset_t * tab = <subset_t *>check_allocarray(2 * num_vertices - 1, sizeof(subset_t))
    tab[0] = 0x7ffffffful >> (31 - num_vertices)

    cdef int beg = 0
    cdef int end = 1

    while beg != end:

        s = tab[beg]
        beg += 1

        if not cslots[s]:
            continue

        g.add_edge(bitset_to_vertex_set(s), bitset_to_vertex_set(s & ~cslots[s]))
        g.add_edge(bitset_to_vertex_set(s), bitset_to_vertex_set(cslots[s]))

        tab[end] = s & ~cslots[s]
        end += 1
        tab[end] = cslots[s]
        end += 1

    sig_free(tab)
    return g


cdef bitset_to_vertex_set(subset_t s):
    """
    Return as a Set object the set corresponding to the given subset_t
    variable.
    """
    from sage.rings.integer import Integer
    from sage.sets.set import Set
    from sage.data_structures.bitset import FrozenBitset
    return Set(list(FrozenBitset((Integer(<unsigned int> s).binary())[::-1])))
