# cython: binding=True
r"""
Interface with Cliquer (clique-related problems)

This module defines functions based on Cliquer, an exact
branch-and-bound algorithm developed by Patric R. J. Ostergard and
written by Sampo Niskanen.

AUTHORS:

- Nathann Cohen (2009-08-14): Initial version

- Jeroen Demeyer (2011-05-06): Make cliquer interruptible (:issue:`11252`)

- Nico Van Cleemput (2013-05-27): Handle the empty graph (:issue:`14525`)

REFERENCE:

[NO2003]_

Methods
-------
"""

# ****************************************************************************
#       Copyright (C) 2009 Nathann Cohen <nathann.cohen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.memory cimport sig_free
from cysignals.signals cimport sig_on, sig_off


cdef extern from "cliquer/cl.c":
    cdef int sage_clique_max(graph_t *g, int ** list_of_vertices)
    cdef int sage_all_clique_max(graph_t *g, int ** list_of_vertices)
    cdef int sage_clique_number(graph_t *g)
    cdef int sage_find_all_clique(graph_t *g, int ** list_of_vertices, int min_size, int max_size)


def max_clique(graph):
    """
    Return the vertex set of a maximum complete subgraph.

    .. NOTE::

        Currently only implemented for undirected graphs. Use
        :meth:`~sage.graphs.digraph.DiGraph.to_undirected` to convert a digraph
        to an undirected graph.

    EXAMPLES::

          sage: C = graphs.PetersenGraph()
          sage: from sage.graphs.cliquer import max_clique
          sage: max_clique(C)
          [7, 9]

    TESTS::

        sage: g = Graph()
        sage: g.clique_maximum()
        []
    """
    if not graph.order():
        return []

    cdef int i
    cdef list int_to_vertex = list(graph)
    cdef dict vertex_to_int = {v: i for i, v in enumerate(int_to_vertex)}

    cdef graph_t* g = graph_new(graph.order())
    for u, v in graph.edge_iterator(labels=None):
        GRAPH_ADD_EDGE(g, vertex_to_int[u], vertex_to_int[v])

    cdef int* list_of_vertices
    cdef int size
    sig_on()
    size = sage_clique_max(g, &list_of_vertices)
    sig_off()
    cdef list b = [int_to_vertex[list_of_vertices[i]] for i in range(size)]

    sig_free(list_of_vertices)
    graph_free(g)
    return b


# computes all the maximum clique of a graph and return its list

def all_max_clique(graph):
    """
    Return the vertex sets of *ALL* the maximum complete subgraphs.

    Returns the list of all maximum cliques, with each clique represented by a
    list of vertices. A clique is an induced complete subgraph, and a maximum
    clique is one of maximal order.

    .. NOTE::

        Currently only implemented for undirected graphs. Use
        :meth:`~sage.graphs.digraph.DiGraph.to_undirected` to convert a digraph
        to an undirected graph.

    ALGORITHM:

    This function is based on Cliquer [NO2003]_.

    EXAMPLES::

        sage: graphs.ChvatalGraph().cliques_maximum() # indirect doctest
        [[0, 1], [0, 4], [0, 6], [0, 9], [1, 2], [1, 5], [1, 7], [2, 3],
         [2, 6], [2, 8], [3, 4], [3, 7], [3, 9], [4, 5], [4, 8], [5, 10],
         [5, 11], [6, 10], [6, 11], [7, 8], [7, 11], [8, 10], [9, 10], [9, 11]]
        sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
        sage: G.show(figsize=[2,2])                                                     # needs sage.plot
        sage: G.cliques_maximum()
        [[0, 1, 2], [0, 1, 3]]
        sage: C = graphs.PetersenGraph()
        sage: C.cliques_maximum()
        [[0, 1], [0, 4], [0, 5], [1, 2], [1, 6], [2, 3], [2, 7], [3, 4],
         [3, 8], [4, 9], [5, 7], [5, 8], [6, 8], [6, 9], [7, 9]]
        sage: C = Graph('DJ{')
        sage: C.cliques_maximum()
        [[1, 2, 3, 4]]

    TESTS::

        sage: g = Graph()
        sage: g.cliques_maximum()
        [[]]
    """
    if not graph.order():
        return [[]]

    cdef int i
    cdef list int_to_vertex = list(graph)
    cdef dict vertex_to_int = {v: i for i, v in enumerate(int_to_vertex)}

    cdef graph_t* g = graph_new(graph.order())
    for u, v in graph.edge_iterator(labels=None):
        GRAPH_ADD_EDGE(g, vertex_to_int[u], vertex_to_int[v])

    cdef int* list_of_vertices
    cdef int size
    sig_on()
    size = sage_all_clique_max(g, &list_of_vertices)
    sig_off()
    cdef list b = []
    cdef list c = []
    for i in range(size):
        if list_of_vertices[i] != -1:
            c.append(int_to_vertex[list_of_vertices[i]])
        else:
            b.append(c)
            c = []

    sig_free(list_of_vertices)
    graph_free(g)

    return sorted(b)


def all_cliques(graph, min_size=0, max_size=0):
    r"""
    Iterator over the cliques in ``graph``.

    A clique is an induced complete subgraph. This method is an iterator over
    all the cliques with size in between ``min_size`` and ``max_size``. By
    default, this method returns only maximum cliques. Each yielded clique is
    represented by a list of vertices.

    .. NOTE::

        Currently only implemented for undirected graphs. Use
        :meth:`~sage.graphs.digraph.DiGraph.to_undirected` to convert a digraph
        to an undirected graph.

    INPUT:

    - ``min_size`` -- integer (default: 0); minimum size of reported cliques.
      When set to 0 (default), this method searches for maximum cliques. In such
      case, parameter ``max_size`` must also be set to 0.

    - ``max_size`` -- integer (default: 0); maximum size of reported cliques.
      When set to 0 (default), the maximum size of the cliques is unbounded.
      When ``min_size`` is set to 0, this parameter must be set to 0.

    ALGORITHM:

    This function is based on Cliquer [NO2003]_.

    EXAMPLES::

        sage: G = graphs.CompleteGraph(5)
        sage: list(sage.graphs.cliquer.all_cliques(G))
        [[0, 1, 2, 3, 4]]
        sage: list(sage.graphs.cliquer.all_cliques(G, 2, 3))
        [[3, 4],
         [2, 3],
         [2, 3, 4],
         [2, 4],
         [1, 2],
         [1, 2, 3],
         [1, 2, 4],
         [1, 3],
         [1, 3, 4],
         [1, 4],
         [0, 1],
         [0, 1, 2],
         [0, 1, 3],
         [0, 1, 4],
         [0, 2],
         [0, 2, 3],
         [0, 2, 4],
         [0, 3],
         [0, 3, 4],
         [0, 4]]
        sage: G.delete_edge([1,3])
        sage: list(sage.graphs.cliquer.all_cliques(G))
        [[0, 2, 3, 4], [0, 1, 2, 4]]

    TESTS::

        sage: G = graphs.CompleteGraph(3)
        sage: list(sage.graphs.cliquer.all_cliques(G, 2, 3))
        [[1, 2], [0, 1], [0, 1, 2], [0, 2]]

        sage: G = graphs.EmptyGraph()
        sage: list(sage.graphs.cliquer.all_cliques(G, 2, 3))
        []

        sage: G = Graph([(0, 1), (0, 1)], multiedges=True)
        sage: list(sage.graphs.cliquer.all_cliques(G, 2, 2))
        [[0, 1]]

        sage: list(sage.graphs.cliquer.all_cliques(G, 0, 2))
        Traceback (most recent call last):
        ...
        ValueError: max_size > 0 is incompatible with min_size == 0

    .. TODO::

        Use the re-entrant functionality of Cliquer [NO2003]_ to avoid storing
        all cliques.
    """
    if not min_size and max_size > 0:
        raise ValueError("max_size > 0 is incompatible with min_size == 0")
    if not graph:
        return

    cdef int i
    cdef list int_to_vertex = list(graph)
    cdef dict vertex_to_int = {v: i for i, v in enumerate(int_to_vertex)}

    cdef graph_t* g = graph_new(graph.order())
    for u, v in graph.edge_iterator(labels=None):
        GRAPH_ADD_EDGE(g, vertex_to_int[u], vertex_to_int[v])

    cdef int* list_of_vertices
    cdef int size = 0
    cdef list c
    try:
        try:
            sig_on()
            size = sage_find_all_clique(g, &list_of_vertices, min_size, max_size)
            sig_off()
        finally:
            graph_free(g)
        c = []
        for i in range(size):
            if list_of_vertices[i] != -1:
                c.append(int_to_vertex[list_of_vertices[i]])
            else:
                yield c
                c = []
    finally:
        if list_of_vertices:
            # We free ``list_of_vertices``,
            # but only if previous computations weren't interrupted before
            # allocating memory for ``list_of_vertices``.
            sig_free(list_of_vertices)


# computes the clique number of a graph

def clique_number(graph):
    """
    Return the size of the largest clique of the graph (clique number).

    .. NOTE::

        Currently only implemented for undirected graphs. Use
        :meth:`~sage.graphs.digraph.DiGraph.to_undirected` to convert a digraph
        to an undirected graph.

    EXAMPLES::

        sage: C = Graph('DJ{')
        sage: C.clique_number()
        4
        sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
        sage: G.show(figsize=[2,2])                                                     # needs sage.plot
        sage: G.clique_number()
        3

    TESTS::

        sage: g = Graph()
        sage: g.clique_number()
        0
    """
    if not graph.order():
        return 0

    cdef int i
    cdef dict vertex_to_int = {v: i for i, v in enumerate(graph)}

    cdef graph_t* g = graph_new(graph.order())
    for u, v in graph.edge_iterator(labels=None):
        GRAPH_ADD_EDGE(g, vertex_to_int[u], vertex_to_int[v])

    cdef int c
    sig_on()
    c = sage_clique_number(g)
    graph_free(g)
    sig_off()
    return c
