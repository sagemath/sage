r"""
Interface with Cliquer (clique-related problems)

This module defines functions based on Cliquer, an exact
branch-and-bound algorithm developed by Patric R. J. Ostergard and
written by Sampo Niskanen.

AUTHORS:

- Nathann Cohen (2009-08-14): Initial version

- Jeroen Demeyer (2011-05-06): Make cliquer interruptible (:trac:`11252`)

- Nico Van Cleemput (2013-05-27): Handle the empty graph (:trac:`14525`)

REFERENCE:

.. [NisOst2003] Sampo Niskanen and Patric R. J. Ostergard,
  "Cliquer User's  Guide, Version 1.0,"
  Communications Laboratory, Helsinki University of Technology,
  Espoo, Finland, Tech. Rep. T48, 2003.

Methods
-------
"""


#*****************************************************************************
#       Copyright (C) 2009 Nathann Cohen <nathann.cohen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


include "cysignals/signals.pxi"
include "cysignals/memory.pxi"


cdef extern from "sage/graphs/cliquer/cl.c":
     cdef int sage_clique_max(graph_t *g, int ** list)
     cdef int sage_all_clique_max(graph_t *g, int ** list)
     cdef int sage_clique_number(graph_t *g)


def max_clique(graph):
    """
    Returns the vertex set of a maximum complete subgraph.

    Currently only implemented for undirected graphs. Use
    to_undirected to convert a digraph to an undirected graph.

    EXAMPLES::

          sage: C=graphs.PetersenGraph()
          sage: max_clique(C)
          [7, 9]

    TEST::

        sage: g = Graph()
        sage: g.clique_maximum()
        []
    """
    if graph.order() == 0:
        return []

    graph,d = graph.relabel(inplace=False, return_map=True)
    d_inv = {}
    for v in d:
        d_inv[d[v]] = v

    cdef graph_t *g
    g=graph_new(graph.order())
    for e in graph.edge_iterator():
        (u,v,w)=e
        GRAPH_ADD_EDGE(g,u,v)

    cdef int* list
    cdef int size
    sig_on()
    size = sage_clique_max(g, &list)
    sig_off()
    b = []
    cdef int i
    for i in range(size):
        b.append(list[i])

    sig_free(list)
    graph_free(g)
    return list_composition(b,d_inv)


# computes all the maximum clique of a graph and return its list

def all_max_clique(graph):
    """
    Returns the vertex sets of *ALL* the maximum complete subgraphs.

    Returns the list of all maximum cliques, with each clique represented by a
    list of vertices. A clique is an induced complete subgraph, and a maximum
    clique is one of maximal order.

    .. NOTE::

       Currently only implemented for undirected graphs. Use to_undirected
       to convert a digraph to an undirected graph.

    ALGORITHM:

    This function is based on Cliquer [NisOst2003]_.

    EXAMPLES::

        sage: graphs.ChvatalGraph().cliques_maximum() # indirect doctest
        [[0, 1], [0, 4], [0, 6], [0, 9], [1, 2], [1, 5], [1, 7], [2, 3],
         [2, 6], [2, 8], [3, 4], [3, 7], [3, 9], [4, 5], [4, 8], [5, 10],
         [5, 11], [6, 10], [6, 11], [7, 8], [7, 11], [8, 10], [9, 10], [9, 11]]
        sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
        sage: G.show(figsize=[2,2])
        sage: G.cliques_maximum()
        [[0, 1, 2], [0, 1, 3]]
        sage: C=graphs.PetersenGraph()
        sage: C.cliques_maximum()
        [[0, 1], [0, 4], [0, 5], [1, 2], [1, 6], [2, 3], [2, 7], [3, 4],
         [3, 8], [4, 9], [5, 7], [5, 8], [6, 8], [6, 9], [7, 9]]
        sage: C = Graph('DJ{')
        sage: C.cliques_maximum()
        [[1, 2, 3, 4]]

    TEST::

        sage: g = Graph()
        sage: g.cliques_maximum()
        [[]]
    """
    if graph.order() == 0:
        return [[]]

    graph,d = graph.relabel(inplace=False, return_map=True)
    d_inv = {}
    for v in d:
        d_inv[d[v]] = v

    cdef graph_t *g
    g=graph_new(graph.order())

    for e in graph.edge_iterator():
        (u,v,w)=e
        GRAPH_ADD_EDGE(g,u,v)

    cdef int* list
    cdef int size
    sig_on()
    size = sage_all_clique_max(g, &list)
    sig_off()
    b = []
    c=[]
    cdef int i
    for i in range(size):
        if(list[i]!=-1):
            c.append(list[i])
        else:
            b.append(list_composition(c,d_inv))
            c=[]

    sig_free(list)
    graph_free(g)

    return sorted(b)


#computes the clique number of a graph

def clique_number(graph):
    """
    Returns the size of the largest clique of the graph (clique
    number).

    Currently only implemented for undirected graphs. Use
    to_undirected to convert a digraph to an undirected graph.

    EXAMPLES::

        sage: C = Graph('DJ{')
        sage: clique_number(C)
        4
        sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
        sage: G.show(figsize=[2,2])
        sage: clique_number(G)
        3

    TEST::

        sage: g = Graph()
        sage: g.clique_number()
        0
    """
    if graph.order() == 0:
        return 0

    graph=graph.relabel(inplace=False)
    cdef graph_t *g
    g=graph_new(graph.order())

    for e in graph.edge_iterator():
        (u,v,w)=e
        GRAPH_ADD_EDGE(g,u,v)

    cdef int c
    sig_on()
    c = sage_clique_number(g)
    graph_free(g)
    sig_off()
    return c


def list_composition(a,b):
    """
    Composes a list ``a`` with a map ``b``.

    EXAMPLES::

        sage: from sage.graphs.cliquer import list_composition
        sage: list_composition([1,3,'a'], {'a':'b', 1:2, 2:3, 3:4})
        [2, 4, 'b']

    """
    value=[]
    for i in a:
        value.append(b[i])
    return value


