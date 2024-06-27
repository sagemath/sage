# cython: binding=True
r"""
Asteroidal triples

This module contains the following function:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`is_asteroidal_triple_free` | Test if the input graph is asteroidal triple-free

Definition
----------

Three independent vertices of a graph form an *asteroidal triple* if every two
of them are connected by a path avoiding the neighborhood of the third one. A
graph is *asteroidal triple-free* (*AT-free*, for short) if it contains no
asteroidal triple [LB1962]_.

Use ``graph_classes.AT_free.description()`` to get some known properties of
AT-free graphs, or visit `this page
<http://www.graphclasses.org/classes/gc_61.html>`_.


Algorithm
---------

This module implements the  *Straightforward algorithm* recalled in [Koh2004]_ and
due to [LB1962]_ for testing if a graph is AT-free or not. This algorithm has time
complexity in `O(n^3)` and space complexity in `O(n^2)`.

This algorithm uses the *connected structure* of the graph, stored into a
`n\times n` matrix `M`. This matrix is such that `M[u][v]==0` if `v\in
(\{u\}\cup N(u))`, and otherwise `M[u][v]` is the unique identifier (a strictly
positive integer) of the connected component of `G\setminus(\{u\}\cup N(u))` to
which `v` belongs. This connected structure can be computed in time `O(n(n+m))`
using `n` BFS.

Now, a triple `u, v, w\in V` is an asteroidal triple if and only if it satisfies
`M[u][v]==M[u][w]` and `M[v][u]==M[v][w]` and `M[w][u]==M[w][v]`, assuming all
these values are positive. Indeed, if `M[u][v]==M[u][w]`, `v` and `w` are in the
same connected component of `G\setminus(\{u\}\cup N(u))`, and so there is a path
between `v` and `w` avoiding the neighborhood of `u`. The algorithm iterates
over all triples.


Functions
---------
"""

# ****************************************************************************
#       Copyright (C) 2015 David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from libc.stdint cimport uint32_t
from cysignals.signals cimport sig_on, sig_off
from memory_allocator cimport MemoryAllocator

from sage.data_structures.bitset_base cimport *
from sage.graphs.base.static_sparse_graph cimport short_digraph, init_short_digraph, free_short_digraph


def is_asteroidal_triple_free(G, certificate=False):
    """
    Test if the input graph is asteroidal triple-free.

    An independent set of three vertices such that each pair is joined by a path
    that avoids the neighborhood of the third one is called an *asteroidal
    triple*. A graph is asteroidal triple-free (AT-free) if it contains no
    asteroidal triples. See the :mod:`module's documentation
    <sage.graphs.asteroidal_triples>` for more details.

    This method returns ``True`` if the graph is AT-free and ``False`` otherwise.

    INPUT:

    - ``G`` -- a Graph

    - ``certificate`` -- boolean (default: ``False``); by default, this method
      returns ``True`` if the graph is asteroidal triple-free and ``False``
      otherwise. When ``certificate==True``, this method returns in addition a
      list of three vertices forming an asteroidal triple if such a triple is
      found, and the empty list otherwise.

    EXAMPLES:

    The complete graph is AT-free, as well as its line graph::

        sage: G = graphs.CompleteGraph(5)
        sage: G.is_asteroidal_triple_free()
        True
        sage: G.is_asteroidal_triple_free(certificate=True)
        (True, [])
        sage: LG = G.line_graph()
        sage: LG.is_asteroidal_triple_free()
        True
        sage: LLG = LG.line_graph()
        sage: LLG.is_asteroidal_triple_free()
        False

    The PetersenGraph is not AT-free::

        sage: from sage.graphs.asteroidal_triples import *
        sage: G = graphs.PetersenGraph()
        sage: G.is_asteroidal_triple_free()
        False
        sage: G.is_asteroidal_triple_free(certificate=True)
        (False, [0, 2, 6])

    TESTS:

    Giving anything else than a Graph::

        sage: from sage.graphs.asteroidal_triples import is_asteroidal_triple_free
        sage: is_asteroidal_triple_free(DiGraph())
        Traceback (most recent call last):
        ...
        ValueError: The first parameter must be a Graph.
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("The first parameter must be a Graph.")

    cdef uint32_t n = <uint32_t>G.order()
    cdef int i

    # ==> Trivial cases
    if n < 3:
        return (True, []) if certificate else True

    # ==> Initialize some data structures for is_asteroidal_triple_free_C
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef uint32_t* waiting_list = <uint32_t*> mem.allocarray(n, sizeof(uint32_t))
    cdef uint32_t* _connected_structure = <uint32_t*> mem.calloc(n * n, sizeof(uint32_t))
    cdef uint32_t** connected_structure = <uint32_t**> mem.allocarray(n, sizeof(uint32_t*))

    # Copying the whole graph to obtain the list of neighbors quicker than by
    # calling out_neighbors. This data structure is well documented in the
    # module sage.graphs.base.static_sparse_graph
    cdef list int_to_vertex = list(G)
    cdef short_digraph sd
    init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_vertex,
                       sort_neighbors=False)

    cdef bitset_t seen
    bitset_init(seen, n)

    connected_structure[0] = _connected_structure
    for i in range(n-1):
        connected_structure[i+1] = connected_structure[i] + n

    cdef list ret = []

    # ==> call is_asteroidal_triple_free_C

    try:
        sig_on()
        ret = is_asteroidal_triple_free_C(n, sd, connected_structure, waiting_list, seen)
        sig_off()

    finally:
        # Release memory
        bitset_free(seen)
        free_short_digraph(sd)

    # ==> We return the result

    if certificate:
        if ret:
            return False, [int_to_vertex[i] for i in ret]
        return True, []

    return False if ret else True


cdef list is_asteroidal_triple_free_C(uint32_t n,
                                      short_digraph sd,
                                      uint32_t** connected_structure,
                                      uint32_t* waiting_list,
                                      bitset_t seen):
    """
    INPUT:

    - ``n`` -- integer; number of points in the graph

    - ``sd`` -- ``short_digraph``; a graph on ``n`` points. This data structure
      is well documented in the module
      :mod:`sage.graphs.base.static_sparse_graph`

    - ``connected_structure`` -- bidimensional array of size `n\times n` used to
      store the connected structure of the graph. All its cells must initially
      be set to 0.

    - ``waiting_list`` -- an array of size `n` to be used for BFS

    - ``seen`` -- a bitset of size `n`

    ALGORITHM:

    See the module's documentation.
    """
    cdef uint32_t waiting_beginning = 0
    cdef uint32_t waiting_end = 0
    cdef uint32_t idx_cc = 0
    cdef uint32_t source, u, v, w
    cdef uint32_t* p_tmp
    cdef uint32_t* end

    # ==> We build the connected structure

    # We run n different BFS taking each vertex as a source
    for source in range(n):

        # The source is forbidden and seen
        bitset_clear(seen)
        bitset_add(seen, source)

        # The neighbors of the source are forbidden and seen
        p_tmp = sd.neighbors[source]
        end = sd.neighbors[source + 1]
        # Iterating over all the outneighbors u of v
        while p_tmp < end:
            bitset_add(seen, p_tmp[0])
            p_tmp += 1

        # We now search for an unseen vertex
        v = bitset_first_in_complement(seen)
        while v != <uint32_t>-1:
            # and add it to the queue
            waiting_list[0] = v
            waiting_beginning = 0
            waiting_end = 0

            # We start a new connected component
            idx_cc += 1
            bitset_add(seen, v)
            connected_structure[source][v] = idx_cc

            # For as long as there are vertices left to explore in this
            # component
            while waiting_beginning <= waiting_end:

                # We pick the first one
                v = waiting_list[waiting_beginning]
                p_tmp = sd.neighbors[v]
                end = sd.neighbors[v + 1]

                # Iterating over all the outneighbors u of v
                while p_tmp < end:
                    u = p_tmp[0]

                    # If we notice one of these neighbors is not seen yet, we
                    # add it to the queue to be explored later
                    if not bitset_in(seen, u):
                        waiting_end += 1
                        waiting_list[waiting_end] = u
                        bitset_add(seen, u)
                        connected_structure[source][u] = idx_cc

                    p_tmp += 1

                waiting_beginning += 1

            # We search for a possibly unseen vertex
            v = bitset_first_in_complement(seen)

    # ==> Now that we have the component structure of the graph, we search for
    # an asteroidal triple.

    # (Possible improvement) right now, the code fixes u and tries to find v,w
    # in the same connected component of G-N[u] by going over all
    # binomial(n-1,2) pairs of point. It would be faster to:
    #
    # - Iterate on all connected components of G-N[u]
    # - Enumerate all v,w in G-N[u]
    #
    # The list of connected components of G-N[u] can be built from
    # connected_structure in O(n) time.

    for u in range(n - 2):
        for v in range(u + 1, n - 1):
            if connected_structure[u][v]:
                for w in range(v + 1, n):
                    if (connected_structure[u][v] == connected_structure[u][w] and
                            connected_structure[v][u] == connected_structure[v][w] and
                            connected_structure[w][u] == connected_structure[w][v]):
                        # We have found an asteroidal triple
                        return [u, v, w]

    # No asteroidal triple was found
    return []
