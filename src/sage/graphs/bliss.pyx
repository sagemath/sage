r"""
Interface with bliss: graph (iso/auto)morphism

Implemented functions:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`automorphism_group` | Returns the automorphism group of the given (di)graph
    :meth:`canonical_form` | Computes a canonical certificate for the given (di) graph.

AUTHORS:

    - Jernej Azarija

"""

#*****************************************************************************
#       Copyright (C) 2015 Jernej Azarija
#       Copyright (C) 2015 Nathann Cohen <nathann.cohen@gail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from cpython cimport PyObject
from libc.limits cimport LONG_MAX

from sage.misc.decorators import rename_keyword

cdef extern from "bliss/graph.hh" namespace "bliss":

    cdef cppclass Stats:
        pass

    cdef cppclass AbstractGraph:
        pass

    cdef cppclass Graph(AbstractGraph):
        Graph(const unsigned int)
        void add_edge(const unsigned int, const unsigned int)
        void find_automorphisms(Stats&, void (*)(void* , unsigned int,
                    const unsigned int*), void*)
        void change_color(const unsigned int, const unsigned int);
        const unsigned int* canonical_form(Stats&, void (*)(void*,unsigned int,
                    const unsigned int*), void*)

    cdef cppclass Digraph(AbstractGraph):

        Digraph(const unsigned int)
        void add_edge(const unsigned int, const unsigned int)
        void find_automorphisms(Stats&, void (*)(void* , unsigned int,
                    const unsigned int*), void*)
        void change_color(const unsigned int, const unsigned int);
        const unsigned int* canonical_form(Stats&, void (*)(void*,unsigned int,
                    const unsigned int*), void*)
        unsigned int get_hash()


cdef void add_gen(void *user_param, unsigned int n, const unsigned int *aut):
    r"""
    This function is called each time a new generator of the automorphism group
    is found.

    This function is used to append the new generators to a Python list. Its
    main job is to translate a permutation into dijoint cycles.

    INPUT:

    - ``user_param`` (``void *``) -- in the current implementation, it points
      toward a Python object which is a pair
      ``(list_of_current_generators,vert_to_integer_labelling)``.

    - ``n`` (int) -- number of points in the graph.

    - ``aut`` (int *) -- an automorphism of the graph.
    """
    cdef int tmp     = 0
    cdef int marker  = 0
    cdef int cur     = 0
    perm        = []
    done        = [False]*n

    gens, int_to_vertex = <object> <PyObject *> user_param

    while True:
        while cur < n and done[cur]:
            cur+=1
        if cur == n:
            break

        marker = tmp = cur
        cycle = [int_to_vertex[cur]]
        done[cur] = True

        while aut[tmp] != marker:
            tmp = aut[tmp]
            done[tmp] = True
            cycle.append(int_to_vertex[tmp])

        perm.append(tuple(cycle))
    gens.append(perm)

cdef Graph *bliss_graph(G, partition, vert2int, int2vert):
    r"""
    Return a bliss copy of a graph G

    INPUT:

    - ``G`` (a graph)

    - ``partition`` -- a partition of the vertex set.

    - ``vert2int, int2vert`` -- Two empty dictionaries. The entries of the
      dictionary are later set to record the labeling of our graph. They are
      taken as arguments to avoid technicalities of returning Python objects in
      Cython functions.
    """
    cdef Graph *g = new Graph(G.order())

    if g == NULL:
        raise MemoryError("Allocation Failed")

    for i,v in enumerate(G.vertices()):
        vert2int[v] = i
        int2vert[i] = v

    for x,y in G.edges(labels=False):
       g.add_edge(vert2int[x],vert2int[y])

    if partition:
        for i in xrange(1, len(partition)):
            for v in partition[i]:
                g.change_color(vert2int[v], i)
    return g

cdef Digraph *bliss_digraph(G, partition, vert2int, int2vert):
    r"""
    Return a bliss copy of a digraph G

    INPUT:

    - ``G`` (a digraph)

    - ``partition`` -- a partition of the vertex set.

    - ``vert2int, int2vert`` -- Two empty dictionaries. The entries of the
      dictionary are later set to record the labeling of our graph. They are
      taken as arguments to avoid technicalities of returning Python objects in
      Cython functions.
    """
    cdef Digraph *g = new Digraph(G.order())

    for i,v in enumerate(G.vertices()):
        vert2int[v] = i
        int2vert[i] = v

    if g == NULL:
        raise MemoryError("Allocation Failed")

    for x,y in G.edges(labels=False):
        g.add_edge(vert2int[x],vert2int[y])

    if partition:
        for i in xrange(1, len(partition)):
            for v in partition[i]:
                g.change_color(vert2int[v], i)
    return g

def automorphism_group(G, partition=None):
    """
    Computes the automorphism group of ``G`` subject to the coloring ``partition.``

    INPUT:

    - ``G`` -- A graph

    - ``partition`` -- A partition of the vertices of ``G`` into color classes.
      Defaults to ``None``, which is equivalent to a partition of size 1.

    TESTS::

        sage: from sage.graphs.bliss import automorphism_group                  # optional - bliss
        sage: G = graphs.PetersenGraph()                                        # optional - bliss
        sage: automorphism_group(G).is_isomorphic(G.automorphism_group())       # optional - bliss
        True

        sage: G = graphs.HeawoodGraph()                                         # optional - bliss
        sage: p = G.bipartite_sets()                                            # optional - bliss
        sage: A = G.automorphism_group(partition=[list(p[0]), list(p[1])])      # optional - bliss
        sage: automorphism_group(G, partition=p).is_isomorphic(A)               # optional - bliss
        True
    """

    cv = 0
    n = G.order()
    vert2int = {}
    int2vert = {}

    cdef Graph *g = NULL
    cdef Digraph *d = NULL
    cdef Stats s

    gens = []
    data = (gens, int2vert)

    if G.is_directed():
        d = bliss_digraph(G, partition, vert2int, int2vert)
        d.find_automorphisms(s, add_gen, <PyObject *> data)
        del d
    else:
        g = bliss_graph(G, partition, vert2int, int2vert)
        g.find_automorphisms(s, add_gen, <PyObject *> data)
        del g

    from sage.groups.perm_gps.permgroup import PermutationGroup
    return PermutationGroup(gens,domain=G)

cdef void empty_hook(void *user_param , unsigned int n, const unsigned int *aut):
    return

@rename_keyword(deprecation=21111, certify='certificate')
def canonical_form(G, partition=None, return_graph=False, certificate=False):
    """
    Return a canonical label of ``G``

    A canonical label ``canonical_form(G)`` of ``G`` is a (di)graph defined on
    `\{0,...,n-1\}` such that ``G`` is isomorphic to ``H`` if and only if
    ``canonical_form(G)`` is equal to ``canonical_form(H)``.

    INPUT:

    - ``G`` -- A graph or digraph.

    - ``partition`` -- A partition of the vertices of ``G`` into color classes.
        Defaults to ``None``.

    - ``return_graph`` -- If set to ``True``, ``canonical_form`` returns the
        canonical graph of G. Otherwise, it returns its set of edges.

    - ``certificate`` -- If set to ``True`` returns the labeling of G into a
      canonical graph.

    TESTS::

        sage: from sage.graphs.bliss import canonical_form                  # optional - bliss
        sage: G = graphs.PetersenGraph()                                    # optional - bliss
        sage: canonical_form(G)                                             # optional - bliss
        [(2, 0), (2, 1), (3, 0), (4, 1), (5, 3), (5, 4), (6, 0), (6, 4),
         (7, 1), (7, 3), (8, 2), (8, 5), (9, 6), (9, 7), (9, 8)]

        sage: P = graphs.GeneralizedPetersenGraph(5,2)                      # optional - bliss
        sage: Q = graphs.PetersenGraph()                                    # optional - bliss
        sage: canonical_form(P) == canonical_form(Q)                        # optional - bliss
        True

        sage: canonical_form(Graph(15),return_graph=True)                   # optional - bliss
        Graph on 15 vertices
        sage: g = digraphs.RandomTournament(40)                             # optional - bliss
        sage: g.is_isomorphic(canonical_form(g,return_graph=True))          # optional - bliss
        True

        sage: g1 = graphs.RandomGNP(100,.4)                                 # optional - bliss
        sage: r = Permutations(range(100)).random_element()                 # optional - bliss
        sage: g2 = Graph([(r[u],r[v]) for u,v in g1.edges(labels=False)])   # optional - bliss
        sage: g1 = canonical_form(g1,return_graph=True)                     # optional - bliss
        sage: g2 = canonical_form(g2,return_graph=True)                     # optional - bliss
        sage: g2 == g2                                                      # optional - bliss
        True

        sage: g = Graph({1: [2]})
        sage: g_ = canonical_form(g, return_graph=True, certify=True)    # optional - bliss
        doctest...: DeprecationWarning: use the option 'certificate' instead of 'certify'
        See http://trac.sagemath.org/21111 for details.
        sage: 0 in g_[0]                                                 # optional - bliss
        True
    """
    # We need this to convert the numbers from <unsigned int> to
    # <long>. This assertion should be true simply for memory reasons.
    assert <unsigned long>(G.order()) <= <unsigned long>LONG_MAX

    cdef const unsigned int* aut
    cdef Graph* g
    cdef Digraph* d
    cdef Stats s
    cdef dict relabel

    cdef list edges = []
    cdef long e, f

    vert2int = {}

    if G.is_directed():
        d = bliss_digraph(G, partition, vert2int, {})
        aut = d.canonical_form(s, empty_hook, NULL)
        for x,y in G.edges(labels=False):
            e,f = aut[ vert2int[x] ], aut[ vert2int[y] ]
            edges.append( (e,f) )
        relabel = {v: <long>aut[vert2int[v]] for v in G}
        del d
    else:
        g = bliss_graph(G, partition, vert2int, {})
        aut = g.canonical_form(s, empty_hook, NULL)
        for x,y in G.edges(labels=False):
            e,f = aut[ vert2int[x] ], aut[ vert2int[y] ]
            edges.append( (e,f) if e > f else (f,e))
        relabel = {v: <long>aut[vert2int[v]] for v in G}
        del g

    if return_graph:
        if G.is_directed():
            from sage.graphs.graph import DiGraph
            G = DiGraph(edges,loops=G.allows_loops(),multiedges=G.allows_multiple_edges())
        else:
            from sage.graphs.graph import Graph
            G = Graph(edges,loops=G.allows_loops(),multiedges=G.allows_multiple_edges())

        G.add_vertices(vert2int.values())
        return (G, relabel) if certificate else G

    if certificate:
        return sorted(edges),relabel

    return sorted(edges)

