r"""
Weakly chordal graphs

This module deals with everything related to weakly chordal graphs. It currently
contains the following functions:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.graphs.weakly_chordal.is_long_hole_free` | Tests whether ``g`` contains an induced cycle of length at least 5.
    :meth:`~sage.graphs.weakly_chordal.is_long_antihole_free` | Tests whether ``g`` contains an induced anticycle of length at least 5.
    :meth:`~sage.graphs.weakly_chordal.is_weakly_chordal` | Tests whether ``g`` is weakly chordal.

Author:

- Birk Eisermann (initial implementation)
- Nathann Cohen (some doc and optimization)

REFERENCES:

.. [NikolopoulosPalios07] Nikolopoulos, S.D. and Palios, L.
  Detecting holes and antiholes in graphs
  Algorithmica, 2007
  Vol. 47, number 2, pages 119--138
  http://www.cs.uoi.gr/~stavros/C-Papers/C-2004-SODA.pdf



Methods
-------
"""

##############################################################################
#       Copyright (C) 2012 Birk Eisermann <eisermbi@fastmail.fm>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

include "sage/misc/bitset.pxi"

cdef inline int has_edge(bitset_t bs, int u, int v, int n):
    return bitset_in(bs, u*n+v)


def is_long_hole_free(g, certificate=False):
    r"""
    Tests whether ``g`` contains an induced cycle of length at least 5.

    INPUT:

    - ``certificate`` -- boolean (default: ``False``)

      Whether to return a certificate. When ``certificate = True``, then
      the function returns

      * ``(True, [])`` if ``g`` does not contain such a cycle.
        For this case, it is not known how to provide a certificate.
      * ``(False, Hole)`` if ``g`` contains an induced cycle of length at
        least 5. ``Hole`` returns this cycle.

      If ``certificate = False``, the function returns just ``True`` or
      ``False`` accordingly.

    ALGORITHM:

    This algorithm tries to find a cycle in the graph of all induced `P_4` of
    `g`, where two copies `P` and `P'` of `P_4` are adjacent if there exists a
    (not necessarily induced) copy of `P_5=u_1u_2u_3u_4u_5` such that
    `P=u_1u_2u_3u_4` and `P'=u_2u_3u_4u_5`.

    This is done through a depth-first-search. For efficiency, the auxiliary
    graph is constructed on-the-fly and never stored in memory.

    The run time of this algorithm is `O(m^2)` [NikolopoulosPalios07]_ ( where
    `m` is the number of edges of the graph ) .

    EXAMPLES:

    The Petersen Graph contains a hole::

        sage: g = graphs.PetersenGraph()
        sage: g.is_long_hole_free()
        False

    The following graph contains a hole, which we want to display::

        sage: g = graphs.FlowerSnark()
        sage: r,h = g.is_long_hole_free(certificate=True)
        sage: r
        False
        sage: Graph(h).is_isomorphic(graphs.CycleGraph(h.order()))
        True

    TESTS:

    Another graph with vertices 2, ..., 8, 10::

        sage: g = Graph({2:[3,8],3:[2,4],4:[3,8,10],5:[6,10],6:[5,7],7:[6,8],8:[2,4,7,10],10:[4,5,8]})
        sage: r,hole = g.is_long_hole_free(certificate=True)
        sage: r
        False
        sage: hole
        Subgraph of (): Graph on 5 vertices
        sage: hole.is_isomorphic(graphs.CycleGraph(hole.order()))
        True
    """
    cdef int a,b,c,i,u,v,d

    # relabel the graph on 0...n-1
    cdef dict label_id = g.relabel(return_map = True)
    cdef dict id_label = {idd:label for label, idd in label_id.iteritems()}

    # A dense copy of our graph
    cdef bitset_t dense_graph
    cdef int n = g.order()
    bitset_init(dense_graph, n*n)
    bitset_set_first_n(dense_graph, 0)
    for u,v in g.edges(labels = False):
        bitset_add(dense_graph,u*n+v)
        bitset_add(dense_graph,v*n+u)

    InPath = {} #vertices of the current path with their position (InPath[v] = i)
    VisitedP3 = {} #stores triples (u,v,w) which represent visited paths of length 3


    def process(a,b,c,i):
        InPath[c] = i  # c is the (i+1)-th vertex at position i
        VisitedP3[a,b,c] = True
        VisitedP3[c,b,a] = True

        for d in g.neighbor_iterator(c):
            if not has_edge(dense_graph,d,a,n) and not has_edge(dense_graph,d,b,n):
                # a-b-c-d form an induced path P_4

                if d in InPath:
                    # d is already contained in InPath
                    # HOLE FOUND !!!
                    if certificate:
                        j = InPath[d]
                        C = [v for v,vj in InPath.iteritems() if vj >= j]
                        C.sort(key = lambda x: InPath[x])
                        C_index = {u:i for i,u in enumerate(C)}

                        # At this step C[0]C[1]..... is a cycle such that any 4
                        # consecutive vertices induce a P4. C may not be an
                        # induced cycle, so we extract one from it.

                        # To do so, we look for the *shortest* edge C[i]C[j]
                        # between two nonconsecutive vertices of C, where the
                        # length is the difference |i-j|.
                        #
                        # C[i]...C[j] is necessarily an induced cycle.⇧

                        gg = g.subgraph(C)
                        gg.delete_edges(zip(C[:-1],C[1:]))

                        abs = lambda x : x if x>0 else -x
                        dist = lambda X : abs(C_index[X[0]]-C_index[X[1]])

                        u,v = min(gg.edges(labels = False), key = dist)
                        u,v = C_index[u], C_index[v]

                        # Return the answer, and relabel it on-the-fly with the
                        # vertices' real name
                        return False, map(lambda x:id_label[x],C[min(u,v): max(u,v)+1])

                    else:
                        return False, None

                elif not VisitedP3.has_key((b,c,d)):
                    # search for another P_4
                    res, hole_vertices = process(b,c,d,i+1)
                    if not res:
                        return False, hole_vertices

        del InPath[c]
        return True, []


    # main algorithm
    # For all triples u,v,w of vertices such that uvw is a P_3
    for u in g:
        InPath[u] = 0   # u is the first vertex at position 0
        for vv,ww in g.edge_iterator(labels = False):
            for v,w in [(vv,ww),(ww,vv)]:
                if has_edge(dense_graph,u,v,n) and u!=w and not has_edge(dense_graph,u,w,n) and not VisitedP3.has_key((u,v,w)):
                    InPath[v] = 1   # v is the second vertex at position 1
                    res,hole = process(u, v, w, 2)
                    if not res:
                        # We relabel the graph before returning the result
                        g.relabel(id_label)
                        # Free the dense graph
                        bitset_free(dense_graph)

                        if certificate:
                            return False, g.subgraph(hole)
                        else:
                            return False
                    del InPath[v]
        del InPath[u]

    # We relabel the graph before returning the result
    g.relabel(id_label)
    # Free the dense graph
    bitset_free(dense_graph)

    if certificate:
        return True, []
    else:
        return True

def is_long_antihole_free(g, certificate = False):
    r"""
    Tests whether the given graph contains an induced subgraph that is
    isomorphic to the complement of a cycle of length at least 5.

    INPUT:

    - ``certificate`` -- boolean (default: ``False``)

      Whether to return a certificate. When ``certificate = True``, then
      the function returns

      * ``(False, Antihole)`` if ``g`` contains an induced complement
        of a cycle of length at least 5 returned as ``Antihole``.
      * ``(True, [])`` if ``g`` does not contain an induced complement of
        a cycle of length at least 5.
        For this case it is not known how to provide a certificate.

      When ``certificate = False``, the function returns just ``True`` or
      ``False`` accordingly.

    ALGORITHM:

    This algorithm tries to find a cycle in the graph of all induced
    `\overline{P_4}` of `g`, where two copies `\overline{P}` and `\overline{P'}`
    of `\overline{P_4}` are adjacent if there exists a (not necessarily induced)
    copy of `\overline{P_5}=u_1u_2u_3u_4u_5` such that
    `\overline{P}=u_1u_2u_3u_4` and `\overline{P'}=u_2u_3u_4u_5`.

    This is done through a depth-first-search. For efficiency, the auxiliary
    graph is constructed on-the-fly and never stored in memory.

    The run time of this algorithm is `O(m^2)` [NikolopoulosPalios07]_ ( where
    `m` is the number of edges of the graph ) .

    EXAMPLES:

    The Petersen Graph contains an antihole::

        sage: g = graphs.PetersenGraph()
        sage: g.is_long_antihole_free()
        False

    The complement of a cycle is an antihole::

        sage: g = graphs.CycleGraph(6).complement()
        sage: r,a = g.is_long_antihole_free(certificate=True)
        sage: r
        False
        sage: a.complement().is_isomorphic( graphs.CycleGraph(6) )
        True

    TESTS:

    Further tests::

        sage: g = Graph({0:[6,7],1:[7,8],2:[8,9],3:[9,10],4:[10,11],5:[11,6],6:[0,5,7],7:[0,1,6],8:[1,2,9],9:[2,3,8],10:[3,4,11],11:[4,5,10]}).complement()
        sage: r,a = g.is_long_antihole_free(certificate=True)
        sage: r
        False
        sage: a.complement().is_isomorphic( graphs.CycleGraph(9) )
        True
    """
    cdef int a,b,c,i,u,v,d

    # relabel the graph on 0...n-1
    cdef dict label_id = g.relabel(return_map = True)
    cdef dict id_label = {idd:label for label, idd in label_id.iteritems()}

    # A dense copy of our graph
    cdef bitset_t dense_graph
    cdef int n = g.order()
    bitset_init(dense_graph, n*n)
    bitset_set_first_n(dense_graph, 0)
    for u,v in g.edges(labels = False):
        bitset_add(dense_graph,u*n+v)
        bitset_add(dense_graph,v*n+u)

    InPath = {} #vertices of the current path with their position (InPath[v] = i)
    VisitedP3 = {} #stores triples (u,v,w) which represent visited paths of length 3


    def process(a,b,c,k):
        InPath[c] = k  # c is the (i+1)-th vertex at position i
        VisitedP3[a,c,b] = True
        VisitedP3[c,a,b] = True
        for d in g.neighbor_iterator(b):
            if has_edge(dense_graph,d,a,n) and not has_edge(dense_graph,d,c,n):
                if InPath.has_key(d):
                    if certificate:  #calculation of induced cycle in complement
                        j = InPath[d]

                        C = [v for v,vj in InPath.iteritems() if vj >= j]
                        C.sort(key = lambda x: InPath[x])
                        C_index = {u:i for i,u in enumerate(C)}

                        # At this step C[0]C[1]..... is an anticycle such that
                        # any 4 consecutive vertices induce the complement of a
                        # P_4. C may not be an induced anticycle, so we extract one
                        # from it.

                        # To do so, we look for the *shortest* nonedge C[i]C[j]
                        # between two nonconsecutive vertices of C, where the
                        # length is the difference |i-j|.
                        #
                        # C[i]...C[j] is necessarily an induced anticycle.⇧

                        gg = g.subgraph(C).complement()
                        gg.delete_edges(zip(C[:-1],C[1:]))

                        abs = lambda x : x if x>0 else -x
                        dist = lambda X : abs(C_index[X[0]]-C_index[X[1]])

                        u,v = min(gg.edges(labels = False), key = dist)
                        u,v = C_index[u], C_index[v]

                        # Return the answer, and relabel it on-the-fly with the
                        # vertices' real name
                        return False, map(lambda x:id_label[x],C[min(u,v): max(u,v)+1])

                    else:
                        return False, []

                elif not VisitedP3.has_key((b,d,c)):
                    r,antihole = process(b,c,d,k+1)
                    if not r:
                        return False, antihole

        del InPath[c]
        return True, []


    # main algorithm
    # For all triples u,v,w of vertices such that uvw is a complement of P_3
    for u in g:
        InPath[u] = 1
        for v,w in g.edge_iterator(labels = False):
            if not has_edge(dense_graph,u,v,n) and not has_edge(dense_graph,u,w,n) and not VisitedP3.has_key((v,w,u)):
                InPath[v] = 0
                r,antihole = process(v, u, w, 2)
                if not r:
                    # We relabel the graph before returning the result
                    g.relabel(id_label)
                    # Free the dense graph
                    bitset_free(dense_graph)

                    if certificate:
                        return False, g.subgraph(antihole)
                    else:
                        return False
                del InPath[v]
        del InPath[u]

    # We relabel the graph before returning the result
    g.relabel(id_label)
    # Free the dense graph
    bitset_free(dense_graph)

    if certificate:
        return True, []
    else:
        return True

def is_weakly_chordal(g, certificate = False):
    r"""
    Tests whether the given graph is weakly chordal, i.e., the graph and its
    complement have no induced cycle of length at least 5.

    INPUT:

    - ``certificate`` -- Boolean value (default: ``False``) whether to
      return a certificate. If ``certificate = False``, return ``True`` or
      ``False`` according to the graph. If ``certificate = True``, return

      * ``(False, forbidden_subgraph)`` when the graph contains a
        forbidden subgraph H, this graph is returned.
      * ``(True, [])`` when the graph is weakly chordal.
          For this case, it is not known how to provide a certificate.

    ALGORITHM:

    This algorithm checks whether the graph ``g`` or its complement
    contain an induced cycle of length at least 5.

    Using is_long_hole_free() and is_long_antihole_free() yields a run time
    of `O(m^2)` (where `m` is the number of edges of the graph).

    EXAMPLES:

    The Petersen Graph is not weakly chordal and contains a hole::

        sage: g = graphs.PetersenGraph()
        sage: r,s = g.is_weakly_chordal(certificate = True)
        sage: r
        False
        sage: l = len(s.vertices())
        sage: s.is_isomorphic( graphs.CycleGraph(l) )
        True
    """

    if certificate:
        r,forbid_subgr = g.is_long_hole_free(certificate=True)
        if not r:
            return False, forbid_subgr

        return g.is_long_antihole_free(certificate=True)
    else:
        return g.is_long_hole_free() and g.is_long_antihole_free()

