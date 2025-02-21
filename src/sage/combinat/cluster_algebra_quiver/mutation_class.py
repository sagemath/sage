# sage.doctest: needs sage.graphs
r"""
mutation_class

This file contains helper functions for compute the mutation class of a cluster algebra or quiver.

For the compendium on the cluster algebra and quiver package see [MS2011]_

AUTHORS:

- Gregg Musiker
- Christian Stump
"""
# ****************************************************************************
#       Copyright (C) 2011 Gregg Musiker <musiker@math.mit.edu>
#                          Christian Stump <christian.stump@univie.ac.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import time
from sage.groups.perm_gps.partn_ref.refinement_graphs import search_tree, get_orbits
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import infinity
from sage.graphs.digraph import DiGraph
from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import _edge_list_to_matrix


def _principal_part(mat):
    """
    Return the principal part of a matrix.

    INPUT:

    - ``mat`` -- a matrix with at least as many rows as columns

    OUTPUT: the top square part of the matrix ``mat``

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _principal_part
        sage: M = Matrix([[1,2],[3,4],[5,6]]); M                                        # needs sage.modules
        [1 2]
        [3 4]
        [5 6]
        sage: _principal_part(M)                                                        # needs sage.modules
        [1 2]
        [3 4]
    """
    n, m = mat.ncols(), mat.nrows() - mat.ncols()
    if m < 0:
        raise ValueError('the input matrix has more columns than rows')
    elif m == 0:
        return mat
    else:
        return mat.submatrix(0, 0, n, n)


def _digraph_mutate(dg, k, frozen=None):
    """
    Return a digraph obtained from ``dg`` by mutating at vertex ``k``.

    Vertices can be labelled by anything, and frozen vertices must
    be explicitly given.

    INPUT:

    - ``dg`` -- a digraph with integral edge labels with ``n+m`` vertices
    - ``k`` -- the vertex at which ``dg`` is mutated
    - ``frozen`` -- the list of frozen vertices (default: empty list)

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _digraph_mutate
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['A',4]).digraph()                                     # needs sage.modules
        sage: dg.edges(sort=True)                                                       # needs sage.modules
        [(0, 1, (1, -1)), (2, 1, (1, -1)), (2, 3, (1, -1))]
        sage: _digraph_mutate(dg,2).edges(sort=True)                                    # needs sage.modules
        [(0, 1, (1, -1)), (1, 2, (1, -1)), (3, 2, (1, -1))]

    TESTS::

       sage: dg = DiGraph([('a','b',(1,-1)),('c','a',(1,-1))])
       sage: _digraph_mutate(dg,'a').edges(sort=True)
       [('a', 'c', (1, -1)), ('b', 'a', (1, -1)), ('c', 'b', (1, -1))]
       sage: _digraph_mutate(dg,'a',frozen=['b','c']).edges(sort=True)
       [('a', 'c', (1, -1)), ('b', 'a', (1, -1))]

       sage: dg = DiGraph([('a','b',(2,-2)),('c','a',(2,-2)),('b','c',(2,-2))])
       sage: _digraph_mutate(dg,'a').edges(sort=True)
       [('a', 'c', (2, -2)), ('b', 'a', (2, -2)), ('c', 'b', (2, -2))]
    """
    # assert sorted(list(dg)) == list(range(n + m))
    # this is not assumed anymore
    if frozen is None:
        frozen = []

    edge_it = dg.incoming_edge_iterator(dg, True)
    edges = {(v1, v2): label for v1, v2, label in edge_it}
    edge_it = dg.incoming_edge_iterator([k], True)
    in_edges = list(edge_it)
    edge_it = dg.outgoing_edge_iterator([k], True)
    out_edges = list(edge_it)

    in_edges_new = [(v2, v1, (-label[1], -label[0]))
                    for (v1, v2, label) in in_edges]
    out_edges_new = [(v2, v1, (-label[1], -label[0]))
                     for (v1, v2, label) in out_edges]
    diag_edges_new = []
    diag_edges_del = []

    for (v1, v2, label1) in in_edges:
        l11, l12 = label1
        for (w1, w2, label2) in out_edges:
            if v1 in frozen and w2 in frozen:
                continue
            l21, l22 = label2
            if (v1, w2) in edges:
                diag_edges_del.append((v1, w2))
                a, b = edges[(v1, w2)]
                a, b = a + l11 * l21, b - l12 * l22
                diag_edges_new.append((v1, w2, (a, b)))
            elif (w2, v1) in edges:
                diag_edges_del.append((w2, v1))
                a, b = edges[(w2, v1)]
                a, b = b + l11 * l21, a - l12 * l22
                if a < 0:
                    diag_edges_new.append((w2, v1, (b, a)))
                elif a > 0:
                    diag_edges_new.append((v1, w2, (a, b)))
            else:
                a, b = l11 * l21, -l12 * l22
                diag_edges_new.append((v1, w2, (a, b)))

    del_edges = [tuple(ed[:2]) for ed in in_edges + out_edges]
    del_edges += diag_edges_del
    new_edges = in_edges_new + out_edges_new
    new_edges += diag_edges_new
    new_edges += [(v1, v2, edges[(v1, v2)]) for (v1, v2) in edges
                  if (v1, v2) not in del_edges]

    dg_new = DiGraph()
    dg_new.add_vertices(list(dg))
    for v1, v2, label in new_edges:
        dg_new.add_edge(v1, v2, label)

    return dg_new


def _matrix_to_digraph(M):
    """
    Return the digraph obtained from the matrix ``M``.

    In order to generate a quiver, we assume that ``M`` is skew-symmetrizable.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _matrix_to_digraph
        sage: _matrix_to_digraph(matrix(3,[0,1,0,-1,0,-1,0,1,0]))                       # needs sage.modules
        Digraph on 3 vertices
    """
    n = M.ncols()

    dg = DiGraph(sparse=True)
    for i, j in M.nonzero_positions():
        if i >= n:
            a, b = M[i, j], -M[i, j]
        else:
            a, b = M[i, j], M[j, i]
        if a > 0:
            dg._backend.add_edge(i, j, (a, b), True)
        elif i >= n:
            dg._backend.add_edge(j, i, (-a, -b), True)
    for i in range(M.nrows()):
        if i not in dg:
            dg.add_vertex(i)
    return dg


def _dg_canonical_form(dg, frozen=None):
    """
    Turn the digraph ``dg`` into its canonical form, and return the
    corresponding isomorphism and the vertex orbits of the automorphism group.

    The labels of ``dg`` are assumed to be integers between `0` and `n + m - 1`.

    The canonical form has the following additional property: the frozen
    vertices are the final labels.

    .. WARNING:: The input ``dg`` is modified.

    INPUT:

    - ``dg`` -- a directed graph having edge labels (a, b) with a > 0

    - ``frozen`` -- list (default: ``[]``) of frozen vertices

    OUTPUT: dictionary {original label: canonical label}

    - list of orbits of mutable vertices (using canonical labels)

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _dg_canonical_form
        sage: dg = ClusterQuiver(['B',4]).digraph(); dg.edges(sort=True)                # needs sage.modules
        [(0, 1, (1, -1)), (2, 1, (1, -1)), (2, 3, (1, -2))]
        sage: _dg_canonical_form(dg); dg.edges(sort=True)                               # needs sage.modules
        ({0: 0, 1: 3, 2: 1, 3: 2}, [[0], [3], [1], [2]])
        [(0, 3, (1, -1)), (1, 2, (1, -2)), (1, 3, (1, -1))]

    TESTS::

        sage: dg2 = ClusterQuiver(DiGraph({0:[1,2]})).digraph()                         # needs sage.modules
        sage: _dg_canonical_form(dg2); dg2.edges(sort=True)                             # needs sage.modules
        ({0: 0, 1: 1, 2: 2}, [[0], [1, 2]])
        [(0, 1, (1, -1)), (0, 2, (1, -1))]

        sage: dg2 = ClusterQuiver(DiGraph({0:[1,2]})).digraph()                         # needs sage.modules
        sage: _dg_canonical_form(dg2, frozen=[0]); dg2.edges(sort=True)                 # needs sage.modules
        ({0: 2, 1: 0, 2: 1}, [[2], [0, 1]])
        [(2, 0, (1, -1)), (2, 1, (1, -1))]

        sage: dg3 = ClusterQuiver(DiGraph({0:[1,2],1:[3]})).digraph()                   # needs sage.modules
        sage: _dg_canonical_form(dg3, frozen=[0,3]); dg3.edges(sort=True)               # needs sage.modules
        ({0: 2, 1: 1, 2: 0, 3: 3}, [[2], [1], [0], [3]])
        [(1, 3, (1, -1)), (2, 0, (1, -1)), (2, 1, (1, -1))]

        sage: dg3 = ClusterQuiver(DiGraph({2:[1,3],1:[0],3:[4]})).digraph()             # needs sage.modules
        sage: _dg_canonical_form(dg3, frozen=[4,0]); dg3.edges(sort=True)               # needs sage.modules
        ({0: 4, 1: 1, 2: 0, 3: 2, 4: 3}, [[4, 3], [1, 2], [0]])
        [(0, 1, (1, -1)), (0, 2, (1, -1)), (1, 4, (1, -1)), (2, 3, (1, -1))]
    """
    vertices = list(dg)
    n_plus_m = dg.order()
    if frozen is not None:
        partition = [[v for v in vertices if v not in frozen], frozen]
    else:
        frozen = []
        partition = [ vertices ]

    # here dg is modified by inserting new vertices at middle of edges
    partition_add, edges = _graph_without_edge_labels(dg, vertices)

    partition += partition_add
    automorphism_group, _, iso = search_tree(dg, partition=partition, lab=True,
                                             dig=True, certificate=True)
    orbits = get_orbits(automorphism_group, n_plus_m)
    orbits = [[iso[i] for i in orbit] for orbit in orbits]

    removed = []
    for v in iso:
        if v not in vertices:
            removed.append(v)
            v1, _, _ = next(dg._backend.iterator_in_edges([v], True))
            _, w2, _ = next(dg._backend.iterator_out_edges([v], True))
            dg._backend.del_vertex(v)
            add_index = True
            index = 0
            while add_index:
                l = partition_add[index]
                if v in l:
                    dg._backend.add_edge(v1, w2, edges[index], True)
                    add_index = False
                index += 1
    for v in removed:
        del iso[v]

    dg._backend.relabel(iso, True)
    return iso, orbits


def _mutation_class_iter( dg, n, m, depth=infinity, return_dig6=False, show_depth=False, up_to_equivalence=True, sink_source=False ):
    """
    Return an iterator for mutation class of dg with respect to several parameters.

    .. NOTE:: assuming that the frozen vertices start at vertex n.

    INPUT:

    - ``dg`` -- a digraph with `n+m` vertices
    - ``depth`` -- positive integer or infinity specifying (roughly) how many
      steps away from the initial seed to mutate
    - ``return_dig6`` -- indicates whether to convert digraph data to dig6
      string data
    - ``show_depth`` -- if ``True``, indicates that a running count of the
      depth is to be displayed
    - ``up_to_equivalence`` -- if ``True``, only one digraph for each
      graph-isomorphism class is recorded
    - ``sink_source`` -- if ``True``, only mutations at sinks or sources are
      applied

    EXAMPLES::

        sage: # needs sage.modules
        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _mutation_class_iter
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['A',[1,2],1]).digraph()
        sage: itt = _mutation_class_iter(dg, 3,0)
        sage: next(itt)[0].edges(sort=True)
        [(0, 1, (1, -1)), (0, 2, (1, -1)), (1, 2, (1, -1))]
        sage: next(itt)[0].edges(sort=True)
        [(0, 2, (1, -1)), (1, 0, (2, -2)), (2, 1, (1, -1))]
    """
    # assuming that the frozen vertices are at the end (from n to n + m - 1)
    # assert sorted(dg) == list(range(n + m))  # commented out for speed
    mlist = list(range(n, n + m))

    timer = time.time()
    depth_counter = 0
    if up_to_equivalence:
        iso, orbits = _dg_canonical_form( dg, mlist )
        iso_inv = { iso[a]: a for a in iso }

    dig6 = _digraph_to_dig6( dg, hashable=True )
    dig6s = {}
    if up_to_equivalence:
        orbits = [ orbit[0] for orbit in orbits ]
        dig6s[dig6] = [ orbits, [], iso_inv ]
    else:
        dig6s[dig6] = [list(range(n)), [] ]
    if return_dig6:
        yield (dig6, [])
    else:
        yield (dg, [])

    gets_bigger = True
    if show_depth:
        timer2 = time.time()
        dc = str(depth_counter)
        dc += ' ' * (5 - len(dc))
        nr = str(len(dig6s))
        nr += ' ' * (10 - len(nr))
        print(f"Depth: {dc} found: {nr} Time: {timer2 - timer:.2f} s")

    while gets_bigger and depth_counter < depth:
        gets_bigger = False
        for key in list(dig6s):
            mutation_indices = [ i for i in dig6s[key][0] if i < n ]
            if mutation_indices:
                dg = _dig6_to_digraph( key )
            while mutation_indices:
                i = mutation_indices.pop()
                if not sink_source or _dg_is_sink_source( dg, i ):
                    dg_new = _digraph_mutate(dg, i, frozen=mlist)
                    if up_to_equivalence:
                        iso, orbits = _dg_canonical_form( dg_new, mlist )
                        i_new = iso[i]
                        iso_inv = { iso[a]: a for a in iso }
                    else:
                        i_new = i
                    dig6_new = _digraph_to_dig6( dg_new, hashable=True )
                    if dig6_new in dig6s:
                        if i_new in dig6s[dig6_new][0]:
                            dig6s[dig6_new][0].remove(i_new)
                    else:
                        gets_bigger = True
                        if up_to_equivalence:
                            orbits = [ orbit[0] for orbit in orbits if i_new not in orbit ]
                            iso_history = { a: dig6s[key][2][iso_inv[a]] for a in iso }
                            i_history = iso_history[i_new]
                            history = dig6s[key][1] + [i_history]
                            dig6s[dig6_new] = [orbits,history,iso_history]
                        else:
                            orbits = list(range(n))
                            del orbits[i_new]
                            history = dig6s[key][1] + [i_new]
                            dig6s[dig6_new] = [orbits,history]
                        if return_dig6:
                            yield (dig6_new,history)
                        else:
                            yield (dg_new,history)
        depth_counter += 1
        if show_depth and gets_bigger:
            timer2 = time.time()
            dc = str(depth_counter)
            dc += ' ' * (5 - len(dc))
            nr = str(len(dig6s))
            nr += ' ' * (10 - len(nr))
            print(f"Depth: {dc} found: {nr} Time: {timer2 - timer:.2f} s")


def _digraph_to_dig6(dg, hashable=False):
    """
    Return the dig6 and edge data of the digraph dg.

    INPUT:

    - ``dg`` -- a digraph
    - ``hashable`` -- boolean (default: ``False``); if ``True``, the edge
      labels are turned into a dict

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _digraph_to_dig6
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['A',4]).digraph()                                     # needs sage.modules
        sage: _digraph_to_dig6(dg)                                                      # needs sage.modules
        ('COD?', {})
    """
    dig6 = dg.dig6_string()
    D = {}
    for E in dg._backend.iterator_in_edges(dg,True):
        if E[2] != (1,-1):
            D[ (E[0],E[1]) ] = E[2]
    if hashable:
        D = tuple( sorted( D.items() ) )
    return (dig6,D)


def _dig6_to_digraph( dig6 ):
    """
    Return the digraph obtained from the dig6 and edge data.

    INPUT:

    - ``dig6`` -- a pair ``(dig6, edges)`` where ``dig6`` is a string encoding a digraph and ``edges`` is a dict or tuple encoding edges

    EXAMPLES::

        sage: # needs sage.modules
        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _digraph_to_dig6
        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _dig6_to_digraph
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['A',4]).digraph()
        sage: data = _digraph_to_dig6(dg)
        sage: _dig6_to_digraph(data)
        Digraph on 4 vertices
        sage: _dig6_to_digraph(data).edges(sort=True)
        [(0, 1, (1, -1)), (2, 1, (1, -1)), (2, 3, (1, -1))]
    """
    dig6, edges = dig6
    dg = DiGraph( dig6 )
    if not isinstance(edges, dict):
        edges = dict( edges )
    for edge in dg._backend.iterator_in_edges(dg,False):
        if edge in edges:
            dg.set_edge_label( edge[0],edge[1],edges[edge] )
        else:
            dg.set_edge_label( edge[0],edge[1], (1,-1) )
    return dg


def _dig6_to_matrix( dig6 ):
    """
    Return the matrix obtained from the dig6 and edge data.

    INPUT:

    - ``dig6`` -- a pair ``(dig6, edges)`` where ``dig6`` is a string
      encoding a digraph and ``edges`` is a dict or tuple encoding edges

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _digraph_to_dig6, _dig6_to_matrix
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['A',4]).digraph()                                     # needs sage.modules
        sage: data = _digraph_to_dig6(dg)                                               # needs sage.modules
        sage: _dig6_to_matrix(data)                                                     # needs sage.modules
        [ 0  1  0  0]
        [-1  0 -1  0]
        [ 0  1  0  1]
        [ 0  0 -1  0]
    """
    dg = _dig6_to_digraph(dig6)
    return _edge_list_to_matrix(dg.edges(sort=True), list(range(dg.order())), [])


def _dg_is_sink_source( dg, v ):
    """
    Return ``True`` iff the digraph dg has a sink or a source at vertex `v`.

    INPUT:

    - ``dg`` -- a digraph
    - ``v`` -- a vertex of dg

    EXAMPLES::

        sage: # needs sage.modules
        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _dg_is_sink_source
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['A',[1,2],1]).digraph()
        sage: _dg_is_sink_source(dg, 0)
        True
        sage: _dg_is_sink_source(dg, 1)
        True
        sage: _dg_is_sink_source(dg, 2)
        False
    """
    in_edges = list(dg._backend.iterator_in_edges([v],True))
    out_edges = list(dg._backend.iterator_out_edges([v],True))
    return not ( in_edges and out_edges )


def _graph_without_edge_labels(dg, vertices):
    """
    Expand the graph, by introducing new vertices in the middle of edges.
    Return the corresponding partition of the new vertices.

    There is one new vertex for each edge with label not equal to ``(1,-1)``.
    These vertices are numbered starting from the smallest available integer.

    Each edge having a non-trivial label is replaced by two consecutive edges,
    passing through one new vertex.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _graph_without_edge_labels
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['B',4]).digraph(); dg.edges(sort=True)                # needs sage.modules
        [(0, 1, (1, -1)), (2, 1, (1, -1)), (2, 3, (1, -2))]
        sage: _graph_without_edge_labels(dg, range(4)); dg.edges(sort=True)             # needs sage.modules
        ([[4]], ((1, -2),))
        [(0, 1, (1, -1)), (2, 1, (1, -1)), (2, 4, (1, -1)), (4, 3, (1, -1))]
    """
    vertices = list(vertices)
    edges = dg.edge_iterator(labels=True)
    edge_labels = tuple(sorted({label for _, _, label in edges
                                if label != (1, -1)}))
    edge_partition = [[] for _ in edge_labels]
    i = 0
    while i in vertices:
        i += 1

    # We have to pass the old vertices to the iterator.
    # Otherwise the ``edge_iterator`` will use the ``vertices``
    # of ``dg``, which are not well-defined, if we add vertices along the way.
    for u, v, label in dg.edge_iterator(vertices=vertices, labels=True):
        if label != (1, -1):  # Ignore edges we just added.
            index = edge_labels.index(label)
            edge_partition[index].append(i)
            dg._backend.add_edge(u, i, (1, -1), True)
            dg._backend.add_edge(i, v, (1, -1), True)
            dg._backend.del_edge(u, v, label, True)
            i += 1
            while i in vertices:
                i += 1
    return edge_partition, edge_labels


def _has_two_cycles( dg ):
    """
    Return ``True`` if the input digraph has a 2-cycle and ``False`` otherwise.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _has_two_cycles
        sage: _has_two_cycles(DiGraph([[0,1],[1,0]]))
        True
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: _has_two_cycles(ClusterQuiver(['A',3]).digraph())                         # needs sage.modules
        False
    """
    edge_set = dg.edges(sort=True, labels=False)
    for (v,w) in edge_set:
        if (w,v) in edge_set:
            return True
    return False


def _is_valid_digraph_edge_set( edges, frozen=0 ):
    """
    Return ``True`` if the input data is the edge set of a digraph for a quiver
    (no loops, no 2-cycles, edge-labels of the specified format), and return
    ``False`` otherwise.

    INPUT:

    - ``frozen`` -- integer (default: 0); the number of frozen vertices

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _is_valid_digraph_edge_set
        sage: _is_valid_digraph_edge_set( [[0,1,'a'],[2,3,(1,-1)]] )
        The given digraph has edge labels which are not integral or integral 2-tuples.
        False
        sage: _is_valid_digraph_edge_set( [[0,1,None],[2,3,(1,-1)]] )
        True
        sage: _is_valid_digraph_edge_set( [[0,1,'a'],[2,3,(1,-1)],[3,2,(1,-1)]] )
        The given digraph or edge list contains oriented 2-cycles.
        False
    """
    try:
        dg = DiGraph()
        dg.allow_multiple_edges(True)
        dg.add_edges( edges )

        # checks if the digraph contains loops
        if dg.has_loops():
            print("The given digraph or edge list contains loops")
            return False

        # checks if the digraph contains oriented 2-cycles
        if _has_two_cycles( dg ):
            print("The given digraph or edge list contains oriented 2-cycles.")
            return False

        # checks if all edge labels are 'None', positive integers or tuples of positive integers
        if not all( i is None or ( i in ZZ and i > 0 ) or ( isinstance(i, tuple) and len(i) == 2 and i[0] in ZZ and i[1] in ZZ ) for i in dg.edge_labels() ):
            print("The given digraph has edge labels which are not integral or integral 2-tuples.")
            return False

        # checks if all edge labels for multiple edges are 'None' or positive integers
        if dg.has_multiple_edges():
            for e in set( dg.multiple_edges(labels=False) ):
                if not all( i is None or ( i in ZZ and i > 0 ) for i in dg.edge_label( e[0], e[1] ) ):
                    print("The given digraph or edge list contains multiple edges with non-integral labels.")
                    return False

        n = dg.order() - frozen
        if n < 0:
            print("The number of frozen variables is larger than the number of vertices.")
            return False

        if any(e[0] >= n for e in dg.edges(sort=True, labels=False)):
            print("The given digraph or edge list contains edges within the frozen vertices.")
            return False

        return True
    except Exception:
        print("Could not even build a digraph from the input data.")
        return False
