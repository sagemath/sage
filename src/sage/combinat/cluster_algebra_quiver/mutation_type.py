# sage_setup: distribution = sagemath-graphs
# sage.doctest: needs sage.graphs
r"""
Helper functions for mutation types of quivers

This file contains helper functions for detecting the mutation type of
a cluster algebra or quiver.

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

from pathlib import Path
import pickle

from copy import copy

from sage.misc.cachefunc import cached_function
from sage.misc.flatten import flatten
from sage.graphs.digraph import DiGraph
from sage.combinat.combination import Combinations
from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import QuiverMutationType


def is_mutation_finite(M, nr_of_checks=None):
    r"""
    Use a non-deterministic method by random mutations in various
    directions. Can result in a wrong answer.

    .. WARNING::

        This method modifies the input matrix ``M``!

    INPUT:

    - ``nr_of_checks`` -- number of mutations applied (default: ``None``);
      standard is 500*(number of vertices of self)

    ALGORITHM:

    A quiver is mutation infinite if and only if every edge label (a,-b) satisfy a*b > 4.
    Thus, we apply random mutations in random directions

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_type import is_mutation_finite

        sage: Q = ClusterQuiver(['A',10])                                               # needs sage.modules
        sage: M = Q.b_matrix()                                                          # needs sage.modules
        sage: is_mutation_finite(M)                                                     # needs sage.modules
        (True, None)

        sage: # needs sage.modules
        sage: Q = ClusterQuiver([(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(2,9)])
        sage: M = Q.b_matrix()
        sage: is_mutation_finite(M)  # random
        (False, [9, 6, 9, 8, 9, 4, 0, 4, 5, 2, 1, 0, 1, 0, 7, 1, 9, 2, 5, 7, 8, 6, 3, 0, 2, 5, 4, 2, 6, 9, 2, 7, 3, 5, 3, 7, 9, 5, 9, 0, 2, 7, 9, 2, 4, 2, 1, 6, 9, 4, 3, 5, 0, 8, 2, 9, 5, 3, 7, 0, 1, 8, 3, 7, 2, 7, 3, 4, 8, 0, 4, 9, 5, 2, 8, 4, 8, 1, 7, 8, 9, 1, 5, 0, 8, 7, 4, 8, 9, 8, 0, 7, 4, 7, 1, 2, 8, 6, 1, 3, 9, 3, 9, 1, 3, 2, 4, 9, 5, 1, 2, 9, 4, 8, 5, 3, 4, 6, 8, 9, 2, 5, 9, 4, 6, 2, 1, 4, 9, 6, 0, 9, 8, 0, 4, 7, 9, 2, 1, 6])

    Check that :issue:`19495` is fixed::

        sage: dg = DiGraph(); dg.add_vertex(0); S = ClusterSeed(dg); S                  # needs sage.modules
        A seed for a cluster algebra of rank 1
        sage: S.is_mutation_finite()                                                    # needs sage.modules
        True
    """
    import random
    n = M.ncols()
    if n <= 2:
        return True, None
    if nr_of_checks is None:
        nr_of_checks = 1000 * n
    k = random.randint(0, n - 1)
    path = []
    for i in range(nr_of_checks):
        # avoid mutating back in the same direction
        next_k = random.randint(0, n - 2)
        if next_k >= k:
            next_k += 1
        k = next_k
        M.mutate(k)
        path.append(k)
        for i, j in M.nonzero_positions():
            if i < j and M[i, j] * M[j, i] < -4:
                return False, path
    return True, None


def _triangles(dg):
    """
    Return a list of all oriented triangles in the digraph ``dg``.

    EXAMPLES::

        sage: # needs sage.modules
        sage: from sage.combinat.cluster_algebra_quiver.mutation_type import _triangles
        sage: Q = ClusterQuiver(['A',3])
        sage: _triangles(Q.digraph())
        []
        sage: Q.mutate([0,1])
        sage: _triangles(Q.digraph())
        [([(2, 0), (0, 1), (1, 2)], True)]
        sage: Q2 = ClusterQuiver(['A',[1,2],1])
        sage: _triangles(Q2.digraph())
        [([(1, 2), (1, 0), (2, 0)], False)]
        sage: Q2.mutate(2)
        sage: _triangles(Q2.digraph())
        [([(1, 0), (0, 2), (2, 1)], True)]
    """
    E = dg.edges(sort=True, labels=False)
    V = list(dg)
    trians = []
    flat_trians = []
    for e in E:
        v1, v2 = e
        for v in V:
            if v not in e:
                if (v, v1) in E:
                    if (v2, v) in E:
                        flat_trian = sorted([v,v1,v2])
                        if flat_trian not in flat_trians:
                            flat_trians.append( flat_trian )
                            trians.append( ( [(v,v1),(v1,v2),(v2,v)], True ) )
                    elif (v, v2) in E:
                        flat_trian = sorted([v,v1,v2])
                        if flat_trian not in flat_trians:
                            flat_trians.append( flat_trian )
                            trians.append( ( [(v,v1),(v1,v2),(v,v2)], False ) )
                if (v1, v) in E:
                    if (v2, v) in E:
                        flat_trian = sorted([v,v1,v2])
                        if flat_trian not in flat_trians:
                            flat_trians.append( flat_trian )
                            trians.append( ( [(v1,v),(v1,v2),(v2,v)], False ) )
                    elif (v, v2) in E:
                        flat_trian = sorted([v,v1,v2])
                        if flat_trian not in flat_trians:
                            flat_trians.append( flat_trian )
                            trians.append( ( [(v1,v),(v1,v2),(v,v2)], False ) )
    return trians


def _all_induced_cycles_iter(dg):
    """
    Return an iterator for all induced oriented cycles of length
    greater than or equal to 4 in the digraph ``dg``.

    EXAMPLES::

        sage: # needs sage.modules
        sage: from sage.combinat.cluster_algebra_quiver.mutation_type import _all_induced_cycles_iter
        sage: Q = ClusterQuiver(['A',[6,0],1]); Q
        Quiver on 6 vertices of type ['D', 6]
        sage: next(_all_induced_cycles_iter(Q.digraph()))
        ([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)], True)
        sage: Q.mutate(0)
        sage: next(_all_induced_cycles_iter(Q.digraph()))
        ([(1, 2), (2, 3), (3, 4), (4, 5), (5, 1)], True)
        sage: Q2 = ClusterQuiver(['A',[2,3],1])
        sage: next(_all_induced_cycles_iter(Q2.digraph()))
        ([(1, 0), (1, 2), (3, 2), (3, 4), (4, 0)], False)
    """
    dg_new = DiGraph(dg)
    E = dg_new.edges(sort=True)
    for v1, v2, label in E:
        dg_new.add_edge((v2, v1, label))
    induced_sets = []
    cycle_iter = dg_new.all_cycles_iterator(simple=True)
    for cycle in cycle_iter:
        if len(cycle) > 3:
            cycle_set = set(cycle)
            if not any(cycle_set.issuperset(induced_set) for induced_set in induced_sets):
                induced_sets.append(cycle_set)
                if len(cycle) > 4:
                    sg = dg.subgraph(cycle)
                    is_oriented = True
                    V = list(sg)
                    while is_oriented and V:
                        v = V.pop()
                        if not sg.in_degree(v) == 1:
                            is_oriented = False
                    yield (sg.edges(sort=True, labels=False), is_oriented)

# a debug function


def _false_return(s=False):
    """
    Return 'unknown'.

    Written for potential debugging purposes.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_type import _false_return
        sage: _false_return()
        'unknown'
    """
#    Uncomment these three lines for debugging purposes.
#    if s:
#        print('DEBUG: error %s' % s)
    return 'unknown'


def _reset_dg(dg, vertices, dict_in_out, del_vertices):
    """
    Delete the specified vertices (del_vertices) from the DiGraph dg,
    and the lists vertices and dict_in_out.

    Note that vertices and dict_in_out are the vertices of dg and a
    dictionary of in- and out-degrees that depend on the digraph
    ``dg`` but they are passed through as arguments so the function
    can change their values.

    EXAMPLES::

        sage: # needs sage.modules
        sage: from sage.combinat.cluster_algebra_quiver.mutation_type import _reset_dg
        sage: dg = ClusterQuiver(['A',[2,2],1]).digraph(); dg
        Digraph on 4 vertices
        sage: vertices = list(dg)
        sage: dict_in_out = {}
        sage: for v in vertices: dict_in_out[v] = (dg.in_degree(v), dg.out_degree(v), dg.degree(v))
        sage: _reset_dg(dg,vertices, dict_in_out, [1])
        sage: dg
        Digraph on 3 vertices
        sage: vertices
        [0, 2, 3]
        sage: dict_in_out
        {0: (1, 0, 1), 2: (1, 0, 1), 3: (0, 2, 2)}
    """
    del_vertices = list(set(del_vertices))
    for v in del_vertices:
        if v in dg:
            dg.delete_vertex(v)
        else:
            print(v)
            print(dg.edges(sort=True))
        vertices.remove(v)
        del dict_in_out[v]
    for v in vertices:
        dict_in_out[v] = (dg.in_degree(v), dg.out_degree(v), dg.degree(v))


def _check_special_BC_cases(dg, n, check_letter_list, check_twist_list,
                            hope_letter_list, conn_vert_list=False):
    """
    Test if dg (on at most `n` vertices) is a quiver of type `A` or
    `D` (as given in hope_letter_list) with conn_vert_list (if
    given) as connecting vertices.

    Since this is supposed to be run on a ``dg`` coming from a larger
    quiver where vertices have already been removed (outside of the
    connecting vertices), this program therefore recognizes the type
    of the larger quiver as an `n`-vertex quiver of letter on
    ``check_letter_list`` and twist on ``check_twist_list``.  This
    method is utilized in _connected_mutation_type to test for types
    BC, BB, CC, BD, or CD.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_type import _check_special_BC_cases
        sage: dg = DiGraph(1)
        sage: _check_special_BC_cases(dg, 3, ['BC'], [1], ['A'], [[0]])
        ['BC', 2, 1]
        sage: dg = DiGraph(2); dg.add_edge([0,1])
        sage: _check_special_BC_cases(dg, 4, ['BC'], [1], ['A'], [[0]])
        ['BC', 3, 1]
        sage: dg = DiGraph(2); dg.add_edge([0,1])
        sage: _check_special_BC_cases(dg, 4, ['BB'], [1], ['A'], [[0,1]])
        ['BB', 3, 1]
        sage: _check_special_BC_cases(dg, 4, ['C', 'CD'], [None, None], ['A', 'D'], [ [], [0] ])
        ['C', 4]
        sage: dg.add_edges([[1, 2], [1, 3]])
        sage: _check_special_BC_cases(dg, 4, ['C', 'CD'], [None, None], ['A', 'D'], [ [], [0] ])
        ['CD', 3, 1]
    """
    # if dg is not connected, mutation type is not recognized.
    if not dg.is_connected():
        return 'unknown'
    # divides into cases depending on whether or not a list 'conn_vert_list' of connecting vertices is given.
    if conn_vert_list:
        mut_type = _connected_mutation_type_AAtildeD( dg, ret_conn_vert=True )
        # when 'conn_vert_list' is given, the output of _connected_mutation_type_AAtildeD is
        # either 'unknown' or a pair (mut_type, conn_verts).  Then, it is tested if the vertices can be glued together as desired.
        if not mut_type == 'unknown':
            mut_type, conn_verts = mut_type
    else:
        # when conn_vert_list == False, the output of _connected_mutation_type _AAtildeD is simply 'unknown' or the mutation type.
        # no 'connecting vertices' need to be computed.
        mut_type = _connected_mutation_type_AAtildeD( dg, ret_conn_vert=False )
        conn_verts = []
    # when the mutation type is recognized, program now tries more specifically to figure out 'letter' and 'twist'
    if not mut_type == 'unknown':
        for i in range( len( check_letter_list ) ):
            check_letter = check_letter_list[i]
            check_twist = check_twist_list[i]
            hope_letter = hope_letter_list[i]
            if conn_vert_list:
                conn_vert = set(conn_vert_list[i])
            else:
                conn_vert = set()
            # Now, tries to connect up the quiver components (keeping in mind ['D',3] - ['A',3] equivalence)
            if hope_letter == 'D' and mut_type._letter == 'A' and mut_type._rank == 3 and not mut_type._twist:
                hope_letter = 'A'
                if conn_vert_list:
                    conn_verts = list(set(dg).difference(conn_verts))
            if mut_type._letter == hope_letter and not mut_type._twist and conn_vert.issubset(conn_verts):
                if len(check_letter) > 1:
                    check_twist = 1
                if check_twist:
                    n -= 1
                return QuiverMutationType([check_letter, n, check_twist])
    return 'unknown'


def _connected_mutation_type(dg):
    """
    Assuming that ``dg`` is a connected digraph, checks the mutation
    type of ``dg`` as a valued quiver.

    EXAMPLES::

        sage: # needs sage.modules
        sage: from sage.combinat.cluster_algebra_quiver.mutation_type import _connected_mutation_type
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['A',3]).digraph(); _connected_mutation_type( dg )
        ['A', 3]
        sage: dg = ClusterQuiver(['D',7]).digraph(); _connected_mutation_type( dg )
        ['D', 7]
        sage: dg = ClusterQuiver(['BC',4,1]).digraph(); _connected_mutation_type( dg )
        ['BC', 4, 1]
    """
    dg = DiGraph( dg )
    # defining some shorthands
    n = dg.order()
    edges = dg.edges(sort=True)
    vertices = list(dg)
    # initializing lists of the edges with labels (2,-1) or (1,-2); (4,-1) or (1,-4); or (2,-2), respectively
    exc_labels = []
    exc_labels41 = []
    double_edges = []
    # letter = None

    # replacing higher labels by multiple edges.  Multiple edges and acyclic is a sign that quiver is infinite mutation type with the exception of A_tilde where there might be one multiple edge with multiplicity 2.  Multiple edges is at least a sign that the quiver is of 'undetermined finite mutation type'.
    dg.allow_multiple_edges( True )
    for edge in edges:
        label = edge[2]
        if label not in [(1,-1),(2,-2),(1,-2),(2,-1),(4,-1),(1,-4)]:
            # _false_return(i) is a simple function that simply returns 'unknown'.  For debugging purposes, it
            # can also output 'DEBUG: error i' if desired.
            # this command is used many times in this code, something times without the argument i.
            return _false_return(2)
        elif label == (2,-2):
            dg.set_edge_label( edge[0], edge[1], 1 )
            dg.add_edge( edge[0], edge[1], 1 )
            double_edges.append( edge )
            if len( double_edges ) > 1:
                return _false_return()
        elif label == (1,-1):
            dg.set_edge_label( edge[0], edge[1], 1 )
        elif label in [(2,-1),(1,-2)]:
            exc_labels.append( edge )
        elif label in [(1,-4),(4,-1)]:
            exc_labels41.append( edge )

    # creating a dictionary of in-, out- and total degrees
    dict_in_out = {}
    for v in vertices:
        dict_in_out[v] = (dg.in_degree(v), dg.out_degree(v), dg.degree(v))

    if len( exc_labels ) + len( exc_labels41 ) + len( double_edges ) > 4:
        return _false_return()

    # test for the labels (4,-1) and (1,-4) which can only appear in affine type BC
    if exc_labels41:
        # tests a two-vertex quiver to see if it is of type ['BC',1,1]
        if len(exc_labels41) == 1 and dict_in_out[exc_labels41[0][0]][2] == dict_in_out[exc_labels41[0][1]][2] == 1:
            return QuiverMutationType(['BC',1,1])
        # test if quiver contains a triangle T with edges [ (0, 1, (2, -1)), (2, 0, (2,-1)), (1, 2, (1, -4)) ] or [ (0, 1, (1, -2)), (2, 0, (1,-2)), (1, 2, (4, -1)) ].
        if len(exc_labels41) == 1 and len(exc_labels) == 2:
            bool2 = exc_labels41[0][2] == (4,-1) and exc_labels[0][2] == exc_labels[1][2] == (1,-2)
            bool3 = exc_labels41[0][2] == (1,-4) and exc_labels[0][2] == exc_labels[1][2] == (2,-1)
            if bool2 or bool3:
                v1, v2, label = exc_labels41[0]
                label1, label2 = exc_labels
                # delete the two vertices associated to the edge with label (1,-4) or (4,-1) and test if the rest of the quiver is of type A.
                # the third vertex of the triangle T should be a connecting_vertex.
                if label1[1] == label2[0] and label2[1] == v1 and v2 == label1[0] and dict_in_out[v1][2] == dict_in_out[v2][2] == 2:
                    _reset_dg( dg, vertices, dict_in_out, [v1,v2] )
                    return _check_special_BC_cases( dg, n, ['BC'], [1], ['A'], [[label1[1]]] )
                elif label1[0] == label2[1] and label1[1] == v1 and v2 == label2[0] and dict_in_out[v1][2] == dict_in_out[v2][2] == 2:
                    _reset_dg( dg, vertices, dict_in_out, [v1,v2] )
                    return _check_special_BC_cases( dg, n, ['BC'], [1], ['A'], [[label1[0]]] )
                else:
                    return _false_return()
            else:
                return _false_return()
        else:
            return _false_return()

    # the program now performs further tests in the case that there are no edges of type (1,-4) nor (4,-1)

    # first test for affine type C: if there are 4 exceptional labels, test if both belong to triangles with leaves
    if len( exc_labels ) == 4:
        exc_labels12 = [labl for labl in exc_labels if labl[2] == (1, -2)]
        exc_labels21 = [labl for labl in exc_labels if labl[2] == (2, -1)]
        # check that we have two labels of one kind and one label of the other
        if len( exc_labels12 ) != 2 or len( exc_labels21 ) != 2:
            return _false_return()

        label121 = exc_labels12[0]
        label122 = exc_labels12[1]
        label211 = exc_labels21[0]
        label212 = exc_labels21[1]

        # affine type B
        if label211[1] == label121[0] and label212[1] == label122[0]:
            pass
        elif label212[1] == label121[0] and label211[1] == label122[0]:
            label211, label212 = label212, label211
        # affine type C
        elif label121[1] == label211[0] and label122[1] == label212[0]:
            pass
        elif label122[1] == label211[0] and label121[1] == label212[0]:
            label211, label212 = label212, label211
        # affine type BC
        elif label121[1] == label211[0] and label212[1] == label122[0]:
            pass
        elif label121[1] == label212[0] and label211[1] == label122[0]:
            label211, label212 = label212, label211
        elif label122[1] == label211[0] and label212[1] == label121[0]:
            label121, label122 = label122, label121
        elif label122[1] == label212[0] and label211[1] == label121[0]:
            pass
        else:
            return _false_return()

        # tests for which configuration the two (1,-2) and two (2,-1) edges are in.
        bool1 = dg.has_edge(label121[1],label211[0],1) and dict_in_out[label211[1]][0] == dict_in_out[label211[1]][1] == 1
        bool2 = dg.has_edge(label122[1],label212[0],1) and dict_in_out[label212[1]][0] == dict_in_out[label212[1]][1] == 1
        bool12 = not ( label121[1] == label122[1] and label211[0] == label212[0] )
        bool3 = dg.has_edge(label211[1],label121[0],1) and dict_in_out[label121[1]][0] == dict_in_out[label121[1]][1] == 1
        bool4 = dg.has_edge(label212[1],label122[0],1) and dict_in_out[label122[1]][0] == dict_in_out[label122[1]][1] == 1
        bool34 = not ( label211[1] == label212[1] and label121[0] == label122[0] )
        bool5 = dg.has_edge(label211[1],label121[0],1) and dict_in_out[label121[1]][0] == dict_in_out[label121[1]][1] == 1
        bool6 = dg.has_edge(label122[1],label212[0],1) and dict_in_out[label212[1]][0] == dict_in_out[label212[1]][1] == 1
        bool56 = not ( label211[1] == label122[1] and label121[0] == label212[0] )
        bool7 = dg.has_edge(label212[1],label122[0],1) and dict_in_out[label122[1]][0] == dict_in_out[label122[1]][1] == 1
        bool8 = dg.has_edge(label121[1],label211[0],1) and dict_in_out[label211[1]][0] == dict_in_out[label211[1]][1] == 1
        bool78 = not ( label212[1] == label121[1] and label122[0] == label211[0] )

        nb1 = len( set(dg.neighbors(label121[1])).intersection(dg.neighbors(label211[0])) ) <= 1
        nb2 = len( set(dg.neighbors(label122[1])).intersection(dg.neighbors(label212[0])) ) <= 1
        nb3 = len( set(dg.neighbors(label211[1])).intersection(dg.neighbors(label121[0])) ) <= 1
        nb4 = len( set(dg.neighbors(label212[1])).intersection(dg.neighbors(label122[0])) ) <= 1

        if bool1 and bool2 and bool12 and nb1 and nb2:
            v1,v2 = label211[1],label212[1]
            _reset_dg( dg, vertices, dict_in_out, [v1,v2] )
            return _check_special_BC_cases( dg, n, ['CC'], [1], ['A'] )
        if bool3 and bool4 and bool34 and nb3 and nb4:
            v1,v2 = label121[1],label122[1]
            _reset_dg( dg, vertices, dict_in_out, [v1,v2] )
            return _check_special_BC_cases( dg, n, ['BB'], [1], ['A'] )
        elif bool5 and bool6 and bool56 and nb2 and nb3:
            v1,v2 = label121[1],label212[1]
            _reset_dg( dg, vertices, dict_in_out, [v1,v2] )
            return _check_special_BC_cases( dg, n, ['BC'], [1], ['A'] )
        elif bool7 and bool8 and bool78 and nb1 and nb4:
            v1,v2 = label122[1],label211[1]
            _reset_dg( dg, vertices, dict_in_out, [v1,v2] )
            return _check_special_BC_cases( dg, n, ['BC'], [1], ['A'] )
        else:
            return _false_return()

    # first test for affine type C: if there are three exceptional labels, we must be in both cases below of the same construction
    elif len( exc_labels ) == 3:
        exc_labels12 = [labl for labl in exc_labels if labl[2] == (1, -2)]
        exc_labels21 = [labl for labl in exc_labels if labl[2] == (2, -1)]
        # check that we have two labels of one kind and one label of the other
        if exc_labels12 == [] or exc_labels21 == []:
            return _false_return()
        if len( exc_labels12 ) == 2:
            label1,label2 = exc_labels12
            label3 = exc_labels21[0]
            if dict_in_out[label2[0]][2] == 1 or dict_in_out[label2[1]][2] == 1:
                label1, label2 = label2, label1
            if dict_in_out[label1[0]][2] == 1:
                if label2[1] == label3[0] and dict_in_out[label2[1]][2] == 2 and dg.has_edge(label3[1],label2[0],1):
                    v1,v2 = label3[1],label2[0]
                    _reset_dg( dg, vertices, dict_in_out, [label2[1]] )
                    if len( set(dg.neighbors_out(v2)).intersection(dg.neighbors_in(v1)) ) > 0:
                        return _false_return()
                    elif len( set(dg.neighbors_out(v1)).intersection(dg.neighbors_in(v2)) ) > 0:
                        return _check_special_BC_cases( dg, n, ['BC'],[1],['A'],[[v1,v2]] )
                    else:
                        return _check_special_BC_cases( dg, n, ['BC'],[1],['A'] )
                elif label3[1] == label2[0] and dict_in_out[label3[1]][2] == 2 and dg.has_edge(label2[1],label3[0],1):
                    v1,v2 = label2[1],label3[0]
                    _reset_dg( dg, vertices, dict_in_out, [label3[1]] )
                    if len( set(dg.neighbors_out(v2)).intersection(dg.neighbors_in(v1)) ) > 0:
                        return _false_return()
                    elif len( set(dg.neighbors_out(v1)).intersection(dg.neighbors_in(v2)) ) > 0:
                        return _check_special_BC_cases( dg, n, ['CC'],[1],['A'],[[v1,v2]] )
                    else:
                        return _check_special_BC_cases( dg, n, ['CC'],[1],['A'] )
                else:
                    return _false_return()
            elif dict_in_out[label1[1]][2] == 1:
                if label3[1] == label2[0] and dict_in_out[label3[1]][2] == 2 and dg.has_edge(label2[1],label3[0],1):
                    v1,v2 = label2[1],label3[0]
                    _reset_dg( dg, vertices, dict_in_out, [label3[1]] )
                    if len( set(dg.neighbors_out(v2)).intersection(dg.neighbors_in(v1)) ) > 0:
                        return _false_return()
                    elif len( set(dg.neighbors_out(v1)).intersection(dg.neighbors_in(v2)) ) > 0:
                        return _check_special_BC_cases( dg, n, ['BC'],[1],['A'],[[v1,v2]] )
                    else:
                        return _check_special_BC_cases( dg, n, ['BC'],[1],['A'] )
                elif label2[1] == label3[0] and dict_in_out[label2[1]][2] == 2 and dg.has_edge(label3[1],label2[0],1):
                    v1,v2 = label3[1],label2[0]
                    _reset_dg( dg, vertices, dict_in_out, [label2[1]] )
                    if len( set(dg.neighbors_out(v2)).intersection(dg.neighbors_in(v1)) ) > 0:
                        return _false_return()
                    elif len( set(dg.neighbors_out(v1)).intersection(dg.neighbors_in(v2)) ) > 0:
                        return _check_special_BC_cases( dg, n, ['BB'],[1],['A'],[[v1,v2]] )
                    else:
                        return _check_special_BC_cases( dg, n, ['BB'],[1],['A'] )
                else:
                    return _false_return()
            elif label1[1] == label2[1] == label3[0] and dict_in_out[label1[1]][2] == 3 and dg.has_edge(label3[1],label1[0],1) and dg.has_edge(label3[1],label2[0],1) and dict_in_out[label2[0]][2] == dict_in_out[label1[0]][2] == 2:
                _reset_dg( dg, vertices, dict_in_out, [label1[1]] )
                return _check_special_BC_cases( dg, n, ['BD'],[1],['D'] )
            elif label1[0] == label2[0] == label3[1] and dict_in_out[label1[0]][2] == 3 and dg.has_edge(label1[1],label3[0],1) and dg.has_edge(label2[1],label3[0],1) and dict_in_out[label2[1]][2] == dict_in_out[label1[1]][2] == 2:
                _reset_dg( dg, vertices, dict_in_out, [label1[0]] )
                return _check_special_BC_cases( dg, n, ['CD'],[1],['D'] )
            else:
                return _false_return()
        elif len( exc_labels21 ) == 2:
            label1,label2 = exc_labels21
            label3 = exc_labels12[0]
            if dict_in_out[label2[0]][2] == 1 or dict_in_out[label2[1]][2] == 1:
                label1, label2 = label2, label1
            if dict_in_out[label1[1]][2] == 1:
                if label2[1] == label3[0] and dict_in_out[label2[1]][2] == 2 and dg.has_edge(label3[1],label2[0],1):
                    v1, v2 = label3[1], label2[0]
                    _reset_dg( dg, vertices, dict_in_out, [label2[1]] )
                    if len( set(dg.neighbors_out(v2)).intersection(dg.neighbors_in(v1)) ) > 0:
                        return _false_return()
                    elif len( set(dg.neighbors_out(v1)).intersection(dg.neighbors_in(v2)) ) > 0:
                        return _check_special_BC_cases( dg, n, ['CC'],[1],['A'],[[v1,v2]] )
                    else:
                        return _check_special_BC_cases( dg, n, ['CC'],[1],['A'] )
                elif label3[1] == label2[0] and dict_in_out[label3[1]][2] == 2 and dg.has_edge(label2[1],label3[0],1):
                    v1, v2 = label2[1], label3[0]
                    _reset_dg( dg, vertices, dict_in_out, [label3[1]] )
                    if len( set(dg.neighbors_out(v2)).intersection(dg.neighbors_in(v1)) ) > 0:
                        return _false_return()
                    elif len( set(dg.neighbors_out(v1)).intersection(dg.neighbors_in(v2)) ) > 0:
                        return _check_special_BC_cases( dg, n, ['BC'],[1],['A'],[[v1,v2]] )
                    else:
                        return _check_special_BC_cases( dg, n, ['BC'],[1],['A'] )
                else:
                    return _false_return()
            elif dict_in_out[label1[0]][2] == 1:
                if label3[1] == label2[0] and dict_in_out[label3[1]][2] == 2 and dg.has_edge(label2[1],label3[0],1):
                    v1, v2 = label2[1], label3[0]
                    _reset_dg( dg, vertices, dict_in_out, [label3[1]] )
                    if len( set(dg.neighbors_out(v2)).intersection(dg.neighbors_in(v1)) ) > 0:
                        return _false_return()
                    elif len( set(dg.neighbors_out(v1)).intersection(dg.neighbors_in(v2)) ) > 0:
                        return _check_special_BC_cases( dg, n, ['BB'],[1],['A'],[[v1,v2]] )
                    else:
                        return _check_special_BC_cases( dg, n, ['BB'],[1],['A'] )
                elif label2[1] == label3[0] and dict_in_out[label2[1]][2] == 2 and dg.has_edge(label3[1],label2[0],1):
                    v1, v2 = label3[1], label2[0]
                    _reset_dg( dg, vertices, dict_in_out, [label2[1]] )
                    if len( set(dg.neighbors_out(v2)).intersection(dg.neighbors_in(v1)) ) > 0:
                        return _false_return()
                    elif len( set(dg.neighbors_out(v1)).intersection(dg.neighbors_in(v2)) ) > 0:
                        return _check_special_BC_cases( dg, n, ['BC'],[1],['A'],[[v1,v2]] )
                    else:
                        return _check_special_BC_cases( dg, n, ['BC'],[1],['A'] )
                else:
                    return _false_return()
            elif label1[0] == label2[0] == label3[1] and dict_in_out[label1[0]][2] == 3 and dg.has_edge(label1[1],label3[0],1) and dict_in_out[label1[1]][2] == 2 and dg.has_edge(label2[1],label3[0],1) and dict_in_out[label2[1]][2] == 2:
                _reset_dg( dg, vertices, dict_in_out, [label3[1]] )
                return _check_special_BC_cases( dg, n, ['BD'],[1],['D'] )
            elif label1[1] == label2[1] == label3[0] and dict_in_out[label3[0]][2] == 3 and dg.has_edge(label3[1],label1[0],1) and dict_in_out[label1[0]][2] == 2 and dg.has_edge(label3[1],label2[0],1) and dict_in_out[label2[0]][2] == 2:
                _reset_dg( dg, vertices, dict_in_out, [label3[0]] )
                return _check_special_BC_cases( dg, n, ['CD'],[1],['D'] )
            else:
                return _false_return()

    # first test for finite types B and C: if there are two exceptional labels, they must belong to an oriented triangle and the vertex between must be a leaf
    # first test for affine type C: if there are two exceptional labels, they must belong to leaves
    # first test for affine type B: if there are two exceptional labels, they must be...
    elif len( exc_labels ) == 2:
        label1, label2 = exc_labels
        if label1[1] == label2[0]:
            pass
        elif label2[1] == label1[0]:
            label1, label2 = label2, label1
        else:
            # the exceptional case in affine type BC_2 is checked
            if label2[2] == (1, -2) and label1[2] == (2, -1):
                label1, label2 = label2, label1
            if label1[2] == (1, -2) and label2[2] == (2, -1):
                if label1[1] == label2[1] and dict_in_out[label1[1]][2] == 2 and dict_in_out[label1[0]][2] == 1 and dict_in_out[label2[0]][2] == 1:
                    return QuiverMutationType(['BC',2,1])
                elif label1[0] == label2[0] and dict_in_out[label1[0]][2] == 2 and dict_in_out[label1[1]][2] == 1 and dict_in_out[label2[1]][2] == 1:
                    return QuiverMutationType(['BC',2,1])
            # the cases in affine type B/C are checked where the exceptional labels connect to leaves
            v11, v12, label1 = label1
            v21, v22, label2 = label2
            if dict_in_out[v11][2] == 1:
                in_out1 = 'out'
            elif dict_in_out[v12][2] == 1:
                in_out1 = 'in'
            else:
                return _false_return()
            if dict_in_out[v21][2] == 1:
                in_out2 = 'out'
            elif dict_in_out[v22][2] == 1:
                in_out2 = 'in'
            else:
                return _false_return()
            if label1 == label2:
                if in_out1 == in_out2 == 'in':
                    if label1 == (1,-2):
                        return _check_special_BC_cases( dg, n, ['BB'],[1],['A'] )
                    else:
                        return _check_special_BC_cases( dg, n, ['CC'],[1],['A'] )
                elif in_out1 == in_out2 == 'out':
                    if label1 == (1,-2):
                        return _check_special_BC_cases( dg, n, ['CC'],[1],['A'] )
                    else:
                        return _check_special_BC_cases( dg, n, ['BB'],[1],['A'] )
                else:
                    return _check_special_BC_cases( dg, n, ['BC'],[1],['A'] )
            else:
                if in_out1 == in_out2:
                    return _check_special_BC_cases( dg, n, ['BC'],[1],['A'] )
                else:
                    if label1 == (1,-2):
                        if in_out1 == 'in':
                            return _check_special_BC_cases( dg, n, ['BB'],[1],['A'] )
                        else:
                            return _check_special_BC_cases( dg, n, ['CC'],[1],['A'] )
                    else:
                        if in_out1 == 'in':
                            return _check_special_BC_cases( dg, n, ['CC'],[1],['A'] )
                        else:
                            return _check_special_BC_cases( dg, n, ['BB'],[1],['A'] )

        v1,v,label1 = label1
        v,v2,label2 = label2
        if dg.has_multiple_edges():
            if all( edge == (v2,v1,1) for edge in dg.multiple_edges() ):
                if dict_in_out[v2][2] == dict_in_out[v1][2] == 3:
                    _reset_dg( dg, vertices, dict_in_out, [v1,v2] )
                    if label1 == (1,-2) and label2 == (2,-1):
                        return _check_special_BC_cases( dg, n, ['CC'],[1],['A'],[[v]] )
                    elif label1 == (2,-1) and label2 == (1,-2):
                        return _check_special_BC_cases( dg, n, ['BB'],[1],['A'],[[v]] )
                    else:
                        return _false_return()
                elif dict_in_out[v][0] == dict_in_out[v][1] == 1:
                    dg.remove_multiple_edges()
                    dg = DiGraph( dg )
                    _reset_dg( dg, vertices, dict_in_out, [v] )
                    if dict_in_out[v1][0] == dict_in_out[v1][1] == dict_in_out[v2][0] == dict_in_out[v2][1] == 1 and next(dg.neighbor_out_iterator(v1)) == next(dg.neighbor_in_iterator(v2)):
                        if label1 == (2,-1) and label2 == (1,-2):
                            return _check_special_BC_cases( dg, n, ['CD'],[1],['A'] )
                        elif label1 == (1,-2) and label2 == (2,-1):
                            return _check_special_BC_cases( dg, n, ['BD'],[1],['A'] )
                    else:
                        return _false_return()
                else:
                    return _false_return()
            else:
                return _false_return()
        elif not dict_in_out[v][0] == 1 or not dict_in_out[v][1] == 1:
            return _false_return()
        else:
            if dg.has_edge(v2,v1,1):
                nr_same_neighbors = len( set(dg.neighbors_out(v1)).intersection(dg.neighbors_in(v2)) )
                nr_other_neighbors = len( set(dg.neighbors_out(v2)).intersection(dg.neighbors_in(v1)) )
                nr_contained_cycles = len([ cycle for cycle, is_oriented in _all_induced_cycles_iter( dg ) if v1 in flatten(cycle) and v2 in flatten(cycle) ] )
                if nr_same_neighbors + nr_other_neighbors + nr_contained_cycles > 2:
                    return _false_return()
                if label1 == (2,-1) and label2 == (1,-2):
                    if n == 4 and (nr_same_neighbors == 2 or nr_other_neighbors == 1):
                        return QuiverMutationType(['CD',n-1,1])
                    # checks for affine A
                    if nr_same_neighbors + nr_other_neighbors > 1:
                        mt_tmp = _check_special_BC_cases( dg, n, ['C','CD'],[None,None],['A','D'],[[],[v]] )
                    else:
                        _reset_dg( dg, vertices, dict_in_out, [v] )
                        mt_tmp = _check_special_BC_cases( dg, n, ['C','CD'],[None,None],['A','D'] )
                    if mt_tmp == 'unknown':
                        dg.delete_edges([[v2,v1],[v1,v],[v,v2]])
                        dg.add_edges([[v1,v2,1],[v,v1,1],[v2,v,1]])
                        if nr_same_neighbors + nr_other_neighbors > 1:
                            #_reset_dg( dg, vertices, dict_in_out, [v] )
                            return _check_special_BC_cases( dg, n, ['CD'],[None],['D'],[[v]] )
                        else:
                            return _check_special_BC_cases( dg, n, ['CD'],[None],['D'] )
                    else:
                        return mt_tmp
                elif label1 == (1,-2) and label2 == (2,-1):
                    if n == 4 and (nr_same_neighbors == 2 or nr_other_neighbors == 1):
                        return QuiverMutationType(['BD',n-1,1])
                    # checks for affine A
                    if nr_same_neighbors + nr_other_neighbors > 1:
                        mt_tmp = _check_special_BC_cases( dg, n, ['B','BD'],[None,None],['A','D'],[[],[v]] )
                    else:
                        _reset_dg( dg, vertices, dict_in_out, [v] )
                        mt_tmp = _check_special_BC_cases( dg, n, ['B','BD'],[None,None],['A','D'] )
                    if mt_tmp == 'unknown':
                        dg.delete_edges([[v2,v1],[v1,v],[v,v2]])
                        dg.add_edges([[v1,v2,1],[v,v1,1],[v2,v,1]])
                        if nr_same_neighbors + nr_other_neighbors > 1:
                            #_reset_dg( dg, vertices, dict_in_out, [v] )
                            return _check_special_BC_cases( dg, n, ['BD'],[None],['D'],[[v]] )
                        else:
                            return _check_special_BC_cases( dg, n, ['BD'],[None],['D'] )
                    else:
                        return mt_tmp
                else:
                    return _false_return()
            elif dict_in_out[v1][2] == 1 and dict_in_out[v2][2] == 1:
                if label1 == (1,-2) and label2 == (1,-2):
                    return _check_special_BC_cases( dg, n, ['BC'],[1],['A'] )
                elif label1 == (2,-1) and label2 == (2,-1):
                    return _check_special_BC_cases( dg, n, ['BC'],[1],['A'] )
                elif label1 == (1,-2) and label2 == (2,-1):
                    return _check_special_BC_cases( dg, n, ['CC'],[1],['A'] )
                elif label1 == (2,-1) and label2 == (1,-2):
                    return _check_special_BC_cases( dg, n, ['BB'],[1],['A'] )
                else:
                    return _false_return()
            elif dict_in_out[v][0] == dict_in_out[v][1] == 1 and dict_in_out[v1][0] == dict_in_out[v1][1] == 1 and dict_in_out[v2][0] == dict_in_out[v2][1] == 1:
                _reset_dg( dg, vertices, dict_in_out, [v] )
                if n == 4 and ( label1, label2 ) == ( (2,-1), (1,-2) ):
                    return _check_special_BC_cases( dg, n, ['CD'],[1],['A'] )
                elif n > 4 and ( label1, label2 ) == ( (2,-1), (1,-2) ):
                    return _check_special_BC_cases( dg, n, ['CD'],[1],['D'] )
                elif n == 4 and ( label1, label2 ) == ( (1,-2), (2,-1) ):
                    return _check_special_BC_cases( dg, n, ['BD'],[1],['A'] )
                elif n > 4 and ( label1, label2 ) == ( (1,-2), (2,-1) ):
                    return _check_special_BC_cases( dg, n, ['BD'],[1],['D'] )
                else:
                    return _false_return()
            else:
                return _false_return()

    # second tests for finite types B and C: if there is only one exceptional label, it must belong to a leaf
    # also tests for affine type B: this exceptional label must belong to a leaf of a type D quiver
    elif len( exc_labels ) == 1:
        label = exc_labels[0]
        v_out = label[0]
        v_in = label[1]
        label = label[2]
        if label == (1,-2):
            if dict_in_out[ v_in ][0] == 1 and dict_in_out[ v_in ][1] == 0:
                #_reset_dg( dg, vertices, dict_in_out, [v_in] )
                return _check_special_BC_cases( dg, n, ['B','BD'],[None,1],['A','D'],[[v_in],[v_in]] )
            elif dict_in_out[ v_out ][0] == 0 and dict_in_out[ v_out ][1] == 1:
                #_reset_dg( dg, vertices, dict_in_out, [v_out] )
                return _check_special_BC_cases( dg, n, ['C','CD'],[None,1],['A','D'],[[v_out],[v_out]] )
            else:
                return _false_return()
        elif label == (2,-1):
            if dict_in_out[ v_out ][0] == 0 and dict_in_out[ v_out ][1] == 1:
                #_reset_dg( dg, vertices, dict_in_out, [v_out] )
                return _check_special_BC_cases( dg, n, ['B','BD'],[None,1],['A','D'],[[v_out],[v_out]] )
            elif dict_in_out[ v_in ][0] == 1 and dict_in_out[ v_in ][1] == 0:
                #_reset_dg( dg, vertices, dict_in_out, [v_in] )
                return _check_special_BC_cases( dg, n, ['C','CD'],[None,1],['A','D'],[[v_in],[v_in]] )
            else:
                return _false_return()

    # if no edges of type (1,-2) nor (2,-1), then tests for type A, affine A, or D.
    return _connected_mutation_type_AAtildeD(dg)


def _connected_mutation_type_AAtildeD(dg, ret_conn_vert=False):
    """
    Return mutation type of ClusterQuiver(dg) for DiGraph dg if it is
    of type finite A, affine A, or finite D.

    For all other types (including affine D), outputs 'unknown'

    See [BPRS2009]_ and [Vat2008]_ (by Vatne) for theoretical details.

    .. TODO::

        Improve this algorithm to also recognize affine D.

    INPUT:

    - ``ret_conn_vert`` -- boolean (default: ``False``); if ``True``,
      returns 'connecting vertices', technical information that is
      used in the algorithm

    A brief description of the algorithm::

         Looks for a long_cycle (of length >= 4) in the digraph dg.  If there is more than one than the mutation_type is 'unknown'.
         Otherwise, checks if each edge of long_cycle connects to a type A quiver.  If so, then ClusterQuiver(dg) is of type D or affine A.

         If there is no long_cycle, then checks that there are no multiple edges, all triangles are oriented,
         no vertices of valence higher than 4, that vertices of valence 4 are incident to two oriented triangles,
         and that a vertex of valence 3 has exactly two incident arrows as part of an oriented triangle.
         All these checks ensures that ClusterQuiver(dg) is of type A except for three exceptions that are also checked.

    EXAMPLES::

        sage: # needs sage.modules
        sage: from sage.combinat.cluster_algebra_quiver.mutation_type import _connected_mutation_type_AAtildeD
        sage: Q = ClusterQuiver(['A',[7,0],1]); Q.mutate([0,1,4])
        sage: _connected_mutation_type_AAtildeD(Q.digraph(),ret_conn_vert=True)
        [['D', 7], [0, 4]]
        sage: Q2 = ClusterQuiver(['A',[5,2],1]); Q2.mutate([4,5])
        sage: _connected_mutation_type_AAtildeD(Q2.digraph() )
        ['A', [2, 5], 1]
        sage: Q3 = ClusterQuiver(['E',6]); Q3.mutate([5,2,1])
        sage: _connected_mutation_type_AAtildeD(Q3.digraph(),ret_conn_vert=True)
        'unknown'
    """
    # naming the vertices
    vertices = list(dg)
    n = dg.order()

    # Test if ClusterQuiver(dg) is of type D_n Type 1, i.e. A_{n-2} plus two leaves
    if n > 3:
        # check if any vertices have a neighborhood with two leaves.  If so, prune the two leaves and retest on this smaller digraph.
        # note that this step is unnecessary for digraphs with fewer than four vertices.
        for v in vertices:
            dead_neighbors = [ v_n for v_n in dg.neighbors(v) if dg.degree(v_n) == 1 ]
            if len( dead_neighbors ) >= 2:
                dg_tmp = DiGraph( dg )
                dg_tmp.delete_vertices( dead_neighbors[:2] )
                type_tmp = _connected_mutation_type_AAtildeD( dg_tmp, ret_conn_vert=True )
                if type_tmp == 'unknown':
                    return _false_return()
                # if smaller digraph is of finite A type with v as a 'connecting vertex', then glueing back the two leaves yields type finite D.
                if type_tmp[0].letter() == 'A' and type_tmp[0].is_finite():
                    if v in type_tmp[1]:
                        type_tmp[1].remove(v)
                        if n == 4:
                            type_tmp[1].extend(dead_neighbors[:2])
                        if ret_conn_vert:
                            return [ QuiverMutationType( ['D',n] ), type_tmp[1] ]
                        else:
                            return QuiverMutationType( ['D',n] )
                    # note that if v is not a 'connecting vertex' then we make no conclusion either way.
                else:
                    return _false_return(3)
        # Test if ClusterQuiver(dg) is of type D_n Type 2 or 3, i.e. two type A quivers plus a 4-cycle or triangulated square glued together
        # at 'connecting vertices'.

        # Exception 1 (Type 2 of D_n)
        exception_graph1 = DiGraph()
        exception_graph1.add_edges([(0,1),(1,2),(2,3),(3,0)])

        # Exception 2 (Type 3 of D_n)
        exception_graph2 = DiGraph()
        exception_graph2.add_edges([(0,1),(1,2),(0,3),(3,2),(2,0)])

        # Let c_1 be a pair of 2-valent vertices and c_2 be a pair of two other vertices.
        # If together, they make an induced 4-cycle and deleting c_1 yields two connected components,
        # then retest of both components.  If still connected after deleting c_1, then return 'unknown'.

        # If on the other hand, (c1 and c2) is isomorphic to a triangulated square, then
        # delete c1.  This ensures that c2 is an edge of the triangulated square, and we delete
        # it regardless of orientation.  Then check if the digraph has exactly two connected
        # components, and again this testing method is rerun on both components.

        for c1 in Combinations( [ vertex for vertex in vertices if dg.degree(vertex) == 2], 2 ):
            del_vertices = list( vertices )
            del_vertices.remove( c1[0] )
            del_vertices.remove( c1[1] )
            for c2 in Combinations( del_vertices, 2 ):
                comb = c1 + c2
                sg = dg.subgraph( comb )

                # Exception 1 case (4-cycle):
                edges = sg.edges(sort=True, labels=False)
                if (c1[0], c1[1]) not in edges and (c1[1], c1[0]) not in edges and sg.is_isomorphic(exception_graph1):
                    dg_tmp = DiGraph(dg)
                    dg_tmp.delete_vertices(c1)

                    components = dg_tmp.connected_components(sort=False)
                    # if not len(components) == 2:
                    if len(components) != 2:
                        return _false_return(4)
                    else:
                        dg_tmp1 = dg_tmp.subgraph( components[0] )
                        type_tmp1 = _connected_mutation_type_AAtildeD( dg_tmp1, ret_conn_vert=True )
                        dg_tmp2 = dg_tmp.subgraph( components[1] )
                        type_tmp2 = _connected_mutation_type_AAtildeD( dg_tmp2, ret_conn_vert=True )

                        if type_tmp1 == 'unknown' or type_tmp2 == 'unknown':
                            return _false_return()

                        # Assuming that the two components are recognized, initialize this in a format it can be returned as output
                        type_tmp = []
                        type_tmp.append( [ type_tmp1[0], type_tmp2[0] ] )
                        type_tmp[0].sort(key=str)
                        type_tmp.append( type_tmp1[1] + type_tmp2[1] )
                        type_tmp[1].sort(key=str)

                        # Need to make sure the two vertices in c2 are both 'connecting vertices'.
                        if not set(c2).issubset(type_tmp[1]):
                            return _false_return(5)

                        if type_tmp[0][0].letter() == 'A' and type_tmp[0][0].is_finite() and type_tmp[0][1].letter() == 'A' and type_tmp[0][1].is_finite():
                            if ret_conn_vert:
                                type_tmp[1].extend(c1)
                                #type_tmp[1].remove(c2[0])
                                #type_tmp[1].remove(c2[1])
                                return [ QuiverMutationType( ['D',n] ), type_tmp[1] ]
                            else:
                                return QuiverMutationType( ['D',n] )

                # Exception 2 case (triangulated square):
                if sg.is_isomorphic( exception_graph2 ):
                    dg_tmp = DiGraph( dg )
                    dg_tmp.delete_vertices( c1 )
                    if tuple( c2 ) in dg_tmp.edges(sort=True, labels=False):
                        dg_tmp.delete_edge( tuple( c2 ) )
                    else:
                        c2.reverse()
                        dg_tmp.delete_edge( tuple( c2 ) )
                    components = dg_tmp.connected_components(sort=False)
                    if len(components) != 2:
                        return _false_return(7)
                    else:
                        dg_tmp1 = dg_tmp.subgraph( components[0] )
                        type_tmp1 = _connected_mutation_type_AAtildeD( dg_tmp1, ret_conn_vert=True )

                        if type_tmp1 == 'unknown':
                            return _false_return()
                        dg_tmp2 = dg_tmp.subgraph( components[1] )
                        type_tmp2 = _connected_mutation_type_AAtildeD( dg_tmp2, ret_conn_vert=True )

                        # Assuming that the two components are recognized, initialize this in
                        # a format it can be returned as output (just as above)
                        type_tmp = []
                        type_tmp.append( [ type_tmp1[0], type_tmp2[0] ] )
                        type_tmp[0].sort(key=str)
                        type_tmp.append( type_tmp1[1] + type_tmp2[1] )
                        type_tmp[1].sort(key=str)
                        if type_tmp2 == 'unknown':
                            return _false_return()
                        if not set(c2).issubset(type_tmp[1]) and len( set(type_tmp[1]).intersection(c2) ) == 1:
                            return _false_return(5.5)
                        if type_tmp[0][0].letter() == 'A' and type_tmp[0][0].is_finite() and type_tmp[0][1].letter() == 'A' and type_tmp[0][1].is_finite():
                            if ret_conn_vert:
                                type_tmp[1].remove(c2[0])
                                type_tmp[1].remove(c2[1])
                                #type_tmp[1].extend(c1)
                                return [ QuiverMutationType( ['D',n] ), type_tmp[1] ]
                            else:
                                return QuiverMutationType( ['D',n] )

    # The following tests are done regardless of the number of vertices in dg.
    # If there are 1, 2, or 3 vertices in dg, we would have skipped above tests and gone directly here.

    # Initialize a long_cycle.
    long_cycle = False

    # test that there is no triple-edge or higher multiplicity and that there is at most one double-edge.
    if dg.has_multiple_edges():
        multiple_edges = dg.multiple_edges(labels=False)
        if len(multiple_edges) > 2:
            return _false_return(14)
        elif len(multiple_edges) == 2:
            # we think of the double-edge as a long_cycle, an unoriented 2-cycle.
            long_cycle = [multiple_edges, ['A', n - 1, 1]]

    # creating a dictionary of in-, out- and total degrees
    dict_in_out = {}
    for v in vertices:
        dict_in_out[v] = (dg.in_degree(v), dg.out_degree(v), dg.degree(v))

    # computing the absolute degree of dg
    abs_deg = max( [ x[2] for x in list( dict_in_out.values() ) ] )

    # edges = dg.edges(sort=True, labels=False )

    # test that no vertex has valency more than 4
    if abs_deg > 4:
        return _false_return(16)
    else:
        # constructing all oriented and unoriented triangles
        trians = _triangles( dg )
        oriented_trians = [ trian[0] for trian in trians if trian[1] ]
        unoriented_trians = [ trian[0] for trian in trians if not trian[1] ]

        oriented_trian_edges = []
        for oriented_trian in oriented_trians:
            oriented_trian_edges.extend( oriented_trian )

        # test that no edge is in more than two oriented triangles
        multiple_trian_edges = []
        for edge in oriented_trian_edges:
            count = oriented_trian_edges.count(edge)
            if count > 2:
                return _false_return(17)
            elif count == 2:
                multiple_trian_edges.append( edge )
        multiple_trian_edges = list(set(multiple_trian_edges))

        # test that there at most three edges appearing in exactly two oriented triangles
        count = len(multiple_trian_edges)
        if count >= 4:
            return _false_return(321)
        # if two edges appearing in exactly two oriented triangles, test that the two edges together
        # determine a unique triangle
        elif count > 1:
            test_triangles = [[tuple(trian) for trian in oriented_trians
                               if edge in trian]
                              for edge in multiple_trian_edges]
            unique_triangle = set.intersection(*map(set, test_triangles))
            if len(unique_triangle) != 1:
                return _false_return(19)
            else:
                # if a long_cycle had previously been found, this unique oriented triangle is a second long_cycle, a contradiction.
                if long_cycle:
                    return _false_return(20)
                else:
                    unique_triangle = unique_triangle.pop()
                    long_cycle = [ unique_triangle, QuiverMutationType( ['D',n] ) ]
        # if one edge appearing in exactly two oriented triangles, test that it is not a double-edge and then
        # test that either the third or fourth vertices (from the oriented triangles) is of degree 2.
        # Then initializes the long_cycle as this triangle including the degree 2 vertex, as long as no other long_cycles.
        elif count == 1 and not dg.has_multiple_edges() and multiple_trian_edges[0] not in dg.multiple_edges():
            multiple_trian_edge = multiple_trian_edges[0]
            neighbors = list(set(dg.neighbors( multiple_trian_edge[0] )).intersection(dg.neighbors( multiple_trian_edge[1] )))
            if dg.degree( neighbors[0] ) == 2:
                unique_triangle = [ multiple_trian_edge, ( multiple_trian_edge[1], neighbors[0] ), ( neighbors[0], multiple_trian_edge[0] ) ]
            elif dg.degree( neighbors[1] ) == 2:
                unique_triangle = [ multiple_trian_edge, ( multiple_trian_edge[1], neighbors[1] ), ( neighbors[1], multiple_trian_edge[0] ) ]
            else:
                return _false_return(201)

            if long_cycle:
                # if a long_cycle had previously been found, then the specified oriented triangle is a second long_cycle, a contradiction.
                return _false_return(202)
            else:
                long_cycle = [ unique_triangle, QuiverMutationType( ['D',n] ) ]

        # there can be at most 1 unoriented triangle and this triangle is the exceptional circle of type A_tilde
        if unoriented_trians:
            if len(unoriented_trians) == 1:
                if long_cycle:
                    return _false_return(21)
                else:
                    long_cycle = [ unoriented_trians[0], ['A',n-1,1] ]
            else:
                return _false_return(22)

        for v in vertices:
            w = dict_in_out[v]
            if w[2] == 4:
                # if a vertex has valency 4 than the 4 neighboring edges must be contained in 2 oriented triangles
                if w[0] != 2:
                    return _false_return(23)
                else:
                    in_neighbors = dg.neighbors_in( v )
                    out_neighbors = dg.neighbors_out( v )
                    if len( out_neighbors ) == 1:
                        out_neighbors.extend(out_neighbors)
                    if len( in_neighbors ) == 1:
                        in_neighbors.extend(in_neighbors)

                    if (in_neighbors[0], v) not in oriented_trian_edges:
                        return _false_return(24)
                    elif (in_neighbors[1], v) not in oriented_trian_edges:
                        return _false_return(25)
                    elif (v, out_neighbors[0]) not in oriented_trian_edges:
                        return _false_return(26)
                    elif (v, out_neighbors[1]) not in oriented_trian_edges:
                        return _false_return(27)

            # if a vertex has valency 3 than 2 of its neighboring edges must be contained in an oriented triangle and the remaining must not
            elif w[2] == 3:
                if w[0] == 1:
                    in_neighbors = dg.neighbors_in( v )
                    out_neighbors = dg.neighbors_out( v )
                    if (in_neighbors[0],v) not in oriented_trian_edges:
                        return _false_return(28)
                    elif len( out_neighbors ) == 1:
                        if (v,out_neighbors[0]) not in oriented_trian_edges:
                            return _false_return(29)
                    else:
                        if (v,out_neighbors[0]) in oriented_trian_edges and (v,out_neighbors[1]) in oriented_trian_edges:
                            if not long_cycle:
                                return _false_return(30)
                            if not long_cycle[1] == QuiverMutationType(['D',n]):
                                return _false_return(31)
                            if (v, out_neighbors[0]) not in long_cycle[0] and (v, out_neighbors[1]) not in long_cycle[0]:
                                return _false_return(32)
                        if (v,out_neighbors[0]) not in oriented_trian_edges and (v,out_neighbors[1]) not in oriented_trian_edges:
                            return _false_return(33)
                elif w[0] == 2:
                    in_neighbors = dg.neighbors_in( v )
                    out_neighbors = dg.neighbors_out( v )
                    if (v, out_neighbors[0]) not in oriented_trian_edges:
                        return _false_return(34)
                    elif len( in_neighbors ) == 1:
                        if (in_neighbors[0],v) not in oriented_trian_edges:
                            return _false_return(35)
                    else:
                        if (in_neighbors[0],v) in oriented_trian_edges and (in_neighbors[1],v) in oriented_trian_edges:
                            if not long_cycle:
                                return _false_return(36)
                            if not long_cycle[1] == QuiverMutationType(['D',n]):
                                return _false_return(37)
                            if (in_neighbors[0], v) not in long_cycle[0] and (in_neighbors[1], v) not in long_cycle[0]:
                                return _false_return(38)
                        if (in_neighbors[0], v) not in oriented_trian_edges and (in_neighbors[1], v) not in oriented_trian_edges:
                            return _false_return(39)
                else:
                    return _false_return(40)

    # there can exist at most one larger oriented or unoriented induced cycle
    # if it is oriented, we are in finite type D, otherwise we are in affine type A

    # Above code found long_cycles would be an unoriented 2-cycle or an oriented triangle.
    # The method _all_induced_cycles_iter only looks for induced cycles on 4 or more vertices.

    for cycle, is_oriented in _all_induced_cycles_iter( dg ):
        # if there already was a long_cycle and we found another one, then have a contradiction.
        if long_cycle:
            return _false_return(41)
        # otherwise, we obtain cases depending on whether or not the found long_cycle is oriented.
        elif is_oriented:
            long_cycle = [ cycle, QuiverMutationType(['D',n]) ]
        else:
            long_cycle = [ cycle, ['A',n-1,1] ]
    # if we haven't found a "long_cycle", we are in finite type A
    if not long_cycle:
        long_cycle = [[], QuiverMutationType(['A', n])]

    # The 'connected vertices' are now computed.
    # Attention: 0-1-2 in type A_3 has connecting vertices 0 and 2, while in type D_3 it has connecting vertex 1;
    # this is not caught here.
    if ret_conn_vert:
        connecting_vertices = []
        o_trian_verts = flatten( oriented_trian_edges )
        long_cycle_verts = flatten( long_cycle[0] )
        for v in vertices:
            w = dict_in_out[v]
            # if the quiver consists of only one vertex, it is of type A_1 and the vertex is a connecting vertex
            if w[2] == 0:
                connecting_vertices.append( v )
            # if a vertex is a leaf in a type A quiver, it is a connecting vertex
            elif w[2] == 1:
                connecting_vertices.append( v )
            # if a vertex is of valence two and contained in an oriented 3-cycle, it is a connecting vertex
            elif w[0] == 1 and w[1] == 1:
                if v in o_trian_verts and v not in long_cycle_verts:
                    connecting_vertices.append( v )

    # post-parsing 1: if we are in the affine type A case, the two parameters for the non-oriented long cycle are computed
    if isinstance(long_cycle[1], list) and len( long_cycle[1] ) == 3 and long_cycle[1][0] == 'A' and long_cycle[1][2] == 1:
        tmp = list( long_cycle[0] )
        e = tmp.pop()
        cycle = [e]
        v = e[1]
        while tmp:
            e = next(x for x in tmp if v in x)
            if v == e[0]:
                cycle.append(e)
                v = e[1]
            else:
                v = e[0]
            tmp.remove( e )

        tmp = list( cycle )
        if len( long_cycle[0] ) == 2:
            edge = long_cycle[0][0]
            sg = DiGraph( dg )
            sg. delete_vertices(edge)
            connected_components = sg.connected_components(sort=False)
            cycle = []
            if connected_components:
                cycle.append( ( edge[0], edge[1], len( connected_components[0] ) + 1 ) )
            else:
                cycle.append( ( edge[0], edge[1], 1 ) )
        else:
            for edge in tmp:
                sg = DiGraph( dg )
                sg. delete_vertices(edge)
                connected_components = sg.connected_components(sort=False)
                if len( connected_components ) == 2:
                    #if len( list_intersection( [ connected_components[0], list_substract( long_cycle[0], [edge] )[0] ] ) ) > 0:
                    if len( set(connected_components[0]).intersection( set(long_cycle[0]).difference([edge]).pop() ) ) > 0:
                        cycle.remove(edge)
                        cycle.append( (edge[0],edge[1], len( connected_components[1] ) + 1 ) )
                    else:
                        cycle.remove(edge)
                        cycle.append( (edge[0],edge[1], len( connected_components[0] ) + 1 ) )
                else:
                    cycle.remove(edge)
                    cycle.append( (edge[0],edge[1], 1 ) )
        r = sum(x[2] for x in cycle)
        r = max(r, n - r)
        if ret_conn_vert:
            return [ QuiverMutationType( ['A',[r,n-r],1] ), connecting_vertices ]
        else:
            return QuiverMutationType( ['A',[r,n-r],1] )

    # post-parsing 2: if we are in another type, it is returned
    else:
        if ret_conn_vert:
            return [ long_cycle[1], connecting_vertices ]
        else:
            return long_cycle[1]


@cached_function
def load_data(n, user=True):
    r"""
    Load a dict with keys being tuples representing exceptional
    QuiverMutationTypes, and with values being lists or sets
    containing all mutation equivalent quivers as dig6 data.

    We check

    - the data stored by the user (unless ``user=False`` was given)
    - and the data installed by the optional package ``database_mutation_class``.

    INPUT:

    - ``user`` -- boolean (default: ``True``); whether to look at user
      data. If not, only consider the optional package.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_type import load_data
        sage: load_data(2) # random - depends on the data the user has stored
        {('G', 2): [('AO', (((0, 1), (1, -3)),)), ('AO', (((0, 1), (3, -1)),))]}

    TESTS:

    We test data from the ``database_mutation_class`` optional package::

        sage: load_data(2, user=False)      # optional - database_mutation_class
        {('G', 2): [('AO', (((0, 1), (1, -3)),)), ('AO', (((0, 1), (3, -1)),))]}
        sage: D = load_data(3, user=False)  # optional - database_mutation_class
        sage: sorted(D.items())             # optional - database_mutation_class
        [(('G', 2, -1),
          [('BH?', (((1, 2), (1, -3)),)),
           ('BGO', (((2, 1), (3, -1)),)),
           ('BW?', (((0, 1), (3, -1)),)),
           ('BP?', (((0, 1), (1, -3)),)),
           ('BP_', (((0, 1), (1, -3)), ((2, 0), (3, -1)))),
           ('BP_', (((0, 1), (3, -1)), ((1, 2), (1, -3)), ((2, 0), (2, -2))))]),
         (('G', 2, 1),
          [('BH?', (((1, 2), (3, -1)),)),
           ('BGO', (((2, 1), (1, -3)),)),
           ('BW?', (((0, 1), (1, -3)),)),
           ('BP?', (((0, 1), (3, -1)),)),
           ('BKO', (((1, 0), (3, -1)), ((2, 1), (1, -3)))),
           ('BP_', (((0, 1), (2, -2)), ((1, 2), (1, -3)), ((2, 0), (3, -1))))])]
    """
    from sage.env import DOT_SAGE, SAGE_SHARE

    # we check
    # - if the data is stored by the user, and if this is not the case
    # - if the data is stored by the optional package install
    paths = [Path(SAGE_SHARE)]
    if user:
        paths.append(Path(DOT_SAGE))
    data = {}
    for path in paths:
        file = path / 'cluster_algebra_quiver' / f'mutation_classes_{n}.dig6'
        try:
            with open(file, 'rb') as fobj:
                data_new = pickle.load(fobj)
        except (OSError, FileNotFoundError, pickle.UnpicklingError):
            # File does not exist, corrupt pickle, wrong Python version...
            pass
        else:
            data.update(data_new)
    return data


def _mutation_type_from_data(n, dig6, compute_if_necessary=True):
    r"""
    Return the mutation type from the given dig6 data by looking into
    the precomputed mutation types

    Attention: it is assumed that dig6 is the dig6 data of the
    canonical form of the given quiver!

    EXAMPLES::

        sage: # needs sage.modules
        sage: from sage.combinat.cluster_algebra_quiver.mutation_class import _digraph_to_dig6
        sage: from sage.combinat.cluster_algebra_quiver.mutation_type import _mutation_type_from_data
        sage: from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
        sage: dg = ClusterQuiver(['F',4]).canonical_label().digraph()
        sage: dig6 = _digraph_to_dig6(dg,hashable=True); dig6
        ('CCo?', (((1, 3), (2, -1)),))
        sage: _mutation_type_from_data(4,dig6)
        ['F', 4]
    """
    # we try to load the data from a library
    data = load_data(n)
    # if this didn't work, we construct all exceptional quivers with n vertices
    if compute_if_necessary and data == {}:
        from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import save_quiver_data
        save_quiver_data(n, up_to=False, types='Exceptional', verbose=False)
        load_data.clear_cache()
        data = load_data(n)
    # finally, we check if the given quiver is in one of the exceptional mutation classes
    for mutation_type in data:
        if dig6 in data[ mutation_type ]:
            return QuiverMutationType( mutation_type )
    return 'unknown'


def _mutation_type_test(n):
    """
    Test all quivers (of the given types) of rank n to check that
    mutation_type() works.

    Affine type D does not return ``True`` since this test is not implemented.

    EXAMPLES::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_type import _mutation_type_test

        sage: _mutation_type_test(2) # long time
        True ('A', (1, 1), 1)
        True ('A', 2)
        True ('B', 2)
        True ('BC', 1, 1)
        True ('G', 2)

        sage: _mutation_type_test(3) # long time
        True ('A', (2, 1), 1)
        True ('A', 3)
        True ('B', 3)
        True ('BB', 2, 1)
        True ('BC', 2, 1)
        True ('C', 3)
        True ('CC', 2, 1)
        True ('G', 2, -1)
        True ('G', 2, 1)

        sage: _mutation_type_test(4) # not tested
        True ('A', (2, 2), 1)
        True ('A', (3, 1), 1)
        True ('A', 4)
        True ('B', 4)
        True ('BB', 3, 1)
        True ('BC', 3, 1)
        True ('BD', 3, 1)
        True ('C', 4)
        True ('CC', 3, 1)
        True ('CD', 3, 1)
        True ('D', 4)
        True ('F', 4)
        True ('G', 2, (1, 1))
        True ('G', 2, (1, 3))
        True ('G', 2, (3, 3))

        sage: _mutation_type_test(5) # not tested
        True ('A', (3, 2), 1)
        True ('A', (4, 1), 1)
        True ('A', 5)
        True ('B', 5)
        True ('BB', 4, 1)
        True ('BC', 4, 1)
        True ('BD', 4, 1)
        True ('C', 5)
        True ('CC', 4, 1)
        True ('CD', 4, 1)
        False ('D', 4, 1)
        True ('D', 5)
        True ('F', 4, -1)
        True ('F', 4, 1)
    """
    from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import _construct_classical_mutation_classes
    from sage.combinat.cluster_algebra_quiver.mutation_class import _dig6_to_digraph
    from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
    data = _construct_classical_mutation_classes( n )
    keys = data.keys()
    for mutation_type in sorted(keys, key=str):
        mt = QuiverMutationType( mutation_type )
        print(all( ClusterQuiver(_dig6_to_digraph(dig6)).mutation_type() == mt for dig6 in data[mutation_type]), mutation_type)
    from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import _construct_exceptional_mutation_classes
    data = _construct_exceptional_mutation_classes( n )
    keys = data.keys()
    for mutation_type in sorted(keys, key=str):
        mt = QuiverMutationType( mutation_type )
        print(all( ClusterQuiver(_dig6_to_digraph(dig6)).mutation_type() == mt for dig6 in data[mutation_type]), mutation_type)


def _random_tests(mt, k, mut_class=None, nr_mut=5):
    """
    Provide random tests to find bugs in the mutation type methods.

    INPUT:

    - ``mt`` something that can be turned into a ``QuiverMutationType``
    - ``k`` -- integer; the number of tests performed for each quiver of rank ``n``
    - ``mut_class`` -- if given, this mutation class is used
    - ``nr_mut`` -- integer (default: 5); the number of mutations performed before
      testing

    The idea of this random test is to start with a mutation type
    and compute its mutation class (or have this class given). Now,
    every quiver in this mutation class is slightly changed in order
    to obtain a matrix of the same type or something very similar.
    Now, the new type is computed and checked if it stays stable for
    ``nr_mut``'s many mutations.

    TESTS::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_type import _random_tests
        sage: _random_tests( ['A',3], 1)                                                # needs sage.modules
        testing ['A', 3]
    """
    from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
    from sage.combinat.cluster_algebra_quiver.mutation_class import _dig6_to_matrix, _matrix_to_digraph, _digraph_mutate, _edge_list_to_matrix
    import random
    if mut_class is None:
        mut_class = ClusterQuiver(mt).mutation_class(data_type='dig6')
    print("testing " + str(mt))
    for dig6 in mut_class:
        M_const = _dig6_to_matrix( dig6 )
        nz = [ (i,j) for i,j in M_const.nonzero_positions() if i > j ]
        # performing k tests on the matrix M_const
        for i in range(k):
            M = copy( M_const )
            # every pair M[i,j],M[j,i] is possibly changed
            # while the property of being skew-symmetrizable is kept
            for i,j in nz:
                a,b = M[i,j],M[j,i]
                skew_sym = False
                while not skew_sym:
                    ran = random.randint(1, 2)
                    if ran == 1:
                        M[i, j], M[j, i] = -M[j, i], -M[i, j]
                    elif ran == 2:
                        ran2 = random.randint(1, 8)
                        if ran2 == 1:
                            c, d = 1, -1
                        elif ran2 == 2:
                            c, d = 1, -2
                        elif ran2 == 3:
                            c, d = 2, -1
                        elif ran2 == 4:
                            c, d = 1, -3
                        elif ran2 == 5:
                            c, d = 3, -1
                        elif ran2 == 6:
                            c, d = 2, -2
                        elif ran2 == 7:
                            c, d = 1, -4
                        elif ran2 == 8:
                            c, d = 4, -1
                        M[i, j], M[j, i] = c, d
                    if M.is_skew_symmetrizable(positive=True):
                        skew_sym = True
                    else:
                        M[i, j], M[j, i] = a, b
            # we now have a new matrix M
            # and a new digraph db
            dg = _matrix_to_digraph(M)
            mt = _connected_mutation_type(dg)
            mut = -1
            # we perform nr_mut many mutations
            for k in range(nr_mut):
                # while making sure that we do not mutate back
                mut_tmp = mut
                while mut == mut_tmp:
                    mut = random.randint(0, dg.order() - 1)
                dg_new = _digraph_mutate(dg, mut)
                mt_new = _connected_mutation_type(dg_new)
                if mt != mt_new:
                    print("FOUND ERROR!")
                    print(_edge_list_to_matrix(dg.edges(sort=True),
                                               list(range(dg.order())), []))
                    print("has mutation type " + str(mt) + " while it has mutation type " + str(mt_new) + " after mutating at " + str(mut) + ":")
                    print(_edge_list_to_matrix(dg_new.edges(sort=True),
                                               list(range(dg.order())), []))
                    return dg, dg_new
                else:
                    dg = dg_new


def _random_multi_tests(n, k, nr_mut=5):
    """
    Provide multiple random tests to find bugs in the mutation type methods.

    INPUT:

    - ``n`` -- integer; the rank of the mutation types to test
    - ``k`` -- integer; the number of tests performed for each quiver of rank ``n``
    - ``nr_mut`` -- integer (default: 5); the number of mutations performed before testing

    TESTS::

        sage: from sage.combinat.cluster_algebra_quiver.mutation_type import _random_multi_tests
        sage: _random_multi_tests(2,1)  # not tested
        testing ('A', (1, 1), 1)
        testing ('A', 2)
        testing ('B', 2)
        testing ('BC', 1, 1)

        sage: _random_multi_tests(3,1)  # not tested
        testing ('A', (2, 1), 1)
        testing ('A', 3)
        testing ('B', 3)
        testing ('BB', 2, 1)
        testing ('BC', 2, 1)
        testing ('C', 3)
        testing ('CC', 2, 1)

        sage: _random_multi_tests(4,1)  # not tested
        testing ('A', (2, 2), 1)
        testing ('A', (3, 1), 1)
        testing ('A', 4)
        testing ('B', 4)
        testing ('BB', 3, 1)
        testing ('BC', 3, 1)
        testing ('BD', 3, 1)
        testing ('C', 4)
        testing ('CC', 3, 1)
        testing ('CD', 3, 1)
        testing ('D', 4)
    """
    from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import _construct_classical_mutation_classes
    mutation_classes = _construct_classical_mutation_classes(n)
    for mutation_type in sorted(mutation_classes, key=str):
        _random_tests(mutation_type, k,
                      mut_class=mutation_classes[mutation_type], nr_mut=nr_mut)
