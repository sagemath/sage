"""
Graph-theoretic partition backtrack functions

EXAMPLES::

    sage: import sage.groups.perm_gps.partn_ref.refinement_graphs

REFERENCE:

- [1] McKay, Brendan D. *Practical Graph Isomorphism*. Congressus Numerantium,
  Vol. 30 (1981), pp. 45-87.
"""
# ****************************************************************************
#       Copyright (C) 2006 - 2011 Robert L. Miller <rlmillster@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.groups.perm_gps.partn_ref.data_structures cimport *
from sage.data_structures.bitset_base cimport *
from sage.rings.integer cimport Integer
from sage.graphs.base.sparse_graph cimport SparseGraph
from sage.graphs.base.dense_graph cimport DenseGraph, copy_dense_graph
from sage.groups.perm_gps.partn_ref.double_coset cimport double_coset


def isomorphic(G1, G2, partn, ordering2, dig, use_indicator_function, sparse=False):
    """
    Test whether two graphs are isomorphic.

    EXAMPLES::

        sage: from sage.groups.perm_gps.partn_ref.refinement_graphs import isomorphic

        sage: G = Graph(2)
        sage: H = Graph(2)
        sage: isomorphic(G, H, [[0,1]], [0,1], 0, 1)
        {0: 0, 1: 1}
        sage: isomorphic(G, H, [[0,1]], [0,1], 0, 1)
        {0: 0, 1: 1}
        sage: isomorphic(G, H, [[0],[1]], [0,1], 0, 1)
        {0: 0, 1: 1}
        sage: isomorphic(G, H, [[0],[1]], [1,0], 0, 1)
        {0: 1, 1: 0}

        sage: G = Graph(3)
        sage: H = Graph(3)
        sage: isomorphic(G, H, [[0,1,2]], [0,1,2], 0, 1)
        {0: 0, 1: 1, 2: 2}
        sage: G.add_edge(0,1)
        sage: isomorphic(G, H, [[0,1,2]], [0,1,2], 0, 1)
        False
        sage: H.add_edge(1,2)
        sage: isomorphic(G, H, [[0,1,2]], [0,1,2], 0, 1)
        {0: 1, 1: 2, 2: 0}
    """
    cdef PartitionStack *part
    cdef int *output
    cdef int *ordering
    cdef CGraph G
    cdef GraphStruct GS1 = GraphStruct()
    cdef GraphStruct GS2 = GraphStruct()
    cdef GraphStruct GS
    cdef int i, j, n = -1
    cdef list partition, cell
    cdef bint loops = 0

    from sage.graphs.graph import Graph
    from sage.graphs.generic_graph import GenericGraph
    from copy import copy
    which_G = 1
    for G_in in [G1, G2]:
        if which_G == 1:
            GS = GS1
            first = True
        else:
            GS = GS2
            first = False
        if isinstance(G_in, GenericGraph):
            if G_in.has_loops():
                loops = 1
            if n == -1:
                n = G_in.num_verts()
            elif n != G_in.num_verts():
                return False
            if G_in.vertices(sort=True) != list(range(n)):
                G_in = copy(G_in)
                to = G_in.relabel(return_map=True)
                frm = {}
                for v in to:
                    frm[to[v]] = v
                if first:
                    partition = [[to[v] for v in cell] for cell in partn]
            else:
                if first:
                    partition = partn
                to = list(range(n))
                frm = to
            if sparse:
                G = SparseGraph(n)
            else:
                G = DenseGraph(n)
            if G_in.is_directed():
                for i, j in G_in.edge_iterator(labels=False):
                    G.add_arc(i, j)
            else:
                for i, j in G_in.edge_iterator(labels=False):
                    G.add_arc(i, j)
                    G.add_arc(j, i)
        elif isinstance(G_in, CGraph):
            G = <CGraph> G_in
            if n == -1:
                n = G.num_verts
            elif n != <int>G.num_verts:
                return False
            if not loops:
                for i in range(n):
                    if G.has_arc_unsafe(i,i):
                        loops = 1
            to = {}
            for a in G.verts():
                to[a] = a
            frm = to
            if first:
                partition = partn
        else:
            raise TypeError("G must be a Sage graph")
        if first:
            frm1 = frm
        else:
            frm2 = frm
            to2 = to
        GS.G = G
        GS.directed = 1 if dig else 0
        GS.loops = 1
        GS.use_indicator = 1 if use_indicator_function else 0
        which_G += 1

    if n == 0:
        return {}

    part = PS_from_list(partition)
    ordering = <int *> sig_malloc(n * sizeof(int))
    output = <int *> sig_malloc(n * sizeof(int))
    if part is NULL or ordering is NULL or output is NULL:
        PS_dealloc(part)
        sig_free(ordering)
        sig_free(output)
        raise MemoryError
    for i from 0 <= i < n:
        ordering[i] = to2[ordering2[i]]

    GS1.scratch = <int *> sig_malloc((3*n+1) * sizeof(int))
    GS2.scratch = <int *> sig_malloc((3*n+1) * sizeof(int))
    if GS1.scratch is NULL or GS2.scratch is NULL:
        sig_free(GS1.scratch)
        sig_free(GS2.scratch)
        PS_dealloc(part)
        sig_free(ordering)
        raise MemoryError

    cdef bint isomorphic = double_coset(<void *>GS1, <void *>GS2, part, ordering, n, &all_children_are_equivalent, &refine_by_degree, &compare_graphs, NULL, NULL, output)

    PS_dealloc(part)
    sig_free(ordering)
    sig_free(GS1.scratch)
    sig_free(GS2.scratch)
    if isomorphic:
        output_py = dict([[frm1[i], frm2[output[i]]] for i from 0 <= i < n])
    else:
        output_py = False
    sig_free(output)
    return output_py


def search_tree(G_in, partition, lab=True, dig=False, dict_rep=False, certificate=False,
                    verbosity=0, use_indicator_function=True, sparse=True,
                    base=False, order=False):
    """
    Compute canonical labels and automorphism groups of graphs.

    INPUT:

    - ``G_in`` -- a Sage graph
    - ``partition`` -- list of lists representing a partition of the vertices
    - ``lab`` -- if ``True``, compute and return the canonical label in addition to the
      automorphism group
    - ``dig`` -- set to ``True`` for digraphs and graphs with loops;  if ``True``, does not
      use optimizations based on Lemma 2.25 in [1] that are valid only for
      simple graphs
    - ``dict_rep`` -- if ``True``, return a dictionary with keys the vertices of the
      input graph G_in and values elements of the set the permutation group
      acts on.  (The point is that graphs are arbitrarily labelled, often
      0..n-1, and permutation groups always act on 1..n.  This dictionary
      maps vertex labels (such as 0..n-1) to the domain of the permutations.)
    - ``certificate`` -- if ``True``, return the permutation from G to its canonical label
    - ``verbosity`` -- currently ignored
    - ``use_indicator_function`` -- option to turn off indicator function
      (``True`` is generally faster)
    - ``sparse`` -- whether to use sparse or dense representation of the graph
      (ignored if G is already a CGraph - see sage.graphs.base)
    - ``base`` -- whether to return the first sequence of split vertices (used in
      computing the order of the group)
    - ``order`` -- whether to return the order of the automorphism group

    OUTPUT:

    Depends on the options. If more than one thing is returned, they are in a
    tuple in the following order:

    - list of generators in list-permutation format -- always
    - dict -- if dict_rep
    - graph -- if lab
    - dict -- if certificate
    - list -- if base
    - integer -- if order

    EXAMPLES::

        sage: st = sage.groups.perm_gps.partn_ref.refinement_graphs.search_tree
        sage: from sage.graphs.base.dense_graph import DenseGraph
        sage: from sage.graphs.base.sparse_graph import SparseGraph

    Graphs on zero vertices::

        sage: G = Graph()
        sage: st(G, [[]], order=True)
        ([], Graph on 0 vertices, 1)

    Graphs on one vertex::

        sage: G = Graph(1)
        sage: st(G, [[0]], order=True)
        ([], Graph on 1 vertex, 1)

    Graphs on two vertices::

        sage: G = Graph(2)
        sage: st(G, [[0,1]], order=True)
        ([[1, 0]], Graph on 2 vertices, 2)
        sage: st(G, [[0],[1]], order=True)
        ([], Graph on 2 vertices, 1)
        sage: G.add_edge(0,1)
        sage: st(G, [[0,1]], order=True)
        ([[1, 0]], Graph on 2 vertices, 2)
        sage: st(G, [[0],[1]], order=True)
        ([], Graph on 2 vertices, 1)

    Graphs on three vertices::

        sage: G = Graph(3)
        sage: st(G, [[0,1,2]], order=True)
        ([[0, 2, 1], [1, 0, 2]], Graph on 3 vertices, 6)
        sage: st(G, [[0],[1,2]], order=True)
        ([[0, 2, 1]], Graph on 3 vertices, 2)
        sage: st(G, [[0],[1],[2]], order=True)
        ([], Graph on 3 vertices, 1)
        sage: G.add_edge(0,1)
        sage: st(G, [range(3)], order=True)
        ([[1, 0, 2]], Graph on 3 vertices, 2)
        sage: st(G, [[0],[1,2]], order=True)
        ([], Graph on 3 vertices, 1)
        sage: st(G, [[0,1],[2]], order=True)
        ([[1, 0, 2]], Graph on 3 vertices, 2)

    The Dodecahedron has automorphism group of size 120::

        sage: G = graphs.DodecahedralGraph()
        sage: Pi = [range(20)]
        sage: st(G, Pi, order=True)[2]
        120

    The three-cube has automorphism group of size 48::

        sage: G = graphs.CubeGraph(3)
        sage: G.relabel()
        sage: Pi = [G.vertices(sort=False)]
        sage: st(G, Pi, order=True)[2]
        48

    We obtain the same output using different types of Sage graphs::

        sage: G = graphs.DodecahedralGraph()
        sage: GD = DenseGraph(20)
        sage: GS = SparseGraph(20)
        sage: for i,j,_ in G.edge_iterator():
        ....:  GD.add_arc(i,j); GD.add_arc(j,i)
        ....:  GS.add_arc(i,j); GS.add_arc(j,i)
        sage: Pi = [range(20)]
        sage: a,b = st(G, Pi)
        sage: asp,bsp = st(GS, Pi)
        sage: ade,bde = st(GD, Pi)
        sage: bsg = Graph()
        sage: bdg = Graph()
        sage: for i in range(20):
        ....:    for j in range(20):
        ....:        if bsp.has_arc(i,j):
        ....:            bsg.add_edge(i,j)
        ....:        if bde.has_arc(i,j):
        ....:            bdg.add_edge(i,j)
        sage: a, b.graph6_string()
        ([[0, 19, 3, 2, 6, 5, 4, 17, 18, 11, 10, 9, 13, 12, 16, 15, 14, 7, 8, 1],
          [0, 1, 8, 9, 13, 14, 7, 6, 2, 3, 19, 18, 17, 4, 5, 15, 16, 12, 11, 10],
          [1, 8, 9, 10, 11, 12, 13, 14, 7, 6, 2, 3, 4, 5, 15, 16, 17, 18, 19, 0]],
         'S?[PG__OQ@?_?_?P?CO?_?AE?EC?Ac?@O')
        sage: a == asp
        True
        sage: a == ade
        True
        sage: b == bsg
        True
        sage: b == bdg
        True

    Cubes!::

        sage: C = graphs.CubeGraph(1)
        sage: gens, order = st(C, [C.vertices(sort=False)], lab=False, order=True); order
        2
        sage: C = graphs.CubeGraph(2)
        sage: gens, order = st(C, [C.vertices(sort=False)], lab=False, order=True); order
        8
        sage: C = graphs.CubeGraph(3)
        sage: gens, order = st(C, [C.vertices(sort=False)], lab=False, order=True); order
        48
        sage: C = graphs.CubeGraph(4)
        sage: gens, order = st(C, [C.vertices(sort=False)], lab=False, order=True); order
        384
        sage: C = graphs.CubeGraph(5)
        sage: gens, order = st(C, [C.vertices(sort=False)], lab=False, order=True); order
        3840
        sage: C = graphs.CubeGraph(6)
        sage: gens, order = st(C, [C.vertices(sort=False)], lab=False, order=True); order
        46080

    One can also turn off the indicator function (note: this will take longer)::

        sage: D1 = DiGraph({0:[2],2:[0],1:[1]}, loops=True)
        sage: D2 = DiGraph({1:[2],2:[1],0:[0]}, loops=True)
        sage: a,b = st(D1, [D1.vertices(sort=False)], dig=True, use_indicator_function=False)
        sage: c,d = st(D2, [D2.vertices(sort=False)], dig=True, use_indicator_function=False)
        sage: b==d
        True

    This example is due to Chris Godsil::

        sage: HS = graphs.HoffmanSingletonGraph()
        sage: alqs = [Set(c) for c in (HS.complement()).cliques_maximum()]
        sage: Y = Graph([alqs, lambda s,t: len(s.intersection(t))==0])
        sage: Y0,Y1 = Y.connected_components_subgraphs()
        sage: st(Y0, [Y0.vertices(sort=False)])[1] == st(Y1, [Y1.vertices(sort=False)])[1]
        True
        sage: st(Y0, [Y0.vertices(sort=False)])[1] == st(HS, [HS.vertices(sort=False)])[1]
        True
        sage: st(HS, [HS.vertices(sort=False)])[1] == st(Y1, [Y1.vertices(sort=False)])[1]
        True

    Certain border cases need to be tested as well::

        sage: G = Graph('Fll^G')
        sage: a,b,c = st(G, [range(G.num_verts())], order=True); b
        Graph on 7 vertices
        sage: c
        48
        sage: G = Graph(21)
        sage: st(G, [range(G.num_verts())], order=True)[2] == factorial(21)
        True

        sage: G = Graph('^????????????????????{??N??@w??FaGa?PCO@CP?AGa?_QO?Q@G?CcA??cc????Bo????{????F_')
        sage: perm = {3:15, 15:3}
        sage: H = G.relabel(perm, inplace=False)
        sage: st(G, [range(G.num_verts())])[1] == st(H, [range(H.num_verts())])[1]
        True

        sage: st(Graph(':Dkw'), [range(5)], lab=False, dig=True)
        [[4, 1, 2, 3, 0], [0, 2, 1, 3, 4]]
    """
    cdef CGraph G
    cdef int i, j, n
    cdef Integer I
    cdef bint loops
    cdef aut_gp_and_can_lab *output
    cdef PartitionStack *part
    from sage.graphs.graph import Graph
    from sage.graphs.generic_graph import GenericGraph
    from copy import copy
    if isinstance(G_in, GenericGraph):
        loops = G_in.has_loops()
        n = G_in.num_verts()
        if G_in.vertices(sort=False) != list(range(n)):
            G_in = copy(G_in)
            to = G_in.relabel(return_map=True)
            frm = {}
            for v in to:
                frm[to[v]] = v
            partition = [[to[v] for v in cell] for cell in partition]
        else:
            to = dict(enumerate(range(n)))
            frm = to
        if sparse:
            G = SparseGraph(n)
        else:
            G = DenseGraph(n)
        if G_in.is_directed():
            for i, j in G_in.edge_iterator(labels=False):
                G.add_arc(i, j)
        else:
            for i, j in G_in.edge_iterator(labels=False):
                G.add_arc(i, j)
                G.add_arc(j, i)
    elif isinstance(G_in, CGraph):
        G = <CGraph> G_in
        n = G.num_verts
        loops = 0
        for i in range(n):
            if G.has_arc_unsafe(i, i):
                loops = 1
        to = {}
        for a in G.verts():
            to[a] = a
        frm = to
    else:
        raise TypeError("G must be a Sage graph")

    cdef GraphStruct GS = GraphStruct()
    GS.G = G
    GS.directed = 1 if dig else 0
    GS.loops = loops
    GS.use_indicator = 1 if use_indicator_function else 0

    if n == 0:
        return_tuple = [[]]
        if dict_rep:
            return_tuple.append({})
        if lab:
            if isinstance(G_in, GenericGraph):
                G_C = copy(G_in)
            else:
                if isinstance(G, SparseGraph):
                    G_C = SparseGraph(n)
                else:
                    G_C = DenseGraph(n)
            return_tuple.append(G_C)
        if certificate:
            return_tuple.append({})
        if base:
            return_tuple.append([])
        if order:
            return_tuple.append(Integer(1))
        if len(return_tuple) == 1:
            return return_tuple[0]
        else:
            return tuple(return_tuple)

    GS.scratch = <int *> sig_malloc( (3*G.num_verts + 1) * sizeof(int) )
    part = PS_from_list(partition)
    if GS.scratch is NULL or part is NULL:
        PS_dealloc(part)
        sig_free(GS.scratch)
        raise MemoryError

    output = get_aut_gp_and_can_lab(<void *>GS, part, G.num_verts, &all_children_are_equivalent, &refine_by_degree, &compare_graphs, lab, NULL, NULL, NULL)
    sig_free( GS.scratch )
    # prepare output
    list_of_gens = []
    for i from 0 <= i < output.num_gens:
        list_of_gens.append([output.generators[j+i*G.num_verts] for j from 0 <= j < <int>G.num_verts])
    return_tuple = [list_of_gens]
    if dict_rep:
        ddd = {}
        for v in frm:
            ddd[frm[v]] = v if v != 0 else n
        return_tuple.append(ddd)
    if lab:
        if isinstance(G_in, GenericGraph):
            G_C = copy(G_in)
            G_C.relabel([output.relabeling[i] for i from 0 <= i < n])
        else:
            if isinstance(G, SparseGraph):
                G_C = SparseGraph(n)
            else:
                G_C = DenseGraph(n)
            for i from 0 <= i < n:
                for j in G.out_neighbors(i):
                    G_C.add_arc(output.relabeling[i],output.relabeling[j])
        return_tuple.append(G_C)
    if certificate:
        dd = {}
        for i from 0 <= i < <int>G.num_verts:
            dd[frm[i]] = output.relabeling[i]
        return_tuple.append(dd)
    if base:
        return_tuple.append([output.group.base_orbits[i][0] for i from 0 <= i < output.group.base_size])
    if order:
        I = Integer()
        SC_order(output.group, 0, I.value)
        return_tuple.append(I)
    PS_dealloc(part)
    deallocate_agcl_output(output)
    if len(return_tuple) == 1:
        return return_tuple[0]
    else:
        return tuple(return_tuple)


cdef int refine_by_degree(PartitionStack *PS, void *S, int *cells_to_refine_by, int ctrb_len) noexcept:
    r"""
    Refine the input partition by checking degrees of vertices to the given
    cells.

    INPUT:

    - ``PS`` -- a partition stack, whose finest partition is the partition to be
      refined
    - ``S`` -- a graph struct object, which contains scratch space, the graph in
      question, and some flags
    - ``cells_to_refine_by`` -- list of pointers to cells to check degrees against
      in refining the other cells (updated in place). Must be allocated to
      length at least the degree of PS, since the array may grow
    - ``ctrb_len`` -- how many cells in cells_to_refine_by

    OUTPUT:

    An integer invariant under the orbits of $S_n$.  That is, if $\gamma$ is a
    permutation of the vertices, then
    $$ I(G, PS, cells_to_refine_by) = I( \gamma(G), \gamma(PS), \gamma(cells_to_refine_by) ) .$$
    """
    cdef GraphStruct GS = <GraphStruct> S
    cdef CGraph G = GS.G
    cdef int current_cell_against = 0
    cdef int current_cell, i, r
    cdef int first_largest_subcell
    cdef int invariant = 1
    cdef int max_degree
    cdef int *degrees = GS.scratch # length 3n+1
    cdef bint necessary_to_split_cell
    cdef int against_index
    if <int>G.num_verts != PS.degree and PS.depth == 0:
        # should be less verts, then, so place the "nonverts" in separate cell at the end
        current_cell = 0
        while current_cell < PS.degree:
            i = current_cell
            r = 0
            while True:
                if G.has_vertex(PS.entries[i]):
                    degrees[i-current_cell] = 0
                else:
                    r = 1
                    degrees[i-current_cell] = 1
                i += 1
                if PS.levels[i-1] <= PS.depth:
                    break
            if r != 0:
                sort_by_function(PS, current_cell, degrees)
            current_cell = i
    while not PS_is_discrete(PS) and current_cell_against < ctrb_len:
        invariant += 1
        current_cell = 0
        while current_cell < PS.degree:
            invariant += 50
            i = current_cell
            necessary_to_split_cell = 0
            max_degree = 0
            while True:
                degrees[i-current_cell] = degree(PS, G, i, cells_to_refine_by[current_cell_against], 0)
                if degrees[i-current_cell] != degrees[0]:
                    necessary_to_split_cell = 1
                if degrees[i-current_cell] > max_degree:
                    max_degree = degrees[i-current_cell]
                i += 1
                if PS.levels[i-1] <= PS.depth:
                    break
            # now, i points to the next cell (before refinement)
            if necessary_to_split_cell:
                invariant += 10
                first_largest_subcell = sort_by_function(PS, current_cell, degrees)
                invariant += first_largest_subcell + max_degree
                against_index = current_cell_against
                while against_index < ctrb_len:
                    if cells_to_refine_by[against_index] == current_cell:
                        cells_to_refine_by[against_index] = first_largest_subcell
                        break
                    against_index += 1
                r = current_cell
                while True:
                    if r == current_cell or PS.levels[r-1] == PS.depth:
                        if r != first_largest_subcell:
                            cells_to_refine_by[ctrb_len] = r
                            ctrb_len += 1
                    r += 1
                    if r >= i:
                        break
                invariant += (i - current_cell)
            current_cell = i
        if GS.directed:
            # if we are looking at a digraph, also compute
            # the reverse degrees and sort by them
            current_cell = 0
            while current_cell < PS.degree: # current_cell is still a valid cell
                invariant += 20
                i = current_cell
                necessary_to_split_cell = 0
                max_degree = 0
                while True:
                    degrees[i-current_cell] = degree(PS, G, i, cells_to_refine_by[current_cell_against], 1)
                    if degrees[i-current_cell] != degrees[0]:
                        necessary_to_split_cell = 1
                    if degrees[i-current_cell] > max_degree:
                        max_degree = degrees[i-current_cell]
                    i += 1
                    if PS.levels[i-1] <= PS.depth:
                        break
                # now, i points to the next cell (before refinement)
                if necessary_to_split_cell:
                    invariant += 7
                    first_largest_subcell = sort_by_function(PS, current_cell, degrees)
                    invariant += first_largest_subcell + max_degree
                    against_index = current_cell_against
                    while against_index < ctrb_len:
                        if cells_to_refine_by[against_index] == current_cell:
                            cells_to_refine_by[against_index] = first_largest_subcell
                            break
                        against_index += 1
                    against_index = ctrb_len
                    r = current_cell
                    while True:
                        if r == current_cell or PS.levels[r-1] == PS.depth:
                            if r != first_largest_subcell:
                                cells_to_refine_by[against_index] = r
                                against_index += 1
                                ctrb_len += 1
                        r += 1
                        if r >= i:
                            break
                    invariant += (i - current_cell)
                current_cell = i
        current_cell_against += 1
    if GS.use_indicator:
        return invariant
    else:
        return 0

cdef int compare_graphs(int *gamma_1, int *gamma_2, void *S1, void *S2, int degree) noexcept:
    r"""
    Compare gamma_1(S1) and gamma_2(S2).

    Return -1 if gamma_1(S1) < gamma_2(S2), 0 if gamma_1(S1) ==
    gamma_2(S2), 1 if gamma_1(S1) > gamma_2(S2).

    INPUT:

    - ``gamma_1``, ``gamma_2`` -- list permutations (inverse)
    - ``S1``, ``S2`` -- graph struct objects
    """
    cdef size_t i, j
    cdef GraphStruct GS1 = <GraphStruct> S1
    cdef GraphStruct GS2 = <GraphStruct> S2
    cdef CGraph G1 = GS1.G
    cdef CGraph G2 = GS2.G
    if G1.active_vertices.size != G2.active_vertices.size or \
       not bitset_cmp(G1.active_vertices, G2.active_vertices):
        for i from 0 <= i < <size_t>degree:
            if G1.has_vertex(gamma_1[i]) != G2.has_vertex(gamma_2[i]):
                return G1.has_vertex(gamma_1[i]) - G2.has_vertex(gamma_2[i])
    for i from 0 <= i < G1.num_verts:
        for j from 0 <= j < G1.num_verts:
            if G1.has_arc_unsafe(gamma_1[i], gamma_1[j]):
                if not G2.has_arc_unsafe(gamma_2[i], gamma_2[j]):
                    return 1
            elif G2.has_arc_unsafe(gamma_2[i], gamma_2[j]):
                return -1
    return 0

cdef bint all_children_are_equivalent(PartitionStack *PS, void *S) noexcept:
    """
    Return ``True`` if every refinement of the current partition results in the
    same structure.

    .. WARNING::

        Converse does not hold in general!  See Lemma 2.25 of [1] for details.

    INPUT:

    - ``PS`` -- the partition stack to be checked
    - ``S`` -- a graph struct object
    """
    cdef GraphStruct GS = <GraphStruct> S
    if GS.directed or GS.loops:
        return 0
    cdef int i, n = PS.degree
    cdef bint in_cell = 0
    cdef int nontrivial_cells = 0
    cdef int total_cells = PS_num_cells(PS)
    if n <= total_cells + 4:
        return 1
    for i from 0 <= i < n-1:
        if PS.levels[i] <= PS.depth:
            if in_cell:
                nontrivial_cells += 1
            in_cell = 0
        else:
            in_cell = 1
    if in_cell:
        nontrivial_cells += 1
    if n == total_cells + nontrivial_cells:
        return 1
    if n == total_cells + nontrivial_cells + 1:
        return 1
    return 0

cdef inline int degree(PartitionStack *PS, CGraph G, int entry, int cell_index, bint reverse) noexcept:
    """
    Return the number of edges from the vertex corresponding to entry to
    vertices in the cell corresponding to cell_index.

    INPUT:

    - ``PS`` -- the partition stack to be checked
    - ``S`` -- a graph struct object
    - ``entry`` -- the position of the vertex in question in the entries of PS
    - ``cell_index`` -- the starting position of the cell in question in the entries
      of PS
    - ``reverse`` -- whether to check for arcs in the other direction
    """
    cdef int num_arcs = 0
    entry = PS.entries[entry]
    if not reverse:
        while True:
            if G.has_arc_unsafe(PS.entries[cell_index], entry):
                num_arcs += 1
            if PS.levels[cell_index] > PS.depth:
                cell_index += 1
            else:
                break
    else:
        while True:
            if G.has_arc_unsafe(entry, PS.entries[cell_index]):
                num_arcs += 1
            if PS.levels[cell_index] > PS.depth:
                cell_index += 1
            else:
                break
    return num_arcs


def all_labeled_graphs(n):
    """
    Return all labeled graphs on n vertices {0,1,...,n-1}.

    Used in classifying isomorphism types (naive approach), and more
    importantly in benchmarking the search algorithm.

    EXAMPLES::

        sage: from sage.groups.perm_gps.partn_ref.refinement_graphs import all_labeled_graphs
        sage: st = sage.groups.perm_gps.partn_ref.refinement_graphs.search_tree
        sage: Glist = {}
        sage: Giso  = {}
        sage: for n in [1..5]:  # long time (4s on sage.math, 2011)
        ....:  Glist[n] = all_labeled_graphs(n)
        ....:  Giso[n] = []
        ....:  for g in Glist[n]:
        ....:      a, b = st(g, [range(n)])
        ....:      inn = False
        ....:      for gi in Giso[n]:
        ....:          if b == gi:
        ....:              inn = True
        ....:      if not inn:
        ....:          Giso[n].append(b)
        sage: for n in Giso:  # long time
        ....:    print("{} {}".format(n, len(Giso[n])))
        1 1
        2 2
        3 4
        4 11
        5 34
    """
    from sage.graphs.graph import Graph
    TE = []
    for i in range(n):
        for j in range(i):
            TE.append((i, j))
    m = len(TE)
    Glist= []
    for i in range(2**m):
        G = Graph(n)
        b = Integer(i).binary()
        b = '0'*(m-len(b)) + b
        for i in range(m):
            if int(b[i]):
                G.add_edge(TE[i])
        Glist.append(G)
    return Glist


def random_tests(num=10, n_max=60, perms_per_graph=5):
    """
    Test to make sure that C(gamma(G)) == C(G) for random permutations gamma
    and random graphs G, and that isomorphic returns an isomorphism.

    INPUT:

    - ``num`` -- run tests for this many graphs
    - ``n_max`` -- test graphs with at most this many vertices
    - ``perms_per_graph`` -- test each graph with this many random permutations

    DISCUSSION:

    This code generates num random graphs G on at most n_max vertices.  The
    density of edges is chosen randomly between 0 and 1.

    For each graph G generated, we uniformly generate perms_per_graph random
    permutations and verify that the canonical labels of G and the image of G
    under the generated permutation are equal, and that the isomorphic function
    returns an isomorphism.

    TESTS::

        sage: import sage.groups.perm_gps.partn_ref.refinement_graphs
        sage: sage.groups.perm_gps.partn_ref.refinement_graphs.random_tests()  # long time
        All passed: 200 random tests on 20 graphs.
    """
    from sage.misc.prandom import random, randint
    from sage.graphs.graph_generators import GraphGenerators
    from sage.graphs.digraph_generators import DiGraphGenerators
    from sage.combinat.permutation import Permutations
    from copy import copy
    cdef int i, j, num_tests = 0, num_graphs = 0
    GG = GraphGenerators()
    DGG = DiGraphGenerators()
    for mmm in range(num):
        p = random()
        n = randint(1, n_max)
        S = Permutations(n)

        G = GG.RandomGNP(n, p)
        H = copy(G)
        for i from 0 <= i < perms_per_graph:
            G = copy(H)
            G1 = search_tree(G, [G.vertices(sort=False)])[1]
            perm = list(S.random_element())
            perm = [perm[j]-1 for j from 0 <= j < n]
            G.relabel(perm)
            G2 = search_tree(G, [G.vertices(sort=False)])[1]
            if G1 != G2:
                print("search_tree FAILURE: graph6-")
                print(H.graph6_string())
                print(perm)
                return
            isom = isomorphic(G, H, [list(range(n))], list(range(n)), 0, 1)
            if not isom or G.relabel(isom, inplace=False) != H:
                print("isom FAILURE: graph6-")
                print(H.graph6_string())
                print(perm)
                return

        D = DGG.RandomDirectedGNP(n, p)
        D.allow_loops(True)
        for i from 0 <= i < n:
            if random() < p:
                D.add_edge(i,i)
        E = copy(D)
        for i from 0 <= i < perms_per_graph:
            D = copy(E)
            D1 = search_tree(D, [D.vertices(sort=False)], dig=True)[1]
            perm = list(S.random_element())
            perm = [perm[j]-1 for j from 0 <= j < n]
            D.relabel(perm)
            D2 = search_tree(D, [D.vertices(sort=False)], dig=True)[1]
            if D1 != D2:
                print("search_tree FAILURE: dig6-")
                print(E.dig6_string())
                print(perm)
                return
            isom = isomorphic(D, E, [list(range(n))], list(range(n)), 1, 1)
            if not isom or D.relabel(isom, inplace=False) != E:
                print("isom FAILURE: dig6-")
                print(E.dig6_string())
                print(perm)
                print(isom)
                return
        num_tests += 4*perms_per_graph
        num_graphs += 2
    print("All passed: %d random tests on %d graphs." % (num_tests, num_graphs))


def orbit_partition(gamma, list_perm=False):
    r"""
    Assuming that G is a graph on vertices 0,1,...,n-1, and gamma is an
    element of SymmetricGroup(n), returns the partition of the vertex
    set determined by the orbits of gamma, considered as action on the
    set 1,2,...,n where we take 0 = n. In other words, returns the
    partition determined by a cyclic representation of gamma.

    INPUT:

    - ``list_perm`` -- if ``True``, assumes
      ``gamma`` is a list representing the map
      `i \mapsto ``gamma``[i]`

    EXAMPLES::

        sage: # needs sage.groups
        sage: from sage.groups.perm_gps.partn_ref.refinement_graphs import orbit_partition
        sage: G = graphs.PetersenGraph()
        sage: S = SymmetricGroup(10)
        sage: gamma = S('(10,1,2,3,4)(5,6,7)(8,9)')
        sage: orbit_partition(gamma)
        [[1, 2, 3, 4, 0], [5, 6, 7], [8, 9]]
        sage: gamma = S('(10,5)(1,6)(2,7)(3,8)(4,9)')
        sage: orbit_partition(gamma)
        [[1, 6], [2, 7], [3, 8], [4, 9], [5, 0]]
    """
    if list_perm:
        n = len(gamma)
        seen = [1] + [0]*(n-1)
        i = 0
        p = 0
        partition = [[0]]
        while sum(seen) < n:
            if gamma[i] != partition[p][0]:
                partition[p].append(gamma[i])
                i = gamma[i]
                seen[i] = 1
            else:
                for j in range(n):
                    if seen[j] == 0:
                        i = j
                        break
                partition.append([i])
                p += 1
                seen[i] = 1
        return partition
    else:
        n = len(gamma.domain())
        l = []
        for i in range(1,n+1):
            orb = gamma.orbit(i)
            if orb not in l:
                l.append(orb)
        for i in l:
            for j in range(len(i)):
                if i[j] == n:
                    i[j] = 0
        return l


def coarsest_equitable_refinement(CGraph G, list partition, bint directed):
    """
    Return the coarsest equitable refinement of ``partition`` for ``G``.

    This is a helper function for the graph function of the same name.

    DOCTEST (More thorough testing in ``sage/graphs/graph.py``)::

        sage: from sage.groups.perm_gps.partn_ref.refinement_graphs import coarsest_equitable_refinement
        sage: from sage.graphs.base.sparse_graph import SparseGraph
        sage: coarsest_equitable_refinement(SparseGraph(7), [[0], [1,2,3,4], [5,6]], 0)
        [[0], [1, 2, 3, 4], [5, 6]]
    """
    cdef int i, j = 0, k = 0, n = G.num_verts

    # set up partition stack and graph struct
    cdef PartitionStack *nu = PS_new(n, 0)
    for cell in partition:
        for i in cell:
            nu.entries[j] = i
            nu.levels[j] = n
            j += 1
        nu.levels[j-1] = 0
        PS_move_min_to_front(nu, k, j-1)
        k = j

    cdef GraphStruct GS = GraphStruct()
    GS.G = G
    GS.scratch = <int *> sig_malloc((3*n+1) * sizeof(int))
    if GS.scratch is NULL:
        PS_dealloc(nu)
        raise MemoryError
    GS.directed = directed
    GS.use_indicator = 0

    # set up cells to refine by
    cdef int num_cells = len(partition)
    cdef int *alpha = <int *>sig_malloc(n * sizeof(int))
    if alpha is NULL:
        PS_dealloc(nu)
        sig_free(GS.scratch)
        raise MemoryError
    j = 0
    for i from 0 <= i < num_cells:
        alpha[i] = j
        j += len(partition[i])

    # refine, and get the result
    refine_by_degree(nu, <void *>GS, alpha, num_cells)

    eq_part = []
    cell = []
    for i from 0 <= i < n:
        cell.append(nu.entries[i])
        if nu.levels[i] <= 0:
            eq_part.append(cell)
            cell = []

    PS_dealloc(nu)
    sig_free(GS.scratch)
    sig_free(alpha)

    return eq_part


def get_orbits(list gens, int n):
    """
    Compute orbits given a list of generators of a permutation group, in list
    format.

    This is a helper function for automorphism groups of graphs.

    DOCTEST (More thorough testing in ``sage/graphs/graph.py``)::

        sage: from sage.groups.perm_gps.partn_ref.refinement_graphs import get_orbits
        sage: get_orbits([[1,2,3,0,4,5], [0,1,2,3,5,4]], 6)
        [[0, 1, 2, 3], [4, 5]]
    """
    cdef int i, j
    if len(gens) == 0:
        return [[i] for i from 0 <= i < n]

    cdef OrbitPartition *OP = OP_new(n)
    cdef int *perm_ints = <int *> sig_malloc(n * sizeof(int))
    if perm_ints is NULL:
        OP_dealloc(OP)
        raise MemoryError

    for gen in gens:
        for i from 0 <= i < n:
            perm_ints[i] = gen[i]
        OP_merge_list_perm(OP, perm_ints)

    orbit_dict = {}
    for i from 0 <= i < n:
        j = OP_find(OP, i)
        if j in orbit_dict:
            orbit_dict[j].append(i)
        else:
            orbit_dict[j] = [i]

    OP_dealloc(OP)
    sig_free(perm_ints)

    return list(orbit_dict.itervalues())


# Canonical augmentation
from cpython.ref cimport *

# Dense graphs: adding edges

# This implements an augmentation scheme as follows:
# * Seed objects are graphs with n vertices and no edges.
# * Augmentations consist of adding a single edge, or a loop.

cdef void *dg_edge_gen_next(void *data, int *degree, bint *mem_err) noexcept:
    r"""
    The ``next`` function in an edge iterator. The iterator generates unique
    representatives under the action of the automorphism group of the parent
    graph on edges not in the graph, which are to considered for adding to the
    graph.
    """
    cdef dg_edge_gen_data *degd = <dg_edge_gen_data *> data
    cdef GraphStruct graph = <GraphStruct> degd.graph
    cdef subset *edge_candidate
    cdef int u, v, reject
    cdef bint mem_err_sub = 0
    if mem_err[0]:
        (<canonical_generator_data *> degd.edge_iterator.data).mem_err = 1
    while True:
        edge_candidate = <subset *> degd.edge_iterator.next(degd.edge_iterator.data, NULL, &mem_err_sub)
        if edge_candidate is NULL:
            break
        reject = 0
        if bitset_len(&edge_candidate.bits) < (1 if graph.loops else 2):
            reject = 1
        else:
            u = bitset_first(&edge_candidate.bits)
            v = bitset_next(&edge_candidate.bits, u+1)
            if v == -1:
                v = u
            if graph.G.has_arc_unsafe(u, v):
                reject = 1
        if not reject:
            break
    if mem_err_sub:
        mem_err[0] = 1
    return edge_candidate

cdef void *allocate_degd(int degree) noexcept:
    r"""
    Allocate the data part of the iterator over edges to add to the graph.
    """
    cdef dg_edge_gen_data *degd = <dg_edge_gen_data *> sig_malloc(sizeof(dg_edge_gen_data))
    cdef iterator *edge_iterator = allocate_subset_gen(degree, 2)
    if degd is NULL or edge_iterator is NULL:
        sig_free(degd)
        free_subset_gen(edge_iterator)
        return NULL
    edge_iterator = setup_set_gen(edge_iterator, degree, 2)
    if edge_iterator is NULL:
        sig_free(degd)
        return NULL
    degd.edge_iterator = edge_iterator
    return degd

cdef void deallocate_degd(void *data) noexcept:
    r"""
    Deallocate the data part of the iterator over edges to add to the graph.
    """
    cdef dg_edge_gen_data *degd = <dg_edge_gen_data *> data
    free_subset_gen(degd.edge_iterator)
    sig_free(degd)

cdef int gen_children_dg_edge(void *S, aut_gp_and_can_lab *group, iterator *it) noexcept:
    r"""
    Setup an iterator over edges to be added.
    """
    cdef GraphStruct GS = <GraphStruct> S
    cdef int n = GS.G.num_verts
    (<dg_edge_gen_data *> it.data).graph = <void *> GS
    cdef iterator *edge_iterator = setup_set_gen((<dg_edge_gen_data *> it.data).edge_iterator, n, 2)
    if edge_iterator is not NULL:
        start_canonical_generator(group.group, NULL, n, edge_iterator)
    return (edge_iterator is NULL)

cdef void *apply_dg_edge_aug(void *parent, void *aug, void *child, int *degree, bint *mem_err) noexcept:
    r"""
    Apply the augmentation to ``parent`` storing the result in ``child``. Here
    ``aug`` represents an edge to be added.
    """
    cdef GraphStruct GS_child = <GraphStruct> child, GS_par = <GraphStruct> parent
    cdef DenseGraph DG = <DenseGraph> GS_child.G, DG_par = <DenseGraph> GS_par.G
    cdef subset *edge = <subset *> aug
    cdef int u, v

    # copy DG_par edges to DG
    copy_dense_graph(DG, DG_par)

    # add the edge
    u = bitset_first(&edge.bits)
    v = bitset_next(&edge.bits, u+1)
    if v == -1:
        DG.add_arc_unsafe(u, u)
    else:
        DG.add_arc_unsafe(u, v)
        DG.add_arc_unsafe(v, u)

    degree[0] = DG.num_verts
    return <void *> GS_child

cdef void *allocate_dg_edge(int n, bint loops) noexcept:
    r"""
    Allocate an object for this augmentation scheme.
    """
    cdef GraphStruct GS
    cdef DenseGraph G
    cdef int *scratch
    try:
        GS = GraphStruct()
        G  = DenseGraph(n)
        scratch = <int *> sig_malloc((3*n+1) * sizeof(int))
        if scratch is NULL:
            raise MemoryError
    except MemoryError:
        return NULL
    Py_INCREF(GS)
    Py_INCREF(G)
    GS.G = G
    GS.directed = 0
    GS.loops = loops
    GS.use_indicator = 1
    GS.scratch = scratch
    return <void *> GS

cdef void free_dg_edge(void *child) noexcept:
    r"""
    Deallocate an object for this augmentation scheme.
    """
    cdef GraphStruct GS = <GraphStruct> child
    sig_free(GS.scratch)
    Py_DECREF(GS.G)
    Py_DECREF(GS)

cdef void *canonical_dg_edge_parent(void *child, void *parent, int *permutation, int *degree, bint *mem_err) noexcept:
    r"""
    Apply ``permutation`` to ``child``, determine an arbitrary parent by
    deleting the lexicographically largest edge, apply the inverse of
    ``permutation`` to the result and store the result in ``parent``.
    """
    cdef GraphStruct GS_par = <GraphStruct> parent, GS = <GraphStruct> child
    cdef DenseGraph DG_par = <DenseGraph> GS_par.G, DG = <DenseGraph> GS.G
    cdef int u, v, n = DG.num_verts
    cdef int *scratch = GS_par.scratch

    # copy DG edges to DG_par
    copy_dense_graph(DG_par, DG)

    # remove the right edge
    for u from 0 <= u < n:
        scratch[permutation[u]] = u
    for u from n > u >= 0:
        if DG.in_degrees[scratch[u]] != 0:
            break
    for v from u >= v >= 0:
        if DG.has_arc_unsafe(scratch[u], scratch[v]):
            break
    DG_par.del_arc_unsafe(scratch[u], scratch[v])
    if u != v:
        DG_par.del_arc_unsafe(scratch[v], scratch[u])

    degree[0] = n
    return <void *> GS_par

cdef iterator *allocate_dg_edge_gen(int degree, int depth, bint loops) noexcept:
    r"""
    Allocate the iterator for generating graphs.
    """
    cdef iterator *dg_edge_gen = <iterator *> sig_malloc(sizeof(iterator))
    cdef canonical_generator_data *cgd = allocate_cgd(depth, degree)
    if dg_edge_gen is NULL or cgd is NULL:
        sig_free(dg_edge_gen)
        deallocate_cgd(cgd)
        return NULL
    cdef int i, j
    for i from 0 <= i < depth:
        cgd.object_stack[i] = allocate_dg_edge(degree, loops)
        cgd.parent_stack[i] = allocate_dg_edge(degree, loops)
        cgd.iterator_stack[i].data = allocate_degd(degree)
        cgd.iterator_stack[i].next = &dg_edge_gen_next
        if cgd.iterator_stack[i].data is NULL or \
           cgd.object_stack[i]        is NULL or \
           cgd.parent_stack[i]        is NULL:
            for j from 0 <= j <= i:
                deallocate_degd(cgd.iterator_stack[j].data)
                free_dg_edge(cgd.object_stack[j])
                free_dg_edge(cgd.parent_stack[j])
            sig_free(dg_edge_gen)
            deallocate_cgd(cgd)
            return NULL
    dg_edge_gen.data = <void *> cgd
    dg_edge_gen.next = canonical_generator_next
    return dg_edge_gen

cdef void free_dg_edge_gen(iterator *dg_edge_gen) noexcept:
    r"""
    Deallocate the iterator for generating graphs.
    """
    cdef canonical_generator_data *cgd = <canonical_generator_data *> dg_edge_gen.data
    deallocate_cgd(cgd)
    sig_free(dg_edge_gen)


def generate_dense_graphs_edge_addition(int n, bint loops, G=None, depth=None,
                                        bint construct=False,
                                        bint indicate_mem_err=True):
    r"""

    EXAMPLES::

        sage: from sage.groups.perm_gps.partn_ref.refinement_graphs import generate_dense_graphs_edge_addition

    ::

        sage: for n in [0..6]:
        ....:     print(generate_dense_graphs_edge_addition(n,1))
        1
        2
        6
        20
        90
        544
        5096

    ::

        sage: for n in [0..7]:
        ....:     print(generate_dense_graphs_edge_addition(n,0))
        1
        1
        2
        4
        11
        34
        156
        1044
        sage: generate_dense_graphs_edge_addition(8,0)  # long time (about 14 seconds at 2.4 GHz)
        12346
    """
    from sage.graphs.graph import Graph
    cdef iterator *graph_iterator
    cdef DenseGraph DG, ODG
    cdef GraphStruct GS
    if n < 0:
        return [] if construct else Integer(0)
    if n == 0:
        return [Graph(0, implementation='c_graph', sparse=False, loops=loops)] if construct else Integer(1)
    if n == 1:
        if not loops:
            return [Graph(1, implementation='c_graph', sparse=False, loops=loops)] if construct else Integer(1)
        else:
            if construct:
                G = Graph(1, implementation='c_graph', sparse=False, loops=loops)
                (<CGraph>G._backend._cg).add_arc_unsafe(0,0)
                return [G, Graph(1, implementation='c_graph', sparse=False, loops=loops)]
            else:
                return Integer(2)

    if depth is None:
        depth = n*n

    graph_iterator = allocate_dg_edge_gen(n, depth, loops)
    if graph_iterator is NULL:
        raise MemoryError

    GS = (<GraphStruct> (<canonical_generator_data *> graph_iterator.data).object_stack[0])
    if G is not None:
        DG = GS.G
        for u,v in G.edges(sort=False, labels=False):
            DG.add_arc(u,v)
            if u != v:
                DG.add_arc(v,u)

    graph_iterator = setup_canonical_generator(n,
        all_children_are_equivalent,
        refine_by_degree,
        compare_graphs,
        gen_children_dg_edge,
        apply_dg_edge_aug,
        free_dg_edge,
        deallocate_degd,
        free_subset,
        canonical_dg_edge_parent,
        depth, 0, graph_iterator)

    start_canonical_generator(NULL, <void *> GS, n, graph_iterator)

    cdef list out_list
    cdef void *thing
    cdef Integer number
    cdef bint mem_err = 0
    if construct:
        out_list = []
    else:
        number = Integer(0)
    if construct:
        while True:
            thing = graph_iterator.next(graph_iterator.data, NULL, &mem_err)
            if thing is NULL:
                break
            ODG = (<GraphStruct>thing).G
            G = Graph(0, implementation='c_graph', sparse=False)
            DG = DenseGraph(ODG.active_vertices.size, extra_vertices=0)
            copy_dense_graph(DG, ODG)
            G._backend._cg = DG
            out_list.append(G)
    else:
        while True:
            thing = graph_iterator.next(graph_iterator.data, NULL, &mem_err)
            if thing is NULL:
                break
            number += 1

    free_dg_edge_gen(graph_iterator)
    if mem_err:
        if indicate_mem_err:
            raise MemoryError
        else:
            out_list.append(MemoryError())
    if construct:
        return out_list
    else:
        return number


# Dense graphs: adding vertices

# This implements an augmentation scheme as follows:
# * Seed objects are graphs with one vertex and no edges.
# * Augmentations consist of adding a single vertex connected to some subset of
#   the previous vertices.

cdef int gen_children_dg_vert(void *S, aut_gp_and_can_lab *group, iterator *it) noexcept:
    r"""
    Setup an iterator over subsets to join a new vertex to.
    """
    cdef GraphStruct GS = <GraphStruct> S
    cdef int n = GS.G.num_verts
    cdef iterator *subset_iterator = setup_set_gen(it, n, n)
    if subset_iterator is not NULL:
        start_canonical_generator(group.group, NULL, n, subset_iterator)
    return (subset_iterator is NULL)

cdef void *apply_dg_vert_aug(void *parent, void *aug, void *child, int *degree, bint *mem_err) noexcept:
    r"""
    Apply the augmentation to ``parent`` storing the result in ``child``. Here
    ``aug`` represents a subset to join to a new vertex.
    """
    cdef GraphStruct GS_child = <GraphStruct> child, GS_par = <GraphStruct> parent
    cdef DenseGraph DG = <DenseGraph> GS_child.G, DG_par = <DenseGraph> GS_par.G
    cdef subset *set1 = <subset *> aug
    cdef int u, n = DG_par.num_verts

    # copy DG_par edges to DG
    copy_dense_graph(DG, DG_par)
    DG.add_vertex_unsafe(n)

    # add the edges
    u = bitset_first(&set1.bits)
    while u != -1:
        DG.add_arc_unsafe(u, n)
        DG.add_arc_unsafe(n, u)
        u = bitset_next(&set1.bits, u+1)

    degree[0] = n+1
    return <void *> GS_child

cdef void *allocate_dg_vert(int n, int depth) noexcept:
    r"""
    Allocate an object for this augmentation scheme.
    """
    cdef GraphStruct GS
    cdef DenseGraph G
    cdef int *scratch
    try:
        GS = GraphStruct()
        G  = DenseGraph(0, extra_vertices=depth)
        bitset_set_first_n(G.active_vertices, n)
        G.num_verts = n
        scratch = <int *> sig_malloc((3*depth+1) * sizeof(int))
        if scratch is NULL:
            raise MemoryError
    except MemoryError:
        return NULL
    Py_INCREF(GS)
    Py_INCREF(G)
    GS.G = G
    GS.directed = 0
    GS.loops = 0
    GS.use_indicator = 1
    GS.scratch = scratch
    return <void *> GS

cdef void free_dg_vert(void *child) noexcept:
    r"""
    Deallocate an object for this augmentation scheme.
    """
    cdef GraphStruct GS = <GraphStruct> child
    sig_free(GS.scratch)
    Py_DECREF(GS.G)
    Py_DECREF(GS)

cdef void *canonical_dg_vert_parent(void *child, void *parent, int *permutation, int *degree, bint *mem_err) noexcept:
    r"""
    Apply ``permutation`` to ``child``, determines an arbitrary parent by
    deleting the lexicographically largest vertex, apply the inverse of
    ``permutation`` to the result and store the result in ``parent``.
    """
    cdef GraphStruct GS_par = <GraphStruct> parent, GS = <GraphStruct> child
    cdef DenseGraph DG_par = <DenseGraph> GS_par.G, DG = <DenseGraph> GS.G
    cdef int u, n = DG_par.num_verts
    cdef int *scratch = GS.scratch

    # copy DG edges to DG_par
    copy_dense_graph(DG_par, DG)

    # remove the right vertex
    for u from 0 <= u <= n:
        scratch[permutation[u]] = u
    DG_par.del_vertex_unsafe(scratch[n])

    degree[0] = n
    return <void *> GS_par

cdef iterator *allocate_dg_vert_gen(int degree, int depth) noexcept:
    r"""
    Allocate the iterator for generating graphs.
    """
    cdef iterator *dg_vert_gen = <iterator *> sig_malloc(sizeof(iterator))
    cdef canonical_generator_data *cgd = allocate_cgd(depth, degree)
    cdef canonical_generator_data *cgd2
    if dg_vert_gen is NULL or cgd is NULL:
        sig_free(dg_vert_gen)
        deallocate_cgd(cgd)
        return NULL
    cdef int i, j
    for i from 0 <= i < depth:
        cgd.object_stack[i] = allocate_dg_vert(i+degree,depth+degree-1)
        cgd.parent_stack[i] = allocate_dg_vert(i+degree,depth+degree-1)
        if cgd.object_stack[i]        is NULL or \
           cgd.parent_stack[i]        is NULL:
            for j from 0 <= j <= i:
                free_dg_vert(cgd.object_stack[j])
                free_dg_vert(cgd.parent_stack[j])
            sig_free(dg_vert_gen)
            deallocate_cgd(cgd)
            return NULL
    for i from 0 <= i < depth-1:
                               # TODO: in fact, should this not happen in
                               # dg_vert_gen_children!? otherwise iterator[i].data will be NULL
                               # and no problems.....
        if allocate_subset_gen_2(i+degree, i+degree, cgd.iterator_stack+i):
            for j from 0 <= j < depth:
                free_dg_vert(cgd.object_stack[j])
                free_dg_vert(cgd.parent_stack[j])
            for j from 0 <= j < i:
                cgd2 = <canonical_generator_data *> cgd.iterator_stack[j].data
                deallocate_cgd(cgd2)
            sig_free(dg_vert_gen)
            deallocate_cgd(cgd)
            return NULL
    dg_vert_gen.data = <void *> cgd
    dg_vert_gen.next = canonical_generator_next
    return dg_vert_gen

cdef void free_dg_vert_gen(iterator *dg_vert_gen) noexcept:
    r"""
    Deallocate the iterator for generating graphs.
    """
    cdef canonical_generator_data *cgd = <canonical_generator_data *> dg_vert_gen.data
    deallocate_cgd(cgd)
    sig_free(dg_vert_gen)

cdef void free_cgd_2(void *data) noexcept:
    r"""
    A simpler alternative to ``free_cgd``.
    """
    cdef canonical_generator_data *cgd = <canonical_generator_data *> data
    deallocate_cgd(cgd)


def generate_dense_graphs_vert_addition(int n, base_G=None,
                                        bint construct=False,
                                        bint indicate_mem_err=True):
    r"""

    EXAMPLES::

        sage: from sage.groups.perm_gps.partn_ref.refinement_graphs import generate_dense_graphs_vert_addition

    ::

        sage: for n in [0..7]:
        ....:   generate_dense_graphs_vert_addition(n)
        1
        2
        4
        8
        19
        53
        209
        1253
        sage: generate_dense_graphs_vert_addition(8)  # long time
        13599

    TESTS::

        sage: from sage.groups.perm_gps.partn_ref.refinement_graphs import generate_dense_graphs_vert_addition
        sage: generate_dense_graphs_vert_addition(10, base_G=Graph('HEhf^rs'))
        11
    """
    from sage.graphs.graph import Graph
    cdef iterator *graph_iterator
    cdef DenseGraph DG, ODG
    cdef GraphStruct GS
    if n < 2:
        if construct:
            L = []
            if n < 0:
                return L
            L.append(Graph(0, implementation='c_graph', sparse=False))
            if n == 0:
                return L
            L.append(Graph(0, implementation='c_graph', sparse=False))
            L.reverse()
            return L
        else:
            if n < 0:
                return Integer(0)
            if n == 0:
                return Integer(1)
            return Integer(2)

    cdef int start_deg = 1 if base_G is None else base_G.num_verts()
    graph_iterator = allocate_dg_vert_gen(start_deg, n+1-start_deg)
    if graph_iterator is NULL:
        raise MemoryError

    GS = (<GraphStruct> (<canonical_generator_data *> graph_iterator.data).object_stack[0])
    DG = GS.G
    if base_G is not None:
        for v in base_G.vertices(sort=False):
            DG.add_vertex(v)
        for u,v in base_G.edges(sort=False, labels=False):
            DG.add_arc(u,v)
            DG.add_arc(v,u)

    graph_iterator = setup_canonical_generator(start_deg,
        all_children_are_equivalent,
        refine_by_degree,
        compare_graphs,
        gen_children_dg_vert,
        apply_dg_vert_aug,
        free_dg_vert,
        free_cgd_2,
        free_subset,
        canonical_dg_vert_parent,
        n+1-start_deg, 0, graph_iterator)

    start_canonical_generator(NULL, <void *> GS, DG.num_verts, graph_iterator)

    cdef list out_list
    cdef void *thing
    cdef Integer number
    cdef bint mem_err = 0
    if construct:
        out_list = []
    else:
        number = Integer(0)
    if construct:
        while True:
            thing = graph_iterator.next(graph_iterator.data, NULL, &mem_err)
            if thing is NULL:
                break
            ODG = (<GraphStruct>thing).G
            G = Graph(0, implementation='c_graph', sparse=False)
            DG = DenseGraph(ODG.active_vertices.size, extra_vertices=0)
            copy_dense_graph(DG, ODG)
            G._backend._cg = DG
            out_list.append(G)
    else:
        while True:
            thing = graph_iterator.next(graph_iterator.data, NULL, &mem_err)
            if thing is NULL:
                break
            number += 1

    free_dg_vert_gen(graph_iterator)
    if mem_err:
        if indicate_mem_err:
            raise MemoryError
        else:
            out_list.append(MemoryError())
    if construct:
        if base_G is None:
            out_list = [Graph(0, implementation='c_graph', sparse=False)] + out_list
        return out_list
    else:
        if base_G is None:
            number += 1
        return number
