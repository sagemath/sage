# cython: binding=True
# distutils: language = c++
r"""
Decomposition by clique minimal separators

This module implements methods related to the decomposition of a graph by clique
minimal separators. See [TY1984]_ and [BPS2010]_ for more details on the
algorithms.

Methods
-------
"""
# ****************************************************************************
# Copyright (C) 2019 David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from libcpp.vector cimport vector
from libc.stdint cimport uint32_t
from cysignals.signals cimport sig_on, sig_off, sig_check
from memory_allocator cimport MemoryAllocator

from sage.graphs.base.static_sparse_graph cimport short_digraph
from sage.graphs.base.static_sparse_graph cimport init_short_digraph
from sage.graphs.base.static_sparse_graph cimport free_short_digraph
from sage.graphs.base.static_sparse_graph cimport has_edge

from sage.sets.set import Set
from sage.graphs.traversals cimport maximum_cardinality_search_M_short_digraph


def make_tree(atoms, cliques):
    r"""
    Return a tree of atoms and cliques.

    The atoms are the leaves of the tree and the cliques are the internal
    vertices. The number of atoms is the number of cliques plus one.

    As a clique may appear several times in the list ``cliques``, vertices are
    numbered by pairs `(i, S)`, where `0 \leq i < |atoms| + |cliques|` and `S`
    is either an atom or a clique.

    The root of the tree is the only vertex with even or null degree, i.e., 0 if
    ``cliques`` is empty and 2 otherwise. When ``cliques`` is not empty, other
    internal vertices (each of which is a clique) have degree 3, and the
    leaves (vertices of degree 1) are the atoms.

    INPUT:

    - ``atoms`` -- list of atoms

    - ``cliques`` -- list of cliques

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.clique_separators import make_tree
        sage: G = graphs.Grid2dGraph(2, 4)
        sage: A, Sc = G.atoms_and_clique_separators()
        sage: T = make_tree(A, Sc)
        sage: all(u[1] in A for u in T if T.degree(u) == 1)
        True
        sage: all(u[1] in Sc for u in T if T.degree(u) != 1)
        True

    TESTS::

        sage: from sage.graphs.graph_decompositions.clique_separators import make_tree
        sage: make_tree([0], [1, 2])
        Traceback (most recent call last):
        ...
        ValueError: the number of atoms must be the number of cliques plus one
    """
    if (atoms or cliques) and len(atoms) != len(cliques) + 1:
        raise ValueError("the number of atoms must be the number of cliques plus one")

    from sage.graphs.graph import Graph
    T = Graph()
    if cliques:
        # As a clique can appear several times, we number the vertices by
        # pairs (int, Set), with 0 <= int < |atoms| + |cliques|
        T.add_path(list(enumerate(cliques)))
        j = len(cliques)
        T.add_edges((s, (i + j, a)) for s, (i, a) in zip(enumerate(cliques), enumerate(atoms)))
        # We have |atoms| = |cliques| + 1. So |atoms| + |cliques| = 2 * j + 1
        T.add_edge((j - 1, cliques[-1]), (2 * j, atoms[-1]))

    elif atoms:
        # The graph has no clique separator
        T.add_vertex(atoms[0])

    return T


def make_labelled_rooted_tree(atoms, cliques):
    r"""
    Return a :class:`~LabelledRootedTree` of atoms and cliques.

    The atoms are the leaves of the tree and the cliques are the internal
    vertices. The number of atoms is the number of cliques plus one.

    EXAMPLES::

        sage: G = graphs.PathGraph(5)
        sage: ascii_art(G.atoms_and_clique_separators(rooted_tree=True))
          _____{3}_____
         /            /
        {3, 4}   ____{2}____
                /          /
               {2, 3}   __{1}__
                       /      /
                      {0, 1} {1, 2}

    TESTS::

        sage: from sage.graphs.graph_decompositions.clique_separators import make_labelled_rooted_tree
        sage: make_labelled_rooted_tree([0], [1, 2])
        Traceback (most recent call last):
        ...
        ValueError: the number of atoms must be the number of cliques plus one
    """
    from sage.combinat.rooted_tree import LabelledRootedTree
    if not atoms and not cliques:
        return LabelledRootedTree([])

    if len(atoms) != len(cliques) + 1:
        raise ValueError("the number of atoms must be the number of cliques plus one")

    def to_tree(i, n):
        if i < n:
            return LabelledRootedTree([LabelledRootedTree([], label=atoms[i]), to_tree(i + 1, n)],
                                      label=cliques[i])
        return LabelledRootedTree([], label=atoms[i])

    return to_tree(0, len(cliques))


cdef inline bint is_clique(short_digraph sd, vector[int] Hx) noexcept:
    """
    Check if the subgraph sd[Hx] is a clique.

    This is a helper function of ``atoms_and_clique_separators``.
    """
    cdef size_t Hx_size = Hx.size()
    cdef size_t i, j
    cdef int u
    for i in range(Hx_size - 1):
        u = Hx[i]
        for j in range(i + 1, Hx_size):
            if not has_edge(sd, u, Hx[j]):
                return False
    return True


def atoms_and_clique_separators(G, tree=False, rooted_tree=False, separators=False):
    r"""
    Return the atoms of the decomposition of `G` by clique minimal separators.

    Let `G = (V, E)` be a graph. A set `S \subset V` is a clique separator if
    `G[S]` is a clique and the graph `G \setminus S` has at least 2 connected
    components. Let `C \subset V` be the vertices of a connected component of `G
    \setminus S`. The graph `G[C + S]` is an *atom* if it has no clique
    separator.

    This method implements the algorithm proposed in [BPS2010]_, that improves
    upon the algorithm proposed in [TY1984]_, for computing the atoms and the
    clique minimal separators of a graph. This algorithm is based on the
    :meth:`~sage.graphs.traversals.maximum_cardinality_search_M` graph traversal
    and has time complexity in `O(|V|\cdot|E|)`.

    If the graph is not connected, we insert empty separators between the lists
    of separators of each connected components. See the examples below for more
    details.

    INPUT:

    - ``G`` -- a Sage graph

    - ``tree`` -- boolean (default: ``False``); whether to return the result as
      a directed tree in which internal nodes are clique separators and leaves
      are the atoms of the decomposition. Since a clique separator is repeated
      when its removal partition the graph into 3 or more connected components,
      vertices are labels by tuples `(i, S)`, where `S` is the set of vertices
      of the atom or the clique separator, and `0 \leq i \leq |T|`.

    - ``rooted_tree`` -- boolean (default: ``False``); whether to return the
      result as a :class:`~sage.combinat.rooted_tree.LabelledRootedTree`. When
      ``tree`` is ``True``, this parameter is ignored.

    - ``separators`` -- boolean (default: ``False``); whether to also return the
      complete list of separators considered during the execution of the
      algorithm. When ``tree`` or ``rooted_tree`` is ``True``, this parameter is
      ignored.

    OUTPUT:

    - By default, return a tuple `(A, S_c)`, where `A` is the list of atoms of
      the graph in the order of discovery, and `S_c` is the list of clique
      separators, with possible repetitions, in the order the separator has been
      considered. If furthermore ``separators`` is ``True``, return a tuple `(A,
      S_h, S_c)`, where `S_c` is the list of considered separators of the graph
      in the order they have been considered.

    - When ``tree`` is ``True``, format the result as a directed tree

    - When ``rooted_tree`` is ``True`` and ``tree`` is ``False``, format the
      output as a :class:`~sage.combinat.rooted_tree.LabelledRootedTree`

    EXAMPLES:

    Example of [BPS2010]_::

        sage: G = Graph({'a': ['b', 'k'], 'b': ['c'], 'c': ['d', 'j', 'k'],
        ....:            'd': ['e', 'f', 'j', 'k'], 'e': ['g'],
        ....:            'f': ['g', 'j', 'k'], 'g': ['j', 'k'], 'h': ['i', 'j'],
        ....:            'i': ['k'], 'j': ['k']})
        sage: atoms, cliques = G.atoms_and_clique_separators()
        sage: sorted(sorted(a) for a in atoms)
        [['a', 'b', 'c', 'k'],
         ['c', 'd', 'j', 'k'],
         ['d', 'e', 'f', 'g', 'j', 'k'],
         ['h', 'i', 'j', 'k']]
        sage: sorted(sorted(c) for c in cliques)
        [['c', 'k'], ['d', 'j', 'k'], ['j', 'k']]
        sage: T = G.atoms_and_clique_separators(tree=True)
        sage: T.is_tree()
        True
        sage: T.diameter() == len(atoms)
        True
        sage: all(u[1] in atoms for u in T if T.degree(u) == 1)
        True
        sage: all(u[1] in cliques for u in T if T.degree(u) != 1)
        True

    A graph without clique separator::

        sage: G = graphs.CompleteGraph(5)
        sage: G.atoms_and_clique_separators()
        ([{0, 1, 2, 3, 4}], [])
        sage: ascii_art(G.atoms_and_clique_separators(rooted_tree=True))
        {0, 1, 2, 3, 4}

    Graphs with several biconnected components::

        sage: G = graphs.PathGraph(4)
        sage: ascii_art(G.atoms_and_clique_separators(rooted_tree=True))
          ____{2}____
         /          /
        {2, 3}   __{1}__
                /      /
               {1, 2} {0, 1}

        sage: G = graphs.WindmillGraph(3, 4)
        sage: G.atoms_and_clique_separators()
        ([{0, 1, 2}, {0, 3, 4}, {0, 5, 6}, {0, 8, 7}], [{0}, {0}, {0}])
        sage: ascii_art(G.atoms_and_clique_separators(rooted_tree=True))
          ________{0}________
         /                  /
        {0, 1, 2}   _______{0}______
                   /               /
                  {0, 3, 4}   ____{0}___
                             /         /
                            {0, 8, 7} {0, 5, 6}

    When the removal of a clique separator results in `k > 2` connected
    components, this separator is repeated `k - 1` times, but the repetitions
    are not necessarily contiguous::

        sage: G = Graph(2)
        sage: for i in range(5):
        ....:     G.add_cycle([0, 1, G.add_vertex()])
        sage: ascii_art(G.atoms_and_clique_separators(rooted_tree=True))
          _________{0, 1}_____
         /                   /
        {0, 1, 4}   ________{0, 1}_____
                   /                  /
                  {0, 1, 2}   _______{0, 1}___
                             /               /
                            {0, 1, 3}   ____{0, 1}
                                       /         /
                                      {0, 1, 5} {0, 1, 6}

        sage: G = graphs.StarGraph(3)
        sage: G.subdivide_edges(G.edges(sort=False), 2)
        sage: ascii_art(G.atoms_and_clique_separators(rooted_tree=True))
          ______{5}______
         /              /
        {1, 5}   ______{7}______
                /              /
               {2, 7}   ______{9}______
                       /              /
                      {9, 3}   ______{6}______
                              /              /
                             {6, 7}   ______{4}_____
                                     /             /
                                    {4, 5}   _____{0}_____
                                            /            /
                                           {0, 6}   ____{8}____
                                                   /          /
                                                  {8, 9}   __{0}__
                                                          /      /
                                                         {0, 8} {0, 4}

    If the graph is not connected, we insert empty separators between the lists
    of separators of each connected components. For instance, let `G` be a graph
    with 3 connected components. The method returns the list `S_c =
    [S_0,\cdots,S_{i},\ldots, S_{j},\ldots,S_{k-1}]` of `k` clique separators,
    where `i` and `j` are the indexes of the inserted empty separators and `0
    \leq i < j < k - 1`. The method also returns the list `A =
    [A_0,\ldots,S_{k}]` of the `k + 1` atoms, with `k + 1 \geq 3`. The lists of
    atoms and clique separators of each of the connected components are
    respectively `[A_0,\ldots,A_{i}]` and `[S_0,\ldots,S_{i-1}]`,
    `[A_{i+1},\ldots,A_{j}]` and `[S_{i+1},\ldots,S_{j-1}]`, and
    `[A_{j+1},\ldots,A_{k}]` and `[S_{j+1},\ldots,S_{k-1}]`. One can check that
    for each connected component, we get one atom more than clique separators::

        sage: G = graphs.PathGraph(3) * 3
        sage: A, Sc = G.atoms_and_clique_separators()
        sage: A
        [{1, 2}, {0, 1}, {4, 5}, {3, 4}, {8, 7}, {6, 7}]
        sage: Sc
        [{1}, {}, {4}, {}, {7}]
        sage: i , j = [i for i, s in enumerate(Sc) if not s]
        sage: i, j
        (1, 3)
        sage: A[:i+1], Sc[:i]
        ([{1, 2}, {0, 1}], [{1}])
        sage: A[i+1:j+1], Sc[i+1:j]
        ([{4, 5}, {3, 4}], [{4}])
        sage: A[j+1:], Sc[j+1:]
        ([{8, 7}, {6, 7}], [{7}])
        sage: I = [-1, i, j, len(Sc)]
        sage: for i, j in zip(I[:-1], I[1:]):
        ....:     print(A[i+1:j+1], Sc[i+1:j])
        [{1, 2}, {0, 1}] [{1}]
        [{4, 5}, {3, 4}] [{4}]
        [{8, 7}, {6, 7}] [{7}]
        sage: ascii_art(G.atoms_and_clique_separators(rooted_tree=True))
          ______{1}______
         /              /
        {1, 2}   ______{}______
                /             /
               {0, 1}   _____{4}_____
                       /            /
                      {4, 5}   ____{}_____
                              /          /
                             {3, 4}   __{7}__
                                     /      /
                                    {6, 7} {8, 7}

    Loops and multiple edges are ignored::

        sage: G.allow_loops(True)
        sage: G.add_edges([(u, u) for u in G])
        sage: G.allow_multiple_edges(True)
        sage: G.add_edges(G.edges(sort=False))
        sage: ascii_art(G.atoms_and_clique_separators(rooted_tree=True))
          ______{1}______
         /              /
        {1, 2}   ______{}______
                /             /
               {0, 1}   _____{4}_____
                       /            /
                      {4, 5}   ____{}_____
                              /          /
                             {3, 4}   __{7}__
                                     /      /
                                    {6, 7} {8, 7}

    We can check that the returned list of separators is valid::

        sage: G = graphs.RandomGNP(50, .1)
        sage: while not G.is_connected():
        ....:     G = graphs.RandomGNP(50, .1)
        sage: _, separators, _ = G.atoms_and_clique_separators(separators=True)
        sage: for S in separators:
        ....:     H = G.copy()
        ....:     H.delete_vertices(S)
        ....:     if H.is_connected():
        ....:         raise ValueError("something goes wrong")

    TESTS::

        sage: EmptyGraph = Graph()
        sage: EmptyGraph.atoms_and_clique_separators()
        ([], [])
        sage: EmptyGraph.atoms_and_clique_separators(separators=True)
        ([], [], [])
        sage: EmptyGraph.atoms_and_clique_separators(tree=True)
        Graph on 0 vertices
        sage: EmptyGraph.atoms_and_clique_separators(rooted_tree=True)
        None[]
        sage: ascii_art(EmptyGraph.atoms_and_clique_separators(rooted_tree=True))
        None
        sage: I4 = Graph(4)
        sage: I4.atoms_and_clique_separators()
        ([{0}, {1}, {2}, {3}], [{}, {}, {}])
        sage: I4.atoms_and_clique_separators(separators=True)
        ([{0}, {1}, {2}, {3}], [{}, {}, {}], [{}, {}, {}])
        sage: I4.atoms_and_clique_separators(tree=True)
        Graph on 7 vertices
        sage: I4.atoms_and_clique_separators(rooted_tree=True)
        {}[{0}[], {}[{1}[], {}[{3}[], {2}[]]]]
        sage: ascii_art(I4.atoms_and_clique_separators(rooted_tree=True))
          ___{}___
         /       /
        {0}   __{}___
             /      /
            {1}   _{}_
                 /   /
                {3} {2}
    """
    cdef list A = []   # atoms
    cdef list Sh = []  # separators
    cdef list Sc = []  # clique separators
    cdef bint first = True

    if not G.is_connected():
        from sage.graphs.graph import Graph

        for cc in G.connected_components(sort=False):
            g = Graph([cc, G.edge_boundary(cc, cc, False, False)],
                      format='vertices_and_edges',
                      loops=True, multiedges=True)
            res = g.atoms_and_clique_separators(tree=False, rooted_tree=False, separators=separators)

            # Update lists of atoms, separators and clique separators
            A.extend(res[0])
            if separators:
                if not first:
                    Sh.append(Set())
                Sh.extend(res[1])
            if not first:
                Sc.append(Set())
            if separators:
                Sc.extend(res[2])
            else:
                Sc.extend(res[1])
            first = False

        # Format and return the result
        if tree:
            return make_tree(A, Sc)
        elif rooted_tree:
            return make_labelled_rooted_tree(A, Sc)
        elif separators:
            return A, Sh, Sc
        return A, Sc

    cdef int N = G.order()
    cdef list int_to_vertex = list(G)

    # Copying the whole graph to obtain the list of neighbors quicker than by
    # calling out_neighbors. This data structure is well documented in the
    # module sage.graphs.base.static_sparse_graph
    cdef short_digraph sd
    init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_vertex)

    # variables for the manipulation of the short digraph
    cdef uint32_t** p_vertices = sd.neighbors
    cdef uint32_t* p_tmp
    cdef uint32_t* p_end

    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int* alpha = <int*>mem.calloc(N, sizeof(int))
    cdef int* alpha_inv = <int*>mem.calloc(N, sizeof(int))
    cdef bint* X = <bint*>mem.calloc(N, sizeof(bint))
    cdef bint* active = <bint*>mem.calloc(N, sizeof(bint))
    cdef int* waiting_list = <int*>mem.calloc(N, sizeof(int))
    cdef int* seen = <int*>mem.calloc(N, sizeof(int))
    cdef list F = []
    cdef int u, v, waiting_begin, waiting_end

    sig_on()
    maximum_cardinality_search_M_short_digraph(sd, 0, alpha, alpha_inv, F, X)
    sig_off()

    # Instead of building the graph H of the triangulation, extracting the
    # neighbors of vertex alpha[i] and then removing that vertex from H, we
    # directly build the neighborhoods. Note that vertices are removed from the
    # graph of the triangulation in the order of alpha. Hence, at step i, this
    # graph has no vertex u such that alpha_inv[u] < i, and edge (u, v) is
    # removed from the graph when min(alpha_inv[u], alpha_inv[v]) is removed.
    # The neighborhood of x at step i is thus {v in N_H(x) | alpha_inv[v] > i}.
    cdef vector[vector[int]] H
    H.resize(N)
    for u in range(N):
        p_tmp = p_vertices[u]
        p_end = p_vertices[u + 1]
        while p_tmp < p_end:
            v = p_tmp[0]
            p_tmp += 1
            if u < v:
                # We consider edge (u, v) only once, when u > v, and the
                # short_digraph data structure ensures that neighbors are sorted
                break
            if alpha_inv[u] < alpha_inv[v]:
                if X[u]:
                    H[u].push_back(v)
            elif X[v]:
                H[v].push_back(u)
    for u, v in F:
        if alpha_inv[u] < alpha_inv[v]:
            if X[u]:
                H[u].push_back(v)
        elif X[v]:
            H[v].push_back(u)

    # Instead of using a copy Gp of G and removing from it the vertices of the
    # connected component of an atom after its discovery, we use an array of
    # booleans to avoid visiting inactive vertices
    for u in range(N):
        active[u] = True
        seen[u] = -1

    cdef frozenset Sint
    cdef vector[int] Sint_min
    cdef vector[int] Cint
    cdef vector[int] Hx

    for i in range(N):
        sig_check()
        x = alpha[i]
        if X[x] and not H[x].empty():
            Hx = H[x]

            if separators:
                Sh.append(Set(int_to_vertex[u] for u in Hx))

            if is_clique(sd, Hx):
                # The subgraph G[H[x]] = G[S] is a clique.
                # We extract the connected component of Gp - S containing x
                Sint = frozenset(Hx)
                Sint_min.clear()
                Cint.clear()
                Cint.push_back(x)
                seen[x] = x
                waiting_list[0] = x
                waiting_begin = 0
                waiting_end = 0
                while waiting_begin <= waiting_end:

                    u = waiting_list[waiting_begin]
                    waiting_begin += 1
                    p_tmp = p_vertices[u]
                    end = p_vertices[u + 1]

                    while p_tmp < end:
                        v = p_tmp[0]
                        p_tmp += 1

                        if active[v] and seen[v] != x:
                            seen[v] = x
                            if v in Sint:
                                # We keep only the vertices of the clique
                                # separator incident to the connected component
                                # containing x
                                Sint_min.push_back(v)
                            else:
                                Cint.push_back(v)
                                waiting_end += 1
                                waiting_list[waiting_end] = v

                # Store the atom Smin + C and the minimal clique separator Smin
                Smin = Set(int_to_vertex[u] for u in Sint_min)
                A.append(Set(Smin.set().union(int_to_vertex[u] for u in Cint)))
                Sc.append(Smin)

                # "Remove" the vertices of Cint from the graph Gp
                for u in Cint:
                    active[u] = False

    free_short_digraph(sd)
    H.clear()

    # We add the last atom
    if Sc:
        A.append(Set(int_to_vertex[x] for x in range(N) if active[x]))
    elif G:
        # The graph has no clique separator
        A.append(Set(int_to_vertex))

    # Format and return the result
    if tree:
        return make_tree(A, Sc)
    elif rooted_tree:
        return make_labelled_rooted_tree(A, Sc)
    if separators:
        return A, Sh, Sc
    return A, Sc
