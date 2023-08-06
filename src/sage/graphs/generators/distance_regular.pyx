r"""
Database of distance regular graphs

In this module we construct several distance regular graphs
and group them in a function that maps intersection arrays
to graphs.

For a survey on distance-regular graph see [BCN1989]_ or [VDKT2016]_.

EXAMPLES::

    sage: G = graphs.cocliques_HoffmannSingleton()
    sage: G.is_distance_regular()
    True
    sage: H = graphs.distance_regular_graph([15, 14, 10, 3, 1, 5, 12, 15])
    sage: H == G
    True
    sage: G = graphs.distance_regular_graph([27, 10, 1, 1, 10, 27])
    sage: G.is_distance_regular(True)
    ([27, 10, 1, None], [None, 1, 10, 27])

AUTHORS:

- Ivo Maffei (2020-07-28): initial version

"""

# ****************************************************************************
#       Copyright (C) 2020 Ivo Maffei <ivomaffei@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import itertools

from cysignals.signals cimport sig_check

from sage.graphs.graph import Graph
from sage.misc.lazy_import import LazyImport
from sage.graphs.generators.smallgraphs import (FosterGraph, BiggsSmithGraph,
                                                CoxeterGraph, LivingstoneGraph,
                                                WellsGraph, GossetGraph,
                                                HoffmanSingletonGraph,
                                                SimsGewirtzGraph,
                                                HigmanSimsGraph)
from sage.graphs.generators.platonic_solids import DodecahedralGraph
from sage.graphs.strongly_regular_db import strongly_regular_graph

codes = LazyImport('sage.coding', 'codes_catalog', as_name='codes')
libgap = LazyImport('sage.libs.gap.libgap', 'libgap')
Matrix = LazyImport('sage.matrix.constructor', 'Matrix')
VectorSpace = LazyImport('sage.modules.free_module', 'VectorSpace')
vector = LazyImport('sage.modules.free_module_element', 'vector')
GF = LazyImport('sage.rings.finite_rings.finite_field_constructor', 'GF')


def cocliques_HoffmannSingleton():
    r"""
    Return the graph obtained from the cocliques of the Hoffmann-Singleton graph.

    This is a distance-regular graph with intersection array
    `[15, 14, 10, 3; 1, 5, 12, 15]`.

    EXAMPLES::

        sage: G = graphs.cocliques_HoffmannSingleton()
        sage: G.is_distance_regular(True)
        ([15, 14, 10, 3, None], [None, 1, 5, 12, 15])

    REFERENCES:

    The construction of this graph can be found in [BCN1989]_ p. 392.
    """
    from sage.graphs.graph_generators import GraphGenerators

    D = GraphGenerators.HoffmanSingletonGraph()
    DC = D.complement()

    cocliques = [frozenset(c) for c in DC.cliques_maximum()]  # 100 of this

    edges = []
    for c1, c2 in itertools.combinations(cocliques, 2):
        if len(c1.intersection(c2)) == 8:
            edges.append((c1, c2))

    G = Graph(edges, format="list_of_edges")
    return G


def locally_GQ42_distance_transitive_graph():
    r"""
    Return the unique amply regular graph with `\mu = 6` which is locally
    a generalised quadrangle.

    This graph is distance-regular with intersection array
    `[45, 32, 12, 1; 1, 6, 32, 45]`.

    This graph is also distance-transitive.

    EXAMPLES::

        sage: G = graphs.locally_GQ42_distance_transitive_graph()  # optional - internet gap_packages
        sage: G.is_distance_regular(True)  # optional - internet gap_packages
        ([45, 32, 12, 1, None], [None, 1, 6, 32, 45])

    REFERENCES:

    A description of this graph can be found in [BCN1989]_ p.399.
    This construction is due to Dima Pasechnik.
    """
    H = libgap.AtlasGroup("3^2.U4(3).D8", libgap.NrMovedPoints, 756)
    Ns = H.NormalSubgroups()
    for N in Ns:
        if len(N.GeneratorsSmallest()) == 7:  # there is only one
            break

    G = Graph(libgap.Orbit(N, [1, 9], libgap.OnSets), format='list_of_edges')
    G.name("locally GQ(4,2) distance transitive graph")
    return G


def ConwaySmith_for_3S7():
    r"""
    Return the Conway-Smith graph related to `3 Sym(7)`.

    This is a distance-regular graph with intersection array
    `[10, 6, 4, 1; 1, 2, 6, 10]`.

    EXAMPLES::

        sage: G = graphs.ConwaySmith_for_3S7()                                          # needs sage.modules sage.rings.finite_rings sage.rings.number_field
        sage: G.is_distance_regular(True)                                               # needs sage.modules sage.rings.finite_rings sage.rings.number_field
        ([10, 6, 4, 1, None], [None, 1, 2, 6, 10])

    REFERENCES:

    A description and construction of this graph can be found in
    [BCN1989]_ p. 399.
    """
    from sage.rings.number_field.number_field import CyclotomicField

    F = CyclotomicField(3)
    w = F.gen()

    V = VectorSpace(GF(4), 6)
    z2 = GF(4)('z2')  # GF(4) = {0, 1, z2, z2+1}

    W = V.span([(0, 0, 1, 1, 1, 1), (0, 1, 0, 1, z2, z2 + 1), (1, 0, 0, 1, z2 + 1, z2)])
    # we only need the 45 vectors with 2 zero entries
    # we also embed everything into CC

    K = []
    for v in W:
        # check zero entries
        zeros = 0
        for x in v:
            if x.is_zero():
                zeros += 1

        if zeros == 2:
            # send to F and in K
            # z2 -> w
            # z2+1 -> w^2
            vv = []  # new vector
            for x in v:
                if x == z2:
                    vv.append(w)
                elif x == z2 + 1:
                    vv.append(w**2)
                else:
                    vv.append(int(x))

            # now vv is the new vector in F
            vv = vector(F, vv)
            K.append(vv)

    # we need to add other vectors
    for i in range(6):
        # create e_i
        ei = [0, 0, 0, 0, 0, 0]
        ei[i] = 1
        ei = vector(F, ei)

        K.append(2 * ei)
        K.append(2 * w * ei)
        K.append(2 * w**2 * ei)
    # now K is all the 63 vertices

    for v in K:
        v.set_immutable()

    def has_edge(u, v):
        return sum(u[i].conjugate() * v[i] for i in range(6)) == 2

    G = Graph()
    for Ki, Kj in itertools.combinations(K, 2):
        if has_edge(Ki, Kj):
            G.add_edge((Ki, Kj))

    G.name("Conway-Smith graph for 3S7")
    return G


def graph_3O73():
    r"""
    Return the graph related to the group `3 O(7,3)`.

    This graph is distance-regular with intersection array
    `[117, 80, 24, 1; 1, 12, 80, 117]`.

    The graph is also distance transitive with `3.O(7,3)` as automorphism
    group

    EXAMPLES::

        sage: G = graphs.graph_3O73()  # optional - internet gap_packages
        sage: G.is_distance_regular(True)  # optional - internet gap_packages
        ([117, 80, 24, 1, None], [None, 1, 12, 80, 117])

    REFERENCES:

    A description and construction of this graph can be found in
    [BCN1989]_ p. 400.
    """
    group = libgap.AtlasGroup("3.O7(3)", libgap.NrMovedPoints, 1134)
    G = Graph(libgap.Orbit(group, [1, 3], libgap.OnSets), format='list_of_edges')
    G.name("Distance transitive graph with automorphism group 3.O_7(3)")
    return G


def FosterGraph3S6():
    r"""
    Return the Foster graph for `3.Sym(6)`.

    This graph is distance-regular with intersection array
    `[6, 4, 2, 1; 1, 1, 4, 6]`.

    The graph is also distance transitive.

    EXAMPLES::

        sage: G = graphs.FosterGraph3S6()                                               # needs sage.libs.gap
        sage: G.is_distance_regular(True)                                               # needs sage.libs.gap
        ([6, 4, 2, 1, None], [None, 1, 1, 4, 6])

    REFERENCES:

    A description and construction of this graph can be found in
    [BCN1989]_ p. 397.
    """
    a = libgap.eval(("(2,6)(3,5)(4,11)(7,17)(8,16)(9,14)(13,22)(15,25)"
                    "(18,29)(19,28)(20,21)(24,30)(26,35)(27,33)(31,39)"
                     "(34,38)(36,43)(37,40)(42,44)"))
    b = libgap.eval(("(1,2,7,12,4)(3,8,18,20,10)(5,9,19,21,11)(6,13,17,26,15)"
                     "(14,23,28,31,24)(16,22,29,36,27)(25,32,35,42,34)"
                     "(30,37,39,44,38)(33,40,43,45,41)"))

    group = libgap.Group(a, b)

    G = Graph(group.Orbit([1, 7], libgap.OnSets), format='list_of_edges')
    G.name("Foster graph for 3.Sym(6) graph")
    return G


def J2Graph():
    r"""
    Return the distance-transitive graph with automorphism group `J_2`.

    EXAMPLES::

        sage: G = graphs.J2Graph()  # optional - internet gap_packages
        sage: G.is_distance_regular(True) # optional - internet gap_packages
        ([10, 8, 8, 2, None], [None, 1, 1, 4, 5])

    REFERENCES:

    A description and construction of this graph can be found in
    [BCN1989]_ p. 408.
    """
    group = libgap.AtlasGroup("J2", libgap.NrMovedPoints, 315)
    G = Graph(group.Orbit([1, 9], libgap.OnSets), format='list_of_edges')
    G.name("J_2 graph")
    return G


def IvanovIvanovFaradjevGraph():
    r"""
    Return the IvanovIvanovFaradjev graph.

    The graph is distance-transitive with automorphism group `3.M_{22}`.

    EXAMPLES::

        sage: G = graphs.IvanovIvanovFaradjevGraph()  # optional - internet gap_packages
        sage: G.is_distance_regular(True)  # optional - internet gap_packages
        ([7, 6, 4, 4, 4, 1, 1, 1, None], [None, 1, 1, 1, 2, 4, 4, 6, 7])

    REFERENCES:

    A description and construction of this graph can be found in
    [BCN1989]_ p. 369.
    """

    group = libgap.AtlasGroup("3.M22", libgap.NrMovedPoints, 990)
    graph = Graph(group.Orbit([1, 22], libgap.OnSets), format='list_of_edges')

    graph.name("Ivanov-Ivanov-Faradjev Graph")
    return graph


def LargeWittGraph():
    r"""
    Return the large Witt graph.

    This is a distance-regular graph with intersection array
    `[30,28,24;1,3,15]`.

    EXAMPLES::

        sage: g = graphs.LargeWittGraph()                                               # needs sage.modules
        sage: g.is_distance_regular(True)                                               # needs sage.modules
        ([30, 28, 24, None], [None, 1, 3, 15])

    REFERENCES:

    A description of this graph can be found in
    [BCN1989]_ p. 366.
    This construction is taken from
    http://mathworld.wolfram.com/LargeWittGraph.html
    """
    import itertools

    C = codes.GolayCode(GF(2), extended=True)
    vertices = [c for c in C if c.hamming_weight() == 8]

    edges = []
    for v, w in itertools.combinations(vertices, 2):
        if not set(v.support()).intersection(w.support()):
            edges.append((v, w))

    W = Graph(edges, format='list_of_edges')
    W.name("Large Witt graph")
    return W


def TruncatedWittGraph():
    r"""
    Return the truncated Witt graph.

    This builds the large Witt graph, then removes
    all vertices whose codeword start with a 1.

    The graph is distance-regular with intersection array
    `[15,14,12;1,1,9]`.

    EXAMPLES::

         sage: G = graphs.TruncatedWittGraph()  # long time                             # needs sage.modules
         sage: G.is_distance_regular(True)      # long time (due to above)              # needs sage.modules
         ([15, 14, 12, None], [None, 1, 1, 9])

    REFERENCES:

    A description and construction of this graph can be found in
    [BCN1989]_ p. 367.
    """
    # get large witt graph and remove all vertices which start with a 1
    G = LargeWittGraph()
    G.delete_vertices(filter(lambda x: x[0] == 1, G.vertices(sort=False)))

    G.name("Truncated Witt graph")
    return G


def DoublyTruncatedWittGraph():
    r"""
    Return the doubly truncated Witt graph.

    This builds the truncated Witt graph, then removes
    all vertices whose codeword start with a 1.

    The graph is distance-regular with intersection array
    `[7,6,4,4;1,1,1,6]`.

    EXAMPLES::

         sage: G = graphs.DoublyTruncatedWittGraph()                                    # needs sage.modules
         sage: G.is_distance_regular(True)                                              # needs sage.modules
         ([7, 6, 4, 4, None], [None, 1, 1, 1, 6])

    REFERENCES:

    A description and construction of this graph can be found in
    [BCN1989]_ p. 368.
    """
    G = TruncatedWittGraph()
    G.delete_vertices(filter(lambda x: x[1] == 1, G.vertices(sort=False)))

    G.name("Doubly Truncated Witt graph")
    return G


def distance_3_doubly_truncated_Golay_code_graph():
    r"""
    Return a distance-regular graph with intersection array
    `[9, 8, 6, 3; 1, 1, 3, 8]`.

    EXAMPLES::

        sage: G = graphs.distance_3_doubly_truncated_Golay_code_graph()         # long time, needs sage.modules sage.rings.finite_rings
        sage: G.is_distance_regular(True)       # long time (due to above)              # needs sage.modules sage.rings.finite_rings
        ([9, 8, 6, 3, None], [None, 1, 1, 3, 8])

    ALGORITHM:

    Compute the binary Golay code and truncate it twice. Compute its coset graph.
    Take a vertex and compute the set of vertices at distance 3
    from the vertex chosen. This set constitutes the set of vertices of our
    distance-regular graph. Moreover we have an edge `(u,v)` if the coset graph
    contains such edge.

    REFERENCES:

    Description and construction of this graph are taken from [BCN1989]_ p. 364.
    """
    G = codes.GolayCode(GF(2), extended=False).punctured([0, 1]).cosetGraph()
    v = G.vertices(sort=False)[0]
    it = G.breadth_first_search(v, distance=3, report_distance=True)
    vertices = [w for (w, d) in it if d == 3]

    edges = [(a, b) for a, b in itertools.combinations(vertices, 2)
             if G.has_edge((a, b))]

    H = Graph(edges, format='list_of_edges')
    return H


def shortened_00_11_binary_Golay_code_graph():
    r"""
    Return a distance-regular graph with intersection array
    `[21, 20, 16, 6, 2, 1; 1, 2, 6, 16, 20, 21]`.

    EXAMPLES::

        sage: G = graphs.shortened_00_11_binary_Golay_code_graph()      # long time (9 s), needs sage.modules sage.rings.finite_rings
        sage: G.is_distance_regular(True)       # long time                             # needs sage.modules sage.rings.finite_rings
        ([21, 20, 16, 6, 2, 1, None], [None, 1, 2, 6, 16, 20, 21])

    ALGORITHM:

    Compute the binary Golay code. Compute the subcode whose codewords start
    with 00 or 11. Remove the first two entries from all codewords of the newly
    found linear code and compute its coset graph.

    REFERENCES:

    Description and construction of this graph can be found in [BCN1989]_ p. 365.
    """
    from sage.coding.linear_code import LinearCode

    code = codes.GolayCode(GF(2), False)
    C_basis = code.basis()

    # Now special shortening
    v = C_basis[0] + C_basis[1]  # v has 11 at the start
    C_basis = C_basis[2:]
    C_basis.append(v)
    C_basis = list(map(lambda x: x[2:], C_basis))

    code = LinearCode(Matrix(GF(2), C_basis))

    G = code.cosetGraph()
    G.name("Shortened 00 11 binary Golay code")
    return G


def shortened_000_111_extended_binary_Golay_code_graph():
    r"""
    Return a distance-regular graph with intersection array
    `[21, 20, 16, 9, 2, 1; 1, 2, 3, 16, 20, 21]`.

    EXAMPLES::

        sage: G = graphs.shortened_000_111_extended_binary_Golay_code_graph()   # long time (25 s), needs sage.modules sage.rings.finite_rings
        sage: G.is_distance_regular(True)       # long time                             # needs sage.modules sage.rings.finite_rings
        ([21, 20, 16, 9, 2, 1, None], [None, 1, 2, 3, 16, 20, 21])

    ALGORITHM:

    Compute the extended binary Golay code. Compute its subcode whose codewords
    start with 000 or 111. Remove the first 3 entries from all the codewords
    from the new linear code and compute its coset graph.

    REFERENCES:

    Description and construction of this graph can be found in [BCN1989]_ p. 365.
    """
    from sage.coding.linear_code import LinearCode

    code = codes.GolayCode(GF(2))
    C_basis = code.basis()

    # now special shortening
    v = C_basis[0] + C_basis[1] + C_basis[2]  # v has 111 at the start
    C_basis = C_basis[3:]
    C_basis.append(v)
    C_basis = list(map(lambda x: x[3:], C_basis))

    code = LinearCode(Matrix(GF(2), C_basis))

    G = code.cosetGraph()
    G.name("Shortened 000 111 extended binary Golay code")
    return G


def vanLintSchrijverGraph():
    r"""
    Return the van Lint-Schrijver graph.

    The graph is distance-regular with intersection array
    `[6, 5, 5, 4; 1, 1, 2, 6]`.

    EXAMPLES::

         sage: G = graphs.vanLintSchrijverGraph()                                       # needs sage.modules
         sage: G.is_distance_regular(True)                                              # needs sage.modules
         ([6, 5, 5, 4, None], [None, 1, 1, 2, 6])

    REFERENCES:

    For a description of this graph see [BCN1989]_ p. 373.
    """
    from sage.coding.linear_code import LinearCode

    one = vector(GF(3), [1, 1, 1, 1, 1, 1])
    G = LinearCode(Matrix(GF(3), one)).cosetGraph()

    vertices = [v for v in G.vertices(sort=False) if v.dot_product(one) in {1, 2}]
    edges = [(v, w) for v, w in itertools.combinations(vertices, 2)
             if G.has_edge((v, w))]

    H = Graph(edges, format='list_of_edges')
    H.name("Linst-Schrijver graph")
    return H


def LeonardGraph():
    r"""
    Return the Leonard graph.

    The graph is distance-regular with intersection array
    `[12, 11, 10, 7; 1, 2, 5, 12]`.

    EXAMPLES::

         sage: G = graphs.LeonardGraph()                                                # needs sage.combinat sage.modules
         sage: G.is_distance_regular(True)                                              # needs sage.combinat sage.modules
         ([12, 11, 10, 7, None], [None, 1, 2, 5, 12])

    REFERENCES:

    For a description of this graph see [BCN1989]_ p. 371.
    """
    from sage.combinat.matrices.hadamard_matrix import hadamard_matrix

    M = hadamard_matrix(12)
    edges = []
    for i, j, k, l in itertools.product(range(12), repeat=4):
        if i == k or j == l:
            continue
        if M[i, j] * M[i, l] * M[k, j] * M[k, l] == -1:
            edges.append(((i, j), (k, l)))

    D = Graph(edges, format="list_of_edges")
    blocks = [frozenset(cl) for cl in D.cliques_maximum()]

    edges = [(p, b) for b in blocks for p in b]
    G = Graph(edges, format="list_of_edges")
    return G


def UstimenkoGraph(const int m, const int q):
    r"""
    Return the Ustimenko graph with parameters `(m, q)`.

    This is the distance 1 or 2 graph of the dual polar graph `C_{m-1}(q)`.
    The graph is distance-regular with classical with parameters
    `(d,q^2, qbinom(3,1,q) -1, qbinom(m+1,1,q) -1)`

    INPUT:

    - ``m, q`` -- integers; ``q`` must be a prime power and ``m > 1``.

    EXAMPLES::

        sage: G = graphs.UstimenkoGraph(4, 2)                                           # needs sage.libs.gap
        sage: G.is_distance_regular(True)                                               # needs sage.libs.gap
        ([70, 32, None], [None, 1, 35])

    REFERENCES:

    See [BCN1989]_ p. 279 or [VDKT2016]_ p. 22.

    TESTS::

        sage: # long time, needs sage.libs.gap
        sage: G = graphs.UstimenkoGraph(5, 2)
        sage: G.order()
        2295
        sage: G.is_distance_regular(True)
        ([310, 224, None], [None, 1, 35])
        sage: G = graphs.UstimenkoGraph(4,3)
        sage: G.is_distance_regular(True)
        ([390, 243, None], [None, 1, 130])
    """
    from sage.graphs.graph_generators import graphs

    G = graphs.SymplecticDualPolarGraph(2*m - 2, q)

    edgesToAdd = []
    for v in G:
        for w in G.neighbor_iterator(v):
            for u in G.neighbor_iterator(w):
                sig_check()
                if u != v and not G.has_edge(u, v):
                    # then u,v are at distance 2
                    edgesToAdd.append((u, v))

    G.add_edges(edgesToAdd)
    G.name(f"Ustimenko graph ({m}, {q})")
    return G


def BilinearFormsGraph(const int d, const int e, const int q):
    r"""
    Return a bilinear forms graph with the given parameters.

    This builds a graph whose vertices are all `d`x`e` matrices over
    `GF(q)`. Two vertices are adjacent if the difference of the two
    matrices has rank 1.

    The graph is distance-regular with classical parameters
    `(\min(d, e), q, q-1 , q^{\max(d, e)}-1)`.

    INPUT:

    - ``d, e`` -- integers; dimension of the matrices
    - ``q`` -- integer; a prime power

    EXAMPLES::

        sage: # needs sage.modules
        sage: G = graphs.BilinearFormsGraph(3, 3, 2)
        sage: G.is_distance_regular(True)
        ([49, 36, 16, None], [None, 1, 6, 28])
        sage: G = graphs.BilinearFormsGraph(3,3,3)      # not tested (20 s)             # needs sage.rings.finite_rings
        sage: G.order()                         # not tested (due to above)             # needs sage.rings.finite_rings
        19683
        sage: G = graphs.BilinearFormsGraph(3, 4, 2)    # long time                     # needs sage.rings.finite_rings
        sage: G.is_distance_regular(True)       # long time                             # needs sage.rings.finite_rings
        ([105, 84, 48, None], [None, 1, 6, 28])

    REFERENCES:

    See [BCN1989]_ pp. 280-282 for a rather detailed discussion, otherwise
    see [VDKT2016]_ p. 21.

    TESTS::

        sage: # needs sage.modules
        sage: G = graphs.BilinearFormsGraph(2,3,2)
        sage: G.is_distance_regular(True)
        ([21, 12, None], [None, 1, 6])
        sage: H = graphs.BilinearFormsGraph(3,2,2)
        sage: H.is_isomorphic(G)
        True
        sage: G = graphs.BilinearFormsGraph(5, 1, 3)
        sage: K = graphs.CompleteGraph(G.order())
        sage: K.is_isomorphic(G)
        True
    """
    from itertools import product as cartprod

    Fq = GF(q)
    Fqelems = list(Fq)
    FqToInt = {x: n for n, x in enumerate(Fqelems)}
    dim = d * e
    matricesOverq = cartprod(range(q), repeat=dim)
    qto = [int(q**jj) for jj in range(dim)]

    rank1Matrices = []
    for u in VectorSpace(Fq, d):
        if u.is_zero() or not u[u.support()[0]].is_one():
            continue

        for v in VectorSpace(Fq, e):
            if v.is_zero():
                continue

            sig_check()
            M = [0] * dim
            for row in range(d):
                for col in range(e):
                    M[e*row + col] = u[row] * v[col]

            rank1Matrices.append(M)

    edges = []
    for m1 in matricesOverq:
        intM1 = 0  # represents vector m1 as integer base q
        for jj in range(dim):
            intM1 += m1[jj] * qto[jj]

        for m2 in rank1Matrices:
            sig_check()
            intM3 = 0
            for jj in range(dim):
                intM3 += FqToInt[Fqelems[m1[jj]] + Fqelems[m2[jj]]] * qto[jj]

            edges.append((intM1, intM3))

    G = Graph(edges, format='list_of_edges')
    G.name("Bilinear forms graphs over F_%d with parameters (%d, %d)" % (q, d, e))
    return G


def AlternatingFormsGraph(const int n, const int q):
    r"""
    Return the alternating forms graph with the given parameters.

    This builds a graph whose vertices are all `n`x`n` skew-symmetric
    matrices over `GF(q)` with zero diagonal. Two vertices are adjacent
    if and only if the difference of the two matrices has rank 2.

    This grap is distance-regular with classical parameters
    `(\lfloor \frac n 2 \rfloor,  q^2, q^2 - 1, q^{2 \lceil \frac n 2 \rceil -1})`.

    INPUT:

    - ``n`` -- integer
    - ``q`` -- a prime power

    EXAMPLES::

        sage: G = graphs.AlternatingFormsGraph(5, 2)  # long time
        sage: G.is_distance_regular(True)  # long time
        ([155, 112, None], [None, 1, 20])

    REFERENCES:

    See [BCN1989]_ pp. 282-284 for a rather detailed discussion, otherwise
    see [VDKT2016]_ p. 22.

    TESTS::

         sage: # needs sage.modules
         sage: G = graphs.AlternatingFormsGraph(6,2)    # not tested (2 min)            # needs sage.rings.finite_rings
         sage: G.order()                        # not tested (because of above)         # needs sage.rings.finite_rings
         32768
         sage: G.is_distance_regular(True)      # not tested (33 min)                   # needs sage.rings.finite_rings
         ([651, 560, 256, None], [None, 1, 20, 336])
         sage: G = graphs.AlternatingFormsGraph(4, 3)
         sage: G.is_distance_regular(True)
         ([260, 162, None], [None, 1, 90])
    """
    # n x n zero-diagonal skew-symmetric matrix
    # can be represented by the upper triangular entries
    # there are n*(n-1) // 2 of them
    size = (n * (n-1)) // 2
    V = VectorSpace(GF(q), size)

    # construct all rank 2 matrices
    rank2Matrices = set()
    Vn = VectorSpace(GF(q), n)
    basis = set(Vn.basis())
    e = [Vn([0]*i + [1] + [0]*(n - i - 1)) for i in range(n)]
    for v in e:
        v.set_immutable()

    scalars = [x for x in GF(q) if not x.is_zero()]
    Vseen = set()
    for v in Vn:
        if v.is_zero() or not v[v.support()[0]].is_one():
            continue
        v.set_immutable()
        # remove from basis e_i s.t. (v[i-1] =) v_i != 0
        i = v.support()[0]
        Ubasis = basis.difference([e[i]])

        for u in Vn.span_of_basis(Ubasis):
            sig_check()
            if u.is_zero() or not u[u.support()[0]].is_one():
                continue
            u.set_immutable()
            if u in Vseen:
                continue

            M = []
            for row in range(n - 1):
                upperRow = [0] * (n - 1 - row)
                for col in range(row + 1, n):
                    upperRow[col - row - 1] = v[row]*u[col] - u[row]*v[col]
                M += upperRow

            for scalar in scalars:
                N = tuple(map(lambda x: scalar * x, M))
                rank2Matrices.add(N)

        Vseen.add(v)

    # now we have all matrices of rank 2
    edges = []
    for m1 in V:
        t1 = tuple(m1)
        for m2 in rank2Matrices:
            sig_check()
            t3 = tuple([t1[i] + m2[i] for i in range(size)])
            edges.append((t1, t3))

    G = Graph(edges, format='list_of_edges')
    G.name("Alternating forms graph on (F_%d)^%d" % (q, n))
    return G


def HermitianFormsGraph(const int n, const int r):
    r"""
    Return the Hermitian forms graph with the given parameters.

    We build a graph whose vertices are all ``n``x``n`` Hermitian matrices
    over ``GF(r^2)``. Two  vertices are adjacent if the difference of the two
    vertices has rank 1.

    This graph is distance-regular with classical parameters
    `(n, - r, - r - 1, - (- r)^d - 1)`.

    INPUT:

    - ``n`` -- integer
    - ``r`` -- a prime power

    EXAMPLES::

        sage: # needs sage.modules sage.rings.finite_rings
        sage: G = graphs.HermitianFormsGraph(2, 2)
        sage: G.is_distance_regular(True)
        ([5, 4, None], [None, 1, 2])
        sage: G = graphs.HermitianFormsGraph(3, 3)      # not tested (2 min)
        sage: G.order()                         # not tested (bacuase of the above)
        19683

    REFERENCES:

    See [BCN1989]_ p. 285 or [VDKT2016]_ p. 22.

    TESTS::

         sage: # needs sage.modules sage.rings.finite_rings
         sage: G = graphs.HermitianFormsGraph(3, 2)
         sage: G.is_distance_regular(True)
         ([21, 20, 16, None], [None, 1, 2, 12])
         sage: G = graphs.HermitianFormsGraph(2, 3)
         sage: G.is_distance_regular(True)
         ([20, 18, None], [None, 1, 6])
    """
    q = r * r
    Fr = GF(r)
    Fq = GF(q)
    i = Fq.gen()
    ir = i**r

    toR = {(a + i*b): (a + ir*b) for a, b in itertools.product(Fr, repeat=2)}

    def build_mat(v, w):
        # get upper diagonal entries
        res = []
        used_v = 0
        used_w = 0
        for row in range(n):
            res += [v[used_v]] + [v[used_v + 1 + j] + i * w[used_w + j]
                                  for j in range(n - 1 - row)]
            used_v += n - row
            used_w += n - 1 - row

        return tuple(res)

    # produce all rank1 matrices
    rank1Matrices = []
    for w1 in VectorSpace(Fr, n):
        if not w1.is_zero():
            # build matrix
            nonZero = 0
            while nonZero < n and w1[nonZero] == 0:
                nonZero += 1

            for w2 in VectorSpace(Fr, n - nonZero - 1):
                # get upper triangular entries
                sig_check()

                v = [w1[nonZero]] + \
                    [w1[nonZero + 1 + j] + i * w2[j]
                     for j in range(n - nonZero - 1)]

                res = []
                for row in range(nonZero):
                    res += [0] * (n - row)

                res += v

                for row in range(1, n - nonZero):
                    factor = toR[v[row]] / v[0]
                    res += list(map(lambda x: factor * x, v[row:]))

                rank1Matrices.append(res)

    Vs = VectorSpace(Fr, (n * (n+1)) // 2)
    Va = VectorSpace(Fr, (n * (n-1)) // 2)

    edges = []
    for a, b in itertools.product(Vs, Va):
        M = build_mat(a, b)
        for R in rank1Matrices:
            N = tuple([M[i] + R[i] for i in range((n * (n+1)) // 2)])
            edges.append((M, N))

    G = Graph(edges, format='list_of_edges')
    G.name(f"Hermitian forms graph on (F_{q})^{n}")
    return G


def DoubleOddGraph(const int n):
    r"""
    Return the double odd graph on `2n+1` points.

    The graph is obtained using the subsets of size `n` and `n+1`
    of `{1, 2, ..., 2n+1}` as vertices. Two vertices are adjacent if one
    is included in the other.

    The graph is distance-transitive.

    INPUT:

    - ``n`` -- integer; must be greater than 0

    EXAMPLES::

         sage: G = graphs.DoubleOddGraph(5)
         sage: G.is_distance_regular(True)
         ([6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, None],
          [None, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6])
         sage: G = graphs.DoubleOddGraph(3)
         sage: G.diameter()
         7
         sage: G.is_distance_regular(True)
         ([4, 3, 3, 2, 2, 1, 1, None], [None, 1, 1, 2, 2, 3, 3, 4])

    REFERENCES:

    See [BCN1989]_ pp. 259-261 or [VDKT2016]_ p. 25.

    TESTS:

    DoubleOddGraph is bipartite double of OddGraph::

         sage: H = graphs.OddGraph(4)
         sage: G1 = graphs.DoubleOddGraph(3)
         sage: vertices = [(x, 0) for x in H] + [(x, 1) for x in H]
         sage: G2 = Graph([vertices, lambda i, j:
         ....: i[1] != j[1] and H.has_edge(i[0], j[0])])
         sage: G2.is_isomorphic(G1)
         True
    """
    from sage.combinat.integer_vector import IntegerVectors

    if n < 1:
        raise ValueError("n must be >= 1")

    cdef list edges, s1
    cdef int i

    # a binary vector of size 2n + 1 represents a set
    edges = []
    for s in IntegerVectors(n, k=2*n + 1, max_part=1):
        s1 = list(s)
        for i in range(2*n + 1):
            sig_check()
            if s1[i] == 0:
                s2 = list(s)  # duplicate list
                s2[i] = 1
                edges.append((tuple(s1), tuple(s2)))

    G = Graph(edges, format='list_of_edges')
    G.name("Bipartite double of Odd graph on a set of %d elements" % (2*n + 1))
    return G


def HalfCube(const int n):
    r"""
    Return the halved cube in `n` dimensions.

    The graph is distance-regular with classical parameters
    `(\lfloor \frac n 2 \rfloor, 1, 2, 2 \lceil \frac n 2 \rceil -1)`.

    INPUT:

    - ``n`` -- integer; must be greater than 2

    EXAMPLES::

        sage: G = graphs.HalfCube(8)
        sage: G.is_distance_regular(True)
        ([28, 15, 6, 1, None], [None, 1, 6, 15, 28])
        sage: G = graphs.HalfCube(4)
        sage: G.is_distance_regular(True)
        ([6, 1, None], [None, 1, 6])

    REFERENCES:

    See [BCN1989]_ pp. 264, 265 or [VDKT2016]_ p. 21.
    This construction can be found on
    :wikipedia:`Halved_cube_graph#Equivalent_constructions`

    TESTS:

    HalfCube is a half of the CubeGraph::

         sage: H = graphs.CubeGraph(8)
         sage: s1, s2 = H.bipartite_sets()
         sage: G1 = Graph([s1, lambda i, j: H.distance(i, j) == 2])
         sage: G2 = graphs.HalfCube(8)
         sage: G1.is_isomorphic(G2)
         True
    """
    from math import cos, sin, pi

    if n < 2:
        raise ValueError("the dimension must be n > 1")

    cdef int u, uu, v, i, j
    cdef list E = []
    cdef dict pos = {}  # dictionary of positions
    cdef float theta = pi / (n - 1)
    cdef list cosi = [<float>cos(i*theta) for i in range(n - 1)]
    cdef list sini = [<float>sin(i*theta) for i in range(n - 1)]

    for u in range(2**(n - 1)):
        sig_check()
        pos[u] = (sum(((u >> (n-2-i)) & 1) * cosi[i] for i in range(n - 1)),
                  sum(((u >> (n-2-i)) & 1) * sini[i] for i in range(n - 1)))

        for i in range(n - 1):
            uu = u ^ (1 << i)
            if u < uu:
                E.append((u, uu))
            for j in range(i + 1, n - 1):
                v = uu ^ (1 << j)
                if u < v:
                    E.append((u, v))

    G = Graph([range(2**(n - 1)), E], format='vertices_and_edges')
    G.set_pos(pos)
    G.name("Half %d Cube" % n)
    return G


def GrassmannGraph(const int q, const int n, const int input_e):
    r"""
    Return the Grassmann graph with parameters `(q, n, e)`.

    This builds the Grassmann graph `J_q(n,e)`. That is, for a vector
    space `V = \mathbb F(q)^n` the output is the graph on the subspaces
    of dimension `e` where two subspaces are adjacent if their intersection
    has dimension `e-1`.

    This graph is distance-regular with classical parameters
    `(\min(e, n-e), q, q, \genfrac {[}{]} {0pt} {} {n-e+1} 1 _q -1)`

    INPUT:

    - ``q`` -- a prime power
    - ``n, e`` -- integers with ``n > e+1``

    EXAMPLES::

        sage: G = graphs.GrassmannGraph(2, 4, 2)                                        # needs sage.modules sage.rings.finite_rings
        sage: G.is_distance_regular(True)                                               # needs sage.modules sage.rings.finite_rings
        ([18, 8, None], [None, 1, 9])

    REFERENCES:

    See [BCN1989]_ pp. 268-272 or [VDKT2016]_ p. 21.

    TESTS::

        sage: # needs sage.modules sage.rings.finite_rings
        sage: G = graphs.GrassmannGraph(2, 6, 3)        # long time
        sage: G.is_distance_regular(True)       # long time
        ([98, 72, 32, None], [None, 1, 9, 49])
        sage: G = graphs.GrassmannGraph(3, 4, 2)
        sage: G.is_distance_regular(True)
        ([48, 27, None], [None, 1, 16])
    """
    from sage.combinat.designs import design_catalog as designs

    if n <= input_e + 1:
        raise ValueError(f"Impossible parameters n <= e+1 ({n} > {input_e + 1})")

    e = input_e
    if n < 2 * input_e:
        e = n - input_e

    PG = designs.ProjectiveGeometryDesign(n - 1, e - 1, q)
    # we want the intersection graph
    # the size of the intersection must be (q^{e-1} - 1) / (q-1)
    size = (q**(e - 1) - 1) // (q - 1)
    G = PG.intersection_graph([size])
    G.name("Grassmann graph J_%d(%d, %d)" % (q, n, e))
    return G


def DoubleGrassmannGraph(const int q, const int e):
    r"""
    Return the bipartite double of the distance-`e` graph of the Grassmann graph `J_q(n,e)`.

    This graph can also be descirbed as follows:
    Let `V` be the vector space of dimension `n` over `GF(q)`.
    The vertex set is the set of `e+1` or `e` subspaces of `V`.
    Two vertices are adjacent if one subspace is contained in the other.

    This graph is distance-transitive.

    INPUT:

    - ``q`` -- a prime power
    - ``e`` -- integer

    EXAMPLES::

        sage: G = graphs.DoubleGrassmannGraph(2,1)                                      # needs sage.modules
        sage: G.diameter()                                                              # needs sage.modules
        3
        sage: G.is_distance_regular(True)                                               # needs sage.modules
        ([3, 2, 2, None], [None, 1, 1, 3])


    REFERENCES:

    See [BCN1989]_ pp. 272, 273 or [VDKT2016]_ p. 25.

    TESTS::

         sage: # needs sage.modules
         sage: G = graphs.DoubleGrassmannGraph(5,1)
         sage: G.order()
         62
         sage: G.is_distance_regular(True)
         ([6, 5, 5, None], [None, 1, 1, 6])
         sage: G = graphs.DoubleGrassmannGraph(3, 2)    # long time                     # needs sage.rings.finite_rings
         sage: G.order()                        # long time                             # needs sage.rings.finite_rings
         2420
         sage: G.is_distance_regular(True)      # long time                             # needs sage.rings.finite_rings
         ([13, 12, 12, 9, 9, None], [None, 1, 1, 4, 4, 13])
    """
    n = 2*e + 1
    V = VectorSpace(GF(q), n)

    edges = []
    for W in V.subspaces(e + 1):
        Wbasis = frozenset(W.basis())
        for U in W.subspaces(e):
            sig_check()
            Ubasis = frozenset(U.basis())
            edges.append((Wbasis, Ubasis))

    G = Graph(edges, format='list_of_edges')
    G.name("Double Grassmann graph (%d, %d, %d)" % (n, e, q))
    return G


def is_from_GQ_spread(list arr):
    r"""
    Return a pair `(s, t)` if the graph obtained from a GQ of order `(s, t)`
    with a spread has the intersection array passed. We also require that such
    GQ can be built by Sage.
    If no such pair exists, then return ``False``.

    INPUT:

    - ``arr`` -- list; an intersection array

    EXAMPLES::

         sage: from sage.graphs.generators.distance_regular import \
         ....: is_from_GQ_spread, graph_from_GQ_spread
         sage: is_from_GQ_spread([125, 120, 1, 1, 24, 125])                             # needs sage.libs.pari
         (5, 25)
         sage: G = graph_from_GQ_spread(5, 25)                                          # needs sage.libs.pari
         sage: G.is_distance_regular(True)                                              # needs sage.libs.pari
         ([125, 120, 1, None], [None, 1, 24, 125])

    REFERENCES:

    The graphs we are looking for are antipodal covers of complete graphs.
    See [BCN1989]_ pp. 385, 386 for a discussion on these particular case.

    TESTS::

         sage: from sage.graphs.generators.distance_regular import \
         ....: is_from_GQ_spread
         sage: is_from_GQ_spread([343, 336, 1, 1, 48, 343])                             # needs sage.libs.pari
         (7, 49)
         sage: is_from_GQ_spread([343, 336, 1, 2, 48, 343])                             # needs sage.libs.pari
         False

    Check that we don't get ``True`` for inexisting GQs::

         sage: from sage.graphs.generators.distance_regular import \
         ....: is_from_GQ_spread
         sage: s = 5
         sage: t = 6
         sage: [s * t, s * (t-1), 1, 1, t - 1, s * t]
         [30, 25, 1, 1, 5, 30]
         sage: is_from_GQ_spread([30, 25, 1, 1, 5, 30])                                 # needs sage.libs.pari
         False
    """
    from sage.combinat.designs import design_catalog as designs

    if len(arr) != 6:
        return False

    t = arr[4] + 1
    if t <= 1:  # avoid division by 0
        return False

    s = arr[1] // (t-1)
    if s == 1 and t == 1:  # in this case we don't get a connected graph
        return False

    if arr != [s * t, s * (t-1), 1, 1, t - 1, s * t]:
        return False

    # check Sage can build it (it may not exist)
    if designs.generalised_quadrangle_with_spread(s, t, existence=True) \
       is not True:
        return False

    return (s, t)


def graph_from_GQ_spread(const int s, const int t):
    r"""
    Return the point graph of the generalised quadrangle with
    order `(s, t)` after removing one of its spreads.

    These graphs are antipodal covers of complete graphs and, in particular,
    distance-regular graphs of diameter 3.

    INPUT:

    - ``s, t`` -- integers; order of the generalised quadrangle

    EXAMPLES::

         sage: from sage.graphs.generators.distance_regular import \
         ....: graph_from_GQ_spread
         sage: G = graph_from_GQ_spread(4, 16)                                          # needs sage.libs.pari
         sage: G.is_distance_regular(True)                                              # needs sage.libs.pari
         ([64, 60, 1, None], [None, 1, 15, 64])

    REFERENCES:

    The graphs constructed here follow [BCN1989]_ pp. 385, 386.

    TESTS::

         sage: from sage.graphs.generators.distance_regular import \
         ....: graph_from_GQ_spread, is_from_GQ_spread
         sage: is_from_GQ_spread([64, 60, 1, 1, 15, 64])                                # needs sage.libs.pari
         (4, 16)
         sage: graph_from_GQ_spread(*is_from_GQ_spread([27, 24, 1, 1, 8, 27]))          # needs sage.libs.pari
         Graph on 112 vertices
         sage: _.is_distance_regular(True)                                              # needs sage.libs.pari
         ([27, 24, 1, None], [None, 1, 8, 27])
    """
    from sage.combinat.designs import design_catalog as designs

    (GQ, S) = designs.generalised_quadrangle_with_spread(s, t, check=False)

    edges = []
    for b in GQ.blocks():
        if b in S:  # skip blocks in spread
            continue
        for p1, p2 in itertools.combinations(b, 2):
            sig_check()
            edges.append((p1, p2))

    return Graph(edges, format="list_of_edges")


def GeneralisedDodecagonGraph(const int s, const int t):
    r"""
    Return the point-graph of a generalised dodecagon of order `(s,t)`.

    INPUT:

    - ``s, t`` -- integers; order of the generalised dodecagon

    EXAMPLES::

        sage: # optional - gap_packages internet
        sage: G = graphs.GeneralisedDodecagonGraph(1, 5)
        sage: G.is_distance_regular(True)
        ([6, 5, 5, 5, 5, 5, None], [None, 1, 1, 1, 1, 1, 6])
        sage: H = graphs.GeneralisedDodecagonGraph(5, 1)
        sage: H.order()
        23436
        sage: H.is_distance_regular(True)       # not tested (6 min)
        ([10, 5, 5, 5, 5, 5, None], [None, 1, 1, 1, 1, 1, 2])

    .. NOTE::

        This function indirectly uses the GAP's AtlasRep package.
        Thus you may need an internet connection and the optional Sage's
        package ``gap_packages``.

    REFERENCES:

    See [BCN1989]_ pp. 200-205 for a discussion of distance-regular graphs from
    generalised polygons.

    TESTS:

    Test all graphs of order `(1, q)`::

        sage: # optional - gap_packages internet
        sage: G = graphs.GeneralisedDodecagonGraph(1, 4)
        sage: G.is_distance_regular(True)
        ([5, 4, 4, 4, 4, 4, None], [None, 1, 1, 1, 1, 1, 5])
        sage: G = graphs.GeneralisedDodecagonGraph(1, 3)
        sage: G.is_distance_regular(True)
        ([4, 3, 3, 3, 3, 3, None], [None, 1, 1, 1, 1, 1, 4])
        sage: G = graphs.GeneralisedDodecagonGraph(1, 2)
        sage: G.is_distance_regular(True)
        ([3, 2, 2, 2, 2, 2, None], [None, 1, 1, 1, 1, 1, 3])
        sage: G = graphs.GeneralisedDodecagonGraph(1, 1)
        sage: G.is_distance_regular(True)
        ([2, 1, 1, 1, 1, 1, None], [None, 1, 1, 1, 1, 1, 2])

    Now test all graphs of order `(q, 1)`::

        sage: # optional - gap_packages internet
        sage: G = graphs.GeneralisedDodecagonGraph(4, 1)
        sage: G.is_distance_regular(True)
        ([8, 4, 4, 4, 4, 4, None], [None, 1, 1, 1, 1, 1, 2])
        sage: G = graphs.GeneralisedDodecagonGraph(3, 1)
        sage: G.is_distance_regular(True)
        ([6, 3, 3, 3, 3, 3, None], [None, 1, 1, 1, 1, 1, 2])
        sage: G = graphs.GeneralisedDodecagonGraph(2, 1)
        sage: G.is_distance_regular(True)
        ([4, 2, 2, 2, 2, 2, None], [None, 1, 1, 1, 1, 1, 2])
    """
    from sage.arith.misc import is_prime_power

    cdef int q = 0
    cdef int orderType = 0

    # decide the type of graph
    if s == 1:  # (1, q)
        q = t
    elif t == 1:  # (q, 1)
        q = s
        orderType = 1
    else:
        raise ValueError(
            f"Generalised dodecagon of order ({s}, {t}) does not exist")

    if q == 1:  # order (1, 1)
        from sage.graphs.generators.basic import CycleGraph
        return CycleGraph(12)

    if not is_prime_power(q):
        raise ValueError(
            f"No generalised dodecagon of order ({s}, {t}) is known")

    if orderType == 0:
        # incidence graph of hexagon (q,q)
        H = GeneralisedHexagonGraph(q, q)
        lines = _extract_lines(H)

        edges = []
        for l in lines:
            for p in l:
                sig_check()
                edges.append((p, l))

        G = Graph(edges, format='list_of_edges')
        G.name("Generalised dodecagon of order (1, %d)" % q)
        return G

    else:  # orderType == 1
        # dual
        H = GeneralisedDodecagonGraph(t, s)
        G = _line_graph_generalised_polygon(H)
        G.name("Generalised dodecagon of order (%s, %d)" % (s, t))
        return G


def GeneralisedOctagonGraph(const int s, const int t):
    r"""
    Return the point-graph of a generalised octagon of order `(s,t)`.

    INPUT:

    - ``s, t`` -- integers; order of the generalised octagon

    EXAMPLES::

        sage: # needs sage.libs.gap
        sage: G = graphs.GeneralisedOctagonGraph(1, 4)
        sage: G.is_distance_regular(True)
        ([5, 4, 4, 4, None], [None, 1, 1, 1, 5])
        sage: G = graphs.GeneralisedOctagonGraph(2, 4)  # optional - gap_packages internet
        sage: G.is_distance_regular(True)  # optional - gap_packages internet
        ([10, 8, 8, 8, None], [None, 1, 1, 1, 5])
        sage: G = graphs.GeneralisedOctagonGraph(5, 1)
        sage: G.is_distance_regular(True)
        ([10, 5, 5, 5, None], [None, 1, 1, 1, 2])

    .. NOTE::

        This function uses the GAP's AtlasRep package to build the graphs
        of order `(2, 4)` or `(4, 2)`. For those graphs you need an internet
        connection and Sage's optional package ``gap_packages``.

    REFERENCES:

    See [BCN1989]_ pp. 200-205 for a discussion of distance-regular graphs from
    generalised polygons.

    TESTS::

        sage: G = graphs.GeneralisedOctagonGraph(8, 64)                                 # needs sage.libs.gap
        Traceback (most recent call last):
        ...
        NotImplementedError: Graph would be too big
        sage: G = graphs.GeneralisedOctagonGraph(4, 16)                                 # needs sage.libs.gap
        Traceback (most recent call last):
        ...
        ValueError: generalised octagons of order (q, q^2) are known only for odd powers q of 2
    """
    from sage.arith.misc import is_prime_power
    from sage.libs.gap.libgap import libgap
    from sage.graphs.strongly_regular_db import strongly_regular_graph

    cdef int q = 0
    cdef int orderType = 0

    if s == 1:  # (1, q)
        q = t
    elif t == 1:  # (q, 1)
        q = s
        orderType = 1
    elif s**2 == t:  # (q, q^2)
        q = s
        (p, k) = is_prime_power(q, get_data=True)

        if p != 2 or k % 2 != 1:
            raise ValueError(("generalised octagons of order (q, q^2) "
                              "are known only for odd powers q of 2"))
        orderType = 2
    elif t**2 == s:  # (q^2, q)
        q = t
        orderType = 1
    else:
        raise ValueError(f"No generalised octagon of order ({s}, {t}) is known")

    if q == 1:  # order (1, 1)
        from sage.graphs.generators.basic import CycleGraph
        return CycleGraph(8)

    if not is_prime_power(q):
        raise ValueError(f"No generalised octagon of order ({s}, {t}) is known")

    if orderType == 0:
        # incidence graph of generalised quadrangle (q, q)

        H = strongly_regular_graph((q+1) * (q*q + 1), q * (q+1), q - 1, q + 1,
                                   check=False)

        lines = _extract_lines(H)

        edges = []
        for l in lines:
            for p in l:
                sig_check()
                edges.append((p, l))

        G = Graph(edges, format='list_of_edges')
        G.name("Generalised octagon of order (1, %d)" % q)
        return G

    elif orderType == 1:
        # dual
        H = GeneralisedOctagonGraph(t, s)
        G = _line_graph_generalised_polygon(H)
        G.name("Generalised octagon of order(%d, %d)" % (s, t))
        return G
    else:
        if q == 2:
            group = libgap.AtlasGroup("2F4(2)", libgap.NrMovedPoints, 1755)
            G = Graph(libgap.Orbit(group, [1, 73], libgap.OnSets),
                      format='list_of_edges')
            G.name("Generalised octagon of order (2, 4)")
            return G
        else:
            raise NotImplementedError("Graph would be too big")


def GeneralisedHexagonGraph(const int s, const int t):
    r"""
    Return the point-graph of a generalised hexagon of order `(s,t)`.

    INPUT:

    - ``s, t`` -- integers; order of the generalised hexagon

    EXAMPLES::

        sage: # needs sage.libs.gap
        sage: G = graphs.GeneralisedHexagonGraph(5, 5)  # optional - gap_packages internet
        sage: G.is_distance_regular(True)  # optional - gap_packages internet
        ([30, 25, 25, None], [None, 1, 1, 6])
        sage: G = graphs.GeneralisedHexagonGraph(7, 1)
        sage: G.is_distance_regular(True)
        ([14, 7, 7, None], [None, 1, 1, 2])
        sage: graphs.GeneralisedHexagonGraph(1, 1)
        Cycle graph: Graph on 6 vertices

    .. NOTE::

        This function uses the GAP's AtlasRep package to build GHs
        of order `(q, q)`, `(q, q^3)` or `(q^3, q)`. For those graphs you need
        an internet connection and Sage's optional package ``gap_packages``.

    REFERENCES:

    See [BCN1989]_ pp. 200-205 for a discussion of distance-regular graphs from
    generalised polygons.

    TESTS::

        sage: # optional - gap_packages internet
        sage: G = graphs.GeneralisedHexagonGraph(4, 4)
        sage: G.is_distance_regular(True)
        ([20, 16, 16, None], [None, 1, 1, 5])
        sage: G = graphs.GeneralisedHexagonGraph(3, 3)
        sage: G.is_distance_regular(True)
        ([12, 9, 9, None], [None, 1, 1, 4])
        sage: G = graphs.GeneralisedHexagonGraph(2, 2)
        sage: G.is_distance_regular(True)
        ([6, 4, 4, None], [None, 1, 1, 3])
        sage: G = graphs.GeneralisedHexagonGraph(2, 8)
        sage: G.is_distance_regular(True)
        ([18, 16, 16, None], [None, 1, 1, 9])
    """
    from sage.arith.misc import is_prime_power
    from sage.libs.gap.libgap import libgap
    from sage.combinat.designs import design_catalog as designs

    cdef int q = 0
    cdef int orderType = 0

    if s == 1:  # (1, q)
        q = t
    elif t == 1:  # (q, 1)
        q = s
        orderType = 1
    elif s == t:  # (q, q)
        q = s
        orderType = 2
    elif s**3 == t:  # (q, q^3)
        q = s
        orderType = 3
    elif t**3 == s:  # (q^3, q)
        q = t
        orderType = 1
    else:
        raise ValueError(f"No generalised hexagon of order ({s}, {t}) is known")

    if q == 1:  # order (1, 1)
        from sage.graphs.generators.basic import CycleGraph
        return CycleGraph(6)

    if not is_prime_power(q):
        raise ValueError(f"No generalised hexagon of order ({s}, {t}) is known")

    if orderType == 0:
        # incidence graph of generalised 3-gon of order (q, q)
        PG2 = designs.ProjectiveGeometryDesign(2, 1, q)

        edges = []
        for l in PG2.blocks():
            for p in l:
                sig_check()
                edges.append((p, tuple(l)))

        G = Graph(edges, format='list_of_edges')
        G.name("Generalised hexagon of order (1, %d)" % q)
        return G

    elif orderType == 1:
        # dual graph
        H = GeneralisedHexagonGraph(t, s)
        G = _line_graph_generalised_polygon(H)
        G.name("Generalised hexagon of order(%d, %d)" % (s, t))
        return G

    elif orderType == 2:
        # we use the group G2(q)
        # if q == 2, then G2(2) is isomorphic to U3(3).2
        if q == 2:
            group = libgap.AtlasGroup("U3(3).2", libgap.NrMovedPoints, 63)
            G = Graph(libgap.Orbit(group, [1, 19], libgap.OnSets),
                      format='list_of_edges')
            G.name("Generalised hexagon of order (%d, %d)" % (q, q))
            return G

        elif q == 3:  # we don't have permutation representation; so we build it
            matrixRep = libgap.AtlasGroup("G2(3)", libgap.Position, 7)
            e1 = vector(GF(3), [1, 0, 0, 0, 0, 0, 0])
            orb = libgap.Orbit(matrixRep, e1, libgap.OnLines)
            group = libgap.Action(matrixRep, orb, libgap.OnLines)

            # now group is our permutation representation
            G = Graph(libgap.Orbit(group, [1, 52], libgap.OnSets),
                      format='list_of_edges')
            G.name("Generalised hexagon of order (%d, %d)" % (q, q))
            return G

        elif q <= 5:
            n = 1365 if q == 4 else 3906
            p = 43 if q == 4 else 185
            group = libgap.AtlasGroup("G2(%d)" % q, libgap.NrMovedPoints, n)

            G = Graph(libgap.Orbit(group, [1, p], libgap.OnSets),
                      format='list_of_edges')
            G.name("Generalised hexagon of order (%d, %d)" % (q, q))
            return G

        else:
            raise NotImplementedError("Graph would be too big")

    elif orderType == 3:
        if q > 3:
            raise NotImplementedError("Graph would be too big")

        movedPoints = 819 if q == 2 else 26572
        group = libgap.AtlasGroup("3D4(%d)" % q, libgap.NrMovedPoints, movedPoints)

        G = Graph(libgap.Orbit(group, [1, 2], libgap.OnSets),
                  format='list_of_edges')
        G.name("Generalised hexagon of order (%d, %d)" % (q, q**3))
        return G


def _extract_lines(G):
    r"""
    Return the set of lines from the point-graph of a generalised polygon.

    In particular, given a graph `G` we return the set of singular lines:
    Let `(x,y)` be an edge, then `\{x,y\}^\bot^\bot` is a singular line.
    We define `x^\bot =` neighbours of `x` and `x` for `x` a vertex and
    `S^\bot =` intersection of `x^\bot` for all `x \in S`.

    INPUT:

    - ``G`` -- a graph

    OUTPUT:

    A list of tuples where each tuple represent a line via the vertices
    contained in that line.

    EXAMPLES::

        sage: # needs sage.libs.gap
        sage: from sage.graphs.generators.distance_regular import _extract_lines
        sage: G = graphs.GeneralisedHexagonGraph(1, 8)
        sage: lines = _extract_lines(G)
        sage: len(lines)
        657
        sage: type(lines)
        <class 'list'>
        sage: line = lines[0]
        sage: type(line)
        <class 'tuple'>
        sage: line[0] in G  # the elements in the line are vertices
        True

    REFERENCES:

    See [BCN1989]_ pp. 200-205 for a discussion of distance-regular graphs from
    generalised polygons. See also [BCN1989]_ pp. 28, 29 for some theory about
    singular lines.
    """
    lines = []
    edges = set(G.edges(labels=False, sort=False))

    while edges:
        (x, y) = edges.pop()

        # compute line
        botX = set(G.neighbors(x, closed=True))
        botY = set(G.neighbors(y, closed=True))
        bot1 = botX.intersection(botY)

        b = bot1.pop()
        bot2 = frozenset(G.neighbors(b, closed=True))
        for v in bot1:
            sig_check()
            s = frozenset(G.neighbors(v, closed=True))
            bot2 = bot2.intersection(s)

        # now bot2 is a line
        lines.append(tuple(bot2))  # we need tuple or GAP will complain later

        # remove already handled edges
        for u, v in itertools.product(bot2, repeat=2):
            sig_check()
            try:
                edges.remove((u, v))
            except KeyError:
                pass  # ignore this
    # end while edges

    return lines


def _line_graph_generalised_polygon(H):
    r"""
    Return the line-graph of the generalised polygon whose point-graph is `H`.

    In particular, return the line-graph of the incidence structure defined
    by the singular lines of the graph `H`.

    See also :func:`sage.graphs.generators.distance_regular._extract_lines`.

    INPUT:

    - ``H`` -- a graph

    EXAMPLES::

        sage: # needs sage.libs.gap
        sage: from sage.graphs.generators.distance_regular import \
        ....: _line_graph_generalised_polygon
        sage: G = graphs.GeneralisedHexagonGraph(1, 8)
        sage: H = _line_graph_generalised_polygon(G)
        sage: H.is_distance_regular(True)
        ([16, 8, 8, None], [None, 1, 1, 2])
        sage: G = graphs.GeneralisedHexagonGraph(3, 3) # optional - gap_packages internet
        sage: H = _line_graph_generalised_polygon(G)   # optional - gap_packages internet
        sage: G.is_isomorphic(H)                       # optional - gap_packages internet
        True

    REFERENCES:

    See [BCN1989]_ pp. 200-205 for a discussion of distance-regular graphs from
    generalised polygons. See also [BCN1989]_ pp. 28, 29 for some theory about
    singular lines.
    """
    lines = _extract_lines(H)

    # get a map (point -> all lines incident to point)
    vToLines = {v: [] for v in H}
    for l in lines:
        for p in l:
            sig_check()
            vToLines[p].append(l)

    edges = []
    for v in vToLines:
        lines = vToLines[v]
        for l1, l2 in itertools.combinations(lines, 2):
            sig_check()
            edges.append((l1, l2))

    return Graph(edges, format="list_of_edges")


def _intersection_array_from_graph(G):
    r"""
    Return the intersection array of the graph `G`.
    If `G` is not distance-regular, then return ``False``.

    This is a simple wrapper around
    :meth:`sage.graphs.distances_all_pairs.is_distance_regular` to return a list
    instead of a pair of lists

    INPUT:

    - G -- a graph

    EXAMPLES::

        sage: from sage.graphs.generators.distance_regular import \
        ....: _intersection_array_from_graph
        sage: _intersection_array_from_graph(graphs.FosterGraph())                      # needs networkx
        [3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3]
        sage: graphs.FosterGraph().is_distance_regular(True)                            # needs networkx
        ([3, 2, 2, 2, 2, 1, 1, 1, None], [None, 1, 1, 1, 1, 2, 2, 2, 3])
        sage: graphs.DartGraph().is_distance_regular()
        False
        sage: _intersection_array_from_graph(graphs.DartGraph())
        False

    TESTS::

        sage: from sage.graphs.generators.distance_regular import \
        ....: _intersection_array_from_graph
        sage: _intersection_array_from_graph(Graph())
        []
        sage: _intersection_array_from_graph(Graph(3))
        []
        sage: _intersection_array_from_graph(graphs.CompleteGraph(7))
        [6, 1]
    """
    t = G.is_distance_regular(True)
    if t is False:
        return False

    return t[0][:-1] + t[1][1:]


cdef enum ClassicalParametersGraph:
    NonExisting = 0,
    Johnson,
    Hamming,
    HalvedCube,
    UnitaryDualPolar,
    HermitianForms,
    GeneralisedHexagon,
    Grassmann,
    OrthogonalDualPolar1,
    SymplecticDualPolar,
    OrthogonalDualPolar2,
    UnitaryDualPolar1,
    UnitaryDualPolar2,
    Ustimenko,
    BilinearForms,
    AlternatingForms,
    LieE77,
    AffineE6


def is_classical_parameters_graph(list array):
    r"""
    Return a tuple of parameters representing the array given. If such no tuple
    can be produced, it returns ``False``.

    Given an intersection array, if it represents a family of  distance-regular
    graphs with classical parameters, then this function  returns a tuple
    consisting of the  parameters `(d, b, \alpha, \beta)` and a fourth parameter
    which is the enum ``CalssicalParametersGraph`` indicating the family with
    the given itersection array.
    If the array doesn't belong to any classical parameter graph, then this
    function returns ``False``.
    If the array belongs to a sporadic graph rather than a family of graphs,
    then the function returns ``False``. This is to reduce the overlap with
    ``sage.graphs.generators.distance_regular._sporadic_graph_database``.


    .. NOTE::

        The array given as an input is expected to be an intersection array.
        If this is not the case, then some exception may be raised.

    INPUT:

    - ``array`` -- list; an intersection array

    OUTPUT:

    ``False`` or a tuple ``(d, b, alpha, beta, gamma)``.

    EXAMPLES::

        sage: from sage.graphs.generators.distance_regular import \
        ....: is_classical_parameters_graph
        sage: G = graphs.HammingGraph(5, 4)
        sage: G.is_distance_regular(True)
        ([15, 12, 9, 6, 3, None], [None, 1, 2, 3, 4, 5])
        sage: is_classical_parameters_graph([15, 12, 9, 6, 3, 1, 2, 3, 4, 5])           # needs sage.combinat
        (5, 1, 0, 3, 2)

    REFERENCES:

    See [BCN1989]_ chapter 9 for a discussion of distance-regular graphs with
    classical parameters. See [BCN1989]_ chapter 6.2 for a method to compute
    the classical parameters of a graph. See also [VDKT2016]_ section 3.1.1.

    TESTS::

        sage: from sage.graphs.generators.distance_regular import \
        ....: is_classical_parameters_graph
        sage: is_classical_parameters_graph([68, 64, 1, 17])  # srg not drg             # needs sage.combinat
        False
        sage: G = graphs.GossetGraph() # sporadic classical parameters graph
        sage: G.is_distance_regular(True)
        ([27, 10, 1, None], [None, 1, 10, 27])
        sage: is_classical_parameters_graph([27, 10, 1, 1, 10, 27])                     # needs sage.combinat
        False
    """
    from sage.misc.functional import log
    from sage.rings.integer_ring import ZZ
    from sage.arith.misc import is_prime_power
    from sage.combinat.q_analogues import q_binomial

    def integral_log(const int x, const int b):
        # compute log_b(x) if is not a positive integer, return -1
        if x <= 0:
            return -1
        k = log(x, b)
        if k in ZZ and k > 0:
            return int(k)
        return -1

    def check_parameters(int d, int b, int alpha, int beta, list arr):
        bs = [(q_binomial(d, 1, b) - q_binomial(i, 1, b)) *
              (beta - alpha * q_binomial(i, 1, b)) for i in range(d)]
        cs = [q_binomial(i, 1, b) * (1 + alpha*q_binomial(i-1, 1, b))
              for i in range(1, d+1)]
        return arr == bs + cs

    def b_(i):
        if i == d:
            return 0
        return array[i]

    def c_(i):
        if i == 0:
            return 0
        return array[d - 1 + i]

    def a_(i):
        return b_(0) - b_(i) - c_(i)

    if len(array) % 2:
        return False

    d = len(array) // 2
    if d < 3:
        return False

    # if array is classical parameters, then we have 2 cases:
    # 1. a_i != \lambda c_i for 2<= i <= d
    # 2. a_i == \lambda c_i for 0<= i <= d
    if all(a_(i) == a_(1) * c_(i) for i in range(2, d+1)):
        case = 2
    elif all(a_(i) != a_(1) * c_(i) for i in range(2, d+1)):
        case = 1
    else:
        return False

    if case == 1:
        # b = (a_2 c_3 - c_2 a_3) / (a_1 c_3 - a_3)
        num = a_(2) * c_(3) - c_(2) * a_(3)
        den = a_(1) * c_(3) - a_(3)
        b = num / den
        if b in {0, -1}:
            return False

        alpha = c_(2) / (1 + b) - 1
        beta = b_(0) / q_binomial(d, 1, b)

    else:  # case 2
        # b is either c_2 - 1 or - a_1 - 1
        # try c_2 - 1
        valid = True
        b = c_(2) - 1
        if b in {0, -1}:
            valid = False
        else:
            alpha = c_(2) / (1 + b) - 1
            beta = b_(0) / q_binomial(d, 1, b)
            valid = check_parameters(d, b, alpha, beta, array)

        if not valid:  # must have -a_1 - 1
            b = -a_(1) - 1
            if b in {0, -1}:
                return False  # even this case is invalid

            alpha = c_(2) / (1 + b) - 1
            beta = a_(1) + 1 - alpha*(q_binomial(d, 1, b) - 1)

    if not check_parameters(d, b, alpha, beta, array):
        return False

    gamma = ClassicalParametersGraph.NonExisting

    if b == 1:
        if alpha == 1 and beta >= d:  # since beta+d = n >= 2*d
            # Johnson Graph
            gamma = ClassicalParametersGraph.Johnson
        elif alpha == 0:
            # Hamming Graph
            gamma = ClassicalParametersGraph.Hamming
        elif alpha == 2 and (beta == 2*d + 1 or beta == 2*d - 1):
            # Halved cube graph
            gamma = ClassicalParametersGraph.HalvedCube
        else:
            return False  # no other (unbounbded) drg exists with b = 1

    elif b < 0 and is_prime_power(-b):
        if (alpha + 1 == (1 + b*b) / (1 + b) and
                beta + 1 == (1 - b**(d + 1)) / (1 + b)):
            # U(2d,r)
            gamma = ClassicalParametersGraph.UnitaryDualPolar1
        elif alpha + 1 == b and beta + 1 == - (b**d):
            gamma = ClassicalParametersGraph.HermitianForms
        elif (d == 3 and alpha + 1 == 1 / (1+b) and
              beta + 1 == q_binomial(3, 1, -b)):
            gamma = ClassicalParametersGraph.GeneralisedHexagon
        else:
            return False

    if gamma != ClassicalParametersGraph.NonExisting:
        return (d, b, alpha, beta, gamma)

    # all remaining cases need b to be a prime power
    (p, k) = is_prime_power(b, get_data=True)
    if k == 0:  # b not a prime power
        return False
    r = p**(k//2)  # will be used later

    if alpha == b and integral_log((beta+1) * (b-1) + 1, b) >= d + 1:
        # we checked that beta + 1 = (b^(n-d+1) - 1)/(b - 1) for n >= 2d
        # Grassmann graph
        gamma = ClassicalParametersGraph.Grassmann

    elif alpha == 0 and beta * beta in {1, b, b * b, b**3, b**4}:
        # checked beta in {b^0, b^(0.5), b, b^(1.5), b^2}
        # dual polar graphs
        if beta == 1:
            gamma = ClassicalParametersGraph.OrthogonalDualPolar1
        elif beta == b:
            gamma = ClassicalParametersGraph.SymplecticDualPolar
        elif beta == b * b:
            gamma = ClassicalParametersGraph.OrthogonalDualPolar2
        else:
            if k % 2 == 1:
                return False  # we need b a square
            if beta == r**3:
                gamma = ClassicalParametersGraph.UnitaryDualPolar1
            elif beta == r:
                gamma = ClassicalParametersGraph.UnitaryDualPolar2

    elif (k % 2 == 0 and alpha + 1 == q_binomial(3, 1, r) and
          beta + 1 in {q_binomial(2*d + 2, 1, r),
                       q_binomial(2*d + 1, 1, r)}):
        gamma = ClassicalParametersGraph.Ustimenko

    elif alpha + 1 == b and integral_log(beta + 1, b) >= d:
        gamma = ClassicalParametersGraph.BilinearForms

    elif (k % 2 == 0 and alpha + 1 == b and
          beta + 1 in {r**(2*d - 1), r**(2*d + 1)}):
        gamma = ClassicalParametersGraph.AlternatingForms

    elif (d == 3 and k % 4 == 0 and alpha + 1 == q_binomial(5, 1, p**(k//4)) and
          beta + 1 == q_binomial(10, 1, p**(k//4))):
        gamma = ClassicalParametersGraph.LieE77

    elif d == 3 and k % 4 == 0 and alpha + 1 == b and beta + 1 == (p**(k//4))**9:
        gamma = ClassicalParametersGraph.AffineE6

    if gamma == ClassicalParametersGraph.NonExisting:
        return False
    return (d, b, alpha, beta, gamma)


def graph_with_classical_parameters(int d, int b, alpha_in, beta_in, int gamma):
    r"""
    Return the graph with the classical parameters given.

    The last parameter ``gamma`` is meant to be an element of the enum
    ``ClassicalParametersGraph`` used to identify the family of graphs to
    construct.
    In particular this function doesn't build any sporadic graph.
    To build such a graph use
    :func:`sage.graphs.generators.distance_regular.distance_regular_graph`.

    INPUT:

    - ``d, b, alpha_in, beta_in`` -- numbers; the parameters of the graph;
      ``d`` and ``b`` must be integers

    - ``gamma`` -- element of the enum ``ClassicalParametersGraph``

    EXAMPLES::

        sage: from sage.graphs.generators.distance_regular import *
        sage: graph_with_classical_parameters(3, 1, 1, 3, 1)
        Johnson graph with parameters 6,3: Graph on 20 vertices

    The last parameter is very important as it takes precedence.
    This function will not check that the other four parameters match the correct
    family. Use
    :func:`sage.graphs.generators.distance_regular.is_classical_parameters_graph`
    to check the parameters::

        sage: from sage.graphs.generators.distance_regular import *
        sage: graph_with_classical_parameters(3, 1, 1, 3, 2)
        Hamming Graph with parameters 3,4: Graph on 64 vertices
        sage: G = _; G.is_distance_regular(True)
        ([9, 6, 3, None], [None, 1, 2, 3])
        sage: is_classical_parameters_graph([9, 6, 3, 1, 2, 3])                         # needs sage.combinat
        (3, 1, 0, 3, 2)

    Two families of graphs are not implemented yet::

        sage: from sage.graphs.generators.distance_regular import *
        sage: graph_with_classical_parameters(3, 16, 15, 511, 17)
        Traceback (most recent call last):
        ...
        NotImplementedError: Graph would be too big
        sage: graph_with_classical_parameters(3, 16, 30, 1022, 16)
        Traceback (most recent call last):
        ...
        NotImplementedError: Graph would be too big

    REFERENCES:

    See [BCN1989]_ chapter 9 for a discussion of distance-regular graphs with
    classical parameters. See also [VDKT2016]_ section 3.1.1.

    TESTS::

        sage: graph_with_classical_parameters(3, 1, 2, 3, 3)
        Half 4 Cube: Graph on 8 vertices
        sage: graph_with_classical_parameters(3, 2, 0, 2, 9)                            # needs sage.libs.gap
        Symplectic Dual Polar Graph DSp(6, 2): Graph on 135 vertices
        sage: graph_with_classical_parameters(3, 2, 2, 14, 7)  # long time
        Grassmann graph J_2(6, 3): Graph on 1395 vertices
        sage: graph_with_classical_parameters(3, -2, -2, 6, 6) # optional - gap_packages internet
        Generalised hexagon of order (2, 8): Graph on 819 vertices
    """
    from sage.rings.rational import Rational
    from sage.misc.functional import sqrt, log
    from sage.graphs.generators.families import JohnsonGraph, HammingGraph
    from sage.graphs.generators.classical_geometries import \
        UnitaryDualPolarGraph, OrthogonalDualPolarGraph, SymplecticDualPolarGraph

    alpha = Rational(alpha_in)
    beta = Rational(beta_in)
    if alpha.is_integer():
        alpha = int(alpha)
    if beta.is_integer():
        beta = int(beta)

    if gamma == ClassicalParametersGraph.Johnson:
        return JohnsonGraph(beta + d, d)

    elif gamma == ClassicalParametersGraph.Hamming:
        return HammingGraph(d, beta + 1)

    elif gamma == ClassicalParametersGraph.HalvedCube:
        a = 0 if beta == 2*d + 1 else 1
        return HalfCube(beta + a)

    elif gamma == ClassicalParametersGraph.UnitaryDualPolar:
        return UnitaryDualPolarGraph(2 * d, -b)

    elif gamma == ClassicalParametersGraph.HermitianForms:
        return HermitianFormsGraph(d, (-b)**2)

    elif gamma == ClassicalParametersGraph.GeneralisedHexagon:
        q = -b
        return GeneralisedHexagonGraph(q, q**3)

    elif gamma == ClassicalParametersGraph.Grassmann:
        n = int(log((beta + 1) * (b - 1) + 1, b)) + d - 1
        return GrassmannGraph(b, n, d)

    elif gamma == ClassicalParametersGraph.OrthogonalDualPolar1:
        return OrthogonalDualPolarGraph(1, d, b)

    elif gamma == ClassicalParametersGraph.SymplecticDualPolar:
        return SymplecticDualPolarGraph(2 * d, b)

    elif gamma == ClassicalParametersGraph.OrthogonalDualPolar2:
        return OrthogonalDualPolarGraph(-1, d, b)

    elif gamma == ClassicalParametersGraph.UnitaryDualPolar1:
        r = int(sqrt(b))
        return UnitaryDualPolarGraph(2*d + 1, r)

    elif gamma == ClassicalParametersGraph.UnitaryDualPolar2:
        r = int(sqrt(b))
        return UnitaryDualPolarGraph(2 * d, r)

    elif gamma == ClassicalParametersGraph.Ustimenko:
        q = int(sqrt(b))
        m = int(log((beta+1) * (q-1) + 1, q)) - 1
        UstimenkoGraph(m, q)

    elif gamma == ClassicalParametersGraph.BilinearForms:
        e = int(log(beta + 1, b))
        return BilinearFormsGraph(d, e, b)

    elif gamma == ClassicalParametersGraph.AlternatingForms:
        q = int(sqrt(b))
        a = 0 if beta + 1 == q**(2*d - 1) else 1
        return AlternatingFormsGraph(2*d + a, q)

    elif (gamma == ClassicalParametersGraph.LieE77 or
          gamma == ClassicalParametersGraph.AffineE6):
        raise NotImplementedError("Graph would be too big")

    raise ValueError("Incorrect family of graphs")


def is_pseudo_partition_graph(list arr):
    r"""
    Return `(m, a)` if the intersection array given satisfies:
    `b_i = (m - i)(1 + a(m - 1 - i))` for `0 \leq i < d`
    `c_i = i(1 + a(i - 1))` for `0 \leq i < d`
    `c_d = (2d + 2 - m) d (1 + a(d - 1))` where `d` is the diameter of the graph.

    If such pair `(m, a)` doesn't exist or the diameter is less than 3, then
    this function returns ``False``.

    These graphs are called pseudo partition graphs in [BCN1989]_ chapter 6.3.

    INPUT:

    - ``arr`` -- list; intersection array

    OUTPUT:

    A pair `(m, a)` of integers or ``False`` if such pair doesn't exist.

    EXAMPLES::

        sage: from sage.graphs.generators.distance_regular import *
        sage: is_pseudo_partition_graph([36, 25, 16, 1, 4, 18])
        (6, 1)
        sage: pseudo_partition_graph(6, 1)  # long time
        Folded Johnson graph with parameters 12,6: Graph on 462 vertices
        sage: _.is_distance_regular(True)  # long time
        ([36, 25, 16, None], [None, 1, 4, 18])

    REFERENCE:

    See [BCN1989]_ pp. 197, 198 or [VDKT2016]_ pp. 38, 39.

    TESTS::

        sage: from sage.graphs.generators.distance_regular import *
        sage: is_pseudo_partition_graph([])
        False
        sage: is_pseudo_partition_graph([36, 25, 16, 1, 0, 18])
        False
        sage: is_pseudo_partition_graph([217, 156, 105, 1, 12, 33])
        (7, 5)
        sage: pseudo_partition_graph(7, 5)
        Traceback (most recent call last):
        ...
        ValueError: No known graph exists
    """
    d = len(arr)
    if d % 2 != 0:
        return False

    d = d // 2

    if d < 3:
        return False

    # c_2 = 2 (1+a)
    c2 = arr[d+1]
    if c2 % 2 != 0:
        return False
    a = c2//2 - 1

    cd = arr[2*d - 1]
    K = d * (1+a * (d-1))
    if cd % K != 0:
        return False

    gamma = cd // K
    m = 2*d + 2 - gamma

    # we must have m = 2*d or 2*d +1
    if m not in {2*d, 2*d + 1}:
        return False

    newArr = [(m-i) * (1 + a * (m-1-i)) for i in range(d)] + \
             [i * (1 + a * (i-1)) for i in range(1, d)] + \
             [(2*d + 2 - m) * d * (1 + a * (d-1))]

    if arr == newArr:
        return (m, a)

    return False


def pseudo_partition_graph(int m, int a):
    r"""
    Return a pseudo partition graph with the given parameters.

    A graph is a pseudo partition graph if it is distance-regular with
    diameter at least 3 and whose intersection numbers satisfy:
    `b_i = (m - i)(1 + a(m - 1 - i))` for `0 \leq i < d`
    `c_i = i(1 + a(i - 1))` for `0 \leq i < d`
    `c_d = (2d + 2 - m) d (1 + a(d - 1))` where `d` is the diameter of the graph.

    INPUT:

    - ``m, a`` -- integers; parameters of the graph

    EXAMPLES::

        sage: from sage.graphs.generators.distance_regular import *
        sage: pseudo_partition_graph(6, 1)  # long time
        Folded Johnson graph with parameters 12,6: Graph on 462 vertices

    Not all graphs built with this function are pseudo partition graphs as
    intended by
    :func:`sage.graphs.generators.distance_regular.is_pseudo_partition_graph`,
    since that function requires the diameter to be at least 3::

        sage: from sage.graphs.generators.distance_regular import *
        sage: pseudo_partition_graph(3, 1)
        Folded Johnson graph with parameters 6,3: Graph on 10 vertices
        sage: G=_; G.is_distance_regular(True)
        ([9, None], [None, 1])
        sage: is_pseudo_partition_graph([9, 1])
        False

    REFERENCES:

    See [BCN1989]_ pp. 197, 198 or [VDKT2016]_ pp. 38, 39 for a discussion of
    known pseudo partition graphs.

    TESTS::

        sage: from sage.graphs.generators.distance_regular import *
        sage: pseudo_partition_graph(3, 3)
        Traceback (most recent call last):
        ...
        ValueError: No known graph exists
        sage: pseudo_partition_graph(8, 0).is_distance_regular(True)
        ([8, 7, 6, 5, None], [None, 1, 2, 3, 8])
        sage: pseudo_partition_graph(6, 2).is_distance_regular(True)
        ([66, 45, 28, None], [None, 1, 6, 30])
    """
    from sage.graphs.generators.families import JohnsonGraph, FoldedCubeGraph
    from sage.graphs.bipartite_graph import BipartiteGraph

    if a == 0:
        return FoldedCubeGraph(m)
    elif a == 1:
        return JohnsonGraph(2 * m, m).folded_graph()
    elif a == 2:
        return BipartiteGraph(FoldedCubeGraph(2 * m)).project_left()

    raise ValueError("No known graph exists")


cdef enum NearPolygonGraph:
    RegularPolygon = 0,
    GeneralisedPolygon,
    OddGraph,
    DoubleOdd,
    DoubleGrassmann,
    FoldedCube,
    HammingGraph,
    DualPolarGraph


def is_near_polygon(array):
    r"""
    Return a tuple of parameters which identify the near polygon graph with
    the given intersection array. If such tuple doesn't exist, return ``False``.

    Note that ``array`` may be the intersection array of a near polygon, but if
    such graph has diameter less than 3, then this function will return
    ``False``.

    INPUT:

    - ``array`` -- list; intersection array

    OUTPUT:

    The tuple has the form ``(id, params)`` where ``id`` is a value of the
    enum `NearPolygonGraph` which identify a family of graphs and ``params``
    are all parameters needed to construct the final graph.

    EXAMPLES::

        sage: from sage.graphs.generators.distance_regular import (
        ....: is_near_polygon, near_polygon_graph)
        sage: is_near_polygon([7, 6, 6, 5, 5, 4, 1, 1, 2, 2, 3, 3])                     # needs sage.combinat
        (2, 7)
        sage: near_polygon_graph(2, 7)
        Odd Graph with parameter 7: Graph on 1716 vertices
        sage: _.is_distance_regular(True)
        ([7, 6, 6, 5, 5, 4, None], [None, 1, 1, 2, 2, 3, 3])

    REFERENCES:

    See [BCN1989]_ pp. 198-206 for some theory about near polygons as well as
    a list of known examples.

    TESTS::

        sage: # needs sage.combinat sage.libs.pari
        sage: from sage.graphs.generators.distance_regular import (
        ....: is_near_polygon, near_polygon_graph)
        sage: is_near_polygon([7, 6, 6, 4, 4, 1, 1, 3, 3, 7])
        (4, (2, 2))
        sage: near_polygon_graph(4, (2, 2))
        Double Grassmann graph (5, 2, 2): Graph on 310 vertices
        sage: near_polygon_graph(*is_near_polygon([3, 2, 2, 1, 1, 3]))                  # needs sage.rings.finite_rings
        Generalised hexagon of order (1, 2): Graph on 14 vertices
        sage: is_near_polygon([16, 12, 8, 4, 1, 2, 3, 4])
        (6, (4, 5))
        sage: is_near_polygon([])
        False
        sage: is_near_polygon([25, 16, 9, 4, 1, 1, 4, 9, 16, 25])  # JohnsonGraph
        False
    """
    from sage.arith.misc import is_prime_power
    from sage.combinat.q_analogues import q_binomial
    from sage.misc.functional import log

    if len(array) % 2 != 0:
        return False

    d = len(array) // 2

    if d < 3:
        return False

    k = array[0]
    l = k - array[1] - 1

    if l < 0:
        return False

    if (any(array[i] != k - (l + 1) * array[d - 1 + i] for i in range(1, d)) or
            k < (l+1) * array[2*d - 1]):
        return False

    # additional checks
    if k < (l + 1) * array[2*d - 1] or k % (l + 1) != 0:
        return False

    # check if it is known example
    if (k == (l + 1) * array[2*d - 1] and
            all(array[d + i] == 1 for i in range(d - 1)) and
            (l + 1 > 1 or array[2*d - 1] - 1 > 1)):  # last 2 reject regular polygons
        # generalised polygon
        s = l + 1
        t = array[2*d - 1] - 1

        if ((d == 3 and (s == 1 or t == 1) and is_prime_power(s * t)) or
            (d, s, t) in {(3, 2, 2), (3, 3, 3), (3, 4, 4), (3, 5, 5), (3, 2, 8),
                          (3, 8, 2), (3, 3, 27), (3, 27, 3), (4, 2, 4), (4, 4, 2),
                          (6, 1, 2), (6, 1, 3), (6, 1, 4), (6, 1, 5), (6, 2, 1),
                          (6, 3, 1), (6, 4, 1), (6, 5, 1)}):
            return (NearPolygonGraph.GeneralisedPolygon, (d, s, t))

        if d == 4 and (s == 1 or t == 1):
            q = s * t
            if strongly_regular_graph((q + 1) * (q*q + 1), q * (q + 1), q - 1, q + 1,
                                      existence=True):
                return (NearPolygonGraph.GeneralisedPolygon, (d, s, t))

        # otherwise not known generalised polygon
        return False

    n = 2 * d if k == (l + 1) * array[2*d - 1] else 2*d + 1

    if (k == 2 and l == 0 and all(array[d + i] == 1 for i in range(d - 1)) and
            array[2*d - 1] in {1, 2}):
        return (NearPolygonGraph.RegularPolygon, 2*d + 2 - array[2*d - 1])

    if (l == 0 and k == d + 1 and n == 2*d + 1 and
            all(array[d + i] == (i + 2) // 2 for i in range(d))):
        return (NearPolygonGraph.OddGraph, d + 1)

    if (l == 0 and k == n and all(array[d - 1 + i] == i for i in range(1, d))
            and array[2*d - 1] == d * (2*d + 2 - n)):
        return (NearPolygonGraph.FoldedCube, k)

    if (l == 0 and n == 2 * d and d % 2 == 1 and (d - 1) // 2 + 1 == k and
            all(array[d - 1 + i] == (i + 1) // 2 for i in range(1, d + 1))):
        return (NearPolygonGraph.DoubleOdd, k - 1)

    if (l == 0 and n == 2 * d and d % 2 == 1 and
            is_prime_power(array[d + 2] - 1) and
            all(array[d - 1 + i] == q_binomial((i + 1) // 2, 1, array[d + 2] - 1)
                for i in range(1, d + 1)) and
            k == q_binomial((d - 1) // 2 + 1, 1, array[d + 2] - 1)):
        return (NearPolygonGraph.DoubleGrassmann, (array[d + 2] - 1, (d - 1) // 2))

    if (n == 2 * d and k == (l+1) * d and
            all(array[d - 1 + i] == i for i in range(1, d + 1))):
        return (NearPolygonGraph.HammingGraph, (d, l + 2))

    if (n == 2 * d and is_prime_power(array[d + 1] - 1) and
            (l + 1) in [(array[d + 1] - 1) ** e for e in [0, 0.5, 1, 1.5, 2]] and
            k == (l + 1) * q_binomial(d, 1, array[d + 1] - 1) and
            all(array[d - 1 + i] == q_binomial(i, 1, array[d + 1] - 1)
                for i in range(1, d + 1))):
        return (NearPolygonGraph.DualPolarGraph, (d, array[d + 1] - 1,
                log(l + 1, array[d + 1] - 1)))

    # otherwise we don't know the near polygon
    return False


def near_polygon_graph(family, params):
    r"""
    Return the near polygon graph with the given parameters.

    The input is expected to be the result of the function
    :func:`sage.graphs.generators.distance_regular.is_near_polygon`.

    INPUT:

    - ``family`` -- int; an element of the enum ``NearPolygonGraph``.

    - ``params`` -- int or tuple; the parameters needed to construct a graph
      of the family ``family``.

    EXAMPLES::

        sage: from sage.graphs.generators.distance_regular import is_near_polygon, near_polygon_graph
        sage: near_polygon_graph(*is_near_polygon(                                      # needs sage.combinat
        ....:     [6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6]))
        Bipartite double of Odd graph on a set of 11 elements: Graph on 924 vertices
        sage: G=_; G.is_distance_regular(True)                                          # needs sage.combinat
        ([6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, None],
         [None, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6])

    REFERENCES:

    See [BCN1989]_ pp. 198-206 for some theory about near polygons as well as
    a list of known examples.

    TESTS::

        sage: # needs sage.combinat
        sage: near_polygon_graph(12, 9)
        Traceback (most recent call last):
        ...
        ValueError: No known near polygons with the given parameters
        sage: is_near_polygon([2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2])
        (0, 12)
        sage: near_polygon_graph((0, 12))
        Traceback (most recent call last):
        ...
        TypeError: ...near_polygon_graph() takes exactly 2 positional arguments (1 given)
        sage: near_polygon_graph(0, 12)
        Cycle graph: Graph on 12 vertices
        sage: near_polygon_graph(*is_near_polygon([8, 7, 6, 5, 1, 2, 3, 8]))
        Folded Cube Graph: Graph on 128 vertices
    """

    if family == NearPolygonGraph.RegularPolygon:
        from sage.graphs.generators.basic import CycleGraph
        return CycleGraph(params)

    if family == NearPolygonGraph.GeneralisedPolygon:
        d, s, t = params
        if d == 3:
            return GeneralisedHexagonGraph(s, t)
        if d == 4:
            return GeneralisedOctagonGraph(s, t)
        if d == 6:
            return GeneralisedDodecagonGraph(s, t)

    if family == NearPolygonGraph.OddGraph:
        from sage.graphs.generators.families import OddGraph
        return OddGraph(params)

    if family == NearPolygonGraph.DoubleOdd:
        return DoubleOddGraph(params)

    if family == NearPolygonGraph.DoubleGrassmann:
        return DoubleGrassmannGraph(*params)

    if family == NearPolygonGraph.FoldedCube:
        from sage.graphs.generators.families import FoldedCubeGraph
        return FoldedCubeGraph(params)

    if family == NearPolygonGraph.HammingGraph:
        from sage.graphs.generators.families import HammingGraph
        return HammingGraph(*params)

    if family == NearPolygonGraph.DualPolarGraph:
        from sage.graphs.generators.classical_geometries import (
            UnitaryDualPolarGraph,
            OrthogonalDualPolarGraph,
            SymplecticDualPolarGraph)

        d, q, e = params
        if e == 0:
            return OrthogonalDualPolarGraph(1, d, q)
        if e == 0.5:
            return UnitaryDualPolarGraph(2 * d, int(q**0.5))
        if e == 1:
            return SymplecticDualPolarGraph(2 * d, q)
        if e == 1.5:
            return UnitaryDualPolarGraph(2*d + 1, int(q**0.5))
        if e == 2:
            return OrthogonalDualPolarGraph(-1, d, q)

    raise ValueError("No known near polygons with the given parameters")


# dictionary intersection_array (as tuple)  -> construction
# of spordaic distance-regular graphs
_sporadic_graph_database = {
    (3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3): FosterGraph,
    (7, 6, 4, 4, 4, 1, 1, 1, 1, 1, 1, 2, 4, 4, 6, 7): IvanovIvanovFaradjevGraph,
    (3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3): BiggsSmithGraph,
    (22, 21, 20, 16, 6, 2, 1, 1, 2, 6, 16, 20, 21, 22): lambda:
    codes.GolayCode(GF(2), False).punctured([0]).cosetGraph().bipartite_double(),
    (23, 22, 21, 20, 3, 2, 1, 1, 2, 3, 20, 21, 22, 23): lambda:
    codes.GolayCode(GF(2), False).cosetGraph().bipartite_double(),
    (21, 20, 16, 6, 2, 1, 1, 2, 6, 16, 20, 21):
    shortened_00_11_binary_Golay_code_graph,
    (21, 20, 16, 9, 2, 1, 1, 2, 3, 16, 20, 21):
    shortened_000_111_extended_binary_Golay_code_graph,
    (22, 21, 20, 3, 2, 1, 1, 2, 3, 20, 21, 22): lambda:
    codes.GolayCode(GF(2), extended=False).shortened([0]).cosetGraph(),
    (3, 2, 1, 1, 1, 1, 1, 1, 2, 3): DodecahedralGraph,
    (22, 20, 18, 2, 1, 1, 2, 9, 20, 22): lambda:
    codes.GolayCode(GF(3)).shortened([0]).cosetGraph(),
    (7, 6, 6, 1, 1, 1, 1, 6, 6, 7): lambda:
    HoffmanSingletonGraph().bipartite_double(),
    (10, 9, 8, 2, 1, 1, 2, 8, 9, 10): lambda:
    SimsGewirtzGraph().bipartite_double(),
    (16, 15, 12, 4, 1, 1, 4, 12, 15, 16): lambda:
    strongly_regular_graph(77, 16, 0, check=False).bipartite_double(),
    (22, 21, 16, 6, 1, 1, 6, 16, 21, 22): lambda:
    HigmanSimsGraph().bipartite_double(),
    (3, 2, 2, 1, 1, 1, 1, 2): CoxeterGraph,
    (6, 5, 5, 4, 1, 1, 2, 6): vanLintSchrijverGraph,
    (7, 6, 4, 4, 1, 1, 1, 6): DoublyTruncatedWittGraph,
    (9, 8, 6, 3, 1, 1, 3, 8): distance_3_doubly_truncated_Golay_code_graph,
    (10, 8, 8, 2, 1, 1, 4, 5): J2Graph,
    (11, 10, 6, 1, 1, 1, 5, 11): LivingstoneGraph,
    (5, 4, 1, 1, 1, 1, 4, 5): WellsGraph,
    (6, 4, 2, 1, 1, 1, 4, 6): FosterGraph3S6,
    (10, 6, 4, 1, 1, 2, 6, 10):  ConwaySmith_for_3S7,
    (20, 18, 4, 1, 1, 2, 18, 20): lambda:
    codes.GolayCode(GF(3), extended=False).shortened([0]).cosetGraph(),
    (45, 32, 12, 1, 1, 6, 32, 45): locally_GQ42_distance_transitive_graph,
    (117, 80, 24, 1, 1, 12, 80, 117): graph_3O73,
    (22, 21, 20, 1, 2, 6): lambda:
    codes.GolayCode(GF(2), extended=False).punctured([0]).cosetGraph(),
    (23, 22, 21, 1, 2, 3): lambda:
    codes.GolayCode(GF(2), extended=False).cosetGraph(),
    (24, 23, 22, 21, 1, 2, 3, 24): lambda: codes.GolayCode(GF(2)).cosetGraph(),
    (12, 11, 10, 7, 1, 2, 5, 12): LeonardGraph,
    (15, 14, 10, 3, 1, 5, 12, 15): cocliques_HoffmannSingleton,
    (27, 10, 1, 1, 10, 27): GossetGraph,
    (30, 28, 24, 1, 3, 15): LargeWittGraph,
    (15, 14, 12, 1, 1, 9): TruncatedWittGraph,
    (24, 22, 20, 1, 2, 12): lambda: codes.GolayCode(GF(3)).cosetGraph(),
    (21, 20, 16, 1, 2, 12): lambda:
    codes.GolayCode(GF(2), extended=False).punctured([0, 1]).cosetGraph()
}

_infinite_families_database = [
    (is_classical_parameters_graph, graph_with_classical_parameters),
    (is_pseudo_partition_graph, pseudo_partition_graph),
    (is_near_polygon, near_polygon_graph),
    (is_from_GQ_spread, graph_from_GQ_spread),
]


def distance_regular_graph(list arr, existence=False, check=True):
    r"""
    Return a distance-regular graph with the intersection array given.

    INPUT:

    - ``arr`` -- list; intersection array of the graph

    - ``existence`` -- boolean (optional); instead of building the graph return:

      - ``True`` - if a graph with the given intersection array exists;

      - ``False`` - if there is no graph with the given intersection array;

      - ``Unknown`` - if Sage doesn't know if such a graph exists.

    - ``check`` -- boolean (optional); if ``True``, then checks that the result
      of this function has the given intersection array. Default: ``True``

    EXAMPLES::

        sage: graphs.distance_regular_graph([21,20,16,1,2,12], existence=True)
        True
        sage: G = graphs.distance_regular_graph([12,11,10,7,1,2,5,12], check=False)     # needs sage.combinat sage.modules
        sage: G.is_distance_regular(True)                                               # needs sage.combinat sage.modules
        ([12, 11, 10, 7, None], [None, 1, 2, 5, 12])

    REFERENCES:

    See [BCN1989]_ and [VDKT2016]_.

    TESTS::

        sage: graphs.distance_regular_graph([3, 2, 2, 1, 1, 1, 1, 2, 2, 3],             # needs sage.combinat
        ....:                               existence=True)
        True
        sage: graphs.distance_regular_graph([3, 2, 2, 1, 2, 1, 1, 2, 2, 3],
        ....:                               existence=True)
        False
        sage: graphs.distance_regular_graph([18, 16, 16, 1, 1, 9])  # optional - internet gap_packages
        Generalised hexagon of order (2, 8): Graph on 819 vertices

        sage: # needs sage.combinat
        sage: graphs.distance_regular_graph([14, 12, 10, 8, 6, 4, 2,
        ....:                                1, 2, 3, 4, 5, 6, 7])
        Hamming Graph with parameters 7,3: Graph on 2187 vertices
        sage: graphs.distance_regular_graph([66, 45, 28, 1, 6, 30])
        Graph on 1024 vertices
        sage: graphs.distance_regular_graph([6,5,5,5,1,1,1,6])
        Generalised octagon of order (1, 5): Graph on 312 vertices
        sage: graphs.distance_regular_graph([64, 60, 1, 1, 15, 64], check=True)
        Graph on 325 vertices
    """
    from sage.misc.unknown import Unknown
    from sage.categories.sets_cat import EmptySetError

    # check if drg module is installed
    try:
        import drg
        from drg import InfeasibleError
        drgModule = True
    except ModuleNotFoundError:
        drgModule = False

    def result(G):
        if check:
            array = _intersection_array_from_graph(G)
            if array != arr:
                raise RuntimeError(("Sage built the wrong distance-regular "
                                    f"graph; expected {arr}, result {array}"))
        return G

    def is_iterable(obj):
        try:
            iter(obj)
            return True
        except TypeError:
            return False

    n = len(arr)
    d = n // 2

    # check that arr makes sense:
    if drgModule:
        try:
            parameters = drg.DRGParameters(arr[:d], arr[d:])
        except (AssertionError, InfeasibleError, TypeError) as err:
            if existence:
                return False
            raise EmptySetError(("No distance-regular graphs with "
                                 f"parameters {arr} exists; error: {err}"))
    else:
        # basic checks
        if len(arr) % 2 == 1 or any([i <= 0 for i in arr]) or \
           any([x != int(x) for x in arr]) or \
           any([(arr[i] - arr[i + 1]) < 0 for i in range(d - 1)]) or \
           any([(arr[d + i + 1] - arr[d + i]) < 0 for i in range(d - 1)]):
            if existence:
                return False
            raise EmptySetError(("No distance-regular graphs with "
                                 f"parameters {arr} exists"))

    # handle diameter < 3
    if d == 1 and arr[1] == 1:
        if existence:
            return True
        from sage.graphs.generators.basic import CompleteGraph
        return result(CompleteGraph(arr[0] + 1))

    if d == 2:
        k = arr[0]
        mu = arr[3]
        l = k - arr[1] - 1  # a1 = k - b1 - c1
        v = (k * (k - l - 1)) // mu + k + 1

        if existence:
            return strongly_regular_graph(v, k, l, mu, existence=True)
        return result(strongly_regular_graph(v, k, l, mu))

    t = tuple(arr)
    if t in _sporadic_graph_database:
        if existence:
            return True
        return result(_sporadic_graph_database[t]())

    for (f, g) in _infinite_families_database:
        t = f(arr)
        if t is not False:
            if existence:
                return True

            G = g(*t) if is_iterable(t) else g(t)
            return result(G)

    # now try drg feasibility
    if drgModule:
        try:
            parameters.check_feasible()
        except (InfeasibleError, TypeError, AssertionError) as err:
            if existence:
                return False
            raise EmptySetError(("No distance-regular graphs with "
                                 f"parameters {arr} exists; reason: {err}"))

    if existence:
        return Unknown
    raise RuntimeError(
        f"No distance-regular graph with intersection array {arr} known")
