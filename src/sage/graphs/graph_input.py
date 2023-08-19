r"""
Functions for reading/building graphs/digraphs

This module gathers functions needed to build a graph from any other data.

.. NOTE::

    This is an **internal** module of Sage. All features implemented here are
    made available to end-users through the constructors of :class:`Graph` and
    :class:`DiGraph`.

Note that because they are called by the constructors of :class:`Graph` and
:class:`DiGraph`, most of these functions modify a graph inplace.

{INDEX_OF_FUNCTIONS}

Functions
---------

"""
from sage.cpython.string import bytes_to_str
from sage.misc.rest_index_of_methods import gen_rest_table_index
import sys


def from_graph6(G, g6_string):
    r"""
    Fill ``G`` with the data of a graph6 string.

    INPUT:

    - ``G`` -- a graph

    - ``g6_string`` -- a graph6 string

    EXAMPLES::

        sage: from sage.graphs.graph_input import from_graph6
        sage: g = Graph()
        sage: from_graph6(g, 'IheA@GUAo')
        sage: g.is_isomorphic(graphs.PetersenGraph())
        True
    """
    from .generic_graph_pyx import length_and_string_from_graph6, binary_string_from_graph6

    if isinstance(g6_string, bytes):
        g6_string = bytes_to_str(g6_string)
    elif not isinstance(g6_string, str):
        raise ValueError('if input format is graph6, then g6_string must be a string')
    n = g6_string.find('\n')
    if n == -1:
        n = len(g6_string)
    ss = g6_string[:n]
    n, s = length_and_string_from_graph6(ss)
    m = binary_string_from_graph6(s, n)
    expected = n * (n - 1) // 2 + (6 - n * (n - 1) // 2) % 6
    if len(m) > expected:
        raise RuntimeError("the string (%s) seems corrupt: for n = %d, the string is too long" % (ss, n))
    elif len(m) < expected:
        raise RuntimeError("the string (%s) seems corrupt: for n = %d, the string is too short" % (ss, n))
    G.add_vertices(range(n))
    k = 0
    for i in range(n):
        for j in range(i):
            if m[k] == '1':
                G._backend.add_edge(i, j, None, False)
            k += 1


def from_sparse6(G, g6_string):
    r"""
    Fill ``G`` with the data of a sparse6 string.

    INPUT:

    - ``G`` -- a graph

    - ``g6_string`` -- a sparse6 string

    EXAMPLES::

        sage: from sage.graphs.graph_input import from_sparse6
        sage: g = Graph()
        sage: from_sparse6(g, ':I`ES@obGkqegW~')
        sage: g.is_isomorphic(graphs.PetersenGraph())
        True
    """
    from .generic_graph_pyx import length_and_string_from_graph6, int_to_binary_string

    if isinstance(g6_string, bytes):
        g6_string = bytes_to_str(g6_string)
    elif not isinstance(g6_string, str):
        raise ValueError('if input format is graph6, then g6_string must be a string')

    n = g6_string.find('\n')
    if n == -1:
        n = len(g6_string)
    s = g6_string[:n]
    n, s = length_and_string_from_graph6(s[1:])
    if not n:
        edges = []
    else:
        from sage.rings.integer_ring import ZZ
        k = int((ZZ(n) - 1).nbits())
        ords = [ord(i) for i in s]
        if any(o > 126 or o < 63 for o in ords):
            raise RuntimeError("the string seems corrupt: valid characters are \n" + ''.join([chr(i) for i in range(63, 127)]))
        bits = ''.join([int_to_binary_string(o-63).zfill(6) for o in ords])
        if not k:
            b = [int(x) for x in bits]
            x = [0] * len(b)
        else:
            b = []
            x = []
            for i in range(0, len(bits)-k, k+1):
                b.append(int(bits[i:i+1], 2))
                x.append(int(bits[i+1:i+k+1], 2))
        v = 0
        edges = []
        for i in range(len(b)):
            v += b[i]  # +1 if b[i] == 1 else 0
            if x[i] > v:
                v = x[i]
            else:
                if v < n:
                    edges.append((x[i], v))
    G.add_vertices(range(n))
    G.add_edges(edges)


def from_dig6(G, dig6_string):
    r"""
    Fill ``G`` with the data of a dig6 string.

    INPUT:

    - ``G`` -- a graph

    - ``dig6_string`` -- a dig6 string

    EXAMPLES::

        sage: from sage.graphs.graph_input import from_dig6
        sage: g = DiGraph()
        sage: from_dig6(g, digraphs.Circuit(10).dig6_string())
        sage: g.is_isomorphic(digraphs.Circuit(10))
        True
    """
    from .generic_graph_pyx import length_and_string_from_graph6, binary_string_from_dig6
    if isinstance(dig6_string, bytes):
        dig6_string = bytes_to_str(dig6_string)
    elif not isinstance(dig6_string, str):
        raise ValueError('if input format is dig6, then dig6_string must be a string')
    n = dig6_string.find('\n')
    if n == -1:
        n = len(dig6_string)
    ss = dig6_string[:n]
    n, s = length_and_string_from_graph6(ss)
    m = binary_string_from_dig6(s, n)
    expected = n**2
    if len(m) > expected:
        raise RuntimeError("the string (%s) seems corrupt: for n = %d, the string is too long" % (ss, n))
    elif len(m) < expected:
        raise RuntimeError("the string (%s) seems corrupt: for n = %d, the string is too short" % (ss, n))
    G.add_vertices(range(n))
    k = 0
    for i in range(n):
        for j in range(n):
            if m[k] == '1':
                G._backend.add_edge(i, j, None, True)
            k += 1


def from_seidel_adjacency_matrix(G, M):
    r"""
    Fill ``G`` with the data of a Seidel adjacency matrix.

    INPUT:

    - ``G`` -- a graph

    - ``M`` -- a Seidel adjacency matrix

    EXAMPLES::

        sage: from sage.graphs.graph_input import from_seidel_adjacency_matrix
        sage: g = Graph()
        sage: sam = graphs.PetersenGraph().seidel_adjacency_matrix()                    # needs sage.modules
        sage: from_seidel_adjacency_matrix(g, sam)                                      # needs sage.modules
        sage: g.is_isomorphic(graphs.PetersenGraph())                                   # needs sage.modules
        True
    """
    from sage.structure.element import is_Matrix
    from sage.rings.integer_ring import ZZ
    assert is_Matrix(M)

    if M.base_ring() != ZZ:
        try:
            M = M.change_ring(ZZ)
        except TypeError:
            raise ValueError("the adjacency matrix of a Seidel graph must" +
                             " have only 0,1,-1 integer entries")

    if M.is_sparse():
        entries = set(M[i, j] for i, j in M.nonzero_positions())
    else:
        entries = set(M.list())

    if any(e < -1 or e > 1 for e in entries):
        raise ValueError("the adjacency matrix of a Seidel graph must" +
                         " have only 0,1,-1 integer entries")
    if any(i == j for i, j in M.nonzero_positions()):
        raise ValueError("the adjacency matrix of a Seidel graph must" +
                         " have 0s on the main diagonal")
    if not M.is_symmetric():
        raise ValueError("the adjacency matrix of a Seidel graph must be symmetric")

    G.add_vertices(range(M.nrows()))
    G.add_edges((i, j) for i, j in M.nonzero_positions() if i <= j and M[i, j] < 0)


def from_adjacency_matrix(G, M, loops=False, multiedges=False, weighted=False):
    r"""
    Fill ``G`` with the data of an adjacency matrix.

    INPUT:

    - ``G`` -- a :class:`Graph` or :class:`DiGraph`

    - ``M`` -- an adjacency matrix

    - ``loops``, ``multiedges``, ``weighted`` -- booleans (default: ``False``);
      whether to consider the graph as having loops, multiple edges, or weights

    EXAMPLES::

        sage: from sage.graphs.graph_input import from_adjacency_matrix
        sage: g = Graph()
        sage: from_adjacency_matrix(g, graphs.PetersenGraph().adjacency_matrix())       # needs sage.modules
        sage: g.is_isomorphic(graphs.PetersenGraph())                                   # needs sage.modules
        True
    """
    from sage.structure.element import is_Matrix
    from sage.rings.integer_ring import ZZ
    assert is_Matrix(M)
    # note: the adjacency matrix might be weighted and hence not
    # necessarily consists of integers
    if not weighted and M.base_ring() != ZZ:
        try:
            M = M.change_ring(ZZ)
        except TypeError:
            if weighted is False:
                raise ValueError("the adjacency matrix of a non-weighted graph" +
                                 " must have only nonnegative integer entries")
            weighted = True

    if M.is_sparse():
        entries = set(M[i, j] for i, j in M.nonzero_positions())
    else:
        entries = set(M.list())

    if not weighted and any(e < 0 for e in entries):
        if weighted is False:
            raise ValueError("the adjacency matrix of a non-weighted graph" +
                             " must have only nonnegative integer entries")
        weighted = True
        if multiedges is None:
            multiedges = False
    if weighted is None:
        weighted = False

    if multiedges is None:
        multiedges = ((not weighted) and any(e != 0 and e != 1 for e in entries))

    if not loops and any(M[i, i] for i in range(M.nrows())):
        if loops is False:
            raise ValueError("the adjacency matrix of a non-weighted graph" +
                             " must have zeroes on the diagonal")
        loops = True
    if loops is None:
        loops = False
    G.allow_loops(loops, check=False)
    G.allow_multiple_edges(multiedges, check=False)
    G.add_vertices(range(M.nrows()))
    if G.is_directed():
        pairs = M.nonzero_positions()
    else:
        pairs = ((i, j) for i, j in M.nonzero_positions() if i <= j)
    if weighted:
        G.add_edges((i, j, M[i][j]) for i, j in pairs)
    elif multiedges:
        G.add_edges((i, j) for i, j in pairs for _ in range(int(M[i][j])))
    else:
        G.add_edges((i, j) for i, j in pairs)
    G._weighted = weighted


def from_incidence_matrix(G, M, loops=False, multiedges=False, weighted=False):
    r"""
    Fill ``G`` with the data of an incidence matrix.

    INPUT:

    - ``G`` -- a graph

    - ``M`` -- an incidence matrix

    - ``loops``, ``multiedges``, ``weighted`` -- booleans (default: ``False``);
      whether to consider the graph as having loops, multiple edges, or weights

    EXAMPLES::

        sage: from sage.graphs.graph_input import from_incidence_matrix
        sage: g = Graph()
        sage: from_incidence_matrix(g, graphs.PetersenGraph().incidence_matrix())       # needs sage.modules
        sage: g.is_isomorphic(graphs.PetersenGraph())                                   # needs sage.modules
        True
    """
    from sage.structure.element import is_Matrix
    assert is_Matrix(M)

    oriented = any(M[pos] < 0 for pos in M.nonzero_positions(copy=False))

    positions = []
    for i in range(M.ncols()):
        NZ = M.nonzero_positions_in_column(i)
        if len(NZ) == 1:
            if oriented:
                raise ValueError("column {} of the (oriented) incidence "
                                 "matrix contains only one nonzero value".format(i))
            elif M[NZ[0], i] != 2:
                raise ValueError("each column of a non-oriented incidence "
                                 "matrix must sum to 2, but column {} does not".format(i))
            if loops is None:
                loops = True
            positions.append((NZ[0], NZ[0]))
        elif (len(NZ) != 2 or
              (oriented and not ((M[NZ[0], i] == +1 and M[NZ[1], i] == -1) or
                                 (M[NZ[0], i] == -1 and M[NZ[1], i] == +1))) or
              (not oriented and (M[NZ[0], i] != 1 or M[NZ[1], i] != 1))):
            msg = "there must be one or two nonzero entries per column in an incidence matrix, "
            msg += "got entries {} in column {}".format([M[j, i] for j in NZ], i)
            raise ValueError(msg)
        else:
            positions.append(tuple(NZ))

    if weighted is None:
        G._weighted = False
    if multiedges is None:
        total = len(positions)
        multiedges = len(set(positions)) < total
    G.allow_loops(False if loops is None else loops, check=False)
    G.allow_multiple_edges(multiedges, check=False)
    G.add_vertices(range(M.nrows()))
    G.add_edges(positions)


def from_oriented_incidence_matrix(G, M, loops=False, multiedges=False, weighted=False):
    r"""
    Fill ``G`` with the data of an *oriented* incidence matrix.

    An oriented incidence matrix is the incidence matrix of a directed graph, in
    which each non-loop edge corresponds to a `+1` and a `-1`, indicating its
    source and destination.

    INPUT:

    - ``G`` -- a :class:`DiGraph`

    - ``M`` -- an incidence matrix

    - ``loops``, ``multiedges``, ``weighted`` -- booleans (default: ``False``);
      whether to consider the graph as having loops, multiple edges, or weights

    .. NOTE:: ``weighted`` is currently ignored.

    EXAMPLES::

        sage: from sage.graphs.graph_input import from_oriented_incidence_matrix
        sage: g = DiGraph()
        sage: im = digraphs.Circuit(10).incidence_matrix()                              # needs sage.modules
        sage: from_oriented_incidence_matrix(g, im)                                     # needs sage.modules
        sage: g.is_isomorphic(digraphs.Circuit(10))                                     # needs sage.modules
        True

    TESTS:

    Fix bug reported in :trac:`22985`::

        sage: DiGraph(matrix ([[1,0,0,1],[0,0,1,1],[0,0,1,1]]).transpose())             # needs sage.modules
        Traceback (most recent call last):
        ...
        ValueError: each column represents an edge: -1 goes to 1

    Handle incidence matrix containing a column with only zeros (:trac:`29275`)::

        sage: m = Matrix([[0,1],[0,-1],[0,0]]); m                                       # needs sage.modules
        [ 0  1]
        [ 0 -1]
        [ 0  0]
        sage: G = DiGraph(m, format='incidence_matrix')                                 # needs sage.modules
        sage: list(G.edges(sort=True, labels=False))                                    # needs sage.modules
        [(1, 0)]

    Handle incidence matrix [[1],[-1]] (:trac:`29275`)::

        sage: m = Matrix([[1],[-1]]); m                                                 # needs sage.modules
        [ 1]
        [-1]
        sage: G = DiGraph(m, format='incidence_matrix')                                 # needs sage.modules
        sage: list(G.edges(sort=True, labels=False))                                    # needs sage.modules
        [(1, 0)]
    """
    from sage.structure.element import is_Matrix
    assert is_Matrix(M)

    positions = []
    for c in M.columns():
        NZ = c.nonzero_positions()
        if not NZ:
            continue
        if len(NZ) != 2:
            raise ValueError("there must be two nonzero entries (-1 & 1) per column")
        L = sorted([c[i] for i in NZ])
        if L != [-1, 1]:
            raise ValueError("each column represents an edge: -1 goes to 1")
        if c[NZ[0]] == -1:
            positions.append(tuple(NZ))
        else:
            positions.append((NZ[1], NZ[0]))
    if multiedges is None:
        total = len(positions)
        multiedges = len(set(positions)) < total
    G.allow_loops(bool(loops), check=False)
    G.allow_multiple_edges(multiedges, check=False)
    G.add_vertices(range(M.nrows()))
    G.add_edges(positions)


def from_dict_of_dicts(G, M, loops=False, multiedges=False, weighted=False, convert_empty_dict_labels_to_None=False):
    r"""
    Fill ``G`` with the data of a dictionary of dictionaries.

    INPUT:

    - ``G`` -- a graph

    - ``M`` -- a dictionary of dictionaries

    - ``loops``, ``multiedges``, ``weighted`` -- booleans (default: ``False``);
      whether to consider the graph as having loops, multiple edges, or weights

    - ``convert_empty_dict_labels_to_None`` -- booleans (default: ``False``);
      whether to adjust for empty dicts instead of ``None`` in NetworkX default
      edge labels

    EXAMPLES::

        sage: from sage.graphs.graph_input import from_dict_of_dicts
        sage: g = Graph()
        sage: from_dict_of_dicts(g, graphs.PetersenGraph().to_dictionary(edge_labels=True))
        sage: g.is_isomorphic(graphs.PetersenGraph())
        True

    TESTS:

    :trac:`32831` is fixed::

        sage: DiGraph({0: {}, 1: {}, 2: {}, 3: {}, 4: {}})
        Digraph on 5 vertices
    """
    if any(not isinstance(M[u], dict) for u in M):
        raise ValueError("input dict must be a consistent format")

    if not loops:
        if any(u in neighb for u, neighb in M.items()):
            if loops is False:
                u = next(u for u, neighb in M.items() if u in neighb)
                raise ValueError("the graph was built with loops=False but input M has a loop at {}".format(u))
            loops = True
        if loops is None:
            loops = False
    if weighted is None:
        G._weighted = False
    input_multiedges = multiedges
    if multiedges is not False:
        if not all(isinstance(M[u][v], list) for u in M for v in M[u]):
            if multiedges:
                raise ValueError("dict of dicts for multigraph must be in the format {v: {u: list}}")
            multiedges = False
        if multiedges is None and M:
            multiedges = True

    G.allow_loops(loops, check=False)
    G.allow_multiple_edges(multiedges, check=False)
    verts = set().union(M.keys(), *M.values())
    G.add_vertices(verts)
    if convert_empty_dict_labels_to_None:
        def relabel(x):
            return x if x != {} else None
    else:
        def relabel(x):
            return x

    is_directed = G.is_directed()
    if not is_directed and multiedges:
        v_to_id = {v: i for i, v in enumerate(verts)}
        for u in M:
            for v in M[u]:
                if v_to_id[u] <= v_to_id[v] or v not in M or u not in M[v] or u == v:
                    for label in M[u][v]:
                        G._backend.add_edge(u, v, relabel(label), False)
    elif multiedges:
        for u in M:
            for v in M[u]:
                for label in M[u][v]:
                    G._backend.add_edge(u, v, relabel(label), is_directed)
    else:
        for u in M:
            for v in M[u]:
                G._backend.add_edge(u, v, relabel(M[u][v]), is_directed)
    if not G.size() and input_multiedges is not True:
        G.allow_multiple_edges(False, check=False)


def from_dict_of_lists(G, D, loops=False, multiedges=False, weighted=False):
    r"""
    Fill ``G`` with the data of a dictionary of lists.

    INPUT:

    - ``G`` -- a :class:`Graph` or :class:`DiGraph`

    - ``D`` -- a dictionary of lists

    - ``loops``, ``multiedges``, ``weighted`` -- booleans (default: ``False``);
      whether to consider the graph as having loops, multiple edges, or weights

    EXAMPLES::

        sage: from sage.graphs.graph_input import from_dict_of_lists
        sage: g = Graph()
        sage: from_dict_of_lists(g, graphs.PetersenGraph().to_dictionary())
        sage: g.is_isomorphic(graphs.PetersenGraph())
        True
    """
    verts = set().union(D.keys(), *D.values())
    if not loops:
        if any(u in neighb for u, neighb in D.items()):
            if loops is False:
                u = next(u for u, neighb in D.items() if u in neighb)
                raise ValueError("the graph was built with loops=False but input D has a loop at {}".format(u))
            loops = True
        if loops is None:
            loops = False
    if weighted is None:
        G._weighted = False
    if not multiedges:
        for u in D:
            if len(set(D[u])) != len(D[u]):
                if multiedges is False:
                    v = next((v for v in D[u] if D[u].count(v) > 1))
                    raise ValueError("non-multigraph got several edges (%s, %s)" % (u, v))
                multiedges = True
                break
        if multiedges is None:
            multiedges = False
    G.allow_loops(loops, check=False)
    G.allow_multiple_edges(multiedges, check=False)
    G.add_vertices(verts)

    is_directed = G.is_directed()
    if not is_directed and multiedges:
        v_to_id = {v: i for i, v in enumerate(verts)}
        for u in D:
            for v in D[u]:
                if (v_to_id[u] <= v_to_id[v] or
                        v not in D or u not in D[v] or u == v):
                    G._backend.add_edge(u, v, None, False)
    else:
        for u in D:
            for v in D[u]:
                G._backend.add_edge(u, v, None, is_directed)


def from_networkx_graph(G, gnx, weighted=None, loops=None, multiedges=None,
                        convert_empty_dict_labels_to_None=None):
    r"""
    Fill `G` with the data of a NetworkX (di)graph.

    INPUT:

    - ``G`` -- a :class:`Graph` or :class:`DiGraph`

    - ``gnx`` -- a NetworkX ``Graph``, ``MultiGraph``, ``DiGraph`` or
      ``MultiDiGraph``

    - ``weighted`` -- boolean (default: ``None``); whether graph thinks of
      itself as weighted or not. See
      :meth:`~sage.graphs.generic_graph.GenericGraph.weighted`.

    - ``loops`` -- boolean (default: ``None``); whether to allow loops

    - ``multiedges`` -- boolean (default: ``None``); whether to allow multiple
      edges

    - ``convert_empty_dict_labels_to_None`` -- boolean (default: ``None``);
      whether to replace the default edge labels used by NetworkX (empty
      dictionaries) by ``None``, the default Sage edge label. When set to
      ``False``, empty dictionaries are not converted to ``None``.

    EXAMPLES:

    Feeding a :class:`Graph` with a NetworkX ``Graph``::

        sage: # needs networkx
        sage: from sage.graphs.graph_input import from_networkx_graph
        sage: import networkx
        sage: G = Graph()
        sage: _ = gnx = networkx.Graph()
        sage: _ = gnx.add_edge(0, 1)
        sage: _ = gnx.add_edge(1, 2)
        sage: from_networkx_graph(G, gnx)
        sage: G.edges(sort=True, labels=False)
        [(0, 1), (1, 2)]

    Feeding a :class:`Graph` with a NetworkX ``MultiGraph``::

        sage: # needs networkx
        sage: G = Graph()
        sage: gnx = networkx.MultiGraph()
        sage: _ = gnx.add_edge(0, 1)
        sage: _ = gnx.add_edge(0, 1)
        sage: from_networkx_graph(G, gnx)
        sage: G.edges(sort=True, labels=False)
        [(0, 1), (0, 1)]
        sage: G = Graph()
        sage: from_networkx_graph(G, gnx, multiedges=False)
        sage: G.edges(sort=True, labels=False)
        [(0, 1)]

    When feeding a :class:`Graph` `G` with a NetworkX ``DiGraph`` `D`, `G` has
    one edge `(u, v)` whenever `D` has arc `(u, v)` or `(v, u)` or both::

        sage: # needs networkx
        sage: G = Graph()
        sage: D = networkx.DiGraph()
        sage: _ = D.add_edge(0, 1)
        sage: from_networkx_graph(G, D)
        sage: G.edges(sort=True, labels=False)
        [(0, 1)]
        sage: G = Graph()
        sage: _ = D.add_edge(1, 0)
        sage: from_networkx_graph(G, D)
        sage: G.edges(sort=True, labels=False)
        [(0, 1)]

    When feeding a :class:`Graph` `G` with a NetworkX ``MultiDiGraph`` `D`, the
    number of edges between `u` and `v` in `G` is the maximum between the number
    of arcs `(u, v)` and the number of arcs `(v, u)` in D`::

        sage: # needs networkx
        sage: G = Graph()
        sage: D = networkx.MultiDiGraph()
        sage: _ = D.add_edge(0, 1)
        sage: _ = D.add_edge(1, 0)
        sage: _ = D.add_edge(1, 0)
        sage: D.edges()
        OutMultiEdgeDataView([(0, 1), (1, 0), (1, 0)])
        sage: from_networkx_graph(G, D)
        sage: G.edges(sort=True, labels=False)
        [(0, 1), (0, 1)]

    Feeding a :class:`DiGraph` with a NetworkX ``DiGraph``::

        sage: # needs networkx
        sage: from sage.graphs.graph_input import from_networkx_graph
        sage: import networkx
        sage: G = DiGraph()
        sage: _ = gnx = networkx.DiGraph()
        sage: _ = gnx.add_edge(0, 1)
        sage: _ = gnx.add_edge(1, 2)
        sage: from_networkx_graph(G, gnx)
        sage: G.edges(sort=True, labels=False)
        [(0, 1), (1, 2)]

    Feeding a :class:`DiGraph` with a NetworkX ``MultiDiGraph``::

        sage: # needs networkx
        sage: G = DiGraph()
        sage: gnx = networkx.MultiDiGraph()
        sage: _ = gnx.add_edge(0, 1)
        sage: _ = gnx.add_edge(0, 1)
        sage: from_networkx_graph(G, gnx)
        sage: G.edges(sort=True, labels=False)
        [(0, 1), (0, 1)]
        sage: G = DiGraph()
        sage: from_networkx_graph(G, gnx, multiedges=False)
        sage: G.edges(sort=True, labels=False)
        [(0, 1)]

    When feeding a :class:`DiGraph` `G` with a NetworkX ``Graph`` `H`, `G` has
    both arcs `(u, v)` and `(v, u)` if `G` has edge `(u, v)`::

        sage: # needs networkx
        sage: G = DiGraph()
        sage: H = networkx.Graph()
        sage: _ = H.add_edge(0, 1)
        sage: from_networkx_graph(G, H)
        sage: G.edges(labels=False, sort=True)
        [(0, 1), (1, 0)]

    When feeding a :class:`DiGraph` `G` with a NetworkX ``MultiGraph`` `H`, `G`
    has `k` arcs `(u, v)` and `k` arcs `(v, u)` if `H` has `k` edges `(u, v)`,
    unless parameter ``multiedges`` is set to ``False``::

        sage: # needs networkx
        sage: G = DiGraph()
        sage: H = networkx.MultiGraph()
        sage: _ = H.add_edge(0, 1)
        sage: _ = H.add_edge(0, 1)
        sage: _ = H.add_edge(0, 1)
        sage: H.edges()
        MultiEdgeDataView([(0, 1), (0, 1), (0, 1)])
        sage: from_networkx_graph(G, H)
        sage: G.edges(labels=False, sort=True)
        [(0, 1), (0, 1), (0, 1), (1, 0), (1, 0), (1, 0)]
        sage: G = DiGraph()
        sage: from_networkx_graph(G, H, multiedges=False)
        sage: G.edges(labels=False, sort=True)
        [(0, 1), (1, 0)]

    TESTS:

    The first parameter must be a :class:`Graph` or :class:`DiGraph`::

        sage: from sage.graphs.graph_input import from_networkx_graph
        sage: from_networkx_graph("foo", "bar")
        Traceback (most recent call last):
        ...
        ValueError: the first parameter must a Sage Graph or DiGraph

    The second parameter must be a NetworkX ``Graph``, ``MultiGraph``,
      ``DiGraph`` or ``MultiDiGraph``::

        sage: from sage.graphs.graph_input import from_networkx_graph
        sage: from_networkx_graph(Graph(), "bar")                                       # needs networkx
        Traceback (most recent call last):
        ...
        ValueError: the second parameter must be a NetworkX (Multi)(Di)Graph
    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, (Graph, DiGraph)):
        raise ValueError("the first parameter must a Sage Graph or DiGraph")
    import networkx
    if not isinstance(gnx, (networkx.Graph, networkx.DiGraph)):
        raise ValueError("the second parameter must be a NetworkX (Multi)(Di)Graph")

    if G.is_directed() != gnx.is_directed():
        if gnx.is_directed():
            gnx = gnx.to_undirected()
        else:
            gnx = gnx.to_directed()

    if weighted is None:
        if multiedges is None:
            multiedges = gnx.is_multigraph()
        if loops is None:
            loops = any(u == v for u, v in gnx.edges())

    G.allow_loops(loops, check=False)
    G.allow_multiple_edges(multiedges, check=False)
    G.add_vertices(gnx.nodes())
    G.set_vertices(gnx.nodes(data=True))
    if convert_empty_dict_labels_to_None is not False:
        def r(label):
            return None if label == {} else label
        G.add_edges((u, v, r(ll)) for u, v, ll in gnx.edges(data=True))
    else:
        G.add_edges(gnx.edges(data=True))


__doc__ = __doc__.format(INDEX_OF_FUNCTIONS=gen_rest_table_index(sys.modules[__name__]))
