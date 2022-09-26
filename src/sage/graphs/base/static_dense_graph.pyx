# cython: binding=True
r"""
Static dense graphs

This module gathers everything which is related to static dense graphs, i.e. :

- The vertices are integer from `0` to `n-1`
- No labels on vertices/edges
- No multiple edges
- No addition/removal of vertices

This being said, it is technically possible to add/remove edges. The data
structure does not mind at all.

It is all based on the binary matrix data structure described in
``data_structures/binary_matrix.pxd``, which is almost a copy of the bitset data
structure. The only difference is that it differentiates the rows (the vertices)
instead of storing the whole data in a long bitset, and we can use that.

For an overview of graph data structures in sage, see
:mod:`~sage.graphs.base.overview`.

Index
-----

**Cython functions**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    ``dense_graph_init`` | Fill a binary matrix with the information from a Sage (di)graph.

**Python functions**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`is_strongly_regular` | Check whether the graph is strongly regular
    :meth:`is_triangle_free` | Check whether `G` is triangle free
    :meth:`triangles_count` | Return the number of triangles containing `v`, for every `v`
    :meth:`connected_subgraph_iterator` | Iterator over the induced connected subgraphs of order at most `k`

Functions
---------
"""
from sage.data_structures.binary_matrix cimport *
from cysignals.signals cimport sig_on, sig_off, sig_check


cdef dict dense_graph_init(binary_matrix_t m, g, translation=None, force_undirected=False):
    r"""
    Fill a binary matrix with the information from a Sage (di)graph.

    INPUT:

    - ``binary_matrix_t m`` -- the binary matrix to be filled

    - ``g`` -- a graph or digraph

    - ``translation`` -- (default: `None``); several options for this parameter
      used to specify the mapping from vertices to integers:

      - ``True``, ``False``, ``None`` -- the `i`-th vertex in the binary matrix
        corresponds to vertex ``g.vertices(sort=True)[i]``.
        When set to ``True``, a dictionary encoding the mapping from the
        vertices of `g` to integers in `(0, \dots, n-1)` is returned.

      - a ``list`` -- defines a mapping from integers in `(0, \dots, n-1)` to
        the vertices in the graph `g`. So the `i`-th vertex in the binary matrix
        corresponds to vertex ``translation[i]``.
        **Beware that no checks are made that this input is correct**.

      - a ``dict`` -- defines a mapping from the vertices of the graph to
        integers in `(0, \dots, n-1)`. So the `i`-th vertex in the binary matrix
        corresponds to ``translation[v]`` for some vertex `v \in g`.
        **Beware that no checks are made that this input is correct**.

    - ``force_undirected`` -- boolean (default: ``False``); whether to consider
      the graph as undirected or not
    """
    cdef dict d_translation = {}
    from sage.graphs.graph import Graph
    cdef bint is_undirected = isinstance(g, Graph) or force_undirected
    cdef int n = g.order()
    cdef int i, j

    binary_matrix_init(m, n, n)

    if translation and translation is not True:
        if isinstance(translation, list):
            # We are given a mapping from integers to vertices
            d_translation = {v: i for i, v in enumerate(translation)}
        elif isinstance(translation, dict):
            # We are given a mapping vertices to integers
            d_translation = translation

    # If the vertices are 0...n-1, let's avoid an unnecessary dictionary
    if not d_translation and set(g.vertex_iterator()) == set(range(n)):
        if translation is True:
            d_translation = {i: i for i in range(n)}

        for i, j in g.edge_iterator(labels=False):
            binary_matrix_set1(m, i, j)
            if is_undirected:
                binary_matrix_set1(m, j, i)
    else:
        if not d_translation:
            d_translation = {v: i for i, v in enumerate(g.vertices(sort=True))}

        for u, v in g.edge_iterator(labels=False):
            binary_matrix_set1(m, d_translation[u], d_translation[v])
            if is_undirected:
                binary_matrix_set1(m, d_translation[v], d_translation[u])

    if translation is True:
        return d_translation


def is_strongly_regular(g, parameters=False):
    r"""
    Check whether the graph is strongly regular.

    A simple graph `G` is said to be strongly regular with parameters
    `(n, k, \lambda, \mu)` if and only if:

    * `G` has `n` vertices

    * `G` is `k`-regular

    * Any two adjacent vertices of `G` have `\lambda` common neighbors

    * Any two non-adjacent vertices of `G` have `\mu` common neighbors

    By convention, the complete graphs, the graphs with no edges and the empty
    graph are not strongly regular.

    See the :wikipedia:`Strongly regular graph`.

    INPUT:

    - ``parameters`` -- boolean (default: ``False``); whether to return the
      quadruple `(n, k, \lambda, \mu)`. If ``parameters = False`` (default),
      this method only returns ``True`` and ``False`` answers.
      If ``parameters = True``, the ``True`` answers are replaced by quadruples
      `(n, k, \lambda, \mu)`. See definition above.

    EXAMPLES:

    Petersen's graph is strongly regular::

        sage: g = graphs.PetersenGraph()
        sage: g.is_strongly_regular()
        True
        sage: g.is_strongly_regular(parameters=True)
        (10, 3, 0, 1)

    And Clebsch's graph is too::

        sage: g = graphs.ClebschGraph()
        sage: g.is_strongly_regular()
        True
        sage: g.is_strongly_regular(parameters=True)
        (16, 5, 0, 2)

    But Chvatal's graph is not::

        sage: g = graphs.ChvatalGraph()
        sage: g.is_strongly_regular()
        False

    Complete graphs are not strongly regular. (:trac:`14297`) ::

        sage: g = graphs.CompleteGraph(5)
        sage: g.is_strongly_regular()
        False

    Completements of complete graphs are not strongly regular::

        sage: g = graphs.CompleteGraph(5).complement()
        sage: g.is_strongly_regular()
        False

    The empty graph is not strongly regular::

        sage: g = graphs.EmptyGraph()
        sage: g.is_strongly_regular()
        False

    If the input graph has loops or multiedges an exception is raised::

        sage: Graph([(1,1),(2,2)],loops=True).is_strongly_regular()
        Traceback (most recent call last):
        ...
        ValueError: This method is not known to work on graphs with
        loops. Perhaps this method can be updated to handle them, but in the
        meantime if you want to use it please disallow loops using
        allow_loops().

        sage: Graph([(1,2),(1,2)],multiedges=True).is_strongly_regular()
        Traceback (most recent call last):
        ...
        ValueError: This method is not known to work on graphs with
        multiedges. Perhaps this method can be updated to handle them, but in
        the meantime if you want to use it please disallow multiedges using
        allow_multiple_edges().
    """
    g._scream_if_not_simple()
    cdef binary_matrix_t m
    cdef bitset_t b_tmp
    cdef int n = g.order()
    cdef int inter
    cdef int i, j, l, k

    if not g.order() or not g.size():  # no vertices or no edges
        return False

    if g.is_clique():
        return False

    cdef list degree = g.degree()
    k = degree[0]
    if any(d != k for d in degree):
        return False

    bitset_init(b_tmp, n)

    # m is now our copy of the graph
    dense_graph_init(m, g, translation={v: i for i, v in enumerate(g)})

    cdef int llambda = -1
    cdef int mu = -1

    for i in range(n):
        for j in range(i + 1, n):

            # The intersection of the common neighbors of i and j is a AND of
            # their respective rows. A popcount then returns its cardinality.
            bitset_and(b_tmp, m.rows[i], m.rows[j])
            inter = bitset_len(b_tmp)

            # Check that this cardinality is correct according to the values of lambda and mu
            if binary_matrix_get(m, i, j):
                if llambda == -1:
                    llambda = inter
                elif llambda != inter:
                    binary_matrix_free(m)
                    bitset_free(b_tmp)
                    return False
            else:
                if mu == -1:
                    mu = inter
                elif mu != inter:
                    binary_matrix_free(m)
                    bitset_free(b_tmp)
                    return False

    binary_matrix_free(m)
    bitset_free(b_tmp)

    if parameters:
        return (n, k, llambda, mu)
    else:
        return True


def is_triangle_free(G, certificate=False):
    r"""
    Check whether `G` is triangle free.

    INPUT:

    - ``G`` -- a Sage graph

    - ``certificate`` -- boolean (default: ``False``); whether to return a
      triangle if one is found

    EXAMPLES::

        sage: from sage.graphs.base.static_dense_graph import is_triangle_free
        sage: is_triangle_free(graphs.PetersenGraph())
        True
        sage: K4 = graphs.CompleteGraph(4)
        sage: is_triangle_free(K4)
        False
        sage: b, certif = is_triangle_free(K4, certificate=True)
        sage: K4.subgraph(certif).is_clique()
        True

    TESTS::

        sage: from sage.graphs.base.static_dense_graph import is_triangle_free
        sage: is_triangle_free(Graph())
        True
        sage: is_triangle_free(Graph(), certificate=True)
        (True, [])
    """
    G._scream_if_not_simple()
    cdef int n = G.order()
    if n < 3:
        return (True, []) if certificate else True

    cdef binary_matrix_t g
    cdef list int_to_vertex = list(G)
    dense_graph_init(g, G, translation=int_to_vertex)

    cdef mp_size_t i, j, k
    for i in range(n):
        j = bitset_next(g.rows[i], i + 1)
        while j != -1:
            if bitset_are_disjoint(g.rows[i], g.rows[j]):
                j = bitset_next(g.rows[i], j + 1)
            else:
                if certificate:
                    # Search for a common neighbor
                    k = bitset_first(g.rows[i])
                    while bitset_not_in(g.rows[j], k):
                        k = bitset_next(g.rows[j], k + 1)
                    certif = [int_to_vertex[i], int_to_vertex[j], int_to_vertex[k]]

                binary_matrix_free(g)
                return (False, certif) if certificate else False

    binary_matrix_free(g)
    return (True, []) if certificate else True


def triangles_count(G):
    r"""
    Return the number of triangles containing `v`, for every `v`.

    INPUT:

    - ``G`` -- a simple Sage graph

    EXAMPLES::

        sage: from sage.graphs.base.static_dense_graph import triangles_count
        sage: triangles_count(graphs.PetersenGraph())
        {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0}
        sage: sum(triangles_count(graphs.CompleteGraph(15)).values()) == 3 * binomial(15, 3)
        True
    """
    from sage.rings.integer import Integer
    G._scream_if_not_simple()
    cdef int n = G.order()

    cdef uint64_t * count = <uint64_t *> check_calloc(n, sizeof(uint64_t))

    cdef binary_matrix_t g
    cdef list int_to_vertex = list(G)
    dense_graph_init(g, G, translation=int_to_vertex)

    cdef bitset_t b_tmp
    bitset_init(b_tmp, n)

    cdef int i, j
    cdef uint64_t tmp_count = 0

    for i in range(n):
        for j in range(i + 1, n):
            if not bitset_in(g.rows[i], j):
                continue
            bitset_and(b_tmp, g.rows[i], g.rows[j])
            tmp_count = bitset_len(b_tmp)
            count[i] += tmp_count
            count[j] += tmp_count

    ans = {v: Integer(count[i] // 2) for i, v in enumerate(int_to_vertex)}

    bitset_free(b_tmp)
    binary_matrix_free(g)
    sig_free(count)

    return ans


def connected_subgraph_iterator(G, k=None, bint vertices_only=False):
    r"""
    Iterator over the induced connected subgraphs of order at most `k`.

    This method implements a iterator over the induced connected subgraphs of
    the input (di)graph. An induced subgraph of a graph is another graph, formed
    from a subset of the vertices of the graph and all of the edges connecting
    pairs of vertices in that subset (:wikipedia:`Induced_subgraph`).

    As for method :meth:`sage.graphs.generic_graph.connected_components`, edge
    orientation is ignored. Hence, the directed graph with a single arc `0 \to
    1` is considered connected.

    INPUT:

    - ``G`` -- a :class:`Graph` or a :class:`DiGraph`; loops and multiple edges
      are allowed

    - ``k`` -- (optional) integer; maximum order of the connected subgraphs to
      report; by default, the method iterates over all connected subgraphs
      (equivalent to ``k == n``)

    - ``vertices_only`` -- boolean (default: ``False``); whether to return
      (Di)Graph or list of vertices

    EXAMPLES::

        sage: G = DiGraph([(1, 2), (2, 3), (3, 4), (4, 2)])
        sage: list(G.connected_subgraph_iterator())
        [Subgraph of (): Digraph on 1 vertex,
         Subgraph of (): Digraph on 2 vertices,
         Subgraph of (): Digraph on 3 vertices,
         Subgraph of (): Digraph on 4 vertices,
         Subgraph of (): Digraph on 3 vertices,
         Subgraph of (): Digraph on 1 vertex,
         Subgraph of (): Digraph on 2 vertices,
         Subgraph of (): Digraph on 3 vertices,
         Subgraph of (): Digraph on 2 vertices,
         Subgraph of (): Digraph on 1 vertex,
         Subgraph of (): Digraph on 2 vertices,
         Subgraph of (): Digraph on 1 vertex]
        sage: list(G.connected_subgraph_iterator(vertices_only=True))
        [[1], [1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 4],
         [2], [2, 3], [2, 3, 4], [2, 4], [3], [3, 4], [4]]
        sage: list(G.connected_subgraph_iterator(k=2))
        [Subgraph of (): Digraph on 1 vertex,
         Subgraph of (): Digraph on 2 vertices,
         Subgraph of (): Digraph on 1 vertex,
         Subgraph of (): Digraph on 2 vertices,
         Subgraph of (): Digraph on 2 vertices,
         Subgraph of (): Digraph on 1 vertex,
         Subgraph of (): Digraph on 2 vertices,
         Subgraph of (): Digraph on 1 vertex]
        sage: list(G.connected_subgraph_iterator(k=2, vertices_only=True))
        [[1], [1, 2], [2], [2, 3], [2, 4], [3], [3, 4], [4]]

        sage: G = DiGraph([(1, 2), (2, 1)])
        sage: list(G.connected_subgraph_iterator())
        [Subgraph of (): Digraph on 1 vertex,
         Subgraph of (): Digraph on 2 vertices,
         Subgraph of (): Digraph on 1 vertex]
        sage: list(G.connected_subgraph_iterator(vertices_only=True))
        [[1], [1, 2], [2]]

    TESTS:

    The Path Graph of order `n` has `n (n + 1) / 2` connected subgraphs::

        sage: G = graphs.PathGraph(10)
        sage: len(list(G.connected_subgraph_iterator(vertices_only=True)))
        55
        sage: len(list(G.connected_subgraph_iterator(vertices_only=False)))
        55

    The Complete Graph of order `n` has `2^n - 1` connected subgraphs::

        sage: G = graphs.CompleteGraph(5)
        sage: len(list(G.connected_subgraph_iterator(vertices_only=False))) == 2**G.order() - 1
        True
        sage: G = graphs.CompleteGraph(6)
        sage: len(list(G.connected_subgraph_iterator(vertices_only=True))) == 2**G.order() - 1
        True

    Checks that it works with general graphs and corner cases::

        sage: G = DiGraph([(1, 2), (1, 2)], multiedges=True)
        sage: len(list(G.connected_subgraph_iterator()))
        3
        sage: len(list(G.connected_subgraph_iterator(k=0)))
        0
        sage: len(list(G.connected_subgraph_iterator(vertices_only=True)))
        3

        sage: G = Graph([(1, 2), (1, 1)], loops=True)
        sage: len(list(G.connected_subgraph_iterator(vertices_only=False)))
        3
        sage: G = Graph([(1, 2), (1, 2), (1, 1)], loops=True, multiedges=True)
        sage: len(list(G.connected_subgraph_iterator()))
        3
        sage: len(list(G.connected_subgraph_iterator(k=1)))
        2
        sage: len(list(G.connected_subgraph_iterator(vertices_only=True)))
        3
    """
    cdef Py_ssize_t mk = G.order() if k is None else k
    cdef Py_ssize_t n = G.order()
    if not n or mk < 1:
        return

    cdef list int_to_vertex = list(G)
    cdef binary_matrix_t DG
    sig_on()
    dense_graph_init(DG, G, translation=int_to_vertex, force_undirected=True)

    # We use a stack of bitsets. We need 3 bitsets per level with at most n + 1
    # levels, so 3 * n + 3 bitsets. We also need 1 bitset that we create at the
    # same time, so we need overall 3 * n + 4 bitsets
    cdef binary_matrix_t stack
    binary_matrix_init(stack, 3 * n + 4, n)
    sig_off()

    cdef bitset_t current  # current subset of vertices
    cdef bitset_t left     # remaining vertices to consider
    cdef bitset_t boundary  # neighbors of the current subset
    # candidate vertices for extending the current subset, i.e., vertices that
    # are both in left and in boundary
    cdef bitset_t candidates = stack.rows[3 * n + 3]

    cdef Py_ssize_t level
    cdef Py_ssize_t u, v, a

    # We first generate subsets containing vertex 0, then the subsets containing
    # vertex 1 but not vertex 0 since we have already generated all subsets
    # containing 0, then subsets containing 2 but not 0 or 1, etc.
    for u in range(n):
        sig_check()

        if vertices_only:
            yield [int_to_vertex[u]]
        else:
            yield G.subgraph([int_to_vertex[u]])

        # We initialize the loop with vertices u in current, {u+1, ..., n-1}
        # in left, and N(u) in boundary
        bitset_clear(stack.rows[0])
        bitset_add(stack.rows[0], u)
        bitset_set_first_n(stack.rows[1], u + 1)
        bitset_complement(stack.rows[1], stack.rows[1])
        bitset_copy(stack.rows[2], DG.rows[u])
        level = 0

        while level >= 0:
            sig_check()

            # We take the values at the top of the stack
            current = stack.rows[level]
            left = stack.rows[level + 1]
            boundary = stack.rows[level + 2]

            bitset_and(candidates, left, boundary)

            # Search for a candidate vertex v
            v = bitset_next(candidates, u + 1)
            if v >= 0 and bitset_len(current) < mk:
                # We select vertex v
                bitset_discard(left, v)

                # Since we have not modified 'level', the bitsets for iterating
                # without v are already at the top of the stack, with correct
                # values

                # We also build the subset with v and so we add values at the
                # top of the stack
                level += 3
                bitset_copy(stack.rows[level], current)
                bitset_add(stack.rows[level], v)
                bitset_copy(stack.rows[level + 1], left)
                bitset_union(stack.rows[level + 2], boundary, DG.rows[v])

                # We yield that new subset
                if vertices_only:
                    yield [int_to_vertex[a] for a in range(u, n)
                           if bitset_in(stack.rows[level], a)]
                else:
                    yield G.subgraph([int_to_vertex[a] for a in range(u, n)
                                      if bitset_in(stack.rows[level], a)])

            else:
                # We cannot extend the current subset, either due to a lack of
                # candidate (v == -1), or because the current subset has maximum
                # size. We pop
                level -= 3

    sig_on()
    binary_matrix_free(stack)
    binary_matrix_free(DG)
    sig_off()
