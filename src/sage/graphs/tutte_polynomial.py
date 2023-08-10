r"""
Tutte polynomial

This module implements a deletion-contraction algorithm for computing
the Tutte polynomial as described in the paper [HPR2010]_.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`tutte_polynomial` | Computes the Tutte polynomial of the input graph

Authors:

- Mike Hansen (06-2013), Implemented the algorithm.
- Jernej Azarija (06-2013), Tweaked the code, added documentation

Definition
-----------

Given a graph `G`, with `n` vertices and `m` edges and `k(G)`
connected components we define the Tutte polynomial of `G` as

.. MATH::

    \sum_H (x-1) ^{k(H) - c} (y-1)^{k(H) - |E(H)|-n}

where the sum ranges over all induced subgraphs `H` of `G`.

Functions
---------
"""

from contextlib import contextmanager
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.misc_c import prod
from sage.rings.integer_ring import ZZ
from sage.misc.decorators import sage_wraps

######################
# Graph Modification #
######################


@contextmanager
def removed_multiedge(G, unlabeled_edge):
    r"""
    A context manager which removes an edge with multiplicity from the
    graph `G` and restores it upon exiting.

    EXAMPLES::

        sage: from sage.graphs.tutte_polynomial import removed_multiedge
        sage: G = Graph(multiedges=True)
        sage: G.add_edges([(0,1,'a'),(0,1,'b')])
        sage: G.edges(sort=True)
        [(0, 1, 'a'), (0, 1, 'b')]
        sage: with removed_multiedge(G,(0,1)) as Y:
        ....:     G.edges(sort=True)
        []
        sage: G.edges(sort=True)
        [(0, 1, 'a'), (0, 1, 'b')]
    """
    u, v = unlabeled_edge
    edges = G.edge_boundary([u], [v], labels=True)
    G.delete_multiedge(u, v)
    try:
        yield
    finally:
        G.add_edges(edges)


@contextmanager
def removed_edge(G, edge):
    r"""
    A context manager which removes an edge from the graph `G` and
    restores it upon exiting.

    EXAMPLES::

        sage: from sage.graphs.tutte_polynomial import removed_edge
        sage: G = Graph()
        sage: G.add_edge(0,1)
        sage: G.edges(sort=True)
        [(0, 1, None)]
        sage: with removed_edge(G,(0,1)) as Y:
        ....:     G.edges(sort=True); G.vertices(sort=True)
        []
        [0, 1]
        sage: G.edges(sort=True)
        [(0, 1, None)]
    """
    G.delete_edge(edge)
    try:
        yield
    finally:
        G.add_edge(edge)


@contextmanager
def contracted_edge(G, unlabeled_edge):
    r"""
    Delete the first vertex in the edge, and make all the edges that
    went from it go to the second vertex.

    EXAMPLES::

        sage: from sage.graphs.tutte_polynomial import contracted_edge
        sage: G = Graph(multiedges=True)
        sage: G.add_edges([(0,1,'a'),(1,2,'b'),(0,3,'c')])
        sage: G.edges(sort=True)
        [(0, 1, 'a'), (0, 3, 'c'), (1, 2, 'b')]
        sage: with contracted_edge(G,(0,1)) as Y:
        ....:     G.edges(sort=True); G.vertices(sort=True)
        [(1, 2, 'b'), (1, 3, 'c')]
        [1, 2, 3]
        sage: G.edges(sort=True)
        [(0, 1, 'a'), (0, 3, 'c'), (1, 2, 'b')]
    """
    v1, v2 = unlabeled_edge
    loops = G.allows_loops()

    v1_edges = G.edges_incident(v1)
    G.delete_vertex(v1)
    added_edges = []

    for start, end, label in v1_edges:
        other_vertex = start if start != v1 else end
        edge = (other_vertex, v2, label)
        if loops or other_vertex != v2:
            G.add_edge(edge)
        added_edges.append(edge)

    try:
        yield
    finally:
        for edge in added_edges:
            G.delete_edge(edge)
        for edge in v1_edges:
            G.add_edge(edge)


@contextmanager
def removed_loops(G):
    r"""
    A context manager which removes all the loops in the graph `G`.
    It yields a list of the loops, and restores the loops upon
    exiting.

    EXAMPLES::

        sage: from sage.graphs.tutte_polynomial import removed_loops
        sage: G = Graph(multiedges=True, loops=True)
        sage: G.add_edges([(0,1,'a'),(1,2,'b'),(0,0,'c')])
        sage: G.edges(sort=True)
        [(0, 0, 'c'), (0, 1, 'a'), (1, 2, 'b')]
        sage: with removed_loops(G) as Y:
        ....:     G.edges(sort=True); G.vertices(sort=True); Y
        [(0, 1, 'a'), (1, 2, 'b')]
        [0, 1, 2]
        [(0, 0, 'c')]
        sage: G.edges(sort=True)
        [(0, 0, 'c'), (0, 1, 'a'), (1, 2, 'b')]
    """
    loops = G.loops()
    for edge in loops:
        G.delete_edge(edge)
    try:
        yield loops
    finally:
        for edge in loops:
            G.add_edge(edge)


def underlying_graph(G):
    r"""
    Given a graph `G` with multi-edges, returns a graph where all the
    multi-edges are replaced with a single edge.

    EXAMPLES::

        sage: from sage.graphs.tutte_polynomial import underlying_graph
        sage: G = Graph(multiedges=True)
        sage: G.add_edges([(0,1,'a'),(0,1,'b')])
        sage: G.edges(sort=True)
        [(0, 1, 'a'), (0, 1, 'b')]
        sage: underlying_graph(G).edges(sort=True)
        [(0, 1, None)]
    """
    from sage.graphs.graph import Graph
    g = Graph()
    g.allow_loops(True)
    for edge in set(G.edges(sort=False, labels=False)):
        g.add_edge(edge)
    return g


def edge_multiplicities(G):
    r"""
    Return the dictionary of multiplicities of the edges in the
    graph `G`.

    EXAMPLES::

        sage: from sage.graphs.tutte_polynomial import edge_multiplicities
        sage: G = Graph({1: [2,2,3], 2: [2], 3: [4,4], 4: [2,2,2]})
        sage: sorted(edge_multiplicities(G).items())
        [((1, 2), 2), ((1, 3), 1), ((2, 2), 1), ((2, 4), 3), ((3, 4), 2)]
    """
    d = {}
    for edge in G.edges(sort=False, labels=False):
        d[edge] = d.setdefault(edge, 0) + 1
    return d

########
# Ears #
########


class Ear():
    r"""
    An ear is a sequence of vertices

    Here is the definition from [HPR2010]_:

    An ear in a graph is a path `v_1 - v_2 - \dots - v_n - v_{n+1}`
    where `d(v_1) > 2`, `d(v_{n+1}) > 2` and
    `d(v_2) = d(v_3) = \dots = d(v_n) = 2`.

    A cycle is viewed as a special ear where `v_1 = v_{n+1}` and the
    restriction on the degree of this vertex is lifted.

    INPUT:
    """
    def __init__(self, graph, end_points, interior, is_cycle):
        """
        EXAMPLES::

            sage: G = graphs.PathGraph(4)
            sage: G.add_edges([(0,4),(0,5),(3,6),(3,7)])
            sage: from sage.graphs.tutte_polynomial import Ear
            sage: E = Ear(G,[0,3],[1,2],False)
        """
        self.end_points = end_points
        self.interior = interior
        self.is_cycle = is_cycle
        self.graph = graph

    @property
    def s(self):
        """
        Returns the number of distinct edges in this ear.

        EXAMPLES::

            sage: G = graphs.PathGraph(4)
            sage: G.add_edges([(0,4),(0,5),(3,6),(3,7)])
            sage: from sage.graphs.tutte_polynomial import Ear
            sage: E = Ear(G,[0,3],[1,2],False)
            sage: E.s
            3
        """
        return len(self.interior) + 1

    @property
    def vertices(self):
        """
        Returns the vertices of this ear.

        EXAMPLES::

            sage: G = graphs.PathGraph(4)
            sage: G.add_edges([(0,4),(0,5),(3,6),(3,7)])
            sage: from sage.graphs.tutte_polynomial import Ear
            sage: E = Ear(G,[0,3],[1,2],False)
            sage: E.vertices
            [0, 1, 2, 3]
        """
        return sorted(self.end_points + self.interior)

    @lazy_attribute
    def unlabeled_edges(self):
        """
        Returns the edges in this ear.

        EXAMPLES::

            sage: G = graphs.PathGraph(4)
            sage: G.add_edges([(0,4),(0,5),(3,6),(3,7)])
            sage: from sage.graphs.tutte_polynomial import Ear
            sage: E = Ear(G,[0,3],[1,2],False)
            sage: E.unlabeled_edges
            [(0, 1), (1, 2), (2, 3)]
        """
        return self.graph.edges_incident(vertices=self.interior, labels=False)

    @staticmethod
    def find_ear(g):
        """
        Finds the first ear in a graph.

        EXAMPLES::

            sage: G = graphs.PathGraph(4)
            sage: G.add_edges([(0,4),(0,5),(3,6),(3,7)])
            sage: from sage.graphs.tutte_polynomial import Ear
            sage: E = Ear.find_ear(G)
            sage: E.s
            3
            sage: E.unlabeled_edges
            [(0, 1), (1, 2), (2, 3)]
            sage: E.vertices
            [0, 1, 2, 3]
        """
        degree_two_vertices = [v for v, degree
                               in g.degree_iterator(labels=True)
                               if degree == 2]
        subgraph = g.subgraph(degree_two_vertices)
        for component in subgraph.connected_components(sort=False):
            edges = g.edges_incident(vertices=component, labels=True)
            all_vertices = sorted(set(sum([e[:2] for e in edges], ())))
            if len(all_vertices) < 3:
                continue
            end_points = [v for v in all_vertices if v not in component]
            if not end_points:
                end_points = [component[0]]

            ear_is_cycle = end_points[0] == end_points[-1]

            if ear_is_cycle:
                for e in end_points:
                    if e in component:
                        component.remove(e)

            return Ear(g, end_points, component, ear_is_cycle)

    @contextmanager
    def removed_from(self, G):
        r"""
        A context manager which removes the ear from the graph `G`.

        EXAMPLES::

            sage: G = graphs.PathGraph(4)
            sage: G.add_edges([(0,4),(0,5),(3,6),(3,7)])
            sage: len(G.edges(sort=True))
            7
            sage: from sage.graphs.tutte_polynomial import Ear
            sage: E = Ear.find_ear(G)
            sage: with E.removed_from(G) as Y:
            ....:     G.edges(sort=True)
            [(0, 4, None), (0, 5, None), (3, 6, None), (3, 7, None)]
            sage: len(G.edges(sort=True))
            7
        """
        deleted_edges = []
        for edge in G.edges_incident(vertices=self.interior, labels=True):
            G.delete_edge(edge)
            deleted_edges.append(edge)
        for v in self.interior:
            G.delete_vertex(v)

        try:
            yield
        finally:
            for edge in deleted_edges:
                G.add_edge(edge)

##################
# Edge Selection #
##################


class EdgeSelection():
    pass


class VertexOrder(EdgeSelection):
    def __init__(self, order):
        """
        EXAMPLES::

            sage: from sage.graphs.tutte_polynomial import VertexOrder
            sage: A = VertexOrder([4,6,3,2,1,7])
            sage: A.order
            [4, 6, 3, 2, 1, 7]
            sage: A.inverse_order
            {1: 4, 2: 3, 3: 2, 4: 0, 6: 1, 7: 5}
        """
        self.order = list(order)
        self.inverse_order = dict([reversed(_) for _ in enumerate(order)])

    def __call__(self, graph):
        """
        EXAMPLES::

            sage: from sage.graphs.tutte_polynomial import VertexOrder
            sage: A = VertexOrder([4,0,3,2,1,5])
            sage: G = graphs.PathGraph(6)
            sage: A(G)
            (3, 4, None)
        """
        for v in self.order:
            edges = graph.edges_incident([v])
            if edges:
                edges.sort(key=lambda x: self.inverse_order[x[0] if x[0] != v else x[1]])
                return edges[0]
        raise RuntimeError("no edges left to select")


class MinimizeSingleDegree(EdgeSelection):
    def __call__(self, graph):
        """
        EXAMPLES::

            sage: from sage.graphs.tutte_polynomial import MinimizeSingleDegree
            sage: G = graphs.PathGraph(6)
            sage: MinimizeSingleDegree()(G)
            (0, 1, None)
        """
        degrees = list(graph.degree_iterator(labels=True))
        degrees.sort(key=lambda x: x[1])  # Sort by degree
        for v, degree in degrees:
            for e in graph.edges_incident([v], labels=True):
                return e
        raise RuntimeError("no edges left to select")


class MinimizeDegree(EdgeSelection):
    def __call__(self, graph):
        """
        EXAMPLES::

            sage: from sage.graphs.tutte_polynomial import MinimizeDegree
            sage: G = graphs.PathGraph(6)
            sage: MinimizeDegree()(G)
            (0, 1, None)
        """
        degrees = dict(graph.degree_iterator(labels=True))
        edges = graph.edges(labels=True, sort=False)
        if edges:
            return min(edges, key=lambda x: degrees[x[0]] + degrees[x[1]])
        raise RuntimeError("no edges left to select")


class MaximizeDegree(EdgeSelection):
    def __call__(self, graph):
        """
        EXAMPLES::

            sage: from sage.graphs.tutte_polynomial import MaximizeDegree
            sage: G = graphs.PathGraph(6)
            sage: MaximizeDegree()(G)
            (1, 2, None)
        """
        degrees = dict(graph.degree_iterator(labels=True))
        edges = graph.edges(labels=True, sort=False)
        if edges:
            return max(edges, key=lambda x: degrees[x[0]] + degrees[x[1]])
        raise RuntimeError("no edges left to select")


###########
# Caching #
###########


def _cache_key(G):
    """
    Return the key used to cache the result for the graph G

    This is used by the decorator :func:`_cached`.

    EXAMPLES::

        sage: from sage.graphs.tutte_polynomial import _cache_key
        sage: G = graphs.DiamondGraph()
        sage: print(_cache_key(G))
        ((0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    """
    return tuple(G.canonical_label().edges(labels=False, sort=True))


def _cached(func):
    """
    Wrapper used to cache results of the function `func`

    This uses the function :func:`_cache_key`.

    EXAMPLES::

        sage: from sage.graphs.tutte_polynomial import tutte_polynomial
        sage: G = graphs.PetersenGraph()
        sage: T = tutte_polynomial(G)  #indirect doctest
        sage: tutte_polynomial(G)(1,1)  #indirect doctest
        2000
    """
    @sage_wraps(func)
    def wrapper(G, *args, **kwds):
        cache = kwds.setdefault('cache', {})
        key = _cache_key(G)
        if key in cache:
            return cache[key]
        result = func(G, *args, **kwds)
        cache[key] = result
        return result
    wrapper.original_func = func
    return wrapper


####################
# Tutte Polynomial #
####################

@_cached
def tutte_polynomial(G, edge_selector=None, cache=None):
    r"""
    Return the Tutte polynomial of the graph `G`.

    INPUT:

    - ``edge_selector`` (optional; method) this argument allows the user
      to specify his own heuristic for selecting edges used in the deletion
      contraction recurrence

    - ``cache`` -- (optional; dict) a dictionary to cache the Tutte
      polynomials generated in the recursive process.  One will be
      created automatically if not provided.

    EXAMPLES:

    The Tutte polynomial of any tree of order `n` is `x^{n-1}`::

        sage: all(T.tutte_polynomial() == x**9 for T in graphs.trees(10))               # needs sage.symbolic
        True

    The Tutte polynomial of the Petersen graph is::

        sage: P = graphs.PetersenGraph()
        sage: P.tutte_polynomial()
        x^9 + 6*x^8 + 21*x^7 + 56*x^6 + 12*x^5*y + y^6 + 114*x^5 + 70*x^4*y
        + 30*x^3*y^2 + 15*x^2*y^3 + 10*x*y^4 + 9*y^5 + 170*x^4 + 170*x^3*y
        + 105*x^2*y^2 + 65*x*y^3 + 35*y^4 + 180*x^3 + 240*x^2*y + 171*x*y^2
        + 75*y^3 + 120*x^2 + 168*x*y + 84*y^2 + 36*x + 36*y

    The Tutte polynomial of a connected graph `G` evaluated at (1,1) is the number of
    spanning trees of `G`::

        sage: G = graphs.RandomGNP(10,0.6)
        sage: while not G.is_connected():
        ....:     G = graphs.RandomGNP(10,0.6)
        sage: G.tutte_polynomial()(1,1) == G.spanning_trees_count()                     # needs sage.modules
        True

    Given that `T(x,y)` is the Tutte polynomial of a graph `G` with
    `n` vertices and `c` connected components, then `(-1)^{n-c} x^k
    T(1-x,0)` is the chromatic polynomial of `G`. ::

        sage: G = graphs.OctahedralGraph()
        sage: T = G.tutte_polynomial()
        sage: R = PolynomialRing(ZZ, 'x')
        sage: R((-1)^5*x*T(1-x,0)).factor()                                             # needs sage.symbolic
        (x - 2) * (x - 1) * x * (x^3 - 9*x^2 + 29*x - 32)
        sage: G.chromatic_polynomial().factor()
        (x - 2) * (x - 1) * x * (x^3 - 9*x^2 + 29*x - 32)

    TESTS:

    Providing an external cache::

        sage: cache = {}
        sage: _ = graphs.RandomGNP(7,.5).tutte_polynomial(cache=cache)
        sage: len(cache) > 0
        True

    Verify that :trac:`18366` is fixed::

        sage: g = Graph(multiedges=True)
        sage: g.add_edges([(0,1,1),(1,5,2),(5,3,3),(5,2,4),(2,4,5),(0,2,6),(0,3,7),(0,4,8),(0,5,9)])
        sage: g.tutte_polynomial()(1,1)
        52
        sage: g.spanning_trees_count()
        52
    """
    R = ZZ['x, y']
    if G.num_edges() == 0:
        return R.one()

    G = G.relabel(inplace=False, immutable=False)  # making sure the vertices are integers
    G.allow_loops(True)
    G.allow_multiple_edges(True)

    if edge_selector is None:
        edge_selector = MinimizeSingleDegree()
    x, y = R.gens()
    return _tutte_polynomial_internal(G, x, y, edge_selector, cache=cache)


@_cached
def _tutte_polynomial_internal(G, x, y, edge_selector, cache=None):
    """
    Does the recursive computation of the Tutte polynomial.

    INPUT:

    - ``G`` -- the graph
    - ``x,y`` -- the variables `x,y` respectively
    - ``edge_selector`` -- the heuristic for selecting edges used in the
      deletion contraction recurrence

    TESTS::

        sage: P = graphs.CycleGraph(5)
        sage: P.tutte_polynomial() # indirect doctest
        x^4 + x^3 + x^2 + x + y
    """
    if not G.num_edges():
        return x.parent().one()

    def recursive_tp(graph=None):
        """
        The recursive call -- used so that we do not have to specify
        the same arguments everywhere.
        """
        if graph is None:
            graph = G
        return _tutte_polynomial_internal(graph, x, y, edge_selector, cache=cache)

    # Remove loops
    with removed_loops(G) as loops:
        if loops:
            return y**len(loops) * recursive_tp()

    uG = underlying_graph(G)
    em = edge_multiplicities(G)
    d = list(em.values())

    def yy(start, end):
        return sum(y**i for i in range(start, end+1))

    # Lemma 1
    if G.is_forest():
        return prod(x + yy(1, d_i-1) for d_i in d)

    # Theorem 1: from Haggard, Pearce, Royle 2008
    blocks, cut_vertices = G.blocks_and_cut_vertices()
    if len(blocks) > 1:
        return prod([recursive_tp(G.subgraph(block)) for block in blocks])

    components = G.connected_components_number()
    edge = edge_selector(G)
    unlabeled_edge = edge[:2]

    with removed_edge(G, edge):
        if G.connected_components_number() > components:
            with contracted_edge(G, unlabeled_edge):
                return x*recursive_tp()

    ##################################
    # We are in the biconnected case #
    ##################################

    # Theorem 4: from Haggard, Pearce, and Royle Note that the formula
    # at http://homepages.ecs.vuw.ac.nz/~djp/files/TOMS10.pdf is
    # slightly incorrect.  The initial sum should only go to n-2
    # instead of n (allowing for the last part of the recursion).
    # Additionally, the first operand of the final product should be
    # (x+y^{1...(d_n+d_{n-1}-1)}) instead of just (x+y^(d_n+d_{n-1}-1)
    if uG.num_verts() == uG.num_edges():  # G is a multi-cycle
        n = len(d)
        result = 0
        for i in range(n - 2):
            term = (prod((x + yy(1, d_j-1)) for d_j in d[i+1:]) *
                    prod((yy(0, d_k-1)) for d_k in d[:i]))
            result += term
        # The last part of the recursion
        result += (x + yy(1, d[-1] + d[-2] - 1))*prod(yy(0, d_i-1)
                                                      for d_i in d[:-2])
        return result

    # Theorem 3 from Haggard, Pearce, and Royle, adapted to multi-ears
    ear = Ear.find_ear(uG)
    if ear is not None:
        if (ear.is_cycle and ear.vertices == G.vertices(sort=True)):
            # The graph is an ear (cycle) We should never be in this
            # case since we check for multi-cycles above
            return y + sum(x**i for i in range(1, ear.s))
        else:
            with ear.removed_from(G):
                # result = sum(x^i for i in range(ear.s)) #single ear case
                result = sum((prod(x + yy(1, em[e]-1) for e in ear.unlabeled_edges[i+1:])
                              * prod(yy(0, em[e]-1) for e in ear.unlabeled_edges[:i]))
                             for i in range(len(ear.unlabeled_edges)))
                result *= recursive_tp()

                with contracted_edge(G, [ear.end_points[0],
                                         ear.end_points[-1]]):
                    result += prod(yy(0, em[e]-1)
                                   for e in ear.unlabeled_edges)*recursive_tp()

            return result

    # Theorem 2
    if len(em) == 1:  # the graph is just a multiedge
        return x + sum(y**i for i in range(1, em[unlabeled_edge]))
    else:
        with removed_multiedge(G, unlabeled_edge):
            result = recursive_tp()
            with contracted_edge(G, unlabeled_edge):
                result += sum(y**i for i in range(em[unlabeled_edge]))*recursive_tp()
        return result
