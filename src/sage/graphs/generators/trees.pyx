r"""
Generation of trees

This module gathers methods related to the generation of trees.
All these methods appear in :mod:`sage.graphs.graph_generators`.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`BalancedTree` | Return the perfectly balanced tree of height `h \geq 1`, whose root has degree `r \geq 2`.
    :meth:`~sage.graphs.generators.degree_sequence.DegreeSequenceTree` | Return a tree with the given degree sequence.
    :meth:`FibonacciTree` | Return the graph of the Fibonacci Tree `F_{n}`.
    :meth:`RandomLobster` | Return a random lobster.
    :meth:`RandomTree` | Return a random tree on `n` nodes numbered `0` through `n-1`.
    :meth:`RandomTreePowerlaw` | Return a tree with a power law degree distribution, or ``False`` on failure.
    :meth:`trees` | Return a generator of the distinct trees on a fixed number of vertices.
    :meth:`nauty_gentreeg` | Return a generator which creates non-isomorphic trees from nauty's gentreeg program.

There are different ways to enumerate all non-isomorphic trees of order `n`.
We can use the implementation of the algorithm proposed in [WROM1986]_::

    sage: gen = graphs.trees(10)
    sage: T = next(gen); T
    Graph on 10 vertices
    sage: T.is_tree()
    True

We can also use nauty's gentreeg::

    sage: gen = graphs.nauty_gentreeg("10")
    sage: T = next(gen); T
    Graph on 10 vertices
    sage: T.is_tree()
    True

Note that nauty's gentreeg can only be used to generate trees with
`1 \leq n \leq 128` vertices::

    sage: next(graphs.nauty_gentreeg("0"))
    Traceback (most recent call last):
    ...
    ValueError: wrong format of parameter options
    sage: next(graphs.nauty_gentreeg("1"))
    Graph on 1 vertex
    sage: next(graphs.nauty_gentreeg("128"))
    Graph on 128 vertices
    sage: next(graphs.nauty_gentreeg("129"))
    Traceback (most recent call last):
    ...
    ValueError: wrong format of parameter options

We don't have this limitation with method :meth:`trees`::

    sage: next(graphs.trees(0))
    Graph on 0 vertices
    sage: next(graphs.trees(129))
    Graph on 129 vertices
    sage: next(graphs.trees(1000))
    Graph on 1000 vertices

Nauty's gentreeg can be used to generate trees with bounded maximum degree::

    sage: gen = graphs.nauty_gentreeg("8 -D3")
    sage: all(max(g.degree()) <= 3 for g in gen)
    True
    sage: len(list(graphs.nauty_gentreeg("8 -D2")))
    1

Nauty's gentreeg can be used to generate trees with bounded diameter::

    sage: all(g.diameter() == 3 for g in graphs.nauty_gentreeg("8 -Z3"))
    True
    sage: len(list(graphs.nauty_gentreeg("8 -Z3")))
    3

The number of trees on the first few vertex counts agrees with :oeis:`A000055`::

    sage: [len(list(graphs.trees(i))) for i in range(15)]
    [1, 1, 1, 1, 2, 3, 6, 11, 23, 47, 106, 235, 551, 1301, 3159]
    sage: [len(list(graphs.nauty_gentreeg(str(i)))) for i in range(1, 15)]
    [1, 1, 1, 2, 3, 6, 11, 23, 47, 106, 235, 551, 1301, 3159]

Methods
-------
"""

# ****************************************************************************
#       Copyright (C) 2009 Ryan Dingman
#                     2010 Harald Schilly
#                     2010 Yann Laigle-Chapuy
#                     2010 Edward Scheinerman
#                     2025 David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sys
from libc.limits cimport INT_MAX
from cysignals.memory cimport check_allocarray, sig_free
from sage.graphs.base.sparse_graph cimport SparseGraph
from sage.graphs.base.sparse_graph cimport SparseGraphBackend

from sage.graphs.graph import Graph
from sage.misc.randstate import current_randstate
from sage.misc.randstate import set_random_seed
from sage.misc.prandom import randint


def BalancedTree(r, h):
    r"""
    Return the perfectly balanced tree of height `h \geq 1`,
    whose root has degree `r \geq 2`.

    The number of vertices of this graph is `1 + r + r^2 + \cdots + r^h`, that
    is, `\frac{r^{h+1} - 1}{r - 1}`. The number of edges is one less than the
    number of vertices.

    INPUT:

    - ``r`` -- positive integer `\geq 2`; the degree of the root node

    - ``h`` -- positive integer `\geq 1`; the height of the balanced tree

    OUTPUT:

    The perfectly balanced tree of height `h \geq 1` and whose root has
    degree `r \geq 2`.

    EXAMPLES:

    A balanced tree whose root node has degree `r = 2`, and of height
    `h = 1`, has order 3 and size 2::

        sage: G = graphs.BalancedTree(2, 1); G
        Balanced tree: Graph on 3 vertices
        sage: G.order()
        3
        sage: G.size()
        2
        sage: r = 2; h = 1
        sage: v = 1 + r
        sage: v; v - 1
        3
        2

    Plot a balanced tree of height 5, whose root node has degree `r = 3`::

        sage: G = graphs.BalancedTree(3, 5)
        sage: G.plot()                          # long time                             # needs sage.plot
        Graphics object consisting of 728 graphics primitives

    A tree is bipartite. If its vertex set is finite, then it is planar::

        sage: # needs networkx
        sage: r = randint(2, 5); h = randint(1, 7)
        sage: T = graphs.BalancedTree(r, h)
        sage: T.is_bipartite()
        True
        sage: T.is_planar()
        True
        sage: v = (r^(h + 1) - 1) / (r - 1)
        sage: T.order() == v
        True
        sage: T.size() == v - 1
        True

    TESTS:

    Normally we would only consider balanced trees whose root node
    has degree `r \geq 2`, but the construction degenerates
    gracefully::

        sage: graphs.BalancedTree(1, 10)
        Balanced tree: Graph on 11 vertices

    Similarly, we usually want the tree must have height `h \geq 1`
    but the algorithm also degenerates gracefully here::

        sage: graphs.BalancedTree(3, 0)
        Balanced tree: Graph on 1 vertex

    The construction is the same as the one of networkx::

        sage: # needs networkx
        sage: import networkx
        sage: r = randint(2, 4); h = randint(1, 5)
        sage: T = graphs.BalancedTree(r, h)
        sage: N = Graph(networkx.balanced_tree(r, h), name="Balanced tree")
        sage: T.is_isomorphic(N)
        True
    """
    # Compute the number of vertices per level of the tree
    order = [r**l for l in range(h + 1)]
    # Compute the first index of the vertices of a level
    begin = [0]
    begin.extend(begin[-1] + val for val in order)
    # The number of vertices of the tree is the first index of level h + 1
    T = Graph(begin[-1], name="Balanced tree")

    # Add edges of the r-ary tree
    for level in range(h):
        start = begin[level + 1]
        for u in range(begin[level], begin[level + 1]):
            T.add_edges((u, v) for v in range(start, start + r))
            start += r
    return T


def FibonacciTree(n):
    r"""
    Return the graph of the Fibonacci Tree `F_{n}`.

    The Fibonacci tree `F_{n}` is recursively defined as the tree
    with a root vertex and two attached child trees `F_{n-1}` and
    `F_{n-2}`, where `F_{1}` is just one vertex and `F_{0}` is empty.

    INPUT:

    - ``n`` -- the recursion depth of the Fibonacci Tree

    EXAMPLES::

        sage: g = graphs.FibonacciTree(3)                                               # needs sage.libs.pari
        sage: g.is_tree()                                                               # needs sage.libs.pari
        True

    ::

        sage: l1 = [graphs.FibonacciTree(_).order() + 1 for _ in range(6)]              # needs sage.libs.pari
        sage: l2 = list(fibonacci_sequence(2,8))                                        # needs sage.libs.pari
        sage: l1 == l2                                                                  # needs sage.libs.pari
        True

    AUTHORS:

    - Harald Schilly and Yann Laigle-Chapuy (2010-03-25)
    """
    T = Graph(name="Fibonacci-Tree-%d" % n)
    if n == 1:
        T.add_vertex(0)
    if n < 2:
        return T

    from sage.combinat.combinat import fibonacci_sequence
    F = list(fibonacci_sequence(n + 2))
    s = 1.618 ** (n / 1.618 - 1.618)
    pos = {}

    def fib(level, node, y):
        pos[node] = (node, y)
        if level < 2:
            return
        level -= 1
        y -= s
        diff = F[level]
        T.add_edge(node, node - diff)
        if level == 1:  # only one child
            pos[node - diff] = (node, y)
            return
        T.add_edge(node, node + diff)
        fib(level, node - diff, y)
        fib(level - 1, node + diff, y)

    T.add_vertices(range(sum(F[:-1])))
    fib(n, F[n + 1] - 1, 0)
    T.set_pos(pos)

    return T


def Caterpillar(spine):
    r"""
    Return the caterpillar tree with given spine sequence.

    A caterpillar tree consists of leaves attached to a path (the "spine").

    INPUT:

    - ``spine`` -- list of nonnegative integers in the form
       `[a_1, a_2, \dots, a_n]`, where `a_i` is the number of leaves adjacent
       to the `i`-th vertex on the spine (except for the first and last vertex,
       which have `a_1 + 1` and `a_n + 1` leaf-neighbors, respectively)

    OUTPUT:

    A caterpillar tree of diameter `n+1` on `n + 2 + \sum_{i=1}^n a_i` vertices,
    `n` of which are not leaves.

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm if the returned graph does not have too many
    vertices. The spine vertices are positioned on a straight line together
    with two leaves at its ends. Every edge in the drawing has unit length.

    EXAMPLES:

    Caterpillars with all-zero spine sequence are paths::

        sage: graphs.Caterpillar([]).is_isomorphic(graphs.PathGraph(2))
        True
        sage: graphs.Caterpillar([0]).is_isomorphic(graphs.PathGraph(3))
        True
        sage: graphs.Caterpillar([0, 0]).is_isomorphic(graphs.PathGraph(4))
        True

    Caterpillars with singleton spine are stars::

        sage: graphs.Caterpillar([1]).is_isomorphic(graphs.StarGraph(3))
        True
        sage: graphs.Caterpillar([2]).is_isomorphic(graphs.StarGraph(4))
        True
        sage: graphs.Caterpillar([3]).is_isomorphic(graphs.StarGraph(5))
        True

    Distinct spine sequences can yield isomorphic caterpillars::

        sage: graphs.Caterpillar([1,1,2]).is_isomorphic(graphs.Caterpillar([2,1,1]))
        True

    TESTS:

    Generated graphs have diameter ``len(spine) + 1``::

        sage: graphs.Caterpillar([7]).diameter()
        2
        sage: graphs.Caterpillar([2,2,2,2]).diameter()
        5
        sage: graphs.Caterpillar([0,1,1,0]).diameter()
        5
    """
    spine = list(spine)
    cdef int spine_len = len(spine)
    cdef int n_vertices = spine_len + 2 + sum(spine)
    T = Graph(n_vertices, name=f"Caterpillar({','.join(map(str, spine))})")

    # add spine
    for i in range(spine_len - 1):
        T._backend.add_edge(i, i + 1, None, False)

    # add a leaf at both ends of the spine
    T._backend.add_edge(spine_len + 1, 0, None, False)
    if spine:
        T._backend.add_edge(spine_len - 1, spine_len, None, False)

    # add leaves
    cdef int v = spine_len + 2
    for i, d in enumerate(spine):
        for j in range(d):
            T._backend.add_edge(i, v + j, None, False)
        v += d

    # add embedding
    cdef int max_leaves = max(spine, default=0)
    if (spine_len < 10 and max_leaves < 3) or (spine_len < 6 and max_leaves < 7):
        T._pos = {spine_len + 1: (-1, 0), spine_len: (spine_len, 0)}
        radius = 0.3
        v = spine_len + 2
        for x, d in enumerate(spine):
            T._pos[x] = (x, 0)
            mid = v + d // 2
            T._line_embedding(range(v, mid), first=(x - radius, 1), last=(x + radius, 1))
            T._line_embedding(range(mid, v + d), first=(x - radius, -1), last=(x + radius, -1))
            v += d

    return T


def RandomLobster(n, p, q, seed=None):
    r"""
    Return a random lobster.

    A lobster is a tree that reduces to a caterpillar when pruning all
    leaf vertices. A caterpillar is a tree that reduces to a path when
    pruning all leaf vertices (`q=0`).

    INPUT:

    - ``n`` -- expected number of vertices in the backbone

    - ``p`` -- probability of adding an edge to the backbone

    - ``q`` -- probability of adding an edge (claw) to the arms

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    EXAMPLES:

    We check a random graph with 12 backbone
    nodes and probabilities `p = 0.7` and `q = 0.3`::

        sage: # needs networkx
        sage: G = graphs.RandomLobster(12, 0.7, 0.3)
        sage: leaves = [v for v in G.vertices(sort=False) if G.degree(v) == 1]
        sage: G.delete_vertices(leaves)                                 # caterpillar
        sage: leaves = [v for v in G.vertices(sort=False) if G.degree(v) == 1]
        sage: G.delete_vertices(leaves)                                 # path
        sage: s = G.degree_sequence()
        sage: if G:
        ....:     if G.n_vertices() == 1:
        ....:         assert s == [0]
        ....:     else:
        ....:         assert s[-2:] == [1, 1]
        ....:     assert all(d == 2 for d in s[:-2])

    ::

        sage: G = graphs.RandomLobster(9, .6, .3)                                       # needs networkx
        sage: G.show()                          # long time                             # needs networkx sage.plot
    """
    if seed is None:
        seed = int(current_randstate().long_seed() % sys.maxsize)
    import networkx
    return Graph(networkx.random_lobster(n, p, q, seed=seed))


def RandomTree(n, seed=None):
    r"""
    Return a random tree on `n` nodes numbered `0` through `n-1`.

    By Cayley's theorem, there are `n^{n-2}` trees with vertex
    set `\{0,1,\dots,n-1\}`. This constructor chooses one of these uniformly
    at random.

    ALGORITHM:

    The algorithm works by generating an `(n-2)`-long
    random sequence of numbers chosen independently and uniformly
    from `\{0,1,\dots,n-1\}` and then applies an inverse
    Prufer transformation.

    INPUT:

    - ``n`` -- number of vertices in the tree

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    EXAMPLES::

        sage: G = graphs.RandomTree(10)
        sage: G.is_tree()
        True
        sage: G.show()                          # long time                             # needs sage.plot

    TESTS:

    Ensuring that we encounter no unexpected surprise ::

        sage: all( graphs.RandomTree(10).is_tree()
        ....:      for i in range(100) )
        True

    Random tree with one and zero vertices::

        sage: graphs.RandomTree(0)
        Graph on 0 vertices
        sage: graphs.RandomTree(1)
        Graph on 1 vertex
    """
    g = Graph(n)
    if n <= 1:
        return g

    if seed is not None:
        set_random_seed(seed)

    # create random Prufer code
    code = [randint(0, n - 1) for i in range(n - 2)]

    # We count the number of symbols of each type.
    # count[k] is the number of times k appears in code
    #
    # (count[k] is set to -1 when the corresponding vertex is not
    # available anymore)
    count = [0] * n
    for k in code:
        count[k] += 1

    # We use a heap to store vertices for which count[k] == 0 and get the vertex
    # with smallest index
    from heapq import heapify, heappop, heappush
    zeros = [x for x in range(n) if not count[x]]
    heapify(zeros)

    for s in code:
        x = heappop(zeros)
        g.add_edge(x, s)
        count[x] = -1
        count[s] -= 1
        if not count[s]:
            heappush(zeros, s)

    # Adding as an edge the last two available vertices
    g.add_edge(zeros)

    return g


def RandomTreePowerlaw(n, gamma=3, tries=1000, seed=None):
    """
    Return a tree with a power law degree distribution, or ``False`` on failure.

    From the NetworkX documentation: a trial power law degree sequence is chosen
    and then elements are swapped with new elements from a power law
    distribution until the sequence makes a tree (size = order - 1).

    INPUT:

    - ``n`` -- number of vertices

    - ``gamma`` -- exponent of power law distribution

    - ``tries`` -- number of attempts to adjust sequence to make a tree

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    EXAMPLES:

    We check that the generated graph is a tree::

        sage: G = graphs.RandomTreePowerlaw(10, 3)                                      # needs networkx
        sage: G.is_tree()                                                               # needs networkx
        True
        sage: G.order(), G.size()                                                       # needs networkx
        (10, 9)

    ::

        sage: G = graphs.RandomTreePowerlaw(15, 2)                                      # needs networkx
        sage: if G:                             # random output         # long time, needs networkx sage.plot
        ....:     G.show()
    """
    if seed is None:
        seed = int(current_randstate().long_seed() % sys.maxsize)
    import networkx
    try:
        return Graph(networkx.random_powerlaw_tree(n, gamma, seed=seed, tries=tries))
    except networkx.NetworkXError:
        return False


def trees(n):
    r"""
    Return a generator of the distinct trees on a fixed number of vertices.

    INPUT:

    - ``n`` -- integer; the size of the trees created

    OUTPUT:

    A generator which creates an exhaustive, duplicate-free listing of the
    connected free (unlabeled) trees with ``n`` vertices. A tree is a graph
    with no cycles.

    ALGORITHM:

    Uses the algorithm that generates each new tree in constant time described
    in [WROM1986]_. See the documentation for, and implementation of, the
    :mod:`sage.graphs.generators.trees` module.

    EXAMPLES:

    We create an iterator, then loop over its elements::

        sage: tree_iterator = graphs.trees(7)
        sage: for T in tree_iterator:
        ....:     print(T.degree_sequence())
        [2, 2, 2, 2, 2, 1, 1]
        [3, 2, 2, 2, 1, 1, 1]
        [3, 2, 2, 2, 1, 1, 1]
        [4, 2, 2, 1, 1, 1, 1]
        [3, 3, 2, 1, 1, 1, 1]
        [3, 3, 2, 1, 1, 1, 1]
        [4, 3, 1, 1, 1, 1, 1]
        [3, 2, 2, 2, 1, 1, 1]
        [4, 2, 2, 1, 1, 1, 1]
        [5, 2, 1, 1, 1, 1, 1]
        [6, 1, 1, 1, 1, 1, 1]

    The number of trees on the first few vertex counts.
    This is sequence A000055 in Sloane's OEIS::

        sage: [len(list(graphs.trees(i))) for i in range(15)]
        [1, 1, 1, 1, 2, 3, 6, 11, 23, 47, 106, 235, 551, 1301, 3159]
    """
    from sage.graphs.generators.trees import TreeIterator
    return iter(TreeIterator(n))


cdef class TreeIterator:
    r"""
    This class iterates over all trees with n vertices (up to isomorphism).

    EXAMPLES::

        sage: from sage.graphs.generators.trees import TreeIterator
        sage: def check_trees(n):
        ....:     trees = []
        ....:     for t in TreeIterator(n):
        ....:         if not t.is_tree():
        ....:             return False
        ....:         if t.n_vertices() != n:
        ....:             return False
        ....:         if t.n_edges() != n - 1:
        ....:             return False
        ....:         for tree in trees:
        ....:             if tree.is_isomorphic(t):
        ....:                 return False
        ....:         trees.append(t)
        ....:     return True
        sage: check_trees(10)
        True

    ::

        sage: from sage.graphs.generators.trees import TreeIterator
        sage: count = 0
        sage: for t in TreeIterator(15):
        ....:     count += 1
        sage: count
        7741
    """

    def __init__(self, int n):
        r"""
        Initialize an iterator over all trees with `n` vertices.

        EXAMPLES::

            sage: from sage.graphs.generators.trees import TreeIterator
            sage: t = TreeIterator(100) # indirect doctest
            sage: print(t)
            Iterator over all trees with 100 vertices
        """
        self.n = n
        self.l = NULL
        self.current_level_sequence = NULL
        self.first_time = 1

    def __dealloc__(self):
        r"""
        EXAMPLES::
            sage: from sage.graphs.generators.trees import TreeIterator
            sage: t = TreeIterator(100)
            sage: t = None # indirect doctest
        """
        sig_free(self.l)
        sig_free(self.current_level_sequence)

    def __str__(self) -> str:
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.graphs.generators.trees import TreeIterator
            sage: t = TreeIterator(100)
            sage: print(t)  # indirect doctest
            Iterator over all trees with 100 vertices
        """
        return "Iterator over all trees with %s vertices" % (self.n)

    def __iter__(self):
        r"""
        Return an iterator over all the trees with `n` vertices.

        EXAMPLES::

            sage: from sage.graphs.generators.trees import TreeIterator
            sage: t = TreeIterator(4)
            sage: list(iter(t))
            [Graph on 4 vertices, Graph on 4 vertices]
        """
        return self

    def __next__(self):
        r"""
        Return the next tree with `n` vertices.

        EXAMPLES::

            sage: from sage.graphs.generators.trees import TreeIterator
            sage: T = TreeIterator(5)
            sage: [t for t in T] # indirect doctest
            [Graph on 5 vertices, Graph on 5 vertices, Graph on 5 vertices]

        TESTS:

        This used to be broken for trees with no vertices
        and was fixed in :issue:`13719`::

            sage: from sage.graphs.generators.trees import TreeIterator
            sage: T = TreeIterator(0)
            sage: [t for t in T] # indirect doctest
            [Graph on 0 vertices]
        """
        if not self.first_time and not self.q:
            raise StopIteration

        if self.first_time == 1:
            self.first_time = 0
            if self.n:
                self.l = <int *>check_allocarray(self.n, sizeof(int))
                self.current_level_sequence = <int *>check_allocarray(self.n, sizeof(int))

                self.generate_first_level_sequence()
            else:
                self.q = 0
        else:
            self.generate_next_level_sequence()

        cdef int i
        cdef int vertex1
        cdef int vertex2
        cdef object G

        G = Graph(self.n, sparse=True)
        cdef SparseGraph SG = (<SparseGraphBackend?> G._backend)._cg

        for i in range(2, self.n + 1):
            vertex1 = i - 1
            vertex2 = self.current_level_sequence[i - 1] - 1
            SG.add_arc_unsafe(vertex1, vertex2)

        return G

    cdef int generate_first_level_sequence(self) noexcept:
        r"""
        Generates the level sequence representing the first tree with `n` vertices
        """
        cdef int i
        cdef int k

        k = (self.n // 2) + 1

        if self.n == 4:
            self.p = 3
        else:
            self.p = self.n
        self.q = self.n - 1
        self.h1 = k
        self.h2 = self.n
        if self.n % 2:
            self.c = INT_MAX  # oo
        else:
            self.c = self.n + 1

        self.r = k

        for i in range(1, k + 1):
            self.l[i - 1] = i
        for i in range(k + 1, self.n + 1):
            self.l[i - 1] = i - k + 1
        for i in range(self.n):
            self.current_level_sequence[i] = i
        if self.n > 2:
            self.current_level_sequence[k] = 1
        if self.n <= 3:
            self.q = 0

        return 0

    cdef int generate_next_level_sequence(self) noexcept:
        r"""
        Generates the level sequence representing the next tree with `n` vertices
        """
        cdef int i
        cdef int fixit = 0

        cdef int needr = 0
        cdef int needc = 0
        cdef int needh2 = 0

        cdef int n = self.n
        cdef int p = self.p
        cdef int q = self.q
        cdef int h1 = self.h1
        cdef int h2 = self.h2
        cdef int c = self.c
        cdef int r = self.r
        cdef int *l = self.l
        cdef int *w = self.current_level_sequence

        if c == n + 1 or p == h2 and (l[h1 - 1] == l[h2 - 1] + 1 and n - h2 > r - h1 or l[h1 - 1] == l[h2 - 1] and n - h2 + 1 < r - h1):
            if l[r - 1] > 3:
                p = r
                q = w[r - 1]
                if h1 == r:
                    h1 = h1 - 1
                fixit = 1
            else:
                p = r
                r = r - 1
                q = 2

        if p <= h1:
            h1 = p - 1
        if p <= r:
            needr = 1
        elif p <= h2:
            needh2 = 1
        elif l[h2 - 1] == l[h1 - 1] - 1 and n - h2 == r - h1:
            if p <= c:
                needc = 1
        else:
            c = INT_MAX

        cdef int oldp = p
        cdef int delta = q - p
        cdef int oldlq = l[q - 1]
        cdef int oldwq = w[q - 1]
        p = INT_MAX

        for i in range(oldp, n + 1):
            l[i - 1] = l[i - 1 + delta]
            if l[i - 1] == 2:
                w[i - 1] = 1
            else:
                p = i
                if l[i - 1] == oldlq:
                    q = oldwq
                else:
                    q = w[i - 1 + delta] - delta
                w[i - 1] = q
            if needr == 1 and l[i - 1] == 2:
                needr = 0
                needh2 = 1
                r = i - 1
            if needh2 == 1 and l[i - 1] <= l[i - 2] and i > r + 1:
                needh2 = 0
                h2 = i - 1
                if l[h2 - 1] == l[h1 - 1] - 1 and n - h2 == r - h1:
                    needc = 1
                else:
                    c = INT_MAX
            if needc == 1:
                if l[i - 1] != l[h1 - h2 + i - 1] - 1:
                    needc = 0
                    c = i
                else:
                    c = i + 1

        if fixit == 1:
            r = n - h1 + 1
            for i in range(r + 1, n + 1):
                l[i - 1] = i - r + 1
                w[i - 1] = i - 1
            w[r] = 1
            h2 = n
            p = n
            q = p - 1
            c = INT_MAX
        else:
            if p == INT_MAX:
                if l[oldp - 2] != 2:
                    p = oldp - 1
                else:
                    p = oldp - 2
                q = w[p - 1]
            if needh2 == 1:
                h2 = n
                if l[h2 - 1] == l[h1 - 1] - 1 and h1 == r:
                    c = n + 1
                else:
                    c = INT_MAX

        self.p = p
        self.q = q
        self.h1 = h1
        self.h2 = h2
        self.c = c
        self.r = r
        self.l = l
        self.current_level_sequence = w

        return 0


def nauty_gentreeg(options='', debug=False):
    r"""
    Return a generator which creates non-isomorphic trees from nauty's gentreeg
    program.

    INPUT:

    - ``options`` -- string (default: ``""``); a string passed to ``gentreeg``
      as if it was run at a system command line. At a minimum, you *must* pass
      the number of vertices you desire. Sage expects the graphs to be in
      nauty's "sparse6" format, do not set an option to change this default or
      results will be unpredictable.

    - ``debug`` -- boolean (default: ``False``); if ``True`` the first line of
      ``gentreeg``'s output to standard error is captured and the first call to
      the generator's ``next()`` function will return this line as a string. A
      line leading with ">A" indicates a successful initiation of the program
      with some information on the arguments, while a line beginning with ">E"
      indicates an error with the input.

    The possible options, obtained as output of ``gentreeg -help``::

           n            : the number of vertices. Must be in range 1..128
        res/mod         : only generate subset res out of subsets 0..mod-1
          -D<int>       : an upper bound for the maximum degree
          -Z<int>:<int> : bounds on the diameter
          -q            : suppress auxiliary output

    Options which cause ``gentreeg`` to use an output format different than the
    sparse6 format are not listed above (-p, -l, -u) as they will confuse the
    creation of a Sage graph. The res/mod option can be useful when using the
    output in a routine run several times in parallel.

    OUTPUT:

    A generator which will produce the graphs as Sage graphs. These will be
    simple graphs: no loops, no multiple edges, no directed edges.

    .. SEEALSO::

        :meth:`trees` -- another generator of trees

    EXAMPLES:

    The generator can be used to construct trees for testing, one at a time
    (usually inside a loop). Or it can be used to create an entire list all at
    once if there is sufficient memory to contain it::

        sage: gen = graphs.nauty_gentreeg("4")
        sage: next(gen)
        Graph on 4 vertices
        sage: next(gen)
        Graph on 4 vertices
        sage: next(gen)
        Traceback (most recent call last):
        ...
        StopIteration

    The number of trees on the first few vertex counts. This agrees with
    :oeis:`A000055`::

        sage: [len(list(graphs.nauty_gentreeg(str(i)))) for i in range(1, 15)]
        [1, 1, 1, 2, 3, 6, 11, 23, 47, 106, 235, 551, 1301, 3159]

    The ``debug`` switch can be used to examine ``gentreeg``'s reaction to the
    input in the ``options`` string.  We illustrate success. (A failure will be
    a string beginning with ">E".)  Passing the "-q" switch to ``gentreeg`` will
    suppress the indicator of a successful initiation, and so the first returned
    value might be an empty string if ``debug`` is ``True``::

        sage: gen = graphs.nauty_gentreeg("4", debug=True)
        sage: print(next(gen))
        >A ...gentreeg ...
        sage: gen = graphs.nauty_gentreeg("4 -q", debug=True)
        sage: next(gen)
        ''

    TESTS:

    The number `n` of vertices must be in range 1..128::

        sage: list(graphs.nauty_gentreeg("0", debug=False))
        Traceback (most recent call last):
        ...
        ValueError: wrong format of parameter options
        sage: list(graphs.nauty_gentreeg("0", debug=True))
        ['>E gentreeg: n must be in the range 1..128\n']
        sage: list(graphs.nauty_gentreeg("200", debug=True))
        ['>E gentreeg: n must be in the range 1..128\n']

    Wrong input::

        sage: list(graphs.nauty_gentreeg("3 -x", debug=False))
        Traceback (most recent call last):
        ...
        ValueError: wrong format of parameter options
        sage: list(graphs.nauty_gentreeg("3 -x", debug=True))
        ['>E Usage: ...gentreeg [-D#] [-Z#:#] [-ulps] [-q] n... [res/mod] ...
        sage: list(graphs.nauty_gentreeg("3", debug=True))
        ['>A ...gentreeg ...\n', Graph on 3 vertices]
    """
    import shlex
    import subprocess
    from sage.features.nauty import NautyExecutable
    gen_path = NautyExecutable("gentreeg").absolute_filename()
    with subprocess.Popen(shlex.quote(gen_path) + " {0}".format(options), shell=True,
                          stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, close_fds=True,
                          encoding='latin-1') as sp:
        msg = sp.stderr.readline()
        if debug:
            yield msg
        elif msg.startswith('>E'):
            raise ValueError('wrong format of parameter options')
        gen = sp.stdout
        while True:
            try:
                s = next(gen)
            except StopIteration:
                # Exhausted list of graphs from nauty geng
                return
            G = Graph(s[:-1], format='sparse6', loops=False, multiedges=False)
            yield G
