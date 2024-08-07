# cython: binding=True
r"""
Tree decompositions

This module implements tree-decomposition methods.

A tree-decomposition of a graph `G = (V, E)` is a pair `(X, T)`, where `X=\{X_1,
X_2, \ldots, X_t\}` is a family of subsets of `V`, usually called *bags*, and
`T` is a tree of order `t` whose nodes are the subsets `X_i` satisfying the
following properties:

- The union of all sets `X_i` equals `V`. That is, each vertex of the graph `G`
  is associated with at least one tree node.

- For every edge `(v, w)` in the graph, there is a subset `X_i` that contains
  both `v` and `w`. That is, each edge of the graph `G` appears in a tree node.

- The nodes associated with vertex `v \in V` form a connected subtree of
  `T`. That is, if `X_i` and `X_j` both contain a vertex `v \in V`, then all
  nodes `X_k` of the tree in the (unique) path between `X_i` and `X_j` contain
  `v` as well, and we have `X_i \cap X_j \subseteq X_k`.

The *width* of a tree decomposition is the size of the largest set `X_i` minus
one, i.e., `\max_{X_i \in X} |X_i| - 1`, and the *treewidth* `tw(G)` of a graph
`G` is the minimum width among all possible tree decompositions of `G`. Observe
that, the size of the largest set is diminished by one in order to make the
treewidth of a tree equal to one.

The *length* of a tree decomposition, as proposed in [DG2006]_, is the maximum
*diameter* in `G` of its bags, where the diameter of a bag `X_i` is the largest
distance in `G` between the vertices in `X_i` (i.e., `\max_{u, v \in X_i}
\dist_G(u, v)`). The *treelength* `tl(G)` of a graph `G` is the minimum length
among all possible tree decompositions of `G`.

While deciding whether a graph has treelength 1 can be done in linear time
(equivalent to deciding if the graph is chordal), deciding if it has treelength
at most `k` for any fixed constant `k \leq 2` is NP-complete [Lokshtanov2009]_.

Treewidth and treelength are different measures of tree-likeness. In particular,
trees have treewidth and treelength 1::

    sage: T = graphs.RandomTree(20)
    sage: T.treewidth()
    1
    sage: T.treelength()
    1

The treewidth of a cycle is 2 and its treelength is `\lceil n/3 \rceil`::

    sage: [graphs.CycleGraph(n).treewidth() for n in range(3, 11)]
    [2, 2, 2, 2, 2, 2, 2, 2]
    sage: [graphs.CycleGraph(n).treelength() for n in range(3, 11)]
    [1, 2, 2, 2, 3, 3, 3, 4]

The treewidth of a clique is `n-1` and its treelength is 1::

    sage: [graphs.CompleteGraph(n).treewidth() for n in range(3, 11)]
    [2, 3, 4, 5, 6, 7, 8, 9]
    sage: [graphs.CompleteGraph(n).treelength() for n in range(3, 11)]
    [1, 1, 1, 1, 1, 1, 1, 1]


.. SEEALSO::

    - :wikipedia:`Tree_decomposition`
    - :wikipedia:`Treewidth`


**This module contains the following methods**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`treewidth` | Compute the treewidth of `G` (and provide a decomposition).
    :meth:`treelength` | Compute the treelength of `G` (and provide a decomposition).
    :meth:`make_nice_tree_decomposition` | Return a *nice* tree decomposition (TD) of the TD `tree_decomp`.
    :meth:`label_nice_tree_decomposition` | Return a nice tree decomposition with nodes labelled accordingly.
    :meth:`is_valid_tree_decomposition` | Check whether `T` is a valid tree-decomposition for `G`.
    :meth:`reduced_tree_decomposition` | Return a reduced tree-decomposition of `T`.
    :meth:`width_of_tree_decomposition` | Return the width of the tree decomposition `T` of `G`.
    :meth:`length_of_tree_decomposition` | Return the length of the tree decomposition `T` of `G`.


.. TODO:

    - Approximation of treelength based on :meth:`~sage.graphs.graph.Graph.lex_M`
    - Approximation of treelength based on BFS Layering


Methods
-------
"""
# ****************************************************************************
#       Copyright (C) 2020 David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.sets.set import Set
from sage.misc.cachefunc import cached_function
from sage.sets.disjoint_set import DisjointSet
from sage.rings.infinity import Infinity
from sage.graphs.distances_all_pairs cimport c_distances_all_pairs
from cysignals.memory cimport sig_calloc, sig_free

from itertools import chain
from collections.abc import Iterable


def is_valid_tree_decomposition(G, T):
    r"""
    Check whether `T` is a valid tree-decomposition for `G`.

    INPUT:

    - ``G`` -- a sage Graph

    - ``T`` -- a tree decomposition, i.e., a tree whose vertices are the bags
      (subsets of vertices) of the decomposition

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import is_valid_tree_decomposition
        sage: K = graphs.CompleteGraph(4)
        sage: T = Graph()
        sage: T.add_vertex(Set(K))
        sage: is_valid_tree_decomposition(K, T)
        True

        sage: G = graphs.RandomGNP(10, .2)
        sage: T = G.treewidth(certificate=True)
        sage: is_valid_tree_decomposition(G, T)
        True

    The union of the bags is the set of vertices of `G`::

        sage: G = graphs.PathGraph(4)
        sage: T = G.treewidth(certificate=True)
        sage: _ = G.add_vertex()
        sage: is_valid_tree_decomposition(G, T)
        False

    Each edge of `G` is contained in a bag::

        sage: G = graphs.PathGraph(4)
        sage: T = G.treewidth(certificate=True)
        sage: G.add_edge(0, 3)
        sage: is_valid_tree_decomposition(G, T)
        False

    The bags containing a vertex `v` form a subtree of `T`::

        sage: G = graphs.PathGraph(4)
        sage: X1, X2, X3 = Set([0, 1]), Set([1, 2]), Set([2, 3])
        sage: T = Graph([(X1, X3), (X3, X2)])
        sage: is_valid_tree_decomposition(G, T)
        False

    TESTS:

    Check that both parameters are sage graphs::

        sage: is_valid_tree_decomposition("foo", Graph())
        Traceback (most recent call last):
        ...
        ValueError: the first parameter must be a sage Graph
        sage: is_valid_tree_decomposition(Graph(), "foo")
        Traceback (most recent call last):
        ...
        ValueError: the second parameter must be a sage Graph

    Check that `T` is a tree::

        sage: is_valid_tree_decomposition(Graph(), Graph(2))
        Traceback (most recent call last):
        ...
        ValueError: the second parameter must be a tree

    The vertices of `T` must be iterables::

        sage: is_valid_tree_decomposition(Graph(1), Graph(1))
        Traceback (most recent call last):
        ...
        ValueError: the vertices of T must be iterables

    Small cases::

        sage: is_valid_tree_decomposition(Graph(), Graph())
        True
        sage: is_valid_tree_decomposition(Graph(1), Graph())
        False
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("the first parameter must be a sage Graph")
    if not isinstance(T, Graph):
        raise ValueError("the second parameter must be a sage Graph")
    if not T:
        return not G
    elif not T.is_tree():
        raise ValueError("the second parameter must be a tree")

    for X in T:
        if not isinstance(X, Iterable):
            raise ValueError("the vertices of T must be iterables")

    # 1. The union of the bags equals V
    if set(G) != set(chain(*T)):
        return False

    # 2. Each edge of G is contained in a bag
    vertex_to_bags = {u: set() for u in G}
    for Xi in T:
        for u in Xi:
            vertex_to_bags[u].add(Xi)

    for u, v in G.edge_iterator(labels=False):
        if all(v not in Xi for Xi in vertex_to_bags[u]):
            return False

    # 3. The bags containing a vertex v form a connected subset of T
    for X in vertex_to_bags.values():
        D = DisjointSet(X)
        for Xi in X:
            for Xj in T.neighbor_iterator(Xi):
                if Xj in X:
                    D.union(Xi, Xj)
        if D.number_of_subsets() > 1:
            return False

    return True


def reduced_tree_decomposition(T):
    r"""
    Return a reduced tree-decomposition of `T`.

    We merge all edges between two sets `S` and `S'` where `S` is a subset of
    `S'`. To do so, we use a simple union-find data structure to record merge
    operations and the good sets.

    .. WARNING::

        This method assumes that the vertices of the input tree `T` are hashable
        and have attribute ``issuperset``, e.g., ``frozenset`` or
        :class:`~sage.sets.set.Set_object_enumerated`.

    INPUT:

    - ``T`` -- a tree-decomposition

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import reduced_tree_decomposition
        sage: from sage.graphs.graph_decompositions.tree_decomposition import is_valid_tree_decomposition
        sage: G = graphs.PathGraph(3)
        sage: T = Graph()
        sage: T.add_path([Set([0]), Set([0, 1]), Set([1]), Set([1, 2]),  Set([2])])
        sage: T.order()
        5
        sage: is_valid_tree_decomposition(G, T)
        True
        sage: T2 = reduced_tree_decomposition(T)
        sage: is_valid_tree_decomposition(G, T2)
        True
        sage: T2.order()
        2

    TESTS::

        sage: G = graphs.PathGraph(3)
        sage: T = G.treewidth(certificate=True)
        sage: is_valid_tree_decomposition(G, T)
        True
        sage: T == reduced_tree_decomposition(T)
        True
        sage: G = Graph(1)
        sage: T = G.treewidth(certificate=True)
        sage: T.order()
        1
        sage: T == reduced_tree_decomposition(T)
        True
    """
    if T.order() < 2:
        return T

    def get_ancestor(ancestor, u):
        if ancestor[u] == u:
            return u
        ancestor[u] = get_ancestor(ancestor, ancestor[u])
        return ancestor[u]

    ancestor = {u: u for u in T}
    for u, v in T.edge_iterator(labels=False):
        u = get_ancestor(ancestor, u)
        v = get_ancestor(ancestor, v)
        if u == v:
            continue
        elif u.issuperset(v):
            ancestor[v] = u
        elif v.issuperset(u):
            ancestor[u] = v

    from sage.graphs.graph import Graph
    H = Graph(multiedges=False, name="Reduced tree-decomposition of {}".format(T.name()))
    for u, v in T.edge_iterator(labels=False):
        u = get_ancestor(ancestor, u)
        v = get_ancestor(ancestor, v)
        if u != v:
            H.add_edge(u, v)
    return H


def width_of_tree_decomposition(G, T, check=True):
    r"""
    Return the width of the tree decomposition `T` of `G`.

    The width of a tree-decomposition is the size of the largest bag minus
    1. The empty graph and a graph of order 1 have treewidth 0.

    INPUT:

    - ``G`` -- a sage Graph

    - ``T`` -- a tree-decomposition for `G`

    - ``check`` -- boolean (default: ``True``); whether to check that the
      tree-decomposition `T` is valid for `G`

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import width_of_tree_decomposition
        sage: G = graphs.PathGraph(3)
        sage: T = G.treewidth(certificate=True)
        sage: width_of_tree_decomposition(G, T, check=True)
        1

    TESTS::

        sage: G = Graph()
        sage: T = G.treewidth(certificate=True)
        sage: width_of_tree_decomposition(G, T, check=True)
        0
        sage: width_of_tree_decomposition(Graph(1), T, check=True)
        Traceback (most recent call last):
        ...
        ValueError: the tree-decomposition is not valid for this graph
    """
    if check and not is_valid_tree_decomposition(G, T):
        raise ValueError("the tree-decomposition is not valid for this graph")

    if T:
        return max(0, max(len(u) for u in T) - 1)
    return 0


def _from_tree_decompositions_of_atoms_to_tree_decomposition(T_atoms, cliques):
    r"""
    Return a tree-decomposition formed by the tree-decompositions of the atoms.

    This is a helper method to avoid duplicated code.

    This method builds the tree decomposition of the graph by connecting the
    tree decompositions of its atoms. This is done in an order that is
    consistent with the order of the atoms and cliques returned by method
    :meth:`~Graph.atoms_and_clique_separators`. More precisely, the first clique
    separates the first atom from the rest of the graph (call G1 this part of
    the graph), the second clique separates (in G1) the second atom from the
    rest of the graph G1, etc. So we merge the tree decompositions in the
    reverse order of the atoms.

    INPUT:

    - ``T_atoms`` -- list of tree-decompositions

    - ``cliques`` -- list of clique separators

    EXAMPLES:

    Indirect doctest::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import is_valid_tree_decomposition
        sage: G = graphs.Grid2dGraph(2, 3)
        sage: T = G.treewidth(algorithm='sage', certificate=True)
        sage: is_valid_tree_decomposition(G, T)
        True

    TESTS::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import _from_tree_decompositions_of_atoms_to_tree_decomposition
        sage: _from_tree_decompositions_of_atoms_to_tree_decomposition([], [Set(range(2))])
        Traceback (most recent call last):
        ...
        ValueError: the number of cliques must be one less than the number of tree-decompositions of atoms
    """
    if len(T_atoms) != len(cliques) + 1:
        raise ValueError("the number of cliques must be one less than the "
                         "number of tree-decompositions of atoms")

    T = T_atoms[-1]
    for i in range(len(cliques) - 1, -1, -1):
        A = T_atoms[i]
        C = cliques[i]

        # We search for a vertex in A and T containing clique C
        ua, ut = None, None
        for u in A:
            if u.issuperset(C):
                ua = u
                break
        for u in T:
            if u.issuperset(C):
                ut = u
                break
        if ua and ut:
            # We merge T and A
            T.add_vertices(A)
            T.add_edges(A.edges(sort=False))
            T.add_edge(ua, ut)
        else:
            # This should never happen
            raise RuntimeError("something goes wrong. Please report the issue "
                               "to sage-devel@googlegroups.com")

    return T


def treewidth(g, k=None, kmin=None, certificate=False, algorithm=None, nice=False):
    r"""
    Compute the treewidth of `g` (and provide a decomposition).

    INPUT:

    - ``g`` -- a Sage Graph

    - ``k`` -- integer (default: ``None``); indicates the width to be
      considered. When ``k`` is an integer, the method checks that the graph has
      treewidth `\leq k`. If ``k`` is ``None`` (default), the method computes
      the optimal treewidth.

    - ``kmin`` -- integer (default: ``None``); when specified, search for a
      tree-decomposition of width at least ``kmin``. This parameter is useful
      when the graph can be decomposed into atoms.  This parameter is ignored
      when ``k`` is not ``None`` or when ``algorithm == 'tdlib'``.

    - ``certificate`` -- boolean (default: ``False``); whether to return the
      tree-decomposition itself

    - ``algorithm`` -- whether to use ``'sage'`` or ``'tdlib'`` (requires the
      installation of the :ref:`spkg_sagemath_tdlib` package). The default behaviour is to use
      'tdlib' if it is available, and Sage's own algorithm when it is not.

    - ``nice`` -- boolean (default: ``False``); whether or not to return the
      nice tree decomposition, provided ``certificate`` is ``True``

    OUTPUT:

    ``g.treewidth()`` returns treewidth of the graph ``g``.

    When ``k`` is specified, it returns ``False`` if there is no tree
    decomposition of width `\leq k`, and ``True`` otherwise.

    When ``certificate=True``, the tree decomposition is returned.

    When ``nice=True``, the nice tree decomposition is returned.

    ALGORITHM:

    This function virtually explores the graph of all pairs ``(vertex_cut, cc)``,
    where ``vertex_cut`` is a vertex cut of the graph of cardinality `\leq k+1`,
    and ``connected_component`` is a connected component of the graph induced by
    ``G-vertex_cut``.

    We deduce that the pair ``(vertex_cut, cc)`` is feasible with treewidth `k`
    if ``cc`` is empty, or if a vertex ``v`` from ``vertex_cut`` can be replaced
    with a vertex from ``cc``, such that the pair ``(vertex_cut+v,cc-v)`` is
    feasible.

    .. NOTE::

        The implementation would be much faster if ``cc``, the argument of the
        recursive function, was a bitset. It would also be very nice to not copy
        the graph in order to compute connected components, for this is really a
        waste of time.

    .. SEEALSO::

        :meth:`~sage.graphs.graph_decompositions.vertex_separation.path_decomposition`
        computes the pathwidth of a graph. See also the
        :mod:`~sage.graphs.graph_decompositions.vertex_separation` module.

    EXAMPLES:

    The PetersenGraph has treewidth 4::

        sage: petersen = graphs.PetersenGraph()
        sage: petersen.treewidth(algorithm='sage')
        4
        sage: petersen.treewidth(algorithm='sage', certificate=True)
        Tree decomposition: Graph on 6 vertices

    The PetersenGraph has treewidth 4 (with ``tdlib``)::

        sage: petersen = graphs.PetersenGraph()
        sage: petersen.treewidth(algorithm='tdlib')                    # optional - tdlib
        4
        sage: petersen.treewidth(algorithm='tdlib', certificate=True)  # optional - tdlib
        Tree decomposition: Graph on 6 vertices

    Nice tree decomposition of the PetersenGraph has 28 nodes::

        sage: petersen = graphs.PetersenGraph()
        sage: petersen.treewidth(algorithm='sage', certificate=True, nice=True)
        Nice tree decomposition of Tree decomposition: Graph on 28 vertices
        sage: petersen.treewidth(algorithm='tdlib', certificate=True, nice=True)  # optional - tdlib
        Nice tree decomposition of Tree decomposition: Graph on 28 vertices

    The treewidth of a 2-dimensional grid is its smallest side::

        sage: graphs.Grid2dGraph(2,5).treewidth()
        2
        sage: graphs.Grid2dGraph(3,5).treewidth()
        3

    When parameter ``kmin`` is specified, the method searches for a
    tree-decomposition of width at least ``kmin``::

        sage: g = graphs.PetersenGraph()
        sage: g.treewidth()
        4
        sage: g.treewidth(kmin=2, algorithm='sage')
        4
        sage: g.treewidth(kmin=g.order(), certificate=True, algorithm='sage')
        Tree decomposition: Graph on 1 vertex

    TESTS::

        sage: g = graphs.PathGraph(3)
        sage: g.treewidth()
        1
        sage: g = 2*graphs.PathGraph(3)
        sage: g.treewidth()
        1
        sage: g.treewidth(certificate=True)
        Tree decomposition: Graph on 4 vertices
        sage: g.treewidth(2)
        True
        sage: g.treewidth(1)
        True
        sage: Graph(1).treewidth()
        0
        sage: Graph(0).treewidth()
        -1
        sage: g = graphs.PetersenGraph()
        sage: g.treewidth(k=2)
        False
        sage: g.treewidth(k=6)
        True
        sage: g.treewidth(k=6, kmin=2)
        True
        sage: g.treewidth(certificate=True).is_tree()
        True
        sage: g.treewidth(k=3, certificate=True)
        False
        sage: T = g.treewidth(k=4, certificate=True)
        sage: T
        Tree decomposition: Graph on 6 vertices
        sage: from sage.graphs.graph_decompositions.tree_decomposition import is_valid_tree_decomposition
        sage: is_valid_tree_decomposition(g, T)
        True

    All edges do appear (:issue:`17893`)::

        sage: from itertools import combinations
        sage: g = graphs.PathGraph(10)
        sage: td = g.treewidth(certificate=True)
        sage: for bag in td:
        ....:    g.delete_edges(list(combinations(bag,2)))
        sage: g.size()
        0

    :issue:`19358`::

        sage: g = Graph()
        sage: for i in range(3):
        ....:     for j in range(2):
        ....:         g.add_path([i,(i,j),(i+1)%3])
        sage: g.treewidth()
        2

    The decomposition is a tree (:issue:`23546`)::

        sage: g = Graph({0:[1,2], 3:[4,5]})
        sage: t = g.treewidth(certificate=True)
        sage: t.is_tree()
        True
        sage: vertices = set()
        sage: for s in t:
        ....:     vertices = vertices.union(s)
        sage: vertices == set(g)
        True

    Check that the use of atoms and clique separators is correct
    (:issue:`30993`)::

        sage: g = 2 * graphs.Grid2dGraph(2, 3)
        sage: g.treewidth(algorithm='sage')
        2
        sage: g.treewidth(algorithm='sage', certificate=True)
        Tree decomposition: Graph on 8 vertices
        sage: g.treewidth(algorithm='sage', certificate=True, kmin=4)
        Tree decomposition: Graph on 4 vertices

    Check that :issue:`38159` is fixed ::

        sage: G = Graph('I~~}vPlr_')
        sage: G.treewidth(algorithm='sage') == G.treewidth(algorithm='tdlib')  # optional - tdlib
        True

    Trivially true::

        sage: graphs.PetersenGraph().treewidth(k=35)
        True
        sage: graphs.PetersenGraph().treewidth(k=35, certificate=True)
        Tree decomposition: Graph on 1 vertex

    Bad input::

        sage: graphs.PetersenGraph().treewidth(k=-3)
        Traceback (most recent call last):
        ...
        ValueError: k(=-3) must be a nonnegative integer
    """
    # Check Input
    if algorithm not in [None, 'tdlib', 'sage']:
        raise ValueError("'algorithm' must be equal to 'tdlib', 'sage', or None")

    try:
        import sage.graphs.graph_decompositions.tdlib as tdlib
        tdlib_found = True
    except ImportError:
        tdlib_found = False

    if algorithm is None:
        algorithm = 'tdlib' if tdlib_found else 'sage'

    if (k is not None) and k < 0:
        raise ValueError(f"k(={k}) must be a nonnegative integer")

    # Silly cases
    from sage.graphs.graph import Graph
    if not g.order():
        if certificate:
            return Graph()
        return -1 if k is None else True

    if k is not None and k >= g.order() - 1:
        if certificate:
            return Graph({Set(g): []}, name="Tree decomposition")
        return True

    kmin = 0 if kmin is None else kmin
    if k is None and kmin >= g.order() - 1:
        if certificate:
            return Graph({Set(g): []}, name="Tree decomposition")
        return kmin

    # TDLIB
    if algorithm == 'tdlib':
        if not tdlib_found:
            from sage.features import FeatureNotPresentError
            from sage.features.tdlib import Tdlib
            raise FeatureNotPresentError(Tdlib())

        tree_decomp = tdlib.treedecomposition_exact(g, -1 if k is None else k)
        width = tdlib.get_width(tree_decomp)

        if certificate:
            if k is None or width <= k:
                return make_nice_tree_decomposition(g, tree_decomp) if nice else tree_decomp
            return False
        return width if k is None else width <= k

    # The treewidth of a graph is the maximum over its atoms. So, we decompose
    # the graph by clique minimal separators, compute the treewidth of each of
    # its atoms, and combine the results.
    # This decomposition also deals with disconnected cases.
    atoms, cliques = g.atoms_and_clique_separators()
    if cliques:
        # If we do not need the tree decomposition
        if not certificate:
            if k is None:
                for a in atoms:
                    kmin = max(kmin, g.subgraph(a).treewidth(algorithm=algorithm, kmin=kmin))
                return kmin
            elif max(len(c) for c in cliques) - 1 > k:
                return False
            else:
                return all(g.subgraph(a).treewidth(algorithm=algorithm, k=k) for a in atoms)

        # Otherwise, compute the tree decomposition of each atom
        T = []
        for a in atoms:
            ga = g.subgraph(a)
            Ta = ga.treewidth(algorithm=algorithm, certificate=True, kmin=kmin)
            kmin = max(kmin, width_of_tree_decomposition(ga, Ta, check=False))
            T.append(Ta)

        # Merge the resulting trees
        tree_decomp = _from_tree_decompositions_of_atoms_to_tree_decomposition(T, cliques)

        return make_nice_tree_decomposition(g, tree_decomp) if nice else tree_decomp

    # Forcing k to be defined
    if k is None:
        for i in range(max(kmin, g.clique_number() - 1, min(g.degree())), g.order()):
            ans = g.treewidth(algorithm=algorithm, k=i, certificate=certificate, nice=nice)
            if ans:
                return ans if certificate else i

    # This is the recursion described in the method's documentation. All
    # computations are cached, and depends on the pair ``cut,
    # connected_component`` only.
    #
    # It returns either a boolean or the corresponding tree-decomposition, as a
    # list of edges between vertex cuts (as it is done for the complete
    # tree-decomposition at the end of the main function.
    @cached_function
    def rec(cut, cc):
        # Easy cases
        if len(cut) > k:
            return False
        if len(cc) + len(cut) <= k + 1:
            return [(cut, cut.union(cc))] if certificate else True

        # We explore all possible extensions of the cut
        for v in cc:

            # New cuts and connected components, with v respectively added
            # and removed
            cutv = cut.union([v])
            ccv = cc.difference([v])

            # The values returned by the recursive calls.
            sons = []

            # Removing v may have disconnected cc. We iterate on its
            # connected components
            for cci in g.subgraph(ccv).connected_components(sort=False):

                # The recursive subcalls. We remove on-the-fly the vertices
                # from the cut which play no role in separating the
                # connected component from the rest of the graph.
                reduced_cut = frozenset([x for x in cutv
                                         if any(xx in cci for xx in g.neighbor_iterator(x))])
                son = rec(reduced_cut, frozenset(cci))
                if not son:
                    break

                if certificate:
                    sons.extend(son)
                    sons.append((cut, cutv))
                    sons.append((cutv, reduced_cut))

            # Weird Python syntax which is useful once in a lifetime : if break
            # was never called in the loop above, we return "sons".
            else:
                return sons if certificate else True

        return False

    # Main call to rec function, i.e. rec({v}, V-{v})
    V = list(g)
    v = frozenset([V.pop()])
    TD = rec(v, frozenset(V))

    if TD is False:
        return False

    if not certificate:
        return True

    # Building the tree-decomposition graph. Its vertices are cuts of the
    # decomposition, and there is an edge from a cut C1 to a cut C2 if C2 is an
    # immediate subcall of C1
    tree_decomp = Graph()
    tree_decomp.add_edges(((Set(x), Set(y)) for x, y in TD), loops=False)

    # The tree-decomposition may contain a lot of useless nodes.
    # We merge all edges between two sets S, S' where S is a subset of S'
    tree_decomp = reduced_tree_decomposition(tree_decomp)

    tree_decomp.name("Tree decomposition")
    if nice:
        tree_decomp = make_nice_tree_decomposition(g, tree_decomp)

    return tree_decomp


def make_nice_tree_decomposition(graph, tree_decomp):
    r"""
    Return a *nice* tree decomposition (TD) of the TD ``tree_decomp``.

    See page 161 of [CFKLMPPS15]_ for a description of the nice tree decomposition.

    A *nice* TD `NT` is a rooted tree with four types of nodes:

    - *Leaf* nodes have no children and bag size 1;
    - *Introduce* nodes have one child: If `v \in NT` is an introduce node and
      `w \in NT` its child, then `Bag(v) = Bag(w) \cup \{ x \}`, where `x` is the
      introduced node;
    - *Forget* nodes have one child: If `v \in NT` is a forget node and
      `w \in NT` its child, then `Bag(v) = Bag(w) \setminus \{ x \}`, where `x` is the
      forgotten node;
    - *Join* nodes have two children, both identical to the parent.

    INPUT:

    - ``graph`` -- a Sage graph

    - ``tree_decomp`` -- a tree decomposition

    OUTPUT: a nice tree decomposition

    .. WARNING::

        This method assumes that the vertices of the input tree `tree_decomp`
        are hashable and have attribute ``issuperset``, e.g., ``frozenset`` or
        :class:`~sage.sets.set.Set_object_enumerated_with_category`.

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import make_nice_tree_decomposition
        sage: petersen = graphs.PetersenGraph()
        sage: petersen_TD = petersen.treewidth(certificate=True)
        sage: make_nice_tree_decomposition(petersen, petersen_TD)
        Nice tree decomposition of Tree decomposition: Graph on 28 vertices

    ::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import make_nice_tree_decomposition
        sage: cherry = graphs.CompleteBipartiteGraph(1, 2)
        sage: cherry_TD = cherry.treewidth(certificate=True)
        sage: make_nice_tree_decomposition(cherry, cherry_TD)
        Nice tree decomposition of Tree decomposition: Graph on 7 vertices

    ::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import make_nice_tree_decomposition
        sage: bip_one_four = graphs.CompleteBipartiteGraph(1, 4)
        sage: bip_one_four_TD = bip_one_four.treewidth(certificate=True)
        sage: make_nice_tree_decomposition(bip_one_four, bip_one_four_TD)
        Nice tree decomposition of Tree decomposition: Graph on 15 vertices

    Check that :issue:`36843` is fixed::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import make_nice_tree_decomposition
        sage: triangle = graphs.CompleteGraph(3)
        sage: triangle_TD = triangle.treewidth(certificate=True)
        sage: make_nice_tree_decomposition(triangle, triangle_TD)
        Nice tree decomposition of Tree decomposition: Graph on 7 vertices

    ::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import make_nice_tree_decomposition
        sage: graph = graphs.CompleteBipartiteGraph(2, 5)
        sage: graph_TD = graph.treewidth(certificate=True)
        sage: make_nice_tree_decomposition(graph, graph_TD)
        Nice tree decomposition of Tree decomposition: Graph on 25 vertices

    ::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import make_nice_tree_decomposition
        sage: empty_graph = graphs.EmptyGraph()
        sage: tree_decomp = empty_graph.treewidth(certificate=True)
        sage: len(make_nice_tree_decomposition(empty_graph, tree_decomp))
        0

    ::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import make_nice_tree_decomposition
        sage: singleton = graphs.CompleteGraph(1)
        sage: tree_decomp = singleton.treewidth(certificate=True)
        sage: make_nice_tree_decomposition(singleton, tree_decomp)
        Nice tree decomposition of Tree decomposition: Graph on 3 vertices

    ::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import make_nice_tree_decomposition
        sage: an_edge = graphs.CompleteGraph(2)
        sage: tree_decomp = an_edge.treewidth(certificate=True)
        sage: make_nice_tree_decomposition(an_edge, tree_decomp)
        Nice tree decomposition of Tree decomposition: Graph on 5 vertices
    """
    if not is_valid_tree_decomposition(graph, tree_decomp):
        raise ValueError("input must be a valid tree decomposition for this graph")

    name = f"Nice tree decomposition of {tree_decomp.name()}"
    from sage.graphs.graph import Graph
    if not tree_decomp:
        return Graph(name=name)

    # Step 1: Ensure the tree is directed and has a root
    # Choose a root and orient the edges from root-to-leaves direction
    #
    # Testing <= 1 for the special case when one bag containing all vertices
    leaves = [u for u in tree_decomp if tree_decomp.degree(u) <= 1]

    from sage.graphs.digraph import DiGraph
    if len(leaves) == 1:
        root = leaves[0]
        directed_tree = DiGraph(tree_decomp)
    else:
        root = leaves.pop()

        directed_tree = DiGraph(tree_decomp.breadth_first_search(start=root, edges=True),
                                format='list_of_edges')

    # Relabel the graph in range (0, |tree_decomp| - 1)
    bags_to_int = directed_tree.relabel(inplace=True, return_map=True)
    # Get the new name of the root node
    root = bags_to_int[root]
    # Force bags to be of type Set to simplify code
    bag = {ui: Set(u) for u, ui in bags_to_int.items()}

    # Step 2: Add the root node and the leaf nodes, with empty bags
    # To each leaf node of `directed_tree`, we add a child with empty bag.
    # We also add a new root with empty bag.
    root, old_root = directed_tree.add_vertex(), root
    directed_tree.add_edge(root, old_root)
    bag[root] = Set()
    for vi, u in enumerate(leaves, start=root + 1):
        directed_tree.add_edge(bags_to_int[u], vi)
        bag[vi] = Set()

    # Step 3: Ensure that each node of directed_tree has at most 2 children.
    # If a node has more than 2 children, introduce new nodes to
    # make sure each node has at most 2 children:
    #
    # If v has k > 2 children (w_1, w_2, ..., w_k), we disconnect (w_1, ..., w_{k-1})
    # from v, and introduce k - 2 new nodes (u_1, u_2, ..., u_{k-2}).
    # We then let w_i be the children of u_i for 1 <= i <= k - 2.
    # We also let w_{k-1} be the second child of u_{k-2}, and
    # u_i the second child of u_{i-1}.
    # Finally, we let u_1 the second child of u.
    # Each node u_i has the same bag as u.

    # We need to call list(...) since we modify directed_tree
    for ui in list(directed_tree):
        if directed_tree.out_degree(ui) > 2:
            children = directed_tree.neighbors_out(ui)
            children.pop() # one vertex remains a child of ui

            directed_tree.delete_edges((ui, vi) for vi in children)

            new_nodes = [directed_tree.add_vertex() for _ in range(len(children) - 1)]

            directed_tree.add_edge(ui, new_nodes[0])
            directed_tree.add_path(new_nodes)
            directed_tree.add_edges(zip(new_nodes, children))
            directed_tree.add_edge(new_nodes[-1], children[-1])

            bag.update((vi, bag[ui]) for vi in new_nodes)

    # Step 4: If current vertex v has two children w1 and w2,
    # then bag[v] == bag[w1] == bag[w2]
    for current_node in list(directed_tree):
        if directed_tree.out_degree(current_node) < 2:
            continue
        for neighbor in directed_tree.neighbor_out_iterator(current_node):
            if bag[current_node] != bag[neighbor]:
                directed_tree.delete_edge(current_node, neighbor)
                new_node = directed_tree.add_vertex()
                directed_tree.add_path([current_node, new_node, neighbor])
                bag[new_node] = bag[current_node]

    # Step 5: If the node v has only one child, then it is either an introduce
    # node or a forget node.
    def add_path_of_intro_nodes(u, v):
        """
        Replace the arc (u, v) by a path of introduce nodes.
        """
        if len(bag[u]) + 1 == len(bag[v]):
            return

        diff = list(bag[v] - bag[u])
        diff.pop()

        last_node = u
        for w in diff:
            new_node = directed_tree.add_vertex()
            bag[new_node] = bag[last_node].union(Set((w,)))
            directed_tree.add_edge(last_node, new_node)
            last_node = new_node

        directed_tree.add_edge(last_node, v)
        directed_tree.delete_edge(u, v)

    def add_path_of_forget_nodes(u, v):
        """
        Replace the arc (u, v) by a path of forget nodes.
        """
        if len(bag[v]) + 1 == len(bag[u]):
            return

        diff = list(bag[u] - bag[v])
        diff.pop()

        last_node = u
        for w in diff:
            new_node = directed_tree.add_vertex()
            bag[new_node] = bag[last_node] - {w}
            directed_tree.add_edge(last_node, new_node)
            last_node = new_node

        directed_tree.add_edge(last_node, v)
        directed_tree.delete_edge(u, v)

    for ui in list(directed_tree):
        if directed_tree.out_degree(ui) != 1:
            continue

        vi = next(directed_tree.neighbor_out_iterator(ui))
        bag_ui, bag_vi = bag[ui], bag[vi]

        # Merge the nodes if the two bags are the same
        if bag_ui == bag_vi:
            if directed_tree.in_degree(ui) == 1:
                parent = next(directed_tree.neighbor_in_iterator(ui))
                directed_tree.add_edge(parent, vi)
            else:
                root = vi
            directed_tree.delete_vertex(ui)

        # Add paths of intro / forget nodes accordingly

        elif bag_ui.issubset(bag_vi):
            add_path_of_intro_nodes(ui, vi)

        elif bag_vi.issubset(bag_ui):
            add_path_of_forget_nodes(ui, vi)

        # Handle the case when the two nodes are not related in any way above
        else:
            wi = directed_tree.add_vertex()
            bag[wi] = bag[ui] & bag[vi]
            directed_tree.add_path([ui, wi, vi])
            directed_tree.delete_edge(ui, vi)
            add_path_of_forget_nodes(ui, wi)
            add_path_of_intro_nodes(wi, vi)

    # Return the nice tree decomposition after the processing
    nice_tree_decomp = Graph(directed_tree, name=name)

    bfs_ordering = nice_tree_decomp.breadth_first_search(start=root)
    relabeling = {u: (i, bag[u]) for i, u in enumerate(bfs_ordering)}
    nice_tree_decomp.relabel(inplace=True, perm=relabeling)

    return nice_tree_decomp


def label_nice_tree_decomposition(nice_TD, root, directed=False):
    r"""
    Return a nice tree decomposition with nodes labelled accordingly.

    INPUT:

    - ``nice_TD`` -- a nice tree decomposition

    - ``root`` -- the root of the nice tree decomposition

    - ``directed`` -- boolean (default: ``False``); whether to return the nice
      tree decomposition as a directed graph rooted at vertex ``root`` or as an
      undirected graph

    OUTPUT: a nice tree decomposition with nodes labelled

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import make_nice_tree_decomposition, label_nice_tree_decomposition
        sage: claw = graphs.CompleteBipartiteGraph(1, 3)
        sage: claw_TD = claw.treewidth(certificate=True)
        sage: nice_TD = make_nice_tree_decomposition(claw, claw_TD)
        sage: root = sorted(nice_TD.vertices())[0]
        sage: label_TD = label_nice_tree_decomposition(nice_TD, root, directed=True)
        sage: label_TD.name()
        'Labelled Nice tree decomposition of Tree decomposition'
        sage: for node in sorted(label_TD):  # random
        ....:     print(node, label_TD.get_vertex(node))
        (0, {}) forget
        (1, {0}) forget
        (2, {0, 1}) intro
        (3, {0}) forget
        (4, {0, 3}) intro
        (5, {0}) forget
        (6, {0, 2}) intro
        (7, {2}) intro
        (8, {}) leaf
    """
    from sage.graphs.digraph import DiGraph
    from sage.graphs.graph import Graph

    directed_TD = DiGraph(nice_TD.breadth_first_search(start=root, edges=True),
                          format='list_of_edges',
                          name='Labelled {}'.format(nice_TD.name()))

    # The loop starts from the root node
    # We assume the tree decomposition is valid and nice,
    # hence saving time on checking.
    for node in directed_TD:
        out_deg = directed_TD.out_degree(node)

        if out_deg == 2:
            directed_TD.set_vertex(node, 'join')
        elif out_deg == 1:
            current_bag = node[1]
            child_bag = directed_TD.neighbors_out(node)[0][1]

            if len(current_bag) == len(child_bag) + 1:
                directed_TD.set_vertex(node, 'intro')
            else:
                directed_TD.set_vertex(node, 'forget')
        else:
            directed_TD.set_vertex(node, 'leaf')

    if directed:
        return directed_TD
    return Graph(directed_TD, name=nice_TD.name())


#
# Treelength
#

def length_of_tree_decomposition(G, T, check=True):
    r"""
    Return the length of the tree decomposition `T` of `G`.

    The *length* of a tree decomposition, as proposed in [DG2006]_, is the
    maximum *diameter* in `G` of its bags, where the diameter of a bag `X_i` is
    the largest distance in `G` between the vertices in `X_i` (i.e., `\max_{u, v
    \in X_i} \dist_G(u, v)`). See the documentation of the
    :mod:`~sage.graphs.graph_decompositions.tree_decomposition` module for more
    details.

    INPUT:

    - ``G`` -- a graph

    - ``T`` -- a tree-decomposition for `G`

    - ``check`` -- boolean (default: ``True``); whether to check that the
      tree-decomposition `T` is valid for `G`

    EXAMPLES:

    Trees and cliques have treelength 1::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import length_of_tree_decomposition
        sage: G = graphs.CompleteGraph(5)
        sage: tl, T = G.treelength(certificate=True)
        sage: tl
        1
        sage: length_of_tree_decomposition(G, T, check=True)
        1
        sage: G = graphs.RandomTree(20)
        sage: tl, T = G.treelength(certificate=True)
        sage: tl
        1
        sage: length_of_tree_decomposition(G, T, check=True)
        1

    The Petersen graph has treelength 2::

        sage: G = graphs.PetersenGraph()
        sage: tl, T = G.treelength(certificate=True)
        sage: tl
        2
        sage: length_of_tree_decomposition(G, T)
        2

    When a tree-decomposition has a single bag containing all vertices of a
    graph, the length of this tree-decomposition is the diameter of the graph::

        sage: G = graphs.Grid2dGraph(2, 5)
        sage: G.treelength()
        2
        sage: G.diameter()
        5
        sage: T = Graph({Set(G): []})
        sage: length_of_tree_decomposition(G, T)
        5

    TESTS::

        sage: G = Graph()
        sage: _, T = G.treelength(certificate=True)
        sage: length_of_tree_decomposition(G, T, check=True)
        0
        sage: length_of_tree_decomposition(Graph(1), T, check=True)
        Traceback (most recent call last):
        ...
        ValueError: the tree-decomposition is not valid for this graph
    """
    if check and not is_valid_tree_decomposition(G, T):
        raise ValueError("the tree-decomposition is not valid for this graph")

    cdef unsigned int n = G.order()

    if n < 2:
        return 0
    if any(len(bag) == n for bag in T):
        return G.diameter()

    cdef unsigned int i, j

    # We map vertices to integers in range 0..n-1
    cdef list int_to_vertex = list(G)
    cdef dict vertex_to_int = {u: i for i, u in enumerate(int_to_vertex)}

    # We compute the distance matrix.
    cdef unsigned short * c_distances = c_distances_all_pairs(G, vertex_list=int_to_vertex)
    cdef unsigned short ** distances = <unsigned short **>sig_calloc(n, sizeof(unsigned short *))
    for i in range(n):
        distances[i] = c_distances + i * n

    # We now compute the maximum lengths of the bags
    from itertools import combinations
    cdef list bag_int
    cdef unsigned short dij
    cdef unsigned short length = 0
    for bag in T:
        bag_int = [vertex_to_int[u] for u in bag]
        for i, j in combinations(bag_int, 2):
            dij = distances[i][j]
            if dij > length:
                length = dij

    sig_free(c_distances)
    sig_free(distances)

    return length


def treelength_lowerbound(G):
    r"""
    Return a lower bound on the treelength of `G`.

    See [DG2006]_ for more details.

    INPUT:

    - ``G`` -- a sage Graph

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import treelength_lowerbound
        sage: G = graphs.PetersenGraph()
        sage: treelength_lowerbound(G)
        1
        sage: G.treelength()
        2
        sage: G = graphs.CycleGraph(5)
        sage: treelength_lowerbound(G)
        2
        sage: G.treelength()
        2

    TESTS::

        sage: treelength_lowerbound(Graph())
        0
    """
    if G.is_cycle():
        from sage.arith.misc import integer_ceil as ceil
        return int(ceil(G.order() / 3))

    lowerbound = 0
    girth = G.girth()
    if girth is not Infinity:
        lowerbound = max(lowerbound, girth / 3)

    return int(lowerbound)


cdef class TreelengthConnected:
    r"""
    Compute the treelength of a connected graph (and provide a decomposition).

    This class implements an algorithm for computing the treelength of a
    connected graph that virtually explores the graph of all pairs
    ``(vertex_cut, connected_component)``, where ``vertex_cut`` is a vertex cut
    of the graph of length `\leq k`, and ``connected_component`` is a connected
    component of the graph induced by ``G - vertex_cut``.

    We deduce that the pair ``(vertex_cut, connected_component)`` is feasible
    with treelength `k` if ``connected_component`` is empty, or if a vertex
    ``v`` from ``vertex_cut`` can be replaced with a vertex from
    ``connected_component``, such that the pair ``(vertex_cut + v,
    connected_component - v)`` is feasible.

    INPUT:

    - ``G`` -- a sage Graph

    - ``k`` -- integer (default: ``None``); indicates the length to be
      considered. When `k` is an integer, the method checks that the graph has
      treelength `\leq k`. If `k` is ``None`` (default), the method computes the
      optimal treelength.

    - ``certificate`` -- boolean (default: ``False``); whether to also compute
      the tree-decomposition itself

    OUTPUT:

    ``TreelengthConnected(G)`` returns the treelength of `G`. When `k` is
    specified, it returns ``False`` when no tree-decomposition of length
    `\leq k` exists or ``True`` otherwise. When ``certificate=True``, the
    tree-decomposition is also returned.

    EXAMPLES:

    A clique has treelength 1::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import TreelengthConnected
        sage: TreelengthConnected(graphs.CompleteGraph(3)).get_length()
        1
        sage: TC = TreelengthConnected(graphs.CompleteGraph(4), certificate=True)
        sage: TC.get_length()
        1
        sage: TC.get_tree_decomposition()
        Tree decomposition of Complete graph: Graph on 1 vertex

    A cycle has treelength `\lceil n/3 \rceil`::

        sage: TreelengthConnected(graphs.CycleGraph(6)).get_length()
        2
        sage: TreelengthConnected(graphs.CycleGraph(7)).get_length()
        3
        sage: TreelengthConnected(graphs.CycleGraph(7), k=3).is_less_than_k()
        True
        sage: TreelengthConnected(graphs.CycleGraph(7), k=2).is_less_than_k()
        False

    TESTS:

    The input graph must be connected::

        sage: TreelengthConnected(Graph(2))
        Traceback (most recent call last):
        ...
        ValueError: the graph is not connected

    The parameter `k` must be nonnegative::

        sage: TreelengthConnected(Graph(1), k=-1)
        Traceback (most recent call last):
        ...
        ValueError: k (= -1) must be a nonnegative integer

    Parameter ``certificate`` must be ``True`` to get a tree decomposition::

        sage: TreelengthConnected(Graph(1), certificate=False).get_tree_decomposition()
        Traceback (most recent call last):
        ...
        ValueError: parameter 'certificate' has not been set to True

    When parameter `k` is specified and ``certificate`` is ``True``, the
    computed tree decomposition is any valid tree decomposition with length at
    most `k`. However, this tree decomposition exists only if the treelength of
    `G` is at most `k` (i.e., `tl(G) \leq k`)::

        sage: G = graphs.Grid2dGraph(2, 3)
        sage: TC = TreelengthConnected(G, k=2, certificate=True)
        sage: TC.is_less_than_k()
        True
        sage: TC.get_tree_decomposition()
        Tree decomposition of 2D Grid Graph for [2, 3]: Graph on 3 vertices
        sage: TC = TreelengthConnected(G, k=1, certificate=True)
        sage: TC.is_less_than_k()
        False
        sage: TC.get_tree_decomposition()
        Traceback (most recent call last):
        ...
        ValueError: no tree decomposition with length <= 1 was found
    """

    def __init__(self, G, k=None, certificate=False):
        r"""
        Initialize this object and compute the treelength of `G`.

        INPUT:

        - ``G`` -- a sage Graph

        - ``k`` -- integer (default: ``None``); indicates the length to be
          considered. When `k` is an integer, the method checks that the graph
          has treelength `\leq k`. If `k` is ``None`` (default), the method
          computes the optimal treelength.

        - ``certificate`` -- boolean (default: ``False``); whether to compute
          the tree-decomposition itself

        TESTS::

            sage: from sage.graphs.graph_decompositions.tree_decomposition import TreelengthConnected
            sage: G = graphs.CycleGraph(4)
            sage: TreelengthConnected(G).get_length()
            2
        """
        if k is not None and k < 0:
            raise ValueError("k (= {}) must be a nonnegative integer".format(k))
        G._scream_if_not_simple()
        if not G.is_connected():
            raise ValueError("the graph is not connected")

        self.certificate = certificate
        self.k_is_defined = k is not None
        self.k = k if self.k_is_defined else 0

        if certificate:
            from sage.graphs.graph import Graph
            self.name = "Tree decomposition of {}".format(G.name())

        self.n = G.order()
        self.distances = NULL  # used in the destructor

        # Trivial cases
        if self.n <= 1 or (self.k_is_defined and self.n <= k):
            if certificate:
                if self.n:
                    self.tree = Graph({Set(G): []}, format='dict_of_lists', name=self.name)
                else:
                    self.tree = Graph(name=self.name)
            self.length = 0 if self.n <= 1 else G.diameter(algorithm='DHV')
            self.leq_k = True  # We know that k is nonnegative
            return

        if self.k_is_defined and not k:
            # We have at least 2 vertices and 1 edges, so tl >= 1
            self.leq_k = False
            return

        if G.is_clique():
            if certificate:
                self.tree = Graph({Set(G): []}, format='dict_of_lists', name=self.name)
            self.length = 1
            self.leq_k = True
            return

        cdef unsigned int i, j

        # If the vertices are not labeled 0..n-1, we relabel the graph. This
        # way, the labeling of the vertices matches the rows and columns of the
        # distance matrix
        if set(G) == set(range(self.n)):
            graph = G
            self.perm_inv = dict()
        else:
            graph, perm = G.relabel(inplace=False, return_map=True)
            self.perm_inv = {i: u for u, i in perm.items()}

        # We compute the distance matrix.
        self.c_distances = c_distances_all_pairs(graph, vertex_list=list(range(self.n)))
        self.distances = <unsigned short **>sig_calloc(self.n, sizeof(unsigned short *))
        for i in range(self.n):
            self.distances[i] = self.c_distances + i * self.n

        # and the diameter of the graph
        self.diameter = 0
        for i in range(self.n):
            for j in range(i, self.n):
                self.diameter = max(self.diameter, self.distances[i][j])

        if self.k_is_defined and k >= self.diameter:
            # All vertices fit in one bag
            if certificate:
                self.tree = Graph({Set(G): []}, format='dict_of_lists', name=self.name)
            self.length = self.diameter
            self.leq_k = True
            return

        # Forcing k to be defined
        if not self.k_is_defined:
            for i in range(treelength_lowerbound(graph), self.diameter + 1):
                ans = self._treelength(graph, i)
                if ans:
                    self.length = i
                    return

        # If k is defined
        ans = self._treelength(graph, k)
        if ans:
            self.length = k
            self.leq_k = True
        else:
            self.leq_k = False

    def __dealloc__(self):
        r"""
        Destroy the object.

        TESTS::

            sage: from sage.graphs.graph_decompositions.tree_decomposition import TreelengthConnected
            sage: G = graphs.CycleGraph(4)
            sage: TreelengthConnected(G).get_length()
            2
        """
        if self.distances:
            sig_free(self.c_distances)
            sig_free(self.distances)

    cdef bint _treelength(self, g, k) noexcept:
        r"""
        Check whether the treelength of `g` is at most `k`.

        INPUT:

        - ``g`` -- a sage Graph

        - ``k`` -- integer; indicates the length to be considered

        TESTS::

            sage: from sage.graphs.graph_decompositions.tree_decomposition import TreelengthConnected
            sage: G = graphs.CycleGraph(4)
            sage: TreelengthConnected(G, k=2).is_less_than_k()
            True
        """

        # This is the recursion described in the method's documentation. All
        # computations are cached, and depends on the pair ``cut,
        # connected_component`` only.
        #
        # It returns either a boolean or the corresponding tree-decomposition,
        # as a list of edges between vertex cuts (used to build the complete
        # tree-decomposition at the end of the _treelength method).
        @cached_function
        def rec(cut, cc):
            cdef int v
            cdef frozenset reduced_cut

            if len(cc) == 1:
                [v] = cc
                # We identify the neighbors of v in cut
                reduced_cut = cut.intersection(g.neighbor_iterator(v))
                # We can form a new bag with its closed neighborhood, and this
                # bag has diameter at most 2. Furthermore, if k == 1, we know
                # that the bag cut has diameter <= 1, and so the new bag has
                # diameter 1
                if self.certificate:
                    if cut == reduced_cut:
                        return [(cut, cut.union(cc))]
                    # We need to forget some vertices
                    return [(cut, reduced_cut), (reduced_cut, reduced_cut.union(cc))]

                return True

            # We explore all possible extensions of the cut
            cdef frozenset cutv
            cdef frozenset ccv
            cdef frozenset cci
            cdef frozenset reduced_cuti
            cdef list sons
            cdef int x

            for v in cc:

                # We know that the cut has diameter <= k. So we check is adding
                # v to the cut does not make its diameter > k
                if any(self.distances[v][x] > k for x in cut):
                    continue
                # We add v to the cut and remove it from cc
                cutv = cut.union([v])
                ccv = cc.difference([v])

                # The values returned by the recursive calls.
                sons = []

                # Removing v may have disconnected cc. We iterate on its
                # connected components
                for _cci in g.subgraph(ccv).connected_components(sort=False):
                    cci = frozenset(_cci)

                    # The recursive subcalls. We remove on-the-fly the vertices
                    # from the cut which play no role in separating the
                    # connected component from the rest of the graph. That is,
                    # we identify the vertices of cutv with a neighbor in cci
                    reduced_cuti = frozenset([x for x in cutv
                                              if any(xx in cci for xx in g.neighbor_iterator(x))])
                    if not reduced_cuti:
                        # This should not happen
                        break

                    # and we do a recursive call
                    son = rec(reduced_cuti, cci)
                    if not son:
                        break

                    # We get a valid decomposition of cci
                    if self.certificate:
                        # We connect cut, cutv, reduced_cci and son
                        sons.append((cut, cutv))
                        if v in reduced_cuti:
                            sons.append((cutv, reduced_cuti))
                        else:
                            reduced_cutv = reduced_cuti.union([v])
                            sons.append((cutv, reduced_cutv))
                            sons.append((reduced_cutv, reduced_cuti))
                        sons.extend(son)

                # Weird Python syntax which is useful once in a lifetime : if break
                # was never called in the loop above, we return "sons".
                else:
                    return sons if self.certificate else True

            return False

        # Main call to rec function, i.e. rec({v}, V-{v})
        cdef list V = list(g)
        cdef frozenset v = frozenset([V.pop()])
        TD = rec(v, frozenset(V))

        if TD is False:
            return False

        if not self.certificate:
            return True

        # Building the tree-decomposition graph. Its vertices are cuts of the
        # decomposition, and there is an edge from a cut C1 to a cut C2 if C2 is an
        # immediate subcall of C1. If needed, the vertices are relabeled.
        if self.perm_inv:
            def good_label(x):
                return Set([self.perm_inv[i] for i in x])
        else:
            def good_label(x):
                return Set(x)

        from sage.graphs.graph import Graph
        T = Graph([(good_label(x), good_label(y)) for x, y in TD if x != y],
                  format='list_of_edges')
        self.tree = reduced_tree_decomposition(T)
        self.tree.name(self.name)
        return True

    def get_tree_decomposition(self):
        """
        Return the tree-decomposition.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.tree_decomposition import TreelengthConnected
            sage: G = graphs.CycleGraph(4)
            sage: TreelengthConnected(G, certificate=True).get_tree_decomposition()
            Tree decomposition of Cycle graph: Graph on 2 vertices
            sage: G.diameter()
            2
            sage: TreelengthConnected(G, k=2, certificate=True).get_tree_decomposition()
            Tree decomposition of Cycle graph: Graph on 1 vertex
            sage: TreelengthConnected(G, k=1, certificate=True).get_tree_decomposition()
            Traceback (most recent call last):
            ...
            ValueError: no tree decomposition with length <= 1 was found

        TESTS::

            sage: G = graphs.CycleGraph(4)
            sage: TreelengthConnected(G, certificate=False).get_tree_decomposition()
            Traceback (most recent call last):
            ...
            ValueError: parameter 'certificate' has not been set to True
        """
        if self.certificate:
            if self.k_is_defined and not self.leq_k:
                raise ValueError("no tree decomposition with length <= {} was found".format(self.k))
            return self.tree
        else:
            raise ValueError("parameter 'certificate' has not been set to True")

    def get_length(self):
        """
        Return the length of the tree decomposition.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.tree_decomposition import TreelengthConnected
            sage: G = graphs.CycleGraph(4)
            sage: TreelengthConnected(G).get_length()
            2
            sage: TreelengthConnected(G, k=2).get_length()
            2
            sage: TreelengthConnected(G, k=1).get_length()
            Traceback (most recent call last):
            ...
            ValueError: no tree decomposition with length <= 1 was found

        TESTS::

            sage: TreelengthConnected(Graph()).get_length()
            0
        """
        if self.k_is_defined and not self.leq_k:
            raise ValueError("no tree decomposition with length <= {} was found".format(self.k))
        return self.length

    def is_less_than_k(self):
        """
        Return whether a tree decomposition with length at most `k` was found.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.tree_decomposition import TreelengthConnected
            sage: G = graphs.CycleGraph(4)
            sage: TreelengthConnected(G, k=1).is_less_than_k()
            False
            sage: TreelengthConnected(G, k=2).is_less_than_k()
            True
            sage: TreelengthConnected(G).is_less_than_k()
            Traceback (most recent call last):
            ...
            ValueError: parameter 'k' has not been specified

        TESTS::

            sage: TreelengthConnected(Graph(), k=1).is_less_than_k()
            True
        """
        if self.k_is_defined:
            return self.leq_k
        raise ValueError("parameter 'k' has not been specified")


def treelength(G, k=None, certificate=False):
    r"""
    Compute the treelength of `G` (and provide a decomposition).

    The *length* of a tree decomposition, as proposed in [DG2006]_, is the
    maximum *diameter* in `G` of its bags, where the diameter of a bag `X_i` is
    the largest distance in `G` between the vertices in `X_i` (i.e., `\max_{u, v
    \in X_i} \dist_G(u, v)`). The *treelength* `tl(G)` of a graph `G` is the
    minimum length among all possible tree decompositions of `G`.
    See the documentation of the
    :mod:`~sage.graphs.graph_decompositions.tree_decomposition` module for more
    details.

    INPUT:

    - ``G`` -- a sage Graph

    - ``k`` -- integer (default: ``None``); indicates the length to be
      considered. When `k` is an integer, the method checks that the graph has
      treelength `\leq k`. If `k` is ``None`` (default), the method computes the
      optimal treelength.

    - ``certificate`` -- boolean (default: ``False``); whether to also return
      the tree-decomposition itself

    OUTPUT:

    ``G.treelength()`` returns the treelength of `G`. When `k` is specified, it
    returns ``False`` when no tree-decomposition of length `\leq k` exists or
    ``True`` otherwise. When ``certificate=True``, the tree-decomposition is
    also returned.

    ALGORITHM:

    This method virtually explores the graph of all pairs ``(vertex_cut,
    connected_component)``, where ``vertex_cut`` is a vertex cut of the graph of
    length `\leq k`, and ``connected_component`` is a connected component of the
    graph induced by ``G - vertex_cut``.

    We deduce that the pair ``(vertex_cut, connected_component)`` is feasible
    with treelength `k` if ``connected_component`` is empty, or if a vertex
    ``v`` from ``vertex_cut`` can be replaced with a vertex from
    ``connected_component``, such that the pair ``(vertex_cut + v,
    connected_component - v)`` is feasible.

    In practice, this method decomposes the graph by its clique minimal
    separators into atoms, computes the treelength of each of atom and returns
    the maximum value over all the atoms. Indeed, we have that `tl(G) = \max_{X
    \in A} tl(G[X])` where `A` is the set of atoms of the decomposition by
    clique separators of `G`. When ``certificate == True``, the
    tree-decompositions of the atoms are connected to each others by adding
    edges with respect to the clique separators.

    .. SEEALSO::

        - :meth:`treewidth` computes the treewidth of a graph.
        - :meth:`~sage.graphs.graph_decompositions.vertex_separation.path_decomposition`
          computes the pathwidth of a graph.
        - module :mod:`~sage.graphs.graph_decompositions.vertex_separation`.
        - :meth:`~sage.graphs.graph_decompositions.clique_separators.atoms_and_clique_separators`

    EXAMPLES:

    The PetersenGraph has treelength 2::

        sage: G = graphs.PetersenGraph()
        sage: G.treelength()
        2

    Disconnected graphs have infinite treelength::

        sage: G = Graph(2)
        sage: G.treelength()
        +Infinity
        sage: G.treelength(k=+Infinity)
        True
        sage: G.treelength(k=2)
        False
        sage: G.treelength(certificate=True)
        Traceback (most recent call last):
        ...
        ValueError: the tree decomposition of a disconnected graph is not defined

    Chordal graphs have treelength 1::

        sage: G = graphs.RandomChordalGraph(30)
        sage: while not G.is_connected():
        ....:     G = graphs.RandomChordalGraph(30)
        sage: G.treelength()
        1

    Cycles have treelength `\lceil n/3 \rceil`::

        sage: [graphs.CycleGraph(n).treelength() for n in range(3, 11)]
        [1, 2, 2, 2, 3, 3, 3, 4]

    TESTS:

    Check that the decomposition by clique separators is valid::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import TreelengthConnected
        sage: from sage.graphs.graph_decompositions.tree_decomposition import is_valid_tree_decomposition
        sage: G = graphs.StarGraph(3)
        sage: G.subdivide_edges(G.edges(sort=False), 2)
        sage: G = G.cartesian_product(graphs.CycleGraph(3))
        sage: tl, T = G.treelength(certificate=True)
        sage: tl == TreelengthConnected(G).get_length()
        True
        sage: is_valid_tree_decomposition(G, T)
        True

    Corner cases::

        sage: Graph().treelength()
        0
        sage: Graph().treelength(certificate=True)
        (0, Tree decomposition: Graph on 0 vertices)
        sage: Graph(1).treelength()
        0
        sage: Graph(1).treelength(k=0)
        True
        sage: Graph(1).treelength(certificate=True)
        (0, Tree decomposition: Graph on 1 vertex)
        sage: Graph(1).treelength(k=0, certificate=True)
        (True, Tree decomposition: Graph on 1 vertex)
        sage: G = graphs.PathGraph(2)
        sage: G.treelength()
        1
        sage: G.treelength(k=0)
        False
        sage: G.treelength(certificate=True)
        (1, Tree decomposition of Path graph: Graph on 1 vertex)
        sage: G.treelength(certificate=True, k=0)
        (False, None)
        sage: G.treelength(certificate=True, k=1)
        (True, Tree decomposition of Path graph: Graph on 1 vertex)
        sage: G.treelength(certificate=True, k=0)
        (False, None)
        sage: G.treelength(k=-1)
        Traceback (most recent call last):
        ...
        ValueError: k(=-1) must be a nonnegative integer
    """
    if G.is_directed():
        raise ValueError("this method is defined for undirected graphs only")
    if k is not None and k < 0:
        raise ValueError("k(={}) must be a nonnegative integer".format(k))

    cdef str name = "Tree decomposition"
    if G.name():
        name += " of {}".format(G.name())

    # Corner cases
    from sage.graphs.graph import Graph
    if G.order() <= 1:
        answer = 0 if k is None else True
        if certificate:
            if G:
                answer = answer, Graph({Set(G): []}, format='dict_of_lists', name=name)
            else:
                answer = answer, Graph(name=name)
        return answer
    if not G.is_connected():
        if certificate:
            raise ValueError("the tree decomposition of a disconnected graph is not defined")
        elif k is None:
            return +Infinity
        else:
            return k is Infinity
    if k == 0:
        return (False, None) if certificate else False
    if not certificate and G.is_chordal():
        return 1 if k is None else True

    # We decompose the graph by clique minimal separators into atoms and solve
    # the problem on each of them
    atoms, cliques = G.atoms_and_clique_separators()

    if not cliques:
        # We have a single atom
        TC = TreelengthConnected(G, k=k, certificate=certificate)
        if certificate:
            if k is None:
                return TC.get_length(), TC.get_tree_decomposition()
            elif TC.is_less_than_k():
                return True, TC.get_tree_decomposition()
            else:
                return False, None
        if k is None:
            return TC.get_length()
        return TC.is_less_than_k()

    # As some atoms might be isomorphic, we use a dictionary keyed by immutable
    # copies of canonical graphs to store intermediate results.
    cdef dict data = dict()
    cdef list result = []
    cdef int tl = 1  # The graph is connected and of order at least 2
    cdef dict certif_inv
    cdef dict perm

    for atom in atoms:

        ga = G.subgraph(atom)
        if ga.is_clique():
            if certificate:
                result.append(Graph({Set(atom): []}, format='dict_of_lists'))
            continue

        gc, certif = ga.canonical_label(certificate=True)
        gci = gc.copy(immutable=True)

        if gci in data:
            # We already solved the problem for an isomorphic atom.
            if certificate:
                # We deduce the solution for this atom
                certif_inv = {i: u for u, i in certif.items()}
                perm = {u: Set([certif_inv[i] for i in u]) for u in data[gci]}
                result.append(data[gci].relabel(perm=perm, inplace=False, immutable=False))
            continue

        # We solve the problem for this atom and store the result
        TC = TreelengthConnected(gci, k=k, certificate=certificate)
        if certificate:
            T = TC.get_tree_decomposition()
            data[gci] = T
            certif_inv = {i: u for u, i in certif.items()}
            perm = {u: Set([certif_inv[i] for i in u]) for u in T}
            result.append(T.relabel(perm=perm, inplace=False, immutable=False))
        if k is None:
            tl = max(tl, TC.get_length())
        elif not TC.is_less_than_k():
            return False if not certificate else (False, None)

    if not certificate:
        if k is None:
            return tl
        return True

    # We now build the tree decomposition of the graph by connecting the tree
    # decompositions of its atoms.
    T = _from_tree_decompositions_of_atoms_to_tree_decomposition(result, cliques)

    # The tree-decomposition may contain a lot of useless nodes.
    # We merge all edges between two sets S,S' where S is a subset of S'
    T = reduced_tree_decomposition(T)
    T.name(name)
    if k is None:
        return tl, T
    return True, T
