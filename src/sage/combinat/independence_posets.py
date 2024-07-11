"""
Independence Poset

This is an implementation of independence posets.

It includes functions to compute maximal orthogonal pairs,
tight orthogonal pairs, rowmotion, labelling,
and generating Galois graphs of lattices.

REFERENCES:

- [TW2017]_

- [TW2018]_

AUTHORS:

- Nathan Williams (2022-03-24): 1.0
- Deep Desai (2022-03-24): 1.0
- Ariba Tahsin (2022-03-24): 1.0
- Ryan Lofton (2022-03-24): 1.0
- Ramesh Kanakala (2022-03-24): 1.0
- Savitha Venkatesh (2022-03-24): 1.0
"""
# ****************************************************************************
#       Copyright (C) 2013 Nathan Williams <Nathan.Williams1@utdallas.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.combinat.posets.posets import Poset
from sage.combinat.subset import Subsets
from sage.graphs.digraph import DiGraph
from sage.sets.set import Set


##############################################################
# code for maximal orthogonal pairs
##############################################################


def orthogonal_pairs(G):
    """
    Iterate over all orthogonal pairs of `G`.

    An orthogonal pair is a pair
    of sets `(X, Y)` such that there is no edge
    from a node in `X` to a node in `Y`.

    INPUT:

    - `G` -- an acyclic directed graph

    OUTPUT:

    all orthogonal pairs of `G``

    .. NOTE::

        This is a slow, brute force implementation.

    EXAMPLES:

    The orthogonal pairs for a small graph ::

        sage: from sage.combinat.independence_posets import orthogonal_pairs
        sage: G = DiGraph([[1, 2], [(1, 2)]])
        sage: list(orthogonal_pairs(G))
        [({}, {}), ({}, {1}), ({}, {2}), ({}, {1, 2}),
         ({1}, {}), ({2}, {}), ({2}, {1}), ({1, 2}, {})]

    The orthogonal pairs for a complete graph ::

        sage: G = digraphs.TransitiveTournament(3)
        sage: list(orthogonal_pairs(G))
        [({}, {}), ({}, {0}), ({}, {1}), ({}, {2}), ({}, {0, 1}), ({}, {0, 2}),
         ({}, {1, 2}), ({}, {0, 1, 2}), ({0}, {}), ({1}, {}), ({1}, {0}),
         ({2}, {}), ({2}, {0}), ({2}, {1}), ({2}, {0, 1}), ({0, 1}, {}),
         ({0, 2}, {}), ({1, 2}, {}), ({1, 2}, {0}), ({0, 1, 2}, {})]
    """
    vertices = Set(G)
    for s in Subsets(vertices):
        bads = Set(j for i in s for j in G.neighbor_out_iterator(i))
        for t in Subsets(vertices.difference(s + bads)):
            yield (s, t)


def poset_of_orthogonal_pairs(G):
    """
    Return the poset of orthogonal pairs.

    INPUT:

    - ``G`` -- an acyclic directed graph

    OUTPUT:

    a poset of orthogonal pairs of g with the inclusion
    relations on both elements

    EXAMPLES:

    The poset of orthogonal pairs for a small graph ::

        sage: from sage.combinat.independence_posets import poset_of_orthogonal_pairs
        sage: G = DiGraph([[1, 2], [(1, 2)]])
        sage: poset_of_orthogonal_pairs(G)
        Finite poset containing 8 elements

    The poset of orthogonal pairs for a completely connected graph ::

        sage: G = digraphs.TransitiveTournament(3)
        sage: poset_of_orthogonal_pairs(G)
        Finite poset containing 20 elements
    """
    def inc(a, b) -> bool:
        return a[0].issubset(b[0]) and a[1].issubset(b[1])

    return Poset([list(orthogonal_pairs(G)), inc])


def mops(G) -> list:
    """
    Return the maximal orthogonal pairs in ``G``.

    INPUT:

    - ``G`` -- an acyclic directed graph

    OUTPUT:

    list of pairs of sets `(X, Y)` containing the maximal elements in a poset
    of orthogonal pairs with respect to the condition that
    no edges run from `X` to `Y`.

    EXAMPLES:

    The maximal orthogonal pairs for a connected two-vertex graph ::

        sage: from sage.combinat.independence_posets import mops
        sage: G = DiGraph([[1, 2], [(1, 2)]])
        sage: mops(G)
        [({}, {1, 2}), ({2}, {1}), ({1, 2}, {})]

    The maximal orthogonal pairs for a completely connected 3 vertex graph ::

        sage: G = digraphs.TransitiveTournament(3)
        sage: mops(G)
        [({}, {0, 1, 2}), ({2}, {0, 1}), ({1, 2}, {0}), ({0, 1, 2}, {})]
    """
    return poset_of_orthogonal_pairs(G).maximal_elements()


def poset_of_mops(G):
    """
    Return the poset of maximal orthogonal pairs in ``G``.

    INPUT:

    - ``G`` -- an acyclic directed graph

    OUTPUT:

    a poset of the maximal orthogonal pairs of ``G``
    with the inclusion relation for only the first element

    EXAMPLES:

    The poset of maximal orthogonal pairs for a connected two-vertex graph ::

        sage: from sage.combinat.independence_posets import poset_of_mops
        sage: G = DiGraph([[1, 2], [(1, 2)]])
        sage: poset_of_mops(G)
        Finite poset containing 3 elements

    The poset of maximal orthogonal pairs for
    a completely connected 3 vertex graph ::

        sage: G = digraphs.TransitiveTournament(3)
        sage: poset_of_mops(G)
        Finite poset containing 4 elements
    """
    return Poset([mops(G), lambda a, b: a[0].issubset(b[0])])


def complete_mop(G, S):
    r"""
    Return a list of maximal orthogonal pairs in ``G`` extending ``S``.

    INPUT:

    - ``G`` -- an acyclic directed graph

    - ``S`` -- a subset of the graph ``G`` of type ``Set``

    OUTPUT:

    the list of maximal elements in the poset
    of orthogonal pairs `(X, Y)` such that `X \subset S`

    EXAMPLES:

    Example of the ``complete_mop`` ::

        sage: from sage.combinat.independence_posets import complete_mop
        sage: G = DiGraph([[1..5], lambda i, j: i < j])
        sage: S = Set([2..4])
        sage: complete_mop(G, S)
        [({2, 3, 4}, {1})]

    .. NOTE::

        This is a slow, brute force implementation.
    """
    edges = set(G.edges(labels=False))
    vertices = Set(G)
    sets = [(S, t) for t in Subsets(vertices.difference(S))
            if not any((i, j) in edges for i in S for j in t)]

    def inc(a, b) -> bool:
        return a[0].issubset(b[0]) and a[1].issubset(b[1])

    return [u for u in sets if not any(inc(u, v) for v in sets if v != u)]


##############################################################
# code for left modular lattices
##############################################################


def join_irr_label(L, c, j):
    r"""
    Return the label of the join-irreducible element ``j``.

    INPUT:

    ``L`` -- a left-modular lattice

    ``c`` -- a left-modular chain in ``L``

    ``j`` -- join-irreducible element

    OUTPUT:

    The label is the minimum index `i` such that `j < x_i`, where `x_i \in c`.

    EXAMPLES::

        sage: from sage.combinat.independence_posets import join_irr_label
        sage: L = LatticePoset(([1,2,3,4,5,6],
        ....:     [[1,2],[1,3],[3,4],[4,5],[2,5],[2,6],[6,5],[2,4]]))
        sage: chain = [1,3,4,5]
        sage: join_irr_label(L, chain, 5)
        3
    """
    for i in range(len(c)):
        if L.le(j, c[i]):
            return i


def meet_irr_label(L, c, m):
    r"""
    Return the label of the meet-irreducible element ``m``.

    INPUT:

    ``L`` -- a left-modular lattice

    ``c`` -- a left-modular chain in ``L``

    ``m`` -- meet-irreducible element

    OUTPUT:

    The label is the maximum index `i` such that `m > x_{i-1}`, where `x_{i-1} \in c`.

    EXAMPLES::

        sage: from sage.combinat.independence_posets import meet_irr_label
        sage: L = LatticePoset(([1,2,3,4,5,6],
        ....:     [[1,2],[1,3],[3,4],[4,5],[2,5],[2,6],[6,5],[2,4]]))
        sage: chain = [1,3,4,5]
        sage: meet_irr_label(L, chain, 5)
        4
    """
    for i in range(len(c) - 1, -1, -1):
        if L.ge(m, c[i]):
            return i + 1


def downward_labels(L, c, z):
    r"""
    Return the list of downward labels of ``z`` in the lattice ``L``.

    INPUT:

    - ``L`` -- a trim lattice

    - ``c`` -- a left-modular chain in ``L``

    - ``z`` -- an element in the lattice ``L``

    OUTPUT:

    the list of downward labels for z, defined as
    `\min\{\beta_J(j) | j \in J, y \vee j = z\}` where
    `J` is the set of all join-irreducibles in `L` and
    `\beta_J(j)` is the join-irreducible label of `j`
    for all `y \in L, y \leq z`

    EXAMPLES::

        sage: from sage.combinat.independence_posets import downward_labels
        sage: L = LatticePoset(([1,2,3,4,5,6],
        ....:     [[1,2],[1,3],[3,4],[4,5],[2,5],[2,6],[6,5],[2,4]]))
        sage: chain = [1,2,4,5]
        sage: downward_labels(L, chain, 4)
        [2, 1]
        sage: L = LatticePoset(([1,2,3,4,5,6],
        ....:     [[1,2],[1,3],[3,4],[4,5],[2,5],[2,6],[6,5],[2,4]]))
        sage: chain = [1,2,4,5]
        sage: downward_labels(L, chain, 5)
        [2, 3]
    """
    J = L.join_irreducibles()
    JL = {j: join_irr_label(L, c, j) for j in J}
    cov = L.lower_covers(z)
    return [min(JL[j] for j in J if L.join([y, j]) == z) for y in cov]


def upward_labels(L, c, z):
    r"""
    INPUT:

    ``L`` -- a trim lattice L

    ``c`` -- a left-modular chain c

    ``z`` -- an element z in the lattice ``L``

    OUTPUT:

    the list of upward labels for ``z``, defined as
    `\max\{\beta_M(m) | m \in M, z \wedge m = y\}` where `M` is
    the set of all meet-irreducibles in `L` and `\beta_M(m)` is the
    meet-irreducible label of `m` for all `y \in L, y \geq z`

    EXAMPLES::

        sage: from sage.combinat.independence_posets import upward_labels
        sage: L = LatticePoset(([1,2,3,4,5,6],
        ....:     [[1,2],[1,3],[3,4],[4,5],[2,5],[2,6],[6,5],[2,4]]))
        sage: chain = [1,2,4,5]
        sage: upward_labels(L, chain, 1)
        [1, 2]
        sage: L = LatticePoset(([1,2,3,4,5,6],
        ....:     [[1,2],[1,3],[3,4],[4,5],[2,5],[2,6],[6,5],[2,4]]))
        sage: chain = [1,2,4,5]
        sage: upward_labels(L, chain, 2)
        [3, 2]
    """
    M = L.meet_irreducibles()
    ML = {m: meet_irr_label(L, c, m) for m in M}
    cov = L.upper_covers(z)
    return [max(ML[m] for m in M if L.meet([y, m]) == z) for y in cov]


def cover_label(L, c, edge):
    r"""
    INPUT:

    ``L`` -- a trim lattice L

    ``c`` -- a left-modular chain c

    ``edge`` -- a covering relation ``[x,y]`` where `x \leq y`

    OUTPUT:

    the label for the covering relation ``edge`` in ``L``

    EXAMPLES::

        sage: from sage.combinat.independence_posets import cover_label
        sage: L = LatticePoset(([1,2,3,4,5,6],
        ....:     [[1,2],[1,3],[3,4],[4,5],[2,5],[2,6],[6,5],[2,4]]))
        sage: chain = [1,2,4,5]
        sage: cover_label(L, chain, [2, 4])
        2
        sage: L = LatticePoset(([1,2,3,4,5,6],
        ....:     [[1,2],[1,3],[3,4],[4,5],[2,5],[2,6],[6,5],[2,4]]))
        sage: chain = [1,2,4,5]
        sage: cover_label(L, chain, [2, 6])
        3
    """
    x, y = edge
    J = L.join_irreducibles()
    JL = {j: join_irr_label(L, c, j) for j in J}
    return min(JL[j] for j in J if L.join([x, j]) == y)


##############################################################
# tight orthogonal pairs (tops)
##############################################################
def complete_top(G, u):
    """
    INPUT:

    - ``G`` -- a finite acyclic directed graph

    - ``u`` -- an independent set in ``G``

    OUTPUT:

    returns a tuple (``m``, ``u``) independent sets in ``G``
    which are orthogonal (there is no edge from ``m`` to ``u`` in ``G``)
    and tight (if an element in ``m`` is increased (removed and replaced
    by a larger element) or any element of ``u`` is decreased,
    or a new element is added into either ``m`` or ``u``
    then the result will no longer be an orthogonal pair)

    EXAMPLES::

        sage: from sage.combinat.independence_posets import complete_top
        sage: G = DiGraph([[1, 2, 3, 4], [[4,3], [3,2], [2,1]]])
        sage: complete_top(G, Set([2, 4]))
        ({1}, {2, 4})
        sage: complete_top(G, Set([1]))
        ({4}, {1})
        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9],
        ....:     [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6],
        ....:      [4,7], [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: complete_top(G, Set([1,4]))
        ({8, 3, 7}, {1, 4})
        sage: complete_top(G, Set([5, 6, 9]))
        ({1, 4}, {9, 5, 6})
    """
    m = Set()
    for k in G.topological_sort():
        ngk_m = any(g in m for g in G.neighbor_in_iterator(k))
        nlk_u = any(g in u for g in G.neighbor_out_iterator(k))
        if not ngk_m and not nlk_u and k not in u:
            m = m.union(Set([k]))
    return (m, u)


def top_to_mop(G, du):
    """
    Transform a tight orthogonal pair into a maximal orthogonal pair.

    INPUT:

    - ``G`` -- a finite acyclic directed graph

    - ``du`` -- a tight orthogonal pair (``d``, ``u``)

    OUTPUT:

    a maximal orthogonal pair corresponding to the given tight orthogonal pair (``d``, ``u``)

    EXAMPLES::

        sage: from sage.combinat.independence_posets import complete_top, top_to_mop
        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9],
        ....:     [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6],
        ....:      [4,7], [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: du = complete_top(G, Set([4, 7]))
        sage: top_to_mop(G, du)
        ({1, 3, 6, 8, 9}, {4, 7})
    """
    X, Y = du
    setG = Set(G)
    for b in setG.difference(Y):
        if not any(g in Y for g in G.neighbor_out_iterator(b)):
            X = X.union(Set([b]))
    for a in setG.difference(X):
        if not any(g in X for g in G.neighbor_in_iterator(a)):
            Y = Y.union(Set([a]))
    return (X, Y)


def maximal_top(G):
    """
    Return the maximal tight orthogonal pair.

    INPUT:

    - ``G`` -- a finite acyclic directed graph

    OUTPUT:

    a tuple (``m``, ``u``) where ``m`` is the set which is
    tight and orthogonal to an empty ``u``.

    EXAMPLES::

        sage: from sage.combinat.independence_posets import maximal_top
        sage: G = DiGraph([[1, 2, 3, 4], [[4,3], [3,2], [2,1]]])
        sage: maximal_top(G)
        ({2, 4}, {})
        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9],
        ....:     [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6],
        ....:      [4,7], [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: maximal_top(G)
        ({1, 4, 5, 6, 9}, {})
    """
    return complete_top(G, Set())


def minimal_top(G):
    """
    Return the minimal tight orthogonal pair.

    INPUT:

    - ``G`` -- a finite acyclic directed graph

    OUTPUT:

    a tuple (``m``, ``u``) where ``u`` is the set which
    is tight and orthogonal to an empty ``m``.

    EXAMPLES::

        sage: from sage.combinat.independence_posets import minimal_top
        sage: G = DiGraph([[1, 2, 3, 4], [[4,3], [3,2], [2,1]]])
        sage: minimal_top(G)
        ({}, {1, 3})
        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9],
        ....:     [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6],
        ....:      [4,7], [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: minimal_top(G)
        ({}, {1, 4, 5, 6, 9})
    """
    m = Set()
    for k in reversed(G.topological_sort()):
        if not any(g in m for g in G.neighbor_out_iterator(k)):
            m = m.union(Set([k]))
    return (Set(), m)


def flip(G, du, j):
    r"""
    Perform a flip if `j \in D` or `j \in U`.

    INPUT:

    - ``G`` -- a finite acyclic directed graph

    - ``du`` -- a tight orthogonal pair (D, U)

    - ``j`` -- an element `j` of ``G``

    OUTPUT:

    a new tight orthogonal pair ``(d, u)``

    EXAMPLES:

    Flipping on an element in U ::

        sage: from sage.combinat.independence_posets import complete_top, flip
        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9],
        ....:     [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6], [4,7],
        ....:      [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: du = complete_top(G, Set([1,4]))
        sage: flip(G, du, 4)
        ({8, 3, 4}, {2})

    Flipping on an element in D ::

        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9],
        ....:     [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6], [4,7],
        ....:      [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: du = complete_top(G, Set([1,4]))
        sage: flip(G, du, 3)
        ({9, 5, 6}, {3, 4})

    Flipping on an element in neither D nor U ::

        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9],
        ....:     [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6], [4,7],
        ....:      [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: du = complete_top(G, Set([1,4]))
        sage: flip(G, du, 5)
        ({8, 3, 7}, {1, 4})
    """
    a, b = du
    if j not in b and j not in a:
        return du
    P = Poset(G)
    if j in b:
        a2 = Set([k for k in a if not P.ge(k, j)] + [j])
        b2 = Set(k for k in b if not P.le(k, j))
    elif j in a:
        a2 = Set(k for k in a if not P.ge(k, j))
        b2 = Set([k for k in b if not P.le(k, j)] + [j])

    for k in P.linear_extension():
        if not P.le(k, j) \
           and not Set(G.neighbor_out_iterator(k)).intersection(b) \
           and not Set(G.neighbor_in_iterator(k)).intersection(a2) \
                and (k not in b2):
            a2 = a2.union(Set([k]))
    for k in reversed(P.linear_extension()):
        if not P.ge(k, j) \
           and not Set(G.neighbor_in_iterator(k)).intersection(a) \
           and not Set(G.neighbor_out_iterator(k)).intersection(b2) \
                and (k not in a2):
            b2 = b2.union(Set([k]))
    return (a2, b2)


def flip_up(G, du, j):
    r"""
    Perform a flip only if `j \in U`.

    INPUT:

    - ``G`` -- a finite acyclic directed graph

    - ``du`` -- a tight orthogonal pair (D, U)

    - ``j`` -- an element j in U

    OUTPUT:

    a new tight orthogonal pair ``(d, u)``

    EXAMPLES:

    Flipping on an element in U ::

        sage: from sage.combinat.independence_posets import complete_top, flip_up
        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9],
        ....:     [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6], [4,7],
        ....:      [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: du = complete_top(G, Set([1,4]))
        sage: flip_up(G, du, 4)
        ({8, 3, 4}, {2})

    Flipping on an element in D ::

        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9],
        ....:     [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6], [4,7],
        ....:      [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: du = complete_top(G, Set([1,4]))
        sage: flip_up(G, du, 3)
        ({8, 3, 7}, {1, 4})

    Flipping on an element in neither D nor U ::

        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9],
        ....:     [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6], [4,7],
        ....:      [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: du = complete_top(G, Set([1,4]))
        sage: flip_up(G, du, 5)
        ({8, 3, 7}, {1, 4})
    """
    if j not in du[1]:
        return du
    return flip(G, du, j)


def flips(G, du, s):
    r"""
    Perform a sequence of flips.

    INPUT:

    - ``G`` -- a directed acyclic graph

    - ``du`` -- a tight orthogonal pair (D, U)

    - ``s`` -- a sequence of elements `j` in ``G``

    OUTPUT:

    a new top ``(d, u)`` after flipping on each element in ``s``

    This applies `\displaystyle{ \left( \prod_{j \in s} \text{flip}_j \right)(D, U) }`.

    EXAMPLES::

        sage: from sage.combinat.independence_posets import complete_top, flip, flips
        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9],
        ....:     [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6], [4,7],
        ....:      [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: du = complete_top(G, Set([1, 4]))
        sage: a = flip(G, du, 4)
        sage: b = flip(G, a, 2)
        sage: flip(G, b, 3)
        ({2, 6, 7}, {3})
        sage: flips(G, du, [4, 2, 3])
        ({2, 6, 7}, {3})
    """
    s_iter = s if isinstance(s, list) else [s]
    dv = tuple(du)
    for i in s_iter:
        dv = flip(G, dv, i)
    return dv


def make_flip_poset(G):
    """
    Return the independence poset corresponding to the graph ``G``.

    INPUT:

    - ``G`` -- a finite directed acyclic graph

    EXAMPLES::

        sage: from sage.combinat.independence_posets import make_flip_poset
        sage: G = DiGraph([[1, 2, 3, 4], [[4, 3], [3, 2], [2, 1]]])
        sage: make_flip_poset(G)
        Finite poset containing 8 elements

    The below graph generates a Tamari lattice with 14 elements ::

        sage: G = DiGraph([[1, 2, 3, 4, 5, 6],
        ....:     [[6, 5], [6, 4], [5, 4], [5, 3], [5, 2],
        ....:      [4, 2], [4, 1], [3, 2], [2, 1]]])
        sage: make_flip_poset(G)
        Finite poset containing 14 elements
    """
    edges = []
    els = Set()
    last = []
    new = [minimal_top(G)]
    while new:
        last = new
        els = els.union(Set(last))
        new = []
        for a in last:
            for i in a[1]:
                b = flip(G, a, i)
                edges.append([a, b])
                if b not in new:
                    new = new + [b]
    return Poset((list(els), edges))


def rowmotion(L, x):
    """
    Perform the rowmotion of element `x` in lattice `L`.

    INPUT:

    `L` -- trim lattice

    `x` -- element to apply rowmotion on

    OUTPUT:

    Rowmotion (with respect with a particular labeling on the lattice `L`)
    on an element `x` is the unique element `y` such that the set of labels
    of upper covers `s > y` equals the set of labels of lower covers `s < x`

    EXAMPLES:

    Example of rowmotion on an independence poset ::

        sage: from sage.combinat.independence_posets import make_flip_poset, rowmotion
        sage: G = DiGraph([[1, 2, 3, 4], [[4, 3], [3, 2], [2, 1]]])
        sage: L = make_flip_poset(G)
        sage: x = L.list()[3]
        sage: rowmotion(L, x)
        ({1, 3}, {4})
    """
    lower_x = x[0]
    return next(elt for elt in L if elt[1] == lower_x)


def toggle(G, x):
    """
    Toggle the vertex ``x`` in the graph ``G``.

    INPUT:

    ``G`` -- a finite directed acyclic graph

    ``x`` -- an element in ``G``

    OUTPUT:

    the directed acyclic graph given by reversing all edges incident to ``x`` in ``G``

    EXAMPLES:

    Example of a usage of toggle.
    We take a graph, ``G`` and compute the independence poset
    on ``G`` and ``tog_3(G)`` ::

        sage: from sage.combinat.independence_posets import toggle, make_flip_poset
        sage: G = DiGraph([[1, 2, 3], [[3, 2], [2, 1]]])
        sage: G_3 = toggle(G, 3)
        sage: P_3 = make_flip_poset(G_3)
        sage: P_3
        Finite poset containing 5 elements
    """
    H = DiGraph(G)
    for i in H.edges():
        if x in i:
            H.reverse_edge(i)
    return H


def toggles(G, s):
    """
    Toggle all elements of ``s`` in ``G``.

    INPUT:

    - ``G`` -- a finite directed acyclic graph

    - ``s`` -- a sequence of elements

    OUTPUT:

    This toggles all elements `x` in the graph `G` (i.e. this reverses all edges
    incident to `x`) as they appear in the sequence `s`.

    EXAMPLES:

    Example of a usage of toggles. We take a graph, `G` and compute the
    independence poset on `G` and `tog_4(tog_2(G))` in two different ways:
    once manually using ``toggle`` and once using ``toggles`` ::

        sage: from sage.combinat.independence_posets import toggle, toggles, make_flip_poset
        sage: G = DiGraph([[1, 2, 3, 4, 5, 6], [[6, 5], [6, 4], [5, 4], [5, 3], [5, 2], [4, 2], [4, 1], [3, 2], [2, 1]]])
        sage: P = make_flip_poset(G)
        sage: G_2 = toggle(G, 2)
        sage: G_24_1 = toggle(G_2, 4)
        sage: G_24_2 = toggles(G, [2, 4])
        sage: G_24_1 == G_24_2
        True

    ``toggles`` will also accept a singular element for ``s``;
    in which case it will be the same thing as a ``toggle`` ::

        sage: G = DiGraph([[1, 2, 3], [[3, 2], [2, 1]]])
        sage: P_1 = toggles(G, 2)
        sage: P_2 = toggle(G, 2)
        sage: P_1 == P_2
        True
    """
    s_iter = s if isinstance(s, list) else [s]
    H = DiGraph(G)
    for i in s_iter:
        H = toggle(H, i)
    return H
