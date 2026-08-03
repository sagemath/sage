r"""
The fundamental theorem of finite semidistributive lattices

This is a result shown in [RST2024]_, stating that a finite lattice is
semidistributive if and only if it is isomorphic to some lattice of
maximal orthogonal pairs of a finite two-acyclic factorization system.

This module implements the two key concepts related to the theorem:
the lattice of maximal orthogonal pairs, and two-acyclic factorization
systems.
"""
from sage.combinat.posets.lattices import LatticePoset
from sage.graphs.digraph import DiGraph


def right_orthogonal(G, X):
    r"""
    Return the right orthogonal of a given set of vertices.

    If `X` is a set of vertices of `G`, the right orthogonal of `X` is defined as
    the set of vertices of `G` that are not in `X` or out-neighbors of any vertices
    in `X`.

    This set is defined on page 3 of [RST2024]_, but left nameless - we have opted
    for "right orthogonal" as a placeholder term.

    INPUT:

    - ``G`` -- DiGraph

    - ``X`` -- Set; a subset of vertices

    OUTPUT: The right orthogonal of ``X`` as a set of vertices.

    EXAMPLES:

    Consider the following graph with 2 vertices and 1 edge. `0` points to `1`,
    so the right orthogonal of {0} is empty. However, `1` does not point to `0`,
    so the the right orthogonal of {1} contains `0`::

        sage: from sage.combinat.posets.semidistributivity import right_orthogonal
        sage: G = DiGraph([(0, 1)])
        sage: right_orthogonal(G, {1})
        {0}
        sage: right_orthogonal(G, {0})
        set()

    A directional 4-cycle::

        sage: G = DiGraph([(0, 1), (1, 2), (2, 3), (3, 0)])
        sage: right_orthogonal(G, {0})
        {2, 3}
        sage: right_orthogonal(G, {1, 2})
        {0}
        sage: right_orthogonal(G, {1, 2, 3})
        set()

    A bidirectional 4-cycle. Note how adding edges removes elements from the right
    orthogonal sets::

        sage: G = DiGraph(graphs.CycleGraph(4))
        sage: right_orthogonal(G, {0})
        {2}
        sage: right_orthogonal(G, {1, 2})
        set()

    .. SEEALSO::

        - :meth:`neighbors_out()`
        - Dual function: :meth:`left_orthogonal()`

    REFERENCE:

    - [RST2024]_

    TESTS:

    The right orthogonal of the empty set is always the entire graph::

        sage: for i in range(20):                                                      # needs sage.numerical.mip
        ....:     G = digraphs.RandomDirectedGNP(10, .3)
        ....:     assert right_orthogonal(G, set()) == set(G)

    The right orthogonal of the entire graph is always the empty set::

        sage: for i in range(20):                                                      # needs sage.numerical.mip
        ....:     G = digraphs.RandomDirectedGNP(10, .3)
        ....:     assert right_orthogonal(G, set(G)) == set()
    """
    l = set(G)
    for x in X:
        l.discard(x)
        for y in G.neighbor_out_iterator(x):
            l.discard(y)
    return l

def left_orthogonal(G, X):
    r"""
    Return the left orthogonal of a given set of vertices.

    If `X` is a set of vertices of `G`, the left orthogonal of `X` is defined as
    the set of vertices of `G` that are not in `X` or in-neighbors of any vertices
    in `X`.

    This set is defined on page 3 of [RST2024]_, but left nameless - we have opted
    for "left orthogonal" as a placeholder term.

    INPUT:

    - ``G`` -- DiGraph

    - ``X`` -- Set; a subset of vertices

    OUTPUT: The left orthogonal of ``X`` as a set of vertices.

    EXAMPLES:

    Consider the following graph with 2 vertices and 1 edge. `0` points to `1`,
    so the left orthogonal of {1} is empty. However, `1` does not point to `0`,
    so the the left orthogonal of {0} contains `1`::

        sage: from sage.combinat.posets.semidistributivity import left_orthogonal
        sage: G = DiGraph([(0, 1)])
        sage: left_orthogonal(G, {0})
        {1}
        sage: left_orthogonal(G, {1})
        set()

    A directional 4-cycle::

        sage: G = DiGraph([(0, 1), (1, 2), (2, 3), (3, 0)])
        sage: left_orthogonal(G, {0})
        {1, 2}
        sage: left_orthogonal(G, {1, 2})
        {3}
        sage: left_orthogonal(G, {1, 2, 3})
        set()

    A bidirectional 4-cycle. Note how adding edges removes elements from the left
    orthogonal sets::

        sage: G = DiGraph(graphs.CycleGraph(4))
        sage: left_orthogonal(G, {0})
        {2}
        sage: left_orthogonal(G, {1, 2})
        set()

    .. SEEALSO::

        - :meth:`neighbors_in()`
        - Dual function: :meth:`right_orthogonal()`

    REFERENCE:

    - [RST2024]_

    TESTS:

    The left orthogonal of the empty set is always the entire graph::

        sage: for i in range(20):                                                      # needs sage.numerical.mip
        ....:     G = digraphs.RandomDirectedGNP(10, .3)
        ....:     assert left_orthogonal(G, set()) == set(G)

    The left orthogonal of the entire graph is always the empty set::

        sage: for i in range(20):                                                      # needs sage.numerical.mip
        ....:     G = digraphs.RandomDirectedGNP(10, .3)
        ....:     assert left_orthogonal(G, set(G)) == set()
    """
    l = set(G)
    for x in X:
        l.discard(x)
        for y in G.neighbor_in_iterator(x):
            l.discard(y)
    return l

def maximal_orthogonal_pairs_lattice(G, labels="pair"):
    r"""
    Return the lattice of maximal orthogonal pairs of ``G``.

    A maximal orthogonal pair is a pair `(X, Y)`, where X and Y are sets of
    vertices, such that `X` is the left orthogonal of `Y`, and `Y` is the right
    orthogonal of `X`. Equivalently, `X` and `Y` are disjoint with no arrow from
    any element of `X` to any element of `Y`, and maximal with this property.

    The set of all maximal orthogonal pairs, ordered by inclusion on the first
    component (or, equivalently, reverse inclusion on the second) forms a lattice.
    This lattice is sometimes written as `Pairs(\rightarrow)`.

    INPUT:

    - ``G`` -- DiGraph

    - ``labels`` -- string; either "left", "right", or "pair".
        Since each side of the pair determines the other, the elements of the lattice
        can be labeled by the first component, the second, or both.

        #. ``left`` -- label each maximal orthogonal pair by its first component.
        #. ``right`` -- label each maximal orthogonal pair by its second component.
        #. ``pair`` -- label each maximal orthogonal pair by itself.

    OUTPUT: The lattice of maximal orthogonal pairs of ``G``. The elements of the
    lattice are labeled according to the above rules. Note that the type of each
    component of a pair is ``frozenset`` and not ``set`` as it must be hashable.

    ALGORITHM:

    Using each pair already calculated, find remaining pairs by adding/removing
    vertices from each component, constructing the lattice in the process. Start
    at `({}, G)`, which corresponds to the bottom of the lattice, and iterate
    while moving up.

    We make use of two tricks to efficiently calculate both sides of each pair:

    - if `(L, R)` is a maximal orthogonal pair, and we want to add some vertex `x`
        to `L`, the right orthogonal of `L \cup {x}` can be obtained by removing
        the out-neighbors of `x` from `R`. Additionally, the resulting set is sure
        to be the right component of some pair.

    - if `R` is the right component of some pair, the left component `L` is the
        union of all sets whose right orthogonal is `R`. Combining this with the
        previous insight means we can merge all elements that give the same `R` to
        obtain the left orthogonal of the new pair.

    EXAMPLES::

        sage: from sage.combinat.posets.semidistributivity import maximal_orthogonal_pairs_lattice
        sage: G = DiGraph([(0, 1), (1, 2), (2, 3), (3, 0)])
        sage: maximal_orthogonal_pairs_lattice(G)
        Finite lattice containing 10 elements
        sage: list(maximal_orthogonal_pairs_lattice(G))
        [(frozenset(), frozenset({0, 1, 2, 3})),
        (frozenset({3}), frozenset({1, 2})),
        (frozenset({2}), frozenset({0, 1})),
        (frozenset({2, 3}), frozenset({1})),
        (frozenset({1}), frozenset({0, 3})),
        (frozenset({1, 2}), frozenset({0})),
        (frozenset({0}), frozenset({2, 3})),
        (frozenset({0, 3}), frozenset({2})),
        (frozenset({0, 1}), frozenset({3})),
        (frozenset({0, 1, 2, 3}), frozenset())]

        sage: from sage.combinat.posets.semidistributivity import right_orthogonal
        sage: from sage.combinat.posets.semidistributivity import left_orthogonal
        sage: L = maximal_orthogonal_pairs_lattice(G, labels="left")
        sage: R = maximal_orthogonal_pairs_lattice(G, labels="right")
        sage: all([right_orthogonal(G, L[i]) == set(R[i]) for i in (0..9)])
        True
        sage: all([left_orthogonal(G, R[i]) == set(L[i]) for i in (0..9)])
        True

    .. SEEALSO::

        - :meth:`right_orthogonal()`
        - :meth:`left_orthogonal()`

    REFERENCES:

    - [RST2024]_

    - [Muh2021]_

    - [TW2018]_

    TESTS::

        sage: G = DiGraph()
        sage: maximal_orthogonal_pairs_lattice(G)
        Finite lattice containing 1 elements
        sage: list(maximal_orthogonal_pairs_lattice(G))
        [(frozenset(), frozenset())]
    """
    Pairs = DiGraph()
    Pairs.add_vertex(frozenset(G))
    # dictionary of pairs, where
    # the first component is indexed by the second
    # for example, pairs[second_term] should give first_term
    pairs = {frozenset(G): frozenset()}
    next_pairs = [frozenset(G)]
    while next_pairs:
        new_pairs = []
        for rt in next_pairs:
            covering_pairs = []
            for x in G:
                if x in pairs[rt]:
                    continue
                # calculate the new right orthogonal by removing vertices
                new_rt = rt.difference(G.neighbors_out(x))
                new_rt = new_rt.difference([x])
                covering_pairs.append(new_rt)
                if new_rt in pairs:
                    # merge all left components with the same right orthogonal
                    pairs[new_rt] = pairs[new_rt].union(pairs[rt])
                else:
                    pairs[new_rt] = pairs[rt]
                    new_pairs.append(new_rt)
                pairs[new_rt] = pairs[new_rt].union(frozenset({x}))
            # generate the upper covers
            for new_rt in covering_pairs:
                if rt != new_rt:
                    Pairs.add_edge(rt, new_rt)
        next_pairs = new_pairs
    L = LatticePoset(Pairs)
    if labels == "left":
        return L.relabel(pairs)
    if labels != "right":
        return L.relabel(lambda v: (pairs[v], v))
    return L

def surjective_edges(G, loops=False):
    r"""
    Return the list of surjective edges of ``G``.

    An edge `xy` is said to be surjective if, for all `z` such that `yz` is an
    edge, `xz` is also an edge. Since any loop is trivially surjective, they are
    not counted by default.

    For more information, see :meth:`is_two_acyclic_factorization_system()`.

    INPUT:

    - ``G`` -- DiGraph

    - ``loops`` -- boolean (default: ``False``)
    ; whether to count loops

    EXAMPLES::

        sage: from sage.combinat.posets.semidistributivity import surjective_edges
        sage: G = DiGraph([(0, 1), (0, 2), (1, 3), (2, 3)], loops=True)
        sage: surjective_edges(G)
        [(1, 3), (2, 3)]
        sage: G.add_edges([(3, 1), (2, 1), (0, 0)])
        sage: surjective_edges(G)
        [(2, 1), (2, 3)]
        sage: surjective_edges(G, loops=True)
        [(0, 0), (2, 1), (2, 3)]

    .. SEEALSO::

        - Dual function: :meth:`injective_edges()`
        - :meth:`is_two_acyclic_factorization_system()`

    REFERENCE:

    - [RST2024]_

    TESTS::

        sage: G = DiGraph()
        sage: surjective_edges(G)
        []
        sage: surjective_edges(G, loops=True)
        []
    """
    E = []
    for x, y in G.edge_iterator(labels=False):
        if x != y:
            if all(G.has_edge(x, z) for z in G.neighbors_out(y)):
                E.append((x, y))
        elif loops:
            E.append((x, y))
    return E

def injective_edges(G, loops=False):
    r"""
    Return the list of injective edges of ``G``.

    An edge `yz` is said to be injective if, for all `x` such that `xy` is an
    edge, `xz` is also an edge. Since any loop is trivially injective, they are
    not counted by default.

    For more information, see :meth:`is_two_acyclic_factorization_system()`.

    INPUT:

    - ``G`` -- DiGraph

    - ``loops`` -- boolean (default: ``False``); whether to count loops

    EXAMPLES::

        sage: from sage.combinat.posets.semidistributivity import injective_edges
        sage: G = DiGraph([(0, 1), (0, 2), (1, 3), (2, 3)], loops=True)
        sage: injective_edges(G)
        [(0, 1), (0, 2)]
        sage: G.add_edges([(2, 0), (2, 1), (3, 3)])
        sage: injective_edges(G)
        [(0, 1), (2, 1)]
        sage: injective_edges(G, loops=True)
        [(0, 1), (2, 1), (3, 3)]

    .. SEEALSO::

        - Dual function: :meth:`surjective_edges()`
        - :meth:`is_two_acyclic_factorization_system()`

    REFERENCE:

    - [RST2024]_

    TESTS::

        sage: G = DiGraph()
        sage: injective_edges(G)
        []
        sage: injective_edges(G, loops=True)
        []
    """
    E = []
    for y, z in G.edge_iterator(labels=False):
        if y != z:
            if all(G.has_edge(x, z) for x in G.neighbors_in(y)):
                E.append((y, z))
        elif loops:
            E.append((y, z))
    return E

def is_two_acyclic_factorization_system(G, certificate=False):
    r"""
    Return whether ``G`` forms a two-acyclic factorization system.

    Given a DiGraph, consider the following binary reflexive relation defined on
    the set of its vertices:

    .. MATH::

        x \rightarrow y \iff x = y \text{or} xy \text{is an edge}

    This is equivalent to viewing ``G`` as the representation of some arbitrary
    binary relation. Now let `\twoheadrightarrow` and `\hookrightarrow` denote the
    previous relation restricted to surjective and injective edges, respectively.
    The triple `(rightarrow, \twoheadrightarrow, \hookrightarrow)` is said to be
    two-acyclic factorization system if:

    #. for all `x, z`, `x \rightarrow z` if and only if there exists some `y` such
        that `x \twoheadrightarrow y` and `y \hookrightarrow z`. If this holds we say
        that `(rightarrow, \twoheadrightarrow, \hookrightarrow)` is a factorization
        system.

    #. for all `x, y`, we do not have `x \twoheadrightarrow y \twoheadrightarrow x`
        or `x \hookrightarrow y \hookrightarrow x` without `x = y`. If this holds we
        say that the factorization system obeys the order condition.

    #. for all `x, y`, we do not have `x \twoheadrightarrow y \hookrightarrow x`
        without `x = y`. If this holds we say that the factorization system obeys the
        brick condition.

    Two-acyclic factorization systems play a key role in the Fundamental Theorem of
    Finite Semidistributive Lattices: as shown in [RST2024]_, checking whether a
    DiGraph forms a two-acyclic factorization system is equivalent to checking if
    the lattice of maximal orthogonal pairs is semidistributive.

    INPUT:

    - ``G`` -- DiGraph

    - ``certificate`` -- boolean (default: ``False``); whether to return a
        certificate

    OUTPUT:

    * If ``certificate=False``, return a boolean value.

    * If ``certificate=True``, and ``G`` is a two-acyclic factorization system,
      return a pair ``(True, None)``.

      * If ``G`` is not a two-acyclic factorization system, return a pair ``(False,
        (edge, error))``, where ``edge`` is an edge of the graph where some condition
        fails and ``error`` is one of the following strings:

        #. ``not_factorization_system`` -- ``G`` does not form a factorization system,
           ``edge`` cannot be factorized as a surjective edge followed by an injective

        #. ``surj_not_order`` -- the surjective edge ``edge`` fails the order condition

        #. ``inj_not_order`` -- the injective edge ``edge`` fails the order condition

        #. ``not_brick`` -- ``edge`` fails the brick condition

    EXAMPLES:

    Examples of graphs failing one or more conditions::

        sage: from sage.combinat.posets.semidistributivity import is_two_acyclic_factorization_system
        sage: G = DiGraph([(0, 1), (1, 0), (0, 0), (1, 1)], loops=True)
        sage: is_two_acyclic_factorization_system(G, certificate=True)            # (0, 1), (1, 0) are both surjective
        (False, ((0, 1), 'surj_not_order'))

        sage: G.add_edge(0, 2)                                                    # make (1, 0) no longer surjective
        sage: is_two_acyclic_factorization_system(G, certificate=True)            # (0, 1) surjective, (1, 0) injective
        (False, ((0, 1), 'not_brick'))

        sage: G.add_edge(1, 3)                                                    # make (0, 1) no longer surjective
        sage: is_two_acyclic_factorization_system(G, certificate=True)            # (0, 1), (1, 0) are both injective
        (False, ((0, 1), 'inj_not_order'))

        sage: G.add_edge(4, 1)                                                    # make (1, 0) no longer injective
        sage: is_two_acyclic_factorization_system(G, certificate=True)            # surjective edges are (0, 2), (1, 3)
        ....:                                                                     # injective edges are (0, 1), (4, 1)
        ....:                                                                     # so no factorization exists for (1, 0)
        (False, ((1, 0), 'not_factorization_system'))

    Example of a two-acyclic factorization system, from Figure 1. of [RST2024]_ ::

        sage: G = DiGraph([('b', 'a'), ('c', 'b'), ('c', 'a'), ('d', 'c'), ('e', 'd'), ('e', 'c'),
        ....:              ('e', 'b'), ('f', 'e'), ('f', 'c'), ('f', 'b'), ('f', 'a'), ('g', 'f'),
        ....:              ('g', 'e'), ('g', 'b'), ('g', 'a')])
        sage: is_two_acyclic_factorization_system(G)
        True

    Note that by the fundamental theorem of finite semidistributive lattices, its lattice of
    maximal orthogonal pairs must be semidistributive::

        sage: from sage.combinat.posets.semidistributivity import maximal_orthogonal_pairs_lattice
        sage: L = maximal_orthogonal_pairs_lattice(G)
        sage: L.is_semidistributive()
        True

    .. SEEALSO::

        - :meth:`surjective_edges()`
        - :meth:`injective_edges()`
        - :meth:`maximal_orthogonal_pairs_lattice()`

    REFERENCE:

    - [RST2024]_

    TESTS::

        sage: G = DiGraph()
        sage: is_two_acyclic_factorization_system(G)
        True
    """
    surjEdges = DiGraph(surjective_edges(G))
    injEdges = DiGraph(injective_edges(G))
    # the set of all factorizable edges
    Mult = set()
    for x, y in surjEdges.edge_iterator(labels=False):
        # since we demand the relations be reflexive
        # all surjective edges can be factorized as: xy => x ->> y -> y
        Mult.add((x, y))
        if surjEdges.has_edge(y, x):
            # since loops are not counted, we can assume x != y
            # so we have x ->> y ->> x with x != y, violating the order condition
            if certificate:
                return (False, ((x, y), "surj_not_order"))
            return False
        if y in injEdges:
            for z in injEdges.neighbor_out_iterator(y):
                if z == x:
                    # we have x ->> y -> x with x != y, violating the brick condition
                    if certificate:
                        return (False, ((x, y), "not_brick"))
                    return False
                # we have x ->> y -> z, so (x, z) can be factorized
                Mult.add((x, z))
    for x, y in injEdges.edge_iterator(labels=False):
        # all injective edges can be factorized as: xy => x ->> x -> y
        Mult.add((x, y))
        if injEdges.has_edge(y, x):
            # we have x -> y -> x with x != y, violating the order condition
            if certificate:
                return (False, ((x, y), "inj_not_order"))
            return False
    for x, y in G.edge_iterator(labels=False):
        if x != y and (x, y) not in Mult:
            # in a factorization system, all edges should be factorizable
            if certificate:
                return (False, ((x, y), "not_factorization_system"))
            return False
    # at this point all conditions should have been checked
    if certificate:
        return (True, None)
    return True
