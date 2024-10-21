# sage_setup: distribution = sagemath-graphs
# cython: binding=True
r"""
Cutwidth

This module implements several algorithms to compute the cutwidth of a graph and
the corresponding ordering of the vertices. It also implements tests functions
for evaluation the width of a linear ordering (or layout).

Given an ordering
`v_1,\cdots, v_n` of the vertices of `V(G)`, its *cost* is defined as:

.. MATH::

    c(v_1, ..., v_n) = \max_{1\leq i \leq n-1} c'(\{v_1, ..., v_i\})

Where

.. MATH::

    c'(S) = |\{(u,w)\in E(G)\mid u\in S\text{ and }w\in V(G)\backslash S\}|

The *cutwidth* of a graph `G` is equal to the minimum cost of an ordering of its
vertices.


**This module contains the following methods**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`cutwidth` | Return the cutwidth of the graph and the corresponding vertex ordering.
    :meth:`cutwidth_dyn` | Compute the cutwidth of `G` using an exponential time and space algorithm based on dynamic programming
    :meth:`cutwidth_MILP` | Compute the cutwidth of `G` and the optimal ordering of its vertices using an MILP formulation
    :meth:`width_of_cut_decomposition` | Return the width of the cut decomposition induced by the linear ordering `L` of the vertices of `G`


Exponential algorithm for cutwidth
----------------------------------

In order to find an optimal ordering of the vertices for the vertex separation,
this algorithm tries to save time by computing the function `c'(S)` **at most
once** once for each of the sets `S\subseteq V(G)`. These values are stored in
an array of size `2^n` where reading the value of `c'(S)` or updating it can be
done in constant time.

Assuming that we can compute the cost of a set `S` and remember it, finding an
optimal ordering is an easy task. Indeed, we can think of the sequence `v_1,
..., v_n` of vertices as a sequence of *sets* `\{v_1\}, \{v_1,v_2\}, ...,
\{v_1,...,v_n\}`, whose cost is precisely `\max c'(\{v_1\}), c'(\{v_1,v_2\}),
... , c'(\{v_1,...,v_n\})`. Hence, when considering the digraph on the `2^n`
sets `S\subseteq V(G)` where there is an arc from `S` to `S'` if `S'=S\cap
\{v\}` for some `v` (that is, if the sets `S` and `S'` can be consecutive in a
sequence), an ordering of the vertices of `G` corresponds to a *path* from
`\emptyset` to `\{v_1,...,v_n\}`. In this setting, checking whether there exists
a ordering of cost less than `k` can be achieved by checking whether there
exists a directed path `\emptyset` to `\{v_1,...,v_n\}` using only sets of cost
less than `k`. This is just a depth-first-search, for each `k`.

**Lazy evaluation of** `c'`

In the previous algorithm, most of the time is actually spent on the computation
of `c'(S)` for each set `S\subseteq V(G)` -- i.e. `2^n` computations of
neighborhoods. This can be seen as a huge waste of time when noticing that it is
useless to know that the value `c'(S)` for a set `S` is less than `k` if all the
paths leading to `S` have a cost greater than `k`. For this reason, the value of
`c'(S)` is computed lazily during the depth-first search. Explanation :

When the depth-first search discovers a set of size less than `k`, the costs of
its out-neighbors (the potential sets that could follow it in the optimal
ordering) are evaluated. When an out-neighbor is found that has a cost smaller
than `k`, the depth-first search continues with this set, which is explored with
the hope that it could lead to a path toward `\{v_1,...,v_n\}`. On the other
hand, if an out-neighbour has a cost larger than `k` it is useless to attempt to
build a cheap sequence going though this set, and the exploration stops
there. This way, a large number of sets will never be evaluated and *a lot* of
computational time is saved this way.

Besides, some improvement is also made by "improving" the values found by
`c'`. Indeed, `c'(S)` is a lower bound on the cost of a sequence containing the
set `S`, but if all out-neighbors of `S` have a cost of `c'(S) + 5` then one
knows that having `S` in a sequence means a total cost of at least `c'(S) +
5`. For this reason, for each set `S` we store the value of `c'(S)`, and replace
it by `\max (c'(S), \min_{\text{next}})` (where `\min_{\text{next}}` is the
minimum of the costs of the out-neighbors of `S`) once the costs of these
out-neighbors have been evaluated by the algorithm.

This algorithm and its implementation are very similar to
:meth:`sage.graphs.graph_decompositions.vertex_separation.vertex_separation_exp`.
The main difference is in the computation of `c'(S)`. See the :mod:`vertex
separation module's documentation
<sage.graphs.graph_decompositions.vertex_separation>` for more details on this
algorithm.

.. NOTE::

    Because of its current implementation, this algorithm only works on graphs
    on strictly less than 32 vertices. This can be changed to 64 if necessary,
    but 32 vertices already require 4GB of memory.


MILP formulation for the cutwidth
---------------------------------

We describe a mixed integer linear program (MILP) for determining an
optimal layout for the cutwidth of `G`.

**Variables:**

- `x_v^k` -- variable set to 1 if vertex `v` is placed in the ordering at
  position `i` with `i\leq k`, and 0 otherwise

- `y_{u,v}^{k}` -- variable set to 1 if one of `u` or `v` is at a position
  `i\leq k` and the other is at a position `j>k`, and so we have to count edge
  `uv` at position `k`. Otherwise, `y_{u,v}^{k}=0`. The value of `y_{u,v}^{k}`
  is a xor of the values of `x_u^k` and `x_v^k`.

- `z` -- objective value to minimize. It is equal to the maximum over all
  position `k` of the number of edges with one extremity at position at most `k`
  and the other at position strictly more than `k`, that is `\sum_{uv\in
  E}y_{u,v}^{k}`.


**MILP formulation:**

.. MATH::
    :nowrap:

    \begin{alignat*}{2}
    \intertext{Minimize:}
    &z&\\
    \intertext{Subject to:}
    \sum_{i=0}^{k-1}x_v^i &\leq k*x_v^{k} & \forall v\in V,\ k\in[1,n-1] \quad(1)\\
    x_v^n & =1  & \quad \forall v\in V \quad(2)\\
    \sum_{v\in V}x_v^k & = k+1 &\quad  \forall k\in[0,n-1] \quad(3)\\
    x_u^k - x_v^k & \leq y_{u,v}^k &\quad  \forall uv\in E,\ \forall k\in[0,n-1]  \quad(4)\\
    x_v^k - x_u^k & \leq y_{u,v}^k &\quad  \forall uv\in E,\ \forall k\in[0,n-1] \quad(5)\\
    \sum_{uv\in E}y_{u,v}^k &\leq z &\quad  \forall k\in[0,n-1] \quad(6)\\
    0 \leq z &\leq |E|
    \end{alignat*}

Constraints (1)-(3) ensure that all vertices have a distinct position.
Constraints (4)-(5) force variable `y_{u,v}^k` to 1 if the edge is in the cut.
Constraint (6) count the number of edges starting at position at most `k` and
ending at a position strictly larger than `k`.

This formulation corresponds to method :meth:`cutwidth_MILP`.


Authors
-------

- David Coudert (2015-06): Initial version


Methods
-------
"""

# ****************************************************************************
#       Copyright (C) 2015 David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from libc.stdint cimport uint8_t
from libc.string cimport memset
from cysignals.memory cimport check_allocarray, sig_free
from cysignals.signals cimport sig_check

from sage.graphs.graph_decompositions.fast_digraph cimport FastDigraph, popcount32
from sage.graphs.graph_decompositions.vertex_separation cimport find_order
from sage.graphs.graph_decompositions.vertex_separation import is_valid_ordering
from sage.rings.integer_ring import ZZ


################################################################################
# Measurement function of the width of some layout
################################################################################

def width_of_cut_decomposition(G, L):
    r"""
    Return the width of the cut decomposition induced by the linear ordering
    `L` of the vertices of `G`.

    If `G` is an instance of :mod:`Graph <sage.graphs.graph>`, this function
    returns the width `cw_L(G)` of the cut decomposition induced by the linear
    ordering `L` of the vertices of `G`.

    .. MATH::

        cw_L(G) =  \max_{0\leq i< |V|-1} |\{(u,w)\in E(G)\mid u\in L[:i]\text{ and }w\in V(G)\setminus L[:i]\}|

    INPUT:

    - ``G`` -- a Graph

    - ``L`` -- a linear ordering of the vertices of ``G``

    EXAMPLES:

    Cut decomposition of a Cycle graph::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: G = graphs.CycleGraph(6)
        sage: L = G.vertices(sort=False)
        sage: cutwidth.width_of_cut_decomposition(G, L)
        2

    Cut decomposition of a Path graph::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: P = graphs.PathGraph(6)
        sage: cutwidth.width_of_cut_decomposition(P, [0, 1, 2, 3, 4, 5])
        1
        sage: cutwidth.width_of_cut_decomposition(P, [5, 0, 1, 2, 3, 4])
        2
        sage: cutwidth.width_of_cut_decomposition(P, [0, 2, 4, 1, 3, 5])
        5

    TESTS:

    Giving a wrong linear ordering::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: cutwidth.width_of_cut_decomposition(Graph(), ['a','b'])
        Traceback (most recent call last):
        ...
        ValueError: the input linear vertex ordering L is not valid for G
    """
    if not is_valid_ordering(G, L):
        raise ValueError("the input linear vertex ordering L is not valid for G")
    elif G.order() <= 1:
        return 0

    cdef int i, x, y
    cdef dict position = {u: i for i, u in enumerate(L)}

    # We count for each position `i` the number of edges going from vertices at
    # positions in `0..i` to vertices at positions in `i+1..n-1`, for each
    # `x\leq i<n-1`.
    cdef list cpt = [0] * G.order()
    for u, v in G.edge_iterator(labels=None):
        x, y = position[u], position[v]
        if x > y:
            x, y = y, x
        # Edge (u,v) contributes 1 to the number of edges going from vertices at
        # positions `0..i` to vertices at positions `i+1..n-1` for each `x\leq
        # i < n-1`.
        for i in range(x, y):
            cpt[i] += 1

    # The width of L is the maximum computed value.
    return max(cpt)


################################################################################
# Front end method for cutwidth
################################################################################

def cutwidth(G, algorithm='exponential', cut_off=0, solver=None, verbose=False,
             *, integrality_tolerance=1e-3):
    r"""
    Return the cutwidth of the graph and the corresponding vertex ordering.

    INPUT:

    - ``G`` -- a Graph or a DiGraph

    - ``algorithm`` -- string (default: ``'exponential'``); algorithm to use
      among:

      - ``exponential`` -- use an exponential time and space algorithm based on
        dynamic programming. This algorithm only works on graphs with strictly
        less than 32 vertices.

      - ``MILP`` -- use a mixed integer linear programming formulation. This
        algorithm has no size restriction but could take a very long time

    - ``cut_off`` -- integer (default: 0); used to stop the search as soon as a
      solution with width at most ``cut_off`` is found, if any. If this bound
      cannot be reached, the best solution found is returned.

    - ``solver`` -- string (default: ``None``); specifies a Mixed Integer Linear
      Programming (MILP) solver to be used. If set to ``None``, the default one
      is used. For more information on MILP solvers and which default solver is
      used, see the method :meth:`solve
      <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
      :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- boolean (default: ``False``); whether to display
      information on the computations

    - ``integrality_tolerance`` -- float; parameter for use with MILP solvers
      over an inexact base ring; see
      :meth:`MixedIntegerLinearProgram.get_values`.

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.

    EXAMPLES:

    Cutwidth of a Complete Graph::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth
        sage: G = graphs.CompleteGraph(5)
        sage: cw,L = cutwidth(G); cw
        6
        sage: K = graphs.CompleteGraph(6)
        sage: cw,L = cutwidth(K); cw
        9
        sage: cw,L = cutwidth(K+K); cw
        9

    The cutwidth of a `p\times q` Grid Graph with `p\leq q` is `p+1`::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth
        sage: G = graphs.Grid2dGraph(3,3)
        sage: cw,L = cutwidth(G); cw
        4
        sage: G = graphs.Grid2dGraph(3,5)
        sage: cw,L = cutwidth(G); cw
        4

    TESTS:

    Comparison of algorithms::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth
        sage: for i in range(2):                # long time                             # needs sage.numerical.mip
        ....:     G = graphs.RandomGNP(7, 0.3)
        ....:     ve, le = cutwidth(G, algorithm='exponential')
        ....:     vm, lm = cutwidth(G, algorithm='MILP', solver='GLPK')
        ....:     if ve != vm:
        ....:        raise ValueError("Something goes wrong!")

    Given a wrong algorithm::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth
        sage: cutwidth(graphs.PathGraph(2), algorithm='SuperFast')
        Traceback (most recent call last):
        ...
        ValueError: algorithm "SuperFast" has not been implemented yet, please contribute

    Given anything else than a Graph::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth
        sage: cutwidth(range(4))
        Traceback (most recent call last):
        ...
        ValueError: the first parameter must be a Graph

    Giving a wrong type cut off::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth
        sage: cutwidth(Graph(), cut_off='toto')
        Traceback (most recent call last):
        ...
        ValueError: the specified cut off parameter must be an integer

    Cutwidth of a graph with one edge (:issue:`32131`)::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth
        sage: G = Graph([(0, 1)])
        sage: cutwidth(G, algorithm='exponential')
        (1, [0, 1])
        sage: cutwidth(G, algorithm='MILP', solver='GLPK')                              # needs sage.numerical.mip
        (1, [0, 1])

    Cutwidth of a disconnected graph::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth
        sage: G = Graph(5)
        sage: G.add_edge(2, 3)
        sage: cutwidth(G, algorithm='exponential')
        (1, [2, 3, 0, 1, 4])
        sage: cutwidth(G, algorithm='MILP', solver='GLPK')                              # needs sage.numerical.mip
        (1, [2, 3, 0, 1, 4])
    """
    from sage.graphs.graph import Graph

    if not isinstance(G, Graph):
        raise ValueError('the first parameter must be a Graph')

    if cut_off not in ZZ:
        raise ValueError("the specified cut off parameter must be an integer")
    elif G.size() <= cut_off:
        # We have a trivial solution
        return width_of_cut_decomposition(G, list(G)), list(G)

    cdef list CC
    if not G.is_connected():
        CC = G.connected_components_subgraphs()
    else:
        CC = [G]

    # If the graph has several connected components. We solve the problem on
    # each of them and concatenate the partial orderings. The cutwidth is the
    # maximum over all these subgraphs.
    cdef int cw = 0
    cdef list L = []
    cdef int cwH
    cdef list LH
    cdef int this_cut_off = cut_off

    for H in CC:

        if H.size() <= this_cut_off:
            # We can directly add the vertices to the solution
            L.extend(H)

        else:
            # We have a connected graph and we call the desired algorithm
            if algorithm == "exponential":
                cwH, LH = cutwidth_dyn(H, lower_bound=this_cut_off)

            elif algorithm == "MILP":
                cwH, LH = cutwidth_MILP(H, lower_bound=this_cut_off, solver=solver,
                                        verbose=verbose, integrality_tolerance=integrality_tolerance)

            else:
                raise ValueError('algorithm "{}" has not been implemented yet, please contribute'.format(algorithm))

            # We update the cutwidth and ordering
            cw = max(cw, cwH)
            L.extend(LH)
            this_cut_off = max(cw, this_cut_off)

    return cw, L


################################################################################
# Dynamic Programming algorithm for cutwidth
################################################################################

def cutwidth_dyn(G, lower_bound=0):
    r"""
    Dynamic programming algorithm for the cutwidth of a Graph.

    This function uses dynamic programming algorithm for determining an optimal
    layout for the cutwidth of `G`. See the :mod:`module's documentation
    <sage.graphs.graph_decompositions.cutwidth>` for more details on this
    method.

    INPUT:

    - ``G`` -- a Graph

    - ``lower_bound`` -- integer (default: 0); the algorithm returns immediately
      if it finds a solution lower or equal to ``lower_bound`` (in which case it
      may not be optimal).

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.

    .. NOTE::

        Because of its current implementation, this algorithm only works on
        graphs on strictly less than 32 vertices. This can be changed to 63 if
        necessary, but 32 vertices already require 4GB of memory.

    TESTS:

    Giving anything else than a Graph::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: cutwidth.cutwidth_dyn([])
        Traceback (most recent call last):
        ...
        ValueError: the parameter must be a Graph

    Giving a too large Graph::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: cutwidth.cutwidth_dyn(graphs.PathGraph(40))
        Traceback (most recent call last):
        ...
        ValueError: the graph should have at most 31 vertices

    Giving a wrong type lower bound::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: cutwidth.cutwidth_dyn(Graph(), lower_bound='toto')
        Traceback (most recent call last):
        ...
        ValueError: the specified lower bound must be an integer
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("the parameter must be a Graph")

    if G.order() >= 32:
        raise ValueError("the graph should have at most 31 vertices")

    if lower_bound not in ZZ:
        raise ValueError("the specified lower bound must be an integer")

    cdef FastDigraph g = FastDigraph(G)

    cdef unsigned int mem = 1 << g.n
    cdef uint8_t* neighborhoods = <uint8_t*> check_allocarray(mem, sizeof(uint8_t))

    memset(neighborhoods, <uint8_t> -1, mem)

    cdef int i, k
    cdef list order

    try:
        for k in range(lower_bound, G.size() + 1):
            for i in range(g.n):
                sig_check()
                if exists(g, neighborhoods, 0, 0, i, k) <= k:
                    order = find_order(g, neighborhoods, k)
                    return k, [g.int_to_vertices[i] for i in order]

        order = find_order(g, neighborhoods, k)
        return k, [g.int_to_vertices[i] for i in order]
    finally:
        sig_free(neighborhoods)


cdef inline int exists(FastDigraph g, uint8_t* neighborhoods, int S, int cost_S, int v, int k) noexcept:
    r"""
    Check whether an ordering with the given cost `k` exists, and updates data
    in the neighborhoods array at the same time. See the module's documentation.

    INPUT:

    - ``g`` -- a FastDiGraph

    - ``neighborhoods`` -- an array of size `2^(g.n)`; stores for each subset
      `X\subseteq V` of vertices of the graph the number of edges from `X` to
      `V\setminus X`

    - ``S`` -- integer; encodes the predecessor subset of vertices (from which
      is issued the current call)

    - ``cost_S`` -- integer; the number of edges from `S` to `V\setminus S`

    - ``v`` -- integer; a vertex such that the current subset of vertices is
      `current==S\cup\{v\}`

    - ``k`` -- integer; the maximum admissible cost for a solution
    """
    cdef int current = S | 1 << v
    # If this is true, it means the set has not been evaluated yet
    if neighborhoods[current] == <uint8_t> -1:
        # The number of edges from `current` to `V\setminus current` is the
        # number of edges from `S` to `V\setminus S`, minus the number of edges
        # from `S` to vertex `v`, plus the number of edges from `v` to
        # `V\setminus (S\cup \{v\})`. This can be computed adding the degree of
        # `c` to `cost_S`, and then removing twice the number of edges from `S`
        # to `v`.
        neighborhoods[current] = cost_S + g.degree[v] - 2 * popcount32(S & g.graph[v])

    # If the cost of this set is too high, there is no point in going further.
    # Same thing if the current set is the whole vertex set.
    if neighborhoods[current] > k or (current == (1 << g.n) - 1):
        return neighborhoods[current]

    # Minimum of the costs of the outneighbors, initialized with large constant.
    cdef int mini = (<uint8_t> -1)

    cdef int i

    # For each possible extension of the current set with a vertex, check whether
    # there exists a cheap path toward {1..n}, and update the cost.
    for i in range(g.n):
        if (current >> i) & 1:  # if i in S
            continue

        mini = min(mini, exists(g, neighborhoods, current, neighborhoods[current], i, k))

        # We have found a path !
        if mini <= k:
            return mini

    # Updating the cost of the current set with the minimum of the cost of its
    # outneighbors.
    neighborhoods[current] = mini

    return neighborhoods[current]


################################################################################
# MILP formulations for cutwidth
################################################################################

def cutwidth_MILP(G, lower_bound=0, solver=None, verbose=0,
                  *, integrality_tolerance=1e-3):
    r"""
    MILP formulation for the cutwidth of a Graph.

    This method uses a mixed integer linear program (MILP) for determining an
    optimal layout for the cutwidth of `G`. See the :mod:`module's documentation
    <sage.graphs.graph_decompositions.cutwidth>` for more details on this MILP
    formulation.

    INPUT:

    - ``G`` -- a Graph

    - ``lower_bound`` -- integer (default: 0); the algorithm searches for a
      solution with cost larger or equal to ``lower_bound``. If the given bound
      is larger than the optimal solution the returned solution might not be
      optimal. If the given bound is too high, the algorithm might not be able
      to find a feasible solution.

    - ``solver`` -- string (default: ``None``); specifies a Mixed Integer Linear
      Programming (MILP) solver to be used. If set to ``None``, the default one
      is used. For more information on MILP solvers and which default solver is
      used, see the method :meth:`solve
      <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
      :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- integer (default: 0); sets the level of verbosity. Set
      to 0 by default, which means quiet.

    - ``integrality_tolerance`` -- float; parameter for use with MILP solvers
      over an inexact base ring; see
      :meth:`MixedIntegerLinearProgram.get_values`.

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.

    EXAMPLES:

    Cutwidth of a Cycle graph::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: G = graphs.CycleGraph(5)
        sage: cw, L = cutwidth.cutwidth_MILP(G); cw                                     # needs sage.numerical.mip
        2
        sage: cw == cutwidth.width_of_cut_decomposition(G, L)                           # needs sage.numerical.mip
        True
        sage: cwe, Le = cutwidth.cutwidth_dyn(G); cwe
        2

    Cutwidth of a Complete graph::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: G = graphs.CompleteGraph(4)
        sage: cw, L = cutwidth.cutwidth_MILP(G); cw                                     # needs sage.numerical.mip
        4
        sage: cw == cutwidth.width_of_cut_decomposition(G, L)                           # needs sage.numerical.mip
        True

    Cutwidth of a Path graph::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: G = graphs.PathGraph(3)
        sage: cw, L = cutwidth.cutwidth_MILP(G); cw                                     # needs sage.numerical.mip
        1
        sage: cw == cutwidth.width_of_cut_decomposition(G, L)                           # needs sage.numerical.mip
        True

    TESTS:

    Comparison with exponential algorithm::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: for i in range(2):                # long time                             # needs sage.numerical.mip
        ....:     G = graphs.RandomGNP(7, 0.3)
        ....:     ve, le = cutwidth.cutwidth_dyn(G)
        ....:     vm, lm = cutwidth.cutwidth_MILP(G, solver='GLPK')
        ....:     if ve != vm:
        ....:        print("The solution is not optimal!")

    Giving a too large lower bound::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth_MILP
        sage: G = graphs.CycleGraph(3)
        sage: cutwidth_MILP(G, lower_bound=G.size()+1)                                  # needs sage.numerical.mip
        Traceback (most recent call last):
        ...
        MIPSolverException: ...

    Giving anything else than a Graph::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth_MILP
        sage: cutwidth_MILP([])
        Traceback (most recent call last):
        ...
        ValueError: the first input parameter must be a Graph
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("the first input parameter must be a Graph")

    from sage.numerical.mip import MixedIntegerLinearProgram
    p = MixedIntegerLinearProgram(maximization=False, solver=solver)

    # Declaration of variables.
    x = p.new_variable(binary=True, nonnegative=True)
    y = p.new_variable(binary=True, nonnegative=True)
    z = p.new_variable(integer=True, nonnegative=True)

    N = G.order()

    # All vertices at different positions
    for v in G:
        for k in range(N - 1):
            p.add_constraint(p.sum(x[v, i] for i in range(k)) <= k * x[v, k])
        p.add_constraint(x[v, N - 1] == 1)
    for k in range(N):
        p.add_constraint(p.sum(x[v, k] for v in G) == k + 1)

    # Edge uv counts at position i if one of u or v is placed at a position in
    # [0,i] and the other is placed at a position in [i+1,n].
    for u, v in G.edge_iterator(labels=None):
        for i in range(N):
            p.add_constraint(x[u, i] - x[v, i] <= y[u, v, i])
            p.add_constraint(x[v, i] - x[u, i] <= y[u, v, i])

    # Lower bound on the solution
    p.add_constraint(lower_bound <= z['z'])

    # Objective
    p.add_constraint(z['z'] <= G.size())
    for i in range(N):
        p.add_constraint(p.sum(y[u, v, i] for u, v in G.edge_iterator(labels=None)) <= z['z'])

    p.set_objective(z['z'])

    _ = p.solve(log=verbose)

    # We now extract the ordering and the cost of the solution
    val_x = p.get_values(x, convert=bool, tolerance=integrality_tolerance)
    cdef int cw = p.get_values(z['z'], convert=True, tolerance=integrality_tolerance)
    cdef list seq = []
    cdef set to_see = set(G)
    for k in range(N):
        for u in to_see:
            if val_x[u, k]:
                seq.append(u)
                to_see.discard(u)
                break

    return cw, seq
