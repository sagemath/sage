# cython: binding=True
r"""
Vertex separation

This module implements several algorithms to compute the vertex separation of a
digraph and the corresponding ordering of the vertices. It also implements tests
functions for evaluation the width of a linear ordering.

Given an ordering
`v_1,\cdots, v_n` of the vertices of `V(G)`, its *cost* is defined as:

.. MATH::

    c(v_1, ..., v_n) = \max_{1\leq i \leq n} c'(\{v_1, ..., v_i\})

Where

.. MATH::

    c'(S) = |N^+_G(S)\backslash S|

The *vertex separation* of a digraph `G` is equal to the minimum cost of an
ordering of its vertices.

**Vertex separation and pathwidth**

The vertex separation is defined on a digraph, but one can obtain from a graph
`G` a digraph `D` with the same vertex set, and in which each edge `uv` of `G`
is replaced by two edges `uv` and `vu` in `D`. The vertex separation of `D` is
equal to the pathwidth of `G`, and the corresponding ordering of the vertices of
`D`, also called a *layout*, encodes an optimal path-decomposition of `G`.
This is a result of Kinnersley [Kin1992]_ and Bodlaender [Bod1998]_.


**This module contains the following methods**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`pathwidth` | Compute the pathwidth of ``self`` (and provides a decomposition)
    :meth:`path_decomposition` | Return the pathwidth of the given graph and the ordering of the vertices resulting in a corresponding path decomposition
    :meth:`vertex_separation` | Return an optimal ordering of the vertices and its cost for vertex-separation
    :meth:`vertex_separation_exp` | Compute the vertex separation of `G` using an exponential time and space algorithm
    :meth:`vertex_separation_MILP` | Compute the vertex separation of `G` and the optimal ordering of its vertices using an MILP formulation
    :meth:`vertex_separation_BAB` | Compute the vertex separation of `G` and the optimal ordering of its vertices using a branch and bound algorithm
    :meth:`lower_bound` | Return a lower bound on the vertex separation of `G`
    :meth:`is_valid_ordering` | Test if the linear vertex ordering `L` is valid for (di)graph `G`
    :meth:`width_of_path_decomposition` | Return the width of the path decomposition induced by the linear ordering `L` of the vertices of `G`
    :meth:`linear_ordering_to_path_decomposition`| Return the path decomposition encoded in the ordering `L`


Exponential algorithm for vertex separation
-------------------------------------------

In order to find an optimal ordering of the vertices for the vertex separation,
this algorithm tries to save time by computing the function `c'(S)` **at most
once** once for each of the sets `S\subseteq V(G)`. These values are stored in
an array of size `2^n` where reading the value of `c'(S)` or updating it can be
done in constant (and small) time.

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

.. NOTE::

    Because of its current implementation, this algorithm only works on graphs
    on less than 32 vertices. This can be changed to 64 if necessary, but 32
    vertices already require 4GB of memory. Running it on 64 bits is not
    expected to be doable by the computers of the next decade `:-D`

**Lower bound on the vertex separation**

One can obtain a lower bound on the vertex separation of a graph in exponential
time but *small* memory by computing once the cost of each set `S`. Indeed, the
cost of a sequence `v_1, ..., v_n` corresponding to sets `\{v_1\}, \{v_1,v_2\},
..., \{v_1,...,v_n\}` is

.. MATH::

    \max c'(\{v_1\}),c'(\{v_1,v_2\}),...,c'(\{v_1,...,v_n\})\geq\max c'_1,...,c'_n

where `c_i` is the minimum cost of a set `S` on `i` vertices. Evaluating the
`c_i` can take time (and in particular more than the previous exact algorithm),
but it does not need much memory to run.


MILP formulation for the vertex separation
------------------------------------------

We describe below a mixed integer linear program (MILP) for determining an
optimal layout for the vertex separation of `G`, which is an improved version of
the formulation proposed in [SP2010]_. It aims at building a sequence `S_t` of
sets such that an ordering `v_1, ..., v_n` of the vertices correspond to
`S_0=\{v_1\}, S_2=\{v_1,v_2\}, ..., S_{n-1}=\{v_1,...,v_n\}`.

**Variables:**


- `y_v^t` -- Variable set to 1 if `v\in S_t`, and 0 otherwise. The order of
  `v` in the layout is the smallest `t` such that `y_v^t==1`.

- `u_v^t` -- Variable set to 1 if `v\not \in S_t` and `v` has an in-neighbor in
  `S_t`. It is set to 0 otherwise.

- `x_v^t` -- Variable set to 1 if either `v\in S_t` or if `v` has an in-neighbor
  in `S_t`. It is set to 0 otherwise.

- `z` -- Objective value to minimize. It is equal to the maximum over all step
  `t` of the number of vertices such that `u_v^t==1`.

**MILP formulation:**

.. MATH::
    :nowrap:

    \begin{alignat}{2}
    \text{Minimize:}
    &z&\\
    \text{Such that:}
    x_v^t &\leq x_v^{t+1}& \forall v\in V,\ 0\leq t\leq n-2\\
    y_v^t &\leq y_v^{t+1}& \forall v\in V,\ 0\leq t\leq n-2\\
    y_v^t &\leq x_w^t& \forall v\in V,\ \forall w\in N^+(v),\ 0\leq t\leq n-1\\
    \sum_{v \in V} y_v^{t} &= t+1& 0\leq t\leq n-1\\
    x_v^t-y_v^t&\leq u_v^t & \forall v \in V,\ 0\leq t\leq n-1\\
    \sum_{v \in V} u_v^t &\leq z& 0\leq t\leq n-1\\
    0 \leq x_v^t &\leq 1& \forall v\in V,\ 0\leq t\leq n-1\\
    0 \leq u_v^t &\leq 1& \forall v\in V,\ 0\leq t\leq n-1\\
    y_v^t &\in \{0,1\}& \forall v\in V,\ 0\leq t\leq n-1\\
    0 \leq z &\leq n&
    \end{alignat}

The vertex separation of `G` is given by the value of `z`, and the order of
vertex `v` in the optimal layout is given by the smallest `t` for which
`y_v^t==1`.


Branch and Bound algorithm for the vertex separation
----------------------------------------------------

We describe below the principle of a branch and bound algorithm (BAB) for
determining an optimal ordering for the vertex separation of `G`, as proposed in
[CMN2014]_.

**Greedy steps:**

Let us denote `{\cal L}(S)` the set of all possible orderings of the vertices in
`S`, and let `{\cal L}_P(S)\subseteq {\cal L}(S)` be the orderings starting with
a prefix `P`. Let also `c(L)` be the cost of the ordering `L\in{\cal L}(V)` as
defined above.

Given a digraph `D=(V,A)`, a set `S\subset V`, and a prefix `P`, it has been
proved in [CMN2014]_ that `\min_{L\in{\cal L}_P(V)} c(L) = \min_{L\in{\cal
L}_{P+v}(V)} c(L)` holds in two (non exhaustive) cases:

.. MATH::

    \text{or} \begin{cases}
    N^+(v)\subseteq S\cup N^+(S)\\
    v\in N^+(S)\text{ and }N^+(v)\setminus(S\cup N^+(S)) = \{w\}
    \end{cases}

In other words, if we find a vertex `v` satisfying the above conditions, the
best possible ordering with prefix `P` has the same cost as the best possible
ordering with prefix `P+v`. So we can greedily extend the prefix with vertices
satisfying the conditions which results in a significant reduction of the search
space.


**The algorithm:**

Given the current prefix `P` and the current upper bound `UB` (either an input
upper bound or the cost of the best solution found so far), apply the following
steps:

- Extend the prefix `P` into a prefix `P'` using the greedy steps as described
  above.

- Sort the vertices `v\in V\setminus P'` by increasing values of `|N^+(P+v)|`,
  and prune the vertices with a value larger or equal to `UB`. Let `\Delta` be
  the resulting sorted list.

- Repeat with prefix `P'+v` for all `v\in\Delta` and keep the best found
  solution.

If a lower bound is passed to the algorithm, it will stop as soon as a solution
with cost equal to that lower bound is found.


**Storing prefixes:**

If for a prefix `P` we have `c(P)<\min_{L\in{\cal L}_P(V)} c(L)=C`, then for any
permutation `P'` of `P` we have `\min_{L\in{\cal L}_{P'}(V)} c(L)\geq C`.

Thus, given such a prefix `P` there is no need to explore any of the orderings
starting with one of its permutations. To do so, we store `P` (as a set of
vertices) to cut branches later. See [CMN2014]_ for more details.

Since the number of stored sets can get very large, one can control the maximum
length and the maximum number of stored prefixes.


Authors
-------

- Nathann Cohen (2011-10): Initial version and exact exponential algorithm

- David Coudert (2012-04): MILP formulation and tests functions

- David Coudert (2015-01): BAB formulation and tests functions


Methods
-------
"""

# ****************************************************************************
#       Copyright (C) 2011 Nathann Cohen <nathann.cohen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from libc.string cimport memset
from cysignals.memory cimport check_malloc, sig_malloc, sig_free
from cysignals.signals cimport sig_check, sig_on, sig_off

from sage.graphs.graph_decompositions.fast_digraph cimport FastDigraph, compute_out_neighborhood_cardinality, popcount32
from libc.stdint cimport uint8_t
from sage.data_structures.binary_matrix cimport *
from sage.graphs.base.static_dense_graph cimport dense_graph_init


###############
# Lower Bound #
###############

def lower_bound(G):
    r"""
    Return a lower bound on the vertex separation of `G`.

    INPUT:

    - ``G`` -- a Graph or a DiGraph

    OUTPUT:

    A lower bound on the vertex separation of `D` (see the module's
    documentation).

    .. NOTE::

        This method runs in exponential time but has no memory constraint.


    EXAMPLES:

    On a circuit::

        sage: from sage.graphs.graph_decompositions.vertex_separation import lower_bound
        sage: g = digraphs.Circuit(6)
        sage: lower_bound(g)
        1

    TESTS:

    Given anything else than a Graph or a DiGraph::

        sage: from sage.graphs.graph_decompositions.vertex_separation import lower_bound
        sage: lower_bound(range(2))
        Traceback (most recent call last):
        ...
        ValueError: the parameter must be a Graph or a DiGraph

    Given a too large graph::

        sage: from sage.graphs.graph_decompositions.vertex_separation import lower_bound
        sage: lower_bound(graphs.PathGraph(50))
        Traceback (most recent call last):
        ...
        ValueError: the (di)graph can have at most 31 vertices

    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, Graph) and not isinstance(G, DiGraph):
        raise ValueError("the parameter must be a Graph or a DiGraph")

    if G.order() >= 32:
        raise ValueError("the (di)graph can have at most 31 vertices")

    cdef FastDigraph FD = FastDigraph(G)
    cdef int * g = FD.graph
    cdef unsigned int n = <unsigned int>FD.n

    # minimums[i] is means to store the value of c'_{i+1}
    cdef uint8_t* minimums = <uint8_t*> check_malloc(n * sizeof(uint8_t))
    cdef unsigned int i

    # They are initialized to n
    for i in range(n):
        minimums[i] = n

    cdef uint8_t tmp, tmp_count

    # We go through all sets
    for i in range(1, <unsigned int> (1 << n)):
        tmp_count = <uint8_t> popcount32(i)
        tmp = <uint8_t> compute_out_neighborhood_cardinality(FD, i)

        # And update the costs
        minimums[tmp_count-1] = minimum(minimums[tmp_count-1], tmp)

    # We compute the maximum of all those values
    for i in range(1, n):
        minimums[0] = maximum(minimums[0], minimums[i])

    cdef int lb = minimums[0]

    sig_free(minimums)

    return lb


###################################################################
# Method for turning an ordering to a path decomposition and back #
###################################################################

def linear_ordering_to_path_decomposition(G, L):
    """
    Return the path decomposition encoded in the ordering L

    INPUT:

    - ``G`` -- a Graph

    - ``L`` -- a linear ordering for G

    OUTPUT:

    A path graph whose vertices are the bags of the path decomposition.

    EXAMPLES:

    The bags of an optimal path decomposition of a path-graph have two vertices
    each::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation
        sage: from sage.graphs.graph_decompositions.vertex_separation import linear_ordering_to_path_decomposition
        sage: g = graphs.PathGraph(5)
        sage: pw, L = vertex_separation(g, algorithm = "BAB"); pw
        1
        sage: h = linear_ordering_to_path_decomposition(g, L)
        sage: sorted(h, key=str)
        [{0, 1}, {1, 2}, {2, 3}, {3, 4}]
        sage: sorted(h.edge_iterator(labels=None), key=str)
        [({0, 1}, {1, 2}), ({1, 2}, {2, 3}), ({2, 3}, {3, 4})]

    Giving a non-optimal linear ordering::

        sage: g = graphs.PathGraph(5)
        sage: L = [1, 4, 0, 2, 3]
        sage: from sage.graphs.graph_decompositions.vertex_separation import width_of_path_decomposition
        sage: width_of_path_decomposition(g, L)
        3
        sage: h = linear_ordering_to_path_decomposition(g, L)
        sage: h.vertices(sort=True)
        [{0, 2, 3, 4}, {0, 1, 2}]

    The bags of the path decomposition of a cycle have three vertices each::

        sage: g = graphs.CycleGraph(6)
        sage: pw, L = vertex_separation(g, algorithm = "BAB"); pw
        2
        sage: h = linear_ordering_to_path_decomposition(g, L)
        sage: sorted(h, key=str)
        [{0, 1, 5}, {1, 2, 5}, {2, 3, 4}, {2, 4, 5}]
        sage: sorted(h.edge_iterator(labels=None), key=str)
        [({0, 1, 5}, {1, 2, 5}), ({1, 2, 5}, {2, 4, 5}), ({2, 4, 5}, {2, 3, 4})]


    TESTS::

        sage: linear_ordering_to_path_decomposition(Graph(), [])
        Graph on 0 vertices
        sage: linear_ordering_to_path_decomposition(DiGraph(), [])
        Traceback (most recent call last):
        ...
        ValueError: the first parameter must be a Graph
        sage: g = graphs.CycleGraph(6)
        sage: linear_ordering_to_path_decomposition(g, list(range(7)))
        Traceback (most recent call last):
        ...
        ValueError: the input linear vertex ordering L is not valid for G
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("the first parameter must be a Graph")
    if not G:
        return Graph()
    if not is_valid_ordering(G, L):
        raise ValueError("the input linear vertex ordering L is not valid for G")

    cdef set seen = set()     # already treated vertices
    cdef set covered = set()  # vertices in the neighborhood of seen but not in seen
    cdef list bags = list()   # The bags of the path decomposition

    # We build the bags of the path-decomposition, and avoid adding useless bags
    for u in L:
        seen.add(u)
        covered.update(G.neighbor_iterator(u))
        covered.difference_update(seen)
        new_bag = covered.union([u])
        if bags:
            if new_bag.issubset(bags[-1]):
                continue
            if new_bag.issuperset(bags[-1]):
                bags.pop()

        bags.append(new_bag)

    # We now build a graph whose vertices are bags
    from sage.sets.set import Set
    H = Graph()
    H.add_path([Set(bag) for bag in bags])
    return H


##################################################################
# Front end methods for path decomposition and vertex separation #
##################################################################

def pathwidth(self, k=None, certificate=False, algorithm="BAB", verbose=False,
              max_prefix_length=20, max_prefix_number=10**6, *, solver=None):
    r"""
    Compute the pathwidth of ``self`` (and provides a decomposition)

    INPUT:

    - ``k`` -- integer (default: ``None``); the width to be considered. When
      ``k`` is an integer, the method checks that the graph has pathwidth
      `\leq k`. If ``k`` is ``None`` (default), the method computes the optimal
      pathwidth.

    - ``certificate`` -- boolean (default: ``False``); whether to return the
      path-decomposition itself

    - ``algorithm`` -- string (default: ``"BAB"``); algorithm to use among:

      - ``"BAB"`` -- Use a branch-and-bound algorithm. This algorithm has no
        size restriction but could take a very long time on large graphs. It can
        also be used to test is the input graph has pathwidth `\leq k`, in which
        cas it will return the first found solution with width `\leq k` is
        ``certificate==True``.

      - ``exponential`` -- Use an exponential time and space algorithm. This
        algorithm only works of graphs on less than 32 vertices.

      - ``MILP`` -- Use a mixed integer linear programming formulation. This
        algorithm has no size restriction but could take a very long time.

    - ``verbose`` -- boolean (default: ``False``); whether to display
      information on the computations

    - ``max_prefix_length`` -- integer (default: 20); limits the length of the
      stored prefixes to prevent storing too many prefixes. This parameter is
      used only when ``algorithm=="BAB"``.

    - ``max_prefix_number`` -- integer (default: 10**6); upper bound on the
      number of stored prefixes used to prevent using too much memory. This
      parameter is used only when ``algorithm=="BAB"``.

    - ``solver`` -- string (default: ``None``); specify a Mixed Integer Linear
      Programming (MILP) solver to be used. If set to ``None``, the default one
      is used. For more information on MILP solvers and which default solver is
      used, see the method :meth:`solve
      <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
      :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    OUTPUT:

    Return the pathwidth of ``self``. When ``k`` is specified, it returns
    ``False`` when no path-decomposition of width `\leq k` exists or ``True``
    otherwise. When ``certificate=True``, the path-decomposition is also
    returned.

    .. SEEALSO::

        * :meth:`Graph.treewidth` -- computes the treewidth of a graph
        * :meth:`~sage.graphs.graph_decompositions.vertex_separation.vertex_separation`
          -- computes the vertex separation of a (di)graph

    EXAMPLES:

    The pathwidth of a cycle is equal to 2::

        sage: g = graphs.CycleGraph(6)
        sage: g.pathwidth()
        2
        sage: pw, decomp = g.pathwidth(certificate=True)
        sage: sorted(decomp, key=str)
        [{0, 1, 5}, {1, 2, 5}, {2, 3, 4}, {2, 4, 5}]

    The pathwidth of a Petersen graph is 5::

        sage: g = graphs.PetersenGraph()
        sage: g.pathwidth()
        5
        sage: g.pathwidth(k=2)
        False
        sage: g.pathwidth(k=6)
        True
        sage: g.pathwidth(k=6, certificate=True)
        (True, Graph on 5 vertices)

    TESTS:

    Given anything else than a Graph::

        sage: from sage.graphs.graph_decompositions.vertex_separation import pathwidth
        sage: pathwidth(DiGraph())
        Traceback (most recent call last):
        ...
        ValueError: the parameter must be a Graph

    Given a wrong algorithm::

        sage: from sage.graphs.graph_decompositions.vertex_separation import pathwidth
        sage: pathwidth(Graph(), algorithm="SuperFast")
        Traceback (most recent call last):
        ...
        ValueError: algorithm "SuperFast" has not been implemented yet, please contribute

    Using a specific solver::

        sage: g = graphs.PetersenGraph()
        sage: g.pathwidth(solver='SCIP')  # optional - pyscipopt
        5
    """
    from sage.graphs.graph import Graph
    if not isinstance(self, Graph):
        raise ValueError("the parameter must be a Graph")

    pw, L = vertex_separation(self, algorithm=algorithm, verbose=verbose,
                              cut_off=k, upper_bound=None if k is None else (k+1),
                              max_prefix_length=max_prefix_length,
                              max_prefix_number=max_prefix_number,
                              solver=solver)

    if k is None:
        return (pw, linear_ordering_to_path_decomposition(self, L)) if certificate else pw
    if pw < 0:
        # no solution found
        return (False, Graph()) if certificate else False
    return (pw <= k, linear_ordering_to_path_decomposition(self, L)) if certificate else pw <= k


def path_decomposition(G, algorithm="BAB", cut_off=None, upper_bound=None, verbose=False,
                       max_prefix_length=20, max_prefix_number=10**6):
    r"""
    Return the pathwidth of the given graph and the ordering of the vertices
    resulting in a corresponding path decomposition.

    INPUT:

    - ``G`` -- a Graph

    - ``algorithm`` -- string (default: ``"BAB"``); algorithm to use among:

      - ``"BAB"`` -- Use a branch-and-bound algorithm. This algorithm has no
        size restriction but could take a very long time on large graphs. It can
        also be used to test is the input (di)graph has vertex separation at
        most ``upper_bound`` or to return the first found solution with vertex
        separation less or equal to a ``cut_off`` value.

      - ``exponential`` -- Use an exponential time and space algorithm. This
        algorithm only works of graphs on less than 32 vertices.

      - ``MILP`` -- Use a mixed integer linear programming formulation. This
        algorithm has no size restriction but could take a very long time.

    - ``upper_bound`` -- integer (default: ``None``); parameter used by the
      ``"BAB"`` algorithm. If specified, the algorithm searches for a solution
      with ``width < upper_bound``. It helps cutting branches.  However, if the
      given upper bound is too low, the algorithm may not be able to find a
      solution.

    - ``cut_off`` -- integer (default: ``None``); parameter used by the
      ``"BAB"`` algorithm. This bound allows us to stop the search as soon as a
      solution with width at most ``cut_off`` is found, if any. If this bound
      cannot be reached, the best solution found is returned, unless a too low
      ``upper_bound`` is given.

    - ``verbose`` -- boolean (default: ``False``); whether to display
      information on the computations

    - ``max_prefix_length`` -- integer (default: 20); limits the length of the
      stored prefixes to prevent storing too many prefixes. This parameter is
      used only when ``algorithm=="BAB"``.

    - ``max_prefix_number`` -- integer (default: 10**6); upper bound on the
      number of stored prefixes used to prevent using too much memory. This
      parameter is used only when ``algorithm=="BAB"``.

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.

    .. SEEALSO::

        * :meth:`Graph.treewidth` -- computes the treewidth of a graph

    EXAMPLES:

    The pathwidth of a cycle is equal to 2::

        sage: from sage.graphs.graph_decompositions.vertex_separation import path_decomposition
        sage: g = graphs.CycleGraph(6)
        sage: pw, L = path_decomposition(g, algorithm = "BAB"); pw
        2
        sage: pw, L = path_decomposition(g, algorithm = "exponential"); pw
        2
        sage: pw, L = path_decomposition(g, algorithm="MILP"); pw                       # needs sage.numerical.mip
        2

    TESTS:

    Given anything else than a Graph::

        sage: from sage.graphs.graph_decompositions.vertex_separation import path_decomposition
        sage: path_decomposition(DiGraph())
        Traceback (most recent call last):
        ...
        ValueError: the parameter must be a Graph

    Given a wrong algorithm::

        sage: from sage.graphs.graph_decompositions.vertex_separation import path_decomposition
        sage: path_decomposition(Graph(), algorithm="SuperFast")
        Traceback (most recent call last):
        ...
        ValueError: algorithm "SuperFast" has not been implemented yet, please contribute

    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("the parameter must be a Graph")

    return vertex_separation(G, algorithm=algorithm, cut_off=cut_off, upper_bound=upper_bound,
                             verbose=verbose, max_prefix_length=max_prefix_length,
                             max_prefix_number=max_prefix_number)


def vertex_separation(G, algorithm="BAB", cut_off=None, upper_bound=None, verbose=False,
                      max_prefix_length=20, max_prefix_number=10**6,
                      *, solver=None, integrality_tolerance=1e-3):
    r"""
    Return an optimal ordering of the vertices and its cost for
    vertex-separation.

    INPUT:

    - ``G`` -- a Graph or a DiGraph

    - ``algorithm`` -- string (default: ``"BAB"``); algorithm to use among:

      - ``"BAB"`` -- Use a branch-and-bound algorithm. This algorithm has no
        size restriction but could take a very long time on large graphs. It can
        also be used to test is the input (di)graph has vertex separation at
        most ``upper_bound`` or to return the first found solution with vertex
        separation less or equal to a ``cut_off`` value.

      - ``exponential`` -- Use an exponential time and space algorithm. This
        algorithm only works of graphs on less than 32 vertices.

      - ``MILP`` -- Use a mixed integer linear programming formulation. This
        algorithm has no size restriction but could take a very long time.

    - ``upper_bound`` -- integer (default: ``None``); parameter used by the
      ``"BAB"`` algorithm. If specified, the algorithm searches for a solution
      with ``width < upper_bound``. It helps cutting branches.  However, if the
      given upper bound is too low, the algorithm may not be able to find a
      solution.

    - ``cut_off`` -- integer (default: ``None``); parameter used by the
      ``"BAB"`` algorithm. This bound allows us to stop the search as soon as a
      solution with width at most ``cut_off`` is found, if any. If this bound
      cannot be reached, the best solution found is returned, unless a too low
      ``upper_bound`` is given.

    - ``verbose`` -- boolean (default: ``False``); whether to display
      information on the computations

    - ``max_prefix_length`` -- integer (default: 20); limits the length of the
      stored prefixes to prevent storing too many prefixes. This parameter is
      used only when ``algorithm=="BAB"``.

    - ``max_prefix_number`` -- integer (default: 10**6); upper bound on the
      number of stored prefixes used to prevent using too much memory. This
      parameter is used only when ``algorithm=="BAB"``.

    - ``solver`` -- string (default: ``None``); specify a Mixed Integer Linear
      Programming (MILP) solver to be used. If set to ``None``, the default one
      is used. For more information on MILP solvers and which default solver is
      used, see the method :meth:`solve
      <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
      :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``integrality_tolerance`` -- float; parameter for use with MILP solvers
      over an inexact base ring; see
      :meth:`MixedIntegerLinearProgram.get_values`.

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.

    EXAMPLES:

    Comparison of methods::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation

        sage: # needs sage.combinat
        sage: G = digraphs.DeBruijn(2,3)
        sage: vs,L = vertex_separation(G, algorithm="BAB"); vs
        2
        sage: vs,L = vertex_separation(G, algorithm="exponential"); vs
        2
        sage: vs,L = vertex_separation(G, algorithm="MILP"); vs                         # needs sage.numerical.mip
        2

        sage: G = graphs.Grid2dGraph(3,3)
        sage: vs,L = vertex_separation(G, algorithm="BAB"); vs
        3
        sage: vs,L = vertex_separation(G, algorithm="exponential"); vs
        3
        sage: vs,L = vertex_separation(G, algorithm="MILP"); vs                         # needs sage.numerical.mip
        3

    Digraphs with multiple strongly connected components::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation
        sage: D = digraphs.Path(8)
        sage: print(vertex_separation(D))
        (0, [7, 6, 5, 4, 3, 2, 1, 0])
        sage: D = digraphs.RandomDirectedAcyclicGraph(10, .5)
        sage: vs,L = vertex_separation(D); vs
        0
        sage: K4 = DiGraph( graphs.CompleteGraph(4) )
        sage: D = K4+K4
        sage: D.add_edge(0, 4)
        sage: print(vertex_separation(D))
        (3, [4, 5, 6, 7, 0, 1, 2, 3])
        sage: D = K4+K4+K4
        sage: D.add_edge(0, 4)
        sage: D.add_edge(0, 8)
        sage: print(vertex_separation(D))
        (3, [10, 11, 8, 9, 4, 5, 6, 7, 0, 1, 2, 3])

    Using a specific MILP solver::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation
        sage: G = graphs.PetersenGraph()
        sage: vs, L = vertex_separation(G, algorithm="MILP", solver="SCIP"); vs     # optional - pyscipopt, needs sage.numerical.mip
        5

    TESTS:

    Given a wrong algorithm::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation
        sage: vertex_separation(Graph(), algorithm="SuperFast")
        Traceback (most recent call last):
        ...
        ValueError: algorithm "SuperFast" has not been implemented yet, please contribute

    Given anything else than a Graph or a DiGraph::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation
        sage: vertex_separation(range(4))
        Traceback (most recent call last):
        ...
        ValueError: the parameter must be a Graph or a DiGraph
    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph

    CC = []
    if isinstance(G, Graph):
        if not G.is_connected():
            # We decompose the graph into connected components.
            CC = G.connected_components()

    elif isinstance(G, DiGraph):
        if not G.is_strongly_connected():
            # We decompose the digraph into strongly connected components and
            # arrange them in the inverse order of the topological sort of the
            # digraph of the strongly connected components.
            scc_digraph = G.strongly_connected_components_digraph()
            CC = scc_digraph.topological_sort()[::-1]

    else:
        raise ValueError('the parameter must be a Graph or a DiGraph')

    if cut_off is None:
        cut_off = 0

    if CC:
        # The graph has several (strongly) connected components. We solve the
        # problem on each of them and order partial solutions in the same order
        # than in list CC. The vertex separation is the maximum over all these
        # subgraphs.
        vs, L = 0, []
        for V in CC:

            if len(V) == 1:
                # We can directly add this vertex to the solution
                L.extend(V)

            else:
                # We build the (strongly) connected subgraph and do a recursive
                # call to get its vertex separation and corresponding ordering
                H = G.subgraph(V)
                vsH, LH = vertex_separation(H, algorithm=algorithm,
                                            cut_off=cut_off,
                                            upper_bound=upper_bound,
                                            verbose=verbose,
                                            max_prefix_length=max_prefix_length,
                                            max_prefix_number=max_prefix_number,
                                            solver=solver,
                                            integrality_tolerance=integrality_tolerance)

                if vsH == -1:
                    # We have not been able to find a solution. This case
                    # happens when a too low upper bound is given.
                    return -1, []

                # We update the vertex separation and ordering
                vs = max(vs, vsH)
                L.extend(LH)

                # We also update the cut_off parameter that could speed up
                # resolution for other components (used when algorithm=="BAB")
                cut_off = max(cut_off, vs)

        return vs, L

    # We have a (strongly) connected graph and we call the desired algorithm
    if algorithm == "exponential":
        return vertex_separation_exp(G, verbose=verbose)

    elif algorithm == "MILP":
        return vertex_separation_MILP(G, solver=solver, verbose=verbose,
                                      integrality_tolerance=integrality_tolerance)

    elif algorithm == "BAB":
        return vertex_separation_BAB(G, cut_off=cut_off, upper_bound=upper_bound, verbose=verbose,
                                     max_prefix_length=max_prefix_length, max_prefix_number=max_prefix_number)

    else:
        raise ValueError('algorithm "{}" has not been implemented yet, please contribute'.format(algorithm))


################################
# Exact exponential algorithms #
################################

def vertex_separation_exp(G, verbose=False):
    r"""
    Return an optimal ordering of the vertices and its cost for
    vertex-separation.

    INPUT:

    - ``G`` -- a Graph or a DiGraph

    - ``verbose`` -- boolean (default: ``False``); whether to display
      information on the computations

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.

    .. NOTE::

        Because of its current implementation, this algorithm only works on
        graphs on less than 32 vertices. This can be changed to 54 if necessary,
        but 32 vertices already require 4GB of memory.

    EXAMPLES:

    The vertex separation of a circuit is equal to 1::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation_exp
        sage: g = digraphs.Circuit(6)
        sage: vertex_separation_exp(g)
        (1, [0, 1, 2, 3, 4, 5])

    TESTS:

    Given anything else than a Graph or a DiGraph::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation_exp
        sage: vertex_separation_exp(range(3))
        Traceback (most recent call last):
        ...
        ValueError: the parameter must be a Graph or a DiGraph

    Graphs with non-integer vertices::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation_exp
        sage: D = digraphs.DeBruijn(2,3)                                                # needs sage.combinat
        sage: vertex_separation_exp(D)                                                  # needs sage.combinat
        (2, ['000', '001', '100', '010', '101', '011', '110', '111'])

    Given a too large graph::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation_exp
        sage: vertex_separation_exp(graphs.PathGraph(50))
        Traceback (most recent call last):
        ...
        ValueError: the graph should have at most 31 vertices
    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, Graph) and not isinstance(G, DiGraph):
        raise ValueError("the parameter must be a Graph or a DiGraph")

    if G.order() >= 32:
        raise ValueError("the graph should have at most 31 vertices")

    cdef FastDigraph g = FastDigraph(G)

    if verbose:
        print("Memory allocation")
        g.print_adjacency_matrix()

    cdef unsigned int mem = 1 << g.n
    cdef uint8_t * neighborhoods = <uint8_t *>check_malloc(mem)

    memset(neighborhoods, <uint8_t> -1, mem)

    cdef int i, k
    for k in range(g.n):
        if verbose:
            print("Looking for a strategy of cost", str(k))

        sig_check()
        if exists(g, neighborhoods, 0, k) <= k:
            break

    if verbose:
        print("... Found !")
        print("Now computing the ordering")

    cdef list order = find_order(g, neighborhoods, k)

    sig_free(neighborhoods)

    return k, [g.int_to_vertices[i] for i in order]


###############################################################################
# Actual algorithm, breadth-first search and updates of the costs of the sets #
###############################################################################

cdef inline int exists(FastDigraph g, uint8_t* neighborhoods, int current, int cost):
    """
    Check whether an ordering with the given cost exists, and updates data in
    the neighborhoods array at the same time. See the module's documentation.
    """
    # If this is true, it means the set has not been evaluated yet
    if neighborhoods[current] == <uint8_t> -1:
        neighborhoods[current] = compute_out_neighborhood_cardinality(g, current)

    # If the cost of this set is too high, there is no point in going further.
    # Same thing if the current set is the whole vertex set.
    if neighborhoods[current] > cost or (current == (1 << g.n) - 1):
        return neighborhoods[current]

    # Minimum of the costs of the outneighbors
    cdef int mini = g.n

    cdef int i
    cdef int next_set

    for i in range(g.n):
        if (current >> i) & 1:
            continue

        # For each of the out-neighbors next_set of current
        next_set = current | 1 << i

        # Check whether there exists a cheap path toward {1..n}, and updated the
        # cost.
        mini = minimum(mini, exists(g, neighborhoods, next_set, cost))

        # We have found a path !
        if mini <= cost:
            return mini

    # Updating the cost of the current set with the minimum of the cost of its
    # outneighbors.
    neighborhoods[current] = mini

    return neighborhoods[current]


cdef list find_order(FastDigraph g, uint8_t* neighborhoods, int cost):
    """
    Return the ordering once we are sure it exists
    """
    cdef list ordering = []
    cdef int current = 0
    cdef int n = g.n
    cdef int i

    while n:
        # We look for n vertices
        for i in range(g.n):
            if (current >> i) & 1:
                continue

            # Find the next set with small cost (we know it exists)
            next_set = current | 1 << i
            if neighborhoods[next_set] <= cost:
                ordering.append(i)
                current = next_set
                break

        # One less to find
        n -= 1

    return ordering


# Min/Max functions

cdef inline int minimum(int a, int b):
    if a < b:
        return a
    else:
        return b


cdef inline int maximum(int a, int b):
    if a > b:
        return a
    else:
        return b


#################################################################
# Function for testing the validity of a linear vertex ordering #
#################################################################

def is_valid_ordering(G, L):
    r"""
    Test if the linear vertex ordering `L` is valid for (di)graph `G`.

    A linear ordering `L` of the vertices of a (di)graph `G` is valid if all
    vertices of `G` are in `L`, and if `L` contains no other vertex and no
    duplicated vertices.

    INPUT:

    - ``G`` -- a Graph or a DiGraph.

    - ``L`` -- an ordered list of the vertices of ``G``.


    OUTPUT:

    Returns ``True`` if `L` is a valid vertex ordering for `G`, and ``False``
    otherwise.


    EXAMPLES:

    Path decomposition of a cycle::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = graphs.CycleGraph(6)
        sage: L = G.vertices(sort=True)
        sage: vertex_separation.is_valid_ordering(G, L)
        True
        sage: vertex_separation.is_valid_ordering(G, [1,2])
        False

    TESTS:

    Giving anything else than a Graph or a DiGraph::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: vertex_separation.is_valid_ordering(2, [])
        Traceback (most recent call last):
        ...
        ValueError: the input parameter must be a Graph or a DiGraph

    Giving anything else than a list::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = graphs.CycleGraph(6)
        sage: vertex_separation.is_valid_ordering(G, {})
        Traceback (most recent call last):
        ...
        ValueError: the second parameter must be of type 'list'
    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, Graph) and not isinstance(G, DiGraph):
        raise ValueError("the input parameter must be a Graph or a DiGraph")
    if not isinstance(L, list):
        raise ValueError("the second parameter must be of type 'list'")

    return set(L) == set(G)


####################################################################
# Measurement functions of the widths of some graph decompositions #
####################################################################

def width_of_path_decomposition(G, L):
    r"""
    Return the width of the path decomposition induced by the linear ordering
    `L` of the vertices of `G`.

    If `G` is an instance of :mod:`Graph <sage.graphs.graph>`, this function
    returns the width `pw_L(G)` of the path decomposition induced by the linear
    ordering `L` of the vertices of `G`. If `G` is a :mod:`DiGraph
    <sage.graphs.digraph>`, it returns instead the width `vs_L(G)` of the
    directed path decomposition induced by the linear ordering `L` of the
    vertices of `G`, where

    .. MATH::

        vs_L(G) & =  \max_{0\leq i< |V|-1} | N^+(L[:i])\setminus L[:i] |\\
        pw_L(G) & =  \max_{0\leq i< |V|-1} | N(L[:i])\setminus L[:i] |\\

    INPUT:

    - ``G`` -- a Graph or a DiGraph

    - ``L`` -- a linear ordering of the vertices of ``G``

    EXAMPLES:

    Path decomposition of a cycle::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = graphs.CycleGraph(6)
        sage: L = G.vertices(sort=True)
        sage: vertex_separation.width_of_path_decomposition(G, L)
        2

    Directed path decomposition of a circuit::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = digraphs.Circuit(6)
        sage: L = G.vertices(sort=True)
        sage: vertex_separation.width_of_path_decomposition(G, L)
        1

    TESTS:

    Path decomposition of a BalancedTree::

        sage: # needs networkx
        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = graphs.BalancedTree(3,2)
        sage: pw, L = vertex_separation.path_decomposition(G)
        sage: pw == vertex_separation.width_of_path_decomposition(G, L)
        True
        sage: L.reverse()
        sage: pw == vertex_separation.width_of_path_decomposition(G, L)
        False

    Directed path decomposition of a circuit::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = digraphs.Circuit(8)
        sage: vs, L = vertex_separation.vertex_separation(G)
        sage: vs == vertex_separation.width_of_path_decomposition(G, L)
        True
        sage: L = [0,4,6,3,1,5,2,7]
        sage: vs == vertex_separation.width_of_path_decomposition(G, L)
        False

    Giving a wrong linear ordering::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = Graph()
        sage: vertex_separation.width_of_path_decomposition(G, ['a','b'])
        Traceback (most recent call last):
        ...
        ValueError: the input linear vertex ordering L is not valid for G
    """
    if not is_valid_ordering(G, L):
        raise ValueError("the input linear vertex ordering L is not valid for G")

    neighbors = G.neighbors_out if G.is_directed() else G.neighbors

    cdef int vsL = 0
    cdef set S = set()
    cdef set neighbors_of_S_in_V_minus_S = set()

    for u in L:

        # We remove u from the neighbors of S
        neighbors_of_S_in_V_minus_S.discard(u)

        # We add vertex u to the set S
        S.add(u)

        # We add the (out-)neighbors of u to the neighbors of S
        for v in neighbors(u):
            if v not in S:
                neighbors_of_S_in_V_minus_S.add(v)

        # We update the cost of the vertex separation
        vsL = max(vsL, len(neighbors_of_S_in_V_minus_S))

    return vsL


##########################################
# MILP formulation for vertex separation #
##########################################

def _vertex_separation_MILP_formulation(G, integrality=False, solver=None):
    r"""
    MILP formulation of the vertex separation of `G` and the optimal ordering of its vertices.

    This MILP is an improved version of the formulation proposed in [SP2010]_. See the
    :mod:`module's documentation <sage.graphs.graph_decompositions.vertex_separation>` for
    more details on this MILP formulation.

    INPUT:

    - ``G`` -- a Graph or a DiGraph

    - ``integrality`` -- boolean (default: ``False``); specify if variables
      `x_v^t` and `u_v^t` must be integral or if they can be relaxed. This has
      no impact on the validity of the solution, but it is sometimes faster to
      solve the problem using binary variables only.

    - ``solver`` -- string (default: ``None``); specify a Mixed Integer Linear
      Programming (MILP) solver to be used. If set to ``None``, the default one
      is used. For more information on MILP solvers and which default solver is
      used, see the method :meth:`solve
      <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
      :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    OUTPUT:

    - the :class:`~sage.numerical.mip.MixedIntegerLinearProgram`

    - :class:`sage.numerical.mip.MIPVariable` objects ``x``, ``u``, ``y``, ``z``.

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.vertex_separation import _vertex_separation_MILP_formulation
        sage: G = digraphs.DeBruijn(2,3)                                                # needs sage.combinat
        sage: p, x, u, y, z = _vertex_separation_MILP_formulation(G)                    # needs sage.combinat sage.numerical.mip
        sage: p                                                                         # needs sage.combinat sage.numerical.mip
        Mixed Integer Program (minimization, 193 variables, 449 constraints)
    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, Graph) and not isinstance(G, DiGraph):
        raise ValueError("the first input parameter must be a Graph or a DiGraph")

    from sage.numerical.mip import MixedIntegerLinearProgram
    p = MixedIntegerLinearProgram(maximization=False, solver=solver)

    # Declaration of variables.
    x = p.new_variable(binary=integrality, nonnegative=True)
    u = p.new_variable(binary=integrality, nonnegative=True)
    y = p.new_variable(binary=True)
    z = p.new_variable(integer=True, nonnegative=True)

    N = G.order()
    V = list(G)
    neighbors_out = G.neighbors_out if G.is_directed() else G.neighbors

    # (2) x[v,t] <= x[v,t+1]   for all v in V, and for t:=0..N-2
    # (3) y[v,t] <= y[v,t+1]   for all v in V, and for t:=0..N-2
    for v in V:
        for t in range(N - 1):
            p.add_constraint(x[v, t] - x[v, t + 1] <= 0)
            p.add_constraint(y[v, t] - y[v, t + 1] <= 0)

    # (4) y[v,t] <= x[w,t]  for all v in V, for all w in N^+(v), and for all t:=0..N-1
    for v in V:
        for w in neighbors_out(v):
            for t in range(N):
                p.add_constraint(y[v, t] - x[w, t] <= 0)

    # (5) sum_{v in V} y[v,t] == t+1 for t:=0..N-1
    for t in range(N):
        p.add_constraint(p.sum(y[v, t] for v in V) == t + 1)

    # (6) u[v,t] >= x[v,t]-y[v,t]    for all v in V, and for all t:=0..N-1
    for v in V:
        for t in range(N):
            p.add_constraint(x[v, t] - y[v, t] - u[v, t] <= 0)

    # (7) z >= sum_{v in V} u[v,t]   for all t:=0..N-1
    for t in range(N):
        p.add_constraint(p.sum(u[v, t] for v in V) - z['z'] <= 0)

    # (8)(9) 0 <= x[v,t] and u[v,t] <= 1
    if not integrality:
        for v in V:
            for t in range(N):
                p.add_constraint(x[v, t], min=0, max=1)
                p.add_constraint(u[v, t], min=0, max=1)

    # (10) y[v,t] in {0,1}
    # already declared

    # (11) 0 <= z <= |V|
    p.add_constraint(z['z'] <= N)

    #  (1) Minimize z
    p.set_objective(z['z'])

    return p, x, u, y, z


def vertex_separation_MILP(G, integrality=False, solver=None, verbose=0,
                           *, integrality_tolerance=1e-3):
    r"""
    Compute the vertex separation of `G` and the optimal ordering of its
    vertices using an MILP formulation.

    This function uses a mixed integer linear program (MILP) for determining an
    optimal layout for the vertex separation of `G`. This MILP is an improved
    version of the formulation proposed in [SP2010]_. See the :mod:`module's
    documentation <sage.graphs.graph_decompositions.vertex_separation>` for more
    details on this MILP formulation.

    INPUT:

    - ``G`` -- a Graph or a DiGraph

    - ``integrality`` -- boolean (default: ``False``); specify if variables
      `x_v^t` and `u_v^t` must be integral or if they can be relaxed. This has
      no impact on the validity of the solution, but it is sometimes faster to
      solve the problem using binary variables only.

    - ``solver`` -- string (default: ``None``); specify a Mixed Integer Linear
      Programming (MILP) solver to be used. If set to ``None``, the default one
      is used. For more information on MILP solvers and which default solver is
      used, see the method :meth:`solve
      <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
      :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- integer (default: ``0``); sets the level of verbosity. Set
      to 0 by default, which means quiet.

    - ``integrality_tolerance`` -- float; parameter for use with MILP solvers
      over an inexact base ring; see
      :meth:`MixedIntegerLinearProgram.get_values`.

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.

    EXAMPLES:

    Vertex separation of a De Bruijn digraph::

        sage: # needs sage.combinat
        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = digraphs.DeBruijn(2,3)
        sage: vs, L = vertex_separation.vertex_separation_MILP(G); vs                   # needs sage.numerical.mip
        2
        sage: vs == vertex_separation.width_of_path_decomposition(G, L)                 # needs sage.numerical.mip
        True
        sage: vse, Le = vertex_separation.vertex_separation(G); vse
        2

    The vertex separation of a circuit is 1::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = digraphs.Circuit(6)
        sage: vs, L = vertex_separation.vertex_separation_MILP(G); vs                   # needs sage.numerical.mip
        1

    TESTS:

    Comparison with exponential algorithm::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: for i in range(10):                                                       # needs sage.numerical.mip
        ....:     G = digraphs.RandomDirectedGNP(10, 0.2)
        ....:     ve, le = vertex_separation.vertex_separation(G)
        ....:     vm, lm = vertex_separation.vertex_separation_MILP(G)
        ....:     if ve != vm:
        ....:        raise ValueError("the solution is not optimal")

    Comparison with different values of the integrality parameter::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: for i in range(10):               # long time (11s on sage.math, 2012), needs sage.numerical.mip
        ....:     G = digraphs.RandomDirectedGNP(10, 0.2)
        ....:     va, la = vertex_separation.vertex_separation_MILP(G, integrality=False)
        ....:     vb, lb = vertex_separation.vertex_separation_MILP(G, integrality=True)
        ....:     if va != vb:
        ....:        raise ValueError("the integrality parameter changes the result")

    Giving anything else than a Graph or a DiGraph::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: vertex_separation.vertex_separation_MILP([])                              # needs sage.numerical.mip
        Traceback (most recent call last):
        ...
        ValueError: the first input parameter must be a Graph or a DiGraph
    """
    from sage.numerical.mip import MIPSolverException

    p, _, _, y, z = _vertex_separation_MILP_formulation(G, integrality=integrality, solver=solver)
    N = G.order()
    V = list(G)

    try:
        _ = p.solve(log=verbose)
    except MIPSolverException:
        if integrality:
            raise ValueError("unbounded or unexpected error")
        else:
            raise ValueError("unbounded or unexpected error, try with 'integrality = True'")

    taby = p.get_values(y, convert=bool, tolerance=integrality_tolerance)
    vs = p.get_values(z, convert=True, tolerance=integrality_tolerance)['z']
    # since exactly one vertex is processed per step, we can reconstruct the sequence
    seq = []
    seen = set()
    for t in range(N):
        for v in V:
            if taby[v, t] and v not in seen:
                seq.append(v)
                seen.add(v)
                break

    return vs, seq


##########################################
# Branch and Bound for vertex separation #
##########################################

def vertex_separation_BAB(G,
                          cut_off=None,
                          upper_bound=None,
                          max_prefix_length=20,
                          max_prefix_number=10**6,
                          verbose=False):
    r"""
    Branch and Bound algorithm for the vertex separation.

    This method implements the branch and bound algorithm for the vertex
    separation of directed graphs and the pathwidth of undirected graphs
    proposed in [CMN2014]_. The implementation is valid for both Graph and
    DiGraph. See the documentation of the
    :mod:`~sage.graphs.graph_decompositions.vertex_separation` module.

    INPUT:

    - ``G`` -- a Graph or a DiGraph.

    - ``cut_off`` -- integer (default: ``None``); bound to consider in the
      branch and bound algorithm. This allows us to stop the search as soon as a
      solution with width at most ``cut_off`` is found, if any. If this bound
      cannot be reached, the best solution found is returned, unless a too low
      ``upper_bound`` is given.

    - ``upper_bound`` -- integer (default: ``None``); if specified, the
      algorithm searches for a solution with ``width < upper_bound``. It helps
      cutting branches.  However, if the given upper bound is too low, the
      algorithm may not be able to find a solution.

    - ``max_prefix_length`` -- integer (default: 20); limits the length of the
      stored prefixes to prevent storing too many prefixes

    - ``max_prefix_number`` -- integer (default: 10**6); upper bound on the
      number of stored prefixes used to prevent using too much memory

    - ``verbose`` -- boolean (default: ``False``); display some information when
      set to ``True``

    OUTPUT:

    - ``width`` -- the computed vertex separation

    - ``seq`` -- an ordering of the vertices of width ``width``

    EXAMPLES:

    The algorithm is valid for the vertex separation::

        sage: from sage.graphs.graph_decompositions import vertex_separation as VS
        sage: D = digraphs.RandomDirectedGNP(15, .2)
        sage: vb, seqb = VS.vertex_separation_BAB(D)
        sage: vd, seqd = VS.vertex_separation_exp(D)
        sage: vb == vd
        True
        sage: vb == VS.width_of_path_decomposition(D, seqb)
        True

    The vertex separation of a `N\times N` grid is `N`::

        sage: from sage.graphs.graph_decompositions import vertex_separation as VS
        sage: G = graphs.Grid2dGraph(4,4)
        sage: vs, seq = VS.vertex_separation_BAB(G); vs
        4
        sage: vs == VS.width_of_path_decomposition(G, seq)
        True

    The vertex separation of a `N\times M` grid with `N<M` is `N`::

        sage: from sage.graphs.graph_decompositions import vertex_separation as VS
        sage: G = graphs.Grid2dGraph(3,5)
        sage: vs, seq = VS.vertex_separation_BAB(G); vs
        3
        sage: vs == VS.width_of_path_decomposition(G, seq)
        True

    The vertex separation of circuit of order `N\geq 2` is 1::

        sage: from sage.graphs.graph_decompositions import vertex_separation as VS
        sage: D = digraphs.Circuit(10)
        sage: vs, seq = VS.vertex_separation_BAB(D); vs
        1
        sage: vs == VS.width_of_path_decomposition(D, seq)
        True

    The vertex separation of cycle of order `N\geq 3` is 2::

        sage: from sage.graphs.graph_decompositions import vertex_separation as VS
        sage: G = graphs.CycleGraph(10)
        sage: vs, seq = VS.vertex_separation_BAB(G); vs
        2

    The vertex separation of ``MycielskiGraph(5)`` is 10::

        sage: from sage.graphs.graph_decompositions import vertex_separation as VS
        sage: G = graphs.MycielskiGraph(5)
        sage: vs, seq = VS.vertex_separation_BAB(G); vs
        10

    Searching for any solution with width less or equal to ``cut_off``::

        sage: from sage.graphs.graph_decompositions import vertex_separation as VS
        sage: G = graphs.MycielskiGraph(5)
        sage: VS.vertex_separation_BAB(G, cut_off=11)[0] <= 11
        True
        sage: VS.vertex_separation_BAB(G, cut_off=10)[0] <= 10
        True
        sage: VS.vertex_separation_BAB(G, cut_off=9)[0] <= 9
        False

    Testing for the existence of a solution with width strictly less than ``upper_bound``::

        sage: from sage.graphs.graph_decompositions import vertex_separation as VS
        sage: G = graphs.MycielskiGraph(5)
        sage: vs, seq = VS.vertex_separation_BAB(G, upper_bound=11); vs
        10
        sage: vs, seq = VS.vertex_separation_BAB(G, upper_bound=10); vs
        -1
        sage: vs, seq = VS.vertex_separation_BAB(G, cut_off=11, upper_bound=10); vs
        -1

    Changing the parameters of the prefix storage::

        sage: from sage.graphs.graph_decompositions import vertex_separation as VS
        sage: G = graphs.MycielskiGraph(5)
        sage: vs, seq = VS.vertex_separation_BAB(G, max_prefix_length=0); vs
        10
        sage: vs, seq = VS.vertex_separation_BAB(G, max_prefix_number=5); vs
        10
        sage: vs, seq = VS.vertex_separation_BAB(G, max_prefix_number=0); vs
        10

    TESTS:

    Giving anything else than a Graph or a DiGraph::

        sage: from sage.graphs.graph_decompositions import vertex_separation as VS
        sage: VS.vertex_separation_BAB(range(5))
        Traceback (most recent call last):
        ...
        ValueError: the input parameter must be a Graph or a DiGraph

    Giving an empty Graph or DiGraph::

        sage: from sage.graphs.graph_decompositions import vertex_separation as VS
        sage: VS.vertex_separation_BAB(Graph())
        (0, [])
        sage: VS.vertex_separation_BAB(DiGraph())
        (0, [])

    Giving a too low upper bound::

        sage: from sage.graphs.graph_decompositions import vertex_separation as VS
        sage: VS.vertex_separation_BAB(digraphs.Circuit(3), upper_bound=0)
        Traceback (most recent call last):
        ...
        ValueError: the input upper bound must be at least 1
    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, DiGraph) and not isinstance(G, Graph):
        raise ValueError("the input parameter must be a Graph or a DiGraph")

    cdef int n = G.order()
    if not n:
        return 0, []

    cut_off = 0 if cut_off is None else cut_off
    upper_bound = n if upper_bound is None else upper_bound
    if upper_bound < 1:
        raise ValueError("the input upper bound must be at least 1")

    # ==> Allocate and initialize some data structures

    # We use a binary matrix to store the (di)graph. This way the neighborhood
    # of a vertex is stored in one bitset.
    cdef binary_matrix_t H
    cdef int i
    cdef list int_to_vertex = list(G)
    cdef dict vertex_to_int = {v: i for i, v in enumerate(int_to_vertex)}
    dense_graph_init(H, G, translation=vertex_to_int)

    # We need 2 bitsets here + 3 per call to vertex_separation_BAB_C, so overall
    # 3*n + 2. We use another binary matrix as a pool of bitsets.
    cdef binary_matrix_t bm_pool
    binary_matrix_init(bm_pool, 3 * n + 2, n)

    cdef int * prefix = <int *>sig_malloc(n * sizeof(int))
    cdef int * positions = <int *>sig_malloc(n * sizeof(int))
    if not prefix or not positions:
        sig_free(prefix)
        sig_free(positions)
        binary_matrix_free(H)
        binary_matrix_free(bm_pool)
        raise MemoryError("unable to allocate data structures")

    cdef list best_seq = list(range(n))
    for i in range(n):
        prefix[i] = i
        positions[i] = i

    cdef int width = upper_bound
    cdef list order = []
    cdef set prefix_storage = set()

    try:
        # ==> Call the cython method
        sig_on()
        width = vertex_separation_BAB_C(H=H,
                                        n=n,
                                        prefix=prefix,
                                        positions=positions,
                                        best_seq=best_seq,
                                        level=0,
                                        b_prefix=bm_pool.rows[3 * n],
                                        b_prefix_and_neighborhood=bm_pool.rows[3 * n + 1],
                                        cut_off=cut_off,
                                        upper_bound=upper_bound,
                                        current_cost=0,
                                        bm_pool=bm_pool,
                                        prefix_storage=prefix_storage,
                                        max_prefix_length=max_prefix_length,
                                        max_prefix_number=max_prefix_number,
                                        verbose=verbose)

        sig_off()

        # ==> Build the final ordering
        order = [int_to_vertex[best_seq[i]] for i in range(n)]

    finally:
        if verbose:
            print('Stored prefixes: {}'.format(len(prefix_storage)))
        sig_free(prefix)
        sig_free(positions)
        binary_matrix_free(H)
        binary_matrix_free(bm_pool)

    return (width if width < upper_bound else -1), order


cdef inline _my_invert_positions(int *prefix, int *positions, int pos_a, int pos_b):
    """
    Permute vertices at positions ``pos_a`` and ``pos_b`` in array ``prefix``,
    and record the new positions in array ``positions``.
    """
    if pos_a != pos_b:
        positions[prefix[pos_a]], positions[prefix[pos_b]] = positions[prefix[pos_b]], positions[prefix[pos_a]]
        prefix[pos_a], prefix[pos_b] = prefix[pos_b], prefix[pos_a]


cdef int vertex_separation_BAB_C(binary_matrix_t H,
                                 int n,
                                 int *prefix,
                                 int *positions,
                                 list best_seq,
                                 int level,
                                 bitset_t b_prefix,
                                 bitset_t b_prefix_and_neighborhood,
                                 int cut_off,
                                 int upper_bound,
                                 int current_cost,
                                 binary_matrix_t bm_pool,
                                 set prefix_storage,
                                 int max_prefix_length,
                                 int max_prefix_number,
                                 bint verbose):
    r"""
    Branch and Bound algorithm for the process number and the vertex separation.

    INPUT:

    - ``H`` -- a binary matrix storing the adjacency of the (di)graph

    - ``n`` -- integer; the number of vertices of the (di)graph

    - ``prefix`` -- array of `n` integers; contains a permutation of the
      vertices. The vertices forming the current prefix under consideration are
      stored in cells ``[0, level - 1]``.

    - ``positions`` -- array of `n` integers; associates to each vertex its
      index in array ``prefix``

    - ``best_seq`` -- array of `n` integers; stores the best ordering found so
      far

    - ``level`` -- integer; specifies the length of the current prefix

    - ``b_prefix`` -- bitset of size `n`; records the vertices in the current
      prefix (in cells ``[0,level-1]``)

    - ``b_prefix_and_neighborhood`` -- bitset of size `n`; records the
      vertices in the current prefix and the vertices in its neighborhood

    - ``cut_off`` -- integer; bound to consider in the branch and bound
      algorithm. This allows us to stop the search as soon as a solution with
      width at most ``cut_off`` is found, if any.

    - ``upper_bound`` -- integer; the algorithm searches for a solution with
      ``width < upper_bound``. It helps cutting branches. Each time a new
      solution is found, the upper bound is reduced.

    - ``bm_pool`` -- a binary matrix with `3*n+2` rows of size `n`; each rows is
      a bitset of size `n`. This data structure is used as a pool of
      initialized bitsets. Each call of this method needs 3 bitsets for local
      operations, so it uses rows ``[3 * level, 3 * level + 2]``.

    - ``prefix_storage`` -- set; used to store prefixes

    - ``max_prefix_length`` -- integer; maximum length of the stored prefixes to
      prevent storing too many prefixes

    - ``max_prefix_number`` -- integer; upper bound on the number of stored
      prefixes used to prevent using too much memory

    - ``verbose`` -- boolean (default: ``False``); display some information when
      set to ``True``
    """
    cdef int i

    # ==> Test termination

    if level == n:
        if current_cost < upper_bound:
            for i in range(n):
                best_seq[i] = prefix[i]
            if verbose:
                print("New upper bound: {}".format(current_cost))

        return current_cost

    cdef int delta_i, j, v, select_it
    cdef list delta = list()
    cdef int loc_level = level

    # ==> Allocate local data structures

    cdef bitset_s *loc_b_prefix = bm_pool.rows[3 * level]
    cdef bitset_s *loc_b_pref_and_neigh = bm_pool.rows[3 * level + 1]
    cdef bitset_s *b_tmp = bm_pool.rows[3 * level + 2]
    bitset_copy(loc_b_prefix, b_prefix)
    bitset_copy(loc_b_pref_and_neigh, b_prefix_and_neighborhood)

    # ==> Greedy steps
    #
    # We extend the current prefix with all vertices u such that either
    # (i) All out-neighbors of u are in the prefix or in its out-neighborhood
    # (ii) or u is an out-neighbor of the prefix and all but one of its
    #      out-neighbors are in the prefix or in its out-neighborhood.

    select_it = 0
    i = loc_level
    while i < n:

        j = prefix[i]

        if bitset_issubset(H.rows[j], loc_b_pref_and_neigh):
            # (i) Vertex j is such that all its out-neighbors are in the prefix
            # or in its out-neighborhood (so in loc_b_pref_and_neigh).
            bitset_add(loc_b_pref_and_neigh, j)
            select_it = 1

        elif bitset_in(loc_b_pref_and_neigh, j) and not bitset_in(loc_b_prefix, j):
            bitset_difference(b_tmp, H.rows[j], loc_b_pref_and_neigh)
            if bitset_len(b_tmp) == 1:
                # (ii) Vertex j is an out-neighbor of the prefix and all but one
                # of its out-neighbors are in the prefix or in its
                # out-neighborhood.
                v = bitset_first(b_tmp)
                bitset_add(loc_b_pref_and_neigh, v)
                select_it = 1

        if select_it:
            # We add j to the prefix and update neighborhoods
            _my_invert_positions(prefix, positions, i, loc_level)
            loc_level += 1
            bitset_add(loc_b_prefix, j)
            select_it = 0
            # We search for vertices that can now be selected
            i = loc_level
        else:
            i += 1

    # ==> Test termination
    #
    if loc_level == n:
        if current_cost < upper_bound:
            for i in range(n):
                best_seq[i] = prefix[i]
            if verbose:
                print("New upper bound: {}".format(current_cost))

        return current_cost

    # ==> Test if the prefix is in prefix_storage
    #
    # The set S of vertices of a prefix P is in prefix_storage if the branch
    # with prefix P is such that c(P)<\min_{L\in{\cal L}_P(V)} c(L). In such
    # case, there is no need to continue exploration for the current branch.
    cdef frozenset frozen_prefix

    if loc_level <= max_prefix_length:
        frozen_prefix = frozenset(prefix[i] for i in range(loc_level))
        if frozen_prefix in prefix_storage:
            return upper_bound

    # ==> Sort and Prune
    #
    # We compute for each remaining vertex v a lower bound on the width of any
    # ordering with prefix prefix+v
    for i in range(loc_level, n):
        j = prefix[i]
        bitset_union(b_tmp, loc_b_pref_and_neigh, H.rows[j])
        bitset_difference(b_tmp, b_tmp, loc_b_prefix)
        bitset_discard(b_tmp, j)
        delta_i = bitset_len(b_tmp)
        if delta_i < upper_bound:
            delta.append((delta_i, j))

    delta.sort()

    # ==> Recursion
    for delta_i, i in delta:

        delta_i = max(current_cost, delta_i)

        if delta_i >= upper_bound:
            break

        # We extend the current prefix with vertex i and explore the branch
        bitset_union(b_tmp, loc_b_pref_and_neigh,  H.rows[i])
        bitset_discard(b_tmp, i)
        _my_invert_positions(prefix, positions, positions[i], loc_level)
        bitset_add(loc_b_prefix, i)

        cost_i = vertex_separation_BAB_C(H=H,
                                         n=n,
                                         prefix=prefix,
                                         positions=positions,
                                         best_seq=best_seq,
                                         level=loc_level + 1,
                                         b_prefix=loc_b_prefix,
                                         b_prefix_and_neighborhood=b_tmp,
                                         cut_off=cut_off,
                                         upper_bound=upper_bound,
                                         current_cost=delta_i,
                                         bm_pool=bm_pool,
                                         prefix_storage=prefix_storage,
                                         max_prefix_length=max_prefix_length,
                                         max_prefix_number=max_prefix_number,
                                         verbose=verbose)

        bitset_discard(loc_b_prefix, i)

        if cost_i < upper_bound:
            upper_bound = cost_i
            if upper_bound <= cut_off:
                # We are satisfied with current solution.
                break

    # ==> Update prefix_storage
    #
    # If the prefix P is such that c(P)<\min_{L\in{\cal L}_P(V)} c(L), no other
    # prefix P' on the same set S=V(P) of vertices can lead to a better
    # solution.
    if (loc_level <= max_prefix_length
            and current_cost < upper_bound
            and len(prefix_storage) < max_prefix_number):
        prefix_storage.add(frozen_prefix)

    return upper_bound
