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
This is a result of Kinnersley [Kin92]_ and Bodlaender [Bod98]_.


**This module contains the following methods**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`path_decomposition` | Returns the pathwidth of the given graph and the ordering of the vertices resulting in a corresponding path decomposition
    :meth:`vertex_separation` | Returns an optimal ordering of the vertices and its cost for vertex-separation
    :meth:`vertex_separation_exp` | Computes the vertex separation of `G` using an exponential time and space algorithm
    :meth:`vertex_separation_MILP` | Computes the vertex separation of `G` and the optimal ordering of its vertices using an MILP formulation
    :meth:`vertex_separation_BAB` | Computes the vertex separation of `G` and the optimal ordering of its vertices using a branch and bound algorithm
    :meth:`lower_bound` | Returns a lower bound on the vertex separation of `G`
    :meth:`is_valid_ordering` | Test if the linear vertex ordering `L` is valid for (di)graph `G`
    :meth:`width_of_path_decomposition` | Returns the width of the path decomposition induced by the linear ordering `L` of the vertices of `G`


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
the formulation proposed in [SP10]_. It aims at building a sequence `S_t` of
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
[CMN14]_.

**Greedy steps:**

Let us denote `{\cal L}(S)` the set of all possible orderings of the vertices in
`S`, and let `{\cal L}_P(S)\subseteq {\cal L}(S)` be the orderings starting with
a prefix `P`. Let also `c(L)` be the cost of the ordering `L\in{\cal L}(V)` as
defined above.

Given a digraph `D=(V,A)`, a set `S\subset V`, and a prefix `P`, it has been
proved in [CMN14]_ that `\min_{L\in{\cal L}_P(V)} c(L) = \min_{L\in{\cal
L}_{P+v}(V)} c(L)` holds in two (non exhaustive) cases:

.. MATH::

    \text{or} \begin{cases}
    N^+(v)\subseteq S\cup N^+(S)\\
    v\in N^+(S)\text{ and }N^+(v)\setminus(S\cup N^+(S)) = \{w\}
    \end{cases}

In other words, if we find a vertex `v` satisfying the above conditions, the best
possible ordering with prefix `P` has the same cost as the best possible
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
vertices) to cut branches later. See [CMN14]_ for more details.

Since the number of stored sets can get very large, one can control the maximum
length and the maximum number of stored prefixes.


REFERENCES
----------

.. [Bod98] *A partial k-arboretum of graphs with bounded treewidth*, Hans
  L. Bodlaender, Theoretical Computer Science 209(1-2):1-45, 1998.

.. [Kin92] *The vertex separation number of a graph equals its path-width*,
  Nancy G. Kinnersley, Information Processing Letters 42(6):345-350, 1992.

.. [SP10] *Lightpath Reconfiguration in WDM networks*, Fernando Solano and
  Michal Pioro, IEEE/OSA Journal of Optical Communication and Networking
  2(12):1010-1021, 2010.

.. [CMN14] *Experimental Evaluation of a Branch and Bound Algorithm for
  computing Pathwidth*, David Coudert, Dorian Mazauric, and Nicolas Nisse. In
  Symposium on Experimental Algorithms (SEA), volume 8504 of LNCS, Copenhagen,
  Denmark, pages 46-58, June 2014,
  http://hal.inria.fr/hal-00943549/document

Authors
-------

- Nathann Cohen (2011-10): Initial version and exact exponential algorithm

- David Coudert (2012-04): MILP formulation and tests functions

- David Coudert (2015-01): BAB formulation and tests functions


Methods
-------
"""
from __future__ import print_function

include "cysignals/memory.pxi"
include "cysignals/signals.pxi"
include 'sage/ext/cdefs.pxi'
from sage.graphs.graph_decompositions.fast_digraph cimport FastDigraph, compute_out_neighborhood_cardinality, popcount32
from libc.stdint cimport uint8_t, int8_t
include "sage/data_structures/binary_matrix.pxi"
from sage.graphs.base.static_dense_graph cimport dense_graph_init

#*****************************************************************************
#          Copyright (C) 2011 Nathann Cohen <nathann.cohen@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

###############
# Lower Bound #
###############

def lower_bound(G):
    r"""
    Returns a lower bound on the vertex separation of `G`

    INPUT:

    - ``G`` -- a Graph or a DiGraph

    OUTPUT:

    A lower bound on the vertex separation of `D` (see the module's
    documentation).

    .. NOTE::

        This method runs in exponential time but has no memory constraint.


    EXAMPLE:

    On a circuit::

        sage: from sage.graphs.graph_decompositions.vertex_separation import lower_bound
        sage: g = digraphs.Circuit(6)
        sage: lower_bound(g)
        1

    TEST:

    Given anything else than a Graph or a DiGraph::

        sage: from sage.graphs.graph_decompositions.vertex_separation import lower_bound
        sage: lower_bound(range(2))
        Traceback (most recent call last):
        ...
        ValueError: The parameter must be a Graph or a DiGraph.

    Given a too large graph::

        sage: from sage.graphs.graph_decompositions.vertex_separation import lower_bound
        sage: lower_bound(graphs.PathGraph(50))
        Traceback (most recent call last):
        ...
        ValueError: The (di)graph can have at most 31 vertices.

    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, Graph) and not isinstance(G, DiGraph):
        raise ValueError("The parameter must be a Graph or a DiGraph.")

    if G.order() >= 32:
        raise ValueError("The (di)graph can have at most 31 vertices.")

    cdef FastDigraph FD = FastDigraph(G)
    cdef int * g = FD.graph
    cdef int n = FD.n

    # minimums[i] is means to store the value of c'_{i+1}
    minimums = <uint8_t *> sig_malloc(sizeof(uint8_t)* n)
    cdef unsigned int i

    # They are initialized to n
    for 0<= i< n:
        minimums[i] = n

    cdef uint8_t tmp, tmp_count

    # We go through all sets
    for 1<= i< <unsigned int> (1<<n):
        tmp_count = <uint8_t> popcount32(i)
        tmp = <uint8_t> compute_out_neighborhood_cardinality(FD, i)

        # And update the costs
        minimums[tmp_count-1] = minimum(minimums[tmp_count-1], tmp)

    # We compute the maximum of all those values
    for 1<= i< n:
        minimums[0] = maximum(minimums[0], minimums[i])

    cdef int min = minimums[0]

    sig_free(minimums)

    return min

##################################################################
# Front end methods for path decomposition and vertex separation #
##################################################################

def path_decomposition(G, algorithm = "BAB", cut_off=None, upper_bound=None, verbose = False,
                       max_prefix_length=20, max_prefix_number=10**6):
    r"""
    Returns the pathwidth of the given graph and the ordering of the vertices
    resulting in a corresponding path decomposition.

    INPUT:

    - ``G`` -- a Graph

    - ``algorithm`` -- (default: ``"BAB"``) Specify the algorithm to use among

      - ``"BAB"`` -- Use a branch-and-bound algorithm. This algorithm has no
        size restriction but could take a very long time on large graphs. It can
        also be used to test is the input (di)graph has vertex separation at
        most ``upper_bound`` or to return the first found solution with vertex
        separation less or equal to a ``cut_off`` value.

      - ``exponential`` -- Use an exponential time and space algorithm. This
        algorithm only works of graphs on less than 32 vertices.

      - ``MILP`` -- Use a mixed integer linear programming formulation. This
        algorithm has no size restriction but could take a very long time.

    - ``upper_bound`` -- (default: ``None``) This is parameter is used by the
      ``"BAB"`` algorithm. If specified, the algorithm searches for a solution
      with ``width < upper_bound``. It helps cutting branches.  However, if the
      given upper bound is too low, the algorithm may not be able to find a
      solution.

    - ``cut_off`` -- (default: None) This is parameter is used by the ``"BAB"``
      algorithm. This bound allows us to stop the search as soon as a solution
      with width at most ``cut_off`` is found, if any. If this bound cannot be
      reached, the best solution found is returned, unless a too low
      ``upper_bound`` is given.

    - ``verbose`` (boolean) -- whether to display information on the
      computations.

    - ``max_prefix_length`` -- (default: 20) limits the length of the stored
      prefixes to prevent storing too many prefixes. This parameter is used only
      when ``algorithm=="BAB"``.

    - ``max_prefix_number`` -- (default: 10**6) upper bound on the number of
      stored prefixes used to prevent using too much memory. This parameter is
      used only when ``algorithm=="BAB"``.

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.

    .. SEEALSO::

        * :meth:`Graph.treewidth` -- computes the treewidth of a graph

    EXAMPLE:

    The pathwidth of a cycle is equal to 2::

        sage: from sage.graphs.graph_decompositions.vertex_separation import path_decomposition
        sage: g = graphs.CycleGraph(6)
        sage: pw, L = path_decomposition(g, algorithm = "BAB"); pw
        2
        sage: pw, L = path_decomposition(g, algorithm = "exponential"); pw
        2
        sage: pw, L = path_decomposition(g, algorithm = "MILP"); pw
        2

    TEST:

    Given anything else than a Graph::

        sage: from sage.graphs.graph_decompositions.vertex_separation import path_decomposition
        sage: path_decomposition(DiGraph())
        Traceback (most recent call last):
        ...
        ValueError: The parameter must be a Graph.

    Given a wrong algorithm::

        sage: from sage.graphs.graph_decompositions.vertex_separation import path_decomposition
        sage: path_decomposition(Graph(), algorithm="SuperFast")
        Traceback (most recent call last):
        ...
        ValueError: Algorithm "SuperFast" has not been implemented yet. Please contribute.

    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("The parameter must be a Graph.")

    return vertex_separation(G, algorithm=algorithm, cut_off=cut_off, upper_bound=upper_bound,
                             verbose=verbose, max_prefix_length=max_prefix_length,
                             max_prefix_number=max_prefix_number)


def vertex_separation(G, algorithm = "BAB", cut_off=None, upper_bound=None, verbose = False,
                      max_prefix_length=20, max_prefix_number=10**6):
    r"""
    Returns an optimal ordering of the vertices and its cost for
    vertex-separation.

    INPUT:

    - ``G`` -- a Graph or a DiGraph

    - ``algorithm`` -- (default: ``"BAB"``) Specify the algorithm to use among

      - ``"BAB"`` -- Use a branch-and-bound algorithm. This algorithm has no
        size restriction but could take a very long time on large graphs. It can
        also be used to test is the input (di)graph has vertex separation at
        most ``upper_bound`` or to return the first found solution with vertex
        separation less or equal to a ``cut_off`` value.

      - ``exponential`` -- Use an exponential time and space algorithm. This
        algorithm only works of graphs on less than 32 vertices.

      - ``MILP`` -- Use a mixed integer linear programming formulation. This
        algorithm has no size restriction but could take a very long time.

    - ``upper_bound`` -- (default: ``None``) This is parameter is used by the
      ``"BAB"`` algorithm. If specified, the algorithm searches for a solution
      with ``width < upper_bound``. It helps cutting branches.  However, if the
      given upper bound is too low, the algorithm may not be able to find a
      solution.

    - ``cut_off`` -- (default: None) This is parameter is used by the ``"BAB"``
      algorithm. This bound allows us to stop the search as soon as a solution
      with width at most ``cut_off`` is found, if any. If this bound cannot be
      reached, the best solution found is returned, unless a too low
      ``upper_bound`` is given.

    - ``verbose`` (boolean) -- whether to display information on the
      computations.

    - ``max_prefix_length`` -- (default: 20) limits the length of the stored
      prefixes to prevent storing too many prefixes. This parameter is used only
      when ``algorithm=="BAB"``.

    - ``max_prefix_number`` -- (default: 10**6) upper bound on the number of
      stored prefixes used to prevent using too much memory. This parameter is
      used only when ``algorithm=="BAB"``.

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.

    EXAMPLES:

    Comparison of methods::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation
        sage: G = digraphs.DeBruijn(2,3)
        sage: vs,L = vertex_separation(G, algorithm="BAB"); vs
        2
        sage: vs,L = vertex_separation(G, algorithm="exponential"); vs
        2
        sage: vs,L = vertex_separation(G, algorithm="MILP"); vs
        2
        sage: G = graphs.Grid2dGraph(3,3)
        sage: vs,L = vertex_separation(G, algorithm="BAB"); vs
        3
        sage: vs,L = vertex_separation(G, algorithm="exponential"); vs
        3
        sage: vs,L = vertex_separation(G, algorithm="MILP"); vs
        3

    Digraphs with multiple strongly connected components::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation
        sage: D = digraphs.Path(8)
        sage: print(vertex_separation(D))
        (0, [7, 6, 5, 4, 3, 2, 1, 0])
        sage: D = DiGraph( random_DAG(30) )
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
        (3, [8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3])

    TESTS:

    Given a wrong algorithm::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation
        sage: vertex_separation(Graph(), algorithm="SuperFast")
        Traceback (most recent call last):
        ...
        ValueError: Algorithm "SuperFast" has not been implemented yet. Please contribute.

    Given anything else than a Graph or a DiGraph::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation
        sage: vertex_separation(range(4))
        Traceback (most recent call last):
        ...
        ValueError: The parameter must be a Graph or a DiGraph.
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
        raise ValueError('The parameter must be a Graph or a DiGraph.')


    if CC:
        # The graph has several (strongly) connected components. We solve the
        # problem on each of them and order partial solutions in the same order
        # than in list CC. The vertex separation is the maximum over all these
        # subgraphs.
        vs, L = 0, []
        for V in CC:

            if len(V)==1:
                # We can directly add this vertex to the solution
                L.extend(V)

            else:
                # We build the (strongly) connected subgraph and do a recursive
                # call to get its vertex separation and corresponding ordering
                H = G.subgraph(V)
                vsH,LH = vertex_separation(H, algorithm      = algorithm,
                                           cut_off           = cut_off,
                                           upper_bound       = upper_bound,
                                           verbose           = verbose,
                                           max_prefix_length = max_prefix_length,
                                           max_prefix_number = max_prefix_number)

                if vsH==-1:
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
        return vertex_separation_exp(G, verbose = verbose)

    elif algorithm == "MILP":
        return vertex_separation_MILP(G, verbosity = (1 if verbose else 0))

    elif algorithm == "BAB":
        return vertex_separation_BAB(G, cut_off=cut_off, upper_bound=upper_bound, verbose=verbose,
                                     max_prefix_length=max_prefix_length, max_prefix_number = max_prefix_number)

    else:
        raise ValueError('Algorithm "{}" has not been implemented yet. Please contribute.'.format(algorithm))


################################
# Exact exponential algorithms #
################################

def vertex_separation_exp(G, verbose = False):
    r"""
    Returns an optimal ordering of the vertices and its cost for
    vertex-separation.

    INPUT:

    - ``G`` -- a Graph or a DiGraph.

    - ``verbose`` (boolean) -- whether to display information on the
      computations.

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.

    .. NOTE::

        Because of its current implementation, this algorithm only works on
        graphs on less than 32 vertices. This can be changed to 54 if necessary,
        but 32 vertices already require 4GB of memory.

    EXAMPLE:

    The vertex separation of a circuit is equal to 1::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation_exp
        sage: g = digraphs.Circuit(6)
        sage: vertex_separation_exp(g)
        (1, [0, 1, 2, 3, 4, 5])

    TEST:

    Given anything else than a Graph or a DiGraph::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation_exp
        sage: vertex_separation_exp(range(3))
        Traceback (most recent call last):
        ...
        ValueError: The parameter must be a Graph or a DiGraph.

    Graphs with non-integer vertices::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation_exp
        sage: D=digraphs.DeBruijn(2,3)
        sage: vertex_separation_exp(D)
        (2, ['000', '001', '100', '010', '101', '011', '110', '111'])

    Given a too large graph::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation_exp
        sage: vertex_separation_exp(graphs.PathGraph(50))
        Traceback (most recent call last):
        ...
        ValueError: The graph should have at most 31 vertices !
    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, Graph) and not isinstance(G, DiGraph):
        raise ValueError("The parameter must be a Graph or a DiGraph.")

    if G.order() >= 32:
        raise ValueError("The graph should have at most 31 vertices !")

    cdef FastDigraph g = FastDigraph(G)

    if verbose:
        print("Memory allocation")
        g.print_adjacency_matrix()

    sig_on()

    cdef unsigned int mem = 1 << g.n
    cdef uint8_t * neighborhoods = <uint8_t *> sig_malloc(mem)

    if neighborhoods == NULL:
        sig_off()
        raise MemoryError("Error allocating memory. I just tried to allocate "+str(mem>>10)+"MB, could that be too much ?")

    memset(neighborhoods, <uint8_t> -1, mem)

    cdef int i,j , k
    for k in range(g.n):
        if verbose:
            print("Looking for a strategy of cost", str(k))

        if exists(g, neighborhoods, 0, k) <= k:
            break

    if verbose:
        print("... Found !")
        print("Now computing the ordering")

    cdef list order = find_order(g, neighborhoods, k)

    sig_free(neighborhoods)
    sig_off()

    return k, list( g.int_to_vertices[i] for i in order )

##############################################################################
# Actual algorithm, breadh-first search and updates of the costs of the sets #
##############################################################################

# Check whether an ordering with the given cost exists, and updates data in the
# neighborhoods array at the same time. See the module's documentation

cdef inline int exists(FastDigraph g, uint8_t * neighborhoods, int current, int cost):

    # If this is true, it means the set has not been evaluated yet
    if neighborhoods[current] == <uint8_t>-1:
        neighborhoods[current] = compute_out_neighborhood_cardinality(g, current)

    # If the cost of this set is too high, there is no point in going further.
    # Same thing if the current set is the whole vertex set.
    if neighborhoods[current] > cost or (current == (1<<g.n)-1):
        return neighborhoods[current]

    # Minimum of the costs of the outneighbors
    cdef int mini = g.n

    cdef int i
    cdef int next_set


    for i in range(g.n):
        if (current >> i)&1:
            continue

        # For each of the out-neighbors next_set of current
        next_set = current | 1<<i

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

# Returns the ordering once we are sure it exists
cdef list find_order(FastDigraph g, uint8_t * neighborhoods, int cost):
    cdef list ordering = []
    cdef int current = 0
    cdef int n = g.n
    cdef int i

    while n:
        # We look for n vertices

        for i in range(g.n):
            if (current >> i)&1:
                continue

            # Find the next set with small cost (we know it exists)
            next_set = current | 1<<i
            if neighborhoods[next_set] <= cost:
                ordering.append(i)
                current = next_set
                break

        # One less to find
        n -= 1

    return ordering

# Min/Max functions

cdef inline int minimum(int a, int b):
    if a<b:
        return a
    else:
        return b

cdef inline int maximum(int a, int b):
    if a>b:
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


    EXAMPLE:

    Path decomposition of a cycle::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = graphs.CycleGraph(6)
        sage: L = [u for u in G.vertices()]
        sage: vertex_separation.is_valid_ordering(G, L)
        True
        sage: vertex_separation.is_valid_ordering(G, [1,2])
        False

    TEST:

    Giving anything else than a Graph or a DiGraph::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: vertex_separation.is_valid_ordering(2, [])
        Traceback (most recent call last):
        ...
        ValueError: The input parameter must be a Graph or a DiGraph.

    Giving anything else than a list::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = graphs.CycleGraph(6)
        sage: vertex_separation.is_valid_ordering(G, {})
        Traceback (most recent call last):
        ...
        ValueError: The second parameter must be of type 'list'.
    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, Graph) and not isinstance(G, DiGraph):
        raise ValueError("The input parameter must be a Graph or a DiGraph.")
    if not isinstance(L, list):
        raise ValueError("The second parameter must be of type 'list'.")

    return set(L) == set(G.vertices())


####################################################################
# Measurement functions of the widths of some graph decompositions #
####################################################################

def width_of_path_decomposition(G, L):
    r"""
    Returns the width of the path decomposition induced by the linear ordering
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
        sage: L = [u for u in G.vertices()]
        sage: vertex_separation.width_of_path_decomposition(G, L)
        2

    Directed path decomposition of a circuit::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = digraphs.Circuit(6)
        sage: L = [u for u in G.vertices()]
        sage: vertex_separation.width_of_path_decomposition(G, L)
        1

    TESTS:

    Path decomposition of a BalancedTree::

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
        ValueError: The input linear vertex ordering L is not valid for G.
    """
    if not is_valid_ordering(G, L):
        raise ValueError("The input linear vertex ordering L is not valid for G.")

    neighbors = G.neighbors_out if G.is_directed() else G.neighbors

    vsL = 0
    S = set()
    neighbors_of_S_in_V_minus_S = set()

    for u in L:

        # We remove u from the neighbors of S
        neighbors_of_S_in_V_minus_S.discard(u)

        # We add vertex u to the set S
        S.add(u)

        # We add the (out-)neighbors of u to the neighbors of S
        for v in neighbors(u):
            if (not v in S):
                neighbors_of_S_in_V_minus_S.add(v)

        # We update the cost of the vertex separation
        vsL = max( vsL, len(neighbors_of_S_in_V_minus_S) )

    return vsL


##########################################
# MILP formulation for vertex separation #
##########################################

def vertex_separation_MILP(G, integrality = False, solver = None, verbosity = 0):
    r"""
    Computes the vertex separation of `G` and the optimal ordering of its
    vertices using an MILP formulation.

    This function uses a mixed integer linear program (MILP) for determining an
    optimal layout for the vertex separation of `G`. This MILP is an improved
    version of the formulation proposed in [SP10]_. See the :mod:`module's
    documentation <sage.graphs.graph_decompositions.vertex_separation>` for more
    details on this MILP formulation.

    INPUT:

    - ``G`` -- a Graph or a DiGraph

    - ``integrality`` -- (default: ``False``) Specify if variables `x_v^t` and
      `u_v^t` must be integral or if they can be relaxed. This has no impact on
      the validity of the solution, but it is sometimes faster to solve the
      problem using binary variables only.

    - ``solver`` -- (default: ``None``) Specify a Linear Program (LP) solver to
      be used. If set to ``None``, the default one is used. For more information
      on LP solvers and which default solver is used, see the method
      :meth:`solve<sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the
      class
      :class:`MixedIntegerLinearProgram<sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- integer (default: ``0``). Sets the level of verbosity. Set
      to 0 by default, which means quiet.

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.

    EXAMPLE:

    Vertex separation of a De Bruijn digraph::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = digraphs.DeBruijn(2,3)
        sage: vs, L = vertex_separation.vertex_separation_MILP(G); vs
        2
        sage: vs == vertex_separation.width_of_path_decomposition(G, L)
        True
        sage: vse, Le = vertex_separation.vertex_separation(G); vse
        2

    The vertex separation of a circuit is 1::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: G = digraphs.Circuit(6)
        sage: vs, L = vertex_separation.vertex_separation_MILP(G); vs
        1

    TESTS:

    Comparison with exponential algorithm::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: for i in range(10):
        ...       G = digraphs.RandomDirectedGNP(10, 0.2)
        ...       ve, le = vertex_separation.vertex_separation(G)
        ...       vm, lm = vertex_separation.vertex_separation_MILP(G)
        ...       if ve != vm:
        ...          print("The solution is not optimal!")

    Comparison with different values of the integrality parameter::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: for i in range(10):  # long time (11s on sage.math, 2012)
        ....:     G = digraphs.RandomDirectedGNP(10, 0.2)
        ....:     va, la = vertex_separation.vertex_separation_MILP(G, integrality=False)
        ....:     vb, lb = vertex_separation.vertex_separation_MILP(G, integrality=True)
        ....:     if va != vb:
        ....:        print("The integrality parameter changes the result!")

    Giving anything else than a Graph or a DiGraph::

        sage: from sage.graphs.graph_decompositions import vertex_separation
        sage: vertex_separation.vertex_separation_MILP([])
        Traceback (most recent call last):
        ...
        ValueError: The first input parameter must be a Graph or a DiGraph.
    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, Graph) and not isinstance(G, DiGraph):
        raise ValueError("The first input parameter must be a Graph or a DiGraph.")

    from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException
    p = MixedIntegerLinearProgram( maximization = False, solver = solver )

    # Declaration of variables.
    x = p.new_variable(binary=integrality, nonnegative=True)
    u = p.new_variable(binary=integrality, nonnegative=True)
    y = p.new_variable(binary=True)
    z = p.new_variable(integer=True, nonnegative=True)

    N = G.num_verts()
    V = G.vertices()
    neighbors_out = G.neighbors_out if G.is_directed() else G.neighbors

    # (2) x[v,t] <= x[v,t+1]   for all v in V, and for t:=0..N-2
    # (3) y[v,t] <= y[v,t+1]   for all v in V, and for t:=0..N-2
    for v in V:
        for t in xrange(N-1):
            p.add_constraint( x[v,t] - x[v,t+1] <= 0 )
            p.add_constraint( y[v,t] - y[v,t+1] <= 0 )

    # (4) y[v,t] <= x[w,t]  for all v in V, for all w in N^+(v), and for all t:=0..N-1
    for v in V:
        for w in neighbors_out(v):
            for t in xrange(N):
                p.add_constraint( y[v,t] - x[w,t] <= 0 )

    # (5) sum_{v in V} y[v,t] == t+1 for t:=0..N-1
    for t in xrange(N):
        p.add_constraint( p.sum([ y[v,t] for v in V ]) == t+1 )

    # (6) u[v,t] >= x[v,t]-y[v,t]    for all v in V, and for all t:=0..N-1
    for v in V:
        for t in xrange(N):
            p.add_constraint( x[v,t] - y[v,t] - u[v,t] <= 0 )

    # (7) z >= sum_{v in V} u[v,t]   for all t:=0..N-1
    for t in xrange(N):
        p.add_constraint( p.sum([ u[v,t] for v in V ]) - z['z'] <= 0 )

    # (8)(9) 0 <= x[v,t] and u[v,t] <= 1
    if not integrality:
        for v in V:
            for t in xrange(N):
                p.add_constraint( 0 <= x[v,t] <= 1 )
                p.add_constraint( 0 <= u[v,t] <= 1 )

    # (10) y[v,t] in {0,1}
    p.set_binary( y )

    # (11) 0 <= z <= |V|
    p.add_constraint( z['z'] <= N )

    #  (1) Minimize z
    p.set_objective( z['z'] )

    try:
        obj = p.solve( log=verbosity )
    except MIPSolverException:
        if integrality:
            raise ValueError("Unbounded or unexpected error")
        else:
            raise ValueError("Unbounded or unexpected error. Try with 'integrality = True'.")

    taby = p.get_values( y )
    tabz = p.get_values( z )
    # since exactly one vertex is processed per step, we can reconstruct the sequence
    seq = []
    for t in xrange(N):
        for v in V:
            if (taby[v,t] > 0) and (not v in seq):
                seq.append(v)
                break
    vs = int(round( tabz['z'] ))

    return vs, seq

##########################################
# Branch and Bound for vertex separation #
##########################################

def vertex_separation_BAB(G,
                          cut_off               = None,
                          upper_bound           = None,
                          max_prefix_length     = 20,
                          max_prefix_number     = 10**6,
                          verbose               = False):
    r"""
    Branch and Bound algorithm for the vertex separation.

    This method implements the branch and bound algorithm for the vertex
    separation of directed graphs and the pathwidth of undirected graphs
    proposed in [CMN14]_. The implementation is valid for both Graph and
    DiGraph. See the documentation of the
    :mod:`~sage.graphs.graph_decompositions.vertex_separation` module.

    INPUT:

    - ``G`` -- a Graph or a DiGraph.

    - ``cut_off`` -- (default: None) bound to consider in the branch and  bound
      algorithm. This allows us to stop the search as soon as a solution with
      width at most ``cut_off`` is found, if any. If this bound cannot be
      reached, the best solution found is returned, unless a too low
      ``upper_bound`` is given.

    - ``upper_bound`` -- (default: None) if specified, the algorithm searches
      for a solution with ``width < upper_bound``. It helps cutting branches.
      However, if the given upper bound is too low, the algorithm may not be
      able to find a solution.

    - ``max_prefix_length`` -- (default: 20) limits the length of the stored
      prefixes to prevent storing too many prefixes.

    - ``max_prefix_number`` -- (default: 10**6) upper bound on the number of
      stored prefixes used to prevent using too much memory.

    - ``verbose`` -- (default: False) display some information when set to True.

    OUTPUT:

    - ``width`` -- the computed vertex separation

    - ``seq`` -- an ordering of the vertices of width ``width``.


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
        sage: vs, seq = VS.vertex_separation_BAB(G, cut_off=11); vs
        11
        sage: vs, seq = VS.vertex_separation_BAB(G, cut_off=10); vs
        10
        sage: vs, seq = VS.vertex_separation_BAB(G, cut_off=9); vs
        10

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
        ValueError: The input parameter must be a Graph or a DiGraph.

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
        ValueError: The input upper bound must be at least 1.
    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, DiGraph) and not isinstance(G, Graph):
        raise ValueError("The input parameter must be a Graph or a DiGraph.")

    cdef int n = G.order()
    if n==0:
        return 0, []

    cut_off = 0 if cut_off is None else cut_off
    upper_bound = n if upper_bound is None else upper_bound
    if upper_bound < 1:
        raise ValueError("The input upper bound must be at least 1.")

    # ==> Allocate and initialize some data structures

    # We use a binary matrix to store the (di)graph. This way the neighborhoud
    # of a vertex is stored in one bitset.
    cdef binary_matrix_t H
    cdef dict vertex_to_int = dense_graph_init(H, G, translation = True)
    cdef int i
    cdef dict int_to_vertex = dict((i, v) for v,i in vertex_to_int.iteritems())

    # We need 2 bitsets here + 3 per call to vertex_separation_BAB_C, so overall
    # 3*n + 2. We use another binary matrix as a pool of bitsets.
    cdef binary_matrix_t bm_pool
    binary_matrix_init(bm_pool, 3*n+2, n)

    cdef int * prefix    = <int *>sig_malloc(n * sizeof(int))
    cdef int * positions = <int *>sig_malloc(n * sizeof(int))
    if prefix==NULL or positions==NULL:
        sig_free(prefix)
        sig_free(positions)
        binary_matrix_free(H)
        binary_matrix_free(bm_pool)
        raise MemoryError("Unable to allocate data strutures.")

    cdef list best_seq = range(n)
    for i in range(n):
        prefix[i] = i
        positions[i] = i

    cdef int width = upper_bound
    cdef list order = list()
    cdef set prefix_storage = set()

    try:
        # ==> Call the cython method
        sig_on()
        width = vertex_separation_BAB_C(H                         = H,
                                        n                         = n,
                                        prefix                    = prefix,
                                        positions                 = positions,
                                        best_seq                  = best_seq,
                                        level                     = 0,
                                        b_prefix                  = bm_pool.rows[3*n],
                                        b_prefix_and_neighborhood = bm_pool.rows[3*n+1],
                                        cut_off                   = cut_off,
                                        upper_bound               = upper_bound,
                                        current_cost              = 0,
                                        bm_pool                   = bm_pool,
                                        prefix_storage            = prefix_storage,
                                        max_prefix_length         = max_prefix_length,
                                        max_prefix_number         = max_prefix_number,
                                        verbose                   = verbose)

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

    return (width if width<upper_bound else -1), order

cdef inline _my_invert_positions(int *prefix, int *positions, int pos_a, int pos_b):
    """
    Permute vertices at positions ``pos_a`` and ``pos_b`` in array ``prefix``,
    and record the new positions in array ``positions``.
    """
    if pos_a!=pos_b:
        positions[prefix[pos_a]],positions[prefix[pos_b]] = positions[prefix[pos_b]],positions[prefix[pos_a]]
        prefix[pos_a], prefix[pos_b] = prefix[pos_b], prefix[pos_a]


cdef int vertex_separation_BAB_C(binary_matrix_t H,
                                 int             n,
                                 int *           prefix,
                                 int *           positions,
                                 list            best_seq,
                                 int             level,
                                 bitset_t        b_prefix,
                                 bitset_t        b_prefix_and_neighborhood,
                                 int             cut_off,
                                 int             upper_bound,
                                 int             current_cost,
                                 binary_matrix_t bm_pool,
                                 set             prefix_storage,
                                 int             max_prefix_length,
                                 int             max_prefix_number,
                                 bint            verbose):
    r"""
    Branch and Bound algorithm for the process number and the vertex separation.

    INPUT:

    - ``H`` -- a binary matrix storing the adjacency of the (di)graph

    - ``n`` -- the number of vertices of the (di)graph

    - ``prefix`` -- array of ``n`` integers containing a permutation of the
      vertices. The vertices forming the current prefix under consideration are
      stored in cells ``[0,level-1]``.

    - ``positions`` -- array of ``n`` integers associating to each vertex its
      index in array ``prefix``.

    - ``best_seq`` -- array of ``n`` integers storing the best ordering found so
      far.

    - ``level`` -- an integer specifying the length of the current prefix.

    - ``b_prefix`` -- a bitset of size ``n`` recording the vertices in the
      current prefix (in cells ``[0,level-1]``).

    - ``b_prefix_and_neighborhood`` -- a bitset of size ``n`` recording the
      vertices in the current prefix and the vertices in its neighborhood.

    - ``cut_off`` -- (default: None) bound to consider in the branch and  bound
      algorithm. This allows us to stop the search as soon as a solution with
      width at most ``cut_off`` is found, if any.

    - ``upper_bound`` -- the algorithm searches for a solution with ``width <
      upper_bound``. It helps cutting branches. Each time a new solution is
      found, the upper bound is reduced.

    - ``bm_pool`` -- a binary matrix used with ``3*n+2`` rows of size
      ``n``. Each rows is a bitset of size ``n``. This data structure is used as
      a pool of initialized bitsets. Each call of this method needs 3 bitsets
      for local operations, so it uses rows ``[3*level,3*level+2]``.

    - ``prefix_storage`` -- set used to store prefixes.

    - ``max_prefix_length`` -- maximum length of the stored prefixes to prevent
      storing too many prefixes.

    - ``max_prefix_number`` -- upper bound on the number of stored prefixes used
      to prevent using too much memory.

    - ``verbose`` -- (default: False) display some information when set to True.
    """
    cdef int i

    # ==> Test termination

    if level==n:
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

    cdef bitset_s *loc_b_prefix         = bm_pool.rows[3*level]
    cdef bitset_s *loc_b_pref_and_neigh = bm_pool.rows[3*level+1]
    cdef bitset_s *b_tmp                = bm_pool.rows[3*level+2]
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
    while i<n:

        j = prefix[i]

        if bitset_issubset(H.rows[j], loc_b_pref_and_neigh):
            # (i) Vertex j is such that all its out-neighbors are in the prefix
            # or in its out-neighborhood (so in loc_b_pref_and_neigh).
            bitset_add(loc_b_pref_and_neigh, j)
            select_it = 1

        elif bitset_in(loc_b_pref_and_neigh, j) and not bitset_in(loc_b_prefix, j):
            bitset_difference(b_tmp, H.rows[j], loc_b_pref_and_neigh)
            if bitset_len(b_tmp)==1:
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
    if loc_level==n:
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

    if loc_level<=max_prefix_length:
        frozen_prefix = frozenset([prefix[i] for i in range(loc_level)])
        if frozen_prefix in prefix_storage:
            return upper_bound


    # ==> Sort and Prune
    #
    # We compute for each remaining vertex v a lower bound on the width of any
    # ordering with prefix prefix+v
    for i from loc_level <= i < n:
        j = prefix[i]
        bitset_union(b_tmp, loc_b_pref_and_neigh, H.rows[j])
        bitset_difference(b_tmp, b_tmp, loc_b_prefix)
        bitset_discard(b_tmp, j)
        delta_i = bitset_len(b_tmp)
        if delta_i < upper_bound:
            delta.append( (delta_i, j) )

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

        cost_i = vertex_separation_BAB_C(H                         = H,
                                         n                         = n,
                                         prefix                    = prefix,
                                         positions                 = positions,
                                         best_seq                  = best_seq,
                                         level                     = loc_level+1,
                                         b_prefix                  = loc_b_prefix,
                                         b_prefix_and_neighborhood = b_tmp,
                                         cut_off                   = cut_off,
                                         upper_bound               = upper_bound,
                                         current_cost              = delta_i,
                                         bm_pool                   = bm_pool,
                                         prefix_storage            = prefix_storage,
                                         max_prefix_length         = max_prefix_length,
                                         max_prefix_number         = max_prefix_number,
                                         verbose                   = verbose)

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
    if loc_level<=max_prefix_length and current_cost<upper_bound and len(prefix_storage)<max_prefix_number:
        prefix_storage.add(frozen_prefix)

    return upper_bound
