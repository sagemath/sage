# cython: binding=True
# distutils: language = c++

"""
Graph coloring

This module gathers all methods related to graph coloring. Here is what it can
do :

**Proper vertex coloring**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`all_graph_colorings` | Compute all `n`-colorings a graph
    :meth:`first_coloring` | Return the first vertex coloring found
    :meth:`number_of_n_colorings` | Compute the number of `n`-colorings of a graph
    :meth:`numbers_of_colorings` | Compute the number of colorings of a graph
    :meth:`chromatic_number` | Return the chromatic number of the graph
    :meth:`vertex_coloring` | Compute vertex colorings and chromatic numbers

**Fractional relaxations**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`fractional_chromatic_number` | Return the fractional chromatic number of the graph
    :meth:`fractional_chromatic_index` | Return the fractional chromatic index of the graph

**Other colorings**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`grundy_coloring` | Compute Grundy numbers and Grundy colorings
    :meth:`b_coloring` | Compute b-chromatic numbers and b-colorings
    :meth:`edge_coloring` | Compute chromatic index and edge colorings
    :meth:`round_robin` | Compute a round-robin coloring of the complete graph on `n` vertices
    :meth:`linear_arboricity` | Compute the linear arboricity of the given graph
    :meth:`acyclic_edge_coloring` | Compute an acyclic edge coloring of the current graph



AUTHORS:

- Tom Boothby (2008-02-21): Initial version
- Carlo Hamalainen (2009-03-28): minor change: switch to C++ DLX solver
- Nathann Cohen (2009-10-24): Coloring methods using linear programming

Methods
-------
"""

# ****************************************************************************
#           Copyright (C) 2008 Tom Boothby <boothby@u.washington.edu>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy
from sage.combinat.matrices.dlxcpp import DLXCPP
from sage.plot.colors import rainbow
from libcpp.vector cimport vector
from libcpp.pair cimport pair

from sage.numerical.mip import MixedIntegerLinearProgram
from sage.numerical.mip import MIPSolverException
from sage.graphs.independent_sets import IndependentSets


def all_graph_colorings(G, n, count_only=False, hex_colors=False, vertex_color_dict=False):
    r"""
    Compute all `n`-colorings of a graph.

    This method casts the graph coloring problem into an exact cover problem,
    and passes this into an implementation of the Dancing Links algorithm
    described by Knuth (who attributes the idea to Hitotumatu and Noshita).

    INPUT:

    * ``G`` -- a graph

    * ``n`` -- a positive integer; the number of colors

    * ``count_only`` -- boolean (default: ``False``); when set to ``True``, it
      returns 1 for each coloring

    * ``hex_colors`` -- boolean (default: ``False``); when set to ``False``,
      colors are labeled [0, 1, ..., `n - 1`], otherwise the RGB Hex labeling
      is used

    * ``vertex_color_dict`` -- boolean (default: ``False``); when set to
      ``True``, it returns a dictionary ``{vertex: color}``, otherwise it
      returns a dictionary ``{color: [list of vertices]}``

    The construction works as follows. Columns:

    * The first `|V|` columns correspond to a vertex -- a `1` in this column
      indicates that this vertex has a color.

    * After those `|V|` columns, we add `n*|E|` columns -- a `1` in these
      columns indicate that a particular edge is incident to a vertex with a
      certain color.

    Rows:

    * For each vertex, add `n` rows; one for each color `c`. Place a `1` in the
      column corresponding to the vertex, and a `1` in the appropriate column
      for each edge incident to the vertex, indicating that that edge is
      incident to the color `c`.

    * If `n > 2`, the above construction cannot be exactly covered since each
      edge will be incident to only two vertices (and hence two colors) - so we
      add `n*|E|` rows, each one containing a `1` for each of the `n*|E|`
      columns. These get added to the cover solutions "for free" during the
      backtracking.

    Note that this construction results in `n*|V| + 2*n*|E| + n*|E|` entries in
    the matrix.  The Dancing Links algorithm uses a sparse representation, so if
    the graph is simple, `|E| \leq |V|^2` and `n <= |V|`, this construction runs
    in `O(|V|^3)` time.  Back-conversion to a coloring solution is a simple scan
    of the solutions, which will contain `|V| + (n-2)*|E|` entries, so runs in
    `O(|V|^3)` time also. For most graphs, the conversion will be much faster
    -- for example, a planar graph will be transformed for `4`-coloring in
    linear time since `|E| = O(|V|)`.

    REFERENCES:

    http://www-cs-staff.stanford.edu/~uno/papers/dancing-color.ps.gz

    EXAMPLES::

        sage: from sage.graphs.graph_coloring import all_graph_colorings
        sage: G = Graph({0: [1, 2, 3], 1: [2]})
        sage: n = 0
        sage: for C in all_graph_colorings(G, 3, hex_colors=True):
        ....:     parts = [C[k] for k in C]
        ....:     for P in parts:
        ....:         l = len(P)
        ....:         for i in range(l):
        ....:             for j in range(i + 1, l):
        ....:                 if G.has_edge(P[i], P[j]):
        ....:                     raise RuntimeError("Coloring Failed.")
        ....:     n+=1
        sage: print("G has %s 3-colorings." % n)
        G has 12 3-colorings.


    TESTS::

        sage: G = Graph({0: [1, 2, 3], 1: [2]})
        sage: for C in all_graph_colorings(G, 0):
        ....:     print(C)
        sage: for C in all_graph_colorings(G, -1):
        ....:     print(C)
        Traceback (most recent call last):
        ...
        ValueError: n must be non-negative
        sage: G = Graph({0: [1], 1: [2]})
        sage: for c in all_graph_colorings(G, 2, vertex_color_dict=True):
        ....:     print(c)
        {0: 0, 1: 1, 2: 0}
        {0: 1, 1: 0, 2: 1}
        sage: for c in all_graph_colorings(G, 2, hex_colors=True):
        ....:     print(sorted(c.items()))
        [('#00ffff', [1]), ('#ff0000', [0, 2])]
        [('#00ffff', [0, 2]), ('#ff0000', [1])]
        sage: for c in all_graph_colorings(G, 2, hex_colors=True, vertex_color_dict=True):
        ....:     print(c)
        {0: '#ff0000', 1: '#00ffff', 2: '#ff0000'}
        {0: '#00ffff', 1: '#ff0000', 2: '#00ffff'}
        sage: for c in all_graph_colorings(G, 2, vertex_color_dict=True):
        ....:     print(c)
        {0: 0, 1: 1, 2: 0}
        {0: 1, 1: 0, 2: 1}
        sage: for c in all_graph_colorings(G, 2, count_only=True, vertex_color_dict=True):
        ....:     print(c)
        1
        1
    """
    G._scream_if_not_simple(allow_multiple_edges=True)

    if not n:
        return
    if n < 0:
        raise ValueError("n must be non-negative")

    cdef list V = list(G)

    cdef int nV = G.order()
    cdef int nE = G.size()

    cdef vector[pair[int, vector[int]]] ones
    cdef dict Vd = {}
    cdef dict colormap = {}
    cdef int k = 0
    for i, v in enumerate(V):
        Vd[v] = i
        for c in range(n):
            ones.push_back((k, [i]))
            colormap[k] = (v, c)
            k += 1

    cdef int kk = nV
    cdef int v0, v1
    for e in G.edges(labels=False, sort=False):
        v0 = n * Vd[e[0]]
        v1 = n * Vd[e[1]]
        for c in range(n):
            ones[v0].second.push_back(kk + c)
            ones[v1].second.push_back(kk + c)
            v0 += 1
            v1 += 1
        kk += n

    if n > 2:
        for i in range(n * nE):
            ones.push_back((k + i, [nV + i]))

    cdef list colors = rainbow(n)
    cdef dict color_dict = {col: i for i, col in enumerate(colors)}

    cdef list ones_second = [ones[i].second for i in range(len(ones))]
    cdef dict coloring = {}

    try:
        for a in DLXCPP(ones_second):
            if count_only:
                yield 1
                continue
            coloring = {}
            if vertex_color_dict:
                for x in a:
                    if x in colormap:
                        v,c = colormap[x]
                        if hex_colors:
                            coloring[v] = colors[c]
                        else:
                            coloring[v] = color_dict[colors[c]]
            else:
                for x in a:
                    if x in colormap:
                        v,c = colormap[x]
                        if hex_colors:
                            if colors[c] in coloring:
                                coloring[colors[c]].append(v)
                            else:
                                coloring[colors[c]] = [v]
                        else:
                            if color_dict[colors[c]] in coloring:
                                coloring[color_dict[colors[c]]].append(v)
                            else:
                                coloring[color_dict[colors[c]]] = [v]
            yield coloring
    except RuntimeError:
        raise RuntimeError("too much recursion, Graph coloring failed")

cpdef first_coloring(G, n=0, hex_colors=False):
    r"""
    Return the first vertex coloring found.

    If a natural number `n` is provided, returns the first found coloring with
    at least `n` colors.

    INPUT:

    -  ``n`` -- integer (default: 0); the minimal number of colors to try

    - ``hex_colors`` -- boolean (default: ``False``); when set to ``True``, the
      partition returned is a dictionary whose keys are colors and whose values
      are the color classes (ideal for plotting)

    EXAMPLES::

        sage: from sage.graphs.graph_coloring import first_coloring
        sage: G = Graph({0: [1, 2, 3], 1: [2]})
        sage: sorted(first_coloring(G, 3))
        [[0], [1, 3], [2]]
    """
    G._scream_if_not_simple(allow_multiple_edges=True)
    cdef int o = G.order()
    for m in range(n, o + 1):
        for C in all_graph_colorings(G, m, hex_colors=True):
            if hex_colors:
                return C
            else:
                return list(C.values())

cpdef number_of_n_colorings(G, n):
    r"""
    Compute the number of `n`-colorings of a graph

    INPUT:

    - ``G`` -- a graph

    - ``n`` -- a positive integer; the number of colors

    EXAMPLES::

        sage: from sage.graphs.graph_coloring import number_of_n_colorings
        sage: G = Graph({0: [1, 2, 3], 1: [2]})
        sage: number_of_n_colorings(G, 3)
        12
    """
    # Take care of the stupid stuff
    if not n:
        return int(G.order() == 0)
    if n == 1:
        return int(G.size() == 0)
    if n < 0:
        # negative colors?? what does that even mean?
        return 0

    m = 0
    for C in all_graph_colorings(G, n, count_only=True):
        m += 1
    return m

cpdef numbers_of_colorings(G):
    r"""
    Compute the number of colorings of a graph.

    Return the number of `n`-colorings of the graph `G` for all `n` from `0` to
    `|V|`.

    EXAMPLES::

        sage: from sage.graphs.graph_coloring import numbers_of_colorings
        sage: G = Graph({0: [1, 2, 3], 1: [2]})
        sage: numbers_of_colorings(G)
        [0, 0, 0, 12, 72]
    """
    cdef int o = G.order()
    cdef list answer = [number_of_n_colorings(G, n) for n in range(o + 1)]
    return answer

cpdef chromatic_number(G):
    r"""
    Return the chromatic number of the graph.

    The chromatic number is the minimal number of colors needed to color the
    vertices of the graph `G`.

    EXAMPLES::

        sage: from sage.graphs.graph_coloring import chromatic_number
        sage: G = Graph({0: [1, 2, 3], 1: [2]})
        sage: chromatic_number(G)
        3

        sage: G = graphs.PetersenGraph()
        sage: G.chromatic_number()
        3
    """
    G._scream_if_not_simple(allow_multiple_edges=True)
    cdef int o = G.order()
    cdef int m
    if not o:
        return 0
    if not G.size():
        return 1
    elif G.is_bipartite(): # can we do it in linear time?
        return 2
    else: # counting cliques is faster than our brute-force method...
        m = G.clique_number()
    if m >= o - 1: # marginal improvement... if there's an o-1 clique and not an o clique, don't waste our time coloring.
        return m
    for n in range(m, o + 1):
        for C in all_graph_colorings(G, n):
            return n


def vertex_coloring(g, k=None, value_only=False, hex_colors=False, solver=None, verbose=0):
    r"""
    Compute Vertex colorings and chromatic numbers.

    This function can compute the chromatic number of the given graph or test
    its `k`-colorability.

    See the :wikipedia:`Graph_coloring` for further details on graph coloring.

    INPUT:

    - ``g`` -- a graph.

    - ``k`` -- integer (default: ``None``); tests whether the graph is
      `k`-colorable.  The function returns a partition of the vertex set in `k`
      independent sets if possible and ``False`` otherwise.

    - ``value_only`` -- boolean (default: ``False``):

      - When set to ``True``, only the chromatic number is returned.

      - When set to ``False`` (default), a partition of the vertex set into
        independent sets is returned if possible.

    - ``hex_colors`` -- boolean (default: ``False``); when set to ``True``, the
      partition returned is a dictionary whose keys are colors and whose values
      are the color classes (ideal for plotting).

    - ``solver`` -- (default: ``None``); specify a Linear Program (LP) solver to
      be used. If set to ``None``, the default one is used. For more information
      on LP solvers and which default solver is used, see the method
      :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the
      class :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- integer (default: ``0``); sets the level of verbosity. Set
      to 0 by default, which means quiet.


    OUTPUT:

    - If ``k=None`` and ``value_only=False``, then return a partition of the
      vertex set into the minimum possible of independent sets.

    - If ``k=None`` and ``value_only=True``, return the chromatic number.

    - If ``k`` is set and ``value_only=None``, return ``False`` if the graph is
      not `k`-colorable, and a partition of the vertex set into `k` independent
      sets otherwise.

    - If ``k`` is set and ``value_only=True``, test whether the graph is
      `k`-colorable, and return ``True`` or ``False`` accordingly.

    EXAMPLES::

       sage: from sage.graphs.graph_coloring import vertex_coloring
       sage: g = graphs.PetersenGraph()
       sage: vertex_coloring(g, value_only=True)
       3

    TESTS:

    Empty graph::

       sage: from sage.graphs.graph_coloring import vertex_coloring
       sage: empty = Graph()
       sage: vertex_coloring(empty, value_only=True)
       0
       sage: vertex_coloring(empty, hex_colors=True)
       {}
       sage: vertex_coloring(empty)
       []
    """
    g._scream_if_not_simple(allow_multiple_edges=True)
    from sage.plot.colors import rainbow
    cdef list colorings, value
    cdef set vertices
    cdef list deg
    cdef list neighbors
    cdef list classes

    # If k is None, tries to find an optimal coloring
    if k is None:
        # No need to start a linear program if the graph is an
        # independent set, is bipartite, or is empty.
        # - Empty graph
        if not g.order():
            if value_only:
                return 0
            elif hex_colors:
                return dict()
            else:
                return []
        # - Independent set
        if not g.size():
            if value_only:
                return 1
            elif hex_colors:
                return {rainbow(1)[0]: list(g)}
            else:
                return [list(g)]
        # - Bipartite set
        if g.is_bipartite():
            if value_only:
                return 2
            elif hex_colors:
                return dict(zip(rainbow(2), g.bipartite_sets()))
            else:
                return g.bipartite_sets()

        # - No need to try any k smaller than the maximum clique in
        # - the graph No need to try k less than |G|/alpha(G), as each color
        #   class is at most alpha(G)
        # - max, because we know it is not bipartite
        from math import ceil
        k = int(max([3, g.clique_number(), ceil(g.order() / len(g.independent_set()))]))

        while True:
            # tries to color the graph, increasing k each time it fails.
            tmp = vertex_coloring(g, k=k, value_only=value_only,
                                  hex_colors=hex_colors, solver=solver, verbose=verbose)
            if tmp is not False:
                if value_only:
                    return k
                else:
                    return tmp
            k += 1
    else:
        # Is the graph empty?
        # If the graph is empty, something should be returned.
        # This is not so stupid, as the graph could be emptied
        # by the test of degeneracy.
        if not g.order():
            if value_only:
                return True
            elif hex_colors:
                return {color: [] for color in rainbow(k)}
            else:
                return [[] for i in range(k)]
        # Is the graph connected?
        # This is not so stupid, as the graph could be disconnected
        # by the test of degeneracy (as previously).
        if not g.is_connected():
            if value_only:
                for component in g.connected_components():
                    tmp = vertex_coloring(g.subgraph(component), k=k,
                                          value_only=value_only,
                                          hex_colors=hex_colors,
                                          solver=solver, verbose=verbose)
                    if tmp is False:
                        return False
                return True
            colorings = []
            for component in g.connected_components():
                tmp = vertex_coloring(g.subgraph(component), k=k,
                                      value_only=value_only,
                                      hex_colors=False,
                                      solver=solver, verbose=verbose)
                if tmp is False:
                    return False
                colorings.append(tmp)
            value = [[] for color in range(k)]
            for color in range(k):
                for component in colorings:
                    value[color].extend(component[color])
            if hex_colors:
                return dict(zip(rainbow(k), value))
            else:
                return value

        # Degeneracy
        # Vertices whose degree is less than k are of no importance in
        # the coloring.
        if min(g.degree()) < k:
            vertices = set(g)
            deg = []
            tmp = [v for v in vertices if g.degree(v) < k]
            while tmp:
                v = tmp.pop(0)
                neighbors = list(set(g.neighbor_iterator(v)) & vertices)
                if v in vertices and len(neighbors) < k:
                    vertices.remove(v)
                    tmp.extend(neighbors)
                    deg.append(v)
            if value_only:
                return vertex_coloring(g.subgraph(list(vertices)), k=k,
                                       value_only=value_only,
                                       hex_colors=hex_colors,
                                       solver=solver, verbose=verbose)
            value = vertex_coloring(g.subgraph(list(vertices)), k=k,
                                    value_only=value_only,
                                    hex_colors=False,
                                    solver=solver, verbose=verbose)
            if value is False:
                return False
            while deg:
                for classe in value:
                    if set(classe).isdisjoint(g.neighbor_iterator(deg[-1])):
                        classe.append(deg[-1])
                        deg.pop(-1)
                        break
            if hex_colors:
                return dict(zip(rainbow(k), value))
            else:
                return value

        p = MixedIntegerLinearProgram(maximization=True, solver=solver)
        color = p.new_variable(binary=True)

        # a vertex has exactly one color
        for v in g:
            p.add_constraint(p.sum(color[v,i] for i in range(k)), min=1, max=1)

        # adjacent vertices have different colors
        for u, v in g.edge_iterator(labels=None):
            for i in range(k):
                p.add_constraint(color[u,i] + color[v,i], max=1)

        # The first vertex is colored with 1. It costs nothing to say
        # it, and it can help.
        p.add_constraint(color[next(g.vertex_iterator()),0],  max=1, min=1)

        try:
            if value_only:
                p.solve(objective_only=True, log=verbose)
                return True
            else:
                p.solve(log=verbose)
        except MIPSolverException:
            return False

        color = p.get_values(color)
        # builds the color classes
        classes = [[] for i in range(k)]

        for v in g:
            for i in range(k):
                if color[v,i] == 1:
                    classes[i].append(v)
                    break

        if hex_colors:
            return dict(zip(rainbow(len(classes)), classes))
        else:
            return classes

# Fractional relaxations
def fractional_chromatic_number(G, solver='PPL', verbose=0,
                                check_components=True, check_bipartite=True):
    r"""
    Return the fractional chromatic number of the graph.

    Fractional coloring is a relaxed version of vertex coloring with several
    equivalent definitions, such as the optimum value in a linear relaxation of
    the integer program that gives the usual chromatic number. It is also equal
    to the fractional clique number by LP-duality.

    ALGORITHM:

    The fractional chromatic number is computed via the usual Linear Program.
    The LP solved by sage is essentially,

    .. MATH::

        \mbox{Minimize : }&\sum_{I\in \mathcal{I}(G)} x_{I}\\
        \mbox{Such that : }&\\
        &\forall v\in V(G), \sum_{I\in \mathcal{I}(G),\, v\in I}x_{v}\geq 1\\
        &\forall I\in \mathcal{I}(G), x_{I} \geq 0

    where `\mathcal{I}(G)` is the set of maximal independent sets of `G` (see
    Section 2.1 of [CFKPR2010]_ to know why it is sufficient to consider maximal
    independent sets). As optional optimisations, we construct the LP on each
    biconnected component of `G` (and output the maximum value), and avoid using
    the LP if G is bipartite (as then the output must be 1 or 2).

    .. NOTE::

        Computing the fractional chromatic number can be very slow. Since the
        variables of the LP are independent sets, in general the LP has size
        exponential in the order of the graph. In the current implementation a
        list of all maximal independent sets is created and stored, which can be
        both slow and memory-hungry.

    INPUT:

    - ``G`` -- a graph

    - ``solver`` -- (default: ``"PPL"``); specify a Linear Program (LP) solver
      to be used. If set to ``None``, the default one is used. For more
      information on LP solvers and which default solver is used, see the method
      :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the
      class :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

      .. NOTE::

          The default solver used here is ``"PPL"`` which provides exact
          results, i.e. a rational number, although this may be slower that
          using other solvers.

    - ``verbose`` -- integer (default: `0`); sets the level of verbosity of
      the LP solver

    - ``check_components`` -- boolean (default: ``True``); whether the method is
      called on each biconnected component of `G`

    - ``check_bipartite`` -- boolean (default: ``True``); whether the graph is
      checked for bipartiteness. If the graph is bipartite then we can avoid
      creating and solving the LP.

    EXAMPLES:

    The fractional chromatic number of a `C_5` is `5/2`::

        sage: g = graphs.CycleGraph(5)
        sage: g.fractional_chromatic_number()
        5/2

    TESTS::

        sage: G = graphs.RandomGNP(20, .2)
        sage: a = G.fractional_chromatic_number(check_components=True)
        sage: b = G.fractional_chromatic_number(check_components=False)
        sage: a == b
        True
    """
    G._scream_if_not_simple()

    if not G.order():
        return 0
    if not G.size():
        # The fractional chromatic number of an independent set is 1
        return 1
    if check_bipartite and G.is_bipartite():
        # at this point we've already ascertained g.size() > 0
        # so g is a bipartite graph with at least one edge
        return 2
    if check_components:
        return max(fractional_chromatic_number(G.subgraph(b), solver=solver,
                                               verbose=verbose,
                                               check_components=False,
                                               check_bipartite=check_bipartite)
                   for b in G.blocks_and_cut_vertices()[0])

    Is = [frozenset(I) for I in IndependentSets(G, maximal=True)]

    # Initialize LP for fractional chromatic number, we want to minimize the
    # total weight
    p = MixedIntegerLinearProgram(solver=solver, maximization=False)

    # One nonnegative variable per maximal independent set
    w = p.new_variable(nonnegative=True)

    # the objective is the sum of weights of the independent sets
    p.set_objective(p.sum(w[I] for I in Is))

    # such that each vertex gets total weight at least 1
    for v in G:
        p.add_constraint(p.sum(w[I] for I in Is if v in I), min=1)

    obj = p.solve(log=verbose)

    return obj

def fractional_chromatic_index(G, solver="PPL", verbose_constraints=False, verbose=0):
    r"""
    Return the fractional chromatic index of the graph.

    The fractional chromatic index is a relaxed version of edge-coloring. An
    edge coloring of a graph being actually a covering of its edges into the
    smallest possible number of matchings, the fractional chromatic index of a
    graph `G` is the smallest real value `\chi_f(G)` such that there exists a
    list of matchings `M_1, \ldots, M_k` of `G` and coefficients `\alpha_1,
    \ldots, \alpha_k` with the property that each edge is covered by the
    matchings in the following relaxed way

    .. MATH::

        \forall e \in E(G), \sum_{e \in M_i} \alpha_i \geq 1.

    For more information, see the :wikipedia:`Fractional_coloring`.

    ALGORITHM:

    The fractional chromatic index is computed through Linear Programming
    through its dual. The LP solved by sage is actually:

    .. MATH::

        \mbox{Maximize : }&\sum_{e\in E(G)} r_{e}\\
        \mbox{Such that : }&\\
        &\forall M\text{ matching }\subseteq G, \sum_{e\in M}r_{v}\leq 1\\

    INPUT:

    - ``G`` -- a graph

    - ``solver`` -- (default: ``"PPL"``); specify a Linear Program (LP) solver
      to be used. If set to ``None``, the default one is used. For more
      information on LP solvers and which default solver is used, see the method
      :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the
      class :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

      .. NOTE::

          The default solver used here is ``"PPL"`` which provides exact
          results, i.e. a rational number, although this may be slower that
          using other solvers. Be aware that this method may loop endlessly when
          using some non exact solvers as reported in :trac:`23658` and
          :trac:`23798`.

    - ``verbose_constraints`` -- boolean (default: ``False``); whether to
      display which constraints are being generated

    - ``verbose`` -- integer (default: `0`); sets the level of verbosity of the
      LP solver

    EXAMPLES:

    The fractional chromatic index of a `C_5` is `5/2`::

        sage: g = graphs.CycleGraph(5)
        sage: g.fractional_chromatic_index()
        5/2

    TESTS:

    Issue reported in :trac:`23658` and :trac:`23798` with non exact
    solvers::

        sage: g = graphs.PetersenGraph()
        sage: g.fractional_chromatic_index(solver='GLPK')  # known bug (#23798)
        3.0
        sage: g.fractional_chromatic_index(solver='PPL')
        3
    """
    G._scream_if_not_simple()

    if not G.order():
        return 0
    if not G.size():
        return 1

    frozen_edges = [frozenset(e) for e in G.edges(labels=False, sort=False)]

    # Initialize LP for maximum weight matching
    M = MixedIntegerLinearProgram(solver=solver, constraint_generation=True)

    # One variable per edge
    b = M.new_variable(binary=True, nonnegative=True)

    # We want to select at most one incident edge per vertex (matching)
    for u in G:
        M.add_constraint(M.sum(b[frozenset(e)] for e in G.edges_incident(u, labels=False)) <= 1)

    #
    # Initialize LP for fractional chromatic number
    p = MixedIntegerLinearProgram(solver=solver, constraint_generation=True)

    # One variable per edge
    r = p.new_variable(nonnegative=True)

    # We want to maximize the sum of weights on the edges
    p.set_objective(p.sum(r[fe] for fe in frozen_edges))

    # Each edge being by itself a matching, its weight cannot be more than 1
    for fe in frozen_edges:
        p.add_constraint(r[fe] <= 1)

    obj = p.solve(log=verbose)

    while True:

        # Update the weights of edges for the matching problem
        M.set_objective(M.sum(p.get_values(r[fe]) * b[fe] for fe in frozen_edges))

        # If the maximum matching has weight at most 1, we are done !
        if M.solve(log=verbose) <= 1:
            break

        # Otherwise, we add a new constraint
        matching = [fe for fe in frozen_edges if M.get_values(b[fe]) == 1]
        p.add_constraint(p.sum(r[fe] for fe in matching) <= 1)
        if verbose_constraints:
            print("Adding a constraint on matching : {}".format(matching))

        # And solve again
        obj = p.solve(log=verbose)

    # Accomplished !
    return obj

def grundy_coloring(g, k, value_only=True, solver=None, verbose=0):
    r"""
    Compute Grundy numbers and Grundy colorings.

    The method computes the worst-case of a first-fit coloring with less than
    `k` colors.

    Definition:

    A first-fit coloring is obtained by sequentially coloring the vertices of a
    graph, assigning them the smallest color not already assigned to one of its
    neighbors. The result is clearly a proper coloring, which usually requires
    much more colors than an optimal vertex coloring of the graph, and heavily
    depends on the ordering of the vertices.

    The number of colors required by the worst-case application of this
    algorithm on a graph `G` is called the Grundy number, written `\Gamma (G)`.

    Equivalent formulation:

    Equivalently, a Grundy coloring is a proper vertex coloring such that any
    vertex colored with `i` has, for every `j < i`, a neighbor colored with
    `j`. This can define a Linear Program, which is used here to compute the
    Grundy number of a graph.

    .. NOTE::

       This method computes a grundy coloring using at *MOST* `k` colors. If
       this method returns a value equal to `k`, it cannot be assumed that `k`
       is equal to `\Gamma(G)`. Meanwhile, if it returns any value `k' < k`,
       this is a certificate that the Grundy number of the given graph is `k'`.

       As `\Gamma(G)\leq \Delta(G)+1`, it can also be assumed that `\Gamma(G) =
       k` if ``grundy_coloring(g, k)`` returns `k` when `k = \Delta(G) +1`.

    INPUT:

    - ``k`` -- integer;  maximum number of colors

    - ``value_only`` -- boolean (default: ``True``); when set to ``True``, only
      the number of colors is returned. Otherwise, the pair ``(nb_colors,
      coloring)`` is returned, where ``coloring`` is a dictionary associating
      its color (integer) to each vertex of the graph.

    - ``solver`` -- (default: ``None``); specify a Linear Program (LP) solver to
      be used. If set to ``None``, the default one is used. For more information
      on LP solvers and which default solver is used, see the method
      :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the
      class :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- integer (default: ``0``); sets the level of verbosity. Set
      to 0 by default, which means quiet.

    ALGORITHM:

    Integer Linear Program.

    EXAMPLES:

    The Grundy number of a `P_4` is equal to 3::

        sage: from sage.graphs.graph_coloring import grundy_coloring
        sage: g = graphs.PathGraph(4)
        sage: grundy_coloring(g, 4)
        3

    The Grundy number of the PetersenGraph is equal to 4::

        sage: g = graphs.PetersenGraph()
        sage: grundy_coloring(g, 5)
        4

    It would have been sufficient to set the value of ``k`` to 4 in
    this case, as `4 = \Delta(G)+1`.
    """
    g._scream_if_not_simple(allow_multiple_edges=True)

    p = MixedIntegerLinearProgram(solver=solver)

    # b[v,i] is set to 1 if and only if v is colored with i
    b = p.new_variable(binary=True)

    # is_used[i] is set to 1 if and only if color[i] is used by some vertex
    is_used = p.new_variable(binary=True)

    # Each vertex is in exactly one color class
    for v in g:
        p.add_constraint(p.sum(b[v,i] for i in range(k)), max=1, min=1)

    # Two adjacent vertices have different colors
    for u,v in g.edge_iterator(labels=None):
        for i in range(k):
            p.add_constraint(b[v,i] + b[u,i], max=1)

    # The following constraints ensure that if v is colored with i, then it has
    # a neighbor colored with j for every j<i

    for i in range(k):
        for j in range(i):
            for v in g:

                # If b[v,i] == 0, then the following constraint is always
                # satisfied, as a sum of binary variables is always positive.
                # If it is equal to 1, then at least one of the other variables
                # must be set to 1 too.

                p.add_constraint(p.sum(b[u,j] for u in g.neighbor_iterator(v)) - b[v,i], min=0)

    # is_used[i] can be set to 1 only if the color is used
    for i in range(k):
        p.add_constraint(p.sum(b[v,i] for v in g) - is_used[i], min=0)

    # Trying to use as many colors as possible
    p.set_objective(p.sum(is_used[i] for i in range(k)))

    try:
        obj = p.solve(log=verbose, objective_only=value_only)
        from sage.rings.integer import Integer
        obj = Integer(obj)

    except MIPSolverException:
        raise ValueError("this graph cannot be colored with k colors")

    if value_only:
        return obj

    # Building the dictionary associating its color to every vertex

    b = p.get_values(b)
    cdef dict coloring = {}

    for v in g:
        for i in range(k):
            if b[v,i] == 1:
                coloring[v] = i
                break

    return obj, coloring


def b_coloring(g, k, value_only=True, solver=None, verbose=0):
    r"""
    Compute b-chromatic numbers and b-colorings.

    This function computes a b-coloring with at most `k` colors that maximizes
    the number of colors, if such a coloring exists.

    Definition :

    Given a proper coloring of a graph `G` and a color class `C` such that none
    of its vertices have neighbors in all the other color classes, one can
    eliminate color class `C` assigning to each of its elements a missing color
    in its neighborhood.

    Let a b-vertex be a vertex with neighbors in all other colorings. Then, one
    can repeat the above procedure until a coloring is obtained where every
    color class contains a b-vertex, in which case none of the color classes can
    be eliminated with the same idea. So, one can define a b-coloring as a
    proper coloring where each color class has a b-vertex.

    In the worst case, after successive applications of the above procedure, one
    get a proper coloring that uses a number of colors equal to the the
    b-chromatic number of `G` (denoted `\chi_b(G)`): the maximum `k` such that
    `G` admits a b-coloring with `k` colors.

    A useful upper bound for calculating the b-chromatic number is the
    following. If `G` admits a b-coloring with `k` colors, then there are `k`
    vertices of degree at least `k - 1` (the b-vertices of each color
    class). So, if we set `m(G) = \max \{k | \text{there are } k \text{ vertices
    of degree at least } k - 1 \}`, we have that `\chi_b(G) \leq m(G)`.


    .. NOTE::

       This method computes a b-coloring that uses at *MOST* `k` colors. If this
       method returns a value equal to `k`, it cannot be assumed that `k` is
       equal to `\chi_b(G)`. Meanwhile, if it returns any value `k' < k`, this
       is a certificate that the Grundy number of the given graph is `k'`.

       As `\chi_b(G)\leq m(G)`, it can be assumed that `\chi_b(G) = k` if
       ``b_coloring(g, k)`` returns `k` when `k = m(G)`.

    INPUT:

    - ``k`` -- integer; maximum number of colors

    - ``value_only`` -- boolean (default: ``True``); when set to ``True``, only
      the number of colors is returned. Otherwise, the pair ``(nb_colors,
      coloring)`` is returned, where ``coloring`` is a dictionary associating
      its color (integer) to each vertex of the graph.

    - ``solver`` -- (default: ``None``); specify a Linear Program (LP) solver to
      be used. If set to ``None``, the default one is used. For more information
      on LP solvers and which default solver is used, see the method
      :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the
      class :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- integer (default: ``0``); sets the level of verbosity. Set
      to 0 by default, which means quiet.

    ALGORITHM:

    Integer Linear Program.

    EXAMPLES:

    The b-chromatic number of a `P_5` is equal to 3::

        sage: from sage.graphs.graph_coloring import b_coloring
        sage: g = graphs.PathGraph(5)
        sage: b_coloring(g, 5)
        3

    The b-chromatic number of the Petersen Graph is equal to 3::

        sage: g = graphs.PetersenGraph()
        sage: b_coloring(g, 5)
        3

    It would have been sufficient to set the value of ``k`` to 4 in this case,
    as `4 = m(G)`.
    """
    g._scream_if_not_simple(allow_multiple_edges=True)

    # Calculate the upper bound m(G)
    # To do so, it takes the list of degrees in non-increasing order and
    # computes the largest i, such that the ith degree on the list is at least
    # i - 1 (note that in the code we need to take in consideration that the
    # indices of the list starts with 0)

    cdef list deg = g.degree_sequence()
    cdef int n = g.order()
    for i in range(n):
        if deg[i] < i:
            break
    if i != n - 1:
        m = i
    else:
        m = n

    # In case the k specified by the user is greater than m(G), make k = m(G)
    if k > m:
        k = m


    p = MixedIntegerLinearProgram(solver=solver)

    # color[v,i] is set to 1 if and only if v is colored i
    color = p.new_variable(binary=True)

    # b[v,i] is set to 1 if and only if v is a b-vertex from color class i
    b = p.new_variable(binary=True)

    # is_used[i] is set to 1 if and only if color[i] is used by some vertex
    is_used = p.new_variable(binary=True)

    # Each vertex is in exactly one class
    for v in g:
        p.add_constraint(p.sum(color[v,i] for i in range(k)), min=1, max=1)

    # Adjacent vertices have distinct colors
    for u, v in g.edge_iterator(labels=None):
        for i in range(k):
            p.add_constraint(color[u,i] + color[v,i], max=1)

    # The following constraints ensure that if v is a b-vertex of color i
    # then it has a neighbor colored j for every j != i

    for v in g:
        for i in range(k):
            for j in range(k):
                if j != i:
                    # If v is not a b-vertex of color i, the constraint
                    # is always satisfied, since the only possible
                    # negative term in this case is -is_used[j] which is
                    # cancelled by + 1. If v is a b-vertex of color i
                    # then we MUST have sum(color[w,j] for w in g.neighbors(v))
                    # valued at least 1, which means that v has a neighbour in
                    # color j, as desired.
                    p.add_constraint(p.sum(color[w,j] for w in g.neighbor_iterator(v))
                                         - b[v,i] + 1 - is_used[j], min=0)

    # if color i is used, there is a vertex colored i
    for i in range(k):
        p.add_constraint(p.sum(color[v,i] for v in g) - is_used[i], min=0)

    # if there is a vertex colored with color i, then i is used
    for v in g:
        for i in range(k):
            p.add_constraint(color[v,i] - is_used[i], max=0)


    # a color class is used if and only if it has one b-vertex
    for i in range(k):
       p.add_constraint(p.sum(b[w,i] for w in g) - is_used[i], min=0, max=0)


    # We want to maximize the number of used colors
    p.set_objective(p.sum(is_used[i] for i in range(k)))


    try:
        obj = p.solve(log=verbose, objective_only=value_only)
        from sage.rings.integer import Integer
        obj = Integer(obj)

    except MIPSolverException:
        raise ValueError("this graph cannot be colored with k colors")

    if value_only:
        return obj


    # Building the dictionary associating its color to every vertex

    c = p.get_values(color)
    cdef dict coloring = {}

    for v in g:
        for i in range(k):
            if c[v,i] == 1:
                coloring[v] = i
                break

    return obj, coloring

def edge_coloring(g, value_only=False, vizing=False, hex_colors=False, solver=None, verbose=0):
    r"""
    Compute chromatic index and edge colorings.

    INPUT:

    - ``g`` -- a graph.

    - ``value_only`` -- boolean (default: ``False``):

      - When set to ``True``, only the chromatic index is returned

      - When set to ``False``, a partition of the edge set into matchings is
        returned if possible

    - ``vizing`` -- boolean (default: ``False``):

      - When set to ``True``, tries to find a `\Delta + 1`-edge-coloring, where
        `\Delta` is equal to the maximum degree in the graph

      - When set to ``False``, tries to find a `\Delta`-edge-coloring, where
        `\Delta` is equal to the maximum degree in the graph. If impossible,
        tries to find and returns a `\Delta + 1`-edge-coloring. This implies
        that ``value_only=False``

    - ``hex_colors`` -- boolean (default: ``False``); when set to ``True``, the
      partition returned is a dictionary whose keys are colors and whose values
      are the color classes (ideal for plotting)

    - ``solver`` -- (default: ``None``); specify a Linear Program (LP) solver to
      be used. If set to ``None``, the default one is used. For more information
      on LP solvers and which default solver is used, see the method
      :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the
      class :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- integer (default: ``0``); sets the level of verbosity. Set
      to 0 by default, which means quiet.

    OUTPUT:

    In the following, `\Delta` is equal to the maximum degree in the graph
    ``g``.

    - If ``vizing=True`` and ``value_only=False``, return a partition of the
      edge set into `\Delta + 1` matchings.

    - If ``vizing=False`` and ``value_only=True``, return the chromatic index.

    - If ``vizing=False`` and ``value_only=False``, return a partition of the
      edge set into the minimum number of matchings.

    - If ``vizing=True`` and ``value_only=True``, should return something, but
      mainly you are just trying to compute the maximum degree of the graph, and
      this is not the easiest way. By Vizing's theorem, a graph has a chromatic
      index equal to `\Delta` or to `\Delta + 1`.

    .. NOTE::

       In a few cases, it is possible to find very quickly the chromatic index
       of a graph, while it remains a tedious job to compute a corresponding
       coloring. For this reason, ``value_only = True`` can sometimes be much
       faster, and it is a bad idea to compute the whole coloring if you do not
       need it !

    .. SEEALSO::

        - :wikipedia:`Edge_coloring` for further details on edge coloring
        - :meth:`~Graph.chromatic_index`
        - :meth:`~Graph.fractional_chromatic_index`
        - :meth:`~Graph.chromatic_number`
        - :meth:`sage.graphs.graph_coloring.vertex_coloring`

    EXAMPLES:

    The Petersen graph has chromatic index 4::

       sage: from sage.graphs.graph_coloring import edge_coloring
       sage: g = graphs.PetersenGraph()
       sage: edge_coloring(g, value_only=True, solver='GLPK')
       4
       sage: edge_coloring(g, value_only=False, solver='GLPK')
       [[(0, 1), (2, 3), (4, 9), (5, 7), (6, 8)],
        [(0, 4), (1, 2), (3, 8), (6, 9)],
        [(0, 5), (2, 7)],
        [(1, 6), (3, 4), (5, 8), (7, 9)]]
       sage: edge_coloring(g, value_only=False, hex_colors=True, solver='GLPK')
       {'#00ffff': [(0, 5), (2, 7)],
        '#7f00ff': [(1, 6), (3, 4), (5, 8), (7, 9)],
        '#7fff00': [(0, 4), (1, 2), (3, 8), (6, 9)],
        '#ff0000': [(0, 1), (2, 3), (4, 9), (5, 7), (6, 8)]}

    Complete graphs are colored using the linear-time round-robin coloring::

       sage: from sage.graphs.graph_coloring import edge_coloring
       sage: len(edge_coloring(graphs.CompleteGraph(20)))
       19

    The chromatic index of a non connected graph is the maximum over its
    connected components::

       sage: g = graphs.CompleteGraph(4) + graphs.CompleteGraph(10)
       sage: edge_coloring(g, value_only=True)
       9

    TESTS:

    Graph without edge::

       sage: g = Graph(2)
       sage: edge_coloring(g)
       []
       sage: edge_coloring(g, value_only=True)
       0
       sage: edge_coloring(g, hex_colors=True)
       {}
    """
    g._scream_if_not_simple()
    from sage.plot.colors import rainbow

    if not g.order() or not g.size():
        if value_only:
            return 0
        return dict() if hex_colors else list()

    if vizing:
        value_only = False

    # The chromatic index of g is the maximum value over its connected
    # components, and the edge coloring is the union of the edge
    # coloring of its connected components
    cdef list L = [g] if g.is_connected() else g.connected_components_subgraphs()
    cdef int chi = 0
    cdef list classes = [], vertices
    cdef list values
    for h in L:

        if not h.size():
            continue

        # We get the vertex of maximum degree and its degree
        Delta,X = max([(d, v) for v,d in h.degree_iterator(labels=True)], key=lambda x: x[0])

        if value_only:
            if Delta + 1 <= chi:
                continue
            if h.is_overfull():
                chi = max(chi, Delta + 1)
                continue

        if h.is_clique():
            if value_only:
                chi = max(chi, h.order() if h.order() % 2 else (h.order() - 1))
                continue
            vertices = list(h)
            r = round_robin(h.order())
            # create missing color classes, if any
            for i in range(len(classes), max(r.edge_labels()) + 1):
                classes.append([])
            # add edges to classes
            r.relabel(perm=vertices, inplace=True)
            for u, v, c in r.edge_iterator():
                classes[c].append((u, v))
            continue

        # Vizing's coloring uses Delta + 1 colors. Otherwise, we try both.
        values = [Delta + 1] if vizing else [Delta, Delta + 1]

        for k in values:
            p = MixedIntegerLinearProgram(maximization=True, solver=solver)
            color = p.new_variable(binary=True)
            # A vertex cannot have two incident edges with the same color.
            for v in h:
                for i in range(k):
                    p.add_constraint(p.sum(color[frozenset((u,v)),i] for u in h.neighbor_iterator(v)) <= 1)
            # An edge must have a color
            for u,v in h.edge_iterator(labels=False):
                p.add_constraint(p.sum(color[frozenset((u,v)),i] for i in range(k)) == 1)
            # We color the edges of the vertex of maximum degree
            for i,v in enumerate(h.neighbor_iterator(X)):
                p.add_constraint( color[frozenset((v,X)),i] == 1 )
            try:
                p.solve(objective_only=value_only, log=verbose)
                break
            except MIPSolverException:
                if k == Delta + 1:
                    raise RuntimeError("Something is wrong! Certainly a problem in the"
                                           " algorithm... please contact sage-devel@googlegroups.com")
                # The coloring fails with Delta colors
                if value_only:
                    k = k + 1
                    break

        if value_only:
            chi = max(chi, k)
        else:
            # create missing color classes, if any
            for i in range(len(classes), k):
                classes.append([])
            # add edges to color classes
            color = p.get_values(color)
            for e in h.edge_iterator(labels=False):
                fe = frozenset(e)
                for i in range(k):
                    if color[fe,i] == 1:
                        classes[i].append(e)
                        break

    if value_only:
        return chi
    # if needed, builds a dictionary from the color classes adding colors
    if hex_colors:
        return dict(zip(rainbow(len(classes)), classes))
    else:
        return classes

def round_robin(n):
    r"""
    Compute a round-robin coloring of the complete graph on `n` vertices.

    A round-robin coloring of the complete graph `G` on `2n` vertices
    (`V = [0, \dots, 2n - 1]`) is a proper coloring of its edges such that
    the edges with color `i` are all the `(i + j, i - j)` plus the
    edge `(2n - 1, i)`.

    If `n` is odd, one obtain a round-robin coloring of the complete graph
    through the round-robin coloring of the graph with `n + 1` vertices.

    INPUT:

    - ``n`` -- the number of vertices in the complete graph

    OUTPUT:

    - A :meth:`~sage.graphs.graph_generators.GraphGenerators.CompleteGraph` with
      labelled edges such that the label of each edge is its color.

    EXAMPLES::

        sage: from sage.graphs.graph_coloring import round_robin
        sage: round_robin(3).edges()
        [(0, 1, 2), (0, 2, 1), (1, 2, 0)]

    ::

        sage: round_robin(4).edges()
        [(0, 1, 2), (0, 2, 1), (0, 3, 0), (1, 2, 0), (1, 3, 1), (2, 3, 2)]


    For higher orders, the coloring is still proper and uses the expected
    number of colors::

        sage: g = round_robin(9)
        sage: sum(Set(e[2] for e in g.edges_incident(v)).cardinality() for v in g) == 2 * g.size()
        True
        sage: Set(e[2] for e in g.edge_iterator()).cardinality()
        9

    ::

        sage: g = round_robin(10)
        sage: sum(Set(e[2] for e in g.edges_incident(v)).cardinality() for v in g) == 2 * g.size()
        True
        sage: Set(e[2] for e in g.edge_iterator()).cardinality()
        9
    """
    if n <= 1:
        raise ValueError("there must be at least two vertices in the graph")
    def my_mod(x, y):
        return x - y * (x // y)
    if not n % 2:
        from sage.graphs.generators.basic import CompleteGraph
        g = CompleteGraph(n)
        for i in range(n - 1):
            g.set_edge_label(n - 1, i, i)
            for j in range(1, (n - 1) // 2 + 1):
                g.set_edge_label(my_mod(i - j, n - 1), my_mod(i + j, n - 1), i)
        return g
    else:
        g = round_robin(n + 1)
        g.delete_vertex(n)
        return g

def linear_arboricity(g, plus_one=None, hex_colors=False, value_only=False, solver=None, verbose=0):
    r"""
    Compute the linear arboricity of the given graph.

    The linear arboricity of a graph `G` is the least number `la(G)` such that
    the edges of `G` can be partitioned into linear forests (i.e. into forests
    of paths).

    Obviously, `la(G)\geq \left\lceil \frac{\Delta(G)}{2} \right\rceil`.

    It is conjectured in [Aki1980]_ that `la(G)\leq \left\lceil
    \frac{\Delta(G)+1}{2} \right\rceil`.

    INPUT:

    - ``plus_one`` -- integer (default: ``None``); whether to use `\left\lceil
      \frac{\Delta(G)}{2} \right\rceil` or `\left\lceil \frac{\Delta(G)+1}{2}
      \right\rceil` colors.

      - If ``0``, computes a decomposition of `G` into `\left\lceil
        \frac{\Delta(G)}{2} \right\rceil` forests of paths

      - If ``1``, computes a decomposition of `G` into `\left\lceil
        \frac{\Delta(G)+1}{2} \right\rceil` colors, which is the conjectured
        general bound.

      - If ``plus_one = None`` (default), computes a decomposition using the
        least possible number of colors.

    - ``hex_colors`` -- boolean (default: ``False``):

      - If ``hex_colors = True``, the function returns a dictionary associating
        to each color a list of edges (meant as an argument to the
        ``edge_colors`` keyword of the ``plot`` method).

      - If ``hex_colors = False`` (default value), returns a list of graphs
        corresponding to each color class.

    - ``value_only`` -- boolean (default: ``False``):

      - If ``value_only = True``, only returns the linear arboricity as an
        integer value.

      - If ``value_only = False``, returns the color classes according to the
        value of ``hex_colors``

    - ``solver`` -- (default: ``None``); specify a Linear Program (LP) solver to
      be used. If set to ``None``, the default one is used. For more information
      on LP solvers and which default solver is used, see the method
      :meth:`~sage.numerical.mip.MixedIntegerLinearProgram.solve` of the class
      :class:`~sage.numerical.mip.MixedIntegerLinearProgram`.

    - ``verbose`` -- integer (default: ``0``); sets the level of verbosity of
      the LP solver. Set to 0 by default, which means quiet.

    ALGORITHM:

    Linear Programming

    COMPLEXITY:

    NP-Hard

    EXAMPLES:

    Obviously, a square grid has a linear arboricity of 2, as the set of
    horizontal lines and the set of vertical lines are an admissible partition::

        sage: from sage.graphs.graph_coloring import linear_arboricity
        sage: g = graphs.Grid2dGraph(4, 4)
        sage: g1,g2 = linear_arboricity(g)

    Each graph is of course a forest::

        sage: g1.is_forest() and g2.is_forest()
        True

    Of maximum degree 2::

        sage: max(g1.degree()) <= 2 and max(g2.degree()) <= 2
        True

    Which constitutes a partition of the whole edge set::

        sage: all((g1.has_edge(e) or g2.has_edge(e)) for e in g.edge_iterator(labels=None))
        True

    TESTS:

    Asking for the value of the linear arboricity only (:trac:`24991`)::

        sage: from sage.graphs.graph_coloring import linear_arboricity
        sage: sorted(linear_arboricity(G, value_only=True) for G in graphs(4))
        [0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]

    Test parameter ``hex_color`` (:trac:`26228`)::

        sage: from sage.graphs.graph_coloring import linear_arboricity
        sage: g = graphs.Grid2dGraph(4, 4)
        sage: d = linear_arboricity(g, hex_colors=True)
        sage: sorted(d)
        ['#00ffff', '#ff0000']
    """
    g._scream_if_not_simple()
    from sage.rings.integer import Integer

    if plus_one is None:
        try:
            return linear_arboricity(g,
                                     plus_one=0,
                                     value_only=value_only,
                                     hex_colors=hex_colors,
                                     solver=solver,
                                     verbose=verbose)
        except ValueError:
            return linear_arboricity(g,
                                     plus_one=1,
                                     value_only=value_only,
                                     hex_colors=hex_colors,
                                     solver=solver,
                                     verbose=verbose)
    elif plus_one == 1:
        k = (Integer(1 + max(g.degree())) / 2).ceil()
    elif not plus_one:
        k = (Integer(max(g.degree())) / 2).ceil()
    else:
        raise ValueError("plus_one must be equal to 0,1, or to None!")

    from sage.plot.colors import rainbow

    p = MixedIntegerLinearProgram(solver=solver)

    # c is a boolean value such that c[i,(u,v)] = 1 if and only if (u,v) is
    # colored with i
    c = p.new_variable(binary=True)

    # relaxed value
    r = p.new_variable(nonnegative=True)

    MAD = 1 - 1 / (Integer(g.order()) * 2)

    # Partition of the edges
    for u,v in g.edge_iterator(labels=None):
        p.add_constraint(p.sum(c[i,frozenset((u,v))] for i in range(k)), max=1, min=1)

    for i in range(k):

        # r greater than c
        for u,v in g.edge_iterator(labels=None):
            p.add_constraint(r[i,(u,v)] + r[i,(v,u)] - c[i,frozenset((u,v))], max=0, min=0)


        # Maximum degree 2
        for u in g:
            p.add_constraint(p.sum(c[i,frozenset((u,v))] for v in g.neighbor_iterator(u)), max=2)

            # no cycles
            p.add_constraint(p.sum(r[i,(u,v)] for v in g.neighbor_iterator(u)),max=MAD)

    try:
        p.solve(objective_only=value_only, log=verbose)
        if value_only:
            return k

    except MIPSolverException:
        if plus_one:
            raise RuntimeError("It looks like you have found a counterexample to a very old conjecture. Please do not loose it ! Please publish it, and send a post to sage-devel to warn us. We implore you!")
        else:
            raise ValueError("this graph cannot be colored with the given number of colors")

    c = p.get_values(c)

    cdef list answer

    if hex_colors:
        answer = [[] for i in range(k)]
        def add(uv, i):
            return answer[i].append(uv)
    else:
        gg = copy(g)
        gg.delete_edges(g.edge_iterator())
        answer = [copy(gg) for i in range(k)]
        def add(uv, i):
            return answer[i].add_edge(uv)

    for i in range(k):
        for u,v in g.edge_iterator(labels=None):
            if c[i,frozenset((u,v))]  == 1:
                add((u,v),i)

    if hex_colors:
        return dict(zip(rainbow(len(answer)), answer))
    else:
        return answer


def acyclic_edge_coloring(g, hex_colors=False, value_only=False, k=0, solver=None, verbose=0):
    r"""
    Compute an acyclic edge coloring of the current graph.

    An edge coloring of a graph is a assignment of colors to the edges of a
    graph such that :

    - the coloring is proper (no adjacent edges share a color)
    - For any two colors `i,j`, the union of the edges colored with `i` or `j`
      is a forest.

    The least number of colors such that such a coloring exists for a graph `G`
    is written `\chi'_a(G)`, also called the acyclic chromatic index of `G`.

    It is conjectured that this parameter cannot be too different from the
    obvious lower bound `\Delta(G) \leq \chi'_a(G)`, `\Delta(G)` being the
    maximum degree of `G`, which is given by the first of the two constraints.
    Indeed, it is conjectured that `\Delta(G)\leq \chi'_a(G)\leq \Delta(G) + 2`.

    INPUT:

    - ``hex_colors`` -- boolean (default: ``False``):

      - If ``hex_colors = True``, the function returns a dictionary associating
        to each color a list of edges (meant as an argument to the
        ``edge_colors`` keyword of the ``plot`` method).

      - If ``hex_colors = False`` (default value), returns a list of graphs
        corresponding to each color class.

    - ``value_only`` -- boolean (default: ``False``):

      - If ``value_only = True``, only returns the acyclic chromatic index as an
        integer value

      - If ``value_only = False``, returns the color classes according to the
        value of ``hex_colors``

    - ``k`` -- integer; the number of colors to use.

      - If ``k > 0``, computes an acyclic edge coloring using `k` colors.

      - If ``k = 0`` (default), computes a coloring of `G` into `\Delta(G) + 2`
        colors, which is the conjectured general bound.

      - If ``k = None``, computes a decomposition using the least possible
        number of colors.

    - ``solver`` -- (default: ``None``); specify a Linear Program (LP) solver to
      be used. If set to ``None``, the default one is used. For more information
      on LP solvers and which default solver is used, see the method
      :meth:`~sage.numerical.mip.MixedIntegerLinearProgram.solve` of the class
      :class:`~sage.numerical.mip.MixedIntegerLinearProgram`.

    - ``verbose`` -- integer (default: ``0``); sets the level of verbosity of
      the LP solver. Set to 0 by default, which means quiet.

    ALGORITHM:

    Linear Programming

    EXAMPLES:

    The complete graph on 8 vertices cannot be acyclically edge-colored with
    less `\Delta + 1` colors, but it can be colored with `\Delta + 2 = 9`::

        sage: from sage.graphs.graph_coloring import acyclic_edge_coloring
        sage: g = graphs.CompleteGraph(8)
        sage: colors = acyclic_edge_coloring(g)

    Each color class is of course a matching ::

        sage: all(max(gg.degree()) <= 1 for gg in colors)
        True

    These matchings being a partition of the edge set::

        sage: all(any(gg.has_edge(e) for gg in colors) for e in g.edge_iterator(labels=False))
        True

    Besides, the union of any two of them is a forest ::

        sage: all(g1.union(g2).is_forest() for g1 in colors for g2 in colors)
        True

    If one wants to acyclically color a cycle on `4` vertices, at least 3 colors
    will be necessary. The function raises an exception when asked to color it
    with only 2::

        sage: g = graphs.CycleGraph(4)
        sage: acyclic_edge_coloring(g, k=2)
        Traceback (most recent call last):
        ...
        ValueError: this graph cannot be colored with the given number of colors

    The optimal coloring give us `3` classes::

        sage: colors = acyclic_edge_coloring(g, k=None)
        sage: len(colors)
        3

    TESTS:

    Ticket :trac:`24991` is fixed::

        sage: from sage.graphs.graph_coloring import acyclic_edge_coloring
        sage: sorted(acyclic_edge_coloring(G, value_only=True) for G in graphs(4))
        [2, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5]

    Test parameter ``hex_color`` (:trac:`26228`)::

        sage: from sage.graphs.graph_coloring import acyclic_edge_coloring
        sage: g = graphs.CompleteGraph(4)
        sage: d = acyclic_edge_coloring(g, hex_colors=True)
        sage: sorted(d)
        ['#0066ff', '#00ff66', '#cbff00', '#cc00ff', '#ff0000']

    The acyclic chromatic index of a graph without edge is 0 (:trac:`27079`)::

        sage: from sage.graphs.graph_coloring import acyclic_edge_coloring
        sage: g = Graph(3)
        sage: acyclic_edge_coloring(g, k=None, value_only=True)
        0
        sage: acyclic_edge_coloring(g, k=None, hex_colors=True)
        {}
        sage: acyclic_edge_coloring(g, k=None, hex_colors=False)
        []

    Empty graph  (:trac:`27079`)::

        sage: from sage.graphs.graph_coloring import acyclic_edge_coloring
        sage: acyclic_edge_coloring(Graph(), k=None, value_only=True)
        0
    """
    g._scream_if_not_simple(allow_multiple_edges=True)

    from sage.rings.integer import Integer
    from sage.combinat.subset import Subsets
    from sage.plot.colors import rainbow

    if not g.order() or not g.size():
        if k == 0:
            k = 2
        if value_only:
            return 0 if k is None else k
        else:
            if k is None:
                return {} if hex_colors else []
            return {c: [] for c in rainbow(k)} if hex_colors else [copy(g) for _ in range(k)]

    if k is None:
        k = max(g.degree())

        while True:
            try:
                return acyclic_edge_coloring(g,
                                             value_only=value_only,
                                             hex_colors=hex_colors,
                                             k=k,
                                             solver=solver,
                                             verbose=verbose)
            except ValueError:
                k += 1

        raise RuntimeError("this should not happen, please report a bug!")

    elif not k:
        k = max(g.degree()) + 2

    p = MixedIntegerLinearProgram(solver=solver)

    # c is a binary variable such that c[i,(u,v)] = 1 if and only if (u,v) is
    # colored with i
    c = p.new_variable(binary=True)

    # relaxed value
    r = p.new_variable(nonnegative=True)

    def E(x, y):
        return frozenset((x,y))

    MAD = 1 - 1/(Integer(g.order()) * 2)

    # Partition of the edges: each edge is assigned a unique color
    for u,v in g.edge_iterator(labels=None):
        p.add_constraint(p.sum(c[i,E(u,v)] for i in range(k)), max=1, min=1)


    for i in range(k):

        # Maximum degree 1
        for u in g:
            p.add_constraint(p.sum(c[i,E(u,v)] for v in g.neighbor_iterator(u)), max=1)

    for i,j in Subsets(range(k), 2):
        # r is greater than c
        for u in g:
            p.add_constraint(p.sum(r[(i,j),(u,v)] for v in g.neighbor_iterator(u)), max=MAD)

        # r greater than c
        for u,v in g.edge_iterator(labels=None):
            p.add_constraint(r[(i,j),(u,v)] + r[(i,j),(v,u)] - c[i,E(u,v)] - c[j,E(u,v)], max=0, min=0)

    p.set_objective(None)

    try:
        p.solve(objective_only=value_only, log=verbose)
        if value_only:
            return k

    except MIPSolverException:
        if k == max(g.degree()) + 2:
            raise RuntimeError("It looks like you have found a counterexample to "
                               "a very old conjecture. Please do not loose it ! "
                               "Please publish it, and send a post to sage-devel "
                               "to warn us. We implore you!")
        else:
            raise ValueError("this graph cannot be colored with the given number of colors")

    c = p.get_values(c)

    if hex_colors:
        answer = [[] for i in range(k)]
        def add(uv, i):
            return answer[i].append(uv)
    else:
        gg = copy(g)
        gg.delete_edges(g.edge_iterator())
        answer = [copy(gg) for i in range(k)]
        def add(uv, i):
            return answer[i].add_edge(uv)

    for i in range(k):
        for u,v in g.edge_iterator(labels=None):
            if c[i,E(u,v)] == 1:
                add((u,v), i)

    if hex_colors:
        return dict(zip(rainbow(len(answer)), answer))
    else:
        return answer


cdef class Test:
    r"""
    This class performs randomized testing for all_graph_colorings.

    Since everything else in this file is derived from all_graph_colorings, this
    is a pretty good randomized tester for the entire file.  Note that for a
    graph `G`, ``G.chromatic_polynomial()`` uses an entirely different
    algorithm, so we provide a good, independent test.
    """

    def random(self, tests=1000):
        r"""
        Call ``self.random_all_graph_colorings()``.

        In the future, if other methods are added, it should call them, too.

        TESTS::

            sage: from sage.graphs.graph_coloring import Test
            sage: Test().random(1)
        """
        self.random_all_graph_colorings(tests)

    def random_all_graph_colorings(self, tests=2):
        r"""
        Verify the results of ``all_graph_colorings()`` in three ways:

        #. all colorings are unique

        #. number of m-colorings is `P(m)` (where `P` is the chromatic
           polynomial of the graph being tested)

        #. colorings are valid -- that is, that no two vertices of the same
           color share an edge.

        TESTS::

            sage: from sage.graphs.graph_coloring import Test
            sage: Test().random_all_graph_colorings(1)
        """
        from sage.graphs.generators.random import RandomGNP
        cdef set S
        cdef list parts
        for _ in range(tests):
            G = RandomGNP(10, .5)
            Q = G.chromatic_polynomial()
            chi = G.chromatic_number()

            S = set()

            for C in all_graph_colorings(G, chi):
                parts = [C[k] for k in C]
                for P in parts:
                    l = len(P)
                    for i in range(l):
                        for j in range(i+1, l):
                            if G.has_edge(P[i], P[j]):
                                raise RuntimeError("coloring failed")

                # make the dict into a frozenset for quick uniqueness checking
                S.add(frozenset((k, frozenset(C[k])) for k in C))

            if len(S) != Q(chi):
                raise RuntimeError("incorrect number of unique colorings")
