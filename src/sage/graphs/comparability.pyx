# cython: binding=True
r"""
Comparability and permutation graphs

This module implements method related to :wikipedia:`Comparability_graph` and
:wikipedia:`Permutation_graph`, that is, for the moment, only recognition
algorithms.

Most of the information found here can also be found in [ST1994]_ or [Sha1997]_.

The following methods are implemented in this module

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~is_comparability_MILP` | Check whether the graph is a comparability graph (MILP)
    :meth:`~greedy_is_comparability` | Check whether the graph is a comparability graph (greedy algorithm)
    :meth:`~greedy_is_comparability_with_certificate` | Check whether the graph is a comparability graph and returns certificates (greedy algorithm)
    :meth:`~is_comparability` | Check whether the graph is a comparability graph
    :meth:`~is_permutation` | Check whether the graph is a permutation graph.
    :meth:`~is_transitive` | Check whether the digraph is transitive.

Author:

- Nathann Cohen 2012-04

Graph classes
-------------

**Comparability graphs**

A graph is a comparability graph if it can be obtained from a poset by adding an
edge between any two elements that are comparable. Co-comparability graph are
complements of such graphs, i.e. graphs built from a poset by adding an edge
between any two incomparable elements.

For more information on comparability graphs, see the
:wikipedia:`Comparability_graph`.

**Permutation graphs**

Definitions:

- A permutation `\pi = \pi_1\pi_2\dots\pi_n` defines a graph on `n` vertices
  such that `i\sim j` when `\pi` reverses `i` and `j` (i.e. when `i<j` and
  `\pi_j < \pi_i`. A graph is a permutation graph whenever it can be built
  through this construction.

- A graph is a permutation graph if it can be build from two parallel lines are
  the intersection graph of segments intersecting both lines.

- A graph is a permutation graph if it is both a comparability graph and a
  co-comparability graph.

For more information on permutation graphs, see the
:wikipedia:`Permutation_graph`.


Recognition algorithm for comparability graphs
----------------------------------------------

**Greedy algorithm**

This algorithm attempts to build a transitive orientation of a given graph `G`,
that is an orientation `D` such that for any directed `uv`-path of `D` there
exists in `D` an edge `uv`. This already determines a notion of equivalence
between some edges of `G` :

  In `G`, two edges `uv` and `uv'` (incident to a common vertex `u`) such that
  `vv'\not\in G` need necessarily be oriented *the same way* (that is that they
  should either both *leave* or both *enter* `u`). Indeed, if one enters `G`
  while the other leaves it, these two edges form a path of length two, which is
  not possible in any transitive orientation of `G` as `vv'\not\in G`.

Hence, we can say that in this case a *directed edge* `uv` is equivalent to a
*directed edge* `uv'` (to mean that if one belongs to the transitive
orientation, the other one must be present too) in the same way that `vu` is
equivalent to `v'u`. We can thus define equivalence classes on oriented edges,
to represent set of edges that imply each other. We can thus define `C^G_{uv}`
to be the equivalence class in `G` of the oriented edge `uv`.

Of course, if there exists a transitive orientation of a graph `G`, then no edge
`uv` implies its contrary `vu`, i.e. it is necessary to ensure that `\forall
uv\in G, vu\not\in C^G_{uv}`. The key result on which the greedy algorithm is
built is the following (see [ST1994]_):

  **Theorem** -- The following statements are equivalent :

  - `G` is a comparability graph
  - `\forall uv\in G, vu\not\in C^G_{uv}`
  - The edges of `G` can be partitioned into `B_1,...,B_k` where `B_i` is the
    equivalence class of some oriented edge in `G-B_1-\dots-B_{i-1}`

Hence, ensuring that a graph is a comparability graph can be done by checking
that no equivalence class is contradictory. Building the orientation, however,
requires to build equivalence classes step by step until an orientation has been
found for all of them.

**Mixed Integer Linear Program**

A MILP formulation is available to check the other methods for correction. It is
easily built :

  To each edge are associated two binary variables (one for each possible
  direction).  We then ensure that each triangle is transitively oriented, and
  that each pair of incident edges `uv, uv'` such that `vv'\not\in G` do not
  create a 2-path.

Here is the formulation:

.. MATH::

    \mbox{Maximize : }&\mbox{Nothing}\\
    \mbox{Such that : }&\\
    &\forall uv\in G\\
    &\cdot o_{uv}+o_{vu} = 1\\
    &\forall u\in G, \forall v,v'\in N(v)\text{ such that }vv'\not\in G\\
    &\cdot o_{uv} + o_{v'u} - o_{v'v} \leq 1\\
    &\cdot o_{uv'} + o_{vu} - o_{vv'} \leq 1\\
    &\forall u\in G, \forall v,v'\in N(v)\text{ such that }vv'\in G\\
    &\cdot o_{uv} + o_{v'u} \leq  1\\
    &\cdot o_{uv'} + o_{vu} \leq  1\\
    &o_{uv}\text{ is a binary variable}\\

.. NOTE::

  The MILP formulation is usually much slower than the greedy algorithm. This
  MILP has been implemented to check the results of the greedy algorithm that
  has been implemented to check the results of a faster algorithm which has not
  been implemented yet.

Certificates
------------

**Comparability graphs**

The *yes*-certificates that a graph is a comparability graphs are transitive
orientations of it. The *no*-certificates, on the other hand, are odd cycles of
such graph. These odd cycles have the property that around each vertex `v` of
the cycle its two incident edges must have the same orientation (toward `v`, or
outward `v`) in any transitive orientation of the graph. This is impossible
whenever the cycle has odd length. Explanations are given in the "Greedy
algorithm" part of the previous section.

**Permutation graphs**

Permutation graphs are precisely the intersection of comparability graphs and
co-comparability graphs. Hence, negative certificates are precisely negative
certificates of comparability or co-comparability. Positive certificates are a
pair of permutations that can be used through
:meth:`~sage.graphs.graph_generators.GraphGenerators.PermutationGraph` (whose
documentation says more about what these permutations represent).

Implementation details
----------------------

**Test that the equivalence classes are not self-contradictory**

This is done by a call to :meth:`Graph.is_bipartite`, and here is how :

   Around a vertex `u`, any two edges `uv, uv'` such that `vv'\not\in G` are
   equivalent. Hence, the equivalence class of edges around a vertex are
   precisely the connected components of the complement of the graph induced by
   the neighbors of `u`.

   In each equivalence class (around a given vertex `u`), the edges should all
   have the same orientation, i.e. all should go toward `u` at the same time, or
   leave it at the same time. To represent this, we create a graph with vertices
   for all equivalent classes around all vertices of `G`, and link `(v, C)` to
   `(u, C')` if `u\in C` and `v\in C'`.

   A bipartite coloring of this graph with colors 0 and 1 tells us that the
   edges of an equivalence class `C` around `u` should be directed toward `u` if
   `(u, C)` is colored with `0`, and outward if `(u, C)` is colored with `1`.

   If the graph is not bipartite, this is the proof that some equivalence class
   is self-contradictory !


.. NOTE::

    The greedy algorithm implemented here is just there to check the correction
    of more complicated ones, and it is reaaaaaaaaaaaalllly bad whenever you
    look at it with performance in mind.

Methods
-------
"""

# ****************************************************************************
#       Copyright (C) 2012 Nathann Cohen <nathann.cohen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.memory cimport sig_free
from sage.graphs.distances_all_pairs cimport c_distances_all_pairs
from copy import copy


#####################
# Greedy Algorithms #
#####################

def greedy_is_comparability(g, no_certificate=False, equivalence_class=False):
    r"""
    Check whether the graph is a comparability graph (greedy algorithm).

    This method only returns no-certificates.

    To understand how this method works, please consult the documentation of the
    :mod:`comparability module <sage.graphs.comparability>`.

    INPUT:

    - ``g`` -- a graph

    - ``no_certificate`` -- whether to return a *no*-certificate when the graph
      is not a comparability graph. This certificate is an odd cycle of edges,
      each of which implies the next. It is set to ``False`` by default.

    - ``equivalence_class`` -- whether to return an equivalence class
      if the graph is a comparability graph

    OUTPUT:

    - If the graph is a comparability graph and ``no_certificate = False``, this
      method returns ``True`` or ``(True, an_equivalence_class)`` according to
      the value of ``equivalence_class``.

    - If the graph is *not* a comparability graph, this method returns ``False``
      or ``(False, odd_cycle)`` according to the value of ``no_certificate``.

    EXAMPLES:

    The Petersen Graph is not transitively orientable::

        sage: from sage.graphs.comparability import greedy_is_comparability as is_comparability
        sage: g = graphs.PetersenGraph()
        sage: is_comparability(g)
        False
        sage: is_comparability(g, no_certificate=True)
        (False, [2, 1, 0, 4, 3, 2])

    But the Bull graph is::

        sage: g = graphs.BullGraph()
        sage: is_comparability(g)
        True

    TESTS:

    Check that the method is working even when vertices are of incomparable
    types::

        sage: from sage.graphs.comparability import greedy_is_comparability
        sage: G = Graph([('a', 1), (1, 2), (2, 3)])
        sage: greedy_is_comparability(G, equivalence_class=True)
        (True, [('a', 1), (2, 1), (2, 3)])
    """
    cdef int i, j

    # Each vertex can partition its neighbors into equivalence classes
    cdef dict equivalence_classes = {}
    for v in g:
        equivalence_classes[v] = g.subgraph(vertices=g.neighbors(v)).complement().connected_components(sort=False)

    # We build a graph h with one vertex per (vertex of g + equivalence class)
    from sage.graphs.graph import Graph
    h = Graph()
    h.add_vertices([(v, i) for v in g for i in range(len(equivalence_classes[v]))])

    # We add an edge between two vertices of h if they represent
    # opposed equivalence classes

    for u, v in g.edge_iterator(labels=False):

        for i, s in enumerate(equivalence_classes[v]):
            if u in s:
                break

        for j, s in enumerate(equivalence_classes[u]):
            if v in s:
                break

        h.add_edge((v, i), (u, j))

    # Is it a comparability graph ?

    cdef int isit
    isit, certif = h.is_bipartite(certificate=True)

    if isit:
        if equivalence_class:
            # We use a mapping between vertices and integers to deal with
            # vertices of different types
            int_to_vertex = list(g)
            vertex_to_int = {u: i for i, u in enumerate(int_to_vertex)}

            # Returning the largest equivalence class
            cc = max(h.connected_components(sort=False), key=len)

            edges = []
            for v, sid in cc:
                s = equivalence_classes[v][sid]

                # For each edge we pick the good orientations
                vi = vertex_to_int[v]
                if certif[v, sid] == 1:
                    edges.extend((vi, vertex_to_int[vv]) for vv in s)
                else:
                    edges.extend((vertex_to_int[vv], vi) for vv in s)

            # We return the value but take care of removing edges that were
            # added twice.
            edges = [(int_to_vertex[u], int_to_vertex[v]) for u, v in sorted(set(edges))]
            return True, edges

        return True

    if no_certificate:
        cycle = [v for v, _ in certif]
        cycle.append(cycle[0])
        return False, cycle
    return False


def greedy_is_comparability_with_certificate(g, certificate=False):
    r"""
    Check whether the graph is a comparability graph and returns
    certificates (greedy algorithm).

    This method can return certificates of both *yes* and *no* answers.

    To understand how this method works, please consult the documentation of the
    :mod:`comparability module <sage.graphs.comparability>`.

    INPUT:

    - ``g`` -- a graph

    - ``certificate`` -- boolean; whether to return a
      certificate. *Yes*-answers the certificate is a transitive orientation of
      `G`, and a *no* certificates is an odd cycle of sequentially forcing
      edges.

    EXAMPLES:

    The 5-cycle or the Petersen Graph are not transitively orientable::

        sage: from sage.graphs.comparability import greedy_is_comparability_with_certificate as is_comparability
        sage: is_comparability(graphs.CycleGraph(5), certificate=True)
        (False, [2, 1, 0, 4, 3, 2])
        sage: g = graphs.PetersenGraph()
        sage: is_comparability(g)
        False
        sage: is_comparability(g, certificate=True)
        (False, [2, 1, 0, 4, 3, 2])

    But the Bull graph is::

        sage: g = graphs.BullGraph()
        sage: is_comparability(g)
        True
        sage: is_comparability(g, certificate = True)
        (True, Digraph on 5 vertices)
        sage: is_comparability(g, certificate = True)[1].is_transitive()
        True

    TESTS:

    Check that the method is working even when vertices are of incomparable
    types::

        sage: from sage.graphs.comparability import greedy_is_comparability_with_certificate
        sage: G = Graph([('a', 1), (1, 2), (2, 3)])
        sage: greedy_is_comparability_with_certificate(G, certificate=True)
        (True, Digraph on 4 vertices)
    """
    isit, certif = greedy_is_comparability(g, no_certificate=True, equivalence_class=True)
    if not isit:
        if certificate:
            return False, certif
        return False

    elif not certificate:
        return True

    gg = copy(g)
    from sage.graphs.digraph import DiGraph
    h = DiGraph()
    h.add_vertices(gg)

    for u, v in certif:
        gg.delete_edge(u, v)
        h.add_edge(u, v)

    # While there are some edges left to be oriented
    while gg.size():

        # We take an equivalence class and orient it
        isit, certif = greedy_is_comparability(gg, no_certificate=True, equivalence_class=True)

        # Then remove it from the former graph
        for u, v in certif:
            gg.delete_edge(u, v)
            h.add_edge(u, v)

    return True, h


###################
# Integer Program #
###################

def is_comparability_MILP(g, certificate=False, solver=None, verbose=0):
    r"""
    Check whether the graph is a comparability graph (MILP).

    INPUT:

    - ``g`` -- a graph

    - ``certificate`` -- boolean (default: ``False``); whether to return a
      certificate for yes instances. This method cannot return negative
      certificates.

    - ``solver`` -- string (default: ``None``); specifies a Mixed Integer Linear
      Programming (MILP) solver to be used. If set to ``None``, the default one
      is used. For more information on MILP solvers and which default solver is
      used, see the method :meth:`solve
      <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
      :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- integer (default: 0); sets the level of verbosity. Set
      to 0 by default, which means quiet.

    EXAMPLES:

    The 5-cycle or the Petersen Graph are not transitively orientable::

        sage: from sage.graphs.comparability import is_comparability_MILP as is_comparability
        sage: is_comparability(graphs.CycleGraph(5), certificate=True)                  # needs sage.numerical.mip
        (False, None)
        sage: g = graphs.PetersenGraph()
        sage: is_comparability(g, certificate=True)                                     # needs sage.numerical.mip
        (False, None)

    But the Bull graph is::

        sage: g = graphs.BullGraph()
        sage: is_comparability(g)                                                       # needs sage.numerical.mip
        True
        sage: is_comparability(g, certificate=True)                                     # needs sage.numerical.mip
        (True, Digraph on 5 vertices)
        sage: is_comparability(g, certificate=True)[1].is_transitive()                  # needs sage.numerical.mip
        True
    """
    from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException
    p = MixedIntegerLinearProgram(solver=solver)
    o = p.new_variable(binary=True)

    for u, v in g.edge_iterator(labels=False):
        p.add_constraint(o[u, v] + o[v, u] == 1)

    from itertools import combinations
    for u in g:
        for v, vv in combinations(g.neighbors(u), 2):

            # If there is an edge between v and vv, we must be sure it is in the
            # good direction when v-u-vv is a directed path
            if g.has_edge(v, vv):
                p.add_constraint(o[u, v] + o[vv, u] - o[vv, v] <= 1)
                p.add_constraint(o[u, vv] + o[v, u] - o[v, vv] <= 1)

            # If there is no edge, there are only two orientations possible (see
            # the module's documentation about edges which imply each other)
            else:
                p.add_constraint(o[u, v] + o[vv, u] <= 1)
                p.add_constraint(o[u, vv] + o[v, u] <= 1)

    try:
        p.solve(log=verbose)
        if not certificate:
            return True

        # Building the transitive orientation
        from sage.graphs.digraph import DiGraph
        d = DiGraph()
        d.add_vertices(g)

        tol = 0 if p.base_ring().is_exact() else 1e-6
        o = p.get_values(o, convert=True, tolerance=tol)
        for u, v in g.edge_iterator(labels=False):
            if o[u, v]:
                d.add_edge(u, v)
            else:
                d.add_edge(v, u)

        return True, d

    except MIPSolverException:
        if certificate:
            return False, None
        return False


###############
# Empty shell #
###############

def is_comparability(g, algorithm='greedy', certificate=False, check=True,
                     solver=None, verbose=0):
    r"""
    Check whether the graph is a comparability graph.

    INPUT:

    - ``g`` -- a graph

    - ``algorithm`` -- string (default: ``'greedy'``); choose the implementation
      used to do the test

      - ``'greedy'`` -- a greedy algorithm (see the documentation of the
        :mod:`comparability module <sage.graphs.comparability>`)

      - ``'MILP'`` -- a Mixed Integer Linear Program formulation of the
        problem. Beware, for this implementation is unable to return negative
        certificates ! When ``certificate = True``, negative certificates are
        always equal to ``None``. ``True`` certificates are valid, though.

    - ``certificate`` -- boolean (default: ``False``); whether to return a
      certificate. *Yes*-answers the certificate is a transitive orientation of
      `G`, and a *no* certificates is an odd cycle of sequentially forcing
      edges.

    - ``check`` -- boolean (default: ``True``); whether to check that the
      yes-certificates are indeed transitive. As it is very quick
      compared to the rest of the operation, it is enabled by default.

    - ``solver`` -- string (default: ``None``); specifies a Mixed Integer Linear
      Programming (MILP) solver to be used. If set to ``None``, the default one
      is used. For more information on MILP solvers and which default solver is
      used, see the method :meth:`solve
      <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
      :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- integer (default: 0); sets the level of verbosity. Set
      to 0 by default, which means quiet.

    EXAMPLES::

        sage: from sage.graphs.comparability import is_comparability
        sage: g = graphs.PetersenGraph()
        sage: is_comparability(g)
        False
        sage: is_comparability(graphs.CompleteGraph(5), certificate=True)
        (True, Digraph on 5 vertices)

    TESTS:

    Let us ensure that no exception is raised when we go over all small graphs::

        sage: from sage.graphs.comparability import is_comparability
        sage: [len([g for g in graphs(i) if is_comparability(g, certificate=True)[0]]) for i in range(7)]
        [1, 1, 2, 4, 11, 33, 144]
    """
    g._scream_if_not_simple()
    if not g.size():
        if certificate:
            from sage.graphs.digraph import DiGraph
            return True, DiGraph(g)
        return True

    if algorithm == "greedy":
        comparability_test = greedy_is_comparability_with_certificate(g, certificate=certificate)
    elif algorithm == "MILP":
        comparability_test = is_comparability_MILP(g, certificate=certificate,
                                                   solver=solver, verbose=verbose)

    if not certificate:
        return comparability_test

    # Checking that the orientation found is indeed transitive. No
    # reason why it should not, but no reason why we should not check
    # anyway :-p
    isit, certif = comparability_test

    if check and isit and (not certif.is_transitive()):
        raise ValueError("Looks like there is a bug somewhere. The "
                         "algorithm thinks that the orientation is "
                         "transitive, but we just checked and it is not. "
                         "Please report the bug on sage-devel, and give "
                         "us the graph that made this method fail !")

    return isit, certif


def is_permutation(g, algorithm='greedy', certificate=False, check=True,
                   solver=None, verbose=0):
    r"""
    Check whether the graph is a permutation graph.

    For more information on permutation graphs, refer to the documentation of
    the :mod:`comparability module <sage.graphs.comparability>`.

    INPUT:

    - ``g`` -- a graph

    - ``algorithm`` -- string (default: ``'greedy'``); choose the implementation
      used for the subcalls to :meth:`is_comparability`

      - ``'greedy'`` -- a greedy algorithm (see the documentation of the
        :mod:`comparability module <sage.graphs.comparability>`)

      - ``'MILP'`` -- a Mixed Integer Linear Program formulation of the
        problem. Beware, for this implementation is unable to return negative
        certificates ! When ``certificate = True``, negative certificates are
        always equal to ``None``. ``True`` certificates are valid, though.

    - ``certificate`` -- boolean (default: ``False``); whether to return a
      certificate for the answer given. For ``True`` answers the certificate is
      a permutation, for ``False`` answers it is a no-certificate for the test
      of comparability or co-comparability.

    - ``check`` -- boolean (default: ``True``); whether to check that the
      permutations returned indeed create the expected Permutation graph. Pretty
      cheap compared to the rest, hence a good investment. It is enabled by
      default.

    - ``solver`` -- string (default: ``None``); specifies a Mixed Integer Linear
      Programming (MILP) solver to be used. If set to ``None``, the default one
      is used. For more information on MILP solvers and which default solver is
      used, see the method :meth:`solve
      <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
      :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- integer (default: 0); sets the level of verbosity. Set
      to 0 by default, which means quiet.

    .. NOTE::

        As the ``True`` certificate is a :class:`Permutation` object, the
        segment intersection model of the permutation graph can be visualized
        through a call to :meth:`Permutation.show
        <sage.combinat.permutation.Permutation.show>`.

    EXAMPLES:

    A permutation realizing the bull graph::

        sage: from sage.graphs.comparability import is_permutation
        sage: g = graphs.BullGraph()
        sage: _ , certif = is_permutation(g, certificate=True)
        sage: h = graphs.PermutationGraph(*certif)
        sage: h.is_isomorphic(g)
        True

    Plotting the realization as an intersection graph of segments::

        sage: true, perm = is_permutation(g, certificate=True)
        sage: p1 = Permutation([nn+1 for nn in perm[0]])
        sage: p2 = Permutation([nn+1 for nn in perm[1]])
        sage: p = p2 * p1.inverse()
        sage: p.show(representation='braid')                                            # needs sage.plot

    TESTS:

    Trying random permutations, first with the greedy algorithm::

        sage: from sage.graphs.comparability import is_permutation
        sage: for i in range(20):
        ....:     p = Permutations(10).random_element()
        ....:     g1 = graphs.PermutationGraph(p)
        ....:     isit, certif = is_permutation(g1, certificate=True)
        ....:     if not isit:
        ....:        print("Something is wrong here !!")
        ....:        break
        ....:     g2 = graphs.PermutationGraph(*certif)
        ....:     if not g1.is_isomorphic(g2):
        ....:        print("Something is wrong here !!")
        ....:        break

    Then with MILP::

        sage: from sage.graphs.comparability import is_permutation
        sage: for i in range(20):                                                       # needs sage.numerical.mip
        ....:     p = Permutations(10).random_element()
        ....:     g1 = graphs.PermutationGraph(p)
        ....:     isit, certif = is_permutation(g1, algorithm='MILP', certificate=True)
        ....:     if not isit:
        ....:        print("Something is wrong here !!")
        ....:        break
        ....:     g2 = graphs.PermutationGraph(*certif)
        ....:     if not g1.is_isomorphic(g2):
        ....:        print("Something is wrong here !!")
        ....:        break
    """
    if not certificate:
        # No certificate... A piece of cake
        return (is_comparability(g, algorithm=algorithm, solver=solver, verbose=verbose) and
                is_comparability(g.complement(), algorithm=algorithm, solver=solver, verbose=verbose))

    # First poset, we stop if it fails
    isit, certif = is_comparability(g, algorithm=algorithm, certificate=True,
                                    solver=solver, verbose=verbose)
    if not isit:
        return False, certif

    # Second poset
    isit, co_certif = is_comparability(g.complement(), algorithm=algorithm, certificate=True,
                                       solver=solver, verbose=verbose)
    if not isit:
        return False, co_certif

    # Building the two orderings
    tmp = list(co_certif.edges(labels=False, sort=False))
    for u, v in certif.edge_iterator(labels=False):
        co_certif.add_edge(v, u)
    certif.add_edges(tmp)

    ordering = certif.topological_sort()
    co_ordering = co_certif.topological_sort()

    # Try to build the Permutation graph from the permutations, just to make
    # sure nothing weird happened !
    if check:
        from sage.graphs.graph_generators import GraphGenerators
        pg = GraphGenerators().PermutationGraph(ordering, co_ordering)
        if not pg.is_isomorphic(g):
            raise ValueError("There is a mistake somewhere ! It looks like "
                             "the Permutation Graph model computed does "
                             "not match the input graph !")

    return True, (ordering, co_ordering)


def is_transitive(g, certificate=False):
    r"""
    Test whether the digraph is transitive.

    A digraph is transitive if for any pair of vertices `u,v\in G` linked by a
    `uv`-path the edge `uv` belongs to `G`.

    INPUT:

    - ``certificate`` -- whether to return a certificate for negative answers

      - If ``certificate = False`` (default), this method returns ``True`` or
        ``False`` according to the graph.

      - If ``certificate = True``, this method either returns ``True`` answers
        or yield a pair of vertices `uv` such that there exists a `uv`-path in
        `G` but `uv\not\in G`.

    EXAMPLES::

        sage: digraphs.Circuit(4).is_transitive()
        False
        sage: digraphs.Circuit(4).is_transitive(certificate=True)
        (0, 2)
        sage: digraphs.RandomDirectedGNP(30,.2).is_transitive()
        False
        sage: D = digraphs.DeBruijn(5, 2)                                               # needs sage.combinat
        sage: D.is_transitive()                                                         # needs sage.combinat
        False
        sage: cert = D.is_transitive(certificate=True)                                  # needs sage.combinat
        sage: D.has_edge(*cert)                                                         # needs sage.combinat
        False
        sage: bool(D.shortest_path(*cert))                                              # needs sage.combinat
        True
        sage: digraphs.RandomDirectedGNP(20,.2).transitive_closure().is_transitive()    # needs networkx
        True
    """
    cdef int n = g.order()

    if n <= 2:
        return True

    cdef list int_to_vertex = list(g)
    cdef unsigned short * distances = c_distances_all_pairs(g, vertex_list=int_to_vertex)
    cdef unsigned short * c_distances = distances

    cdef int i, j

    # Only 3 distances can appear in the matrix of all distances : 0, 1, and
    # infinity. Anything else is a proof of nontransitivity !

    for j in range(n):
        for i in range(n):
            if c_distances[i] != <unsigned short> -1 and c_distances[i] > 1:
                sig_free(distances)
                if certificate:

                    return int_to_vertex[j], int_to_vertex[i]
                else:
                    return False

        c_distances += n

    sig_free(distances)
    return True
