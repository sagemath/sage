# cython: binding=True
r"""
Products of graphs

This module gathers everything related to graph products. At the moment it
contains an implementation of a recognition algorithm for graphs that can be
written as a Cartesian product of smaller ones.

Author:

- Nathann Cohen (May 2012 -- coded while watching the election of Francois
  Hollande on TV)

Cartesian product of graphs -- the recognition problem
------------------------------------------------------

First, a definition:

  **Definition** The Cartesian product of two graphs `G` and `H`, denoted
  `G\square H`, is a graph defined on the pairs `(g, h)\in V(G)\times V(H)`.

  Two elements `(g, h),(g', h')\in V(G\square H)` are adjacent in `G\square H`
  if and only if :

  - `g=g'` and `hh'\in H`; or
  - `h=h'` and `gg'\in G`

Two remarks follow :

#. The Cartesian product is commutative

#. Any edge `uv` of a graph `G_1 \square \cdots \square G_k` can be given a
   color `i` corresponding to the unique index `i` such that `u_i` and `v_i`
   differ.

The problem that is of interest to us in the present module is the following:

  **Recognition problem** Given a graph `G`, can we guess whether there exist
  graphs `G_1, ..., G_k` such that `G=G_1\square \cdots \square G_k` ?

This problem can actually be solved, and the resulting factorization is
unique. What is explained below can be found in the book *Handbook of Product
Graphs* [HIK2011]_.

Everything is actually based on simple observations. Given a graph `G`, finding
out whether `G` can be written as the product of several graphs can be attempted
by trying to color its edges according to some rules. Indeed, if we are to color
the edges of `G` in such a way that each color class represents a factor of `G`,
we must ensure several things.

  **Remark 1** In any cycle of `G` no color can appear exactly once.

  Indeed, if only one edge `uv` of a cycle were labelled with color `i`, it
  would mean that:

  #. The only difference between `u` and `v` lies in their `i` th coordinate

  #. It is possible to go from `u` to `v` by changing only coordinates
     different from the `i` th

  A contradiction indeed.

  .. image:: ../../../media/cycle.png

  That means that, for instance, the edges of a triangle necessarily have the
  same color.

  **Remark 2** If two consecutive edges `u_1u_2` and `u_2u_3` have different
  colors, there necessarily exists a unique vertex `u_4` different from `u_2`
  and incident to both `u_1` and `u_3`.

  In this situation, opposed edges necessarily have the same colors because of
  the previous remark.

  .. image:: ../../../media/square.png

  **1st criterion** : As a corollary, we know that:

  #. If two vertices `u,v` have a *unique* common neighbor `x`, then `ux` and
     `xv` have the same color.

  #. If two vertices `u, v` have more that two common neighbors `x_1, ...,
     x_k` then all edges between the `x_i` and the vertices of `u,v` have the
     same color. This is also a consequence of the first remark.

  **2nd criterion** : if two edges `uv` and `u'v'` of the product graph
  `G\square H` are such that `d(u,u')+d(v,v')\neq d(u,v') + d(v,u')` then the
  two edges `uv` and `u'v'` necessarily have the same color.

    This is a consequence of the fact that for any two vertices `u,v` of
    `G\square H` (where `u=(u_G,u_H)` and `v=(v_G,v_H)`), we have `d(u,v) =
    d_G(u_G,v_G)+d_H(u_H,v_H)`. Indeed, a shortest path from `u` to `v` in
    `G\square H` contains the information of a shortest path from `u_G` to `v_G`
    in `G`, and a shortest path from `u_H` to `v_H` in `H`.

The algorithm
^^^^^^^^^^^^^

The previous remarks tell us that some edges are in some way equivalent to some
others, i.e. that their colors are equal. In order to compute the coloring we
are looking for, we therefore build a graph on the *edges* of a graph `G`,
linking two edges whenever they are found to be equivalent according to the
previous remarks.

All that is left to do is to compute the connected components of this new graph,
as each of them representing the edges of a factor. Of course, only one
connected component indicates that the graph has no factorization.

Then again, please refer to [HIK2011]_ for any technical question.

To Do
^^^^^

This implementation is made at Python level, and some parts of the algorithm
could be rewritten in Cython to save time. Especially when enumerating all pairs
of edges and computing their distances. This can easily be done in C with the
functions from the :mod:`sage.graphs.distances_all_pairs` module.

Methods
-------
"""

#******************************************************************************
#          Copyright (C) 2012 Nathann Cohen <nathann.cohen@gmail.com>         *
#                                                                             *
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)*
#                         http://www.gnu.org/licenses/                        *
#******************************************************************************

from copy import copy


def is_cartesian_product(g, certificate=False, relabeling=False):
    r"""
    Test whether the graph is a Cartesian product.

    INPUT:

    - ``certificate`` -- boolean (default: ``False``); if ``certificate =
      False`` (default) the method only returns ``True`` or ``False``
      answers. If ``certificate = True``, the ``True`` answers are replaced by
      the list of the factors of the graph.

    - ``relabeling`` -- boolean (default: ``False``); if ``relabeling = True``
      (implies ``certificate = True``), the method also returns a dictionary
      associating to each vertex its natural coordinates as a vertex of a
      product graph. If `g` is not a Cartesian product, ``None`` is returned
      instead.

    .. SEEALSO::

        - :meth:`sage.graphs.generic_graph.GenericGraph.cartesian_product`

        - :mod:`~sage.graphs.graph_decompositions.graph_products` -- a module on
          graph products.

    .. NOTE::

        This algorithm may run faster whenever the graph's vertices are integers
        (see :meth:`~sage.graphs.generic_graph.GenericGraph.relabel`). Give it a
        try if it is too slow !

    EXAMPLES:

    The Petersen graph is prime::

        sage: from sage.graphs.graph_decompositions.graph_products import is_cartesian_product
        sage: g = graphs.PetersenGraph()
        sage: is_cartesian_product(g)
        False

    A 2d grid is the product of paths::

        sage: g = graphs.Grid2dGraph(5,5)
        sage: p1, p2 = is_cartesian_product(g, certificate = True)
        sage: p1.is_isomorphic(graphs.PathGraph(5))
        True
        sage: p2.is_isomorphic(graphs.PathGraph(5))
        True

    Forgetting the graph's labels, then finding them back::

        sage: g.relabel()
        sage: b,D = g.is_cartesian_product(g, relabeling=True)
        sage: b
        True
        sage: D  # random isomorphism
        {0: (20, 0), 1: (20, 1), 2: (20, 2), 3: (20, 3), 4: (20, 4),
         5: (15, 0), 6: (15, 1), 7: (15, 2), 8: (15, 3), 9: (15, 4),
         10: (10, 0), 11: (10, 1), 12: (10, 2), 13: (10, 3), 14: (10, 4),
         15: (5, 0), 16: (5, 1), 17: (5, 2), 18: (5, 3), 19: (5, 4),
         20: (0, 0), 21: (0, 1), 22: (0, 2), 23: (0, 3), 24: (0, 4)}

    And of course, we find the factors back when we build a graph from a
    product::

        sage: g = graphs.PetersenGraph().cartesian_product(graphs.CycleGraph(3))
        sage: g1, g2 = is_cartesian_product(g, certificate = True)
        sage: any( x.is_isomorphic(graphs.PetersenGraph()) for x in [g1,g2])
        True
        sage: any( x.is_isomorphic(graphs.CycleGraph(3)) for x in [g1,g2])
        True

    TESTS:

    Wagner's Graph (:trac:`13599`)::

        sage: g = graphs.WagnerGraph()
        sage: g.is_cartesian_product()
        False

    Empty and one-element graph (:trac:`19546`)::

        sage: Graph().is_cartesian_product()
        False
        sage: Graph({0:[]}).is_cartesian_product()
        False
    """
    g._scream_if_not_simple()
    if g.is_directed():
        raise NotImplementedError("recognition of Cartesian product is not implemented for directed graphs")
    if relabeling:
        certificate = True

    from sage.rings.integer import Integer

    if not g.is_connected():
        raise NotImplementedError("recognition of Cartesian product is not implemented for disconnected graphs")

    # Of course the number of vertices of g cannot be prime !
    if g.order() <= 3 or Integer(g.order()).is_prime():
        return (False, None) if relabeling else False

    from sage.graphs.graph import Graph

    # As we need the vertices of g to be linearly ordered, we copy the graph and
    # relabel it
    cdef list int_to_vertex = list(g)
    cdef dict vertex_to_int = {vert: i for i, vert in enumerate(int_to_vertex)}
    g_int = g.relabel(perm=vertex_to_int, inplace=False)

    # Reorder the vertices of an edge
    r = lambda x,y: (x, y) if x < y else (y, x)

    cdef int x, y, u, v
    cdef set un, intersect

    # The equivalence graph on the edges of g
    h = Graph()
    h.add_vertices(r(x,y) for x, y in g_int.edge_iterator(labels=False))

    # For all pairs of vertices u,v of G, according to their number of common
    # neighbors... See the module's documentation !
    for u in g_int:
        un = set(g_int.neighbor_iterator(u))
        for v in g_int.breadth_first_search(u):

            # u and v are different
            if u == v:
                continue

            # List of common neighbors
            intersect = un & set(g_int.neighbor_iterator(v))

            # If u and v have no neighbors and uv is not an edge then their
            # distance is at least 3. As we enumerate the vertices in a
            # breadth-first search, it means that we already checked all the
            # vertices at distance less than two from u, and we are done with
            # this loop !
            if not intersect:
                if g_int.has_edge(u, v):
                    continue
                else:
                    break

            # If uv is an edge
            if g_int.has_edge(u, v):
                h.add_path([r(u, x) for x in intersect] + [r(v, x) for x in intersect])

            # Only one common neighbor
            elif len(intersect) == 1:
                x = intersect.pop()
                h.add_edge(r(u, x), r(v, x))

            # Exactly 2 neighbors
            elif len(intersect) == 2:
                x, y = intersect
                h.add_edge(r(u, x), r(v, y))
                h.add_edge(r(v, x), r(u, y))
            # More
            else:
                h.add_path([r(u,x) for x in intersect] + [r(v,x) for x in intersect])

    # Edges uv and u'v' such that d(u,u')+d(v,v') != d(u,v')+d(v,u') are also
    # equivalent

    cdef list edges = list(g_int.edges(labels=False, sort=False))
    cdef dict d = g_int.distance_all_pairs()
    cdef int uu, vv
    for i, (u, v) in enumerate(edges):
        du = d[u]
        dv = d[v]
        for j in range(i + 1, g_int.size()):
            uu, vv = edges[j]
            if du[uu] + dv[vv] != du[vv] + dv[uu]:
                h.add_edge(r(u, v), r(uu, vv))

    # Gathering the connected components, relabeling the vertices on-the-fly
    edges = [[(int_to_vertex[u], int_to_vertex[v]) for u, v in cc]
             for cc in h.connected_components()]

    #Print the graph, distinguishing the edges according to their color classes
    #
    #from sage.plot.colors import rainbow
    #g.show(edge_colors = dict(zip(rainbow(len(edges)),edges)))

    # Only one connected component ?
    if len(edges) == 1:
        return (False, None) if relabeling else False

    # Building the list of factors
    cdef list factors = []
    for cc in edges:
        tmp = Graph()
        tmp.add_edges(cc)
        factors.append(tmp.subgraph(vertices = tmp.connected_components()[0]))

    # Computing the product of these graphs
    answer = factors[0]
    for i in range(1,len(factors)):
        answer = answer.cartesian_product(factors[i])

    # Checking that the resulting graph is indeed isomorphic to what we have.
    isiso, dictt = g.is_isomorphic(answer, certificate=True)
    if not isiso:
        raise ValueError("something weird happened during the algorithm... "
                         "Please report the bug and give us the graph instance"
                         " that made it fail !")
    if relabeling:
        return isiso, dictt
    if certificate:
        return factors
    else:
        return True
