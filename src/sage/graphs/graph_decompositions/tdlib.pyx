# distutils: language = c++
# sage_setup: distribution = sagemath-tdlib

r"""
Interface with TdLib (algorithms for tree decompositions)

This module defines functions based on TdLib, a library that implements
algorithms for tree decompositions written by Lukas Larisch.

**Definition** :

A `tree decomposition` of a graph `G` is a pair `(T, \beta)` consisting of a
tree T and a function `\beta: V(T) \rightarrow 2^{V(G)}` associating with each
node `t \in V(T)` a set of vertices `\beta (t) \subseteq V(G)` such that

* (T1) for every edge `e \in E(G)` there is a node `t \in V(T)` with `e
  \subseteq \beta (t)`, and

* (T2) for all `v \in V(G)` the set `\beta^{-1} := \{t \in V(T): v \in \beta
  (t)\}` is non-empty and connected in T.

The width of `(T, \beta)` is defined as `max\{|\beta (t)|-1: t \in V(T) \}`.
The treewidth of G is defined as the minimum width over all tree decompositions
of `G`.

**Some known results** :

- Trees have treewidth 1

- Cycles have treewidth 2

- Series-parallel graphs have treewidth at most 2

- Cliques must be contained in some bag of a tree decomposition

Computing the treewidth or a tree decomposition of a given graph is NP-hard in
general.

**This module contains the following functions** :

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`treedecomposition_exact` | Compute a tree decomposition of exact width
    :meth:`get_width` | Return the width of a given tree decomposition


AUTHOR: Lukas Larisch (10-25-2015): Initial version

REFERENCE:

- [ST1993]_

- [AP1986]_

- [Bod1993]_

Methods
-------
"""

from libcpp.vector cimport vector
from cysignals.signals cimport sig_on, sig_off

from sage.sets.set import Set
from sage.graphs.graph import Graph

cdef extern from "sage_tdlib.cpp":
    int sage_exact_decomposition(vector[unsigned int] &V_G, vector[unsigned int] &E_G, vector[vector[int]] &V_T, vector[unsigned int] &E_T, int lb)


##############################################################
# ########## GRAPH/DECOMPOSITION ENCODING/DECODING ###########
# the following will be used implicitly to do the translation
# between Sage graph encoding and BGL graph encoding.

cdef make_tdlib_graph(G, vertex_to_int, vector[unsigned int] &V, vector[unsigned int] &E):
    for i in range(G.order()):
        V.push_back(i)

    for u, v in G.edge_iterator(labels=False):
        E.push_back(vertex_to_int[u])
        E.push_back(vertex_to_int[v])


cdef make_sage_decomp(G, vector[vector[int]] &V, vector[unsigned int] &E, int_to_vertex):
    cdef int i, j
    for i in range(len(V)):
        G.add_vertex(Set([int_to_vertex[j] for j in V[i]]))

    for i in range(0, len(E), 2):
        G.add_edge(Set([int_to_vertex[j] for j in V[E[i]]]), Set([int_to_vertex[j] for j in V[E[i+1]]]))


##############################################################
# ########## EXACT ALGORITHMS ################################

def treedecomposition_exact(G, lb=-1):
    r"""
    Compute a tree decomposition of exact width.

    The returned tree decomposition is exact iff the given lower bound is not
    greater than the treewidth of the input graph. Otherwise a tree
    decomposition of a width than matches the given lower bound will be
    computed.

    INPUT:

    - ``G`` -- a generic graph

    - ``lb`` -- integer (default: -1); a lower bound to the treewidth of G,
      e.g. computed by lower_bound

    OUTPUT:

    A tree decomposition of ``G`` of ``tw(G)``, if the lower bound was not
    greater than ``tw(G)``, otherwise a tree decomposition of ``width = lb``.

    .. WARNING::

        The computation can take a lot of time for a graph `G` on more than
        about 30 vertices and `tw(G) > 3`.

    EXAMPLES::

        sage: # optional - tdlib
        sage: import sage.graphs.graph_decompositions.tdlib as tdlib
        sage: G = graphs.HouseGraph()
        sage: T = tdlib.treedecomposition_exact(G)
        sage: T.show(vertex_size=2000)

    TESTS::

        sage: # optional - tdlib
        sage: import sage.graphs.graph_decompositions.tdlib as tdlib
        sage: G = graphs.HouseGraph()
        sage: T = tdlib.treedecomposition_exact(G)
        sage: G = graphs.PetersenGraph()
        sage: T = tdlib.treedecomposition_exact(G)

    """
    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cdef list int_to_vertex = list(G)
    cdef dict vertex_to_int = {v: i for i, v in enumerate(G)}

    make_tdlib_graph(G, vertex_to_int, V_G, E_G)

    cdef int c_lb = lb

    sig_on()
    sage_exact_decomposition(V_G, E_G, V_T, E_T, c_lb)
    sig_off()

    T = Graph(name="Tree decomposition")
    make_sage_decomp(T, V_T, E_T, int_to_vertex)

    return T


def get_width(T):
    """
    Return the width of a given tree decomposition.

    The width of a tree decompositions is the maximal size of a bag minus one.

    INPUT:

    - ``T`` -- a tree decomposition

    OUTPUT:

    - The width of ``T``

    EXAMPLES::

        sage: # optional - tdlib
        sage: import sage.graphs.graph_decompositions.tdlib as tdlib
        sage: G = graphs.PetersenGraph()
        sage: T = tdlib.treedecomposition_exact(G)
        sage: tdlib.get_width(T)
        4
    """
    return (max(len(x) for x in T) - 1) if T else -1
