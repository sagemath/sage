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

**This module containes the following functions** :

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`treedecomposition_exact` | Computes a tree decomposition of exact width
    :meth:`get_width` | Returns the width of a given tree decomposition


AUTHOR: Lukas Larisch (10-25-2015): Initial version

REFERENCE:

.. [ST93] \P. D. Seymour and Robin Thomas,
   Graph searching and a min-max theorem for tree-width,
   J. Comb. Theory Ser. B 58, 1 (May 1993), 22-33.

.. [AP86] \S. Arnborg, A. Proskurowski,
   Characterization and Recognition of Partial 3-Trees,
   SIAM Journal of Alg. and Discrete Methods,
   Vol. 7, pp. 305-314, 1986

.. [Bodlaender93] \H. L. Bodlaender,
   A Tourist Guide through Treewidth, Acta Cybern. 1993

Methods
-------
"""

from libcpp.vector cimport vector
from cysignals.signals cimport sig_on, sig_off

from sage.sets.set import Set
from sage.graphs.graph import Graph

cdef extern from "tdlib/sage_tdlib.cpp":
     int sage_exact_decomposition(vector[unsigned int] &V_G, vector[unsigned int] &E_G, vector[vector[int]] &V_T, vector[unsigned int] &E_T, int lb)

##############################################################
############ GRAPH/DECOMPOSITION ENCODING/DECODING ###########
#the following will be used implicitly do the translation
#between Sage graph encoding and BGL graph encoding.

cdef make_tdlib_graph(G, vector[unsigned int] &V, vector[unsigned int] &E):
    for v in G.vertices():
        V.push_back(v)

    for u,v in G.edges(labels=False):
        E.push_back(u)
        E.push_back(v)

cdef make_sage_decomp(G, vector[vector[int]] &V, vector[unsigned int] &E, label_map):
    for i in range(0, len(V)):
        G.add_vertex(Set([label_map[j] for j in V[i]]))

    for i in range(0, len(E), 2):
        G.add_edge(Set(V[E[i]]), Set(V[E[i+1]]))

##############################################################
############ EXACT ALGORITHMS ################################

def treedecomposition_exact(G, lb=-1):
    r"""
    Computes a tree decomposition of exact width, iff the given lower bound
    is not greater than the treewidth of the input graph. Otherwise
    a tree decomposition of a width than matches the given lower bound
    will be computed.

    INPUT:

    - ``G`` -- a generic graph

    - ``lb`` -- a lower bound to the treewidth of G, e.g. computed by lower_bound (default: ``'-1'``)

    OUTPUT:

    A tree decomposition of ``G`` of ``tw(G)``, if the lower bound was not
    greater than ``tw(G)``, otherwise a tree decomposition of ``width = lb``.

..  WARNING::

    The computation can take a lot of time for a graph `G` on more than about 30
    vertices and `tw(G) > 3`.

    EXAMPLES::

        sage: import sage.graphs.graph_decompositions.tdlib as tdlib # optional - tdlib
        sage: G = graphs.HouseGraph()                                # optional - tdlib
        sage: T = tdlib.treedecomposition_exact(G)                   # optional - tdlib
        sage: T.show(vertex_size=2000)                               # optional - tdlib

    TESTS::

        sage: import sage.graphs.graph_decompositions.tdlib as tdlib # optional - tdlib
        sage: G = graphs.HouseGraph()                                # optional - tdlib
        sage: T = tdlib.treedecomposition_exact(G)                   # optional - tdlib
        sage: G = graphs.PetersenGraph()                             # optional - tdlib
        sage: T = tdlib.treedecomposition_exact(G)                   # optional - tdlib

    """
    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    V = G.vertices()
    G_int = G.relabel(inplace=False)

    make_tdlib_graph(G_int, V_G, E_G)

    cdef int c_lb = lb

    sig_on()
    sage_exact_decomposition(V_G, E_G, V_T, E_T, c_lb)
    sig_off()

    T = Graph()
    T.name("Tree decomposition")
    make_sage_decomp(T, V_T, E_T, V)

    return T

def get_width(T):
    """
    Returns the width (maximal size of a bag minus one) of a given tree decomposition.

    INPUT:

    - ``T`` -- a tree decomposition

    OUTPUT:

    - The width of ``T``

    EXAMPLES::
        sage: import sage.graphs.graph_decompositions.tdlib as tdlib # optional - tdlib
        sage: G = graphs.PetersenGraph()                             # optional - tdlib
        sage: T = tdlib.treedecomposition_exact(G)                   # optional - tdlib
        sage: tdlib.get_width(T)                                     # optional - tdlib
        4
    """
    return max(len(x) for x in T)-1 if len(T) > 0 else -1
