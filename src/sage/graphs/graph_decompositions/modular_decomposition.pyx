# distutils: language = c++
# distutils: extra_compile_args = -std=c++11
r"""
Modular Decomposition

This module implements the function for computing the modular decomposition
of undirected graphs.

AUTHORS:

- Lokesh Jain (2017): first implementation of the linear time algorithm of
  D. Corneil, M. Habib, C. Paul and M. Tedder [TCHP2008]_

- David Einstein (2018): added the algorithm of M. Habib and M. Maurer [HM1979]_

- Cyril Bouvier (2024): second implementation of the linear time algorithm
  of D. Corneil, M. Habib, C. Paul and M. Tedder [TCHP2008]_
"""
# ****************************************************************************
#       Copyright (C) 2017 Lokesh Jain <lokeshj1703@gmail.com>
#                     2018 David Einstein <deinst@gmail.com>
#                     2024 Cyril Bouvier <cyril.bouvier@lirmm.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from cython.operator cimport dereference as deref

from enum import IntEnum

from sage.graphs.base.c_graph cimport CGraph, CGraphBackend
from sage.graphs.graph_decompositions.slice_decomposition cimport \
        extended_lex_BFS
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.misc.lazy_import import lazy_import
from sage.misc.random_testing import random_testing


################################################################################
#                     Corneil-Habib-Paul-Tedder algorithm                      #
################################################################################
def corneil_habib_paul_tedder_algorithm(G):
    r"""
    Compute the modular decomposition by the algorithm of Corneil, Habib, Paul
    and Tedder.

    INPUT:

    - ``G`` -- the graph for which modular decomposition tree needs to be
      computed

    OUTPUT: an object of type Node representing the modular decomposition tree
    of the graph G

    This function compute the modular decomposition of the given graph by the
    algorithm of Corneil, Habib, Paul and Tedder [TCHP2008]_. It is a recursive,
    linear-time algorithm that first computes the slice decomposition of the
    graph (via the extended lexBFS algorithm) and then computes the modular
    decomposition by calling itself recursively on the slices of the previously
    computed slice decomposition.

    .. SEEALSO::

        * :mod:`~sage.graphs.graph_decompositions.slice_decomposition` --
          compute a slice decomposition of the simple undirect graph

    This function should not be used directly, it should be called via the
    ``modular_decomposition`` method of ``Graph`` with the parameter
    ``algorithm='corneil_habib_paul_tedder'``.

    This functions assumes that ``graph`` is a object of the class ``Graph`` and
    is a simple graph.

    TESTS::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: recreate_decomposition(15, corneil_habib_paul_tedder_algorithm,
        ....:                                                         3, 4, 0.2)
        sage: recreate_decomposition(10, corneil_habib_paul_tedder_algorithm,
        ....:                                                         4, 5, 0.2)
        sage: recreate_decomposition(3, corneil_habib_paul_tedder_algorithm,
        ....:                                                         6, 5, 0.2)
    """
    cdef CGraphBackend Gbackend = <CGraphBackend> G._backend
    cdef CGraph cg = Gbackend.cg()

    cdef vector[int] sigma
    cdef vector[vector[int]] lex_label
    cdef vector[size_t] xslice_len

    # Compute the slice decomposition using the extended lexBFS algorithm
    extended_lex_BFS(cg, sigma, NULL, -1, NULL, &xslice_len, &lex_label)

    cdef SDData SD
    SD.set_from_data(0, sigma.data(), xslice_len.data(), lex_label.data())

    MD = corneil_habib_paul_tedder_inner(SD)

    r = md_tree_node_to_md_tree(MD, Gbackend)
    dealloc_md_tree_nodes_recursively(MD)
    return r


cdef object _md_tree_node_to_md_tree_inner_rec(const md_tree_node *n,
                                               CGraphBackend Gb):
    """
    Utility function for :func:`md_tree_node_to_md_tree`.
    """
    cdef md_tree_node *c
    if deref(n).is_leaf():
        return Node.create_leaf(Gb.vertex_label(deref(n).vertex))
    else:
        if deref(n).is_series():
            node = Node(NodeType.SERIES)
        elif deref(n).is_parallel():
            node = Node(NodeType.PARALLEL)
        else:  # is_prime
            node = Node(NodeType.PRIME)
        node.children.extend(
            _md_tree_node_to_md_tree_inner_rec(c, Gb)
                                                  for c in deref(n).children)
        return node


cdef object md_tree_node_to_md_tree(const md_tree_node *n, CGraphBackend Gb):
    """
    This function converts a modular decomposition tree (given as a pointer to a
    md_tree_node) into an object of the Python class Node (which is then
    converted into the correct output by :meth:`Graph.modular_decomposition`).

    The graph backend is needed to convert the int stored into the md_tree_node
    into the corresponding vertex of the Python graph.

    This function deals with the case of an empty tree and then delegates the
    actual conversion to :func:`_md_tree_node_to_md_tree_inner_rec`.

    TESTS:

    Indirect doctests::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: corneil_habib_paul_tedder_algorithm(Graph(1))
        NORMAL [0]
        sage: corneil_habib_paul_tedder_algorithm(Graph(2))
        PARALLEL [NORMAL [0], NORMAL [1]]
        sage: corneil_habib_paul_tedder_algorithm(graphs.CompleteGraph(3))
        SERIES [NORMAL [1], NORMAL [2], NORMAL [0]]
    """
    if n == NULL:
        return Node(NodeType.EMPTY)
    else:
        return _md_tree_node_to_md_tree_inner_rec(n, Gb)


################################################################################
class NodeType(IntEnum):
    """
    NodeType is an enumeration class used to define the various types of nodes
    in modular decomposition tree.

    The various node types defined are

    - ``PARALLEL`` -- indicates the node is a parallel module

    - ``SERIES`` -- indicates the node is a series module

    - ``PRIME`` -- indicates the node is a prime module

    - ``EMPTY`` -- indicates a empty tree

    - ``NORMAL`` -- indicates the node is normal containing a vertex
    """
    PRIME = 0
    SERIES = 1
    PARALLEL = 2
    NORMAL = 3
    EMPTY = -1

    def __repr__(self) -> str:
        r"""
        Return a string representation of a ``NodeType`` object.

        TESTS::

            sage: from sage.graphs.graph_decompositions.modular_decomposition import NodeType
            sage: repr(NodeType.PARALLEL)
            'PARALLEL'
            sage: str(NodeType.PRIME)
            'PRIME'
        """
        return self.name

    __str__ = __repr__


class Node:
    """
    Node class stores information about the node type.

    Node type can be ``PRIME``, ``SERIES``, ``PARALLEL``, ``NORMAL`` or
    ``EMPTY``.

    - ``node_type`` -- is of type NodeType and specifies the type of node
    """
    def __init__(self, node_type):
        r"""
        Create a node with the given node type.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.modular_decomposition import *
            sage: n = Node(NodeType.SERIES); n.node_type
            SERIES
            sage: n.children
            []
        """
        self.node_type = node_type
        self.children = []

    def is_prime(self):
        r"""
        Check whether ``self`` is a prime node.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.modular_decomposition import *
            sage: n = Node(NodeType.PRIME)
            sage: n.children.append(Node.create_leaf(1))
            sage: n.children.append(Node.create_leaf(2))
            sage: n.is_prime()
            True
            sage: (n.children[0].is_prime(), n.children[1].is_prime())
            (False, False)
        """
        return self.node_type == NodeType.PRIME

    def is_series(self):
        r"""
        Check whether ``self`` is a series node.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.modular_decomposition import *
            sage: n = Node(NodeType.SERIES)
            sage: n.children.append(Node.create_leaf(1))
            sage: n.children.append(Node.create_leaf(2))
            sage: n.is_series()
            True
            sage: (n.children[0].is_series(), n.children[1].is_series())
            (False, False)
        """
        return self.node_type == NodeType.SERIES

    def is_empty(self):
        r"""
        Check whether ``self`` is an empty node.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.modular_decomposition import *
            sage: Node(NodeType.EMPTY).is_empty()
            True
            sage: Node.create_leaf(1).is_empty()
            False
        """
        return self.node_type == NodeType.EMPTY

    def is_leaf(self):
        r"""
        Check whether ``self`` is a leaf.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.modular_decomposition import *
            sage: n = Node(NodeType.PRIME)
            sage: n.children.append(Node.create_leaf(1))
            sage: n.children.append(Node.create_leaf(2))
            sage: n.is_leaf()
            False
            sage: all(c.is_leaf() for c in n.children)
            True
        """
        return self.node_type == NodeType.NORMAL

    def __repr__(self):
        r"""
        Return a string representation of the node.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.modular_decomposition import *
            sage: n = Node(NodeType.PRIME)
            sage: n.children.append(Node.create_leaf(1))
            sage: n.children.append(Node.create_leaf(2))
            sage: str(n)
            'PRIME [NORMAL [1], NORMAL [2]]'
        """
        return f"{self.node_type} {self.children}"

    def __eq__(self, other):
        r"""
        Compare two nodes for equality.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.modular_decomposition import *
            sage: n1 = Node(NodeType.PRIME)
            sage: n2 = Node(NodeType.PRIME)
            sage: n3 = Node(NodeType.SERIES)
            sage: n1 == n2
            True
            sage: n1 == n3
            False
        """
        return (self.node_type == other.node_type and
                self.children == other.children)

    @classmethod
    def create_leaf(cls, v):
        """
        Return Node object that is a leaf corresponding to the vertex ``v``

        INPUT:

        - ``vertex`` -- vertex number

        OUTPUT: a node object representing the vertex with ``node_type`` set as
        ``NodeType.NORMAL``

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.modular_decomposition import Node
            sage: node = Node.create_leaf(2)
            sage: node
            NORMAL [2]
        """
        node = cls(NodeType.NORMAL)
        node.children.append(v)
        return node


def print_md_tree(root):
    """
    Print the modular decomposition tree.

    INPUT:

    - ``root`` -- root of the modular decomposition tree

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: print_md_tree(habib_maurer_algorithm(graphs.IcosahedralGraph()))
        PRIME
         3
         4
         7
         9
         11
         1
         5
         8
         0
         2
         6
         10
    """

    def recursive_print_md_tree(root, level):
        """
        Print the modular decomposition tree at root.

        INPUT:

        - ``root`` -- root of the modular decomposition tree

        - ``level`` -- indicates the depth of root in the original modular
          decomposition tree
        """
        if root.node_type != NodeType.NORMAL:
            print("{}{}".format(level, str(root.node_type)))
            for tree in root.children:
                recursive_print_md_tree(tree, level + " ")
        else:
            print("{}{}".format(level, str(root.children[0])))

    recursive_print_md_tree(root, "")


# =============================================================================
# Habib Maurer algorithm
# =============================================================================

def gamma_classes(graph):
    """
    Partition the edges of the graph into Gamma classes.

    Two distinct edges are Gamma related if they share a vertex but are not
    part of a triangle.  A Gamma class of edges is a collection of edges such
    that any edge in the class can be reached from any other by a chain of
    Gamma related edges (that are also in the class).

    The two important properties of the Gamma class

    * The vertex set corresponding to a Gamma class is a module
    * If the graph is not fragile (neither it or its complement is
      disconnected) then there is exactly one class that visits all the
      vertices of the graph, and this class consists of just the edges that
      connect the maximal strong modules of that graph.

    EXAMPLES:

    The gamma_classes of the octahedral graph are the three 4-cycles
    corresponding to the slices through the center of the octahedron::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import gamma_classes
        sage: g = graphs.OctahedralGraph()
        sage: sorted(gamma_classes(g), key=str)
        [frozenset({0, 1, 4, 5}), frozenset({0, 2, 3, 5}), frozenset({1, 2, 3, 4})]

    TESTS:

    Ensure that the returned vertex sets from some random graphs are modules::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import test_gamma_modules
        sage: test_gamma_modules(2, 10, 0.5)
    """
    from itertools import chain
    from sage.sets.disjoint_set import DisjointSet

    pieces = DisjointSet(frozenset(e) for e in graph.edge_iterator(labels=False))
    for v in graph:
        neighborhood = graph.subgraph(vertices=graph.neighbors(v))
        for component in neighborhood.complement().connected_components(sort=False):
            v1 = component[0]
            e = frozenset([v1, v])
            for vi in component[1:]:
                ei = frozenset([vi, v])
                pieces.union(e, ei)
    return {frozenset(chain.from_iterable(loe)): loe for loe in pieces}


def habib_maurer_algorithm(graph, g_classes=None):
    """
    Compute the modular decomposition by the algorithm of Habib and Maurer.

    Compute the modular decomposition of the given graph by the algorithm of
    Habib and Maurer [HM1979]_ . If the graph is disconnected or its complement
    is disconnected return a tree with a ``PARALLEL`` or ``SERIES`` node at the
    root and children being the modular decomposition of the subgraphs induced
    by the components. Otherwise, the root is ``PRIME`` and the modules are
    identified by having identical neighborhoods in the gamma class that spans
    the vertices of the subgraph (exactly one is guaranteed to exist). The gamma
    classes only need to be computed once, as the algorithm computes the the
    classes for the current root and each of the submodules. See also [BM1983]_
    for an equivalent algorithm described in greater detail.

    This function should not be used directly, it should be called via the
    ``modular_decomposition`` method of ``Graph`` with the parameter
    ``algorithm='habib_maurer'``.

    This functions assumes that ``graph`` is a object of the class ``Graph``, is
    a simple graph and has at least 1 vertex.

    INPUT:

    - ``graph`` -- the graph for which modular decomposition tree needs to be
      computed

    - ``g_classes`` -- dictionary (default: ``None``); a dictionary whose values
      are the gamma classes of the graph, and whose keys are a frozenset of the
      vertices corresponding to the class. Used internally.

    OUTPUT: the modular decomposition tree of the graph

    EXAMPLES:

    The Icosahedral graph is Prime::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: print_md_tree(habib_maurer_algorithm(graphs.IcosahedralGraph()))
        PRIME
         3
         4
         7
         9
         11
         1
         5
         8
         0
         2
         6
         10

    The Octahedral graph is not Prime::

        sage: print_md_tree(habib_maurer_algorithm(graphs.OctahedralGraph()))
        SERIES
         PARALLEL
          0
          5
         PARALLEL
          1
          4
         PARALLEL
          2
          3

    Tetrahedral Graph is Series::

        sage: print_md_tree(habib_maurer_algorithm(graphs.TetrahedralGraph()))
        SERIES
         0
         1
         2
         3

    Modular Decomposition tree containing both parallel and series modules::

        sage: d = {2:[4,3,5], 1:[4,3,5], 5:[3,2,1,4], 3:[1,2,5], 4:[1,2,5]}
        sage: g = Graph(d)
        sage: print_md_tree(habib_maurer_algorithm(g))
        SERIES
         PARALLEL
          1
          2
         PARALLEL
          3
          4
         5

    Graph from Marc Tedder implementation of modular decomposition::

        sage: d = {1:[5,4,3,24,6,7,8,9,2,10,11,12,13,14,16,17], 2:[1],
        ....:       3:[24,9,1], 4:[5,24,9,1], 5:[4,24,9,1], 6:[7,8,9,1],
        ....:       7:[6,8,9,1], 8:[6,7,9,1], 9:[6,7,8,5,4,3,1], 10:[1],
        ....:       11:[12,1], 12:[11,1], 13:[14,16,17,1], 14:[13,17,1],
        ....:       16:[13,17,1], 17:[13,14,16,18,1], 18:[17], 24:[5,4,3,1]}
        sage: g = Graph(d)
        sage: test_modular_decomposition(habib_maurer_algorithm(g), g)
        True

    Tetrahedral Graph is Series::

        sage: print_md_tree(habib_maurer_algorithm(graphs.TetrahedralGraph()))
        SERIES
         0
         1
         2
         3

    Modular Decomposition tree containing both parallel and series modules::

        sage: d = {2:[4,3,5], 1:[4,3,5], 5:[3,2,1,4], 3:[1,2,5], 4:[1,2,5]}
        sage: g = Graph(d)
        sage: print_md_tree(habib_maurer_algorithm(g))
        SERIES
         PARALLEL
          1
          2
         PARALLEL
          3
          4
         5

    TESTS:

    Ensure that a random graph and an isomorphic graph have identical modular
    decompositions. ::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import permute_decomposition
        sage: permute_decomposition(2, habib_maurer_algorithm, 20, 0.5)                 # needs sage.groups
    """
    if graph.order() == 1:
        root = Node.create_leaf(next(graph.vertex_iterator()))
        return root

    elif not graph.is_connected():
        root = Node(NodeType.PARALLEL)
        root.children = [habib_maurer_algorithm(graph.subgraph(vertices=sg), g_classes)
                         for sg in graph.connected_components(sort=False)]
        return root

    g_comp = graph.complement()
    if g_comp.is_connected():
        from collections import defaultdict
        root = Node(NodeType.PRIME)
        if g_classes is None:
            g_classes = gamma_classes(graph)
        vertex_set = frozenset(graph)
        edges = [tuple(e) for e in g_classes[vertex_set]]
        sub = graph.subgraph(edges=edges)
        d = defaultdict(list)
        for v in sub:
            for v1 in sub.neighbor_iterator(v):
                d[v1].append(v)
        d1 = defaultdict(list)
        for k, v in d.items():
            d1[frozenset(v)].append(k)
        root.children = [habib_maurer_algorithm(graph.subgraph(vertices=sg), g_classes)
                         for sg in d1.values()]
        return root

    root = Node(NodeType.SERIES)
    root.children = [habib_maurer_algorithm(graph.subgraph(vertices=sg), g_classes)
                     for sg in g_comp.connected_components(sort=False)]
    return root


################################################################################
#                   Exported modular_decomposition function                    #
################################################################################
def modular_decomposition(G, algorithm=None):
    r"""
    Return the modular decomposition of the current graph.

    This function should not be used directly, it should be called via the
    ``modular_decomposition`` method of ``Graph``.

    TESTS::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *

        sage: modular_decomposition(Graph())
        EMPTY []

        sage: modular_decomposition(Graph(1))
        NORMAL [0]

        sage: modular_decomposition(DiGraph())
        Traceback (most recent call last):
        ...
        TypeError: the input must be an undirected Sage graph

        sage: modular_decomposition(Graph(5, loops=True))
        Traceback (most recent call last):
        ...
        ValueError: This method is not known to work on graphs with loops...

        sage: modular_decomposition(Graph(5, multiedges=True))
        Traceback (most recent call last):
        ...
        ValueError: This method is not known to work on graphs with multiedges...

        sage: modular_decomposition(Graph(), algorithm='silly walk')
        Traceback (most recent call last):
        ...
        ValueError: unknown algorithm "silly walk"
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise TypeError("the input must be an undirected Sage graph")
    G._scream_if_not_simple()

    if algorithm is None:
        algorithm = "corneil_habib_paul_tedder"

    if algorithm not in ("habib_maurer", "corneil_habib_paul_tedder"):
        raise ValueError(f'unknown algorithm "{algorithm}"')

    if not G.order():
        return Node(NodeType.EMPTY)
    if G.order() == 1:
        D = Node(NodeType.NORMAL)
        D.children.append(next(G.vertex_iterator()))
        return D
    if algorithm == "habib_maurer":
        return habib_maurer_algorithm(G)
    return corneil_habib_paul_tedder_algorithm(G)


# ============================================================================
# Below functions are implemented to test the modular decomposition tree
# ============================================================================
# Function implemented for testing
def test_modular_decomposition(tree_root, graph):
    """
    Test the input modular decomposition tree using recursion.

    INPUT:

    - ``tree_root`` -- root of the modular decomposition tree to be tested

    - ``graph`` -- graph whose modular decomposition tree needs to be tested

    OUTPUT: ``True`` if input tree is a modular decomposition else ``False``

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = graphs.HexahedralGraph()
        sage: test_modular_decomposition(modular_decomposition(g), g)
        True
    """
    if tree_root.node_type != NodeType.NORMAL:
        for module in tree_root.children:
            if not test_module(module, graph):
                # test whether modules pass the defining
                # characteristics of modules
                return False
            if not test_modular_decomposition(module,
                                              graph.subgraph(get_vertices(module))):
                # recursively test the modular decomposition subtrees
                return False

        if not test_maximal_modules(tree_root, graph):
            # test whether the mdoules are maximal in nature
            return False

    return True


# Function implemented for testing
def test_maximal_modules(tree_root, graph):
    r"""
    Test the maximal nature of modules in a modular decomposition tree.

    Suppose the module `M = [M_1, M_2, \cdots, n]` is the input modular
    decomposition tree. Algorithm forms pairs like `(M_1, M_2), (M_1, M_3),
    \cdots, (M_1, M_n)`; `(M_2, M_3), (M_2, M_4), \cdots, (M_2, M_n)`; `\cdots`
    and so on and tries to form a module using the pair. If the module formed
    has same type as `M` and is of type ``SERIES`` or ``PARALLEL`` then the
    formed module is not considered maximal. Otherwise it is considered maximal
    and `M` is not a modular decomposition tree.

    INPUT:

    - ``tree_root`` -- modular decomposition tree whose modules are tested for
      maximal nature

    - ``graph`` -- graph whose modular decomposition tree is tested

    OUTPUT:

    ``True`` if all modules at first level in the modular decomposition tree
    are maximal in nature

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = graphs.HexahedralGraph()
        sage: test_maximal_modules(modular_decomposition(g), g)
        True
    """
    if tree_root.node_type != NodeType.NORMAL:
        for index, module in enumerate(tree_root.children):
            for other_index in range(index + 1, len(tree_root.children)):

                # compute the module formed using modules at index and
                # other_index
                module_formed = form_module(index, other_index,
                                            tree_root, graph)

                if module_formed[0]:
                    # Module formed and the parent of the formed module
                    # should not both be of type SERIES or PARALLEL
                    mod_type = get_module_type(graph.subgraph(module_formed[1]))
                    if (mod_type == tree_root.node_type and
                            (tree_root.node_type == NodeType.PARALLEL or
                             tree_root.node_type == NodeType.SERIES)):
                        continue
                    return False
    return True


def get_vertices(component_root):
    """
    Compute the list of vertices in the (co)component.

    INPUT:

    - ``component_root`` -- root of the (co)component whose vertices need to be
      returned as a list

    OUTPUT:

    list of vertices in the (co)component

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: forest = Node(NodeType.PRIME)
        sage: forest.children = [Node.create_leaf(2), Node.create_leaf(0),
        ....:                    Node.create_leaf(3), Node.create_leaf(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [Node.create_leaf(4),
        ....:                         Node.create_leaf(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [Node.create_leaf(6),
        ....:                           Node.create_leaf(7)]
        sage: forest.children.insert(1, series_node)
        sage: forest.children.insert(3, parallel_node)
        sage: get_vertices(forest)
        [2, 4, 5, 0, 6, 7, 3, 1]
    """
    vertices = []

    # inner recursive function to recurse over the elements in the
    # ``component``
    def recurse_component(node, vertices):
        if node.node_type == NodeType.NORMAL:
            vertices.append(node.children[0])
            return
        for child in node.children:
            recurse_component(child, vertices)

    recurse_component(component_root, vertices)
    return vertices


# Function implemented for testing
def get_module_type(graph):
    """
    Return the module type of the root of the modular decomposition tree of
    ``graph``.

    INPUT:

    - ``graph`` -- input sage graph

    OUTPUT:

    ``PRIME`` if graph is PRIME, ``PARALLEL`` if graph is PARALLEL and
    ``SERIES`` if graph is of type SERIES

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import get_module_type
        sage: g = graphs.HexahedralGraph()
        sage: get_module_type(g)
        PRIME
    """
    if not graph.is_connected():
        return NodeType.PARALLEL
    elif graph.complement().is_connected():
        return NodeType.PRIME
    return NodeType.SERIES


# Function implemented for testing
def form_module(index, other_index, tree_root, graph):
    r"""
    Forms a module out of the modules in the module pair.

    Let `M_1` and `M_2` be the input modules. Let `V` be the set of vertices in
    these modules. Suppose `x` is a neighbor of subset of the vertices in `V`
    but not all the vertices and `x` does not belong to `V`. Then the set of
    modules also include the module which contains `x`. This process is repeated
    until a module is formed and the formed module if subset of `V` is returned.

    INPUT:

    - ``index`` -- first module in the module pair

    - ``other_index`` -- second module in the module pair

    - ``tree_root`` -- modular decomposition tree which contains the modules
      in the module pair

    - ``graph`` -- graph whose modular decomposition tree is created

    OUTPUT:

    ``[module_formed, vertices]`` where ``module_formed`` is ``True`` if
    module is formed else ``False`` and ``vertices`` is a list of vertices
    included in the formed module

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = graphs.HexahedralGraph()
        sage: tree_root = modular_decomposition(g)
        sage: form_module(0, 2, tree_root, g)
        [False, {0, 1, 2, 3, 4, 5, 6, 7}]
    """
    vertices = set(get_vertices(tree_root.children[index]))
    vertices.update(get_vertices(tree_root.children[other_index]))

    # stores all neighbors which are common for all vertices in V
    common_neighbors = set()

    # stores all neighbors of vertices in V which are outside V
    all_neighbors = set()

    while True:
        # remove vertices from all_neighbors and common_neighbors
        all_neighbors.difference_update(vertices)
        common_neighbors.difference_update(vertices)

        for v in vertices:
            # stores the neighbors of v which are outside the set of vertices
            neighbor_list = set(graph.neighbors(v))
            neighbor_list.difference_update(vertices)

            # update all_neighbors and common_neighbors using the
            # neighbor_list
            all_neighbors.update(neighbor_list)
            common_neighbors.intersection_update(neighbor_list)

        if all_neighbors == common_neighbors:  # indicates a module is formed

            # module formed covers the entire graph
            if len(vertices) == graph.order():
                return [False, vertices]

            return [True, vertices]

        # add modules containing uncommon neighbors into the formed module
        for v in (all_neighbors - common_neighbors):
            for index in range(len(tree_root.children)):
                if v in get_vertices(tree_root.children[index]):
                    vertices.update(get_vertices(tree_root.children[index]))
                    break


# Function implemented for testing
def test_module(module, graph):
    """
    Test whether input module is actually a module.

    INPUT:

    - ``module`` -- module which needs to be tested

    - ``graph`` -- input sage graph which contains the module

    OUTPUT: ``True`` if input module is a module by definition else ``False``

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = graphs.HexahedralGraph()
        sage: tree_root = modular_decomposition(g)
        sage: test_module(tree_root, g)
        True
        sage: test_module(tree_root.children[0], g)
        True
    """
    # A single vertex is a module
    if module.node_type == NodeType.NORMAL:
        return True

    # vertices contained in module
    vertices_in_module = get_vertices(module)

    # vertices outside module
    vertices_outside = list(set(graph.vertices(sort=False)) - set(vertices_in_module))

    # Nested module with only one child
    if module.node_type != NodeType.NORMAL and len(module.children) == 1:
        return False

    # If children of SERIES module are all SERIES modules
    if module.node_type == NodeType.SERIES:
        if children_node_type(module, NodeType.SERIES):
            return False

    # If children of PARALLEL module are all PARALLEL modules
    if module.node_type == NodeType.PARALLEL:
        if children_node_type(module, NodeType.PARALLEL):
            return False

    # check the module by definition. Vertices in a module should all either
    # be connected or disconnected to any vertex outside module
    for v in vertices_outside:
        if not either_connected_or_not_connected(v, vertices_in_module, graph):
            return False
    return True


# Function implemented for testing
def children_node_type(module, node_type):
    """
    Check whether the node type of the children of ``module`` is ``node_type``.

    INPUT:

    - ``module`` -- module which is tested

    - ``node_type`` -- input node_type

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = graphs.OctahedralGraph()
        sage: tree_root = modular_decomposition(g)
        sage: print_md_tree(tree_root)
        SERIES
         PARALLEL
          2
          3
         PARALLEL
          1
          4
         PARALLEL
          0
          5
        sage: children_node_type(tree_root, NodeType.SERIES)
        False
        sage: children_node_type(tree_root, NodeType.PARALLEL)
        True
    """
    return all(node.node_type == node_type for node in module.children)


# Function implemented for testing
def either_connected_or_not_connected(v, vertices_in_module, graph):
    """
    Check whether ``v`` is connected or disconnected to all vertices in the
    module.

    INPUT:

    - ``v`` -- vertex tested

    - ``vertices_in_module`` -- list containing vertices in the module

    - ``graph`` -- graph to which the vertices belong

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = graphs.OctahedralGraph()
        sage: print_md_tree(modular_decomposition(g))
        SERIES
         PARALLEL
          2
          3
         PARALLEL
          1
          4
         PARALLEL
          0
          5
        sage: either_connected_or_not_connected(2, [1, 4], g)
        True
        sage: either_connected_or_not_connected(2, [3, 4], g)
        False
    """
    # marks whether vertex v is connected to first vertex in the module
    connected = graph.has_edge(vertices_in_module[0], v)

    # if connected is True then all vertices in module should be connected to
    # v else all should be disconnected
    return all(graph.has_edge(u, v) == connected for u in vertices_in_module)


def tree_to_nested_tuple(root):
    r"""
    Convert a modular decomposition tree to a nested tuple.

    INPUT:

    - ``root`` -- the root of the modular decomposition tree

    OUTPUT:

    A tuple whose first element is the type of the root of the tree and whose
    subsequent nodes are either vertex labels in the case of leaves or tuples
    representing the child subtrees.

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = graphs.OctahedralGraph()
        sage: tree_to_nested_tuple(modular_decomposition(g))
        (SERIES, [(PARALLEL, [2, 3]), (PARALLEL, [1, 4]), (PARALLEL, [0, 5])])
    """
    if root.node_type == NodeType.NORMAL:
        return root.children[0]
    else:
        return (root.node_type, [tree_to_nested_tuple(x) for x in root.children])


def nested_tuple_to_tree(nest):
    r"""
    Turn a tuple representing the modular decomposition into a tree.

    INPUT:

    - ``nest`` -- a nested tuple of the form returned by
      :meth:`tree_to_nested_tuple`

    OUTPUT: the root node of a modular decomposition tree

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: tree = (NodeType.SERIES, 1, 2, (NodeType.PARALLEL, 3, 4))
        sage: print_md_tree(nested_tuple_to_tree(tree))
        SERIES
         1
         2
         PARALLEL
          3
          4
    """
    if not isinstance(nest, tuple):
        return Node.create_leaf(nest)

    root = Node(nest[0])
    root.children = [nested_tuple_to_tree(n) for n in nest[1:]]
    return root


def equivalent_trees(root1, root2):
    r"""
    Check that two modular decomposition trees are the same.

    Verify that the structure of the trees is the same. Two leaves are
    equivalent if they represent the same vertex, two internal nodes are
    equivalent if they have the same nodes type and the same number of children
    and there is a matching between the children such that each pair of
    children is a pair of equivalent subtrees.

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: t1 = nested_tuple_to_tree((NodeType.SERIES, 1, 2,
        ....:             (NodeType.PARALLEL, 3, 4)))
        sage: t2 = nested_tuple_to_tree((NodeType.SERIES,
        ....:             (NodeType.PARALLEL, 4, 3), 2, 1))
        sage: equivalent_trees(t1, t2)
        True
    """
    # internal definition
    def node_id(root):
        return (root.node_type, frozenset(get_vertices(root)))

    if root1.node_type != root2.node_type:
        return False

    if len(root1.children) != len(root2.children):
        return False

    if root1.node_type == NodeType.NORMAL:
        return root1.children[0] == root2.children[0]

    child_map = {}
    for node in root2.children:
        child_map[node_id(node)] = node

    for node in root1.children:
        id = node_id(node)
        if id not in child_map:
            return False
        if not equivalent_trees(node, child_map[id]):
            return False

    return True


def relabel_tree(root, perm):
    r"""
    Relabel the leaves of a tree according to a dictionary.

    INPUT:

    - ``root`` -- the root of the tree

    - ``perm`` -- a function, dictionary, list, permutation, or ``None``
      representing the relabeling. See
      :meth:`~sage.graphs.generic_graph.GenericGraph.relabel` for description of
      the permutation input.

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: tuple_tree = (NodeType.SERIES, 1, 2, (NodeType.PARALLEL, 3, 4))
        sage: tree = nested_tuple_to_tree(tuple_tree)
        sage: print_md_tree(relabel_tree(tree, (4,3,2,1)))
        SERIES
         4
         3
         PARALLEL
          2
          1
    """
    # If perm is not a dictionary, we build one !
    if perm is None:

        # vertices() returns a sorted list:
        # this guarantees consistent relabeling
        perm = {v: i for i, v in enumerate(get_vertices(root))}

    elif isinstance(perm, dict):
        from copy import copy
        # If all vertices do not have a new label, the code will touch the
        # dictionary. Let us keep the one we received from the user clean !
        perm = copy(perm)

    elif isinstance(perm, (list, tuple)):
        perm = dict(zip(sorted(get_vertices(root)), perm))

    elif isinstance(perm, PermutationGroupElement):
        n = len(get_vertices(root))
        ddict = {}
        for i in range(1, n):
            ddict[i] = perm(i) % n
        if n > 0:
            ddict[0] = perm(n) % n
        perm = ddict

    elif callable(perm):
        perm = {i: perm(i) for i in get_vertices(root)}

    else:
        raise TypeError("type of perm is not supported for relabeling")

    if root.node_type == NodeType.NORMAL:
        return Node.create_leaf(perm[root.children[0]])
    else:
        new_root = Node(root.node_type)
        new_root.children = [relabel_tree(child, perm) for child in root.children]
        return new_root


# =============================================================================
# Random tests
# =============================================================================

@random_testing
def test_gamma_modules(trials, vertices, prob, verbose=False):
    r"""
    Verify that the vertices of each gamma class of a random graph are modules
    of that graph.

    INPUT:

    - ``trials`` -- the number of trials to run

    - ``vertices`` -- the size of the graph to use

    - ``prob`` -- the probability that any given edge is in the graph
      See :meth:`~sage.graphs.generators.random.RandomGNP` for more details

    - ``verbose`` -- print information on each trial

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: test_gamma_modules(3, 7, 0.5)
    """
    from sage.graphs.generators.random import RandomGNP
    for _ in range(trials):
        g = RandomGNP(vertices, prob)
        if verbose:
            print(g.graph6_string())
        g_classes = gamma_classes(g)
        for module in g_classes.keys():
            m_list = list(module)
            for v in g:
                if v not in module:
                    assert either_connected_or_not_connected(v, m_list, g)
        if verbose:
            print("Passes!")


@random_testing
def permute_decomposition(trials, algorithm, vertices, prob, verbose=False):
    r"""
    Check that a graph and its permuted relabeling have the same modular
    decomposition.

    We generate a ``trials`` random graphs and then generate an isomorphic graph
    by relabeling the original graph. We then verify

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: permute_decomposition(30, habib_maurer_algorithm, 10, 0.5)
    """
    from sage.graphs.generators.random import RandomGNP
    from sage.combinat.permutation import Permutations
    for _ in range(trials):
        g1 = RandomGNP(vertices, prob)
        random_perm = Permutations(list(g1)).random_element()
        g2 = g1.relabel(perm=random_perm, inplace=False)
        if verbose:
            print(g1.graph6_string())
            print(random_perm)
        t1 = algorithm(g1)
        t2 = algorithm(g2)
        assert test_modular_decomposition(t1, g1)
        assert test_modular_decomposition(t2, g2)
        t1p = relabel_tree(t1, random_perm)
        assert equivalent_trees(t1p, t2)
        if verbose:
            print("Passes!")


def random_md_tree(max_depth, max_fan_out, leaf_probability):
    r"""
    Create a random MD tree.

    INPUT:

    - ``max_depth`` -- the maximum depth of the tree

    - ``max_fan_out`` -- the maximum number of children a node can have
      (must be >=4 as a prime node must have at least 4 vertices)

    - ``leaf_probability`` -- the probability that a subtree is a leaf

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: set_random_seed(0)
        sage: tree_to_nested_tuple(random_md_tree(2, 5, 0.5))
        (PRIME, [0, 1, (PRIME, [2, 3, 4, 5, 6]), 7, (PARALLEL, [8, 9, 10])])
    """

    from sage.misc.prandom import choice, randint, random

    if max_fan_out < 4:
        raise ValueError("max_fan_out must be at least 4")

    # Internal function
    def rand_md_tree(max_depth, parent_type):
        r"""
        Create the subtrees of a node.

        A child of a node cannot have the same node type as its parent if its
        parent's node type is either PARALLEL or SERIES.  Also its ``max_depth``
        is one less than its parent's.
        """
        if random() < leaf_probability or max_depth == 1:
            root = Node.create_leaf(current_leaf[0])
            current_leaf[0] += 1
            return root
        if parent_type == NodeType.PRIME:
            node_type = choice([NodeType.PRIME, NodeType.SERIES, NodeType.PARALLEL])
        elif parent_type == NodeType.SERIES:
            node_type = choice([NodeType.PRIME, NodeType.PARALLEL])
        else:
            node_type = choice([NodeType.PRIME, NodeType.SERIES])
        if node_type == NodeType.PRIME:
            num_children = randint(4, max_fan_out)
        else:
            num_children = randint(2, max_fan_out)
        root = Node(node_type)
        root.children = [rand_md_tree(max_depth - 1, node_type)
                         for _ in range(num_children)]
        return root

    # a hack around python2's lack of 'nonlocal'
    current_leaf = [0]
    node_type = choice([NodeType.PRIME, NodeType.SERIES, NodeType.PARALLEL])
    num_children = randint(4, max_fan_out)
    root = Node(node_type)
    root.children = [rand_md_tree(max_depth, node_type)
                     for _ in range(num_children)]
    return root


def md_tree_to_graph(root, prime_node_generator=None):
    r"""
    Create a graph having the given MD tree.

    For the prime nodes, the parameter ``prime_node_generator`` is called with
    the number of vertices as the only argument. If it is ``None``, the path
    graph is used (it is prime when the length is 4 or more).

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: from sage.graphs.graph_generators import graphs
        sage: tup1 = (NodeType.PRIME, 1, (NodeType.SERIES, 2, 3),
        ....:        (NodeType.PARALLEL, 4, 5), 6)
        sage: tree1 = nested_tuple_to_tree(tup1)
        sage: g1 = md_tree_to_graph(tree1)
        sage: g2 = Graph({1: [2, 3], 2: [1, 3, 4, 5], 3: [1, 2, 4, 5],
        ....:             4: [2, 3, 6], 5: [2, 3, 6], 6: [4, 5]})
        sage: g1.is_isomorphic(g2)
        True

        sage: G = md_tree_to_graph(Node(NodeType.EMPTY))
        sage: G.is_isomorphic(Graph())
        True

        sage: tree = Node(NodeType.SERIES)
        sage: tree.children.extend(Node.create_leaf(i) for i in range(5))
        sage: G = md_tree_to_graph(tree)
        sage: G.is_isomorphic(graphs.CompleteGraph(5))
        True

        sage: tree = Node(NodeType.PRIME)
        sage: tree.children.extend(Node.create_leaf(i) for i in range(5))
        sage: png = lambda n: (graphs.PathGraph if n == 4 else graphs.CycleGraph)(n)
        sage: G = md_tree_to_graph(tree, prime_node_generator=png)
        sage: G.is_isomorphic(graphs.CycleGraph(5))
        True
    """
    from itertools import product, combinations
    from sage.graphs.graph import Graph

    if prime_node_generator is None:
        from sage.graphs.graph_generators import graphs
        prime_node_generator = graphs.PathGraph

    def tree_to_vertices_and_edges(root):
        r"""
        Give the list of vertices and edges of the graph having the given md tree.
        """
        if root.is_leaf():
            return (root.children, [])
        children_ve = [tree_to_vertices_and_edges(child) for child in root.children]
        vertices = [v for vs, es in children_ve for v in vs]
        edges = [e for vs, es in children_ve for e in es]
        vertex_lists = [vs for vs, es in children_ve]
        if root.is_prime():
            G = prime_node_generator(len(vertex_lists))
            G.relabel(range(len(vertex_lists)))
            for i1, i2 in G.edge_iterator(labels=False):
                edges.extend(product(vertex_lists[i1], vertex_lists[i2]))
        elif root.is_series():
            for vs1, vs2 in combinations(vertex_lists, 2):
                edges.extend(product(vs1, vs2))
        # else: no edge to be created for PARALLEL nodes
        return (vertices, edges)

    if root.is_empty():
        return Graph()
    vs, es = tree_to_vertices_and_edges(root)
    return Graph([vs, es], format='vertices_and_edges')


@random_testing
def recreate_decomposition(trials, algorithm, max_depth, max_fan_out,
                           leaf_probability, verbose=False):
    r"""
    Verify that we can recreate a random MD tree.

    We create a random MD tree, then create a graph having that decomposition,
    then find a modular decomposition for that graph, and verify that the two
    modular decomposition trees are equivalent.

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: recreate_decomposition(3, habib_maurer_algorithm, 4, 6, 0.5,
        ....:                         verbose=False)
    """
    for _ in range(trials):
        rand_tree = random_md_tree(max_depth, max_fan_out, leaf_probability)
        if verbose:
            print_md_tree(rand_tree)
        graph = md_tree_to_graph(rand_tree)
        if verbose:
            print(graph.graph6_string())
            print(graph.to_dictionary())
        reconstruction = algorithm(graph)
        if verbose:
            print_md_tree(reconstruction)
        assert equivalent_trees(rand_tree, reconstruction)
        if verbose:
            print("Passes!")
