# -*- coding: utf-8 -*-
r"""
Modular Decomposition

This module implements the function for computing the modular decomposition
of undirected graphs.
"""
# ****************************************************************************
#       Copyright (C) 2017 Lokesh Jain <lokeshj1703@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from collections import deque
from enum import Enum


class NodeType(Enum):
    """
    NodeType is an enumeration class used to define the various types of nodes
    in modular decomposition tree.

    The various node types defined are

    - ``PARALLEL`` -- indicates the node is a parallel module

    - ``SERIES`` -- indicates the node is a series module

    - ``PRIME`` -- indicates the node is a prime module

    - ``FOREST`` -- indicates a forest containing trees

    - ``NORMAL`` -- indicates the node is normal containing a vertex
    """
    PRIME = 0
    SERIES = 1
    PARALLEL = 2
    NORMAL = 3
    FOREST = -1

    def __repr__(self):
        r"""
        String representation of this node type.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.modular_decomposition import NodeType
            sage: repr(NodeType.PARALLEL)
            'PARALLEL'
        """
        return self.name

    def __str__(self):
        """
        String representation of this node type.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.modular_decomposition import NodeType
            sage: str(NodeType.PARALLEL)
            'PARALLEL'
        """
        return repr(self)


class NodeSplit(Enum):
    """
    Enumeration class used to specify the split that has occurred at the node or
    at any of its descendants.

    ``NodeSplit`` is defined for every node in modular decomposition tree and is
    required during the refinement and promotion phase of modular decomposition
    tree computation. Various node splits defined are

    - ``LEFT_SPLIT`` -- indicates a left split has occurred

    - ``RIGHT_SPLIT`` -- indicates a right split has occurred

    - ``BOTH_SPLIT`` -- indicates both left and right split have occurred

    - ``NO_SPLIT`` -- indicates no split has occurred
    """
    LEFT_SPLIT = 1
    RIGHT_SPLIT = 2
    BOTH_SPLIT = 3
    NO_SPLIT = 0


class VertexPosition(Enum):
    """
    Enumeration class used to define position of a vertex w.r.t source in
    modular decomposition.

    For computing modular decomposition of connected graphs a source vertex is
    chosen. The position of vertex is w.r.t this source vertex. The various
    positions defined are

    - ``LEFT_OF_SOURCE`` -- indicates vertex is to left of source and is a
      neighbour of source vertex

    - ``RIGHT_OF_SOURCE`` -- indicates vertex is to right of source and is
      connected to but not a neighbour of source vertex

    - ``SOURCE`` -- indicates vertex is source vertex
    """
    LEFT_OF_SOURCE = -1
    RIGHT_OF_SOURCE = 1
    SOURCE = 0


class Node:
    """
    Node class stores information about the node type, node split and index of
    the node in the parent tree.

    Node type can be ``PRIME``, ``SERIES``, ``PARALLEL``, ``NORMAL`` or
    ``FOREST``. Node split can be ``NO_SPLIT``, ``LEFT_SPLIT``, ``RIGHT_SPLIT``
    or ``BOTH_SPLIT``. A node is split in the refinement phase and the split
    used is propagated to the ancestors.

    - ``node_type`` -- is of type NodeType and specifies the type of node

    - ``node_split`` -- is of type NodeSplit and specifies the type of splits
      which have occurred in the node and its descendants

    - ``index_in_root`` -- specifies the index of the node in the forest
      obtained after promotion phase

    - ``comp_num`` -- specifies the number given to nodes in a (co)component
      before refinement

    - ``is_separated`` -- specifies whether a split has occurred with the node
      as the root
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
        self.node_split = NodeSplit.NO_SPLIT
        self.index_in_root = -1
        self.comp_num = -1
        self.is_separated = False
        self.children = []

    def set_node_split(self, node_split):
        """
        Add node_split to the node split of self.

        ``LEFT_SPLIT`` and ``RIGHT_SPLIT`` can exist together in ``self`` as
        ``BOTH_SPLIT``.

        INPUT:

        - ``node_split`` -- node_split to be added to self

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.modular_decomposition import *
            sage: node = Node(NodeType.PRIME)
            sage: node.set_node_split(NodeSplit.LEFT_SPLIT)
            sage: node.node_split == NodeSplit.LEFT_SPLIT
            True
            sage: node.set_node_split(NodeSplit.RIGHT_SPLIT)
            sage: node.node_split == NodeSplit.BOTH_SPLIT
            True
            sage: node = Node(NodeType.PRIME)
            sage: node.set_node_split(NodeSplit.BOTH_SPLIT)
            sage: node.node_split == NodeSplit.BOTH_SPLIT
            True
        """
        if self.node_split == NodeSplit.NO_SPLIT:
            self.node_split = node_split
        elif ((self.node_split == NodeSplit.LEFT_SPLIT and
                       node_split == NodeSplit.RIGHT_SPLIT) or
                  (self.node_split == NodeSplit.RIGHT_SPLIT and
                           node_split == NodeSplit.LEFT_SPLIT)):
            self.node_split = NodeSplit.BOTH_SPLIT

    def has_left_split(self):
        """
        Check whether ``self`` has ``LEFT_SPLIT``.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.modular_decomposition import *
            sage: node = Node(NodeType.PRIME)
            sage: node.set_node_split(NodeSplit.LEFT_SPLIT)
            sage: node.has_left_split()
            True
            sage: node = Node(NodeType.PRIME)
            sage: node.set_node_split(NodeSplit.BOTH_SPLIT)
            sage: node.has_left_split()
            True
        """
        return (self.node_split == NodeSplit.LEFT_SPLIT or
                self.node_split == NodeSplit.BOTH_SPLIT)

    def has_right_split(self):
        """
        Check whether ``self`` has ``RIGHT_SPLIT``.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.modular_decomposition import *
            sage: node = Node(NodeType.PRIME)
            sage: node.set_node_split(NodeSplit.RIGHT_SPLIT)
            sage: node.has_right_split()
            True
            sage: node = Node(NodeType.PRIME)
            sage: node.set_node_split(NodeSplit.BOTH_SPLIT)
            sage: node.has_right_split()
            True
        """
        return (self.node_split == NodeSplit.RIGHT_SPLIT or
                self.node_split == NodeSplit.BOTH_SPLIT)

    def __repr__(self):
        r"""
        Return a string representation of the node.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.modular_decomposition import *
            sage: n = Node(NodeType.PRIME)
            sage: n.children.append(create_normal_node(1))
            sage: n.children.append(create_normal_node(2))
            sage: str(n)
            'PRIME [NORMAL [1], NORMAL [2]]'
        """
        if self.node_type == NodeType.SERIES:
            s = "SERIES "
        elif self.node_type == NodeType.PARALLEL:
            s = "PARALLEL "
        elif self.node_type == NodeType.PRIME:
            s = "PRIME "
        elif self.node_type == NodeType.FOREST:
            s = "FOREST "
        else:
            s = "NORMAL "

        s += str(self.children)
        return s

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
                self.node_split == other.node_split and
                self.index_in_root == other.index_in_root and
                self.comp_num == other.comp_num and
                self.is_separated == other.is_separated and
                self.children == other.children)


def modular_decomposition(graph):
    """
    Compute the modular decomposition tree of ``graph``.

    The tree structure is represented in form of nested lists. A tree node is
    an object of type Node. The Node object further contains a list of its
    children

    INPUT:

    - ``graph`` -- the graph for which modular decomposition tree needs to be
      computed

    OUTPUT:

    A nested list representing the modular decomposition tree computed for the
    graph

    EXAMPLES:

    The Icosahedral graph is Prime::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: print_md_tree(modular_decomposition(graphs.IcosahedralGraph()))
        PRIME
         5
         7
         11
         1
         8
         0
         9
         4
         10
         6
         2
         3

    The Octahedral graph is not Prime::

        sage: print_md_tree(modular_decomposition(graphs.OctahedralGraph()))
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

    Tetrahedral Graph is Series::

        sage: print_md_tree(modular_decomposition(graphs.TetrahedralGraph()))
        SERIES
              3
              2
              1
              0

    Modular Decomposition tree containing both parallel and series modules::

        sage: d = {2:[4,3,5], 1:[4,3,5], 5:[3,2,1,4], 3:[1,2,5], 4:[1,2,5]}
        sage: g = Graph(d)
        sage: print_md_tree(modular_decomposition(g))
        SERIES
              5
              PARALLEL
                3
                4
              PARALLEL
                1
                2

    TESTS:

    Bad Input::

        sage: g = DiGraph()
        sage: modular_decomposition(g)
        Traceback (most recent call last):
        ...
        ValueError: Graph must be undirected

    Empty Graph is Prime::

        sage: g = Graph()
        sage: modular_decomposition(g)
        PRIME []

    Graph from Marc Tedder implementation of modular decomposition::

        sage: d = {1:[5,4,3,24,6,7,8,9,2,10,11,12,13,14,16,17], 2:[1],
        ....:       3:[24,9,1], 4:[5,24,9,1], 5:[4,24,9,1], 6:[7,8,9,1],
        ....:       7:[6,8,9,1], 8:[6,7,9,1], 9:[6,7,8,5,4,3,1], 10:[1],
        ....:       11:[12,1], 12:[11,1], 13:[14,16,17,1], 14:[13,17,1],
        ....:       16:[13,17,1], 17:[13,14,16,18,1], 18:[17], 24:[5,4,3,1]}
        sage: g = Graph(d)
        sage: test_modular_decomposition(modular_decomposition(g), g)
        True

    Graph from the :wikipedia:`Modular_decomposition`::

        sage: d2 = {1:[2,3,4], 2:[1,4,5,6,7], 3:[1,4,5,6,7], 4:[1,2,3,5,6,7],
        ....:       5:[2,3,4,6,7], 6:[2,3,4,5,8,9,10,11],
        ....:       7:[2,3,4,5,8,9,10,11], 8:[6,7,9,10,11], 9:[6,7,8,10,11],
        ....:       10:[6,7,8,9], 11:[6,7,8,9]}
        sage: g = Graph(d2)
        sage: test_modular_decomposition(modular_decomposition(g), g)
        True
    """
    if graph.is_directed():
        raise ValueError("Graph must be undirected")

    if not graph.order():  # Empty Graph
        return create_prime_node()

    if graph.order() == 1:  # Single vertex graph
        root = create_normal_node(next(graph.vertex_iterator()))
        return root

    if not graph.is_connected():

        # Parallel case: the tree contains the MD trees of its connected
        # components as subtrees
        components = graph.connected_components()
        root = create_parallel_node()
        for component in components:
            root.children.append(modular_decomposition(graph.subgraph(component)))
        return root
    elif graph.complement().is_connected(): # Prime Graph
        root = create_prime_node()
    else:
        root = create_series_node()  # Series Graph

    bfs_generator = graph.breadth_first_search(next(graph.vertex_iterator()),
                                               report_distance=True)

    prev_level_distance = -1  # used as a demarker for different levels in bfs
    prev_level_list = []  # stores the vertices in previous level

    # dictionary stores the distance of vertices from the SOURCE
    vertex_dist = {}

    # dictionary stores the position of vertices w.r.t SOURCE
    vertex_status = {}
    vertex_status[next(graph.vertex_iterator())] = VertexPosition.SOURCE

    # Different distances from the source vertex are considered
    # as different levels in the algorithm
    for (vertex, distance) in bfs_generator:
        vertex_dist[vertex] = distance

        # Mark the neighbours of source as LEFT_OF_SOURCE as they appear to
        # left of source in the forest, other vertices are marked as
        # RIGHT_OF_SOURCE
        if distance == 1:
            vertex_status[vertex] = VertexPosition.LEFT_OF_SOURCE
        elif distance != 0:
            vertex_status[vertex] = VertexPosition.RIGHT_OF_SOURCE

        if distance != prev_level_distance:  # On start of new level in BFS
            if prev_level_list:
                # MD Tree is computed for each level and added to the forest
                root.children.append(modular_decomposition(graph.subgraph(prev_level_list)))
            prev_level_list = []
            prev_level_distance = distance
        prev_level_list.append(vertex)

    # The last level is left out in the above loop
    root.children.append(modular_decomposition(graph.subgraph(prev_level_list)))

    # The MD tree for the neighbours of source marked as LEFT_OF_SOURCE
    # are placed left of Source in the forest. root.children[1] is required to
    # be source and root.children[0] is required to be the MD tree for the
    # neighbours therefore, the first two elements in the list are replaced
    root.children[0], root.children[1] = root.children[1], root.children[0]

    root.node_type = NodeType.FOREST
    clear_node_split_info(root)
    number_cocomponents(root, vertex_status)
    number_components(root, vertex_status)
    refine(graph, root, vertex_dist, vertex_status)
    promote_left(root)
    promote_right(root)
    promote_child(root)
    assembly(graph, root, vertex_status, vertex_dist)

    if root.node_type == NodeType.FOREST:
        return root.children[0]
    else:
        return root


def number_components(root, vertex_status):
    """
    Number the components to the right of SOURCE vertex in the forest input to
    the assembly phase

    INPUT:

    - ``root`` -- the forest which contains the components and cocomponents

    - ``vertex_status`` -- dictionary which stores the position of vertex
      w.r.t SOURCE

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                           create_normal_node(7)]
        sage: forest.children.append(series_node)
        sage: forest.children.append(parallel_node)
        sage: vertex_status = {2: VertexPosition.LEFT_OF_SOURCE,
        ....:                  3: VertexPosition.SOURCE,
        ....:                  1: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  4: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  5: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  6: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  7: VertexPosition.RIGHT_OF_SOURCE}
        sage: number_components(forest, vertex_status)
        sage: forest.children[-1].children[0].comp_num
        2
        sage: forest.children[-1].children[1].comp_num
        3

    TESTS::

        sage: (forest.children[-1].children[0].comp_num == 2 and
        ....:  forest.children[-1].children[1].comp_num == 3)
        True
        sage: (forest.children[-2].children[0].comp_num == 1 and
        ....:  forest.children[-2].children[1].comp_num == 1)
        True
    """
    comp_num = 0
    flag = False

    if not root:  # root is empty
        return ValueError("Input forest {} is empty".format(root))

    for node in root.children:

        # flag set to True after source vertex is encountered
        if (node.node_type == NodeType.NORMAL and
            vertex_status[node.children[0]] == VertexPosition.SOURCE):
            flag = True
            continue

        if not flag:  # Cocomponents are skipped
            continue

        comp_num += recursively_number_parts(node, comp_num, NodeType.PARALLEL)


def number_cocomponents(root, vertex_status):
    """
    Number the cocomponents to the left of SOURCE vertex in the forest input to
    the assembly phase

    INPUT:

    - ``root`` -- the forest which contains the cocomponents and components

    - ``vertex_status`` -- dictionary which stores the position of vertex
      w.r.t SOURCE

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                           create_normal_node(7)]
        sage: forest.children.insert(1, series_node)
        sage: forest.children.insert(2, parallel_node)
        sage: vertex_status = {2: VertexPosition.LEFT_OF_SOURCE,
        ....:                  3: VertexPosition.SOURCE,
        ....:                  1: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  4: VertexPosition.LEFT_OF_SOURCE,
        ....:                  5: VertexPosition.LEFT_OF_SOURCE,
        ....:                  6: VertexPosition.LEFT_OF_SOURCE,
        ....:                  7: VertexPosition.LEFT_OF_SOURCE}
        sage: number_cocomponents(forest, vertex_status)
        sage: forest.children[1].children[0].comp_num
        1
        sage: forest.children[1].children[1].comp_num
        2

    TESTS::

        sage: (forest.children[1].children[0].comp_num and
        ....:  forest.children[1].children[1].comp_num == 2)
        True
        sage: (forest.children[2].children[0].comp_num == 3 and
        ....:  forest.children[2].children[1].comp_num == 3)
        True
    """
    cocomp_num = 0
    for node in root.children:
        # Only cocomponents are numbered
        if (node.node_type == NodeType.NORMAL and
            vertex_status[node.children[0]] == VertexPosition.SOURCE):
            break
        cocomp_num += recursively_number_parts(node, cocomp_num, NodeType.SERIES)


def recursively_number_parts(part_root, part_num, by_type):
    """
    Recursively number the nodes in the (co)components(parts).

    If the ``node_type`` of ``part_root`` is same as ``by_type`` then
    ``part_num`` is incremented for subtree at each child of ``part_root`` else
    part is numbered by ``part_num``.

    INPUT:

    - ``part_root`` -- root of the part to be numbered

    - ``part_num`` -- input number to be used as reference for numbering
      the (co)components

    - ``by_type`` -- type which determines how numbering is done

    OUTPUT:

    The value incremented to ``part_num``.

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: recursively_number_parts(series_node, 1, NodeType.SERIES)
        2
        sage: series_node.comp_num
        1
        sage: series_node.children[0].comp_num
        1
        sage: series_node.children[1].comp_num
        2

    TESTS::

        sage: (series_node.comp_num == 1 and
        ....:  series_node.children[0].comp_num == 1 and
        ....:  series_node.children[1].comp_num == 2)
        True
    """
    # inner function
    def number_subtree(subtree_root, number):
        """
        set the ``comp_num`` for all the nodes in the subtree to ``number``

        INPUT:

        - ``subtree_root`` -- root of the subtree to be numbered

        - ``number`` -- number assigned to the subtree
        """
        subtree_root.comp_num = number
        if subtree_root.node_type != NodeType.NORMAL:
            for child in subtree_root.children:
                number_subtree(child, number)

    orig_part_num = part_num

    if part_root.node_type == by_type:
        # if node_type is same as tree's node_type then cocomp_num is
        # incremented before assigning to each subtree
        part_root.comp_num = part_num
        for child in part_root.children:
            number_subtree(child, part_num)
            part_num += 1
    else:
        # entire tree is numbered by cocomp_num
        number_subtree(part_root, part_num)
        part_num += 1
    return part_num - orig_part_num


def assembly(graph, root, vertex_status, vertex_dist):
    """
    Assemble the forest obtained after the promotion phase into a modular
    decomposition tree.

    INPUT:

    - ``graph`` -- graph whose MD tree is to be computed

    - ``root`` -- Forest which would be assembled into a MD tree

    - ``vertex_status`` -- Dictionary which stores the position of vertex with
      respect to the source

    - ``vertex_dist`` -- Dictionary which stores the distance of vertex from
      source vertex

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(6, 1)
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                           create_normal_node(7)]
        sage: forest.children.insert(1, series_node)
        sage: forest.children.insert(3, parallel_node)
        sage: vertex_status = {2: VertexPosition.LEFT_OF_SOURCE,
        ....:                  3: VertexPosition.SOURCE,
        ....:                  1: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  4: VertexPosition.LEFT_OF_SOURCE,
        ....:                  5: VertexPosition.LEFT_OF_SOURCE,
        ....:                  6: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  7: VertexPosition.RIGHT_OF_SOURCE}
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 3}
        sage: forest.children[0].comp_num = 1
        sage: forest.children[1].comp_num = 1
        sage: forest.children[1].children[0].comp_num = 1
        sage: forest.children[1].children[1].comp_num = 1
        sage: number_components(forest, vertex_status)
        sage: assembly(g, forest, vertex_status, vertex_dist)
        sage: forest.children
        [PRIME [NORMAL [2], SERIES [NORMAL [4], NORMAL [5]], NORMAL [3],
                PARALLEL [NORMAL [6], NORMAL [7]], NORMAL [1]]]

        sage: g.add_edge(4, 2)
        sage: g.add_edge(5, 2)
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                           create_normal_node(7)]
        sage: forest.children.insert(1, series_node)
        sage: forest.children.insert(3, parallel_node)
        sage: number_cocomponents(forest, vertex_status)
        sage: assembly(g, forest, vertex_status, vertex_dist)
        sage: forest.children
        [PRIME [NORMAL [2], SERIES [NORMAL [4], NORMAL [5], NORMAL [3]],
                PARALLEL [NORMAL [6], NORMAL [7]], NORMAL [1]]]
    """
    # Maps index to the mu computed for the (co)component at the index
    mu = {}

    # Stores index in the forest containing the source vertex
    source_index = -1

    # Maps index to list of vertices in the (co)component at the index
    vertices_in_component = {}

    # comp_num of parent should be equal to comp_num of its first child
    update_comp_num(root)

    for index, component in enumerate(root.children):

        if (component.node_type == NodeType.NORMAL and
            vertex_status[component.children[0]] == VertexPosition.SOURCE):
            source_index = root.children.index(component)

        vertices_in_component[index] = get_vertices(component)
        component.index_in_root = index

    # compute mu values for (co)components
    for index, component in enumerate(root.children):
        if index < source_index:
            mu[index] = compute_mu_for_co_component(graph, index,
                                                    source_index, root,
                                                    vertices_in_component)
        elif index > source_index:
            mu[index] = compute_mu_for_component(graph, index,
                                                 source_index, root,
                                                 vertices_in_component)

    mu[source_index] = root.children[source_index]

    # stores the leftmost cocomponent included in the module containing
    # source_index
    left = root.children[source_index]

    # stores the rightmost component included in the module containing
    # source_index
    right = root.children[source_index]

    while len(root.children) != 1:
        # source_index is changed every time a new module is formed therefore
        # updated. left or right are also updated every time module is formed.

        # First series module is attempted
        result, source_index = check_series(root, left, right,
                                              source_index, mu)
        if result:
            left = root.children[source_index].children[0]
            continue

        # If series module can't be formed, parallel is tried
        result, source_index = check_parallel(graph, root, left, right,
                                                source_index, mu, vertex_dist,
                                                vertices_in_component)
        if result:
            right = root.children[source_index].children[-1]
            continue

        # Finally a prime module is formed if both
        # series and parallel can not be created
        result, source_index = check_prime(graph, root, left, right,
                                             source_index, mu, vertex_dist,
                                             vertices_in_component)
        if result:
            if root.children[source_index].children[0].index_in_root != -1:
                left = root.children[source_index].children[0]
            if root.children[source_index].children[-1].index_in_root != -1:
                right = root.children[source_index].children[-1]


def update_comp_num(node):
    """
    Set the ``comp_num`` of ``node`` to the ``comp_num`` of its first child.

    INPUT:

    - ``node`` -- node whose comp_num needs to be updated

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.comp_num = 2
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: series_node.children[0].comp_num = 3
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                           create_normal_node(7)]
        sage: forest.children.insert(0, series_node)
        sage: forest.children.insert(3, parallel_node)
        sage: update_comp_num(forest)
        sage: series_node.comp_num
        3
        sage: forest.comp_num
        2
    """
    if node.node_type != NodeType.NORMAL:
        node.comp_num = node.children[0].comp_num
        for child in node.children:
            update_comp_num(child)


def check_prime(graph, root, left, right,
                source_index, mu, vertex_dist,
                vertices_in_component):
    """
    Assemble the forest to create a prime module.

    INPUT:

    - ``root`` -- forest which needs to be assembled

    - ``left`` -- the leftmost fragment of the last module

    - ``right`` -- the rightmost fragment of the last module

    - ``source_index`` -- index of the tree containing the source vertex

    - ``mu`` -- dictionary which maps the (co)components with their mu values

    - ``vertex_dist`` -- dictionary which stores the distance of vertex from
      source vertex

    - ``vertices_in_component`` -- dictionary which stores a list of various
      vertices in a (co)component

    OUTPUT:

    ``[module_formed, source_index]`` where ``module_formed`` is ``True`` if
    module is formed else ``False`` and ``source_index`` is the index of the
    new module which contains the source vertex

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(6, 1)
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                           create_normal_node(7)]
        sage: forest.children.insert(1, series_node)
        sage: forest.children.insert(3, parallel_node)
        sage: vertex_status = {2: VertexPosition.LEFT_OF_SOURCE,
        ....:                  3: VertexPosition.SOURCE,
        ....:                  1: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  4: VertexPosition.LEFT_OF_SOURCE,
        ....:                  5: VertexPosition.LEFT_OF_SOURCE,
        ....:                  6: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  7: VertexPosition.RIGHT_OF_SOURCE}
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 3}
        sage: source_index = 2
        sage: vertices_in_component = {}
        sage: mu = {}
        sage: left = right = forest.children[2]
        sage: for index, component in enumerate(forest.children):
        ....:     vertices_in_component[index] = get_vertices(component)
        ....:     component.index_in_root = index
        sage: for index, component in enumerate(forest.children):
        ....:     if index < source_index:
        ....:         mu[index] = compute_mu_for_co_component(g, index,
        ....:                                           source_index, forest,
        ....:                                           vertices_in_component)
        ....:     elif index > source_index:
        ....:         mu[index] = compute_mu_for_component(g, index,
        ....:                                           source_index, forest,
        ....:                                           vertices_in_component)
        sage: forest.children[0].comp_num = 1
        sage: forest.children[1].comp_num = 1
        sage: forest.children[1].children[0].comp_num = 1
        sage: forest.children[1].children[1].comp_num = 1
        sage: number_components(forest, vertex_status)
        sage: check_prime(g, forest, left, right,
        ....:              source_index, mu, vertex_dist,
        ....:              vertices_in_component)
        [True, 0]
        sage: forest.children
        [PRIME [NORMAL [2], SERIES [NORMAL [4], NORMAL [5]], NORMAL [3],
                PARALLEL [NORMAL [6], NORMAL [7]], NORMAL [1]]]
    """
    # stores the index of rightmost component included in the prime module
    new_right_index = source_index
    if source_index + 1 < len(root.children):
        new_right_index += 1

    # stores the index of leftmost component included in the prime module
    new_left_index = source_index
    if source_index >= 1:
        new_left_index -= 1

    # stores the indices of the cocomponents included in the prime module
    # the cocomponents are extracted one by one from left_queue for adding
    # more components
    left_queue = deque()

    # stores the indices of the components included in the prime module
    # the components are extracted one by one from right_queue for adding
    # more cocomponents
    right_queue = deque()

    if new_left_index != source_index:
        left_queue.append(new_left_index)
    if new_right_index != source_index:
        right_queue.append(new_right_index)

    while left_queue or right_queue:

        if left_queue:

            # cocomponent indices extracted from the queue
            left_index = left_queue.popleft()

            # more components added based on the below condition
            while (new_right_index < len(root.children) - 1 and
                   root.children[new_right_index].index_in_root < mu[left_index].index_in_root):
                new_right_index += 1
                right_queue.append(new_right_index)

            # cocomponent added while cocomponent at left_index
            # has cocomponent to its left with same comp_num
            while has_left_cocomponent_fragment(root, left_index):
                if left_index >= 1:
                    left_index -= 1
                    if new_left_index > left_index:
                        left_queue.append(left_index)
                    new_left_index = min(left_index, new_left_index)

        if right_queue:

            # component indices extracted from the queue
            right_index = right_queue.popleft()

            # more cocomponents added based on the below condition
            while (new_left_index > 0 and
                   root.children[new_left_index].index_in_root > mu[right_index].index_in_root):
                new_left_index -= 1
                left_queue.append(new_left_index)

            # component is added while component at right_index
            # has component to its right with same comp_num
            # or has a connected component with vertices at different
            # level from the source vertex
            while (has_right_component_fragment(root, right_index) or
                   has_right_layer_neighbor(graph, root,
                                            right_index, vertex_dist,
                                            vertices_in_component)):

                if has_right_layer_neighbor(graph, root,
                                            right_index, vertex_dist,
                                            vertices_in_component):
                    new_left_index = 0
                    new_right_index = len(root.children) - 1
                    break

                if right_index + 1 < len(root.children):
                    right_index += 1
                    if new_right_index < right_index:
                        right_queue.append(right_index)
                    new_right_index = max(right_index, new_right_index)

    node = create_prime_node()

    # vertices or modules are added in the prime_module
    for temp in range(new_left_index, new_right_index + 1):
        node.children.append(root.children[temp])

    # list elements included in the prime module
    # are removed from the forest
    root.children[new_left_index:new_right_index + 1] = []

    #insert the newly created prime module in the forest
    root.children.insert(new_left_index, node)

    return [True, new_left_index]


def check_parallel(graph, root, left, right,
                   source_index, mu, vertex_dist,
                   vertices_in_component):
    """
    Assemble the forest to create a parallel module.

    INPUT:

    - ``root`` -- forest which needs to be assembled

    - ``left`` -- the leftmost fragment of the last module

    - ``right`` -- the rightmost fragment of the last module

    - ``source_index`` -- index of the tree containing the source vertex

    - ``mu`` -- dictionary which maps the (co)components with their mu values

    - ``vertex_dist`` -- dictionary which stores the distance of vertex from
      source vertex

    - ``vertices_in_component`` -- dictionary which stores a list of various
      vertices in a (co)component

    OUTPUT:

    ``[module_formed, source_index]`` where ``module_formed`` is ``True`` if
    module is formed else ``False`` and ``source_index`` is the index of the
    new module which contains the source vertex

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(2, 1)
        sage: g.add_edge(4, 1)
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                   create_normal_node(7), create_normal_node(1)]
        sage: forest.children.insert(1, series_node)
        sage: forest.children.append(parallel_node)
        sage: vertex_status = {2: VertexPosition.LEFT_OF_SOURCE,
        ....:                  3: VertexPosition.SOURCE,
        ....:                  1: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  4: VertexPosition.LEFT_OF_SOURCE,
        ....:                  5: VertexPosition.LEFT_OF_SOURCE,
        ....:                  6: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  7: VertexPosition.RIGHT_OF_SOURCE}
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 2}
        sage: source_index = 2
        sage: vertices_in_component = {}
        sage: mu = {}
        sage: left = right = forest.children[2]
        sage: for index, component in enumerate(forest.children):
        ....:     vertices_in_component[index] = get_vertices(component)
        ....:     component.index_in_root = index
        sage: for index, component in enumerate(forest.children):
        ....:     if index < source_index:
        ....:         mu[index] = compute_mu_for_co_component(g, index,
        ....:                                           source_index, forest,
        ....:                                           vertices_in_component)
        ....:     elif index > source_index:
        ....:         mu[index] = compute_mu_for_component(g, index,
        ....:                                           source_index, forest,
        ....:                                           vertices_in_component)
        sage: number_components(forest, vertex_status)
        sage: check_parallel(g, forest, left, right,
        ....:              source_index, mu, vertex_dist,
        ....:              vertices_in_component)
        [True, 2]
        sage: forest.children
        [NORMAL [2],
         SERIES [NORMAL [4], NORMAL [5]],
         PARALLEL [NORMAL [3], NORMAL [6], NORMAL [7], NORMAL [1]]]
    """
    # stores the index of rightmost component included in the parallel module
    new_right_index = source_index

    while new_right_index + 1 < len(root.children):

        # component at new_right_index + 1 is added only if it doesn't have
        # a component to its right with same comp_num
        if has_right_component_fragment(root, new_right_index + 1):
            break

        # component at new_right_index + 1 is added only if it doesn't have a
        # connected component to its right with vertices at different level
        # from its vertices
        if has_right_layer_neighbor(graph, root, new_right_index + 1,
                                    vertex_dist, vertices_in_component):
            break

        # stores the index in root of new component to be added in the
        # parallel module
        i = root.children[new_right_index + 1].index_in_root

        # condition for adding more components in the parallel module
        if mu[i].index_in_root >= left.index_in_root:
            new_right_index += 1
        else:
            break

    # if new_right_index > source_index then only parallel
    # module can be formed
    if source_index != new_right_index:
        node = create_parallel_node()
        for temp in range(source_index, new_right_index + 1):

            # if module X to be included in the new parallel module Y is also
            # parallel then children of X and not X are included in Y
            if root.children[temp].node_type == NodeType.PARALLEL:
                for child in root.children[temp].children:
                    node.children.append(child)
                    child.index_in_root = root.children[temp].index_in_root
            else:
                node.children.append(root.children[temp])

        # list elements included in the parallel module are removed from the
        # forest
        root.children[source_index:new_right_index + 1] = []

        # insert the newly created parallel module into the forest
        root.children.insert(source_index, node)

        return [True, source_index]

    # no parallel module was formed
    return [False, source_index]


def check_series(root, left, right, source_index, mu):
    """
    Assemble the forest to create a series module.

    INPUT:

    - ``root`` -- forest which needs to be assembled

    - ``left`` -- The leftmost fragment of the last module

    - ``right`` -- The rightmost fragment of the last module

    - ``source_index`` -- index of the tree containing the source vertex

    - ``mu`` -- dictionary which maps the (co)components with their mu values

    - ``vertex_dist`` -- dictionary which stores the distance of vertex from
      source vertex

    - ``vertices_in_component`` -- dictionary which stores a list of various
      vertices in a (co)component

    OUTPUT:

    ``[module_formed, source_index]`` where ``module_formed`` is ``True`` if
    module is formed else ``False`` and ``source_index`` is the index of the
    new module which contains the source vertex

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(2, 1)
        sage: g.add_edge(6, 1)
        sage: g.add_edge(4, 2)
        sage: g.add_edge(5, 2)
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                           create_normal_node(7)]
        sage: forest.children.insert(1, series_node)
        sage: forest.children.insert(3, parallel_node)
        sage: vertex_status = {2: VertexPosition.LEFT_OF_SOURCE,
        ....:                  3: VertexPosition.SOURCE,
        ....:                  1: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  4: VertexPosition.LEFT_OF_SOURCE,
        ....:                  5: VertexPosition.LEFT_OF_SOURCE,
        ....:                  6: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  7: VertexPosition.RIGHT_OF_SOURCE}
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 3}
        sage: source_index = 2
        sage: vertices_in_component = {}
        sage: mu = {}
        sage: left = right = forest.children[2]
        sage: for index, component in enumerate(forest.children):
        ....:     vertices_in_component[index] = get_vertices(component)
        ....:     component.index_in_root = index
        sage: for index, component in enumerate(forest.children):
        ....:     if index < source_index:
        ....:         mu[index] = compute_mu_for_co_component(g, index,
        ....:                                           source_index, forest,
        ....:                                           vertices_in_component)
        ....:     elif index > source_index:
        ....:         mu[index] = compute_mu_for_component(g, index,
        ....:                                           source_index, forest,
        ....:                                           vertices_in_component)
        sage: number_cocomponents(forest, vertex_status)
        sage: number_components(forest, vertex_status)
        sage: check_series(forest, left, right,
        ....:              source_index, mu)
        [True, 1]
        sage: forest.children
        [NORMAL [2],
         SERIES [NORMAL [4], NORMAL [5], NORMAL [3]],
         PARALLEL [NORMAL [6], NORMAL [7]],
         NORMAL [1]]
    """
    # stores the index of leftmost component included in the parallel module
    new_left_index = source_index

    while new_left_index > 0:

        # cocomponent at new_left_index - 1 is added only if it doesn't have
        # a cocomponent to its left with same comp_num
        if has_left_cocomponent_fragment(root, new_left_index - 1):
            break

        # stores the index in root of new cocomponent to be added in the
        # series module
        i = root.children[new_left_index - 1].index_in_root

        # condition for adding more cocomponents in the series module
        if mu[i].index_in_root <= right.index_in_root:
            new_left_index -= 1
        else:
            break

    # if new_left_index < source_index then only series module can be formed
    if source_index != new_left_index:
        node = create_series_node()
        for temp in range(new_left_index, source_index + 1):

            if root.children[temp].node_type == NodeType.SERIES:
                # if module X to be included in the new series module Y is
                # also series then children of X and not X are included in Y
                for child in root.children[temp].children:
                    child.index_in_root = root.children[temp].index_in_root
                    node.children.append(child)
            else:
                node.children.append(root.children[temp])

        # list elements included in the series module
        # are removed from the forest
        root.children[new_left_index:source_index + 1] = []

        # insert the newly created series module into the forest
        root.children.insert(new_left_index, node)

        return [True, new_left_index]

    # no series module could be formed
    return [False, new_left_index]


def has_left_cocomponent_fragment(root, cocomp_index):
    """
    Check whether cocomponent at ``cocomp_index`` has a cocomponent to its left
    with same ``comp_num``.

    INPUT:

    - ``root`` -- the forest to which cocomponent belongs

    - ``cocomp_index`` -- index at which cocomponent is present in root

    OUTPUT:

    ``True`` if cocomponent at ``cocomp_index`` has a cocomponent to its left with
    same ``comp_num``, and ``False`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                           create_normal_node(7)]
        sage: forest.children.insert(1, series_node)
        sage: forest.children.insert(3, parallel_node)
        sage: forest.children[0].comp_num = 1
        sage: forest.children[1].comp_num = 1
        sage: forest.children[1].children[0].comp_num = 1
        sage: forest.children[1].children[1].comp_num = 1
        sage: has_left_cocomponent_fragment(forest, 1)
        True
        sage: has_left_cocomponent_fragment(forest, 0)
        False
    """
    return any(root.children[index].comp_num == root.children[cocomp_index].comp_num
                   for index in range(cocomp_index))

def has_right_component_fragment(root, comp_index):
    """
    Check whether component at ``comp_index`` has a component to its right with
    same ``comp_num``.

    INPUT:

    - ``root`` -- the forest to which component belongs

    - ``comp_index`` -- index at which component is present in root

    OUTPUT:

    ``True`` if component at ``comp_index`` has a component to its right with
    same ``comp_num``, and ``False`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                           create_normal_node(7)]
        sage: forest.children.insert(1, series_node)
        sage: forest.children.insert(3, parallel_node)
        sage: forest.children[3].comp_num = 1
        sage: forest.children[4].comp_num = 1
        sage: has_right_component_fragment(forest, 3)
        True
        sage: has_right_component_fragment(forest, 4)
        False
    """
    return any(root.children[index].comp_num == root.children[comp_index].comp_num
                   for index in range(comp_index + 1, len(root.children)))

def has_right_layer_neighbor(graph, root, comp_index,
                             vertex_dist, vertices_in_component):
    """
    Check whether component at ``comp_index`` has a connected component to its
    right with vertices at different level from the source vertex.

    INPUT:

    - ``root`` -- the forest to which component belongs

    - ``comp_index`` -- index at which component is present in root

    - ``vertex_dist`` -- dictionary which stores the distance of vertex from
      source vertex

    - ``vertices_in_component`` -- dictionary which stores a list of various
      vertices in a (co)component

    OUTPUT:

    ``True`` if component at ``comp_index`` has a right layer neighbor, and
    ``False`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(2, 1)
        sage: g.add_edge(6, 1)
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                           create_normal_node(7)]
        sage: forest.children.insert(1, series_node)
        sage: forest.children.insert(3, parallel_node)
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 3}
        sage: vertices_in_component = {}
        sage: for index, component in enumerate(forest.children):
        ....:     vertices_in_component[index] = get_vertices(component)
        ....:     component.index_in_root = index
        sage: has_right_layer_neighbor(g, forest, 3, vertex_dist,
        ....:                          vertices_in_component)
        True

    """
    for index in range(comp_index + 1, len(root.children)):

        # check vertex in component at index has different level from vertex
        # in component at comp_index and are connected to each other
        if ((vertex_dist[get_vertex_in(root.children[index])] >
                 vertex_dist[get_vertex_in(root.children[comp_index])]
             ) and
            (is_component_connected(graph, root.children[index].index_in_root,
                                    root.children[comp_index].index_in_root,
                                    vertices_in_component)
            )):
            return True

    return False


def get_vertex_in(node):
    """
    Return the first vertex encountered in the depth-first traversal of the
    tree rooted at node

    INPUT:

    - ``tree`` -- input modular decomposition tree

    OUTPUT:

    Return the first vertex encountered in recursion

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: forest.children.insert(1, series_node)
        sage: get_vertex_in(forest)
        2
    """
    while node.node_type != NodeType.NORMAL:
        node = node.children[0]
    return node.children[0]

def compute_mu_for_co_component(graph, component_index, source_index,
                                root, vertices_in_component):
    """
    Compute the mu value for co-component

    INPUT:

    - ``graph`` -- Graph whose MD tree needs to be computed

    - ``component_index`` -- index of the co-component

    - ``source_index`` -- index of the source in the forest

    - ``root`` -- the forest which needs to be assembled into a MD tree

    - ``vertices_in_component`` -- dictionary which maps index i to list of
      vertices in the tree at index i in the forest

    OUTPUT:

    The mu value (component in the forest) for the co-component

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(2, 1)
        sage: g.add_edge(6, 1)
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                           create_normal_node(7)]
        sage: forest.children.insert(1, series_node)
        sage: forest.children.insert(3, parallel_node)
        sage: vertices_in_component = {}
        sage: for index, component in enumerate(forest.children):
        ....:     vertices_in_component[index] = get_vertices(component)
        sage: compute_mu_for_co_component(g, 0, 2, forest,
        ....:                             vertices_in_component)
        NORMAL [1]
        sage: compute_mu_for_co_component(g, 1, 2, forest,
        ....:                             vertices_in_component)
        NORMAL [3]
    """
    for index in range(len(root.children) - 1, source_index, -1):
        if is_component_connected(graph, component_index,
                                  index, vertices_in_component):
            return root.children[index]

    # return the default value
    return root.children[source_index]

def compute_mu_for_component(graph, component_index, source_index,
                             root, vertices_in_component):
    """
    Compute the mu value for component

    INPUT:

    - ``graph`` -- Graph whose MD tree needs to be computed

    - ``component_index`` -- index of the component

    - ``source_index`` -- index of the source in the forest

    - ``root`` -- the forest which needs to be assembled into a MD tree

    - ``vertices_in_component`` -- dictionary which maps index i to list of
      vertices in the tree at the index i in the forest

    OUTPUT:

    The mu value (co-component in the forest) for the component

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(6, 1)
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                           create_normal_node(7)]
        sage: forest.children.insert(1, series_node)
        sage: forest.children.insert(3, parallel_node)
        sage: vertices_in_component = {}
        sage: for index, component in enumerate(forest.children):
        ....:     vertices_in_component[index] = get_vertices(component)
        sage: compute_mu_for_component(g, 3, 2, forest,
        ....:                          vertices_in_component)
        SERIES [NORMAL [4], NORMAL [5]]
        sage: compute_mu_for_component(g, 4, 2, forest,
        ....:                          vertices_in_component)
        NORMAL [2]
    """
    # default mu value for a component
    mu_for_component = root.children[0]

    for index in range(source_index):
        if (mu_for_component == root.children[index] and
            is_component_connected(graph, component_index,
                                   index, vertices_in_component)):
            mu_for_component = root.children[index + 1]

    # return the default value
    return mu_for_component


def is_component_connected(graph, index1, index2, vertices_in_component):
    """
    Check whether the two specified (co)components are connected.

    INPUT:

    - ``graph`` -- Graph whose MD tree needs to be computed

    - ``index1`` -- index of the first (co)component

    - ``index2`` -- index of the second (co)component

    - ``vertices_in_component`` -- dictionary which maps index i to list of
      vertices in the tree at the index i in the forest

    OUTPUT:

    ``True`` if the (co)components are connected else ``False``

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(6, 1)
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                           create_normal_node(7)]
        sage: forest.children.insert(1, series_node)
        sage: forest.children.insert(3, parallel_node)
        sage: vertices_in_component = {}
        sage: for index, component in enumerate(forest.children):
        ....:     vertices_in_component[index] = get_vertices(component)
        sage: is_component_connected(g, 0, 1, vertices_in_component)
        False
        sage: is_component_connected(g, 0, 3, vertices_in_component)
        True
    """
    V1 = vertices_in_component[index1]
    V2 = frozenset(vertices_in_component[index2])

    return any(u in V2 for v in V1 for u in graph.neighbor_iterator(v))

def get_vertices(component_root):
    """
    Compute the list of vertices in the (co)component

    INPUT:

    - ``component_root`` -- root of the (co)component whose vertices need to be
      returned as a list

    OUTPUT:

    list of vertices in the (co)component

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                           create_normal_node(7)]
        sage: forest.children.insert(1, series_node)
        sage: forest.children.insert(3, parallel_node)
        sage: get_vertices(forest)
        [2, 4, 5, 3, 6, 7, 1]
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

def promote_left(root):
    """
    Perform the promotion phase on the forest root.

    If child and parent both are marked by LEFT_SPLIT then child is removed
    and placed just before the parent

    INPUT:

    - ``root`` -- The forest which needs to be promoted

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(4, 7)
        sage: g.add_edge(2, 1)
        sage: g.add_edge(6, 1)
        sage: g.add_edge(4, 2)
        sage: g.add_edge(5, 2)
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:           create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                           create_normal_node(7)]
        sage: forest.children.insert(1, series_node)
        sage: forest.children.insert(3, parallel_node)
        sage: vertex_status = {2: VertexPosition.LEFT_OF_SOURCE,
        ....:                  3: VertexPosition.SOURCE,
        ....:                  1: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  4: VertexPosition.LEFT_OF_SOURCE,
        ....:                  5: VertexPosition.LEFT_OF_SOURCE,
        ....:                  6: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  7: VertexPosition.RIGHT_OF_SOURCE}
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 3}
        sage: x = {u for u in g.neighbor_iterator(2)
        ....:            if vertex_dist[u] != vertex_dist[2]}
        sage: maximal_subtrees_with_leaves_in_x(forest, 2, x, vertex_status,
        ....:                                   False, 0)
        sage: promote_left(forest)
        sage: forest
        FOREST [NORMAL [2], SERIES [NORMAL [4], NORMAL [5]], NORMAL [3],
                PARALLEL [NORMAL [6]], PARALLEL [NORMAL [7]],
                PARALLEL [], NORMAL [1]]
    """
    q = deque()

    # q has [parent, child] elements as parent needs to be modified
    for child in root.children:
        q.append([root, child])

    while q:

        parent, child = q.popleft()

        if child.node_type == NodeType.NORMAL:
            continue

        # stores the elements to be removed from the child
        to_remove = []

        # stores the index of child in parent list
        index = parent.children.index(child)

        for grand_child in child.children:

            # if tree and child both have LEFT_SPLIT then tree from
            # child is inserted just before child in the parent
            if grand_child.has_left_split() and child.has_left_split():
                parent.children.insert(index, grand_child)
                index += 1
                to_remove.append(grand_child)
                q.append([parent, grand_child])
            else:
                q.append([child, grand_child])

        for grand_child in to_remove:
            child.children.remove(grand_child)


def promote_right(root):
    """
    Perform the promotion phase on the forest root.

    If child and parent both are marked by RIGHT_SPLIT then child is removed
    and placed just after the parent

    INPUT:

    - ``root`` -- the forest which needs to be promoted

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(4, 7)
        sage: g.add_edge(2, 1)
        sage: g.add_edge(6, 1)
        sage: g.add_edge(4, 2)
        sage: g.add_edge(5, 2)
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                           create_normal_node(7)]
        sage: forest.children.insert(1, series_node)
        sage: forest.children.insert(3, parallel_node)
        sage: vertex_status = {2: VertexPosition.LEFT_OF_SOURCE,
        ....:                  3: VertexPosition.SOURCE,
        ....:                  1: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  4: VertexPosition.LEFT_OF_SOURCE,
        ....:                  5: VertexPosition.LEFT_OF_SOURCE,
        ....:                  6: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  7: VertexPosition.RIGHT_OF_SOURCE}
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 3}
        sage: refine(g, forest, vertex_dist, vertex_status)
        sage: promote_right(forest)
        sage: forest
        FOREST [NORMAL [2], SERIES [SERIES [NORMAL [4]], SERIES [NORMAL [5]]],
                NORMAL [3], PARALLEL [], PARALLEL [NORMAL [7]],
                PARALLEL [NORMAL [6]], NORMAL [1]]
    """
    q = deque()

    # q has [parent, child] elements as parent needs to be modified
    for child in root.children:
        q.append([root, child])

    while q:

        parent, child = q.popleft()

        if child.node_type == NodeType.NORMAL:
            continue

        # stores the elements to be removed from the child
        to_remove = []

        # stores the index of child in parent list
        index = parent.children.index(child)

        for grand_child in child.children:

            # if tree and child both have RIGHT_SPLIT then tree from
            # child is inserted just after child in the parent
            if grand_child.has_right_split() and child.has_right_split():
                parent.children.insert(index + 1, grand_child)
                to_remove.append(grand_child)
                q.append([parent, grand_child])
            else:
                q.append([child, grand_child])

        for grand_child in to_remove:
            child.children.remove(grand_child)


def promote_child(root):
    """
    Perform the promotion phase on the forest `root`.

    If marked parent has no children it is removed, if it has one child then
    it is replaced by its child

    INPUT:

    - ``root`` -- the forest which needs to be promoted

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(4, 7)
        sage: g.add_edge(2, 1)
        sage: g.add_edge(6, 1)
        sage: g.add_edge(4, 2)
        sage: g.add_edge(5, 2)
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                           create_normal_node(7)]
        sage: forest.children.insert(1, series_node)
        sage: forest.children.insert(3, parallel_node)
        sage: vertex_status = {2: VertexPosition.LEFT_OF_SOURCE,
        ....:                  3: VertexPosition.SOURCE,
        ....:                  1: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  4: VertexPosition.LEFT_OF_SOURCE,
        ....:                  5: VertexPosition.LEFT_OF_SOURCE,
        ....:                  6: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  7: VertexPosition.RIGHT_OF_SOURCE}
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 3}
        sage: refine(g, forest, vertex_dist, vertex_status)
        sage: promote_right(forest)
        sage: promote_child(forest)
        sage: forest
        FOREST [NORMAL [2], SERIES [NORMAL [4], NORMAL [5]], NORMAL [3],
                NORMAL [7], NORMAL [6], NORMAL [1]]
    """
    q = deque()

    # q has [parent, child] elements as parent needs to be modified
    for child in root.children:
        q.append([root, child])

    while q:

        parent, child = q.popleft()

        if child.node_type == NodeType.NORMAL:
            continue

        # if child node itself has only one child
        if (len(child.children) == 1 and
            (child.node_split != NodeSplit.NO_SPLIT or
             child.node_type == NodeType.FOREST)):
            # replace child node by its own child

            grand_child = child.children[0]
            index = parent.children.index(child)
            parent.children.insert(index, grand_child)
            parent.children.remove(child)
            q.append([parent, grand_child])
        # if child node has no children
        elif ((not child.children) and child.node_split != NodeSplit.NO_SPLIT):
            # remove the child node
            parent.children.remove(child)
        else:
            for grand_child in child.children:
                q.append([child, grand_child])


def clear_node_split_info(root):
    """
    Set the node_split of nodes to NO_SPLIT

    INPUT:

    - ``root`` -- the forest which needs to be cleared of split information

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: series_node.children[0].node_split = NodeSplit.LEFT_SPLIT
        sage: series_node.node_split = NodeSplit.RIGHT_SPLIT
        sage: forest.children.insert(1, series_node)
        sage: clear_node_split_info(forest)
        sage: series_node.node_split == NodeSplit.NO_SPLIT
        True
        sage: series_node.children[0].node_split == NodeSplit.NO_SPLIT
        True
    """
    root.node_split = NodeSplit.NO_SPLIT

    if root.node_type != NodeType.NORMAL:
        for node in root.children:
            clear_node_split_info(node)


def refine(graph, root, vertex_dist, vertex_status):
    """
    Refine the forest based on the active edges

    INPUT:

    - ``graph`` -- graph whose MD tree needs to be computed

    - ``root`` -- the forest which needs to be assembled into a MD tree

    - ``vertex_dist`` -- dictionary mapping the vertex with distance from the
      source

    - ``vertex_status`` -- dictionary mapping the vertex to the position w.r.t.
      source

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(2, 1)
        sage: g.add_edge(6, 1)
        sage: g.add_edge(4, 2)
        sage: g.add_edge(5, 2)
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                           create_normal_node(7)]
        sage: forest.children.insert(1, series_node)
        sage: forest.children.insert(3, parallel_node)
        sage: vertex_status = {2: VertexPosition.LEFT_OF_SOURCE,
        ....:                  3: VertexPosition.SOURCE,
        ....:                  1: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  4: VertexPosition.LEFT_OF_SOURCE,
        ....:                  5: VertexPosition.LEFT_OF_SOURCE,
        ....:                  6: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  7: VertexPosition.RIGHT_OF_SOURCE}
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 3}
        sage: refine(g, forest, vertex_dist, vertex_status)
        sage: forest
        FOREST [NORMAL [2], SERIES [NORMAL [4], NORMAL [5]], NORMAL [3],
               PARALLEL [PARALLEL [NORMAL [6]], PARALLEL [NORMAL [7]]],
               NORMAL [1]]
    """
    x_used = []

    # active edges of each vertex in the graph is used to refine the forest
    for v in graph.vertices():
        if v in vertex_status and vertex_status[v] == VertexPosition.SOURCE:
            continue

        # set of vertices connected through active edges to v
        x = {u for u in graph.neighbor_iterator(v)
                        if vertex_dist[u] != vertex_dist[v]}

        if x not in x_used:
            x_used.append(x)
            maximal_subtrees_with_leaves_in_x(root, v, x,
                                              vertex_status, False, 0)

    get_child_splits(root)


def get_child_splits(root):
    """
    Add the node_split of children to the parent node

    INPUT:

    - ``root`` -- input modular decomposition tree

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: series_node.children[0].node_split = NodeSplit.LEFT_SPLIT
        sage: series_node.node_split = NodeSplit.RIGHT_SPLIT
        sage: forest.children.insert(1, series_node)
        sage: get_child_splits(forest)
        sage: series_node.node_split == NodeSplit.BOTH_SPLIT
        True
        sage: forest.node_split == NodeSplit.BOTH_SPLIT
        True
    """
    if root.node_type != NodeType.NORMAL:
        for node in root.children:
            get_child_splits(node)
            root.set_node_split(node.node_split)


def maximal_subtrees_with_leaves_in_x(root, v, x, vertex_status,
                                      tree_left_of_source, level):
    """
    Refine the forest based on the active edges(x) of vertex v

    INPUT:

    - ``root`` -- the forest which needs to be assembled into a MD tree

    - ``v`` -- the vertex used to refine

    - ``x`` -- set of vertices connected to v and at different distance
      from source compared to v

    - ``vertex_status`` -- dictionary mapping the vertex to the position
      w.r.t source

    - ``tree_left_of_source`` -- flag indicating whether tree is

    - ``level`` -- indicates the recursion level, 0 for root

    OUTPUT:

    ``[contained_in_x, split]`` where ``contained_in_x`` is ``True`` if all
    vertices in root are subset of x else ``False`` and ``split`` is the
    split which occurred at any node in root

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = Graph()
        sage: g.add_vertices([1, 2, 3, 4, 5, 6, 7])
        sage: g.add_edge(2, 3)
        sage: g.add_edge(4, 3)
        sage: g.add_edge(5, 3)
        sage: g.add_edge(2, 6)
        sage: g.add_edge(2, 7)
        sage: g.add_edge(2, 1)
        sage: g.add_edge(6, 1)
        sage: g.add_edge(4, 2)
        sage: g.add_edge(5, 2)
        sage: forest = Node(NodeType.FOREST)
        sage: forest.children = [create_normal_node(2),
        ....:                    create_normal_node(3), create_normal_node(1)]
        sage: series_node = Node(NodeType.SERIES)
        sage: series_node.children = [create_normal_node(4),
        ....:                         create_normal_node(5)]
        sage: parallel_node = Node(NodeType.PARALLEL)
        sage: parallel_node.children = [create_normal_node(6),
        ....:                           create_normal_node(7)]
        sage: forest.children.insert(1, series_node)
        sage: forest.children.insert(3, parallel_node)
        sage: vertex_status = {2: VertexPosition.LEFT_OF_SOURCE,
        ....:                  3: VertexPosition.SOURCE,
        ....:                  1: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  4: VertexPosition.LEFT_OF_SOURCE,
        ....:                  5: VertexPosition.LEFT_OF_SOURCE,
        ....:                  6: VertexPosition.RIGHT_OF_SOURCE,
        ....:                  7: VertexPosition.RIGHT_OF_SOURCE}
        sage: vertex_dist = {2: 1, 4: 1, 5: 1, 3: 0, 6: 2, 7: 2, 1: 3}
        sage: x = {u for u in g.neighbor_iterator(2)
        ....:            if vertex_dist[u] != vertex_dist[2]}
        sage: maximal_subtrees_with_leaves_in_x(forest, 2, x, vertex_status,
        ....:                                   False, 0)
        sage: forest
        FOREST [NORMAL [2], SERIES [NORMAL [4], NORMAL [5]], NORMAL [3],
                PARALLEL [NORMAL [6], NORMAL [7]], NORMAL [1]]
        sage: x = {u for u in g.neighbor_iterator(1)
        ....:            if vertex_dist[u] != vertex_dist[1]}
        sage: maximal_subtrees_with_leaves_in_x(forest, 1, x, vertex_status,
        ....:                                   False, 0)
        sage: forest
        FOREST [NORMAL [2], SERIES [NORMAL [4], NORMAL [5]], NORMAL [3],
                PARALLEL [PARALLEL [NORMAL [6]], PARALLEL [NORMAL [7]]],
                NORMAL [1]]
    """
    def update_node_info(node, node_type, node_split, comp_num, subtree_list):
        """
        Set the various fields for a tree node and update its subtrees

        - ``node`` -- node whose fields need to be updated
        - ``node_type`` -- node_type to be set
        - ``node_split`` -- node_split to be set
        - ``comp_num`` -- comp_num to be set
        - ``subtree_list`` -- list containing the subtrees

        """
        node.node_type = node_type
        node.node_split = node_split
        node.comp_num = comp_num
        node.children = subtree_list

    return_split = NodeSplit.NO_SPLIT  # initialize split to NO_SPLIT

    # all trees in a forest are refined using x
    if root.node_type == NodeType.FOREST:

        # indicates whether tree is left of source, True if left of source
        left_flag = True

        for node in root.children:
            if (node.node_type == NodeType.NORMAL and
                node.children[0] in vertex_status and
                vertex_status[node.children[0]] == VertexPosition.SOURCE):
                left_flag = False
            subtree_result = maximal_subtrees_with_leaves_in_x(node, v, x,
                                                               vertex_status,
                                                               left_flag,
                                                               level)
            if subtree_result:
                # Mark the ancestors
                root.set_node_split(subtree_result[1])

    # handles the prime, series and parallel cases
    elif root.node_type != NodeType.NORMAL:

        flag = True  # indicates the entire root is contained in x
        split_flag = False  # indicates a split is required
        Ta = []  # contains subtrees with leaves in x
        Tb = []  # contains subtrees with leaves not in x

        for node in root.children:

            # refines the children of root
            subtree_result = maximal_subtrees_with_leaves_in_x(node, v, x,
                                                               vertex_status,
                                                               tree_left_of_source,
                                                               level + 1)

            if subtree_result:
                flag = flag and subtree_result[0]

                # add the node split of children to root
                root.set_node_split(subtree_result[1])

                if subtree_result[0]:
                    Ta.append(node)
                    split_flag = True
                else:
                    Tb.append(node)

        if root.node_type == NodeType.PRIME:
            # mark all the children of prime nodes
            for node in root.children:
                node.set_node_split(root.node_split)

        if flag:
            # return if all subtrees are in x, no split required
            return [True, root.node_split]
        elif split_flag:  # split required`

            split = NodeSplit.LEFT_SPLIT

            # if v is right of source and tree is also right of source then
            # RIGHT_SPLIT
            if (vertex_status[v] == VertexPosition.RIGHT_OF_SOURCE and
                not tree_left_of_source):
                split = NodeSplit.RIGHT_SPLIT

            # add the split to root node_split
            root.set_node_split(split)

            if root.node_type == NodeType.PRIME:
                # mark all the children of prime nodes
                for node in root.children:
                    node.set_node_split(split)
                return [False, split]

            if root.is_separated:
                # if root has already been split then further split not
                # required
                return [flag, root.node_split]

            node_type = root.node_type
            root.is_separated = True

            # root[1] would now contain Ta and Tb
            root.children = []

            # add two nodes for Ta and Tb
            a = create_parallel_node()
            update_node_info(a, node_type, root.node_split,
                             Ta[0].comp_num, Ta)
            b = create_parallel_node()
            update_node_info(b, node_type, root.node_split,
                             Tb[0].comp_num, Tb)
            root.children.append(a)
            root.children.append(b)

        return_split = root.node_split
        return [flag, return_split]
    # root is a vertex and is contained in x
    elif root.children[0] in x:
        return [True, root.node_split]
    # root is a vertex and is not contained in x
    else:
        return [False, root.node_split]


def create_prime_node():
    """
    Return a prime node with no children

    OUTPUT:

    A node object with node_type set as NodeType.PRIME

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import create_prime_node
        sage: node = create_prime_node()
        sage: node
        PRIME []
    """
    return Node(NodeType.PRIME)


def create_parallel_node():
    """
    Return a parallel node with no children

    OUTPUT:

    A node object with node_type set as NodeType.PARALLEL

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import create_parallel_node
        sage: node = create_parallel_node()
        sage: node
        PARALLEL []
    """
    return Node(NodeType.PARALLEL)


def create_series_node():
    """
    Return a series node with no children

    OUTPUT:

    A node object with node_type set as NodeType.SERIES

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import create_series_node
        sage: node = create_series_node()
        sage: node
        SERIES []
    """
    return Node(NodeType.SERIES)


def create_normal_node(vertex):
    """
    Return a normal node with no children

    INPUT:

    - ``vertex`` -- vertex number

    OUTPUT:

    A node object representing the vertex with node_type set as NodeType.NORMAL

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import create_normal_node
        sage: node = create_normal_node(2)
        sage: node
        NORMAL [2]
    """
    node = Node(NodeType.NORMAL)
    node.children.append(vertex)
    return node

def print_md_tree(root):
    """
    Print the modular decomposition tree

    INPUT:

    - ``root`` -- root of the modular decomposition tree

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: print_md_tree(modular_decomposition(graphs.IcosahedralGraph()))
        PRIME
         5
         7
         11
         1
         8
         0
         9
         4
         10
         6
         2
         3
    """

    def recursive_print_md_tree(root, level):
        """
        Print the modular decomposition tree at root

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
            print("{}{}".format(level,str(root.children[0])))

    recursive_print_md_tree(root, "")



#==============================================================================
#  Habib Maurer algorithm
#==============================================================================

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
        for component in neighborhood.complement().connected_components():
            v1 = component[0]
            e = frozenset([v1, v])
            for vi in component[1:]:
                ei = frozenset([vi, v])
                pieces.union(e, ei)
    return {frozenset(chain.from_iterable(loe)): loe for loe in pieces}


def habib_maurer_algorithm(graph, g_classes=None):
    """
    Compute the modular decomposition by the algorithm of Habib and Maurer

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

    INPUT:

    - ``graph`` -- the graph for which modular decomposition tree needs to be
      computed

    - ``g_classes`` -- dictionary (default: ``None``); a dictionary whose values
      are the gamma classes of the graph, and whose keys are a frozenset of the
      vertices corresponding to the class. Used internally.

    OUTPUT:

    The modular decomposition tree of the graph.

    EXAMPLES:

    The Icosahedral graph is Prime::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: print_md_tree(habib_maurer_algorithm(graphs.IcosahedralGraph()))
        PRIME
         1
         5
         7
         8
         11
         0
         2
         6
         3
         9
         4
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

    Graph from the :wikipedia:`Modular_decomposition`::

        sage: d2 = {1:[2,3,4], 2:[1,4,5,6,7], 3:[1,4,5,6,7], 4:[1,2,3,5,6,7],
        ....:       5:[2,3,4,6,7], 6:[2,3,4,5,8,9,10,11],
        ....:       7:[2,3,4,5,8,9,10,11], 8:[6,7,9,10,11], 9:[6,7,8,10,11],
        ....:       10:[6,7,8,9], 11:[6,7,8,9]}
        sage: g = Graph(d2)
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

    Bad Input::

        sage: g = DiGraph()
        sage: habib_maurer_algorithm(g)
        Traceback (most recent call last):
        ...
        ValueError: Graph must be undirected

    Empty Graph is Prime::

        sage: g = Graph()
        sage: habib_maurer_algorithm(g)
        PRIME []


    Ensure that a random graph and an isomorphic graph have identical modular
    decompositions. ::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import permute_decomposition
        sage: permute_decomposition(2, habib_maurer_algorithm, 20, 0.5)
    """
    if graph.is_directed():
        raise ValueError("Graph must be undirected")

    if not graph.order():
        return create_prime_node()

    if graph.order() == 1:
        root = create_normal_node(next(graph.vertex_iterator()))
        return root

    elif not graph.is_connected():
        root = create_parallel_node()
        root.children = [habib_maurer_algorithm(graph.subgraph(vertices=sg), g_classes)
                             for sg in graph.connected_components()]
        return root

    g_comp = graph.complement()
    if g_comp.is_connected():
        from collections import defaultdict
        root = create_prime_node()
        if g_classes is None:
            g_classes = gamma_classes(graph)
        vertex_set = frozenset(graph)
        edges = [tuple(e) for e in g_classes[vertex_set] ]
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

    root = create_series_node()
    root.children = [habib_maurer_algorithm(graph.subgraph(vertices=sg), g_classes)
                     for sg in g_comp.connected_components()]
    return root


#=============================================================================
# Below functions are implemented to test the modular decomposition tree
#=============================================================================

# Function implemented for testing
def test_modular_decomposition(tree_root, graph):
    """
    Test the input modular decomposition tree using recursion.

    INPUT:

    - ``tree_root`` -- root of the modular decomposition tree to be tested

    - ``graph`` -- graph whose modular decomposition tree needs to be tested

    OUTPUT:

    ``True`` if input tree is a modular decomposition else ``False``

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
                    if ((get_module_type(graph.subgraph(module_formed[1])) ==
                             tree_root.node_type
                         ) and
                        (tree_root.node_type == NodeType.PARALLEL or
                             tree_root.node_type == NodeType.SERIES
                        )):
                        continue
                    return False
    return True


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
    Test whether input module is actually a module

    INPUT:

    - ``module`` -- module which needs to be tested

    - ``graph`` -- input sage graph which contains the module

    OUTPUT:

    ``True`` if input module is a module by definition else ``False``

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
    Check whether the node type of the childrens of ``module`` is ``node_type``.

    INPUT:

    - ``module`` -- module which is tested

    - ``node_type`` -- input node_type

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: g = graphs.OctahedralGraph()
        sage: tree_root = modular_decomposition(g)
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

    OUTPUT:

    The root node of a modular decomposition tree.

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
        return create_normal_node(nest)

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
    Relabel the leaves of a tree according to a dictionary

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
    from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
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
        perm = {i: perm(i) for i in get_vertices(root) }

    else:
        raise TypeError("type of perm is not supported for relabeling")

    if root.node_type == NodeType.NORMAL:
        return create_normal_node(perm[root.children[0]])
    else:
        new_root = Node(root.node_type)
        new_root.children = [relabel_tree(child, perm) for child in root.children]
        return new_root


#==============================================================================
#   Random tests
#==============================================================================

from sage.misc.random_testing import random_testing
@random_testing
def test_gamma_modules(trials, vertices, prob, verbose=False):
    r"""
    Verify that the vertices of each gamma class of a random graph are modules
    of that graph.

    INPUT:

    - ``trials`` -- the number of trials to run

    - ``vertices`` -- the size of the graph to use

    - ``prob`` -- the probability that any given edge is in the graph.
      See :meth:`~sage.graphs.generators.random.RandomGNP` for more details.

    - ``verbose`` -- print information on each trial.

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
                    assert(either_connected_or_not_connected(v, m_list, g))
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
        assert(test_modular_decomposition(t1, g1))
        assert(test_modular_decomposition(t2, g2))
        t1p = relabel_tree(t1, random_perm)
        assert(equivalent_trees(t1p, t2))
        if verbose:
            print("Passses!")

def random_md_tree(max_depth, max_fan_out, leaf_probability):
    r"""
    Create a random MD tree.

    INPUT:

    - ``max_depth`` -- the maximum depth of the tree.

    - ``max_fan_out`` -- the maximum number of children a node can have
      (must be >=4 as a prime node must have at least 4 vertices).

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

    #internal function
    def rand_md_tree(max_depth, parent_type):
        r"""
        Create the subtrees of a node.

        A child of a node cannot have the same node type as its parent if its
        parent's node type is either PARALLEL or SERIES.  Also its ``max_depth``
        is one less than its parent's.
        """
        if random() < leaf_probability or max_depth == 1:
            root = create_normal_node(current_leaf[0])
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

def md_tree_to_graph(root):
    r"""
    Create a graph having the given MD tree.

    For the prime nodes we use that every path of length 4 or more is prime.

    TODO: accept a function that generates prime graphs as a parameter and
    use that in the prime nodes.

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.modular_decomposition import *
        sage: tup1 = (NodeType.PRIME, 1, (NodeType.SERIES, 2, 3),
        ....:        (NodeType.PARALLEL, 4, 5), 6)
        sage: tree1 = nested_tuple_to_tree(tup1)
        sage: g1 = md_tree_to_graph(tree1)
        sage: g2 = Graph({1: [2, 3], 2: [1, 3, 4, 5], 3: [1, 2, 4, 5],
        ....:             4: [2, 3, 6], 5: [2, 3, 6], 6: [4, 5]})
        sage: g1.is_isomorphic(g2)
        True
    """
    from itertools import product, combinations
    from sage.graphs.graph import Graph
    def tree_to_vertices_and_edges(root):
        r"""
        Give the list of vertices and edges of the graph having the given md tree.
        """
        if root.node_type == NodeType.NORMAL:
            return (root.children, [])
        children_ve = [tree_to_vertices_and_edges(child) for child in root.children]
        vertices = [v for vs, es in children_ve for v in vs]
        edges = [e for vs, es in children_ve for e in es]
        vertex_lists = [vs for vs, es in children_ve]
        if root.node_type == NodeType.PRIME:
            for vs1, vs2 in zip(vertex_lists, vertex_lists[1:]):
                for v1, v2 in product(vs1, vs2):
                    edges.append((v1, v2))
        elif root.node_type == NodeType.SERIES:
            for vs1, vs2 in combinations(vertex_lists, 2):
                for v1, v2 in product(vs1, vs2):
                    edges.append((v1, v2))
        return (vertices, edges)

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
        assert(equivalent_trees(rand_tree, reconstruction))
        if verbose:
            print("Passes!")
