r"""
Cographs

A cograph is a `P_4`-free graph, that is a graph without induced path of order
4. Any cograph may be constructed, starting from the single vertex graph, by a
sequence of join and disjoint union operations. See the :wikipedia:`Cograph` for
more details on this graph class, and :oeis:`A000084` to know the number of
cographs of order `n \geq 1`.

This module implements the following methods concerning cographs:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`cographs` | Return an iterator over the cographs of order `n`.

Methods
-------
"""
# ****************************************************************************
#       Copyright (C) 2017-2024 Marianna Spyrakou
#                          2024 David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.combinat.partitions import AccelAsc_next


class CoTree:
    r"""
    Generic cotree node.

    This data structure is used for the generation of cographs in
    :meth:`cographs`.
    """
    def __init__(self, name='root'):
        r"""
        Initialize a cotree.

        INPUT:

        - ``name`` -- either an operation ('U' or 'J') or the size of the
          subtree rooted at this node

        EXAMPLES::

            sage: from sage.graphs.cographs import CoTree
            sage: CoTree(1)
            ( 1 )
            sage: CoTree()
            ( root )
        """
        self.name = name
        self.children = []
        self.info = None
        self.parent = None

    def __str__(self):
        r"""
        Return a string representation of ``self``.

        The subtree of a node is inside brackets and starts with the operation
        to perform between the children ('J' for join or 'U' for disjoint
        union).

        EXAMPLES::

            sage: from sage.graphs.cographs import CoTree
            sage: CoTree(1)
            ( 1 )
            sage: CoTree('J')
            [ J ]
            sage: next(graphs.cographs(4, as_graph=False))  # indirect doctest
            [ J ( 0 ) ( 1 ) ( 2 ) ( 3 ) ]
        """
        first = '[' if self.name in ['J', 'U'] else '('
        last = ' ]' if self.name in ['J', 'U'] else ' )'
        s = f"{first} {self.name}"
        for child in self.children:
            s += f' {child}'
        s += last
        return s

    __repr__ = __str__

    def add_child(self, node):
        r"""
        Add cotree ``node`` in the list of children of ``self``.

        INPUT:

        - ``node`` -- a CoTree

        EXAMPLES::

            sage: from sage.graphs.cographs import CoTree
            sage: T = CoTree('J')
            sage: T.add_child(CoTree(1))
            sage: T
            [ J ( 1 ) ]
        """
        assert isinstance(node, CoTree)
        self.children.append(node)
        node.parent = self

    def copy_tree(self, T):
        r"""
        Make `T` a copy of ``self``.

        INPUT:

        - ``T`` -- a CoTree

        EXAMPLES::

            sage: from sage.graphs.cographs import CoTree
            sage: T = CoTree('J')
            sage: T.add_child(CoTree(1))
            sage: T.add_child(CoTree(2))
            sage: T
            [ J ( 1 ) ( 2 ) ]
            sage: B =  CoTree('U')
            sage: T.copy_tree(B)
            sage: B
            [ U ( 1 ) ( 2 ) ]
        """
        for i, child in enumerate(self.children):
            T.add_child(CoTree(child.name))
            child.copy_tree(T.children[i])

    def reset_info(self):
        r"""
        Reset parameter ``info`` from all nodes of ``self``.

        EXAMPLES::

            sage: from sage.graphs.cographs import CoTree
            sage: T = CoTree(1)
            sage: B = CoTree(2)
            sage: T.add_child(B)
            sage: C = CoTree(3)
            sage: B.add_child(C)
            sage: C.info = 'info'
            sage: T.reset_info()
            sage: C.info is None
            True
        """
        for child in self.children:
            child.reset_info()
        self.info = None


def rebuild_node(u, P):
    r"""
    Replace the subtree rooted at `u` by a subtree induced by partition `P`.

    This is a helper method to method :meth:`cographs`.

    INPUT:

    - ``u`` -- a ``CoTree``

    - ``P`` -- a partition encoding the new structure of the children of `u`

    EXAMPLES::

        sage: next(graphs.cographs(3, as_graph=True)).vertices()  # indirect doctest
        [0, 1, 2]
    """
    if P is None:
        print('P is None')
    u.children = []  # delete the subtree rooted at u
    for value in P:
        this_child = CoTree(value)
        u.add_child(this_child)
        if value > 1:
            for j in range(value):
                this_child.add_child(CoTree(1))


def find_pivot(T):
    r"""
    Search for a pivot node in `T`.

    This is a helper method to method :meth:`cographs`.

    A node in `T` is a ``pivot`` if it is not a leaf, it does not induce a
    maximum partition and it is the first such node in the inverse postorder
    traversal.

    INPUT:

    - ``T`` -- a ``CoTree``

    EXAMPLES::

        sage: next(graphs.cographs(3, as_graph=True)).vertices()  # indirect doctest
        [0, 1, 2]
    """
    for child in reversed(T.children):
        pivot = find_pivot(child)
        if pivot is not None:
            return pivot

    # Check if T is a pivot
    i = T.name
    if (i != 1 and ((i//2 != T.children[0].name) or
                    (i//2 + i % 2 != T.children[1].name))):
        T.info = 'p'  # pivot mark
        return T
    return None


def next_tree(T):
    r"""
    Check if there is another tree after `T`.

    This is a helper method to method :meth:`cographs`.

    This methods returns ``True`` if there is a tree after `T`, and if so it
    modifies the input tree `T` that becomes the next tree.

    INPUT:

    - ``T`` -- a ``CoTree``

    EXAMPLES::

        sage: next(graphs.cographs(3, as_graph=True)).vertices()  # indirect doctest
        [0, 1, 2]
    """
    pivot = find_pivot(T)
    if pivot is None:
        return False

    # Find the next partition induced by the subtree pivot.
    # Partitions are represented in ascending order (i.e., `P_i \leq P_{i+1}`)
    # and we search for the next partition in lexicographic order.
    partition = [c.name for c in pivot.children]
    P = AccelAsc_next(partition)
    # and rebuild the subtree of pivot accordingly
    rebuild_node(pivot, P)

    x = pivot
    while x.parent is not None:
        ancestor = x.parent
        # Modify the bigger siblings of x (i.e., the siblings placed after
        # x in the list of children of ancestor)
        is_bigger_sibling = False
        for y, this_child in enumerate(ancestor.children):
            if this_child.info == 'p':
                is_bigger_sibling = True
            elif is_bigger_sibling:  # true only for bigger siblings of x
                if x.name == this_child.name:
                    temp = CoTree(x.name)
                    x.copy_tree(temp)
                    ancestor.children[y] = temp  # copy subtree T(x) in T(y)
                    temp.parent = ancestor
                else:
                    P = [1] * ancestor.children[y].name
                    rebuild_node(ancestor.children[y], P)

        # Make the parent of x (if any) the new pivot
        x.info = None  # reset the pivot mark
        ancestor.info = 'p'  # parent gets the pivot mark
        x = ancestor

    return True


def cographs(n, as_graph=True, immutable=False):
    r"""
    Return an iterator over the cographs of order `n`.

    A cograph is a `P_4`-free graph, that is a graph without induced path of
    order 4. Any cograph may be constructed, starting from the single vertex
    graph, by a sequence of :meth:`sage.graphs.graph.Graph.join` and
    :meth:`sage.graphs.generic_graph.GenericGraph.disjoint_union` operations.
    See the :wikipedia:`Cograph` for more details.

    This method implements the generator of all cographs of order `n` proposed
    in [JPD2018]_. The algorithm generates one by one every cotree with `n`
    nodes, and each cotree is generated by using its previous cotree. The time
    to construct the first cotree is `O(n)` and the time spent between two
    consecutive outputs is `O(n)`. Hence, the overall complexity of the
    algorithm is `O(n*M_n)`, where `n` is the number of nodes and `M_n` is the
    total number of cographs with `n` nodes (see :oeis:`A000084`).

    INPUT:

    - ``n`` -- integer larger or equal to 1

    - ``as_graph`` -- boolean (default: ``True``); whether to return graphs or
      the tree data structure encoding the graph

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph. This parameter is used only when ``as_graph
      is True``.

    EXAMPLES:

    The output can be either cotrees or graphs::

        sage: for t in graphs.cographs(3, as_graph=False):
        ....:     print(t)
        [ J ( 0 ) ( 1 ) ( 2 ) ]
        [ J [ U ( 0 ) ( 1 ) ( 2 ) ] ]
        [ J ( 0 ) [ U ( 1 ) ( 2 ) ] ]
        [ J [ U ( 0 ) [ J ( 1 ) ( 2 ) ] ] ]
        sage: for g in graphs.cographs(3, as_graph=True):
        ....:     print(g.edges(labels=False, sort=True))
        [(0, 1), (0, 2), (1, 2)]
        []
        [(0, 1), (0, 2)]
        [(1, 2)]

    Check that we agree with :oeis:`A000084`::

        sage: [sum(1 for _ in graphs.cographs(n, as_graph=False))
        ....:  for n in range(1, 8)]
        [1, 2, 4, 10, 24, 66, 180]

    TESTS::

        sage: g = next(graphs.cographs(2, as_graph=True, immutable=False))
        sage: g.is_immutable()
        False
        sage: g = next(graphs.cographs(2, as_graph=True, immutable=True))
        sage: g.is_immutable()
        True
        sage: next(graphs.cographs(0))
        Traceback (most recent call last):
        ...
        ValueError: parameter n must be at least >= 1
    """
    if n < 1:
        raise ValueError('parameter n must be at least >= 1')
    if as_graph:
        def func(T):
            return tree_to_graph(T, immutable=immutable)
    else:
        def func(T):
            return T

    if n == 1:
        T = CoTree('J')
        B = CoTree('U')
        B.add_child(CoTree(0))
        T.add_child(B)
        yield func(T)
        return

    # Construct the minimum tree
    T = CoTree(n)
    for j in range(n):
        T.add_child(CoTree(1))

    while True:
        # T corresponds to 2 cotrees: one with 'U' root and one with 'J' root
        tree_J = CoTree(T.name)
        T.copy_tree(tree_J)
        tree_U = CoTree(T.name)
        T.copy_tree(tree_U)
        # tree_J has root 'J'
        change_label(tree_J, True, [0])
        # tree 0 has root 'U'
        change_label(tree_U, False, [0])
        tree_UU = CoTree('J')
        tree_UU.add_child(tree_U)
        yield func(tree_J)
        yield func(tree_UU)
        if not next_tree(T):
            break


def change_label(tree, status, counter):
    r"""
    Set the names of the nodes of ``tree``.

    This is a helper method to method :meth:`cographs`.

    The input ``tree`` is such that each node has as label its number of
    children.  This method changes the label of each node so that a parallel
    node is labeled 'U', a series node is labeled 'J' and a leaf node gets a
    unique number.

    INPUT:

    - ``tree`` -- the tree to relabel

    - ``status`` -- boolean; used to to detect series (``True``) and parallel
      (``False``) internal nodes

    - ``counter`` -- list; the first integer of the list is used to assign a
      unique number to the leaves of the tree

    EXAMPLES::

        sage: next(graphs.cographs(4, as_graph=True)).vertices()  # indirect doctest
        [0, 1, 2, 3]
    """
    if tree.name != 1:
        tree.name = 'J' if status else 'U'
    else:
        tree.name = counter[0]
        counter[0] += 1
    for child in tree.children:
        if child is not None:
            change_label(child, not status, counter)


def tree_to_graph(tree, immutable=False):
    r"""
    Return the cograph represented by ``tree``.

    This is a helper method to method :meth:`cographs`.

    EXAMPLES::

        sage: for t in graphs.cographs(3, as_graph=True):  # indirect doctest
        ....:     print(t.edges(labels=False, sort=True))
        [(0, 1), (0, 2), (1, 2)]
        []
        [(0, 1), (0, 2)]
        [(1, 2)]
    """
    from sage.graphs.graph import Graph
    g = Graph()
    _tree_to_graph_rec(tree, g)
    return g.copy(immutable=True) if immutable else g


def _tree_to_graph_rec(tree, g):
    r"""
    Add recursively one by one the vertices and edges of ``tree`` to ``g``.

    This is a helper method to method :meth:`tree_to_graph`.

    EXAMPLES::

        sage: for t in graphs.cographs(3, as_graph=True):  # indirect doctest
        ....:     print(t.edges(labels=False, sort=True))
        [(0, 1), (0, 2), (1, 2)]
        []
        [(0, 1), (0, 2)]
        [(1, 2)]
    """
    for child in tree.children:
        _tree_to_graph_rec(child, g)
    if tree.name not in ['J', 'U']:
        tree.info = 'v'
        g.add_vertex(tree.name)
        _find_neighbors(tree, g)
        tree.reset_info()


def _find_neighbors(tree, g):
    r"""
    Identify the neighbors of node ``tree`` in ``g`` and add needed edges.

    This is a helper method to method :meth:`tree_to_graph`.

    EXAMPLES::

        sage: for t in graphs.cographs(3, as_graph=True):  # indirect doctest
        ....:     print(t.edges(labels=False, sort=True))
        [(0, 1), (0, 2), (1, 2)]
        []
        [(0, 1), (0, 2)]
        [(1, 2)]
    """
    ancestor = tree.parent
    ancestor.info = 'v'
    while ancestor is not None:
        if ancestor.name == 'J':
            for sibling in ancestor.children:
                if sibling != tree and sibling.info != 'v':
                    _add_edge(tree, sibling, g)
        elif ancestor.name == 'U':
            ancestor.info = 'v'
        ancestor = ancestor.parent


def _add_edge(u, v, g):
    r"""
    Add an edge in `g` between the nodes associated to the cotrees `u` and `v`.

    This is a helper method to method :meth:`tree_to_graph`.

    If `v` is not a leaf, then the method add edges between `u` and all leaves
    of the subtree rooted at `v`.

    The method assumes that `u` is a leaf (not tested).

    EXAMPLES::

        sage: for t in graphs.cographs(3, as_graph=True):  # indirect doctest
        ....:     print(t.edges(labels=False, sort=True))
        [(0, 1), (0, 2), (1, 2)]
        []
        [(0, 1), (0, 2)]
        [(1, 2)]
    """
    if v.name in ['J', 'U']:
        for child in v.children:
            _add_edge(u, child, g)
    else:
        g.add_edge(u.name, v.name)
