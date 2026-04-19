r"""
Tamari blossoming trees

This module implements the blossoming trees in bijection with Tamari intervals.

These blossoming trees are trees with half-edges bicolored following some local
rules, each node has two buds, and each edge has two legs. Buds are matched
with legs in a planar way, leaving only two dangling buds. The coloring can be
replaced by marking one of the dangling buds. The blossoming tree is
represented as a plane tree, in which a bud is exactly a leaf. We take the
convention that the root bud is a dangling bud with red half-edges next to it
in the counter-clockwise order.

REFERENCES:

[FFN2025]_

AUTHORS:

- Wenjie Fang (2026): initial implementation
"""

# ****************************************************************************
#       Copyright (C) 2026 Wenjie Fang <fwjmath@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from collections.abc import Iterator
from typing import Self
from math import acos, cos, sin, pi

from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.binary_tree import from_tamari_sorting_tuple, BinaryTree
from sage.combinat.dyck_word import DyckWords
from sage.combinat.integer_lists.invlex import IntegerListsLex
from sage.combinat.interval_posets import (TamariIntervalPoset,
                                           TamariIntervalPosets)
from sage.combinat.ordered_tree import OrderedTree, LabelledOrderedTree
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.latex import latex
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.prandom import uniform, randrange
from sage.plot.arc import arc
from sage.plot.arrow import arrow2d
from sage.plot.bezier_path import bezier_path
from sage.plot.circle import circle
from sage.plot.graphics import Graphics
from sage.plot.line import line
from sage.plot.polygon import polygon2d
from sage.rings.integer import Integer
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.sets.positive_integers import PositiveIntegers
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation


class TamariBlossomingTree(Element, UniqueRepresentation,
                           metaclass=InheritComparisonClasscallMetaclass):
    r"""
    The class of bicolored blossoming trees, which are in bijection with
    intervals in the Tamari lattice.

    A (bicolored) blossoming tree is a plane unrooted tree satisfying:

    - Each edge consists of two half-edges, one colored red and the other blue.
    - Each node has two extra unpaired uncolored half-edges called **buds**,
      which separate cyclically its other adjacent half-edges into two
      monochromatic parts. In other words, for each node, when reading color of
      adjacent half-edges in the cyclic order, colors change exactly twice, and
      the buds are placed at such position of changes.

    The size of a bicolored blossoming tree is the number of edges (not
    counting buds). They are in bijection with intervals in the Tamari lattice
    formed by binary trees of the same size (the number of internal nodes).

    As a convention, although a blossoming tree is unrooted, we represent it by
    a rooted plane tree with the root as the node with the dangling bud with red
    half-edges next to it in the counter-clockwise order. As a consequence, all
    internal nodes of the rooted plane tree have two buds, except for the root
    which only has one, separating blue half-edges on the left and red ones on
    the right.

    For usage, it is the best to use conversion functions provided by this
    class, instead of its constructor, which has less flexibility.

    .. NOTE::

        The metaclass
        :class:`~sage.misc.inherit_comparison.InheritComparisonClasscallMetaclass`
        comes from a conjoint of the metaclass
        :class:`~sage.misc.inherit_comparison.InheritComparisonMetaclass`
        imposed by the base class
        :class:`~sage.structure.unique_representation.UniqueRepresentation`
        and the metaclass
        :class:`~sage.misc.classcall_metaclass.ClasscallMetaclass`
        imposed by :class:`~sage.structure.element.Element`.

    EXAMPLES::

        sage: T = [[], [[[], []], [], []]]
        sage: TB = TamariBlossomingTree.from_plane_tree(T)
        sage: TB
        Tamari blossoming tree [[], [[], []], [[], []]] of size 2
        sage: B1, B2 = TB.to_tamari()
        sage: B1, B2
        ([[., .], .], [[., .], .])
        sage: TB is TamariBlossomingTree.from_tamari(B1, B2)
        True
        sage: TB.is_modern()
        True
        sage: TB.is_synchronized()
        True
    """

    @staticmethod
    def __checkbuds(tree):
        r"""
        Internal function. Check recursively whether every node has two buds. We
        do not suppose ``tree`` to be of type
        :class:`~sage.combinat.ordered_tree.OrderedTree`.

        TESTS::

            sage: TamariBlossomingTree(OrderedTree([]))
            Traceback (most recent call last):
            ...
            ValueError: not a blossoming tree, bud count incorrect
            sage: TamariBlossomingTree(OrderedTree([[], [[]]]))
            Traceback (most recent call last):
            ...
            ValueError: not a blossoming tree, bud count incorrect
            sage: TamariBlossomingTree(OrderedTree([[[], [], [[], []]], []]))
            Traceback (most recent call last):
            ...
            ValueError: not a blossoming tree, bad matching
            sage: T = OrderedTree([[], [[], [], [[], []]]])
            sage: TamariBlossomingTree(T)._tree
            1[2[], 3[4[], 5[], 6[7[], 8[]]]]
        """
        if not tree:
            return
        if len([x for x in tree if not x]) != 2:
            raise ValueError('not a blossoming tree, bud count incorrect')
        for st in tree:
            TamariBlossomingTree.__checkbuds(st)

    def __get_meandric_order(self) -> tuple[list[int], list[int]]:
        r"""
        Internal function. Compute the order of nodes and edges in the meandric
        representation. In the intermediate steps, we have:

        - Buds are represented by the label of its node and 1, 2 as the order
          of buds of the same node.
        - Legs are represented by the label of both its ends, and 1, 2 as the
          order of legs of the same node.
        """
        def aux(tree, budleg, budcnt):
            """
            This auxiliary function computes recursively a list of buds and legs
            with the representations:

            - A bud is represented by ``((r,), k)`` where ``r`` is the label of
              its node, and ``k`` its order in prefix order.
            - A leg is represented by ``((r1, r2), i)``, where ``r1`` and ``r2``
              are the labels of the nodes adjacent to the edge of the leg, and
              ``i`` presents the side it is on.
            """
            rlabel = tree.label()
            for st in tree:
                if not st:  # bud
                    if rlabel not in budcnt:
                        budcnt[rlabel] = 0
                    budcnt[rlabel] += 1
                    budleg.append(((rlabel,), budcnt[rlabel]))
                else:  # edge, two legs
                    budleg.append(((rlabel, st.label()), 1))
                    aux(st, budleg, budcnt)
                    budleg.append(((rlabel, st.label()), 2))

        # first compute the list of buds and legs
        budleg = [((self._tree.label(),), 1)]
        budcnt = {self._tree.label(): 1}
        aux(self._tree, budleg, budcnt)

        # then match them up
        matching, stack = [], []
        for halfedge in budleg:
            if len(halfedge[0]) == 1:  # bud
                stack.append(halfedge)
            else:  # leg
                matching.append((stack.pop(), halfedge))
        if len(stack) != 2:  # should never happen as budcount already checked
            raise ValueError('error during matching: incorrect matching')

        # get it in a dictionary
        matchdict = {}
        for bud, leg in matching:
            matchdict[bud] = leg
            matchdict[leg] = bud

        # go through the dictionary to get the path, thus the order
        curnode, curord = (self._tree.label(),), 1
        morder = [curnode]
        while True:
            curord = 3 - curord  # possible values are 1 and 2
            if (curnode, curord) not in matchdict:
                break
            curedge, curord = matchdict[curnode, curord]
            morder.append(curedge)
            curord = 3 - curord
            curnode, curord = matchdict[curedge, curord]
            morder.append(curnode)

        # A last check, but should never happen
        if len(morder) != self._size * 2 + 1:
            raise ValueError('error during matching: no Hamiltonian path')

        # compute both node order and edge order
        norder = [morder[i][0] for i in range(0, len(morder), 2)]
        eorder = [morder[i] for i in range(1, len(morder), 2)]

        # Finally, we return the order
        return norder, eorder

    @staticmethod
    def __classcall_private__(cls, tree):
        """
        Internal function serving two purposes:

        - Capture any instanciation of such that the constructed instance is
          always an element of
          :class:`~sage.combinat.tamari_blossoming_tree.TamariBlossomingTrees_all`
          so that the relation Parent/Element is unambiguous. We thus need to
          construct the object using the ``element_class`` method of
          :class:`~sage.combinat.tamari_blossoming_tree.TamariBlossomingTrees_all`.
          Note that all instances of
          :class:`~sage.combinat.tamari_blossoming_tree.TamariBlossomingTree`
          will have
          :class:`~sage.combinat.tamari_blossoming_tree.TamariBlossomingTrees_all`
          as parent.
        - Uniformize the input for cached representation. As we accept nested
          lists representing trees, which is not hashable, we need to convert
          it to :class:`~sage.combinat.ordered_tree.OrderedTree`, which is
          hashable, so that we may use a cached representation.

        TESTS::

            sage: T = OrderedTree([[], [[], [], [[], []]]])
            sage: B = TamariBlossomingTree(T)
            sage: B.parent()
            Tamari blossoming trees
            sage: B.parent() is TamariBlossomingTrees()
            True
            sage: B2 = TamariBlossomingTrees()(T)
            sage: B2 is B
            True
            sage: B3 = TamariBlossomingTrees(2)(T)
            sage: B3 is B
            True
            sage: B.parent() is B2.parent()
            True
            sage: B.parent() is B3.parent()
            True
        """
        A = TamariBlossomingTrees_all()
        if not isinstance(tree, OrderedTree):
            tree = OrderedTree(tree)
        return A.element_class(A, tree)

    def __init__(self, parent, tree: OrderedTree) -> None:
        r"""
        Initialize a Tamari blossoming tree with a plane tree.

        We consider the root as a bud, so every internal node has two leaves,
        except the root which has only one. We also check that the root is
        really a dangling node (i.e., adjacent to a dangling bud), and raise an
        error otherwise.

        INPUT:

        - ``tree`` -- a plane tree of the type
          :class:`~sage.combinat.ordered_tree.OrderedTree`.

        EXAMPLES::

            sage: T = OrderedTree([[], [[], [], [[], []]]])
            sage: TamariBlossomingTree(T)._tree
            1[2[], 3[4[], 5[], 6[7[], 8[]]]]

        TESTS::

            sage: TamariBlossomingTree(OrderedTree([]))
            Traceback (most recent call last):
            ...
            ValueError: not a blossoming tree, bud count incorrect
            sage: TamariBlossomingTree(OrderedTree([[], [[]]]))
            Traceback (most recent call last):
            ...
            ValueError: not a blossoming tree, bud count incorrect
            sage: TamariBlossomingTree(OrderedTree([[[], [], [[], []]], []]))
            Traceback (most recent call last):
            ...
            ValueError: not a blossoming tree, bad matching
        """
        def matching_word(tree):
            """
            Internal function. Return an iterator of the matching word with
            buds as 1 and legs as 0. We do not suppose ``tree`` to be of type
            :class:`~sage.combinat.ordered_tree.OrderedTree`. We do not count
            the root here.
            """
            accu = []
            for t in tree:
                if not t:  # a bud
                    yield 1
                else:  # not a bud, but an edge to the next subtree
                    yield -1
                    yield from matching_word(t)
                    yield -1

        # case of rebuilding the element
        if isinstance(tree, TamariBlossomingTree):
            return
        # check root leaves
        if len([x for x in tree if not x]) != 1:
            raise ValueError('not a blossoming tree, bud count incorrect')

        # check for all nodes
        for st in tree:
            TamariBlossomingTree.__checkbuds(st)

        # check for matching (whether the root is a dangling bud)
        ht = 0
        for e in matching_word(tree):
            ht += e
            if ht < 0:
                raise ValueError('not a blossoming tree, bad matching')

        # all tests passed, construct the object
        # the tree, with the canonical labeling
        self._tree = OrderedTree(tree).canonical_labelling()
        # size is given by the number of edges in the tree (excluding buds)
        self._size = (self._tree.node_number() - 1) // 3
        # the meandric order of nodes
        self._node_order, self._edge_order = self.__get_meandric_order()
        # element/parent class
        Element.__init__(self, parent)

    def _repr_(self) -> str:
        r"""
        Return a string representing the blossoming tree and its size.

        TESTS::

            sage: T = OrderedTree([[], [[], [], [[], []]]])
            sage: TamariBlossomingTree(T)
            Tamari blossoming tree [[], [[], [], [[], []]]] of size 2
        """
        s = f'Tamari blossoming tree {OrderedTree(self._tree)}'
        s += f' of size {self._size}'
        return s

    def _ascii_art_(self):
        r"""
        Return the ascii art of the blossoming tree, using that of OrderedTree.

        TESTS::

            sage: T = OrderedTree([[], [[], [], [[], []]]])
            sage: ascii_art(TamariBlossomingTree(T))
              __o___
             /     /
            o   __o___
               / /   /
              o o   o_
                   / /
                  o o
        """
        return OrderedTree(self._tree)._ascii_art_()

    def _unicode_art_(self):
        r"""
        Return the unicode art of the blossoming tree, using that of
        OrderedTree.

        TESTS::

            sage: T = OrderedTree([[], [[], [], [[], []]]])
            sage: unicode_art(TamariBlossomingTree(T))
            ╭──o─╮
            │    │
            o ╭─┬o─╮
              │ │  │
              o o ╭o╮
                  │ │
                  o o
        """
        return OrderedTree(self._tree)._unicode_art_()

    def __hash__(self) -> int:
        r"""
        Return the hash value of the blossoming tree.

        TESTS::

            sage: TBT5 = TamariBlossomingTrees(5)
            sage: len(set(hash(T) for T in TBT5)) == TBT5.cardinality()
            True
        """
        return self._tree.__hash__()

    def __eq__(self, other) -> bool:
        r"""
        Equality test, which only needs to compare the underlying trees.

        .. NOTE::

            This delegates the comparison to
            :class:`~sage.combinat.ordered_tree.LabelledOrderedTree`

        TESTS::

            sage: T1 = [[], [[[], []], [], []]]
            sage: TB1 = TamariBlossomingTree.from_plane_tree(T1)
            sage: T2 = OrderedTree([[], [[], []], [[], []]])
            sage: TB2 = TamariBlossomingTree(T2)
            sage: TB1 == TB2
            True
        """
        if not isinstance(other, TamariBlossomingTree):
            return False
        return self._tree == other._tree

    def __ne__(self, other) -> bool:
        r"""
        Test for inequality, uses ``__eq__``.

        TESTS::

            sage: T = OrderedTree([[], [[], []], [[], []]])
            sage: TB1 = TamariBlossomingTree.from_plane_tree(T)
            sage: TB2 = TamariBlossomingTree(T)
            sage: TB1 != TB2  # should be different, even from the same tree
            True
        """
        return not self == other

    def size(self) -> Integer:
        r"""
        Return the size of the Tamari blossoming tree.

        OUTPUT:

        The size of the current blossoming tree, which is the number of edges
        excluding buds. This convention agrees with that of
        :class:`~sage.combinat.interval_posets.TamariIntervalPoset`.

        EXAMPLES::

            sage: T = OrderedTree([[], [[], [], [[], []]]])
            sage: T = TamariBlossomingTree(T)
            sage: T.size()
            2
            sage: T.size() == T.to_TIP().size()
            True
        """
        return Integer(self._size)

    def to_plane_tree(self) -> OrderedTree:
        r"""
        Return the blossoming tree as an
        :class:`~sage.combinat.ordered_tree.OrderedTree`.

        The buds simply become leaves.

        EXAMPLES::

            sage: T = OrderedTree([[], [[], [], [[], []]]])
            sage: T == TamariBlossomingTree(T).to_plane_tree()
            True
        """
        return OrderedTree(self._tree)

    def to_tamari(self) -> tuple[BinaryTree, BinaryTree]:
        r"""
        Return the Tamari interval in bijection with ``self``, under the form
        of a pair of binary trees.

        OUTPUT:

        A pair of binary trees comparable in the Tamari lattice (thus a Tamari
        interval)

        EXAMPLES::

            sage: B0 = BinaryTree()
            sage: T0 = OrderedTree([[]])
            sage: TamariBlossomingTree(T0).to_tamari() == (B0, B0)
            True
            sage: T = OrderedTree([[], [[], [], [[], []]]])
            sage: TamariBlossomingTree(T).to_tamari()
            ([., [., .]], [., [., .]])
            sage: Tl = OrderedTree([[], [[], [], [[], []], [[[], []], [], []]],
            ....:                   [[], [], [[[], []], [], [[], []], []]]])
            sage: B3, B4 = TamariBlossomingTree(Tl).to_tamari()
            sage: B3
            [[[[., [., [., .]]], .], [[., .], .]], .]
            sage: B4
            [[., [[., [., .]], .]], [., [[., .], .]]]
        """
        def from_dual_bracket_vector(dvec):
            """
            This function converts dual bracket vectors to binary trees
            """
            if not dvec:
                return BinaryTree()
            ridx = len(dvec) - 1
            while ridx != dvec[ridx]:
                ridx -= 1
            ltree = from_dual_bracket_vector(dvec[:ridx])
            rtree = from_dual_bracket_vector(dvec[ridx + 1:])
            return BinaryTree([ltree, rtree])

        # get the orders of nodes and edges
        norder = self._node_order
        eorder = self._edge_order

        # get the bracket vector (lower) and the dual bracket vector (higher)
        bvec, dvec = [], []
        for i, eoi in enumerate(eorder):
            idx = sorted(norder.index(x) for x in eoi)
            bvec.append(idx[1] - 1 - i)
            dvec.append(i - idx[0])

        # get the trees
        ltree = from_tamari_sorting_tuple(bvec)
        rtree = from_dual_bracket_vector(dvec)
        return ltree, rtree

    @staticmethod
    def from_tamari(ltree, htree) -> Self:
        r"""
        Return the blossoming tree corresponding to the given Tamari interval,
        given as a pair of binary trees (not necessarily of type
        :class:`~sage.combinat.binary_tree.BinaryTree`).

        INPUT:

        - ``ltree`` -- the lower binary tree in the given Tamari interval.
        - ``htree`` -- the upper binary tree in the given Tamari interval.

        OUTPUT:

        A blossoming tree in bijection with the given Tamari interval

        EXAMPLES::

            sage: B0 = BinaryTree()
            sage: TB0 = TamariBlossomingTree(OrderedTree([[]]))
            sage: TamariBlossomingTree.from_tamari(B0, B0) == TB0
            True
            sage: Tl = OrderedTree([[], [[], [], [[], []], [[[], []], [], []]],
            ....:                   [[], [], [[[], []], [], [[], []], []]]])
            sage: TlB = TamariBlossomingTree(Tl)
            sage: B3 = [[[[None, [None, []]], None], [[], None]], None]
            sage: B4 = [[None, [[None, []], None]], [None, [[], None]]]
            sage: TamariBlossomingTree.from_tamari(B3, B4) == TlB
            True
            sage: TamariBlossomingTree.from_tamari(B4, B3)
            Traceback (most recent call last):
            ...
            ValueError: not a Tamari interval
        """
        def traversal(node, parent, cycord):
            # internal function, which go through the tree given by cycord
            # we provide parent to know where to cut
            if node == -1:  # a bud
                return []
            children = cycord[node]
            if parent in cycord[node]:
                pidx = children.index(parent)
                children = children[pidx + 1:] + children[:pidx]
            return [traversal(x, node, cycord) for x in children]

        # initialization and verification
        btl, bth = BinaryTree(ltree), BinaryTree(htree)
        if not btl.tamari_lequal(bth):
            raise ValueError('not a Tamari interval')

        # compute the bracket vector and the dual bracket vector
        bvec = btl.tamari_sorting_tuple()[0]
        dvec = bth.tamari_sorting_tuple(reverse=True)[0][::-1]
        n = len(bvec)

        # construct the edges between nodes, with their order
        upper = [[] for _ in range(n + 1)]  # neighbors by upper arcs
        lower = [[] for _ in range(n + 1)]  # neighbors by lower arcs
        for i in range(n):
            a, b = i - dvec[i], bvec[i] + i + 1
            upper[a].append(b)  # counter-clockwise
            lower[b].append(a)  # clockwise
        # edges in counterclockwise order (left to right for trees)
        # buds are represented by -1
        cycord = [[-1] + lower[i][::-1] + [-1] + upper[i] for i in range(n + 1)]
        # get rid of the first bud (ugly pop...)
        cycord[0].pop(0)
        # traversal for the plane tree
        # 2 saying that it is the root (so parent inexistent)
        ptree = traversal(0, -2, cycord)
        return TamariBlossomingTree(OrderedTree(ptree))

    def to_TIP(self) -> TamariIntervalPoset:
        r"""
        Return the corresponding Tamari interval-poset in bijection with the
        current instance of blossoming tree.

        OUTPUT:

        An object of type
        :class:`~sage.combinat.interval_posets.TamariIntervalPoset`
        representing the corresponding Tamari interval poset

        EXAMPLES::

            sage: Tl = OrderedTree([[], [[], [], [[], []], [[[], []], [], []]],
            ....:                   [[], [], [[[], []], [], [[], []], []]]])
            sage: TlB = TamariBlossomingTree(Tl)
            sage: tip = TlB.to_TIP()
            sage: B3 = [[[[None, [None, []]], None], [[], None]], None]
            sage: B4 = [[None, [[None, []], None]], [None, [[], None]]]
            sage: tip.lower_binary_tree() == BinaryTree(B3)
            True
            sage: tip.upper_binary_tree() == BinaryTree(B4)
            True
        """
        t1, t2 = self.to_tamari()
        return TamariIntervalPosets.from_binary_trees(t1, t2)

    @staticmethod
    def from_TIP(tip) -> Self:
        r"""
        Construct the blossoming tree in bijection with a given Tamari
        interval-poset.

        INPUT:

        - ``tip`` -- a
          :class:`~sage.combinat.interval_posets.TamariIntervalPoset` object
          representing a Tamari interval

        EXAMPLES::

            sage: B3 = [[[[None, [None, []]], None], [[], None]], None]
            sage: B4 = [[None, [[None, []], None]], [None, [[], None]]]
            sage: B3, B4 = BinaryTree(B3), BinaryTree(B4)
            sage: tip = TamariIntervalPosets.from_binary_trees(B3, B4)
            sage: Tl = OrderedTree([[], [[], [], [[], []], [[[], []], [], []]],
            ....:                   [[], [], [[[], []], [], [[], []], []]]])
            sage: TlB = TamariBlossomingTree(Tl)
            sage: TlB == TamariBlossomingTree.from_TIP(tip)
            True
            sage: B0 = BinaryTree()
            sage: tip2 = TamariIntervalPosets.from_binary_trees(B0, B0)
            sage: TlB == TamariBlossomingTree.from_TIP(tip2)
            False
        """
        t1, t2 = tip.lower_binary_tree(), tip.upper_binary_tree()
        return TamariBlossomingTree.from_tamari(t1, t2)

    @staticmethod
    def binary_tree_plot(btree) -> Graphics:
        r"""
        Utility function for plotting binary trees in the "upside-down
        Chapoton" way, i.e., leaves are drawn on a horizontal line, but the
        root still on top.

        INPUT:

        - ``btree`` -- a binary tree, not necessarily of type BinaryTree

        OUTPUT:

        A plot of ``btree``, with leaves on a horizontal line, and edges all of
        slope `1` or `-1`. Labels are ignored.

        EXAMPLES::

            sage: B3 = [[[[None, [None, []]], None], [[], None]], None]
            sage: g = TamariBlossomingTree.binary_tree_plot(B3)

        .. PLOT::

            B3 = [[[[None, [None, []]], None], [[], None]], None]
            g = TamariBlossomingTree.binary_tree_plot(B3)
            sphinx_plot(g)
        """
        # auxiliary function to compute coordinates of internal nodes
        def aux(t, a, b, points):
            if t.is_empty():
                return
            # current point
            points.append(((a + b) / 2, (b - a) / 2))
            # recursive calls
            lt, rt = tuple(t)
            lsize = lt.node_number()
            aux(lt, a, a + lsize, points)
            aux(rt, a + lsize + 1, b, points)
            return

        # initialization
        bt = BinaryTree(btree)
        n = bt.node_number()
        G = Graphics()
        G.set_aspect_ratio(1)

        # get node positions
        points = []
        aux(bt, 0, n, points)

        # plot, first leaves, then nodes and edges
        for i in range(n + 1):
            G += circle((i, 0), 0.05, fill=True, zorder=2)
        for x, y in points:
            G += circle((x, y), 0.2, fill=True, facecolor='white', zorder=2)
            G += line([(x, y), (x + y, 0)], zorder=1, thickness=1)
            G += line([(x, y), (x - y, 0)], zorder=1, thickness=1)
        G.axes(show=False)
        return G

    def tamari_dual(self) -> Self:
        r"""
        Return the current blossoming tree with colors exchanged.

        This is equivalent to taking the dual in the Tamari lattice for the
        corresponding interval. **Not to be confused** with the mirror
        symmetric of blossoming trees, which is implemented in
        :meth:`reflection`.

        .. NOTE::

            We use a feature of :meth:`from_plane_tree` instead of using the
            usual
            :meth:`~sage.combinat.binary_tree.BinaryTree.left_right_symmetry`
            in :class:`~sage.combinat.binary_tree.BinaryTree` to avoid
            exceeding the limit on the number of recursions when tree size is
            large.

        EXAMPLES::

            sage: B3 = [[[[None, [None, []]], None], [[], None]], None]
            sage: B4 = [[None, [[None, []], None]], [None, [[], None]]]
            sage: B3, B4 = BinaryTree(B3), BinaryTree(B4)
            sage: TB = TamariBlossomingTree.from_tamari(B3, B4)
            sage: B4d, B3d = TB.tamari_dual().to_tamari()
            sage: B4d == B4.left_right_symmetry()
            True
            sage: B3d == B3.left_right_symmetry()
            True
        """
        # use the fact that from_plane_tree picks the bud with the opposite
        # color of the root, so exchanges the color (thus duality)
        return TamariBlossomingTree.from_plane_tree(self._tree)

    def plot_meandric(self, semicircle=True, arrow=True) -> Graphics:
        r"""
        Plot the meandric tree of ``self``.

        The meandric tree is the blossoming tree drawn in a way such that the
        drawing is planar, with buds all on the same horizontal line, which
        separates the red half-edges above and the blue ones below. For more
        details, see [FFN2025]_.

        INPUT:

        - ``semicircle`` -- optional, sets whether the arcs are drawn as
          semicircles or a Bezier arc. Default: ``True``, which draws
          semicircles.
        - ``arrow`` -- optional, sets whether draw horizontal arrows tips.
          Default: ``True``.

        EXAMPLES::

            sage: Tl = OrderedTree([[], [[], [], [[], []], [[[], []], [], []]],
            ....:                   [[], [], [[[], []], [], [[], []], []]]])
            sage: g = TamariBlossomingTree(Tl).plot_meandric()

        .. PLOT::

            Tl = OrderedTree([[], [[], [], [[], []], [[[], []], [], []]],
                              [[], [], [[[], []], [], [[], []], []]]])
            g = TamariBlossomingTree(Tl).plot_meandric()
            sphinx_plot(g)
        """
        def sqnode(x, y):
            """
            Draw a white square node (middle of segments) at position (x, y).
            """
            diam = 0.1
            return polygon2d([[x - diam, y - diam], [x + diam, y - diam],
                              [x + diam, y + diam], [x - diam, y + diam]],
                             edgecolor='black', rgbcolor='white', zorder=2)

        def cirnode(x, y):
            """
            Draw a black circle node at position (x, y).
            """
            return circle([x, y], 0.15, fill=True, edgecolor='black',
                          facecolor='black', zorder=2)

        def semicir(x1, x2, isupper):
            """
            Draw a semi-circle from (x1, 0) to (x2, 0), with ``isupper``
            indicating whether it is in the upper or lower plane.
            """
            sec = (0, pi) if isupper else (pi, 2 * pi)
            color = 'blue' if isupper else 'red'
            return arc([(x1 + x2) / 2, 0], (x2 - x1) / 2, sector=sec, zorder=1,
                       rgbcolor=color)

        def bezierarc(x1, x2, isupper):
            """
            Draw a Bezier arc from (x1, 0) to (x2, 0), with ``isupper``
            indicating whether it is in the upper or lower plane.
            """
            cp1 = [x1 * 0.7 + x2 * 0.3, (x2 - x1) * 0.6]
            cp2 = [x1 * 0.3 + x2 * 0.7, (x2 - x1) * 0.6]
            if not isupper:
                cp1[1], cp2[1] = -cp1[1], -cp2[1]
            cp1 = tuple(cp1)
            cp2 = tuple(cp2)
            color = 'blue' if isupper else 'red'
            return bezier_path([[(x1, 0), cp1, cp2, (x2, 0)]], zorder=1,
                               rgbcolor=color)

        # initialization
        G = Graphics()
        G.set_aspect_ratio(1)
        arcfct = semicir if semicircle else bezierarc

        # draw the vertices, black circle for tree node, white squares for edges
        n = self._size
        for i in range(2 * n + 1):  # tree nodes
            if i % 2 == 0:
                G += cirnode(i, 0)
                if arrow:
                    G += arrow2d((i, 0), (i + 0.6, 0), rgbcolor='black',
                                 width=1, arrowsize=2)
                    G += arrow2d((i, 0), (i - 0.6, 0), rgbcolor='black',
                                 width=1, arrowsize=2)
                else:
                    G += line([(i, 0), (i + 0.6, 0)], rgbcolor='black')
                    G += line([(i, 0), (i - 0.6, 0)], rgbcolor='black')
            else:
                G += sqnode(i, 0)

        # draw the arcs (tree edges), depending on options
        norder, eorder = self._node_order, self._edge_order
        for i in range(n):
            nidx1, nidx2 = eorder[i]
            k, m = sorted((norder.index(nidx1), norder.index(nidx2)))
            G += arcfct(k * 2, i * 2 + 1, True)   # upper arc
            G += arcfct(i * 2 + 1, m * 2, False)  # lower arc
        G.axes(show=False)
        return G

    def _latex_(self) -> str:
        r"""
        Return latex code for the meandric drawing of the current blossoming
        tree.

        EXAMPLES::

            sage: T = OrderedTree([[], [[], [], [[], []]]])
            sage: T = TamariBlossomingTree(T)
            sage: type(latex(T))
            <class 'sage.misc.latex.LatexExpr'>
        """
        n = self._size
        tikz = []
        latex.add_package_to_preamble_if_available('tikz')
        # header
        tikz.append('\\begin{tikzpicture}\n')
        # zorder=0
        # arrows
        tikz.append(f'\\foreach \\x in {{0, 2, ..., {2 * n}}} '
                    '\\draw[-latex, thick] (\\x, 0) -- ++(-0.6, 0);\n')
        tikz.append(f'\\foreach \\x in {{0, 2, ..., {2 * n}}} '
                    '\\draw[-latex, thick] (\\x, 0) -- ++(0.6, 0);\n')
        # zorder=1
        # tree edges
        norder, eorder = self._node_order, self._edge_order
        for i in range(n):
            nidx1, nidx2 = eorder[i]
            k, m = sorted((norder.index(nidx1), norder.index(nidx2)))
            # upper arc
            tikz.append(f'\\draw[blue, very thick] ({k * 2}, 0) arc '
                        f'(180:0:{i - k + 0.5});\n')
            # lower arc
            tikz.append(f'\\draw[red, very thick] ({m * 2}, 0) arc '
                        f'(0:-180:{m - i - 0.5});\n')
        # zorder=2
        # trees nodes, which are circles
        tikz.append(f'\\foreach \\x in {{0,  2, ..., {2 * n}}} \\filldraw[black]'
                    ' (\\x, 0) circle (0.15);\n')
        # square nodes for edges, which are squares
        tikz.append(f'\\foreach \\x in {{1, 3, ..., {2 * n - 1}}}'
                    '\\node[draw=black, fill=white, '
                    'very thick, minimum size=0.2] '
                    '(square) at (\\x, 0) {};\n')
        # ending
        tikz.append('\\end{tikzpicture}\n')
        return ''.join(tikz)

    @staticmethod
    def __find_dangling_bud(tree: LabelledOrderedTree) -> list[int]:
        r"""
        Internal function. Return the dangling buds of ``tree``, in the order
        of counterclockwise order starting from the root. In ``tree`` we assume
        that there is a bud for the root (that is not in ``tree``), labeled 0
        (which should not be in canonical labeling of ``tree``). We also assume
        that the tree is valid. Note that the root is not necessarily a dangling
        bud.

        TESTS::

            sage: T1 = OrderedTree([[], [[], [], [[], []]]])
            sage: T1 = TamariBlossomingTree(T1)
            sage: T2 = OrderedTree([[[], [], [[], []]], []])
            sage: T2 = TamariBlossomingTree.from_plane_tree(T2)
            sage: T1 == T2
            True
        """
        def aux(t, buds, dyck):
            for st in t:
                if not st:  # bud
                    buds.append(st.label())
                    dyck.append(1)
                else:  # edge and subtree
                    dyck.append(-1)
                    aux(st, buds, dyck)
                    dyck.append(-1)

        buds = [0]  # the root bud
        dyck = [1]  # the root bud
        aux(tree, buds, dyck)

        # find the latest lowest point
        lowest, lidx, height = 0, 0, 0
        for i, step in enumerate(dyck):
            height += step
            if height <= lowest:
                lowest, lidx = height, i + 1
        # find the corresponding bud
        budcnt = dyck[:lidx].count(1)
        # splice to get a walk above 0, go for the latest point with height 1
        ndyck = dyck[lidx:] + dyck[:lidx]
        oneidx, height = 0, 0
        for i, step in enumerate(ndyck):
            height += step
            if height == 1 and step < 0:
                oneidx = i
        oneidx += lidx
        if oneidx > len(ndyck):
            oneidx -= len(ndyck)
        budcnt2 = dyck[:oneidx].count(1)
        return [buds[budcnt], buds[budcnt2]]

    @staticmethod
    def __get_cycle_order(t: LabelledOrderedTree) -> list[int]:
        r"""
        Internal function. Return the cyclic order (in the sens of maps) of the
        given tree.

        More precisely, this function returns a dictionary with
        node labels as keys and a (cyclic) list of its neighbors in
        counterclockwise order. The root bud is labeled 0, under the assumption
        that 0 is not present in the canonical labeling.
        """
        def aux(tree, parent, cycord):
            cycord[tree.label()] = [parent] + [st.label() for st in tree]
            for st in tree:
                aux(st, tree.label(), cycord)

        cycord = {0: [t.label()]}
        aux(t, 0, cycord)
        return cycord

    @staticmethod
    def from_plane_tree(tree) -> Self:
        r"""
        Return the blossoming tree corresponding to the given tree.

        We suppose that the root of the tree is a bud. Comparing to the
        constructor, we do not fail when the buds are not matching, but tries
        to find the correct bud.

        We assume that the root, which is a bud, has red half-edges next to it
        in counter-clockwise order (so the left one). We then find the unpaired
        bud with the opposite property, to simplify the reflection operation.

        .. NOTE::

            This function is a thin wrapper of the internal function
            ``_from_plane_tree``, which provides two extra functionalities
            by which end users should not be concerned.

        INPUT:

        - ``tree`` -- a plane tree with two buds on each node (one for the
          root).

        OUTPUT:

        An object of type
        :class:`~sage.combinat.tamari_blossoming_tree.TamariBlossomingTree`
        representing the blossoming tree given by ``tree``.

        EXAMPLES::

            sage: pt = [[], [[[], []], [], []]]
            sage: B1, B2 = TamariBlossomingTree.from_plane_tree(pt).to_tamari()
            sage: B1 == BinaryTree([[], None])
            True
            sage: B1 == B2
            True

        TESTS::

            sage: TamariBlossomingTree.from_plane_tree([[], []])
            Traceback (most recent call last):
            ...
            ValueError: not a blossoming tree, bud count incorrect
            sage: TamariBlossomingTree(OrderedTree([[], [[[], []], [], []]]))
            Traceback (most recent call last):
            ...
            ValueError: not a blossoming tree, bad matching
        """
        return TamariBlossomingTree._from_plane_tree(tree)

    @staticmethod
    def _from_plane_tree(tree, skip_check=False, random_bud=False) -> Self:
        r"""
        Return the blossoming tree corresponding to the given tree.

        We suppose that the root of the tree is a bud. Comparing to the
        constructor, we do not fail when the buds are not matching, but tries
        to find the correct bud.

        We assume that the root, which is a bud, has red half-edges next to it
        in counter-clockwise order (so the left one). We then find the unpaired
        bud with the opposite property, to simplify the reflection operation.

        .. NOTE::

            This function comes with two extra options that are not meant to be
            used by end users, with which we may skip the check for bud
            conditions (each node has two buds), and we may pick one of the
            dangling buds at random as the root to simplify the random
            generation.

        INPUT:

        - ``tree`` -- a plane tree with two buds on each node (one for the
          root).
        - ``skip_check`` -- skip checking for bud conditions. Default:
          ``False``.
        - ``random_bud`` -- choose a bud at random, instead of the first.
          Default: ``False``.

        OUTPUT:

        An object of type
        :class:`~sage.combinat.tamari_blossoming_tree.TamariBlossomingTree`
        representing the blossoming tree given by ``tree``.

        TESTS::

            sage: TamariBlossomingTree._from_plane_tree([[], []])
            Traceback (most recent call last):
            ...
            ValueError: not a blossoming tree, bud count incorrect
            sage: TamariBlossomingTree._from_plane_tree([[], []],
            ....:                                       skip_check=True)
            Traceback (most recent call last):
            ...
            ValueError: invalid blossoming tree: bud colors
            sage: pt = [[], [[[], []], [], []]]
            sage: res = {TamariBlossomingTree._from_plane_tree(pt,
            ....:                                              random_bud=True)
            ....:        for _ in range(100)}
            sage: len(res) <= 2
            True
            sage: while len(res) < 2:
            ....:     T = TamariBlossomingTree._from_plane_tree(pt,
            ....:                                               random_bud=True)
            ....:     res.add(T)
            sage: len(res)
            2
        """
        def traverse(node, parent, cycord):
            """
            Internal function, construct a plane tree out of the cycle order.
            The parameter ``parent`` is for knowing where to cut.
            """
            pidx = cycord[node].index(parent)
            stnodes = cycord[node][pidx + 1:] + cycord[node][:pidx]
            return [traverse(stn, node, cycord) for stn in stnodes]

        # check buds
        if not skip_check:
            if len([x for x in tree if not x]) != 1:
                raise ValueError('not a blossoming tree, bud count incorrect')
            for st in tree:
                TamariBlossomingTree.__checkbuds(st)

        tree = OrderedTree(tree).canonical_labelling()
        dangling = TamariBlossomingTree.__find_dangling_bud(tree)

        # compute bud color
        cycord = TamariBlossomingTree.__get_cycle_order(tree)
        dcolor = [0, 0]
        for i in range(2):
            bud = dangling[i]
            if bud == 0:
                dcolor[i] = 0
                continue
            color = 0
            curpos = bud
            while curpos != 0:  # going up to the root
                prevpos = curpos
                curpos = cycord[curpos][0]
                pidx = cycord[curpos].index(prevpos)
                for sibling in cycord[curpos][pidx + 1:]:
                    if len(cycord[sibling]) == 1:  # a bud
                        color = 1 - color
                color = 1 - color  # going to the opposite half-edge
            dcolor[i] = 1 - color  # accounting for the root bud

        # select against colors
        if sum(dcolor) != 1:
            raise ValueError('invalid blossoming tree: bud colors')
        didx = dcolor.index(1)  # select the opposite color
        if random_bud:
            didx = randrange(2)
        dangling = dangling[didx]
        rroot = cycord[dangling][0]  # the only neighbor of a bud is the root
        rtree = traverse(rroot, dangling, cycord)
        return TamariBlossomingTree(OrderedTree(rtree))

    def reflection(self) -> Self:
        r"""
        Return the blossoming tree that is the mirror image of the current
        blossoming tree.

        Note that the dangling buds change in general, so the
        root dangling bud will be the one with the correct color.

        OUTPUT:

        A blossoming tree that is the mirror image of the current one. Not to be
        confused with the Tamari dual computed by :meth:`tamari_dual`.

        EXAMPLES::

            sage: pt = [[], [[[], []], [], []]]
            sage: TB = TamariBlossomingTree.from_plane_tree(pt)
            sage: TB == TB.reflection().reflection()
            True
            sage: TB.reflection().to_plane_tree()
            [[], [[], []], [[], []]]
        """
        tree = self.to_plane_tree().left_right_symmetry()
        return TamariBlossomingTree._from_plane_tree(tree, skip_check=True)

    def plot_blossoming(self, aspect=1.0, layout='tree') -> Graphics:
        r"""
        Plot the blossoming tree of ``self``, using the plot of
        :class:`~sage.combinat.ordered_tree.OrderedTree`, but with buds better
        spaced.

        INPUT:

        - ``aspect`` -- ratio of aspect, default to ``1.0``
        - ``layout`` -- the algorithm for layout, with three possible options:
            * ``'tree'``: uses ``layout_tree`` of undirected graph in sage
            * ``'planar'``: uses ``layout_planar`` of undirected graph in sage

        OUTPUT:

        A plot of the current blossoming tree.

        EXAMPLES::

            sage: Tl = OrderedTree([[], [[], [], [[], []], [[[], []], [], []]],
            ....:                   [[], [], [[[], []], [], [[], []], []]]])
            sage: g = TamariBlossomingTree(Tl).plot_blossoming()

        .. PLOT::

            Tl = OrderedTree([[], [[], [], [[], []], [[[], []], [], []]],
                              [[], [], [[[], []], [], [[], []], []]]])
            g = TamariBlossomingTree(Tl).plot_blossoming()
            sphinx_plot(g)
        """
        def euclid_dist(p1, p2):
            return sum([(p1[i] - p2[i]) ** 2 for i in range(2)]) ** 0.5

        def rad_dir(p1, p2):
            plen = euclid_dist(p1, p2)
            p = [(p2[i] - p1[i]) / plen for i in range(2)]
            r = acos(p[0])
            if p[1] < 0:
                r = -r
            return r

        def shift_point(p, rad, dist):
            px = p[0] + cos(rad) * dist
            py = p[1] + sin(rad) * dist
            return (px, py)

        def plot_bud(origp, rad, m, bud, dbuds):
            p2 = shift_point(pos[rn], rad, m)
            w = 1
            color = 'red'
            if bud == dbuds[0]:
                w = 3
                color = 'darkgreen'
            elif bud == dbuds[1]:
                w = 2
                color = 'darkgreen'
            return line([pos[rn], p2], rgbcolor=color, thickness=w, zorder=1)

        # construct the graph and the embedding
        t = LabelledOrderedTree([self._tree], label=0)
        cycord = TamariBlossomingTree.__get_cycle_order(self._tree)
        degs = [len(cycord[x]) for x in cycord]
        n = len(cycord)
        realnodes = [i for i in range(n) if degs[i] > 1]
        g = t.to_undirected_graph()
        # using clockwise direction, instead of counterclockwise for trees
        embed = {x: cycord[x][::-1] for x in cycord}
        g.set_embedding(embed)
        # use appropriate layout algorithm
        pos = None
        if layout == 'tree':
            pos = g.layout_tree()
        elif layout == 'planar':
            pos = g.layout_planar(on_embedding=embed, external_face=(0, 1))
        else:
            raise ValueError('invalid option for parameter "layout"')

        # Initialize the graphic
        G = Graphics()
        G.set_aspect_ratio(aspect)

        # Normalize node positions
        xcoords = [pos[i][0] for i in realnodes]
        minx, maxx = min(xcoords), max(xcoords)
        ycoords = [pos[i][1] for i in realnodes]
        miny, maxy = min(ycoords), max(ycoords)
        multfact = aspect / (maxy - miny) * (maxx - minx)
        for node in realnodes:
            e = pos[node]
            pos[node] = (e[0], (e[1] - miny) * multfact + miny)

        # compute min edge for scaling (excluding buds)
        minedge = None
        for n1 in realnodes:
            for n2 in [i for i in cycord[n1] if degs[i] > 1]:
                edlen = euclid_dist(pos[n1], pos[n2])
                if minedge is None or minedge > edlen:
                    minedge = edlen

        # draw the real nodes (not the buds)
        for node in realnodes:
            G += circle(pos[node], 0.02 * edlen, fill=True, zorder=2)

        # draw the real edges
        for n1 in realnodes:
            for n2 in [i for i in cycord[n1] if degs[i] > 1 and i > n1]:
                G += line([pos[n1], pos[n2]], zorder=1, thickness=1)

        # draw the buds
        dbuds = TamariBlossomingTree.__find_dangling_bud(self._tree)
        budlen = minedge * 0.3
        for rn in realnodes:
            ncnt = len(cycord[rn])
            budidx = [i for i in range(ncnt) if degs[cycord[rn][i]] == 1]
            if budidx[1] - budidx[0] in [1, ncnt - 1]:  # two consecutive buds
                rad1 = rad_dir(pos[rn], pos[cycord[rn][budidx[0] - 1]])
                eidx2 = budidx[1] + 1 - ncnt  # using negative index for cyclic
                rad2 = rad_dir(pos[rn], pos[cycord[rn][eidx2]])
                if rad2 <= rad1:
                    rad2 += pi * 2
                # trisection of angle
                rbuds = [rad1 + (rad2 - rad1) / 3 * (1 + i) for i in range(2)]
                for i in range(2):
                    G += plot_bud(pos[rn], rbuds[i], budlen,
                                  cycord[rn][budidx[i]], dbuds)
            else:  # two non-consecutive buds, we put each one in the middle
                for i in range(2):
                    rad1 = rad_dir(pos[rn], pos[cycord[rn][budidx[i] - 1]])
                    eidx2 = budidx[i] + 1 - ncnt
                    rad2 = rad_dir(pos[rn], pos[cycord[rn][eidx2]])
                    if rad2 <= rad1:
                        rad2 += pi * 2
                    rbud = (rad1 + rad2) / 2
                    G += plot_bud(pos[rn], rbud, budlen,
                                  cycord[rn][budidx[i]], dbuds)

        # output
        G.axes(show=False)
        return G

    @staticmethod
    def _binary_tree_arcs(btree: BinaryTree) -> list[tuple[int]]:
        """
        Internal function. Return the list of arcs in the smooth drawing of a
        given binary tree.

        INPUT:

        - ``btree`` -- a binary tree. Both nested lists and tree objects in
          Sagemath are accepted.

        OUTPUT:

        The list of arcs in the smooth drawing, represented by leaves on its
        both ends.

        TESTS::

            sage: B = BinaryTree([[[[None, [None, []]], None], [[], None]],
            ....:                 None])
            sage: sorted(TamariBlossomingTree._binary_tree_arcs(B))
            [(0, 3), (0, 4), (0, 7), (0, 8), (1, 3), (2, 3), (5, 6), (5, 7)]
        """
        def aux(bt, offset, arcs):
            if not bt:
                return
            arcs.append((offset, offset + bt.node_number()))
            stlist = list(bt)
            aux(stlist[0], offset, arcs)
            aux(stlist[1], offset + stlist[0].node_number() + 1, arcs)

        arclist = []
        aux(btree, 0, arclist)
        return arclist

    @staticmethod
    def binary_tree_smooth_drawing(btree, color='blue') -> Graphics:
        r"""
        Plot the smooth drawing of a binary tree, which is obtained by
        converting each internal node of a binary tree and its edges into an
        arc linking its leftmost and rightmost leaves. Such a drawing is planar.
        See [FFN2025]_ for more details.

        INPUT:

        - ``btree`` -- a binary tree, not necessarily of type
          :class:`~sage.combinat.binary_tree.BinaryTree`.
        - ``color`` -- color of the arcs. Default: ``blue``.

        OUTPUT:

        The smooth drawing of a binary tree with the given color

        EXAMPLES::

            sage: B3 = [[[[None, [None, []]], None], [[], None]], None]
            sage: g = TamariBlossomingTree.binary_tree_smooth_drawing(B3)

        .. PLOT::

            B3 = [[[[None, [None, []]], None], [[], None]], None]
            g = TamariBlossomingTree.binary_tree_smooth_drawing(B3)
            sphinx_plot(g)
        """
        # initialization
        bt = BinaryTree(btree)
        G = Graphics()
        G.set_aspect_ratio(0.5)

        # plot the arcs
        for e in TamariBlossomingTree._binary_tree_arcs(bt):
            G += arc([(e[0] + e[1]) / 2, 0], (e[1] - e[0]) / 2, sector=(0, pi),
                     rgbcolor=color)
        G.axes(show=False)
        return G

    def smooth_drawing(self) -> Graphics:
        r"""
        Plot the smooth drawing of ``self``.

        The smooth drawing of a blossoming tree is the combination of the smooth
        drawings of the two binary trees of the Tamari interval represented by
        the current blossoming tree, with that of the larger element above and
        that of the smaller one below, reflected by the horizontal axis.

        See [FFN2025]_ for more details.

        EXAMPLES::

            sage: Tl = OrderedTree([[], [[], [], [[], []], [[[], []], [], []]],
            ....:                   [[], [], [[[], []], [], [[], []], []]]])
            sage: g = TamariBlossomingTree(Tl).smooth_drawing()

        .. PLOT::

            Tl = OrderedTree([[], [[], [], [[], []], [[[], []], [], []]],
                              [[], [], [[[], []], [], [[], []], []]]])
            g = TamariBlossomingTree(Tl).smooth_drawing()
            sphinx_plot(g)
        """
        def cirnode(x, y):
            return circle([x, y], 0.1, fill=True, edgecolor='black',
                          facecolor='black', zorder=2)

        def semicir(x1, x2, isupper):
            sec = (0, pi) if isupper else (pi, 2 * pi)
            color = 'blue' if isupper else 'red'
            return arc([(x1 + x2) / 2, 0], (x2 - x1) / 2, sector=sec, zorder=1,
                       rgbcolor=color)

        # initialization
        G = Graphics()
        G.set_aspect_ratio(1)

        # draw the vertices, black circle for tree node, white squares for edges
        n = self._size
        for i in range(n + 1):  # tree nodes
            G += cirnode(i, 0)

        # draw the arcs according to upper and lower trees
        trees = self.to_tamari()
        for i in range(2):
            for e in TamariBlossomingTree._binary_tree_arcs(trees[i]):
                G += semicir(e[0], e[1], i == 1)  # 0 is lower, 1 is upper
        G.axes(show=False)
        return G

    def is_synchronized(self) -> bool:
        r"""
        Return whether the Tamari interval presented by the current blossoming
        tree is synchronized, i.e., two buds of the same node are adjacent.

        OUTPUT:

        ``True`` if the blossoming tree is synchronized, otherwise ``False``

        EXAMPLES::

            sage: T1 = [[], [[], [], [[], []]]]
            sage: TamariBlossomingTree.from_plane_tree(T1).is_synchronized()
            True
            sage: T2 = [[], [[], [[], []], []]]
            sage: TamariBlossomingTree.from_plane_tree(T2).is_synchronized()
            False
        """
        def aux(tree, isroot=False):
            """
            Check synchronized condition on subtree
            """
            # get indices of buds
            idx = [i for i, st in enumerate(tree) if not st]
            # bud number, and consecutive check
            if isroot:
                # at the root: one bud at the beginning or the end
                if len(idx) != 1 or (idx[0] != 0 and idx[0] != len(tree) - 1):
                    return False
            else:
                # otherwise: two buds consecutive
                if len(idx) != 2 or idx[1] - idx[0] != 1:
                    return False
            for st in tree:
                if st and not aux(st):  # an internal node failing the test
                    return False
            return True
        return aux(self._tree, isroot=True)

    def is_modern(self) -> bool:
        r"""
        Return whether the Tamari interval associated to the current blossoming
        tree is modern, using the function
        :meth:`~sage.combinat.interval_posets.TamariIntervalPoset.is_modern` in
        :class:`~sage.combinat.interval_posets.TamariIntervalPoset`.

        OUTPUT:

        ``True`` if the blossoming tree is modern, and ``False`` otherwise.

        EXAMPLES::

            sage: T1 = [[], [[], [], [[], []]]]
            sage: TamariBlossomingTree.from_plane_tree(T1).is_modern()
            True
            sage: T2 = [[], [[[], [], [[], []]], [], []]]
            sage: TamariBlossomingTree.from_plane_tree(T2).is_modern()
            False
        """
        return self.to_TIP().is_modern()


class _RandomPath:
    r"""
    This class contains static functions related to the generation of random
    positive lattice paths with steps `(1, k - 1)` and `(1, -1)` of length
    `kn`, with `n` and `k` given in parameters.
    """

    @staticmethod
    def gen_comb(n: int, k: int) -> list[int]:
        r"""
        Generate a random combination of `n` elements among `kn + 1` elements
        using a random approach, which is faster than the unranking approach.

        INPUT:

        - ``n`` -- a positive integer
        - ``k`` -- a strictly positive integer at least 2

        OUTPUT:

        A list of elements in the random combination of ``n`` elements among
        integers from ``1`` to ``kn+1``, not necessary in order.

        EXAMPLES::

            sage: from sage.combinat.tamari_blossoming_tree import _RandomPath
            sage: l = _RandomPath.gen_comb(3, 2)
            sage: len(l)
            3
            sage: max(l) < 7 and min(l) >= 0
            True
            sage: res = {tuple(sorted(_RandomPath.gen_comb(3, 2)))
            ....:        for _ in range(1000)}
            sage: len(res) <= binomial(7, 3)
            True
            sage: while len(res) < binomial(7, 3):
            ....:     res.add(tuple(sorted(_RandomPath.gen_comb(3, 2))))
            sage: len(res) == binomial(7, 3)
            True
        """
        # test validity
        if n < 0 or k <= 1:
            raise ValueError("invalid parameter")
        # get a random set with each element appearing with prob 1/k
        # the size of the set is close to n, with sqrt(n) standard deviation
        # better than unranking in terms of performance
        s: list[int] = []   # the random set
        cs: list[int] = []  # its complement
        for i in range(k * n + 1):
            if randrange(k) == 1:
                s.append(i)
            else:
                cs.append(i)
        cnt: int = len(s)
        if cnt > n:  # First case: too many elements
            # remove randomly until getting the correct number
            while cnt > n:
                s[randrange(cnt)] = s[-1]
                s.pop()
                cnt -= 1
        else:
            # add elements randomly until getting the correct number
            while cnt < n:
                id: int = randrange(k * n + 1 - cnt)
                s.append(cs[id])
                cs[id] = cs[-1]
                cs.pop()
                cnt += 1
        return s

    @staticmethod
    def cutting(cardlist: list[float], size: int) -> list[(float, int)]:
        """
        Utility function to generate the cutting ratio list according to a
        list of (relative) count of objects of sizes from ``0`` to ``size``.
        The cutting ratio list tells us the probability to generate pairs of
        each size separation.

        More precisely, we suppose that ``cardlist[i]`` is the number of objects
        of size ``i``. Then the output is a list of pairs of elements ``(p, j)``
        such that ``p`` is the number of pairs of objects of total size
        ``size`` with the first object of size ``j``.

        This function is provided for memoization of the cutting ratio list.

        .. NOTE::

            In its usage, ``cardlist`` is sometimes renormalized with respect to
            the exponential growth, i.e., ``cardlist[n]`` may be the number of
            objects of size `n` divided by `c^n`, where `c` is the growth rate
            of ``cardlist[n]``. This does not affect the validity of the
            result, as all probabilities concern objects of the same size.

        INPUT:

        - ``cardlist`` -- a list of the cardinality of objects of sizes from
          ``0`` to ``size``
        - ``size`` -- the size of elements to generate

        OUTPUT:

        A list of pairs ``(p, j)`` of first element size and its probability.
        Note that, for more efficient random generation, the list is sorted with
        decreasing probability.

        EXAMPLES::

            sage: from sage.combinat.tamari_blossoming_tree import _RandomPath
            sage: [(int(x), y)
            ....:  for x, y in _RandomPath.cutting([1.0, 1.0, 2.0], 2)]
            [(2, 0), (2, 2), (1, 1)]
        """
        # check list size
        if len(cardlist) != size + 1:
            raise ValueError("invalid parameter: l does not have correct size.")
        cutting: list[(float, int)] = []
        for i in range(size + 1):
            cutting.append((cardlist[i] * cardlist[size - i], i))
        # sort with decreasing probability for efficient random generation
        cutting.sort(key=lambda x: x[0], reverse=True)
        return cutting

    @staticmethod
    def _comb_to_path(n: int, k: int, uset: list[int]) -> list[int]:
        """
        Utility function to convert a subset ``uset`` of integers from ``1`` to
        `nk` of size `k` into a lattice path with up steps `(1, k - 1)` and
        `(1, -1)`, while staying always weakly above the x-axis and returning to
        the x-axis at the end. They are counted by the Fuss--Catalan numbers,
        and they are generated using the cyclic lemma. More precisely, we use
        the elements in ``uset`` as positions of up step to obtain a lattice
        path of total length `nk + 1`, then find the last lowest point of the
        path and rotate it there to obtain the path we want, while removing the
        last down step.

        We note that this function assume the validity of inputs, as it should
        only be called by other internal functions, which always give inputs of
        the correct form.

        INPUT:

        - ``n`` -- the number of up steps in the path. There must be `(k - 1) n`
          down steps in the same path.
        - ``k`` -- the parameter regulating the up steps. For Dyck paths, we
          have `k = 2`.
        - ``uset`` -- the subset of positions of up steps in the non-rotated
          path.

        OUTPUT:

        A lattice path of steps `(1, k - 1)` and `(1, -1)` ending on the
        `x`-axis while staying weakly above it, represented by variation in
        `y`-coordinates.

        TESTS::

            sage: from sage.combinat.tamari_blossoming_tree import _RandomPath
            sage: _RandomPath._comb_to_path(4, 3, [1, 3, 7, 8])
            [2, -1, 2, -1, -1, -1, 2, 2, -1, -1, -1, -1]
            sage: _RandomPath._comb_to_path(4, 3, [3, 4, 10, 12])
            [2, -1, 2, -1, -1, -1, 2, 2, -1, -1, -1, -1]
        """
        path = [-1] * (k * n + 1)
        for e in uset:
            path[e] = k - 1
        # find last lowest point
        lidx, minh, height = 0, 0, 0
        for i, step in enumerate(path):
            height += step
            if height < minh:
                lidx, minh = i + 1, height
        # rotate for the path
        path = path[lidx:] + path[:lidx]
        # remove last step
        path.pop()
        return path

    @staticmethod
    def gen_path(n: int, k: int) -> list[int]:
        r"""
        Return a uniformly random lattice path staying weakly above the x-axis
        with `kn` steps , `n` of them up steps `(1, k - 1)` and others down
        steps `(1, -1)`. Such lattice paths are sometimes also called "Raney
        paths", and are counted by Fuss--Catalan numbers. For more details, see
        the documentation of ``_comb_to_path``.

        INPUT:

        - ``n`` -- the size of the lattice path (the length will be `kn`).
        - ``k`` -- the parameter corresponding to the arity of the `k`-ary trees
          that are in bijection with the generated random lattice path.

        OUTPUT:

        A uniformly random lattice path represented by variation in
        `y`-coordinates.

        TESTS::

            sage: from sage.combinat.tamari_blossoming_tree import _RandomPath
            sage: res = {tuple(_RandomPath.gen_path(3, 3)) for _ in range(1000)}
            sage: len(res) <= binomial(10, 3) / 10
            True
            sage: while len(res) * 10 < binomial(10, 3):
            ....:     res.add(tuple(_RandomPath.gen_path(3, 3)))
            sage: len(res) == binomial(10, 3) / 10
            True
        """
        uset = _RandomPath.gen_comb(n, k)
        return _RandomPath._comb_to_path(n, k, uset)


class TamariBlossomingTrees(UniqueRepresentation, Parent):
    """
    Factory class for Tamari blossoming trees. See
    :class:`~sage.combinat.tamari_blossoming_tree.TamariBlossomingTree` for
    details of the elements.

    INPUT:

    - ``size`` -- size of Tamari blossoming trees in the set (optional,
      default: None, standing for Tamari blossoming trees of all sizes)

    OUTPUT:

    An object representing the set of Tamari blossoming trees of the given size.

    EXAMPLES::

        sage: from sage.combinat.tamari_blossoming_tree import (
        ....:     TamariBlossomingTrees)
        sage: TamariBlossomingTrees()
        Tamari blossoming trees
        sage: TamariBlossomingTrees(10)
        Tamari blossoming trees of size 10

    .. NOTE::

        This is a factory class whose constructor returns instances of
        subclasses corresponding to call parameters.
    """

    def __classcall_private__(cls, size=None):
        """
        Return the appropriate subclass according to the type of the parameter.

        INPUT:

        - ``size`` -- size of Tamari blossoming trees in the set (optional,
          default: None, standing for Tamari blossoming trees of all sizes)

        We use the private version here to avoid infinite recursion when
        initializing
        :class:`~sage.combinat.tamari_blossoming_tree.TamariBlossomingTrees_all`.

        TESTS::

            sage: from sage.combinat.tamari_blossoming_tree import (
            ....: TamariBlossomingTrees_all, TamariBlossomingTrees_size)
            sage: isinstance(TamariBlossomingTrees(4),
            ....:            TamariBlossomingTrees_size)
            True
            sage: isinstance(TamariBlossomingTrees(),
            ....:            TamariBlossomingTrees_all)
            True
            sage: isinstance(TamariBlossomingTrees(),
            ....:            TamariBlossomingTrees_size)
            False
            sage: isinstance(TamariBlossomingTrees(4),
            ....:            TamariBlossomingTrees_all)
            False
            sage: TamariBlossomingTrees() is TamariBlossomingTrees_all()
            True
            sage: TamariBlossomingTrees(-3)
            Traceback (most recent call last):
            ...
            ValueError: size must be a strictly positive integer
            sage: TamariBlossomingTrees(2.6)
            Traceback (most recent call last):
            ...
            ValueError: size must be a strictly positive integer
        """
        if size is None:
            return TamariBlossomingTrees_all()
        if not isinstance(size, (Integer, int)) or size <= 0:
            raise ValueError('size must be a strictly positive integer')
        return TamariBlossomingTrees_size(size)


class TamariBlossomingTrees_all(DisjointUnionEnumeratedSets,
                                TamariBlossomingTrees):
    """
    The enumerated set of all Tamari blossoming trees.
    """

    Element = TamariBlossomingTree

    def __init__(self):
        """
        TESTS::

            sage: TBTA = TamariBlossomingTrees()
            sage: TBTA.cardinality()
            +Infinity
            sage: TBTA is TamariBlossomingTrees()
            True
            sage: it = iter(TBTA)
            sage: [next(it) for _ in range(4)]
            [Tamari blossoming tree [[], [[], []]] of size 1,
             Tamari blossoming tree [[], [[], [[], []], []]] of size 2,
             Tamari blossoming tree [[], [[], [], [[], []]]] of size 2,
             Tamari blossoming tree [[], [[], []], [[], []]] of size 2]
            sage: next(it).parent()
            Tamari blossoming trees
            sage: TestSuite(TBTA).run()  # long time (5s)
        """
        DisjointUnionEnumeratedSets.__init__(
            self, Family(PositiveIntegers(), TamariBlossomingTrees_size),
            facade=True, keepkey=False
        )

    def _repr_(self) -> str:
        """
        TESTS::

            sage: from sage.combinat.tamari_blossoming_tree import (
            ....:     TamariBlossomingTrees)
            sage: TamariBlossomingTrees()
            Tamari blossoming trees
        """
        return 'Tamari blossoming trees'

    def __contains__(self, elem) -> bool:
        """
        TEST::

            sage: from sage.combinat.tamari_blossoming_tree import (
            ....:     TamariBlossomingTrees)
            sage: TBT4 = TamariBlossomingTrees(4)
            sage: TBT4.random_element() in TamariBlossomingTrees()
            True
        """
        return isinstance(elem, self.element_class)


class TamariBlossomingTrees_size(TamariBlossomingTrees):
    """
    The enumerated set of Tamari blossoming trees of a given size.
    """

    def __init__(self, size):
        """
        TESTS::

            sage: TBTS = TamariBlossomingTrees(5)
            sage: TBTS is TamariBlossomingTrees(5)
            True
            sage: for i in range(1, 6):
            ....:     TestSuite(TamariBlossomingTrees(i)).run()
        """
        super().__init__(category=(FiniteEnumeratedSets(),))
        self._size = size

    def _repr_(self) -> str:
        """
        TESTS::

            sage: TamariBlossomingTrees(8)
            Tamari blossoming trees of size 8
        """
        return f'Tamari blossoming trees of size {self._size}'

    @lazy_attribute
    def _parent_for(self):
        """
        The parent of elements generated by ``self``. In our case, it should be
        :class:`~sage.combinat.tamari_blossoming_tree.TamariBlossomingTrees_all`
        so that the parent of every Tamari blossoming tree is the same.

        TESTS::

            sage: TamariBlossomingTrees(8)._parent_for
            Tamari blossoming trees
        """
        return TamariBlossomingTrees_all()

    @lazy_attribute
    def element_class(self):
        """
        Overriding the original method of the base class to make the element
        class uniform across all implementation of sets of Tamari blossoming
        trees.

        TESTS::

            sage: TBTS = TamariBlossomingTrees(5)
            sage: TBTS.element_class is TamariBlossomingTrees().element_class
            True
        """
        return self._parent_for.element_class

    def _element_constructor_(self, tree) -> TamariBlossomingTree:
        """
        Internal element constructor. As a side-effect of our choice to
        uniformize the parent class of all Tamari blossoming trees, this
        constructor is needed for an idempotent initialization of elements
        produced by this class.

        TESTS::

            sage: TBT5 = TamariBlossomingTrees(5)
            sage: T = TBT5.random_element()
            sage: T == TBT5(T)
            True
            sage: TBT5([[]])
            Traceback (most recent call last):
            ...
            ValueError: the given tree does not have the correct size
        """
        A = TamariBlossomingTrees_all()
        if isinstance(tree, TamariBlossomingTree):
            tree = tree._tree
        if not isinstance(tree, OrderedTree):
            tree = OrderedTree(tree)
        T = self.element_class(A, tree)
        if T.size() != self._size:
            raise ValueError('the given tree does not have the correct size')
        return T

    def __contains__(self, elem) -> bool:
        """
        TESTS::

            sage: TBT4 = TamariBlossomingTrees(4)
            sage: TBT4.random_element() in TBT4
            True
            sage: TBT4.random_element() in TamariBlossomingTrees(5)
            False
        """
        return (isinstance(elem, self.element_class)
                and elem.size() == self._size)

    def cardinality(self) -> Integer:
        r"""
        Return the cardinality of ``self``, i.e., the number of Tamari
        blossoming trees of size `n` used to initialize the instance.

        The formula, which is the same as that for Tamari intervals, is given
        in [Cha2008]_:

        .. MATH::

            \frac{2}{n (n + 1)} \binom{4n + 1}{n - 1}

        EXAMPLES::

            sage: from sage.combinat.tamari_blossoming_tree import (
            ....:     TamariBlossomingTrees)
            sage: [len(TamariBlossomingTrees(n)) for n in range(1, 6)]
            [1, 3, 13, 68, 399]
        """
        n = self._size
        return (2 * Integer(4 * n + 1).binomial(n - 1)) // (n * (n + 1))

    @staticmethod
    def __path_to_tree(path: list[int]) -> TamariBlossomingTree:
        """
        Internal function. Return a plane tree as a nested list corresponding
        to the lattice path given by ``path``.

        The lattice path we consider here should be with steps `(1, 3)` and
        `(1, -1)`, and is represented as a list of the y-coordinate variation.
        Furthermore, it must intersect with the x-axis only at both ends, so
        that it can be considered as a pair of paths by last return
        decomposition. For more information, see :meth:`random_element`.

        INPUT:

        - ``path`` -- a lattice path represented by the height variation
          between steps

        OUTPUT:

        A plane tree as nested list in bijection with the given path.

        TESTS::

            sage: from sage.combinat.tamari_blossoming_tree import (
            ....:     TamariBlossomingTrees)
            sage: TBS5 = TamariBlossomingTrees(5)
            sage: len(set(x for x in TBS5)) == TBS5.cardinality()
            True
        """
        # start with a stack of sub-trees in construction
        stack = [[0, []]]
        for step in path:
            if step == 3:  # new node
                stack.append([0, []])
            else:  # depending on type
                if stack[-1][0] < 2:  # add bud
                    stack[-1][0] += 1
                    stack[-1][1].append([])
                else:  # subtree completed
                    subtree = stack.pop()[1]
                    stack[-1][1].append(subtree)
        # Get the tree (list of subtrees for the moment)
        stack = stack[-1][1][0]
        # pop the last bud, which is always the last child
        stack.pop()
        return stack

    def random_element(self) -> TamariBlossomingTree:
        r"""
        Generate a uniformly random Tamari blossoming tree of a given size.

        OUTPUT:

        A uniformly random Tamari blossoming tree, obtained through random
        generation of lattice path.

        ALGORITHM:

        Let `A` be the set of rooted plane trees such that each internal node
        has two buds. Then a tree in `A` can be decomposed at the root into
        three sequences of sub-trees in `A`, separated by the two buds of the
        root. A Tamari blossoming tree represented as a rooted plane tree can
        thus be decomposed at the root as two sequences of sub-trees separated
        by the only bud (the other bud hangs above the root).

        It is clear that lattice paths with steps `(1, 3)` and `(1, -1)`
        returning to the x-axis while staying always weakly above it are in
        bijection with sequences of trees in `A` by separation at contacts. The
        parts between contacts are the same family of lattices walks with the
        extra condition of never touching the x-axis except at both ends. They
        are in bijection with trees in `A` by decomposition at the last
        returning to the height 3, 2 and 1.

        We perform random generation by considering a blossoming tree as a pair
        of lattice paths of a given total size, and the precomputation consists
        of storing the relative probability of how the total size is split.

        EXAMPLES::

            sage: from sage.combinat.tamari_blossoming_tree import (
            ....:     TamariBlossomingTrees)
            sage: TB = TamariBlossomingTrees(100).random_element()
            sage: TB
            Tamari blossoming tree ... of size 100
            sage: TB in TamariBlossomingTrees()
            True
            sage: TB in TamariBlossomingTrees(100)
            True
        """
        # Compute lazily the cutting
        if not hasattr(self, '_cutting'):
            # compute the size of trees, given by binomial(4n + 1, n) / (4n + 1)
            # normalized by dividing the growth factor 4^4 / 3^3
            # precision is enough, as the rest grows as n^(-3/2)
            l: list[float] = [1.0]  # no need to use numerical_approx with prec
            for i in range(1, self._size + 1):
                nextitem = l[-1] * (4 * i - 1) * (4 * i - 2) * (4 * i - 3) / 64
                nextitem /= (3 * i + 1) * i * (3 * i - 1) / 9
                l.append(nextitem)
                # counting for generation
            self._cutting = _RandomPath.cutting(l, self._size)
            self._cutting_sum = sum([x[0] for x in self._cutting])

        # get the correct size separation
        cnt = uniform(0, self._cutting_sum)
        s1 = -1
        for e in self._cutting:
            if cnt < e[0]:
                s1 = e[1]
                break
            else:
                cnt -= e[0]
        s2 = self._size - s1
        # generate the lattice path based on the pair of random paths
        p1 = _RandomPath.gen_path(s1, 4)
        p2 = _RandomPath.gen_path(s2, 4)
        path = [3] + p1 + [-1] + p2 + [-1, -1]
        # convert lattice path to blossoming tree
        tree = TamariBlossomingTrees_size.__path_to_tree(path)
        return TamariBlossomingTree._from_plane_tree(tree, skip_check=True,
                                                     random_bud=True)

    def __iter__(self) -> Iterator[TamariBlossomingTree]:
        """
        The generation is done using
        :class:`~sage.combinat.integer_lists.invlex.IntegerListsLex` to generate
        lattice paths in bijection with Tamari blossoming trees. See
        :meth:`random_element` for details of the bijection.

        .. NOTE::

            This method of iteration losses a factor of `n`, as only two leaves
            can be dangling, and all other leaves fail to be the root, so that
            each blossoming tree is generated only once. The current timing is
            therefore slower than using conversion from
            :class:`~sage.combinat.interval_posets.TamariIntervalPosets_size`.

        TESTS::

            sage: from sage.combinat.tamari_blossoming_tree import (
            ....:     TamariBlossomingTrees)
            sage: TBS = TamariBlossomingTrees
            sage: all(len(set(TBS(n))) == TBS(n).cardinality()
            ....:     for n in range(1, 6))
            True
            sage: all(len(list(TBS(n))) == TBS(n).cardinality()
            ....:     for n in range(1, 6))
            True
        """
        def ballot(m):
            if m == 0:
                yield []
                return
            for ss in IntegerListsLex(length=m * 3, floor=lambda x: x // 3 + 1,
                                      ceiling=lambda x: m, min_slope=0,
                                      check=False):
                accu = [3] * ss[0] + [-1]
                for i in range(1, len(ss)):
                    accu.extend([3] * (ss[i] - ss[i - 1]))
                    accu.append(-1)
                yield accu
            return

        for m1 in range(self._size + 1):
            m2 = self._size - m1
            for p1 in ballot(m1):
                for p2 in ballot(m2):
                    path = [3] + p1 + [-1] + p2 + [-1, -1]
                    tree = TamariBlossomingTrees_size.__path_to_tree(path)
                    try:
                        yield TamariBlossomingTree(tree)
                    except ValueError:  # bad luck, the root is not dangling
                        pass
        return

    def random_synchronized_element(self) -> TamariBlossomingTree:
        r"""
        Return a uniformly random synchronized blossoming tree in the current
        set of Tamari blossoming trees of a given size.

        OUTPUT:

        A uniformly random synchronized blossoming tree in the current instance.

        .. NOTE::

            This is a thin wrapper of the random generation in
            :class:`~sage.combinat.tamari_blossoming_tree.SynchronizedBlossomingTreeFactory`,
            which automatically cache the instance.

        EXAMPLES::

            sage: B = TamariBlossomingTrees(100).random_synchronized_element()
            sage: B
            Tamari blossoming tree ... of size 100
            sage: B.is_synchronized()
            True
        """
        if not hasattr(self, 'SBTF'):
            self.SBTF = SynchronizedBlossomingTreeFactory(self._size)
        return self.SBTF.random_element()

    def random_modern_element(self) -> TamariBlossomingTree:
        r"""
        Return a uniformly random modern blossoming tree in the current set of
        Tamari blossoming trees of a given size.

        OUTPUT:

        A uniformly random modern blossoming tree in the current instance.

        .. NOTE::

            This is a thin wrapper of the random generation facility in
            :class:`~sage.combinat.tamari_blossoming_tree.ModernBlossomingTreeFactory`,
            which automatically cache the instance, thus also the
            precomputation within.

        EXAMPLES::

            sage: B = TamariBlossomingTrees(100).random_modern_element()
            sage: B
            Tamari blossoming tree ... of size 100
            sage: B.is_modern()
            True
        """
        if not hasattr(self, 'MBTF'):
            self.MBTF = ModernBlossomingTreeFactory(self._size)
        return self.MBTF.random_element()


class SynchronizedBlossomingTreeFactory(SageObject, UniqueRepresentation):
    r"""
    This class is for uniform random generation of synchronized blossoming
    trees, which are in bijection with modern Tamari intervals, of a given
    size. No precomputation is needed here, but we keep the same convention as
    :class:`~sage.combinat.tamari_blossoming_tree.ModernBlossomingTreeFactory`
    for the methods.

    EXAMPLES::

        sage: from sage.combinat.tamari_blossoming_tree import (
        ....:     SynchronizedBlossomingTreeFactory)
        sage: SBTF = SynchronizedBlossomingTreeFactory(100)
        sage: SBTF
        Random generator of synchronized blossoming trees of size 100
        sage: SBTF.random_element()
        Tamari blossoming tree ... of size 100
        sage: SBTF.random_element().is_synchronized()
        True
    """

    def __init__(self, size: int) -> None:
        r"""
        Initialization of the generator by the size of synchronized blossoming
        trees to generate.

        INPUT:

        - ``size`` -- the size of the synchronized blossoming tree to generate.

        EXAMPLES::

            sage: from sage.combinat.tamari_blossoming_tree import (
            ....:     SynchronizedBlossomingTreeFactory)
            sage: SynchronizedBlossomingTreeFactory(-3)
            Traceback (most recent call last):
            ...
            ValueError: invalid parameter size
            sage: SBTF = SynchronizedBlossomingTreeFactory(100)
            sage: SBTF
            Random generator of synchronized blossoming trees of size 100
        """
        if size <= 0:
            raise ValueError('invalid parameter size')
        self._size = size

    def _repr_(self) -> str:
        """
        Return a string representing the instance of random generator and the
        size of objects that it is going to generate.

        TESTS::

            sage: from sage.combinat.tamari_blossoming_tree import (
            ....:     SynchronizedBlossomingTreeFactory)
            sage: SynchronizedBlossomingTreeFactory(17)
            Random generator of synchronized blossoming trees of size 17
        """
        s = 'Random generator of synchronized blossoming trees of size'
        s += f' {self._size}'
        return s

    def random_element(self) -> TamariBlossomingTree:
        r"""
        Generate a random synchronized blossoming tree of a given size.

        OUTPUT:

        A uniformly random synchronized blossoming tree, obtained through
        random generation of lattice paths.

        ALGORITHM:

        According to [FFN2025]_, a Tamari blossoming tree is synchronized if and
        only if the two buds of the same node is always side by side. In such a
        case, we may identify them, and what we need to generate is a plane tree
        where each node has one extra bud. Now, root such a plane tree with buds
        by one of the bud, we only need to generate a sequence of plane trees
        with a bud at each internal node. Let `A` be the set of such trees. A
        tree in `A` can be decomposed at the root as two sequences of trees in
        `A`, separated by the bud.

        It is clear that lattice paths with steps `(1, 2)` and `(1, -1)`
        starting and ending on the x-axis while staying weakly above it (also
        called 2-Dyck paths in some literature) are in bijection with such
        sequence of plane trees. Such a lattice path can be considered a
        sequence of primitive paths (touching only the x-axis at both ends),
        and each of them can be decomposed into two 2-Dyck paths by the last
        returning to the height 2 and 1. It is clear that `A` is in bijection
        with primitive paths by isomorphism of recursive decomposition.

        We note that no cutting is needed here, as we simply generate one
        sequence instead of a pair of them. This is why precomputing is not
        needed here.

        EXAMPLES::

            sage: from sage.combinat.tamari_blossoming_tree import (
            ....:     SynchronizedBlossomingTreeFactory)
            sage: B = SynchronizedBlossomingTreeFactory(100).random_element()
            sage: B
            Tamari blossoming tree ... of size 100
            sage: B.is_synchronized()
            True
        """
        path = _RandomPath.gen_path(self._size, 3)
        stack = [[0, []]]
        for step in path:
            if step == 2:  # new nodes
                stack.append([0, []])
            else:  # depending on type
                if stack[-1][0] == 0:  # add two buds
                    stack[-1][0] = 1
                    stack[-1][1].append([])
                    stack[-1][1].append([])
                else:  # subtree completed
                    subtree = stack.pop()[1]
                    stack[-1][1].append(subtree)
        tree = stack[-1][1]
        tree.append([])  # add the extra bud besides the root
        return TamariBlossomingTree._from_plane_tree(tree, skip_check=True,
                                                     random_bud=True)


class ModernBlossomingTreeFactory(SageObject, UniqueRepresentation):
    r"""
    This class is for uniform random generation of blossoming trees associated
    with modern Tamari intervals of a given size. As some precomputation is
    needed, it is the best practice to keep the same instance when generating
    multiple modern blossoming trees. For one-shot generation, we also provide a
    static method.

    EXAMPLES::

        sage: from sage.combinat.tamari_blossoming_tree import (
        ....:     ModernBlossomingTreeFactory)
        sage: MBTF = ModernBlossomingTreeFactory(100)
        sage: MBTF
        Random generator of modern blossoming trees of size 100
        sage: MBTF.random_element()
        Tamari blossoming tree ... of size 100
        sage: MBTF.random_element().is_modern()
        True
    """

    def __init__(self: Self, size: int) -> None:
        r"""
        Initialize a random generator of modern blossoming trees of a given size
        along with the needed precomputation. See the documentation of
        :meth:`random_element` for more precise information on the algorithm
        and the related precomputation.

        INPUT:

        - ``size`` -- the size of modern blossoming trees to generate, which is
          the number of internal edges (not including buds)

        EXAMPLES::

            sage: from sage.combinat.tamari_blossoming_tree import (
            ....:     ModernBlossomingTreeFactory)
            sage: MBTF = ModernBlossomingTreeFactory(100)
            sage: MBTF
            Random generator of modern blossoming trees of size 100
            sage: ModernBlossomingTreeFactory(-4)
            Traceback (most recent call last):
            ...
            ValueError: invalid parameter size.
        """
        if size <= 0:
            raise ValueError('invalid parameter size.')
        self._size = size
        # compute the size of trees
        # two parts, each given by the series
        # 1 + C(z) = 1 + \sum_{n \geq 1} \frac{2^{n-1}}{n+1} \binom{2n}{n} z^n
        # growth rate 8
        l: list[float] = [8.0, 1.0]  # float suffices as for other families
        for i in range(2, size + 1):
            l.append(l[-1] * (i - 0.5) / (i + 1))  # growth rate divided out
        # counting for generation
        self.cutting = _RandomPath.cutting(l, size)
        self.cutting_sum = sum([x[0] for x in self.cutting])

    def _repr_(self) -> str:
        """
        Return a string representing the instance of random generator and the
        size of objects that it is going to generate.

        TESTS::

            sage: from sage.combinat.tamari_blossoming_tree import (
            ....:     ModernBlossomingTreeFactory)
            sage: ModernBlossomingTreeFactory(16)
            Random generator of modern blossoming trees of size 16
        """
        return (f'Random generator of modern blossoming trees'
                f' of size {self._size}')

    def random_element(self) -> TamariBlossomingTree:
        r"""
        Generate a uniformly random modern blossoming tree of a given size.

        OUTPUT:

        A uniformly random modern blossoming tree obtained using bijection with
        lattice paths.

        ALGORITHM:

        According to Section 5.5 of [FFN2025]_, the generating function of
        modern blossoming trees can be written as `(1 + C)^2`, with:

        - `C = \frac{A}{1 - B}`
        - `B = \frac{z(1 + C)}{1 - B}`
        - `A = \frac{z(1 + C)^2}{1 - B} = B (1 + C)`

        Each series counts a certain family of trees with buds.

        By solving it, we know that

        - `C(z)` is the series of Dyck paths with weight 2 on every non-initial
          up-step
        - `B(z)` is the series of Dyck paths of C without touching the x-axis in
          middle
        - `A(z)` is the series of Dyck paths with weight 2 on every up-step
          except the first and the second one on the x-axis (there may be only
          one such step)

        We may then generate a modern blossoming tree by generating Dyck paths
        with the correct weights, then interpreting such weights in the proper
        way to obtain a modern blossoming tree according to the recursive
        decomposition. Details of the bijection for each family can be found in
        the documentation of the corresponding sub-functions.

        EXAMPLES::

            sage: from sage.combinat.tamari_blossoming_tree import (
            ....:     ModernBlossomingTreeFactory)
            sage: B = ModernBlossomingTreeFactory(100).random_element()
            sage: B
            Tamari blossoming tree ... of size 100
            sage: B.is_modern()
            True
        """
        def genC(dtree: OrderedTree) -> list[OrderedTree]:
            r"""
            Generate a forest counted by the series `1 + C(z)`, which
            is a sequence of B-trees followed by an A-tree, given a Dyck path
            and the colors of the up-steps on the x-axis. We represent the Dyck
            path by a plane tree so that the colors can be generated on the fly.
            """
            if not dtree:  # empty tree
                return []
            if len(dtree) == 1 and not dtree[0]:  # simple tree
                return [OrderedTree([[], []])]
            idx = len(dtree) - 1  # index of cutting
            while idx > 0:  # never check the first
                if randrange(2) == 1:  # bad color, B stops here
                    break
                idx -= 1
            idx += 1
            lasttree = genA(OrderedTree(dtree[:idx]))
            treelist = [lasttree]
            treelist.extend([genB(dtree[i]) for i in range(idx, len(dtree))])
            return treelist

        def genB(dtree: OrderedTree) -> OrderedTree:
            r"""
            Generate a `B`-tree, i.e., a tree counted by `B(z)`, given a Dyck
            path counted by `B(z)` (colors generated on the fly). The Dyck path
            (without the initial and final step) is represented by a plane tree.
            """
            if not dtree:  # empty tree
                return OrderedTree([[], []])
            idx = 0  # index of cutting
            while idx < len(dtree):
                if randrange(2) == 1:  # bad color, B stops here
                    break
                idx += 1
            # first sequence: a sequence of B
            treelist = [genB(dtree[i]) for i in range(idx)]
            # the first bud
            treelist.append(OrderedTree([]))
            # second sequence: a sequence given by C
            treelist.extend(genC(OrderedTree(dtree[idx:])))
            # the second bud
            treelist.append(OrderedTree([]))
            return OrderedTree(treelist)

        def genA(dtree: OrderedTree) -> OrderedTree:
            r"""
            Generate an `A`-tree, i.e. a tree counted by `A(z)`, given a Dyck
            path counted by `A(z)` with colors given on the fly (so no color
            for the two up-steps starting on the x-axis). The Dyck path is
            again given as a plane tree.
            """
            if not dtree:  # empty tree, should not happen!
                raise ValueError('internal error on genA')
            # first part: same as B, and there is already a separating bud
            treelist = list(genB(dtree[0]))
            # second part: a sequence from C
            treelist.extend(genC(OrderedTree(dtree[1:])))
            return OrderedTree(treelist)

        # get the size separation for (1 + C)^2
        cnt = uniform(0, self.cutting_sum)
        s1, s2 = -1, -1
        for e in self.cutting:
            if cnt < e[0]:
                s1 = e[1]
                break
            else:
                cnt -= e[0]
        s2 = self._size - s1
        # first C sequence
        l1 = genC(DyckWords(s1).random_element().to_ordered_tree())
        # a bud
        l1.append(OrderedTree([]))
        # second C sequence
        l2 = genC(DyckWords(s2).random_element().to_ordered_tree())
        l1.extend(l2)
        tree = OrderedTree(l1)
        return TamariBlossomingTree._from_plane_tree(tree, skip_check=True,
                                                     random_bud=True)
