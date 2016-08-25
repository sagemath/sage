# -*- coding: utf-8 -*-
"""
Binary Trees

This module deals with binary trees as mathematical (in particular immutable)
objects.

.. NOTE::

    If you need the data-structure for example to represent sets or hash
    tables with AVL trees, you should have a look at :mod:`sage.misc.sagex_ds`.

AUTHORS:

- Florent Hivert (2010-2011): initial implementation.
- Adrien Boussicault (2015): Hook statistics.

REFERENCES:

.. [LodayRonco] Jean-Louis Loday and Maria O. Ronco.
   *Hopf algebra of the planar binary trees*,
   Advances in Mathematics, volume 139, issue 2,
   10 November 1998, pp. 293-309.
   http://www.sciencedirect.com/science/article/pii/S0001870898917595

.. [HNT05] Florent Hivert, Jean-Christophe Novelli, and Jean-Yves Thibon.
   *The algebra of binary search trees*,
   :arxiv:`math/0401089v2`.

.. [CP12] Gregory Chatel, Viviane Pons.
   *Counting smaller trees in the Tamari order*,
   :arxiv:`1212.0751v1`.
"""
#*****************************************************************************
#       Copyright (C) 2010 Florent Hivert <Florent.Hivert@univ-rouen.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
# python3
from __future__ import division, absolute_import

from sage.structure.list_clone import ClonableArray
from sage.combinat.abstract_tree import (AbstractClonableTree,
                                         AbstractLabelledClonableTree)
from sage.combinat.ordered_tree import LabelledOrderedTrees
from sage.rings.integer import Integer
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.lazy_attribute import lazy_attribute, lazy_class_attribute
from sage.combinat.combinatorial_map import combinatorial_map
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.misc.cachefunc import cached_method
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets


class BinaryTree(AbstractClonableTree, ClonableArray):
    """
    Binary trees.

    Binary trees here mean ordered (a.k.a. plane) finite binary
    trees, where "ordered" means that the children of each node are
    ordered.

    Binary trees contain nodes and leaves, where each node has two
    children while each leaf has no children. The number of leaves
    of a binary tree always equals the number of nodes plus `1`.

    INPUT:

    - ``children`` -- ``None`` (default) or a list, tuple or iterable of
      length `2` of binary trees or convertible objects. This corresponds
      to the standard recursive definition of a binary tree as either a
      leaf or a pair of binary trees. Syntactic sugar allows leaving out
      all but the outermost calls of the ``BinaryTree()`` constructor, so
      that, e. g., ``BinaryTree([BinaryTree(None),BinaryTree(None)])`` can
      be shortened to ``BinaryTree([None,None])``. It is also allowed to
      abbreviate ``[None, None]`` by ``[]``.

    - ``check`` -- (default: ``True``) whether check for binary should be
      performed or not.

    EXAMPLES::

        sage: BinaryTree()
        .
        sage: BinaryTree(None)
        .
        sage: BinaryTree([])
        [., .]
        sage: BinaryTree([None, None])
        [., .]
        sage: BinaryTree([None, []])
        [., [., .]]
        sage: BinaryTree([[], None])
        [[., .], .]
        sage: BinaryTree("[[], .]")
        [[., .], .]
        sage: BinaryTree([None, BinaryTree([None, None])])
        [., [., .]]

        sage: BinaryTree([[], None, []])
        Traceback (most recent call last):
        ...
        ValueError: this is not a binary tree

    TESTS::

        sage: t1 = BinaryTree([[None, [[],[[], None]]],[[],[]]])
        sage: t2 = BinaryTree([[[],[]],[]])
        sage: with t1.clone() as t1c:
        ....:     t1c[1,1,1] = t2
        sage: t1 == t1c
        False
    """
    __metaclass__ = InheritComparisonClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        """
        Ensure that binary trees created by the enumerated sets and directly
        are the same and that they are instances of :class:`BinaryTree`.

        TESTS::

            sage: from sage.combinat.binary_tree import BinaryTrees_all
            sage: issubclass(BinaryTrees_all().element_class, BinaryTree)
            True
            sage: t0 = BinaryTree([[],[[], None]])
            sage: t0.parent()
            Binary trees
            sage: type(t0)
            <class 'sage.combinat.binary_tree.BinaryTrees_all_with_category.element_class'>

            sage: t1 = BinaryTrees()([[],[[], None]])
            sage: t1.parent() is t0.parent()
            True
            sage: type(t1) is type(t0)
            True

            sage: t1 = BinaryTrees(4)([[],[[], None]])
            sage: t1.parent() is t0.parent()
            True
            sage: type(t1) is type(t0)
            True
        """
        return cls._auto_parent.element_class(cls._auto_parent, *args, **opts)

    @lazy_class_attribute
    def _auto_parent(cls):
        """
        The automatic parent of the elements of this class.

        When calling the constructor of an element of this class, one needs a
        parent. This class attribute specifies which parent is used.

        EXAMPLES::

            sage: BinaryTree._auto_parent
            Binary trees
            sage: BinaryTree([None, None]).parent()
            Binary trees
        """
        return BinaryTrees_all()

    def __init__(self, parent, children = None, check = True):
        """
        TESTS::

            sage: BinaryTree([None, None]).parent()
            Binary trees
            sage: BinaryTree("[., [., [., .]]]")
            [., [., [., .]]]
            sage: BinaryTree("[.,.,.]")
            Traceback (most recent call last):
            ...
            ValueError: this is not a binary tree
            sage: all(BinaryTree(repr(bt)) == bt for i in range(6) for bt in BinaryTrees(i))
            True
        """
        if (isinstance(children, str)):  # if the input is the repr of a binary tree
            children = children.replace(".", "None")
            from ast import literal_eval
            children = literal_eval(children)
        if children is None:
            children = []
        elif (children == [] or children == ()
              or isinstance(children, (Integer, int))):
            children = [None, None]
        if (children.__class__ is self.__class__ and
                children.parent() == parent):
            children = list(children)
        else:
            children = [self.__class__(parent, x) for x in children]
        ClonableArray.__init__(self, parent, children, check=check)

    def check(self):
        """
        Check that ``self`` is a binary tree.

        EXAMPLES::

            sage: BinaryTree([[], []])     # indirect doctest
            [[., .], [., .]]
            sage: BinaryTree([[], [], []]) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: this is not a binary tree
            sage: BinaryTree([[]])         # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: this is not a binary tree
        """
        if not (not self or len(self) == 2):
            raise ValueError("this is not a binary tree")

    def _repr_(self):
        """
        TESTS::

            sage: t1 = BinaryTree([[], None]); t1  # indirect doctest
            [[., .], .]
            sage: BinaryTree([[None, t1], None])   # indirect doctest
            [[., [[., .], .]], .]
        """
        if not self:
            return "."
        else:
            return super(BinaryTree, self)._repr_()

    def _ascii_art_( self ):
        r"""
        TESTS::

            sage: ascii_art(BinaryTree())
            <BLANKLINE>
            sage: ascii_art(BinaryTree([]))
            o
            sage: for bt in BinaryTrees(3):
            ....:     print(ascii_art(bt))
            o
             \
              o
               \
                o
            o
             \
              o
             /
            o
              o
             / \
            o   o
              o
             /
            o
             \
              o
                o
               /
              o
             /
            o
            sage: ascii_art(BinaryTree([None,[]]))
            o
             \
              o
            sage: ascii_art(BinaryTree([None,[None,[]]]))
            o
             \
              o
               \
                o
            sage: ascii_art(BinaryTree([None,[[],None]]))
            o
             \
              o
             /
            o
            sage: ascii_art(BinaryTree([None,[[[],[]],[]]]))
               o
                \
                _o_
               /   \
              o     o
             / \
            o   o
            sage: ascii_art(BinaryTree([None,[[None,[[],[]]],None]]))
            o
             \
              o
             /
            o
             \
              o
             / \
            o   o
            sage: ascii_art(BinaryTree([[],None]))
              o
             /
            o
            sage: ascii_art(BinaryTree([[[[],None], None],None]))
                  o
                 /
                o
               /
              o
             /
            o
            sage: ascii_art(BinaryTree([[[],[]],None]))
                o
               /
              o
             / \
            o   o
            sage: ascii_art(BinaryTree([[[None,[]],[[[],None],None]], None]))
                   o
                  /
              ___o___
             /       \
            o         o
             \       /
              o     o
                   /
                  o
            sage: ascii_art(BinaryTree([[None,[[],[]]],None]))
              o
             /
            o
             \
              o
             / \
            o   o
            sage: ascii_art(BinaryTree([[],[]]))
              o
             / \
            o   o
            sage: ascii_art(BinaryTree([[],[[],None]]))
              _o_
             /   \
            o     o
                 /
                o
            sage: ascii_art(BinaryTree([[None,[]],[[[],None],None]]))
              ___o___
             /       \
            o         o
             \       /
              o     o
                   /
                  o
            sage: ascii_art(BinaryTree([[[],[]],[[],None]]))
                __o__
               /     \
              o       o
             / \     /
            o   o   o
            sage: ascii_art(BinaryTree([[[],[]],[[],[]]]))
                __o__
               /     \
              o       o
             / \     / \
            o   o   o   o
            sage: ascii_art(BinaryTree([[[[],[]],[[],[]]],[]]))
                    ___o___
                   /       \
                __o__       o
               /     \
              o       o
             / \     / \
            o   o   o   o
            sage: ascii_art(BinaryTree([[],[[[[],[]],[[],[]]],[]]]))
              _____o______
             /            \
            o           ___o___
                       /       \
                    __o__       o
                   /     \
                  o       o
                 / \     / \
                o   o   o   o
        """
        node_to_str = lambda bt: str(bt.label()) if hasattr(bt, "label") else "o"

        if self.is_empty():
            from sage.typeset.ascii_art import empty_ascii_art
            return empty_ascii_art

        from sage.typeset.ascii_art import AsciiArt
        if self[0].is_empty() and self[1].is_empty():
            bt_repr = AsciiArt( [node_to_str(self)] )
            bt_repr._root = 1
            return bt_repr
        if self[0].is_empty():
            node = node_to_str(self)
            rr_tree = self[1]._ascii_art_()
            if rr_tree._root > 2:
                f_line = " " ** Integer( rr_tree._root - 3 ) + node
                s_line = " " ** Integer( len( node ) + rr_tree._root - 3 ) + "\\"
                t_repr = AsciiArt( [f_line, s_line] ) * rr_tree
                t_repr._root = rr_tree._root - 2
            else:
                f_line = node
                s_line = " " + "\\"
                t_line = " " ** Integer( len( node ) + 1 )
                t_repr = AsciiArt( [f_line, s_line] ) * ( AsciiArt( [t_line] ) + rr_tree )
                t_repr._root = rr_tree._root
            t_repr._baseline = t_repr._h - 1
            return t_repr
        if self[1].is_empty():
            node = node_to_str(self)
            lr_tree = self[0]._ascii_art_()
            f_line = " " ** Integer( lr_tree._root + 1 ) + node
            s_line = " " ** Integer( lr_tree._root ) + "/"
            t_repr = AsciiArt( [f_line, s_line] ) * lr_tree
            t_repr._root = lr_tree._root + 2
            t_repr._baseline = t_repr._h - 1
            return t_repr
        node = node_to_str(self)
        lr_tree = self[0]._ascii_art_()
        rr_tree = self[1]._ascii_art_()
        nb_ = lr_tree._l - lr_tree._root + rr_tree._root - 1
        nb_L = nb_ // 2
        nb_R = nb_L + ( 1 if nb_ % 2 == 1 else 0 )
        f_line = " " ** Integer( lr_tree._root + 1 ) + "_" ** Integer( nb_L ) + node
        f_line += "_" ** Integer( nb_R )
        s_line = " " ** Integer( lr_tree._root ) + "/" + " " ** Integer( len( node ) + rr_tree._root - 1 + ( lr_tree._l - lr_tree._root ) ) + "\\"
        t_repr = AsciiArt( [f_line, s_line] ) * ( lr_tree + AsciiArt( [" " ** Integer( len( node ) + 2 )] ) + rr_tree )
        t_repr._root = lr_tree._root + nb_L + 2
        t_repr._baseline = t_repr._h - 1
        return t_repr

    def is_empty(self):
        """
        Return whether ``self`` is empty.

        The notion of emptiness employed here is the one which defines
        a binary tree to be empty if its root is a leaf. There is
        precisely one empty binary tree.

        EXAMPLES::

            sage: BinaryTree().is_empty()
            True
            sage: BinaryTree([]).is_empty()
            False
            sage: BinaryTree([[], None]).is_empty()
            False
        """
        return not self

    def graph(self, with_leaves=True):
        """
        Convert ``self`` to a digraph. By default, this graph contains
        both nodes and leaves, hence is never empty. To obtain a graph
        which contains only the nodes, the ``with_leaves`` optional
        keyword variable has to be set to ``False``.

        INPUT:

        - ``with_leaves`` -- (default: ``True``) a Boolean, determining
          whether the resulting graph will be formed from the leaves
          and the nodes of ``self`` (if ``True``), or only from the
          nodes of ``self`` (if ``False``)

        EXAMPLES::

            sage: t1 = BinaryTree([[], None])
            sage: t1.graph()
            Digraph on 5 vertices
            sage: t1.graph(with_leaves=False)
            Digraph on 2 vertices

            sage: t1 = BinaryTree([[], [[], None]])
            sage: t1.graph()
            Digraph on 9 vertices
            sage: t1.graph().edges()
            [(0, 1, None), (0, 4, None), (1, 2, None), (1, 3, None), (4, 5, None), (4, 8, None), (5, 6, None), (5, 7, None)]
            sage: t1.graph(with_leaves=False)
            Digraph on 4 vertices
            sage: t1.graph(with_leaves=False).edges()
            [(0, 1, None), (0, 2, None), (2, 3, None)]

            sage: t1 = BinaryTree()
            sage: t1.graph()
            Digraph on 1 vertex
            sage: t1.graph(with_leaves=False)
            Digraph on 0 vertices

            sage: BinaryTree([]).graph()
            Digraph on 3 vertices
            sage: BinaryTree([]).graph(with_leaves=False)
            Digraph on 1 vertex

            sage: t1 = BinaryTree([[], [[], []]])
            sage: t1.graph(with_leaves=False)
            Digraph on 5 vertices
            sage: t1.graph(with_leaves=False).edges()
            [(0, 1, None), (0, 2, None), (2, 3, None), (2, 4, None)]
        """
        from sage.graphs.graph import DiGraph

        if with_leaves:   # We want leaves and nodes.

            # Special treatment for the case when self is empty.
            # In this case, rec(self, 0) would give a false result.
            if not self:
                return DiGraph({0: []})

            res = DiGraph()
            # The edge set of res will be built up step by step using the
            # following function:

            def rec(tr, idx):
                if not tr:  # tr is a leaf.
                    return
                else:  # tr is a node.
                    nbl = 2 * tr[0].node_number() + 1
                    res.add_edges([[idx, idx + 1], [idx, idx + 1 + nbl]])
                    rec(tr[0], idx + 1)
                    rec(tr[1], idx + nbl + 1)
            rec(self, 0)
            return res

        else:   # We want only the nodes.

            # Special treatment for the case when self has only 1 node.
            # In this case, the general DiGraph construction would
            # falsely yield an empty graph (since it adds nodes only
            # implicitly by adding edges).
            if self.node_number() == 1:
                return DiGraph({0: []})

            res = DiGraph()
            # The edge set of res will be built up step by step using the
            # following function:

            def rec(tr, idx):
                if not tr:  # tr is a leaf.
                    return
                else:  # tr is a node.
                    nbl = tr[0].node_number()
                    if nbl > 0:
                        res.add_edge([idx, idx + 1])
                        rec(tr[0], idx + 1)
                    if tr[1].node_number() > 0:
                        res.add_edge([idx, idx + nbl + 1])
                        rec(tr[1], idx + nbl + 1)
            rec(self, 0)
            return res

    def canonical_labelling(self, shift=1):
        r"""
        Return a labelled version of ``self``.

        The canonical labelling of a binary tree is a certain labelling of the
        nodes (not the leaves) of the tree.
        The actual canonical labelling is currently unspecified. However, it
        is guaranteed to have labels in `1...n` where `n` is the number of
        nodes of the tree. Moreover, two (unlabelled) trees compare as equal if
        and only if their canonical labelled trees compare as equal.

        EXAMPLES::

            sage: BinaryTree().canonical_labelling()
            .
            sage: BinaryTree([]).canonical_labelling()
            1[., .]
            sage: BinaryTree([[[], [[], None]], [[], []]]).canonical_labelling()
            5[2[1[., .], 4[3[., .], .]], 7[6[., .], 8[., .]]]
        """
        LTR = self.parent().labelled_trees()
        if self:
            sz0 = self[0].node_number()
            return LTR([self[0].canonical_labelling(shift),
                        self[1].canonical_labelling(shift + 1 + sz0)],
                       label=shift + sz0)
        else:
            return LTR(None)

    def show(self, with_leaves=False):
        """
        Show the binary tree ``show``, with or without leaves depending
        on the Boolean keyword variable ``with_leaves``.

        .. WARNING::

            Left and right children might get interchanged in
            the actual picture. Moreover, for a labelled binary
            tree, the labels shown in the picture are not (in
            general) the ones given by the labelling!

            Use :meth:`_latex_`, ``view``,
            :meth:`_ascii_art_` or ``pretty_print`` for more
            faithful representations of the data of the tree.

        TESTS::

            sage: t1 = BinaryTree([[], [[], None]])
            sage: t1.show()
        """
        try:
            self.graph(with_leaves=with_leaves).show(layout='tree', tree_root=0, tree_orientation="down")
        except RuntimeError:
            # This is for the border case BinaryTree().show().
            self.graph(with_leaves=with_leaves).show()

    def make_node(self, child_list = [None, None]):
        """
        Modify ``self`` so that it becomes a node with children ``child_list``.

        INPUT:

        - ``child_list`` -- a pair of binary trees (or objects convertible to)

        .. NOTE:: ``self`` must be in a mutable state.

        .. SEEALSO::

            :meth:`make_leaf <sage.combinat.binary_tree.BinaryTree.make_leaf>`

        EXAMPLES::

            sage: t = BinaryTree()
            sage: t.make_node([None, None])
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: with t.clone() as t1:
            ....:     t1.make_node([None, None])
            sage: t, t1
            (., [., .])
            sage: with t.clone() as t:
            ....:     t.make_node([BinaryTree(), BinaryTree(), BinaryTree([])])
            Traceback (most recent call last):
            ...
            ValueError: the list must have length 2
            sage: with t1.clone() as t2:
            ....:     t2.make_node([t1, t1])
            sage: with t2.clone() as t3:
            ....:     t3.make_node([t1, t2])
            sage: t1, t2, t3
            ([., .], [[., .], [., .]], [[., .], [[., .], [., .]]])
        """
        self._require_mutable()
        child_lst = [self.__class__(self.parent(), x) for x in child_list]
        if not(len(child_lst) == 2):
            raise ValueError("the list must have length 2")
        self.__init__(self.parent(), child_lst, check=False)

    def make_leaf(self):
        """
        Modify ``self`` so that it becomes a leaf (i. e., an empty tree).

        .. NOTE:: ``self`` must be in a mutable state.

        .. SEEALSO::

            :meth:`make_node <sage.combinat.binary_tree.BinaryTree.make_node>`

        EXAMPLES::

            sage: t = BinaryTree([None, None])
            sage: t.make_leaf()
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: with t.clone() as t1:
            ....:     t1.make_leaf()
            sage: t, t1
            ([., .], .)
        """
        self._require_mutable()
        self.__init__(self.parent(), None)

    def _to_dyck_word_rec(self, usemap="1L0R"):
        r"""
        EXAMPLES::

            sage: BinaryTree()._to_dyck_word_rec()
            []
            sage: BinaryTree([])._to_dyck_word_rec()
            [1, 0]
            sage: BinaryTree([[[], [[], None]], [[], []]])._to_dyck_word_rec()
            [1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0]
            sage: BinaryTree([[None,[]],None])._to_dyck_word_rec("L1R0")
            [1, 1, 0, 0, 1, 0]
            sage: BinaryTree([[[], [[], None]], [[], []]])._to_dyck_word_rec("L1R0")
            [1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0]
        """
        if self:
            w = []
            for l in usemap:
                if l == "L":
                    w += self[0]._to_dyck_word_rec(usemap)
                elif l == "R":
                    w += self[1]._to_dyck_word_rec(usemap)
                elif l == "1":
                    w += [1]
                elif l == "0":
                    w += [0]
            return w
        else:
            return []

    @combinatorial_map(name="to the Tamari corresponding Dyck path")
    def to_dyck_word_tamari(self):
        r"""
        Return the Dyck word associated with ``self`` in consistency with
        the Tamari order on Dyck words and binary trees.

        The bijection is defined recursively as follows:

        - a leaf is associated with an empty Dyck word;

        - a tree with children `l,r` is associated with the Dyck word
          `T(l) 1 T(r) 0`.

        EXAMPLES::

            sage: BinaryTree().to_dyck_word_tamari()
            []
            sage: BinaryTree([]).to_dyck_word_tamari()
            [1, 0]
            sage: BinaryTree([[None,[]],None]).to_dyck_word_tamari()
            [1, 1, 0, 0, 1, 0]
            sage: BinaryTree([[[], [[], None]], [[], []]]).to_dyck_word_tamari()
            [1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0]
        """
        return self.to_dyck_word("L1R0")

    def tamari_interval(self, other):
        r"""
        Return the Tamari interval between ``self`` and ``other`` as a
        :class:`~sage.combinat.interval_posets.TamariIntervalPoset`.

        A "Tamari interval" is an interval in the Tamari poset.
        See :meth:`tamari_lequal` for the definition of the Tamari poset.

        INPUT:

        - ``other`` -- a binary tree greater or equal to ``self``
          in the Tamari order

        EXAMPLES::

            sage: bt = BinaryTree([[None, [[], None]], None])
            sage: ip = bt.tamari_interval(BinaryTree([None, [[None, []], None]])); ip
            The Tamari interval of size 4 induced by relations [(2, 4), (3, 4), (3, 1), (2, 1)]
            sage: ip.lower_binary_tree()
            [[., [[., .], .]], .]
            sage: ip.upper_binary_tree()
            [., [[., [., .]], .]]
            sage: ip.interval_cardinality()
            4
            sage: ip.number_of_tamari_inversions()
            2
            sage: list(ip.binary_trees())
            [[., [[., [., .]], .]],
             [[., [., [., .]]], .],
             [., [[[., .], .], .]],
             [[., [[., .], .]], .]]
            sage: bt.tamari_interval(BinaryTree([[None,[]],[]]))
            Traceback (most recent call last):
            ...
            ValueError: The two binary trees are not comparable on the Tamari lattice.

        TESTS:

        Setting ``other`` equal to ``bt`` gives an interval consisting of
        just one element::

            sage: ip = bt.tamari_interval(bt)
            sage: ip
            The Tamari interval of size 4 induced by relations [(1, 4), (2, 3), (3, 4), (3, 1), (2, 1)]
            sage: list(ip.binary_trees())
            [[[., [[., .], .]], .]]

        Empty trees work well::

            sage: bt = BinaryTree()
            sage: ip = bt.tamari_interval(bt)
            sage: ip
            The Tamari interval of size 0 induced by relations []
            sage: list(ip.binary_trees())
            [.]
        """
        from sage.combinat.interval_posets import TamariIntervalPosets
        return TamariIntervalPosets.from_binary_trees(self, other)

    def tamari_join(self, other):
        r"""
        Return the join of the binary trees ``self`` and ``other``
        (of equal size) in the `n`-th Tamari poset (where `n` is
        the size of these trees).

        The `n`-th Tamari poset (defined in :meth:`tamari_lequal`)
        is known to be a lattice, and the map from the `n`-th
        symmetric group `S_n` to the `n`-th Tamari poset defined
        by sending every permutation `p \in S_n` to the binary
        search tree of `p` (more precisely, to
        ``p.binary_search_tree_shape()``) is a lattice
        homomorphism. (See Theorem 6.2 in [Read04]_.)

        .. SEEALSO::

            :meth:`tamari_lequal`, :meth:`tamari_meet`.

        AUTHORS:

        Viviane Pons and Darij Grinberg, 18 June 2014.

        EXAMPLES::

            sage: a = BinaryTree([None, [None, []]])
            sage: b = BinaryTree([None, [[], None]])
            sage: c = BinaryTree([[None, []], None])
            sage: d = BinaryTree([[[], None], None])
            sage: e = BinaryTree([[], []])
            sage: a.tamari_join(c) == a
            True
            sage: b.tamari_join(c) == b
            True
            sage: c.tamari_join(e) == a
            True
            sage: d.tamari_join(e) == e
            True
            sage: e.tamari_join(b) == a
            True
            sage: e.tamari_join(a) == a
            True

        ::

            sage: b1 = BinaryTree([None, [[[], None], None]])
            sage: b2 = BinaryTree([[[], None], []])
            sage: b1.tamari_join(b2)
            [., [[., .], [., .]]]
            sage: b3 = BinaryTree([[], [[], None]])
            sage: b1.tamari_join(b3)
            [., [., [[., .], .]]]
            sage: b2.tamari_join(b3)
            [[., .], [., [., .]]]

        The universal property of the meet operation is
        satisfied::

            sage: def test_uni_join(p, q):
            ....:     j = p.tamari_join(q)
            ....:     if not p.tamari_lequal(j):
            ....:         return False
            ....:     if not q.tamari_lequal(j):
            ....:         return False
            ....:     for r in p.tamari_greater():
            ....:         if q.tamari_lequal(r) and not j.tamari_lequal(r):
            ....:             return False
            ....:     return True
            sage: all( test_uni_join(p, q) for p in BinaryTrees(3) for q in BinaryTrees(3) )
            True
            sage: p = BinaryTrees(6).random_element(); q = BinaryTrees(6).random_element(); test_uni_join(p, q)
            True

        Border cases::

            sage: b = BinaryTree(None)
            sage: b.tamari_join(b)
            .
            sage: b = BinaryTree([])
            sage: b.tamari_join(b)
            [., .]

        REFERENCES:

        .. [Read04] Nathan Reading.
           *Cambrian Lattices*.
           :arxiv:`math/0402086v2`.
        """
        # We use Reading's result that the projection from the symmetric
        # group is a lattice homomorphism.
        a = self.to_132_avoiding_permutation()
        b = other.to_132_avoiding_permutation()
        return a.permutohedron_join(b).binary_search_tree_shape(left_to_right=False)

    def tamari_meet(self, other, side="right"):
        r"""
        Return the meet of the binary trees ``self`` and ``other``
        (of equal size) in the `n`-th Tamari poset (where `n` is
        the size of these trees).

        The `n`-th Tamari poset (defined in :meth:`tamari_lequal`)
        is known to be a lattice, and the map from the `n`-th
        symmetric group `S_n` to the `n`-th Tamari poset defined
        by sending every permutation `p \in S_n` to the binary
        search tree of `p` (more precisely, to
        ``p.binary_search_tree_shape()``) is a lattice
        homomorphism. (See Theorem 6.2 in [Read04]_.)

        .. SEEALSO::

            :meth:`tamari_lequal`, :meth:`tamari_join`.

        AUTHORS:

        Viviane Pons and Darij Grinberg, 18 June 2014.

        EXAMPLES::

            sage: a = BinaryTree([None, [None, []]])
            sage: b = BinaryTree([None, [[], None]])
            sage: c = BinaryTree([[None, []], None])
            sage: d = BinaryTree([[[], None], None])
            sage: e = BinaryTree([[], []])
            sage: a.tamari_meet(c) == c
            True
            sage: b.tamari_meet(c) == c
            True
            sage: c.tamari_meet(e) == d
            True
            sage: d.tamari_meet(e) == d
            True
            sage: e.tamari_meet(b) == d
            True
            sage: e.tamari_meet(a) == e
            True

        ::

            sage: b1 = BinaryTree([None, [[[], None], None]])
            sage: b2 = BinaryTree([[[], None], []])
            sage: b1.tamari_meet(b2)
            [[[[., .], .], .], .]
            sage: b3 = BinaryTree([[], [[], None]])
            sage: b1.tamari_meet(b3)
            [[[[., .], .], .], .]
            sage: b2.tamari_meet(b3)
            [[[[., .], .], .], .]

        The universal property of the meet operation is
        satisfied::

            sage: def test_uni_meet(p, q):
            ....:     m = p.tamari_meet(q)
            ....:     if not m.tamari_lequal(p):
            ....:         return False
            ....:     if not m.tamari_lequal(q):
            ....:         return False
            ....:     for r in p.tamari_smaller():
            ....:         if r.tamari_lequal(q) and not r.tamari_lequal(m):
            ....:             return False
            ....:     return True
            sage: all( test_uni_meet(p, q) for p in BinaryTrees(3) for q in BinaryTrees(3) )
            True
            sage: p = BinaryTrees(6).random_element(); q = BinaryTrees(6).random_element(); test_uni_meet(p, q)
            True

        Border cases::

            sage: b = BinaryTree(None)
            sage: b.tamari_meet(b)
            .
            sage: b = BinaryTree([])
            sage: b.tamari_meet(b)
            [., .]
        """
        # We use Reading's result that the projection from the symmetric
        # group is a lattice homomorphism.
        a = self.to_132_avoiding_permutation()
        b = other.to_132_avoiding_permutation()
        return a.permutohedron_meet(b).binary_search_tree_shape(left_to_right=False)

    @combinatorial_map(name="to Dyck paths: up step, left tree, down step, right tree")
    def to_dyck_word(self, usemap="1L0R"):
        r"""
        Return the Dyck word associated with ``self`` using the given map.

        INPUT:

        - ``usemap`` -- a string, either ``1L0R``, ``1R0L``, ``L1R0``, ``R1L0``

        The bijection is defined recursively as follows:

        - a leaf is associated to the empty Dyck Word

        - a tree with children `l,r` is associated with the Dyck word
          described by ``usemap`` where `L` and `R` are respectively the
          Dyck words associated with the trees `l` and `r`.

        EXAMPLES::

            sage: BinaryTree().to_dyck_word()
            []
            sage: BinaryTree([]).to_dyck_word()
            [1, 0]
            sage: BinaryTree([[[], [[], None]], [[], []]]).to_dyck_word()
            [1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0]
            sage: BinaryTree([[None,[]],None]).to_dyck_word()
            [1, 1, 0, 1, 0, 0]
            sage: BinaryTree([[None,[]],None]).to_dyck_word("1R0L")
            [1, 0, 1, 1, 0, 0]
            sage: BinaryTree([[None,[]],None]).to_dyck_word("L1R0")
            [1, 1, 0, 0, 1, 0]
            sage: BinaryTree([[None,[]],None]).to_dyck_word("R1L0")
            [1, 1, 0, 1, 0, 0]
            sage: BinaryTree([[None,[]],None]).to_dyck_word("R10L")
            Traceback (most recent call last):
            ...
            ValueError: R10L is not a correct map

        TESTS::

            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
            sage: bt == bt.to_dyck_word().to_binary_tree()
            True
            sage: bt == bt.to_dyck_word("1R0L").to_binary_tree("1R0L")
            True
            sage: bt == bt.to_dyck_word("L1R0").to_binary_tree("L1R0")
            True
            sage: bt == bt.to_dyck_word("R1L0").to_binary_tree("R1L0")
            True
        """
        from sage.combinat.dyck_word import DyckWord
        if usemap not in ["1L0R", "1R0L", "L1R0", "R1L0"]:
            raise ValueError("%s is not a correct map" % usemap)
        return DyckWord(self._to_dyck_word_rec(usemap))

    def _to_ordered_tree(self, bijection="left", root=None):
        r"""
        Internal recursive method to obtain an ordered tree from a binary
        tree.

        EXAMPLES::

            sage: bt = BinaryTree([[],[]])
            sage: bt._to_ordered_tree()
            [[], [[]]]
            sage: bt._to_ordered_tree(bijection="right")
            [[[]], []]
            sage: bt._to_ordered_tree(bijection="none")
            Traceback (most recent call last):
            ...
            ValueError: the bijection argument should be either left or right
            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
            sage: bt._to_ordered_tree()
            [[], [[], []], [[], [[]]]]
            sage: bt._to_ordered_tree(bijection="right")
            [[[[]], [[]]], [[]], []]
        """
        close_root = False
        if(root is None):
            from sage.combinat.ordered_tree import OrderedTree
            root = OrderedTree().clone()
            close_root = True
        if(self):
            left, right = self[0], self[1]
            if bijection == "left":
                root = left._to_ordered_tree(bijection=bijection, root=root)
                root.append(right._to_ordered_tree(bijection=bijection,
                                                   root=None))
            elif bijection == "right":
                root.append(left._to_ordered_tree(bijection=bijection,
                                                  root=None))
                root = right._to_ordered_tree(bijection=bijection,
                                              root=root)
            else:
                raise ValueError("the bijection argument should be either left or right")
        if close_root:
            root.set_immutable()
        return root

    @combinatorial_map(name="To ordered tree, left child = left brother")
    def to_ordered_tree_left_branch(self):
        r"""
        Return an ordered tree of size `n+1` by the following recursive rule:

        - if `x` is the left child of `y`, `x` becomes the left brother
          of `y`
        - if `x` is the right child of `y`, `x` becomes the last child
          of `y`

        EXAMPLES::

            sage: bt = BinaryTree([[],[]])
            sage: bt.to_ordered_tree_left_branch()
            [[], [[]]]
            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
            sage: bt.to_ordered_tree_left_branch()
            [[], [[], []], [[], [[]]]]
        """
        return self._to_ordered_tree()

    @combinatorial_map(name="To ordered tree, right child = right brother")
    def to_ordered_tree_right_branch(self):
        r"""
        Return an ordered tree of size `n+1` by the following recursive rule:

        - if `x` is the right child of `y`, `x` becomes the right brother
          of `y`
        - if `x` is the left child of `y`, `x` becomes the first child
          of `y`

        EXAMPLES::

            sage: bt = BinaryTree([[],[]])
            sage: bt.to_ordered_tree_right_branch()
            [[[]], []]
            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
            sage: bt.to_ordered_tree_right_branch()
            [[[[]], [[]]], [[]], []]
        """
        return self._to_ordered_tree(bijection="right")

    def _postfix_word(self, left_first = True, start = 1):
        r"""
        Internal recursive method to obtain a postfix canonical read of the
        binary tree.

        EXAMPLES::

            sage: bt = BinaryTree([[],[]])
            sage: bt._postfix_word()
            [1, 3, 2]
            sage: bt._postfix_word(left_first=False)
            [3, 1, 2]
            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
            sage: bt._postfix_word()
            [1, 3, 4, 2, 6, 8, 7, 5]
            sage: bt._postfix_word(left_first=False)
            [8, 6, 7, 3, 4, 1, 2, 5]
        """
        if not self:
            return []
        left = self[0]._postfix_word(left_first, start)
        label = start + self[0].node_number()
        right = self[1]._postfix_word(left_first, start = label + 1)
        if left_first:
            left.extend(right)
            left.append(label)
            return left
        else:
            right.extend(left)
            right.append(label)
            return right

    @combinatorial_map(name="To 312 avoiding permutation")
    def to_312_avoiding_permutation(self):
        r"""
        Return a 312-avoiding permutation corresponding to the binary tree.

        The linear extensions of a binary tree form an interval of the weak
        order called the sylvester class of the tree. This permutation is
        the minimal element of this sylvester class.

        EXAMPLES::

            sage: bt = BinaryTree([[],[]])
            sage: bt.to_312_avoiding_permutation()
            [1, 3, 2]
            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
            sage: bt.to_312_avoiding_permutation()
            [1, 3, 4, 2, 6, 8, 7, 5]

        TESTS::

            sage: bt = BinaryTree([[],[]])
            sage: bt == bt.to_312_avoiding_permutation().binary_search_tree_shape(left_to_right=False)
            True
            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
            sage: bt == bt.to_312_avoiding_permutation().binary_search_tree_shape(left_to_right=False)
            True
        """
        from sage.combinat.permutation import Permutation
        return Permutation(self._postfix_word())

    @combinatorial_map(name="To complete tree")
    def as_ordered_tree(self, with_leaves=True):
        r"""
        Return the same tree seen as an ordered tree. By default, leaves
        are transformed into actual nodes, but this can be avoided by
        setting the optional variable ``with_leaves`` to ``False``.

        EXAMPLES::

            sage: bt = BinaryTree([]); bt
            [., .]
            sage: bt.as_ordered_tree()
            [[], []]
            sage: bt.as_ordered_tree(with_leaves = False)
            []
            sage: bt = bt.canonical_labelling(); bt
            1[., .]
            sage: bt.as_ordered_tree()
            1[None[], None[]]
            sage: bt.as_ordered_tree(with_leaves=False)
            1[]
        """
        if with_leaves:
            children = [child.as_ordered_tree(with_leaves) for child in self]
        else:
            if not self:
                raise ValueError("The empty binary tree cannot be made into an ordered tree with with_leaves = False")
            children = [child.as_ordered_tree(with_leaves) for child in self if not child.is_empty()]
        if self in LabelledBinaryTrees():
            from sage.combinat.ordered_tree import LabelledOrderedTree
            return LabelledOrderedTree(children, label = self.label())
        else:
            from sage.combinat.ordered_tree import OrderedTree
            return OrderedTree(children)

    @combinatorial_map(name="To graph")
    def to_undirected_graph(self, with_leaves=False):
        r"""
        Return the undirected graph obtained from the tree nodes and edges.

        Leaves are ignored by default, but one can set ``with_leaves`` to
        ``True`` to obtain the graph of the complete tree.

        INPUT:

        - ``with_leaves`` -- (default: ``False``) a Boolean, determining
          whether the resulting graph will be formed from the leaves
          and the nodes of ``self`` (if ``True``), or only from the
          nodes of ``self`` (if ``False``)

        EXAMPLES::

            sage: bt = BinaryTree([])
            sage: bt.to_undirected_graph()
            Graph on 1 vertex
            sage: bt.to_undirected_graph(with_leaves=True)
            Graph on 3 vertices

            sage: bt = BinaryTree()
            sage: bt.to_undirected_graph()
            Graph on 0 vertices
            sage: bt.to_undirected_graph(with_leaves=True)
            Graph on 1 vertex

        If the tree is labelled, we use its labelling to label the graph.
        Otherwise, we use the graph canonical labelling which means that
        two different trees can have the same graph.

        EXAMPLES::

            sage: bt = BinaryTree([[],[None,[]]])
            sage: bt.canonical_labelling().to_undirected_graph() == bt.to_undirected_graph()
            False
            sage: BinaryTree([[],[]]).to_undirected_graph() == BinaryTree([[[],None],None]).to_undirected_graph()
            True
        """
        if (not with_leaves) and (not self):
            # this case needs extra care :(
            from sage.graphs.graph import Graph
            return Graph([])
        return self.as_ordered_tree(with_leaves).to_undirected_graph()

    def to_tilting(self):
        """
        Transform a binary tree into a tilting object.

        Let `t` be a binary tree with `n` nodes. There exists a unique
        depiction of `t` (above the diagonal) such that all leaves are
        regularly distributed on the diagonal line from `(0,0)` to
        `(n,n)` and all edges are either horizontal or vertical. This
        method provides the coordinates of this depiction, with the
        root as the top-left vertex.

        OUTPUT:

        a list of pairs of integers.

        Every vertex of the binary tree is mapped to a pair of
        integers. The conventions are the following. The root has
        coordinates `(0, n)` where `n` is the node number.
        If a vertex is the left (right) son of
        another vertex, they share the first (second) coordinate.

        EXAMPLES::

            sage: t = BinaryTrees(1)[0]
            sage: t.to_tilting()
            [(0, 1)]

            sage: for t in BinaryTrees(2):
            ....:     print(t.to_tilting())
            [(1, 2), (0, 2)]
            [(0, 1), (0, 2)]

            sage: from sage.combinat.abstract_tree import from_hexacode
            sage: t = from_hexacode('2020222002000', BinaryTrees())
            sage: print(t.to_tilting())
            [(0, 1), (2, 3), (4, 5), (6, 7), (4, 7), (8, 9), (10, 11),
            (8, 11), (4, 11), (12, 13), (4, 13), (2, 13), (0, 13)]

            sage: t2 = DyckWord([1,1,1,1,0,1,1,0,0,0,1,1,0,1,0,1,1,0,1,1,0,0,0,0,0,0]).to_binary_tree()
            sage: len(t2.to_tilting()) == t2.node_number()
            True
        """
        if not self:
            return []
        return self._to_tilting_rec()[0]

    def _to_tilting_rec(self, shift=0):
        """
        Auxiliary method for :meth:`to_tilting`.

        INPUT:

        ``shift`` -- an integer (default 0)

        OUTPUT:

        list of tilting coordinates and number of leaves

        EXAMPLES::

            sage: all(t._to_tilting_rec()[1] == 4 for t in BinaryTrees(3))
            True
        """
        u, v = self
        if u:
            resu, N_u = u._to_tilting_rec(shift=shift)
        else:
            resu = []
            N_u = 1
        if v:
            tilt_v, N_v = v._to_tilting_rec(shift=shift + N_u)
            resu.extend(tilt_v)
        else:
            N_v = 1
        N = N_u + N_v
        resu.append((shift, shift + N - 1))
        return (resu, N)
    
    @combinatorial_map(name="To poset")
    def to_poset(self, with_leaves=False, root_to_leaf=False):
        r"""
        Return the poset obtained by interpreting the tree as a Hasse
        diagram.

        The default orientation is from leaves to root but you can
        pass ``root_to_leaf=True`` to obtain the inverse orientation.

        Leaves are ignored by default, but one can set ``with_leaves`` to
        ``True`` to obtain the poset of the complete tree.

        INPUT:

        - ``with_leaves`` -- (default: ``False``) a Boolean, determining
          whether the resulting poset will be formed from the leaves
          and the nodes of ``self`` (if ``True``), or only from the
          nodes of ``self`` (if ``False``)
        - ``root_to_leaf`` -- (default: ``False``) a Boolean,
          determining whether the poset orientation should be from root
          to leaves (if ``True``) or from leaves to root (if ``False``).

        EXAMPLES::

            sage: bt = BinaryTree([])
            sage: bt.to_poset()
            Finite poset containing 1 elements
            sage: bt.to_poset(with_leaves=True)
            Finite poset containing 3 elements
            sage: P1 = bt.to_poset(with_leaves=True)
            sage: len(P1.maximal_elements())
            1
            sage: len(P1.minimal_elements())
            2
            sage: bt = BinaryTree([])
            sage: P2 = bt.to_poset(with_leaves=True,root_to_leaf=True)
            sage: len(P2.maximal_elements())
            2
            sage: len(P2.minimal_elements())
            1

        If the tree is labelled, we use its labelling to label the poset.
        Otherwise, we use the poset canonical labelling::

            sage: bt = BinaryTree([[],[None,[]]]).canonical_labelling()
            sage: bt
            2[1[., .], 3[., 4[., .]]]
            sage: bt.to_poset().cover_relations()
            [[4, 3], [3, 2], [1, 2]]

        Let us check that the empty binary tree is correctly handled::

            sage: bt = BinaryTree()
            sage: bt.to_poset()
            Finite poset containing 0 elements
            sage: bt.to_poset(with_leaves=True)
            Finite poset containing 1 elements
        """
        if (not with_leaves) and (not self):
            # this case needs extra care :(
            from sage.combinat.posets.posets import Poset
            return Poset({})
        return self.as_ordered_tree(with_leaves).to_poset(root_to_leaf)

    @combinatorial_map(name="To 132 avoiding permutation")
    def to_132_avoiding_permutation(self):
        r"""
        Return a 132-avoiding permutation corresponding to the binary tree.

        The linear extensions of a binary tree form an interval of the weak
        order called the sylvester class of the tree. This permutation is
        the maximal element of this sylvester class.

        EXAMPLES::

            sage: bt = BinaryTree([[],[]])
            sage: bt.to_132_avoiding_permutation()
            [3, 1, 2]
            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
            sage: bt.to_132_avoiding_permutation()
            [8, 6, 7, 3, 4, 1, 2, 5]

        TESTS::

            sage: bt = BinaryTree([[],[]])
            sage: bt == bt.to_132_avoiding_permutation().binary_search_tree_shape(left_to_right=False)
            True
            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
            sage: bt == bt.to_132_avoiding_permutation().binary_search_tree_shape(left_to_right=False)
            True
        """
        from sage.combinat.permutation import Permutation
        return Permutation(self._postfix_word(left_first=False))

    def left_children_node_number(self, direction='left'):
        r"""
        Return the number of nodes which are left children in ``self``.

        Every node (except the root) is either the left child or the
        right child of its parent node. The total number of nodes
        is `1` plus the number of left-children nodes plus the number of
        right-children nodes.

        INPUT:

        - ``direction`` -- either ``'left'`` (default) or ``'right'`` ; if
          set to ``'right'``, instead count nodes that are right children

        EXAMPLES::

            sage: bt = BinaryTree([[None,[[],[]]],[None,[[],None]]])
            sage: ascii_art(bt)
              __o__
             /     \
            o       o
             \       \
              o       o
             / \     /
            o   o   o
            sage: bt.left_children_node_number('left')
            3
            sage: bt.left_children_node_number('right')
            4

            sage: all([5 == 1 + bt.left_children_node_number()
            ....:             + bt.left_children_node_number('right')
            ....:     for bt in BinaryTrees(5)])
            True

        TESTS::

            sage: BinaryTree([[],None]).left_children_node_number()
            1
            sage: BinaryTree([None,[]]).left_children_node_number()
            0
            sage: BinaryTree([]).left_children_node_number()
            0
            sage: BinaryTree().left_children_node_number()
            0

            sage: BinaryTree([[],None]).left_children_node_number('right')
            0
            sage: BinaryTree([None,[]]).left_children_node_number('right')
            1
            sage: BinaryTree([]).left_children_node_number('right')
            0
            sage: BinaryTree().left_children_node_number('right')
            0
        """
        if self.is_empty():
            return 0
        res = 0
        if self[0]:
            if direction == 'left':
                res += 1
            res += self[0].left_children_node_number(direction)
        if self[1]:
            if direction == 'right':
                res += 1
            res += self[1].left_children_node_number(direction)
        return res

    @combinatorial_map(order = 2, name="Left-right symmetry")
    def left_right_symmetry(self):
        r"""
        Return the left-right symmetrized tree of ``self``.

        EXAMPLES::

            sage: BinaryTree().left_right_symmetry()
            .
            sage: BinaryTree([]).left_right_symmetry()
            [., .]
            sage: BinaryTree([[],None]).left_right_symmetry()
            [., [., .]]
            sage: BinaryTree([[None, []],None]).left_right_symmetry()
            [., [[., .], .]]
        """
        if not self:
            return BinaryTree()
        tree = [self[1].left_right_symmetry(), self[0].left_right_symmetry()]
        if(not self in LabelledBinaryTrees()):
            return BinaryTree(tree)
        return LabelledBinaryTree(tree, label = self.label())

    @combinatorial_map(order=2, name="Left border symmetry")
    def left_border_symmetry(self):
        r"""
        Return the tree where a symmetry has been applied recursively on
        all left borders. If a tree is made of three trees `[T_1, T_2,
        T_3]` on its left border, it becomes `[T_3', T_2', T_1']` where
        same symmetry has been applied to `T_1, T_2, T_3`.

        EXAMPLES::

            sage: BinaryTree().left_border_symmetry()
            .
            sage: BinaryTree([]).left_border_symmetry()
            [., .]
            sage: BinaryTree([[None,[]],None]).left_border_symmetry()
            [[., .], [., .]]
            sage: BinaryTree([[None,[None,[]]],None]).left_border_symmetry()
            [[., .], [., [., .]]]
            sage: bt = BinaryTree([[None,[None,[]]],None]).canonical_labelling()
            sage: bt
            4[1[., 2[., 3[., .]]], .]
            sage: bt.left_border_symmetry()
            1[4[., .], 2[., 3[., .]]]
        """
        if not self:
            return BinaryTree()
        border = []
        labelled = self in LabelledBinaryTrees()
        labels = []
        t = self
        while(t):
            border.append(t[1].left_border_symmetry())
            if labelled:
                labels.append(t.label())
            t = t[0]
        tree = BinaryTree()
        for r in border:
            if labelled:
                tree = LabelledBinaryTree([tree, r], label=labels.pop(0))
            else:
                tree = BinaryTree([tree, r])
        return tree

    def canopee(self):
        """
        Return the canopee of ``self``.

        The *canopee* of a non-empty binary tree `T` with `n` internal nodes is
        the list `l` of `0` and `1` of length `n-1` obtained by going along the
        leaves of `T` from left to right except the two extremal ones, writing
        `0` if the leaf is a right leaf and `1` if the leaf is a left leaf.

        EXAMPLES::

            sage: BinaryTree([]).canopee()
            []
            sage: BinaryTree([None, []]).canopee()
            [1]
            sage: BinaryTree([[], None]).canopee()
            [0]
            sage: BinaryTree([[], []]).canopee()
            [0, 1]
            sage: BinaryTree([[[], [[], None]], [[], []]]).canopee()
            [0, 1, 0, 0, 1, 0, 1]

        The number of pairs `(t_1, t_2)` of binary trees of size `n` such that
        the canopee of `t_1` is the complementary of the canopee of `t_2` is
        also the number of Baxter permutations (see [DG94]_, see
        also :oeis:`A001181`). We check this in small cases::

            sage: [len([(u,v) for u in BinaryTrees(n) for v in BinaryTrees(n)
            ....:       if map(lambda x:1-x, u.canopee()) == v.canopee()])
            ....:    for n in range(1, 5)]
            [1, 2, 6, 22]

        Here is a less trivial implementation of this::

            sage: from sage.sets.finite_set_map_cy import fibers
            sage: from sage.misc.all import attrcall
            sage: def baxter(n):
            ....:     f = fibers(lambda t: tuple(t.canopee()),
            ....:                   BinaryTrees(n))
            ....:     return sum(len(f[i])*len(f[tuple(1-x for x in i)])
            ....:                for i in f)
            sage: [baxter(n) for n in range(1, 7)]
            [1, 2, 6, 22, 92, 422]

        TESTS::

            sage: t = BinaryTree().canopee()
            Traceback (most recent call last):
            ...
            ValueError: canopee is only defined for non empty binary trees

        REFERENCES:

        .. [DG94] \S. Dulucq and O. Guibert. Mots de piles, tableaux
           standards et permutations de Baxter, proceedings of
           Formal Power Series and Algebraic Combinatorics, 1994.
        """
        if not self:
            raise ValueError("canopee is only defined for non empty binary trees")
        res = []

        def add_leaf_rec(tr):
            for i in range(2):
                if tr[i]:
                    add_leaf_rec(tr[i])
                else:
                    res.append(1 - i)
        add_leaf_rec(self)
        return res[1:-1]

    def in_order_traversal_iter(self):
        """
        The depth-first infix-order traversal iterator for the binary
        tree ``self``.

        This method iters each vertex (node and leaf alike) of the given
        binary tree following the depth-first infix order traversal
        algorithm.

        The *depth-first infix order traversal algorithm* iterates
        through a binary tree as follows::

            iterate through the left subtree (by the depth-first infix
                order traversal algorithm);
            yield the root;
            iterate through the right subtree (by the depth-first infix
                order traversal algorithm).

        For example on the following binary tree `T`, where we denote
        leaves by `a, b, c, \ldots` and nodes by `1, 2, 3, \ldots`::

            |     ____3____          |
            |    /         \         |
            |   1          __7__     |
            |  / \        /     \    |
            | a   2      _5_     8   |
            |    / \    /   \   / \  |
            |   b   c  4     6 h   i |
            |         / \   / \      |
            |        d   e f   g     |

        the depth-first infix-order traversal algorithm iterates through
        the vertices of `T` in the following order:
        `a,1,b,2,c,3,d,4,e,5,f,6,g,7,h,8,i`.

        See :meth:`in_order_traversal` for a version of this algorithm
        which not only iterates through, but actually does something at
        the vertices of tree.

        TESTS::

            sage: b = BinaryTree([[],[[],[]]]); ascii_art([b])
            [   _o_     ]
            [  /   \    ]
            [ o     o   ]
            [      / \  ]
            [     o   o ]
            sage: ascii_art(list(b.in_order_traversal_iter()))
            [ , o, ,   _o_    , , o, ,   o  , , o,  ]
            [         /   \             / \         ]
            [        o     o           o   o        ]
            [             / \                       ]
            [            o   o                      ]
            sage: ascii_art(filter(lambda node: node.label() is not None,
            ....:     b.canonical_labelling().in_order_traversal_iter()))
            [ 1,   _2_    , 3,   4  , 5 ]
            [     /   \         / \     ]
            [    1     4       3   5    ]
            [         / \               ]
            [        3   5              ]

            sage: list(BinaryTree(None).in_order_traversal_iter())
            [.]
        """
        if self.is_empty():
            yield self
            return
        # TODO:: PYTHON 3
        # yield from self[0].in_order_traversal_iter()
        for left_subtree in self[0].in_order_traversal_iter():
            yield left_subtree
        yield self
        # TODO:: PYTHON 3
        # yield from self[1].in_order_traversal_iter()
        for right_subtree in self[1].in_order_traversal_iter():
            yield right_subtree

    def in_order_traversal(self, node_action=None, leaf_action=None):
        r"""
        Explore the binary tree ``self`` using the depth-first infix-order
        traversal algorithm, executing the ``node_action`` function
        whenever traversing a node and executing the ``leaf_action``
        function whenever traversing a leaf.

        In more detail, what this method does to a tree `T` is the
        following::

            if the root of `T` is a node:
                apply in_order_traversal to the left subtree of `T`
                    (with the same node_action and leaf_action);
                apply node_action to the root of `T`;
                apply in_order_traversal to the right subtree of `T`
                    (with the same node_action and leaf_action);
            else:
                apply leaf_action to the root of `T`.

        For example on the following binary tree `T`, where we denote
        leaves by `a, b, c, \ldots` and nodes by `1, 2, 3, \ldots`::

            |     ____3____          |
            |    /         \         |
            |   1          __7__     |
            |  / \        /     \    |
            | a   2      _5_     8   |
            |    / \    /   \   / \  |
            |   b   c  4     6 h   i |
            |         / \   / \      |
            |        d   e f   g     |

        this method first applies ``leaf_action`` to `a`, then applies
        ``node_action`` to `1`, then ``leaf_action`` to `b`, then
        ``node_action`` to `2`, etc., with the vertices being traversed
        in the order `a,1,b,2,c,3,d,4,e,5,f,6,g,7,h,8,i`.

        See :meth:`in_order_traversal_iter` for a version of this
        algorithm which only iterates through the vertices rather than
        applying any function to them.

        INPUT:

        - ``node_action`` -- (optional) a function which takes a node in input
          and does something during the exploration
        - ``leaf_action`` -- (optional) a function which takes a leaf in input
          and does something during the exploration

        TESTS::

            sage: nb_leaf = 0
            sage: def l_action(_):
            ....:    global nb_leaf
            ....:    nb_leaf += 1
            sage: nb_node = 0
            sage: def n_action(_):
            ....:    global nb_node
            ....:    nb_node += 1

            sage: BinaryTree().in_order_traversal(n_action, l_action)
            sage: nb_leaf, nb_node
            (1, 0)

            sage: nb_leaf, nb_node = 0, 0
            sage: b = BinaryTree([[],[[],[]]]); b
            [[., .], [[., .], [., .]]]
            sage: b.in_order_traversal(n_action, l_action)
            sage: nb_leaf, nb_node
            (6, 5)
            sage: nb_leaf, nb_node = 0, 0
            sage: b = b.canonical_labelling()
            sage: b.in_order_traversal(n_action, l_action)
            sage: nb_leaf, nb_node
            (6, 5)
            sage: l = []
            sage: b.in_order_traversal(lambda node: l.append( node.label() ))
            sage: l
            [1, 2, 3, 4, 5]

            sage: leaf = 'a'
            sage: l = []
            sage: def l_action(_):
            ....:    global leaf, l
            ....:    l.append(leaf)
            ....:    leaf = chr( ord(leaf)+1 )
            sage: n_action = lambda node: l.append( node.label() )
            sage: b = BinaryTree([[None,[]],[[[],[]],[]]]).\
            ....:     canonical_labelling()
            sage: b.in_order_traversal(n_action, l_action)
            sage: l
            ['a', 1, 'b', 2, 'c', 3, 'd', 4, 'e', 5, 'f', 6, 'g', 7, 'h', 8,
             'i']
        """
        if leaf_action is None:
            leaf_action = lambda x: None
        if node_action is None:
            node_action = lambda x: None

        for node in self.in_order_traversal_iter():
            if node.is_empty():
                leaf_action(node)
            else:
                node_action(node)

    def tamari_lequal(self, t2):
        r"""
        Return ``True`` if ``self`` is less or equal to another binary
        tree ``t2`` (of the same size as ``self``) in the Tamari order.

        The Tamari order on binary trees of size `n` is the partial order
        on the set of all binary trees of size `n` generated by the
        following requirement:  If a binary tree `T'` is obtained by
        right rotation (see :meth:`right_rotate`) from a binary tree `T`,
        then `T < T'`.
        This not only is a well-defined partial order, but actually is
        a lattice structure on the set of binary trees of size `n`, and
        is a quotient of the weak order on the `n`-th symmetric group
        (also known as the right permutohedron order, see
        :meth:`~sage.combinat.permutation.Permutation.permutohedron_lequal`).
        See [CP12]_. The set of binary trees of size `n` equipped with
        the Tamari order is called the `n`-th Tamari poset.

        The Tamari order can equivalently be defined as follows:

        If `T` and `S` are two binary trees of size `n`, then the
        following four statements are equivalent:

        - We have `T \leq S` in the Tamari order.

        - There exist elements `t` and `s` of the Sylvester classes
          (:meth:`sylvester_class`) of `T` and `S`, respectively,
          such that `t \leq s` in the weak order on the symmetric
          group.

        - The 132-avoiding permutation corresponding to `T` (see
          :meth:`to_132_avoiding_permutation`) is `\leq` to the
          132-avoiding permutation corresponding to `S` in the weak
          order on the symmetric group.

        - The 312-avoiding permutation corresponding to `T` (see
          :meth:`to_312_avoiding_permutation`) is `\leq` to the
          312-avoiding permutation corresponding to `S` in the weak
          order on the symmetric group.

        .. SEEALSO::

            :meth:`tamari_smaller`, :meth:`tamari_greater`,
            :meth:`tamari_pred`, :meth:`tamari_succ`,
            :meth:`tamari_interval`

        EXAMPLES:

        This tree::

            |       o    |
            |      / \   |
            |     o   o  |
            |    /       |
            |   o        |
            |  / \       |
            | o   o      |

        is Tamari-`\leq` to the following tree::

            |     _o_     |
            |    /   \    |
            |   o     o   |
            |  / \     \  |
            | o   o     o |

        Checking this::

            sage: b = BinaryTree([[[[], []], None], []])
            sage: c = BinaryTree([[[],[]],[None,[]]])
            sage: b.tamari_lequal(c)
            True

        TESTS::

            sage: for T in BinaryTrees(4):
            ....:     for S in T.tamari_smaller():
            ....:         if S != T and T.tamari_lequal(S):
            ....:             print("FAILURE")
            ....:         if not S.tamari_lequal(T):
            ....:             print("FAILURE")
        """
        self_perm = self.to_312_avoiding_permutation()
        t2_perm = t2.to_312_avoiding_permutation()
        return self_perm.permutohedron_lequal(t2_perm)

    def tamari_greater(self):
        r"""
        The list of all trees greater or equal to ``self`` in the Tamari
        order.

        This is the order filter of the Tamari order generated by ``self``.

        See :meth:`tamari_lequal` for the definition of the Tamari poset.

        .. SEEALSO::

            :meth:`tamari_smaller`

        EXAMPLES:

        For example, the tree::

            |     __o__   |
            |    /     \  |
            |   o       o |
            |  / \     /  |
            | o   o   o   |

        has these trees greater or equal to it::

            |o          , o        , o        , o        ,  o       ,   o      ,|
            | \            \          \          \           \           \      |
            |  o            o          o           o         _o_        __o__   |
            |   \            \          \           \       /   \      /     \  |
            |    o            o          o          _o_    o     o    o       o |
            |     \            \        / \        /   \    \     \    \     /  |
            |      o            o      o   o      o     o    o     o    o   o   |
            |       \            \          \          /                        |
            |        o            o          o        o                         |
            |         \          /                                              |
            |          o        o                                               |
            <BLANKLINE>
            |   o        ,   o      ,   _o_      ,   _o__     ,   __o__    ,   ___o___  ,|
            |  / \          / \        /   \        /    \       /     \      /       \  |
            | o   o        o   o      o     o      o     _o_    o       o    o         o |
            |      \            \          / \          /   \    \       \    \       /  |
            |       o            o        o   o        o     o    o       o    o     o   |
            |        \            \            \            /      \            \        |
            |         o            o            o          o        o            o       |
            |          \          /                                                      |
            |           o        o                                                       |
            <BLANKLINE>
            |     _o_    ,     __o__  |
            |    /   \        /     \ |
            |   o     o      o       o|
            |  / \     \    / \     / |
            | o   o     o  o   o   o  |

        TESTS::

            sage: B = BinaryTree
            sage: b = B([None, B([None, B([None, B([])])])]);b
            [., [., [., [., .]]]]
            sage: b.tamari_greater()
            [[., [., [., [., .]]]]]
            sage: b = B([B([B([B([]), None]), None]), None]);b
            [[[[., .], .], .], .]
            sage: b.tamari_greater()
            [[., [., [., [., .]]]], [., [., [[., .], .]]],
            [., [[., .], [., .]]], [., [[., [., .]], .]],
            [., [[[., .], .], .]], [[., .], [., [., .]]],
            [[., .], [[., .], .]], [[., [., .]], [., .]],
            [[., [., [., .]]], .], [[., [[., .], .]], .],
            [[[., .], .], [., .]], [[[., .], [., .]], .],
            [[[., [., .]], .], .], [[[[., .], .], .], .]]
        """
        from sage.combinat.tools import transitive_ideal
        return transitive_ideal(lambda x: x.tamari_succ(), self)

    def tamari_pred(self):
        r"""
        Compute the list of predecessors of ``self`` in the Tamari poset.

        This list is computed by performing all left rotates possible on
        its nodes.

        See :meth:`tamari_lequal` for the definition of the Tamari poset.

        EXAMPLES:

        For this tree::

            |     __o__   |
            |    /     \  |
            |   o       o |
            |  / \     /  |
            | o   o   o   |

        the list is::

            |        o ,       _o_   |
            |       /         /   \  |
            |     _o_        o     o |
            |    /   \      /     /  |
            |   o     o    o     o   |
            |  / \        /          |
            | o   o      o           |

        TESTS::

            sage: B = BinaryTree
            sage: b = B([B([B([B([]), None]), None]), None]);b
            [[[[., .], .], .], .]
            sage: b.tamari_pred()
            []
            sage: b = B([None, B([None, B([None, B([])])])]);b
            [., [., [., [., .]]]]
            sage: b.tamari_pred()
            [[[., .], [., [., .]]], [., [[., .], [., .]]], [., [., [[., .], .]]]]
        """
        res = []
        if self.is_empty():
            return []
        if not self[1].is_empty():
            res.append(self.left_rotate())
        B = self.parent()._element_constructor_
        return (res +
                [B([g, self[1]]) for g in self[0].tamari_pred()] +
                [B([self[0], d]) for d in self[1].tamari_pred()])

    def tamari_smaller(self):
        r"""
        The list of all trees smaller or equal to ``self`` in the Tamari
        order.

        This is the order ideal of the Tamari order generated by ``self``.

        See :meth:`tamari_lequal` for the definition of the Tamari poset.

        .. SEEALSO::

            :meth:`tamari_greater`

        EXAMPLES:

        The tree::

            |     __o__   |
            |    /     \  |
            |   o       o |
            |  / \     /  |
            | o   o   o   |

        has these trees smaller or equal to it::

            |    __o__  ,       _o_  ,        o ,         o,         o,           o |
            |   /     \        /   \         /           /          /            /  |
            |  o       o      o     o      _o_          o          o            o   |
            | / \     /      /     /      /   \        / \        /            /    |
            |o   o   o      o     o      o     o      o   o      o            o     |
            |              /            / \          /          /            /      |
            |             o            o   o        o          o            o       |
            |                                      /          / \          /        |
            |                                     o          o   o        o         |
            |                                                            /          |
            |                                                           o           |

        TESTS::

            sage: B = BinaryTree
            sage: b = B([None, B([None, B([None, B([])])])]);b
            [., [., [., [., .]]]]
            sage: b.tamari_smaller()
            [[., [., [., [., .]]]], [., [., [[., .], .]]],
            [., [[., .], [., .]]], [., [[., [., .]], .]],
            [., [[[., .], .], .]], [[., .], [., [., .]]],
            [[., .], [[., .], .]], [[., [., .]], [., .]],
            [[., [., [., .]]], .], [[., [[., .], .]], .],
            [[[., .], .], [., .]], [[[., .], [., .]], .],
            [[[., [., .]], .], .], [[[[., .], .], .], .]]
            sage: b = B([B([B([B([]), None]), None]), None]);b
            [[[[., .], .], .], .]
            sage: b.tamari_smaller()
            [[[[[., .], .], .], .]]
        """
        from sage.combinat.tools import transitive_ideal
        return transitive_ideal(lambda x: x.tamari_pred(), self)

    def tamari_succ(self):
        r"""
        Compute the list of successors of ``self`` in the Tamari poset.

        This is the list of all trees obtained by a right rotate of
        one of its nodes.

        See :meth:`tamari_lequal` for the definition of the Tamari poset.

        EXAMPLES:

        The list of successors of::

            |     __o__   |
            |    /     \  |
            |   o       o |
            |  / \     /  |
            | o   o   o   |

        is::

            |   _o__     ,   ___o___  ,     _o_     |
            |  /    \       /       \      /   \    |
            | o     _o_    o         o    o     o   |
            |      /   \    \       /    / \     \  |
            |     o     o    o     o    o   o     o |
            |          /      \                     |
            |         o        o                    |

        TESTS::

            sage: B = BinaryTree
            sage: b = B([B([B([B([]), None]), None]), None]);b
            [[[[., .], .], .], .]
            sage: b.tamari_succ()
            [[[[., .], .], [., .]], [[[., .], [., .]], .], [[[., [., .]], .], .]]

            sage: b = B([])
            sage: b.tamari_succ()
            []

            sage: b = B([[],[]])
            sage: b.tamari_succ()
            [[., [., [., .]]]]
        """
        res = []
        if self.is_empty():
            return []
        B = self.parent()._element_constructor_
        if not self[0].is_empty():
            res.append(self.right_rotate())
        return (res +
                [B([g, self[1]]) for g in self[0].tamari_succ()] +
                [B([self[0], d]) for d in self[1].tamari_succ()])

    def single_edge_cut_shapes(self):
        """
        Return the list of possible single-edge cut shapes for the binary tree.

        This is used in :meth:`sage.combinat.interval_posets.TamariIntervalPoset.is_new`.

        OUTPUT:

        a list of triples `(m, i, n)` of integers

        This is a list running over all inner edges (i.e., edges
        joining two non-leaf vertices) of the binary tree. The removal
        of each inner edge defines two binary trees (connected
        components), the root-tree and the sub-tree. Thus, to every
        inner edge, we can assign three positive integers:
        `m` is the node number of the root-tree `R`, and `n` is
        the node number of the sub-tree `S`. The integer `i` is the
        index of the leaf of `R` on which `S` is grafted to obtain the
        original tree. The leaves of `R` are numbered starting from
        `1` (from left to right), hence `1 \leq i \leq m+1`.

        In fact, each of `m` and `n` determines the other, as the
        total node number of `R` and `S` is the node number of
        ``self``.

        EXAMPLES::

            sage: BT = BinaryTrees(3)
            sage: [t.single_edge_cut_shapes() for t in BT]
            [[(2, 3, 1), (1, 2, 2)],
            [(2, 2, 1), (1, 2, 2)],
            [(2, 1, 1), (2, 3, 1)],
            [(2, 2, 1), (1, 1, 2)],
            [(2, 1, 1), (1, 1, 2)]]

            sage: BT = BinaryTrees(2)
            sage: [t.single_edge_cut_shapes() for t in BT]
            [[(1, 2, 1)], [(1, 1, 1)]]

            sage: BT = BinaryTrees(1)
            sage: [t.single_edge_cut_shapes() for t in BT]
            [[]]
        """
        resu = []
        left, right = list(self)
        L = left.node_number()
        R = right.node_number()
        if L:
            resu += [(m + R + 1, i, n)
                     for m, i, n in left.single_edge_cut_shapes()]
            resu += [(R + 1, 1, L)]
        if R:
            resu += [(m + L + 1, i + L + 1, n)
                     for m, i, n in right.single_edge_cut_shapes()]
            resu += [(L + 1, L + 2, R)]
        return resu

    def comb(self, side='left'):
        r"""
        Return the comb of a tree.

        There are two combs in a binary tree: a left comb and a right comb.

        Consider all the vertices of the leftmost (resp. rightmost) branch of
        the root. The left (resp. right) comb is the list of right (resp. left)
        subtrees of each of these vertices.

        INPUT:

        - ``side`` -- (default: 'left') set to 'left' to obtain a left
          comb, and to 'right' to obtain a right comb.

        OUTPUT:

        A list of binary trees.

        EXAMPLES::

            sage: BT = BinaryTree( '.' )
            sage: [BT.comb('left'), BT.comb('right')]
            [[], []]
            sage: BT = BinaryTree( '[.,.]' )
            sage: [BT.comb('left'), BT.comb('right')]
            [[], []]
            sage: BT = BinaryTree( '[[[.,.], .], [.,.]]' )
            sage: BT.comb('left')
            [., .]
            sage: BT.comb('right')
            [.]
            sage: BT = BinaryTree( '[[[[., [., .]], .], [[., .], [[[., .], [., .]], [., .]]]], [., [[[., .], [[[., .], [., .]], .]], .]]]' )
            sage: ascii_art(BT)
                    ________o________
                   /                 \
                __o__                 o
               /     \                 \
              o     __o___              o
             /     /      \            /
            o     o       _o_       __o__
             \           /   \     /     \
              o         o     o   o       o
                       / \               /
                      o   o             o
                                       / \
                                      o   o
            sage: BT.comb('left')
            [[[., .], [[[., .], [., .]], [., .]]], ., [., .]]
            sage: ascii_art(BT.comb('left'))
            [   __o___     , , o ]
            [  /      \          ]
            [ o       _o_        ]
            [        /   \       ]
            [       o     o      ]
            [      / \           ]
            [     o   o          ]
            sage: BT.comb('right')
            [., [[., .], [[[., .], [., .]], .]]]
            sage: ascii_art(BT.comb('right'))
            [ ,   __o__   ]
            [    /     \  ]
            [   o       o ]
            [          /  ]
            [         o   ]
            [        / \  ]
            [       o   o ]
        """

        def _comb(side):
            if self.is_empty():
                return []
            tree = self[side]
            res = []
            while not tree.is_empty():
                res.append(tree[1 - side])
                tree = tree[side]
            return res
        if side == 'left':
            return _comb(0)
        elif side == 'right':
            return _comb(1)

    def hook_number(self):
        r"""
        Return the number of hooks.

        Recalling that a branch is a path from a vertex of the tree to a leaf,
        the leftmost (resp. rightmost) branch of a vertex `v` is the branch from
        `v` made only of left (resp. right) edges.

        The hook of a vertex `v` is a set of vertices formed by the
        union of `{v}`, and the vertices of its leftmost and rightmost branches.

        There is a unique way to partition the set of vertices in hooks.
        The number of hooks in such a partition is the hook number of the tree.

        We can obtain this partition recursively by extracting the root's hook
        and iterating the processus on each tree of the remaining forest.

        EXAMPLES::

            sage: BT = BinaryTree( '.' )
            sage: BT.hook_number()
            0
            sage: BT = BinaryTree( '[.,.]' )
            sage: BT.hook_number()
            1
            sage: BT = BinaryTree( '[[[.,.], .], [.,.]]' ); ascii_art(BT)
                o
               / \
              o   o
             /
            o
            sage: BT.hook_number()
            1
            sage: BT = BinaryTree( '[[[[., [., .]], .], [[., .], [[[., .], [., .]], [., .]]]], [., [[[., .], [[[., .], [., .]], .]], .]]]' )
            sage: ascii_art(BT)
                    ________o________
                   /                 \
                __o__                 o
               /     \                 \
              o     __o___              o
             /     /      \            /
            o     o       _o_       __o__
             \           /   \     /     \
              o         o     o   o       o
                       / \               /
                      o   o             o
                                       / \
                                      o   o
            sage: BT.hook_number()
            6
        """
        if self.is_empty():
            return 0
        return 1 + sum(t.hook_number()
                       for t in self.comb('left') + self.comb('right'))

    def twisting_number(self):
        r"""
        Return a pair (number of maximal left branches, number of maximal right
        branches).

        Recalling that a branch of a vertex `v` is a path from a vertex of the
        tree to a leaf, a left (resp. right) branch is a branch made only of
        left (resp. right) edges. The length of a branch is the number of edges
        composing it. A left (resp. right) branch is maximal if it is not
        included in a strictly longer left (resp. right) branch.


        OUTPUT :

        A list of two integers.

        EXAMPLES::

            sage: BT = BinaryTree( '.' )
            sage: BT.twisting_number()
            [0, 0]
            sage: BT = BinaryTree( '[.,.]' )
            sage: BT.twisting_number()
            [0, 0]
            sage: BT = BinaryTree( '[[[.,.], .], [.,.]]' ); ascii_art(BT)
                o
               / \
              o   o
             /
            o
            sage: BT.twisting_number()
            [1, 1]
            sage: BT = BinaryTree( '[[[[., [., .]], .], [[., .], [[[., .], [., .]], [., .]]]], [., [[[., .], [[[., .], [., .]], .]], .]]]' )
            sage: ascii_art(BT)
                    ________o________
                   /                 \
                __o__                 o
               /     \                 \
              o     __o___              o
             /     /      \            /
            o     o       _o_       __o__
             \           /   \     /     \
              o         o     o   o       o
                       / \               /
                      o   o             o
                                       / \
                                      o   o
            sage: BT.twisting_number()
            [5, 6]
            sage: BT = BinaryTree( '[.,[[[.,.],.],.]]' ); ascii_art(BT)
              o
               \
                o
               /
              o
             /
            o
            sage: BT.twisting_number()
            [1, 1]
        """
        tn = [0, 0]
        if self.node_number() <= 1:
            return tn

        L = self.comb('left')
        if len(L):
            tn[0] += 1
            for h in L:
                tw = BinaryTree([None, h]).twisting_number()
                tn[0] += tw[0]
                tn[1] += tw[1]

        R = self.comb('right')
        if len(R):
            tn[1] += 1
            for l in R:
                tw = BinaryTree([l, None]).twisting_number()
                tn[0] += tw[0]
                tn[1] += tw[1]
        return tn

    def q_hook_length_fraction(self, q=None, q_factor=False):
        r"""
        Compute the ``q``-hook length fraction of the binary tree ``self``,
        with an additional "q-factor" if desired.

        If `T` is a (plane) binary tree and `q` is a polynomial
        indeterminate over some ring, then the `q`-hook length fraction
        `h_{q} (T)` of `T` is defined by

        .. MATH::

            h_{q} (T)
            = \frac{[\lvert T \rvert]_q!}{\prod_{t \in T}
            [\lvert T_t \rvert]_q},

        where the product ranges over all nodes `t` of `T`, where `T_t`
        denotes the subtree of `T` consisting of `t` and its all
        descendants, and where for every tree `S`, we denote by
        `\lvert S \rvert` the number of nodes of `S`. While this
        definition only shows that `h_{q} (T)` is a rational function
        in `T`, it is in fact easy to show that `h_{q} (T)` is
        actually a polynomial in `T`, and thus makes sense when any
        element of a commutative ring is substituted for `q`.
        This can also be explicitly seen from the following recursive
        formula for `h_{q} (T)`:

        .. MATH::

            h_{q} (T)
            = \binom{ \lvert T \rvert - 1 }{ \lvert T_1 \rvert }_q
            h_{q} (T_1) h_{q} (T_2),

        where `T` is any nonempty binary tree, and `T_1` and `T_2` are
        the two child trees of the root of `T`, and where
        `\binom{a}{b}_q` denotes a `q`-binomial coefficient.

        A variation of the `q`-hook length fraction is the following
        "`q`-hook length fraction with `q`-factor":

        .. MATH::

            f_{q} (T)
            = h_{q} (T) \cdot
            \prod_{t \in T} q^{\lvert T_{\mathrm{right}(t)} \rvert},

        where for every node `t`, we denote by `\mathrm{right}(t)` the
        right child of `t`.
        This `f_{q} (T)` differs from `h_{q} (T)` only in a
        multiplicative factor, which is a power of `q`.

        When `q = 1`, both `f_{q} (T)` and `h_{q} (T)` equal the number
        of permutations whose binary search tree (see [HNT05]_ for the
        definition) is `T` (after dropping the labels). For example,
        there are `20` permutations which give a binary tree of the
        following shape::

            |     __o__   |
            |    /     \  |
            |   o       o |
            |  / \     /  |
            | o   o   o   |

        by the binary search insertion algorithm, in accordance with
        the fact that this tree satisfies `f_{1} (T) = 20`.

        When `q` is considered as a polynomial indeterminate,
        `f_{q} (T)` is the generating function for all permutations
        whose binary search tree is `T` (after dropping the labels)
        with respect to the number of inversions (i. e., the Coxeter
        length) of the permutations.

        Objects similar to `h_{q} (T)` also make sense for general
        ordered forests (rather than just binary trees), see e. g.
        [BW88]_, Theorem 9.1.

        INPUT:

        - ``q`` -- a ring element which is to be substituted as `q`
          into the `q`-hook length fraction (by default, this is
          set to be the indeterminate `q` in the polynomial ring
          `\ZZ[q]`)

        - ``q_factor`` -- a Boolean (default: ``False``) which
          determines whether to compute `h_{q} (T)` or to
          compute `f_{q} (T)` (namely, `h_{q} (T)` is obtained when
          ``q_factor == False``, and `f_{q} (T)` is obtained when
          ``q_factor == True``)

        REFERENCES:

        .. [BW88] Anders Bjoerner, Michelle L. Wachs,
           *Generalized quotients in Coxeter groups*.
           Transactions of the American Mathematical Society,
           vol. 308, no. 1, July 1988.
           http://www.ams.org/journals/tran/1988-308-01/S0002-9947-1988-0946427-X/S0002-9947-1988-0946427-X.pdf

        EXAMPLES:

        Let us start with a simple example. Actually, let us start
        with the easiest possible example -- the binary tree with
        only one vertex (which is a leaf)::

            sage: b = BinaryTree()
            sage: b.q_hook_length_fraction()
            1
            sage: b.q_hook_length_fraction(q_factor=True)
            1

        Nothing different for a tree with one node and two leaves::

            sage: b = BinaryTree([]); b
            [., .]
            sage: b.q_hook_length_fraction()
            1
            sage: b.q_hook_length_fraction(q_factor=True)
            1

        Let us get to a more interesting tree::

            sage: b = BinaryTree([[[],[]],[[],None]]); b
            [[[., .], [., .]], [[., .], .]]
            sage: b.q_hook_length_fraction()(q=1)
            20
            sage: b.q_hook_length_fraction()
            q^7 + 2*q^6 + 3*q^5 + 4*q^4 + 4*q^3 + 3*q^2 + 2*q + 1
            sage: b.q_hook_length_fraction(q_factor=True)
            q^10 + 2*q^9 + 3*q^8 + 4*q^7 + 4*q^6 + 3*q^5 + 2*q^4 + q^3
            sage: b.q_hook_length_fraction(q=2)
            465
            sage: b.q_hook_length_fraction(q=2, q_factor=True)
            3720
            sage: q = PolynomialRing(ZZ, 'q').gen()
            sage: b.q_hook_length_fraction(q=q**2)
            q^14 + 2*q^12 + 3*q^10 + 4*q^8 + 4*q^6 + 3*q^4 + 2*q^2 + 1

        Let us check the fact that `f_{q} (T)` is the generating function
        for all permutations whose binary search tree is `T` (after
        dropping the labels) with respect to the number of inversions of
        the permutations::

            sage: def q_hook_length_fraction_2(T):
            ....:     P = PolynomialRing(ZZ, 'q')
            ....:     q = P.gen()
            ....:     res = P.zero()
            ....:     for w in T.sylvester_class():
            ....:         res += q ** Permutation(w).length()
            ....:     return res
            sage: def test_genfun(i):
            ....:     return all( q_hook_length_fraction_2(T)
            ....:                 == T.q_hook_length_fraction(q_factor=True)
            ....:                 for T in BinaryTrees(i) )
            sage: test_genfun(4)
            True
        """
        from sage.combinat.q_analogues import q_binomial

        if q is None:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            from sage.rings.integer_ring import ZZ
            basering = PolynomialRing(ZZ, 'q')
            q = basering.gen()
        else:
            basering = q.base_ring()

        if q_factor:
            def product_of_subtrees(b):
                if b.is_empty():
                    return basering.one()
                b0 = b[0]
                b1 = b[1]
                return q_binomial(b.node_number() - 1, b0.node_number(), q=q) * \
                    product_of_subtrees(b0) * product_of_subtrees(b1) * \
                    q ** (b1.node_number())
        else:
            def product_of_subtrees(b):
                if b.is_empty():
                    return basering.one()
                b0 = b[0]
                b1 = b[1]
                return q_binomial(b.node_number() - 1, b0.node_number(), q=q) * \
                    product_of_subtrees(b0) * product_of_subtrees(b1)

        return product_of_subtrees(self)

    @combinatorial_map(name="Right rotate")
    def right_rotate(self):
        r"""
        Return the result of right rotation applied to the binary
        tree ``self``.

        Right rotation on binary trees is defined as follows:
        Let `T` be a binary tree such that the left child of the
        root of `T` is a node. Let `C` be the right
        child of the root of `T`, and let `A` and `B` be the
        left and right children of the left child of the root
        of `T`. (Keep in mind that nodes of trees are identified
        with the subtrees consisting of their descendants.)
        Then, the right rotation of `T` is the binary tree in
        which the left child of the root is `A`, whereas
        the right child of the root is a node whose left and right
        children are `B` and `C`. In pictures::

            |     *                      *     |
            |    / \                    / \    |
            |   *   C -right-rotate->  A   *   |
            |  / \                        / \  |
            | A   B                      B   C |

        where asterisks signify a single node each (but `A`, `B`
        and `C` might be empty).

        For example, ::

            |     o                     _o_   |
            |    /                     /   \  |
            |   o    -right-rotate->  o     o |
            |  / \                         /  |
            | o   o                       o   |
            <BLANKLINE>
            |       __o__                         _o__      |
            |      /     \                       /    \     |
            |     o       o  -right-rotate->    o     _o_   |
            |    / \                           /     /   \  |
            |   o   o                         o     o     o |
            |  /     \                               \      |
            | o       o                               o     |

        Right rotation is the inverse operation to left rotation
        (:meth:`left_rotate`).

        The right rotation operation introduced here is the one defined
        in Definition 2.1 of [CP12]_.

        .. SEEALSO::

            :meth:`left_rotate`

        EXAMPLES::

            sage: b = BinaryTree([[[],[]], None]); ascii_art([b])
            [     o ]
            [    /  ]
            [   o   ]
            [  / \  ]
            [ o   o ]
            sage: ascii_art([b.right_rotate()])
            [   _o_   ]
            [  /   \  ]
            [ o     o ]
            [      /  ]
            [     o   ]
            sage: b = BinaryTree([[[[],None],[None,[]]], []]); ascii_art([b])
            [       __o__   ]
            [      /     \  ]
            [     o       o ]
            [    / \        ]
            [   o   o       ]
            [  /     \      ]
            [ o       o     ]
            sage: ascii_art([b.right_rotate()])
            [     _o__      ]
            [    /    \     ]
            [   o     _o_   ]
            [  /     /   \  ]
            [ o     o     o ]
            [        \      ]
            [         o     ]
        """
        B = self.parent()._element_constructor_
        return B([self[0][0], B([self[0][1], self[1]])])

    @combinatorial_map(name="Left rotate")
    def left_rotate(self):
        r"""
        Return the result of left rotation applied to the binary
        tree ``self``.

        Left rotation on binary trees is defined as follows:
        Let `T` be a binary tree such that the right child of the
        root of `T` is a node. Let `A` be the left
        child of the root of `T`, and let `B` and `C` be the
        left and right children of the right child of the root
        of `T`. (Keep in mind that nodes of trees are identified
        with the subtrees consisting of their descendants.)
        Then, the left rotation of `T` is the binary tree in
        which the right child of the root is `C`, whereas
        the left child of the root is a node whose left and right
        children are `A` and `B`. In pictures::

            |   *                        *   |
            |  / \                      / \  |
            | A   *  -left-rotate->    *   C |
            |    / \                  / \    |
            |   B   C                A   B   |

        where asterisks signify a single node each (but `A`, `B`
        and `C` might be empty).

        For example, ::

            |   _o_                        o |
            |  /   \                      /  |
            | o     o  -left-rotate->    o   |
            |      /                    / \  |
            |     o                    o   o |
            <BLANKLINE>
            |       __o__                            o |
            |      /     \                          /  |
            |     o       o  -left-rotate->        o   |
            |    / \                              /    |
            |   o   o                            o     |
            |  /     \                          / \    |
            | o       o                        o   o   |
            |                                 /     \  |
            |                                o       o |

        Left rotation is the inverse operation to right rotation
        (:meth:`right_rotate`).

        .. SEEALSO::

            :meth:`right_rotate`

        EXAMPLES::

            sage: b = BinaryTree([[],[[],None]]); ascii_art([b])
            [   _o_   ]
            [  /   \  ]
            [ o     o ]
            [      /  ]
            [     o   ]
            sage: ascii_art([b.left_rotate()])
            [     o ]
            [    /  ]
            [   o   ]
            [  / \  ]
            [ o   o ]
            sage: b.left_rotate().right_rotate() == b
            True
        """
        B = self.parent()._element_constructor_
        return B([B([self[0], self[1][0]]), self[1][1]])

    @combinatorial_map(name="Over operation on Binary Trees")
    def over(self, bt):
        r"""
        Return ``self`` over ``bt``, where "over" is the ``over``
        (`/`) operation.

        If `T` and `T'` are two binary trees, then `T` over `T'`
        (written `T / T'`) is defined as the tree obtained by grafting
        `T'` on the rightmost leaf of `T`. More precisely, `T / T'` is
        defined by identifying the root of the `T'` with the rightmost
        leaf of `T`. See section 4.5 of [HNT05]_.

        If `T` is empty, then `T / T' = T'`.

        The definition of this "over" operation goes back to
        Loday-Ronco [LodRon0102066]_ (Definition 2.2), but it is
        denoted by `\backslash` and called the "under" operation there.
        In fact, trees in sage have their root at the top, contrary to
        the trees in [LodRon0102066]_ which are growing upwards. For
        this reason, the names of the over and under operations are
        swapped, in order to keep a graphical meaning.
        (Our notation follows that of section 4.5 of [HNT05]_.)

        .. SEEALSO::

            :meth:`under`

        EXAMPLES:

        Showing only the nodes of a binary tree, here is an
        example for the over operation::

            |   o       __o__       _o_         |
            |  / \  /  /     \  =  /   \        |
            | o   o   o       o   o     o       |
            |          \     /           \      |
            |           o   o           __o__   |
            |                          /     \  |
            |                         o       o |
            |                          \     /  |
            |                           o   o   |

        A Sage example::

            sage: b1 = BinaryTree([[],[[],[]]])
            sage: b2 = BinaryTree([[None, []],[]])
            sage: ascii_art((b1, b2, b1/b2))
            (   _o_    ,   _o_  ,   _o_           )
            (  /   \      /   \    /   \          )
            ( o     o    o     o  o     o_        )
            (      / \    \            /  \       )
            (     o   o    o          o    o      )
            (                               \     )
            (                               _o_   )
            (                              /   \  )
            (                             o     o )
            (                              \      )
            (                               o     )

        TESTS::

            sage: b1 = BinaryTree([[],[]]); ascii_art([b1])
            [   o   ]
            [  / \  ]
            [ o   o ]
            sage: b2 = BinaryTree([[None,[]],[[],None]]); ascii_art([b2])
            [   __o__   ]
            [  /     \  ]
            [ o       o ]
            [  \     /  ]
            [   o   o   ]
            sage: ascii_art([b1.over(b2)])
            [   _o_         ]
            [  /   \        ]
            [ o     o       ]
            [        \      ]
            [       __o__   ]
            [      /     \  ]
            [     o       o ]
            [      \     /  ]
            [       o   o   ]

        The same in the labelled case::

            sage: b1 = b1.canonical_labelling()
            sage: b2 = b2.canonical_labelling()
            sage: ascii_art([b1.over(b2)])
            [   _2_         ]
            [  /   \        ]
            [ 1     3       ]
            [        \      ]
            [       __3__   ]
            [      /     \  ]
            [     1       5 ]
            [      \     /  ]
            [       2   4   ]
        """
        B = self.parent()._element_constructor_
        if self.is_empty():
            return bt
        if hasattr(self, "label"):
            lab = self.label()
            return B([self[0], self[1].over(bt)], lab)
        else:
            return B([self[0], self[1].over(bt)])

    __div__ = over
    __truediv__ = over

    @combinatorial_map(name="Under operation on Binary Trees")
    def under(self, bt):
        r"""
        Return ``self`` under ``bt``, where "under" is the ``under``
        (`\backslash`) operation.

        If `T` and `T'` are two binary trees, then `T` under `T'`
        (written `T \backslash T'`) is defined as the tree obtained
        by grafting `T` on the leftmost leaf of `T'`. More precisely,
        `T \backslash T'` is defined by identifying the root of `T`
        with the leftmost leaf of `T'`.

        If `T'` is empty, then `T \backslash T' = T`.

        The definition of this "under" operation goes back to
        Loday-Ronco [LodRon0102066]_ (Definition 2.2), but it is
        denoted by `/` and called the "over" operation there. In fact,
        trees in sage have their root at the top, contrary to the trees
        in [LodRon0102066]_ which are growing upwards. For this reason,
        the names of the over and under operations are swapped, in
        order to keep a graphical meaning.
        (Our notation follows that of section 4.5 of [HNT05]_.)

        .. SEEALSO::

            :meth:`over`

        EXAMPLES:

        Showing only the nodes of a binary tree, here is an
        example for the under operation::

            sage: b1 = BinaryTree([[],[]])
            sage: b2 = BinaryTree([None,[]])
            sage: ascii_art((b1, b2, b1 \ b2))
            (   o  , o  ,     _o_   )
            (  / \    \      /   \  )
            ( o   o    o    o     o )
            (              / \      )
            (             o   o     )

        TESTS::

            sage: b1 = BinaryTree([[],[[None,[]],None]]); ascii_art([b1])
            [   _o_   ]
            [  /   \  ]
            [ o     o ]
            [      /  ]
            [     o   ]
            [      \  ]
            [       o ]
            sage: b2 = BinaryTree([[],[None,[]]]); ascii_art([b2])
            [   o     ]
            [  / \    ]
            [ o   o   ]
            [      \  ]
            [       o ]
            sage: ascii_art([b1.under(b2)])
            [        o_     ]
            [       /  \    ]
            [      o    o   ]
            [     /      \  ]
            [   _o_       o ]
            [  /   \        ]
            [ o     o       ]
            [      /        ]
            [     o         ]
            [      \        ]
            [       o       ]

        The same in the labelled case::

            sage: b1 = b1.canonical_labelling()
            sage: b2 = b2.canonical_labelling()
            sage: ascii_art([b1.under(b2)])
            [        2_     ]
            [       /  \    ]
            [      1    3   ]
            [     /      \  ]
            [   _2_       4 ]
            [  /   \        ]
            [ 1     5       ]
            [      /        ]
            [     3         ]
            [      \        ]
            [       4       ]
        """
        B = self.parent()._element_constructor_
        if bt.is_empty():
            return self
        lab = None
        if hasattr(bt, "label"):
            lab = bt.label()
            return B([self.under(bt[0]), bt[1]], lab)
        else:
            return B([self.under(bt[0]), bt[1]])

    _backslash_ = under

    def sylvester_class(self, left_to_right=False):
        r"""
        Iterate over the sylvester class corresponding to the binary tree
        ``self``.

        The sylvester class of a tree `T` is the set of permutations
        `\sigma` whose right-to-left binary search tree (a notion defined
        in [HNT05]_, Definition 7) is `T` after forgetting the labels.
        This is an equivalence class of the sylvester congruence (the
        congruence on words which holds two words `uacvbw` and `ucavbw`
        congruent whenever `a`, `b`, `c` are letters satisfying
        `a \leq b < c`, and extends by transitivity) on the symmetric
        group.

        For example the following tree's sylvester class consists of the
        permutations `(1,3,2)` and `(3,1,2)`::

            [   o   ]
            [  / \  ]
            [ o   o ]

        (only the nodes are drawn here).

        The right-to-left binary search tree of a word is constructed by
        an RSK-like insertion algorithm which proceeds as follows: Start
        with an empty labelled binary tree, and read the word from right
        to left. Each time a letter is read from the word, insert this
        letter in the existing tree using binary search tree insertion
        (:meth:`~sage.combinat.binary_tree.LabelledBinaryTree.binary_search_insert`).
        This is what the
        :meth:`~sage.combinat.permutation.Permutation.binary_search_tree`
        method computes if it is given the keyword
        ``left_to_right=False``.

        Here are two more descriptions of the sylvester class of a binary
        search tree:

        - The sylvester class of a binary search tree `T` is the set of
          all linear extensions of the poset corresponding to `T` (that
          is, of the poset whose Hasse diagram is `T`, with the root on
          top), provided that the nodes of `T` are labelled with
          `1, 2, \ldots, n` in a binary-search-tree way (i.e., every left
          descendant of a node has a label smaller than that of the node,
          and every right descendant of a node has a label higher than
          that of the node).

        - The sylvester class of a binary search tree `T` (with vertex
          labels `1, 2, \ldots, n`) is the interval `[u, v]` in the right
          permutohedron order
          (:meth:`~sage.combinat.permutation.Permutation.permutohedron_lequal`),
          where `u` is the 312-avoiding permutation corresponding to `T`
          (:meth:`to_312_avoiding_permutation`), and where `v` is the
          132-avoiding permutation corresponding to `T`
          (:meth:`to_132_avoiding_permutation`).

        If the optional keyword variable ``left_to_right`` is set to
        ``True``, then the *left* sylvester class of ``self`` is
        returned instead. This is the set of permutations `\sigma` whose
        left-to-right binary search tree (that is, the result of the
        :meth:`~sage.combinat.permutation.Permutation.binary_search_tree`
        with ``left_to_right`` set to ``True``) is ``self``. It is an
        equivalence class of the left sylvester congruence.

        .. WARNING::

            This method yields the elements of the sylvester class as
            raw lists, not as permutations!

        EXAMPLES:

        Verifying the claim that the right-to-left binary search trees of
        the permutations in the sylvester class of a tree `t` all equal
        `t`::

            sage: def test_bst_of_sc(n, left_to_right):
            ....:     for t in BinaryTrees(n):
            ....:         for p in t.sylvester_class(left_to_right=left_to_right):
            ....:             p_per = Permutation(p)
            ....:             tree = p_per.binary_search_tree(left_to_right=left_to_right)
            ....:             if not BinaryTree(tree) == t:
            ....:                 return False
            ....:     return True
            sage: test_bst_of_sc(4, False)
            True
            sage: test_bst_of_sc(5, False)   # long time
            True

        The same with the left-to-right version of binary search::

            sage: test_bst_of_sc(4, True)
            True
            sage: test_bst_of_sc(5, True)   # long time
            True

        Checking that the sylvester class is the set of linear extensions
        of the poset of the tree::

            sage: all( sorted(t.canonical_labelling().sylvester_class())
            ....:      == sorted(list(v) for v in t.canonical_labelling().to_poset().linear_extensions())
            ....:      for t in BinaryTrees(4) )
            True

        TESTS::

            sage: list(BinaryTree([[],[]]).sylvester_class())
            [[1, 3, 2], [3, 1, 2]]
            sage: bt = BinaryTree([[[],None],[[],[]]])
            sage: l = list(bt.sylvester_class()); l
            [[1, 2, 4, 6, 5, 3],
             [1, 4, 2, 6, 5, 3],
             [1, 4, 6, 2, 5, 3],
             [1, 4, 6, 5, 2, 3],
             [4, 1, 2, 6, 5, 3],
             [4, 1, 6, 2, 5, 3],
             [4, 1, 6, 5, 2, 3],
             [4, 6, 1, 2, 5, 3],
             [4, 6, 1, 5, 2, 3],
             [4, 6, 5, 1, 2, 3],
             [1, 2, 6, 4, 5, 3],
             [1, 6, 2, 4, 5, 3],
             [1, 6, 4, 2, 5, 3],
             [1, 6, 4, 5, 2, 3],
             [6, 1, 2, 4, 5, 3],
             [6, 1, 4, 2, 5, 3],
             [6, 1, 4, 5, 2, 3],
             [6, 4, 1, 2, 5, 3],
             [6, 4, 1, 5, 2, 3],
             [6, 4, 5, 1, 2, 3]]
            sage: len(l) == Integer(bt.q_hook_length_fraction()(q=1))
            True

        Border cases::

            sage: list(BinaryTree().sylvester_class())
            [[]]
            sage: list(BinaryTree([]).sylvester_class())
            [[1]]
        """
        if self.is_empty():
            yield []
            return
        from itertools import product
        from sage.combinat.words.word import Word as W
        from sage.combinat.words.shuffle_product import ShuffleProduct_w1w2 \
            as shuffle

        if left_to_right:
            builder = lambda i, p: [i] + list(p)
        else:
            builder = lambda i, p: list(p) + [i]

        shift = self[0].node_number() + 1
        for l, r in product(self[0].sylvester_class(left_to_right=left_to_right),
                            self[1].sylvester_class(left_to_right=left_to_right)):
            for p in shuffle(W(l), W([shift + ri for ri in r])):
                yield builder(shift, p)

    def is_full(self):
        r"""
        Return ``True`` if ``self`` is full, else return ``False``.

        A full binary tree is a tree in which every node either has two
        child nodes or has two child leaves.

        This is also known as *proper binary tree* or *2-tree* or *strictly
        binary tree*.

        For example::

            |       __o__   |
            |      /     \  |
            |     o       o |
            |    / \        |
            |   o   o       |
            |  /     \      |
            | o       o     |

        is not full but the next one is::

            |         ___o___   |
            |        /       \  |
            |     __o__       o |
            |    /     \        |
            |   o       o       |
            |  / \     / \      |
            | o   o   o   o     |

        EXAMPLES::

            sage: BinaryTree([[[[],None],[None,[]]], []]).is_full()
            False
            sage: BinaryTree([[[[],[]],[[],[]]], []]).is_full()
            True
            sage: ascii_art(filter(lambda bt: bt.is_full(), BinaryTrees(5)))
            [   _o_    ,     _o_   ]
            [  /   \        /   \  ]
            [ o     o      o     o ]
            [      / \    / \      ]
            [     o   o  o   o     ]
        """
        if self.is_empty():
            return True
        if self[0].is_empty() != self[1].is_empty():
            return False
        return self[0].is_full() and self[1].is_full()

    def is_perfect(self):
        r"""
        Return ``True`` if ``self`` is perfect, else return ``False``.

        A perfect binary tree is a full tree in which all leaves are at the
        same depth.

        For example::

            |         ___o___   |
            |        /       \  |
            |     __o__       o |
            |    /     \        |
            |   o       o       |
            |  / \     / \      |
            | o   o   o   o     |

        is not perfect but the next one is::

            |     __o__     |
            |    /     \    |
            |   o       o   |
            |  / \     / \  |
            | o   o   o   o |

        EXAMPLES::

            sage: lst = lambda i: filter(lambda bt: bt.is_perfect(), BinaryTrees(i))
            sage: for i in range(8): ascii_art(lst(i)) # long time
            [  ]
            [ o ]
            [  ]
            [   o   ]
            [  / \  ]
            [ o   o ]
            [  ]
            [  ]
            [  ]
            [     __o__     ]
            [    /     \    ]
            [   o       o   ]
            [  / \     / \  ]
            [ o   o   o   o ]
        """
        return 2 ** self.depth() - 1 == self.node_number()

    def is_complete(self):
        r"""
        Return ``True`` if ``self`` is complete, else return ``False``.

        In a nutshell, a complete binary tree is a perfect binary tree
        except possibly in the last level, with all nodes in the last
        level "flush to the left".

        In more detail:
        A complete binary tree (also called binary heap) is a binary tree in
        which every level, except possibly the last one (the deepest), is
        completely filled. At depth `n`, all nodes must be as far left as
        possible.

        For example::

            |         ___o___   |
            |        /       \  |
            |     __o__       o |
            |    /     \        |
            |   o       o       |
            |  / \     / \      |
            | o   o   o   o     |

        is not complete but the following ones are::

            |     __o__          _o_            ___o___     |
            |    /     \        /   \          /       \    |
            |   o       o      o     o      __o__       o   |
            |  / \     / \    / \          /     \     / \  |
            | o   o   o   o, o   o    ,   o       o   o   o |
            |                            / \     /          |
            |                           o   o   o           |

        EXAMPLES::

            sage: lst = lambda i: filter(lambda bt: bt.is_complete(), BinaryTrees(i))
            sage: for i in range(8): ascii_art(lst(i)) # long time
            [  ]
            [ o ]
            [   o ]
            [  /  ]
            [ o   ]
            [   o   ]
            [  / \  ]
            [ o   o ]
            [     o   ]
            [    / \  ]
            [   o   o ]
            [  /      ]
            [ o       ]
            [     _o_   ]
            [    /   \  ]
            [   o     o ]
            [  / \      ]
            [ o   o     ]
            [     __o__   ]
            [    /     \  ]
            [   o       o ]
            [  / \     /  ]
            [ o   o   o   ]
            [     __o__     ]
            [    /     \    ]
            [   o       o   ]
            [  / \     / \  ]
            [ o   o   o   o ]
        """
        if self.is_empty():
            return True
        # self := L ^ R
        dL = self[0].depth()
        dR = self[1].depth()
        # if L is perfect
        if self[0].is_perfect():
            # if the depth of R == depth of L then R must be complete
            if dL == dR:
                return self[1].is_complete()
            # else R is perfect with depth equals depth of L - 1
            elif dL == dR + 1:
                return self[1].is_perfect()
            return False
        # L is not perfect then R is perfect and the depth of L = the depth of
        # R + 1
        return self[0].is_complete() and self[1].is_perfect() and dL == dR + 1


# Abstract class to serve as a Factory no instance are created.
class BinaryTrees(UniqueRepresentation, Parent):
    """
    Factory for binary trees.

    INPUT:

    - ``size`` -- (optional) an integer

    OUTPUT:

    - the set of all binary trees (of the given ``size`` if specified)

    EXAMPLES::

        sage: BinaryTrees()
        Binary trees

        sage: BinaryTrees(2)
        Binary trees of size 2

    .. NOTE:: this is a factory class whose constructor returns instances of
              subclasses.

    .. NOTE:: the fact that BinaryTrees is a class instead of a simple callable
              is an implementation detail. It could be changed in the future
              and one should not rely on it.
    """
    @staticmethod
    def __classcall_private__(cls, n=None):
        """
        TESTS::

            sage: from sage.combinat.binary_tree import BinaryTrees_all, BinaryTrees_size
            sage: isinstance(BinaryTrees(2), BinaryTrees)
            True
            sage: isinstance(BinaryTrees(), BinaryTrees)
            True
            sage: BinaryTrees(2) is BinaryTrees_size(2)
            True
            sage: BinaryTrees(5).cardinality()
            42
            sage: BinaryTrees() is BinaryTrees_all()
            True
        """
        if n is None:
            return BinaryTrees_all()
        else:
            if not (isinstance(n, (Integer, int)) and n >= 0):
                raise ValueError("n must be a nonnegative integer")
            return BinaryTrees_size(Integer(n))

    @cached_method
    def leaf(self):
        """
        Return a leaf tree with ``self`` as parent.

        EXAMPLES::

            sage: BinaryTrees().leaf()
            .

        TEST::

            sage: (BinaryTrees().leaf() is
            ....:  sage.combinat.binary_tree.BinaryTrees_all().leaf())
            True
        """
        return self(None)

#################################################################
# Enumerated set of all binary trees
#################################################################


class BinaryTrees_all(DisjointUnionEnumeratedSets, BinaryTrees):

    def __init__(self):
        """
        TESTS::

            sage: from sage.combinat.binary_tree import BinaryTrees_all
            sage: B = BinaryTrees_all()
            sage: B.cardinality()
            +Infinity

            sage: it = iter(B)
            sage: (next(it), next(it), next(it), next(it), next(it))
            (., [., .], [., [., .]], [[., .], .], [., [., [., .]]])
            sage: next(it).parent()
            Binary trees
            sage: B([])
            [., .]

            sage: B is BinaryTrees_all()
            True
            sage: TestSuite(B).run() # long time
            """
        DisjointUnionEnumeratedSets.__init__(
            self, Family(NonNegativeIntegers(), BinaryTrees_size),
            facade=True, keepkey = False)

    def _repr_(self):
        """
        TEST::

            sage: BinaryTrees()   # indirect doctest
            Binary trees
        """
        return "Binary trees"

    def __contains__(self, x):
        """
        TESTS::

            sage: S = BinaryTrees()
            sage: 1 in S
            False
            sage: S([]) in S
            True
        """
        return isinstance(x, self.element_class)

    def __call__(self, x=None, *args, **keywords):
        """
        Ensure that ``None`` instead of ``0`` is passed by default.

        TESTS::

            sage: B = BinaryTrees()
            sage: B()
            .
        """
        return super(BinaryTrees, self).__call__(x, *args, **keywords)

    def unlabelled_trees(self):
        """
        Return the set of unlabelled trees associated to ``self``.

        EXAMPLES::

            sage: BinaryTrees().unlabelled_trees()
            Binary trees
        """
        return self

    def labelled_trees(self):
        """
        Return the set of labelled trees associated to ``self``.

        EXAMPLES::

            sage: BinaryTrees().labelled_trees()
            Labelled binary trees
        """
        return LabelledBinaryTrees()

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: B = BinaryTrees()
            sage: B._element_constructor_([])
            [., .]
            sage: B([[],[]]) # indirect doctest
            [[., .], [., .]]
            sage: B(None)    # indirect doctest
            .
        """
        return self.element_class(self, *args, **keywords)

    Element = BinaryTree


#################################################################
# Enumerated set of binary trees of a given size
#################################################################


class BinaryTrees_size(BinaryTrees):
    """
    The enumerated sets of binary trees of given size

    TESTS::

        sage: from sage.combinat.binary_tree import BinaryTrees_size
        sage: for i in range(6): TestSuite(BinaryTrees_size(i)).run()
    """
    def __init__(self, size):
        """
        TESTS::

            sage: S = BinaryTrees(3)
            sage: S == loads(dumps(S))
            True

            sage: S is BinaryTrees(3)
            True
        """
        super(BinaryTrees_size, self).__init__(category = FiniteEnumeratedSets())
        self._size = size

    def _repr_(self):
        """
        TESTS::

            sage: BinaryTrees(3)   # indirect doctest
            Binary trees of size 3
        """
        return "Binary trees of size %s" % (self._size)

    def __contains__(self, x):
        """
        TESTS::

            sage: S = BinaryTrees(3)
            sage: 1 in S
            False
            sage: S([[],[]]) in S
            True
        """
        return isinstance(x, self.element_class) and x.node_number() == self._size

    def _an_element_(self):
        """
        TESTS::

            sage: BinaryTrees(0).an_element()  # indirect doctest
            .
        """
        return self.first()

    def cardinality(self):
        """
        The cardinality of ``self``

        This is a Catalan number.

        TESTS::

            sage: BinaryTrees(0).cardinality()
            1
            sage: BinaryTrees(5).cardinality()
            42
        """
        from .combinat import catalan_number
        return catalan_number(self._size)

    def random_element(self):
        r"""
        Return a random ``BinaryTree`` with uniform probability.

        This method generates a random ``DyckWord`` and then uses a
        bijection between Dyck words and binary trees.

        EXAMPLES::

            sage: BinaryTrees(5).random_element() # random
            [., [., [., [., [., .]]]]]
            sage: BinaryTrees(0).random_element()
            .
            sage: BinaryTrees(1).random_element()
            [., .]

        TESTS::

            sage: all([BinaryTrees(10).random_element() in BinaryTrees(10) for i in range(20)])
            True
        """
        from sage.combinat.dyck_word import CompleteDyckWords_size
        return CompleteDyckWords_size(self._size).random_element().to_binary_tree()

    def __iter__(self):
        """
        A basic generator.

        .. TODO:: could be optimized.

        TESTS::

            sage: BinaryTrees(0).list()
            [.]
            sage: BinaryTrees(1).list()
            [[., .]]
            sage: BinaryTrees(3).list()
            [[., [., [., .]]], [., [[., .], .]], [[., .], [., .]], [[., [., .]], .], [[[., .], .], .]]
        """
        if self._size == 0:
            yield self._element_constructor_()
        else:
            for i in range(0, self._size):
                for lft in self.__class__(i):
                    for rgt in self.__class__(self._size - 1 - i):
                        yield self._element_constructor_([lft, rgt])

    @lazy_attribute
    def _parent_for(self):
        """
        The parent of the elements generated by ``self``.

        TESTS::

            sage: S = BinaryTrees(3)
            sage: S._parent_for
            Binary trees
        """
        return BinaryTrees_all()

    @lazy_attribute
    def element_class(self):
        """
        TESTS::

            sage: S = BinaryTrees(3)
            sage: S.element_class
            <class 'sage.combinat.binary_tree.BinaryTrees_all_with_category.element_class'>
            sage: S.first().__class__ == BinaryTrees().first().__class__
            True
        """
        return self._parent_for.element_class

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: S = BinaryTrees(0)
            sage: S([])   # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: wrong number of nodes
            sage: S(None)   # indirect doctest
            .

            sage: S = BinaryTrees(1)   # indirect doctest
            sage: S([])
            [., .]
        """
        res = self.element_class(self._parent_for, *args, **keywords)
        if res.node_number() != self._size:
            raise ValueError("wrong number of nodes")
        return res


class LabelledBinaryTree(AbstractLabelledClonableTree, BinaryTree):
    """
    Labelled binary trees.

    A labelled binary tree is a binary tree (see :class:`BinaryTree` for
    the meaning of this) with a label assigned to each node.
    The labels need not be integers, nor are they required to be distinct.
    ``None`` can be used as a label.

    .. WARNING::

        While it is possible to assign values to leaves (not just nodes)
        using this class, these labels are disregarded by various
        methods such as
        :meth:`~sage.combinat.abstract_tree.AbstractLabelledTree.labels`,
        :meth:`~sage.combinat.abstract_tree.AbstractLabelledTree.map_labels`,
        and (ironically)
        :meth:`~sage.combinat.abstract_tree.AbstractLabelledTree.leaf_labels`.

    INPUT:

    - ``children`` -- ``None`` (default) or a list, tuple or iterable of
      length `2` of labelled binary trees or convertible objects. This
      corresponds to the standard recursive definition of a labelled
      binary tree as being either a leaf, or a pair of:

      - a pair of labelled binary trees,
      - and a label.

      (The label is specified in the keyword variable ``label``; see
      below.)

      Syntactic sugar allows leaving out all but the outermost calls
      of the ``LabelledBinaryTree()`` constructor, so that, e. g.,
      ``LabelledBinaryTree([LabelledBinaryTree(None),LabelledBinaryTree(None)])``
      can be shortened to ``LabelledBinaryTree([None,None])``. However,
      using this shorthand, it is impossible to label any vertex of
      the tree other than the root (because there is no way to pass a
      ``label`` variable without calling ``LabelledBinaryTree``
      explicitly).

      It is also allowed to abbreviate ``[None, None]`` by ``[]`` if
      one does not want to label the leaves (which one should not do
      anyway!).

    - ``label`` -- (default: ``None``) the label to be put on the root
      of this tree.

    - ``check`` -- (default: ``True``) whether checks should be
      performed or not.

    .. TODO::

        It is currently not possible to use ``LabelledBinaryTree()``
        as a shorthand for ``LabelledBinaryTree(None)`` (in analogy to
        similar syntax in the ``BinaryTree`` class).

    EXAMPLES::

        sage: LabelledBinaryTree(None)
        .
        sage: LabelledBinaryTree(None, label="ae")    # not well supported
        'ae'
        sage: LabelledBinaryTree([])
        None[., .]
        sage: LabelledBinaryTree([], label=3)    # not well supported
        3[., .]
        sage: LabelledBinaryTree([None, None])
        None[., .]
        sage: LabelledBinaryTree([None, None], label=5)
        5[., .]
        sage: LabelledBinaryTree([None, []])
        None[., None[., .]]
        sage: LabelledBinaryTree([None, []], label=4)
        4[., None[., .]]
        sage: LabelledBinaryTree([[], None])
        None[None[., .], .]
        sage: LabelledBinaryTree("[[], .]", label=False)
        False[None[., .], .]
        sage: LabelledBinaryTree([None, LabelledBinaryTree([None, None], label=4)], label=3)
        3[., 4[., .]]
        sage: LabelledBinaryTree([None, BinaryTree([None, None])], label=3)
        3[., None[., .]]

        sage: LabelledBinaryTree([[], None, []])
        Traceback (most recent call last):
        ...
        ValueError: this is not a binary tree

        sage: LBT = LabelledBinaryTree
        sage: t1 = LBT([[LBT([], label=2), None], None], label=4); t1
        4[None[2[., .], .], .]

    TESTS::

        sage: t1 = LabelledBinaryTree([[None, [[],[[], None]]],[[],[]]])
        sage: t2 = LabelledBinaryTree([[[],[]],[]])
        sage: with t1.clone() as t1c:
        ....:     t1c[1,1,1] = t2
        sage: t1 == t1c
        False

    We check for :trac:`16314`::

        sage: t1 = LBT([ LBT([LBT([], label=2),
        ....:                 LBT([], label=5)], label=6),
        ....:            None], label=4); t1
        4[6[2[., .], 5[., .]], .]
        sage: class Foo(LabelledBinaryTree):
        ....:     pass
        sage: t2 = Foo(t1.parent(), t1); t2
        4[6[2[., .], 5[., .]], .]
        sage: t2.label()
        4
        sage: t2[0].label()
        6
        sage: t2.__class__, t2[0].__class__
        (<class '__main__.Foo'>, <class '__main__.Foo'>)
    """
    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        """
        Ensure that trees created by the sets and directly are the same and
        that they are instances of :class:`LabelledTree`.

        TESTS::

            sage: issubclass(LabelledBinaryTrees().element_class, LabelledBinaryTree)
            True
            sage: t0 = LabelledBinaryTree([[],[[], None]], label = 3)
            sage: t0.parent()
            Labelled binary trees
            sage: type(t0)
            <class 'sage.combinat.binary_tree.LabelledBinaryTrees_with_category.element_class'>
        """
        return cls._auto_parent.element_class(cls._auto_parent, *args, **opts)

    @lazy_class_attribute
    def _auto_parent(cls):
        """
        The automatic parent of the elements of this class.

        When calling the constructor of an element of this class, one needs a
        parent. This class attribute specifies which parent is used.

        EXAMPLES::

            sage: LabelledBinaryTree._auto_parent
            Labelled binary trees
            sage: LabelledBinaryTree([], label = 3).parent()
            Labelled binary trees
        """
        return LabelledBinaryTrees()

    def _repr_(self):
        """
        TESTS::

            sage: LBT = LabelledBinaryTree
            sage: t1 = LBT([[LBT([], label=2), None], None], label=4); t1
            4[None[2[., .], .], .]
            sage: LBT([[],[[], None]], label = 3)   # indirect doctest
            3[None[., .], None[None[., .], .]]
        """
        if not self:
            if self._label is not None:
                return repr(self._label)
            else:
                return "."
        else:
            return "%s%s" % (self._label, self[:])

    def binary_search_insert(self, letter):
        """
        Return the result of inserting a letter ``letter`` into the
        right strict binary search tree ``self``.

        INPUT:

        - ``letter`` -- any object comparable with the labels of ``self``

        OUTPUT:

        The right strict binary search tree ``self`` with ``letter``
        inserted into it according to the binary search insertion
        algorithm.

        .. NOTE:: ``self`` is supposed to be a binary search tree.
                  This is not being checked!

        A right strict binary search tree is defined to be a labelled
        binary tree such that for each node `n` with label `x`,
        every descendant of the left child of `n` has a label `\leq x`,
        and every descendant of the right child of `n` has a label
        `> x`. (Here, only nodes count as descendants, and every node
        counts as its own descendant too.) Leaves are assumed to have
        no labels.

        Given a right strict binary search tree `t` and a letter `i`,
        the result of inserting `i` into `t` (denoted `Ins(i, t)` in
        the following) is defined recursively as follows:

        - If `t` is empty, then `Ins(i, t)` is the tree with one node
          only, and this node is labelled with `i`.

        - Otherwise, let `j` be the label of the root of `t`. If
          `i > j`, then `Ins(i, t)` is obtained by replacing the
          right child of `t` by `Ins(i, r)` in `t`, where `r` denotes
          the right child of `t`. If `i \leq j`, then `Ins(i, t)` is
          obtained by replacing the left child of `t` by `Ins(i, l)`
          in `t`, where `l` denotes the left child of `t`.

        See, for example, [HNT05]_ for properties of this algorithm.

        .. WARNING::

            If `t` is nonempty, then inserting `i` into `t` does not
            change the root label of `t`. Hence, as opposed to
            algorithms like Robinson-Schensted-Knuth, binary
            search tree insertion involves no bumping.

        EXAMPLES:

        The example from Fig. 2 of [HNT05]_::

            sage: LBT = LabelledBinaryTree
            sage: x = LBT(None)
            sage: x
            .
            sage: x = x.binary_search_insert("b"); x
            b[., .]
            sage: x = x.binary_search_insert("d"); x
            b[., d[., .]]
            sage: x = x.binary_search_insert("e"); x
            b[., d[., e[., .]]]
            sage: x = x.binary_search_insert("a"); x
            b[a[., .], d[., e[., .]]]
            sage: x = x.binary_search_insert("b"); x
            b[a[., b[., .]], d[., e[., .]]]
            sage: x = x.binary_search_insert("d"); x
            b[a[., b[., .]], d[d[., .], e[., .]]]
            sage: x = x.binary_search_insert("a"); x
            b[a[a[., .], b[., .]], d[d[., .], e[., .]]]
            sage: x = x.binary_search_insert("c"); x
            b[a[a[., .], b[., .]], d[d[c[., .], .], e[., .]]]

        Other examples::

            sage: LBT = LabelledBinaryTree
            sage: LBT(None).binary_search_insert(3)
            3[., .]
            sage: LBT([], label = 1).binary_search_insert(3)
            1[., 3[., .]]
            sage: LBT([], label = 3).binary_search_insert(1)
            3[1[., .], .]
            sage: res = LBT(None)
            sage: for i in [3,1,5,2,4,6]:
            ....:     res = res.binary_search_insert(i)
            sage: res
            3[1[., 2[., .]], 5[4[., .], 6[., .]]]
        """
        LT = self.parent()._element_constructor_
        if not self:
            return LT([], label = letter)
        else:
            if letter <= self.label():
                fils = self[0].binary_search_insert(letter)
                return LT([fils, self[1]], label=self.label())
            else:
                fils = self[1].binary_search_insert(letter)
                return LT([self[0], fils], label=self.label())

    def semistandard_insert(self, letter):
        """
        Return the result of inserting a letter ``letter`` into the
        semistandard tree ``self`` using the bumping algorithm.

        INPUT:

        - ``letter`` -- any object comparable with the labels of ``self``

        OUTPUT:

        The semistandard tree ``self`` with ``letter`` inserted into it
        according to the bumping algorithm.

        .. NOTE:: ``self`` is supposed to be a semistandard tree.
                  This is not being checked!

        A semistandard tree is defined to be a labelled binary tree
        such that for each node `n` with label `x`, every descendant of
        the left child of `n` has a label `> x`, and every descendant
        of the right child of `n` has a label `\geq x`. (Here, only
        nodes count as descendants, and every node counts as its own
        descendant too.) Leaves are assumed to have no labels.

        Given a semistandard tree `t` and a letter `i`, the result of
        inserting `i` into `t` (denoted `Ins(i, t)` in the following)
        is defined recursively as follows:

        - If `t` is empty, then `Ins(i, t)` is the tree with one node
          only, and this node is labelled with `i`.

        - Otherwise, let `j` be the label of the root of `t`. If
          `i \geq j`, then `Ins(i, t)` is obtained by replacing the
          right child of `t` by `Ins(i, r)` in `t`, where `r` denotes
          the right child of `t`. If `i < j`, then `Ins(i, t)` is
          obtained by replacing the label at the root of `t` by `i`,
          and replacing the left child of `t` by `Ins(j, l)`
          in `t`, where `l` denotes the left child of `t`.

        This algorithm is similar to the Robinson-Schensted-Knuth
        insertion algorithm for semistandard Young tableaux.

        AUTHORS:

        - Darij Grinberg (10 Nov 2013).

        EXAMPLES::

            sage: LBT = LabelledBinaryTree
            sage: x = LBT(None)
            sage: x
            .
            sage: x = x.semistandard_insert("b"); x
            b[., .]
            sage: x = x.semistandard_insert("d"); x
            b[., d[., .]]
            sage: x = x.semistandard_insert("e"); x
            b[., d[., e[., .]]]
            sage: x = x.semistandard_insert("a"); x
            a[b[., .], d[., e[., .]]]
            sage: x = x.semistandard_insert("b"); x
            a[b[., .], b[d[., .], e[., .]]]
            sage: x = x.semistandard_insert("d"); x
            a[b[., .], b[d[., .], d[e[., .], .]]]
            sage: x = x.semistandard_insert("a"); x
            a[b[., .], a[b[d[., .], .], d[e[., .], .]]]
            sage: x = x.semistandard_insert("c"); x
            a[b[., .], a[b[d[., .], .], c[d[e[., .], .], .]]]

        Other examples::

            sage: LBT = LabelledBinaryTree
            sage: LBT(None).semistandard_insert(3)
            3[., .]
            sage: LBT([], label = 1).semistandard_insert(3)
            1[., 3[., .]]
            sage: LBT([], label = 3).semistandard_insert(1)
            1[3[., .], .]
            sage: res = LBT(None)
            sage: for i in [3,1,5,2,4,6]:
            ....:     res = res.semistandard_insert(i)
            sage: res
            1[3[., .], 2[5[., .], 4[., 6[., .]]]]
        """
        LT = self.parent()._element_constructor_
        if not self:
            return LT([], label = letter)
        else:
            root_label = self.label()
            if letter < root_label:
                fils = self[0].semistandard_insert(root_label)
                return LT([fils, self[1]], label=letter)
            else:
                fils = self[1].semistandard_insert(letter)
                return LT([self[0], fils], label=root_label)

    def right_rotate(self):
        r"""
        Return the result of right rotation applied to the labelled
        binary tree ``self``.

        Right rotation on labelled binary trees is defined as
        follows: Let `T` be a labelled binary tree such that the
        left child of the root of `T` is a node. Let
        `C` be the right child of the root of `T`, and let `A`
        and `B` be the left and right children of the left child
        of the root of `T`. (Keep in mind that nodes of trees are
        identified with the subtrees consisting of their
        descendants.) Furthermore, let `y` be the label at the
        root of `T`, and `x` be the label at the left child of the
        root of `T`.
        Then, the right rotation of `T` is the labelled binary
        tree in which the root is labelled `x`, the left child of
        the root is `A`, whereas the right child of the root is a
        node labelled `y` whose left and right children are `B`
        and `C`. In pictures::

            |     y                      x     |
            |    / \                    / \    |
            |   x   C -right-rotate->  A   y   |
            |  / \                        / \  |
            | A   B                      B   C |

        Right rotation is the inverse operation to left rotation
        (:meth:`left_rotate`).

        TESTS::

            sage: LB = LabelledBinaryTree
            sage: b = LB([LB([LB([],"A"), LB([],"B")],"x"),LB([],"C")], "y"); b
            y[x[A[., .], B[., .]], C[., .]]
            sage: b.right_rotate()
            x[A[., .], y[B[., .], C[., .]]]
        """
        B = self.parent()._element_constructor_
        s0 = self[0]
        return B([s0[0], B([s0[1], self[1]], self.label())], s0.label())

    def left_rotate(self):
        r"""
        Return the result of left rotation applied to the labelled
        binary tree ``self``.

        Left rotation on labelled binary trees is defined as
        follows: Let `T` be a labelled binary tree such that the
        right child of the root of `T` is a node. Let
        `A` be the left child of the root of `T`, and let `B`
        and `C` be the left and right children of the right child
        of the root of `T`. (Keep in mind that nodes of trees are
        identified with the subtrees consisting of their
        descendants.) Furthermore, let `x` be the label at the
        root of `T`, and `y` be the label at the right child of the
        root of `T`.
        Then, the left rotation of `T` is the labelled binary tree
        in which the root is labelled `y`, the right child of the
        root is `C`, whereas the left child of the root is a node
        labelled `x` whose left and right children are `A` and `B`.
        In pictures::

           |     y                    x     |
           |    / \                  / \    |
           |   x   C <-left-rotate- A   y   |
           |  / \                      / \  |
           | A   B                    B   C |

        Left rotation is the inverse operation to right rotation
        (:meth:`right_rotate`).

        TESTS::

            sage: LB = LabelledBinaryTree
            sage: b = LB([LB([LB([],"A"), LB([],"B")],"x"),LB([],"C")], "y"); b
            y[x[A[., .], B[., .]], C[., .]]
            sage: b == b.right_rotate().left_rotate()
            True
        """
        B = self.parent()._element_constructor_
        s1 = self[1]
        return B([B([self[0], s1[0]], self.label()), s1[1]], s1.label())

    def heap_insert(self, l):
        r"""
        Return the result of inserting a letter ``l`` into the binary
        heap (tree) ``self``.

        A binary heap is a labelled complete binary tree such that for
        each node, the label at the node is greater or equal to the
        label of each of its child nodes. (More precisely, this is
        called a max-heap.)

        For example::

            |     _7_   |
            |    /   \  |
            |   5     6 |
            |  / \      |
            | 3   4     |

        is a binary heap.

        See :wikipedia:`Binary_heap#Insert` for a description of how to
        insert a letter into a binary heap. The result is another binary
        heap.

        INPUT:

        - ``letter`` -- any object comparable with the labels of ``self``

        .. NOTE::

            ``self`` is assumed to be a binary heap (tree). No check is
            performed.

        TESTS::

            sage: h = LabelledBinaryTree(None)
            sage: h = h.heap_insert(3); ascii_art([h])
            [ 3 ]
            sage: h = h.heap_insert(4); ascii_art([h])
            [   4 ]
            [  /  ]
            [ 3   ]
            sage: h = h.heap_insert(6); ascii_art([h])
            [   6   ]
            [  / \  ]
            [ 3   4 ]
            sage: h = h.heap_insert(2); ascii_art([h])
            [     6   ]
            [    / \  ]
            [   3   4 ]
            [  /      ]
            [ 2       ]
            sage: ascii_art([h.heap_insert(5)])
            [     _6_   ]
            [    /   \  ]
            [   5     4 ]
            [  / \      ]
            [ 2   3     ]
        """
        B = self.parent()._element_constructor_
        if self.is_empty():
            return B([], l)

        if self.label() < l:
            label_root = l
            label_insert = self.label()
        else:
            label_root = self.label()
            label_insert = l
        L, R = self
        dL = L.depth()
        dR = R.depth()
        # if depth of L is greater than the depth of R
        if dL > dR:
            # if L is perfect we insert in R
            if L.is_perfect():
                return B([L, R.heap_insert(label_insert)], label_root)
            # we insert in L
            return B([L.heap_insert(label_insert), R], label_root)
        # else ==> dL == dR
        # if R is perfect we have to insert on the leftmost leaf
        if R.is_perfect():
            # ## TODO:: can be optimized...
            return B([L.heap_insert(label_insert), R], label_root)
        # else we insert on the right
        return B([L, R.heap_insert(label_insert)], label_root)

    _UnLabelled = BinaryTree


class LabelledBinaryTrees(LabelledOrderedTrees):
    """
    This is a parent stub to serve as a factory class for trees with various
    labels constraints.
    """
    def _repr_(self):
        """
        TESTS::

            sage: LabelledBinaryTrees()   # indirect doctest
            Labelled binary trees
        """
        return "Labelled binary trees"

    def _an_element_(self):
        """
        Return a labelled binary tree.

        EXAMPLE::

            sage: LabelledBinaryTrees().an_element()   # indirect doctest
            toto[42[3[., .], 3[., .]], 5[None[., .], None[., .]]]
        """
        LT = self._element_constructor_
        t = LT([], label=3)
        t1 = LT([t, t], label=42)
        t2 = LT([[], []], label=5)
        return LT([t1, t2], label="toto")

    def unlabelled_trees(self):
        """
        Return the set of unlabelled trees associated to ``self``.

        EXAMPLES::

            sage: LabelledBinaryTrees().unlabelled_trees()
            Binary trees

        This is used to compute the shape::

            sage: t = LabelledBinaryTrees().an_element().shape(); t
            [[[., .], [., .]], [[., .], [., .]]]
            sage: t.parent()
            Binary trees

        TESTS::

            sage: t = LabelledBinaryTrees().an_element()
            sage: t.canonical_labelling()
            4[2[1[., .], 3[., .]], 6[5[., .], 7[., .]]]
        """
        return BinaryTrees_all()

    def labelled_trees(self):
        """
        Return the set of labelled trees associated to ``self``.

        EXAMPLES::

            sage: LabelledBinaryTrees().labelled_trees()
            Labelled binary trees
        """
        return self

    Element = LabelledBinaryTree



################################################################
# Interface attempt with species...
#
# Kept here for further reference when species will be improved
################################################################
# from sage.combinat.species.library import (
#     CombinatorialSpecies, EmptySetSpecies, SingletonSpecies)

# BT = CombinatorialSpecies()
# F =  EmptySetSpecies()
# N =  SingletonSpecies()
# BT.define(F+N*(BT*BT))
# # b3 = BT.isotypes(range(3))
# # tr = b3.list()[1]

# def BTsp_to_bintrees(bt):
#     """
#     sage: from sage.combinat.binary_tree import BT, BTsp_to_bintrees
#     sage: BTsp_to_bintrees(BT.isotypes(range(5))[0])
#     [., [., [., [., [., .]]]]]
#     sage: def spls(size):
#     ....:     return map(BTsp_to_bintrees, BT.isotypes(range(size)).list())
#     sage: spls(3)
#     [[., [., [., .]]], [., [[., .], .]], [[., .], [., .]], [[., [., .]], .], [[[., .], .], .]]
#     sage: all(spls(i) == BinaryTrees(i).list() for i in range(5))
#     True
#     """
#     if len(bt) == 0:
#         return BinaryTree()
#     elif len(bt) == 2 and len(bt[1]) == 2:
#         return BinaryTree(map(BTsp_to_bintrees, list(bt[1])))
#     else:
#         raise ValueError
