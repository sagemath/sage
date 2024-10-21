# sage_setup: distribution = sagemath-graphs
"""
Ordered Rooted Trees

AUTHORS:

- Florent Hivert (2010-2011): initial revision
- Frédéric Chapoton (2010): contributed some methods
"""
# ****************************************************************************
#       Copyright (C) 2010 Florent Hivert <Florent.Hivert@univ-rouen.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import itertools

from sage.structure.list_clone import ClonableArray, ClonableList
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.lazy_attribute import lazy_class_attribute
from sage.misc.lazy_import import lazy_import
from sage.combinat.abstract_tree import (AbstractClonableTree,
                                         AbstractLabelledClonableTree)
from sage.combinat.combinatorial_map import combinatorial_map
from sage.misc.cachefunc import cached_method
from sage.categories.sets_cat import Sets, EmptySetError
from sage.rings.integer import Integer
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.rings.infinity import Infinity

lazy_import('sage.combinat.dyck_word', 'CompleteDyckWords_size')


class OrderedTree(AbstractClonableTree, ClonableList,
        metaclass=InheritComparisonClasscallMetaclass):
    """
    The class of (ordered rooted) trees.

    An ordered tree is constructed from a node, called the root, on which one
    has grafted a possibly empty list of trees. There is a total order on the
    children of a node which is given by the order of the elements in the
    list. Note that there is no empty ordered tree (so the smallest ordered
    tree consists of just one node).

    INPUT:

    One can create a tree from any list (or more generally iterable) of trees
    or objects convertible to a tree. Alternatively a string is also
    accepted. The syntax is the same as for printing: children are grouped by
    square brackets.

    EXAMPLES::

        sage: x = OrderedTree([])
        sage: x1 = OrderedTree([x,x])
        sage: x2 = OrderedTree([[],[]])
        sage: x1 == x2
        True
        sage: tt1 = OrderedTree([x,x1,x2])
        sage: tt2 = OrderedTree([[], [[], []], x2])
        sage: tt1 == tt2
        True

        sage: OrderedTree([]) == OrderedTree()
        True

    TESTS::

        sage: x1.__hash__() == x2.__hash__()
        True
        sage: tt1.__hash__() == tt2.__hash__()
        True

    Trees are usually immutable. However they inherit from
    :class:`sage.structure.list_clone.ClonableList`, so that they can be
    modified using the clone protocol. Let us now see what this means.

    Trying to modify a non-mutable tree raises an error::

        sage: tt1[1] = tt2
        Traceback (most recent call last):
        ...
        ValueError: object is immutable; please change a copy instead.

    Here is the correct way to do it::

        sage: with tt2.clone() as tt2:
        ....:     tt2[1] = tt1
        sage: tt2
        [[], [[], [[], []], [[], []]], [[], []]]

    It is also possible to append a child to a tree::

        sage: with tt2.clone() as tt3:
        ....:     tt3.append(OrderedTree([]))
        sage: tt3
        [[], [[], [[], []], [[], []]], [[], []], []]

    Or to insert a child in a tree::

        sage: with tt2.clone() as tt3:
        ....:     tt3.insert(2, OrderedTree([]))
        sage: tt3
        [[], [[], [[], []], [[], []]], [], [[], []]]

    We check that ``tt1`` is not modified and that everything is correct with
    respect to equality::

        sage: tt1
        [[], [[], []], [[], []]]
        sage: tt1 == tt2
        False
        sage: tt1.__hash__() == tt2.__hash__()
        False

    TESTS::

        sage: tt1bis = OrderedTree(tt1)
        sage: with tt1.clone() as tt1:
        ....:     tt1[1] = tt1bis
        sage: tt1
        [[], [[], [[], []], [[], []]], [[], []]]
        sage: tt1 == tt2
        True
        sage: tt1.__hash__() == tt2.__hash__()
        True
        sage: len(tt1)
        3
        sage: tt1[2]
        [[], []]
        sage: tt1[3]
        Traceback (most recent call last):
        ...
        IndexError: list index out of range
        sage: tt1[1:2]
        [[[], [[], []], [[], []]]]

    Various tests involving construction, equality and hashing::

        sage: OrderedTree() == OrderedTree()
        True
        sage: t1 = OrderedTree([[],[[]]])
        sage: t2 = OrderedTree([[],[[]]])
        sage: t1 == t2
        True
        sage: t2 = OrderedTree(t1)
        sage: t1 == t2
        True
        sage: t1 = OrderedTree([[],[[]]])
        sage: t2 = OrderedTree([[[]],[]])
        sage: t1 == t2
        False

        sage: t1 = OrderedTree([[],[[]]])
        sage: t2 = OrderedTree([[],[[]]])
        sage: t1.__hash__() == t2.__hash__()
        True
        sage: t2 = OrderedTree([[[]],[]])
        sage: t1.__hash__() == t2.__hash__()
        False
        sage: OrderedTree().__hash__() == OrderedTree([]).__hash__()
        True
        sage: tt1 = OrderedTree([t1,t2,t1])
        sage: tt2 = OrderedTree([t1, [[[]],[]], t1])
        sage: tt1.__hash__() == tt2.__hash__()
        True

    Check that the hash value is correctly updated after modification::

        sage: with tt2.clone() as tt2:
        ....:     tt2[1,1] = tt1
        sage: tt1.__hash__() == tt2.__hash__()
        False
    """
    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        """
        Ensure that trees created by the enumerated sets and directly
        are the same and that they are instances of :class:`OrderedTree`

        TESTS::

            sage: issubclass(OrderedTrees().element_class, OrderedTree)
            True
            sage: t0 = OrderedTree([[],[[], []]])
            sage: t0.parent()
            Ordered trees
            sage: type(t0)
            <class 'sage.combinat.ordered_tree.OrderedTrees_all_with_category.element_class'>

            sage: t1 = OrderedTrees()([[],[[], []]])
            sage: t1.parent() is t0.parent()
            True
            sage: type(t1) is type(t0)
            True

            sage: t1 = OrderedTrees(4)([[],[[]]])
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

            sage: OrderedTree([[],[[]]])._auto_parent
            Ordered trees
            sage: OrderedTree([[],[[]]]).parent()
            Ordered trees

        .. NOTE::

            It is possible to bypass the automatic parent mechanism using::

                sage: t1 = OrderedTree.__new__(OrderedTree, Parent(), [])
                sage: t1.__init__(Parent(), [])
                sage: t1
                []
                sage: t1.parent()
                <sage.structure.parent.Parent object at ...>
        """
        return OrderedTrees_all()

    def __init__(self, parent=None, children=None, check=True):
        """
        TESTS::

            sage: t1 = OrderedTrees(4)([[],[[]]])
            sage: TestSuite(t1).run()
            sage: OrderedTrees()("[]") # indirect doctest
            []
            sage: all(OrderedTree(repr(tr)) == tr for i in range(6) for tr in OrderedTrees(i))
            True
        """
        if children is None:
            children = []
        if isinstance(children, str):
            children = eval(children)
        if (children.__class__ is self.__class__ and
                children.parent() == parent):
            children = list(children)
        else:
            children = [self.__class__(parent, x) for x in children]
        ClonableArray.__init__(self, parent, children, check=check)

    def is_empty(self):
        """
        Return if ``self`` is the empty tree.

        For ordered trees, this always returns ``False``.

        .. NOTE:: this is different from ``bool(t)`` which returns whether
                  ``t`` has some child or not.

        EXAMPLES::

            sage: t = OrderedTrees(4)([[],[[]]])
            sage: t.is_empty()
            False
            sage: bool(t)
            True
        """
        return False

    def _to_binary_tree_rec(self, bijection='left'):
        r"""
        Internal recursive method to obtain a binary tree from an ordered
        tree.

        See :meth:`to_binary_tree_left_branch` and
        :meth:`to_binary_tree_right_branch` for what it does.

        EXAMPLES::

            sage: T = OrderedTree([[],[]])
            sage: T._to_binary_tree_rec()
            [[., .], .]
            sage: T._to_binary_tree_rec(bijection='right')
            [., [., .]]
            sage: T = OrderedTree([[], [[], []], [[], [[]]]])
            sage: T._to_binary_tree_rec()
            [[[., .], [[., .], .]], [[., .], [., .]]]
            sage: T._to_binary_tree_rec(bijection='right')
            [., [[., [., .]], [[., [[., .], .]], .]]]
        """
        from sage.combinat.binary_tree import BinaryTree
        root = BinaryTree()
        if bijection == "left":
            for child in self:
                root = BinaryTree([root, child._to_binary_tree_rec(bijection)])
        elif bijection == "right":
            children = list(self)
            children.reverse()
            for child in children:
                root = BinaryTree([child._to_binary_tree_rec(bijection), root])
        else:
            raise ValueError("the bijection argument should be either "
                             "left or right")
        return root

    @combinatorial_map(name="To binary tree, left brother = left child")
    def to_binary_tree_left_branch(self):
        r"""
        Return a binary tree of size `n-1` (where `n` is the size of `t`,
        and where `t` is ``self``) obtained from `t` by the following
        recursive rule:

        - if `x` is the left brother of `y` in `t`, then `x` becomes the
          left child of `y`;
        - if `x` is the last child of `y` in `t`, then `x` becomes the
          right child of `y`,

        and removing the root of `t`.

        EXAMPLES::

            sage: T = OrderedTree([[],[]])
            sage: T.to_binary_tree_left_branch()
            [[., .], .]
            sage: T = OrderedTree([[], [[], []], [[], [[]]]])
            sage: T.to_binary_tree_left_branch()
            [[[., .], [[., .], .]], [[., .], [., .]]]

        TESTS::

            sage: T = OrderedTree([[],[]])
            sage: T == T.to_binary_tree_left_branch().to_ordered_tree_left_branch()
            True
            sage: T = OrderedTree([[], [[], []], [[], [[]]]])
            sage: T == T.to_binary_tree_left_branch().to_ordered_tree_left_branch()
            True
        """
        return self._to_binary_tree_rec()

    @combinatorial_map(name="To parallelogram polyomino")
    def to_parallelogram_polyomino(self, bijection=None):
        r"""
        Return a polyomino parallelogram.

        INPUT:

        - ``bijection`` -- (default: ``'Boussicault-Socci'``) is the name of the
          bijection to use. Possible values are ``'Boussicault-Socci'``,
          ``'via dyck and Delest-Viennot'``.

        EXAMPLES::

            sage: # needs sage.combinat sage.modules
            sage: T = OrderedTree([[[], [[], [[]]]], [], [[[],[]]], [], []])
            sage: T.to_parallelogram_polyomino(bijection='Boussicault-Socci')
            [[0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1],
             [1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0]]
            sage: T = OrderedTree( [] )
            sage: T.to_parallelogram_polyomino()
            [[1], [1]]
            sage: T = OrderedTree( [[]] )
            sage: T.to_parallelogram_polyomino()
            [[0, 1], [1, 0]]
            sage: T = OrderedTree( [[],[]] )
            sage: T.to_parallelogram_polyomino()
            [[0, 1, 1], [1, 1, 0]]
            sage: T = OrderedTree( [[[]]] )
            sage: T.to_parallelogram_polyomino()
            [[0, 0, 1], [1, 0, 0]]
        """
        if (bijection is None) or (bijection == 'Boussicault-Socci'):
            return self._to_parallelogram_polyomino_Boussicault_Socci()
        if bijection == 'via dyck and Delest-Viennot':
            raise NotImplementedError
        raise ValueError('unknown bijection')

    def _to_parallelogram_polyomino_Boussicault_Socci(self):
        r"""
        Return the polyomino parallelogram using the Boussicault-Socci
        bijection.

        EXAMPLES::

            sage: # needs sage.combinat sage.modules
            sage: T = OrderedTree(
            ....:     [[[], [[], [[]]]], [], [[[],[]]], [], []]
            ....: )
            sage: T._to_parallelogram_polyomino_Boussicault_Socci()
            [[0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1], [1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0]]
            sage: T = OrderedTree( [] )
            sage: T._to_parallelogram_polyomino_Boussicault_Socci()
            [[1], [1]]
            sage: T = OrderedTree( [[]] )
            sage: T._to_parallelogram_polyomino_Boussicault_Socci()
            [[0, 1], [1, 0]]
            sage: T = OrderedTree( [[],[]] )
            sage: T._to_parallelogram_polyomino_Boussicault_Socci()
            [[0, 1, 1], [1, 1, 0]]
            sage: T = OrderedTree( [[[]]] )
            sage: T._to_parallelogram_polyomino_Boussicault_Socci()
            [[0, 0, 1], [1, 0, 0]]
        """
        from sage.combinat.parallelogram_polyomino import ParallelogramPolyomino
        if self.node_number() == 1:
            return ParallelogramPolyomino([[1], [1]])
        upper_nodes = []
        lower_nodes = []
        w_coordinate = {}
        w_coordinate[()] = 0
        h_coordinate = {}

        cpt = 0
        for h in range(0, self.depth(), 2):
            for node in self.paths_at_depth(h):
                h_coordinate[node] = cpt
                lower_nodes.append(node)
                cpt += 1

        cpt = 0
        for h in range(1, self.depth(), 2):
            for node in self.paths_at_depth(h):
                w_coordinate[node] = cpt
                upper_nodes.append(node)
                cpt += 1

        def W(path):
            if path in w_coordinate:
                return w_coordinate[path]
            else:
                return w_coordinate[path[:-1]]

        def H(path):
            if path in h_coordinate:
                return h_coordinate[path]
            else:
                return h_coordinate[path[:-1]]

        lower_path = []
        for i in range(1, len(lower_nodes)):
            lower_path.append(0)
            lower_path += [1] * (W(lower_nodes[i]) - W(lower_nodes[i - 1]))
        lower_path.append(0)
        lower_path += [1] * (self.node_number() - len(lower_path))

        upper_path = []
        for i in range(1, len(upper_nodes)):
            upper_path.append(1)
            upper_path += [0] * (H(upper_nodes[i]) - H(upper_nodes[i - 1]))
        upper_path.append(1)
        upper_path += [0] * (self.node_number() - len(upper_path))

        return ParallelogramPolyomino([lower_path, upper_path])

    @combinatorial_map(name="To binary tree, right brother = right child")
    def to_binary_tree_right_branch(self):
        r"""
        Return a binary tree of size `n-1` (where `n` is the size of `t`,
        and where `t` is ``self``) obtained from `t` by the following
        recursive rule:

        - if `x` is the right brother of `y` in `t`, then`x` becomes the
          right child of `y`;
        - if `x` is the first child of `y` in `t`, then `x` becomes the
          left child of `y`,

        and removing the root of `t`.

        EXAMPLES::

            sage: T = OrderedTree([[],[]])
            sage: T.to_binary_tree_right_branch()
            [., [., .]]
            sage: T = OrderedTree([[], [[], []], [[], [[]]]])
            sage: T.to_binary_tree_right_branch()
            [., [[., [., .]], [[., [[., .], .]], .]]]

        TESTS::

            sage: T = OrderedTree([[],[]])
            sage: T == T.to_binary_tree_right_branch().to_ordered_tree_right_branch()
            True
            sage: T = OrderedTree([[], [[], []], [[], [[]]]])
            sage: T == T.to_binary_tree_right_branch().to_ordered_tree_right_branch()
            True
        """
        return self._to_binary_tree_rec(bijection='right')

    @combinatorial_map(name="To Dyck path")
    def to_dyck_word(self):
        r"""
        Return the Dyck path corresponding to ``self`` where the maximal
        height of the Dyck path is the depth of ``self`` .

        EXAMPLES::

            sage: T = OrderedTree([[],[]])
            sage: T.to_dyck_word()                                                      # needs sage.combinat
            [1, 0, 1, 0]
            sage: T = OrderedTree([[],[[]]])
            sage: T.to_dyck_word()                                                      # needs sage.combinat
            [1, 0, 1, 1, 0, 0]
            sage: T = OrderedTree([[], [[], []], [[], [[]]]])
            sage: T.to_dyck_word()                                                      # needs sage.combinat
            [1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0]
        """
        word = []
        for child in self:
            word.append(1)
            word.extend(child.to_dyck_word())
            word.append(0)
        from sage.combinat.dyck_word import DyckWord
        return DyckWord(word)

    @combinatorial_map(name="To graph")
    def to_undirected_graph(self):
        r"""
        Return the undirected graph obtained from the tree nodes and edges.

        The graph is endowed with an embedding, so that it will be displayed
        correctly.

        EXAMPLES::

            sage: t = OrderedTree([])
            sage: t.to_undirected_graph()
            Graph on 1 vertex
            sage: t = OrderedTree([[[]],[],[]])
            sage: t.to_undirected_graph()
            Graph on 5 vertices

        If the tree is labelled, we use its labelling to label the graph. This
        will fail if the labels are not all distinct.
        Otherwise, we use the graph canonical labelling which means that
        two different trees can have the same graph.

        EXAMPLES::

            sage: t = OrderedTree([[[]],[],[]])
            sage: t.canonical_labelling().to_undirected_graph()
            Graph on 5 vertices

        TESTS::

            sage: t.canonical_labelling().to_undirected_graph() == t.to_undirected_graph()
            False
            sage: OrderedTree([[],[]]).to_undirected_graph() == OrderedTree([[[]]]).to_undirected_graph()
            True
            sage: OrderedTree([[],[],[]]).to_undirected_graph() == OrderedTree([[[[]]]]).to_undirected_graph()
            False
        """
        from sage.graphs.graph import Graph
        g = Graph()
        if self in LabelledOrderedTrees():
            relabel = False
        else:
            self = self.canonical_labelling()
            relabel = True
        roots = [self]
        g.add_vertex(name=self.label())
        emb = {self.label(): []}
        while roots:
            node = roots.pop()
            children = reversed([child.label() for child in node])
            emb[node.label()].extend(children)
            for child in node:
                g.add_vertex(name=child.label())
                emb[child.label()] = [node.label()]
                g.add_edge(child.label(), node.label())
                roots.append(child)
        g.set_embedding(emb)
        if relabel:
            g = g.canonical_label()
        return g

    @combinatorial_map(name="To poset")
    def to_poset(self, root_to_leaf=False):
        r"""
        Return the poset obtained by interpreting the tree as a Hasse
        diagram. The default orientation is from leaves to root but you can
        pass ``root_to_leaf=True`` to obtain the inverse orientation.

        INPUT:

        - ``root_to_leaf`` -- boolean (default: ``False``); ``True`` if the
          poset orientation should be from root to leaves

        EXAMPLES::

            sage: t = OrderedTree([])
            sage: t.to_poset()
            Finite poset containing 1 elements
            sage: p = OrderedTree([[[]],[],[]]).to_poset()
            sage: p.height(), p.width()                                                 # needs networkx
            (3, 3)

        If the tree is labelled, we use its labelling to label the poset.
        Otherwise, we use the poset canonical labelling::

            sage: t = OrderedTree([[[]],[],[]]).canonical_labelling().to_poset()
            sage: t.height(), t.width()                                                 # needs networkx
            (3, 3)
        """
        if self in LabelledOrderedTrees():
            relabel = False
        else:
            self = self.canonical_labelling()
            relabel = True
        relations = []
        elements = [self.label()]
        roots = [self]
        while roots:
            node = roots.pop()
            for child in node:
                elements.append(child.label())
                relations.append((node.label(), child.label())
                                 if root_to_leaf else (child.label(),
                                                       node.label()))
                roots.append(child)
        from sage.combinat.posets.posets import Poset
        p = Poset([elements, relations])
        if relabel:
            p = p.canonical_label()
        return p

    @combinatorial_map(order=2, name="Left-right symmetry")
    def left_right_symmetry(self):
        r"""
        Return the symmetric tree of ``self``.

        The symmetric tree `s(T)` of an ordered tree `T` is
        defined as follows:
        If `T` is an ordered tree with children `C_1, C_2, \ldots, C_k`
        (listed from left to right), then the symmetric tree `s(T)` of
        `T` is the ordered tree with children
        `s(C_k), s(C_{k-1}), \ldots, s(C_1)` (from left to right).

        EXAMPLES::

            sage: T = OrderedTree([[],[[]]])
            sage: T.left_right_symmetry()
            [[[]], []]
            sage: T = OrderedTree([[], [[], []], [[], [[]]]])
            sage: T.left_right_symmetry()
            [[[[]], []], [[], []], []]
        """
        children = [c.left_right_symmetry() for c in self]
        children.reverse()
        return OrderedTree(children)

    def plot(self):
        r"""
        Plot the tree ``self``.

        .. WARNING::

            For a labelled tree, this will fail unless all labels are
            distinct. For unlabelled trees, some arbitrary labels are chosen.
            Use :meth:`_latex_`, ``view``,
            :meth:`_ascii_art_` or ``pretty_print`` for more
            faithful representations of the data of the tree.

        EXAMPLES::

            sage: p = OrderedTree([[[]],[],[]])
            sage: ascii_art(p)
              _o__
             / / /
            o o o
            |
            o
            sage: p.plot()                                                              # needs sage.plot
            Graphics object consisting of 10 graphics primitives

        .. PLOT::

            P = OrderedTree([[[]],[],[]]).plot()
            sphinx_plot(P)

        Now a labelled example::

            sage: g = OrderedTree([[],[[]],[]]).canonical_labelling()
            sage: ascii_art(g)
              _1__
             / / /
            2 3 5
              |
              4
            sage: g.plot()                                                              # needs sage.plot
            Graphics object consisting of 10 graphics primitives

        .. PLOT::

            P = OrderedTree([[],[[]],[]]).canonical_labelling().plot()
            sphinx_plot(P)
        """
        try:
            root = self.label()
            g = self.to_undirected_graph()
        except AttributeError:
            root = 1
            g = self.canonical_labelling().to_undirected_graph()
        return g.plot(layout='tree', tree_root=root, tree_orientation='down')

    def sort_key(self):
        """
        Return a tuple of nonnegative integers encoding the ordered
        tree ``self``.

        The first entry of the tuple is the number of children of the
        root. Then the rest of the tuple is the concatenation of the
        tuples associated to these children (we view the children of
        a tree as trees themselves) from left to right.

        This tuple characterizes the tree uniquely, and can be used to
        sort the ordered trees.

        .. NOTE::

            By default, this method does not encode any extra
            structure that ``self`` might have -- e.g., if you were
            to define a class ``EdgeColoredOrderedTree`` which
            implements edge-colored trees and which inherits from
            :class:`OrderedTree`, then the :meth:`sort_key` method
            it would inherit would forget about the colors of the
            edges (and thus would not characterize edge-colored
            trees uniquely). If you want to preserve extra data,
            you need to override this method or use a new method.
            For instance, on the :class:`LabelledOrderedTree`
            subclass, this method is overridden by a slightly
            different method, which encodes not only the numbers
            of children of the nodes of ``self``, but also their
            labels.
            Be careful with using overridden methods, however:
            If you have (say) a class ``BalancedTree`` which
            inherits from :class:`OrderedTree` and which encodes
            balanced trees, and if you have another class
            ``BalancedLabelledOrderedTree`` which inherits both
            from ``BalancedOrderedTree`` and from
            :class:`LabelledOrderedTree`, then (depending on the MRO)
            the default :meth:`sort_key` method on
            ``BalancedLabelledOrderedTree`` (unless manually
            overridden) will be taken either from ``BalancedTree``
            or from :class:`LabelledOrderedTree`, and in the former
            case will ignore the labelling!

        EXAMPLES::

            sage: RT = OrderedTree
            sage: RT([[],[[]]]).sort_key()
            (2, 0, 1, 0)
            sage: RT([[[]],[]]).sort_key()
            (2, 1, 0, 0)
        """
        l = len(self)
        if l == 0:
            return (0,)
        resu = [l] + [u for t in self for u in t.sort_key()]
        return tuple(resu)

    @cached_method
    def normalize(self, inplace=False):
        r"""
        Return the normalized tree of ``self``.

        INPUT:

        - ``inplace`` -- boolean (default: ``False``); if ``True``,
          then ``self`` is modified and nothing returned. Otherwise
          the normalized tree is returned.

        The normalization of an ordered tree `t` is an ordered tree `s`
        which has the property that `t` and `s` are isomorphic as
        *unordered* rooted trees, and that if two ordered trees `t` and
        `t'` are isomorphic as *unordered* rooted trees, then the
        normalizations of `t` and `t'` are identical. In other words,
        normalization is a map from the set of ordered trees to itself
        which picks a representative from every equivalence class with
        respect to the relation of "being isomorphic as unordered
        trees", and maps every ordered tree to the representative
        chosen from its class.

        This map proceeds recursively by first normalizing every
        subtree, and then sorting the subtrees according to the value
        of the :meth:`sort_key` method.

        Consider the quotient map `\pi` that sends a planar rooted tree to
        the associated unordered rooted tree. Normalization is the
        composite `s \circ \pi`, where `s` is a section of `\pi`.

        EXAMPLES::

            sage: OT = OrderedTree
            sage: ta = OT([[],[[]]])
            sage: tb = OT([[[]],[]])
            sage: ta.normalize() == tb.normalize()
            True
            sage: ta == tb
            False

        An example with inplace normalization::

            sage: OT = OrderedTree
            sage: ta = OT([[],[[]]])
            sage: tb = OT([[[]],[]])
            sage: ta.normalize(inplace=True); ta
            [[], [[]]]
            sage: tb.normalize(inplace=True); tb
            [[], [[]]]
        """
        if not inplace:
            with self.clone() as res:
                resl = res._get_list()
                for i in range(len(resl)):
                    resl[i] = resl[i].normalize()
                resl.sort(key=lambda t: t.sort_key())
            return res
        else:
            resl = self._get_list()
            for i in range(len(resl)):
                resl[i] = resl[i].normalize()
            resl.sort(key=lambda t: t.sort_key())


# Abstract class to serve as a Factory no instance are created.
class OrderedTrees(UniqueRepresentation, Parent):
    """
    Factory for ordered trees.

    INPUT:

    - ``size`` -- integer (optional)

    OUTPUT:

    - the set of all ordered trees (of the given ``size`` if specified)

    EXAMPLES::

        sage: OrderedTrees()
        Ordered trees

        sage: OrderedTrees(2)
        Ordered trees of size 2

    .. NOTE:: this is a factory class whose constructor returns instances of
              subclasses.

    .. NOTE:: the fact that OrderedTrees is a class instead of a simple callable
              is an implementation detail. It could be changed in the future
              and one should not rely on it.
    """
    @staticmethod
    def __classcall_private__(cls, n=None):
        """
        TESTS::

            sage: from sage.combinat.ordered_tree import OrderedTrees_all, OrderedTrees_size
            sage: isinstance(OrderedTrees(2), OrderedTrees)
            True
            sage: isinstance(OrderedTrees(), OrderedTrees)
            True
            sage: OrderedTrees(2) is OrderedTrees_size(2)
            True
            sage: OrderedTrees(5).cardinality()
            14
            sage: OrderedTrees() is OrderedTrees_all()
            True
        """
        if n is None:
            return OrderedTrees_all()
        else:
            if not (isinstance(n, (Integer, int)) and n >= 0):
                raise ValueError("n must be a nonnegative integer")
            return OrderedTrees_size(Integer(n))

    @cached_method
    def leaf(self):
        """
        Return a leaf tree with ``self`` as parent.

        EXAMPLES::

            sage: OrderedTrees().leaf()
            []

        TESTS::

            sage: (OrderedTrees().leaf() is
            ....:     sage.combinat.ordered_tree.OrderedTrees_all().leaf())
            True
        """
        return self([])


class OrderedTrees_all(DisjointUnionEnumeratedSets, OrderedTrees):
    """
    The set of all ordered trees.

    EXAMPLES::

        sage: OT = OrderedTrees(); OT
        Ordered trees
        sage: OT.cardinality()
        +Infinity
    """

    def __init__(self):
        """
        TESTS::

            sage: from sage.combinat.ordered_tree import OrderedTrees_all
            sage: B = OrderedTrees_all()
            sage: B.cardinality()
            +Infinity

            sage: it = iter(B)
            sage: (next(it), next(it), next(it), next(it), next(it))
            ([], [[]], [[], []], [[[]]], [[], [], []])
            sage: next(it).parent()
            Ordered trees
            sage: B([])
            []

            sage: B is OrderedTrees_all()
            True
            sage: TestSuite(B).run() # long time
        """
        DisjointUnionEnumeratedSets.__init__(
            self, Family(NonNegativeIntegers(), OrderedTrees_size),
            facade=True, keepkey=False)

    def _repr_(self):
        """
        TESTS::

            sage: OrderedTrees()   # indirect doctest
            Ordered trees
        """
        return "Ordered trees"

    def __contains__(self, x):
        """
        TESTS::

            sage: T = OrderedTrees()
            sage: 1 in T
            False
            sage: T([]) in T
            True
        """
        return isinstance(x, self.element_class)

    def unlabelled_trees(self):
        """
        Return the set of unlabelled trees associated to ``self``.

        EXAMPLES::

            sage: OrderedTrees().unlabelled_trees()
            Ordered trees
        """
        return self

    def labelled_trees(self):
        """
        Return the set of labelled trees associated to ``self``.

        EXAMPLES::

            sage: OrderedTrees().labelled_trees()
            Labelled ordered trees
        """
        return LabelledOrderedTrees()

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: T = OrderedTrees()
            sage: T([])     # indirect doctest
            []
        """
        return self.element_class(self, *args, **keywords)

    Element = OrderedTree


from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.composition import Compositions
#################################################################
# Enumerated set of binary trees of a given size
#################################################################


class OrderedTrees_size(OrderedTrees):
    """
    The enumerated sets of binary trees of a given size.

    EXAMPLES::

        sage: S = OrderedTrees(3); S
        Ordered trees of size 3
        sage: S.cardinality()
        2
        sage: S.list()
        [[[], []], [[[]]]]
    """

    def __init__(self, size):
        """
        TESTS::

            sage: from sage.combinat.ordered_tree import OrderedTrees_size
            sage: TestSuite(OrderedTrees_size(0)).run()
            sage: for i in range(6): TestSuite(OrderedTrees_size(i)).run()
        """
        super().__init__(category=FiniteEnumeratedSets())
        self._size = size

    def _repr_(self):
        """
        TESTS::

            sage: OrderedTrees(3)   # indirect doctest
            Ordered trees of size 3
        """
        return "Ordered trees of size {}".format(self._size)

    def __contains__(self, x):
        """
        TESTS::

            sage: T = OrderedTrees(3)
            sage: 1 in T
            False
            sage: T([[],[]]) in T
            True
        """
        return isinstance(x, self.element_class) and x.node_number() == self._size

    def _an_element_(self):
        """
        TESTS::

            sage: OrderedTrees(3).an_element()   # indirect doctest
            [[], []]
        """
        if self._size == 0:
            raise EmptySetError
        return self.first()

    def cardinality(self):
        """
        The cardinality of ``self``.

        This is a Catalan number.

        TESTS::

            sage: OrderedTrees(0).cardinality()
            0
            sage: OrderedTrees(1).cardinality()
            1
            sage: OrderedTrees(6).cardinality()
            42
        """
        if self._size == 0:
            return Integer(0)
        else:
            from .combinat import catalan_number
            return catalan_number(self._size - 1)

    def random_element(self):
        """
        Return a random ``OrderedTree`` with uniform probability.

        This method generates a random ``DyckWord`` and then uses a
        bijection between Dyck words and ordered trees.

        EXAMPLES::

            sage: OrderedTrees(5).random_element()  # random                            # needs sage.combinat
            [[[], []], []]
            sage: OrderedTrees(0).random_element()
            Traceback (most recent call last):
            ...
            EmptySetError: there are no ordered trees of size 0
            sage: OrderedTrees(1).random_element()                                      # needs sage.combinat
            []

        TESTS::

            sage: all(OrderedTrees(10).random_element() in OrderedTrees(10)             # needs sage.combinat
            ....:     for i in range(20))
            True
        """
        if self._size == 0:
            raise EmptySetError("there are no ordered trees of size 0")
        return CompleteDyckWords_size(self._size - 1).random_element().to_ordered_tree()

    def __iter__(self):
        """
        A basic generator.

        .. TODO:: could be optimized.

        TESTS::

            sage: OrderedTrees(0).list()
            []
            sage: OrderedTrees(1).list()
            [[]]
            sage: OrderedTrees(2).list()
            [[[]]]
            sage: OrderedTrees(3).list()
            [[[], []], [[[]]]]
            sage: OrderedTrees(4).list()
            [[[], [], []], [[], [[]]], [[[]], []], [[[], []]], [[[[]]]]]
        """
        if self._size == 0:
            return
        for c in Compositions(self._size - 1):
            for lst in itertools.product(*[self.__class__(i) for i in c]):
                yield self._element_constructor_(lst)

    @lazy_attribute
    def _parent_for(self):
        """
        Return the parent of the element generated by ``self``.

        TESTS::

            sage: OrderedTrees(3)._parent_for
            Ordered trees
        """
        return OrderedTrees_all()

    @lazy_attribute
    def element_class(self):
        """
        The class of the element of ``self``.

        EXAMPLES::

            sage: from sage.combinat.ordered_tree import OrderedTrees_size, OrderedTrees_all
            sage: S = OrderedTrees_size(3)
            sage: S.element_class is OrderedTrees().element_class
            True
            sage: S.first().__class__ == OrderedTrees_all().first().__class__
            True
        """
        return self._parent_for.element_class

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: S = OrderedTrees(0)
            sage: S([])   # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: wrong number of nodes

            sage: S = OrderedTrees(1)   # indirect doctest
            sage: S([])
            []
        """
        res = self.element_class(self._parent_for, *args, **keywords)
        if res.node_number() != self._size:
            raise ValueError("wrong number of nodes")
        return res


class LabelledOrderedTree(AbstractLabelledClonableTree, OrderedTree):
    """
    Labelled ordered trees.

    A labelled ordered tree is an ordered tree with a label attached at each
    node.

    INPUT:

    - ``children`` -- list or tuple or more generally any iterable
      of trees or object convertible to trees
    - ``label`` -- any Sage object (default: ``None``)

    EXAMPLES::

        sage: x = LabelledOrderedTree([], label = 3); x
        3[]
        sage: LabelledOrderedTree([x, x, x], label = 2)
        2[3[], 3[], 3[]]
        sage: LabelledOrderedTree((x, x, x), label = 2)
        2[3[], 3[], 3[]]
        sage: LabelledOrderedTree([[],[[], []]], label = 3)
        3[None[], None[None[], None[]]]
    """
    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        """
        Ensure that trees created by the sets and directly are the same and
        that they are instances of :class:`LabelledOrderedTree`

        TESTS::

            sage: issubclass(LabelledOrderedTrees().element_class, LabelledOrderedTree)
            True
            sage: t0 = LabelledOrderedTree([[],[[], []]], label = 3)
            sage: t0.parent()
            Labelled ordered trees
            sage: type(t0)
            <class 'sage.combinat.ordered_tree.LabelledOrderedTrees_with_category.element_class'>
        """
        return cls._auto_parent.element_class(cls._auto_parent, *args, **opts)

    @lazy_class_attribute
    def _auto_parent(cls):
        """
        The automatic parent of the elements of this class.

        When calling the constructor of an element of this class, one needs a
        parent. This class attribute specifies which parent is used.

        EXAMPLES::

            sage: LabelledOrderedTree._auto_parent
            Labelled ordered trees
            sage: LabelledOrderedTree([], label = 3).parent()
            Labelled ordered trees
        """
        return LabelledOrderedTrees()

    _UnLabelled = OrderedTree

    __hash__ = ClonableArray.__hash__

    @combinatorial_map(order=2, name="Left-right symmetry")
    def left_right_symmetry(self):
        r"""
        Return the symmetric tree of ``self``.

        The symmetric tree `s(T)` of a labelled ordered tree `T` is
        defined as follows:
        If `T` is a labelled ordered tree with children
        `C_1, C_2, \ldots, C_k` (listed from left to right), then the
        symmetric tree `s(T)` of `T` is a labelled ordered tree with
        children `s(C_k), s(C_{k-1}), \ldots, s(C_1)` (from left to
        right), and with the same root label as `T`.

        .. NOTE::

            If you have a subclass of :meth:`LabelledOrderedTree`
            which also inherits from another subclass of
            :meth:`OrderedTree` which does not come with a labelling,
            then (depending on the method resolution order) it might
            happen that this method gets overridden by an
            implementation from that other subclass, and thus forgets
            about the labels. In this case you need to manually
            override this method on your subclass.

        EXAMPLES::

            sage: L2 = LabelledOrderedTree([], label=2)
            sage: L3 = LabelledOrderedTree([], label=3)
            sage: T23 = LabelledOrderedTree([L2, L3], label=4)
            sage: T23.left_right_symmetry()
            4[3[], 2[]]
            sage: T223 = LabelledOrderedTree([L2, T23], label=17)
            sage: T223.left_right_symmetry()
            17[4[3[], 2[]], 2[]]
            sage: T223.left_right_symmetry().left_right_symmetry() == T223
            True
        """
        children = [c.left_right_symmetry() for c in self]
        children.reverse()
        return LabelledOrderedTree(children, label=self.label())

    def sort_key(self):
        """
        Return a tuple of nonnegative integers encoding the labelled
        tree ``self``.

        The first entry of the tuple is a pair consisting of the
        number of children of the root and the label of the root. Then
        the rest of the tuple is the concatenation of the tuples
        associated to these children (we view the children of
        a tree as trees themselves) from left to right.

        This tuple characterizes the labelled tree uniquely, and can
        be used to sort the labelled ordered trees provided that the
        labels belong to a type which is totally ordered.

        .. WARNING::

            This method overrides :meth:`OrderedTree.sort_key`
            and returns a result different from what the latter
            would return, as it wants to encode the whole labelled
            tree including its labelling rather than just the
            unlabelled tree. Therefore, be careful with using this
            method on subclasses of :class:`LabelledOrderedTree`;
            under some circumstances they could inherit it from
            another superclass instead of from :class:`OrderedTree`,
            which would cause the method to forget the labelling.
            See the docstring of :meth:`OrderedTree.sort_key`.

        EXAMPLES::

            sage: L2 = LabelledOrderedTree([], label=2)
            sage: L3 = LabelledOrderedTree([], label=3)
            sage: T23 = LabelledOrderedTree([L2, L3], label=4)
            sage: T23.sort_key()
            ((2, 4), (0, 2), (0, 3))
            sage: T32 = LabelledOrderedTree([L3, L2], label=5)
            sage: T32.sort_key()
            ((2, 5), (0, 3), (0, 2))
            sage: T23322 = LabelledOrderedTree([T23, T32, L2], label=14)
            sage: T23322.sort_key()
            ((3, 14), (2, 4), (0, 2), (0, 3), (2, 5), (0, 3), (0, 2), (0, 2))
        """
        l = len(self)
        if l == 0:
            return ((0, self.label()),)
        resu = [(l, self.label())] + [u for t in self for u in t.sort_key()]
        return tuple(resu)


class LabelledOrderedTrees(UniqueRepresentation, Parent):
    """
    This is a parent stub to serve as a factory class for trees with various
    label constraints.

    EXAMPLES::

        sage: LOT = LabelledOrderedTrees(); LOT
        Labelled ordered trees
        sage: x = LOT([], label = 3); x
        3[]
        sage: x.parent() is LOT
        True
        sage: y = LOT([x, x, x], label = 2); y
        2[3[], 3[], 3[]]
        sage: y.parent() is LOT
        True
    """

    def __init__(self, category=None):
        """
        TESTS::

            sage: TestSuite(LabelledOrderedTrees()).run()
        """
        if category is None:
            category = Sets()
        Parent.__init__(self, category=category)

    def _repr_(self):
        """
        TESTS::

            sage: LabelledOrderedTrees()   # indirect doctest
            Labelled ordered trees
        """
        return "Labelled ordered trees"

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: LabelledOrderedTrees().cardinality()
            +Infinity
        """
        return Infinity

    def _an_element_(self):
        """
        Return a labelled ordered tree.

        EXAMPLES::

            sage: LabelledOrderedTrees().an_element()   # indirect doctest
            toto[3[], 42[3[], 3[]], 5[None[]]]
        """
        LT = self._element_constructor_
        t = LT([], label=3)
        t1 = LT([t, t], label=42)
        t2 = LT([[]], label=5)
        return LT([t, t1, t2], label='toto')

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: T = LabelledOrderedTrees()
            sage: T([], label=2)     # indirect doctest
            2[]
        """
        return self.element_class(self, *args, **keywords)

    def unlabelled_trees(self):
        """
        Return the set of unlabelled trees associated to ``self``.

        This is the set of ordered trees, since ``self`` is the set of
        labelled ordered trees.

        EXAMPLES::

            sage: LabelledOrderedTrees().unlabelled_trees()
            Ordered trees
        """
        return OrderedTrees_all()

    def labelled_trees(self):
        """
        Return the set of labelled trees associated to ``self``.

        This is precisely ``self``, because ``self`` already is the set
        of labelled ordered trees.

        EXAMPLES::

            sage: LabelledOrderedTrees().labelled_trees()
            Labelled ordered trees
            sage: LOT = LabelledOrderedTrees()
            sage: x = LOT([], label = 3)
            sage: y = LOT([x, x, x], label = 2)
            sage: y.canonical_labelling()
            1[2[], 3[], 4[]]
        """
        return self

    Element = LabelledOrderedTree
