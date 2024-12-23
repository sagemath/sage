r"""
Abstract Recursive Trees

The purpose of this class is to help implement trees with a specific structure
on the children of each node. For instance, one could want to define a tree in
which each node sees its children as linearly (see the :mod:`Ordered Trees
<sage.combinat.ordered_tree>` module) or cyclically ordered.

**Tree structures**

Conceptually, one can define a tree structure from any object that can contain
others. Indeed, a list can contain lists which contain lists which contain
lists, and thus define a tree ... The same can be done with sets, or any kind
of iterable objects.

While any iterable is sufficient to encode trees, it can prove useful to have
other methods available like isomorphism tests (see next section), conversions
to DiGraphs objects (see :meth:`~.AbstractLabelledTree.as_digraph`) or
computation of the number of automorphisms constrained by the structure on
children. Providing such methods is the whole purpose of the
:class:`AbstractTree` class.

As a result, the :class:`AbstractTree` class is not meant to be
instantiated, but extended. It is expected that classes extending this one may
also inherit from classes representing iterables, for instance
:class:`~sage.structure.list_clone.ClonableArray` or :class:`~sage.structure.list_clone.ClonableList`

**Constrained Trees**

The tree built from a specific container will reflect the properties of the
container. Indeed, if ``A`` is an iterable class whose elements are linearly
ordered, a class ``B`` extending both of :class:`AbstractTree` and ``A`` will
be such that the children of a node will be linearly ordered. If ``A`` behaves
like a set (i.e. if there is no order on the elements it contains), then two
trees will be considered as equal if one can be obtained from the other
through permutations between the children of a same node (see next section).

**Paths and ID**

It is expected that each element of a set of children should be identified by
its index in the container. This way, any node of the tree can be identified
by a word describing a path from the root node.

**Canonical labellings**

Equality between instances of classes extending both :class:`AbstractTree`
and ``A`` is entirely defined by the equality defined on the elements of
``A``. A canonical labelling of such a tree, however, should be such that
two trees ``a`` and ``b`` satisfying ``a == b`` have the same canonical
labellings. On the other hand, the canonical labellings of trees ``a`` and
``b`` satisfying ``a != b`` are expected to be different.

For this reason, the values returned by the :meth:`canonical_labelling
<AbstractTree.canonical_labelling>` method heavily
depend on the data structure used for a node's children and **should be**
**overridden** by most of the classes extending :class:`AbstractTree` if it is
incoherent with the data structure.

**Authors**

- Florent Hivert (2010-2011): initial revision
- Frédéric Chapoton (2011): contributed some methods
"""
import itertools

from sage.structure.list_clone import ClonableArray
from sage.rings.integer import Integer
from sage.misc.misc_c import prod

# Unfortunately Cython forbids multiple inheritance. Therefore, we do not
# inherit from SageObject to be able to inherit from Element or a subclass
# of it later.


class AbstractTree:
    """
    Abstract Tree.

    There is no data structure defined here, as this class is meant to be
    extended, not instantiated.

    .. rubric:: How should this class be extended?

    A class extending :class:`AbstractTree
    <sage.combinat.abstract_tree.AbstractTree>` should respect several
    assumptions:

    * For a tree ``T``, the call ``iter(T)`` should return an iterator on the
      children of the root ``T``.

    * The :meth:`canonical_labelling
      <AbstractTree.canonical_labelling>` method
      should return the same value for trees that are considered equal (see the
      "canonical labellings" section in the documentation of the
      :class:`AbstractTree <sage.combinat.abstract_tree.AbstractTree>` class).

    * For a tree ``T`` the call ``T.parent().labelled_trees()`` should return
      a parent for labelled trees of the same kind: for example,

      - if ``T`` is a binary tree, ``T.parent()`` is ``BinaryTrees()`` and
        ``T.parent().labelled_trees()`` is ``LabelledBinaryTrees()``

      - if ``T`` is an ordered tree, ``T.parent()`` is ``OrderedTrees()`` and
        ``T.parent().labelled_trees()`` is ``LabelledOrderedTrees()``

    TESTS::

        sage: TestSuite(OrderedTree()).run()
        sage: TestSuite(BinaryTree()).run()
    """

    def pre_order_traversal_iter(self):
        r"""
        The depth-first pre-order traversal iterator.

        This method iters each node following the depth-first pre-order
        traversal algorithm (recursive implementation). The algorithm
        is::

            yield the root (in the case of binary trees, if it is not
                a leaf);
            then explore each subtree (by the algorithm) from the
                leftmost one to the rightmost one.

        EXAMPLES:

        For example, on the following binary tree `b`::

            |   ___3____      |
            |  /        \     |
            | 1         _7_   |
            |  \       /   \  |
            |   2     5     8 |
            |        / \      |
            |       4   6     |

        (only the nodes shown), the depth-first pre-order traversal
        algorithm explores `b` in the following order of nodes:
        `3,1,2,7,5,4,6,8`.

        Another example::

            |     __1____ |
            |    /  /   / |
            |   2  6   8_ |
            |   |  |  / / |
            |   3_ 7 9 10 |
            |  / /        |
            | 4 5         |

        The algorithm explores this labelled tree in the following
        order: `1,2,3,4,5,6,7,8,9,10`.

        TESTS::

            sage: b = BinaryTree([[None,[]],[[[],[]],[]]]).canonical_labelling()
            sage: ascii_art([b])
            [   ___3____      ]
            [  /        \     ]
            [ 1         _7_   ]
            [  \       /   \  ]
            [   2     5     8 ]
            [        / \      ]
            [       4   6     ]
            sage: [n.label() for n in b.pre_order_traversal_iter()]
            [3, 1, 2, 7, 5, 4, 6, 8]

            sage: t = OrderedTree([[[[],[]]],[[]],[[],[]]]).canonical_labelling()
            sage: ascii_art([t])
            [     __1____   ]
            [    /  /   /   ]
            [   2  6   8_   ]
            [   |  |  / /   ]
            [   3_ 7 9 10   ]
            [  / /          ]
            [ 4 5           ]
            sage: [n.label() for n in t.pre_order_traversal_iter()]
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

            sage: [n for n in BinaryTree(None).pre_order_traversal_iter()]
            []

        The following test checks that things do not go wrong if some among
        the descendants of the tree are equal or even identical::

            sage: u = BinaryTree(None)
            sage: v = BinaryTree([u, u])
            sage: w = BinaryTree([v, v])
            sage: t = BinaryTree([w, w])
            sage: t.node_number()
            7
            sage: l = [1 for i in t.pre_order_traversal_iter()]
            sage: len(l)
            7
        """
        if self.is_empty():
            return
        yield self
        yield from itertools.chain(*[c.pre_order_traversal_iter()
                                     for c in self])

    def iterative_pre_order_traversal(self, action=None):
        r"""
        Run the depth-first pre-order traversal algorithm (iterative
        implementation) and subject every node encountered
        to some procedure ``action``. The algorithm is::

            manipulate the root with function `action` (in the case
                of a binary tree, only if the root is not a leaf);
            then explore each subtree (by the algorithm) from the
                leftmost one to the rightmost one.

        INPUT:

        - ``action`` -- (optional) a function which takes a node as
          input, and does something during the exploration

        OUTPUT:

        ``None``. (This is *not* an iterator.)

        .. SEEALSO::

            - :meth:`~sage.combinat.abstract_tree.AbstractTree.pre_order_traversal_iter()`
            - :meth:`~sage.combinat.abstract_tree.AbstractTree.pre_order_traversal()`

        TESTS::

            sage: l = []
            sage: b = BinaryTree([[None,[]],[[[],[]],[]]]).canonical_labelling()
            sage: b
            3[1[., 2[., .]], 7[5[4[., .], 6[., .]], 8[., .]]]
            sage: b.iterative_pre_order_traversal(lambda node: l.append(node.label()))
            sage: l
            [3, 1, 2, 7, 5, 4, 6, 8]

            sage: t = OrderedTree([[[[],[]]],[[]],[[],[]]]).canonical_labelling()
            sage: t
            1[2[3[4[], 5[]]], 6[7[]], 8[9[], 10[]]]
            sage: l = []
            sage: t.iterative_pre_order_traversal(lambda node: l.append(node.label()))
            sage: l
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            sage: l = []

            sage: BinaryTree().canonical_labelling().\
            ....:     pre_order_traversal(lambda node: l.append(node.label()))
            sage: l
            []
            sage: OrderedTree([]).canonical_labelling().\
            ....:     iterative_pre_order_traversal(lambda node: l.append(node.label()))
            sage: l
            [1]

        The following test checks that things do not go wrong if some among
        the descendants of the tree are equal or even identical::

            sage: u = BinaryTree(None)
            sage: v = BinaryTree([u, u])
            sage: w = BinaryTree([v, v])
            sage: t = BinaryTree([w, w])
            sage: t.node_number()
            7
            sage: l = []
            sage: t.iterative_pre_order_traversal(lambda node: l.append(1))
            sage: len(l)
            7
        """
        if self.is_empty():
            return
        if action is None:
            def action(x):
                return None
        stack = []
        stack.append(self)
        while stack:
            node = stack.pop()
            action(node)
            for i in range(len(node)):
                subtree = node[-i - 1]
                if not subtree.is_empty():
                    stack.append(subtree)

    def pre_order_traversal(self, action=None):
        r"""
        Run the depth-first pre-order traversal algorithm (recursive
        implementation) and subject every node encountered
        to some procedure ``action``. The algorithm is::

            manipulate the root with function `action` (in the case
                of a binary tree, only if the root is not a leaf);
            then explore each subtree (by the algorithm) from the
                leftmost one to the rightmost one.

        INPUT:

        - ``action`` -- (optional) a function which takes a node as
          input, and does something during the exploration

        OUTPUT:

        ``None``. (This is *not* an iterator.)

        EXAMPLES:

        For example, on the following binary tree `b`::

            |   ___3____      |
            |  /        \     |
            | 1         _7_   |
            |  \       /   \  |
            |   2     5     8 |
            |        / \      |
            |       4   6     |

        the depth-first pre-order traversal algorithm explores `b` in the
        following order of nodes: `3,1,2,7,5,4,6,8`.

        Another example::

            |     __1____ |
            |    /  /   / |
            |   2  6   8_ |
            |   |  |  / / |
            |   3_ 7 9 10 |
            |  / /        |
            | 4 5         |

        The algorithm explores this tree in the following order:
        `1,2,3,4,5,6,7,8,9,10`.

        .. SEEALSO::

            - :meth:`~sage.combinat.abstract_tree.AbstractTree.pre_order_traversal_iter()`
            - :meth:`~sage.combinat.abstract_tree.AbstractTree.iterative_pre_order_traversal()`

        TESTS::

            sage: l = []
            sage: b = BinaryTree([[None,[]],[[[],[]],[]]]).canonical_labelling()
            sage: b
            3[1[., 2[., .]], 7[5[4[., .], 6[., .]], 8[., .]]]
            sage: b.pre_order_traversal(lambda node: l.append(node.label()))
            sage: l
            [3, 1, 2, 7, 5, 4, 6, 8]
            sage: li = []
            sage: b.iterative_pre_order_traversal(lambda node: li.append(node.label()))
            sage: l == li
            True

            sage: t = OrderedTree([[[[],[]]],[[]],[[],[]]]).canonical_labelling()
            sage: t
            1[2[3[4[], 5[]]], 6[7[]], 8[9[], 10[]]]
            sage: l = []
            sage: t.pre_order_traversal(lambda node: l.append(node.label()))
            sage: l
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            sage: li = []
            sage: t.iterative_pre_order_traversal(lambda node: li.append(node.label()))
            sage: l == li
            True

            sage: l = []
            sage: BinaryTree().canonical_labelling().\
            ....:     pre_order_traversal(lambda node: l.append(node.label()))
            sage: l
            []
            sage: OrderedTree([]).canonical_labelling().\
            ....:     pre_order_traversal(lambda node: l.append(node.label()))
            sage: l
            [1]

        The following test checks that things do not go wrong if some among
        the descendants of the tree are equal or even identical::

            sage: u = BinaryTree(None)
            sage: v = BinaryTree([u, u])
            sage: w = BinaryTree([v, v])
            sage: t = BinaryTree([w, w])
            sage: t.node_number()
            7
            sage: l = []
            sage: t.pre_order_traversal(lambda node: l.append(1))
            sage: len(l)
            7
        """
        if action is None:
            def action(x):
                return None
        for node in self.pre_order_traversal_iter():
            action(node)

    def post_order_traversal_iter(self):
        r"""
        The depth-first post-order traversal iterator.

        This method iters each node following the depth-first post-order
        traversal algorithm (recursive implementation). The algorithm
        is::

            explore each subtree (by the algorithm) from the
                leftmost one to the rightmost one;
            then yield the root (in the case of binary trees, only if
                it is not a leaf).

        EXAMPLES:

        For example on the following binary tree `b`::

            |   ___3____      |
            |  /        \     |
            | 1         _7_   |
            |  \       /   \  |
            |   2     5     8 |
            |        / \      |
            |       4   6     |

        (only the nodes are shown), the depth-first post-order traversal
        algorithm explores `b` in the following order of nodes:
        `2,1,4,6,5,8,7,3`.

        For another example, consider the labelled tree::

            |     __1____ |
            |    /  /   / |
            |   2  6   8_ |
            |   |  |  / / |
            |   3_ 7 9 10 |
            |  / /        |
            | 4 5         |

        The algorithm explores this tree in the following order:
        `4,5,3,2,7,6,9,10,8,1`.

        TESTS::

            sage: b = BinaryTree([[None,[]],[[[],[]],[]]]).canonical_labelling()
            sage: ascii_art([b])
            [   ___3____      ]
            [  /        \     ]
            [ 1         _7_   ]
            [  \       /   \  ]
            [   2     5     8 ]
            [        / \      ]
            [       4   6     ]
            sage: [node.label() for node in b.post_order_traversal_iter()]
            [2, 1, 4, 6, 5, 8, 7, 3]

            sage: t = OrderedTree([[[[],[]]],[[]],[[],[]]]).canonical_labelling()
            sage: ascii_art([t])
            [     __1____   ]
            [    /  /   /   ]
            [   2  6   8_   ]
            [   |  |  / /   ]
            [   3_ 7 9 10   ]
            [  / /          ]
            [ 4 5           ]
            sage: [node.label() for node in t.post_order_traversal_iter()]
            [4, 5, 3, 2, 7, 6, 9, 10, 8, 1]

            sage: [node.label() for node in BinaryTree().canonical_labelling().\
            ....:     post_order_traversal_iter()]
            []
            sage: [node.label() for node in OrderedTree([]).\
            ....:     canonical_labelling().post_order_traversal_iter()]
            [1]

        The following test checks that things do not go wrong if some among
        the descendants of the tree are equal or even identical::

            sage: u = BinaryTree(None)
            sage: v = BinaryTree([u, u])
            sage: w = BinaryTree([v, v])
            sage: t = BinaryTree([w, w])
            sage: t.node_number()
            7
            sage: l = [1 for i in t.post_order_traversal_iter()]
            sage: len(l)
            7
        """
        if self.is_empty():
            return
        yield from itertools.chain(*[c.post_order_traversal_iter()
                                     for c in self])
        yield self

    def post_order_traversal(self, action=None):
        r"""
        Run the depth-first post-order traversal algorithm (recursive
        implementation) and subject every node encountered
        to some procedure ``action``. The algorithm is::

            explore each subtree (by the algorithm) from the
                leftmost one to the rightmost one;
            then manipulate the root with function `action` (in the
                case of a binary tree, only if the root is not a leaf).

        INPUT:

        - ``action`` -- (optional) a function which takes a node as
          input, and does something during the exploration

        OUTPUT:

        ``None``. (This is *not* an iterator.)

        .. SEEALSO::

            - :meth:`~sage.combinat.abstract_tree.AbstractTree.post_order_traversal_iter()`
            - :meth:`~sage.combinat.abstract_tree.AbstractTree.iterative_post_order_traversal()`

        TESTS::

            sage: l = []
            sage: b = BinaryTree([[None,[]],[[[],[]],[]]]).canonical_labelling()
            sage: b
            3[1[., 2[., .]], 7[5[4[., .], 6[., .]], 8[., .]]]
            sage: b.post_order_traversal(lambda node: l.append(node.label()))
            sage: l
            [2, 1, 4, 6, 5, 8, 7, 3]

            sage: t = OrderedTree([[[[],[]]],[[]],[[],[]]]).\
            ....:     canonical_labelling(); t
            1[2[3[4[], 5[]]], 6[7[]], 8[9[], 10[]]]
            sage: l = []
            sage: t.post_order_traversal(lambda node: l.append(node.label()))
            sage: l
            [4, 5, 3, 2, 7, 6, 9, 10, 8, 1]

            sage: l = []
            sage: BinaryTree().canonical_labelling().\
            ....:     post_order_traversal(lambda node: l.append(node.label()))
            sage: l
            []
            sage: OrderedTree([]).canonical_labelling().\
            ....:     post_order_traversal(lambda node: l.append(node.label()))
            sage: l
            [1]

        The following test checks that things do not go wrong if some among
        the descendants of the tree are equal or even identical::

            sage: u = BinaryTree(None)
            sage: v = BinaryTree([u, u])
            sage: w = BinaryTree([v, v])
            sage: t = BinaryTree([w, w])
            sage: t.node_number()
            7
            sage: l = []
            sage: t.post_order_traversal(lambda node: l.append(1))
            sage: len(l)
            7
        """
        if action is None:
            def action(x):
                return None
        for node in self.post_order_traversal_iter():
            action(node)

    def iterative_post_order_traversal(self, action=None):
        r"""
        Run the depth-first post-order traversal algorithm (iterative
        implementation) and subject every node encountered
        to some procedure ``action``. The algorithm is::

            explore each subtree (by the algorithm) from the
                leftmost one to the rightmost one;
            then manipulate the root with function `action` (in the
                case of a binary tree, only if the root is not a leaf).

        INPUT:

        - ``action`` -- (optional) a function which takes a node as
          input, and does something during the exploration

        OUTPUT:

        ``None``. (This is *not* an iterator.)

        .. SEEALSO::

            - :meth:`~sage.combinat.abstract_tree.AbstractTree.post_order_traversal_iter()`

        TESTS::

            sage: l = []
            sage: b = BinaryTree([[None,[]],[[[],[]],[]]]).canonical_labelling()
            sage: b
            3[1[., 2[., .]], 7[5[4[., .], 6[., .]], 8[., .]]]
            sage: b.iterative_post_order_traversal(lambda node: l.append(node.label()))
            sage: l
            [2, 1, 4, 6, 5, 8, 7, 3]

            sage: t = OrderedTree([[[[],[]]],[[]],[[],[]]]).canonical_labelling()
            sage: t
            1[2[3[4[], 5[]]], 6[7[]], 8[9[], 10[]]]
            sage: l = []
            sage: t.iterative_post_order_traversal(lambda node: l.append(node.label()))
            sage: l
            [4, 5, 3, 2, 7, 6, 9, 10, 8, 1]

            sage: l = []
            sage: BinaryTree().canonical_labelling().\
            ....:     iterative_post_order_traversal(
            ....:         lambda node: l.append(node.label()))
            sage: l
            []
            sage: OrderedTree([]).canonical_labelling().\
            ....:     iterative_post_order_traversal(
            ....:         lambda node: l.append(node.label()))
            sage: l
            [1]

        The following test checks that things do not go wrong if some among
        the descendants of the tree are equal or even identical::

            sage: u = BinaryTree(None)
            sage: v = BinaryTree([u, u])
            sage: w = BinaryTree([v, v])
            sage: t = BinaryTree([w, w])
            sage: t.node_number()
            7
            sage: l = []
            sage: t.iterative_post_order_traversal(lambda node: l.append(1))
            sage: len(l)
            7
        """
        if self.is_empty():
            return
        if action is None:
            def action(x):
                return None
        stack = [self]
        while stack:
            node = stack[-1]
            if node is not None:
                # A "None" on the stack means that the node right before
                # it on the stack has already been "exploded" into
                # subtrees, and should not be exploded again, but instead
                # should be manipulated and removed from the stack.
                stack.append(None)
                for i in range(len(node)):
                    subtree = node[-i - 1]
                    if not subtree.is_empty():
                        stack.append(subtree)
            else:
                stack.pop()
                node = stack.pop()
                action(node)

    def breadth_first_order_traversal(self, action=None):
        r"""
        Run the breadth-first post-order traversal algorithm
        and subject every node encountered to some procedure
        ``action``. The algorithm is::

            queue <- [ root ];
            while the queue is not empty:
                node <- pop( queue );
                manipulate the node;
                prepend to the queue the list of all subtrees of
                    the node (from the rightmost to the leftmost).

        INPUT:

        - ``action`` -- (optional) a function which takes a node as
          input, and does something during the exploration

        OUTPUT:

        ``None``. (This is *not* an iterator.)

        EXAMPLES:

        For example, on the following binary tree `b`::

            |   ___3____      |
            |  /        \     |
            | 1         _7_   |
            |  \       /   \  |
            |   2     5     8 |
            |        / \      |
            |       4   6     |

        the breadth-first order traversal algorithm explores `b` in the
        following order of nodes: `3,1,7,2,5,8,4,6`.

        TESTS::

            sage: b = BinaryTree([[None,[]],[[[],[]],[]]]).canonical_labelling()
            sage: l = []
            sage: b.breadth_first_order_traversal(lambda node: l.append(node.label()))
            sage: l
            [3, 1, 7, 2, 5, 8, 4, 6]

            sage: t = OrderedTree([[[[],[]]],[[]],[[],[]]]).canonical_labelling()
            sage: t
            1[2[3[4[], 5[]]], 6[7[]], 8[9[], 10[]]]
            sage: l = []
            sage: t.breadth_first_order_traversal(lambda node: l.append(node.label()))
            sage: l
            [1, 2, 6, 8, 3, 7, 9, 10, 4, 5]

            sage: l = []
            sage: BinaryTree().canonical_labelling().\
            ....:     breadth_first_order_traversal(
            ....:         lambda node: l.append(node.label()))
            sage: l
            []
            sage: OrderedTree([]).canonical_labelling().\
            ....:     breadth_first_order_traversal(
            ....:         lambda node: l.append(node.label()))
            sage: l
            [1]
        """
        if self.is_empty():
            return
        if action is None:
            def action(x):
                return None
        queue = []
        queue.append(self)
        while queue:
            node = queue.pop()
            action(node)
            for subtree in node:
                if not subtree.is_empty():
                    queue.insert(0, subtree)

    def paths_at_depth(self, depth, path=[]):
        r"""
        Return a generator for all paths at a fixed depth.

        This iterates over all paths for nodes that are at the given depth.

        Here the root is considered to have depth 0.

        INPUT:

        - ``depth`` -- integer
        - ``path`` -- (optional) list; given path used in the recursion

        .. WARNING::

            The ``path`` option should not be used directly.

        .. SEEALSO::

            :meth:`paths`, :meth:`paths_to_the_right`, :meth:`node_number_at_depth`

        EXAMPLES::

            sage: T = OrderedTree([[[], [[], [[]]]], [], [[[],[]]], [], []])
            sage: ascii_art(T)
                 ______o_______
                /    /   /  / /
              _o__  o   o  o o
             /   /      |
            o   o_      o_
               / /     / /
              o o     o o
                |
                o
            sage: list(T.paths_at_depth(0))
            [()]
            sage: list(T.paths_at_depth(2))
            [(0, 0), (0, 1), (2, 0)]
            sage: list(T.paths_at_depth(4))
            [(0, 1, 1, 0)]
            sage: list(T.paths_at_depth(5))
            []

            sage: T2 = OrderedTree([])
            sage: list(T2.paths_at_depth(0))
            [()]
        """
        if not depth:
            yield tuple(path)
        else:
            for i in range(len(self)):
                yield from self[i].paths_at_depth(depth - 1, path + [i])

    def node_number_at_depth(self, depth):
        r"""
        Return the number of nodes at a given depth.

        This counts all nodes that are at the given depth.

        Here the root is considered to have depth 0.

        INPUT:

        - ``depth`` -- integer

        .. SEEALSO::

            :meth:`node_number`, :meth:`node_number_to_the_right`, :meth:`paths_at_depth`

        EXAMPLES::

            sage: T = OrderedTree([[[], [[]]], [[], [[[]]]], []])
            sage: ascii_art(T)
                ___o____
               /    /  /
              o_   o_ o
             / /  / /
            o o  o o
              |    |
              o    o
                   |
                   o
            sage: [T.node_number_at_depth(i) for i in range(6)]
            [1, 3, 4, 2, 1, 0]

        TESTS:

        Check that the empty tree has no nodes (:issue:`29134`)::

            sage: T = BinaryTree()
            sage: T
            .
            sage: T.is_empty()
            True
            sage: [T.node_number_at_depth(i) for i in range(3)]
            [0, 0, 0]
        """
        if self.is_empty():
            return Integer(0)
        if depth == 0:
            return Integer(1)
        return sum(son.node_number_at_depth(depth - 1) for son in self)

    def paths_to_the_right(self, path):
        r"""
        Return a generator of paths for all nodes at the same
        depth and to the right of the node identified by ``path``.

        This iterates over the paths for nodes that are at the same
        depth as the given one, and strictly to its right.

        INPUT:

        - ``path`` -- any path in the tree

        .. SEEALSO::

            :meth:`paths`, :meth:`paths_at_depth`, :meth:`node_number_to_the_right`

        EXAMPLES::

            sage: T = OrderedTree([[[], [[]]], [[], [[[]]]], []])
            sage: ascii_art(T)
                ___o____
               /    /  /
              o_   o_ o
             / /  / /
            o o  o o
              |    |
              o    o
                   |
                   o
            sage: g = T.paths_to_the_right(())
            sage: list(g)
            []

            sage: g = T.paths_to_the_right((0,))
            sage: list(g)
            [(1,), (2,)]

            sage: g = T.paths_to_the_right((0,1))
            sage: list(g)
            [(1, 0), (1, 1)]

            sage: g = T.paths_to_the_right((0,1,0))
            sage: list(g)
            [(1, 1, 0)]

            sage: g = T.paths_to_the_right((1,2))
            sage: list(g)
            []
        """
        depth = len(path)
        if (not depth) or path[0] >= len(self):
            return
        for i in range(path[0] + 1, len(self)):
            for p in self[i].paths_at_depth(depth - 1, path=[i]):
                yield p
        for p in self[path[0]].paths_to_the_right(path[1:]):
            yield tuple([path[0]] + list(p))

    def node_number_to_the_right(self, path):
        r"""
        Return the number of nodes at the same depth and to the right of
        the node identified by ``path``.

        This counts the nodes that are at the same depth as the given
        one, and strictly to its right.

        .. SEEALSO::

            :meth:`node_number`, :meth:`node_number_at_depth`, :meth:`paths_to_the_right`

        EXAMPLES::

            sage: T = OrderedTree([[[], [[]]], [[], [[[]]]], []])
            sage: ascii_art(T)
                ___o____
               /    /  /
              o_   o_ o
             / /  / /
            o o  o o
              |    |
              o    o
                   |
                   o
            sage: T.node_number_to_the_right(())
            0
            sage: T.node_number_to_the_right((0,))
            2
            sage: T.node_number_to_the_right((0,1))
            2
            sage: T.node_number_to_the_right((0,1,0))
            1

            sage: T = OrderedTree([])
            sage: T.node_number_to_the_right(())
            0
        """
        depth = len(path)
        if depth == 0:
            return Integer(0)
        result = sum(son.node_number_at_depth(depth - 1)
                     for son in self[path[0] + 1:])
        if path[0] < len(self) and path[0] >= 0:
            result += self[path[0]].node_number_to_the_right(path[1:])
        return result

    def subtrees(self):
        """
        Return a generator for all nonempty subtrees of ``self``.

        The number of nonempty subtrees of a tree is its number of
        nodes. (The word "nonempty" makes a difference only in the
        case of binary trees. For ordered trees, for example, all
        trees are nonempty.)

        EXAMPLES::

            sage: list(OrderedTree([]).subtrees())
            [[]]
            sage: list(OrderedTree([[],[[]]]).subtrees())
            [[[], [[]]], [], [[]], []]

            sage: list(OrderedTree([[],[[]]]).canonical_labelling().subtrees())
            [1[2[], 3[4[]]], 2[], 3[4[]], 4[]]

            sage: list(BinaryTree([[],[[],[]]]).subtrees())
            [[[., .], [[., .], [., .]]], [., .], [[., .], [., .]], [., .], [., .]]

            sage: v = BinaryTree([[],[]])
            sage: list(v.canonical_labelling().subtrees())
            [2[1[., .], 3[., .]], 1[., .], 3[., .]]

        TESTS::

            sage: t = OrderedTree([[], [[], [[], []], [[], []]], [[], []]])
            sage: t.node_number() == len(list(t.subtrees()))
            True
            sage: list(BinaryTree().subtrees())
            []
            sage: bt = BinaryTree([[],[[],[]]])
            sage: bt.node_number() == len(list(bt.subtrees()))
            True
        """
        return self.pre_order_traversal_iter()

    def paths(self):
        """
        Return a generator for all paths to nodes of ``self``.

        OUTPUT:

        This method returns a list of sequences of integers. Each of these
        sequences represents a path from the root node to some node. For
        instance, `(1, 3, 2, 5, 0, 3)` represents the node obtained by
        choosing the 1st child of the root node (in the ordering returned
        by ``iter``), then the 3rd child of its child, then the 2nd child
        of the latter, etc. (where the labelling of the children is
        zero-based).

        The root element is represented by the empty tuple ``()``.

        .. SEEALSO::

            :meth:`paths_at_depth`, :meth:`paths_to_the_right`

        EXAMPLES::

            sage: list(OrderedTree([]).paths())
            [()]
            sage: list(OrderedTree([[],[[]]]).paths())
            [(), (0,), (1,), (1, 0)]

            sage: list(BinaryTree([[],[[],[]]]).paths())
            [(), (0,), (1,), (1, 0), (1, 1)]

        TESTS::

            sage: t = OrderedTree([[], [[], [[], []], [[], []]], [[], []]])
            sage: t.node_number() == len(list(t.paths()))
            True
            sage: list(BinaryTree().paths())
            []
            sage: bt = BinaryTree([[],[[],[]]])
            sage: bt.node_number() == len(list(bt.paths()))
            True
        """
        if not self.is_empty():
            yield ()
            for i, t in enumerate(self):
                for p in t.paths():
                    yield (i,) + p

    def node_number(self):
        """
        Return the number of nodes of ``self``.

        .. SEEALSO::

            :meth:`node_number_at_depth`, :meth:`node_number_to_the_right`

        EXAMPLES::

            sage: OrderedTree().node_number()
            1
            sage: OrderedTree([]).node_number()
            1
            sage: OrderedTree([[],[]]).node_number()
            3
            sage: OrderedTree([[],[[]]]).node_number()
            4
            sage: OrderedTree([[], [[], [[], []], [[], []]], [[], []]]).node_number()
            13

        EXAMPLES::

            sage: BinaryTree(None).node_number()
            0
            sage: BinaryTree([]).node_number()
            1
            sage: BinaryTree([[], None]).node_number()
            2
            sage: BinaryTree([[None, [[], []]], None]).node_number()
            5
        """
        if self.is_empty():
            return Integer(0)
        else:
            return sum((i.node_number() for i in self), Integer(1))

    def depth(self):
        """
        Return the depth of ``self``.

        EXAMPLES::

            sage: OrderedTree().depth()
            1
            sage: OrderedTree([]).depth()
            1
            sage: OrderedTree([[],[]]).depth()
            2
            sage: OrderedTree([[],[[]]]).depth()
            3
            sage: OrderedTree([[], [[], [[], []], [[], []]], [[], []]]).depth()
            4

            sage: BinaryTree().depth()
            0
            sage: BinaryTree([[],[[],[]]]).depth()
            3
        """
        if self:
            return Integer(1 + max(i.depth() for i in self))
        else:
            return Integer(0 if self.is_empty() else 1)

    def _ascii_art_(self):
        r"""
        TESTS::

            sage: t = OrderedTree([])
            sage: ascii_art(t)
            o
            sage: t = OrderedTree([[]])
            sage: aa = ascii_art(t);aa
            o
            |
            o
            sage: aa.get_baseline()
            2
            sage: tt1 = OrderedTree([[],[[],[],[[[[]]]]],[[[],[],[],[]]]])
            sage: ascii_art(tt1)
              _____o_______
             /    /       /
            o   _o__     o
               / / /     |
              o o o    __o___
                  |   / / / /
                  o  o o o o
                  |
                  o
                  |
                  o
            sage: ascii_art(tt1.canonical_labelling())
              ______1_______
             /    /        /
            2   _3__      10
               / / /      |
              4 5 6    ___11____
                  |   /  /  /  /
                  7  12 13 14 15
                  |
                  8
                  |
                  9
            sage: ascii_art(OrderedTree([[],[[]]]))
                  o_
                 / /
                o o
                  |
                  o
            sage: t = OrderedTree([[[],[[[],[]]],[[]]],[[[[[],[]]]]],[[],[]]])
            sage: ascii_art(t)
                  _____o_______
                 /       /    /
              __o____   o    o_
             /   /  /   |   / /
            o   o  o    o  o o
                |  |    |
                o_ o    o
               / /      |
              o o       o_
                       / /
                      o o
            sage: ascii_art(t.canonical_labelling())
                  ______1________
                 /       /      /
              __2____   10     16_
             /   /  /   |     /  /
            3   4  8    11   17 18
                |  |    |
                5_ 9    12
               / /      |
              6 7       13_
                       /  /
                      14 15
        """
        def node_to_str(t):
            return str(t.label()) if hasattr(t, "label") else "o"

        if self.is_empty():
            from sage.typeset.ascii_art import empty_ascii_art
            return empty_ascii_art

        from sage.typeset.ascii_art import AsciiArt
        if len(self) == 0:
            t_repr = AsciiArt([node_to_str(self)])
            t_repr._root = 1
            return t_repr
        if len(self) == 1:
            repr_child = self[0]._ascii_art_()
            sep = AsciiArt([" " * (repr_child._root - 1)])
            t_repr = AsciiArt([node_to_str(self)])
            t_repr._root = 1
            repr_root = (sep + t_repr) * (sep + AsciiArt(["|"]))
            t_repr = repr_root * repr_child
            t_repr._root = repr_child._root
            t_repr._baseline = t_repr._h - 1
            return t_repr
        # General case
        l_repr = [subtree._ascii_art_() for subtree in self]
        acc = l_repr.pop(0)
        whitesep = acc._root + 1
        lf_sep = " " * (acc._root + 1) + "_" * (acc._l - acc._root)
        ls_sep = " " * (acc._root) + "/" + " " * (acc._l - acc._root)
        while l_repr:
            t_repr = l_repr.pop(0)
            acc += AsciiArt([" "]) + t_repr
            if len(l_repr) == 0:
                lf_sep += "_" * (t_repr._root + 1)
            else:
                lf_sep += "_" * (t_repr._l + 1)
            ls_sep += " " * (t_repr._root) + "/" + " " * (t_repr._l - t_repr._root)
        mid = whitesep + (len(lf_sep) - whitesep) // 2
        node = node_to_str(self)
        t_repr = AsciiArt([lf_sep[:mid - 1] + node + lf_sep[mid + len(node) - 1:], ls_sep]) * acc
        t_repr._root = mid
        t_repr._baseline = t_repr._h - 1
        return t_repr

    def _unicode_art_(self):
        r"""
        TESTS::

            sage: t = OrderedTree([])
            sage: unicode_art(t)
            o
            sage: t = OrderedTree([[]])
            sage: aa = unicode_art(t);aa
            o
            │
            o
            sage: aa.get_baseline()
            2
            sage: tt1 = OrderedTree([[],[[],[],[[[[]]]]],[[[],[],[],[]]]])
            sage: unicode_art(tt1)
            ╭───┬─o────╮
            │   │      │
            o ╭─o─╮    o
              │ │ │    │
              o o o ╭─┬o┬─╮
                  │ │ │ │ │
                  o o o o o
                  │
                  o
                  │
                  o
            sage: unicode_art(tt1.canonical_labelling())
            ╭───┬──1─────╮
            │   │        │
            2 ╭─3─╮      10
              │ │ │      │
              4 5 6 ╭──┬11┬──╮
                  │ │  │  │  │
                  7 12 13 14 15
                  │
                  8
                  │
                  9
            sage: unicode_art(OrderedTree([[],[[]]]))
            ╭o╮
            │ │
            o o
              │
              o
            sage: t = OrderedTree([[[],[[[],[]]],[[]]],[[[[[],[]]]]],[[],[]]])
            sage: unicode_art(t)
               ╭────o┬───╮
               │     │   │
            ╭──o──╮  o  ╭o╮
            │  │  │  │  │ │
            o  o  o  o  o o
               │  │  │
              ╭o╮ o  o
              │ │    │
              o o   ╭o╮
                    │ │
                    o o
            sage: unicode_art(t.canonical_labelling())
               ╭──────1─────╮
               │      │     │
            ╭──2──╮   10  ╭16╮
            │  │  │   │   │  │
            3  4  8   11  17 18
               │  │   │
              ╭5╮ 9   12
              │ │     │
              6 7   ╭13╮
                    │  │
                    14 15
        """

        def node_to_str(t):
            if hasattr(t, "label"):
                return str(t.label())
            else:
                return "o"
        # other possible choices for nodes would be u"█ ▓ ░ ╋ ╬"

        if self.is_empty():
            from sage.typeset.unicode_art import empty_unicode_art
            return empty_unicode_art

        from sage.typeset.unicode_art import UnicodeArt
        if not len(self):
            t_repr = UnicodeArt([node_to_str(self)])
            t_repr._root = 0
            return t_repr

        if len(self) == 1:
            repr_child = self[0]._unicode_art_()
            sep = UnicodeArt([" " * repr_child._root])
            t_repr = UnicodeArt([node_to_str(self)])
            repr_root = (sep + t_repr) * (sep + UnicodeArt(["│"]))
            t_repr = repr_root * repr_child
            t_repr._root = repr_child._root
            t_repr._baseline = t_repr._h - 1
            return t_repr

        # General case
        l_repr = [subtree._unicode_art_() for subtree in self]
        acc = l_repr.pop(0)
        whitesep = acc._root
        lf_sep = " " * whitesep + "╭" + "─" * (acc._l - acc._root)
        ls_sep = " " * whitesep + "│" + " " * (acc._l - acc._root)
        while l_repr:
            tr = l_repr.pop(0)
            acc += UnicodeArt([" "]) + tr
            if not len(l_repr):
                lf_sep += "─" * (tr._root) + "╮"
                ls_sep += " " * (tr._root) + "│"
            else:
                lf_sep += "─" * (tr._root) + "┬" + "─" * (tr._l - tr._root)
                ls_sep += " " * (tr._root) + "│" + " " * (tr._l - tr._root)
        mid = whitesep + (len(lf_sep) - whitesep) // 2
        node = node_to_str(self)
        lf_sep = (lf_sep[:mid - len(node) // 2] + node +
                  lf_sep[mid + len(node) - len(node) // 2:])
        t_repr = UnicodeArt([lf_sep, ls_sep]) * acc
        t_repr._root = mid
        t_repr._baseline = t_repr._h - 1
        return t_repr

    def canonical_labelling(self, shift=1):
        """
        Return a labelled version of ``self``.

        The actual canonical labelling is currently unspecified. However, it
        is guaranteed to have labels in `1...n` where `n` is the number of
        nodes of the tree. Moreover, two (unlabelled) trees compare as equal if
        and only if their canonical labelled trees compare as equal.

        EXAMPLES::

            sage: t = OrderedTree([[], [[], [[], []], [[], []]], [[], []]])
            sage: t.canonical_labelling()
            1[2[], 3[4[], 5[6[], 7[]], 8[9[], 10[]]], 11[12[], 13[]]]

            sage: BinaryTree([]).canonical_labelling()
            1[., .]
            sage: BinaryTree([[],[[],[]]]).canonical_labelling()
            2[1[., .], 4[3[., .], 5[., .]]]

        TESTS::

            sage: BinaryTree().canonical_labelling()
            .
        """
        LTR = self.parent().labelled_trees()
        liste = []
        deca = 1
        for subtree in self:
            liste += [subtree.canonical_labelling(shift + deca)]
            deca += subtree.node_number()
        return LTR._element_constructor_(liste, label=shift)

    def to_hexacode(self):
        r"""
        Transform a tree into a hexadecimal string.

        The definition of the hexacode is recursive. The first letter is
        the valence of the root as a hexadecimal (up to 15), followed by
        the concatenation of the hexacodes of the subtrees.

        This method only works for trees where every vertex has
        valency at most 15.

        See :func:`from_hexacode` for the reverse transformation.

        EXAMPLES::

            sage: from sage.combinat.abstract_tree import from_hexacode
            sage: LT = LabelledOrderedTrees()
            sage: from_hexacode('2010', LT).to_hexacode()
            '2010'
            sage: LT.an_element().to_hexacode()
            '3020010'
            sage: t = from_hexacode('a0000000000000000', LT)
            sage: t.to_hexacode()
            'a0000000000'

            sage: OrderedTrees(6).an_element().to_hexacode()
            '500000'

        TESTS::

            sage: one = LT([], label='@')
            sage: LT([one for _ in range(15)], label='@').to_hexacode()
            'f000000000000000'
            sage: LT([one for _ in range(16)], label='@').to_hexacode()
            Traceback (most recent call last):
            ...
            ValueError: the width of the tree is too large
        """
        if len(self) > 15:
            raise ValueError("the width of the tree is too large")
        if self.node_number() == 1:
            return "0"
        return ("%x" % len(self)) + "".join(u.to_hexacode() for u in self)

    def tree_factorial(self):
        r"""
        Return the tree-factorial of ``self``.

        Definition:

        The tree-factorial `T!` of a tree `T` is the product `\prod_{v\in
        T}\#\mbox{children}(v)`.

        EXAMPLES::

            sage: LT = LabelledOrderedTrees()
            sage: t = LT([LT([],label=6),LT([],label=1)],label=9)
            sage: t.tree_factorial()
            3

            sage: BinaryTree([[],[[],[]]]).tree_factorial()
            15

        TESTS::

            sage: BinaryTree().tree_factorial()
            1
        """
        nb = self.node_number()
        if nb <= 1:
            return Integer(1)
        return nb * prod(s.tree_factorial() for s in self)

    def _latex_(self):
        r"""
        Generate `\LaTeX` output which can be easily modified.

        TESTS::

            sage: latex(BinaryTree([[[],[]],[[],None]]))
            { \newcommand{\nodea}{\node[draw,circle] (a) {$$}
            ;}\newcommand{\nodeb}{\node[draw,circle] (b) {$$}
            ;}\newcommand{\nodec}{\node[draw,circle] (c) {$$}
            ;}\newcommand{\noded}{\node[draw,circle] (d) {$$}
            ;}\newcommand{\nodee}{\node[draw,circle] (e) {$$}
            ;}\newcommand{\nodef}{\node[draw,circle] (f) {$$}
            ;}\begin{tikzpicture}[auto]
            \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                     \&         \&         \& \nodea  \&         \&         \&         \\
                     \& \nodeb  \&         \&         \&         \& \nodee  \&         \\
             \nodec  \&         \& \noded  \&         \& \nodef  \&         \&         \\
            };
            <BLANKLINE>
            \path[ultra thick, red] (b) edge (c) edge (d)
                (e) edge (f)
                (a) edge (b) edge (e);
            \end{tikzpicture}}
        """
        #######################################################################
        # load tikz in the preamble for *view*
        from sage.misc.latex import latex
        latex.add_package_to_preamble_if_available("tikz")
        #######################################################################
        # latex environment : TikZ
        begin_env = "\\begin{tikzpicture}[auto]\n"
        end_env = "\\end{tikzpicture}"
        # it uses matrix trick to place each node
        matrix_begin = "\\matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\\&]{\n"
        matrix_end = "\\\\\n};\n"
        # a basic path to each edges
        path_begin = "\\path[ultra thick, red] "
        path_end = ";\n"
        # to make a pretty output, it creates one LaTeX command for
        # each node
        cmd = "\\node"
        new_cmd1 = "\\newcommand{" + cmd
        new_cmd2 = "}{\\node[draw,circle] ("
        new_cmd3 = ") {$"
        new_cmd4 = "$}\n;}"
        # some variables to simplify code
        sep = "\\&"
        space = " " * 9
        sepspace = sep + space
        spacesep = space + sep

        def node_to_str(node):
            return " " + node + " " * (len(space) - 1 - len(node))
        # # TODO:: modify how to create nodes --> new_cmd : \\node[...] in create_node
        num = [0]

        def resolve(self):
            nodes = []
            matrix = []
            edges = []

            def create_node(self):
                r"""
                create a name (infixe reading)
                 -> ex: b
                create a new command:
                 -> ex: \newcommand{\nodeb}{\node[draw,circle] (b) {$$};
                return the name and the command to build:
                  . the matrix
                  . and the edges
                """
                name = "".join(chr(ord(x) + 49) for x in str(num[0]))
                node = cmd + name
                nodes.append((name,
                    (str(self.label()) if hasattr(self, "label") else ""))
                )
                num[0] += 1
                return node, name

            def empty_tree():
                r"""
                TESTS::

                    sage: t = BinaryTree()
                    sage: print(latex(t))
                    { \begin{tikzpicture}[auto]
                    \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                             \\
                    };
                    \end{tikzpicture}}
                """
                matrix.append(space)

            def one_node_tree(self):
                r"""
                TESTS::

                    sage: t = BinaryTree([]); print(latex(t))
                    { \newcommand{\nodea}{\node[draw,circle] (a) {$$}
                    ;}\begin{tikzpicture}[auto]
                    \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                     \nodea  \\
                    };
                    \end{tikzpicture}}
                    sage: t = OrderedTree([]); print(latex(t))
                    { \newcommand{\nodea}{\node[draw,circle] (a) {$$}
                    ;}\begin{tikzpicture}[auto]
                    \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                     \nodea  \\
                    };
                    \end{tikzpicture}}
                """
                node, _ = create_node(self)
                matrix.append(node_to_str(node))

            def concat_matrix(mat, mat2):
                lmat = len(mat)
                lmat2 = len(mat2)
                for i in range(max(lmat, lmat2)):
                    # mat[i] --> n & n & ...
                    # mat2[i] -> n' & n' & ...
                    # ==> n & n & ... & n' & n' & ...
                    try:
                        mat[i] += sep + mat2[i]
                    except IndexError:
                        if i >= lmat:
                            if i != 0:
                                # mat[i] does not exist but
                                # mat[0] has k "&"
                                # mat2[i] -> n' & n' & ...
                                # ==> (_ &)*k+1 n' & n' & ...
                                nb_of_and = mat[0].count(sep) - mat2[0].count(sep)
                                mat.append(spacesep * (nb_of_and) + mat2[i])
                            else:
                                # mat is empty
                                # mat2[i] -> n' & n' & ...
                                # ==> mat2
                                mat.extend(mat2)
                                return
                        else:
                            # mat[i] -> n & n & ...
                            # mat2[i] does not exist but mat2[0] exists
                            # # and has k "&"
                            # NOTE:: i != 0 because that is a no-empty subtree.
                            # ==> n & n & ... (& _)*k+1
                            nb_of_and = mat2[0].count(sep)
                            mat[i] += sepspace * (nb_of_and + 1)

            def tmp(subtree, edge, nodes, edges, matrix):
                if not subtree.is_empty():
                    # # create representation of the subtree
                    nodes_st, matrix_st, edges_st = resolve(subtree)
                    # # add its nodes to the "global" nodes set
                    nodes.extend(nodes_st)
                    # # create a new edge between the root and the subtree
                    edge.append(nodes_st[0][0])
                    # # add the subtree edges to the "global" edges set
                    edges.extend(edges_st)
                    # # build a new matrix by concatenation
                    concat_matrix(matrix, matrix_st)
                else:
                    concat_matrix(matrix, [space])

            def pair_nodes_tree(self, nodes, edges, matrix):
                r"""
                TESTS::

                    sage: t = OrderedTree([[[],[]],[[],[]]]).\
                    ....:     canonical_labelling(); print(latex(t))
                    { \newcommand{\nodea}{\node[draw,circle] (a) {$1$}
                    ;}\newcommand{\nodeb}{\node[draw,circle] (b) {$2$}
                    ;}\newcommand{\nodec}{\node[draw,circle] (c) {$3$}
                    ;}\newcommand{\noded}{\node[draw,circle] (d) {$4$}
                    ;}\newcommand{\nodee}{\node[draw,circle] (e) {$5$}
                    ;}\newcommand{\nodef}{\node[draw,circle] (f) {$6$}
                    ;}\newcommand{\nodeg}{\node[draw,circle] (g) {$7$}
                    ;}\begin{tikzpicture}[auto]
                    \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                             \&         \&         \& \nodea  \&         \&         \&         \\
                             \& \nodeb  \&         \&         \&         \& \nodee  \&         \\
                     \nodec  \&         \& \noded  \&         \& \nodef  \&         \& \nodeg  \\
                    };
                    <BLANKLINE>
                    \path[ultra thick, red] (b) edge (c) edge (d)
                        (e) edge (f) edge (g)
                        (a) edge (b) edge (e);
                    \end{tikzpicture}}
                    sage: t = BinaryTree([[],[[],[]]]); print(latex(t))
                    { \newcommand{\nodea}{\node[draw,circle] (a) {$$}
                    ;}\newcommand{\nodeb}{\node[draw,circle] (b) {$$}
                    ;}\newcommand{\nodec}{\node[draw,circle] (c) {$$}
                    ;}\newcommand{\noded}{\node[draw,circle] (d) {$$}
                    ;}\newcommand{\nodee}{\node[draw,circle] (e) {$$}
                    ;}\begin{tikzpicture}[auto]
                    \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                             \& \nodea  \&         \&         \&         \\
                     \nodeb  \&         \&         \& \nodec  \&         \\
                             \&         \& \noded  \&         \& \nodee  \\
                    };
                    <BLANKLINE>
                    \path[ultra thick, red] (c) edge (d) edge (e)
                        (a) edge (b) edge (c);
                    \end{tikzpicture}}
                """
                # build all subtree matrices.
                node, name = create_node(self)
                edge = [name]
                split = len(self) // 2
                # the left part
                for i in range(split):
                    tmp(self[i], edge, nodes, edges, matrix)
                # # prepare the root line
                nb_of_and = matrix[0].count(sep)
                # the middle
                for i in range(len(matrix)):
                    matrix[i] += sepspace
                # the right part
                for i in range(split, len(self)):
                    tmp(self[i], edge, nodes, edges, matrix)

                # # create the root line
                root_line = (spacesep * (nb_of_and + 1) + node_to_str(node) +
                    sepspace * (matrix[0].count(sep) - nb_of_and - 1))
                matrix.insert(0, root_line)
                # add edges from the root
                edges.append(edge)

            def odd_nodes_tree(self, nodes, edges, matrix):
                r"""
                TESTS::

                    sage: t = OrderedTree([[]]).canonical_labelling()
                    sage: print(latex(t))
                    { \newcommand{\nodea}{\node[draw,circle] (a) {$1$}
                    ;}\newcommand{\nodeb}{\node[draw,circle] (b) {$2$}
                    ;}\begin{tikzpicture}[auto]
                    \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                     \nodea  \\
                     \nodeb  \\
                    };
                    <BLANKLINE>
                    \path[ultra thick, red] (a) edge (b);
                    \end{tikzpicture}}
                    sage: t = OrderedTree([[[],[]]]).canonical_labelling(); print(latex(t))
                    { \newcommand{\nodea}{\node[draw,circle] (a) {$1$}
                    ;}\newcommand{\nodeb}{\node[draw,circle] (b) {$2$}
                    ;}\newcommand{\nodec}{\node[draw,circle] (c) {$3$}
                    ;}\newcommand{\noded}{\node[draw,circle] (d) {$4$}
                    ;}\begin{tikzpicture}[auto]
                    \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                             \& \nodea  \&         \\
                             \& \nodeb  \&         \\
                     \nodec  \&         \& \noded  \\
                    };
                    <BLANKLINE>
                    \path[ultra thick, red] (b) edge (c) edge (d)
                        (a) edge (b);
                    \end{tikzpicture}}
                    sage: t = OrderedTree([[[],[],[]]]).canonical_labelling(); print(latex(t))
                    { \newcommand{\nodea}{\node[draw,circle] (a) {$1$}
                    ;}\newcommand{\nodeb}{\node[draw,circle] (b) {$2$}
                    ;}\newcommand{\nodec}{\node[draw,circle] (c) {$3$}
                    ;}\newcommand{\noded}{\node[draw,circle] (d) {$4$}
                    ;}\newcommand{\nodee}{\node[draw,circle] (e) {$5$}
                    ;}\begin{tikzpicture}[auto]
                    \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                             \& \nodea  \&         \\
                             \& \nodeb  \&         \\
                     \nodec  \& \noded  \& \nodee  \\
                    };
                    <BLANKLINE>
                    \path[ultra thick, red] (b) edge (c) edge (d) edge (e)
                        (a) edge (b);
                    \end{tikzpicture}}
                    sage: t = OrderedTree([[[],[],[]],[],[]]).canonical_labelling(); print(latex(t))
                    { \newcommand{\nodea}{\node[draw,circle] (a) {$1$}
                    ;}\newcommand{\nodeb}{\node[draw,circle] (b) {$2$}
                    ;}\newcommand{\nodec}{\node[draw,circle] (c) {$3$}
                    ;}\newcommand{\noded}{\node[draw,circle] (d) {$4$}
                    ;}\newcommand{\nodee}{\node[draw,circle] (e) {$5$}
                    ;}\newcommand{\nodef}{\node[draw,circle] (f) {$6$}
                    ;}\newcommand{\nodeg}{\node[draw,circle] (g) {$7$}
                    ;}\begin{tikzpicture}[auto]
                    \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                             \&         \&         \& \nodea  \&         \\
                             \& \nodeb  \&         \& \nodef  \& \nodeg  \\
                     \nodec  \& \noded  \& \nodee  \&         \&         \\
                    };
                    <BLANKLINE>
                    \path[ultra thick, red] (b) edge (c) edge (d) edge (e)
                        (a) edge (b) edge (f) edge (g);
                    \end{tikzpicture}}
                """
                # build all subtree matrices.
                node, name = create_node(self)
                edge = [name]
                split = len(self) // 2
                # the left part
                for i in range(split):
                    tmp(self[i], edge, nodes, edges, matrix)
                # # prepare the root line
                if matrix:
                    nb_of_and = matrix[0].count(sep)
                    sizetmp = len(matrix[0])
                else:
                    nb_of_and = 0
                    sizetmp = 0
                # the middle
                tmp(self[split], edge, nodes, edges, matrix)
                nb_of_and += matrix[0][sizetmp:].split("node")[0].count(sep)

                # the right part
                for i in range(split + 1, len(self)):
                    tmp(self[i], edge, nodes, edges, matrix)

                # # create the root line
                root_line = (spacesep * (nb_of_and) + node_to_str(node) +
                    sepspace * (matrix[0].count(sep) - nb_of_and))
                matrix.insert(0, root_line)
                # add edges from the root
                edges.append(edge)
            if self.is_empty():
                empty_tree()
            elif len(self) == 0 or all(subtree.is_empty()
                                       for subtree in self):
                one_node_tree(self)
            elif not len(self) % 2:
                pair_nodes_tree(self, nodes, edges, matrix)
            else:
                odd_nodes_tree(self, nodes, edges, matrix)
            return nodes, matrix, edges

        nodes, matrix, edges = resolve(self)

        def make_cmd(nodes):
            cmds = []
            for name, label in nodes:
                cmds.append(new_cmd1 + name + new_cmd2 +
                    name + new_cmd3 +
                    label + new_cmd4)
            return cmds

        def make_edges(edges):
            all_paths = []
            for edge in edges:
                path = "(" + edge[0] + ")"
                for i in range(1, len(edge)):
                    path += " edge (%s)" % edge[i]
                all_paths.append(path)
            return all_paths
        return ("{ " +
            "".join(make_cmd(nodes)) +
            begin_env +
                (matrix_begin +
                    "\\\\ \n".join(matrix) +
                matrix_end +
                ("\n" +
                path_begin +
                    "\n\t".join(make_edges(edges)) +
                path_end if edges else "")
                if matrix else "") +
            end_env +
            "}")


class AbstractClonableTree(AbstractTree):
    """
    Abstract Clonable Tree.

    An abstract class for trees with clone protocol (see
    :mod:`~sage.structure.list_clone`). It is expected that classes extending
    this one may also inherit from classes like :class:`~sage.structure.list_clone.ClonableArray` or
    :class:`~sage.structure.list_clone.ClonableList` depending whether one
    wants to build trees where adding a child is allowed.

    .. NOTE:: Due to the limitation of Cython inheritance, one cannot inherit
       here from :class:`~sage.structure.list_clone.ClonableElement`, because
       it would prevent us from later inheriting from
       :class:`~sage.structure.list_clone.ClonableArray` or
       :class:`~sage.structure.list_clone.ClonableList`.

    .. rubric:: How should this class be extended ?

    A class extending :class:`AbstractClonableTree
    <sage.combinat.abstract_tree.AbstractClonableTree>` should satisfy the
    following assumptions:

    * An instantiable class extending :class:`AbstractClonableTree
      <sage.combinat.abstract_tree.AbstractClonableTree>` should also extend
      the :class:`ClonableElement <sage.structure.list_clone.ClonableElement>`
      class or one of its subclasses generally, at least :class:`ClonableArray
      <sage.structure.list_clone.ClonableArray>`.

    * To respect the Clone protocol, the :meth:`AbstractClonableTree.check`
      method should be overridden by the new class.

    See also the assumptions in :class:`AbstractTree`.
    """

    def check(self):
        """
        Check that ``self`` is a correct tree.

        This method does nothing. It is implemented here because many
        extensions of :class:`AbstractClonableTree
        <sage.combinat.abstract_tree.AbstractClonableTree>` also extend
        :class:`sage.structure.list_clone.ClonableElement`, which requires it.

        It should be overridden in subclasses in order to check that the
        characterizing property of the respective kind of tree holds (eg: two
        children for binary trees).

        EXAMPLES::

            sage: OrderedTree([[],[[]]]).check()
            sage: BinaryTree([[],[[],[]]]).check()
        """
        pass

    def __setitem__(self, idx, value):
        """
        Substitute a subtree.

        .. NOTE::

            The tree ``self`` must be in a mutable state. See
            :mod:`sage.structure.list_clone` for more details about
            mutability.  The default implementation here assume that the
            container of the node implement a method ``_setitem`` with signature
            `self._setitem(idx, value)`. It is usually provided by inheriting
            from :class:`~sage.structure.list_clone.ClonableArray`.

        INPUT:

        - ``idx`` -- a valid path in ``self`` identifying a node

        - ``value`` -- the tree to be substituted

        EXAMPLES:

        Trying to modify a non mutable tree raises an error::

            sage: x = OrderedTree([])
            sage: x[0] =  OrderedTree([[]])
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.

        Here is the correct way to do it::

            sage: x = OrderedTree([[],[[]]])
            sage: with x.clone() as x:
            ....:     x[0] = OrderedTree([[]])
            sage: x
            [[[]], [[]]]

        One can also substitute at any depth::

            sage: y = OrderedTree(x)
            sage: with x.clone() as x:
            ....:     x[0,0] = OrderedTree([[]])
            sage: x
            [[[[]]], [[]]]
            sage: y
            [[[]], [[]]]
            sage: with y.clone() as y:
            ....:     y[(0,)] = OrderedTree([])
            sage: y
            [[], [[]]]

        This works for binary trees as well::

            sage: bt = BinaryTree([[],[[],[]]]); bt
            [[., .], [[., .], [., .]]]
            sage: with bt.clone() as bt1:
            ....:     bt1[0,0] = BinaryTree([[[], []], None])
            sage: bt1
            [[[[[., .], [., .]], .], .], [[., .], [., .]]]

        TESTS::

            sage: x = OrderedTree([])
            sage: with x.clone() as x:
            ....:     x[0] = OrderedTree([[]])
            Traceback (most recent call last):
            ...
            IndexError: list assignment index out of range

            sage: x = OrderedTree([]); x = OrderedTree([x,x]); x = OrderedTree([x,x]); x = OrderedTree([x,x])
            sage: with x.clone() as x:
            ....:     x[0,0] = OrderedTree()
            sage: x
            [[[], [[], []]], [[[], []], [[], []]]]
        """
        if not isinstance(value, self.__class__):
            raise TypeError('the given value is not a tree')
        if isinstance(idx, tuple):
            self.__setitem_rec__(idx, 0, value)
        else:
            self._setitem(idx, value)

    def __setitem_rec__(self, idx, i, value):
        """
        TESTS::

            sage: x = OrderedTree([[[], []],[[]]])
            sage: with x.clone() as x:              # indirect doctest
            ....:     x[0,1] = OrderedTree([[[]]])
            sage: x
            [[[], [[[]]]], [[]]]
        """
        if i == len(idx) - 1:
            self._setitem(idx[-1], value)
        else:
            with self[idx[i]].clone() as child:
                child.__setitem_rec__(idx, i + 1, value)
            self[idx[i]] = child

    def __getitem__(self, idx):
        """
        Return the ``idx``-th child of ``self`` (which is a subtree) if
        ``idx`` is an integer, or the ``idx[n-1]``-th child of the
        ``idx[n-2]``-th child of the ... of the ``idx[0]``-th child of
        ``self`` if ``idx`` is a list (or iterable) of length `n`.

        The indexing of the children is zero-based.

        INPUT:

        - ``idx`` -- integer; or a valid path in ``self`` identifying a node

        .. NOTE::

            The default implementation here assumes that the container of the
            node inherits from
            :class:`~sage.structure.list_clone.ClonableArray`.

        EXAMPLES::

            sage: x = OrderedTree([[],[[]]])
            sage: x[1,0]
            []
            sage: x = OrderedTree([[],[[]]])
            sage: x[()]
            [[], [[]]]
            sage: x[(0,)]
            []
            sage: x[0,0]
            Traceback (most recent call last):
            ...
            IndexError: list index out of range

            sage: u = BinaryTree(None)
            sage: v = BinaryTree([u, u])
            sage: w = BinaryTree([u, v])
            sage: t = BinaryTree([v, w])
            sage: z = BinaryTree([w, t])
            sage: z[0,1]
            [., .]
            sage: z[0,0]
            .
            sage: z[1]
            [[., .], [., [., .]]]
            sage: z[1,1]
            [., [., .]]
            sage: z[1][1,1]
            [., .]
        """
        if isinstance(idx, slice):
            return ClonableArray.__getitem__(self, idx)
        try:
            i = int(idx)
        except TypeError:
            res = self
            # idx is supposed to be an iterable of ints
            for i in idx:
                res = ClonableArray._getitem(res, i)
            return res
        else:
            return ClonableArray._getitem(self, i)


class AbstractLabelledTree(AbstractTree):
    """
    Abstract Labelled Tree.

    Typically a class for labelled trees is constructed by inheriting from
    a class for unlabelled trees and :class:`AbstractLabelledTree`.

    .. rubric:: How should this class be extended ?

    A class extending :class:`AbstractLabelledTree
    <sage.combinat.abstract_tree.AbstractLabelledTree>` should respect the
    following assumptions:

    * For a labelled tree ``T`` the call ``T.parent().unlabelled_trees()``
      should return a parent for unlabelled trees of the same kind: for
      example,

      - if ``T`` is a binary labelled tree, ``T.parent()`` is
        ``LabelledBinaryTrees()`` and ``T.parent().unlabelled_trees()`` is
        ``BinaryTrees()``

      - if ``T`` is an ordered labelled tree, ``T.parent()`` is
        ``LabelledOrderedTrees()`` and ``T.parent().unlabelled_trees()`` is
        ``OrderedTrees()``

    * In the same vein, the class of ``T`` should contain an attribute
      ``_UnLabelled`` which should be the class for the corresponding
      unlabelled trees.

    See also the assumptions in :class:`AbstractTree`.

    .. SEEALSO:: :class:`AbstractTree`
    """

    def __init__(self, parent, children, label=None, check=True):
        """
        TESTS::

            sage: LabelledOrderedTree([])
            None[]
            sage: LabelledOrderedTree([], 3)
            3[]
            sage: LT = LabelledOrderedTree
            sage: t = LT([LT([LT([], label=42), LT([], 21)])], label=1)
            sage: t
            1[None[42[], 21[]]]
            sage: LabelledOrderedTree(OrderedTree([[],[[],[]],[]]))
            None[None[], None[None[], None[]], None[]]

        We test that inheriting from `LabelledOrderedTree` allows construction from a
        `LabelledOrderedTree` (:issue:`16314`)::

            sage: LBTS = LabelledOrderedTrees()
            sage: class Foo(LabelledOrderedTree):
            ....:     def bar(self):
            ....:         print("bar called")
            sage: foo = Foo(LBTS, [], label=1); foo
            1[]
            sage: foo1 = LBTS([LBTS([], label=21)], label=42); foo1
            42[21[]]
            sage: foo2 = Foo(LBTS, foo1); foo2
            42[21[]]
            sage: foo2[0]
            21[]
            sage: foo2.__class__
            <class '__main__.Foo'>
            sage: foo2[0].__class__
            <class '__main__.Foo'>
            sage: foo2.bar()
            bar called
            sage: foo2.label()
            42
        """
        # We must initialize the label before the subtrees to allows rooted
        # trees canonization. Indeed it needs that ``self``._hash_() is working
        # at the end of the call super().__init__(...)
        if isinstance(children, AbstractLabelledTree):
            if label is None:
                self._label = children._label
            else:
                self._label = label
        else:
            self._label = label
        super().__init__(parent, children, check=check)

    def _repr_(self):
        """
        Return the string representation of ``self``.

        TESTS::

            sage: LabelledOrderedTree([])            # indirect doctest
            None[]
            sage: LabelledOrderedTree([], label=3)   # indirect doctest
            3[]
            sage: LabelledOrderedTree([[],[[]]])     # indirect doctest
            None[None[], None[None[]]]
            sage: LabelledOrderedTree([[],LabelledOrderedTree([[]], label=2)], label=3)
            3[None[], 2[None[]]]
        """
        return "%s%s" % (self._label, self[:])

    def label(self, path=None):
        """
        Return the label of ``self``.

        INPUT:

        - ``path`` -- ``None`` (default) or a path (list or tuple of
          children index in the tree)

        OUTPUT: the label of the subtree indexed by ``path``

        EXAMPLES::

            sage: t = LabelledOrderedTree([[],[]], label = 3)
            sage: t.label()
            3
            sage: t[0].label()
            sage: t = LabelledOrderedTree([LabelledOrderedTree([], 5),[]], label = 3)
            sage: t.label()
            3
            sage: t[0].label()
            5
            sage: t[1].label()
            sage: t.label([0])
            5
        """
        if path is None:
            return self._label
        else:
            tr = self
            for i in path:
                tr = tr[i]
            return tr._label

    def labels(self):
        """
        Return the list of labels of ``self``.

        EXAMPLES::

            sage: LT = LabelledOrderedTree
            sage: t = LT([LT([],label='b'),LT([],label='c')],label='a')
            sage: t.labels()
            ['a', 'b', 'c']

            sage: LBT = LabelledBinaryTree
            sage: LBT([LBT([],label=1),LBT([],label=4)],label=2).labels()
            [2, 1, 4]
        """
        return [t.label() for t in self.subtrees()]

    def leaf_labels(self):
        """
        Return the list of labels of the leaves of ``self``.

        In case of a labelled binary tree, these "leaves" are not actually
        the leaves of the binary trees, but the nodes whose both children
        are leaves!

        EXAMPLES::

            sage: LT = LabelledOrderedTree
            sage: t = LT([LT([],label='b'),LT([],label='c')],label='a')
            sage: t.leaf_labels()
            ['b', 'c']

            sage: LBT = LabelledBinaryTree
            sage: bt = LBT([LBT([],label='b'),LBT([],label='c')],label='a')
            sage: bt.leaf_labels()
            ['b', 'c']
            sage: LBT([], label='1').leaf_labels()
            ['1']
            sage: LBT(None).leaf_labels()
            []
        """
        return [t.label() for t in self.subtrees() if t.node_number() == 1]

    def __eq__(self, other):
        """
        Test if ``self`` is equal to ``other``.

        TESTS::

            sage  LabelledOrderedTree() == LabelledOrderedTree()
            True
            sage  LabelledOrderedTree([]) == LabelledOrderedTree()
            False
            sage: t1 = LabelledOrderedTree([[],[[]]])
            sage: t2 = LabelledOrderedTree([[],[[]]])
            sage: t1 == t2
            True
            sage: t2 = LabelledOrderedTree(t1)
            sage: t1 == t2
            True
            sage: t1 = LabelledOrderedTree([[],[[]]])
            sage: t2 = LabelledOrderedTree([[[]],[]])
            sage: t1 == t2
            False
        """
        return super().__eq__(other) and self._label == other._label

    def _hash_(self):
        """
        Return the hash value for ``self``.

        TESTS::

            sage: t1 = LabelledOrderedTree([[],[[]]], label = 1); t1hash = t1.__hash__()
            sage: LabelledOrderedTree([[],[[]]], label = 1).__hash__() == t1hash
            True
            sage: LabelledOrderedTree([[[]],[]], label = 1).__hash__() == t1hash
            False
            sage: LabelledOrderedTree(t1, label = 1).__hash__() == t1hash
            True
            sage: LabelledOrderedTree([[],[[]]], label = 25).__hash__() == t1hash
            False
            sage: LabelledOrderedTree(t1, label = 25).__hash__() == t1hash
            False

            sage: LabelledBinaryTree([[],[[],[]]], label = 25).__hash__() #random
            8544617749928727644

        We check that the hash value depends on the value of the labels of the
        subtrees::

            sage: LBT = LabelledBinaryTree
            sage: t1 = LBT([], label = 1)
            sage: t2 = LBT([], label = 2)
            sage: t3 = LBT([], label = 3)
            sage: t12 = LBT([t1, t2], label = "a")
            sage: t13 = LBT([t1, t3], label = "a")
            sage: t12.__hash__() != t13.__hash__()
            True
        """
        return self._UnLabelled._hash_(self) ^ hash(self._label)

    def shape(self):
        """
        Return the unlabelled tree associated to ``self``.

        EXAMPLES::

            sage: t = LabelledOrderedTree([[],[[]]], label = 25).shape(); t
            [[], [[]]]

            sage: LabelledBinaryTree([[],[[],[]]], label = 25).shape()
            [[., .], [[., .], [., .]]]

            sage: LRT = LabelledRootedTree
            sage: tb = LRT([],label='b')
            sage: LRT([tb, tb], label='a').shape()
            [[], []]

        TESTS::

            sage: t.parent()
            Ordered trees
            sage: type(t)
            <class 'sage.combinat.ordered_tree.OrderedTrees_all_with_category.element_class'>
        """
        TR = self.parent().unlabelled_trees()
        if not self:
            return TR.leaf()
        else:
            return TR._element_constructor_([i.shape() for i in self])

    def as_digraph(self):
        """
        Return a directed graph version of ``self``.

        .. WARNING::

            At this time, the output makes sense only if ``self`` is a
            labelled binary tree with no repeated labels and no ``None``
            labels.

        EXAMPLES::

           sage: LT = LabelledOrderedTrees()
           sage: t1 = LT([LT([],label=6),LT([],label=1)],label=9)
           sage: t1.as_digraph()
           Digraph on 3 vertices

           sage: t = BinaryTree([[None, None],[[],None]])
           sage: lt = t.canonical_labelling()
           sage: lt.as_digraph()
           Digraph on 4 vertices
        """
        from sage.graphs.digraph import DiGraph
        resu = {self.label():
                [t.label() for t in self if not t.is_empty()]}
        resu = DiGraph(resu, format='dict_of_lists')
        for t in self:
            if not t.is_empty():
                resu = resu.union(t.as_digraph())
        return resu


class AbstractLabelledClonableTree(AbstractLabelledTree,
                                   AbstractClonableTree):
    """
    Abstract Labelled Clonable Tree.

    This class takes care of modification for the label by the clone protocol.

    .. NOTE:: Due to the limitation of Cython inheritance, one cannot inherit
       here from :class:`~sage.structure.list_clone.ClonableArray`, because it would prevent us to
       inherit later from :class:`~sage.structure.list_clone.ClonableList`.
    """

    def set_root_label(self, label):
        """
        Set the label of the root of ``self``.

        INPUT:

        - ``label`` -- any Sage object

        OUTPUT: none, ``self`` is modified in place

        .. NOTE::

            ``self`` must be in a mutable state. See
            :mod:`sage.structure.list_clone` for more details about
            mutability.

        EXAMPLES::

            sage: t = LabelledOrderedTree([[],[[],[]]])
            sage: t.set_root_label(3)
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: with t.clone() as t:
            ....:     t.set_root_label(3)
            sage: t.label()
            3
            sage: t
            3[None[], None[None[], None[]]]

        This also works for binary trees::

            sage: bt = LabelledBinaryTree([[],[]])
            sage: bt.set_root_label(3)
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: with bt.clone() as bt:
            ....:     bt.set_root_label(3)
            sage: bt.label()
            3
            sage: bt
            3[None[., .], None[., .]]

        TESTS::

            sage: with t.clone() as t:
            ....:     t[0] = LabelledOrderedTree(t[0], label = 4)
            sage: t
            3[4[], None[None[], None[]]]
            sage: with t.clone() as t:
            ....:     t[1,0] = LabelledOrderedTree(t[1,0], label = 42)
            sage: t
            3[4[], None[42[], None[]]]
        """
        self._require_mutable()
        self._label = label

    def set_label(self, path, label):
        """
        Change the label of subtree indexed by ``path`` to ``label``.

        INPUT:

        - ``path`` -- ``None`` (default) or a path (list or tuple of children
          index in the tree)

        - ``label`` -- any sage object

        OUTPUT: nothing, ``self`` is modified in place

        .. NOTE::

            ``self`` must be in a mutable state. See
            :mod:`sage.structure.list_clone` for more details about
            mutability.

        EXAMPLES::

            sage: t = LabelledOrderedTree([[],[[],[]]])
            sage: t.set_label((0,), 4)
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: with t.clone() as t:
            ....:     t.set_label((0,), 4)
            sage: t
            None[4[], None[None[], None[]]]
            sage: with t.clone() as t:
            ....:     t.set_label((1,0), label = 42)
            sage: t
            None[4[], None[42[], None[]]]

        .. TODO::

            Do we want to implement the following syntactic sugar::

                with t.clone() as tt:
                    tt.labels[1,2] = 3 ?
        """
        self._require_mutable()
        path = tuple(path)
        if path == ():
            self._label = label
        else:
            with self[path[0]].clone() as child:
                child.set_label(path[1:], label)
            self[path[0]] = child

    def map_labels(self, f):
        """
        Apply the function `f` to the labels of ``self``.

        This method returns a copy of ``self`` on which the function `f` has
        been applied on all labels (a label `x` is replaced by `f(x)`).

        EXAMPLES::

            sage: LT = LabelledOrderedTree
            sage: t = LT([LT([],label=1),LT([],label=7)],label=3); t
            3[1[], 7[]]
            sage: t.map_labels(lambda z:z+1)
            4[2[], 8[]]

            sage: LBT = LabelledBinaryTree
            sage: bt = LBT([LBT([],label=1),LBT([],label=4)],label=2); bt
            2[1[., .], 4[., .]]
            sage: bt.map_labels(lambda z:z+1)
            3[2[., .], 5[., .]]
        """
        if self.is_empty():
            return self
        return self.parent()([t.map_labels(f) for t in self],
                             label=f(self.label()))


def from_hexacode(ch, parent=None, label='@'):
    r"""
    Transform a hexadecimal string into a tree.

    INPUT:

    - ``ch`` -- a hexadecimal string

    - ``parent`` -- kind of trees to be produced. If ``None``, this will
      be ``LabelledOrderedTrees``

    - ``label`` -- a label (default: ``'@'``) to be used for every vertex
      of the tree

    See :meth:`AbstractTree.to_hexacode` for the description of the encoding

    See :func:`_from_hexacode_aux` for the actual code

    EXAMPLES::

        sage: from sage.combinat.abstract_tree import from_hexacode
        sage: from_hexacode('12000', LabelledOrderedTrees())
        @[@[@[], @[]]]
        sage: from_hexacode('12000')
        @[@[@[], @[]]]

        sage: from_hexacode('1200', LabelledOrderedTrees())
        @[@[@[], @[]]]

    It can happen that only a prefix of the word is used::

        sage: from_hexacode('a'+14*'0', LabelledOrderedTrees())
        @[@[], @[], @[], @[], @[], @[], @[], @[], @[], @[]]

    One can choose the label::

        sage: from_hexacode('1200', LabelledOrderedTrees(), label='o')
        o[o[o[], o[]]]

    One can also create other kinds of trees::

        sage: from_hexacode('1200', OrderedTrees())
        [[[], []]]
    """
    if parent is None:
        from sage.combinat.ordered_tree import LabelledOrderedTrees
        parent = LabelledOrderedTrees()
    return _from_hexacode_aux(ch, parent, label)[0]


def _from_hexacode_aux(ch, parent, label='@'):
    r"""
    Transform a hexadecimal string into a tree and a remainder string.

    INPUT:

    - ``ch`` -- a hexadecimal string

    - ``parent`` -- kind of trees to be produced

    - ``label`` -- a label (default: ``'@'``) to be used for every vertex
      of the tree

    This method is used in :func:`from_hexacode`

    EXAMPLES::

        sage: from sage.combinat.abstract_tree import _from_hexacode_aux
        sage: _from_hexacode_aux('12000', LabelledOrderedTrees())
        (@[@[@[], @[]]], '0')

        sage: _from_hexacode_aux('1200', LabelledOrderedTrees())
        (@[@[@[], @[]]], '')

        sage: _from_hexacode_aux('1200', OrderedTrees())
        ([[[], []]], '')

        sage: _from_hexacode_aux('a00000000000000', LabelledOrderedTrees())
        (@[@[], @[], @[], @[], @[], @[], @[], @[], @[], @[]], '0000')
    """
    Trees = parent
    width = int(ch[0], 16)  # hexadecimal input
    remainder = ch[1:]
    if width == 0:
        return (Trees([], label), remainder)
    branches = {}
    for i in range(width):
        tree, remainder = _from_hexacode_aux(remainder, parent, label)
        branches[i] = tree
    return (Trees(branches.values(), label), remainder)
