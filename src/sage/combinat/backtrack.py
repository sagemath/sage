r"""
Backtracking

This library contains generic tools for constructing large sets whose
elements can be enumerated by exploring a search space with a (lazy)
tree or graph structure.

- :class:`GenericBacktracker`: Depth first search through a tree
  described by a ``children`` function, with branch pruning, etc.

Deprecated classes (use :func:`RecursivelyEnumeratedSet` instead):

- :class:`SearchForest`: Depth and breadth first
  search through a tree described by a ``children`` function.

- :class:`TransitiveIdeal`: Depth first search through a
  graph described by a ``neighbours`` relation.

- :class:`TransitiveIdealGraded`: Breadth first search
  through a graph described by a ``neighbours`` relation.

Deprecation details:

- ``SearchForest(seeds, succ)`` keeps the same behavior as before
  :trac:`6637` and is now the same as ``RecursivelyEnumeratedSet(seeds,
  succ, structure='forest', enumeration='depth')``.

- ``TransitiveIdeal(succ, seeds)`` keeps the same behavior as before
  :trac:`6637` and is now the same as ``RecursivelyEnumeratedSet(seeds,
  succ, structure=None, enumeration='naive')``.

- ``TransitiveIdealGraded(succ, seeds, max_depth)`` keeps the same behavior
  as before :trac:`6637` and is now the same as
  ``RecursivelyEnumeratedSet(seeds, succ, structure=None,
  enumeration='breadth', max_depth=max_depth)``.

.. todo::

    - For now the code of :class:`SearchForest` is still in
      ``sage/combinat/backtrack.py``.  It should be moved in
      ``sage/sets/recursively_enumerated_set.pyx`` into a class named
      :class:`RecursivelyEnumeratedSet_forest` in a later ticket.

    - Deprecate ``TransitiveIdeal`` and ``TransitiveIdealGraded``.

    - Once the deprecation has been there for enough time: delete
      ``TransitiveIdeal`` and ``TransitiveIdealGraded``.

"""
#*****************************************************************************
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>,
#                     2009 Nicolas M. Thiery <nthiery at users.sf.net>
#                     2010 Nicolas Borie <nicolas.borie at math.u-psud.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.monoids import Monoids
from sage.structure.parent import Parent
from sage.misc.prandom import randint
from sage.misc.abstract_method import abstract_method
from sage.categories.commutative_additive_semigroups import (
        CommutativeAdditiveSemigroups)
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer_ring import ZZ
from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_generic


def _imap_and_filter_none(function, iterable):
    r"""
    Return an iterator over the elements ``function(x)``, where ``x``
    iterates through ``iterable``, such that ``function(x)`` is not
    ``None``.

    EXAMPLES::

        sage: from sage.combinat.backtrack import _imap_and_filter_none
        sage: p = _imap_and_filter_none(lambda x: x if is_prime(x) else None, range(15))
        sage: [next(p), next(p), next(p), next(p), next(p), next(p)]
        [2, 3, 5, 7, 11, 13]
        sage: p = _imap_and_filter_none(lambda x: x+x, ['a','b','c','d','e'])
        sage: [next(p), next(p), next(p), next(p), next(p)]
        ['aa', 'bb', 'cc', 'dd', 'ee']
    """
    for x in iterable:
        x = function(x)
        if x is not None:
            yield x


class GenericBacktracker(object):
    r"""
    A generic backtrack tool for exploring a search space organized as a tree,
    with branch pruning, etc.

    See also :class:`SearchForest` and :class:`TransitiveIdeal` for
    handling simple special cases.
    """
    def __init__(self, initial_data, initial_state):
        r"""
        EXAMPLES::

            sage: from sage.combinat.backtrack import GenericBacktracker
            sage: p = GenericBacktracker([], 1)
            sage: loads(dumps(p))
            <sage.combinat.backtrack.GenericBacktracker object at 0x...>
        """
        self._initial_data = initial_data
        self._initial_state = initial_state

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.permutation import PatternAvoider
            sage: p = PatternAvoider(Permutations(4), [[1,3,2]])
            sage: len(list(p))
            14
        """
        #Initialize the stack of generators with the initial data.
        #The generator in stack[i] is a generator for the i^th level
        #of the search tree.
        stack = []
        stack.append(self._rec(self._initial_data, self._initial_state))

        done = False
        while not done:
            #Try to get the next object in this level
            try:
                obj, state, yld = next(stack[-1])
            except StopIteration:
                #If there are no more, go back up the tree
                #We also need to check if we've exhausted all
                #possibilities
                stack.pop()
                done = len(stack) == 0
                continue

            #If the return state is None, then obj is a leaf
            #of the search tree.  If yld is True, then obj
            #should be yielded.
            if yld is True:
                yield obj
            if state is not None:
                stack.append( self._rec(obj, state) )

def search_forest_iterator(roots, children, algorithm='depth'):
    r"""
    Return an iterator on the nodes of the forest having the given
    roots, and where ``children(x)`` returns the children of the node ``x``
    of the forest.  Note that every node of the tree is returned,
    not simply the leaves.

    INPUT:

    - ``roots`` -- a list (or iterable)
    - ``children`` -- a function returning a list (or iterable)
    - ``algorithm`` -- ``'depth'`` or ``'breadth'`` (default: ``'depth'``)

    EXAMPLES:

    We construct the prefix tree of binary sequences of length at most
    three, and enumerate its nodes::

        sage: from sage.combinat.backtrack import search_forest_iterator
        sage: list(search_forest_iterator([[]], lambda l: [l+[0], l+[1]]
        ....:                                   if len(l) < 3 else []))
        [[], [0], [0, 0], [0, 0, 0], [0, 0, 1], [0, 1], [0, 1, 0],
         [0, 1, 1], [1], [1, 0], [1, 0, 0], [1, 0, 1], [1, 1], [1, 1, 0], [1, 1, 1]]

    By default, the nodes are iterated through by depth first search.
    We can instead use a breadth first search (increasing depth)::

        sage: list(search_forest_iterator([[]], lambda l: [l+[0], l+[1]]
        ....:                                   if len(l) < 3 else [],
        ....:                             algorithm='breadth'))
        [[],
         [0], [1],
         [0, 0], [0, 1], [1, 0], [1, 1],
         [0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1],
         [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]]

    This allows for iterating trough trees of infinite depth::

        sage: it = search_forest_iterator([[]], lambda l: [l+[0], l+[1]], algorithm='breadth')
        sage: [ next(it) for i in range(16) ]
        [[],
         [0], [1], [0, 0], [0, 1], [1, 0], [1, 1],
         [0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1],
         [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1],
         [0, 0, 0, 0]]

    Here is an iterator through the prefix tree of sequences of
    letters in `0,1,2` without repetitions, sorted by length; the
    leaves are therefore permutations::

        sage: list(search_forest_iterator([[]], lambda l: [l + [i] for i in range(3) if i not in l],
        ....:                             algorithm='breadth'))
        [[],
         [0], [1], [2],
         [0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1],
         [0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]
    """
    # Little trick: the same implementation handles both depth and
    # breadth first search. Setting position to -1 makes a depth search
    # (you ask the children for the last node you met). Setting
    # position on 0 makes a breadth search (enumerate all the
    # descendants of a node before going on to the next father)
    if algorithm == 'depth':
        position = -1
    else:
        position = 0

    # Invariant:
    #  - for breadth first search: stack[i] contains an iterator over the nodes
    #    of depth ``i`` in the tree
    #  - for depth first search: stack[i] contains an iterator over the children
    #    of the node at depth ``i-1`` in the current branch (assuming a virtual
    #    father of all roots at depth ``-1``)
    stack = [iter(roots)]
    while stack:
        try:
            node = next(stack[position])
        except StopIteration:
            # If there are no more, go back up the tree
            # We also need to check if we've exhausted all
            # possibilities
            stack.pop(position)
            continue

        yield node
        stack.append( iter(children(node)) )

class SearchForest(Parent):
    r"""
    The enumerated set of the nodes of the forest having the given
    ``roots``, and where ``children(x)`` returns the children of the
    node ``x`` of the forest.

    See also :class:`GenericBacktracker`, :class:`TransitiveIdeal`,
    and :class:`TransitiveIdealGraded`.

    INPUT:

    - ``roots`` -- a list (or iterable)
    - ``children`` -- a function returning a list (or iterable, or iterator)
    - ``post_process`` -- a function defined over the nodes of the
      forest (default: no post processing)
    - ``algorithm`` -- ``'depth'`` or ``'breadth'`` (default: ``'depth'``)
    - ``category`` -- a category (default: :class:`EnumeratedSets`)

    The option ``post_process`` allows for customizing the nodes that
    are actually produced. Furthermore, if ``f(x)`` returns ``None``,
    then ``x`` won't be output at all.

    EXAMPLES:

    We construct the set of all binary sequences of length at most
    three, and list them::

        sage: from sage.combinat.backtrack import SearchForest
        sage: S = SearchForest( [[]],
        ....:     lambda l: [l+[0], l+[1]] if len(l) < 3 else [],
        ....:     category=FiniteEnumeratedSets())
        sage: S.list()
        [[],
         [0], [0, 0], [0, 0, 0], [0, 0, 1], [0, 1], [0, 1, 0], [0, 1, 1],
         [1], [1, 0], [1, 0, 0], [1, 0, 1], [1, 1], [1, 1, 0], [1, 1, 1]]

    ``SearchForest`` needs to be explicitly told that the set is
    finite for the following to work::

        sage: S.category()
        Category of finite enumerated sets
        sage: S.cardinality()
        15

    We proceed with the set of all lists of letters in ``0,1,2``
    without repetitions, ordered by increasing length (i.e. using a
    breadth first search through the tree)::

        sage: from sage.combinat.backtrack import SearchForest
        sage: tb = SearchForest( [[]],
        ....:       lambda l: [l + [i] for i in range(3) if i not in l],
        ....:       algorithm = 'breadth',
        ....:       category=FiniteEnumeratedSets())
        sage: tb[0]
        []
        sage: tb.cardinality()
        16
        sage: list(tb)
        [[],
         [0], [1], [2],
         [0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1],
         [0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]

    For infinite sets, this option should be set carefully to ensure
    that all elements are actually generated. The following example
    builds the set of all ordered pairs `(i,j)` of nonnegative
    integers such that `j\leq 1`::

        sage: from sage.combinat.backtrack import SearchForest
        sage: I = SearchForest([(0,0)],
        ....:                  lambda l: [(l[0]+1, l[1]), (l[0], 1)]
        ....:                            if l[1] == 0 else [(l[0], l[1]+1)])

    With a depth first search, only the elements of the form `(i,0)`
    are generated::

        sage: depth_search = I.depth_first_search_iterator()
        sage: [next(depth_search) for i in range(7)]
        [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0), (6, 0)]

    Using instead breadth first search gives the usual anti-diagonal
    iterator::

        sage: breadth_search = I.breadth_first_search_iterator()
        sage: [next(breadth_search) for i in range(15)]
        [(0, 0),
         (1, 0), (0, 1),
         (2, 0), (1, 1), (0, 2),
         (3, 0), (2, 1), (1, 2), (0, 3),
         (4, 0), (3, 1), (2, 2), (1, 3), (0, 4)]

    .. rubric:: Deriving subclasses

    The class of a parent `A` may derive from :class:`SearchForest` so
    that `A` can benefit from enumeration tools. As a running example,
    we consider the problem of enumerating integers whose binary
    expansion have at most three nonzero digits. For example, `3 =
    2^1 + 2^0` has two nonzero digits. `15 = 2^3 + 2^2 + 2^1 + 2^0`
    has four nonzero digits. In fact, `15` is the smallest integer
    which is not in the enumerated set.

    To achieve this, we use ``SearchForest`` to enumerate binary tuples
    with at most three nonzero digits, apply a post processing to
    recover the corresponding integers, and discard tuples finishing
    by zero.

    A first approach is to pass the ``roots`` and ``children``
    functions as arguments to :meth:`SearchForest.__init__`::

        sage: from sage.combinat.backtrack import SearchForest
        sage: class A(UniqueRepresentation, SearchForest):
        ....:     def __init__(self):
        ....:         SearchForest.__init__(self, [()],
        ....:             lambda x : [x+(0,), x+(1,)] if sum(x) < 3 else [],
        ....:             lambda x : sum(x[i]*2^i for i in range(len(x))) if sum(x) != 0 and x[-1] != 0 else None,
        ....:             algorithm = 'breadth',
        ....:             category=InfiniteEnumeratedSets())
        sage: MyForest = A(); MyForest
        An enumerated set with a forest structure
        sage: MyForest.category()
        Category of infinite enumerated sets
        sage: p = iter(MyForest)
        sage: [next(p) for i in range(30)]
        [1, 2, 3, 4, 6, 5, 7, 8, 12, 10, 14, 9, 13, 11, 16, 24, 20, 28, 18, 26, 22, 17, 25, 21, 19, 32, 48, 40, 56, 36]

    An alternative approach is to implement ``roots`` and ``children``
    as methods of the subclass (in fact they could also be attributes
    of `A`). Namely, ``A.roots()`` must return an iterable containing
    the enumeration generators, and ``A.children(x)`` must return an
    iterable over the children of `x`. Optionally, `A` can have a
    method or attribute such that ``A.post_process(x)`` returns the
    desired output for the node ``x`` of the tree::

        sage: from sage.combinat.backtrack import SearchForest
        sage: class A(UniqueRepresentation, SearchForest):
        ....:     def __init__(self):
        ....:         SearchForest.__init__(self, algorithm = 'breadth',
        ....:                               category=InfiniteEnumeratedSets())
        ....:
        ....:     def roots(self):
        ....:         return [()]
        ....:
        ....:     def children(self, x):
        ....:         if sum(x) < 3:
        ....:             return [x+(0,), x+(1,)]
        ....:         else:
        ....:             return []
        ....:
        ....:     def post_process(self, x):
        ....:         if sum(x) == 0 or x[-1] == 0:
        ....:             return None
        ....:         else:
        ....:             return sum(x[i]*2^i for i in range(len(x)))
        sage: MyForest = A(); MyForest
        An enumerated set with a forest structure
        sage: MyForest.category()
        Category of infinite enumerated sets
        sage: p = iter(MyForest)
        sage: [next(p) for i in range(30)]
        [1, 2, 3, 4, 6, 5, 7, 8, 12, 10, 14, 9, 13, 11, 16, 24, 20, 28, 18, 26, 22, 17, 25, 21, 19, 32, 48, 40, 56, 36]

    .. warning::

        A :class:`SearchForest` instance is picklable if and only if
        the input functions are themselves picklable. This excludes
        anonymous or interactively defined functions::

            sage: def children(x):
            ....:     return [x+1]
            sage: S = SearchForest( [1], children, category=InfiniteEnumeratedSets())
            sage: dumps(S)
            Traceback (most recent call last):
            ...
            PicklingError: Can't pickle <...function...>: attribute lookup ... failed

        Let us now fake ``children`` being defined in a Python module::

            sage: import __main__
            sage: __main__.children = children
            sage: S = SearchForest( [1], children, category=InfiniteEnumeratedSets())
            sage: loads(dumps(S))
            An enumerated set with a forest structure
    """
    def __init__(self, roots = None, children = None, post_process = None,
                 algorithm = 'depth', facade = None, category=None):
        r"""
        TESTS::

            sage: from sage.combinat.backtrack import SearchForest
            sage: S = SearchForest(NN, lambda x : [], lambda x: x^2 if x.is_prime() else None)
            sage: S.category()
            Category of enumerated sets
        """
        if roots is not None:
            self._roots = roots
        if children is not None:
            self.children = children
        if post_process is not None:
            self.post_process = post_process
        self._algorithm = algorithm
        Parent.__init__(self, facade = facade, category = EnumeratedSets().or_subcategory(category))

    __len__ = None

    def _repr_(self):
        r"""
        TESTS::

            sage: from sage.combinat.backtrack import SearchForest
            sage: SearchForest( [1], lambda x: [x+1])
            An enumerated set with a forest structure
        """
        return "An enumerated set with a forest structure"

    def roots(self):
        r"""
        Return an iterable over the roots of ``self``.

        EXAMPLES::

            sage: from sage.combinat.backtrack import SearchForest
            sage: I = SearchForest([(0,0)], lambda l: [(l[0]+1, l[1]), (l[0], 1)] if l[1] == 0 else [(l[0], l[1]+1)])
            sage: [i for i in I.roots()]
            [(0, 0)]
            sage: I = SearchForest([(0,0),(1,1)], lambda l: [(l[0]+1, l[1]), (l[0], 1)] if l[1] == 0 else [(l[0], l[1]+1)])
            sage: [i for i in I.roots()]
            [(0, 0), (1, 1)]
        """
        return self._roots

    @abstract_method
    def children(self, x):
        r"""
        Return the children of the element ``x``

        The result can be a list, an iterable, an iterator, or even a
        generator.

        EXAMPLES::

            sage: from sage.combinat.backtrack import SearchForest
            sage: I = SearchForest([(0,0)], lambda l: [(l[0]+1, l[1]), (l[0], 1)] if l[1] == 0 else [(l[0], l[1]+1)])
            sage: [i for i in I.children((0,0))]
            [(1, 0), (0, 1)]
            sage: [i for i in I.children((1,0))]
            [(2, 0), (1, 1)]
            sage: [i for i in I.children((1,1))]
            [(1, 2)]
            sage: [i for i in I.children((4,1))]
            [(4, 2)]
            sage: [i for i in I.children((4,0))]
            [(5, 0), (4, 1)]
        """

    def __iter__(self):
        r"""
        Return an iterator over the elements of ``self``.

        EXAMPLES::

            sage: from sage.combinat.backtrack import SearchForest
            sage: def children(l):
            ....:      return [l+[0], l+[1]]
            ....:
            sage: C = SearchForest(([],), children)
            sage: f = C.__iter__()
            sage: next(f)
            []
            sage: next(f)
            [0]
            sage: next(f)
            [0, 0]
        """
        iter = search_forest_iterator(self.roots(),
                                      self.children,
                                      algorithm = self._algorithm)
        if hasattr(self, "post_process"):
            iter = _imap_and_filter_none(self.post_process, iter)
        return iter

    def depth_first_search_iterator(self):
        r"""
        Return a depth first search iterator over the elements of ``self``

        EXAMPLES::

            sage: from sage.combinat.backtrack import SearchForest
            sage: f = SearchForest([[]],
            ....:                  lambda l: [l+[0], l+[1]] if len(l) < 3 else [])
            sage: list(f.depth_first_search_iterator())
            [[], [0], [0, 0], [0, 0, 0], [0, 0, 1], [0, 1], [0, 1, 0], [0, 1, 1], [1], [1, 0], [1, 0, 0], [1, 0, 1], [1, 1], [1, 1, 0], [1, 1, 1]]
        """
        return iter(self)

    def breadth_first_search_iterator(self):
        r"""
        Return a breadth first search iterator over the elements of ``self``

        EXAMPLES::

            sage: from sage.combinat.backtrack import SearchForest
            sage: f = SearchForest([[]],
            ....:                  lambda l: [l+[0], l+[1]] if len(l) < 3 else [])
            sage: list(f.breadth_first_search_iterator())
            [[], [0], [1], [0, 0], [0, 1], [1, 0], [1, 1], [0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]]
            sage: S = SearchForest([(0,0)],
            ....: lambda x : [(x[0], x[1]+1)] if x[1] != 0 else [(x[0]+1,0), (x[0],1)],
            ....: post_process = lambda x: x if ((is_prime(x[0]) and is_prime(x[1])) and ((x[0] - x[1]) == 2)) else None)
            sage: p = S.breadth_first_search_iterator()
            sage: [next(p), next(p), next(p), next(p), next(p), next(p), next(p)]
            [(5, 3), (7, 5), (13, 11), (19, 17), (31, 29), (43, 41), (61, 59)]
        """
        iter = search_forest_iterator(self.roots(), self.children, algorithm='breadth')
        if hasattr(self, "post_process"):
            iter = _imap_and_filter_none(self.post_process, iter)
        return iter

    def _elements_of_depth_iterator_rec(self, depth=0):
        r"""
        Return an iterator over the elements of ``self`` of given depth.
        An element of depth `n` can be obtained applying `n` times the
        children function from a root. This function is not affected
        by post processing.

        EXAMPLES::

            sage: from sage.combinat.backtrack import SearchForest
            sage: I = SearchForest([(0,0)], lambda l: [(l[0]+1, l[1]), (l[0], 1)] if l[1] == 0 else [(l[0], l[1]+1)])
            sage: list(I._elements_of_depth_iterator_rec(8))
            [(8, 0), (7, 1), (6, 2), (5, 3), (4, 4), (3, 5), (2, 6), (1, 7), (0, 8)]
            sage: I = SearchForest([[]], lambda l: [l+[0], l+[1]] if len(l) < 3 else [])
            sage: list(I._elements_of_depth_iterator_rec(0))
            [[]]
            sage: list(I._elements_of_depth_iterator_rec(1))
            [[0], [1]]
            sage: list(I._elements_of_depth_iterator_rec(2))
            [[0, 0], [0, 1], [1, 0], [1, 1]]
            sage: list(I._elements_of_depth_iterator_rec(3))
            [[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]]
            sage: list(I._elements_of_depth_iterator_rec(4))
            []
        """
        if depth == 0:
            for node in self.roots():
                yield node
        else:
            for father in self._elements_of_depth_iterator_rec(depth - 1):
                for node in self.children(father):
                    yield node

    def elements_of_depth_iterator(self, depth=0):
        r"""
        Return an iterator over the elements of ``self`` of given depth.
        An element of depth `n` can be obtained applying `n` times the
        children function from a root.

        EXAMPLES::

            sage: from sage.combinat.backtrack import SearchForest
            sage: S = SearchForest([(0,0)] ,
            ....:        lambda x : [(x[0], x[1]+1)] if x[1] != 0 else [(x[0]+1,0), (x[0],1)],
            ....:        post_process = lambda x: x if ((is_prime(x[0]) and is_prime(x[1]))
            ....:                                        and ((x[0] - x[1]) == 2)) else None)
            sage: p = S.elements_of_depth_iterator(8)
            sage: next(p)
            (5, 3)
            sage: S = SearchForest(NN, lambda x : [],
            ....:                      lambda x: x^2 if x.is_prime() else None)
            sage: p = S.elements_of_depth_iterator(0)
            sage: [next(p), next(p), next(p), next(p), next(p)]
            [4, 9, 25, 49, 121]
        """
        iter = self._elements_of_depth_iterator_rec(depth)
        if hasattr(self, "post_process"):
            iter = _imap_and_filter_none(self.post_process, iter)
        return iter

    def __contains__(self, elt):
        r"""
        Return ``True`` if ``elt`` is in ``self``.

        .. warning::

           This is achieved by iterating through the elements until
           ``elt`` is found. In particular, this method will never
           stop when ``elt`` is not in ``self`` and ``self`` is
           infinite.

        EXAMPLES::

            sage: from sage.combinat.backtrack import SearchForest
            sage: S = SearchForest( [[]], lambda l: [l+[0], l+[1]] if len(l) < 3 else [], category=FiniteEnumeratedSets())
            sage: [4] in S
            False
            sage: [1] in S
            True
            sage: [1,1,1,1] in S
            False
            sage: all(S.__contains__(i) for i in iter(S))
            True
            sage: S = SearchForest([1], lambda x: [x+1], category=InfiniteEnumeratedSets())
            sage: 1 in S
            True
            sage: 732 in S
            True
            sage: -1 in S # not tested : Will never stop

        The algorithm uses a random enumeration of the nodes of the
        forest. This choice was motivated by examples in which both
        depth first search and breadth first search failed. The
        following example enumerates all ordered pairs of nonnegative
        integers, starting from an infinite set of roots, where each
        roots has an infinite number of children::

            sage: from sage.combinat.backtrack import SearchForest
            sage: S = SearchForest(Family(NN, lambda x : (x, 0)),
            ....: lambda x : Family(PositiveIntegers(), lambda y : (x[0], y)) if x[1] == 0 else [])
            sage: p = S.depth_first_search_iterator()
            sage: [next(p), next(p), next(p), next(p), next(p), next(p), next(p)]
            [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6)]
            sage: p = S.breadth_first_search_iterator()
            sage: [next(p), next(p), next(p), next(p), next(p), next(p), next(p)]
            [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0), (6, 0)]
            sage: (0,0) in S
            True
            sage: (1,1) in S
            True
            sage: (10,10) in S
            True
            sage: (42,18) in S
            True

        We now consider the same set of all ordered pairs of
        nonnegative integers but constructed in a different way. There
        still are infinitely many roots, but each node has a single
        child. From each root starts an infinite branch of breadth
        `1`::

            sage: S = SearchForest(Family(NN, lambda x : (x, 0)) , lambda x : [(x[0], x[1]+1)])
            sage: p = S.depth_first_search_iterator()
            sage: [next(p), next(p), next(p), next(p), next(p), next(p), next(p)]
            [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6)]
            sage: p = S.breadth_first_search_iterator()
            sage: [next(p), next(p), next(p), next(p), next(p), next(p), next(p)]
            [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0), (6, 0)]
            sage: (0,0) in S
            True
            sage: (1,1) in S
            True
            sage: (10,10) in S
            True
            sage: (37,11) in S
            True
        """
        stack = [iter(self.roots())]
        while stack:
            position = randint(0,len(stack)-1)
            try:
                node = next(stack[position])
            except StopIteration:
                stack.pop(position)
                continue

            if node == elt:
                return True
            stack.append( iter(self.children(node)) )
        return False

    def map_reduce(self, map_function = None,
                   reduce_function = None,
                   reduce_init = None):
        r"""
        Apply a Map/Reduce algorithm on ``self``

        INPUT:

        - ``map_function`` -- a function from the element of ``self`` to some
          set with a reduce operation (e.g.: a monoid). The default value is
          the constant function ``1``.

        - ``reduce_function`` -- the reduce function (e.g.: the addition of a
          monoid). The default value is ``+``.

        - ``reduce_init`` -- the initialisation of the reduction (e.g.: the
          neutral element of the monoid). The default value is ``0``.

        .. note::

            the effect of the default values is to compute the cardinality
            of ``self``.

        EXAMPLES::

            sage: seeds = [([i],i, i) for i in range(1,10)]
            sage: def succ(t):
            ....:     list, sum, last = t
            ....:     return [(list + [i], sum + i, i) for i in range(1, last)]
            sage: F = RecursivelyEnumeratedSet(seeds, succ,
            ....:                       structure='forest', enumeration='depth')

            sage: y = var('y')
            sage: def map_function(t):
            ....:     li, sum, _ = t
            ....:     return y ^ sum
            sage: reduce_function = lambda x,y: x + y
            sage: F.map_reduce(map_function, reduce_function, 0)
            y^45 + y^44 + y^43 + 2*y^42 + 2*y^41 + 3*y^40 + 4*y^39 + 5*y^38 + 6*y^37 + 8*y^36 + 9*y^35 + 10*y^34 + 12*y^33 + 13*y^32 + 15*y^31 + 17*y^30 + 18*y^29 + 19*y^28 + 21*y^27 + 21*y^26 + 22*y^25 + 23*y^24 + 23*y^23 + 23*y^22 + 23*y^21 + 22*y^20 + 21*y^19 + 21*y^18 + 19*y^17 + 18*y^16 + 17*y^15 + 15*y^14 + 13*y^13 + 12*y^12 + 10*y^11 + 9*y^10 + 8*y^9 + 6*y^8 + 5*y^7 + 4*y^6 + 3*y^5 + 2*y^4 + 2*y^3 + y^2 + y

        Here is an example with the default values::

            sage: F.map_reduce()
            511

        .. SEEALSO:: :mod:`sage.parallel.map_reduce`
        """
        import sage.parallel.map_reduce
        return sage.parallel.map_reduce.RESetMapReduce(
            forest = self,
            map_function = map_function,
            reduce_function = reduce_function,
            reduce_init = reduce_init).run()


class PositiveIntegerSemigroup(UniqueRepresentation, SearchForest):
    r"""
    The commutative additive semigroup of positive integers.

    This class provides an example of algebraic structure which
    inherits from :class:`SearchForest`. It builds the positive
    integers a la Peano, and endows it with its natural commutative
    additive semigroup structure.

    EXAMPLES::

        sage: from sage.combinat.backtrack import PositiveIntegerSemigroup
        sage: PP = PositiveIntegerSemigroup()
        sage: PP.category()
        Join of Category of monoids and Category of commutative additive semigroups and Category of infinite enumerated sets and Category of facade sets
        sage: PP.cardinality()
        +Infinity
        sage: PP.one()
        1
        sage: PP.an_element()
        1
        sage: some_elements = list(PP.some_elements()); some_elements
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100]

    TESTS::

        sage: from sage.combinat.backtrack import PositiveIntegerSemigroup
        sage: PP = PositiveIntegerSemigroup()

    We factor out the long test from the ``TestSuite``::

        sage: TestSuite(PP).run(skip='_test_enumerated_set_contains')
        sage: PP._test_enumerated_set_contains()  # long time
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.combinat.backtrack import PositiveIntegerSemigroup
            sage: PP = PositiveIntegerSemigroup()
        """
        SearchForest.__init__(self, facade = ZZ, category=(InfiniteEnumeratedSets(), CommutativeAdditiveSemigroups(), Monoids()))

    def roots(self):
        r"""
        Return the single root of ``self``.

        EXAMPLES::

            sage: from sage.combinat.backtrack import PositiveIntegerSemigroup
            sage: PP = PositiveIntegerSemigroup()
            sage: list(PP.roots())
            [1]
        """
        return [ZZ(1)]

    def children(self, x):
        r"""
        Return the single child ``x+1`` of the integer ``x``

        EXAMPLES::

            sage: from sage.combinat.backtrack import PositiveIntegerSemigroup
            sage: PP = PositiveIntegerSemigroup()
            sage: list(PP.children(1))
            [2]
            sage: list(PP.children(42))
            [43]
        """
        return [ZZ(x+1)]

    def one(self):
        r"""
        Return the unit of ``self``.

        EXAMPLES::

            sage: from sage.combinat.backtrack import PositiveIntegerSemigroup
            sage: PP = PositiveIntegerSemigroup()
            sage: PP.one()
            1
        """
        return self.first()

class TransitiveIdeal(RecursivelyEnumeratedSet_generic):
    r"""
    Generic tool for constructing ideals of a relation.

    INPUT:

    - ``relation`` -- a function (or callable) returning a list (or iterable)
    - ``generators`` -- a list (or iterable)

    Returns the set `S` of elements that can be obtained by repeated
    application of ``relation`` on the elements of ``generators``.

    Consider ``relation`` as modeling a directed graph (possibly with
    loops, cycles, or circuits). Then `S` is the ideal generated by
    ``generators`` under this relation.

    Enumerating the elements of `S` is achieved by depth first search
    through the graph. The time complexity is `O(n+m)` where `n` is
    the size of the ideal, and `m` the number of edges in the
    relation. The memory complexity is the depth, that is the maximal
    distance between a generator and an element of `S`.

    See also :class:`SearchForest` and :class:`TransitiveIdealGraded`.

    EXAMPLES::

        sage: from sage.combinat.backtrack import TransitiveIdeal
        sage: [i for i in TransitiveIdeal(lambda i: [i+1] if i<10 else [], [0])]
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        sage: [i for i in TransitiveIdeal(lambda i: [mod(i+1,3)], [0])]
        [0, 1, 2]
        sage: [i for i in TransitiveIdeal(lambda i: [mod(i+2,3)], [0])]
        [0, 2, 1]
        sage: [i for i in TransitiveIdeal(lambda i: [mod(i+2,10)], [0])]
        [0, 2, 4, 6, 8]
        sage: sorted(i for i in TransitiveIdeal(lambda i: [mod(i+3,10),mod(i+5,10)], [0]))
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: sorted(i for i in TransitiveIdeal(lambda i: [mod(i+4,10),mod(i+6,10)], [0]))
        [0, 2, 4, 6, 8]
        sage: [i for i in TransitiveIdeal(lambda i: [mod(i+3,9)], [0,1])]
        [0, 1, 3, 4, 6, 7]

        sage: [p for p in TransitiveIdeal(lambda x:[x],[Permutation([3,1,2,4]), Permutation([2,1,3,4])])]
        [[2, 1, 3, 4], [3, 1, 2, 4]]

    We now illustrate that the enumeration is done lazily, by depth first
    search::

        sage: C = TransitiveIdeal(lambda x: [x-1, x+1], (-10, 0, 10))
        sage: f = C.__iter__()
        sage: [ next(f) for i in range(6) ]
        [0, 1, 2, 3, 4, 5]

    We compute all the permutations of 3::

        sage: sorted(p for p in TransitiveIdeal(attrcall("permutohedron_succ"), [Permutation([1,2,3])]))
        [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]

    We compute all the permutations which are larger than [3,1,2,4],
    [2,1,3,4] in the right permutohedron::

        sage: sorted(p for p in TransitiveIdeal(attrcall("permutohedron_succ"),
        ....:     [Permutation([3,1,2,4]), Permutation([2,1,3,4])]))
        [[2, 1, 3, 4], [2, 1, 4, 3], [2, 3, 1, 4], [2, 3, 4, 1],
         [2, 4, 1, 3], [2, 4, 3, 1], [3, 1, 2, 4], [3, 1, 4, 2],
         [3, 2, 1, 4], [3, 2, 4, 1], [3, 4, 1, 2], [3, 4, 2, 1],
         [4, 2, 1, 3], [4, 2, 3, 1], [4, 3, 1, 2], [4, 3, 2, 1]]

    Using TransitiveIdeal people have been using the ``__contains__``
    method provided from the ``__iter__`` method. We need to make sure that
    this continues to work::

        sage: T = TransitiveIdeal(lambda a:[a+7,a+5], [0])
        sage: 12 in T
        True

    """
    def __init__(self, succ, generators):
        r"""
        TESTS::

            sage: from sage.combinat.backtrack import TransitiveIdeal
            sage: C = TransitiveIdeal(factor, (1, 2, 3))
            sage: C._succ
            <function factor at ...>
            sage: C._generators
            (1, 2, 3)
            sage: loads(dumps(C))   # should test for equality with C, but equality is not implemented
        """
        RecursivelyEnumeratedSet_generic.__init__(self, seeds=generators, successors=succ, enumeration='naive')
        self._generators = self._seeds
        self._succ = self.successors

    def __iter__(self):
        r"""
        Return an iterator on the elements of ``self``.

        TESTS::

            sage: from sage.combinat.backtrack import TransitiveIdeal
            sage: C = TransitiveIdeal(lambda x: [1,2], ())
            sage: list(C) # indirect doctest
            []

            sage: C = TransitiveIdeal(lambda x: [1,2], (1,))
            sage: list(C) # indirect doctest
            [1, 2]

            sage: C = TransitiveIdeal(lambda x: [], (1,2))
            sage: list(C) # indirect doctest
            [1, 2]

        """
        return self.naive_search_iterator()

class TransitiveIdealGraded(RecursivelyEnumeratedSet_generic):
    r"""
    Generic tool for constructing ideals of a relation.

    INPUT:

    - ``relation`` -- a function (or callable) returning a list (or iterable)

    - ``generators`` -- a list (or iterable)

    - ``max_depth`` -- (Default: infinity) Specifies the maximal depth to
      which elements are computed

    Return the set `S` of elements that can be obtained by repeated
    application of ``relation`` on the elements of ``generators``.

    Consider ``relation`` as modeling a directed graph (possibly with
    loops, cycles, or circuits). Then `S` is the ideal generated by
    ``generators`` under this relation.

    Enumerating the elements of `S` is achieved by breadth first search
    through the graph; hence elements are enumerated by increasing
    distance from the generators. The time complexity is `O(n+m)`
    where `n` is the size of the ideal, and `m` the number of edges in
    the relation. The memory complexity is the depth, that is the
    maximal distance between a generator and an element of `S`.

    See also :class:`SearchForest` and :class:`TransitiveIdeal`.

    EXAMPLES::

        sage: from sage.combinat.backtrack import TransitiveIdealGraded
        sage: [i for i in TransitiveIdealGraded(lambda i: [i+1] if i<10 else [], [0])]
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    We now illustrate that the enumeration is done lazily, by breadth first search::

        sage: C = TransitiveIdealGraded(lambda x: [x-1, x+1], (-10, 0, 10))
        sage: f = C.__iter__()

    The elements at distance 0 from the generators::

        sage: sorted([ next(f) for i in range(3) ])
        [-10, 0, 10]

    The elements at distance 1 from the generators::

        sage: sorted([ next(f) for i in range(6) ])
        [-11, -9, -1, 1, 9, 11]

    The elements at distance 2 from the generators::

        sage: sorted([ next(f) for i in range(6) ])
        [-12, -8, -2, 2, 8, 12]

    The enumeration order between elements at the same distance is not specified.

    We compute all the permutations which are larger than [3,1,2,4] or
    [2,1,3,4] in the permutohedron::

        sage: sorted(p for p in TransitiveIdealGraded(attrcall("permutohedron_succ"),
        ....:     [Permutation([3,1,2,4]), Permutation([2,1,3,4])]))
        [[2, 1, 3, 4], [2, 1, 4, 3], [2, 3, 1, 4], [2, 3, 4, 1],
         [2, 4, 1, 3], [2, 4, 3, 1], [3, 1, 2, 4], [3, 1, 4, 2],
         [3, 2, 1, 4], [3, 2, 4, 1], [3, 4, 1, 2], [3, 4, 2, 1],
         [4, 2, 1, 3], [4, 2, 3, 1], [4, 3, 1, 2], [4, 3, 2, 1]]
    """
    def __init__(self, succ, generators, max_depth=float("inf")):
        r"""
        TESTS::

            sage: from sage.combinat.backtrack import TransitiveIdealGraded
            sage: C = TransitiveIdealGraded(factor, (1, 2, 3))
            sage: C._succ
            <function factor at ...>
            sage: C._generators
            (1, 2, 3)
            sage: loads(dumps(C))   # should test for equality with C, but equality is not implemented
        """
        RecursivelyEnumeratedSet_generic.__init__(self, seeds=generators, successors=succ, enumeration='breadth', max_depth=max_depth)
        self._generators = self._seeds
        self._succ = self.successors

    def __iter__(self):
        r"""
        Return an iterator on the elements of ``self``.

        TESTS::

            sage: from sage.combinat.backtrack import TransitiveIdealGraded
            sage: C = TransitiveIdealGraded(lambda x: [1,2], ())
            sage: list(C) # indirect doctest
            []

            sage: C = TransitiveIdealGraded(lambda x: [1,2], (1,))
            sage: list(C) # indirect doctest
            [1, 2]

            sage: C = TransitiveIdealGraded(lambda x: [], (1,2))
            sage: list(C) # indirect doctest
            [1, 2]

        ::

            sage: fn = lambda i: [i+1] if i<10 else []
            sage: C = TransitiveIdealGraded(fn, [0], max_depth=1)
            sage: list(C)
            [0, 1]
        """
        return self.breadth_first_search_iterator()

