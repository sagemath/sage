r"""
Recursively Enumerated Sets

A set `S` is called recursively enumerable if there is an algorithm that
enumerates the members of `S`. We consider here the recursively enumerated
sets that are described by some ``seeds`` and a successor function
``successors``.  The successor function may have some structure (symmetric,
graded, forest) or not. The elements of a set having a symmetric, graded or
forest structure can be enumerated uniquely without keeping all of them in
memory. Many kinds of iterators are provided in this module: depth first
search, breadth first search and elements of given depth.

See :wikipedia:`Recursively_enumerable_set`.

See the documentation of :func:`RecursivelyEnumeratedSet` below for the
description of the inputs.

AUTHORS:

- Sébastien Labbé, April 2014, at Sage Days 57, Cernay-la-ville

EXAMPLES:

.. RUBRIC:: No hypothesis on the structure

What we mean by "no hypothesis" is that the set is not known
to be a forest, symmetric, or graded. However, it may have other
structure, such as not containing an oriented cycle, that does not
help with the enumeration.

In this example, the seed is 0 and the successor function is either ``+2``
or ``+3``. This is the set of nonnegative linear combinations of 2 and 3::

    sage: succ = lambda a:[a+2,a+3]
    sage: C = RecursivelyEnumeratedSet([0], succ)
    sage: C
    A recursively enumerated set (breadth first search)

Breadth first search::

    sage: it = C.breadth_first_search_iterator()
    sage: [next(it) for _ in range(10)]
    [0, 2, 3, 4, 5, 6, 7, 8, 9, 10]

Depth first search::

    sage: it = C.depth_first_search_iterator()
    sage: [next(it) for _ in range(10)]
    [0, 3, 6, 9, 12, 15, 18, 21, 24, 27]

.. RUBRIC:: Symmetric structure

The origin ``(0, 0)`` as seed and the upper, lower, left and right lattice
point as successor function. This function is symmetric since `p` is a
successor of `q` if and only if `q` is a successor or `p`::

    sage: succ = lambda a: [(a[0]-1,a[1]), (a[0],a[1]-1), (a[0]+1,a[1]), (a[0],a[1]+1)]
    sage: seeds = [(0,0)]
    sage: C = RecursivelyEnumeratedSet(seeds, succ, structure='symmetric', enumeration='depth')
    sage: C
    A recursively enumerated set with a symmetric structure (depth first search)

In this case, depth first search is the default enumeration for iteration::

    sage: it_depth = iter(C)
    sage: [next(it_depth) for _ in range(10)]
    [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9)]

Breadth first search::

    sage: it_breadth = C.breadth_first_search_iterator()
    sage: [next(it_breadth) for _ in range(13)]
    [(0, 0),
     (-1, 0), (0, -1), (1, 0), (0, 1),
     (-2, 0), (-1, -1), (-1, 1), (0, -2), (1, -1), (2, 0), (1, 1), (0, 2)]

Levels (elements of given depth)::

    sage: sorted(C.graded_component(0))
    [(0, 0)]
    sage: sorted(C.graded_component(1))
    [(-1, 0), (0, -1), (0, 1), (1, 0)]
    sage: sorted(C.graded_component(2))
    [(-2, 0), (-1, -1), (-1, 1), (0, -2), (0, 2), (1, -1), (1, 1), (2, 0)]

.. RUBRIC:: Graded structure

Identity permutation as seed and ``permutohedron_succ`` as successor
function::

    sage: succ = attrcall("permutohedron_succ")
    sage: seed = [Permutation([1..5])]
    sage: R = RecursivelyEnumeratedSet(seed, succ, structure='graded')
    sage: R
    A recursively enumerated set with a graded structure (breadth first search)

Depth first search iterator::

    sage: it_depth = R.depth_first_search_iterator()
    sage: [next(it_depth) for _ in range(5)]
    [[1, 2, 3, 4, 5],
     [1, 2, 3, 5, 4],
     [1, 2, 5, 3, 4],
     [1, 2, 5, 4, 3],
     [1, 5, 2, 4, 3]]

Breadth first search iterator::

    sage: it_breadth = R.breadth_first_search_iterator()
    sage: [next(it_breadth) for _ in range(5)]
    [[1, 2, 3, 4, 5],
     [2, 1, 3, 4, 5],
     [1, 3, 2, 4, 5],
     [1, 2, 4, 3, 5],
     [1, 2, 3, 5, 4]]

Elements of given depth iterator::

    sage: sorted(R.elements_of_depth_iterator(9))
    [[4, 5, 3, 2, 1], [5, 3, 4, 2, 1], [5, 4, 2, 3, 1], [5, 4, 3, 1, 2]]
    sage: list(R.elements_of_depth_iterator(10))
    [[5, 4, 3, 2, 1]]

Graded components (set of elements of the same depth)::

    sage: # needs sage.combinat
    sage: sorted(R.graded_component(0))
    [[1, 2, 3, 4, 5]]
    sage: sorted(R.graded_component(1))
    [[1, 2, 3, 5, 4], [1, 2, 4, 3, 5], [1, 3, 2, 4, 5], [2, 1, 3, 4, 5]]
    sage: sorted(R.graded_component(9))
    [[4, 5, 3, 2, 1], [5, 3, 4, 2, 1], [5, 4, 2, 3, 1], [5, 4, 3, 1, 2]]
    sage: sorted(R.graded_component(10))
    [[5, 4, 3, 2, 1]]

.. RUBRIC:: Forest structure (Example 1)

The set of words over the alphabet `\{a,b\}` can be generated from the
empty word by appending the letter `a` or `b` as a successor function. This set
has a forest structure::

    sage: seeds = ['']
    sage: succ = lambda w: [w+'a', w+'b']
    sage: C = RecursivelyEnumeratedSet(seeds, succ, structure='forest')
    sage: C
    An enumerated set with a forest structure

Depth first search iterator::

    sage: it = C.depth_first_search_iterator()
    sage: [next(it) for _ in range(6)]
    ['', 'a', 'aa', 'aaa', 'aaaa', 'aaaaa']

Breadth first search iterator::

    sage: it = C.breadth_first_search_iterator()
    sage: [next(it) for _ in range(6)]
    ['', 'a', 'b', 'aa', 'ab', 'ba']

This example was provided by Florent Hivert.

How to define a set using those classes?

Only two things are necessary to define a set using a
:class:`RecursivelyEnumeratedSet` object (the other
classes being very similar):

.. MATH::

    \begin{array}{ccc}
                      & \emptyset \\
        \hfil\swarrow & \downarrow & \searrow\hfil\\
                    a & b & c \\
        \begin{array}{ccc}
             \swarrow & \downarrow & \searrow \\
                   aa & ab & ac \\
        \end{array}   &
        \begin{array}{ccc}
             \swarrow & \downarrow & \searrow \\
                   ba & bb & bc \\
        \end{array}   &
        \begin{array}{ccc}
             \swarrow & \downarrow & \searrow \\
                   ca & cb & cc \\
        \end{array}
    \end{array}

For the previous example, the two necessary pieces of information are:

- the initial element ``""``;

- the function::

      lambda x: [x + letter for letter in ['a', 'b', 'c']

This would actually describe an **infinite** set, as such rules describe
"all words" on 3 letters. Hence, it is a good idea to replace the function by::

    lambda x: [x + letter for letter in ['a', 'b', 'c']] if len(x) < 2 else []

or even::

    sage: def children(x):
    ....:     if len(x) < 2:
    ....:         for letter in ['a', 'b', 'c']:
    ....:             yield x+letter

We can then create the :class:`RecursivelyEnumeratedSet` object with either::

    sage: S = RecursivelyEnumeratedSet([''],
    ....:     lambda x: [x + letter for letter in ['a', 'b', 'c']]
    ....:               if len(x) < 2 else [],
    ....:     structure='forest', enumeration='depth',
    ....:     category=FiniteEnumeratedSets())
    sage: S.list()
    ['', 'a', 'aa', 'ab', 'ac', 'b', 'ba', 'bb', 'bc', 'c', 'ca', 'cb', 'cc']

or::

    sage: S = RecursivelyEnumeratedSet([''], children,
    ....:     structure='forest', enumeration='depth',
    ....:     category=FiniteEnumeratedSets())
    sage: S.list()
    ['', 'a', 'aa', 'ab', 'ac', 'b', 'ba', 'bb', 'bc', 'c', 'ca', 'cb', 'cc']

.. RUBRIC:: Forest structure (Example 2)

This example was provided by Florent Hivert.

Here is a little more involved example. We want to iterate through all
permutations of a given set `S`. One solution is to take elements of `S` one
by one and insert them at every position. So a node of the generating tree
contains two pieces of information:

- the list ``lst`` of already inserted element;
- the set ``st`` of the yet to be inserted element.

We want to generate a permutation only if ``st`` is empty (leaves on the
tree). Also suppose for the sake of the example, that instead of list we want
to generate tuples. This selection of some nodes and final mapping of a
function to the element is done by the ``post_process = f`` argument. The
convention is that the generated elements are the ``s := f(n)``, except when
``s`` not ``None`` when no element is generated at all. Here is the code::

    sage: def children(node):
    ....:     (lst, st) = node
    ....:     st = set(st) # make a copy
    ....:     if st:
    ....:        el = st.pop()
    ....:        for i in range(len(lst) + 1):
    ....:            yield (lst[0:i] + [el] + lst[i:], st)
    sage: list(children(([1,2], {3,7,9})))
    [([9, 1, 2], {3, 7}), ([1, 9, 2], {3, 7}), ([1, 2, 9], {3, 7})]
    sage: def post_process(node):
    ....:     (l, s) = node
    ....:     return tuple(l) if not s else None
    sage: S = RecursivelyEnumeratedSet( [([], {1,3,6,8})],
    ....:     children, post_process=post_process,
    ....:     structure='forest', enumeration='depth',
    ....:     category=FiniteEnumeratedSets())
    sage: S.list()
    [(6, 3, 1, 8), (3, 6, 1, 8), (3, 1, 6, 8), (3, 1, 8, 6), (6, 1, 3, 8),
     (1, 6, 3, 8), (1, 3, 6, 8), (1, 3, 8, 6), (6, 1, 8, 3), (1, 6, 8, 3),
     (1, 8, 6, 3), (1, 8, 3, 6), (6, 3, 8, 1), (3, 6, 8, 1), (3, 8, 6, 1),
     (3, 8, 1, 6), (6, 8, 3, 1), (8, 6, 3, 1), (8, 3, 6, 1), (8, 3, 1, 6),
     (6, 8, 1, 3), (8, 6, 1, 3), (8, 1, 6, 3), (8, 1, 3, 6)]
    sage: S.cardinality()
    24
"""

# ****************************************************************************
#       Copyright (C) 2014 Sebastien Labbe <slabqc at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.parent cimport Parent
from sage.categories.enumerated_sets import EnumeratedSets
from sage.misc.abstract_method import abstract_method
from sage.misc.prandom import randint
from collections import deque


def RecursivelyEnumeratedSet(seeds, successors, structure=None,
            enumeration=None, max_depth=float("inf"), post_process=None,
            facade=None, category=None):
    r"""
    Return a recursively enumerated set.

    A set `S` is called recursively enumerable if there is an algorithm that
    enumerates the members of `S`. We consider here the recursively
    enumerated sets that are described by some ``seeds`` and a successor
    function ``successors``.

    Let `U` be a set and ``successors`` `:U \to 2^U` be a successor function
    associating to each element of `U` a subset of `U`. Let ``seeds`` be a
    subset of `U`. Let `S\subseteq U` be the set of elements of `U` that
    can be reached from a seed by applying recursively the ``successors``
    function. This class provides different kinds of iterators (breadth first,
    depth first, elements of given depth, etc.) for the elements of `S`.

    See :wikipedia:`Recursively_enumerable_set`.

    INPUT:

    - ``seeds`` -- list (or iterable) of hashable objects
    - ``successors`` -- function (or callable) returning a list (or iterable) of
      hashable objects
    - ``structure`` -- string (default: ``None``); structure of the
      set, possible values are:

      - ``None`` -- nothing is known about the structure of the set
      - ``'forest'`` -- if the ``successors`` function generates a *forest*, that
        is, each element can be reached uniquely from a seed
      - ``'graded'`` -- if the ``successors`` function is *graded*, that is, all
        paths from a seed to a given element have equal length
      - ``'symmetric'`` -- if the relation is *symmetric*, that is,
        ``y in successors(x)`` if and only if ``x in successors(y)``

    - ``enumeration`` -- ``'depth'``, ``'breadth'``, ``'naive'`` or ``None``
      (default: ``None``); the default enumeration for the
      ``__iter__`` function
    - ``max_depth`` -- integer (default: ``float("inf")``); limit
      the search to a certain depth, currently works only for breadth first
      search
    - ``post_process`` -- (default: ``None``) for forest only
    - ``facade`` -- (default: ``None``)
    - ``category`` -- (default: ``None``)

    EXAMPLES:

    A recursive set with no other information::

        sage: f = lambda a: [a+3, a+5]
        sage: C = RecursivelyEnumeratedSet([0], f)
        sage: C
        A recursively enumerated set (breadth first search)
        sage: it = iter(C)
        sage: [next(it) for _ in range(10)]
        [0, 3, 5, 6, 8, 10, 9, 11, 13, 15]

    A recursive set with a forest structure::

        sage: f = lambda a: [2*a,2*a+1]
        sage: C = RecursivelyEnumeratedSet([1], f, structure='forest'); C
        An enumerated set with a forest structure
        sage: it = C.depth_first_search_iterator()
        sage: [next(it) for _ in range(7)]
        [1, 2, 4, 8, 16, 32, 64]
        sage: it = C.breadth_first_search_iterator()
        sage: [next(it) for _ in range(7)]
        [1, 2, 3, 4, 5, 6, 7]

    A recursive set given by a symmetric relation::

        sage: f = lambda a: [a-1,a+1]
        sage: C = RecursivelyEnumeratedSet([10, 15], f, structure='symmetric')
        sage: C
        A recursively enumerated set with a symmetric structure (breadth first search)
        sage: it = iter(C)
        sage: [next(it) for _ in range(7)]
        [10, 15, 9, 11, 14, 16, 8]

    A recursive set given by a graded relation::

        sage: # needs sage.symbolic
        sage: def f(a):
        ....:     return [a + 1, a + I]
        sage: C = RecursivelyEnumeratedSet([0], f, structure='graded'); C
        A recursively enumerated set with a graded structure (breadth first search)
        sage: it = iter(C)
        sage: [next(it) for _ in range(7)]
        [0, 1, I, 2, I + 1, 2*I, 3]

    .. WARNING::

        If you do not set a good structure, you might obtain bad results,
        like elements generated twice::

            sage: f = lambda a: [a-1,a+1]
            sage: C = RecursivelyEnumeratedSet([0], f, structure='graded')
            sage: it = iter(C)
            sage: [next(it) for _ in range(7)]
            [0, -1, 1, -2, 0, 2, -3]

    TESTS:

    The successors method is an attribute::

        sage: R = RecursivelyEnumeratedSet([1], lambda x: [x+1, x-1])
        sage: R.successors(4)
        [5, 3]

    ::

        sage: C = RecursivelyEnumeratedSet((1, 2, 3), factor)
        sage: C.successors
        <function factor at ...>
        sage: C._seeds
        (1, 2, 3)
    """
    if structure is None:
        if enumeration is None:
            enumeration = 'breadth'
        return RecursivelyEnumeratedSet_generic(seeds, successors,
                enumeration, max_depth, facade=facade, category=category)
    if structure == 'symmetric':
        if enumeration is None:
            enumeration = 'breadth'
        return RecursivelyEnumeratedSet_symmetric(seeds, successors,
                enumeration, max_depth, facade=facade, category=category)
    if structure == 'forest':
        if enumeration is None:
            enumeration = 'depth'
        return RecursivelyEnumeratedSet_forest(roots=seeds, children=successors,
                algorithm=enumeration, post_process=post_process,
                facade=facade, category=category)
    if structure == 'graded':
        if enumeration is None:
            enumeration = 'breadth'
        return RecursivelyEnumeratedSet_graded(seeds, successors, enumeration,
                max_depth, facade=facade, category=category)

    raise ValueError("Unknown value for structure (={})".format(structure))


cdef class RecursivelyEnumeratedSet_generic(Parent):
    r"""
    A generic recursively enumerated set.

    For more information, see :func:`RecursivelyEnumeratedSet`.

    EXAMPLES::

        sage: f = lambda a:[a+1]

    Different structure for the sets::

        sage: RecursivelyEnumeratedSet([0], f, structure=None)
        A recursively enumerated set (breadth first search)
        sage: RecursivelyEnumeratedSet([0], f, structure='graded')
        A recursively enumerated set with a graded structure (breadth first search)
        sage: RecursivelyEnumeratedSet([0], f, structure='symmetric')
        A recursively enumerated set with a symmetric structure (breadth first search)
        sage: RecursivelyEnumeratedSet([0], f, structure='forest')
        An enumerated set with a forest structure

    Different default enumeration algorithms::

        sage: RecursivelyEnumeratedSet([0], f, enumeration='breadth')
        A recursively enumerated set (breadth first search)
        sage: RecursivelyEnumeratedSet([0], f, enumeration='naive')
        A recursively enumerated set (naive search)
        sage: RecursivelyEnumeratedSet([0], f, enumeration='depth')
        A recursively enumerated set (depth first search)
    """
    def __init__(self, seeds, successors,
                 enumeration='depth', max_depth=float("inf"),
                 post_process=None, facade=None, category=None):
        r"""
        TESTS::

            sage: f = lambda a: [a+3, a+5]
            sage: C = RecursivelyEnumeratedSet([0], f)
            sage: C
            A recursively enumerated set (breadth first search)
        """
        assert enumeration in ['naive', 'depth', 'breadth'], \
                    "unknown enumeration(={})".format(enumeration)

        self._seeds = seeds
        self.successors = successors
        self._enumeration = enumeration
        self._max_depth = max_depth

        if post_process is not None:
            self.post_process = post_process
        self._graded_component = None
        Parent.__init__(self, facade=facade, category=EnumeratedSets().or_subcategory(category))

    def __reduce__(self):
        r"""
        Return a tuple of three elements:

        - The function :func:`RecursivelyEnumeratedSet`
        - Arguments for the function :func:`RecursivelyEnumeratedSet`
        - The actual state of ``self``.

        EXAMPLES::

            sage: C = RecursivelyEnumeratedSet((1, 2, 3), factor)
            sage: loads(dumps(C))
            A recursively enumerated set (breadth first search)
        """
        try:
            pp = self.post_process
        except AttributeError:
            pp = None

        classname = self.__class__.__name__
        if classname.startswith('RecursivelyEnumeratedSet_graded'):
            struct = 'graded'
        elif classname.startswith('RecursivelyEnumeratedSet_symmetric'):
            struct = 'symmetric'
        elif classname.startswith('RecursivelyEnumeratedSet_forest'):
            struct = 'forest'
        elif classname.startswith('RecursivelyEnumeratedSet_generic'):
            struct = None

        args = (self._seeds, self.successors, struct,
                self._enumeration, self._max_depth, pp)
        return (RecursivelyEnumeratedSet, args, self.__getstate__())

    def __getstate__(self):
        r"""
        Get the current state of ``self``. Used in pickling.

        EXAMPLES::

            sage: C = RecursivelyEnumeratedSet((1, 2, 3), factor)
            sage: C.__getstate__()
            (None,)
        """
        return (self._graded_component, )

    def __setstate__(self, l):
        r"""
        Set the state of ``self``. Used in pickling.

        INPUT:

        - ``l`` -- the state in the pickle

        EXAMPLES::

            sage: C = RecursivelyEnumeratedSet((1, 2, 3), factor)
            sage: C.__setstate__(C.__getstate__())
        """
        self._graded_component = l[0]
        # Since github issue #21312, the graded component iterator is not used
        # anymore but maybe some previously pickled object still have it
        # self._graded_component_it = l[1]

    def __len__(self):
        """
        Disable ``__len__()`` from :class:`Parent` :issue:`12955`.

        Because Python assumes ``__len__()`` is fast and we cannot
        have a fast default implementation.

        EXAMPLES::

            sage: f = lambda a: [a+3, a+5]
            sage: C = RecursivelyEnumeratedSet([0], f)
            sage: len(C)
            Traceback (most recent call last):
            ...
            TypeError: cannot compute length of A recursively enumerated set (breadth first search)
        """
        raise TypeError(f"cannot compute length of {self}")

    def __iter__(self):
        r"""
        Iterate on the elements of ``self``.

        The enumeration is done depth first or breadth first depending on
        the value of ``self._enumeration``.

        EXAMPLES::

            sage: f = lambda a: [a+3, a+5]
            sage: it_naive = iter(RecursivelyEnumeratedSet([0], f, enumeration='naive'))
            sage: it_depth = iter(RecursivelyEnumeratedSet([0], f, enumeration='depth'))
            sage: it_breadth = iter(RecursivelyEnumeratedSet([0], f, enumeration='breadth'))
            sage: sorted([next(it_naive) for _ in range(10)])
            [0, 3, 5, 6, 8, 9, 10, 11, 12, 13]
            sage: [next(it_depth) for _ in range(10)]
            [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
            sage: [next(it_breadth) for _ in range(10)]
            [0, 3, 5, 6, 8, 10, 9, 11, 13, 15]
        """
        if self._enumeration == 'naive':
            return self.naive_search_iterator()
        elif self._enumeration == 'breadth':
            return self.breadth_first_search_iterator(max_depth=self._max_depth)
        elif self._enumeration == 'depth':
            return self.depth_first_search_iterator()

        raise ValueError("unknown value for enumeration(={})".format(self._enumeration))

    def __contains__(self, elt):
        r"""
        Return ``True`` if ``elt`` is in ``self``.

        .. WARNING::

           This is achieved by iterating through the elements using the
           default enumeration until ``elt`` is found. In particular, this
           method will never stop when ``elt`` is not in ``self`` and
           ``self`` is infinite or when ``elt`` is in ``self`` but the
           enumeration is not appropriate.

        EXAMPLES::

            sage: f = lambda a:[a+3,a+5]
            sage: R = RecursivelyEnumeratedSet([0], f)
            sage: R
            A recursively enumerated set (breadth first search)
            sage: 8 in R
            True

        ::

            sage: R = RecursivelyEnumeratedSet([0], f, enumeration='depth')
            sage: R
            A recursively enumerated set (depth first search)
            sage: it = iter(R)
            sage: [next(it) for _ in range(6)]
            [0, 5, 10, 15, 20, 25]
            sage: 8 in R     # (should return True) not tested: does not terminate
            sage: 7 in R     # (should return False) not tested: does not terminate
        """
        return any(node == elt for node in self)

    def _repr_(self):
        r"""
        TESTS::

            sage: f = lambda x: [x-1, x+1]
            sage: RecursivelyEnumeratedSet([1], f, structure=None)
            A recursively enumerated set (breadth first search)

        ::

            sage: RecursivelyEnumeratedSet([1], f, structure='graded')
            A recursively enumerated set with a graded structure (breadth first search)

        ::

            sage: RecursivelyEnumeratedSet([1], f, structure='symmetric')
            A recursively enumerated set with a symmetric structure (breadth first search)

        When ``max_depth`` is set::

            sage: RecursivelyEnumeratedSet([1], f, structure='symmetric', max_depth=4)
            A recursively enumerated set with a symmetric structure (breadth
            first search) with max_depth=4
        """
        L = ["A recursively enumerated set"]
        classname = self.__class__.__name__
        if classname.startswith('RecursivelyEnumeratedSet_graded'):
            L.append("with a graded structure")
        elif classname.startswith('RecursivelyEnumeratedSet_symmetric'):
            L.append("with a symmetric structure")
        elif classname.startswith('RecursivelyEnumeratedSet_forest'):
            L.append("with a forest structure")

        if self._enumeration in ['depth', 'breadth']:
            L.append("({} first search)".format(self._enumeration))
        else:
            L.append("({} search)".format(self._enumeration))

        if not self._max_depth == float('inf'):
            L.append("with max_depth={}".format(self._max_depth))
        return " ".join(L)

    cpdef seeds(self):
        r"""
        Return an iterable over the seeds of ``self``.

        EXAMPLES::

            sage: R = RecursivelyEnumeratedSet([1], lambda x: [x + 1, x - 1])
            sage: R.seeds()
            [1]
        """
        return self._seeds

    # using this in a .pyx file makes sage crash at startup
    # @abstract_method
    # def successors(self, x):
    #     r"""
    #     Return the successors of the element ``x``
    #
    #     OUTPUT:
    #
    #         an iterable
    #
    #     EXAMPLES::
    #
    #         sage: R = RecursivelyEnumeratedSet([1], lambda x: [x+1, x-1])
    #         sage: R.successors(4)
    #         [5, 3]
    #     """

    def graded_component_iterator(self):
        r"""
        Iterate over the graded components of ``self``.

        A graded component is a set of elements of the same depth.

        It is currently implemented only for graded or symmetric structure.

        OUTPUT: an iterator of sets

        EXAMPLES::

            sage: f = lambda a: [a+3, a+5]
            sage: C = RecursivelyEnumeratedSet([0], f)
            sage: it = C.graded_component_iterator()    # todo: not implemented
        """
        raise NotImplementedError("graded_component_iterator method currently"
                                  " implemented only for graded or symmetric structure")

    cpdef graded_component(self, depth):
        r"""
        Return the graded component of given depth.

        This method caches each lower graded component.

        A graded component is a set of elements of the same depth where the
        depth of an element is its minimal distance to a root.

        It is currently implemented only for graded or symmetric structure.

        INPUT:

        - ``depth`` -- integer

        OUTPUT: set

        EXAMPLES::

            sage: f = lambda a: [a+3, a+5]
            sage: C = RecursivelyEnumeratedSet([0], f)
            sage: C.graded_component(0)
            Traceback (most recent call last):
            ...
            NotImplementedError: graded_component_iterator method currently implemented only for graded or symmetric structure
        """
        raise NotImplementedError("graded_component_iterator method currently"
                                  " implemented only for graded or symmetric structure")

    def elements_of_depth_iterator(self, depth):
        r"""
        Iterate over the elements of ``self`` of given depth.

        An element of depth `n` can be obtained by applying the
        successor function `n` times to a seed.

        INPUT:

        - ``depth`` -- integer

        OUTPUT: an iterator

        EXAMPLES::

            sage: f = lambda a: [a-1, a+1]
            sage: S = RecursivelyEnumeratedSet([5, 10], f, structure='symmetric')
            sage: it = S.elements_of_depth_iterator(2)
            sage: sorted(it)
            [3, 7, 8, 12]
        """
        return iter(self.graded_component(depth))

    def breadth_first_search_iterator(self, max_depth=None):
        r"""
        Iterate on the elements of ``self`` (breadth first).

        This code remembers every element generated.

        The elements are guaranteed to be enumerated in the order in which they
        are first visited (left-to-right traversal).

        INPUT:

        - ``max_depth`` -- (default: ``self._max_depth``) specifies the
          maximal depth to which elements are computed

        EXAMPLES::

            sage: f = lambda a: [a+3, a+5]
            sage: C = RecursivelyEnumeratedSet([0], f)
            sage: it = C.breadth_first_search_iterator()
            sage: [next(it) for _ in range(10)]
            [0, 3, 5, 6, 8, 10, 9, 11, 13, 15]
        """
        if max_depth is None:
            max_depth = self._max_depth
        current_level = self._seeds
        known = set(current_level)
        if max_depth >= 0:
            yield from current_level
        depth = 0
        while current_level and depth < max_depth:
            next_level = []
            for x in current_level:
                for y in self.successors(x):
                    if y is None or y in known:
                        continue
                    yield y
                    next_level.append(y)
                    known.add(y)
            current_level = next_level
            depth += 1

    def _breadth_first_search_iterator_using_queue(self):
        r"""
        Iterate on the elements of ``self`` (breadth first).

        This code remembers every element generated and uses python
        queues. It is 3 times slower than the other one.

        See :wikipedia:`Breadth-first_search`.

        EXAMPLES::

            sage: f = lambda a: [a+3, a+5]
            sage: C = RecursivelyEnumeratedSet([0], f)
            sage: it = C._breadth_first_search_iterator_using_queue()
            sage: [next(it) for _ in range(10)]
            [0, 3, 5, 6, 8, 10, 9, 11, 13, 15]
        """
        cdef set known
        known = set(self._seeds)
        q = deque(self._seeds)
        while q:
            x = q.popleft()
            yield x
            for y in self.successors(x):
                if y is None or y in known:
                    continue
                q.append(y)
                known.add(y)

    def naive_search_iterator(self):
        r"""
        Iterate on the elements of ``self`` (in no particular order).

        This code remembers every element generated.

        TESTS:

        We compute all the permutations of 3::

            sage: # needs sage.combinat
            sage: seeds = [Permutation([1,2,3])]
            sage: succ = attrcall("permutohedron_succ")
            sage: R = RecursivelyEnumeratedSet(seeds, succ)
            sage: sorted(R.naive_search_iterator())
            [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
        """
        cdef set known, todo
        known = set(self._seeds)
        todo = known.copy()
        while todo:
            x = todo.pop()
            yield x
            for y in self.successors(x):
                if y is None or y in known:
                    continue
                todo.add(y)
                known.add(y)

    def depth_first_search_iterator(self):
        r"""
        Iterate on the elements of ``self`` (depth first).

        This code remembers every element generated.

        The elements are traversed right-to-left, so the last element returned
        by the successor function is visited first.

        See :wikipedia:`Depth-first_search`.

        EXAMPLES::

            sage: f = lambda a: [a+3, a+5]
            sage: C = RecursivelyEnumeratedSet([0], f)
            sage: it = C.depth_first_search_iterator()
            sage: [next(it) for _ in range(10)]
            [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
        """
        cdef list stack
        cdef set known
        stack = list(self._seeds)
        known = set()
        while stack:
            x = stack.pop()
            if x is None or x in known:
                continue
            yield x
            known.add(x)
            for y in self.successors(x):
                stack.append(y)

    def to_digraph(self, max_depth=None, loops=True, multiedges=True):
        r"""
        Return the directed graph of the recursively enumerated set.

        INPUT:

        - ``max_depth`` -- (default: ``self._max_depth``) specifies the
          maximal depth for which outgoing edges of elements are computed
        - ``loops`` -- boolean (default: ``True``); option for the digraph
        - ``multiedges`` -- boolean (default: ``True``); option of the digraph

        OUTPUT: a directed graph

        .. WARNING::

            If the set is infinite, this will loop forever unless ``max_depth``
            is finite.

        EXAMPLES::

            sage: child = lambda i: [(i+3) % 10, (i+8) % 10]
            sage: R = RecursivelyEnumeratedSet([0], child)
            sage: R.to_digraph()                                                        # needs sage.graphs
            Looped multi-digraph on 10 vertices

        Digraph of a recursively enumerated set with a symmetric structure of
        infinite cardinality using ``max_depth`` argument::

            sage: succ = lambda a: [(a[0]-1,a[1]), (a[0],a[1]-1), (a[0]+1,a[1]), (a[0],a[1]+1)]
            sage: seeds = [(0,0)]
            sage: C = RecursivelyEnumeratedSet(seeds, succ, structure='symmetric')
            sage: C.to_digraph(max_depth=3)                                             # needs sage.graphs
            Looped multi-digraph on 41 vertices

        The ``max_depth`` argument can be given at the creation of the set::

            sage: C = RecursivelyEnumeratedSet(seeds, succ, structure='symmetric',
            ....:                              max_depth=2)
            sage: C.to_digraph()                                                        # needs sage.graphs
            Looped multi-digraph on 25 vertices

        Digraph of a recursively enumerated set with a graded structure::

            sage: f = lambda a: [a+1, a+I]
            sage: C = RecursivelyEnumeratedSet([0], f, structure='graded')
            sage: C.to_digraph(max_depth=4)                                             # needs sage.graphs
            Looped multi-digraph on 21 vertices
        """
        successors = self.successors
        it = self.breadth_first_search_iterator(max_depth=max_depth)
        E = [(u, v) for u in it for v in successors(u)]
        from sage.graphs.digraph import DiGraph
        return DiGraph(E, format='list_of_edges', loops=loops,
                       multiedges=multiedges)


cdef class RecursivelyEnumeratedSet_symmetric(RecursivelyEnumeratedSet_generic):
    r"""
    Generic tool for constructing ideals of a symmetric relation.

    INPUT:

    - ``seeds`` -- list (or iterable) of hashable objects
    - ``successors`` -- function (or callable) returning a list (or iterable)
    - ``enumeration`` -- ``'depth'``, ``'breadth'`` or ``None`` (default: ``None``)
    - ``max_depth`` -- integer (default: ``float("inf")``)

    EXAMPLES::

        sage: f = lambda a: [a-1,a+1]
        sage: C = RecursivelyEnumeratedSet([0], f, structure='symmetric')
        sage: C
        A recursively enumerated set with a symmetric structure (breadth first search)
        sage: it = iter(C)
        sage: [next(it) for _ in range(7)]
        [0, -1, 1, -2, 2, -3, 3]

    TESTS:

    Do not use lambda functions for saving purposes::

        sage: f = lambda a: [a-1,a+1]
        sage: C = RecursivelyEnumeratedSet([0], f, structure='symmetric')
        sage: loads(dumps(C))
        Traceback (most recent call last):
        ...
        PicklingError: ...

    This works in the command line but apparently not as a doctest::

        sage: def f(a): return [a-1,a+1]
        sage: C = RecursivelyEnumeratedSet([0], f, structure='symmetric')
        sage: loads(dumps(C))
        Traceback (most recent call last):
        ...
        PicklingError: ...
    """

    def breadth_first_search_iterator(self, max_depth=None):
        r"""
        Iterate on the elements of ``self`` (breadth first).

        This iterator makes use of the graded structure by remembering only
        the last two graded components since the structure is symmetric.

        The elements are guaranteed to be enumerated in the order in which they
        are first visited (left-to-right traversal).

        INPUT:

        - ``max_depth`` -- (default: ``self._max_depth``) specifies the
          maximal depth to which elements are computed

        EXAMPLES::

            sage: f = lambda a: [(a[0]-1,a[1]), (a[0],a[1]-1), (a[0]+1,a[1]), (a[0],a[1]+1)]
            sage: C = RecursivelyEnumeratedSet([(0,0)], f, structure='symmetric')
            sage: s = list(C.breadth_first_search_iterator(max_depth=2)); s
            [(0, 0),
             (-1, 0), (0, -1), (1, 0), (0, 1),
             (-2, 0), (-1, -1), (-1, 1), (0, -2), (1, -1), (2, 0), (1, 1), (0, 2)]

        This iterator is used by default for symmetric structure::

            sage: it = iter(C)
            sage: s == [next(it) for _ in range(13)]
            True

        TESTS:

        Check that :issue:`28674` is fixed::

            sage: D = RecursivelyEnumeratedSet([(0,0)], f)
            sage: s == list(D.breadth_first_search_iterator(max_depth=2))
            True
        """
        cdef list C
        cdef set set_A, set_B
        cdef int depth
        if max_depth is None:
            max_depth = self._max_depth

        set_A = set()
        B = self._seeds
        set_B = set(B)
        if max_depth >= 0:
            yield from B
        depth = 0
        while B and depth < max_depth:
            C = list()
            set_C = set()
            for x in B:
                for y in self.successors(x):
                    if y is None or y in set_C or y in set_A or y in set_B:
                        continue
                    yield y
                    C.append(y)
                    set_C.add(y)
            set_A = set_B
            set_B = set_C
            B = C
            depth += 1

    def graded_component_iterator(self):
        r"""
        Iterate over the graded components of ``self``.

        A graded component is a set of elements of the same depth.

        The enumeration remembers only the last two graded components
        generated since the structure is symmetric.

        OUTPUT: an iterator of sets

        EXAMPLES::

            sage: f = lambda a: [a-1, a+1]
            sage: S = RecursivelyEnumeratedSet([10], f, structure='symmetric')
            sage: it = S.graded_component_iterator()
            sage: [sorted(next(it)) for _ in range(5)]
            [[10], [9, 11], [8, 12], [7, 13], [6, 14]]

        Starting with two generators::

            sage: f = lambda a: [a-1, a+1]
            sage: S = RecursivelyEnumeratedSet([5, 10], f, structure='symmetric')
            sage: it = S.graded_component_iterator()
            sage: [sorted(next(it)) for _ in range(5)]
            [[5, 10], [4, 6, 9, 11], [3, 7, 8, 12], [2, 13], [1, 14]]

        Gaussian integers::

            sage: # needs sage.symbolic
            sage: def f(a):
            ....:     return [a + 1, a + I]
            sage: S = RecursivelyEnumeratedSet([0], f, structure='symmetric')
            sage: it = S.graded_component_iterator()
            sage: [sorted(next(it)) for _ in range(7)]
            [[0],
             [I, 1],
             [2*I, I + 1, 2],
             [3*I, 2*I + 1, I + 2, 3],
             [4*I, 3*I + 1, 2*I + 2, I + 3, 4],
             [5*I, 4*I + 1, 3*I + 2, 2*I + 3, I + 4, 5],
             [6*I, 5*I + 1, 4*I + 2, 3*I + 3, 2*I + 4, I + 5, 6]]

        TESTS:

        Note that interrupting the computation (``KeyboardInterrupt`` for
        instance) breaks the iterator::

            sage: # needs sage.symbolic
            sage: def f(a):
            ....:     sleep(0.05r)
            ....:     return [a - 1, a + 1]
            sage: C = RecursivelyEnumeratedSet([0], f, structure='symmetric')
            sage: it = C.graded_component_iterator()
            sage: next(it)
            {0}
            sage: next(it)
            {-1, 1}
            sage: from cysignals.alarm import alarm
            sage: alarm(0.02); next(it)
            Traceback (most recent call last):
            ...
            AlarmInterrupt
            sage: next(it)
            Traceback (most recent call last):
            ...
            StopIteration
        """
        cdef set A, B
        A = set()
        B = set(self._seeds)
        while B:
            yield B
            A, B = B, self._get_next_graded_component(A, B)

    cpdef graded_component(self, depth):
        r"""
        Return the graded component of given depth.

        This method caches each lower graded component. See
        :meth:`graded_component_iterator` to generate each graded component
        without caching the previous ones.

        A graded component is a set of elements of the same depth where the
        depth of an element is its minimal distance to a root.

        INPUT:

        - ``depth`` -- integer

        OUTPUT: set

        EXAMPLES::

            sage: f = lambda a: [a-1,a+1]
            sage: C = RecursivelyEnumeratedSet([10, 15], f, structure='symmetric')
            sage: for i in range(5): sorted(C.graded_component(i))
            [10, 15]
            [9, 11, 14, 16]
            [8, 12, 13, 17]
            [7, 18]
            [6, 19]

        TESTS:

        We make sure that :issue:`21312` is fixed::

            sage: def f(a):
            ....:    sleep(0.1r)
            ....:    return [a - 1, a + 1]
            sage: C = RecursivelyEnumeratedSet([0], f, structure='symmetric')
            sage: from cysignals.alarm import alarm
            sage: alarm(0.45); C.graded_component(10)
            Traceback (most recent call last):
            ...
            AlarmInterrupt
            sage: C.graded_component(1)
            {-1, 1}
            sage: C.graded_component(2)
            {-2, 2}
            sage: C.graded_component(3)
            {-3, 3}
            sage: C.graded_component(4)
            {-4, 4}
            sage: C.graded_component(5)
            {-5, 5}
        """
        cdef set A, B, C
        if self._graded_component is None:
            A = set()
            B = set(self._seeds)
            C = self._get_next_graded_component(A, B)
            self._graded_component = [B, C]
        while len(self._graded_component) <= depth:
            A = self._graded_component[-2]
            B = self._graded_component[-1]
            C = self._get_next_graded_component(A, B)
            self._graded_component.append(C)
        return self._graded_component[depth]

    cdef set _get_next_graded_component(self, set A, set B):
        r"""
        Return the set of elements of depth `n+1`.

        INPUT:

        - ``A`` -- set, the set of elements of depth n-1
        - ``B`` -- set, the set of elements of depth n

        OUTPUT: ``C``; the set of elements of depth n+1

        .. TODO::

            Can :class:`collections.OrderedDict` help maintain the
            breadth first search enumeration for each graded component?

        EXAMPLES::

            sage: f = lambda a: [a-1, a+1]
            sage: S = RecursivelyEnumeratedSet([5, 10], f, structure='symmetric')
            sage: it = S.graded_component_iterator()
            sage: [sorted(next(it)) for _ in range(3)] # indirect doctest
            [[5, 10], [4, 6, 9, 11], [3, 7, 8, 12]]
        """
        cdef set C
        C = set()
        for x in B:
            for y in self.successors(x):
                if (y is None or y in A or y in B):
                    continue
                C.add(y)
        return C


cdef class RecursivelyEnumeratedSet_graded(RecursivelyEnumeratedSet_generic):
    r"""
    Generic tool for constructing ideals of a graded relation.

    INPUT:

    - ``seeds`` -- list (or iterable) of hashable objects
    - ``successors`` -- function (or callable) returning a list (or iterable)
    - ``enumeration`` -- ``'depth'``, ``'breadth'`` or ``None`` (default: ``None``)
    - ``max_depth`` -- integer (default: ``float("inf")``)

    EXAMPLES::

        sage: f = lambda a: [(a[0]+1,a[1]), (a[0],a[1]+1)]
        sage: C = RecursivelyEnumeratedSet([(0,0)], f, structure='graded', max_depth=3)
        sage: C
        A recursively enumerated set with a graded structure (breadth first
        search) with max_depth=3
        sage: list(C)
        [(0, 0),
         (1, 0), (0, 1),
         (2, 0), (1, 1), (0, 2),
         (3, 0), (2, 1), (1, 2), (0, 3)]
    """
    def breadth_first_search_iterator(self, max_depth=None):
        r"""
        Iterate on the elements of ``self`` (breadth first).

        This iterator makes use of the graded structure by remembering only
        the elements of the current depth.

        The elements are guaranteed to be enumerated in the order in which they
        are first visited (left-to-right traversal).

        INPUT:

        - ``max_depth`` -- (default: ``self._max_depth``) specifies the
          maximal depth to which elements are computed

        EXAMPLES::

            sage: f = lambda a: [(a[0]+1,a[1]), (a[0],a[1]+1)]
            sage: C = RecursivelyEnumeratedSet([(0,0)], f, structure='graded')
            sage: list(C.breadth_first_search_iterator(max_depth=3))
            [(0, 0),
             (1, 0), (0, 1),
             (2, 0), (1, 1), (0, 2),
             (3, 0), (2, 1), (1, 2), (0, 3)]
        """
        cdef list next_level
        cdef set set_next_level
        cdef int depth
        if max_depth is None:
            max_depth = self._max_depth
        current_level = self._seeds
        if max_depth >= 0:
            yield from current_level
        depth = 0
        while current_level and depth < max_depth:
            next_level = list()
            set_next_level = set()

            for x in current_level:
                for y in self.successors(x):
                    if y is None or y in set_next_level:
                        continue
                    yield y
                    next_level.append(y)
                    set_next_level.add(y)
            current_level = next_level
            depth += 1

    def graded_component_iterator(self):
        r"""
        Iterate over the graded components of ``self``.

        A graded component is a set of elements of the same depth.

        The algorithm remembers only the current graded component generated
        since the structure is graded.

        OUTPUT: an iterator of sets

        EXAMPLES::

            sage: f = lambda a: [(a[0]+1,a[1]), (a[0],a[1]+1)]
            sage: C = RecursivelyEnumeratedSet([(0,0)], f, structure='graded', max_depth=3)
            sage: it = C.graded_component_iterator()
            sage: for _ in range(4): sorted(next(it))
            [(0, 0)]
            [(0, 1), (1, 0)]
            [(0, 2), (1, 1), (2, 0)]
            [(0, 3), (1, 2), (2, 1), (3, 0)]

        TESTS:

        Make sure that :issue:`20225` is fixed::

            sage: child = lambda k:[2*k,2*k+1] if k<8 else []
            sage: root = [0]
            sage: R = RecursivelyEnumeratedSet(root, child, structure='graded')
            sage: it = R.graded_component_iterator()
            sage: for _ in range(7): next(it)
            {0}
            {1}
            {2, 3}
            {4, 5, 6, 7}
            {8, 9, 10, 11, 12, 13, 14, 15}
            set()
            set()
        """
        cdef set B
        B = set(self._seeds)
        while True:
            yield B
            B = self._get_next_graded_component(B)

    cpdef graded_component(self, depth):
        r"""
        Return the graded component of given depth.

        This method caches each lower graded component. See
        :meth:`graded_component_iterator` to generate each graded component
        without caching the previous ones.

        A graded component is a set of elements of the same depth where the
        depth of an element is its minimal distance to a root.

        INPUT:

        - ``depth`` -- integer

        OUTPUT: set

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: def f(a):
            ....:     return [a + 1, a + I]
            sage: C = RecursivelyEnumeratedSet([0], f, structure='graded')
            sage: for i in range(5): sorted(C.graded_component(i))
            [0]
            [I, 1]
            [2*I, I + 1, 2]
            [3*I, 2*I + 1, I + 2, 3]
            [4*I, 3*I + 1, 2*I + 2, I + 3, 4]

        TESTS:

        We make sure that :issue:`21312` is fixed::

            sage: # needs sage.symbolic
            sage: def f(a):
            ....:    sleep(0.1r)
            ....:    return [a + 1, a + I]
            sage: C = RecursivelyEnumeratedSet([0], f, structure='graded')
            sage: from cysignals.alarm import alarm
            sage: alarm(0.45); C.graded_component(10)
            Traceback (most recent call last):
            ...
            AlarmInterrupt
            sage: C.graded_component(2)
            {2*I, I + 1, 2}
            sage: C.graded_component(3)
            {3*I, 2*I + 1, I + 2, 3}
            sage: C.graded_component(4)
            {4*I, 3*I + 1, 2*I + 2, I + 3, 4}
        """
        cdef set B, C
        if self._graded_component is None:
            B = set(self._seeds)
            self._graded_component = [B]
        while len(self._graded_component) <= depth:
            B = self._graded_component[-1]
            C = self._get_next_graded_component(B)
            self._graded_component.append(C)
        return self._graded_component[depth]

    cdef set _get_next_graded_component(self, set B):
        r"""
        Return the set of elements of depth `n+1`.

        INPUT:

        - ``B`` -- set, the set of elements of depth `n`

        OUTPUT: ``C``; the set of elements of depth `n+1`

        .. TODO::

            Can :class:`collections.OrderedDict` help maintain the
            breadth first search enumeration for each graded component?

        EXAMPLES::

            sage: f = lambda a: [(a[0]+1,a[1]), (a[0],a[1]+1)]
            sage: C = RecursivelyEnumeratedSet([(0,0)], f, structure='graded')
            sage: it = C.graded_component_iterator()
            sage: [sorted(next(it)) for _ in range(2)] # indirect doctest
            [[(0, 0)], [(0, 1), (1, 0)]]
        """
        cdef set C
        C = set()
        for x in B:
            for y in self.successors(x):
                if (y is None or y in B):
                    continue
                C.add(y)
        return C


def _imap_and_filter_none(function, iterable):
    r"""
    Return an iterator over the elements ``function(x)``, where ``x``
    iterates through ``iterable``, such that ``function(x)`` is not
    ``None``.

    EXAMPLES::

        sage: from sage.sets.recursively_enumerated_set import _imap_and_filter_none
        sage: p = _imap_and_filter_none(lambda x: x if is_prime(x) else None, range(15))
        sage: [next(p), next(p), next(p), next(p), next(p), next(p)]
        [2, 3, 5, 7, 11, 13]
        sage: p = _imap_and_filter_none(lambda x: x + x, ['a','b','c','d','e'])
        sage: [next(p), next(p), next(p), next(p), next(p)]
        ['aa', 'bb', 'cc', 'dd', 'ee']
    """
    for x in iterable:
        x = function(x)
        if x is not None:
            yield x


def search_forest_iterator(roots, children, algorithm='depth'):
    r"""
    Return an iterator on the nodes of the forest having the given
    roots, and where ``children(x)`` returns the children of the node ``x``
    of the forest.  Note that every node of the tree is returned,
    not simply the leaves.

    INPUT:

    - ``roots`` -- list (or iterable)
    - ``children`` -- a function returning a list (or iterable)
    - ``algorithm`` -- ``'depth'`` or ``'breadth'`` (default: ``'depth'``)

    EXAMPLES:

    We construct the prefix tree of binary sequences of length at most
    three, and enumerate its nodes::

        sage: from sage.sets.recursively_enumerated_set import search_forest_iterator
        sage: list(search_forest_iterator([[]], lambda l: [l + [0], l + [1]]
        ....:                                   if len(l) < 3 else []))
        [[], [0], [0, 0], [0, 0, 0], [0, 0, 1], [0, 1], [0, 1, 0],
         [0, 1, 1], [1], [1, 0], [1, 0, 0], [1, 0, 1], [1, 1], [1, 1, 0], [1, 1, 1]]

    By default, the nodes are iterated through by depth first search.
    We can instead use a breadth first search (increasing depth)::

        sage: list(search_forest_iterator([[]], lambda l: [l + [0], l + [1]]
        ....:                                   if len(l) < 3 else [],
        ....:                             algorithm='breadth'))
        [[],
         [0], [1],
         [0, 0], [0, 1], [1, 0], [1, 1],
         [0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1],
         [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]]

    This allows for iterating through trees of infinite depth::

        sage: it = search_forest_iterator([[]], lambda l: [l + [0], l + [1]],
        ....:                             algorithm='breadth')
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
    # breadth first search. Setting position to -1 results in a depth search
    # (you ask the children for the last node you met). Setting
    # position on 0 results in a breadth search (enumerate all the
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


class RecursivelyEnumeratedSet_forest(Parent):
    r"""
    The enumerated set of the nodes of the forest having the given
    ``roots``, and where ``children(x)`` returns the children of the
    node ``x`` of the forest.

    See also :class:`sage.combinat.backtrack.GenericBacktracker`,
    :class:`RecursivelyEnumeratedSet_graded`, and
    :class:`RecursivelyEnumeratedSet_symmetric`.

    INPUT:

    - ``roots`` -- list (or iterable)
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

        sage: from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest
        sage: S = RecursivelyEnumeratedSet_forest( [[]],
        ....:     lambda l: [l + [0], l + [1]] if len(l) < 3 else [],
        ....:     category=FiniteEnumeratedSets())
        sage: S.list()
        [[],
         [0], [0, 0], [0, 0, 0], [0, 0, 1], [0, 1], [0, 1, 0], [0, 1, 1],
         [1], [1, 0], [1, 0, 0], [1, 0, 1], [1, 1], [1, 1, 0], [1, 1, 1]]

    ``RecursivelyEnumeratedSet_forest`` needs to be explicitly told that the set is
    finite for the following to work::

        sage: S.category()
        Category of finite enumerated sets
        sage: S.cardinality()
        15

    We proceed with the set of all lists of letters in ``0,1,2``
    without repetitions, ordered by increasing length (i.e. using a
    breadth first search through the tree)::

        sage: from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest
        sage: tb = RecursivelyEnumeratedSet_forest( [[]],
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

        sage: from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest
        sage: I = RecursivelyEnumeratedSet_forest([(0,0)],
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

    The class of a parent `A` may derive from :class:`RecursivelyEnumeratedSet_forest` so
    that `A` can benefit from enumeration tools. As a running example,
    we consider the problem of enumerating integers whose binary
    expansion have at most three nonzero digits. For example, `3 =
    2^1 + 2^0` has two nonzero digits. `15 = 2^3 + 2^2 + 2^1 + 2^0`
    has four nonzero digits. In fact, `15` is the smallest integer
    which is not in the enumerated set.

    To achieve this, we use ``RecursivelyEnumeratedSet_forest`` to enumerate binary tuples
    with at most three nonzero digits, apply a post processing to
    recover the corresponding integers, and discard tuples finishing
    by zero.

    A first approach is to pass the ``roots`` and ``children``
    functions as arguments to :meth:`RecursivelyEnumeratedSet_forest.__init__`::

        sage: from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest
        sage: class A(UniqueRepresentation, RecursivelyEnumeratedSet_forest):
        ....:     def __init__(self):
        ....:         RecursivelyEnumeratedSet_forest.__init__(self, [()],
        ....:             lambda x: [x + (0,), x + (1,)] if sum(x) < 3 else [],
        ....:             lambda x: sum(x[i]*2^i for i in range(len(x)))
        ....:                           if sum(x) != 0 and x[-1] != 0 else None,
        ....:             algorithm='breadth',
        ....:             category=InfiniteEnumeratedSets())
        sage: MyForest = A(); MyForest
        An enumerated set with a forest structure
        sage: MyForest.category()
        Category of infinite enumerated sets
        sage: p = iter(MyForest)
        sage: [next(p) for i in range(30)]
        [1, 2, 3, 4, 6, 5, 7, 8, 12, 10, 14, 9, 13, 11, 16, 24,
         20, 28, 18, 26, 22, 17, 25, 21, 19, 32, 48, 40, 56, 36]

    An alternative approach is to implement ``roots`` and ``children``
    as methods of the subclass (in fact they could also be attributes
    of `A`). Namely, ``A.roots()`` must return an iterable containing
    the enumeration generators, and ``A.children(x)`` must return an
    iterable over the children of `x`. Optionally, `A` can have a
    method or attribute such that ``A.post_process(x)`` returns the
    desired output for the node ``x`` of the tree::

        sage: from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest
        sage: class A(UniqueRepresentation, RecursivelyEnumeratedSet_forest):
        ....:     def __init__(self):
        ....:         RecursivelyEnumeratedSet_forest.__init__(self, algorithm='breadth',
        ....:                               category=InfiniteEnumeratedSets())
        ....:     def roots(self):
        ....:         return [()]
        ....:     def children(self, x):
        ....:         if sum(x) < 3:
        ....:             return [x + (0,), x + (1,)]
        ....:         else:
        ....:             return []
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
        [1, 2, 3, 4, 6, 5, 7, 8, 12, 10, 14, 9, 13, 11, 16, 24,
         20, 28, 18, 26, 22, 17, 25, 21, 19, 32, 48, 40, 56, 36]

    .. warning::

        A :class:`RecursivelyEnumeratedSet_forest` instance is picklable if and only if
        the input functions are themselves picklable. This excludes
        anonymous or interactively defined functions::

            sage: def children(x):
            ....:     return [x + 1]
            sage: S = RecursivelyEnumeratedSet_forest([1], children, category=InfiniteEnumeratedSets())
            sage: dumps(S)
            Traceback (most recent call last):
            ...
            PicklingError: Can't pickle <...function...>: attribute lookup ... failed

        Let us now fake ``children`` being defined in a Python module::

            sage: import __main__
            sage: __main__.children = children
            sage: S = RecursivelyEnumeratedSet_forest([1], children, category=InfiniteEnumeratedSets())
            sage: loads(dumps(S))
            An enumerated set with a forest structure
    """
    def __init__(self, roots = None, children = None, post_process = None,
                 algorithm = 'depth', facade = None, category=None):
        r"""
        TESTS::

            sage: from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest
            sage: S = RecursivelyEnumeratedSet_forest(NN, lambda x: [], lambda x: x^2 if x.is_prime() else None)
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

            sage: from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest
            sage: RecursivelyEnumeratedSet_forest( [1], lambda x: [x+1])
            An enumerated set with a forest structure
        """
        return "An enumerated set with a forest structure"

    def roots(self):
        r"""
        Return an iterable over the roots of ``self``.

        EXAMPLES::

            sage: from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest
            sage: I = RecursivelyEnumeratedSet_forest([(0,0)], lambda l: [(l[0]+1, l[1]), (l[0], 1)] if l[1] == 0 else [(l[0], l[1]+1)])
            sage: [i for i in I.roots()]
            [(0, 0)]
            sage: I = RecursivelyEnumeratedSet_forest([(0,0),(1,1)], lambda l: [(l[0]+1, l[1]), (l[0], 1)] if l[1] == 0 else [(l[0], l[1]+1)])
            sage: [i for i in I.roots()]
            [(0, 0), (1, 1)]
        """
        return self._roots

    @abstract_method
    def children(self, x):
        r"""
        Return the children of the element ``x``.

        The result can be a list, an iterable, an iterator, or even a
        generator.

        EXAMPLES::

            sage: from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest
            sage: I = RecursivelyEnumeratedSet_forest([(0,0)], lambda l: [(l[0]+1, l[1]), (l[0], 1)] if l[1] == 0 else [(l[0], l[1]+1)])
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

            sage: from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest
            sage: def children(l):
            ....:      return [l + [0], l + [1]]
            sage: C = RecursivelyEnumeratedSet_forest(([],), children)
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
        Return a depth first search iterator over the elements of ``self``.

        EXAMPLES::

            sage: from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest
            sage: f = RecursivelyEnumeratedSet_forest([[]],
            ....:                  lambda l: [l + [0], l + [1]] if len(l) < 3 else [])
            sage: list(f.depth_first_search_iterator())
            [[], [0], [0, 0], [0, 0, 0], [0, 0, 1], [0, 1], [0, 1, 0], [0, 1, 1],
                 [1], [1, 0], [1, 0, 0], [1, 0, 1], [1, 1], [1, 1, 0], [1, 1, 1]]
        """
        return iter(self)

    def breadth_first_search_iterator(self):
        r"""
        Return a breadth first search iterator over the elements of ``self``.

        EXAMPLES::

            sage: from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest
            sage: f = RecursivelyEnumeratedSet_forest([[]],
            ....:                  lambda l: [l+[0], l+[1]] if len(l) < 3 else [])
            sage: list(f.breadth_first_search_iterator())
            [[], [0], [1], [0, 0], [0, 1], [1, 0], [1, 1], [0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]]
            sage: S = RecursivelyEnumeratedSet_forest([(0,0)],
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
        An element of depth `n` can be obtained by applying the
        children function `n` times from a root. This function is not affected
        by post processing.

        EXAMPLES::

            sage: from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest
            sage: I = RecursivelyEnumeratedSet_forest([(0,0)], lambda l: [(l[0]+1, l[1]), (l[0], 1)] if l[1] == 0 else [(l[0], l[1]+1)])
            sage: list(I._elements_of_depth_iterator_rec(8))
            [(8, 0), (7, 1), (6, 2), (5, 3), (4, 4), (3, 5), (2, 6), (1, 7), (0, 8)]
            sage: I = RecursivelyEnumeratedSet_forest([[]], lambda l: [l+[0], l+[1]] if len(l) < 3 else [])
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
            yield from self.roots()
        else:
            for father in self._elements_of_depth_iterator_rec(depth - 1):
                yield from self.children(father)

    def elements_of_depth_iterator(self, depth=0):
        r"""
        Return an iterator over the elements of ``self`` of given depth.
        An element of depth `n` can be obtained by applying the
        children function `n` times from a root.

        EXAMPLES::

            sage: from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest
            sage: S = RecursivelyEnumeratedSet_forest([(0,0)] ,
            ....:        lambda x : [(x[0], x[1]+1)] if x[1] != 0 else [(x[0]+1,0), (x[0],1)],
            ....:        post_process = lambda x: x if ((is_prime(x[0]) and is_prime(x[1]))
            ....:                                        and ((x[0] - x[1]) == 2)) else None)
            sage: p = S.elements_of_depth_iterator(8)
            sage: next(p)
            (5, 3)
            sage: S = RecursivelyEnumeratedSet_forest(NN, lambda x : [],
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

            sage: from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest
            sage: S = RecursivelyEnumeratedSet_forest([[]], lambda l: [l + [0], l + [1]] if len(l) < 3 else [],
            ....:                                     category=FiniteEnumeratedSets())
            sage: [4] in S
            False
            sage: [1] in S
            True
            sage: [1,1,1,1] in S
            False
            sage: all(S.__contains__(i) for i in iter(S))
            True
            sage: S = RecursivelyEnumeratedSet_forest([1], lambda x: [x + 1], category=InfiniteEnumeratedSets())
            sage: 1 in S
            True
            sage: 732 in S
            True
            sage: -1 in S  # not tested : Will never stop

        The algorithm uses a random enumeration of the nodes of the
        forest. This choice was motivated by examples in which both
        depth first search and breadth first search failed. The
        following example enumerates all ordered pairs of nonnegative
        integers, starting from an infinite set of roots, where each
        root has an infinite number of children::

            sage: from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest
            sage: S = RecursivelyEnumeratedSet_forest(
            ....:         Family(NN, lambda x: (x, 0)),
            ....:         lambda x: Family(PositiveIntegers(), lambda y: (x[0], y)) if x[1] == 0 else [])
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

            sage: S = RecursivelyEnumeratedSet_forest(Family(NN, lambda x: (x, 0)),
            ....:                                     lambda x: [(x[0], x[1] + 1)])
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
        Apply a Map/Reduce algorithm on ``self``.

        INPUT:

        - ``map_function`` -- a function from the element of ``self`` to some
          set with a reduce operation (e.g.: a monoid). The default value is
          the constant function ``1``.

        - ``reduce_function`` -- the reduce function (e.g.: the addition of a
          monoid); the default value is ``+``

        - ``reduce_init`` -- the initialisation of the reduction (e.g.: the
          neutral element of the monoid); the default value is ``0``

        .. NOTE::

            the effect of the default values is to compute the cardinality
            of ``self``.

        EXAMPLES::

            sage: seeds = [([i], i, i) for i in range(1, 10)]
            sage: def succ(t):
            ....:     list, sum, last = t
            ....:     return [(list + [i], sum + i, i) for i in range(1, last)]
            sage: F = RecursivelyEnumeratedSet(seeds, succ,
            ....:         structure='forest', enumeration='depth')

            sage: # needs sage.symbolic
            sage: y = var('y')
            sage: def map_function(t):
            ....:     li, sum, _ = t
            ....:     return y ^ sum
            sage: def reduce_function(x, y):
            ....:     return x + y
            sage: F.map_reduce(map_function, reduce_function, 0)
            y^45 + y^44 + y^43 + 2*y^42 + 2*y^41 + 3*y^40 + 4*y^39 + 5*y^38 + 6*y^37
            + 8*y^36 + 9*y^35 + 10*y^34 + 12*y^33 + 13*y^32 + 15*y^31 + 17*y^30
            + 18*y^29 + 19*y^28 + 21*y^27 + 21*y^26 + 22*y^25 + 23*y^24 + 23*y^23
            + 23*y^22 + 23*y^21 + 22*y^20 + 21*y^19 + 21*y^18 + 19*y^17 + 18*y^16
            + 17*y^15 + 15*y^14 + 13*y^13 + 12*y^12 + 10*y^11 + 9*y^10 + 8*y^9 + 6*y^8
            + 5*y^7 + 4*y^6 + 3*y^5 + 2*y^4 + 2*y^3 + y^2 + y

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
