r"""
Recursively enumerated set

A set `S` is called recursively enumerable if there is an algorithm that
enumerates the members of `S`. We consider here the recursively enumerated
sets that are described by some ``seeds`` and a successor function
``successors``.  The successor function may have some structure (symmetric,
graded, forest) or not. The elements of a set having a symmetric, graded or
forest structure can be enumerated uniquely without keeping all of them in
memory. Many kinds of iterators are provided in this module: depth first
search, breadth first search or elements of given depth.

See :wikipedia:`Recursively_enumerable_set`.

See documentation of :func:`RecursivelyEnumeratedSet` below for the
description of the inputs.

AUTHORS:

- Sebastien Labbe, April 2014, at Sage Days 57, Cernay-la-ville

EXAMPLES:

Forest structure
----------------

The set of words over the alphabet `\{a,b\}` can be generated from the
empty word by appending letter `a` or `b` as a successor function. This set
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

Symmetric structure
-------------------

The origin ``(0, 0)`` as seed and the upper, lower, left and right lattice
point as successor function. This function is symmetric::

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
    sage: [next(it_breadth) for _ in range(10)]
    [(0, 0), (0, 1), (0, -1), (1, 0), (-1, 0), (-1, 1), (-2, 0), (0, 2), (2, 0), (-1, -1)]

Levels (elements of given depth)::

    sage: sorted(C.graded_component(0))
    [(0, 0)]
    sage: sorted(C.graded_component(1))
    [(-1, 0), (0, -1), (0, 1), (1, 0)]
    sage: sorted(C.graded_component(2))
    [(-2, 0), (-1, -1), (-1, 1), (0, -2), (0, 2), (1, -1), (1, 1), (2, 0)]

Graded structure
----------------

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
     [1, 3, 2, 4, 5],
     [1, 2, 4, 3, 5],
     [2, 1, 3, 4, 5],
     [1, 2, 3, 5, 4]]

Elements of given depth iterator::

    sage: list(R.elements_of_depth_iterator(9))
    [[5, 3, 4, 2, 1], [4, 5, 3, 2, 1], [5, 4, 2, 3, 1], [5, 4, 3, 1, 2]]
    sage: list(R.elements_of_depth_iterator(10))
    [[5, 4, 3, 2, 1]]

Graded components (set of elements of the same depth)::

    sage: sorted(R.graded_component(0))
    [[1, 2, 3, 4, 5]]
    sage: sorted(R.graded_component(1))
    [[1, 2, 3, 5, 4], [1, 2, 4, 3, 5], [1, 3, 2, 4, 5], [2, 1, 3, 4, 5]]
    sage: sorted(R.graded_component(9))
    [[4, 5, 3, 2, 1], [5, 3, 4, 2, 1], [5, 4, 2, 3, 1], [5, 4, 3, 1, 2]]
    sage: sorted(R.graded_component(10))
    [[5, 4, 3, 2, 1]]

No hypothesis on the structure
------------------------------

By "no hypothesis" is meant neither a forest, neither symmetric neither
graded, it may have other structure like not containing oriented cycle but
this does not help for enumeration.

In this example, the seed is 0 and the successor function is either ``+2``
or ``+3``. This is the set of non negative linear combinations of 2 and 3::

    sage: succ = lambda a:[a+2,a+3]
    sage: C = RecursivelyEnumeratedSet([0], succ)
    sage: C
    A recursively enumerated set (breadth first search)

Breadth first search::

    sage: it = C.breadth_first_search_iterator()
    sage: [next(it) for _ in range(10)]
    [0, 2, 3, 4, 5, 6, 8, 9, 7, 10]

Depth first search::

    sage: it = C.depth_first_search_iterator()
    sage: [next(it) for _ in range(10)]
    [0, 3, 6, 9, 12, 15, 18, 21, 24, 27]

"""
#*****************************************************************************
#       Copyright (C) 2014 Sebastien Labbe <slabqc at gmail.com>
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
from sage.structure.parent cimport Parent
from sage.categories.enumerated_sets import EnumeratedSets
from sage.combinat.backtrack import SearchForest
#from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from collections import deque

def RecursivelyEnumeratedSet(seeds, successors, structure=None,
            enumeration=None, max_depth=float("inf"), post_process=None,
            facade=None, category=None):
    r"""
    Return a recursively enumerated set.

    A set `S` is called recursively enumerable if there is an algorithm that
    enumerates the members of `S`. We consider here the recursively
    enumerated set that are described by some ``seeds`` and a successor
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
    - ``structure`` -- string (optional, default: ``None``), structure of the
      set, possible values are:

      - ``None`` -- nothing is known about the structure of the set.
      - ``'forest'`` -- if the ``successors`` function generates a *forest*, that
        is, each element can be reached uniquely from a seed.
      - ``'graded'`` -- if the ``successors`` function is *graded*, that is, all
        paths from a seed to a given element have equal length.
      - ``'symmetric'`` -- if the relation is *symmetric*, that is,
        ``y in successors(x)`` if and only if ``x in successors(y)``

    - ``enumeration`` -- ``'depth'``, ``'breadth'``, ``'naive'`` or ``None``
      (optional, default: ``None``). The default enumeration for the
      ``__iter__`` function.
    - ``max_depth`` -- integer (optional, default: ``float("inf")``), limit
      the search to a certain depth, currently works only for breadth first
      search
    - ``post_process`` -- (optional, default: ``None``), for forest only
    - ``facade`` -- (optional, default: ``None``)
    - ``category`` -- (optional, default: ``None``)

    EXAMPLES:

    A recursive set with no other information::

        sage: f = lambda a: [a+3, a+5]
        sage: C = RecursivelyEnumeratedSet([0], f)
        sage: C
        A recursively enumerated set (breadth first search)
        sage: it = iter(C)
        sage: [next(it) for _ in range(10)]
        [0, 3, 5, 8, 10, 6, 9, 11, 13, 15]

    A recursive set with a forest structure::

        sage: f = lambda a: [2*a,2*a+1]
        sage: C = RecursivelyEnumeratedSet([1], f, structure='forest')
        sage: C
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
        [10, 15, 16, 9, 11, 14, 8]

    A recursive set given by a graded relation::

        sage: f = lambda a: [a+1, a+I]
        sage: C = RecursivelyEnumeratedSet([0], f, structure='graded')
        sage: C
        A recursively enumerated set with a graded structure (breadth first search)
        sage: it = iter(C)
        sage: [next(it) for _ in range(7)]
        [0, 1, I, I + 1, 2, 2*I, I + 2]

    .. WARNING::

        If you do not set the good structure, you might obtain bad results,
        like elements generated twice::

            sage: f = lambda a: [a-1,a+1]
            sage: C = RecursivelyEnumeratedSet([0], f, structure='graded')
            sage: it = iter(C)
            sage: [next(it) for _ in range(7)]
            [0, 1, -1, 0, 2, -2, 1]

    TESTS:

    The succesors method is an attribute::

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
        if enumeration is None: enumeration = 'breadth'
        return RecursivelyEnumeratedSet_generic(seeds, successors,
                enumeration, max_depth, facade=facade, category=category)
    if structure == 'symmetric':
        if enumeration is None: enumeration = 'breadth'
        return RecursivelyEnumeratedSet_symmetric(seeds, successors,
                enumeration, max_depth, facade=facade, category=category)
    if structure == 'forest':
        if enumeration is None: enumeration = 'depth'
        return RecursivelyEnumeratedSet_forest(roots=seeds, children=successors,
                algorithm=enumeration, post_process=post_process,
                facade=facade, category=category)
    if structure == 'graded':
        if enumeration is None: enumeration = 'breadth'
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
        # Since trac ticket #21312, the graded component iterator is not used
        # anymore but maybe some previously pickled object still have it
        # self._graded_component_it = l[1]

    def __len__(self):
        """
        Disable ``__len__()`` from :class:`Parent` :trac:`12955`.

        Because Python assumes ``__len__()`` is fast and we cannot
        have a fast default implementation.

        EXAMPLES::

            sage: f = lambda a: [a+3, a+5]
            sage: C = RecursivelyEnumeratedSet([0], f)
            sage: len(C)
            Traceback (most recent call last):
            ...
            TypeError: 'NoneType' object cannot be interpreted as an index
        """
        return None

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
            sage: [next(it_naive) for _ in range(10)]
            [0, 3, 8, 11, 5, 6, 9, 10, 12, 13]
            sage: [next(it_depth) for _ in range(10)]
            [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
            sage: [next(it_breadth) for _ in range(10)]
            [0, 3, 5, 8, 10, 6, 9, 11, 13, 15]
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

            sage: R = RecursivelyEnumeratedSet([1], lambda x: [x+1, x-1])
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

        OUTPUT:

        An iterator of sets.

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

        OUTPUT:

        A set.

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

        An element of depth `n` can be obtained applying `n` times the
        successor function to a seed.

        INPUT:

        - ``depth`` -- integer

        OUTPUT:

        An iterator.

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

        INPUT:

        - ``max_depth`` -- (Default: ``None``) specifies the maximal depth
          to which elements are computed; if ``None``, the value of
          ``self._max_depth`` is used

        EXAMPLES::

            sage: f = lambda a: [a+3, a+5]
            sage: C = RecursivelyEnumeratedSet([0], f)
            sage: it = C.breadth_first_search_iterator()
            sage: [next(it) for _ in range(10)]
            [0, 3, 5, 8, 10, 6, 9, 11, 13, 15]
        """
        if max_depth is None:
            max_depth = self._max_depth
        current_level = self._seeds
        known = set(current_level)
        depth = 0
        while current_level and depth <= max_depth:
            next_level = set()
            for x in current_level:
                yield x
                for y in self.successors(x):
                    if y is None or y in known:
                        continue
                    next_level.add(y)
                    known.add(y)
            current_level = next_level
            depth += 1

    def _breadth_first_search_iterator_from_graded_component_iterator(self, max_depth=None):
        r"""
        Iterate on the elements of ``self`` (breadth first).

        This code remembers only elements needed by the graded component iterator
        to generate the next graded component.

        This method is the default breadth first search iterator when the
        structure is symmetric or graded.

        INPUT:

        - ``max_depth`` -- (Default: ``None``) specifies the maximal depth
          to which elements are computed; if ``None``, the value of
          ``self._max_depth`` is used

        .. NOTE::

            Calling next in this iterator will be either quite slow or very fast
            since it generates the whole graded component before yielding the
            elements of each graded component.

        EXAMPLES::

            sage: f = lambda a: [(a[0]+1,a[1]), (a[0],a[1]+1)]
            sage: C = RecursivelyEnumeratedSet([(0,0)], f, structure='graded')
            sage: it = C._breadth_first_search_iterator_from_graded_component_iterator(max_depth=3)
            sage: list(it)
            [(0, 0), (0, 1), (1, 0), (2, 0), (1, 1), (0, 2)]

        This iterator is used by default for symmetric structure::

            sage: f = lambda a: [a-1,a+1]
            sage: S = RecursivelyEnumeratedSet([10], f, structure='symmetric')
            sage: it = iter(S)
            sage: [next(it) for _ in range(7)]
            [10, 9, 11, 8, 12, 13, 7]
        """
        if max_depth is None:
            max_depth = self._max_depth
        it = self.graded_component_iterator()
        cdef int i = 0
        while i < max_depth:
            graded_component = next(it)
            for a in graded_component:
                yield a
            i += 1

    def _breadth_first_search_iterator_using_queue(self):
        r"""
        Iterate on the elements of ``self`` (breadth first).

        This code remembers every elements generated and uses python
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

        This code remembers every elements generated.

        TESTS:

        We compute all the permutations of 3::

            sage: seeds = [Permutation([1,2,3])]
            sage: succ = attrcall("permutohedron_succ")
            sage: R = RecursivelyEnumeratedSet(seeds, succ)
            sage: list(R.naive_search_iterator())
            [[1, 2, 3], [2, 1, 3], [1, 3, 2], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
        """
        cdef set known, todo
        known = set(self._seeds)
        todo = known.copy()
        while todo:
            x = todo.pop()
            yield x
            for y in self.successors(x):
                if y == None or y in known:
                    continue
                todo.add(y)
                known.add(y)

    def depth_first_search_iterator(self):
        r"""
        Iterate on the elements of ``self`` (depth first).

        This code remembers every elements generated.

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

        - ``max_depth`` -- (default: ``None``) specifies the maximal depth
          for which outgoing edges of elements are computed; if ``None``, the
          value of ``self._max_depth`` is used
        - ``loops`` -- (default: ``True``) option for the digraph
        - ``multiedges`` -- (default: ``True``) option of the digraph

        OUTPUT:

        A directed graph

        .. WARNING::

            If the set is infinite, this will loop forever unless ``max_depth``
            is finite.

        EXAMPLES::

            sage: child = lambda i: [(i+3) % 10, (i+8) % 10]
            sage: R = RecursivelyEnumeratedSet([0], child)
            sage: R.to_digraph()
            Looped multi-digraph on 10 vertices

        Digraph of an recursively enumerated set with a symmetric structure of
        infinite cardinality using ``max_depth`` argument::

            sage: succ = lambda a: [(a[0]-1,a[1]), (a[0],a[1]-1), (a[0]+1,a[1]), (a[0],a[1]+1)]
            sage: seeds = [(0,0)]
            sage: C = RecursivelyEnumeratedSet(seeds, succ, structure='symmetric')
            sage: C.to_digraph(max_depth=4)
            Looped multi-digraph on 41 vertices

        The ``max_depth`` argument can be given at the creation of the set::

            sage: C = RecursivelyEnumeratedSet(seeds, succ, structure='symmetric', max_depth=3)
            sage: C.to_digraph()
            Looped multi-digraph on 25 vertices

        Digraph of an recursively enumerated set with a graded structure::

            sage: f = lambda a: [a+1, a+I]
            sage: C = RecursivelyEnumeratedSet([0], f, structure='graded')
            sage: C.to_digraph(max_depth=4)
            Looped multi-digraph on 21 vertices
        """
        successors = self.successors
        it = self.breadth_first_search_iterator(max_depth=max_depth)
        E = [(u,v) for u in it for v in successors(u)]
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
        [0, 1, -1, 2, -2, 3, -3]

    TESTS:

    Do not use lambda functions for saving purposes::

        sage: f = lambda a: [a-1,a+1]
        sage: C = RecursivelyEnumeratedSet([0], f, structure='symmetric')
        sage: loads(dumps(C))
        Traceback (most recent call last):
        ...
        PicklingError: Can't pickle <type 'function'>: attribute lookup __builtin__.function failed

    This works in the command line but apparently not as a doctest::

        sage: def f(a): return [a-1,a+1]
        sage: C = RecursivelyEnumeratedSet([0], f, structure='symmetric')
        sage: loads(dumps(C))
        Traceback (most recent call last):
        ...
        PicklingError: Can't pickle <type 'function'>: attribute lookup __builtin__.function failed
    """
    breadth_first_search_iterator = RecursivelyEnumeratedSet_generic._breadth_first_search_iterator_from_graded_component_iterator

    def graded_component_iterator(self):
        r"""
        Iterate over the graded components of ``self``.

        A graded component is a set of elements of the same depth.

        The enumeration remembers only the last two graded components
        generated since the structure is symmetric.

        OUTPUT:

        An iterator of sets.

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

            sage: f = lambda a: [a+1, a+I]
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

            sage: def f(a):
            ....:     sleep(0.05)
            ....:     return [a-1,a+1]
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
        cdef set A,B
        A = set()
        B = set(self._seeds)
        while B:
            yield B
            A,B = B, self._get_next_graded_component(A, B)

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

        OUTPUT:

        A set.

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

        We make sure that :trac:`21312` is fixed::

            sage: def f(a):
            ....:    sleep(0.1)
            ....:    return [a-1,a+1]
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

        OUTPUT:

        - ``C`` -- set, the set of elements of depth n+1

        .. TODO::

            Can ``collections.OrderedDict`` can help maintain the breadth
            first search enumeration for each graded component?

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
        sage: sorted(C)
        [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0),
         (1, 1), (1, 2), (2, 0), (2, 1), (3, 0)]
    """
    def breadth_first_search_iterator(self, max_depth=None):
        r"""
        Iterate on the elements of ``self`` (breadth first).

        This iterator make use of the graded structure by remembering only
        the elements of the current depth.

        INPUT:

        - ``max_depth`` -- (Default: ``None``) Specifies the maximal depth
          to which elements are computed. If None, the value of
          ``self._max_depth`` is used.

        EXAMPLES::

            sage: f = lambda a: [(a[0]+1,a[1]), (a[0],a[1]+1)]
            sage: C = RecursivelyEnumeratedSet([(0,0)], f, structure='graded')
            sage: it = C.breadth_first_search_iterator(max_depth=3)
            sage: list(it)
            [(0, 0), (0, 1), (1, 0), (2, 0), (1, 1),
             (0, 2), (3, 0), (1, 2), (0, 3), (2, 1)]
        """
        cdef set next_level
        cdef int depth
        if max_depth is None:
            max_depth = self._max_depth
        current_level = self._seeds
        depth = 0
        while current_level and depth <= max_depth:
            next_level = set()
            for x in current_level:
                yield x
                for y in self.successors(x):
                    if y is None or y in next_level:
                        continue
                    next_level.add(y)
            current_level = next_level
            depth += 1

    def graded_component_iterator(self):
        r"""
        Iterate over the graded components of ``self``.

        A graded component is a set of elements of the same depth.

        The algorithm remembers only the current graded component generated
        since the structure is graded.

        OUTPUT:

        An iterator of sets.

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

        Make sure that :trac:`20225` is fixed::

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

        OUTPUT:

        A set.

        EXAMPLES::

            sage: f = lambda a: [a+1, a+I]
            sage: C = RecursivelyEnumeratedSet([0], f, structure='graded')
            sage: for i in range(5): sorted(C.graded_component(i))
            [0]
            [I, 1]
            [2*I, I + 1, 2]
            [3*I, 2*I + 1, I + 2, 3]
            [4*I, 3*I + 1, 2*I + 2, I + 3, 4]

        TESTS:

        We make sure that :trac:`21312` is fixed::

            sage: def f(a):
            ....:    sleep(0.1)
            ....:    return [a+1, a+I]
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

        OUTPUT:

        - ``C`` -- set, the set of elements of depth `n+1`

        .. TODO::

            Can ``collections.OrderedDict`` can help maintain the breadth
            first search enumeration for each graded component?

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

RecursivelyEnumeratedSet_forest = SearchForest
