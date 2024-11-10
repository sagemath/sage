r"""
Disjoint-set data structure

The main entry point is :func:`DisjointSet` which chooses the appropriate
type to return. For more on the data structure, see :func:`DisjointSet`.

This module defines a class for mutable partitioning of a set, which
cannot be used as a key of a dictionary, a vertex of a graph, etc. For
immutable partitioning see :class:`SetPartition`.

AUTHORS:

- Sébastien Labbé (2008) - Initial version.
- Sébastien Labbé (2009-11-24) - Pickling support
- Sébastien Labbé (2010-01) - Inclusion into sage (:issue:`6775`).
- Giorgos Mousa (2024-04-22): Optimize

EXAMPLES:

Disjoint set of integers from ``0`` to ``n - 1``::

    sage: s = DisjointSet(6)
    sage: s
    {{0}, {1}, {2}, {3}, {4}, {5}}
    sage: s.union(2, 4)
    sage: s.union(1, 3)
    sage: s.union(5, 1)
    sage: s
    {{0}, {1, 3, 5}, {2, 4}}
    sage: s.find(3)
    1
    sage: s.find(5)
    1
    sage: list(map(s.find, range(6)))
    [0, 1, 2, 1, 2, 1]

Disjoint set of hashables objects::

    sage: d = DisjointSet('abcde')
    sage: d
    {{'a'}, {'b'}, {'c'}, {'d'}, {'e'}}
    sage: d.union('a', 'b')
    sage: d.union('b', 'c')
    sage: d.union('c', 'd')
    sage: d
    {{'a', 'b', 'c', 'd'}, {'e'}}
    sage: d.find('c')
    'a'
"""

# ****************************************************************************
#       Copyright (C) 2009 Sébastien Labbé <slabqc at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer cimport Integer
from sage.structure.sage_object cimport SageObject
from cpython.object cimport PyObject_RichCompare
from sage.groups.perm_gps.partn_ref.data_structures cimport *
from sage.misc.lazy_import import LazyImport

SetPartition = LazyImport('sage.combinat.set_partition', 'SetPartition')


cpdef DisjointSet(arg):
    r"""
    Construct a disjoint set where each element of ``arg`` is in its
    own set. If ``arg`` is an integer, then the disjoint set returned is
    made of the integers from ``0`` to ``arg - 1``.

    A disjoint-set data structure (sometimes called union-find data structure)
    is a data structure that keeps track of a partitioning of a set into a
    number of separate, nonoverlapping sets. It performs two operations:

    - :meth:`~sage.sets.disjoint_set.DisjointSet_of_hashables.find` --
      Determine which set a particular element is in.

    - :meth:`~sage.sets.disjoint_set.DisjointSet_of_hashables.union` --
      Combine or merge two sets into a single set.

    REFERENCES:

    - :wikipedia:`Disjoint-set_data_structure`

    INPUT:

    - ``arg`` -- nonnegative integer or an iterable of hashable objects

    EXAMPLES:

    From a nonnegative integer::

        sage: DisjointSet(5)
        {{0}, {1}, {2}, {3}, {4}}

    From an iterable::

        sage: DisjointSet('abcde')
        {{'a'}, {'b'}, {'c'}, {'d'}, {'e'}}
        sage: DisjointSet(range(6))
        {{0}, {1}, {2}, {3}, {4}, {5}}
        sage: DisjointSet(['yi', 45, 'cheval'])
        {{'cheval'}, {'yi'}, {45}}

    From a set partition (see :issue:`38693`)::

        sage: SP = SetPartition(DisjointSet(5))
        sage: DisjointSet(SP)
        {{0}, {1}, {2}, {3}, {4}}
        sage: DisjointSet(SP) == DisjointSet(5)
        True

    TESTS::

        sage: DisjointSet(0)
        {}
        sage: DisjointSet('')
        {}
        sage: DisjointSet([])
        {}

    The argument must be a nonnegative integer::

        sage: DisjointSet(-1)
        Traceback (most recent call last):
        ...
        ValueError: arg must be a nonnegative integer (-1 given)

    or an iterable::

        sage: DisjointSet(4.3)                                                          # needs sage.rings.real_mpfr
        Traceback (most recent call last):
        ...
        TypeError: 'sage.rings.real_mpfr.RealLiteral' object is not iterable

    Moreover, the iterable must consist of hashable::

        sage: DisjointSet([{}, {}])
        Traceback (most recent call last):
        ...
        TypeError: unhashable type: 'dict'
    """
    if isinstance(arg, (Integer, int)):
        if arg < 0:
            raise ValueError('arg must be a nonnegative integer (%s given)' % arg)
        return DisjointSet_of_integers(arg)
    elif isinstance(arg, SetPartition):
        return DisjointSet(arg.base_set())
    else:
        return DisjointSet_of_hashables(arg)

cdef class DisjointSet_class(SageObject):
    r"""
    Common class and methods for :class:`DisjointSet_of_integers` and
    :class:`DisjointSet_of_hashables`.
    """
    def _repr_(self):
        r"""
        Return ``self`` as a unique ``str``.

        EXAMPLES::

            sage: e = DisjointSet(5)
            sage: e.union(2, 4)
            sage: e._repr_()
            '{{0}, {1}, {2, 4}, {3}}'
            sage: e = DisjointSet(5)
            sage: e.union(4, 2)
            sage: e._repr_()
            '{{0}, {1}, {2, 4}, {3}}'

        ::

            sage: e = DisjointSet(range(5))
            sage: e.union(2, 4)
            sage: e._repr_()
            '{{0}, {1}, {2, 4}, {3}}'
            sage: e = DisjointSet(range(5))
            sage: e.union(4, 2)
            sage: e._repr_()
            '{{0}, {1}, {2, 4}, {3}}'
        """
        res = []
        for l in (<dict?>self.root_to_elements_dict()).itervalues():
            l.sort()
            res.append('{%s}' % ', '.join(repr(u) for u in l))
        res.sort()
        return '{%s}' % ', '.join(res)

    def __iter__(self):
        """
        Iterate over elements of the set.

        EXAMPLES::

            sage: d = DisjointSet(4)
            sage: d.union(2, 0)
            sage: sorted(d)
            [[0, 2], [1], [3]]
            sage: d = DisjointSet('abc')
            sage: sorted(d)
            [['a'], ['b'], ['c']]
        """
        return iter((<dict?>self.root_to_elements_dict()).itervalues())

    def __richcmp__(self, other, int op):
        r"""
        Compare the disjoint sets ``self`` and ``other``.

        EXAMPLES::

            sage: d = DisjointSet(5)
            sage: d == d
            True

        ::

            sage: e = DisjointSet(5)
            sage: e == d
            True

        ::

            sage: d.union(0, 3)
            sage: d.union(3, 4)
            sage: e.union(4, 0)
            sage: e.union(3, 0)
            sage: e == d
            True

        ::

            sage: DisjointSet(3) == DisjointSet(5)
            False
            sage: DisjointSet(3) == 4
            False

        ::

            sage: d = DisjointSet('abcde')
            sage: e = DisjointSet('abcde')
            sage: d.union('a', 'b')
            sage: d.union('b', 'c')
            sage: e.union('c', 'a')
            sage: e == d
            False
            sage: e.union('a', 'b')
            sage: e == d
            True
        """
        from sage.sets.set import Set
        s = Set(map(Set, self.root_to_elements_dict().values()))
        try:
            t = Set(map(Set, other.root_to_elements_dict().values()))
        except AttributeError:
            return NotImplemented
        return PyObject_RichCompare(s, t, op)

    def __dealloc__(self):
        r"""
        Deallocate ``self`` (i.e. the ``self._nodes``).

        EXAMPLES::

            sage: d = DisjointSet(5)
            sage: del d
            sage: d = DisjointSet('abc')
            sage: del d
        """
        OP_dealloc(self._nodes)

    def __reduce__(self):
        r"""
        Return a tuple of three elements:

        - The function :func:`DisjointSet`
        - Arguments for the function :func:`DisjointSet`
        - The actual state of ``self``.

        EXAMPLES::

            sage: d = DisjointSet(5)
            sage: d.__reduce__()
            (<built-in function DisjointSet>, (5,), [0, 1, 2, 3, 4])

        ::

            sage: d.union(2, 4)
            sage: d.union(1, 3)
            sage: d.__reduce__()
            (<built-in function DisjointSet>, (5,), [0, 1, 2, 1, 2])
        """
        return DisjointSet, (self._nodes.degree,), self.__getstate__()

    cpdef cardinality(self):
        r"""
        Return the number of elements in ``self``, *not* the number of subsets.

        EXAMPLES::

            sage: d = DisjointSet(5)
            sage: d.cardinality()
            5
            sage: d.union(2, 4)
            sage: d.cardinality()
            5
            sage: d = DisjointSet(range(5))
            sage: d.cardinality()
            5
            sage: d.union(2, 4)
            sage: d.cardinality()
            5
        """
        return self._nodes.degree

    cpdef number_of_subsets(self):
        r"""
        Return the number of subsets in ``self``.

        EXAMPLES::

            sage: d = DisjointSet(5)
            sage: d.number_of_subsets()
            5
            sage: d.union(2, 4)
            sage: d.number_of_subsets()
            4
            sage: d = DisjointSet(range(5))
            sage: d.number_of_subsets()
            5
            sage: d.union(2, 4)
            sage: d.number_of_subsets()
            4
        """
        return self._nodes.num_cells

cdef class DisjointSet_of_integers(DisjointSet_class):
    r"""
    Disjoint set of integers from ``0`` to ``n-1``.

    EXAMPLES::

        sage: d = DisjointSet(5)
        sage: d
        {{0}, {1}, {2}, {3}, {4}}
        sage: d.union(2, 4)
        sage: d.union(0, 2)
        sage: d
        {{0, 2, 4}, {1}, {3}}
        sage: d.find(2)
        2

    TESTS::

        sage: a = DisjointSet(5)
        sage: a == loads(dumps(a))
        True

    ::

        sage: a.union(3, 4)
        sage: a == loads(dumps(a))
        True
    """
    def __init__(self, n):
        r"""
        Construct the ``DisjointSet`` where each element (integers from ``0``
        to ``n-1``) is in its own set.

        INPUT:

        - ``n`` -- nonnegative integer

        EXAMPLES::

            sage: DisjointSet(6)
            {{0}, {1}, {2}, {3}, {4}, {5}}
            sage: DisjointSet(1)
            {{0}}
            sage: DisjointSet(0)
            {}
        """
        self._nodes = OP_new(n)

    def __getstate__(self):
        r"""
        Return a list of the parent of each node from ``0`` to ``n-1``.

        EXAMPLES::

            sage: d = DisjointSet(5)
            sage: d.__getstate__()
            [0, 1, 2, 3, 4]
            sage: d.union(2, 3)
            sage: d.__getstate__()
            [0, 1, 2, 2, 4]
            sage: d.union(3, 0)
            sage: d.__getstate__()
            [2, 1, 2, 2, 4]

        Other parents are obtained when the operations are done in a
        different order::

            sage: d = DisjointSet(5)
            sage: d.union(0, 3)
            sage: d.__getstate__()
            [0, 1, 2, 0, 4]
            sage: d.union(2, 0)
            sage: d.__getstate__()
            [0, 1, 0, 0, 4]
        """
        cdef Py_ssize_t card = self._nodes.degree
        cdef list l = [None] * card
        cdef int i
        for i in range(card):
            l[i] = self._nodes.parent[i]
        return l

    def __setstate__(self, l):
        r"""
        Merge the nodes ``i`` and ``l[i]`` (using union) for ``i`` in
        ``range(len(l))``.

        INPUT:

        - ``l`` -- list of nodes

        EXAMPLES::

            sage: d = DisjointSet(5)
            sage: d.__setstate__([0, 1, 2, 3, 4])
            sage: d
            {{0}, {1}, {2}, {3}, {4}}

        ::

            sage: d = DisjointSet(5)
            sage: d.__setstate__([1, 2, 3, 4, 0])
            sage: d
            {{0, 1, 2, 3, 4}}

        ::

            sage: d = DisjointSet(5)
            sage: d.__setstate__([1, 1, 1])
            sage: d
            {{0, 1, 2}, {3}, {4}}

        ::

            sage: d = DisjointSet(5)
            sage: d.__setstate__([3, 3, 3])
            sage: d
            {{0, 1, 2, 3}, {4}}
        """
        cdef int i, parent
        for i, parent in enumerate(l):
            self.union(parent, i)

    cpdef int find(self, int i) except -1:
        r"""
        Return the representative of the set that ``i`` currently belongs to.

        INPUT:

        - ``i`` -- element in ``self``

        EXAMPLES::

            sage: e = DisjointSet(5)
            sage: e.union(4, 2)
            sage: e
            {{0}, {1}, {2, 4}, {3}}
            sage: e.find(2)
            4
            sage: e.find(4)
            4
            sage: e.union(1, 3)
            sage: e
            {{0}, {1, 3}, {2, 4}}
            sage: e.find(1)
            1
            sage: e.find(3)
            1
            sage: e.union(3, 2)
            sage: e
            {{0}, {1, 2, 3, 4}}
            sage: [e.find(i) for i in range(5)]
            [0, 1, 1, 1, 1]
            sage: e.find(2**10)
            Traceback (most recent call last):
            ...
            ValueError: i must be between 0 and 4 (1024 given)

        .. NOTE::

            This method performs input checks. To avoid them you may directly
            use :meth:`~sage.groups.perm_gps.partn_ref.data_structures.OP_find`.
        """
        card = self._nodes.degree
        if i < 0 or i >= card:
            raise ValueError('i must be between 0 and %s (%s given)' % (card - 1, i))
        return OP_find(self._nodes, i)

    cpdef void union(self, int i, int j) except *:
        r"""
        Combine the set of ``i`` and the set of ``j`` into one.

        All elements in those two sets will share the same representative
        that can be retrieved using
        :meth:`~sage.sets.disjoint_set.DisjointSet_of_integers.find`.

        INPUT:

        - ``i`` -- element in ``self``
        - ``j`` -- element in ``self``

        EXAMPLES::

            sage: d = DisjointSet(5)
            sage: d
            {{0}, {1}, {2}, {3}, {4}}
            sage: d.union(0, 1)
            sage: d
            {{0, 1}, {2}, {3}, {4}}
            sage: d.union(2, 4)
            sage: d
            {{0, 1}, {2, 4}, {3}}
            sage: d.union(1, 4)
            sage: d
            {{0, 1, 2, 4}, {3}}
            sage: d.union(1, 5)
            Traceback (most recent call last):
            ...
            ValueError: j must be between 0 and 4 (5 given)

        .. NOTE::

            This method performs input checks. To avoid them you may directly
            use :meth:`~sage.groups.perm_gps.partn_ref.data_structures.OP_join`.
        """
        cdef int card = self._nodes.degree
        if i < 0 or i >= card:
            raise ValueError('i must be between 0 and %s (%s given)' % (card - 1, i))
        if j < 0 or j >= card:
            raise ValueError('j must be between 0 and %s (%s given)' % (card - 1, j))
        OP_join(self._nodes, i, j)

    def make_set(self):
        r"""
        Add a new element into a new set containing only the new element.

        According to :wikipedia:`Disjoint-set_data_structure#Making_new_sets` the
        `make_set` operation adds a new element into a new set containing only
        the new element. The new set is added at the end of `self`.

        EXAMPLES::
            sage: d = DisjointSet(5)
            sage: d.union(1, 2)
            sage: d.union(0, 1)
            sage: d.make_set()
            sage: d
            {{0, 1, 2}, {3}, {4}, {5}}
            sage: d.find(1)
            1
        """
        OP_make_set(self._nodes)

    cpdef root_to_elements_dict(self):
        r"""
        Return the dictionary where the keys are the roots of ``self`` and the
        values are the elements in the same set as the root.

        EXAMPLES::

            sage: d = DisjointSet(5)
            sage: sorted(d.root_to_elements_dict().items())
            [(0, [0]), (1, [1]), (2, [2]), (3, [3]), (4, [4])]
            sage: d.union(2, 3)
            sage: sorted(d.root_to_elements_dict().items())
            [(0, [0]), (1, [1]), (2, [2, 3]), (4, [4])]
            sage: d.union(3, 0)
            sage: sorted(d.root_to_elements_dict().items())
            [(1, [1]), (2, [0, 2, 3]), (4, [4])]
            sage: d
            {{0, 2, 3}, {1}, {4}}
        """
        cdef dict s = {}
        cdef int i, o
        for i in range(self._nodes.degree):
            o = OP_find(self._nodes, i)
            if o not in s:
                s[o] = []
            s[o].append(i)
        return s

    cpdef element_to_root_dict(self):
        r"""
        Return the dictionary where the keys are the elements of ``self`` and
        the values are their representative inside a list.

        EXAMPLES::

            sage: d = DisjointSet(5)
            sage: d.union(2, 3)
            sage: d.union(4, 1)
            sage: e = d.element_to_root_dict()
            sage: e
            {0: 0, 1: 4, 2: 2, 3: 2, 4: 4}
            sage: WordMorphism(e)                                                       # needs sage.combinat
            WordMorphism: 0->0, 1->4, 2->2, 3->2, 4->4
        """
        cdef dict d = {}
        cdef int i
        for i in range(self._nodes.degree):
            d[i] = OP_find(self._nodes, i)
        return d

    cpdef to_digraph(self):
        r"""
        Return the current digraph of ``self`` where `(a, b)` is an oriented
        edge if `b` is the parent of `a`.

        EXAMPLES::

            sage: d = DisjointSet(5)
            sage: d.union(2, 3)
            sage: d.union(4, 1)
            sage: d.union(3, 4)
            sage: d
            {{0}, {1, 2, 3, 4}}
            sage: g = d.to_digraph(); g                                                 # needs sage.graphs
            Looped digraph on 5 vertices
            sage: g.edges(sort=True)                                                    # needs sage.graphs
            [(0, 0, None), (1, 2, None), (2, 2, None), (3, 2, None), (4, 2, None)]

        The result depends on the ordering of the union::

            sage: d = DisjointSet(5)
            sage: d.union(1, 2)
            sage: d.union(1, 3)
            sage: d.union(1, 4)
            sage: d
            {{0}, {1, 2, 3, 4}}
            sage: d.to_digraph().edges(sort=True)                                       # needs sage.graphs
            [(0, 0, None), (1, 1, None), (2, 1, None), (3, 1, None), (4, 1, None)]
        """
        cdef dict d = {i: [self._nodes.parent[i]] for i in range(self._nodes.degree)}
        from sage.graphs.digraph import DiGraph
        return DiGraph(d)

cdef class DisjointSet_of_hashables(DisjointSet_class):
    r"""
    Disjoint set of hashables.

    EXAMPLES::

        sage: d = DisjointSet('abcde')
        sage: d
        {{'a'}, {'b'}, {'c'}, {'d'}, {'e'}}
        sage: d.union('a', 'c')
        sage: d
        {{'a', 'c'}, {'b'}, {'d'}, {'e'}}
        sage: d.find('a')
        'a'

    TESTS::

        sage: a = DisjointSet('abcdef')
        sage: a == loads(dumps(a))
        True

    ::

        sage: a.union('a', 'c')
        sage: a == loads(dumps(a))
        True
    """
    def __init__(self, iterable):
        r"""
        Construct the trivial disjoint set where each element is in its own set.

        INPUT:

        - ``iterable`` -- iterable of hashable objects

        EXAMPLES::

            sage: DisjointSet('abcde')
            {{'a'}, {'b'}, {'c'}, {'d'}, {'e'}}
            sage: DisjointSet(range(6))
            {{0}, {1}, {2}, {3}, {4}, {5}}
            sage: DisjointSet(['yi', 45, 'cheval'])
            {{'cheval'}, {'yi'}, {45}}
            sage: DisjointSet(set([0, 1, 2, 3, 4]))
            {{0}, {1}, {2}, {3}, {4}}
        """
        cdef int i
        self._int_to_el = []
        self._el_to_int = {}
        for i, e in enumerate(iterable):
            self._int_to_el.append(e)
            self._el_to_int[e] = i
        self._nodes = OP_new(len(self._int_to_el))

    def __reduce__(self):
        r"""
        Return a tuple of three elements:

        - The function :func:`DisjointSet`
        - Arguments for the function :func:`DisjointSet`
        - The actual state of ``self``.

        EXAMPLES::

            sage: DisjointSet(range(5))
            {{0}, {1}, {2}, {3}, {4}}
            sage: d = _
            sage: d.__reduce__()
            (<built-in function DisjointSet>,
             ([0, 1, 2, 3, 4],),
             [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)])

         ::

            sage: d.union(2, 4)
            sage: d.union(1, 3)
            sage: d.__reduce__()
            (<built-in function DisjointSet>,
             ([0, 1, 2, 3, 4],),
             [(0, 0), (1, 1), (2, 2), (3, 1), (4, 2)])
        """
        return DisjointSet, (self._int_to_el,), self.__getstate__()

    def __getstate__(self):
        r"""
        Return a list of pairs (``n``, parent of ``n``) for each node ``n``.

        EXAMPLES::

            sage: d = DisjointSet('abcde')
            sage: d.__getstate__()
            [('a', 'a'), ('b', 'b'), ('c', 'c'), ('d', 'd'), ('e', 'e')]
            sage: d.union('c', 'd')
            sage: d.__getstate__()
            [('a', 'a'), ('b', 'b'), ('c', 'c'), ('d', 'c'), ('e', 'e')]
            sage: d.union('d', 'a')
            sage: d.__getstate__()
            [('a', 'c'), ('b', 'b'), ('c', 'c'), ('d', 'c'), ('e', 'e')]

        Other parents are obtained when the operations are done in a
        different order::

            sage: d = DisjointSet('abcde')
            sage: d.union('d', 'c')
            sage: d.__getstate__()
            [('a', 'a'), ('b', 'b'), ('c', 'd'), ('d', 'd'), ('e', 'e')]
        """
        cdef int card = self._nodes.degree
        cdef list l = [None] * card
        cdef int i
        for i in range(card):
            l[i] = self._int_to_el[self._nodes.parent[i]]
        return list(zip(self._int_to_el, l))

    def __setstate__(self, l):
        r"""
        Merge the nodes ``a`` and ``b`` for each pair of nodes
        ``(a, b)`` in ``l``.

        INPUT:

        - ``l`` -- list of pair of nodes

        EXAMPLES::

            sage: d = DisjointSet('abcde')
            sage: d.__setstate__([('a', 'a'), ('b', 'b'), ('c', 'c')])
            sage: d
            {{'a'}, {'b'}, {'c'}, {'d'}, {'e'}}

        ::

            sage: d = DisjointSet('abcde')
            sage: d.__setstate__([('a', 'b'), ('b', 'c'), ('c', 'd'), ('d', 'e')])
            sage: d
            {{'a', 'b', 'c', 'd', 'e'}}
        """
        for a, b in l:
            self.union(a, b)

    cpdef find(self, e):
        r"""
        Return the representative of the set that ``e`` currently belongs to.

        INPUT:

        - ``e`` -- element in ``self``

        EXAMPLES::

            sage: e = DisjointSet(range(5))
            sage: e.union(4, 2)
            sage: e
            {{0}, {1}, {2, 4}, {3}}
            sage: e.find(2)
            4
            sage: e.find(4)
            4
            sage: e.union(1,3)
            sage: e
            {{0}, {1, 3}, {2, 4}}
            sage: e.find(1)
            1
            sage: e.find(3)
            1
            sage: e.union(3, 2)
            sage: e
            {{0}, {1, 2, 3, 4}}
            sage: [e.find(i) for i in range(5)]
            [0, 1, 1, 1, 1]
            sage: e.find(5)
            Traceback (most recent call last):
            ...
            KeyError: 5
        """
        cdef int i = <int> self._el_to_int[e]
        cdef int r = <int> OP_find(self._nodes, i)
        return self._int_to_el[r]

    cpdef void union(self, e, f) noexcept:
        r"""
        Combine the set of ``e`` and the set of ``f`` into one.

        All elements in those two sets will share the same representative
        that can be retrieved using
        :meth:`~sage.sets.disjoint_set.DisjointSet_of_hashables.find`.

        INPUT:

        - ``e`` -- element in ``self``
        - ``f`` -- element in ``self``

        EXAMPLES::

            sage: e = DisjointSet('abcde')
            sage: e
            {{'a'}, {'b'}, {'c'}, {'d'}, {'e'}}
            sage: e.union('a', 'b')
            sage: e
            {{'a', 'b'}, {'c'}, {'d'}, {'e'}}
            sage: e.union('c', 'e')
            sage: e
            {{'a', 'b'}, {'c', 'e'}, {'d'}}
            sage: e.union('b', 'e')
            sage: e
            {{'a', 'b', 'c', 'e'}, {'d'}}
            sage: e.union('a', 2**10)
            KeyError: 1024
            ...
        """
        cdef int i = <int> self._el_to_int[e]
        cdef int j = <int> self._el_to_int[f]
        OP_join(self._nodes, i, j)

    def make_set(self, new_elt=None):
        r"""
        Add a new element into a new set containing only the new element.

        According to :wikipedia:`Disjoint-set_data_structure#Making_new_sets`
        the `make_set` operation adds a new element into a new set containing
        only the new element. The new set is added at the end of `self`.

        INPUT:

        - ``new_elt`` -- (optional) element to add. If `None`, then an integer
          is added.

        EXAMPLES::

            sage: e = DisjointSet('abcde')
            sage: e.union('d', 'c')
            sage: e.union('c', 'e')
            sage: e.make_set('f')
            sage: e
            {{'a'}, {'b'}, {'c', 'd', 'e'}, {'f'}}
            sage: e.union('f', 'b')
            sage: e
            {{'a'}, {'b', 'f'}, {'c', 'd', 'e'}}
            sage: e.make_set('e'); e
            {{'a'}, {'b', 'f'}, {'c', 'd', 'e'}}
            sage: e.make_set(); e
            {{'a'}, {'b', 'f'}, {'c', 'd', 'e'}, {6}}
        """
        if new_elt is None:
            new_elt = self._nodes.degree
        if new_elt not in self._int_to_el:
            d = self._nodes.degree
            self._int_to_el.append(new_elt)
            self._el_to_int[new_elt] = d
            OP_make_set(self._nodes)

    cpdef root_to_elements_dict(self):
        r"""
        Return the dictionary where the keys are the roots of ``self`` and the
        values are the elements in the same set.

        EXAMPLES::

            sage: d = DisjointSet(range(5))
            sage: d.union(2, 3)
            sage: d.union(4, 1)
            sage: e = d.root_to_elements_dict()
            sage: sorted(e.items())
            [(0, [0]), (2, [2, 3]), (4, [1, 4])]
        """
        cdef dict s = {}
        for e in self._int_to_el:
            r = self.find(e)
            if r not in s:
                s[r] = []
            s[r].append(e)
        return s

    cpdef element_to_root_dict(self):
        r"""
        Return the dictionary where the keys are the elements of ``self`` and
        the values are their representative inside a list.

        EXAMPLES::

            sage: d = DisjointSet(range(5))
            sage: d.union(2, 3)
            sage: d.union(4, 1)
            sage: e = d.element_to_root_dict()
            sage: sorted(e.items())
            [(0, 0), (1, 4), (2, 2), (3, 2), (4, 4)]
            sage: WordMorphism(e)                                                       # needs sage.combinat
            WordMorphism: 0->0, 1->4, 2->2, 3->2, 4->4
        """
        cdef dict d = {}
        for a in self._int_to_el:
            d[a] = self.find(a)
        return d

    cpdef to_digraph(self):
        r"""
        Return the current digraph of ``self`` where `(a, b)` is an oriented
        edge if `b` is the parent of `a`.

        EXAMPLES::

            sage: d = DisjointSet(range(5))
            sage: d.union(2, 3)
            sage: d.union(4, 1)
            sage: d.union(3, 4)
            sage: d
            {{0}, {1, 2, 3, 4}}
            sage: g = d.to_digraph()
            sage: g                                                                     # needs sage.graphs
            Looped digraph on 5 vertices
            sage: g.edges(sort=True)                                                    # needs sage.graphs
            [(0, 0, None), (1, 2, None), (2, 2, None), (3, 2, None), (4, 2, None)]

        The result depends on the ordering of the union::

            sage: d = DisjointSet(range(5))
            sage: d.union(1, 2)
            sage: d.union(1, 3)
            sage: d.union(1, 4)
            sage: d
            {{0}, {1, 2, 3, 4}}
            sage: d.to_digraph().edges(sort=True)                                       # needs sage.graphs
            [(0, 0, None), (1, 1, None), (2, 1, None), (3, 1, None), (4, 1, None)]
        """
        cdef dict d = {}
        cdef int i
        for i in range(self._nodes.degree):
            e = self._int_to_el[i]
            p = self._int_to_el[self._nodes.parent[i]]
            d[e] = [p]
        from sage.graphs.digraph import DiGraph
        return DiGraph(d)
