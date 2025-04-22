r"""
Disjoint-set data structure

The main entry point is :func:`DisjointSet` which chooses the appropriate
type to return. For more on the data structure, see :func:`DisjointSet`.

This module defines a class for mutable partitioning of a set, which
cannot be used as a key of a dictionary, a vertex of a graph, etc. For
immutable partitioning see :class:`SetPartition`.

Behavior based on input ``arg``:

* ``arg=None`` (or omitted): Returns an empty dynamic disjoint set.
  Elements can be added later using the ``union`` method.
  ``dynamic=True`` is implied.
* ``arg`` is a non-negative integer `n`: Returns a static disjoint
  set containing integers `0` to `n-1`. Elements *cannot* be added
  later via ``union`` or ``make_set``. The ``dynamic`` flag is ignored.
* ``arg`` is an iterable (e.g., list, string, set):
    * If ``dynamic=False`` (default): Returns a static disjoint set
      containing the elements from the iterable. ``union`` will raise a
      ``KeyError`` if called with elements not initially present.
      ``make_set`` will raise a ``TypeError``.
    * If ``dynamic=True``: Returns a dynamic disjoint set containing
      the elements from the iterable. ``union`` will automatically add
      new elements if they are not present. ``make_set`` can also be used.
* ``arg`` is a `SetPartition`: Creates a disjoint set from its base
  set, respecting the ``dynamic`` flag.

Static vs. Dynamic:

* Static sets (default for iterables, always for integers) are
  optimized for cases where the ground set of elements is fixed. They
  provide slightly better performance and enforce the fixed-set
  constraint by raising errors if you attempt to add elements after
  creation (via ``union`` or ``make_set``).
* Dynamic sets (``dynamic=True`` or created empty) are designed for
  building the set structure incrementally. The ``union`` operation
  conveniently adds any missing elements it encounters.

Use ``dynamic=True`` when the set of elements is not known in advance or
when you prefer the convenience of implicit element addition during unions.
Otherwise, use the default static behavior for performance and stricter
set definition.

You can check if an element is part of the ground set using the ``in``
operator (e.g., ``if element in my_disjoint_set: ...``).


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

    sage: d = DisjointSet('abcde')  # Static by default now
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
from sage.sets.set import Set
from cysignals.memory cimport sig_calloc, sig_free

SetPartition = LazyImport('sage.combinat.set_partition', 'SetPartition')


cpdef DisjointSet(arg=None, dynamic=False):
    r"""
    Construct a disjoint set data structure (union-find).

    This function acts as a factory, returning the appropriate specialized
    DisjointSet object based on the input arguments and the ``dynamic`` flag.

    REFERENCES:

    - :wikipedia:`Disjoint-set_data_structure`

    INPUT:

    - ``arg`` -- (default: ``None``) Nonnegative integer, an iterable of
      hashable objects, a ``SetPartition``, or ``None``.
    - ``dynamic`` -- boolean (default: ``False``); If ``arg`` is an iterable
      or ``None``, determines whether to create a dynamic or static disjoint
      set. Ignored if ``arg`` is an integer.

    EXAMPLES::

    Empty (dynamic) set::

        sage: D = DisjointSet()
        sage: D
        {}
        sage: 1 in D
        False
        sage: D.union(1, 2)
        sage: D
        {{1, 2}}
        sage: 1 in D
        True
        sage: D.make_set(3)  # Explicitly add 3
        sage: D
        {{1, 2}, {3}}
        sage: isinstance(D, sage.sets.disjoint_set.DynamicDisjointSet_of_hashables)
        True

    Static set from integer::

        sage: S_int = DisjointSet(5)
        sage: S_int
        {{0}, {1}, {2}, {3}, {4}}
        sage: 3 in S_int
        True
        sage: 6 in S_int
        False
        sage: S_int.union(1, 6)
        Traceback (most recent call last):
        ...
        ValueError: j must be between 0 and 4 (6 given)
        sage: S_int.make_set()
        Traceback (most recent call last):
        ...
        TypeError: Cannot add elements to a static DisjointSet_of_integers.
        sage: isinstance(S_int, sage.sets.disjoint_set.DisjointSet_of_integers)
        True

    Static set from iterable (default)::

        sage: S_static = DisjointSet('abc')
        sage: S_static
        {{'a'}, {'b'}, {'c'}}
        sage: 'a' in S_static
        True
        sage: 'd' in S_static
        False
        sage: S_static.union('a', 'd')
        Traceback (most recent call last):
        ...
        KeyError: "Element 'd' not found in static DisjointSet. Use DisjointSet(..., dynamic=True) to allow adding elements via union."
        sage: S_static.make_set('d')
        Traceback (most recent call last):
        ...
        TypeError: Cannot add elements to a static DisjointSet (created with dynamic=False). Use dynamic=True if elements need to be added after creation.
        sage: isinstance(S_static, sage.sets.disjoint_set.DisjointSet_of_hashables)
        True

    Dynamic set from iterable::

        sage: S_dyn = DisjointSet('abc', dynamic=True)
        sage: S_dyn
        {{'a'}, {'b'}, {'c'}}
        sage: 'd' in S_dyn
        False
        sage: S_dyn.union('a', 'd')  # 'd' is added automatically
        sage: S_dyn
        {{'a', 'd'}, {'b'}, {'c'}}
        sage: 'd' in S_dyn
        True
        sage: S_dyn.union('e', 'f')  # 'e' and 'f' are added
        sage: S_dyn
        {{'a', 'd'}, {'b'}, {'c'}, {'e', 'f'}}
        sage: isinstance(S_dyn, sage.sets.disjoint_set.DynamicDisjointSet_of_hashables)
        True

    From SetPartition::

        sage: SP = SetPartition([[1,3], [2]])
        sage: D_static = DisjointSet(SP)  # Static by default
        sage: D_static
        {{1, 3}, {2}}
        sage: D_static.make_set(4)
        Traceback (most recent call last):
        ...
        TypeError: Cannot add elements to a static DisjointSet (created with dynamic=False). Use dynamic=True if elements need to be added after creation.
        sage: D_dyn = DisjointSet(SP, dynamic=True)  # Dynamic
        sage: D_dyn
        {{1, 3}, {2}}
        sage: D_dyn.make_set(4)
        sage: D_dyn
        {{1, 3}, {2}, {4}}


    TESTS::

        sage: DisjointSet(0)  # Integer case
        {}
        sage: DisjointSet('')  # Static hashable
        {}
        sage: DisjointSet([], dynamic=True)  # Dynamic hashable
        {}
        sage: DisjointSet(-1)
        Traceback (most recent call last):
        ...
        ValueError: arg must be a nonnegative integer (-1 given)
        sage: DisjointSet(4.3)  # needs sage.rings.real_mpfr
        Traceback (most recent call last):
        ...
        TypeError: Cannot create DisjointSet from 'RealLiteral'; input must be None, int, iterable, or SetPartition
        sage: DisjointSet([{}, {}])
        Traceback (most recent call last):
        ...
        TypeError: Input elements must be hashable: unhashable type: 'dict'

        # Test unpickling argument format
        sage: D_sp = DisjointSet( (['a','b'], False) )  # Static from pickle args
        sage: D_sp
        {{'a'}, {'b'}}
        sage: isinstance(D_sp, sage.sets.disjoint_set.DisjointSet_of_hashables)
        True
        sage: D_dp = DisjointSet( (['a','b'], True) )  # Dynamic from pickle args
        sage: D_dp
        {{'a'}, {'b'}}
        sage: isinstance(D_dp, sage.sets.disjoint_set.DynamicDisjointSet_of_hashables)
        True

    """
    # Handle None -> Empty Dynamic Set
    if arg is None:
        return DynamicDisjointSet_of_hashables([])

    # Handle Integer -> Static Integer Set
    if isinstance(arg, (Integer, int)):
        if arg < 0:
            raise ValueError('arg must be a nonnegative integer (%s given)' % arg)
        return DisjointSet_of_integers(arg)

    # Handle SetPartition
    if isinstance(arg, SetPartition):
        base = arg.base_set()
        elements = [] if base is None else list(base) 
        if dynamic:
            ds = DynamicDisjointSet_of_hashables(elements)
        else:
            ds = DisjointSet_of_hashables(elements)

        for block in arg:
            if len(block) > 1:
                block_list = list(block)
                first = block_list[0]
                for i in range(1, len(block_list)):
                    # Note: For dynamic=True, union handles adding if needed (though shouldn't be needed here)
                    # For dynamic=False, elements must exist (which they do)
                    ds.union(first, block_list[i])
        return ds
        
    # Handle unpickling args: (elements, dynamic_flag)
    if isinstance(arg, tuple) and len(arg) == 2 and isinstance(arg[1], bool):
        elements, dynamic_flag = arg
        # Ensure elements is iterable (should be list from getstate)
        try:
            iter(elements)
        except TypeError:
             raise TypeError("Invalid arguments for DisjointSet unpickling: elements part is not iterable") from None
        if dynamic_flag:
            return DynamicDisjointSet_of_hashables(elements)
        return DisjointSet_of_hashables(elements)

    # Handle general iterable -> Static or Dynamic Hashable Set
    try:
        it = iter(arg)
    except TypeError:
         raise TypeError(f"Cannot create DisjointSet from '{type(arg).__name__}'; input must be None, int, iterable, or SetPartition") from None

    if dynamic:
        return DynamicDisjointSet_of_hashables(arg)
    return DisjointSet_of_hashables(arg)

cdef class DisjointSet_class(SageObject):
    r"""
    Base class for DisjointSet implementations.
    """ 

    def __cinit__(self):
        self._nodes = NULL

    def __dealloc__(self):
        if self._nodes is not NULL:
            OP_dealloc(self._nodes)

    cpdef Py_ssize_t cardinality(self) noexcept:
        r"""
        Return the number of elements in the ground set.
        """
        return self._nodes.degree if self._nodes is not NULL else 0

    cpdef Py_ssize_t number_of_subsets(self) noexcept:
        r"""
        Return the number of disjoint subsets currently formed.
        """
        return self._nodes.num_cells if self._nodes is not NULL else 0

cdef class DisjointSet_of_integers(DisjointSet_class):
    r"""
    Static disjoint set of integers from ``0`` to ``n-1``.
    Elements cannot be added after creation.

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

    def __init__(self, Py_ssize_t n):
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
        if n < 0:
             raise ValueError("Argument must be a non-negative integer")
        self._nodes = OP_new(n)
        if self._nodes == NULL:
            raise MemoryError("Failed to allocate OrbitPartition structure")

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
        if self._nodes is NULL:
            return []

        cdef Py_ssize_t card = self._nodes.degree
        cdef list l = [None] * card
        cdef Py_ssize_t i
        for i in range(card):
            l[i] = self._nodes.parent[i]
        return l

    def __setstate__(self, l):
        r"""
        Set the internal parent array directly from the state list ``l``
        after validating the state and checking for cycles.

        This method will assume the object has already been initialized with the
        correct size via ``__reduce__`` and the factory and it directly
        manipulates the underlying C structure for efficiency during
        unpickling, but only after ensuring the provided state is valid
        and acyclic.

        INPUT:

        - ``l`` -- list of parent indices. ``l[i]`` is the parent of ``i``.
          The length of ``l`` must match ``self.cardinality()``. The implied
          parent graph must be a forest (acyclic).

        EXAMPLES::

            sage: d = DisjointSet(5)
            sage: state = [1, 1, 1, 3, 4]  # Represents unions (0,1), (2,1)
            sage: d.__setstate__(state)
            sage: d  # Represents {{0, 1, 2}, {3}, {4}} - check roots
            {{0, 1, 2}, {3}, {4}}
            sage: d.find(0), d.find(1), d.find(2), d.find(3), d.find(4)
            (1, 1, 1, 3, 4)

            sage: d = DisjointSet(5)
            sage: state = [1, 2, 3, 4, 4]  # Chain 0->1->2->3->4 (root)
            sage: d.__setstate__(state)
            sage: d  # Represents {{0, 1, 2, 3, 4}}
            {{0, 1, 2, 3, 4}}
            sage: len(list(d))  # Check number of subsets
            1

            sage: d = DisjointSet(5)
            sage: d.__setstate__([1, 1, 1])
            Traceback (most recent call last):
            ...
            ValueError: State list length 3 does not match cardinality 5

            sage: d = DisjointSet(5)
            sage: d.__setstate__([0, 1, 2, 5, 4])
            Traceback (most recent call last):
            ...
            ValueError: Invalid parent index 5 for element 3

            sage: d = DisjointSet(5)
            sage: state = [1, 2, 3, 4, 0]  # Cycle 0->1->2->3->4->0
            sage: d.__setstate__(state)
            Traceback (most recent call last):
            ...
            ValueError: Cycle detected in parent list starting at element 0

            sage: d = DisjointSet(4)
            sage: state = [1, 0, 3, 2]  # Cycles 0<->1 and 2<->3
            sage: d.__setstate__(state)
            Traceback (most recent call last):
            ...
            ValueError: Cycle detected in parent list starting at element 0

        TESTS::

            sage: d = DisjointSet(0)
            sage: d.__setstate__([])
            sage: d
            {}
        """
        cdef Py_ssize_t n = self.cardinality()
        if len(l) != n:
            raise ValueError(f"State list length {len(l)} does not match cardinality {n}")

        # Ensuring _nodes is allocated if n > 0
        if self._nodes is NULL and n > 0:
             # This case shouldn't ideally happen if i haven't done any blunder in pickling!
             raise RuntimeError("Internal error: _nodes not allocated during __setstate__")
        elif n == 0:
             return 

        # Validate all indices first
        cdef Py_ssize_t i, j
        cdef int parent_idx 
        for i in range(n):
            parent_idx = l[i]
            if not (0 <= parent_idx < n):
                raise ValueError(f"Invalid parent index {parent_idx} for element {i}")

        # Cycle detection with C array
        cdef int *visited_int = <int *> sig_calloc(n, sizeof(int))  # Use int array
        if not visited_int:
            raise MemoryError("Failed to allocate memory for cycle detection")
        cdef int visited_flag = 1  # Use incrementing flag instead of resetting
        try:
            for i in range(n):
                j = i
                while j != l[j]:
                    if visited_int[j] == visited_flag:
                        # Cycle detected!
                        raise ValueError(f"Cycle detected in parent list starting at element {i}")
                    if visited_int[j] != 0:
                        # Already visited in a previous traversal (part of a valid tree), stop early
                        break
                    visited_int[j] = visited_flag
                    j = l[j]
                if visited_int[j] == 0:  # Mark root if not visited before
                     visited_int[j] = visited_flag
                elif visited_int[j] != visited_flag:
                     # Root was visited by another traversal's flag, valid tree
                     pass
                # Else: visited_int[j] == visited_flag means we started at a root, which is fine.

                visited_flag += 1  # Increment flag for next starting node

            # If no cycles detected, apply the state
            for i in range(n):
                self._nodes.parent[i] = l[i]
        finally:
            # Ensuring allocated memory is always freed in the end
            sig_free(visited_int)


    def __reduce__(self):
        return DisjointSet, (self.cardinality(),), self.__getstate__()

    def __contains__(self, item):
        r"""
        Check if ``item`` is an integer within the valid range ``0`` to ``n-1``.

        EXAMPLES::

            sage: D = DisjointSet(5)
            sage: 3 in D
            True
            sage: 5 in D
            False
            sage: -1 in D
            False
            sage: 'a' in D
            False
        """
        return isinstance(item, (int, Integer)) and 0 <= item < self.cardinality()

    cpdef int find(self, int i) except? -1:
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
        cdef Py_ssize_t card = self._nodes.degree
        if not (0 <= i < card): 
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
        if self._nodes is NULL:
            raise ValueError("Disjoint set is empty")

        cdef Py_ssize_t card = self._nodes.degree
        if not (0 <= i < card):
            raise ValueError('i must be between 0 and %s (%s given)' % (card - 1, i))
        if not (0 <= j < card):
            raise ValueError('j must be between 0 and %s (%s given)' % (card - 1, j))
        OP_join(self._nodes, i, j)

    def make_set(self):
        r"""
        Raise an error: elements cannot be added to static integer sets."""
        raise TypeError("Cannot add elements to a static DisjointSet_of_integers.")

    # Representation and Iteration Methods
    def _repr_(self):
        res = []
        for l in (<dict?>self.root_to_elements_dict()).itervalues():
            l.sort()
            res.append('{%s}' % ', '.join(map(str, l)))
        res.sort()
        return '{%s}' % ', '.join(res)

    def __iter__(self):
        return iter((<dict?>self.root_to_elements_dict()).itervalues())

    def __richcmp__(self, other, int op):
        s = Set(map(frozenset, self.root_to_elements_dict().values()))
        try:
            if isinstance(other, DisjointSet_of_integers):
                 t = Set(map(frozenset, other.root_to_elements_dict().values()))
                 return PyObject_RichCompare(s, t, op)
            elif isinstance(other, DisjointSet_class):
                 t = Set(map(frozenset, other.root_to_elements_dict().values()))
                 # This comparison make sense in the case where elements are integers
                 return PyObject_RichCompare(s, t, op)
            else:
                 return NotImplemented
        except AttributeError: 
            return NotImplemented


    cpdef dict root_to_elements_dict(self):
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
        if self._nodes is NULL:
            return {} 

        cdef dict s = {}
        cdef Py_ssize_t i
        cdef int o
        for i in range(self._nodes.degree):
            o = OP_find(self._nodes, i)
            if o not in s:
                s[o] = []
            s[o].append(i)
        return s

    cpdef dict element_to_root_dict(self):
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
        if self._nodes is NULL:
            return {}

        cdef dict d = {}
        cdef Py_ssize_t i
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
        if self._nodes is NULL:
            from sage.graphs.digraph import DiGraph
            return DiGraph()

        cdef dict d = {i: [self._nodes.parent[i]] for i in range(self._nodes.degree)}
        from sage.graphs.digraph import DiGraph 
        return DiGraph(d)

cdef class DisjointSet_of_hashables(DisjointSet_class):
    r"""
    Static disjoint set of hashable objects.

    Elements must be provided at creation and cannot be added later.
    Use ``DisjointSet(..., dynamic=True)`` for dynamic sets.

    EXAMPLES::

        sage: # Static creation (default via factory)
        sage: d = DisjointSet('abcde')
        sage: d
        {{'a'}, {'b'}, {'c'}, {'d'}, {'e'}}
        sage: d.union('a', 'c')
        sage: d
        {{'a', 'c'}, {'b'}, {'d'}, {'e'}}
        sage: d.find('a')
        'a'
        sage: 'c' in d
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
        cdef Py_ssize_t i = 0
        self._int_to_el = []
        self._el_to_int = {}
        try:
            for e in iterable:
                hash(e)  # Checking hashability
                if e in self._el_to_int:
                     raise ValueError(f"Duplicate element found in input iterable: {e!r}")
                self._int_to_el.append(e)
                self._el_to_int[e] = i
                i += 1
        except TypeError as te:
             raise TypeError(f"Input elements must be hashable: {te}") from te

        self._nodes = OP_new(len(self._int_to_el))
        if self._nodes == NULL:
            raise MemoryError("Failed to allocate OrbitPartition structure")

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
             ([0, 1, 2, 3, 4], False),
             [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)])

         ::

            sage: d.union(2, 4)
            sage: d.union(1, 3)
            sage: d.__reduce__()
            (<built-in function DisjointSet>,
             ([0, 1, 2, 3, 4], False),
             [(0, 0), (1, 1), (2, 2), (3, 1), (4, 2)])
        """
        return DisjointSet, (self._int_to_el, False), self.__getstate__()

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
        if self._nodes is NULL:
            return []

        cdef Py_ssize_t card = self._nodes.degree
        cdef list parent_elements = [None] * card
        cdef Py_ssize_t i, parent_idx
        for i in range(card):
            parent_idx = self._nodes.parent[i]
            if 0 <= parent_idx < len(self._int_to_el):
                 parent_elements[i] = self._int_to_el[parent_idx]
            else:
                 raise IndexError(f"Internal state error: Invalid parent index {parent_idx} for element index {i}")
        return list(zip(self._int_to_el, parent_elements))

    def __setstate__(self, state_list):
        r"""
        Sets the internal parent array directly from the state list
        after validating the state and checking for cycles.

        This method assumes the object has already been initialized with the
        correct elements via ``__reduce__`` and the factory and it directly
        manipulates the underlying C structure for better efficiency during
        unpickling, but only after ensuring the provided state is valid
        and acyclic.

        INPUT:

        - ``state_list`` -- List of ``(element, parent_element)`` pairs.
          The length and order of elements must match the object's internal
          element list ``_int_to_el``. The implied parent graph must be a
          forest (acyclic).

        EXAMPLES::

            sage: d = DisjointSet('abcde')
            sage: state = [('a', 'a'), ('b', 'b'), ('c', 'c'), ('d', 'd'), ('e', 'e')]
            sage: d.__setstate__(state)
            sage: d
            {{'a'}, {'b'}, {'c'}, {'d'}, {'e'}}

            sage: d = DisjointSet('abcde')
            sage: state = [('a', 'b'), ('b', 'b'), ('c', 'd'), ('d', 'd'), ('e', 'e')]
            sage: d.__setstate__(state)
            sage: d  # Represents {{'a', 'b'}, {'c', 'd'}, {'e'}}
            {{'a', 'b'}, {'c', 'd'}, {'e'}}
            sage: d.find('a')
            'b'
            sage: d.find('c')
            'd'

            sage: d = DisjointSet('abcde')
            sage: state = [('a', 'b'), ('b', 'c'), ('c', 'd'), ('d', 'e'), ('e', 'e')]
            sage: d.__setstate__(state)
            sage: d  # Represents {{'a', 'b', 'c', 'd', 'e'}}
            {{'a', 'b', 'c', 'd', 'e'}}
            sage: len(list(d))
            1

            sage: d = DisjointSet('abcde')
            sage: d.__setstate__([('a', 'a'), ('b', 'b'), ('c', 'c')])
            Traceback (most recent call last):
            ...
            ValueError: State list length 3 does not match number of elements 5

            sage: d = DisjointSet('abc')
            sage: d.__setstate__([('a', 'a'), ('x', 'b'), ('c', 'c')])
            Traceback (most recent call last):
            ...
            ValueError: State elements do not match at index 1: expected 'b', got 'x'

            sage: d = DisjointSet('abc')
            sage: d.__setstate__([('a', 'a'), ('b', 'x'), ('c', 'c')])
            Traceback (most recent call last):
            ...
            ValueError: Invalid parent element 'x' not found in set

            sage: d = DisjointSet('abc')
            sage: state = [('a', 'b'), ('b', 'c'), ('c', 'a')]  # Cycle a->b->c->a
            sage: d.__setstate__(state)
            Traceback (most recent call last):
            ...
            ValueError: Cycle detected in parent list at element 'a'

        TESTS::
            sage: d = DisjointSet([])
            sage: d.__setstate__([])
            sage: d
            {}
        """
        cdef Py_ssize_t n = len(self._int_to_el)
        if len(state_list) != n:
            raise ValueError(f"State list length {len(state_list)} does not match number of elements {n}")

        if self._nodes is NULL and n > 0:
             raise RuntimeError("Internal error: _nodes not allocated during __setstate__")
        elif n == 0: 
             return 

        # First Pass: Validate state and build temporary parent index array 
        cdef list temp_parent = [0] * n 
        cdef Py_ssize_t i, j, p_i 
        cdef object e, parent_e, expected_e

        for i, (e, parent_e) in enumerate(state_list):
            expected_e = self._int_to_el[i]
            if e != expected_e:
                raise ValueError(f"State elements do not match at index {i}: expected {expected_e!r}, got {e!r}")

            # Checking if parent element exists and get its index
            try:
                p_i = self._el_to_int[parent_e]
            except KeyError:
                raise ValueError(f"Invalid parent element {parent_e!r} not found in set") from None

            # Store parent index in temporary list
            temp_parent[i] = p_i

        # Second Pass: Check for cycles using the temporary parent index array
        cdef int *visited_int = <int *> sig_calloc(n, sizeof(int))  # Use int array
        if not visited_int:
            raise MemoryError("Failed to allocate memory for cycle detection")
        cdef int visited_flag = 1  # Use incrementing flag
        try:
            for i in range(n):
                j = i
                while j != temp_parent[j]:
                    if visited_int[j] == visited_flag:
                        # Cycle detected!
                        raise ValueError(f"Cycle detected in parent list at element {self._int_to_el[i]!r}")
                    if visited_int[j] != 0:
                        break
                    visited_int[j] = visited_flag
                    j = temp_parent[j]
                # Check root if loop terminated normally (and mark if first time seen)
                if visited_int[j] == 0:  # Mark root if not visited before
                     visited_int[j] = visited_flag
                elif visited_int[j] != visited_flag:
                     pass

                visited_flag += 1  # Increment flag for next starting node

            # If no cycles detected, apply the state
            for i in range(n):
                self._nodes.parent[i] = temp_parent[i]
        finally:
            # Ensuring allocated memory is always freed in the end
            sig_free(visited_int)



    def __contains__(self, item):
        r"""
        Check if ``item`` is one of the elements in this disjoint set.

        EXAMPLES::

            sage: D = DisjointSet('abc')
            sage: 'a' in D
            True
            sage: 'd' in D
            False
            sage: D_int = DisjointSet([1, 2])
            sage: 1 in D_int
            True
            sage: 3 in D_int
            False
        """
        return item in self._el_to_int

    cpdef find(self, e):
        r"""
        Return the representative of the set containing element ``e``.
        Raises KeyError if ``e`` is not in the set.

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
        cdef int i, r
        try:
            i = <int> self._el_to_int[e]
        except KeyError:
            raise KeyError(e)

        if self._nodes is NULL:
            raise RuntimeError("Internal state error: _nodes is NULL but element map is not empty")

        r = <int> OP_find(self._nodes, i)
        if 0 <= r < len(self._int_to_el):
             return self._int_to_el[r]
        else:
             raise IndexError(f"Internal state error: Invalid representative index {r}")

    cpdef void union(self, e, f) except *:
        r"""
        Combine the sets containing ``e`` and ``f``.
        Raises KeyError if ``e`` or ``f`` are not present in the set.

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
            Traceback (most recent call last):
            ...
            KeyError: 'Element 1024 not found in static DisjointSet. Use DisjointSet(..., dynamic=True) to allow adding elements via union.'
        """
        cdef int i, j
        try:
            i = <int> self._el_to_int[e]
        except KeyError:
            raise KeyError(f"Element {e!r} not found in static DisjointSet. Use DisjointSet(..., dynamic=True) to allow adding elements via union.") from None
        try:
            j = <int> self._el_to_int[f]
        except KeyError:
            raise KeyError(f"Element {f!r} not found in static DisjointSet. Use DisjointSet(..., dynamic=True) to allow adding elements via union.") from None
        
        if self._nodes is NULL:
            raise RuntimeError("Internal state error: _nodes is NULL but element map is not empty")
        OP_join(self._nodes, i, j)

    def make_set(self, new_elt=None):
        r"""
        Raise an error: elements cannot be added to static hashable sets.
        """
        raise TypeError("Cannot add elements to a static DisjointSet (created with dynamic=False). Use dynamic=True if elements need to be added after creation.")

    def _repr_(self):
        res = []
        for l in (<dict?>self.root_to_elements_dict()).itervalues():
            try:
                 l_sorted = sorted(l)
            except TypeError:
                 l_sorted = l 
            res.append('{%s}' % ', '.join(map(repr, l_sorted)))
        res.sort()
        return '{%s}' % ', '.join(res)

    def __iter__(self):
        return iter((<dict?>self.root_to_elements_dict()).itervalues())

    def __richcmp__(self, other, int op):
        s = Set(map(frozenset, self.root_to_elements_dict().values()))
        try:
            if isinstance(other, (DisjointSet_of_hashables, DynamicDisjointSet_of_hashables, DisjointSet_of_integers)):
                 t = Set(map(frozenset, other.root_to_elements_dict().values()))
                 return PyObject_RichCompare(s, t, op)
            else:
                 return NotImplemented
        except AttributeError:
            return NotImplemented

    cpdef dict root_to_elements_dict(self):
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
        if self._nodes is NULL and self._int_to_el:
            raise RuntimeError("Internal state error: _nodes is NULL but element list is not empty")
        cdef dict s = {}
        cdef object r, e
        for e in self._int_to_el:
            r = self.find(e) 
            if r not in s:
                s[r] = []
            s[r].append(e)
        return s

    cpdef dict element_to_root_dict(self):
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
        if self._nodes is NULL and self._int_to_el:
            raise RuntimeError("Internal state error: _nodes is NULL but element list is not empty")
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
        if self._nodes is NULL:
            from sage.graphs.digraph import DiGraph
            return DiGraph()

        cdef dict d = {}
        cdef Py_ssize_t i, p_idx
        cdef object e, p
        from sage.graphs.digraph import DiGraph

        for i in range(self._nodes.degree):
            if i < len(self._int_to_el):
                 e = self._int_to_el[i]
                 p_idx = self._nodes.parent[i]
                 if 0 <= p_idx < len(self._int_to_el):
                     p = self._int_to_el[p_idx]
                     d[e] = [p]
                 else:
                     raise IndexError(f"Internal state error: Invalid parent index {p_idx} for element index {i}")
            else:
                 raise IndexError(f"Internal state error: Invalid element index {i}")
        return DiGraph(d)

cdef class DynamicDisjointSet_of_hashables(DisjointSet_of_hashables):  
    r"""
    Dynamic disjoint set of hashable objects.

    Elements can be provided at creation or added implicitly via `union`
    or explicitly via `make_set`.
    """

    def __init__(self, iterable=None):
        r"""
        Construct the dynamic disjoint set, optionally from an iterable.

        INPUT:

        - ``iterable`` -- (default: ``None``) An iterable of hashable
          elements, or ``None`` to create an empty dynamic set.

        EXAMPLES::

            sage: D_empty = DisjointSet(dynamic=True)
            sage: D_empty
            {}
            sage: isinstance(D_empty, sage.sets.disjoint_set.DynamicDisjointSet_of_hashables)
            True

            sage: D_abc = DisjointSet(['a', 'b', 'c'], dynamic=True)
            sage: D_abc
            {{'a'}, {'b'}, {'c'}}
            sage: isinstance(D_abc, sage.sets.disjoint_set.DynamicDisjointSet_of_hashables)
            True

            sage: D = DisjointSet(['a', 'b'], dynamic=True)
            sage: D.union('a', 'b')
            sage: factory, args, state = D.__reduce__()
            sage: factory is DisjointSet
            True
            sage: args  # Note the True flag
            (['a', 'b'], True)
            sage: isinstance(state, list) and len(state) == 2
            True
            sage: isinstance(state[0], tuple) and len(state[0]) == 2
            True

        TESTS::

            sage: DisjointSet(['a', 'b', 'a'], dynamic=True)
            Traceback (most recent call last):
            ...
            ValueError: Duplicate element found in input iterable: 'a'

            sage: DisjointSet([['a']], dynamic=True)
            Traceback (most recent call last):
            ...
            TypeError: Input elements must be hashable: unhashable type: 'list'

        """
        self._int_to_el = []
        self._el_to_int = {}
        cdef Py_ssize_t i = 0

        elements_to_process = [] if iterable is None else iterable

        try:
            for e in elements_to_process:
                hash(e) 
                if e in self._el_to_int:
                     raise ValueError(f"Duplicate element found in input iterable: {e!r}")
                self._int_to_el.append(e)
                self._el_to_int[e] = i
                i += 1
        except TypeError as te:
             raise TypeError(f"Input elements must be hashable: {te}") from te

        self._nodes = OP_new(len(self._int_to_el))
        if self._nodes == NULL:
            raise MemoryError("Failed to allocate OrbitPartition structure")

    # Override __reduce__ to set dynamic_flag=True
    def __reduce__(self):
        return DisjointSet, (self._int_to_el, True), self.__getstate__()

    # Override union to add elements
    cpdef void union(self, e, f) except *:
        r"""
        Combine the sets containing ``e`` and ``f``.

        If ``e`` or ``f`` are not already present in the set, they are
        automatically added via ``make_set`` before performing the union.
        """
        cdef int i, j
        if e not in self:
            self.make_set(e)
        if f not in self:
            self.make_set(f)

        i = <int> self._el_to_int[e]
        j = <int> self._el_to_int[f]
        OP_join(self._nodes, i, j)

    # Override make_set to provide the working implementation
    def make_set(self, new_elt=None):
        r"""
        Add ``new_elt`` to the disjoint set as a new singleton subset.

        If ``new_elt`` is already present, this operation does nothing.
        If ``new_elt`` is ``None``, a unique integer element is generated
        and added.

        INPUT:

        - ``new_elt`` -- (optional) A hashable element to add. If `None`,
          an integer is generated.

        EXAMPLES::

            sage: D = DisjointSet(dynamic=True)
            sage: D.make_set('a')
            sage: D
            {{'a'}}
            sage: D.make_set('b')
            sage: D
            {{'a'}, {'b'}}
            sage: D.make_set('a')  # Already present, no change
            sage: D
            {{'a'}, {'b'}}
            sage: D.make_set()  # Add generated integer
            sage: D
            {{'a'}, {'b'}, {2}}
            sage: D.make_set()
            sage: D
            {{'a'}, {'b'}, {2}, {3}}
        """
        cdef Py_ssize_t new_index, current_degree, i
        cdef OrbitPartition *new_nodes

        if new_elt is None:
            current_card = self.cardinality()
            new_elt = current_card
            while new_elt in self._el_to_int:
                 new_elt += 1
        else:
             # Checking hashability
             try:
                 hash(new_elt)
             except TypeError as te:
                 raise TypeError(f"New element must be hashable: {te}") from te

        if new_elt not in self._el_to_int:
            self._int_to_el.append(new_elt)
            new_index = len(self._int_to_el) - 1  # NEW
            self._el_to_int[new_elt] = new_index

            current_degree = 0 if self._nodes is NULL else self._nodes.degree
            new_nodes = OP_new(current_degree + 1)
            if new_nodes == NULL:
                # Cleans up mappings if allocation fails
                self._el_to_int.pop(new_elt)
                self._int_to_el.pop()
                raise MemoryError("Failed to allocate OrbitPartition structure in make_set")

            if self._nodes is not NULL:
                for i in range(current_degree):
                    new_nodes.parent[i] = self._nodes.parent[i]
                    new_nodes.rank[i] = self._nodes.rank[i]
                # Increments cell count for the new singleton set
                new_nodes.num_cells = self._nodes.num_cells + 1
                # Deallocates old nodes after successful copy and setup
                OP_dealloc(self._nodes)
            else:
                new_nodes.num_cells = 1

            # Initialize the new element as a singleton set in the new structure
            new_nodes.parent[new_index] = new_index
            new_nodes.rank[new_index] = 0

            # Update self._nodes to point to the new structure
            self._nodes = new_nodes
