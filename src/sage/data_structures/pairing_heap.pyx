# distutils: language = c++
r"""
Pairing Heap

This module proposes several implementations of the pairing heap data structure
[FSST1986]_. See the :wikipedia:`Pairing_heap` for more information on this
min-heap data structure.

- :class:`PairingHeap_of_n_integers`: a pairing heap data structure with fixed
  capacity `n`. Its items are integers in the range `0..n-1`. Values can be of
  any type equipped with a comparison method (``<=``).

- :class:`PairingHeap_of_n_hashables`: a pairing heap data structure with fixed
  capacity `n`. Its items can be of any hashable type. Values can be of any type
  equipped with a comparison method (``<=``).

- ``PairingHeap``: interface to a pairing heap data structure written in C++.
  The advantages of this data structure are that: its capacity is unbounded;
  items can be of any hashable type equipped with a hashing method that can be
  supported by ``std::unordered_map``; values can be of any specified type
  equipped with a comparison method (``<=``). This data structure is for
  internal use and therefore cannot be accessed from a shell.

EXAMPLES:

Pairing heap of `n` integers in the range `0..n-1`::

    sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
    sage: P = PairingHeap_of_n_integers(10); P
    PairingHeap_of_n_integers: capacity 10, size 0
    sage: P.push(1, 3)
    sage: P.push(2, 2)
    sage: P
    PairingHeap_of_n_integers: capacity 10, size 2
    sage: P.top()
    (2, 2)
    sage: P.decrease(1, 1)
    sage: P.top()
    (1, 1)
    sage: P.pop()
    sage: P.top()
    (2, 2)

    sage: P = PairingHeap_of_n_integers(10)
    sage: P.push(1, (2, 'a'))
    sage: P.push(2, (2, 'b'))
    sage: P.top()
    (1, (2, 'a'))

Pairing heap of `n` hashables::

    sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_hashables
    sage: P = PairingHeap_of_n_hashables(10); P
    PairingHeap_of_n_hashables: capacity 10, size 0
    sage: P.push(1, 3)
    sage: P.push('b', 2)
    sage: P.push((1, 'abc'), 4)
    sage: P.top()
    ('b', 2)
    sage: P.decrease((1, 'abc'), 1)
    sage: P.top()
    ((1, 'abc'), 1)
    sage: P.pop()
    sage: P.top()
    ('b', 2)

    sage: # needs sage.graphs
    sage: P = PairingHeap_of_n_hashables(10)
    sage: P.push(('a', 1), (2, 'b'))
    sage: P.push(2, (2, 'a'))
    sage: g = Graph(2, immutable=True)
    sage: P.push(g, (3, 'z'))
    sage: P.top()
    (2, (2, 'a'))
    sage: P.decrease(g, (1, 'z'))
    sage: P.top()
    (Graph on 2 vertices, (1, 'z'))
    sage: while P:
    ....:     print(P.top())
    ....:     P.pop()
    (Graph on 2 vertices, (1, 'z'))
    (2, (2, 'a'))
    (('a', 1), (2, 'b'))

AUTHORS:

- David Coudert (2024) - Initial version.

"""
# ******************************************************************************
# Copyright (C) 2024 David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ******************************************************************************


from libcpp.pair cimport pair
from cpython.ref cimport PyObject, Py_INCREF, Py_XDECREF
from cysignals.signals cimport sig_on, sig_off, sig_check
from cysignals.memory cimport check_allocarray, sig_free
from sage.data_structures.bitset_base cimport (bitset_init, bitset_free,
                                               bitset_add, bitset_remove,
                                               bitset_in)
from sage.misc.prandom import shuffle


# ==============================================================================
# Class PairingHeap_class
# ==============================================================================

cdef class PairingHeap_class:
    r"""
    Common class and methods for :class:`PairingHeap_of_n_integers` and
    :class:`PairingHeap_of_n_hashables`.
    """
    def __repr__(self):
        r"""
        Return a string representing ``self``.

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
            sage: P = PairingHeap_of_n_integers(5); P
            PairingHeap_of_n_integers: capacity 5, size 0
            sage: P.push(1, 2)
            sage: P
            PairingHeap_of_n_integers: capacity 5, size 1
        """
        return f"{type(self).__name__}: capacity {self.n}, size {len(self)}"

    def __bool__(self):
        r"""
        Check whether ``self`` is not empty.

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
            sage: P = PairingHeap_of_n_integers(5)
            sage: 'not empty' if P else 'empty'
            'empty'
            sage: P.push(1, 2)
            sage: 'not empty' if P else 'empty'
            'not empty'
        """
        return self.root != NULL

    cpdef bint empty(self) noexcept:
        r"""
        Check whether the heap is empty or not.

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
            sage: P = PairingHeap_of_n_integers(5)
            sage: P.empty()
            True
            sage: P.push(1, 2)
            sage: P.empty()
            False
        """
        return self.root == NULL

    cpdef bint full(self) noexcept:
        r"""
        Check whether the heap is full or not.

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
            sage: P = PairingHeap_of_n_integers(2)
            sage: P.full()
            False
            sage: P.push(0, 2)
            sage: P.push(1, 3)
            sage: P.full()
            True
        """
        return self.n == self.number_of_items

    def __len__(self):
        r"""
        Return the number of items in the heap.

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
            sage: P = PairingHeap_of_n_integers(5)
            sage: len(P)
            0
            sage: P.push(1, 2)
            sage: len(P)
            1
        """
        return self.number_of_items

    size = __len__

    cpdef tuple top(self):
        r"""
        Return the top pair (item, value) of the heap.

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
            sage: P = PairingHeap_of_n_integers(5)
            sage: P.push(1, 2)
            sage: P.top()
            (1, 2)
            sage: P.push(3, 1)
            sage: P.top()
            (3, 1)

            sage: P = PairingHeap_of_n_integers(3)
            sage: P.top()
            Traceback (most recent call last):
            ...
            ValueError: trying to access the top of an empty heap
        """
        raise NotImplementedError()

    cpdef object top_value(self):
        r"""
        Return the value of the top item of the heap.

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
            sage: P = PairingHeap_of_n_integers(5)
            sage: P.push(1, 2)
            sage: P.top()
            (1, 2)
            sage: P.top_value()
            2

            sage: P = PairingHeap_of_n_integers(3)
            sage: P.top_value()
            Traceback (most recent call last):
            ...
            ValueError: trying to access the top of an empty heap
        """
        if not self:
            raise ValueError("trying to access the top of an empty heap")
        return <object>self.root.value

    cpdef void pop(self) noexcept:
        r"""
        Remove the top item from the heap.

        If the heap is already empty, we do nothing.

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
            sage: P = PairingHeap_of_n_integers(5); P
            PairingHeap_of_n_integers: capacity 5, size 0
            sage: P.push(1, 2); P
            PairingHeap_of_n_integers: capacity 5, size 1
            sage: P.pop(); P
            PairingHeap_of_n_integers: capacity 5, size 0
            sage: P.pop(); P
            PairingHeap_of_n_integers: capacity 5, size 0
        """
        raise NotImplementedError()

# ==============================================================================
# Class PairingHeap_of_n_integers
# ==============================================================================

cdef class PairingHeap_of_n_integers(PairingHeap_class):
    r"""
    Pairing Heap for items in range [0..n - 1].

    EXAMPLES::

        sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
        sage: P = PairingHeap_of_n_integers(5); P
        PairingHeap_of_n_integers: capacity 5, size 0
        sage: P.push(1, 3)
        sage: P.push(2, 2)
        sage: P
        PairingHeap_of_n_integers: capacity 5, size 2
        sage: P.top()
        (2, 2)
        sage: P.decrease(1, 1)
        sage: P.top()
        (1, 1)
        sage: P.pop()
        sage: P.top()
        (2, 2)

    TESTS::

        sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
        sage: P = PairingHeap_of_n_integers(0)
        Traceback (most recent call last):
        ...
        ValueError: the capacity of the heap must be strictly positive
        sage: P = PairingHeap_of_n_integers(1); P
        PairingHeap_of_n_integers: capacity 1, size 0
        sage: P = PairingHeap_of_n_integers(5)
        sage: P.push(11, 3)
        Traceback (most recent call last):
        ...
        ValueError: item must be in range 0..4
    """
    def __init__(self, size_t n):
        r"""
        Construct the ``PairingHeap_of_n_integers`` where items are integers
        from ``0`` to ``n-1``.

        INPUT:

        - ``n`` -- strictly positive integer; the maximum number of items in the heap

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
            sage: P = PairingHeap_of_n_integers(5); P
            PairingHeap_of_n_integers: capacity 5, size 0
            sage: P.push(1, 2)
            sage: P
            PairingHeap_of_n_integers: capacity 5, size 1
            sage: P.push(2, 3)
            sage: P
            PairingHeap_of_n_integers: capacity 5, size 2
            sage: P.pop()
            sage: P
            PairingHeap_of_n_integers: capacity 5, size 1
            sage: P.push(10, 1)
            Traceback (most recent call last):
            ...
            ValueError: item must be in range 0..4
            sage: P = PairingHeap_of_n_integers(0)
            Traceback (most recent call last):
            ...
            ValueError: the capacity of the heap must be strictly positive
            sage: P = PairingHeap_of_n_integers(1); P
            PairingHeap_of_n_integers: capacity 1, size 0
        """
        if not n:
            raise ValueError("the capacity of the heap must be strictly positive")
        self.n = n
        self.root = NULL
        self.nodes = <PairingHeapNodePy *>check_allocarray(n, sizeof(PairingHeapNodePy))
        bitset_init(self.active, n)
        self.number_of_items = 0

    def __dealloc__(self):
        """
        Deallocate ``self``.

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
            sage: P = PairingHeap_of_n_integers(5)
            sage: del P
        """
        sig_free(self.nodes)
        bitset_free(self.active)

    cpdef void push(self, size_t item, object value) except *:
        r"""
        Insert an item into the heap with specified value (priority).

        INPUT:

        - ``item`` -- non negative integer; the item to consider

        - ``value`` -- the value associated with ``item``

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
            sage: P = PairingHeap_of_n_integers(5)
            sage: P.push(1, 2)
            sage: P.top()
            (1, 2)
            sage: P.push(3, 1)
            sage: P.top()
            (3, 1)

        TESTS::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
            sage: P = PairingHeap_of_n_integers(5)
            sage: P.push(1, 2)
            sage: P.push(1, 2)
            Traceback (most recent call last):
            ...
            ValueError: 1 is already in the heap
            sage: P.push(11, 2)
            Traceback (most recent call last):
            ...
            ValueError: item must be in range 0..4
        """
        if item >= self.n:
            raise ValueError(f"item must be in range 0..{self.n - 1}")
        if item in self:
            raise ValueError(f"{item} is already in the heap")

        cdef PairingHeapNodePy * p = self.nodes + item
        Py_INCREF(value)
        p.value = <PyObject *>value
        p.prev = p.next = p.child = NULL
        if self.root == NULL:
            self.root = p
        else:
            self.root = PairingHeapNodePy._merge(self.root, p)
        bitset_add(self.active, item)
        self.number_of_items += 1

    cpdef tuple top(self):
        r"""
        Return the top pair (item, value) of the heap.

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
            sage: P = PairingHeap_of_n_integers(5)
            sage: P.push(1, 2)
            sage: P.top()
            (1, 2)
            sage: P.push(3, 1)
            sage: P.top()
            (3, 1)

            sage: P = PairingHeap_of_n_integers(3)
            sage: P.top()
            Traceback (most recent call last):
            ...
            ValueError: trying to access the top of an empty heap
        """
        if not self:
            raise ValueError("trying to access the top of an empty heap")
        return self.root - self.nodes, <object>self.root.value

    cpdef size_t top_item(self) except *:
        r"""
        Return the top item of the heap.

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
            sage: P = PairingHeap_of_n_integers(5)
            sage: P.push(1, 2)
            sage: P.top()
            (1, 2)
            sage: P.top_item()
            1

            sage: P = PairingHeap_of_n_integers(3)
            sage: P.top_item()
            Traceback (most recent call last):
            ...
            ValueError: trying to access the top of an empty heap
        """
        if not self:
            raise ValueError("trying to access the top of an empty heap")
        return self.root - self.nodes

    cpdef void pop(self) noexcept:
        r"""
        Remove the top item from the heap.

        If the heap is already empty, we do nothing.

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
            sage: P = PairingHeap_of_n_integers(5); P
            PairingHeap_of_n_integers: capacity 5, size 0
            sage: P.push(1, 2); P
            PairingHeap_of_n_integers: capacity 5, size 1
            sage: P.push(2, 3); P
            PairingHeap_of_n_integers: capacity 5, size 2
            sage: P.pop(); P
            PairingHeap_of_n_integers: capacity 5, size 1
            sage: P.pop(); P
            PairingHeap_of_n_integers: capacity 5, size 0
            sage: P.pop(); P
            PairingHeap_of_n_integers: capacity 5, size 0
        """
        if not self:
            return
        cdef size_t item = self.top_item()
        Py_XDECREF(<PyObject *>self.nodes[item].value)
        bitset_remove(self.active, item)
        self.number_of_items -= 1
        self.root = PairingHeapNodePy._pair(self.root.child)

    cpdef void decrease(self, size_t item, object new_value) except *:
        r"""
        Decrease the value of specified item.

        This method is more permissive than it should as it can also be used to
        push an item in the heap.

        INPUT:

        - ``item`` -- non negative integer; the item to consider

        - ``new_value`` -- the new value for ``item``

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
            sage: P = PairingHeap_of_n_integers(5)
            sage: 3 in P
            False
            sage: P.decrease(3, 33)
            sage: 3 in P
            True
            sage: P.top()
            (3, 33)
            sage: P.push(1, 10)
            sage: P.top()
            (1, 10)
            sage: P.decrease(3, 7)
            sage: P.top()
            (3, 7)

        TESTS::

            sage: P = PairingHeap_of_n_integers(5)
            sage: P.push(1, 3)
            sage: P.decrease(1, 2)
            sage: P.decrease(1, 2)
            Traceback (most recent call last):
            ...
            ValueError: the new value must be less than the current value
        """
        if item >= self.n:
            raise ValueError(f"item must be in range 0..{self.n - 1}")
        cdef PairingHeapNodePy * p
        if bitset_in(self.active, item):
            p = self.nodes + item
            if <object>p.value <= new_value:
                raise ValueError("the new value must be less than the current value")
            Py_XDECREF(<PyObject *>p.value)
            Py_INCREF(new_value)
            p.value = <PyObject *>new_value
            if p.prev != NULL:
                PairingHeapNodePy._unlink(p)
                self.root = PairingHeapNodePy._merge(self.root, p)
        else:
            self.push(item, new_value)

    def __contains__(self, size_t item):
        r"""
        Check whether the specified item is in the heap.

        INPUT:

        - ``item`` -- non negative integer; the item to consider

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
            sage: P = PairingHeap_of_n_integers(5)
            sage: 3 in P
            False
            sage: P.push(3, 33)
            sage: 3 in P
            True
            sage: 100 in P
            False
        """
        if item >= self.n:
            return False
        return bitset_in(self.active, item)

    contains = __contains__

    cpdef object value(self, size_t item):
        r"""
        Return the value associated with the item.

        INPUT:

        - ``item`` -- non negative integer; the item to consider

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
            sage: P = PairingHeap_of_n_integers(5)
            sage: P.push(3, 33)
            sage: P.push(1, 10)
            sage: P.value(3)
            33
            sage: P.value(7)
            Traceback (most recent call last):
            ...
            ValueError: 7 is not in the heap
        """
        if item not in self:
            raise ValueError(f"{item} is not in the heap")
        return <object>self.nodes[item].value


# ==============================================================================
# Class PairingHeap_of_n_hashables
# ==============================================================================

cdef class PairingHeap_of_n_hashables(PairingHeap_class):
    r"""
    Pairing Heap for ``n`` hashable items.

    EXAMPLES::

        sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_hashables
        sage: P = PairingHeap_of_n_hashables(5); P
        PairingHeap_of_n_hashables: capacity 5, size 0
        sage: P.push(1, 3)
        sage: P.push('abc', 2)
        sage: P
        PairingHeap_of_n_hashables: capacity 5, size 2
        sage: P.top()
        ('abc', 2)
        sage: P.decrease(1, 1)
        sage: P.top()
        (1, 1)
        sage: P.pop()
        sage: P.top()
        ('abc', 2)

        sage: P = PairingHeap_of_n_hashables(5)
        sage: P.push(1, (2, 3))
        sage: P.push('a', (2, 2))
        sage: P.push('b', (3, 3))
        sage: P.push('c', (2, 1))
        sage: P.top()
        ('c', (2, 1))
        sage: P.push(Graph(2, immutable=True), (1, 7))
        sage: P.top()
        (Graph on 2 vertices, (1, 7))
        sage: P.decrease('b', (1, 5))
        sage: P.top()
        ('b', (1, 5))

    TESTS::

        sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_hashables
        sage: P = PairingHeap_of_n_hashables(0)
        Traceback (most recent call last):
        ...
        ValueError: the capacity of the heap must be strictly positive
        sage: P = PairingHeap_of_n_hashables(1); P
        PairingHeap_of_n_hashables: capacity 1, size 0
        sage: P.push(11, 3)
        sage: P.push(12, 4)
        Traceback (most recent call last):
        ...
        ValueError: the heap is full

        sage: P = PairingHeap_of_n_hashables(10)
        sage: P.push(1, 'John')
        sage: P.push(4, 42)
        Traceback (most recent call last):
        ...
        TypeError: unsupported operand parent(s) for >=: 'Integer Ring' and '<class 'str'>'
    """
    def __init__(self, size_t n):
        r"""
        Construct the ``PairingHeap_of_n_hashables``.

        This pairing heap has a maximum capacity of `n` items and each item is a
        hashable object.

        INPUT:

        - ``n`` -- strictly positive integer; the maximum number of items in the heap

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_hashables
            sage: P = PairingHeap_of_n_hashables(2); P
            PairingHeap_of_n_hashables: capacity 2, size 0
            sage: P.push(1, 2)
            sage: P
            PairingHeap_of_n_hashables: capacity 2, size 1
            sage: P.push(2, 3)
            sage: P
            PairingHeap_of_n_hashables: capacity 2, size 2
            sage: P.full()
            True
            sage: P.push(10, 1)
            Traceback (most recent call last):
            ...
            ValueError: the heap is full
            sage: P.pop()
            sage: P
            PairingHeap_of_n_hashables: capacity 2, size 1
            sage: P.push(10, 1)

    TESTS::

            sage: P = PairingHeap_of_n_hashables(0)
            Traceback (most recent call last):
            ...
            ValueError: the capacity of the heap must be strictly positive
            sage: P = PairingHeap_of_n_hashables(1); P
            PairingHeap_of_n_hashables: capacity 1, size 0
        """
        if not n:
            raise ValueError("the capacity of the heap must be strictly positive")
        self.n = n
        self.root = NULL
        self.nodes = <PairingHeapNodePy *>check_allocarray(n, sizeof(PairingHeapNodePy))
        self.number_of_items = 0
        self._int_to_item = [None] * n
        self._item_to_int = dict()
        self.free_idx = list(range(n))

    def __dealloc__(self):
        """
        Deallocate ``self``.

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_integers
            sage: P = PairingHeap_of_n_integers(5)
            sage: del P
        """
        sig_free(self.nodes)

    cpdef void push(self, object item, object value) except *:
        r"""
        Insert an item into the heap with specified value (priority).

        INPUT:

        - ``item`` -- non negative integer; the item to consider

        - ``value`` -- the value associated with ``item``

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_hashables
            sage: P = PairingHeap_of_n_hashables(5)
            sage: P.push(1, 2)
            sage: P.top()
            (1, 2)
            sage: P.push(3, 1)
            sage: P.top()
            (3, 1)

        TESTS::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_hashables
            sage: P = PairingHeap_of_n_hashables(2)
            sage: P.push(1, 2)
            sage: P.push(1, 2)
            Traceback (most recent call last):
            ...
            ValueError: 1 is already in the heap
            sage: P.push(11, 2)
            sage: P.push(7, 5)
            Traceback (most recent call last):
            ...
            ValueError: the heap is full
        """
        if item in self:
            raise ValueError(f"{item} is already in the heap")
        if self.full():
            raise ValueError("the heap is full")

        cdef size_t idx = self.free_idx.pop()
        self._int_to_item[idx] = item
        self._item_to_int[item] = idx
        cdef PairingHeapNodePy * p = self.nodes + idx
        Py_INCREF(value)
        p.value = <PyObject *>value
        p.prev = p.next = p.child = NULL
        if self.root == NULL:
            self.root = p
        else:
            self.root = PairingHeapNodePy._merge(self.root, p)
        self.number_of_items += 1

    cpdef tuple top(self):
        r"""
        Return the top pair (item, value) of the heap.

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_hashables
            sage: P = PairingHeap_of_n_hashables(5)
            sage: P.push(1, 2)
            sage: P.top()
            (1, 2)
            sage: P.push(3, 1)
            sage: P.top()
            (3, 1)

            sage: P = PairingHeap_of_n_hashables(3)
            sage: P.top()
            Traceback (most recent call last):
            ...
            ValueError: trying to access the top of an empty heap
        """
        if not self:
            raise ValueError("trying to access the top of an empty heap")
        cdef size_t idx = self.root - self.nodes
        return self._int_to_item[idx], <object>self.root.value

    cpdef object top_item(self):
        r"""
        Return the top item of the heap.

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_hashables
            sage: P = PairingHeap_of_n_hashables(5)
            sage: P.push(1, 2)
            sage: P.top()
            (1, 2)
            sage: P.top_item()
            1

            sage: P = PairingHeap_of_n_hashables(3)
            sage: P.top_item()
            Traceback (most recent call last):
            ...
            ValueError: trying to access the top of an empty heap
        """
        if not self:
            raise ValueError("trying to access the top of an empty heap")
        cdef size_t idx = self.root - self.nodes
        return self._int_to_item[idx]

    cpdef void pop(self) noexcept:
        r"""
        Remove the top item from the heap.

        If the heap is already empty, we do nothing.

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_hashables
            sage: P = PairingHeap_of_n_hashables(5); len(P)
            0
            sage: P.push(1, 2); len(P)
            1
            sage: P.push(2, 3); len(P)
            2
            sage: P.pop(); len(P)
            1
            sage: P.pop(); len(P)
            0
            sage: P.pop(); len(P)
            0
        """
        if not self:
            return
        cdef object item = self.top_item()
        cdef size_t idx = self._item_to_int[item]
        Py_XDECREF(<PyObject *>self.nodes[idx].value)
        self.free_idx.append(idx)
        del self._item_to_int[item]
        self.number_of_items -= 1
        self.root = PairingHeapNodePy._pair(self.root.child)

    cpdef void decrease(self, object item, object new_value) except *:
        r"""
        Decrease the value of specified item.

        This method is more permissive than it should as it can also be used to
        push an item in the heap.

        INPUT:

        - ``item`` -- the item to consider

        - ``new_value`` -- the new value for ``item``

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_hashables
            sage: P = PairingHeap_of_n_hashables(5)
            sage: 3 in P
            False
            sage: P.decrease(3, 33)
            sage: 3 in P
            True
            sage: P.top()
            (3, 33)
            sage: P.push(1, 10)
            sage: P.top()
            (1, 10)
            sage: P.decrease(3, 7)
            sage: P.top()
            (3, 7)

        TESTS::

            sage: P = PairingHeap_of_n_hashables(5)
            sage: P.push(1, 3)
            sage: P.decrease(1, 2)
            sage: P.decrease(1, 2)
            Traceback (most recent call last):
            ...
            ValueError: the new value must be less than the current value
        """
        cdef PairingHeapNodePy * p
        cdef size_t idx
        if item in self:
            idx = self._item_to_int[item]
            p = self.nodes + idx
            if <object>p.value <= new_value:
                raise ValueError("the new value must be less than the current value")
            Py_XDECREF(<PyObject *>p.value)
            Py_INCREF(new_value)
            p.value = <PyObject *>new_value
            if p.prev != NULL:
                PairingHeapNodePy._unlink(p)
                self.root = PairingHeapNodePy._merge(self.root, p)
        else:
            self.push(item, new_value)

    def __contains__(self, object item):
        r"""
        Check whether the specified item is in the heap.

        INPUT:

        - ``item`` -- the item to consider

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_hashables
            sage: P = PairingHeap_of_n_hashables(5)
            sage: 3 in P
            False
            sage: P.push(3, 33)
            sage: 3 in P
            True
            sage: 100 in P
            False
        """
        return item in self._item_to_int

    contains = __contains__

    cpdef object value(self, object item):
        r"""
        Return the value associated with the item.

        INPUT:

        - ``item`` -- the item to consider

        EXAMPLES::

            sage: from sage.data_structures.pairing_heap import PairingHeap_of_n_hashables
            sage: P = PairingHeap_of_n_hashables(5)
            sage: P.push(3, 33)
            sage: P.push(1, 10)
            sage: P.value(3)
            33
            sage: P.value(7)
            Traceback (most recent call last):
            ...
            ValueError: 7 is not in the heap
        """
        if item not in self:
            raise ValueError(f"{item} is not in the heap")
        cdef size_t idx = self._item_to_int[item]
        return <object>self.nodes[idx].value


# ==============================================================================
# Methods to check the validity of the pairing heaps
# ==============================================================================

def _test_PairingHeap_from_C(n=100):
    r"""
    Test :class:`~sage.data_structures.pairing_heap.PairingHeap`.

    INPUT:

    - ``n`` -- a strictly positive integer (default: 100); the maximum capacity
      of the heap

    TESTS::

        sage: from sage.data_structures.pairing_heap import _test_PairingHeap_from_C
        sage: _test_PairingHeap_from_C(100)
    """
    from sage.misc.prandom import randint, shuffle
    sig_on()
    cdef PairingHeap[size_t, size_t] PH = PairingHeap[size_t, size_t]()
    sig_off()

    # Initialize a list of tuples (value, item) randomly ordered
    items = list(range(n))
    values = list(range(n))
    shuffle(items)
    shuffle(values)
    cdef list Lref = list(zip(values, items))

    for value, item in Lref:
        PH.push(item, value)
        sig_check()

    L = []
    while not PH.empty():
        item, value = PH.top()
        L.append((value, item))
        PH.pop()
        sig_check()

    if L != sorted(Lref):
        raise ValueError('the order is not good')

    # Test decrease key operations. We first push items in the heap with an
    # excess of k in the value. Then we decrease the keys in a random order by
    # random values until returning to the origianl values. We finally check the
    # validity of the resulting ordering.
    k = 10
    dec = {item: k for item in items}
    shuffle(Lref)
    for value, item in Lref:
        PH.push(item, value + k)
        sig_check()

    L = list(items)
    while L:
        i = randint(0, len(L) - 1)
        item = L[i]
        d = randint(1, dec[item])
        dec[item] -= d
        if not dec[item]:
            L[i] = L[-1]
            L.pop()
        PH.decrease(item, PH.value(item) - d)
        sig_check()

    L = []
    while not PH.empty():
        item, value = PH.top()
        L.append((value, item))
        PH.pop()
        sig_check()

    if L != sorted(Lref):
        raise ValueError('the order is not good')

    # Different cost function
    from sage.functions.trig import sin, cos
    sig_on()
    cdef PairingHeap[size_t, pair[size_t, size_t]] HH = PairingHeap[size_t, pair[size_t, size_t]]()
    sig_off()

    for i in range(n):
        HH.push(i, (sin(i), cos(i)))
        sig_check()

    L = []
    while not HH.empty():
        L.append(HH.top())
        HH.pop()
        sig_check()

    for (u, cu), (v, cv) in zip(L, L[1:]):
        if cu > cv:
            print(u, cu, v, cv)

    # We finally show that an error is raised when trying to access the top of
    # an empty heap
    try:
        _ = HH.top()
        print("something goes wrong, the error has not been raised")
    except ValueError as msg:
        # The error has been properly handled
        pass

    try:
        _ = HH.top_item()
        print("something goes wrong, the error has not been raised")
    except ValueError as msg:
        # The error has been properly handled
        pass

    try:
        _ = HH.top_value()
        print("something goes wrong, the error has not been raised")
    except ValueError as msg:
        # The error has been properly handled
        pass

    # Or to get the value associated to an item that is not in the heap
    try:
        _ = HH.value(123)
        print("something goes wrong, the error has not been raised")
    except ValueError as msg:
        # The error has been properly handled
        pass


def _test_PairingHeap_of_n_integers(n=100):
    r"""
    Test :class:`~sage.data_structures.pairing_heap.PairingHeap_of_n_integers`.

    INPUT:

    - ``n`` -- a strictly positive integer (default: 100); the maximum capacity
      of the heap

    TESTS::

        sage: from sage.data_structures.pairing_heap import _test_PairingHeap_of_n_integers
        sage: _test_PairingHeap_of_n_integers(100)
    """
    from sage.misc.prandom import randint, shuffle

    sig_on()
    cdef PairingHeap_of_n_integers P = PairingHeap_of_n_integers(n)
    sig_off()

    # Initialize a list of tuples (value, item) randomly ordered
    cdef list items = list(range(n))
    cdef list values = list(range(n))
    shuffle(items)
    shuffle(values)
    cdef list Lref = list(zip(values, items))

    cdef int value, item
    for value, item in Lref:
        P.push(item, value)
        sig_check()

    cdef list L = []
    while P:
        item, value = P.top()
        L.append((value, item))
        P.pop()
        sig_check()

    if L != sorted(Lref):
        raise ValueError('the order is not good')

    # Test decrease key operations. We first push items in the heap with an
    # excess of k in the value. Then we decrease the keys in a random order by
    # random values until returning to the origianl values. We finally check the
    # validity of the resulting ordering.
    cdef int k = 10
    cdef list dec = [k] * n
    shuffle(Lref)
    for value, item in Lref:
        P.push(item, value + k)
        sig_check()

    L = list(items)
    while L:
        i = randint(0, len(L) - 1)
        item = L[i]
        d = randint(1, dec[item])
        dec[item] -= d
        if not dec[item]:
            L[i] = L[-1]
            L.pop()
        P.decrease(item, P.value(item) - d)
        sig_check()

    L = []
    while P:
        item, value = P.top()
        L.append((value, item))
        P.pop()
        sig_check()

    if L != sorted(Lref):
        raise ValueError('the order is not good')

    # We finally show that an error is raised when trying to access the top of
    # an empty heap
    try:
        _ = P.top()
        print("something goes wrong, the error has not been raised")
    except ValueError as msg:
        # The error has been properly handled
        pass

    try:
        _ = P.top_item()
        print("something goes wrong, the error has not been raised")
    except ValueError as msg:
        # The error has been properly handled
        pass

    try:
        _ = P.top_value()
        print("something goes wrong, the error has not been raised")
    except ValueError as msg:
        # The error has been properly handled
        pass

    # Or to get the value associated to an item that is not in the heap
    try:
        _ = P.value(123)
        print("something goes wrong, the error has not been raised")
    except ValueError as msg:
        # The error has been properly handled
        pass


def _test_PairingHeap_of_n_hashables(n=100):
    r"""
    Test :class:`~sage.data_structures.pairing_heap.PairingHeap_of_n_hashables`.

    INPUT:

    - ``n`` -- a strictly positive integer (default: 100); the maximum capacity
      of the heap

    TESTS::

        sage: from sage.data_structures.pairing_heap import _test_PairingHeap_of_n_hashables
        sage: _test_PairingHeap_of_n_hashables(100)
    """
    from sage.misc.prandom import randint, shuffle

    sig_on()
    cdef PairingHeap_of_n_hashables P = PairingHeap_of_n_hashables(n)
    sig_off()

    # Initialize a list of tuples (value, item) randomly ordered
    cdef list items = [(str(i), i) for i in range(n)]
    cdef list values = list(range(n))
    shuffle(items)
    shuffle(values)
    cdef list Lref = list(zip(values, items))

    for value, item in Lref:
        P.push(item, value)
        sig_check()

    cdef list L = []
    while P:
        item, value = P.top()
        L.append((value, item))
        P.pop()
        sig_check()

    if L != sorted(Lref):
        raise ValueError('the order is not good')

    # Test decrease key operations. We first push items in the heap with an
    # excess of k in the value. Then we decrease the keys in a random order by
    # random values until returning to the origianl values. We finally check the
    # validity of the resulting ordering.
    cdef int k = 10
    cdef dict dec = {item: k for item in items}
    shuffle(Lref)
    for value, item in Lref:
        P.push(item, value + k)
        sig_check()

    L = list(items)
    while L:
        i = randint(0, len(L) - 1)
        item = L[i]
        d = randint(1, dec[item])
        dec[item] -= d
        if not dec[item]:
            L[i] = L[-1]
            L.pop()
        P.decrease(item, P.value(item) - d)
        sig_check()

    L = []
    while P:
        item, value = P.top()
        L.append((value, item))
        P.pop()
        sig_check()

    if L != sorted(Lref):
        raise ValueError('the order is not good')

    # We finally show that an error is raised when trying to access the top of
    # an empty heap
    try:
        _ = P.top()
        print("something goes wrong, the error has not been raised")
    except ValueError as msg:
        # The error has been properly handled
        pass

    try:
        _ = P.top_item()
        print("something goes wrong, the error has not been raised")
    except ValueError as msg:
        # The error has been properly handled
        pass

    try:
        _ = P.top_value()
        print("something goes wrong, the error has not been raised")
    except ValueError as msg:
        # The error has been properly handled
        pass

    # Or to get the value associated to an item that is not in the heap
    try:
        _ = P.value(123)
        print("something goes wrong, the error has not been raised")
    except ValueError as msg:
        # The error has been properly handled
        pass


def compare_heaps(n=100, verbose=False):
    r"""
    Check that the heaps behave the same.

    This method selects a list of instructions: push items in some order,
    decrease the values of the items in some order, extract all items in order.
    Then it applies the same instructions to a ``PairingHeap``, a
    :class:`PairingHeap_of_n_integers` and a
    :class:`PairingHeap_of_n_hashables`. It checks that all heaps report the
    list of items in the same order and it measures the running time.

    INPUT:

    - ``n`` -- a strictly positive integer (default: 100); the maximum capacity
      of the heap

    - ``verbose`` -- boolean (default: ``False``); whether to display
      information about the running times

    EXAMPLES::

        sage: from sage.data_structures.pairing_heap import compare_heaps
        sage: compare_heaps(n=100)
        sage: compare_heaps(n=100, verbose=True)  # random
        PairingHeap_of_n_integers: 7.800000000024454e-05
        PairingHeap_of_n_hashables: 9.400000000026054e-05
        PairingHeap (C++): 6.899999999987472e-05
        sage: compare_heaps(1000000, verbose=True)  # not tested (long time), random
        PairingHeap_of_n_integers: 1.5106779999999995
        PairingHeap_of_n_hashables: 4.998040000000001
        PairingHeap (C++): 1.7841750000000012
    """
    from sage.misc.prandom import shuffle
    from sage.misc.timing import cputime

    items = list(range(n))
    values = list(range(n))
    shuffle(items)
    shuffle(values)
    Lref = list(zip(values, items))
    k = 10
    dec_order = list(items)
    shuffle(dec_order)

    t = cputime()
    sig_on()
    cdef PairingHeap_of_n_integers P = PairingHeap_of_n_integers(n)
    sig_off()
    for value, item in Lref:
        P.push(item, value + k)
        sig_check()
    for item in dec_order:
        P.decrease(item, P.value(item) - k)
        sig_check()
    LP = []
    while P:
        LP.append(P.top())
        P.pop()
        sig_check()
    t = cputime(t)
    if verbose:
        print(f"PairingHeap_of_n_integers: {t}")

    t = cputime()
    sig_on()
    cdef PairingHeap_of_n_hashables Q = PairingHeap_of_n_hashables(n)
    sig_off()
    for value, item in Lref:
        Q.push(item, value + k)
        sig_check()
    for item in dec_order:
        Q.decrease(item, Q.value(item) - k)
        sig_check()
    LQ = []
    while Q:
        LQ.append(Q.top())
        Q.pop()
        sig_check()
    t = cputime(t)
    if verbose:
        print(f"PairingHeap_of_n_hashables: {t}")

    t = cputime()
    sig_on()
    cdef PairingHeap[size_t, size_t] PH = PairingHeap[size_t, size_t]()
    sig_off()
    for value, item in Lref:
        PH.push(item, value + k)
        sig_check()
    for item in dec_order:
        PH.decrease(item, PH.value(item) - k)
        sig_check()
    LPH = []
    while not PH.empty():
        LPH.append(PH.top())
        PH.pop()
        sig_check()
    t = cputime(t)
    if verbose:
        print(f"PairingHeap (C++): {t}")

    if LPH != LP or LP != LQ:
        print('something goes wrong')
