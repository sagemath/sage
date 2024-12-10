# ******************************************************************************
# Copyright (C) 2024 David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ******************************************************************************

# ==============================================================================
# Interface to pairing heap data structure from ./pairing_heap.h
# ==============================================================================

from libcpp.pair cimport pair
from cpython cimport PyObject

cdef extern from "./pairing_heap.h" namespace "pairing_heap":
    cdef cppclass PairingHeap[TypeOfItem, TypeOfValue]:
        PairingHeap() except +
        PairingHeap(PairingHeap[TypeOfItem, TypeOfValue]) except +
        bint empty()
        void push(TypeOfItem, TypeOfValue) except +
        pair[TypeOfItem, TypeOfValue] top() except +
        TypeOfItem top_item() except +
        TypeOfValue top_value() except +
        void pop() except +
        void decrease(TypeOfItem, TypeOfValue) except +
        bint contains(TypeOfItem)
        TypeOfValue value(TypeOfItem) except +

    cdef cppclass PairingHeapNodePy:
        PyObject * value           # value associated with the item
        PairingHeapNodePy * prev   # Previous sibling of the node or parent
        PairingHeapNodePy * next   # Next sibling of the node
        PairingHeapNodePy * child  # First child of the node

        @staticmethod
        PairingHeapNodePy * _merge(PairingHeapNodePy * a, PairingHeapNodePy * b) except +
        @staticmethod
        PairingHeapNodePy * _pair(PairingHeapNodePy * p) except +
        @staticmethod
        void _link(PairingHeapNodePy * a, PairingHeapNodePy * b)
        @staticmethod
        void _unlink(PairingHeapNodePy * p)


# ==============================================================================
# Pairing heap data structure with fixed capacity n
# ==============================================================================

from sage.data_structures.bitset_base cimport bitset_t

cdef class PairingHeap_class:
    cdef size_t n                   # maximum number of items
    cdef PairingHeapNodePy * root   # pointer to the top of the heap
    cdef PairingHeapNodePy * nodes  # array of size n to store items
    cdef size_t number_of_items     # number of active items
    cpdef bint empty(self) noexcept
    cpdef bint full(self) noexcept
    cpdef tuple top(self)
    cpdef object top_value(self)
    cpdef void pop(self) noexcept


cdef class PairingHeap_of_n_integers(PairingHeap_class):
    cdef bitset_t active  # bitset to identify active items
    cpdef void push(self, size_t item, object value) except *
    cpdef size_t top_item(self) except *
    cpdef void decrease(self, size_t item, object new_value) except *
    cpdef object value(self, size_t item)


cdef class PairingHeap_of_n_hashables(PairingHeap_class):
    cdef list _int_to_item  # mapping from integers to items
    cdef dict _item_to_int  # mapping from items to integers
    cdef list free_idx      # list of free indexes
    cpdef void push(self, object item, object value) except *
    cpdef object top_item(self)
    cpdef void decrease(self, object item, object new_value) except *
    cpdef object value(self, object item)
