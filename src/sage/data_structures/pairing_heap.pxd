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

cdef extern from "./pairing_heap.h" namespace "pairing_heap":
    cdef cppclass PairingHeap[TypeOfItem, TypeOfValue]:
        PairingHeap() except +
        PairingHeap(PairingHeap[TypeOfItem, TypeOfValue]) except +
        bint empty()
        void reset()
        void push(TypeOfItem, TypeOfValue) except +
        pair[TypeOfItem, TypeOfValue] top() except +
        TypeOfItem top_item() except +
        TypeOfValue top_value() except +
        void pop() except +
        void decrease(TypeOfItem, TypeOfValue) except +
        bint contains(TypeOfItem)
        TypeOfValue value(TypeOfItem) except +

# ==============================================================================
# Pairing heap data structure with fixed capacity n
# ==============================================================================

from sage.data_structures.bitset_base cimport bitset_t

ctypedef struct PairingHeapNode:
    void * value             # value associated with the item
    PairingHeapNode * prev   # Previous sibling of the node or parent
    PairingHeapNode * succ   # Next sibling of the node
    PairingHeapNode * child  # First child of the node


cdef bint _compare(PairingHeapNode * a, PairingHeapNode * b) except *
cdef PairingHeapNode * _pair(PairingHeapNode * p) except *
cdef PairingHeapNode * _merge(PairingHeapNode * a, PairingHeapNode * b) except *
cdef _link(PairingHeapNode * a, PairingHeapNode * b) except *
cdef _unlink(PairingHeapNode * p) except *


cdef class PairingHeap_class:
    cdef str name                 # name of the data structure
    cdef size_t n                 # maximum number of items
    cdef PairingHeapNode * root   # pointer to the top of the heap
    cdef PairingHeapNode * nodes  # array of size n to store items
    cdef bitset_t active          # bitset to identify active items
    cdef size_t number_of_items   # number of active items
    cpdef bint empty(self) noexcept
    cpdef bint full(self) noexcept


cdef class PairingHeap_of_n_integers(PairingHeap_class):
    cpdef void push(self, size_t item, object value) except *
    cpdef tuple top(self) except *
    cpdef size_t top_item(self) except *
    cpdef object top_value(self) except *
    cpdef void pop(self) noexcept
    cpdef void decrease(self, size_t item, object new_value) except *
    cpdef object value(self, size_t item) except *


cdef class PairingHeap_of_n_hashables(PairingHeap_class):
    cdef list _int_to_item  # mapping from integers to items
    cdef dict _item_to_int  # mapping from items to integers
    cpdef void push(self, object item, object value) except *
    cpdef tuple top(self) except *
    cpdef object top_item(self) except *
    cpdef object top_value(self) except *
    cpdef void pop(self) noexcept
    cpdef void decrease(self, object item, object new_value) except *
    cpdef object value(self, object item) except *
