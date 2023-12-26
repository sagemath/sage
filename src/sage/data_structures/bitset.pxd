#*****************************************************************************
#     Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.data_structures.bitset_base cimport bitset_t

# Python layer over bitset_t
cdef class FrozenBitset:
    cdef bitset_t _bitset
    cdef FrozenBitset _new(self,long int capacity) noexcept
    cpdef FrozenBitset _larger_capacity_(self, long size) noexcept
    cpdef long capacity(self) noexcept
    cpdef bint isempty(self) noexcept
    cpdef bint issubset(self, FrozenBitset other) except -1
    cpdef bint issuperset(self, FrozenBitset other) except -1
    cpdef bint isdisjoint(self, FrozenBitset other) except -1
    cpdef _union(self, FrozenBitset other) noexcept
    cpdef intersection(self, FrozenBitset other) noexcept
    cpdef difference(self, FrozenBitset other) noexcept
    cpdef symmetric_difference(self, FrozenBitset other) noexcept
    cpdef complement(self) noexcept
    cpdef __copy__(self) noexcept

cdef class Bitset(FrozenBitset):
    cpdef __copy__(self) noexcept
    cpdef update(self, FrozenBitset other) noexcept
    cpdef intersection_update(self, FrozenBitset other) noexcept
    cpdef difference_update(self, FrozenBitset other) noexcept
    cpdef symmetric_difference_update(self, FrozenBitset other) noexcept
    cpdef add(self, unsigned long n) noexcept
    cpdef remove(self, unsigned long n) noexcept
    cpdef discard(self, unsigned long n) noexcept
    cpdef pop(self) noexcept
    cpdef clear(self) noexcept

