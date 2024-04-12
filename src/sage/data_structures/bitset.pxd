# sage_setup: distribution = sagemath-categories
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
    cdef FrozenBitset _new(self,long int capacity)
    cpdef FrozenBitset _larger_capacity_(self, long size)
    cpdef long capacity(self) noexcept
    cpdef bint isempty(self) noexcept
    cpdef bint issubset(self, FrozenBitset other) except -1
    cpdef bint issuperset(self, FrozenBitset other) except -1
    cpdef bint isdisjoint(self, FrozenBitset other) except -1
    cpdef _union(self, FrozenBitset other)
    cpdef intersection(self, FrozenBitset other)
    cpdef difference(self, FrozenBitset other)
    cpdef symmetric_difference(self, FrozenBitset other)
    cpdef complement(self)
    cpdef __copy__(self)

cdef class Bitset(FrozenBitset):
    cpdef __copy__(self)
    cpdef update(self, FrozenBitset other)
    cpdef intersection_update(self, FrozenBitset other)
    cpdef difference_update(self, FrozenBitset other)
    cpdef symmetric_difference_update(self, FrozenBitset other)
    cpdef add(self, unsigned long n)
    cpdef remove(self, unsigned long n)
    cpdef discard(self, unsigned long n)
    cpdef pop(self)
    cpdef clear(self)

