"""
Cython bitset types
"""
#*****************************************************************************
#     Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# This file declares the bitset types. The implementation of the basic
# Cython type bitset_t (inline functions) is in bitset.pxi, the
# implementation of the Python later is in bitset.pyx. The latter file
# also contains all doctests.

cdef extern from *:
    # constant literals
    int index_shift "(sizeof(unsigned long)==8 ? 6 : 5)"
    unsigned long offset_mask "(sizeof(unsigned long)==8 ? 0x3F : 0x1F)"

    # Given an element index n in a set, (n >> index_shift) gives the
    # corresponding limb number, while (n & offset_mask) gives the bit
    # number inside of the limb.

cdef struct bitset_s:
    # The size of a bitset B counts the maximum number of bits that B can
    # hold. This size is independent of how many elements of B are toggled to
    # 1. For example, say B is the bitset 1001. Then B has size 4, with the
    # first and fourth elements toggled to 1, reading from left to right.
    # We can also think of the size of a bitset as its capacity.
    long size

    # A limb is that part of a bitset that can fit into an unsigned long
    # (typically, 32 bits on a 32-bit machine and 64 bits on a 64-bit
    # machine). This counts the number of limbs to represent a bitset.
    # If a bitset has size <= n, then the whole bitset fits into a limb
    # and we only require one limb to represent the bitset. However, if
    # the bitset has size > n, we require more than one limb to
    # represent the bitset. For example, if a limb is 64 bits in length
    # and the bitset has size 96 bits, then we require at most two limbs
    # to represent the bitset.
    #
    # NOTE: we assume that a limb corresponds to an MPIR limb (this
    # assumption is always true in practice).
    long limbs

    # The individual bits of a bitset.
    unsigned long *bits

ctypedef bitset_s bitset_t[1]


# Python layer over bitset_t
cdef class FrozenBitset:
    cdef bitset_t _bitset
    cdef FrozenBitset _new(self,long int capacity)
    cpdef FrozenBitset _larger_capacity_(self, long size)
    cpdef long capacity(self)
    cpdef bint isempty(self)
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

