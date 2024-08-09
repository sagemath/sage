# distutils: depends = bitset_intrinsics.h
# distutils: libraries = gmp
"""
A fast bitset datatype in Cython

Operations between bitsets are only guaranteed to work if the bitsets
have the same size, with the exception of ``bitset_realloc``.  Similarly, you
should not try to access elements of a bitset beyond the size.

AUTHORS:

- Robert Bradshaw (2008)
- Rudi Pendavingh, Stefan van Zwam (2013-06-06): added functions map, lex_cmp,
  pickle, unpickle
- Jeroen Demeyer (2014-09-05): use mpn_* functions from MPIR in the
  implementation (:issue:`13352` and :issue:`16937`)
- Simon King (2014-10-28): ``bitset_rshift`` and ``bitset_lshift`` respecting
  the size of the given bitsets (:issue:`15820`)
"""

#*****************************************************************************
#     Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# This file declares the bitset types and the (inline functions).
# Few functions that are not inline are in the `pyx` file.
# The python wrapper and all the doctests are in bitset.pyx/bitset.pxd.

from libc.string cimport strlen
from cysignals.memory cimport check_calloc, check_allocarray, check_reallocarray, sig_malloc, sig_free
from memory_allocator  cimport MemoryAllocator
from memory_allocator.memory_allocator cimport align
from cython.operator import preincrement as preinc

from sage.cpython.string cimport str_to_bytes, bytes_to_str
from sage.libs.gmp.mpn cimport *
from sage.libs.gmp.types cimport *
from sage.data_structures.sparse_bitset cimport sparse_bitset_t


cdef extern from *:
    # Given an element index n in a set, (n >> index_shift) gives the
    # corresponding limb number.
    int index_shift "(sizeof(mp_limb_t) == 8 ? 6 : 5)"


cdef struct bitset_s:
    # The size of a bitset B counts the maximum number of bits that B can
    # hold. This size is independent of how many elements of B are toggled to
    # 1. For example, say B is the bitset 1001. Then B has size 4, with the
    # first and fourth elements toggled to 1, reading from left to right.
    # We can also think of the size of a bitset as its capacity.
    mp_bitcnt_t size

    # A limb is that part of a bitset that can fit into an mp_limb_t
    # (typically, 32 bits on a 32-bit machine and 64 bits on a 64-bit
    # machine). This counts the number of limbs to represent a bitset.
    # If a bitset has size <= n, then the whole bitset fits into a limb
    # and we only require one limb to represent the bitset. However, if
    # the bitset has size > n, we require more than one limb to
    # represent the bitset. For example, if a limb is 64 bits in length
    # and the bitset has size 96 bits, then we require at most two limbs
    # to represent the bitset.
    #
    # NOTE: some code assumes that mp_limb_t is an unsigned long
    # (this assumption is always true in practice).
    mp_size_t limbs

    # The individual bits of a bitset.
    mp_limb_t* bits

ctypedef bitset_s bitset_t[1]

ctypedef fused fused_bitset_t:
    bitset_t
    sparse_bitset_t

# Doctests for the functions in this file are in sage/data_structures/bitset.pyx

#############################################################################
# Functions that can be optimized by intrinsics.
#############################################################################

cdef extern from "bitset_intrinsics.h":
    cdef const mp_bitcnt_t LIMB_SIZE
    cdef const mp_bitcnt_t ALIGNMENT

    cdef mp_bitcnt_t _set_non_zero(
            mp_limb_t* bits, mp_bitcnt_t* non_zero_chunks, mp_bitcnt_t limbs) nogil

    # Bitset Comparison
    ctypedef enum cmpop_t:
        EQUAL
        SUBSET
        DISJOINT
    cdef bint _bitset_isempty(mp_limb_t* bits, mp_bitcnt_t limbs) nogil
    cdef bint _bitset_cmp(mp_limb_t* a, mp_limb_t* b, mp_bitcnt_t limbs, cmpop_t) nogil

    # Bitset Comparison for sparse bitsets (only subset and disjoint)
    cdef bint _sparse_bitset_cmp(
            mp_limb_t* a, mp_bitcnt_t* a_non_zero_chunks, mp_bitcnt_t a_n_non_zero_chunks,
            mp_limb_t* b, cmpop_t) nogil

    # Bitset Searching
    cdef long _bitset_first_in_limb(mp_limb_t limb) nogil
    cdef long _bitset_first_in_limb_nonzero(mp_limb_t) nogil
    cdef long _bitset_len(mp_limb_t* bits, mp_bitcnt_t limbs) nogil

    # Bitset Arithmetic
    ctypedef enum operation_t:
        AND
        OR
        ANDNOT
        XOR
    cdef void _bitset_operation(
            mp_limb_t* dst, mp_limb_t* a, mp_limb_t* b, mp_bitcnt_t limbs, operation_t) nogil

    # Bitset Arithmetic for sparse bitsets.
    cdef mp_bitcnt_t _sparse_bitset_operation(
            mp_limb_t* dst, mp_bitcnt_t* dst_non_zero_chunks, mp_limb_t* a,
            mp_limb_t* b, mp_bitcnt_t limbs, operation_t) nogil

#############################################################################
# Creating limb patterns
#############################################################################
#
# NOTE: In all functions in this section, the index n is interpreted
# modulo GMP_LIMB_BITS, the number of bits in a limb.
#
cdef inline mp_limb_t limb_one_set_bit(mp_bitcnt_t n) noexcept:
    """
    Return a limb with only bit n set.
    """
    return (<mp_limb_t>1) << (n % GMP_LIMB_BITS)

cdef inline mp_limb_t limb_one_zero_bit(mp_bitcnt_t n) noexcept:
    """
    Return a limb with all bits set, except for bit n.
    """
    return ~((<mp_limb_t>1) << (n % GMP_LIMB_BITS))

cdef inline mp_limb_t limb_lower_bits_down(mp_bitcnt_t n) noexcept:
    """
    Return a limb with the lower n bits set, where n is interpreted
    in [0 .. GMP_LIMB_BITS-1].
    """
    return ((<mp_limb_t>1) << (n % GMP_LIMB_BITS)) - 1

cdef inline mp_limb_t limb_lower_bits_up(mp_bitcnt_t n) noexcept:
    """
    Return a limb with the lower n bits set, where n is interpreted
    in [1 .. GMP_LIMB_BITS].
    """
    return (<mp_limb_t>(-1)) >> ((<unsigned int>(-n)) % GMP_LIMB_BITS)

#############################################################################
# Bitset Initalization
#############################################################################
cdef inline bint bitset_init(fused_bitset_t bits, mp_bitcnt_t size) except -1:
    """
    Allocate an empty bitset of size ``size``.

    Size must be at least 1.
    """
    if size <= 0:
        raise ValueError("bitset capacity must be greater than 0")

    cdef mp_bitcnt_t extra

    bits.size = size
    if fused_bitset_t is bitset_t:
        bits.limbs = (size - 1) // (8 * LIMB_SIZE) + 1
        bits.bits = <mp_limb_t*>check_calloc(bits.limbs, LIMB_SIZE)
    else:
        bits.limbs = ((size - 1) // (8*ALIGNMENT) + 1) * (ALIGNMENT//LIMB_SIZE)
        extra = (ALIGNMENT + LIMB_SIZE - 2) // LIMB_SIZE
        bits.mem = check_calloc(bits.limbs + extra, LIMB_SIZE)
        bits.bits = <mp_limb_t*> align(bits.mem, ALIGNMENT)
        bits.non_zero_chunks_are_initialized = False
        bits.non_zero_chunks = <mp_bitcnt_t*> check_allocarray((bits.limbs*LIMB_SIZE) // ALIGNMENT, sizeof(mp_bitcnt_t))

cdef inline bint bitset_check_alignment(fused_bitset_t bits) noexcept:
    """
    Return whether the bitset is aligned correctly.
    """
    cdef size_t address = <size_t> bits.bits
    return address == (address & ~(ALIGNMENT - 1))

cdef inline int bitset_realloc(bitset_t bits, mp_bitcnt_t size) except -1:
    """
    Reallocate a bitset to size ``size``. If reallocation is larger, new bitset
    does not contain any of the extra bits.
    """
    cdef mp_size_t limbs_old = bits.limbs
    cdef mp_bitcnt_t size_old = bits.size
    if size_old == size:
        return 0
    if size <= 0:
        raise ValueError("bitset capacity must be greater than 0")

    cdef mp_size_t limbs_new = (size - 1) // (8 * LIMB_SIZE) + 1
    bits.bits = <mp_limb_t*>check_reallocarray(bits.bits, limbs_new, LIMB_SIZE)
    bits.size = size
    bits.limbs = limbs_new

    if bits.limbs > limbs_old:
        # Zero any extra limbs
        mpn_zero(bits.bits + limbs_old, bits.limbs - limbs_old)
    elif bits.size < size_old:
        # Zero removed bits
        bitset_fix(bits)

cdef inline void bitset_free(fused_bitset_t bits) noexcept:
    """
    Deallocate the memory in bits.
    """
    if fused_bitset_t is bitset_t:
        sig_free(bits.bits)
    else:
        sig_free(bits.mem)
        sig_free(bits.non_zero_chunks)

cdef inline void bitset_clear(fused_bitset_t bits) noexcept:
    """
    Remove all elements from the set.
    """
    mpn_zero(bits.bits, bits.limbs)
    if fused_bitset_t is sparse_bitset_t:
        bits.non_zero_chunks_are_initialized = False

cdef inline void bitset_zero(fused_bitset_t bits) noexcept:
    """
    Remove all elements from the set.

    This function is the same as bitset_clear(bits).
    """
    mpn_zero(bits.bits, bits.limbs)
    if fused_bitset_t is sparse_bitset_t:
        bits.non_zero_chunks_are_initialized = False

cdef inline void bitset_copy(fused_bitset_t dst, fused_bitset_t src) noexcept:
    """
    Copy the bitset src over to the bitset dst, overwriting dst.

    We assume ``dst.limbs == src.limbs``.
    """
    mpn_copyi(dst.bits, src.bits, src.limbs)
    if fused_bitset_t is sparse_bitset_t:
        dst.non_zero_chunks_are_initialized = False

cdef inline void bitset_copy_flex(fused_bitset_t dst, fused_bitset_t src) noexcept:
    """
    Copy the bitset src over to the bitset dst, overwriting dst.

    Set additional limbs in dst to zero.

    We assume ``dst.limbs >= src.limbs``.
    """
    mpn_copyi(dst.bits, src.bits, src.limbs)
    mpn_zero(dst.bits + src.limbs, dst.limbs - src.limbs)
    if fused_bitset_t is sparse_bitset_t:
        dst.non_zero_chunks_are_initialized = False

cdef inline void bitset_fix(fused_bitset_t bits) noexcept:
    """
    Clear upper bits in upper limb which should be zero.
    """
    bits.bits[bits.limbs - 1] &= limb_lower_bits_up(bits.size)

cdef inline void sparse_bitset_set_non_zero(sparse_bitset_t bits) noexcept nogil:
    """
    Set the nonzero chunks of ``bits``.
    """
    bits.n_non_zero_chunks = _set_non_zero(bits.bits, bits.non_zero_chunks, bits.limbs)
    bits.non_zero_chunks_are_initialized = True

#############################################################################
# Bitset Comparison
#############################################################################

cdef inline bint mpn_equal_bits(mp_srcptr b1, mp_srcptr b2, mp_bitcnt_t n) noexcept:
    """
    Return ``True`` iff the first n bits of *b1 and *b2 agree.
    """
    cdef mp_size_t nlimbs = n // GMP_LIMB_BITS
    cdef mp_limb_t mask = limb_lower_bits_down(n)
    if nlimbs > 0 and mpn_cmp(b1, b2, nlimbs) != 0:
        return False
    if mask == 0:
        return True

    cdef mp_limb_t b1h = b1[nlimbs]
    cdef mp_limb_t b2h = b2[nlimbs]
    return (b1h ^ b2h) & mask == 0

cdef inline bint mpn_equal_bits_shifted(mp_srcptr b1, mp_srcptr b2, mp_bitcnt_t n, mp_bitcnt_t offset) noexcept:
    """
    Return ``True`` iff the first n bits of *b1 and the bits ranging from
    offset to offset+n of *b2 agree.
    """
    cdef mp_bitcnt_t bit_offset = offset % GMP_LIMB_BITS
    cdef mp_size_t i2 = offset//GMP_LIMB_BITS
    if bit_offset==0:
        return mpn_equal_bits(b1, b2 + i2, n)
    cdef mp_size_t neg_bit_offset = GMP_LIMB_BITS-bit_offset
    # limbs of b1 to be considered
    cdef mp_size_t nlimbs = n // GMP_LIMB_BITS
    # bits of an additional limb of b1 to be considered
    cdef mp_limb_t tmp_limb
    cdef mp_size_t i1
    for i1 in range(nlimbs):
        tmp_limb = (b2[i2] >> bit_offset)
        tmp_limb |= (b2[preinc(i2)] << neg_bit_offset)
        if tmp_limb != b1[i1]:
            return False
    cdef mp_limb_t mask = limb_lower_bits_down(n)
    if mask == 0:
        return True

    cdef mp_limb_t b1h = b1[nlimbs]
    tmp_limb = (b2[i2] >> bit_offset)
    if (n%GMP_LIMB_BITS)+bit_offset > GMP_LIMB_BITS:
        # Need bits from the next limb of b2
        tmp_limb |= (b2[preinc(i2)] << neg_bit_offset)
    return (b1h ^ tmp_limb) & mask == 0

cdef inline bint bitset_isempty(fused_bitset_t bits) noexcept nogil:
    """
    Test whether bits is empty.  Return ``True`` (i.e., 1) if the set is
    empty, ``False`` (i.e., 0) otherwise.
    """
    return _bitset_isempty(bits.bits, bits.limbs)

cdef inline bint bitset_is_zero(fused_bitset_t bits) noexcept:
    """
    Test whether bits is empty (i.e., zero).  Return ``True`` (1) if
    the set is empty, ``False`` (0) otherwise.

    This function is the same as bitset_is_empty(bits).
    """
    return bitset_isempty(bits)

cdef inline bint bitset_eq(fused_bitset_t a, fused_bitset_t b) noexcept:
    """
    Compare bitset a and b.  Return ``True`` (i.e., 1) if the sets are
    equal, and ``False`` (i.e., 0) otherwise.

    We assume ``a.limbs >= b.limbs``.
    """
    return _bitset_cmp(a.bits, b.bits, b.limbs, EQUAL)

cdef inline int bitset_cmp(fused_bitset_t a, fused_bitset_t b) noexcept:
    """
    Compare bitsets a and b.  Return 0 if the two sets are
    identical, and consistently return -1 or 1 for two sets that are
    not equal.

    We assume ``a.limbs >= b.limbs``.
    """
    return mpn_cmp(a.bits, b.bits, b.limbs)

cdef inline int bitset_lex_cmp(fused_bitset_t a, fused_bitset_t b) noexcept:
    """
    Compare bitsets ``a`` and ``b`` using lexicographical ordering.

    In this order `a < b` if, for some `k`, the first `k` elements from
    `[0 ... n-1]` are in `a` if and only if they are in `b`, and the
    `(k+1)`st element is in `b` but not ``a``. So `1010 < 1011` and
    `1010111 < 1011000`.

    We assume ``a.limbs == b.limbs``.

    INPUT:

    - ``a`` -- a bitset
    - ``b`` -- a bitset, assumed to have the same size as ``a``

    OUTPUT:

    Return ``0`` if the two sets are identical, return ``1`` if ``a > b``,
    and return ``-1`` if ``a < b``.
    """
    cdef long i = bitset_first_diff(a, b)
    if i == -1:
        return 0
    if bitset_in(a, i):
        return 1
    else:
        return -1

cdef inline bint bitset_issubset(fused_bitset_t a, fused_bitset_t b) noexcept nogil:
    """
    Test whether a is a subset of b (i.e., every element in a is also
    in b).

    We assume ``a.limbs <= b.limbs``.
    """
    if fused_bitset_t is sparse_bitset_t and a.non_zero_chunks_are_initialized:
        return _sparse_bitset_cmp(a.bits, a.non_zero_chunks, a.n_non_zero_chunks, b.bits, SUBSET)
    else:
        return _bitset_cmp(a.bits, b.bits, a.limbs, SUBSET)

cdef inline bint bitset_issuperset(fused_bitset_t a, fused_bitset_t b) noexcept nogil:
    """
    Test whether a is a superset of b (i.e., every element in b is also
    in a).

    We assume ``a.limbs >= b.limbs``.
    """
    return bitset_issubset(b, a)

cdef inline bint bitset_are_disjoint(fused_bitset_t a, fused_bitset_t b) noexcept:
    """
    Test whether ``a`` and ``b`` have an empty intersection.

    We assume ``a.limbs <= b.limbs``.
    """
    if fused_bitset_t is sparse_bitset_t and a.non_zero_chunks_are_initialized:
        return _sparse_bitset_cmp(a.bits, a.non_zero_chunks, a.n_non_zero_chunks, b.bits, DISJOINT)
    else:
        return _bitset_cmp(a.bits, b.bits, a.limbs, DISJOINT)


#############################################################################
# Bitset Bit Manipulation
#############################################################################

cdef inline bint bitset_in(fused_bitset_t bits, mp_bitcnt_t n) noexcept:
    """
    Check if n is in bits.  Return ``True`` (i.e., 1) if n is in the
    set, ``False`` (i.e., 0) otherwise.
    """
    return (bits.bits[n >> index_shift] >> (n % GMP_LIMB_BITS)) & 1

cdef inline bint bitset_check(fused_bitset_t bits, mp_bitcnt_t n) noexcept:
    """
    Check if n is in bits.  Return ``True`` (i.e., 1) if n is in the
    set, ``False`` (i.e., 0) otherwise.

    This function is the same as bitset_in(bits, n).
    """
    return bitset_in(bits, n)

cdef inline bint bitset_not_in(fused_bitset_t bits, mp_bitcnt_t n) noexcept:
    """
    Check if n is not in bits.  Return ``True`` (i.e., 1) if n is not in the
    set, ``False`` (i.e., 0) otherwise.
    """
    return not bitset_in(bits, n)

cdef inline bint bitset_remove(fused_bitset_t bits, mp_bitcnt_t n) except -1:
    """
    Remove ``n`` from ``bits``.

    This raises a :exc:`KeyError` if ``n`` is not contained in ``bits``.
    """
    if not bitset_in(bits, n):
        raise KeyError(n)
    bitset_discard(bits, n)
    if fused_bitset_t is sparse_bitset_t:
        bits.non_zero_chunks_are_initialized = False

cdef inline void bitset_discard(fused_bitset_t bits, mp_bitcnt_t n) noexcept:
    """
    Remove n from bits.
    """
    bits.bits[n >> index_shift] &= limb_one_zero_bit(n)
    if fused_bitset_t is sparse_bitset_t:
        bits.non_zero_chunks_are_initialized = False

cdef inline void bitset_unset(fused_bitset_t bits, mp_bitcnt_t n) noexcept:
    """
    Remove n from bits.

    This function is the same as bitset_discard(bits, n).
    """
    bitset_discard(bits, n)

cdef inline void bitset_add(fused_bitset_t bits, mp_bitcnt_t n) noexcept:
    """
    Add n to bits.
    """
    bits.bits[n >> index_shift] |= limb_one_set_bit(n)
    if fused_bitset_t is sparse_bitset_t:
        bits.non_zero_chunks_are_initialized = False

cdef inline void bitset_set(fused_bitset_t bits, mp_bitcnt_t n) noexcept:
    """
    Add n to bits.

    This function is the same as bitset_add(bits, n).
    """
    bitset_add(bits, n)

cdef inline void bitset_set_to(bitset_t bits, mp_bitcnt_t n, bint b) noexcept:
    """
    If b is True, add n to bits.  If b is False, remove n from bits.
    """
    bitset_unset(bits, n)
    bits.bits[n >> index_shift] |= (<mp_limb_t>b) << (n % GMP_LIMB_BITS)

cdef inline void bitset_flip(fused_bitset_t bits, mp_bitcnt_t n) noexcept:
    """
    If n is in bits, remove n from bits.  If n is not in bits, add n
    to bits.
    """
    bits.bits[n >> index_shift] ^= limb_one_set_bit(n)
    if fused_bitset_t is sparse_bitset_t:
        bits.non_zero_chunks_are_initialized = False

cdef inline void bitset_set_first_n(fused_bitset_t bits, mp_bitcnt_t n) noexcept:
    """
    Set exactly the first n bits.
    """
    cdef mp_size_t i
    cdef mp_size_t index = n >> index_shift
    for i in range(index):
        bits.bits[i] = -1
    if index < bits.limbs:
        bits.bits[index] = limb_lower_bits_down(n)
    for i in range(index+1, bits.limbs):
        bits.bits[i] = 0
    if fused_bitset_t is sparse_bitset_t:
        bits.non_zero_chunks_are_initialized = False

#############################################################################
# Bitset Searching
#############################################################################

cdef inline long bitset_first(fused_bitset_t a) noexcept:
    """
    Calculate the index of the first element in the set. If the set
    is empty, returns -1.
    """
    cdef mp_size_t i
    for i in range(a.limbs):
        if a.bits[i]:
            return (i << index_shift) | _bitset_first_in_limb_nonzero(a.bits[i])
    return -1

cdef inline long bitset_first_in_complement(fused_bitset_t a) noexcept:
    """
    Calculate the index of the first element not in the set. If the set
    is full, returns -1.
    """
    cdef mp_size_t i
    cdef mp_bitcnt_t j
    for i in range(a.limbs):
        if ~a.bits[i]:
            j = (i << index_shift) | _bitset_first_in_limb_nonzero(~a.bits[i])
            if j >= a.size:
                return -1
            return <mp_size_t>j
    return -1

cdef inline long bitset_pop(fused_bitset_t a) except -1:
    """
    Remove and return an arbitrary element from the set.

    This raises a :exc:`KeyError` if the set is empty.
    """
    cdef long i = bitset_first(a)
    if i == -1:
        raise KeyError('pop from an empty set')
    bitset_discard(a, i)
    return i

cdef inline long bitset_first_diff(fused_bitset_t a, fused_bitset_t b) noexcept:
    """
    Calculate the index of the first difference between a and b.  If a
    and b are equal, then return -1.

    We assume ``a.limbs == b.limbs``.
    """
    cdef mp_size_t i
    for i in range(a.limbs):
        if a.bits[i] != b.bits[i]:
            return (i << index_shift) | _bitset_first_in_limb_nonzero(a.bits[i] ^ b.bits[i])
    return -1

cdef inline long bitset_next(fused_bitset_t a, mp_bitcnt_t n) noexcept:
    """
    Calculate the index of the next element in the set, starting at
    (and including) n.  Return -1 if there are no elements from n
    onwards.
    """
    if n >= a.size:
        return -1
    cdef mp_size_t i = n >> index_shift
    cdef mp_limb_t limb = a.bits[i] & ~limb_lower_bits_down(n)
    cdef long ret = _bitset_first_in_limb(limb)
    if ret != -1:
        return (i << index_shift) | ret
    for i in range((n >> index_shift) + 1, a.limbs):
        if a.bits[i]:
            return (i << index_shift) | _bitset_first_in_limb_nonzero(a.bits[i])
    return -1

cdef inline long bitset_next_diff(fused_bitset_t a, fused_bitset_t b, mp_bitcnt_t n) noexcept:
    """
    Calculate the index of the next element that differs between a and
    b, starting at (and including) n.  Return -1 if there are no
    elements differing between a and b from n onwards.

    We assume ``a.limbs == b.limbs``.
    """
    if n >= a.size:
        return -1
    cdef mp_size_t i = n >> index_shift
    cdef mp_limb_t limb = (a.bits[i] ^ b.bits[i]) & ~limb_lower_bits_down(n)
    cdef long ret = _bitset_first_in_limb(limb)
    if ret != -1:
        return (i << index_shift) | ret
    for i in range((n >> index_shift) + 1, a.limbs):
        if a.bits[i] != b.bits[i]:
            return (i << index_shift) | _bitset_first_in_limb(a.bits[i] ^ b.bits[i])
    return -1

cdef inline long bitset_len(fused_bitset_t bits) noexcept nogil:
    """
    Calculate the number of items in the set (i.e., the number of nonzero bits).
    """
    return _bitset_len(bits.bits, bits.limbs)

cdef inline long bitset_hash(fused_bitset_t bits) noexcept:
    """
    Calculate a (very naive) hash function.

    This function should not depend on the size of the bitset, only on
    the items in the bitset.
    """
    cdef mp_limb_t hash = 0
    cdef mp_size_t i
    for i in range(bits.limbs):
        hash += bits.bits[i]
    return hash

#############################################################################
# Bitset Arithmetic
#############################################################################

cdef inline void bitset_complement(fused_bitset_t r, fused_bitset_t a) noexcept:
    """
    Set r to be the complement of a, overwriting r.

    We assume ``r.limbs == a.limbs``.
    """
    mpn_com(r.bits, a.bits, a.limbs)
    bitset_fix(r)
    if fused_bitset_t is sparse_bitset_t:
        r.non_zero_chunks_are_initialized = False

cdef inline void bitset_not(fused_bitset_t r, fused_bitset_t a) noexcept:
    """
    Set r to be the complement of a, overwriting r.

    We assume ``r.limbs == a.limbs``.

    This function is the same as bitset_complement(r, a).
    """
    bitset_complement(r, a)
    if fused_bitset_t is sparse_bitset_t:
        r.non_zero_chunks_are_initialized = False

cdef inline void bitset_intersection(fused_bitset_t r, fused_bitset_t a, fused_bitset_t b) noexcept nogil:
    """
    Set r to the intersection of a and b, overwriting r.

    We assume ``a.limbs >= r.limbs == b.limbs``.
    """
    _bitset_operation(r.bits, a.bits, b.bits, b.limbs, AND)
    if fused_bitset_t is sparse_bitset_t:
        r.non_zero_chunks_are_initialized = False

cdef inline void sparse_bitset_intersection(sparse_bitset_t r, fused_bitset_t a, fused_bitset_t b) noexcept nogil:
    """
    Set r to the intersection of a and b, overwriting r.

    Also set the nonzero positions of ``r``.

    We assume ``a.limbs >= r.limbs == b.limbs``.
    """
    r.n_non_zero_chunks = _sparse_bitset_operation(r.bits, r.non_zero_chunks, a.bits, b.bits, b.limbs, AND)
    r.non_zero_chunks_are_initialized = True

cdef inline void bitset_and(fused_bitset_t r, fused_bitset_t a, fused_bitset_t b) noexcept:
    """
    Set r to the intersection of a and b, overwriting r.

    We assume ``a.limbs >= r.limbs == b.limbs``.

    This function is the same as bitset_intersection(r, a, b).
    """
    bitset_intersection(r, a, b)

cdef inline void bitset_union(fused_bitset_t r, fused_bitset_t a, fused_bitset_t b) noexcept nogil:
    """
    Set r to the union of a and b, overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.
    """
    _bitset_operation(r.bits, a.bits, b.bits, b.limbs, OR)
    if fused_bitset_t is sparse_bitset_t:
        r.non_zero_chunks_are_initialized = False

cdef inline void sparse_bitset_union(sparse_bitset_t r, fused_bitset_t a, fused_bitset_t b) noexcept nogil:
    """
    Set r to the union of a and b, overwriting r.

    Also set the nonzero positions of ``r``.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.
    """
    r.n_non_zero_chunks = _sparse_bitset_operation(r.bits, r.non_zero_chunks, a.bits, b.bits, b.limbs, OR)
    r.non_zero_chunks_are_initialized = True

cdef inline void bitset_or(fused_bitset_t r, fused_bitset_t a, fused_bitset_t b) noexcept:
    """
    Set r to the union of a and b, overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.

    This function is the same as bitset_union(r, a, b).
    """
    bitset_union(r, a, b)

cdef inline void bitset_difference(fused_bitset_t r, fused_bitset_t a, fused_bitset_t b) noexcept:
    """
    Set r to the difference of a and b (i.e., things in a that are not
    in b), overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.
    """
    _bitset_operation(r.bits, a.bits, b.bits, b.limbs, ANDNOT)
    if fused_bitset_t is sparse_bitset_t:
        r.non_zero_chunks_are_initialized = False

cdef inline void sparse_bitset_difference(sparse_bitset_t r, fused_bitset_t a, fused_bitset_t b) noexcept:
    """
    Set r to the difference of a and b (i.e., things in a that are not
    in b), overwriting r.

    Also set the nonzero positions of ``r``.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.
    """
    r.n_non_zero_chunks = _sparse_bitset_operation(r.bits, r.non_zero_chunks, a.bits, b.bits, b.limbs, ANDNOT)
    r.non_zero_chunks_are_initialized = True

cdef inline void bitset_symmetric_difference(fused_bitset_t r, fused_bitset_t a, fused_bitset_t b) noexcept:
    """
    Set r to the symmetric difference of a and b, overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.
    """
    _bitset_operation(r.bits, a.bits, b.bits, b.limbs, XOR)
    if fused_bitset_t is sparse_bitset_t:
        r.non_zero_chunks_are_initialized = False

cdef inline void sparse_bitset_symmetric_difference(sparse_bitset_t r, fused_bitset_t a, fused_bitset_t b) noexcept:
    """
    Set r to the symmetric difference of a and b, overwriting r.

    Also set the nonzero positions of ``r``.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.
    """
    r.n_non_zero_chunks = _sparse_bitset_operation(r.bits, r.non_zero_chunks, a.bits, b.bits, b.limbs, XOR)
    r.non_zero_chunks_are_initialized = True

cdef inline void bitset_xor(fused_bitset_t r, fused_bitset_t a, fused_bitset_t b) noexcept:
    """
    Set r to the symmetric difference of a and b, overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.

    This function is the same as bitset_symmetric_difference(r, a, b).
    """
    bitset_symmetric_difference(r, a, b)

cdef inline void bitset_rshift(fused_bitset_t r, fused_bitset_t a, mp_bitcnt_t n) noexcept:
    """
    Shift the bitset ``a`` right by ``n`` bits and store the result in
    ``r``.

    There are no assumptions on the sizes of ``a`` and ``r``.  Bits which are
    shifted outside of the resulting bitset are discarded.
    """
    if n >= a.size:
        mpn_zero(r.bits, r.limbs)
        return

    # Number of limbs on the right of a which will totally be shifted out
    cdef mp_size_t nlimbs = n >> index_shift
    # Number of limbs to be shifted assuming r is large enough
    cdef mp_size_t shifted_limbs = a.limbs - nlimbs
    # Number of bits to shift additionally
    cdef mp_bitcnt_t nbits = n % GMP_LIMB_BITS

    if shifted_limbs < r.limbs:
        if nbits:
            mpn_rshift(r.bits, a.bits + nlimbs, shifted_limbs, nbits)
        else:
            mpn_copyi(r.bits, a.bits + nlimbs, shifted_limbs)

        # Clear top limbs (note that r.limbs - shifted_limbs >= 1)
        mpn_zero(r.bits + (r.limbs - nlimbs), r.limbs - shifted_limbs)
    else:
        # Number of limbs to shift is r.limbs
        if nbits:
            mpn_rshift(r.bits, a.bits + nlimbs, r.limbs, nbits)
            if shifted_limbs > r.limbs:
                # Add the additional bits from top limb of a
                r.bits[r.limbs-1] |= a.bits[r.limbs+nlimbs] << (GMP_LIMB_BITS - nbits)
        else:
            mpn_copyi(r.bits, a.bits + nlimbs, r.limbs)

        # Clear bits outside bitset in top limb
        bitset_fix(r)

    if fused_bitset_t is sparse_bitset_t:
        r.non_zero_chunks_are_initialized = False

cdef inline void bitset_lshift(fused_bitset_t r, fused_bitset_t a, mp_bitcnt_t n) noexcept:
    """
    Shift the bitset ``a`` left by ``n`` bits and store the result in
    ``r``.

    There are no assumptions on the sizes of ``a`` and ``r``.  Bits which are
    shifted outside of the resulting bitset are discarded.
    """
    if n >= r.size:
        mpn_zero(r.bits, r.limbs)
        return

    # Number of limbs on the right of r which will totally be zeroed
    cdef mp_size_t nlimbs = n >> index_shift
    # Number of limbs to be shifted assuming a is large enough
    cdef mp_size_t shifted_limbs = r.limbs - nlimbs
    # Number of bits to shift additionally
    cdef mp_bitcnt_t nbits = n % GMP_LIMB_BITS

    cdef mp_limb_t out = 0
    if shifted_limbs > a.limbs:
        if nbits:
            out = mpn_lshift(r.bits + nlimbs, a.bits, a.limbs, nbits)
        else:
            mpn_copyd(r.bits + nlimbs, a.bits, a.limbs)

        # Clear top limbs (note that shifted_limbs - a.limbs >= 1)
        mpn_zero(r.bits + a.limbs + nlimbs, shifted_limbs - a.limbs)
        # Store extra limb shifted in from a
        r.bits[nlimbs+a.limbs] = out
    else:
        if nbits:
            mpn_lshift(r.bits + nlimbs, a.bits, shifted_limbs, nbits)
        else:
            mpn_copyd(r.bits + nlimbs, a.bits, shifted_limbs)

        # Clear bits outside bitset in top limb
        bitset_fix(r)

    # Clear bottom limbs
    mpn_zero(r.bits, nlimbs)

    if fused_bitset_t is sparse_bitset_t:
        r.non_zero_chunks_are_initialized = False

cdef inline int bitset_map(fused_bitset_t r, fused_bitset_t a, m) except -1:
    """
    Fill bitset ``r`` so ``r == {m[i] for i in a}``.

    We assume ``m`` has a dictionary interface such that
    ``m[i]`` is an integer in ``[0 ... n-1]`` for all ``i`` in ``a``,
    where ``n`` is the capacity of ``r``.
    """
    if fused_bitset_t is sparse_bitset_t:
        r.non_zero_chunks_are_initialized = False
    cdef long i
    bitset_clear(r)
    i = bitset_first(a)
    while i >= 0:
        bitset_add(r, <mp_bitcnt_t> m[i])
        i = bitset_next(a, i + 1)
    return 0

#############################################################################
# Hamming Weights
#############################################################################

cdef inline long bitset_hamming_weight(fused_bitset_t a) noexcept:
    return bitset_len(a)

#############################################################################
# Bitset Conversion
#############################################################################

cdef char* bitset_chars(char* s, fused_bitset_t bits, char zero=*, char one=*) noexcept

cdef int bitset_from_char(bitset_t bits, char* s, char zero=*, char one=*) except -1

cdef int bitset_from_str(bitset_t bits, object s, char zero=*, char one=*) except -1

cdef bitset_string(fused_bitset_t bits)

cdef bitset_bytes(fused_bitset_t bits)

cdef list bitset_list(fused_bitset_t bits)

cdef bitset_pickle(bitset_t bs)

cdef bitset_unpickle(bitset_t bs, tuple input)
