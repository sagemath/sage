"""
A fast bitset datatype in Cython.

Operations between bitsets are only guaranteed to work if the bitsets
have the same size, with the exception of ``bitset_realloc``.  Similarly, you
should not try to access elements of a bitset beyond the size.

AUTHORS:

- Robert Bradshaw (2008)
- Rudi Pendavingh, Stefan van Zwam (2013-06-06): added functions map, lex_cmp,
  pickle, unpickle
- Jeroen Demeyer (2014-09-05): use mpn_* functions from MPIR in the
  implementation (:trac:`13352` and :trac:`16937`)
- Simon King (2014-10-28): ``bitset_rshift`` and ``bitset_lshift`` respecting
  the size of the given bitsets (:trac:`15820`)
"""

#*****************************************************************************
#     Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "cysignals/memory.pxi"
from libc.string cimport strlen
from sage.libs.gmp.mpn cimport *
from sage.data_structures.bitset cimport *
from cython.operator import preincrement as preinc

# Doctests for the functions in this file are in sage/data_structures/bitset.pyx

#############################################################################
# Creating limb patterns
#############################################################################
#
# NOTE: In all functions in this section, the index n is interpreted
# modulo GMP_LIMB_BITS, the number of bits in a limb.
#
cdef inline mp_limb_t limb_one_set_bit(mp_bitcnt_t n):
    """
    Return a limb with only bit n set.
    """
    return (<mp_limb_t>1) << (n % GMP_LIMB_BITS)

cdef inline mp_limb_t limb_one_zero_bit(mp_bitcnt_t n):
    """
    Return a limb with all bits set, except for bit n.
    """
    return ~((<mp_limb_t>1) << (n % GMP_LIMB_BITS))

cdef inline mp_limb_t limb_lower_bits_down(mp_bitcnt_t n):
    """
    Return a limb with the lower n bits set, where n is interpreted
    in [0 .. GMP_LIMB_BITS-1].
    """
    return ((<mp_limb_t>1) << (n % GMP_LIMB_BITS)) - 1

cdef inline mp_limb_t limb_lower_bits_up(mp_bitcnt_t n):
    """
    Return a limb with the lower n bits set, where n is interpreted
    in [1 .. GMP_LIMB_BITS].
    """
    return (<mp_limb_t>(-1)) >> ((<unsigned int>(-n)) % GMP_LIMB_BITS)

#############################################################################
# Bitset Initalization
#############################################################################
cdef inline bint bitset_init(bitset_t bits, mp_bitcnt_t size) except -1:
    """
    Allocate an empty bitset of size ``size``.

    Size must be at least 1.
    """
    if size <= 0:
        raise ValueError("bitset capacity must be greater than 0")

    bits.size = size
    bits.limbs = (size - 1) / (8 * sizeof(mp_limb_t)) + 1
    bits.bits = <mp_limb_t*>check_calloc(bits.limbs, sizeof(mp_limb_t))

cdef inline bint bitset_realloc(bitset_t bits, mp_bitcnt_t size) except -1:
    """
    Reallocate a bitset to size size. If reallocation is larger, new bitset
    does not contain any of the extra bits.
    """
    cdef mp_size_t limbs_old = bits.limbs
    cdef mp_bitcnt_t size_old = bits.size
    if size_old == size:
        return 0
    if size <= 0:
        raise ValueError("bitset capacity must be greater than 0")

    bits.limbs = (size - 1) / (8 * sizeof(mp_limb_t)) + 1
    tmp = <mp_limb_t*>sig_realloc(bits.bits, bits.limbs * sizeof(mp_limb_t))
    if tmp != NULL:
        bits.bits = tmp
    else:
        bits.limbs = limbs_old
        raise MemoryError
    bits.size = size

    if bits.limbs > limbs_old:
        # Zero any extra limbs
        mpn_zero(bits.bits + limbs_old, bits.limbs - limbs_old)
    elif bits.size < size_old:
        # Zero removed bits
        bitset_fix(bits)

cdef inline void bitset_free(bitset_t bits):
    """
    Deallocate the memory in bits.
    """
    sig_free(bits.bits)

cdef inline void bitset_clear(bitset_t bits):
    """
    Remove all elements from the set.
    """
    mpn_zero(bits.bits, bits.limbs)

cdef inline void bitset_zero(bitset_t bits):
    """
    Remove all elements from the set.

    This function is the same as bitset_clear(bits).
    """
    mpn_zero(bits.bits, bits.limbs)

cdef inline void bitset_copy(bitset_t dst, bitset_t src):
    """
    Copy the bitset src over to the bitset dst, overwriting dst.

    We assume ``dst.limbs == src.limbs``.
    """
    mpn_copyi(dst.bits, src.bits, src.limbs)

cdef inline void bitset_fix(bitset_t bits):
    """
    Clear upper bits in upper limb which should be zero.
    """
    bits.bits[bits.limbs - 1] &= limb_lower_bits_up(bits.size)

#############################################################################
# Bitset Comparison
#############################################################################

cdef inline bint mpn_equal_bits(mp_srcptr b1, mp_srcptr b2, mp_bitcnt_t n):
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

cdef bint mpn_equal_bits_shifted(mp_srcptr b1, mp_srcptr b2, mp_bitcnt_t n, mp_bitcnt_t offset):
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
    for i1 from 0 <= i1 < nlimbs:
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

cdef inline bint bitset_isempty(bitset_t bits):
    """
    Test whether bits is empty.  Return True (i.e., 1) if the set is
    empty, False (i.e., 0) otherwise.
    """
    # First check lowest limb
    if bits.bits[0]:
        return False
    if bits.limbs == 1:
        return True
    # Compare bits to itself shifted by 1 limb. If these compare equal,
    # all limbs must be 0.
    return mpn_cmp(bits.bits+1, bits.bits, bits.limbs-1) == 0

cdef inline bint bitset_is_zero(bitset_t bits):
    """
    Test whether bits is empty (i.e., zero).  Return True (1) if
    the set is empty, False (0) otherwise.

    This function is the same as bitset_is_empty(bits).
    """
    return bitset_isempty(bits)

cdef inline bint bitset_eq(bitset_t a, bitset_t b):
    """
    Compare bitset a and b.  Return True (i.e., 1) if the sets are
    equal, and False (i.e., 0) otherwise.

    We assume ``a.limbs >= b.limbs``.
    """
    return mpn_cmp(a.bits, b.bits, b.limbs) == 0

cdef inline int bitset_cmp(bitset_t a, bitset_t b):
    """
    Compare bitsets a and b.  Returns 0 if the two sets are
    identical, and consistently return -1 or 1 for two sets that are
    not equal.

    We assume ``a.limbs >= b.limbs``.
    """
    return mpn_cmp(a.bits, b.bits, b.limbs)

cdef inline int bitset_lex_cmp(bitset_t a, bitset_t b):
    """
    Compare bitsets ``a`` and ``b`` using lexicographical ordering.

    In this order `a < b` if, for some `k`, the first `k` elements from
    `[0 ... n-1]` are in `a` if and only if they are in `b`, and the
    `(k+1)`st element is in `b` but not ``a``. So `1010 < 1011` and
    `1010111 < 1011000`.

    We assume ``a.limbs == b.limbs``.

    INPUT:

    - ``a`` -- a bitset
    - ``b`` -- a bitset, assumed to have the same size as ``a``.

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

cdef inline bint bitset_issubset(bitset_t a, bitset_t b):
    """
    Test whether a is a subset of b (i.e., every element in a is also
    in b).

    We assume ``a.limbs <= b.limbs``.
    """
    cdef mp_size_t i
    for i from 0 <= i < a.limbs:
        if (a.bits[i] & ~b.bits[i]) != 0:
            return False
    return True

cdef inline bint bitset_issuperset(bitset_t a, bitset_t b):
    """
    Test whether a is a superset of b (i.e., every element in b is also
    in a).

    We assume ``a.limbs >= b.limbs``.
    """
    return bitset_issubset(b, a)


cdef inline bint bitset_are_disjoint(bitset_t a, bitset_t b):
    """
    Tests whether ``a`` and ``b`` have an empty intersection.

    We assume ``a.limbs <= b.limbs``.
    """
    cdef mp_size_t i
    for i from 0 <= i < a.limbs:
        if (a.bits[i]&b.bits[i]) != 0:
            return False
    return True


#############################################################################
# Bitset Bit Manipulation
#############################################################################

cdef inline bint bitset_in(bitset_t bits, mp_bitcnt_t n):
    """
    Check if n is in bits.  Return True (i.e., 1) if n is in the
    set, False (i.e., 0) otherwise.
    """
    return (bits.bits[n >> index_shift] >> (n % GMP_LIMB_BITS)) & 1

cdef inline bint bitset_check(bitset_t bits, mp_bitcnt_t n):
    """
    Check if n is in bits.  Return True (i.e., 1) if n is in the
    set, False (i.e., 0) otherwise.

    This function is the same as bitset_in(bits, n).
    """
    return bitset_in(bits, n)

cdef inline bint bitset_not_in(bitset_t bits, mp_bitcnt_t n):
    """
    Check if n is not in bits.  Return True (i.e., 1) if n is not in the
    set, False (i.e., 0) otherwise.
    """
    return not bitset_in(bits, n)

cdef inline bint bitset_remove(bitset_t bits, mp_bitcnt_t n) except -1:
    """
    Remove n from bits.  Raise KeyError if n is not contained in bits.
    """
    if not bitset_in(bits, n):
        raise KeyError(n)
    bitset_discard(bits, n)

cdef inline void bitset_discard(bitset_t bits, mp_bitcnt_t n):
    """
    Remove n from bits.
    """
    bits.bits[n >> index_shift] &= limb_one_zero_bit(n)

cdef inline void bitset_unset(bitset_t bits, mp_bitcnt_t n):
    """
    Remove n from bits.

    This function is the same as bitset_discard(bits, n).
    """
    bitset_discard(bits, n)


cdef inline void bitset_add(bitset_t bits, mp_bitcnt_t n):
    """
    Add n to bits.
    """
    bits.bits[n >> index_shift] |= limb_one_set_bit(n)

cdef inline void bitset_set(bitset_t bits, mp_bitcnt_t n):
    """
    Add n to bits.

    This function is the same as bitset_add(bits, n).
    """
    bitset_add(bits, n)

cdef inline void bitset_set_to(bitset_t bits, mp_bitcnt_t n, bint b):
    """
    If b is True, add n to bits.  If b is False, remove n from bits.
    """
    bitset_unset(bits, n)
    bits.bits[n >> index_shift] |= (<mp_limb_t>b) << (n % GMP_LIMB_BITS)

cdef inline void bitset_flip(bitset_t bits, mp_bitcnt_t n):
    """
    If n is in bits, remove n from bits.  If n is not in bits, add n
    to bits.
    """
    bits.bits[n >> index_shift] ^= limb_one_set_bit(n)

cdef inline void bitset_set_first_n(bitset_t bits, mp_bitcnt_t n):
    """
    Set exactly the first n bits.
    """
    cdef mp_size_t i
    cdef mp_size_t index = n >> index_shift
    for i from 0 <= i < index:
        bits.bits[i] = -1
    if index < bits.limbs:
        bits.bits[index] = limb_lower_bits_down(n)
    for i from index < i < bits.limbs:
        bits.bits[i] = 0

#############################################################################
# Bitset Searching
#############################################################################

cdef inline long _bitset_first_in_limb_nonzero(mp_limb_t limb):
    """
    Given a non-zero limb of a bitset, return the index of the first
    nonzero bit.
    """
    return mpn_scan1(&limb, 0)

cdef inline long _bitset_first_in_limb(mp_limb_t limb):
    """
    Given a limb of a bitset, return the index of the first nonzero
    bit. If there are no bits set in the limb, return -1.
    """
    if limb == 0:
        return -1
    return mpn_scan1(&limb, 0)

cdef inline long bitset_first(bitset_t a):
    """
    Calculate the index of the first element in the set. If the set
    is empty, returns -1.
    """
    cdef mp_size_t i
    for i from 0 <= i < a.limbs:
        if a.bits[i]:
            return (i << index_shift) | _bitset_first_in_limb_nonzero(a.bits[i])
    return -1

cdef inline long bitset_first_in_complement(bitset_t a):
    """
    Calculate the index of the first element not in the set. If the set
    is full, returns -1.
    """
    cdef mp_size_t i, j
    for i from 0 <= i < a.limbs:
        if ~a.bits[i]:
            j = (i << index_shift) | _bitset_first_in_limb_nonzero(~a.bits[i])
            if j >= a.size:
                j = -1
            return j
    return -1

cdef inline long bitset_pop(bitset_t a) except -1:
    """
    Remove and return an arbitrary element from the set. Raise
    KeyError if the set is empty.
    """
    cdef long i = bitset_first(a)
    if i == -1:
        raise KeyError('pop from an empty set')
    bitset_discard(a, i)
    return i

cdef inline long bitset_first_diff(bitset_t a, bitset_t b):
    """
    Calculate the index of the first difference between a and b.  If a
    and b are equal, then return -1.

    We assume ``a.limbs == b.limbs``.
    """
    cdef mp_size_t i
    for i from 0 <= i < a.limbs:
        if a.bits[i] != b.bits[i]:
            return (i << index_shift) | _bitset_first_in_limb_nonzero(a.bits[i] ^ b.bits[i])
    return -1

cdef inline long bitset_next(bitset_t a, mp_bitcnt_t n):
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
    for i from (n >> index_shift) < i < a.limbs:
        if a.bits[i]:
            return (i << index_shift) | _bitset_first_in_limb_nonzero(a.bits[i])
    return -1

cdef inline long bitset_next_diff(bitset_t a, bitset_t b, mp_bitcnt_t n):
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
    for i from (n >> index_shift) < i < a.limbs:
        if a.bits[i] != b.bits[i]:
            return (i << index_shift) | _bitset_first_in_limb(a.bits[i] ^ b.bits[i])
    return -1

cdef inline long bitset_len(bitset_t bits):
    """
    Calculate the number of items in the set (i.e., the number of nonzero bits).
    """
    return mpn_popcount(bits.bits, bits.limbs)

cdef inline long bitset_hash(bitset_t bits):
    """
    Calculate a (very naive) hash function.

    This function should not depend on the size of the bitset, only on
    the items in the bitset.
    """
    cdef mp_limb_t hash = 0
    cdef mp_size_t i
    for i from 0 <= i < bits.limbs:
        hash += bits.bits[i]
    return hash

#############################################################################
# Bitset Arithmetic
#############################################################################

cdef inline void bitset_complement(bitset_t r, bitset_t a):
    """
    Set r to be the complement of a, overwriting r.

    We assume ``r.limbs == a.limbs``.
    """
    mpn_com(r.bits, a.bits, a.limbs)
    bitset_fix(r)

cdef inline void bitset_not(bitset_t r, bitset_t a):
    """
    Set r to be the complement of a, overwriting r.

    We assume ``r.limbs == a.limbs``.

    This function is the same as bitset_complement(r, a).
    """
    bitset_complement(r, a)

cdef inline void bitset_intersection(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the intersection of a and b, overwriting r.

    We assume ``a.limbs >= r.limbs == b.limbs``.
    """
    mpn_and_n(r.bits, a.bits, b.bits, b.limbs)

cdef inline void bitset_and(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the intersection of a and b, overwriting r.

    We assume ``a.limbs >= r.limbs == b.limbs``.

    This function is the same as bitset_intersection(r, a, b).
    """
    mpn_and_n(r.bits, a.bits, b.bits, b.limbs)

cdef inline void bitset_union(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the union of a and b, overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.
    """
    mpn_ior_n(r.bits, a.bits, b.bits, b.limbs)

cdef inline void bitset_or(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the union of a and b, overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.

    This function is the same as bitset_union(r, a, b).
    """
    mpn_ior_n(r.bits, a.bits, b.bits, b.limbs)

cdef inline void bitset_difference(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the difference of a and b (i.e., things in a that are not
    in b), overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.
    """
    mpn_andn_n(r.bits, a.bits, b.bits, b.limbs)

cdef inline void bitset_symmetric_difference(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the symmetric difference of a and b, overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.
    """
    mpn_xor_n(r.bits, a.bits, b.bits, b.limbs)

cdef inline void bitset_xor(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the symmetric difference of a and b, overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.

    This function is the same as bitset_symmetric_difference(r, a, b).
    """
    mpn_xor_n(r.bits, a.bits, b.bits, b.limbs)


cdef void bitset_rshift(bitset_t r, bitset_t a, mp_bitcnt_t n):
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

cdef void bitset_lshift(bitset_t r, bitset_t a, mp_bitcnt_t n):
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


cdef int bitset_map(bitset_t r, bitset_t a, m) except -1:
    """
    Fill bitset ``r`` so ``r == {m[i] for i in a}``.

    We assume ``m`` has a dictionary interface such that
    ``m[i]`` is an integer in ``[0 ... n-1]`` for all ``i`` in ``a``,
    where ``n`` is the capacity of ``r``.
    """
    cdef long i
    bitset_clear(r)
    i = bitset_first(a)
    while i >= 0:
        bitset_add(r, m[i])
        i = bitset_next(a, i + 1)
    return 0

#############################################################################
# Hamming Weights
#############################################################################

cdef inline long bitset_hamming_weight(bitset_t a):
    return bitset_len(a)

#############################################################################
# Bitset Conversion
#############################################################################

cdef char* bitset_chars(char* s, bitset_t bits, char zero=c'0', char one=c'1'):
    """
    Return a string representation of the bitset in s, using zero for
    the character representing the items not in the bitset and one for
    the character representing the items in the bitset.

    The string is both stored in s and returned.  If s is NULL, then a
    new string is allocated.
    """
    cdef long i
    if s == NULL:
        s = <char *>sig_malloc(bits.size + 1)
    for i from 0 <= i < bits.size:
        s[i] = one if bitset_in(bits, i) else zero
    s[bits.size] = 0
    return s

cdef int bitset_from_str(bitset_t bits, char* s, char zero=c'0', char one=c'1') except -1:
    """
    Initialize a bitset with a set derived from the character string
    s, where one represents the character indicating set membership.
    """
    bitset_init(bits, strlen(s))
    cdef long i
    for i from 0 <= i < bits.size:
        bitset_set_to(bits, i, s[i] == one)
    return 0

cdef bitset_string(bitset_t bits):
    """
    Return a python string representing the bitset.
    """
    cdef char* s = bitset_chars(NULL, bits)
    cdef object py_s
    py_s = s
    sig_free(s)
    return py_s

cdef list bitset_list(bitset_t bits):
    """
    Return a list of elements in the bitset.
    """
    cdef list elts = []
    cdef long elt = bitset_first(bits)
    while elt >= 0:
        elts.append(elt)
        elt = bitset_next(bits, elt + 1)
    return elts

cdef bitset_pickle(bitset_t bs):
    """
    Convert ``bs`` to a reasonably compact Python structure.

    Useful for pickling objects using bitsets as internal data structure.
    To ensure this works on 32-bit and 64-bit machines, the size of a long
    is stored too.
    """
    version = 0
    data = []
    for i from 0 <= i < bs.limbs:
        data.append(bs.bits[i])
    return (version, bs.size, bs.limbs, sizeof(unsigned long), tuple(data))

cdef bitset_unpickle(bitset_t bs, tuple input):
    """
    Convert the data into a bitset.

    Companion of ``bitset_pickle()``. Assumption: ``bs`` has been initialized.
    """
    version, size, limbs, longsize, data = input
    if version != 0:
        raise TypeError("bitset was saved with newer version of Sage. Please upgrade.")
    if bs.size != size:
        bitset_realloc(bs, size)
    if sizeof(unsigned long) == longsize and bs.limbs == limbs:
        for i from 0 <= i < bs.limbs:
            bs.bits[i] = data[i]
    else:
        storage = 8 * longsize  # number of elements encoded in one limb
        adder = 0
        bitset_clear(bs)
        for i from 0 <= i < limbs:
            for j from 0 <= j < storage:
                if (data[i] >> j) & 1:
                    bitset_add(bs, j + adder)
            adder += storage
