r"""
Sequences of bounded integers

This module provides :class:`BoundedIntegerSequence`, which implements
sequences of bounded integers and is for many (but not all) operations faster
than representing the same sequence as a Python :class:`tuple`.

The underlying data structure is similar to :class:`~sage.misc.bitset.Bitset`,
which means that certain operations are implemented by using fast shift
operations from MPIR.  The following boilerplate functions can be
cimported in Cython modules:

- ``cdef bint biseq_init(biseq_t R, mp_size_t l, mp_size_t itemsize) except -1``

  Allocate memory for a bounded integer sequence of length ``l`` with
  items fitting in ``itemsize`` bits.

- ``cdef inline void biseq_dealloc(biseq_t S)``

  Deallocate the memory used by ``S``.

- ``cdef bint biseq_init_copy(biseq_t R, biseq_t S)``

  Initialize ``R`` as a copy of ``S``.

- ``cdef tuple biseq_pickle(biseq_t S)``

  Return a triple ``(bitset_data, itembitsize, length)`` defining ``S``.

- ``cdef bint biseq_unpickle(biseq_t R, tuple bitset_data, mp_bitcnt_t itembitsize, mp_size_t length) except -1``

  Initialise ``R`` from data returned by ``biseq_pickle``.

- ``cdef bint biseq_init_list(biseq_t R, list data, size_t bound) except -1``

  Convert a list to a bounded integer sequence, which must not be allocated.

- ``cdef inline Py_hash_t biseq_hash(biseq_t S)``

  Hash value for ``S``.

- ``cdef inline bint biseq_richcmp(biseq_t S1, biseq_t S2, int op)``

  Comparison of ``S1`` and ``S2``. This takes into account the bound, the
  length, and the list of items of the two sequences.

- ``cdef bint biseq_init_concat(biseq_t R, biseq_t S1, biseq_t S2) except -1``

  Concatenate ``S1`` and ``S2`` and write the result to ``R``. Does not test
  whether the sequences have the same bound!

- ``cdef inline bint biseq_startswith(biseq_t S1, biseq_t S2)``

  Is ``S1=S2+something``? Does not check whether the sequences have the same
  bound!

- ``cdef mp_size_t biseq_contains(biseq_t S1, biseq_t S2, mp_size_t start) except -2``

  Return the position in ``S1`` of ``S2`` as a subsequence of
  ``S1[start:]``, or ``-1`` if ``S2`` is not a subsequence. Does not check
  whether the sequences have the same bound!

- ``cdef mp_size_t biseq_starswith_tail(biseq_t S1, biseq_t S2, mp_size_t start) except -2:``

  Return the smallest number ``i`` such that the bounded integer sequence
  ``S1`` starts with the sequence ``S2[i:]``, where ``start <= i <
  S1.length``, or return ``-1`` if no such ``i`` exists.

- ``cdef mp_size_t biseq_index(biseq_t S, size_t item, mp_size_t start) except -2``

  Return the position in ``S`` of the item in ``S[start:]``, or ``-1`` if
  ``S[start:]`` does not contain the item.

- ``cdef size_t biseq_getitem(biseq_t S, mp_size_t index)``

  Return ``S[index]``, without checking margins.

- ``cdef size_t biseq_getitem_py(biseq_t S, mp_size_t index)``

  Return ``S[index]`` as Python ``int``, without checking margins.

- ``cdef biseq_inititem(biseq_t S, mp_size_t index, size_t item)``

  Set ``S[index] = item``, without checking margins and assuming that ``S[index]``
  has previously been zero.

- ``cdef inline void biseq_clearitem(biseq_t S, mp_size_t index)``

  Set ``S[index] = 0``, without checking margins.

- ``cdef bint biseq_init_slice(biseq_t R, biseq_t S, mp_size_t start, mp_size_t stop, mp_size_t step) except -1``

  Initialise ``R`` with ``S[start:stop:step]``.

AUTHORS:

- Simon King, Jeroen Demeyer (2014-10): initial version (:issue:`15820`)
"""
# ****************************************************************************
#       Copyright (C) 2014 Simon King <simon.king@uni-jena.de>
#       Copyright (C) 2014 Jeroen Demeyer <jdemeyer@cage.ugennt.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_check, sig_on, sig_off
from sage.data_structures.bitset_base cimport *

from cpython.long cimport PyLong_FromSize_t
from cpython.slice cimport PySlice_GetIndicesEx
from sage.libs.flint.flint cimport FLINT_BIT_COUNT as BIT_COUNT
from sage.structure.richcmp cimport richcmp_not_equal, rich_to_bool

cimport cython

###################
#  Boilerplate
#   cdef functions
###################

#
# (De)allocation, copying
#


@cython.overflowcheck
cdef bint biseq_init(biseq_t R, mp_size_t l, mp_bitcnt_t itemsize) except -1:
    """
    Allocate memory for a bounded integer sequence of length ``l`` with
    items fitting in ``itemsize`` bits.
    """
    cdef mp_bitcnt_t totalbitsize
    if l:
        totalbitsize = l * itemsize
    else:
        totalbitsize = 1
    bitset_init(R.data, totalbitsize)
    R.length = l
    R.itembitsize = itemsize
    R.mask_item = limb_lower_bits_up(itemsize)

cdef inline void biseq_dealloc(biseq_t S) noexcept:
    """
    Deallocate the memory used by ``S``.
    """
    bitset_free(S.data)

cdef bint biseq_init_copy(biseq_t R, biseq_t S) except -1:
    """
    Initialize ``R`` as a copy of ``S``.
    """
    biseq_init(R, S.length, S.itembitsize)
    sig_on()
    bitset_copy(R.data, S.data)
    sig_off()

#
# Pickling
#

cdef tuple biseq_pickle(biseq_t S):
    return (bitset_pickle(S.data), S.itembitsize, S.length)

cdef bint biseq_unpickle(biseq_t R, tuple bitset_data, mp_bitcnt_t itembitsize, mp_size_t length) except -1:
    biseq_init(R, length, itembitsize)
    sig_on()
    bitset_unpickle(R.data, bitset_data)
    sig_off()
    return 1

#
# Conversion
#

cdef bint biseq_init_list(biseq_t R, list data, size_t bound) except -1:
    """
    Convert a list into a bounded integer sequence and write the result
    into ``R``, which must not be initialised.

    INPUT:

    - ``data`` -- list of integers

    - ``bound`` -- a number which is the maximal value of an item
    """
    cdef mp_size_t index = 0
    cdef size_t item_c

    biseq_init(R, len(data), BIT_COUNT(bound|<size_t>1))

    for item in data:
        sig_check()
        item_c = item
        if item_c > bound:
            raise OverflowError("list item {!r} larger than {}".format(item, bound) )
        biseq_inititem(R, index, item_c)
        index += 1

cdef inline Py_hash_t biseq_hash(biseq_t S) noexcept:
    return S.itembitsize*(<Py_hash_t>1073807360)+bitset_hash(S.data)

cdef inline bint biseq_richcmp(biseq_t S1, biseq_t S2, int op) noexcept:
    if S1.itembitsize != S2.itembitsize:
        return richcmp_not_equal(S1.itembitsize, S2.itembitsize, op)
    if S1.length != S2.length:
        return richcmp_not_equal(S1.length, S2.length, op)
    return rich_to_bool(op, bitset_cmp(S1.data, S2.data))

#
# Arithmetics
#

cdef bint biseq_init_concat(biseq_t R, biseq_t S1, biseq_t S2) except -1:
    """
    Concatenate two bounded integer sequences ``S1`` and ``S2``.

    ASSUMPTION:

    - The two sequences must have equivalent bounds, i.e., the items on the
      sequences must fit into the same number of bits.

    OUTPUT:

    The result is written into ``R``, which must not be initialised
    """
    biseq_init(R, S1.length + S2.length, S1.itembitsize)
    sig_on()
    bitset_lshift(R.data, S2.data, S1.length * S1.itembitsize)
    bitset_or(R.data, R.data, S1.data)
    sig_off()


cdef inline bint biseq_startswith(biseq_t S1, biseq_t S2) except -1:
    """
    Test if bounded integer sequence ``S1`` starts with bounded integer
    sequence ``S2``.

    ASSUMPTION:

    - The two sequences must have equivalent bounds, i.e., the items on the
      sequences must fit into the same number of bits. This condition is not
      tested.
    """
    if S2.length > S1.length:
        return False
    if S2.length == 0:
        return True
    sig_on()
    ret = mpn_equal_bits(S1.data.bits, S2.data.bits, S2.data.size)
    sig_off()
    return ret


cdef mp_size_t biseq_index(biseq_t S, size_t item, mp_size_t start) except -2:
    """
    Return the position in ``S`` of an item in ``S[start:]``, or -1 if
    ``S[start:]`` does not contain the item.
    """
    cdef mp_size_t index
    sig_on()
    for index from start <= index < S.length:
        if biseq_getitem(S, index) == item:
            sig_off()
            return index
    sig_off()
    return -1


cdef inline size_t biseq_getitem(biseq_t S, mp_size_t index) noexcept:
    """
    Get item ``S[index]``, without checking margins.
    """
    cdef mp_bitcnt_t limb_index, bit_index
    bit_index = (<mp_bitcnt_t>index) * S.itembitsize
    limb_index = bit_index // GMP_LIMB_BITS
    bit_index %= GMP_LIMB_BITS

    cdef mp_limb_t out
    out = (S.data.bits[limb_index]) >> bit_index
    if bit_index + S.itembitsize > GMP_LIMB_BITS:
        # Our item is stored using 2 limbs, add the part from the upper limb
        out |= (S.data.bits[limb_index+1]) << (GMP_LIMB_BITS - bit_index)
    return out & S.mask_item

cdef biseq_getitem_py(biseq_t S, mp_size_t index):
    """
    Get item ``S[index]`` as a Python ``int``, without
    checking margins.
    """
    cdef size_t out = biseq_getitem(S, index)
    return PyLong_FromSize_t(out)

cdef inline void biseq_inititem(biseq_t S, mp_size_t index, size_t item) noexcept:
    """
    Set ``S[index] = item``, without checking margins.

    Note that it is assumed that ``S[index] == 0`` before the assignment.
    """
    cdef mp_bitcnt_t limb_index, bit_index
    bit_index = (<mp_bitcnt_t>index) * S.itembitsize
    limb_index = bit_index // GMP_LIMB_BITS
    bit_index %= GMP_LIMB_BITS

    S.data.bits[limb_index] |= (item << bit_index)
    # Have some bits been shifted out of bound?
    if bit_index + S.itembitsize > GMP_LIMB_BITS:
        # Our item is stored using 2 limbs, add the part from the upper limb
        S.data.bits[limb_index+1] |= (item >> (GMP_LIMB_BITS - bit_index))

cdef inline void biseq_clearitem(biseq_t S, mp_size_t index) noexcept:
    """
    Set ``S[index] = 0``, without checking margins.

    In contrast to ``biseq_inititem``, the previous content of ``S[index]``
    will be erased.
    """
    cdef mp_bitcnt_t limb_index, bit_index
    bit_index = (<mp_bitcnt_t>index) * S.itembitsize
    limb_index = bit_index // GMP_LIMB_BITS
    bit_index %= GMP_LIMB_BITS

    S.data.bits[limb_index] &= ~(S.mask_item << bit_index)
    # Have some bits been shifted out of bound?
    if bit_index + S.itembitsize > GMP_LIMB_BITS:
        # Our item is stored using 2 limbs, add the part from the upper limb
        S.data.bits[limb_index+1] &= ~(S.mask_item >> (GMP_LIMB_BITS - bit_index))

cdef bint biseq_init_slice(biseq_t R, biseq_t S, mp_size_t start, mp_size_t stop, mp_size_t step) except -1:
    """
    Create the slice ``S[start:stop:step]`` as bounded integer sequence
    and write the result to ``R``, which must not be initialised.
    """
    cdef mp_size_t length = 0
    if step > 0:
        if stop > start:
            length = ((stop-start-1)//step)+1
    else:
        if stop < start:
            length = ((stop-start+1)//step)+1
    biseq_init(R, length, S.itembitsize)

    if not length:
        return 0

    if step == 1:
        # Slicing essentially boils down to a shift operation.
        sig_on()
        bitset_rshift(R.data, S.data, start*S.itembitsize)
        sig_off()
        return 0

    # In the general case, we move item by item.
    cdef mp_size_t src_index = start
    cdef mp_size_t tgt_index
    sig_on()
    for tgt_index in range(length):
        biseq_inititem(R, tgt_index, biseq_getitem(S, src_index))
        src_index += step
    sig_off()


cdef mp_size_t biseq_contains(biseq_t S1, biseq_t S2, mp_size_t start) except -2:
    """
    Test if the bounded integer sequence ``S1[start:]`` contains a
    sub-sequence ``S2``.

    INPUT:

    - ``S1``, ``S2`` -- two bounded integer sequences
    - ``start`` -- integer; start index

    OUTPUT:

    The smallest index ``i >= start`` such that ``S1[i:]`` starts with
    ``S2``, or ``-1`` if ``S1[start:]`` does not contain ``S2``.

    ASSUMPTION:

    - The two sequences must have equivalent bounds, i.e., the items on the
      sequences must fit into the same number of bits. This condition is not
      tested.
    """
    if S2.length == 0:
        return start
    cdef mp_size_t index
    sig_on()
    for index from start <= index <= S1.length-S2.length:
        if mpn_equal_bits_shifted(S2.data.bits, S1.data.bits,
                S2.length*S2.itembitsize, index*S2.itembitsize):
            sig_off()
            return index
    sig_off()
    return -1

cdef mp_size_t biseq_startswith_tail(biseq_t S1, biseq_t S2, mp_size_t start) except -2:
    """
    Return the smallest index ``i`` such that the bounded integer sequence
    ``S1`` starts with the sequence ``S2[i:]``, where ``start <= i <
    S2.length``.

    INPUT:

    - ``S1``, ``S2`` -- two bounded integer sequences
    - ``start`` -- integer; start index

    OUTPUT:

    The smallest index ``i >= start`` such that ``S1`` starts with ``S2[i:],
    or ``-1`` if no such ``i < S2.length`` exists.

    ASSUMPTION:

    - The two sequences must have equivalent bounds, i.e., the items on the
      sequences must fit into the same number of bits. This condition is not
      tested.
    """
    # Increase start if S1 is too short to contain S2[start:]
    if S1.length < S2.length - start:
        start = S2.length - S1.length
    cdef mp_size_t index
    sig_on()
    for index from start <= index < S2.length:
        if mpn_equal_bits_shifted(S1.data.bits, S2.data.bits,
                (S2.length - index)*S2.itembitsize, index*S2.itembitsize):
            sig_off()
            return index
    sig_off()
    return -1


###########################################
# A cdef class that wraps the above, and
# behaves like a tuple

from sage.rings.integer cimport smallInteger


cdef class BoundedIntegerSequence:
    """
    A sequence of nonnegative uniformly bounded integers.

    INPUT:

    - ``bound`` -- nonnegative integer. When zero, a :exc:`ValueError`
      will be raised. Otherwise, the given bound is replaced by the
      power of two that is at least the given bound.
    - ``data`` -- list of integers

    EXAMPLES:

    We showcase the similarities and differences between bounded integer
    sequences and lists respectively tuples.

    To distinguish from tuples or lists, we use pointed brackets for the
    string representation of bounded integer sequences::

        sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
        sage: S = BoundedIntegerSequence(21, [2, 7, 20]); S
        <2, 7, 20>

    Each bounded integer sequence has a bound that is a power of two, such
    that all its item are less than this bound::

        sage: S.bound()
        32
        sage: BoundedIntegerSequence(16, [2, 7, 20])
        Traceback (most recent call last):
        ...
        OverflowError: list item 20 larger than 15

    Bounded integer sequences are iterable, and we see that we can recover the
    originally given list::

        sage: L = [randint(0,31) for i in range(5000)]
        sage: S = BoundedIntegerSequence(32, L)
        sage: list(L) == L
        True

    Getting items and slicing works in the same way as for lists::

        sage: n = randint(0,4999)
        sage: S[n] == L[n]
        True
        sage: m = randint(0,1000)
        sage: n = randint(3000,4500)
        sage: s = randint(1, 7)
        sage: list(S[m:n:s]) == L[m:n:s]
        True
        sage: list(S[n:m:-s]) == L[n:m:-s]
        True

    The :meth:`index` method works different for bounded integer sequences and
    tuples or lists. If one asks for the index of an item, the behaviour is
    the same. But we can also ask for the index of a sub-sequence::

        sage: L.index(L[200]) == S.index(L[200])
        True
        sage: S.index(S[100:2000])    # random
        100

    Similarly, containment tests work for both items and sub-sequences::

        sage: S[200] in S
        True
        sage: S[200:400] in S
        True
        sage: S[200]+S.bound() in S
        False

    Bounded integer sequences are immutable, and thus copies are
    identical. This is the same for tuples, but of course not for lists::

        sage: T = tuple(S)
        sage: copy(T) is T
        True
        sage: copy(S) is S
        True
        sage: copy(L) is L
        False

    Concatenation works in the same way for lists, tuples and bounded
    integer sequences::

        sage: M = [randint(0,31) for i in range(5000)]
        sage: T = BoundedIntegerSequence(32, M)
        sage: list(S+T)==L+M
        True
        sage: list(T+S)==M+L
        True
        sage: (T+S == S+T) == (M+L == L+M)
        True

    However, comparison works different for lists and bounded integer
    sequences. Bounded integer sequences are first compared by bound, then by
    length, and eventually by *reverse* lexicographical ordering::

        sage: S = BoundedIntegerSequence(21, [4,1,6,2,7,20,9])
        sage: T = BoundedIntegerSequence(51, [4,1,6,2,7,20])
        sage: S < T   # compare by bound, not length
        True
        sage: T < S
        False
        sage: S.bound() < T.bound()
        True
        sage: len(S) > len(T)
        True

    ::

        sage: T = BoundedIntegerSequence(21, [0,0,0,0,0,0,0,0])
        sage: S < T    # compare by length, not lexicographically
        True
        sage: T < S
        False
        sage: list(T) < list(S)
        True
        sage: len(T) > len(S)
        True

    ::

        sage: T = BoundedIntegerSequence(21, [4,1,5,2,8,20,9])
        sage: T > S   # compare by reverse lexicographic ordering...
        True
        sage: S > T
        False
        sage: len(S) == len(T)
        True
        sage: list(S) > list(T) # direct lexicographic ordering is different
        True

    TESTS:

    We test against various corner cases::

        sage: BoundedIntegerSequence(16, [2, 7, -20])
        Traceback (most recent call last):
        ...
        OverflowError: can...t convert negative value to size_t
        sage: BoundedIntegerSequence(1, [0, 0, 0])
        <0, 0, 0>
        sage: BoundedIntegerSequence(1, [0, 1, 0])
        Traceback (most recent call last):
        ...
        OverflowError: list item 1 larger than 0
        sage: BoundedIntegerSequence(0, [0, 1, 0])
        Traceback (most recent call last):
        ...
        ValueError: positive bound expected
        sage: BoundedIntegerSequence(2, [])
        <>
        sage: BoundedIntegerSequence(2, []) == BoundedIntegerSequence(4, []) # The bounds differ
        False
        sage: BoundedIntegerSequence(16, [2, 7, 4])[1:1]
        <>
    """
    def __cinit__(self, *args, **kwds):
        """
        Allocate memory for underlying data.

        INPUT:

        - ``bound`` -- nonnegative integer
        - ``data`` -- ignored

        .. WARNING::

            If ``bound=0`` then no allocation is done.  Hence, this should
            only be done internally, when calling :meth:`__new__` without :meth:`__init__`.

        EXAMPLES::

            sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
            sage: BoundedIntegerSequence(21, [4,1,6,2,7,20,9])  # indirect doctest
            <4, 1, 6, 2, 7, 20, 9>
        """
        # In __init__, we'll raise an error if the bound is 0.
        self.data.data.bits = NULL

    def __dealloc__(self):
        """
        Free the memory from underlying data.

        EXAMPLES::

            sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(21, [4,1,6,2,7,20,9])
            sage: del S     # indirect doctest
        """
        biseq_dealloc(self.data)

    def __init__(self, bound, data):
        """
        INPUT:

        - ``bound`` -- positive integer; the given bound is replaced by
          the next power of two that is greater than the given bound

        - ``data`` -- list of nonnegative integers; all less than
          ``bound``

        EXAMPLES::

            sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
            sage: L = [randint(0,26) for i in range(5000)]
            sage: S = BoundedIntegerSequence(57, L)   # indirect doctest
            sage: list(S) == L
            True
            sage: S = BoundedIntegerSequence(11, [4,1,6,2,7,4,9]); S
            <4, 1, 6, 2, 7, 4, 9>
            sage: S.bound()
            16

        Non-positive bounds or bounds which are too large result in errors::

            sage: BoundedIntegerSequence(-1, L)
            Traceback (most recent call last):
            ...
            ValueError: positive bound expected
            sage: BoundedIntegerSequence(0, L)
            Traceback (most recent call last):
            ...
            ValueError: positive bound expected
            sage: BoundedIntegerSequence(2^64+1, L)
            Traceback (most recent call last):
            ...
            OverflowError: ... int too large to convert...

        We are testing the corner case of the maximal possible bound::

            sage: S = BoundedIntegerSequence(2*(sys.maxsize+1), [8, 8, 26, 18, 18, 8, 22, 4, 17, 22, 22, 7, 12, 4, 1, 7, 21, 7, 10, 10])
            sage: S
            <8, 8, 26, 18, 18, 8, 22, 4, 17, 22, 22, 7, 12, 4, 1, 7, 21, 7, 10, 10>

        Items that are too large::

            sage: BoundedIntegerSequence(100, [2^256])
            Traceback (most recent call last):
            ...
            OverflowError: ... int too large to convert...
            sage: BoundedIntegerSequence(100, [100])
            Traceback (most recent call last):
            ...
            OverflowError: list item 100 larger than 99

        Bounds that are too large::

            sage: BoundedIntegerSequence(2^256, [200])
            Traceback (most recent call last):
            ...
            OverflowError: ... int too large to convert...
        """
        if bound <= 0:
            raise ValueError("positive bound expected")
        biseq_init_list(self.data, data, bound-1)

    def __copy__(self):
        """
        :class:`BoundedIntegerSequence` is immutable, copying returns ``self``.

        EXAMPLES::

            sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(21, [4,1,6,2,7,20,9])
            sage: copy(S) is S
            True
        """
        return self

    def __reduce__(self):
        """
        Pickling of :class:`BoundedIntegerSequence`.

        EXAMPLES::

            sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
            sage: L = [randint(0,26) for i in range(5000)]
            sage: S = BoundedIntegerSequence(32, L)
            sage: loads(dumps(S)) == S    # indirect doctest
            True

        TESTS:

        The discussion at :issue:`15820` explains why the following is a good test::

            sage: X = BoundedIntegerSequence(21, [4,1,6,2,7,2,3])
            sage: S = BoundedIntegerSequence(21, [0,0,0,0,0,0,0])
            sage: loads(dumps(X+S))
            <4, 1, 6, 2, 7, 2, 3, 0, 0, 0, 0, 0, 0, 0>
            sage: loads(dumps(X+S)) == X+S
            True
            sage: T = BoundedIntegerSequence(21, [0,4,0,1,0,6,0,2,0,7,0,2,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
            sage: T[1::2]
            <4, 1, 6, 2, 7, 2, 3, 0, 0, 0, 0, 0, 0, 0>
            sage: T[1::2] == X+S
            True
            sage: loads(dumps(X[1::2])) == X[1::2]
            True
        """
        return NewBISEQ, biseq_pickle(self.data)

    def __len__(self):
        """
        EXAMPLES::

            sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
            sage: L = [randint(0,26) for i in range(5000)]
            sage: S = BoundedIntegerSequence(57, L)   # indirect doctest
            sage: len(S) == len(L)
            True
        """
        return self.data.length

    def __bool__(self):
        """
        A bounded integer sequence is nonzero if and only if its length is nonzero.

        EXAMPLES::

            sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(13, [0,0,0])
            sage: bool(S)
            True
            sage: bool(S[1:1])
            False
        """
        return self.data.length!=0

    def __repr__(self):
        """
        String representation.

        To distinguish it from Python tuples or lists, we use pointed brackets
        as delimiters.

        EXAMPLES::

            sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
            sage: BoundedIntegerSequence(21, [4,1,6,2,7,20,9])   # indirect doctest
            <4, 1, 6, 2, 7, 20, 9>
            sage: BoundedIntegerSequence(21, [0,0]) + BoundedIntegerSequence(21, [0,0])
            <0, 0, 0, 0>
        """
        return "<" + ", ".join(str(x) for x in self) + ">"

    def bound(self):
        """
        Return the bound of this bounded integer sequence.

        All items of this sequence are nonnegative integers less than the
        returned bound. The bound is a power of two.

        EXAMPLES::

            sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(21, [4,1,6,2,7,20,9])
            sage: T = BoundedIntegerSequence(51, [4,1,6,2,7,20,9])
            sage: S.bound()
            32
            sage: T.bound()
            64
        """
        return smallInteger(1) << self.data.itembitsize

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
            sage: L = [randint(0,26) for i in range(5000)]
            sage: S = BoundedIntegerSequence(27, L)
            sage: list(S) == L   # indirect doctest
            True

        TESTS::

            sage: list(BoundedIntegerSequence(1, []))
            []

        The discussion at :issue:`15820` explains why this is a good test::

            sage: S = BoundedIntegerSequence(21, [0,0,0,0,0,0,0])
            sage: X = BoundedIntegerSequence(21, [4,1,6,2,7,2,3])
            sage: list(X)
            [4, 1, 6, 2, 7, 2, 3]
            sage: list(X+S)
            [4, 1, 6, 2, 7, 2, 3, 0, 0, 0, 0, 0, 0, 0]
            sage: list(BoundedIntegerSequence(21, [0,0]) + BoundedIntegerSequence(21, [0,0]))
            [0, 0, 0, 0]
        """
        cdef mp_size_t index
        for index in range(self.data.length):
            yield biseq_getitem_py(self.data, index)

    def __getitem__(self, index):
        """
        Get single items or slices.

        EXAMPLES::

            sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(21, [4,1,6,2,7,20,9])
            sage: S[2]
            6
            sage: S[1::2]
            <1, 2, 20>
            sage: S[-1::-2]
            <9, 7, 6, 4>

        TESTS::

            sage: S = BoundedIntegerSequence(10^8, list(range(9)))
            sage: S[-1]
            8
            sage: S[8]
            8
            sage: S[9]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
            sage: S[-10]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
            sage: S[2^63]
            Traceback (most recent call last):
            ...
            OverflowError: ... int too large to convert to ...

        ::

            sage: S[-1::-2]
            <8, 6, 4, 2, 0>
            sage: S[1::2]
            <1, 3, 5, 7>

        ::

            sage: L = [randint(0,26) for i in range(5000)]
            sage: S = BoundedIntegerSequence(27, L)
            sage: S[1234] == L[1234]
            True
            sage: list(S[100:2000:3]) == L[100:2000:3]
            True
            sage: list(S[3000:10:-7]) == L[3000:10:-7]
            True
            sage: S[:] == S
            True
            sage: S[:] is S
            True

        ::

            sage: S = BoundedIntegerSequence(21, [0,0,0,0,0,0,0])
            sage: X = BoundedIntegerSequence(21, [4,1,6,2,7,2,3])
            sage: (X+S)[6]
            3
            sage: (X+S)[10]
            0
            sage: (X+S)[12:]
            <0, 0>

        ::

            sage: S[2:2] == X[4:2]
            True

        ::

            sage: S = BoundedIntegerSequence(6, [3, 5, 3, 1, 5, 2, 2, 5, 3, 3, 4])
            sage: S[10]
            4

        ::

            sage: B = BoundedIntegerSequence(27, [8, 8, 26, 18, 18, 8, 22, 4, 17, 22, 22, 7, 12, 4, 1, 7, 21, 7, 10, 10])
            sage: B[8:]
            <17, 22, 22, 7, 12, 4, 1, 7, 21, 7, 10, 10>

        ::

            sage: B1 = BoundedIntegerSequence(8, [0,7])
            sage: B2 = BoundedIntegerSequence(8, [2,1,4])
            sage: B1[0:1]+B2
            <0, 2, 1, 4>
        """
        cdef BoundedIntegerSequence out
        cdef Py_ssize_t start, stop, step, slicelength
        if isinstance(index, slice):
            PySlice_GetIndicesEx(index, self.data.length, &start, &stop, &step, &slicelength)
            if start==0 and stop==self.data.length and step==1:
                return self
            out = BoundedIntegerSequence.__new__(BoundedIntegerSequence, 0, None)
            biseq_init_slice(out.data, self.data, start, stop, step)
            return out
        cdef Py_ssize_t ind
        try:
            ind = index
        except TypeError:
            raise TypeError("Sequence index must be integer or slice")
        if ind < 0:
            ind += self.data.length
        if ind < 0 or ind >= self.data.length:
            raise IndexError("index out of range")
        return biseq_getitem_py(self.data, ind)

    def __contains__(self, other):
        """
        Tells whether this bounded integer sequence contains an item or a sub-sequence.

        EXAMPLES::

            sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(21, [4,1,6,2,7,20,9])
            sage: 6 in S
            True
            sage: BoundedIntegerSequence(21, [2, 7, 20]) in S
            True

        The bound of the sequences matters::

            sage: BoundedIntegerSequence(51, [2, 7, 20]) in S
            False

        ::

            sage: 6+S.bound() in S
            False
            sage: S.index(6+S.bound())
            Traceback (most recent call last):
            ...
            ValueError: 38 is not in sequence

        TESTS:

        The discussion at :issue:`15820` explains why the following are good tests::

            sage: X = BoundedIntegerSequence(21, [4,1,6,2,7,2,3])
            sage: S = BoundedIntegerSequence(21, [0,0,0,0,0,0,0])
            sage: loads(dumps(X+S))
            <4, 1, 6, 2, 7, 2, 3, 0, 0, 0, 0, 0, 0, 0>
            sage: loads(dumps(X+S)) == X+S
            True
            sage: T = BoundedIntegerSequence(21, [0,4,0,1,0,6,0,2,0,7,0,2,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
            sage: T[3::2]==(X+S)[1:]
            True
            sage: T[3::2] in X+S
            True
            sage: T = BoundedIntegerSequence(21, [0,4,0,1,0,6,0,2,0,7,0,2,0,3,0,0,0,16,0,0,0,0,0,0,0,0,0,0])
            sage: T[3::2] in (X+S)
            False

        ::

            sage: S1 = BoundedIntegerSequence(4, [1,3])
            sage: S2 = BoundedIntegerSequence(4, [0])
            sage: S2 in S1
            False

        ::

            sage: B = BoundedIntegerSequence(27, [8, 8, 26, 18, 18, 8, 22, 4, 17, 22, 22, 7, 12, 4, 1, 7, 21, 7, 10, 10])
            sage: B.index(B[8:])
            8

        ::

            sage: -1 in B
            False
        """
        if not isinstance(other, BoundedIntegerSequence):
            try:
                return biseq_index(self.data, other, 0) >= 0
            except OverflowError:
                return False
        cdef BoundedIntegerSequence right = other
        if self.data.itembitsize!=right.data.itembitsize:
            return False
        return biseq_contains(self.data, right.data, 0) >= 0

    cpdef list list(self):
        """
        Convert this bounded integer sequence to a list.

        NOTE:

        A conversion to a list is also possible by iterating over the
        sequence.

        EXAMPLES::

            sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
            sage: L = [randint(0,26) for i in range(5000)]
            sage: S = BoundedIntegerSequence(32, L)
            sage: S.list() == list(S) == L
            True

        The discussion at :issue:`15820` explains why the following is a good test::

            sage: (BoundedIntegerSequence(21, [0,0]) + BoundedIntegerSequence(21, [0,0])).list()
            [0, 0, 0, 0]
        """
        cdef mp_size_t i
        return [biseq_getitem_py(self.data, i) for i in range(self.data.length)]

    cpdef bint startswith(self, BoundedIntegerSequence other) noexcept:
        """
        Tells whether ``self`` starts with a given bounded integer sequence

        EXAMPLES::

            sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
            sage: L = [randint(0,26) for i in range(5000)]
            sage: S = BoundedIntegerSequence(27, L)
            sage: L0 = L[:1000]
            sage: T = BoundedIntegerSequence(27, L0)
            sage: S.startswith(T)
            True
            sage: L0[-1] = (L0[-1] + 1) % 27
            sage: T = BoundedIntegerSequence(27, L0)
            sage: S.startswith(T)
            False
            sage: L0[-1] = (L0[-1] - 1) % 27
            sage: L0[0] = (L0[0] + 1) % 27
            sage: T = BoundedIntegerSequence(27, L0)
            sage: S.startswith(T)
            False
            sage: L0[0] = (L0[0] - 1) % 27

        The bounds of the sequences must be compatible, or :meth:`startswith`
        returns ``False``::

            sage: T = BoundedIntegerSequence(51, L0)
            sage: S.startswith(T)
            False
        """
        if self.data.itembitsize != other.data.itembitsize:
            return False
        return biseq_startswith(self.data, other.data)

    def index(self, other):
        """
        The index of a given item or sub-sequence of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(21, [4,1,6,2,6,20,9,0])
            sage: S.index(6)
            2
            sage: S.index(5)
            Traceback (most recent call last):
            ...
            ValueError: 5 is not in sequence
            sage: S.index(BoundedIntegerSequence(21, [6, 2, 6]))
            2
            sage: S.index(BoundedIntegerSequence(21, [6, 2, 7]))
            Traceback (most recent call last):
            ...
            ValueError: not a sub-sequence

        The bound of (sub-)sequences matters::

            sage: S.index(BoundedIntegerSequence(51, [6, 2, 6]))
            Traceback (most recent call last):
            ...
            ValueError: not a sub-sequence
            sage: S.index(0)
            7
            sage: S.index(S.bound())
            Traceback (most recent call last):
            ...
            ValueError: 32 is not in sequence

        TESTS::

            sage: S = BoundedIntegerSequence(10^9, [2, 2, 2, 1, 2, 4, 3, 3, 3, 2, 2, 0])
            sage: S[11]
            0
            sage: S.index(0)
            11

        ::

            sage: S.index(-3)
            Traceback (most recent call last):
            ...
            ValueError: -3 is not in sequence
            sage: S.index(2^100)
            Traceback (most recent call last):
            ...
            ValueError: 1267650600228229401496703205376 is not in sequence
            sage: S.index("hello")
            Traceback (most recent call last):
            ...
            TypeError: an integer is required
        """
        cdef mp_size_t out
        if not isinstance(other, BoundedIntegerSequence):
            try:
                out = biseq_index(self.data, other, 0)
            except OverflowError:
                out = -1
            if out >= 0:
                return out
            raise ValueError("{!r} is not in sequence".format(other))

        cdef BoundedIntegerSequence right = other
        if self.data.itembitsize != right.data.itembitsize:
            out = -1
        else:
            out = biseq_contains(self.data, right.data, 0)
        if out >= 0:
            return out
        raise ValueError("not a sub-sequence")

    def __add__(self, other):
        """
        Concatenation of bounded integer sequences.

        NOTE:

        There is no coercion happening, as bounded integer sequences are not
        considered to be elements of an object.

        EXAMPLES::

            sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(21, [4,1,6,2,7,20,9])
            sage: T = BoundedIntegerSequence(21, [4,1,6,2,8,15])
            sage: S+T
            <4, 1, 6, 2, 7, 20, 9, 4, 1, 6, 2, 8, 15>
            sage: T+S
            <4, 1, 6, 2, 8, 15, 4, 1, 6, 2, 7, 20, 9>
            sage: S in S+T
            True
            sage: T in S+T
            True
            sage: BoundedIntegerSequence(21, [4,1,6,2,7,20,9,4]) in S+T
            True
            sage: T+list(S)
            Traceback (most recent call last):
            ...
            TypeError: Cannot convert list to sage.data_structures.bounded_integer_sequences.BoundedIntegerSequence
            sage: T+None
            Traceback (most recent call last):
            ...
            TypeError: cannot concatenate bounded integer sequence and None

        TESTS:

        The discussion at :issue:`15820` explains why the following are good tests::

            sage: BoundedIntegerSequence(21, [0,0]) + BoundedIntegerSequence(21, [0,0])
            <0, 0, 0, 0>
            sage: B1 = BoundedIntegerSequence(2^30, [10^9+1, 10^9+2])
            sage: B2 = BoundedIntegerSequence(2^30, [10^9+3, 10^9+4])
            sage: B1 + B2
            <1000000001, 1000000002, 1000000003, 1000000004>
        """
        cdef BoundedIntegerSequence myself, right, out
        if other is None or self is None:
            raise TypeError('cannot concatenate bounded integer sequence and None')
        myself = self  # may result in a type error
        right = other  # --"--
        if right.data.itembitsize != myself.data.itembitsize:
            raise ValueError("can only concatenate bounded integer sequences of compatible bounds")
        out = BoundedIntegerSequence.__new__(BoundedIntegerSequence, 0, None)
        biseq_init_concat(out.data, myself.data, right.data)
        return out

    cpdef BoundedIntegerSequence maximal_overlap(self, BoundedIntegerSequence other):
        """
        Return ``self``'s maximal trailing sub-sequence that ``other`` starts with.

        Return ``None`` if there is no overlap.

        EXAMPLES::

            sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
            sage: X = BoundedIntegerSequence(21, [4,1,6,2,7,2,3])
            sage: S = BoundedIntegerSequence(21, [0,0,0,0,0,0,0])
            sage: T = BoundedIntegerSequence(21, [2,7,2,3,0,0,0,0,0,0,0,1])
            sage: (X+S).maximal_overlap(T)
            <2, 7, 2, 3, 0, 0, 0, 0, 0, 0, 0>
            sage: print((X+S).maximal_overlap(BoundedIntegerSequence(21, [2,7,2,3,0,0,0,0,0,1])))
            None
            sage: (X+S).maximal_overlap(BoundedIntegerSequence(21, [0,0]))
            <0, 0>
            sage: B1 = BoundedIntegerSequence(4,[1,2,3,2,3,2,3])
            sage: B2 = BoundedIntegerSequence(4,[2,3,2,3,2,3,1])
            sage: B1.maximal_overlap(B2)
            <2, 3, 2, 3, 2, 3>
        """
        cdef mp_size_t i = biseq_startswith_tail(other.data, self.data, 0)
        if i==-1:
            return None
        return self[i:]

    def __richcmp__(self, other, op):
        """
        Comparison of bounded integer sequences.

        We compare, in this order:

        - The bound of ``self`` and ``other``

        - The length of ``self`` and ``other``

        - Reverse lexicographical ordering, i.e., the sequences' items
          are compared starting with the last item.

        EXAMPLES:

        Comparison by bound::

            sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(21, [4,1,6,2,7,20,9])
            sage: T = BoundedIntegerSequence(51, [4,1,6,2,7,20,9])
            sage: S < T
            True
            sage: T < S
            False
            sage: list(T) == list(S)
            True

        Comparison by length::

            sage: T = BoundedIntegerSequence(21, [0,0,0,0,0,0,0,0])
            sage: S < T
            True
            sage: T < S
            False
            sage: list(T) < list(S)
            True
            sage: len(T) > len(S)
            True

        Comparison by *reverse* lexicographical ordering::

            sage: T = BoundedIntegerSequence(21, [4,1,5,2,8,20,9])
            sage: T > S
            True
            sage: S > T
            False
            sage: list(S)> list(T)
            True
        """
        cdef BoundedIntegerSequence right
        cdef BoundedIntegerSequence left
        if other is None or self is None:
            return NotImplemented
        try:
            right = other
            left = self
        except TypeError:
            return NotImplemented
        return biseq_richcmp(left.data, right.data, op)

    def __hash__(self):
        """
        The hash takes into account the content and the bound of the sequence.

        EXAMPLES::

            sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
            sage: S = BoundedIntegerSequence(21, [4,1,6,2,7,20,9])
            sage: T = BoundedIntegerSequence(51, [4,1,6,2,7,20,9])
            sage: S == T
            False
            sage: list(S) == list(T)
            True
            sage: S.bound() == T.bound()
            False
            sage: hash(S) == hash(T)
            False
            sage: T = BoundedIntegerSequence(31, [4,1,6,2,7,20,9])
            sage: T.bound() == S.bound()
            True
            sage: hash(S) == hash(T)
            True
        """
        cdef Py_hash_t h = biseq_hash(self.data)
        if h == -1:
            return 0
        return h

cpdef BoundedIntegerSequence NewBISEQ(tuple bitset_data, mp_bitcnt_t itembitsize, mp_size_t length):
    """
    Helper function for unpickling of :class:`BoundedIntegerSequence`.

    EXAMPLES::

        sage: from sage.data_structures.bounded_integer_sequences import BoundedIntegerSequence
        sage: L = [randint(0,26) for i in range(5000)]
        sage: S = BoundedIntegerSequence(32, L)
        sage: loads(dumps(S)) == S    # indirect doctest
        True

    TESTS:

    We test a corner case::

        sage: S = BoundedIntegerSequence(8,[])
        sage: S
        <>
        sage: loads(dumps(S)) == S
        True

    And another one::

        sage: S = BoundedIntegerSequence(2*sys.maxsize, [8, 8, 26, 18, 18, 8, 22, 4, 17, 22, 22, 7, 12, 4, 1, 7, 21, 7, 10, 10])
        sage: loads(dumps(S))
        <8, 8, 26, 18, 18, 8, 22, 4, 17, 22, 22, 7, 12, 4, 1, 7, 21, 7, 10, 10>
    """
    cdef BoundedIntegerSequence out = BoundedIntegerSequence.__new__(BoundedIntegerSequence)
    biseq_unpickle(out.data, bitset_data, itembitsize, length)
    return out


def _biseq_stresstest():
    """
    This function creates many bounded integer sequences and manipulates them
    in various ways, in order to try to detect random memory corruptions.

    This runs forever and must be interrupted (this means that
    interrupting is also checked).

    TESTS::

        sage: from sage.data_structures.bounded_integer_sequences import _biseq_stresstest
        sage: alarm(1); _biseq_stresstest()  # long time
        Traceback (most recent call last):
        ...
        AlarmInterrupt
    """
    cdef int branch
    cdef Py_ssize_t x, y, z
    from sage.misc.prandom import randint
    cdef list L = [BoundedIntegerSequence(6, [randint(0,5) for z in range(randint(4,10))]) for y in range(100)]
    cdef BoundedIntegerSequence S, T
    while True:
        branch = randint(0,4)
        if branch == 0:
            L[randint(0,99)] = L[randint(0,99)]+L[randint(0,99)]
        elif branch == 1:
            x = randint(0,99)
            if len(L[x]):
                y = randint(0,len(L[x])-1)
                z = randint(y,len(L[x])-1)
                L[randint(0,99)] = L[x][y:z]
            else:
                L[x] = BoundedIntegerSequence(6, [randint(0,5) for z in range(randint(4,10))])
        elif branch == 2:
            t = list(L[randint(0,99)])
            t = repr(L[randint(0,99)])
            t = L[randint(0,99)].list()
        elif branch == 3:
            x = randint(0,99)
            if len(L[x]):
                y = randint(0,len(L[x])-1)
                t = L[x][y]
                try:
                    t = L[x].index(t)
                except ValueError:
                    raise ValueError("{} should be in {} (bound {}) at position {}".format(t,L[x],L[x].bound(),y))
            else:
                L[x] = BoundedIntegerSequence(6, [randint(0,5) for z in range(randint(4,10))])
        elif branch == 4:
            S = L[randint(0,99)]
            T = L[randint(0,99)]
            biseq_startswith(S.data,T.data)
            biseq_contains(S.data, T.data, 0)
            biseq_startswith_tail(S.data, T.data, 0)
