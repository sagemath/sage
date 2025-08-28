from sage.libs.gmp.types cimport *
from sage.data_structures.bitset cimport *

# A biseq (bounded integer sequence) is a sequence of nonnegative
# integers, each fitting in "itembitsize" bits. We store the sequence
# in a bitset of length at least length*itembitsize.
ctypedef struct biseq_s:
    # A bitset comprising at least length*itembitsize bits is used to
    # store the items of this sequence. Note that bitset_t stores the
    # data in a field mp_limb_t *bits, where mp_limb_t is defined in
    # sage.libs.gmp.types.
    #
    # Normally, the "at least" in the first sentence above will be
    # an equality, but we need the "at least" for the special case
    # where the length is zero: then we will use a bitset of length 1.
    # All extra bits (not in the first length*itembitsize) should be
    # zero.
    bitset_t data

    # Number of items in this sequence.
    mp_size_t length

    # Bitsize (ranging from 1 to GMP_LIMB_BITS) of one item of this
    # sequence. Note: Each item is a nonnegative integer, and all
    # items of this sequence satisfy an upper bound. We do not store
    # the exact bound for the items of this sequence, but store the
    # bitsize that is sufficient to store one item.
    mp_bitcnt_t itembitsize

    # If L is a limb, then "L & mask_item" extracts the first
    # itembitsize bits of it, i.e., it extracts one item.
    size_t mask_item

ctypedef biseq_s biseq_t[1]

cdef bint biseq_init(biseq_t R, mp_size_t l, mp_bitcnt_t itemsize) except -1
cdef void biseq_dealloc(biseq_t S) noexcept
cdef bint biseq_init_copy(biseq_t R, biseq_t S) except -1
cdef tuple biseq_pickle(biseq_t S)
cdef bint biseq_unpickle(biseq_t R, tuple bitset_data, mp_bitcnt_t itembitsize, mp_size_t length) except -1
cdef bint biseq_init_list(biseq_t R, list data, size_t bound) except -1
cdef Py_hash_t biseq_hash(biseq_t S) noexcept
cdef bint biseq_richcmp(biseq_t S1, biseq_t S2, int op) noexcept
cdef bint biseq_init_concat(biseq_t R, biseq_t S1, biseq_t S2) except -1
cdef bint biseq_startswith(biseq_t S1, biseq_t S2) except -1
cdef mp_size_t biseq_contains(biseq_t S1, biseq_t S2, mp_size_t start) except -2
cdef mp_size_t biseq_startswith_tail(biseq_t S1, biseq_t S2, mp_size_t start) except -2
cdef mp_size_t biseq_index(biseq_t S, size_t item, mp_size_t start) except -2
cdef size_t biseq_getitem(biseq_t S, mp_size_t index) noexcept
cdef biseq_getitem_py(biseq_t S, mp_size_t index)
cdef void biseq_inititem(biseq_t S, mp_size_t index, size_t item) noexcept
cdef void biseq_clearitem(biseq_t S, mp_size_t index) noexcept
cdef bint biseq_init_slice(biseq_t R, biseq_t S, mp_size_t start, mp_size_t stop, mp_size_t step) except -1

cdef class BoundedIntegerSequence:
    cdef biseq_t data
    cpdef bint startswith(self, BoundedIntegerSequence other) noexcept
    cpdef list list(self)
    cpdef BoundedIntegerSequence maximal_overlap(self, BoundedIntegerSequence other)

cpdef BoundedIntegerSequence NewBISEQ(tuple bitset_data, mp_bitcnt_t itembitsize, mp_size_t length)
