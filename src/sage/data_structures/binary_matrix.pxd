r"""
A binary matrix datatype in Cython

It's almost a copy of the bitset datatype, but allows a differentiation of the
rows of the matrix. That's the only advantage compared to storing all the rows
of a matrix in a loooooonng bitset.

A ``binary_matrix_t`` structure contains:

- ``mp_bitcnt_t n_cols`` -- number of columns

- ``mp_bitcnt_t n_rows`` -- number of rows

- ``bitset_t * rows`` -- ``rows[i]`` points toward a block of type ``bitset_t``
  containing the bits of row `i`
"""

from sage.data_structures.bitset_base cimport *

cdef struct binary_matrix_s:
    mp_bitcnt_t n_cols
    mp_bitcnt_t n_rows
    bitset_t* rows

ctypedef binary_matrix_s[1] binary_matrix_t

cdef inline int binary_matrix_init(binary_matrix_t m, mp_bitcnt_t n_rows, mp_bitcnt_t n_cols) except -1:
    r"""
    Allocate an empty binary matrix with ``n_rows`` rows and ``n_cols`` columns.
    """
    cdef mp_bitcnt_t i

    m.n_cols = n_cols
    m.n_rows = n_rows
    m.rows = <bitset_t *>sig_malloc(n_rows * sizeof(bitset_t))
    if m.rows == NULL:
        raise MemoryError

    for i in range(n_rows):
        bitset_init(m.rows[i], n_cols)

cdef inline int binary_matrix_realloc(binary_matrix_t m, mp_bitcnt_t n_rows, mp_bitcnt_t n_cols) except -1:
    r"""
    Reallocate a binary matrix to have ``n_rows`` rows and ``n_cols`` columns.

    Possibly new positions do not contain extra bits.
    """
    cdef mp_bitcnt_t i
    for i in range(n_rows, m.n_rows):
        bitset_free(m.rows[i])

    m.rows = <bitset_t *>check_reallocarray(m.rows, n_rows, sizeof(bitset_t))

    for i in range(min(n_rows, m.n_rows)):
        bitset_realloc(m.rows[i], n_cols)

    for i in range(m.n_rows, n_rows):
        bitset_init(m.rows[i], n_cols)

    m.n_cols = n_cols
    m.n_rows = n_rows

cdef inline void binary_matrix_free(binary_matrix_t m) noexcept:
    r"""
    Free the memory allocated by the matrix
    """
    cdef mp_bitcnt_t i

    for i in range(m.n_rows):
        bitset_free(m.rows[i])
    sig_free(m.rows)

cdef inline void binary_matrix_copy(binary_matrix_t dst, binary_matrix_t src) noexcept:
    """
    Copy the binary matrix src over to the binary matrix dst, overwriting dst.

    We assume that ``dst.n_rows == src.n_rows`` and ``dst.n_cols == src.n_cols``.
    """
    cdef mp_bitcnt_t i
    for i in range(dst.n_rows):
        bitset_copy(dst.rows[i], src.rows[i])

cdef inline void binary_matrix_fill(binary_matrix_t m, bint bit) noexcept:
    r"""
    Fill the whole matrix with a bit
    """
    cdef mp_bitcnt_t i

    if bit: # set the matrix to 1
        for i in range(m.n_rows):
            bitset_set_first_n(m.rows[i], m.n_cols)
    else:
        for i in range(m.n_rows):
            bitset_clear(m.rows[i])

cdef inline void binary_matrix_complement(binary_matrix_t m) noexcept:
    r"""
    Complement all of the matrix' bits.
    """
    cdef mp_bitcnt_t i
    for i in range(m.n_rows):
        bitset_complement(m.rows[i], m.rows[i])

cdef inline void binary_matrix_set1(binary_matrix_t m, mp_bitcnt_t row, mp_bitcnt_t col) noexcept:
    r"""
    Set an entry to 1
    """
    bitset_add(m.rows[row], col)

cdef inline void binary_matrix_set0(binary_matrix_t m, mp_bitcnt_t row, mp_bitcnt_t col) noexcept:
    r"""
    Set an entry to 0
    """
    bitset_discard(m.rows[row], col)

cdef inline void binary_matrix_set(binary_matrix_t m, mp_bitcnt_t row, mp_bitcnt_t col, bint value) noexcept:
    r"""
    Set an entry
    """
    bitset_set_to(m.rows[row],col,value)

cdef inline bint binary_matrix_get(binary_matrix_t m, mp_bitcnt_t row, mp_bitcnt_t col) noexcept:
    r"""
    Return the value of a given entry
    """
    return bitset_in(m.rows[row], col)

cdef inline binary_matrix_print(binary_matrix_t m):
    r"""
    Print the binary matrix
    """
    cdef mp_bitcnt_t i,j
    import sys
    for i in range(m.n_rows):
        for j in range(m.n_cols):
            sys.stdout.write("1" if binary_matrix_get(m, i, j) else ".",)
        print("")
