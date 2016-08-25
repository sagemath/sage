r"""
A binary matrix datatype in Cython

It's almost a copy of the bitset datatype, but allows a differentiation of the
rows of the matrix. That's the only advantage compared to storing all the rows
of a matrix in a loooooonng bitset.

A ``binary_matrix_t`` structure contains :

- ``Py_ssize_t n_cols`` -- number of columns

- ``Py_ssize_t n_rows`` -- number of rows

- ``bitset_t * rows`` -- ``rows[i]`` points toward a block of type ``bitset_t``
  containing the bits of row `i`.

"""

from sage.data_structures.binary_matrix cimport *
include 'sage/data_structures/bitset.pxi'


cdef inline binary_matrix_init(binary_matrix_t m, Py_ssize_t n_rows, Py_ssize_t n_cols):
    r"""
    Allocate an empty binary matrix with ``n_rows`` rows and ``n_cols`` columns.
    """
    cdef Py_ssize_t i

    m.n_cols = n_cols
    m.n_rows = n_rows
    m.rows = <bitset_t *>sig_malloc(n_rows * sizeof(bitset_t))
    if m.rows == NULL:
        raise MemoryError

    for i in range(n_rows):
        bitset_init(m.rows[i], n_cols)

cdef inline binary_matrix_free(binary_matrix_t m):
    r"""
    Frees the memory allocated by the matrix
    """
    cdef Py_ssize_t i

    for i in range(m.n_rows):
        bitset_free(m.rows[i])
    sig_free(m.rows)

cdef inline binary_matrix_fill(binary_matrix_t m, bint bit):
    r"""
    Fill the whole matrix with a bit
    """
    cdef Py_ssize_t i

    if bit: # set the matrix to 1
        for i in range(m.n_rows):
            bitset_set_first_n(m.rows[i], m.n_cols)
    else:
        for i in range(m.n_rows):
            bitset_clear(m.rows[i])

cdef inline binary_matrix_complement(binary_matrix_t m):
    r"""
    Complements all of the matrix' bits.
    """
    cdef Py_ssize_t i
    for i in range(m.n_rows):
        bitset_complement(m.rows[i], m.rows[i])

cdef inline binary_matrix_set1(binary_matrix_t m, Py_ssize_t row, Py_ssize_t col):
    r"""
    Sets an entry to 1
    """
    bitset_add(m.rows[row], col)

cdef inline binary_matrix_set0(binary_matrix_t m, Py_ssize_t row, Py_ssize_t col):
    r"""
    Sets an entry to 0
    """
    bitset_discard(m.rows[row], col)

cdef inline binary_matrix_set(binary_matrix_t m, Py_ssize_t row, Py_ssize_t col, bint value):
    r"""
    Sets an entry
    """
    bitset_set_to(m.rows[row],col,value)

cdef inline bint binary_matrix_get(binary_matrix_t m, Py_ssize_t row, Py_ssize_t col):
    r"""
    Returns the value of a given entry
    """
    return bitset_in(m.rows[row], col)

cdef inline binary_matrix_print(binary_matrix_t m):
    r"""
    Prints the binary matrix
    """
    cdef Py_ssize_t i,j
    import sys
    for i in range(m.n_rows):
        for j in range(m.n_cols):
            sys.stdout.write("1" if binary_matrix_get(m, i, j) else ".",)
        print("")
