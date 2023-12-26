from sage.matrix.matrix_sparse cimport Matrix_sparse


cdef class Matrix_generic_sparse(Matrix_sparse):
    cdef dict _entries
    cdef object _zero
