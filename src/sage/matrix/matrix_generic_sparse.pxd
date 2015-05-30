cimport matrix_sparse

cdef class Matrix_generic_sparse(matrix_sparse.Matrix_sparse):
    cdef dict _entries
    cdef object _zero

