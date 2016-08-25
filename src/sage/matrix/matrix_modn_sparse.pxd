# distutils: language = c 

from .matrix_sparse cimport Matrix_sparse

include 'sage/modules/vector_modn_sparse_h.pxi'

cdef class Matrix_modn_sparse(Matrix_sparse):
    cdef c_vector_modint* rows
    cdef public int p
    cdef swap_rows_c(self, Py_ssize_t n1, Py_ssize_t n2)

    cdef _init_linbox(self)
