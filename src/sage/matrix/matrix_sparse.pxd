from sage.matrix.matrix cimport Matrix

cdef class Matrix_sparse(Matrix):
    cdef dict _multiply_classical_entries(self, Matrix_sparse right)
