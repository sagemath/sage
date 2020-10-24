from ..eclib cimport mat

cdef class Matrix:
    cdef mat* M

cdef class MatrixFactory:
    cdef new_matrix(self, mat M)




