from sage.libs.eclib cimport mat, scalar, scalar_to_mpz
from sage.rings.integer cimport Integer

cdef class Matrix:
    cdef mat* M

cdef class MatrixFactory:
    cdef new_matrix(self, mat M)

cdef inline Integer scalar_to_integer(const scalar& s):
    cdef Integer res = Integer()
    scalar_to_mpz(res.value, s)
    return res
