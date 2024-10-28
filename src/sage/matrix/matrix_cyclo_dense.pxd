from sage.libs.gmp.types cimport mpz_t
from sage.matrix.matrix_dense cimport Matrix_dense
from sage.matrix.matrix_rational_dense cimport Matrix_rational_dense

cdef class Matrix_cyclo_dense(Matrix_dense):

    # Matrix over ZZ that stores elements
    cdef Matrix_rational_dense _matrix

    # Degree of base cyclotomic field
    cdef int _degree
    cdef int _n

    cdef _randomize_rational_column_unsafe(Matrix_cyclo_dense self,
        Py_ssize_t col, mpz_t nump1, mpz_t denp1, distribution=?)

