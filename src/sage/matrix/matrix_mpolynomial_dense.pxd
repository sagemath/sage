# sage_setup: distribution = sagemath-singular
from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense

from sage.libs.singular.decl cimport ideal

cdef class Matrix_mpolynomial_dense(Matrix_generic_dense):
    pass
