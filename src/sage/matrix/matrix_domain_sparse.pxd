# sage_setup: distribution = sagemath-modules
from .matrix cimport Matrix

cdef class Matrix_domain_sparse(Matrix):
    pass
