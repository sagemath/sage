# sage_setup: distribution = sagemath-modules

from sage.modules.vector_double_dense cimport Vector_double_dense


cdef class Vector_real_double_dense(Vector_double_dense):
    pass
