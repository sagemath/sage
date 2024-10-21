# sage_setup: distribution = sagemath-modules

from sage.modules.vector_numpy_dense cimport Vector_numpy_dense


cdef class Vector_numpy_integer_dense(Vector_numpy_dense):

    pass
