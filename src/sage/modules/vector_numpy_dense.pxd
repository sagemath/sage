# sage_setup: distribution = sagemath-modules

cimport numpy

from sage.modules.free_module_element cimport FreeModuleElement


cdef class Vector_numpy_dense(FreeModuleElement):
    cdef object _numpy_dtype
    cdef int _numpy_dtypeint
    cdef object _python_dtype
    cdef object _sage_dtype
    cdef object _sage_vector_dtype
    cdef numpy.ndarray _vector_numpy
    cdef Vector_numpy_dense _new(self, numpy.ndarray vector_numpy)
    cdef _replace_self_with_numpy(self, numpy.ndarray numpy_array)
