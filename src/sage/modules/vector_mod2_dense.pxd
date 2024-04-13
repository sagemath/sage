# sage_setup: distribution = sagemath-linbox
from sage.modules.free_module_element cimport FreeModuleElement
from sage.libs.m4ri cimport mzd_t


cdef class Vector_mod2_dense(FreeModuleElement):
    cdef mzd_t* _entries
    cdef object _base_ring

    cdef _new_c(self)
    cdef _init(self, Py_ssize_t degree, parent)
