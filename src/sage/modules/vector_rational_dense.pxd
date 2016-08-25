from .free_module_element cimport FreeModuleElement
from sage.libs.gmp.types cimport mpq_t

cdef class Vector_rational_dense(FreeModuleElement):
        cdef mpq_t* _entries
        cdef _new_c(self)
        cdef _init(self, Py_ssize_t degree, parent)

