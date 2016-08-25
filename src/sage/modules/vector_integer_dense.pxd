from .free_module_element cimport FreeModuleElement
from sage.libs.gmp.types cimport mpz_t

cdef class Vector_integer_dense(FreeModuleElement):
        cdef mpz_t* _entries
        cdef _new_c(self)
        cdef _init(self, Py_ssize_t degree, parent)

