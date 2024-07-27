# sage_setup: distribution = sagemath-modules
from sage.structure.element cimport Element, ModuleElement

cdef class IndexedFreeModuleElement(ModuleElement):
    cdef public dict _monomial_coefficients
    cdef long _hash
    cdef bint _hash_set

    cpdef _add_(self, other)
    cpdef _sub_(self, other)
    cpdef _neg_(self)

    cpdef dict monomial_coefficients(self, bint copy=*)
