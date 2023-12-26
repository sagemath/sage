from sage.structure.element cimport Element, ModuleElement

cdef class IndexedFreeModuleElement(ModuleElement):
    cdef public dict _monomial_coefficients
    cdef long _hash
    cdef bint _hash_set

    cpdef _add_(self, other) noexcept
    cpdef _sub_(self, other) noexcept
    cpdef _neg_(self) noexcept

    cpdef dict monomial_coefficients(self, bint copy=*) noexcept
