from sage.structure.element cimport AlgebraElement, ModuleElement

cdef class LaurentSeries(AlgebraElement):
    cdef ModuleElement __u
    cdef long __n

    cdef _normalize(self) noexcept
    cpdef _add_(self, other) noexcept
    cpdef _mul_(self, other) noexcept

