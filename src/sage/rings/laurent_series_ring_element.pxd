from sage.structure.element cimport AlgebraElement, ModuleElement

cdef class LaurentSeries(AlgebraElement):
    cdef ModuleElement __u
    cdef long __n

    cdef _normalize(self)
    cpdef _add_(self, other)
    cpdef _mul_(self, other)

