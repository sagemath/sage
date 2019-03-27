from sage.structure.element cimport AlgebraElement, ModuleElement

cdef class LaurentSeries(AlgebraElement):
    cpdef ModuleElement __u
    cdef long __n

    cdef __normalize(self)
    cpdef _add_(self, other)
    cpdef _mul_(self, other)

