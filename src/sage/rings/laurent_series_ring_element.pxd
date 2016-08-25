from sage.structure.element cimport AlgebraElement, ModuleElement

cdef class LaurentSeries(AlgebraElement):
    cpdef ModuleElement __u
    cdef long __n
