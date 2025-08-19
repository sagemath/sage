from sage.structure.element cimport AlgebraElement
from sage.rings.power_series_ring_element cimport PowerSeries

cdef class LaurentSeries(AlgebraElement):
    cdef PowerSeries __u
    cdef long __n

    cdef _normalize(self)
    cpdef _add_(self, other)
    cpdef _mul_(self, other)
