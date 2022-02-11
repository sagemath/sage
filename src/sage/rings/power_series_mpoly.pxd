from sage.structure.element cimport ModuleElement
from .power_series_ring_element cimport PowerSeries

cdef class PowerSeries_mpoly(PowerSeries):
    cdef ModuleElement __f
    cdef object __poly
    cdef object __list
    cdef bint _truncated
