from sage.structure.element cimport ModuleElement
from sage.rings.power_series_ring_element cimport PowerSeries

cdef class PowerSeries_mpoly(PowerSeries):
    cdef ModuleElement __f
    cdef object _poly
    cdef object __list
    cdef bint _truncated
