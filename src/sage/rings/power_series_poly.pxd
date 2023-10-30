from sage.rings.power_series_ring_element cimport PowerSeries
from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.categories.action cimport Action

cdef class PowerSeries_poly(PowerSeries):
    cdef Polynomial __f

cdef class BaseRingFloorDivAction(Action):
    pass

