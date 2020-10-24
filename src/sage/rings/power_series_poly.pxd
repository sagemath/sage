from .power_series_ring_element cimport PowerSeries
from sage.rings.polynomial.polynomial_element cimport Polynomial


cdef class PowerSeries_poly(PowerSeries):
    cdef Polynomial __f
