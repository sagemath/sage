from cypari2.gen cimport Gen as pari_gen
from sage.rings.power_series_ring_element cimport PowerSeries

cdef class PowerSeries_pari(PowerSeries):
    cdef pari_gen g
