# sage_setup: distribution = sagemath-categories
from sage.structure.element cimport AlgebraElement, ModuleElement
from sage.rings.laurent_series_ring_element cimport LaurentSeries

cdef class PuiseuxSeries(AlgebraElement):
     cdef LaurentSeries _l
     cdef size_t _e
