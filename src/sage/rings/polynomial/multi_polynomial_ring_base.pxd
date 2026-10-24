from sage.rings.ring cimport Ring
from sage.structure.parent cimport Parent

cdef class MPolynomialRing_base(Ring):
    cdef object _ngens
    cdef object _term_order
    cdef public object _indices
    cdef public object _has_singular
    cdef public object _magma_gens
    cdef public dict _magma_cache


cdef class BooleanPolynomialRing_base(MPolynomialRing_base):
    pass
