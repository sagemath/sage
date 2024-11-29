from sage.rings.ring cimport CommutativeRing, Ring
from sage.structure.parent cimport Parent

cdef class MPolynomialRing_base(CommutativeRing):
    cdef object _ngens
    cdef object _term_order
    cdef public object _indices
    cdef public object _has_singular
    cdef public object _magma_gens
    cdef public dict _magma_cache

    cdef _coerce_c_impl(self, x)


cdef class BooleanPolynomialRing_base(MPolynomialRing_base):
    pass
