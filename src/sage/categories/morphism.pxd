from sage.structure.element cimport Element
from .map cimport Map


cdef class Morphism(Map):
    pass

cdef class SetMorphism(Morphism):
    cdef object _function
    cpdef bint _eq_c_impl(left, Element right)
