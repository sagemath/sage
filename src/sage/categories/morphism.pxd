# sage_setup: distribution = sagemath-objects
from sage.structure.element cimport Element
from sage.categories.map cimport Map


cdef class Morphism(Map):
    pass

cdef class SetMorphism(Morphism):
    cdef object _function
    cpdef bint _eq_c_impl(left, Element right) noexcept

cdef class SetIsomorphism(SetMorphism):
    cdef object _inverse
