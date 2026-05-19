from sage.structure.parent_gens cimport ParentWithGens

cpdef bint _is_Field(x) except -2

cdef class Ring(ParentWithGens):
    cdef public object _zero_element
    cdef public object _one_element


cdef class CommutativeRing(Ring):
    pass

cdef class IntegralDomain(Ring):
    pass

cdef class DedekindDomain(Ring):
    pass

cdef class PrincipalIdealDomain(Ring):
    pass

cdef class Field(Ring):
    pass

cdef class Algebra(Ring):
    pass

cdef class CommutativeAlgebra(Ring):
    pass
