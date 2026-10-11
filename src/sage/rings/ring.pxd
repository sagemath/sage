from sage.structure.parent cimport Parent

cpdef bint _is_Field(x) except -2

cdef class Ring(Parent):
    cdef public object _zero_element
    cdef public object _one_element
    cdef public object _gens
    cdef public object _latex_names

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
