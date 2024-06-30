from sage.structure.element cimport CommutativeRingElement
from sage.structure.sage_object cimport SageObject

cdef class FiniteRingElement(CommutativeRingElement):
    pass

cdef class FinitePolyExtElement(FiniteRingElement):
    pass

cdef class Cache_base(SageObject):
    cpdef FinitePolyExtElement fetch_int(self, number)

