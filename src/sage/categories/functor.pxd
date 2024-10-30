# sage_setup: distribution = sagemath-objects
from sage.structure.sage_object cimport SageObject

cdef class Functor(SageObject):
    cdef __weakref__
    cdef object __domain
    cdef object __codomain
