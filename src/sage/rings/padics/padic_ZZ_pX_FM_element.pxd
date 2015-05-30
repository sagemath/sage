include "sage/libs/ntl/decl.pxi"

from sage.rings.padics.padic_ZZ_pX_element cimport pAdicZZpXElement
from sage.structure.element cimport RingElement, ModuleElement

cdef class pAdicZZpXFMElement(pAdicZZpXElement):
    cdef ZZ_pX_c value
    cdef pAdicZZpXFMElement _new_c(self)
    cdef pAdicZZpXFMElement _lshift_c(self, long n)
    cdef pAdicZZpXFMElement _rshift_c(self, long n)

    cpdef pAdicZZpXFMElement unit_part(self)
