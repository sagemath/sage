from sage.libs.gsl.types cimport gsl_complex

cimport sage.structure.element
cimport sage.rings.ring
cimport sage.structure.element
from sage.structure.element cimport RingElement, ModuleElement
from cypari2.types cimport GEN


cdef class ComplexDoubleField_class(sage.rings.ring.Field):
    pass

cdef class ComplexDoubleElement(sage.structure.element.FieldElement):
    cdef gsl_complex _complex
    cdef ComplexDoubleElement _new_c(self, gsl_complex x)
    cpdef _add_(self, other)
    cpdef _mul_(self, other)

cdef ComplexDoubleElement new_ComplexDoubleElement()
