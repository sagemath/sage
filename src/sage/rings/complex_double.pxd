from sage.libs.gsl.types cimport gsl_complex

cimport sage.structure.element
cimport sage.rings.ring
cimport sage.structure.element
from sage.structure.element cimport RingElement, ModuleElement
from sage.libs.pari.types cimport GEN


cdef class ComplexDoubleField_class(sage.rings.ring.Field):
    pass

cdef class ComplexDoubleElement(sage.structure.element.FieldElement):
    cdef gsl_complex _complex
    cdef GEN _gen(self)
    cdef ComplexDoubleElement _new_c(self, gsl_complex x)

cdef ComplexDoubleElement new_ComplexDoubleElement()
