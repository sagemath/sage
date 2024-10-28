from sage.libs.gsl.types cimport gsl_complex

cimport sage.structure.element
cimport sage.rings.abc


cdef class ComplexDoubleField_class(sage.rings.abc.ComplexDoubleField):
    pass


cdef class ComplexDoubleElement(sage.structure.element.FieldElement):
    cdef gsl_complex _complex
    cdef ComplexDoubleElement _new_c(self, gsl_complex x)
    cpdef _add_(self, other)
    cpdef _mul_(self, other)
    cpdef _pow_(self, other)


cdef ComplexDoubleElement new_ComplexDoubleElement()
