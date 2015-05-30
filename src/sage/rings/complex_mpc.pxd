from sage.libs.mpc cimport *

cimport sage.rings.ring

cdef class MPComplexNumber(sage.structure.element.FieldElement):
    cdef mpc_t value
    cdef char init
    cdef MPComplexNumber _new(self)

cdef class MPComplexField_class(sage.rings.ring.Field):
    cdef readonly int __prec
    cdef mpc_rnd_t __rnd
    cdef object __rnd_str
    cdef object __real_field
    cdef object __imag_field
    cdef MPComplexNumber _new(self)
    cpdef _an_element_(self)
