# sage_setup: distribution = sagemath-modules
from sage.libs.mpc.types cimport mpc_t, mpc_rnd_t

cimport sage.structure.element
from sage.rings.ring cimport Field

cdef class MPComplexNumber(sage.structure.element.FieldElement):
    cdef mpc_t value
    cdef char init
    cdef MPComplexNumber _new(self)
    cpdef _add_(self, other)
    cpdef _mul_(self, other)

cdef class MPComplexField_class(Field):
    cdef readonly int _prec
    cdef mpc_rnd_t __rnd
    cdef object __rnd_str
    cdef object __real_field
    cdef object __imag_field
    cdef MPComplexNumber _new(self)
    cpdef _an_element_(self)
