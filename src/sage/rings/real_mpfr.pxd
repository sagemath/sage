from sage.libs.mpfr cimport *

cimport sage.rings.ring
cimport sage.structure.element

cdef extern from "pari/pari.h":
    ctypedef long* GEN


cdef class RealNumber(sage.structure.element.RingElement)  # forward decl

cdef class RealField_class(sage.rings.ring.Field):
    cdef int __prec
    cdef bint sci_not
    cdef mpfr_rnd_t rnd
    cdef object rnd_str
    cdef RealNumber _new(self)


cdef class RealNumber(sage.structure.element.RingElement):
    cdef mpfr_t value
    cdef char init
    cdef RealNumber _new(self)
    cdef _set(self, x, int base)
    cdef _set_from_GEN_REAL(self, GEN g)
    cdef RealNumber abs(RealNumber self)

cpdef RealField(int prec=*, int sci_not=*, rnd=*)
