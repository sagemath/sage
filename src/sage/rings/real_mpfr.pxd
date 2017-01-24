from sage.libs.mpfr cimport *

cimport sage.rings.ring
cimport sage.structure.element
from sage.libs.cypari2.types cimport GEN


cdef class RealNumber(sage.structure.element.RingElement)  # forward decl

cdef class RealField_class(sage.rings.ring.Field):
    cdef int __prec
    cdef bint sci_not
    cdef mpfr_rnd_t rnd
    cdef object rnd_str
    cdef inline RealNumber _new(self):
        """Return a new real number with parent ``self``."""
        return <RealNumber>(RealNumber.__new__(RealNumber, self))

cdef class RealNumber(sage.structure.element.RingElement):
    cdef mpfr_t value
    cdef inline RealNumber _new(self):
        """Return a new real number with same parent as ``self``."""
        return <RealNumber>(RealNumber.__new__(RealNumber, self._parent))
    cdef _set(self, x, int base)
    cdef _set_from_GEN_REAL(self, GEN g)
    cdef RealNumber abs(RealNumber self)

cpdef RealField(int prec=*, int sci_not=*, rnd=*)
