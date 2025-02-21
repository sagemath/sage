from sage.libs.mpfr.types cimport mpfr_rnd_t, mpfr_t, mpfr_prec_t

cimport sage.rings.abc
cimport sage.structure.element

cdef class RealNumber(sage.structure.element.RingElement)  # forward decl

cdef class RealField_class(sage.rings.abc.RealField):
    cdef mpfr_prec_t _prec
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
    cpdef _add_(self, other)
    cpdef _mul_(self, other)
    cpdef _mod_(self, right)
    cdef _set(self, x, int base)
    cdef RealNumber abs(RealNumber self)

cpdef RealField(mpfr_prec_t prec=*, int sci_not=*, rnd=*)
