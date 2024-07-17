from sage.libs.mpfr.types cimport mpfr_prec_t
from sage.libs.mpfi.types cimport mpfi_t

cimport sage.structure.element
from sage.rings.real_mpfi cimport RealIntervalFieldElement, RealIntervalField_class


cdef class ComplexIntervalFieldElement(sage.structure.element.FieldElement):
    cdef mpfi_t __re
    cdef mpfi_t __im
    cdef mpfr_prec_t _prec
    cdef object _multiplicative_order

    cdef inline ComplexIntervalFieldElement _new(self):
        """
        Quickly create a new complex interval with the same parent as
        ``self``.
        """
        cdef type t = type(self)
        cdef object _multiplicative_order = None
        return t.__new__(t, self._parent)

    cdef inline RealIntervalFieldElement _new_real(self):
        """
        Quickly create a new real interval with the same precision as
        ``self``.
        """
        P = <RealIntervalField_class>(self._parent.real_field())
        return P._new()
