from sage.libs.mpfi cimport *

cimport sage.structure.element
from .real_mpfi cimport RealIntervalFieldElement

cdef class ComplexIntervalFieldElement(sage.structure.element.FieldElement):
    cdef mpfi_t __re
    cdef mpfi_t __im
    cdef mp_prec_t _prec

    cdef RealIntervalFieldElement abs_c(ComplexIntervalFieldElement self)
    cdef RealIntervalFieldElement norm_c(ComplexIntervalFieldElement self)

    cdef inline ComplexIntervalFieldElement _new(self):
        """
        Quickly create a new complex interval with the same parent as
        ``self``.
        """
        cdef type t = type(self)
        return t.__new__(t, self._parent)
