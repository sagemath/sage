# sage_setup: distribution = sagemath-ntl
from sage.libs.gmp.types cimport mpq_t
from sage.rings.padics.padic_ZZ_pX_element cimport pAdicZZpXElement
from sage.structure.element cimport RingElement, ModuleElement
from sage.libs.ntl.types cimport ZZ_pX_c
from sage.libs.ntl.ntl_ZZ_pX cimport ntl_ZZ_pX
from sage.rings.padics.padic_ZZ_pX_CR_element cimport pAdicZZpXCRElement

cdef class pAdicZZpXCAElement(pAdicZZpXElement):
    cdef ZZ_pX_c value
    cdef long absprec

    cdef bint _set_prec_both_with_ordp(self, long ordp, long absprec, long relprec) except -1
    cdef int _set(self, ZZ_pX_c* value, long absprec) except -1
    cdef int _set_from_mpq_part2(self, mpq_t x) except -1

    cpdef pAdicZZpXCRElement to_fraction_field(self)
    cdef pAdicZZpXCAElement _new_c(self, long absprec)
    cdef pAdicZZpXCAElement _lshift_c(self, long n)
    cdef pAdicZZpXCAElement _rshift_c(self, long n)
    cpdef pAdicZZpXCAElement unit_part(self)
    cpdef _ntl_rep_abs(self)
    cpdef ntl_ZZ_pX _ntl_rep(self)

    cpdef pAdicZZpXCAElement lift_to_precision(self, absprec=*)
