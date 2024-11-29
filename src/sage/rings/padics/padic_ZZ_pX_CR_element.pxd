from sage.libs.gmp.types cimport mpz_t, mpq_t
from sage.rings.padics.padic_ZZ_pX_element cimport pAdicZZpXElement
from sage.structure.element cimport RingElement, ModuleElement
from sage.libs.ntl.types cimport ZZ_pX_c, ZZX_c
from sage.libs.ntl.ntl_ZZ_pX cimport ntl_ZZ_pX

cdef class pAdicZZpXCRElement(pAdicZZpXElement):
    cdef ZZ_pX_c unit
    cdef long ordp
    cdef long relprec
    # relprec = 0 if self.unit is not constructed.  This includes _exact_zeros and _inexact_zeros
    # relprec > 0 if self.unit is constructed and self is normalized
    # relprec < 0 if self.unit is constructed but self is not normalized.  The actual relprec is -self.relprec

    cdef int _set(self, ZZ_pX_c* unit, long ordp, long relprec) except -1
    cdef int _set_from_mpq_part1(self, mpz_t num_unit, mpz_t den_unit, mpq_t x) except -1
    cdef int _set_from_mpq_part2(self, mpz_t num_unit, mpz_t den_unit) except -1
    cdef int _set_from_ZZX_part1(self, ZZX_c poly, long absprec, long relprec) except -1
    cdef int _set_from_ZZ_pX_part1(self, ZZ_pX_c* poly) except -1
    cdef int _set_from_ZZ_pX_part2(self, ZZ_pX_c* poly) except -1

    cdef pAdicZZpXCRElement _new_c(self, long relprec)
    cdef int _internal_lshift(self, long shift) except -1
    cdef int _normalize(self) except -1
    cdef pAdicZZpXCRElement _lshift_c(self, long n)
    cdef pAdicZZpXCRElement _rshift_c(self, long n)
    cpdef pAdicZZpXCRElement unit_part(self)
    cpdef ntl_ZZ_pX _ntl_rep_unnormalized(self)
    cpdef _ntl_rep_abs(self)
    cpdef ntl_ZZ_pX _ntl_rep(self)

    cpdef pAdicZZpXCRElement lift_to_precision(self, absprec=*)
