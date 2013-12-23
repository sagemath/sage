include "sage/ext/cdefs.pxi"

cimport sage.rings.padics.padic_base_generic_element
from sage.rings.padics.padic_base_generic_element cimport pAdicBaseGenericElement

cimport sage.structure.element
from sage.structure.element cimport CommutativeRingElement, RingElement, ModuleElement, Element

cimport sage.rings.integer
from sage.rings.integer cimport Integer

cimport sage.rings.rational
from sage.rings.rational cimport Rational

cimport sage.rings.padics.pow_computer
from sage.rings.padics.pow_computer cimport PowComputer_class

from sage.libs.pari.gen cimport gen as pari_gen
from sage.libs.pari.pari_instance cimport PariInstance

cdef class pAdicCappedRelativeElement(pAdicBaseGenericElement):
    cdef mpz_t unit    #An exact zero is indicated by unit < 0
    cdef long ordp
    cdef long relprec
    cdef bint _normalized

    cdef int _set_zero(pAdicCappedRelativeElement self, absprec) except -1
    cdef int _set_prec(pAdicCappedRelativeElement self, long relprec) except -1
    cdef int _set(pAdicCappedRelativeElement self, long ordp, mpz_t unit, long relprec) except -1
    cdef int _set_from_CR_rel(pAdicCappedRelativeElement self, pAdicCappedRelativeElement other, long relprec) except -1
    cdef int _set_from_CR_both(pAdicCappedRelativeElement self, pAdicCappedRelativeElement other, long absprec, long relprec) except -1

    cdef pAdicCappedRelativeElement _new_c(pAdicCappedRelativeElement self)
    cdef int _normalize(pAdicCappedRelativeElement self) except -1

    cdef int _set_from_Integer( pAdicCappedRelativeElement self, Integer x, absprec, relprec) except -1
    cdef int _set_from_Rational( pAdicCappedRelativeElement self, Rational x, absprec, relprec) except -1

    cpdef ModuleElement _neg_(self)
    cpdef RingElement _floordiv_c_impl(self, RingElement right)
    cdef pAdicCappedRelativeElement _lshift_c(pAdicCappedRelativeElement self, long shift)
    cdef pAdicCappedRelativeElement _rshift_c(pAdicCappedRelativeElement self, long shift)
    cdef pAdicCappedRelativeElement lift_to_precision_c(pAdicCappedRelativeElement self, long absprec)
    cdef pari_gen _to_gen(pAdicCappedRelativeElement self)
    cdef lift_c(self)
    cdef teichmuller_list(pAdicCappedRelativeElement self)
    cpdef pAdicCappedRelativeElement unit_part(self)

