# sage_setup: distribution = sagemath-pari
from cpython cimport array

cimport sage.structure.element
from sage.libs.gmp.types cimport mpz_t, mpq_t
from sage.structure.element cimport Element, RingElement
from sage.rings.padics.local_generic_element cimport LocalGenericElement
from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational

cpdef gauss_table(long long p, int f, int prec, bint use_longs)

cdef class pAdicGenericElement(LocalGenericElement):
    cdef long valuation_c(self) noexcept
    cpdef val_unit(self)

    cdef int _set_from_Integer(self, Integer x, absprec, relprec) except -1
    cdef int _set_from_mpz(self, mpz_t x) except -1
    cdef int _set_from_mpz_rel(self, mpz_t x, long relprec) except -1
    cdef int _set_from_mpz_abs(self, mpz_t value, long absprec) except -1
    cdef int _set_from_mpz_both(self, mpz_t x, long absprec, long relprec) except -1

    cdef int _set_from_Rational(self, Rational x, absprec, relprec) except -1
    cdef int _set_from_mpq(self, mpq_t x) except -1
    cdef int _set_from_mpq_rel(self, mpq_t x, long relprec) except -1
    cdef int _set_from_mpq_abs(self, mpq_t value, long absprec) except -1
    cdef int _set_from_mpq_both(self, mpq_t x, long absprec, long relprec) except -1

    cdef int _pshift_self(self, long shift) except -1

    cdef int _cmp_units(left, pAdicGenericElement right) except -2

    cdef int _set_inexact_zero(self, long absprec) except -1
    cdef int _set_exact_zero(self) except -1

    cpdef bint _is_exact_zero(self) except -1
    cpdef bint _is_inexact_zero(self) except -1
    cpdef bint _is_zero_rep(self) except -1

    cdef bint _set_prec_abs(self, long absprec) except -1
    cdef bint _set_prec_rel(self, long relprec) except -1
    cdef bint _set_prec_both(self, long absprec, long relprec) except -1

    cpdef abs(self, prec=*)
    cpdef _mod_(self, right)
    cpdef _floordiv_(self, right)
    cpdef bint _is_base_elt(self, p) except -1
