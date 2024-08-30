from sage.libs.singular.decl cimport ring, poly, number, intvec, bigintmat
from sage.libs.singular.function cimport Resolution

from sage.rings.rational cimport Rational
from sage.structure.element cimport Element
from sage.rings.integer cimport Integer
from sage.rings.finite_rings.integer_mod cimport IntegerMod_abstract
from sage.rings.finite_rings.element_givaro cimport Cache_givaro
from sage.rings.finite_rings.element_givaro cimport FiniteField_givaroElement as FFgivE
from sage.rings.finite_rings.element_ntl_gf2e cimport Cache_ntl_gf2e
from sage.rings.finite_rings.element_ntl_gf2e cimport FiniteField_ntl_gf2eElement as FFgf2eE
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomialRing_libsingular

from sage.rings.number_field.number_field_base cimport NumberField

# ======================================
# Conversion from Singular to Sage types
# ======================================

cdef Rational si2sa_QQ(number (*), number **, ring (*))
cdef Integer  si2sa_ZZ(number (*),ring (*))

cdef FFgivE   si2sa_GFqGivaro(number *n, ring *_ring, Cache_givaro cache)
cdef FFgf2eE  si2sa_GFqNTLGF2E(number *n, ring *_ring, Cache_ntl_gf2e cache)
cdef object   si2sa_GFq_generic(number *n, ring *_ring, object base)
cdef object   si2sa_ZZmod(number *n, ring *_ring, object base)

cdef object   si2sa_NF(number *n, ring *_ring, object base)

cdef object si2sa_intvec(intvec *v)
cdef object si2sa_bigintvec(bigintmat *v)

# dispatches to all the above.
cdef object si2sa(number *n, ring *_ring, object base)

cdef list singular_monomial_exponents(poly *p, ring *r)
cpdef list si2sa_resolution(Resolution res)
cpdef tuple si2sa_resolution_graded(Resolution res, tuple degrees)

# ======================================
# Conversion from Sage to Singular types
# ======================================

cdef number *sa2si_QQ(Rational ,ring (*)) noexcept
cdef number *sa2si_ZZ(Integer d, ring *_ring) noexcept

cdef number *sa2si_GFqGivaro(int exp ,ring (*)) noexcept
cdef number *sa2si_GFqNTLGF2E(FFgf2eE elem, ring *_ring) noexcept
cdef number *sa2si_GFq_generic(object vector, ring *_ring) noexcept
cdef number *sa2si_ZZmod(IntegerMod_abstract d, ring *_ring) noexcept

cdef number *sa2si_NF(object element, ring *_ring) noexcept

# dispatches to all the above.
cdef number *sa2si(Element elem, ring * _ring) noexcept

# ==============
# Initialisation
# ==============

cdef int overflow_check(unsigned long e, ring *_ring) except -1

cdef init_libsingular()
