from sage.libs.singular.decl cimport skStrategy, ring

from sage.rings.polynomial.multi_polynomial_libsingular cimport \
    MPolynomialRing_libsingular, MPolynomial_libsingular
from sage.structure.sage_object cimport SageObject
from sage.rings.polynomial.plural cimport NCPolynomialRing_plural, NCPolynomial_plural

cdef class GroebnerStrategy(SageObject):
    cdef skStrategy *_strat
    cdef ring *_parent_ring
    cdef MPolynomialRing_libsingular _parent
    cdef object _ideal

    cpdef MPolynomial_libsingular normal_form(self, MPolynomial_libsingular p)

cdef class NCGroebnerStrategy(SageObject):
    cdef skStrategy *_strat
    cdef NCPolynomialRing_plural _parent
    cdef object _ideal

    cpdef NCPolynomial_plural normal_form(self, NCPolynomial_plural p)
