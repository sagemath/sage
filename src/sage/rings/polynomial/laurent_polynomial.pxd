from sage.rings.quotient_ring_element import QuotientRingElement
from sage.structure.element cimport CommutativeAlgebraElement, ModuleElement, RingElement, Element
from sage.rings.polynomial.polydict cimport ETuple, PolyDict
from sage.rings.polynomial.multi_polynomial cimport MPolynomial
from sage.rings.polynomial.polynomial_element cimport Polynomial

cdef class LaurentPolynomial_generic(CommutativeAlgebraElement):
    pass

cdef class LaurentPolynomial_univariate(LaurentPolynomial_generic):
    cpdef ModuleElement __u
    cdef long __n

cdef class LaurentPolynomial_mpair(LaurentPolynomial_generic):
    cdef ETuple _mon
    cdef MPolynomial _poly
    cdef PolyDict _prod
    cdef _new_c(self)

