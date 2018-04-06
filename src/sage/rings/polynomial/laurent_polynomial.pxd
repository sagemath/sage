from sage.structure.element cimport CommutativeAlgebraElement, ModuleElement, RingElement, Element
from sage.rings.polynomial.polydict cimport ETuple, PolyDict
from sage.rings.polynomial.multi_polynomial cimport MPolynomial


cdef class LaurentPolynomial(CommutativeAlgebraElement):
    cpdef _add_(self, other)
    cpdef _mul_(self, other)
    cpdef _floordiv_(self, other)
    cpdef long number_of_terms(self) except -1


cdef class LaurentPolynomial_univariate(LaurentPolynomial):
    cpdef ModuleElement __u
    cdef long __n


cdef class LaurentPolynomial_mpair(LaurentPolynomial):
    cdef ETuple _mon
    cdef MPolynomial _poly
    cdef PolyDict _prod
    cdef _new_c(self)
