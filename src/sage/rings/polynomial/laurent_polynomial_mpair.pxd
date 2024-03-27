from sage.rings.polynomial.laurent_polynomial cimport LaurentPolynomial
from sage.rings.polynomial.multi_polynomial cimport MPolynomial
from sage.rings.polynomial.polydict cimport ETuple, PolyDict


cdef class LaurentPolynomial_mpair(LaurentPolynomial):
    cdef ETuple _mon
    cdef MPolynomial _poly
    cdef PolyDict _prod
    cdef _compute_polydict(self)
    cdef _normalize(self, i=*)
    cpdef rescale_vars(self, dict d, h=*, new_ring=*)
    cpdef toric_coordinate_change(self, M, h=*, new_ring=*)
    cpdef toric_substitute(self, v, v1, a, h=*, new_ring=*)
