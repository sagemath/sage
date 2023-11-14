from sage.rings.polynomial.laurent_polynomial cimport LaurentPolynomial
from sage.rings.polynomial.multi_polynomial cimport MPolynomial
from sage.rings.polynomial.polydict cimport ETuple, PolyDict


cdef class LaurentPolynomial_mpair(LaurentPolynomial):
    cdef ETuple _mon
    cdef MPolynomial _poly
    cdef PolyDict _prod
    cdef _compute_polydict(self) noexcept
    cdef _normalize(self, i=*) noexcept
    cpdef rescale_vars(self, dict d, h=*, new_ring=*) noexcept
    cpdef toric_coordinate_change(self, M, h=*, new_ring=*) noexcept
    cpdef toric_substitute(self, v, v1, a, h=*, new_ring=*) noexcept
