from sage.libs.flint.types cimport fmpz_poly_t

cdef class Apply:
    cdef fmpz_poly_t f, g, ff, gg
    cdef int apply_to_monomial_flint(self, fmpz_poly_t ans, int i, int j, int a, int b, int c, int d) except -1
