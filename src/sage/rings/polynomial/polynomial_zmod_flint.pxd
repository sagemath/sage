from sage.libs.flint.types cimport nmod_poly_t, nmod_poly_struct, fmpz_poly_t
from sage.structure.parent cimport Parent
from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint

ctypedef nmod_poly_struct celement
ctypedef unsigned long cparent

include "polynomial_template_header.pxi"

cdef cparent get_cparent(parent) except? 0

cdef class Polynomial_zmod_flint(Polynomial_template):
    cdef Polynomial_template _new(self)
    cdef int _set_list(self, x) except -1
    cdef int _set_fmpz_poly(self, fmpz_poly_t) except -1
    cpdef _mul_trunc(self, Polynomial_zmod_flint other, length)
    cpdef _mul_trunc_opposite(self, Polynomial_zmod_flint other, length)
    cpdef rational_reconstruct(self, m, n_deg=?, d_deg=?)
