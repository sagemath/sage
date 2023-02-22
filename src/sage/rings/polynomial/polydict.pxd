cdef class PolyDict:
    cdef dict __repn

    cdef PolyDict _new(self, dict pdict)
    cpdef remove_zeros(self)


cdef class ETuple:
    cdef size_t _length
    cdef size_t _nonzero
    cdef int *_data

    cpdef int unweighted_degree(self) except *
    cpdef int weighted_degree(self, tuple w) except *
    cpdef int unweighted_quotient_degree(self, ETuple other) except *
    cpdef int weighted_quotient_degree(self, ETuple other, tuple w) except *

    cpdef ETuple eadd(self, ETuple other)
    cpdef ETuple esub(self, ETuple other)
    cpdef ETuple emul(self, int factor)
    cpdef ETuple emin(self, ETuple other)
    cpdef ETuple emax(self, ETuple other)
    cpdef ETuple eadd_p(self, int other, size_t pos)
    cpdef ETuple eadd_scaled(self, ETuple other, int scalar)
    cpdef int dotprod(self, ETuple other)
    cpdef ETuple escalar_div(self, int n)
    cdef ETuple divide_by_gcd(self, ETuple other)
    cdef ETuple divide_by_var(self, size_t index)
    cdef bint divides(self, ETuple other)
    cpdef bint is_constant(self)
    cpdef bint is_multiple_of(self, int n)
    cpdef list nonzero_positions(self, bint sort=*)
    cpdef common_nonzero_positions(self, ETuple other, bint sort=*)
    cpdef list nonzero_values(self, bint sort=*)
    cpdef ETuple reversed(self)
    cdef ETuple _new(self)
    cdef int get_exp(self, size_t i)

cpdef int gen_index(PolyDict x)
cpdef ETuple monomial_exponent(PolyDict p)
