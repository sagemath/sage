cdef class PolyDict:
    cdef dict __repn

    cdef PolyDict _new(self, dict pdict) noexcept
    cpdef remove_zeros(self, zero_test=*) noexcept


cdef class ETuple:
    cdef size_t _length
    cdef size_t _nonzero
    cdef int *_data

    cdef size_t get_position(self, size_t i, size_t start, size_t end) noexcept
    cdef ETuple _new(self) noexcept
    cdef int get_exp(self, size_t i) noexcept

    cpdef int unweighted_degree(self) except *
    cpdef int weighted_degree(self, tuple w) except *
    cpdef int unweighted_quotient_degree(self, ETuple other) except *
    cpdef int weighted_quotient_degree(self, ETuple other, tuple w) except *

    cpdef ETuple eadd(self, ETuple other) noexcept
    cpdef ETuple esub(self, ETuple other) noexcept
    cpdef ETuple emul(self, int factor) noexcept
    cpdef ETuple emin(self, ETuple other) noexcept
    cpdef ETuple emax(self, ETuple other) noexcept
    cpdef ETuple eadd_p(self, int other, size_t pos) noexcept
    cpdef ETuple eadd_scaled(self, ETuple other, int scalar) noexcept
    cpdef int dotprod(self, ETuple other) except *
    cpdef ETuple escalar_div(self, int n) noexcept
    cpdef ETuple divide_by_gcd(self, ETuple other) noexcept
    cpdef ETuple divide_by_var(self, size_t pos) noexcept
    cpdef bint divides(self, ETuple other) except *
    cpdef bint is_constant(self) noexcept
    cpdef bint is_multiple_of(self, int n) except *
    cpdef list nonzero_positions(self, bint sort=*) noexcept
    cpdef common_nonzero_positions(self, ETuple other, bint sort=*) noexcept
    cpdef list nonzero_values(self, bint sort=*) noexcept
    cpdef ETuple reversed(self) noexcept

cpdef int gen_index(PolyDict x) noexcept
cpdef ETuple monomial_exponent(PolyDict p) noexcept
