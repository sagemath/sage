# sage_setup: distribution = sagemath-categories
cdef class PolyDict:
    cdef dict __repn

    cdef PolyDict _new(self, dict pdict)
    cpdef remove_zeros(self, zero_test=*)


cdef class ETuple:
    cdef size_t _length
    cdef size_t _nonzero
    cdef int *_data

    cdef size_t get_position(self, size_t i, size_t start, size_t end) noexcept
    cdef ETuple _new(self)
    cdef int get_exp(self, size_t i) noexcept

    # need a cdef version for function pointers
    cdef int _unweighted_degree(self) except *
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
    cpdef int dotprod(self, ETuple other) except *
    cpdef ETuple escalar_div(self, int n)
    cpdef ETuple divide_by_gcd(self, ETuple other)
    cpdef ETuple divide_by_var(self, size_t pos)
    cpdef bint divides(self, ETuple other) except *
    cpdef bint is_constant(self) noexcept
    cpdef bint is_multiple_of(self, int n) except *
    cpdef list nonzero_positions(self, bint sort=*)
    cpdef common_nonzero_positions(self, ETuple other, bint sort=*)
    cpdef list nonzero_values(self, bint sort=*)
    cpdef ETuple reversed(self)

cpdef int gen_index(PolyDict x) noexcept
cpdef ETuple monomial_exponent(PolyDict p)
