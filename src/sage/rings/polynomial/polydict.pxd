cdef class PolyDict:
    cdef dict __repn
    cdef object __zero

cdef class ETuple:
    cdef size_t _length
    cdef size_t _nonzero
    cdef int *_data

    cpdef size_t unweighted_degree(self)
    cdef size_t weighted_degree(self, tuple w)
    cdef size_t unweighted_quotient_degree(self, ETuple other)
    cdef size_t weighted_quotient_degree(self, ETuple other, tuple w)
    cpdef ETuple eadd(ETuple self, ETuple other)
    cpdef ETuple esub(ETuple self, ETuple other)
    cpdef ETuple emul(ETuple self, int factor)
    cpdef ETuple emin(ETuple self, ETuple other)
    cpdef ETuple emax(ETuple self, ETuple other)
    cpdef ETuple eadd_p(ETuple self, int other, int pos)
    cpdef ETuple eadd_scaled(ETuple self, ETuple other, int scalar)
    cpdef int dotprod(ETuple self, ETuple other)
    cpdef ETuple escalar_div(ETuple self, int n)
    cdef ETuple divide_by_gcd(self, ETuple other)
    cdef ETuple divide_by_var(self, size_t index)
    cdef bint divides(self, ETuple other)
    cpdef bint is_constant(ETuple self)
    cpdef bint is_multiple_of(ETuple self, int n)
    cpdef list nonzero_positions(ETuple self, bint sort=*)
    cpdef common_nonzero_positions(ETuple self, ETuple other, bint sort=*)
    cpdef list nonzero_values(ETuple self, bint sort=*)
    cpdef ETuple reversed(ETuple self)
    cdef ETuple _new(ETuple self)
    cdef size_t get_exp(ETuple self, int i)

