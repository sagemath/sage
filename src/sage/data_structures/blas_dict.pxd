cpdef int iaxpy(a, dict X, dict Y, bint remove_zeros=*, bint factor_on_left=*) except -1
cpdef dict axpy(a, dict X, dict Y, bint factor_on_left=*) noexcept
cpdef dict negate(dict D) noexcept
cpdef dict scal(a, dict D, bint factor_on_left=*) noexcept
cpdef dict add(dict D, dict D2) noexcept
cpdef dict sum(dict_iter) noexcept
cpdef dict linear_combination(dict_factor_iter, bint factor_on_left=*) noexcept
cpdef dict sum_of_monomials(monomials, scalar) noexcept
cpdef dict sum_of_terms(index_coeff_pairs) noexcept
cdef dict remove_zeros(dict D) noexcept
cpdef dict convert_remove_zeroes(dict D, R) noexcept
