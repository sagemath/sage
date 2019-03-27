cpdef int iaxpy(a, dict X, dict Y, bint remove_zeros=*, bint factor_on_left=*) except -1
cpdef dict axpy(a, dict X, dict Y, bint factor_on_left=*)
cpdef dict negate(dict D)
cpdef dict scal(a, dict D, bint factor_on_left=*)
cpdef dict add(dict D, dict D2)
cpdef dict sum(dict_iter)
cpdef dict linear_combination(dict_factor_iter, bint factor_on_left=*)

