from .generic_backend cimport GenericBackend
from sage.matrix.matrix2 cimport Matrix

cdef class MatrixBackend(GenericBackend):

    cdef Matrix objective_function
    cdef Matrix G_matrix
    cdef str prob_name
    cdef bint is_maximize

    cdef Matrix row_lower_bound
    cdef list row_lower_bound_indicator
    cdef Matrix row_upper_bound
    cdef list row_upper_bound_indicator
    cdef Matrix col_lower_bound
    cdef list col_lower_bound_indicator
    cdef Matrix col_upper_bound
    cdef list col_upper_bound_indicator

    cdef list row_name_var
    cdef list col_name_var

    cdef list is_integer
    cdef list is_binary
    cdef list is_continuous

    cdef object _base_ring
