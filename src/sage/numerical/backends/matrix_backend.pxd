from .generic_backend cimport GenericBackend

cdef class MatrixBackend(GenericBackend):

    cdef list objective_function
    cdef list G_matrix
    cdef str prob_name
    cdef bint is_maximize

    cdef list row_lower_bound
    cdef list row_upper_bound
    cdef list col_lower_bound
    cdef list col_upper_bound

    cdef list row_name_var
    cdef list col_name_var

    cdef list is_integer
    cdef list is_binary
    cdef list is_continuous

    cdef object _base_ring
