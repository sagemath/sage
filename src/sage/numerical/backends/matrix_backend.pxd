from .generic_backend cimport GenericBackend
#from sage.matrix.matrix2 cimport Matrix

cdef class MatrixBackend(GenericBackend):

    cdef object objective_function
    cdef object G_matrix
    cdef str prob_name
    cdef bint is_maximize

    cdef object row_lower_bound
    cdef object row_upper_bound
    cdef object col_lower_bound
    cdef object col_upper_bound

    cdef list row_name_var
    cdef list col_name_var

    cdef list is_integer
    cdef list is_binary
    cdef list is_continuous

    cdef object _base_ring
