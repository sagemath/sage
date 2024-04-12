# sage_setup: distribution = sagemath-polyhedra
from sage.numerical.backends.generic_sdp_backend cimport GenericSDPBackend


cdef class MatrixSDPBackend(GenericSDPBackend):

    cdef list objective_function
    cdef list coeffs_matrix
    cdef bint is_maximize

    cdef list row_name_var
    cdef list col_name_var
    cdef str name

    cdef object _base_ring
