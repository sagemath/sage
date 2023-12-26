# sage_setup: distribution = sagemath-categories
cdef class Envelope():
    cdef readonly sign
    cdef f
    cdef f_limit_start
    cdef list precomputed
    cdef readonly max_part
    cdef readonly min_slope, max_slope

cdef class IntegerListsBackend():
    cdef readonly min_sum, max_sum
    cdef readonly min_length, max_length
    cdef readonly min_part, max_part
    cdef readonly min_slope, max_slope
    cdef readonly Envelope floor, ceiling
    cdef public dict _cached_methods  # Support cached_method
