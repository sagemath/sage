# sage_setup: distribution = sagemath-categories
cdef class SymmetricReductionStrategy:
    cdef list _lm
    cdef list _lengths
    cdef object _min_lm
    cdef int _tail
    cdef object _R
    cdef object _parent
