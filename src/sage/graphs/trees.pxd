# sage_setup: distribution = sagemath-graphs

cdef class TreeIterator:
    cdef int vertices
    cdef int first_time
    cdef int p
    cdef int q
    cdef int h1
    cdef int h2
    cdef int c
    cdef int r
    cdef int *l
    cdef int *current_level_sequence
    cdef int generate_first_level_sequence(self) noexcept
    cdef int generate_next_level_sequence(self) noexcept
