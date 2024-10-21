# sage_setup: distribution = sagemath-categories
cdef class ContainChildren():
    cdef int parentpid
    cdef int exitcode, exceptcode
    cdef bint silent
