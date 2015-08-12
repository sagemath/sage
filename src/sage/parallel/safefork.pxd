cdef class ContainChildren(object):
    cdef int parentpid
    cdef int exitcode, exceptcode
    cdef bint silent
