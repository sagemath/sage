cdef class _LazyString(object):
    cdef func
    cdef args
    cdef kwargs
    cdef value(self)
    cpdef update_lazy_string(self, args, kwds)
