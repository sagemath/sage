# sage_setup: distribution = sagemath-objects
cdef class _LazyString():
    cdef func
    cdef args
    cdef kwargs
    cdef val(self)
    cpdef update_lazy_string(self, args, kwds)
