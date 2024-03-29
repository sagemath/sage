cdef class _LazyString():
    cdef func
    cdef args
    cdef kwargs
    cdef val(self) noexcept
    cpdef update_lazy_string(self, args, kwds) noexcept
