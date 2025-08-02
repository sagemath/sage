cdef class LazyImport:
    cdef readonly object _object   # The actual object if imported, None otherwise
    cdef object _module
    cdef object _name
    cdef object _as_name
    cdef object _namespace
    cdef bint _at_startup
    cdef object _deprecation
    cdef object _feature

    cdef inline get_object(self)
    cpdef _get_object(self)
