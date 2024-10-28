# sage_setup: distribution = sagemath-objects
cdef class WeakValueDictionary(dict):
    cdef __weakref__
    cdef callback
    cdef int _guard_level
    cdef list _pending_removals

    cdef int _set_item(self, key, value) except -1
    cdef int _enter_iter(self) except -1
    cdef int _exit_iter(self) except -1


cdef class CachedWeakValueDictionary(WeakValueDictionary):
    cdef tuple cache
    cdef Py_ssize_t cache_index
