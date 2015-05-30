cdef class CachedFunction(object):
    cdef public str __name__
    cdef public str __module__
    cdef object _argument_fixer
    cdef public object _fix_to_pos
    cdef public object f
    # cache is not always of type "dict"!
    cdef public object cache
    cdef tuple _default_key
    cdef bint is_classmethod
    cdef argfix_init(self)
    cpdef get_cache(self)
    cpdef clear_cache(self)
    cdef object key

cdef class CachedMethod

cdef class CachedMethodCaller(CachedFunction):
    cdef public object _instance
    cdef public bint _inst_in_key
    cdef public CachedMethod _cachedmethod

cdef class CachedMethodCallerNoArgs(CachedFunction):
    cdef public object _instance

cdef class CachedMethod(object):
    cdef str _cache_name
    cdef public str __name__
    cdef public str __module__
    cdef CachedFunction _cachedfunc
    cdef int nargs
    cpdef dict _get_instance_cache(self, inst)

