# sage_setup: distribution = sagemath-objects
from sage.categories.map cimport Map


cdef class DefaultConvertMap(Map):
    pass


cdef class DefaultConvertMap_unique(DefaultConvertMap):
    pass


cdef class NamedConvertMap(Map):
    cdef readonly method_name


cdef class TryMap(Map):
    cdef Map _map_p
    cdef Map _map_b
    cdef _error_types


cdef class CallableConvertMap(Map):
    cdef bint _parent_as_first_arg
    cdef _func


cdef Map CCallableConvertMap(domain, codomain, void* func, name)
