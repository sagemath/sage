from sage.categories.action cimport Action
from sage.categories.map cimport Map

cdef class DefaultConvertMap(Map):
    cdef public bint _force_use

cdef class DefaultConvertMap_unique(DefaultConvertMap):
    pass

cdef class NamedConvertMap(Map):
    cdef readonly method_name
    cdef public bint _force_use

cdef class TryMap(Map):
    cdef Map _map_p
    cdef Map _map_b
    cdef _error_types

cdef Map CCallableConvertMap(domain, codomain, void* func, name)
