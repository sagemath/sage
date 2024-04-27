from sage.structure.sage_object cimport SageObject

cdef class Function(SageObject):
    cdef unsigned int _serial
    cdef int _nargs
    cdef object _name
    cdef object _alt_name
    cdef object _latex_name
    cdef object _conversions
    cdef object _evalf_params_first
    cdef _is_registered(self)
    cdef _register_function(self)

cdef class BuiltinFunction(Function):
    cdef object _preserved_arg
    cdef _is_registered(self)

cdef class GinacFunction(BuiltinFunction):
    cdef object _ginac_name
    cdef _is_registered(self)
    cdef _register_function(self)

cdef class SymbolicFunction(Function):
    # cache hash value
    cdef long _hash_(self) except -1
    cdef bint __hinit
    cdef long __hcache
    cdef _is_registered(self)
