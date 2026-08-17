from sage.structure.sage_object cimport SageObject

cdef class Function(SageObject):
    cdef unsigned int _serial
    cdef int _nargs
    cdef str _name
    cdef str _alt_name
    cdef str _latex_name
    cdef object _conversions
    cdef object _evalf_params_first
    cdef _is_registered(self)
    cdef _register_function(self)

cdef class BuiltinFunction(Function):
    cdef int _preserved_arg  # 0 if none, otherwise in [1.._nargs], see function.pyx
    cdef _is_registered(self)

cdef class GinacFunction(BuiltinFunction):
    cdef str _ginac_name
    cdef _is_registered(self)
    cdef _register_function(self)

cdef class SymbolicFunction(Function):
    # cache hash value
    cdef Py_hash_t _hash_(self) except -1
    cdef bint __hinit
    cdef Py_hash_t __hcache
    cdef _is_registered(self)
