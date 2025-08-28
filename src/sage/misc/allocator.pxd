from cpython.object cimport *

cdef hook_tp_functions_type(object t, newfunc tp_new, destructor tp_dealloc, bint useGC)

cdef hook_tp_functions(object global_dummy, newfunc tp_new, destructor tp_dealloc, bint useGC)
