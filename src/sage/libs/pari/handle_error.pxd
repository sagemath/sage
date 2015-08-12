from .types cimport GEN

cdef void _pari_init_error_handling()
cdef int _pari_err_handle(GEN E) except 0
cdef void _pari_err_recover(long errnum)
