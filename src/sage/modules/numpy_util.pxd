from libc.stdint cimport uintptr_t
from sage.libs.m4ri cimport *

cpdef int set_mzd_from_numpy(uintptr_t entries_addr, Py_ssize_t degree, x) except -1
