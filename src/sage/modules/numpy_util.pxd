from libc.stdint cimport uintptr_t
from sage.libs.m4ri cimport *

cpdef int set_mzd_from_numpy(uintptr_t entries_addr, Py_ssize_t degree, x) except -1
# Note: we don't actually need ``cimport`` to work, which means this header file is not used in practice
# neither do we need ``cpdef`` (``def`` is sufficient)
