from libc.stdint cimport uintptr_t
from sage.libs.m4ri cimport *
from sage.matrix.matrix_mod2_dense cimport Matrix_mod2_dense

cpdef int set_matrix_mod2_from_numpy(Matrix_mod2_dense a, b) except -1

cpdef int set_mzd_from_numpy(uintptr_t entries_addr, Py_ssize_t degree, x) except -1
# Note: we don't actually need ``cimport`` to work, which means this header file is not used in practice
# neither do we need ``cpdef`` (``def`` is sufficient)
