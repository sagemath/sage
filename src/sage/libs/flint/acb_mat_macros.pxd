# Macros from acb_mat.h
# See https://github.com/flintlib/flint/issues/152

from .types cimport *

cdef extern from "flint_wrap.h":
    acb_ptr acb_mat_entry(acb_mat_t mat, slong i, slong j)
    # Macro giving a pointer to the entry at row *i* and column *j*.

    slong acb_mat_nrows(acb_mat_t)
    # Returns the number of rows of the matrix.

    slong acb_mat_ncols(acb_mat_t)
    # Returns the number of columns of the matrix.
