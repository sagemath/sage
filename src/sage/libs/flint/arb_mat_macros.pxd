# Macros from arb_mat.h
# See https://github.com/flintlib/flint/issues/152

from .types cimport *

cdef extern from "flint_wrap.h":
    arb_ptr arb_mat_entry(arb_mat_t mat, slong i, slong j)
    # Macro giving a pointer to the entry at row *i* and column *j*.

    slong arb_mat_nrows(arb_mat_t)
    # Returns the number of rows of the matrix.

    slong arb_mat_ncols(arb_mat_t)
    # Returns the number of columns of the matrix.
