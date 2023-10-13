# Macros from fmpq_mat.h
# See https://github.com/flintlib/flint/issues/152

from .types cimport *

cdef extern from "flint_wrap.h":
    fmpq * fmpq_mat_entry(fmpq_mat_t mat, slong i, slong j)
    # Macro giving a pointer to the entry at row *i* and column *j*.

    slong fmpq_mat_nrows(fmpq_mat_t)
    # Returns the number of rows of the matrix.

    slong fmpq_mat_ncols(fmpq_mat_t)
    # Returns the number of columns of the matrix.
