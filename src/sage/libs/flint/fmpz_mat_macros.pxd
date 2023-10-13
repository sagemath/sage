# Macros from fmpz_mat.h
# See https://github.com/flintlib/flint/issues/152

from .types cimport *

cdef extern from "flint_wrap.h":
    fmpz * fmpz_mat_entry(fmpz_mat_t mat, slong i, slong j)
    # Macro giving a pointer to the entry at row *i* and column *j*.

    slong fmpz_mat_nrows(fmpz_mat_t)
    # Returns the number of rows of the matrix.

    slong fmpz_mat_ncols(fmpz_mat_t)
    # Returns the number of columns of the matrix.
