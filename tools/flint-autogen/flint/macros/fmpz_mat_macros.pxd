# Macros from fmpz_mat.h
# See https://github.com/flintlib/flint/issues/1529

from .types cimport *

cdef extern from "flint_wrap.h":
    fmpz * fmpz_mat_entry(fmpz_mat_t mat, slong i, slong j)
    slong fmpz_mat_nrows(fmpz_mat_t)
    slong fmpz_mat_ncols(fmpz_mat_t)
