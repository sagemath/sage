# Macros from fmpq_mat.h
# See https://github.com/flintlib/flint/issues/1529

from .types cimport *

cdef extern from "flint_wrap.h":
    fmpq * fmpq_mat_entry(fmpq_mat_t mat, slong i, slong j)
    slong fmpq_mat_nrows(fmpq_mat_t)
    slong fmpq_mat_ncols(fmpq_mat_t)
