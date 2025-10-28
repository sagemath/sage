# Macros from acb_mat.h
# See https://github.com/flintlib/flint/issues/1529

from .types cimport *

cdef extern from "flint_wrap.h":
    acb_ptr acb_mat_entry(acb_mat_t mat, slong i, slong j)
    slong acb_mat_nrows(acb_mat_t)
    slong acb_mat_ncols(acb_mat_t)
