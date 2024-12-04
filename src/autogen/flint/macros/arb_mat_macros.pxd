# Macros from arb_mat.h
# See https://github.com/flintlib/flint/issues/1529

from .types cimport *

cdef extern from "flint_wrap.h":
    arb_ptr arb_mat_entry(arb_mat_t mat, slong i, slong j)
    slong arb_mat_nrows(arb_mat_t)
    slong arb_mat_ncols(arb_mat_t)
