include "sage/ext/cdefs.pxi"
include "fplll.pxi"

cdef class FP_LLL:
    cdef object fp_map
    cdef ZZ_mat[mpz_t] *_lattice
