# sage_setup: distribution = sagemath-ntl
ctypedef double celement

include "matrix_modn_dense_template_header.pxi"

from sage.rings.finite_rings.integer_mod cimport IntegerMod_abstract

cdef class Matrix_modn_dense_double(Matrix_modn_dense_template):
    cdef IntegerMod_abstract _get_template
    cdef bint _fits_int32
