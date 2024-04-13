# sage_setup: distribution = sagemath-ntl
ctypedef float celement

include "matrix_modn_dense_template_header.pxi"

from sage.rings.finite_rings.integer_mod cimport IntegerMod_int

cdef class Matrix_modn_dense_float(Matrix_modn_dense_template):
    cdef IntegerMod_int _get_template
