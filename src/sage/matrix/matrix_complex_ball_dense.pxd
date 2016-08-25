from sage.libs.arb.types cimport acb_mat_t
from .matrix_dense cimport Matrix_dense
from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense
from sage.structure.parent cimport Parent

cdef void matrix_to_acb_mat(acb_mat_t target, source)
cdef Matrix_generic_dense acb_mat_to_matrix(
    acb_mat_t source, Parent CIF)

cdef class Matrix_complex_ball_dense(Matrix_dense):
    cdef acb_mat_t value
