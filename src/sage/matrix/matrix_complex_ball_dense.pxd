from sage.libs.arb.types cimport acb_mat_t
from sage.matrix.matrix_dense cimport Matrix_dense
from sage.structure.parent cimport Parent

# Helper function declarations
cdef void matrix_to_acb_mat(acb_mat_t target, source) except +
cdef Matrix_generic_dense acb_mat_to_matrix(
    acb_mat_t source, Parent CIF) except +

# Class definition
cdef class Matrix_complex_ball_dense(Matrix_dense):
    cdef acb_mat_t value
    
    # Constructor declaration
    cdef Matrix_complex_ball_dense _new(self, Py_ssize_t nrows, Py_ssize_t ncols) except +
    
    # Exposed method declaration
    cpdef _pow_int(self, n) except +
