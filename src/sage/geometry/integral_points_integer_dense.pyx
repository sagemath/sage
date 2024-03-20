# sage_setup: distribution = sagemath-linbox
#cython: wraparound=False, boundscheck=False

from sage.modules.vector_integer_dense cimport Vector_integer_dense as VectorClass
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense as MatrixClass

include "integral_points.pxi"
