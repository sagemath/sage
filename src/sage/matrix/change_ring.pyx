# Importing necessary modules
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.matrix_real_double_dense cimport Matrix_real_double_dense
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.rings.real_double import RDF

def integer_to_real_double_dense(Matrix_integer_dense A):
    """
    Fast conversion of a matrix over the integers to a matrix with
    real double entries.

    INPUT:
        A -- a dense matrix over the integers

    OUTPUT:
        -- a dense real double matrix

    EXAMPLES::

        sage: a = matrix(ZZ,2,3,[-2,-5,3,4,8,1030339830489349908])
        sage: a.change_ring(RDF)
        [                  -2.0                   -5.0                    3.0]
        [                   4.0                    8.0 1.0303398304893499e+18]
        sage: import sage.matrix.change_ring
        sage: sage.matrix.change_ring.integer_to_real_double_dense(a)
        [                  -2.0                   -5.0                    3.0]
        [                   4.0                    8.0 1.0303398304893499e+18]
    """
    cdef Py_ssize_t i, j
    cdef Matrix_real_double_dense M
    
    # Creating a MatrixSpace with real double entries
    S = MatrixSpace(RDF, A._nrows, A._ncols, sparse=False)
    M = S()
    
    # Copying elements from the integer matrix to the real double matrix
    for i in range(A._nrows):
        for j in range(A._ncols):
            M[i, j] = RDF(A[i, j])
    
    return M
