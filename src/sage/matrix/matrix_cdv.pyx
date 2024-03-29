from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense
from sage.structure.element cimport RingElement
from sage.rings.infinity import Infinity

# Function definition
cpdef hessenbergize_cdvf(Matrix_generic_dense H) except *:
    """
    Compute the Hessenberg form of the given matrix over a complete discrete
    valuation field.

    This function assumes that H is a square matrix over a complete discrete
    valuation field. The pivot on each column is always chosen with maximal
    relative precision, ensuring the numerical stability of the algorithm.

    INPUT:
    - H: Matrix over a complete discrete valuation field

    OUTPUT:
    None (H is modified in-place)

    EXAMPLES:
    sage: K = Qp(5, print_mode="digits", prec=5)
    sage: H = matrix(K, 3, 3, range(9))
    sage: H.hessenbergize()
    sage: H
    [        0  ...00010  ...00002]
    [ ...00003  ...00024 ...000010]
    [ ...00000  ...44440  ...44443]
    """
    