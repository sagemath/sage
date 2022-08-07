"Evaluation"


def QFEvaluateVector(Q, v):
    """
    Evaluate this quadratic form Q on a vector or matrix of elements
    coercible to the base ring of the quadratic form.  If a vector
    is given then the output will be the ring element Q(v), but if a
    matrix is given then the output will be the quadratic form Q'
    which in matrix notation is given by:

    .. MATH::

            Q' = v^t * Q * v.

    Note: This is a Python wrapper for the fast evaluation routine
    QFEvaluateVector_cdef().  This routine is for internal use and is
    called more conveniently as Q(M).

    INPUT:

    - Q -- QuadraticForm over a base ring R
    - v -- a tuple or list (or column matrix) of Q.dim() elements of R

    OUTPUT:

        an element of R

    EXAMPLES::

        sage: from sage.quadratic_forms.quadratic_form__evaluate import QFEvaluateVector
        sage: Q = QuadraticForm(ZZ, 4, range(10)); Q
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 0 1 2 3 ]
        [ * 4 5 6 ]
        [ * * 7 8 ]
        [ * * * 9 ]
        sage: QFEvaluateVector(Q, (1,0,0,0))
        0
        sage: QFEvaluateVector(Q, (1,0,1,0))
        9

    """
    return QFEvaluateVector_cdef(Q, v)



cdef QFEvaluateVector_cdef(Q, v):
    """
    Routine to quickly evaluate a quadratic form Q on a vector v.  See
    the Python wrapper function QFEvaluate() above for details.

    """
    # If we are passed a matrix A, return the quadratic form Q(A(x))
    # (In matrix notation: A^t * Q * A)
    n = Q.dim()

    tmp_val = Q.base_ring()(0)
    for i from 0 <= i < n:
        for j from i <= j < n:
            tmp_val += Q[i,j] * v[i] * v[j]

    # Return the value (over R)
    return Q.base_ring().coerce(tmp_val)



def QFEvaluateMatrix(Q, M, Q2):
    """
    Evaluate this quadratic form Q on a matrix M of elements coercible
    to the base ring of the quadratic form, which in matrix notation
    is given by:

            Q2 = M^t * Q * M.

    Note: This is a Python wrapper for the fast evaluation routine
    QFEvaluateMatrix_cdef().  This routine is for internal use and is
    called more conveniently as Q(M).  The inclusion of Q2 as an
    argument is to avoid having to create a QuadraticForm here, which
    for now creates circular imports.

    INPUT:

    - Q -- QuadraticForm over a base ring R
    - M -- a Q.dim() x Q2.dim() matrix of elements of R

    OUTPUT:

    - Q2 -- a QuadraticForm over R

    EXAMPLES::

        sage: from sage.quadratic_forms.quadratic_form__evaluate import QFEvaluateMatrix
        sage: Q = QuadraticForm(ZZ, 4, range(10)); Q
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 0 1 2 3 ]
        [ * 4 5 6 ]
        [ * * 7 8 ]
        [ * * * 9 ]
        sage: Q2 = QuadraticForm(ZZ, 2)
        sage: M = Matrix(ZZ, 4, 2, [1,0,0,0, 0,1,0,0]); M
        [1 0]
        [0 0]
        [0 1]
        [0 0]
        sage: QFEvaluateMatrix(Q, M, Q2)
        Quadratic form in 2 variables over Integer Ring with coefficients:
        [ 0 2 ]
        [ * 7 ]

    """
    return QFEvaluateMatrix_cdef(Q, M, Q2)


cdef QFEvaluateMatrix_cdef(Q, M, Q2):
    """
    Routine to quickly evaluate a quadratic form Q on a matrix M.  See
    the Python wrapper function QFEvaluateMatrix() above for details.

    """
    # Create the new quadratic form
    n = Q.dim()
    m = Q2.dim()

    # TODO: Check the dimensions of M are compatible with those of Q and Q2

    # Evaluate Q(M) into Q2
    for k from 0 <= k < m:
        for l from k <= l < m:
            tmp_sum = Q2.base_ring()(0)
            for i from 0 <= i < n:
                for j from i <= j < n:
                    if (k == l):
                        tmp_sum += Q[i,j] * (M[i,k] * M[j,l])
                    else:
                        tmp_sum += Q[i,j] * (M[i,k] * M[j,l] + M[i,l] * M[j,k])
            Q2[k,l] = tmp_sum
    return Q2
