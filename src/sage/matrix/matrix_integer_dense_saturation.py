"""
Saturation over ZZ
"""

from copy import copy

from sage.arith.misc import binomial, GCD as gcd
from sage.matrix import matrix_integer_dense_hnf
from sage.matrix.constructor import identity_matrix, random_matrix
from sage.misc.randstate import current_randstate
from sage.misc.verbose import verbose
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.rings.integer_ring import ZZ


def p_saturation(A, p, proof=True):
    """
    INPUT:

    - ``A`` -- a matrix over ZZ
    - ``p`` -- a prime
    - ``proof`` -- boolean (default: ``True``)

    OUTPUT:

    The p-saturation of the matrix A, i.e., a new matrix in Hermite form
    whose row span a ZZ-module that is p-saturated.

    EXAMPLES::

        sage: from sage.matrix.matrix_integer_dense_saturation import p_saturation
        sage: A = matrix(ZZ, 2, 2, [3,2,3,4]); B = matrix(ZZ, 2,3,[1,2,3,4,5,6])
        sage: A.det()
        6
        sage: C = A*B; C
        [11 16 21]
        [19 26 33]
        sage: C2 = p_saturation(C, 2); C2
        [ 1  8 15]
        [ 0  9 18]
        sage: C2.index_in_saturation()
        9
        sage: C3 = p_saturation(C, 3); C3
        [ 1  0 -1]
        [ 0  2  4]
        sage: C3.index_in_saturation()
        2
    """
    tm = verbose("%s-saturating a %sx%s matrix" % (p, A.nrows(), A.ncols()))
    H = A.hermite_form(include_zero_rows=False, proof=proof)
    while True:
        if p == 2:
            A = H.change_ring(GF(p))
        else:
            try:
                # Faster than change_ring
                A = H._reduce(p)
            except OverflowError:
                # fall back to generic GF(p) matrices
                A = H.change_ring(GF(p))
        assert A.nrows() <= A.ncols()
        K = A.kernel()
        if K.dimension() == 0:
            verbose("done saturating", tm)
            return H
        B = K.basis_matrix().lift()
        C = ((B * H) / p).change_ring(ZZ)
        H = H.stack(C).hermite_form(include_zero_rows=False, proof=proof)
    verbose("done saturating", tm)

def random_sublist_of_size(k, n):
    """
    INPUT:

    - ``k`` -- integer
    - ``n`` -- integer

    OUTPUT: a randomly chosen sublist of ``range(k)`` of size `n`

    EXAMPLES::

        sage: import sage.matrix.matrix_integer_dense_saturation as s
        sage: l = s.random_sublist_of_size(10, 3)
        sage: len(l)
        3
        sage: l_check = [-1] + l + [10]
        sage: all(l_check[i] < l_check[i+1] for i in range(4))
        True
        sage: l = s.random_sublist_of_size(10, 7)
        sage: len(l)
        7
        sage: l_check = [-1] + l + [10]
        sage: all(l_check[i] < l_check[i+1] for i in range(8))
        True
    """
    if n > k:
        raise ValueError("n must be <= len(v)")
    if n == k:
        return list(range(k))
    if n >= k // 2 + 5:
        # use complement
        w = random_sublist_of_size(k, k - n)
        m = set(w)
        w = [z for z in range(k) if z not in m]
        assert len(w) == n
        return w

    randrange = current_randstate().python_random().randrange

    w = set()
    while len(w) < n:
        z = randrange(k)
        if z not in w:
            w.add(z)
    return sorted(w)


def solve_system_with_difficult_last_row(B, A):
    """
    Solve the matrix equation ``B*Z = A`` when the last row of `B`
    contains huge entries.

    INPUT:

    - ``B`` -- a square n x n nonsingular matrix with painful big bottom row
    - ``A`` -- an n x k matrix

    OUTPUT: the unique solution to ``B*Z = As``

    EXAMPLES::

        sage: from sage.matrix.matrix_integer_dense_saturation import solve_system_with_difficult_last_row
        sage: B = matrix(ZZ, 3, [1,2,3, 3,-1,2,939239082,39202803080,2939028038402834]); A = matrix(ZZ,3,2,[1,2,4,3,-1,0])
        sage: X = solve_system_with_difficult_last_row(B, A); X
        [  290668794698843/226075992027744         468068726971/409557956572]
        [-226078357385539/1582531944194208       1228691305937/2866905696004]
        [      2365357795/1582531944194208           -17436221/2866905696004]
        sage: B*X == A
        True
    """
    # See the comments in the function of the same name in matrix_integer_dense_hnf.py.
    # This function is just a generalization of that one to A a matrix.
    C = copy(B)
    while True:
        C[C.nrows()-1] = random_matrix(ZZ,1,C.ncols()).row(0)
        try:
            X = C.solve_right(A)
        except ValueError:
            verbose("Try difficult solve again with different random vector")
        else:
            break
    D = B.matrix_from_rows(range(C.nrows()-1))
    N = D._rational_kernel_flint()
    if N.ncols() != 1:
        verbose("Difficult solve quickly failed.  Using direct approach.")
        return B.solve_right(A)

    tm = verbose("Recover correct linear combinations")
    k = N.matrix_from_columns([0])

    # The sought for solution Z to B*Z = A is some linear combination
    #       Z = X + alpha*k
    #  Let w be the last row of B; then Z satisfies
    #       w * Z = A'
    # where A' is the last row of A.  Thus
    #       w * (X + alpha*k) = A'
    # so    w * X + alpha*w*k = A'
    # so    alpha*w*k  = A' - w*X.
    w = B[-1]  # last row of B
    A_prime = A[-1]  # last row of A
    lhs = w*k
    rhs = A_prime - w * X

    if lhs[0] == 0:
        verbose("Difficult solve quickly failed.  Using direct approach.")
        return B.solve_right(A)

    for i in range(X.ncols()):
        alpha = rhs[i] / lhs[0]
        X.set_column(i, (X.matrix_from_columns([i]) + alpha*k).list())
    verbose("Done getting linear combinations.", tm)
    return X

def saturation(A, proof=True, p=0, max_dets=5):
    r"""
    Compute a saturation matrix of `A`.

    INPUT:

    - ``A`` -- a matrix over `\ZZ`
    - ``proof`` -- boolean (default: ``True``)
    - ``p`` -- integer (default: 0); if not 0 only guarantees that output is
      `p`-saturated
    - ``max_dets`` -- integer (default: 4); max number of dets of submatrices to
      compute

    OUTPUT: matrix; saturation of the matrix `A`

    EXAMPLES::

        sage: from sage.matrix.matrix_integer_dense_saturation import saturation
        sage: A = matrix(ZZ, 2, 2, [3,2,3,4]); B = matrix(ZZ, 2,3,[1,2,3,4,5,6]); C = A*B
        sage: C
        [11 16 21]
        [19 26 33]
        sage: C.index_in_saturation()
        18
        sage: S = saturation(C); S
        [11 16 21]
        [-2 -3 -4]
        sage: S.index_in_saturation()
        1
        sage: saturation(C, proof=False)
        [11 16 21]
        [-2 -3 -4]
        sage: saturation(C, p=2)
        [11 16 21]
        [-2 -3 -4]
        sage: saturation(C, p=2, max_dets=1)
        [11 16 21]
        [-2 -3 -4]
    """
    # Find a submatrix of full rank and instead saturate that matrix.
    r = A.rank()
    if A.is_square() and r == A.nrows():
        return identity_matrix(ZZ, r)
    if A.nrows() > r:
        P = []
        while len(P) < r:
            P = matrix_integer_dense_hnf.probable_pivot_rows(A)
        A = A.matrix_from_rows(P)

    # Factor out all common factors from all rows, just in case.
    A = copy(A)
    A._factor_out_common_factors_from_each_row()

    if A.nrows() <= 1:
        return A

    A, zero_cols = A._delete_zero_columns()

    if max_dets > 0:
        # Take the GCD of at most num_dets randomly chosen determinants.
        nr = A.nrows()
        nc = A.ncols()
        d = 0
        trials = min(binomial(nc, nr), max_dets)
        already_tried = []
        while len(already_tried) < trials:
            v = random_sublist_of_size(nc, nr)
            tm = verbose('saturation -- checking det condition on submatrix')
            d = gcd(d, A.matrix_from_columns(v).determinant(proof=proof))
            verbose('saturation -- got det down to %s' % d, tm)
            if gcd(d, p) == 1:
                return A._insert_zero_columns(zero_cols)
            already_tried.append(v)

        if gcd(d, p) == 1:
            # already p-saturated
            return A._insert_zero_columns(zero_cols)

        # Factor and p-saturate at each p.
        # This is not a good algorithm, because all the HNF's in it are really slow!
        #
        #tm = verbose('factoring gcd %s of determinants'%d)
        #limit = 2**31-1
        #F = d.factor(limit = limit)
        #D = [p for p, e in F if p <= limit]
        #B = [n for n, e in F if n > limit]  # all big factors -- there will only be at most one
        #assert len(B) <= 1
        #C = B[0]
        #for p in D:
        #    A = p_saturation(A, p=p, proof=proof)

    # This is a really simple but powerful algorithm.
    # FACT: If A is a matrix of full rank, then hnf(transpose(A))^(-1)*A is a saturation of A.
    # To make this practical we use solve_system_with_difficult_last_row, since the
    # last column of HNF's are typically the only really big ones.
    B = A.transpose().hermite_form(include_zero_rows=False, proof=proof)
    B = B.transpose()

    # Now compute B^(-1) * A
    C = solve_system_with_difficult_last_row(B, A)
    return C.change_ring(ZZ)._insert_zero_columns(zero_cols)

def index_in_saturation(A, proof=True):
    r"""
    The index of A in its saturation.

    INPUT:

    - ``A`` -- matrix over `\ZZ`

    - ``proof`` -- boolean (``True`` or ``False``)

    OUTPUT: integer

    EXAMPLES::

        sage: from sage.matrix.matrix_integer_dense_saturation import index_in_saturation
        sage: A = matrix(ZZ, 2, 2, [3,2,3,4]); B = matrix(ZZ, 2,3,[1,2,3,4,5,6]); C = A*B; C
        [11 16 21]
        [19 26 33]
        sage: index_in_saturation(C)
        18
        sage: W = C.row_space()
        sage: S = W.saturation()
        sage: W.index_in(S)
        18

    For any zero matrix the index in its saturation is 1 (see :issue:`13034`)::

        sage: m = matrix(ZZ, 3)
        sage: m
        [0 0 0]
        [0 0 0]
        [0 0 0]
        sage: m.index_in_saturation()
        1
        sage: m = matrix(ZZ, 2, 3)
        sage: m
        [0 0 0]
        [0 0 0]
        sage: m.index_in_saturation()
        1

    TESTS::

        sage: zero = matrix(ZZ, [[]])
        sage: zero.index_in_saturation()
        1
    """
    r = A.rank()
    if r == 0:
        return ZZ.one()
    if r < A.nrows():
        A = A.hermite_form(proof=proof, include_zero_rows=False)
    if A.is_square():
        return abs(A.determinant(proof=proof))
    A = A.transpose()
    A = A.hermite_form(proof=proof, include_zero_rows=False)
    return abs(A.determinant(proof=proof))
