"""
Modular algorithm to compute Hermite normal forms of integer matrices

AUTHORS:

- Clement Pernet and William Stein (2008-02-07): initial version
"""

from copy import copy

from sage.arith.misc import CRT_list, previous_prime
from sage.matrix.constructor import identity_matrix, matrix, random_matrix
from sage.misc.timing import cputime
from sage.misc.verbose import verbose
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR


def max_det_prime(n):
    """
    Return the largest prime so that it is reasonably efficient to
    compute modulo that prime with n x n matrices in LinBox.

    INPUT:

    - ``n`` -- positive integer

    OUTPUT: a prime number

    EXAMPLES::

        sage: from sage.matrix.matrix_integer_dense_hnf import max_det_prime
        sage: max_det_prime(10000)
        8388593
        sage: max_det_prime(1000)
        8388593
        sage: max_det_prime(10)
        8388593
    """
    # See #14032: LinBox now uses a constant bound of 2^23.
    # This is the largest prime less than that bound.
    return Integer(8388593)


def det_from_modp_and_divisor(A, d, p, z_mod, moduli, z_so_far=ZZ(1), N_so_far=ZZ(1)):
    """
    This is used for internal purposes for computing determinants
    quickly (with the hybrid `p`-adic / multimodular algorithm).

    INPUT:

    - ``A`` -- a square matrix
    - ``d`` -- a divisor of the determinant of A
    - ``p`` -- a prime
    - ``z_mod`` -- values of det/d (mod ...)
    - ``moduli`` -- the moduli so far
    - ``z_so_far`` -- for a modulus p in the list moduli,
      (z_so_far mod p) is the determinant of A modulo p
    - ``N_so_far`` -- N_so_far is the product over the primes in the list moduli

    OUTPUT:

    - A triple (det bound, new z_so_far, new N_so_far).

    EXAMPLES::

        sage: a = matrix(ZZ, 3, [6, 1, 2, -56, -2, -1, -11, 2, -3])
        sage: factor(a.det())
        -1 * 13 * 29
        sage: d = 13
        sage: import sage.matrix.matrix_integer_dense_hnf as matrix_integer_dense_hnf
        sage: matrix_integer_dense_hnf.det_from_modp_and_divisor(a, d, 97, [], [])
        (-377, -29, 97)
        sage: a.det()
        -377
    """
    tm = verbose("Multimodular stage of det calculation -- using p = %s" % p, level=2)
    z = A.mod(p).det() / d
    z = z.lift()
    z_mod.append(z)
    moduli.append(p)
    z = CRT_list([z_so_far, z], [N_so_far, p])
    N = N_so_far*p

    if z > N // 2:
        z -= N
    verbose("Finished multimodular det for p = %s" % p, tm, level=2)
    return (d * z, z, N)


def det_given_divisor(A, d, proof=True, stabilize=2):
    """
    Given a divisor d of the determinant of A, compute the determinant of A.

    INPUT:

    - ``A`` -- square integer matrix
    - ``d`` -- nonzero integer that is assumed to divide the determinant of A
    - ``proof`` -- boolean (default: ``True``); compute det modulo enough primes
      so that the determinant is computed provably correctly (via the
      Hadamard bound).  It would be VERY hard for ``det()`` to fail even
      when ``proof`` is ``False``.
    - ``stabilize`` -- integer (default: 2); if proof = False, then compute
      the determinant modulo `p` until ``stabilize`` successive modulo
      determinant computations stabilize.

    OUTPUT: integer; determinant

    EXAMPLES::

        sage: import sage.matrix.matrix_integer_dense_hnf as matrix_integer_dense_hnf
        sage: a = matrix(ZZ,3,[-1, -1, -1, -20, 4, 1, -1, 1, 2])
        sage: matrix_integer_dense_hnf.det_given_divisor(a, 3)
        -30
        sage: matrix_integer_dense_hnf.det_given_divisor(a, 3, proof=False)
        -30
        sage: matrix_integer_dense_hnf.det_given_divisor(a, 3, proof=False, stabilize=1)
        -30
        sage: a.det()
        -30

    Here we illustrate proof=False giving a wrong answer::

        sage: p = matrix_integer_dense_hnf.max_det_prime(2)
        sage: q = previous_prime(p)
        sage: a = matrix(ZZ, 2, [p, 0, 0, q])
        sage: p * q
        70368442188091
        sage: matrix_integer_dense_hnf.det_given_divisor(a, 1, proof=False, stabilize=2)
        0

    This still works, because we do not work modulo primes that divide
    the determinant bound, which is found using a `p`-adic algorithm::

        sage: a.det(proof=False, stabilize=2)
        70368442188091

    3 primes is enough::

        sage: matrix_integer_dense_hnf.det_given_divisor(a, 1, proof=False, stabilize=3)
        70368442188091
        sage: matrix_integer_dense_hnf.det_given_divisor(a, 1, proof=False, stabilize=5)
        70368442188091
        sage: matrix_integer_dense_hnf.det_given_divisor(a, 1, proof=True)
        70368442188091

    TESTS::

        sage: m = diagonal_matrix(ZZ, 68, [2]*66 + [1,1])
        sage: m.det()
        73786976294838206464
    """
    p = max_det_prime(A.nrows())
    z_mod = []
    moduli = []
    assert d != 0
    z_so_far = 1
    N_so_far = 1
    if proof:
        N = 1
        B = (2 * 10**A.hadamard_bound()) // d + 1
        dd = d
        # bad verbose statement, since computing the log overflows!
        est = int(RR(B).log() / RR(p).log()) + 1
        cnt = 1
        verbose("Multimodular det -- need to use about %s primes." % est,
                level=1)
        while N < B:
            if d % p != 0:
                tm = cputime()
                dd, z_so_far, N_so_far = det_from_modp_and_divisor(A, d, p, z_mod, moduli, z_so_far, N_so_far)
                N *= p
                verbose("computed det mod p=%s which is %s (of about %s)" % (p, cnt, est), tm)
            p = previous_prime(p)
            cnt += 1
        return dd
    else:
        val = []
        while True:
            if d % p:
                tm = cputime()
                dd, z_so_far, N_so_far = det_from_modp_and_divisor(A, d, p, z_mod, moduli, z_so_far, N_so_far)
                verbose("computed det mod %s" % p, tm)
                val.append(dd)
                if len(val) >= stabilize and len(set(val[-stabilize:])) == 1:
                    return val[-1]
            p = previous_prime(p)


def det_padic(A, proof=True, stabilize=2):
    """
    Return the determinant of A, computed using a `p`-adic/multimodular
    algorithm.

    INPUT:

    - ``A`` -- a square matrix

    - ``proof`` -- boolean

    - ``stabilize`` -- (default: 2) if proof False, number of successive primes so that
      CRT det must stabilize

    EXAMPLES::

        sage: import sage.matrix.matrix_integer_dense_hnf as h
        sage: a = matrix(ZZ, 3, [1..9])
        sage: h.det_padic(a)
        0
        sage: a = matrix(ZZ, 3, [1,2,5,-7,8,10,192,5,18])
        sage: h.det_padic(a)
        -3669
        sage: a.determinant(algorithm='ntl')
        -3669
    """
    if not A.is_square():
        raise ValueError("A must be a square matrix")
    r = A.rank()
    if r < A.nrows():
        return ZZ.zero()
    v = random_matrix(ZZ, A.nrows(), 1)
    d = A._solve_right_nonsingular_square(v, check_rank=False).denominator()
    return det_given_divisor(A, d, proof=proof, stabilize=stabilize)


def double_det(A, b, c, proof):
    """
    Compute the determinants of the stacked integer matrices
    A.stack(b) and A.stack(c).

    INPUT:

    - ``A`` -- an (n-1) x n matrix
    - ``b`` -- a 1 x n matrix
    - ``c`` -- a 1 x n matrix
    - ``proof`` -- whether or not to compute the det modulo enough times to
      provably compute the determinant

    OUTPUT: a pair of two integers

    EXAMPLES::

        sage: from sage.matrix.matrix_integer_dense_hnf import double_det
        sage: A = matrix(ZZ, 2, 3, [1,2,3, 4,-2,5])
        sage: b = matrix(ZZ, 1, 3, [1,-2,5])
        sage: c = matrix(ZZ, 1, 3, [8,2,10])
        sage: A.stack(b).det()
        -48
        sage: A.stack(c).det()
        42
        sage: double_det(A, b, c, False)
        (-48, 42)
    """
    # We use the "two for the price of one" algorithm, which I made up. (William Stein)

    # This is a clever trick!  First we transpose everything.  Then
    # we use that if [A|b]*v = c then [A|c]*w = b with w easy to write down!
    # In fact w is got from v by dividing all entries by -v[n], where n is the
    # number of rows of v, and *also* dividing the last entry of w by v[n] again.
    # See this as an algebra exercise where you have to think of matrix vector
    # multiply as "linear combination of columns".
    A = A.transpose()
    b = b.transpose()
    c = c.transpose()
    t = verbose('starting double det')
    B = A.augment(b)
    v = B.solve_right(-c)

    db = det_given_divisor(B, v.denominator(), proof=proof)

    n = v.nrows()
    vn = v[n - 1, 0]
    w = (-1 / vn) * v
    w[n - 1] = w[n - 1] / vn
    dc = det_given_divisor(A.augment(c), w.denominator(), proof=proof)

    verbose('finished double det', t)

    return (db, dc)


def add_column_fallback(B, a, proof):
    """
    Simplistic version of add_column, in case the powerful clever one
    fails (e.g., B is singular).

    INPUT:

    - ``B`` -- a square matrix (may be singular)
    - ``a`` -- an n x 1 matrix, where B has n rows
    - ``proof`` -- boolean; whether to prove result correct

    OUTPUT: x; a vector such that ``H' = H_B.augment(x)`` is the HNF of
    ``A = B.augment(a)``

    EXAMPLES::

        sage: B = matrix(ZZ,3, [-1, -1, 1, -3, 8, -2, -1, -1, -1])
        sage: a = matrix(ZZ,3,1, [1,2,3])
        sage: import sage.matrix.matrix_integer_dense_hnf as matrix_integer_dense_hnf
        sage: matrix_integer_dense_hnf.add_column_fallback(B, a, True)
        [-3]
        [-7]
        [-2]
        sage: matrix_integer_dense_hnf.add_column_fallback(B, a, False)
        [-3]
        [-7]
        [-2]
        sage: B.augment(a).hermite_form()
        [ 1  1  1 -3]
        [ 0 11  1 -7]
        [ 0  0  2 -2]
    """
    tt = verbose('add column fallback...')
    W = B.augment(matrix(ZZ,B.nrows(),a.list()))
    H, _ = hnf(W, proof)
    C = H.matrix_from_columns([H.ncols()-1])
    verbose('finished add column fallback', tt)
    return C


def solve_system_with_difficult_last_row(B, a):
    """
    Solve ``B*x = a`` when the last row of `B` contains huge entries using
    a clever trick that reduces the problem to solve ``C*x = a`` where `C`
    is `B` but with the last row replaced by something small, along
    with one easy null space computation.  The latter are both solved
    `p`-adically.

    INPUT:

    - ``B`` -- a square n x n nonsingular matrix with painful big bottom row
    - ``a`` -- an n x 1 column matrix

    OUTPUT: the unique solution to ``B*x = a``

    EXAMPLES::

        sage: from sage.matrix.matrix_integer_dense_hnf import solve_system_with_difficult_last_row
        sage: B = matrix(ZZ, 3, [1,2,4, 3,-4,7, 939082,2930982,132902384098234])
        sage: a = matrix(ZZ,3,1, [1,2,5])
        sage: z = solve_system_with_difficult_last_row(B, a)
        sage: z
        [ 106321906985474/132902379815497]
        [132902385037291/1329023798154970]
        [        -5221794/664511899077485]
        sage: B*z
        [1]
        [2]
        [5]
    """
    # Here's how:
    # 1. We make a copy of B but with the last *nasty* row of B replaced
    #    by a random very nice row.
    C = copy(B)
    while True:
        C[C.nrows()-1] = random_matrix(ZZ,1,C.ncols()).row(0)
        # 2. Then we find the unique solution to C * x = a
        try:
            x = C.solve_right(a)
        except ValueError:
            verbose("Try difficult solve again with different random vector")
        else:
            break

    # 3. We next delete the last row of B and find a basis vector k
    #    for the 1-dimensional kernel.
    D = B.matrix_from_rows(range(C.nrows()-1))
    N = D._rational_kernel_iml()
    if N.ncols() != 1:
        verbose("Try difficult solve again with different random vector")
        return solve_system_with_difficult_last_row(B, a)

    k = N.matrix_from_columns([0])

    # 4. The sought for solution z to B*z = a is some linear combination
    #
    #       z = x + alpha*k
    #
    #  of x and k, where k is the above fixed basis for the kernel of D.
    #  Setting w to be the last row of B, this column vector z satisfies
    #
    #       w * z = a'
    #
    # where a' is the last entry of a.  Thus
    #
    #       w * (x + alpha*k) = a'
    #
    # so    w * x + alpha*w*k = a'
    # so    alpha*w*k  = a' - w*x.

    w = B[-1]  # last row of B
    a_prime = a[-1]
    lhs = w*k
    rhs = a_prime - w * x

    if lhs[0] == 0:
        verbose("Try difficult solve again with different random vector")
        return solve_system_with_difficult_last_row(B, a)

    alpha = rhs[0] / lhs[0]
    z = x + alpha*k
    return z


def add_column(B, H_B, a, proof):
    """
    The add column procedure.

    INPUT:

    - ``B`` -- a square matrix (may be singular)
    - ``H_B`` -- the Hermite normal form of B
    - ``a`` -- an n x 1 matrix, where B has n rows
    - ``proof`` -- boolean; whether to prove result correct, in case we use fallback method

    OUTPUT:

    - x -- a vector such that H' = H_B.augment(x) is the HNF of A = B.augment(a)

    EXAMPLES::

        sage: B = matrix(ZZ, 3, 3, [1,2,5, 0,-5,3, 1,1,2])
        sage: H_B = B.echelon_form()
        sage: a = matrix(ZZ, 3, 1, [1,8,-2])
        sage: import sage.matrix.matrix_integer_dense_hnf as hnf
        sage: x = hnf.add_column(B, H_B, a, True); x
        [18]
        [ 3]
        [23]
        sage: H_B.augment(x)
        [ 1  0 17 18]
        [ 0  1  3  3]
        [ 0  0 18 23]
        sage: B.augment(a).echelon_form()
        [ 1  0 17 18]
        [ 0  1  3  3]
        [ 0  0 18 23]
    """
    verbose('starting add_column')

    if B.rank() < B.nrows():
        return add_column_fallback(B, a, proof)
    else:
        z = solve_system_with_difficult_last_row(B, a)

    zd, d = z._clear_denom()
    x = H_B * zd
    if d != 1:
        for i in range(x.nrows()):
            x[i, 0] = x[i, 0] / d

    return x


def add_row(A, b, pivots, include_zero_rows):
    """
    The add row procedure.

    INPUT:

    - ``A`` -- a matrix in Hermite normal form with n column
    - ``b`` -- an n x 1 row matrix
    - ``pivots`` -- sorted list of integers; the pivot positions of A

    OUTPUT:

    - ``H`` -- the Hermite normal form of A.stack(b)
    - ``new_pivots`` -- the pivot columns of H

    EXAMPLES::

        sage: import sage.matrix.matrix_integer_dense_hnf as hnf
        sage: A = matrix(ZZ, 2, 3, [-21, -7, 5, 1,20,-7])
        sage: b = matrix(ZZ, 1,3, [-1,1,-1])
        sage: hnf.add_row(A, b, A.pivots(), True)
        (
        [ 1  6 29]
        [ 0  7 28]
        [ 0  0 46], [0, 1, 2]
        )
        sage: A.stack(b).echelon_form()
        [ 1  6 29]
        [ 0  7 28]
        [ 0  0 46]
    """
    t = verbose('add hnf row')
    H, pivs = A._add_row_and_maintain_echelon_form(b.row(0), pivots)
    if include_zero_rows and H.nrows() != A.nrows() + 1:
        H = H.matrix_from_rows(range(A.nrows() + 1))
    verbose('finished add hnf row', t)
    return H, pivs


def pivots_of_hnf_matrix(H):
    """
    Return the pivot columns of a matrix H assumed to be in HNF.

    INPUT:

    - ``H`` -- a matrix that must be HNF

    OUTPUT: list of pivots

    EXAMPLES::

        sage: H = matrix(ZZ, 3, 5, [1, 0, 0, 45, -36, 0, 1, 0, 131, -107, 0, 0, 0, 178, -145]); H
        [   1    0    0   45  -36]
        [   0    1    0  131 -107]
        [   0    0    0  178 -145]
        sage: import sage.matrix.matrix_integer_dense_hnf as matrix_integer_dense_hnf
        sage: matrix_integer_dense_hnf.pivots_of_hnf_matrix(H)
        [0, 1, 3]
    """
    pivots = []
    r = -1
    for j in range(H.ncols()):
        # Find first nonzero position (counting from bottom) in the j-th column
        for i in reversed(range(H.nrows())):
            if H[i, j]:
                if i > r:
                    pivots.append(j)
                    r = i
                else:
                    break
    return pivots


def hnf_square(A, proof):
    """
    INPUT:

    - ``A`` -- a nonsingular n x n matrix over the integers

    OUTPUT: the Hermite normal form of A

    EXAMPLES::

        sage: import sage.matrix.matrix_integer_dense_hnf as hnf
        sage: A = matrix(ZZ, 3, [-21, -7, 5, 1,20,-7, -1,1,-1])
        sage: hnf.hnf_square(A, False)
        [ 1  6 29]
        [ 0  7 28]
        [ 0  0 46]
        sage: A.echelon_form()
        [ 1  6 29]
        [ 0  7 28]
        [ 0  0 46]
    """
    n = A.nrows()
    m = A.ncols()
    if n != m:
        raise ValueError("A must be square.")

    # Small cases -- do not use this algorithm
    if n <= 3:
        return A.echelon_form(algorithm='pari')

    if A.rank() < A.nrows():
        raise ValueError("matrix must have full rank")

    t = verbose("starting slicings")
    B = A.matrix_from_rows(range(m-2)).matrix_from_columns(range(n-1))
    c = A.matrix_from_rows([m-2]).matrix_from_columns(range(n-1))
    d = A.matrix_from_rows([m-1]).matrix_from_columns(range(n-1))
    b = A.matrix_from_columns([n-1]).matrix_from_rows(range(m-2))
    verbose("done slicing", t)

    try:
        d1, d2 = double_det(B, c, d, proof=proof)
    except (ValueError, ZeroDivisionError):
        d1 = B.stack(c).det(proof=proof)
        d2 = B.stack(d).det(proof=proof)

    g, k, l = d1._xgcd(d2, minimal=True)

    W = B.stack(k * c + l * d)
    verbose("submatrix det: g=%s" % g)
    CUTOFF = 2147483647  # 2^31-1
    if g == 0:
        # Big trouble -- matrix is not invertible
        # Since we have no good conditioning code at present,
        # in this case we just fall back to using pari.
        H = W.echelon_form(algorithm='pari')
    elif 2 * g > CUTOFF:
        # Unlikely that g will be large on even slightly random input
        # if it is, we fallback to the traditional algorithm.
        # A nasty example is A = n*random_matrix(ZZ,m), where
        # this algorithm gets killed.  This is not random input though.
        f = W.gcd()
        g = g / (f**W.nrows())
        if 2 * g <= CUTOFF:
            verbose("Found common factor of %s -- dividing out; get new g = %s" % (f, g))
            W0 = (W / f).change_ring(ZZ)
            H = W0._hnf_mod(2 * g)
            H *= f
        else:
            verbose("Falling back to PARI HNF since input matrix is ill conditioned for p-adic hnf algorithm.")
            # We need more clever preconditioning?
            # It is important to *not* just do the submatrix, since
            # the whole rest of the algorithm will likely be very slow in
            # weird cases where the det is large.
            # E.g., matrix all of whose rows but 1 are multiplied by some
            # fixed scalar n.
            raise NotImplementedError("fallback to PARI")
            # H = W.hermite_form(algorithm='pari')
    else:
        H = W._hnf_mod(2 * g)

    x = add_column(W, H, b.stack(matrix(1, 1,
                                        [k*A[m-2,m-1] + l*A[m-1,m-1]])),
                   proof)
    Hprime = H.augment(x)
    pivots = pivots_of_hnf_matrix(Hprime)

    Hprime, pivots = add_row(Hprime, A.matrix_from_rows([m - 2]),
                             pivots, include_zero_rows=False)
    Hprime, pivots = add_row(Hprime, A.matrix_from_rows([m - 1]),
                             pivots, include_zero_rows=False)
    return Hprime.matrix_from_rows(range(m))


def interleave_matrices(A, B, cols1, cols2):
    """
    INPUT:

    - A, B -- matrices with the same number of rows
    - cols1, cols2 -- disjoint lists of integers

    OUTPUT:

    construct a new matrix C by sticking the columns
    of A at the positions specified by cols1 and the
    columns of B at the positions specified by cols2.

    EXAMPLES::

        sage: A = matrix(ZZ, 2, [1,2,3,4]); B = matrix(ZZ, 2, [-1,5,2,3])
        sage: A
        [1 2]
        [3 4]
        sage: B
        [-1  5]
        [ 2  3]
        sage: import sage.matrix.matrix_integer_dense_hnf as hnf
        sage: hnf.interleave_matrices(A, B, [1,3], [0,2])
        [-1  1  5  2]
        [ 2  3  3  4]
    """
    D = A.augment(B)
    w = cols1 + cols2
    v = [w.index(i) for i in range(len(cols1) + len(cols2))]
    return D.matrix_from_columns(v)


def probable_pivot_rows(A):
    """
    Return rows of A that are very likely to be pivots.

    This really finds the pivots of A modulo a random prime.

    INPUT:

    - ``A`` -- a matrix

    OUTPUT: a tuple of integers

    EXAMPLES::

        sage: import sage.matrix.matrix_integer_dense_hnf as matrix_integer_dense_hnf
        sage: a = matrix(ZZ,3,[0, -1, -1, 0, -20, 1, 0, 1, 2])
        sage: a
        [  0  -1  -1]
        [  0 -20   1]
        [  0   1   2]
        sage: matrix_integer_dense_hnf.probable_pivot_rows(a)
        (0, 1)
    """
    return probable_pivot_columns(A.transpose())


def probable_pivot_columns(A):
    """
    INPUT:

    - ``A`` -- a matrix

    OUTPUT: a tuple of integers

    EXAMPLES::

        sage: import sage.matrix.matrix_integer_dense_hnf as matrix_integer_dense_hnf
        sage: a = matrix(ZZ,3,[0, -1, -1, 0, -20, 1, 0, 1, 2])
        sage: a
        [  0  -1  -1]
        [  0 -20   1]
        [  0   1   2]
        sage: matrix_integer_dense_hnf.probable_pivot_columns(a)
        (1, 2)
    """
    p = ZZ.random_element(10007, 46000).next_prime()
    return A._reduce(p).pivots()


def ones(H, pivots):
    """
    Find all 1 pivot columns of the matrix H in Hermite form, along
    with the corresponding rows, and also the non 1 pivot columns and
    non-pivot rows.  Here a 1 pivot column is a pivot column so that
    the leading bottom entry is 1.

    INPUT:

    - ``H`` -- matrix in Hermite form
    - ``pivots`` -- list of integers (all pivot positions of H)

    OUTPUT:

    4-tuple of integer lists: onecol, onerow, non_oneol, non_onerow

    EXAMPLES::

        sage: H = matrix(ZZ, 3, 5, [1, 0, 0, 45, -36, 0, 1, 0, 131, -107, 0, 0, 0, 178, -145]); H
        [   1    0    0   45  -36]
        [   0    1    0  131 -107]
        [   0    0    0  178 -145]
        sage: import sage.matrix.matrix_integer_dense_hnf as matrix_integer_dense_hnf
        sage: matrix_integer_dense_hnf.ones(H, [0,1,3])
        ([0, 1], [0, 1], [2], [2])
    """
    # Find the "onecol" pivot columns of H, i.e., the columns
    # that contain exactly one "1" entry and all other entries 0.
    onecol = []
    onerow = []
    i = 0
    for c in pivots:
        if H[i, c] == 1:
            onecol.append(c)
            onerow.append(i)
        i += 1
    onecol_set = set(onecol)
    non_onerow = [j for j in range(len(pivots)) if j not in onerow]
    non_onecol = [j for j in range(H.ncols()) if j not in onecol_set][:len(non_onerow)]
    return onecol, onerow, non_onecol, non_onerow


def extract_ones_data(H, pivots):
    """
    Compute ones data and corresponding submatrices of H.

    This is used to optimized the :func:`add_row` function.

    INPUT:

    - ``H`` -- a matrix in HNF
    - ``pivots`` -- list of all pivot column positions of H

    OUTPUT:

    C, D, E, onecol, onerow, non_onecol, non_onerow
    where onecol, onerow, non_onecol, non_onerow are as for
    the ones function, and C, D, E are matrices:

    - ``C`` -- submatrix of all non-onecol columns and onecol rows
    - ``D`` -- all non-onecol columns and other rows
    - ``E`` -- inverse of D

    If D is not invertible or there are 0 or more than 2 non onecols,
    then C, D, and E are set to None.

    EXAMPLES::

        sage: H = matrix(ZZ, 3, 4, [1, 0, 0, 7, 0, 1, 5, 2, 0, 0, 6, 6])
        sage: import sage.matrix.matrix_integer_dense_hnf as matrix_integer_dense_hnf
        sage: matrix_integer_dense_hnf.extract_ones_data(H, [0,1,2])
        (
        [0]
        [5], [6], [1/6], [0, 1], [0, 1], [2], [2]
        )

    Here we get None's since the (2,2) position submatrix is not invertible.
        sage: H = matrix(ZZ, 3, 5, [1, 0, 0, 45, -36, 0, 1, 0, 131, -107, 0, 0, 0, 178, -145]); H
        [   1    0    0   45  -36]
        [   0    1    0  131 -107]
        [   0    0    0  178 -145]
        sage: import sage.matrix.matrix_integer_dense_hnf as matrix_integer_dense_hnf
        sage: matrix_integer_dense_hnf.extract_ones_data(H, [0,1,3])
        (None, None, None, [0, 1], [0, 1], [2], [2])
    """
    onecol, onerow, non_onecol, non_onerow = ones(H, pivots)
    verbose('extract_ones -- got submatrix of size %s' % len(non_onecol))
    if len(non_onecol) in [1, 2]:
        # Extract submatrix of all non-onecol columns and onecol rows
        C = H.matrix_from_rows_and_columns(onerow, non_onecol)
        # Extract submatrix of all non-onecol columns and other rows
        D = H.matrix_from_rows_and_columns(non_onerow, non_onecol).transpose()
        tt = verbose("extract ones -- INVERT %s x %s" % (len(non_onerow), len(non_onecol)), level=1)
        try:
            E = D**(-1)
        except ZeroDivisionError:
            C = D = E = None
        verbose("done inverting", tt, level=1)
        return C, D, E, onecol, onerow, non_onecol, non_onerow
    else:
        return None, None, None, onecol, onerow, non_onecol, non_onerow


def is_in_hnf_form(H, pivots):
    """
    Return whether the matrix ``H`` is in Hermite normal form
    with given pivot columns.

    INPUT:

    - ``H`` -- matrix
    - ``pivots`` -- sorted list of integers

    OUTPUT: boolean

    EXAMPLES::

        sage: a = matrix(ZZ,3,5,[-2, -6, -3, -17, -1, 2, -1, -1, -2, -1, -2, -2, -6, 9, 2])
        sage: import sage.matrix.matrix_integer_dense_hnf as matrix_integer_dense_hnf
        sage: matrix_integer_dense_hnf.is_in_hnf_form(a,range(3))
        False
        sage: e = a.hermite_form(); p = a.pivots()
        sage: matrix_integer_dense_hnf.is_in_hnf_form(e, p)
        True
    """
    tt = verbose('testing if matrix is in HNF')
    r = 0
    pivots_set = set(pivots)
    for j in range(H.ncols()):
        if j in pivots_set:
            for i in range(r + 1, H.nrows()):
                if H[i, j]:
                    verbose('not HNF because nonzeros below pivot position',tt)
                    return False
            for i in range(r):
                if H[i, j] < 0 or H[i, j] >= H[r, j]:
                    verbose('not HNF because negative or too big above pivot position',tt)
                    return False
            r += 1
        else:
            for i in range(r, H.nrows()):
                if H[i, j]:
                    verbose('not HNF nonzero in wrong place in nonpivot column',tt)
                    return False
    verbose('done verifying in HNF -- yes', tt)
    return True


def probable_hnf(A, include_zero_rows, proof):
    """
    Return the HNF of A or raise an exception if something involving
    the randomized nature of the algorithm goes wrong along the way.

    Calling this function again a few times should result it in it
    working, at least if proof=True.

    INPUT:

    - ``A`` -- a matrix
    - ``include_zero_rows`` -- boolean
    - ``proof`` -- boolean

    OUTPUT:

    the Hermite normal form of A.
    cols -- pivot columns

    EXAMPLES::

        sage: a = matrix(ZZ,4,3,[-1, -1, -1, -20, 4, 1, -1, 1, 2,1,2,3])
        sage: import sage.matrix.matrix_integer_dense_hnf as matrix_integer_dense_hnf
        sage: matrix_integer_dense_hnf.probable_hnf(a, True, True)
        (
        [1 0 0]
        [0 1 0]
        [0 0 1]
        [0 0 0], [0, 1, 2]
        )
        sage: matrix_integer_dense_hnf.probable_hnf(a, False, True)
        (
        [1 0 0]
        [0 1 0]
        [0 0 1], [0, 1, 2]
        )
        sage: matrix_integer_dense_hnf.probable_hnf(a, False, False)
        (
        [1 0 0]
        [0 1 0]
        [0 0 1], [0, 1, 2]
        )
    """
    # Find left-most full rank submatrix by working modulo a prime
    rows = list(probable_pivot_rows(A))
    B = A.matrix_from_rows(rows)
    cols = list(probable_pivot_columns(B))
    C = B.matrix_from_columns(cols)
    # Now C is a submatrix of A that has full rank and is square.

    # We compute the HNF of C, which is a square nonsingular matrix.
    try:
        H = hnf_square(C, proof=proof)
    except NotImplementedError:
        # raise
        # this signals that we must fallback to PARI
        verbose("generic random modular HNF algorithm failed -- we fall back to PARI")
        H = A.hermite_form(algorithm='pari', include_zero_rows=include_zero_rows, proof=proof)
        return H, H.pivots()

    # The transformation matrix to HNF is the unique
    # matrix U such that U * C = H, i.e., U = H*C^(-1).

    if len(cols) < B.ncols():
        # We compute the HNF of B by multiplying the matrix D
        # got from the columns not in C by U:
        # We want to compute X = U*D.  But U = H*C^(-1),
        # so X = U*D = H*C^(-1)*D.
        # So C*H^(-1)*X = D

        # find y s.t C*y = D
        #   H^(-1)*X = y ===> X = H*y
        #
        cols_set = set(cols)
        cols2 = [i for i in range(B.ncols()) if i not in cols_set]
        D = B.matrix_from_columns(cols2)
        Y = C.solve_right(D)
        H2 = H * Y
        H2 = H2.change_ring(ZZ)

        # The HNF of B is got by assembling together
        # the matrices H and H2.
        H = interleave_matrices(H, H2, cols, cols2)

    pivots = pivots_of_hnf_matrix(H)

    # Now H is the HNF of the matrix B.
    # Finally we add all remaining rows of A to H using
    # the add_row function.

    C, D, E, onecol, onerow, non_onecol, non_onerow = extract_ones_data(H, cols)
    if not proof and len(non_onecol) == 0:
        # Identity matrix -- done
        verbose("hnf -- got identity matrix -- early abort (0)")
        if include_zero_rows:
            H = pad_zeros(H, A.nrows())
        return H, pivots

    rows_set = set(rows)
    for i in range(A.nrows()):
        if i not in rows_set:
            v = A.matrix_from_rows([i])
            if v == 0:
                continue
            if E is None:
                H, pivots = add_row(H, v, pivots, include_zero_rows=False)
                C, D, E, onecol, onerow, non_onecol, non_onerow = extract_ones_data(H, pivots)
                if not proof and len(non_onecol) == 0:
                    # Identity matrix -- done
                    verbose("hnf -- got identity matrix -- early abort (1)")
                    if include_zero_rows:
                        H = pad_zeros(H, A.nrows())
                    return H, pivots
            else:
                z = A.matrix_from_rows_and_columns([i], non_onecol)
                w = A.matrix_from_rows_and_columns([i], onecol)
                tt = verbose("checking denom (%s x %s)" % (D.nrows(),
                                                           D.ncols()))
                Y = (z - w * C).transpose()
                k = E * Y
                verbose("done checking denom",tt)
                if k.denominator() != 1:
                    H, pivots = add_row(H, v, pivots, include_zero_rows=False)
                    D = H.matrix_from_rows_and_columns(non_onerow, non_onecol).transpose()
                nn = ones(H, pivots)
                if not proof and len(nn[2]) == 0:
                    verbose("hnf -- got identity matrix -- early abort (2)")
                    if include_zero_rows:
                        H = pad_zeros(H, A.nrows())
                    return H, pivots

    if include_zero_rows:
        H = pad_zeros(H, A.nrows())
    return H, pivots


def pad_zeros(A, nrows):
    """
    Add zeros to the bottom of A so that the resulting matrix has nrows.

    INPUT:

    - ``A`` -- a matrix
    - ``nrows`` -- integer that is at least as big as the number of rows of A

    OUTPUT: a matrix with nrows rows

    EXAMPLES::

        sage: import sage.matrix.matrix_integer_dense_hnf as matrix_integer_dense_hnf
        sage: a = matrix(ZZ, 2, 4, [1, 0, 0, 7, 0, 1, 5, 2])
        sage: matrix_integer_dense_hnf.pad_zeros(a, 4)
        [1 0 0 7]
        [0 1 5 2]
        [0 0 0 0]
        [0 0 0 0]
        sage: matrix_integer_dense_hnf.pad_zeros(a, 2)
        [1 0 0 7]
        [0 1 5 2]
    """
    nz = nrows - A.nrows()
    if nz == 0:
        return A
    if nz < 0:
        return A.matrix_from_rows(range(nrows))
    return A.stack(matrix(ZZ, nz, A.ncols()))


def hnf(A, include_zero_rows=True, proof=True):
    """
    Return the Hermite Normal Form of a general integer matrix A,
    along with the pivot columns.

    INPUT:

    - ``A`` -- an n x m matrix A over the integers
    - ``include_zero_rows`` -- boolean (default: ``True``); whether or not to include zero
      rows in the output matrix
    - ``proof`` -- whether or not to prove the result correct

    OUTPUT: tuple of:

    - ``matrix`` -- the Hermite normal form of A
    - ``pivots`` -- the pivot column positions of A

    EXAMPLES::

        sage: import sage.matrix.matrix_integer_dense_hnf as matrix_integer_dense_hnf
        sage: a = matrix(ZZ,3,5,[-2, -6, -3, -17, -1, 2, -1, -1, -2, -1, -2, -2, -6, 9, 2])
        sage: matrix_integer_dense_hnf.hnf(a)
        (
        [   2    0   26  -75  -10]
        [   0    1   27  -73   -9]
        [   0    0   37 -106  -13], [0, 1, 2]
        )
        sage: matrix_integer_dense_hnf.hnf(a.transpose())
        (
        [1 0 0]
        [0 1 0]
        [0 0 1]
        [0 0 0]
        [0 0 0], [0, 1, 2]
        )
        sage: matrix_integer_dense_hnf.hnf(a.transpose(), include_zero_rows=False)
        (
        [1 0 0]
        [0 1 0]
        [0 0 1], [0, 1, 2]
        )
    """
    if A.nrows() <= 1:
        np = A.nonzero_positions()
        if not np:
            pivots = []
            if not include_zero_rows:
                A = A.new_matrix(0)  # 0 rows
        else:
            i,j = np[0]
            if A[i,j] < 0:
                A = -A
            pivots = [j]
        return A, pivots

    if not proof:
        H, pivots = probable_hnf(A, include_zero_rows=include_zero_rows,
                                 proof=False)
        if not include_zero_rows and len(pivots) > H.nrows():
            return H.matrix_from_rows(range(len(pivots))), pivots

    while True:
        try:
            H, pivots = probable_hnf(A, include_zero_rows=include_zero_rows,
                                 proof=True)
        except ValueError:
            verbose("The attempt failed since the pivots must have been wrong. We try again.")
            continue

        if is_in_hnf_form(H, pivots):
            if not include_zero_rows and len(pivots) > H.nrows():
                H = H.matrix_from_rows(range(len(pivots)))
            return H, pivots
        verbose("After attempt the return matrix is not in HNF form since pivots must have been wrong. We try again.")


def hnf_with_transformation(A, proof=True):
    """
    Compute the HNF H of A along with a transformation matrix U
    such that U*A = H.

    INPUT:

    - ``A`` -- an n x m matrix A over the integers
    - ``proof`` -- whether or not to prove the result correct

    OUTPUT: tuple of:

    - ``matrix`` -- the Hermite normal form H of A
    - ``U`` -- a unimodular matrix such that U * A = H

    EXAMPLES::

        sage: import sage.matrix.matrix_integer_dense_hnf as matrix_integer_dense_hnf
        sage: A = matrix(ZZ, 2, [1, -5, -10, 1, 3, 197]); A
        [  1  -5 -10]
        [  1   3 197]
        sage: H, U = matrix_integer_dense_hnf.hnf_with_transformation(A)
        sage: H
        [  1   3 197]
        [  0   8 207]
        sage: U
        [ 0  1]
        [-1  1]
        sage: U*A
        [  1   3 197]
        [  0   8 207]
    """
    # All we do is augment the input matrix with the identity matrix of the appropriate rank on the right.
    C = A.augment(identity_matrix(ZZ, A.nrows()))
    H, _ = hnf(C, include_zero_rows=True, proof=proof)
    U = H.matrix_from_columns(range(A.ncols(), H.ncols()))
    H2 = H.matrix_from_columns(range(A.ncols()))
    return H2, U


def hnf_with_transformation_tests(n=10, m=5, trials=10):
    """
    Use this to randomly test that hnf with transformation matrix
    is working.

    EXAMPLES::

        sage: from sage.matrix.matrix_integer_dense_hnf import hnf_with_transformation_tests
        sage: hnf_with_transformation_tests(n=15, m=10, trials=10)
        0 1 2 3 4 5 6 7 8 9
    """
    for i in range(trials):
        print(i, end=" ")
        A = random_matrix(ZZ, n, m)
        H, U = hnf_with_transformation(A)
        assert H == U * A
        H, U = hnf_with_transformation(A, proof=False)
        assert H == U * A


# -----------------------------------------------------------------------------
# Code for testing and benchmarking
# ------------------------------------------------------------------------------


def benchmark_hnf(nrange, bits=4):
    """
    Run benchmark program.

    EXAMPLES::

        sage: import sage.matrix.matrix_integer_dense_hnf as hnf
        sage: hnf.benchmark_hnf([10,25],32)
        ('sage', 10, 32, ...),
        ('sage', 25, 32, ...),
    """
    b = 2**bits
    for n in nrange:
        a = random_matrix(ZZ, n, x=-b, y=b)
        t = cputime()
        h, _ = hnf(a, proof=False)
        tm = cputime(t)
        print('%s,' % (('sage', n, bits, tm),))


def benchmark_magma_hnf(nrange, bits=4):
    """
    EXAMPLES::

        sage: import sage.matrix.matrix_integer_dense_hnf as hnf
        sage: hnf.benchmark_magma_hnf([50,100],32)     # optional - magma
        ('magma', 50, 32, ...),
        ('magma', 100, 32, ...),
    """
    from sage.interfaces.magma import magma
    b = 2**bits
    for n in nrange:
        a = magma('MatrixAlgebra(IntegerRing(),%s)![Random(%s,%s) : i in [1..%s]]' % (n, -b, b, n**2))
        t = magma.cputime()
        a.EchelonForm()
        tm = magma.cputime(t)
        print('%s,' % (('magma', n, bits, tm),))


def sanity_checks(times=50, n=8, m=5, proof=True, stabilize=2,
                  check_using_magma=True):
    """
    Run random sanity checks on the modular `p`-adic HNF with tall and wide matrices
    both dense and sparse.

    INPUT:

    - ``times`` -- number of times to randomly try matrices with each shape
    - ``n`` -- number of rows
    - ``m`` -- number of columns
    - ``proof`` -- test with proof true
    - ``stabilize`` -- parameter to pass to hnf algorithm when proof is False
    - ``check_using_magma`` -- if ``True`` use Magma instead of PARI to check
      correctness of computed HNF's. Since PARI's HNF is buggy and slow (as of
      2008-02-16 non-pivot entries sometimes are not normalized to be
      nonnegative) the default is Magma.

    EXAMPLES::

        sage: import sage.matrix.matrix_integer_dense_hnf as matrix_integer_dense_hnf
        sage: matrix_integer_dense_hnf.sanity_checks(times=5, check_using_magma=False)
        small 8 x 5
        0 1 2 3 4  (done)
        big 8 x 5
        0 1 2 3 4  (done)
        small 5 x 8
        0 1 2 3 4  (done)
        big 5 x 8
        0 1 2 3 4  (done)
        sparse 8 x 5
        0 1 2 3 4  (done)
        sparse 5 x 8
        0 1 2 3 4  (done)
        ill conditioned -- 1000*A -- 8 x 5
        0 1 2 3 4  (done)
        ill conditioned -- 1000*A but one row -- 8 x 5
        0 1 2 3 4  (done)
    """
    if check_using_magma:
        from sage.interfaces.magma import magma

    def __do_check(v):
        """
        This is used internally by the sanity check code.
        """
        for i, a in enumerate(v):
            global sanity
            sanity = a
            print(i, end=" ")
            if check_using_magma:
                if magma(hnf(a)[0]) != magma(a).EchelonForm():
                    print("bug computing hnf of a matrix")
                    print('a = matrix(ZZ, %s, %s, %s)' % (a.nrows(), a.ncols(),
                                                          a.list()))
                    return
            else:
                if hnf(a)[0] != a.echelon_form(algorithm='pari'):
                    print("bug computing hnf of a matrix")
                    print('a = matrix(ZZ, %s, %s, %s)' % (a.nrows(), a.ncols(),
                                                          a.list()))
                    return
        print(" (done)")
    print("small %s x %s" % (n, m))
    __do_check([random_matrix(ZZ, n, m, x=-1, y=1) for _ in range(times)])
    print("big %s x %s" % (n, m))
    __do_check([random_matrix(ZZ, n, m, x=-2**32, y=2**32)
                for _ in range(times)])

    print("small %s x %s" % (m, n))
    __do_check([random_matrix(ZZ, m, n, x=-1, y=1) for _ in range(times)])
    print("big %s x %s" % (m, n))
    __do_check([random_matrix(ZZ, m, n, x=-2**32,y=2**32)
                for _ in range(times)])

    print("sparse %s x %s" % (n, m))
    __do_check([random_matrix(ZZ, n, m, density=0.1) for _ in range(times)])
    print("sparse %s x %s" % (m, n))
    __do_check([random_matrix(ZZ, m, n, density=0.1) for _ in range(times)])

    print("ill conditioned -- 1000*A -- %s x %s" % (n, m))
    __do_check([1000*random_matrix(ZZ, n, m, x=-1, y=1) for _ in range(times)])

    print("ill conditioned -- 1000*A but one row -- %s x %s" % (n, m))
    v = []
    for _ in range(times):
        a = 1000 * random_matrix(ZZ, n, m, x=-1, y=1)
        a[a.nrows() - 1] = a[a.nrows() - 1] / 1000
        v.append(a)
    __do_check(v)
