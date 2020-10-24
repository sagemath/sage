"""
Misc matrix algorithms

Code goes here mainly when it needs access to the internal structure
of several classes, and we want to avoid circular cimports.

NOTE: The whole problem of avoiding circular imports -- the reason for
existence of this file -- is now a non-issue, since some bugs in
Cython were fixed.  Probably all this code should be moved into the
relevant classes and this file deleted.
"""

from cysignals.signals cimport sig_check

from sage.ext.mod_int cimport *
from sage.libs.gmp.mpz cimport *
from sage.libs.gmp.mpq cimport *
from sage.libs.mpfr cimport *

from sage.libs.flint.fmpz cimport fmpz_set_mpz, fmpz_one
from sage.libs.flint.fmpq cimport fmpq_set_mpq, fmpq_canonicalise
from sage.libs.flint.fmpq_mat cimport fmpq_mat_entry_num, fmpq_mat_entry_den, fmpq_mat_entry

from sage.arith.rational_reconstruction cimport mpq_rational_reconstruction

from sage.data_structures.binary_search cimport *
from sage.modules.vector_integer_sparse cimport *
from sage.modules.vector_rational_sparse cimport *
from sage.modules.vector_modn_sparse cimport *

from .matrix0 cimport Matrix
from .matrix_integer_dense cimport Matrix_integer_dense
from .matrix_integer_sparse cimport Matrix_integer_sparse
from .matrix_rational_dense cimport Matrix_rational_dense
from .matrix_rational_sparse cimport Matrix_rational_sparse

from sage.rings.integer_ring   import ZZ
from sage.rings.rational_field import QQ

from sage.rings.integer cimport Integer
from sage.arith.all import previous_prime, CRT_basis

from sage.rings.real_mpfr import  is_RealField
from sage.rings.real_mpfr cimport RealNumber


from sage.misc.verbose import verbose, get_verbose

def matrix_integer_dense_rational_reconstruction(Matrix_integer_dense A, Integer N):
    """
    Given a matrix over the integers and an integer modulus, do
    rational reconstruction on all entries of the matrix, viewed as
    numbers mod N.  This is done efficiently by assuming there is a
    large common factor dividing the denominators.

    INPUT:

        A -- matrix
        N -- an integer

    EXAMPLES::

        sage: B = ((matrix(ZZ, 3,4, [1,2,3,-4,7,2,18,3,4,3,4,5])/3)%500).change_ring(ZZ)
        sage: sage.matrix.misc.matrix_integer_dense_rational_reconstruction(B, 500)
        [ 1/3  2/3    1 -4/3]
        [ 7/3  2/3    6    1]
        [ 4/3    1  4/3  5/3]

    TESTS:

    Check that :trac:`9345` is fixed::

        sage: A = random_matrix(ZZ, 3)
        sage: sage.matrix.misc.matrix_integer_dense_rational_reconstruction(A, 0)
        Traceback (most recent call last):
        ...
        ZeroDivisionError: The modulus cannot be zero
    """
    if not N:
        raise ZeroDivisionError("The modulus cannot be zero")
    cdef Matrix_rational_dense R
    R = Matrix_rational_dense.__new__(Matrix_rational_dense,
                                      A.parent().change_ring(QQ), 0,0,0)

    cdef mpz_t a, bnd, other_bnd, denom, tmp
    cdef mpq_t qtmp
    cdef Integer _bnd
    cdef Py_ssize_t i, j
    cdef int do_it

    mpz_init_set_si(denom, 1)
    mpz_init(a)
    mpz_init(tmp)
    mpz_init(other_bnd)
    mpq_init(qtmp)

    _bnd = (N//2).isqrt()
    mpz_init_set(bnd, _bnd.value)
    mpz_sub(other_bnd, N.value, bnd)

    for i in range(A._nrows):
        for j in range(A._ncols):
            sig_check()
            A.get_unsafe_mpz(i, j, a)
            if mpz_cmp_ui(denom, 1) != 0:
                mpz_mul(a, a, denom)
            mpz_fdiv_r(a, a, N.value)
            do_it = 0
            if mpz_cmp(a, bnd) <= 0:
                do_it = 1
            elif mpz_cmp(a, other_bnd) >= 0:
                mpz_sub(a, a, N.value)
                do_it = 1
            if do_it:
                fmpz_set_mpz(fmpq_mat_entry_num(R._matrix, i, j), a)
                if mpz_cmp_ui(denom, 1) != 0:
                    fmpz_set_mpz(fmpq_mat_entry_den(R._matrix, i, j), denom)
                    fmpq_canonicalise(fmpq_mat_entry(R._matrix, i, j))
                else:
                    fmpz_one(fmpq_mat_entry_den(R._matrix, i, j))
            else:
                # Otherwise have to do it the hard way
                A.get_unsafe_mpz(i, j, tmp)
                mpq_rational_reconstruction(qtmp, tmp, N.value)
                mpz_lcm(denom, denom, mpq_denref(qtmp))
                fmpq_set_mpq(fmpq_mat_entry(R._matrix, i, j), qtmp)

    mpz_clear(denom)
    mpz_clear(a)
    mpz_clear(tmp)
    mpz_clear(other_bnd)
    mpz_clear(bnd)
    mpq_clear(qtmp)

    return R


def matrix_integer_sparse_rational_reconstruction(Matrix_integer_sparse A, Integer N):
    """
    Given a sparse matrix over the integers and an integer modulus, do
    rational reconstruction on all entries of the matrix, viewed as
    numbers mod N.

    EXAMPLES::

        sage: A = matrix(ZZ, 3, 4, [(1/3)%500, 2, 3, (-4)%500, 7, 2, 2, 3, 4, 3, 4, (5/7)%500], sparse=True)
        sage: sage.matrix.misc.matrix_integer_sparse_rational_reconstruction(A, 500)
        [1/3   2   3  -4]
        [  7   2   2   3]
        [  4   3   4 5/7]

    TESTS:

    Check that :trac:`9345` is fixed::

        sage: A = random_matrix(ZZ, 3, sparse=True)
        sage: sage.matrix.misc.matrix_integer_sparse_rational_reconstruction(A, 0)
        Traceback (most recent call last):
        ...
        ZeroDivisionError: The modulus cannot be zero
    """
    if not N:
        raise ZeroDivisionError("The modulus cannot be zero")
    cdef Matrix_rational_sparse R
    R = Matrix_rational_sparse.__new__(Matrix_rational_sparse,
                                      A.parent().change_ring(QQ), 0,0,0)

    cdef mpq_t t
    cdef mpz_t a, bnd, other_bnd, denom
    cdef Integer _bnd
    cdef Py_ssize_t i, j
    cdef int do_it
    cdef mpz_vector* A_row
    cdef mpq_vector* R_row

    mpq_init(t)
    mpz_init_set_si(denom, 1)
    mpz_init(a)
    mpz_init(other_bnd)

    _bnd = (N//2).isqrt()
    mpz_init_set(bnd, _bnd.value)
    mpz_sub(other_bnd, N.value, bnd)

    for i in range(A._nrows):
        sig_check()
        A_row = &A._matrix[i]
        R_row = &R._matrix[i]
        reallocate_mpq_vector(R_row, A_row.num_nonzero)
        R_row.num_nonzero = A_row.num_nonzero
        R_row.degree = A_row.degree
        for j in range(A_row.num_nonzero):
            sig_check()
            mpz_set(a, A_row.entries[j])
            if mpz_cmp_ui(denom, 1) != 0:
                mpz_mul(a, a, denom)
            mpz_fdiv_r(a, a, N.value)
            do_it = 0
            if mpz_cmp(a, bnd) <= 0:
                do_it = 1
            elif mpz_cmp(a, other_bnd) >= 0:
                mpz_sub(a, a, N.value)
                do_it = 1
            if do_it:
                mpz_set(mpq_numref(t), a)
                if mpz_cmp_ui(denom, 1) != 0:
                    mpz_set(mpq_denref(t), denom)
                    mpq_canonicalize(t)
                else:
                    mpz_set_si(mpq_denref(t), 1)
                mpq_set(R_row.entries[j], t)
                R_row.positions[j] = A_row.positions[j]
            else:
                # Otherwise have to do it the hard way
                mpq_rational_reconstruction(t, A_row.entries[j], N.value)
                mpq_set(R_row.entries[j], t)
                R_row.positions[j] = A_row.positions[j]
                mpz_lcm(denom, denom, mpq_denref(t))

    mpq_clear(t)
    mpz_clear(denom)
    mpz_clear(a)
    mpz_clear(other_bnd)
    mpz_clear(bnd)

    return R


def matrix_rational_echelon_form_multimodular(Matrix self, height_guess=None, proof=None):
    """
    Returns reduced row-echelon form using a multi-modular
    algorithm.  Does not change self.

    REFERENCE: Chapter 7 of Stein's "Explicitly Computing Modular Forms".

    INPUT:

    - height_guess -- integer or None
    - proof -- boolean or None (default: None, see proof.linear_algebra or
      sage.structure.proof). Note that the global Sage default is proof=True

    OUTPUT: a pair consisting of a matrix in echelon form and a tuple of pivot
    positions.

    ALGORITHM:

    The following is a modular algorithm for computing the echelon
    form.  Define the height of a matrix to be the max of the
    absolute values of the entries.

    Given Matrix A with n columns (self).

     0. Rescale input matrix A to have integer entries.  This does
        not change echelon form and makes reduction modulo lots of
        primes significantly easier if there were denominators.
        Henceforth we assume A has integer entries.

     1. Let c be a guess for the height of the echelon form.  E.g.,
        c=1000, e.g., if matrix is very sparse and application is to
        computing modular symbols.

     2. Let M = n * c * H(A) + 1,
        where n is the number of columns of A.

     3. List primes p_1, p_2, ..., such that the product of
        the p_i is at least M.

     4. Try to compute the rational reconstruction CRT echelon form
        of A mod the product of the p_i.  If rational
        reconstruction fails, compute 1 more echelon forms mod the
        next prime, and attempt again.  Make sure to keep the
        result of CRT on the primes from before, so we don't have
        to do that computation again.  Let E be this matrix.

     5. Compute the denominator d of E.
        Attempt to prove that result is correct by checking that

              H(d*E)*ncols(A)*H(A) < (prod of reduction primes)

        where H denotes the height.   If this fails, do step 4 with
        a few more primes.

    EXAMPLES::

        sage: A = matrix(QQ, 3, 7, [1..21])
        sage: from sage.matrix.misc import matrix_rational_echelon_form_multimodular
        sage: E, pivots = matrix_rational_echelon_form_multimodular(A)
        sage: E
        [ 1  0 -1 -2 -3 -4 -5]
        [ 0  1  2  3  4  5  6]
        [ 0  0  0  0  0  0  0]
        sage: pivots
        (0, 1)

        sage: A = matrix(QQ, 3, 4, [0,0] + [1..9] + [-1/2^20])
        sage: E, pivots = matrix_rational_echelon_form_multimodular(A)
        sage: E
        [                1                 0                 0 -10485761/1048576]
        [                0                 1                 0  27262979/4194304]
        [                0                 0                 1                 2]
        sage: pivots
        (0, 1, 2)

        sage: A.echelon_form()
        [                1                 0                 0 -10485761/1048576]
        [                0                 1                 0  27262979/4194304]
        [                0                 0                 1                 2]
        sage: A.pivots()
        (0, 1, 2)
    """
    if proof is None:
        from sage.structure.proof.proof import get_flag
        proof = get_flag(proof, "linear_algebra")

    verbose("Multimodular echelon algorithm on %s x %s matrix"%(self._nrows, self._ncols), caller_name="multimod echelon")
    cdef Matrix E
    if self._nrows == 0 or self._ncols == 0:
        return self, ()

    B, _ = self._clear_denom()

    height = self.height()
    if height_guess is None:
        height_guess = 10000000*(height+100)
    tm = verbose("height_guess = %s"%height_guess, level=2, caller_name="multimod echelon")

    if proof:
        M = self._ncols * height_guess * height  +  1
    else:
        M = height_guess + 1

    if self.is_sparse():
        from sage.matrix.matrix_modn_sparse import MAX_MODULUS
        p = MAX_MODULUS + 1
    else:
        from sage.matrix.matrix_modn_dense_double import MAX_MODULUS
        p = MAX_MODULUS + 1
    t = None
    X = []
    best_pivots = []
    prod = 1
    problem = 0
    lifts = {}
    while True:
        p = previous_prime(p)
        while prod < M:
            problem = problem + 1
            if problem > 50:
                verbose("echelon multi-modular possibly not converging?", caller_name="multimod echelon")
            t = verbose("echelon modulo p=%s (%.2f%% done)"%(
                       p, 100*float(len(str(prod))) / len(str(M))), level=2, caller_name="multimod echelon")

            # We use denoms=False, since we made self integral by calling clear_denom above.
            A = B._mod_int(p)
            t = verbose("time to reduce matrix mod p:",t, level=2, caller_name="multimod echelon")
            A.echelonize()
            t = verbose("time to put reduced matrix in echelon form:",t, level=2, caller_name="multimod echelon")

            # a worthwhile check / shortcut.
            if self._nrows >= self._ncols and self._nrows == len(A.pivots()):
                verbose("done: the echelon form mod p is the identity matrix and possibly some 0 rows", caller_name="multimod echelon")
                E = self.parent()(0)
                one = self.base_ring().one()
                for i in range(self._nrows):
                    E.set_unsafe(i, i, one)
                return E, tuple(range(self._nrows))

            c = cmp_pivots(best_pivots, A.pivots())
            if c <= 0:
                best_pivots = A.pivots()
                X.append(A)
                prod = prod * p
            else:
                # do not save A since it is bad.
                verbose("Excluding this prime (bad pivots).", caller_name="multimod echelon")
            t = verbose("time for pivot compare", t, level=2, caller_name="multimod echelon")
            p = previous_prime(p)
        # Find set of best matrices.
        Y = []
        # recompute product, since may drop bad matrices
        prod = 1
        t = verbose("now comparing pivots and dropping any bad ones", level=2, t=t, caller_name="multimod echelon")
        for i in range(len(X)):
            if cmp_pivots(best_pivots, X[i].pivots()) <= 0:
                p = X[i].base_ring().order()
                if p not in lifts:
                    t0 = verbose("Lifting a good matrix", level=2, caller_name = "multimod echelon")
                    lift = X[i].lift()
                    lifts[p] = (lift, p)
                    verbose("Finished lift", level=2, caller_name= "multimod echelon", t=t0)
                Y.append(lifts[p])
                prod = prod * X[i].base_ring().order()
        verbose("finished comparing pivots", level=2, t=t, caller_name="multimod echelon")
        try:
            if not Y:
                raise ValueError("not enough primes")
            t = verbose("start crt linear combination", level=2, caller_name="multimod echelon")
            a = CRT_basis([w[1] for w in Y])
            t = verbose('got crt basis', level=2, t=t, caller_name="multimod echelon")

            # take the linear combination of the lifts of the elements
            # of Y times coefficients in a
            L = a[0]*(Y[0][0])
            assert Y[0][0].is_sparse() == L.is_sparse()
            for j in range(1,len(Y)):
                L += a[j]*(Y[j][0])
            verbose("time to take linear combination of matrices over ZZ is",t, level=2, caller_name="multimod echelon")
            t = verbose("now doing rational reconstruction", level=2, caller_name="multimod echelon")
            E = L.rational_reconstruction(prod)
            L = 0  # free memory
            verbose('rational reconstruction completed', t, level=2, caller_name="multimod echelon")
        except ValueError as msg:
            verbose(msg, level=2)
            verbose("Not enough primes to do CRT lift; redoing with several more primes.", level=2, caller_name="multimod echelon")
            M = prod * p*p*p
            continue

        if not proof:
            verbose("Not checking validity of result (since proof=False).", level=2, caller_name="multimod echelon")
            break
        d   = E.denominator()
        hdE = int((d*E).height())
        if hdE * self.ncols() * height < prod:
            verbose("Validity of result checked.", level=2, caller_name="multimod echelon")
            break
        verbose("Validity failed; trying again with more primes.", level=2, caller_name="multimod echelon")
        M = prod * p*p*p
    #end while
    verbose("total time",tm, level=2, caller_name="multimod echelon")
    return E, tuple(best_pivots)


def cmp_pivots(x, y):
    """
    Compare two sequences of pivot columns.

    If x is shorter than y, return -1, i.e., x < y, "not as good".
    If x is longer than y, then x > y, so "better" and return +1.
    If the length is the same, then x is better, i.e., x > y
    if the entries of x are correspondingly <= those of y with
    one being strictly less.

    INPUT:

    - x, y -- lists or tuples of integers

    EXAMPLES:

    We illustrate each of the above comparisons. ::

        sage: sage.matrix.misc.cmp_pivots([1,2,3], [4,5,6,7])
        -1
        sage: sage.matrix.misc.cmp_pivots([1,2,3,5], [4,5,6])
        1
        sage: sage.matrix.misc.cmp_pivots([1,2,4], [1,2,3])
        -1
        sage: sage.matrix.misc.cmp_pivots([1,2,3], [1,2,3])
        0
        sage: sage.matrix.misc.cmp_pivots([1,2,3], [1,2,4])
        1
    """
    x = tuple(x)
    y = tuple(y)
    if len(x) < len(y):
        return -1
    if len(x) > len(y):
        return 1
    if x < y:
        return 1
    elif x == y:
        return 0
    else:
        return -1


def hadamard_row_bound_mpfr(Matrix A):
    """
    Given a matrix A with entries that coerce to RR, compute the row
    Hadamard bound on the determinant.

    INPUT:

        A -- a matrix over RR

    OUTPUT:

        integer -- an integer n such that the absolute value of the
                   determinant of this matrix is at most $10^n$.

    EXAMPLES:

    We create a very large matrix, compute the row Hadamard bound,
    and also compute the row Hadamard bound of the transpose, which
    happens to be sharp. ::

        sage: a = matrix(ZZ, 2, [2^10000,3^10000,2^50,3^19292])
        sage: import sage.matrix.misc
        sage: sage.matrix.misc.hadamard_row_bound_mpfr(a.change_ring(RR))
        13976
        sage: len(str(a.det()))
        12215
        sage: sage.matrix.misc.hadamard_row_bound_mpfr(a.transpose().change_ring(RR))
        12215

    Note that in the above example using RDF would overflow::

        sage: b = a.change_ring(RDF)
        sage: b._hadamard_row_bound()
        Traceback (most recent call last):
        ...
        OverflowError: cannot convert float infinity to integer
    """
    if not is_RealField(A.base_ring()):
        raise TypeError("A must have base field an mpfr real field.")

    cdef RealNumber a, b
    cdef mpfr_t s, d, pr
    cdef Py_ssize_t i, j

    mpfr_init(s)
    mpfr_init(d)
    mpfr_init(pr)
    mpfr_set_si(d, 0, MPFR_RNDU)

    for i in range(A._nrows):
        mpfr_set_si(s, 0, MPFR_RNDU)
        for j in range(A._ncols):
            sig_check()
            a = A.get_unsafe(i, j)
            mpfr_mul(pr, a.value, a.value, MPFR_RNDU)
            mpfr_add(s, s, pr, MPFR_RNDU)
        mpfr_log10(s, s, MPFR_RNDU)
        mpfr_add(d, d, s, MPFR_RNDU)
    b = a._new()
    mpfr_set(b.value, d, MPFR_RNDU)
    b /= 2
    mpfr_clear(s)
    mpfr_clear(d)
    mpfr_clear(pr)
    return b.ceil()
