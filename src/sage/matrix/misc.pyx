"""
Misc matrix algorithms
"""

from cysignals.signals cimport sig_check

from sage.arith.misc import CRT_basis, next_prime, previous_prime
from sage.arith.rational_reconstruction cimport mpq_rational_reconstruction
from sage.data_structures.binary_search cimport *
from sage.ext.mod_int cimport *
from sage.libs.gmp.mpq cimport *
from sage.libs.gmp.mpz cimport *
from sage.misc.lazy_import import LazyImport
from sage.misc.lazy_string import lazy_string
from sage.misc.verbose import verbose
from sage.modules.vector_integer_sparse cimport *
from sage.modules.vector_modn_sparse cimport *
from sage.modules.vector_rational_sparse cimport *
from sage.rings.integer cimport Integer
from sage.rings.rational_field import QQ

from sage.matrix.matrix0 cimport Matrix
from sage.matrix.matrix_integer_sparse cimport Matrix_integer_sparse
from sage.matrix.matrix_rational_sparse cimport Matrix_rational_sparse

matrix_integer_dense_rational_reconstruction = \
  LazyImport('sage.matrix.misc_flint', 'matrix_integer_dense_rational_reconstruction',
             deprecation=35758)
hadamard_row_bound_mpfr = \
  LazyImport('sage.matrix.misc_mpfr', 'hadamard_row_bound_mpfr',
             deprecation=35758)


def matrix_integer_sparse_rational_reconstruction(Matrix_integer_sparse A, Integer N):
    r"""
    Given a sparse matrix over the integers and an integer modulus, do
    rational reconstruction on all entries of the matrix, viewed as
    numbers mod `N`.

    EXAMPLES::

        sage: A = matrix(ZZ, 3, 4, [(1/3)%500, 2, 3, (-4)%500, 7, 2, 2, 3, 4, 3, 4, (5/7)%500], sparse=True)
        sage: from sage.matrix.misc import matrix_integer_sparse_rational_reconstruction
        sage: matrix_integer_sparse_rational_reconstruction(A, 500)
        [1/3   2   3  -4]
        [  7   2   2   3]
        [  4   3   4 5/7]

    TESTS:

    Check that :issue:`9345` is fixed::

        sage: A = random_matrix(ZZ, 3, sparse=True)
        sage: sage.matrix.misc.matrix_integer_sparse_rational_reconstruction(A, 0)
        Traceback (most recent call last):
        ...
        ZeroDivisionError: The modulus cannot be zero

    Check that :issue:`42533` is fixed: the GMP temporaries must be released
    even when reconstruction of an entry fails, which is the normal signal
    that not enough primes have been used yet::

        sage: import resource
        sage: A = matrix(ZZ, 1, 1, {(0, 0): 12345}, sparse=True)

    The entry is genuinely not reconstructible modulo 13 (it reduces to 8,
    while the bound on numerator and denominator is 2), so every call below
    really does take the failing path::

        sage: matrix_integer_sparse_rational_reconstruction(A, 13)
        Traceback (most recent call last):
        ...
        ValueError: rational reconstruction does not exist

    ::

        sage: def leak(N):
        ....:     before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        ....:     for _ in range(N):
        ....:         try:
        ....:             matrix_integer_sparse_rational_reconstruction(A, 13)
        ....:         except ValueError:
        ....:             pass
        ....:     after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        ....:     return (after - before) * 1024   # ru_maxrss is in kilobytes

    Loop (at most 30 times) until we have 6 consecutive zeros when calling
    ``leak(10000)``. Before the fix each failing call leaked about 113 bytes,
    so 10000 iterations grew the resident set by more than a megabyte::

        sage: zeros = 0
        sage: for i in range(30):  # long time
        ....:     n = leak(10000)
        ....:     print("Leaked {} bytes".format(n))
        ....:     if n == 0:
        ....:         zeros += 1
        ....:         if zeros >= 6:
        ....:             break
        ....:     else:
        ....:         zeros = 0
        Leaked...
        Leaked 0 bytes
        Leaked 0 bytes
        Leaked 0 bytes
        Leaked 0 bytes
        Leaked 0 bytes

    Since ``ru_maxrss`` is a high-water mark, the printed lines alone would
    also match a run that exhausted the 30 iterations without ever reaching
    six consecutive zeros, so check that explicitly::

        sage: zeros >= 6  # long time
        True
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
    mpz_init(bnd)

    try:
        _bnd = (N//2).isqrt()
        mpz_set(bnd, _bnd.value)
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
    finally:
        mpq_clear(t)
        mpz_clear(denom)
        mpz_clear(a)
        mpz_clear(other_bnd)
        mpz_clear(bnd)

    return R


def _multimodular_echelon_progress(p, prod, M):
    """
    Return a progress message for modular echelon form.

    TESTS::

        sage: from sage.matrix.misc import _multimodular_echelon_progress
        sage: _multimodular_echelon_progress(7, 10, 100)
        'echelon modulo p=7 (66.67% done)'
    """
    return "echelon modulo p=%s (%.2f%% done)" % (
        p, 100*float(len(str(prod))) / len(str(M)))


_sparse_prime_product_cache = {}


def _sparse_prime_product(max_modulus):
    r"""
    Return the product of all primes below ``max_modulus``, cached.

    This is the largest modulus the optimized sparse modular backend can
    reach, so it is what decides whether that backend can finish a given
    reconstruction at all.  It depends on nothing but ``max_modulus``, hence
    the cache: recomputing it costs about a millisecond, which would dominate
    the echelon form of a small matrix.

    TESTS::

        sage: from sage.matrix.misc import _sparse_prime_product
        sage: _sparse_prime_product(10)
        210
        sage: _sparse_prime_product(10) is _sparse_prime_product(10)
        True
        sage: from sage.matrix.matrix_modn_sparse import MAX_MODULUS
        sage: _sparse_prime_product(MAX_MODULUS).nbits()
        66444
    """
    cached = _sparse_prime_product_cache.get(max_modulus)
    if cached is not None:
        return cached
    from sage.arith.misc import prime_range
    from sage.misc.misc_c import prod as product
    cached = Integer(product(prime_range(max_modulus)))
    _sparse_prime_product_cache[max_modulus] = cached
    return cached


def _multimodular_echelon_next_modulus(M, prod):
    r"""
    Return the target modulus to use after a failed reconstruction attempt.

    ``M`` grows geometrically, but the result always exceeds ``prod``: a
    retry must buy at least one new prime, otherwise the next pass would
    reconstruct from the very same primes and fail in the same way.

    TESTS::

        sage: from sage.matrix.misc import _multimodular_echelon_next_modulus
        sage: _multimodular_echelon_next_modulus(2^60, 2^61) == 2^73
        True

    Growth alone can leave ``M`` behind ``prod``; then take the smallest
    modulus that requires another prime::

        sage: _multimodular_echelon_next_modulus(2, 10^6)
        1000001
    """
    M = M << (M.bit_length() // 5 + 1)
    if M <= prod:
        M = prod + 1
    return M


def _multimodular_echelon_candidate_has_shape(Matrix E, pivots):
    r"""
    Return whether ``E`` has reduced-echelon shape for ``pivots``.

    This checks the structural invariant used by the exact row-space
    certificate in :func:`matrix_rational_echelon_form_multimodular`.

    TESTS::

        sage: from sage.matrix.misc import _multimodular_echelon_candidate_has_shape
        sage: _multimodular_echelon_candidate_has_shape(
        ....:     matrix(QQ, [[0, 1]]), (1,))
        True
        sage: _multimodular_echelon_candidate_has_shape(
        ....:     matrix(QQ, [[46337, 1]], sparse=True), (1,))
        False
        sage: _multimodular_echelon_candidate_has_shape(
        ....:     matrix(QQ, [[1, 2], [0, 1]]), (0, 1))
        False
        sage: _multimodular_echelon_candidate_has_shape(
        ....:     matrix(QQ, [[0, 1], [1, 0]]), (1, 0))
        False
    """
    cdef Py_ssize_t i, j, pivot
    cdef Py_ssize_t r = len(pivots)
    if (r > E.nrows()
            or any(pivot < 0 or pivot >= E.ncols() for pivot in pivots)
            or any(pivots[i] >= pivots[i + 1] for i in range(r - 1))):
        return False
    if (not E[r:].is_zero()
            or not E[:r].matrix_from_columns(pivots).is_one()):
        return False
    if E.is_sparse():
        for i, j in E.nonzero_positions(copy=False):
            if i < r and j < pivots[i]:
                return False
    else:
        for i in range(r):
            pivot = pivots[i]
            for j in range(pivot):
                if not E.get_is_zero_unsafe(i, j):
                    return False
    return True


def matrix_rational_echelon_form_multimodular(Matrix self, height_guess=None, proof=None):
    """
    Return reduced row-echelon form using a multi-modular
    algorithm.  Does not change ``self``.

    REFERENCE: Chapter 7 of Stein's "Explicitly Computing Modular Forms".

    INPUT:

    - ``height_guess`` -- integer or ``None``; an estimate for `H(dE)`, where
      `E` is the output echelon form and `d` is its denominator
    - ``proof`` -- boolean or ``None`` (default: ``None``, see
      ``proof.linear_algebra`` or ``sage.structure.proof``). Note that the
      global Sage default is proof=True.  With ``proof=False``, the algorithm
      uses an explicitly supplied ``height_guess`` to attempt reconstruction
      early, subject to backend fallback, and checks each candidate first with
      the usual height certificate; if that fails, it checks the
      reduced-echelon shape and an exact row-space identity.  The early
      attempt can reduce the number of modular images when the output is much
      smaller than the input, but failed
      reconstructions and the exact check can also make it slower

    OUTPUT: a pair consisting of a matrix in echelon form and a tuple of pivot
    positions.

    .. NOTE::

        A sparse input that is already close to dense, and wide enough that
        its echelon form can be far taller than itself, is delegated to dense
        FLINT: the optimized sparse modular backend is capped at small moduli
        and cannot always reach the required prime product.  The result is
        still returned sparse, but ``height_guess`` and ``proof`` have no
        effect on that path, since FLINT is exact.

    ALGORITHM:

    The following is a modular algorithm for computing the echelon
    form.  Define the height of a matrix to be the max of the
    absolute values of the entries.

    Given Matrix A with n columns (self).

     0. Rescale each row of the input matrix A independently to a primitive
        integer row.  This does not change echelon form and makes reduction
        modulo lots of primes significantly easier if there were
        denominators.  Henceforth we assume A has integer entries.

     1. Let c be a guess for H(dE), where E is the output echelon form
        and d is its denominator.  E.g., c=1000 if the matrix is very
        sparse and the application is to computing modular symbols.

     2. Let M = n * c * H(A) + 1, where n is the number of columns of A.
        If ``proof=False`` and c was supplied explicitly, start instead with
        M = c + 1 in order to attempt rational reconstruction earlier.

     3. List primes p_1, p_2, ..., such that the product of
        the p_i is at least M.  Sparse matrices ordinarily start with the
        primes supported by the optimized sparse modular backend.  If M
        already exceeds the product of that complete finite range, start
        genuinely sparse inputs directly with larger primes and generic
        sparse matrices; otherwise switch to those matrices if the optimized
        range is exhausted later.  When a wide input is already at least
        three-quarters dense and a heuristic based on the target and a
        Hadamard bound indicates that the complete optimized prime range may
        be insufficient, use dense FLINT instead.

     4. Try to compute the rational reconstruction CRT echelon form
        of A mod the product of the p_i.  If rational
        reconstruction fails, compute more echelon forms modulo subsequent
        primes, and attempt again.  Make sure to keep the
        result of CRT on the primes from before, so we don't have
        to do that computation again.  Let E be this matrix.

     5. Compute the denominator d of E.  Attempt to prove that the result is
        correct by checking that

              H(d*E)*ncols(A)*H(A) < (prod of reduction primes)

        where H denotes the height.  With ``proof=False``, if this certificate
        fails, check the exact identity

              d*A = A[:, P] * (d*E)[:r]

        where P is the tuple of r pivot columns, after explicitly checking
        that E has reduced-echelon shape for P.  If the applicable check
        fails, do step 4 with more primes.

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

    A small benchmark, showing that flint fraction-free multimodular algorithm
    is always faster than the fraction-free multimodular algorithm implemented in Python::

        sage: import copy
        sage: def benchmark(num_row, num_col, entry_size, timeout=2, integer_coefficient=True):
        ....:     A = matrix(QQ, [[
        ....:         randint(1, 2^entry_size) if integer_coefficient else ZZ(randint(1, 2^entry_size))/randint(1, 2^entry_size)
        ....:         for col in range(num_col)] for row in range(num_row)])
        ....:     data=[]
        ....:     for algorithm in ("flint:fflu", "flint:multimodular", "padic", "multimodular"):
        ....:         # classical is too slow
        ....:         B = copy.copy(A)
        ....:         t = walltime()
        ....:         alarm(timeout)
        ....:         try:
        ....:             B.echelonize(algorithm=algorithm)
        ....:         except AlarmInterrupt:
        ....:             pass
        ....:         finally:
        ....:             cancel_alarm()
        ....:         data.append((round(walltime(t), 4), algorithm))
        ....:     return sorted(data)
        sage: benchmark(20, 20, 10000)  # long time
        [...'flint:multimodular'...'multimodular'...'flint:fflu'...]
        sage: benchmark(39, 40, 200)  # long time
        [...'flint:multimodular'...'flint:fflu'...'multimodular'...]

    In older versions of flint
    before this `issue <https://github.com/flintlib/flint/issues/2129>`_
    is fixed, ``algorithm='flint'`` (automatic choice) may be slower than
    ``algorithm='flint:multimodular'``.

    In this case, there are more columns than rows, which means the resulting
    matrix has height much higher than the input matrix. We check that the function
    does not take too long::

        sage: A = matrix(QQ, [[randint(1, 2^500) for col in range(40)] for row in range(20)])
        sage: t = walltime()
        sage: A.echelonize(algorithm="multimodular")  # long time
        sage: t = walltime(t)  # long time
        sage: (t < 10, t)  # long time
        (True, ...)

    TESTS:

    Check that the correctness bound accounts for the integer matrix obtained
    after clearing denominators (:issue:`42411`)::

        sage: entries = [[1000001000, -1/501000500, 3], [1, -1, 1]]
        sage: expected = matrix(QQ, entries).echelon_form(algorithm='flint')
        sage: matrix(QQ, entries).echelon_form(algorithm='multimodular') == expected
        True
        sage: matrix(QQ, entries, sparse=True).echelon_form() == expected
        True

    An early reconstruction is checked exactly even when the standard
    validity bound is disabled::

        sage: output_height = (expected.denominator() * expected).height()
        sage: matrix(QQ, entries, sparse=True).echelon_form(
        ....:     height_guess=output_height, proof=False) == expected
        True

    First calibrate the instrumentation on a case that needs the exact
    fallback, then check that the inexpensive height certificate bypasses
    it::

        sage: import sage.matrix.misc as matrix_misc
        sage: A = matrix(QQ, [[10^6, 1, 0], [10^6 - 1, 1, 0]])
        sage: expected = A.echelon_form(algorithm='flint')
        sage: original_shape_check = matrix_misc._multimodular_echelon_candidate_has_shape
        sage: shape_checks = []
        sage: def count_shape_check(*args):
        ....:     shape_checks.append(None)
        ....:     return original_shape_check(*args)
        sage: matrix_misc._multimodular_echelon_candidate_has_shape = count_shape_check
        sage: try:
        ....:     D = 43*53*61
        ....:     fallback = matrix(QQ, [[D, D/43, D/53, D/61, 0],
        ....:                            [0,    0,    0,    0, 1]], sparse=True)
        ....:     fallback_expected = fallback.dense_matrix().echelon_form(
        ....:         algorithm='flint')
        ....:     fallback_result = fallback.echelon_form(
        ....:         algorithm='multimodular', height_guess=1, proof=False)
        ....:     hook_is_active = bool(shape_checks)
        ....:     shape_checks.clear()
        ....:     result = A.__copy__().echelon_form(
        ....:         algorithm='multimodular', proof=False)
        ....: finally:
        ....:     matrix_misc._multimodular_echelon_candidate_has_shape = original_shape_check
        sage: (result == expected, fallback_result == fallback_expected)
        (True, True)
        sage: (hook_is_active, shape_checks)
        (True, [])

    A premature rational reconstruction is rejected by the exact check even
    when ``proof=False``::

        sage: A = matrix(QQ, 2, 3,
        ....:            [-17, 4/5, 25/13, -2/3, 19/5, -17/11], sparse=True)
        sage: expected = A.dense_matrix().echelon_form(algorithm='flint')
        sage: output_height = (expected.denominator() * expected).height()
        sage: A.echelon_form(algorithm='multimodular',
        ....:                height_guess=output_height, proof=False) == expected
        True

    The full identity is necessary.  Rational reconstruction may reuse a
    shared denominator, so even a pivot residue of one need not reconstruct
    as one::

        sage: N = ZZ(499)
        sage: C = matrix(QQ, [[1, 1/7, 1/11, 1/13, 0],
        ....:                 [0,   0,    0,    0, 1]])
        sage: L = matrix(ZZ, 2, 5, [x % N for x in C.list()], sparse=True)
        sage: candidate = L.rational_reconstruction(N)
        sage: candidate[1, 4]
        3/1001
        sage: B0 = matrix(ZZ, [[1001, 143, 91, 77, 0],
        ....:                  [   0,   0,  0,  0, 1]], sparse=True)
        sage: d = candidate.denominator()
        sage: dE = (d*candidate).change_ring(ZZ)
        sage: P = (0, 4); nonpivots = (1, 2, 3)
        sage: (d * B0.matrix_from_columns(nonpivots)
        ....:  == B0.matrix_from_columns(P)
        ....:     * dE.matrix_from_columns(nonpivots))
        True
        sage: d*B0 == B0.matrix_from_columns(P) * dE
        False

    This situation can occur in the multimodular algorithm itself.  The first
    sparse reduction prime is 46337, and the following shared denominator is
    congruent to 8 modulo that prime::

        sage: D = 43*53*61
        sage: A = matrix(QQ, [[D, D/43, D/53, D/61, 0],
        ....:                 [0,    0,    0,    0, 1]], sparse=True)
        sage: expected = A.dense_matrix().echelon_form(algorithm='flint')
        sage: A.echelon_form(algorithm='multimodular',
        ....:                height_guess=1, proof=False) == expected
        True

    An exact row-space check also rejects a reconstruction obtained only from
    bad-pivot primes::

        sage: A = matrix(QQ, [[46337, 19453]], sparse=True)
        sage: expected = A.dense_matrix().echelon_form(algorithm='flint')
        sage: A.echelon_form(algorithm='multimodular',
        ....:                height_guess=46336, proof=False) == expected
        True

    Check both ends of the optimized sparse prime range.  A target reached
    only after using the final prime must be reconstructed immediately; if
    the target is still not reached, the computation continues with generic
    sparse matrices over larger primes::

        sage: import sage.matrix.misc as matrix_misc
        sage: import sage.matrix.matrix_modn_sparse as modn_sparse
        sage: original_max_modulus = modn_sparse.MAX_MODULUS
        sage: original_previous_prime = matrix_misc.previous_prime
        sage: original_next_prime = matrix_misc.next_prime
        sage: previous_prime_inputs = []
        sage: next_prime_inputs = []
        sage: def record_previous_prime(p):
        ....:     previous_prime_inputs.append(p)
        ....:     return original_previous_prime(p)
        sage: def record_next_prime(p, proof=None):
        ....:     next_prime_inputs.append(
        ....:         (p > modn_sparse.MAX_MODULUS, proof))
        ....:     return original_next_prime(p, proof=proof)
        sage: modn_sparse.MAX_MODULUS = 4
        sage: matrix_misc.previous_prime = record_previous_prime
        sage: matrix_misc.next_prime = record_next_prime
        sage: try:
        ....:     A = matrix(QQ, [[1, 10, 0]], sparse=True)
        ....:     expected = A.dense_matrix().echelon_form(algorithm='flint')
        ....:     continued = A.echelon_form(
        ....:         algorithm='multimodular', height_guess=1, proof=False)
        ....:     continuation_trace = (list(previous_prime_inputs),
        ....:                           list(next_prime_inputs))
        ....:     previous_prime_inputs.clear()
        ....:     next_prime_inputs.clear()
        ....:     direct_input = matrix(QQ, [[1, 10, 0]], sparse=True)
        ....:     direct = direct_input.echelon_form(
        ....:         algorithm='multimodular', height_guess=8, proof=False)
        ....:     direct_trace = (list(previous_prime_inputs),
        ....:                     list(next_prime_inputs))
        ....:     previous_prime_inputs.clear()
        ....:     next_prime_inputs.clear()
        ....:     dense_input = [[1, 1, 1], [1, 0, 1]]
        ....:     dense_fallback = matrix(QQ, dense_input,
        ....:                             sparse=True).echelon_form(
        ....:         algorithm='multimodular', height_guess=1, proof=True)
        ....:     dense_trace = (list(previous_prime_inputs),
        ....:                    list(next_prime_inputs))
        ....:     previous_prime_inputs.clear()
        ....:     next_prime_inputs.clear()
        ....:     early = matrix(QQ, dense_input, sparse=True).echelon_form(
        ....:         algorithm='multimodular', height_guess=1, proof=False)
        ....:     early_trace = (list(previous_prime_inputs),
        ....:                    list(next_prime_inputs))
        ....:     previous_prime_inputs.clear()
        ....:     next_prime_inputs.clear()
        ....:     boundary = matrix(QQ, [[1, 1]], sparse=True).echelon_form(
        ....:         algorithm='multimodular', height_guess=2, proof=True)
        ....:     boundary_trace = (list(previous_prime_inputs),
        ....:                       list(next_prime_inputs))
        ....: finally:
        ....:     matrix_misc.next_prime = original_next_prime
        ....:     matrix_misc.previous_prime = original_previous_prime
        ....:     modn_sparse.MAX_MODULUS = original_max_modulus
        sage: (continued == expected, continued.is_sparse(),
        ....:  continuation_trace)
        (True, True, ([5, 3], [(True, True)]))
        sage: (direct == expected, direct.is_sparse(), direct_trace)
        (True, True, ([], [(True, True)]))
        sage: dense_fallback
        [1 0 1]
        [0 1 0]
        sage: (dense_fallback.is_sparse(), dense_trace)
        (True, ([], []))
        sage: (early == dense_fallback, early_trace)
        (True, ([5], []))
        sage: (boundary, boundary_trace)
        ([1 1], ([5, 3], []))

    Check that independent row scaling keeps the performance-critical bound
    small when rows have distinct denominators.  Counting calls to
    ``previous_prime`` avoids a timing-dependent test::

        sage: import sage.matrix.misc as matrix_misc
        sage: denominators = [next_prime(10^4 + 100*i) for i in range(80)]
        sage: entries = {(i, i): 1/q for i, q in enumerate(denominators)}
        sage: entries.update({(i, 119): 1 for i in range(80)})
        sage: A = matrix(QQ, 80, 120, entries)
        sage: expected = matrix(QQ, 80, 120)
        sage: for i, q in enumerate(denominators):
        ....:     expected[i, i] = 1
        ....:     expected[i, 119] = q
        sage: original_previous_prime = matrix_misc.previous_prime
        sage: def run(proof):
        ....:     primes = []
        ....:     def count_previous_prime(p):
        ....:         primes.append(p)
        ....:         return original_previous_prime(p)
        ....:     matrix_misc.previous_prime = count_previous_prime
        ....:     try:
        ....:         result = A.__copy__().echelon_form(
        ....:             algorithm='multimodular', proof=proof)
        ....:     finally:
        ....:         matrix_misc.previous_prime = original_previous_prime
        ....:     return result, len(primes)
        sage: proved, proved_primes = run(True)
        sage: unproved, unproved_primes = run(False)
        sage: proved == unproved == expected
        True
        sage: all(0 < count < 10
        ....:     for count in (proved_primes, unproved_primes))
        True
    """
    if proof is None:
        from sage.structure.proof.proof import get_flag
        proof = get_flag(proof, "linear_algebra")

    verbose("Multimodular echelon algorithm on %s x %s matrix" % (self._nrows, self._ncols), caller_name="multimod echelon")
    cdef Matrix E, dE
    cdef Matrix_integer_sparse B_sparse
    cdef bint default_height_guess = height_guess is None
    cdef bint dense_flint_candidate = False
    cdef bint generic_sparse_moduli = False
    cdef bint sparse_input = self.is_sparse()
    cdef bint use_dense_flint = False
    cdef Py_ssize_t i, nonzero_count
    if self._nrows == 0 or self._ncols == 0:
        return self, ()

    B, height = self._clear_denom_rowwise()
    if not height:
        return self.parent()(0), ()

    if height_guess is None:
        # Base the heuristic on the primitive integer matrix actually used by
        # the modular computation, so denominator clearing is accounted for.
        height_guess = 10000000*(height+100)
    tm = verbose("height_guess = %s" % height_guess, level=2, caller_name="multimod echelon")

    cdef Integer M
    from sage.arith.misc import integer_floor as floor
    if proof or default_height_guess:
        # The modular computation is done with denominators cleared, so H(A)
        # in the reconstruction bound is the height of this integer matrix.
        M = floor(max(1, self._ncols * height_guess * height + 1))
    else:
        # An explicit height guess opts into early reconstruction.  Each
        # candidate is checked below.
        M = floor(max(2, height_guess + 1))

    if sparse_input:
        from sage.matrix.matrix_modn_sparse import MAX_MODULUS
        # The generic sparse fallback below is deliberately retained for
        # genuinely sparse matrices.  For a matrix whose storage is already
        # close to dense, however, FLINT is both faster and memory-appropriate
        # when the Hadamard reconstruction bound can exceed the complete
        # optimized sparse prime range.
        B_sparse = B
        nonzero_count = 0
        for i in range(self._nrows):
            nonzero_count += B_sparse._matrix[i].num_nonzero
        # Only a wide input can have an echelon form whose height is much
        # larger than its own, which is what exhausts the prime range; a
        # matrix with at least as many rows as columns either meets the
        # full-rank shortcut below or has its output bounded by ``ncols``.
        dense_flint_candidate = (
            self._nrows < self._ncols
            and 4 * nonzero_count >= 3 * self.nrows() * self.ncols())
        sparse_prime_product = _sparse_prime_product(MAX_MODULUS)
        if M > sparse_prime_product:
            if dense_flint_candidate:
                use_dense_flint = True
            else:
                from sage.rings.finite_rings.finite_field_constructor import GF
                generic_sparse_moduli = True
        elif dense_flint_candidate and (proof or default_height_guess):
            # Compare the Hadamard bound 2 * r^r * H(B)^(2r) against the prime
            # product.  Materializing it outright would build an integer of
            # about 2*r*H(B).nbits() bits, which for a wide input with large
            # entries dwarfs the matrix itself, so settle it from bit lengths
            # whenever they are conclusive.  H(B)^(2r) has between
            # 2r*(nbits(H)-1)+1 and 2r*nbits(H) bits; only when those straddle
            # the prime product do we evaluate exactly, and there the value is
            # by construction no larger than that product.
            rank_bound = Integer(self._nrows)
            head = 2 * rank_bound**rank_bound
            threshold_bits = sparse_prime_product.nbits()
            low_bits = head.nbits() + 2 * rank_bound * (height.nbits() - 1)
            high_bits = head.nbits() + 2 * rank_bound * height.nbits()
            if low_bits > threshold_bits:
                use_dense_flint = True
            elif high_bits < threshold_bits:
                use_dense_flint = False
            else:
                use_dense_flint = (head * height**(2 * rank_bound)
                                   >= sparse_prime_product)
        if use_dense_flint:
            verbose("Using dense FLINT because the optimized sparse "
                    "prime range may be insufficient.", level=2,
                    caller_name="multimod echelon")
            B_sparse = None
            B = None
            E = self.dense_matrix()
            E.echelonize(algorithm='flint:multimodular')
            pivots = E.pivots()
            result = E.sparse_matrix()
            # ``echelonize`` caches the dense matrix as its own echelon form,
            # so it would only be reclaimed by the cyclic collector.
            E.clear_cache()
            return result, pivots
        if generic_sparse_moduli:
            p = 1 << 255
        else:
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
        while prod < M:
            if generic_sparse_moduli:
                p = next_prime(p, proof=True)
            elif sparse_input and p <= 2:
                # Matrix_modn_sparse is limited to small C-int moduli.  Once
                # those primes are exhausted, keep the computation sparse but
                # use fewer, larger primes to amortize the generic backend.
                from sage.rings.finite_rings.finite_field_constructor import GF
                generic_sparse_moduli = True
                p = next_prime(1 << 255, proof=True)
            else:
                try:
                    p = previous_prime(p)
                except ValueError:
                    raise RuntimeError("ran out of primes in multimodular "
                                       "echelon form")
            problem = problem + 1
            if problem > 50:
                verbose("echelon multi-modular possibly not converging?", caller_name="multimod echelon")
            t = verbose(lazy_string(_multimodular_echelon_progress,
                                    p, prod, M),
                        level=2, caller_name="multimod echelon")

            # We use denoms=False, since the rows of B are integral.
            if generic_sparse_moduli:
                A = B.change_ring(GF(p))
            else:
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
        # Find set of best matrices.
        Y = []
        # recompute product, since may drop bad matrices
        prod = 1
        t = verbose("now comparing pivots and dropping any bad ones", level=2, t=t, caller_name="multimod echelon")
        for i in range(len(X)):
            if cmp_pivots(best_pivots, X[i].pivots()) <= 0:
                q = X[i].base_ring().order()
                if q not in lifts:
                    t0 = verbose("Lifting a good matrix", level=2, caller_name="multimod echelon")
                    lift = X[i].lift()
                    lifts[q] = (lift, q)
                    verbose("Finished lift", level=2, caller_name="multimod echelon", t=t0)
                Y.append(lifts[q])
                prod = prod * X[i].base_ring().order()
        verbose("finished comparing pivots", level=2, t=t, caller_name="multimod echelon")
        if prod < M:
            continue
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
            verbose("Not enough primes to do CRT lift; redoing with more "
                    "primes.", level=2,
                    caller_name="multimod echelon")
            M = _multimodular_echelon_next_modulus(M, prod)
            continue

        d = E.denominator()
        dE = d*E
        hdE = int(dE.height())
        if hdE * self.ncols() * height < prod:
            verbose("Validity checked using the height bound.", level=2,
                    caller_name="multimod echelon")
            break

        if not proof:
            r = len(best_pivots)
            if _multimodular_echelon_candidate_has_shape(E, best_pivots):
                dE_integer = dE[:r].change_ring(B.base_ring())
                if d*B == B.matrix_from_columns(best_pivots) * dE_integer:
                    verbose("Validity of early reconstruction checked.",
                            level=2, caller_name="multimod echelon")
                    break
            verbose("Early reconstruction failed the row-space check; "
                    "trying more primes.", level=2,
                    caller_name="multimod echelon")
            M = _multimodular_echelon_next_modulus(M, prod)
            continue

        verbose("Validity failed; trying again with more primes.", level=2, caller_name="multimod echelon")
        M = _multimodular_echelon_next_modulus(M, prod)
    #end while
    verbose("total time",tm, level=2, caller_name="multimod echelon")
    return E, tuple(best_pivots)


def cmp_pivots(x, y):
    r"""
    Compare two sequences of pivot columns.

    If `x` is shorter than `y`, return `-1`, i.e., `x < y`, "not as good".
    If `x` is longer than `y`, then `x > y`, so "better" and return `+1`.
    If the length is the same, then `x` is better, i.e., `x > y`
    if the entries of `x` are correspondingly `\leq` those of `y` with
    one being strictly less.

    INPUT:

    - ``x``, ``y`` -- lists or tuples of integers

    EXAMPLES:

    We illustrate each of the above comparisons. ::

        sage: from sage.matrix.misc import cmp_pivots
        sage: cmp_pivots([1,2,3], [4,5,6,7])
        -1
        sage: cmp_pivots([1,2,3,5], [4,5,6])
        1
        sage: cmp_pivots([1,2,4], [1,2,3])
        -1
        sage: cmp_pivots([1,2,3], [1,2,3])
        0
        sage: cmp_pivots([1,2,3], [1,2,4])
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
    if x == y:
        return 0
    return -1
