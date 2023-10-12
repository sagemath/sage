# distutils: libraries = flint
# distutils: depends = flint/ca_mat.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    ca_ptr ca_mat_entry_ptr(ca_mat_t mat, long i, long j)
    # Returns a pointer to the entry at row *i* and column *j*.
    # Equivalent to :macro:`ca_mat_entry` but implemented as a function.

    void ca_mat_init(ca_mat_t mat, long r, long c, ca_ctx_t ctx)
    # Initializes the matrix, setting it to the zero matrix with *r* rows
    # and *c* columns.

    void ca_mat_clear(ca_mat_t mat, ca_ctx_t ctx)
    # Clears the matrix, deallocating all entries.

    void ca_mat_swap(ca_mat_t mat1, ca_mat_t mat2, ca_ctx_t ctx)
    # Efficiently swaps *mat1* and *mat2*.

    void ca_mat_window_init(ca_mat_t window, const ca_mat_t mat, long r1, long c1, long r2, long c2, ca_ctx_t ctx)
    # Initializes *window* to a window matrix into the submatrix of *mat*
    # starting at the corner at row *r1* and column *c1* (inclusive) and ending
    # at row *r2* and column *c2* (exclusive).

    void ca_mat_window_clear(ca_mat_t window, ca_ctx_t ctx)
    # Frees the window matrix.

    void ca_mat_set(ca_mat_t dest, const ca_mat_t src, ca_ctx_t ctx)
    void ca_mat_set_fmpz_mat(ca_mat_t dest, const fmpz_mat_t src, ca_ctx_t ctx)
    void ca_mat_set_fmpq_mat(ca_mat_t dest, const fmpq_mat_t src, ca_ctx_t ctx)
    # Sets *dest* to *src*. The operands must have identical dimensions.

    void ca_mat_set_ca(ca_mat_t mat, const ca_t c, ca_ctx_t ctx)
    # Sets *mat* to the matrix with the scalar *c* on the main diagonal
    # and zeros elsewhere.

    void ca_mat_transfer(ca_mat_t res, ca_ctx_t res_ctx, const ca_mat_t src, ca_ctx_t src_ctx)
    # Sets *res* to *src* where the corresponding context objects *res_ctx* and
    # *src_ctx* may be different.
    # This operation preserves the mathematical value represented by *src*,
    # but may result in a different internal representation depending on the
    # settings of the context objects.

    void ca_mat_randtest(ca_mat_t mat, flint_rand_t state, long depth, long bits, ca_ctx_t ctx)
    # Sets *mat* to a random matrix with entries having complexity up to
    # *depth* and *bits* (see :func:`ca_randtest`).

    void ca_mat_randtest_rational(ca_mat_t mat, flint_rand_t state, long bits, ca_ctx_t ctx)
    # Sets *mat* to a random rational matrix with entries up to *bits* bits in size.

    void ca_mat_randops(ca_mat_t mat, flint_rand_t state, long count, ca_ctx_t ctx)
    # Randomizes *mat* in-place by performing elementary row or column operations.
    # More precisely, at most count random additions or subtractions of distinct
    # rows and columns will be performed. This leaves the rank (and for square matrices,
    # the determinant) unchanged.

    void ca_mat_print(const ca_mat_t mat, ca_ctx_t ctx)
    # Prints *mat* to standard output. The entries are printed on separate lines.

    void ca_mat_printn(const ca_mat_t mat, long digits, ca_ctx_t ctx)
    # Prints a decimal representation of *mat* with precision specified by *digits*.
    # The entries are comma-separated with square brackets and comma separation
    # for the rows.

    void ca_mat_zero(ca_mat_t mat, ca_ctx_t ctx)
    # Sets all entries in *mat* to zero.

    void ca_mat_one(ca_mat_t mat, ca_ctx_t ctx)
    # Sets the entries on the main diagonal of *mat* to one, and
    # all other entries to zero.

    void ca_mat_ones(ca_mat_t mat, ca_ctx_t ctx)
    # Sets all entries in *mat* to one.

    void ca_mat_pascal(ca_mat_t mat, int triangular, ca_ctx_t ctx)
    # Sets *mat* to a Pascal matrix, whose entries are binomial coefficients.
    # If *triangular* is 0, constructs a full symmetric matrix
    # with the rows of Pascal's triangle as successive antidiagonals.
    # If *triangular* is 1, constructs the upper triangular matrix with
    # the rows of Pascal's triangle as columns, and if *triangular* is -1,
    # constructs the lower triangular matrix with the rows of Pascal's
    # triangle as rows.

    void ca_mat_stirling(ca_mat_t mat, int kind, ca_ctx_t ctx)
    # Sets *mat* to a Stirling matrix, whose entries are Stirling numbers.
    # If *kind* is 0, the entries are set to the unsigned Stirling numbers
    # of the first kind. If *kind* is 1, the entries are set to the signed
    # Stirling numbers of the first kind. If *kind* is 2, the entries are
    # set to the Stirling numbers of the second kind.

    void ca_mat_hilbert(ca_mat_t mat, ca_ctx_t ctx)
    # Sets *mat* to the Hilbert matrix, which has entries `A_{i,j} = 1/(i+j+1)`.

    void ca_mat_dft(ca_mat_t mat, int type, ca_ctx_t ctx)
    # Sets *mat* to the DFT (discrete Fourier transform) matrix of order *n*
    # where *n* is the smallest dimension of *mat* (if *mat* is not square,
    # the matrix is extended periodically along the larger dimension).
    # The *type* parameter selects between four different versions
    # of the DFT matrix (in which `\omega = e^{2\pi i/n}`):
    # * Type 0 -- entries `A_{j,k} = \omega^{-jk}`
    # * Type 1 -- entries `A_{j,k} = \omega^{jk} / n`
    # * Type 2 -- entries `A_{j,k} = \omega^{-jk} / \sqrt{n}`
    # * Type 3 -- entries `A_{j,k} = \omega^{jk} / \sqrt{n}`
    # The type 0 and 1 matrices are inverse pairs, and similarly for the
    # type 2 and 3 matrices.

    truth_t ca_mat_check_equal(const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
    # Compares *A* and *B* for equality.

    truth_t ca_mat_check_is_zero(const ca_mat_t A, ca_ctx_t ctx)
    # Tests if *A* is the zero matrix.

    truth_t ca_mat_check_is_one(const ca_mat_t A, ca_ctx_t ctx)
    # Tests if *A* has ones on the main diagonal and zeros elsewhere.

    void ca_mat_transpose(ca_mat_t res, const ca_mat_t A, ca_ctx_t ctx)
    # Sets *res* to the transpose of *A*.

    void ca_mat_conj(ca_mat_t res, const ca_mat_t A, ca_ctx_t ctx)
    # Sets *res* to the entrywise complex conjugate of *A*.

    void ca_mat_conj_transpose(ca_mat_t res, const ca_mat_t A, ca_ctx_t ctx)
    # Sets *res* to the conjugate transpose (Hermitian transpose) of *A*.

    void ca_mat_neg(ca_mat_t res, const ca_mat_t A, ca_ctx_t ctx)
    # Sets *res* to the negation of *A*.

    void ca_mat_add(ca_mat_t res, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
    # Sets *res* to the sum of *A* and *B*.

    void ca_mat_sub(ca_mat_t res, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
    # Sets *res* to the difference of *A* and *B*.

    void ca_mat_mul_classical(ca_mat_t res, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
    void ca_mat_mul_same_nf(ca_mat_t res, const ca_mat_t A, const ca_mat_t B, ca_field_t K, ca_ctx_t ctx)
    void ca_mat_mul(ca_mat_t res, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
    # Sets *res* to the matrix product of *A* and *B*.
    # The *classical* version uses classical multiplication.
    # The *same_nf* version assumes (not checked) that both *A* and *B*
    # have coefficients in the same simple algebraic number field *K*
    # or in `\mathbb{Q}`.
    # The default version chooses an algorithm automatically.

    void ca_mat_mul_si(ca_mat_t B, const ca_mat_t A, long c, ca_ctx_t ctx)
    void ca_mat_mul_fmpz(ca_mat_t B, const ca_mat_t A, const fmpz_t c, ca_ctx_t ctx)
    void ca_mat_mul_fmpq(ca_mat_t B, const ca_mat_t A, const fmpq_t c, ca_ctx_t ctx)
    void ca_mat_mul_ca(ca_mat_t B, const ca_mat_t A, const ca_t c, ca_ctx_t ctx)
    # Sets *B* to *A* multiplied by the scalar *c*.

    void ca_mat_div_si(ca_mat_t B, const ca_mat_t A, long c, ca_ctx_t ctx)
    void ca_mat_div_fmpz(ca_mat_t B, const ca_mat_t A, const fmpz_t c, ca_ctx_t ctx)
    void ca_mat_div_fmpq(ca_mat_t B, const ca_mat_t A, const fmpq_t c, ca_ctx_t ctx)
    void ca_mat_div_ca(ca_mat_t B, const ca_mat_t A, const ca_t c, ca_ctx_t ctx)
    # Sets *B* to *A* divided by the scalar *c*.

    void ca_mat_add_ca(ca_mat_t B, const ca_mat_t A, const ca_t c, ca_ctx_t ctx)
    void ca_mat_sub_ca(ca_mat_t B, const ca_mat_t A, const ca_t c, ca_ctx_t ctx)
    # Sets *B* to *A* plus or minus the scalar *c* (interpreted as a diagonal matrix).

    void ca_mat_addmul_ca(ca_mat_t B, const ca_mat_t A, const ca_t c, ca_ctx_t ctx)
    void ca_mat_submul_ca(ca_mat_t B, const ca_mat_t A, const ca_t c, ca_ctx_t ctx)
    # Sets the matrix *B* to *B* plus (or minus) the matrix *A* multiplied by the scalar *c*.

    void ca_mat_sqr(ca_mat_t B, const ca_mat_t A, ca_ctx_t ctx)
    # Sets *B* to the square of *A*.

    void ca_mat_pow_ui_binexp(ca_mat_t B, const ca_mat_t A, unsigned long exp, ca_ctx_t ctx)
    # Sets *B* to *A* raised to the power *exp*, evaluated using
    # binary exponentiation.

    void _ca_mat_ca_poly_evaluate(ca_mat_t res, ca_srcptr poly, long len, const ca_mat_t A, ca_ctx_t ctx)
    void ca_mat_ca_poly_evaluate(ca_mat_t res, const ca_poly_t poly, const ca_mat_t A, ca_ctx_t ctx)
    # Sets *res* to `f(A)` where *f* is the polynomial given by *poly*
    # and *A* is a square matrix. Uses the Paterson-Stockmeyer algorithm.

    truth_t ca_mat_find_pivot(long * pivot_row, ca_mat_t mat, long start_row, long end_row, long column, ca_ctx_t ctx)
    # Attempts to find a nonzero entry in *mat* with column index *column*
    # and row index between *start_row* (inclusive) and *end_row* (exclusive).
    # If the return value is ``T_TRUE``, such an element exists,
    # and *pivot_row* is set to the row index.
    # If the return value is ``T_FALSE``, no such element exists
    # (all entries in this part of the column are zero).
    # If the return value is ``T_UNKNOWN``, it is unknown whether such
    # an element exists (zero certification failed).
    # This function is destructive: any elements that are nontrivially
    # zero but can be certified zero will be overwritten by exact zeros.

    int ca_mat_lu_classical(long * rank, long * P, ca_mat_t LU, const ca_mat_t A, int rank_check, ca_ctx_t ctx)
    int ca_mat_lu_recursive(long * rank, long * P, ca_mat_t LU, const ca_mat_t A, int rank_check, ca_ctx_t ctx)
    int ca_mat_lu(long * rank, long * P, ca_mat_t LU, const ca_mat_t A, int rank_check, ca_ctx_t ctx)
    # Computes a generalized LU decomposition `A = PLU` of a given
    # matrix *A*, writing the rank of *A* to *rank*.
    # If *A* is a nonsingular square matrix, *LU* will be set to
    # a unit diagonal lower triangular matrix *L* and an upper
    # triangular matrix *U* (the diagonal of *L* will not be stored
    # explicitly).
    # If *A* is an arbitrary matrix of rank *r*, *U* will be in row
    # echelon form having *r* nonzero rows, and *L* will be lower
    # triangular but truncated to *r* columns, having implicit ones on
    # the *r* first entries of the main diagonal. All other entries will
    # be zero.
    # If a nonzero value for ``rank_check`` is passed, the function
    # will abandon the output matrix in an undefined state and set
    # the rank to 0 if *A* is detected to be rank-deficient.
    # The algorithm can fail if it fails to certify that a pivot
    # element is zero or nonzero, in which case the correct rank
    # cannot be determined.
    # The return value is 1 on success and 0 on failure. On failure,
    # the data in the output variables
    # ``rank``, ``P`` and ``LU`` will be meaningless.
    # The *classical* version uses iterative Gaussian elimination.
    # The *recursive* version uses a block recursive algorithm
    # to take advantage of fast matrix multiplication.

    int ca_mat_fflu(long * rank, long * P, ca_mat_t LU, ca_t den, const ca_mat_t A, int rank_check, ca_ctx_t ctx)
    # Similar to :func:`ca_mat_lu`, but computes a fraction-free
    # LU decomposition using the Bareiss algorithm.
    # The denominator is written to *den*.
    # Note that despite being "fraction-free", this algorithm may
    # introduce fractions due to incomplete symbolic simplifications.

    truth_t ca_mat_nonsingular_lu(long * P, ca_mat_t LU, const ca_mat_t A, ca_ctx_t ctx)
    # Wrapper for :func:`ca_mat_lu`.
    # If *A* can be proved to be invertible/nonsingular, returns ``T_TRUE`` and sets *P* and *LU* to a LU decomposition `A = PLU`.
    # If *A* can be proved to be singular, returns ``T_FALSE``.
    # If *A* cannot be proved to be either singular or nonsingular, returns ``T_UNKNOWN``.
    # When the return value is ``T_FALSE`` or ``T_UNKNOWN``, the
    # LU factorization is not completed and the values of
    # *P* and *LU* are arbitrary.

    truth_t ca_mat_nonsingular_fflu(long * P, ca_mat_t LU, ca_t den, const ca_mat_t A, ca_ctx_t ctx)
    # Wrapper for :func:`ca_mat_fflu`.
    # Similar to :func:`ca_mat_nonsingular_lu`, but computes a fraction-free
    # LU decomposition using the Bareiss algorithm.
    # The denominator is written to *den*.
    # Note that despite being "fraction-free", this algorithm may
    # introduce fractions due to incomplete symbolic simplifications.

    truth_t ca_mat_inv(ca_mat_t X, const ca_mat_t A, ca_ctx_t ctx)
    # Determines if the square matrix *A* is nonsingular, and if successful,
    # sets `X = A^{-1}` and returns ``T_TRUE``.
    # Returns ``T_FALSE`` if *A* is singular, and ``T_UNKNOWN`` if the
    # rank of *A* cannot be determined.

    truth_t ca_mat_nonsingular_solve_adjugate(ca_mat_t X, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
    truth_t ca_mat_nonsingular_solve_fflu(ca_mat_t X, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
    truth_t ca_mat_nonsingular_solve_lu(ca_mat_t X, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
    truth_t ca_mat_nonsingular_solve(ca_mat_t X, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
    # Determines if the square matrix *A* is nonsingular, and if successful,
    # solves `AX = B` and returns ``T_TRUE``.
    # Returns ``T_FALSE`` if *A* is singular, and ``T_UNKNOWN`` if the
    # rank of *A* cannot be determined.

    void ca_mat_solve_tril_classical(ca_mat_t X, const ca_mat_t L, const ca_mat_t B, int unit, ca_ctx_t ctx)
    void ca_mat_solve_tril_recursive(ca_mat_t X, const ca_mat_t L, const ca_mat_t B, int unit, ca_ctx_t ctx)
    void ca_mat_solve_tril(ca_mat_t X, const ca_mat_t L, const ca_mat_t B, int unit, ca_ctx_t ctx)
    void ca_mat_solve_triu_classical(ca_mat_t X, const ca_mat_t U, const ca_mat_t B, int unit, ca_ctx_t ctx)
    void ca_mat_solve_triu_recursive(ca_mat_t X, const ca_mat_t U, const ca_mat_t B, int unit, ca_ctx_t ctx)
    void ca_mat_solve_triu(ca_mat_t X, const ca_mat_t U, const ca_mat_t B, int unit, ca_ctx_t ctx)
    # Solves the lower triangular system `LX = B` or the upper triangular system
    # `UX = B`, respectively. It is assumed (not checked) that the diagonal
    # entries are nonzero. If *unit* is set, the main diagonal of *L* or *U*
    # is taken to consist of all ones, and in that case the actual entries on
    # the diagonal are not read at all and can contain other data.
    # The *classical* versions perform the computations iteratively while the
    # *recursive* versions perform the computations in a block recursive
    # way to benefit from fast matrix multiplication. The default versions
    # choose an algorithm automatically.

    void ca_mat_solve_fflu_precomp(ca_mat_t X, const long * perm, const ca_mat_t A, const ca_t den, const ca_mat_t B, ca_ctx_t ctx)
    void ca_mat_solve_lu_precomp(ca_mat_t X, const long * P, const ca_mat_t LU, const ca_mat_t B, ca_ctx_t ctx)
    # Solves `AX = B` given the precomputed nonsingular LU decomposition `A = PLU`
    # or fraction-free LU decomposition with denominator *den*.
    # The matrices `X` and `B` are allowed to be aliased with each other,
    # but `X` is not allowed to be aliased with `LU`.

    int ca_mat_rank(long * rank, const ca_mat_t A, ca_ctx_t ctx)
    # Computes the rank of the matrix *A*. If successful, returns 1 and
    # writes the rank to ``rank``. If unsuccessful, returns 0.

    int ca_mat_rref_fflu(long * rank, ca_mat_t R, const ca_mat_t A, ca_ctx_t ctx)
    int ca_mat_rref_lu(long * rank, ca_mat_t R, const ca_mat_t A, ca_ctx_t ctx)
    int ca_mat_rref(long * rank, ca_mat_t R, const ca_mat_t A, ca_ctx_t ctx)
    # Computes the reduced row echelon form (rref) of a given matrix.
    # On success, sets *R* to the rref of *A*, writes the rank to
    # *rank*, and returns 1. On failure to certify the correct rank,
    # returns 0, leaving the data in *rank* and *R* meaningless.
    # The *fflu* version computes a fraction-free LU decomposition and
    # then converts the output ro rref form. The *lu* version computes a
    # regular LU decomposition and then converts the output to rref form.
    # The default version uses an automatic algorithm choice and may
    # implement additional methods for special cases.

    int ca_mat_right_kernel(ca_mat_t X, const ca_mat_t A, ca_ctx_t ctx)
    # Sets *X* to a basis of the right kernel (nullspace) of *A*.
    # The output matrix *X* will be resized in-place to have a number
    # of columns equal to the nullity of *A*.
    # Returns 1 on success. On failure, returns 0 and leaves the data
    # in *X* meaningless.

    void ca_mat_trace(ca_t trace, const ca_mat_t mat, ca_ctx_t ctx)
    # Sets *trace* to the sum of the entries on the main diagonal of *mat*.

    void ca_mat_det_berkowitz(ca_t det, const ca_mat_t A, ca_ctx_t ctx)
    int ca_mat_det_lu(ca_t det, const ca_mat_t A, ca_ctx_t ctx)
    int ca_mat_det_bareiss(ca_t det, const ca_mat_t A, ca_ctx_t ctx)
    void ca_mat_det_cofactor(ca_t det, const ca_mat_t A, ca_ctx_t ctx)
    void ca_mat_det(ca_t det, const ca_mat_t A, ca_ctx_t ctx)
    # Sets *det* to the determinant of the square matrix *A*.
    # Various algorithms are available:
    # * The *berkowitz* version uses the division-free Berkowitz algorithm
    # performing `O(n^4)` operations. Since no zero tests are required, it
    # is guaranteed to succeed.
    # * The *cofactor* version performs cofactor expansion. This is currently
    # only supported for matrices up to size 4.
    # * The *lu* and *bareiss* versions use rational LU decomposition
    # and fraction-free LU decomposition (Bareiss algorithm) respectively,
    # requiring `O(n^3)` operations. These algorithms can fail if zero
    # certification fails (see :func:`ca_mat_nonsingular_lu`); they
    # return 1 for success and 0 for failure.
    # Note that the Bareiss algorithm, despite being "fraction-free",
    # may introduce fractions due to incomplete symbolic simplifications.
    # The default function chooses an algorithm automatically.
    # It will, in addition, recognize trivially rational and integer
    # matrices and evaluate those determinants using
    # :type:`fmpq_mat_t` or :type:`fmpz_mat_t`.
    # The various algorithms can produce different symbolic
    # forms of the same determinant. Which algorithm performs better
    # depends strongly and sometimes
    # unpredictably on the structure of the matrix.

    void ca_mat_adjugate_cofactor(ca_mat_t adj, ca_t det, const ca_mat_t A, ca_ctx_t ctx)
    void ca_mat_adjugate_charpoly(ca_mat_t adj, ca_t det, const ca_mat_t A, ca_ctx_t ctx)
    void ca_mat_adjugate(ca_mat_t adj, ca_t det, const ca_mat_t A, ca_ctx_t ctx)
    # Sets *adj* to the adjuate matrix of *A* and *det* to the determinant
    # of *A*, both computed simultaneously.
    # The *cofactor* version uses cofactor expansion.
    # The *charpoly* version computes and
    # evaluates the characteristic polynomial.
    # The default version uses an automatic algorithm choice.

    void _ca_mat_charpoly_berkowitz(ca_ptr cp, const ca_mat_t mat, ca_ctx_t ctx)
    void ca_mat_charpoly_berkowitz(ca_poly_t cp, const ca_mat_t mat, ca_ctx_t ctx)
    int _ca_mat_charpoly_danilevsky(ca_ptr cp, const ca_mat_t mat, ca_ctx_t ctx)
    int ca_mat_charpoly_danilevsky(ca_poly_t cp, const ca_mat_t mat, ca_ctx_t ctx)
    void _ca_mat_charpoly(ca_ptr cp, const ca_mat_t mat, ca_ctx_t ctx)
    void ca_mat_charpoly(ca_poly_t cp, const ca_mat_t mat, ca_ctx_t ctx)
    # Sets *poly* to the characteristic polynomial of *mat* which must be
    # a square matrix. If the matrix has *n* rows, the underscore method
    # requires space for `n + 1` output coefficients.
    # The *berkowitz* version uses a division-free algorithm
    # requiring `O(n^4)` operations.
    # The *danilevsky* version only performs `O(n^3)` operations, but
    # performs divisions and needs to check for zero which can fail.
    # This version returns 1 on success and 0 on failure.
    # The default version chooses an algorithm automatically.

    int ca_mat_companion(ca_mat_t mat, const ca_poly_t poly, ca_ctx_t ctx)
    # Sets *mat* to the companion matrix of *poly*.
    # This function verifies that the leading coefficient of *poly*
    # is provably nonzero and that the output matrix has the right size,
    # returning 1 on success.
    # It returns 0 if the leading coefficient of *poly* cannot be
    # proved nonzero or if the size of the output matrix does not match.

    int ca_mat_eigenvalues(ca_vec_t lambda, unsigned long * exp, const ca_mat_t mat, ca_ctx_t ctx)
    # Attempts to compute all complex eigenvalues of the given matrix *mat*.
    # On success, returns 1 and sets *lambda* to the distinct eigenvalues
    # with corresponding multiplicities in *exp*.
    # The eigenvalues are returned in arbitrary order.
    # On failure, returns 0 and leaves the values in *lambda* and *exp*
    # arbitrary.
    # This function effectively computes the characteristic polynomial
    # and then calls :type:`ca_poly_roots`.

    truth_t ca_mat_diagonalization(ca_mat_t D, ca_mat_t P, const ca_mat_t A, ca_ctx_t ctx)
    # Matrix diagonalization: attempts to compute a diagonal matrix *D*
    # and an invertible matrix *P* such that `A = PDP^{-1}`.
    # Returns ``T_TRUE`` if *A* is diagonalizable and the computation
    # succeeds, ``T_FALSE`` if *A* is provably not diagonalizable,
    # and ``T_UNKNOWN`` if it is unknown whether *A* is diagonalizable.
    # If the return value is not ``T_TRUE``, the values in *D* and *P*
    # are arbitrary.

    int ca_mat_jordan_blocks(ca_vec_t lambda, long * num_blocks, long * block_lambda, long * block_size, const ca_mat_t A, ca_ctx_t ctx)
    # Computes the blocks of the Jordan canonical form of *A*.
    # On success, returns 1 and sets *lambda* to the unique eigenvalues
    # of *A*, sets *num_blocks* to the number of Jordan blocks,
    # entry *i* of *block_lambda* to the index of the eigenvalue
    # in Jordan block *i*, and entry *i* of *block_size* to the size
    # of Jordan block *i*. On failure, returns 0, leaving arbitrary
    # values in the output variables.
    # The user should allocate space in *block_lambda* and *block_size*
    # for up to *n* entries where *n* is the size of the matrix.
    # The Jordan form is unique up to the ordering of blocks, which
    # is arbitrary.

    void ca_mat_set_jordan_blocks(ca_mat_t mat, const ca_vec_t lambda, long num_blocks, long * block_lambda, long * block_size, ca_ctx_t ctx)
    # Sets *mat* to the concatenation of the Jordan blocks
    # given in *lambda*, *num_blocks*, *block_lambda* and *block_size*.
    # See :func:`ca_mat_jordan_blocks` for an explanation of these
    # variables.

    int ca_mat_jordan_transformation(ca_mat_t mat, const ca_vec_t lambda, long num_blocks, long * block_lambda, long * block_size, const ca_mat_t A, ca_ctx_t ctx)
    # Given the precomputed Jordan block decomposition
    # (*lambda*, *num_blocks*, *block_lambda*, *block_size*) of the
    # square matrix *A*, computes the corresponding transformation
    # matrix *P* such that `A = P J P^{-1}`.
    # On success, writes *P* to *mat* and returns 1. On failure,
    # returns 0, leaving the value of *mat* arbitrary.

    int ca_mat_jordan_form(ca_mat_t J, ca_mat_t P, const ca_mat_t A, ca_ctx_t ctx)
    # Computes the Jordan decomposition `A = P J P^{-1}` of the given
    # square matrix *A*. The user can pass *NULL* for the output
    # variable *P*, in which case only *J* is computed.
    # On success, returns 1. On failure, returns 0, leaving the values
    # of *J* and *P* arbitrary.
    # This function is a convenience wrapper around
    # :func:`ca_mat_jordan_blocks`, :func:`ca_mat_set_jordan_blocks` and
    # :func:`ca_mat_jordan_transformation`. For computations with
    # the Jordan decomposition, it is often better to use those
    # methods directly since they give direct access to the
    # spectrum and block structure.

    int ca_mat_exp(ca_mat_t res, const ca_mat_t A, ca_ctx_t ctx)
    # Matrix exponential: given a square matrix *A*, sets *res* to
    # `e^A` and returns 1 on success. If unsuccessful, returns 0,
    # leaving the values in *res* arbitrary.
    # This function uses Jordan decomposition. The matrix exponential
    # always exists, but computation can fail if computing the Jordan
    # decomposition fails.

    truth_t ca_mat_log(ca_mat_t res, const ca_mat_t A, ca_ctx_t ctx)
    # Matrix logarithm: given a square matrix *A*, sets *res* to a
    # logarithm `\log(A)` and returns ``T_TRUE`` on success.
    # If *A* can be proved to have no logarithm, returns ``T_FALSE``.
    # If the existence of a logarithm cannot be proved, returns
    # ``T_UNKNOWN``.
    # This function uses the Jordan decomposition, and the branch of
    # the matrix logarithm is defined by taking the principal values
    # of the logarithms of all eigenvalues.
