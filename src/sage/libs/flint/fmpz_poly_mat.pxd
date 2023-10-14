# distutils: libraries = flint
# distutils: depends = flint/fmpz_poly_mat.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void fmpz_poly_mat_init(fmpz_poly_mat_t mat, slong rows, slong cols)
    # Initialises a matrix with the given number of rows and columns for use.

    void fmpz_poly_mat_init_set(fmpz_poly_mat_t mat, const fmpz_poly_mat_t src)
    # Initialises a matrix ``mat`` of the same dimensions as ``src``,
    # and sets it to a copy of ``src``.

    void fmpz_poly_mat_clear(fmpz_poly_mat_t mat)
    # Frees all memory associated with the matrix. The matrix must be
    # reinitialised if it is to be used again.

    slong fmpz_poly_mat_nrows(const fmpz_poly_mat_t mat)
    # Returns the number of rows in ``mat``.

    slong fmpz_poly_mat_ncols(const fmpz_poly_mat_t mat)
    # Returns the number of columns in ``mat``.

    fmpz_poly_struct * fmpz_poly_mat_entry(const fmpz_poly_mat_t mat, slong i, slong j)
    # Gives a reference to the entry at row ``i`` and column ``j``.
    # The reference can be passed as an input or output variable to any
    # ``fmpz_poly`` function for direct manipulation of the matrix element.
    # No bounds checking is performed.

    void fmpz_poly_mat_set(fmpz_poly_mat_t mat1, const fmpz_poly_mat_t mat2)
    # Sets ``mat1`` to a copy of ``mat2``.

    void fmpz_poly_mat_swap(fmpz_poly_mat_t mat1, fmpz_poly_mat_t mat2)
    # Swaps ``mat1`` and ``mat2`` efficiently.

    void fmpz_poly_mat_swap_entrywise(fmpz_poly_mat_t mat1, fmpz_poly_mat_t mat2)
    # Swaps two matrices by swapping the individual entries rather than swapping
    # the contents of the structs.

    void fmpz_poly_mat_print(const fmpz_poly_mat_t mat, const char * x)
    # Prints the matrix ``mat`` to standard output, using the
    # variable ``x``.

    void fmpz_poly_mat_randtest(fmpz_poly_mat_t mat, flint_rand_t state, slong len, flint_bitcnt_t bits)
    # This is equivalent to applying ``fmpz_poly_randtest`` to all entries
    # in the matrix.

    void fmpz_poly_mat_randtest_unsigned(fmpz_poly_mat_t mat, flint_rand_t state, slong len, flint_bitcnt_t bits)
    # This is equivalent to applying ``fmpz_poly_randtest_unsigned`` to
    # all entries in the matrix.

    void fmpz_poly_mat_randtest_sparse(fmpz_poly_mat_t A, flint_rand_t state, slong len, flint_bitcnt_t bits, float density)
    # Creates a random matrix with the amount of nonzero entries given
    # approximately by the ``density`` variable, which should be a fraction
    # between 0 (most sparse) and 1 (most dense).
    # The nonzero entries will have random lengths between 1 and ``len``.

    void fmpz_poly_mat_zero(fmpz_poly_mat_t mat)
    # Sets ``mat`` to the zero matrix.

    void fmpz_poly_mat_one(fmpz_poly_mat_t mat)
    # Sets ``mat`` to the unit or identity matrix of given shape,
    # having the element 1 on the main diagonal and zeros elsewhere.
    # If ``mat`` is nonsquare, it is set to the truncation of a unit matrix.

    bint fmpz_poly_mat_equal(const fmpz_poly_mat_t mat1, const fmpz_poly_mat_t mat2)
    # Returns nonzero if ``mat1`` and ``mat2`` have the same shape and
    # all their entries agree, and returns zero otherwise.

    bint fmpz_poly_mat_is_zero(const fmpz_poly_mat_t mat)
    # Returns nonzero if all entries in ``mat`` are zero, and returns
    # zero otherwise.

    bint fmpz_poly_mat_is_one(const fmpz_poly_mat_t mat)
    # Returns nonzero if all entries of ``mat`` on the main diagonal
    # are the constant polynomial 1 and all remaining entries are zero,
    # and returns zero otherwise. The matrix need not be square.

    bint fmpz_poly_mat_is_empty(const fmpz_poly_mat_t mat)
    # Returns a non-zero value if the number of rows or the number of
    # columns in ``mat`` is zero, and otherwise returns
    # zero.

    bint fmpz_poly_mat_is_square(const fmpz_poly_mat_t mat)
    # Returns a non-zero value if the number of rows is equal to the
    # number of columns in ``mat``, and otherwise returns zero.

    slong fmpz_poly_mat_max_bits(const fmpz_poly_mat_t A)
    # Returns the maximum number of bits among the coefficients of the
    # entries in ``A``, or the negative of that value if any
    # coefficient is negative.

    slong fmpz_poly_mat_max_length(const fmpz_poly_mat_t A)
    # Returns the maximum polynomial length among all the entries in ``A``.

    void fmpz_poly_mat_transpose(fmpz_poly_mat_t B, const fmpz_poly_mat_t A)
    # Sets `B` to `A^t`.

    void fmpz_poly_mat_evaluate_fmpz(fmpz_mat_t B, const fmpz_poly_mat_t A, const fmpz_t x)
    # Sets the ``fmpz_mat_t`` ``B`` to ``A`` evaluated entrywise
    # at the point ``x``.

    void fmpz_poly_mat_scalar_mul_fmpz_poly(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, const fmpz_poly_t c)
    # Sets ``B`` to ``A`` multiplied entrywise by the polynomial ``c``.

    void fmpz_poly_mat_scalar_mul_fmpz(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, const fmpz_t c)
    # Sets ``B`` to ``A`` multiplied entrywise by the integer ``c``.

    void fmpz_poly_mat_add(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B)
    # Sets ``C`` to the sum of ``A`` and ``B``.
    # All matrices must have the same shape. Aliasing is allowed.

    void fmpz_poly_mat_sub(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B)
    # Sets ``C`` to the sum of ``A`` and ``B``.
    # All matrices must have the same shape. Aliasing is allowed.

    void fmpz_poly_mat_neg(fmpz_poly_mat_t B, const fmpz_poly_mat_t A)
    # Sets ``B`` to the negation of ``A``.
    # The matrices must have the same shape. Aliasing is allowed.

    void fmpz_poly_mat_mul(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B)
    # Sets ``C`` to the matrix product of ``A`` and ``B``.
    # The matrices must have compatible dimensions for matrix multiplication.
    # Aliasing is allowed. This function automatically chooses between
    # classical and KS multiplication.

    void fmpz_poly_mat_mul_classical(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B)
    # Sets ``C`` to the matrix product of ``A`` and ``B``,
    # computed using the classical algorithm. The matrices must have
    # compatible dimensions for matrix multiplication. Aliasing is allowed.

    void fmpz_poly_mat_mul_KS(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B)
    # Sets ``C`` to the matrix product of ``A`` and ``B``,
    # computed using Kronecker segmentation. The matrices must have
    # compatible dimensions for matrix multiplication. Aliasing is allowed.

    void fmpz_poly_mat_mullow(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B, slong len)
    # Sets ``C`` to the matrix product of ``A`` and ``B``,
    # truncating each entry in the result to length ``len``.
    # Uses classical matrix multiplication. The matrices must have
    # compatible dimensions for matrix multiplication. Aliasing is allowed.

    void fmpz_poly_mat_sqr(fmpz_poly_mat_t B, const fmpz_poly_mat_t A)
    # Sets ``B`` to the square of ``A``, which must be a square matrix.
    # Aliasing is allowed. This function automatically chooses between
    # classical and KS squaring.

    void fmpz_poly_mat_sqr_classical(fmpz_poly_mat_t B, const fmpz_poly_mat_t A)
    # Sets ``B`` to the square of ``A``, which must be a square matrix.
    # Aliasing is allowed. This function uses direct formulas for very small
    # matrices, and otherwise classical matrix multiplication.

    void fmpz_poly_mat_sqr_KS(fmpz_poly_mat_t B, const fmpz_poly_mat_t A)
    # Sets ``B`` to the square of ``A``, which must be a square matrix.
    # Aliasing is allowed. This function uses Kronecker segmentation.

    void fmpz_poly_mat_sqrlow(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, slong len)
    # Sets ``B`` to the square of ``A``, which must be a square matrix,
    # truncating all entries to length ``len``.
    # Aliasing is allowed. This function uses direct formulas for very small
    # matrices, and otherwise classical matrix multiplication.

    void fmpz_poly_mat_pow(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, ulong exp)
    # Sets ``B`` to ``A`` raised to the power ``exp``, where ``A``
    # is a square matrix. Uses exponentiation by squaring. Aliasing is allowed.

    void fmpz_poly_mat_pow_trunc(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, ulong exp, slong len)
    # Sets ``B`` to ``A`` raised to the power ``exp``, truncating
    # all entries to length ``len``, where ``A`` is a square matrix.
    # Uses exponentiation by squaring. Aliasing is allowed.

    void fmpz_poly_mat_prod(fmpz_poly_mat_t res, fmpz_poly_mat_t * const factors, slong n)
    # Sets ``res`` to the product of the ``n`` matrices given in
    # the vector ``factors``, all of which must be square and of the
    # same size. Uses binary splitting.

    slong fmpz_poly_mat_find_pivot_any(const fmpz_poly_mat_t mat, slong start_row, slong end_row, slong c)
    # Attempts to find a pivot entry for row reduction.
    # Returns a row index `r` between ``start_row`` (inclusive) and
    # ``stop_row`` (exclusive) such that column `c` in ``mat`` has
    # a nonzero entry on row `r`, or returns -1 if no such entry exists.
    # This implementation simply chooses the first nonzero entry
    # it encounters. This is likely to be a nearly optimal choice if all
    # entries in the matrix have roughly the same size, but can lead to
    # unnecessary coefficient growth if the entries vary in size.

    slong fmpz_poly_mat_find_pivot_partial(const fmpz_poly_mat_t mat, slong start_row, slong end_row, slong c)
    # Attempts to find a pivot entry for row reduction.
    # Returns a row index `r` between ``start_row`` (inclusive) and
    # ``stop_row`` (exclusive) such that column `c` in ``mat`` has
    # a nonzero entry on row `r`, or returns -1 if no such entry exists.
    # This implementation searches all the rows in the column and
    # chooses the nonzero entry of smallest degree. If there are several
    # entries with the same minimal degree, it chooses the entry with
    # the smallest coefficient bit bound. This heuristic typically reduces
    # coefficient growth when the matrix entries vary in size.

    slong fmpz_poly_mat_fflu(fmpz_poly_mat_t B, fmpz_poly_t den, slong * perm, const fmpz_poly_mat_t A, int rank_check)
    # Uses fraction-free Gaussian elimination to set (``B``, ``den``) to a
    # fraction-free LU decomposition of ``A`` and returns the
    # rank of ``A``. Aliasing of ``A`` and ``B`` is allowed.
    # Pivot elements are chosen with ``fmpz_poly_mat_find_pivot_partial``.
    # If ``perm`` is non-``NULL``, the permutation of
    # rows in the matrix will also be applied to ``perm``.
    # If ``rank_check`` is set, the function aborts and returns 0 if the
    # matrix is detected not to have full rank without completing the
    # elimination.
    # The denominator ``den`` is set to `\pm \operatorname{det}(A)`, where
    # the sign is decided by the parity of the permutation. Note that the
    # determinant is not generally the minimal denominator.

    slong fmpz_poly_mat_rref(fmpz_poly_mat_t B, fmpz_poly_t den, const fmpz_poly_mat_t A)
    # Sets (``B``, ``den``) to the reduced row echelon form of
    # ``A`` and returns the rank of ``A``. Aliasing of ``A`` and
    # ``B`` is allowed.
    # The denominator ``den`` is set to `\pm \operatorname{det}(A)`.
    # Note that the determinant is not generally the minimal denominator.

    void fmpz_poly_mat_trace(fmpz_poly_t trace, const fmpz_poly_mat_t mat)
    # Computes the trace of the matrix, i.e. the sum of the entries on
    # the main diagonal. The matrix is required to be square.

    void fmpz_poly_mat_det(fmpz_poly_t det, const fmpz_poly_mat_t A)
    # Sets ``det`` to the determinant of the square matrix ``A``. Uses
    # a direct formula, fraction-free LU decomposition, or interpolation,
    # depending on the size of the matrix.

    void fmpz_poly_mat_det_fflu(fmpz_poly_t det, const fmpz_poly_mat_t A)
    # Sets ``det`` to the determinant of the square matrix ``A``.
    # The determinant is computed by performing a fraction-free LU
    # decomposition on a copy of ``A``.

    void fmpz_poly_mat_det_interpolate(fmpz_poly_t det, const fmpz_poly_mat_t A)
    # Sets ``det`` to the determinant of the square matrix ``A``.
    # The determinant is computed by determining a bound `n` for its length,
    # evaluating the matrix at `n` distinct points, computing the determinant
    # of each integer matrix, and forming the interpolating polynomial.

    slong fmpz_poly_mat_rank(const fmpz_poly_mat_t A)
    # Returns the rank of ``A``. Performs fraction-free LU decomposition
    # on a copy of ``A``.

    int fmpz_poly_mat_inv(fmpz_poly_mat_t Ainv, fmpz_poly_t den, const fmpz_poly_mat_t A)
    # Sets (``Ainv``, ``den``) to the inverse matrix of ``A``.
    # Returns 1 if ``A`` is nonsingular and 0 if ``A`` is singular.
    # Aliasing of ``Ainv`` and ``A`` is allowed.
    # More precisely, ``det`` will be set to the determinant of ``A``
    # and ``Ainv`` will be set to the adjugate matrix of ``A``.
    # Note that the determinant is not necessarily the minimal denominator.
    # Uses fraction-free LU decomposition, followed by solving for
    # the identity matrix.

    slong fmpz_poly_mat_nullspace(fmpz_poly_mat_t res, const fmpz_poly_mat_t mat)
    # Computes the right rational nullspace of the matrix ``mat`` and
    # returns the nullity.
    # More precisely, assume that ``mat`` has rank `r` and nullity `n`.
    # Then this function sets the first `n` columns of ``res``
    # to linearly independent vectors spanning the nullspace of ``mat``.
    # As a result, we always have rank(``res``) `= n`, and
    # ``mat`` `\times` ``res`` is the zero matrix.
    # The computed basis vectors will not generally be in a reduced form.
    # In general, the polynomials in each column vector in the result
    # will have a nontrivial common GCD.

    int fmpz_poly_mat_solve(fmpz_poly_mat_t X, fmpz_poly_t den, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B)
    # Solves the equation `AX = B` for nonsingular `A`. More precisely, computes
    # (``X``, ``den``) such that `AX = B \times \operatorname{den}`.
    # Returns 1 if `A` is nonsingular and 0 if `A` is singular.
    # The computed denominator will not generally be minimal.
    # Uses fraction-free LU decomposition followed by fraction-free
    # forward and back substitution.

    int fmpz_poly_mat_solve_fflu(fmpz_poly_mat_t X, fmpz_poly_t den, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B)
    # Solves the equation `AX = B` for nonsingular `A`. More precisely, computes
    # (``X``, ``den``) such that `AX = B \times \operatorname{den}`.
    # Returns 1 if `A` is nonsingular and 0 if `A` is singular.
    # The computed denominator will not generally be minimal.
    # Uses fraction-free LU decomposition followed by fraction-free
    # forward and back substitution.

    void fmpz_poly_mat_solve_fflu_precomp(fmpz_poly_mat_t X, const slong * perm, const fmpz_poly_mat_t FFLU, const fmpz_poly_mat_t B)
    # Performs fraction-free forward and back substitution given a precomputed
    # fraction-free LU decomposition and corresponding permutation.
