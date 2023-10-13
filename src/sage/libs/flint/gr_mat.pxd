# distutils: libraries = flint
# distutils: depends = flint/gr_mat.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    gr_ptr gr_mat_entry_ptr(gr_mat_t mat, long i, long j, gr_ctx_t ctx)
    # Function returning a pointer to the entry at row *i* and column
    # *j* of the matrix *mat*. The indices must be in bounds.

    void gr_mat_init(gr_mat_t mat, long rows, long cols, gr_ctx_t ctx)
    # Initializes *mat* to a matrix with the given number of rows and
    # columns.

    int gr_mat_init_set(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
    # Initializes *res* to a copy of the matrix *mat*.

    void gr_mat_clear(gr_mat_t mat, gr_ctx_t ctx)
    # Clears the matrix.

    void gr_mat_swap(gr_mat_t mat1, gr_mat_t mat2, gr_ctx_t ctx)
    # Swaps *mat1* and *mat12* efficiently.

    int gr_mat_swap_entrywise(gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
    # Performs a deep swap of *mat1* and *mat2*, swapping the individual
    # entries rather than the top-level structures.

    void gr_mat_window_init(gr_mat_t window, const gr_mat_t mat, long r1, long c1, long r2, long c2, gr_ctx_t ctx)
    # Initializes *window* to a window matrix into the submatrix of *mat*
    # starting at the corner at row *r1* and column *c1* (inclusive) and ending
    # at row *r2* and column *c2* (exclusive).
    # The indices must be within bounds.

    void gr_mat_window_clear(gr_mat_t window, gr_ctx_t ctx)
    # Frees the window matrix.

    int gr_mat_write(gr_stream_t out, const gr_mat_t mat, gr_ctx_t ctx)
    # Write *mat* to the stream *out*.

    int gr_mat_print(const gr_mat_t mat, gr_ctx_t ctx)
    # Prints *mat* to standard output.

    truth_t gr_mat_equal(const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
    # Returns whether *mat1* and *mat2* are equal.

    truth_t gr_mat_is_zero(const gr_mat_t mat, gr_ctx_t ctx)
    truth_t gr_mat_is_one(const gr_mat_t mat, gr_ctx_t ctx)
    truth_t gr_mat_is_neg_one(const gr_mat_t mat, gr_ctx_t ctx)
    # Returns whether *mat* respectively is the zero matrix or
    # the scalar matrix with 1 or -1 on the main diagonal.

    truth_t gr_mat_is_scalar(const gr_mat_t mat, gr_ctx_t ctx)
    # Returns whether *mat* is a scalar matrix, being a diagonal matrix
    # with identical elements on the main diagonal.

    int gr_mat_zero(gr_mat_t res, gr_ctx_t ctx)
    # Sets *res* to the zero matrix.

    int gr_mat_one(gr_mat_t res, gr_ctx_t ctx)
    # Sets *res* to the scalar matrix with 1 on the main diagonal
    # and zero elsewhere.

    int gr_mat_set(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_set_fmpz_mat(gr_mat_t res, const fmpz_mat_t mat, gr_ctx_t ctx)
    int gr_mat_set_fmpq_mat(gr_mat_t res, const fmpq_mat_t mat, gr_ctx_t ctx)
    # Sets *res* to the value of *mat*.

    int gr_mat_set_scalar(gr_mat_t res, gr_srcptr c, gr_ctx_t ctx)
    int gr_mat_set_ui(gr_mat_t res, unsigned long c, gr_ctx_t ctx)
    int gr_mat_set_si(gr_mat_t res, long c, gr_ctx_t ctx)
    int gr_mat_set_fmpz(gr_mat_t res, const fmpz_t c, gr_ctx_t ctx)
    int gr_mat_set_fmpq(gr_mat_t res, const fmpq_t c, gr_ctx_t ctx)
    # Set *res* to the scalar matrix with *c* on the main diagonal
    # and zero elsewhere.

    int gr_mat_concat_horizontal(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)

    int gr_mat_concat_vertical(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)

    int gr_mat_transpose(gr_mat_t B, const gr_mat_t A, gr_ctx_t ctx)
    # Sets *B* to the transpose of *A*.

    int gr_mat_swap_rows(gr_mat_t mat, long * perm, long r, long s, gr_ctx_t ctx)
    # Swaps rows ``r`` and ``s`` of ``mat``.  If ``perm`` is non-``NULL``, the
    # permutation of the rows will also be applied to ``perm``.

    int gr_mat_swap_cols(gr_mat_t mat, long * perm, long r, long s, gr_ctx_t ctx)
    # Swaps columns ``r`` and ``s`` of ``mat``.  If ``perm`` is non-``NULL``, the
    # permutation of the columns will also be applied to ``perm``.

    int gr_mat_invert_rows(gr_mat_t mat, long * perm, gr_ctx_t ctx)
    # Swaps rows ``i`` and ``r - i`` of ``mat`` for ``0 <= i < r/2``, where
    # ``r`` is the number of rows of ``mat``. If ``perm`` is non-``NULL``, the
    # permutation of the rows will also be applied to ``perm``.

    int gr_mat_invert_cols(gr_mat_t mat, long * perm, gr_ctx_t ctx)
    # Swaps columns ``i`` and ``c - i`` of ``mat`` for ``0 <= i < c/2``, where
    # ``c`` is the number of columns of ``mat``. If ``perm`` is non-``NULL``, the
    # permutation of the columns will also be applied to ``perm``.

    truth_t gr_mat_is_empty(const gr_mat_t mat, gr_ctx_t ctx)
    # Returns whether *mat* is an empty matrix, having either zero
    # rows or zero column. This predicate is always decidable (even if
    # the underlying ring is not computable), returning
    # ``T_TRUE`` or ``T_FALSE``.

    truth_t gr_mat_is_square(const gr_mat_t mat, gr_ctx_t ctx)
    # Returns whether *mat* is a square matrix, having the same number
    # of rows as columns (not the same thing as being a perfect square!).
    # This predicate is always decidable (even if the underlying ring
    # is not computable), returning ``T_TRUE`` or ``T_FALSE``.

    int gr_mat_neg(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)

    int gr_mat_add(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)

    int gr_mat_sub(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)

    int gr_mat_mul_classical(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
    int gr_mat_mul_strassen(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
    int gr_mat_mul_generic(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
    int gr_mat_mul(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
    # Matrix multiplication. The default function can be overloaded by specific rings;
    # otherwise, it falls back to :func:`gr_mat_mul_generic` which currently
    # only performs classical multiplication.

    int gr_mat_sqr(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)

    int gr_mat_add_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)
    int gr_mat_sub_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)
    int gr_mat_mul_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)
    int gr_mat_addmul_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)
    int gr_mat_submul_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)
    int gr_mat_div_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)

    int _gr_mat_gr_poly_evaluate(gr_mat_t res, gr_srcptr poly, long len, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_gr_poly_evaluate(gr_mat_t res, const gr_poly_t poly, const gr_mat_t mat, gr_ctx_t ctx)
    # Sets *res* to the matrix obtained by evaluating the
    # scalar polynomial *poly* with matrix argument *mat*.

    truth_t gr_mat_is_upper_triangular(const gr_mat_t mat, gr_ctx_t ctx)
    truth_t gr_mat_is_lower_triangular(const gr_mat_t mat, gr_ctx_t ctx)
    # Returns whether *mat* is upper (respectively lower) triangular, having
    # zeros everywhere below (respectively above) the main diagonal.
    # The matrix need not be square.

    truth_t gr_mat_is_diagonal(const gr_mat_t mat, gr_ctx_t ctx)
    # Returns whether *mat* is a diagonal matrix, having zeros everywhere
    # except on the main diagonal.
    # The matrix need not be square.

    int gr_mat_mul_diag(gr_mat_t res, const gr_mat_t A, const gr_vec_t D, gr_ctx_t ctx)
    int gr_mat_diag_mul(gr_mat_t res, const gr_vec_t D, const gr_mat_t A, gr_ctx_t ctx)
    # Set *res* to the product `AD` or `DA` respectively, where `D` is
    # a diagonal matrix represented as a vector of entries.

    int gr_mat_find_nonzero_pivot_large_abs(long * pivot_row, gr_mat_t mat, long start_row, long end_row, long column, gr_ctx_t ctx)
    int gr_mat_find_nonzero_pivot_generic(long * pivot_row, gr_mat_t mat, long start_row, long end_row, long column, gr_ctx_t ctx)
    int gr_mat_find_nonzero_pivot(long * pivot_row, gr_mat_t mat, long start_row, long end_row, long column, gr_ctx_t ctx)
    # Attempts to find a nonzero element in column number *column*
    # of the matrix *mat* in a row between *start_row* (inclusive)
    # and *end_row* (exclusive).
    # On success, sets ``pivot_row`` to the row index and returns
    # ``GR_SUCCESS``. If no nonzero pivot element exists, returns ``GR_DOMAIN``.
    # If no nonzero pivot element exists and zero-testing fails for some
    # element, returns the flag ``GR_UNABLE``.
    # This function may be destructive: any elements that are nontrivially
    # zero but can be certified zero may be overwritten by exact zeros.

    int gr_mat_lu_classical(long * rank, long * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
    int gr_mat_lu_recursive(long * rank, long * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
    int gr_mat_lu(long * rank, long * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
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
    # This currently only works as expected for square matrices.
    # The algorithm can fail if it fails to certify that a pivot
    # element is zero or nonzero, in which case the correct rank
    # cannot be determined. It can also fail if a pivot element
    # is not invertible. In these cases the ``GR_UNABLE`` and/or
    # ``GR_DOMAIN`` flags will be returned. On failure,
    # the data in the output variables
    # ``rank``, ``P`` and ``LU`` will be meaningless.
    # The *classical* version uses iterative Gaussian elimination.
    # The *recursive* version uses a block recursive algorithm
    # to take advantage of fast matrix multiplication.

    int gr_mat_fflu(long * rank, long * P, gr_mat_t LU, gr_ptr den, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
    # Similar to :func:`gr_mat_lu`, but computes a fraction-free
    # LU decomposition using the Bareiss algorithm.
    # The denominator is written to *den*.

    int gr_mat_nonsingular_solve_tril_classical(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
    int gr_mat_nonsingular_solve_tril_recursive(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
    int gr_mat_nonsingular_solve_tril(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx)
    int gr_mat_nonsingular_solve_triu_classical(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx)
    int gr_mat_nonsingular_solve_triu_recursive(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx)
    int gr_mat_nonsingular_solve_triu(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx)
    # Solves the lower triangular system `LX = B` or the upper triangular system
    # `UX = B`, respectively. Division by the the diagonal entries must
    # be possible; if not a division fails, ``GR_DOMAIN`` is returned
    # even if the system is solvable.
    # If *unit* is set, the main diagonal of *L* or *U*
    # is taken to consist of all ones, and in that case the actual entries on
    # the diagonal are not read at all and can contain other data.
    # The *classical* versions perform the computations iteratively while the
    # *recursive* versions perform the computations in a block recursive
    # way to benefit from fast matrix multiplication. The default versions
    # choose an algorithm automatically.

    int gr_mat_nonsingular_solve_fflu(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
    int gr_mat_nonsingular_solve_lu(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
    int gr_mat_nonsingular_solve(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
    # Solves `AX = B`. If *A* is not invertible,
    # returns ``GR_DOMAIN`` even if the system has a solution.

    int gr_mat_nonsingular_solve_fflu_precomp(gr_mat_t X, const long * perm, const gr_mat_t LU, const gr_mat_t B, gr_ctx_t ctx)
    int gr_mat_nonsingular_solve_lu_precomp(gr_mat_t X, const long * perm, const gr_mat_t LU, const gr_mat_t B, gr_ctx_t ctx)
    # Solves `AX = B` given a precomputed FFLU or LU factorization of *A*.

    int gr_mat_nonsingular_solve_den_fflu(gr_mat_t X, gr_ptr den, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
    int gr_mat_nonsingular_solve_den(gr_mat_t X, gr_ptr den, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
    # Solves `AX = B` over the fraction field of the present ring
    # (assumed to be an integral domain), returning `X` with
    # an implied denominator *den*.
    # If *A* is not invertible over the fraction field, returns
    # ``GR_DOMAIN`` even if the system has a solution.

    int gr_mat_solve_field(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
    # Solves `AX = B` where *A* is not necessarily square and not necessarily
    # invertible. Assuming that the ring is a field, a return value of
    # ``GR_DOMAIN`` indicates that the system has no solution.
    # If there are multiple solutions, an arbitrary solution is returned.

    int gr_mat_det_fflu(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_det_berkowitz(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_det_lu(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_det_cofactor(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_det_generic_field(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
    int gr_mat_det_generic_integral_domain(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
    int gr_mat_det_generic(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
    int gr_mat_det(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
    # Sets *res* to the determinant of the square matrix *mat*.
    # Various algorithms are available:
    # * The *berkowitz* version uses the division-free Berkowitz algorithm
    # performing `O(n^4)` operations. Since no zero tests are required, it
    # is guaranteed to succeed if the ring arithmetic succeeds.
    # * The *cofactor* version performs cofactor expansion. This is currently
    # only supported for matrices up to size 4, and for larger
    # matrices returns the ``GR_UNABLE`` flag.
    # * The *lu* and *fflu* versions use rational LU decomposition
    # and fraction-free LU decomposition (Bareiss algorithm) respectively,
    # requiring `O(n^3)` operations. These algorithms can fail if zero
    # certification or inversion fails, in which case the ``GR_UNABLE``
    # flag is returned.
    # * The *generic*, *generic_field* and *generic_integral_domain*
    # versions choose an appropriate algorithm for a generic ring
    # depending on the availability of division.
    # * The *default* method can be overloaded.
    # If the matrix is not square, ``GR_DOMAIN`` is returned.

    int gr_mat_trace(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
    # Sets *res* to the trace (sum of entries on the main diagonal) of
    # the square matrix *mat*.
    # If the matrix is not square, ``GR_DOMAIN`` is returned.

    int gr_mat_rank_fflu(long * rank, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_rank_lu(long * rank, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_rank(long * rank, const gr_mat_t mat, gr_ctx_t ctx)
    # Sets *res* to the rank of *mat*.
    # The default method returns ``GR_DOMAIN`` if the element ring
    # is not an integral domain, in which case the usual rank is
    # not well-defined. The *fflu* and *lu* variants currently do
    # not check the element domain, and simply return this flag if they
    # encounter an impossible inverse in the execution of the
    # respective algorithms.

    int gr_mat_rref_lu(long * rank, gr_mat_t R, const gr_mat_t A, gr_ctx_t ctx)
    int gr_mat_rref_fflu(long * rank, gr_mat_t R, const gr_mat_t A, gr_ctx_t ctx)
    int gr_mat_rref(long * rank, gr_mat_t R, const gr_mat_t A, gr_ctx_t ctx)
    # Sets *R* to the reduced row echelon form of *A*, also setting
    # *rank* to its rank.

    int gr_mat_rref_den_fflu(long * rank, gr_mat_t R, gr_ptr den, const gr_mat_t A, gr_ctx_t ctx)
    int gr_mat_rref_den(long * rank, gr_mat_t R, gr_ptr den, const gr_mat_t A, gr_ctx_t ctx)
    # Like *rref*, but computes the reduced row echelon multiplied
    # by a common (not necessarily minimal) denominator which is written
    # to *den*. This can be used to compute the rref over an integral
    # domain which is not a field.

    int gr_mat_nullspace(gr_mat_t X, const gr_mat_t A, gr_ctx_t ctx)
    # Sets *X* to a basis for the (right) nullspace of *A*.
    # On success, the output matrix will be resized to the correct
    # number of columns.
    # The basis is not guaranteed to be presented in a
    # canonical or minimal form.
    # If the ring is not a field, this is implied to compute a nullspace
    # basis over the fraction field. The result may be meaningless
    # if the ring is not an integral domain.

    int gr_mat_inv(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
    # Sets *res* to the inverse of *mat*, computed by solving
    # `A A^{-1} = I`.
    # Returns ``GR_DOMAIN`` if it can be determined that *mat* is not
    # invertible over the present ring (warning: this may not work
    # over non-integral domains). If invertibility cannot be proved,
    # returns ``GR_UNABLE``.
    # To compute the inverse over the fraction field, one may use
    # :func:`gr_mat_nonsingular_solve_den` or :func:`gr_mat_adjugate`.

    int gr_mat_adjugate_charpoly(gr_mat_t adj, gr_ptr det, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_adjugate_cofactor(gr_mat_t adj, gr_ptr det, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_adjugate(gr_mat_t adj, gr_ptr det, const gr_mat_t mat, gr_ctx_t ctx)
    # Sets *adj* to the adjugate matrix of *mat*, simultaneously
    # setting *det* to the determinant of *mat*. We have
    # `\operatorname{adj}(A) A = A \operatorname{adj}(A) = \det(A) I`,
    # and `A^{-1} = \operatorname{adj}(A) / \det(A)` when *A*
    # is invertible.
    # The *cofactor* version uses cofactor expansion, requiring the
    # evaluation of `n^2` determinants.
    # The *charpoly* version computes and then evaluates the
    # characteristic polynomial, requiring `O(n^{1/2})`
    # matrix multiplications plus `O(n^3)` or `O(n^4)` operations
    # for the characteristic polynomial itself depending on the
    # algorithm used.

    int _gr_mat_charpoly(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_charpoly(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx)
    # Computes the characteristic polynomial using a default
    # algorithm choice. The
    # underscore method assumes that *res* is a preallocated
    # array of `n + 1` coefficients.

    int _gr_mat_charpoly_berkowitz(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_charpoly_berkowitz(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx)
    # Sets *res* to the characteristic polynomial of the square matrix
    # *mat*, computed using the division-free Berkowitz algorithm.
    # The number of operations is `O(n^4)` where *n* is the
    # size of the matrix.

    int _gr_mat_charpoly_danilevsky_inplace(gr_ptr res, gr_mat_t mat, gr_ctx_t ctx)
    int _gr_mat_charpoly_danilevsky(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_charpoly_danilevsky(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx)
    int _gr_mat_charpoly_gauss(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_charpoly_gauss(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx)
    int _gr_mat_charpoly_householder(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_charpoly_householder(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx)
    # Sets *res* to the characteristic polynomial of the square matrix
    # *mat*, computed using the Danilevsky algorithm,
    # Hessenberg reduction using Gaussian elimination,
    # and Hessenberg reduction using Householder reflections.
    # The number of operations of each method is `O(n^3)` where *n* is the
    # size of the matrix. The *inplace* version overwrites the input matrix.
    # These methods require divisions and can therefore fail when the
    # ring is not a field. They also require zero tests.
    # The *householder* version also requires square roots.
    # The flags ``GR_UNABLE`` or ``GR_DOMAIN`` are returned when
    # an impossible division or square root
    # is encountered or when a comparison cannot be performed.

    int _gr_mat_charpoly_faddeev(gr_ptr res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_charpoly_faddeev(gr_poly_t res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx)
    int _gr_mat_charpoly_faddeev_bsgs(gr_ptr res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_charpoly_faddeev_bsgs(gr_poly_t res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx)
    # Sets *res* to the characteristic polynomial of the square matrix
    # *mat*, computed using the Faddeev-LeVerrier algorithm.
    # If the optional output argument *adj* is not *NULL*, it is
    # set to the adjugate matrix, which is computed free of charge.
    # The *bsgs* version uses a baby-step giant-step strategy,
    # also known as the Preparata-Sarwate algorithm.
    # This reduces the complexity from `O(n^4)` to `O(n^{3.5})` operations
    # at the cost of requiring `n^{0.5}` temporary matrices to be
    # stored.
    # This method requires divisions by small integers and can
    # therefore fail (returning the ``GR_UNABLE`` or ``GR_DOMAIN`` flags)
    # in finite characteristic or when the underlying ring does
    # not implement a division algorithm.

    int _gr_mat_charpoly_from_hessenberg(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_charpoly_from_hessenberg(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx)
    # Sets *res* to the characteristic polynomial of the square matrix
    # *mat*, which is assumed to be in Hessenberg form (this is
    # currently not checked).

    int gr_mat_minpoly_field(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx)
    # Compute the minimal polynomial of the matrix *mat*.
    # The algorithm assumes that the coefficient ring is a field.

    int gr_mat_apply_row_similarity(gr_mat_t M, long r, gr_ptr d, gr_ctx_t ctx)
    # Applies an elementary similarity transform to the `n\times n` matrix `M`
    # in-place.
    # If `P` is the `n\times n` identity matrix the zero entries of whose row
    # `r` (`0`-indexed) have been replaced by `d`, this transform is equivalent
    # to `M = P^{-1}MP`.
    # Similarity transforms preserve the determinant, characteristic polynomial
    # and minimal polynomial.

    int gr_mat_eigenvalues(gr_vec_t lmbda, gr_vec_t mult, const gr_mat_t mat, int flags, gr_ctx_t ctx)
    int gr_mat_eigenvalues_other(gr_vec_t lmbda, gr_vec_t mult, const gr_mat_t mat, gr_ctx_t mat_ctx, int flags, gr_ctx_t ctx)
    # Finds all eigenvalues of the given matrix in the ring defined by *ctx*,
    # storing the eigenvalues without duplication in *lambda* (a vector with
    # elements of type ``ctx``) and the corresponding multiplicities in
    # *mult* (a vector with elements of type ``fmpz``).
    # The interface is essentially the same as that of
    # :func:`gr_poly_roots`; see its documentation for details.

    int gr_mat_diagonalization_precomp(gr_vec_t D, gr_mat_t L, gr_mat_t R, const gr_mat_t A, const gr_vec_t eigenvalues, const gr_vec_t mult, gr_ctx_t ctx)
    int gr_mat_diagonalization_generic(gr_vec_t D, gr_mat_t L, gr_mat_t R, const gr_mat_t A, int flags, gr_ctx_t ctx)
    int gr_mat_diagonalization(gr_vec_t D, gr_mat_t L, gr_mat_t R, const gr_mat_t A, int flags, gr_ctx_t ctx)
    # Computes a diagonalization `LAR = D` given a square matrix `A`,
    # where `D` is a diagonal matrix (returned as a vector) of the eigenvalues
    # repeated according to their multiplicities,
    # `L` is a matrix of left eigenvectors,
    # and `R` is a matrix of right eigenvectors,
    # normalized such that `L = R^{-1}`.
    # This implies that `A = RDL = RDR^{-1}`.
    # Either `L` or `R` (or both) can be set to ``NULL`` to omit computing
    # the respective matrix.
    # If the matrix has entries in a field then a return flag
    # of ``GR_DOMAIN`` indicates that the matrix is non-diagonalizable
    # over this field.
    # The *precomp* version requires as input a precomputed set of eigenvalues
    # with corresponding multiplicities, which can be computed
    # with :func:`gr_mat_eigenvalues`.

    int gr_mat_set_jordan_blocks(gr_mat_t mat, const gr_vec_t lmbda, long num_blocks, long * block_lambda, long * block_size, gr_ctx_t ctx)
    int gr_mat_jordan_blocks(gr_vec_t lmbda, long * num_blocks, long * block_lambda, long * block_size, const gr_mat_t A, gr_ctx_t ctx)
    int gr_mat_jordan_transformation(gr_mat_t mat, const gr_vec_t lmbda, long num_blocks, long * block_lambda, long * block_size, const gr_mat_t A, gr_ctx_t ctx)
    int gr_mat_jordan_form(gr_mat_t J, gr_mat_t P, const gr_mat_t A, gr_ctx_t ctx)

    int gr_mat_exp_jordan(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)
    int gr_mat_exp(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)

    int gr_mat_log_jordan(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)
    int gr_mat_log(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)

    truth_t gr_mat_is_hessenberg(const gr_mat_t mat, gr_ctx_t ctx)
    # Returns whether *mat* is in upper Hessenberg form.

    int gr_mat_hessenberg_gauss(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_hessenberg_householder(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
    int gr_mat_hessenberg(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
    # Sets *res* to an upper Hessenberg form of *mat*.
    # The *gauss* version uses Gaussian elimination.
    # The *householder* version uses Householder reflections.
    # These methods require divisions and zero testing
    # and can therefore fail (returning ``GR_UNABLE`` or ``GR_DOMAIN``)
    # when the ring is not a field.
    # The *householder* version additionally requires complex
    # conjugation and the ability to compute square roots.

    int gr_mat_randtest(gr_mat_t res, flint_rand_t state, gr_ctx_t ctx)
    # Sets *res* to a random matrix. The distribution is nonuniform.

    int gr_mat_randops(gr_mat_t mat, flint_rand_t state, long count, gr_ctx_t ctx)
    # Randomises *mat* in-place by performing elementary row or column
    # operations. More precisely, at most *count* random additions or
    # subtractions of distinct rows and columns will be performed.

    int gr_mat_randpermdiag(int * parity, gr_mat_t mat, flint_rand_t state, gr_ptr diag, long n, gr_ctx_t ctx)
    # Sets mat to a random permutation of the diagonal matrix with *n* leading entries given by
    # the vector ``diag``. Returns ``GR_DOMAIN`` if the main diagonal of ``mat``
    # does not have room for at least *n* entries.
    # The parity (0 or 1) of the permutation is written to ``parity``.

    int gr_mat_randrank(gr_mat_t mat, flint_rand_t state, long rank, gr_ctx_t ctx)
    # Sets ``mat`` to a random sparse matrix with the given rank, having exactly as many
    # non-zero elements as the rank. The matrix can be transformed into a dense matrix
    # with unchanged rank by subsequently calling :func:`gr_mat_randops`.
    # This operation only makes sense over integral domains (currently not checked).

    int gr_mat_ones(gr_mat_t res, gr_ctx_t ctx)
    # Sets all entries in *res* to one.

    int gr_mat_pascal(gr_mat_t res, int triangular, gr_ctx_t ctx)
    # Sets *res* to a Pascal matrix, whose entries are binomial coefficients.
    # If *triangular* is 0, constructs a full symmetric matrix
    # with the rows of Pascal's triangle as successive antidiagonals.
    # If *triangular* is 1, constructs the upper triangular matrix with
    # the rows of Pascal's triangle as columns, and if *triangular* is -1,
    # constructs the lower triangular matrix with the rows of Pascal's
    # triangle as rows.

    int gr_mat_stirling(gr_mat_t res, int kind, gr_ctx_t ctx)
    # Sets *res* to a Stirling matrix, whose entries are Stirling numbers.
    # If *kind* is 0, the entries are set to the unsigned Stirling numbers
    # of the first kind. If *kind* is 1, the entries are set to the signed
    # Stirling numbers of the first kind. If *kind* is 2, the entries are
    # set to the Stirling numbers of the second kind.

    int gr_mat_hilbert(gr_mat_t res, gr_ctx_t ctx)
    # Sets *res* to the Hilbert matrix, which has entries `1/(i+j+1)`
    # for `i, j \ge 0`.

    int gr_mat_hadamard(gr_mat_t res, gr_ctx_t ctx)
    # If possible, sets *res* to a Hadamard matrix of the provided size
    # and returns ``GR_SUCCESS``. Returns ``GR_DOMAIN``
    # if no Hadamard matrix of the given size exists,
    # and ``GR_UNABLE`` if the implementation does
    # not know how to construct a Hadamard matrix of the given
    # size.
    # A Hadamard matrix of size *n* can only exist if *n* is 0, 1, 2,
    # or a multiple of 4. It is not known whether a
    # Hadamard matrix exists for every size that is a multiple of 4.
    # This function uses the Paley construction, which
    # succeeds for all *n* of the form `n = 2^e` or `n = 2^e (q + 1)` where
    # *q* is an odd prime power. Orders *n* for which Hadamard matrices are
    # known to exist but for which this construction fails are
    # 92, 116, 156, ... (OEIS A046116).

    int gr_mat_reduce_row(long * column, gr_mat_t A, long * P, long * L, long m, gr_ctx_t ctx)
    # Reduce row n of the matrix `A`, assuming the prior rows are in Gauss
    # form. However those rows may not be in order. The entry `i` of the array
    # `P` is the row of `A` which has a pivot in the `i`-th column. If no such
    # row exists, the entry of `P` will be `-1`. The function sets *column* to the column
    # in which the `n`-th row has a pivot after reduction. This will always be
    # chosen to be the first available column for a pivot from the left. This
    # information is also updated in `P`. Entry `i` of the array `L` contains the
    # number of possibly nonzero columns of `A` row `i`. This speeds up reduction
    # in the case that `A` is chambered on the right. Otherwise the entries of
    # `L` can all be set to the number of columns of `A`. We require the entries
    # of `L` to be monotonic increasing.
