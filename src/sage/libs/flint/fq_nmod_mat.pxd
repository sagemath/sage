# distutils: libraries = flint
# distutils: depends = flint/fq_nmod_mat.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void fq_nmod_mat_init(fq_nmod_mat_t mat, slong rows, slong cols, const fq_nmod_ctx_t ctx)
    # Initialises ``mat`` to a ``rows``-by-``cols`` matrix with
    # coefficients in `\mathbf{F}_{q}` given by ``ctx``. All elements
    # are set to zero.

    void fq_nmod_mat_init_set(fq_nmod_mat_t mat, const fq_nmod_mat_t src, const fq_nmod_ctx_t ctx)
    # Initialises ``mat`` and sets its dimensions and elements to
    # those of ``src``.

    void fq_nmod_mat_clear(fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    # Clears the matrix and releases any memory it used. The matrix
    # cannot be used again until it is initialised. This function must be
    # called exactly once when finished using an :type:`fq_nmod_mat_t` object.

    void fq_nmod_mat_set(fq_nmod_mat_t mat, const fq_nmod_mat_t src, const fq_nmod_ctx_t ctx)
    # Sets ``mat`` to a copy of ``src``. It is assumed
    # that ``mat`` and ``src`` have identical dimensions.

    fq_nmod_struct * fq_nmod_mat_entry(const fq_nmod_mat_t mat, slong i, slong j)
    # Directly accesses the entry in ``mat`` in row `i` and column `j`,
    # indexed from zero. No bounds checking is performed.

    void fq_nmod_mat_entry_set(fq_nmod_mat_t mat, slong i, slong j, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    # Sets the entry in ``mat`` in row `i` and column `j` to ``x``.

    slong fq_nmod_mat_nrows(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    # Returns the number of rows in ``mat``.

    slong fq_nmod_mat_ncols(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    # Returns the number of columns in ``mat``.

    void fq_nmod_mat_swap(fq_nmod_mat_t mat1, fq_nmod_mat_t mat2, const fq_nmod_ctx_t ctx)
    # Swaps two matrices. The dimensions of ``mat1`` and ``mat2``
    # are allowed to be different.

    void fq_nmod_mat_swap_entrywise(fq_nmod_mat_t mat1, fq_nmod_mat_t mat2, const fq_nmod_ctx_t ctx)
    # Swaps two matrices by swapping the individual entries rather than swapping
    # the contents of the structs.

    void fq_nmod_mat_zero(fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    # Sets all entries of ``mat`` to 0.

    void fq_nmod_mat_one(fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    # Sets all diagonal entries of ``mat`` to 1 and all other entries to 0.

    void fq_nmod_mat_swap_rows(fq_nmod_mat_t mat, slong * perm, slong r, slong s, const fq_nmod_ctx_t ctx)
    # Swaps rows ``r`` and ``s`` of ``mat``.  If ``perm`` is non-``NULL``, the
    # permutation of the rows will also be applied to ``perm``.

    void fq_nmod_mat_swap_cols(fq_nmod_mat_t mat, slong * perm, slong r, slong s, const fq_nmod_ctx_t ctx)
    # Swaps columns ``r`` and ``s`` of ``mat``.  If ``perm`` is non-``NULL``, the
    # permutation of the columns will also be applied to ``perm``.

    void fq_nmod_mat_invert_rows(fq_nmod_mat_t mat, slong * perm, const fq_nmod_ctx_t ctx)
    # Swaps rows ``i`` and ``r - i`` of ``mat`` for ``0 <= i < r/2``, where
    # ``r`` is the number of rows of ``mat``. If ``perm`` is non-``NULL``, the
    # permutation of the rows will also be applied to ``perm``.

    void fq_nmod_mat_invert_cols(fq_nmod_mat_t mat, slong * perm, const fq_nmod_ctx_t ctx)
    # Swaps columns ``i`` and ``c - i`` of ``mat`` for ``0 <= i < c/2``, where
    # ``c`` is the number of columns of ``mat``. If ``perm`` is non-``NULL``, the
    # permutation of the columns will also be applied to ``perm``.

    void fq_nmod_mat_set_nmod_mat(fq_nmod_mat_t mat1, const nmod_mat_t mat2, const fq_nmod_ctx_t ctx)
    # Sets the matrix ``mat1`` to the matrix ``mat2``.

    void fq_nmod_mat_set_fmpz_mod_mat(fq_nmod_mat_t mat1, const fmpz_mod_mat_t mat2, const fq_nmod_ctx_t ctx)
    # Sets the matrix ``mat1`` to the matrix ``mat2``.

    void fq_nmod_mat_concat_vertical(fq_nmod_mat_t res, const fq_nmod_mat_t mat1, const fq_nmod_mat_t mat2, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to vertical concatenation of (``mat1``, ``mat2``) in that order. Matrix dimensions : ``mat1`` : `m \times n`, ``mat2`` : `k \times n`, ``res`` : `(m + k) \times n`.

    void fq_nmod_mat_concat_horizontal(fq_nmod_mat_t res, const fq_nmod_mat_t mat1, const fq_nmod_mat_t mat2, const fq_nmod_ctx_t ctx)
    # Sets ``res`` to horizontal concatenation of (``mat1``, ``mat2``) in that order. Matrix dimensions : ``mat1`` : `m \times n`, ``mat2`` : `m \times k`, ``res``  : `m \times (n + k)`.

    int fq_nmod_mat_print_pretty(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    # Pretty-prints ``mat`` to ``stdout``. A header is printed
    # followed by the rows enclosed in brackets.

    int fq_nmod_mat_fprint_pretty(FILE * file, const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    # Pretty-prints ``mat`` to ``file``. A header is printed
    # followed by the rows enclosed in brackets.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    int fq_nmod_mat_print(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    # Prints ``mat`` to ``stdout``. A header is printed followed
    # by the rows enclosed in brackets.

    int fq_nmod_mat_fprint(FILE * file, const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    # Prints ``mat`` to ``file``. A header is printed followed by
    # the rows enclosed in brackets.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    void fq_nmod_mat_window_init(fq_nmod_mat_t window, const fq_nmod_mat_t mat, slong r1, slong c1, slong r2, slong c2, const fq_nmod_ctx_t ctx)
    # Initializes the matrix ``window`` to be an ``r2 - r1`` by
    # ``c2 - c1`` submatrix of ``mat`` whose ``(0,0)`` entry
    # is the ``(r1, c1)`` entry of ``mat``.  The memory for the
    # elements of ``window`` is shared with ``mat``.

    void fq_nmod_mat_window_clear(fq_nmod_mat_t window, const fq_nmod_ctx_t ctx)
    # Clears the matrix ``window`` and releases any memory that it
    # uses.  Note that the memory to the underlying matrix that
    # ``window`` points to is not freed.

    void fq_nmod_mat_randtest(fq_nmod_mat_t mat, flint_rand_t state, const fq_nmod_ctx_t ctx)
    # Sets the elements of ``mat`` to random elements of
    # `\mathbf{F}_{q}`, given by ``ctx``.

    int fq_nmod_mat_randpermdiag(fq_nmod_mat_t mat, flint_rand_t state, fq_nmod_struct * diag, slong n, const fq_nmod_ctx_t ctx)
    # Sets ``mat`` to a random permutation of the diagonal matrix
    # with `n` leading entries given by the vector ``diag``. It is
    # assumed that the main diagonal of ``mat`` has room for at
    # least `n` entries.
    # Returns `0` or `1`, depending on whether the permutation is even
    # or odd respectively.

    void fq_nmod_mat_randrank(fq_nmod_mat_t mat, flint_rand_t state, slong rank, const fq_nmod_ctx_t ctx)
    # Sets ``mat`` to a random sparse matrix with the given rank,
    # having exactly as many non-zero elements as the rank, with the
    # non-zero elements being uniformly random elements of
    # `\mathbf{F}_{q}`.
    # The matrix can be transformed into a dense matrix with unchanged
    # rank by subsequently calling :func:`fq_nmod_mat_randops`.

    void fq_nmod_mat_randops(fq_nmod_mat_t mat, slong count, flint_rand_t state, const fq_nmod_ctx_t ctx)
    # Randomises ``mat`` by performing elementary row or column
    # operations. More precisely, at most ``count`` random additions
    # or subtractions of distinct rows and columns will be performed.
    # This leaves the rank (and for square matrices, determinant)
    # unchanged.

    void fq_nmod_mat_randtril(fq_nmod_mat_t mat, flint_rand_t state, int unit, const fq_nmod_ctx_t ctx)
    # Sets ``mat`` to a random lower triangular matrix. If
    # ``unit`` is 1, it will have ones on the main diagonal,
    # otherwise it will have random nonzero entries on the main
    # diagonal.

    void fq_nmod_mat_randtriu(fq_nmod_mat_t mat, flint_rand_t state, int unit, const fq_nmod_ctx_t ctx)
    # Sets ``mat`` to a random upper triangular matrix. If
    # ``unit`` is 1, it will have ones on the main diagonal,
    # otherwise it will have random nonzero entries on the main
    # diagonal.

    bint fq_nmod_mat_equal(const fq_nmod_mat_t mat1, const fq_nmod_mat_t mat2, const fq_nmod_ctx_t ctx)
    # Returns nonzero if mat1 and mat2 have the same dimensions and elements,
    # and zero otherwise.

    bint fq_nmod_mat_is_zero(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    # Returns a non-zero value if all entries ``mat`` are zero, and
    # otherwise returns zero.

    bint fq_nmod_mat_is_one(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    # Returns a non-zero value if all entries ``mat`` are zero except the
    # diagonal entries which must be one, otherwise returns zero.

    bint fq_nmod_mat_is_empty(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    # Returns a non-zero value if the number of rows or the number of
    # columns in ``mat`` is zero, and otherwise returns zero.

    bint fq_nmod_mat_is_square(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    # Returns a non-zero value if the number of rows is equal to the
    # number of columns in ``mat``, and otherwise returns zero.

    void fq_nmod_mat_add(fq_nmod_mat_t C, const fq_nmod_mat_t A, const fq_nmod_mat_t B,  const fq_nmod_ctx_t ctx)
    # Computes `C = A + B`. Dimensions must be identical.

    void fq_nmod_mat_sub(fq_nmod_mat_t C, const fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
    # Computes `C = A - B`. Dimensions must be identical.

    void fq_nmod_mat_neg(fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
    # Sets `B = -A`. Dimensions must be identical.

    void fq_nmod_mat_mul(fq_nmod_mat_t C, const fq_nmod_mat_t A, const fq_nmod_mat_t B,  const fq_nmod_ctx_t ctx)
    # Sets `C = AB`. Dimensions must be compatible for matrix
    # multiplication. Aliasing is allowed. This function automatically chooses
    # between classical and KS multiplication.

    void fq_nmod_mat_mul_classical(fq_nmod_mat_t C, const fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
    # Sets `C = AB`. Dimensions must be compatible for matrix multiplication.
    # `C` is not allowed to be aliased with `A` or `B`. Uses classical
    # matrix multiplication.

    void fq_nmod_mat_mul_KS(fq_nmod_mat_t C, const fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
    # Sets `C = AB`. Dimensions must be compatible for matrix
    # multiplication.  `C` is not allowed to be aliased with `A` or
    # `B`. Uses Kronecker substitution to perform the multiplication
    # over the integers.

    void fq_nmod_mat_submul(fq_nmod_mat_t D, const fq_nmod_mat_t C, const fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
    # Sets `D = C + AB`. `C` and `D` may be aliased with each other but
    # not with `A` or `B`.

    void fq_nmod_mat_mul_vec(fq_nmod_struct * c, const fq_nmod_mat_t A, const fq_nmod_struct * b, slong blen, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_mul_vec_ptr(fq_nmod_struct * const * c, const fq_nmod_mat_t A, const fq_nmod_struct * const * b, slong blen, const fq_nmod_ctx_t ctx)
    # Compute a matrix-vector product of ``A`` and ``(b, blen)`` and store the result in ``c``.
    # The vector ``(b, blen)`` is either truncated or zero-extended to the number of columns of ``A``.
    # The number entries written to ``c`` is always equal to the number of rows of ``A``.

    void fq_nmod_mat_vec_mul(fq_nmod_struct * c, const fq_nmod_struct * a, slong alen, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_vec_mul_ptr(fq_nmod_struct * const * c, const fq_nmod_struct * const * a, slong alen, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
    # Compute a vector-matrix product of ``(a, alen)`` and ``B`` and and store the result in ``c``.
    # The vector ``(a, alen)`` is either truncated or zero-extended to the number of rows of ``B``.
    # The number entries written to ``c`` is always equal to the number of columns of ``B``.

    int fq_nmod_mat_inv(fq_nmod_mat_t B, fq_nmod_mat_t A, const fq_nmod_ctx_t ctx)
    # Sets `B = A^{-1}` and returns `1` if `A` is invertible. If `A` is singular,
    # returns `0` and sets the elements of `B` to undefined values.
    # `A` and `B` must be square matrices with the same dimensions.

    slong fq_nmod_mat_lu(slong * P, fq_nmod_mat_t A, int rank_check, const fq_nmod_ctx_t ctx)
    # Computes a generalised LU decomposition `LU = PA` of a given
    # matrix `A`, returning the rank of `A`.
    # If `A` is a nonsingular square matrix, it will be overwritten with
    # a unit diagonal lower triangular matrix `L` and an upper
    # triangular matrix `U` (the diagonal of `L` will not be stored
    # explicitly).
    # If `A` is an arbitrary matrix of rank `r`, `U` will be in row
    # echelon form having `r` nonzero rows, and `L` will be lower
    # triangular but truncated to `r` columns, having implicit ones on
    # the `r` first entries of the main diagonal. All other entries will
    # be zero.
    # If a nonzero value for ``rank_check`` is passed, the function
    # will abandon the output matrix in an undefined state and return 0
    # if `A` is detected to be rank-deficient.
    # This function calls ``fq_nmod_mat_lu_recursive``.

    slong fq_nmod_mat_lu_classical(slong * P, fq_nmod_mat_t A, int rank_check, const fq_nmod_ctx_t ctx)
    # Computes a generalised LU decomposition `LU = PA` of a given
    # matrix `A`, returning the rank of `A`. The behavior of this
    # function is identical to that of ``fq_nmod_mat_lu``. Uses Gaussian
    # elimination.

    slong fq_nmod_mat_lu_recursive(slong * P, fq_nmod_mat_t A, int rank_check, const fq_nmod_ctx_t ctx)
    # Computes a generalised LU decomposition `LU = PA` of a given
    # matrix `A`, returning the rank of `A`. The behavior of this
    # function is identical to that of ``fq_nmod_mat_lu``. Uses recursive
    # block decomposition, switching to classical Gaussian elimination
    # for sufficiently small blocks.

    slong fq_nmod_mat_rref(fq_nmod_mat_t A, const fq_nmod_ctx_t ctx)
    # Puts `A` in reduced row echelon form and returns the rank of `A`.
    # The rref is computed by first obtaining an unreduced row echelon
    # form via LU decomposition and then solving an additional
    # triangular system.

    slong fq_nmod_mat_reduce_row(fq_nmod_mat_t A, slong * P, slong * L, slong n, const fq_nmod_ctx_t ctx)
    # Reduce row n of the matrix `A`, assuming the prior rows are in Gauss
    # form. However those rows may not be in order. The entry `i` of the array
    # `P` is the row of `A` which has a pivot in the `i`-th column. If no such
    # row exists, the entry of `P` will be `-1`. The function returns the column
    # in which the `n`-th row has a pivot after reduction. This will always be
    # chosen to be the first available column for a pivot from the left. This
    # information is also updated in `P`. Entry `i` of the array `L` contains the
    # number of possibly nonzero columns of `A` row `i`. This speeds up reduction
    # in the case that `A` is chambered on the right. Otherwise the entries of
    # `L` can all be set to the number of columns of `A`. We require the entries
    # of `L` to be monotonic increasing.

    void fq_nmod_mat_solve_tril(fq_nmod_mat_t X, const fq_nmod_mat_t L, const fq_nmod_mat_t B, int unit, const fq_nmod_ctx_t ctx)
    # Sets `X = L^{-1} B` where `L` is a full rank lower triangular
    # square matrix. If ``unit`` = 1, `L` is assumed to have ones on
    # its main diagonal, and the main diagonal will not be read.  `X`
    # and `B` are allowed to be the same matrix, but no other aliasing
    # is allowed. Automatically chooses between the classical and
    # recursive algorithms.

    void fq_nmod_mat_solve_tril_classical(fq_nmod_mat_t X, const fq_nmod_mat_t L, const fq_nmod_mat_t B, int unit, const fq_nmod_ctx_t ctx)
    # Sets `X = L^{-1} B` where `L` is a full rank lower triangular
    # square matrix. If ``unit`` = 1, `L` is assumed to have ones on
    # its main diagonal, and the main diagonal will not be read.  `X`
    # and `B` are allowed to be the same matrix, but no other aliasing
    # is allowed. Uses forward substitution.

    void fq_nmod_mat_solve_tril_recursive(fq_nmod_mat_t X, const fq_nmod_mat_t L, const fq_nmod_mat_t B, int unit, const fq_nmod_ctx_t ctx)
    # Sets `X = L^{-1} B` where `L` is a full rank lower triangular
    # square matrix. If ``unit`` = 1, `L` is assumed to have ones on
    # its main diagonal, and the main diagonal will not be read.  `X`
    # and `B` are allowed to be the same matrix, but no other aliasing
    # is allowed.
    # Uses the block inversion formula
    # .. math ::
    # \begin{pmatrix} A & 0 \\ C & D \end{pmatrix}^{-1}
    # \begin{pmatrix} X \\ Y \end{pmatrix} =
    # \begin{pmatrix} A^{-1} X \\ D^{-1} ( Y - C A^{-1} X ) \end{pmatrix}
    # to reduce the problem to matrix multiplication and triangular
    # solving of smaller systems.

    void fq_nmod_mat_solve_triu(fq_nmod_mat_t X, const fq_nmod_mat_t U, const fq_nmod_mat_t B, int unit, const fq_nmod_ctx_t ctx)
    # Sets `X = U^{-1} B` where `U` is a full rank upper triangular
    # square matrix. If ``unit`` = 1, `U` is assumed to have ones on
    # its main diagonal, and the main diagonal will not be read.  `X`
    # and `B` are allowed to be the same matrix, but no other aliasing
    # is allowed. Automatically chooses between the classical and
    # recursive algorithms.

    void fq_nmod_mat_solve_triu_classical(fq_nmod_mat_t X, const fq_nmod_mat_t U, const fq_nmod_mat_t B, int unit, const fq_nmod_ctx_t ctx)
    # Sets `X = U^{-1} B` where `U` is a full rank upper triangular
    # square matrix. If ``unit`` = 1, `U` is assumed to have ones on
    # its main diagonal, and the main diagonal will not be read.  `X`
    # and `B` are allowed to be the same matrix, but no other aliasing
    # is allowed. Uses forward substitution.

    void fq_nmod_mat_solve_triu_recursive(fq_nmod_mat_t X, const fq_nmod_mat_t U, const fq_nmod_mat_t B, int unit, const fq_nmod_ctx_t ctx)
    # Sets `X = U^{-1} B` where `U` is a full rank upper triangular
    # square matrix. If ``unit`` = 1, `U` is assumed to have ones on
    # its main diagonal, and the main diagonal will not be read.  `X`
    # and `B` are allowed to be the same matrix, but no other aliasing
    # is allowed.
    # Uses the block inversion formula
    # .. math ::
    # \begin{pmatrix} A & B \\ 0 & D \end{pmatrix}^{-1}
    # \begin{pmatrix} X \\ Y \end{pmatrix} =
    # \begin{pmatrix} A^{-1} (X - B D^{-1} Y) \\ D^{-1} Y \end{pmatrix}
    # to reduce the problem to matrix multiplication and triangular
    # solving of smaller systems.

    int fq_nmod_mat_solve(fq_nmod_mat_t X, const fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
    # Solves the matrix-matrix equation `AX = B`.
    # Returns `1` if `A` has full rank; otherwise returns `0` and sets the
    # elements of `X` to undefined values.
    # The matrix `A` must be square.

    int fq_nmod_mat_can_solve(fq_nmod_mat_t X, const fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
    # Solves the matrix-matrix equation `AX = B` over `Fq`.
    # Returns `1` if a solution exists; otherwise returns `0` and sets the
    # elements of `X` to zero. If more than one solution exists, one of the
    # valid solutions is given.
    # There are no restrictions on the shape of `A` and it may be singular.

    void fq_nmod_mat_similarity(fq_nmod_mat_t M, slong r, fq_nmod_t d, const fq_nmod_ctx_t ctx)
    # Applies a similarity transform to the `n\times n` matrix `M` in-place.
    # If `P` is the `n\times n` identity matrix the zero entries of whose row
    # `r` (`0`-indexed) have been replaced by `d`, this transform is equivalent
    # to `M = P^{-1}MP`.
    # Similarity transforms preserve the determinant, characteristic polynomial
    # and minimal polynomial.
    # The value `d` is required to be reduced modulo the modulus of the entries
    # in the matrix.

    void fq_nmod_mat_charpoly_danilevsky(fq_nmod_poly_t p, const fq_nmod_mat_t M, const fq_nmod_ctx_t ctx)
    # Compute the characteristic polynomial `p` of the matrix `M`. The matrix
    # is assumed to be square.

    void fq_nmod_mat_charpoly(fq_nmod_poly_t p, const fq_nmod_mat_t M, const fq_nmod_ctx_t ctx)
    # Compute the characteristic polynomial `p` of the matrix `M`. The matrix
    # is required to be square, otherwise an exception is raised.

    void fq_nmod_mat_minpoly(fq_nmod_poly_t p, const fq_nmod_mat_t M, const fq_nmod_ctx_t ctx)
    # Compute the minimal polynomial `p` of the matrix `M`. The matrix
    # is required to be square, otherwise an exception is raised.
