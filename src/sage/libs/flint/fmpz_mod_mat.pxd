# distutils: libraries = flint
# distutils: depends = flint/fmpz_mod_mat.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    fmpz * fmpz_mod_mat_entry(const fmpz_mod_mat_t mat, slong i, slong j)
    # Return a reference to the element at row ``i`` and column ``j`` of ``mat``.

    void fmpz_mod_mat_set_entry(fmpz_mod_mat_t mat, slong i, slong j, const fmpz_t val)
    # Set the entry at row ``i`` and column ``j`` of ``mat`` to ``val``.

    void fmpz_mod_mat_init(fmpz_mod_mat_t mat, slong rows, slong cols, const fmpz_t n)
    # Initialise ``mat`` as a matrix with the given number of ``rows`` and
    # ``cols`` and modulus ``n``.

    void fmpz_mod_mat_init_set(fmpz_mod_mat_t mat, const fmpz_mod_mat_t src)
    # Initialise ``mat`` and set it equal to the matrix ``src``, including the
    # number of rows and columns and the modulus.

    void fmpz_mod_mat_clear(fmpz_mod_mat_t mat)
    # Clear ``mat`` and release any memory it used.

    slong fmpz_mod_mat_nrows(const fmpz_mod_mat_t mat)

    slong fmpz_mod_mat_ncols(const fmpz_mod_mat_t mat)
    # Return the number of columns of ``mat``.

    void _fmpz_mod_mat_set_mod(fmpz_mod_mat_t mat, const fmpz_t n)
    # Set the modulus of the matrix ``mat`` to ``n``.

    void fmpz_mod_mat_one(fmpz_mod_mat_t mat)
    # Set ``mat`` to the identity matrix (ones down the diagonal).

    void fmpz_mod_mat_zero(fmpz_mod_mat_t mat)
    # Set ``mat`` to the zero matrix.

    void fmpz_mod_mat_swap(fmpz_mod_mat_t mat1, fmpz_mod_mat_t mat2)
    # Efficiently swap the matrices ``mat1`` and ``mat2``.

    void fmpz_mod_mat_swap_entrywise(fmpz_mod_mat_t mat1, fmpz_mod_mat_t mat2)
    # Swaps two matrices by swapping the individual entries rather than swapping
    # the contents of the structs.

    bint fmpz_mod_mat_is_empty(const fmpz_mod_mat_t mat)
    # Return `1` if ``mat`` has either zero rows or columns.

    bint fmpz_mod_mat_is_square(const fmpz_mod_mat_t mat)
    # Return `1` if ``mat`` has the same number of rows and columns.

    void _fmpz_mod_mat_reduce(fmpz_mod_mat_t mat)
    # Reduce all the entries of ``mat`` by the modulus ``n``. This function is
    # only needed internally.

    void fmpz_mod_mat_randtest(fmpz_mod_mat_t mat, flint_rand_t state)
    # Generate a random matrix with the existing dimensions and entries in
    # `[0, n)` where ``n`` is the modulus.

    void fmpz_mod_mat_window_init(fmpz_mod_mat_t window, const fmpz_mod_mat_t mat, slong r1, slong c1, slong r2, slong c2)
    # Initializes the matrix ``window`` to be an ``r2 - r1`` by
    # ``c2 - c1`` submatrix of ``mat`` whose ``(0, 0)`` entry
    # is the ``(r1, c1)`` entry of ``mat``. The memory for the
    # elements of ``window`` is shared with ``mat``.

    void fmpz_mod_mat_window_clear(fmpz_mod_mat_t window)
    # Clears the matrix ``window`` and releases any memory that it
    # uses. Note that the memory to the underlying matrix that
    # ``window`` points to is not freed.

    void fmpz_mod_mat_concat_horizontal(fmpz_mod_mat_t res, const fmpz_mod_mat_t mat1, const fmpz_mod_mat_t mat2)
    # Sets ``res`` to vertical concatenation of (``mat1``, ``mat2``)                            in that order. Matrix dimensions : ``mat1`` : `m \times n`,                               ``mat2`` : `k \times n`, ``res`` : `(m + k) \times n`.

    void fmpz_mod_mat_concat_vertical(fmpz_mod_mat_t res, const fmpz_mod_mat_t mat1, const fmpz_mod_mat_t mat2)
    # Sets ``res`` to horizontal concatenation of (``mat1``, ``mat2``)
    # in that order. Matrix dimensions : ``mat1`` : `m \times n`,
    # ``mat2`` : `m \times k`, ``res``  : `m \times (n + k)`.

    void fmpz_mod_mat_print_pretty(const fmpz_mod_mat_t mat)
    # Prints the given matrix to ``stdout``.  The format is an
    # opening square bracket then on each line a row of the matrix, followed
    # by a closing square bracket. Each row is written as an opening square
    # bracket followed by a space separated list of coefficients followed
    # by a closing square bracket.

    bint fmpz_mod_mat_is_zero(const fmpz_mod_mat_t mat)
    # Return `1` if ``mat`` is the zero matrix.

    void fmpz_mod_mat_set(fmpz_mod_mat_t B, const fmpz_mod_mat_t A)
    # Set ``B`` to equal ``A``.

    void fmpz_mod_mat_transpose(fmpz_mod_mat_t B, const fmpz_mod_mat_t A)
    # Set ``B`` to the transpose of ``A``.

    void fmpz_mod_mat_set_fmpz_mat(fmpz_mod_mat_t A, const fmpz_mat_t B)
    # Set ``A`` to the matrix ``B`` reducing modulo the modulus of ``A``.

    void fmpz_mod_mat_get_fmpz_mat(fmpz_mat_t A, const fmpz_mod_mat_t B)
    # Set ``A`` to a lift of ``B``.

    void fmpz_mod_mat_add(fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B)
    # Set ``C`` to `A + B`.

    void fmpz_mod_mat_sub(fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B)
    # Set ``C`` to `A - B`.

    void fmpz_mod_mat_neg(fmpz_mod_mat_t B, const fmpz_mod_mat_t A)
    # Set ``B`` to `-A`.

    void fmpz_mod_mat_scalar_mul_si(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, slong c)
    # Set ``B`` to `cA` where ``c`` is a constant.

    void fmpz_mod_mat_scalar_mul_ui(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, ulong c)
    # Set ``B`` to `cA` where ``c`` is a constant.

    void fmpz_mod_mat_scalar_mul_fmpz(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, fmpz_t c)
    # Set ``B`` to `cA` where ``c`` is a constant.

    void fmpz_mod_mat_mul(fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B)
    # Set ``C`` to ``A\times B``. The number of rows of ``B`` must match the
    # number of columns of ``A``.

    void _fmpz_mod_mat_mul_classical_threaded_pool_op(fmpz_mod_mat_t D, const fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B, int op, thread_pool_handle * threads, slong num_threads)
    # Set ``D`` to ``A\times B + op*C`` where ``op`` is ``+1``, ``-1`` or ``0``.

    void _fmpz_mod_mat_mul_classical_threaded_op(fmpz_mod_mat_t D, const fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B, int op)
    # Set ``D`` to ``A\times B + op*C`` where ``op`` is ``+1``, ``-1`` or ``0``.

    void fmpz_mod_mat_mul_classical_threaded(fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B)
    # Set ``C`` to ``A\times B``. The number of rows of ``B`` must match the
    # number of columns of ``A``.

    void fmpz_mod_mat_sqr(fmpz_mod_mat_t B, const fmpz_mod_mat_t A)
    # Set ``B`` to ``A^2``. The matrix ``A`` must be square.

    void fmpz_mod_mat_mul_fmpz_vec(fmpz * c, const fmpz_mod_mat_t A, const fmpz * b, slong blen)
    void fmpz_mod_mat_mul_fmpz_vec_ptr(fmpz * const * c, const fmpz_mod_mat_t A, const fmpz * const * b, slong blen)
    # Compute a matrix-vector product of ``A`` and ``(b, blen)`` and store the result in ``c``.
    # The vector ``(b, blen)`` is either truncated or zero-extended to the number of columns of ``A``.
    # The number entries written to ``c`` is always equal to the number of rows of ``A``.

    void fmpz_mod_mat_fmpz_vec_mul(fmpz * c, const fmpz * a, slong alen, const fmpz_mod_mat_t B)
    void fmpz_mod_mat_fmpz_vec_mul_ptr(fmpz * const * c, const fmpz * const * a, slong alen, const fmpz_mod_mat_t B)
    # Compute a vector-matrix product of ``(a, alen)`` and ``B`` and and store the result in ``c``.
    # The vector ``(a, alen)`` is either truncated or zero-extended to the number of rows of ``B``.
    # The number entries written to ``c`` is always equal to the number of columns of ``B``.

    void fmpz_mod_mat_trace(fmpz_t trace, const fmpz_mod_mat_t mat)
    # Set ``trace`` to the trace of the matrix ``mat``.

    slong fmpz_mod_mat_rref(slong * perm, fmpz_mod_mat_t mat)
    # Uses Gauss-Jordan elimination to set ``mat`` to its reduced row echelon
    # form and returns the rank of ``mat``.
    # If ``perm`` is non-``NULL``, the permutation of
    # rows in the matrix will also be applied to ``perm``.
    # The modulus is assumed to be prime.

    void fmpz_mod_mat_strong_echelon_form(fmpz_mod_mat_t mat)
    # Transforms `mat` into the strong echelon form of `mat`. The Howell form and the
    # strong echelon form are equal up to permutation of the rows, see
    # [FieHof2014]_ for a definition of the strong echelon form and the
    # algorithm used here.
    # `mat` must have at least as many rows as columns.

    slong fmpz_mod_mat_howell_form(fmpz_mod_mat_t mat)
    # Transforms `mat` into the Howell form of `mat`.  For a definition of the
    # Howell form see [StoMul1998]_. The Howell form is computed by first
    # putting `mat` into strong echelon form and then ordering the rows.
    # `mat` must have at least as many rows as columns.

    int fmpz_mod_mat_inv(fmpz_mod_mat_t B, fmpz_mod_mat_t A)
    # Sets `B = A^{-1}` and returns `1` if `A` is invertible. If `A` is singular,
    # returns `0` and sets the elements of `B` to undefined values.
    # `A` and `B` must be square matrices with the same dimensions.
    # The modulus is assumed to be prime.

    slong fmpz_mod_mat_lu(slong * P, fmpz_mod_mat_t A, int rank_check)
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
    # The modulus is assumed to be prime.

    void fmpz_mod_mat_solve_tril(fmpz_mod_mat_t X, const fmpz_mod_mat_t L, const fmpz_mod_mat_t B, int unit)
    # Sets `X = L^{-1} B` where `L` is a full rank lower triangular
    # square matrix. If ``unit`` = 1, `L` is assumed to have ones on
    # its main diagonal, and the main diagonal will not be read.  `X`
    # and `B` are allowed to be the same matrix, but no other aliasing
    # is allowed. Automatically chooses between the classical and
    # recursive algorithms.
    # The modulus is assumed to be prime.

    void fmpz_mod_mat_solve_triu(fmpz_mod_mat_t X, const fmpz_mod_mat_t U, const fmpz_mod_mat_t B, int unit)
    # Sets `X = U^{-1} B` where `U` is a full rank upper triangular
    # square matrix. If ``unit`` = 1, `U` is assumed to have ones on
    # its main diagonal, and the main diagonal will not be read.  `X`
    # and `B` are allowed to be the same matrix, but no other aliasing
    # is allowed. Automatically chooses between the classical and
    # recursive algorithms.
    # The modulus is assumed to be prime.

    int fmpz_mod_mat_solve(fmpz_mod_mat_t X, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B)
    # Solves the matrix-matrix equation `AX = B`.
    # Returns `1` if `A` has full rank; otherwise returns `0` and sets the
    # elements of `X` to undefined values.
    # The matrix `A` must be square.
    # The modulus is assumed to be prime.

    int fmpz_mod_mat_can_solve(fmpz_mod_mat_t X, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B)
    # Solves the matrix-matrix equation `AX = B` over `Fp`.
    # Returns `1` if a solution exists; otherwise returns `0` and sets the
    # elements of `X` to zero. If more than one solution exists, one of the
    # valid solutions is given.
    # There are no restrictions on the shape of `A` and it may be singular.
    # The modulus is assumed to be prime.

    void fmpz_mod_mat_similarity(fmpz_mod_mat_t M, slong r, fmpz_t d)
    # Applies a similarity transform to the `n\times n` matrix `M` in-place.
    # If `P` is the `n\times n` identity matrix the zero entries of whose row
    # `r` (`0`-indexed) have been replaced by `d`, this transform is equivalent
    # to `M = P^{-1}MP`.
    # Similarity transforms preserve the determinant, characteristic polynomial
    # and minimal polynomial.
    # The value `d` is required to be reduced modulo the modulus of the entries
    # in the matrix.
    # The modulus is assumed to be prime.

    void fmpz_mod_mat_charpoly(fmpz_mod_poly_t p, const fmpz_mod_mat_t M, const fmpz_mod_ctx_t ctx)
    # Compute the characteristic polynomial `p` of the matrix `M`. The matrix
    # is required to be square, otherwise an exception is raised.

    void fmpz_mod_mat_minpoly(fmpz_mod_poly_t p, const fmpz_mod_mat_t M, const fmpz_mod_ctx_t ctx)
    # Compute the minimal polynomial `p` of the matrix `M`. The matrix
    # is required to be square, otherwise an exception is raised.
    # The modulus is assumed to be prime.
