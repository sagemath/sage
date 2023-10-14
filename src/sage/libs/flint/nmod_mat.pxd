# distutils: libraries = flint
# distutils: depends = flint/nmod_mat.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void nmod_mat_init(nmod_mat_t mat, slong rows, slong cols, mp_limb_t n)
    # Initialises ``mat`` to a ``rows``-by-``cols`` matrix with
    # coefficients modulo `n`, where `n` can be any nonzero integer that
    # fits in a limb. All elements are set to zero.

    void nmod_mat_init_set(nmod_mat_t mat, const nmod_mat_t src)
    # Initialises ``mat`` and sets its dimensions, modulus and elements
    # to those of ``src``.

    void nmod_mat_clear(nmod_mat_t mat)
    # Clears the matrix and releases any memory it used. The matrix
    # cannot be used again until it is initialised. This function must be
    # called exactly once when finished using an ``nmod_mat_t`` object.

    void nmod_mat_set(nmod_mat_t mat, const nmod_mat_t src)
    # Sets ``mat`` to a copy of ``src``. It is assumed
    # that ``mat`` and ``src`` have identical dimensions.

    void nmod_mat_swap(nmod_mat_t mat1, nmod_mat_t mat2)
    # Exchanges ``mat1`` and ``mat2``.

    void nmod_mat_swap_entrywise(nmod_mat_t mat1, nmod_mat_t mat2)
    # Swaps two matrices by swapping the individual entries rather than swapping
    # the contents of the structs.

    MACRO nmod_mat_entry(nmod_mat_t mat, slong i, slong j)
    # Directly accesses the entry in ``mat`` in row `i` and column `j`,
    # indexed from zero. No bounds checking is performed. This macro can be
    # used both for reading and writing coefficients.

    mp_limb_t nmod_mat_get_entry(const nmod_mat_t mat, slong i, slong j)
    # Get the entry at row `i` and column `j` of the matrix ``mat``.

    mp_limb_t * nmod_mat_entry_ptr(const nmod_mat_t mat, slong i, slong j)
    # Return a pointer to the entry at row `i` and column `j` of the matrix
    # ``mat``.

    void nmod_mat_set_entry(nmod_mat_t mat, slong i, slong j, mp_limb_t x)
    # Set the entry at row `i` and column `j` of the matrix ``mat`` to
    # ``x``.

    slong nmod_mat_nrows(const nmod_mat_t mat)
    # Returns the number of rows in ``mat``.

    slong nmod_mat_ncols(const nmod_mat_t mat)
    # Returns the number of columns in ``mat``.

    void nmod_mat_zero(nmod_mat_t mat)
    # Sets all entries of the matrix ``mat`` to zero.

    bint nmod_mat_is_zero(const nmod_mat_t mat)
    # Returns `1` if all entries of the matrix ``mat`` are zero.

    void nmod_mat_window_init(nmod_mat_t window, const nmod_mat_t mat, slong r1, slong c1, slong r2, slong c2)
    # Initializes the matrix ``window`` to be an ``r2 - r1`` by
    # ``c2 - c1`` submatrix of ``mat`` whose ``(0,0)`` entry
    # is the ``(r1, c1)`` entry of ``mat``. The memory for the
    # elements of ``window`` is shared with ``mat``.

    void nmod_mat_window_clear(nmod_mat_t window)
    # Clears the matrix ``window`` and releases any memory that it
    # uses. Note that the memory to the underlying matrix that
    # ``window`` points to is not freed.

    void nmod_mat_concat_vertical(nmod_mat_t res, const nmod_mat_t mat1, const nmod_mat_t mat2)
    # Sets ``res`` to vertical concatenation of (`mat1`, ``mat2``) in that order. Matrix dimensions : ``mat1`` : `m \times n`, ``mat2`` : `k \times n`, ``res`` : `(m + k) \times n`.

    void nmod_mat_concat_horizontal(nmod_mat_t res, const nmod_mat_t mat1, const nmod_mat_t mat2)
    # Sets ``res`` to horizontal concatenation of (``mat1``, ``mat2``) in that order. Matrix dimensions : ``mat1`` : `m \times n`, ``mat2`` : `m \times k`, ``res``  : `m \times (n + k)`.

    void nmod_mat_print_pretty(const nmod_mat_t mat)
    # Pretty-prints ``mat`` to ``stdout``. A header is printed followed
    # by the rows enclosed in brackets. Each column is right-aligned to the
    # width of the modulus written in decimal, and the columns are separated by
    # spaces.
    # For example::
    # <2 x 3 integer matrix mod 2903>
    # [   0    0 2607]
    # [ 622    0    0]

    int nmod_mat_fprint_pretty(FILE * file, const nmod_mat_t mat)
    # Same as ``nmod_mat_print_pretty`` but printing to ``file``.

    int nmod_mat_print(const nmod_mat_t mat)
    # Currently, same as ``nmod_mat_print_pretty``.

    int nmod_mat_fprint(FILE * f, const nmod_mat_t mat)
    # Currently, same as ``nmod_mat_fprint_pretty``.

    void nmod_mat_randtest(nmod_mat_t mat, flint_rand_t state)
    # Sets the elements to a random matrix with entries between `0` and `m-1`
    # inclusive, where `m` is the modulus of ``mat``. A sparse matrix is
    # generated with increased probability.

    void nmod_mat_randfull(nmod_mat_t mat, flint_rand_t state)
    # Sets the element to random numbers likely to be close to the modulus
    # of the matrix. This is used to test potential overflow-related bugs.

    int nmod_mat_randpermdiag(nmod_mat_t mat, flint_rand_t state, mp_srcptr diag, slong n)
    # Sets ``mat`` to a random permutation of the diagonal matrix
    # with `n` leading entries given by the vector ``diag``. It is
    # assumed that the main diagonal of ``mat`` has room for at
    # least `n` entries.
    # Returns `0` or `1`, depending on whether the permutation is even
    # or odd respectively.

    void nmod_mat_randrank(nmod_mat_t mat, flint_rand_t state, slong rank)
    # Sets ``mat`` to a random sparse matrix with the given rank,
    # having exactly as many non-zero elements as the rank, with the
    # non-zero elements being uniformly random integers between `0`
    # and `m-1` inclusive, where `m` is the modulus of ``mat``.
    # The matrix can be transformed into a dense matrix with unchanged
    # rank by subsequently calling :func:`nmod_mat_randops`.

    void nmod_mat_randops(nmod_mat_t mat, slong count, flint_rand_t state)
    # Randomises ``mat`` by performing elementary row or column
    # operations. More precisely, at most ``count`` random additions
    # or subtractions of distinct rows and columns will be performed.
    # This leaves the rank (and for square matrices, determinant)
    # unchanged.

    void nmod_mat_randtril(nmod_mat_t mat, flint_rand_t state, int unit)
    # Sets ``mat`` to a random lower triangular matrix. If ``unit`` is 1,
    # it will have ones on the main diagonal, otherwise it will have random
    # nonzero entries on the main diagonal.

    void nmod_mat_randtriu(nmod_mat_t mat, flint_rand_t state, int unit)
    # Sets ``mat`` to a random upper triangular matrix. If ``unit`` is 1,
    # it will have ones on the main diagonal, otherwise it will have random
    # nonzero entries on the main diagonal.

    bint nmod_mat_equal(const nmod_mat_t mat1, const nmod_mat_t mat2)
    # Returns nonzero if ``mat1`` and ``mat2`` have the same dimensions and elements,
    # and zero otherwise. The moduli are ignored.

    bint nmod_mat_is_zero_row(const nmod_mat_t mat, slong i)
    # Returns a non-zero value if row `i` of ``mat`` is zero.

    void nmod_mat_transpose(nmod_mat_t B, const nmod_mat_t A)
    # Sets `B` to the transpose of `A`. Dimensions must be compatible.
    # `B` and `A` may be the same object if and only if the matrix is square.

    void nmod_mat_swap_rows(nmod_mat_t mat, slong * perm, slong r, slong s)
    # Swaps rows ``r`` and ``s`` of ``mat``.  If ``perm`` is non-``NULL``, the
    # permutation of the rows will also be applied to ``perm``.

    void nmod_mat_swap_cols(nmod_mat_t mat, slong * perm, slong r, slong s)
    # Swaps columns ``r`` and ``s`` of ``mat``.  If ``perm`` is non-``NULL``, the
    # permutation of the columns will also be applied to ``perm``.

    void nmod_mat_invert_rows(nmod_mat_t mat, slong * perm)
    # Swaps rows ``i`` and ``r - i`` of ``mat`` for ``0 <= i < r/2``, where
    # ``r`` is the number of rows of ``mat``. If ``perm`` is non-``NULL``, the
    # permutation of the rows will also be applied to ``perm``.

    void nmod_mat_invert_cols(nmod_mat_t mat, slong * perm)
    # Swaps columns ``i`` and ``c - i`` of ``mat`` for ``0 <= i < c/2``, where
    # ``c`` is the number of columns of ``mat``. If ``perm`` is non-``NULL``, the
    # permutation of the columns will also be applied to ``perm``.

    void nmod_mat_permute_rows(nmod_mat_t mat, const slong * perm_act, slong * perm_store)
    # Permutes rows of the matrix ``mat`` according to permutation ``perm_act``
    # and, if ``perm_store`` is not ``NULL``, apply the same permutation to it.

    void nmod_mat_add(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
    # Computes `C = A + B`. Dimensions must be identical.

    void nmod_mat_sub(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
    # Computes `C = A - B`. Dimensions must be identical.

    void nmod_mat_neg(nmod_mat_t A, const nmod_mat_t B)
    # Sets `B = -A`. Dimensions must be identical.

    void nmod_mat_scalar_mul(nmod_mat_t B, const nmod_mat_t A, mp_limb_t c)
    # Sets `B = cA`, where the scalar `c` is assumed to be reduced
    # modulo the modulus. Dimensions of `A` and `B` must be identical.

    void nmod_mat_scalar_addmul_ui(nmod_mat_t dest, const nmod_mat_t X, const nmod_mat_t Y, const mp_limb_t b)
    # Sets `dest = X + bY`, where the scalar `b` is assumed to be reduced
    # modulo the modulus. Dimensions of dest, X and Y must be identical.
    # dest can be aliased with X or Y.

    void nmod_mat_scalar_mul_fmpz(nmod_mat_t res, const nmod_mat_t M, const fmpz_t c)
    # Sets `B = cA`, where the scalar `c` is of type ``fmpz_t``. Dimensions of `A`
    # and `B` must be identical.

    void nmod_mat_mul(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
    # Sets `C = AB`. Dimensions must be compatible for matrix multiplication.
    # Aliasing is allowed. This function automatically chooses between classical
    # and Strassen multiplication.

    void _nmod_mat_mul_classical_op(nmod_mat_t D, const nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, int op)

    void nmod_mat_mul_classical(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
    # Sets `C = AB`. Dimensions must be compatible for matrix multiplication.
    # `C` is not allowed to be aliased with `A` or `B`. Uses classical
    # matrix multiplication, creating a temporary transposed copy of `B`
    # to improve memory locality if the matrices are large enough,
    # and packing several entries of `B` into each word if the modulus
    # is very small.

    void _nmod_mat_mul_classical_threaded_pool_op(nmod_mat_t D, const nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, int op, thread_pool_handle * threads, slong num_threads)
    # Multithreaded version of ``_nmod_mat_mul_classical``.

    void _nmod_mat_mul_classical_threaded_op(nmod_mat_t D, const nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, int op)
    # Multithreaded version of ``_nmod_mat_mul_classical``.

    void nmod_mat_mul_classical_threaded(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
    # Multithreaded version of ``nmod_mat_mul_classical``.

    void nmod_mat_mul_strassen(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
    # Sets `C = AB`. Dimensions must be compatible for matrix multiplication.
    # `C` is not allowed to be aliased with `A` or `B`. Uses Strassen
    # multiplication (the Strassen-Winograd variant).

    int nmod_mat_mul_blas(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
    # Tries to set `C = AB` using BLAS and returns `1` for success and `0` for failure. Dimensions must be compatible for matrix multiplication.

    void nmod_mat_addmul(nmod_mat_t D, const nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
    # Sets `D = C + AB`. `C` and `D` may be aliased with each other but
    # not with `A` or `B`. Automatically selects between classical
    # and Strassen multiplication.

    void nmod_mat_submul(nmod_mat_t D, const nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
    # Sets `D = C + AB`. `C` and `D` may be aliased with each other but
    # not with `A` or `B`.

    void nmod_mat_mul_nmod_vec(mp_limb_t * c, const nmod_mat_t A, const mp_limb_t * b, slong blen)
    void nmod_mat_mul_nmod_vec_ptr(mp_limb_t * const * c, const nmod_mat_t A, const mp_limb_t * const * b, slong blen)
    # Compute a matrix-vector product of ``A`` and ``(b, blen)`` and store the result in ``c``.
    # The vector ``(b, blen)`` is either truncated or zero-extended to the number of columns of ``A``.
    # The number entries written to ``c`` is always equal to the number of rows of ``A``.

    void nmod_mat_nmod_vec_mul(mp_limb_t * c, const mp_limb_t * a, slong alen, const nmod_mat_t B)
    void nmod_mat_nmod_vec_mul_ptr(mp_limb_t * const * c, const mp_limb_t * const * a, slong alen, const nmod_mat_t B)
    # Compute a vector-matrix product of ``(a, alen)`` and ``B`` and and store the result in ``c``.
    # The vector ``(a, alen)`` is either truncated or zero-extended to the number of rows of ``B``.
    # The number entries written to ``c`` is always equal to the number of columns of ``B``.

    void _nmod_mat_pow(nmod_mat_t dest, const nmod_mat_t mat, ulong pow)

    void nmod_mat_pow(nmod_mat_t dest, const nmod_mat_t mat, ulong pow)
    # Sets `dest = mat^{pow}`. ``dest`` and ``mat`` may be aliased. Implements

    mp_limb_t nmod_mat_trace(const nmod_mat_t mat)
    # Computes the trace of the matrix, i.e. the sum of the entries on
    # the main diagonal. The matrix is required to be square.

    mp_limb_t nmod_mat_det_howell(const nmod_mat_t A)
    # Returns the determinant of `A`.

    mp_limb_t nmod_mat_det(const nmod_mat_t A)
    # Returns the determinant of `A`.

    slong nmod_mat_rank(const nmod_mat_t A)
    # Returns the rank of `A`. The modulus of `A` must be a prime number.

    int nmod_mat_inv(nmod_mat_t B, const nmod_mat_t A)
    # Sets `B = A^{-1}` and returns `1` if `A` is invertible.
    # If `A` is singular, returns `0` and sets the elements of
    # `B` to undefined values.
    # `A` and `B` must be square matrices with the same dimensions
    # and modulus. The modulus must be prime.

    void nmod_mat_solve_tril(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit)
    # Sets `X = L^{-1} B` where `L` is a full rank lower triangular square
    # matrix. If ``unit`` = 1, `L` is assumed to have ones on its
    # main diagonal, and the main diagonal will not be read.
    # `X` and `B` are allowed to be the same matrix, but no other
    # aliasing is allowed. Automatically chooses between the classical and
    # recursive algorithms.

    void nmod_mat_solve_tril_classical(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit)
    # Sets `X = L^{-1} B` where `L` is a full rank lower triangular square
    # matrix. If ``unit`` = 1, `L` is assumed to have ones on its
    # main diagonal, and the main diagonal will not be read.
    # `X` and `B` are allowed to be the same matrix, but no other
    # aliasing is allowed. Uses forward substitution.

    void nmod_mat_solve_tril_recursive(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit)
    # Sets `X = L^{-1} B` where `L` is a full rank lower triangular square
    # matrix. If ``unit`` = 1, `L` is assumed to have ones on its
    # main diagonal, and the main diagonal will not be read.
    # `X` and `B` are allowed to be the same matrix, but no other
    # aliasing is allowed.
    # Uses the block inversion formula
    # .. math ::
    # \begin{pmatrix} A & 0 \\ C & D \end{pmatrix}^{-1}
    # \begin{pmatrix} X \\ Y \end{pmatrix} =
    # \begin{pmatrix} A^{-1} X \\ D^{-1} ( Y - C A^{-1} X ) \end{pmatrix}
    # to reduce the problem to matrix multiplication and triangular solving
    # of smaller systems.

    void nmod_mat_solve_triu(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit)
    # Sets `X = U^{-1} B` where `U` is a full rank upper triangular square
    # matrix. If ``unit`` = 1, `U` is assumed to have ones on its
    # main diagonal, and the main diagonal will not be read.
    # `X` and `B` are allowed to be the same matrix, but no other
    # aliasing is allowed. Automatically chooses between the classical and
    # recursive algorithms.

    void nmod_mat_solve_triu_classical(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit)
    # Sets `X = U^{-1} B` where `U` is a full rank upper triangular square
    # matrix. If ``unit`` = 1, `U` is assumed to have ones on its
    # main diagonal, and the main diagonal will not be read.
    # `X` and `B` are allowed to be the same matrix, but no other
    # aliasing is allowed. Uses forward substitution.

    void nmod_mat_solve_triu_recursive(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit)
    # Sets `X = U^{-1} B` where `U` is a full rank upper triangular square
    # matrix. If ``unit`` = 1, `U` is assumed to have ones on its
    # main diagonal, and the main diagonal will not be read.
    # `X` and `B` are allowed to be the same matrix, but no other
    # aliasing is allowed.
    # Uses the block inversion formula
    # .. math ::
    # \begin{pmatrix} A & B \\ 0 & D \end{pmatrix}^{-1}
    # \begin{pmatrix} X \\ Y \end{pmatrix} =
    # \begin{pmatrix} A^{-1} (X - B D^{-1} Y) \\ D^{-1} Y \end{pmatrix}
    # to reduce the problem to matrix multiplication and triangular solving
    # of smaller systems.

    int nmod_mat_solve(nmod_mat_t X, const nmod_mat_t A, const nmod_mat_t B)
    # Solves the matrix-matrix equation `AX = B` over `\mathbb{Z} / p \mathbb{Z}` where `p`
    # is the modulus of `X` which must be a prime number. `X`, `A`, and `B`
    # should have the same moduli.
    # Returns `1` if `A` has full rank; otherwise returns `0` and sets the
    # elements of `X` to undefined values.
    # The matrix `A` must be square.

    int nmod_mat_can_solve_inner(slong * rank, slong * perm, slong * pivots, nmod_mat_t X, const nmod_mat_t A, const nmod_mat_t B)
    # As for :func:`nmod_mat_can_solve` except that if `rank` is not `NULL` the
    # value it points to will be set to the rank of `A`. If `perm` is not `NULL`
    # then it must be a valid initialised permutation whose length is the number
    # of rows of `A`. After the function call it will be set to the row
    # permutation given by LU decomposition of `A`. If `pivots` is not `NULL`
    # then it must an initialised vector. Only the first `*rank` of these will be
    # set by the function call. They are set to the columns of the pivots chosen
    # by the LU decomposition of `A`.

    int nmod_mat_can_solve(nmod_mat_t X, const nmod_mat_t A, const nmod_mat_t B)
    # Solves the matrix-matrix equation `AX = B` over `\mathbb{Z} / p \mathbb{Z}` where `p`
    # is the modulus of `X` which must be a prime number. `X`, `A`, and `B`
    # should have the same moduli.
    # Returns `1` if a solution exists; otherwise returns `0` and sets the
    # elements of `X` to zero. If more than one solution exists, one of the
    # valid solutions is given.
    # There are no restrictions on the shape of `A` and it may be singular.

    int nmod_mat_solve_vec(mp_ptr x, const nmod_mat_t A, mp_srcptr b)
    # Solves the matrix-vector equation `Ax = b` over `\mathbb{Z} / p \mathbb{Z}` where `p`
    # is the modulus of `A` which must be a prime number.
    # Returns `1` if `A` has full rank; otherwise returns `0` and sets the
    # elements of `x` to undefined values.

    slong nmod_mat_lu(slong * P, nmod_mat_t A, int rank_check)
    slong nmod_mat_lu_classical(slong * P, nmod_mat_t A, int rank_check)
    slong nmod_mat_lu_classical_delayed(slong * P, nmod_mat_t A, int rank_check)
    slong nmod_mat_lu_recursive(slong * P, nmod_mat_t A, int rank_check)
    # Computes a generalised LU decomposition `LU = PA` of a given
    # matrix `A`, returning the rank of `A`.
    # If `A` is a nonsingular square matrix, it will be overwritten with
    # a unit diagonal lower triangular matrix `L` and an upper triangular
    # matrix `U` (the diagonal of `L` will not be stored explicitly).
    # If `A` is an arbitrary matrix of rank `r`, `U` will be in row echelon
    # form having `r` nonzero rows, and `L` will be lower triangular
    # but truncated to `r` columns, having implicit ones on the `r` first
    # entries of the main diagonal. All other entries will be zero.
    # If a nonzero value for ``rank_check`` is passed, the
    # function will abandon the output matrix in an undefined state and
    # return 0 if `A` is detected to be rank-deficient.
    # The *classical* version uses direct Gaussian elimination.
    # The *classical_delayed* version also uses Gaussian elimination,
    # but performs delayed modular reductions.
    # The *recursive* version uses block recursive decomposition.
    # The default function chooses an algorithm automatically.

    slong nmod_mat_rref(nmod_mat_t A)
    # Puts `A` in reduced row echelon form and returns the rank of `A`.
    # The rref is computed by first obtaining an unreduced row echelon
    # form via LU decomposition and then solving an additional
    # triangular system.

    slong nmod_mat_reduce_row(nmod_mat_t A, slong * P, slong * L, slong n)
    # Reduce row n of the matrix `A`, assuming the prior rows are in Gauss
    # form. However those rows may not be in order. The entry `i` of the array
    # `P` is the row of `A` which has a pivot in the `i`-th column. If no such
    # row exists, the entry of `P` will be `-1`. The function returns the column
    # in which the `n`-th row has a pivot after reduction. This will always be
    # chosen to be the first available column for a pivot from the left. This
    # information is also updated in `P`. Entry `i` of the array `L` contains the
    # number of possibly nonzero columns of `A` row `i`. This speeds up reduction
    # in the case that `A` is chambered on the right. Otherwise the entries of `L`
    # can all be set to the number of columns of `A`. We require the entries of
    # `L` to be monotonic increasing.

    slong nmod_mat_nullspace(nmod_mat_t X, const nmod_mat_t A)
    # Computes the nullspace of `A` and returns the nullity.
    # More precisely, this function sets `X` to a maximum rank matrix
    # such that `AX = 0` and returns the rank of `X`. The columns of
    # `X` will form a basis for the nullspace of `A`.
    # `X` must have sufficient space to store all basis vectors
    # in the nullspace.
    # This function computes the reduced row echelon form and then reads
    # off the basis vectors.

    void nmod_mat_similarity(nmod_mat_t M, slong r, ulong d)
    # Applies a similarity transform to the `n\times n` matrix `M` in-place.
    # If `P` is the `n\times n` identity matrix the zero entries of whose row
    # `r` (`0`-indexed) have been replaced by `d`, this transform is equivalent
    # to `M = P^{-1}MP`.
    # Similarity transforms preserve the determinant, characteristic polynomial
    # and minimal polynomial.
    # The value `d` is required to be reduced modulo the modulus of the entries
    # in the matrix.

    void nmod_mat_charpoly_berkowitz(nmod_poly_t p, const nmod_mat_t M)
    void nmod_mat_charpoly_danilevsky(nmod_poly_t p, const nmod_mat_t M)
    void nmod_mat_charpoly(nmod_poly_t p, const nmod_mat_t M)
    # Compute the characteristic polynomial `p` of the matrix `M`. The matrix
    # is required to be square, otherwise an exception is raised.
    # The *danilevsky* algorithm assumes that the modulus is prime.

    void nmod_mat_minpoly(nmod_poly_t p, const nmod_mat_t M)
    # Compute the minimal polynomial `p` of the matrix `M`. The matrix
    # is required to be square, otherwise an exception is raised.

    void nmod_mat_strong_echelon_form(nmod_mat_t A)
    # Puts `A` into strong echelon form. The Howell form and the strong echelon
    # form are equal up to permutation of the rows, see [FieHof2014]_ for a
    # definition of the strong echelon form and the algorithm used here.
    # Note that [FieHof2014]_ defines strong echelon form as a lower left normal form,
    # while the implemented version returns an upper right normal form,
    # agreeing with the definition of Howell form in [StoMul1998]_.
    # `A` must have at least as many rows as columns.

    slong nmod_mat_howell_form(nmod_mat_t A)
    # Puts `A` into Howell form and returns the number of non-zero rows.
    # For a definition of the Howell form see [StoMul1998]_. The Howell form
    # is computed by first putting `A` into strong echelon form and then ordering
    # the rows.
    # `A` must have at least as many rows as columns.
