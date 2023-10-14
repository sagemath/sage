# distutils: libraries = flint
# distutils: depends = flint/fmpz_mat.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void fmpz_mat_init(fmpz_mat_t mat, slong rows, slong cols)
    # Initialises a matrix with the given number of rows and columns for use.

    void fmpz_mat_clear(fmpz_mat_t mat)
    # Clears the given matrix.

    void fmpz_mat_set(fmpz_mat_t mat1, const fmpz_mat_t mat2)
    # Sets ``mat1`` to a copy of ``mat2``. The dimensions of
    # ``mat1`` and ``mat2`` must be the same.

    void fmpz_mat_init_set(fmpz_mat_t mat, const fmpz_mat_t src)
    # Initialises the matrix ``mat`` to the same size as ``src`` and
    # sets it to a copy of ``src``.

    slong fmpz_mat_nrows(const fmpz_mat_t mat)
    slong fmpz_mat_ncols(const fmpz_mat_t mat)
    # Returns respectively the number of rows and columns of the matrix.

    void fmpz_mat_swap(fmpz_mat_t mat1, fmpz_mat_t mat2)
    # Swaps two matrices. The dimensions of ``mat1`` and ``mat2``
    # are allowed to be different.

    void fmpz_mat_swap_entrywise(fmpz_mat_t mat1, fmpz_mat_t mat2)
    # Swaps two matrices by swapping the individual entries rather than swapping
    # the contents of the structs.

    fmpz * fmpz_mat_entry(const fmpz_mat_t mat, slong i, slong j)
    # Returns a reference to the entry of ``mat`` at row `i` and column `j`.
    # This reference can be passed as an input or output variable to any
    # function in the ``fmpz`` module for direct manipulation.
    # Both `i` and `j` must not exceed the dimensions of the matrix.
    # This function is implemented as a macro.

    void fmpz_mat_zero(fmpz_mat_t mat)
    # Sets all entries of ``mat`` to 0.

    void fmpz_mat_one(fmpz_mat_t mat)
    # Sets ``mat`` to the unit matrix, having ones on the main diagonal
    # and zeroes elsewhere. If ``mat`` is nonsquare, it is set to the
    # truncation of a unit matrix.

    void fmpz_mat_swap_rows(fmpz_mat_t mat, slong * perm, slong r, slong s)
    # Swaps rows ``r`` and ``s`` of ``mat``.  If ``perm`` is non-``NULL``, the
    # permutation of the rows will also be applied to ``perm``.

    void fmpz_mat_swap_cols(fmpz_mat_t mat, slong * perm, slong r, slong s)
    # Swaps columns ``r`` and ``s`` of ``mat``.  If ``perm`` is non-``NULL``, the
    # permutation of the columns will also be applied to ``perm``.

    void fmpz_mat_invert_rows(fmpz_mat_t mat, slong * perm)
    # Swaps rows ``i`` and ``r - i`` of ``mat`` for ``0 <= i < r/2``, where
    # ``r`` is the number of rows of ``mat``. If ``perm`` is non-``NULL``, the
    # permutation of the rows will also be applied to ``perm``.

    void fmpz_mat_invert_cols(fmpz_mat_t mat, slong * perm)
    # Swaps columns ``i`` and ``c - i`` of ``mat`` for ``0 <= i < c/2``, where
    # ``c`` is the number of columns of ``mat``. If ``perm`` is non-``NULL``, the
    # permutation of the columns will also be applied to ``perm``.

    void fmpz_mat_window_init(fmpz_mat_t window, const fmpz_mat_t mat, slong r1, slong c1, slong r2, slong c2)
    # Initializes the matrix ``window`` to be an ``r2 - r1`` by
    # ``c2 - c1`` submatrix of ``mat`` whose ``(0,0)`` entry
    # is the ``(r1, c1)`` entry of ``mat``. The memory for the
    # elements of ``window`` is shared with ``mat``.

    void fmpz_mat_window_clear(fmpz_mat_t window)
    # Clears the matrix ``window`` and releases any memory that it
    # uses. Note that the memory to the underlying matrix that
    # ``window`` points to is not freed.

    void fmpz_mat_randbits(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits)
    # Sets the entries of ``mat`` to random signed integers whose absolute
    # values have the given number of binary bits.

    void fmpz_mat_randtest(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits)
    # Sets the entries of ``mat`` to random signed integers whose
    # absolute values have a random number of bits up to the given number
    # of bits inclusive.

    void fmpz_mat_randintrel(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits)
    # Sets ``mat`` to be a random *integer relations* matrix, with
    # signed entries up to the given number of bits.
    # The number of columns of ``mat`` must be equal to one more than
    # the number of rows. The format of the matrix is a set of random integers
    # in the left hand column and an identity matrix in the remaining square
    # submatrix.

    void fmpz_mat_randsimdioph(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits, flint_bitcnt_t bits2)
    # Sets ``mat`` to a random *simultaneous diophantine* matrix.
    # The matrix must be square. The top left entry is set to ``2^bits2``.
    # The remainder of that row is then set to signed random integers of the
    # given number of binary bits. The remainder of the first column is zero.
    # Running down the rest of the diagonal are the values ``2^bits`` with
    # all remaining entries zero.

    void fmpz_mat_randntrulike(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits, ulong q)
    # Sets a square matrix ``mat`` of even dimension to a random
    # *NTRU like* matrix.
    # The matrix is broken into four square submatrices. The top left submatrix
    # is set to the identity. The bottom left submatrix is set to the zero
    # matrix. The bottom right submatrix is set to `q` times the identity matrix.
    # Finally the top right submatrix has the following format. A random vector
    # `h` of length `r/2` is created, with random signed entries of the given
    # number of bits. Then entry `(i, j)` of the submatrix is set to
    # `h[i + j \bmod{r/2}]`.

    void fmpz_mat_randntrulike2(fmpz_mat_t mat, flint_rand_t state, flint_bitcnt_t bits, ulong q)
    # Sets a square matrix ``mat`` of even dimension to a random
    # *NTRU like* matrix.
    # The matrix is broken into four square submatrices. The top left submatrix
    # is set to `q` times the identity matrix. The top right submatrix is set to
    # the zero matrix. The bottom right submatrix is set to the identity matrix.
    # Finally the bottom left submatrix has the following format. A random vector
    # `h` of length `r/2` is created, with random signed entries of the given
    # number of bits. Then entry `(i, j)` of the submatrix is set to
    # `h[i + j \bmod{r/2}]`.

    void fmpz_mat_randajtai(fmpz_mat_t mat, flint_rand_t state, double alpha)
    # Sets a square matrix ``mat`` to a random *ajtai* matrix.
    # The diagonal entries `(i, i)` are set to a random entry in the range
    # `[1, 2^{b-1}]` inclusive where `b = \lfloor(2 r - i)^\alpha\rfloor` for some
    # double parameter `\alpha`. The entries below the diagonal in column `i`
    # are set to a random entry in the range `(-2^b + 1, 2^b - 1)` whilst the
    # entries to the right of the diagonal in row `i` are set to zero.

    int fmpz_mat_randpermdiag(fmpz_mat_t mat, flint_rand_t state, const fmpz * diag, slong n)
    # Sets ``mat`` to a random permutation of the rows and columns of a
    # given diagonal matrix. The diagonal matrix is specified in the form of
    # an array of the `n` initial entries on the main diagonal.
    # The return value is `0` or `1` depending on whether the permutation is
    # even or odd.

    void fmpz_mat_randrank(fmpz_mat_t mat, flint_rand_t state, slong rank, flint_bitcnt_t bits)
    # Sets ``mat`` to a random sparse matrix with the given rank,
    # having exactly as many non-zero elements as the rank, with the
    # nonzero elements being random integers of the given bit size.
    # The matrix can be transformed into a dense matrix with unchanged
    # rank by subsequently calling :func:`fmpz_mat_randops`.

    void fmpz_mat_randdet(fmpz_mat_t mat, flint_rand_t state, const fmpz_t det)
    # Sets ``mat`` to a random sparse matrix with minimal number of
    # nonzero entries such that its determinant has the given value.
    # Note that the matrix will be zero if ``det`` is zero.
    # In order to generate a non-zero singular matrix, the function
    # :func:`fmpz_mat_randrank` can be used.
    # The matrix can be transformed into a dense matrix with unchanged
    # determinant by subsequently calling :func:`fmpz_mat_randops`.

    void fmpz_mat_randops(fmpz_mat_t mat, flint_rand_t state, slong count)
    # Randomises ``mat`` by performing elementary row or column operations.
    # More precisely, at most ``count`` random additions or subtractions of
    # distinct rows and columns will be performed. This leaves the rank
    # (and for square matrices, the determinant) unchanged.

    int fmpz_mat_fprint(FILE * file, const fmpz_mat_t mat)
    # Prints the given matrix to the stream ``file``.  The format is
    # the number of rows, a space, the number of columns, two spaces, then
    # a space separated list of coefficients, one row after the other.
    # In case of success, returns a positive value;  otherwise, returns
    # a non-positive value.

    int fmpz_mat_fprint_pretty(FILE * file, const fmpz_mat_t mat)
    # Prints the given matrix to the stream ``file``.  The format is an
    # opening square bracket, then on each line a row of the matrix, followed
    # by a closing square bracket. Each row is written as an opening square
    # bracket followed by a space separated list of coefficients followed
    # by a closing square bracket.
    # In case of success, returns a positive value;  otherwise, returns
    # a non-positive value.

    int fmpz_mat_print(const fmpz_mat_t mat)
    # Prints the given matrix to the stream ``stdout``.  For further
    # details, see :func:`fmpz_mat_fprint`.

    int fmpz_mat_print_pretty(const fmpz_mat_t mat)
    # Prints the given matrix to ``stdout``.  For further details,
    # see :func:`fmpz_mat_fprint_pretty`.

    int fmpz_mat_fread(FILE* file, fmpz_mat_t mat)
    # Reads a matrix from the stream ``file``, storing the result
    # in ``mat``.  The expected format is the number of rows, a
    # space, the number of columns, two spaces, then a space separated
    # list of coefficients, one row after the other.
    # In case of success, returns a positive number.  In case of failure,
    # returns a non-positive value.

    int fmpz_mat_read(fmpz_mat_t mat)
    # Reads a matrix from ``stdin``, storing the result
    # in ``mat``.
    # In case of success, returns a positive number.  In case of failure,
    # returns a non-positive value.

    int fmpz_mat_equal(const fmpz_mat_t mat1, const fmpz_mat_t mat2)
    # Returns a non-zero value if ``mat1`` and ``mat2`` have
    # the same dimensions and entries, and zero otherwise.

    int fmpz_mat_is_zero(const fmpz_mat_t mat)
    # Returns a non-zero value if all entries ``mat`` are zero, and
    # otherwise returns zero.

    int fmpz_mat_is_one(const fmpz_mat_t mat)
    # Returns a non-zero value if ``mat`` is the unit matrix or the truncation
    # of a unit matrix, and otherwise returns zero.

    int fmpz_mat_is_empty(const fmpz_mat_t mat)
    # Returns a non-zero value if the number of rows or the number of
    # columns in ``mat`` is zero, and otherwise returns
    # zero.

    int fmpz_mat_is_square(const fmpz_mat_t mat)
    # Returns a non-zero value if the number of rows is equal to the
    # number of columns in ``mat``, and otherwise returns zero.

    int fmpz_mat_is_zero_row(const fmpz_mat_t mat, slong i)
    # Returns a non-zero value if row `i` of ``mat`` is zero.

    int fmpz_mat_equal_col(fmpz_mat_t M, slong m, slong n)
    # Returns `1` if columns `m` and `n` of the matrix `M` are equal, otherwise
    # returns `0`.

    int fmpz_mat_equal_row(fmpz_mat_t M, slong m, slong n)
    # Returns `1` if rows `m` and `n` of the matrix `M` are equal, otherwise
    # returns `0`.

    void fmpz_mat_transpose(fmpz_mat_t B, const fmpz_mat_t A)
    # Sets `B` to `A^T`, the transpose of `A`. Dimensions must be compatible.
    # `A` and `B` are allowed to be the same object if `A` is a square matrix.

    void fmpz_mat_concat_vertical(fmpz_mat_t res, const fmpz_mat_t mat1, const fmpz_mat_t mat2)
    # Sets ``res`` to vertical concatenation of (``mat1``, ``mat2``)
    # in that order. Matrix dimensions: ``mat1``: `m \times n`,
    # ``mat2``: `k \times n`, ``res``: `(m + k) \times n`.

    void fmpz_mat_concat_horizontal(fmpz_mat_t res, const fmpz_mat_t mat1, const fmpz_mat_t mat2)
    # Sets ``res`` to horizontal concatenation of (``mat1``, ``mat2``)
    # in that order. Matrix dimensions: ``mat1``: `m \times n`,
    # ``mat2``: `m \times k`, ``res``: `m \times (n + k)`.

    void fmpz_mat_get_nmod_mat(nmod_mat_t Amod, const fmpz_mat_t A)
    # Sets the entries of ``Amod`` to the entries of ``A`` reduced
    # by the modulus of ``Amod``.

    void fmpz_mat_set_nmod_mat(fmpz_mat_t A, const nmod_mat_t Amod)
    # Sets the entries of ``Amod`` to the residues in ``Amod``,
    # normalised to the interval `-m/2 <= r < m/2` where `m` is the modulus.

    void fmpz_mat_set_nmod_mat_unsigned(fmpz_mat_t A, const nmod_mat_t Amod)
    # Sets the entries of ``Amod`` to the residues in ``Amod``,
    # normalised to the interval `0 <= r < m` where `m` is the modulus.

    void fmpz_mat_CRT_ui(fmpz_mat_t res, const fmpz_mat_t mat1, const fmpz_t m1, const nmod_mat_t mat2, int sign)
    # Given ``mat1`` with entries modulo ``m`` and ``mat2``
    # with modulus `n`, sets ``res`` to the CRT reconstruction modulo `mn`
    # with entries satisfying `-mn/2 <= c < mn/2` (if sign = 1)
    # or `0 <= c < mn` (if sign = 0).

    void fmpz_mat_multi_mod_ui_precomp(nmod_mat_t * residues, slong nres, const fmpz_mat_t mat, const fmpz_comb_t comb, fmpz_comb_temp_t temp)
    # Sets each of the ``nres`` matrices in ``residues`` to ``mat`` reduced modulo
    # the modulus of the respective matrix, given precomputed ``comb`` and
    # ``comb_temp`` structures.
    # Note: ``fmpz.h`` must be included **before** ``fmpz_mat.h`` in order for
    # this function to be declared.

    void fmpz_mat_multi_mod_ui(nmod_mat_t * residues, slong nres, const fmpz_mat_t mat)
    # Sets each of the ``nres`` matrices in ``residues`` to ``mat``
    # reduced modulo the modulus of the respective matrix.
    # This function is provided for convenience purposes.
    # For reducing or reconstructing multiple integer matrices over the same
    # set of moduli, it is faster to use ``fmpz_mat_multi_mod_precomp``.

    void fmpz_mat_multi_CRT_ui_precomp(fmpz_mat_t mat, nmod_mat_t * const residues, slong nres, const fmpz_comb_t comb, fmpz_comb_temp_t temp, int sign)
    # Reconstructs ``mat`` from its images modulo the ``nres`` matrices in
    # ``residues``, given precomputed ``comb`` and ``comb_temp`` structures.
    # Note: ``fmpz.h`` must be included **before** ``fmpz_mat.h`` in order for
    # this function to be declared.

    void fmpz_mat_multi_CRT_ui(fmpz_mat_t mat, nmod_mat_t * const residues, slong nres, int sign)
    # Reconstructs ``mat`` from its images modulo the ``nres`` matrices
    # in ``residues``.
    # This function is provided for convenience purposes.
    # For reducing or reconstructing multiple integer matrices over the same
    # set of moduli, it is faster to use :func:`fmpz_mat_multi_CRT_ui_precomp`.

    void fmpz_mat_add(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
    # Sets ``C`` to the elementwise sum `A + B`. All inputs must
    # be of the same size. Aliasing is allowed.

    void fmpz_mat_sub(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
    # Sets ``C`` to the elementwise difference `A - B`. All inputs must
    # be of the same size. Aliasing is allowed.

    void fmpz_mat_neg(fmpz_mat_t B, const fmpz_mat_t A)
    # Sets ``B`` to the elementwise negation of ``A``. Both inputs
    # must be of the same size. Aliasing is allowed.

    void fmpz_mat_scalar_mul_si(fmpz_mat_t B, const fmpz_mat_t A, slong c)
    void fmpz_mat_scalar_mul_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c)
    void fmpz_mat_scalar_mul_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c)
    # Set ``B = A*c`` where ``A`` is an ``fmpz_mat_t`` and ``c``
    # is a scalar respectively of type ``slong``, ``ulong``,
    # or ``fmpz_t``. The dimensions of ``A`` and ``B`` must
    # be compatible.

    void fmpz_mat_scalar_addmul_si(fmpz_mat_t B, const fmpz_mat_t A, slong c)
    void fmpz_mat_scalar_addmul_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c)
    void fmpz_mat_scalar_addmul_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c)
    # Set ``B = B + A*c`` where ``A`` is an ``fmpz_mat_t`` and ``c``
    # is a scalar respectively of type ``slong``, ``ulong``,
    # or ``fmpz_t``. The dimensions of ``A`` and ``B`` must
    # be compatible.

    void fmpz_mat_scalar_submul_si(fmpz_mat_t B, const fmpz_mat_t A, slong c)
    void fmpz_mat_scalar_submul_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c)
    void fmpz_mat_scalar_submul_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c)
    # Set ``B = B - A*c`` where ``A`` is an ``fmpz_mat_t`` and ``c``
    # is a scalar respectively of type ``slong``, ``ulong``,
    # or ``fmpz_t``. The dimensions of ``A`` and ``B`` must
    # be compatible.

    void fmpz_mat_scalar_addmul_nmod_mat_ui(fmpz_mat_t B, const nmod_mat_t A, ulong c)
    void fmpz_mat_scalar_addmul_nmod_mat_fmpz(fmpz_mat_t B, const nmod_mat_t A, const fmpz_t c)
    # Set ``B = B + A*c`` where ``A`` is an ``nmod_mat_t`` and ``c``
    # is a scalar respectively of type ``ulong`` or ``fmpz_t``.
    # The dimensions of ``A`` and ``B`` must be compatible.

    void fmpz_mat_scalar_divexact_si(fmpz_mat_t B, const fmpz_mat_t A, slong c)
    void fmpz_mat_scalar_divexact_ui(fmpz_mat_t B, const fmpz_mat_t A, ulong c)
    void fmpz_mat_scalar_divexact_fmpz(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t c)
    # Set ``A = B / c``, where ``B`` is an ``fmpz_mat_t`` and ``c``
    # is a scalar respectively of type ``slong``, ``ulong``,
    # or ``fmpz_t``, which is assumed to divide all elements of
    # ``B`` exactly.

    void fmpz_mat_scalar_mul_2exp(fmpz_mat_t B, const fmpz_mat_t A, ulong exp)
    # Set the matrix ``B`` to the matrix ``A``, of the same dimensions,
    # multiplied by `2^{exp}`.

    void fmpz_mat_scalar_tdiv_q_2exp(fmpz_mat_t B, const fmpz_mat_t A, ulong exp)
    # Set the matrix ``B`` to the matrix ``A``, of the same dimensions,
    # divided by `2^{exp}`, rounding down towards zero.

    void fmpz_mat_scalar_smod(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t P)
    # Set the matrix ``B`` to the matrix ``A``, of the same dimensions,
    # with each entry reduced modulo `P` in the symmetric moduli system. We
    # require `P > 0`.

    void fmpz_mat_mul(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
    # Sets ``C`` to the matrix product `C = A B`. The matrices must have
    # compatible dimensions for matrix multiplication. Aliasing
    # is allowed.
    # This function automatically switches between classical and
    # multimodular multiplication, based on a heuristic comparison of
    # the dimensions and entry sizes.

    void fmpz_mat_mul_classical(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
    # Sets ``C`` to the matrix product `C = A B` computed using
    # classical matrix algorithm.
    # The matrices must have compatible dimensions for matrix multiplication.
    # No aliasing is allowed.

    void fmpz_mat_mul_strassen(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
    # Sets `C = AB`. Dimensions must be compatible for matrix multiplication.
    # `C` is not allowed to be aliased with `A` or `B`. Uses Strassen
    # multiplication (the Strassen-Winograd variant).

    void _fmpz_mat_mul_multi_mod(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B, int sign, flint_bitcnt_t bits)
    void fmpz_mat_mul_multi_mod(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
    # Sets ``C`` to the matrix product `C = AB` computed using a multimodular
    # algorithm. `C` is computed modulo several small prime numbers
    # and reconstructed using the Chinese Remainder Theorem. This generally
    # becomes more efficient than classical multiplication for large matrices.
    # The absolute value of the elements of `C` should be `< 2^{\text{bits}}`,
    # and ``sign`` should be `0` if the entries of `C` are known to be nonnegative
    # and `1` otherwise. The function
    # :func:`fmpz_mat_mul_multi_mod` calculates a rigorous bound automatically.
    # If the default bound is too pessimistic, :func:`_fmpz_mat_mul_multi_mod`
    # can be used with a custom bound.
    # The matrices must have compatible dimensions for matrix multiplication.
    # No aliasing is allowed.

    int fmpz_mat_mul_blas(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
    # Tries to set `C = AB` using BLAS and returns `1` for success and `0` for failure.
    # Dimensions must be compatible for matrix multiplication. No aliasing is allowed.
    # This function currently will fail if the matrices are empty, their dimensions are too large, or their max bits size is over one million bits.

    void fmpz_mat_mul_fft(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
    # Aliasing is allowed.

    void fmpz_mat_sqr(fmpz_mat_t B, const fmpz_mat_t A)
    # Sets ``B`` to the square of the matrix ``A``, which must be
    # a square matrix. Aliasing is allowed.
    # The function calls :func:`fmpz_mat_mul` for dimensions less than 12 and
    # calls :func:`fmpz_mat_sqr_bodrato` for cases in which the latter is faster.

    void fmpz_mat_sqr_bodrato(fmpz_mat_t B, const fmpz_mat_t A)
    # Sets ``B`` to the square of the matrix ``A``, which must be
    # a square matrix. Aliasing is allowed.
    # The Bodrato algorithm is described in [Bodrato2010]_.
    # It is highly efficient for squaring matrices which satisfy both the
    # following conditions: (a) large elements,  (b) dimensions less than 150.

    void fmpz_mat_pow(fmpz_mat_t B, const fmpz_mat_t A, ulong e)
    # Sets ``B`` to the matrix ``A`` raised to the power ``e``,
    # where ``A`` must be a square matrix. Aliasing is allowed.

    void _fmpz_mat_mul_small(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
    # This internal function sets `C` to the matrix product `C = A B` computed
    # using classical matrix algorithm assuming that all entries of `A` and `B`
    # are small, that is, have bits `\le FLINT\_BITS - 2`. No aliasing is allowed.

    void _fmpz_mat_mul_double_word(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
    # This function is only for internal use and assumes that either:
    # - the entries of `A` and `B` are all nonnegative and strictly less than `2^{2*FLINT\_BITS}`, or
    # - the entries of `A` and `B` are all strictly less than `2^{2*FLINT\_BITS - 1}` in absolute value.

    void fmpz_mat_mul_fmpz_vec(fmpz * c, const fmpz_mat_t A, const fmpz * b, slong blen)
    void fmpz_mat_mul_fmpz_vec_ptr(fmpz * const * c, const fmpz_mat_t A, const fmpz * const * b, slong blen)
    # Compute a matrix-vector product of ``A`` and ``(b, blen)`` and store the result in ``c``.
    # The vector ``(b, blen)`` is either truncated or zero-extended to the number of columns of ``A``.
    # The number of entries written to ``c`` is always equal to the number of rows of ``A``.

    void fmpz_mat_fmpz_vec_mul(fmpz * c, const fmpz * a, slong alen, const fmpz_mat_t B)
    void fmpz_mat_fmpz_vec_mul_ptr(fmpz * const * c, const fmpz * const * a, slong alen, const fmpz_mat_t B)
    # Compute a vector-matrix product of ``(a, alen)`` and ``B`` and store the result in ``c``.
    # The vector ``(a, alen)`` is either truncated or zero-extended to the number of rows of ``B``.
    # The number of entries written to ``c`` is always equal to the number of columns of ``B``.

    int fmpz_mat_inv(fmpz_mat_t Ainv, fmpz_t den, const fmpz_mat_t A)
    # Sets (``Ainv``, ``den``) to the inverse matrix of ``A``.
    # Returns 1 if ``A`` is nonsingular and 0 if ``A`` is singular.
    # Aliasing of ``Ainv`` and ``A`` is allowed.
    # The denominator is not guaranteed to be minimal, but is guaranteed
    # to be a divisor of the determinant of ``A``.
    # This function uses a direct formula for matrices of size two or less,
    # and otherwise solves for the identity matrix using
    # fraction-free LU decomposition.

    void fmpz_mat_kronecker_product(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
    # Sets ``C`` to the Kronecker product of ``A`` and ``B``.

    void fmpz_mat_content(fmpz_t mat_gcd, const fmpz_mat_t A)
    # Sets ``mat_gcd`` as the gcd of all the elements of the matrix ``A``.
    # Returns 0 if the matrix is empty.

    void fmpz_mat_trace(fmpz_t trace, const fmpz_mat_t mat)
    # Computes the trace of the matrix, i.e. the sum of the entries on
    # the main diagonal. The matrix is required to be square.

    void fmpz_mat_det(fmpz_t det, const fmpz_mat_t A)
    # Sets ``det`` to the determinant of the square matrix `A`.
    # The matrix of dimension `0 \times 0` is defined to have determinant 1.
    # This function automatically chooses between :func:`fmpz_mat_det_cofactor`,
    # :func:`fmpz_mat_det_bareiss`, :func:`fmpz_mat_det_modular` and
    # :func:`fmpz_mat_det_modular_accelerated`
    # (with ``proved`` = 1), depending on the size of the matrix
    # and its entries.

    void fmpz_mat_det_cofactor(fmpz_t det, const fmpz_mat_t A)
    # Sets ``det`` to the determinant of the square matrix `A`
    # computed using direct cofactor expansion. This function only
    # supports matrices up to size `4 \times 4`.

    void fmpz_mat_det_bareiss(fmpz_t det, const fmpz_mat_t A)
    # Sets ``det`` to the determinant of the square matrix `A`
    # computed using the Bareiss algorithm. A copy of the input matrix is
    # row reduced using fraction-free Gaussian elimination, and the
    # determinant is read off from the last element on the main
    # diagonal.

    void fmpz_mat_det_modular(fmpz_t det, const fmpz_mat_t A, int proved)
    # Sets ``det`` to the determinant of the square matrix `A`
    # (if ``proved`` = 1), or a probabilistic value for the
    # determinant (``proved`` = 0), computed using a multimodular
    # algorithm.
    # The determinant is computed modulo several small primes and
    # reconstructed using the Chinese Remainder Theorem.
    # With ``proved`` = 1, sufficiently many primes are chosen
    # to satisfy the bound computed by ``fmpz_mat_det_bound``.
    # With ``proved`` = 0, the determinant is considered determined
    # if it remains unchanged modulo several consecutive primes
    # (currently if their product exceeds `2^{100}`).

    void fmpz_mat_det_modular_accelerated(fmpz_t det, const fmpz_mat_t A, int proved)
    # Sets ``det`` to the determinant of the square matrix `A`
    # (if ``proved`` = 1), or a probabilistic value for the
    # determinant (``proved`` = 0), computed using a multimodular
    # algorithm.
    # This function uses the same basic algorithm as ``fmpz_mat_det_modular``,
    # but instead of computing `\det(A)` directly, it generates a divisor `d`
    # of `\det(A)` and then computes `x = \det(A) / d` modulo several
    # small primes not dividing `d`. This typically accelerates the
    # computation by requiring fewer primes for large matrices, since `d`
    # with high probability will be nearly as large as the determinant.
    # This trick is described in [AbbottBronsteinMulders1999]_.

    void fmpz_mat_det_modular_given_divisor(fmpz_t det, const fmpz_mat_t A, const fmpz_t d, int proved)
    # Given a positive divisor `d` of `\det(A)`, sets ``det`` to the
    # determinant of the square matrix `A` (if ``proved`` = 1), or a
    # probabilistic value for the determinant (``proved`` = 0), computed
    # using a multimodular algorithm.

    void fmpz_mat_det_bound(fmpz_t bound, const fmpz_mat_t A)
    # Sets ``bound`` to a nonnegative integer `B` such that
    # `|\det(A)| \le B`. Assumes `A` to be a square matrix.
    # The bound is computed from the Hadamard inequality
    # `|\det(A)| \le \prod \|a_i\|_2` where the product is taken
    # over the rows `a_i` of `A`.

    void fmpz_mat_det_bound_nonzero(fmpz_t bound, const fmpz_mat_t A)
    # As per ``fmpz_mat_det_bound()`` but excludes zero columns. For use with
    # non-square matrices.

    void fmpz_mat_det_divisor(fmpz_t d, const fmpz_mat_t A)
    # Sets `d` to some positive divisor of the determinant of the given
    # square matrix `A`, if the determinant is nonzero. If `|\det(A)| = 0`,
    # `d` will always be set to zero.
    # A divisor is obtained by solving `Ax = b` for an arbitrarily chosen
    # right-hand side `b` using Dixon's algorithm and computing the least
    # common multiple of the denominators in `x`. This yields a divisor `d`
    # such that `|\det(A)| / d` is tiny with very high probability.

    void fmpz_mat_similarity(fmpz_mat_t A, slong r, fmpz_t d)
    # Applies a similarity transform to the `n\times n` matrix `M` in-place.
    # If `P` is the `n\times n` identity matrix the zero entries of whose row
    # `r` (`0`-indexed) have been replaced by `d`, this transform is equivalent
    # to `M = P^{-1}MP`.
    # Similarity transforms preserve the determinant, characteristic polynomial
    # and minimal polynomial.

    void _fmpz_mat_charpoly_berkowitz(fmpz * cp, const fmpz_mat_t mat)
    # Sets ``(cp, n+1)`` to the characteristic polynomial of
    # an `n \times n` square matrix.

    void fmpz_mat_charpoly_berkowitz(fmpz_poly_t cp, const fmpz_mat_t mat)
    # Computes the characteristic polynomial of length `n + 1` of
    # an `n \times n` square matrix. Uses an `O(n^4)` algorithm based on the
    # method of Berkowitz.

    void _fmpz_mat_charpoly_modular(fmpz * cp, const fmpz_mat_t mat)
    # Sets ``(cp, n+1)`` to the characteristic polynomial of
    # an `n \times n` square matrix.

    void fmpz_mat_charpoly_modular(fmpz_poly_t cp, const fmpz_mat_t mat)
    # Computes the characteristic polynomial of length `n + 1` of
    # an `n \times n` square matrix. Uses a modular method based on an `O(n^3)`
    # method over `\mathbb{Z}/n\mathbb{Z}`.

    void _fmpz_mat_charpoly(fmpz * cp, const fmpz_mat_t mat)
    # Sets ``(cp, n+1)`` to the characteristic polynomial of
    # an `n \times n` square matrix.

    void fmpz_mat_charpoly(fmpz_poly_t cp, const fmpz_mat_t mat)
    # Computes the characteristic polynomial of length `n + 1` of
    # an `n \times n` square matrix.

    slong _fmpz_mat_minpoly_modular(fmpz * cp, const fmpz_mat_t mat)
    # Sets ``(cp, n+1)`` to the modular polynomial of
    # an `n \times n` square matrix and returns its length.

    void fmpz_mat_minpoly_modular(fmpz_poly_t cp, const fmpz_mat_t mat)
    # Computes the minimal polynomial of an `n \times n` square matrix.
    # Uses a modular method based on an average time `O(n^3)`, worst case
    # `O(n^4)` method over `\mathbb{Z}/n\mathbb{Z}`.

    slong _fmpz_mat_minpoly(fmpz * cp, const fmpz_mat_t mat)
    # Sets ``cp`` to the minimal polynomial of an `n \times n` square
    # matrix and returns its length.

    void fmpz_mat_minpoly(fmpz_poly_t cp, const fmpz_mat_t mat)
    # Computes the minimal polynomial of an `n \times n` square matrix.

    slong fmpz_mat_rank(const fmpz_mat_t A)
    # Returns the rank, that is, the number of linearly independent columns
    # (equivalently, rows), of `A`. The rank is computed by row reducing
    # a copy of `A`.

    int fmpz_mat_col_partition(slong * part, fmpz_mat_t M, int short_circuit)
    # Returns the number `p` of distinct columns of `M` (or `0` if the flag
    # ``short_circuit`` is set and this number is greater than the number
    # of rows of `M`). The entries of array ``part`` are set to values in
    # `[0, p)` such that two entries of part are equal iff the corresponding
    # columns of `M` are equal. This function is used in van Hoeij polynomial
    # factoring.

    int fmpz_mat_solve(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B)
    # Solves the equation `AX = B` for nonsingular `A`. More precisely, computes
    # (``X``, ``den``) such that `AX = B \times \operatorname{den}`.
    # Returns 1 if `A` is nonsingular and 0 if `A` is singular.
    # The computed denominator will not generally be minimal.
    # This function uses Cramer's rule for small systems and
    # fraction-free LU decomposition followed by fraction-free forward
    # and back substitution for larger systems.
    # Note that for very large systems, it is faster to compute a modular
    # solution using ``fmpz_mat_solve_dixon``.

    int fmpz_mat_solve_fflu(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B)
    # Solves the equation `AX = B` for nonsingular `A`. More precisely, computes
    # (``X``, ``den``) such that `AX = B \times \operatorname{den}`.
    # Returns 1 if `A` is nonsingular and 0 if `A` is singular.
    # The computed denominator will not generally be minimal.
    # Uses fraction-free LU decomposition followed by fraction-free
    # forward and back substitution.

    int fmpz_mat_solve_fflu_precomp(fmpz_mat_t X, const slong * perm, const fmpz_mat_t FFLU, const fmpz_mat_t B)
    # Performs fraction-free forward and back substitution given a precomputed
    # fraction-free LU decomposition and corresponding permutation. If no
    # impossible division is encountered, the function returns `1`. This does not
    # mean the system has a solution, however a return value of `0` can only
    # occur if the system is insoluble.
    # If the return value is `1` and `r` is the rank of the matrix `A` whose FFLU
    # we have, then the first `r` rows of `p(A)y = p(b)d` hold, where `d` is the
    # denominator of the FFLU. The remaining rows must be checked by the caller.

    int fmpz_mat_solve_cramer(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B)
    # Solves the equation `AX = B` for nonsingular `A`. More precisely, computes
    # (``X``, ``den``) such that `AX = B \times \operatorname{den}`.
    # Returns 1 if `A` is nonsingular and 0 if `A` is singular.
    # Uses Cramer's rule. Only systems of size up to `3 \times 3` are allowed.

    void fmpz_mat_solve_bound(fmpz_t N, fmpz_t D, const fmpz_mat_t A, const fmpz_mat_t B)
    # Assuming that `A` is nonsingular, computes integers `N` and `D`
    # such that the reduced numerators and denominators `n/d` in
    # `A^{-1} B` satisfy the bounds `0 \le |n| \le N` and `0 \le d \le D`.

    int fmpz_mat_solve_dixon(fmpz_mat_t X, fmpz_t M, const fmpz_mat_t A, const fmpz_mat_t B)
    # Solves `AX = B` given a nonsingular square matrix `A` and a matrix `B` of
    # compatible dimensions, using a modular algorithm. In particular,
    # Dixon's p-adic lifting algorithm is used (currently a non-adaptive version).
    # This is generally the preferred method for large dimensions.
    # More precisely, this function computes an integer `M` and an integer
    # matrix `X` such that `AX = B \bmod M` and such that all the reduced
    # numerators and denominators of the elements `x = p/q` in the full
    # solution satisfy `2|p|q < M`. As such, the explicit rational solution
    # matrix can be recovered uniquely by passing the output of this
    # function to ``fmpq_mat_set_fmpz_mat_mod``.
    # A nonzero value is returned if `A` is nonsingular. If `A` is singular,
    # zero is returned and the values of the output variables will be
    # undefined.
    # Aliasing between input and output matrices is allowed.

    void _fmpz_mat_solve_dixon_den(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B, const nmod_mat_t Ainv, mp_limb_t p, const fmpz_t N, const fmpz_t D)
    # Solves the equation `AX = B` for nonsingular `A`. More precisely, computes
    # (``X``, ``den``) such that `AX = B \times \operatorname{den}` using a
    # ``p``-adic algorithm for the supplied prime ``p``. The values ``N`` and
    # ``D`` are absolute value bounds for the numerator and denominator of the
    # solution.
    # Uses the Dixon lifting algorithm with early termination once the lifting
    # stabilises.

    int fmpz_mat_solve_dixon_den(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B)
    # Solves the equation `AX = B` for nonsingular `A`. More precisely, computes
    # (``X``, ``den``) such that `AX = B \times \operatorname{den}`.
    # Returns 1 if `A` is nonsingular and 0 if `A` is singular.
    # The computed denominator will not generally be minimal.
    # Uses the Dixon lifting algorithm with early termination once the lifting
    # stabilises.

    int fmpz_mat_solve_multi_mod_den(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B)
    # Solves the equation `AX = B` for nonsingular `A`. More precisely, computes
    # (``X``, ``den``) such that `AX = B \times \operatorname{den}`.
    # Returns 1 if `A` is nonsingular and 0 if `A` is singular.
    # The computed denominator will not generally be minimal.
    # Uses a Chinese remainder algorithm with early termination once the lifting
    # stabilises.

    int fmpz_mat_can_solve_multi_mod_den(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B)
    # Returns `1` if the system `AX = B` can be solved. If so it computes
    # (``X``, ``den``) such that `AX = B \times \operatorname{den}`. The
    # computed denominator will not generally be minimal.
    # Uses a Chinese remainder algorithm.
    # Note that the matrices `A` and `B` may have any shape as long as they have
    # the same number of rows.

    int fmpz_mat_can_solve_fflu(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B)
    # Returns `1` if the system `AX = B` can be solved. If so it computes
    # (``X``, ``den``) such that `AX = B \times \operatorname{den}`. The
    # computed denominator will not generally be minimal.
    # Uses a fraction free LU decomposition algorithm.
    # Note that the matrices `A` and `B` may have any shape as long as they have
    # the same number of rows.

    int fmpz_mat_can_solve(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A, const fmpz_mat_t B)
    # Returns `1` if the system `AX = B` can be solved. If so it computes
    # (``X``, ``den``) such that `AX = B \times \operatorname{den}`. The
    # computed denominator will not generally be minimal.
    # Note that the matrices `A` and `B` may have any shape as long as they have
    # the same number of rows.

    slong fmpz_mat_find_pivot_any(const fmpz_mat_t mat, slong start_row, slong end_row, slong c)
    # Attempts to find a pivot entry for row reduction.
    # Returns a row index `r` between ``start_row`` (inclusive) and
    # ``stop_row`` (exclusive) such that column `c` in ``mat`` has
    # a nonzero entry on row `r`, or returns -1 if no such entry exists.
    # This implementation simply chooses the first nonzero entry
    # it encounters. This is likely to be a nearly optimal choice if all
    # entries in the matrix have roughly the same size, but can lead to
    # unnecessary coefficient growth if the entries vary in size.

    slong fmpz_mat_fflu(fmpz_mat_t B, fmpz_t den, slong * perm, const fmpz_mat_t A, int rank_check)
    # Uses fraction-free Gaussian elimination to set (``B``, ``den``) to a
    # fraction-free LU decomposition of ``A`` and returns the
    # rank of ``A``. Aliasing of ``A`` and ``B`` is allowed.
    # Pivot elements are chosen with ``fmpz_mat_find_pivot_any``.
    # If ``perm`` is non-``NULL``, the permutation of
    # rows in the matrix will also be applied to ``perm``.
    # If ``rank_check`` is set, the function aborts and returns 0 if the
    # matrix is detected not to have full rank without completing the
    # elimination.
    # The denominator ``den`` is set to `\pm \operatorname{det}(S)` where
    # `S` is an appropriate submatrix of `A` (`S = A` if `A` is square)
    # and the sign is decided by the parity of the permutation. Note that the
    # determinant is not generally the minimal denominator.
    # The fraction-free LU decomposition is defined in [NakTurWil1997]_.

    slong fmpz_mat_rref(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A)
    # Sets (``B``, ``den``) to the reduced row echelon form of ``A``
    # and returns the rank of ``A``. Aliasing of ``A`` and ``B``
    # is allowed.
    # The algorithm used chooses between ``fmpz_mat_rref_fflu`` and
    # ``fmpz_mat_rref_mul`` based on the dimensions of the input matrix.

    slong fmpz_mat_rref_fflu(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A)
    # Sets (``B``, ``den``) to the reduced row echelon form of ``A``
    # and returns the rank of ``A``. Aliasing of ``A`` and ``B``
    # is allowed.
    # The algorithm proceeds by first computing a row echelon form using
    # ``fmpz_mat_fflu``. Letting the upper part of this matrix be
    # `(U | V) P` where `U` is full rank upper triangular and `P` is a
    # permutation matrix, we obtain the rref by setting `V` to `U^{-1} V`
    # using back substitution. Scaling each completed row in the back
    # substitution to the denominator ``den``, we avoid introducing
    # new fractions. This strategy is equivalent to the fraction-free
    # Gauss-Jordan elimination in [NakTurWil1997]_, but faster since
    # only the part `V` corresponding to the null space has to be updated.
    # The denominator ``den`` is set to `\pm \operatorname{det}(S)` where
    # `S` is an appropriate submatrix of `A` (`S = A` if `A` is square).
    # Note that the determinant is not generally the minimal denominator.

    slong fmpz_mat_rref_mul(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A)
    # Sets (``B``, ``den``) to the reduced row echelon form of ``A``
    # and returns the rank of ``A``. Aliasing of ``A`` and ``B``
    # is allowed.
    # The algorithm works by computing the reduced row echelon form of ``A``
    # modulo a prime `p` using ``nmod_mat_rref``. The pivot columns and rows
    # of this matrix will then define a non-singular submatrix of ``A``,
    # nonsingular solving and matrix multiplication can then be used to determine
    # the reduced row echelon form of the whole of ``A``. This procedure is
    # described in [Stein2007]_.

    int fmpz_mat_is_in_rref_with_rank(const fmpz_mat_t A, const fmpz_t den, slong rank)
    # Checks that the matrix `A/den` is in reduced row echelon form of rank
    # ``rank``, returns 1 if so and 0 otherwise.

    slong fmpz_mat_rref_mod(slong * perm, fmpz_mat_t A, const fmpz_t p)
    # Uses fraction-free Gauss-Jordan elimination to set ``A``
    # to its reduced row echelon form and returns the rank of ``A``.
    # All computations are done modulo p.
    # Pivot elements are chosen with ``fmpz_mat_find_pivot_any``.
    # If ``perm`` is non-``NULL``, the permutation of
    # rows in the matrix will also be applied to ``perm``.

    void fmpz_mat_strong_echelon_form_mod(fmpz_mat_t A, const fmpz_t mod)
    # Transforms `A` such that `A` modulo ``mod`` is the strong echelon form
    # of the input matrix modulo ``mod``. The Howell form and the strong
    # echelon form are equal up to permutation of the rows, see [FieHof2014]_
    # for a definition of the strong echelon form and the algorithm used here.
    # `A` must have at least as many rows as columns.

    slong fmpz_mat_howell_form_mod(fmpz_mat_t A, const fmpz_t mod)
    # Transforms `A` such that `A` modulo ``mod`` is the Howell form of the
    # input matrix modulo ``mod``.
    # For a definition of the Howell form see [StoMul1998]_. The Howell form
    # is computed by first putting `A` into strong echelon form and then ordering
    # the rows.
    # `A` must have at least as many rows as columns.

    slong fmpz_mat_nullspace(fmpz_mat_t B, const fmpz_mat_t A)
    # Computes a basis for the right rational nullspace of `A` and returns
    # the dimension of the nullspace (or nullity). `B` is set to a matrix with
    # linearly independent columns and maximal rank such that `AB = 0`
    # (i.e. `Ab = 0` for each column `b` in `B`), and the rank of `B` is
    # returned.
    # In general, the entries in `B` will not be minimal: in particular,
    # the pivot entries in `B` will generally differ from unity.
    # `B` must be allocated with sufficient space to represent the result
    # (at most `n \times n` where `n` is the number of columns of `A`).

    slong fmpz_mat_rref_fraction_free(slong * perm, fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A)
    # Computes an integer matrix ``B`` and an integer ``den`` such that
    # ``B / den`` is the unique row reduced echelon form (RREF) of ``A``
    # and returns the rank, i.e. the number of nonzero rows in ``B``.
    # Aliasing of ``B`` and ``A`` is allowed, with an in-place
    # computation being more efficient. The size of ``B`` must be
    # the same as that of ``A``.
    # The permutation order will be written to ``perm`` unless this
    # argument is ``NULL``. That is, row ``i`` of the output matrix will
    # correspond to row ``perm[i]`` of the input matrix.
    # The denominator will always be a divisor of the determinant of (some
    # submatrix of) `A`, but is not guaranteed to be minimal or canonical in
    # any other sense.

    void fmpz_mat_hnf(fmpz_mat_t H, const fmpz_mat_t A)
    # Computes an integer matrix ``H`` such that ``H`` is the unique (row)
    # Hermite normal form of ``A``. The algorithm used is selected from the
    # implementations in FLINT to be the one most likely to be optimal, based on
    # the characteristics of the input matrix.
    # Aliasing of ``H`` and ``A`` is allowed. The size of ``H`` must be
    # the same as that of ``A``.

    void fmpz_mat_hnf_transform(fmpz_mat_t H, fmpz_mat_t U, const fmpz_mat_t A)
    # Computes an integer matrix ``H`` such that ``H`` is the unique (row)
    # Hermite normal form of ``A`` along with the transformation matrix
    # ``U`` such that `UA = H`. The algorithm used is selected from the
    # implementations in FLINT as per ``fmpz_mat_hnf``.
    # Aliasing of ``H`` and ``A`` is allowed. The size of ``H`` must be
    # the same as that of ``A`` and ``U`` must be square of \compatible
    # dimension (having the same number of rows as ``A``).

    void fmpz_mat_hnf_classical(fmpz_mat_t H, const fmpz_mat_t A)
    # Computes an integer matrix ``H`` such that ``H`` is the unique (row)
    # Hermite normal form of ``A``. The algorithm used is straightforward and
    # is described, for example, in [Algorithm 2.4.4] [Coh1996]_.
    # Aliasing of ``H`` and ``A`` is allowed. The size of ``H`` must be
    # the same as that of ``A``.

    void fmpz_mat_hnf_xgcd(fmpz_mat_t H, const fmpz_mat_t A)
    # Computes an integer matrix ``H`` such that ``H`` is the unique (row)
    # Hermite normal form of ``A``. The algorithm used is an improvement on the
    # basic algorithm and uses extended gcds to speed up computation, this method
    # is described, for example, in [Algorithm 2.4.5] [Coh1996]_.
    # Aliasing of ``H`` and ``A`` is allowed. The size of ``H`` must be
    # the same as that of ``A``.

    void fmpz_mat_hnf_modular(fmpz_mat_t H, const fmpz_mat_t A, const fmpz_t D)
    # Computes an integer matrix ``H`` such that ``H`` is the unique (row)
    # Hermite normal form of the `m\times n` matrix ``A``, where ``A`` is
    # assumed to be of rank `n` and ``D`` is known to be a positive multiple of
    # the determinant of the non-zero rows of ``H``. The algorithm used here is
    # due to Domich, Kannan and Trotter [DomKanTro1987]_ and is also described
    # in [Algorithm 2.4.8] [Coh1996]_.
    # Aliasing of ``H`` and ``A`` is allowed. The size of ``H`` must be
    # the same as that of ``A``.

    void fmpz_mat_hnf_modular_eldiv(fmpz_mat_t A, const fmpz_t D)
    # Transforms the `m\times n` matrix ``A`` into Hermite normal form,
    # where ``A`` is assumed to be of rank `n` and ``D`` is known to be a
    # positive multiple of the largest elementary divisor of ``A``.
    # The algorithm used here is described in [FieHof2014]_.

    void fmpz_mat_hnf_minors(fmpz_mat_t H, const fmpz_mat_t A)
    # Computes an integer matrix ``H`` such that ``H`` is the unique (row)
    # Hermite normal form of the `m\times n` matrix ``A``, where ``A`` is
    # assumed to be of rank `n`. The algorithm used here is due to Kannan and
    # Bachem [KanBac1979]_ and takes the principal minors to Hermite normal
    # form in turn.
    # Aliasing of ``H`` and ``A`` is allowed. The size of ``H`` must be
    # the same as that of ``A``.

    void fmpz_mat_hnf_pernet_stein(fmpz_mat_t H, const fmpz_mat_t A, flint_rand_t state)
    # Computes an integer matrix ``H`` such that ``H`` is the unique (row)
    # Hermite normal form of the `m\times n` matrix ``A``. The algorithm used
    # here is due to Pernet and Stein [PernetStein2010]_.
    # Aliasing of ``H`` and ``A`` is allowed. The size of ``H`` must be
    # the same as that of ``A``.

    int fmpz_mat_is_in_hnf(const fmpz_mat_t A)
    # Checks that the given matrix is in Hermite normal form, returns 1 if so and
    # 0 otherwise.

    void fmpz_mat_snf(fmpz_mat_t S, const fmpz_mat_t A)
    # Computes an integer matrix ``S`` such that ``S`` is the unique Smith
    # normal form of ``A``. The algorithm used is selected from the
    # implementations in FLINT to be the one most likely to be optimal, based on
    # the characteristics of the input matrix.
    # Aliasing of ``S`` and ``A`` is allowed. The size of ``S`` must be
    # the same as that of ``A``.

    void fmpz_mat_snf_diagonal(fmpz_mat_t S, const fmpz_mat_t A)
    # Computes an integer matrix ``S`` such that ``S`` is the unique Smith
    # normal form of the diagonal matrix ``A``. The algorithm used simply takes
    # gcds of pairs on the diagonal in turn until the Smith form is obtained.
    # Aliasing of ``S`` and ``A`` is allowed. The size of ``S`` must be
    # the same as that of ``A``.

    void fmpz_mat_snf_kannan_bachem(fmpz_mat_t S, const fmpz_mat_t A)
    # Computes an integer matrix ``S`` such that ``S`` is the unique Smith
    # normal form of the diagonal matrix ``A``. The algorithm used here is due
    # to Kannan and Bachem [KanBac1979]_
    # Aliasing of ``S`` and ``A`` is allowed. The size of ``S`` must be
    # the same as that of ``A``.

    void fmpz_mat_snf_iliopoulos(fmpz_mat_t S, const fmpz_mat_t A, const fmpz_t mod)
    # Computes an integer matrix ``S`` such that ``S`` is the unique Smith
    # normal form of the nonsingular `n\times n` matrix ``A``. The algorithm
    # used is due to Iliopoulos [Iliopoulos1989]_.
    # Aliasing of ``S`` and ``A`` is allowed. The size of ``S`` must be
    # the same as that of ``A``.

    int fmpz_mat_is_in_snf(const fmpz_mat_t A)
    # Checks that the given matrix is in Smith normal form, returns 1 if so and 0
    # otherwise.

    void fmpz_mat_gram(fmpz_mat_t B, const fmpz_mat_t A)
    # Sets ``B`` to the Gram matrix of the `m`-dimensional lattice ``L`` in
    # `n`-dimensional Euclidean space `R^n` spanned by the rows of
    # the `m \times n` matrix ``A``. Dimensions must be compatible.
    # ``A`` and ``B`` are allowed to be the same object if ``A`` is a
    # square matrix.

    int fmpz_mat_is_hadamard(const fmpz_mat_t H)
    # Returns nonzero iff `H` is a Hadamard matrix, meaning
    # that it is a square matrix, only has entries that are `\pm 1`,
    # and satisfies `H^T = n H^{-1}` where `n` is the matrix size.

    int fmpz_mat_hadamard(fmpz_mat_t H)
    # Attempts to set the matrix `H` to a Hadamard matrix, returning 1 if
    # successful and 0 if unsuccessful.
    # A Hadamard matrix of size `n` can only exist if `n` is 1, 2,
    # or a multiple of 4. It is not known whether a
    # Hadamard matrix exists for every size that is a multiple of 4.
    # This function uses the Paley construction, which
    # succeeds for all `n` of the form `n = 2^e` or `n = 2^e (q + 1)` where
    # `q` is an odd prime power. Orders `n` for which Hadamard matrices are
    # known to exist but for which this construction fails are
    # 92, 116, 156, ... (OEIS A046116).

    int fmpz_mat_get_d_mat(d_mat_t B, const fmpz_mat_t A)
    # Sets the entries of ``B`` as doubles corresponding to the entries of
    # ``A``, rounding down towards zero if the latter cannot be represented
    # exactly. The return value is -1 if any entry of ``A`` is too large to
    # fit in the normal range of a double, and 0 otherwise.

    int fmpz_mat_get_d_mat_transpose(d_mat_t B, const fmpz_mat_t A)
    # Sets the entries of ``B`` as doubles corresponding to the entries of
    # the transpose of ``A``, rounding down towards zero if the latter cannot
    # be represented exactly. The return value is -1 if any entry of ``A`` is
    # too large to fit in the normal range of a double, and 0 otherwise.

    void fmpz_mat_chol_d(d_mat_t R, const fmpz_mat_t A)
    # Computes ``R``, the Cholesky factor of a symmetric, positive definite
    # matrix ``A`` using the Cholesky decomposition process. (Sets ``R``
    # such that `A = RR^{T}` where ``R`` is a lower triangular matrix.)

    int fmpz_mat_is_reduced(const fmpz_mat_t A, double delta, double eta)
    int fmpz_mat_is_reduced_gram(const fmpz_mat_t A, double delta, double eta)
    # Returns a non-zero value if the basis ``A`` is LLL-reduced with factor
    # (``delta``, ``eta``), and otherwise returns zero.
    # The second version assumes ``A`` is the Gram matrix of the basis.

    int fmpz_mat_is_reduced_with_removal(const fmpz_mat_t A, double delta, double eta, const fmpz_t gs_B, int newd)
    int fmpz_mat_is_reduced_gram_with_removal(const fmpz_mat_t A, double delta, double eta, const fmpz_t gs_B, int newd)
    # Returns a non-zero value if the basis ``A`` is LLL-reduced with factor
    # (``delta``, ``eta``) for each of the first ``newd`` vectors and the squared
    # Gram-Schmidt length of each of the remaining `i`-th vectors
    # (where `i \ge` ``newd``) is greater than ``gs_B``, and otherwise returns zero.
    # The second version assumes ``A`` is the Gram matrix of the basis.

    void fmpz_mat_lll_original(fmpz_mat_t A, const fmpq_t delta, const fmpq_t eta)
    # Takes a basis `x_1, x_2, \ldots, x_m` of the lattice `L \subset R^n` (as
    # the rows of a `m \times n` matrix ``A``). The output is a (``delta``,
    # ``eta``)-reduced basis `y_1, y_2, \ldots, y_m` of the lattice `L` (as
    # the rows of the same `m \times n` matrix ``A``).

    void fmpz_mat_lll_storjohann(fmpz_mat_t A, const fmpq_t delta, const fmpq_t eta)
    # Takes a basis `x_1, x_2, \ldots, x_m` of the lattice `L \subset R^n` (as
    # the rows of a `m \times n` matrix ``A``). The output is an (``delta``,
    # ``eta``)-reduced basis `y_1, y_2, \ldots, y_m` of the lattice `L` (as
    # the rows of the same `m \times n` matrix ``A``). Uses a modified version of
    # LLL, which has better complexity in terms of the lattice dimension,
    # introduced by Storjohann.
    # See "Faster Algorithms for Integer Lattice Basis Reduction." Technical
    # Report 249. Zurich, Switzerland: Department Informatik, ETH. July 30,
    # 1996.

from .fmpz_mat_macros cimport *
