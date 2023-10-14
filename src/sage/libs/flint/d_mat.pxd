# distutils: libraries = flint
# distutils: depends = flint/d_mat.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void d_mat_init(d_mat_t mat, slong rows, slong cols)
    # Initialises a matrix with the given number of rows and columns for use.

    void d_mat_clear(d_mat_t mat)
    # Clears the given matrix.

    void d_mat_set(d_mat_t mat1, const d_mat_t mat2)
    # Sets ``mat1`` to a copy of ``mat2``. The dimensions of
    # ``mat1`` and ``mat2`` must be the same.

    void d_mat_swap(d_mat_t mat1, d_mat_t mat2)
    # Swaps two matrices. The dimensions of ``mat1`` and ``mat2``
    # are allowed to be different.

    void d_mat_swap_entrywise(d_mat_t mat1, d_mat_t mat2)
    # Swaps two matrices by swapping the individual entries rather than swapping
    # the contents of the structs.

    double d_mat_entry(d_mat_t mat, slong i, slong j)
    # Returns the entry of ``mat`` at row `i` and column `j`.
    # Both `i` and `j` must not exceed the dimensions of the matrix.
    # This function is implemented as a macro.

    double d_mat_get_entry(const d_mat_t mat, slong i, slong j)
    # Returns the entry of ``mat`` at row `i` and column `j`.
    # Both `i` and `j` must not exceed the dimensions of the matrix.

    double * d_mat_entry_ptr(const d_mat_t mat, slong i, slong j)
    # Returns a pointer to the entry of ``mat`` at row `i` and column
    # `j`. Both `i` and `j` must not exceed the dimensions of the matrix.

    void d_mat_zero(d_mat_t mat)
    # Sets all entries of ``mat`` to 0.

    void d_mat_one(d_mat_t mat)
    # Sets ``mat`` to the unit matrix, having ones on the main diagonal
    # and zeroes elsewhere. If ``mat`` is nonsquare, it is set to the
    # truncation of a unit matrix.

    void d_mat_randtest(d_mat_t mat, flint_rand_t state, slong minexp, slong maxexp)
    # Sets the entries of ``mat`` to random signed numbers with exponents
    # between ``minexp`` and ``maxexp`` or zero.

    void d_mat_print(const d_mat_t mat)
    # Prints the given matrix to the stream ``stdout``.

    int d_mat_equal(const d_mat_t mat1, const d_mat_t mat2)
    # Returns a non-zero value if ``mat1`` and ``mat2`` have
    # the same dimensions and entries, and zero otherwise.

    int d_mat_approx_equal(const d_mat_t mat1, const d_mat_t mat2, double eps)
    # Returns a non-zero value if ``mat1`` and ``mat2`` have
    # the same dimensions and entries within ``eps`` of each other,
    # and zero otherwise.

    int d_mat_is_zero(const d_mat_t mat)
    # Returns a non-zero value if all entries ``mat`` are zero, and
    # otherwise returns zero.

    int d_mat_is_approx_zero(const d_mat_t mat, double eps)
    # Returns a non-zero value if all entries ``mat`` are zero to within
    # ``eps`` and otherwise returns zero.

    int d_mat_is_empty(const d_mat_t mat)
    # Returns a non-zero value if the number of rows or the number of
    # columns in ``mat`` is zero, and otherwise returns
    # zero.

    int d_mat_is_square(const d_mat_t mat)
    # Returns a non-zero value if the number of rows is equal to the
    # number of columns in ``mat``, and otherwise returns zero.

    void d_mat_transpose(d_mat_t B, const d_mat_t A)
    # Sets `B` to `A^T`, the transpose of `A`. Dimensions must be compatible.
    # `A` and `B` are allowed to be the same object if `A` is a square matrix.

    void d_mat_mul_classical(d_mat_t C, const d_mat_t A, const d_mat_t B)
    # Sets ``C`` to the matrix product `C = A B`. The matrices must have
    # compatible dimensions for matrix multiplication (an exception is raised
    # otherwise). Aliasing is allowed.

    void d_mat_gso(d_mat_t B, const d_mat_t A)
    # Takes a subset of `R^m` `S = {a_1, a_2, \ldots, a_n}` (as the columns of
    # a `m \times n` matrix ``A``) and generates an orthonormal set
    # `S' = {b_1, b_2, \ldots, b_n}` (as the columns of the `m \times n` matrix
    # ``B``) that spans the same subspace of `R^m` as `S`.
    # This uses an algorithm of Schwarz-Rutishauser. See pp. 9 of
    # https://people.inf.ethz.ch/gander/papers/qrneu.pdf

    void d_mat_qr(d_mat_t Q, d_mat_t R, const d_mat_t A)
    # Computes the `QR` decomposition of a matrix ``A`` using the Gram-Schmidt
    # process. (Sets ``Q`` and ``R`` such that `A = QR` where ``R`` is
    # an upper triangular matrix and ``Q`` is an orthogonal matrix.)
    # This uses an algorithm of Schwarz-Rutishauser. See pp. 9 of
    # https://people.inf.ethz.ch/gander/papers/qrneu.pdf
