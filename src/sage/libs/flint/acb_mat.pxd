# distutils: libraries = flint
# distutils: depends = flint/acb_mat.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void acb_mat_init(acb_mat_t mat, slong r, slong c)
    # Initializes the matrix, setting it to the zero matrix with *r* rows
    # and *c* columns.

    void acb_mat_clear(acb_mat_t mat)
    # Clears the matrix, deallocating all entries.

    slong acb_mat_allocated_bytes(const acb_mat_t x)
    # Returns the total number of bytes heap-allocated internally by this object.
    # The count excludes the size of the structure itself. Add
    # ``sizeof(acb_mat_struct)`` to get the size of the object as a whole.

    void acb_mat_window_init(acb_mat_t window, const acb_mat_t mat, slong r1, slong c1, slong r2, slong c2)
    # Initializes *window* to a window matrix into the submatrix of *mat*
    # starting at the corner at row *r1* and column *c1* (inclusive) and ending
    # at row *r2* and column *c2* (exclusive).

    void acb_mat_window_clear(acb_mat_t window)
    # Frees the window matrix.

    void acb_mat_set(acb_mat_t dest, const acb_mat_t src)

    void acb_mat_set_fmpz_mat(acb_mat_t dest, const fmpz_mat_t src)

    void acb_mat_set_round_fmpz_mat(acb_mat_t dest, const fmpz_mat_t src, slong prec)

    void acb_mat_set_fmpq_mat(acb_mat_t dest, const fmpq_mat_t src, slong prec)

    void acb_mat_set_arb_mat(acb_mat_t dest, const arb_mat_t src)

    void acb_mat_set_round_arb_mat(acb_mat_t dest, const arb_mat_t src, slong prec)
    # Sets *dest* to *src*. The operands must have identical dimensions.

    void acb_mat_randtest(acb_mat_t mat, flint_rand_t state, slong prec, slong mag_bits)
    # Sets *mat* to a random matrix with up to *prec* bits of precision
    # and with exponents of width up to *mag_bits*.

    void acb_mat_randtest_eig(acb_mat_t mat, flint_rand_t state, acb_srcptr E, slong prec)
    # Sets *mat* to a random matrix with the prescribed eigenvalues
    # supplied as the vector *E*. The output matrix is required to be
    # square. We generate a random unitary matrix via a matrix
    # exponential, and then evaluate an inverse Schur decomposition.

    void acb_mat_printd(const acb_mat_t mat, slong digits)
    # Prints each entry in the matrix with the specified number of decimal digits.

    void acb_mat_fprintd(FILE * file, const acb_mat_t mat, slong digits)
    # Prints each entry in the matrix with the specified number of decimal
    # digits to the stream *file*.

    int acb_mat_equal(const acb_mat_t mat1, const acb_mat_t mat2)
    # Returns whether the matrices have the same dimensions and identical
    # intervals as entries.

    int acb_mat_overlaps(const acb_mat_t mat1, const acb_mat_t mat2)
    # Returns whether the matrices have the same dimensions
    # and each entry in *mat1* overlaps with the corresponding entry in *mat2*.

    int acb_mat_contains(const acb_mat_t mat1, const acb_mat_t mat2)

    int acb_mat_contains_fmpz_mat(const acb_mat_t mat1, const fmpz_mat_t mat2)

    int acb_mat_contains_fmpq_mat(const acb_mat_t mat1, const fmpq_mat_t mat2)
    # Returns whether the matrices have the same dimensions and each entry
    # in *mat2* is contained in the corresponding entry in *mat1*.

    int acb_mat_eq(const acb_mat_t mat1, const acb_mat_t mat2)
    # Returns whether *mat1* and *mat2* certainly represent the same matrix.

    int acb_mat_ne(const acb_mat_t mat1, const acb_mat_t mat2)
    # Returns whether *mat1* and *mat2* certainly do not represent the same matrix.

    int acb_mat_is_real(const acb_mat_t mat)
    # Returns whether all entries in *mat* have zero imaginary part.

    int acb_mat_is_empty(const acb_mat_t mat)
    # Returns whether the number of rows or the number of columns in *mat* is zero.

    int acb_mat_is_square(const acb_mat_t mat)
    # Returns whether the number of rows is equal to the number of columns in *mat*.

    int acb_mat_is_exact(const acb_mat_t mat)
    # Returns whether all entries in *mat* have zero radius.

    int acb_mat_is_zero(const acb_mat_t mat)
    # Returns whether all entries in *mat* are exactly zero.

    int acb_mat_is_finite(const acb_mat_t mat)
    # Returns whether all entries in *mat* are finite.

    int acb_mat_is_triu(const acb_mat_t mat)
    # Returns whether *mat* is upper triangular; that is, all entries
    # below the main diagonal are exactly zero.

    int acb_mat_is_tril(const acb_mat_t mat)
    # Returns whether *mat* is lower triangular; that is, all entries
    # above the main diagonal are exactly zero.

    int acb_mat_is_diag(const acb_mat_t mat)
    # Returns whether *mat* is a diagonal matrix; that is, all entries
    # off the main diagonal are exactly zero.

    void acb_mat_zero(acb_mat_t mat)
    # Sets all entries in mat to zero.

    void acb_mat_one(acb_mat_t mat)
    # Sets the entries on the main diagonal to ones,
    # and all other entries to zero.

    void acb_mat_ones(acb_mat_t mat)
    # Sets all entries in the matrix to ones.

    void acb_mat_indeterminate(acb_mat_t mat)
    # Sets all entries in the matrix to indeterminate (NaN).

    void acb_mat_dft(acb_mat_t mat, int type, slong prec)
    # Sets *mat* to the DFT (discrete Fourier transform) matrix of order *n*
    # where *n* is the smallest dimension of *mat* (if *mat* is not square,
    # the matrix is extended periodically along the larger dimension).
    # Here, we use the normalized DFT matrix
    # .. math ::
    # A_{j,k} = \frac{\omega^{jk}}{\sqrt{n}}, \quad \omega = e^{-2\pi i/n}.
    # The *type* parameter is currently ignored and should be set to 0.
    # In the future, it might be used to select a different convention.

    void acb_mat_transpose(acb_mat_t dest, const acb_mat_t src)
    # Sets *dest* to the exact transpose *src*. The operands must have
    # compatible dimensions. Aliasing is allowed.

    void acb_mat_conjugate_transpose(acb_mat_t dest, const acb_mat_t src)
    # Sets *dest* to the conjugate transpose of *src*. The operands must have
    # compatible dimensions. Aliasing is allowed.

    void acb_mat_conjugate(acb_mat_t dest, const acb_mat_t src)
    # Sets *dest* to the elementwise complex conjugate of *src*.

    void acb_mat_bound_inf_norm(mag_t b, const acb_mat_t A)
    # Sets *b* to an upper bound for the infinity norm (i.e. the largest
    # absolute value row sum) of *A*.

    void acb_mat_frobenius_norm(arb_t res, const acb_mat_t A, slong prec)
    # Sets *res* to the Frobenius norm (i.e. the square root of the sum
    # of squares of entries) of *A*.

    void acb_mat_bound_frobenius_norm(mag_t res, const acb_mat_t A)
    # Sets *res* to an upper bound for the Frobenius norm of *A*.

    void acb_mat_neg(acb_mat_t dest, const acb_mat_t src)
    # Sets *dest* to the exact negation of *src*. The operands must have
    # the same dimensions.

    void acb_mat_add(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec)
    # Sets res to the sum of *mat1* and *mat2*. The operands must have the same dimensions.

    void acb_mat_sub(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec)
    # Sets *res* to the difference of *mat1* and *mat2*. The operands must have
    # the same dimensions.

    void acb_mat_mul_classical(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec)

    void acb_mat_mul_threaded(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec)

    void acb_mat_mul_reorder(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec)

    void acb_mat_mul(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec)
    # Sets *res* to the matrix product of *mat1* and *mat2*. The operands must have
    # compatible dimensions for matrix multiplication.
    # The *classical* version performs matrix multiplication in the trivial way.
    # The *threaded* version performs classical multiplication but splits the
    # computation over the number of threads returned by *flint_get_num_threads()*.
    # The *reorder* version reorders the data and performs one to four real
    # matrix multiplications via :func:`arb_mat_mul`.
    # The default version chooses an algorithm automatically.

    void acb_mat_mul_entrywise(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec)
    # Sets *res* to the entrywise product of *mat1* and *mat2*.
    # The operands must have the same dimensions.

    void acb_mat_sqr_classical(acb_mat_t res, const acb_mat_t mat, slong prec)

    void acb_mat_sqr(acb_mat_t res, const acb_mat_t mat, slong prec)
    # Sets *res* to the matrix square of *mat*. The operands must both be square
    # with the same dimensions.

    void acb_mat_pow_ui(acb_mat_t res, const acb_mat_t mat, ulong exp, slong prec)
    # Sets *res* to *mat* raised to the power *exp*. Requires that *mat*
    # is a square matrix.

    void acb_mat_approx_mul(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, slong prec)
    # Approximate matrix multiplication. The input radii are ignored and
    # the output matrix is set to an approximate floating-point result.
    # For performance reasons, the radii in the output matrix will *not*
    # necessarily be written (zeroed), but will remain zero if they
    # are already zeroed in *res* before calling this function.

    void acb_mat_scalar_mul_2exp_si(acb_mat_t B, const acb_mat_t A, slong c)
    # Sets *B* to *A* multiplied by `2^c`.

    void acb_mat_scalar_addmul_si(acb_mat_t B, const acb_mat_t A, slong c, slong prec)

    void acb_mat_scalar_addmul_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, slong prec)

    void acb_mat_scalar_addmul_arb(acb_mat_t B, const acb_mat_t A, const arb_t c, slong prec)

    void acb_mat_scalar_addmul_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, slong prec)
    # Sets *B* to `B + A \times c`.

    void acb_mat_scalar_mul_si(acb_mat_t B, const acb_mat_t A, slong c, slong prec)

    void acb_mat_scalar_mul_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, slong prec)

    void acb_mat_scalar_mul_arb(acb_mat_t B, const acb_mat_t A, const arb_t c, slong prec)

    void acb_mat_scalar_mul_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, slong prec)
    # Sets *B* to `A \times c`.

    void acb_mat_scalar_div_si(acb_mat_t B, const acb_mat_t A, slong c, slong prec)

    void acb_mat_scalar_div_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, slong prec)

    void acb_mat_scalar_div_arb(acb_mat_t B, const acb_mat_t A, const arb_t c, slong prec)

    void acb_mat_scalar_div_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, slong prec)
    # Sets *B* to `A / c`.

    int acb_mat_lu_classical(slong * perm, acb_mat_t LU, const acb_mat_t A, slong prec)

    int acb_mat_lu_recursive(slong * perm, acb_mat_t LU, const acb_mat_t A, slong prec)

    int acb_mat_lu(slong * perm, acb_mat_t LU, const acb_mat_t A, slong prec)
    # Given an `n \times n` matrix `A`, computes an LU decomposition `PLU = A`
    # using Gaussian elimination with partial pivoting.
    # The input and output matrices can be the same, performing the
    # decomposition in-place.
    # Entry `i` in the permutation vector perm is set to the row index in
    # the input matrix corresponding to row `i` in the output matrix.
    # The algorithm succeeds and returns nonzero if it can find `n` invertible
    # (i.e. not containing zero) pivot entries. This guarantees that the matrix
    # is invertible.
    # The algorithm fails and returns zero, leaving the entries in `P` and `LU`
    # undefined, if it cannot find `n` invertible pivot elements.
    # In this case, either the matrix is singular, the input matrix was
    # computed to insufficient precision, or the LU decomposition was
    # attempted at insufficient precision.
    # The *classical* version uses Gaussian elimination directly while
    # the *recursive* version performs the computation in a block recursive
    # way to benefit from fast matrix multiplication. The default version
    # chooses an algorithm automatically.

    void acb_mat_solve_tril_classical(acb_mat_t X, const acb_mat_t L, const acb_mat_t B, int unit, slong prec)

    void acb_mat_solve_tril_recursive(acb_mat_t X, const acb_mat_t L, const acb_mat_t B, int unit, slong prec)

    void acb_mat_solve_tril(acb_mat_t X, const acb_mat_t L, const acb_mat_t B, int unit, slong prec)

    void acb_mat_solve_triu_classical(acb_mat_t X, const acb_mat_t U, const acb_mat_t B, int unit, slong prec)

    void acb_mat_solve_triu_recursive(acb_mat_t X, const acb_mat_t U, const acb_mat_t B, int unit, slong prec)

    void acb_mat_solve_triu(acb_mat_t X, const acb_mat_t U, const acb_mat_t B, int unit, slong prec)
    # Solves the lower triangular system `LX = B` or the upper triangular system
    # `UX = B`, respectively. If *unit* is set, the main diagonal of *L* or *U*
    # is taken to consist of all ones, and in that case the actual entries on
    # the diagonal are not read at all and can contain other data.
    # The *classical* versions perform the computations iteratively while the
    # *recursive* versions perform the computations in a block recursive
    # way to benefit from fast matrix multiplication. The default versions
    # choose an algorithm automatically.

    void acb_mat_solve_lu_precomp(acb_mat_t X, const slong * perm, const acb_mat_t LU, const acb_mat_t B, slong prec)
    # Solves `AX = B` given the precomputed nonsingular LU decomposition `A = PLU`.
    # The matrices `X` and `B` are allowed to be aliased with each other,
    # but `X` is not allowed to be aliased with `LU`.

    int acb_mat_solve(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec)

    int acb_mat_solve_lu(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec)

    int acb_mat_solve_precond(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec)
    # Solves `AX = B` where `A` is a nonsingular `n \times n` matrix
    # and `X` and `B` are `n \times m` matrices.
    # If `m > 0` and `A` cannot be inverted numerically (indicating either that
    # `A` is singular or that the precision is insufficient), the values in the
    # output matrix are left undefined and zero is returned. A nonzero return
    # value guarantees that `A` is invertible and that the exact solution
    # matrix is contained in the output.
    # Three algorithms are provided:
    # * The *lu* version performs LU decomposition directly in ball arithmetic.
    # This is fast, but the bounds typically blow up exponentially with *n*,
    # even if the system is well-conditioned. This algorithm is usually
    # the best choice at very high precision.
    # * The *precond* version computes an approximate inverse to precondition
    # the system. This is usually several times slower than direct LU
    # decomposition, but the bounds do not blow up with *n* if the system is
    # well-conditioned. This algorithm is usually
    # the best choice for large systems at low to moderate precision.
    # * The default version selects between *lu* and *precomp* automatically.
    # The automatic choice should be reasonable most of the time, but users
    # may benefit from trying either *lu* or *precond* in specific applications.
    # For example, the *lu* solver often performs better for ill-conditioned
    # systems where use of very high precision is unavoidable.

    int acb_mat_inv(acb_mat_t X, const acb_mat_t A, slong prec)
    # Sets `X = A^{-1}` where `A` is a square matrix, computed by solving
    # the system `AX = I`.
    # If `A` cannot be inverted numerically (indicating either that
    # `A` is singular or that the precision is insufficient), the values in the
    # output matrix are left undefined and zero is returned.
    # A nonzero return value guarantees that the matrix is invertible
    # and that the exact inverse is contained in the output.

    void acb_mat_det_lu(acb_t det, const acb_mat_t A, slong prec)

    void acb_mat_det_precond(acb_t det, const acb_mat_t A, slong prec)

    void acb_mat_det(acb_t det, const acb_mat_t A, slong prec)
    # Sets *det* to the determinant of the matrix *A*.
    # The *lu* version uses Gaussian elimination with partial pivoting. If at
    # some point an invertible pivot element cannot be found, the elimination is
    # stopped and the magnitude of the determinant of the remaining submatrix
    # is bounded using Hadamard's inequality.
    # The *precond* version computes an approximate LU factorization of *A*
    # and multiplies by the inverse *L* and *U* martices as preconditioners
    # to obtain a matrix close to the identity matrix [Rum2010]_. An enclosure
    # for this determinant is computed using Gershgorin circles. This is about
    # four times slower than direct Gaussian elimination, but much more
    # numerically stable.
    # The default version automatically selects between the *lu* and *precond*
    # versions and additionally handles small or triangular matrices
    # by direct formulas.

    void acb_mat_approx_solve_triu(acb_mat_t X, const acb_mat_t U, const acb_mat_t B, int unit, slong prec)

    void acb_mat_approx_solve_tril(acb_mat_t X, const acb_mat_t L, const acb_mat_t B, int unit, slong prec)

    int acb_mat_approx_lu(slong * P, acb_mat_t LU, const acb_mat_t A, slong prec)

    void acb_mat_approx_solve_lu_precomp(acb_mat_t X, const slong * perm, const acb_mat_t A, const acb_mat_t B, slong prec)

    int acb_mat_approx_solve(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec)

    int acb_mat_approx_inv(acb_mat_t X, const acb_mat_t A, slong prec)
    # These methods perform approximate solving *without any error control*.
    # The radii in the input matrices are ignored, the computations are done
    # numerically with floating-point arithmetic (using ordinary
    # Gaussian elimination and triangular solving, accelerated through
    # the use of block recursive strategies for large matrices), and the
    # output matrices are set to the approximate floating-point results with
    # zeroed error bounds.

    void _acb_mat_charpoly(acb_ptr poly, const acb_mat_t mat, slong prec)

    void acb_mat_charpoly(acb_poly_t poly, const acb_mat_t mat, slong prec)
    # Sets *poly* to the characteristic polynomial of *mat* which must be
    # a square matrix. If the matrix has *n* rows, the underscore method
    # requires space for `n + 1` output coefficients.
    # Employs a division-free algorithm using `O(n^4)` operations.

    void _acb_mat_companion(acb_mat_t mat, acb_srcptr poly, slong prec)

    void acb_mat_companion(acb_mat_t mat, const acb_poly_t poly, slong prec)
    # Sets the *n* by *n* matrix *mat* to the companion matrix of the polynomial
    # *poly* which must have degree *n*.
    # The underscore method reads `n + 1` input coefficients.

    void acb_mat_exp_taylor_sum(acb_mat_t S, const acb_mat_t A, slong N, slong prec)
    # Sets *S* to the truncated exponential Taylor series `S = \sum_{k=0}^{N-1} A^k / k!`.
    # See :func:`arb_mat_exp_taylor_sum` for implementation notes.

    void acb_mat_exp(acb_mat_t B, const acb_mat_t A, slong prec)
    # Sets *B* to the exponential of the matrix *A*, defined by the Taylor series
    # .. math ::
    # \exp(A) = \sum_{k=0}^{\infty} \frac{A^k}{k!}.
    # The function is evaluated as `\exp(A/2^r)^{2^r}`, where `r` is chosen
    # to give rapid convergence of the Taylor series.
    # Error bounds are computed as for :func:`arb_mat_exp`.

    void acb_mat_trace(acb_t trace, const acb_mat_t mat, slong prec)
    # Sets *trace* to the trace of the matrix, i.e. the sum of entries on the
    # main diagonal of *mat*. The matrix is required to be square.

    void _acb_mat_diag_prod(acb_t res, const acb_mat_t mat, slong a, slong b, slong prec)

    void acb_mat_diag_prod(acb_t res, const acb_mat_t mat, slong prec)
    # Sets *res* to the product of the entries on the main diagonal of *mat*.
    # The underscore method computes the product of the entries between
    # index *a* inclusive and *b* exclusive (the indices must be in range).

    void acb_mat_get_mid(acb_mat_t B, const acb_mat_t A)
    # Sets the entries of *B* to the exact midpoints of the entries of *A*.

    void acb_mat_add_error_mag(acb_mat_t mat, const mag_t err)
    # Adds *err* in-place to the radii of the entries of *mat*.

    int acb_mat_approx_eig_qr(acb_ptr E, acb_mat_t L, acb_mat_t R, const acb_mat_t A, const mag_t tol, slong maxiter, slong prec)
    # Computes floating-point approximations of all the *n* eigenvalues
    # (and optionally eigenvectors) of the
    # given *n* by *n* matrix *A*. The approximations of the
    # eigenvalues are written to the vector *E*, in no particular order.
    # If *L* is not *NULL*, approximations of the corresponding left
    # eigenvectors are written to the rows of *L*. If *R* is not *NULL*,
    # approximations of the corresponding
    # right eigenvectors are written to the columns of *R*.
    # The parameters *tol* and *maxiter* can be used to control the
    # target numerical error and the maximum number of iterations
    # allowed before giving up. Passing *NULL* and 0 respectively results
    # in default values being used.
    # Uses the implicitly shifted QR algorithm with reduction
    # to Hessenberg form.
    # No guarantees are made about the accuracy of the output. A nonzero
    # return value indicates that the QR iteration converged numerically,
    # but this is only a heuristic termination test and does not imply
    # any statement whatsoever about error bounds.
    # The output may also be accurate even if this function returns zero.

    void acb_mat_eig_global_enclosure(mag_t eps, const acb_mat_t A, acb_srcptr E, const acb_mat_t R, slong prec)
    # Given an *n* by *n* matrix *A*, a length-*n* vector *E*
    # containing approximations of the eigenvalues of *A*,
    # and an *n* by *n* matrix *R* containing approximations of
    # the corresponding right eigenvectors, computes a rigorous bound
    # `\varepsilon` such that every eigenvalue `\lambda` of *A* satisfies
    # `|\lambda - \hat \lambda_k| \le \varepsilon`
    # for some `\hat \lambda_k` in *E*.
    # In other words, the union of the balls
    # `B_k = \{z : |z - \hat \lambda_k| \le \varepsilon\}` is guaranteed to
    # be an enclosure of all eigenvalues of *A*.
    # Note that there is no guarantee that each ball `B_k` can be
    # identified with a single eigenvalue: it is possible that some
    # balls contain several eigenvalues while other balls contain
    # no eigenvalues. In other words, this method is not powerful enough
    # to compute isolating balls for the individual eigenvalues (or even
    # for clusters of eigenvalues other than the whole spectrum).
    # Nevertheless, in practice the balls `B_k` will represent
    # eigenvalues one-to-one with high probability if the
    # given approximations are good.
    # The output can be used to certify
    # that all eigenvalues of *A* lie in some region of the complex
    # plane (such as a specific half-plane, strip, disk, or annulus)
    # without the need to certify the individual eigenvalues.
    # The output is easily converted into lower or upper bounds
    # for the absolute values or real or imaginary parts
    # of the spectrum, and with high probability these bounds will be tight.
    # Using :func:`acb_add_error_mag` and :func:`acb_union`, the output
    # can also be converted to a single :type:`acb_t` enclosing
    # the whole spectrum of *A* in a rectangle, but note that to
    # test whether a condition holds for all eigenvalues of *A*, it
    # is typically better to iterate over the individual balls `B_k`.
    # This function implements the fast algorithm in Theorem 1 in
    # [Miy2010]_ which extends the Bauer-Fike theorem. Approximations
    # *E* and *R* can, for instance, be computed using
    # :func:`acb_mat_approx_eig_qr`.
    # No assumptions are made about the structure of *A* or the
    # quality of the given approximations.

    void acb_mat_eig_enclosure_rump(acb_t lmbda, acb_mat_t J, acb_mat_t R, const acb_mat_t A, const acb_t lambda_approx, const acb_mat_t R_approx, slong prec)
    # Given an *n* by *n* matrix  *A* and an approximate
    # eigenvalue-eigenvector pair *lambda_approx* and *R_approx* (where
    # *R_approx* is an *n* by 1 matrix), computes an enclosure
    # *lambda* guaranteed to contain at least one of the eigenvalues of *A*,
    # along with an enclosure *R* for a corresponding right eigenvector.
    # More generally, this function can handle clustered (or repeated)
    # eigenvalues. If *R_approx* is an *n* by *k*
    # matrix containing approximate eigenvectors for a presumed cluster
    # of *k* eigenvalues near *lambda_approx*,
    # this function computes an enclosure *lambda*
    # guaranteed to contain at
    # least *k* eigenvalues of *A* along with a matrix *R* guaranteed to
    # contain a basis for the *k*-dimensional invariant subspace
    # associated with these eigenvalues.
    # Note that for multiple eigenvalues, determining the individual eigenvectors is
    # an ill-posed problem; describing an enclosure of the invariant subspace
    # is the best we can hope for.
    # For `k = 1`, it is guaranteed that `AR - R \lambda` contains
    # the zero matrix. For `k > 2`, this cannot generally be guaranteed
    # (in particular, *A* might not diagonalizable).
    # In this case, we can still compute an approximately diagonal
    # *k* by *k* interval matrix `J \approx \lambda I` such that `AR - RJ`
    # is guaranteed to contain the zero matrix.
    # This matrix has the property that the Jordan canonical form of
    # (any exact matrix contained in) *A* has a *k* by *k* submatrix
    # equal to the Jordan canonical form of
    # (some exact matrix contained in) *J*.
    # The output *J* is optional (the user can pass *NULL* to omit it).
    # The algorithm follows section 13.4 in [Rum2010]_, corresponding
    # to the ``verifyeig()`` routine in INTLAB.
    # The initial approximations can, for instance, be computed using
    # :func:`acb_mat_approx_eig_qr`.
    # No assumptions are made about the structure of *A* or the
    # quality of the given approximations.

    int acb_mat_eig_simple_rump(acb_ptr E, acb_mat_t L, acb_mat_t R, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, slong prec)

    int acb_mat_eig_simple_vdhoeven_mourrain(acb_ptr E, acb_mat_t L, acb_mat_t R, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, slong prec)

    int acb_mat_eig_simple(acb_ptr E, acb_mat_t L, acb_mat_t R, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, slong prec)
    # Computes all the eigenvalues (and optionally corresponding
    # eigenvectors) of the given *n* by *n* matrix *A*.
    # Attempts to prove that *A* has *n* simple (isolated)
    # eigenvalues, returning 1 if successful
    # and 0 otherwise. On success, isolating complex intervals for the
    # eigenvalues are written to the vector *E*, in no particular order.
    # If *L* is not *NULL*, enclosures of the corresponding left
    # eigenvectors are written to the rows of *L*. If *R* is not *NULL*,
    # enclosures of the corresponding
    # right eigenvectors are written to the columns of *R*.
    # The left eigenvectors are normalized so that `L = R^{-1}`.
    # This produces a diagonalization `LAR = D` where *D* is the
    # diagonal matrix with the entries in *E* on the diagonal.
    # The user supplies approximations *E_approx* and *R_approx*
    # of the eigenvalues and the right eigenvectors.
    # The initial approximations can, for instance, be computed using
    # :func:`acb_mat_approx_eig_qr`.
    # No assumptions are made about the structure of *A* or the
    # quality of the given approximations.
    # Two algorithms are implemented:
    # * The *rump* version calls :func:`acb_mat_eig_enclosure_rump` repeatedly
    # to certify eigenvalue-eigenvector pairs one by one. The iteration is
    # stopped to return non-success if a new eigenvalue overlaps with
    # previously computed one. Finally, *L* is computed by a matrix inversion.
    # This has complexity `O(n^4)`.
    # * The *vdhoeven_mourrain* version uses the algorithm in [HM2017]_ to
    # certify all eigenvalues and eigenvectors in one step. This has
    # complexity `O(n^3)`.
    # The default version currently uses *vdhoeven_mourrain*.
    # By design, these functions terminate instead of attempting to
    # compute eigenvalue clusters if some eigenvalues cannot be isolated.
    # To compute all eigenvalues of a matrix allowing for overlap,
    # :func:`acb_mat_eig_multiple_rump` may be used as a fallback,
    # or :func:`acb_mat_eig_multiple` may be used in the first place.

    int acb_mat_eig_multiple_rump(acb_ptr E, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, slong prec)

    int acb_mat_eig_multiple(acb_ptr E, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, slong prec)
    # Computes all the eigenvalues of the given *n* by *n* matrix *A*.
    # On success, the output vector *E* contains *n* complex intervals,
    # each representing one eigenvalue of *A* with the correct
    # multiplicities in case of overlap.
    # The output intervals are either disjoint or identical, and
    # identical intervals are guaranteed to be grouped consecutively.
    # Each complete run of *k* identical intervals thus represents a cluster of
    # exactly *k* eigenvalues which could not be separated from each
    # other at the current precision, but which could be isolated
    # from the other `n - k` eigenvalues of the matrix.
    # The user supplies approximations *E_approx* and *R_approx*
    # of the eigenvalues and the right eigenvectors.
    # The initial approximations can, for instance, be computed using
    # :func:`acb_mat_approx_eig_qr`.
    # No assumptions are made about the structure of *A* or the
    # quality of the given approximations.
    # The *rump* algorithm groups approximate eigenvalues that are close
    # and calls :func:`acb_mat_eig_enclosure_rump` repeatedly to validate
    # each cluster. The complexity is `O(m n^3)` for *m* clusters.
    # The default version, as currently implemented, first attempts to
    # call :func:`acb_mat_eig_simple_vdhoeven_mourrain` hoping that the
    # eigenvalues are actually simple. It then uses the *rump* algorithm as
    # a fallback.

from .acb_mat_macros cimport *
