# distutils: libraries = flint
# distutils: depends = flint/arb_mat.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void arb_mat_init(arb_mat_t mat, long r, long c)
    # Initializes the matrix, setting it to the zero matrix with *r* rows
    # and *c* columns.

    void arb_mat_clear(arb_mat_t mat)
    # Clears the matrix, deallocating all entries.

    long arb_mat_allocated_bytes(const arb_mat_t x)
    # Returns the total number of bytes heap-allocated internally by this object.
    # The count excludes the size of the structure itself. Add
    # ``sizeof(arb_mat_struct)`` to get the size of the object as a whole.

    void arb_mat_window_init(arb_mat_t window, const arb_mat_t mat, long r1, long c1, long r2, long c2)
    # Initializes *window* to a window matrix into the submatrix of *mat*
    # starting at the corner at row *r1* and column *c1* (inclusive) and ending
    # at row *r2* and column *c2* (exclusive).

    void arb_mat_window_clear(arb_mat_t window)
    # Frees the window matrix.

    void arb_mat_set(arb_mat_t dest, const arb_mat_t src)

    void arb_mat_set_fmpz_mat(arb_mat_t dest, const fmpz_mat_t src)

    void arb_mat_set_round_fmpz_mat(arb_mat_t dest, const fmpz_mat_t src, long prec)

    void arb_mat_set_fmpq_mat(arb_mat_t dest, const fmpq_mat_t src, long prec)
    # Sets *dest* to *src*. The operands must have identical dimensions.

    void arb_mat_randtest(arb_mat_t mat, flint_rand_t state, long prec, long mag_bits)
    # Sets *mat* to a random matrix with up to *prec* bits of precision
    # and with exponents of width up to *mag_bits*.

    void arb_mat_printd(const arb_mat_t mat, long digits)
    # Prints each entry in the matrix with the specified number of decimal digits.

    void arb_mat_fprintd(FILE * file, const arb_mat_t mat, long digits)
    # Prints each entry in the matrix with the specified number of decimal
    # digits to the stream *file*.

    int arb_mat_equal(const arb_mat_t mat1, const arb_mat_t mat2)
    # Returns whether the matrices have the same dimensions
    # and identical intervals as entries.

    int arb_mat_overlaps(const arb_mat_t mat1, const arb_mat_t mat2)
    # Returns whether the matrices have the same dimensions
    # and each entry in *mat1* overlaps with the corresponding entry in *mat2*.

    int arb_mat_contains(const arb_mat_t mat1, const arb_mat_t mat2)

    int arb_mat_contains_fmpz_mat(const arb_mat_t mat1, const fmpz_mat_t mat2)

    int arb_mat_contains_fmpq_mat(const arb_mat_t mat1, const fmpq_mat_t mat2)
    # Returns whether the matrices have the same dimensions and each entry
    # in *mat2* is contained in the corresponding entry in *mat1*.

    int arb_mat_eq(const arb_mat_t mat1, const arb_mat_t mat2)
    # Returns whether *mat1* and *mat2* certainly represent the same matrix.

    int arb_mat_ne(const arb_mat_t mat1, const arb_mat_t mat2)
    # Returns whether *mat1* and *mat2* certainly do not represent the same matrix.

    int arb_mat_is_empty(const arb_mat_t mat)
    # Returns whether the number of rows or the number of columns in *mat* is zero.

    int arb_mat_is_square(const arb_mat_t mat)
    # Returns whether the number of rows is equal to the number of columns in *mat*.

    int arb_mat_is_exact(const arb_mat_t mat)
    # Returns whether all entries in *mat* have zero radius.

    int arb_mat_is_zero(const arb_mat_t mat)
    # Returns whether all entries in *mat* are exactly zero.

    int arb_mat_is_finite(const arb_mat_t mat)
    # Returns whether all entries in *mat* are finite.

    int arb_mat_is_triu(const arb_mat_t mat)
    # Returns whether *mat* is upper triangular; that is, all entries
    # below the main diagonal are exactly zero.

    int arb_mat_is_tril(const arb_mat_t mat)
    # Returns whether *mat* is lower triangular; that is, all entries
    # above the main diagonal are exactly zero.

    int arb_mat_is_diag(const arb_mat_t mat)
    # Returns whether *mat* is a diagonal matrix; that is, all entries
    # off the main diagonal are exactly zero.

    void arb_mat_zero(arb_mat_t mat)
    # Sets all entries in mat to zero.

    void arb_mat_one(arb_mat_t mat)
    # Sets the entries on the main diagonal to ones,
    # and all other entries to zero.

    void arb_mat_ones(arb_mat_t mat)
    # Sets all entries in the matrix to ones.

    void arb_mat_indeterminate(arb_mat_t mat)
    # Sets all entries in the matrix to indeterminate (NaN).

    void arb_mat_hilbert(arb_mat_t mat, long prec)
    # Sets *mat* to the Hilbert matrix, which has entries `A_{j,k} = 1/(j+k+1)`.

    void arb_mat_pascal(arb_mat_t mat, int triangular, long prec)
    # Sets *mat* to a Pascal matrix, whose entries are binomial coefficients.
    # If *triangular* is 0, constructs a full symmetric matrix
    # with the rows of Pascal's triangle as successive antidiagonals.
    # If *triangular* is 1, constructs the upper triangular matrix with
    # the rows of Pascal's triangle as columns, and if *triangular* is -1,
    # constructs the lower triangular matrix with the rows of Pascal's
    # triangle as rows.
    # The entries are computed using recurrence relations.
    # When the dimensions get large, some precision loss is possible; in that
    # case, the user may wish to create the matrix at slightly higher precision
    # and then round it to the final precision.

    void arb_mat_stirling(arb_mat_t mat, int kind, long prec)
    # Sets *mat* to a Stirling matrix, whose entries are Stirling numbers.
    # If *kind* is 0, the entries are set to the unsigned Stirling numbers
    # of the first kind. If *kind* is 1, the entries are set to the signed
    # Stirling numbers of the first kind. If *kind* is 2, the entries are
    # set to the Stirling numbers of the second kind.
    # The entries are computed using recurrence relations.
    # When the dimensions get large, some precision loss is possible; in that
    # case, the user may wish to create the matrix at slightly higher precision
    # and then round it to the final precision.

    void arb_mat_dct(arb_mat_t mat, int type, long prec)
    # Sets *mat* to the DCT (discrete cosine transform) matrix of order *n*
    # where *n* is the smallest dimension of *mat* (if *mat* is not square,
    # the matrix is extended periodically along the larger dimension).
    # There are many different conventions for defining DCT matrices; here,
    # we use the normalized "DCT-II" transform matrix
    # .. math ::
    # A_{j,k} = \sqrt{\frac{2}{n}} \cos\left(\frac{\pi j}{n} \left(k+\frac{1}{2}\right)\right)
    # which satisfies `A^{-1} = A^T`.
    # The *type* parameter is currently ignored and should be set to 0.
    # In the future, it might be used to select a different convention.

    void arb_mat_transpose(arb_mat_t dest, const arb_mat_t src)
    # Sets *dest* to the exact transpose *src*. The operands must have
    # compatible dimensions. Aliasing is allowed.

    void arb_mat_bound_inf_norm(mag_t b, const arb_mat_t A)
    # Sets *b* to an upper bound for the infinity norm (i.e. the largest
    # absolute value row sum) of *A*.

    void arb_mat_frobenius_norm(arb_t res, const arb_mat_t A, long prec)
    # Sets *res* to the Frobenius norm (i.e. the square root of the sum
    # of squares of entries) of *A*.

    void arb_mat_bound_frobenius_norm(mag_t res, const arb_mat_t A)
    # Sets *res* to an upper bound for the Frobenius norm of *A*.

    void arb_mat_neg(arb_mat_t dest, const arb_mat_t src)
    # Sets *dest* to the exact negation of *src*. The operands must have
    # the same dimensions.

    void arb_mat_add(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, long prec)
    # Sets res to the sum of *mat1* and *mat2*. The operands must have the same dimensions.

    void arb_mat_sub(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, long prec)
    # Sets *res* to the difference of *mat1* and *mat2*. The operands must have
    # the same dimensions.

    void arb_mat_mul_classical(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, long prec)

    void arb_mat_mul_threaded(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, long prec)

    void arb_mat_mul_block(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, long prec)

    void arb_mat_mul(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, long prec)
    # Sets *res* to the matrix product of *mat1* and *mat2*. The operands must have
    # compatible dimensions for matrix multiplication.
    # The *classical* version performs matrix multiplication in the trivial way.
    # The *block* version decomposes the input matrices into one or several
    # blocks of uniformly scaled matrices and multiplies
    # large blocks via *fmpz_mat_mul*. It also invokes
    # :func:`_arb_mat_addmul_rad_mag_fast` for the radius matrix multiplications.
    # The *threaded* version performs classical multiplication but splits the
    # computation over the number of threads returned by *flint_get_num_threads()*.
    # The default version chooses an algorithm automatically.

    void arb_mat_mul_entrywise(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, long prec)
    # Sets *C* to the entrywise product of *A* and *B*.
    # The operands must have the same dimensions.

    void arb_mat_sqr_classical(arb_mat_t B, const arb_mat_t A, long prec)

    void arb_mat_sqr(arb_mat_t res, const arb_mat_t mat, long prec)

    void arb_mat_pow_ui(arb_mat_t res, const arb_mat_t mat, unsigned long exp, long prec)
    # Sets *res* to *mat* raised to the power *exp*. Requires that *mat*
    # is a square matrix.

    void _arb_mat_addmul_rad_mag_fast(arb_mat_t C, mag_srcptr A, mag_srcptr B, long ar, long ac, long bc)
    # Helper function for matrix multiplication.
    # Adds to the radii of *C* the matrix product of the matrices represented
    # by *A* and *B*, where *A* is a linear array of coefficients in row-major
    # order and *B* is a linear array of coefficients in column-major order.
    # This function assumes that all exponents are small and is unsafe
    # for general use.

    void arb_mat_approx_mul(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, long prec)
    # Approximate matrix multiplication. The input radii are ignored and
    # the output matrix is set to an approximate floating-point result.
    # The radii in the output matrix will *not* necessarily be zeroed.

    void arb_mat_scalar_mul_2exp_si(arb_mat_t B, const arb_mat_t A, long c)
    # Sets *B* to *A* multiplied by `2^c`.

    void arb_mat_scalar_addmul_si(arb_mat_t B, const arb_mat_t A, long c, long prec)

    void arb_mat_scalar_addmul_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, long prec)

    void arb_mat_scalar_addmul_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, long prec)
    # Sets *B* to `B + A \times c`.

    void arb_mat_scalar_mul_si(arb_mat_t B, const arb_mat_t A, long c, long prec)

    void arb_mat_scalar_mul_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, long prec)

    void arb_mat_scalar_mul_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, long prec)
    # Sets *B* to `A \times c`.

    void arb_mat_scalar_div_si(arb_mat_t B, const arb_mat_t A, long c, long prec)

    void arb_mat_scalar_div_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, long prec)

    void arb_mat_scalar_div_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, long prec)
    # Sets *B* to `A / c`.

    int arb_mat_lu_classical(long * perm, arb_mat_t LU, const arb_mat_t A, long prec)

    int arb_mat_lu_recursive(long * perm, arb_mat_t LU, const arb_mat_t A, long prec)

    int arb_mat_lu(long * perm, arb_mat_t LU, const arb_mat_t A, long prec)
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

    void arb_mat_solve_tril_classical(arb_mat_t X, const arb_mat_t L, const arb_mat_t B, int unit, long prec)

    void arb_mat_solve_tril_recursive(arb_mat_t X, const arb_mat_t L, const arb_mat_t B, int unit, long prec)

    void arb_mat_solve_tril(arb_mat_t X, const arb_mat_t L, const arb_mat_t B, int unit, long prec)

    void arb_mat_solve_triu_classical(arb_mat_t X, const arb_mat_t U, const arb_mat_t B, int unit, long prec)

    void arb_mat_solve_triu_recursive(arb_mat_t X, const arb_mat_t U, const arb_mat_t B, int unit, long prec)

    void arb_mat_solve_triu(arb_mat_t X, const arb_mat_t U, const arb_mat_t B, int unit, long prec)
    # Solves the lower triangular system `LX = B` or the upper triangular system
    # `UX = B`, respectively. If *unit* is set, the main diagonal of *L* or *U*
    # is taken to consist of all ones, and in that case the actual entries on
    # the diagonal are not read at all and can contain other data.
    # The *classical* versions perform the computations iteratively while the
    # *recursive* versions perform the computations in a block recursive
    # way to benefit from fast matrix multiplication. The default versions
    # choose an algorithm automatically.

    void arb_mat_solve_lu_precomp(arb_mat_t X, const long * perm, const arb_mat_t LU, const arb_mat_t B, long prec)
    # Solves `AX = B` given the precomputed nonsingular LU decomposition `A = PLU`.
    # The matrices `X` and `B` are allowed to be aliased with each other,
    # but `X` is not allowed to be aliased with `LU`.

    int arb_mat_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, long prec)

    int arb_mat_solve_lu(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, long prec)

    int arb_mat_solve_precond(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, long prec)
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
    # the system [HS1967]_. This is usually several times slower than direct LU
    # decomposition, but the bounds do not blow up with *n* if the system is
    # well-conditioned. This algorithm is usually
    # the best choice for large systems at low to moderate precision.
    # * The default version selects between *lu* and *precomp* automatically.
    # The automatic choice should be reasonable most of the time, but users
    # may benefit from trying either *lu* or *precond* in specific applications.
    # For example, the *lu* solver often performs better for ill-conditioned
    # systems where use of very high precision is unavoidable.

    int arb_mat_solve_preapprox(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, const arb_mat_t R, const arb_mat_t T, long prec)
    # Solves `AX = B` where `A` is a nonsingular `n \times n` matrix
    # and `X` and `B` are `n \times m` matrices, given an approximation
    # `R` of the matrix inverse of `A`, and given the approximation `T`
    # of the solution `X`.
    # If `m > 0` and `A` cannot be inverted numerically (indicating either that
    # `A` is singular or that the precision is insufficient, or that `R` is
    # not a close enough approximation of the inverse of `A`), the values in the
    # output matrix are left undefined and zero is returned. A nonzero return
    # value guarantees that `A` is invertible and that the exact solution
    # matrix is contained in the output.

    int arb_mat_inv(arb_mat_t X, const arb_mat_t A, long prec)
    # Sets `X = A^{-1}` where `A` is a square matrix, computed by solving
    # the system `AX = I`.
    # If `A` cannot be inverted numerically (indicating either that
    # `A` is singular or that the precision is insufficient), the values in the
    # output matrix are left undefined and zero is returned.
    # A nonzero return value guarantees that the matrix is invertible
    # and that the exact inverse is contained in the output.

    void arb_mat_det_lu(arb_t det, const arb_mat_t A, long prec)

    void arb_mat_det_precond(arb_t det, const arb_mat_t A, long prec)

    void arb_mat_det(arb_t det, const arb_mat_t A, long prec)
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

    void arb_mat_approx_solve_triu(arb_mat_t X, const arb_mat_t U, const arb_mat_t B, int unit, long prec)

    void arb_mat_approx_solve_tril(arb_mat_t X, const arb_mat_t L, const arb_mat_t B, int unit, long prec)

    int arb_mat_approx_lu(long * P, arb_mat_t LU, const arb_mat_t A, long prec)

    void arb_mat_approx_solve_lu_precomp(arb_mat_t X, const long * perm, const arb_mat_t A, const arb_mat_t B, long prec)

    int arb_mat_approx_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, long prec)

    int arb_mat_approx_inv(arb_mat_t X, const arb_mat_t A, long prec)
    # These methods perform approximate solving *without any error control*.
    # The radii in the input matrices are ignored, the computations are done
    # numerically with floating-point arithmetic (using ordinary
    # Gaussian elimination and triangular solving, accelerated through
    # the use of block recursive strategies for large matrices), and the
    # output matrices are set to the approximate floating-point results with
    # zeroed error bounds.
    # Approximate solutions are useful for computing preconditioning matrices
    # for certified solutions. Some users may also find these methods useful
    # for doing ordinary numerical linear algebra in applications where
    # error bounds are not needed.

    int _arb_mat_cholesky_banachiewicz(arb_mat_t A, long prec)

    int arb_mat_cho(arb_mat_t L, const arb_mat_t A, long prec)
    # Computes the Cholesky decomposition of *A*, returning nonzero iff
    # the symmetric matrix defined by the lower triangular part of *A*
    # is certainly positive definite.
    # If a nonzero value is returned, then *L* is set to the lower triangular
    # matrix such that `A = L * L^T`.
    # If zero is returned, then either the matrix is not symmetric positive
    # definite, the input matrix was computed to insufficient precision,
    # or the decomposition was attempted at insufficient precision.
    # The underscore method computes *L* from *A* in-place, leaving the
    # strict upper triangular region undefined.

    void arb_mat_solve_cho_precomp(arb_mat_t X, const arb_mat_t L, const arb_mat_t B, long prec)
    # Solves `AX = B` given the precomputed Cholesky decomposition `A = L L^T`.
    # The matrices *X* and *B* are allowed to be aliased with each other,
    # but *X* is not allowed to be aliased with *L*.

    int arb_mat_spd_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, long prec)
    # Solves `AX = B` where *A* is a symmetric positive definite matrix
    # and *X* and *B* are `n \times m` matrices, using Cholesky decomposition.
    # If `m > 0` and *A* cannot be factored using Cholesky decomposition
    # (indicating either that *A* is not symmetric positive definite or that
    # the precision is insufficient), the values in the
    # output matrix are left undefined and zero is returned. A nonzero return
    # value guarantees that the symmetric matrix defined through the lower
    # triangular part of *A* is invertible and that the exact solution matrix
    # is contained in the output.

    void arb_mat_inv_cho_precomp(arb_mat_t X, const arb_mat_t L, long prec)
    # Sets `X = A^{-1}` where `A` is a symmetric positive definite matrix
    # whose Cholesky decomposition *L* has been computed with
    # :func:`arb_mat_cho`.
    # The inverse is calculated using the method of [Kri2013]_ which is more
    # efficient than solving `AX = I` with :func:`arb_mat_solve_cho_precomp`.

    int arb_mat_spd_inv(arb_mat_t X, const arb_mat_t A, long prec)
    # Sets `X = A^{-1}` where *A* is a symmetric positive definite matrix.
    # It is calculated using the method of [Kri2013]_ which computes fewer
    # intermediate results than solving `AX = I` with :func:`arb_mat_spd_solve`.
    # If *A* cannot be factored using Cholesky decomposition
    # (indicating either that *A* is not symmetric positive definite or that
    # the precision is insufficient), the values in the
    # output matrix are left undefined and zero is returned.  A nonzero return
    # value guarantees that the symmetric matrix defined through the lower
    # triangular part of *A* is invertible and that the exact inverse
    # is contained in the output.

    int _arb_mat_ldl_inplace(arb_mat_t A, long prec)

    int _arb_mat_ldl_golub_and_van_loan(arb_mat_t A, long prec)

    int arb_mat_ldl(arb_mat_t res, const arb_mat_t A, long prec)
    # Computes the `LDL^T` decomposition of *A*, returning nonzero iff
    # the symmetric matrix defined by the lower triangular part of *A*
    # is certainly positive definite.
    # If a nonzero value is returned, then *res* is set to a lower triangular
    # matrix that encodes the `L * D * L^T` decomposition of *A*.
    # In particular, `L` is a lower triangular matrix with ones on its diagonal
    # and whose strictly lower triangular region is the same as that of *res*.
    # `D` is a diagonal matrix with the same diagonal as that of *res*.
    # If zero is returned, then either the matrix is not symmetric positive
    # definite, the input matrix was computed to insufficient precision,
    # or the decomposition was attempted at insufficient precision.
    # The underscore methods compute *res* from *A* in-place, leaving the
    # strict upper triangular region undefined.
    # The default method uses algorithm 4.1.2 from [GVL1996]_.

    void arb_mat_solve_ldl_precomp(arb_mat_t X, const arb_mat_t L, const arb_mat_t B, long prec)
    # Solves `AX = B` given the precomputed `A = LDL^T` decomposition
    # encoded by *L*.  The matrices *X* and *B* are allowed to be aliased
    # with each other, but *X* is not allowed to be aliased with *L*.

    void arb_mat_inv_ldl_precomp(arb_mat_t X, const arb_mat_t L, long prec)
    # Sets `X = A^{-1}` where `A` is a symmetric positive definite matrix
    # whose `LDL^T` decomposition encoded by *L* has been computed with
    # :func:`arb_mat_ldl`.
    # The inverse is calculated using the method of [Kri2013]_ which is more
    # efficient than solving `AX = I` with :func:`arb_mat_solve_ldl_precomp`.

    void _arb_mat_charpoly(arb_ptr poly, const arb_mat_t mat, long prec)

    void arb_mat_charpoly(arb_poly_t poly, const arb_mat_t mat, long prec)
    # Sets *poly* to the characteristic polynomial of *mat* which must be
    # a square matrix. If the matrix has *n* rows, the underscore method
    # requires space for `n + 1` output coefficients.
    # Employs a division-free algorithm using `O(n^4)` operations.

    void _arb_mat_companion(arb_mat_t mat, arb_srcptr poly, long prec)

    void arb_mat_companion(arb_mat_t mat, const arb_poly_t poly, long prec)
    # Sets the *n* by *n* matrix *mat* to the companion matrix of the polynomial
    # *poly* which must have degree *n*.
    # The underscore method reads `n + 1` input coefficients.

    void arb_mat_exp_taylor_sum(arb_mat_t S, const arb_mat_t A, long N, long prec)
    # Sets *S* to the truncated exponential Taylor series `S = \sum_{k=0}^{N-1} A^k / k!`.
    # Uses rectangular splitting to compute the sum using `O(\sqrt{N})`
    # matrix multiplications. The recurrence relation for factorials
    # is used to get scalars that are small integers instead of full
    # factorials. As in [Joh2014b]_, all divisions are postponed to
    # the end by computing partial factorials of length `O(\sqrt{N})`.
    # The scalars could be reduced by doing more divisions, but this
    # appears to be slower in most cases.

    void arb_mat_exp(arb_mat_t B, const arb_mat_t A, long prec)
    # Sets *B* to the exponential of the matrix *A*, defined by the Taylor series
    # .. math ::
    # \exp(A) = \sum_{k=0}^{\infty} \frac{A^k}{k!}.
    # The function is evaluated as `\exp(A/2^r)^{2^r}`, where `r` is chosen
    # to give rapid convergence.
    # The elementwise error when truncating the Taylor series after *N*
    # terms is bounded by the error in the infinity norm, for which we have
    # .. math ::
    # \left\|\exp(2^{-r}A) - \sum_{k=0}^{N-1}
    # \frac{\left(2^{-r} A\right)^k}{k!} \right\|_{\infty} =
    # \left\|\sum_{k=N}^{\infty} \frac{\left(2^{-r} A\right)^k}{k!}\right\|_{\infty} \le
    # \sum_{k=N}^{\infty} \frac{(2^{-r} \|A\|_{\infty})^k}{k!}.
    # We bound the sum on the right using :func:`mag_exp_tail`.
    # Truncation error is not added to entries whose values are determined
    # by the sparsity structure of `A`.

    void arb_mat_trace(arb_t trace, const arb_mat_t mat, long prec)
    # Sets *trace* to the trace of the matrix, i.e. the sum of entries on the
    # main diagonal of *mat*. The matrix is required to be square.

    void _arb_mat_diag_prod(arb_t res, const arb_mat_t mat, long a, long b, long prec)

    void arb_mat_diag_prod(arb_t res, const arb_mat_t mat, long prec)
    # Sets *res* to the product of the entries on the main diagonal of *mat*.
    # The underscore method computes the product of the entries between
    # index *a* inclusive and *b* exclusive (the indices must be in range).

    void arb_mat_entrywise_is_zero(fmpz_mat_t dest, const arb_mat_t src)
    # Sets each entry of *dest* to indicate whether the corresponding
    # entry of *src* is certainly zero.
    # If the entry of *src* at row `i` and column `j` is zero according to
    # :func:`arb_is_zero` then the entry of *dest* at that row and column
    # is set to one, otherwise that entry of *dest* is set to zero.

    void arb_mat_entrywise_not_is_zero(fmpz_mat_t dest, const arb_mat_t src)
    # Sets each entry of *dest* to indicate whether the corresponding
    # entry of *src* is not certainly zero.
    # This the complement of :func:`arb_mat_entrywise_is_zero`.

    long arb_mat_count_is_zero(const arb_mat_t mat)
    # Returns the number of entries of *mat* that are certainly zero
    # according to :func:`arb_is_zero`.

    long arb_mat_count_not_is_zero(const arb_mat_t mat)
    # Returns the number of entries of *mat* that are not certainly zero.

    void arb_mat_get_mid(arb_mat_t B, const arb_mat_t A)
    # Sets the entries of *B* to the exact midpoints of the entries of *A*.

    void arb_mat_add_error_mag(arb_mat_t mat, const mag_t err)
    # Adds *err* in-place to the radii of the entries of *mat*.
