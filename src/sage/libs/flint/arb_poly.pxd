# distutils: libraries = flint
# distutils: depends = flint/arb_poly.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void arb_poly_init(arb_poly_t poly)
    # Initializes the polynomial for use, setting it to the zero polynomial.

    void arb_poly_clear(arb_poly_t poly)
    # Clears the polynomial, deallocating all coefficients and the
    # coefficient array.

    void arb_poly_fit_length(arb_poly_t poly, long len)
    # Makes sure that the coefficient array of the polynomial contains at
    # least *len* initialized coefficients.

    void _arb_poly_set_length(arb_poly_t poly, long len)
    # Directly changes the length of the polynomial, without allocating or
    # deallocating coefficients. The value should not exceed the allocation length.

    void _arb_poly_normalise(arb_poly_t poly)
    # Strips any trailing coefficients which are identical to zero.

    long arb_poly_allocated_bytes(const arb_poly_t x)
    # Returns the total number of bytes heap-allocated internally by this object.
    # The count excludes the size of the structure itself. Add
    # ``sizeof(arb_poly_struct)`` to get the size of the object as a whole.

    long arb_poly_length(const arb_poly_t poly)
    # Returns the length of *poly*, i.e. zero if *poly* is
    # identically zero, and otherwise one more than the index
    # of the highest term that is not identically zero.

    long arb_poly_degree(const arb_poly_t poly)
    # Returns the degree of *poly*, defined as one less than its length.
    # Note that if one or several leading coefficients are balls
    # containing zero, this value can be larger than the true
    # degree of the exact polynomial represented by *poly*,
    # so the return value of this function is effectively
    # an upper bound.

    int arb_poly_is_zero(const arb_poly_t poly)

    int arb_poly_is_one(const arb_poly_t poly)

    int arb_poly_is_x(const arb_poly_t poly)
    # Returns 1 if *poly* is exactly the polynomial 0, 1 or *x*
    # respectively. Returns 0 otherwise.

    void arb_poly_zero(arb_poly_t poly)

    void arb_poly_one(arb_poly_t poly)
    # Sets *poly* to the constant 0 respectively 1.

    void arb_poly_set(arb_poly_t dest, const arb_poly_t src)
    # Sets *dest* to a copy of *src*.

    void arb_poly_set_round(arb_poly_t dest, const arb_poly_t src, long prec)
    # Sets *dest* to a copy of *src*, rounded to *prec* bits.

    void arb_poly_set_trunc(arb_poly_t dest, const arb_poly_t src, long n)

    void arb_poly_set_trunc_round(arb_poly_t dest, const arb_poly_t src, long n, long prec)
    # Sets *dest* to a copy of *src*, truncated to length *n* and rounded to *prec* bits.

    void arb_poly_set_coeff_si(arb_poly_t poly, long n, long c)

    void arb_poly_set_coeff_arb(arb_poly_t poly, long n, const arb_t c)
    # Sets the coefficient with index *n* in *poly* to the value *c*.
    # We require that *n* is nonnegative.

    void arb_poly_get_coeff_arb(arb_t v, const arb_poly_t poly, long n)
    # Sets *v* to the value of the coefficient with index *n* in *poly*.
    # We require that *n* is nonnegative.

    void _arb_poly_shift_right(arb_ptr res, arb_srcptr poly, long len, long n)

    void arb_poly_shift_right(arb_poly_t res, const arb_poly_t poly, long n)
    # Sets *res* to *poly* divided by `x^n`, throwing away the lower coefficients.
    # We require that *n* is nonnegative.

    void _arb_poly_shift_left(arb_ptr res, arb_srcptr poly, long len, long n)

    void arb_poly_shift_left(arb_poly_t res, const arb_poly_t poly, long n)
    # Sets *res* to *poly* multiplied by `x^n`.
    # We require that *n* is nonnegative.

    void arb_poly_truncate(arb_poly_t poly, long n)
    # Truncates *poly* to have length at most *n*, i.e. degree
    # strictly smaller than *n*. We require that *n* is nonnegative.

    long arb_poly_valuation(const arb_poly_t poly)
    # Returns the degree of the lowest term that is not exactly zero in *poly*.
    # Returns -1 if *poly* is the zero polynomial.

    void arb_poly_set_fmpz_poly(arb_poly_t poly, const fmpz_poly_t src, long prec)

    void arb_poly_set_fmpq_poly(arb_poly_t poly, const fmpq_poly_t src, long prec)

    void arb_poly_set_si(arb_poly_t poly, long src)
    # Sets *poly* to *src*, rounding the coefficients to *prec* bits.

    void arb_poly_printd(const arb_poly_t poly, long digits)
    # Prints the polynomial as an array of coefficients, printing each
    # coefficient using *arb_printd*.

    void arb_poly_fprintd(FILE * file, const arb_poly_t poly, long digits)
    # Prints the polynomial as an array of coefficients to the stream *file*,
    # printing each coefficient using *arb_fprintd*.

    void arb_poly_randtest(arb_poly_t poly, flint_rand_t state, long len, long prec, long mag_bits)
    # Creates a random polynomial with length at most *len*.

    int arb_poly_contains(const arb_poly_t poly1, const arb_poly_t poly2)

    int arb_poly_contains_fmpz_poly(const arb_poly_t poly1, const fmpz_poly_t poly2)

    int arb_poly_contains_fmpq_poly(const arb_poly_t poly1, const fmpq_poly_t poly2)
    # Returns nonzero iff *poly1* contains *poly2*.

    int arb_poly_equal(const arb_poly_t A, const arb_poly_t B)
    # Returns nonzero iff *A* and *B* are equal as polynomial balls, i.e. all
    # coefficients have equal midpoint and radius.

    int _arb_poly_overlaps(arb_srcptr poly1, long len1, arb_srcptr poly2, long len2)

    int arb_poly_overlaps(const arb_poly_t poly1, const arb_poly_t poly2)
    # Returns nonzero iff *poly1* overlaps with *poly2*. The underscore
    # function requires that *len1* ist at least as large as *len2*.

    int arb_poly_get_unique_fmpz_poly(fmpz_poly_t z, const arb_poly_t x)
    # If *x* contains a unique integer polynomial, sets *z* to that value and returns
    # nonzero. Otherwise (if *x* represents no integers or more than one integer),
    # returns zero, possibly partially modifying *z*.

    void _arb_poly_majorant(arb_ptr res, arb_srcptr poly, long len, long prec)

    void arb_poly_majorant(arb_poly_t res, const arb_poly_t poly, long prec)
    # Sets *res* to an exact real polynomial whose coefficients are
    # upper bounds for the absolute values of the coefficients in *poly*,
    # rounded to *prec* bits.

    void _arb_poly_add(arb_ptr C, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long prec)
    # Sets *{C, max(lenA, lenB)}* to the sum of *{A, lenA}* and *{B, lenB}*.
    # Allows aliasing of the input and output operands.

    void arb_poly_add(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, long prec)

    void arb_poly_add_si(arb_poly_t C, const arb_poly_t A, long B, long prec)
    # Sets *C* to the sum of *A* and *B*.

    void _arb_poly_sub(arb_ptr C, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long prec)
    # Sets *{C, max(lenA, lenB)}* to the difference of *{A, lenA}* and *{B, lenB}*.
    # Allows aliasing of the input and output operands.

    void arb_poly_sub(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, long prec)
    # Sets *C* to the difference of *A* and *B*.

    void arb_poly_add_series(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, long len, long prec)
    # Sets *C* to the sum of *A* and *B*, truncated to length *len*.

    void arb_poly_sub_series(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, long len, long prec)
    # Sets *C* to the difference of *A* and *B*, truncated to length *len*.

    void arb_poly_neg(arb_poly_t C, const arb_poly_t A)
    # Sets *C* to the negation of *A*.

    void arb_poly_scalar_mul_2exp_si(arb_poly_t C, const arb_poly_t A, long c)
    # Sets *C* to *A* multiplied by `2^c`.

    void arb_poly_scalar_mul(arb_poly_t C, const arb_poly_t A, const arb_t c, long prec)
    # Sets *C* to *A* multiplied by *c*.

    void arb_poly_scalar_div(arb_poly_t C, const arb_poly_t A, const arb_t c, long prec)
    # Sets *C* to *A* divided by *c*.

    void _arb_poly_mullow_classical(arb_ptr C, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long n, long prec)

    void _arb_poly_mullow_block(arb_ptr C, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long n, long prec)

    void _arb_poly_mullow(arb_ptr C, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long n, long prec)
    # Sets *{C, n}* to the product of *{A, lenA}* and *{B, lenB}*, truncated to
    # length *n*. The output is not allowed to be aliased with either of the
    # inputs. We require `\mathrm{lenA} \ge \mathrm{lenB} > 0`,
    # `n > 0`, `\mathrm{lenA} + \mathrm{lenB} - 1 \ge n`.
    # The *classical* version uses a plain loop. This has good numerical
    # stability but gets slow for large *n*.
    # The *block* version decomposes the product into several
    # subproducts which are computed exactly over the integers.
    # It first attempts to find an integer `c`
    # such that `A(2^c x)` and `B(2^c x)` have slowly varying
    # coefficients, to reduce the number of blocks.
    # The scaling factor `c` is chosen in a quick, heuristic way
    # by picking the first and last nonzero terms in each polynomial.
    # If the indices in `A` are `a_2, a_1` and the log-2 magnitudes
    # are `e_2, e_1`, and the indices in `B` are `b_2, b_1`
    # with corresponding magnitudes `f_2, f_1`, then we compute
    # `c` as the weighted arithmetic mean of the slopes,
    # rounded to the nearest integer:
    # .. math ::
    # c = \left\lfloor
    # \frac{(e_2 - e_1) + (f_2 + f_1)}{(a_2 - a_1) + (b_2 - b_1)}
    # + \frac{1}{2}
    # \right \rfloor.
    # This strategy is used because it is simple. It is not optimal
    # in all cases, but will typically give good performance when
    # multiplying two power series with a similar decay rate.
    # The default algorithm chooses the *classical* algorithm for
    # short polynomials and the *block* algorithm for long polynomials.
    # If the input pointers are identical (and the lengths are the same),
    # they are assumed to represent the same polynomial, and its
    # square is computed.

    void arb_poly_mullow_classical(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, long n, long prec)

    void arb_poly_mullow_ztrunc(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, long n, long prec)

    void arb_poly_mullow_block(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, long n, long prec)

    void arb_poly_mullow(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, long n, long prec)
    # Sets *C* to the product of *A* and *B*, truncated to length *n*.
    # If the same variable is passed for *A* and *B*, sets *C* to the square
    # of *A* truncated to length *n*.

    void _arb_poly_mul(arb_ptr C, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long prec)
    # Sets *{C, lenA + lenB - 1}* to the product of *{A, lenA}* and *{B, lenB}*.
    # The output is not allowed to be aliased with either of the
    # inputs. We require `\mathrm{lenA} \ge \mathrm{lenB} > 0`.
    # This function is implemented as a simple wrapper for :func:`_arb_poly_mullow`.
    # If the input pointers are identical (and the lengths are the same),
    # they are assumed to represent the same polynomial, and its
    # square is computed.

    void arb_poly_mul(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, long prec)
    # Sets *C* to the product of *A* and *B*.
    # If the same variable is passed for *A* and *B*, sets *C* to the
    # square of *A*.

    void _arb_poly_inv_series(arb_ptr Q, arb_srcptr A, long Alen, long len, long prec)
    # Sets *{Q, len}* to the power series inverse of *{A, Alen}*. Uses Newton iteration.

    void arb_poly_inv_series(arb_poly_t Q, const arb_poly_t A, long n, long prec)
    # Sets *Q* to the power series inverse of *A*, truncated to length *n*.

    void _arb_poly_div_series(arb_ptr Q, arb_srcptr A, long Alen, arb_srcptr B, long Blen, long n, long prec)
    # Sets *{Q, n}* to the power series quotient of *{A, Alen}* by *{B, Blen}*.
    # Uses Newton iteration followed by multiplication.

    void arb_poly_div_series(arb_poly_t Q, const arb_poly_t A, const arb_poly_t B, long n, long prec)
    # Sets *Q* to the power series quotient *A* divided by *B*, truncated to length *n*.

    void _arb_poly_div(arb_ptr Q, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long prec)

    void _arb_poly_rem(arb_ptr R, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long prec)

    void _arb_poly_divrem(arb_ptr Q, arb_ptr R, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long prec)

    int arb_poly_divrem(arb_poly_t Q, arb_poly_t R, const arb_poly_t A, const arb_poly_t B, long prec)
    # Performs polynomial division with remainder, computing a quotient `Q` and
    # a remainder `R` such that `A = BQ + R`. The implementation reverses the
    # inputs and performs power series division.
    # If the leading coefficient of `B` contains zero (or if `B` is identically
    # zero), returns 0 indicating failure without modifying the outputs.
    # Otherwise returns nonzero.

    void _arb_poly_div_root(arb_ptr Q, arb_t R, arb_srcptr A, long len, const arb_t c, long prec)
    # Divides `A` by the polynomial `x - c`, computing the quotient `Q` as well
    # as the remainder `R = f(c)`.

    void _arb_poly_taylor_shift(arb_ptr g, const arb_t c, long n, long prec)
    void arb_poly_taylor_shift(arb_poly_t g, const arb_poly_t f, const arb_t c, long prec)
    # Sets *g* to the Taylor shift `f(x+c)`.
    # The underscore methods act in-place on *g* = *f* which has length *n*.

    void _arb_poly_compose(arb_ptr res, arb_srcptr poly1, long len1, arb_srcptr poly2, long len2, long prec)
    void arb_poly_compose(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, long prec)
    # Sets *res* to the composition `h(x) = f(g(x))` where `f` is given by
    # *poly1* and `g` is given by *poly2*.
    # The underscore method does not support aliasing of the output
    # with either input polynomial.

    void _arb_poly_compose_series(arb_ptr res, arb_srcptr poly1, long len1, arb_srcptr poly2, long len2, long n, long prec)
    void arb_poly_compose_series(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, long n, long prec)
    # Sets *res* to the power series composition `h(x) = f(g(x))` truncated
    # to order `O(x^n)` where `f` is given by *poly1* and `g` is given by *poly2*.
    # Wraps :func:`_gr_poly_compose_series` which chooses automatically
    # between various algorithms.
    # We require that the constant term in `g(x)` is exactly zero.
    # The underscore method does not support aliasing of the output
    # with either input polynomial.

    void _arb_poly_revert_series_lagrange(arb_ptr h, arb_srcptr f, long flen, long n, long prec)

    void arb_poly_revert_series_lagrange(arb_poly_t h, const arb_poly_t f, long n, long prec)

    void _arb_poly_revert_series_newton(arb_ptr h, arb_srcptr f, long flen, long n, long prec)

    void arb_poly_revert_series_newton(arb_poly_t h, const arb_poly_t f, long n, long prec)

    void _arb_poly_revert_series_lagrange_fast(arb_ptr h, arb_srcptr f, long flen, long n, long prec)

    void arb_poly_revert_series_lagrange_fast(arb_poly_t h, const arb_poly_t f, long n, long prec)

    void _arb_poly_revert_series(arb_ptr h, arb_srcptr f, long flen, long n, long prec)

    void arb_poly_revert_series(arb_poly_t h, const arb_poly_t f, long n, long prec)
    # Sets `h` to the power series reversion of `f`, i.e. the expansion
    # of the compositional inverse function `f^{-1}(x)`,
    # truncated to order `O(x^n)`, using respectively
    # Lagrange inversion, Newton iteration, fast Lagrange inversion,
    # and a default algorithm choice.
    # We require that the constant term in `f` is exactly zero and that the
    # linear term is nonzero. The underscore methods assume that *flen*
    # is at least 2, and do not support aliasing.

    void _arb_poly_evaluate_horner(arb_t y, arb_srcptr f, long len, const arb_t x, long prec)

    void arb_poly_evaluate_horner(arb_t y, const arb_poly_t f, const arb_t x, long prec)

    void _arb_poly_evaluate_rectangular(arb_t y, arb_srcptr f, long len, const arb_t x, long prec)

    void arb_poly_evaluate_rectangular(arb_t y, const arb_poly_t f, const arb_t x, long prec)

    void _arb_poly_evaluate(arb_t y, arb_srcptr f, long len, const arb_t x, long prec)

    void arb_poly_evaluate(arb_t y, const arb_poly_t f, const arb_t x, long prec)
    # Sets `y = f(x)`, evaluated respectively using Horner's rule,
    # rectangular splitting, and an automatic algorithm choice.

    void _arb_poly_evaluate_acb_horner(acb_t y, arb_srcptr f, long len, const acb_t x, long prec)

    void arb_poly_evaluate_acb_horner(acb_t y, const arb_poly_t f, const acb_t x, long prec)

    void _arb_poly_evaluate_acb_rectangular(acb_t y, arb_srcptr f, long len, const acb_t x, long prec)

    void arb_poly_evaluate_acb_rectangular(acb_t y, const arb_poly_t f, const acb_t x, long prec)

    void _arb_poly_evaluate_acb(acb_t y, arb_srcptr f, long len, const acb_t x, long prec)

    void arb_poly_evaluate_acb(acb_t y, const arb_poly_t f, const acb_t x, long prec)
    # Sets `y = f(x)` where `x` is a complex number, evaluating the
    # polynomial respectively using Horner's rule,
    # rectangular splitting, and an automatic algorithm choice.

    void _arb_poly_evaluate2_horner(arb_t y, arb_t z, arb_srcptr f, long len, const arb_t x, long prec)

    void arb_poly_evaluate2_horner(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, long prec)

    void _arb_poly_evaluate2_rectangular(arb_t y, arb_t z, arb_srcptr f, long len, const arb_t x, long prec)

    void arb_poly_evaluate2_rectangular(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, long prec)

    void _arb_poly_evaluate2(arb_t y, arb_t z, arb_srcptr f, long len, const arb_t x, long prec)

    void arb_poly_evaluate2(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, long prec)
    # Sets `y = f(x), z = f'(x)`, evaluated respectively using Horner's rule,
    # rectangular splitting, and an automatic algorithm choice.
    # When Horner's rule is used, the only advantage of evaluating the
    # function and its derivative simultaneously is that one does not have
    # to generate the derivative polynomial explicitly.
    # With the rectangular splitting algorithm, the powers can be reused,
    # making simultaneous evaluation slightly faster.

    void _arb_poly_evaluate2_acb_horner(acb_t y, acb_t z, arb_srcptr f, long len, const acb_t x, long prec)

    void arb_poly_evaluate2_acb_horner(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, long prec)

    void _arb_poly_evaluate2_acb_rectangular(acb_t y, acb_t z, arb_srcptr f, long len, const acb_t x, long prec)

    void arb_poly_evaluate2_acb_rectangular(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, long prec)

    void _arb_poly_evaluate2_acb(acb_t y, acb_t z, arb_srcptr f, long len, const acb_t x, long prec)

    void arb_poly_evaluate2_acb(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, long prec)
    # Sets `y = f(x), z = f'(x)`, evaluated respectively using Horner's rule,
    # rectangular splitting, and an automatic algorithm choice.

    void _arb_poly_product_roots(arb_ptr poly, arb_srcptr xs, long n, long prec)

    void arb_poly_product_roots(arb_poly_t poly, arb_srcptr xs, long n, long prec)
    # Generates the polynomial `(x-x_0)(x-x_1)\cdots(x-x_{n-1})`.

    void _arb_poly_product_roots_complex(arb_ptr poly, arb_srcptr r, long rn, acb_srcptr c, long cn, long prec)

    void arb_poly_product_roots_complex(arb_poly_t poly, arb_srcptr r, long rn, acb_srcptr c, long cn, long prec)
    # Generates the polynomial
    # .. math ::
    # \left(\prod_{i=0}^{rn-1} (x-r_i)\right) \left(\prod_{i=0}^{cn-1} (x-c_i)(x-\bar{c_i})\right)
    # having *rn* real roots given by the array *r* and having `2cn` complex roots
    # in conjugate pairs given by the length-*cn* array *c*.
    # Either *rn* or *cn* or both may be zero.
    # Note that only one representative from each complex conjugate pair
    # is supplied (unless a pair is supposed to
    # be repeated with higher multiplicity).
    # To construct a polynomial from complex roots where the conjugate pairs
    # have not been distinguished, use :func:`acb_poly_product_roots` instead.

    arb_ptr * _arb_poly_tree_alloc(long len)
    # Returns an initialized data structured capable of representing a
    # remainder tree (product tree) of *len* roots.

    void _arb_poly_tree_free(arb_ptr * tree, long len)
    # Deallocates a tree structure as allocated using *_arb_poly_tree_alloc*.

    void _arb_poly_tree_build(arb_ptr * tree, arb_srcptr roots, long len, long prec)
    # Constructs a product tree from a given array of *len* roots. The tree
    # structure must be pre-allocated to the specified length using
    # :func:`_arb_poly_tree_alloc`.

    void _arb_poly_evaluate_vec_iter(arb_ptr ys, arb_srcptr poly, long plen, arb_srcptr xs, long n, long prec)

    void arb_poly_evaluate_vec_iter(arb_ptr ys, const arb_poly_t poly, arb_srcptr xs, long n, long prec)
    # Evaluates the polynomial simultaneously at *n* given points, calling
    # :func:`_arb_poly_evaluate` repeatedly.

    void _arb_poly_evaluate_vec_fast_precomp(arb_ptr vs, arb_srcptr poly, long plen, arb_ptr * tree, long len, long prec)

    void _arb_poly_evaluate_vec_fast(arb_ptr ys, arb_srcptr poly, long plen, arb_srcptr xs, long n, long prec)

    void arb_poly_evaluate_vec_fast(arb_ptr ys, const arb_poly_t poly, arb_srcptr xs, long n, long prec)
    # Evaluates the polynomial simultaneously at *n* given points, using
    # fast multipoint evaluation.

    void _arb_poly_interpolate_newton(arb_ptr poly, arb_srcptr xs, arb_srcptr ys, long n, long prec)

    void arb_poly_interpolate_newton(arb_poly_t poly, arb_srcptr xs, arb_srcptr ys, long n, long prec)
    # Recovers the unique polynomial of length at most *n* that interpolates
    # the given *x* and *y* values. This implementation first interpolates in the
    # Newton basis and then converts back to the monomial basis.

    void _arb_poly_interpolate_barycentric(arb_ptr poly, arb_srcptr xs, arb_srcptr ys, long n, long prec)

    void arb_poly_interpolate_barycentric(arb_poly_t poly, arb_srcptr xs, arb_srcptr ys, long n, long prec)
    # Recovers the unique polynomial of length at most *n* that interpolates
    # the given *x* and *y* values. This implementation uses the barycentric
    # form of Lagrange interpolation.

    void _arb_poly_interpolation_weights(arb_ptr w, arb_ptr * tree, long len, long prec)

    void _arb_poly_interpolate_fast_precomp(arb_ptr poly, arb_srcptr ys, arb_ptr * tree, arb_srcptr weights, long len, long prec)

    void _arb_poly_interpolate_fast(arb_ptr poly, arb_srcptr xs, arb_srcptr ys, long len, long prec)

    void arb_poly_interpolate_fast(arb_poly_t poly, arb_srcptr xs, arb_srcptr ys, long n, long prec)
    # Recovers the unique polynomial of length at most *n* that interpolates
    # the given *x* and *y* values, using fast Lagrange interpolation.
    # The precomp function takes a precomputed product tree over the
    # *x* values and a vector of interpolation weights as additional inputs.

    void _arb_poly_derivative(arb_ptr res, arb_srcptr poly, long len, long prec)
    # Sets *{res, len - 1}* to the derivative of *{poly, len}*.
    # Allows aliasing of the input and output.

    void arb_poly_derivative(arb_poly_t res, const arb_poly_t poly, long prec)
    # Sets *res* to the derivative of *poly*.

    void _arb_poly_nth_derivative(arb_ptr res, arb_srcptr poly, unsigned long n, long len, long prec)
    # Sets *{res, len - n}* to the nth derivative of *{poly, len}*. Does
    # nothing if *len <= n*. Allows aliasing of the input and output.

    void arb_poly_nth_derivative(arb_poly_t res, const arb_poly_t poly, unsigned long n, long prec)
    # Sets *res* to the nth derivative of *poly*.

    void _arb_poly_integral(arb_ptr res, arb_srcptr poly, long len, long prec)
    # Sets *{res, len}* to the integral of *{poly, len - 1}*.
    # Allows aliasing of the input and output.

    void arb_poly_integral(arb_poly_t res, const arb_poly_t poly, long prec)
    # Sets *res* to the integral of *poly*.

    void _arb_poly_borel_transform(arb_ptr res, arb_srcptr poly, long len, long prec)

    void arb_poly_borel_transform(arb_poly_t res, const arb_poly_t poly, long prec)
    # Computes the Borel transform of the input polynomial, mapping `\sum_k a_k x^k`
    # to `\sum_k (a_k / k!) x^k`. The underscore method allows aliasing.

    void _arb_poly_inv_borel_transform(arb_ptr res, arb_srcptr poly, long len, long prec)

    void arb_poly_inv_borel_transform(arb_poly_t res, const arb_poly_t poly, long prec)
    # Computes the inverse Borel transform of the input polynomial, mapping `\sum_k a_k x^k`
    # to `\sum_k a_k k! x^k`. The underscore method allows aliasing.

    void _arb_poly_binomial_transform_basecase(arb_ptr b, arb_srcptr a, long alen, long len, long prec)

    void arb_poly_binomial_transform_basecase(arb_poly_t b, const arb_poly_t a, long len, long prec)

    void _arb_poly_binomial_transform_convolution(arb_ptr b, arb_srcptr a, long alen, long len, long prec)

    void arb_poly_binomial_transform_convolution(arb_poly_t b, const arb_poly_t a, long len, long prec)

    void _arb_poly_binomial_transform(arb_ptr b, arb_srcptr a, long alen, long len, long prec)

    void arb_poly_binomial_transform(arb_poly_t b, const arb_poly_t a, long len, long prec)
    # Computes the binomial transform of the input polynomial, truncating
    # the output to length *len*.
    # The binomial transform maps the coefficients `a_k` in the input polynomial
    # to the coefficients `b_k` in the output polynomial via
    # `b_n = \sum_{k=0}^n (-1)^k {n \choose k} a_k`.
    # The binomial transform is equivalent to the power series composition
    # `f(x) \to (1-x)^{-1} f(x/(x-1))`, and is its own inverse.
    # The *basecase* version evaluates coefficients one by one from the
    # definition, generating the binomial coefficients by a recurrence
    # relation.
    # The *convolution* version uses the identity
    # `T(f(x)) = B^{-1}(e^x B(f(-x)))` where `T` denotes the binomial
    # transform operator and `B` denotes the Borel transform operator.
    # This only costs a single polynomial multiplication, plus some
    # scalar operations.
    # The default version automatically chooses an algorithm.
    # The underscore methods do not support aliasing, and assume that
    # the lengths are nonzero.

    void _arb_poly_graeffe_transform(arb_ptr b, arb_srcptr a, long len, long prec)

    void arb_poly_graeffe_transform(arb_poly_t b, const arb_poly_t a, long prec)
    # Computes the Graeffe transform of input polynomial.
    # The Graeffe transform `G` of a polynomial `P` is defined through the
    # equation `G(x^2) = \pm P(x)P(-x)`.
    # The sign is given by `(-1)^d`, where `d = deg(P)`.
    # The Graeffe transform has the property that its roots are exactly the
    # squares of the roots of P.
    # The underscore method assumes that *a* and *b* are initialized,
    # *a* is of length *len*, and *b* is of length at least *len*.
    # Both methods allow aliasing.

    void _arb_poly_pow_ui_trunc_binexp(arb_ptr res, arb_srcptr f, long flen, unsigned long exp, long len, long prec)
    # Sets *{res, len}* to *{f, flen}* raised to the power *exp*, truncated
    # to length *len*. Requires that *len* is no longer than the length
    # of the power as computed without truncation (i.e. no zero-padding is performed).
    # Does not support aliasing of the input and output, and requires
    # that *flen* and *len* are positive.
    # Uses binary exponentiation.

    void arb_poly_pow_ui_trunc_binexp(arb_poly_t res, const arb_poly_t poly, unsigned long exp, long len, long prec)
    # Sets *res* to *poly* raised to the power *exp*, truncated to length *len*.
    # Uses binary exponentiation.

    void _arb_poly_pow_ui(arb_ptr res, arb_srcptr f, long flen, unsigned long exp, long prec)
    # Sets *res* to *{f, flen}* raised to the power *exp*. Does not
    # support aliasing of the input and output, and requires that
    # *flen* is positive.

    void arb_poly_pow_ui(arb_poly_t res, const arb_poly_t poly, unsigned long exp, long prec)
    # Sets *res* to *poly* raised to the power *exp*.

    void _arb_poly_pow_series(arb_ptr h, arb_srcptr f, long flen, arb_srcptr g, long glen, long len, long prec)
    # Sets *{h, len}* to the power series `f(x)^{g(x)} = \exp(g(x) \log f(x))` truncated
    # to length *len*. This function detects special cases such as *g* being an
    # exact small integer or `\pm 1/2`, and computes such powers more
    # efficiently. This function does not support aliasing of the output
    # with either of the input operands. It requires that all lengths
    # are positive, and assumes that *flen* and *glen* do not exceed *len*.

    void arb_poly_pow_series(arb_poly_t h, const arb_poly_t f, const arb_poly_t g, long len, long prec)
    # Sets *h* to the power series `f(x)^{g(x)} = \exp(g(x) \log f(x))` truncated
    # to length *len*. This function detects special cases such as *g* being an
    # exact small integer or `\pm 1/2`, and computes such powers more
    # efficiently.

    void _arb_poly_pow_arb_series(arb_ptr h, arb_srcptr f, long flen, const arb_t g, long len, long prec)
    # Sets *{h, len}* to the power series `f(x)^g = \exp(g \log f(x))` truncated
    # to length *len*. This function detects special cases such as *g* being an
    # exact small integer or `\pm 1/2`, and computes such powers more
    # efficiently. This function does not support aliasing of the output
    # with either of the input operands. It requires that all lengths
    # are positive, and assumes that *flen* does not exceed *len*.

    void arb_poly_pow_arb_series(arb_poly_t h, const arb_poly_t f, const arb_t g, long len, long prec)
    # Sets *h* to the power series `f(x)^g = \exp(g \log f(x))` truncated
    # to length *len*.

    void _arb_poly_sqrt_series(arb_ptr g, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_sqrt_series(arb_poly_t g, const arb_poly_t h, long n, long prec)
    # Sets *g* to the power series square root of *h*, truncated to length *n*.
    # Uses division-free Newton iteration for the reciprocal square root,
    # followed by a multiplication.
    # The underscore method does not support aliasing of the input and output
    # arrays. It requires that *hlen* and *n* are greater than zero.

    void _arb_poly_rsqrt_series(arb_ptr g, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_rsqrt_series(arb_poly_t g, const arb_poly_t h, long n, long prec)
    # Sets *g* to the reciprocal power series square root of *h*, truncated to length *n*.
    # Uses division-free Newton iteration.
    # The underscore method does not support aliasing of the input and output
    # arrays. It requires that *hlen* and *n* are greater than zero.

    void _arb_poly_log_series(arb_ptr res, arb_srcptr f, long flen, long n, long prec)

    void arb_poly_log_series(arb_poly_t res, const arb_poly_t f, long n, long prec)
    # Sets *res* to the power series logarithm of *f*, truncated to length *n*.
    # Uses the formula `\log(f(x)) = \int f'(x) / f(x) dx`, adding the logarithm of the
    # constant term in *f* as the constant of integration.
    # The underscore method supports aliasing of the input and output
    # arrays. It requires that *flen* and *n* are greater than zero.

    void _arb_poly_log1p_series(arb_ptr res, arb_srcptr f, long flen, long n, long prec)

    void arb_poly_log1p_series(arb_poly_t res, const arb_poly_t f, long n, long prec)
    # Computes the power series `\log(1+f)`, with better accuracy when the constant term of *f* is small.

    void _arb_poly_atan_series(arb_ptr res, arb_srcptr f, long flen, long n, long prec)

    void arb_poly_atan_series(arb_poly_t res, const arb_poly_t f, long n, long prec)

    void _arb_poly_asin_series(arb_ptr res, arb_srcptr f, long flen, long n, long prec)

    void arb_poly_asin_series(arb_poly_t res, const arb_poly_t f, long n, long prec)

    void _arb_poly_acos_series(arb_ptr res, arb_srcptr f, long flen, long n, long prec)

    void arb_poly_acos_series(arb_poly_t res, const arb_poly_t f, long n, long prec)
    # Sets *res* respectively to the power series inverse tangent,
    # inverse sine and inverse cosine of *f*, truncated to length *n*.
    # Uses the formulas
    # .. math ::
    # \tan^{-1}(f(x)) = \int f'(x) / (1+f(x)^2) dx,
    # \sin^{-1}(f(x)) = \int f'(x) / (1-f(x)^2)^{1/2} dx,
    # \cos^{-1}(f(x)) = -\int f'(x) / (1-f(x)^2)^{1/2} dx,
    # adding the inverse
    # function of the constant term in *f* as the constant of integration.
    # The underscore methods supports aliasing of the input and output
    # arrays. They require that *flen* and *n* are greater than zero.

    void _arb_poly_exp_series_basecase(arb_ptr f, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_exp_series_basecase(arb_poly_t f, const arb_poly_t h, long n, long prec)

    void _arb_poly_exp_series(arb_ptr f, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_exp_series(arb_poly_t f, const arb_poly_t h, long n, long prec)
    # Sets `f` to the power series exponential of `h`, truncated to length `n`.
    # The basecase version uses a simple recurrence for the coefficients,
    # requiring `O(nm)` operations where `m` is the length of `h`.
    # The main implementation uses Newton iteration, starting from a small
    # number of terms given by the basecase algorithm. The complexity
    # is `O(M(n))`. Redundant operations in the Newton iteration are
    # avoided by using the scheme described in [HZ2004]_.
    # The underscore methods support aliasing and allow the input to be
    # shorter than the output, but require the lengths to be nonzero.

    void _arb_poly_sin_cos_series(arb_ptr s, arb_ptr c, arb_srcptr h, long hlen, long n, long prec)
    void arb_poly_sin_cos_series(arb_poly_t s, arb_poly_t c, const arb_poly_t h, long n, long prec)
    # Sets *s* and *c* to the power series sine and cosine of *h*, computed
    # simultaneously.
    # The underscore method supports aliasing and requires the lengths to be nonzero.

    void _arb_poly_sin_series(arb_ptr s, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_sin_series(arb_poly_t s, const arb_poly_t h, long n, long prec)

    void _arb_poly_cos_series(arb_ptr c, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_cos_series(arb_poly_t c, const arb_poly_t h, long n, long prec)
    # Respectively evaluates the power series sine or cosine. These functions
    # simply wrap :func:`_arb_poly_sin_cos_series`. The underscore methods
    # support aliasing and require the lengths to be nonzero.

    void _arb_poly_tan_series(arb_ptr g, arb_srcptr h, long hlen, long len, long prec)

    void arb_poly_tan_series(arb_poly_t g, const arb_poly_t h, long n, long prec)
    # Sets *g* to the power series tangent of *h*.
    # For small *n* takes the quotient of the sine and cosine as computed
    # using the basecase algorithm. For large *n*, uses Newton iteration
    # to invert the inverse tangent series. The complexity is `O(M(n))`.
    # The underscore version does not support aliasing, and requires
    # the lengths to be nonzero.

    void _arb_poly_sin_cos_pi_series(arb_ptr s, arb_ptr c, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_sin_cos_pi_series(arb_poly_t s, arb_poly_t c, const arb_poly_t h, long n, long prec)

    void _arb_poly_sin_pi_series(arb_ptr s, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_sin_pi_series(arb_poly_t s, const arb_poly_t h, long n, long prec)

    void _arb_poly_cos_pi_series(arb_ptr c, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_cos_pi_series(arb_poly_t c, const arb_poly_t h, long n, long prec)

    void _arb_poly_cot_pi_series(arb_ptr c, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_cot_pi_series(arb_poly_t c, const arb_poly_t h, long n, long prec)
    # Compute the respective trigonometric functions of the input
    # multiplied by `\pi`.

    void _arb_poly_sinh_cosh_series_basecase(arb_ptr s, arb_ptr c, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_sinh_cosh_series_basecase(arb_poly_t s, arb_poly_t c, const arb_poly_t h, long n, long prec)

    void _arb_poly_sinh_cosh_series_exponential(arb_ptr s, arb_ptr c, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_sinh_cosh_series_exponential(arb_poly_t s, arb_poly_t c, const arb_poly_t h, long n, long prec)

    void _arb_poly_sinh_cosh_series(arb_ptr s, arb_ptr c, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_sinh_cosh_series(arb_poly_t s, arb_poly_t c, const arb_poly_t h, long n, long prec)

    void _arb_poly_sinh_series(arb_ptr s, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_sinh_series(arb_poly_t s, const arb_poly_t h, long n, long prec)

    void _arb_poly_cosh_series(arb_ptr c, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_cosh_series(arb_poly_t c, const arb_poly_t h, long n, long prec)
    # Sets *s* and *c* respectively to the hyperbolic sine and cosine of the
    # power series *h*, truncated to length *n*.
    # The implementations mirror those for sine and cosine, except that
    # the *exponential* version computes both functions using the exponential
    # function instead of the hyperbolic tangent.

    void _arb_poly_sinc_series(arb_ptr s, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_sinc_series(arb_poly_t s, const arb_poly_t h, long n, long prec)
    # Sets *c* to the sinc function of the power series *h*, truncated
    # to length *n*.

    void _arb_poly_sinc_pi_series(arb_ptr s, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_sinc_pi_series(arb_poly_t s, const arb_poly_t h, long n, long prec)
    # Compute the sinc function of the input multiplied by `\pi`.

    void _arb_poly_lambertw_series(arb_ptr res, arb_srcptr z, long zlen, int flags, long len, long prec)

    void arb_poly_lambertw_series(arb_poly_t res, const arb_poly_t z, int flags, long len, long prec)
    # Sets *res* to the Lambert W function of the power series *z*.
    # If *flags* is 0, the principal branch is computed; if *flags* is 1,
    # the second real branch `W_{-1}(z)` is computed.
    # The underscore method allows aliasing, but assumes that the lengths are nonzero.

    void _arb_poly_gamma_series(arb_ptr res, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_gamma_series(arb_poly_t res, const arb_poly_t h, long n, long prec)

    void _arb_poly_rgamma_series(arb_ptr res, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_rgamma_series(arb_poly_t res, const arb_poly_t h, long n, long prec)

    void _arb_poly_lgamma_series(arb_ptr res, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_lgamma_series(arb_poly_t res, const arb_poly_t h, long n, long prec)

    void _arb_poly_digamma_series(arb_ptr res, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_digamma_series(arb_poly_t res, const arb_poly_t h, long n, long prec)
    # Sets *res* to the series expansion of `\Gamma(h(x))`, `1/\Gamma(h(x))`,
    # or `\log \Gamma(h(x))`, `\psi(h(x))`, truncated to length *n*.
    # These functions first generate the Taylor series at the constant
    # term of *h*, and then call :func:`_arb_poly_compose_series`.
    # The Taylor coefficients are generated using the Riemann zeta function
    # if the constant term of *h* is a small integer,
    # and with Stirling's series otherwise.
    # The underscore methods support aliasing of the input and output
    # arrays, and require that *hlen* and *n* are greater than zero.

    void _arb_poly_rising_ui_series(arb_ptr res, arb_srcptr f, long flen, unsigned long r, long trunc, long prec)

    void arb_poly_rising_ui_series(arb_poly_t res, const arb_poly_t f, unsigned long r, long trunc, long prec)
    # Sets *res* to the rising factorial `(f) (f+1) (f+2) \cdots (f+r-1)`, truncated
    # to length *trunc*. The underscore method assumes that *flen*, *r* and *trunc*
    # are at least 1, and does not support aliasing. Uses binary splitting.

    void arb_poly_zeta_series(arb_poly_t res, const arb_poly_t s, const arb_t a, int deflate, long n, long prec)
    # Sets *res* to the Hurwitz zeta function `\zeta(s,a)` where `s` a power
    # series and `a` is a constant, truncated to length *n*.
    # To evaluate the usual Riemann zeta function, set `a = 1`.
    # If *deflate* is nonzero, evaluates `\zeta(s,a) + 1/(1-s)`, which
    # is well-defined as a limit when the constant term of `s` is 1.
    # In particular, expanding `\zeta(s,a) + 1/(1-s)` with `s = 1+x`
    # gives the Stieltjes constants
    # .. math ::
    # \sum_{k=0}^{n-1} \frac{(-1)^k}{k!} \gamma_k(a) x^k.
    # If `a = 1`, this implementation uses the reflection formula if the midpoint
    # of the constant term of `s` is negative.

    void _arb_poly_riemann_siegel_theta_series(arb_ptr res, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_riemann_siegel_theta_series(arb_poly_t res, const arb_poly_t h, long n, long prec)
    # Sets *res* to the series expansion of the Riemann-Siegel theta
    # function
    # .. math ::
    # \theta(h) = \arg \left(\Gamma\left(\frac{2ih+1}{4}\right)\right) - \frac{\log \pi}{2} h
    # where the argument of the gamma function is chosen continuously
    # as the imaginary part of the log gamma function.
    # The underscore method does not support aliasing of the input
    # and output arrays, and requires that the lengths are greater
    # than zero.

    void _arb_poly_riemann_siegel_z_series(arb_ptr res, arb_srcptr h, long hlen, long n, long prec)

    void arb_poly_riemann_siegel_z_series(arb_poly_t res, const arb_poly_t h, long n, long prec)
    # Sets *res* to the series expansion of the Riemann-Siegel Z-function
    # .. math ::
    # Z(h) = e^{i\theta(h)} \zeta(1/2+ih).
    # The zeros of the Z-function on the real line precisely
    # correspond to the imaginary parts of the zeros of
    # the Riemann zeta function on the critical line.
    # The underscore method supports aliasing of the input
    # and output arrays, and requires that the lengths are greater
    # than zero.

    void _arb_poly_root_bound_fujiwara(mag_t bound, arb_srcptr poly, long len)

    void arb_poly_root_bound_fujiwara(mag_t bound, arb_poly_t poly)
    # Sets *bound* to an upper bound for the magnitude of all the complex
    # roots of *poly*. Uses Fujiwara's bound
    # .. math ::
    # 2 \max \left\{\left|\frac{a_{n-1}}{a_n}\right|,
    # \left|\frac{a_{n-2}}{a_n}\right|^{1/2},
    # \cdots,
    # \left|\frac{a_1}{a_n}\right|^{1/(n-1)},
    # \left|\frac{a_0}{2a_n}\right|^{1/n}
    # \right\}
    # where `a_0, \ldots, a_n` are the coefficients of *poly*.

    void _arb_poly_newton_convergence_factor(arf_t convergence_factor, arb_srcptr poly, long len, const arb_t convergence_interval, long prec)
    # Given an interval `I` specified by *convergence_interval*, evaluates a bound
    # for `C = \sup_{t,u \in I} \frac{1}{2} |f''(t)| / |f'(u)|`,
    # where `f` is the polynomial defined by the coefficients *{poly, len}*.
    # The bound is obtained by evaluating `f'(I)` and `f''(I)` directly.
    # If `f` has large coefficients, `I` must be extremely precise in order to
    # get a finite factor.

    int _arb_poly_newton_step(arb_t xnew, arb_srcptr poly, long len, const arb_t x, const arb_t convergence_interval, const arf_t convergence_factor, long prec)
    # Performs a single step with Newton's method.
    # The input consists of the polynomial `f` specified by the coefficients
    # *{poly, len}*, an interval `x = [m-r, m+r]` known to contain a single root of `f`,
    # an interval `I` (*convergence_interval*) containing `x` with an
    # associated bound (*convergence_factor*) for
    # `C = \sup_{t,u \in I} \frac{1}{2} |f''(t)| / |f'(u)|`,
    # and a working precision *prec*.
    # The Newton update consists of setting
    # `x' = [m'-r', m'+r']` where `m' = m - f(m) / f'(m)`
    # and `r' = C r^2`. The expression `m - f(m) / f'(m)` is evaluated
    # using ball arithmetic at a working precision of *prec* bits, and the
    # rounding error during this evaluation is accounted for in the output.
    # We now check that `x' \in I` and `m' < m`. If both conditions are
    # satisfied, we set *xnew* to `x'` and return nonzero.
    # If either condition fails, we set *xnew* to `x` and return zero,
    # indicating that no progress was made.

    void _arb_poly_newton_refine_root(arb_t r, arb_srcptr poly, long len, const arb_t start, const arb_t convergence_interval, const arf_t convergence_factor, long eval_extra_prec, long prec)
    # Refines a precise estimate of a polynomial root to high precision
    # by performing several Newton steps, using nearly optimally
    # chosen doubling precision steps.
    # The inputs are defined as for *_arb_poly_newton_step*, except for
    # the precision parameters: *prec* is the target accuracy and
    # *eval_extra_prec* is the estimated number of guard bits that need
    # to be added to evaluate the polynomial accurately close to the root
    # (typically, if the polynomial has large coefficients of alternating
    # signs, this needs to be approximately the bit size of the coefficients).

    void _arb_poly_swinnerton_dyer_ui(arb_ptr poly, unsigned long n, long trunc, long prec)

    void arb_poly_swinnerton_dyer_ui(arb_poly_t poly, unsigned long n, long prec)
    # Computes the Swinnerton-Dyer polynomial `S_n`, which has degree `2^n`
    # and is the rational minimal polynomial of the sum
    # of the square roots of the first *n* prime numbers.
    # If *prec* is set to zero, a precision is chosen automatically such
    # that :func:`arb_poly_get_unique_fmpz_poly` should be successful.
    # Otherwise a working precision of *prec* bits is used.
    # The underscore version accepts an additional *trunc* parameter. Even
    # when computing a truncated polynomial, the array *poly* must have room for
    # `2^n + 1` coefficients, used as temporary space.
