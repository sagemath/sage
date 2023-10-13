# distutils: libraries = flint
# distutils: depends = flint/gr.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    long gr_ctx_sizeof_elem(gr_ctx_t ctx)
    # Return ``sizeof(type)`` where ``type`` is the underlying C
    # type for elements of *ctx*.

    int gr_ctx_clear(gr_ctx_t ctx)
    # Clears the context object *ctx*, freeing any memory
    # allocated by this object.
    # Some context objects may require that no elements are cleared after calling
    # this method, and may leak memory if not all elements have
    # been cleared when calling this method.
    # If *ctx* is derived from a base ring, the base ring context
    # may also be required to stay alive until after this
    # method is called.

    int gr_ctx_write(gr_stream_t out, gr_ctx_t ctx)
    int gr_ctx_print(gr_ctx_t ctx)
    int gr_ctx_println(gr_ctx_t ctx)
    int gr_ctx_get_str(char ** s, gr_ctx_t ctx)
    # Writes a description of the structure *ctx* to the stream *out*,
    # prints it to *stdout*, or sets *s* to a pointer to
    # a heap-allocated string of the description (the user must free
    # the string with ``flint_free``).
    # The *println* version prints a trailing newline.

    int gr_ctx_set_gen_name(gr_ctx_t ctx, const char * s)
    int gr_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)
    # Set the name of the generator (univariate polynomial ring,
    # finite field, etc.) or generators (multivariate).
    # The name is used when printing and may be used to choose
    # coercions.

    void gr_init(gr_ptr res, gr_ctx_t ctx)
    # Initializes *res* to a valid variable and sets it to the
    # zero element of the ring *ctx*.

    void gr_clear(gr_ptr res, gr_ctx_t ctx)
    # Clears *res*, freeing any memory allocated by this object.

    void gr_swap(gr_ptr x, gr_ptr y, gr_ctx_t ctx)
    # Swaps *x* and *y* efficiently.

    void gr_set_shallow(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    # Sets *res* to a shallow copy of *x*, copying the struct data.

    gr_ptr gr_heap_init(gr_ctx_t ctx)
    # Return a pointer to a single new heap-allocated element of *ctx*
    # set to 0.

    void gr_heap_clear(gr_ptr x, gr_ctx_t ctx)
    # Free the single heap-allocated element *x* of *ctx* which should
    # have been created with :func:`gr_heap_init`.

    gr_ptr gr_heap_init_vec(long len, gr_ctx_t ctx)
    # Return a pointer to a new heap-allocated vector of *len*
    # initialized elements.

    void gr_heap_clear_vec(gr_ptr x, long len, gr_ctx_t ctx)
    # Clear the *len* elements in the heap-allocated vector *len* and
    # free the vector itself.

    int gr_randtest(gr_ptr res, flint_rand_t state, gr_ctx_t ctx)
    # Sets *res* to a random element of the domain *ctx*.
    # The distribution is determined by the implementation.
    # Typically the distribution is non-uniform in order to
    # find corner cases more easily in test code.

    int gr_randtest_not_zero(gr_ptr res, flint_rand_t state, gr_ctx_t ctx)
    # Sets *res* to a random nonzero element of the domain *ctx*.
    # This operation will fail and return ``GR_DOMAIN`` in the zero ring.

    int gr_randtest_small(gr_ptr res, flint_rand_t state, gr_ctx_t ctx)
    # Sets *res* to a "small" element of the domain *ctx*.
    # This is suitable for randomized testing where a "large" argument
    # could result in excessive computation time.

    int gr_write(gr_stream_t out, gr_srcptr x, gr_ctx_t ctx)
    int gr_print(gr_srcptr x, gr_ctx_t ctx)
    int gr_println(gr_srcptr x, gr_ctx_t ctx)
    int gr_get_str(char ** s, gr_srcptr x, gr_ctx_t ctx)
    # Writes a description of the element *x* to the stream *out*,
    # or prints it to *stdout*, or sets *s* to a pointer to
    # a heap-allocated string of the description (the user must free
    # the string with ``flint_free``). The *println* version prints a
    # trailing newline.

    int gr_set_str(gr_ptr res, const char * x, gr_ctx_t ctx)
    # Sets *res* to the string description in *x*.

    int gr_write_n(gr_stream_t out, gr_srcptr x, long n, gr_ctx_t ctx)
    int gr_get_str_n(char ** s, gr_srcptr x, long n, gr_ctx_t ctx)
    # String conversion where real and complex numbers may be rounded
    # to *n* digits.

    int gr_set(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    # Sets *res* to a copy of the element *x*.

    int gr_set_other(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
    # Sets *res* to the element *x* of the structure *x_ctx* which
    # may be different from *ctx*. This returns the ``GR_DOMAIN`` flag
    # if *x* is not an element of *ctx* or cannot be converted
    # unambiguously to *ctx*.  The ``GR_UNABLE`` flag is returned
    # if the conversion is not implemented.

    int gr_set_ui(gr_ptr res, unsigned long x, gr_ctx_t ctx)
    int gr_set_si(gr_ptr res, long x, gr_ctx_t ctx)
    int gr_set_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx)
    int gr_set_fmpq(gr_ptr res, const fmpq_t x, gr_ctx_t ctx)
    int gr_set_d(gr_ptr res, double x, gr_ctx_t ctx)
    # Sets *res* to the value *x*. If no reasonable conversion to the
    # domain *ctx* is possible, returns ``GR_DOMAIN``.

    int gr_get_si(long * res, gr_srcptr x, gr_ctx_t ctx)
    int gr_get_ui(unsigned long * res, gr_srcptr x, gr_ctx_t ctx)
    int gr_get_fmpz(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
    int gr_get_fmpq(fmpq_t res, gr_srcptr x, gr_ctx_t ctx)
    int gr_get_d(double * res, gr_srcptr x, gr_ctx_t ctx)
    # Sets *res* to the value *x*. This returns the ``GR_DOMAIN`` flag
    # if *x* cannot be converted to the target type.
    # For floating-point output types, the output may be rounded.

    int gr_set_fmpz_2exp_fmpz(gr_ptr res, const fmpz_t x, const fmpz_t y, gr_ctx_t ctx)
    int gr_get_fmpz_2exp_fmpz(fmpz_t res1, fmpz_t res2, gr_srcptr x, gr_ctx_t ctx)
    # Set or retrieve a dyadic number.

    int gr_get_fexpr(fexpr_t res, gr_srcptr x, gr_ctx_t ctx)
    int gr_get_fexpr_serialize(fexpr_t res, gr_srcptr x, gr_ctx_t ctx)
    # Sets *res* to a symbolic expression representing *x*.
    # The *serialize* version may generate a representation of the
    # internal representation which is not intended to be human-readable.

    int gr_set_fexpr(gr_ptr res, fexpr_vec_t inputs, gr_vec_t outputs, const fexpr_t x, gr_ctx_t ctx)
    # Sets *res* to the evaluation of the expression *x* in the
    # given ring or structure.
    # The user must provide vectors *inputs* and *outputs* which
    # may be empty initially and which may be used as scratch space
    # during evaluation. Non-empty vectors may be given to map symbols
    # to predefined values.

    int gr_zero(gr_ptr res, gr_ctx_t ctx)
    int gr_one(gr_ptr res, gr_ctx_t ctx)
    int gr_neg_one(gr_ptr res, gr_ctx_t ctx)
    # Sets *res* to the ring element 0, 1 or -1.

    int gr_gen(gr_ptr res, gr_ctx_t ctx)
    # Sets *res* to a generator of this domain. The meaning of
    # "generator" depends on the domain.

    int gr_gens(gr_vec_t res, gr_ctx_t ctx)
    # Sets *res* to a vector containing the generators of this domain
    # where this makes sense, for example in a multivariate polynomial
    # ring.

    truth_t gr_is_zero(gr_srcptr x, gr_ctx_t ctx)
    truth_t gr_is_one(gr_srcptr x, gr_ctx_t ctx)
    truth_t gr_is_neg_one(gr_srcptr x, gr_ctx_t ctx)
    # Returns whether *x* is equal to the ring element 0, 1 or -1,
    # respectively.

    truth_t gr_equal(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
    # Returns whether the elements *x* and *y* are equal.

    truth_t gr_is_integer(gr_srcptr x, gr_ctx_t ctx)
    # Returns whether *x* represents an integer.

    truth_t gr_is_rational(gr_srcptr x, gr_ctx_t ctx)
    # Returns whether *x* represents a rational number.

    int gr_neg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    # Sets *res* to `-x`.

    int gr_add(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
    int gr_add_ui(gr_ptr res, gr_srcptr x, unsigned long y, gr_ctx_t ctx)
    int gr_add_si(gr_ptr res, gr_srcptr x, long y, gr_ctx_t ctx)
    int gr_add_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
    int gr_add_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
    int gr_add_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
    int gr_other_add(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)
    # Sets *res* to `x + y`.

    int gr_sub(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
    int gr_sub_ui(gr_ptr res, gr_srcptr x, unsigned long y, gr_ctx_t ctx)
    int gr_sub_si(gr_ptr res, gr_srcptr x, long y, gr_ctx_t ctx)
    int gr_sub_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
    int gr_sub_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
    int gr_sub_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
    int gr_other_sub(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)
    # Sets *res* to `x - y`.

    int gr_mul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
    int gr_mul_ui(gr_ptr res, gr_srcptr x, unsigned long y, gr_ctx_t ctx)
    int gr_mul_si(gr_ptr res, gr_srcptr x, long y, gr_ctx_t ctx)
    int gr_mul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
    int gr_mul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
    int gr_mul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
    int gr_other_mul(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)
    # Sets *res* to `x \cdot y`.

    int gr_addmul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
    int gr_addmul_ui(gr_ptr res, gr_srcptr x, unsigned long y, gr_ctx_t ctx)
    int gr_addmul_si(gr_ptr res, gr_srcptr x, long y, gr_ctx_t ctx)
    int gr_addmul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
    int gr_addmul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
    int gr_addmul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
    # Sets *res* to `\mathrm{res } + x \cdot y`.
    # Rings may override the default
    # implementation to perform this operation in one step without
    # allocating a temporary variable, without intermediate rounding, etc.

    int gr_submul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
    int gr_submul_ui(gr_ptr res, gr_srcptr x, unsigned long y, gr_ctx_t ctx)
    int gr_submul_si(gr_ptr res, gr_srcptr x, long y, gr_ctx_t ctx)
    int gr_submul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
    int gr_submul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
    int gr_submul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
    # Sets *res* to `\mathrm{res } - x \cdot y`.
    # Rings may override the default
    # implementation to perform this operation in one step without
    # allocating a temporary variable, without intermediate rounding, etc.

    int gr_mul_two(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    # Sets *res* to `2x`. The default implementation adds *x*
    # to itself.

    int gr_sqr(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    # Sets *res* to `x ^ 2`. The default implementation multiplies *x*
    # with itself.

    int gr_mul_2exp_si(gr_ptr res, gr_srcptr x, long y, gr_ctx_t ctx)
    int gr_mul_2exp_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
    # Sets *res* to `x \cdot 2^y`. This may perform `x \cdot 2^{-y}`
    # when *y* is negative, allowing exact division by powers of two
    # even if `2^{y}` is not representable.

    int gr_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
    int gr_div_ui(gr_ptr res, gr_srcptr x, unsigned long y, gr_ctx_t ctx)
    int gr_div_si(gr_ptr res, gr_srcptr x, long y, gr_ctx_t ctx)
    int gr_div_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
    int gr_div_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
    int gr_div_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
    int gr_other_div(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)
    # Sets *res* to the quotient `x / y` if such an element exists
    # in the present ring. Returns the flag ``GR_DOMAIN`` if no such
    # quotient exists.
    # Returns the flag ``GR_UNABLE`` if the implementation is unable
    # to perform the computation.
    # When the ring is not a field, the definition of division may
    # vary depending on the ring. A ring implementation may define
    # `x / y = x y^{-1}` and return ``GR_DOMAIN`` when `y^{-1}` does not
    # exist; alternatively, it may attempt to solve the equation
    # `q y = x` (which, for example, gives the usual exact
    # division in `\mathbb{Z}`).

    truth_t gr_is_invertible(gr_srcptr x, gr_ctx_t ctx)
    # Returns whether *x* has a multiplicative inverse in the present ring,
    # i.e. whether *x* is a unit.

    int gr_inv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    # Sets *res* to the multiplicative inverse of *x* in the present ring,
    # if such an element exists.
    # Returns the flag ``GR_DOMAIN`` if *x* is not invertible, or
    # ``GR_UNABLE`` if the implementation is unable to perform
    # the computation.

    truth_t gr_divides(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
    # Returns whether *x* divides *y*.

    int gr_divexact(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
    int gr_divexact_ui(gr_ptr res, gr_srcptr x, unsigned long y, gr_ctx_t ctx)
    int gr_divexact_si(gr_ptr res, gr_srcptr x, long y, gr_ctx_t ctx)
    int gr_divexact_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
    int gr_divexact_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
    int gr_other_divexact(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)
    # Sets *res* to the quotient `x / y`, assuming that this quotient
    # is exact in the present ring.
    # Rings may optimize this operation by not verifying that the
    # division is possible. If the division is not actually exact, the
    # implementation may set *res* to a nonsense value and still
    # return the ``GR_SUCCESS`` flag.

    int gr_euclidean_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
    int gr_euclidean_rem(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
    int gr_euclidean_divrem(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
    # In a Euclidean ring, these functions perform some version of Euclidean
    # division with remainder, where the choice of quotient is
    # implementation-defined. For example, it is standard to use
    # the round-to-floor quotient in `\mathbb{Z}` and a round-to-nearest quotient in `\mathbb{Z}[i]`.
    # In non-Euclidean rings, these functions may implement some generalization of
    # Euclidean division with remainder.

    int gr_pow(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
    int gr_pow_ui(gr_ptr res, gr_srcptr x, unsigned long y, gr_ctx_t ctx)
    int gr_pow_si(gr_ptr res, gr_srcptr x, long y, gr_ctx_t ctx)
    int gr_pow_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
    int gr_pow_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
    int gr_pow_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
    int gr_other_pow(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)
    # Sets *res* to the power `x ^ y`, the interpretation of which
    # depends on the ring when `y \not \in \mathbb{Z}`.
    # Returns the flag ``GR_DOMAIN`` if this power cannot be assigned
    # a meaningful value in the present ring, or ``GR_UNABLE`` if
    # the implementation is unable to perform the computation.
    # For subrings of `\mathbb{C}`, it is implied that the principal
    # power `x^y = \exp(y \log(x))` is computed for `x \ne 0`.
    # Default implementations of the powering methods support raising
    # elements to integer powers using a generic implementation of
    # exponentiation by squaring. Particular rings
    # should override these methods with faster versions or
    # to support more general notions of exponentiation when possible.

    truth_t gr_is_square(gr_srcptr x, gr_ctx_t ctx)
    # Returns whether *x* is a perfect square in the present ring.

    int gr_sqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    int gr_rsqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    # Sets *res* to a square root of *x* (respectively reciprocal
    # square root) in the present ring, if such an element exists.
    # Returns the flag ``GR_DOMAIN`` if *x* is not a perfect square
    # (also for zero, when computing the reciprocal square root), or
    # ``GR_UNABLE`` if the implementation is unable to perform
    # the computation.

    int gr_gcd(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
    # Sets *res* to a greatest common divisor (GCD) of *x* and *y*.
    # Since the GCD is unique only up to multiplication by a unit,
    # an implementation-defined representative is chosen.

    int gr_lcm(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
    # Sets *res* to a least common multiple (LCM) of *x* and *y*.
    # Since the LCM is unique only up to multiplication by a unit,
    # an implementation-defined representative is chosen.

    int gr_factor(gr_ptr c, gr_vec_t factors, gr_vec_t exponents, gr_srcptr x, int flags, gr_ctx_t ctx)
    # Given `x \in R`, computes a factorization
    # `x = c {f_1}^{e_1} \ldots {f_n}^{e_n}`
    # where `f_k` will be irreducible or prime (depending on `R`).
    # The prefactor `c` stores a unit, sign, or coefficient, e.g.\ the
    # sign `-1`, `0` or `+1` in `\mathbb{Z}`, or a sign multiplied
    # by the coefficient content in `\mathbb{Z}[x]`.
    # Note that this function outputs `c` as an element of the
    # same ring as the input: for example, in `\mathbb{Z}[x]`,
    # `c` will be a constant polynomial rather than an
    # element of the coefficient ring.
    # The exponents `e_k` are output as a vector of ``fmpz`` elements.
    # The factors `f_k` are guaranteed to be distinct,
    # but they are not guaranteed to be sorted in any particular
    # order.

    int gr_numerator(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    int gr_denominator(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    # Return a numerator `p` and denominator `q` such that `x = p/q`.
    # For typical fraction fields, the denominator will be minimal
    # and canonical.
    # However, some rings may return an arbitrary denominator as long
    # as the numerator matches.
    # The default implementations simply return `p = x` and `q = 1`.

    int gr_floor(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    int gr_ceil(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    int gr_trunc(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    int gr_nint(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    # In the real and complex numbers, sets *res* to the integer closest
    # to *x*, respectively rounding towards minus infinity, plus infinity,
    # zero, or the nearest integer (with tie-to-even).

    int gr_abs(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    # Sets *res* to the absolute value of *x*, which maybe defined
    # both in complex rings and in any ordered ring.

    int gr_i(gr_ptr res, gr_ctx_t ctx)
    # Sets *res* to the imaginary unit.

    int gr_conj(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    int gr_re(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    int gr_im(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    int gr_sgn(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    int gr_csgn(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    int gr_arg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
    # These methods may return the flag ``GR_DOMAIN`` (or ``GR_UNABLE``)
    # when the ring is not a subring of the real or complex numbers.

    int gr_pos_inf(gr_ptr res, gr_ctx_t ctx)
    int gr_neg_inf(gr_ptr res, gr_ctx_t ctx)
    int gr_uinf(gr_ptr res, gr_ctx_t ctx)
    int gr_undefined(gr_ptr res, gr_ctx_t ctx)
    int gr_unknown(gr_ptr res, gr_ctx_t ctx)

    int gr_cmp(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
    int gr_cmp_other(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
    # Sets *res* to -1, 0 or 1 according to whether *x* is less than,
    # equal or greater than the absolute value of *y*.
    # This may return ``GR_DOMAIN`` if the ring is not an ordered ring.

    int gr_cmpabs(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
    int gr_cmpabs_other(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
    # Sets *res* to -1, 0 or 1 according to whether the absolute value
    # of *x* is less than, equal or greater than the absolute value of *y*.
    # This may return ``GR_DOMAIN`` if the ring is not an ordered ring.

    int gr_ctx_fq_prime(fmpz_t p, gr_ctx_t ctx)

    int gr_ctx_fq_degree(long * deg, gr_ctx_t ctx)

    int gr_ctx_fq_order(fmpz_t q, gr_ctx_t ctx)

    int gr_fq_frobenius(gr_ptr res, gr_srcptr x, long e, gr_ctx_t ctx)

    int gr_fq_multiplicative_order(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)

    int gr_fq_norm(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)

    int gr_fq_trace(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)

    truth_t gr_fq_is_primitive(gr_srcptr x, gr_ctx_t ctx)

    int gr_fq_pth_root(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
