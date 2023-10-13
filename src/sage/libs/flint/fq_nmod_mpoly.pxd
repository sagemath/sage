# distutils: libraries = flint
# distutils: depends = flint/fq_nmod_mpoly.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void fq_nmod_mpoly_ctx_init(fq_nmod_mpoly_ctx_t ctx, long nvars, const ordering_t ord, const fq_nmod_ctx_t fqctx)
    # Initialise a context object for a polynomial ring with the given number of variables and the given ordering.
    # It will have coefficients in the finite field *fqctx*.
    # The possibilities for the ordering are ``ORD_LEX``, ``ORD_DEGLEX`` and ``ORD_DEGREVLEX``.

    long fq_nmod_mpoly_ctx_nvars(const fq_nmod_mpoly_ctx_t ctx)
    # Return the number of variables used to initialize the context.

    ordering_t fq_nmod_mpoly_ctx_ord(const fq_nmod_mpoly_ctx_t ctx)
    # Return the ordering used to initialize the context.

    void fq_nmod_mpoly_ctx_clear(fq_nmod_mpoly_ctx_t ctx)
    # Release any space allocated by an *ctx*.

    void fq_nmod_mpoly_init(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Initialise *A* for use with the given an initialised context object. Its value is set to zero.

    void fq_nmod_mpoly_init2(fq_nmod_mpoly_t A, long alloc, const fq_nmod_mpoly_ctx_t ctx)
    # Initialise *A* for use with the given an initialised context object. Its value is set to zero.
    # It is allocated with space for *alloc* terms and at least ``MPOLY_MIN_BITS`` bits for the exponents.

    void fq_nmod_mpoly_init3(fq_nmod_mpoly_t A, long alloc, flint_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx)
    # Initialise *A* for use with the given an initialised context object. Its value is set to zero.
    # It is allocated with space for *alloc* terms and *bits* bits for the exponents.

    void fq_nmod_mpoly_fit_length(fq_nmod_mpoly_t A, long len, const fq_nmod_mpoly_ctx_t ctx)
    # Ensure that *A* has space for at least *len* terms.

    void fq_nmod_mpoly_realloc(fq_nmod_mpoly_t A, long alloc, const fq_nmod_mpoly_ctx_t ctx)
    # Reallocate *A* to have space for *alloc* terms.
    # Assumes the current length of the polynomial is not greater than *alloc*.

    void fq_nmod_mpoly_clear(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Release any space allocated for *A*.

    char * fq_nmod_mpoly_get_str_pretty(const fq_nmod_mpoly_t A, const char ** x, const fq_nmod_mpoly_ctx_t ctx)
    # Return a string, which the user is responsible for cleaning up, representing *A*, given an array of variable strings *x*.

    int fq_nmod_mpoly_fprint_pretty(FILE * file, const fq_nmod_mpoly_t A, const char ** x, const fq_nmod_mpoly_ctx_t ctx)
    # Print a string representing *A* to *file*.

    int fq_nmod_mpoly_print_pretty(const fq_nmod_mpoly_t A, const char ** x, const fq_nmod_mpoly_ctx_t ctx)
    # Print a string representing *A* to ``stdout``.

    int fq_nmod_mpoly_set_str_pretty(fq_nmod_mpoly_t A, const char * str, const char ** x, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to the polynomial in the null-terminates string *str* given an array *x* of variable strings.
    # If parsing *str* fails, *A* is set to zero, and `-1` is returned. Otherwise, `0`  is returned.
    # The operations ``+``, ``-``, ``*``, and ``/`` are permitted along with integers and the variables in *x*. The character ``^`` must be immediately followed by the (integer) exponent.
    # If any division is not exact, parsing fails.

    void fq_nmod_mpoly_gen(fq_nmod_mpoly_t A, long var, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to the variable of index *var*, where `var = 0` corresponds to the variable with the most significance with respect to the ordering.

    int fq_nmod_mpoly_is_gen(const fq_nmod_mpoly_t A, long var, const fq_nmod_mpoly_ctx_t ctx)
    # If `var \ge 0`, return `1` if *A* is equal to the `var`-th generator, otherwise return `0`.
    # If `var < 0`, return `1` if the polynomial is equal to any generator, otherwise return `0`.

    void fq_nmod_mpoly_set(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to *B*.

    int fq_nmod_mpoly_equal(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
    # Return `1` if *A* is equal to *B*, else return `0`.

    void fq_nmod_mpoly_swap(fq_nmod_mpoly_t A, fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
    # Efficiently swap *A* and *B*.

    int fq_nmod_mpoly_is_fq_nmod(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Return `1` if *A* is a constant, else return `0`.

    void fq_nmod_mpoly_get_fq_nmod(fq_nmod_t c, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Assuming that *A* is a constant, set *c* to this constant.
    # This function throws if *A* is not a constant.

    void fq_nmod_mpoly_set_fq_nmod(fq_nmod_mpoly_t A, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx)
    void fq_nmod_mpoly_set_ui(fq_nmod_mpoly_t A, unsigned long c, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to the constant *c*.

    void fq_nmod_mpoly_set_fq_nmod_gen(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to the constant given by :func:`fq_nmod_gen`.

    void fq_nmod_mpoly_zero(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to the constant `0`.

    void fq_nmod_mpoly_one(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to the constant `1`.

    int fq_nmod_mpoly_equal_fq_nmod(const fq_nmod_mpoly_t A, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx)
    # Return `1` if *A* is equal to the constant *c*, else return `0`.

    int fq_nmod_mpoly_is_zero(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Return `1` if *A* is the constant `0`, else return `0`.

    int fq_nmod_mpoly_is_one(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Return `1` if *A* is the constant `1`, else return `0`.

    int fq_nmod_mpoly_degrees_fit_si(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Return `1` if the degrees of *A* with respect to each variable fit into an ``slong``, otherwise return `0`.

    void fq_nmod_mpoly_degrees_fmpz(fmpz ** degs, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    void fq_nmod_mpoly_degrees_si(long * degs, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Set *degs* to the degrees of *A* with respect to each variable.
    # If *A* is zero, all degrees are set to `-1`.

    void fq_nmod_mpoly_degree_fmpz(fmpz_t deg, const fq_nmod_mpoly_t A, long var, const fq_nmod_mpoly_ctx_t ctx)
    long fq_nmod_mpoly_degree_si(const fq_nmod_mpoly_t A, long var, const fq_nmod_mpoly_ctx_t ctx)
    # Either return or set *deg* to the degree of *A* with respect to the variable of index *var*.
    # If *A* is zero, the degree is defined to be `-1`.

    int fq_nmod_mpoly_total_degree_fits_si(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Return `1` if the total degree of *A* fits into an ``slong``, otherwise return `0`.

    void fq_nmod_mpoly_total_degree_fmpz(fmpz_t tdeg, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    long fq_nmod_mpoly_total_degree_si(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Either return or set *tdeg* to the total degree of *A*.
    # If *A* is zero, the total degree is defined to be `-1`.

    void fq_nmod_mpoly_used_vars(int * used, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # For each variable index `i`, set ``used[i]`` to nonzero if the variable of index `i` appears in *A* and to zero otherwise.

    void fq_nmod_mpoly_get_coeff_fq_nmod_monomial(fq_nmod_t c, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t M, const fq_nmod_mpoly_ctx_t ctx)
    # Assuming that *M* is a monomial, set *c* to the coefficient of the corresponding monomial in *A*.
    # This function throws if *M* is not a monomial.

    void fq_nmod_mpoly_set_coeff_fq_nmod_monomial(fq_nmod_mpoly_t A, const fq_nmod_t c, const fq_nmod_mpoly_t M, const fq_nmod_mpoly_ctx_t ctx)
    # Assuming that *M* is a monomial, set the coefficient of the corresponding monomial in *A* to *c*.
    # This function throws if *M* is not a monomial.

    void fq_nmod_mpoly_get_coeff_fq_nmod_fmpz(fq_nmod_t c, const fq_nmod_mpoly_t A, fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx)
    void fq_nmod_mpoly_get_coeff_fq_nmod_ui(fq_nmod_t c, const fq_nmod_mpoly_t A, const unsigned long * exp, const fq_nmod_mpoly_ctx_t ctx)
    # Set *c* to the coefficient of the monomial with exponent vector *exp*.

    void fq_nmod_mpoly_set_coeff_fq_nmod_fmpz(fq_nmod_mpoly_t A, const fq_nmod_t c, fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx)
    void fq_nmod_mpoly_set_coeff_fq_nmod_ui(fq_nmod_mpoly_t A, const fq_nmod_t c, const unsigned long * exp, const fq_nmod_mpoly_ctx_t ctx)
    # Set the coefficient of the monomial with exponent *exp* to *c*.

    void fq_nmod_mpoly_get_coeff_vars_ui(fq_nmod_mpoly_t C, const fq_nmod_mpoly_t A, const long * vars, const unsigned long * exps, long length, const fq_nmod_mpoly_ctx_t ctx)
    # Set *C* to the coefficient of *A* with respect to the variables in *vars* with powers in the corresponding array *exps*.
    # Both *vars* and *exps* point to array of length *length*. It is assumed that `0 < length \le nvars(A)` and that the variables in *vars* are distinct.

    int fq_nmod_mpoly_cmp(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
    # Return `1` (resp. `-1`, or `0`) if *A* is after (resp. before, same as) *B* in some arbitrary but fixed total ordering of the polynomials.
    # This ordering agrees with the usual ordering of monomials when *A* and *B* are both monomials.

    int fq_nmod_mpoly_is_canonical(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Return `1` if *A* is in canonical form. Otherwise, return `0`.
    # To be in canonical form, all of the terms must have nonzero coefficients, and the terms must be sorted from greatest to least.

    long fq_nmod_mpoly_length(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Return the number of terms in *A*.
    # If the polynomial is in canonical form, this will be the number of nonzero coefficients.

    void fq_nmod_mpoly_resize(fq_nmod_mpoly_t A, long new_length, const fq_nmod_mpoly_ctx_t ctx)
    # Set the length of *A* to ``new_length``.
    # Terms are either deleted from the end, or new zero terms are appended.

    void fq_nmod_mpoly_get_term_coeff_fq_nmod(fq_nmod_t c, const fq_nmod_mpoly_t A, long i, const fq_nmod_mpoly_ctx_t ctx)
    # Set *c* to the coefficient of the term of index *i*.

    void fq_nmod_mpoly_set_term_coeff_ui(fq_nmod_mpoly_t A, long i, unsigned long c, const fq_nmod_mpoly_ctx_t ctx)
    # Set the coefficient of the term of index *i* to *c*.

    int fq_nmod_mpoly_term_exp_fits_si(const fq_nmod_mpoly_t A, long i, const fq_nmod_mpoly_ctx_t ctx)
    int fq_nmod_mpoly_term_exp_fits_ui(const fq_nmod_mpoly_t A, long i, const fq_nmod_mpoly_ctx_t ctx)
    # Return `1` if all entries of the exponent vector of the term of index `i` fit into an ``slong`` (resp. a ``ulong``). Otherwise, return `0`.

    void fq_nmod_mpoly_get_term_exp_fmpz(fmpz ** exp, const fq_nmod_mpoly_t A, long i, const fq_nmod_mpoly_ctx_t ctx)
    void fq_nmod_mpoly_get_term_exp_ui(unsigned long * exp, const fq_nmod_mpoly_t A, long i, const fq_nmod_mpoly_ctx_t ctx)
    void fq_nmod_mpoly_get_term_exp_si(long * exp, const fq_nmod_mpoly_t A, long i, const fq_nmod_mpoly_ctx_t ctx)
    # Set *exp* to the exponent vector of the term of index *i*.
    # The ``_ui`` (resp. ``_si``) version throws if any entry does not fit into a ``ulong`` (resp. ``slong``).

    unsigned long fq_nmod_mpoly_get_term_var_exp_ui(const fq_nmod_mpoly_t A, long i, long var, const fq_nmod_mpoly_ctx_t ctx)
    long fq_nmod_mpoly_get_term_var_exp_si(const fq_nmod_mpoly_t A, long i, long var, const fq_nmod_mpoly_ctx_t ctx)
    # Return the exponent of the variable *var* of the term of index *i*.
    # This function throws if the exponent does not fit into a ``ulong`` (resp. ``slong``).

    void fq_nmod_mpoly_set_term_exp_fmpz(fq_nmod_mpoly_t A, long i, fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx)
    void fq_nmod_mpoly_set_term_exp_ui(fq_nmod_mpoly_t A, long i, const unsigned long * exp, const fq_nmod_mpoly_ctx_t ctx)
    # Set the exponent of the term of index *i* to *exp*.

    void fq_nmod_mpoly_get_term(fq_nmod_mpoly_t M, const fq_nmod_mpoly_t A, long i, const fq_nmod_mpoly_ctx_t ctx)
    # Set *M* to the term of index *i* in *A*.

    void fq_nmod_mpoly_get_term_monomial(fq_nmod_mpoly_t M, const fq_nmod_mpoly_t A, long i, const fq_nmod_mpoly_ctx_t ctx)
    # Set *M* to the monomial of the term of index *i* in *A*. The coefficient of *M* will be one.

    void fq_nmod_mpoly_push_term_fq_nmod_fmpz(fq_nmod_mpoly_t A, const fq_nmod_t c, fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx)
    void fq_nmod_mpoly_push_term_fq_nmod_ffmpz(fq_nmod_mpoly_t A, const fq_nmod_t c, const fmpz * exp, const fq_nmod_mpoly_ctx_t ctx)
    void fq_nmod_mpoly_push_term_fq_nmod_ui(fq_nmod_mpoly_t A, const fq_nmod_t c, const unsigned long * exp, const fq_nmod_mpoly_ctx_t ctx)
    # Append a term to *A* with coefficient *c* and exponent vector *exp*.
    # This function runs in constant average time.

    void fq_nmod_mpoly_sort_terms(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Sort the terms of *A* into the canonical ordering dictated by the ordering in *ctx*.
    # This function simply reorders the terms: It does not combine like terms, nor does it delete terms with coefficient zero.
    # This function runs in linear time in the bit size of *A*.

    void fq_nmod_mpoly_combine_like_terms(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Combine adjacent like terms in *A* and delete terms with coefficient zero.
    # If the terms of *A* were sorted to begin with, the result will be in canonical form.
    # This function runs in linear time in the bit size of *A*.

    void fq_nmod_mpoly_reverse(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to the reversal of *B*.

    void fq_nmod_mpoly_randtest_bound(fq_nmod_mpoly_t A, flint_rand_t state, long length, unsigned long exp_bound, const fq_nmod_mpoly_ctx_t ctx)
    # Generate a random polynomial with length up to *length* and exponents in the range ``[0, exp_bound - 1]``.
    # The exponents of each variable are generated by calls to  ``n_randint(state, exp_bound)``.

    void fq_nmod_mpoly_randtest_bounds(fq_nmod_mpoly_t A, flint_rand_t state, long length, unsigned long *exp_bounds, const fq_nmod_mpoly_ctx_t ctx)
    # Generate a random polynomial with length up to *length* and exponents in the range ``[0, exp_bounds[i] - 1]``.
    # The exponents of the variable of index *i* are generated by calls to ``n_randint(state, exp_bounds[i])``.

    void fq_nmod_mpoly_randtest_bits(fq_nmod_mpoly_t A, flint_rand_t state, long length, mp_limb_t exp_bits, const fq_nmod_mpoly_ctx_t ctx)
    # Generate a random polynomial with length up to *length* and exponents whose packed form does not exceed the given bit count.

    void fq_nmod_mpoly_add_fq_nmod(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_t C, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to `B + c`.

    void fq_nmod_mpoly_sub_fq_nmod(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_t C, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to `B - c`.

    void fq_nmod_mpoly_add(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_t C, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to `B + C`.

    void fq_nmod_mpoly_sub(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_t C, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to `B - C`.

    void fq_nmod_mpoly_neg(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to `-B`.

    void fq_nmod_mpoly_scalar_mul_fq_nmod(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to `B \times c`.

    void fq_nmod_mpoly_make_monic(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to *B* divided by the leading coefficient of *B*.
    # This throws if *B* is zero.

    void fq_nmod_mpoly_derivative(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, long var, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to the derivative of *B* with respect to the variable of index *var*.

    void fq_nmod_mpoly_evaluate_all_fq_nmod(fq_nmod_t ev, const fq_nmod_mpoly_t A, fq_nmod_struct * const *  vals, const fq_nmod_mpoly_ctx_t ctx)
    # Set *ev* the evaluation of *A* where the variables are replaced by the corresponding elements of the array *vals*.

    void fq_nmod_mpoly_evaluate_one_fq_nmod(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, long var, const fq_nmod_t val, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to the evaluation of *B* where the variable of index *var* is replaced by *val*.

    int fq_nmod_mpoly_compose_fq_nmod_poly(fq_nmod_poly_t A, const fq_nmod_mpoly_t B, fq_nmod_poly_struct * const * C, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to the evaluation of *B* where the variables are replaced by the corresponding elements of the array *C*.
    # The context object of *B* is *ctxB*.
    # Return `1` for success and `0` for failure.

    int fq_nmod_mpoly_compose_fq_nmod_mpoly(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, fq_nmod_mpoly_struct * const * C, const fq_nmod_mpoly_ctx_t ctxB, const fq_nmod_mpoly_ctx_t ctxAC)
    # Set *A* to the evaluation of *B* where the variables are replaced by the corresponding elements of the array *C*.
    # Both *A* and the elements of *C* have context object *ctxAC*, while *B* has context object *ctxB*.
    # Neither *A* nor *B* is allowed to alias any other polynomial.
    # Return `1` for success and `0` for failure.

    void fq_nmod_mpoly_compose_fq_nmod_mpoly_gen(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const long * c, const fq_nmod_mpoly_ctx_t ctxB, const fq_nmod_mpoly_ctx_t ctxAC)
    # Set *A* to the evaluation of *B* where the variable of index *i* in *ctxB* is replaced by the variable of index ``c[i]`` in *ctxAC*.
    # The length of the array *C* is the number of variables in *ctxB*.
    # If any ``c[i]`` is negative, the corresponding variable of *B* is replaced by zero. Otherwise, it is expected that ``c[i]`` is less than the number of variables in *ctxAC*.

    void fq_nmod_mpoly_mul(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_t C, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to *B* times *C*.

    int fq_nmod_mpoly_pow_fmpz(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fmpz_t k, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to `B` raised to the *k*-th power.
    # Return `1` for success and `0` for failure.

    int fq_nmod_mpoly_pow_ui(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, unsigned long k, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to `B` raised to the *k*-th power.
    # Return `1` for success and `0` for failure.

    int fq_nmod_mpoly_divides(fq_nmod_mpoly_t Q, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
    # If *A* is divisible by *B*, set *Q* to the exact quotient and return `1`. Otherwise, set *Q* to zero and return `0`.

    void fq_nmod_mpoly_div(fq_nmod_mpoly_t Q, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
    # Set *Q* to the quotient of *A* by *B*, discarding the remainder.

    void fq_nmod_mpoly_divrem(fq_nmod_mpoly_t Q, fq_nmod_mpoly_t R, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
    # Set *Q* and *R* to the quotient and remainder of *A* divided by *B*.

    void fq_nmod_mpoly_divrem_ideal(fq_nmod_mpoly_struct ** Q, fq_nmod_mpoly_t R, const fq_nmod_mpoly_t A, fq_nmod_mpoly_struct * const * B, long len, const fq_nmod_mpoly_ctx_t ctx)
    # This function is as per :func:`fq_nmod_mpoly_divrem` except that it takes an array of divisor polynomials *B* and it returns an array of quotient polynomials *Q*.
    # The number of divisor (and hence quotient) polynomials, is given by *len*.

    void fq_nmod_mpoly_term_content(fq_nmod_mpoly_t M, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Set *M* to the GCD of the terms of *A*.
    # If *A* is zero, *M* will be zero. Otherwise, *M* will be a monomial with coefficient one.

    int fq_nmod_mpoly_content_vars(fq_nmod_mpoly_t g, const fq_nmod_mpoly_t A, long * vars, long vars_length, const fq_nmod_mpoly_ctx_t ctx)
    # Set *g* to the GCD of the coefficients of *A* when viewed as a polynomial in the variables *vars*.
    # Return `1` for success and `0` for failure. Upon success, *g* will be independent of the variables *vars*.

    int fq_nmod_mpoly_gcd(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
    # Try to set *G* to the monic GCD of *A* and *B*. The GCD of zero and zero is defined to be zero.
    # If the return is `1` the function was successful. Otherwise the return is  `0` and *G* is left untouched.

    int fq_nmod_mpoly_gcd_cofactors(fq_nmod_mpoly_t G, fq_nmod_mpoly_t Abar, fq_nmod_mpoly_t Bbar, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
    # Do the operation of :func:`fq_nmod_mpoly_gcd` and also compute `Abar = A/G` and `Bbar = B/G` if successful.

    int fq_nmod_mpoly_gcd_brown(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
    int fq_nmod_mpoly_gcd_hensel(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
    int fq_nmod_mpoly_gcd_zippel(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
    # Try to set *G* to the GCD of *A* and *B* using various algorithms.

    int fq_nmod_mpoly_resultant(fq_nmod_mpoly_t R, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, long var, const fq_nmod_mpoly_ctx_t ctx)
    # Try to set *R* to the resultant of *A* and *B* with respect to the variable of index *var*.

    int fq_nmod_mpoly_discriminant(fq_nmod_mpoly_t D, const fq_nmod_mpoly_t A, long var, const fq_nmod_mpoly_ctx_t ctx)
    # Try to set *D* to the discriminant of *A* with respect to the variable of index *var*.

    int fq_nmod_mpoly_sqrt(fq_nmod_mpoly_t Q, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # If `Q^2=A` has a solution, set `Q` to a solution and return `1`, otherwise return `0` and set `Q` to zero.

    int fq_nmod_mpoly_is_square(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Return `1` if *A* is a perfect square, otherwise return `0`.

    int fq_nmod_mpoly_quadratic_root(fq_nmod_mpoly_t Q, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
    # If `Q^2+AQ=B` has a solution, set `Q` to a solution and return `1`, otherwise return `0`.

    void fq_nmod_mpoly_univar_init(fq_nmod_mpoly_univar_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Initialize *A*.

    void fq_nmod_mpoly_univar_clear(fq_nmod_mpoly_univar_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Clear *A*.

    void fq_nmod_mpoly_univar_swap(fq_nmod_mpoly_univar_t A, fq_nmod_mpoly_univar_t B, const fq_nmod_mpoly_ctx_t ctx)
    # Swap *A* and `B`.

    void fq_nmod_mpoly_to_univar(fq_nmod_mpoly_univar_t A, const fq_nmod_mpoly_t B, long var, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to a univariate form of *B* by pulling out the variable of index *var*.
    # The coefficients of *A* will still belong to the content *ctx* but will not depend on the variable of index *var*.

    void fq_nmod_mpoly_from_univar(fq_nmod_mpoly_t A, const fq_nmod_mpoly_univar_t B, long var, const fq_nmod_mpoly_ctx_t ctx)
    # Set *A* to the normal form of *B* by putting in the variable of index *var*.
    # This function is undefined if the coefficients of *B* depend on the variable of index *var*.

    int fq_nmod_mpoly_univar_degree_fits_si(const fq_nmod_mpoly_univar_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Return `1` if the degree of *A* with respect to the main variable fits an ``slong``. Otherwise, return `0`.

    long fq_nmod_mpoly_univar_length(const fq_nmod_mpoly_univar_t A, const fq_nmod_mpoly_ctx_t ctx)
    # Return the number of terms in *A* with respect to the main variable.

    long fq_nmod_mpoly_univar_get_term_exp_si(fq_nmod_mpoly_univar_t A, long i, const fq_nmod_mpoly_ctx_t ctx)
    # Return the exponent of the term of index *i* of *A*.

    void fq_nmod_mpoly_univar_get_term_coeff(fq_nmod_mpoly_t c, const fq_nmod_mpoly_univar_t A, long i, const fq_nmod_mpoly_ctx_t ctx)
    void fq_nmod_mpoly_univar_swap_term_coeff(fq_nmod_mpoly_t c, fq_nmod_mpoly_univar_t A, long i, const fq_nmod_mpoly_ctx_t ctx)
    # Set (resp. swap) *c* to (resp. with) the coefficient of the term of index *i* of *A*.
