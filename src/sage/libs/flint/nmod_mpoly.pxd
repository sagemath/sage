# distutils: libraries = flint
# distutils: depends = flint/nmod_mpoly.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void nmod_mpoly_ctx_init(nmod_mpoly_ctx_t ctx, slong nvars, const ordering_t ord, mp_limb_t n)
    # Initialise a context object for a polynomial ring with the given number of variables and the given ordering.
    # It will have coefficients modulo *n*. Setting `n = 0` will give undefined behavior.
    # The possibilities for the ordering are ``ORD_LEX``, ``ORD_DEGLEX`` and ``ORD_DEGREVLEX``.

    slong nmod_mpoly_ctx_nvars(const nmod_mpoly_ctx_t ctx)
    # Return the number of variables used to initialize the context.

    ordering_t nmod_mpoly_ctx_ord(const nmod_mpoly_ctx_t ctx)
    # Return the ordering used to initialize the context.

    mp_limb_t nmod_mpoly_ctx_modulus(const nmod_mpoly_ctx_t ctx)
    # Return the modulus used to initialize the context.

    void nmod_mpoly_ctx_clear(nmod_mpoly_ctx_t ctx)
    # Release any space allocated by *ctx*.

    void nmod_mpoly_init(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Initialise *A* for use with the given an initialised context object. Its value is set to zero.

    void nmod_mpoly_init2(nmod_mpoly_t A, slong alloc, const nmod_mpoly_ctx_t ctx)
    # Initialise *A* for use with the given an initialised context object. Its value is set to zero.
    # It is allocated with space for *alloc* terms and at least ``MPOLY_MIN_BITS`` bits for the exponent widths.

    void nmod_mpoly_init3(nmod_mpoly_t A, slong alloc, flint_bitcnt_t bits, const nmod_mpoly_ctx_t ctx)
    # Initialise *A* for use with the given an initialised context object. Its value is set to zero.
    # It is allocated with space for *alloc* terms and *bits* bits for the exponents.

    void nmod_mpoly_fit_length(nmod_mpoly_t A, slong len, const nmod_mpoly_ctx_t ctx)
    # Ensure that *A* has space for at least *len* terms.

    void nmod_mpoly_realloc(nmod_mpoly_t A, slong alloc, const nmod_mpoly_ctx_t ctx)
    # Reallocate *A* to have space for *alloc* terms.
    # Assumes the current length of the polynomial is not greater than *alloc*.

    void nmod_mpoly_clear(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Release any space allocated for *A*.

    char * nmod_mpoly_get_str_pretty(const nmod_mpoly_t A, const char ** x, const nmod_mpoly_ctx_t ctx)
    # Return a string, which the user is responsible for cleaning up, representing *A*, given an array of variable strings *x*.

    int nmod_mpoly_fprint_pretty(FILE * file, const nmod_mpoly_t A, const char ** x, const nmod_mpoly_ctx_t ctx)
    # Print a string representing *A* to *file*.

    int nmod_mpoly_print_pretty(const nmod_mpoly_t A, const char ** x, const nmod_mpoly_ctx_t ctx)
    # Print a string representing *A* to ``stdout``.

    int nmod_mpoly_set_str_pretty(nmod_mpoly_t A, const char * str, const char ** x, const nmod_mpoly_ctx_t ctx)
    # Set *A* to the polynomial in the null-terminates string *str* given an array *x* of variable strings.
    # If parsing *str* fails, *A* is set to zero, and `-1` is returned. Otherwise, `0`  is returned.
    # The operations ``+``, ``-``, ``*``, and ``/`` are permitted along with integers and the variables in *x*. The character ``^`` must be immediately followed by the (integer) exponent.
    # If any division is not exact, parsing fails.

    void nmod_mpoly_gen(nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)
    # Set *A* to the variable of index *var*, where `var = 0` corresponds to the variable with the most significance with respect to the ordering.

    bint nmod_mpoly_is_gen(const nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)
    # If `var \ge 0`, return `1` if *A* is equal to the `var`-th generator, otherwise return `0`.
    # If `var < 0`, return `1` if the polynomial is equal to any generator, otherwise return `0`.

    void nmod_mpoly_set(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    # Set *A* to *B*.

    bint nmod_mpoly_equal(const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    # Return `1` if *A* is equal to *B*, else return `0`.

    void nmod_mpoly_swap(nmod_mpoly_t A, nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    # Efficiently swap *A* and *B*.

    bint nmod_mpoly_is_ui(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Return `1` if *A* is a constant, else return `0`.

    ulong nmod_mpoly_get_ui(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Assuming that *A* is a constant, return this constant.
    # This function throws if *A* is not a constant.

    void nmod_mpoly_set_ui(nmod_mpoly_t A, ulong c, const nmod_mpoly_ctx_t ctx)
    # Set *A* to the constant *c*.

    void nmod_mpoly_zero(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Set *A* to the constant `0`.

    void nmod_mpoly_one(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Set *A* to the constant `1`.

    bint nmod_mpoly_equal_ui(const nmod_mpoly_t A, ulong c, const nmod_mpoly_ctx_t ctx)
    # Return `1` if *A* is equal to the constant *c*, else return `0`.

    bint nmod_mpoly_is_zero(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Return `1` if *A* is the constant `0`, else return `0`.

    bint nmod_mpoly_is_one(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Return `1` if *A* is the constant `1`, else return `0`.

    int nmod_mpoly_degrees_fit_si(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Return `1` if the degrees of *A* with respect to each variable fit into an ``slong``, otherwise return `0`.

    void nmod_mpoly_degrees_fmpz(fmpz ** degs, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_degrees_si(slong * degs, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Set *degs* to the degrees of *A* with respect to each variable.
    # If *A* is zero, all degrees are set to `-1`.

    void nmod_mpoly_degree_fmpz(fmpz_t deg, const nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)
    slong nmod_mpoly_degree_si(const nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)
    # Either return or set *deg* to the degree of *A* with respect to the variable of index *var*.
    # If *A* is zero, the degree is defined to be `-1`.

    int nmod_mpoly_total_degree_fits_si(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Return `1` if the total degree of *A* fits into an ``slong``, otherwise return `0`.

    void nmod_mpoly_total_degree_fmpz(fmpz_t tdeg, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    slong nmod_mpoly_total_degree_si(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Either return or set *tdeg* to the total degree of *A*.
    # If *A* is zero, the total degree is defined to be `-1`.

    void nmod_mpoly_used_vars(int * used, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # For each variable index *i*, set ``used[i]`` to nonzero if the variable of index *i* appears in *A* and to zero otherwise.

    ulong nmod_mpoly_get_coeff_ui_monomial(const nmod_mpoly_t A, const nmod_mpoly_t M, const nmod_mpoly_ctx_t ctx)
    # Assuming that *M* is a monomial, return the coefficient of the corresponding monomial in *A*.
    # This function throws if *M* is not a monomial.

    void nmod_mpoly_set_coeff_ui_monomial(nmod_mpoly_t A, ulong c, const nmod_mpoly_t M, const nmod_mpoly_ctx_t ctx)
    # Assuming that *M* is a monomial, set the coefficient of the corresponding monomial in *A* to *c*.
    # This function throws if *M* is not a monomial.

    ulong nmod_mpoly_get_coeff_ui_fmpz(const nmod_mpoly_t A, fmpz * const * exp, const nmod_mpoly_ctx_t ctx)
    ulong nmod_mpoly_get_coeff_ui_ui(const nmod_mpoly_t A, const ulong * exp, const nmod_mpoly_ctx_t ctx)
    # Return the coefficient of the monomial with exponent *exp*.

    void nmod_mpoly_set_coeff_ui_fmpz(nmod_mpoly_t A, ulong c, fmpz * const * exp, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_set_coeff_ui_ui(nmod_mpoly_t A, ulong c, const ulong * exp, const nmod_mpoly_ctx_t ctx)
    # Set the coefficient of the monomial with exponent *exp* to `c`.

    void nmod_mpoly_get_coeff_vars_ui(nmod_mpoly_t C, const nmod_mpoly_t A, const slong * vars, const ulong * exps, slong length, const nmod_mpoly_ctx_t ctx)
    # Set *C* to the coefficient of *A* with respect to the variables in *vars* with powers in the corresponding array *exps*.
    # Both *vars* and *exps* point to array of length *length*. It is assumed that ``0 < length \le nvars(A)`` and that the variables in *vars* are distinct.

    int nmod_mpoly_cmp(const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    # Return `1` (resp. `-1`, or `0`) if *A* is after (resp. before, same as) *B* in some arbitrary but fixed total ordering of the polynomials.
    # This ordering agrees with the usual ordering of monomials when *A* and *B* are both monomials.

    mp_limb_t * nmod_mpoly_term_coeff_ref(nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
    # Return a reference to the coefficient of index *i* of *A*.

    bint nmod_mpoly_is_canonical(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Return `1` if *A* is in canonical form. Otherwise, return `0`.
    # To be in canonical form, all of the terms must have nonzero coefficients, and the terms must be sorted from greatest to least.

    slong nmod_mpoly_length(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Return the number of terms in *A*.
    # If the polynomial is in canonical form, this will be the number of nonzero coefficients.

    void nmod_mpoly_resize(nmod_mpoly_t A, slong new_length, const nmod_mpoly_ctx_t ctx)
    # Set the length of *A* to ``new_length``.
    # Terms are either deleted from the end, or new zero terms are appended.

    ulong nmod_mpoly_get_term_coeff_ui(const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
    # Return the coefficient of the term of index *i*.

    void nmod_mpoly_set_term_coeff_ui(nmod_mpoly_t A, slong i, ulong c, const nmod_mpoly_ctx_t ctx)
    # Set the coefficient of the term of index *i* to *c*.

    int nmod_mpoly_term_exp_fits_si(const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_term_exp_fits_ui(const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
    # Return `1` if all entries of the exponent vector of the term of index *i* fit into an ``slong`` (resp. a ``ulong``). Otherwise, return `0`.

    void nmod_mpoly_get_term_exp_fmpz(fmpz ** exp, const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_get_term_exp_ui(ulong * exp, const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)

    void nmod_mpoly_get_term_exp_si(slong * exp, const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
    # Set *exp* to the exponent vector of the term of index *i*.
    # The ``_ui`` (resp. ``_si``) version throws if any entry does not fit into a ``ulong`` (resp. ``slong``).

    ulong nmod_mpoly_get_term_var_exp_ui(const nmod_mpoly_t A, slong i, slong var, const nmod_mpoly_ctx_t ctx)
    slong nmod_mpoly_get_term_var_exp_si(const nmod_mpoly_t A, slong i, slong var, const nmod_mpoly_ctx_t ctx)
    # Return the exponent of the variable *var* of the term of index *i*.
    # This function throws if the exponent does not fit into a ``ulong`` (resp. ``slong``).

    void nmod_mpoly_set_term_exp_fmpz(nmod_mpoly_t A, slong i, fmpz * const * exp, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_set_term_exp_ui(nmod_mpoly_t A, slong i, const ulong * exp, const nmod_mpoly_ctx_t ctx)
    # Set the exponent of the term of index *i* to *exp*.

    void nmod_mpoly_get_term(nmod_mpoly_t M, const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
    # Set *M* to the term of index *i* in *A*.

    void nmod_mpoly_get_term_monomial(nmod_mpoly_t M, const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
    # Set *M* to the monomial of the term of index *i* in *A*. The coefficient of *M* will be one.

    void nmod_mpoly_push_term_ui_fmpz(nmod_mpoly_t A, ulong c, fmpz * const * exp, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_push_term_ui_ffmpz(nmod_mpoly_t A, ulong c, const fmpz * exp, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_push_term_ui_ui(nmod_mpoly_t A, ulong c, const ulong * exp, const nmod_mpoly_ctx_t ctx)
    # Append a term to *A* with coefficient *c* and exponent vector *exp*.
    # This function runs in constant average time.

    void nmod_mpoly_sort_terms(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Sort the terms of *A* into the canonical ordering dictated by the ordering in *ctx*.
    # This function simply reorders the terms: It does not combine like terms, nor does it delete terms with coefficient zero.
    # This function runs in linear time in the bit size of *A*.

    void nmod_mpoly_combine_like_terms(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Combine adjacent like terms in *A* and delete terms with coefficient zero.
    # If the terms of *A* were sorted to begin with, the result will be in canonical form.
    # This function runs in linear time in the bit size of *A*.

    void nmod_mpoly_reverse(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    # Set *A* to the reversal of *B*.

    void nmod_mpoly_randtest_bound(nmod_mpoly_t A, flint_rand_t state, slong length, ulong exp_bound, const nmod_mpoly_ctx_t ctx)
    # Generate a random polynomial with length up to *length* and exponents in the range ``[0, exp_bound - 1]``.
    # The exponents of each variable are generated by calls to  ``n_randint(state, exp_bound)``.

    void nmod_mpoly_randtest_bounds(nmod_mpoly_t A, flint_rand_t state, slong length, ulong * exp_bounds, const nmod_mpoly_ctx_t ctx)
    # Generate a random polynomial with length up to *length* and exponents in the range ``[0, exp_bounds[i] - 1]``.
    # The exponents of the variable of index *i* are generated by calls to ``n_randint(state, exp_bounds[i])``.

    void nmod_mpoly_randtest_bits(nmod_mpoly_t A, flint_rand_t state, slong length, mp_limb_t exp_bits, const nmod_mpoly_ctx_t ctx)
    # Generate a random polynomial with length up to *length* and exponents whose packed form does not exceed the given bit count.

    void nmod_mpoly_add_ui(nmod_mpoly_t A, const nmod_mpoly_t B, ulong c, const nmod_mpoly_ctx_t ctx)
    # Set *A* to `B + c`.

    void nmod_mpoly_sub_ui(nmod_mpoly_t A, const nmod_mpoly_t B, ulong c, const nmod_mpoly_ctx_t ctx)
    # Set *A* to `B - c`.

    void nmod_mpoly_add(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
    # Set *A* to `B + C`.

    void nmod_mpoly_sub(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
    # Set *A* to `B - C`.

    void nmod_mpoly_neg(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    # Set *A* to `-B`.

    void nmod_mpoly_scalar_mul_ui(nmod_mpoly_t A, const nmod_mpoly_t B, ulong c, const nmod_mpoly_ctx_t ctx)
    # Set *A* to `B \times c`.

    void nmod_mpoly_make_monic(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    # Set *A* to *B* divided by the leading coefficient of *B*.
    # This throws if *B* is zero or the leading coefficient is not invertible.

    void nmod_mpoly_derivative(nmod_mpoly_t A, const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx)
    # Set *A* to the derivative of *B* with respect to the variable of index *var*.

    ulong nmod_mpoly_evaluate_all_ui(const nmod_mpoly_t A, const ulong * vals, const nmod_mpoly_ctx_t ctx)
    # Return the evaluation of *A* where the variables are replaced by the corresponding elements of the array *vals*.

    void nmod_mpoly_evaluate_one_ui(nmod_mpoly_t A, const nmod_mpoly_t B, slong var, ulong val, const nmod_mpoly_ctx_t ctx)
    # Set *A* to the evaluation of *B* where the variable of index *var* is replaced by *val*.

    int nmod_mpoly_compose_nmod_poly(nmod_poly_t A, const nmod_mpoly_t B, nmod_poly_struct * const * C, const nmod_mpoly_ctx_t ctx)
    # Set *A* to the evaluation of *B* where the variables are replaced by the corresponding elements of the array *C*.
    # The context object of *B* is *ctxB*.
    # Return `1` for success and `0` for failure.

    int nmod_mpoly_compose_nmod_mpoly_geobucket(nmod_mpoly_t A, const nmod_mpoly_t B, nmod_mpoly_struct * const * C, const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC)
    int nmod_mpoly_compose_nmod_mpoly_horner(nmod_mpoly_t A, const nmod_mpoly_t B, nmod_mpoly_struct * const * C, const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC)
    int nmod_mpoly_compose_nmod_mpoly(nmod_mpoly_t A, const nmod_mpoly_t B, nmod_mpoly_struct * const * C, const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC)
    # Set *A* to the evaluation of *B* where the variables are replaced by the corresponding elements of the array *C*.
    # Both *A* and the elements of *C* have context object *ctxAC*, while *B* has context object *ctxB*.
    # Neither of *A* and *B* is allowed to alias any other polynomial.
    # Return `1` for success and `0` for failure.
    # The main method attempts to perform the calculation using matrices and chooses heuristically between the ``geobucket`` and ``horner`` methods if needed.

    void nmod_mpoly_compose_nmod_mpoly_gen(nmod_mpoly_t A, const nmod_mpoly_t B, const slong * c, const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC)
    # Set *A* to the evaluation of *B* where the variable of index *i* in *ctxB* is replaced by the variable of index ``c[i]`` in *ctxAC*.
    # The length of the array *C* is the number of variables in *ctxB*.
    # If any ``c[i]`` is negative, the corresponding variable of *B* is replaced by zero. Otherwise, it is expected that ``c[i]`` is less than the number of variables in *ctxAC*.

    void nmod_mpoly_mul(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
    # Set *A* to `B \times C`.

    void nmod_mpoly_mul_johnson(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_mul_heap_threaded(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
    # Set *A* to `B \times C` using Johnson's heap-based method.
    # The first version always uses one thread.

    int nmod_mpoly_mul_array(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_mul_array_threaded(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
    # Try to set *A* to `B \times C` using arrays.
    # If the return is `0`, the operation was unsuccessful. Otherwise, it was successful, and the return is `1`.
    # The first version always uses one thread.

    int nmod_mpoly_mul_dense(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
    # Try to set *A* to `B \times C` using univariate arithmetic.
    # If the return is `0`, the operation was unsuccessful. Otherwise, it was successful and the return is `1`.

    int nmod_mpoly_pow_fmpz(nmod_mpoly_t A, const nmod_mpoly_t B, const fmpz_t k, const nmod_mpoly_ctx_t ctx)
    # Set *A* to *B* raised to the *k*-th power.
    # Return `1` for success and `0` for failure.

    int nmod_mpoly_pow_ui(nmod_mpoly_t A, const nmod_mpoly_t B, ulong k, const nmod_mpoly_ctx_t ctx)
    # Set *A* to *B* raised to the *k*-th power.
    # Return `1` for success and `0` for failure.

    int nmod_mpoly_divides(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    # If *A* is divisible by *B*, set *Q* to the exact quotient and return `1`. Otherwise, set *Q* to zero and return `0`.
    # Note that the function ``nmod_mpoly_div`` below may be faster if the quotient is known to be exact.

    void nmod_mpoly_div(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    # Set *Q* to the quotient of *A* by *B*, discarding the remainder.

    void nmod_mpoly_divrem(nmod_mpoly_t Q, nmod_mpoly_t R, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    # Set *Q* and *R* to the quotient and remainder of *A* divided by *B*.

    void nmod_mpoly_divrem_ideal(nmod_mpoly_struct ** Q, nmod_mpoly_t R, const nmod_mpoly_t A, nmod_mpoly_struct * const * B, slong len, const nmod_mpoly_ctx_t ctx)
    # This function is as per :func:`nmod_mpoly_divrem` except that it takes an array of divisor polynomials *B* and it returns an array of quotient polynomials *Q*.
    # The number of divisor (and hence quotient) polynomials, is given by *len*.

    int nmod_mpoly_divides_dense(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    # Try to do the operation of ``nmod_mpoly_divides`` using univariate arithmetic.
    # If the return is `-1`, the operation was unsuccessful. Otherwise, it was successful and the return is `0` or `1`.

    int nmod_mpoly_divides_monagan_pearce(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    # Do the operation of ``nmod_mpoly_divides`` using the algorithm of Michael Monagan and Roman Pearce.

    int nmod_mpoly_divides_heap_threaded(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    # Do the operation of ``nmod_mpoly_divides`` using a heap and multiple threads.
    # This function should only be called once ``global_thread_pool`` has been initialized.

    void nmod_mpoly_term_content(nmod_mpoly_t M, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Set *M* to the GCD of the terms of *A*.
    # If *A* is zero, *M* will be zero. Otherwise, *M* will be a monomial with coefficient one.

    int nmod_mpoly_content_vars(nmod_mpoly_t g, const nmod_mpoly_t A, slong * vars, slong vars_length, const nmod_mpoly_ctx_t ctx)
    # Set *g* to the GCD of the coefficients of *A* when viewed as a polynomial in the variables *vars*.
    # Return `1` for success and `0` for failure. Upon success, *g* will be independent of the variables *vars*.

    int nmod_mpoly_gcd(nmod_mpoly_t G, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    # Try to set *G* to the monic GCD of *A* and *B*. The GCD of zero and zero is defined to be zero.
    # If the return is `1` the function was successful. Otherwise the return is  `0` and *G* is left untouched.

    int nmod_mpoly_gcd_cofactors(nmod_mpoly_t G, nmod_mpoly_t Abar, nmod_mpoly_t Bbar, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    # Do the operation of :func:`nmod_mpoly_gcd` and also compute `Abar = A/G` and `Bbar = B/G` if successful.

    int nmod_mpoly_gcd_brown(nmod_mpoly_t G, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_gcd_hensel(nmod_mpoly_t G, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_gcd_zippel(nmod_mpoly_t G, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    # Try to set *G* to the GCD of *A* and *B* using various algorithms.

    int nmod_mpoly_resultant(nmod_mpoly_t R, const nmod_mpoly_t A, const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx)
    # Try to set *R* to the resultant of *A* and *B* with respect to the variable of index *var*.

    int nmod_mpoly_discriminant(nmod_mpoly_t D, const nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)
    # Try to set *D* to the discriminant of *A* with respect to the variable of index *var*.

    int nmod_mpoly_sqrt(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # If `Q^2=A` has a solution, set *Q* to a solution and return `1`, otherwise return `0` and set *Q* to zero.

    bint nmod_mpoly_is_square(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    # Return `1` if *A* is a perfect square, otherwise return `0`.

    int nmod_mpoly_quadratic_root(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    # If `Q^2+AQ=B` has a solution, set *Q* to a solution and return `1`, otherwise return `0`.

    void nmod_mpoly_univar_init(nmod_mpoly_univar_t A, const nmod_mpoly_ctx_t ctx)
    # Initialize *A*.

    void nmod_mpoly_univar_clear(nmod_mpoly_univar_t A, const nmod_mpoly_ctx_t ctx)
    # Clear *A*.

    void nmod_mpoly_univar_swap(nmod_mpoly_univar_t A, nmod_mpoly_univar_t B, const nmod_mpoly_ctx_t ctx)
    # Swap *A* and *B*.

    void nmod_mpoly_to_univar(nmod_mpoly_univar_t A, const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx)
    # Set *A* to a univariate form of *B* by pulling out the variable of index *var*.
    # The coefficients of *A* will still belong to the content *ctx* but will not depend on the variable of index *var*.

    void nmod_mpoly_from_univar(nmod_mpoly_t A, const nmod_mpoly_univar_t B, slong var, const nmod_mpoly_ctx_t ctx)
    # Set *A* to the normal form of *B* by putting in the variable of index *var*.
    # This function is undefined if the coefficients of *B* depend on the variable of index *var*.

    int nmod_mpoly_univar_degree_fits_si(const nmod_mpoly_univar_t A, const nmod_mpoly_ctx_t ctx)
    # Return `1` if the degree of *A* with respect to the main variable fits an ``slong``. Otherwise, return `0`.

    slong nmod_mpoly_univar_length(const nmod_mpoly_univar_t A, const nmod_mpoly_ctx_t ctx)
    # Return the number of terms in *A* with respect to the main variable.

    slong nmod_mpoly_univar_get_term_exp_si(nmod_mpoly_univar_t A, slong i, const nmod_mpoly_ctx_t ctx)
    # Return the exponent of the term of index *i* of *A*.

    void nmod_mpoly_univar_get_term_coeff(nmod_mpoly_t c, const nmod_mpoly_univar_t A, slong i, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_univar_swap_term_coeff(nmod_mpoly_t c, nmod_mpoly_univar_t A, slong i, const nmod_mpoly_ctx_t ctx)
    # Set (resp. swap) *c* to (resp. with) the coefficient of the term of index *i* of *A*.

    void nmod_mpoly_pow_rmul(nmod_mpoly_t A, const nmod_mpoly_t B, ulong k, const nmod_mpoly_ctx_t ctx)
    # Set *A* to *B* raised to the *k*-th power using repeated multiplications.

    void nmod_mpoly_div_monagan_pearce(nmod_mpoly_t polyq, const nmod_mpoly_t poly2, const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx)
    # Set ``polyq`` to the quotient of ``poly2`` by ``poly3``,
    # discarding the remainder (with notional remainder coefficients reduced
    # modulo the leading coefficient of ``poly3``). Implements "Polynomial
    # division using dynamic arrays, heaps and packed exponents" by Michael
    # Monagan and Roman Pearce. This function is exceptionally efficient if the
    # division is known to be exact.

    void nmod_mpoly_divrem_monagan_pearce(nmod_mpoly_t q, nmod_mpoly_t r, const nmod_mpoly_t poly2, const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx)
    # Set ``polyq`` and ``polyr`` to the quotient and remainder of
    # ``poly2`` divided by ``poly3``, (with remainder coefficients reduced
    # modulo the leading coefficient of ``poly3``). Implements "Polynomial
    # division using dynamic arrays, heaps and packed exponents" by Michael
    # Monagan and Roman Pearce.

    void nmod_mpoly_divrem_ideal_monagan_pearce(nmod_mpoly_struct ** q, nmod_mpoly_t r, const nmod_mpoly_t poly2, nmod_mpoly_struct * const * poly3, slong len, const nmod_mpoly_ctx_t ctx)
    # This function is as per ``nmod_mpoly_divrem_monagan_pearce`` except
    # that it takes an array of divisor polynomials ``poly3``, and it returns
    # an array of quotient polynomials ``q``. The number of divisor (and hence
    # quotient) polynomials, is given by *len*. The function computes
    # polynomials `q_i = q[i]` such that ``poly2`` is
    # `r + \sum_{i=0}^{\mbox{len - 1}} q_ib_i`, where `b_i =` ``poly3[i]``.
