# distutils: libraries = flint
# distutils: depends = flint/fmpz_mod_mpoly.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void fmpz_mod_mpoly_ctx_init(fmpz_mod_mpoly_ctx_t ctx, slong nvars, const ordering_t ord, const fmpz_t p)
    # Initialise a context object for a polynomial ring modulo *n* with *nvars* variables and ordering *ord*.
    # The possibilities for the ordering are ``ORD_LEX``, ``ORD_DEGLEX`` and ``ORD_DEGREVLEX``.

    slong fmpz_mod_mpoly_ctx_nvars(const fmpz_mod_mpoly_ctx_t ctx)
    # Return the number of variables used to initialize the context.

    ordering_t fmpz_mod_mpoly_ctx_ord(const fmpz_mod_mpoly_ctx_t ctx)
    # Return the ordering used to initialize the context.

    void fmpz_mod_mpoly_ctx_get_modulus(fmpz_t n, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *n* to the modulus used to initialize the context.

    void fmpz_mod_mpoly_ctx_clear(fmpz_mod_mpoly_ctx_t ctx)
    # Release up any space allocated by an *ctx*.

    void fmpz_mod_mpoly_init(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Initialise *A* for use with the given an initialised context object. Its value is set to zero.

    void fmpz_mod_mpoly_init2(fmpz_mod_mpoly_t A, slong alloc, const fmpz_mod_mpoly_ctx_t ctx)
    # Initialise *A* for use with the given an initialised context object. Its value is set to zero.
    # It is allocated with space for *alloc* terms and at least ``MPOLY_MIN_BITS`` bits for the exponents.

    void fmpz_mod_mpoly_init3(fmpz_mod_mpoly_t A, slong alloc, flint_bitcnt_t bits, const fmpz_mod_mpoly_ctx_t ctx)
    # Initialise *A* for use with the given an initialised context object. Its value is set to zero.
    # It is allocated with space for *alloc* terms and *bits* bits for the exponents.

    void fmpz_mod_mpoly_clear(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Release any space allocated for *A*.

    char * fmpz_mod_mpoly_get_str_pretty(const fmpz_mod_mpoly_t A, const char ** x, const fmpz_mod_mpoly_ctx_t ctx)
    # Return a string, which the user is responsible for cleaning up, representing *A*, given an array of variable strings *x*.

    int fmpz_mod_mpoly_fprint_pretty(FILE * file, const fmpz_mod_mpoly_t A, const char ** x, const fmpz_mod_mpoly_ctx_t ctx)
    # Print a string representing *A* to *file*.

    int fmpz_mod_mpoly_print_pretty(const fmpz_mod_mpoly_t A, const char ** x, const fmpz_mod_mpoly_ctx_t ctx)
    # Print a string representing *A* to ``stdout``.

    int fmpz_mod_mpoly_set_str_pretty(fmpz_mod_mpoly_t A, const char * str, const char ** x, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to the polynomial in the null-terminates string *str* given an array *x* of variable strings.
    # If parsing *str* fails, *A* is set to zero, and `-1` is returned. Otherwise, `0`  is returned.
    # The operations ``+``, ``-``, ``*``, and ``/`` are permitted along with integers and the variables in *x*. The character ``^`` must be immediately followed by the (integer) exponent.
    # If any division is not exact, parsing fails.

    void fmpz_mod_mpoly_gen(fmpz_mod_mpoly_t A, slong var, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to the variable of index *var*, where `var = 0` corresponds to the variable with the most significance with respect to the ordering.

    bint fmpz_mod_mpoly_is_gen(const fmpz_mod_mpoly_t A, slong var, const fmpz_mod_mpoly_ctx_t ctx)
    # If `var \ge 0`, return `1` if *A* is equal to the `var`-th generator, otherwise return `0`.
    # If `var < 0`, return `1` if the polynomial is equal to any generator, otherwise return `0`.

    void fmpz_mod_mpoly_set(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to *B*.

    bint fmpz_mod_mpoly_equal(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
    # Return `1` if *A* is equal to *B*, else return `0`.

    void fmpz_mod_mpoly_swap(fmpz_mod_mpoly_t poly1, fmpz_mod_mpoly_t poly2, const fmpz_mod_mpoly_ctx_t ctx)
    # Efficiently swap *A* and *B*.

    bint fmpz_mod_mpoly_is_fmpz(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Return `1` if *A* is a constant, else return `0`.

    void fmpz_mod_mpoly_get_fmpz(fmpz_t c, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Assuming that *A* is a constant, set *c* to this constant.
    # This function throws if *A* is not a constant.

    void fmpz_mod_mpoly_set_fmpz(fmpz_mod_mpoly_t A, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_set_ui(fmpz_mod_mpoly_t A, ulong c, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_set_si(fmpz_mod_mpoly_t A, slong c, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to the constant *c*.

    void fmpz_mod_mpoly_zero(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to the constant `0`.

    void fmpz_mod_mpoly_one(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to the constant `1`.

    bint fmpz_mod_mpoly_equal_fmpz(const fmpz_mod_mpoly_t A, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
    bint fmpz_mod_mpoly_equal_ui(const fmpz_mod_mpoly_t A, ulong c, const fmpz_mod_mpoly_ctx_t ctx)
    bint fmpz_mod_mpoly_equal_si(const fmpz_mod_mpoly_t A, slong c, const fmpz_mod_mpoly_ctx_t ctx)
    # Return `1` if *A* is equal to the constant *c*, else return `0`.

    bint fmpz_mod_mpoly_is_zero(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Return `1` if *A* is the constant `0`, else return `0`.

    bint fmpz_mod_mpoly_is_one(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Return `1` if *A* is the constant `1`, else return `0`.

    int fmpz_mod_mpoly_degrees_fit_si(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Return `1` if the degrees of *A* with respect to each variable fit into an ``slong``, otherwise return `0`.

    void fmpz_mod_mpoly_degrees_fmpz(fmpz ** degs, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_degrees_si(slong * degs, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *degs* to the degrees of *A* with respect to each variable.
    # If *A* is zero, all degrees are set to `-1`.

    void fmpz_mod_mpoly_degree_fmpz(fmpz_t deg, const fmpz_mod_mpoly_t A, slong var, const fmpz_mod_mpoly_ctx_t ctx)
    slong fmpz_mod_mpoly_degree_si(const fmpz_mod_mpoly_t A, slong var, const fmpz_mod_mpoly_ctx_t ctx)
    # Either return or set *deg* to the degree of *A* with respect to the variable of index *var*.
    # If *A* is zero, the degree is defined to be `-1`.

    int fmpz_mod_mpoly_total_degree_fits_si(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Return `1` if the total degree of *A* fits into an ``slong``, otherwise return `0`.

    void fmpz_mod_mpoly_total_degree_fmpz(fmpz_t tdeg, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    slong fmpz_mod_mpoly_total_degree_si(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Either return or set *tdeg* to the total degree of *A*.
    # If *A* is zero, the total degree is defined to be `-1`.

    void fmpz_mod_mpoly_used_vars(int * used, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # For each variable index *i*, set ``used[i]`` to nonzero if the variable of index *i* appears in *A* and to zero otherwise.

    void fmpz_mod_mpoly_get_coeff_fmpz_monomial(fmpz_t c, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t M, const fmpz_mod_mpoly_ctx_t ctx)
    # Assuming that *M* is a monomial, set *c* to the coefficient of the corresponding monomial in *A*.
    # This function throws if *M* is not a monomial.

    void fmpz_mod_mpoly_set_coeff_fmpz_monomial(fmpz_mod_mpoly_t A, const fmpz_t c, const fmpz_mod_mpoly_t M, const fmpz_mod_mpoly_ctx_t ctx)
    # Assuming that *M* is a monomial, set the coefficient of the corresponding monomial in *A* to *c*.
    # This function throws if *M* is not a monomial.

    void fmpz_mod_mpoly_get_coeff_fmpz_fmpz(fmpz_t c, const fmpz_mod_mpoly_t A, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_get_coeff_fmpz_ui(fmpz_t c, const fmpz_mod_mpoly_t A, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *c* to the coefficient of the monomial with exponent vector *exp*.

    void fmpz_mod_mpoly_set_coeff_fmpz_fmpz(fmpz_mod_mpoly_t A, const fmpz_t c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_set_coeff_ui_fmpz(fmpz_mod_mpoly_t A, ulong c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_set_coeff_si_fmpz(fmpz_mod_mpoly_t A, slong c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_set_coeff_fmpz_ui(fmpz_mod_mpoly_t A, const fmpz_t c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_set_coeff_ui_ui(fmpz_mod_mpoly_t A, ulong c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_set_coeff_si_ui(fmpz_mod_mpoly_t A, slong c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
    # Set the coefficient of the monomial with exponent vector *exp* to *c*.

    void fmpz_mod_mpoly_get_coeff_vars_ui(fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_t A, const slong * vars, const ulong * exps, slong length, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *C* to the coefficient of *A* with respect to the variables in *vars* with powers in the corresponding array *exps*.
    # Both *vars* and *exps* point to array of length *length*. It is assumed that `0 < length \le nvars(A)` and that the variables in *vars* are distinct.

    int fmpz_mod_mpoly_cmp(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
    # Return `1` (resp. `-1`, or `0`) if *A* is after (resp. before, same as) *B* in some arbitrary but fixed total ordering of the polynomials.
    # This ordering agrees with the usual ordering of monomials when *A* and *B* are both monomials.

    bint fmpz_mod_mpoly_is_canonical(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Return `1` if *A* is in canonical form. Otherwise, return `0`.
    # To be in canonical form, all of the terms must have nonzero coefficient, and the terms must be sorted from greatest to least.

    slong fmpz_mod_mpoly_length(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Return the number of terms in *A*.
    # If the polynomial is in canonical form, this will be the number of nonzero coefficients.

    void fmpz_mod_mpoly_resize(fmpz_mod_mpoly_t A, slong new_length, const fmpz_mod_mpoly_ctx_t ctx)
    # Set the length of *A* to ``new_length``.
    # Terms are either deleted from the end, or new zero terms are appended.

    void fmpz_mod_mpoly_get_term_coeff_fmpz(fmpz_t c, const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *c* to the coefficient of the term of index *i*.

    void fmpz_mod_mpoly_set_term_coeff_fmpz(fmpz_mod_mpoly_t A, slong i, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_set_term_coeff_ui(fmpz_mod_mpoly_t A, slong i, ulong c, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_set_term_coeff_si(fmpz_mod_mpoly_t A, slong i, slong c, const fmpz_mod_mpoly_ctx_t ctx)
    # Set the coefficient of the term of index *i* to *c*.

    int fmpz_mod_mpoly_term_exp_fits_si(const fmpz_mod_mpoly_t poly, slong i, const fmpz_mod_mpoly_ctx_t ctx)
    int fmpz_mod_mpoly_term_exp_fits_ui(const fmpz_mod_mpoly_t poly, slong i, const fmpz_mod_mpoly_ctx_t ctx)
    # Return `1` if all entries of the exponent vector of the term of index *i* fit into an ``slong`` (resp. a ``ulong``). Otherwise, return `0`.

    void fmpz_mod_mpoly_get_term_exp_fmpz(fmpz ** exp, const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_get_term_exp_ui(ulong * exp, const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_get_term_exp_si(slong * exp, const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *exp* to the exponent vector of the term of index *i*.
    # The ``_ui`` (resp. ``_si``) version throws if any entry does not fit into a ``ulong`` (resp. ``slong``).

    ulong fmpz_mod_mpoly_get_term_var_exp_ui(const fmpz_mod_mpoly_t A, slong i, slong var, const fmpz_mod_mpoly_ctx_t ctx)
    slong fmpz_mod_mpoly_get_term_var_exp_si(const fmpz_mod_mpoly_t A, slong i, slong var, const fmpz_mod_mpoly_ctx_t ctx)
    # Return the exponent of the variable *var* of the term of index *i*.
    # This function throws if the exponent does not fit into a ``ulong`` (resp. ``slong``).

    void fmpz_mod_mpoly_set_term_exp_fmpz(fmpz_mod_mpoly_t A, slong i, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_set_term_exp_ui(fmpz_mod_mpoly_t A, slong i, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
    # Set the exponent vector of the term of index *i* to *exp*.

    void fmpz_mod_mpoly_get_term(fmpz_mod_mpoly_t M, const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *M* to the term of index *i* in *A*.

    void fmpz_mod_mpoly_get_term_monomial(fmpz_mod_mpoly_t M, const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *M* to the monomial of the term of index *i* in *A*. The coefficient of *M* will be one.

    void fmpz_mod_mpoly_push_term_fmpz_fmpz(fmpz_mod_mpoly_t A, const fmpz_t c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_push_term_fmpz_ffmpz(fmpz_mod_mpoly_t A, const fmpz_t c, const fmpz * exp, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_push_term_ui_fmpz(fmpz_mod_mpoly_t A, ulong c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_push_term_ui_ffmpz(fmpz_mod_mpoly_t A, ulong c, const fmpz * exp, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_push_term_si_fmpz(fmpz_mod_mpoly_t A, slong c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_push_term_si_ffmpz(fmpz_mod_mpoly_t A, slong c, const fmpz * exp, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_push_term_fmpz_ui(fmpz_mod_mpoly_t A, const fmpz_t c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_push_term_ui_ui(fmpz_mod_mpoly_t A, ulong c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_push_term_si_ui(fmpz_mod_mpoly_t A, slong c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
    # Append a term to *A* with coefficient *c* and exponent vector *exp*.
    # This function runs in constant average time.

    void fmpz_mod_mpoly_sort_terms(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Sort the terms of *A* into the canonical ordering dictated by the ordering in *ctx*.
    # This function simply reorders the terms: It does not combine like terms, nor does it delete terms with coefficient zero.
    # This function runs in linear time in the size of *A*.

    void fmpz_mod_mpoly_combine_like_terms(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Combine adjacent like terms in *A* and delete terms with coefficient zero.
    # If the terms of *A* were sorted to begin with, the result will be in canonical form.
    # This function runs in linear time in the size of *A*.

    void fmpz_mod_mpoly_reverse(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to the reversal of *B*.

    void fmpz_mod_mpoly_randtest_bound(fmpz_mod_mpoly_t A, flint_rand_t state, slong length, ulong exp_bound, const fmpz_mod_mpoly_ctx_t ctx)
    # Generate a random polynomial with length up to *length* and exponents in the range ``[0, exp_bound - 1]``.
    # The exponents of each variable are generated by calls to ``n_randint(state, exp_bound)``.

    void fmpz_mod_mpoly_randtest_bounds(fmpz_mod_mpoly_t A, flint_rand_t state, slong length, ulong * exp_bounds, const fmpz_mod_mpoly_ctx_t ctx)
    # Generate a random polynomial with length up to *length* and exponents in the range ``[0, exp_bounds[i] - 1]``.
    # The exponents of the variable of index *i* are generated by calls to ``n_randint(state, exp_bounds[i])``.

    void fmpz_mod_mpoly_randtest_bits(fmpz_mod_mpoly_t A, flint_rand_t state, slong length, mp_limb_t exp_bits, const fmpz_mod_mpoly_ctx_t ctx)
    # Generate a random polynomial with length up to *length* and exponents whose packed form does not exceed the given bit count.

    void fmpz_mod_mpoly_add_fmpz(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_add_ui(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, ulong c, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_add_si(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, slong c, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to `B + c`.

    void fmpz_mod_mpoly_sub_fmpz(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_sub_ui(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, ulong c, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_sub_si(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, slong c, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to `B - c`.

    void fmpz_mod_mpoly_add(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to `B + C`.

    void fmpz_mod_mpoly_sub(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to `B - C`.

    void fmpz_mod_mpoly_neg(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to `-B`.

    void fmpz_mod_mpoly_scalar_mul_fmpz(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_scalar_mul_ui(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, ulong c, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_scalar_mul_si(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, slong c, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to `B \times c`.

    void fmpz_mod_mpoly_scalar_addmul_fmpz(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C, const fmpz_t d, const fmpz_mod_mpoly_ctx_t ctx)
    # Sets *A* to `B + C \times d`.

    void fmpz_mod_mpoly_make_monic(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to *B* divided by the leading coefficient of *B*. This throws if *B* is zero or the leading coefficient is not invertible.

    void fmpz_mod_mpoly_derivative(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to the derivative of *B* with respect to the variable of index *var*.

    void fmpz_mod_mpoly_evaluate_all_fmpz(fmpz_t eval, const fmpz_mod_mpoly_t A, fmpz * const * vals, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *ev* to the evaluation of *A* where the variables are replaced by the corresponding elements of the array *vals*.

    void fmpz_mod_mpoly_evaluate_one_fmpz(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, slong var, const fmpz_t val, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to the evaluation of *B* where the variable of index *var* is replaced by *val*.
    # Return `1` for success and `0` for failure.

    int fmpz_mod_mpoly_compose_fmpz_poly(fmpz_poly_t A, const fmpz_mod_mpoly_t B, fmpz_poly_struct * const * C, const fmpz_mod_mpoly_ctx_t ctxB)
    # Set *A* to the evaluation of *B* where the variables are replaced by the corresponding elements of the array *C*.
    # The context object of *B* is *ctxB*.
    # Return `1` for success and `0` for failure.

    int fmpz_mod_mpoly_compose_fmpz_mod_mpoly_geobucket(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, fmpz_mod_mpoly_struct * const * C, const fmpz_mod_mpoly_ctx_t ctxB, const fmpz_mod_mpoly_ctx_t ctxAC)
    int fmpz_mod_mpoly_compose_fmpz_mod_mpoly(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, fmpz_mod_mpoly_struct * const * C, const fmpz_mod_mpoly_ctx_t ctxB, const fmpz_mod_mpoly_ctx_t ctxAC)
    # Set *A* to the evaluation of *B* where the variables are replaced by the corresponding elements of the array *C*.
    # Both *A* and the elements of *C* have context object *ctxAC*, while *B* has context object *ctxB*.
    # The length of the array *C* is the number of variables in *ctxB*.
    # Neither *A* nor *B* is allowed to alias any other polynomial.
    # Return `1` for success and `0` for failure.
    # The main method attempts to perform the calculation using matrices and chooses heuristically between the ``geobucket`` and ``horner`` methods if needed.

    void fmpz_mod_mpoly_compose_fmpz_mod_mpoly_gen(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const slong * c, const fmpz_mod_mpoly_ctx_t ctxB, const fmpz_mod_mpoly_ctx_t ctxAC)
    # Set *A* to the evaluation of *B* where the variable of index *i* in *ctxB* is replaced by the variable of index ``c[i]`` in *ctxAC*.
    # The length of the array *C* is the number of variables in *ctxB*.
    # If any ``c[i]`` is negative, the corresponding variable of *B* is replaced by zero. Otherwise, it is expected that ``c[i]`` is less than the number of variables in *ctxAC*.

    void fmpz_mod_mpoly_mul(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to `B \times C`.

    void fmpz_mod_mpoly_mul_johnson(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to `B \times C` using Johnson's heap-based method.

    int fmpz_mod_mpoly_mul_dense(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx)
    # Try to set *A* to `B \times C` using dense arithmetic.
    # If the return is `0`, the operation was unsuccessful. Otherwise, it was successful and the return is `1`.

    int fmpz_mod_mpoly_pow_fmpz(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_t k, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to *B* raised to the `k`-th power.
    # Return `1` for success and `0` for failure.

    int fmpz_mod_mpoly_pow_ui(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, ulong k, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to *B* raised to the `k`-th power.
    # Return `1` for success and `0` for failure.

    int fmpz_mod_mpoly_divides(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
    # If *A* is divisible by *B*, set *Q* to the exact quotient and return `1`. Otherwise, set *Q* to zero and return `0`.

    void fmpz_mod_mpoly_div(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *Q* to the quotient of *A* by *B*, discarding the remainder.

    void fmpz_mod_mpoly_divrem(fmpz_mod_mpoly_t Q, fmpz_mod_mpoly_t R, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *Q* and *R* to the quotient and remainder of *A* divided by *B*.

    void fmpz_mod_mpoly_divrem_ideal(fmpz_mod_mpoly_struct ** Q, fmpz_mod_mpoly_t R, const fmpz_mod_mpoly_t A, fmpz_mod_mpoly_struct * const * B, slong len, const fmpz_mod_mpoly_ctx_t ctx)
    # This function is as per :func:`fmpz_mod_mpoly_divrem` except that it takes an array of divisor polynomials *B* and it returns an array of quotient polynomials *Q*.
    # The number of divisor (and hence quotient) polynomials, is given by *len*.

    void fmpz_mod_mpoly_term_content(fmpz_mod_mpoly_t M, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *M* to the GCD of the terms of *A*.
    # If *A* is zero, *M* will be zero. Otherwise, *M* will be a monomial with coefficient one.

    int fmpz_mod_mpoly_content_vars(fmpz_mod_mpoly_t g, const fmpz_mod_mpoly_t A, slong * vars, slong vars_length, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *g* to the GCD of the coefficients of *A* when viewed as a polynomial in the variables *vars*.
    # Return `1` for success and `0` for failure. Upon success, *g* will be independent of the variables *vars*.

    int fmpz_mod_mpoly_gcd(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
    # Try to set *G* to the monic GCD of *A* and *B*. The GCD of zero and zero is defined to be zero.
    # If the return is `1` the function was successful. Otherwise the return is `0` and *G* is left untouched.

    int fmpz_mod_mpoly_gcd_cofactors(fmpz_mod_mpoly_t G, fmpz_mod_mpoly_t Abar, fmpz_mod_mpoly_t Bbar, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
    # Do the operation of :func:`fmpz_mod_mpoly_gcd` and also compute `Abar = A/G` and `Bbar = B/G` if successful.

    int fmpz_mod_mpoly_gcd_brown(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
    int fmpz_mod_mpoly_gcd_hensel(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
    int fmpz_mod_mpoly_gcd_subresultant(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
    int fmpz_mod_mpoly_gcd_zippel(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
    int fmpz_mod_mpoly_gcd_zippel2(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
    # Try to set *G* to the GCD of *A* and *B* using various algorithms.

    int fmpz_mod_mpoly_resultant(fmpz_mod_mpoly_t R, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx)
    # Try to set *R* to the resultant of *A* and *B* with respect to the variable of index *var*.

    int fmpz_mod_mpoly_discriminant(fmpz_mod_mpoly_t D, const fmpz_mod_mpoly_t A, slong var, const fmpz_mod_mpoly_ctx_t ctx)
    # Try to set *D* to the discriminant of *A* with respect to the variable of index *var*.

    int fmpz_mod_mpoly_sqrt(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # If `Q^2=A` has a solution, set *Q* to a solution and return `1`, otherwise return `0` and set *Q* to zero.

    bint fmpz_mod_mpoly_is_square(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Return `1` if *A* is a perfect square, otherwise return `0`.

    int fmpz_mod_mpoly_quadratic_root(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
    # If `Q^2+AQ=B` has a solution, set *Q* to a solution and return `1`, otherwise return `0`.

    void fmpz_mod_mpoly_univar_init(fmpz_mod_mpoly_univar_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Initialize *A*.

    void fmpz_mod_mpoly_univar_clear(fmpz_mod_mpoly_univar_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Clear *A*.

    void fmpz_mod_mpoly_univar_swap(fmpz_mod_mpoly_univar_t A, fmpz_mod_mpoly_univar_t B, const fmpz_mod_mpoly_ctx_t ctx)
    # Swap *A* and *B*.

    void fmpz_mod_mpoly_to_univar(fmpz_mod_mpoly_univar_t A, const fmpz_mod_mpoly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to a univariate form of *B* by pulling out the variable of index *var*.
    # The coefficients of *A* will still belong to the content *ctx* but will not depend on the variable of index *var*.

    void fmpz_mod_mpoly_from_univar(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_univar_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx)
    # Set *A* to the normal form of *B* by putting in the variable of index *var*.
    # This function is undefined if the coefficients of *B* depend on the variable of index *var*.

    int fmpz_mod_mpoly_univar_degree_fits_si(const fmpz_mod_mpoly_univar_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Return `1` if the degree of *A* with respect to the main variable fits an ``slong``. Otherwise, return `0`.

    slong fmpz_mod_mpoly_univar_length(const fmpz_mod_mpoly_univar_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # Return the number of terms in *A* with respect to the main variable.

    slong fmpz_mod_mpoly_univar_get_term_exp_si(fmpz_mod_mpoly_univar_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
    # Return the exponent of the term of index *i* of *A*.

    void fmpz_mod_mpoly_univar_get_term_coeff(fmpz_mod_mpoly_t c, const fmpz_mod_mpoly_univar_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_univar_swap_term_coeff(fmpz_mod_mpoly_t c, fmpz_mod_mpoly_univar_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
    # Set (resp. swap) *c* to (resp. with) the coefficient of the term of index *i* of *A*.

    void fmpz_mod_mpoly_univar_set_coeff_ui(fmpz_mod_mpoly_univar_t Ax, ulong e, const fmpz_mod_mpoly_t c, const fmpz_mod_mpoly_ctx_t ctx)
    # Set the coefficient of `X^e` in *Ax* to *c*.

    int fmpz_mod_mpoly_univar_resultant(fmpz_mod_mpoly_t R, const fmpz_mod_mpoly_univar_t Ax, const fmpz_mod_mpoly_univar_t Bx, const fmpz_mod_mpoly_ctx_t ctx)
    # Try to set *R* to the resultant of *Ax* and *Bx*.

    int fmpz_mod_mpoly_univar_discriminant(fmpz_mod_mpoly_t D, const fmpz_mod_mpoly_univar_t Ax, const fmpz_mod_mpoly_ctx_t ctx)
    # Try to set *D* to the discriminant of *Ax*.

    void fmpz_mod_mpoly_inflate(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz * shift, const fmpz * stride, const fmpz_mod_mpoly_ctx_t ctx)
    # Apply the function ``e -> shift[v] + stride[v]*e`` to each exponent ``e`` corresponding to the variable ``v``.
    # It is assumed that each shift and stride is not negative.

    void fmpz_mod_mpoly_deflate(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz * shift, const fmpz * stride, const fmpz_mod_mpoly_ctx_t ctx)
    # Apply the function ``e -> (e - shift[v])/stride[v]`` to each exponent ``e`` corresponding to the variable ``v``.
    # If any ``stride[v]`` is zero, the corresponding numerator ``e - shift[v]`` is assumed to be zero, and the quotient is defined as zero.
    # This allows the function to undo the operation performed by :func:`fmpz_mod_mpoly_inflate` when possible.

    void fmpz_mod_mpoly_deflation(fmpz * shift, fmpz * stride, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    # For each variable `v` let `S_v` be the set of exponents appearing on `v`.
    # Set ``shift[v]`` to `\operatorname{min}(S_v)` and set ``stride[v]`` to `\operatorname{gcd}(S-\operatorname{min}(S_v))`.
    # If *A* is zero, all shifts and strides are set to zero.
