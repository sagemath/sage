# distutils: libraries = flint
# distutils: depends = flint/fexpr.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void fexpr_init(fexpr_t expr) noexcept
    # Initializes *expr* for use. Its value is set to the atomic
    # integer 0.

    void fexpr_clear(fexpr_t expr) noexcept
    # Clears *expr*, freeing its allocated memory.

    fexpr_ptr _fexpr_vec_init(slong len) noexcept
    # Returns a heap-allocated vector of *len* initialized expressions.

    void _fexpr_vec_clear(fexpr_ptr vec, slong len) noexcept
    # Clears the *len* expressions in *vec* and frees *vec* itself.

    void fexpr_fit_size(fexpr_t expr, slong size) noexcept
    # Ensures that *expr* has room for *size* words.

    void fexpr_set(fexpr_t res, const fexpr_t expr) noexcept
    # Sets *res* to the a copy of *expr*.

    void fexpr_swap(fexpr_t a, fexpr_t b) noexcept
    # Swaps *a* and *b* efficiently.

    slong fexpr_depth(const fexpr_t expr) noexcept
    # Returns the depth of *expr* as a symbolic expression tree.

    slong fexpr_num_leaves(const fexpr_t expr) noexcept
    # Returns the number of leaves (atoms, counted with repetition)
    # in the expression *expr*.

    slong fexpr_size(const fexpr_t expr) noexcept
    # Returns the number of words in the internal representation
    # of *expr*.

    slong fexpr_size_bytes(const fexpr_t expr) noexcept
    # Returns the number of bytes in the internal representation
    # of *expr*. The count excludes the size of the structure itself.
    # Add ``sizeof(fexpr_struct)`` to get the size of the object as a
    # whole.

    slong fexpr_allocated_bytes(const fexpr_t expr) noexcept
    # Returns the number of allocated bytes in the internal
    # representation of *expr*. The count excludes the size of the
    # structure itself. Add ``sizeof(fexpr_struct)`` to get the size of
    # the object as a whole.

    bint fexpr_equal(const fexpr_t a, const fexpr_t b) noexcept
    # Checks if *a* and *b* are exactly equal as expressions.

    bint fexpr_equal_si(const fexpr_t expr, slong c) noexcept

    bint fexpr_equal_ui(const fexpr_t expr, ulong c) noexcept
    # Checks if *expr* is an atomic integer exactly equal to *c*.

    ulong fexpr_hash(const fexpr_t expr) noexcept
    # Returns a hash of the expression *expr*.

    int fexpr_cmp_fast(const fexpr_t a, const fexpr_t b) noexcept
    # Compares *a* and *b* using an ordering based on the internal
    # representation, returning -1, 0 or 1. This can be used, for
    # instance, to maintain sorted arrays of expressions for binary
    # search; the sort order has no mathematical significance.

    bint fexpr_is_integer(const fexpr_t expr) noexcept
    # Returns whether *expr* is an atomic integer

    bint fexpr_is_symbol(const fexpr_t expr) noexcept
    # Returns whether *expr* is an atomic symbol.

    bint fexpr_is_string(const fexpr_t expr) noexcept
    # Returns whether *expr* is an atomic string.

    bint fexpr_is_atom(const fexpr_t expr) noexcept
    # Returns whether *expr* is any atom.

    void fexpr_zero(fexpr_t res) noexcept
    # Sets *res* to the atomic integer 0.

    bint fexpr_is_zero(const fexpr_t expr) noexcept
    # Returns whether *expr* is the atomic integer 0.

    bint fexpr_is_neg_integer(const fexpr_t expr) noexcept
    # Returns whether *expr* is any negative atomic integer.

    void fexpr_set_si(fexpr_t res, slong c) noexcept
    void fexpr_set_ui(fexpr_t res, ulong c) noexcept
    void fexpr_set_fmpz(fexpr_t res, const fmpz_t c) noexcept
    # Sets *res* to the atomic integer *c*.

    int fexpr_get_fmpz(fmpz_t res, const fexpr_t expr) noexcept
    # Sets *res* to the atomic integer in *expr*. This aborts
    # if *expr* is not an atomic integer.

    void fexpr_set_symbol_builtin(fexpr_t res, slong id) noexcept
    # Sets *res* to the builtin symbol with internal index *id*
    # (see :ref:`fexpr-builtin`).

    bint fexpr_is_builtin_symbol(const fexpr_t expr, slong id) noexcept
    # Returns whether *expr* is the builtin symbol with index *id*
    # (see :ref:`fexpr-builtin`).

    bint fexpr_is_any_builtin_symbol(const fexpr_t expr) noexcept
    # Returns whether *expr* is any builtin symbol
    # (see :ref:`fexpr-builtin`).

    void fexpr_set_symbol_str(fexpr_t res, const char * s) noexcept
    # Sets *res* to the symbol given by *s*.

    char * fexpr_get_symbol_str(const fexpr_t expr) noexcept
    # Returns the symbol in *expr* as a string. The string must
    # be freed with :func:`flint_free`.
    # This aborts if *expr* is not an atomic symbol.

    void fexpr_set_string(fexpr_t res, const char * s) noexcept
    # Sets *res* to the atomic string *s*.

    char * fexpr_get_string(const fexpr_t expr) noexcept
    # Assuming that *expr* is an atomic string, returns a copy of this
    # string. The string must be freed with :func:`flint_free`.

    void fexpr_write(calcium_stream_t stream, const fexpr_t expr) noexcept
    # Writes *expr* to *stream*.

    void fexpr_print(const fexpr_t expr) noexcept
    # Prints *expr* to standard output.

    char * fexpr_get_str(const fexpr_t expr) noexcept
    # Returns a string representation of *expr*. The string must
    # be freed with :func:`flint_free`.
    # Warning: string literals appearing in expressions
    # are currently not escaped.

    void fexpr_write_latex(calcium_stream_t stream, const fexpr_t expr, ulong flags) noexcept
    # Writes the LaTeX representation of *expr* to *stream*.

    void fexpr_print_latex(const fexpr_t expr, ulong flags) noexcept
    # Prints the LaTeX representation of *expr* to standard output.

    char * fexpr_get_str_latex(const fexpr_t expr, ulong flags) noexcept
    # Returns a string of the LaTeX representation of *expr*. The string
    # must be freed with :func:`flint_free`.
    # Warning: string literals appearing in expressions
    # are currently not escaped.

    slong fexpr_nargs(const fexpr_t expr) noexcept
    # Returns the number of arguments *n* in the function call
    # `f(e_1,\ldots,e_n)` represented
    # by *expr*. If *expr* is an atom, returns -1.

    void fexpr_func(fexpr_t res, const fexpr_t expr) noexcept
    # Assuming that *expr* represents a function call
    # `f(e_1,\ldots,e_n)`, sets *res* to the function expression *f*.

    void fexpr_view_func(fexpr_t view, const fexpr_t expr) noexcept
    # As :func:`fexpr_func`, but sets *view* to a shallow view
    # instead of copying the expression.
    # The variable *view* must not be initialized before use or
    # cleared after use, and *expr* must not be modified or cleared
    # as long as *view* is in use.

    void fexpr_arg(fexpr_t res, const fexpr_t expr, slong i) noexcept
    # Assuming that *expr* represents a function call
    # `f(e_1,\ldots,e_n)`, sets *res* to the argument `e_{i+1}`.
    # Note that indexing starts from 0.
    # The index must be in bounds, with `0 \le i < n`.

    void fexpr_view_arg(fexpr_t view, const fexpr_t expr, slong i) noexcept
    # As :func:`fexpr_arg`, but sets *view* to a shallow view
    # instead of copying the expression.
    # The variable *view* must not be initialized before use or
    # cleared after use, and *expr* must not be modified or cleared
    # as long as *view* is in use.

    void fexpr_view_next(fexpr_t view) noexcept
    # Assuming that *view* is a shallow view of a function argument `e_i`
    # in a function call `f(e_1,\ldots,e_n)`, sets *view* to
    # a view of the next argument `e_{i+1}`.
    # This function can be called when *view* refers to the last argument
    # `e_n`, provided that *view* is not used afterwards.
    # This function can also be called when *view* refers to the function *f*,
    # in which case it will make *view* point to `e_1`.

    bint fexpr_is_builtin_call(const fexpr_t expr, slong id) noexcept
    # Returns whether *expr* has the form `f(\ldots)` where *f* is
    # a builtin function defined by *id* (see :ref:`fexpr-builtin`).

    bint fexpr_is_any_builtin_call(const fexpr_t expr) noexcept
    # Returns whether *expr* has the form `f(\ldots)` where *f* is
    # any builtin function (see :ref:`fexpr-builtin`).

    void fexpr_call0(fexpr_t res, const fexpr_t f) noexcept
    void fexpr_call1(fexpr_t res, const fexpr_t f, const fexpr_t x1) noexcept
    void fexpr_call2(fexpr_t res, const fexpr_t f, const fexpr_t x1, const fexpr_t x2) noexcept
    void fexpr_call3(fexpr_t res, const fexpr_t f, const fexpr_t x1, const fexpr_t x2, const fexpr_t x3) noexcept
    void fexpr_call4(fexpr_t res, const fexpr_t f, const fexpr_t x1, const fexpr_t x2, const fexpr_t x3, const fexpr_t x4) noexcept
    void fexpr_call_vec(fexpr_t res, const fexpr_t f, fexpr_srcptr args, slong len) noexcept
    # Creates the function call `f(x_1,\ldots,x_n)`.
    # The *vec* version takes the arguments as an array *args*
    # and *n* is given by *len*.
    # Warning: aliasing between inputs and outputs is not implemented.

    void fexpr_call_builtin1(fexpr_t res, slong f, const fexpr_t x1) noexcept
    void fexpr_call_builtin2(fexpr_t res, slong f, const fexpr_t x1, const fexpr_t x2) noexcept
    # Creates the function call `f(x_1,\ldots,x_n)`, where *f* defines
    # a builtin symbol.

    bint fexpr_contains(const fexpr_t expr, const fexpr_t x) noexcept
    # Returns whether *expr* contains the expression *x* as a subexpression
    # (this includes the case where *expr* and *x* are equal).

    int fexpr_replace(fexpr_t res, const fexpr_t expr, const fexpr_t x, const fexpr_t y) noexcept
    # Sets *res* to the expression *expr* with all occurrences of the subexpression
    # *x* replaced by the expression *y*. Returns a boolean value indicating whether
    # any replacements have been performed.
    # Aliasing is allowed between *res* and *expr* but not between *res*
    # and *x* or *y*.

    int fexpr_replace2(fexpr_t res, const fexpr_t expr, const fexpr_t x1, const fexpr_t y1, const fexpr_t x2, const fexpr_t y2) noexcept
    # Like :func:`fexpr_replace`, but simultaneously replaces *x1* by *y1*
    # and *x2* by *y2*.

    int fexpr_replace_vec(fexpr_t res, const fexpr_t expr, const fexpr_vec_t xs, const fexpr_vec_t ys) noexcept
    # Sets *res* to the expression *expr* with all occurrences of the
    # subexpressions given by entries in *xs* replaced by the corresponding
    # expressions in *ys*. It is required that *xs* and *ys* have the same length.
    # Returns a boolean value indicating whether any replacements
    # have been performed.
    # Aliasing is allowed between *res* and *expr* but not between *res*
    # and the entries of *xs* or *ys*.

    void fexpr_set_fmpq(fexpr_t res, const fmpq_t x) noexcept
    # Sets *res* to the rational number *x*. This creates an atomic
    # integer if the denominator of *x* is one, and otherwise creates a
    # division expression.

    void fexpr_set_arf(fexpr_t res, const arf_t x) noexcept
    void fexpr_set_d(fexpr_t res, double x) noexcept
    # Sets *res* to an expression for the value of the
    # floating-point number *x*. NaN is represented
    # as ``Undefined``. For a regular value, this creates an atomic integer
    # or a rational fraction if the exponent is small, and otherwise
    # creates an expression of the form ``Mul(m, Pow(2, e))``.

    void fexpr_set_re_im_d(fexpr_t res, double x, double y) noexcept
    # Sets *res* to an expression for the complex number with real part
    # *x* and imaginary part *y*.

    void fexpr_neg(fexpr_t res, const fexpr_t a) noexcept
    void fexpr_add(fexpr_t res, const fexpr_t a, const fexpr_t b) noexcept
    void fexpr_sub(fexpr_t res, const fexpr_t a, const fexpr_t b) noexcept
    void fexpr_mul(fexpr_t res, const fexpr_t a, const fexpr_t b) noexcept
    void fexpr_div(fexpr_t res, const fexpr_t a, const fexpr_t b) noexcept
    void fexpr_pow(fexpr_t res, const fexpr_t a, const fexpr_t b) noexcept
    # Constructs an arithmetic expression with given arguments.
    # No simplifications whatsoever are performed.

    bint fexpr_is_arithmetic_operation(const fexpr_t expr) noexcept
    # Returns whether *expr* is of the form `f(e_1,\ldots,e_n)`
    # where *f* is one of the arithmetic operators ``Pos``, ``Neg``,
    # ``Add``, ``Sub``, ``Mul``, ``Div``.

    void fexpr_arithmetic_nodes(fexpr_vec_t nodes, const fexpr_t expr) noexcept
    # Sets *nodes* to a vector of subexpressions of *expr* such that *expr*
    # is an arithmetic expression with *nodes* as leaves.
    # More precisely, *expr* will be constructed out of nested application
    # the arithmetic operators
    # ``Pos``, ``Neg``, ``Add``, ``Sub``, ``Mul``, ``Div`` with
    # integers and expressions in *nodes* as leaves.
    # Powers ``Pow`` with an atomic integer exponent are also allowed.
    # The nodes are output without repetition but are not automatically sorted in
    # a canonical order.

    int fexpr_get_fmpz_mpoly_q(fmpz_mpoly_q_t res, const fexpr_t expr, const fexpr_vec_t vars, const fmpz_mpoly_ctx_t ctx) noexcept
    # Sets *res* to the expression *expr* as a formal rational
    # function of the subexpressions in *vars*.
    # The vector *vars* must have the same length as the number of
    # variables specified in *ctx*.
    # To build *vars* automatically for a given expression,
    # :func:`fexpr_arithmetic_nodes` may be used.
    # Returns 1 on success and 0 on failure. Failure can occur for the
    # following reasons:
    # * A subexpression is encountered that cannot be interpreted
    # as an arithmetic operation and does not appear (exactly) in *vars*.
    # * Overflow (too many terms or too large exponent).
    # * Division by zero (a zero denominator is encountered).
    # It is important to note that this function views *expr* as
    # a formal rational function with *vars* as formal indeterminates.
    # It does thus not check for algebraic relations between *vars*
    # and can implicitly divide by zero if *vars* are not algebraically
    # independent.

    void fexpr_set_fmpz_mpoly(fexpr_t res, const fmpz_mpoly_t poly, const fexpr_vec_t vars, const fmpz_mpoly_ctx_t ctx) noexcept
    void fexpr_set_fmpz_mpoly_q(fexpr_t res, const fmpz_mpoly_q_t frac, const fexpr_vec_t vars, const fmpz_mpoly_ctx_t ctx) noexcept
    # Sets *res* to an expression for the multivariate polynomial *poly*
    # (or rational function *frac*),
    # using the expressions in *vars* as the variables. The length
    # of *vars* must agree with the number of variables in *ctx*.
    # If *NULL* is passed for *vars*, a default choice of symbols
    # is used.

    int fexpr_expanded_normal_form(fexpr_t res, const fexpr_t expr, ulong flags) noexcept
    # Sets *res* to *expr* converted to expanded normal form viewed
    # as a formal rational function with its non-arithmetic subexpressions
    # as terminal nodes.
    # This function first computes nodes with :func:`fexpr_arithmetic_nodes`,
    # sorts the nodes, evaluates to a rational function with
    # :func:`fexpr_get_fmpz_mpoly_q`, and then converts back to an
    # expression with :func:`fexpr_set_fmpz_mpoly_q`.
    # Optional *flags* are reserved for future use.

    void fexpr_vec_init(fexpr_vec_t vec, slong len) noexcept
    # Initializes *vec* to a vector of length *len*. All entries
    # are set to the atomic integer 0.

    void fexpr_vec_clear(fexpr_vec_t vec) noexcept
    # Clears the vector *vec*.

    void fexpr_vec_print(const fexpr_vec_t vec) noexcept
    # Prints *vec* to standard output.

    void fexpr_vec_swap(fexpr_vec_t x, fexpr_vec_t y) noexcept
    # Swaps *x* and *y* efficiently.

    void fexpr_vec_fit_length(fexpr_vec_t vec, slong len) noexcept
    # Ensures that *vec* has space for *len* entries.

    void fexpr_vec_set(fexpr_vec_t dest, const fexpr_vec_t src) noexcept
    # Sets *dest* to a copy of *src*.

    void fexpr_vec_append(fexpr_vec_t vec, const fexpr_t expr) noexcept
    # Appends *expr* to the end of the vector *vec*.

    slong fexpr_vec_insert_unique(fexpr_vec_t vec, const fexpr_t expr) noexcept
    # Inserts *expr* without duplication into vec, returning its
    # position. If this expression already exists, *vec* is unchanged.
    # If this expression does not exist in *vec*, it is appended.

    void fexpr_vec_set_length(fexpr_vec_t vec, slong len) noexcept
    # Sets the length of *vec* to *len*, truncating or zero-extending as needed.

    void _fexpr_vec_sort_fast(fexpr_ptr vec, slong len) noexcept
    # Sorts the *len* entries in *vec* using
    # the comparison function :func:`fexpr_cmp_fast`.
