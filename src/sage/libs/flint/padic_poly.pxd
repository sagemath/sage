# distutils: libraries = flint
# distutils: depends = flint/padic_poly.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void padic_poly_init(padic_poly_t poly)
    # Initialises ``poly`` for use, setting its length to zero.
    # The precision of the polynomial is set to ``PADIC_DEFAULT_PREC``.
    # A corresponding call to :func:`padic_poly_clear` must be made
    # after finishing with the :type:`padic_poly_t` to free the memory
    # used by the polynomial.

    void padic_poly_init2(padic_poly_t poly, slong alloc, slong prec)
    # Initialises ``poly`` with space for at least ``alloc`` coefficients
    # and sets the length to zero.  The allocated coefficients are all set to
    # zero.  The precision is set to ``prec``.

    void padic_poly_realloc(padic_poly_t poly, slong alloc, const fmpz_t p)
    # Reallocates the given polynomial to have space for ``alloc``
    # coefficients.  If ``alloc`` is zero the polynomial is cleared
    # and then reinitialised.  If the current length is greater than
    # ``alloc`` the polynomial is first truncated to length ``alloc``.

    void padic_poly_fit_length(padic_poly_t poly, slong len)
    # If ``len`` is greater than the number of coefficients currently
    # allocated, then the polynomial is reallocated to have space for at
    # least ``len`` coefficients.  No data is lost when calling this
    # function.
    # The function efficiently deals with the case where ``fit_length`` is
    # called many times in small increments by at least doubling the number
    # of allocated coefficients when length is larger than the number of
    # coefficients currently allocated.

    void _padic_poly_set_length(padic_poly_t poly, slong len)
    # Demotes the coefficients of ``poly`` beyond ``len`` and sets
    # the length of ``poly`` to ``len``.
    # Note that if the current length is greater than ``len`` the
    # polynomial may no slonger be in canonical form.

    void padic_poly_clear(padic_poly_t poly)
    # Clears the given polynomial, releasing any memory used.  It must
    # be reinitialised in order to be used again.

    void _padic_poly_normalise(padic_poly_t poly)
    # Sets the length of ``poly`` so that the top coefficient is non-zero.
    # If all coefficients are zero, the length is set to zero.  This function
    # is mainly used internally, as all functions guarantee normalisation.

    void _padic_poly_canonicalise(fmpz *poly, slong *v, slong len, const fmpz_t p)
    void padic_poly_canonicalise(padic_poly_t poly, const fmpz_t p)
    # Brings the polynomial ``poly`` into canonical form,
    # assuming that it is normalised already.  Does *not*
    # carry out any reduction.

    void padic_poly_reduce(padic_poly_t poly, const padic_ctx_t ctx)
    # Reduces the polynomial ``poly`` modulo `p^N`, assuming
    # that it is in canonical form already.

    void padic_poly_truncate(padic_poly_t poly, slong n, const fmpz_t p)
    # Truncates the polynomial to length at most~`n`.

    slong padic_poly_degree(const padic_poly_t poly)
    # Returns the degree of the polynomial ``poly``.

    slong padic_poly_length(const padic_poly_t poly)
    # Returns the length of the polynomial ``poly``.

    slong padic_poly_val(const padic_poly_t poly)
    # Returns the valuation of the polynomial ``poly``,
    # which is defined to be the minimum valuation of all
    # its coefficients.
    # The valuation of the zero polynomial is~`0`.
    # Note that this is implemented as a macro and can be
    # used as either a ``lvalue`` or a ``rvalue``.

    slong padic_poly_prec(padic_poly_t poly)
    # Returns the precision of the polynomial ``poly``.
    # Note that this is implemented as a macro and can be
    # used as either a ``lvalue`` or a ``rvalue``.
    # Note that increasing the precision might require
    # a call to :func:`padic_poly_reduce`.

    void padic_poly_randtest(padic_poly_t f, flint_rand_t state, slong len, const padic_ctx_t ctx)
    # Sets `f` to a random polynomial of length at most ``len``
    # with entries reduced modulo `p^N`.

    void padic_poly_randtest_not_zero(padic_poly_t f, flint_rand_t state, slong len, const padic_ctx_t ctx)
    # Sets `f` to a non-zero random polynomial of length at most ``len``
    # with entries reduced modulo `p^N`.

    void padic_poly_randtest_val(padic_poly_t f, flint_rand_t state, slong val, slong len, const padic_ctx_t ctx)
    # Sets `f` to a random polynomial of length at most ``len``
    # with at most the prescribed valuation ``val`` and entries
    # reduced modulo `p^N`.
    # Specifically, we aim to set the valuation to be exactly equal
    # to ``val``, but do not check for additional cancellation
    # when creating the coefficients.

    void padic_poly_set_padic(padic_poly_t poly, const padic_t x, const padic_ctx_t ctx)
    # Sets the polynomial ``poly`` to the `p`-adic number `x`,
    # reduced to the precision of the polynomial.

    void padic_poly_set(padic_poly_t poly1, const padic_poly_t poly2, const padic_ctx_t ctx)
    # Sets the polynomial ``poly1`` to the polynomial ``poly2``,
    # reduced to the precision of ``poly1``.

    void padic_poly_set_si(padic_poly_t poly, slong x, const padic_ctx_t ctx)
    # Sets the polynomial ``poly`` to the ``signed slong``
    # integer `x` reduced to the precision of the polynomial.

    void padic_poly_set_ui(padic_poly_t poly, ulong x, const padic_ctx_t ctx)
    # Sets the polynomial ``poly`` to the ``unsigned slong``
    # integer `x` reduced to the precision of the polynomial.

    void padic_poly_set_fmpz(padic_poly_t poly, const fmpz_t x, const padic_ctx_t ctx)
    # Sets the polynomial ``poly`` to the integer `x`
    # reduced to the precision of the polynomial.

    void padic_poly_set_fmpq(padic_poly_t poly, const fmpq_t x, const padic_ctx_t ctx)
    # Sets the polynomial ``poly`` to the value of the rational `x`,
    # reduced to the precision of the polynomial.

    void padic_poly_set_fmpz_poly(padic_poly_t rop, const fmpz_poly_t op, const padic_ctx_t ctx)
    # Sets the polynomial ``rop`` to the integer polynomial ``op``
    # reduced to the precision of the polynomial.

    void padic_poly_set_fmpq_poly(padic_poly_t rop, const fmpq_poly_t op, const padic_ctx_t ctx)
    # Sets the polynomial ``rop`` to the value of the rational
    # polynomial ``op``, reduced to the precision of the polynomial.

    int padic_poly_get_fmpz_poly(fmpz_poly_t rop, const padic_poly_t op, const padic_ctx_t ctx)
    # Sets the integer polynomial ``rop`` to the value of the `p`-adic
    # polynomial ``op`` and returns `1` if the polynomial is `p`-adically
    # integral.  Otherwise, returns `0`.

    void padic_poly_get_fmpq_poly(fmpq_poly_t rop, const padic_poly_t op, const padic_ctx_t ctx)
    # Sets ``rop`` to the rational polynomial corresponding to
    # the `p`-adic polynomial ``op``.

    void padic_poly_zero(padic_poly_t poly)
    # Sets ``poly`` to the zero polynomial.

    void padic_poly_one(padic_poly_t poly)
    # Sets ``poly`` to the constant polynomial `1`,
    # reduced to the precision of the polynomial.

    void padic_poly_swap(padic_poly_t poly1, padic_poly_t poly2)
    # Swaps the two polynomials ``poly1`` and ``poly2``,
    # including their precisions.
    # This is done efficiently by swapping pointers.

    void padic_poly_get_coeff_padic(padic_t c, const padic_poly_t poly, slong n, const padic_ctx_t ctx)
    # Sets `c` to the coefficient of `x^n` in the polynomial,
    # reduced modulo the precision of `c`.

    void padic_poly_set_coeff_padic(padic_poly_t f, slong n, const padic_t c, const padic_ctx_t ctx)
    # Sets the coefficient of `x^n` in the polynomial `f` to `c`,
    # reduced to the precision of the polynomial `f`.
    # Note that this operation can take linear time in the length
    # of the polynomial.

    bint padic_poly_equal(const padic_poly_t poly1, const padic_poly_t poly2)
    # Returns whether the two polynomials ``poly1`` and ``poly2``
    # are equal.

    bint padic_poly_is_zero(const padic_poly_t poly)
    # Returns whether the polynomial ``poly`` is the zero polynomial.

    bint padic_poly_is_one(const padic_poly_t poly)
    # Returns whether the polynomial ``poly`` is equal
    # to the constant polynomial~`1`, taking the precision
    # of the polynomial into account.

    void _padic_poly_add(fmpz *rop, slong *rval, slong N, const fmpz *op1, slong val1, slong len1, slong N1, const fmpz *op2, slong val2, slong len2, slong N2, const padic_ctx_t ctx)
    # Sets ``(rop, *val, FLINT_MAX(len1, len2)`` to the sum of
    # ``(op1, val1, len1)`` and ``(op2, val2, len2)``.
    # Assumes that the input is reduced and guarantees that this is
    # also the case for the output.
    # Assumes that `\min\{v_1, v_2\} < N`.
    # Supports aliasing between the output and input arguments.

    void padic_poly_add(padic_poly_t f, const padic_poly_t g, const padic_poly_t h, const padic_ctx_t ctx)
    # Sets `f` to the sum `g + h`.

    void _padic_poly_sub(fmpz *rop, slong *rval, slong N, const fmpz *op1, slong val1, slong len1, slong N1, const fmpz *op2, slong val2, slong len2, slong N2, const padic_ctx_t ctx)
    # Sets ``(rop, *val, FLINT_MAX(len1, len2)`` to the difference of
    # ``(op1, val1, len1)`` and ``(op2, val2, len2)``.
    # Assumes that the input is reduced and guarantees that this is
    # also the case for the output.
    # Assumes that `\min\{v_1, v_2\} < N`.
    # Support aliasing between the output and input arguments.

    void padic_poly_sub(padic_poly_t f, const padic_poly_t g, const padic_poly_t h, const padic_ctx_t ctx)
    # Sets `f` to the difference `g - h`.

    void padic_poly_neg(padic_poly_t f, const padic_poly_t g, const padic_ctx_t ctx)
    # Sets `f` to `-g`.

    void _padic_poly_scalar_mul_padic(fmpz *rop, slong *rval, slong N, const fmpz *op, slong val, slong len, const padic_t c, const padic_ctx_t ctx)
    # Sets ``(rop, *rval, len)`` to ``(op, val, len)`` multiplied
    # by the scalar `c`.
    # The result will only be correctly reduced if the polynomial
    # is non-zero.  Otherwise, the array ``(rop, len)`` will be
    # set to zero but the valuation ``*rval`` might be wrong.

    void padic_poly_scalar_mul_padic(padic_poly_t rop, const padic_poly_t op, const padic_t c, const padic_ctx_t ctx)
    # Sets the polynomial ``rop`` to the product of the
    # polynomial ``op`` and the `p`-adic number `c`,
    # reducing the result modulo `p^N`.

    void _padic_poly_mul(fmpz *rop, slong *rval, slong N, const fmpz *op1, slong val1, slong len1, const fmpz *op2, slong val2, slong len2, const padic_ctx_t ctx)
    # Sets ``(rop, *rval, len1 + len2 - 1)`` to the product of
    # ``(op1, val1, len1)`` and ``(op2, val2, len2)``.
    # Assumes that the resulting valuation ``*rval``, which is
    # the sum of the valuations ``val1`` and ``val2``, is less
    # than the precision~`N` of the context.
    # Assumes that ``len1 >= len2 > 0``.

    void padic_poly_mul(padic_poly_t res, const padic_poly_t poly1, const padic_poly_t poly2, const padic_ctx_t ctx)
    # Sets the polynomial ``res`` to the product of the two polynomials
    # ``poly1`` and ``poly2``, reduced modulo `p^N`.

    void _padic_poly_pow(fmpz *rop, slong *rval, slong N, const fmpz *op, slong val, slong len, ulong e, const padic_ctx_t ctx)
    # Sets the polynomial ``(rop, *rval, e (len - 1) + 1)`` to the
    # polynomial ``(op, val, len)`` raised to the power~`e`.
    # Assumes that `e > 1` and ``len > 0``.
    # Does not support aliasing between the input and output arguments.

    void padic_poly_pow(padic_poly_t rop, const padic_poly_t op, ulong e, const padic_ctx_t ctx)
    # Sets the polynomial ``rop`` to the polynomial ``op`` raised
    # to the power~`e`, reduced to the precision in ``rop``.
    # In the special case `e = 0`, sets ``rop`` to the constant
    # polynomial one reduced to the precision of ``rop``.
    # Also note that when `e = 1`, this operation sets ``rop`` to
    # ``op`` and then reduces ``rop``.
    # When the valuation of the input polynomial is negative,
    # this results in a loss of `p`-adic precision.  Suppose
    # that the input polynomial is given to precision~`N` and
    # has valuation~`v < 0`.  The result then has valuation
    # `e v < 0` but is only correct to precision `N + (e - 1) v`.

    void padic_poly_inv_series(padic_poly_t g, const padic_poly_t f, slong n, const padic_ctx_t ctx)
    # Computes the power series inverse `g` of `f` modulo `X^n`,
    # where `n \geq 1`.
    # Given the polynomial `f \in \mathbf{Q}[X] \subset \mathbf{Q}_p[X]`,
    # there exists a unique polynomial `f^{-1} \in \mathbf{Q}[X]` such that
    # `f f^{-1} = 1` modulo `X^n`.  This function sets `g` to `f^{-1}`
    # reduced modulo `p^N`.
    # Assumes that the constant coefficient of `f` is non-zero.
    # Moreover, assumes that the valuation of the constant coefficient
    # of `f` is minimal among the coefficients of `f`.
    # Note that the result `g` is zero if and only if  `- \operatorname{ord}_p(f) \geq N`.

    void _padic_poly_derivative(fmpz *rop, slong *rval, slong N, const fmpz *op, slong val, slong len, const padic_ctx_t ctx)
    # Sets ``(rop, rval)`` to the derivative of ``(op, val)`` reduced
    # modulo `p^N`.
    # Supports aliasing of the input and the output parameters.

    void padic_poly_derivative(padic_poly_t rop, const padic_poly_t op, const padic_ctx_t ctx)
    # Sets ``rop`` to the derivative of ``op``, reducing the
    # result modulo the precision of ``rop``.

    void padic_poly_shift_left(padic_poly_t rop, const padic_poly_t op, slong n, const padic_ctx_t ctx)
    # Notationally, sets the polynomial ``rop`` to the polynomial ``op``
    # multiplied by `x^n`, where `n \geq 0`, and reduces the result.

    void padic_poly_shift_right(padic_poly_t rop, const padic_poly_t op, slong n, const padic_ctx_t ctx)
    # Notationally, sets the polynomial ``rop`` to the polynomial
    # ``op`` after floor division by `x^n`, where `n \geq 0`, ensuring
    # the result is reduced.

    void _padic_poly_evaluate_padic(fmpz_t u, slong *v, slong N, const fmpz *poly, slong val, slong len, const fmpz_t a, slong b, const padic_ctx_t ctx)
    void padic_poly_evaluate_padic(padic_t y, const padic_poly_t poly, const padic_t a, const padic_ctx_t ctx)
    # Sets the `p`-adic number ``y`` to ``poly`` evaluated at `a`,
    # reduced in the given context.
    # Suppose that the polynomial can be written as `F(X) = p^w f(X)`
    # with `\operatorname{ord}_p(f) = 1`, that `\operatorname{ord}_p(a) = b` and that both are
    # defined to precision~`N`.  Then `f` is defined to precision
    # `N-w` and so `f(a)` is defined to precision `N-w` when `a` is
    # integral and `N-w+(n-1)b` when `b < 0`, where `n = \deg(f)`.  Thus,
    # `y = F(a)` is defined to precision `N` when `a` is integral and
    # `N+(n-1)b` when `b < 0`.

    void _padic_poly_compose(fmpz *rop, slong *rval, slong N, const fmpz *op1, slong val1, slong len1, const fmpz *op2, slong val2, slong len2, const padic_ctx_t ctx)
    # Sets ``(rop, *rval, (len1-1)*(len2-1)+1)`` to the composition
    # of the two input polynomials, reducing the result modulo `p^N`.
    # Assumes that ``len1`` is non-zero.
    # Does not support aliasing.

    void padic_poly_compose(padic_poly_t rop, const padic_poly_t op1, const padic_poly_t op2, const padic_ctx_t ctx)
    # Sets ``rop`` to the composition of ``op1`` and ``op2``,
    # reducing the result in the given context.
    # To be clear about the order of composition, let `f(X)` and `g(X)`
    # denote the polynomials ``op1`` and ``op2``, respectively.
    # Then ``rop`` is set to `f(g(X))`.

    void _padic_poly_compose_pow(fmpz *rop, slong *rval, slong N, const fmpz *op, slong val, slong len, slong k, const padic_ctx_t ctx)
    # Sets ``(rop, *rval, (len - 1)*k + 1)`` to the composition of
    # ``(op, val, len)`` and the monomial `x^k`, where `k \geq 1`.
    # Assumes that ``len`` is positive.
    # Supports aliasing between the input and output polynomials.

    void padic_poly_compose_pow(padic_poly_t rop, const padic_poly_t op, slong k, const padic_ctx_t ctx)
    # Sets ``rop`` to the composition of ``op`` and the monomial `x^k`,
    # where `k \geq 1`.
    # Note that no reduction takes place.

    int padic_poly_debug(const padic_poly_t poly)
    # Prints the data defining the `p`-adic polynomial ``poly``
    # in a simple format useful for debugging purposes.
    # In the current implementation, always returns `1`.

    int _padic_poly_fprint(FILE *file, const fmpz *poly, slong val, slong len, const padic_ctx_t ctx)
    int padic_poly_fprint(FILE *file, const padic_poly_t poly, const padic_ctx_t ctx)
    # Prints a simple representation of the polynomial ``poly``
    # to the stream ``file``.
    # A non-zero polynomial is represented by the number of coefficients,
    # two spaces, followed by a list of the coefficients, which are printed
    # in a way depending on the print mode,
    # In the ``PADIC_TERSE`` mode, the coefficients are printed as
    # rational numbers.
    # The ``PADIC_SERIES`` mode is currently not supported and will
    # raise an abort signal.
    # In the ``PADIC_VAL_UNIT`` mode, the coefficients are printed
    # in the form `p^v u`.
    # The zero polynomial is represented by ``"0"``.
    # In the current implementation, always returns `1`.

    int _padic_poly_print(const fmpz *poly, slong val, slong len, const padic_ctx_t ctx)
    int padic_poly_print(const padic_poly_t poly, const padic_ctx_t ctx)
    # Prints a simple representation of the polynomial ``poly``
    # to ``stdout``.
    # In the current implementation, always returns `1`.

    int _padic_poly_fprint_pretty(FILE *file, const fmpz *poly, slong val, slong len, const char *var, const padic_ctx_t ctx)
    int padic_poly_fprint_pretty(FILE *file, const padic_poly_t poly, const char *var, const padic_ctx_t ctx)
    int _padic_poly_print_pretty(const fmpz *poly, slong val, slong len, const char *var, const padic_ctx_t ctx)
    int padic_poly_print_pretty(const padic_poly_t poly, const char *var, const padic_ctx_t ctx)

    bint _padic_poly_is_canonical(const fmpz *op, slong val, slong len, const padic_ctx_t ctx)
    bint padic_poly_is_canonical(const padic_poly_t op, const padic_ctx_t ctx)
    bint _padic_poly_is_reduced(const fmpz *op, slong val, slong len, slong N, const padic_ctx_t ctx)
    bint padic_poly_is_reduced(const padic_poly_t op, const padic_ctx_t ctx)
