# distutils: libraries = flint
# distutils: depends = flint/padic.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    fmpz * padic_unit(const padic_t op)
    # Returns the unit part of the `p`-adic number as a FLINT integer, which
    # can be used as an operand for the ``fmpz`` functions.

    long padic_val(const padic_t op)
    # Returns the valuation part of the `p`-adic number.
    # Note that this function is implemented as a macro and that
    # the expression ``padic_val(op)`` can be used as both an
    # *lvalue* and an *rvalue*.

    long padic_get_val(const padic_t op)
    # Returns the valuation part of the `p`-adic number.

    long padic_prec(const padic_t op)
    # Returns the precision of the `p`-adic number.
    # Note that this function is implemented as a macro and that
    # the expression ``padic_prec(op)`` can be used as both an
    # *lvalue* and an *rvalue*.

    long padic_get_prec(const padic_t op)
    # Returns the precision of the `p`-adic number.

    void padic_ctx_init(padic_ctx_t ctx, const fmpz_t p, long min, long max, padic_print_mode mode)
    # Initialises the context ``ctx`` with the given data.
    # Assumes that `p` is a prime.  This is not verified but the subsequent
    # behaviour is undefined if `p` is a composite number.
    # Assumes that ``min`` and ``max`` are non-negative and that
    # ``min`` is at most ``max``, raising an ``abort`` signal
    # otherwise.
    # Assumes that the printing mode is one of ``PADIC_TERSE``,
    # ``PADIC_SERIES``, or ``PADIC_VAL_UNIT``.  Using the example
    # `x = 7^{-1} 12` in `\mathbf{Q}_7`, these behave as follows:
    # In ``PADIC_TERSE`` mode, a `p`-adic number is printed
    # in the same way as a rational number, e.g. ``12/7``.
    # In ``PADIC_SERIES`` mode, a `p`-adic number is printed
    # digit by digit, e.g. ``5*7^-1 + 1``.
    # In ``PADIC_VAL_UNIT`` mode, a `p`-adic number is
    # printed showing the valuation and unit parts separately,
    # e.g. ``12*7^-1``.

    void padic_ctx_clear(padic_ctx_t ctx)
    # Clears all memory that has been allocated as part of the context.

    int _padic_ctx_pow_ui(fmpz_t rop, unsigned long e, const padic_ctx_t ctx)
    # Sets ``rop`` to `p^e` as efficiently as possible, where
    # ``rop`` is expected to be an uninitialised ``fmpz_t``.
    # If the return value is non-zero, it is the responsibility of
    # the caller to clear the returned integer.

    void padic_init(padic_t rop)
    # Initialises the `p`-adic number with the precision set to
    # ``PADIC_DEFAULT_PREC``, which is defined as `20`.

    void padic_init2(padic_t rop, long N)
    # Initialises the `p`-adic number ``rop`` with precision `N`.

    void padic_clear(padic_t rop)
    # Clears all memory used by the `p`-adic number ``rop``.

    void _padic_canonicalise(padic_t rop, const padic_ctx_t ctx)
    # Brings the `p`-adic number ``rop`` into canonical form.
    # That is to say, ensures that either `u = v = 0` or
    # `p \nmid u`.  There is no reduction modulo a power
    # of `p`.

    void _padic_reduce(padic_t rop, const padic_ctx_t ctx)
    # Given a `p`-adic number ``rop`` in canonical form,
    # reduces it modulo `p^N`.

    void padic_reduce(padic_t rop, const padic_ctx_t ctx)
    # Ensures that the `p`-adic number ``rop`` is reduced.

    void padic_randtest(padic_t rop, flint_rand_t state, const padic_ctx_t ctx)
    # Sets ``rop`` to a random `p`-adic number modulo `p^N` with valuation
    # in the range `[- \lceil N/10\rceil, N)`, `[N - \lceil -N/10\rceil, N)`, or `[-10, 0)`
    # as `N` is positive, negative or zero, whenever ``rop`` is non-zero.

    void padic_randtest_not_zero(padic_t rop, flint_rand_t state, const padic_ctx_t ctx)
    # Sets ``rop`` to a random non-zero `p`-adic number modulo `p^N`,
    # where the range of the valuation is as for the function
    # :func:`padic_randtest`.

    void padic_randtest_int(padic_t rop, flint_rand_t state, const padic_ctx_t ctx)
    # Sets ``rop`` to a random `p`-adic integer modulo `p^N`.
    # Note that whenever `N \leq 0`, ``rop`` is set to zero.

    void padic_set(padic_t rop, const padic_t op, const padic_ctx_t ctx)
    # Sets ``rop`` to the `p`-adic number ``op``.

    void padic_set_si(padic_t rop, long op, const padic_ctx_t ctx)
    # Sets the `p`-adic number ``rop`` to the
    # ``slong`` integer ``op``.

    void padic_set_ui(padic_t rop, unsigned long op, const padic_ctx_t ctx)
    # Sets the `p`-adic number ``rop`` to the ``ulong``
    # integer ``op``.

    void padic_set_fmpz(padic_t rop, const fmpz_t op, const padic_ctx_t ctx)
    # Sets the `p`-adic number ``rop`` to the integer ``op``.

    void padic_set_fmpq(padic_t rop, const fmpq_t op, const padic_ctx_t ctx)
    # Sets ``rop`` to the rational ``op``.

    void padic_set_mpz(padic_t rop, const mpz_t op, const padic_ctx_t ctx)
    # Sets the `p`-adic number ``rop`` to the MPIR integer ``op``.

    void padic_set_mpq(padic_t rop, const mpq_t op, const padic_ctx_t ctx)
    # Sets ``rop`` to the MPIR rational ``op``.

    void padic_get_fmpz(fmpz_t rop, const padic_t op, const padic_ctx_t ctx)
    # Sets the integer ``rop`` to the exact `p`-adic integer ``op``.
    # If ``op`` is not a `p`-adic integer, raises an ``abort`` signal.

    void padic_get_fmpq(fmpq_t rop, const padic_t op, const padic_ctx_t ctx)
    # Sets the rational ``rop`` to the `p`-adic number ``op``.

    void padic_get_mpz(mpz_t rop, const padic_t op, const padic_ctx_t ctx)
    # Sets the MPIR integer ``rop`` to the `p`-adic integer ``op``.
    # If ``op`` is not a `p`-adic integer, raises an ``abort`` signal.

    void padic_get_mpq(mpq_t rop, const padic_t op, const padic_ctx_t ctx)
    # Sets the MPIR rational ``rop`` to the value of ``op``.

    void padic_swap(padic_t op1, padic_t op2)
    # Swaps the two `p`-adic numbers ``op1`` and ``op2``.
    # Note that this includes swapping the precisions.  In particular, this
    # operation is not equivalent to swapping ``op1`` and ``op2``
    # using :func:`padic_set` and an auxiliary variable whenever the
    # precisions of the two elements are different.

    void padic_zero(padic_t rop)
    # Sets the `p`-adic number ``rop`` to zero.

    void padic_one(padic_t rop)
    # Sets the `p`-adic number ``rop`` to one, reduced modulo the
    # precision of ``rop``.

    int padic_is_zero(const padic_t op)
    # Returns whether ``op`` is equal to zero.

    int padic_is_one(const padic_t op)
    # Returns whether ``op`` is equal to one, that is, whether
    # `u = 1` and `v = 0`.

    int padic_equal(const padic_t op1, const padic_t op2)
    # Returns whether ``op1`` and ``op2`` are equal, that is,
    # whether `u_1 = u_2` and `v_1 = v_2`.

    long * _padic_lifts_exps(long *n, long N)
    # Given a positive integer `N` define the sequence
    # `a_0 = N, a_1 = \lceil a_0/2\rceil, \dotsc, a_{n-1} = \lceil a_{n-2}/2\rceil = 1`.
    # Then `n = \lceil\log_2 N\rceil + 1`.
    # This function sets `n` and allocates and returns the array `a`.

    void _padic_lifts_pows(fmpz *pow, const long *a, long n, const fmpz_t p)
    # Given an array `a` as computed above, this function
    # computes the corresponding powers of `p`, that is,
    # ``pow[i]`` is equal to `p^{a_i}`.

    void padic_add(padic_t rop, const padic_t op1, const padic_t op2, const padic_ctx_t ctx)
    # Sets ``rop`` to the sum of ``op1`` and ``op2``.

    void padic_sub(padic_t rop, const padic_t op1, const padic_t op2, const padic_ctx_t ctx)
    # Sets ``rop`` to the difference of ``op1`` and ``op2``.

    void padic_neg(padic_t rop, const padic_t op, const padic_ctx_t ctx)
    # Sets ``rop`` to the additive inverse of ``op``.

    void padic_mul(padic_t rop, const padic_t op1, const padic_t op2, const padic_ctx_t ctx)
    # Sets ``rop`` to the product of ``op1`` and ``op2``.

    void padic_shift(padic_t rop, const padic_t op, long v, const padic_ctx_t ctx)
    # Sets ``rop`` to the product of ``op`` and `p^v`.

    void padic_div(padic_t rop, const padic_t op1, const padic_t op2, const padic_ctx_t ctx)
    # Sets ``rop`` to the quotient of ``op1`` and ``op2``.

    void _padic_inv_precompute(padic_inv_t S, const fmpz_t p, long N)
    # Pre-computes some data and allocates temporary space for
    # `p`-adic inversion using Hensel lifting.

    void _padic_inv_clear(padic_inv_t S)
    # Frees the memory used by `S`.

    void _padic_inv_precomp(fmpz_t rop, const fmpz_t op, const padic_inv_t S)
    # Sets ``rop`` to the inverse of ``op`` modulo `p^N`,
    # assuming that ``op`` is a unit and `N \geq 1`.
    # In the current implementation, allows aliasing, but this might
    # change in future versions.
    # Uses some data `S` precomputed by calling the function
    # :func:`_padic_inv_precompute`.  Note that this object
    # is not declared ``const`` and in fact it carries a field
    # providing temporary work space.  This allows repeated calls of
    # this function to avoid repeated memory allocations, as used
    # e.g. by the function :func:`padic_log`.

    void _padic_inv(fmpz_t rop, const fmpz_t op, const fmpz_t p, long N)
    # Sets ``rop`` to the inverse of ``op`` modulo `p^N`,
    # assuming that ``op`` is a unit and `N \geq 1`.
    # In the current implementation, allows aliasing, but this might
    # change in future versions.

    void padic_inv(padic_t rop, const padic_t op, const padic_ctx_t ctx)
    # Computes the inverse of ``op`` modulo `p^N`.
    # Suppose that ``op`` is given as `x = u p^v`.
    # Raises an ``abort`` signal if `v < -N`.  Otherwise,
    # computes the inverse of `u` modulo `p^{N+v}`.
    # This function employs Hensel lifting of an inverse modulo `p`.

    int padic_sqrt(padic_t rop, const padic_t op, const padic_ctx_t ctx)
    # Returns whether ``op`` is a `p`-adic square.  If this is
    # the case, sets ``rop`` to one of the square roots;  otherwise,
    # the value of ``rop`` is undefined.
    # We have the following theorem:
    # Let `u \in \mathbf{Z}^{\times}`.  Then `u` is a
    # square if and only if `u \bmod p` is a square in
    # `\mathbf{Z} / p \mathbf{Z}`, for `p > 2`, or if
    # `u \bmod 8` is a square in `\mathbf{Z} / 8 \mathbf{Z}`,
    # for `p = 2`.

    void padic_pow_si(padic_t rop, const padic_t op, long e, const padic_ctx_t ctx)
    # Sets ``rop`` to ``op`` raised to the power `e`,
    # which is defined as one whenever `e = 0`.
    # Assumes that some computations involving `e` and the
    # valuation of ``op`` do not overflow in the ``slong``
    # range.
    # Note that if the input `x = p^v u` is defined modulo `p^N`
    # then `x^e = p^{ev} u^e` is defined modulo `p^{N + (e - 1) v}`,
    # which is a precision loss in case `v < 0`.

    long _padic_exp_bound(long v, long N, const fmpz_t p)
    # Returns an integer `i` such that for all `j \geq i` we have
    # `\operatorname{ord}_p(x^j / j!) \geq N`, where `\operatorname{ord}_p(x) = v`.
    # When `p` is a word-sized prime,
    # returns `\left\lceil \frac{(p-1)N - 1}{(p-1)v - 1}\right\rceil`.
    # Otherwise, returns `\lceil N/v\rceil`.
    # Assumes that `v < N`.  Moreover, `v` has to be at least `2` or `1`,
    # depending on whether `p` is `2` or odd.

    void _padic_exp_rectangular(fmpz_t rop, const fmpz_t u, long v, const fmpz_t p, long N)
    void _padic_exp_balanced(fmpz_t rop, const fmpz_t u, long v, const fmpz_t p, long N)
    void _padic_exp(fmpz_t rop, const fmpz_t u, long v, const fmpz_t p, long N)
    # Sets ``rop`` to the `p`-exponential function evaluated at
    # `x = p^v u`, reduced modulo `p^N`.
    # Assumes that `x \neq 0`, that `\operatorname{ord}_p(x) < N` and that
    # `\exp(x)` converges, that is, that `\operatorname{ord}_p(x)` is at least
    # `2` or `1` depending on whether the prime `p` is `2` or odd.
    # Supports aliasing between ``rop`` and `u`.

    int padic_exp(padic_t y, const padic_t x, const padic_ctx_t ctx)
    # Returns whether the `p`-adic exponential function converges at
    # the `p`-adic number `x`, and if so sets `y` to its value.
    # The `p`-adic exponential function is defined by the usual series
    # .. math ::
    # \exp_p(x) = \sum_{i = 0}^{\infty} \frac{x^i}{i!}
    # but this only converges only when `\operatorname{ord}_p(x) > 1 / (p - 1)`.  For
    # elements `x \in \mathbf{Q}_p`, this means that `\operatorname{ord}_p(x) \geq 1`
    # when `p \geq 3` and `\operatorname{ord}_2(x) \geq 2` when `p = 2`.

    int padic_exp_rectangular(padic_t y, const padic_t x, const padic_ctx_t ctx)
    # Returns whether the `p`-adic exponential function converges at
    # the `p`-adic number `x`, and if so sets `y` to its value.
    # Uses a rectangular splitting algorithm to evaluate the series
    # expression of `\exp(x) \bmod{p^N}`.

    int padic_exp_balanced(padic_t y, const padic_t x, const padic_ctx_t ctx)
    # Returns whether the `p`-adic exponential function converges at
    # the `p`-adic number `x`, and if so sets `y` to its value.
    # Uses a balanced approach, balancing the size of chunks of `x`
    # with the valuation and hence the rate of convergence, which
    # results in a quasi-linear algorithm in `N`, for fixed `p`.

    long _padic_log_bound(long v, long N, const fmpz_t p)
    # Returns `b` such that for all `i \geq b` we have
    # .. math ::
    # i v - \operatorname{ord}_p(i) \geq N
    # where `v \geq 1`.
    # Assumes that `1 \leq v < N` or `2 \leq v < N` when `p` is
    # odd or `p = 2`, respectively, and also that `N < 2^{f-2}`
    # where `f` is ``FLINT_BITS``.

    void _padic_log(fmpz_t z, const fmpz_t y, long v, const fmpz_t p, long N)
    void _padic_log_rectangular(fmpz_t z, const fmpz_t y, long v, const fmpz_t p, long N)
    void _padic_log_satoh(fmpz_t z, const fmpz_t y, long v, const fmpz_t p, long N)
    void _padic_log_balanced(fmpz_t z, const fmpz_t y, long v, const fmpz_t p, long N)
    # Computes
    # .. math ::
    # z = - \sum_{i = 1}^{\infty} \frac{y^i}{i} \pmod{p^N},
    # reduced modulo `p^N`.
    # Note that this can be used to compute the `p`-adic logarithm
    # via the equation
    # .. math ::
    # \log(x) & = \sum_{i=1}^{\infty} (-1)^{i-1} \frac{(x-1)^i}{i} \\
    # & = - \sum_{i=1}^{\infty} \frac{(1-x)^i}{i}.
    # Assumes that `y = 1 - x` is non-zero and that `v = \operatorname{ord}_p(y)`
    # is at least `1` when `p` is odd and at least `2` when `p = 2`
    # so that the series converges.
    # Assumes that `v < N`, and hence in particular `N \geq 2`.
    # Does not support aliasing between `y` and `z`.

    int padic_log(padic_t rop, const padic_t op, const padic_ctx_t ctx)
    # Returns whether the `p`-adic logarithm function converges at
    # the `p`-adic number ``op``, and if so sets ``rop`` to its
    # value.
    # The `p`-adic logarithm function is defined by the usual series
    # .. math ::
    # \log_p(x) = \sum_{i=1}^{\infty} (-1)^{i-1} \frac{(x-1)^i}{i}
    # but this only converges when `\operatorname{ord}_p(x - 1)` is at least `2`
    # or `1` when `p = 2` or `p > 2`, respectively.

    int padic_log_rectangular(padic_t rop, const padic_t op, const padic_ctx_t ctx)
    # Returns whether the `p`-adic logarithm function converges at
    # the `p`-adic number ``op``, and if so sets ``rop`` to its
    # value.
    # Uses a rectangular splitting algorithm to evaluate the series
    # expression of `\log(x) \bmod{p^N}`.

    int padic_log_satoh(padic_t rop, const padic_t op, const padic_ctx_t ctx)
    # Returns whether the `p`-adic logarithm function converges at
    # the `p`-adic number ``op``, and if so sets ``rop`` to its
    # value.
    # Uses an algorithm based on a result of Satoh, Skjernaa and Taguchi
    # that `\operatorname{ord}_p\bigl(a^{p^k} - 1\bigr) > k`, which implies that
    # .. math ::
    # \log(a) \equiv p^{-k} \Bigl( \log\bigl(a^{p^k}\bigr) \pmod{p^{N+k}}
    # \Bigr) \pmod{p^N}.

    int padic_log_balanced(padic_t rop, const padic_t op, const padic_ctx_t ctx)
    # Returns whether the `p`-adic logarithm function converges at
    # the `p`-adic number ``op``, and if so sets ``rop`` to its
    # value.

    void _padic_teichmuller(fmpz_t rop, const fmpz_t op, const fmpz_t p, long N)
    # Computes the Teichm\"uller lift of the `p`-adic unit ``op``,
    # assuming that `N \geq 1`.
    # Supports aliasing between ``rop`` and ``op``.

    void padic_teichmuller(padic_t rop, const padic_t op, const padic_ctx_t ctx)
    # Computes the Teichm\"uller lift of the `p`-adic unit ``op``.
    # If ``op`` is a `p`-adic integer divisible by `p`, sets ``rop``
    # to zero, which satisfies `t^p - t = 0`, although it is clearly not
    # a `(p-1)`-st root of unity.
    # If ``op`` has negative valuation, raises an ``abort`` signal.

    unsigned long padic_val_fac_ui_2(unsigned long n)
    # Computes the `2`-adic valuation of `n!`.
    # Note that since `n` fits into an ``ulong``, so does
    # `\operatorname{ord}_2(n!)` since `\operatorname{ord}_2(n!) \leq (n - 1) / (p - 1) = n - 1`.

    unsigned long padic_val_fac_ui(unsigned long n, const fmpz_t p)
    # Computes the `p`-adic valuation of `n!`.
    # Note that since `n` fits into an ``ulong``, so does
    # `\operatorname{ord}_p(n!)` since `\operatorname{ord}_p(n!) \leq (n - 1) / (p - 1)`.

    void padic_val_fac(fmpz_t rop, const fmpz_t op, const fmpz_t p)
    # Sets ``rop`` to the `p`-adic valuation of the factorial
    # of ``op``, assuming that ``op`` is non-negative.

    char * padic_get_str(char * str, const padic_t op, const padic_ctx_t ctx)
    # Returns the string representation of the `p`-adic number ``op``
    # according to the printing mode set in the context.
    # If ``str`` is ``NULL`` then a new block of memory is allocated
    # and a pointer to this is returned.  Otherwise, it is assumed that
    # the string ``str`` is large enough to hold the representation and
    # it is also the return value.

    int _padic_fprint(FILE * file, const fmpz_t u, long v, const padic_ctx_t ctx)
    int padic_fprint(FILE * file, const padic_t op, const padic_ctx_t ctx)
    # Prints the string representation of the `p`-adic number ``op``
    # to the stream ``file``.
    # In the current implementation, always returns `1`.

    int _padic_print(const fmpz_t u, long v, const padic_ctx_t ctx)
    int padic_print(const padic_t op, const padic_ctx_t ctx)
    # Prints the string representation of the `p`-adic number ``op``
    # to the stream ``stdout``.
    # In the current implementation, always returns `1`.

    void padic_debug(const padic_t op)
    # Prints debug information about ``op`` to the stream ``stdout``,
    # in the format ``"(u v N)"``.
