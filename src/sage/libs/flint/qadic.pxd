# distutils: libraries = flint
# distutils: depends = flint/qadic.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void qadic_ctx_init(qadic_ctx_t ctx, const fmpz_t p, long d, long min, long max, const char *var, padic_print_mode mode)
    # Initialises the context ``ctx`` with prime `p`, extension degree `d`,
    # variable name ``var`` and printing mode ``mode``. The defining polynomial
    # is chosen as a Conway polynomial if possible and otherwise as a random
    # sparse polynomial.
    # Stores powers of `p` with exponents between ``min`` (inclusive) and
    # ``max`` exclusive.  Assumes that ``min`` is at most ``max``.
    # Assumes that `p` is a prime.
    # Assumes that the string ``var`` is a null-terminated string
    # of length at least one.
    # Assumes that the printing mode is one of ``PADIC_TERSE``,
    # ``PADIC_SERIES``, or ``PADIC_VAL_UNIT``.
    # This function also carries out some relevant precomputation for
    # arithmetic in `\mathbf{Q}_p / (p^N)` such as powers of `p` close
    # to `p^N`.

    void qadic_ctx_init_conway(qadic_ctx_t ctx, const fmpz_t p, long d, long min, long max, const char *var, padic_print_mode mode)
    # Initialises the context ``ctx`` with prime `p`, extension degree `d`,
    # variable name ``var`` and printing mode ``mode``. The defining polynomial
    # is chosen as a Conway polynomial, hence has restrictions on the
    # prime and the degree.
    # Stores powers of `p` with exponents between ``min`` (inclusive) and
    # ``max`` exclusive.  Assumes that ``min`` is at most ``max``.
    # Assumes that `p` is a prime.
    # Assumes that the string ``var`` is a null-terminated string
    # of length at least one.
    # Assumes that the printing mode is one of ``PADIC_TERSE``,
    # ``PADIC_SERIES``, or ``PADIC_VAL_UNIT``.
    # This function also carries out some relevant precomputation for
    # arithmetic in `\mathbf{Q}_p / (p^N)` such as powers of `p` close
    # to `p^N`.

    void qadic_ctx_clear(qadic_ctx_t ctx)
    # Clears all memory that has been allocated as part of the context.

    long qadic_ctx_degree(const qadic_ctx_t ctx)
    # Returns the extension degree.

    void qadic_ctx_print(const qadic_ctx_t ctx)
    # Prints the data from the given context.

    void qadic_init(qadic_t rop)
    # Initialises the element ``rop``, setting its value to `0`.

    void qadic_init2(qadic_t rop, long prec)
    # Initialises the element ``rop`` with the given output precision,
    # setting the value to `0`.

    void qadic_clear(qadic_t rop)
    # Clears the element ``rop``.

    void _fmpz_poly_reduce(fmpz *R, long lenR, const fmpz *a, const long *j, long len)
    # Reduces a polynomial ``(R, lenR)`` modulo a sparse monic
    # polynomial `f(X) = \sum_{i} a_{i} X^{j_{i}}` of degree at
    # least `2`.
    # Assumes that the array `j` of positive length ``len`` is
    # sorted in ascending order.
    # Allows zero-padding in ``(R, lenR)``.

    void _fmpz_mod_poly_reduce(fmpz *R, long lenR, const fmpz *a, const long *j, long len, const fmpz_t p)
    # Reduces a polynomial ``(R, lenR)`` modulo a sparse monic
    # polynomial `f(X) = \sum_{i} a_{i} X^{j_{i}}` of degree at
    # least `2` in `\mathbf{Z}/(p)`, where `p` is typically a prime
    # power.
    # Assumes that the array `j` of positive length ``len`` is
    # sorted in ascending order.
    # Allows zero-padding in ``(R, lenR)``.

    void qadic_reduce(qadic_t rop, const qadic_ctx_t ctx)
    # Reduces ``rop`` modulo `f(X)` and `p^N`.

    long qadic_val(const qadic_t op)
    # Returns the valuation of ``op``.

    long qadic_prec(const qadic_t op)
    # Returns the precision of ``op``.

    void qadic_randtest(qadic_t rop, flint_rand_t state, const qadic_ctx_t ctx)
    # Generates a random element of `\mathbf{Q}_q`.

    void qadic_randtest_not_zero(qadic_t rop, flint_rand_t state, const qadic_ctx_t ctx)
    # Generates a random non-zero element of `\mathbf{Q}_q`.

    void qadic_randtest_val(qadic_t rop, flint_rand_t state, long v, const qadic_ctx_t ctx)
    # Generates a random element of `\mathbf{Q}_q` with prescribed
    # valuation ``val``.
    # Note that if `v \geq N` then the element is necessarily zero.

    void qadic_randtest_int(qadic_t rop, flint_rand_t state, const qadic_ctx_t ctx)
    # Generates a random element of `\mathbf{Q}_q` with non-negative valuation.

    void qadic_set(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    # Sets ``rop`` to ``op``.

    void qadic_zero(qadic_t rop)
    # Sets ``rop`` to zero.

    void qadic_one(qadic_t rop)
    # Sets ``rop`` to one, reduced in the given context.
    # Note that if the precision `N` is non-positive then ``rop``
    # is actually set to zero.

    void qadic_gen(qadic_t rop, const qadic_ctx_t ctx)
    # Sets ``rop`` to the generator `X` for the extension
    # when `N > 0`, and zero otherwise.  If the extension degree
    # is one, raises an abort signal.

    void qadic_set_ui(qadic_t rop, unsigned long op, const qadic_ctx_t ctx)
    # Sets ``rop`` to the integer ``op``, reduced in the
    # context.

    int qadic_get_padic(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    # If the element ``op`` lies in `\mathbf{Q}_p`, sets ``rop``
    # to its value and returns `1`;  otherwise, returns `0`.

    int qadic_is_zero(const qadic_t op)
    # Returns whether ``op`` is equal to zero.

    int qadic_is_one(const qadic_t op)
    # Returns whether ``op`` is equal to one in the given
    # context.

    int qadic_equal(const qadic_t op1, const qadic_t op2)
    # Returns whether ``op1`` and ``op2`` are equal.

    void qadic_add(qadic_t rop, const qadic_t op1, const qadic_t op2, const qadic_ctx_t ctx)
    # Sets ``rop`` to the sum of ``op1`` and ``op2``.
    # Assumes that both ``op1`` and ``op2`` are reduced in the
    # given context and ensures that ``rop`` is, too.

    void qadic_sub(qadic_t rop, const qadic_t op1, const qadic_t op2, const qadic_ctx_t ctx)
    # Sets ``rop`` to the difference of ``op1`` and ``op2``.
    # Assumes that both ``op1`` and ``op2`` are reduced in the
    # given context and ensures that ``rop`` is, too.

    void qadic_neg(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    # Sets ``rop`` to the negative of ``op``.
    # Assumes that ``op`` is reduced in the given context and
    # ensures that ``rop`` is, too.

    void qadic_mul(qadic_t rop, const qadic_t op1, const qadic_t op2, const qadic_ctx_t ctx)
    # Sets ``rop`` to the product of ``op1`` and ``op2``,
    # reducing the output in the given context.

    void _qadic_inv(fmpz *rop, const fmpz *op, long len, const fmpz *a, const long *j, long lena, const fmpz_t p, long N)
    # Sets ``(rop, d)`` to the inverse of ``(op, len)``
    # modulo `f(X)` given by ``(a,j,lena)`` and `p^N`.
    # Assumes that ``(op,len)`` has valuation `0`, that is,
    # that it represents a `p`-adic unit.
    # Assumes that ``len`` is at most `d`.
    # Does not support aliasing.

    void qadic_inv(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    # Sets ``rop`` to the inverse of ``op``, reduced in the given context.

    void _qadic_pow(fmpz *rop, const fmpz *op, long len, const fmpz_t e, const fmpz *a, const long *j, long lena, const fmpz_t p)
    # Sets ``(rop, 2*d-1)`` to ``(op,len)`` raised to the power `e`,
    # reduced modulo `f(X)` given by ``(a, j, lena)`` and `p`, which
    # is expected to be a prime power.
    # Assumes that `e \geq 0` and that ``len`` is positive and at most `d`.
    # Although we require that ``rop`` provides space for
    # `2d - 1` coefficients, the output will be reduced modulo
    # `f(X)`, which is a polynomial of degree `d`.
    # Does not support aliasing.

    void qadic_pow(qadic_t rop, const qadic_t op, const fmpz_t e, const qadic_ctx_t ctx)
    # Sets ``rop`` the ``op`` raised to the power `e`.
    # Currently assumes that `e \geq 0`.
    # Note that for any input ``op``, ``rop`` is set to one in the
    # given context whenever `e = 0`.

    int qadic_sqrt(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    # Return ``1`` if the input is a square (to input precision). If so, set
    # ``rop`` to a square root (truncated to output precision).

    void _qadic_exp_rectangular(fmpz *rop, const fmpz *op, long v, long len, const fmpz *a, const long *j, long lena, const fmpz_t p, long N, const fmpz_t pN)
    # Sets ``(rop, 2*d - 1)`` to the exponential of ``(op, v, len)``
    # reduced modulo `p^N`, assuming that the series converges.
    # Assumes that ``(op, v, len)`` is non-zero.
    # Does not support aliasing.

    int qadic_exp_rectangular(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    # Returns whether the exponential series converges at ``op``
    # and sets ``rop`` to its value reduced modulo in the given
    # context.

    void _qadic_exp_balanced(fmpz *rop, const fmpz *x, long v, long len, const fmpz *a, const long *j, long lena, const fmpz_t p, long N, const fmpz_t pN)
    # Sets ``(rop, d)`` to the exponential of ``(op, v, len)``
    # reduced modulo `p^N`, assuming that the series converges.
    # Assumes that ``len`` is in `[1,d)` but supports zero padding,
    # including the special case when ``(op, len)`` is zero.
    # Supports aliasing between ``rop`` and ``op``.

    int qadic_exp_balanced(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    # Returns whether the exponential series converges at ``op``
    # and sets ``rop`` to its value reduced modulo in the given
    # context.

    void _qadic_exp(fmpz *rop, const fmpz *op, long v, long len, const fmpz *a, const long *j, long lena, const fmpz_t p, long N, const fmpz_t pN)
    # Sets ``(rop, 2*d - 1)`` to the exponential of ``(op, v, len)``
    # reduced modulo `p^N`, assuming that the series converges.
    # Assumes that ``(op, v, len)`` is non-zero.
    # Does not support aliasing.

    int qadic_exp(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    # Returns whether the exponential series converges at ``op``
    # and sets ``rop`` to its value reduced modulo in the given
    # context.
    # The exponential series converges if the valuation of ``op``
    # is at least `2` or `1` when `p` is even or odd, respectively.

    void _qadic_log_rectangular(fmpz *z, const fmpz *y, long v, long len, const fmpz *a, const long *j, long lena, const fmpz_t p, long N, const fmpz_t pN)
    # Computes
    # .. math ::
    # z = - \sum_{i = 1}^{\infty} \frac{y^i}{i} \pmod{p^N}.
    # Note that this can be used to compute the `p`-adic logarithm
    # via the equation
    # .. math ::
    # \log(x) & = \sum_{i=1}^{\infty} (-1)^{i-1} \frac{(x-1)^i}{i} \\
    # & = - \sum_{i=1}^{\infty} \frac{(1-x)^i}{i}.
    # Assumes that `y = 1 - x` is non-zero and that `v = \operatorname{ord}_p(y)`
    # is at least `1` when `p` is odd and at least `2` when `p = 2`
    # so that the series converges.
    # Assumes that `y` is reduced modulo `p^N`.
    # Assumes that `v < N`, and in particular `N \geq 2`.
    # Supports aliasing between `y` and `z`.

    int qadic_log_rectangular(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    # Returns whether the `p`-adic logarithm function converges at
    # ``op``, and if so sets ``rop`` to its value.

    void _qadic_log_balanced(fmpz *z, const fmpz *y, long len, const fmpz *a, const long *j, long lena, const fmpz_t p, long N, const fmpz_t pN)
    # Computes `(z, d)` as
    # .. math ::
    # z = - \sum_{i = 1}^{\infty} \frac{y^i}{i} \pmod{p^N}.
    # Assumes that `v = \operatorname{ord}_p(y)` is at least `1` when `p` is odd and
    # at least `2` when `p = 2` so that the series converges.
    # Supports aliasing between `z` and `y`.

    int qadic_log_balanced(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    # Returns whether the `p`-adic logarithm function converges at
    # ``op``, and if so sets ``rop`` to its value.

    void _qadic_log(fmpz *z, const fmpz *y, long v, long len, const fmpz *a, const long *j, long lena, const fmpz_t p, long N, const fmpz_t pN)
    # Computes `(z, d)` as
    # .. math ::
    # z = - \sum_{i = 1}^{\infty} \frac{y^i}{i} \pmod{p^N}.
    # Note that this can be used to compute the `p`-adic logarithm
    # via the equation
    # .. math ::
    # \log(x) & = \sum_{i=1}^{\infty} (-1)^{i-1} \frac{(x-1)^i}{i} \\
    # & = - \sum_{i=1}^{\infty} \frac{(1-x)^i}{i}.
    # Assumes that `y = 1 - x` is non-zero and that `v = \operatorname{ord}_p(y)`
    # is at least `1` when `p` is odd and at least `2` when `p = 2`
    # so that the series converges.
    # Assumes that `(y, d)` is reduced modulo `p^N`.
    # Assumes that `v < N`, and hence in particular `N \geq 2`.
    # Supports aliasing between `z` and `y`.

    int qadic_log(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    # Returns whether the `p`-adic logarithm function converges at
    # ``op``, and if so sets ``rop`` to its value.
    # The `p`-adic logarithm function is defined by the usual series
    # .. math ::
    # \log_p(x) = \sum_{i=1}^{\infty} (-1)^{i-1} \frac{(x-1)^i}{i}
    # but this only converges when `\operatorname{ord}_p(x)` is at least `2` or `1`
    # when `p = 2` or `p > 2`, respectively.

    void _qadic_frobenius_a(fmpz *rop, long e, const fmpz *a, const long *j, long lena, const fmpz_t p, long N)
    # Computes `\sigma^e(X) \bmod{p^N}` where `X` is such that
    # `\mathbf{Q}_q \cong \mathbf{Q}_p[X]/(f(X))`.
    # Assumes that the precision `N` is at least `2` and that the
    # extension is non-trivial, i.e. `d \geq 2`.
    # Assumes that `0 < e < d`.
    # Sets ``(rop, 2*d-1)``, although the actual length of the
    # output will be at most `d`.

    void _qadic_frobenius(fmpz *rop, const fmpz *op, long len, long e, const fmpz *a, const long *j, long lena, const fmpz_t p, long N)
    # Sets ``(rop, 2*d-1)`` to `\Sigma` evaluated at ``(op, len)``.
    # Assumes that ``len`` is positive but at most `d`.
    # Assumes that `0 < e < d`.
    # Does not support aliasing.

    void qadic_frobenius(qadic_t rop, const qadic_t op, long e, const qadic_ctx_t ctx)
    # Evaluates the homomorphism `\Sigma^e` at ``op``.
    # Recall that `\mathbf{Q}_q / \mathbf{Q}_p` is Galois with Galois group
    # `\langle \Sigma \rangle \cong \langle \sigma \rangle`, which is also
    # isomorphic to `\mathbf{Z}/d\mathbf{Z}`, where
    # `\sigma \in \operatorname{Gal}(\mathbf{F}_q/\mathbf{F}_p)` is the Frobenius element
    # `\sigma \colon x \mapsto x^p` and `\Sigma` is its lift to
    # `\operatorname{Gal}(\mathbf{Q}_q/\mathbf{Q}_p)`.
    # This functionality is implemented as ``GaloisImage()`` in Magma.

    void _qadic_teichmuller(fmpz *rop, const fmpz *op, long len, const fmpz *a, const long *j, long lena, const fmpz_t p, long N)
    # Sets ``(rop, d)`` to the Teichmüller lift of ``(op, len)``
    # modulo `p^N`.
    # Does not support aliasing.

    void qadic_teichmuller(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    # Sets ``rop`` to the Teichmüller lift of ``op`` to the
    # precision given in the context.
    # For a unit ``op``, this is the unique `(q-1)`\th root of unity
    # which is congruent to ``op`` modulo `p`.
    # Sets ``rop`` to zero if ``op`` is zero in the given context.
    # Raises an exception if the valuation of ``op`` is negative.

    void _qadic_trace(fmpz_t rop, const fmpz *op, long len, const fmpz *a, const long *j, long lena, const fmpz_t pN)
    void qadic_trace(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    # Sets ``rop`` to the trace of ``op``.
    # For an element `a \in \mathbf{Q}_q`, multiplication by `a` defines
    # a `\mathbf{Q}_p`-linear map on `\mathbf{Q}_q`.  We define the trace
    # of `a` as the trace of this map.  Equivalently, if `\Sigma` generates
    # `\operatorname{Gal}(\mathbf{Q}_q / \mathbf{Q}_p)` then the trace of `a` is equal to
    # `\sum_{i=0}^{d-1} \Sigma^i (a)`.

    void _qadic_norm(fmpz_t rop, const fmpz *op, long len, const fmpz *a, const long *j, long lena, const fmpz_t p, long N)
    # Sets ``rop`` to the norm of the element ``(op,len)``
    # in `\mathbf{Z}_q` to precision `N`, where ``len`` is at
    # least one.
    # The result will be reduced modulo `p^N`.
    # Note that whenever ``(op,len)`` is a unit, so is its norm.
    # Thus, the output ``rop`` of this function will typically
    # not have to be canonicalised or reduced by the caller.

    void qadic_norm(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    # Computes the norm of ``op`` to the given precision.
    # Algorithm selection is automatic depending on the input.

    void qadic_norm_analytic(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    # Whenever ``op`` has valuation greater than `(p-1)^{-1}`, this
    # routine computes its norm ``rop`` via
    # .. math ::
    # \operatorname{Norm} (x) = \exp \Bigl( \bigl( \operatorname{Trace} \log (x) \bigr) \Bigr).
    # In the special case that ``op`` lies in `\mathbf{Q}_p`, returns
    # its norm as `\operatorname{Norm}(x) = x^d`, where `d` is the extension degree.
    # Otherwise, raises an ``abort`` signal.
    # The complexity of this implementation is quasi-linear in `d` and `N`,
    # and polynomial in `\log p`.

    void qadic_norm_resultant(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    # Sets ``rop`` to the norm of ``op``, using the formula
    # .. math ::
    # \operatorname{Norm}(x) = \ell(f)^{-\deg(a)} \operatorname{Res}(f(X), a(X)),
    # where `\mathbf{Q}_q \cong \mathbf{Q}_p[X] / (f(X))`, `\ell(f)` is the
    # leading coefficient of `f(X)`, and `a(X) \in \mathbf{Q}_p[X]` denotes
    # the same polynomial as `x`.
    # The complexity of the current implementation is given by
    # `\mathcal{O}(d^4 M(N \log p))`, where `M(n)` denotes the
    # complexity of multiplying to `n`-bit integers.

    int qadic_fprint_pretty(FILE *file, const qadic_t op, const qadic_ctx_t ctx)
    # Prints a pretty representation of ``op`` to ``file``.
    # In the current implementation, always returns `1`.  The return code is
    # part of the function's signature to allow for a later implementation to
    # return the number of characters printed or a non-positive error code.

    int qadic_print_pretty(const qadic_t op, const qadic_ctx_t ctx)
    # Prints a pretty representation of ``op`` to ``stdout``.
    # In the current implementation, always returns `1`.  The return code is
    # part of the function's signature to allow for a later implementation to
    # return the number of characters printed or a non-positive error code.
