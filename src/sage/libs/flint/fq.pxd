# distutils: libraries = flint
# distutils: depends = flint/fq.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void fq_ctx_init(fq_ctx_t ctx, const fmpz_t p, slong d, const char *var)
    # Initialises the context for prime `p` and extension degree `d`,
    # with name ``var`` for the generator.  By default, it will try
    # use a Conway polynomial; if one is not available, a random
    # irreducible polynomial will be used.
    # Assumes that `p` is a prime.
    # Assumes that the string ``var`` is a null-terminated string
    # of length at least one.

    int _fq_ctx_init_conway(fq_ctx_t ctx, const fmpz_t p, slong d, const char *var)
    # Attempts to initialise the context for prime `p` and extension
    # degree `d`, with name ``var`` for the generator using a Conway
    # polynomial for the modulus.
    # Returns `1` if the Conway polynomial is in the database for the
    # given size and the initialization is successful; otherwise,
    # returns `0`.
    # Assumes that `p` is a prime.
    # Assumes that the string ``var`` is a null-terminated string
    # of length at least one.

    void fq_ctx_init_conway(fq_ctx_t ctx, const fmpz_t p, slong d, const char *var)
    # Initialises the context for prime `p` and extension degree `d`,
    # with name ``var`` for the generator using a Conway polynomial
    # for the modulus.
    # Assumes that `p` is a prime.
    # Assumes that the string ``var`` is a null-terminated string
    # of length at least one.

    void fq_ctx_init_modulus(fq_ctx_t ctx, const fmpz_mod_poly_t modulus, const fmpz_mod_ctx_t ctxp, const char *var)
    # Initialises the context for given ``modulus`` with name
    # ``var`` for the generator.
    # Assumes that ``modulus`` is an irreducible polynomial over the finite field `\mathbf{F}_{p}` in ``ctxp``.
    # Assumes that the string ``var`` is a null-terminated string
    # of length at least one.

    void fq_ctx_clear(fq_ctx_t ctx)
    # Clears all memory that has been allocated as part of the context.

    const fmpz_mod_poly_struct* fq_ctx_modulus(const fq_ctx_t ctx)
    # Returns a pointer to the modulus in the context.

    slong fq_ctx_degree(const fq_ctx_t ctx)
    # Returns the degree of the field extension
    # `[\mathbf{F}_{q} : \mathbf{F}_{p}]`, which
    # is equal to `\log_{p} q`.

    const fmpz * fq_ctx_prime(const fq_ctx_t ctx)
    # Returns a pointer to the prime `p` in the context.

    void fq_ctx_order(fmpz_t f, const fq_ctx_t ctx)
    # Sets `f` to be the size of the finite field.

    int fq_ctx_fprint(FILE * file, const fq_ctx_t ctx)
    # Prints the context information to ``file``. Returns 1 for a
    # success and a negative number for an error.

    void fq_ctx_print(const fq_ctx_t ctx)
    # Prints the context information to ``stdout``.

    void fq_ctx_randtest(fq_ctx_t ctx, flint_rand_t state)
    # Initializes ``ctx`` to a random finite field.  Assumes that
    # ``fq_ctx_init`` has not been called on ``ctx`` already.

    void fq_ctx_randtest_reducible(fq_ctx_t ctx, flint_rand_t state)
    # Initializes ``ctx`` to a random extension of a prime field.
    # The modulus may or may not be irreducible.  Assumes that
    # ``fq_ctx_init`` has not been called on ``ctx`` already.

    void fq_init(fq_t rop, const fq_ctx_t ctx)
    # Initialises the element ``rop``, setting its value to `0`.

    void fq_init2(fq_t rop, const fq_ctx_t ctx)
    # Initialises ``poly`` with at least enough space for it to be an element
    # of ``ctx`` and sets it to `0`.

    void fq_clear(fq_t rop, const fq_ctx_t ctx)
    # Clears the element ``rop``.

    void _fq_sparse_reduce(fmpz *R, slong lenR, const fq_ctx_t ctx)
    # Reduces ``(R, lenR)`` modulo the polynomial `f` given by the
    # modulus of ``ctx``.

    void _fq_dense_reduce(fmpz *R, slong lenR, const fq_ctx_t ctx)
    # Reduces ``(R, lenR)`` modulo the polynomial `f` given by the
    # modulus of ``ctx`` using Newton division.

    void _fq_reduce(fmpz *r, slong lenR, const fq_ctx_t ctx)
    # Reduces ``(R, lenR)`` modulo the polynomial `f` given by the
    # modulus of ``ctx``.  Does either sparse or dense reduction
    # based on ``ctx->sparse_modulus``.

    void fq_reduce(fq_t rop, const fq_ctx_t ctx)
    # Reduces the polynomial ``rop`` as an element of
    # `\mathbf{F}_p[X] / (f(X))`.

    void fq_add(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
    # Sets ``rop`` to the sum of ``op1`` and ``op2``.

    void fq_sub(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
    # Sets ``rop`` to the difference of ``op1`` and ``op2``.

    void fq_sub_one(fq_t rop, const fq_t op1, const fq_ctx_t ctx)
    # Sets ``rop`` to the difference of ``op1`` and `1`.

    void fq_neg(fq_t rop, const fq_t op, const fq_ctx_t ctx)
    # Sets ``rop`` to the negative of ``op``.

    void fq_mul(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
    # Sets ``rop`` to the product of ``op1`` and ``op2``,
    # reducing the output in the given context.

    void fq_mul_fmpz(fq_t rop, const fq_t op, const fmpz_t x, const fq_ctx_t ctx)
    # Sets ``rop`` to the product of ``op`` and `x`,
    # reducing the output in the given context.

    void fq_mul_si(fq_t rop, const fq_t op, slong x, const fq_ctx_t ctx)
    # Sets ``rop`` to the product of ``op`` and `x`,
    # reducing the output in the given context.

    void fq_mul_ui(fq_t rop, const fq_t op, ulong x, const fq_ctx_t ctx)
    # Sets ``rop`` to the product of ``op`` and `x`,
    # reducing the output in the given context.

    void fq_sqr(fq_t rop, const fq_t op, const fq_ctx_t ctx)
    # Sets ``rop`` to the square of ``op``,
    # reducing the output in the given context.

    void fq_div(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
    # Sets ``rop`` to the quotient of ``op1`` and ``op2``,
    # reducing the output in the given context.

    void _fq_inv(fmpz *rop, const fmpz *op, slong len, const fq_ctx_t ctx)
    # Sets ``(rop, d)`` to the inverse of the non-zero element
    # ``(op, len)``.

    void fq_inv(fq_t rop, const fq_t op, const fq_ctx_t ctx)
    # Sets ``rop`` to the inverse of the non-zero element ``op``.

    void fq_gcdinv(fq_t f, fq_t inv, const fq_t op, const fq_ctx_t ctx)
    # Sets ``inv`` to be the inverse of ``op`` modulo the modulus
    # of ``ctx``.  If ``op`` is not invertible, then ``f`` is
    # set to a factor of the modulus; otherwise, it is set to one.

    void _fq_pow(fmpz *rop, const fmpz *op, slong len, const fmpz_t e, const fq_ctx_t ctx)
    # Sets ``(rop, 2*d-1)`` to ``(op,len)`` raised to the power `e`,
    # reduced modulo `f(X)`, the modulus of ``ctx``.
    # Assumes that `e \geq 0` and that ``len`` is positive and at most `d`.
    # Although we require that ``rop`` provides space for
    # `2d - 1` coefficients, the output will be reduced modulo
    # `f(X)`, which is a polynomial of degree `d`.
    # Does not support aliasing.

    void fq_pow(fq_t rop, const fq_t op, const fmpz_t e, const fq_ctx_t ctx)
    # Sets ``rop`` the ``op`` raised to the power `e`.
    # Currently assumes that `e \geq 0`.
    # Note that for any input ``op``, ``rop`` is set to `1`
    # whenever `e = 0`.

    void fq_pow_ui(fq_t rop, const fq_t op, const ulong e, const fq_ctx_t ctx)
    # Sets ``rop`` the ``op`` raised to the power `e`.
    # Currently assumes that `e \geq 0`.
    # Note that for any input ``op``, ``rop`` is set to `1`
    # whenever `e = 0`.

    int fq_sqrt(fq_t rop, const fq_t op1, const fq_ctx_t ctx)
    # Sets ``rop`` to the square root of ``op1`` if it is a square, and return
    # `1`, otherwise return `0`.

    void fq_pth_root(fq_t rop, const fq_t op1, const fq_ctx_t ctx)
    # Sets ``rop`` to a `p^{th}` root root of ``op1``.  Currently,
    # this computes the root by raising ``op1`` to `p^{d-1}` where
    # `d` is the degree of the extension.

    int fq_is_square(const fq_t op, const fq_ctx_t ctx)
    # Return ``1`` if ``op`` is a square.

    int fq_fprint_pretty(FILE *file, const fq_t op, const fq_ctx_t ctx)
    # Prints a pretty representation of ``op`` to ``file``.
    # In the current implementation, always returns `1`.  The return code is
    # part of the function's signature to allow for a later implementation to
    # return the number of characters printed or a non-positive error code.

    int fq_print_pretty(const fq_t op, const fq_ctx_t ctx)
    # Prints a pretty representation of ``op`` to ``stdout``.
    # In the current implementation, always returns `1`.  The return code is
    # part of the function's signature to allow for a later implementation to
    # return the number of characters printed or a non-positive error code.

    int fq_fprint(FILE * file, const fq_t op, const fq_ctx_t ctx)
    # Prints a representation of ``op`` to ``file``.
    # For further details on the representation used, see
    # :func:`fmpz_mod_poly_fprint`.

    void fq_print(const fq_t op, const fq_ctx_t ctx)
    # Prints a representation of ``op`` to ``stdout``.
    # For further details on the representation used, see
    # :func:`fmpz_mod_poly_print`.

    char * fq_get_str(const fq_t op, const fq_ctx_t ctx)
    # Returns the plain FLINT string representation of the element
    # ``op``.

    char * fq_get_str_pretty(const fq_t op, const fq_ctx_t ctx)
    # Returns a pretty representation of the element ``op`` using the
    # null-terminated string ``x`` as the variable name.

    void fq_randtest(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)
    # Generates a random element of `\mathbf{F}_q`.

    void fq_randtest_not_zero(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)
    # Generates a random non-zero element of `\mathbf{F}_q`.

    void fq_randtest_dense(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)
    # Generates a random element of `\mathbf{F}_q` which has an
    # underlying polynomial with dense coefficients.

    void fq_rand(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)
    # Generates a high quality random element of `\mathbf{F}_q`.

    void fq_rand_not_zero(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)
    # Generates a high quality non-zero random element of `\mathbf{F}_q`.

    void fq_set(fq_t rop, const fq_t op, const fq_ctx_t ctx)
    # Sets ``rop`` to ``op``.

    void fq_set_si(fq_t rop, const slong x, const fq_ctx_t ctx)
    # Sets ``rop`` to ``x``, considered as an element of
    # `\mathbf{F}_p`.

    void fq_set_ui(fq_t rop, const ulong x, const fq_ctx_t ctx)
    # Sets ``rop`` to ``x``, considered as an element of
    # `\mathbf{F}_p`.

    void fq_set_fmpz(fq_t rop, const fmpz_t x, const fq_ctx_t ctx)
    # Sets ``rop`` to ``x``, considered as an element of
    # `\mathbf{F}_p`.

    void fq_swap(fq_t op1, fq_t op2, const fq_ctx_t ctx)
    # Swaps the two elements ``op1`` and ``op2``.

    void fq_zero(fq_t rop, const fq_ctx_t ctx)
    # Sets ``rop`` to zero.

    void fq_one(fq_t rop, const fq_ctx_t ctx)
    # Sets ``rop`` to one, reduced in the given context.

    void fq_gen(fq_t rop, const fq_ctx_t ctx)
    # Sets ``rop`` to a generator for the finite field.
    # There is no guarantee this is a multiplicative generator of
    # the finite field.

    int fq_get_fmpz(fmpz_t rop, const fq_t op, const fq_ctx_t ctx)
    # If ``op`` has a lift to the integers, return `1` and set ``rop`` to the lift in `[0,p)`.
    # Otherwise, return `0` and leave `rop` undefined.

    void fq_get_fmpz_poly(fmpz_poly_t a, const fq_t b, const fq_ctx_t ctx)

    void fq_get_fmpz_mod_poly(fmpz_mod_poly_t a, const fq_t b, const fq_ctx_t ctx)
    # Set ``a`` to a representative of ``b`` in ``ctx``.
    # The representatives are taken in `(\mathbb{Z}/p\mathbb{Z})[x]/h(x)` where `h(x)` is the defining polynomial in ``ctx``.

    void fq_set_fmpz_poly(fq_t a, const fmpz_poly_t b, const fq_ctx_t ctx)

    void fq_set_fmpz_mod_poly(fq_t a, const fmpz_mod_poly_t b, const fq_ctx_t ctx)
    # Set ``a`` to the element in ``ctx`` with representative ``b``.
    # The representatives are taken in `(\mathbb{Z}/p\mathbb{Z})[x]/h(x)` where `h(x)` is the defining polynomial in ``ctx``.

    void fq_get_fmpz_mod_mat(fmpz_mod_mat_t col, const fq_t a, const fq_ctx_t ctx)
    # Convert ``a`` to a column vector of length ``degree(ctx)``.

    void fq_set_fmpz_mod_mat(fq_t a, const fmpz_mod_mat_t col, const fq_ctx_t ctx)
    # Convert a column vector ``col`` of length ``degree(ctx)`` to an element of ``ctx``.

    int fq_is_zero(const fq_t op, const fq_ctx_t ctx)
    # Returns whether ``op`` is equal to zero.

    int fq_is_one(const fq_t op, const fq_ctx_t ctx)
    # Returns whether ``op`` is equal to one.

    int fq_equal(const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
    # Returns whether ``op1`` and ``op2`` are equal.

    int fq_is_invertible(const fq_t op, const fq_ctx_t ctx)
    # Returns whether ``op`` is an invertible element.

    int fq_is_invertible_f(fq_t f, const fq_t op, const fq_ctx_t ctx)
    # Returns whether ``op`` is an invertible element.  If it is not,
    # then ``f`` is set of a factor of the modulus.

    void _fq_trace(fmpz_t rop, const fmpz *op, slong len, const fq_ctx_t ctx)
    # Sets ``rop`` to the trace of the non-zero element ``(op, len)``
    # in `\mathbf{F}_{q}`.

    void fq_trace(fmpz_t rop, const fq_t op, const fq_ctx_t ctx)
    # Sets ``rop`` to the trace of ``op``.
    # For an element `a \in \mathbf{F}_q`, multiplication by `a` defines
    # a `\mathbf{F}_p`-linear map on `\mathbf{F}_q`.  We define the
    # trace of `a` as the trace of this map.  Equivalently, if `\Sigma`
    # generates `\operatorname{Gal}(\mathbf{F}_q / \mathbf{F}_p)` then the trace of
    # `a` is equal to `\sum_{i=0}^{d-1} \Sigma^i (a)`, where `d =
    # \log_{p} q`.

    void _fq_norm(fmpz_t rop, const fmpz *op, slong len, const fq_ctx_t ctx)
    # Sets ``rop`` to the norm of the non-zero element ``(op, len)``
    # in `\mathbf{F}_{q}`.

    void fq_norm(fmpz_t rop, const fq_t op, const fq_ctx_t ctx)
    # Computes the norm of ``op``.
    # For an element `a \in \mathbf{F}_q`, multiplication by `a` defines
    # a `\mathbf{F}_p`-linear map on `\mathbf{F}_q`.  We define the norm
    # of `a` as the determinant of this map.  Equivalently, if `\Sigma` generates
    # `\operatorname{Gal}(\mathbf{F}_q / \mathbf{F}_p)` then the trace of `a` is equal to
    # `\prod_{i=0}^{d-1} \Sigma^i (a)`, where
    # `d = \text{dim}_{\mathbf{F}_p}(\mathbf{F}_q)`.
    # Algorithm selection is automatic depending on the input.

    void _fq_frobenius(fmpz *rop, const fmpz *op, slong len, slong e, const fq_ctx_t ctx)
    # Sets ``(rop, 2d-1)`` to the image of ``(op, len)`` under the
    # Frobenius operator raised to the e-th power, assuming that neither
    # ``op`` nor ``e`` are zero.

    void fq_frobenius(fq_t rop, const fq_t op, slong e, const fq_ctx_t ctx)
    # Evaluates the homomorphism `\Sigma^e` at ``op``.
    # Recall that `\mathbf{F}_q / \mathbf{F}_p` is Galois with Galois group
    # `\langle \sigma \rangle`, which is also isomorphic to
    # `\mathbf{Z}/d\mathbf{Z}`, where
    # `\sigma \in \operatorname{Gal}(\mathbf{F}_q/\mathbf{F}_p)` is the Frobenius element
    # `\sigma \colon x \mapsto x^p`.

    int fq_multiplicative_order(fmpz * ord, const fq_t op, const fq_ctx_t ctx)
    # Computes the order of ``op`` as an element of the
    # multiplicative group of ``ctx``.
    # Returns 0 if ``op`` is 0, otherwise it returns 1 if ``op``
    # is a generator of the multiplicative group, and -1 if it is not.
    # This function can also be used to check primitivity of a generator of
    # a finite field whose defining polynomial is not primitive.

    int fq_is_primitive(const fq_t op, const fq_ctx_t ctx)
    # Returns whether ``op`` is primitive, i.e., whether it is a
    # generator of the multiplicative group of ``ctx``.

    void fq_bit_pack(fmpz_t f, const fq_t op, flint_bitcnt_t bit_size, const fq_ctx_t ctx)
    # Packs ``op`` into bitfields of size ``bit_size``, writing the
    # result to ``f``.

    void fq_bit_unpack(fq_t rop, const fmpz_t f, flint_bitcnt_t bit_size, const fq_ctx_t ctx)
    # Unpacks into ``rop`` the element with coefficients packed into
    # fields of size ``bit_size`` as represented by the integer
    # ``f``.
