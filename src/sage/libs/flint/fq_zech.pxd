# distutils: libraries = flint
# distutils: depends = flint/fq_zech.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void fq_zech_ctx_init(fq_zech_ctx_t ctx, const fmpz_t p, long d, const char *var)
    # Initialises the context for prime `p` and extension degree `d`,
    # with name ``var`` for the generator.  By default, it will try
    # use a Conway polynomial; if one is not available, a random
    # primitive polynomial will be used.
    # Assumes that `p` is a prime and :math:`p^d < 2^{\mathtt{FLINT\_BITS}}`.
    # Assumes that the string ``var`` is a null-terminated string
    # of length at least one.

    int _fq_zech_ctx_init_conway(fq_zech_ctx_t ctx, const fmpz_t p, long d, const char *var)
    # Attempts to initialise the context for prime `p` and extension
    # degree `d`, with name ``var`` for the generator using a Conway
    # polynomial for the modulus.
    # Returns `1` if the Conway polynomial is in the database for the
    # given size and the initialization is successful; otherwise,
    # returns `0`.
    # Assumes that `p` is a prime and `p^d < 2^\mathtt{FLINT\_BITS}`.
    # Assumes that the string ``var`` is a null-terminated string
    # of length at least one.

    void fq_zech_ctx_init_conway(fq_zech_ctx_t ctx, const fmpz_t p, long d, const char *var)
    # Initialises the context for prime `p` and extension degree `d`,
    # with name ``var`` for the generator using a Conway polynomial
    # for the modulus.
    # Assumes that `p` is a prime and `p^d < 2^\mathtt{FLINT\_BITS}`.
    # Assumes that the string ``var`` is a null-terminated string
    # of length at least one.

    void fq_zech_ctx_init_random(fq_zech_ctx_t ctx, const fmpz_t p, long d, const char *var)
    # Initialises the context for prime `p` and extension degree `d`,
    # with name ``var`` for the generator using a random primitive
    # polynomial.
    # Assumes that `p` is a prime and `p^d < 2^\mathtt{FLINT\_BITS}`.
    # Assumes that the string ``var`` is a null-terminated string
    # of length at least one.

    void fq_zech_ctx_init_modulus(fq_zech_ctx_t ctx, const nmod_poly_t modulus, const char *var)
    # Initialises the context for given ``modulus`` with name
    # ``var`` for the generator.
    # Assumes that ``modulus`` is an primitive polynomial over
    # `\mathbf{F}_{p}`. An exception is raised if a non-primitive modulus is
    # detected.
    # Assumes that the string ``var`` is a null-terminated string
    # of length at least one.

    int fq_zech_ctx_init_modulus_check(fq_zech_ctx_t ctx, const nmod_poly_t modulus, const char *var)
    # As per the previous function, but returns `0` if the modulus was not
    # primitive and `1` if the context was successfully initialised with the
    # given modulus. No exception is raised.

    void fq_zech_ctx_init_fq_nmod_ctx(fq_zech_ctx_t ctx, fq_nmod_ctx_t ctxn)
    # Initializes the context ``ctx`` to be the Zech representation
    # for the finite field given by ``ctxn``.

    int fq_zech_ctx_init_fq_nmod_ctx_check(fq_zech_ctx_t ctx, fq_nmod_ctx_t ctxn)
    # As per the previous function but returns `0` if a non-primitive modulus is
    # detected. Returns `0` if the Zech representation was successfully
    # initialised.

    void fq_zech_ctx_clear(fq_zech_ctx_t ctx)
    # Clears all memory that has been allocated as part of the context.

    const nmod_poly_struct* fq_zech_ctx_modulus(const fq_zech_ctx_t ctx)
    # Returns a pointer to the modulus in the context.

    long fq_zech_ctx_degree(const fq_zech_ctx_t ctx)
    # Returns the degree of the field extension
    # `[\mathbf{F}_{q} : \mathbf{F}_{p}]`, which
    # is equal to `\log_{p} q`.

    fmpz * fq_zech_ctx_prime(const fq_zech_ctx_t ctx)
    # Returns a pointer to the prime `p` in the context.

    void fq_zech_ctx_order(fmpz_t f, const fq_zech_ctx_t ctx)
    # Sets `f` to be the size of the finite field.

    mp_limb_t fq_zech_ctx_order_ui(const fq_zech_ctx_t ctx)
    # Returns the size of the finite field.

    int fq_zech_ctx_fprint(FILE * file, const fq_zech_ctx_t ctx)
    # Prints the context information to {\tt{file}}. Returns 1 for a
    # success and a negative number for an error.

    void fq_zech_ctx_print(const fq_zech_ctx_t ctx)
    # Prints the context information to {\tt{stdout}}.

    void fq_zech_ctx_randtest(fq_zech_ctx_t ctx, flint_rand_t state)
    # Initializes ``ctx`` to a random finite field.  Assumes that
    # ``fq_zech_ctx_init`` has not been called on ``ctx`` already.

    void fq_zech_ctx_randtest_reducible(fq_zech_ctx_t ctx, flint_rand_t state)
    # Since the Zech logarithm representation does not work with a
    # non-irreducible modulus, does the same as
    # ``fq_zech_ctx_randtest``.

    void fq_zech_init(fq_zech_t rop, const fq_zech_ctx_t ctx)
    # Initialises the element ``rop``, setting its value to `0`.

    void fq_zech_init2(fq_zech_t rop, const fq_zech_ctx_t ctx)
    # Initialises ``poly`` with at least enough space for it to be an element
    # of ``ctx`` and sets it to `0`.

    void fq_zech_clear(fq_zech_t rop, const fq_zech_ctx_t ctx)
    # Clears the element ``rop``.

    void _fq_zech_sparse_reduce(mp_ptr R, long lenR, const fq_zech_ctx_t ctx)
    # Reduces ``(R, lenR)`` modulo the polynomial `f` given by the
    # modulus of ``ctx``.

    void _fq_zech_dense_reduce(mp_ptr R, long lenR, const fq_zech_ctx_t ctx)
    # Reduces ``(R, lenR)`` modulo the polynomial `f` given by the
    # modulus of ``ctx`` using Newton division.

    void _fq_zech_reduce(mp_ptr r, long lenR, const fq_zech_ctx_t ctx)
    # Reduces ``(R, lenR)`` modulo the polynomial `f` given by the
    # modulus of ``ctx``.  Does either sparse or dense reduction
    # based on ``ctx->sparse_modulus``.

    void fq_zech_reduce(fq_zech_t rop, const fq_zech_ctx_t ctx)
    # Reduces the polynomial ``rop`` as an element of
    # `\mathbf{F}_p[X] / (f(X))`.

    void fq_zech_add(fq_zech_t rop, const fq_zech_t op1, const fq_zech_t op2, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to the sum of ``op1`` and ``op2``.

    void fq_zech_sub(fq_zech_t rop, const fq_zech_t op1, const fq_zech_t op2, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to the difference of ``op1`` and ``op2``.

    void fq_zech_sub_one(fq_zech_t rop, const fq_zech_t op1, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to the difference of ``op1`` and `1`.

    void fq_zech_neg(fq_zech_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to the negative of ``op``.

    void fq_zech_mul(fq_zech_t rop, const fq_zech_t op1, const fq_zech_t op2, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to the product of ``op1`` and ``op2``,
    # reducing the output in the given context.

    void fq_zech_mul_fmpz(fq_zech_t rop, const fq_zech_t op, const fmpz_t x, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to the product of ``op`` and `x`,
    # reducing the output in the given context.

    void fq_zech_mul_si(fq_zech_t rop, const fq_zech_t op, long x, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to the product of ``op`` and `x`,
    # reducing the output in the given context.

    void fq_zech_mul_ui(fq_zech_t rop, const fq_zech_t op, unsigned long x, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to the product of ``op`` and `x`,
    # reducing the output in the given context.

    void fq_zech_sqr(fq_zech_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to the square of ``op``,
    # reducing the output in the given context.

    void fq_zech_div(fq_zech_t rop, const fq_zech_t op1, const fq_zech_t op2, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to the quotient of ``op1`` and ``op2``,
    # reducing the output in the given context.

    void _fq_zech_inv(mp_ptr *rop, mp_srcptr *op, long len, const fq_zech_ctx_t ctx)
    # Sets ``(rop, d)`` to the inverse of the non-zero element
    # ``(op, len)``.

    void fq_zech_inv(fq_zech_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to the inverse of the non-zero element ``op``.

    void fq_zech_gcdinv(fq_zech_t f, fq_zech_t inv, const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Sets ``inv`` to be the inverse of ``op`` modulo the modulus
    # of ``ctx`` and sets ``f`` to one.  Since the modulus for
    # ``ctx`` is always irreducible, ``op`` is always invertible.

    void _fq_zech_pow(fmpz *rop, const fmpz *op, long len, const fmpz_t e, const fmpz * a, const long *j, long lena, const fmpz_t p)
    # Sets ``(rop, 2*d-1)`` to ``(op,len)`` raised to the power `e`,
    # reduced modulo `f(X)`, the modulus of ``ctx``.
    # Assumes that `e \geq 0` and that ``len`` is positive and at most `d`.
    # Although we require that ``rop`` provides space for
    # `2d - 1` coefficients, the output will be reduced modulo
    # `f(X)`, which is a polynomial of degree `d`.
    # Does not support aliasing.

    void fq_zech_pow(fq_zech_t rop, const fq_zech_t op, const fmpz_t e, const fq_zech_ctx_t ctx)
    # Sets ``rop`` the ``op`` raised to the power `e`.
    # Currently assumes that `e \geq 0`.
    # Note that for any input ``op``, ``rop`` is set to `1`
    # whenever `e = 0`.

    void fq_zech_pow_ui(fq_zech_t rop, const fq_zech_t op, const unsigned long e, const fq_zech_ctx_t ctx)
    # Sets ``rop`` the ``op`` raised to the power `e`.
    # Currently assumes that `e \geq 0`.
    # Note that for any input ``op``, ``rop`` is set to `1`
    # whenever `e = 0`.

    int fq_zech_sqrt(fq_zech_t rop, const fq_zech_t op1, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to the square root of ``op1`` if it is a square, and return
    # `1`, otherwise return `0`.

    void fq_zech_pth_root(fq_zech_t rop, const fq_zech_t op1, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to a `p^{th}` root root of ``op1``.  Currently,
    # this computes the root by raising ``op1`` to `p^{d-1}` where
    # `d` is the degree of the extension.

    int fq_zech_is_square(const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Return ``1`` if ``op`` is a square.

    int fq_zech_fprint_pretty(FILE *file, const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Prints a pretty representation of ``op`` to ``file``.
    # In the current implementation, always returns `1`.  The return code is
    # part of the function's signature to allow for a later implementation to
    # return the number of characters printed or a non-positive error code.

    void fq_zech_print_pretty(const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Prints a pretty representation of ``op`` to ``stdout``.
    # In the current implementation, always returns `1`.  The return code is
    # part of the function's signature to allow for a later implementation to
    # return the number of characters printed or a non-positive error code.

    int fq_zech_fprint(FILE * file, const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Prints a representation of ``op`` to ``file``.

    void fq_zech_print(const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Prints a representation of ``op`` to ``stdout``.

    char * fq_zech_get_str(const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Returns the plain FLINT string representation of the element
    # ``op``.

    char * fq_zech_get_str_pretty(const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Returns a pretty representation of the element ``op`` using the
    # null-terminated string ``x`` as the variable name.

    void fq_zech_randtest(fq_zech_t rop, flint_rand_t state, const fq_zech_ctx_t ctx)
    # Generates a random element of `\mathbf{F}_q`.

    void fq_zech_randtest_not_zero(fq_zech_t rop, flint_rand_t state, const fq_zech_ctx_t ctx)
    # Generates a random non-zero element of `\mathbf{F}_q`.

    void fq_zech_randtest_dense(fq_zech_t rop, flint_rand_t state, const fq_zech_ctx_t ctx)
    # Generates a random element of `\mathbf{F}_q` which has an
    # underlying polynomial with dense coefficients.

    void fq_zech_rand(fq_zech_t rop, flint_rand_t state, const fq_zech_ctx_t ctx)
    # Generates a high quality random element of `\mathbf{F}_q`.

    void fq_zech_rand_not_zero(fq_zech_t rop, flint_rand_t state, const fq_zech_ctx_t ctx)
    # Generates a high quality non-zero random element of `\mathbf{F}_q`.

    void fq_zech_set(fq_zech_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to ``op``.

    void fq_zech_set_si(fq_zech_t rop, const long x, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to ``x``, considered as an element of
    # `\mathbf{F}_p`.

    void fq_zech_set_ui(fq_zech_t rop, const unsigned long x, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to ``x``, considered as an element of
    # `\mathbf{F}_p`.

    void fq_zech_set_fmpz(fq_zech_t rop, const fmpz_t x, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to ``x``, considered as an element of
    # `\mathbf{F}_p`.

    void fq_zech_swap(fq_zech_t op1, fq_zech_t op2, const fq_zech_ctx_t ctx)
    # Swaps the two elements ``op1`` and ``op2``.

    void fq_zech_zero(fq_zech_t rop, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to zero.

    void fq_zech_one(fq_zech_t rop, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to one, reduced in the given context.

    void fq_zech_gen(fq_zech_t rop, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to a generator for the finite field.
    # There is no guarantee this is a multiplicative generator of
    # the finite field.

    int fq_zech_get_fmpz(fmpz_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
    # If ``op`` has a lift to the integers, return `1` and set ``rop`` to the lift in `[0,p)`.
    # Otherwise, return `0` and leave `rop` undefined.

    void fq_zech_get_fq_nmod(fq_nmod_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to the ``fq_nmod_t`` element corresponding to ``op``.

    void fq_zech_set_fq_nmod(fq_zech_t rop, const fq_nmod_t op, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to the ``fq_zech_t`` element corresponding to ``op``.

    void fq_zech_get_nmod_poly(nmod_poly_t a, const fq_zech_t b, const fq_zech_ctx_t ctx)
    # Set ``a`` to a representative of ``b`` in ``ctx``.
    # The representatives are taken in `(\mathbb{Z}/p\mathbb{Z})[x]/h(x)` where `h(x)` is the defining polynomial in ``ctx``.

    void fq_zech_set_nmod_poly(fq_zech_t a, const nmod_poly_t b, const fq_zech_ctx_t ctx)
    # Set ``a`` to the element in ``ctx`` with representative ``b``.
    # The representatives are taken in `(\mathbb{Z}/p\mathbb{Z})[x]/h(x)` where `h(x)` is the defining polynomial in ``ctx``.

    void fq_zech_get_nmod_mat(nmod_mat_t col, const fq_zech_t a, const fq_zech_ctx_t ctx)
    # Convert ``a`` to a column vector of length ``degree(ctx)``.

    void fq_zech_set_nmod_mat(fq_zech_t a, const nmod_mat_t col, const fq_zech_ctx_t ctx)
    # Convert a column vector ``col`` of length ``degree(ctx)`` to
    # an element of ``ctx``.

    int fq_zech_is_zero(const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Returns whether ``op`` is equal to zero.

    int fq_zech_is_one(const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Returns whether ``op`` is equal to one.

    int fq_zech_equal(const fq_zech_t op1, const fq_zech_t op2, const fq_zech_ctx_t ctx)
    # Returns whether ``op1`` and ``op2`` are equal.

    int fq_zech_is_invertible(const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Returns whether ``op`` is an invertible element.

    int fq_zech_is_invertible_f(fq_zech_t f, const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Returns whether ``op`` is an invertible element.  If it is not,
    # then ``f`` is set of a factor of the modulus.  Since the
    # modulus for an ``fq_zech_ctx_t`` is always irreducible, then
    # any non-zero ``op`` will be invertible.

    void fq_zech_trace(fmpz_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Sets ``rop`` to the trace of ``op``.
    # For an element `a \in \mathbf{F}_q`, multiplication by `a` defines
    # a `\mathbf{F}_p`-linear map on `\mathbf{F}_q`.  We define the
    # trace of `a` as the trace of this map.  Equivalently, if `\Sigma`
    # generates `\operatorname{Gal}(\mathbf{F}_q / \mathbf{F}_p)` then the trace of
    # `a` is equal to `\sum_{i=0}^{d-1} \Sigma^i (a)`, where `d =
    # \log_{p} q`.

    void fq_zech_norm(fmpz_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Computes the norm of ``op``.
    # For an element `a \in \mathbf{F}_q`, multiplication by `a` defines
    # a `\mathbf{F}_p`-linear map on `\mathbf{F}_q`.  We define the norm
    # of `a` as the determinant of this map.  Equivalently, if `\Sigma` generates
    # `\operatorname{Gal}(\mathbf{F}_q / \mathbf{F}_p)` then the trace of `a` is equal to
    # `\prod_{i=0}^{d-1} \Sigma^i (a)`, where
    # `d = \text{dim}_{\mathbf{F}_p}(\mathbf{F}_q)`.
    # Algorithm selection is automatic depending on the input.

    void fq_zech_frobenius(fq_zech_t rop, const fq_zech_t op, long e, const fq_zech_ctx_t ctx)
    # Evaluates the homomorphism `\Sigma^e` at ``op``.
    # Recall that `\mathbf{F}_q / \mathbf{F}_p` is Galois with Galois group
    # `\langle \sigma \rangle`, which is also isomorphic to
    # `\mathbf{Z}/d\mathbf{Z}`, where
    # `\sigma \in \operatorname{Gal}(\mathbf{F}_q/\mathbf{F}_p)` is the Frobenius element
    # `\sigma \colon x \mapsto x^p`.

    int fq_zech_multiplicative_order(fmpz * ord, const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Computes the order of ``op`` as an element of the
    # multiplicative group of ``ctx``.
    # Returns 0 if ``op`` is 0, otherwise it returns 1 if ``op``
    # is a generator of the multiplicative group, and -1 if it is not.
    # Note that ``ctx`` must already correspond to a finite field defined by
    # a primitive polynomial and so this function cannot be used to check
    # primitivity of the generator, but can be used to check that other elements
    # are primitive.

    int fq_zech_is_primitive(const fq_zech_t op, const fq_zech_ctx_t ctx)
    # Returns whether ``op`` is primitive, i.e., whether it is a
    # generator of the multiplicative group of ``ctx``.

    void fq_zech_bit_pack(fmpz_t f, const fq_zech_t op, flint_bitcnt_t bit_size, const fq_zech_ctx_t ctx)
    # Packs ``op`` into bitfields of size ``bit_size``, writing the
    # result to ``f``.

    void fq_zech_bit_unpack(fq_zech_t rop, const fmpz_t f, flint_bitcnt_t bit_size, const fq_zech_ctx_t ctx)
    # Unpacks into ``rop`` the element with coefficients packed into
    # fields of size ``bit_size`` as represented by the integer
    # ``f``.
