# distutils: libraries = flint
# distutils: depends = flint/fq_nmod.h

from sage.libs.flint.types cimport fq_nmod_ctx_t, fq_nmod_t, fmpz_t, nmod_poly_t, slong, ulong

# flint/fq_nmod.h
cdef extern from "flint_wrap.h":
    void fq_nmod_ctx_init(fq_nmod_ctx_t ctx, const fmpz_t p, slong d, const char *var)
    void fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx, const fmpz_t p, slong d, const char *var)
    void fq_nmod_ctx_init_modulus(fq_nmod_ctx_t ctx,
                         nmod_poly_t modulus,
                         const char *var)

    void fq_nmod_ctx_clear(fq_nmod_ctx_t ctx)

    fmpz_t fq_nmod_ctx_prime(fq_nmod_ctx_t ctx)
    slong fq_nmod_ctx_degree(const fq_nmod_ctx_t ctx)
    void fq_nmod_ctx_order(fmpz_t f, const fq_nmod_ctx_t ctx)

    void fq_nmod_ctx_print(const fq_nmod_ctx_t ctx)

    #  Memory managment

    void fq_nmod_init(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
    void fq_nmod_clear(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
    void fq_nmod_reduce(fq_nmod_t rop, const fq_nmod_ctx_t ctx)

    #  Basic arithmetic

    void fq_nmod_add(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)
    void fq_nmod_sub(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)
    void fq_nmod_sub_one(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx)
    void fq_nmod_neg(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx)
    void fq_nmod_mul(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)
    void fq_nmod_mul_fmpz(fq_nmod_t rop, const fq_nmod_t op, const fmpz_t x, const fq_nmod_ctx_t ctx)
    void fq_nmod_mul_si(fq_nmod_t rop, const fq_nmod_t op, slong x, const fq_nmod_ctx_t ctx)
    void fq_nmod_mul_ui(fq_nmod_t rop, const fq_nmod_t op, ulong x, const fq_nmod_ctx_t ctx)
    void fq_nmod_sqr(fq_nmod_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
    void fq_nmod_inv(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx)
    void fq_nmod_gcdinv(fq_nmod_t rop, fq_nmod_t inv, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
    void fq_nmod_pow(fq_nmod_t rop, const fq_nmod_t op1, const fmpz_t e, const fq_nmod_ctx_t ctx)
    void fq_nmod_pow_ui(fq_nmod_t rop, const fq_nmod_t op, const ulong e, const fq_nmod_ctx_t ctx)
    void fq_nmod_pth_root(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx)

    #  Comparison

    int fq_nmod_equal(const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)
    int fq_nmod_is_zero(const fq_nmod_t op, const fq_nmod_ctx_t ctx)
    int fq_nmod_is_one(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    #  Assignments and conversions

    void fq_nmod_set(fq_nmod_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
    void fq_nmod_set_fmpz(fq_nmod_t rop, const fmpz_t x, const fq_nmod_ctx_t ctx)
    void fq_nmod_set_ui(fq_nmod_t rop, const ulong x, const fq_nmod_ctx_t ctx)
    void fq_nmod_set_si(fq_nmod_t rop, const slong x, const fq_nmod_ctx_t ctx)
    void fq_nmod_set_coeff_fmpz(fq_nmod_t rop, const fmpz_t x, const ulong n, const fq_nmod_ctx_t ctx)
    void fq_nmod_swap(fq_nmod_t op1, fq_nmod_t op2, const fq_nmod_ctx_t ctx)
    void fq_nmod_zero(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
    void fq_nmod_one(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
    void fq_nmod_gen(fq_nmod_t rop, const fq_nmod_ctx_t ctx)

    #  Output

    void fq_nmod_print(const fq_nmod_t op, const fq_nmod_ctx_t ctx)
    int fq_nmod_print_pretty(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    char * fq_nmod_get_str(const fq_nmod_t op, const fq_nmod_ctx_t ctx)
    char * fq_nmod_get_str_pretty(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

    #  Special functions

    void fq_nmod_trace(fmpz_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
    void fq_nmod_frobenius(fq_nmod_t rop, const fq_nmod_t op, slong e, const fq_nmod_ctx_t ctx)
    void fq_nmod_norm(fmpz_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
