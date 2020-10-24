# distutils: libraries = flint
# distutils: depends = flint/fq.h

from sage.libs.flint.types cimport fq_ctx_t, fq_t, fmpz_t, fmpz_mod_poly_t, slong, ulong

# flint/fq.h
cdef extern from "flint_wrap.h":
    void fq_ctx_init(fq_ctx_t ctx, const fmpz_t p, slong d, const char *var)
    void fq_ctx_init_conway(fq_ctx_t ctx, const fmpz_t p, slong d, const char *var)
    void fq_ctx_init_modulus(fq_ctx_t ctx,
                         fmpz_mod_poly_t modulus,
                         const char *var)

    void fq_ctx_clear(fq_ctx_t ctx)

    fmpz_t fq_ctx_prime(fq_ctx_t ctx)
    slong fq_ctx_degree(const fq_ctx_t ctx)
    void fq_ctx_order(fmpz_t f, const fq_ctx_t ctx)

    void fq_ctx_print(const fq_ctx_t ctx)

    #  Memory managment

    void fq_init(fq_t rop, const fq_ctx_t ctx)
    void fq_clear(fq_t rop, const fq_ctx_t ctx)
    void fq_reduce(fq_t rop, const fq_ctx_t ctx)

    #  Basic arithmetic

    void fq_add(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
    void fq_sub(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
    void fq_sub_one(fq_t rop, const fq_t op1, const fq_ctx_t ctx)
    void fq_neg(fq_t rop, const fq_t op1, const fq_ctx_t ctx)
    void fq_mul(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
    void fq_mul_fmpz(fq_t rop, const fq_t op, const fmpz_t x, const fq_ctx_t ctx)
    void fq_mul_si(fq_t rop, const fq_t op, slong x, const fq_ctx_t ctx)
    void fq_mul_ui(fq_t rop, const fq_t op, ulong x, const fq_ctx_t ctx)
    void fq_div(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
    void fq_sqr(fq_t rop, const fq_t op, const fq_ctx_t ctx)
    void fq_inv(fq_t rop, const fq_t op1, const fq_ctx_t ctx)
    void fq_gcdinv(fq_t rop, fq_t inv, const fq_t op, const fq_ctx_t ctx)
    void fq_pow(fq_t rop, const fq_t op1, const fmpz_t e, const fq_ctx_t ctx)
    void fq_pow_ui(fq_t rop, const fq_t op, const ulong e, const fq_ctx_t ctx)
    void fq_pth_root(fq_t rop, const fq_t op1, const fq_ctx_t ctx)

    #  Comparison

    int fq_equal(const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
    int fq_is_zero(const fq_t op, const fq_ctx_t ctx)
    int fq_is_one(const fq_t op, const fq_ctx_t ctx)

    #  Assignments and conversions

    void fq_set(fq_t rop, const fq_t op, const fq_ctx_t ctx)
    void fq_set_fmpz(fq_t rop, const fmpz_t x, const fq_ctx_t ctx)
    void fq_set_ui(fq_t rop, const ulong x, const fq_ctx_t ctx)
    void fq_set_si(fq_t rop, const slong x, const fq_ctx_t ctx)
    void fq_set_coeff_fmpz(fq_t rop, const fmpz_t x, const ulong n, const fq_ctx_t ctx)
    void fq_swap(fq_t op1, fq_t op2, const fq_ctx_t ctx)
    void fq_zero(fq_t rop, const fq_ctx_t ctx)
    void fq_one(fq_t rop, const fq_ctx_t ctx)
    void fq_gen(fq_t rop, const fq_ctx_t ctx)

    #  Output

    void fq_print(const fq_t op, const fq_ctx_t ctx)
    int fq_print_pretty(const fq_t op, const fq_ctx_t ctx)

    char * fq_get_str(const fq_t op, const fq_ctx_t ctx)
    char * fq_get_str_pretty(const fq_t op, const fq_ctx_t ctx)

    #  Special functions

    void fq_trace(fmpz_t rop, const fq_t op, const fq_ctx_t ctx)
    void fq_frobenius(fq_t rop, const fq_t op, slong e, const fq_ctx_t ctx)
    void fq_norm(fmpz_t rop, const fq_t op, const fq_ctx_t ctx)
