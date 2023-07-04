# distutils: libraries = flint
# distutils: depends = flint/fmpz_mod.h

# flint/fmpz_mod.h
cdef extern from "flint_wrap.h":
    void fmpz_mod_ctx_init(fmpz_mod_ctx_t ctx, const fmpz_t n)
    void fmpz_mod_ctx_init_ui(fmpz_mod_ctx_t ctx, ulong n)
    void fmpz_mod_ctx_clear(fmpz_mod_ctx_t ctx)

    const fmpz * fmpz_mod_ctx_modulus(const fmpz_mod_ctx_t ctx)

    void fmpz_mod_ctx_set_modulus(fmpz_mod_ctx_t ctx, const fmpz_t n)
    void fmpz_mod_ctx_set_modulus_ui(fmpz_mod_ctx_t ctx, ulong n)

    int fmpz_mod_is_canonical(const fmpz_t a, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_assert_canonical(const fmpz_t a, const fmpz_mod_ctx_t ctx)

    int fmpz_mod_is_one(const fmpz_t a, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_add(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_sub(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_neg(fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_mul(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_inv(fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx)

    int fmpz_mod_divides(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_pow_ui(fmpz_t a, const fmpz_t b, ulong pow, const fmpz_mod_ctx_t ctx)
    int fmpz_mod_pow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t pow, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_discrete_log_pohlig_hellman_init(fmpz_mod_discrete_log_pohlig_hellman_t L)

    void fmpz_mod_discrete_log_pohlig_hellman_clear(fmpz_mod_discrete_log_pohlig_hellman_t L)

    double fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(
                    fmpz_mod_discrete_log_pohlig_hellman_t L,
                    const fmpz_t p)

    void fmpz_mod_discrete_log_pohlig_hellman_run(
                    fmpz_t x,
                    const fmpz_mod_discrete_log_pohlig_hellman_t L,
                    const fmpz_t y)

    const fmpz * fmpz_mod_discrete_log_pohlig_hellman_primitive_root(
                    fmpz_mod_discrete_log_pohlig_hellman_t L)

    int fmpz_next_smooth_prime(fmpz_t a, const fmpz_t b)
