# sage_setup: distribution = sagemath-flint
# distutils: libraries = flint
# distutils: depends = flint/fmpq.h

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport mpq_t
from sage.libs.flint.types cimport fmpz_t, fmpq_t, flint_rand_t, mp_bitcnt_t, fmpz
from sage.libs.mpfr.types cimport mpfr_t, mpfr_rnd_t

# flint/fmpq.h
cdef extern from "flint_wrap.h":
    fmpz * fmpq_numref(fmpq_t)
    fmpz * fmpq_denref(fmpq_t)
    void fmpq_init(fmpq_t)
    void fmpq_clear(fmpq_t)
    void fmpq_one(fmpq_t)
    void fmpq_zero(fmpq_t)
    bint fmpq_is_zero(fmpq_t)
    bint fmpq_is_one(fmpq_t)
    int fmpq_sgn(const fmpq_t x)
    void fmpq_set(fmpq_t dest, const fmpq_t src)
    void fmpq_swap(fmpq_t op1, fmpq_t op2)
    void fmpq_neg(fmpq_t dest, const fmpq_t src)
    void fmpq_abs(fmpq_t dest, const fmpq_t src)
    int fmpq_cmp(const fmpq_t x, const fmpq_t y)
    void fmpq_canonicalise(fmpq_t res)
    int fmpq_is_canonical(const fmpq_t x)
    void fmpq_set_si(fmpq_t res, long p, unsigned long q)
    void fmpq_set_fmpz_frac(fmpq_t res, const fmpz_t p, const fmpz_t q)
    void fmpq_set_mpq(fmpq_t dest, const mpq_t src)
    void fmpq_get_mpq(mpq_t dest, const fmpq_t src)
    int fmpq_get_mpfr(mpfr_t r, const fmpq_t x, mpfr_rnd_t rnd)
    void flint_mpq_init_set_readonly(mpq_t z, const fmpq_t f)
    void flint_mpq_clear_readonly(mpq_t z)
    void fmpq_init_set_readonly(fmpq_t f, const mpq_t z)
    void fmpq_clear_readonly(fmpq_t f)
    char * fmpq_get_str(char * str, int b, const fmpq_t x)
    void fmpq_fprint(FILE * file, const fmpq_t x)
    void fmpq_print(const fmpq_t x)
    void fmpq_randtest(fmpq_t res, flint_rand_t state, mp_bitcnt_t bits)
    void fmpq_randtest_not_zero(fmpq_t res, flint_rand_t state, mp_bitcnt_t bits)
    void fmpq_randbits(fmpq_t res, flint_rand_t state, mp_bitcnt_t bits)
    void fmpq_add(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    void fmpq_add_si(fmpq_t res, const fmpq_t op1, long c)
    void fmpq_add_fmpz(fmpq_t res, const fmpq_t op1, const fmpz_t c)
    void fmpq_sub(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    void fmpq_sub_si(fmpq_t res, const fmpq_t op1, long c)
    void fmpq_sub_fmpz(fmpq_t res, const fmpq_t op1, const fmpz_t c)
    void fmpq_mul(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    void fmpq_mul_fmpz(fmpq_t res, const fmpq_t op, const fmpz_t x)
    void fmpq_pow_si(fmpq_t rop, const fmpq_t op, long e)
    void fmpq_addmul(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    void fmpq_submul(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    void fmpq_inv(fmpq_t dest, const fmpq_t src)
    void fmpq_div(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    void fmpq_div_fmpz(fmpq_t res, const fmpq_t op, const fmpz_t x)
    void fmpq_mul_2exp(fmpq_t res, const fmpq_t x, mp_bitcnt_t exp)
    void fmpq_div_2exp(fmpq_t res, const fmpq_t x, mp_bitcnt_t exp)
    int fmpq_mod_fmpz(fmpz_t res, const fmpq_t x, const fmpz_t mod)
    void fmpq_gcd(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    int fmpq_reconstruct_fmpz(fmpq_t res, const fmpz_t a, const fmpz_t m)
    int fmpq_reconstruct_fmpz_2(fmpq_t res, const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D)
    mp_bitcnt_t fmpq_height_bits(const fmpq_t x)
    void fmpq_height(fmpz_t height, const fmpq_t x)
    void fmpq_next_calkin_wilf(fmpq_t res, const fmpq_t x)
    void fmpq_next_signed_calkin_wilf(fmpq_t res, const fmpq_t x)
    void fmpq_next_minimal(fmpq_t res, const fmpq_t x)
    void fmpq_next_signed_minimal(fmpq_t res, const fmpq_t x)
    long fmpq_get_cfrac(fmpz * c, fmpq_t rem, const fmpq_t x, long n)
    void fmpq_set_cfrac(fmpq_t x, const fmpz * c, long n)
    long fmpq_cfrac_bound(const fmpq_t x)
    void fmpq_dedekind_sum_naive(fmpq_t s, const fmpz_t h, const fmpz_t k)
    void fmpq_dedekind_sum_coprime_large(fmpq_t s, const fmpz_t h, const fmpz_t k)
    double fmpq_dedekind_sum_coprime_d(double h, double k)
    void fmpq_dedekind_sum_coprime(fmpq_t s, const fmpz_t h, const fmpz_t k)
    void fmpq_dedekind_sum(fmpq_t s, const fmpz_t h, const fmpz_t k)
    void fmpq_harmonic_ui(fmpq_t x, unsigned long n)
