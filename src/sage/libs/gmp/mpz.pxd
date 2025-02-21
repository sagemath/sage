# sage_setup: distribution = sagemath-objects
# distutils: libraries = gmp

from sage.libs.gmp.types cimport *
from libc.stdio cimport FILE

from libc.stdint cimport intmax_t, uintmax_t

cdef extern from "gmp.h":

    ### Integer Functions ###

    # Initialization Functions
    void mpz_init (mpz_t integer)
    void mpz_init2 (mpz_t integer, unsigned long n)
    void mpz_clear (mpz_t integer)
    void mpz_realloc2 (mpz_t integer, unsigned long n)

    # Assignment Functions
    void mpz_set (mpz_t rop, mpz_t op)
    void mpz_set_ui (mpz_t rop, unsigned long int op)
    void mpz_set_si (mpz_t rop, signed long int op)
    void mpz_set_ux (mpz_t rop, uintmax_t op)
    void mpz_set_sx (mpz_t rop, intmax_t op)
    void mpz_set_d (mpz_t rop, double op)
    void mpz_set_q (mpz_t rop, mpq_t op)
    void mpz_set_f (mpz_t rop, mpf_t op)
    int mpz_set_str (mpz_t rop, char *str, int base)
    void mpz_swap (mpz_t rop1, mpz_t rop2)

    # Combined Initialization and Assignment Functions
    void mpz_init_set (mpz_t rop, mpz_t op)
    void mpz_init_set_ui (mpz_t rop, unsigned long int op)
    void mpz_init_set_si (mpz_t rop, signed long int op)
    void mpz_init_set_ux (mpz_t rop, uintmax_t op)
    void mpz_init_set_sx (mpz_t rop, intmax_t op)
    void mpz_init_set_d (mpz_t rop, double op)
    int mpz_init_set_str (mpz_t rop, char *str, int base)

    # Conversion Functions
    unsigned long int mpz_get_ui (mpz_t op)
    signed long int mpz_get_si (mpz_t op)
    uintmax_t mpz_get_ux (mpz_t op)
    intmax_t mpz_get_sx (mpz_t op)
    double mpz_get_d (mpz_t op)
    double mpz_get_d_2exp (long int *exp, mpz_t op)
    char * mpz_get_str (char *str, int base, mpz_t op)

    # Arithmetic Functions
    void mpz_add (mpz_t rop, mpz_t op1, mpz_t op2)
    void mpz_add_ui (mpz_t rop, mpz_t op1, unsigned long int op2)
    void mpz_sub (mpz_t rop, mpz_t op1, mpz_t op2)
    void mpz_sub_ui (mpz_t rop, mpz_t op1, unsigned long int op2)
    void mpz_ui_sub (mpz_t rop, unsigned long int op1, mpz_t op2)
    void mpz_mul (mpz_t rop, mpz_t op1, mpz_t op2)
    void mpz_mul_si (mpz_t rop, mpz_t op1, long int op2)
    void mpz_mul_ui (mpz_t rop, mpz_t op1, unsigned long int op2)
    void mpz_addmul (mpz_t rop, mpz_t op1, mpz_t op2)
    void mpz_addmul_ui (mpz_t rop, mpz_t op1, unsigned long int op2)
    void mpz_submul (mpz_t rop, mpz_t op1, mpz_t op2)
    void mpz_submul_ui (mpz_t rop, mpz_t op1, unsigned long int op2)
    void mpz_mul_2exp (mpz_t rop, mpz_t op1, unsigned long int op2)
    void mpz_neg (mpz_t rop, mpz_t op)
    void mpz_abs (mpz_t rop, mpz_t op)

    # Division Functions
    void mpz_cdiv_q (mpz_t q, mpz_t n, mpz_t d)
    void mpz_cdiv_r (mpz_t r, mpz_t n, mpz_t d)
    void mpz_cdiv_qr (mpz_t q, mpz_t r, mpz_t n, mpz_t d)
    unsigned long int mpz_cdiv_q_ui (mpz_t q, mpz_t n, unsigned long int d)
    unsigned long int mpz_cdiv_r_ui (mpz_t r, mpz_t n, unsigned long int d)
    unsigned long int mpz_cdiv_qr_ui (mpz_t q, mpz_t r, mpz_t n, unsigned long int d)
    unsigned long int mpz_cdiv_ui (mpz_t n, unsigned long int d)
    void mpz_cdiv_q_2exp (mpz_t q, mpz_t n, unsigned long int b)
    void mpz_cdiv_r_2exp (mpz_t r, mpz_t n, unsigned long int b)
    void mpz_fdiv_q (mpz_t q, mpz_t n, mpz_t d)
    void mpz_fdiv_r (mpz_t r, mpz_t n, mpz_t d)
    void mpz_fdiv_qr (mpz_t q, mpz_t r, mpz_t n, mpz_t d)
    unsigned long int mpz_fdiv_q_ui (mpz_t q, mpz_t n, unsigned long int d)
    unsigned long int mpz_fdiv_r_ui (mpz_t r, mpz_t n, unsigned long int d)
    unsigned long int mpz_fdiv_qr_ui (mpz_t q, mpz_t r, mpz_t n, unsigned long int d)
    unsigned long int mpz_fdiv_ui (mpz_t n, unsigned long int d)
    void mpz_fdiv_q_2exp (mpz_t q, mpz_t n, unsigned long int b)
    void mpz_fdiv_r_2exp (mpz_t r, mpz_t n, unsigned long int b)
    void mpz_tdiv_q (mpz_t q, mpz_t n, mpz_t d)
    void mpz_tdiv_r (mpz_t r, mpz_t n, mpz_t d)
    void mpz_tdiv_qr (mpz_t q, mpz_t r, mpz_t n, mpz_t d)
    unsigned long int mpz_tdiv_q_ui (mpz_t q, mpz_t n, unsigned long int d)
    unsigned long int mpz_tdiv_r_ui (mpz_t r, mpz_t n, unsigned long int d)
    unsigned long int mpz_tdiv_qr_ui (mpz_t q, mpz_t r, mpz_t n, unsigned long int d)
    unsigned long int mpz_tdiv_ui (mpz_t n, unsigned long int d)
    void mpz_tdiv_q_2exp (mpz_t q, mpz_t n, unsigned long int b)
    void mpz_tdiv_r_2exp (mpz_t r, mpz_t n, unsigned long int b)
    void mpz_mod (mpz_t r, mpz_t n, mpz_t d)
    unsigned long int mpz_mod_ui (mpz_t r, mpz_t n, unsigned long int d)
    void mpz_divexact (mpz_t q, mpz_t n, mpz_t d)
    void mpz_divexact_ui (mpz_t q, mpz_t n, unsigned long d)
    bint mpz_divisible_p (mpz_t n, mpz_t d)
    bint mpz_divisible_ui_p (mpz_t n, unsigned long int d)
    bint mpz_divisible_2exp_p (mpz_t n, unsigned long int b)
    bint mpz_congruent_p (mpz_t n, mpz_t c, mpz_t d)
    bint mpz_congruent_ui_p (mpz_t n, unsigned long int c, unsigned long int d)
    bint mpz_congruent_2exp_p (mpz_t n, mpz_t c, unsigned long int b)

    # Exponentiation Functions
    void mpz_powm (mpz_t rop, mpz_t base, mpz_t exp, mpz_t mod)
    void mpz_powm_ui (mpz_t rop, mpz_t base, unsigned long int exp, mpz_t mod)
    void mpz_pow_ui (mpz_t rop, mpz_t base, unsigned long int exp)
    void mpz_ui_pow_ui (mpz_t rop, unsigned long int base, unsigned long int exp)

    # Root Extraction Functions
    int mpz_root (mpz_t rop, mpz_t op, unsigned long int n)
    void mpz_rootrem (mpz_t root, mpz_t rem, mpz_t u, unsigned long int n)
    void mpz_sqrt (mpz_t rop, mpz_t op)
    void mpz_sqrtrem (mpz_t rop1, mpz_t rop2, mpz_t op)
    bint mpz_perfect_power_p (mpz_t op)
    bint mpz_perfect_square_p (mpz_t op)

    # Number Theoretic Functions
    bint mpz_probab_prime_p (mpz_t n, int reps)
    void mpz_nextprime (mpz_t rop, mpz_t op)
    void mpz_gcd (mpz_t rop, mpz_t op1, mpz_t op2)
    unsigned long int mpz_gcd_ui (mpz_t rop, mpz_t op1, unsigned long int op2)
    void mpz_gcdext (mpz_t g, mpz_t s, mpz_t t, mpz_t a, mpz_t b)
    void mpz_lcm (mpz_t rop, mpz_t op1, mpz_t op2)
    void mpz_lcm_ui (mpz_t rop, mpz_t op1, unsigned long op2)
    int mpz_invert (mpz_t rop, mpz_t op1, mpz_t op2)
    int mpz_jacobi (mpz_t a, mpz_t b)
    int mpz_legendre (mpz_t a, mpz_t p)
    int mpz_kronecker (mpz_t a, mpz_t b)
    int mpz_kronecker_si (mpz_t a, long b)
    int mpz_kronecker_ui (mpz_t a, unsigned long b)
    int mpz_si_kronecker (long a, mpz_t b)
    int mpz_ui_kronecker (unsigned long a, mpz_t b)
    unsigned long int mpz_remove (mpz_t rop, mpz_t op, mpz_t f)
    void mpz_fac_ui (mpz_t rop, unsigned long int op)
    void mpz_2fac_ui (mpz_t rop, unsigned long int op)
    void mpz_mfac_uiui (mpz_t rop, unsigned long int n, unsigned long int m)
    void mpz_bin_ui (mpz_t rop, mpz_t n, unsigned long int k)
    void mpz_bin_uiui (mpz_t rop, unsigned long int n, unsigned long int k)
    void mpz_fib_ui (mpz_t fn, unsigned long int n)
    void mpz_fib2_ui (mpz_t fn, mpz_t fnsub1, unsigned long int n)
    void mpz_lucnum_ui (mpz_t ln, unsigned long int n)
    void mpz_lucnum2_ui (mpz_t ln, mpz_t lnsub1, unsigned long int n)

    # Comparison Functions
    int mpz_cmp (mpz_t op1, mpz_t op2)
    int mpz_cmp_d (mpz_t op1, double op2)
    int mpz_cmp_si (mpz_t op1, signed long int op2)
    int mpz_cmp_ui (mpz_t op1, unsigned long int op2)
    int mpz_cmpabs (mpz_t op1, mpz_t op2)
    int mpz_cmpabs_d (mpz_t op1, double op2)
    int mpz_cmpabs_ui (mpz_t op1, unsigned long int op2)
    int mpz_sgn (mpz_t op)

    # Logical and Bit Manipulation Functions
    void mpz_and (mpz_t rop, mpz_t op1, mpz_t op2)
    void mpz_ior (mpz_t rop, mpz_t op1, mpz_t op2)
    void mpz_xor (mpz_t rop, mpz_t op1, mpz_t op2)
    void mpz_com (mpz_t rop, mpz_t op)
    unsigned long int mpz_popcount (mpz_t op)
    unsigned long int mpz_hamdist (mpz_t op1, mpz_t op2)
    unsigned long int mpz_scan0 (mpz_t op, unsigned long int starting_bit)
    unsigned long int mpz_scan1 (mpz_t op, unsigned long int starting_bit)
    void mpz_setbit (mpz_t rop, unsigned long int bit_index)
    void mpz_clrbit (mpz_t rop, unsigned long int bit_index)
    void mpz_combit (mpz_t rop, unsigned long int bit_index)
    int mpz_tstbit (mpz_t op, unsigned long int bit_index)

    # Input and Output Functions
    size_t mpz_out_str (FILE *stream, int base, mpz_t op)
    size_t mpz_inp_str (mpz_t rop, FILE *stream, int base)
    size_t mpz_out_raw (FILE *stream, mpz_t op)
    size_t mpz_inp_raw (mpz_t rop, FILE *stream)

    # Random Number Functions
    void mpz_urandomb (mpz_t rop, gmp_randstate_t state, unsigned long int n)
    void mpz_urandomm (mpz_t rop, gmp_randstate_t state, mpz_t n)
    void mpz_rrandomb (mpz_t rop, gmp_randstate_t state, unsigned long int n)
    void mpz_random (mpz_t rop, mp_size_t max_size)
    void mpz_random2 (mpz_t rop, mp_size_t max_size)

    # Integer Import and Export
    void mpz_import (mpz_t rop, size_t count, int order, int size, int endian, size_t nails, void *op)
    void * mpz_export (void *rop, size_t *countp, int order, int size, int endian, size_t nails, mpz_t op)

    # Miscellaneous Functions
    bint mpz_fits_ulong_p (mpz_t op)
    bint mpz_fits_slong_p (mpz_t op)
    bint mpz_fits_uint_p (mpz_t op)
    bint mpz_fits_sint_p (mpz_t op)
    bint mpz_fits_ushort_p (mpz_t op)
    bint mpz_fits_sshort_p (mpz_t op)
    bint mpz_odd_p (mpz_t op)
    bint mpz_even_p (mpz_t op)
    size_t mpz_sizeinbase (mpz_t op, int base)

    # Special Functions
    void * _mpz_realloc (mpz_t integer, mp_size_t new_alloc)
    mp_limb_t mpz_getlimbn (mpz_t op, mp_size_t n)
    size_t mpz_size (mpz_t op)
    mpz_srcptr mpz_roinit_n (mpz_t x, const mp_limb_t *xp, mp_size_t xs)
