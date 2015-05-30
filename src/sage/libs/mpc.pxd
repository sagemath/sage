from sage.libs.mpfr cimport *

cdef extern from "mpc.h":
    ctypedef struct __mpc_struct:
        mpfr_t re
        mpfr_t im
    ctypedef __mpc_struct mpc_t[1]
    ctypedef __mpc_struct* mpc_ptr
    ctypedef __mpc_struct* mpc_srcptr

    ctypedef enum mpc_rnd_t:
        MPC_RNDNN = 0
        MPC_RNDZN = 1
        MPC_RNDUN = 2
        MPC_RNDDN = 3
        MPC_RNDNZ = 16
        MPC_RNDZZ = 17
        MPC_RNDUZ = 18
        MPC_RNDDZ = 19
        MPC_RNDNU = 32
        MPC_RNDZU = 33
        MPC_RNDUU = 34
        MPC_RNDDU = 35
        MPC_RNDND = 48
        MPC_RNDZD = 49
        MPC_RNDUD = 50
        MPC_RNDDD = 51

    # Memory management
    void mpc_init (mpc_t)
    void mpc_init2 (mpc_t, mp_prec_t)
    void mpc_clear (mpc_t)

    # Precision accessors
    mp_prec_t mpc_get_prec (mpc_t)
    void mpc_set_prec (mpc_t, mp_prec_t)

    # Set real part to given value and imaginary part to +0
    int  mpc_set_ui (mpc_t, unsigned long int, mpc_rnd_t)
    int  mpc_set_si (mpc_t, long int, mpc_rnd_t)
    int  mpc_set_z (mpc_t, mpz_t, mpc_rnd_t)
    int  mpc_set_d (mpc_t, double, mpc_rnd_t)
    int  mpc_set_fr (mpc_t, mpfr_t, mpc_rnd_t)
    # Set value
    int  mpc_set (mpc_t, mpc_t, mpc_rnd_t)
    int  mpc_set_ui_ui (mpc_t, unsigned long int, unsigned long int, mpc_rnd_t)
    int  mpc_set_si_si (mpc_t, long int, long int, mpc_rnd_t)
    int  mpc_set_d_d (mpc_t, double, double, mpc_rnd_t)
    int  mpc_set_fr_fr (mpc_t, mpfr_t, mpfr_t, mpc_rnd_t)
    void mpc_set_nan(mpc_t)
    void mpc_swap(mpc_t, mpc_t)

    # Comparisons
    int  mpc_cmp (mpc_t, mpc_t)
    int  mpc_cmp_si_si (mpc_t, long int, long int)
    int  mpc_cmp_si (mpc_t, long int)

    # Projection
    int mpc_real (mpfr_t rop, mpc_t op, mpfr_rnd_t rnd)
    int mpc_imag (mpfr_t rop, mpc_t op, mpfr_rnd_t rnd)
    mpfr_t mpc_realref (mpc_t op)
    mpfr_t mpc_imagref (mpc_t op)
    int mpc_arg (mpfr_t rop, mpc_t op, mpfr_rnd_t rnd)
    int mpc_proj (mpc_t rop, mpc_t op, mpc_rnd_t rnd)

    # Arithmetic
    int  mpc_neg (mpc_t, mpc_t, mpc_rnd_t)
    int  mpc_conj (mpc_t, mpc_t, mpc_rnd_t)

    int  mpc_add (mpc_t, mpc_t, mpc_t, mpc_rnd_t)
    int  mpc_sub (mpc_t, mpc_t, mpc_t, mpc_rnd_t)
    int  mpc_mul (mpc_t, mpc_t, mpc_t, mpc_rnd_t)
    int  mpc_div (mpc_t, mpc_t, mpc_t, mpc_rnd_t)
    int  mpc_sqr (mpc_t, mpc_t, mpc_rnd_t)

    int  mpc_add_fr (mpc_t, mpc_t, mpfr_t, mpc_rnd_t)
    int  mpc_add_ui (mpc_t, mpc_t, unsigned long int, mpc_rnd_t)
    int  mpc_add_si (mpc_t, mpc_t, unsigned long, mpc_rnd_t)
    int  mpc_sub_fr (mpc_t, mpc_t, mpfr_t, mpc_rnd_t)
    int  mpc_fr_sub (mpc_t, mpfr_t, mpc_t, mpc_rnd_t)
    int  mpc_sub_ui (mpc_t, mpc_t, unsigned long int, mpc_rnd_t)
    int  mpc_ui_sub (mpc_t, unsigned long int, mpc_t, mpc_rnd_t)
    int  mpc_ui_ui_sub (mpc_t, unsigned long int, unsigned long int, mpc_t, mpc_rnd_t)
    int  mpc_mul_fr (mpc_t, mpc_t, mpfr_t, mpc_rnd_t)
    int  mpc_mul_ui (mpc_t, mpc_t, unsigned long int, mpc_rnd_t)
    int  mpc_mul_si (mpc_t, mpc_t, long int, mpc_rnd_t)
    int  mpc_mul_i  (mpc_t, mpc_t, int, mpc_rnd_t)
    int  mpc_ui_div (mpc_t, unsigned long int, mpc_t, mpc_rnd_t)
    int  mpc_div_fr (mpc_t, mpc_t, mpfr_t, mpc_rnd_t)
    int  mpc_fr_div (mpc_t, mpfr_t, mpc_t, mpc_rnd_t)
    int  mpc_div_ui (mpc_t, mpc_t, unsigned long int, mpc_rnd_t)
    int  mpc_ui_div (mpc_t, unsigned long int, mpc_t, mpc_rnd_t)
    int  mpc_div_2ui (mpc_t, mpc_t, unsigned long int, mpc_rnd_t)
    int  mpc_div_2si (mpc_t, mpc_t, long int, mpc_rnd_t)
    int  mpc_mul_2ui (mpc_t, mpc_t, unsigned long int, mpc_rnd_t)
    int  mpc_mul_2si (mpc_t, mpc_t, long int, mpc_rnd_t)
    #
    int  mpc_abs (mpfr_t, mpc_t, mp_rnd_t)
    int  mpc_norm (mpfr_t, mpc_t, mp_rnd_t)


    # Power functions and logarithm
    int mpc_sqrt (mpc_t, mpc_t, mpc_rnd_t)
    int mpc_exp (mpc_t, mpc_t, mpc_rnd_t)
    int mpc_log (mpc_t, mpc_t, mpc_rnd_t)
    int mpc_pow (mpc_t rop, mpc_t op1, mpc_t op2, mpc_rnd_t rnd)
    int mpc_pow_si (mpc_t rop, mpc_t op1, long op2, mpc_rnd_t rnd)
    int mpc_pow_ui (mpc_t rop, mpc_t op1, unsigned long op2, mpc_rnd_t rnd)
    int mpc_pow_z (mpc_t rop, mpc_t op1, mpz_t, mpc_rnd_t rnd)
    int mpc_pow_d (mpc_t rop, mpc_t op1, double, mpc_rnd_t rnd)
    int mpc_pow_fr (mpc_t rop, mpc_t op1, mpfr_t op2, mpc_rnd_t rnd)

    # Trigonometric functions
    void mpc_sin (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_cos (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_tan (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_sinh (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_cosh (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_tanh (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_asin (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_acos (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_atan (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_asinh (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_acosh (mpc_t, mpc_t, mpc_rnd_t)
    void mpc_atanh (mpc_t, mpc_t, mpc_rnd_t)

    # Random Function
    int  mpc_urandom (mpc_t, gmp_randstate_t)

    # utility
    int MPC_INEX_RE(int)
    int MPC_INEX_IM(int)
