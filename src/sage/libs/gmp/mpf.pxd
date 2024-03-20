# sage_setup: distribution = sagemath-objects
# distutils: libraries = gmp

from sage.libs.gmp.types cimport *

cdef extern from "gmp.h":

    ### Floating-point Functions ###

    # Initialization Functions
    void mpf_set_default_prec (unsigned long int prec)
    unsigned long int mpf_get_default_prec ()
    void mpf_init (mpf_t x)
    void mpf_init2 (mpf_t x, unsigned long int prec)
    void mpf_clear (mpf_t x)
    unsigned long int mpf_get_prec (mpf_t op)
    void mpf_set_prec (mpf_t rop, unsigned long int prec)
    void mpf_set_prec_raw (mpf_t rop, unsigned long int prec)

    # Assignment Functions
    void mpf_set (mpf_t rop, mpf_t op)
    void mpf_set_ui (mpf_t rop, unsigned long int op)
    void mpf_set_si (mpf_t rop, signed long int op)
    void mpf_set_d (mpf_t rop, double op)
    void mpf_set_z (mpf_t rop, mpz_t op)
    void mpf_set_q (mpf_t rop, mpq_t op)
    int mpf_set_str (mpf_t rop, char *str, int base)
    void mpf_swap (mpf_t rop1, mpf_t rop2)

    # Combined Initialization and Assignment Functions
    void mpf_init_set (mpf_t rop, mpf_t op)
    void mpf_init_set_ui (mpf_t rop, unsigned long int op)
    void mpf_init_set_si (mpf_t rop, signed long int op)
    void mpf_init_set_d (mpf_t rop, double op)
    int mpf_init_set_str (mpf_t rop, char *str, int base)

    # Conversion Functions
    double mpf_get_d (mpf_t op)
    double mpf_get_d_2exp (signed long int *exp, mpf_t op)
    long mpf_get_si (mpf_t op)
    unsigned long mpf_get_ui (mpf_t op)
    char * mpf_get_str (char *str, mp_exp_t *expptr, int base, size_t n_digits, mpf_t op)

    # Arithmetic Functions
    void mpf_add (mpf_t rop, mpf_t op1, mpf_t op2)
    void mpf_add_ui (mpf_t rop, mpf_t op1, unsigned long int op2)
    void mpf_sub (mpf_t rop, mpf_t op1, mpf_t op2)
    void mpf_ui_sub (mpf_t rop, unsigned long int op1, mpf_t op2)
    void mpf_sub_ui (mpf_t rop, mpf_t op1, unsigned long int op2)
    void mpf_mul (mpf_t rop, mpf_t op1, mpf_t op2)
    void mpf_mul_ui (mpf_t rop, mpf_t op1, unsigned long int op2)
    void mpf_div (mpf_t rop, mpf_t op1, mpf_t op2)
    void mpf_ui_div (mpf_t rop, unsigned long int op1, mpf_t op2)
    void mpf_div_ui (mpf_t rop, mpf_t op1, unsigned long int op2)
    void mpf_sqrt (mpf_t rop, mpf_t op)
    void mpf_sqrt_ui (mpf_t rop, unsigned long int op)
    void mpf_pow_ui (mpf_t rop, mpf_t op1, unsigned long int op2)
    void mpf_neg (mpf_t rop, mpf_t op)
    void mpf_abs (mpf_t rop, mpf_t op)
    void mpf_mul_2exp (mpf_t rop, mpf_t op1, unsigned long int op2)
    void mpf_div_2exp (mpf_t rop, mpf_t op1, unsigned long int op2)

    # Comparison Functions
    int mpf_cmp (mpf_t op1, mpf_t op2)
    int mpf_cmp_d (mpf_t op1, double op2)
    int mpf_cmp_ui (mpf_t op1, unsigned long int op2)
    int mpf_cmp_si (mpf_t op1, signed long int op2)
    int mpf_eq (mpf_t op1, mpf_t op2, unsigned long int op3)
    void mpf_reldiff (mpf_t rop, mpf_t op1, mpf_t op2)
    int mpf_sgn (mpf_t op)

    # Input and Output Functions
    # size_t mpf_out_str (file *stream, int base, size_t n_digits, mpf_t op)
    # size_t mpf_inp_str (mpf_t rop, file *stream, int base)

    # Miscellaneous Functions
    void mpf_ceil (mpf_t rop, mpf_t op)
    void mpf_floor (mpf_t rop, mpf_t op)
    void mpf_trunc (mpf_t rop, mpf_t op)
    bint mpf_integer_p (mpf_t op)
    bint mpf_fits_ulong_p (mpf_t op)
    bint mpf_fits_slong_p (mpf_t op)
    bint mpf_fits_uint_p (mpf_t op)
    bint mpf_fits_sint_p (mpf_t op)
    bint mpf_fits_ushort_p (mpf_t op)
    bint mpf_fits_sshort_p (mpf_t op)
    void mpf_urandomb (mpf_t rop, gmp_randstate_t state, unsigned long int nbits)
    void mpf_random2 (mpf_t rop, mp_size_t max_size, mp_exp_t exp)
