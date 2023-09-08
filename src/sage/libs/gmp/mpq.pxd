# sage_setup: distribution = sagemath-objects
# distutils: libraries = gmp

from sage.libs.gmp.types cimport *

cdef extern from "gmp.h":

    ### Rational Functions ###

    void mpq_canonicalize (mpq_t op)

    # Initialization and Assignment Functions
    void mpq_init (mpq_t dest_rational)
    void mpq_clear (mpq_t rational_number)
    void mpq_set (mpq_t rop, mpq_t op)
    void mpq_set_z (mpq_t rop, mpz_t op)
    void mpq_set_ui (mpq_t rop, unsigned long int op1, unsigned long int op2)
    void mpq_set_si (mpq_t rop, signed long int op1, unsigned long int op2)
    int mpq_set_str (mpq_t rop, char *str, int base)
    void mpq_swap (mpq_t rop1, mpq_t rop2)

    # Conversion Functions
    double mpq_get_d (mpq_t op)
    void mpq_set_d (mpq_t rop, double op)
    void mpq_set_f (mpq_t rop, mpf_t op)
    char * mpq_get_str (char *str, int base, mpq_t op)

    # Arithmetic Functions
    void mpq_add (mpq_t sum, mpq_t addend1, mpq_t addend2)
    void mpq_sub (mpq_t difference, mpq_t minuend, mpq_t subtrahend)
    void mpq_mul (mpq_t product, mpq_t multiplier, mpq_t multiplicand)
    void mpq_mul_2exp (mpq_t rop, mpq_t op1, unsigned long int op2)
    void mpq_div (mpq_t quotient, mpq_t dividend, mpq_t divisor)
    void mpq_div_2exp (mpq_t rop, mpq_t op1, unsigned long int op2)
    void mpq_neg (mpq_t negated_operand, mpq_t operand)
    void mpq_abs (mpq_t rop, mpq_t op)
    void mpq_inv (mpq_t inverted_number, mpq_t number)

    # Comparison Functions
    int mpq_cmp (mpq_t op1, mpq_t op2)
    int mpq_cmp_ui (mpq_t op1, unsigned long int num2, unsigned long int den2)
    int mpq_cmp_si (mpq_t op1, long int num2, unsigned long int den2)
    int mpq_cmp_z (const mpq_t op1, const mpz_t op2)
    int mpq_sgn (mpq_t op)
    int mpq_equal (mpq_t op1, mpq_t op2)

    # Applying Integer Functions to Rationals
    mpz_t mpq_numref (mpq_t op)
    mpz_t mpq_denref (mpq_t op)
    void mpq_get_num (mpz_t numerator, mpq_t rational)
    void mpq_get_den (mpz_t denominator, mpq_t rational)
    void mpq_set_num (mpq_t rational, mpz_t numerator)
    void mpq_set_den (mpq_t rational, mpz_t denominator)

    # Input and Output Functions
    # size_t mpq_out_str (file *stream, int base, mpq_t op)
    # size_t mpq_inp_str (mpq_t rop, file *stream, int base)

