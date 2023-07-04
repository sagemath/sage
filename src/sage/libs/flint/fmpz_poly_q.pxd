# distutils: libraries = flint
# distutils: depends = flint/qadic.h

from sage.libs.flint.types cimport *

# flint/qadic.h
cdef extern from "flint_wrap.h":
    #* Accessing numerator and denominator ***************************************/
    # macros
    fmpz_poly_struct* fmpz_poly_q_numref(fmpz_poly_q_t op)
    fmpz_poly_struct* fmpz_poly_q_denref(fmpz_poly_q_t op)

    void fmpz_poly_q_canonicalise(fmpz_poly_q_t rop)
    int fmpz_poly_q_is_canonical(const fmpz_poly_q_t op)

    #* Memory management *********************************************************/
    void fmpz_poly_q_init(fmpz_poly_q_t rop)
    void fmpz_poly_q_clear(fmpz_poly_q_t rop)

    #* Randomisation *************************************************************/
    void fmpz_poly_q_randtest(fmpz_poly_q_t poly, flint_rand_t state,
                          long len1, mp_bitcnt_t bits1, 
                          long len2, mp_bitcnt_t bits2)
    void fmpz_poly_q_randtest_not_zero(fmpz_poly_q_t poly, flint_rand_t state, 
                                   long len1, mp_bitcnt_t bits1, 
                                   long len2, mp_bitcnt_t bits2)

    #* Assignment ****************************************************************/
    void fmpz_poly_q_set(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
    void fmpz_poly_q_set_si(fmpz_poly_q_t rop, long op)
    void fmpz_poly_q_swap(fmpz_poly_q_t op1, fmpz_poly_q_t op2)
    void fmpz_poly_q_zero(fmpz_poly_q_t rop)
    void fmpz_poly_q_one(fmpz_poly_q_t rop)
    void fmpz_poly_q_neg(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
    void fmpz_poly_q_inv(fmpz_poly_q_t rop, const fmpz_poly_q_t op)

    #* Comparison ****************************************************************/
    int fmpz_poly_q_is_zero(const fmpz_poly_q_t op)
    int fmpz_poly_q_is_one(const fmpz_poly_q_t op)
    int fmpz_poly_q_equal(const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)

    #* Addition and subtraction **************************************************/
    void fmpz_poly_q_add_in_place(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
    void fmpz_poly_q_sub_in_place(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
    void fmpz_poly_q_add(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
    void fmpz_poly_q_sub(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
    void fmpz_poly_q_addmul(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
    void fmpz_poly_q_submul(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)

    #* Scalar multiplication and division ****************************************/
    void fmpz_poly_q_scalar_mul_si(fmpz_poly_q_t rop, const fmpz_poly_q_t op, long x)
    void fmpz_poly_q_scalar_mul_mpz(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const mpz_t x)
    void fmpz_poly_q_scalar_mul_mpq(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const mpq_t x)
    void fmpz_poly_q_scalar_div_si(fmpz_poly_q_t rop, const fmpz_poly_q_t op, long x)
    void fmpz_poly_q_scalar_div_mpz(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const mpz_t x)
    void fmpz_poly_q_scalar_div_mpq(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const mpq_t x)

    #* Multiplication and division ***********************************************/
    void fmpz_poly_q_mul(fmpz_poly_q_t rop, 
                     const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
    void fmpz_poly_q_div(fmpz_poly_q_t rop, 
                     const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)

    #* Powering ******************************************************************/
    void fmpz_poly_q_pow(fmpz_poly_q_t rop, const fmpz_poly_q_t op, unsigned long exp)

    #* Derivative ****************************************************************/
    void fmpz_poly_q_derivative(fmpz_poly_q_t rop, const fmpz_poly_q_t op)

    #* Evaluation ****************************************************************/
    int fmpz_poly_q_evaluate(mpq_t rop, const fmpz_poly_q_t f, const mpq_t a)

    #* Input and output **********************************************************/
    int fmpz_poly_q_set_str(fmpz_poly_q_t rop, const char *s)
    char * fmpz_poly_q_get_str(const fmpz_poly_q_t op)
    char * fmpz_poly_q_get_str_pretty(const fmpz_poly_q_t op, const char *x)
    int fmpz_poly_q_print(const fmpz_poly_q_t op)
    int fmpz_poly_q_print_pretty(const fmpz_poly_q_t op, const char *x)
