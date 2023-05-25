# distutils: libraries = flint
# distutils: depends = flint/qadic.h

from sage.libs.flint.types cimport *

# flint/qadic.h
cdef extern from "flint_wrap.h":
    # macros
    long padic_val(padic_t x)
    long padic_prec(padic_t x)

    fmpz * padic_unit(const padic_t x)
    slong padic_get_val(const padic_t x)
    slong padic_get_prec(const padic_t x)

    #* Context *******************************************************************/
    void padic_ctx_init(padic_ctx_t ctx, const fmpz_t p, slong min, slong max,
                    padic_print_mode mode)
    void padic_ctx_clear(padic_ctx_t ctx)

    int _padic_ctx_pow_ui(fmpz_t rop, ulong e, const padic_ctx_t ctx)

    void padic_ctx_pow_ui(fmpz_t rop, ulong e, const padic_ctx_t ctx)

    #* Memory management *********************************************************/
    void padic_init(padic_t rop)
    void padic_init2(padic_t rop, slong N)
    void padic_clear(padic_t rop)
    void _padic_canonicalise(padic_t rop, const padic_ctx_t ctx)
    void _padic_reduce(padic_t rop, const padic_ctx_t ctx)
    void padic_reduce(padic_t rop, const padic_ctx_t ctx)

    #* Randomisation *************************************************************/
    void padic_randtest(padic_t rop, flint_rand_t state, const padic_ctx_t ctx)
    void padic_randtest_not_zero(padic_t rop, flint_rand_t state,
                             const padic_ctx_t ctx)
    void padic_randtest_int(padic_t rop, flint_rand_t state,
                        const padic_ctx_t ctx)

    #* Assignments and conversions ***********************************************/
    void padic_set(padic_t rop, const padic_t op, const padic_ctx_t ctx)
    void padic_set_si(padic_t rop, slong op, const padic_ctx_t ctx)
    void padic_set_ui(padic_t rop, ulong op, const padic_ctx_t ctx)
    void padic_set_fmpz(padic_t rop, const fmpz_t op, const padic_ctx_t ctx)
    void padic_set_fmpq(padic_t rop, const fmpq_t op, const padic_ctx_t ctx)
    void padic_set_mpz(padic_t rop, const mpz_t op, const padic_ctx_t ctx)
    void padic_set_mpq(padic_t rop, const mpq_t op, const padic_ctx_t ctx)
    void padic_get_fmpz(fmpz_t rop, const padic_t op, const padic_ctx_t ctx)
    void padic_get_fmpq(fmpq_t rop, const padic_t op, const padic_ctx_t ctx)
    void padic_get_mpz(mpz_t rop, const padic_t op, const padic_ctx_t ctx)
    void padic_get_mpq(mpq_t rop, const padic_t op, const padic_ctx_t ctx)
    void padic_swap(padic_t op1, padic_t op2)
    void padic_zero(padic_t rop)
    void padic_one(padic_t rop)

    #* Comparison ****************************************************************/
    int padic_is_zero(const padic_t op)
    int padic_is_one(const padic_t op)
    int padic_equal(const padic_t op1, const padic_t op2)

    #* Arithmetic operations *****************************************************/
    slong * _padic_lifts_exps(slong *n, slong N)
    void _padic_lifts_pows(fmpz *pow, const slong *a, slong n, const fmpz_t p)
    void padic_add(padic_t rop, const padic_t op1, const padic_t op2,
               const padic_ctx_t ctx)
    void padic_sub(padic_t rop, const padic_t op1, const padic_t op2,
               const padic_ctx_t ctx)
    void padic_neg(padic_t rop, const padic_t op, const padic_ctx_t ctx)
    void padic_mul(padic_t rop, const padic_t op1, const padic_t op2,
               const padic_ctx_t ctx)
    void padic_shift(padic_t rop, const padic_t op, slong v, const padic_ctx_t ctx)
    void padic_div(padic_t rop, const padic_t op1, const padic_t op2,
               const padic_ctx_t ctx)
    void _padic_inv_precompute(padic_inv_t S, const fmpz_t p, slong N)
    void _padic_inv_clear(padic_inv_t S)
    void _padic_inv_precomp(fmpz_t rop, const fmpz_t op, const padic_inv_t S)
    void _padic_inv(fmpz_t rop, const fmpz_t op, const fmpz_t p, slong N)
    void padic_inv(padic_t rop, const padic_t op, const padic_ctx_t ctx)
    int padic_sqrt(padic_t rop, const padic_t op, const padic_ctx_t ctx)
    void padic_pow_si(padic_t rop, const padic_t op, slong e,
                  const padic_ctx_t ctx)

    #* Exponential ***************************************************************/
    slong _padic_exp_bound(slong v, slong N, const fmpz_t p)
    void _padic_exp(fmpz_t rop, const fmpz_t u, slong v, const fmpz_t p, slong N)
    void _padic_exp_rectangular(fmpz_t rop, const fmpz_t u, slong v, const fmpz_t p, slong N)
    void _padic_exp_balanced(fmpz_t rop, const fmpz_t u, slong v, const fmpz_t p, slong N)
    int padic_exp(padic_t rop, const padic_t op, const padic_ctx_t ctx)
    int padic_exp_rectangular(padic_t rop, const padic_t op, const padic_ctx_t ctx)
    int padic_exp_balanced(padic_t rop, const padic_t op, const padic_ctx_t ctx)

    #* Logarithm *****************************************************************/
    slong _padic_log_bound(slong v, slong N, const fmpz_t p)
    void _padic_log(fmpz_t z, const fmpz_t y, slong v, const fmpz_t p, slong N)
    void _padic_log_rectangular(fmpz_t z, const fmpz_t y, slong v, const fmpz_t p, slong N)
    void _padic_log_satoh(fmpz_t z, const fmpz_t y, slong v, const fmpz_t p, slong N)
    void _padic_log_balanced(fmpz_t z, const fmpz_t y, slong v, const fmpz_t p, slong N)
    int padic_log(padic_t rop, const padic_t op, const padic_ctx_t ctx)
    int padic_log_rectangular(padic_t rop, const padic_t op, const padic_ctx_t ctx)
    int padic_log_satoh(padic_t rop, const padic_t op, const padic_ctx_t ctx)
    int padic_log_balanced(padic_t rop, const padic_t op, const padic_ctx_t ctx)

    #* Special functions *********************************************************/
    void _padic_teichmuller(fmpz_t rop, const fmpz_t op, const fmpz_t p, slong N)
    void padic_teichmuller(padic_t rop, const padic_t op, const padic_ctx_t ctx)
    ulong padic_val_fac_ui_2(ulong N)
    ulong padic_val_fac_ui(ulong N, const fmpz_t p)
    void padic_val_fac(fmpz_t rop, const fmpz_t op, const fmpz_t p)

    #* Input and output **********************************************************/
    char * padic_get_str(char * str, const padic_t op, const padic_ctx_t ctx)
    int _padic_fprint(FILE * file, const fmpz_t u, slong v, const padic_ctx_t ctx)
    int padic_fprint(FILE * file, const padic_t op, const padic_ctx_t ctx)

    int _padic_print(const fmpz_t u, slong v, const padic_ctx_t ctx)
    int padic_print(const padic_t op, const padic_ctx_t ctx)
    void padic_debug(const padic_t op)
