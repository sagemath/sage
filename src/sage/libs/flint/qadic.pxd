# distutils: libraries = flint
# distutils: depends = flint/qadic.h

from sage.libs.flint.types cimport *

# flint/qadic.h
cdef extern from "flint_wrap.h":
    long qadic_val(const qadic_t op)
    long qadic_prec(const qadic_t op)

    void qadic_ctx_init_conway(qadic_ctx_t ctx, 
                           const fmpz_t p, long d, long min, long max, 
                           const char *var, padic_print_mode mode)
    void qadic_ctx_clear(qadic_ctx_t ctx)
    long qadic_ctx_degree(const qadic_ctx_t ctx)
    void qadic_ctx_print(const qadic_ctx_t ctx)

    #* Memory management *********************************************************/
    void qadic_init(qadic_t x)
    void qadic_init2(qadic_t rop, long prec)
    void qadic_clear(qadic_t x)
    void _fmpz_poly_reduce(fmpz *R, long lenR, const fmpz *a, const long *j, long len)
    void  _fmpz_mod_poly_reduce(fmpz *R, long lenR, 
                      const fmpz *a, const long *j, long len, const fmpz_t p)
    void qadic_reduce(qadic_t x, const qadic_ctx_t ctx)

    #* Randomisation *************************************************************/
    void qadic_randtest(qadic_t x, flint_rand_t state, const qadic_ctx_t ctx)
    void qadic_randtest_not_zero(qadic_t x, flint_rand_t state, const qadic_ctx_t ctx)
    void qadic_randtest_val(qadic_t x, flint_rand_t state, long val, 
                   const qadic_ctx_t ctx)
    void qadic_randtest_int(qadic_t x, flint_rand_t state, const qadic_ctx_t ctx)

    #* Assignments and conversions ***********************************************/
    void qadic_zero(qadic_t op)
    void qadic_one(qadic_t op)
    void qadic_gen(qadic_t x, const qadic_ctx_t ctx)
    void qadic_set_ui(qadic_t rop, unsigned long op, const qadic_ctx_t ctx)
    int qadic_get_padic(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    void qadic_set(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    void qadic_set_fmpz_poly(qadic_t rop, const fmpz_poly_t op, 
                         const qadic_ctx_t ctx)

    #* Comparison ****************************************************************/
    int qadic_is_zero(const qadic_t op)
    int qadic_is_one(const qadic_t op)
    int qadic_equal(const qadic_t op1, const qadic_t op2)

    #* Basic arithmetic **********************************************************/
    void qadic_add(qadic_t x, const qadic_t y, const qadic_t z, const qadic_ctx_t ctx)
    void qadic_sub(qadic_t x, const qadic_t y, const qadic_t z, const qadic_ctx_t ctx)
    void qadic_neg(qadic_t x, const qadic_t y, const qadic_ctx_t ctx)
    void qadic_mul(qadic_t x, const qadic_t y, const qadic_t z,
                          const qadic_ctx_t ctx)
    void _qadic_inv(fmpz *rop, const fmpz *op, long len, 
                const fmpz *a, const long *j, long lena, 
                const fmpz_t p, long N)
    void qadic_inv(qadic_t x, const qadic_t y, const qadic_ctx_t ctx)
    void _qadic_pow(fmpz *rop, const fmpz *op, long len, const fmpz_t e, 
                const fmpz *a, const long *j, long lena, 
                const fmpz_t p)
    void qadic_pow(qadic_t x, const qadic_t y, const fmpz_t e, const qadic_ctx_t ctx)

    #* Special functions *********************************************************/
    void _qadic_exp_rectangular(fmpz *rop, const fmpz *op, long v, long len, 
                            const fmpz *a, const long *j, long lena, 
                            const fmpz_t p, long N, const fmpz_t pN)
    int qadic_exp_rectangular(qadic_t rop, const qadic_t op, 
                          const qadic_ctx_t ctx)
    void _qadic_exp_balanced(fmpz *rop, const fmpz *op, long v, long len, 
                         const fmpz *a, const long *j, long lena, 
                         const fmpz_t p, long N, const fmpz_t pN)
    int qadic_exp_balanced(qadic_t rop, const qadic_t op, 
                       const qadic_ctx_t ctx)
    void _qadic_exp(fmpz *rop, const fmpz *op, long v, long len, 
                           const fmpz *a, const long *j, long lena, 
                           const fmpz_t p, long N, const fmpz_t pN)
    int qadic_exp(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    void _qadic_log_rectangular(fmpz *z, const fmpz *y, long v, long len, 
                            const fmpz *a, const long *j, long lena, 
                            const fmpz_t p, long N, const fmpz_t pN)
    int qadic_log_rectangular(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    void _qadic_log_balanced(fmpz *z, const fmpz *y, long len, 
                         const fmpz *a, const long *j, long lena, 
                         const fmpz_t p, long N, const fmpz_t pN)
    int qadic_log_balanced(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    void _qadic_log(fmpz *z, const fmpz *y, long v, long len, 
                const fmpz *a, const long *j, long lena, 
                const fmpz_t p, long N, const fmpz_t pN)
    int qadic_log(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    void _qadic_frobenius(fmpz *rop, const fmpz *op, long len, long e, 
                  const fmpz *a, const long *j, long lena, 
                  const fmpz_t p, long N)
    void qadic_frobenius(qadic_t rop, const qadic_t op, long e, const qadic_ctx_t ctx)
    void _qadic_teichmuller(fmpz *rop, const fmpz *op, long len, 
                        const fmpz *a, const long *j, long lena, 
                        const fmpz_t p, long N)
    void qadic_teichmuller(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    void _qadic_trace(fmpz_t rop, const fmpz *op, long len, 
                  const fmpz *a, const long *j, long lena, const fmpz_t pN)
    void qadic_trace(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    void _qadic_norm_resultant(fmpz_t rop, const fmpz *op, long len, 
                           const fmpz *a, const long *j, long lena, 
                           const fmpz_t p, long N)
    void _qadic_norm_analytic(fmpz_t rop, const fmpz *y, long v, long len, 
                          const fmpz *a, const long *j, long lena, 
                          const fmpz_t p, long N)
    void _qadic_norm(fmpz_t rop, const fmpz *op, long len, 
                 const fmpz *a, const long *j, long lena, 
                 const fmpz_t p, long N)
    void qadic_norm(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    void qadic_norm_analytic(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    void qadic_norm_resultant(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
    int qadic_sqrt(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)

    #* Output ********************************************************************/
    int qadic_fprint_pretty(FILE *file, const qadic_t op, const qadic_ctx_t ctx)
    int qadic_print_pretty(const qadic_t op, const qadic_ctx_t ctx)
    int qadic_debug(const qadic_t op)
