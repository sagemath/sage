# distutils: libraries = flint
# distutils: depends = flint/padic_poly.h

from sage.libs.flint.types cimport *

# flint/padic_poly.h
cdef extern from "flint_wrap.h":
    #*  Helper functions  ********************************************************/
    long _fmpz_vec_ord_p(const fmpz *vec, long len, const fmpz_t p)

    #*  Memory management  *******************************************************/
    void padic_poly_init(padic_poly_t poly)
    void padic_poly_init2(padic_poly_t poly, long alloc, long prec)
    void padic_poly_clear(padic_poly_t poly)
    void padic_poly_realloc(padic_poly_t f, long alloc, const fmpz_t p)
    void padic_poly_fit_length(padic_poly_t f, long len)
    void _padic_poly_set_length(padic_poly_t poly, long len)
    void _padic_poly_normalise(padic_poly_t f)
    void _padic_poly_canonicalise(fmpz *poly, long *v, long len, const fmpz_t p)
    void padic_poly_canonicalise(padic_poly_t poly, const fmpz_t p)
    void padic_poly_reduce(padic_poly_t f, const padic_ctx_t ctx)
    void padic_poly_truncate(padic_poly_t poly, long n, const fmpz_t p)

    #*  Polynomial parameters  ***************************************************/
    long padic_poly_degree(const padic_poly_t poly)
    long padic_poly_length(const padic_poly_t poly)
    long padic_poly_val(const padic_poly_t poly)
    # macros
    long padic_poly_val(padic_poly_t poly)
    long padic_poly_prec(padic_poly_t poly)

    #*  Randomisation  ***********************************************************/
    void padic_poly_randtest(padic_poly_t f, flint_rand_t state, 
                         long len, const padic_ctx_t ctx)
    void padic_poly_randtest_not_zero(padic_poly_t f, flint_rand_t state,
                                  long len, const padic_ctx_t ctx)
    void padic_poly_randtest_val(padic_poly_t f, flint_rand_t state, 
                             long val, long len, const padic_ctx_t ctx)

    #*  Assignment and basic manipulation  ***************************************/
    void padic_poly_set(padic_poly_t f, const padic_poly_t g, 
                    const padic_ctx_t ctx)
    void padic_poly_set_padic(padic_poly_t poly, const padic_t x, 
                          const padic_ctx_t ctx)
    void padic_poly_set_si(padic_poly_t poly, long x, const padic_ctx_t ctx)
    void padic_poly_set_ui(padic_poly_t poly, unsigned long x, const padic_ctx_t ctx)
    void padic_poly_set_fmpz(padic_poly_t poly, const fmpz_t x, 
                         const padic_ctx_t ctx)
    void padic_poly_set_fmpq(padic_poly_t poly, const fmpq_t x, 
                         const padic_ctx_t ctx)
    void padic_poly_set_fmpz_poly(padic_poly_t rop, const fmpz_poly_t op, 
                              const padic_ctx_t ctx)
    void padic_poly_set_fmpq_poly(padic_poly_t rop, 
                              const fmpq_poly_t op, const padic_ctx_t ctx)
    int padic_poly_get_fmpz_poly(fmpz_poly_t rop, const padic_poly_t op, 
                             const padic_ctx_t ctx)
    void padic_poly_get_fmpq_poly(fmpq_poly_t rop, 
                              const padic_poly_t op, const padic_ctx_t ctx)
    void padic_poly_zero(padic_poly_t poly)
    void padic_poly_one(padic_poly_t poly)
    void padic_poly_swap(padic_poly_t poly1, padic_poly_t poly2)

    #*  Getting and setting coefficients  ****************************************/
    void padic_poly_get_coeff_padic(padic_t c, const padic_poly_t poly, long n, 
                                           const padic_ctx_t ctx)
    void padic_poly_set_coeff_padic(padic_poly_t f, long n, const padic_t c, 
                                                const padic_ctx_t ctx)

    #*  Comparison  **************************************************************/
    int padic_poly_equal(const padic_poly_t f, const padic_poly_t g)
    int padic_poly_is_zero(const padic_poly_t poly)
    int padic_poly_is_one(const padic_poly_t poly)

    #*  Addition and subtraction  ************************************************/
    void _padic_poly_add(fmpz *rop, long *rval, long N, 
                     const fmpz *op1, long val1, long len1, long N1, 
                     const fmpz *op2, long val2, long len2, long N2, 
                     const padic_ctx_t ctx)
    void padic_poly_add(padic_poly_t f, 
                    const padic_poly_t g, const padic_poly_t h, 
                    const padic_ctx_t ctx)
    void _padic_poly_sub(fmpz *rop, long *rval, long N, 
                     const fmpz *op1, long val1, long len1, long N1, 
                     const fmpz *op2, long val2, long len2, long N2, 
                     const padic_ctx_t ctx)
    void padic_poly_sub(padic_poly_t f, 
                    const padic_poly_t g, const padic_poly_t h, 
                    const padic_ctx_t ctx)
    void padic_poly_neg(padic_poly_t f, const padic_poly_t g, 
                    const padic_ctx_t ctx)

#*  Scalar multiplication and division  **************************************/
    void _padic_poly_scalar_mul_padic(fmpz *rop, long *rval, long N, 
                                  const fmpz *op, long val, long len, 
                                  const padic_t c, const padic_ctx_t ctx)
    void padic_poly_scalar_mul_padic(padic_poly_t rop, const padic_poly_t op, 
                                 const padic_t c, const padic_ctx_t ctx)

#*  Multiplication  **********************************************************/
    void _padic_poly_mul(fmpz *rop, long *rval, long N, 
                     const fmpz *op1, long val1, long len1, 
                     const fmpz *op2, long val2, long len2, 
                     const padic_ctx_t ctx)
    void padic_poly_mul(padic_poly_t f, 
                    const padic_poly_t g, const padic_poly_t h, 
                    const padic_ctx_t ctx)

#*  Powering  ****************************************************************/
    void _padic_poly_pow(fmpz *rop, long *rval, long N, 
                     const fmpz *op, long val, long len, unsigned long e,
                     const padic_ctx_t ctx)
    void padic_poly_pow(padic_poly_t rop, const padic_poly_t op, unsigned long e, 
                    const padic_ctx_t ctx)

#*  Series inversion  ********************************************************/
    void padic_poly_inv_series(padic_poly_t Qinv, const padic_poly_t Q, long n, 
                           const padic_ctx_t ctx)

#*  Derivative  **************************************************************/
    void _padic_poly_derivative(fmpz *rop, long *rval, long N, 
                            const fmpz *op, long val, long len, 
                            const padic_ctx_t ctx)
    void padic_poly_derivative(padic_poly_t rop, 
                           const padic_poly_t op, const padic_ctx_t ctx)

#*  Shifting  ****************************************************************/
    void padic_poly_shift_left(padic_poly_t rop, const padic_poly_t op, long n, 
                           const padic_ctx_t ctx)
    void padic_poly_shift_right(padic_poly_t rop, const padic_poly_t op, long n, 
                            const padic_ctx_t ctx)

#*  Evaluation  **************************************************************/
    void _padic_poly_evaluate_padic(fmpz_t u, long *v, long N, 
                                const fmpz *poly, long val, long len, 
                                const fmpz_t a, long b, const padic_ctx_t ctx)
    void padic_poly_evaluate_padic(padic_t y, const padic_poly_t poly, 
                                          const padic_t x, const padic_ctx_t ctx)

#*  Composition  *************************************************************/
    void _padic_poly_compose(fmpz *rop, long *rval, long N, 
                         const fmpz *op1, long val1, long len1, 
                         const fmpz *op2, long val2, long len2, 
                         const padic_ctx_t ctx)
    void padic_poly_compose(padic_poly_t rop, 
                        const padic_poly_t op1, const padic_poly_t op2, 
                        const padic_ctx_t ctx)
    void _padic_poly_compose_pow(fmpz *rop, long *rval, long N, 
                             const fmpz *op, long val, long len, long k, 
                             const padic_ctx_t ctx)
    void padic_poly_compose_pow(padic_poly_t rop, const padic_poly_t op, long k, 
                            const padic_ctx_t ctx)

    #*  Input and output  ********************************************************/
    int padic_poly_debug(const padic_poly_t poly)
    int _padic_poly_fprint(FILE *file, const fmpz *poly, long val, long len, 
                       const padic_ctx_t ctx)
    int padic_poly_fprint(FILE *file, const padic_poly_t poly, 
                      const padic_ctx_t ctx)
    int _padic_poly_print(const fmpz *poly, long val, long len, 
                      const padic_ctx_t ctx)
    int padic_poly_print(const padic_poly_t poly, const padic_ctx_t ctx)
    int _padic_poly_fprint_pretty(FILE *file, 
                              const fmpz *poly, long val, long len, 
                              const char *var, 
                              const padic_ctx_t ctx)
    int padic_poly_fprint_pretty(FILE *file, 
                             const padic_poly_t poly, const char *var, 
                             const padic_ctx_t ctx)
    int _padic_poly_print_pretty(FILE *file, 
                             const fmpz *poly, long val, long len, 
                             const char *var, 
                             const padic_ctx_t ctx)
    int padic_poly_print_pretty(const padic_poly_t poly, const char *var, 
                            const padic_ctx_t ctx)

    #*  Testing  *****************************************************************/
    int _padic_poly_is_canonical(const fmpz *op, long val, long len, 
                             const padic_ctx_t ctx)
    int padic_poly_is_canonical(const padic_poly_t op, const padic_ctx_t ctx)
    int _padic_poly_is_reduced(const fmpz *op, long val, long len, long N, 
                           const padic_ctx_t ctx)
    int padic_poly_is_reduced(const padic_poly_t op, const padic_ctx_t ctx)
