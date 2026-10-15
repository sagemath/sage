#pragma once
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>

#define DEBUG 0

#if DEBUG

/* Wrapper functions for FLINT commands that don't show up in profiling data. */
void fmpz_fmma_wrapper(fmpz_t f, const fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d) {
  fmpz_fmma(f, a, b, c, d);
}

void fmpz_fmms_wrapper(fmpz_t f, const fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d) {
  fmpz_fmms(f, a, b, c, d);
}

void _fmpz_vec_sub_wrapper(fmpz *res, const fmpz *vec1, const fmpz *vec2, slong len2) {
  _fmpz_vec_sub(res, vec1, vec2, len2);
}

void _fmpz_vec_scalar_divexact_fmpz_wrapper(fmpz *vec1, const fmpz *vec2, slong len2, const fmpz_t x) {
  _fmpz_vec_scalar_divexact_fmpz(vec1, vec2, len2, x);
}

#else
  #define fmpz_fmma_wrapper(f, a, b, c, d) fmpz_fmma(f, a, b, c, d)
  #define fmpz_fmms_wrapper(f, a, b, c, d) fmpz_fmms(f, a, b, c, d)
  #define _fmpz_vec_sub_wrapper(res, vec1, vec2, len2) _fmpz_vec_sub(res, vec1, vec2, len2)
  #define _fmpz_vec_scalar_divexact_fmpz_wrapper(vec1, vec2, len2, x) _fmpz_vec_scalar_divexact_fmpz(vec1, vec2, len2, x)
#endif

inline int is_mpz(fmpz f);

void fmpz_opt_mul(fmpz_t res, const fmpz_t a, const fmpz_t b);
void fmpz_opt_addmul(fmpz_t res, const fmpz_t a, const fmpz_t b, const fmpz_t c);
void fmpz_fmid(fmpz_t res, const fmpz_t a, const fmpz_t b);
void fmpz_cmid(fmpz_t res, const fmpz_t a, const fmpz_t b);
void fmpz_sqrt_c(fmpz_t res, const fmpz_t a);
void fmpz_div_q(fmpz_t res, const fmpz_t a, const fmpz_t b, int r);
void fmpq_floor_ceil_quad(fmpz_t res, int r, const fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d, const fmpz_t q);

void fmpz_fmms_divexact(fmpz *res, const fmpz *a, const fmpz *b, const fmpz *c, const fmpz *d, const fmpz_t e);
void _fmpz_vec_fmms(fmpz *res, const fmpz *a, const fmpz *b, const fmpz *c, const fmpz *d, int n);
void _fmpz_vec_scalar_fmma(fmpz *res, const fmpz *a, const fmpz_t b, const fmpz *c, const fmpz_t d, int n);
void _fmpz_vec_scalar_fmms(fmpz *res, const fmpz *a, const fmpz_t b, const fmpz *c, const fmpz_t d, int n);
void _fmpz_vec_scalar_fmms_divexact(fmpz *res, const fmpz *a, const fmpz *b, const fmpz *c, const fmpz *d, int n, const fmpz *e);
void _fmpz_vec_fmms_divexact(fmpz *res, const fmpz *a, const fmpz *b, const fmpz *c, const fmpz *d, const fmpz_t e, int n);
void _fmpz_vec_scalar_fmma_one(fmpz *res, const fmpz *a, const fmpz_t b, const fmpz *c, int n);
void _fmpz_vec_scalar_fmms_one(fmpz *res, const fmpz *a, const fmpz_t b, const fmpz *c, int n);

void hankel_determinant_direct(fmpz_t res, const fmpz *seq, int n);
int hankel_determinant_condensation(fmpz_t res, const fmpz *seq, int n, fmpz *w);

int _fmpz_poly_all_real_roots(const fmpz *poly, int n, fmpz *f0, fmpz *f1,
                              int force_squarefree, const fmpz_t a, const fmpz_t b);


