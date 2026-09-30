/*
  FLINT utility functions.

#*****************************************************************************
#       Copyright (C) 2019-2026 Kiran S. Kedlaya <kskedl@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

*/

#include <stdlib.h>
#include "flint-utils.h"

#define ROTATE(x, y, z) { x = y; y = z; z = x; }
/*****
  Arithmetic functions

  As with FLINT library functions, aliasing is allowed unless specified.
*****/

inline int is_mpz(fmpz f) {
  return(COEFF_IS_MPZ(f));
}

/* Set res to a*b, interpreting NULL for a as 1. */
void fmpz_opt_mul(fmpz_t res, const fmpz_t a, const fmpz_t b) {
  if (a) fmpz_mul(res, a, b); else fmpz_set(res, b);
}

/* Set res to a*b + c, interpreting NULL for a as 1.
   Aliasing not allowed between res and c unless a == NULL. */
void fmpz_opt_addmul(fmpz_t res, const fmpz_t a, const fmpz_t b, const fmpz_t c) {
  if (a) {
    fmpz_mul(res, a, b);
    fmpz_add(res, res, c);
  }
  else fmpz_add(res, b, c);
}

/* Set res to floor((a+b)/2). */
void fmpz_fmid(fmpz_t res, const fmpz_t a, const fmpz_t b) {
  fmpz_add(res, a, b);
  fmpz_fdiv_q_ui(res, res, 2);
}

/* Set res to ceil((a+b)/2). */
void fmpz_cmid(fmpz_t res, const fmpz_t a, const fmpz_t b) {
  fmpz_add(res, a, b);
  fmpz_cdiv_q_ui(res, res, 2);
}

/* Set res to ceil(sqrt(a)). For the floor, use FLINT's built-in fmpz_sqrt instead. */
void fmpz_sqrt_c(fmpz_t res, const fmpz_t a) {
  int s = fmpz_root(res, a, 2);
  if (!s) fmpz_add_ui(res, res, 1);
}

/* Set res to floor(a/b) if r == 0 and ceil(r/b) if r == 1. Aliasing allowed. */
void fmpz_div_q(fmpz_t res, const fmpz_t a, const fmpz_t b, int r) {
  if (r) fmpz_cdiv_q(res, a, b); else fmpz_fdiv_q(res, a, b);
}

/* Set res to the floor (if r==0) or ceiling (if r==1) of (a/b + c sqrt(q))/d).
   No aliasing allowed. b and d must be positive.
   If b is NULL we interpret it as 1. If c is NULL we interpret it as 0. */
void fmpq_floor_ceil_quad(fmpz_t res, int r, const fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d, const fmpz_t q) {
  const fmpz *tmp;

  if (c) {
    int s = fmpz_sgn(c) > 0;
    fmpz_mul(res, c, c);
    fmpz_mul(res, res, q);
    if (r ^ s) fmpz_sqrt(res, res); else fmpz_sqrt_c(res, res);
    if (!s) fmpz_neg(res, res);
    if (b) fmpz_mul(res, res, b);
    fmpz_add(res, res, a);
    tmp = res;
  } else tmp = a;
  if (b) { fmpz_div_q(res, tmp, b, r); tmp = res; }
  fmpz_div_q(res, tmp, d, r);
}

/* Compute (a*b-c*d)/e assuming that division is exact.
   Aliasing allowed except between res and e. */
void fmpz_fmms_divexact(fmpz *res, const fmpz *a, const fmpz *b, const fmpz *c, const fmpz *d, const fmpz_t e) {
  fmpz_fmms(res, a, b, c, d);
  fmpz_divexact(res, res, e);
}

/* Compute the vector res with res[i] = a[i]*b[i].
   Aliasing allowed except for *partial* overlap between res and a, b. */
void _fmpz_vec_mul(fmpz *res, const fmpz *a, const fmpz *b, int n) {
  for (int i=0; i<n; i++) fmpz_mul(res+i, a+i, b+i);
}

/* Compute the vector res with res[i] = a[i]*b[i] - c[i]*d[i].
   Aliasing allowed except for *partial* overlap between res and a, b, c, or d. */
void _fmpz_vec_fmms(fmpz *res, const fmpz *a, const fmpz *b, const fmpz *c, const fmpz *d, int n) {
  for (int i=0; i<n; i++) fmpz_fmms(res+i, a+i, b+i, c+i, d+i);
}

/* Compute the vector res with res[i] = a[i]*b + c[i]*d.
   No aliasing allowed between res and b or d.
   No partial overlap allowed between res and a or c. */
void _fmpz_vec_scalar_fmma(fmpz *res, const fmpz *a, const fmpz_t b, const fmpz *c, const fmpz_t d, int n) {
  for (int i=0; i<n; i++) fmpz_fmma(res+i, a+i, b, c+i, d);
}

/* Compute the vector res with res[i] = a[i]*b - c[i]*d.
   No aliasing allowed between res and b or d.
   No partial overlap allowed between res and a or c. */
void _fmpz_vec_scalar_fmms(fmpz *res, const fmpz *a, const fmpz_t b, const fmpz *c, const fmpz_t d, int n) {
  for (int i=0; i<n; i++) fmpz_fmms(res+i, a+i, b, c+i, d);
}

/* Compute the vector res with res[i] = (a[i]*b - c[i]*d)/e assuming that division is exact.
   No aliasing allowed between res and b, d, or e.
   No partial overlap allowed between res and a or c. */
void _fmpz_vec_scalar_fmms_divexact(fmpz *res, const fmpz *a, const fmpz *b, const fmpz *c, const fmpz *d, int n, const fmpz *e) {
  _fmpz_vec_scalar_fmms(res, a, b, c, d, n);
  _fmpz_vec_scalar_divexact_fmpz_wrapper(res, res, n, e);
 }

/* Compute the vector res with res[i] = (a[i]*b[i] - c[i]*d[i])/e[i] assuming that division is exact.
   No partial overlap allowed between res and a, b, c, d, or e. */
void _fmpz_vec_fmms_divexact(fmpz *res, const fmpz *a, const fmpz *b, const fmpz *c, const fmpz *d, const fmpz_t e, int n) {
  for (int i=0; i<n; i++) fmpz_fmms_divexact(res+i, a+i, b+i, c+i, d+i, e+i);
}

/* Compute the vector res with res[i] = a[i]*b + c[i].
   Aliasing not allowed between res and b. */
void _fmpz_vec_scalar_fmma_one(fmpz *res, const fmpz *a, const fmpz_t b, const fmpz *c, int n) {
  _fmpz_vec_scalar_mul_fmpz(res, a, n, b);
  _fmpz_vec_add(res, res, c, n);
}

/* Compute the vector res with res[i] = a[i]*b - c[i].
   Aliasing not allowed between res and b. */
void _fmpz_vec_scalar_fmms_one(fmpz *res, const fmpz *a, const fmpz_t b, const fmpz *c, int n) {
  _fmpz_vec_scalar_mul_fmpz(res, a, n, b);
  _fmpz_vec_sub(res, res, c, n);
}

/*
  Compute the Hankel determinant associated to the sequence seq of odd length n
  by a direct application of FLINT's determinant function.
*/

void hankel_determinant_direct(fmpz_t res, const fmpz *seq, int n) {
  fmpz_mat_t mat;
  int s = n/2 + 1;
  fmpz_mat_init(mat, s, s);
  for (int i=0; i<s; i++)
    for (int j=0; j<s; j++)
      fmpz_set(fmpz_mat_entry(mat, i, j), seq+i+j);
  fmpz_mat_det(res, mat);
  fmpz_mat_clear(mat);
  return;
}

/*
  Attempt to compute the Hankel determinant associated to the sequence seq of odd length n
  using Dodgson condensation. Returns 1 for success, 0 for failure.

  This function assumes that {w, 2*n-5} is scratch space.
*/

int hankel_determinant_condensation(fmpz_t res, const fmpz *seq, int n, fmpz *w) {
  if (n == 1) { fmpz_set(res, seq); return(1); }

  int i, n1 = n-2;
  fmpz *f0 = w; // Length n-2
  fmpz *f1 = w+n1; // Length max(0,n-4)
  fmpz *t = (fmpz *)seq;

  _fmpz_vec_fmms(f0, t, t+2, t+1, t+1, n1);
  while (n1 > 1) {
    n1 -= 2;
    for (i=0; i<n1; i++) if (fmpz_is_zero(t+i+2)) return(0); // Failure because of zero division
    _fmpz_vec_fmms_divexact(f1, f0, f0+2, f0+1, f0+1, t+2, n1);
    ROTATE(t, f0, f1)
  }
  fmpz_set(res, f0);
  return(1);
}

/*
    Use a subresultant sequence to test whether a given polynomial has
    real roots. Note that this test has an early abort mechanism:
    having real roots means that the sign sequence has the maximal number
    of sign changes, so the test aborts if any sign change is missed.

    This function assumes that:
        - {poly, n} is a normalized vector with n >= 2 and positive leading coefficient;
        - {f0, n-1} and {f1, n-1} are scratch space.

    We add a (if b is NULL) or a*b (otherwise) to the constant term before testing.

    Based on code by Sebastian Pancratz from the FLINT repository (plus the Ducos variation).
*/

#define sgn_criterion (sgn == 1 || (force_squarefree && sgn == 0))
int _fmpz_poly_all_real_roots(const fmpz *poly, int n, fmpz *f0, fmpz *f1,
                              int force_squarefree, const fmpz_t a, const fmpz_t b) {
  if (n <= 2) return(1);  // Constant or linear polynomial

  /* Put the updated constant term of poly in f1. */
  fmpz_opt_addmul(f1, b, a, poly);

  int sgn;
  if (n == 3) { // Quadratic case, compute discriminant
    fmpz_mul_ui(f1, f1, 4);
    fmpz_fmms_wrapper(f1, f1, poly+2, poly+1, poly+1);
    sgn = fmpz_sgn(f1);
    return(!sgn_criterion);
  }

  /* Introduce some unallocated pointers.
     To gloss these, imagine we are trying to set f0 to be the
     pseudoremainder of f1 modulo f2. */
  fmpz *t; // Receives the leading coefficient of f0
  fmpz *lead1; // Holds the leading coefficient of f1
  fmpz *lead11; // Holds the subleading coefficient of f1
  fmpz *lead2; // Holds the leading coefficient of f2
  fmpz *content; // Holds a value that will be divided (twice) out of f0

  /* At this point deg(poly) = n-1, deg(f1) = n-2.
     Compute the next leading coefficient; instead of the Ducos variation,
     we remove one factor of the leading coefficient of poly. */

  t = f1+n-3; lead1 = f0+n-2; lead11 = lead1-1; lead2 = (fmpz *)poly+n-1;
  fmpz_mul_ui(lead1, lead2, n-1);
  fmpz_mul_ui(t, lead1, 2);
  fmpz_mul_ui(lead11, lead2-1, n-2);
  fmpz_fmms_wrapper(t, t, lead2-2, lead11, lead2-1);

  /* If we miss any one sign change, we cannot have enough. */
  sgn = fmpz_sgn(t);
  if (sgn_criterion) return(0);

  /* Set f0 := deriv(poly). */
  _fmpz_poly_derivative(f0, poly, n-2);

  /* Set f1 to the pseudoremainder of poly modulo f0, again removing one factor of lead2.
     We use t+1 = f1+n-2 as scratch. */

  _fmpz_vec_scalar_fmms(f1+1, poly+1, lead1, f0, lead2, n-4);
  fmpz_mul(f1, f1, lead1);
  fmpz_set_ui(t+1, n-1);
  _fmpz_vec_scalar_fmms(f1, f1, t+1, f0, lead2-1, n-3);

  /* If not forcing squarefree but sgn == 0, we win iff f1 = 0. */
  if (!sgn && !force_squarefree) return (_fmpz_vec_is_zero(f1, n-3));

  /* At this point deg(f0) = n-2, deg(f1) = n-3.
     Start computing the next leading coefficient (up to a positive factor). */

  content = lead1; t = f0+n-4; lead1 = f1+n-3; lead11 = lead1-1;
  fmpz_fmms_wrapper(t, t, lead1, t+1, lead11);
  fmpz_divexact_ui(t, t, n-1);

  _fmpz_vec_scalar_mul_fmpz(f1, f1, n-2, lead2); // Put back in a factor for later

  for (n -= 3; n > 1; n--) {
    /* At this point deg(f0) = n+1, deg(f1) = n.
       Finish computing the next leading coefficient (up to a positive factor). */
    fmpz_sub(t, t, lead1-2);
    if (fmpz_sgn(t) * fmpz_sgn(lead1) > 0) return(0);
    fmpz_fmma_wrapper(t, t, lead1, lead11, lead11);

    /* If we miss any one sign change, we cannot have enough. */
    sgn = fmpz_sgn(t);
    if (sgn_criterion) return(0);

    /* Replace f0 by its pseudoremainder modulo f1, 
       leaving f0[n], f0[n+1] intact as well as f0[n-1] (except to remove content). */
    _fmpz_vec_scalar_fmms_divexact(f0, f0, lead1, f1, t+1, n-1, content);
    _fmpz_vec_sub_wrapper(f0+1, f0+1, f1, n-2);
    _fmpz_vec_scalar_fmma(f0, f0, lead1, f1, lead11, n-1);

    /* If not forcing squarefree but sgn == 0, we win iff f0 = 0. */
    if (!sgn && !force_squarefree) return (_fmpz_vec_is_zero(f0, n-1));

    /* Divide f0 by content. */
    _fmpz_vec_scalar_divexact_fmpz_wrapper(f0, f0, n, content);

    /* Swap f0 with f1 at the pointer level. */
    ROTATE(t, f0, f1)

    /* At this point deg(f0) = n, deg(f1) = n-1.
       Start computing the next leading coefficient (up to a positive factor). */
    content = lead1; t = f0+n-2; lead1 = f1+n-1; lead11 = lead1-1;
    fmpz_fmms_divexact(t, t, lead1, lead11, t+1, content);
  }

  /* Handle the case where the pseudoremainder is a scalar. */
  if (fmpz_sgn(t) * fmpz_sgn(lead1) > 0) return(0);
  fmpz_fmma_wrapper(t, t, lead1, lead11, lead11);
  sgn = fmpz_sgn(t);
  return (!sgn_criterion);
}

