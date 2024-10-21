/* sage_setup: distribution = sagemath-flint
 */
/*
  Low-level code to exhaust over trees of Weil polynomials.
  This code does not implement parallelism; see the Cython wrapper.

  TODO: check for memory leaks.
  TODO: try the Routh-Hurwitz criterion.

#*****************************************************************************
#       Copyright (C) 2019 Kiran S. Kedlaya <kskedl@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

*/

#include "power_sums.h"

/* Check for OpenMP at runtime.
*/
int has_openmp() {
  #if defined(_OPENMP)
  return(1);
  #endif
  return(0);
}

/*
    Use a subresultant (Sturm-Habicht) sequence to test whether a given
    polynomial has all real roots. Note that this test has an early abort
    mechanism: having all real roots means that the sign sequence has
    the maximal number of sign changes, so the test aborts as soon
    as a sign change is missed.

    This function assumes that:
        - {poly, n} is a normalized vector with n >= 2
        - {w, 2*n+1} is scratch space.
    If a and b are not NULL, we add a*b to the constant term before testing.

    Based on code by Sebastian Pancratz from the FLINT repository.
    TODO: compare with floating-point interval arithmetic.
*/

int _fmpz_poly_all_real_roots(fmpz *poly, long n, fmpz *w, int force_squarefree,
			      const fmpz_t a, const fmpz_t b) {
  fmpz *f0     = w + 0*n;
  fmpz *f1     = w + 1*n;
  fmpz *c      = w + 2*n;
  fmpz *d      = w + 2*n+1;
  fmpz *t;

  _fmpz_vec_set(f0, poly, n);
  /* Sanitize input so that n = deg(f0). */
  while ((n > 2) && fmpz_is_zero(f0+n-1))
    n--;
  if (n <= 2) return(1);
  if (a != NULL && b != NULL) fmpz_addmul(f0, a, b);
  _fmpz_poly_derivative(f1, f0, n);
  n--;
  int sgn0_l = fmpz_sgn(f0+n);

  while (1) {
    /* At this point deg(f0) = n, deg(f1) = n-1.
       We explicitly compute the pseudoremainder of f0 modulo f1:
       f0 := f1[n-1]*f0 - f0[n]*x*f1
       f0 := f0[n-1]*f1 - f1[n-1]*f0
    */
    fmpz_set(c, f0+n);
    _fmpz_vec_scalar_mul_fmpz(f0, f0, n, f1+n-1);
    _fmpz_vec_scalar_submul_fmpz(f0+1, f1, n-1, c);
    n--;
    fmpz_set(c, f0+n);
    fmpz_neg(d, f1+n);
    _fmpz_vec_scalar_mul_fmpz(f0, f0, n, d);
    _fmpz_vec_scalar_addmul_fmpz(f0, f1, n, c);

    /* If f0 = 0, we win unless we are insisting on squarefree. */
    if (!force_squarefree && _fmpz_vec_is_zero(f0, n)) return(1);

    /* If we miss any one sign change, we cannot have enough. */
    if (fmpz_sgn(f0+n-1) != sgn0_l) return(0);

    if (n==1) return(1); /* If f0 is a scalar, it is nonzero and we win. */

    /* Extract content from f0; in practice, this seems to do better than
       an explicit subresultant computation. */
    _fmpz_vec_content(c, f0, n);
    _fmpz_vec_scalar_divexact_fmpz(f0, f0, n, c);

    /* Swap f0 with f1. */
    t = f0; f0 = f1; f1 = t;
  }
}



/* Set res to floor(a). */
void fmpq_floor(fmpz_t res, const fmpq_t a) {
  fmpz_fdiv_q(res, fmpq_numref(a), fmpq_denref(a));
};

/* Set res to ceil(a). */
void fmpq_ceil(fmpz_t res, const fmpq_t a) {
  fmpz_cdiv_q(res, fmpq_numref(a), fmpq_denref(a));
};

void fmpz_sqrt_f(fmpz_t res, const fmpz_t a) {
  fmpz_sqrt(res, a);
}

void fmpz_sqrt_c(fmpz_t res, const fmpz_t a) {
  int s = fmpz_is_square(a);
  fmpz_sqrt(res, a);
  if (!s) fmpz_add_ui(res, res, 1);
}

/* Set res to floor(a + b sqrt(q)).
   For efficiency, we do not assume a and b are canonical;
   we must thus be careful about signs. */
void fmpq_floor_quad(fmpz_t res, fmpq_t a,
		     fmpq_t b, const fmpz_t q) {
  if (b==NULL) fmpq_floor(res, a);
  else {
    fmpz *anum = fmpq_numref(a);
    fmpz *aden = fmpq_denref(a);
    int aden_s = fmpz_sgn(aden);
    fmpz *bnum = fmpq_numref(b);
    int bnum_s = fmpz_sgn(bnum);
    fmpz *bden = fmpq_denref(b);
    int bden_s = fmpz_sgn(bden);

    fmpz_mul(res, aden, bnum);
    fmpz_mul(res, res, res);
    fmpz_mul(res, res, q);
    if (bnum_s*bden_s >= 0) fmpz_sqrt_f(res, res);
    else {
      fmpz_sqrt_c(res, res);
      fmpz_neg(res, res);
    }
    fmpz_mul_si(res, res, aden_s*bden_s);
    fmpz_addmul(res, anum, bden);
    if (bden_s > 0) fmpz_fdiv_q(res, res, aden);
    else fmpz_cdiv_q(res, res, aden);
    fmpz_fdiv_q(res, res, bden);
  }
}

/* Set res to ceil(a + b sqrt(q)). */
void fmpq_ceil_quad(fmpz_t res, fmpq_t a,
		     fmpq_t b, const fmpz_t q) {
  if (b==NULL) fmpq_ceil(res, a);
  else {
    fmpz *anum = fmpq_numref(a);
    fmpz *aden = fmpq_denref(a);
    int aden_s = fmpz_sgn(aden);
    fmpz *bnum = fmpq_numref(b);
    int bnum_s = fmpz_sgn(bnum);
    fmpz *bden = fmpq_denref(b);
    int bden_s = fmpz_sgn(bden);

    fmpz_mul(res, aden, bnum);
    fmpz_mul(res, res, res);
    fmpz_mul(res, res, q);
    if (bnum_s*bden_s >= 0) fmpz_sqrt_c(res, res);
    else {
      fmpz_sqrt_f(res, res);
      fmpz_neg(res, res);
    }
    fmpz_mul_si(res, res, aden_s*bden_s);
    fmpz_addmul(res, anum, bden);
    if (bden_s > 0) fmpz_cdiv_q(res, res, aden);
    else fmpz_fdiv_q(res, res, aden);
    fmpz_cdiv_q(res, res, bden);
  }
}

/* Memory allocation and initialization. */
ps_static_data_t *ps_static_init(int d, fmpz_t q, int coeffsign, fmpz_t lead,
				 int cofactor, fmpz *modlist, long node_limit,
				 int force_squarefree) {
  int i, j, k, l;
  ps_static_data_t *st_data;
  fmpz_poly_t pol;
  fmpz_t m, const1;
  fmpq *k1;

  fmpz_poly_init(pol);
  fmpz_init(m);
  fmpz_init_set_ui(const1, 1);

  st_data = (ps_static_data_t *)malloc(sizeof(ps_static_data_t));

  st_data->d = d;
  st_data->sign = coeffsign;
  fmpz_init(st_data->q);
  fmpz_set(st_data->q, q);
  st_data->node_limit = node_limit;
  st_data->force_squarefree = force_squarefree;

  fmpz_init(st_data->lead);
  fmpz_set(st_data->lead, lead);

  st_data->cofactor = _fmpz_vec_init(3);
  switch (cofactor) {
  case 0: /* Cofactor 1 */
    fmpz_set_si(st_data->cofactor, 1);
    fmpz_set_si(st_data->cofactor+1, 0);
    fmpz_set_si(st_data->cofactor+2, 0);
    break;

  case 1: /* Cofactor x+sqrt(q) */
    fmpz_set(st_data->cofactor, st_data->q);
    fmpz_sqrt(st_data->cofactor, st_data->cofactor);
    fmpz_set_si(st_data->cofactor+1, 1);
    fmpz_set_si(st_data->cofactor+2, 0);
    break;

  case 2:  /* Cofactor x-sqrt(q) */
    fmpz_set(st_data->cofactor, st_data->q);
    fmpz_sqrt(st_data->cofactor, st_data->cofactor);
    fmpz_neg(st_data->cofactor, st_data->cofactor);
    fmpz_set_si(st_data->cofactor+1, 1);
    fmpz_set_si(st_data->cofactor+2, 0);
    break;

  case 3: /* Cofactor x^2-q */
    fmpz_neg(st_data->cofactor, st_data->q);
    fmpz_set_si(st_data->cofactor+1, 0);
    fmpz_set_si(st_data->cofactor+2, 1);
    break;
  }

  st_data->modlist = _fmpz_vec_init(d+1);
  st_data->f = _fmpq_vec_init(d+1);
  for (i=0; i<=d; i++) {
    fmpz_set(st_data->modlist+i, modlist+d-i);
    fmpq_set_si(st_data->f+i, d-i, 1);
    fmpq_div_fmpz(st_data->f+i, st_data->f+i, st_data->lead);
    /* In order to apply power sums and Descartes' rule of signs
       when the modulus is 0, we must pretend that the modulus is 1. */
    if (!fmpz_is_zero(st_data->modlist+i))
      fmpq_mul_fmpz(st_data->f+i, st_data->f+i, st_data->modlist+i);
  }

  fmpz_mat_init(st_data->binom_mat, d+1, d+1);
  for (i=0; i<=d; i++)
    for (j=0; j<=d; j++)
      fmpz_bin_uiui(fmpz_mat_entry(st_data->binom_mat, i, j), i, j);

  st_data->hausdorff_mats = (fmpq_mat_t *)malloc((d+1)*sizeof(fmpq_mat_t));
  for (i=0; i<=d; i++) {

    fmpq_mat_init(st_data->hausdorff_mats[i], 2*d+2, d+1);
    fmpq_mat_zero(st_data->hausdorff_mats[i]);

    for (j=0; j<=i; j++)
      for (k=0; k<=i; k++) {
	// The coefficient of t^k in (t-2 sqrt(q))^j (t+2 sqrt(q))^{i-j}, rounding down the exponent of q.
	if ((i-k)%2==0)
	  k1 = fmpq_mat_entry(st_data->hausdorff_mats[i], 2*j, k);
	else
	  k1 = fmpq_mat_entry(st_data->hausdorff_mats[i], 2*j+1, k);
	for (l=0; l<=j; l++) if (k-l >=0 && k-l<=i-j) {
	    fmpz_mul(m, fmpz_mat_entry(st_data->binom_mat, j, l),
		     fmpz_mat_entry(st_data->binom_mat, i-j, k-l));
	    if ((j-l)%2==1) fmpz_neg(m, m);
	    fmpq_add_fmpz(k1, k1, m);
	  }
	fmpq_mul_2exp(k1, k1, i-k);
	for (l=0; l<(i-k)/2; l++) fmpq_mul_fmpz(k1, k1, q);
      }
  }

  st_data->sum_mats = (fmpq_mat_t *)malloc((d+1)*sizeof(fmpq_mat_t));
  for (i=0; i<=d; i++) {

    fmpq_mat_init(st_data->sum_mats[i], 1, d+1);
    fmpq_mat_zero(st_data->sum_mats[i]);

    arith_chebyshev_t_polynomial(pol, i);
    for (j=0; j<=d; j++) {

      /* Coefficients of 2*(i-th Chebyshev polynomial)(x/2).
         If q != 1, the coeff of x^j is multiplied by q^{floor(i-j)/2}. */
      if (j <= i) {
	k1 = fmpq_mat_entry(st_data->sum_mats[i], 0, j);
	fmpq_set_fmpz_frac(k1, fmpz_poly_get_coeff_ptr(pol, j), const1);
	fmpz_mul_2exp(m, const1, j);
	fmpq_div_fmpz(k1, k1, m);
	fmpz_set_ui(m, 2);
	fmpq_mul_fmpz(k1, k1, m);
	if (!fmpz_is_one(st_data->q) && i%2==j%2) {
	  fmpz_set(m, st_data->q);
	  fmpz_pow_ui(m, m, (i-j)/2);
	  fmpq_mul_fmpz(k1, k1, m);
	}
      }

    }
  }

  fmpz_poly_clear(pol);
  fmpz_clear(m);
  fmpz_clear(const1);

  return(st_data);
}

ps_dynamic_data_t *ps_dynamic_init(int d, fmpz_t q, fmpz *coefflist) {
  ps_dynamic_data_t *dy_data;
  int i;

  dy_data = (ps_dynamic_data_t *)malloc(sizeof(ps_dynamic_data_t));
  dy_data->d = d;
  dy_data->q_is_1 = fmpz_is_one(q);

  /* Initialize mutable quantities */
  dy_data->n = d;
  dy_data->node_count = 0;
  dy_data->ascend = 0;
  dy_data->pol = _fmpz_vec_init(d+1);
  dy_data->sympol = _fmpz_vec_init(2*d+3);
  if (coefflist != NULL) {
    dy_data->flag = 1; // Activate this process
    for (i=0; i<=d; i++)
      fmpz_set(dy_data->pol+i, coefflist+i);
  } else dy_data->flag = 0;

  fmpq_mat_init(dy_data->power_sums, d+1, 1);
  fmpq_set_si(fmpq_mat_entry(dy_data->power_sums, 0, 0), d, 1);
  fmpq_mat_init(dy_data->hankel_mat, d/2+1, d/2+1);
  fmpq_mat_init(dy_data->hankel_dets, d/2+1, 1);
  fmpq_set_si(fmpq_mat_entry(dy_data->hankel_dets, 0, 0), d, 1);
  fmpq_mat_init(dy_data->hausdorff_prod, 2*d+2, 1);
  fmpq_mat_init(dy_data->hausdorff_sums1, d+1, d+1);
  fmpq_mat_init(dy_data->hausdorff_sums2, d+1, d+1);

  dy_data->upper = _fmpz_vec_init(d+1);

  /* Allocate scratch space */
  fmpq_mat_init(dy_data->sum_prod, 1, 1);
  dy_data->wlen = 3*d+10;
  dy_data->w = _fmpz_vec_init(dy_data->wlen);
  dy_data->w2len = 5;
  dy_data->w2 = _fmpq_vec_init(dy_data->w2len);
  return(dy_data);
}

/* Split off a subtree.
   The first process gives up on the current branch, up to the first coefficient that is not uniquely specified;
   the remaining work is yielded to the second process, which may in turn be split immediately.
*/
void ps_dynamic_split(ps_dynamic_data_t *dy_data, ps_dynamic_data_t *dy_data2) {
  if ((dy_data == NULL) || (dy_data->flag <= 0) || dy_data2->flag) return;

  int i, d = dy_data->d, n = dy_data->n, ascend = dy_data->ascend;

  for (i=d; i>n+ascend; i--)
    if (fmpz_cmp(dy_data->pol+i, dy_data->upper+i) <0) {
      dy_data2->n = n;
      dy_data2->ascend = ascend;
      _fmpz_vec_set(dy_data2->pol, dy_data->pol, d+1);
      _fmpz_vec_set(dy_data2->upper, dy_data->upper, d+1);
      fmpq_mat_set(dy_data2->power_sums, dy_data->power_sums);
      fmpq_mat_set(dy_data2->hankel_dets, dy_data->hankel_dets);
      if (dy_data->q_is_1) {
	fmpq_mat_set(dy_data2->hausdorff_sums1, dy_data->hausdorff_sums1);
	fmpq_mat_set(dy_data2->hausdorff_sums2, dy_data->hausdorff_sums2);
      }
      fmpz_set(dy_data2->upper+i, dy_data2->pol+i);
      dy_data->ascend = i-n;
      dy_data2->flag = 1; // This process can now itself be split.
      return;
  }
  return;
}

/* Memory deallocation. */
void ps_static_clear(ps_static_data_t *st_data) {
  if (st_data == NULL) return;
  int i, d = st_data->d;
  fmpz_clear(st_data->lead);
  fmpz_clear(st_data->q);
  _fmpz_vec_clear(st_data->cofactor, 3);
  fmpz_mat_clear(st_data->binom_mat);
  _fmpq_vec_clear(st_data->f, d+1);
  _fmpz_vec_clear(st_data->modlist, d+1);
  for (i=0; i<=d; i++)  {
    fmpq_mat_clear(st_data->hausdorff_mats[i]);
    fmpq_mat_clear(st_data->sum_mats[i]);
  }
  free(st_data->hausdorff_mats);
  free(st_data->sum_mats);
  free(st_data);
}

void ps_dynamic_clear(ps_dynamic_data_t *dy_data) {
  if (dy_data == NULL) return;
  int d = dy_data->d;
  _fmpz_vec_clear(dy_data->pol, d+1);
  _fmpz_vec_clear(dy_data->sympol, 2*d+3);
  _fmpz_vec_clear(dy_data->upper, d+1);
  fmpq_mat_clear(dy_data->power_sums);
  fmpq_mat_clear(dy_data->sum_prod);
  fmpq_mat_clear(dy_data->hankel_mat);
  fmpq_mat_clear(dy_data->hankel_dets);
  fmpq_mat_clear(dy_data->hausdorff_prod);
  fmpq_mat_clear(dy_data->hausdorff_sums1);
  fmpq_mat_clear(dy_data->hausdorff_sums2);
  _fmpz_vec_clear(dy_data->w, dy_data->wlen);
  _fmpq_vec_clear(dy_data->w2, dy_data->w2len);
  free(dy_data);
}

/* Subroutines to adjust lower and upper bounds within set_range_from_power_sums.
   These use t0z, t0q, t4q as persistent scratch space.
   The pair (val1, val2) stands for val1 + val2*sqrt(q);
   passing NULL for val2 is a faster variant of passing 0.

  Usage: if g is a monic linear function of the k-th power sum, then
  set_upper(g) or change_upper(g) imposes the condition g >= 0;
  set_lower(g) or change_lower(g) imposes the condition g <= 0.

*/

#define STATE lower, upper, q, f, t0z, t0q, t4q
#define STATE_DECLARE fmpz_t lower, fmpz_t upper, fmpz_t q, fmpq_t f, fmpz_t t0z, fmpq_t t0q, fmpq_t t4q

void set_lower(const fmpq_t val1, const fmpq_t val2, STATE_DECLARE) {
  fmpq_div(t0q, val1, f);
  if (val2==NULL) fmpq_ceil(lower, t0q);
  else {
    fmpq_div(t4q, val2, f);
    fmpq_ceil_quad(lower, t0q, t4q, q);
  }
}

void set_upper(const fmpq_t val1, const fmpq_t val2, STATE_DECLARE) {
  fmpq_div(t0q, val1, f);
  if (val2==NULL) fmpq_floor(upper, t0q);
  else {
    fmpq_div(t4q, val2, f);
    fmpq_floor_quad(upper, t0q, t4q, q);
  }
}

void change_lower(const fmpq_t val1, const fmpq_t val2, STATE_DECLARE) {
  fmpq_div(t0q, val1, f);
  if (val2==NULL) fmpq_ceil(t0z, t0q);
  else {
    fmpq_div(t4q, val2, f);
    fmpq_ceil_quad(t0z, t0q, t4q, q);
  }
  if (fmpz_cmp(t0z, lower) > 0) fmpz_set(lower, t0z);
}

void change_upper(const fmpq_t val1, const fmpq_t val2, STATE_DECLARE) {
  fmpq_div(t0q, val1, f);
  if (val2==NULL) fmpq_floor(t0z, t0q);
  else {
    fmpq_div(t4q, val2, f);
    fmpq_floor_quad(t0z, t0q, t4q, q);
  }
  if (fmpz_cmp(t0z, upper) < 0) fmpz_set(upper, t0z);
}

void change_lower_strict(const fmpq_t val1, const fmpq_t val2, STATE_DECLARE) {
  fmpq_div(t0q, val1, f);
  if (val2==NULL) fmpq_floor(t0z, t0q);
  else {
    fmpq_div(t4q, val2, f);
    fmpq_floor_quad(t0z, t0q, t4q, q);
  }
  fmpz_add_ui(t0z, t0z, 1);
  if (fmpz_cmp(t0z, lower) > 0) fmpz_set(lower, t0z);
}

void change_upper_strict(const fmpq_t val1, const fmpq_t val2, STATE_DECLARE) {
  fmpq_div(t0q, val1, f);
  if (val2==NULL) fmpq_ceil(t0z, t0q);
  else {
    fmpq_div(t4q, val2, f);
    fmpq_ceil_quad(t0z, t0q, t4q, q);
  }
  fmpz_sub_ui(t0z, t0z, 1);
  if (fmpz_cmp(t0z, upper) < 0) fmpz_set(upper, t0z);
}

/* Impose the condition that val1*val3 >= val2, assuming that val1 is a linear 
   monic function of the k-th power sum and val2, val3 do not depend on this sum. */
void impose_quadratic_condition(const fmpq_t val1, const fmpq_t val2,
				const fmpq_t val3, STATE_DECLARE) {
  int s = fmpq_sgn(val3);
  if (s) {
    fmpq_mul(t0q, val2, val2);
    fmpq_div(t0q, t0q, val3);
    fmpq_sub(t0q, val1, t0q);
    if (s>0) change_upper(t0q, NULL, STATE);
    else change_lower(t0q, NULL, STATE);
  }
}

/* The following is the key subroutine: given some initial coefficients, compute
   a lower and upper bound for the next coefficient. Return 1 iff the resulting
   interval is nonempty.

*/
int set_range_from_power_sums(ps_static_data_t *st_data,
			      ps_dynamic_data_t *dy_data) {
  int i, j, r;
  int d = st_data->d;
  int n = dy_data->n;
  int k = d+1-n;
  int q_is_1 = dy_data->q_is_1;
  fmpz *modulus = st_data->modlist+n-1;
  fmpz *pol = dy_data->pol;
  fmpz *q = st_data->q;
  fmpq *f = (fmpq *)(st_data->f+n-1);
  fmpq *t;

  /* Allocate temporary variables from persistent scratch space. */
  fmpz *tpol = dy_data->w;
  fmpz *tpol2 = dy_data->w+d+1;

  fmpz *t0z = dy_data->w+3*d+5;
  fmpz *t1z = dy_data->w+3*d+6;
  fmpz *t2z = dy_data->w+3*d+7;
  fmpz *lower = dy_data->w+3*d+8;
  fmpz *upper = dy_data->w+3*d+9;

  fmpq *t0q = dy_data->w2;
  fmpq *t1q = dy_data->w2+1;
  fmpq *t2q = dy_data->w2+2;
  fmpq *t3q = dy_data->w2+3;
  fmpq *t4q = dy_data->w2+4;

  /* If k>d, no further coefficients to bound. */
  if (k>d) return(1);

  /* Update power_sums[k]. */
  t = fmpq_mat_entry(dy_data->power_sums, k, 0);
  fmpq_set_si(t, -k, 1);
  fmpq_mul_fmpz(t, t, pol+d-k);
  for (i=1; i<k; i++) {
    fmpq_set_si(t0q, -1, 1);
    fmpq_mul_fmpz(t0q, t0q, pol+d-i);
    fmpq_addmul(t, t0q, fmpq_mat_entry(dy_data->power_sums, k-i, 0));
  }
  fmpq_div_fmpz(t, t, pol+d);

  /* Condition: the k-th symmetrized power sum must lie in [-2*sqrt(q), 2*sqrt(q)]. */
  fmpq_mat_mul(dy_data->sum_prod, st_data->sum_mats[k], dy_data->power_sums);
  t = fmpq_mat_entry(dy_data->sum_prod, 0, 0);

  fmpq_set_si(t1q, 2*d, 1);
  if (!q_is_1) {
    fmpz_pow_ui(t0z, q, k/2);
    fmpq_mul_fmpz(t1q, t1q, t0z);
  }
  if (k%2==0) {
    fmpq_sub(t0q, t, t1q);
    set_lower(t0q, NULL, STATE);
    fmpq_add(t0q, t, t1q);
    set_upper(t0q, NULL, STATE);
  } else {
    set_upper(t, t1q, STATE);
    fmpq_neg(t1q, t1q);
    set_lower(t, t1q, STATE);
    }

  /* Compute the divided (n-1)-st derivative of pol, answer in tpol. */
  for (i=0; i<=k; i++)
    fmpz_mul(tpol+i, fmpz_mat_entry(st_data->binom_mat, n-1+i, n-1), pol+n-1+i);

  /* Condition: Descartes' rule of signs applies at -2*sqrt(q), +2*sqrt(q).
   This is only a new condition for the evaluations at these points. */

  fmpq_set_si(t3q, -k, 1);
  fmpq_div_fmpz(t3q, t3q, pol+d);

  for (i=0; 2*i <= k; i++) fmpz_mul_2exp(tpol2+i, tpol+2*i, 2*i);
  _fmpz_poly_evaluate_fmpz(t0z, tpol2, (k+2) / 2, q);
  fmpq_mul_fmpz(t1q, t3q, t0z);

  for (i=0; 2*i+1 <= k; i++) fmpz_mul_2exp(tpol2+i, tpol+2*i+1, 2*i+1);
  _fmpz_poly_evaluate_fmpz(t0z, tpol2, (k+1) / 2, q);
  fmpq_mul_fmpz(t2q, t3q, t0z);

  /* If checking for squarefree, shear endpoints off the range. */
  if (st_data->force_squarefree) {
    change_lower_strict(t1q, t2q, STATE);
    fmpq_neg(t2q, t2q);
    if (k%2==1) change_upper_strict(t1q, t2q, STATE);
    else change_lower_strict(t1q, t2q, STATE);
  }
  else {
    change_lower(t1q, t2q, STATE);
    fmpq_neg(t2q, t2q);
    if (k%2==1) change_upper(t1q, t2q, STATE);
    else change_lower(t1q, t2q, STATE);
  }
  if (fmpz_cmp(lower, upper) > 0) return(0);

  /* Update Hankel matrices. */
  if (k%2==0) {
    fmpq_mat_one(dy_data->hankel_mat);
    for (i=0; i<=k/2; i++)
      for (j=0; j<=k/2; j++)
	fmpq_set(fmpq_mat_entry(dy_data->hankel_mat, i, j),
		 fmpq_mat_entry(dy_data->power_sums, i+j, 0));
    fmpq_mat_det(t0q, dy_data->hankel_mat);
    t = fmpq_mat_entry(dy_data->hankel_dets, k/2-1, 0);
    fmpq_set(fmpq_mat_entry(dy_data->hankel_dets, k/2, 0), t0q);
  }

  /* If modulus==0, then return 1 iff [lower, upper] contains 0
     and Rolle's theorem is satisfied.
   */
  if (fmpz_is_zero(modulus)) {
    if ((fmpz_sgn(lower) > 0) || (fmpz_sgn(upper) < 0) ||
	!_fmpz_poly_all_real_roots(tpol, k+1, tpol2, st_data->force_squarefree,
				   NULL, NULL)) return(0);
    fmpz_zero(lower);
    fmpz_zero(upper);
    return(1);
  } else
    if (fmpz_cmp(lower, upper) > 0) return(0);

  /* Condition: nonnegativity of the Hankel determinant.
     TODO: reimplement this as a subresultant. */
  if (k%2==0) {
    if (fmpq_sgn(t) > 0) {
      fmpq_div(t0q, t0q, t);
      change_upper(t0q, NULL, STATE);
      }
    else if (st_data->force_squarefree || fmpq_sgn(t0q)) return(0);
    else change_upper(fmpq_mat_entry(dy_data->power_sums, k, 0), NULL, STATE);
    if (fmpz_cmp(lower, upper) > 0) return(0);
  }

  /* Condition: the Hausdorff moment criterion for having roots in [-2, 2]. 
     This might be redundant given the Hankel and Descartes criteria. */
  fmpq_mat_mul(dy_data->hausdorff_prod, st_data->hausdorff_mats[k], dy_data->power_sums);
  for (i=0; i<=k; i++) {
    fmpq_set(t1q, fmpq_mat_entry(dy_data->hausdorff_prod, 2*i, 0));
    fmpq_set(t2q, fmpq_mat_entry(dy_data->hausdorff_prod, 2*i+1, 0));
    if (i%2==0) change_upper(t1q, t2q, STATE);
    else change_lower(t1q, t2q, STATE);
    if (q_is_1) {
      fmpq_set(fmpq_mat_entry(dy_data->hausdorff_sums1, k, i), t1q);
      fmpq_set(fmpq_mat_entry(dy_data->hausdorff_sums2, k, i), t2q);
    }
  }
  if (fmpz_cmp(lower, upper) > 0) return(0);

  /* Condition: log convexity based on Cauchy-Schwarz. */
  /* TODO: extend to q != 1 without losing too much efficiency. */
  if (q_is_1) {
    for (i=0; i<=k-2; i++) {
      fmpq_add(t1q, fmpq_mat_entry(dy_data->hausdorff_sums1, k, i),
	       fmpq_mat_entry(dy_data->hausdorff_sums2, k, i));
      fmpq_add(t2q, fmpq_mat_entry(dy_data->hausdorff_sums1, k-1, i),
	       fmpq_mat_entry(dy_data->hausdorff_sums2, k-1, i));
      fmpq_add(t3q, fmpq_mat_entry(dy_data->hausdorff_sums1, k-2, i),
	       fmpq_mat_entry(dy_data->hausdorff_sums2, k-2, i));
      impose_quadratic_condition(t1q, t2q, t3q, STATE);
    }
    for (i=2; i<=k; i++) {
      fmpq_add(t1q, fmpq_mat_entry(dy_data->hausdorff_sums1, k, i),
	       fmpq_mat_entry(dy_data->hausdorff_sums2, k, i));
      fmpq_add(t2q, fmpq_mat_entry(dy_data->hausdorff_sums1, k-1, i-1),
	       fmpq_mat_entry(dy_data->hausdorff_sums2, k-1, i-1));
      fmpq_add(t3q, fmpq_mat_entry(dy_data->hausdorff_sums1, k-2, i-2),
	       fmpq_mat_entry(dy_data->hausdorff_sums2, k-2, i-2));
      impose_quadratic_condition(t1q, t2q, t3q, STATE);
    }
  }
  r = fmpz_cmp(lower, upper);
  if (r>0) return(0);

  /* Check the Rolle condition at the midpoint. If it holds, perform a binary
     search on the left endpoint; otherwise, do a linear search. */
  if (r) {
    fmpz_add(t0z, lower, upper);
    fmpz_fdiv_q_2exp(t0z, t0z, 1);
    r = _fmpz_poly_all_real_roots(tpol, k+1, tpol2, st_data->force_squarefree, t0z, modulus);
  }
  if (r) {
    fmpz_set(t2z, t0z);
    while (fmpz_cmp(lower, t0z)) {
      fmpz_add(t1z, lower, t0z);
      fmpz_fdiv_q_2exp(t1z, t1z, 1);
      r = _fmpz_poly_all_real_roots(tpol, k+1, tpol2, st_data->force_squarefree, t1z, modulus);
      if (r) fmpz_set(t0z, t1z);
      else fmpz_add_ui(lower, t1z, 1);
    }
  } else {
    r = _fmpz_poly_all_real_roots(tpol, k+1, tpol2, st_data->force_squarefree, lower, modulus);
    while (!r) {
      fmpz_add_ui(lower, lower, 1);
      if (fmpz_cmp(lower, upper) > 0) return(0);
      r = _fmpz_poly_all_real_roots(tpol, k+1, tpol2, st_data->force_squarefree, lower, modulus);
    }
    if (fmpz_cmp(lower, t0z)<0) fmpz_sub_ui(upper, t0z, 1);
    fmpz_set(t2z, lower);
  }
  /* Now do a binary search on the right endpoint. */
  while (fmpz_cmp(t2z, upper)) {
      fmpz_add(t1z, t2z, upper);
      fmpz_cdiv_q_2exp(t1z, t1z, 1);
      r = _fmpz_poly_all_real_roots(tpol, k+1, tpol2, st_data->force_squarefree, t1z, modulus);
      if (r) fmpz_set(t2z, t1z);
      else fmpz_sub_ui(upper, t1z, 1);
  }

  /* Set the new upper bound. */
  fmpz_mul(upper, upper, modulus);
  fmpz_add(dy_data->upper+n-1, pol+n-1, upper);

  /* Set the new polynomial value. */
  fmpz_addmul(pol+n-1, lower, modulus);

  /* Correct the k-th power sum and related quantities. */
  t1q = fmpq_mat_entry(dy_data->power_sums, k, 0);
  fmpq_mul_fmpz(t0q, f, lower);
  fmpq_sub(t1q, t1q, t0q);
  if (q_is_1) for (i=0; i<=k; i++) {
      t1q = fmpq_mat_entry(dy_data->hausdorff_sums1, k, i);
      fmpq_sub(t1q, t1q, t0q);
    }
  if (k%2==0) {
    t1q = fmpq_mat_entry(dy_data->hankel_dets, k/2, 0);
    fmpq_submul(t1q, fmpq_mat_entry(dy_data->hankel_dets, k/2-1, 0), t0q);
  }
  return(1);
}

/* Increment the current moving counter and update stored data to match. */
void step_forward(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int n) {
  int d = st_data->d, k = d-n;
  fmpz *pol = dy_data->pol;
  fmpq *tq = fmpq_mat_entry(dy_data->power_sums, k, 0);
  int j;

  fmpz_add(pol+n, pol+n, st_data->modlist+n);
  fmpq_sub(tq, tq, st_data->f+n);
  if (dy_data->q_is_1) for (j=0; j<=k; j++) {
      tq = fmpq_mat_entry(dy_data->hausdorff_sums1, k, j);
      fmpq_sub(tq, tq, st_data->f+n);
    }
  if (k%2==0)
    fmpq_submul(fmpq_mat_entry(dy_data->hankel_dets, k/2, 0),
		st_data->f+n, fmpq_mat_entry(dy_data->hankel_dets, k/2-1, 0));
}

/* Return value sent back in dy_data->flag:
   1: in process
   2: found a solution
   0: tree exhausted
   -1: maximum number of nodes reached
*/

void next_pol(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int max_steps) {
  int d = st_data->d;
  int node_limit = st_data->node_limit;
  fmpz *modlist = st_data->modlist;

  int ascend = dy_data->ascend;
  int n = dy_data->n;
  int q_is_1 = dy_data->q_is_1;
  long node_count = dy_data->node_count;
  fmpz *upper = dy_data->upper;
  fmpz *pol = dy_data->pol;
  fmpz *sympol = dy_data->sympol;
  fmpz *temp = dy_data->w;

  int i, j, flag = 1, count_steps = 0;

  if (dy_data==NULL || !dy_data->flag) return; // No work assigned to this process
  if (n>d) return;

  dy_data->flag = 0; // Prevent work-stealing while this process is running

  while ((flag==1) && (count_steps <= max_steps)) {
    count_steps += 1;
    if (ascend) { // Ascend the tree and step forward as needed.
      n += ascend;
      if (n>d) flag = 0; // This process is complete.
      else {
	ascend = (fmpz_is_zero(modlist+n) || (fmpz_cmp(pol+n, upper+n) >= 0));
	if (!ascend) step_forward(st_data, dy_data, n);
      }
    } else if (n < 0) { // Return a solution.
      _fmpz_vec_zero(sympol, 2*d+3);
      for (i=0; i<=d; i++) {
	fmpz_one(temp);
	for (j=0; j<=i; j++) {
	  fmpz_addmul(sympol+d+i-2*j, pol+i, temp);
	  if (j<i) {
	    fmpz_mul(temp, temp, st_data->q);
	    fmpz_mul_si(temp, temp, i-j);
	    fmpz_divexact_si(temp, temp, j+1);
	  }
	}
      }
      _fmpz_vec_scalar_mul_si(sympol, sympol, 2*d+1, st_data->sign);
      _fmpz_poly_mul_KS(sympol, sympol, 2*d+1, st_data->cofactor, 3);
      ascend = 1;
      flag = 2;
    } else { // Compute children of the current node.
      dy_data->n = n;
      ascend = !set_range_from_power_sums(st_data, dy_data);
      n -= 1;
      if (ascend) {
	node_count += 1;
	if (node_limit != -1 && node_count >= node_limit) flag = -1;
      }
    }
  }
  dy_data->ascend = ascend;
  dy_data->n = n;
  dy_data->node_count = node_count;
  dy_data->flag = flag;
}
