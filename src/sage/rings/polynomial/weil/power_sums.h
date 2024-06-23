/* sage_setup: distribution = sagemath-flint
 */
#pragma once
#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/arith.h>

typedef struct ps_static_data {
  int d, sign, force_squarefree;
  long node_limit;
  fmpz_t a, b, lead, q;
  fmpz_mat_t binom_mat;
  fmpz *cofactor;
  fmpz *modlist;
  fmpq_mat_t *hausdorff_mats;
  fmpq_mat_t *sum_mats;
  fmpq *f;
} ps_static_data_t;

typedef struct ps_dynamic_data {
  int d, n, ascend, flag, q_is_1;
  long node_count;
  fmpq_mat_t power_sums, sum_prod, hankel_mat, hankel_dets,
    hausdorff_prod, hausdorff_sums1, hausdorff_sums2;
  fmpz *pol, *sympol, *upper;

  /* Scratch space */
  fmpz *w;
  long wlen; /* = 4*d+10 */
  fmpq *w2;
  long w2len; /* = 5 */
} ps_dynamic_data_t;

int has_openmp();
ps_static_data_t *ps_static_init(int d, fmpz_t q, int coeffsign, fmpz_t lead,
				 int cofactor, fmpz *modlist, long node_limit,
				 int force_squarefree);
ps_dynamic_data_t *ps_dynamic_init(int d, fmpz_t q, fmpz *coefflist);
void ps_static_clear(ps_static_data_t *st_data);
void ps_dynamic_clear(ps_dynamic_data_t *dy_data);
void extract_pol(int *Q, ps_dynamic_data_t *dy_data);
void ps_dynamic_split(ps_dynamic_data_t *dy_data, ps_dynamic_data_t *dy_data2);
void step_forward(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int n);
void next_pol(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int max_steps);
