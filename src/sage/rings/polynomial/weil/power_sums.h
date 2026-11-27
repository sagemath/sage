#pragma once
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>

#if defined(_OPENMP)
  #include <omp.h>
#endif

typedef struct ps_static_data {
  int d, force_squarefree;
  long node_limit;
  fmpz_t q, q_sqrt, c0, c1;
  int q_is_1, q_is_square, lead_is_1, num_constraints;
  int *constraint_lens;
  fmpz *modlist, *binom_mat, *sum_mats, *eval_pm2_mats, *ranges, *lead_pows, *constraints;
  fmpz_mat_t pol_to_sym;
} ps_static_data_t;

typedef struct ps_dynamic_data {
  int d, n, ascend, flag;
  long node_count;
  fmpz *pol, *sympol, *upper, *power_sums_num, *hankel_dets;

  /* Scratch space */
  fmpz *w;
  long wlen; /* = 3*d+9 */
} ps_dynamic_data_t;

int num_threads();
int is_mpz(fmpz f);

ps_static_data_t *ps_static_init(int d, const fmpz_t q, const fmpz_t lead, const fmpz *modlist, int num_constraints, const fmpz *constraints, long node_limit, int force_squarefree);
ps_dynamic_data_t *ps_dynamic_init(int d, fmpz *coefflist);
void ps_static_clear(ps_static_data_t *st_data);
void ps_dynamic_clear(ps_dynamic_data_t *dy_data);
void ps_cleanup(int n);

int ascend_step_forward(const ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int n);
int reduce_range_from_rolle(const ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int n);
int set_range_from_power_sums(const ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int n);
void reciprocal_transform(const ps_static_data_t *st_data, ps_dynamic_data_t *dy_data);

int ps_dynamic_split(const ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, ps_dynamic_data_t *dy_data2);
int ps_next_pol(const ps_static_data_t *st_data, ps_dynamic_data_t *dy_data, int max_steps);

