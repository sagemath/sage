/* sage_setup: distribution = sagemath-modules
 */
/******************************************************************************
*
*                        DGS - Discrete Gaussian Samplers
*
* Copyright (c) 2014, Martin Albrecht  <martinralbrecht+dgs@googlemail.com>
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
* The views and conclusions contained in the software and documentation are
* those of the authors and should not be interpreted as representing official
* policies, either expressed or implied, of the FreeBSD Project.
******************************************************************************/

#include "dgs.h"
#include <assert.h>
#include <stdlib.h>
#include <limits.h>

/** SIGMA2 **/

dgs_disc_gauss_sigma2p_t *dgs_disc_gauss_sigma2p_init(void) {
  dgs_disc_gauss_sigma2p_t *self = (dgs_disc_gauss_sigma2p_t*)calloc(sizeof(dgs_disc_gauss_sigma2p_t),1);
  if (!self) dgs_die("out of memory");
  self->B = dgs_bern_uniform_init(0);
  return self;
}

void dgs_disc_gauss_sigma2p_mp_call(mpz_t rop, dgs_disc_gauss_sigma2p_t *self, gmp_randstate_t state) {
  while(1) {
    if (!dgs_bern_uniform_call(self->B, state)) {
      mpz_set_ui(rop, 0);
      return;
    }
    int dobreak = 0;
    for(unsigned long i=1; ;i++) {
      for(size_t j=0; j<2*i-2; j++) {
        if(dgs_bern_uniform_call(self->B, state)) {
          dobreak = 1;
          break;
        }
      }
      if (__DGS_LIKELY(dobreak))
        break;
      if (!dgs_bern_uniform_call(self->B, state)) {
        mpz_set_ui(rop, i);
        return;
      }
    }
  }
}

long dgs_disc_gauss_sigma2p_dp_call(dgs_disc_gauss_sigma2p_t *self) {
  while(1) {
    if (!dgs_bern_uniform_call_libc(self->B)) {
      return 0;
    }
    int dobreak = 0;
    for(unsigned long i=1; ;i++) {
      for(size_t j=0; j<2*i-2; j++) {
        if(dgs_bern_uniform_call_libc(self->B)) {
          dobreak = 1;
          break;
        }
      }
      if (__DGS_LIKELY(dobreak))
        break;
      if (!dgs_bern_uniform_call_libc(self->B)) {
        return i;
      }
    }
  }
}


void dgs_disc_gauss_sigma2p_clear(dgs_disc_gauss_sigma2p_t *self) {
  assert(self != NULL);
  if (self->B) dgs_bern_uniform_clear(self->B);
  free(self);
}

/** GENERAL SIGMA :: INIT **/

static inline void _dgs_disc_gauss_mp_init_f(mpfr_t f, const mpfr_t sigma) {
  mpfr_init2(f, mpfr_get_prec(sigma));
  mpfr_set(f, sigma, MPFR_RNDN);
  mpfr_sqr(f, f, MPFR_RNDN); // f = σ²
  mpfr_mul_ui(f, f, 2, MPFR_RNDN); // f = 2 σ²
  mpfr_ui_div(f, 1, f, MPFR_RNDN);  // f = 1/(2 σ²)
  mpfr_neg(f, f, MPFR_RNDN); // f = -1/(2 σ²)
}

static inline void _dgs_disc_gauss_mp_init_upper_bound(mpz_t upper_bound,
                                                       mpz_t upper_bound_minus_one,
                                                       mpz_t two_upper_bound_minus_one,
                                                       mpfr_t sigma, size_t tailcut) {
  mpfr_t tmp;
  mpfr_init2(tmp, mpfr_get_prec(sigma));
  mpz_init(upper_bound);
  mpz_init(upper_bound_minus_one);
  mpz_init(two_upper_bound_minus_one);
  mpfr_mul_ui(tmp, sigma, tailcut, MPFR_RNDN); // tmp = σ·τ
  mpfr_add_ui(tmp, tmp, 1, MPFR_RNDN); // tmp = σ·τ + 1
  mpfr_get_z(upper_bound, tmp, MPFR_RNDU); // upper_bound = ⌈σ·τ + 1⌉
  mpz_sub_ui(upper_bound_minus_one, upper_bound, 1); // upper_bound - 1 = ⌈σ·τ⌉
  mpz_mul_ui(two_upper_bound_minus_one, upper_bound, 2);
  mpz_sub_ui(two_upper_bound_minus_one, two_upper_bound_minus_one, 1); // 2·upper_bound - 1
  mpfr_clear(tmp);
}

static inline void _dgs_disc_gauss_mp_init_bexp(dgs_disc_gauss_mp_t *self, mpfr_t sigma, mpz_t upper_bound) {
  mpfr_init2(self->f, mpfr_get_prec(sigma));
  mpfr_set(self->f, sigma, MPFR_RNDN); // f = σ
  mpfr_sqr(self->f, self->f, MPFR_RNDN); // f = σ²
  mpfr_mul_ui(self->f, self->f, 2, MPFR_RNDN); // f = 2 σ²
  size_t l = 2*mpz_sizeinbase(upper_bound, 2);
  self->Bexp = dgs_bern_exp_mp_init(self->f, l);
}

dgs_disc_gauss_mp_t *dgs_disc_gauss_mp_init(mpfr_t sigma, mpfr_t c, size_t tau, dgs_disc_gauss_alg_t algorithm) {
  if (mpfr_cmp_ui(sigma,0)<= 0)
    dgs_die("sigma must be > 0");
  if (tau == 0)
    dgs_die("tau must be > 0");

  mpfr_prec_t prec = mpfr_get_prec(sigma);
  if (mpfr_get_prec(c) > prec)
    prec = mpfr_get_prec(c);

  dgs_disc_gauss_mp_t *self = (dgs_disc_gauss_mp_t*)calloc(sizeof(dgs_disc_gauss_mp_t),1);
  if (!self) dgs_die("out of memory");

  mpz_init(self->x);
  mpz_init(self->x2);
  mpz_init(self->k);
  mpfr_init2(self->y, prec);
  mpfr_init2(self->z, prec);

  mpfr_init2(self->sigma, prec);
  mpfr_set(self->sigma, sigma, MPFR_RNDN);

  mpfr_init2(self->c, prec);
  mpfr_set(self->c, c, MPFR_RNDN);
  mpz_init(self->c_z);
  mpfr_get_z(self->c_z, c, MPFR_RNDN);
  mpfr_init2(self->c_r, prec);
  mpfr_sub_z(self->c_r, self->c, self->c_z, MPFR_RNDN);

  self->tau = tau;

  switch(algorithm) {

  case DGS_DISC_GAUSS_UNIFORM_ONLINE: {
    _dgs_disc_gauss_mp_init_upper_bound(self->upper_bound,
                                        self->upper_bound_minus_one,
                                        self->two_upper_bound_minus_one,
                                        self->sigma, self->tau);

    self->call = dgs_disc_gauss_mp_call_uniform_online;
    _dgs_disc_gauss_mp_init_f(self->f, self->sigma);

   break;
  }
  case DGS_DISC_GAUSS_UNIFORM_TABLE: {
    _dgs_disc_gauss_mp_init_upper_bound(self->upper_bound,
                                        self->upper_bound_minus_one,
                                        self->two_upper_bound_minus_one,
                                        self->sigma, self->tau);

    self->B = dgs_bern_uniform_init(0);
    _dgs_disc_gauss_mp_init_f(self->f, sigma);

    if (mpfr_zero_p(self->c_r)) { /* c is an integer */
      self->call = dgs_disc_gauss_mp_call_uniform_table;
      if (mpz_cmp_ui(self->upper_bound, ULONG_MAX/sizeof(mpfr_t))>0){
        dgs_disc_gauss_mp_clear(self);
        dgs_die("integer overflow");
      }
      self->rho = (mpfr_t*)malloc(sizeof(mpfr_t)*mpz_get_ui(self->upper_bound));
      if (!self->rho){
        dgs_disc_gauss_mp_clear(self);
        dgs_die("out of memory");
      }

      mpfr_t x_;
      mpfr_init2(x_, prec);
      for(unsigned long x=0; x<mpz_get_ui(self->upper_bound); x++) {
        mpfr_set_ui(x_, x, MPFR_RNDN);
        mpfr_sqr(x_, x_, MPFR_RNDN);
        mpfr_mul(x_, x_, self->f, MPFR_RNDN);
        mpfr_exp(x_, x_, MPFR_RNDN);
        mpfr_init2(self->rho[x], prec);
        mpfr_set(self->rho[x], x_, MPFR_RNDN);
      }
      mpfr_div_ui(self->rho[0],self->rho[0], 2, MPFR_RNDN);
      mpfr_clear(x_);

    } else { /* c is not an integer, we need a bigger table as our nice symmetry is lost */
      self->call = dgs_disc_gauss_mp_call_uniform_table_offset;
      if (mpz_cmp_ui(self->two_upper_bound_minus_one, ULONG_MAX/sizeof(mpfr_t)) > 0){
        dgs_disc_gauss_mp_clear(self);
        dgs_die("integer overflow");
      }
      // we need a bigger table
      self->rho = (mpfr_t*)malloc(sizeof(mpfr_t)*mpz_get_ui(self->two_upper_bound_minus_one));
      if (!self->rho){
        dgs_disc_gauss_mp_clear(self);
        dgs_die("out of memory");
      }

      mpfr_t x_;
      mpfr_init2(x_, prec);
      long absmax = mpz_get_ui(self->upper_bound) - 1;
      for(long x=-absmax; x<=absmax; x++) {
        mpfr_set_si(x_, x, MPFR_RNDN);
        mpfr_sub(x_, x_, self->c_r, MPFR_RNDN);
        mpfr_sqr(x_, x_, MPFR_RNDN);
        mpfr_mul(x_, x_, self->f, MPFR_RNDN);
        mpfr_exp(x_, x_, MPFR_RNDN);
        mpfr_init2(self->rho[x+absmax], prec);
        mpfr_set(self->rho[x+absmax], x_, MPFR_RNDN);
      }
      mpfr_clear(x_);
    }
  }
    break;

  case DGS_DISC_GAUSS_UNIFORM_LOGTABLE: {
    self->call = dgs_disc_gauss_mp_call_uniform_logtable;
    _dgs_disc_gauss_mp_init_upper_bound(self->upper_bound,
                                        self->upper_bound_minus_one,
                                        self->two_upper_bound_minus_one,
                                        self->sigma, self->tau);

    if (!mpfr_zero_p(self->c_r)) {
      dgs_disc_gauss_mp_clear(self);
      dgs_die("algorithm DGS_DISC_GAUSS_UNIFORM_LOGTABLE requires c%1 == 0");
    }

    _dgs_disc_gauss_mp_init_bexp(self, self->sigma, self->upper_bound);
   break;
  }

  case DGS_DISC_GAUSS_SIGMA2_LOGTABLE: {
    self->call = dgs_disc_gauss_mp_call_sigma2_logtable;

    if (!mpfr_zero_p(self->c_r)) {
      dgs_disc_gauss_mp_clear(self);
      dgs_die("algorithm DGS_DISC_GAUSS_SIGMA2_LOGTABLE requires c%1 == 0");
    }

    mpfr_t tmp;
    mpfr_init2(tmp, prec);

    mpfr_t sigma2;
    mpfr_init2(sigma2, prec);
    mpfr_set_ui(sigma2, 2, MPFR_RNDN); // 2
    mpfr_log(sigma2, sigma2, MPFR_RNDN); //log₂ 2
    mpfr_mul_ui(sigma2, sigma2, 2, MPFR_RNDN); //2·log₂ 2
    mpfr_ui_div(sigma2, 1, sigma2, MPFR_RNDN); //1/(2·log₂ 2)
    mpfr_sqrt(sigma2, sigma2, MPFR_RNDN); //σ₂ = sqrt(1/(2·log₂ 2))
    mpfr_div(tmp, sigma, sigma2, MPFR_RNDN);
    mpfr_get_z(self->k, tmp, MPFR_RNDN);
    mpfr_mul_z(self->sigma, sigma2, self->k, MPFR_RNDN); //k·σ₂
    mpfr_clear(sigma2);
    mpfr_clear(tmp);

    _dgs_disc_gauss_mp_init_upper_bound(self->upper_bound,
                                        self->upper_bound_minus_one,
                                        self->two_upper_bound_minus_one,
                                        self->sigma, self->tau);

    _dgs_disc_gauss_mp_init_bexp(self, self->sigma, self->upper_bound);
    self->B = dgs_bern_uniform_init(0);
    self->D2 = dgs_disc_gauss_sigma2p_init();
    break;
  }

  default:
    free(self);
    dgs_die("unknown algorithm %d", algorithm);
  }
  return self;
}

/** GENERAL SIGMA :: CALL **/

void dgs_disc_gauss_mp_call_uniform_table(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state) {
  unsigned long x;
  do {
    mpz_urandomm(self->x, state, self->upper_bound);
    x = mpz_get_ui(self->x);
    mpfr_urandomb(self->y, state);
  } while (mpfr_cmp(self->y, self->rho[x]) >= 0);

  mpz_set_ui(rop, x);
  if(dgs_bern_uniform_call(self->B, state))
    mpz_neg(rop, rop);
  mpz_add(rop, rop, self->c_z);
}

 void dgs_disc_gauss_mp_call_uniform_table_offset(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state) {
  unsigned long x;
  do {
    mpz_urandomm(self->x, state, self->two_upper_bound_minus_one);
    x = mpz_get_ui(self->x);
    mpfr_urandomb(self->y, state);
  } while (mpfr_cmp(self->y, self->rho[x]) >= 0);

  mpz_set_ui(rop, x);
  mpz_sub(rop, rop, self->upper_bound_minus_one);
  mpz_add(rop, rop, self->c_z);
}

void dgs_disc_gauss_mp_call_uniform_online(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state) {
  do {
    mpz_urandomm(self->x, state, self->two_upper_bound_minus_one);
    mpz_sub(self->x, self->x, self->upper_bound_minus_one);
    mpfr_set_z(self->z, self->x, MPFR_RNDN);
    mpfr_sub(self->z, self->z, self->c_r, MPFR_RNDN);
    mpfr_mul(self->z, self->z, self->z, MPFR_RNDN);
    mpfr_mul(self->z, self->z, self->f, MPFR_RNDN);
    mpfr_exp(self->z, self->z, MPFR_RNDN);
    mpfr_urandomb(self->y, state);
  } while (mpfr_cmp(self->y, self->z) >= 0);

  mpz_set(rop, self->x);
  mpz_add(rop, rop, self->c_z);
}


void dgs_disc_gauss_mp_call_uniform_logtable(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state) {
  do {
    mpz_urandomm(self->x, state, self->two_upper_bound_minus_one);
    mpz_sub(self->x, self->x, self->upper_bound_minus_one);
    mpz_mul(self->x2, self->x, self->x);
  } while (dgs_bern_exp_mp_call(self->Bexp, self->x2, state) == 0);
  mpz_set(rop, self->x);
  mpz_add(rop, rop, self->c_z);
}

void dgs_disc_gauss_mp_call_sigma2_logtable(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state) {
  do {
    do {
      dgs_disc_gauss_sigma2p_mp_call(self->x, self->D2, state);
      mpz_urandomm(self->y_z, state, self->k);
      mpz_mul(self->x2, self->k, self->x);
      mpz_mul_ui(self->x2, self->x2, 2);
      mpz_add(self->x2, self->x2, self->y_z);
      mpz_mul(self->x2, self->x2, self->y_z);
    } while (dgs_bern_exp_mp_call(self->Bexp, self->x2, state) == 0);
    mpz_mul(rop, self->k, self->x);
    mpz_add(rop, rop, self->y_z);
    if (mpz_sgn(rop) == 0) {
      if (dgs_bern_uniform_call(self->B, state))
        break;
    } else {
      break;
    }
  } while (1);
  if(dgs_bern_uniform_call(self->B, state))
    mpz_neg(rop, rop);
  mpz_add(rop, rop, self->c_z);
}

/** GENERAL SIGMA :: CLEAR **/

void dgs_disc_gauss_mp_clear(dgs_disc_gauss_mp_t *self) {
  mpfr_clear(self->sigma);
  if (self->B) dgs_bern_uniform_clear(self->B);
  if (self->Bexp) dgs_bern_exp_mp_clear(self->Bexp);
  if (self->D2) dgs_disc_gauss_sigma2p_clear(self->D2);
  mpz_clear(self->x);
  mpz_clear(self->x2);
  mpz_clear(self->k);
  mpfr_clear(self->y);
  mpfr_clear(self->f);
  mpfr_clear(self->z);
  mpfr_clear(self->c);
  mpfr_clear(self->c_r);
  mpz_clear(self->c_z);
  if (self->rho) {
    for(unsigned long x=0; x<mpz_get_ui(self->upper_bound); x++) {
      mpfr_clear(self->rho[x]);
    }
    free(self->rho);
  }
  free(self);
}
