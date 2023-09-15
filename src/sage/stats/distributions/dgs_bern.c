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
#include <math.h>

/*
 * balanced Bernoulli distribution, machine-precision version
 */

dgs_bern_uniform_t* dgs_bern_uniform_init(size_t length) {
  if (length == 0) {
    length = DGS_BERN_UNIFORM_DEFAULT_LENGTH;
  }
  assert(length <= DGS_BERN_UNIFORM_MAX_LENGTH);

  dgs_bern_uniform_t *self = (dgs_bern_uniform_t *)malloc(sizeof(dgs_bern_uniform_t));
  if (!self) dgs_die("out of memory");
  self->length = length;

  self->count = self->length;
  mpz_init(self->tmp);
  return self;
}


void dgs_bern_uniform_clear(dgs_bern_uniform_t *self) {
  mpz_clear(self->tmp);
  free(self);
}

/*
 * Bernoulli distribution, multi-precision version
 */

dgs_bern_mp_t* dgs_bern_mp_init(mpfr_t p) {
  /* we allow 0 here for low precision */
  assert((mpfr_cmp_d(p, 0.0) >= 0) && (mpfr_cmp_d(p, 1.0) < 0));

  dgs_bern_mp_t *self = malloc(sizeof(dgs_bern_mp_t));
  if (!self) dgs_die("out of memory");

  mpfr_init_set(self->p, p, MPFR_RNDN);
  mpfr_init2(self->tmp, mpfr_get_prec(p));
  return self;
}

long dgs_bern_mp_call(dgs_bern_mp_t *self, gmp_randstate_t state) {
  mpfr_urandomb(self->tmp, state);

  if (mpfr_cmp(self->tmp, self->p)<0) {
    return 1;
  } else {
    return 0;
  }
}

void dgs_bern_mp_clear(dgs_bern_mp_t *self) {
  mpfr_clear(self->tmp);
  mpfr_clear(self->p);
  free(self);
}

/*
 * Bernoulli with p = exp(-x/f) for integers x, multi-precision version
 */

dgs_bern_exp_mp_t* dgs_bern_exp_mp_init(mpfr_t f, size_t l) {
  dgs_bern_exp_mp_t *self = (dgs_bern_exp_mp_t *)malloc(sizeof(dgs_bern_exp_mp_t));
  if (!self) dgs_die("out of memory");

  /* l == 0, means we use the precision of f to decide l */
  if (l == 0)
    l = SIZE_MAX;

  self->l = DGS_BERN_EXP_ALLOC_BLOCK_SIZE;
  self->p = (mpfr_t*)malloc(sizeof(mpfr_t)*self->l);
  if (!self->p) dgs_die("out of memory");
  self->B = (dgs_bern_mp_t**)malloc(sizeof(dgs_bern_mp_t)*self->l);
  if (!self->B) dgs_die("out of memory");

  mpfr_t tmp, tmp2;
  mpfr_init2(tmp2, mpfr_get_prec(f));
  mpfr_init_set(tmp, f, MPFR_RNDN); // f
  mpfr_pow_si(tmp, tmp, -1, MPFR_RNDN); // 1/f
  mpfr_neg(tmp, tmp, MPFR_RNDN); // -1/f

  for(size_t i=0; i<l; i++) {
    mpfr_exp(tmp2, tmp, MPFR_RNDN);
    if (mpfr_zero_p(tmp2)) {
      self->l = i+1;
      break;
    }
    if (i%DGS_BERN_EXP_ALLOC_BLOCK_SIZE == 0 && i!=0) {
      self->l += DGS_BERN_EXP_ALLOC_BLOCK_SIZE;
      self->l = (l>self->l) ? self->l : l;
      self->p = realloc(self->p, sizeof(mpfr_t)*self->l);
      if(!self->p) dgs_die("out of memory");
      self->B = realloc(self->B, sizeof(dgs_bern_exp_mp_t)*self->l);
      if(!self->B) dgs_die("out of memory");
    }

    mpfr_init_set(self->p[i], tmp2, MPFR_RNDN);
    self->B[i] = dgs_bern_mp_init(self->p[i]);

    mpfr_mul_ui(tmp, tmp, 2, MPFR_RNDN);
  }
  if (l < self->l)
    self->l = l;
  mpfr_clear(tmp);
  mpfr_clear(tmp2);
  return self;
}

long dgs_bern_exp_mp_call(dgs_bern_exp_mp_t *self, mpz_t x, gmp_randstate_t state) {
  assert(mpz_sgn(x) >= 0);
  long int start = (mpz_sizeinbase(x, 2) < self->l) ? mpz_sizeinbase(x, 2) : self->l;

  for(long int i=start-1; i>=0; i--) {
    if (mpz_tstbit(x, i)) {
      if (dgs_bern_mp_call(self->B[i], state) == 0) {
        return 0;
      }
    }
  }
  return 1;
}

void dgs_bern_exp_mp_clear(dgs_bern_exp_mp_t *self) {
  if(!self)
    return;

  for(size_t i=0; i<self->l; i++) {
    mpfr_clear(self->p[i]);
    dgs_bern_mp_clear(self->B[i]);
  }
  if(self->p)
    free(self->p);
  if(self->B)
    free(self->B);
  free(self);
}

/*
 * Bernoulli distribution, double-precision version
 */

dgs_bern_dp_t* dgs_bern_dp_init(double p) {
  /* we allow 0 here for low precision */
  assert((p >= 0) && (p < 1));

  dgs_bern_dp_t *self = malloc(sizeof(dgs_bern_dp_t));
  if (!self) dgs_die("out of memory");

  self->p = p;
  return self;
}

long dgs_bern_dp_call(dgs_bern_dp_t *self) {
  double c = drand48();
  if (c<self->p)
    return 1;
  else
    return 0;
}

void dgs_bern_dp_clear(dgs_bern_dp_t *self) {
  free(self);
}

/*
 * Bernoulli with p = exp(-x/f) for integers x, multi-precision version
 */

dgs_bern_exp_dp_t* dgs_bern_exp_dp_init(double f, size_t l) {
  dgs_bern_exp_dp_t *self = (dgs_bern_exp_dp_t *)malloc(sizeof(dgs_bern_exp_dp_t));
  if (!self) dgs_die("out of memory");

  /* l == 0, means we use the precision of f to decide l */
  if (l == 0)
    l = SIZE_MAX;

  self->l = DGS_BERN_EXP_ALLOC_BLOCK_SIZE;
  self->p = (double*)malloc(sizeof(double)*self->l);
  if (!self->p) dgs_die("out of memory");
  self->B = (dgs_bern_dp_t**)malloc(sizeof(dgs_bern_dp_t)*self->l);
  if (!self->B) dgs_die("out of memory");

  double tmp = -1.0/f;
  double tmp2;

  for(size_t i=0; i<l; i++) {
    tmp2 = exp(tmp);
    if (tmp2 == 0.0) {
      self->l = i+1;
      break;
    }
    if (i%DGS_BERN_EXP_ALLOC_BLOCK_SIZE == 0 && i!=0) {
      self->l += DGS_BERN_EXP_ALLOC_BLOCK_SIZE;
      self->l = (l>self->l) ? self->l : l;
      self->p = realloc(self->p, sizeof(double)*self->l);
      if(!self->p) dgs_die("out of memory");
      self->B = realloc(self->B, sizeof(dgs_bern_exp_dp_t)*self->l);
      if(!self->B) dgs_die("out of memory");
    }

    self->p[i] = tmp2;
    self->B[i] = dgs_bern_dp_init(self->p[i]);

    tmp = 2*tmp;
  }
  if (l < self->l)
    self->l = l;
  return self;
}

long dgs_bern_exp_dp_call(dgs_bern_exp_dp_t *self, long x) {
  if (x == 0)
    return 1;
  assert(x >= 0);
  long start = self->l;
  for(long i=start-1; i>=0; i--) {
    if (x & (1L<<i)) {
      if (dgs_bern_dp_call(self->B[i]) == 0) {
        return 0;
      }
    }
  }
  return 1;
}

void dgs_bern_exp_dp_clear(dgs_bern_exp_dp_t *self) {
  if(!self)
    return;

  for(size_t i=0; i<self->l; i++) {
    dgs_bern_dp_clear(self->B[i]);
  }
  if(self->p)
    free(self->p);
  if(self->B)
    free(self->B);
  free(self);
}
