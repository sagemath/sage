/* sage_setup: distribution = sagemath-modules
 */
/**
   Bernoulli samplers.

   .. author:: Martin Albrecht <martinralbrecht+dgs@googlemail.com>
 */

/******************************************************************************
*
*                      DGS - Discrete Gaussian Samplers
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

#ifndef DGS_BERN__H
#define DGS_BERN__H

#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include <gmp.h>
#include <mpfr.h>

#include "dgs_misc.h"

/**
   Number of bits sampled at once in ``dgs_bern_uniform_t``

   .. NOTE::

      We want machine independent behaviour at least for the MP functions, so we pick 32.
*/

#define DGS_BERN_UNIFORM_DEFAULT_LENGTH 32

/**
   Maximum number of bits sampled at once in ``dgs_bern_uniform_t``
*/

#define DGS_BERN_UNIFORM_MAX_LENGTH (sizeof(unsigned long)*8)

/**
   Number of blocks allocated in one go in ``dgs_bern_exp_mp_t`` and ``dgs_bern_exp_dp_t``
*/

#define DGS_BERN_EXP_ALLOC_BLOCK_SIZE 16

/**
   Balanced Bernoulli distribution.

   .. NOTE::

       This implementation is faster than ``dgs_bern_mp_t``

*/

typedef struct {

  /**
     Number of bits we sample in each go.
  */

  size_t   length;

  /**
     Number of bits we have consumed yet.
  */

  size_t   count;

  /**
     We sample to this ``mpz_t``.
  */

  mpz_t    tmp;

  /**
     We store the pool of random bits here.
  */

  unsigned long pool;
} dgs_bern_uniform_t;

/**
   Construct a new uniform bit sampler.

   :param length: number of bits to sample from GMP (or 0 for automatic choice)

   .. NOTE::

       Clear with ``dgs_bern_uniform_clear()``.

*/

dgs_bern_uniform_t* dgs_bern_uniform_init(size_t length);

/**
   Sample a new uniformly random bit.

   :param self: Bernoulli state
   :param state: GMP randstate used as randomness source

*/

static inline long dgs_bern_uniform_call(dgs_bern_uniform_t *self, gmp_randstate_t state) {
  assert(self != NULL);
  assert(state != NULL);

  if (__DGS_UNLIKELY(self->count == self->length)) {
    mpz_urandomb(self->tmp, state, self->length);
    self->pool = mpz_get_ui(self->tmp);
    self->count = 0;
  }

  long b = (long)(self->pool & 1);
  self->pool >>= 1;
  self->count++;
  return b;
}

/**
   Sample a new uniformly random bit using libc ``random()``.

   :param self: Bernoulli state

 */

static inline long dgs_bern_uniform_call_libc(dgs_bern_uniform_t *self) {
  assert(self != NULL);
  if (__DGS_UNLIKELY(self->count == self->length)) {
    self->pool = _dgs_randomb_libc(self->length);
    self->count = 0;
  }

  long b = (long)(self->pool & 1);
  self->pool >>= 1;
  self->count++;
  return b;
}

/**
   Clear cache of random bits.

   :param self: Bernoulli state

 */

static inline void dgs_bern_uniform_flush_cache(dgs_bern_uniform_t *self) {
  self->count = self->length;
}

/**
   Clear uniformly random bit sampler.

   :param self: Bernoulli state

*/

void dgs_bern_uniform_clear(dgs_bern_uniform_t *self);

/**
   Multi-precision Bernoulli distribution.
 */

typedef struct {
  /**
     Return 1 with probability `p`
  */
  mpfr_t p;
  mpfr_t tmp; //< used internally
} dgs_bern_mp_t;

/**
   Create new Bernoulli sampler.

   :param p: sampler will return 1 with probability `p` and 0 with probability `1-p`.

   .. NOTE::

       Clear with ``dgs_bern_mp_clear()``.

 */

dgs_bern_mp_t *dgs_bern_mp_init(mpfr_t p);

/**
   Return 1 with probability `p`.

   :param self: Bernoulli state
   :param state: GMP randstate used as randomness source

 */

long dgs_bern_mp_call(dgs_bern_mp_t *self, gmp_randstate_t state);

/**
   Clear Bernoulli sampler.

   :param self: Bernoulli state

 */

void dgs_bern_mp_clear(dgs_bern_mp_t *self);

/**
   double precision Bernoulli distribution.

 */

typedef struct {
  /**
     Return 1 with probability `p`
  */
  double p;
} dgs_bern_dp_t;

/**

   Create new Bernoulli sampler.

   :param p: sampler will return 1 with probability `p` and 0 with probability `1-p`.

   .. NOTE::

       Clear with ``dgs_bern_dp_clear()``.

 */

dgs_bern_dp_t *dgs_bern_dp_init(double p);


/**
   Return 1 with probability `p`.

   :param self: Bernoulli state

  .. NOTE::

      Uses libc ``random()``.

 */

long dgs_bern_dp_call(dgs_bern_dp_t *self);

/**
   Clear Bernoulli sampler.

   :param self: Bernoulli state
 */

void dgs_bern_dp_clear(dgs_bern_dp_t *self);


/**
   Multi-precision Bernoulli distribution with `p = exp(-x/f)` for positive integers `x`.
 */

typedef struct {
  /**
     We support positive `x` up to `2^l-1`
   */

  size_t l;

  /**
     Probabilities for Bernoulli sub-samplers.
   */

  mpfr_t *p;

  /**
     Bernoulli sub-samplers.
  */

  dgs_bern_mp_t **B;

} dgs_bern_exp_mp_t;

/**
   Create new family of Bernoulli samples.

   :param f: samplers return 1 with probability `exp(-x/f)`
   :param l: we support inputs to call `x` up to `2^l-1`

   .. NOTE::

       Clear ``dgs_bern_exp_mp_clear()``.

 */

dgs_bern_exp_mp_t* dgs_bern_exp_mp_init(mpfr_t f, size_t l);


/**
   Return 1 with probability `exp(-x/f)`.

   :param self: Bernoulli state
   :param x: integer with `0 < x < 2^l`
   :param state: GMP randstate used as randomness source

 */
long dgs_bern_exp_mp_call(dgs_bern_exp_mp_t *self, mpz_t x, gmp_randstate_t state);


/**
   Clear Bernoulli sampler family.

   :param self: Bernoulli state

 */

void dgs_bern_exp_mp_clear(dgs_bern_exp_mp_t *self);

/**
   double-precision Bernoulli distribution with `p = exp(-x/f)` for positive integers `x`.
 */

typedef struct {
  /**
     We support positive `x` up to `2^l-1`
   */

  size_t l;

  /**
     Probabilities for Bernoulli sub-samplers.
   */

  double *p;

  /**
     Bernoulli sub-samplers.
  */

  dgs_bern_dp_t **B;

} dgs_bern_exp_dp_t;

/**
   Create new family of Bernoulli samplers.

   :param f: samplers return 1 with probability `exp(-x/f)`
   :param l: we support inputs to call `x` up to `2^l-1`

   .. NOTE::

       Clear ``dgs_bern_exp_mp_clear()``.

 */

dgs_bern_exp_dp_t* dgs_bern_exp_dp_init(double f, size_t l);

/**
   Return 1 with probability `exp(-x/f)`.

   :param self: Bernoulli state
   :param x: integer with `0 < x < 2^l`
   :param state: GMP randstate used as randomness source

 */

long dgs_bern_exp_dp_call(dgs_bern_exp_dp_t *self, long x);

/**
   Clear Bernoulli sampler family.

   :param self: Bernoulli state

 */

void dgs_bern_exp_dp_clear(dgs_bern_exp_dp_t *self);

#endif //DGS_BERN__H
