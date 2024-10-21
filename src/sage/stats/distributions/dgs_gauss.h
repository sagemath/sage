/* sage_setup: distribution = sagemath-modules
 */
/**
   Discrete Gaussians over the Integers.

   A discrete Gaussian distribution on the Integers is a distribution where the
   integer `x` is sampled with probability proportional to `exp(-(x-c)²/(2σ²))`.
   It is denoted by `D_{σ,c}` where `σ` is the width parameter (close to the
   standard deviation) and `c` is the center.

   AVAILABLE ALGORITHMS:

   - ``DGS_DISC_GAUSS_UNIFORM_TABLE`` -- classical rejection sampling, sampling
     from the uniform distribution and accepted with probability proportional to
     `\exp(-(x-c)²/(2σ²))` where `\exp(-(x-c)²/(2σ²))` is precomputed and
     stored in a table. Any real-valued `c` is supported.

 - ``DGS_DISC_GAUSS_UNIFORM_LOGTABLE`` -- samples are drawn from a uniform
   distribution and accepted with probability proportional to
   `\exp(-(x-c)²/(2σ²))` where `\exp(-(x-c)²/(2σ²))` is computed using
   logarithmically many calls to Bernoulli distributions. Only integer-valued
   `c` are supported.

 - ``DGS_DISC_GAUSS_UNIFORM_ONLINE`` -- samples are drawn from a uniform
   distribution and accepted with probability proportional to
   `\exp(-(x-c)²/(2σ²))` where `\exp(-(x-c)²/(2σ²))` is computed in each
   invocation. Typically this is very slow. Any real-valued `c` is accepted.

  - ``DGS_DISC_SIGMA2_LOGTABLE`` -- samples are drawn from an easily samplable
    distribution with `σ = k·σ₂` where `σ₂ := \sqrt{1/(2\log 2)}` and
    accepted  with probability proportional to `\exp(-(x-c)²/(2σ²))` where
    `\exp(-(x-c)²/(2σ²))` is computed using logarithmically many calls to
    Bernoulli distributions (but no calls to `\exp`). Note that this sampler
    adjusts sigma to match `σ₂·k` for some integer `k`.  Only integer-valued
    `c` are supported.

  AVAILABLE PRECISIONS:

  - ``mp`` -- multi-precision using MPFR, cf. ``dgs_gauss_mp.c``

  - ``dp`` -- double precision using machine doubles, cf. ``dgs_gauss_dp.c``

  For readers unfamiliar with the implemented algorithms it makes sense to start
  with ``dgs_gauss_dp.c`` which implements the same algorithms as
  ``dgs_gauss_mp.c`` should be easier to read.

  TYPICAL USAGE::

      dgs_disc_gauss_dp_t *D = dgs_disc_gauss_dp_init(<sigma>, <c>, <tau>, <algorithm>);
      D->call(D); // as often as needed
      dgs_disc_gauss_dp_clear(D);

   .. author:: Martin R. Albrecht <martinralbrecht+dgs@googlemail.com>

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

#ifndef DGS_GAUSS__H
#define DGS_GAUSS__H

#include "dgs_bern.h"

/** UTILITY FUNCTIONS **/

/**
   We consider a double ``x`` an integer if ``fmod(x,1.0) <= DGS_DISC_GAUSS_INTEGER_CUTOFF``

   .. NOTE::

       it is okay put 0.0 here as for typical inputs the above inequality holds
       exactly
*/

#define DGS_DISC_GAUSS_INTEGER_CUTOFF 0.0

/**
   We consider two doubles ``x`` and ``y`` equal if ``abs(x-y) <= DGS_DISC_GAUSS_EQUAL_DIFF``

   .. NOTE::

       the value picked here is somewhat arbitrary
*/

#define DGS_DISC_GAUSS_EQUAL_DIFF 0.001


typedef enum {
  DGS_DISC_GAUSS_UNIFORM_ONLINE    = 0x1, //<call dgs_disc_gauss_mp_call_uniform_online
  DGS_DISC_GAUSS_UNIFORM_TABLE     = 0x2, //<call dgs_disc_gauss_mp_call_uniform_table
  DGS_DISC_GAUSS_UNIFORM_LOGTABLE  = 0x3, //<call dgs_disc_gauss_mp_call_uniform_logtable
  DGS_DISC_GAUSS_SIGMA2_LOGTABLE   = 0x7, //<call dgs_disc_gauss_mp_call_sigma2_logtable
} dgs_disc_gauss_alg_t;


/**
   Discrete Gaussian `D_{σ₂,0}` with `σ₂ := sqrt(1/(2·log(2)))`.

   Return integer `x` with probability

   `ρ_{σ,c}(x) = exp(-(x-c)²/(2σ₂²))/exp(-(\ZZ-c)²/(2σ₂²))`

   where `exp(-(\ZZ-c)²/(2σ₂²)) ≈ \sum_{i=-τσ₂}^{τσ₂} exp(-(i-c)²/(2σ₂²))` is the
   probability for all of the integers.

*/

typedef struct {
  dgs_bern_uniform_t *B; //< Bernoulli sub-samplers
} dgs_disc_gauss_sigma2p_t;

/**
   Create a new sampler for `D_{σ₂,0}`.
*/

dgs_disc_gauss_sigma2p_t *dgs_disc_gauss_sigma2p_init(void);

/**
   Return an ``mpz_t`` sampled from `D_{σ₂,0}`.

   :param rop: target value.
   :param self: discrete Gaussian sampler.
   :param state: entropy pool.

*/

void dgs_disc_gauss_sigma2p_mp_call(mpz_t rop, dgs_disc_gauss_sigma2p_t *self, gmp_randstate_t state);

/**
   Return a ``long`` sampled from `D_{σ₂,0}`.

   :param rop: target value.
   :param self: discrete Gaussian sampler.
   :param state: entropy pool.

*/

long dgs_disc_gauss_sigma2p_dp_call(dgs_disc_gauss_sigma2p_t *self);

/**
   Free `D_{σ₂,0}` sampler.

   :param self: discrete Gaussian sampler state.

*/

void dgs_disc_gauss_sigma2p_clear(dgs_disc_gauss_sigma2p_t *self);

/**
   Double-precision Discrete Gaussians `D_{σ,c}`

   Return integer `x` with probability

   `ρ_{σ,c}(x) = exp(-(x-c)²/(2σ²))/exp(-(\ZZ-c)²/(2σ²))`

   where `exp(-(\ZZ-c)²/(2σ²)) ≈ \sum_{i=-τσ}^{τσ} exp(-(i-c)²/(2σ²))` is the
   probability for all of the integers.
 */

struct _dgs_disc_gauss_mp_t;

typedef struct _dgs_disc_gauss_dp_t {

  /**
     The width parameter `σ`, i.e. samples are accepted with probability
     proportional to `\exp(-(x-c)²/(2σ²))`
  */

  double sigma;

  /**
     The mean of the distribution `c`. The value of `c` does not have to be an
     integer. However, some algorithms only support integer-valued `c`.
  */

  double c;
  double c_r;  //< `c_r := c % 1`
  long   c_z; //< c_z := c - (c_r)

  /**
     Cutoff `τ`, samples outside the range `(⌊c⌉-⌈στ⌉,...,⌊c⌉+⌈στ⌉)` are
     considered to have probability zero. This bound applies to algorithms
     which sample from the uniform distribution.
  */

  size_t tau;

  dgs_disc_gauss_alg_t algorithm;  //<  which algorithm to use

  /**
     We use a uniform Bernoulli to decide signs.
   */

  dgs_bern_uniform_t *B;

  /**
     To realize rejection sampling, we call `B_{exp(-(x·x)/(2σ²))}` and accept
     if it returns 1.

     Used when ``DGS_DISC_GAUSS_UNIFORM_LOGTABLE`` or
     ``DGS_DISC_GAUSS_SIGMA2_LOGTABLE`` is set.
   */

  dgs_bern_exp_dp_t *Bexp;

  /**
     `\D_{σ₂,0}` which is easily sampable`

     Used when ``DGS_DISC_GAUSS_SIGMA2_LOGTABLE`` is set.
  */

  dgs_disc_gauss_sigma2p_t *D2;

  /**
   Return an ``long`` sampled from this sampler

   :param rop: target value.
   :param self: discrete Gaussian sampler.
   :param state: entropy pool.

  */

  long (*call)(struct _dgs_disc_gauss_dp_t *self);

  /**
   We sample ``x`` with ``abs(x) < upper_bound`` in
   ``DGS_DISC_GAUSS_UNIFORM_ONLINE``, ``DGS_DISC_GAUSS_UNIFORM_TABLE`` and
   ``DGS_DISC_GAUSS_UNIFORM_LOGTABLE``.
   */

  long upper_bound;

  /**
   We sample ``x`` with ``abs(x) <= upper_bound - 1`` in
   ``DGS_DISC_GAUSS_UNIFORM_ONLINE``, ``DGS_DISC_GAUSS_UNIFORM_TABLE`` and
   ``DGS_DISC_GAUSS_UNIFORM_LOGTABLE``.
  */

  long upper_bound_minus_one;

  /**
     There are ``2*upper_bound -1`` elements in the range
     ``-upper_bound+1,...,upper_bound-1``.
  */

  long two_upper_bound_minus_one;

  /**
     The multiplier `k` when we sample from `D_{k·σ₂,c}` in
     ``DGS_DISC_GAUSS_SIGMA2_LOGTABLE``.
  */

  long k;


  /**
   Precomputed `-1/(2σ²)`.
  */

  double f;

  /**
     Precomputed values for `exp(-(x-2)²/(2σ²))` in
     ``DGS_DISC_GAUSS_UNIFORM_TABLE``
  */

  double *rho;
} dgs_disc_gauss_dp_t;

/**
 Create a new double-precision discrete Gaussian sampler.

 :param sigma: width parameter `σ`
 :param c: center `c`
 :param tau: cutoff `τ`
 :param algorithm: algorithm to use.

*/

dgs_disc_gauss_dp_t *dgs_disc_gauss_dp_init(double sigma, double c, size_t tau, dgs_disc_gauss_alg_t algorithm);

/**
   Sample from ``dgs_disc_gauss_dp_t`` by rejection sampling using the uniform distribution

   :param self: Discrete Gaussian sampler

 */

long dgs_disc_gauss_dp_call_uniform_online(dgs_disc_gauss_dp_t *self);

/**
   Sample from ``dgs_disc_gauss_dp_t`` by rejection sampling using the uniform
   distribution and tabulated ``exp()`` evaluations.

   :param self: discrete Gaussian sampler

   .. NOTE::

      `c` must be an integer in this algorithm
 */

long dgs_disc_gauss_dp_call_uniform_table(dgs_disc_gauss_dp_t *self);

/**
   Sample from ``dgs_disc_gauss_dp_t`` by rejection sampling using the uniform
   distribution and tabulated ``exp()`` evaluations.

   :param self: discrete Gaussian sampler

   .. NOTE::

      This function makes no assumptions about `c` but requires more resources
      than ``dgs_disc_gauss_dp_call_uniform_table()``.
 */

long dgs_disc_gauss_dp_call_uniform_table_offset(dgs_disc_gauss_dp_t *self);


/**
  Sample from ``dgs_disc_gauss_dp_t`` by rejection sampling using the uniform
  distribution replacing all ``exp()`` calls with calls to Bernoulli
  distributions.

  :param self: discrete Gaussian sampler

  .. NOTE::

      `c` must be an integer in this algorithm
 */

long dgs_disc_gauss_dp_call_uniform_logtable(dgs_disc_gauss_dp_t *self);

/**
   Sample from ``dgs_disc_gauss_dp_t`` by rejection sampling using the
   ``D_{k·σ₂,0}` distribution replacing all ``exp()`` calls with calls to
   Bernoulli distributions.

   :param self: discrete Gaussian sampler

   .. NOTE::

      `c` must be an integer in this algorithm

 */

long dgs_disc_gauss_dp_call_sigma2_logtable(dgs_disc_gauss_dp_t *self);

/**
   The uniform Bernoulli sampler which is used to decide signs caches bits for
   performance reasons. This functions clears this cache of random bits.

   :param self: discrete Gaussian sampler

 */

static inline void dgs_disc_gauss_dp_flush_cache(dgs_disc_gauss_dp_t *self) {
  self->B->count = self->B->length;
}

/**
   Free memory.

   :param self: discrete Gaussian sadpler

 */

void dgs_disc_gauss_dp_clear(dgs_disc_gauss_dp_t *self);


/**
   Multi-precision Discrete Gaussians `D_{σ,c}`

   Return integer `x` with probability

   `ρ_{σ,c}(x) = exp(-(x-c)²/(2σ²))/exp(-(\ZZ-c)²/(2σ²))`

   where `exp(-(\ZZ-c)²/(2σ²)) ≈ \sum_{i=-τσ}^{τσ} exp(-(i-c)²/(2σ^²))` is the
   probability for all of the integers.

*/

typedef struct _dgs_disc_gauss_mp_t {

  /**
      The width parameter `σ`, i.e. samples are accepted with probability
      proportional to `\exp(-(x-c)²/(2σ²))`
   */

  mpfr_t sigma;

  /**
     The mean of the distribution `c`. The value of `c` does not have to be an
     integer. However, some algorithms only support integer-valued `c`.
  */

  mpfr_t c;

  mpfr_t c_r; //< `c_r := c % 1`
  mpz_t c_z;  //< c_z := c - (c_r)

  /**
     Cutoff `τ`, samples outside the range `(⌊c⌉-⌈στ⌉,...,⌊c⌉+⌈στ⌉)` are
     considered to have probability zero. This bound applies to algorithms
     which sample from the uniform distribution.
  */

  size_t tau;

  dgs_disc_gauss_alg_t algorithm; //<  which algorithm to use

  /**
     We use a uniform Bernoulli to decide signs.
   */

  dgs_bern_uniform_t *B;

  /**
     To realize rejection sampling, we call `B_{exp(-(x·x)/(2σ²))}` and accept
     if it returns 1.

     Used when ``DGS_DISC_GAUSS_UNIFORM_LOGTABLE`` or
     ``DGS_DISC_GAUSS_SIGMA2_LOGTABLE`` is set.
   */

  dgs_bern_exp_mp_t *Bexp;


  /**
     `D_{σ₂,0}` which is easily sampable`

     Used when ``DGS_DISC_GAUSS_SIGMA2_LOGTABLE`` is set.
  */

  dgs_disc_gauss_sigma2p_t *D2;

  /**
   Return an ``mpz_t`` sampled from this sampler

   :param rop: target value.
   :param self: discrete Gaussian sampler.
   :param state: entropy pool.

  */

  void (*call)(mpz_t rop, struct _dgs_disc_gauss_mp_t *self, gmp_randstate_t state);

  /**
   We sample ``x`` with ``abs(x) < upper_bound`` in
   ``DGS_DISC_GAUSS_UNIFORM_ONLINE``, ``DGS_DISC_GAUSS_UNIFORM_TABLE`` and
   ``DGS_DISC_GAUSS_UNIFORM_LOGTABLE``.
   */

  mpz_t upper_bound;

  /**
   We sample ``x`` with ``abs(x) <= upper_bound - 1`` in
   ``DGS_DISC_GAUSS_UNIFORM_ONLINE``, ``DGS_DISC_GAUSS_UNIFORM_TABLE`` and
   ``DGS_DISC_GAUSS_UNIFORM_LOGTABLE``.
   */

  mpz_t upper_bound_minus_one;

  /**
     There are ``2*upper_bound -1`` elements in the range
     ``-upper_bound+1,...,upper_bound-1``.
   */

  mpz_t two_upper_bound_minus_one;

  /**
     The multiplier `k` when we sample from `D_{k·σ₂,c}` in
     ``DGS_DISC_GAUSS_SIGMA2_LOGTABLE``.
  */

  mpz_t k;

  /**
   Precomputed `-1/(2σ²)`.
  */

  mpfr_t f;

  mpz_t x; //< space for temporary integer
  mpz_t y_z; //< space for temporary integer
  mpz_t x2; // space for temporary integer
  mpfr_t y; // space for temporary rational number
  mpfr_t z; // space for temporary rational number

  /**
     Precomputed values for `exp(-(x-c)²/(2σ²))` in
     ``DGS_DISC_GAUSS_UNIFORM_TABLE``
  */

  mpfr_t *rho;

} dgs_disc_gauss_mp_t;

dgs_disc_gauss_mp_t *dgs_disc_gauss_mp_init(mpfr_t sigma, mpfr_t c, size_t tau, dgs_disc_gauss_alg_t algorithm);

/**
   Sample from ``dgs_disc_gauss_mp_t`` by rejection sampling using the uniform
   distribution and tabulated ``exp()`` evaluations.

   :param self: discrete Gaussian sampler

   .. NOTE::

      `c` must be an integer in this algorithm

 */

void dgs_disc_gauss_mp_call_uniform_table(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);

/**
   Sample from ``dgs_disc_gauss_mp_t`` by rejection sampling using the uniform
   distribution and tabulated ``exp()`` evaluations.

   :param self: discrete Gaussian sampler

   .. NOTE::

      This function makes no assumptions about `c` but requires more resources
      than ``dgs_disc_gauss_dp_call_uniform_table()``.

 */

void dgs_disc_gauss_mp_call_uniform_table_offset(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);

/**
  Sample from ``dgs_disc_gauss_mp_t`` by rejection sampling using the uniform
  distribution replacing all ``exp()`` calls with call to Bernoulli distributions.

  :param self: discrete Gaussian sampler

  .. NOTE::

     `c` must be an integer in this algorithm
 */

void dgs_disc_gauss_mp_call_uniform_logtable(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);

/**
  Sample from ``dgs_disc_gauss_mp_t`` by rejection sampling using the uniform distribution.

  :param self: discrete Gaussian sampler

 */

void dgs_disc_gauss_mp_call_uniform_online(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);

/**
  Sample from ``dgs_disc_gauss_mp_t`` by rejection sampling using the `D_{k·σ₂,0}`
  distribution replacing all ``exp()`` calls with call to Bernoulli distributions.

  :param self: Discrete Gaussian sampler

  .. NOTE::

     `c` must be an integer in this algorithm.
 */

void dgs_disc_gauss_mp_call_sigma2_logtable(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);

/**
   Clear cache of random bits.

   :param self: discrete Gaussian sampler

 */

static inline void dgs_disc_gauss_mp_flush_cache(dgs_disc_gauss_mp_t *self) {
  self->B->count = self->B->length;
}

/**
   Free memory.

   :param self: discrete Gaussian sadpler

 */

void dgs_disc_gauss_mp_clear(dgs_disc_gauss_mp_t *self);

#endif //DGS_GAUSS__H
