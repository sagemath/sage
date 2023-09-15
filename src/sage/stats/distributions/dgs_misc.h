/* sage_setup: distribution = sagemath-modules
 */
/**
 * \file dgs_misc.h
 *
 * \author Martin Albrecht <martinralbrecht+dgs@googlemail.com>
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
#ifndef DGS_MISC__H
#define DGS_MISC__H

#include <stddef.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <stdarg.h>

/**
 * \brief Macro to help with branch prediction.
 */

#define __DGS_LIKELY(cond)    __builtin_expect ((cond) != 0, 1)

/**
 * \brief Macro to help with branch prediction.
 */

#define __DGS_UNLIKELY(cond)  __builtin_expect ((cond) != 0, 0)


static int const dgs_radix = sizeof(unsigned long)<<3;
static unsigned long const dgs_ffff = -1;

#define __DGS_LSB_BITMASK(n) (dgs_ffff >> (dgs_radix - (n)) % dgs_radix)


static inline unsigned long _dgs_randomb_libc(size_t nbits) {
  size_t n = __DGS_LSB_BITMASK(nbits);
  assert(((RAND_MAX | (RAND_MAX >> 1)) == RAND_MAX));
  if (__DGS_LIKELY(n <= RAND_MAX))
    return random() & n;
  assert(RAND_MAX >= __DGS_LSB_BITMASK(22));
  unsigned long pool = (((unsigned long)random()) << 0) ^ (((unsigned long)random()) << 22) ^ (((unsigned long)random()) << 44);
  return pool & n;
}

static inline unsigned long _dgs_randomm_libc(unsigned long n) {
  assert(n < RAND_MAX);
  long r;
  unsigned long k = RAND_MAX/n;
  do {
    r = random();
  } while (r >= k*n);
  return r%n;
}

static inline void dgs_die(const char *msg, ...) {
  va_list lst;
  va_start(lst, msg);
  vfprintf(stderr, msg, lst);
  fprintf(stderr, "\n");
  va_end(lst);
  abort();
}

#endif //DGS_MISC__H
